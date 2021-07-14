#include "Strategy.h"

Plan::Plan() : Plan(Default()) {}
Plan::Plan(const Plan& other)
    : _entries{CloneEntries(other)}
    , _exits{CloneExits(other)}
{
}
Plan::Plan(Plan&& other) noexcept
    : _entries{std::move(other._entries)}
    , _exits{std::move(other._exits)}
{
}
auto Plan::operator=(const Plan& other)  //
    -> Plan&
{
  _entries = CloneEntries(other);
  _exits   = CloneExits(other);
  return *this;
}
auto Plan::operator=(Plan&& other) noexcept  //
    -> Plan&
{
  _entries = std::move(other._entries);
  _exits   = std::move(other._exits);
  return *this;
}
Plan::Plan(const std::vector<std::unique_ptr<StrategyBase>>& entries,
           const std::vector<std::unique_ptr<StrategyBase>>& exits)
    : _entries{CloneVector(entries)}
    , _exits{CloneVector(exits)}
{
}

auto Plan::entries()  //
    -> std::vector<std::unique_ptr<StrategyBase>>&
{
  return _entries;
}
auto Plan::exits()  //
    -> std::vector<std::unique_ptr<StrategyBase>>&
{
  return _exits;
}
auto Plan::entries() const  //
    -> const std::vector<std::unique_ptr<StrategyBase>>&
{
  return _entries;
}
auto Plan::exits() const  //
    -> const std::vector<std::unique_ptr<StrategyBase>>&
{
  return _exits;
}

auto Plan::CloneEntries(const Plan& other)  //
    -> std::vector<std::unique_ptr<StrategyBase>>
{
  return CloneVector(other._entries);
}
auto Plan::CloneExits(const Plan& other)  //
    -> std::vector<std::unique_ptr<StrategyBase>>
{
  return CloneVector(other._exits);
}
auto Plan::Name() const  //
    -> std::string
{
  auto ret = std::string{"Plan{Enter{"};
  if (!_entries.empty())
  {
    ret += _entries[0]->Name();
    for (std::size_t i = 1; i < _entries.size(); ++i)
    {
      ret += "," + _entries[i]->Name();
    }
  }
  ret += "},Exit{";
  if (!_exits.empty())
  {
    ret += _exits[0]->Name();
    for (std::size_t i = 1; i < _exits.size(); ++i)
    {
      ret += "," + _exits[i]->Name();
    }
  }
  ret += "}}";
  return ret;
}
auto Plan::CloneVector(
    const std::vector<std::unique_ptr<StrategyBase>>& vec)  //
    -> std::vector<std::unique_ptr<StrategyBase>>
{
  auto ret = std::vector<std::unique_ptr<StrategyBase>>{};
  for (auto& v : vec)
  {
    ret.push_back(v->Clone());
  }
  return ret;
}
auto Plan::Default()  //
    -> Plan
{
  auto entries = std::vector<std::unique_ptr<StrategyBase>>{};
  auto exits   = std::vector<std::unique_ptr<StrategyBase>>{};
  entries.push_back(
      std::make_unique<
          Strategy<Identity<Metric::Open>, Identity<Metric::Close>>>());
  entries.push_back(
      std::make_unique<
          Strategy<Identity<Metric::Open>, Identity<Metric::Close>>>());
  exits.push_back(std::make_unique<
                  Strategy<Identity<Metric::Open>, Identity<Metric::Close>>>());
  exits.push_back(std::make_unique<
                  Strategy<Identity<Metric::Open>, Identity<Metric::Close>>>());
  return Plan(entries, exits);
}

auto Plan::StrategiesRun(
    const std::vector<std::unique_ptr<StrategyBase>>& strategies,
    const History&                                    candles) const  //
    -> std::vector<std::map<Candle, Position, Candle::Compare::OpenTime>>
{
  auto ret
      = std::vector<std::map<Candle, Position, Candle::Compare::OpenTime>>{};
  for (auto& s : strategies)
  {
    ret.push_back(s->Evaluate(candles));
  }
  return ret;
}
auto Plan::EntriesRun(const History& candles) const  //
    -> std::vector<std::map<Candle, Position, Candle::Compare::OpenTime>>
{
  return StrategiesRun(_entries, candles);
}
auto Plan::ExitsRun(const History& candles) const  //
    -> std::vector<std::map<Candle, Position, Candle::Compare::OpenTime>>
{
  return StrategiesRun(_exits, candles);
}
auto Plan::StrategiesRunJoin(
    const std::vector<std::map<Candle, Position, Candle::Compare::OpenTime>>&
        runs) const  //
    -> std::map<Candle, Position, Candle::Compare::OpenTime>
{
  auto ret = std::map<Candle, Position, Candle::Compare::OpenTime>{};
  for (auto& run : runs)
  {
    for (auto& [key, value] : run)
    {
      if (value == Position::Short)
      {
        ret[key] = Position::Short;
        continue;
      }
      if (!ret.contains(key))
      {
        ret[key] = Position::Long;
      }
    }
  }
  return ret;
}
auto Plan::JoinEntriesExits(
    const std::map<Candle, Position, Candle::Compare::OpenTime>& entries,
    const std::map<Candle, Position, Candle::Compare::OpenTime>& exits)
    const  //
    -> std::map<Candle, Position, Candle::Compare::OpenTime>
{
  auto ret = std::map<Candle, Position, Candle::Compare::OpenTime>{};
  auto pos = Position::Short;
  for (auto& [key, value] : entries)
  {
    if (pos == Position::Short && value == Position::Long)
    {
      pos = Position::Long;
    }
    else if (pos == Position::Long && exits.contains(key)
             && exits.at(key) == Position::Short)
    {
      pos = Position::Short;
    }
    ret[key] = pos;
  }
  return ret;
}

auto Plan::Run(const History& candles) const  //
    -> std::map<Candle, Position, Candle::Compare::OpenTime>
{
  return JoinEntriesExits(StrategiesRunJoin(EntriesRun(candles)),
                          StrategiesRunJoin(ExitsRun(candles)));
}
auto Plan::Assess(const History& candles,
                  double         loss) const  //
    -> std::map<Candle, double, Candle::Compare::OpenTime>
{
  auto ret          = std::map<Candle, double, Candle::Compare::OpenTime>{};
  auto run          = Run(candles);
  auto entry_candle = static_cast<const Candle*>(nullptr);
  auto exit_candle  = static_cast<const Candle*>(nullptr);
  auto last_pos     = Position::Short;
  if (ret.empty())
  {
    return {};
  }
  ret[run.begin()->first] = 1.0;
  for (auto& [candle, position] : run)
  {
    if (last_pos == Position::Short && position == Position::Long)
    {
      last_pos     = Position::Long;
      entry_candle = &candle;
    }
    else if (last_pos == Position::Long && position == Position::Short)
    {
      last_pos    = Position::Short;
      exit_candle = &candle;
      ret[candle] = exit_candle->close / entry_candle->close - loss;
    }
  }
  return ret;
}
auto Plan::AssessTotal(const History& candles,
                       double         loss) const  //
    -> std::map<Candle, double, Candle::Compare::OpenTime>
{
  auto assess = Assess(candles, loss);
  auto ret    = std::map<Candle, double, Candle::Compare::OpenTime>{};
  auto total  = 1.0;
  for (auto& candle : candles)
  {
    if (assess.contains(candle))
    {
      total *= assess[candle];
    }
    ret[candle] = total;
  }
  return ret;
}
auto Plan::AssessFinalTotal(const History& candles,
                            double         loss) const  //
    -> double
{
  auto assessment = AssessTotal(candles);
  if (assessment.empty())
  {
    return 0.0;
  }
  return assessment.rbegin()->second;
}
auto Plan::AssessWithEntry(double         entry_value,
                           const History& candles,
                           double         loss) const  //
    -> std::map<Candle, double, Candle::Compare::OpenTime>
{
  auto ret        = std::map<Candle, double, Candle::Compare::OpenTime>{};
  auto assessment = AssessTotal(candles, loss);
  for (auto& [key, value] : assessment)
  {
    ret[key] = value * entry_value;
  }
  return ret;
}