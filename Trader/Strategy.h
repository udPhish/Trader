#pragma once
#include <array>
#include <chrono>
#include <functional>
#include <map>
#include <set>
#include <vector>

#include "libudph/Class/Interface.h"

#include "Candle.h"
#include "Indicator.h"
#include "Logger.h"
#include "Strategies.h"

enum ComparisonType
{
  Less,
  LessOrEqual,
  Greater,
  GreaterOrEqual,
  Equal,
  NotEqual
};
template<int Comp>
struct Comparison;
template<>
struct Comparison<ComparisonType::Less>
{
  using Type = std::less<double>;
  inline static constexpr auto Name()  //
      -> const char*
  {
    return "Less";
  }
};
template<>
struct Comparison<ComparisonType::LessOrEqual>
{
  using Type = std::less_equal<double>;
  inline static constexpr auto Name()  //
      -> const char*
  {
    return "LessOrEqual";
  }
};
template<>
struct Comparison<ComparisonType::Greater>
{
  using Type = std::greater<double>;
  inline static constexpr auto Name()  //
      -> const char*
  {
    return "Greater";
  }
};
template<>
struct Comparison<ComparisonType::GreaterOrEqual>
{
  using Type = std::greater_equal<double>;
  inline static constexpr auto Name()  //
      -> const char*
  {
    return "GreaterOrEqual";
  }
};
template<>
struct Comparison<ComparisonType::Equal>
{
  using Type = std::equal_to<double>;
  inline static constexpr auto Name()  //
      -> const char*
  {
    return "Equal";
  }
};
template<>
struct Comparison<ComparisonType::NotEqual>
{
  using Type = std::not_equal_to<double>;
  inline static constexpr auto Name()  //
      -> const char*
  {
    return "NotEqual";
  }
};

struct StrategyBase;
template<class _Leaf>
struct UD::Interface::Traits::InheritDescriptor<StrategyBase, _Leaf>
    : public UD::Interface::Traits::Tag<StrategyBase,
                                        UD::Tag::Interface::Modifier::Abstract>
{
};

struct StrategyBase : public UD::Interface::Interface<StrategyBase>
{
  virtual ~StrategyBase()           = default;
  StrategyBase()                    = default;
  StrategyBase(const StrategyBase&) = default;
  StrategyBase(StrategyBase&&)      = default;
  auto operator=(const StrategyBase&) -> StrategyBase& = default;
  auto operator=(StrategyBase&&) -> StrategyBase& = default;

  [[nodiscard]] virtual auto Name() const  //
      -> std::string = 0;
  [[nodiscard]] virtual auto Difference(const History& candles) const  //
      -> std::map<Candle, double, Candle::Compare::OpenTime> = 0;
  [[nodiscard]] virtual auto Compare(const History& candles) const  //
      -> std::map<Candle, bool, Candle::Compare::OpenTime> = 0;
  [[nodiscard]] virtual auto Evaluate(const History& candles) const  //
      -> std::map<Candle, Position, Candle::Compare::OpenTime> = 0;
  [[nodiscard]] virtual auto RunFirst(const History& candles) const  //
      -> std::map<Candle, double, Candle::Compare::OpenTime> = 0;
  [[nodiscard]] virtual auto RunSecond(const History& candles) const  //
      -> std::map<Candle, double, Candle::Compare::OpenTime> = 0;
};
template<class T>
struct StrategyDerived;

template<class _Leaf, class T>
struct UD::Interface::Traits::InheritDescriptor<StrategyDerived<T>, _Leaf>
    : public UD::Interface::Traits::Inherit<_Leaf, StrategyBase>
    , public UD::Interface::Traits::Tag<StrategyDerived<T>,
                                        UD::Tag::Interface::Modifier::Abstract>
{
};
template<class T>
struct StrategyDerived : public UD::Interface::Interface<StrategyDerived<T>>
{
  ~StrategyDerived() override                 = default;
  StrategyDerived()                           = default;
  StrategyDerived(const StrategyDerived&)     = default;
  StrategyDerived(StrategyDerived&&) noexcept = default;
  auto operator=(const StrategyDerived&) -> StrategyDerived& = default;
  auto operator=(StrategyDerived&&) noexcept -> StrategyDerived& = default;

  [[nodiscard]] auto Name() const  //
      -> std::string override
  {
    return T::Name();
  }
  [[nodiscard]] auto Difference(const History& candles) const  //
      -> std::map<Candle, double, Candle::Compare::OpenTime> override
  {
    return T::Difference(candles);
  }
  [[nodiscard]] auto Compare(const History& candles) const  //
      -> std::map<Candle, bool, Candle::Compare::OpenTime> override
  {
    return T::Compare(candles);
  }
  [[nodiscard]] auto Evaluate(const History& candles) const  //
      -> std::map<Candle, Position, Candle::Compare::OpenTime> override
  {
    return T::Evaluate(candles);
  }
  [[nodiscard]] auto RunFirst(const History& candles) const  //
      -> std::map<Candle, double, Candle::Compare::OpenTime> override
  {
    return T::RunFirst(candles);
  }
  [[nodiscard]] auto RunSecond(const History& candles) const  //
      -> std::map<Candle, double, Candle::Compare::OpenTime> override
  {
    return T::RunSecond(candles);
  }
};
template<class Ind1, class Ind2, int _ComparisonType>
// requires Indicator<Ind1> && Indicator<Ind2>
struct Strategy;

template<class _Leaf, class Ind1, class Ind2, int _ComparisonType>
struct UD::Interface::Traits::
    InheritDescriptor<Strategy<Ind1, Ind2, _ComparisonType>,
                                         _Leaf>
    : public UD::Interface::Traits::
          Inherit<_Leaf, StrategyDerived<Strategy<Ind1, Ind2, _ComparisonType>>>
{
};
template<class Ind1, class Ind2, int _ComparisonType = ComparisonType::Less>
// requires Indicator<Ind1> && Indicator<Ind2>
struct Strategy
    : public UD::Interface::Interface<Strategy<Ind1, Ind2, _ComparisonType>>
{
  ~Strategy() override          = default;
  Strategy()                    = default;
  Strategy(const Strategy&)     = default;
  Strategy(Strategy&&) noexcept = default;
  auto operator=(const Strategy&) -> Strategy& = default;
  auto operator=(Strategy&&) noexcept -> Strategy& = default;

  Strategy(Ind1, Ind2) {}
  static auto Name()  //
      -> std::string
  {
    return std::string(Comparison<_ComparisonType>::Name()) + "<" + Ind1::Name()
         + "," + Ind2::Name() + ">";
  }
  static auto Difference(const History& candles)  //
      -> std::map<Candle, double, Candle::Compare::OpenTime>
  {
    std::map<Candle, double, Candle::Compare::OpenTime> ret;
    auto ind1 = Ind1::Run(candles);
    auto ind2 = Ind2::Run(candles);
    for (auto& candle : candles)
    {
      ret[candle] = ind1[candle] - ind2[candle];
    }
    return ret;
  }
  static auto Compare(const History& candles)  //
      -> std::map<Candle, bool, Candle::Compare::OpenTime>
  {
    std::map<Candle, bool, Candle::Compare::OpenTime> ret;
    for (auto& [key, value] : Difference(candles))
    {
      ret[key] = typename Comparison<_ComparisonType>::Type{}(value, 0.0);
    }
    return ret;
  }
  static auto Evaluate(const History& candles)  //
      -> std::map<Candle, Position, Candle::Compare::OpenTime>
  {
    std::map<Candle, Position, Candle::Compare::OpenTime> ret;
    for (auto& [key, value] : Compare(candles))
    {
      ret[key] = value ? Position::Long : Position::Short;
    }
    return ret;
  }
  static auto RunFirst(const History& candles)  //
      -> std::map<Candle, double, Candle::Compare::OpenTime>
  {
    return Ind1::Run(candles);
  }
  static auto RunSecond(const History& candles)  //
      -> std::map<Candle, double, Candle::Compare::OpenTime>
  {
    return Ind2::Run(candles);
  }
};

struct Plan
{
  std::vector<std::unique_ptr<StrategyBase>> _entries = {};
  std::vector<std::unique_ptr<StrategyBase>> _exits   = {};

  virtual ~Plan() = default;
  Plan();
  Plan(const Plan& other);
  Plan(Plan&& other) noexcept;
  auto operator=(const Plan& other)  //
      -> Plan&;
  auto operator=(Plan&& other) noexcept  //
      -> Plan&;
  Plan(const std::vector<std::unique_ptr<StrategyBase>>& entries,
       const std::vector<std::unique_ptr<StrategyBase>>& exits);

  [[nodiscard]] auto entries()  //
      -> std::vector<std::unique_ptr<StrategyBase>>&;
  [[nodiscard]] auto exits()  //
      -> std::vector<std::unique_ptr<StrategyBase>>&;
  [[nodiscard]] auto entries() const  //
      -> const std::vector<std::unique_ptr<StrategyBase>>&;
  [[nodiscard]] auto exits() const  //
      -> const std::vector<std::unique_ptr<StrategyBase>>&;

  auto CloneEntries(const Plan& other)
      -> std::vector<std::unique_ptr<StrategyBase>>;
  auto CloneExits(const Plan& other)
      -> std::vector<std::unique_ptr<StrategyBase>>;
  template<class T>
  void AddEntry()
  {
    _entries.push_back(std::make_unique<T>());
  }
  template<class T>
  void AddExit()
  {
    _exits.push_back(std::make_unique<T>());
  }
  [[nodiscard]] auto Name() const  //
      -> std::string;
  static auto CloneVector(
      const std::vector<std::unique_ptr<StrategyBase>>& vec)  //
      -> std::vector<std::unique_ptr<StrategyBase>>;
  static auto Default()  //
      -> Plan;

 private:
  [[nodiscard]] auto StrategiesRun(
      const std::vector<std::unique_ptr<StrategyBase>>& strategies,
      const History&                                    candles) const  //
      -> std::vector<std::map<Candle, Position, Candle::Compare::OpenTime>>;
  [[nodiscard]] auto EntriesRun(const History& candles) const  //
      -> std::vector<std::map<Candle, Position, Candle::Compare::OpenTime>>;
  [[nodiscard]] auto ExitsRun(const History& candles) const  //
      -> std::vector<std::map<Candle, Position, Candle::Compare::OpenTime>>;
  [[nodiscard]] auto StrategiesRunJoin(
      const std::vector<std::map<Candle, Position, Candle::Compare::OpenTime>>&
          runs) const  //
      -> std::map<Candle, Position, Candle::Compare::OpenTime>;
  [[nodiscard]] auto JoinEntriesExits(
      const std::map<Candle, Position, Candle::Compare::OpenTime>& entries,
      const std::map<Candle, Position, Candle::Compare::OpenTime>& exits)
      const  //
      -> std::map<Candle, Position, Candle::Compare::OpenTime>;

 public:
  [[nodiscard]] auto Run(const History& candles) const  //
      -> std::map<Candle, Position, Candle::Compare::OpenTime>;
  [[nodiscard]] auto Assess(const History& candles,
                            double         loss = 0.01) const  //
      -> std::map<Candle, double, Candle::Compare::OpenTime>;
  [[nodiscard]] auto AssessTotal(const History& candles,
                                 double         loss = 0.01) const  //
      -> std::map<Candle, double, Candle::Compare::OpenTime>;
  [[nodiscard]] auto AssessFinalTotal(const History& candles,
                                      double         loss = 0.01) const  //
      -> double;
  [[nodiscard]] auto AssessWithEntry(double         entry_value,
                                     const History& candles,
                                     double         loss = 0.01) const  //
      -> std::map<Candle, double, Candle::Compare::OpenTime>;
};