#include "Strategy.h"

Strategy::Strategy() {}
Strategy::Strategy(std::string name, Function function)
    : name{name}, function{function} {}

Plan::Plan(std::vector<Strategy> entries, std::vector<Strategy> exits)
    : entries{entries}, exits{exits} {}
Plan::CandlePositionMap Plan::Run(History candles) const {
  std::vector<CandlePositionMap> entry_results;
  std::vector<CandlePositionMap> exit_results;

  for (auto& entry : entries) entry_results.push_back(entry.function(candles));
  for (auto& exit : exits) exit_results.push_back(exit.function(candles));

  CandlePositionMap ret;
  Position pos = Position::Short;
  for (auto& candle : candles) {
    if (pos == Position::Short) {
      pos = Position::Long;
      for (auto& entry_result : entry_results) {
        if (entry_result.at(candle) == Position::Short) {
          pos = Position::Short;
          break;
        }
      }
    } else if (pos == Position::Long) {
      pos = Position::Short;
      for (auto& exit_result : exit_results) {
        if (exit_result.at(candle) == Position::Long) {
          pos = Position::Long;
          break;
        }
      }
    }
    ret[candle] = pos;
  }
  return ret;
}
Plan::CandlePositionMap Plan::operator()(History candles) const {
  return Run(candles);
}
std::string Plan::Name() const {
  std::string n;
  for (auto& entry : entries) n += entry.name;
  n += "->";
  for (auto& exit : exits) n += exit.name;
  return n;
}

std::map<Candle, double, Candle::Compare::OpenTime> Plan::Assess(
    History candles) const {
  CandlePositionMap positions = Run(candles);
  std::map<Candle, double, Candle::Compare::OpenTime> ret;
  Position pos = Position::Short;
  for (auto [it, entry] = std::pair{candles.begin(), candles.begin()}; it != candles.end();
       ++it) {
    if (pos == Position::Short) {
      if (positions.at(*it) == Position::Long) {
        entry = it;
        pos = Position::Long;
      }
    } else if (pos == Position::Long) {
      if (positions.at(*it) == Position::Short) {
        ret[*entry] = (it->close - entry->close);
        pos = Position::Short;
      }
    }
  }
  return ret;
}
double Plan::AssessTotal(History candles) const {
  std::map<Candle, double, Candle::Compare::OpenTime> assessment =
      Assess(candles);
  double ret = 0;
  for (auto& p : assessment) {
    ret += p.second;
  }
  return ret;
}
double Plan::AssessGain(History candles) const {
  std::map<Candle, double, Candle::Compare::OpenTime> assessment =
      Assess(candles);
  double ret = 0;
  for (auto& p : assessment) {
    if (p.second > 0) {
      ret += p.second;
    }
  }
  return ret;
}
double Plan::AssessLoss(History candles) const {
  std::map<Candle, double, Candle::Compare::OpenTime> assessment =
      Assess(candles);
  double ret = 0;
  for (auto& p : assessment) {
    if (p.second < 0) {
      ret += p.second;
    }
  }
  return ret;
}