#pragma once
#include <array>
#include <chrono>
#include <functional>
#include <map>
#include <set>
#include <vector>

#include "Candle.h"
#include "Strategies.h"
#include "Logger.h"

struct Strategy {
  struct Compare {
    struct Name {
      bool operator()(const Strategy& lhs, const Strategy& rhs) const {
        return lhs.name < rhs.name;
      }
    };
  };
  using CandlePositionMap =
      std::map<Candle, Position, Candle::Compare::OpenTime>;
  using Function = std::function<CandlePositionMap(History)>;

  std::string name;
  Function function;
  Strategy();
  Strategy(std::string name, Function function);
};
struct Plan {
  using CandlePositionMap =
      std::map<Candle, Position, Candle::Compare::OpenTime>;
  using StrategyFunction = std::function<CandlePositionMap(History)>;

  std::vector<Strategy> entries;
  std::vector<Strategy> exits;
  Plan(std::vector<Strategy> entries, std::vector<Strategy> exits);
  CandlePositionMap Run(History candles) const;
  std::map<Candle, double, Candle::Compare::OpenTime> Assess(
      History candles) const;
  double AssessTotal(History candles) const;
  double AssessGain(History candles) const;
  double AssessLoss(History candles) const;
  CandlePositionMap operator()(History candles) const;
  std::string Name() const;
};