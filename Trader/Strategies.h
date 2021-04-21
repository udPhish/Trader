#pragma once

#include <map>
#include <vector>
#include <chrono>

#include "Candle.h"
enum class Position { Long, Short };

// class TrueRange{
//  double Eval(std::array<Candle> candles) override {
//    if (candles.size() < 1) throw;
//    if (candles.size() == 1) return candles[0].high - candles[0].low;
//    return std::max(candles[0].high - candles[0].low,
//                    std::abs(candles[0].high - candles[1].close),
//                    std::abs(candles[0].low - candles[1].close));
//  }
//};

template <std::size_t length>
std::vector<double> SimpleMovingAverage(std::vector<double> values) {
  for (std::size_t i = 0; i < values.size(); ++i) {
    double sum = 0;
    std::size_t end = std::min(values.size(), i + length);
    for (std::size_t j = i; j < end; ++j) {
      sum += values[j];
    }
    values[i] = sum / (end - i);
  }
  return values;
}
std::vector<double> TrueRange(History candles);
template <std::size_t length>
std::vector<double> AverageTrueRange(History candles) {
  return SimpleMovingAverage<length>(TrueRange(candles));
}

template <std::size_t length>
std::map<Candle, Position, Candle::Compare::OpenTime> IsAboveAverageTrueRange(
    History candles) {
}
template <std::size_t length>
std::map<Candle, Position, Candle::Compare::OpenTime> IsAboveSMA(
    History candles) {
  auto closes{candles.ExtractClose()};
  auto sma = SimpleMovingAverage<length>(closes);
  std::map<Candle, Position, Candle::Compare::OpenTime> ret;
  for (auto [i, it] = std::pair(std::size_t(0), candles.begin());
       i < closes.size(); ++i, ++it) {
    if (closes[i] > sma[i])
      ret[*it] = Position::Short;
    else
      ret[*it] = Position::Long;
  }
  return ret;
}

// class ATR : public Indicator<double, length> {
//  double Eval(std::vector<Candle> candles) override {
//    for (auto& candle : candles) {
//    }
//  }
//};
// std::map<std::chrono::milliseconds, Strategy::Position> Strategy::Run(
//    std::set<Candle> candles) {
//  std::map<std::chrono::milliseconds, Strategy::Position> strat;
//
//  for (auto it = candles.begin(); it != candles.end(); ++it) {
//    double ave = it->close;
//    std::size_t i = 1;
//    if (it != candles.begin()) {
//      for (; i < 10; ++i) {
//        ave += std::prev(it, i)->close;
//        if (std::prev(it, i) == candles.begin()) break;
//      }
//    }
//    ave /= i;
//    if (ave < it->close)
//      strat[it->open_time] = Position::Long;
//    else
//      strat[it->open_time] = Position::Short;
//  }
//  return strat;
//}