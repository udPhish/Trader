#pragma once

#include <map>

#include "Candle.h"

template <class T>
concept Indicator = requires(T t, History hist) {
  { T::Run(hist) }
  ->std::same_as<std::map<Candle, double, Candle::Compare::OpenTime>>;
  { T::Name() }
  ->std::same_as<std::string>;
};
//struct Indicator {
//  virtual std::map<Candle, double, Candle::Compare::OpenTime> operator()(const History& candles) = 0;
//  virtual std::string name() = 0;
//};
