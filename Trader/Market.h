#pragma once
#include <chrono>
#include <set>
#include <string>
#include <string_view>
#include <unordered_map>

#include "Candle.h"

struct Market {
  enum class TimeFrame { MINUTE, HOUR, DAY, WEEK, MONTH };

  std::string asset;
  std::unordered_map<TimeFrame, History> candles;

  Market();
  Market(std::string asset);

  friend class boost::serialization::access;
  template <class Archive>
  void serialize(Archive& ar, const unsigned int version) {
    ar& asset;
    ar& candles;
  }
};
