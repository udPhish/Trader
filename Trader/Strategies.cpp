#include "Strategies.h"
auto TrueRange(History candles) -> std::vector<double> {
  std::vector<Candle> vec_candles{candles.begin(), candles.end()};
  std::vector<double> ret;
  if (candles.size() < 1) {
    return ret;
  }
  for (std::size_t i = 0; i < vec_candles.size() - 1; ++i) {
    ret.push_back(std::max(vec_candles[i].high - vec_candles[i].low, std::max(std::abs(vec_candles[i].high - vec_candles[i + 1].close),
                                                                              std::abs(vec_candles[i].low - vec_candles[i + 1].close))));
  }
  ret.push_back(vec_candles[vec_candles.size() - 1].high - vec_candles[vec_candles.size() - 1].low);
  return ret;
}