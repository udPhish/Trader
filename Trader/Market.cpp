#include "Market.h"

Market::Market() = default;
Market::Market(std::string asset) : asset{std::move(asset)} {
  candles.insert({{Market::TimeFrame::MINUTE, std::set<Candle, Candle::Compare::OpenTime>()},
                  {Market::TimeFrame::HOUR, std::set<Candle, Candle::Compare::OpenTime>()},
                  {Market::TimeFrame::DAY, std::set<Candle, Candle::Compare::OpenTime>()},
                  {Market::TimeFrame::WEEK, std::set<Candle, Candle::Compare::OpenTime>()},
                  {Market::TimeFrame::MONTH, std::set<Candle, Candle::Compare::OpenTime>()}});
}