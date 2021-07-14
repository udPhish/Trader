#pragma once

#include <chrono>
#include <map>
#include <set>
#include <vector>

#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/serialization/set.hpp>

#include "Boost_Serialization_chrono.h"

struct Candle
{
  struct Compare
  {
    struct OpenTime
    {
      bool operator()(const Candle& lhs, const Candle& rhs) const
      {
        return lhs.open_time < rhs.open_time;
      }
    };
    struct CloseTime
    {
      bool operator()(const Candle& lhs, const Candle& rhs) const
      {
        return lhs.close_time < rhs.close_time;
      }
    };
    struct Open
    {
      bool operator()(const Candle& lhs, const Candle& rhs) const
      {
        return lhs.open < rhs.open;
      }
    };
    struct High
    {
      bool operator()(const Candle& lhs, const Candle& rhs) const
      {
        return lhs.high < rhs.high;
      }
    };
    struct Low
    {
      bool operator()(const Candle& lhs, const Candle& rhs) const
      {
        return lhs.low < rhs.low;
      }
    };
    struct Close
    {
      bool operator()(const Candle& lhs, const Candle& rhs) const
      {
        return lhs.close < rhs.close;
      }
    };
    struct Volume
    {
      bool operator()(const Candle& lhs, const Candle& rhs) const
      {
        return lhs.volume < rhs.volume;
      }
    };
  };

  std::chrono::milliseconds open_time;
  double                    open;
  double                    high;
  double                    low;
  double                    close;
  double                    volume;
  std::chrono::milliseconds close_time;
  double                    quote_asset_volume;
  std::size_t               number_of_trades;
  double                    taker_buy_base_asset_volume;
  double                    taker_buy_quote_asset_volume;
  double                    ignore;
  // Candle(){}
  // Candle(std::chrono::milliseconds open_time, long double open, long double high,
  //       long double low, long double close, long double volume,
  //       std::chrono::milliseconds close_time, long double quote_asset_volume,
  //       std::size_t number_of_trades, long double taker_buy_base_asset_volume,
  //       long double taker_buy_quote_asset_volume, long double ignore)
  //    : open_time{open_time},
  //      open{open},
  //      high{high},
  //      low{low},
  //      close{close},
  //      volume{volume},
  //      close_time{close_time},
  //      quote_asset_volume{quote_asset_volume},
  //      number_of_trades{number_of_trades},
  //      taker_buy_base_asset_volume{taker_buy_base_asset_volume},
  //      taker_buy_quote_asset_volume{taker_buy_quote_asset_volume},
  //      ignore{ignore} {}

  friend class boost::serialization::access;
  template<class Archive>
  void serialize(Archive& ar, const unsigned int version)
  {
    ar& open_time;
    ar& open;
    ar& high;
    ar& low;
    ar& close;
    ar& volume;
    ar& close_time;
    ar& quote_asset_volume;
    ar& number_of_trades;
    ar& taker_buy_base_asset_volume;
    ar& taker_buy_quote_asset_volume;
    ar& ignore;
  }
};
class History
{
  using CandleSet = std::set<Candle, Candle::Compare::OpenTime>;
  CandleSet m_candles;

 public:
  History() noexcept {}
  History(CandleSet candles)
      : m_candles{candles} {}

  auto begin()
  {
    return m_candles.begin();
  }
  auto rbegin()
  {
    return m_candles.rbegin();
  }
  auto end()
  {
    return m_candles.end();
  }
  auto rend()
  {
    return m_candles.rend();
  }

  auto begin() const
  {
    return m_candles.begin();
  }
  auto rbegin() const
  {
    return m_candles.rbegin();
  }
  auto end() const
  {
    return m_candles.end();
  }
  auto rend() const
  {
    return m_candles.rend();
  }

  auto size() const
  {
    return m_candles.size();
  }
  auto empty() const
  {
    return m_candles.empty();
  }

  template<class T>
  std::map<Candle, T, Candle::Compare::OpenTime> Extract(T Candle::*member) const
  {
    std::map<Candle, T, Candle::Compare::OpenTime> ret;
    for (auto& candle : m_candles)
      ret.insert({candle, candle.*member});
    return ret;
  }
  auto ExtractOpenTime() const
  {
    return Extract(&Candle::open_time);
  }
  auto ExtractCloseTime() const
  {
    return Extract(&Candle::close_time);
  }
  auto ExtractOpen() const
  {
    return Extract(&Candle::open);
  }
  auto ExtractClose() const
  {
    return Extract(&Candle::close);
  }
  auto ExtractHigh() const
  {
    return Extract(&Candle::high);
  }
  auto ExtractLow() const
  {
    return Extract(&Candle::low);
  }
  auto ExtractVolume() const
  {
    return Extract(&Candle::volume);
  }

  friend class boost::serialization::access;
  template<class Archive>
  void serialize(Archive& ar, const unsigned int version)
  {
    ar& m_candles;
  }
};