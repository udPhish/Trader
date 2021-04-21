#pragma once

#include <string>

#include "API.h"
#include "Exchange.h"
#include "Logger.h"

enum class BinanceRequest { PING, TIME, KLINES };
class APIBinance {
  static const std::string HOST;
  static const std::string PORT;
  static const std::string PK;

  API m_api;

 public:
  APIBinance();
  void Connect();
  void Disconnect();

  void Update(Exchange& exchange);
  void Update(Market& market);
  void Update(Market& market, Market::TimeFrame timeframe);
  void UpdateNoConnect(Exchange& exchange);
  void UpdateNoConnect(Market& market);
  void UpdateNoConnect(Market& market, Market::TimeFrame timeframe);
  static std::set<Candle, Candle::Compare::OpenTime> ParseJSON(
      std::string data);

  template <BinanceRequest request_type>
  boost::beast::http::response<boost::beast::http::string_body> Sync();
  template <>
  boost::beast::http::response<boost::beast::http::string_body>
  Sync<BinanceRequest::PING>() {
    return m_api.SyncPost("api/v3/ping", boost::beast::http::verb::get, "", "",
                          {{"X-MBX-APIKEY", PK}});
  }
  template <>
  boost::beast::http::response<boost::beast::http::string_body>
  Sync<BinanceRequest::TIME>() {
    return m_api.SyncPost("api/v3/time", boost::beast::http::verb::get, "", "",
                          {{"X-MBX-APIKEY", PK}});
  }
  template <BinanceRequest request_type>
  boost::beast::http::response<boost::beast::http::string_body> Sync(
      std::string symbol, std::string interval, std::string limit);
  template <>
  boost::beast::http::response<boost::beast::http::string_body>
  Sync<BinanceRequest::KLINES>(std::string symbol, std::string interval,
                               std::string limit) {
    return m_api.SyncPost(
        "api/v3/klines", boost::beast::http::verb::get,
        "symbol=" + symbol + "&interval=" + interval + "&limit=" + limit, "",
        {{"X-MBX-APIKEY", PK}});
  }
  template <BinanceRequest request_type>
  boost::beast::http::response<boost::beast::http::string_body> Sync(
      std::string symbol, Market::TimeFrame timeframe, std::string limit);
  template <>
  boost::beast::http::response<boost::beast::http::string_body>
  Sync<BinanceRequest::KLINES>(std::string symbol, Market::TimeFrame timeframe,
                               std::string limit) {
    std::string timeframe_string;
    switch (timeframe) {
      case Market::TimeFrame::MINUTE:
        timeframe_string = "1m";
        break;
      case Market::TimeFrame::HOUR:
        timeframe_string = "1h";
        break;
      case Market::TimeFrame::DAY:
        timeframe_string = "1d";
        break;
      case Market::TimeFrame::WEEK:
        timeframe_string = "1w";
        break;
      case Market::TimeFrame::MONTH:
        timeframe_string = "1M";
        break;
    }
    if (timeframe_string.size() <= 0)
      Logger::Log<Logger::Level::FatalError>("Timeframe undefined");
    return Sync<BinanceRequest::KLINES>(symbol, timeframe_string, limit);
  }
};
