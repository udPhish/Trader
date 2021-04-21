#include "APIBinance.h"

void APIBinance::Connect() { m_api.Connect(); }
void APIBinance::Disconnect() { m_api.Disconnect(); }
APIBinance::APIBinance() : m_api(HOST, PORT) {}

void APIBinance::Update(Exchange& exchange) {
  Logger::Log<Logger::Level::Message>(wxString("Updating ")
                                      << exchange.markets.size() << " markets");
  Connect();
  UpdateNoConnect(exchange);
  Disconnect();
}
void APIBinance::Update(Market& market) {
  Logger::Log<Logger::Level::Message>(wxString("Updating ") << market.asset);
  Connect();
  UpdateNoConnect(market);
  Disconnect();
}
void APIBinance::Update(Market& market, Market::TimeFrame timeframe) {
  Logger::Log<Logger::Level::Message>(wxString("Updating(timeframe) ")
                                      << market.asset);
  Connect();
  UpdateNoConnect(market, timeframe);
  Disconnect();
}
void APIBinance::UpdateNoConnect(Exchange& exchange) {
  for (auto& market : exchange.markets) {
    Logger::Log<Logger::Level::Message>(wxString() << market.first);
    UpdateNoConnect(market.second);
  }
}
void APIBinance::UpdateNoConnect(Market& market) {
  for (auto& timeframe : market.candles) {
    UpdateNoConnect(market, timeframe.first);
  }
}
void APIBinance::UpdateNoConnect(Market& market, Market::TimeFrame timeframe) {
  auto response = Sync<BinanceRequest::KLINES>(market.asset, timeframe, "1000");
  if (response.result_int() == 429) {
    Logger::Log(wxString("Wait request received (429)... Waiting for ")
                << response.base().at("Retry-After").to_string()
                << " seconds.");

    Logger::Log<Logger::Level::FatalError>("Wait request received (429)");
  }
  Logger::Log<Logger::Level::Status>(
      wxString() << "Weight: "
                 << response.base().at("X-MBX-USED-WEIGHT").to_string());
  market.candles[timeframe] = ParseJSON(response.body());
}

std::set<Candle, Candle::Compare::OpenTime> APIBinance::ParseJSON(
    std::string data) {
  boost::json::array json_candles =
      boost::json::parse(std::string(data)).as_array();

  std::set<Candle, Candle::Compare::OpenTime> candles;
  for (auto& json_candle : json_candles) {
    auto open_time =
        std::chrono::duration_cast<std::chrono::duration<int, std::ratio<60>>>(
            std::chrono::milliseconds(json_candle.at(0).as_int64()));
    double open = wxAtof(json_candle.at(1).as_string().c_str());
    double high = wxAtof(json_candle.at(2).as_string().c_str());
    double low = wxAtof(json_candle.at(3).as_string().c_str());
    double close = wxAtof(json_candle.at(4).as_string().c_str());
    double volume = wxAtof(json_candle.at(5).as_string().c_str());
    auto close_time =
        std::chrono::duration_cast<std::chrono::duration<int, std::ratio<60>>>(
            std::chrono::milliseconds(json_candle.at(6).as_int64()));
    double quote_asset_volume = wxAtof(json_candle.at(7).as_string().c_str());
    std::size_t number_of_trades = json_candle.at(8).as_int64();
    double taker_buy_base_asset_volume =
        wxAtof(json_candle.at(9).as_string().c_str());
    double taker_buy_quote_asset_volume =
        wxAtof(json_candle.at(10).as_string().c_str());
    double ignore = wxAtof(json_candle.at(11).as_string().c_str());
    candles.emplace(Candle{open_time, open, high, low, close, volume,
                           close_time, quote_asset_volume, number_of_trades,
                           taker_buy_base_asset_volume,
                           taker_buy_quote_asset_volume, ignore});
  }
  return candles;
}
const std::string APIBinance::HOST = "api.binance.com";
const std::string APIBinance::PORT = "https";
const std::string APIBinance::PK =
    "bGmYjY5Qx47BBQODtyZ7jl4B6JtL2hogRW6NW3dlU1DHxMe4UyvENhMnVy4sXyS1";