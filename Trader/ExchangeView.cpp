#include "ExchangeView.h"

ExchangeView::ExchangeView(wxWindow* parent, Exchange& exchange)
    : wxWindow(parent, wxID_ANY),
      m_list_markets{new wxListBox{this, wxID_ANY, wxDefaultPosition, wxDefaultSize, 0, nullptr, wxLB_SINGLE | wxLB_SORT}},
      m_choice_timeframe{new wxChoice{this, wxID_ANY}},
      m_exchange{exchange} {
  auto* outer = new wxBoxSizer(wxVERTICAL);

  m_choice_timeframe->Append({"Minute", "Hour", "Day", "Week", "Month"});
  m_choice_timeframe->Select(0);

  if (m_exchange.markets.empty()) {
    Logger::Log<Logger::Level::FatalError>("No markets provided");
  }
  for (auto& it : m_exchange.markets) {
    m_list_markets->Append(it.first);
  }
  m_list_markets->Select(0);

  outer->Add(m_choice_timeframe, wxSizerFlags().Expand());
  outer->Add(m_list_markets, wxSizerFlags().Expand().Proportion(1));
  SetSizer(outer);
  SetAutoLayout(true);
}

auto ExchangeView::SelectedCandles() -> History& {
  wxString market_selection = m_list_markets->GetStringSelection();
  wxString timeframe_selection = m_choice_timeframe->GetStringSelection();
  if (market_selection == "" || timeframe_selection == "") {
    Logger::Log<Logger::Level::FatalError>("Invalid market selection");
  }
  Market::TimeFrame timeframe;
  if (timeframe_selection == "Minute") {
    timeframe = Market::TimeFrame::MINUTE;
  } else if (timeframe_selection == "Hour") {
    timeframe = Market::TimeFrame::HOUR;
  } else if (timeframe_selection == "Day") {
    timeframe = Market::TimeFrame::DAY;
  } else if (timeframe_selection == "Week") {
    timeframe = Market::TimeFrame::WEEK;
  } else if (timeframe_selection == "Month") {
    timeframe = Market::TimeFrame::MONTH;
  } else {
    timeframe = Market::TimeFrame::DAY;
  }

  return m_exchange.markets.at(market_selection.ToStdString()).candles.at(timeframe);
}