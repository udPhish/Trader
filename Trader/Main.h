#pragma once

#define WIN32_LEAN_AND_MEAN
#include <array>
#include <chrono>
#include <set>
#include <sstream>
#include <unordered_map>

#include <wx/wx.h>

#include <wx/file.h>
#include <wx/listbox.h>

#include "openssl/ssl.h"

//#include "openssl/hmac.h"

#include "APIBinance.h"
#include "Candle.h"
#include "CandleView.h"
#include "ExchangeView.h"
#include "Logger.h"
#include "OpenGLCanvas.h"
#include "PlanView.h"
#include "Strategy.h"
#include "Util.h"

class Main : public wxFrame
{
  enum MenuItem
  {
    PING,
    TIME,
    KLINES,
    UPDATE_EXCHANGE,
    _COUNT
  };
  OpenGLCanvas* m_gl_canvas;
  wxButton*     m_query_button;

  Exchange       m_exchange;
  ExchangeView*  m_exchange_view;
  CandleView*    m_candle_view;
  PlanStatsView* m_plan_stats_view;

  PlanView* m_plan_view;

  Logger* m_logger;

  std::string       GetTimestamp();
  std::string       Encrypt(const std::string& key, const std::string& data);
  void              OnMenuRequest(wxCommandEvent& event);
  void              OnSelection(wxCommandEvent& event);
  void              OnPlanSelection(wxCommandEvent& event);
  void              OnQuery(wxCommandEvent& event);
  void              WriteFile(const wxString& name, const wxString& contents);
  wxString          ReadFile(const wxString& name);
  void              RequestPing();
  void              RequestTime();
  void              RequestKlines();
  void              RequestUpdateExchangeNoRefresh(Exchange& exchange);
  void              RequestUpdateExchange(Exchange& exchange);
  void              Parse(const wxString& candle_data);
  std::vector<Plan> PlanList();
  template<class T, class U, class V, class... Ts>
  std::vector<std::unique_ptr<StrategyBase>> CreateAllStrategyCombinations()
  {
    std::vector<std::unique_ptr<StrategyBase>> strats;
    strats.push_back(std::make_unique<Strategy<T, U>>());
    strats.push_back(std::make_unique<Strategy<U, T>>());
    std::vector<std::unique_ptr<StrategyBase>> t_strats
        = CreateAllStrategyCombinations<T, V, Ts...>();
    std::vector<std::unique_ptr<StrategyBase>> u_strats
        = CreateAllStrategyCombinations<U, V, Ts...>();
    for (auto& ptr : t_strats)
      strats.push_back(std::move(ptr));
    for (auto& ptr : u_strats)
      strats.push_back(std::move(ptr));
    return strats;
  }
  template<class T, class U>
  std::vector<std::unique_ptr<StrategyBase>> CreateAllStrategyCombinations()
  {
    std::vector<std::unique_ptr<StrategyBase>> strats;
    strats.push_back(std::make_unique<Strategy<T, U>>());
    strats.push_back(std::make_unique<Strategy<U, T>>());
    return strats;
  }

 public:
  Main();
  ~Main();
};
