#pragma once

#define WIN32_LEAN_AND_MEAN
#include <chrono>
#include <set>
#include <sstream>
#include <unordered_map>
#include <array>

#include <wx/wx.h>

#include <wx/file.h>
#include <wx/listbox.h>

#include "boost/beast.hpp"
#include "boost/json.hpp"
#include "openssl/ssl.h"

#include "boost/asio/ssl/stream.hpp"

//#include "openssl/hmac.h"

#include "APIBinance.h"
#include "Candle.h"
#include "ExchangeView.h"
#include "Logger.h"
#include "OpenGLCanvas.h"
#include "PlanView.h"
#include "Strategy.h"
#include "Util.h"

class Main : public wxFrame {
  enum MenuItem { PING, TIME, KLINES, UPDATE_EXCHANGE, _COUNT };
  OpenGLCanvas* m_gl_canvas;
  wxButton* m_query_button;

  Exchange m_exchange;
  ExchangeView* m_exchange_view;

  PlanView* m_plan_view;

  Logger* m_logger;

  std::string GetTimestamp();
  std::string Encrypt(std::string key, std::string data);
  void OnMenuRequest(wxCommandEvent& event);
  void OnSelection(wxCommandEvent& event);
  void OnPlanSelection(wxCommandEvent& event);
  void OnQuery(wxCommandEvent& event);
  void WriteFile(wxString name, wxString contents);
  wxString ReadFile(wxString name);
  void RequestPing();
  void RequestTime();
  void RequestKlines();
  void RequestUpdateExchange(Exchange& exchange);
  void Parse(const wxString& candle_data);
  std::vector<Plan> PlanList();

 public:
  Main();
  ~Main();
};
