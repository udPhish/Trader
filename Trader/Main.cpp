#include "Main.h"

Main::Main() : wxFrame(nullptr, wxID_ANY, "title") {
  m_logger = new Logger(this);

  wxBoxSizer* outer_sizer = new wxBoxSizer(wxVERTICAL);
  wxBoxSizer* inner_sizer = new wxBoxSizer(wxHORIZONTAL);

  wxMenuBar* menu_bar = new wxMenuBar();

  wxMenu* menu_request = new wxMenu();
  menu_request->Append(MenuItem::PING, "Ping", "Ping server");
  menu_request->Append(MenuItem::TIME, "Time", "Get server time");
  menu_request->Append(MenuItem::KLINES, "Klines", "Get candles");
  menu_request->Append(MenuItem::UPDATE_EXCHANGE, "UpdateExchange",
                       "Update entire exchange");

  menu_bar->Append(menu_request, "Request");

  SetMenuBar(menu_bar);

  wxStatusBar* status_bar = new wxStatusBar(this);
  SetStatusBar(status_bar);

  Bind(wxEVT_COMMAND_MENU_SELECTED, &Main::OnMenuRequest, this, MenuItem::PING,
       MenuItem::_COUNT);

  wxGLAttributes gl_attributes;
   gl_attributes.PlatformDefaults().RGBA().DoubleBuffer().Depth(16).EndList();
  //gl_attributes.PlatformDefaults().Defaults().EndList();

  if (!wxGLCanvas::IsDisplaySupported(gl_attributes))
    wxMessageBox("Unsupported");

  SetMinSize(wxSize(800, 600));
  try {
    m_exchange = Load<Exchange>(std::string("exchange.dat"));
  } catch (...) {
    Logger::Log<Logger::Level::Message>(
        "Load Exchange failed, generating new Exchange.");
    try {
      Save(m_exchange, std::string("exchange.dat"));
    } catch (...) {
      Logger::Log<Logger::Level::FatalError>("Exchange generation failed.");
    }
  }
  m_exchange.AddAssets({"BTCUSDT", "ETHUSDT", "ETHBTC"});
  m_exchange_view = new ExchangeView(this, m_exchange);

  m_plan_view = new PlanView(this, PlanList());


  m_gl_canvas = new OpenGLCanvas(this, &m_exchange_view->SelectedCandles(),
                                 &m_plan_view->SelectedPlan(), gl_attributes);
  //m_gl_canvas->ShouldTest(); //TODO: Remove
  outer_sizer->Add(inner_sizer, wxSizerFlags().Expand().Proportion(4));
  outer_sizer->Add(m_logger, wxSizerFlags().Expand());

  inner_sizer->Add(m_exchange_view, wxSizerFlags().Expand().Proportion(0));
  inner_sizer->Add(m_plan_view, wxSizerFlags().Expand().Proportion(0));
  inner_sizer->Add(m_gl_canvas, wxSizerFlags().Expand().Proportion(1));

  SetSizer(outer_sizer);
  SetAutoLayout(true);

  m_exchange_view->Bind(wxEVT_LISTBOX, &Main::OnSelection, this);
  m_exchange_view->Bind(wxEVT_CHOICE, &Main::OnSelection, this);
  m_plan_view->Bind(wxEVT_LISTBOX, &Main::OnPlanSelection, this);

  m_gl_canvas->InitWorld();
}

Main::~Main() {
  try {
    Save(m_exchange, std::string("exchange.dat"));
  } catch (...) {
    Logger::Log<Logger::Level::FatalError>("Unable to save Exchange.");
  }
}

std::string Main::GetTimestamp() {
  return (std::stringstream()
          << std::chrono::duration_cast<std::chrono::milliseconds>(
                 std::chrono::system_clock::now().time_since_epoch())
                 .count())
      .str();
}
std::string Main::Encrypt(std::string key, std::string data) {
  std::string result;
  static char res_hexstring[64];
  int result_len = 32;
  std::string signature;

  // result = std::string(reinterpret_cast<char*>(HMAC(
  //    EVP_sha256(), key.c_str(), strlen((char*)key.c_str()),
  //              const_cast<unsigned char*>(
  //                  reinterpret_cast<const unsigned char*>(data.c_str())),
  //              strlen((char*)data.c_str()), NULL, NULL)));
  for (int i = 0; i < result_len; i++) {
    sprintf(&(res_hexstring[i * 2]), "%02x", result[i]);
  }

  for (int i = 0; i < 64; i++) {
    signature += res_hexstring[i];
  }

  return signature;
}

void Main::OnMenuRequest(wxCommandEvent& event) {
  switch (event.GetId()) {
    case MenuItem::PING:
      RequestPing();
      break;
    case MenuItem::TIME:
      RequestTime();
      break;
    case MenuItem::KLINES:
      RequestKlines();
      break;
    case MenuItem::UPDATE_EXCHANGE:
      RequestUpdateExchange(m_exchange);
      break;
    default:
      break;
  }
}
void Main::OnSelection(wxCommandEvent& event) {
  m_gl_canvas->UpdateCandles(&m_exchange_view->SelectedCandles());
}
void Main::OnPlanSelection(wxCommandEvent& event) {
  m_gl_canvas->UpdatePlan(&m_plan_view->SelectedPlan());
  Logger::Log(wxString("Total: ") << m_plan_view->SelectedPlan().AssessTotal(
                                         m_exchange_view->SelectedCandles())
                                  << ", Gain: "
                                  << m_plan_view->SelectedPlan().AssessGain(
                                         m_exchange_view->SelectedCandles())
                                  << ", Loss: "
                                  << m_plan_view->SelectedPlan().AssessLoss(
                                         m_exchange_view->SelectedCandles()));
}
void Main::OnQuery(wxCommandEvent& event) {}
void Main::WriteFile(wxString name, wxString contents) {
  wxFile query_file = wxFile(name, wxFile::write);
  query_file.Write(contents);
}
wxString Main::ReadFile(wxString name) {
  if (wxFile::Exists(name)) {
    wxFile query_file = wxFile(name, wxFile::read);
    wxString contents;
    query_file.ReadAll(&contents);
    return contents;
  }
  return wxString();
}

void Main::RequestPing() {
  APIBinance api;

  api.Connect();

  boost::beast::http::response<boost::beast::http::string_body> response =
      api.Sync<BinanceRequest::PING>();

  api.Disconnect();

  wxMessageBox(response.body());
}

void Main::RequestTime() {
  APIBinance api;

  api.Connect();

  boost::beast::http::response<boost::beast::http::string_body> response =
      api.Sync<BinanceRequest::TIME>();

  api.Disconnect();

  WriteFile("query.txt", response.body());

  wxMessageBox(response.body());
}
void Main::RequestKlines() {
  APIBinance api;

  api.Connect();

  boost::beast::http::response<boost::beast::http::string_body> response =
      api.Sync<BinanceRequest::KLINES>("BTCUSDT", "1d", "10");

  api.Disconnect();

  wxMessageBox(response.body());
}
void Main::RequestUpdateExchange(Exchange& exchange) {
  APIBinance api;
  api.Update(exchange);
  m_gl_canvas->Refresh();
  m_gl_canvas->Update();
}
std::vector<Plan> Main::PlanList() {
  std::vector<Strategy> strategies;
  strategies.push_back({"IsAboveSMA<2>", IsAboveSMA<2>});
  strategies.push_back({"IsAboveSMA<4>", IsAboveSMA<4>});
  strategies.push_back({"IsAboveSMA<8>", IsAboveSMA<8>});
  strategies.push_back({"IsAboveSMA<16>", IsAboveSMA<16>});

  std::vector<std::vector<Strategy>> strategy_combinations =
      Combinations(strategies, 2);

  std::vector<Plan> plans;
  for (auto& entry_strategies : strategy_combinations) {
    for (auto& exit_strategies : strategy_combinations) {
      plans.push_back({entry_strategies, exit_strategies});
    }
  }
  return plans;
}