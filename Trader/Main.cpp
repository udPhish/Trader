#include "Main.h"

Main::Main()
    : wxFrame(nullptr, wxID_ANY, "title")
    , m_query_button{nullptr}
    , m_plan_view{nullptr}
{
  m_logger          = new Logger(this);
  m_candle_view     = new CandleView(this);
  m_plan_stats_view = new PlanStatsView(this);

  auto* outer_sizer     = new wxBoxSizer(wxVERTICAL);
  auto* inner_sizer     = new wxBoxSizer(wxHORIZONTAL);
  auto* log_sizer       = new wxBoxSizer(wxHORIZONTAL);
  auto* selection_sizer = new wxBoxSizer(wxHORIZONTAL);

  auto* menu_bar = new wxMenuBar();

  auto* menu_request = new wxMenu();
  menu_request->Append(MenuItem::PING, "Ping", "Ping server");
  menu_request->Append(MenuItem::TIME, "Time", "Get server time");
  menu_request->Append(MenuItem::KLINES, "Klines", "Get candles");
  menu_request->Append(MenuItem::UPDATE_EXCHANGE,
                       "UpdateExchange",
                       "Update entire exchange");

  menu_bar->Append(menu_request, "Request");

  SetMenuBar(menu_bar);

  auto* status_bar = new wxStatusBar(this);
  SetStatusBar(status_bar);

  Bind(wxEVT_COMMAND_MENU_SELECTED,
       &Main::OnMenuRequest,
       this,
       MenuItem::PING,
       MenuItem::_COUNT);

  wxGLAttributes gl_attributes;
  gl_attributes.PlatformDefaults().RGBA().DoubleBuffer().Depth(16).EndList();
  // gl_attributes.PlatformDefaults().Defaults().EndList();

  if (!wxGLCanvas::IsDisplaySupported(gl_attributes))
  {
    wxMessageBox("Unsupported");
  }

  SetMinSize(wxSize(800, 600));
  try
  {
    m_exchange = Load<Exchange>(std::string("exchange.dat"));
  }
  catch (...)
  {
    Logger::Log<Logger::Level::Message>(
        "Load Exchange failed, generating new Exchange.");
    try
    {
      Save(m_exchange, std::string("exchange.dat"));
    }
    catch (...)
    {
      Logger::Log<Logger::Level::FatalError>("Exchange generation failed.");
    }
  }
  if (m_exchange.markets.empty())
  {
    m_exchange.AddAssets({"BTCUSDT", "ETHUSDT", "ETHBTC"});
    RequestUpdateExchangeNoRefresh(m_exchange);
  }
  m_exchange_view = new ExchangeView(this, m_exchange);

  m_plan_view
      = new PlanView(this, m_exchange_view->SelectedCandles(), PlanList());

  m_gl_canvas = new OpenGLCanvas(this,
                                 m_plan_stats_view,
                                 m_candle_view,
                                 &m_exchange_view->SelectedCandles(),
                                 gl_attributes);
  m_gl_canvas->UpdatePlan(m_plan_view->SelectedPlan());
  m_gl_canvas->UpdateStrategy(m_plan_view->SelectedStrategy());
  m_gl_canvas->InitWorld();

  outer_sizer->Add(inner_sizer, wxSizerFlags().Expand().Proportion(4));
  outer_sizer->Add(log_sizer, wxSizerFlags().Expand().Proportion(0));

  log_sizer->Add(m_candle_view, wxSizerFlags().Expand().Proportion(1));
  log_sizer->Add(m_plan_stats_view, wxSizerFlags().Expand().Proportion(1));
  log_sizer->Add(m_logger, wxSizerFlags().Expand().Proportion(2));

  inner_sizer->Add(selection_sizer, wxSizerFlags().Expand().Proportion(2));
  selection_sizer->Add(m_exchange_view, wxSizerFlags().Expand().Proportion(0));
  selection_sizer->Add(m_plan_view, wxSizerFlags().Expand().Proportion(1));
  inner_sizer->Add(m_gl_canvas, wxSizerFlags().Expand().Proportion(2));

  SetSizer(outer_sizer);
  SetAutoLayout(true);

  m_exchange_view->Bind(wxEVT_LISTBOX, &Main::OnSelection, this);
  m_exchange_view->Bind(wxEVT_CHOICE, &Main::OnSelection, this);
  m_plan_view->Bind(wxEVT_LISTBOX, &Main::OnPlanSelection, this);
}

Main::~Main()
{
  try
  {
    Save(m_exchange, std::string("exchange.dat"));
  }
  catch (...)
  {
    Logger::Log<Logger::Level::FatalError>("Unable to save Exchange.");
  }
}

auto Main::GetTimestamp() -> std::string
{
  return (std::stringstream()
          << std::chrono::duration_cast<std::chrono::milliseconds>(
                 std::chrono::system_clock::now().time_since_epoch())
                 .count())
      .str();
}
auto Main::Encrypt(const std::string& key, const std::string& data)
    -> std::string
{
  std::string result;
  static char res_hexstring[64];
  int         result_len = 32;
  std::string signature;

  // result = std::string(reinterpret_cast<char*>(HMAC(
  //    EVP_sha256(), key.c_str(), strlen((char*)key.c_str()),
  //              const_cast<unsigned char*>(
  //                  reinterpret_cast<const unsigned char*>(data.c_str())),
  //              strlen((char*)data.c_str()), NULL, NULL)));
  for (int i = 0; i < result_len; i++)
  {
    sprintf(&(res_hexstring[i * 2]), "%02x", result[i]);
  }

  for (int i = 0; i < 64; i++)
  {
    signature += res_hexstring[i];
  }

  return signature;
}

void Main::OnMenuRequest(wxCommandEvent& event)
{
  switch (event.GetId())
  {
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
void Main::OnSelection(wxCommandEvent& event)
{
  auto& candles = m_exchange_view->SelectedCandles();
  m_gl_canvas->UpdateCandles(&candles);
  m_plan_view->Reset(candles, PlanList());
  m_gl_canvas->InitWorld();
}
void Main::OnPlanSelection(wxCommandEvent& event)
{
  m_gl_canvas->UpdatePlan(m_plan_view->SelectedPlan());
  m_gl_canvas->UpdateStrategy(m_plan_view->SelectedStrategy());
  m_gl_canvas->InitWorldSaveView();
  // Logger::Log(wxString("Total: ") <<
  // m_plan_view->SelectedPlan().AssessTotal(m_exchange_view->SelectedCandles())
  //                                << ", Gain: " <<
  //                                m_plan_view->SelectedPlan().AssessGain(m_exchange_view->SelectedCandles())
  //                                << ", Loss: " <<
  //                                m_plan_view->SelectedPlan().AssessLoss(m_exchange_view->SelectedCandles()));
}
void Main::OnQuery(wxCommandEvent& event) {}
void Main::WriteFile(const wxString& name, const wxString& contents)
{
  wxFile query_file = wxFile(name, wxFile::write);
  query_file.Write(contents);
}
auto Main::ReadFile(const wxString& name) -> wxString
{
  if (wxFile::Exists(name))
  {
    wxFile   query_file = wxFile(name, wxFile::read);
    wxString contents;
    query_file.ReadAll(&contents);
    return contents;
  }
  return wxString();
}

void Main::RequestPing()
{
  APIBinance api;

  api.Connect();

  auto response = api.Sync<BinanceRequest::PING>();

  api.Disconnect();

  wxMessageBox(response.body());
}

void Main::RequestTime()
{
  APIBinance api;

  api.Connect();

  auto response = api.Sync<BinanceRequest::TIME>();

  api.Disconnect();

  WriteFile("query.txt", response.body());

  wxMessageBox(response.body());
}
void Main::RequestKlines()
{
  APIBinance api;

  api.Connect();

  auto response = api.Sync<BinanceRequest::KLINES>("BTCUSDT", "1d", "10");

  api.Disconnect();

  wxMessageBox(response.body());
}
void Main::RequestUpdateExchangeNoRefresh(Exchange& exchange)
{
  APIBinance api;
  api.Update(exchange);
}
void Main::RequestUpdateExchange(Exchange& exchange)
{
  RequestUpdateExchangeNoRefresh(exchange);
  m_gl_canvas->Refresh();
  m_gl_canvas->Update();
}
auto Main::PlanList() -> std::vector<Plan>
{
  auto strategies = std::vector<std::unique_ptr<StrategyBase>>{
      CreateAllStrategyCombinations<SimpleMovingAverage<8, Metric::Close>,
                                    Identity<Metric::Close>,
                                    SimpleMovingAverage<8, Metric::Open>,
                                    Identity<Metric::Open>>()};
  auto strategy_combinations = Combinations(strategies, 2);

  auto plans = std::vector<Plan>{};
  for (auto& entry_strategies : strategy_combinations)
  {
    for (auto& exit_strategies : strategy_combinations)
    {
      plans.push_back({entry_strategies, exit_strategies});
    }
  }
  auto strategies2 = std::vector<std::unique_ptr<StrategyBase>>{
      CreateAllStrategyCombinations<SimpleMovingAverage<8, Metric::Volume>,
                                    Identity<Metric::Volume>>()};
  auto strategy_combinations2 = Combinations(strategies2, 2);

  for (auto& entry_strategies : strategy_combinations2)
  {
    for (auto& exit_strategies : strategy_combinations2)
    {
      plans.push_back({entry_strategies, exit_strategies});
    }
  }
  return plans;
}