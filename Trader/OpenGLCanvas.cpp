#include "OpenGLCanvas.h"

#include <iostream>

OpenGLCanvas::OpenGLCanvas(wxFrame*              parent,
                           PlanStatsView*        plan_stats_view,
                           CandleView*           candle_view,
                           const History*        candles,
                           const wxGLAttributes& attributes)
    : m_candles(candles)
    , m_context{std::make_unique<wxGLContext>(this)}
    , m_ind1_id{0}
    , m_ind2_id{0}
    , m_marker{}
    , m_marker_id{0}
    , m_candle_view{candle_view}
    , m_plan_stats_view{plan_stats_view}
    , m_plan()
    , m_strategy{nullptr}
    , wxGLCanvas{parent,
                 attributes,
                 wxID_ANY,
                 wxDefaultPosition,
                 wxDefaultSize,
                 wxFULL_REPAINT_ON_RESIZE}
    , m_world{}
    , m_dragging{false}
    , m_old_x{0}
    , m_view{&m_world.CreateView(
          UD::World::World<double, 2>::Rectangle{{0.0, 0.0}, {0.0, 0.0}})}
    , m_strat_view{&m_world.CreateView(
          UD::World::World<double, 2>::Rectangle{{0.0, 0.0}, {0.0, 0.0}})}
    , m_strat_id{0}
    , m_mouse_x_id{0}
    , m_mouse_y_id{0}
    , m_selected_entity_id{0}
    , m_text_marker_pos{new wxStaticText(this,
                                         wxID_ANY,
                                         wxT(""),
                                         wxDefaultPosition,
                                         wxDefaultSize,
                                         wxTRANSPARENT_WINDOW)}
{
  Bind(wxEVT_SIZE, &OpenGLCanvas::OnSize, this);
  Bind(wxEVT_PAINT, &OpenGLCanvas::OnPaint, this);

  Bind(wxEVT_MOUSEWHEEL, &OpenGLCanvas::OnMouseWheel, this);
  Bind(wxEVT_KEY_DOWN, &OpenGLCanvas::OnKeyDown, this);
  Bind(wxEVT_MOTION, &OpenGLCanvas::OnMouseMotion, this);
  Bind(wxEVT_LEAVE_WINDOW, &OpenGLCanvas::OnMouseLeaveWindow, this);
  Bind(wxEVT_RIGHT_UP, &OpenGLCanvas::OnMouseRightUp, this);
  Bind(wxEVT_LEFT_UP, &OpenGLCanvas::OnMouseLeftUp, this);
  Bind(wxEVT_LEFT_DOWN, &OpenGLCanvas::OnMouseLeftDown, this);

  SetBackgroundStyle(wxBG_STYLE_CUSTOM);
}
OpenGLCanvas::OpenGLCanvas(wxFrame*              parent,
                           PlanStatsView*        plan_stats_view,
                           CandleView*           candle_view,
                           const wxGLAttributes& attributes)
    : OpenGLCanvas{parent, plan_stats_view, candle_view, nullptr, attributes}
{
}

void OpenGLCanvas::UpdateCandles(const History* candles)
{
  m_candles = candles;
}
void OpenGLCanvas::UpdatePlan(Plan plan)
{
  m_plan = std::move(plan);
}
void OpenGLCanvas::UpdateStrategy(std::unique_ptr<StrategyBase> strategy)
{
  m_strategy = std::move(strategy);
}
void OpenGLCanvas::PrepareViewport()
{
  PrepareGL();
  wxSize frame = GetSize() * GetContentScaleFactor();
  glViewport(0, 0, frame.GetX() / 2, frame.GetY() / 2);
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  glOrtho(0, frame.GetWidth(), frame.GetHeight(), 0, 0, 1);
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();
}
void OpenGLCanvas::PrepareGL()
{
  glEnable(GL_TEXTURE_2D);
  glEnable(GL_COLOR_MATERIAL);
  glEnable(GL_BLEND);
  glDisable(GL_DEPTH_TEST);
  glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
}

void OpenGLCanvas::OnSize(wxSizeEvent& event)
{
  // event.Skip();
  // wxGLCanvas::SetCurrent(*m_context);

  // wxSize frame = GetSize();  //*GetContentScaleFactor();
  // glViewport(0, 0, frame.GetX()/2, frame.GetY()/2);
  Refresh();
  // event.Skip();
}
auto OpenGLCanvas::CreateEntityCandles(const History& candles,
                                       const Plan&    plan)  //
    -> std::vector<UD::World::World<double, 2>::ID>
{
  auto plan_results = plan.Run(candles);
  auto entities     = std::vector<UD::World::World<double, 2>::ID>{};
  for (auto& candle : candles)
  {
    entities.push_back(m_world.Add<EntityCandle>(candle, plan_results[candle]));
  }
  return entities;
}
auto OpenGLCanvas::MultiplyAllPreviousElements(
    const std::map<Candle, double, Candle::Compare::OpenTime>& candles)  //
    -> std::map<Candle, double, Candle::Compare::OpenTime>
{
  auto ret = std::map<Candle, double, Candle::Compare::OpenTime>{};

  double sum = 1.0;
  for (auto& [key, val] : candles)
  {
    sum *= val;
    ret[key] = sum;
  }
  return ret;
}
auto OpenGLCanvas::AccumulateAllPreviousElements(
    const std::map<Candle, double, Candle::Compare::OpenTime>& candles)  //
    -> std::map<Candle, double, Candle::Compare::OpenTime>
{
  auto ret = std::map<Candle, double, Candle::Compare::OpenTime>{};

  double sum = 0;
  for (auto& [key, val] : candles)
  {
    sum += val;
    ret[key] = sum;
  }
  return ret;
}
auto OpenGLCanvas::InitView(
    const std::vector<UD::World::World<double, 2>::ID>& entities)  //
    -> UD::World::World<double, 2>::View*
{
  m_world.ClearViews();
  std::size_t num_entities_to_view = 50;
  auto*       view                 = &m_world.CreateView(
      UD::World::World<double, 2>::Rectangle{{0.0, 0.0}, {1.0, 1.0}});
  view->FitY(entities);
  if (entities.size() < num_entities_to_view)
  {
    view->FitX(entities);
  }
  else
  {
    std::vector<UD::World::World<double, 2>::ID> es;
    for (std::size_t i = entities.size() - num_entities_to_view;
         i < entities.size();
         ++i)
      es.push_back(entities[i]);
    view->FitX(es);
  }
  return view;
}
auto OpenGLCanvas::CreateMarkerMesh()  //
    -> UD::World::VerticalLine<double, 2>
{
  using Point = UD::World::World<double, 2>::Point;
  Point a{m_marker};
  a.colour() = UD::Colour::White();
  return {a, m_view->rect().height()};
}
auto OpenGLCanvas::CreateMarker()  //
    -> UD::World::World<double, 2>::ID
{
  return m_world.Add<UD::World::World<double, 2>::Entity>(CreateMarkerMesh());
}
auto OpenGLCanvas::CreateStratEntity()  //
    -> UD::World::World<double, 2>::ID
{
  const double loss = 0.01;

  auto vm          = UD::World::VariableMesh<double, 2>();
  auto col_white   = UD::Colour::White();
  auto val_results = m_plan.AssessTotal(*m_candles, loss);
  auto prev_val    = val_results.begin()->second;
  auto sum_rate    = 1.0;
  vm.push_back(
      {{static_cast<double>(val_results.begin()->first.open_time.count()
                            + (val_results.begin()->first.close_time.count()
                               - val_results.begin()->first.open_time.count())
                                  / 2.0),
        sum_rate},
       col_white});
  for (auto& [key, val] : val_results)
  {
    if (val != prev_val)
    {
      prev_val = val;
      auto t   = key.open_time.count()
             + (key.close_time.count() - key.open_time.count()) / 2.0;
      auto col = col_white;
      sum_rate *= val;
      double bounded_sum_rate = std::min(2.0, sum_rate);

      col.blue() /= 2;
      col.red() /= 2;
      col.red() *= 2 - bounded_sum_rate;
      col.green() /= 2;
      col.green() *= bounded_sum_rate;
      vm.push_back({{static_cast<double>(t), val}, col});
    }
  }
  if (vm.size() < 1)
  {
    vm.push_back({{0.0, 0.0}, col_white});
  }

  return m_world.Add<UD::World::World<double, 2>::Entity>(vm);
}
auto OpenGLCanvas::CreateInd1Entity(
    std::map<Candle, Position, Candle::Compare::OpenTime> strategy_results)  //
    -> UD::World::World<double, 2>::ID
{
  auto col_blue    = UD::Colour::Blue();
  auto col_green   = UD::Colour::Green();
  auto col_red     = UD::Colour::Red();
  auto col_long    = col_blue + col_green;
  auto col_short   = col_blue + col_red;
  auto vm_ind1     = UD::World::VariableMesh<double, 2>();
  auto ind_results = m_strategy->RunFirst(*m_candles);

  for (auto& [key, val] : ind_results)
  {
    auto t = key.open_time.count()
           + (key.close_time.count() - key.open_time.count()) / 2;
    vm_ind1.push_back({{static_cast<double>(t), val},
                       (strategy_results.at(key) == Position::Long)
                           ? UD::Colour::Colour<>(col_long)
                           : UD::Colour::Colour<>(col_short)});
  }
  return m_world.Add<UD::World::World<double, 2>::Entity>(vm_ind1);
}
auto OpenGLCanvas::CreateInd2Entity(
    std::map<Candle, Position, Candle::Compare::OpenTime> strategy_results)  //
    -> UD::World::World<double, 2>::ID
{
  auto col_green    = UD::Colour::Green();
  auto col_red      = UD::Colour::Red();
  auto col_baseline = col_red + col_green;

  auto vm_ind2     = UD::World::VariableMesh<double, 2>();
  auto ind_results = m_strategy->RunSecond(*m_candles);

  for (auto& [key, val] : ind_results)
  {
    auto t = key.open_time.count()
           + (key.close_time.count() - key.open_time.count()) / 2;
    vm_ind2.push_back({{static_cast<double>(t), val}, col_baseline});
  }
  return m_world.Add<UD::World::World<double, 2>::Entity>(vm_ind2);
}
auto OpenGLCanvas::CreateIndEntities()  //
    -> std::pair<UD::World::World<double, 2>::ID,
                 UD::World::World<double, 2>::ID>
{
  auto strategy_results = m_strategy->Evaluate(*m_candles);
  return {CreateInd1Entity(strategy_results),
          CreateInd2Entity(strategy_results)};
}
void OpenGLCanvas::InitWorldSaveView()
{
  if (!m_candles)
    return;
  if (m_candles->size() <= 0)
    return;

  m_world.ClearEntities();
  auto entities = CreateEntityCandles(*m_candles, m_plan);

  m_strat_id = CreateStratEntity();
  m_world.Get<UD::World::World<double, 2>::Entity>(m_strat_id)
      .ReBaseY(m_view->bottom(), m_view->top());

   auto ind_ids = CreateIndEntities();
   m_ind1_id    = ind_ids.first;
   m_ind2_id    = ind_ids.second;

  InitMouseEntity();
  m_marker_id = CreateMarker();
  Refresh();
}
void OpenGLCanvas::InitWorld()
{
  if (!m_candles)
    return;
  if (m_candles->size() <= 0)
    return;

  m_world.ClearEntities();
  auto entities = CreateEntityCandles(*m_candles, m_plan);

  m_view = InitView(entities);

  m_strat_id = CreateStratEntity();
  m_world.Get<UD::World::World<double, 2>::Entity>(m_strat_id)
      .ReBaseY(m_view->bottom(), m_view->top());

   auto ind_ids = CreateIndEntities();
   m_ind1_id    = ind_ids.first;
   m_ind2_id    = ind_ids.second;

  InitMouseEntity();
  m_marker_id = CreateMarker();
  Refresh();
}
void OpenGLCanvas::OnPaint(wxPaintEvent& evt)
{
  if (!IsShown())
    return;

  auto dc = wxPaintDC(this);

  wxGLCanvas::SetCurrent(*m_context);

  auto frame            = GetSize();
  auto visible_entities = m_view->VisibleX({m_mouse_x_id,
                                            m_mouse_y_id,
                                            m_strat_id,
                                            m_ind1_id,
                                            m_ind2_id,
                                            m_marker_id});
  m_view->FitY(visible_entities);

  if (m_strat_id)
  {
    auto& strat_entity
        = m_world.Get<UD::World::World<double, 2>::Entity>(m_strat_id);
    strat_entity.ReBaseY(m_view->bottom(), m_view->top());
  }

  UpdateMouseEntity();
  UpdateMarkerEntity();
  m_view->Draw(
      {{0.0, 0.0},
       {static_cast<double>(frame.GetX()), static_cast<double>(frame.GetY())}});

  glFlush();
  SwapBuffers();
}

void OpenGLCanvas::OnMouseWheel(wxMouseEvent& evt)
{
  if (evt.GetWheelRotation() > 0)
  {
    m_view->ZoomWidth(0.8);
  }
  else
  {
    m_view->ZoomWidth(1.2);
  }
  Refresh();
}
void OpenGLCanvas::OnKeyDown(wxKeyEvent& evt) {}
auto OpenGLCanvas::GetCandleAtTime(std::chrono::milliseconds time)  //
    -> Candle
{
  for (auto& candle : *m_candles)
    if (time > candle.open_time && time < candle.close_time)
      return candle;
  throw std::exception();
}
auto OpenGLCanvas::GetPlanValueAtTime(std::chrono::milliseconds time)  //
    -> double
{
  for (auto& [candle, value] : m_plan.AssessTotal(*m_candles))
    if (time > candle.open_time && time < candle.close_time)
      return value;
  throw std::exception();
}
auto OpenGLCanvas::MousePosition()  //
    -> UD::Tensor::Vector<2, double>
{
  UD::Tensor::Vector<2, double> ret;
  auto                          mouselocation = wxGetMousePosition();
  ret.x() = mouselocation.x - this->GetScreenPosition().x;
  ret.y() = GetSize().GetY() - (mouselocation.y - this->GetScreenPosition().y);
  return ret;
}
auto OpenGLCanvas::Frame()  //
    -> UD::World::World<double, 2>::Rectangle
{
  return {{0.0, 0.0},
          {static_cast<double>(GetSize().GetX()),
           static_cast<double>(GetSize().GetY())}};
}
auto OpenGLCanvas::MouseWorldPosition()  //
    -> UD::Tensor::Vector<2, double>
{
  return m_view->FramePositionToWorldPosition(Frame(), MousePosition());
}
auto OpenGLCanvas::CreateMouseXAxisMesh()  //
    -> UD::World::HorizontalLine<double, 2>
{
  auto world_position
      = UD::World::World<double, 2>::Point{MouseWorldPosition()};
  world_position.colour() = UD::Colour::White();
  return {world_position, m_view->rect().width()};
}
auto OpenGLCanvas::CreateMouseYAxisMesh()  //
    -> UD::World::VerticalLine<double, 2>
{
  auto world_position
      = UD::World::World<double, 2>::Point{MouseWorldPosition()};
  world_position.colour() = UD::Colour::White();
  return {world_position, m_view->rect().height()};
}
void OpenGLCanvas::InitMouseEntity()
{
  m_mouse_x_id = m_world.Add<UD::World::World<double, 2>::Entity>(
      CreateMouseXAxisMesh());
  m_mouse_y_id = m_world.Add<UD::World::World<double, 2>::Entity>(
      CreateMouseYAxisMesh());
}
void OpenGLCanvas::UpdateMouseEntity()
{
  m_world.Get<UD::World::World<double, 2>::Entity>(m_mouse_x_id)
      .ExchangeMesh(CreateMouseXAxisMesh());
  m_world.Get<UD::World::World<double, 2>::Entity>(m_mouse_y_id)
      .ExchangeMesh(CreateMouseYAxisMesh());
}
void OpenGLCanvas::UpdateMarkerEntity()
{
  m_world.Get<UD::World::World<double, 2>::Entity>(m_marker_id)
      .ExchangeMesh(CreateMarkerMesh());
}
void OpenGLCanvas::UpdateTextMarkerPos()
{
  auto pos = m_view->WorldPositionToFramePosition(Frame(), m_marker);
  pos.y()  = Frame().height() - pos.y();
  m_text_marker_pos->SetLabel(
      wxString("(") << MillisecondsToString<std::chrono::minutes>(
          std::chrono::milliseconds(static_cast<long long>(m_marker.x())))
                    << "," << m_marker.y() << ")");
  m_text_marker_pos->SetPosition({(int)pos.x(), (int)pos.y()});
}
void OpenGLCanvas::OnMouseMotion(wxMouseEvent& evt)
{
  // m_world.Get<UD::World::World<double, 2>::Entity>(m_mouse_x_id).mesh() =
  // xaxis; m_world.Get<UD::World::World<double,
  // 2>::Entity>(m_mouse_y_id).mesh() = yaxis;
  if (m_dragging && evt.Dragging())
  {
    wxSize frame = GetSize();
    double x     = evt.GetX() - m_old_x;
    m_old_x      = evt.GetX();

    double diff = x / static_cast<double>(frame.GetX());

    double ref_scale = m_view->right() - m_view->left();
    double ref_diff  = ref_scale * diff;
    m_view->position().x() -= ref_diff;
  }
  Refresh();
}
void OpenGLCanvas::OnMouseLeaveWindow(wxMouseEvent& evt)
{
  m_dragging = false;
}
void OpenGLCanvas::OnMouseRightUp(wxMouseEvent& evt)
{
  m_marker = m_view->FramePositionToWorldPosition(Frame(), MousePosition());
  m_selected_entity_id = m_world.Intersect(m_marker);
  Candle candle;
  try
  {
    candle = GetCandleAtTime(
        std::chrono::milliseconds(static_cast<long long>(m_marker.x())));
  }
  catch (...)
  {
  }
  m_candle_view->Update(candle);
  m_plan_stats_view->Update(candle, *m_candles, m_plan);
}
void OpenGLCanvas::OnMouseLeftUp(wxMouseEvent& evt)
{
  // UpdateTextMarkerPos();
  // try {
  //  m_candle_view->Update(GetCandleAtTime(std::chrono::milliseconds(static_cast<long
  //  long>(m_marker.x()))));
  //} catch (...) {
  //}
  m_dragging = false;
  Refresh();
}
void OpenGLCanvas::OnMouseLeftDown(wxMouseEvent& evt)
{
  m_dragging = true;
  m_old_x    = evt.GetX();
}

const std::size_t OpenGLCanvas::MAX_ZOOM       = 1000;
const std::size_t OpenGLCanvas::MIN_ZOOM       = 100;
const std::size_t OpenGLCanvas::ZOOM_INCREMENT = 10;
