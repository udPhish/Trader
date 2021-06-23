#include "OpenGLCanvas.h"

#include <iostream>

OpenGLCanvas::OpenGLCanvas(wxFrame* parent, const History* candles,
                           const Plan* plan, const wxGLAttributes& attributes)
    : m_candles(candles),
      m_plan(plan),
      wxGLCanvas(parent, attributes, wxID_ANY, wxDefaultPosition, wxDefaultSize,
                 wxFULL_REPAINT_ON_RESIZE),
      m_position{static_cast<int>(MAX_ZOOM) / 4},
      m_zoom{MAX_ZOOM / 2},
      m_testing{false},
      m_world{},
      m_dragging{false},
      m_old_x{0},
      m_view{&m_world.CreateView(
          UD::World::World<double, 2>::Rectangle{{0.0, 0.0}, {0.0, 0.0}})},
      m_strat_view{&m_world.CreateView(
          UD::World::World<double, 2>::Rectangle{{0.0, 0.0}, {0.0, 0.0}})},
      m_strat_id{},
      m_mouse_x_id{},
      m_mouse_y_id{} {
  Bind(wxEVT_SIZE, &OpenGLCanvas::OnSize, this);
  Bind(wxEVT_PAINT, &OpenGLCanvas::OnPaint, this);

  Bind(wxEVT_MOUSEWHEEL, &OpenGLCanvas::OnMouseWheel, this);
  Bind(wxEVT_KEY_DOWN, &OpenGLCanvas::OnKeyDown, this);
  Bind(wxEVT_MOTION, &OpenGLCanvas::OnMouseMotion, this);
  Bind(wxEVT_LEAVE_WINDOW, &OpenGLCanvas::OnMouseLeaveWindow, this);
  Bind(wxEVT_LEFT_UP, &OpenGLCanvas::OnMouseLeftUp, this);
  Bind(wxEVT_LEFT_DOWN, &OpenGLCanvas::OnMouseLeftDown, this);

  // Bind(wxEVT_LEFT_DOWN, &OpenGLCanvas::OnMouseLeftDown, this);
  // Bind(wxEVT_LEFT_UP, &OpenGLCanvas::OnMouseLeftUp, this);
  // Bind(wxEVT_MOTION, &OpenGLCanvas::OnMouseMotion, this);

  // TODO: Switch to the following when move to Shaders:
  // m_context_attributes.PlatformDefaults()
  //    .CoreProfile()
  //    .OGLVersion(3, 2)
  //    .EndList();
  // m_context =
  //    std::make_unique<wxGLContext>(this, nullptr, &m_context_attributes);

  m_context = std::make_unique<wxGLContext>(this);

  // wxGLContextAttrs context_attributes;
  // context_attributes.PlatformDefaults()
  //    .CoreProfile()
  //    .OGLVersion(3, 2)
  //    .EndList();
  // m_context = std::make_unique<wxGLContext>(this, nullptr,
  // &context_attributes);
  // To avoid flashing on MSW
  SetBackgroundStyle(wxBG_STYLE_CUSTOM);
}
OpenGLCanvas::OpenGLCanvas(wxFrame* parent, const wxGLAttributes& attributes)
    : OpenGLCanvas{parent, nullptr, nullptr, attributes} {}

OpenGLCanvas::~OpenGLCanvas() {}

void OpenGLCanvas::ShouldTest(bool should) { m_testing = should; }
void OpenGLCanvas::UpdateCandles(const History* candles) {
  m_candles = candles;
  InitWorld();
  Refresh();
}
void OpenGLCanvas::UpdatePlan(const Plan* plan) {
  m_plan = plan;
  InitWorld();
  Refresh();
}
void OpenGLCanvas::PrepareViewport() {
  PrepareGL();
  wxSize frame = GetSize() * GetContentScaleFactor();
  // glViewport(0, 0, frame.GetX(), frame.GetY());
  glViewport(0, 0, frame.GetX() / 2, frame.GetY() / 2);
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  glOrtho(0, frame.GetWidth(), frame.GetHeight(), 0, 0, 1);
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();
}
void OpenGLCanvas::PrepareGL() {
  glEnable(GL_TEXTURE_2D);
  glEnable(GL_COLOR_MATERIAL);
  glEnable(GL_BLEND);
  glDisable(GL_DEPTH_TEST);
  glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
}

void OpenGLCanvas::OnSize(wxSizeEvent& event) {
  // event.Skip();
  // wxGLCanvas::SetCurrent(*m_context);

  // wxSize frame = GetSize();  //*GetContentScaleFactor();
  // glViewport(0, 0, frame.GetX()/2, frame.GetY()/2);
  Refresh();
  // event.Skip();
}
void OpenGLCanvas::InitWorld() {
  if (!m_candles) return;
  if (m_candles->size() <= 0) return;
  m_world.Clear();
  m_view = &m_world.CreateView(
      UD::World::World<double, 2>::Rectangle{{0.0, 0.0}, {0.0, 0.0}});
  m_strat_view = &m_world.CreateView(
      UD::World::World<double, 2>::Rectangle{{0.0, 0.0}, {0.0, 0.0}});

  auto plan_results = m_plan->Run(*m_candles);
  auto val_results = m_plan->Assess(*m_candles);
  std::vector<UD::World::World<double, 2>::ID> entities;
  for (auto& candle : *m_candles) {
    entities.push_back(m_world.Add<EntityCandle>(candle, plan_results[candle]));
  }

  auto& temp_ent = *m_world.entities.begin()->second;
  double max_height = temp_ent.mesh().MaxY() + temp_ent.position().y();
  double min_height = temp_ent.mesh().MinY() + temp_ent.position().y();
  for (auto& ent : m_world.entities) {
    double maxy = ent.second->mesh().MaxY() + ent.second->position().y();
    double miny = ent.second->mesh().MinY() + ent.second->position().y();
    max_height = maxy > max_height ? maxy : max_height;
    min_height = miny < min_height ? miny : min_height;
  }
  auto vm = UD::World::VariableMesh<double, 2>();
  auto col_white = UD::Colour::White();
  double min_val = val_results.begin()->second;
  double max_val = min_val;
  double cumulative = 0;
  std::vector<double> valvec;
  for (auto& res : val_results) {
    cumulative += res.second;
    if (cumulative < min_val) min_val = cumulative;
    if (cumulative > max_val) max_val = cumulative;
    valvec.push_back(cumulative);
  }
  for (auto& v : valvec) {
    auto basedv = (v - min_val);
    auto ratev = (basedv / (max_val - min_val));
    auto newrate = ratev * (max_height - min_height);
    auto newv = newrate + min_height;

    v = ((v - min_val) / (max_val - min_val)) * (max_height - min_height) +
        min_height;
  }

  UD::Math::Type::ULong icount = 0;
  for (auto& res : val_results) {
    auto t = res.first.open_time.count() +
             (res.first.close_time.count() - res.first.open_time.count()) / 2;
    vm.push_back({{static_cast<double>(t), valvec[icount]}, col_white});
    ++icount;
  }
  m_strat_id = m_world.Add<UD::World::World<double, 2>::Entity>(vm);
  entities.push_back(m_strat_id);
  m_view->Fit(entities);
  // std::vector<UD::World::World<double, 2>::ID> strats;
  // strats.push_back(m_strat_id);
  // m_strat_view->Fit(strats);

  auto mouselocation = wxGetMousePosition();
  double mousex = mouselocation.x - this->GetScreenPosition().x;
  double mousey = mouselocation.y - this->GetScreenPosition().y;

  wxSize framesize = GetSize();

  mousex =
      m_view->left() + (mousex / framesize.GetX()) * m_view->rect().width();
  mousey = m_view->top() - (mousey / framesize.GetY()) * m_view->rect().width();

  auto xaxis = UD::World::StaticMesh<double, 2, 3>{
      UD::World::World<double, 2>::Point{{m_view->left(), mousey}, col_white},
      UD::World::World<double, 2>::Point{
          {m_view->left() + m_view->rect().width() / 2, mousey}, col_white},
      UD::World::World<double, 2>::Point{{m_view->right(), mousey}, col_white}};
  auto yaxis = UD::World::StaticMesh<double, 2, 3>{
      UD::World::World<double, 2>::Point{{mousex, m_view->bottom()}, col_white},
      UD::World::World<double, 2>::Point{
          {mousex, m_view->bottom() + m_view->rect().height() / 2}, col_white},
      UD::World::World<double, 2>::Point{{mousex, m_view->top()}, col_white}};
  m_mouse_x_id = m_world.Add<UD::World::World<double, 2>::Entity>(xaxis);
  m_mouse_y_id = m_world.Add<UD::World::World<double, 2>::Entity>(yaxis);
  // entities.push_back(m_mouse_x_id);
  // entities.push_back(m_mouse_y_id);
}
void OpenGLCanvas::OnPaint(wxPaintEvent& evt) {
  if (!IsShown()) return;

  wxPaintDC dc(this);

  wxGLCanvas::SetCurrent(*m_context);

  wxSize frame = GetSize();
  auto& strat_entity =
      m_world.Get<UD::World::World<double, 2>::Entity>(m_strat_id);

  m_view->FitYExcluding(m_view->VisibleX(),
                        {m_mouse_x_id, m_mouse_y_id, m_strat_id});
  strat_entity.ScaleYTo(m_view->bottom(), m_view->top());
  strat_entity.ReCenterY(m_view->center().y());

  m_view->Draw(
      {{0.0, 0.0},  // static_cast<double>(frame.GetY()) / 3.0},
       {static_cast<double>(frame.GetX()), static_cast<double>(frame.GetY())}});

  // m_strat_view->Draw({{0.0, 0.0},
  //                   {static_cast<double>(frame.GetX()),
  //                    static_cast<double>(frame.GetY())/3.0 }},
  //                  {m_strat_id, m_mouse_x_id});

  glFlush();
  SwapBuffers();
}

void OpenGLCanvas::OnMouseWheel(wxMouseEvent& evt) {
  if (evt.GetWheelRotation() > 0) {
    m_view->ZoomWidth(0.8);
  } else {
    m_view->ZoomWidth(1.2);
  }
  Refresh();
  // if (m_candles->empty()) return;
  // if (evt.GetWheelRotation() > 0) {
  //  if (m_zoom >= MIN_ZOOM + ZOOM_INCREMENT) {
  //    m_zoom -= ZOOM_INCREMENT;
  //  }
  //} else {
  //  if (m_zoom <= MAX_ZOOM - ZOOM_INCREMENT) {
  //    m_zoom += ZOOM_INCREMENT;
  //    m_position = std::min(
  //        static_cast<int>(m_candles->size()) - (static_cast<int>(m_zoom) /
  //        2), m_position);
  //    m_position = std::max(static_cast<int>(m_zoom) / 2, m_position);
  //  }
  //}
  // Refresh();
}
void OpenGLCanvas::OnKeyDown(wxKeyEvent& evt) {
  // if (m_candles->empty()) return;
  // std::size_t amount = 1;
  // if (wxGetKeyState(WXK_SHIFT)) amount = m_zoom;
  // switch (evt.GetKeyCode()) {
  //  case wxKeyCode::WXK_LEFT:
  //    m_position += amount;
  //    m_position = std::min(
  //        static_cast<int>(m_candles->size()) - (static_cast<int>(m_zoom) /
  //        2), m_position);
  //    break;
  //  case wxKeyCode::WXK_RIGHT:
  //    m_position -= amount;
  //    m_position = std::max(static_cast<int>(m_zoom) / 2, m_position);
  //    break;
  //}
  // Refresh();
}
void OpenGLCanvas::OnMouseLeftUp(wxMouseEvent& evt) { m_dragging = false; }
void OpenGLCanvas::OnMouseLeftDown(wxMouseEvent& evt) {
  m_dragging = true;
  m_old_x = evt.GetX();
}
void OpenGLCanvas::OnMouseMotion(wxMouseEvent& evt) {
  auto col_white = UD::Colour::White();
  auto mouselocation = wxGetMousePosition();
  double mousex = mouselocation.x - this->GetScreenPosition().x;
  double mousey = mouselocation.y - this->GetScreenPosition().y;

  wxSize framesize = GetSize();

  mousex =
      m_view->left() + (mousex / framesize.GetX()) * m_view->rect().width();
  mousey =
      m_view->top() - (mousey / framesize.GetY()) * m_view->rect().height();

  auto xaxis = UD::World::StaticMesh<double, 2, 3>{
      UD::World::World<double, 2>::Point{{m_view->left(), mousey}, col_white},
      UD::World::World<double, 2>::Point{
          {m_view->left() + m_view->rect().width() / 2, mousey}, col_white},
      UD::World::World<double, 2>::Point{{m_view->right(), mousey}, col_white}};
  auto yaxis = UD::World::StaticMesh<double, 2, 3>{
      UD::World::World<double, 2>::Point{{mousex, m_view->bottom()}, col_white},
      UD::World::World<double, 2>::Point{
          {mousex, m_view->bottom() + m_view->rect().height() / 2}, col_white},
      UD::World::World<double, 2>::Point{{mousex, m_view->top()}, col_white}};

  m_world.Get<UD::World::World<double, 2>::Entity>(m_mouse_x_id)
      .ExchangeMesh(xaxis);
  m_world.Get<UD::World::World<double, 2>::Entity>(m_mouse_y_id)
      .ExchangeMesh(yaxis);
  // m_world.Get<UD::World::World<double, 2>::Entity>(m_mouse_x_id).mesh() =
  // xaxis; m_world.Get<UD::World::World<double,
  // 2>::Entity>(m_mouse_y_id).mesh() = yaxis;

  if (m_dragging && evt.Dragging()) {
    wxSize frame = GetSize();
    double x = evt.GetX() - m_old_x;
    m_old_x = evt.GetX();

    double diff = x / static_cast<double>(frame.GetX());

    double ref_scale = m_view->right() - m_view->left();
    double ref_diff = ref_scale * diff;
    m_view->position().x() -= ref_diff;
  }
  Refresh();
}
void OpenGLCanvas::OnMouseLeaveWindow(wxMouseEvent& evt) { m_dragging = false; }

const std::size_t OpenGLCanvas::MAX_ZOOM = 1000;
const std::size_t OpenGLCanvas::MIN_ZOOM = 100;
const std::size_t OpenGLCanvas::ZOOM_INCREMENT = 10;
