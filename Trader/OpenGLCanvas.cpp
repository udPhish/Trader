#include "OpenGLCanvas.h"

#include <iostream>

OpenGLCanvas::OpenGLCanvas(wxFrame* parent, const History* candles,
                           const Plan* plan, const wxGLAttributes& attributes)
    : m_candles(candles),
      m_plan(plan),
      wxGLCanvas(parent, attributes, wxID_ANY, wxDefaultPosition, wxDefaultSize,
                 wxFULL_REPAINT_ON_RESIZE),
      m_position{static_cast<int>(MAX_ZOOM) / 4},
      m_zoom{MAX_ZOOM / 2} {
  Bind(wxEVT_SIZE, &OpenGLCanvas::OnSize, this);
  Bind(wxEVT_PAINT, &OpenGLCanvas::OnPaint, this);

  Bind(wxEVT_MOUSEWHEEL, &OpenGLCanvas::OnMouseWheel, this);
  Bind(wxEVT_KEY_DOWN, &OpenGLCanvas::OnKeyDown, this);
  // Bind(wxEVT_LEFT_DOWN, &OpenGLCanvas::OnMouseLeftDown, this);
  // Bind(wxEVT_LEFT_UP, &OpenGLCanvas::OnMouseLeftUp, this);
  // Bind(wxEVT_MOTION, &OpenGLCanvas::OnMouseMotion, this);

  // wxGLContextAttrs context_attributes;
  // context_attributes.PlatformDefaults()
  //    .CoreProfile()
  //    .OGLVersion(3, 2)
  //    .EndList();
  m_context = std::make_unique<wxGLContext>(this);
  // m_context = std::make_unique<wxGLContext>(this, nullptr,
  // &context_attributes);
  // To avoid flashing on MSW
  // SetBackgroundStyle(wxBG_STYLE_CUSTOM);
}
OpenGLCanvas::OpenGLCanvas(wxFrame* parent, const wxGLAttributes& attributes)
    : OpenGLCanvas{parent, nullptr, nullptr, attributes} {}

OpenGLCanvas::~OpenGLCanvas() {}

void OpenGLCanvas::UpdateCandles(const History* candles) {
  m_candles = candles;
  Refresh();
}
void OpenGLCanvas::UpdatePlan(const Plan* plan) {
  m_plan = plan;
  Refresh();
}
void OpenGLCanvas::PrepareViewport() {
  wxSize frame = GetSize() * GetContentScaleFactor();
  wxGLCanvas::SetCurrent(*m_context);
  // glViewport(0, 0, frame.GetX(), frame.GetY());
  glViewport(0, 0, frame.GetX(), frame.GetY());
  // glScissor(0, 0, frame.GetWidth(), frame.GetHeight());

  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  glOrtho(0, 1, -0.1, 1.1, 0, 1);
  // if (frame.GetX() <= frame.GetY())
  //  glOrtho(0, 1, 0, (double)frame.GetHeight() / frame.GetWidth(), 0, 1);
  // else
  //  glOrtho(0, (double)frame.GetWidth() / frame.GetHeight(), 0, 1, 0, 1);
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();
}
void OpenGLCanvas::PrepareGL() {
  glEnable(GL_TEXTURE_2D);
  glEnable(GL_COLOR_MATERIAL);
  glEnable(GL_BLEND);
  glEnable(GL_SCISSOR_TEST);
  glDisable(GL_DEPTH_TEST);
  glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
}

void OpenGLCanvas::OnSize(wxSizeEvent& event) {
  // event.Skip();
  wxGLCanvas::SetCurrent(*m_context);
  glViewport(100, 100, 100, 100);
  Refresh();
  //event.Skip();
}
void OpenGLCanvas::OnPaint(wxPaintEvent& evt) {
  if (!IsShown()) return;
  if (!m_candles) return;
  if (m_candles->size() <= 0) return;

  wxPaintDC dc(this);

  wxGLCanvas::SetCurrent(*m_context);

  World<2> w;
  for (auto& candle : *m_candles) {
    w.Add<EntityCandle>(candle);
  }
  auto times = m_candles->ExtractOpenTime();
  auto lows = m_candles->ExtractLow();
  auto highs = m_candles->ExtractHigh();
  int min_time = times[0].count();
  int max_time = min_time;
  double lowest = lows[0];
  double highest = highs[0];
  for (std::size_t i = 1; i < times.size(); ++i) {
    min_time = std::min(min_time, times[i].count());
    max_time = std::max(max_time, times[i].count());
    lowest = std::min(lowest, lows[i]);
    highest = std::max(highest, highs[i]);
  }

  World<2>::Rectangle rect = {0, 0, max_time - min_time,
                              highest * 10000 - lowest * 10000};
  World<2>::Vector vec = {
      min_time + (max_time - min_time) / 2,
      lowest * 10000 + (highest * 10000 - lowest * 10000) / 2};
  auto& view = w.CreateView(rect, vec);

  wxSize frame = GetSize();
  Logger::Log(wxString("frame{") << frame.GetX() << "," << frame.GetY() << "}");
  view.Draw({100, 100, 100, 100});
  //view.Draw({0, 0, frame.GetX(), frame.GetY()});
  // auto& w

  // wxGLCanvas::SetCurrent(*m_context);
  // glEnable(GL_DEPTH_TEST);
  // glDepthFunc(GL_LEQUAL);
  // glEnable(GL_BLEND);
  // glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
  // glClearColor(0.0f, 0.0f, 0.0f, 1.0f);
  // glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  // double max_close = 0;
  // double min_close = 0;
  // if (m_candles->size() > 0) {
  //  auto strat = m_plan->Run(*m_candles);
  //  int last_index = m_candles->size() - (m_position - m_zoom / 2);
  //  int first_index = m_candles->size() - (m_position + m_zoom / 2);
  //  if (first_index >= m_candles->size()) return;
  //  if (last_index < 0) return;
  //  auto first = std::next(m_candles->begin(), std::max(first_index, 0));
  //  auto last = std::next(
  //      m_candles->begin(),
  //      std::min(static_cast<std::size_t>(last_index), m_candles->size()));
  //  max_close = first->close;
  //  min_close = first->close;
  //  for (auto it = first; it != last; ++it) {
  //    max_close = std::max(max_close, it->close);
  //    min_close = std::min(min_close, it->close);
  //  }
  //  glBegin(GL_LINES);
  //  auto it = first;
  //  int count = last_index - first_index;
  //  GLfloat prev_x = 0;
  //  GLfloat prev_y = 0;
  //  for (int i = first_index; i < last_index; ++i) {
  //    if (i >= 0 && i < m_candles->size()) {
  //      if (strat.at(*it) == Position::Long) {
  //        glColor4f(0, 1, 0, 1);
  //      } else {
  //        glColor4f(1, 0, 0, 1);
  //      }
  //      glVertex2f(prev_x, prev_y);
  //      prev_x = (static_cast<double>(i) - first_index) / count;
  //      prev_y = (it->close - min_close) / (max_close - min_close);
  //      glVertex2f(prev_x, prev_y);
  //      ++it;
  //    }
  //  }
  //  // for (std::size_t i = 0; i < m_candles->size(); ++i) {
  //  //  auto current = std::next(m_candles->begin(), i);
  //  //  if (strat.at(current->open_time) == Strategy::Position::Long) {
  //  //    glColor4f(0, 1, 0, 1);
  //  //  } else {
  //  //    glColor4f(1, 0, 0, 1);
  //  //  }
  //  //  glVertex2f((double)i / (double)m_candles->size(),
  //  //             (current->close - min_close) / (max_close - min_close));
  //  //}
  //  glEnd();
  //}

  glFlush();
  SwapBuffers();
}

void OpenGLCanvas::OnMouseWheel(wxMouseEvent& evt) {
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
// void OpenGLCanvas::OnMouseLeftUp(wxMouseEvent& evt) {
//  m_dragging = false;
//  evt.Skip();
//}
// void OpenGLCanvas::OnMouseLeftDown(wxMouseEvent& evt) {
//  m_dragging = true;
//  evt.Skip();
//}
// void OpenGLCanvas::OnMouseMotion(wxMouseEvent& evt) {
//  if (evt.LeftIsDown() && evt.Dragging()) {
//  }
//}

const std::size_t OpenGLCanvas::MAX_ZOOM = 1000;
const std::size_t OpenGLCanvas::MIN_ZOOM = 100;
const std::size_t OpenGLCanvas::ZOOM_INCREMENT = 10;
