#pragma once

#include <set>

#include "wx/wx.h"

#include "wx/glcanvas.h"

#include "Candle.h"
#include "EntityCandle.h"
#include "Logger.h"
#include "Strategies.h"
#include "Strategy.h"
#include "Util.h"
#include "World.h"

class OpenGLCanvas : public wxGLCanvas {
  static const std::size_t MAX_ZOOM;
  static const std::size_t MIN_ZOOM;
  static const std::size_t ZOOM_INCREMENT;
  wxGLContextAttrs m_context_attributes;
  std::unique_ptr<wxGLContext> m_context;
  //wxGLContext* m_context;
  const History* m_candles;
  const Plan* m_plan;
  std::size_t m_zoom;
  int m_position;
  bool m_testing;
  World<2> m_world;
  World<2>::View m_view;
  // bool m_dragging;

  void OnSize(wxSizeEvent& event);
  void OnPaint(wxPaintEvent& evt);
  void OnMouseWheel(wxMouseEvent& evt);
  void OnKeyDown(wxKeyEvent& evt);
  // void OnMouseLeftUp(wxMouseEvent& evt);
  // void OnMouseLeftDown(wxMouseEvent& evt);
  // void OnMouseMotion(wxMouseEvent& evt);
  void PrepareViewport();
  void PrepareGL();

 public:
  OpenGLCanvas(wxFrame* parent,  const History* candles,
               const Plan* plan,
               const wxGLAttributes& attributes);
  OpenGLCanvas(wxFrame* parent, const wxGLAttributes& attributes);
  virtual ~OpenGLCanvas();
  void DrawTest();
  void ShouldTest(bool should = true);
  void UpdateCandles( const History* candles);
  void UpdatePlan(const Plan* plan);
};
