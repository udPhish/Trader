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
  std::unique_ptr<wxGLContext> m_context;
  const History* m_candles;
  const Plan* m_plan;
  std::size_t m_zoom;
  int m_position;
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
  void UpdateCandles( const History* candles);
  void UpdatePlan(const Plan* plan);
};
