#pragma once
#include <algorithm>
#include <chrono>
#include <set>
#include <sstream>
#include "Candle.h"
#include "CandleView.h"
#include "EntityCandle.h"
#include "Logger.h"
#include "Math.h"
#include "PlanStatsView.h"
#include "Strategies.h"
#include "Strategy.h"
#include "Util.h"
#include "World.h"
#include "date.h"
#include "wx/glcanvas.h"
#include "wx/wx.h"
class OpenGLCanvas : public wxGLCanvas
{
  static const std::size_t            MAX_ZOOM;
  static const std::size_t            MIN_ZOOM;
  static const std::size_t            ZOOM_INCREMENT;
  wxGLContextAttrs                    m_context_attributes = {};
  const History*                      m_candles            = nullptr;
  Plan                                m_plan               = {};
  std::unique_ptr<wxGLContext>        m_context            = nullptr;
  std::unique_ptr<StrategyBase>       m_strategy           = nullptr;
  UD::World::World<double, 2>         m_world              = {};
  UD::World::World<double, 2>::View*  m_view               = nullptr;
  UD::World::World<double, 2>::View*  m_strat_view         = nullptr;
  UD::World::World<double, 2>::ID     m_strat_id           = 0;
  UD::World::World<double, 2>::ID     m_ind1_id            = 0;
  UD::World::World<double, 2>::ID     m_ind2_id            = 0;
  UD::World::World<double, 2>::ID     m_mouse_x_id         = 0;
  UD::World::World<double, 2>::ID     m_mouse_y_id         = 0;
  UD::World::World<double, 2>::ID     m_marker_id          = 0;
  UD::World::World<double, 2>::ID     m_selected_entity_id = 0;
  UD::World::World<double, 2>::Vector m_marker             = {};
  CandleView*                         m_candle_view        = nullptr;
  PlanStatsView*                      m_plan_stats_view    = nullptr;
  wxStaticText*                       m_text_marker_pos    = nullptr;
  bool                                m_dragging           = false;
  double                              m_old_x              = 0.0;

  void OnSize(wxSizeEvent& event);
  void OnPaint(wxPaintEvent& evt);
  void OnMouseWheel(wxMouseEvent& evt);
  void OnKeyDown(wxKeyEvent& evt);
  void OnMouseLeftUp(wxMouseEvent& evt);
  void OnMouseRightUp(wxMouseEvent& evt);
  void OnMouseLeftDown(wxMouseEvent& evt);
  void OnMouseMotion(wxMouseEvent& evt);
  void OnMouseLeaveWindow(wxMouseEvent& evt);
  void PrepareViewport();
  void PrepareGL();
  void InitMouseEntity();
  void UpdateMouseEntity();
  auto CreateMouseXAxisMesh()  //
      -> UD::World::HorizontalLine<double, 2>;
  auto CreateMouseYAxisMesh()  //
      -> UD::World::VerticalLine<double, 2>;
  auto CreateMarkerMesh()  //
      -> UD::World::VerticalLine<double, 2>;

 public:
  virtual ~OpenGLCanvas() = default;
  OpenGLCanvas(wxFrame*              parent,
               PlanStatsView*        plan_stats_view,
               CandleView*           candle_view,
               const History*        candles,
               const wxGLAttributes& attributes);
  OpenGLCanvas(wxFrame*              parent,
               PlanStatsView*        plan_stats_view,
               CandleView*           candle_view,
               const wxGLAttributes& attributes);
  void InitWorld();
  void InitWorldSaveView();
  void UpdateCandles(const History* candles);
  void UpdatePlan(Plan plan);
  void UpdateStrategy(std::unique_ptr<StrategyBase> strategy);
  void UpdateMarkerEntity();
  void UpdateTextMarkerPos();
  auto MousePosition()  //
      -> UD::Tensor::Vector<2, double>;
  auto Frame()  //
      -> UD::World::World<double, 2>::Rectangle;
  auto MouseWorldPosition()  //
      -> UD::Tensor::Vector<2, double>;
  auto CreateEntityCandles(const History& candles, const Plan& plan)  //
      -> std::vector<UD::World::World<double, 2>::ID>;
  auto CreateStratEntity()  //
      -> UD::World::World<double, 2>::ID;
  auto CreateInd1Entity(
      std::map<Candle, Position, Candle::Compare::OpenTime> strat_results)  //
      -> UD::World::World<double, 2>::ID;
  auto CreateInd2Entity(
      std::map<Candle, Position, Candle::Compare::OpenTime> strat_results)  //
      -> UD::World::World<double, 2>::ID;
  auto CreateIndEntities()  //
      -> std::pair<UD::World::World<double, 2>::ID,
                   UD::World::World<double, 2>::ID>;
  auto InitView(
      const std::vector<UD::World::World<double, 2>::ID>& entities)  //
      -> UD::World::World<double, 2>::View*;
  auto CreateMarker()  //
      -> UD::World::World<double, 2>::ID;
  auto GetCandleAtTime(std::chrono::milliseconds time)  //
      -> Candle;
  auto GetPlanValueAtTime(std::chrono::milliseconds time)  //
      -> double;
  auto AccumulateAllPreviousElements(
      const std::map<Candle, double, Candle::Compare::OpenTime>& candles)  //
      -> std::map<Candle, double, Candle::Compare::OpenTime>;
  auto MultiplyAllPreviousElements(
      const std::map<Candle, double, Candle::Compare::OpenTime>& candles)  //
      -> std::map<Candle, double, Candle::Compare::OpenTime>;
};
