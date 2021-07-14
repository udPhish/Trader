#pragma once
#include <wx/wx.h>

#include <wx/sizer.h>
#include <wx/stattext.h>

#include "Candle.h"
#include "Strategy.h"

struct PlanStatsView : public wxWindow {
  wxStaticText* _previous_position;
  wxStaticText* _current_position;
  wxStaticText* _first_entry_position;
  wxStaticText* _second_entry_position;
  wxStaticText* _first_exit_position;
  wxStaticText* _second_exit_position;
  wxStaticText* _gain_loss;
  wxStaticText* _total;
  wxStaticText* _final_total;
  PlanStatsView(wxWindow* parent);
  void Update(const Candle& candle, const History& candles, const Plan& plan);
};
