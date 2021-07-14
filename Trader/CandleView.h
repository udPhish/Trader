#pragma once
#include <wx/wx.h>

#include <wx/sizer.h>
#include <wx/stattext.h>

#include "Candle.h"
#include "Util.h"

struct CandleView : public wxWindow
{
  wxStaticText* _open_time;
  wxStaticText* _close_time;
  wxStaticText* _open;
  wxStaticText* _close;
  wxStaticText* _high;
  wxStaticText* _low;
  wxStaticText* _volume;
  CandleView(wxWindow* parent);
  void Update(const Candle& candle);
  void UpdateOpenTime(const Candle& candle);
  void UpdateCloseTime(const Candle& candle);
  void UpdateOpen(const Candle& candle);
  void UpdateClose(const Candle& candle);
  void UpdateHigh(const Candle& candle);
  void UpdateLow(const Candle& candle);
  void UpdateVolume(const Candle& candle);
};
