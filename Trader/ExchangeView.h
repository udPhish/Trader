#pragma once
#include <wx/listbox.h>
#include <wx/choice.h>
#include <wx/wx.h>

#include "Exchange.h"
#include "Logger.h"
class ExchangeView : public wxWindow {
  wxListBox* m_list_markets;
  wxChoice* m_choice_timeframe;
  Exchange& m_exchange;
 public:
  ExchangeView(wxWindow* parent, Exchange& exchange);
  History& SelectedCandles();
};
