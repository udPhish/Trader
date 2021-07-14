#include "CandleView.h"

CandleView::CandleView(wxWindow* parent)
    : wxWindow{parent, wxID_ANY}
{
  _open_time  = new wxStaticText(this, wxID_ANY, "value");
  _close_time = new wxStaticText(this, wxID_ANY, "value");
  _open       = new wxStaticText(this, wxID_ANY, "value");
  _close      = new wxStaticText(this, wxID_ANY, "value");
  _high       = new wxStaticText(this, wxID_ANY, "value");
  _low        = new wxStaticText(this, wxID_ANY, "value");
  _volume     = new wxStaticText(this, wxID_ANY, "value");
  auto* s     = new wxFlexGridSizer(7, 2, 2, 10);
  s->Add(new wxStaticText(this, wxID_ANY, "Open Time:"));
  s->Add(_open_time);
  s->Add(new wxStaticText(this, wxID_ANY, "Close Time:"));
  s->Add(_close_time);
  s->Add(new wxStaticText(this, wxID_ANY, "Open:"));
  s->Add(_open);
  s->Add(new wxStaticText(this, wxID_ANY, "Close:"));
  s->Add(_close);
  s->Add(new wxStaticText(this, wxID_ANY, "High:"));
  s->Add(_high);
  s->Add(new wxStaticText(this, wxID_ANY, "Low:"));
  s->Add(_low);
  s->Add(new wxStaticText(this, wxID_ANY, "Volume:"));
  s->Add(_volume);
  s->AddGrowableCol(1);
  SetSizer(s);
  SetAutoLayout(true);
}
void CandleView::UpdateOpenTime(const Candle& candle)
{
  _open_time->SetLabel(MillisecondsToString<std::chrono::minutes>(candle.open_time));
}
void CandleView::UpdateCloseTime(const Candle& candle)
{
  _close_time->SetLabel(MillisecondsToString<std::chrono::minutes>(candle.close_time));
}
void CandleView::UpdateOpen(const Candle& candle)
{
  _open->SetLabel((wxString() << candle.open).ToStdString());
}
void CandleView::UpdateClose(const Candle& candle)
{
  _close->SetLabel((wxString() << candle.close).ToStdString());
}
void CandleView::UpdateHigh(const Candle& candle)
{
  _high->SetLabel((wxString() << candle.high).ToStdString());
}
void CandleView::UpdateLow(const Candle& candle)
{
  _low->SetLabel((wxString() << candle.low).ToStdString());
}
void CandleView::UpdateVolume(const Candle& candle)
{
  _volume->SetLabel((wxString() << candle.volume).ToStdString());
}
void CandleView::Update(const Candle& candle)
{
  UpdateOpenTime(candle);
  UpdateCloseTime(candle);
  UpdateOpen(candle);
  UpdateClose(candle);
  UpdateHigh(candle);
  UpdateLow(candle);
  UpdateVolume(candle);
}