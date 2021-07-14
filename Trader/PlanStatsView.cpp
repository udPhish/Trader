#include "PlanStatsView.h"

PlanStatsView::PlanStatsView(wxWindow* parent) : wxWindow{parent, wxID_ANY}
{
  _previous_position     = new wxStaticText(this, wxID_ANY, "value");
  _current_position      = new wxStaticText(this, wxID_ANY, "value");
  _first_entry_position  = new wxStaticText(this, wxID_ANY, "value");
  _second_entry_position = new wxStaticText(this, wxID_ANY, "value");
  _first_exit_position   = new wxStaticText(this, wxID_ANY, "value");
  _second_exit_position  = new wxStaticText(this, wxID_ANY, "value");
  _gain_loss             = new wxStaticText(this, wxID_ANY, "value");
  _total                 = new wxStaticText(this, wxID_ANY, "value");
  _final_total           = new wxStaticText(this, wxID_ANY, "value");
  auto* pos_title        = new wxStaticText(this, wxID_ANY, "Positions");
  auto* values_title     = new wxStaticText(this, wxID_ANY, "Values");
  pos_title->Center();
  values_title->Center();

  auto* outer_s = new wxBoxSizer(wxVERTICAL);
  outer_s->Add(pos_title, wxSizerFlags().Expand().Proportion(1));

  auto* s = new wxFlexGridSizer(8, 2, 2, 10);
  s->Add(new wxStaticText(this, wxID_ANY, "Previous:"));
  s->Add(_previous_position);
  s->Add(new wxStaticText(this, wxID_ANY, "Current:"));
  s->Add(_current_position);
  s->Add(new wxStaticText(this, wxID_ANY, "--Entry--"));
  s->Add(new wxStaticText(this, wxID_ANY, ""));
  s->Add(new wxStaticText(this, wxID_ANY, "First:"));
  s->Add(_first_entry_position);
  s->Add(new wxStaticText(this, wxID_ANY, "Second:"));
  s->Add(_second_entry_position);
  s->Add(new wxStaticText(this, wxID_ANY, "--Exit--"));
  s->Add(new wxStaticText(this, wxID_ANY, ""));
  s->Add(new wxStaticText(this, wxID_ANY, "First:"));
  s->Add(_first_exit_position);
  s->Add(new wxStaticText(this, wxID_ANY, "Second:"));
  s->Add(_second_exit_position);
  s->AddGrowableCol(1);
  outer_s->Add(s);

  outer_s->Add(values_title);

  auto* s2 = new wxFlexGridSizer(3, 2, 2, 10);
  s2->Add(new wxStaticText(this, wxID_ANY, "Gain:"));
  s2->Add(_gain_loss);
  s2->Add(new wxStaticText(this, wxID_ANY, "Total:"));
  s2->Add(_total);
  s2->Add(new wxStaticText(this, wxID_ANY, "Final Total:"));
  s2->Add(_final_total);
  s2->AddGrowableCol(1);
  outer_s->Add(s2);

  SetSizer(outer_s);
  SetAutoLayout(true);
}
void PlanStatsView::Update(const Candle&  candle,
                           const History& candles,
                           const Plan&    plan)
{
  auto run_results          = plan.Run(candles);
  auto assess_results       = plan.Assess(candles);
  auto assess_total_results = plan.AssessTotal(candles);
  auto first_entry_run      = plan.entries()[0]->Evaluate(candles);
  auto second_entry_run     = plan.entries()[1]->Evaluate(candles);
  auto first_exit_run       = plan.exits()[0]->Evaluate(candles);
  auto second_exit_run      = plan.exits()[1]->Evaluate(candles);

  auto prev_it = run_results.begin();
  for (auto it = run_results.begin(); it != run_results.end(); ++it)
  {
    if (candle.open_time == it->first.open_time)
    {
      break;
    }
    prev_it = it;
  }
  _previous_position->SetLabel(
      run_results[prev_it->first] == Position::Long ? "Long" : "Short");
  _current_position->SetLabel(run_results[candle] == Position::Long ? "Long"
                                                                    : "Short");
  _first_entry_position->SetLabel(
      first_entry_run[candle] == Position::Long ? "Long" : "Short");
  _second_entry_position->SetLabel(
      second_entry_run[candle] == Position::Long ? "Long" : "Short");
  _first_exit_position->SetLabel(
      first_exit_run[candle] == Position::Long ? "Long" : "Short");
  _second_exit_position->SetLabel(
      second_exit_run[candle] == Position::Long ? "Long" : "Short");
  _gain_loss->SetLabel(assess_results.contains(candle)
                           ? wxString() << assess_results[candle]
                           : wxString());
  _total->SetLabel(wxString() << assess_total_results[candle]);
  _final_total->SetLabel(wxString() << plan.AssessFinalTotal(candles));
}