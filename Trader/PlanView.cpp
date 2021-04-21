#include "PlanView.h"

PlanView::PlanView(wxWindow* parent, std::vector<Plan> plans)
    : wxWindow(parent, wxID_ANY),
      m_list_plans{new wxListBox{this, wxID_ANY, wxDefaultPosition,
                                 wxDefaultSize, 0, nullptr,
                                 wxLB_SINGLE | wxLB_SORT}} {
  Add(plans);
  m_list_plans->Select(0);
  wxBoxSizer* s = new wxBoxSizer(wxVERTICAL);

  s->Add(m_list_plans, wxSizerFlags().Expand().Proportion(1));
  SetSizer(s);
  SetAutoLayout(true);
}

Plan& PlanView::SelectedPlan() {
  std::string plan_selection = m_list_plans->GetStringSelection().ToStdString();
  if (plan_selection == "")
    Logger::Log<Logger::Level::FatalError>("Invalid plan selection");
  return Get(plan_selection);
}
void PlanView::Add(Plan plan) {
  m_plans.insert({plan.Name(), plan});
  m_list_plans->Append(plan.Name());
}
void PlanView::Add(std::vector<Plan> plans) {
  for (auto& plan : plans) Add(plan);
}
Plan& PlanView::Get(std::string name) { return m_plans.at(name); }
void PlanView::Clear() {
  m_plans.clear();
  m_list_plans->Clear();
}