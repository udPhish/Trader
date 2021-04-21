#pragma once
#include <set>

#include <wx/listbox.h>
#include <wx/wx.h>

#include "Logger.h"
#include "Strategy.h"

class PlanView : public wxWindow {
  wxListBox* m_list_plans;
  std::unordered_map<std::string, Plan> m_plans;

 public:
  PlanView(wxWindow* parent, std::vector<Plan> plans);
  Plan& SelectedPlan();
  void Add(Plan plan);
  void Add(std::vector<Plan> plans);
  Plan& Get(std::string);
  void Clear();
};