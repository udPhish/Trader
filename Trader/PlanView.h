#pragma once
#include <map>
#include <set>

#include <wx/listbox.h>
#include <wx/thread.h>
#include <wx/wx.h>

#include "Candle.h"
#include "Logger.h"
#include "Strategy.h"

class PlanView
    : public wxWindow
    , public wxThreadHelper
{
  enum ThreadEvent
  {
    Notify,
    Update
  };
  wxListBox*                            m_list_plans          = nullptr;
  wxListBox*                            m_plan_entries        = nullptr;
  wxListBox*                            m_plan_exits          = nullptr;
  History                               m_candles             = {};
  std::unordered_map<std::string, Plan> m_plans               = {};
  std::multimap<double, Plan>           m_best_plans          = {};
  bool                                  m_is_valid            = false;
  wxCriticalSection                     m_critical_best_plans = {};
  wxCriticalSection                     m_critical_plans      = {};
  wxCriticalSection                     m_critical_is_valid   = {};

 protected:
  auto Entry()  //
      -> wxThread::ExitCode override;
  void StartThread();

 public:
  PlanView(wxWindow* parent, History candles, std::vector<Plan> plans);
  auto Get(const std::string& name)  //
      -> Plan;
  auto SelectedPlan()  //
      -> Plan;
  auto SelectedStrategy()  //
      -> std::unique_ptr<StrategyBase>;
  auto GetStrategy(const std::string& name)  //
      -> std::unique_ptr<StrategyBase>;
  void UpdateList();
  void Add(const Plan& plan);
  void Add(const std::vector<Plan>& plans);
  void Clear();
  void UpdateStrategies();
  void OnPlanSelection(wxCommandEvent& event);
  void OnStrategySelection(wxCommandEvent& event);
  void OnEntrySelection(wxCommandEvent& event);
  void OnExitSelection(wxCommandEvent& event);
  void PostUpdate();
  void PostNotify(const wxString& s);
  void Reset(const History& candles, std::vector<Plan> plans);

  void OnThreadUpdate(wxThreadEvent& evt);
  void OnDestroy(wxWindowDestroyEvent& evt);
};