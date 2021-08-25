#include "PlanView.h"

PlanView::PlanView(wxWindow* parent, History candles, std::vector<Plan> plans)
    : wxWindow(parent, wxID_ANY)
    , m_candles{std::move(candles)}
    , m_list_plans{new wxListBox{this,
                                 wxID_ANY,
                                 wxDefaultPosition,
                                 wxDefaultSize,
                                 0,
                                 nullptr,
                                 wxLB_SINGLE}}
    , m_is_valid{true}
    , m_plan_entries{new wxListBox{this,
                                   wxID_ANY,
                                   wxDefaultPosition,
                                   wxDefaultSize,
                                   0,
                                   nullptr,
                                   wxLB_SINGLE | wxLB_SORT}}
    , m_plan_exits{new wxListBox{this,
                                 wxID_ANY,
                                 wxDefaultPosition,
                                 wxDefaultSize,
                                 0,
                                 nullptr,
                                 wxLB_SINGLE | wxLB_SORT}}
{
  m_list_plans->SetMinSize({20, 20});
  m_plan_entries->SetMinSize({10, 10});
  m_plan_exits->SetMinSize({10, 10});
  Add(plans);
  auto* s  = new wxBoxSizer(wxVERTICAL);
  auto* s2 = new wxBoxSizer(wxHORIZONTAL);
  s->Add(m_list_plans, wxSizerFlags().Expand().Proportion(1));
  s->Add(s2, wxSizerFlags().Expand().Proportion(1));
  s2->Add(m_plan_entries, wxSizerFlags().Expand().Proportion(1));
  s2->Add(m_plan_exits, wxSizerFlags().Expand().Proportion(1));
  SetSizer(s);
  SetAutoLayout(true);
  m_list_plans->Bind(wxEVT_LISTBOX, &PlanView::OnPlanSelection, this);
  m_plan_entries->Bind(wxEVT_LISTBOX, &PlanView::OnStrategySelection, this);
  m_plan_exits->Bind(wxEVT_LISTBOX, &PlanView::OnStrategySelection, this);
  m_plan_entries->Bind(wxEVT_LISTBOX, &PlanView::OnEntrySelection, this);
  m_plan_exits->Bind(wxEVT_LISTBOX, &PlanView::OnExitSelection, this);
  Bind(wxEVT_THREAD, &PlanView::OnThreadUpdate, this);
  Bind(wxEVT_DESTROY, &PlanView::OnDestroy, this);
  StartThread();
}

void PlanView::Reset(const History& candles, std::vector<Plan> plans)
{
  if (GetThread() && GetThread()->IsRunning())
  {
    GetThread()->Delete();
  }
  m_candles = candles;
  m_best_plans.clear();
  Add(plans);
  StartThread();
}
auto PlanView::SelectedPlan()  //
    -> Plan
{
  std::string plan_selection = m_list_plans->GetStringSelection().ToStdString();
  if (plan_selection == "")
  {
    Logger::Log<Logger::Level::FatalError>("Invalid plan selection");
  }
  return Get(plan_selection);
}
auto PlanView::SelectedStrategy()  //
    -> std::unique_ptr<StrategyBase>
{
  std::string selection = m_plan_entries->GetStringSelection().ToStdString();
  if (selection == "")
  {
    selection = m_plan_exits->GetStringSelection().ToStdString();
  }
  if (selection == "")
  {
    Logger::Log<Logger::Level::FatalError>("Invalid strategy selection");
  }
  return GetStrategy(selection);
}
void PlanView::UpdateList()
{
  m_list_plans->Clear();
  {
    wxCriticalSectionLocker lock(m_critical_best_plans);
    for (auto& [key, value] : m_best_plans)
    {
      m_list_plans->Insert(value.Name(),0);
    }
  }
  if (!m_list_plans->IsEmpty())
  {
    m_list_plans->Select(0);
    UpdateStrategies();
    m_plan_entries->Select(0);
  }
}
void PlanView::Add(const Plan& plan)
{
  {
    wxCriticalSectionLocker lock(m_critical_best_plans);
    if (m_best_plans.size() < 10)
    {
      m_best_plans.insert({plan.AssessFinalTotal(m_candles), plan});
      UpdateList();
      return;
    }
  }
  wxCriticalSectionLocker lock(m_critical_plans);
  m_plans.insert({plan.Name(), plan});
}
void PlanView::Add(const std::vector<Plan>& plans)
{
  for (auto& plan : plans)
  {
    Add(plan);
  }
}
auto PlanView::Get(const std::string& name)  //
    -> Plan
{
  wxCriticalSectionLocker lock(m_critical_best_plans);
  for (auto& [key, val] : m_best_plans)
  {
    if (val.Name() == name)
    {
      return val;
    }
  }

  throw std::exception();
}
auto PlanView::GetStrategy(const std::string& name)
    -> std::unique_ptr<StrategyBase>
{
  auto selected_plan = SelectedPlan();
  for (auto& strategy : selected_plan.entries())
  {
    if (strategy->Name() == name)
    {
      return strategy->Clone();
    }
  }
  for (auto& strategy : selected_plan.exits())
  {
    if (strategy->Name() == name)
    {
      return strategy->Clone();
    }
  }
  return (*selected_plan.entries().begin())->Clone();
}
void PlanView::Clear()
{
  {
    wxCriticalSectionLocker lock(m_critical_plans);
    m_plans.clear();
  }
  {
    wxCriticalSectionLocker lock(m_critical_best_plans);
    m_best_plans.clear();
  }
  m_list_plans->Clear();
}

void PlanView::UpdateStrategies()
{
  auto plan = SelectedPlan();
  m_plan_entries->Clear();
  m_plan_exits->Clear();

  for (auto& entry : plan.entries())
  {
    m_plan_entries->Append(entry->Name());
  }
  for (auto& exit : plan.exits())
  {
    m_plan_exits->Append(exit->Name());
  }
  m_plan_entries->Select(0);
}
void PlanView::OnPlanSelection(wxCommandEvent& event)
{
  UpdateStrategies();
  event.Skip();
}
void PlanView::OnStrategySelection(wxCommandEvent& event)
{
  event.Skip();
}
void PlanView::OnEntrySelection(wxCommandEvent& event)
{
  m_plan_exits->Select(wxNOT_FOUND);
  event.Skip();
}
void PlanView::OnExitSelection(wxCommandEvent& event)
{
  m_plan_entries->Select(wxNOT_FOUND);
  event.Skip();
}
void PlanView::StartThread()
{
  if (CreateThread(wxTHREAD_JOINABLE) != wxTHREAD_NO_ERROR)
  {
    Logger::Log<Logger::Level::FatalError>(
        "Could not create the worker thread!");
    return;
  }
  if (GetThread()->Run() != wxTHREAD_NO_ERROR)
  {
    Logger::Log<Logger::Level::FatalError>("Could not run the worker thread!");
    return;
  }
}
void PlanView::PostUpdate()
{
  auto* ev = new wxThreadEvent();
  ev->SetInt(ThreadEvent::Update);
  wxQueueEvent(this, ev);
}
void PlanView::PostNotify(const wxString& s)
{
  auto* ev = new wxThreadEvent();
  ev->SetInt(ThreadEvent::Notify);
  ev->SetString(s);
  wxQueueEvent(this, ev);
}
auto PlanView::Entry() -> wxThread::ExitCode
{
  PostNotify("Thread Start");
  const int                             extraction_limit = 100;
  bool                                  should_notify    = false;
  double                                assess_total     = 0;
  std::unordered_map<std::string, Plan> plans;
  while (!GetThread()->TestDestroy())
  {
    GetThread()->Yield();
    wxCriticalSectionLocker lock(m_critical_is_valid);
    if (m_is_valid)
    {
      should_notify = false;
      if (plans.empty())
      {
        wxCriticalSectionLocker lock(m_critical_plans);
        if (m_plans.empty())
        {
          break;
        }
        PostNotify(wxString("Loading Batch (Remaining Plans: ")
                   << m_plans.size() << ")");
        for (std::size_t i = 0; i < extraction_limit; ++i)
        {
          plans.insert(m_plans.extract(m_plans.begin()));
          if (m_plans.empty())
          {
            break;
          }
        }
      }
      else
      {
        auto plan_node = plans.extract(plans.begin());
        // PostNotify(wxString("Testing: ") << plan_node.mapped().Name());
        assess_total = plan_node.mapped().AssessFinalTotal(m_candles);
        {
          wxCriticalSectionLocker lock(m_critical_best_plans);
          if (assess_total > m_best_plans.begin()->first)
          {
            m_best_plans.insert({assess_total, plan_node.mapped()});
            m_best_plans.erase(m_best_plans.begin());
            should_notify = true;
          }
        }
        // PostNotify(wxString("Assessed: ") << assess_total);
      }
      if (should_notify)
      {
        PostNotify(wxString("Found: ") << assess_total);
        PostUpdate();
      }
    }
  }
  PostNotify("Thread End");
  return (wxThread::ExitCode) nullptr;
}
void PlanView::OnThreadUpdate(wxThreadEvent& evt)
{
  switch (evt.GetInt())
  {
    case ThreadEvent::Notify:
      Logger::Log(evt.GetString());
      break;
    case ThreadEvent::Update:
      UpdateList();
      break;
  }
}
void PlanView::OnDestroy(wxWindowDestroyEvent& evt)
{
  {
    wxCriticalSectionLocker lock(m_critical_is_valid);
    m_is_valid = false;
  }

  if (GetThread() && GetThread()->IsRunning())
  {
    GetThread()->Delete();
  }
  evt.Skip();
}