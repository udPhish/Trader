#include "Application.h"

wxIMPLEMENT_APP(Application);

auto Application::OnInit() -> bool
{
  m_main_frame = new Main();
  m_main_frame->Show(true);
  return true;
}
void Application::CleanUp()
{
  wxApp::CleanUp();
}