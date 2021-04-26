#include "Application.h"

wxIMPLEMENT_APP(Application);

bool Application::OnInit() {
  m_main_frame = new Main();
  m_main_frame->Show(true);
  return true;
}