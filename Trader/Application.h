#pragma once
#define WIN32_LEAN_AND_MEAN

#include "wx/wx.h"

#include "Main.h"

class Application : public wxApp
{
 private:
  Main* m_main_frame = nullptr;

 public:
  virtual bool OnInit() wxOVERRIDE;
  virtual void CleanUp() override;
};