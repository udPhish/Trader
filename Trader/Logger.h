#pragma once
#include <ostream>

#include <wx/log.h>
#include <wx/textctrl.h>
#include <wx/wx.h>

class Logger : public wxTextCtrl
{
 public:
  enum Level
  {
    FatalError   = wxLOG_FatalError,
    Error        = wxLOG_Error,
    Warning      = wxLOG_Warning,
    Notification = wxLOG_Message,
    Status       = wxLOG_Status,
    Message      = wxLOG_User
  };
  std::ostream m_stream;
  wxLogStream* m_log;
  wxLogChain   m_log_chain;

 public:
  Logger(wxWindow* parent);

  template<Level = Level::Message>
  static void Log(wxString message);
  template<>
  static void Log<Level::FatalError>(wxString message)
  {
    wxLogFatalError(message);
  }
  template<>
  static void Log<Level::Error>(wxString message)
  {
    wxLogError(message);
  }
  template<>
  static void Log<Level::Warning>(wxString message)
  {
    wxLogWarning(message);
  }
  template<>
  static void Log<Level::Notification>(wxString message)
  {
    wxLogMessage(message);
  }
  template<>
  static void Log<Level::Status>(wxString message)
  {
    wxLogStatus(message);
  }
  template<>
  static void Log<Level::Message>(wxString message)
  {
    wxLogGeneric(Level::Message, message);
  }
};
