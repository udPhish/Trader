#include "Logger.h"
Logger::Logger(wxWindow* parent)
    : wxTextCtrl(parent,
                 wxID_ANY,
                 wxEmptyString,
                 wxDefaultPosition,
                 wxDefaultSize,
                 wxTE_MULTILINE | wxTE_READONLY)
    , m_stream(this)
    , m_log(new wxLogStream(&m_stream))
    , m_log_chain(m_log)
{
  SetBackgroundColour(wxColor("black"));
  SetForegroundColour(wxColor("white"));

  SetMinSize(this->GetFont().GetPixelSize() * 8);

  m_log->SetLogLevel(wxLOG_User);
}