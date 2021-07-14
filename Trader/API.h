#pragma once
#include <sstream>

#include "boost/asio.hpp"
#include "boost/asio/ssl.hpp"
#include "boost/asio/ssl/context.hpp"
#include "boost/asio/ssl/stream.hpp"

#include "boost/beast.hpp"

#include "openssl/ssl.h"

#include "wx/wx.h"
//#include "openssl/hmac.h"

#include "Logger.h"

class API
{
  const std::string m_host;
  const std::string m_port;
  const std::string m_pk;

  boost::asio::io_context                                m_io_context;
  boost::asio::ip::tcp::resolver                         m_resolver;
  boost::asio::ssl::context                              m_ssl_context;
  boost::asio::ssl::stream<boost::asio::ip::tcp::socket> m_ssl_stream;

 public:
  API(std::string host, std::string m_port);
  void Connect();
  void Disconnect();
  auto SyncPost(std::string                                         target,
                boost::beast::http::verb                            action,
                const std::string&                                  query         = "",
                std::string                                         body          = "",
                const std::unordered_map<std::string, std::string>& header_fields = std::unordered_map<std::string, std::string>()) -> boost::beast::http::response<boost::beast::http::string_body>;
};
