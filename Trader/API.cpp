#include "API.h"

void API::Connect()
{
  if (!SSL_set_tlsext_host_name(m_ssl_stream.native_handle(), m_host.c_str()))
  {
    Logger::Log<Logger::Level::FatalError>("Unable to Connect");
  }
  try
  {
    auto const results = m_resolver.resolve(m_host, m_port);
    boost::asio::connect(m_ssl_stream.next_layer(), results);
    m_ssl_stream.handshake(boost::asio::ssl::stream_base::client);
  }
  catch (std::exception& e)
  {
    Logger::Log<Logger::Level::FatalError>("Unable to Connect");
  }
}
void API::Disconnect()
{
  try
  {
    m_ssl_stream.lowest_layer().shutdown(
        boost::asio::socket_base::shutdown_both);
    m_ssl_stream.lowest_layer().close();
  }
  catch (std::exception& e)
  {
    Logger::Log<Logger::Level::FatalError>("Unable to Disconnect");
  }
}
auto API::SyncPost(
    std::string                                         target,
    boost::beast::http::verb                            action,
    const std::string&                                  query,
    std::string                                         body,
    const std::unordered_map<std::string, std::string>& header_fields)
    -> boost::beast::http::response<boost::beast::http::string_body>
{
  boost::beast::http::request<boost::beast::http::string_body> request;

  target = "/" + target;
  if (query.length() > 0)
  {
    target += "?" + query;
  }

  request.method(action);
  request.target(target);
  request.version(11);
  request.set(boost::beast::http::field::host, m_host);
  if (body.length() > 0)
  {
    request.body() = std::move(body);
    request.set(boost::beast::http::field::content_length,
                std::to_string(request.body().length()));
  }
  for (auto& pair : header_fields)
  {
    request.insert(pair.first, pair.second);
  }

  boost::beast::http::write(m_ssl_stream, request);

  boost::beast::flat_buffer                                     buffer;
  boost::beast::http::response<boost::beast::http::string_body> response;

  boost::beast::http::read(m_ssl_stream, buffer, response);

  return response;
}

API::API(std::string host, std::string port)
    : m_host{std::move(host)}
    , m_port{std::move(port)}
    , m_resolver{m_io_context}
    , m_ssl_context{boost::asio::ssl::context::sslv23_client}
    , m_ssl_stream{m_io_context, m_ssl_context}
{
}