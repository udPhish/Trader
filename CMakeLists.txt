# CMAKE VERSION
cmake_minimum_required(VERSION 3.19.0)

set(UDPH_CMAKE_VERSION master)
include(cmake/udph.cmake)

PROJECT_CREATE(Trader)

PROJECT_LOAD_CONAN_PACKAGE(OpenSSL 1.1.1k)
PROJECT_LOAD_GIT_REPOSITORY(Boost "https://github.com/boostorg/boost.git" master)
PROJECT_LOAD_GIT_REPOSITORY(libudph "https://github.com/udPhish/libudph.git" master)
PROJECT_LOAD_GIT_REPOSITORY(wxwidgets "https://github.com/wxWidgets/wxWidgets.git" master)

target_create_exe(Trader)
TARGET_SET_DIR(Trader Trader)

target_append_src(
	Trader
	"API.cpp"
	"APIBinance.cpp"
	"Application.cpp"
	"Candle.cpp"
	"CandleView.cpp"
	"Colour.cpp"
	"EntityCandle.cpp"
	"Exchange.cpp"
	"ExchangeView.cpp"
	"Indicator.cpp"
	"Logger.cpp"
	"Main.cpp"
	"Market.cpp"
	"OpenGLCanvas.cpp"
	"OpenGLFunctions.cpp"
	"PlanStatsView.cpp"
	"PlanView.cpp"
	"Strategies.cpp"
	"Strategy.cpp"
	"Util.cpp"
	"World.cpp"
)
target_append_hdr(
	Trader
	"API.h"
	"APIBinance.h"
	"Application.h"
	"Candle.h"
	"CandleView.h"
	"Colour.h"
	"EntityCandle.h"
	"Exchange.h"
	"ExchangeView.h"
	"Indicator.h"
	"Logger.h"
	"Main.h"
	"Market.h"
	"OpenGLCanvas.h"
	"OpenGLFunctions.h"
	"PlanStatsView.h"
	"PlanView.h"
	"Strategies.h"
	"Strategy.h"
	"Util.h"
	"World.h"
	"Boost_Serialization_chrono.h"
	"date.h"
)

TARGET_APPEND_DEP(Trader OpenSSL::SSL)
TARGET_APPEND_DEP(Trader wxcore)
TARGET_APPEND_DEP(Trader libudph::libudph)
TARGET_APPEND_DEP(Trader Boost::beast)
TARGET_APPEND_DEP(Trader Boost::asio)
TARGET_APPEND_DEP(Trader Boost::serialization)
TARGET_APPEND_DEP(Trader Boost::json)

target_link(Trader)

target_package(Trader)

FINALIZE()