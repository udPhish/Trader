#pragma once

#include <array>
#include <chrono>
#include <filesystem>
#include <fstream>
#include <sstream>

#include "date.h"

#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/tmpdir.hpp>

#include "Logger.h"

// std::array<std::array<double, 4>, 4> OrthoProjection(double left, double
// right,
//                                                     double bottom, double
//                                                     top, double near, double
//                                                     far);
//
// std::array<std::array<double, 4>, 4> PerspectiveProjection(
//    double left, double right, double bottom, double top, double near,
//    double far);

template<class T>
void Save(const T& obj, std::filesystem::path path)
{
  std::ofstream ofs;
  ofs.open(path);
  boost::archive::text_oarchive oa(ofs);
  oa << obj;
}
template<class T>
void Save(const T& obj, std::string filename)
{
  Save(obj,
       std::filesystem::path(std::string(boost::archive::tmpdir()) + "\\"
                             + filename));
}
// template <class T>
// void Save(const T& obj) {
//  Save(obj, boost::archive::tmpdir());
//}
template<class T>
T Load(std::filesystem::path path)
{
  std::ifstream ifs;
  ifs.open(path);
  boost::archive::text_iarchive ia(ifs);
  T                             t;
  ia >> t;
  return t;
}
template<class T>
T Load(std::string filename)
{
  return Load<T>(std::filesystem::path(std::string(boost::archive::tmpdir())
                                       + "\\" + filename));
}
// template <class T>
// T Load() {
//  return Load<T>(boost::archive::tmpdir());
//}
std::vector<std::vector<std::size_t>> CombinationsIndexHelper(
    std::size_t sample_size,
    std::size_t length,
    std::size_t skip);
std::vector<std::vector<std::size_t>> CombinationsIndex(std::size_t sample_size,
                                                        std::size_t length);
template<class T>
std::vector<std::vector<T>> Combinations(const std::vector<T>& values,
                                         std::size_t           length)
{
  std::vector<std::vector<T>> ret;
  for (auto& comb : CombinationsIndex(values.size(), length))
  {
    std::vector<T> vec;
    for (auto& i : comb)
      vec.push_back(values[i]);
    ret.push_back(vec);
  }
  return ret;
}
template<class T>
std::vector<std::vector<std::unique_ptr<T>>> Combinations(
    const std::vector<std::unique_ptr<T>>& values,
    std::size_t                            length)
{
  std::vector<std::vector<std::unique_ptr<T>>> ret;
  for (auto& comb : CombinationsIndex(values.size(), length))
  {
    std::vector<std::unique_ptr<T>> vec;
    for (auto& i : comb)
    {
      vec.push_back(values[i]->Clone());
    }
    ret.push_back(std::move(vec));
  }
  return ret;
}
template<class _Scale>
std::string MillisecondsToString(std::chrono::milliseconds time)
{
  using namespace date;
  auto time_point = floor<_Scale>(
      std::chrono::time_point<std::chrono::system_clock>{} + time);
  std::stringstream ss;
  ss << time_point;
  return ss.str();
}