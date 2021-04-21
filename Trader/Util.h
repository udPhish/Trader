#pragma once

#include <array>
#include <filesystem>
#include <fstream>

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

template <class T>
void Save(const T& obj, std::filesystem::path path) {
  std::ofstream ofs;
  ofs.open(path);
  boost::archive::text_oarchive oa(ofs);
  oa << obj;
}
template <class T>
void Save(const T& obj, std::string filename) {
  Save(obj, std::filesystem::path(std::string(boost::archive::tmpdir()) + "\\" +
                                  filename));
}
// template <class T>
// void Save(const T& obj) {
//  Save(obj, boost::archive::tmpdir());
//}
template <class T>
T Load(std::filesystem::path path) {
  std::ifstream ifs;
  ifs.open(path);
  boost::archive::text_iarchive ia(ifs);
  T t;
  ia >> t;
  return t;
}
template <class T>
T Load(std::string filename) {
  return Load<T>(std::filesystem::path(std::string(boost::archive::tmpdir()) +
                                       "\\" + filename));
}
// template <class T>
// T Load() {
//  return Load<T>(boost::archive::tmpdir());
//}
template <class T>
std::vector<std::vector<T>> Combinations(std::vector<T> values,
                                         std::size_t length) {
  std::vector<std::vector<T>> ret;
  std::vector<std::vector<T>> combs;
  if (length == 1) {
    for(auto& value: values){
      std::vector<T> v;
      v.push_back(value);
      ret.push_back(v);
    }
    return ret;
  }
  for (std::size_t i = 0; i < values.size() - length; ++i) {
    combs = Combinations(std::vector<T>{values.begin() + 1, values.end()},
                         length - 1);
    for (auto& c : combs) {
      std::vector<T> comb;
      comb.push_back(values[i]);
      comb.insert(comb.end(), c.begin(), c.end());
      ret.push_back(comb);
    }
  }
  return ret;
}