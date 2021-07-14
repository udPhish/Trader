#pragma once

#include <unordered_map>
#include <filesystem>

#include <boost/json.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/serialization/unordered_map.hpp>

#include <wx/wx.h>

#include "Market.h"
struct Exchange {
  std::unordered_map<std::string, Market> markets;

  void AddAsset(const std::string& asset);
  void AddAssets(const std::vector<std::string>& assets);

  friend class boost::serialization::access;
  template <class Archive>
  void serialize(Archive& ar, const unsigned int version) {
    ar& markets;
  }
};
