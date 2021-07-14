#include "Exchange.h"
void Exchange::AddAsset(const std::string& asset) {
  markets.insert({asset, Market(asset)});
}
void Exchange::AddAssets(const std::vector<std::string>& assets) {
  for (auto& asset : assets) {
    AddAsset(asset);
  }
}

