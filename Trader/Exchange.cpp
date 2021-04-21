#include "Exchange.h"
void Exchange::AddAsset(std::string asset) {
  markets.insert({asset, Market(asset)});
}
void Exchange::AddAssets(std::vector<std::string> assets) {
  for (auto& asset : assets) {
    AddAsset(asset);
  }
}

