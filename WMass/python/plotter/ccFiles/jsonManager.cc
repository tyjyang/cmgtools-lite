#ifndef JSON_MANAGER_H
#define JSON_MANAGER_H

#include <iostream>
#include <string>
#include <vector>
#include <unordered_map>
#include <utility>      // std::pair
#include "defines.h"

#include "TROOT.h"

using namespace std;

bool isGoodRunLS(const DataType isData, const UInt_t run, const UInt_t lumi) {
  
  // for MC this function always returns true
  if (not isData) return true;

  if ( jsonMap_all.find(run) == jsonMap_all.end() ) return false; // run not found

  auto& validlumis = jsonMap_all.at(run);
  auto match = std::lower_bound(std::begin(validlumis), std::end(validlumis), lumi, 
          [](std::pair<unsigned int, unsigned int>& range, unsigned int val) { return range.second < val; });
  return match->first <= lumi && match->second >= lumi;
}

//==========================================================

#endif
