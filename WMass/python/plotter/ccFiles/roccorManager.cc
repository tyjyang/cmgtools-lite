#ifndef FUNCTIONS_ROCCOR_H
#define FUNCTIONS_ROCCOR_H

#include "TROOT.h"
#include "TFile.h"
#include "TH2.h"
#include "TH3.h"
#include "TF1.h"
#include "TH2Poly.h"
#include "TRandom.h"
#include "TRandom3.h"
#include "TSpline.h"
#include "TCanvas.h"
#include "TGraphAsymmErrors.h"
#include "TLorentzVector.h"
#include "TEfficiency.h"
#include "TVector2.h"
#include "Math/GenVector/LorentzVector.h"
#include "Math/GenVector/PtEtaPhiM4D.h"

#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <cstdlib> //as stdlib.h      
#include <cstdio>
#include <cmath>
#include <array>
#include <string>
#include <vector>
#include <unordered_map>
#include <utility>
#include <algorithm>
#include <boost/algorithm/string/join.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/functional/hash.hpp>
#include <limits>
#include <map>

#include "defines.h"
#include "RoccoR.h"

typedef ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> > PtEtaPhiMVector;

RoccoR* rochester_preVFP = new RoccoR("./RochesterStuff/RoccoR2016aUL.txt");
RoccoR* rochester_postVFP = new RoccoR("./RochesterStuff/RoccoR2016bUL.txt");

Vec_f applyRochester(const DataType isData, const DataEra dtype, const Vec_f& pt, const Vec_f& eta, const Vec_f& phi, const Vec_i& ch, const Vec_i& gen_idx, const Vec_f& gen_pt, const Vec_i& nTrackerLayers) {

  RoccoR* rochester = (dtype >= F_postVFP) ? rochester_postVFP : rochester_preVFP; // will it work?
  
  unsigned int size = pt.size();
  Vec_f res(size, 0.0);
  
  // for data can pass random stuff to gen_idx, gen_pt, and nTrackerLayers  
  if (isData) {

    for (unsigned int i = 0; i < size; ++i) {
      res[i] = pt[i] * rochester->kScaleDT(ch[i], pt[i], eta[i], phi[i]);
    }
    
  } else {

    for (unsigned int i = 0; i < size; ++i) {      
      if(gen_idx.at(i) >= 0) res[i] = pt[i] * rochester->kSpreadMC(ch[i], pt[i], eta[i], phi[i], gen_pt.at(gen_idx.at(i)), 0, 0);
      else res[i] = pt[i] * rochester->kSmearMC(ch[i], pt[i], eta[i], phi[i], nTrackerLayers.at(i), gRandom->Rndm(), 0, 0);
    }        

  }

  return res;
  
}


Vec_f applyRochesterMC(const DataEra dtype, const Vec_f& pt, const Vec_f& eta, const Vec_f& phi, const Vec_i& ch, const Vec_i& gen_idx, const Vec_f& gen_pt, const Vec_i& nTrackerLayers) {

  RoccoR* rochester = (dtype >= F_postVFP) ? rochester_postVFP : rochester_preVFP; // will it work?
  
  unsigned int size = pt.size();
  Vec_f res(size, 0.0);
  
  for (unsigned int i = 0; i < size; ++i) {      
    if(gen_idx.at(i) >= 0) res[i] = pt[i] * rochester->kSpreadMC(ch[i], pt[i], eta[i], phi[i], gen_pt.at(gen_idx.at(i)), 0, 0);
    else res[i] = pt[i] * rochester->kSmearMC(ch[i], pt[i], eta[i], phi[i], nTrackerLayers.at(i), gRandom->Rndm(), 0, 0);
  }        

  return res;
  
}


Vec_f applyRochesterData(const DataEra dtype, const Vec_f& pt, const Vec_f& eta, const Vec_f& phi, const Vec_i& ch) {

  RoccoR* rochester = (dtype >= F_postVFP) ? rochester_postVFP : rochester_preVFP; // will it work?
  
  unsigned int size = pt.size();
  Vec_f res(size, 0.0);
  
  // for data can pass random stuff to gen_idx, gen_pt, and nTrackerLayers  
  for (unsigned int i = 0; i < size; ++i) {
    res[i] = pt[i] * rochester->kScaleDT(ch[i], pt[i], eta[i], phi[i]);
  }
  
  return res;
  
}


#endif
