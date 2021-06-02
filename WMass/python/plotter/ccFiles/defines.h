#ifndef DEFINES_H
#define DEFINES_H

#include <ROOT/RVec.hxx>
#include <ROOT/RVec.hxx>
#include <ROOT/RDataFrame.hxx>
#include <ROOT/RDF/RInterface.hxx>

using Vec_b = ROOT::VecOps::RVec<bool>;
using Vec_d = ROOT::VecOps::RVec<double>;
using Vec_f = ROOT::VecOps::RVec<float>;
using Vec_i = ROOT::VecOps::RVec<int>;
using Vec_ui = ROOT::VecOps::RVec<unsigned int>;

typedef enum {isoTrigPlus=0, isoTrigMinus, isoNotrig, noisoTrigPlus, noisoTrigMinus, noisoNotrig, antiisoTrigPlus, antiisoTrigMinus, antiisoNotrig, trackingAndReco} ScaleFactorType; // to use product
typedef enum {BToH=0, BToF, GToH, B, C, D, E, F, F_preVFP, F_postVFP, G, H} DataEra;
typedef enum {MC=0, Data} DataType;

bool isOddEvent(ULong64_t evt) {

  return (evt%2) ? 1 : 0;       

}

bool isEvenEvent(ULong64_t evt) {

  return (evt%2) ? 0 : 1;       

}


#endif
