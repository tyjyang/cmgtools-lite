#ifndef FUNCTIONS_WMASS_H
#define FUNCTIONS_WMASS_H

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
#include <iostream>
#include <algorithm>
#include <boost/algorithm/string/join.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/functional/hash.hpp>
#include "defines.h"

TF1 * helicityFractionSimple_0 = new TF1("helicityFraction_0", "3./4*(TMath::Sqrt(1-x*x))^2", -1., 1.);
TF1 * helicityFractionSimple_L = new TF1("helicityFraction_L", "3./8.*(1-x)^2"              , -1., 1.);
TF1 * helicityFractionSimple_R = new TF1("helicityFraction_R", "3./8.*(1+x)^2"              , -1., 1.);

TFile *_file_helicityFractionsSimple = NULL;
TH2 * helicityFractionsSimple_0 = NULL;
TH2 * helicityFractionsSimple_L = NULL;
TH2 * helicityFractionsSimple_R = NULL;

string getEnvironmentVariable(const string& env_var_name = "CMSSW_BASE") {

  char* _env_var_ptr = getenv(env_var_name.c_str());
  if (_env_var_ptr == nullptr) {
    cout << "Error: environment variable " << env_var_name << " not found. Exit" << endl;
    exit(EXIT_FAILURE);
  } else {
    string str = string(_env_var_ptr);
    return str;
  }

}

bool etaptIsInBin(const TH2& h,
		  const int& ietabin, const int& iptbin,
		  const float& eta, const float& pt) {
  if (pt >= h.GetYaxis()->GetBinLowEdge(iptbin) and pt < h.GetYaxis()->GetBinUpEdge(iptbin) and \
      eta >= h.GetXaxis()->GetBinLowEdge(ietabin) and eta < h.GetXaxis()->GetBinUpEdge(ietabin)
      )
    return true;
  else
    return false;

}
      
	   
float getValFromTH2(const TH2& h, const float& x, const float& y, const float& sumError=0.0) {
  //std::cout << "x,y --> " << x << "," << y << std::endl;
  int xbin = std::max(1, std::min(h.GetNbinsX(), h.GetXaxis()->FindFixBin(x)));
  int ybin  = std::max(1, std::min(h.GetNbinsY(), h.GetYaxis()->FindFixBin(y)));
  //std::cout << "xbin,ybin --> " << xbin << "," << ybin << std::endl;
  if (sumError)
    return h.GetBinContent(xbin, ybin) + sumError * h.GetBinError(xbin, ybin);
  else
    return h.GetBinContent(xbin, ybin);
}

float getValFromTH2bin(const TH2& h, const int& xbin, const int& ybin, const float& sumError=0.0) {
  if (sumError)
    return h.GetBinContent(xbin, ybin) + sumError * h.GetBinError(xbin, ybin);
  else
    return h.GetBinContent(xbin, ybin);
}

float getRelUncertaintyFromTH2(const TH2& h, const float& x, const float& y, const float valBadRatio = 1.0) {
  //std::cout << "x,y --> " << x << "," << y << std::endl;
  int xbin = std::max(1, std::min(h.GetNbinsX(), h.GetXaxis()->FindFixBin(x)));
  int ybin  = std::max(1, std::min(h.GetNbinsY(), h.GetYaxis()->FindFixBin(y)));
  //std::cout << "xbin,ybin --> " << xbin << "," << ybin << std::endl;
  if (h.GetBinContent(xbin, ybin) != 0.0)
    return h.GetBinError(xbin, ybin)/h.GetBinContent(xbin, ybin);
  else
    return valBadRatio;
}

float getRelUncertaintyFromTH2bin(const TH2& h, const int& xbin, const int& ybin, const float valBadRatio = 1.0) {
  if (h.GetBinContent(xbin, ybin) != 0.0)
    return h.GetBinError(xbin, ybin)/h.GetBinContent(xbin, ybin);
  else
    return valBadRatio;
}

float getAbsUncertaintyFromTH2(const TH2& h, const float& x, const float& y) {
  //std::cout << "x,y --> " << x << "," << y << std::endl;
  int xbin = std::max(1, std::min(h.GetNbinsX(), h.GetXaxis()->FindFixBin(x)));
  int ybin  = std::max(1, std::min(h.GetNbinsY(), h.GetYaxis()->FindFixBin(y)));
  //std::cout << "xbin,ybin --> " << xbin << "," << ybin << std::endl;
  return h.GetBinError(xbin, ybin);
}

float getAbsUncertaintyFromTH2bin(const TH2& h, const int& xbin, const int& ybin) {
  return h.GetBinError(xbin, ybin);
}


float returnChargeVal(float val1, int ch1, float val2, int ch2, ULong64_t evt){

    float retVal = -999.;

    if (evt%2 ) retVal = ch1 > 0 ? val1 : val2; //odd event numbers 
    else        retVal = ch1 > 0 ? val2 : val1; //even event numbers 

    return retVal;

}


float returnPlusVal(float val1, int ch1, float val2, int ch2){
  
  // return value of the lepton with desired charge, without looking at event parity
  // assumes opposite charges
  return (ch1 > 0) ? val1 : val2;

}

float returnMinusVal(float val1, int ch1, float val2, int ch2){
  
  // return value of the lepton with desired charge, without looking at event parity
  // assumes opposite charges
  return (ch1 < 0) ? val1 : val2;

}

float returnChargeValAllEvt(int desiredCharge, float val1, int ch1, float val2, int ch2, float returnForSameCharge = 0.0){
  
  if (ch1*ch2 > 0) {
    return returnForSameCharge;
  } else {
    if (desiredCharge > 0) return (ch1 > 0) ? val1 : val2;
    else                   return (ch1 < 0) ? val1 : val2;
  }

}


float mt_wlike(float pt1, float phi1, int ch1, float pt2, float phi2, int ch2, float met, float phimet, ULong64_t evt) {
  
  //if (ch1 == ch2) return mt_wlike_samesign(pt1,phi1,pt2,phi2,met,phimet);
  float ptL  = 0.0;
  float phiL = 0.0;

  float ptl  = 0.0;
  float phil = 0.0;

  // positive (negative) leptons on odd (even) events
  if (isOddEvent(evt)) {
    if (ch1 > 0) {
      ptL = pt1;
      phiL = phi1;
      ptl = pt2;
      phil = phi2;
    } else {
      ptL = pt2;
      phiL = phi2;
      ptl = pt1;
      phil = phi1;
    }
  } else {
    if (ch1 < 0) {
      ptL = pt1;
      phiL = phi1;
      ptl = pt2;
      phil = phi2;
    } else {
      ptL = pt2;
      phiL = phi2;
      ptl = pt1;
      phil = phi1;
    }
  }

  TVector2 pl = TVector2();
  pl.SetMagPhi(ptl,phil);

  TVector2 metv = TVector2();
  metv.SetMagPhi(met,phimet);
  TVector2 met_wlike = pl+metv;

  return std::sqrt(2*ptL*met_wlike.Mod()*(1-std::cos(phiL-met_wlike.Phi())));

}

float lepPlusMetPt(float pt, float phi, float met, float phimet) {

  TVector2 pl = TVector2();
  pl.SetMagPhi(pt, phi);

  TVector2 met_wlike = TVector2();
  met_wlike.SetMagPhi(met, phimet);
  met_wlike = pl + met_wlike;
  return met_wlike.Mod();
  
}

float mt_wlike_nano(float pt, float phi, float ptOther, float phiOther, float met, float phimet) {
  
  TVector2 pl = TVector2();
  pl.SetMagPhi(ptOther,phiOther);

  TVector2 met_wlike = TVector2();
  met_wlike.SetMagPhi(met,phimet);
  met_wlike = pl + met_wlike;

  return std::sqrt(2*pt*met_wlike.Mod()*(1-std::cos(phi-met_wlike.Phi())));

}

float helicityWeightSimple(float yw, float ptw, float costheta, int pol)
{

  if (!helicityFractionsSimple_0 || !helicityFractionsSimple_L || !helicityFractionsSimple_R) {
    _file_helicityFractionsSimple = new TFile("w-mass-13TeV/fractionReweighting/fractions.root","read");
    helicityFractionsSimple_0 = (TH2F*)(_file_helicityFractionsSimple->Get("fraction0_plus_sym"));
    helicityFractionsSimple_L = (TH2F*)(_file_helicityFractionsSimple->Get("fractionL_plus_sym"));
    helicityFractionsSimple_R = (TH2F*)(_file_helicityFractionsSimple->Get("fractionR_plus_sym"));
  }

  if (std::abs(costheta) > 1.) {
    //std::cout << " found an event with weird cosTheta = " << costheta << std::endl;
    //std::cout << " setting event weight to 0" << std::endl;
    return 0;
  }

  TH2 *hist_f0 = helicityFractionsSimple_0;
  TH2 *hist_fL = helicityFractionsSimple_L;
  TH2 *hist_fR = helicityFractionsSimple_R;

  // float yval  = std::abs(yw) > hist_f0->GetXaxis()->GetXmax() ? hist_f0->GetXaxis()->GetXmax() : yw;
  // float ptval = ptw > hist_f0->GetYaxis()->GetXmax() ? hist_f0->GetYaxis()->GetXmax() : ptw;

  int ywbin = std::max(1, std::min(hist_f0->GetNbinsX(), hist_f0->GetXaxis()->FindBin(yw )));
  int ptbin = std::max(1, std::min(hist_f0->GetNbinsY(), hist_f0->GetYaxis()->FindBin(ptw)));

  float f0 = hist_f0->GetBinContent(ywbin, ptbin);
  float fL = hist_fL->GetBinContent(ywbin, ptbin);
  float fR = hist_fR->GetBinContent(ywbin, ptbin);

  float f0Term = helicityFractionSimple_0->Eval(costheta);
  float fLTerm = helicityFractionSimple_L->Eval(costheta);
  float fRTerm = helicityFractionSimple_R->Eval(costheta);

  if      (pol == 0) return f0*f0Term/(f0*f0Term+fL*fLTerm+fR*fRTerm);
  else if (pol == 1) return fL*fLTerm/(f0*f0Term+fL*fLTerm+fR*fRTerm);
  else if (pol == 2) return fR*fRTerm/(f0*f0Term+fL*fLTerm+fR*fRTerm);
        
  std::cout << "something went wrong in the helicity reweighting" << std::endl;
  return -99999.;

}

float mydeltaPhi(float phi1, float phi2) {
  float result = phi1 - phi2;
  while (result > float(M_PI)) result -= float(2*M_PI);
  while (result <= -float(M_PI)) result += float(2*M_PI);
  return result;
}

float mydeltaR(float eta1, float phi1, float eta2, float phi2) {
  float deta = eta1-eta2;
  float dphi = mydeltaPhi(phi1,phi2);
  return std::sqrt(deta*deta + dphi*dphi);
}


//-------------------

TFile *_file_ratio_FSRoverNoFSR_etaPt_mu = NULL;
TH2D * ratio_FSRoverNoFSR_etaPt_mu = NULL;
TFile *_file_ratio_FSRoverNoFSR_etaPt_el = NULL;
TH2D * ratio_FSRoverNoFSR_etaPt_el = NULL;

// old function, should not be needed this time
float FSRscaleFactorEtaPt(int pdgId, float dresspt, float dresseta) {

  // these histograms are the gen-level xsec before any cut (except the gen_decayId)
  // the yields and ratio of fsr and no-fsr (both on top of W-pt reweighting) are here:
  // muon: http://mciprian.web.cern.ch/mciprian/wmass/13TeV/distribution/FSR_atGenLevel_muon_genBinAnalysis/
  // electron: http://mciprian.web.cern.ch/mciprian/wmass/13TeV/distribution/FSR_atGenLevel_electron_genBinAnalysis/
  // FSR should not change the gen-level xsec, but our fsr weights can artificially do it
  // so, we use the ratio to rescale the fsr weight in bins of gen-level pt and eta
  // the histogram uses the binning for the 2D xsec measurement (as each of them should keep the same xsec as without fsr)
 

  TH2D* hratio_FSR_noFSR = NULL;
  Double_t outlierWgt = 0.0;  // hardcoded below, taken from ratio of yields outside acceptance

  if (abs(pdgId)==11) {

    if (!ratio_FSRoverNoFSR_etaPt_el) {
      _file_ratio_FSRoverNoFSR_etaPt_el = new TFile("w-mass-13TeV/theoryReweighting/FSR_atGenLevel_electron_genEtaPtAnalysis.root","read");
      ratio_FSRoverNoFSR_etaPt_el = (TH2D*) _file_ratio_FSRoverNoFSR_etaPt_el->Get("ratio__Wzpt_el__Wfsr_el");
    }
    hratio_FSR_noFSR = ratio_FSRoverNoFSR_etaPt_el;
    outlierWgt = 0.964796;

  } else if (abs(pdgId)==13) {

     if (!ratio_FSRoverNoFSR_etaPt_mu) {
      _file_ratio_FSRoverNoFSR_etaPt_mu = new TFile("w-mass-13TeV/theoryReweighting/FSR_atGenLevel_muon_genEtaPtAnalysis.root","read");
      ratio_FSRoverNoFSR_etaPt_mu = (TH2D*) _file_ratio_FSRoverNoFSR_etaPt_mu->Get("ratio__Wzpt_mu__Wfsr_mu");
    }
    hratio_FSR_noFSR = ratio_FSRoverNoFSR_etaPt_mu;
    outlierWgt = 0.914322;

  } else {

    return 1.0;

  }

  int etabin = hratio_FSR_noFSR->GetXaxis()->FindFixBin(fabs(dresseta));
  int ptbin = hratio_FSR_noFSR->GetXaxis()->FindFixBin(fabs(dresspt));
  if (ptbin == 0 or ptbin > hratio_FSR_noFSR->GetNbinsY() or etabin == 0 or etabin> hratio_FSR_noFSR->GetNbinsX()) {
    return outlierWgt;
  } else {
    return 1./hratio_FSR_noFSR->GetBinContent(etabin,ptbin); // ratio is FSR/no-FSR, but need the opposite
  }

}

//-------------------

TFile *_file_fsrWeights_simple = NULL;
TH1F * fsrWeights_el = NULL;
TH1F * fsrWeights_mu = NULL;

float fsrPhotosWeightSimple(int pdgId, float dresspt, float barept, bool normToSameGenArea = false, float dresseta = 0) {

  if (!fsrWeights_el || !fsrWeights_mu) {
    _file_fsrWeights_simple = new TFile("w-mass-13TeV/theoryReweighting/photos_rwgt_integrated.root","read");
    fsrWeights_el  = (TH1F*)(_file_fsrWeights_simple->Get("w_e_h_lptBareOverDressed_ratio"));
    fsrWeights_mu  = (TH1F*)(_file_fsrWeights_simple->Get("w_mu_h_lptBareOverDressed_ratio"));
  }
  TH1F *fsrWeights = 0;
  if      (abs(pdgId)==11) fsrWeights = fsrWeights_el;
  else if (abs(pdgId)==13) fsrWeights = fsrWeights_mu;
  else return 1;

  int ratiobin  = std::max(1, std::min(fsrWeights->GetNbinsX(), fsrWeights->GetXaxis()->FindFixBin(barept/dresspt)));
  
  // normfactor is used to keep the total xsec unchanged when using the fsr
  // it was obtained from the gen level cross section for e/mu without any cut (except the gen-pdgid), as the ratio of nomi/fsr (nomi/fsr < 1)
  //Double_t normfactor = (abs(pdgId) == 11) ? 0.964792 : 0.914320; 
  Double_t normfactor = (normToSameGenArea) ? FSRscaleFactorEtaPt(pdgId, dresspt, dresseta) : 1.0; 
  return normfactor * fsrWeights->GetBinContent(ratiobin);

}


TFile *_file_fsrWeights = NULL;
TH3F * fsrWeights_elplus  = NULL;
TH3F * fsrWeights_elminus = NULL;
TH3F * fsrWeights_muplus  = NULL;
TH3F * fsrWeights_muminus = NULL;

float fsrPhotosWeight(int pdgId, float dresseta, float dresspt, float barept) {
  if (!fsrWeights_elplus || !fsrWeights_elminus || !fsrWeights_muplus || !fsrWeights_muminus) {
    _file_fsrWeights = new TFile("w-mass-13TeV/theoryReweighting/photos_rwgt.root","read");
    fsrWeights_elplus  = (TH3F*)(_file_fsrWeights->Get("qed_weights_wp_e"));
    fsrWeights_elminus = (TH3F*)(_file_fsrWeights->Get("qed_weights_wm_e"));
    fsrWeights_muplus  = (TH3F*)(_file_fsrWeights->Get("qed_weights_wp_mu"));
    fsrWeights_muminus = (TH3F*)(_file_fsrWeights->Get("qed_weights_wm_mu"));
  }
  TH3F *fsrWeights = 0;
  if      (abs(pdgId)==11) fsrWeights = ( pdgId>0 ? fsrWeights_elplus : fsrWeights_elminus );
  else if (abs(pdgId)==13) fsrWeights = ( pdgId>0 ? fsrWeights_muplus : fsrWeights_muminus );
  else return 1;

  int etabin = std::max(1, std::min(fsrWeights->GetNbinsX(), fsrWeights->GetXaxis()->FindFixBin(fabs(dresseta))));
  int ptbin  = std::max(1, std::min(fsrWeights->GetNbinsY(), fsrWeights->GetYaxis()->FindFixBin(dresspt)));
  int zbin  = std::max(1, std::min(fsrWeights->GetNbinsZ(), fsrWeights->GetZaxis()->FindFixBin(barept/dresspt)));

  return fsrWeights->GetBinContent(etabin,ptbin,zbin);
}

TFile *_file_dyptWeights = NULL;
TH1F *amcnlody = NULL; 

float dyptWeight(float pt2l, int isZ, bool scaleNormWToGenXsecBeforeCuts = false) {
  if (!amcnlody) {
    _file_dyptWeights = new TFile("w-mass-13TeV/theoryReweighting/zpt_weights.root");
    amcnlody = (TH1F*)(_file_dyptWeights->Get("amcnlo"));
  }
  int ptbin = std::max(1, std::min(amcnlody->GetNbinsX(), amcnlody->GetXaxis()->FindFixBin(pt2l)));
  // for the Z, change the pT *and* normalization to the measured one in data
  // for the W, change only the pT, but scale back to the total xsec of MC@NLO (factor comptued for total xsec in acceptance)
  // when using pt range 26-56 instead of 26-45, there is an additional factor 0.987, so the actual number would be 0.946 
  
  // when using only 0.958 I see that total gen-xsec no-Wpt over with-Wpt is 1.0138013 for W+ and 1.0144267 for W-
  // let's take only 1.014 without distinguishing charges
  //float scaleToMCaNLO = isZ ? 1. : 0.958;
  float scaleW = scaleNormWToGenXsecBeforeCuts ? 0.9714120 : 0.958;  //1.014 * 0.958 
  float scaleToMCaNLO = isZ ? 1. : scaleW;
  // plots are MC/data
  //float scaleToMCaNLO = isZ ? 1. : 0.958;
  return scaleToMCaNLO / amcnlody->GetBinContent(ptbin);
}
//=================
TFile *_file_dyptWeights_MGaMCAtNLO_PowhegMiNNLO = NULL;
TH1F *zptRatio_MGaMCAtNLO_PowhegMiNNLO = NULL; 

float dyptWeight_PowhegMiNNLO(float pt2l, int isZ, bool scaleNormWToGenXsecBeforeCuts = false, bool usePhotos = false) {
  if (!amcnlody) {
    _file_dyptWeights = new TFile("w-mass-13TeV/theoryReweighting/zpt_weights.root");
    amcnlody = (TH1F*)(_file_dyptWeights->Get("amcnlo"));
  }
  int ptbin = std::max(1, std::min(amcnlody->GetNbinsX(), amcnlody->GetXaxis()->FindFixBin(pt2l)));
  // for the Z, change the pT *and* normalization to the measured one in data
  // for the W, change only the pT, but scale back to the total xsec of MC@NLO (factor comptued for total xsec in acceptance)
  // when using pt range 26-56 instead of 26-45, there is an additional factor 0.987, so the actual number would be 0.946 
  
  // when using only 0.958 I see that total gen-xsec no-Wpt over with-Wpt is 1.0138013 for W+ and 1.0144267 for W-
  // let's take only 1.014 without distinguishing charges
  //float scaleToMCaNLO = isZ ? 1. : 0.958;

  if (!zptRatio_MGaMCAtNLO_PowhegMiNNLO) {
    // 2 GeV granularity from 0 to 100 GeV (last bin is overflow)
    _file_dyptWeights_MGaMCAtNLO_PowhegMiNNLO = new TFile("w-mass-13TeV/theoryReweighting/ZptRatio_MGaMCAtNLO_PowhegMiNNLO_dressed_noCuts.root");
    string hname_MGaMCAtNLO_PowhegMiNNLO = "pythia8";
    if (usePhotos) hname_MGaMCAtNLO_PowhegMiNNLO = "pythia8_photos";
    zptRatio_MGaMCAtNLO_PowhegMiNNLO = (TH1F*)(_file_dyptWeights_MGaMCAtNLO_PowhegMiNNLO->Get(hname_MGaMCAtNLO_PowhegMiNNLO.c_str()));
  }
  int ptbinCorr = std::max(1, std::min(zptRatio_MGaMCAtNLO_PowhegMiNNLO->GetNbinsX(), 
				       zptRatio_MGaMCAtNLO_PowhegMiNNLO->GetXaxis()->FindFixBin(pt2l))
			   );

  float scaleW = scaleNormWToGenXsecBeforeCuts ? 0.9714120 : 0.958;  //1.014 * 0.958 // to be revisited with new MC
  float scaleToMCaNLO = isZ ? 1. : scaleW;
  // plots are MC/data
  //float scaleToMCaNLO = isZ ? 1. : 0.958;
  return zptRatio_MGaMCAtNLO_PowhegMiNNLO->GetBinContent(ptbinCorr) * scaleToMCaNLO / amcnlody->GetBinContent(ptbin);
}


//=================

TFile *_file_postfitWeights = NULL;
TH1F *hist_scales_0plus  = NULL;
TH1F *hist_scales_Lplus  = NULL;
TH1F *hist_scales_Rplus  = NULL;
TH1F *hist_scales_0minus = NULL;
TH1F *hist_scales_Lminus = NULL;
TH1F *hist_scales_Rminus = NULL;

float postfitQCDWeight(float pt2l, int pol, int charge, int flav=0) {
  TH1F* hist_scales = NULL;
  if (!_file_postfitWeights) {
    if (flav==0)
      _file_postfitWeights = new TFile("w-helicity-13TeV/theoryReweighting/postfit_wgts_fitlep.root", "read");
    else if (flav==13)
      _file_postfitWeights = new TFile("w-helicity-13TeV/theoryReweighting/postfit_wgts_fitmu.root", "read");
    else if (flav==11)
      _file_postfitWeights = new TFile("w-helicity-13TeV/theoryReweighting/postfit_wgts_fitel.root","read");
    else {
      std::cout << "ERROR! Unknown flavor: " << flav << " Returning 0 weight." << std::endl;
      return 0;
    }
    hist_scales_0plus  = (TH1F*)(_file_postfitWeights->Get("weights_longplus"));
    hist_scales_Lplus  = (TH1F*)(_file_postfitWeights->Get("weights_leftplus"));
    hist_scales_Rplus  = (TH1F*)(_file_postfitWeights->Get("weights_rightplus"));
    hist_scales_0minus = (TH1F*)(_file_postfitWeights->Get("weights_longminus"));
    hist_scales_Lminus = (TH1F*)(_file_postfitWeights->Get("weights_leftminus"));
    hist_scales_Rminus = (TH1F*)(_file_postfitWeights->Get("weights_rightminus"));
  }
  if (charge>0) {
    if (pol==0)      hist_scales = hist_scales_0plus;
    else if (pol==1) hist_scales = hist_scales_Lplus;
    else if (pol==2) hist_scales = hist_scales_Rplus;
    else std::cerr << "ERROR: polarization " << pol << " not defined for postfitQCDWeight()" << std::endl;
  } else {
    if (pol==0)      hist_scales = hist_scales_0minus;
    else if (pol==1) hist_scales = hist_scales_Lminus;
    else if (pol==2) hist_scales = hist_scales_Rminus;
    else std::cerr << "ERROR: polarization " << pol << " not defined for postfitQCDWeight()" << std::endl;
  }
  int ptbin = std::max(1, std::min(hist_scales->GetNbinsX(), hist_scales->GetXaxis()->FindFixBin(pt2l)));
  return hist_scales->GetBinContent(ptbin);
}


//==================================================
TRandom3 *rng = NULL;

float getSmearedVar(float var, float smear, ULong64_t eventNumber, int isData, bool smearOnlyMC=false) {

  if (smearOnlyMC && isData) return var;

  if(!rng) rng = new TRandom3();
  // use eventNumber as seed, otherwise each time the function is called for the same event, the smearer produce a different smeared value
  rng->SetSeed(eventNumber); // make it really random across different jobs    
  return var * ( 1.0 + smear * rng->Gaus());

}

//================================================== 

// Sorry you have to manually keep these consistent
// let's remove BToH, we will never use it
// std::unordered_map<DataEra, std::string> eraNames = { {BToH, "BtoH"}, {BToF, "BtoF"}, {GToH, "GtoH"} };
std::unordered_map<DataEra, std::string> eraNames = { {BToF, "BtoF"}, {GToH, "GtoH"} };
std::unordered_map<DataType, std::string> datatypeNames = { {MC, "MC"}, {Data, "Data"} };
std::unordered_map<ScaleFactorType, std::string> scalefactorNames = { {isoTrigPlus, "isoTrigPlus"}, {isoTrigMinus, "isoTrigMinus"}, {isoNotrig, "isoNotrig"}, {antiisoTrigPlus, "antiisoTrigPlus"}, {antiisoTrigMinus, "antiisoTrigMinus"}, {antiisoNotrig, "antiisoNotrig"} };
  
struct pair_hash
{
    template <class T1, class T2>
    std::size_t operator() (const std::pair<T1, T2> &p) const
    {
      std::size_t seed = 0;
      boost::hash_combine(seed, p.first);
      boost::hash_combine(seed, p.second);
      return seed;
    }
};

std::unordered_map<DataEra, TH1D> hMuonPrefiring = {}; // will store pre and post (only BToF and GToH)
std::unordered_map<std::pair<ScaleFactorType, DataEra>,  TH2D, pair_hash> scaleFactorHist = {};
std::unordered_map<std::pair<ScaleFactorType, DataType>, TH2D, pair_hash> prePostCorrToHist = {};

void initializeScaleFactors(const string& _filename_allSF = "./testMuonSF/scaleFactorProduct_31Mar2021.root") {

  TFile _file_allSF = TFile(_filename_allSF.c_str(), "read");
  if (!_file_allSF.IsOpen())
      std::cerr << "WARNING: Failed to open scaleFactors file " << _filename_allSF << "! No scale factors will be applied\n";

  std::cout << "INFO >>> Initializing histograms for SF from file " << _filename_allSF << std::endl;

  for (auto& corr : scalefactorNames) {
    for (auto& era : eraNames) {
      std::vector<std::string> vars = {"fullSF2D", corr.second, era.second};
      std::string corrname = boost::algorithm::join(vars, "_");
      auto* histptr = static_cast<TH2D*>(_file_allSF.Get(corrname.c_str()));
      if (histptr == nullptr)
	std::cerr << "WARNING: Failed to load correction " << corrname << " in file "
		  << _filename_allSF << "! scale factors for this correction will be set to 1.0";
      histptr->SetDirectory(0);
      DataEra eraVal = era.first;
      ScaleFactorType key = corr.first;
      // std::cout << "Histogram key " << key << " and era " << era.second << std::endl;
      auto corrKey = std::make_pair(key, eraVal);
      scaleFactorHist[corrKey] = *static_cast<TH2D*>(histptr);
    }

    for (auto& datatype : datatypeNames) {
      std::vector<std::string> vars = {"fullSF2D", datatype.second, "preOverPost", corr.second};
      std::string corrname = boost::algorithm::join(vars, "_");
      auto* histptr = static_cast<TH2D*>(_file_allSF.Get(corrname.c_str()));
      if (histptr == nullptr)
	std::cerr << "WARNING: Failed to load correction " << corrname << " in file "
		  << _filename_allSF << "! scale factors for this correction will be set to 1.0";
      histptr->SetDirectory(0);
      DataType typeVal = datatype.first;
      ScaleFactorType key = corr.first;
      // std::cout << "Histogram key " << key << " and era " << era.second << std::endl;
      auto corrKey = std::make_pair(key, typeVal);
      prePostCorrToHist[corrKey] = *static_cast<TH2D*>(histptr);
    }
  
  }
  
  _file_allSF.Close(); // should work since we used TH1D::SetDirectory(0) to detach histogram from file

  std::string _filename_prefiring = "./testMuonSF/muonPrefiring_prePostVFP.root";
  TFile _file_prefiring = TFile(_filename_prefiring.c_str(), "read");
  if (!_file_prefiring.IsOpen())
    std::cerr << "WARNING: Failed to open prefiring file " << _filename_prefiring << "\n";
  std::cout << "INFO >>> Initializing histograms for prefiring from file " << _filename_prefiring << std::endl;
  hMuonPrefiring[BToF] = *(static_cast<TH1D*>(_file_prefiring.Get("muonPrefiring_preVFP")));
  hMuonPrefiring[BToF].SetDirectory(0);
  hMuonPrefiring[GToH] = *(static_cast<TH1D*>(_file_prefiring.Get("muonPrefiring_postVFP")));
  hMuonPrefiring[GToH].SetDirectory(0);
  _file_prefiring.Close();
    
}


float _get_MuonPrefiringSF_singleMuon(float eta, DataEra era = BToF) {

  // no need to care about under/overflow, the prefiring would be 0 there, and the actual range is -2.4000001, 2.4000001
  return 1.0 - hMuonPrefiring[era].GetBinContent(hMuonPrefiring[era].FindFixBin(eta));
  
}

// may just add this in a function, to avoid defining a column
Vec_b prefirableMuon(const Vec_f& pt, const Vec_b& looseId) {

  Vec_b res(pt.size(),false); // initialize to 0
  for (unsigned int i = 0; i < res.size(); ++i) {
    if (pt[i] < 22) continue;
    if (not looseId[i]) continue;
    res[i] = true;
  }
  return res;
  
}

float _get_MuonPrefiringSF(const Vec_f& eta, const Vec_f& pt, const Vec_b& looseId, DataEra era = BToF) {

  // can be called as Muon_eta, Muon_pt, Muon_looseId, no need to use Muon_eta[prefirableMuon]
  float sf = 1.0;
  // get SF = Prod_i( 1 - P_pref[i] )
  // int nBinsX = hMuonPrefiring[era].GetNbinsX(); // not needed if not neglecting under/overflow

  // std::cout << "PREFIRING FOR: " << eraNames[era] << std::endl;

  const TH1D& hprefire = hMuonPrefiring[era];
  for (unsigned int i = 0; i < eta.size(); ++i) {
    if (pt[i] < 22) continue;
    if (not looseId[i]) continue;
    // no need to care about under/overflow, the prefiring would be 0 there, and the actual range is -2.4000001, 2.4000001
    //sf *= (1.0 - hMuonPrefiring[era].GetBinContent(std::max(0, std::min(nBinsX, hMuonPrefiring[era].FindFixBin(eta[i])) )) );
    sf *= (1.0 - hprefire.GetBinContent(hprefire.FindFixBin(eta[i])));
  }
  return sf;  
  
}

Vec_f _get_MuonPrefiringSFvariation(int n_prefireBinNuisance,
				    const Vec_f& eta, const Vec_f& pt, const Vec_b& looseId,
				    DataEra era = BToF
				    ) {

  // this function directly provides the alternative prefiring SF, not the variation on the original one
  // it is supposed to be called instead of the nominal weight
  // it returns a vector used as an event weight to get all variations in the same TH3 (eta-pt-prefireBin)

  Vec_f res(n_prefireBinNuisance, 1.0); // initialize to 1

  int prefireBin = 0;
  float tmpval = 0.0;
  const TH1D& hprefire = hMuonPrefiring[era];

  for (unsigned int i = 0; i < eta.size(); ++i) {

    if (pt[i] < 22) continue;
    if (not looseId[i]) continue;
    prefireBin = hprefire.FindFixBin(eta[i]);
    // fill the vector with the nominal weight in each bin,
    // while the vector bin corresponding to prefireBin will be filled with prefiring probability moved by its uncertainty
    tmpval = res[prefireBin-1];
    res *= (1.0 - hprefire.GetBinContent(prefireBin));
    res[prefireBin-1] = tmpval * (1.0 - hprefire.GetBinContent(prefireBin) - hprefire.GetBinError(prefireBin));

  }
  
  return res;

}

float _get_fullMuonSF(float pt,      float eta,      int charge,
		      float ptOther, float etaOther,
		      DataEra era = BToF,
		      bool isoSF1 = true,
		      bool isoSF2 = true
		      ) {

  // function to get full muon scale factor for  analysis (except prefiring, handled elsewhere)
  // first three arguments are for the triggering muon, second two for the non triggering one
  // isoSF1 and isoSF2 are to use SF for isolation or antiisolation, for triggering and non triggering muons respectively
  // if ptOther < 0, the second lepton is ignored, so we can use this function for Wmass as well (etaOther is not used)
  // may actually use another function for Wmass
  
  //std::cout << "pt,eta       -> " << pt      << "," << eta      << std::endl;
  //std::cout << "pt,eta other -> " << ptOther << "," << etaOther << std::endl;
  // std::cout << "scale factors for " << eraNames[era] << std::endl;
  
  ScaleFactorType sftype = charge > 0 ? isoTrigPlus : isoTrigMinus;
  if (not isoSF1) {
    sftype = charge > 0 ? antiisoTrigPlus : antiisoTrigMinus;
  }

  auto const key = std::make_pair(sftype, era);
  const TH2D& hsf = scaleFactorHist.at(key);
  float sf = getValFromTH2(hsf, eta, pt);
  //std::cout << "scale factor main leg -> " << sf << std::endl;

  if (ptOther > 0.0) {
    sftype = isoSF2 ? isoNotrig : antiisoNotrig;
    auto const keyOther = std::make_pair(sftype, era);
    const TH2D& hsfOther = scaleFactorHist.at(keyOther);
    sf *= getValFromTH2(hsfOther, etaOther, ptOther);
  }
  //std::cout << "final scale factor -> " << sf << std::endl;
  return sf;
}

float _get_fullMuonSF_preOverPost(float pt,      float eta,      int charge,
				  float ptOther, float etaOther,
				  DataType dtype = MC,
				  bool isoSF1 = true,
				  bool isoSF2 = true
				  ) {

  //std::cout <<  "type " << datatypeNames[dtype] << std::endl;
  //std::cout << "pt,eta       -> " << pt      << "," << eta      << std::endl;
  //std::cout << "pt,eta other -> " << ptOther << "," << etaOther << std::endl;

  ScaleFactorType sftype = charge > 0 ? isoTrigPlus : isoTrigMinus;
  if (not isoSF1) {
    sftype = charge > 0 ? antiisoTrigPlus : antiisoTrigMinus;
  }

  auto const key = std::make_pair(sftype, dtype);
  const TH2D& hcorr = prePostCorrToHist.at(key);
  float sf = getValFromTH2(hcorr, eta, pt);
  //std::cout << "scale factor main leg -> " << sf << std::endl;

  if (ptOther > 0.0) {
    sftype = isoSF2 ? isoNotrig : antiisoNotrig;
    auto const keyOther = std::make_pair(sftype, dtype);
    const TH2D& hcorrOther = prePostCorrToHist.at(keyOther);
    sf *= getValFromTH2(hcorrOther, etaOther, ptOther);
  }

  return sf;

}


Vec_f _get_fullMuonSFvariation(int n_tnpBinNuisance,
			       float pt,      float eta, int charge,
			       float ptOther=-1, float etaOther=-1,
			       DataEra era = BToF,
			       bool isoSF1 = true,
			       bool isoSF2 = true
			       ) {

  // this is an helper function to define the effSystvariations for the Wlike analysis
  // idea is to fill again the alternative template for each effStat nuisance, where the
  // nuisance is defined for each single TnP eta-pt bin, and they will be considered as uncorrelated
  // so not as done in SMP-18-012 from the fit to efficiencies with Error function
  //
  // this should replace the nominal SF weight, and return SF for any bin except the specific one corresponding
  // to the nuisance parameter, and SF+err for that one
  // in order to facilitate the usage with the ReplaceWeight functionality to customize event weight per histogram,
  // the input arguments should be the same as in _get_fullMuonSF, and new ones should appear in the beginning
  
  // n_tnpBinNuisance is the number of TnP bins, used to set the size of the RVec that will be returned
  
  // tnpBinNuisance is supposed to start from 1 and be mapped into SF histogram bins as shown below
  //
  // pt | 7 | 8 | 9 |
  //     --- --- ---
  //    | 4 | 5 | 6 |
  //     --- --- ---
  //    | 1 | 2 | 3 |
  //              eta

  // if ptOther < 0 it is assumed only one lepton exists (so this function could also be used for wmass)
  // in that case the values are not used (and etaOther could actually take any value)
  
  ScaleFactorType sftype = charge > 0 ? isoTrigPlus : isoTrigMinus;
  if (not isoSF1) {
    sftype = charge > 0 ? antiisoTrigPlus : antiisoTrigMinus;
  }

  //std::cout << "Entry " << iEntry << ": era " << eraNames[era] << std::endl;
  //std::cout << "pt,eta       -> " << pt      << "," << eta      << std::endl;
  //std::cout << "pt,eta other -> " << ptOther << "," << etaOther << std::endl;

  auto const key = std::make_pair(sftype, era);
  const TH2D& hsf = scaleFactorHist.at(key);
  int nEtaBins = hsf.GetNbinsX();
  int nPtBins  = hsf.GetNbinsY();

  int ietaTnP = std::min(nEtaBins, std::max(1, hsf.GetXaxis()->FindFixBin(eta)));
  int iptTnP  = std::min(nPtBins,  std::max(1, hsf.GetYaxis()->FindFixBin(pt)));
  int tnpBinNuisance = ietaTnP + nEtaBins * (iptTnP - 1);
   // initialize to nominal SF
  Vec_f res(n_tnpBinNuisance, hsf.GetBinContent(ietaTnP, iptTnP));
  // sum error in specific bin
  res[tnpBinNuisance-1] += hsf.GetBinError(ietaTnP, iptTnP);
  
  if (ptOther > 0) {
    ScaleFactorType sftypeOther = isoNotrig;
    if (not isoSF2)
      sftypeOther = antiisoNotrig;
    auto const keyOther = std::make_pair(sftype, era);
    const TH2D& hsfOther = scaleFactorHist.at(keyOther);
    ietaTnP = std::min(nEtaBins, std::max(1, hsfOther.GetXaxis()->FindFixBin(etaOther)));
    iptTnP  = std::min(nPtBins,  std::max(1, hsfOther.GetYaxis()->FindFixBin(ptOther)));
    tnpBinNuisance = ietaTnP + nEtaBins * (iptTnP - 1);
    float tmp = res[tnpBinNuisance-1];
    float sf = hsf.GetBinContent(ietaTnP, iptTnP);
    res *= sf;
    res[tnpBinNuisance-1] = tmp * (sf + hsf.GetBinError(ietaTnP, iptTnP));
  }

  return res;
}


// might become obsolete (better to use next function directly)
Vec_f _get_fullSFvariation_wlike(int n_tnpBinNuisance,
				 float pt,      float eta, int charge,
				 float ptOther=-1, float etaOther=-1,
				 DataEra era = BToF,
				 bool isoSF1 = true,
				 bool isoSF2 = true
				 ) {
  
  // this is an helper function to define the effSystvariations for the Wlike analysis
  // idea is to fill again the alternative template for each effStat nuisance, where the
  // nuisance is defined for each single TnP eta-pt bin, and they will be considered as uncorrelated
  // so not as done in SMP-18-012 from the fit to efficiencies with Error function
  //
  // now, given that the SF weight is already defined for each process in the MCA file, in order to
  // define a histogram-dependent weight variation we have to add a weight which is (1 + epsilon),
  // where epsilon is the relative uncertainty on the SF (let's say on the product of SF for simplicity)
  // this is thus equivalent to using (SF + SF_err) as event weight for the alternate template

  Vec_f res(n_tnpBinNuisance, 1.0); // initialize to 1
  
  // n_tnpBinNuisance is the number of TnP bins, used to set the size of the RVec that will be returned
  
  // tnpBinNuisance is supposed to start from 1 and be mapped into SF histogram bins as shown below
  //
  // pt | 7 | 8 | 9 |
  //     --- --- ---
  //    | 4 | 5 | 6 |
  //     --- --- ---
  //    | 1 | 2 | 3 |
  //              eta

  // if ptOther < 0 it is assumed only one lepton exists (so this function could also be used for wmass)
  // in that case the values are not used (and etaOther could actually take any value)
  
  ScaleFactorType sftype = charge > 0 ? isoTrigPlus : isoTrigMinus;
  if (not isoSF1) {
    sftype = charge > 0 ? antiisoTrigPlus : antiisoTrigMinus;
  }

  //std::cout << "Entry " << iEntry << ": era " << eraNames[era] << std::endl;
  //std::cout << "pt,eta       -> " << pt      << "," << eta      << std::endl;
  //std::cout << "pt,eta other -> " << ptOther << "," << etaOther << std::endl;
  
  const TH2D& hsf = scaleFactorHist.at(std::make_pair(sftype, era));
  int nEtaBins = hsf.GetNbinsX();
  int nPtBins  = hsf.GetNbinsY();

  int ietaTnP = std::min(nEtaBins, std::max(1, hsf.GetXaxis()->FindFixBin(eta)));
  int iptTnP  = std::min(nPtBins,  std::max(1, hsf.GetYaxis()->FindFixBin(pt)));
  int tnpBinNuisance = ietaTnP + nEtaBins * (iptTnP - 1);
  res[tnpBinNuisance-1] = (1.0 + getRelUncertaintyFromTH2bin(hsf, ietaTnP, iptTnP));
  
  if (ptOther > 0) {
    ScaleFactorType sftypeOther = isoNotrig;
    if (not isoSF2)
      sftypeOther = antiisoNotrig;
    const TH2D& hsfOther = scaleFactorHist.at(std::make_pair(sftypeOther, era));
    ietaTnP = std::min(nEtaBins, std::max(1, hsfOther.GetXaxis()->FindFixBin(etaOther)));
    iptTnP  = std::min(nPtBins,  std::max(1, hsfOther.GetYaxis()->FindFixBin(ptOther)));
    tnpBinNuisance = ietaTnP + nEtaBins * (iptTnP - 1);
    res[tnpBinNuisance-1] *= (1.0 + getRelUncertaintyFromTH2bin(hsfOther, ietaTnP, iptTnP));
  }

  return res;
}

double qcdScaleWeight_VptBinned(const double& qcdscale, const double& vpt, const double& ptlow, const double& pthigh) {

  if (vpt >= ptlow and vpt < pthigh)
    return qcdscale;
  else
    return 1.0;
  
}

Vec_f qcdScaleWeight_VptBinned(const Vec_f& qcdscale, const double& vpt, const double& ptlow, const double& pthigh) {

  if (vpt >= ptlow and vpt < pthigh) {
    return qcdscale;
  } else {
    Vec_f res(qcdscale.size(),1.0); // initialize to 1   
    return res;
  }
  
}

int unroll2DTo1D_ptSlices(int pdgid, float pt, float eta){
  float ptmin = 0;
  float etaMax = 0;
  int nEtaBins = 0;
  if (abs(pdgid)==13) {
    ptmin = 26;
    etaMax = 2.4;
    nEtaBins = 48;
  } else {
    ptmin = 30;
    etaMax = 2.5;
    nEtaBins = 50;    
  }
  int etabin = (int) ((eta+etaMax)*10. );
  int ptbin  = (int) (pt-ptmin);
  return (ptbin*nEtaBins + etabin);
}

int unroll2DTo1D_etaSlices(int pdgid, float pt, float eta){
  float ptmin = 0;
  float etaMax = 0;
  int nptbins = 0;
  if (abs(pdgid)==13) {
    ptmin = 26;
    etaMax = 2.4;
    nptbins = 19;
  } else {
    ptmin = 30;
    etaMax = 2.5;
    nptbins = 15;    
  }
  int etabin = (int) ((eta+etaMax)*10. );
  int ptbin  = (int) (pt-ptmin );
  return (ptbin + nptbins * etabin);
}

// for muons, we have independent files for each charge
TFile *_file_effCov_trg_staterr_mu_plus = NULL;
TFile *_file_effCov_trg_staterr_mu_minus = NULL;
TH2F *_hist_relSystErr0_mu_plus = NULL;
TH2F *_hist_relSystErr1_mu_plus = NULL;
TH2F *_hist_relSystErr2_mu_plus = NULL;
TH2F *_hist_relSystErr0_mu_minus = NULL;
TH2F *_hist_relSystErr1_mu_minus = NULL;
TH2F *_hist_relSystErr2_mu_minus = NULL;
TFile *_file_effCov_trg_staterr_el = NULL;
TH2F *_hist_relSystErr0_el = NULL;
TH2F *_hist_relSystErr1_el = NULL;
TH2F *_hist_relSystErr2_el = NULL;

float effSystEtaBins(int inuisance, int pdgId, float eta, float pt, float etamin, float etamax, int charge=1) {
  // for charge, pass a positive or negative integer

  if (inuisance<0 || inuisance>2) {
    std::cout << "ERROR. Nuisance index " << inuisance << " not foreseen for the Erf with 3 parameters. Returning 0 syst." << std::endl;
    return 1.0;
  }

  float ret = 1.0;
  if (eta < etamin || eta > etamax) ret = 1.0;
  else {
    Double_t factor = 1.;
    TH2F *_hist_relSystErr;
    if(abs(pdgId)==11) {
      factor = 2.;
      if(!_file_effCov_trg_staterr_el) {
        _file_effCov_trg_staterr_el = new TFile("../postprocessing/data/leptonSF/new2016_madeSummer2018/systEff_trgel.root","read");
        _hist_relSystErr0_el = (TH2F*)_file_effCov_trg_staterr_el->Get("p0");
        _hist_relSystErr1_el = (TH2F*)_file_effCov_trg_staterr_el->Get("p1");
        _hist_relSystErr2_el = (TH2F*)_file_effCov_trg_staterr_el->Get("p2");
      }
      if      (inuisance==0) _hist_relSystErr = _hist_relSystErr0_el;
      else if (inuisance==1) _hist_relSystErr = _hist_relSystErr1_el;
      else                   _hist_relSystErr = _hist_relSystErr2_el;
    } else {
      factor = sqrt(2.);
      if(!_file_effCov_trg_staterr_mu_plus) {
        _file_effCov_trg_staterr_mu_plus = new TFile("../postprocessing/data/leptonSF/new2016_madeSummer2018/systEff_trgmu_plus_mu.root","read");
        _hist_relSystErr0_mu_plus = (TH2F*)_file_effCov_trg_staterr_mu_plus->Get("p0");
        _hist_relSystErr1_mu_plus = (TH2F*)_file_effCov_trg_staterr_mu_plus->Get("p1");
        _hist_relSystErr2_mu_plus = (TH2F*)_file_effCov_trg_staterr_mu_plus->Get("p2");
      }
      if(!_file_effCov_trg_staterr_mu_minus) {
        _file_effCov_trg_staterr_mu_minus = new TFile("../postprocessing/data/leptonSF/new2016_madeSummer2018/systEff_trgmu_minus_mu.root","read");
        _hist_relSystErr0_mu_minus = (TH2F*)_file_effCov_trg_staterr_mu_minus->Get("p0");
        _hist_relSystErr1_mu_minus = (TH2F*)_file_effCov_trg_staterr_mu_minus->Get("p1");
        _hist_relSystErr2_mu_minus = (TH2F*)_file_effCov_trg_staterr_mu_minus->Get("p2");
      }

      if      (inuisance==0) _hist_relSystErr = (charge > 0) ? _hist_relSystErr0_mu_plus : _hist_relSystErr0_mu_minus;
      else if (inuisance==1) _hist_relSystErr = (charge > 0) ? _hist_relSystErr1_mu_plus : _hist_relSystErr1_mu_minus;
      else                   _hist_relSystErr = (charge > 0) ? _hist_relSystErr2_mu_plus : _hist_relSystErr2_mu_minus;

    }
    
    int etabin = std::max(1, std::min(_hist_relSystErr->GetNbinsX(), _hist_relSystErr->GetXaxis()->FindFixBin(eta)));
    int ptbin  = std::max(1, std::min(_hist_relSystErr->GetNbinsY(), _hist_relSystErr->GetYaxis()->FindFixBin(pt)));
    // add protection for electrons close to gap (events in the gap are already excluded by selection)
    // if(abs(pdgId)==11 and fabs(eta) > 1.4 and fabs(eta) < 1.6 ) {
    //   ret = 1.0
    // } else {    
    ret = 1.0 + factor*_hist_relSystErr->GetBinContent(etabin,ptbin); //blow up the uncertainty
    //}
    // factor = sqrt(2) for muons. 
    // For electrons, need another factor sqrt(2) because the efficiency are made inclusively in charge, but treated as uncorrelated between W+ and W-
  }
  return ret;
}

// some additional older functions, not needed anymore with RDF, but to be checked
bool triggerMatchWlike_nano(int match1, int ch1, int match2, int ch2, ULong64_t evt) {

  if (isOddEvent(evt)) {
    if (ch1 > 0 and match1 > 0) 
      return true;
    else if (ch2 > 0 and match2 > 0) 
      return true;
    else 
      return false;
  } else {
    if (ch1 < 0 and match1 > 0) 
      return true;
    else if (ch2 < 0 and match2 > 0) 
      return true;
    else 
      return false;
  }

}

bool isInAccEtaPt(Float_t eta, Float_t pt, 
		  Float_t etalow, Float_t etahigh, 
		  Float_t ptlow, Float_t pthigh) {

  if (eta > etalow and eta < etahigh and pt > ptlow and pt < pthigh) 
    return 1;
  else 
    return 0;

}

bool isOutAccEtaPt(Float_t eta, Float_t pt, 
		   Float_t etalow, Float_t etahigh, 
		   Float_t ptlow, Float_t pthigh) {

  if (isInAccEtaPt(eta,pt,etalow,etahigh,ptlow,pthigh))
    return 0;
  else
    return 1;		   

}

// return value for gen lepton from Z that was not reconstructed when selecting events in W phase space
// so it is not a generic functions, assume you pass the reco muon pdgId (for first and only reco muon) to get the corresponding gen, assuming matching was done elsewhere and one is already selecting opposite charge leptons elsewhere 
// it also assumes the same pdgId was already selected for reco and gen, so only
float varGenLepFailReco(Int_t recoLepPdgId, Int_t genLep1_pdgId, Float_t genLep1_var, Float_t genLep2_var) {

  if (genLep1_pdgId == recoLepPdgId) 
    return genLep2_var;
  else
    return genLep1_var;

}


float safeRatio(float num, float den, float safe = 1.0) {
  return (den != 0.0) ? (num/den) : safe;
}

int regionIsoMt(bool lowIso, bool lowMt) {

  if      (not lowIso and     lowMt) return 0; // fakes region (failing isolation)
  else if (    lowIso and     lowMt) return 1; // fakes region (passing isolation)
  else if (not lowIso and not lowMt) return 2; // fakes application region
  else if (    lowIso and not lowMt) return 3; // signal region
  return -1;
  
}

#endif
