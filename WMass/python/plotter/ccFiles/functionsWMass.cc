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
#include <string>
#include <vector>
#include <unordered_map>
#include <utility>
#include <iostream>
#include <boost/algorithm/string/join.hpp>
#include <boost/functional/hash.hpp>

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

float getValFromTH2(const TH2& h, const float& x, const float& y) {
  //std::cout << "x,y --> " << x << "," << y << std::endl;
  int xbin = std::max(1, std::min(h.GetNbinsX(), h.GetXaxis()->FindFixBin(x)));
  int ybin  = std::max(1, std::min(h.GetNbinsY(), h.GetYaxis()->FindFixBin(y)));
  //std::cout << "xbin,ybin --> " << xbin << "," << ybin << std::endl;
  return h.GetBinContent(xbin,ybin);
}

bool isOddEvent(ULong64_t evt) {

  return (evt%2) ? 1 : 0;       

}

bool isEvenEvent(ULong64_t evt) {

  return (evt%2) ? 0 : 1;       

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

TRandom3 *rng_mt = NULL;
float mt_wlike_samesign_random(float pt1, float phi1, float pt2, float phi2, float met, float phimet, ULong64_t evt) {

  // for tests
  // randomly select one lepton to compute mt

  TVector2 metv = TVector2();
  metv.SetMagPhi(met,phimet);
  // use same vector and sum met
  TVector2 met_wlike = TVector2();

  Double_t ptL = 0.0;
  Double_t phiL = 0.0;
  if(!rng_mt) rng_mt = new TRandom3();
  // use eventNumber as seed
  rng_mt->SetSeed(evt); 
  if (rng_mt->Rndm() > 0.5) {
    ptL = pt1;
    phiL = phi1;
    met_wlike.SetMagPhi(pt2,phi2);
  } else {
    ptL = pt2;
    phiL = phi2;    
    met_wlike.SetMagPhi(pt1,phi1);
  }

  met_wlike += metv;
  return std::sqrt(2*ptL*met_wlike.Mod()*(1-std::cos(phiL-met_wlike.Phi())));

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

float mt_wlike_nano(float pt, float phi, float ptOther, float phiOther, float met, float phimet) {
  
  TVector2 pl = TVector2();
  pl.SetMagPhi(ptOther,phiOther);

  TVector2 metv = TVector2();
  metv.SetMagPhi(met,phimet);
  TVector2 met_wlike = pl+metv;

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


TFile* _file_prefireMapJets = NULL;
TH2F* prefireMapJets = NULL;

float mydeltaPhi(float phi1, float phi2) {
  float result = phi1 - phi2;
  while (result > float(M_PI)) result -= float(2*M_PI);
  while (result <= -float(M_PI)) result += float(2*M_PI);
  return result;
}

float mydeltaR(float eta1, float phi1, float eta2, float phi2) {
  float deta = std::abs(eta1-eta2);
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

//std::string _filename_allSF = "./testMuonSF/allSFs.root";
std::string _filename_allSF = "./testMuonSF/allSFs_eta0p1.root";

// Sorry you have to manually keep these consistent
typedef enum {BToH=0, BToF, GToH} DataEra;
typedef enum {MC=0, Data} DataType;
std::unordered_map<DataEra, std::string> eraNames = { {BToH, "BtoH"}, {BToF, "BtoF"}, {GToH, "GtoH"} };
std::unordered_map<DataType, std::string> datatypeNames = { {MC, "MC"}, {Data, "Data"} };

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

std::unordered_map<std::pair<std::string, DataEra>, TH2D, pair_hash> corrTypeToHist = {};
std::unordered_map<std::pair<std::string, DataType>, TH2D, pair_hash> prePostCorrToHist = {};
TFile _file_allSF = TFile(_filename_allSF.c_str(), "read");

void initializeScaleFactors() {
    if (!_file_allSF.IsOpen())
        std::cerr << "WARNING: Failed to open scaleFactors file " << _filename_allSF << "! No scale factors will be applied\n";

    std::cout << "INFO >>> Initializing histograms for SF from file " << _filename_allSF << std::endl;
    
    for (auto& era : eraNames) {
      for (auto& corr : {"trigger", "tracking", "idip", "iso", "isonotrig", "antiiso", "antiisonotrig"}) {
            std::vector<std::string> charges = {"both"};
            if (strcmp(corr, "trigger") == 0) {
                charges = {"plus", "minus"};
            }            
            for (auto& charge : charges) {
                std::vector<std::string> vars = {"SF2D", corr, era.second, charge};
                std::string corrname = boost::algorithm::join(vars, "_");
                auto* histptr = static_cast<TH2D*>(_file_allSF.Get(corrname.c_str()));
                if (histptr == nullptr)
                    std::cerr << "WARNING: Failed to load correction " << corrname << " in file "
                            << _filename_allSF << "! scale factors for this correction will be set to 1.0";
		histptr->SetDirectory(0);
                DataEra eraVal = era.first;
		// do not use "both" as charge key for the histograms, keep it simple
		std::string key = corr;
		if (charge != "both") {
		  key += charge;
		}
		// std::cout << "Histogram key " << key << " and era " << era.second << std::endl;
                auto corrKey = std::make_pair(key, eraVal);
                corrTypeToHist[corrKey] = *static_cast<TH2D*>(histptr);
            }
        }
    }

    for (auto& era : datatypeNames) {
      for (auto& corr : {"trigger", "tracking", "idip", "iso", "isonotrig"}) {
            std::vector<std::string> charges = {"both"};
            if (strcmp(corr, "trigger") == 0) {
                charges = {"plus", "minus"};
            }            
            for (auto& charge : charges) {
	      std::vector<std::string> vars = {"SF2D", era.second, "preOverPost", corr, charge};
                std::string corrname = boost::algorithm::join(vars, "_");
                auto* histptr = static_cast<TH2D*>(_file_allSF.Get(corrname.c_str()));
                if (histptr == nullptr)
                    std::cerr << "WARNING: Failed to load correction " << corrname << " in file "
                            << _filename_allSF << "! scale factors for this correction will be set to 1.0";
		histptr->SetDirectory(0);
                DataType typeVal = era.first;
		// do not use "both" as charge key for the histograms, keep it simple
		std::string key = corr;
		if (charge != "both") {
		  key += charge;
		}
		// std::cout << "Histogram key " << key << " and era " << era.second << std::endl;
                auto corrKey = std::make_pair(key, typeVal);
                prePostCorrToHist[corrKey] = *static_cast<TH2D*>(histptr);
            }
        }
    }

    _file_allSF.Close(); // should work since we used TH1D::SetDirectory(0) to detach histogram from file
    
}


// ongoing work
// std::unordered_map<std::string, DataEra>, TH2D, pair_hash> corrTypeToHist = {};;
// TFile _file_allSF_preOverPost = TFile(_filename_allSF_preOverPost.c_str(), "read");
/////////


float _get_AllMuonSF_fast_wlike(const float& pt,      const float& eta, const int& charge,
				const float& ptOther, const float& etaOther,
				DataEra era = BToH, bool noTrackingSF = false//, ULong64_t iEntry = 0
				) {
  if (corrTypeToHist.empty())
      return 1.;

  //std::cout << "Entry " << iEntry << ": era " << eraNames[era] << std::endl;
  //std::cout << "pt,eta       -> " << pt      << "," << eta      << std::endl;
  //std::cout << "pt,eta other -> " << ptOther << "," << etaOther << std::endl;
  float sf = 1.0;
  // not sure there is a more efficient way to compute the sf
  // some elements are common between the 2 leptons, some are not
  std::string triggerSF = charge > 0 ? "triggerplus" : "triggerminus";
  std::vector<std::string> sfnames = {triggerSF, "idip", "iso"};
  std::vector<std::string> sfnamesOther = {      "idip", "isonotrig"};
  if (not noTrackingSF) {
    sfnames.push_back("tracking");
    sfnamesOther.push_back("tracking");
  }
  for (const auto& corr : sfnames) {
    auto key = std::make_pair(corr, era);
    if (corrTypeToHist.find(key) != corrTypeToHist.end()) {
      sf *= getValFromTH2(corrTypeToHist[key],eta,pt);
      //std::cout << "scale factor main leg -> " << sf << std::endl;
    }
  }
  //std::cout << "scale factor main leg -> " << sf << std::endl;
  for (const auto& corr : sfnamesOther) {
    auto key = std::make_pair(corr, era);
    if (corrTypeToHist.find(key) != corrTypeToHist.end())
      sf *= getValFromTH2(corrTypeToHist[key],etaOther,ptOther);
  }
  //std::cout << "final scale factor -> " << sf << std::endl;
  return sf;
}

float _get_AllMuonSF_fast_wmass(const float& pt, const float& eta, const int& charge, DataEra era = BToH, bool noTrackingSF = false) {
  if (corrTypeToHist.empty())
      return 1.;
  
  float sf = 1.0;
  // not sure there is an efficient way to compute the sf
  // some elements are common between the 2 leptons, some are not
  std::string triggerSF = charge > 0 ? "triggerplus" : "triggerminus";
  std::vector<std::string> sfnames = {triggerSF, "idip", "iso"};
  if (not noTrackingSF) sfnames.push_back("tracking");
  for (const auto& corr : sfnames) {
    auto key = std::make_pair(corr, era);
    if (corrTypeToHist.find(key) != corrTypeToHist.end())
      sf *= getValFromTH2(corrTypeToHist[key],eta,pt);
  }
  return sf;
}

// generic function for a single leg, it is supposed to be called by other functions where the list of sf names was defined and passed to this one
float _get_singleMuonSF(const float& pt, const float& eta, const std::vector<std::string> &sfnames, DataEra era = BToH) {
  if (corrTypeToHist.empty())
      return 1.;
  
  float sf = 1.0;
  // not sure there is an efficient way to compute the sf
  // some elements are common between the 2 leptons, some are not
  for (const auto& corr : sfnames) {
    auto key = std::make_pair(corr, era);
    if (corrTypeToHist.find(key) != corrTypeToHist.end())
      sf *= getValFromTH2(corrTypeToHist[key],eta,pt);
  }
  return sf;
}

float _get_muonSF(const float& pt, const float& eta, const int& charge,
		  const bool trigger = true, const bool isolated = true,
		  DataEra era = BToH) {

  // trigger is used to decide whether the trigger sf has to be applied, and which isolation sf to use
  // it is supposed to be true for Wmass analysis, while for the Z Wlike analysis it should be true for the
  // specific lepton that is chosen to mimic the muon from a W, the other one gets no trigger sf (and isonotrig)

  std::vector<std::string> sfnames = {"tracking", "idip"}; // these should be always present
  std::string triggerSF = charge > 0 ? "triggerplus" : "triggerminus";
 
  if (isolated) {
    if (trigger) {
      sfnames.push_back("iso");
      sfnames.push_back(triggerSF);
    } else {
      sfnames.push_back("isonotrig");
    }
  } else {
    if (trigger) {
      sfnames.push_back("antiiso");
      sfnames.push_back(triggerSF);      
    } else {
      sfnames.push_back("antiisonotrig");
    }
  }
  return _get_singleMuonSF(pt, eta, sfnames, era);

}

float _get_AllMuonSF_fast_wlike_preOverPost(const float& pt,      const float& eta, const int& charge,
					    const float& ptOther, const float& etaOther,
					    DataType dtype = MC,
					    bool noTrackingSF = false,  bool noTriggerSF = false//, ULong64_t iEntry = 0
					    ) {
  if (prePostCorrToHist.empty())
      return 1.;

  //std::cout <<  "type " << datatypeNames[dtype] << std::endl;
  //std::cout << "pt,eta       -> " << pt      << "," << eta      << std::endl;
  //std::cout << "pt,eta other -> " << ptOther << "," << etaOther << std::endl;
  float sf = 1.0;
  // not sure there is a more efficient way to compute the sf
  // some elements are common between the 2 leptons, some are not
  std::string triggerSF = charge > 0 ? "triggerplus" : "triggerminus";
  std::vector<std::string> sfnames =      {"idip", "iso"};
  std::vector<std::string> sfnamesOther = {"idip", "isonotrig"};
  if (not noTriggerSF) {
    sfnames.push_back(triggerSF);
  }
  if (not noTrackingSF) {
    sfnames.push_back("tracking");
    sfnamesOther.push_back("tracking");
  }
  for (const auto& corr : sfnames) {
    auto key = std::make_pair(corr, dtype);
    if (prePostCorrToHist.find(key) != prePostCorrToHist.end()) {
      sf *= getValFromTH2(prePostCorrToHist[key],eta,pt);
      //std::cout << "scale factor main leg -> " << sf << std::endl;
    }
  }
  //std::cout << "scale factor main leg -> " << sf << std::endl;
  for (const auto& corr : sfnamesOther) {
    auto key = std::make_pair(corr, dtype);
    if (prePostCorrToHist.find(key) != prePostCorrToHist.end())
      sf *= getValFromTH2(prePostCorrToHist[key],etaOther,ptOther);
  }
  //std::cout << "final scale factor -> " << sf << std::endl;
  return sf;
}

float _get_AllMuonSF_fast_wlike_preOverPost_anyTrig(const float& pt,      const float& eta, const int& charge, const bool& trigMatch,
						    const float& ptOther, const float& etaOther, const bool& trigMatchOther,
						    DataType dtype = MC,
						    bool noTrackingSF = false,  bool noTriggerSF = false//, ULong64_t iEntry = 0
						    ) {
  if (prePostCorrToHist.empty())
      return 1.;

  //std::cout <<  "type " << datatypeNames[dtype] << std::endl;
  //std::cout << "pt,eta       -> " << pt      << "," << eta      << std::endl;
  //std::cout << "pt,eta other -> " << ptOther << "," << etaOther << std::endl;
  float sf = 1.0;
  // not sure there is a more efficient way to compute the sf
  // some elements are common between the 2 leptons, some are not
  std::vector<std::string> sfnames =      {"idip"};
  std::vector<std::string> sfnamesOther = {"idip"};
  if (not noTriggerSF) {
    if (trigMatch) {
      sfnames.push_back(charge > 0 ? "triggerplus" : "triggerminus");
      sfnames.push_back("iso");
      sfnamesOther.push_back("isonotrig");
    } else if (trigMatchOther) {
      sfnamesOther.push_back(charge > 0 ? "triggerminus" : "triggerplus");
      sfnames.push_back("isonotrig");
      sfnamesOther.push_back("iso");
    } else {
      // with the selection this case will never happen
      sfnames.push_back("triggerplus");
      sfnames.push_back("iso");
      sfnamesOther.push_back("isonotrig");
    }
  }
  if (not noTrackingSF) {
    sfnames.push_back("tracking");
    sfnamesOther.push_back("tracking");
  }
  for (const auto& corr : sfnames) {
    auto key = std::make_pair(corr, dtype);
    if (prePostCorrToHist.find(key) != prePostCorrToHist.end()) {
      sf *= getValFromTH2(prePostCorrToHist[key],eta,pt);
      //std::cout << "scale factor main leg -> " << sf << std::endl;
    }
  }
  //std::cout << "scale factor main leg -> " << sf << std::endl;
  for (const auto& corr : sfnamesOther) {
    auto key = std::make_pair(corr, dtype);
    if (prePostCorrToHist.find(key) != prePostCorrToHist.end())
      sf *= getValFromTH2(prePostCorrToHist[key],etaOther,ptOther);
  }
  //std::cout << "final scale factor -> " << sf << std::endl;
  return sf;
}


float _get_AllMuonSF_fast_wmass_preOverPost(const float& pt,      const float& eta, const int& charge,
					    DataType dtype = MC, bool noTrackingSF = false, bool noTriggerSF = false//, ULong64_t iEntry = 0
					    ) {
  if (prePostCorrToHist.empty())
      return 1.;

  //std::cout << "Entry " << iEntry << ": era " << eraNames[era] << std::endl;
  //std::cout << "pt,eta       -> " << pt      << "," << eta      << std::endl;
  //std::cout << "pt,eta other -> " << ptOther << "," << etaOther << std::endl;
  float sf = 1.0;
  // not sure there is a more efficient way to compute the sf
  // some elements are common between the 2 leptons, some are not
  std::string triggerSF = charge > 0 ? "triggerplus" : "triggerminus";
  std::vector<std::string> sfnames = {"idip", "iso"};
  if (not noTriggerSF) {
    sfnames.push_back(triggerSF);
  }
  if (not noTrackingSF) {
    sfnames.push_back("tracking");
  }
  for (const auto& corr : sfnames) {
    auto key = std::make_pair(corr, dtype);
    if (prePostCorrToHist.find(key) != prePostCorrToHist.end()) {
      sf *= getValFromTH2(prePostCorrToHist[key],eta,pt);
      //std::cout << "scale factor main leg -> " << sf << std::endl;
    }
  }
  //std::cout << "final scale factor -> " << sf << std::endl;
  return sf;
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

#endif
