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
// #include "EgammaAnalysis/ElectronTools/src/EnergyScaleCorrection_class.cc"

#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <cstdlib> //as stdlib.h      
#include <cstdio>
#include <string>
#include <vector>

TF1 * helicityFractionSimple_0 = new TF1("helicityFraction_0", "3./4*(TMath::Sqrt(1-x*x))^2", -1., 1.);
TF1 * helicityFractionSimple_L = new TF1("helicityFraction_L", "3./8.*(1-x)^2"              , -1., 1.);
TF1 * helicityFractionSimple_R = new TF1("helicityFraction_R", "3./8.*(1+x)^2"              , -1., 1.);

TFile *_file_helicityFractionsSimple = NULL;
TH2 * helicityFractionsSimple_0 = NULL;
TH2 * helicityFractionsSimple_L = NULL;
TH2 * helicityFractionsSimple_R = NULL;

static string _cmssw_base_ = string(getenv("CMSSW_BASE"));

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

float getValFromTH2(TH2*h = NULL, float x = 0.0, float y = 0.0) {

  if (!h) {
    cout << "Error in getValFromTH2(): uninitialized histogram. Abort" << endl;
    exit(EXIT_FAILURE);
  }
  int xbin = std::max(1, std::min(h->GetNbinsX(), h->GetXaxis()->FindFixBin(x)));
  int ybin  = std::max(1, std::min(h->GetNbinsY(), h->GetYaxis()->FindFixBin(y)));
  return h->GetBinContent(xbin,ybin);

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


// FIXME: better not to use this function, using the average is a bias
float mt_wlike_samesign(float pt1, float phi1, float pt2, float phi2, float met, float phimet) {

  // compute mt with both cases, and return average

  TVector2 metv = TVector2();
  metv.SetMagPhi(met,phimet);
  // use same vector and sum met
  TVector2 met_wlike = TVector2();

  met_wlike.SetMagPhi(pt2,phi2);
  met_wlike += metv;
  Double_t mt_wlike_1 = std::sqrt(2*pt1*met_wlike.Mod()*(1-std::cos(phi1-met_wlike.Phi())));

  met_wlike.SetMagPhi(pt1,phi1);
  met_wlike += metv;
  Double_t mt_wlike_2 = std::sqrt(2*pt2*met_wlike.Mod()*(1-std::cos(phi2-met_wlike.Phi())));

  return 0.5 * (mt_wlike_1 + mt_wlike_2);

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
 

  if (_cmssw_base_ == "") {
    cout << "Setting _cmssw_base_ to environment variable CMSSW_BASE" << endl;
    _cmssw_base_ = getEnvironmentVariable("CMSSW_BASE");
  }

  TH2D* hratio_FSR_noFSR = NULL;
  Double_t outlierWgt = 0.0;  // hardcoded below, taken from ratio of yields outside acceptance

  if (abs(pdgId)==11) {

    if (!ratio_FSRoverNoFSR_etaPt_el) {
      _file_ratio_FSRoverNoFSR_etaPt_el = new TFile(Form("%s/src/CMGTools/WMass/python/plotter/w-mass-13TeV/theoryReweighting/FSR_atGenLevel_electron_genEtaPtAnalysis.root",_cmssw_base_.c_str()),"read");
      ratio_FSRoverNoFSR_etaPt_el = (TH2D*) _file_ratio_FSRoverNoFSR_etaPt_el->Get("ratio__Wzpt_el__Wfsr_el");
    }
    hratio_FSR_noFSR = ratio_FSRoverNoFSR_etaPt_el;
    outlierWgt = 0.964796;

  } else if (abs(pdgId)==13) {

     if (!ratio_FSRoverNoFSR_etaPt_mu) {
      _file_ratio_FSRoverNoFSR_etaPt_mu = new TFile(Form("%s/src/CMGTools/WMass/python/plotter/w-mass-13TeV/theoryReweighting/FSR_atGenLevel_muon_genEtaPtAnalysis.root",_cmssw_base_.c_str()),"read");
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

  if (_cmssw_base_ == "") {
    cout << "Setting _cmssw_base_ to environment variable CMSSW_BASE" << endl;
    _cmssw_base_ = getEnvironmentVariable("CMSSW_BASE");
  }
  
  if (!fsrWeights_el || !fsrWeights_mu) {
    _file_fsrWeights_simple = new TFile(Form("%s/src/CMGTools/WMass/python/plotter/w-mass-13TeV/theoryReweighting/photos_rwgt_integrated.root",_cmssw_base_.c_str()),"read");
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
  if (_cmssw_base_ == "") {
    cout << "Setting _cmssw_base_ to environment variable CMSSW_BASE" << endl;
    _cmssw_base_ = getEnvironmentVariable("CMSSW_BASE");
  }
  if (!fsrWeights_elplus || !fsrWeights_elminus || !fsrWeights_muplus || !fsrWeights_muminus) {
    _file_fsrWeights = new TFile(Form("%s/src/CMGTools/WMass/python/plotter/w-mass-13TeV/theoryReweighting/photos_rwgt.root",_cmssw_base_.c_str()),"read");
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
  if (_cmssw_base_ == "") {
    cout << "Setting _cmssw_base_ to environment variable CMSSW_BASE" << endl;
    _cmssw_base_ = getEnvironmentVariable("CMSSW_BASE");
  }
  if (!amcnlody) {
    _file_dyptWeights = new TFile(Form("%s/src/CMGTools/WMass/python/plotter/w-mass-13TeV/theoryReweighting/zpt_weights.root",_cmssw_base_.c_str()));
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
  if (_cmssw_base_ == "") {
    cout << "Setting _cmssw_base_ to environment variable CMSSW_BASE" << endl;
    _cmssw_base_ = getEnvironmentVariable("CMSSW_BASE");
  }
  if (!amcnlody) {
    _file_dyptWeights = new TFile(Form("%s/src/CMGTools/WMass/python/plotter/w-mass-13TeV/theoryReweighting/zpt_weights.root",_cmssw_base_.c_str()));
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
    _file_dyptWeights_MGaMCAtNLO_PowhegMiNNLO = new TFile(Form("%s/src/CMGTools/WMass/python/plotter/w-mass-13TeV/theoryReweighting/ZptRatio_MGaMCAtNLO_PowhegMiNNLO_dressed_noCuts.root",_cmssw_base_.c_str()));
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
  if (_cmssw_base_ == "") {
    cout << "Setting _cmssw_base_ to environment variable CMSSW_BASE" << endl;
    _cmssw_base_ = getEnvironmentVariable("CMSSW_BASE");
  }
  TH1F* hist_scales = NULL;
  if (!_file_postfitWeights) {
    if (flav==0)
      _file_postfitWeights = new TFile(Form("%s/src/CMGTools/WMass/python/plotter/w-helicity-13TeV/theoryReweighting/postfit_wgts_fitlep.root",_cmssw_base_.c_str()));
    else if (flav==13)
      _file_postfitWeights = new TFile(Form("%s/src/CMGTools/WMass/python/plotter/w-helicity-13TeV/theoryReweighting/postfit_wgts_fitmu.root",_cmssw_base_.c_str()));
    else if (flav==11)
      _file_postfitWeights = new TFile(Form("%s/src/CMGTools/WMass/python/plotter/w-helicity-13TeV/theoryReweighting/postfit_wgts_fitel.root",_cmssw_base_.c_str()));
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



//// ALL THE FOLLOWING PART UNTIL "//// XXXX" IS OBSOLETE, CAN BE REMOVED AT SOME POINT

// ==================
// full muon pT with scale variations
TFile *_file_residualcorr_scaleMu = NULL;
TH2D *_histo_residualcorr_scaleMu = NULL;

float residualScaleMu(float pt, float eta, int isData, const char *fileCorr="../../data/muonscale/scale_correction_nonclosure_mu.root") {
  if(!isData) return 1.;

  if(!_histo_residualcorr_scaleMu) {
    _file_residualcorr_scaleMu = new TFile(fileCorr);
    _histo_residualcorr_scaleMu = (TH2D*)(_file_residualcorr_scaleMu->Get("scales_corrections_2d_muons"));
  }
  
  TH2D *hist = _histo_residualcorr_scaleMu;
  int etabin  = std::max(1, std::min(hist->GetNbinsX(), hist->GetXaxis()->FindFixBin(fabs(eta))));
  int ptbin = std::max(1, std::min(hist->GetNbinsY(), hist->GetYaxis()->FindFixBin(pt)));
  
  float scale = 1. - hist->GetBinContent(etabin,ptbin);
  if (scale < 0) {
    cout << "WARNING in residualScale() function: scale < 0 --> returning 0." << endl;
    return 0;
  } else {
    return scale;
  }

}

float ptScaleUncorr(float pt, float eta, int pdgid, int iPtVar, bool isUp = true) {

  // iPtVar defines the number fo the nuisance. For example, for muons we have 2 of them: 
  // one is flat over all eta and has a 0.3% pt variation, while the other has a non-zero pt variation only for |eta|>2.1 (we assign 0.95%)

  float abseta = fabs(eta);
  float ptsyst = 0.0;

  if (fabs(pdgid) == 11) {

    // electrons, numbers to be discussed further
    if      (iPtVar == 0)  ptsyst = 0.003; //0.3%
    else if (iPtVar == 1) {
        if (abseta >= 1.0) ptsyst = 0.004; // 0.5%
    } else if (iPtVar == 2) {
        if (abseta >= 1.5) ptsyst = 0.0063; // 0.8%
    } else if (iPtVar == 3) {
        if (abseta >= 2.1) ptsyst = 0.008; // 1.0%
    }    

  } else if (fabs(pdgid) == 13) {

    //muons
    if      (iPtVar == 0)  ptsyst = 0.003; //0.3%
    else if (iPtVar == 1) {
      if (abseta >= 2.1) ptsyst = 0.01; // 1.0%
    }

  } else {
    // this should not happen using pdgid
    // cout << "Warning in ptScaleUncorr(): pdg ID not consistent with mu or ele. Returning 0" << endl;
    return 0;    
  }

  if (not isUp) ptsyst *= -1.0; // switch sign, will decrease pt
  return pt * (1 + ptsyst);
  
}

//-------------------
// muons

float ptMuScaleUncorr0Up(float pt, float eta) {
  return ptScaleUncorr(pt, eta, 13, 0, true);
}

float ptMuScaleUncorr0Dn(float pt, float eta) {
  return ptScaleUncorr(pt, eta, 13, 0, false);
}

float ptMuScaleUncorr1Up(float pt, float eta) {
  return ptScaleUncorr(pt, eta, 13, 1, true);
}

float ptMuScaleUncorr1Dn(float pt, float eta) {
  return ptScaleUncorr(pt, eta, 13, 1, false);
}

//// XXXX (PART ABOVE IS OBSOLETE)

//===============================================


//==================================================
TRandom3 *rng = NULL;

float getSmearedVar(float var, float smear, ULong64_t eventNumber, int isData, bool smearOnlyMC=false) {

  if (smearOnlyMC && isData) return var;

  if(!rng) rng = new TRandom3();
  // use eventNumber as seed, otherwise each time the function is called for the same event, the smearer produce a different smeared value
  rng->SetSeed(eventNumber); // make it really random across different jobs    
  return var * ( 1.0 + smear * rng->Gaus());

}

// old function, can probably be removed
float triggerSF_2l(float l11pass, float l12pass, float l21pass, float l22pass, float l1sf, float l2sf, int randomize=0){
  float weight = -999.;

  bool l1pass = (l11pass > -1. || l12pass > -1.);
  bool l2pass = (l21pass > -1. || l22pass > -1.);

  if      ( l1pass && !l2pass) weight = l1sf;
  else if (!l1pass &&  l2pass) weight = l2sf;
  else if ( l1pass &&  l2pass) {
    if   (!randomize) weight = (l1sf+l2sf)/2.;
    else {
      if(!rng) rng = new TRandom3();
      float randy = rng->Uniform(-1.,1.);
      if (randy < 0.) weight = l1sf;
      else            weight = l2sf;
    }
  }

  //else return -999.; // this should never happen

  return weight;
}

//============================     

TFile *_file_recoToSelection_leptonSF_mu = NULL;
TH2F *_histo_recoToSelection_leptonSF_mu = NULL;

float _get_muonSF_recoToSelection(int pdgid, float pt, float eta, bool useBinnedSF = false) {

  if (_cmssw_base_ == "") {
    cout << "Setting _cmssw_base_ to environment variable CMSSW_BASE" << endl;
    _cmssw_base_ = getEnvironmentVariable("CMSSW_BASE");
  }

  string hSFname = "scaleFactor_etaInterpolated";
  if (useBinnedSF) hSFname = "scaleFactorOriginal";

  if (!_histo_recoToSelection_leptonSF_mu) {
    _file_recoToSelection_leptonSF_mu = new TFile(Form("%s/src/CMGTools/WMass/python/postprocessing/data/leptonSF/new2016_madeSummer2018/smoothEfficiency_muons_recoToSel_finerETA.root",_cmssw_base_.c_str()),"read");
    _histo_recoToSelection_leptonSF_mu = (TH2F*)(_file_recoToSelection_leptonSF_mu->Get(hSFname.c_str()));
  }

  if(abs(pdgid)==13) {
    TH2F *histRecoToSelection = _histo_recoToSelection_leptonSF_mu;

    int etabin = std::max(1, std::min(histRecoToSelection->GetNbinsX(), histRecoToSelection->GetXaxis()->FindFixBin(eta)));
    int ptbin  = std::max(1, std::min(histRecoToSelection->GetNbinsY(), histRecoToSelection->GetYaxis()->FindFixBin(pt)));

    return histRecoToSelection->GetBinContent(etabin,ptbin);
  }

  return 0;

}


//============================


TFile *_file_trigger_leptonSF_mu_plus = NULL;
TH2F *_histo_trigger_leptonSF_mu_plus = NULL;
TFile *_file_trigger_leptonSF_mu_minus = NULL;
TH2F *_histo_trigger_leptonSF_mu_minus = NULL;
// following are used to test effstat uncertainties using the stat uncertainty from binned efficiencies only
TFile *_file_triggerBinUncEffStat_leptonSF_mu_plus = NULL;
TH2F *_histo_triggerBinUncEffStat_leptonSF_mu_plus = NULL;
TFile *_file_triggerBinUncEffStat_leptonSF_mu_minus = NULL;
TH2F *_histo_triggerBinUncEffStat_leptonSF_mu_minus = NULL;

float _get_muonSF_selectionToTrigger(int pdgid, float pt, float eta, int charge, 
				     float sumErrorTimesThis = 0.0, bool useStatErrOnly = false,
				     bool useBinnedSF = false) {

  if (_cmssw_base_ == "") {
    cout << "Setting _cmssw_base_ to environment variable CMSSW_BASE" << endl;
    _cmssw_base_ = getEnvironmentVariable("CMSSW_BASE");
  }

  string hSFname = "scaleFactor";
  if (useBinnedSF) hSFname = "scaleFactorOriginal";

  if (!_histo_trigger_leptonSF_mu_plus) {
    _file_trigger_leptonSF_mu_plus = new TFile(Form("%s/src/CMGTools/WMass/python/postprocessing/data/leptonSF/new2016_madeSummer2018/smoothEfficiency_muons_plus_trigger.root",_cmssw_base_.c_str()),"read");
    _histo_trigger_leptonSF_mu_plus = (TH2F*)(_file_trigger_leptonSF_mu_plus->Get(hSFname.c_str()));
  }
  if (!_histo_trigger_leptonSF_mu_minus) {
    _file_trigger_leptonSF_mu_minus = new TFile(Form("%s/src/CMGTools/WMass/python/postprocessing/data/leptonSF/new2016_madeSummer2018/smoothEfficiency_muons_minus_trigger.root",_cmssw_base_.c_str()),"read");
    _histo_trigger_leptonSF_mu_minus = (TH2F*)(_file_trigger_leptonSF_mu_minus->Get(hSFname.c_str()));
  }


  if (useStatErrOnly) {
    // since July 2019 we could take them from the same file as above, as we stored these numbers there, but whatever works is fine)
    if (!_histo_triggerBinUncEffStat_leptonSF_mu_plus) {
      _file_triggerBinUncEffStat_leptonSF_mu_plus = new TFile(Form("%s/src/CMGTools/WMass/python/postprocessing/data/leptonSF/new2016_madeSummer2018/TnPstuff/muon/triggerMuonEffPlus_fromRooFitResult_onlyStatUnc.root",_cmssw_base_.c_str()),"read");
      _histo_triggerBinUncEffStat_leptonSF_mu_plus = (TH2F*)(_file_triggerBinUncEffStat_leptonSF_mu_plus->Get("triggerSF_plus"));
    }
    if (!_histo_triggerBinUncEffStat_leptonSF_mu_minus) {
      _file_triggerBinUncEffStat_leptonSF_mu_minus = new TFile(Form("%s/src/CMGTools/WMass/python/postprocessing/data/leptonSF/new2016_madeSummer2018/TnPstuff/muon/triggerMuonEffMinus_fromRooFitResult_onlyStatUnc.root",_cmssw_base_.c_str()),"read");
      _histo_triggerBinUncEffStat_leptonSF_mu_minus = (TH2F*)(_file_triggerBinUncEffStat_leptonSF_mu_minus->Get("triggerSF_minus"));
    }
  }


  if(abs(pdgid)==13) {
    TH2F *histTrigger = ( charge > 0 ? _histo_trigger_leptonSF_mu_plus : _histo_trigger_leptonSF_mu_minus );

    int etabin = std::max(1, std::min(histTrigger->GetNbinsX(), histTrigger->GetXaxis()->FindFixBin(eta)));
    int ptbin  = std::max(1, std::min(histTrigger->GetNbinsY(), histTrigger->GetYaxis()->FindFixBin(pt)));

    //float out = histTrigger->GetBinContent(etabin,ptbin) + sumErrorTimesThis * histTrigger->GetBinError(etabin,ptbin);
    TH2F *histTriggerBinUnc = nullptr; 
    double binunc = histTrigger->GetBinError(etabin,ptbin);
    if (useStatErrOnly) {
      // see above, this would no longer be needed
      histTriggerBinUnc = ( charge > 0 ? _histo_triggerBinUncEffStat_leptonSF_mu_plus : _histo_triggerBinUncEffStat_leptonSF_mu_minus );
      int etabinUnc = std::max(1, std::min(histTriggerBinUnc->GetNbinsX(), histTriggerBinUnc->GetXaxis()->FindFixBin(eta)));
      int ptbinUnc  = std::max(1, std::min(histTriggerBinUnc->GetNbinsY(), histTriggerBinUnc->GetYaxis()->FindFixBin(pt)));      
      binunc = histTriggerBinUnc->GetBinError(etabinUnc,ptbinUnc);
    }
    // inflate error by sqrt(2.) to account for other SF (here only trigger is considered)
    float out = histTrigger->GetBinContent(etabin,ptbin) + sumErrorTimesThis * TMath::Sqrt(2.) * binunc;

    return out;
  }

  return 0;

}

// return both SF, possibily only for events fulfilling a given condition
float _get_muonSF_TriggerAndIDiso(int pdgid, float pt, float eta, int charge, bool passCut = true) {

  if (passCut) {
    return _get_muonSF_recoToSelection(pdgid,pt,eta)*_get_muonSF_selectionToTrigger(pdgid,pt,eta,charge);
  } else {
    return 1.0;
  }

}

//============================

TFile *_file_trigger_leptonSF_fast_mu_plus = NULL;
TH2F *_histo_trigger_leptonSF_fast_mu_plus = NULL;
TFile *_file_trigger_leptonSF_fast_mu_minus = NULL;
TH2F *_histo_trigger_leptonSF_fast_mu_minus = NULL;
TFile *_file_selection_leptonSF_fast_mu_plus = NULL;
TH2F *_histo_selection_leptonSF_fast_mu_plus = NULL;
TFile *_file_selection_leptonSF_fast_mu_minus = NULL;
TH2F *_histo_selection_leptonSF_fast_mu_minus = NULL;

float _get_muonSF_fast_wlike(float pt1, float eta1, int charge1, float pt2, float eta2, int charge2, ULong64_t event) {

  // to be improved, still experimental
  float pt = 0.0;
  float eta = 0.0;
  int charge = 0.0;

  // positive (negative) leptons on odd (even) events
  if (isOddEvent(event)) {
    if (charge1 > 0) {
      pt = pt1;
      eta = eta1;
      charge = charge1;
    } else {
      pt = pt2;
      eta = eta2;
      charge = charge2;
    }
  } else {
    if (charge1 < 0) {
      pt = pt1;
      eta = eta1;
      charge = charge1;
    } else {
      pt = pt2;
      eta = eta2;
      charge = charge2;
    }
  }

  
  if (_cmssw_base_ == "") {
    cout << "Setting _cmssw_base_ to environment variable CMSSW_BASE" << endl;
    _cmssw_base_ = getEnvironmentVariable("CMSSW_BASE");
  }

  string hSFname = "EGamma_SF2D";

  // trigger
  if (!_histo_trigger_leptonSF_fast_mu_plus) {
    _file_trigger_leptonSF_fast_mu_plus = new TFile(Form("%s/src/CMGTools/WMass/python/plotter/testMuonSF/triggerMuPlus.root",_cmssw_base_.c_str()),"read");
    _histo_trigger_leptonSF_fast_mu_plus = (TH2F*)(_file_trigger_leptonSF_fast_mu_plus->Get(hSFname.c_str()));
  }
  if (!_histo_trigger_leptonSF_fast_mu_minus) {
    _file_trigger_leptonSF_fast_mu_minus = new TFile(Form("%s/src/CMGTools/WMass/python/plotter/testMuonSF/triggerMuMinus.root",_cmssw_base_.c_str()),"read");
    _histo_trigger_leptonSF_fast_mu_minus = (TH2F*)(_file_trigger_leptonSF_fast_mu_minus->Get(hSFname.c_str()));
  }
  // selection
  if (!_histo_selection_leptonSF_fast_mu_plus) {
    _file_selection_leptonSF_fast_mu_plus = new TFile(Form("%s/src/CMGTools/WMass/python/plotter/testMuonSF/selectionMuPlus.root",_cmssw_base_.c_str()),"read");
    _histo_selection_leptonSF_fast_mu_plus = (TH2F*)(_file_selection_leptonSF_fast_mu_plus->Get(hSFname.c_str()));
  }
  if (!_histo_selection_leptonSF_fast_mu_minus) {
    _file_selection_leptonSF_fast_mu_minus = new TFile(Form("%s/src/CMGTools/WMass/python/plotter/testMuonSF/selectionMuMinus.root",_cmssw_base_.c_str()),"read");
    _histo_selection_leptonSF_fast_mu_minus = (TH2F*)(_file_selection_leptonSF_fast_mu_minus->Get(hSFname.c_str()));
  }

  TH2F *histTrigger = ( charge > 0 ? _histo_trigger_leptonSF_fast_mu_plus : _histo_trigger_leptonSF_fast_mu_minus );
  int etabin = std::max(1, std::min(histTrigger->GetNbinsX(), histTrigger->GetXaxis()->FindFixBin(eta)));
  int ptbin  = std::max(1, std::min(histTrigger->GetNbinsY(), histTrigger->GetYaxis()->FindFixBin(pt)));

  float out = histTrigger->GetBinContent(etabin,ptbin);

  TH2F *histSelectionPlus  = _histo_selection_leptonSF_fast_mu_plus;
  TH2F *histSelectionMinus = _histo_selection_leptonSF_fast_mu_minus;
  float ptPlus = 0;
  float ptMinus = 0;
  float etaPlus = 0;
  float etaMinus = 0;
  if (charge1 > 0) {
    ptPlus = pt1;
    etaPlus = eta1;
    ptMinus = pt2;
    etaMinus = eta2;
  } else {
    ptPlus = pt2;
    etaPlus = eta2;
    ptMinus = pt1;
    etaMinus = eta1;
  }

  etabin = std::max(1, std::min(histSelectionPlus->GetNbinsX(), histSelectionPlus->GetXaxis()->FindFixBin(etaPlus)));
  ptbin  = std::max(1, std::min(histSelectionPlus->GetNbinsY(), histSelectionPlus->GetYaxis()->FindFixBin(ptPlus)));
  out *= histSelectionPlus->GetBinContent(etabin,ptbin);
  etabin = std::max(1, std::min(histSelectionMinus->GetNbinsX(), histSelectionMinus->GetXaxis()->FindFixBin(etaMinus)));
  ptbin  = std::max(1, std::min(histSelectionMinus->GetNbinsY(), histSelectionMinus->GetYaxis()->FindFixBin(ptMinus)));
  out *= histSelectionMinus->GetBinContent(etabin,ptbin);

  return out;

}


TFile *_file_allSF = NULL;
// 3 elements, depending on what dataset period one is using for the scale factors
TH2F *_histo_trigger_minus[3] = {NULL, NULL, NULL}; 
TH2F *_histo_trigger_plus[3] = {NULL, NULL, NULL};
TH2F *_histo_iso[3] = {NULL, NULL, NULL};
TH2F *_histo_isonotrig[3] = {NULL, NULL, NULL};
TH2F *_histo_idip[3] = {NULL, NULL, NULL};
TH2F *_histo_tracking[3] = {NULL, NULL, NULL};

// float* getParticleVars(float pt1, float eta1, float pt2, float eta2) {
  
//   static float vars[] = {pt1, eta1, pt2, eta2};
//   std::cout << "1) " << pt1 << "   " << eta1 << "   " << pt2 << "   " << eta2 << std::endl;
//   return vars;

// }


// I think the function cannot have more than 8 arguments, let's see
float _get_AllMuonSF_fast_wlike(float pt1, float eta1, int charge1, int trigMatch1, float pt2, float eta2, int charge2, int trigMatch2, ULong64_t event, int era = 0) {

  // era = 0 to pick SF for all year, 1 for BtoF and 2 for GtoH
  string dataEraForSF = "BtoH";
  if (era == 1) 
    dataEraForSF = "BtoF";
  else if (era == 2 )
    dataEraForSF = "GtoH";

  // to be improved, still experimental
  float pt = 0.0;
  float eta = 0.0;
  int charge = 0.0;
  // positive (negative) leptons on odd (even) events
  if (isOddEvent(event)) {
    if (charge1 > 0) {
      pt = pt1;
      eta = eta1;
      charge = charge1;
    } else {
      pt = pt2;
      eta = eta2;
      charge = charge2;
    }
  } else {
    if (charge1 < 0) {
      pt = pt1;
      eta = eta1;
      charge = charge1;
    } else {
      pt = pt2;
      eta = eta2;
      charge = charge2;
    }
  }

  
  if (_cmssw_base_ == "") {
    cout << "Setting _cmssw_base_ to environment variable CMSSW_BASE" << endl;
    _cmssw_base_ = getEnvironmentVariable("CMSSW_BASE");
  }

  // open the single file we have
  if (!_file_allSF) {
    _file_allSF = new TFile(Form("%s/src/CMGTools/WMass/python/plotter/testMuonSF/allSFs.root",_cmssw_base_.c_str()),"read");
  }

  // trigger
  if (!_histo_trigger_plus[era]) {
    _histo_trigger_plus[era] = (TH2F*)(_file_allSF->Get(Form("SF2D_trigger_%s_plus",dataEraForSF.c_str())));
  }
  if (!_histo_trigger_minus[era]) {
    _histo_trigger_minus[era] = (TH2F*)(_file_allSF->Get(Form("SF2D_trigger_%s_minus",dataEraForSF.c_str())));
  }
  // ID + ip (interaction point, i.e. dz and dxy)
  if (!_histo_idip[era]) {
    _histo_idip[era] = (TH2F*)(_file_allSF->Get(Form("SF2D_idip_%s_both",dataEraForSF.c_str())));
  }
  if (!_histo_tracking[era]) {
    _histo_tracking[era] = (TH2F*)(_file_allSF->Get(Form("SF2D_tracking_%s_both",dataEraForSF.c_str())));
  }
  if (!_histo_iso[era]) {
    _histo_iso[era] = (TH2F*)(_file_allSF->Get(Form("SF2D_iso_%s_both",dataEraForSF.c_str())));
  }
  if (!_histo_isonotrig[era]) {
    _histo_isonotrig[era] = (TH2F*)(_file_allSF->Get(Form("SF2D_isonotrig_%s_both",dataEraForSF.c_str())));
  }


  TH2F *histTrigger = ( charge > 0 ? _histo_trigger_plus[era] : _histo_trigger_minus[era] );
  float sf = getValFromTH2(histTrigger,eta,pt);
  sf *= getValFromTH2(_histo_idip[era],eta1,pt1) * getValFromTH2(_histo_idip[era],eta2,pt2);
  sf *= getValFromTH2(_histo_tracking[era],eta1,pt1) * getValFromTH2(_histo_tracking[era],eta2,pt2);
  // sf *= getValFromTH2(_histo_iso[era],eta1,pt1) * getValFromTH2(_histo_iso[era],eta2,pt2);
  // for SF validation the SF for iso should be applied on the second lepton only if it had passed the trigger
  // so we should not apply this SF, but then also not applying thr isolation cut on that leg, but only on the
  // selected one
  // sf *= getValFromTH2(_histo_iso[era],eta,pt);
  TH2F *histIso1 = (trigMatch1 ? _histo_iso[era] : _histo_isonotrig[era]);
  TH2F *histIso2 = (trigMatch2 ? _histo_iso[era] : _histo_isonotrig[era]);
  sf *= getValFromTH2(histIso1,eta1,pt1) * getValFromTH2(histIso2,eta2,pt2);
  // cout << " sf = " << sf << endl;

  return sf;

}


// details like histograms and location might be updated
float _get_AllMuonSF_fast_wmass(float pt, float eta, int charge, int trigMatch, int era=0) {

  // era = 0 to pick SF for all year, 1 for BtoF and 2 for GtoH
  string dataEraForSF = "BtoH";
  if (era == 1) 
    dataEraForSF = "BtoF";
  else if (era == 2 )
    dataEraForSF = "GtoH";
  
  if (_cmssw_base_ == "") {
    cout << "Setting _cmssw_base_ to environment variable CMSSW_BASE" << endl;
    _cmssw_base_ = getEnvironmentVariable("CMSSW_BASE");
  }

  // open the single file we have
  if (!_file_allSF) {
    _file_allSF = new TFile(Form("%s/src/CMGTools/WMass/python/plotter/testMuonSF/allSFs.root",_cmssw_base_.c_str()),"read");
  }

  // trigger
  if (!_histo_trigger_plus[era]) {
    _histo_trigger_plus[era] = (TH2F*)(_file_allSF->Get(Form("SF2D_trigger_%s_plus",dataEraForSF.c_str())));
  }
  if (!_histo_trigger_minus[era]) {
    _histo_trigger_minus[era] = (TH2F*)(_file_allSF->Get(Form("SF2D_trigger_%s_minus",dataEraForSF.c_str())));
  }
  // ID + ip (interaction point, i.e. dz and dxy)
  if (!_histo_idip[era]) {
    _histo_idip[era] = (TH2F*)(_file_allSF->Get(Form("SF2D_idip_%s_both",dataEraForSF.c_str())));
  }
  if (!_histo_tracking[era]) {
    _histo_tracking[era] = (TH2F*)(_file_allSF->Get(Form("SF2D_tracking_%s_both",dataEraForSF.c_str())));
  }
  if (!_histo_iso[era]) {
    _histo_iso[era] = (TH2F*)(_file_allSF->Get(Form("SF2D_iso_%s_both",dataEraForSF.c_str())));
  }
  if (!_histo_isonotrig[era]) {
    _histo_isonotrig[era] = (TH2F*)(_file_allSF->Get(Form("SF2D_isonotrig_%s_both",dataEraForSF.c_str())));
  }


  TH2F *histTrigger = ( charge > 0 ? _histo_trigger_plus[era] : _histo_trigger_minus[era] );
  float sf = getValFromTH2(histTrigger,eta,pt);
  sf *= getValFromTH2(_histo_idip[era],eta,pt);
  sf *= getValFromTH2(_histo_tracking[era],eta,pt);
  TH2F *histIso = (trigMatch ? _histo_iso[era] : _histo_isonotrig[era]);
  sf *= getValFromTH2(histIso,eta,pt);
  // cout << " sf = " << sf << endl;

  return sf;

}

//============================

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
    if (_cmssw_base_ == "") {
      cout << "Setting _cmssw_base_ to environment variable CMSSW_BASE" << endl;
      _cmssw_base_ = getEnvironmentVariable("CMSSW_BASE");
    }
    TH2F *_hist_relSystErr;
    if(abs(pdgId)==11) {
      factor = 2.;
      if(!_file_effCov_trg_staterr_el) {
        _file_effCov_trg_staterr_el = new TFile(Form("%s/src/CMGTools/WMass/python/postprocessing/data/leptonSF/new2016_madeSummer2018/systEff_trgel.root",_cmssw_base_.c_str()),"read");
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
        _file_effCov_trg_staterr_mu_plus = new TFile(Form("%s/src/CMGTools/WMass/python/postprocessing/data/leptonSF/new2016_madeSummer2018/systEff_trgmu_plus_mu.root",_cmssw_base_.c_str()),"read");
        _hist_relSystErr0_mu_plus = (TH2F*)_file_effCov_trg_staterr_mu_plus->Get("p0");
        _hist_relSystErr1_mu_plus = (TH2F*)_file_effCov_trg_staterr_mu_plus->Get("p1");
        _hist_relSystErr2_mu_plus = (TH2F*)_file_effCov_trg_staterr_mu_plus->Get("p2");
      }
      if(!_file_effCov_trg_staterr_mu_minus) {
        _file_effCov_trg_staterr_mu_minus = new TFile(Form("%s/src/CMGTools/WMass/python/postprocessing/data/leptonSF/new2016_madeSummer2018/systEff_trgmu_minus_mu.root",_cmssw_base_.c_str()),"read");
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


// float _muonTriggerSF_2l_trigMatch(int requiredCharge, // pass positive or negative number, depending on what you want 
// 				  float matchedTrgObjMuPt_l1, float matchedTrgObjMuPt_l2,
// 				  int pdgId1, int pdgId2,  // used to decide which lepton has the required charge
// 				  float pt1, float pt2,
// 				  float eta1, float eta2
// 				  ) {

//   int pdgId = 0;
//   float pt  = 0.0;
//   float eta = 0.0;
//   // muon (negative charge) has positive pdgId, antimuon (postive charge) has negative pdgId
//   // so, product of charge and pdgId_n must be negative to use pdgId_n and not the pther pdgId_n'
//   if (requiredCharge * pdgId1 < 0) {
//     // use lep 1
//     pdgId = pdgId1;
//     pt    = pt1;
//     eta   = eta1;
//     if (matchedTrgObjMuPt_l1 < 0.0) return 0;  // if no match to trigger, discard events
//   } else {
//     // use lep 2
//     pdgId = pdgId2;
//     pt    = pt2;
//     eta   = eta2;
//     if (matchedTrgObjMuPt_l2 < 0.0) return 0;  // if no match to trigger, discard events
//   } 

//   return _get_muonSF_selectionToTrigger(pdgId, pt, eta, requiredCharge);

// }


// bool triggerMatch(int requiredCharge, // pass positive or negative number, depending on what you want 
// 		  float matchedTrgObjMuPt_l1, float matchedTrgObjMuPt_l2,
// 		  int pdgId1, int pdgId2  // used to decide which lepton has the required charge
// 		  ) {

//   int pdgId = 0;
//   // muon (negative charge) has positive pdgId, antimuon (postive charge) has negative pdgId
//   // so, product of charge and pdgId_n must be negative to use pdgId_n and not the pther pdgId_n'
//   if (requiredCharge * pdgId1 < 0) {
//     // use lep 1
//     if (matchedTrgObjMuPt_l1 < 0.0) return 0;  // if no match to trigger, discard events
//     else                            return 1;
//   } else {
//     // use lep 2
//     pdgId = pdgId2;
//     if (matchedTrgObjMuPt_l2 < 0.0) return 0;  // if no match to trigger, discard events
//     else                            return 1;
//   }

// }



// float triggerSFforChargedLepton(int requiredCharge, // pass positive or negative number, depending on what you want 
// 				int pdgid1, int pdgid2,  // used to decide which lepton has the required charge
// 				float pt1, float pt2,
// 				float eta1, float eta2
// 				) {

//   int pdgId = 0;
//   float pt  = 0.0;
//   float eta = 0.0;
  
//   // pdgID > 0 for negative leptons
//   if (requiredCharge * pdgid1 < 0) {
//     pdgId = pdgid1;
//     pt    = pt1;
//     eta   = eta1;
//   } else {
//     pdgId = pdgid2;
//     pt    = pt2;
//     eta   = eta2;
//   }

//   return _get_muonSF_selectionToTrigger(pdgId, pt, eta, requiredCharge);					

// }

// float triggerSFforChargedLeptonMatchingTrigger(int requiredCharge, // pass positive or negative number, depending on what you want 
// 					       float matchedTrgObjMuPt_l1, float matchedTrgObjMuPt_l2,
// 					       int pdgid1, int pdgid2,  // used to decide which lepton has the required charge
// 					       float pt1, float pt2,
// 					       float eta1, float eta2
// 					       ) {

//   float trigMatchPt = 0.0;
//   int pdgId = 0;
//   float pt  = 0.0;
//   float eta = 0.0;

//   float trigMatchPt_other = 0.0;
//   int pdgId_other = 0;
//   float pt_other  = 0.0;
//   float eta_other = 0.0;
  
//   // pdgID > 0 for negative leptons
//   if (requiredCharge * pdgid1 < 0) {
//     pdgId = pdgid1;
//     pt    = pt1;
//     eta   = eta1;    
//     trigMatchPt = matchedTrgObjMuPt_l1;
//     pdgId_other = pdgid2;
//     pt_other    = pt2;
//     eta_other   = eta2;    
//     trigMatchPt_other = matchedTrgObjMuPt_l2;
//   } else {
//     pdgId = pdgid2;
//     trigMatchPt = matchedTrgObjMuPt_l2;
//     pt    = pt2;
//     eta   = eta2;
//     pdgId_other = pdgid1;
//     pt_other    = pt1;
//     eta_other   = eta1;    
//     trigMatchPt_other = matchedTrgObjMuPt_l1;
//   }

//   // try using 1 if both lepton match trigger (efficiency for trigger in 2 lepton phase space is ~ 100%)
//   if (trigMatchPt > 0.0 and trigMatchPt_other > 0.0) {
//     return 1;
//   } else {
//     if (trigMatchPt > 0.0)
//       return _get_muonSF_selectionToTrigger(pdgId, pt, eta, requiredCharge);					
//     else
//       return _get_muonSF_selectionToTrigger(pdgId_other, pt_other, eta_other, -1*requiredCharge);	 
//   }

// }

// float triggerSFforChargedLeptonMatchingTriggerV2(int requiredCharge, // pass positive or negative number, depending on what you want 
// 						 bool matchTrigger_l1, bool matchTrigger_l2,
// 						 int pdgid1, int pdgid2,  // used to decide which lepton has the required charge
// 						 float pt1, float pt2,
// 						 float eta1, float eta2						 
// 						 ) {

//   bool thisLepMatchesTrigger = false;
//   int pdgId = 0;
//   float pt  = 0.0;
//   float eta = 0.0;

//   bool otherLepMatchesTrigger = false;
//   int pdgId_other = 0;
//   float pt_other  = 0.0;
//   float eta_other = 0.0;
  
//   bool isSameSign = (pdgid1*pdgid2 > 0) ? true : false;
  
//   // pdgID > 0 for negative leptons
//   if (isSameSign || requiredCharge * pdgid1 < 0) {
//     pdgId = pdgid1;
//     pt    = pt1;
//     eta   = eta1;    
//     thisLepMatchesTrigger = matchTrigger_l1;
//     pdgId_other = pdgid2;
//     pt_other    = pt2;
//     eta_other   = eta2;    
//     otherLepMatchesTrigger = matchTrigger_l2;
//   } else {
//     pdgId = pdgid2;
//     pt    = pt2;
//     eta   = eta2;
//     thisLepMatchesTrigger = matchTrigger_l2;
//     pdgId_other = pdgid1;
//     pt_other    = pt1;
//     eta_other   = eta1;    
//     otherLepMatchesTrigger = matchTrigger_l1;
//   }

//   // try using 1 if both lepton match trigger (efficiency for trigger in 2 lepton phase space is ~ 100%)
//   if (thisLepMatchesTrigger and otherLepMatchesTrigger) {
//     return 1;
//   } else {
//     int sfCharge = thisLepMatchesTrigger ? requiredCharge : (-1*requiredCharge);
//     if (isSameSign) sfCharge = (pdgid1 > 0) ? -1 : 1; // for same sign uses the charge of the leading for the SF

//     if (thisLepMatchesTrigger)
//       return _get_muonSF_selectionToTrigger(pdgId, pt, eta, sfCharge);					
//     else if (otherLepMatchesTrigger)
//       return _get_muonSF_selectionToTrigger(pdgId_other, pt_other, eta_other, sfCharge);	 
//     else
//       return 0.0;  // this should not happen, but just in case
//   }

// }

// float triggerSFforChargedLeptonMatchingTriggerWlike(bool isOddEvent, // pass positive or negative number, depending on what you want 
// 						    bool matchTrigger_l1, bool matchTrigger_l2,
// 						    int pdgid1, int pdgid2,  // used to decide which lepton has the required charge
// 						    float pt1, float pt2,
// 						    float eta1, float eta2						 
// 						    ) {

//   int requiredCharge = isOddEvent ? 1 : -1;

//   bool thisLepMatchesTrigger = false;
//   int pdgId = 0;
//   float pt  = 0.0;
//   float eta = 0.0;

//   bool otherLepMatchesTrigger = false;
//   int pdgId_other = 0;
//   float pt_other  = 0.0;
//   float eta_other = 0.0;
  
//   bool isSameSign = (pdgid1*pdgid2 > 0) ? true : false;
  
//   // pdgID > 0 for negative leptons
//   if (isSameSign || requiredCharge * pdgid1 < 0) {
//     pdgId = pdgid1;
//     pt    = pt1;
//     eta   = eta1;    
//     thisLepMatchesTrigger = matchTrigger_l1;
//     pdgId_other = pdgid2;
//     pt_other    = pt2;
//     eta_other   = eta2;    
//     otherLepMatchesTrigger = matchTrigger_l2;
//   } else {
//     pdgId = pdgid2;
//     pt    = pt2;
//     eta   = eta2;
//     thisLepMatchesTrigger = matchTrigger_l2;
//     pdgId_other = pdgid1;
//     pt_other    = pt1;
//     eta_other   = eta1;    
//     otherLepMatchesTrigger = matchTrigger_l1;
//   }

//   // consider only case where given lepton matches trigger
//   if (isSameSign) {
//     if (thisLepMatchesTrigger and otherLepMatchesTrigger) {      
//       // charge is known, but what about which eta-pt ?
//       // this is a minor component for Z with SS in MC, might just randomly choose one
//       // use leading for now to be simple, it is really a small component
//       return _get_muonSF_selectionToTrigger(pdgId, pt, eta, -1 * pdgId);					    
//     } else if (thisLepMatchesTrigger) {
//       return _get_muonSF_selectionToTrigger(pdgId, pt, eta, -1 * pdgId);					    
//     } else if (otherLepMatchesTrigger) {
//       return _get_muonSF_selectionToTrigger(pdgId_other, pt_other, eta_other, -1 * pdgId_other);
//     } else {
//       return 0; // should not happen, but just in case
//     }
//   } else {
//     if (thisLepMatchesTrigger) {
//       int sfCharge = requiredCharge;
//       return _get_muonSF_selectionToTrigger(pdgId, pt, eta, sfCharge);					    
//     } else {
//       // if lepton does not match trigger, we decide to reject the event
//       return 0.0; 
//     }
//   }

// }


// bool triggerMatchV2(int requiredCharge, // pass positive or negative number, depending on what you want 
// 		    float matchedTrgObjMuPt_l1, float matchedTrgObjMuPt_l2,
// 		    float matchedTrgObjTkMuPt_l1, float matchedTrgObjTkMuPt_l2,
// 		    int pdgId1, int pdgId2  // used to decide which lepton has the required charge
// 		    ) {

//   // muon (negative charge) has positive pdgId, antimuon (postive charge) has negative pdgId
//   // so, product of charge and pdgId_n must be negative to use pdgId_n and not the pther pdgId_n'
//   if (requiredCharge * pdgId1 < 0) {
//     // use lep 1
//     return (matchedTrgObjMuPt_l1 > 0.0 || matchedTrgObjTkMuPt_l1 > 0.0) ? 1 : 0; 
//   } else {
//     // use lep 2
//     return (matchedTrgObjMuPt_l2 > 0.0 || matchedTrgObjTkMuPt_l2 > 0.0) ? 1 : 0; 
//   }

// }


// bool triggerMatchWlike(bool isOddEvent, // pass positive or negative number, depending on what you want 
// 		       float matchedTrgObjMuPt_l1, float matchedTrgObjMuPt_l2,
// 		       float matchedTrgObjTkMuPt_l1, float matchedTrgObjTkMuPt_l2,
// 		       int pdgId1, int pdgId2  // used to decide which lepton has the required charge
// 		       ) {

//   if (pdgId1 * pdgId2 > 0) {
//     // same sign case
//     return ((matchedTrgObjMuPt_l1 > 0.0 || matchedTrgObjTkMuPt_l1 > 0.0) || (matchedTrgObjMuPt_l2 > 0.0 || matchedTrgObjTkMuPt_l2));
//   }

//   int requiredCharge = isOddEvent ? 1 : -1;
//   // muon (negative charge) has positive pdgId, antimuon (postive charge) has negative pdgId
//   // so, product of charge and pdgId_n must be negative to use pdgId_n and not the pther pdgId_n'
//   if (requiredCharge * pdgId1 < 0) {
//     // use lep 1
//     return (matchedTrgObjMuPt_l1 > 0.0 || matchedTrgObjTkMuPt_l1 > 0.0) ? 1 : 0; 
//   } else {
//     // use lep 2
//     return (matchedTrgObjMuPt_l2 > 0.0 || matchedTrgObjTkMuPt_l2 > 0.0) ? 1 : 0; 
//   }

// }

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
