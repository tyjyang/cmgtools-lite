#ifndef FUNCTIONS_WMASS_H
#define FUNCTIONS_WMASS_H

#include "TROOT.h"
#include "TFile.h"
#include "TH2.h"
#include "TH3.h"
#include "TF1.h"
#include "TH2Poly.h"
#include "TSpline.h"
#include "TCanvas.h"
#include "TGraphAsymmErrors.h"
#include "TLorentzVector.h"
#include "EgammaAnalysis/ElectronTools/src/EnergyScaleCorrection_class.cc"

#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <cstdlib> //as stdlib.h      
#include <cstdio>
#include <string>

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

float returnChargeVal(float val1, int ch1, float val2, int ch2, ULong64_t evt){

    float retVal = -999.;

    if (evt%2 ) retVal = ch1 > 0 ? val1 : val2; //odd event numbers 
    else        retVal = ch1 > 0 ? val2 : val1; //even event numbers 

    return retVal;

}

float mt_wlike(float pt1, float phi1, int ch1, float pt2, float phi2, int ch2, float met, float phimet, ULong64_t evt) {
  float ptL  = returnChargeVal(pt1,ch1,pt2,ch2,evt);
  float phiL = returnChargeVal(phi1,ch1,phi2,ch2,evt);

  float ptl   = returnChargeVal(pt1,ch2,pt2,ch1,evt);
  float phil = returnChargeVal(phi1,ch2,phi2,ch1,evt);
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

float wpt_slope_weight(float wpt, float offset, float slope){
    if(wpt > 20.) return 1.;
    float weight = offset + slope*wpt;
    return weight;
}

float prefireJetsWeight(float eta){
    float feta = fabs(eta);
    if (feta < 2.0) return 0.9986;
    if (feta < 2.2) return 0.9880;
    if (feta < 2.3) return 0.9875;
    if (feta < 2.4) return 0.9871;
    if (feta < 2.5) return 0.9865;
    return 1.;
}

float prefireJetsWeight_2l(float eta1, float eta2) {

  float pf1 = prefireJetsWeight(eta1);
  float pf2 = prefireJetsWeight(eta2);
  return (pf1 + pf2 - pf1*pf2);

}


//-------------------

TFile *_file_ratio_FSRoverNoFSR_etaPt_mu = NULL;
TH2D * ratio_FSRoverNoFSR_etaPt_mu = NULL;
TFile *_file_ratio_FSRoverNoFSR_etaPt_el = NULL;
TH2D * ratio_FSRoverNoFSR_etaPt_el = NULL;

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

TFile *_file_recoToMedium_leptonSF_el = NULL;
TH2F *_histo_recoToMedium_leptonSF_el = NULL;
TFile *_file_recoToLoose_leptonSF_el = NULL;
TH2F *_histo_recoToLoose_leptonSF_el = NULL;
TFile *_file_elereco_leptonSF_gsf = NULL;
TH2F *_histo_elereco_leptonSF_gsf = NULL;

float _get_electronSF_recoToCustomTight(int pdgid, float pt, float eta, float var) {

  if (_cmssw_base_ == "") {
    cout << "Setting _cmssw_base_ to environment variable CMSSW_BASE" << endl;
    _cmssw_base_ = getEnvironmentVariable("CMSSW_BASE");
  }
    
  if (!_histo_elereco_leptonSF_gsf) {
    _file_elereco_leptonSF_gsf = new TFile(Form("%s/src/CMGTools/WMass/python/postprocessing/data/leptonSF/EGM2D_eleGSF.root",_cmssw_base_.c_str()),"read");
    _histo_elereco_leptonSF_gsf = (TH2F*)(_file_elereco_leptonSF_gsf->Get("EGamma_SF2D"));
    _histo_elereco_leptonSF_gsf->Smooth(1,"k3a");
  }

  if (!_histo_recoToMedium_leptonSF_el) {
    _file_recoToMedium_leptonSF_el = new TFile(Form("%s/src/CMGTools/WMass/python/postprocessing/data/leptonSF/EGM2D_eleCutBasedMediumWP.root",_cmssw_base_.c_str()),"read");
    _histo_recoToMedium_leptonSF_el = (TH2F*)(_file_recoToMedium_leptonSF_el->Get("EGamma_SF2D"));
    _histo_recoToMedium_leptonSF_el->Smooth(1,"k3a");
  }

  if (!_histo_recoToLoose_leptonSF_el) { 
    _file_recoToLoose_leptonSF_el = new TFile(Form("%s/src/CMGTools/WMass/python/postprocessing/data/leptonSF/EGM2D_eleCutBasedLooseWP.root",_cmssw_base_.c_str()),"read");
    _histo_recoToLoose_leptonSF_el = (TH2F*)(_file_recoToLoose_leptonSF_el->Get("EGamma_SF2D"));
    _histo_recoToLoose_leptonSF_el->Smooth(1,"k3a");
  }

  if(abs(pdgid)==11) {
    TH2F *histMedium = _histo_recoToMedium_leptonSF_el;
    TH2F *histLoose = _histo_recoToLoose_leptonSF_el;
    int etabin = std::max(1, std::min(histMedium->GetNbinsX(), histMedium->GetXaxis()->FindFixBin(eta)));
    int ptbin  = std::max(1, std::min(histMedium->GetNbinsY(), histMedium->GetYaxis()->FindFixBin(pt)));
    float out = 0;
    if(fabs(eta)<1.479) out = histLoose->GetBinContent(etabin,ptbin)+var*histLoose->GetBinError(etabin,ptbin);
    else out = histMedium->GetBinContent(etabin,ptbin)+var*histMedium->GetBinError(etabin,ptbin);

    TH2F *hist = _histo_elereco_leptonSF_gsf;
    etabin = std::max(1, std::min(hist->GetNbinsX(), hist->GetXaxis()->FindFixBin(eta)));
    ptbin  = std::max(1, std::min(hist->GetNbinsY(), hist->GetYaxis()->FindFixBin(pt)));
    out *= (hist->GetBinContent(etabin,ptbin)+var*(hist->GetBinError(etabin,ptbin) + 0.01*((pt<20) || (pt>80))));

    return out;
  }

  return 0;

}

TFile *_file_trigger_wmass_leptonEff_el = NULL;
TH2F *_histo_trigger_wmass_leptonEff_DATA_el = NULL;
TH2F *_histo_trigger_wmass_leptonEff_MC_el = NULL;

float _get_electronEff_HLT(float pt, float eta, bool isData) {
  if (!_histo_trigger_wmass_leptonEff_DATA_el || !_histo_trigger_wmass_leptonEff_MC_el) {
    _file_trigger_wmass_leptonEff_el = new TFile(Form("%s/src/CMGTools/WMass/python/postprocessing/data/leptonSF/new2016_madeSummer2018/electrons_trigger_pt30to55.root",_cmssw_base_.c_str()),"read");
    _histo_trigger_wmass_leptonEff_DATA_el = (TH2F*)(_file_trigger_wmass_leptonEff_el->Get("EGamma_EffData2D"));
    _histo_trigger_wmass_leptonEff_MC_el = (TH2F*)(_file_trigger_wmass_leptonEff_el->Get("EGamma_EffMC2D")) ;
    //_histo_trigger_wmass_leptonEff_DATA_el->Smooth(1,"k3a");
    //_histo_trigger_wmass_leptonEff_MC_el->Smooth(1,"k3a");
  }
  TH2F *hist = isData ? _histo_trigger_wmass_leptonEff_DATA_el : _histo_trigger_wmass_leptonEff_MC_el;
  int etabin = std::max(1, std::min(hist->GetNbinsX(), hist->GetXaxis()->FindFixBin(eta)));
  int ptbin  = std::max(1, std::min(hist->GetNbinsY(), hist->GetYaxis()->FindFixBin(pt)));
  float out = hist->GetBinContent(etabin,ptbin);
  return out;
}


TFile *_file_trigger_wmass_leptonSF_el = NULL;
TH2F *_histo_trigger_wmass_leptonSF_el = NULL;
TFile *_file_trigger_ee0p1_wmass_leptonSF_el = NULL;
TH2F *_histo_trigger_ee0p1_wmass_leptonSF_el = NULL;
TFile *_file_reco_wmass_leptonSF_el = NULL;
TH2F *_histo_reco_wmass_leptonSF_el = NULL;
TFile *_file_fullid_wmass_leptonSF_el = NULL;
TH2F *_histo_fullid_wmass_leptonSF_el = NULL;
TFile *_file_fullid_ee0p1_wmass_leptonSF_el = NULL;
TH2F *_histo_fullid_ee0p1_wmass_leptonSF_el = NULL;
TFile *_file_cluster_wmass_leptonSF_el = NULL;
TH2F *_histo_cluster_wmass_leptonSF_el = NULL;

float _get_electronSF_anyStep(float pt, float eta, int step, bool geterr=false) {
  if (_cmssw_base_ == "") {
    cout << "Setting _cmssw_base_ to environment variable CMSSW_BASE" << endl;
    _cmssw_base_ = getEnvironmentVariable("CMSSW_BASE");
  }
    
  if (!_histo_trigger_wmass_leptonSF_el) {
    _file_trigger_wmass_leptonSF_el = new TFile(Form("%s/src/CMGTools/WMass/python/postprocessing/data/leptonSF/new2016_madeSummer2018/etaptSmooth_electrons_trigger_30_55_onlyErf.root",_cmssw_base_.c_str()),"read");
    _histo_trigger_wmass_leptonSF_el = (TH2F*)(_file_trigger_wmass_leptonSF_el->Get("scaleFactor"));
  }
  if (!_histo_trigger_ee0p1_wmass_leptonSF_el) {
    _file_trigger_ee0p1_wmass_leptonSF_el = new TFile(Form("%s/src/CMGTools/WMass/python/postprocessing/data/leptonSF/new2016_madeSummer2018/electrons_trigger_endcap0p1.root",_cmssw_base_.c_str()),"read");
    _histo_trigger_ee0p1_wmass_leptonSF_el = (TH2F*)(_file_trigger_ee0p1_wmass_leptonSF_el->Get("scaleFactor"));
  }
  if (!_histo_reco_wmass_leptonSF_el) {
    _file_reco_wmass_leptonSF_el = new TFile(Form("%s/src/CMGTools/WMass/python/postprocessing/data/leptonSF/new2016_madeSummer2018/electrons_reco_pt30to45.root",_cmssw_base_.c_str()),"read");
    _histo_reco_wmass_leptonSF_el = (TH2F*)(_file_reco_wmass_leptonSF_el->Get("EGamma_SF2D"));
  }
  if (!_histo_fullid_wmass_leptonSF_el) {
    _file_fullid_wmass_leptonSF_el = new TFile(Form("%s/src/CMGTools/WMass/python/postprocessing/data/leptonSF/new2016_madeSummer2018/etaptSmooth_electrons_fullID_V2_pt25to55.root",_cmssw_base_.c_str()),"read");
    _histo_fullid_wmass_leptonSF_el = (TH2F*)(_file_fullid_wmass_leptonSF_el->Get("Graph2D_from_scaleFactor_smoothedByGraph"));
  }
  if (!_histo_fullid_ee0p1_wmass_leptonSF_el) {
    _file_fullid_ee0p1_wmass_leptonSF_el = new TFile(Form("%s/src/CMGTools/WMass/python/postprocessing/data/leptonSF/new2016_madeSummer2018/electrons_fullID_V2_endcap0p1.root",_cmssw_base_.c_str()),"read");
    _histo_fullid_ee0p1_wmass_leptonSF_el = (TH2F*)(_file_fullid_ee0p1_wmass_leptonSF_el->Get("EGamma_SF2D"));
  }
  if (!_histo_cluster_wmass_leptonSF_el) {
    _file_cluster_wmass_leptonSF_el = new TFile(Form("%s/src/CMGTools/WMass/python/postprocessing/data/leptonSF/new2016_madeSummer2018/l1EG_eff.root",_cmssw_base_.c_str()),"read");
    _histo_cluster_wmass_leptonSF_el = (TH2F*)(_file_cluster_wmass_leptonSF_el->Get("l1EG_eff"));
  }

  TH2F *hist = 0;
  if (step==1) {
    hist = fabs(eta)<1.566 ? _histo_trigger_wmass_leptonSF_el : _histo_trigger_ee0p1_wmass_leptonSF_el;
  }
  else if (step==2) hist = _histo_reco_wmass_leptonSF_el;
  else if (step==3) {
    hist = fabs(eta)<1.566 ? _histo_fullid_wmass_leptonSF_el : _histo_fullid_ee0p1_wmass_leptonSF_el;
  }
  else if (step==4) hist = _histo_cluster_wmass_leptonSF_el;
  else {
    std::cout << "Step " << step << " of the efficiency corrections not foreseen. Returning 0" << std::endl;
    return 0.;
  }

  int etabin = std::max(1, std::min(hist->GetNbinsX(), hist->GetXaxis()->FindFixBin(eta)));
  int ptbin  = std::max(1, std::min(hist->GetNbinsY(), hist->GetYaxis()->FindFixBin(pt)));
  float out = geterr ? hist->GetBinError(etabin,ptbin) : hist->GetBinContent(etabin,ptbin);
  return out;
}

float eleSF_HLT(float pt, float eta) {
  return _get_electronSF_anyStep(pt,eta,1);
}

float eleSF_HLT_2lfriends(float matchpt1, float sf1, float matchpt2, float sf2) {
  if (matchpt1>-1 && matchpt2>-1) {
    return 0.5*(sf1+sf2);
  } else if (matchpt1>-1) 
    return sf1;
  else
    return sf2;
}

float eleSF_HLT_2l(float matchpt1, float pt1, float eta1, float matchpt2, float pt2, float eta2) {
  if (matchpt1>-1 && matchpt2>-1) {
    float sf1 = _get_electronSF_anyStep(pt1,eta1,1);
    float sf2 = _get_electronSF_anyStep(pt2,eta2,1);
    return 0.5*(sf1+sf2);
  } else if (matchpt1>-1) 
    return _get_electronSF_anyStep(pt1,eta1,1);
  else
    return _get_electronSF_anyStep(pt2,eta2,1);
}

float eleSF_HLT_2lAvg(float matchpt1, float pt1, float eta1, float matchpt2, float pt2, float eta2) {
  if (matchpt1>-1 && matchpt2>-1) {
    double eff1_DATA = _get_electronEff_HLT(pt1,eta1,true);
    double eff2_DATA = _get_electronEff_HLT(pt2,eta2,true);
    float sf1 = _get_electronSF_anyStep(pt1,eta1,1);
    float sf2 = _get_electronSF_anyStep(pt2,eta2,1);
    return (eff1_DATA*sf1 + eff2_DATA*sf2) / (eff1_DATA+eff2_DATA);
  } else if (matchpt1>-1) {
    return _get_electronSF_anyStep(pt1,eta1,1);
  }
  else {
    return _get_electronSF_anyStep(pt2,eta2,1);
  }
}

float eleEff_HLT_2l_DATA(float pt1, float eta1, float pt2, float eta2) {
  float eff1 = _get_electronEff_HLT(pt1,eta1,true);
  float eff2 = _get_electronEff_HLT(pt2,eta2,true);
  return eff1 + eff2 - eff1*eff2;
}

float eleSF_HLT_2lComb(float matchpt1, float pt1, float eta1, float matchpt2, float pt2, float eta2) {
  double eff1_DATA = _get_electronEff_HLT(pt1,eta1,true);
  double eff2_DATA = _get_electronEff_HLT(pt2,eta2,true);
  double eff1_MC = _get_electronEff_HLT(pt1,eta1,false);
  double eff2_MC = _get_electronEff_HLT(pt2,eta2,false);
  double eff_DATA,eff_MC;
  if (matchpt1>-1 && matchpt2>-1) {
    eff_DATA = std::min(1., eff1_DATA + eff2_DATA - eff1_DATA*eff2_DATA);
    eff_MC =   std::min(1., eff1_MC + eff2_MC - eff1_MC*eff2_MC);
  } else if (matchpt1>-1) {
    eff_DATA = eff1_DATA * (1-eff2_DATA);
    eff_MC = eff1_MC * (1-eff2_MC);
    //    return _get_electronSF_anyStep(pt1,eta1,1);
  } else {
    eff_DATA = eff2_DATA * (1-eff1_DATA);
    eff_MC = eff2_MC * (1-eff1_MC);
    //    return _get_electronSF_anyStep(pt2,eta2,1);
  }
  if (eff_DATA>1 || eff_MC>1) 
    std::cout << pt1 << " "<< eta1 << "   2 = " << pt2 << " " << eta2 << "  eff  = " << eff_DATA << " " << eff_MC << std::endl;
  return eff_DATA / eff_MC;
}


float eleSF_GSFReco(float pt, float eta) {
  return _get_electronSF_anyStep(pt,eta,2);
}

float eleSF_FullID(float pt, float eta) {
  return _get_electronSF_anyStep(std::min(pt,float(45.)),eta,3);
}

float eleSF_Clustering(float pt, float eta) {
  return _get_electronSF_anyStep(pt,eta,4);
}

float eleSF_L1Eff(float pt, float eta, bool geterr=false) {
  float sf;
  if (fabs(eta)<1.479 || pt<35) {
    //if (fabs(eta)<1.479) {
    sf = geterr ? 0.0 : 1.0;
  }
  else sf = _get_electronSF_anyStep(pt,eta,4,geterr);
  return sf;
}

float eleSF_L1Eff_2l(float pt1, float eta1, float pt2, float eta2) {
  float eta,pt;
  if (fabs(eta1) > fabs(eta2)) {
    eta = eta1;
    pt = pt1;
  } else {
    eta = eta2;
    pt = pt2;
  }
  if (fabs(eta)<1.479 || pt<35) return 1.;
  return _get_electronSF_anyStep(pt,eta,4);
}


float _lepSF(int pdgId, float pt, float eta, float sf1, float sf2, float sf3, int nSigma=0) {
  float abseta = fabs(eta);
  float syst=0;
  float sf4=1.0; float sf4_err=0.0;
  if (abs(pdgId)==11) {
    sf4 = eleSF_L1Eff(pt,eta);
    sf4_err = eleSF_L1Eff(pt,eta,true);
    if (abseta<1)          syst = 0.006;
    else if (abseta<1.479) syst = 0.008;
    else if (abseta<2)     syst = 0.013;
    else                   syst = 0.016;
    if (abseta>1.479)      syst = hypot(syst,sf4_err);
  } else if (abs(pdgId)==13) {
    if (abseta<1)          syst = 0.002;
    else if (abseta<1.5)   syst = 0.004;
    else                   syst = 0.014;
  }
  double ret = sf1*sf2*sf3*sf4 + nSigma*syst;  
  if (ret <= 0.0) {
    // cout << "Error: returning bad SF: " << ret << "    Please check! Returning 1" << endl;
    // cout << ">>> pdgId, pt, eta, sf1, sf2, sf3, nsigma --> " << pdgId << ", " << pt << ", " << eta << ", " << sf1 << ", " << sf2 << ", " << sf3 << ",  " << nSigma << endl;
    return 1;
  } else {
    return ret;  
  }

}

float lepSF(int pdgId, float pt, float eta, float sf1, float sf2, float sf3) {
  return _lepSF(pdgId,pt,eta,sf1,sf2,sf3,0);
}

float lepSFRelUp(int pdgId, float pt, float eta, float sf1, float sf2, float sf3) {
  return _lepSF(pdgId,pt,eta,sf1,sf2,sf3, 1)/_lepSF(pdgId,pt,eta,sf1,sf2,sf3,0);
}

float lepSFRelDn(int pdgId, float pt, float eta, float sf1, float sf2, float sf3) {
  return _lepSF(pdgId,pt,eta,sf1,sf2,sf3, -1)/_lepSF(pdgId,pt,eta,sf1,sf2,sf3,0);
}

#include "TRandom.h"
TRandom3 *rng = NULL;
EnergyScaleCorrection_class *calibrator = NULL;

float ptCorr(float pt, float eta, float phi, float r9, int run, int isData, ULong64_t eventNumber) {

  if(!calibrator) calibrator = new EnergyScaleCorrection_class("CMGTools/WMass/python/postprocessing/data/leptonScale/el/Legacy2016_07Aug2017_FineEtaR9_ele",0);

  if(!isData) {
    if(!rng) rng = new TRandom3();
    // use eventNumber as seed, otherwise each time the function is called for the same event, the smearer produce a different pt values
    // the better solution would be to have the corrected pt in the friend trees
    rng->SetSeed(eventNumber); // make it really random across different jobs    
  }

  if (isData) return pt *  calibrator->ScaleCorrection(run,fabs(eta)<1.479,r9,fabs(eta),pt);
  else {
    float smear = calibrator->getSmearingSigma(run,fabs(eta)<1.479,r9,fabs(eta),pt,0,0);
    return pt * ( 1.0 + smear * rng->Gaus());
  }
}

TFile *_file_residualcorr_scale = NULL;
TH2D *_histo_residualcorr_scale = NULL;

float residualScale(float pt, float eta, int isData, const char *fileCorr="../postprocessing/data/leptonScale/el/plot_dm_diff.root") {
  if(!isData) return 1.;

  if(!_histo_residualcorr_scale) {
    _file_residualcorr_scale = new TFile(fileCorr);
    _histo_residualcorr_scale = (TH2D*)(_file_residualcorr_scale->Get("histSmooth"));
  }
  
  TH2D *hist = _histo_residualcorr_scale;
  int etabin  = std::max(1, std::min(hist->GetNbinsX(), hist->GetXaxis()->FindFixBin(fabs(eta))));
  int ptbin = std::max(1, std::min(hist->GetNbinsY(), hist->GetYaxis()->FindFixBin(pt)));
  
  const float MZ0 = 91.1876;
  float scale = 1. - hist->GetBinContent(etabin,ptbin)/MZ0/sqrt(2.);
  if (scale < 0) {
    cout << "WARNING in residualScale() function: scale < 0 --> returning 0." << endl;
    return 0;
  } else {
    return scale;
  }

}

float ptElFull(float pt, float eta, int nSigma=0) {

  if (nSigma == 0) return pt;

  // THE FOLLOWING USES STD EGAMMA SYSTEMATICS (W/O RESIDUAL CORRECTIONS, INTEGRATED IN ETA/PT)
  /*
  float relSyst=0.;
  if(fabs(eta)<1.0) relSyst = 0.0015;  
  else if(fabs(eta)<1.479) relSyst = 0.005;  
  else relSyst = 0.01; 
  return (1.+nSigma*relSyst) * pt;
  */
  // the following uses private residual corrections of AN-17-340
  if (_cmssw_base_ == "") {
    cout << "Setting _cmssw_base_ to environment variable CMSSW_BASE" << endl;
    _cmssw_base_ = getEnvironmentVariable("CMSSW_BASE");
  }
  float syst = 1-residualScale(pt,eta,1,Form("%s/src/CMGTools/WMass/python/postprocessing/data/leptonScale/el/plot_dm_diff_closure_smoothWithSpline.root",_cmssw_base_.c_str()));
  return (1. + nSigma*syst) * pt;

}

float ptElFullUp(float pt, float eta) {
  return ptElFull(pt,eta,1);
}

float ptElFullDn(float pt, float eta) {
  return ptElFull(pt,eta,-1);
}

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

//-------------------
// electrons
float ptEleScaleUncorr0Up(float pt, float eta) {
  return ptScaleUncorr(pt, eta, 11, 0, true);
}

float ptEleScaleUncorr0Dn(float pt, float eta) {
  return ptScaleUncorr(pt, eta, 11, 0, false);
}

float ptEleScaleUncorr1Up(float pt, float eta) {
  return ptScaleUncorr(pt, eta, 11, 1, true);
}

float ptEleScaleUncorr1Dn(float pt, float eta) {
  return ptScaleUncorr(pt, eta, 11, 1, false);
}

float ptEleScaleUncorr2Up(float pt, float eta) {
  return ptScaleUncorr(pt, eta, 11, 2, true);
}

float ptEleScaleUncorr2Dn(float pt, float eta) {
  return ptScaleUncorr(pt, eta, 11, 2, false);
}

float ptEleScaleUncorr3Up(float pt, float eta) {
  return ptScaleUncorr(pt, eta, 11, 3, true);
}

float ptEleScaleUncorr3Dn(float pt, float eta) {
  return ptScaleUncorr(pt, eta, 11, 3, false);
}

//===============================================


float ptMuFull(float pt, float eta, int nSigma=0) {

  if (nSigma == 0) return pt;

  // THE FOLLOWING USES THE DERIVED MUON NON-CLOSURE
  if (_cmssw_base_ == "") {
    cout << "Setting _cmssw_base_ to environment variable CMSSW_BASE" << endl;
    _cmssw_base_ = getEnvironmentVariable("CMSSW_BASE");
  }
  float syst = 1-residualScaleMu(pt,eta,1,Form("%s/src/CMGTools/WMass/data/muonscale/scale_correction_nonclosure_mu.root",_cmssw_base_.c_str()));
  //float syst = 0.004; // test
  return (1. + nSigma*syst) * pt;
  
}

float ptMuFullUp(float pt, float eta) {
  return ptMuFull(pt,eta,1);
}

float ptMuFullDn(float pt, float eta) {
  return ptMuFull(pt,eta,-1);
}

//==================================================

float getSmearedVar(float var, float smear, ULong64_t eventNumber, int isData, bool smearOnlyMC=false) {

  if (smearOnlyMC && isData) return var;

  if(!rng) rng = new TRandom3();
  // use eventNumber as seed, otherwise each time the function is called for the same event, the smearer produce a different smeared value
  rng->SetSeed(eventNumber); // make it really random across different jobs    
  return var * ( 1.0 + smear * rng->Gaus());

}

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

//============================

TFile *_file_eleTrigger_EBeta0p1 = NULL;
TH2F *_histo_eleTrigger_EBeta0p1 = NULL;
TFile *_file_eleTrigger_EEeta0p1 = NULL;
TH2F *_histo_eleTrigger_EEeta0p1 = NULL;
TFile *_file_eleID_EBeta0p1 = NULL;
TH2F *_histo_eleID_EBeta0p1 = NULL;
TFile *_file_eleID_EEeta0p1 = NULL;
TH2F *_histo_eleID_EEeta0p1 = NULL;

float _get_electronSF_TriggerAndID(int pdgid, float pt, float eta) {

  if (_cmssw_base_ == "") {
    cout << "Setting _cmssw_base_ to environment variable CMSSW_BASE" << endl;
    _cmssw_base_ = getEnvironmentVariable("CMSSW_BASE");
  }

  // get all histograms 
  if (!_histo_eleTrigger_EBeta0p1) {
    _file_eleTrigger_EBeta0p1 = new TFile(Form("%s/src/CMGTools/WMass/python/postprocessing/data/leptonSF/new2016_madeSummer2018/etaptSmooth_electrons_trigger_30_55_onlyErf.root",_cmssw_base_.c_str()),"read");
    _histo_eleTrigger_EBeta0p1 = (TH2F*)(_file_eleTrigger_EBeta0p1->Get("scaleFactor"));
  }

  if (!_histo_eleTrigger_EEeta0p1) {
    _file_eleTrigger_EEeta0p1 = new TFile(Form("%s/src/CMGTools/WMass/python/postprocessing/data/leptonSF/new2016_madeSummer2018/electrons_trigger_endcap0p1.root",_cmssw_base_.c_str()),"read");
    _histo_eleTrigger_EEeta0p1 = (TH2F*)(_file_eleTrigger_EEeta0p1->Get("scaleFactor"));
  }

  if (!_histo_eleID_EBeta0p1) {
    _file_eleID_EBeta0p1 = new TFile(Form("%s/src/CMGTools/WMass/python/postprocessing/data/leptonSF/new2016_madeSummer2018/TnPstuff/electron/elFullID_V2_granular/smoothEfficiency_electrons_fullID_V2_granular.root",_cmssw_base_.c_str()),"read");
    _histo_eleID_EBeta0p1 = (TH2F*)(_file_eleID_EBeta0p1->Get("scaleFactor"));
  }

  if (!_histo_eleID_EEeta0p1) {
    _file_eleID_EEeta0p1 = new TFile(Form("%s/src/CMGTools/WMass/python/postprocessing/data/leptonSF/new2016_madeSummer2018/TnPstuff/electron/elFullID_V2_endcap0p1/smoothEfficiency_electrons_fulID_V2_endcap0p1.root",_cmssw_base_.c_str()),"read");
    _histo_eleID_EEeta0p1 = (TH2F*)(_file_eleID_EEeta0p1->Get("scaleFactor"));
  }

  int etabinTrigger = 0;
  int etabinID = 0;
  int ptbinTrigger = 0;
  int ptbinID = 0;
  TH2F *histTrigger = NULL;
  TH2F *histID = NULL;

  if (fabs(pdgid) == 11) {
    if (fabs(eta) >= 1.566) {
      histTrigger = _histo_eleTrigger_EEeta0p1;
      histID = _histo_eleID_EEeta0p1;      
    } else {
      histTrigger = _histo_eleTrigger_EBeta0p1;
      histID = _histo_eleID_EBeta0p1;      
    }
    int etabinTrigger = std::max(1, std::min(histTrigger->GetNbinsX(), histTrigger->GetXaxis()->FindFixBin(eta)));
    int ptbinTrigger  = std::max(1, std::min(histTrigger->GetNbinsY(), histTrigger->GetYaxis()->FindFixBin(pt)));
    int etabinID = std::max(1, std::min(histID->GetNbinsX(), histID->GetXaxis()->FindFixBin(eta)));
    int ptbinID  = std::max(1, std::min(histID->GetNbinsY(), histID->GetYaxis()->FindFixBin(pt)));
    return histTrigger->GetBinContent(etabinTrigger,ptbinTrigger) * histID->GetBinContent(etabinID,ptbinID);
  }

  return 0;

}

//float triggerSF_2l_histo(float l1pt, float l1eta, float l11pass, float l12pass, float l2pt, float l2eta, float l21pass, float l22pass){
//  float weight = -999.;
//
//  bool l1pass = (l11pass > -1. || l12pass > -1.);
//  bool l2pass = (l21pass > -1. || l22pass > -1.);
//
//  float l1sf = _get_muonSF_selectionToTrigger(13,l1pt,l1eta);
//  float l2sf = _get_muonSF_selectionToTrigger(13,l2pt,l2eta);
//
//  int randomize = 0;
//
//  if      ( l1pass && !l2pass) weight = l1sf;
//  else if (!l1pass &&  l2pass) weight = l2sf;
//  else if ( l1pass &&  l2pass) {
//    if   (!randomize) weight = (l1sf+l2sf)/2.;
//    else {
//      if(!rng) rng = new TRandom3();
//      float randy = rng->Uniform(-1.,1.);
//      if (randy < 0.) weight = l1sf;
//      else            weight = l2sf;
//    }
//  }
//
//  //else return -999.; // this should never happen
//
//  return weight;
//}

//float triggerSF_2l_byCharge(float l1pt, float l1eta, float l11pass, float l12pass, float l1charge, int charge){
//  float weight = -999.;
//
//  bool l1pass = (l11pass > -1. || l12pass > -1.);
//
//  float l1sf = 1.;
//  if (l1pass) l1sf = _get_muonSF_selectionToTrigger(13,l1pt,l1eta);
//
//  if (l1charge*charge > 0) weight = l1sf;
//  else weight = 1.;
//
//  return weight;
//}
//
//float triggerSF_1l_histo(float l1pt, float l1eta){
//
//  float l1sf = _get_muonSF_selectionToTrigger(13,l1pt,l1eta);
//
//  return l1sf;
//}

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

float lepSFRmuUp(int pdgId, float pt, float eta, int charge, float sf2) {

  // in this function SF1 is obtained inside here, while SF3 (reco SF, which would generally be 1) is treated as prefiring SF
  // at this level we can have electrons, because the cuts are applied afterwards, and this would return bad numbers
  if (fabs(pdgId) != 13) return 0;
  double trigSF = _get_muonSF_selectionToTrigger(pdgId,pt,eta,charge);
  double prefireSF = prefireJetsWeight(eta);
  return _lepSF(pdgId,pt,eta,trigSF,sf2,prefireSF, 1)/_lepSF(pdgId,pt,eta,trigSF,sf2,prefireSF,0);

}

float lepSFRmuDn(int pdgId, float pt, float eta, int charge, float sf2) {

  // in this function SF1 is obtained inside here, while SF3 (reco SF, which would generally be 1) is treated as prefiring SF
  // at this level we can have electrons, because the cuts are applied afterwards, and this would return bad numbers
  if (fabs(pdgId) != 13) return 0;
  double trigSF = _get_muonSF_selectionToTrigger(pdgId,pt,eta,charge);
  double prefireSF = prefireJetsWeight(eta);
  return _lepSF(pdgId,pt,eta,trigSF,sf2,prefireSF, -1)/_lepSF(pdgId,pt,eta,trigSF,sf2,prefireSF,0);

}


float _muonTriggerSF_2l_trigMatch(int requiredCharge, // pass positive or negative number, depending on what you want 
				  float matchedTrgObjMuPt_l1, float matchedTrgObjMuPt_l2,
				  int pdgId1, int pdgId2,  // used to decide which lepton has the required charge
				  float pt1, float pt2,
				  float eta1, float eta2
				  ) {

  int pdgId = 0;
  float pt  = 0.0;
  float eta = 0.0;
  // muon (negative charge) has positive pdgId, antimuon (postive charge) has negative pdgId
  // so, product of charge and pdgId_n must be negative to use pdgId_n and not the pther pdgId_n'
  if (requiredCharge * pdgId1 < 0) {
    // use lep 1
    pdgId = pdgId1;
    pt    = pt1;
    eta   = eta1;
    if (matchedTrgObjMuPt_l1 < 0.0) return 0;  // if no match to trigger, discard events
  } else {
    // use lep 2
    pdgId = pdgId2;
    pt    = pt2;
    eta   = eta2;
    if (matchedTrgObjMuPt_l2 < 0.0) return 0;  // if no match to trigger, discard events
  } 

  return _get_muonSF_selectionToTrigger(pdgId, pt, eta, requiredCharge);

}


bool triggerMatch(int requiredCharge, // pass positive or negative number, depending on what you want 
		  float matchedTrgObjMuPt_l1, float matchedTrgObjMuPt_l2,
		  int pdgId1, int pdgId2  // used to decide which lepton has the required charge
		  ) {

  int pdgId = 0;
  // muon (negative charge) has positive pdgId, antimuon (postive charge) has negative pdgId
  // so, product of charge and pdgId_n must be negative to use pdgId_n and not the pther pdgId_n'
  if (requiredCharge * pdgId1 < 0) {
    // use lep 1
    if (matchedTrgObjMuPt_l1 < 0.0) return 0;  // if no match to trigger, discard events
    else                            return 1;
  } else {
    // use lep 2
    pdgId = pdgId2;
    if (matchedTrgObjMuPt_l2 < 0.0) return 0;  // if no match to trigger, discard events
    else                            return 1;
  }

}



float triggerSFforChargedLepton(int requiredCharge, // pass positive or negative number, depending on what you want 
				int pdgid1, int pdgid2,  // used to decide which lepton has the required charge
				float pt1, float pt2,
				float eta1, float eta2
				) {

  int pdgId = 0;
  float pt  = 0.0;
  float eta = 0.0;
  
  // pdgID > 0 for negative leptons
  if (requiredCharge * pdgid1 < 0) {
    pdgId = pdgid1;
    pt    = pt1;
    eta   = eta1;
  } else {
    pdgId = pdgid2;
    pt    = pt2;
    eta   = eta2;
  }

  return _get_muonSF_selectionToTrigger(pdgId, pt, eta, requiredCharge);					

}

float triggerSFforChargedLeptonMatchingTrigger(int requiredCharge, // pass positive or negative number, depending on what you want 
					       float matchedTrgObjMuPt_l1, float matchedTrgObjMuPt_l2,
					       int pdgid1, int pdgid2,  // used to decide which lepton has the required charge
					       float pt1, float pt2,
					       float eta1, float eta2
					       ) {

  float trigMatchPt = 0.0;
  int pdgId = 0;
  float pt  = 0.0;
  float eta = 0.0;

  float trigMatchPt_other = 0.0;
  int pdgId_other = 0;
  float pt_other  = 0.0;
  float eta_other = 0.0;
  
  // pdgID > 0 for negative leptons
  if (requiredCharge * pdgid1 < 0) {
    pdgId = pdgid1;
    pt    = pt1;
    eta   = eta1;    
    trigMatchPt = matchedTrgObjMuPt_l1;
    pdgId_other = pdgid2;
    pt_other    = pt2;
    eta_other   = eta2;    
    trigMatchPt_other = matchedTrgObjMuPt_l2;
  } else {
    pdgId = pdgid2;
    trigMatchPt = matchedTrgObjMuPt_l2;
    pt    = pt2;
    eta   = eta2;
    pdgId_other = pdgid1;
    pt_other    = pt1;
    eta_other   = eta1;    
    trigMatchPt_other = matchedTrgObjMuPt_l1;
  }

  // try using 1 if both lepton match trigger (efficiency for trigger in 2 lepton phase space is ~ 100%)
  if (trigMatchPt > 0.0 and trigMatchPt_other > 0.0) {
    return 1;
  } else {
    if (trigMatchPt > 0.0)
      return _get_muonSF_selectionToTrigger(pdgId, pt, eta, requiredCharge);					
    else
      return _get_muonSF_selectionToTrigger(pdgId_other, pt_other, eta_other, -1*requiredCharge);	 
  }

}

float triggerSFforChargedLeptonMatchingTriggerV2(int requiredCharge, // pass positive or negative number, depending on what you want 
						 bool matchTrigger_l1, bool matchTrigger_l2,
						 int pdgid1, int pdgid2,  // used to decide which lepton has the required charge
						 float pt1, float pt2,
						 float eta1, float eta2
						 ) {

  bool thisLepMatchesTrigger = false;
  int pdgId = 0;
  float pt  = 0.0;
  float eta = 0.0;

  bool otherLepMatchesTrigger = false;
  int pdgId_other = 0;
  float pt_other  = 0.0;
  float eta_other = 0.0;
  
  // pdgID > 0 for negative leptons
  if (requiredCharge * pdgid1 < 0) {
    pdgId = pdgid1;
    pt    = pt1;
    eta   = eta1;    
    thisLepMatchesTrigger = matchTrigger_l1;
    pdgId_other = pdgid2;
    pt_other    = pt2;
    eta_other   = eta2;    
    otherLepMatchesTrigger = matchTrigger_l2;
  } else {
    pdgId = pdgid2;
    pt    = pt2;
    eta   = eta2;
    thisLepMatchesTrigger = matchTrigger_l2;
    pdgId_other = pdgid1;
    pt_other    = pt1;
    eta_other   = eta1;    
    otherLepMatchesTrigger = matchTrigger_l1;
  }

  // try using 1 if both lepton match trigger (efficiency for trigger in 2 lepton phase space is ~ 100%)
  if (thisLepMatchesTrigger and otherLepMatchesTrigger) {
    return 1;
  } else {
    if (thisLepMatchesTrigger)
      return _get_muonSF_selectionToTrigger(pdgId, pt, eta, requiredCharge);					
    else if (otherLepMatchesTrigger)
      return _get_muonSF_selectionToTrigger(pdgId_other, pt_other, eta_other, -1*requiredCharge);	 
    else
      return 0.0;  // this should not happen, but just in case
  }

}


bool isOddEvent(ULong64_t evt) {

  return (evt%2) ? 1 : 0;       

}

bool isEvenEvent(ULong64_t evt) {

  return (evt%2) ? 0 : 1;       

}

#endif
