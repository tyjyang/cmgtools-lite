#ifndef FAKERATE_CC
#define FAKERATE_CC

#include <TH2.h>
#include <TH2D.h>
#include <TFile.h>
#include <TMath.h>
#include <TF1.h>
#include <cmath>
#include <iostream>
#include <string>
#include <map>
#include <cstdlib> //as stdlib.h         
#include <cstdio>

TH2 * helicityFractions_0 = 0;
TH2 * helicityFractions_L = 0;
TH2 * helicityFractions_R = 0;

TH2 * angular_c = 0;
TH2 * angular_0 = 0;
TH2 * angular_1 = 0;
TH2 * angular_2 = 0;
TH2 * angular_3 = 0;
TH2 * angular_4 = 0;
TH2 * angular_5 = 0;
TH2 * angular_6 = 0;
TH2 * angular_7 = 0;

// Index 0 is for nominal, 1-4 is for variation of FR linear fit parameters (1,2 for up/dn offset, 3,4 for slope)
// 5-6 is used for normalization systematic variation implemented through an eta(-pt) dependent lnN nuisance. 
// Actually index 5 or 6 are used to store the histogram for FR without any variation, to avoid clashes with FRi_mu[0] or FRi_el[0] used when running the nominal FR,
// since the variation is implemented inside fakeRateWeight_promptRateCorr_1l_i_smoothed 
// see for example: w-helicity-13TeV/wmass_e/fakerate-vars/fakeRate-frdata-e-normup.txt 
// Index 7 is the shape variation when awayJetPt > 45 (one could add other systematics, adding other indices)
// Index 8 and 9 accounts for FR obtained subtracting EWK scaled up and down by 1 sigma of their cross section
// Index 10 and 11 to scale  Slope Up/Down with EWK Down/Up (more EWK means lower FR, so it makes sense to decrease the slope to create a band around nominal)
// Keep array index larger than the number of index you will use, so you don't have to increase it everytime you add another index
TH2 * FR_mu = 0;
TH2 * FRi_mu[15] = {0};  
TH2 * FR_el = 0;
TH2 * FRi_el[15] = {0};

// FR for QCD MC, needed not to clash with that on data (above) in case they are used together
TH2 * FR_mu_qcdmc = 0;
TH2 * FRi_mu_qcdmc[5] = {0}; 
TH2 * FR_el_qcdmc = 0;
TH2 * FRi_el_qcdmc[5] = {0};

// prompt rate
TH2 * PR_mu = 0;
TH2 * PRi_mu[15] = {0};
TH2 * PR_el = 0;
TH2 * PRi_el[15] = {0};

TH2 * FRcorrection = 0;
TH2 * FRcorrection_i[5];

// following few lines are obsolete
TH2 * FRnormSyst_el = 0;  // normalization systematics for FR asaf pt and eta (equivalent to lnN nuisance parameter for each template bin)
TH2 * FRnormSyst_el_i[3] = {0}; // will only have up and down, but let's use 3 to include a nominal variation
TH2 * FRnormSyst_mu = 0;  // normalization systematics for FR asaf pt and eta (equivalent to lnN nuisance parameter for each template bin)
TH2 * FRnormSyst_mu_i[3] = {0};

bool loadFRHisto(const std::string &histoName, const std::string file, const char *name) {

  TH2 **histo = 0, **hptr2 = 0;
  TH2 * FR_temp = 0;
  TH2 * FR_normSyst_temp = 0;
  TH2 * PR_temp = 0;
  TH2 * FRcorr_temp = 0;
  if (histoName == "FR_mu")  { histo = & FR_mu;  hptr2 = & FRi_mu[0]; }
  else if (histoName == "FR_mu_qcdmc")  { histo = & FR_mu_qcdmc;  hptr2 = & FRi_mu_qcdmc[0]; }
  else if (histoName == "FR_el")  { histo = & FR_el;  hptr2 = & FRi_el[0]; }
  else if (histoName == "FR_el_qcdmc")  { histo = & FR_el_qcdmc;  hptr2 = & FRi_el_qcdmc[0]; }
  else if (histoName == "PR_el")  { histo = & PR_el;  hptr2 = & PRi_el[0]; }
  else if (histoName == "PR_mu")  { histo = & PR_mu;  hptr2 = & PRi_mu[0]; }
  else if (TString(histoName).BeginsWith("FR_mu_i")) {histo = & FR_temp; hptr2 = & FRi_mu[TString(histoName).ReplaceAll("FR_mu_i","").Atoi()];}
  else if (TString(histoName).BeginsWith("FR_el_i")) {histo = & FR_temp; hptr2 = & FRi_el[TString(histoName).ReplaceAll("FR_el_i","").Atoi()];}
  else if (TString(histoName).Contains("helicityFractions_0")) { histo = & helicityFractions_0; }
  else if (TString(histoName).Contains("helicityFractions_L")) { histo = & helicityFractions_L; }
  else if (TString(histoName).Contains("helicityFractions_R")) { histo = & helicityFractions_R; }
  else if (TString(histoName).Contains("angular_c")) { histo = & angular_c; }
  else if (TString(histoName).Contains("angular_0")) { histo = & angular_0; }
  else if (TString(histoName).Contains("angular_1")) { histo = & angular_1; }
  else if (TString(histoName).Contains("angular_2")) { histo = & angular_2; }
  else if (TString(histoName).Contains("angular_3")) { histo = & angular_3; }
  else if (TString(histoName).Contains("angular_4")) { histo = & angular_4; }
  else if (TString(histoName).Contains("angular_5")) { histo = & angular_5; }
  else if (TString(histoName).Contains("angular_6")) { histo = & angular_6; }
  else if (TString(histoName).Contains("angular_7")) { histo = & angular_7; }
  else if (TString(histoName).BeginsWith("PR_mu_i")) {histo = & PR_temp; hptr2 = & PRi_mu[TString(histoName).ReplaceAll("PR_mu_i","").Atoi()];}
  else if (TString(histoName).BeginsWith("PR_el_i")) {histo = & PR_temp; hptr2 = & PRi_el[TString(histoName).ReplaceAll("PR_el_i","").Atoi()];}
  else if (histoName == "FRnormSyst_el")  { histo = & FRnormSyst_el;  hptr2 = & FRnormSyst_el_i[0];}
  else if (histoName == "FRnormSyst_mu")  { histo = & FRnormSyst_mu;  hptr2 = & FRnormSyst_mu_i[0];}
  else if (TString(histoName).BeginsWith("FRnormSyst_el_i")) {histo = & FR_normSyst_temp; hptr2 = & FRnormSyst_el_i[TString(histoName).ReplaceAll("FRnormSyst_el_i","").Atoi()];}
  else if (TString(histoName).BeginsWith("FRnormSyst_mu_i")) {histo = & FR_normSyst_temp; hptr2 = & FRnormSyst_mu_i[TString(histoName).ReplaceAll("FRnormSyst_mu_i","").Atoi()];}
  else if (histoName == "FR_correction")  { histo = & FRcorrection; hptr2 = & FRcorrection_i[0]; }
  else if (TString(histoName).BeginsWith("FR_correction_i")) {histo = & FRcorr_temp; hptr2 = & FRcorrection_i[TString(histoName).ReplaceAll("FRcorrection_i","").Atoi()];}
  if (histo == 0)  {
    std::cerr << "ERROR: histogram " << histoName << " is not defined in fakeRate.cc." << std::endl;
    return 0;
  }

  TFile *f = TFile::Open(file.c_str());
  if (*histo != 0) {
    if (std::string(name) != (*histo)->GetName()) {
      std::cerr << "WARNING 1: overwriting histogram " << (*histo)->GetName() << std::endl;
    } else {
      TH2* hnew = (TH2*) f->Get(name);
      if (hnew == 0 || hnew->GetNbinsX() != (*histo)->GetNbinsX() || hnew->GetNbinsY() != (*histo)->GetNbinsY()) {
	std::cerr << "WARNING 2: overwriting histogram " << (*histo)->GetName() << std::endl;
      } else {
	bool fail = false;
	for (int ix = 1; ix <= (*histo)->GetNbinsX(); ++ix) {
	  for (int iy = 1; iy <= (*histo)->GetNbinsX(); ++iy) {
	    if ((*histo)->GetBinContent(ix,iy) != hnew->GetBinContent(ix,iy)) {
	      fail = true; break;
	    }
	  }
	}
	if (fail) std::cerr << "WARNING 3: overwriting histogram " << (*histo)->GetName() << std::endl;
      }
    }
    delete *histo;
  }
  if (f->Get(name) == 0) {
    std::cerr << "ERROR: could not find " << name << " in " << file << std::endl;
    *histo = 0;
  } else {
    *histo = (TH2*) f->Get(name)->Clone(name);
    (*histo)->SetDirectory(0);
    if (hptr2) *hptr2 = *histo;
  }
  f->Close();
  return histo != 0;

}

//===============================================================

float getFakeRatenormWeight(float lpt, float leta, int lpdgId, int ivar = 0) {

  if (!ivar) return 1.0; // this is a special case: the "nominal" FRnormSyst TH1 was never created (it would be FR_el or FR_mu)
  // only the variations are considered here
  // ivar == 1: up variation
  // ivar == 2: down variation

  if (FRnormSyst_el_i[ivar] == 0 && FRnormSyst_mu_i[ivar] == 0) {
    std::cout << "Error in getFakeRatenormWeight(): both histograms are 0. Returning 0" << std::endl;	
    return 0;
  } 

  int fid = abs(lpdgId);   
  TH2 *hist = (fid == 11 ? FRnormSyst_el_i[ivar] : FRnormSyst_mu_i[ivar]);
  if (hist == 0) {
    std::cout << "Error in getFakeRatenormWeight(): hist == 0. Returning 0" << std::endl;	
    return 0;
  } 

  // do we need a weight just as a function of eta or pt as well?
  float absleta = std::abs(leta);
  int etabin = std::max(1, std::min(hist->GetNbinsX(), hist->GetXaxis()->FindBin(absleta)));
  int ptbin  = std::max(1, std::min(hist->GetNbinsX(), hist->GetYaxis()->FindBin(lpt)));
  
  // TH1 content is a fraction, like 30%, so we return 1 +/- var
  float var = hist->GetBinContent(etabin, ptbin);
  return (ivar == 1) ? (1 + var) : (1 - var);
  
}

//===============================================================

TFile * _file_EleFR_pt_eta_fitPol2 = nullptr;
TH2D* _histo_EleFR_pt_eta_fitPol2_data_nomi = nullptr;
TH2D* _histo_EleFR_pt_eta_fitPol2_data_ewkUp = nullptr; // used to define a conservative uncertainty
TH2D* _histo_EleFR_pt_eta_fitPol2_data_ewkDn = nullptr; // used to define a conservative uncertainty

//===============================================================

float fakeRateWeight_promptRateCorr_1l_i_smoothed(float lpt, float leta, int lpdgId, bool passWP, int iFR=0, int iPR=0) { //, int expected_pdgId=11) {

  // formula for fake rate including effect of prompt rate
  //
  // Let LNT denote the region passing loose but not tight selection, T the region passing the tight selection.
  // Let p and f denote the prompt and fake lepton rate respectively.
  // Then:
  // N(QCD in T) = f/(p-f) * (p*N(NLT) - (1-p)*N(T))
  // second term is negative by definition of p)
  // If p=1, then N(QCD in T) = f/(1-f) * N(NLT), which is the formula used in function fakeRateWeight_1l_i_smoothed()

  // this is needed to extend lepton FR above 50 GeV, need to use pol2 fit (or binned FR, but pol2 is better)
  if (!_file_EleFR_pt_eta_fitPol2) {
    string cmssw_base_path = "";
    char* _env_var_ptr = getenv("CMSSW_BASE");
    if (_env_var_ptr == nullptr) {
      cout << "Error: environment variable CMSSW_BASE not found. Exit" << endl;
      exit(EXIT_FAILURE);
    } else {
      cmssw_base_path = string(_env_var_ptr);      
    }    
    _file_EleFR_pt_eta_fitPol2 = new TFile(Form("%s/src/CMGTools/WMass/data/fakerate/frAndPrSmoothed_el_frPol2_prPol1.root",cmssw_base_path.c_str()),"read");
    _histo_EleFR_pt_eta_fitPol2_data_nomi = (TH2D*) _file_EleFR_pt_eta_fitPol2->Get("fr_pt_eta_data");  // histo with smoothed FR values and uncertainties (stat.only)
    _histo_EleFR_pt_eta_fitPol2_data_ewkUp = (TH2D*) _file_EleFR_pt_eta_fitPol2->Get("fr_pt_eta_data_subScaledUpEWKMC"); // no uncertainty stored
    _histo_EleFR_pt_eta_fitPol2_data_ewkDn = (TH2D*) _file_EleFR_pt_eta_fitPol2->Get("fr_pt_eta_data_subScaledDownEWKMC"); // no uncertainty stored
    if (!_histo_EleFR_pt_eta_fitPol2_data_nomi|| ! _histo_EleFR_pt_eta_fitPol2_data_ewkUp || !_histo_EleFR_pt_eta_fitPol2_data_ewkDn) {
      cout << "Warning: I couldn't get some histograms from file, with FR smoothed with pol2" << endl;
      exit(EXIT_FAILURE);
    }
  }

  // double fpt = lpt; // what was the purpose of this line?
  double feta = std::fabs(leta); int fid = abs(lpdgId); 

  // int fAbsExpected_pdgId = abs(expected_pdgId);
  // if (fid != fAbsExpected_pdgId) {
  //   return 0;
  // }

  if (FRi_el[iFR] == 0 && FRi_mu[iFR] == 0) {
    // this is the case where the histogram was not loaded correctly (one is 0 because you use the other flavour)
    std::cout << "Error in fakeRateWeight_promptRateCorr_1l_i_smoothed: hist_fr == 0. Returning 0" << std::endl;	
    return 0;
  } 

  if (PRi_el[iPR] == 0 && PRi_mu[iPR] == 0) {
    // as above
    std::cout << "Error in fakeRateWeight_promptRateCorr_1l_i_smoothed: hist_pr == 0. Returning 0" << std::endl;	
    return 0;
  }


  TH2 *hist_fr = (fid == 11 ? FRi_el[iFR] : FRi_mu[iFR]);
  TH2 *hist_pr = (fid == 11 ? PRi_el[iPR] : PRi_mu[iPR]);
  if (hist_fr == 0 || hist_pr == 0) {
    // this is the case where you expect electrons but get a muon, or viceversa
    // Indeed, selection is evaluated as 1 or 0 multiplying the event weight in TTree::Draw(...), so you potentially have all flavours here
    // do not issue warning messages here, unless it is for testing
    //std::cout << "Error in fakeRateWeight_promptRateCorr_1l_i_smoothed: hist_fr == 0. It seems the flavour is not what you expect. Returning 0" << std::endl;	
    return 0;
  }

  Bool_t hasNegativeEta = (hist_fr->GetXaxis()->GetBinLowEdge(1) < 0) ? true : false;
  int etabin = std::max(1, std::min(hist_fr->GetNbinsX(), hist_fr->GetXaxis()->FindBin(hasNegativeEta ? leta : feta)));
  int nFRfitParam = hist_fr->GetNbinsY();
  int nPRfitParam = hist_pr->GetNbinsY();
  // FR
  float p0 = hist_fr->GetBinContent(etabin, 1);
  float p1 = hist_fr->GetBinContent(etabin, 2);
  float p2 = (nFRfitParam > 2) ? hist_fr->GetBinContent(etabin, 3) : 0.0;
  if      (iFR==1) p0 += hist_fr->GetBinError(etabin, 1);
  else if (iFR==2) p0 -= hist_fr->GetBinError(etabin, 1);
  // offset and slope are anticorrelated, when changing slope, have to move the offset as well
  else if (iFR==3 || iFR==10) {
    p1 += hist_fr->GetBinError(etabin, 2);
    p0 -= hist_fr->GetBinError(etabin, 1);
  } else if (iFR==4 || iFR==11) {
    p1 -= hist_fr->GetBinError(etabin, 2);
    p0 += hist_fr->GetBinError(etabin, 1);
  }
  // now PR
  // eta bin is typically the same as for fake rate, but let's allow the possibility that it is different
  etabin = std::max(1, std::min(hist_pr->GetNbinsX(), hist_pr->GetXaxis()->FindBin(hasNegativeEta ? leta : feta)));
  float p0_pr = hist_pr->GetBinContent(etabin, 1);
  float p1_pr = hist_pr->GetBinContent(etabin, 2);
  float p2_pr = (nPRfitParam > 2) ? hist_pr->GetBinContent(etabin, 3) : 0.0;

  if      (iPR==1) p0_pr += hist_pr->GetBinError(etabin, 1);
  else if (iPR==2) p0_pr -= hist_pr->GetBinError(etabin, 1);
  // offset and slope are anticorrelated, when changing slope, have to move the offset as well
  else if (iPR==3) {
    p1_pr += hist_pr->GetBinError(etabin, 2);
    p0_pr -= hist_pr->GetBinError(etabin, 1);
  } else if (iPR==4) {
    p1_pr -= hist_pr->GetBinError(etabin, 2);
    p0_pr += hist_pr->GetBinError(etabin, 1);
  }

  // Marc added the crop at pt=50, but for electrons I think it is not needed
  // I will make so to have it only for muon channel
  if (fid == 13 && lpt > 50.) lpt = 50.; // this helps because after 50 the FR is often not well modeled due to non linear behaviour.
  // could just be the limit of the method, because the ewk contribution gets larger
  // it works because the muon FR rises and then starts decreasing again or gets flat
  // for electrons, the FR decreases with pt, so making it flat above 50 is not good is not good
  // one could use the estimate from pol2 for points above 50 because it works extremely well, but then the pt slope variation must be defined in a different way 
  // one could just plug a temptative uncertainty for points above 50, given that for fakes we keep large regions of eta-pt bins uncorrelated in the fit

  // p2=0 if the saved histogram has only 2 bins
  // in case on is fitting with a straight line and not a parabola, the histogram might still have 3 bins for the parameters, but the last one should be set as 0
  // this saves backward compatibility
  float fr = p0    + p1   *lpt + p2   *lpt*lpt; 
  float pr = p0_pr + p1_pr*lpt + p2_pr*lpt*lpt;

  if (fid == 11 && lpt > 48) {
    Int_t ptbin_tmp = std::max(1, std::min(_histo_EleFR_pt_eta_fitPol2_data_nomi->GetXaxis()->FindFixBin(lpt), _histo_EleFR_pt_eta_fitPol2_data_nomi->GetNbinsX()));
    Int_t etabin_tmp = std::max(1, std::min(_histo_EleFR_pt_eta_fitPol2_data_nomi->GetYaxis()->FindFixBin(leta), _histo_EleFR_pt_eta_fitPol2_data_nomi->GetNbinsY()));
    if (iFR==3 || iFR==10) {
      // get larger value for FR (from EWK down + uncertainty of nominal (stat.uncertainties are basically the same for all the three histograms)
      fr = _histo_EleFR_pt_eta_fitPol2_data_ewkDn->GetBinContent(ptbin_tmp, etabin_tmp) + _histo_EleFR_pt_eta_fitPol2_data_nomi->GetBinError(ptbin_tmp, etabin_tmp);
    } else if (iFR==4 || iFR==11) {
      // get lower value for FR (from EWK Up - uncertainty of nominal (stat.uncertainties are basically the same for all the three histograms)
      fr = _histo_EleFR_pt_eta_fitPol2_data_ewkUp->GetBinContent(ptbin_tmp, etabin_tmp) - _histo_EleFR_pt_eta_fitPol2_data_nomi->GetBinError(ptbin_tmp, etabin_tmp);
    } else {
      fr = _histo_EleFR_pt_eta_fitPol2_data_nomi->GetBinContent(ptbin_tmp, etabin_tmp);    
    }
  }

  if (fid == 13 && pr > 0.98) pr = 0.98; // safety thing, not needed for electrons
  else if (pr > 1.0) pr = 1.0;  // just in case

  // implement an eta-pt dependent lnN nuisance parameter to account for normalization variations
  float FRnormWgt = 1.0; 
  if(fid == 11){
    //if      (iFR==5) FRnormWgt = getFakeRatenormWeight(lpt, feta, fid, 1);
    //else if (iFR==6) FRnormWgt = getFakeRatenormWeight(lpt, feta, fid, 2);
    if      (iFR==5) FRnormWgt = 1.05 + feta*0.06; // from 1.05 to 1.20
    else if (iFR==6) FRnormWgt = 0.95 - feta*0.06;
  }
  else {
    // for muons vary continuously in eta from 5% to 20% between eta = 0 and eta = 2.4
    if      (iFR==5) FRnormWgt = 1.05 + feta*0.0625;
    else if (iFR==6) FRnormWgt = 0.95 - feta*0.0625;
  }

  float weight;

  // for large pt, when using pol2 it can happen that FR > PR, but this was observed for pt > 100 GeV, which is far beyond the range we are interested
  // so in that case the weight can be safely set as 0, because those events are not used in the analysis
  if (pr < fr) {
    //std::cout << "### Error in weight: FR > PR. Please check!" << std::endl;
    //std::cout << " pt: " << lpt << " eta:" << leta << " pdgid: " << lpdgId << std::endl;
    return 0;
  } else if (pr == fr) {
    //std::cout << "### Error in weight: division by 0. Please check" << std::endl;
    //std::cout << " pt: " << lpt << " eta:" << leta << " pdgid: " << lpdgId << std::endl;
    return 0;
  }

  if (passWP) {
    // tight
    // returning a negative weight
    weight = FRnormWgt*fr*(pr-1)/(pr-fr); // pr=1 --> return 0
  } else {
    // not tight (but still loose)
    weight = FRnormWgt*fr*pr/(pr-fr);  // pr=1 --> return fr/(1-fr)
  }

  // if (weight != weight)   std::cout << "weight is NaN" << std::endl;
  // if (fabs(weight) > 10.) std::cout << "event with large weight: " << weight << " pt: " << lpt << " eta:" << leta << " pdgid: " << lpdgId << std::endl;

  return weight;

}

//==============================

float fakeRateWeight_1l_i_smoothed(float lpt, float leta, int lpdgId, bool passWP, int iFR=0) { //, int expected_pdgId=11) {
  if (!passWP) {
    double fpt = lpt; double feta = std::fabs(leta); int fid = abs(lpdgId); 
    // int fAbsExpected_pdgId = abs(expected_pdgId);
    // if (fid != fAbsExpected_pdgId) {
    //   return 0;
    // }
    if (FRi_el[iFR] == 0 and FRi_mu[iFR] == 0) {
      // this is the case where the histogram was not loaded correctly (one is 0 because you use the other flavour)
      std::cout << "Error in fakeRateWeight_1l_i_smoothed: hist == 0. Returning 0" << std::endl;	
      return 0;
    } 
    TH2 *hist = (fid == 11 ? FRi_el[iFR] : FRi_mu[iFR]);
    if (hist == 0) {
      // this is the case where you expect electrons but get a muon, or viceversa
      // Indeed, selection is evaluated as 1 or 0 multiplying the event weight in TTree::Draw(...), so you potentially have all flavours here
      // do not issue warnign mewssages here, unless it is for testing
      //std::cout << "Error in fakeRateWeight_1l_i_smoothed: hist == 0. Returning 0" << std::endl;	
      //std::cout << "pdg ID = " << lpdgId << std::endl;
      return 0;
    }
    Bool_t hasNegativeEta = (hist->GetXaxis()->GetBinLowEdge(1) < 0);
    int etabin = std::max(1, std::min(hist->GetNbinsX(), hist->GetXaxis()->FindBin(hasNegativeEta ? leta : feta)));
    float p0 = hist->GetBinContent(etabin, 1);
    float p1 = hist->GetBinContent(etabin, 2);
    if      (iFR==1) p0 += hist->GetBinError(etabin, 1);
    else if (iFR==2) p0 -= hist->GetBinError(etabin, 1);
    else if (iFR==3) p1 += hist->GetBinError(etabin, 2);
    else if (iFR==4) p1 -= hist->GetBinError(etabin, 2);
    float fr = p0 + p1*lpt;
    return fr/(1-fr);

  } else return 0;
}


float fakeRateWeight_1l_i_smoothed_FRcorr(float lpt, float leta, int lpdgId, bool passWP, int iFR=0, float var=-999, int iFRcorr=0) { //, int expected_pdgId=11) {
  if (!passWP) {
    double fpt = lpt; double feta = std::fabs(leta); int fid = abs(lpdgId); 
    // int fAbsExpected_pdgId = abs(expected_pdgId);
    // if (fid != fAbsExpected_pdgId) {
    //   return 0;
    // }
    if (FRi_el[iFR] == 0 and FRi_mu[iFR] == 0) {
      // this is the case where the histogram was not loaded correctly (one is 0 because you use the other flavour)
      std::cout << "Error in fakeRateWeight_1l_i_smoothed: hist == 0. Returning 0" << std::endl;	
      return 0;
    } 
    TH2 *hist = (fid == 11 ? FRi_el[iFR] : FRi_mu[iFR]);
    if (hist == 0) {
      // this is the case where you expect electrons but get a muon, or viceversa
      // Indeed, selection is evaluated as 1 or 0 multiplying the event weight in TTree::Draw(...), so you potentially have all flavours here
      // do not issue warnign mewssages here, unless it is for testing
      //std::cout << "Error in fakeRateWeight_1l_i_smoothed: hist == 0. Returning 0" << std::endl;	
      //std::cout << "pdg ID = " << lpdgId << std::endl;
      return 0;
    }
    int etabin = std::max(1, std::min(hist->GetNbinsX(), hist->GetXaxis()->FindBin(feta)));
    float p0 = hist->GetBinContent(etabin, 1);
    float p1 = hist->GetBinContent(etabin, 2);
    if      (iFR==1) p0 += hist->GetBinError(etabin, 1);
    else if (iFR==2) p0 -= hist->GetBinError(etabin, 1);
    else if (iFR==3) p1 += hist->GetBinError(etabin, 2);
    else if (iFR==4) p1 -= hist->GetBinError(etabin, 2);
    float fr = p0 + p1*lpt;

    // FR correction (iFRcorr > 0 not supported yet, could be used for variations)
    float FRcorrection = 1;
    TH2 *hist_FRcorr = 0;
    if (var > -999) {
      hist_FRcorr = FRcorrection_i[iFRcorr];
      if (hist_FRcorr == 0) {
	std::cout << "Error in fakeRateWeight_1l_i_smoothed_FRcorr: hist_FRcorr == 0. Applying no correction (i.e. 1)" << std::endl;
	//return 0;
      } else {
	int varbin = std::max(1, std::min(hist_FRcorr->GetNbinsX(), hist_FRcorr->GetXaxis()->FindBin(var))); 
	etabin = std::max(1, std::min(hist_FRcorr->GetNbinsY(), hist_FRcorr->GetYaxis()->FindBin(leta))); 
	FRcorrection = hist_FRcorr->GetBinContent(varbin,etabin); 
      }
    }

    return FRcorrection * fr/(1-fr);

  } else return 0;
}


float fakeRateWeight_1l_i(float lpt, float leta, int lpdgId, bool passWP, int iFR) {
  if (!passWP) {
    double fpt = lpt; double feta = std::fabs(leta); int fid = abs(lpdgId);
    TH2 *hist = (fid == 11 ? FRi_el[iFR] : FRi_mu[iFR]);
    if (hist == 0) return 0;
    int ptbin  = std::max(1, std::min(hist->GetNbinsX(), hist->GetXaxis()->FindBin(fpt)));
    int etabin = std::max(1, std::min(hist->GetNbinsY(), hist->GetYaxis()->FindBin(feta)));
    double fr = hist->GetBinContent(ptbin,etabin);
    return fr/(1-fr);
  } else return 0;
}

float fakeRateWeight_1l(float lpt, float leta, int lpdgId, bool passWP)
{
  return fakeRateWeight_1l_i(lpt, leta, lpdgId, passWP, 0);
}

float fetchFR_i(float l1pt, float l1eta, int l1pdgId, int iFR) 
{
    TH2 *hist1 = (abs(l1pdgId) == 11 ? FRi_el[iFR] : FRi_mu[iFR]);
    if (hist1 == 0) { std::cerr << "ERROR, missing FR for pdgId " << l1pdgId << ", iFR " << iFR << std::endl; }
    int ptbin1  = std::max(1, std::min(hist1->GetNbinsX(), hist1->GetXaxis()->FindBin(l1pt)));
    int etabin1 = std::max(1, std::min(hist1->GetNbinsY(), hist1->GetYaxis()->FindBin(std::abs(l1eta))));
    double fr1 = hist1->GetBinContent(ptbin1,etabin1);
    if (fr1 <= 0)  { std::cerr << "WARNING, FR is " << fr1 << " for " << hist1->GetName() << ", pt " << l1pt << " eta " << l1eta << std::endl; }
    return fr1;
}

   
TF1 * helicityFraction_0 = new TF1("helicityFraction_0", "3./4*(TMath::Sqrt(1-x*x))^2", -1., 1.);
TF1 * helicityFraction_L = new TF1("helicityFraction_L", "3./8.*(1-x)^2"              , -1., 1.);
TF1 * helicityFraction_R = new TF1("helicityFraction_R", "3./8.*(1+x)^2"              , -1., 1.);

float helicityWeight(float yw, float ptw, float costheta, int pol)
{

  float ayw = std::abs(yw);

  if (std::abs(costheta) > 1.) {
    std::cout << " found an event with weird cosTheta = " << costheta << std::endl;
    std::cout << " setting event weight to 0" << std::endl;
    return 0;
  }

  TH2 *hist_f0 = helicityFractions_0;
  TH2 *hist_fL = helicityFractions_L;
  TH2 *hist_fR = helicityFractions_R;

  // float yval  = std::abs(yw) > hist_f0->GetXaxis()->GetXmax() ? hist_f0->GetXaxis()->GetXmax() : yw;
  // float ptval = ptw > hist_f0->GetYaxis()->GetXmax() ? hist_f0->GetYaxis()->GetXmax() : ptw;

  int ywbin = std::max(1, std::min(hist_f0->GetNbinsX(), hist_f0->GetXaxis()->FindBin(ayw)));
  int ptbin = std::max(1, std::min(hist_f0->GetNbinsY(), hist_f0->GetYaxis()->FindBin(ptw)));

  float f0 = hist_f0->GetBinContent(ywbin, ptbin);
  float fL = hist_fL->GetBinContent(ywbin, ptbin);
  float fR = hist_fR->GetBinContent(ywbin, ptbin);

  float f0Term = helicityFraction_0->Eval(costheta);
  float fLTerm = helicityFraction_L->Eval(costheta);
  float fRTerm = helicityFraction_R->Eval(costheta);

  float weight = 0.;
  float max_weight = 4.;

   //std::cout << "normalization of the weights: " << f0*f0Term+fL*fLTerm+fR*fRTerm << std::endl;
   //weight = fR*fRTerm/(f0*f0Term+fL*fLTerm+fR*fRTerm);
   // if (weight >= max_weight) std::cout << " weight for fR larger than 4!!!!" << weight << std::endl;

  if      (pol == 0) return std::min( f0*f0Term/(f0*f0Term+fL*fLTerm+fR*fRTerm), max_weight);
  else if (pol == 1) return std::min( fL*fLTerm/(f0*f0Term+fL*fLTerm+fR*fRTerm), max_weight);
  else if (pol == 2) return std::min( fR*fRTerm/(f0*f0Term+fL*fLTerm+fR*fRTerm), max_weight);
        
  std::cout << "something went wrong in the helicity reweighting" << std::endl;
  return -99999.;

}

// ================================================================================
float angularWeight(float yw, float ptw, float costheta, float phi, int pol)
{

  float ayw = std::abs(yw);

  if (std::abs(costheta) > 1.) {
    std::cout << " found an event with weird cosTheta = " << costheta << std::endl;
    std::cout << " setting event weight to 0" << std::endl;
    return 0;
  }

  TH2 *hist_ac = angular_c;
  TH2 *hist_a0 = angular_0;
  TH2 *hist_a1 = angular_1;
  TH2 *hist_a2 = angular_2;
  TH2 *hist_a3 = angular_3;
  TH2 *hist_a4 = angular_4;
  TH2 *hist_a5 = angular_5;
  TH2 *hist_a6 = angular_6;
  TH2 *hist_a7 = angular_7;

  // float yval  = std::abs(yw) > hist_f0->GetXaxis()->GetXmax() ? hist_f0->GetXaxis()->GetXmax() : yw;
  // float ptval = ptw > hist_f0->GetYaxis()->GetXmax() ? hist_f0->GetYaxis()->GetXmax() : ptw;

  int ywbin = std::max(1, std::min(angular_0->GetNbinsX(), angular_0->GetXaxis()->FindBin(ayw)));
  int ptbin = std::max(1, std::min(angular_0->GetNbinsY(), angular_0->GetYaxis()->FindBin(ptw)));

  float fc = angular_c->GetBinContent(ywbin, ptbin);
  float f0 = angular_0->GetBinContent(ywbin, ptbin);
  float f1 = angular_1->GetBinContent(ywbin, ptbin);
  float f2 = angular_2->GetBinContent(ywbin, ptbin);
  float f3 = angular_3->GetBinContent(ywbin, ptbin);
  float f4 = angular_4->GetBinContent(ywbin, ptbin);
  float f5 = angular_5->GetBinContent(ywbin, ptbin);
  float f6 = angular_6->GetBinContent(ywbin, ptbin);
  float f7 = angular_7->GetBinContent(ywbin, ptbin);

  float theta = TMath::ACos(costheta);

  float fcTerm = (1.+costheta*costheta);
  float f0Term = 0.5*(1.-3.*costheta*costheta);
  float f1Term = TMath::Sin(2.*theta)*TMath::Cos(phi);
  float f2Term = 0.5*TMath::Sin(theta)*TMath::Sin(theta)*TMath::Cos(2.*phi);
  float f3Term = TMath::Sin(theta)*TMath::Cos(phi);
  float f4Term = costheta;
  float f5Term = TMath::Sin(theta)*TMath::Sin(theta)*TMath::Sin(2.*phi);
  float f6Term = TMath::Sin(2.*theta)*TMath::Sin(phi);
  float f7Term = TMath::Sin(theta)*TMath::Sin(phi);


  float weight = 0.;

  float normterm = 3./(16.*TMath::Pi())*(fc*fcTerm + f0*f0Term + f1*f1Term + f2*f2Term + f3*f3Term + f4*f4Term + f5*f5Term + f6*f6Term + f7*f7Term);

  if      (pol == -1) return fc*fcTerm/(normterm);
  else if (pol ==  0) return f0*f0Term/(normterm);
  else if (pol ==  1) return f1*f1Term/(normterm);
  else if (pol ==  2) return f2*f2Term/(normterm);
  else if (pol ==  3) return f3*f3Term/(normterm);
  else if (pol ==  4) return f4*f4Term/(normterm);
  else if (pol ==  5) return f5*f5Term/(normterm);
  else if (pol ==  6) return f6*f6Term/(normterm);
  else if (pol ==  7) return f7*f7Term/(normterm);
        
  std::cout << "something went wrong in the angular coefficient reweighting" << std::endl;
  return -99999.;

}


#endif
