#include "../interface/utility.h"
#include "../../functions.cc"
#include "../../w-mass-13TeV/functionsWMass.cc"

#define CHECK_EVERY_N 10000
#define N_MAX_ENTRIES_PER_SAMPLE 0 // for tests, use number <= 0 to use all events in each sample
#define FOLDER_IN_FILE 0

using namespace std;

//static float intLumi = 30.9 // 35.9 for muons, 30.9 for electrons, measured in 1/fb
//static vector<Double_t> eleEtaBinEdges_double = {0.0, 1.0, 1.479, 2.1, 2.5};
// static vector<Double_t> etaBinEdgesTemplate = {-2.5,-2.4,-2.3,-2.2,-2.1,-2.0,-1.9,-1.8,-1.7,-1.6,-1.566,-1.5,-1.4442,-1.4,-1.3,-1.2,-1.1,-1.0,-0.9,-0.8,-0.7,-0.6,-0.5,-0.4,-0.3,-0.2,-0.1,0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.4442,1.5,1.566,1.6,1.7,1.8,1.9,2.0,2.1,2.2,2.3,2.4,2.5};
// static vector<Double_t> ptBinEdgesTemplate = {30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45};

// static vector<Double_t> genEtaBinEdgesTemplate = {0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.1,2.2,2.3,2.4};
// //static vector<Double_t> genPtBinEdgesTemplate = {26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45};
// static vector<Double_t> genPtBinEdgesTemplate = {30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45};

static vector<Double_t> wptbins = {0.0, 1.971, 2.949, 3.838, 4.733, 5.674, 6.684, 7.781, 8.979, 10.303, 11.777, 13.435, 15.332, 17.525, 20.115, 23.245, 27.173, 32.414, 40.151, 53.858, 13000.0};
static Int_t nPDFweight = 60;


//=============================================================

void fillHistograms(const string& treedir = "./", 
		    const string& outdir = "./", 
		    const Sample& sample = Sample::wjets, 
		    TFile* outputFile = NULL,
		    const bool usePreFSRvar = true,  // if false, use DressedLepton
		    const string& nameMatch = ""
		    ) 
{

  gROOT->SetBatch(kTRUE);
  TH1::SetDefaultSumw2(); //all the following histograms will automatically call TH1::Sumw2()  

  cout << endl;
  cout << "================================================" << endl;
  cout << "Using integrated luminosity L = " << intLumi << endl;
  cout << endl;

  if (outputFile == NULL) {
    cout << "Error: file is NULL. please check. Exit ..." << endl;
    exit(EXIT_FAILURE);
  }

  outputFile->cd();

  bool isMuon = false;
  if (sample == Sample::wenujets) isMuon = false;
  if (sample == Sample::wmunujets) isMuon = true;


  TDirectory *dirSample = NULL;
  string sampleDir = getStringFromEnumSample(sample).c_str();
  cout << "Sample --> " << sampleDir << endl;
  if (FOLDER_IN_FILE) {
    if (outputFile->GetKey(sampleDir.c_str())) dirSample = outputFile->GetDirectory(sampleDir.c_str());
    else dirSample = outputFile->mkdir(sampleDir.c_str());
    dirSample->cd();
  }

  cout << endl;

  //============================
  // BINNING: will have to implement parsing from file

  vector<Double_t> etaBinEdgesTemplate = {-2.4, -2.1, -1.9,-1.7,-1.5,-1.3,-1.2,-1.1,-1.0,-0.9,
					  -0.8,-0.7,-0.6,-0.5,-0.4,-0.3,-0.2,-0.1,0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,
					  1.1,1.2,1.3,1.5,1.7,1.9, 2.1, 2.4};
  vector<Double_t> ptBinEdgesTemplate = {25,26,28,30,31.5,33,34.5,36,37.5,39.0,40.5,42,43.5,45,46.5};
  
  Double_t etaCutReco = etaBinEdgesTemplate.back(); 
  Double_t ptCutReco = ptBinEdgesTemplate.back(); 
  Double_t ptMinReco = ptBinEdgesTemplate[0];

  //=============================


  Int_t netaBins = etaBinEdgesTemplate.size() -1;
  Int_t nptBins = ptBinEdgesTemplate.size() -1;
  Int_t nBinsTemplate = nptBins * netaBins;
  
  TChain* chain = new TChain("tree");
  // INFO: the new friend trees at 13 TeV are inside a "tree_Friend_<sampleName>.root" file, and the tree's name is "Friends"
  // friend trees are still located in a directory called "friends" with respect to base trees
  TChain* friendChain = new TChain("Friends");      
  //TChain* friendChain = NULL;  // leave as NULL if you don't use friend trees
  TChain* SfFriendChain = NULL;
  //if (sampleDir.find("data") == string::npos && sampleDir.find("fake") == string::npos) SfFriendChain = new TChain("sf/t");  
  // leave as NULL if you don't use friend trees

  Bool_t noSumGenWeight = true;
  vector<Double_t> sumGenWeightVector;
  buildChain(chain, sumGenWeightVector, treedir, sample, friendChain, SfFriendChain, nameMatch, noSumGenWeight); 

  // change directory again, when building chain something was messed up
  if (FOLDER_IN_FILE) dirSample->cd();
  else                outputFile->cd();
  //  cout << "check" << endl;

  Bool_t isWsignal = false;
  if (sampleDir.find("wenujets") != string::npos or sampleDir.find("wmunujets") != string::npos) isWsignal = true;

  cout << "Setting TTreeReader and branches" << endl;
  TTreeReader reader (chain);

  TTreeReaderValue<Int_t> isData(reader,"isData");
  TTreeReaderValue<UInt_t> run   (reader,"run");
  TTreeReaderValue<UInt_t> lumi  (reader,"lumi");
  // //TTreeReaderValue<Int_t> nVert  (reader,"nVert");
  // //TTreeReaderValue<Float_t> rho  (reader,"rho");

  // // trigger
  TTreeReaderValue<Int_t>* HLT_SingleElectron = nullptr;
  TTreeReaderValue<Int_t>* HLT_SingleMuon_1 = nullptr;
  TTreeReaderValue<Int_t>* HLT_SingleMuon_2 = nullptr;

  if (isMuon) {
    HLT_SingleMuon_1 = new TTreeReaderValue<Int_t>(reader,"HLT_BIT_HLT_IsoMu24_v");
    HLT_SingleMuon_2 = new TTreeReaderValue<Int_t>(reader,"HLT_BIT_HLT_IsoTkMu24_v");
  } else {
    HLT_SingleElectron = new TTreeReaderValue<Int_t>(reader,"HLT_BIT_HLT_Ele27_WPTight_Gsf_v");
  }

  // // reco met
  // // TTreeReaderValue<Float_t> tkmet    (reader,"met_trkPt");
  // // TTreeReaderValue<Float_t> tkmet_phi(reader,"met_trkPhi");
  TTreeReaderValue<Float_t> pfmet    (reader,"met_pt");
  TTreeReaderValue<Float_t> pfmet_phi(reader,"met_phi");


  // // LepGood branch
  TTreeReaderValue<Int_t> nlep      (reader,"nLepGood");
  TTreeReaderArray<Int_t> lep_pdgId (reader,"LepGood_pdgId");
  TTreeReaderArray<Int_t> lep_charge(reader,"LepGood_charge");  
  TTreeReaderArray<Float_t> lep_calPt(reader,"LepGood_calPt");
  TTreeReaderArray<Float_t> lep_pt  (reader,"LepGood_pt");
  TTreeReaderArray<Float_t> lep_eta (reader,"LepGood_eta");
  TTreeReaderArray<Float_t> lep_phi (reader,"LepGood_phi");
  TTreeReaderArray<Float_t> lep_SF1 (reader,"LepGood_SF1");
  TTreeReaderArray<Float_t> lep_SF2 (reader,"LepGood_SF2");
  TTreeReaderArray<Float_t> lep_SF3 (reader,"LepGood_SF3");

  // // for electronID
  TTreeReaderArray<Int_t> lep_Id  (reader, isMuon ? "LepGood_mediumMuonId" : "LepGood_customId");  
  //TTreeReaderArray<Float_t> lep_Iso (reader, isMuon ? "LepGood_relIso04" : "LepGood_relIso04EA");  // actually, for ele the iso is already in LepGood_customId
  TTreeReaderArray<Int_t> *lep_hltId = nullptr;  
  TTreeReaderArray<Int_t> *lep_tightChargeFix = nullptr;  
  TTreeReaderArray<Int_t> *lep_mcMatchId = nullptr;
  if (not isMuon) {
    lep_hltId          = new TTreeReaderArray<Int_t>(reader,"LepGood_hltId");  
    lep_tightChargeFix = new TTreeReaderArray<Int_t>(reader,"LepGood_tightChargeFix");  
    lep_mcMatchId      = new TTreeReaderArray<Int_t>(reader,"LepGood_mcMatchId");
  }

  // // gen quantities
  TTreeReaderValue<Float_t> nTrueInt(reader,"nTrueInt");
  TTreeReaderValue<Float_t> xsec(reader,"xsec");
  TTreeReaderValue<Float_t> genWeight(reader,"genWeight");

  // // W MC specific branches
  // apparently it is UInt_t for electrons, and Int_t for muons
  //TTreeReaderValue<UInt_t> nGenLep(reader,usePreFSRvar ? "nGenLepPreFSR" : "nGenLepDressed");
  //TTreeReaderValue<Int_t> nGenLep(reader,usePreFSRvar ? "nGenLepPreFSR" : "nGenLepDressed");

  // value instead of array because they have only 1 entry per event
  // this suppresses the following error
  // Error in <TTreeReaderValueBase::CreateProxy()>: The branch nGenLepPreFSR contains data of type unsigned int. It cannot be accessed by a TTreeReaderValue<int>
  TTreeReaderValue<Float_t> GenLep_pdgId( reader, usePreFSRvar ? "GenLepPreFSR_pdgId"  : "GenLepDressed_pdgId"); 
  TTreeReaderValue<Float_t> GenLep_pt( reader, usePreFSRvar ? "GenLepPreFSR_pt"  : "GenLepDressed_pt"); 
  TTreeReaderValue<Float_t> GenLep_eta(reader, usePreFSRvar ? "GenLepPreFSR_eta" : "GenLepDressed_eta");
  TTreeReaderValue<Float_t> GenLep_phi(reader, usePreFSRvar ? "GenLepPreFSR_phi" : "GenLepDressed_phi");

  TTreeReaderValue<Float_t> GenPromptNu_pt(reader, "GenPromptNu_pt");
  TTreeReaderValue<Float_t> GenPromptNu_phi(reader, "GenPromptNu_phi");
  TTreeReaderValue<Float_t> GenLepBare_pt( reader, "GenLepBare_pt"); 

  TTreeReaderValue<Float_t> genw_charge( reader, usePreFSRvar ? "prefsrw_charge"  : "genw_charge");
  TTreeReaderValue<Float_t> genw_decayId(reader, usePreFSRvar ? "prefsrw_decayId" : "genw_decayId");
  TTreeReaderValue<Float_t> genw_pt(reader,      usePreFSRvar ? "prefsrw_pt"      : "genw_pt");
  // // syst weights for W
  TTreeReaderValue<Float_t> mass_80470(reader,"mass_80470"); // mw up, central is 80420
  TTreeReaderValue<Float_t> mass_80370(reader,"mass_80370"); // mw down, central is 80420
  TTreeReaderValue<Float_t> qcd_muRUp(reader,"qcd_muRUp");
  TTreeReaderValue<Float_t> qcd_muRDn(reader,"qcd_muRDn");
  TTreeReaderValue<Float_t> qcd_muFUp(reader,"qcd_muFUp");
  TTreeReaderValue<Float_t> qcd_muFDn(reader,"qcd_muFDn");
  TTreeReaderValue<Float_t> qcd_muRmuFUp(reader,"qcd_muRmuFUp");
  TTreeReaderValue<Float_t> qcd_muRmuFDn(reader,"qcd_muRmuFDn");
  TTreeReaderValue<Float_t> qcd_alphaSUp(reader,"qcd_alphaSUp");
  TTreeReaderValue<Float_t> qcd_alphaSDn(reader,"qcd_alphaSDn");
  // // pdf
  vector< TTreeReaderValue<Float_t>* > pdfwgt;
  for (Int_t i = 1; i <= nPDFweight; ++i) {
    pdfwgt.push_back( new TTreeReaderValue<Float_t>(reader,Form("hessWgt%d",i)) );
  }
  cout << "Branches have been set" << endl;

  //////////////////////////
  // NOTE
  // values of TTreeReaderValue objects are read with * before the variable name, as if they were pointers
  // Pointers to TTreeReaderValue are used with double * (one to access the pointer's content and the other for the TTreeReaderValue convention)
  // TTreeReaderArray variables are used as normal array variables (if they are pointers to TTreeReaderArray.then you still need a * to access the array)
  //////////////////////////

  cout << "Defining histograms" << endl;

  /////////////////////////////////////
  // dummy histogram to easily retrieve information on eta-pt bins
  TH2F* h2_etaPt = new TH2F("h2_etaPt","",nGenEtaBins,genEtaBinEdgesTemplate.data(),nGenPtBins,genPtBinEdgesTemplate.data());
  h2_etaPt->SetDirectory(0); // I don't want to save this histogram in the output file

  // 1 histogram for each bin of the template and for each charge, hence a double vector of histograms
  vector<string> charges = {"plus","minus"}; // 0 for positive, 1 for negative
  vector<string> chargeSigns = {"+","-"}; // 0 for positive, 1 for negative

  vector<TH2F*> h2_charge_eta_pt_SF1;
  vector<TH2F*> h2_charge_eta_pt_SF2;
  vector<TH2F*> h2_charge_eta_pt_SF3;

  vector<TH1F*> h1_charge_eta;
  vector<TH1F*> h1_charge_pt;
  vector<TH2F*> h2_charge_eta_pt_inclusive;

  vector<TH3F*> h3_charge_eta_pt_globalBin;

  // vector<TH3F*> h3_charge_eta_pt_globalBin_lepScaleUp;
  // vector<TH3F*> h3_charge_eta_pt_globalBin_lepScaleDn;
  vector<TH3F*> h3_charge_eta_pt_globalBin_lepEffUp;
  vector<TH3F*> h3_charge_eta_pt_globalBin_lepEffDn;

  Int_t nPtScaleRegions = isMuon ? 2 : 4; 
  // muons or electrons
  vector< vector<TH3F*> > h3_charge_eta_pt_globalBin_lepUncorrScaleUp(charges.size());  // N objects, 1 for each scale uncorrelated part (2 for muons, 4 for electrons)
  vector< vector<TH3F*> > h3_charge_eta_pt_globalBin_lepUncorrScaleDn(charges.size());

  vector< vector<TH3F*> > h3_charge_eta_pt_globalBin_pdf(charges.size());  // will have 60 replicas of pairs of TH3 (the pair is for charge + and -)
  vector<TH3F*> h3_charge_eta_pt_globalBin_alphaSUp;
  //vector<TH3F*> h3_charge_eta_pt_globalBin_wptSlopeUp;
  vector<TH3F*> h3_charge_eta_pt_globalBin_muRUp;
  vector<TH3F*> h3_charge_eta_pt_globalBin_muFUp;
  vector<TH3F*> h3_charge_eta_pt_globalBin_muRmuFUp;
  vector<TH3F*> h3_charge_eta_pt_globalBin_alphaSDn;
  //vector<TH3F*> h3_charge_eta_pt_globalBin_wptSlopeDn;
  vector<TH3F*> h3_charge_eta_pt_globalBin_muRDn;
  vector<TH3F*> h3_charge_eta_pt_globalBin_muFDn;
  vector<TH3F*> h3_charge_eta_pt_globalBin_muRmuFDn;

  vector<TH3F*> h3_charge_eta_pt_globalBin_mWUp;
  vector<TH3F*> h3_charge_eta_pt_globalBin_mWDn;

  vector<TH3F*> h3_charge_eta_pt_globalBin_fsr;

  // will use 10 bins of wptbins array for now
  vector< vector<TH3F*> > h3_charge_eta_pt_globalBin_muRUp_wpt(charges.size());
  vector< vector<TH3F*> > h3_charge_eta_pt_globalBin_muFUp_wpt(charges.size());
  vector< vector<TH3F*> > h3_charge_eta_pt_globalBin_muRmuFUp_wpt(charges.size());
  vector< vector<TH3F*> > h3_charge_eta_pt_globalBin_muRDn_wpt(charges.size());
  vector< vector<TH3F*> > h3_charge_eta_pt_globalBin_muFDn_wpt(charges.size());
  vector< vector<TH3F*> > h3_charge_eta_pt_globalBin_muRmuFDn_wpt(charges.size());


  // ErfParXXEffStatYY  XX = 0,1,2; YY from 1 to 48 or 50
  vector< vector<TH3F*> > h3_charge_eta_pt_globalBin_ErfPar0EffStat(charges.size());
  vector< vector<TH3F*> > h3_charge_eta_pt_globalBin_ErfPar1EffStat(charges.size());
  vector< vector<TH3F*> > h3_charge_eta_pt_globalBin_ErfPar2EffStat(charges.size());

  //vector< vector<TH3F*> > h3_charge_eta_pt_globalBin_BinUncEffStatUp(charges.size());
  //vector< vector<TH3F*> > h3_charge_eta_pt_globalBin_BinUncEffStatDn(charges.size());

  // better to do it as done below
  //Int_t maxEffStat = 2 * ((Int_t) (10 * (genEtaBinEdgesTemplate.back() + 0.0001)));
  //Double_t  etaminForEffStat = -1. * genEtaBinEdgesTemplate.back();

  // the reweighting template has 48 (50) bins for muons from -2.4 to 2.4 (for electrons from -2.5 to 2.5)
  // for signal, we will eventually keep only the nuisances that affect the specific gen_eta bin (in this macro we will define them all and filter them later)
  // if gen_eta for electrons stops at 2.5, then the only relevant EffStat nuisances will go from 2 to 49 rather than from 1 to 50
  // this works whatever gen_eta is used as the upper boundary, meaning that the EffStat 25 and 26 will always be associated to get_eta_bin = 0
  // in general, the gen_eta binning might be coarser than 0.1. This will be handled later by keeping more EffStat nuisances for that bin
  // the match should be based on the reco bins asociated to the gen bins, but generally the template binning (reco) will match the width of the gen_eta binning
  Int_t maxEffStat = isMuon ? 48 : 50;
  Double_t  etaminForEffStat = isMuon ? -2.4 : -2.5;
  cout << "maxEffStat | etaminForEffStat = " << maxEffStat << " | " << etaminForEffStat << endl;
  
  for (UInt_t ch = 0; ch < charges.size(); ++ch) {  
    // for (Int_t bin = 0; bin  < nBinsTemplate; ++bin) {
    //   Int_t ieta = 0;
    //   Int_t ipt = 0;
    //   h2_etaPt->GetBinXYZ(bin+1,ieta,ipt)
    vector<Double_t> globalBin_binning; // temporary vector to be used in TH3 constructor
    for (Int_t bin = 0; bin  <= nGenBinsTemplate; ++bin) 
      globalBin_binning.push_back(0.5+(Double_t)bin);

    h1_charge_pt.push_back(new TH1F(Form("h1_%s_pt",charges[ch].c_str()),"",40,25,65));
    h1_charge_eta.push_back(new TH1F(Form("h1_%s_eta",charges[ch].c_str()),"",50,-2.5,2.5));

    h2_charge_eta_pt_inclusive.push_back(new TH2F(Form("h2_%s_eta_pt_inclusive",charges[ch].c_str()),"",
						  50,-2.5,2.5,
						  40,25,65)	  
					 );

    h2_charge_eta_pt_SF1.push_back(new TH2F(Form("h2_%s_eta_pt_SF1",charges[ch].c_str()),"",
					    50,-2.5,2.5,
					    150,25,55)	  
				   );
    h2_charge_eta_pt_SF2.push_back(new TH2F(Form("h2_%s_eta_pt_SF2",charges[ch].c_str()),"",
					    50,-2.5,2.5,
					    150,25,55)	  
				   );
    h2_charge_eta_pt_SF3.push_back(new TH2F(Form("h2_%s_eta_pt_SF3",charges[ch].c_str()),"",
					    50,-2.5,2.5,
					    150,25,55)	  
				   );
    

    h3_charge_eta_pt_globalBin.push_back(new TH3F(Form("h3_%s_eta_pt_globalBin",charges[ch].c_str()),"",
						  netaBins,etaBinEdgesTemplate.data(),
						  nptBins,ptBinEdgesTemplate.data(),
						  nGenBinsTemplate,globalBin_binning.data())
					 );
    // h3_charge_eta_pt_globalBin_lepScaleUp.push_back(new TH3F(Form("h3_%s_eta_pt_globalBin_lepScaleUp",charges[ch].c_str()),"",
    // 							     netaBins,etaBinEdgesTemplate.data(),
    // 							     nptBins,ptBinEdgesTemplate.data(),
    // 							     nGenBinsTemplate,globalBin_binning.data())
    // 					 );
    // h3_charge_eta_pt_globalBin_lepScaleDn.push_back(new TH3F(Form("h3_%s_eta_pt_globalBin_lepScaleDn",charges[ch].c_str()),"",
    // 							     netaBins,etaBinEdgesTemplate.data(),
    // 							     nptBins,ptBinEdgesTemplate.data(),
    // 							     nGenBinsTemplate,globalBin_binning.data())
    // 					 );
    h3_charge_eta_pt_globalBin_lepEffUp.push_back(new TH3F(Form("h3_%s_eta_pt_globalBin_lepEffUp",charges[ch].c_str()),"",
							   netaBins,etaBinEdgesTemplate.data(),
							   nptBins,ptBinEdgesTemplate.data(),
							   nGenBinsTemplate,globalBin_binning.data())
						  );
    h3_charge_eta_pt_globalBin_lepEffDn.push_back(new TH3F(Form("h3_%s_eta_pt_globalBin_lepEffDn",charges[ch].c_str()),"",
							   netaBins,etaBinEdgesTemplate.data(),
							   nptBins,ptBinEdgesTemplate.data(),
							   nGenBinsTemplate,globalBin_binning.data())
					 );
    for (Int_t iptscale = 0; iptscale < nPtScaleRegions; iptscale++) {
      h3_charge_eta_pt_globalBin_lepUncorrScaleUp[ch].push_back(new TH3F(Form("h3_%s_eta_pt_globalBin_lepUncorrScale%dUp",charges[ch].c_str(),iptscale),"",
									 netaBins,etaBinEdgesTemplate.data(),
									 nptBins,ptBinEdgesTemplate.data(),
									 nGenBinsTemplate,globalBin_binning.data())
								);
      h3_charge_eta_pt_globalBin_lepUncorrScaleDn[ch].push_back(new TH3F(Form("h3_%s_eta_pt_globalBin_lepUncorrScale%dDn",charges[ch].c_str(),iptscale),"",
									 netaBins,etaBinEdgesTemplate.data(),
									 nptBins,ptBinEdgesTemplate.data(),
									 nGenBinsTemplate,globalBin_binning.data())
								);

    }

    h3_charge_eta_pt_globalBin_alphaSUp.push_back(new TH3F(Form("h3_%s_eta_pt_globalBin_alphaSUp",charges[ch].c_str()),"",
							   netaBins,etaBinEdgesTemplate.data(),
							   nptBins,ptBinEdgesTemplate.data(),
							   nGenBinsTemplate,globalBin_binning.data())
					 );
    h3_charge_eta_pt_globalBin_alphaSDn.push_back(new TH3F(Form("h3_%s_eta_pt_globalBin_alphaSDn",charges[ch].c_str()),"",
							   netaBins,etaBinEdgesTemplate.data(),
							   nptBins,ptBinEdgesTemplate.data(),
							   nGenBinsTemplate,globalBin_binning.data())
						  );
    h3_charge_eta_pt_globalBin_mWUp.push_back(new TH3F(Form("h3_%s_eta_pt_globalBin_mWUp",charges[ch].c_str()),"",
						       netaBins,etaBinEdgesTemplate.data(),
						       nptBins,ptBinEdgesTemplate.data(),
						       nGenBinsTemplate,globalBin_binning.data())
					      );
    h3_charge_eta_pt_globalBin_mWDn.push_back(new TH3F(Form("h3_%s_eta_pt_globalBin_mWDn",charges[ch].c_str()),"",
						       netaBins,etaBinEdgesTemplate.data(),
						       nptBins,ptBinEdgesTemplate.data(),
						       nGenBinsTemplate,globalBin_binning.data())
					      );
    h3_charge_eta_pt_globalBin_muRUp.push_back(new TH3F(Form("h3_%s_eta_pt_globalBin_muRUp",charges[ch].c_str()),"",
							   netaBins,etaBinEdgesTemplate.data(),
							   nptBins,ptBinEdgesTemplate.data(),
							   nGenBinsTemplate,globalBin_binning.data())
					       );
    h3_charge_eta_pt_globalBin_muRDn.push_back(new TH3F(Form("h3_%s_eta_pt_globalBin_muRDn",charges[ch].c_str()),"",
							   netaBins,etaBinEdgesTemplate.data(),
							   nptBins,ptBinEdgesTemplate.data(),
							   nGenBinsTemplate,globalBin_binning.data())
					       );
    h3_charge_eta_pt_globalBin_muFUp.push_back(new TH3F(Form("h3_%s_eta_pt_globalBin_muFUp",charges[ch].c_str()),"",
							   netaBins,etaBinEdgesTemplate.data(),
							   nptBins,ptBinEdgesTemplate.data(),
							   nGenBinsTemplate,globalBin_binning.data())
					       );
    h3_charge_eta_pt_globalBin_muFDn.push_back(new TH3F(Form("h3_%s_eta_pt_globalBin_muFDn",charges[ch].c_str()),"",
							   netaBins,etaBinEdgesTemplate.data(),
							   nptBins,ptBinEdgesTemplate.data(),
							   nGenBinsTemplate,globalBin_binning.data())
					       );
    h3_charge_eta_pt_globalBin_muRmuFUp.push_back(new TH3F(Form("h3_%s_eta_pt_globalBin_muRmuFUp",charges[ch].c_str()),"",
							   netaBins,etaBinEdgesTemplate.data(),
							   nptBins,ptBinEdgesTemplate.data(),
							   nGenBinsTemplate,globalBin_binning.data())
						  );
    h3_charge_eta_pt_globalBin_muRmuFDn.push_back(new TH3F(Form("h3_%s_eta_pt_globalBin_muRmuFDn",charges[ch].c_str()),"",
							   netaBins,etaBinEdgesTemplate.data(),
							   nptBins,ptBinEdgesTemplate.data(),
							   nGenBinsTemplate,globalBin_binning.data())
						  );
    // h3_charge_eta_pt_globalBin_wptSlopeUp.push_back(new TH3F(Form("h3_%s_eta_pt_globalBin_wptSlopeUp",charges[ch].c_str()),"",
    // 							   netaBins,etaBinEdgesTemplate.data(),
    // 							   nptBins,ptBinEdgesTemplate.data(),
    // 							   nGenBinsTemplate,globalBin_binning.data())
    // 					 );
    // h3_charge_eta_pt_globalBin_wptSlopeDn.push_back(new TH3F(Form("h3_%s_eta_pt_globalBin_wptSlopeDn",charges[ch].c_str()),"",
    // 							   netaBins,etaBinEdgesTemplate.data(),
    // 							   nptBins,ptBinEdgesTemplate.data(),
    // 							   nGenBinsTemplate,globalBin_binning.data())
    // 					 );

    h3_charge_eta_pt_globalBin_fsr.push_back(new TH3F(Form("h3_%s_eta_pt_globalBin_fsr",charges[ch].c_str()),"",
						      netaBins,etaBinEdgesTemplate.data(),
						      nptBins,ptBinEdgesTemplate.data(),
						      nGenBinsTemplate,globalBin_binning.data())
					     );

    for (Int_t ipdf = 1; ipdf <= nPDFweight; ++ipdf) {
      h3_charge_eta_pt_globalBin_pdf[ch].push_back(new TH3F(Form("h3_%s_eta_pt_globalBin_pdf%d",charges[ch].c_str(),ipdf),"",
							    netaBins,etaBinEdgesTemplate.data(),
							    nptBins,ptBinEdgesTemplate.data(),
							    nGenBinsTemplate,globalBin_binning.data())
						   );
    }

    for (Int_t iwpt = 1; iwpt <= 10; ++iwpt) {     
      h3_charge_eta_pt_globalBin_muRmuFUp_wpt[ch].push_back(new TH3F(Form("h3_%s_eta_pt_globalBin_muRmuF%dUp",charges[ch].c_str(),iwpt),"",
								     netaBins,etaBinEdgesTemplate.data(),
								     nptBins,ptBinEdgesTemplate.data(),
								     nGenBinsTemplate,globalBin_binning.data())
							    );
      h3_charge_eta_pt_globalBin_muRUp_wpt[ch].push_back(new TH3F(Form("h3_%s_eta_pt_globalBin_muR%dUp",charges[ch].c_str(),iwpt),"",
								  netaBins,etaBinEdgesTemplate.data(),
								  nptBins,ptBinEdgesTemplate.data(),
								  nGenBinsTemplate,globalBin_binning.data())
							 );
      h3_charge_eta_pt_globalBin_muFUp_wpt[ch].push_back(new TH3F(Form("h3_%s_eta_pt_globalBin_muF%dUp",charges[ch].c_str(),iwpt),"",
								  netaBins,etaBinEdgesTemplate.data(),
								  nptBins,ptBinEdgesTemplate.data(),
								  nGenBinsTemplate,globalBin_binning.data())
							 );
      h3_charge_eta_pt_globalBin_muRmuFDn_wpt[ch].push_back(new TH3F(Form("h3_%s_eta_pt_globalBin_muRmuF%dDn",charges[ch].c_str(),iwpt),"",
								     netaBins,etaBinEdgesTemplate.data(),
								     nptBins,ptBinEdgesTemplate.data(),
								     nGenBinsTemplate,globalBin_binning.data())
							    );
      h3_charge_eta_pt_globalBin_muRDn_wpt[ch].push_back(new TH3F(Form("h3_%s_eta_pt_globalBin_muR%dDn",charges[ch].c_str(),iwpt),"",
								  netaBins,etaBinEdgesTemplate.data(),
								  nptBins,ptBinEdgesTemplate.data(),
								  nGenBinsTemplate,globalBin_binning.data())
							 );
      h3_charge_eta_pt_globalBin_muFDn_wpt[ch].push_back(new TH3F(Form("h3_%s_eta_pt_globalBin_muF%dDn",charges[ch].c_str(),iwpt),"",
								  netaBins,etaBinEdgesTemplate.data(),
								  nptBins,ptBinEdgesTemplate.data(),
								  nGenBinsTemplate,globalBin_binning.data())
							 );
    }    


    // this corresponds to 3*48 histograms for each signal bin, because we have 3 independent variations for each eta bin
    // actually, this is redundant, because the templates for differential cross section are sparse, and almost only 1 eta bin is populated
    // therefore, we could just define 3 new TH3 and fill them with the erf weight corresponding to that gen-eta
    // or, we could just skip the useless variation in a later step, like in makeTH1FromTH3.py
    // if we define approx 150 histograms we are a little more versatile
    for (Int_t ieff = 0; ieff < maxEffStat; ++ieff) {     
      h3_charge_eta_pt_globalBin_ErfPar0EffStat[ch].push_back(new TH3F(Form("h3_%s_eta_pt_globalBin_ErfPar0EffStat%d",charges[ch].c_str(),ieff+1),"",
								       netaBins,etaBinEdgesTemplate.data(),
								       nptBins,ptBinEdgesTemplate.data(),
								       nGenBinsTemplate,globalBin_binning.data())
							      );
      h3_charge_eta_pt_globalBin_ErfPar1EffStat[ch].push_back(new TH3F(Form("h3_%s_eta_pt_globalBin_ErfPar1EffStat%d",charges[ch].c_str(),ieff+1),"",
								       netaBins,etaBinEdgesTemplate.data(),
								       nptBins,ptBinEdgesTemplate.data(),
								       nGenBinsTemplate,globalBin_binning.data())
							      );
      h3_charge_eta_pt_globalBin_ErfPar2EffStat[ch].push_back(new TH3F(Form("h3_%s_eta_pt_globalBin_ErfPar2EffStat%d",charges[ch].c_str(),ieff+1),"",
								       netaBins,etaBinEdgesTemplate.data(),
								       nptBins,ptBinEdgesTemplate.data(),
								       nGenBinsTemplate,globalBin_binning.data())
							      );
    }


  }
  cout << "Done histograms" << endl;

  ////////////////////
  ////////////////////

  // start event loop                                                                                                                     
                    
  long int nTotal = chain->GetEntries();
  long int nEvents = 0;
  long int nEventsInSample = 0; // events processed for each sample

  Int_t EB0orEE1 = -1;

  Bool_t negativeLeptonHasPassedSelection = false;
  Bool_t positiveLeptonHasPassedSelection = false;

  Double_t wgt = 1.0;
  Double_t wgt_noPt = 1.0; // will store wgt before applying selections dependent on pt, so to use that for the ptscale systematics
  Double_t wgtMuTrigSFup = 1.0;
  Double_t wgtMuTrigSFdn = 1.0;
  // Double_t wgt_ptscaleUp = 1.0;
  // Double_t wgt_ptscaleDn = 1.0;
  Double_t lep1pt = 0.0;
  Double_t ptlow = 0.0;
  Double_t pthigh = 0.0;
  Double_t etalow = 0.0;
  Double_t etahigh = 0.0;
  Double_t lepEffWgtUp = 0.0;
  Double_t lepEffWgtDn = 0.0;
  // Double_t ptLepFullUp = 0.0;
  // Double_t ptLepFullDn = 0.0;
  Double_t lep1calPt = 0.0;

  Bool_t nominalPt_passSelection = false;
  Bool_t scaleUpPt_passSelection = false;
  Bool_t scaleDnPt_passSelection = false;

  vector<Bool_t> scaleUncorrUpPt_passSelection;
  vector<Bool_t> scaleUncorrDnPt_passSelection;
  vector<Double_t> ptUncorrUp;
  vector<Double_t> ptUncorrDn;

  for (Int_t ipt = 0; ipt < nPtScaleRegions; ipt++) {
    scaleUncorrUpPt_passSelection.push_back(false);
    scaleUncorrDnPt_passSelection.push_back(false);
    ptUncorrUp.push_back(0.0);
    ptUncorrDn.push_back(0.0);
  }

  ////////////////////////////////////////////
  // to get correct weight depending on sample in chain
  string currentFile = "";
  Int_t ifile = 0;
  ////////////////////

  Double_t sfTriggerMu = 0.0;
  Double_t sfRecoToSelectionMu = 0.0;
  Double_t sfPrefireMu = 0.0;

  Double_t sumWgt = 9.56169443709e+13;
  // if (treedir == "/eos/cms/store/cmst3/group/wmass/w-helicity-13TeV/trees/SKIMS_muons_latest/")
  //   sumWgt = 9.47291822594e+13;
  // else if (treedir == "/eos/cms/store/cmst3/group/wmass/mciprian/TREES_1LEP_80X_V3_SIGSKIM_WMUNU_FULLSEL_14Feb2019/")
  //   sumWgt = 9.47291822594e+13;
  // else if (treedir == "/u1/mciprian/trees/muon/TREES_1LEP_80X_V3_SIGSKIM_WMUNU_FULLSEL_14Feb2019/")
  //   sumWgt = 9.47291822594e+13;

  if (treedir.find("TREES_1LEP_80X_V3_SIGSKIM_WMUNU_FULLSEL_14Feb2019") != string::npos || 
      treedir.find("SKIMS_muons_latest") != string::npos) {
    sumWgt = 9.47291822594e+13;
    cout << "Warning: setting sumWgt = " << sumWgt << endl;
  }
  
  Double_t intLumiPb = 1000.0 * intLumi;
  Double_t intLumiPbXsecZ = intLumiPb * 2008.4 * 3.; // for Z the xsec in the ntuples is no more valid, it changed
  //Double_t wjets_NLO_wgt_partial = intLumiPb * (3. * 20508.9) / sumWgt;  // WJetsToLNu_, just to speed up, for electrons and muons (same number)
  // use native MC@NLO xsec, not fewz3.1
  Double_t wjets_NLO_wgt_partial = intLumiPb * 60400.0 / sumWgt;  // WJetsToLNu_, just to speed up, for electro

  //Double_t wjets_NLO_wgt_partial = intLumiPb * (3. * 20508.9) / 7.97011396982e+11;  // test with WJetsToLNu_ext2v5_part1 in pccmsrm
  
  // use 9.56169443709e+13 for ntuples in /eos/cms/store/cmst3/group/wmass/mciprian/TREES_1LEP_80X_V3_SIGSKIM_WMUNU_FULLSEL_NOMT_V2/
  // use 9.56169443709e+13 for ntuples in /eos/cms/store/cmst3/group/wmass/w-helicity-13TeV/ntuplesRecoil/TREES_SIGNAL_1l_recoil_fullTrees/
  // use 9.56169443709e+13 for ntuples in /u2/mciprian/TREES_13TeV/muon/signalSkim/
  // use 9.47291822594e+13 with /afs/cern.ch/work/m/mdunser/public/wmassTrees/SKIMS_muons_latest/   // requires reading tree.root.url, might not work, see line below
  // use 9.47291822594e+13 with /eos/cms/store/cmst3/group/wmass/w-helicity-13TeV/trees/SKIMS_muons_latest/
  // 3.54324749853e+13; // obsolete without ext2v5                      

  cout << "Starting loop" << endl;

  //cout << "check -2" << endl;
  while (reader.Next()) {
  
    //cout << "check -1" << endl;
    cout.flush();
    if (nEvents % CHECK_EVERY_N == 0) cout << "\r" << "Analyzing events " << ((Double_t) nEvents)/nTotal*100.0 << " % ";
    //cout << "entry : " << nEvents << endl;
    nEvents++;
    nEventsInSample++;

    if (dynamic_cast<TChain*>(reader.GetTree())->GetFile()->GetName() != currentFile and currentFile != "") { 
      currentFile = dynamic_cast<TChain*>(reader.GetTree())->GetFile()->GetName();                   
      ifile ++;                                                                                      
      nEventsInSample = 1; // reset nEvents when sub sample is changed (useful with N_MAX_ENTRIES_PER_SAMPLE for debugging)
    } else if (dynamic_cast<TChain*>(reader.GetTree())->GetFile()->GetName() != currentFile and currentFile == "") {
      currentFile = dynamic_cast<TChain*>(reader.GetTree())->GetFile()->GetName();                         
    }      
    if (N_MAX_ENTRIES_PER_SAMPLE > 0 && nEventsInSample > N_MAX_ENTRIES_PER_SAMPLE) continue; // this does not skip the tree, but still loop on all events doing nothing

    //if (nEvents > 100) break;  // temporary line for tests
    //checkPoint(currentFile,200,nEvents);

    negativeLeptonHasPassedSelection = false;
    positiveLeptonHasPassedSelection = false;
  
    Double_t absLep1eta = fabs(lep_eta[0]);  
    Double_t absGenEta = fabs(*GenLep_eta);  

    //cout << "check 0" << endl;

    if (isMuon) {

      lep1calPt = ptMuFull(lep_calPt[0],lep_eta[0]);

      // selection electrons
      if (**HLT_SingleMuon_1 < 1 and **HLT_SingleMuon_2 < 1) continue;
      if (*nlep != 1) continue;
      if (lep_Id[0] < 1) continue;
      if (fabs(lep_pdgId[0]) != 13) continue;
      if (absLep1eta >= etaCutReco) continue;
      if (not isFloatEqual(fabs(*genw_decayId), 14.0, 0.001)) continue;
      // genw_charge is saved as float, while lep_charge is integer, protect against floating point precision in the comparison
      if (fabs(*genw_charge - lep_charge[0]) > 0.01) continue;   
      // do not cut here on pt, otherwise the muscale syst is biased, I have to cut on the consistent pt
      //if (lep1calPt < 26 or lep1calPt > ptCutReco) continue;
      //if (mt_2(*pfmet,*pfmet_phi,lep1calPt,lep_phi[0]) < 40) continue;

      // PU reweigthing, trigger scale factors, lepton efficiency scale factors
      // done like this to speed it up
      bool useBinnedSF = false;  // true for tests, otherwise false 
      sfTriggerMu = _get_muonSF_selectionToTrigger(lep_pdgId[0], lep1calPt, lep_eta[0], lep_charge[0], 0.0, false, useBinnedSF);
      wgtMuTrigSFup = _get_muonSF_selectionToTrigger(lep_pdgId[0], lep1calPt, lep_eta[0], lep_charge[0], 1.0, true, useBinnedSF);
      wgtMuTrigSFdn = _get_muonSF_selectionToTrigger(lep_pdgId[0], lep1calPt, lep_eta[0], lep_charge[0], -1.0, true, useBinnedSF);
      // the reco2Selection SF is eta-smoothed
      //sfRecoToSelectionMu = lep_SF2[0]; 
      sfRecoToSelectionMu = _get_muonSF_recoToSelection(lep_pdgId[0], lep1calPt, lep_eta[0], useBinnedSF); // this is also eta-smoothed 
      sfPrefireMu = prefireJetsWeight(lep_eta[0]);
      wgt = sfTriggerMu * sfRecoToSelectionMu * sfPrefireMu; 
      // for muons, get fast weight for efficiency Up/Down
      Double_t syst = 0.0;
      if (absLep1eta<1)          syst = 0.002;
      else if (absLep1eta<1.5)   syst = 0.004;
      else                       syst = 0.014;
      lepEffWgtUp = (wgt + syst)/wgt; // because in the following it will multiply the nominal wgt, which includes the current definition of wgt
      lepEffWgtDn = (wgt - syst)/wgt;
      
      // ptLepFullUp = ptMuFullUp(lep_calPt[0],lep_eta[0]);
      // ptLepFullDn = ptMuFullDn(lep_calPt[0],lep_eta[0]);
      //cout << "pt, up, down: " << lep_calPt[0] << ", " << ptLepFullUp << ", " << ptLepFullDn << endl;

      // if (ptLepFullUp > 26 and ptLepFullUp < ptCutReco and mt_2(*pfmet,*pfmet_phi,ptLepFullUp,lep_phi[0]) > 40)  
      // 	scaleUpPt_passSelection = true;
      // else
      // 	scaleUpPt_passSelection = false;

      // if (ptLepFullDn > 26 and ptLepFullDn < ptCutReco and mt_2(*pfmet,*pfmet_phi,ptLepFullDn,lep_phi[0]) > 40)  
      // 	scaleDnPt_passSelection = true;
      // else
      // 	scaleDnPt_passSelection = false;

    } else {

      lep1calPt = ptElFull(lep_calPt[0],lep_eta[0]);   

      // selection electrons
      if ((*lep_hltId)[0] < 1) continue;
      if (**HLT_SingleElectron != 1) continue;
      if (*nlep != 1) continue;
      if (lep_Id[0] != 1 or (*lep_tightChargeFix)[0] != 2) continue;
      if (absLep1eta > 1.4442 and absLep1eta < 1.566) continue;
      if (absLep1eta > etaCutReco) continue;
      if (fabs(lep_pdgId[0]) != 11) continue;
      if (not isFloatEqual(*genw_decayId, 12.0, 0.001)) continue;
      if ((*lep_mcMatchId)[0] * lep_charge[0] == -24) continue;
      // if (lep1calPt < 30 or lep1calPt > ptCutReco) continue;
      // if (mt_2(*pfmet,*pfmet_phi,lep1calPt,lep_phi[0]) < 40) continue;

      // PU reweigthing, trigger scale factors, lepton efficiency scale factors
      // done like this to speed it up
      //wgt = lepSF(lep_pdgId[0],lep_pt[0],lep_eta[0],lep_SF1[0],lep_SF2[0],lep_SF3[0]); 
      wgt = _get_electronSF_TriggerAndID(lep_pdgId[0],lep_calPt[0],lep_eta[0]) * lep_SF2[0] * eleSF_L1Eff(lep_pt[0],lep_eta[0]);     
      lepEffWgtUp = lepSFRelUp(lep_pdgId[0],lep_pt[0],lep_eta[0],lep_SF1[0],lep_SF2[0],lep_SF3[0]);
      lepEffWgtDn = lepSFRelDn(lep_pdgId[0],lep_pt[0],lep_eta[0],lep_SF1[0],lep_SF2[0],lep_SF3[0]);
      // ptLepFullUp = ptElFullUp(lep_calPt[0],lep_eta[0]);
      // ptLepFullDn = ptElFullDn(lep_calPt[0],lep_eta[0]);

      // if (ptLepFullUp > 30 and ptLepFullUp < ptCutReco and mt_2(*pfmet,*pfmet_phi,ptLepFullUp,lep_phi[0]) > 40)  
      // 	scaleUpPt_passSelection = true;
      // else
      // 	scaleUpPt_passSelection = false;

      // if (ptLepFullDn > 30 and ptLepFullDn < ptCutReco and mt_2(*pfmet,*pfmet_phi,ptLepFullDn,lep_phi[0]) > 40)  
      // 	scaleDnPt_passSelection = true;
      // else
      // 	scaleDnPt_passSelection = false;


    }

    if (lep1calPt > ptMinReco and lep1calPt < ptCutReco and mt_2(*pfmet,*pfmet_phi,lep1calPt,lep_phi[0]) > 40)  
      nominalPt_passSelection = true;
    else
      nominalPt_passSelection = false;

    for (Int_t ipt = 0; ipt < nPtScaleRegions; ipt++) {

      ptUncorrUp[ipt] = ptScaleUncorr(lep_calPt[0],lep_eta[0], lep_pdgId[0], ipt, true);
      ptUncorrDn[ipt] = ptScaleUncorr(lep_calPt[0],lep_eta[0], lep_pdgId[0], ipt, false);

      if (ptUncorrUp[ipt] > ptMinReco and ptUncorrUp[ipt] < ptCutReco and mt_2(*pfmet,*pfmet_phi,ptUncorrUp[ipt],lep_phi[0]) > 40)
	scaleUncorrUpPt_passSelection[ipt] = true;
      else
	scaleUncorrUpPt_passSelection[ipt] = false;

      if (ptUncorrDn[ipt] > ptMinReco and ptUncorrDn[ipt] < ptCutReco and mt_2(*pfmet,*pfmet_phi,ptUncorrDn[ipt],lep_phi[0]) > 40)
	scaleUncorrDnPt_passSelection[ipt] = true;
      else
	scaleUncorrDnPt_passSelection[ipt] = false;
	  
    }

    // add weight for zpt (we use it on W as well)
    // if we used the loop to run on Z, we should use GenLepDressed_pt[1],GenLepDressed_phi[1] instead of GenPromptNu_pt[0],GenPromptNu_phi[0]
    // note that we should always use Dressed lepton (I assume GenLep points to that and not to PreFsr variables, although the difference is generally small)
    wgt *= dyptWeight(pt_2(*GenLep_pt,*GenLep_phi,*GenPromptNu_pt,*GenPromptNu_phi),0, true);
    // try to put common factor together
    wgt *= (wjets_NLO_wgt_partial * *genWeight * puw2016_nTrueInt_36fb(*nTrueInt));

    wgt_noPt = wgt;
    // wgt_ptscaleUp = wgt_noPt * scaleUpPt_passSelection;
    // wgt_ptscaleDn = wgt_noPt * scaleDnPt_passSelection;
    wgt           = wgt_noPt * nominalPt_passSelection;

    // Get charge index: 0 for positive, 1 for negative
    Int_t chargeIndex = 0; 
    Int_t chargeSign = 0; 
    if (lep_pdgId[0] > 0) {
      negativeLeptonHasPassedSelection = true;
      chargeIndex = 1;
      chargeSign = -1;
    } else {
      positiveLeptonHasPassedSelection = true;
      chargeIndex = 0;
      chargeSign = 1;
    }

    //cout << "check 1" << endl;
    
    // Now start filling histograms
    h1_charge_pt[chargeIndex]->Fill(lep1calPt, wgt);
    h1_charge_eta[chargeIndex]->Fill(lep_eta[0], wgt);
    h2_charge_eta_pt_inclusive[chargeIndex]->Fill(lep_eta[0],lep1calPt, wgt);

    if (nominalPt_passSelection) {
      Int_t binx = h2_charge_eta_pt_SF1[chargeIndex]->GetXaxis()->FindFixBin(lep_eta[0]);
      Int_t biny = h2_charge_eta_pt_SF1[chargeIndex]->GetYaxis()->FindFixBin(lep1calPt);
      h2_charge_eta_pt_SF1[chargeIndex]->SetBinContent(binx, biny, isMuon ? sfTriggerMu : lep_SF1[0]);
      h2_charge_eta_pt_SF2[chargeIndex]->SetBinContent(binx, biny, isMuon ? sfRecoToSelectionMu : lep_SF2[0]);
      h2_charge_eta_pt_SF3[chargeIndex]->SetBinContent(binx, biny, lep_SF3[0]);
    }

    // the following is not what we would expect, TH2 is internally built as a TH1 including all the underflow and overflow bins
    // suppose you have a TH2 with 3x2 bins, this is a 2D grip made as 
    //  0  0  0  0  0
    //  0  1  1  1  0
    //  0  1  1  1  0
    //  0  0  0  0  0
    // where 0 is either underflow or overflow, and the number returned by h2_etaPt->FindFixBin(*GenLep_eta, *GenLep_pt) starts from 0 to (NX+2)*(NY+2)-1
    // where NY = 2, NX = 3, therefore it is better to get the global bin by hand using the TH1 methods for X and Y axis
    // Int_t globalBinEtaPt = h2_etaPt->FindFixBin(*GenLep_eta, *GenLep_pt);  // wrong
    // Instead, we want a global bin number going from 1 to 6, where 0 and 7 will be underflow and overflow
    Int_t genEtaBin = h2_etaPt->GetXaxis()->FindFixBin(absGenEta);  // goes from 0 to nGenEtaBins + 1 (underflow and overflow)
    Int_t genPtBin =  h2_etaPt->GetYaxis()->FindFixBin(*GenLep_pt); // goes from 0 to nGenPtBins + 1 (underflow and overflow)
    Int_t globalBinEtaPt = 0;
    if (genEtaBin && genPtBin && (genEtaBin <= nGenEtaBins) && (genPtBin <= nGenPtBins)) 
      globalBinEtaPt = genEtaBin + nGenEtaBins * (genPtBin - 1);  // goes from 1 to nGenEtaBins*nGenPtBins
    // the underflow (bin 0) can be used for the outliers, the overflow for the moment is empty by definition
    // One could devise different ways of separating the under/over-flow bins    

    // fill with reco quantities the bin corresponding to the gen level quantities (it is equivalent to cutting on gen level variables)
    h3_charge_eta_pt_globalBin[chargeIndex]->Fill(lep_eta[0],lep1calPt, globalBinEtaPt, wgt);
    // h3_charge_eta_pt_globalBin_lepScaleUp[chargeIndex]->Fill(lep_eta[0],ptLepFullUp, globalBinEtaPt, wgt_ptscaleUp);
    // h3_charge_eta_pt_globalBin_lepScaleDn[chargeIndex]->Fill(lep_eta[0],ptLepFullDn, globalBinEtaPt, wgt_ptscaleDn);
    for (Int_t ipt = 0; ipt < nPtScaleRegions; ipt++) {
      h3_charge_eta_pt_globalBin_lepUncorrScaleUp[chargeIndex][ipt]->Fill(lep_eta[0],ptUncorrUp[ipt], globalBinEtaPt, wgt_noPt * scaleUncorrUpPt_passSelection[ipt]);
      h3_charge_eta_pt_globalBin_lepUncorrScaleDn[chargeIndex][ipt]->Fill(lep_eta[0],ptUncorrDn[ipt], globalBinEtaPt, wgt_noPt * scaleUncorrDnPt_passSelection[ipt]);
    }
    h3_charge_eta_pt_globalBin_lepEffUp[chargeIndex]->Fill(lep_eta[0],lep1calPt, globalBinEtaPt, lepEffWgtUp * wgt);
    h3_charge_eta_pt_globalBin_lepEffDn[chargeIndex]->Fill(lep_eta[0],lep1calPt, globalBinEtaPt, lepEffWgtDn * wgt);
    // wpt syst
    //h3_charge_eta_pt_globalBin_wptSlopeUp[chargeIndex]->Fill(lep_eta[0],lep1calPt, globalBinEtaPt, wpt_slope_weight(*genw_pt,0.95,0.005) * wgt);
    //h3_charge_eta_pt_globalBin_wptSlopeDn[chargeIndex]->Fill(lep_eta[0],lep1calPt, globalBinEtaPt, wpt_slope_weight(*genw_pt,1.05,-0005) * wgt);
    // qcd syst
    h3_charge_eta_pt_globalBin_muRUp[chargeIndex]->Fill(lep_eta[0],lep1calPt, globalBinEtaPt, *qcd_muRUp * wgt);
    h3_charge_eta_pt_globalBin_muRDn[chargeIndex]->Fill(lep_eta[0],lep1calPt, globalBinEtaPt, *qcd_muRDn * wgt);
    h3_charge_eta_pt_globalBin_muFUp[chargeIndex]->Fill(lep_eta[0],lep1calPt, globalBinEtaPt, *qcd_muFUp * wgt);
    h3_charge_eta_pt_globalBin_muFDn[chargeIndex]->Fill(lep_eta[0],lep1calPt, globalBinEtaPt, *qcd_muFDn * wgt);
    h3_charge_eta_pt_globalBin_muRmuFUp[chargeIndex]->Fill(lep_eta[0],lep1calPt, globalBinEtaPt, *qcd_muRmuFUp * wgt);
    h3_charge_eta_pt_globalBin_muRmuFDn[chargeIndex]->Fill(lep_eta[0],lep1calPt, globalBinEtaPt, *qcd_muRmuFDn * wgt);
    h3_charge_eta_pt_globalBin_alphaSUp[chargeIndex]->Fill(lep_eta[0],lep1calPt, globalBinEtaPt, *qcd_alphaSUp * wgt);
    h3_charge_eta_pt_globalBin_alphaSDn[chargeIndex]->Fill(lep_eta[0],lep1calPt, globalBinEtaPt, *qcd_alphaSDn * wgt);
    h3_charge_eta_pt_globalBin_mWUp[chargeIndex]->Fill(lep_eta[0],lep1calPt, globalBinEtaPt, *mass_80470 * wgt);
    h3_charge_eta_pt_globalBin_mWDn[chargeIndex]->Fill(lep_eta[0],lep1calPt, globalBinEtaPt, *mass_80370 * wgt);
    //Double_t fsrweight = fsrPhotosWeight(*GenLep_pdgId,*GenLep_eta,*GenLep_pt,*GenLepBare_pt);
    Double_t fsrweight = fsrPhotosWeightSimple(*GenLep_pdgId,*GenLep_pt,*GenLepBare_pt, true, *GenLep_eta);
    h3_charge_eta_pt_globalBin_fsr[chargeIndex]->Fill(lep_eta[0],lep1calPt, globalBinEtaPt, fsrweight * wgt);
    
    // pdf syst
    for (Int_t ipdf = 0; ipdf < nPDFweight; ++ipdf) {
      h3_charge_eta_pt_globalBin_pdf[chargeIndex][ipdf]->Fill(lep_eta[0],lep1calPt, globalBinEtaPt, **(pdfwgt[ipdf]) * wgt);
    }

    // cout << "check 2" << endl;

    for (Int_t iwpt = 0; iwpt < 10; ++iwpt) {
      // hardcoded, get 1 every 2 bins instead of defining a different wptbins
      ptlow = wptbins[2*iwpt];
      pthigh = wptbins[2*(iwpt+1)];
      if (*genw_pt > ptlow and *genw_pt < pthigh) {
    	h3_charge_eta_pt_globalBin_muRUp_wpt[chargeIndex][iwpt]   ->Fill(lep_eta[0],lep1calPt, globalBinEtaPt, *qcd_muRUp * wgt);
    	h3_charge_eta_pt_globalBin_muFUp_wpt[chargeIndex][iwpt]   ->Fill(lep_eta[0],lep1calPt, globalBinEtaPt, *qcd_muFUp * wgt);
    	h3_charge_eta_pt_globalBin_muRmuFUp_wpt[chargeIndex][iwpt]->Fill(lep_eta[0],lep1calPt, globalBinEtaPt, *qcd_muRmuFUp * wgt);
    	h3_charge_eta_pt_globalBin_muRDn_wpt[chargeIndex][iwpt]   ->Fill(lep_eta[0],lep1calPt, globalBinEtaPt, *qcd_muRDn * wgt);
    	h3_charge_eta_pt_globalBin_muFDn_wpt[chargeIndex][iwpt]   ->Fill(lep_eta[0],lep1calPt, globalBinEtaPt, *qcd_muFDn * wgt);
    	h3_charge_eta_pt_globalBin_muRmuFDn_wpt[chargeIndex][iwpt]->Fill(lep_eta[0],lep1calPt, globalBinEtaPt, *qcd_muRmuFDn * wgt);     
      } else {
    	h3_charge_eta_pt_globalBin_muRUp_wpt[chargeIndex][iwpt]   ->Fill(lep_eta[0],lep1calPt, globalBinEtaPt, wgt);
    	h3_charge_eta_pt_globalBin_muFUp_wpt[chargeIndex][iwpt]   ->Fill(lep_eta[0],lep1calPt, globalBinEtaPt, wgt);
    	h3_charge_eta_pt_globalBin_muRmuFUp_wpt[chargeIndex][iwpt]->Fill(lep_eta[0],lep1calPt, globalBinEtaPt, wgt);
    	h3_charge_eta_pt_globalBin_muRDn_wpt[chargeIndex][iwpt]   ->Fill(lep_eta[0],lep1calPt, globalBinEtaPt, wgt);
    	h3_charge_eta_pt_globalBin_muFDn_wpt[chargeIndex][iwpt]   ->Fill(lep_eta[0],lep1calPt, globalBinEtaPt, wgt);
    	h3_charge_eta_pt_globalBin_muRmuFDn_wpt[chargeIndex][iwpt]->Fill(lep_eta[0],lep1calPt, globalBinEtaPt, wgt);     
      }
    }

    for (Int_t ieff = 0; ieff < maxEffStat; ++ieff) {     
      etalow = etaminForEffStat + 0.1 * ieff;
      etahigh = etalow + 0.1;
      h3_charge_eta_pt_globalBin_ErfPar0EffStat[chargeIndex][ieff]->Fill(lep_eta[0],lep1calPt, globalBinEtaPt, 
      									 effSystEtaBins(0,lep_pdgId[0],lep_eta[0],lep1calPt,etalow,etahigh,chargeSign) * wgt);
      h3_charge_eta_pt_globalBin_ErfPar1EffStat[chargeIndex][ieff]->Fill(lep_eta[0],lep1calPt, globalBinEtaPt, 
      									 effSystEtaBins(1,lep_pdgId[0],lep_eta[0],lep1calPt,etalow,etahigh,chargeSign) * wgt);
      h3_charge_eta_pt_globalBin_ErfPar2EffStat[chargeIndex][ieff]->Fill(lep_eta[0],lep1calPt, globalBinEtaPt, 
      									 effSystEtaBins(2,lep_pdgId[0],lep_eta[0],lep1calPt,etalow,etahigh,chargeSign) * wgt);
      // test for muons only      
      // if (isMuon) {
      // 	Double_t binunc_scaling_up = 1.0; 
      // 	Double_t binunc_scaling_dn = 1.0; 
      // 	if (lep_eta[0] > etalow && lep_eta[0] < etahigh) {
      // 	  binunc_scaling_up = wgtMuTrigSFup / sfTriggerMu; 
      // 	  binunc_scaling_dn = wgtMuTrigSFdn / sfTriggerMu; 	  	  
      // 	}

      // 	h3_charge_eta_pt_globalBin_BinUncEffStatUp[chargeIndex][ieff]->Fill(lep_eta[0],lep1calPt, globalBinEtaPt, 
      // 									    wgt * binunc_scaling_up);
      // 	h3_charge_eta_pt_globalBin_BinUncEffStatDn[chargeIndex][ieff]->Fill(lep_eta[0],lep1calPt, globalBinEtaPt, 
      // 									    wgt * binunc_scaling_dn);

      // }

    }

  }

  delete h2_etaPt;  // no need to save this as well
  cout << endl;
  cout << "Writing on output file" << endl;

  // if the file is opened in UPDATE mode, the following should overwrite an object if its key inside the file already exists
  outputFile->Write(0,TObject::kOverwrite);

  cout << "End of fillHistograms for " << sampleDir << endl;
  cout << endl;

}


//=============================================================

void loopNanoAOD(//const string& treedir = "/eos/cms/store/cmst3/group/wmass/mciprian/TREES_1LEP_80X_V3_WSKIM_NEW/", 
			 const string& treedir = "/eos/cms/store/cmst3/group/wmass/w-helicity-13TeV/trees/TREES_electrons_1l_V6_TINY/", 
			 const string& outdir = "./", 
			 const string& outfileName = "wmass_varhists.root") {

  createPlotDirAndCopyPhp(outdir);

  gROOT->SetBatch(kTRUE);
  TH1::SetDefaultSumw2(); //all the following histograms will automatically call TH1::Sumw2()  
  cout << endl;


  string cmssw_base = getEnvVariable("CMSSW_BASE");
  cout << "CMSSW_BASE = " << cmssw_base << endl;
  TFile* outputFile = new TFile((outdir + outfileName).c_str(),"RECREATE");
  if (!outputFile || outputFile->IsZombie()) {
    cout << "Error: file not opened. Exit" << endl;
    exit(EXIT_FAILURE);
  }

  outputFile->cd();

  fillHistograms(treedir, outdir, Sample::wmunujets, outputFile);
  
  outputFile->Close();
  delete outputFile;

  if (outdir != "./") {
    cout << "Going to copy this code for future reference in " << outdir << endl;     
    system(Form("cp %s/src/CMGTools/WMass/python/plotter/utilityMacros/src/loopNanoAOD.C %s",cmssw_base.c_str(),outdir.c_str()));
  }  

}
