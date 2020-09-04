#include "../interface/utility.h"
#include "../../functions.cc"
#include "../../w-mass-13TeV/functionsWMass.cc"

#define CHECK_EVERY_N 10000
#define N_MAX_ENTRIES_PER_SAMPLE 0 // for tests, use number <= 0 to use all events in each sample
#define FOLDER_IN_FILE 0
#define MT_CUT_WLIKE 45.0

using namespace std;

//static float intLumi = 30.9 // 35.9 for muons, 30.9 for electrons, measured in 1/fb
//static vector<Double_t> eleEtaBinEdges_double = {0.0, 1.0, 1.479, 2.1, 2.5};
// static vector<Double_t> etaBinEdgesTemplate = {-2.5,-2.4,-2.3,-2.2,-2.1,-2.0,-1.9,-1.8,-1.7,-1.6,-1.566,-1.5,-1.4442,-1.4,-1.3,-1.2,-1.1,-1.0,-0.9,-0.8,-0.7,-0.6,-0.5,-0.4,-0.3,-0.2,-0.1,0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.4442,1.5,1.566,1.6,1.7,1.8,1.9,2.0,2.1,2.2,2.3,2.4,2.5};
// static vector<Double_t> ptBinEdgesTemplate = {30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45};

// static vector<Double_t> genEtaBinEdgesTemplate = {0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.1,2.2,2.3,2.4};
// //static vector<Double_t> genPtBinEdgesTemplate = {26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45};
// static vector<Double_t> genPtBinEdgesTemplate = {30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45};

//static vector<Double_t> wptbins = {0.0, 1.971, 2.949, 3.838, 4.733, 5.674, 6.684, 7.781, 8.979, 10.303, 11.777, 13.435, 15.332, 17.525, 20.115, 23.245, 27.173, 32.414, 40.151, 53.858, 13000.0};
//static Int_t nPDFweight = 60;


//=============================================================

void fillHistograms(const string& treedir = "./", 
		    const string& outdir = "./", 
		    const Sample& sample = Sample::wjets, 
		    TFile* outputFile = NULL,
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

  bool isMuon = true;

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

  // muon
  vector<Double_t> etaBinEdgesTemplateMu = {-2.4,-2.3,-2.2,-2.1,-2.0,-1.9,-1.8,-1.7,-1.6,-1.5,-1.4,-1.3,-1.2,-1.1,-1.0,-0.9,
  					    -0.8,-0.7,-0.6,-0.5,-0.4,-0.3,-0.2,-0.1,0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,
  					    1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.1,2.2,2.3,2.4};
  // vector<Double_t> etaBinEdgesTemplateMu = {-2.4, -2.2, -2.0,-1.9,-1.8,-1.7,-1.6,-1.5,-1.4,-1.3,-1.2,-1.1,-1.0,-0.9,
  // 					    -0.8,-0.7,-0.6,-0.5,-0.4,-0.3,-0.2,-0.1,0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,
  // 					    1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9, 2.0, 2.2, 2.4};
  //vector<Double_t> etaBinEdgesTemplateMu = {-2.4, -2.2, -2.0,-1.8,-1.6,-1.4,-1.2,-1.1,-1.0,-0.9,
  //					    -0.8,-0.7,-0.6,-0.5,-0.4,-0.3,-0.2,-0.1,0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,
  //					    1.1,1.2,1.4,1.6,1.8,2.0, 2.2, 2.4};
  // vector<Double_t> etaBinEdgesTemplateMu = {-2.4, -2.1, -1.9,-1.7,-1.5,-1.3,-1.2,-1.1,-1.0,-0.9,
  // 					    -0.8,-0.7,-0.6,-0.5,-0.4,-0.3,-0.2,-0.1,0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,
  // 					    1.1,1.2,1.3,1.5,1.7,1.9, 2.1, 2.4};
  // vector<Double_t> etaBinEdgesTemplateMu = {-2.4, -2.1, -1.9,-1.7,-1.5,-1.3,-1.2,-1.1,-1.0,-0.9,
  // 					    -0.8,-0.7,-0.6,-0.5,-0.4,-0.3,-0.2,-0.1,0};
  // vector<Double_t> etaBinEdgesTemplateMu = {0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,
  // 					    1.1,1.2,1.3,1.5,1.7,1.9, 2.1, 2.4};
  // vector<Double_t> etaBinEdgesTemplateMu = {-2.4, -2.2, -2.0,-1.8,-1.6,-1.4,-1.2,-1.0,
  //  					    -0.8,-0.6,-0.4,-0.2,0,0.2,0.4,0.6,0.8,1.0,
  // 					    1.2,1.4,1.6,1.8,2.0, 2.2, 2.4};
  vector<Double_t> ptBinEdgesTemplateMu = {26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56};
  //vector<Double_t> ptBinEdgesTemplateMu = {26,28,30,32,34,36,38,40,42,44,46};
  //vector<Double_t> ptBinEdgesTemplateMu = {26,28,30,31.5,33,34.5,36,37.5,39.0,40.5,42,43.5,45,46.5,48,50,52,54,56};
  //vector<Double_t> ptBinEdgesTemplateMu = {25,26,28,30,31.5,33,34.5,36,37.5,39.0,40.5,42,43.5,45,46.5};
  
  vector<Double_t> etaBinEdgesTemplate;
  vector<Double_t> ptBinEdgesTemplate;

  vector<Double_t> *ptr = nullptr;

  ptr = &etaBinEdgesTemplateMu;
  for (UInt_t i = 0; i < ptr->size(); ++i) etaBinEdgesTemplate.push_back(ptr->at(i));
  etaBinEdgesTemplateMu.clear();

  ptr = &ptBinEdgesTemplateMu;
  for (UInt_t i = 0; i < ptr->size(); ++i) ptBinEdgesTemplate.push_back(ptr->at(i));
  ptBinEdgesTemplateMu.clear();
 
  Double_t etaCutReco = etaBinEdgesTemplate.back(); 
  Double_t ptMaxReco = ptBinEdgesTemplate.back(); 
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

  Bool_t noSumGenWeight = true;
  vector<Double_t> sumGenWeightVector;
  buildChain(chain, sumGenWeightVector, treedir, sample, friendChain, nameMatch, noSumGenWeight); 

  // change directory again, when building chain something was messed up
  if (FOLDER_IN_FILE) dirSample->cd();
  else                outputFile->cd();
  //  cout << "check" << endl;

  cout << "Setting TTreeReader and branches" << endl;
  TTreeReader reader (chain);

  //TTreeReaderValue<Int_t> isData(reader,"isData");
  TTreeReaderValue<UInt_t>    run(reader,"run");
  TTreeReaderValue<UInt_t>    lumi(reader,"lumi");
  TTreeReaderValue<ULong64_t> evt(reader,"evt");
  // //TTreeReaderValue<Int_t> nVert  (reader,"nVert");
  // //TTreeReaderValue<Float_t> rho  (reader,"rho");

  // // trigger
  TTreeReaderValue<Int_t>* HLT_SingleMuon_1 = new TTreeReaderValue<Int_t>(reader,"HLT_BIT_HLT_IsoMu24_v");
  TTreeReaderValue<Int_t>* HLT_SingleMuon_2 = new TTreeReaderValue<Int_t>(reader,"HLT_BIT_HLT_IsoTkMu24_v");
  // for trigger match
  TTreeReaderArray<Float_t> lep_trigMatchMuPt  (reader,"LepGood_matchedTrgObjMuPt");
  TTreeReaderArray<Float_t> lep_trigMatchTkMuPt(reader,"LepGood_matchedTrgObjTkMuPt");
  

  // // reco met
  // // TTreeReaderValue<Float_t> tkmet    (reader,"met_trkPt");
  // // TTreeReaderValue<Float_t> tkmet_phi(reader,"met_trkPhi");
  TTreeReaderValue<Float_t> pfmet    (reader,"met_pt");
  TTreeReaderValue<Float_t> pfmet_phi(reader,"met_phi");

  // // LepGood branch
  TTreeReaderValue<Int_t> nlep      (reader,"nLepGood");
  TTreeReaderArray<Int_t> lep_pdgId (reader,"LepGood_pdgId");
  TTreeReaderArray<Int_t> lep_charge(reader,"LepGood_charge");  
  TTreeReaderArray<Float_t> lep_calPt(reader,"LepGood_rocPt");
  //TTreeReaderArray<Float_t> lep_pt  (reader,"LepGood_pt");
  TTreeReaderArray<Float_t> lep_eta (reader,"LepGood_eta");
  TTreeReaderArray<Float_t> lep_phi (reader,"LepGood_phi");

  // muon ID
  TTreeReaderArray<Int_t> lep_Id  (reader, "LepGood_mediumMuonId");  
  // isolation
  TTreeReaderArray<Float_t> lep_iso  (reader, "LepGood_relIso04");  

  TTreeReaderArray<Int_t> *lep_mcMatchId = nullptr;
  // if (not isMuon) {
  //   lep_mcMatchId      = new TTreeReaderArray<Int_t>(reader,"LepGood_mcMatchId");
  // }

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
  TH2F* h2_etaPt = new TH2F("h2_etaPt","",netaBins,etaBinEdgesTemplate.data(),nptBins,ptBinEdgesTemplate.data());
  h2_etaPt->SetDirectory(0); // I don't want to save this histogram in the output file

  vector<string> parities = {"odd","even"}; // to split events

  vector<TH2F*> h2_parity_eta_pt_passIsopassMtpassMz; // pass both
  vector<TH2F*> h2_parity_eta_pt_passIsofailMtfailMz; // fail both
  vector<TH2F*> h2_parity_eta_pt_passIsopassMtfailMz;   // pass mT, fail mZ
  vector<TH2F*> h2_parity_eta_pt_passIsofailMtpassMz;   // fail mT, pass mZ
  vector<TH2F*> h2_parity_Mt_Mz_passIso;   // fail mT, pass mZ

  vector<TH2F*> h2_parity_eta_pt_failIsopassMtpassMz; // pass both
  vector<TH2F*> h2_parity_eta_pt_failIsofailMtfailMz; // fail both
  vector<TH2F*> h2_parity_eta_pt_failIsopassMtfailMz;   // pass mT, fail mZ
  vector<TH2F*> h2_parity_eta_pt_failIsofailMtpassMz;   // fail mT, pass mZ
  vector<TH2F*> h2_parity_Mt_Mz_failIso;   // fail mT, pass mZ

  Int_t nptBins1GeV = (Int_t) (ptMaxReco - ptMinReco + 0.00001);
  
  for (UInt_t ch = 0; ch < parities.size(); ++ch) {  

    h2_parity_eta_pt_passIsopassMtpassMz.push_back(new TH2F(Form("h2_%s_eta_pt_passIsopassMtpassMz",parities[ch].c_str()),"",48,-2.4,2.4,nptBins1GeV,ptMinReco,ptMaxReco));
    h2_parity_eta_pt_passIsofailMtfailMz.push_back(new TH2F(Form("h2_%s_eta_pt_passIsofailMtfailMz",parities[ch].c_str()),"",48,-2.4,2.4,nptBins1GeV,ptMinReco,ptMaxReco));
    h2_parity_eta_pt_passIsopassMtfailMz.push_back(new TH2F(Form("h2_%s_eta_pt_passIsopassMtfailMz",parities[ch].c_str()),"",48,-2.4,2.4,nptBins1GeV,ptMinReco,ptMaxReco));
    h2_parity_eta_pt_passIsofailMtpassMz.push_back(new TH2F(Form("h2_%s_eta_pt_passIsofailMtpassMz",parities[ch].c_str()),"",48,-2.4,2.4,nptBins1GeV,ptMinReco,ptMaxReco));
    h2_parity_Mt_Mz_passIso.push_back(new TH2F(Form("h2_%s_Mt_Mz_passIso",parities[ch].c_str()),"",60,0.0,120.0,60,0.0,120.0));

   h2_parity_eta_pt_failIsopassMtpassMz.push_back(new TH2F(Form("h2_%s_eta_pt_failIsopassMtpassMz",parities[ch].c_str()),"",48,-2.4,2.4,nptBins1GeV,ptMinReco,ptMaxReco));
    h2_parity_eta_pt_failIsofailMtfailMz.push_back(new TH2F(Form("h2_%s_eta_pt_failIsofailMtfailMz",parities[ch].c_str()),"",48,-2.4,2.4,nptBins1GeV,ptMinReco,ptMaxReco));
    h2_parity_eta_pt_failIsopassMtfailMz.push_back(new TH2F(Form("h2_%s_eta_pt_failIsopassMtfailMz",parities[ch].c_str()),"",48,-2.4,2.4,nptBins1GeV,ptMinReco,ptMaxReco));
    h2_parity_eta_pt_failIsofailMtpassMz.push_back(new TH2F(Form("h2_%s_eta_pt_failIsofailMtpassMz",parities[ch].c_str()),"",48,-2.4,2.4,nptBins1GeV,ptMinReco,ptMaxReco));
    h2_parity_Mt_Mz_failIso.push_back(new TH2F(Form("h2_%s_Mt_Mz_failIso",parities[ch].c_str()),"",60,0.0,120.0,60,0.0,120.0));


  }    


  cout << "Done histograms" << endl;

  ////////////////////
  ////////////////////

  // start event loop                                                                                                
  long int nTotal = chain->GetEntries();
  long int nEvents = 0;
  long int nEventsInSample = 0; // events processed for each sample

  Double_t absLep1eta = 0.0;
  Double_t absLep2eta = 0.0;

  ////////////////////////////////////////////
  // to get correct weight depending on sample in chain
  string currentFile = "";
  Int_t ifile = 0;
  ////////////////////

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
    // this does not skip the tree, but still loop on all events doing nothing
    if (N_MAX_ENTRIES_PER_SAMPLE > 0 && nEventsInSample > N_MAX_ENTRIES_PER_SAMPLE) continue; 

    //if (nEvents > 100) break;  // temporary line for tests
    //checkPoint(currentFile,200,nEvents);

    //cout << "check 0" << endl;

    // selection muons (SS leptons)
    if (**HLT_SingleMuon_1 < 1 and **HLT_SingleMuon_2 < 1) continue;
    if (*nlep != 2) continue;
    // now we are sure there are 2 leptons
    absLep1eta = fabs(lep_eta[0]);  
    absLep2eta = fabs(lep_eta[1]);  
    // same sign
    if (lep_pdgId[0]*lep_pdgId[1] != 169) continue; 
    // acceptance
    if (absLep1eta >= etaCutReco or absLep2eta >= etaCutReco) continue;
    if (lep_calPt[0] < ptMinReco or lep_calPt[0] > ptMaxReco) continue;
    if (lep_calPt[1] < ptMinReco or lep_calPt[1] > ptMaxReco) continue;
    // lepton ID
    if (lep_Id[0] < 1 or lep_Id[1] < 1) continue;
    // trigger match (at least one must match to one of the two triggers
    // this loop is for same charge for both leptons, so no matter which one is matched
    if (lep_trigMatchMuPt[0] <= 0.0 and lep_trigMatchTkMuPt[0] <= 0.0 and lep_trigMatchMuPt[1] <= 0.0 and lep_trigMatchTkMuPt[1] <= 0.0) continue;

    // mZ and mT cuts
    Bool_t passMz = false;
    Bool_t passMt = false;
    Bool_t passIso = false;
		    		   
    Double_t zmass = mass_2(lep_calPt[0],lep_eta[0],lep_phi[0],0.1057,lep_calPt[1],lep_eta[1],lep_phi[1],0.1057);
    Double_t mt_wlike = mt_wlike_samesign(lep_calPt[0],lep_phi[0],lep_calPt[1],lep_phi[1],*pfmet,*pfmet_phi);    

    // evaluate but do not cut, so to fill all histograms in same loop
    if (zmass > 60.0 and zmass < 120.0) passMz = true;
    if (lep_iso[0] < 0.15 and lep_iso[1] < 0.15) passIso = true;
    if (mt_wlike > 45.0) passMt = true;  
   
    Int_t parityIndex = 0;
    if (isEvenEvent(*evt)) parityIndex = 1;

    //cout << "check 1" << endl;
    
    // Now start filling histograms
    // fill twice, once for each lepton
    

    if (passIso) {
      h2_parity_Mt_Mz_passIso[parityIndex]->Fill(mt_wlike,zmass);
      h2_parity_Mt_Mz_passIso[parityIndex]->Fill(mt_wlike,zmass);
      if (passMz and passMt) { 
	h2_parity_eta_pt_passIsopassMtpassMz[parityIndex]->Fill(lep_eta[0],lep_calPt[0]);
	h2_parity_eta_pt_passIsopassMtpassMz[parityIndex]->Fill(lep_eta[1],lep_calPt[1]);
      }
      if (not passMz and not passMt) { 
	h2_parity_eta_pt_passIsofailMtfailMz[parityIndex]->Fill(lep_eta[0],lep_calPt[0]);
	h2_parity_eta_pt_passIsofailMtfailMz[parityIndex]->Fill(lep_eta[1],lep_calPt[1]);
      }
      if (passMz and not passMt) { 
	h2_parity_eta_pt_passIsofailMtpassMz[parityIndex]->Fill(lep_eta[0],lep_calPt[0]);
	h2_parity_eta_pt_passIsofailMtpassMz[parityIndex]->Fill(lep_eta[1],lep_calPt[1]);
      }
      if (not passMz and passMt) { 
	h2_parity_eta_pt_passIsopassMtfailMz[parityIndex]->Fill(lep_eta[0],lep_calPt[0]);
	h2_parity_eta_pt_passIsopassMtfailMz[parityIndex]->Fill(lep_eta[1],lep_calPt[1]);
      }
    } else {
      h2_parity_Mt_Mz_failIso[parityIndex]->Fill(mt_wlike,zmass);
      h2_parity_Mt_Mz_failIso[parityIndex]->Fill(mt_wlike,zmass);
      if (passMz and passMt) { 
	h2_parity_eta_pt_failIsopassMtpassMz[parityIndex]->Fill(lep_eta[0],lep_calPt[0]);
	h2_parity_eta_pt_failIsopassMtpassMz[parityIndex]->Fill(lep_eta[1],lep_calPt[1]);
      }
      if (not passMz and not passMt) { 
	h2_parity_eta_pt_failIsofailMtfailMz[parityIndex]->Fill(lep_eta[0],lep_calPt[0]);
	h2_parity_eta_pt_failIsofailMtfailMz[parityIndex]->Fill(lep_eta[1],lep_calPt[1]);
      }
      if (passMz and not passMt) { 
	h2_parity_eta_pt_failIsofailMtpassMz[parityIndex]->Fill(lep_eta[0],lep_calPt[0]);
	h2_parity_eta_pt_failIsofailMtpassMz[parityIndex]->Fill(lep_eta[1],lep_calPt[1]);
      }
      if (not passMz and passMt) { 
	h2_parity_eta_pt_failIsopassMtfailMz[parityIndex]->Fill(lep_eta[0],lep_calPt[0]);
	h2_parity_eta_pt_failIsopassMtfailMz[parityIndex]->Fill(lep_eta[1],lep_calPt[1]);
      }
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

// void fillHistogramsTH3(const string& treedir = "./", 
// 		    const string& outdir = "./", 
// 		    const Sample& sample = Sample::wjets, 
// 		    TFile* outputFile = NULL,
// 		    const string& nameMatch = ""
// 		    ) 
// {

//   gROOT->SetBatch(kTRUE);
//   TH1::SetDefaultSumw2(); //all the following histograms will automatically call TH1::Sumw2()  

//   cout << endl;
//   cout << "================================================" << endl;
//   cout << "Using integrated luminosity L = " << intLumi << endl;
//   cout << endl;

//   if (outputFile == NULL) {
//     cout << "Error: file is NULL. please check. Exit ..." << endl;
//     exit(EXIT_FAILURE);
//   }

//   outputFile->cd();

//   bool isMuon = true;

//   TDirectory *dirSample = NULL;
//   string sampleDir = getStringFromEnumSample(sample).c_str();
//   cout << "Sample --> " << sampleDir << endl;
//   if (FOLDER_IN_FILE) {
//     if (outputFile->GetKey(sampleDir.c_str())) dirSample = outputFile->GetDirectory(sampleDir.c_str());
//     else dirSample = outputFile->mkdir(sampleDir.c_str());
//     dirSample->cd();
//   }

//   cout << endl;

//   //============================
//   // BINNING: will have to implement parsing from file

//   // muon
//   vector<Double_t> etaBinEdgesTemplateMu = {-2.4,-2.3,-2.2,-2.1,-2.0,-1.9,-1.8,-1.7,-1.6,-1.5,-1.4,-1.3,-1.2,-1.1,-1.0,-0.9,
//   					    -0.8,-0.7,-0.6,-0.5,-0.4,-0.3,-0.2,-0.1,0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,
//   					    1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.1,2.2,2.3,2.4};
//   // vector<Double_t> etaBinEdgesTemplateMu = {-2.4, -2.2, -2.0,-1.9,-1.8,-1.7,-1.6,-1.5,-1.4,-1.3,-1.2,-1.1,-1.0,-0.9,
//   // 					    -0.8,-0.7,-0.6,-0.5,-0.4,-0.3,-0.2,-0.1,0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,
//   // 					    1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9, 2.0, 2.2, 2.4};
//   //vector<Double_t> etaBinEdgesTemplateMu = {-2.4, -2.2, -2.0,-1.8,-1.6,-1.4,-1.2,-1.1,-1.0,-0.9,
//   //					    -0.8,-0.7,-0.6,-0.5,-0.4,-0.3,-0.2,-0.1,0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,
//   //					    1.1,1.2,1.4,1.6,1.8,2.0, 2.2, 2.4};
//   // vector<Double_t> etaBinEdgesTemplateMu = {-2.4, -2.1, -1.9,-1.7,-1.5,-1.3,-1.2,-1.1,-1.0,-0.9,
//   // 					    -0.8,-0.7,-0.6,-0.5,-0.4,-0.3,-0.2,-0.1,0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,
//   // 					    1.1,1.2,1.3,1.5,1.7,1.9, 2.1, 2.4};
//   // vector<Double_t> etaBinEdgesTemplateMu = {-2.4, -2.1, -1.9,-1.7,-1.5,-1.3,-1.2,-1.1,-1.0,-0.9,
//   // 					    -0.8,-0.7,-0.6,-0.5,-0.4,-0.3,-0.2,-0.1,0};
//   // vector<Double_t> etaBinEdgesTemplateMu = {0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,
//   // 					    1.1,1.2,1.3,1.5,1.7,1.9, 2.1, 2.4};
//   // vector<Double_t> etaBinEdgesTemplateMu = {-2.4, -2.2, -2.0,-1.8,-1.6,-1.4,-1.2,-1.0,
//   //  					    -0.8,-0.6,-0.4,-0.2,0,0.2,0.4,0.6,0.8,1.0,
//   // 					    1.2,1.4,1.6,1.8,2.0, 2.2, 2.4};
//   vector<Double_t> ptBinEdgesTemplateMu = {26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56};
//   //vector<Double_t> ptBinEdgesTemplateMu = {26,28,30,32,34,36,38,40,42,44,46};
//   //vector<Double_t> ptBinEdgesTemplateMu = {26,28,30,31.5,33,34.5,36,37.5,39.0,40.5,42,43.5,45,46.5,48,50,52,54,56};
//   //vector<Double_t> ptBinEdgesTemplateMu = {25,26,28,30,31.5,33,34.5,36,37.5,39.0,40.5,42,43.5,45,46.5};
  
//   vector<Double_t> etaBinEdgesTemplate;
//   vector<Double_t> ptBinEdgesTemplate;

//   vector<Double_t> *ptr = nullptr;

//   ptr = &etaBinEdgesTemplateMu;
//   for (UInt_t i = 0; i < ptr->size(); ++i) etaBinEdgesTemplate.push_back(ptr->at(i));
//   etaBinEdgesTemplateMu.clear();

//   ptr = &ptBinEdgesTemplateMu;
//   for (UInt_t i = 0; i < ptr->size(); ++i) ptBinEdgesTemplate.push_back(ptr->at(i));
//   ptBinEdgesTemplateMu.clear();
 
//   Double_t etaCutReco = etaBinEdgesTemplate.back(); 
//   Double_t ptMaxReco = ptBinEdgesTemplate.back(); 
//   Double_t ptMinReco = ptBinEdgesTemplate[0];

//   //=============================


//   Int_t netaBins = etaBinEdgesTemplate.size() -1;
//   Int_t nptBins = ptBinEdgesTemplate.size() -1;
//   Int_t nBinsTemplate = nptBins * netaBins;
 
//   TChain* chain = new TChain("tree");
//   // INFO: the new friend trees at 13 TeV are inside a "tree_Friend_<sampleName>.root" file, and the tree's name is "Friends"
//   // friend trees are still located in a directory called "friends" with respect to base trees
//   TChain* friendChain = new TChain("Friends");      
//   //TChain* friendChain = NULL;  // leave as NULL if you don't use friend trees

//   Bool_t noSumGenWeight = true;
//   vector<Double_t> sumGenWeightVector;
//   buildChain(chain, sumGenWeightVector, treedir, sample, friendChain, nameMatch, noSumGenWeight); 

//   // change directory again, when building chain something was messed up
//   if (FOLDER_IN_FILE) dirSample->cd();
//   else                outputFile->cd();
//   //  cout << "check" << endl;

//   cout << "Setting TTreeReader and branches" << endl;
//   TTreeReader reader (chain);

//   //TTreeReaderValue<Int_t> isData(reader,"isData");
//   TTreeReaderValue<UInt_t>    run(reader,"run");
//   TTreeReaderValue<UInt_t>    lumi(reader,"lumi");
//   TTreeReaderValue<ULong64_t> evt(reader,"evt");
//   // //TTreeReaderValue<Int_t> nVert  (reader,"nVert");
//   // //TTreeReaderValue<Float_t> rho  (reader,"rho");

//   // // trigger
//   TTreeReaderValue<Int_t>* HLT_SingleMuon_1 = new TTreeReaderValue<Int_t>(reader,"HLT_BIT_HLT_IsoMu24_v");
//   TTreeReaderValue<Int_t>* HLT_SingleMuon_2 = new TTreeReaderValue<Int_t>(reader,"HLT_BIT_HLT_IsoTkMu24_v");
//   // for trigger match
//   TTreeReaderArray<Float_t> lep_trigMatchMuPt  (reader,"LepGood_matchedTrgObjMuPt");
//   TTreeReaderArray<Float_t> lep_trigMatchTkMuPt(reader,"LepGood_matchedTrgObjTkMuPt");
  

//   // // reco met
//   // // TTreeReaderValue<Float_t> tkmet    (reader,"met_trkPt");
//   // // TTreeReaderValue<Float_t> tkmet_phi(reader,"met_trkPhi");
//   TTreeReaderValue<Float_t> pfmet    (reader,"met_pt");
//   TTreeReaderValue<Float_t> pfmet_phi(reader,"met_phi");

//   // // LepGood branch
//   TTreeReaderValue<Int_t> nlep      (reader,"nLepGood");
//   TTreeReaderArray<Int_t> lep_pdgId (reader,"LepGood_pdgId");
//   TTreeReaderArray<Int_t> lep_charge(reader,"LepGood_charge");  
//   TTreeReaderArray<Float_t> lep_calPt(reader,"LepGood_rocPt");
//   //TTreeReaderArray<Float_t> lep_pt  (reader,"LepGood_pt");
//   TTreeReaderArray<Float_t> lep_eta (reader,"LepGood_eta");
//   TTreeReaderArray<Float_t> lep_phi (reader,"LepGood_phi");

//   // muon ID
//   TTreeReaderArray<Int_t> lep_Id  (reader, "LepGood_mediumMuonId");  
//   // isolation
//   TTreeReaderArray<Float_t> lep_iso  (reader, "LepGood_relIso04");  

//   TTreeReaderArray<Int_t> *lep_mcMatchId = nullptr;
//   // if (not isMuon) {
//   //   lep_mcMatchId      = new TTreeReaderArray<Int_t>(reader,"LepGood_mcMatchId");
//   // }

//   cout << "Branches have been set" << endl;

//   //////////////////////////
//   // NOTE
//   // values of TTreeReaderValue objects are read with * before the variable name, as if they were pointers
//   // Pointers to TTreeReaderValue are used with double * (one to access the pointer's content and the other for the TTreeReaderValue convention)
//   // TTreeReaderArray variables are used as normal array variables (if they are pointers to TTreeReaderArray.then you still need a * to access the array)
//   //////////////////////////

//   cout << "Defining histograms" << endl;

//   /////////////////////////////////////
//   // dummy histogram to easily retrieve information on eta-pt bins
//   TH2F* h2_etaPt = new TH2F("h2_etaPt","",netaBins,etaBinEdgesTemplate.data(),nptBins,ptBinEdgesTemplate.data());
//   h2_etaPt->SetDirectory(0); // I don't want to save this histogram in the output file

//   vector<string> parities = {"odd","even"}; // to split events
//   //vector<TH1F*> h1_charge_eta;
//   //vector<TH1F*> h1_charge_pt;
//   //vector<TH3F*> h2_charge_eta_pt_iso; // inclusive, before mT and mZ cut, can be taken as sum of other ones below
//   vector<TH3F*> h3_parity_eta_pt_iso_passMtpassMz; // pass both
//   vector<TH3F*> h3_parity_eta_pt_iso_failMtfailMz; // fail both
//   vector<TH3F*> h3_parity_eta_pt_iso_passMtfailMz;   // pass mT, fail mZ
//   vector<TH3F*> h3_parity_eta_pt_iso_failMtpassMz;   // fail mT, pass mZ
//   vector<TH3F*> h3_parity_Mt_Mz_iso;   // fail mT, pass mZ

//   Int_t nptBins1GeV = (Int_t) (ptMaxReco - ptMinReco + 0.00001);
  
//   for (UInt_t ch = 0; ch < parities.size(); ++ch) {  

//     //h1_charge_pt.push_back(new TH1F(Form("h1_%s_pt",charges[ch].c_str()),"",nptBins1GeV,ptMinReco,ptMaxReco));
//     //h1_charge_eta.push_back(new TH1F(Form("h1_%s_eta",charges[ch].c_str()),"",48,-2.4,2.4));

//     // h3_charge_eta_pt_iso.push_back(new TH3F(Form("h3_%s_eta_pt_iso",charges[ch].c_str()),"",
//     // 					    48,-2.4,2.4,
//     // 					    nptBins1GeV,ptMinReco,ptMaxReco,
//     // 					    80,0.0,0.8));
//     h3_parity_eta_pt_iso_passMtpassMz.push_back(new TH3F(Form("h3_%s_eta_pt_iso_passMtpassMz",parities[ch].c_str()),"",
// 							 48,-2.4,2.4,
// 							 nptBins1GeV,ptMinReco,ptMaxReco,
// 							 80,0.0,0.8));
//     h3_parity_eta_pt_iso_failMtfailMz.push_back(new TH3F(Form("h3_%s_eta_pt_iso_failMtfailMz",parities[ch].c_str()),"",
// 						     48,-2.4,2.4,
// 							 nptBins1GeV,ptMinReco,ptMaxReco,
// 							 80,0.0,0.8));
//     h3_parity_eta_pt_iso_passMtfailMz.push_back(new TH3F(Form("h3_%s_eta_pt_iso_passMtfailMz",parities[ch].c_str()),"",
// 							 48,-2.4,2.4,
// 							 nptBins1GeV,ptMinReco,ptMaxReco,
// 							 80,0.0,0.8));
//     h3_parity_eta_pt_iso_failMtpassMz.push_back(new TH3F(Form("h3_%s_eta_pt_iso_failMtpassMz",parities[ch].c_str()),"",
// 							 48,-2.4,2.4,
// 							 nptBins1GeV,ptMinReco,ptMaxReco,
// 							 80,0.0,0.8));
//     h3_parity_Mt_Mz_iso.push_back(new TH3F(Form("h3_%s_Mt_Mz_iso",parities[ch].c_str()),"",
// 					   60,0.0,120.0,
// 					   60,0.0,120.0,
// 					   80,0.0,0.8));

//   }    


//   cout << "Done histograms" << endl;

//   ////////////////////
//   ////////////////////

//   // start event loop                                                                                                
//   long int nTotal = chain->GetEntries();
//   long int nEvents = 0;
//   long int nEventsInSample = 0; // events processed for each sample

//   Double_t absLep1eta = 0.0;
//   Double_t absLep2eta = 0.0;

//   ////////////////////////////////////////////
//   // to get correct weight depending on sample in chain
//   string currentFile = "";
//   Int_t ifile = 0;
//   ////////////////////

//   cout << "Starting loop" << endl;

//   //cout << "check -2" << endl;
//   while (reader.Next()) {
  
//     //cout << "check -1" << endl;
//     cout.flush();
//     if (nEvents % CHECK_EVERY_N == 0) cout << "\r" << "Analyzing events " << ((Double_t) nEvents)/nTotal*100.0 << " % ";
//     //cout << "entry : " << nEvents << endl;
//     nEvents++;
//     nEventsInSample++;

//     if (dynamic_cast<TChain*>(reader.GetTree())->GetFile()->GetName() != currentFile and currentFile != "") { 
//       currentFile = dynamic_cast<TChain*>(reader.GetTree())->GetFile()->GetName();                   
//       ifile ++;                                                                                      
//       nEventsInSample = 1; // reset nEvents when sub sample is changed (useful with N_MAX_ENTRIES_PER_SAMPLE for debugging)
//     } else if (dynamic_cast<TChain*>(reader.GetTree())->GetFile()->GetName() != currentFile and currentFile == "") {
//       currentFile = dynamic_cast<TChain*>(reader.GetTree())->GetFile()->GetName();                         
//     }      
//     // this does not skip the tree, but still loop on all events doing nothing
//     if (N_MAX_ENTRIES_PER_SAMPLE > 0 && nEventsInSample > N_MAX_ENTRIES_PER_SAMPLE) continue; 

//     //if (nEvents > 100) break;  // temporary line for tests
//     //checkPoint(currentFile,200,nEvents);

//     //cout << "check 0" << endl;

//     // selection muons (SS leptons)
//     if (**HLT_SingleMuon_1 < 1 and **HLT_SingleMuon_2 < 1) continue;
//     if (*nlep != 2) continue;
//     // now we are sure there are 2 leptons
//     absLep1eta = fabs(lep_eta[0]);  
//     absLep2eta = fabs(lep_eta[1]);  
//     // same sign
//     if (lep_pdgId[0]*lep_pdgId[1] != 169) continue; 
//     // acceptance
//     if (absLep1eta >= etaCutReco or absLep2eta >= etaCutReco) continue;
//     if (lep_calPt[0] < ptMinReco or lep_calPt[0] > ptMaxReco) continue;
//     if (lep_calPt[1] < ptMinReco or lep_calPt[1] > ptMaxReco) continue;
//     // lepton ID
//     if (lep_Id[0] < 1 or lep_Id[1] < 1) continue;
//     // trigger match (at least one must match to one of the two triggers
//     // this loop is for same charge for both leptons, so no matter which one is matched
//     if (lep_trigMatchMuPt[0] <= 0.0 and lep_trigMatchTkMuPt[0] <= 0.0 and lep_trigMatchMuPt[1] <= 0.0 and lep_trigMatchTkMuPt[1] <= 0.0) continue;

//     // mZ and mT cuts
//     Bool_t passMz = false;
//     Bool_t passMt = false;
//     //Bool_t passIso = false;
		    		   
//     Double_t zmass = mass_2(lep_calPt[0],lep_eta[0],lep_phi[0],0.1057,lep_calPt[1],lep_eta[1],lep_phi[1],0.1057);
//     Double_t mt_wlike = mt_wlike_samesign(lep_calPt[0],lep_phi[0],lep_calPt[1],lep_phi[1],*pfmet,*pfmet_phi);    

//     // evaluate but do not cut, so to fill all histograms in same loop
//     if (zmass > 60.0 and zmass < 120.0) passMz = true;
//     //if (lep_iso[0] < 0.15 and lep_iso[1] < 0.15) passIso = true;
//     if (mt_wlike > 45.0) passMt = true;  
   
//     Int_t parityIndex = 0;
//     if (isEvenEvent(*evt)) parityIndex = 1;

//     //cout << "check 1" << endl;
    
//     // Now start filling histograms
//     // fill twice, once for each lepton
//     // h3_charge_eta_pt_iso[parityIndex]->Fill(lep_eta[0],lep_calPt[0], lep_iso[0]);
//     // h3_charge_eta_pt_iso[parityIndex]->Fill(lep_eta[1],lep_calPt[1], lep_iso[1]);
    
//     //Double_t lep1iso = std::max(0,std::min(lep_iso[0]),0.8); // to put large iso in last bin

//     h3_parity_Mt_Mz_iso[parityIndex]->Fill(mt_wlike,zmass,lep_iso[0]);
//     h3_parity_Mt_Mz_iso[parityIndex]->Fill(mt_wlike,zmass,lep_iso[1]);

//     if (passMz and passMt) { 
// 	h3_parity_eta_pt_iso_passMtpassMz[parityIndex]->Fill(lep_eta[0],lep_calPt[0], lep_iso[0]);
// 	h3_parity_eta_pt_iso_passMtpassMz[parityIndex]->Fill(lep_eta[1],lep_calPt[1], lep_iso[1]);
//     }
//     if (not passMz and not passMt) { 
// 	h3_parity_eta_pt_iso_failMtfailMz[parityIndex]->Fill(lep_eta[0],lep_calPt[0], lep_iso[0]);
// 	h3_parity_eta_pt_iso_failMtfailMz[parityIndex]->Fill(lep_eta[1],lep_calPt[1], lep_iso[1]);
//     }
//     if (passMz and not passMt) { 
// 	h3_parity_eta_pt_iso_failMtpassMz[parityIndex]->Fill(lep_eta[0],lep_calPt[0], lep_iso[0]);
// 	h3_parity_eta_pt_iso_failMtpassMz[parityIndex]->Fill(lep_eta[1],lep_calPt[1], lep_iso[1]);
//     }
//     if (not passMz and passMt) { 
// 	h3_parity_eta_pt_iso_passMtfailMz[parityIndex]->Fill(lep_eta[0],lep_calPt[0], lep_iso[0]);
// 	h3_parity_eta_pt_iso_passMtfailMz[parityIndex]->Fill(lep_eta[1],lep_calPt[1], lep_iso[1]);
//     }

//   }

//   delete h2_etaPt;  // no need to save this as well
//   cout << endl;
//   cout << "Writing on output file" << endl;

//   // if the file is opened in UPDATE mode, the following should overwrite an object if its key inside the file already exists
//   outputFile->Write(0,TObject::kOverwrite);

//   cout << "End of fillHistograms for " << sampleDir << endl;
//   cout << endl;

// }


//=============================================================

void loopNtuples_fakes(//const string& treedir = "/eos/cms/store/cmst3/group/wmass/mciprian/TREES_1LEP_80X_V3_WSKIM_NEW/", 
		       const string& treedir = "/eos/cms/store/cmst3/group/wmass/w-mass-13TeV//ntuples/SKIM_2LEP_wlike_mu_V2/", 
		       const string& outdir = "./", 
		       const string& outfileName = "wlike_qcdhists.root",
		       const bool isMuon = true
		       ) {

  createPlotDirAndCopyPhp(outdir);

  gROOT->SetBatch(kTRUE);
  TH1::SetDefaultSumw2(); //all the following histograms will automatically call TH1::Sumw2()  
  cout << endl;


  string cmssw_base = getEnvVariable("CMSSW_BASE");
  cout << "CMSSW_BASE = " << cmssw_base << endl;

  cout << "Loading functions.cc" << endl;
  gROOT->ProcessLine(Form(".L %s/src/CMGTools/WMass/python/plotter/functions.cc+",cmssw_base.c_str())); 
  cout << "Loading functionsWMass.cc" << endl;
  gROOT->ProcessLine(Form(".L %s/src/CMGTools/WMass/python/plotter/w-mass-13TeV/functionsWMass.cc+",cmssw_base.c_str()));

  // //compile some stuff to get functions used in heppy for some reweightings
  // string rootLibraries = string(gSystem->GetLibraries());                                                                                       
  // //cout << "gSystem->GetLibraries():" << endl;
  // //cout << gSystem->GetLibraries() << endl;
  // if (rootLibraries.find("/w-mass-13TeV/functionsWMass_cc.so") == string::npos) 
  //   compileMacro("src/CMGTools/WMass/python/plotter/w-mass-13TeV/functionsWMass.cc");
  // if (rootLibraries.find("/functions_cc.so") == string::npos)
  //   compileMacro("src/CMGTools/WMass/python/plotter/functions.cc");                               

  // gROOT->ProcessLine(Form(".L %s/src/CMGTools/WMass/python/plotter/functions.cc+",cmssw_base.c_str()));
  // gROOT->ProcessLine(Form(".L %s/src/CMGTools/WMass/python/plotter/w-mass-13TeV/functionsWMass.cc+",cmssw_base.c_str()));

  // cout << "PU(0) = " << puw2016_nTrueInt_36fb(0) << endl;
  // cout << "PU(1) = " << puw2016_nTrueInt_36fb(1) << endl;
  // cout << "PU(5) = " << puw2016_nTrueInt_36fb(5) << endl;
  // cout << "PU(30) = " << puw2016_nTrueInt_36fb(30) << endl;
  // cout << "PU(99) = " << puw2016_nTrueInt_36fb(99) << endl;
  // cout << "PU(100) = " << puw2016_nTrueInt_36fb(100) << endl;
  // cout << "TEST SUCCESSFUL" << endl;
  // return;

  TFile* outputFile = new TFile((outdir + outfileName).c_str(),"RECREATE");
  //TFile* outputFile = new TFile((outdir + outfileName).c_str(),"UPDATE");
  if (!outputFile || outputFile->IsZombie()) {
    cout << "Error: file not opened. Exit" << endl;
    exit(EXIT_FAILURE);
  }

  outputFile->cd();

  //fillHistograms(treedir, outdir, Sample::wenujets, outputFile,"WJetsToLNu_NLO_ext2v5_part1");
  //fillHistograms(treedir, outdir, Sample::wmunujets, outputFile, usePreFSRvar, "WJetsToLNu_ext2v5_part1");
  //fillHistograms(treedir, outdir, Sample::wmunujets, outputFile, usePreFSRvar);

  //fillHistograms(treedir, outdir, Sample::data_singleEG, outputFile);
  fillHistograms(treedir, outdir, Sample::data_singleMu, outputFile);
  // if (useFakeRateForElectron) 
  //   fillHistograms(treedir, outdir, Sample::qcd_ele_fake, outputFile); // use fake rate only in SR
  // fillHistograms(treedir, outdir, Sample::qcd_ele, outputFile);
  //  fillHistograms(treedir, outdir, Sample::wenujets, outputFile);
  // fillHistograms(treedir, outdir, Sample::wtaunujets, outputFile);
  // fillHistograms(treedir, outdir, Sample::zjets, outputFile);
  // fillHistograms(treedir, outdir, Sample::top, outputFile);
  // fillHistograms(treedir, outdir, Sample::diboson, outputFile);
  
  outputFile->Close();
  delete outputFile;

  if (outdir != "./") {
    cout << "Going to copy this code for future reference in " << outdir << endl;     
    system(Form("cp %s/src/CMGTools/WMass/python/plotter/utilityMacros/src/loopNtuples_fakes.C %s",cmssw_base.c_str(),outdir.c_str()));
  }  

}
