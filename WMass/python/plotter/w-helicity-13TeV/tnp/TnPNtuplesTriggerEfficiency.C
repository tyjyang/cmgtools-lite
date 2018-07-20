#define TnPNtuplesTriggerEfficiency_cxx
#include "TnPNtuplesTriggerEfficiency.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <iostream>
#include <TLorentzVector.h>

using namespace std;

void TnPNtuplesTriggerEfficiency::setFlavor(int flavor){
  fFlavor = flavor;
}

void TnPNtuplesTriggerEfficiency::setOutfile(TString outfilepath){
  fOutfile = outfilepath;
}

bool TnPNtuplesTriggerEfficiency::isTagLepton(int jj){

  if(fFlavor == 11){
    if (abs(LepGood_pdgId[jj])!=11)                                  return false;
    if (LepGood_calPt[jj]<30)                                        return false;
    if (fabs(LepGood_eta[jj])>1.4442 && fabs(LepGood_eta[jj])<1.566) return false;
    if (fabs(LepGood_eta[jj])>2.5)                                   return false;
    if (LepGood_customId[jj] < 1)                                    return false;
    if (LepGood_hltId[jj] < 1)                                       return false;
    if (LepGood_tightChargeFix[jj] != 2 )                            return false;
  }

  else {
    if (abs(LepGood_pdgId[jj])!=13)   return false;
    if (LepGood_pt[jj]<26)            return false;
    if (fabs(LepGood_eta[jj])> 2.4)   return false;
    if (LepGood_relIso04[jj] > 0.15)  return false;
    if (LepGood_mediumMuonId[jj] < 1) return false;
  }
  return true;

}

void TnPNtuplesTriggerEfficiency::Loop(int maxentries)
{
  // skim
  // HLT_SingleEL : HLT_SingleEl == 1
  // onelep : nLepGood == 1 && abs(LepGood1_pdgId)==11
  // fiducial : abs(LepGood1_eta)<1.4442 || abs(LepGood1_eta)>1.566
  // eleKin : ptElFull(LepGood1_calPt,LepGood1_eta) > 30 && abs(LepGood1_eta)<2.5
  // HLTid : LepGood1_hltId > 0
  // numSel : LepGood1_customId == 1
  
  // -----------------------

  bool doElectrons = (fFlavor == 11);
  if (doElectrons) std::cout << "running on electrons !!! " << std::endl;
  else             std::cout << "running on muons !!! " << std::endl;
  // -----------------------

  // Start of the program
  if (fChain == 0) return;
  Long64_t nentries = fChain->GetEntriesFast();
  Long64_t nbytes = 0, nb = 0;

  // Booking histos and tree with final variables
  bookOutputTree();

  // Vectors to store infos
  vector <float> cand_pt        = {};
  vector <float> cand_eta       = {};
  vector <float> cand_etaSc     = {};
  vector <float> cand_phi       = {};
  vector <float> cand_eleTrgPt  = {};
  vector <float> cand_muTrgPt   = {};
  vector <float> cand_tkMuTrgPt = {};
  vector <int>   cand_matchMC   = {};
  vector <int>   cand_hltSafeId = {};
  vector <int>   cand_customId  = {};
  vector <int>   cand_alsoTag   = {};
  vector <int>   cand_isZero    = {};

  // To compute the lumi weight
  float sigma=1921.8*3.;
  float count_getentries = doElectrons ? 123847915 : 99999999999;    // madgraph, ext1+ext2 //MARC THIS NUMBER IS WRONG FOR MUONS
  float SetLumi=35.9;     


  // Loop over events
  std::cout << "Start looping over events" << endl;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    if (maxentries > 0 && ientry >= maxentries) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    if (!(ientry%250000)) std::cout << ientry << endl;
    thisEntry = ientry;

    // To keep track of the total number of events
    h_entries->Fill(5);  

    // PU weight
    mypuw = 1.;
    if (!isData) mypuw = puw2016_nTrueInt_36fb(nTrueInt);
    
    // Lumi weight for MC only
    totWeight = 1.;
    if (!isData) totWeight = mypuw*genWeight*((sigma*SetLumi)/count_getentries)*(pow(10,3));  

    // Events breakdown  
    h_selection->Fill(0.);

    // 1) analysis cuts: fire the single lepton trigger  
    if ( doElectrons && !HLT_BIT_HLT_Ele27_WPTight_Gsf_v) continue;
    if (!doElectrons && !HLT_BIT_HLT_IsoMu24_v && !HLT_BIT_HLT_IsoTkMu24_v) continue; // require the tag to fire both to avoid bias

    h_selection->Fill(1.);
    
    // 2) at least one good vertex found
    nvtx = nVert;
    if (nvtx<1) continue;
    h_selection->Fill(2.);

    // 3) Gen level match, to be saved as additional info
    int numGenLevel=0;
    bool genLepFound = false;
    bool genPosFound = false;
    TLorentzVector myGenLep(0,0,0,0);  
    TLorentzVector myGenPos(0,0,0,0);  
    if (!isData) {   
      for(int ii=0; ii<nGenPart; ii++){
        int status = GenPart_status[ii];
        int pdgid  = GenPart_pdgId[ii];
        if ( abs(pdgid)==fFlavor && status==23 ) { // check if lepton is correct type and status 23
          if (GenPart_motherId[ii]==23){
            float ptgen  = GenPart_pt[ii];
            float etagen = GenPart_eta[ii];
            float phigen = GenPart_phi[ii];
            if (pdgid>0)  { //marc: is this correct? i though 11 was a negatively charged lepton. shouldn't matter though
              myGenPos.SetPtEtaPhiM(ptgen, etagen, phigen, 0.);
              genPosFound = true;
            }
            if (pdgid<0) {
              myGenLep.SetPtEtaPhiM(ptgen, etagen, phigen, 0.);
              genLepFound = true;
            }
            numGenLevel++;
          } // end check motherid == Z
        } // end if pdgId && status
      } // end loop gen particles
    } // end MC only


    // 4) Tag and probe selection
    std::vector<int> acceptLep;
    
    // full selection, ID+ISO, but not trigger match requirement
    for(int jj=0; jj<nLepGood; jj++){
      if ( isTagLepton(jj) ) acceptLep.push_back(jj);
    }

    bool atLeastOneTag = false;
    // loop on all the ID+ISO leptons
    for (std::vector<int>::const_iterator ilep = acceptLep.begin(); ilep !=acceptLep.end(); ++ilep){ 
      int theOrigIndex = *ilep;

      // kine 
      float lepPt    = fFlavor == 11 ? LepGood_calPt [theOrigIndex] : LepGood_pt[theOrigIndex];     // calibrated pT for electrons
      float lepEta   = LepGood_eta   [theOrigIndex];
      float lepScEta = LepGood_etaSc [theOrigIndex];
      float lepPhi   = LepGood_phi   [theOrigIndex];

      // this lep
      TLorentzVector thisRecoLep(0,0,0,0);
      thisRecoLep.SetPtEtaPhiM(lepPt,lepEta,lepPhi,0);

      // Match with MC truth   
      int matchMC = 0;
      if (!isData) {  
        if(genLepFound && thisRecoLep.DeltaR(myGenLep)<0.3) matchMC = 1;  
        if(genPosFound && thisRecoLep.DeltaR(myGenPos)<0.3) matchMC = 1;  
      } 
      else {
        matchMC = 1;
      } // end check mcMatch

      // HLT Safe ID
      int hltSafeId = LepGood_hltId[theOrigIndex];
      
      // Full ID 
      int customId = LepGood_customId[theOrigIndex];

      // Is this a tag:
      int isThisTag = 0;
      if (fFlavor == 11){
        if (isTagLepton(theOrigIndex) && LepGood_matchedTrgObjElePt[theOrigIndex] > -1.) isThisTag =1;
      }
      else{
        if (isTagLepton(theOrigIndex) && LepGood_matchedTrgObjMuPt[theOrigIndex] > -1. && LepGood_matchedTrgObjTkMuPt[theOrigIndex] > -1.) isThisTag =1; // require muon tag to fire both triggers.
      }

      if (isThisTag) atLeastOneTag = true;

      // Infos to be kept
      cand_pt        . push_back(lepPt);
      cand_eta       . push_back(lepEta);
      cand_etaSc     . push_back(lepScEta);
      cand_phi       . push_back(lepPhi);
      cand_eleTrgPt  . push_back(LepGood_matchedTrgObjElePt[theOrigIndex]);
      cand_muTrgPt   . push_back(LepGood_matchedTrgObjMuPt[theOrigIndex]);
      cand_tkMuTrgPt . push_back(LepGood_matchedTrgObjTkMuPt[theOrigIndex]);
      cand_matchMC   . push_back(matchMC);
      cand_hltSafeId . push_back(hltSafeId);
      cand_customId  . push_back(customId);
      cand_alsoTag   . push_back(isThisTag);

      if (theOrigIndex==0) cand_isZero.push_back(1);
      else cand_isZero.push_back(0);    
      
    } // end leptons in acceptance

    
    //--- 4) at least one tag candidate - should be there by definition of the skim
    if (!atLeastOneTag) {
      // cleaning vectors
      cand_pt        . clear();
      cand_eta       . clear();
      cand_etaSc     . clear();
      cand_phi       . clear();
      cand_eleTrgPt  . clear();
      cand_muTrgPt   . clear();
      cand_tkMuTrgPt . clear();
      cand_matchMC   . clear();
      cand_hltSafeId . clear();
      cand_customId  . clear();
      cand_alsoTag   . clear();
      cand_isZero    . clear();
      continue;
    }
    h_selection->Fill(3.);    
  
    //--- 5) invariant mass and T&P pairs
    // first as tag  
    for(unsigned int iLep1=0; iLep1<cand_pt.size(); ++iLep1) {
      if (!cand_alsoTag[iLep1]) continue;
      TLorentzVector thisLep1(0,0,0,0);   
      thisLep1.SetPtEtaPhiM(cand_pt[iLep1],cand_eta[iLep1],cand_phi[iLep1],0);

      // second as probe 
      for(unsigned int iLep2=0; iLep2<cand_pt.size(); ++iLep2) {
        if (cand_isZero[iLep2]) continue;
        TLorentzVector thisLep2(0,0,0,0);
        thisLep2.SetPtEtaPhiM(cand_pt[iLep2],cand_eta[iLep2],cand_phi[iLep2],0);

        // invariant mass
        pair_mass = (thisLep1+thisLep2).M();
        if (pair_mass<60 || pair_mass>120) continue;
        
        // both matching mc truth?
        mcTrue = cand_matchMC[iLep1] && cand_matchMC[iLep2];
        
        // first as tag, second as probe
        tag_lep_pt          = cand_pt        [iLep1];
        tag_lep_eta         = cand_eta       [iLep1];
        tag_lep_matchMC     = cand_matchMC   [iLep1];

        probe_lep_pt        = cand_pt        [iLep2];
        probe_lep_eta       = cand_eta       [iLep2];
        probe_sc_eta        = cand_etaSc     [iLep2];
        probe_lep_phi       = cand_phi       [iLep2];
        probe_eleTrgPt      = cand_eleTrgPt  [iLep2];
        probe_muTrgPt       = cand_muTrgPt   [iLep2];
        probe_tkMuTrgPt     = cand_tkMuTrgPt [iLep2];
        probe_lep_matchMC   = cand_matchMC   [iLep2];
        probe_lep_hltSafeId = cand_hltSafeId [iLep2];
        probe_lep_customId  = cand_customId  [iLep2];
        probe_lep_alsoTag   = cand_alsoTag   [iLep2];
        
        // Tree filling
        outFile_->cd();
        cddir->cd();  
        outTree_->Fill();

      }  // probes
    }   // tags
  
    // cleaning vectors
    cand_pt        . clear();
    cand_eta       . clear();
    cand_etaSc     . clear();
    cand_phi       . clear();
    cand_matchMC   . clear();
    cand_hltSafeId . clear();
    cand_customId  . clear();
    cand_alsoTag   . clear();
    cand_isZero    . clear();
    cand_eleTrgPt  . clear();
    cand_muTrgPt   . clear();
    cand_tkMuTrgPt . clear();

  }  // Loop over entries

  // Saving output tree and histos
  outFile_    -> cd();
  h_entries   -> Write();
  h_selection -> Write();
  cddir       -> cd();
  outTree_    -> Write();

} // Loop method

// Histos booking
void TnPNtuplesTriggerEfficiency::bookOutputTree() 
{
  outFile_ = new TFile(fOutfile, "RECREATE");    
  outFile_->cd();

  cddir = outFile_->mkdir("IDIsoToHLT");
  cddir->cd();    

  std::cout << "Booking output tree" << endl;
  outTree_ = new TTree("fitter_tree", "fitter_tree");

  outTree_->Branch("tag_lep_pt"          , &tag_lep_pt          , "tag_lep_pt/F");
  outTree_->Branch("tag_lep_eta"         , &tag_lep_eta         , "tag_lep_eta/F");
  outTree_->Branch("tag_lep_matchMC"     , &tag_lep_matchMC     , "tag_lep_matchMC/I");
  outTree_->Branch("probe_lep_pt"        , &probe_lep_pt        , "probe_lep_pt/F");
  outTree_->Branch("probe_lep_eta"       , &probe_lep_eta       , "probe_lep_eta/F");
  outTree_->Branch("probe_sc_eta"        , &probe_sc_eta        , "probe_sc_eta/F");
  outTree_->Branch("probe_lep_phi"       , &probe_lep_phi       , "probe_lep_phi/F");

  outTree_->Branch("probe_eleTrgPt"      , &probe_eleTrgPt      , "probe_eleTrgPt/F");
  outTree_->Branch("probe_muTrgPt"       , &probe_muTrgPt       , "probe_muTrgPt/F");
  outTree_->Branch("probe_tkMuTrgPt"     , &probe_tkMuTrgPt     , "probe_tkMuTrgPt/F");

  outTree_->Branch("probe_lep_matchMC"   , &probe_lep_matchMC   , "probe_lep_matchMC/I");
  outTree_->Branch("probe_lep_hltSafeId" , &probe_lep_hltSafeId , "probe_lep_hltSafeId/I");
  outTree_->Branch("probe_lep_customId"  , &probe_lep_customId  , "probe_lep_customId/I");
  outTree_->Branch("probe_lep_alsoTag"   , &probe_lep_alsoTag   , "probe_lep_alsoTag/I");
  outTree_->Branch("pair_mass"           , &pair_mass           , "pair_mass/F");
  outTree_->Branch("nvtx"                , &nvtx                , "nvtx/I");
  outTree_->Branch("thisEntry"           , &thisEntry           , "thisEntry/I");
  outTree_->Branch("mypuw"               , &mypuw               , "mypuw/F");
  outTree_->Branch("totWeight"           , &totWeight           , "totWeight/F");
  outTree_->Branch("mcTrue"              , &mcTrue              , "mcTrue/I");

  std::cout << "Booking output histos for event breakdown" << endl;
  h_entries   = new TH1F("h_entries"  , "h_entries"  , 10,   0., 10.);
  h_selection = new TH1F("h_selection", "h_selection",  6, -0.5, 5.5);
  h_entries->Sumw2();
  h_selection->Sumw2();
}

// To compute the pileup weight
float TnPNtuplesTriggerEfficiency::puw2016_nTrueInt_36fb(int nTrueInt) 
{ 
  float _puw2016_nTrueInt_36fb[100] = {0.3505407355600995, 0.8996968628890968, 1.100322319466069, 0.9562526765089195, 1.0366251229154624, 1.0713954619016586, 0.7593488199769544, 0.47490309461978414, 0.7059895997695581, 0.8447022252423783, 0.9169159386164522, 1.0248924033173097, 1.0848877947714115, 1.1350984224561655, 1.1589888429954602, 1.169048420382294, 1.1650383018054549, 1.1507200023444994, 1.1152571438041776, 1.0739529436969637, 1.0458014000030829, 1.032500407707141, 1.0391236062781293, 1.041283620738903, 1.0412963370894526, 1.0558823002770783, 1.073481674823461, 1.0887053272606795, 1.1041701696801014, 1.123218903738397, 1.1157169321377927, 1.1052520327174429, 1.0697489590429388, 1.0144652740600584, 0.9402657069968621, 0.857142825520793, 0.7527112615290031, 0.6420618248685722, 0.5324755829715156, 0.4306470627563325, 0.33289171600176093, 0.24686361729094983, 0.17781595237914027, 0.12404411884835284, 0.08487088505600057, 0.056447805688061216, 0.03540829360547507, 0.022412461576677457, 0.013970541270658443, 0.008587896629717911, 0.004986410514292661, 0.00305102303701641, 0.001832072556146534, 0.0011570757619737708, 0.0008992999249003301, 0.0008241241729452477, 0.0008825716073180279, 0.001187003960081393, 0.0016454104270429153, 0.0022514113879764414, 0.003683196037880878, 0.005456695951503178, 0.006165248770884191, 0.007552675218762607, 0.008525338219226993, 0.008654690499815343, 0.006289068906974821, 0.00652551838513972, 0.005139581024893171, 0.005115751962934923, 0.004182527768384693, 0.004317593022028565, 0.0035749335962533355, 0.003773660372937113, 0.002618732319396435, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0};

  if (nTrueInt<100) return _puw2016_nTrueInt_36fb[nTrueInt]; 
  else return 0; 
}
