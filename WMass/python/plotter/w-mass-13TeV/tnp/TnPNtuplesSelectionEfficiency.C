#ifndef TnPNtuplesSelectionEfficiency_cxx
#define TnPNtuplesSelectionEfficiency_cxx

#include "TnPNtuplesBase.C"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <iostream>
#include <TLorentzVector.h>

class TnPNtuplesSelectionEfficiency : public TnPNtuplesBase {
public:

  TnPNtuplesSelectionEfficiency(TTree *tree=0);
  virtual ~TnPNtuplesSelectionEfficiency();  

  void Loop(int maxentries = -1);
  bool isTagLepton(int);
  bool isLeptonInAcceptance(int jj);
  bool hasTriggerMatch(int jj);
  bool trackCleanedFromMuons(int itrk, int imu);
  bool passesFilters();

protected:

};

#endif

#ifdef TnPNtuplesSelectionEfficiency_cxx
using namespace std;

float deltaPhi(float phi1, float phi2) {                                                        
    float result = phi1 - phi2;
    while (result > float(M_PI)) result -= float(2*M_PI);
    while (result <= -float(M_PI)) result += float(2*M_PI);
    return result;
}

float deltaR2(float eta1, float phi1, float eta2, float phi2) {
    float deta = std::abs(eta1-eta2);
    float dphi = deltaPhi(phi1,phi2);
    return deta*deta + dphi*dphi;
}
float deltaR(float eta1, float phi1, float eta2, float phi2) {
    return std::sqrt(deltaR2(eta1,phi1,eta2,phi2));
}


TnPNtuplesSelectionEfficiency::TnPNtuplesSelectionEfficiency(TTree *tree) : TnPNtuplesBase(tree)
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("/eos/user/m/mdunser/w-helicity-13TeV/trees/TREES_2018-07-06-2l_triggerMatch_MUONS/SingleMuon_Run2016C/treeProducerWMass/tree.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("/eos/user/m/mdunser/w-helicity-13TeV/trees/TREES_2018-07-06-2l_triggerMatch_MUONS/SingleMuon_Run2016C/treeProducerWMass/tree.root");
      }
      f->GetObject("tree",tree);

   }
   //if (ftree) tree->AddFriend(ftree); 
   Init(tree);
}

TnPNtuplesSelectionEfficiency::~TnPNtuplesSelectionEfficiency()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

bool TnPNtuplesSelectionEfficiency::passesFilters(){

  if (!Flag_goodVertices) return false;
  if (!Flag_globalTightHalo2016Filter) return false;
  if (!Flag_HBHENoiseFilter) return false;
  if (!Flag_HBHENoiseIsoFilter) return false;
  if (!Flag_EcalDeadCellTriggerPrimitiveFilter) return false;

  return true;

}

bool TnPNtuplesSelectionEfficiency::isLeptonInAcceptance(int jj) {
  if(fFlavor == 11){
    if (Electron_pt[jj]<25)                                            return false;
    if (fabs(Electron_eta[jj])>1.4442 && fabs(Electron_eta[jj])<1.566) return false;
    if (fabs(Electron_eta[jj])>2.5)                                    return false;
  }
  else {
    if (Muon_pt[jj]<25)            return false;
    if (fabs(Muon_eta[jj])> 2.4)   return false;
  }

  return true;
}

bool TnPNtuplesSelectionEfficiency::trackCleanedFromMuons(int itrk, int imu){
  bool muIsGood = true;

  if (!(Muon_isTracker[imu] || Muon_isGlobal[imu])) muIsGood = false;
  if (!(fabs(Muon_dxy[imu]) < 0.2 && fabs(Muon_dz[imu]) < 0.1 )) muIsGood = false;
  if (!((Muon_pfRelIso03_chg[imu]*Muon_pt[imu] < 5. && Muon_pt[imu] < 25.) || Muon_pfRelIso03_chg[imu] < 0.2 ) ) muIsGood = false;

  float dr = deltaR(Muon_eta[imu], Muon_phi[imu], IsoTrack_eta[itrk], IsoTrack_phi[itrk]);

  if (muIsGood && dr < 0.3) return false;

  return true;
}

    //  triggerBits: 1 = TrkIsoVVL, 2 = Iso, 4 = OverlapFilter PFTau, 8 = IsoTkMu, 1024 = 1mu (Mu50) for Muon
    // marc muon_trigs = [ trig for trig in all_trigs if trig.id==13 and ((trig.filterBits>>0 & 1 ) or (trig.filterBits>>1 
bool TnPNtuplesSelectionEfficiency::hasTriggerMatch(int jj){

    bool hasTrigMatch = false;
    for (unsigned int i=0; i<nTrigObj; ++i){
      if (TrigObj_id[i]  != 13 ) continue;
      if (TrigObj_pt[i]   < 24.) continue;
      if (TrigObj_l1pt[i] < 22.) continue;
      if (! (( TrigObj_filterBits[i] & 8) || (TrigObj_l2pt[i] > 10. && (TrigObj_filterBits[i] & 2) )) ) continue;
      if (deltaR(Muon_eta[jj], Muon_phi[jj], TrigObj_eta[i], TrigObj_phi[i]) < 0.3) {
        hasTrigMatch = true;
        break;
      }
    }
    return hasTrigMatch;
}

bool TnPNtuplesSelectionEfficiency::isTagLepton(int jj){

  if(fFlavor == 11){
    if (Electron_pt[jj]<25)                                        return false;
    if (fabs(Electron_eta[jj])>1.4442 && fabs(Electron_eta[jj])<1.566) return false;
    if (fabs(Electron_eta[jj])>2.5)                                   return false;
    // peu importe pour le moment if (Electron_customId[jj] < 1)                                    return false;
    // peu importe pour le moment if (Electron_hltId[jj] < 1)                                       return false;
    // peu importe pour le moment if (Electron_tightChargeFix[jj] != 2 )                            return false;
  }

  else {
    if (Muon_pt[jj]<25)                 return false;
    if (fabs(Muon_eta[jj])> 2.4)        return false;
    if (Muon_pfRelIso04_all[jj] > 0.15) return false;
    if (fabs(Muon_dxy[jj]) > 0.05)      return false;
    if (fabs(Muon_dz[jj]) > 0.20)       return false;
    if (Muon_mediumId[jj] < 1)          return false;
    if (Muon_isGlobal[jj] < 1)          return false;
    //if (!hasTriggerMatch(jj))           return false;
  }
  return true;

}

void TnPNtuplesSelectionEfficiency::Loop(int maxentries)
{

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
  vector <float> cand_pt           = {};
  vector <float> cand_eta          = {};
  vector <float> cand_dxy          = {};
  vector <float> cand_iso          = {};
  vector <float> cand_iso03        = {};
  vector <float> cand_chgiso03     = {};
  vector <float> cand_dz           = {};
  vector <int>   cand_mediumId     = {};
  vector <float> cand_truept       = {};
  vector <float> cand_trueeta      = {};
  vector <float> cand_etaSc        = {};
  vector <float> cand_phi          = {};
  vector <float> cand_charge       = {};
  vector <int>   cand_triggerMatch = {};
  vector <int>   cand_isGlobal     = {};
  vector <int>   cand_isTracker    = {};
  vector <int>   cand_matchMC      = {};
  vector <int>   cand_customId     = {};
  vector <int>   cand_tightCharge  = {};
  vector <int>   cand_fullLepId    = {};
  vector <int>   cand_alsoTag      = {};
  vector <int>   cand_isMuon       = {};
  vector <int>   cand_isZero       = {};

  // To compute the lumi weight
  float sigma=1976.17;
  float count_getentries = doElectrons ? 123847915 : 99999999999;    // madgraph, ext1+ext2 //MARC THIS NUMBER IS WRONG FOR MUONS. doesn't matter, it's a constant
  float SetLumi=35.92546;     

  // Loop over events
  std::cout << "Start looping over events" << endl;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    if (maxentries > 0 && ientry >= maxentries) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    if (!(ientry%250000)) std::cout << ientry << endl;
    thisEntry = ientry;

    bool isData = (run != 1);

    // To keep track of the total number of events
    h_entries->Fill(5);  

    // PU weight
    puw_inc = 1.;
    puw_bf  = 1.;
    puw_gh  = 1.;
    
    if (!isData) {
        puw_inc = puw_2016(Pileup_nTrueInt, 0);
        puw_bf  = puw_2016(Pileup_nTrueInt, 1);
        puw_gh  = puw_2016(Pileup_nTrueInt, 2);
    }
    
    // Lumi weight for MC only
    totWeight = 1.;
    if (!isData) totWeight = puw_inc*genWeight*((sigma*SetLumi)/count_getentries)*(pow(10,3));  


    // selection for the data JSON
    // ============================================
    if (!isGoodRunLS(isData, run, luminosityBlock)) {
        //std::cout << "does not pass the json " << run << "  " << luminosityBlock << std::endl;
        continue;
    }
    _run = run;
    

    // select some MET filters
    // ==================================
    if (!passesFilters()) continue;

    // Events breakdown  
    h_selection->Fill(0.);

    // 1) analysis cuts: fire the single lepton trigger  
    // marc peu importe pour le moment if ( doElectrons && !HLT_BIT_HLT_Ele27_WPTight_Gsf_v) continue;
    if (!doElectrons && !HLT_IsoMu24 && !HLT_IsoTkMu24) continue; // require the event to pass both vetos

    h_selection->Fill(1.);
    
    // 2) at least one good vertex found
    nvtx = PV_npvsGood;
    if (nvtx<1) continue;
    h_selection->Fill(2.);

    // 3) Gen level match, to be saved as additional info
    int numGenLevel=0;
    bool genLepFound = false;
    bool genPosFound = false;
    TLorentzVector myGenLep(0,0,0,0);  
    TLorentzVector myGenPos(0,0,0,0);  
    if (!isData) {   
      for(unsigned int ii=0; ii<nGenPart; ii++){
        int status = GenPart_status[ii];
        int statusFlag = GenPart_statusFlags[ii];
        int pdgid  = GenPart_pdgId[ii];
        // if ( abs(pdgid)==fFlavor && (status==746 || (status == 1 && getStatusFlag(statusFlag,8) )) ) { // check if lepton is correct type and status 23
         if ( abs(pdgid)==fFlavor && (status==746 || (status == 1 && (GenPart_statusFlags[ii] & 1) )) ) { // check if lepton is correct type and status 23
           //if (GenPart_pdgId[GenPart_genPartIdxMother[ii]]==23){
           float ptgen  = GenPart_pt [ii];
           float etagen = GenPart_eta[ii];
           float phigen = GenPart_phi[ii];
           if (pdgid<0)  {
             myGenPos.SetPtEtaPhiM(ptgen, etagen, phigen, 0.);
             genPosFound = true;
           }
           if (pdgid>0) {
             myGenLep.SetPtEtaPhiM(ptgen, etagen, phigen, 0.);
             genLepFound = true;
           }
         numGenLevel++;
        } // end if pdgId && status
      } // end loop gen particles
    } // end MC only


    // 4) Tag and probe selection
    bool atLeastOneTag = false;
    std::vector<int> acceptLep;
    std::vector<int> acceptIsoTracks;
    
    // leptons in the acceptance
    for(unsigned int jj=0; jj<nMuon; jj++){
      if ( isLeptonInAcceptance(jj) ) acceptLep.push_back(jj);
    }

    // isotracks in the acceptance
    // isolated tracks after basic selection 
    //  ( ((pt>5 && (abs(pdgId) == 11 || abs(pdgId) == 13)) || pt > 10) 
    // && (abs(pdgId) < 15 || abs(eta) < 2.5) 
    // && abs(dxy) < 0.2 
    // && abs(dz) < 0.1 
    // && ((pfIsolationDR03().chargedHadronIso < 5 && pt < 25) || pfIsolationDR03().chargedHadronIso/pt < 0.2)) 
    // and lepton veto*

    for(unsigned int jj=0; jj<nIsoTrack; jj++){
      if ( IsoTrack_pt[jj] < 25. ) continue;
      bool isGoodTrack = true;
      for(unsigned int kk=0; kk<nMuon; kk++){
        if (!trackCleanedFromMuons(jj,kk) ) isGoodTrack = false;
      }
      if (isGoodTrack) acceptIsoTracks.push_back(jj);
    }

    // full selection for tags and probes
    //for (unsigned int iLep=0; iLep<acceptLep.size(); iLep++) {
    for (std::vector<int>::const_iterator ilep = acceptLep.begin(); ilep !=acceptLep.end(); ++ilep){ 
      //int theOrigIndex = acceptLep.at(iLep);
      int theOrigIndex = *ilep;

      // kine 
      float lepPt        = Muon_pt                        [theOrigIndex];     // need calibrated pT eventually!!! FIXME
      float lepEta       = Muon_eta                       [theOrigIndex];
      float lepIso       = Muon_pfRelIso04_all            [theOrigIndex];
      float lepIso03     = Muon_pfRelIso03_all            [theOrigIndex];
      float lepChgIso03  = Muon_pfRelIso03_chg            [theOrigIndex];
      float lepDxy       = Muon_dxy                       [theOrigIndex];
      float lepDz        = Muon_dz                        [theOrigIndex];
      int   lepMediumId  = Muon_mediumId                  [theOrigIndex];
      float lepScEta     = Muon_eta                       [theOrigIndex];
      float lepPhi       = Muon_phi                       [theOrigIndex];
      float lepCharge    = Muon_charge                    [theOrigIndex];
      bool  triggerMatch = hasTriggerMatch(theOrigIndex);
      bool  isGlobal     = Muon_isGlobal                  [theOrigIndex];
      bool  isTracker    = Muon_isTracker                 [theOrigIndex];

      // this lep
      TLorentzVector thisRecoLep(0,0,0,0);
      thisRecoLep.SetPtEtaPhiM(lepPt,lepEta,lepPhi,0.105);

      // Match with MC truth   
      int matchMC = 0;
      float truePt = -1;
      float trueEta = -999;
      if (!isData) {  
        if(genLepFound && thisRecoLep.DeltaR(myGenLep)<0.3) {
          matchMC = 1; truePt = myGenLep.Pt(); trueEta = myGenLep.Eta();
        }
        if(genPosFound && thisRecoLep.DeltaR(myGenPos)<0.3) {
          matchMC = 1; truePt = myGenPos.Pt(); trueEta = myGenPos.Eta();
        }
      } 
      else {
        matchMC = 1;
      } // end check mcMatch

      // marc peu importe pour le moment // HLT Safe ID
      // marc peu importe pour le moment int hltSafeId = LepGood_hltId[theOrigIndex];
      // marc peu importe pour le moment 
      // marc peu importe pour le moment // Full ID 
      // marc peu importe pour le moment int customId = LepGood_customId[theOrigIndex];

      // Tight Charge
      int tightCharge = Muon_tightCharge[theOrigIndex];

      // Is this a tag:
      int isThisTag = 0;
      if (doElectrons){
        // marc peu importe pour le moment if (isTagLepton(theOrigIndex) && LepGood_matchedTrgObjElePt[theOrigIndex] > -1.) isThisTag =1;
        continue;
      }
      else {
        if (isTagLepton(theOrigIndex) && triggerMatch ) isThisTag = 1; // marc fix this trigger match! && (LepGood_matchedTrgObjMuPt[theOrigIndex] > -1. || LepGood_matchedTrgObjTkMuPt[theOrigIndex] > -1.) ) isThisTag =1;
      }
      if (isThisTag) atLeastOneTag = true;

      // Infos to be kept
      cand_pt          . push_back(lepPt);
      cand_eta         . push_back(lepEta);
      cand_dxy         . push_back(lepDxy);
      cand_iso         . push_back(lepIso);
      cand_iso03       . push_back(lepIso03);
      cand_chgiso03    . push_back(lepChgIso03);
      cand_dz          . push_back(lepDz);
      cand_mediumId    . push_back(lepMediumId);
      cand_truept      . push_back(truePt);
      cand_trueeta     . push_back(trueEta);
      cand_etaSc       . push_back(lepScEta);
      cand_phi         . push_back(lepPhi);
      cand_charge      . push_back(lepCharge);
      cand_triggerMatch. push_back(triggerMatch);
      cand_isGlobal    . push_back(isGlobal);
      cand_isTracker   . push_back(isTracker);
      cand_matchMC     . push_back(matchMC);
      cand_tightCharge . push_back(tightCharge);
      cand_fullLepId   . push_back(isThisTag);
      cand_alsoTag     . push_back(isThisTag);
      cand_isMuon      . push_back(1);        

      if (theOrigIndex==0) cand_isZero.push_back(1);
      else cand_isZero.push_back(0);    
      
    } // end leptons in acceptance

    // now loop on all the iso tracks in acceptance
    for (std::vector<int>::const_iterator itrk = acceptIsoTracks.begin(); itrk !=acceptIsoTracks.end(); ++itrk){ 
      //int theOrigIndex = acceptLep.at(iLep);
      int theOrigIndex = *itrk;

      // kine 
      float trkPt        = IsoTrack_pt                        [theOrigIndex];     // need calibrated pT eventually!!! FIXME
      float trkEta       = IsoTrack_eta                       [theOrigIndex];
      float trkPhi       = IsoTrack_phi                       [theOrigIndex];
      float trkIso       = IsoTrack_pfRelIso03_all            [theOrigIndex];
      float trkChgIso    = IsoTrack_pfRelIso03_chg            [theOrigIndex];
      float trkDxy       = IsoTrack_dxy                       [theOrigIndex];
      float trkDz        = IsoTrack_dz                        [theOrigIndex];
      bool  isMuon       = 0   ;
      // Infos to be kept
      cand_pt          . push_back(trkPt);
      cand_eta         . push_back(trkEta);
      cand_dxy         . push_back(trkDxy);
      cand_iso         . push_back(trkIso);
      cand_iso03       . push_back(trkIso);
      cand_chgiso03    . push_back(trkChgIso);
      cand_dz          . push_back(trkDz);
      cand_phi         . push_back(trkPhi);
      cand_mediumId    . push_back(-999);
      cand_truept      . push_back(-999);
      cand_trueeta     . push_back(-999);
      cand_etaSc       . push_back(-999);
      cand_charge      . push_back(-999);
      cand_triggerMatch. push_back(-999);
      cand_isGlobal    . push_back(-999);
      cand_isTracker   . push_back(-999);
      cand_matchMC     . push_back(-999);
      cand_tightCharge . push_back(-999);
      cand_fullLepId   . push_back(-999);
      cand_alsoTag     . push_back(0);
      cand_isMuon      . push_back(0);        

      if (theOrigIndex==0) cand_isZero.push_back(1);
      else cand_isZero.push_back(0);    
    }
    // end loop on isotrack candidates

    
    //--- 4)
    if (!atLeastOneTag) {
      // cleaning vectors
      cand_pt          . clear();
      cand_eta         . clear();
      cand_dxy         . clear();
      cand_iso         . clear();
      cand_iso03       . clear();
      cand_chgiso03    . clear();
      cand_dz          . clear();
      cand_mediumId    . clear();
      cand_truept      . clear();
      cand_trueeta     . clear();
      cand_etaSc       . clear();
      cand_phi         . clear();
      cand_charge      . clear();
      cand_triggerMatch. clear();
      cand_isGlobal    . clear();
      cand_isTracker   . clear();
      cand_matchMC     . clear();
      cand_tightCharge . clear();
      cand_fullLepId   . clear();
      cand_alsoTag     . clear();
      cand_isZero      . clear();
      cand_isMuon      . clear();
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
        // if (cand_isZero[iLep2]) continue; # consider 2 candidates as possible probes with same minv, otherwise biases minv for high pt probe
        TLorentzVector thisLep2(0,0,0,0);
        thisLep2.SetPtEtaPhiM(cand_pt[iLep2],cand_eta[iLep2],cand_phi[iLep2],0);

        // invariant mass
        pair_mass = (thisLep1+thisLep2).M();
        if (pair_mass<40 || pair_mass>140) continue;

        // both matching mc truth?
        mcTrue = cand_matchMC[iLep1] && cand_matchMC[iLep2];
        
        // first as tag, second as probe
        tag_pt             = cand_pt          [iLep1];
        tag_eta            = cand_eta         [iLep1];
        tag_matchMC        = cand_matchMC     [iLep1];

        probe_pt           = cand_pt          [iLep2];
        probe_eta          = cand_eta         [iLep2];
        probe_dxy          = cand_dxy         [iLep2];
        probe_iso          = cand_iso         [iLep2];
        probe_iso03        = cand_iso03       [iLep2];
        probe_chgiso03     = cand_chgiso03    [iLep2];
        probe_dz           = cand_dz          [iLep2];
        probe_mediumId     = cand_mediumId    [iLep2];
        probe_truept       = cand_truept      [iLep2];
        probe_trueeta      = cand_trueeta     [iLep2];
        probe_sc_eta       = cand_etaSc       [iLep2];
        probe_phi          = cand_phi         [iLep2];
        probe_charge       = cand_charge      [iLep2];
        probe_triggerMatch = cand_triggerMatch[iLep2];
        probe_isGlobal     = cand_isGlobal    [iLep2];
        probe_isTracker    = cand_isTracker   [iLep2];
        probe_matchMC      = cand_matchMC     [iLep2];
        probe_tightCharge  = cand_tightCharge [iLep2];
        probe_fullLepId    = cand_fullLepId   [iLep2];
        probe_alsoTag      = cand_alsoTag     [iLep2];
        probe_isMuon       = cand_isMuon      [iLep2];
        
        // Tree filling
        outFile_->cd();
        cddir->cd();  
        outTree_->Fill();

      }  // probes
    }   // tags
  
    // cleaning vectors
    cand_pt           .clear();
    cand_eta          .clear();
    cand_dxy          .clear();
    cand_iso          .clear();
    cand_iso03        .clear();
    cand_chgiso03     .clear();
    cand_dz           .clear();
    cand_mediumId     .clear();
    cand_truept       .clear();
    cand_trueeta      .clear();
    cand_etaSc        .clear();
    cand_phi          .clear();
    cand_charge       .clear();
    cand_triggerMatch .clear();
    cand_isGlobal     .clear();
    cand_isTracker    .clear();
    cand_matchMC      .clear();
    cand_customId     .clear();
    cand_tightCharge  .clear();
    cand_fullLepId    .clear();
    cand_alsoTag      .clear();
    cand_isMuon       .clear();
    cand_isZero       .clear();

  }  // Loop over entries

  // Saving output tree and histos
  outFile_    -> cd();
  h_entries   -> Write();
  h_selection -> Write();
  cddir       -> cd();
  outTree_    -> Write();
  outFile_    -> Close();

} // Loop method

#endif
