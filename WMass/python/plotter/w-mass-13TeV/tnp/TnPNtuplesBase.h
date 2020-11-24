//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue Oct 27 16:10:48 2020 by ROOT version 6.10/09
// from TTree Events/Events
// found on file: /eos/cms/store/cmst3/group/wmass/w-mass-13TeV/NanoAOD/DYJetsToMuMu_M-50_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoAODv7/201025_173845/0000/SMP-RunIISummer16NanoAODv7-00336_123.root
//////////////////////////////////////////////////////////

#ifndef TnPNtuplesBase_h
#define TnPNtuplesBase_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TString.h>
#include <TH1F.h>

// Header file for the classes stored in the TTree if any.

class TnPNtuplesBase {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

  // Output tree and histos
  TH1F *h_entries;
  TH1F *h_selection;

  TFile *outFile_;
  TTree* outTree_;
  TDirectory *cddir;

  float tag_lep_pt, tag_lep_eta;
  int tag_lep_matchMC;
  float probe_lep_pt, probe_lep_eta, probe_sc_eta, probe_lep_phi, probe_lep_charge, probe_lep_truept, probe_lep_trueeta;
  float probe_eleTrgPt, probe_muTrgPt, probe_tkMuTrgPt;
  int probe_lep_matchMC;
  int probe_lep_tightCharge, probe_lep_fullLepId, probe_lep_pdgId, probe_triggerMatch;
  int probe_lep_alsoTag;
  float pair_mass;
  int nvtx;
  int thisEntry;
  float mypuw, totWeight;
  int mcTrue;


   // Declaration of leaf types
   UInt_t          run;
   UInt_t          luminosityBlock;
   ULong64_t       event;
   UInt_t          nElectron;
   Float_t         Electron_deltaEtaSC[5];   //[nElectron]
   Float_t         Electron_dr03EcalRecHitSumEt[5];   //[nElectron]
   Float_t         Electron_dr03HcalDepth1TowerSumEt[5];   //[nElectron]
   Float_t         Electron_dr03TkSumPt[5];   //[nElectron]
   Float_t         Electron_dr03TkSumPtHEEP[5];   //[nElectron]
   Float_t         Electron_dxy[5];   //[nElectron]
   Float_t         Electron_dxyErr[5];   //[nElectron]
   Float_t         Electron_dz[5];   //[nElectron]
   Float_t         Electron_dzErr[5];   //[nElectron]
   Float_t         Electron_eInvMinusPInv[5];   //[nElectron]
   Float_t         Electron_energyErr[5];   //[nElectron]
   Float_t         Electron_eta[5];   //[nElectron]
   Float_t         Electron_hoe[5];   //[nElectron]
   Float_t         Electron_ip3d[5];   //[nElectron]
   Float_t         Electron_jetPtRelv2[5];   //[nElectron]
   Float_t         Electron_jetRelIso[5];   //[nElectron]
   Float_t         Electron_mass[5];   //[nElectron]
   Float_t         Electron_miniPFRelIso_all[5];   //[nElectron]
   Float_t         Electron_miniPFRelIso_chg[5];   //[nElectron]
   Float_t         Electron_mvaFall17V1Iso[5];   //[nElectron]
   Float_t         Electron_mvaFall17V1noIso[5];   //[nElectron]
   Float_t         Electron_mvaFall17V2Iso[5];   //[nElectron]
   Float_t         Electron_mvaFall17V2noIso[5];   //[nElectron]
   Float_t         Electron_pfRelIso03_all[5];   //[nElectron]
   Float_t         Electron_pfRelIso03_chg[5];   //[nElectron]
   Float_t         Electron_phi[5];   //[nElectron]
   Float_t         Electron_pt[5];   //[nElectron]
   Float_t         Electron_r9[5];   //[nElectron]
   Float_t         Electron_scEtOverPt[5];   //[nElectron]
   Float_t         Electron_sieie[5];   //[nElectron]
   Float_t         Electron_sip3d[5];   //[nElectron]
   Float_t         Electron_mvaTTH[5];   //[nElectron]
   Int_t           Electron_charge[5];   //[nElectron]
   Int_t           Electron_cutBased[5];   //[nElectron]
   Int_t           Electron_cutBased_Fall17_V1[5];   //[nElectron]
   Int_t           Electron_jetIdx[5];   //[nElectron]
   Int_t           Electron_pdgId[5];   //[nElectron]
   Int_t           Electron_photonIdx[5];   //[nElectron]
   Int_t           Electron_tightCharge[5];   //[nElectron]
   Int_t           Electron_vidNestedWPBitmap[5];   //[nElectron]
   Int_t           Electron_vidNestedWPBitmapHEEP[5];   //[nElectron]
   Bool_t          Electron_convVeto[5];   //[nElectron]
   Bool_t          Electron_cutBased_HEEP[5];   //[nElectron]
   Bool_t          Electron_isPFcand[5];   //[nElectron]
   UChar_t         Electron_lostHits[5];   //[nElectron]
   Bool_t          Electron_mvaFall17V1Iso_WP80[5];   //[nElectron]
   Bool_t          Electron_mvaFall17V1Iso_WP90[5];   //[nElectron]
   Bool_t          Electron_mvaFall17V1Iso_WPL[5];   //[nElectron]
   Bool_t          Electron_mvaFall17V1noIso_WP80[5];   //[nElectron]
   Bool_t          Electron_mvaFall17V1noIso_WP90[5];   //[nElectron]
   Bool_t          Electron_mvaFall17V1noIso_WPL[5];   //[nElectron]
   Bool_t          Electron_mvaFall17V2Iso_WP80[5];   //[nElectron]
   Bool_t          Electron_mvaFall17V2Iso_WP90[5];   //[nElectron]
   Bool_t          Electron_mvaFall17V2Iso_WPL[5];   //[nElectron]
   Bool_t          Electron_mvaFall17V2noIso_WP80[5];   //[nElectron]
   Bool_t          Electron_mvaFall17V2noIso_WP90[5];   //[nElectron]
   Bool_t          Electron_mvaFall17V2noIso_WPL[5];   //[nElectron]
   UChar_t         Electron_seedGain[5];   //[nElectron]
   UInt_t          nFsrPhoton;
   Float_t         FsrPhoton_dROverEt2[3];   //[nFsrPhoton]
   Float_t         FsrPhoton_eta[3];   //[nFsrPhoton]
   Float_t         FsrPhoton_phi[3];   //[nFsrPhoton]
   Float_t         FsrPhoton_pt[3];   //[nFsrPhoton]
   Float_t         FsrPhoton_relIso03[3];   //[nFsrPhoton]
   Int_t           FsrPhoton_muonIdx[3];   //[nFsrPhoton]
   UInt_t          nGenPart;
   Float_t         GenPart_eta[111];   //[nGenPart]
   Float_t         GenPart_mass[111];   //[nGenPart]
   Float_t         GenPart_phi[111];   //[nGenPart]
   Float_t         GenPart_pt[111];   //[nGenPart]
   Int_t           GenPart_genPartIdxMother[111];   //[nGenPart]
   Int_t           GenPart_pdgId[111];   //[nGenPart]
   Int_t           GenPart_status[111];   //[nGenPart]
   Int_t           GenPart_statusFlags[111];   //[nGenPart]
   Float_t         Generator_binvar;
   Float_t         Generator_scalePDF;
   Float_t         Generator_weight;
   Float_t         Generator_x1;
   Float_t         Generator_x2;
   Float_t         Generator_xpdf1;
   Float_t         Generator_xpdf2;
   Int_t           Generator_id1;
   Int_t           Generator_id2;
   Float_t         genWeight;
   Float_t         L1PreFiringWeight_Dn;
   Float_t         L1PreFiringWeight_Nom;
   Float_t         L1PreFiringWeight_Up;
   UInt_t          nMuon;
   Float_t         Muon_dxy[10];   //[nMuon]
   Float_t         Muon_dxyErr[10];   //[nMuon]
   Float_t         Muon_dxybs[10];   //[nMuon]
   Float_t         Muon_dz[10];   //[nMuon]
   Float_t         Muon_dzErr[10];   //[nMuon]
   Float_t         Muon_eta[10];   //[nMuon]
   Float_t         Muon_ip3d[10];   //[nMuon]
   Float_t         Muon_jetPtRelv2[10];   //[nMuon]
   Float_t         Muon_jetRelIso[10];   //[nMuon]
   Float_t         Muon_mass[10];   //[nMuon]
   Float_t         Muon_miniPFRelIso_all[10];   //[nMuon]
   Float_t         Muon_miniPFRelIso_chg[10];   //[nMuon]
   Float_t         Muon_pfRelIso03_all[10];   //[nMuon]
   Float_t         Muon_pfRelIso03_chg[10];   //[nMuon]
   Float_t         Muon_pfRelIso04_all[10];   //[nMuon]
   Float_t         Muon_phi[10];   //[nMuon]
   Float_t         Muon_pt[10];   //[nMuon]
   Float_t         Muon_ptErr[10];   //[nMuon]
   Float_t         Muon_segmentComp[10];   //[nMuon]
   Float_t         Muon_sip3d[10];   //[nMuon]
   Float_t         Muon_softMva[10];   //[nMuon]
   Float_t         Muon_tkRelIso[10];   //[nMuon]
   Float_t         Muon_tunepRelPt[10];   //[nMuon]
   Float_t         Muon_mvaLowPt[10];   //[nMuon]
   Float_t         Muon_mvaTTH[10];   //[nMuon]
   Int_t           Muon_charge[10];   //[nMuon]
   Int_t           Muon_jetIdx[10];   //[nMuon]
   Int_t           Muon_nStations[10];   //[nMuon]
   Int_t           Muon_nTrackerLayers[10];   //[nMuon]
   Int_t           Muon_pdgId[10];   //[nMuon]
   Int_t           Muon_tightCharge[10];   //[nMuon]
   Int_t           Muon_fsrPhotonIdx[10];   //[nMuon]
   UChar_t         Muon_highPtId[10];   //[nMuon]
   Bool_t          Muon_highPurity[10];   //[nMuon]
   Bool_t          Muon_inTimeMuon[10];   //[nMuon]
   Bool_t          Muon_isGlobal[10];   //[nMuon]
   Bool_t          Muon_isPFcand[10];   //[nMuon]
   Bool_t          Muon_isTracker[10];   //[nMuon]
   Bool_t          Muon_looseId[10];   //[nMuon]
   Bool_t          Muon_mediumId[10];   //[nMuon]
   Bool_t          Muon_mediumPromptId[10];   //[nMuon]
   UChar_t         Muon_miniIsoId[10];   //[nMuon]
   UChar_t         Muon_multiIsoId[10];   //[nMuon]
   UChar_t         Muon_mvaId[10];   //[nMuon]
   UChar_t         Muon_mvaLowPtId[10];   //[nMuon]
   UChar_t         Muon_pfIsoId[10];   //[nMuon]
   UChar_t         Muon_puppiIsoId[10];   //[nMuon]
   Bool_t          Muon_softId[10];   //[nMuon]
   Bool_t          Muon_softMvaId[10];   //[nMuon]
   Bool_t          Muon_tightId[10];   //[nMuon]
   UChar_t         Muon_tkIsoId[10];   //[nMuon]
   Bool_t          Muon_triggerIdLoose[10];   //[nMuon]
   Float_t         Pileup_nTrueInt;
   UInt_t          nGenDressedLepton;
   Float_t         GenDressedLepton_eta[4];   //[nGenDressedLepton]
   Float_t         GenDressedLepton_mass[4];   //[nGenDressedLepton]
   Float_t         GenDressedLepton_phi[4];   //[nGenDressedLepton]
   Float_t         GenDressedLepton_pt[4];   //[nGenDressedLepton]
   Int_t           GenDressedLepton_pdgId[4];   //[nGenDressedLepton]
   Bool_t          GenDressedLepton_hasTauAnc[4];   //[nGenDressedLepton]
   UInt_t          nTrigObj;
   Float_t         TrigObj_pt[60];   //[nTrigObj]
   Float_t         TrigObj_eta[60];   //[nTrigObj]
   Float_t         TrigObj_phi[60];   //[nTrigObj]
   Float_t         TrigObj_l1pt[60];   //[nTrigObj]
   Float_t         TrigObj_l1pt_2[60];   //[nTrigObj]
   Float_t         TrigObj_l2pt[60];   //[nTrigObj]
   Int_t           TrigObj_id[60];   //[nTrigObj]
   Int_t           TrigObj_l1iso[60];   //[nTrigObj]
   Int_t           TrigObj_l1charge[60];   //[nTrigObj]
   Int_t           TrigObj_filterBits[60];   //[nTrigObj]
   UInt_t          nOtherPV;
   Float_t         OtherPV_z[3];   //[nOtherPV]
   Float_t         PV_ndof;
   Float_t         PV_x;
   Float_t         PV_y;
   Float_t         PV_z;
   Float_t         PV_chi2;
   Float_t         PV_score;
   Int_t           PV_npvs;
   Int_t           PV_npvsGood;
   Int_t           Electron_genPartIdx[5];   //[nElectron]
   UChar_t         Electron_genPartFlav[5];   //[nElectron]
   Int_t           Muon_genPartIdx[6];   //[nMuon]
   UChar_t         Muon_genPartFlav[6];   //[nMuon]
   UChar_t         Electron_cleanmask[5];   //[nElectron]
   UChar_t         Muon_cleanmask[6];   //[nMuon]
   Bool_t          Flag_HBHENoiseFilter;
   Bool_t          Flag_HBHENoiseIsoFilter;
   Bool_t          Flag_CSCTightHaloFilter;
   Bool_t          Flag_CSCTightHaloTrkMuUnvetoFilter;
   Bool_t          Flag_CSCTightHalo2015Filter;
   Bool_t          Flag_globalTightHalo2016Filter;
   Bool_t          Flag_globalSuperTightHalo2016Filter;
   Bool_t          Flag_HcalStripHaloFilter;
   Bool_t          Flag_hcalLaserEventFilter;
   Bool_t          Flag_EcalDeadCellTriggerPrimitiveFilter;
   Bool_t          Flag_EcalDeadCellBoundaryEnergyFilter;
   Bool_t          Flag_ecalBadCalibFilter;
   Bool_t          Flag_goodVertices;
   Bool_t          Flag_eeBadScFilter;
   Bool_t          Flag_ecalLaserCorrFilter;
   Bool_t          Flag_trkPOGFilters;
   Bool_t          Flag_chargedHadronTrackResolutionFilter;
   Bool_t          Flag_muonBadTrackFilter;
   Bool_t          Flag_BadChargedCandidateFilter;
   Bool_t          Flag_BadPFMuonFilter;
   Bool_t          Flag_BadChargedCandidateSummer16Filter;
   Bool_t          Flag_BadPFMuonSummer16Filter;
   Bool_t          Flag_trkPOG_manystripclus53X;
   Bool_t          Flag_trkPOG_toomanystripclus53X;
   Bool_t          Flag_trkPOG_logErrorTooManyClusters;
   Bool_t          Flag_METFilters;
   Bool_t          HLT_IsoMu24;
   Bool_t          HLT_IsoTkMu24;

   // List of branches
   TBranch    *b_run;
   TBranch    *b_luminosityBlock;
   TBranch    *b_event;
   TBranch    *b_nElectron;
   TBranch    *b_Electron_deltaEtaSC;
   TBranch    *b_Electron_dr03EcalRecHitSumEt;
   TBranch    *b_Electron_dr03HcalDepth1TowerSumEt;
   TBranch    *b_Electron_dr03TkSumPt;
   TBranch    *b_Electron_dr03TkSumPtHEEP;
   TBranch    *b_Electron_dxy;
   TBranch    *b_Electron_dxyErr;
   TBranch    *b_Electron_dz;
   TBranch    *b_Electron_dzErr;
   TBranch    *b_Electron_eInvMinusPInv;
   TBranch    *b_Electron_energyErr;
   TBranch    *b_Electron_eta;
   TBranch    *b_Electron_hoe;
   TBranch    *b_Electron_ip3d;
   TBranch    *b_Electron_jetPtRelv2;
   TBranch    *b_Electron_jetRelIso;
   TBranch    *b_Electron_mass;
   TBranch    *b_Electron_miniPFRelIso_all;
   TBranch    *b_Electron_miniPFRelIso_chg;
   TBranch    *b_Electron_mvaFall17V1Iso;
   TBranch    *b_Electron_mvaFall17V1noIso;
   TBranch    *b_Electron_mvaFall17V2Iso;
   TBranch    *b_Electron_mvaFall17V2noIso;
   TBranch    *b_Electron_pfRelIso03_all;
   TBranch    *b_Electron_pfRelIso03_chg;
   TBranch    *b_Electron_phi;
   TBranch    *b_Electron_pt;
   TBranch    *b_Electron_r9;
   TBranch    *b_Electron_scEtOverPt;
   TBranch    *b_Electron_sieie;
   TBranch    *b_Electron_sip3d;
   TBranch    *b_Electron_mvaTTH;
   TBranch    *b_Electron_charge;
   TBranch    *b_Electron_cutBased;
   TBranch    *b_Electron_cutBased_Fall17_V1;
   TBranch    *b_Electron_jetIdx;
   TBranch    *b_Electron_pdgId;
   TBranch    *b_Electron_photonIdx;
   TBranch    *b_Electron_tightCharge;
   TBranch    *b_Electron_vidNestedWPBitmap;
   TBranch    *b_Electron_vidNestedWPBitmapHEEP;
   TBranch    *b_Electron_convVeto;
   TBranch    *b_Electron_cutBased_HEEP;
   TBranch    *b_Electron_isPFcand;
   TBranch    *b_Electron_lostHits;
   TBranch    *b_Electron_mvaFall17V1Iso_WP80;
   TBranch    *b_Electron_mvaFall17V1Iso_WP90;
   TBranch    *b_Electron_mvaFall17V1Iso_WPL;
   TBranch    *b_Electron_mvaFall17V1noIso_WP80;
   TBranch    *b_Electron_mvaFall17V1noIso_WP90;
   TBranch    *b_Electron_mvaFall17V1noIso_WPL;
   TBranch    *b_Electron_mvaFall17V2Iso_WP80;
   TBranch    *b_Electron_mvaFall17V2Iso_WP90;
   TBranch    *b_Electron_mvaFall17V2Iso_WPL;
   TBranch    *b_Electron_mvaFall17V2noIso_WP80;
   TBranch    *b_Electron_mvaFall17V2noIso_WP90;
   TBranch    *b_Electron_mvaFall17V2noIso_WPL;
   TBranch    *b_Electron_seedGain;
   TBranch    *b_nFsrPhoton;
   TBranch    *b_FsrPhoton_dROverEt2;
   TBranch    *b_FsrPhoton_eta;
   TBranch    *b_FsrPhoton_phi;
   TBranch    *b_FsrPhoton_pt;
   TBranch    *b_FsrPhoton_relIso03;
   TBranch    *b_FsrPhoton_muonIdx;
   TBranch    *b_nGenPart;
   TBranch    *b_GenPart_eta;
   TBranch    *b_GenPart_mass;
   TBranch    *b_GenPart_phi;
   TBranch    *b_GenPart_pt;
   TBranch    *b_GenPart_genPartIdxMother;
   TBranch    *b_GenPart_pdgId;
   TBranch    *b_GenPart_status;
   TBranch    *b_GenPart_statusFlags;
   TBranch    *b_Generator_binvar;
   TBranch    *b_Generator_scalePDF;
   TBranch    *b_Generator_weight;
   TBranch    *b_Generator_x1;
   TBranch    *b_Generator_x2;
   TBranch    *b_Generator_xpdf1;
   TBranch    *b_Generator_xpdf2;
   TBranch    *b_Generator_id1;
   TBranch    *b_Generator_id2;
   TBranch    *b_genWeight;
   TBranch    *b_L1PreFiringWeight_Dn;
   TBranch    *b_L1PreFiringWeight_Nom;
   TBranch    *b_L1PreFiringWeight_Up;
   TBranch    *b_nMuon;
   TBranch    *b_Muon_dxy;
   TBranch    *b_Muon_dxyErr;
   TBranch    *b_Muon_dxybs;
   TBranch    *b_Muon_dz;
   TBranch    *b_Muon_dzErr;
   TBranch    *b_Muon_eta;
   TBranch    *b_Muon_ip3d;
   TBranch    *b_Muon_jetPtRelv2;
   TBranch    *b_Muon_jetRelIso;
   TBranch    *b_Muon_mass;
   TBranch    *b_Muon_miniPFRelIso_all;
   TBranch    *b_Muon_miniPFRelIso_chg;
   TBranch    *b_Muon_pfRelIso03_all;
   TBranch    *b_Muon_pfRelIso03_chg;
   TBranch    *b_Muon_pfRelIso04_all;
   TBranch    *b_Muon_phi;
   TBranch    *b_Muon_pt;
   TBranch    *b_Muon_ptErr;
   TBranch    *b_Muon_segmentComp;
   TBranch    *b_Muon_sip3d;
   TBranch    *b_Muon_softMva;
   TBranch    *b_Muon_tkRelIso;
   TBranch    *b_Muon_tunepRelPt;
   TBranch    *b_Muon_mvaLowPt;
   TBranch    *b_Muon_mvaTTH;
   TBranch    *b_Muon_charge;
   TBranch    *b_Muon_jetIdx;
   TBranch    *b_Muon_nStations;
   TBranch    *b_Muon_nTrackerLayers;
   TBranch    *b_Muon_pdgId;
   TBranch    *b_Muon_tightCharge;
   TBranch    *b_Muon_fsrPhotonIdx;
   TBranch    *b_Muon_highPtId;
   TBranch    *b_Muon_highPurity;
   TBranch    *b_Muon_inTimeMuon;
   TBranch    *b_Muon_isGlobal;
   TBranch    *b_Muon_isPFcand;
   TBranch    *b_Muon_isTracker;
   TBranch    *b_Muon_looseId;
   TBranch    *b_Muon_mediumId;
   TBranch    *b_Muon_mediumPromptId;
   TBranch    *b_Muon_miniIsoId;
   TBranch    *b_Muon_multiIsoId;
   TBranch    *b_Muon_mvaId;
   TBranch    *b_Muon_mvaLowPtId;
   TBranch    *b_Muon_pfIsoId;
   TBranch    *b_Muon_puppiIsoId;
   TBranch    *b_Muon_softId;
   TBranch    *b_Muon_softMvaId;
   TBranch    *b_Muon_tightId;
   TBranch    *b_Muon_tkIsoId;
   TBranch    *b_Muon_triggerIdLoose;
   TBranch    *b_Pileup_nTrueInt;
   TBranch    *b_nGenDressedLepton;
   TBranch    *b_GenDressedLepton_eta;
   TBranch    *b_GenDressedLepton_mass;
   TBranch    *b_GenDressedLepton_phi;
   TBranch    *b_GenDressedLepton_pt;
   TBranch    *b_GenDressedLepton_pdgId;
   TBranch    *b_GenDressedLepton_hasTauAnc;
   TBranch    *b_nTrigObj;
   TBranch    *b_TrigObj_pt;
   TBranch    *b_TrigObj_eta;
   TBranch    *b_TrigObj_phi;
   TBranch    *b_TrigObj_l1pt;
   TBranch    *b_TrigObj_l1pt_2;
   TBranch    *b_TrigObj_l2pt;
   TBranch    *b_TrigObj_id;
   TBranch    *b_TrigObj_l1iso;
   TBranch    *b_TrigObj_l1charge;
   TBranch    *b_TrigObj_filterBits;
   TBranch    *b_nOtherPV;
   TBranch    *b_OtherPV_z;
   TBranch    *b_PV_ndof;
   TBranch    *b_PV_x;
   TBranch    *b_PV_y;
   TBranch    *b_PV_z;
   TBranch    *b_PV_chi2;
   TBranch    *b_PV_score;
   TBranch    *b_PV_npvs;
   TBranch    *b_PV_npvsGood;
   TBranch    *b_Electron_genPartIdx;
   TBranch    *b_Electron_genPartFlav;
   TBranch    *b_Muon_genPartIdx;
   TBranch    *b_Muon_genPartFlav;
   TBranch    *b_Electron_cleanmask;
   TBranch    *b_Muon_cleanmask;
   TBranch    *b_Flag_HBHENoiseFilter;
   TBranch    *b_Flag_HBHENoiseIsoFilter;
   TBranch    *b_Flag_CSCTightHaloFilter;
   TBranch    *b_Flag_CSCTightHaloTrkMuUnvetoFilter;
   TBranch    *b_Flag_CSCTightHalo2015Filter;
   TBranch    *b_Flag_globalTightHalo2016Filter;
   TBranch    *b_Flag_globalSuperTightHalo2016Filter;
   TBranch    *b_Flag_HcalStripHaloFilter;
   TBranch    *b_Flag_hcalLaserEventFilter;
   TBranch    *b_Flag_EcalDeadCellTriggerPrimitiveFilter;
   TBranch    *b_Flag_EcalDeadCellBoundaryEnergyFilter;
   TBranch    *b_Flag_ecalBadCalibFilter;
   TBranch    *b_Flag_goodVertices;
   TBranch    *b_Flag_eeBadScFilter;
   TBranch    *b_Flag_ecalLaserCorrFilter;
   TBranch    *b_Flag_trkPOGFilters;
   TBranch    *b_Flag_chargedHadronTrackResolutionFilter;
   TBranch    *b_Flag_muonBadTrackFilter;
   TBranch    *b_Flag_BadChargedCandidateFilter;
   TBranch    *b_Flag_BadPFMuonFilter;
   TBranch    *b_Flag_BadChargedCandidateSummer16Filter;
   TBranch    *b_Flag_BadPFMuonSummer16Filter;
   TBranch    *b_Flag_trkPOG_manystripclus53X;
   TBranch    *b_Flag_trkPOG_toomanystripclus53X;
   TBranch    *b_Flag_trkPOG_logErrorTooManyClusters;
   TBranch    *b_Flag_METFilters;
   TBranch    *b_HLT_IsoMu24;
   TBranch    *b_HLT_IsoTkMu24;

   TnPNtuplesBase(TTree *tree=0);
   virtual ~TnPNtuplesBase();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop(int maxentries = -1);
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);

  virtual void     setFlavor(int);
  virtual void     setOutfile(TString);
  virtual bool     isTagLepton(int);
  virtual bool     getStatusFlag(Int_t, int);
  virtual void makeMapFromJson(const string myJsonFile) ;
  virtual Bool_t isGoodRunLS(Bool_t isData,UInt_t run, UInt_t lumis);

  std::map< UInt_t, std::vector< std::pair<UInt_t,UInt_t> > > fJsonMap;

  int fFlavor;
  TString fOutfile;

  void bookOutputTree();

  float puw2016_nTrueInt_36fb(int nTrueInt);
};

#endif

#ifdef TnPNtuplesBase_cxx
TnPNtuplesBase::TnPNtuplesBase(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("/eos/cms/store/cmst3/group/wmass/w-mass-13TeV/NanoAOD/DYJetsToMuMu_M-50_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoAODv7/201025_173845/0000/SMP-RunIISummer16NanoAODv7-00336_123.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("/eos/cms/store/cmst3/group/wmass/w-mass-13TeV/NanoAOD/DYJetsToMuMu_M-50_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoAODv7/201025_173845/0000/SMP-RunIISummer16NanoAODv7-00336_123.root");
      }
      f->GetObject("Events",tree);

   }
   // if (ftree) tree->AddFriend(ftree); 
   Init(tree);
}

TnPNtuplesBase::~TnPNtuplesBase()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t TnPNtuplesBase::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t TnPNtuplesBase::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void TnPNtuplesBase::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("run", &run, &b_run);
   fChain->SetBranchAddress("luminosityBlock",                          &luminosityBlock                        , &b_luminosityBlock);
   fChain->SetBranchAddress("event",                                    &event                                  , &b_event);
   fChain->SetBranchAddress("nElectron",                                &nElectron                              , &b_nElectron);
   fChain->SetBranchAddress("Electron_deltaEtaSC",                      &Electron_deltaEtaSC                    , &b_Electron_deltaEtaSC);
   fChain->SetBranchAddress("Electron_dr03EcalRecHitSumEt",             &Electron_dr03EcalRecHitSumEt           , &b_Electron_dr03EcalRecHitSumEt);
   fChain->SetBranchAddress("Electron_dr03HcalDepth1TowerSumEt",        &Electron_dr03HcalDepth1TowerSumEt      , &b_Electron_dr03HcalDepth1TowerSumEt);
   fChain->SetBranchAddress("Electron_dr03TkSumPt",                     &Electron_dr03TkSumPt                   , &b_Electron_dr03TkSumPt);
   fChain->SetBranchAddress("Electron_dr03TkSumPtHEEP",                 &Electron_dr03TkSumPtHEEP               , &b_Electron_dr03TkSumPtHEEP);
   fChain->SetBranchAddress("Electron_dxy",                             &Electron_dxy                           , &b_Electron_dxy);
   fChain->SetBranchAddress("Electron_dxyErr",                          &Electron_dxyErr                        , &b_Electron_dxyErr);
   fChain->SetBranchAddress("Electron_dz",                              &Electron_dz                            , &b_Electron_dz);
   fChain->SetBranchAddress("Electron_dzErr",                           &Electron_dzErr                         , &b_Electron_dzErr);
   fChain->SetBranchAddress("Electron_eInvMinusPInv",                   &Electron_eInvMinusPInv                 , &b_Electron_eInvMinusPInv);
   fChain->SetBranchAddress("Electron_energyErr",                       &Electron_energyErr                     , &b_Electron_energyErr);
   fChain->SetBranchAddress("Electron_eta",                             &Electron_eta                           , &b_Electron_eta);
   fChain->SetBranchAddress("Electron_hoe",                             &Electron_hoe                           , &b_Electron_hoe);
   fChain->SetBranchAddress("Electron_ip3d",                            &Electron_ip3d                          , &b_Electron_ip3d);
   fChain->SetBranchAddress("Electron_jetPtRelv2",                      &Electron_jetPtRelv2                    , &b_Electron_jetPtRelv2);
   fChain->SetBranchAddress("Electron_jetRelIso",                       &Electron_jetRelIso                     , &b_Electron_jetRelIso);
   fChain->SetBranchAddress("Electron_mass",                            &Electron_mass                          , &b_Electron_mass);
   fChain->SetBranchAddress("Electron_miniPFRelIso_all",                &Electron_miniPFRelIso_all              , &b_Electron_miniPFRelIso_all);
   fChain->SetBranchAddress("Electron_miniPFRelIso_chg",                &Electron_miniPFRelIso_chg              , &b_Electron_miniPFRelIso_chg);
   fChain->SetBranchAddress("Electron_mvaFall17V1Iso",                  &Electron_mvaFall17V1Iso                , &b_Electron_mvaFall17V1Iso);
   fChain->SetBranchAddress("Electron_mvaFall17V1noIso",                &Electron_mvaFall17V1noIso              , &b_Electron_mvaFall17V1noIso);
   fChain->SetBranchAddress("Electron_mvaFall17V2Iso",                  &Electron_mvaFall17V2Iso                , &b_Electron_mvaFall17V2Iso);
   fChain->SetBranchAddress("Electron_mvaFall17V2noIso",                &Electron_mvaFall17V2noIso              , &b_Electron_mvaFall17V2noIso);
   fChain->SetBranchAddress("Electron_pfRelIso03_all",                  &Electron_pfRelIso03_all                , &b_Electron_pfRelIso03_all);
   fChain->SetBranchAddress("Electron_pfRelIso03_chg",                  &Electron_pfRelIso03_chg                , &b_Electron_pfRelIso03_chg);
   fChain->SetBranchAddress("Electron_phi",                             &Electron_phi                           , &b_Electron_phi);
   fChain->SetBranchAddress("Electron_pt",                              &Electron_pt                            , &b_Electron_pt);
   fChain->SetBranchAddress("Electron_r9",                              &Electron_r9                            , &b_Electron_r9);
   fChain->SetBranchAddress("Electron_scEtOverPt",                      &Electron_scEtOverPt                    , &b_Electron_scEtOverPt);
   fChain->SetBranchAddress("Electron_sieie",                           &Electron_sieie                         , &b_Electron_sieie);
   fChain->SetBranchAddress("Electron_sip3d",                           &Electron_sip3d                         , &b_Electron_sip3d);
   fChain->SetBranchAddress("Electron_mvaTTH",                          &Electron_mvaTTH                        , &b_Electron_mvaTTH);
   fChain->SetBranchAddress("Electron_charge",                          &Electron_charge                        , &b_Electron_charge);
   fChain->SetBranchAddress("Electron_cutBased",                        &Electron_cutBased                      , &b_Electron_cutBased);
   fChain->SetBranchAddress("Electron_cutBased_Fall17_V1",              &Electron_cutBased_Fall17_V1            , &b_Electron_cutBased_Fall17_V1);
   fChain->SetBranchAddress("Electron_jetIdx",                          &Electron_jetIdx                        , &b_Electron_jetIdx);
   fChain->SetBranchAddress("Electron_pdgId",                           &Electron_pdgId                         , &b_Electron_pdgId);
   fChain->SetBranchAddress("Electron_photonIdx",                       &Electron_photonIdx                     , &b_Electron_photonIdx);
   fChain->SetBranchAddress("Electron_tightCharge",                     &Electron_tightCharge                   , &b_Electron_tightCharge);
   fChain->SetBranchAddress("Electron_vidNestedWPBitmap",               &Electron_vidNestedWPBitmap             , &b_Electron_vidNestedWPBitmap);
   fChain->SetBranchAddress("Electron_vidNestedWPBitmapHEEP",           &Electron_vidNestedWPBitmapHEEP         , &b_Electron_vidNestedWPBitmapHEEP);
   fChain->SetBranchAddress("Electron_convVeto",                        &Electron_convVeto                      , &b_Electron_convVeto);
   fChain->SetBranchAddress("Electron_cutBased_HEEP",                   &Electron_cutBased_HEEP                 , &b_Electron_cutBased_HEEP);
   fChain->SetBranchAddress("Electron_isPFcand",                        &Electron_isPFcand                      , &b_Electron_isPFcand);
   fChain->SetBranchAddress("Electron_lostHits",                        &Electron_lostHits                      , &b_Electron_lostHits);
   fChain->SetBranchAddress("Electron_mvaFall17V1Iso_WP80",             &Electron_mvaFall17V1Iso_WP80           , &b_Electron_mvaFall17V1Iso_WP80);
   fChain->SetBranchAddress("Electron_mvaFall17V1Iso_WP90",             &Electron_mvaFall17V1Iso_WP90           , &b_Electron_mvaFall17V1Iso_WP90);
   fChain->SetBranchAddress("Electron_mvaFall17V1Iso_WPL",              &Electron_mvaFall17V1Iso_WPL            , &b_Electron_mvaFall17V1Iso_WPL);
   fChain->SetBranchAddress("Electron_mvaFall17V1noIso_WP80",           &Electron_mvaFall17V1noIso_WP80         , &b_Electron_mvaFall17V1noIso_WP80);
   fChain->SetBranchAddress("Electron_mvaFall17V1noIso_WP90",           &Electron_mvaFall17V1noIso_WP90         , &b_Electron_mvaFall17V1noIso_WP90);
   fChain->SetBranchAddress("Electron_mvaFall17V1noIso_WPL",            &Electron_mvaFall17V1noIso_WPL          , &b_Electron_mvaFall17V1noIso_WPL);
   fChain->SetBranchAddress("Electron_mvaFall17V2Iso_WP80",             &Electron_mvaFall17V2Iso_WP80           , &b_Electron_mvaFall17V2Iso_WP80);
   fChain->SetBranchAddress("Electron_mvaFall17V2Iso_WP90",             &Electron_mvaFall17V2Iso_WP90           , &b_Electron_mvaFall17V2Iso_WP90);
   fChain->SetBranchAddress("Electron_mvaFall17V2Iso_WPL",              &Electron_mvaFall17V2Iso_WPL            , &b_Electron_mvaFall17V2Iso_WPL);
   fChain->SetBranchAddress("Electron_mvaFall17V2noIso_WP80",           &Electron_mvaFall17V2noIso_WP80         , &b_Electron_mvaFall17V2noIso_WP80);
   fChain->SetBranchAddress("Electron_mvaFall17V2noIso_WP90",           &Electron_mvaFall17V2noIso_WP90         , &b_Electron_mvaFall17V2noIso_WP90);
   fChain->SetBranchAddress("Electron_mvaFall17V2noIso_WPL",            &Electron_mvaFall17V2noIso_WPL          , &b_Electron_mvaFall17V2noIso_WPL);
   fChain->SetBranchAddress("Electron_seedGain",                        &Electron_seedGain                      , &b_Electron_seedGain);
   fChain->SetBranchAddress("nFsrPhoton",                               &nFsrPhoton                             , &b_nFsrPhoton);
   fChain->SetBranchAddress("FsrPhoton_dROverEt2",                      &FsrPhoton_dROverEt2                    , &b_FsrPhoton_dROverEt2);
   fChain->SetBranchAddress("FsrPhoton_eta",                            &FsrPhoton_eta                          , &b_FsrPhoton_eta);
   fChain->SetBranchAddress("FsrPhoton_phi",                            &FsrPhoton_phi                          , &b_FsrPhoton_phi);
   fChain->SetBranchAddress("FsrPhoton_pt",                             &FsrPhoton_pt                           , &b_FsrPhoton_pt);
   fChain->SetBranchAddress("FsrPhoton_relIso03",                       &FsrPhoton_relIso03                     , &b_FsrPhoton_relIso03);
   fChain->SetBranchAddress("FsrPhoton_muonIdx",                        &FsrPhoton_muonIdx                      , &b_FsrPhoton_muonIdx);
   fChain->SetBranchAddress("nGenPart",                                 &nGenPart                               , &b_nGenPart);
   fChain->SetBranchAddress("GenPart_eta",                              &GenPart_eta                            , &b_GenPart_eta);
   fChain->SetBranchAddress("GenPart_mass",                             &GenPart_mass                           , &b_GenPart_mass);
   fChain->SetBranchAddress("GenPart_phi",                              &GenPart_phi                            , &b_GenPart_phi);
   fChain->SetBranchAddress("GenPart_pt",                               &GenPart_pt                             , &b_GenPart_pt);
   fChain->SetBranchAddress("GenPart_genPartIdxMother",                 &GenPart_genPartIdxMother               , &b_GenPart_genPartIdxMother);
   fChain->SetBranchAddress("GenPart_pdgId",                            &GenPart_pdgId                          , &b_GenPart_pdgId);
   fChain->SetBranchAddress("GenPart_status",                           &GenPart_status                         , &b_GenPart_status);
   fChain->SetBranchAddress("GenPart_statusFlags",                      &GenPart_statusFlags                    , &b_GenPart_statusFlags);
   fChain->SetBranchAddress("Generator_binvar",                         &Generator_binvar                       , &b_Generator_binvar);
   fChain->SetBranchAddress("Generator_scalePDF",                       &Generator_scalePDF                     , &b_Generator_scalePDF);
   fChain->SetBranchAddress("Generator_weight",                         &Generator_weight                       , &b_Generator_weight);
   fChain->SetBranchAddress("Generator_x1",                             &Generator_x1                           , &b_Generator_x1);
   fChain->SetBranchAddress("Generator_x2",                             &Generator_x2                           , &b_Generator_x2);
   fChain->SetBranchAddress("Generator_xpdf1",                          &Generator_xpdf1                        , &b_Generator_xpdf1);
   fChain->SetBranchAddress("Generator_xpdf2",                          &Generator_xpdf2                        , &b_Generator_xpdf2);
   fChain->SetBranchAddress("Generator_id1",                            &Generator_id1                          , &b_Generator_id1);
   fChain->SetBranchAddress("Generator_id2",                            &Generator_id2                          , &b_Generator_id2);
   fChain->SetBranchAddress("genWeight",                                &genWeight                              , &b_genWeight);
   fChain->SetBranchAddress("L1PreFiringWeight_Dn",                     &L1PreFiringWeight_Dn                   , &b_L1PreFiringWeight_Dn);
   fChain->SetBranchAddress("L1PreFiringWeight_Nom",                    &L1PreFiringWeight_Nom                  , &b_L1PreFiringWeight_Nom);
   fChain->SetBranchAddress("L1PreFiringWeight_Up",                     &L1PreFiringWeight_Up                   , &b_L1PreFiringWeight_Up);
   fChain->SetBranchAddress("nMuon",                                    &nMuon                                  , &b_nMuon);
   fChain->SetBranchAddress("Muon_dxy",                                 &Muon_dxy                               , &b_Muon_dxy);
   fChain->SetBranchAddress("Muon_dxyErr",                              &Muon_dxyErr                            , &b_Muon_dxyErr);
   fChain->SetBranchAddress("Muon_dxybs",                               &Muon_dxybs                             , &b_Muon_dxybs);
   fChain->SetBranchAddress("Muon_dz",                                  &Muon_dz                                , &b_Muon_dz);
   fChain->SetBranchAddress("Muon_dzErr",                               &Muon_dzErr                             , &b_Muon_dzErr);
   fChain->SetBranchAddress("Muon_eta",                                 &Muon_eta                               , &b_Muon_eta);
   fChain->SetBranchAddress("Muon_ip3d",                                &Muon_ip3d                              , &b_Muon_ip3d);
   fChain->SetBranchAddress("Muon_jetPtRelv2",                          &Muon_jetPtRelv2                        , &b_Muon_jetPtRelv2);
   fChain->SetBranchAddress("Muon_jetRelIso",                           &Muon_jetRelIso                         , &b_Muon_jetRelIso);
   fChain->SetBranchAddress("Muon_mass",                                &Muon_mass                              , &b_Muon_mass);
   fChain->SetBranchAddress("Muon_miniPFRelIso_all",                    &Muon_miniPFRelIso_all                  , &b_Muon_miniPFRelIso_all);
   fChain->SetBranchAddress("Muon_miniPFRelIso_chg",                    &Muon_miniPFRelIso_chg                  , &b_Muon_miniPFRelIso_chg);
   fChain->SetBranchAddress("Muon_pfRelIso03_all",                      &Muon_pfRelIso03_all                    , &b_Muon_pfRelIso03_all);
   fChain->SetBranchAddress("Muon_pfRelIso03_chg",                      &Muon_pfRelIso03_chg                    , &b_Muon_pfRelIso03_chg);
   fChain->SetBranchAddress("Muon_pfRelIso04_all",                      &Muon_pfRelIso04_all                    , &b_Muon_pfRelIso04_all);
   fChain->SetBranchAddress("Muon_phi",                                 &Muon_phi                               , &b_Muon_phi);
   fChain->SetBranchAddress("Muon_pt",                                  &Muon_pt                                , &b_Muon_pt);
   fChain->SetBranchAddress("Muon_ptErr",                               &Muon_ptErr                             , &b_Muon_ptErr);
   fChain->SetBranchAddress("Muon_segmentComp",                         &Muon_segmentComp                       , &b_Muon_segmentComp);
   fChain->SetBranchAddress("Muon_sip3d",                               &Muon_sip3d                             , &b_Muon_sip3d);
   fChain->SetBranchAddress("Muon_softMva",                             &Muon_softMva                           , &b_Muon_softMva);
   fChain->SetBranchAddress("Muon_tkRelIso",                            &Muon_tkRelIso                          , &b_Muon_tkRelIso);
   fChain->SetBranchAddress("Muon_tunepRelPt",                          &Muon_tunepRelPt                        , &b_Muon_tunepRelPt);
   fChain->SetBranchAddress("Muon_mvaLowPt",                            &Muon_mvaLowPt                          , &b_Muon_mvaLowPt);
   fChain->SetBranchAddress("Muon_mvaTTH",                              &Muon_mvaTTH                            , &b_Muon_mvaTTH);
   fChain->SetBranchAddress("Muon_charge",                              &Muon_charge                            , &b_Muon_charge);
   fChain->SetBranchAddress("Muon_jetIdx",                              &Muon_jetIdx                            , &b_Muon_jetIdx);
   fChain->SetBranchAddress("Muon_nStations",                           &Muon_nStations                         , &b_Muon_nStations);
   fChain->SetBranchAddress("Muon_nTrackerLayers",                      &Muon_nTrackerLayers                    , &b_Muon_nTrackerLayers);
   fChain->SetBranchAddress("Muon_pdgId",                               &Muon_pdgId                             , &b_Muon_pdgId);
   fChain->SetBranchAddress("Muon_tightCharge",                         &Muon_tightCharge                       , &b_Muon_tightCharge);
   fChain->SetBranchAddress("Muon_fsrPhotonIdx",                        &Muon_fsrPhotonIdx                      , &b_Muon_fsrPhotonIdx);
   fChain->SetBranchAddress("Muon_highPtId",                            &Muon_highPtId                          , &b_Muon_highPtId);
   fChain->SetBranchAddress("Muon_highPurity",                          &Muon_highPurity                        , &b_Muon_highPurity);
   fChain->SetBranchAddress("Muon_inTimeMuon",                          &Muon_inTimeMuon                        , &b_Muon_inTimeMuon);
   fChain->SetBranchAddress("Muon_isGlobal",                            &Muon_isGlobal                          , &b_Muon_isGlobal);
   fChain->SetBranchAddress("Muon_isPFcand",                            &Muon_isPFcand                          , &b_Muon_isPFcand);
   fChain->SetBranchAddress("Muon_isTracker",                           &Muon_isTracker                         , &b_Muon_isTracker);
   fChain->SetBranchAddress("Muon_looseId",                             &Muon_looseId                           , &b_Muon_looseId);
   fChain->SetBranchAddress("Muon_mediumId",                            &Muon_mediumId                          , &b_Muon_mediumId);
   fChain->SetBranchAddress("Muon_mediumPromptId",                      &Muon_mediumPromptId                    , &b_Muon_mediumPromptId);
   fChain->SetBranchAddress("Muon_miniIsoId",                           &Muon_miniIsoId                         , &b_Muon_miniIsoId);
   fChain->SetBranchAddress("Muon_multiIsoId",                          &Muon_multiIsoId                        , &b_Muon_multiIsoId);
   fChain->SetBranchAddress("Muon_mvaId",                               &Muon_mvaId                             , &b_Muon_mvaId);
   fChain->SetBranchAddress("Muon_mvaLowPtId",                          &Muon_mvaLowPtId                        , &b_Muon_mvaLowPtId);
   fChain->SetBranchAddress("Muon_pfIsoId",                             &Muon_pfIsoId                           , &b_Muon_pfIsoId);
   fChain->SetBranchAddress("Muon_puppiIsoId",                          &Muon_puppiIsoId                        , &b_Muon_puppiIsoId);
   fChain->SetBranchAddress("Muon_softId",                              &Muon_softId                            , &b_Muon_softId);
   fChain->SetBranchAddress("Muon_softMvaId",                           &Muon_softMvaId                         , &b_Muon_softMvaId);
   fChain->SetBranchAddress("Muon_tightId",                             &Muon_tightId                           , &b_Muon_tightId);
   fChain->SetBranchAddress("Muon_tkIsoId",                             &Muon_tkIsoId                           , &b_Muon_tkIsoId);
   fChain->SetBranchAddress("Muon_triggerIdLoose",                      &Muon_triggerIdLoose                    , &b_Muon_triggerIdLoose);
   fChain->SetBranchAddress("Pileup_nTrueInt",                          &Pileup_nTrueInt                        , &b_Pileup_nTrueInt);
   fChain->SetBranchAddress("nGenDressedLepton",                        &nGenDressedLepton                      , &b_nGenDressedLepton);
   fChain->SetBranchAddress("GenDressedLepton_eta",                     &GenDressedLepton_eta                   , &b_GenDressedLepton_eta);
   fChain->SetBranchAddress("GenDressedLepton_mass",                    &GenDressedLepton_mass                  , &b_GenDressedLepton_mass);
   fChain->SetBranchAddress("GenDressedLepton_phi",                     &GenDressedLepton_phi                   , &b_GenDressedLepton_phi);
   fChain->SetBranchAddress("GenDressedLepton_pt",                      &GenDressedLepton_pt                    , &b_GenDressedLepton_pt);
   fChain->SetBranchAddress("GenDressedLepton_pdgId",                   &GenDressedLepton_pdgId                 , &b_GenDressedLepton_pdgId);
   fChain->SetBranchAddress("GenDressedLepton_hasTauAnc",               &GenDressedLepton_hasTauAnc             , &b_GenDressedLepton_hasTauAnc);
   fChain->SetBranchAddress("nTrigObj",                                 &nTrigObj                               , &b_nTrigObj);
   fChain->SetBranchAddress("TrigObj_pt",                               &TrigObj_pt                             , &b_TrigObj_pt);
   fChain->SetBranchAddress("TrigObj_eta",                              &TrigObj_eta                            , &b_TrigObj_eta);
   fChain->SetBranchAddress("TrigObj_phi",                              &TrigObj_phi                            , &b_TrigObj_phi);
   fChain->SetBranchAddress("TrigObj_l1pt",                             &TrigObj_l1pt                           , &b_TrigObj_l1pt);
   fChain->SetBranchAddress("TrigObj_l1pt_2",                           &TrigObj_l1pt_2                         , &b_TrigObj_l1pt_2);
   fChain->SetBranchAddress("TrigObj_l2pt",                             &TrigObj_l2pt                           , &b_TrigObj_l2pt);
   fChain->SetBranchAddress("TrigObj_id",                               &TrigObj_id                             , &b_TrigObj_id);
   fChain->SetBranchAddress("TrigObj_l1iso",                            &TrigObj_l1iso                          , &b_TrigObj_l1iso);
   fChain->SetBranchAddress("TrigObj_l1charge",                         &TrigObj_l1charge                       , &b_TrigObj_l1charge);
   fChain->SetBranchAddress("TrigObj_filterBits",                       &TrigObj_filterBits                     , &b_TrigObj_filterBits);
   fChain->SetBranchAddress("nOtherPV",                                 &nOtherPV                               , &b_nOtherPV);
   fChain->SetBranchAddress("OtherPV_z",                                &OtherPV_z                              , &b_OtherPV_z);
   fChain->SetBranchAddress("PV_ndof",                                  &PV_ndof                                , &b_PV_ndof);
   fChain->SetBranchAddress("PV_x",                                     &PV_x                                   , &b_PV_x);
   fChain->SetBranchAddress("PV_y",                                     &PV_y                                   , &b_PV_y);
   fChain->SetBranchAddress("PV_z",                                     &PV_z                                   , &b_PV_z);
   fChain->SetBranchAddress("PV_chi2",                                  &PV_chi2                                , &b_PV_chi2);
   fChain->SetBranchAddress("PV_score",                                 &PV_score                               , &b_PV_score);
   fChain->SetBranchAddress("PV_npvs",                                  &PV_npvs                                , &b_PV_npvs);
   fChain->SetBranchAddress("PV_npvsGood",                              &PV_npvsGood                            , &b_PV_npvsGood);
   fChain->SetBranchAddress("Electron_genPartIdx",                      &Electron_genPartIdx                    , &b_Electron_genPartIdx);
   fChain->SetBranchAddress("Electron_genPartFlav",                     &Electron_genPartFlav                   , &b_Electron_genPartFlav);
   fChain->SetBranchAddress("Muon_genPartIdx",                          &Muon_genPartIdx                        , &b_Muon_genPartIdx);
   fChain->SetBranchAddress("Muon_genPartFlav",                         &Muon_genPartFlav                       , &b_Muon_genPartFlav);
   fChain->SetBranchAddress("Electron_cleanmask",                       &Electron_cleanmask                     , &b_Electron_cleanmask);
   fChain->SetBranchAddress("Muon_cleanmask",                           &Muon_cleanmask                         , &b_Muon_cleanmask);
   fChain->SetBranchAddress("Flag_HBHENoiseFilter",                     &Flag_HBHENoiseFilter                   , &b_Flag_HBHENoiseFilter);
   fChain->SetBranchAddress("Flag_HBHENoiseIsoFilter",                  &Flag_HBHENoiseIsoFilter                , &b_Flag_HBHENoiseIsoFilter);
   fChain->SetBranchAddress("Flag_CSCTightHaloFilter",                  &Flag_CSCTightHaloFilter                , &b_Flag_CSCTightHaloFilter);
   fChain->SetBranchAddress("Flag_CSCTightHaloTrkMuUnvetoFilter",       &Flag_CSCTightHaloTrkMuUnvetoFilter     , &b_Flag_CSCTightHaloTrkMuUnvetoFilter);
   fChain->SetBranchAddress("Flag_CSCTightHalo2015Filter",              &Flag_CSCTightHalo2015Filter            , &b_Flag_CSCTightHalo2015Filter);
   fChain->SetBranchAddress("Flag_globalTightHalo2016Filter",           &Flag_globalTightHalo2016Filter         , &b_Flag_globalTightHalo2016Filter);
   fChain->SetBranchAddress("Flag_globalSuperTightHalo2016Filter",      &Flag_globalSuperTightHalo2016Filter    , &b_Flag_globalSuperTightHalo2016Filter);
   fChain->SetBranchAddress("Flag_HcalStripHaloFilter",                 &Flag_HcalStripHaloFilter               , &b_Flag_HcalStripHaloFilter);
   fChain->SetBranchAddress("Flag_hcalLaserEventFilter",                &Flag_hcalLaserEventFilter              , &b_Flag_hcalLaserEventFilter);
   fChain->SetBranchAddress("Flag_EcalDeadCellTriggerPrimitiveFilter",  &Flag_EcalDeadCellTriggerPrimitiveFilter, &b_Flag_EcalDeadCellTriggerPrimitiveFilter);
   fChain->SetBranchAddress("Flag_EcalDeadCellBoundaryEnergyFilter",    &Flag_EcalDeadCellBoundaryEnergyFilter  , &b_Flag_EcalDeadCellBoundaryEnergyFilter);
   fChain->SetBranchAddress("Flag_ecalBadCalibFilter",                  &Flag_ecalBadCalibFilter                , &b_Flag_ecalBadCalibFilter);
   fChain->SetBranchAddress("Flag_goodVertices",                        &Flag_goodVertices                      , &b_Flag_goodVertices);
   fChain->SetBranchAddress("Flag_eeBadScFilter",                       &Flag_eeBadScFilter                     , &b_Flag_eeBadScFilter);
   fChain->SetBranchAddress("Flag_ecalLaserCorrFilter",                 &Flag_ecalLaserCorrFilter               , &b_Flag_ecalLaserCorrFilter);
   fChain->SetBranchAddress("Flag_trkPOGFilters",                       &Flag_trkPOGFilters                     , &b_Flag_trkPOGFilters);
   fChain->SetBranchAddress("Flag_chargedHadronTrackResolutionFilter",  &Flag_chargedHadronTrackResolutionFilter, &b_Flag_chargedHadronTrackResolutionFilter);
   fChain->SetBranchAddress("Flag_muonBadTrackFilter",                  &Flag_muonBadTrackFilter                , &b_Flag_muonBadTrackFilter);
   fChain->SetBranchAddress("Flag_BadChargedCandidateFilter",           &Flag_BadChargedCandidateFilter         , &b_Flag_BadChargedCandidateFilter);
   fChain->SetBranchAddress("Flag_BadPFMuonFilter",                     &Flag_BadPFMuonFilter                   , &b_Flag_BadPFMuonFilter);
   fChain->SetBranchAddress("Flag_BadChargedCandidateSummer16Filter",   &Flag_BadChargedCandidateSummer16Filter , &b_Flag_BadChargedCandidateSummer16Filter);
   fChain->SetBranchAddress("Flag_BadPFMuonSummer16Filter",             &Flag_BadPFMuonSummer16Filter           , &b_Flag_BadPFMuonSummer16Filter);
   fChain->SetBranchAddress("Flag_trkPOG_manystripclus53X",             &Flag_trkPOG_manystripclus53X           , &b_Flag_trkPOG_manystripclus53X);
   fChain->SetBranchAddress("Flag_trkPOG_toomanystripclus53X",          &Flag_trkPOG_toomanystripclus53X        , &b_Flag_trkPOG_toomanystripclus53X);
   fChain->SetBranchAddress("Flag_trkPOG_logErrorTooManyClusters",      &Flag_trkPOG_logErrorTooManyClusters    , &b_Flag_trkPOG_logErrorTooManyClusters);
   fChain->SetBranchAddress("Flag_METFilters",                          &Flag_METFilters                        , &b_Flag_METFilters);
   fChain->SetBranchAddress("HLT_IsoMu24",                              &HLT_IsoMu24                            , &b_HLT_IsoMu24);
   fChain->SetBranchAddress("HLT_IsoTkMu24",                            &HLT_IsoTkMu24                          , &b_HLT_IsoTkMu24);
   Notify();
}

Bool_t TnPNtuplesBase::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void TnPNtuplesBase::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t TnPNtuplesBase::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef TnPNtuplesBase_cxx
