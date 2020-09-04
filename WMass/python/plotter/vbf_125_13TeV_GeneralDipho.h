//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Sat Aug 29 17:54:48 2020 by ROOT version 6.12/07
// from TTree vbf_125_13TeV_GeneralDipho/vbf_125_13TeV_GeneralDipho
// found on file: /eos/home-e/escott/HggLegacy/TrainingNtuples/Pass5/2016/VBF.root
//////////////////////////////////////////////////////////

#ifndef vbf_125_13TeV_GeneralDipho_h
#define vbf_125_13TeV_GeneralDipho_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

class vbf_125_13TeV_GeneralDipho {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Int_t           candidate_id;
   Float_t         weight;
   Float_t         dipho_sumpt;
   Float_t         dipho_cosphi;
   Float_t         dipho_mass;
   Float_t         dipho_leadPt;
   Float_t         dipho_leadEt;
   Float_t         dipho_leadEta;
   Float_t         dipho_leadPhi;
   Float_t         dipho_lead_sieie;
   Float_t         dipho_lead_hoe;
   Float_t         dipho_lead_sigmaEoE;
   Float_t         dipho_lead_ptoM;
   Float_t         dipho_leadR9;
   Float_t         dipho_subleadPt;
   Float_t         dipho_subleadEt;
   Float_t         dipho_subleadEta;
   Float_t         dipho_subleadPhi;
   Float_t         dipho_sublead_sieie;
   Float_t         dipho_sublead_hoe;
   Float_t         dipho_sublead_sigmaEoE;
   Float_t         dipho_sublead_ptoM;
   Float_t         dipho_subleadR9;
   Float_t         dipho_leadIDMVA;
   Float_t         dipho_subleadIDMVA;
   Float_t         dipho_lead_elveto;
   Float_t         dipho_sublead_elveto;
   Float_t         result;
   Float_t         dipho_PToM;
   Float_t         sigmarv;
   Float_t         sigmarvDecorr;
   Float_t         sigmawv;
   Float_t         CosPhi;
   Float_t         vtxprob;
   Float_t         pt;
   Float_t         leadSCeta;
   Float_t         subleadSCeta;
   Float_t         dijet_abs_dEta;
   Float_t         dijet_leadEta;
   Float_t         dijet_subleadEta;
   Float_t         dijet_leady;
   Float_t         dijet_subleady;
   Float_t         dijet_LeadJPt;
   Float_t         dijet_SubJPt;
   Float_t         dijet_Zep;
   Float_t         dijet_Mjj;
   Float_t         dipho_PToM;
   Float_t         leadPho_PToM;
   Float_t         sublPho_PToM;
   Float_t         dijet_dipho_dphi_trunc;
   Float_t         dijet_dipho_pt;
   Float_t         dijet_dphi;
   Float_t         dijet_dipho_dphi;
   Float_t         dijet_dPhi_trunc;
   Float_t         cos_dijet_dipho_dphi;
   Float_t         dijet_minDRJetPho;
   Float_t         has3Jet;
   Float_t         dijet_mva;
   Float_t         dipho_dijet_MVA;
   Float_t         dipho_mva;
   Float_t         dijet_dipho_dphi_trunc;
   Float_t         jet1_pt;
   Float_t         jet2_pt;
   Float_t         jet3_pt;
   Float_t         jet1_eta;
   Float_t         jet2_eta;
   Float_t         jet3_eta;
   Float_t         Mjj;
   Float_t         jet1_rawPt;
   Float_t         jet2_rawPt;
   Float_t         jet1_HFHadronEnergyFraction;
   Float_t         jet1_HFEMEnergyFraction;
   Float_t         jet1_HFHadronEnergy;
   Float_t         jet1_HFEMEnergy;
   Float_t         jet1_HFHadronMultiplicity;
   Float_t         jet1_HFEMMultiplicity;
   Float_t         jet2_HFHadronEnergyFraction;
   Float_t         jet2_HFEMEnergyFraction;
   Float_t         jet2_HFHadronEnergy;
   Float_t         jet2_HFEMEnergy;
   Float_t         jet2_HFHadronMultiplicity;
   Float_t         jet2_HFEMMultiplicity;
   Float_t         dijet_nj;
   Float_t         dipho_dijet_ptHjj;
   Float_t         n_jet_30;
   Float_t         dijet_jet1_RMS;
   Float_t         dijet_jet2_RMS;
   Float_t         dijet_jet1_QGL;
   Float_t         dijet_jet2_QGL;
   Float_t         dijet_jet1_pujid_mva;
   Float_t         dijet_jet2_pujid_mva;
   Float_t         dipho_pt;
   Float_t         dijet_pt;
   Float_t         gghMVA_n_rec_jets;
   Float_t         gghMVA_Mjj;
   Float_t         gghMVA_leadEta;
   Float_t         gghMVA_subleadEta;
   Float_t         gghMVA_subsubleadEta;
   Float_t         gghMVA_leadJPt;
   Float_t         gghMVA_SubleadJPt;
   Float_t         gghMVA_SubsubleadJPt;
   Float_t         gghMVA_leadPUMVA;
   Float_t         gghMVA_subleadPUMVA;
   Float_t         gghMVA_subsubleadPUMVA;
   Float_t         gghMVA_leadDeltaPhi;
   Float_t         gghMVA_subleadDeltaPhi;
   Float_t         gghMVA_subsubleadDeltaPhi;
   Float_t         gghMVA_leadDeltaEta;
   Float_t         gghMVA_subleadDeltaEta;
   Float_t         gghMVA_subsubleadDeltaEta;
   Float_t         ggHMVAResult_prob_0J_PTH_0_10;
   Float_t         ggHMVAResult_prob_0J_PTH_GT10;
   Float_t         ggHMVAResult_prob_1J_PTH_0_60;
   Float_t         ggHMVAResult_prob_1J_PTH_60_120;
   Float_t         ggHMVAResult_prob_1J_PTH_120_200;
   Float_t         ggHMVAResult_prob_GE2J_MJJ_0_350_PTH_0_60;
   Float_t         ggHMVAResult_prob_GE3J_MJJ_0_350_PTH_60_120;
   Float_t         ggHMVAResult_prob_GE2J_MJJ_0_350_PTH_120_200;
   Float_t         ggHMVAResult_prob_PTH_GT200;
   Float_t         vbfMvaResult_prob_bkg;
   Float_t         vbfMvaResult_prob_ggH;
   Float_t         vbfMvaResult_prob_VBF;
   Float_t         dijet_minDRJetPho;
   Float_t         dijet_centrality_gg;
   Float_t         dijet_centrality_j3;
   Float_t         dijet_centrality_g;
   Float_t         cosThetaStar;
   Float_t         VH_had_mvascore;
   Float_t         dijet_jet1_match;
   Float_t         dijet_jet2_match;
   Float_t         prompt_pho_1;
   Float_t         prompt_pho_2;
   Float_t         HTXSstage0bin;
   Float_t         HTXSstage1bin;
   Float_t         HTXSstage0bin;
   Float_t         HTXSstage1bin;
   Float_t         HTXSstage1p1bin;
   Float_t         HTXSstage1p1binFine;
   Float_t         HTXSstage1p2bin;
   Float_t         HTXSstage1p2binFine;
   Float_t         HTXSnjets;
   Float_t         HTXSpTH;
   Float_t         HTXSpTV;
   Float_t         rho;
   Int_t           nvtx;
   ULong64_t       event;
   UInt_t          lumi;
   Int_t           processIndex;
   UInt_t          run;
   Int_t           nvtx;
   Float_t         npu;
   Float_t         puweight;

   // List of branches
   TBranch        *b_candidate_id;   //!
   TBranch        *b_weight;   //!
   TBranch        *b_dipho_sumpt;   //!
   TBranch        *b_dipho_cosphi;   //!
   TBranch        *b_dipho_mass;   //!
   TBranch        *b_dipho_leadPt;   //!
   TBranch        *b_dipho_leadEt;   //!
   TBranch        *b_dipho_leadEta;   //!
   TBranch        *b_dipho_leadPhi;   //!
   TBranch        *b_dipho_lead_sieie;   //!
   TBranch        *b_dipho_lead_hoe;   //!
   TBranch        *b_dipho_lead_sigmaEoE;   //!
   TBranch        *b_dipho_lead_ptoM;   //!
   TBranch        *b_dipho_leadR9;   //!
   TBranch        *b_dipho_subleadPt;   //!
   TBranch        *b_dipho_subleadEt;   //!
   TBranch        *b_dipho_subleadEta;   //!
   TBranch        *b_dipho_subleadPhi;   //!
   TBranch        *b_dipho_sublead_sieie;   //!
   TBranch        *b_dipho_sublead_hoe;   //!
   TBranch        *b_dipho_sublead_sigmaEoE;   //!
   TBranch        *b_dipho_sublead_ptoM;   //!
   TBranch        *b_dipho_subleadR9;   //!
   TBranch        *b_dipho_leadIDMVA;   //!
   TBranch        *b_dipho_subleadIDMVA;   //!
   TBranch        *b_dipho_lead_elveto;   //!
   TBranch        *b_dipho_sublead_elveto;   //!
   TBranch        *b_result;   //!
   TBranch        *b_dipho_PToM;   //!
   TBranch        *b_sigmarv;   //!
   TBranch        *b_sigmarvDecorr;   //!
   TBranch        *b_sigmawv;   //!
   TBranch        *b_CosPhi;   //!
   TBranch        *b_vtxprob;   //!
   TBranch        *b_pt;   //!
   TBranch        *b_leadSCeta;   //!
   TBranch        *b_subleadSCeta;   //!
   TBranch        *b_dijet_abs_dEta;   //!
   TBranch        *b_dijet_leadEta;   //!
   TBranch        *b_dijet_subleadEta;   //!
   TBranch        *b_dijet_leady;   //!
   TBranch        *b_dijet_subleady;   //!
   TBranch        *b_dijet_LeadJPt;   //!
   TBranch        *b_dijet_SubJPt;   //!
   TBranch        *b_dijet_Zep;   //!
   TBranch        *b_dijet_Mjj;   //!
   TBranch        *b_dipho_PToM;   //!
   TBranch        *b_leadPho_PToM;   //!
   TBranch        *b_sublPho_PToM;   //!
   TBranch        *b_dijet_dipho_dphi_trunc;   //!
   TBranch        *b_dijet_dipho_pt;   //!
   TBranch        *b_dijet_dphi;   //!
   TBranch        *b_dijet_dipho_dphi;   //!
   TBranch        *b_dijet_dPhi_trunc;   //!
   TBranch        *b_cos_dijet_dipho_dphi;   //!
   TBranch        *b_dijet_minDRJetPho;   //!
   TBranch        *b_has3Jet;   //!
   TBranch        *b_dijet_mva;   //!
   TBranch        *b_dipho_dijet_MVA;   //!
   TBranch        *b_dipho_mva;   //!
   TBranch        *b_dijet_dipho_dphi_trunc;   //!
   TBranch        *b_jet1_pt;   //!
   TBranch        *b_jet2_pt;   //!
   TBranch        *b_jet3_pt;   //!
   TBranch        *b_jet1_eta;   //!
   TBranch        *b_jet2_eta;   //!
   TBranch        *b_jet3_eta;   //!
   TBranch        *b_Mjj;   //!
   TBranch        *b_jet1_rawPt;   //!
   TBranch        *b_jet2_rawPt;   //!
   TBranch        *b_jet1_HFHadronEnergyFraction;   //!
   TBranch        *b_jet1_HFEMEnergyFraction;   //!
   TBranch        *b_jet1_HFHadronEnergy;   //!
   TBranch        *b_jet1_HFEMEnergy;   //!
   TBranch        *b_jet1_HFHadronMultiplicity;   //!
   TBranch        *b_jet1_HFEMMultiplicity;   //!
   TBranch        *b_jet2_HFHadronEnergyFraction;   //!
   TBranch        *b_jet2_HFEMEnergyFraction;   //!
   TBranch        *b_jet2_HFHadronEnergy;   //!
   TBranch        *b_jet2_HFEMEnergy;   //!
   TBranch        *b_jet2_HFHadronMultiplicity;   //!
   TBranch        *b_jet2_HFEMMultiplicity;   //!
   TBranch        *b_dijet_nj;   //!
   TBranch        *b_dipho_dijet_ptHjj;   //!
   TBranch        *b_n_jet_30;   //!
   TBranch        *b_dijet_jet1_RMS;   //!
   TBranch        *b_dijet_jet2_RMS;   //!
   TBranch        *b_dijet_jet1_QGL;   //!
   TBranch        *b_dijet_jet2_QGL;   //!
   TBranch        *b_dijet_jet1_pujid_mva;   //!
   TBranch        *b_dijet_jet2_pujid_mva;   //!
   TBranch        *b_dipho_pt;   //!
   TBranch        *b_dijet_pt;   //!
   TBranch        *b_gghMVA_n_rec_jets;   //!
   TBranch        *b_gghMVA_Mjj;   //!
   TBranch        *b_gghMVA_leadEta;   //!
   TBranch        *b_gghMVA_subleadEta;   //!
   TBranch        *b_gghMVA_subsubleadEta;   //!
   TBranch        *b_gghMVA_leadJPt;   //!
   TBranch        *b_gghMVA_SubleadJPt;   //!
   TBranch        *b_gghMVA_SubsubleadJPt;   //!
   TBranch        *b_gghMVA_leadPUMVA;   //!
   TBranch        *b_gghMVA_subleadPUMVA;   //!
   TBranch        *b_gghMVA_subsubleadPUMVA;   //!
   TBranch        *b_gghMVA_leadDeltaPhi;   //!
   TBranch        *b_gghMVA_subleadDeltaPhi;   //!
   TBranch        *b_gghMVA_subsubleadDeltaPhi;   //!
   TBranch        *b_gghMVA_leadDeltaEta;   //!
   TBranch        *b_gghMVA_subleadDeltaEta;   //!
   TBranch        *b_gghMVA_subsubleadDeltaEta;   //!
   TBranch        *b_ggHMVAResult_prob_0J_PTH_0_10;   //!
   TBranch        *b_ggHMVAResult_prob_0J_PTH_GT10;   //!
   TBranch        *b_ggHMVAResult_prob_1J_PTH_0_60;   //!
   TBranch        *b_ggHMVAResult_prob_1J_PTH_60_120;   //!
   TBranch        *b_ggHMVAResult_prob_1J_PTH_120_200;   //!
   TBranch        *b_ggHMVAResult_prob_GE2J_MJJ_0_350_PTH_0_60;   //!
   TBranch        *b_ggHMVAResult_prob_GE3J_MJJ_0_350_PTH_60_120;   //!
   TBranch        *b_ggHMVAResult_prob_GE2J_MJJ_0_350_PTH_120_200;   //!
   TBranch        *b_ggHMVAResult_prob_PTH_GT200;   //!
   TBranch        *b_vbfMvaResult_prob_bkg;   //!
   TBranch        *b_vbfMvaResult_prob_ggH;   //!
   TBranch        *b_vbfMvaResult_prob_VBF;   //!
   TBranch        *b_dijet_minDRJetPho;   //!
   TBranch        *b_dijet_centrality_gg;   //!
   TBranch        *b_dijet_centrality_j3;   //!
   TBranch        *b_dijet_centrality_g;   //!
   TBranch        *b_cosThetaStar;   //!
   TBranch        *b_VH_had_mvascore;   //!
   TBranch        *b_dijet_jet1_match;   //!
   TBranch        *b_dijet_jet2_match;   //!
   TBranch        *b_prompt_pho_1;   //!
   TBranch        *b_prompt_pho_2;   //!
   TBranch        *b_HTXSstage0bin;   //!
   TBranch        *b_HTXSstage1bin;   //!
   TBranch        *b_HTXSstage0bin;   //!
   TBranch        *b_HTXSstage1bin;   //!
   TBranch        *b_HTXSstage1p1bin;   //!
   TBranch        *b_HTXSstage1p1binFine;   //!
   TBranch        *b_HTXSstage1p2bin;   //!
   TBranch        *b_HTXSstage1p2binFine;   //!
   TBranch        *b_HTXSnjets;   //!
   TBranch        *b_HTXSpTH;   //!
   TBranch        *b_HTXSpTV;   //!
   TBranch        *b_rho;   //!
   TBranch        *b_nvtx;   //!
   TBranch        *b_event;   //!
   TBranch        *b_lumi;   //!
   TBranch        *b_processIndex;   //!
   TBranch        *b_run;   //!
   TBranch        *b_nvtx;   //!
   TBranch        *b_npu;   //!
   TBranch        *b_puweight;   //!

   vbf_125_13TeV_GeneralDipho(TTree *tree=0);
   virtual ~vbf_125_13TeV_GeneralDipho();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef vbf_125_13TeV_GeneralDipho_cxx
vbf_125_13TeV_GeneralDipho::vbf_125_13TeV_GeneralDipho(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("/eos/home-e/escott/HggLegacy/TrainingNtuples/Pass5/2016/VBF.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("/eos/home-e/escott/HggLegacy/TrainingNtuples/Pass5/2016/VBF.root");
      }
      TDirectory * dir = (TDirectory*)f->Get("/eos/home-e/escott/HggLegacy/TrainingNtuples/Pass5/2016/VBF.root:/vbfTagDumper/trees");
      dir->GetObject("vbf_125_13TeV_GeneralDipho",tree);

   }
   Init(tree);
}

vbf_125_13TeV_GeneralDipho::~vbf_125_13TeV_GeneralDipho()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t vbf_125_13TeV_GeneralDipho::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t vbf_125_13TeV_GeneralDipho::LoadTree(Long64_t entry)
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

void vbf_125_13TeV_GeneralDipho::Init(TTree *tree)
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

   fChain->SetBranchAddress("candidate_id", &candidate_id, &b_candidate_id);
   fChain->SetBranchAddress("weight", &weight, &b_weight);
   fChain->SetBranchAddress("dipho_sumpt", &dipho_sumpt, &b_dipho_sumpt);
   fChain->SetBranchAddress("dipho_cosphi", &dipho_cosphi, &b_dipho_cosphi);
   fChain->SetBranchAddress("dipho_mass", &dipho_mass, &b_dipho_mass);
   fChain->SetBranchAddress("dipho_leadPt", &dipho_leadPt, &b_dipho_leadPt);
   fChain->SetBranchAddress("dipho_leadEt", &dipho_leadEt, &b_dipho_leadEt);
   fChain->SetBranchAddress("dipho_leadEta", &dipho_leadEta, &b_dipho_leadEta);
   fChain->SetBranchAddress("dipho_leadPhi", &dipho_leadPhi, &b_dipho_leadPhi);
   fChain->SetBranchAddress("dipho_lead_sieie", &dipho_lead_sieie, &b_dipho_lead_sieie);
   fChain->SetBranchAddress("dipho_lead_hoe", &dipho_lead_hoe, &b_dipho_lead_hoe);
   fChain->SetBranchAddress("dipho_lead_sigmaEoE", &dipho_lead_sigmaEoE, &b_dipho_lead_sigmaEoE);
   fChain->SetBranchAddress("dipho_lead_ptoM", &dipho_lead_ptoM, &b_dipho_lead_ptoM);
   fChain->SetBranchAddress("dipho_leadR9", &dipho_leadR9, &b_dipho_leadR9);
   fChain->SetBranchAddress("dipho_subleadPt", &dipho_subleadPt, &b_dipho_subleadPt);
   fChain->SetBranchAddress("dipho_subleadEt", &dipho_subleadEt, &b_dipho_subleadEt);
   fChain->SetBranchAddress("dipho_subleadEta", &dipho_subleadEta, &b_dipho_subleadEta);
   fChain->SetBranchAddress("dipho_subleadPhi", &dipho_subleadPhi, &b_dipho_subleadPhi);
   fChain->SetBranchAddress("dipho_sublead_sieie", &dipho_sublead_sieie, &b_dipho_sublead_sieie);
   fChain->SetBranchAddress("dipho_sublead_hoe", &dipho_sublead_hoe, &b_dipho_sublead_hoe);
   fChain->SetBranchAddress("dipho_sublead_sigmaEoE", &dipho_sublead_sigmaEoE, &b_dipho_sublead_sigmaEoE);
   fChain->SetBranchAddress("dipho_sublead_ptoM", &dipho_sublead_ptoM, &b_dipho_sublead_ptoM);
   fChain->SetBranchAddress("dipho_subleadR9", &dipho_subleadR9, &b_dipho_subleadR9);
   fChain->SetBranchAddress("dipho_leadIDMVA", &dipho_leadIDMVA, &b_dipho_leadIDMVA);
   fChain->SetBranchAddress("dipho_subleadIDMVA", &dipho_subleadIDMVA, &b_dipho_subleadIDMVA);
   fChain->SetBranchAddress("dipho_lead_elveto", &dipho_lead_elveto, &b_dipho_lead_elveto);
   fChain->SetBranchAddress("dipho_sublead_elveto", &dipho_sublead_elveto, &b_dipho_sublead_elveto);
   fChain->SetBranchAddress("result", &result, &b_result);
   fChain->SetBranchAddress("dipho_PToM", &dipho_PToM, &b_dipho_PToM);
   fChain->SetBranchAddress("sigmarv", &sigmarv, &b_sigmarv);
   fChain->SetBranchAddress("sigmarvDecorr", &sigmarvDecorr, &b_sigmarvDecorr);
   fChain->SetBranchAddress("sigmawv", &sigmawv, &b_sigmawv);
   fChain->SetBranchAddress("CosPhi", &CosPhi, &b_CosPhi);
   fChain->SetBranchAddress("vtxprob", &vtxprob, &b_vtxprob);
   fChain->SetBranchAddress("pt", &pt, &b_pt);
   fChain->SetBranchAddress("leadSCeta", &leadSCeta, &b_leadSCeta);
   fChain->SetBranchAddress("subleadSCeta", &subleadSCeta, &b_subleadSCeta);
   fChain->SetBranchAddress("dijet_abs_dEta", &dijet_abs_dEta, &b_dijet_abs_dEta);
   fChain->SetBranchAddress("dijet_leadEta", &dijet_leadEta, &b_dijet_leadEta);
   fChain->SetBranchAddress("dijet_subleadEta", &dijet_subleadEta, &b_dijet_subleadEta);
   fChain->SetBranchAddress("dijet_leady", &dijet_leady, &b_dijet_leady);
   fChain->SetBranchAddress("dijet_subleady", &dijet_subleady, &b_dijet_subleady);
   fChain->SetBranchAddress("dijet_LeadJPt", &dijet_LeadJPt, &b_dijet_LeadJPt);
   fChain->SetBranchAddress("dijet_SubJPt", &dijet_SubJPt, &b_dijet_SubJPt);
   fChain->SetBranchAddress("dijet_Zep", &dijet_Zep, &b_dijet_Zep);
   fChain->SetBranchAddress("dijet_Mjj", &dijet_Mjj, &b_dijet_Mjj);
//    fChain->SetBranchAddress("dipho_PToM", &dipho_PToM, &b_dipho_PToM);
   fChain->SetBranchAddress("leadPho_PToM", &leadPho_PToM, &b_leadPho_PToM);
   fChain->SetBranchAddress("sublPho_PToM", &sublPho_PToM, &b_sublPho_PToM);
   fChain->SetBranchAddress("dijet_dipho_dphi_trunc", &dijet_dipho_dphi_trunc, &b_dijet_dipho_dphi_trunc);
   fChain->SetBranchAddress("dijet_dipho_pt", &dijet_dipho_pt, &b_dijet_dipho_pt);
   fChain->SetBranchAddress("dijet_dphi", &dijet_dphi, &b_dijet_dphi);
   fChain->SetBranchAddress("dijet_dipho_dphi", &dijet_dipho_dphi, &b_dijet_dipho_dphi);
   fChain->SetBranchAddress("dijet_dPhi_trunc", &dijet_dPhi_trunc, &b_dijet_dPhi_trunc);
   fChain->SetBranchAddress("cos_dijet_dipho_dphi", &cos_dijet_dipho_dphi, &b_cos_dijet_dipho_dphi);
   fChain->SetBranchAddress("dijet_minDRJetPho", &dijet_minDRJetPho, &b_dijet_minDRJetPho);
   fChain->SetBranchAddress("has3Jet", &has3Jet, &b_has3Jet);
   fChain->SetBranchAddress("dijet_mva", &dijet_mva, &b_dijet_mva);
   fChain->SetBranchAddress("dipho_dijet_MVA", &dipho_dijet_MVA, &b_dipho_dijet_MVA);
   fChain->SetBranchAddress("dipho_mva", &dipho_mva, &b_dipho_mva);
//    fChain->SetBranchAddress("dijet_dipho_dphi_trunc", &dijet_dipho_dphi_trunc, &b_dijet_dipho_dphi_trunc);
   fChain->SetBranchAddress("jet1_pt", &jet1_pt, &b_jet1_pt);
   fChain->SetBranchAddress("jet2_pt", &jet2_pt, &b_jet2_pt);
   fChain->SetBranchAddress("jet3_pt", &jet3_pt, &b_jet3_pt);
   fChain->SetBranchAddress("jet1_eta", &jet1_eta, &b_jet1_eta);
   fChain->SetBranchAddress("jet2_eta", &jet2_eta, &b_jet2_eta);
   fChain->SetBranchAddress("jet3_eta", &jet3_eta, &b_jet3_eta);
   fChain->SetBranchAddress("Mjj", &Mjj, &b_Mjj);
   fChain->SetBranchAddress("jet1_rawPt", &jet1_rawPt, &b_jet1_rawPt);
   fChain->SetBranchAddress("jet2_rawPt", &jet2_rawPt, &b_jet2_rawPt);
   fChain->SetBranchAddress("jet1_HFHadronEnergyFraction", &jet1_HFHadronEnergyFraction, &b_jet1_HFHadronEnergyFraction);
   fChain->SetBranchAddress("jet1_HFEMEnergyFraction", &jet1_HFEMEnergyFraction, &b_jet1_HFEMEnergyFraction);
   fChain->SetBranchAddress("jet1_HFHadronEnergy", &jet1_HFHadronEnergy, &b_jet1_HFHadronEnergy);
   fChain->SetBranchAddress("jet1_HFEMEnergy", &jet1_HFEMEnergy, &b_jet1_HFEMEnergy);
   fChain->SetBranchAddress("jet1_HFHadronMultiplicity", &jet1_HFHadronMultiplicity, &b_jet1_HFHadronMultiplicity);
   fChain->SetBranchAddress("jet1_HFEMMultiplicity", &jet1_HFEMMultiplicity, &b_jet1_HFEMMultiplicity);
   fChain->SetBranchAddress("jet2_HFHadronEnergyFraction", &jet2_HFHadronEnergyFraction, &b_jet2_HFHadronEnergyFraction);
   fChain->SetBranchAddress("jet2_HFEMEnergyFraction", &jet2_HFEMEnergyFraction, &b_jet2_HFEMEnergyFraction);
   fChain->SetBranchAddress("jet2_HFHadronEnergy", &jet2_HFHadronEnergy, &b_jet2_HFHadronEnergy);
   fChain->SetBranchAddress("jet2_HFEMEnergy", &jet2_HFEMEnergy, &b_jet2_HFEMEnergy);
   fChain->SetBranchAddress("jet2_HFHadronMultiplicity", &jet2_HFHadronMultiplicity, &b_jet2_HFHadronMultiplicity);
   fChain->SetBranchAddress("jet2_HFEMMultiplicity", &jet2_HFEMMultiplicity, &b_jet2_HFEMMultiplicity);
   fChain->SetBranchAddress("dijet_nj", &dijet_nj, &b_dijet_nj);
   fChain->SetBranchAddress("dipho_dijet_ptHjj", &dipho_dijet_ptHjj, &b_dipho_dijet_ptHjj);
   fChain->SetBranchAddress("n_jet_30", &n_jet_30, &b_n_jet_30);
   fChain->SetBranchAddress("dijet_jet1_RMS", &dijet_jet1_RMS, &b_dijet_jet1_RMS);
   fChain->SetBranchAddress("dijet_jet2_RMS", &dijet_jet2_RMS, &b_dijet_jet2_RMS);
   fChain->SetBranchAddress("dijet_jet1_QGL", &dijet_jet1_QGL, &b_dijet_jet1_QGL);
   fChain->SetBranchAddress("dijet_jet2_QGL", &dijet_jet2_QGL, &b_dijet_jet2_QGL);
   fChain->SetBranchAddress("dijet_jet1_pujid_mva", &dijet_jet1_pujid_mva, &b_dijet_jet1_pujid_mva);
   fChain->SetBranchAddress("dijet_jet2_pujid_mva", &dijet_jet2_pujid_mva, &b_dijet_jet2_pujid_mva);
   fChain->SetBranchAddress("dipho_pt", &dipho_pt, &b_dipho_pt);
   fChain->SetBranchAddress("dijet_pt", &dijet_pt, &b_dijet_pt);
   fChain->SetBranchAddress("gghMVA_n_rec_jets", &gghMVA_n_rec_jets, &b_gghMVA_n_rec_jets);
   fChain->SetBranchAddress("gghMVA_Mjj", &gghMVA_Mjj, &b_gghMVA_Mjj);
   fChain->SetBranchAddress("gghMVA_leadEta", &gghMVA_leadEta, &b_gghMVA_leadEta);
   fChain->SetBranchAddress("gghMVA_subleadEta", &gghMVA_subleadEta, &b_gghMVA_subleadEta);
   fChain->SetBranchAddress("gghMVA_subsubleadEta", &gghMVA_subsubleadEta, &b_gghMVA_subsubleadEta);
   fChain->SetBranchAddress("gghMVA_leadJPt", &gghMVA_leadJPt, &b_gghMVA_leadJPt);
   fChain->SetBranchAddress("gghMVA_SubleadJPt", &gghMVA_SubleadJPt, &b_gghMVA_SubleadJPt);
   fChain->SetBranchAddress("gghMVA_SubsubleadJPt", &gghMVA_SubsubleadJPt, &b_gghMVA_SubsubleadJPt);
   fChain->SetBranchAddress("gghMVA_leadPUMVA", &gghMVA_leadPUMVA, &b_gghMVA_leadPUMVA);
   fChain->SetBranchAddress("gghMVA_subleadPUMVA", &gghMVA_subleadPUMVA, &b_gghMVA_subleadPUMVA);
   fChain->SetBranchAddress("gghMVA_subsubleadPUMVA", &gghMVA_subsubleadPUMVA, &b_gghMVA_subsubleadPUMVA);
   fChain->SetBranchAddress("gghMVA_leadDeltaPhi", &gghMVA_leadDeltaPhi, &b_gghMVA_leadDeltaPhi);
   fChain->SetBranchAddress("gghMVA_subleadDeltaPhi", &gghMVA_subleadDeltaPhi, &b_gghMVA_subleadDeltaPhi);
   fChain->SetBranchAddress("gghMVA_subsubleadDeltaPhi", &gghMVA_subsubleadDeltaPhi, &b_gghMVA_subsubleadDeltaPhi);
   fChain->SetBranchAddress("gghMVA_leadDeltaEta", &gghMVA_leadDeltaEta, &b_gghMVA_leadDeltaEta);
   fChain->SetBranchAddress("gghMVA_subleadDeltaEta", &gghMVA_subleadDeltaEta, &b_gghMVA_subleadDeltaEta);
   fChain->SetBranchAddress("gghMVA_subsubleadDeltaEta", &gghMVA_subsubleadDeltaEta, &b_gghMVA_subsubleadDeltaEta);
   fChain->SetBranchAddress("ggHMVAResult_prob_0J_PTH_0_10", &ggHMVAResult_prob_0J_PTH_0_10, &b_ggHMVAResult_prob_0J_PTH_0_10);
   fChain->SetBranchAddress("ggHMVAResult_prob_0J_PTH_GT10", &ggHMVAResult_prob_0J_PTH_GT10, &b_ggHMVAResult_prob_0J_PTH_GT10);
   fChain->SetBranchAddress("ggHMVAResult_prob_1J_PTH_0_60", &ggHMVAResult_prob_1J_PTH_0_60, &b_ggHMVAResult_prob_1J_PTH_0_60);
   fChain->SetBranchAddress("ggHMVAResult_prob_1J_PTH_60_120", &ggHMVAResult_prob_1J_PTH_60_120, &b_ggHMVAResult_prob_1J_PTH_60_120);
   fChain->SetBranchAddress("ggHMVAResult_prob_1J_PTH_120_200", &ggHMVAResult_prob_1J_PTH_120_200, &b_ggHMVAResult_prob_1J_PTH_120_200);
   fChain->SetBranchAddress("ggHMVAResult_prob_GE2J_MJJ_0_350_PTH_0_60", &ggHMVAResult_prob_GE2J_MJJ_0_350_PTH_0_60, &b_ggHMVAResult_prob_GE2J_MJJ_0_350_PTH_0_60);
   fChain->SetBranchAddress("ggHMVAResult_prob_GE3J_MJJ_0_350_PTH_60_120", &ggHMVAResult_prob_GE3J_MJJ_0_350_PTH_60_120, &b_ggHMVAResult_prob_GE3J_MJJ_0_350_PTH_60_120);
   fChain->SetBranchAddress("ggHMVAResult_prob_GE2J_MJJ_0_350_PTH_120_200", &ggHMVAResult_prob_GE2J_MJJ_0_350_PTH_120_200, &b_ggHMVAResult_prob_GE2J_MJJ_0_350_PTH_120_200);
   fChain->SetBranchAddress("ggHMVAResult_prob_PTH_GT200", &ggHMVAResult_prob_PTH_GT200, &b_ggHMVAResult_prob_PTH_GT200);
   fChain->SetBranchAddress("vbfMvaResult_prob_bkg", &vbfMvaResult_prob_bkg, &b_vbfMvaResult_prob_bkg);
   fChain->SetBranchAddress("vbfMvaResult_prob_ggH", &vbfMvaResult_prob_ggH, &b_vbfMvaResult_prob_ggH);
   fChain->SetBranchAddress("vbfMvaResult_prob_VBF", &vbfMvaResult_prob_VBF, &b_vbfMvaResult_prob_VBF);
//    fChain->SetBranchAddress("dijet_minDRJetPho", &dijet_minDRJetPho, &b_dijet_minDRJetPho);
   fChain->SetBranchAddress("dijet_centrality_gg", &dijet_centrality_gg, &b_dijet_centrality_gg);
   fChain->SetBranchAddress("dijet_centrality_j3", &dijet_centrality_j3, &b_dijet_centrality_j3);
   fChain->SetBranchAddress("dijet_centrality_g", &dijet_centrality_g, &b_dijet_centrality_g);
   fChain->SetBranchAddress("cosThetaStar", &cosThetaStar, &b_cosThetaStar);
   fChain->SetBranchAddress("VH_had_mvascore", &VH_had_mvascore, &b_VH_had_mvascore);
   fChain->SetBranchAddress("dijet_jet1_match", &dijet_jet1_match, &b_dijet_jet1_match);
   fChain->SetBranchAddress("dijet_jet2_match", &dijet_jet2_match, &b_dijet_jet2_match);
   fChain->SetBranchAddress("prompt_pho_1", &prompt_pho_1, &b_prompt_pho_1);
   fChain->SetBranchAddress("prompt_pho_2", &prompt_pho_2, &b_prompt_pho_2);
   fChain->SetBranchAddress("HTXSstage0bin", &HTXSstage0bin, &b_HTXSstage0bin);
   fChain->SetBranchAddress("HTXSstage1bin", &HTXSstage1bin, &b_HTXSstage1bin);
//    fChain->SetBranchAddress("HTXSstage0bin", &HTXSstage0bin, &b_HTXSstage0bin);
//    fChain->SetBranchAddress("HTXSstage1bin", &HTXSstage1bin, &b_HTXSstage1bin);
   fChain->SetBranchAddress("HTXSstage1p1bin", &HTXSstage1p1bin, &b_HTXSstage1p1bin);
   fChain->SetBranchAddress("HTXSstage1p1binFine", &HTXSstage1p1binFine, &b_HTXSstage1p1binFine);
   fChain->SetBranchAddress("HTXSstage1p2bin", &HTXSstage1p2bin, &b_HTXSstage1p2bin);
   fChain->SetBranchAddress("HTXSstage1p2binFine", &HTXSstage1p2binFine, &b_HTXSstage1p2binFine);
   fChain->SetBranchAddress("HTXSnjets", &HTXSnjets, &b_HTXSnjets);
   fChain->SetBranchAddress("HTXSpTH", &HTXSpTH, &b_HTXSpTH);
   fChain->SetBranchAddress("HTXSpTV", &HTXSpTV, &b_HTXSpTV);
   fChain->SetBranchAddress("rho", &rho, &b_rho);
   fChain->SetBranchAddress("nvtx", &nvtx, &b_nvtx);
   fChain->SetBranchAddress("event", &event, &b_event);
   fChain->SetBranchAddress("lumi", &lumi, &b_lumi);
   fChain->SetBranchAddress("processIndex", &processIndex, &b_processIndex);
   fChain->SetBranchAddress("run", &run, &b_run);
//    fChain->SetBranchAddress("nvtx", &nvtx, &b_nvtx);
   fChain->SetBranchAddress("npu", &npu, &b_npu);
   fChain->SetBranchAddress("puweight", &puweight, &b_puweight);
   Notify();
}

Bool_t vbf_125_13TeV_GeneralDipho::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void vbf_125_13TeV_GeneralDipho::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t vbf_125_13TeV_GeneralDipho::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef vbf_125_13TeV_GeneralDipho_cxx
