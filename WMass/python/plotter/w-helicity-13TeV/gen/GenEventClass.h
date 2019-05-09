//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue Apr 16 13:05:22 2019 by ROOT version 6.17/01
// from TTree Events/Events
// found on file: /Users/emanuele/Work/data/cms/wmass/GEN/powh-pythia/wp/GEN_100000.root
//////////////////////////////////////////////////////////

#ifndef GenEventClass_h
#define GenEventClass_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TH1F.h>
#include <TH3F.h>
#include <TString.h>
#include <TLorentzVector.h>

#include <vector>

// Header file for the classes stored in the TTree if any.
#include "TClonesArray.h"

class GenEventClass {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // output histograms
   TH1F *h_lpt, *h_leta, *h_prefsrlpt, *h_prefsrleta, *h_lptDressOverPreFSR;
   TH1F *h_wpt, *h_wy, *h_wmass, *h_genwpt, *h_genwy, *h_genwmass;
   TH1F *h_nfsr, *h_fsrdr_close, *h_fsrpt_close, *h_fsrdr_hard, *h_fsrpt_hard, *h_fsrptfrac_hard;

   TH3F *h3d_fsrdr_hard;

   TFile *outFile_;
   std::vector<TH1F*> histograms;
   
// Fixed size dimensions of array or collections stored in the TTree if any.
   static constexpr Int_t kMaxGenParticle = 1000;

   // Declaration of leaf types
 //baconhep::TGenEventInfo *GenEvtInfo;
   Int_t           id_1;
   Int_t           id_2;
   Double_t        x_1;
   Double_t        x_2;
   Double_t        xPDF_1;
   Double_t        xPDF_2;
   Double_t        scalePDF;
   Float_t         weight;
   //vector<double>  lheweight;
   Int_t           GenParticle_;
   Int_t           GenParticle_parent[kMaxGenParticle];   //[GenParticle_]
   Int_t           GenParticle_pdgId[kMaxGenParticle];   //[GenParticle_]
   Int_t           GenParticle_status[kMaxGenParticle];   //[GenParticle_]
   Float_t         GenParticle_pt[kMaxGenParticle];   //[GenParticle_]
   Float_t         GenParticle_eta[kMaxGenParticle];   //[GenParticle_]
   Float_t         GenParticle_phi[kMaxGenParticle];   //[GenParticle_]
   Float_t         GenParticle_mass[kMaxGenParticle];   //[GenParticle_]
   Float_t         GenParticle_y[kMaxGenParticle];   //[GenParticle_]

   // List of branches
   TBranch        *b_GenEvtInfo_id_1;   //!
   TBranch        *b_GenEvtInfo_id_2;   //!
   TBranch        *b_GenEvtInfo_x_1;   //!
   TBranch        *b_GenEvtInfo_x_2;   //!
   TBranch        *b_GenEvtInfo_xPDF_1;   //!
   TBranch        *b_GenEvtInfo_xPDF_2;   //!
   TBranch        *b_GenEvtInfo_scalePDF;   //!
   TBranch        *b_GenEvtInfo_weight;   //!
   //TBranch        *b_GenEvtInfo_lheweight;   //!
   TBranch        *b_GenParticle_;   //!
   TBranch        *b_GenParticle_parent;   //!
   TBranch        *b_GenParticle_pdgId;   //!
   TBranch        *b_GenParticle_status;   //!
   TBranch        *b_GenParticle_pt;   //!
   TBranch        *b_GenParticle_eta;   //!
   TBranch        *b_GenParticle_phi;   //!
   TBranch        *b_GenParticle_mass;   //!
   TBranch        *b_GenParticle_y;   //!

   GenEventClass(TTree *tree=0);
   virtual ~GenEventClass();
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop(int maxentries = -1);
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
   virtual void     bookHistograms();
   virtual void     writeHistograms();
   virtual void     setOutfile(TString);
   virtual void     setFlavor(int);
   TLorentzVector   getPreFSRLepton();
   std::vector<TLorentzVector> getDressedLeptons(float deltaR=0.1);
   std::vector<TLorentzVector> getFSRPhotons(TLorentzVector fourmom, float deltaR=0.1);
   TLorentzVector   getGenW();
   TLorentzVector   getNeutrino();
   TLorentzVector   getClosestPhoton(std::vector<TLorentzVector> photons,TLorentzVector fourmom);
   TLorentzVector   getHardestPhoton(std::vector<TLorentzVector> photons);
   bool             isPromptFinalStatePhoton(int index);
   bool             isPromptFinalStateLepton(int index);
   
   TString fOutfile;
   int fFlavor;
};

#endif

#ifdef GenEventClass_cxx
GenEventClass::GenEventClass(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("/Users/emanuele/Work/data/cms/wmass/GEN/powh-pythia/wp/GEN_100000.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("/Users/emanuele/Work/data/cms/wmass/GEN/powh-pythia/wp/GEN_100000.root");
      }
      f->GetObject("Events",tree);
   }
   Init(tree);
   fFlavor = 13;
}

GenEventClass::~GenEventClass()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t GenEventClass::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t GenEventClass::LoadTree(Long64_t entry)
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

void GenEventClass::Init(TTree *tree)
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

   fChain->SetBranchAddress("id_1", &id_1, &b_GenEvtInfo_id_1);
   fChain->SetBranchAddress("id_2", &id_2, &b_GenEvtInfo_id_2);
   fChain->SetBranchAddress("x_1", &x_1, &b_GenEvtInfo_x_1);
   fChain->SetBranchAddress("x_2", &x_2, &b_GenEvtInfo_x_2);
   fChain->SetBranchAddress("xPDF_1", &xPDF_1, &b_GenEvtInfo_xPDF_1);
   fChain->SetBranchAddress("xPDF_2", &xPDF_2, &b_GenEvtInfo_xPDF_2);
   fChain->SetBranchAddress("scalePDF", &scalePDF, &b_GenEvtInfo_scalePDF);
   fChain->SetBranchAddress("weight", &weight, &b_GenEvtInfo_weight);
   //fChain->SetBranchAddress("lheweight", &lheweight, &b_GenEvtInfo_lheweight);
   fChain->SetBranchAddress("GenParticle", &GenParticle_, &b_GenParticle_);
   fChain->SetBranchAddress("GenParticle.parent", GenParticle_parent, &b_GenParticle_parent);
   fChain->SetBranchAddress("GenParticle.pdgId", GenParticle_pdgId, &b_GenParticle_pdgId);
   fChain->SetBranchAddress("GenParticle.status", GenParticle_status, &b_GenParticle_status);
   fChain->SetBranchAddress("GenParticle.pt", GenParticle_pt, &b_GenParticle_pt);
   fChain->SetBranchAddress("GenParticle.eta", GenParticle_eta, &b_GenParticle_eta);
   fChain->SetBranchAddress("GenParticle.phi", GenParticle_phi, &b_GenParticle_phi);
   fChain->SetBranchAddress("GenParticle.mass", GenParticle_mass, &b_GenParticle_mass);
   fChain->SetBranchAddress("GenParticle.y", GenParticle_y, &b_GenParticle_y);
   Notify();
}

Bool_t GenEventClass::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void GenEventClass::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
#endif // #ifdef GenEventClass_cxx
