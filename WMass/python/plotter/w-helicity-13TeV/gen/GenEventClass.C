#define GenEventClass_cxx
#include "GenEventClass.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <iostream>
#include <algorithm>


void GenEventClass::Loop(int maxentries)
{
//   In a ROOT session, you can do:
//      root> .L GenEventClass.C
//      root> GenEventClass t
//      root> t.GetEntry(12); // Fill t data members with entry number 12
//      root> t.Show();       // Show values of entry 12
//      root> t.Show(16);     // Read and show values of entry 16
//      root> t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch
   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntries();
   std::cout << "Chain has " << nentries << " events. Will run with maxevents = " << maxentries << std::endl;

   bookHistograms();
   
   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      if (maxentries > 0 && jentry >= maxentries) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;

      if(jentry%10000==0) std::cout << "Processing event " << jentry << std::endl;
      
      std::vector<TLorentzVector> dressedLeptonCollection = getDressedLeptons(TMath::Pi());
      TLorentzVector neutrino = getNeutrino();

      if(dressedLeptonCollection.size()){

        TLorentzVector preFSRLepton = getPreFSRLepton();
        h_prefsrlpt->Fill(preFSRLepton.Pt());
        h_prefsrleta->Fill(preFSRLepton.Eta());
        
        TLorentzVector dressedLepton = dressedLeptonCollection[0];
        h_lpt->Fill(dressedLepton.Pt());
        h_leta->Fill(dressedLepton.Eta());
        h_lptDressOverPreFSR->Fill(dressedLepton.Pt()/preFSRLepton.Pt());
        
        TLorentzVector recw = dressedLepton + neutrino;
        h_wpt->Fill(recw.Pt());
        h_wy->Fill(recw.Rapidity());
        h_wmass->Fill(recw.M());

        TLorentzVector genw = getGenW();
        h_genwpt->Fill(genw.Pt());
        h_genwy->Fill(genw.Rapidity());
        h_genwmass->Fill(genw.M());        

        std::vector<TLorentzVector> photons = getFSRPhotons(dressedLepton,TMath::Pi());
        h_nfsr->Fill(photons.size());
        if(photons.size()>0) {
          TLorentzVector fsr_close = getClosestPhoton(photons,dressedLepton);
          TLorentzVector fsr_hard = getHardestPhoton(photons);
          h_fsrpt_close->Fill(fsr_close.Pt());
          h_fsrdr_close->Fill(fsr_close.DeltaR(dressedLepton));
          h_fsrpt_hard->Fill(fsr_hard.Pt());
          h_fsrdr_hard->Fill(fsr_hard.DeltaR(dressedLepton));
          h_fsrptfrac_hard->Fill(fsr_hard.Pt()/preFSRLepton.Pt());
        } else {
          h_fsrpt_close->Fill(-1);
          h_fsrdr_close->Fill(-1);
          h_fsrpt_hard->Fill(-1);
          h_fsrdr_hard->Fill(-1);
          h_fsrptfrac_hard->Fill(-1);
        }
      }
   }
   writeHistograms();
   outFile_->Close();
}

TLorentzVector GenEventClass::getPreFSRLepton() {
  TLorentzVector lepton;
  for (int igp=0; igp<GenParticle_; ++igp) {
    if (abs(GenParticle_pdgId[igp]) != 13) continue;
    if (!isPromptFinalStateLepton(igp)) continue;
    lepton.SetPtEtaPhiM(GenParticle_pt[igp], GenParticle_eta[igp], GenParticle_phi[igp], std::max(GenParticle_mass[igp],float(0.)));
    break;
  }
  return lepton;
}

std::vector<TLorentzVector> GenEventClass::getDressedLeptons(float cone) {

  std::vector<TLorentzVector> leptons;

  for (int igp=0; igp<GenParticle_; ++igp) {
    if (abs(GenParticle_pdgId[igp]) != 13) continue;
    if (!isPromptFinalStateLepton(igp)) continue;
    TLorentzVector lepton;
    lepton.SetPtEtaPhiM(GenParticle_pt[igp], GenParticle_eta[igp], GenParticle_phi[igp], std::max(GenParticle_mass[igp],float(0.)));
    leptons.push_back(lepton);
  }

  sort(leptons.begin(), leptons.end(),
       [](const TLorentzVector& a, const TLorentzVector& b)
       {return a.Pt()>b.Pt();});
  
  if(leptons.size()==0) return leptons;

  for (int igp=0; igp<GenParticle_; ++igp) {
    if (abs(GenParticle_pdgId[igp]) != 22) continue;
    if (!isPromptFinalStatePhoton(igp)) continue;
    TLorentzVector tmp_photon;
    tmp_photon.SetPtEtaPhiM(GenParticle_pt[igp], GenParticle_eta[igp], GenParticle_phi[igp], std::max(GenParticle_mass[igp],float(0.)));

    for(int il=0; il<(int)leptons.size(); ++il) {
      if(leptons[il].DeltaR(tmp_photon) > cone) continue;
      leptons[il] = leptons[il]+tmp_photon;
      break;
    }
  }

  return leptons;
  
}

TLorentzVector GenEventClass::getGenW() {
  TLorentzVector W;
  for (int igp=0; igp<GenParticle_; ++igp) {
    if (abs(GenParticle_pdgId[igp]) != 24) continue;
    W.SetPtEtaPhiM(GenParticle_pt[igp], GenParticle_eta[igp], GenParticle_phi[igp], std::max(GenParticle_mass[igp],float(0.)));
    break;
  }
  return W;
}

TLorentzVector GenEventClass::getNeutrino() {

    std::vector<TLorentzVector> nus;

    for (int igp=0; igp<GenParticle_; ++igp) {
      if (abs(GenParticle_pdgId[igp]) != 14) continue;
      TLorentzVector nu;
      nu.SetPtEtaPhiM(GenParticle_pt[igp], GenParticle_eta[igp], GenParticle_phi[igp], std::max(GenParticle_mass[igp],float(0.)));
      nus.push_back(nu);
    }

    sort(nus.begin(), nus.end(),
         [](const TLorentzVector& a, const TLorentzVector& b)
         {return a.Pt()>b.Pt();});

    return nus[0];
}

std::vector<TLorentzVector> GenEventClass::getFSRPhotons(TLorentzVector fourmom, float cone) {

  std::vector<TLorentzVector> photons;

  for (int igp=0; igp<GenParticle_; ++igp) {
    if (abs(GenParticle_pdgId[igp]) != 22) continue;
    if (!isPromptFinalStatePhoton(igp)) continue;
    TLorentzVector tmp_photon;
    tmp_photon.SetPtEtaPhiM(GenParticle_pt[igp], GenParticle_eta[igp], GenParticle_phi[igp], std::max(GenParticle_mass[igp],float(0.)));
    float dr = fourmom.DeltaR(tmp_photon);
    if(dr < cone) {
      photons.push_back(tmp_photon);
    }
  }
  return photons;
}

TLorentzVector GenEventClass::getClosestPhoton(std::vector<TLorentzVector> photons, TLorentzVector fourmom) {
  if (photons.size()==0) return TLorentzVector(0,0,0,0);
  if (photons.size()==1) return photons[0];
  sort(photons.begin(), photons.end(),
       [fourmom](const TLorentzVector& a, const TLorentzVector& b)
       {return a.DeltaR(fourmom)<b.DeltaR(fourmom);});
  return photons[0];
}

TLorentzVector GenEventClass::getHardestPhoton(std::vector<TLorentzVector> photons) {
  if (photons.size()==0) return TLorentzVector(0,0,0,0);
  if (photons.size()==1) return photons[0];
  sort(photons.begin(), photons.end(),
       [](const TLorentzVector& a, const TLorentzVector& b)
       {return a.Pt()>b.Pt();});
  return photons[0];
}

bool GenEventClass::isPromptFinalStatePhoton(int index) {
  int motherId = GenParticle_pdgId[GenParticle_parent[index]];
  // with photos the photon can be attached to the lepton or the W/Z
  return (abs(motherId)==13 || abs(motherId)==24);
}

vbool GenEventClass::isPromptFinalStateLepton(int index) {
  int motherId = GenParticle_pdgId[GenParticle_parent[index]];
  // with photos the photon can be attached to the lepton or the W/Z
  return (abs(motherId)==24);
}

void GenEventClass::bookHistograms() {

  outFile_ = new TFile(fOutfile, "RECREATE");    
  outFile_->cd();
  
  h_lpt  = new TH1F("lpt","lepton pt",1000,0,250);              histograms.push_back(h_lpt);
  h_leta = new TH1F("leta","lepton eta",1000,-5.,5.);           histograms.push_back(h_leta);
  h_prefsrlpt  = new TH1F("prefsrlpt","pre-FSR lepton pt",1000,0,250);      histograms.push_back(h_prefsrlpt);
  h_prefsrleta = new TH1F("prefsrleta","pre-FSR lepton eta",1000,-5.,5.);   histograms.push_back(h_prefsrleta);

  h_lptDressOverPreFSR  = new TH1F("lptDressOverPreFSR","dressed lepton pt over preFSR",1000,0.5,1.5);   histograms.push_back(h_lptDressOverPreFSR);

  h_wpt      = new TH1F("wpt","W pt",1000,0,150);              histograms.push_back(h_wpt);  
  h_wy       =  new TH1F("wy","W rapidity",1000,-6,6);         histograms.push_back(h_wy);
  h_wmass    =  new TH1F("wmass","W mass",1000,0,180);         histograms.push_back(h_wmass);
  h_genwpt   = new TH1F("genwpt","gen W pt",1000,0,150);       histograms.push_back(h_genwpt);  
  h_genwy    =  new TH1F("genwy","gen W rapidity",1000,-6,6);  histograms.push_back(h_genwy);
  h_genwmass =  new TH1F("genwmass","gen W mass",1000,0,180);  histograms.push_back(h_genwmass);

  h_nfsr = new TH1F("nfsr","n FSR photons",10,0,10);                     histograms.push_back(h_nfsr);
  h_fsrpt_close = new TH1F("fsrpt_close","FSR pt",5000,0,10);            histograms.push_back(h_fsrpt_close);
  h_fsrdr_close = new TH1F("fsrdr_close","FSR pt",5000,0,TMath::Pi());   histograms.push_back(h_fsrdr_close);
  h_fsrpt_hard = new TH1F("fsrpt_hard","FSR pt",5000,0,10);              histograms.push_back(h_fsrpt_hard);
  h_fsrdr_hard = new TH1F("fsrdr_hard","FSR pt",5000,0,TMath::Pi());     histograms.push_back(h_fsrdr_hard);

  h_fsrptfrac_hard = new TH1F("fsrptfrac_hard","fraction FSR pt / preFSR lepton pt",5000,0,0.1);      histograms.push_back(h_fsrptfrac_hard);


}

void GenEventClass::writeHistograms() {
  outFile_->cd();
  for (int ih=0; ih<(int)histograms.size(); ++ih) {
    histograms[ih]->Write();
  }
}

void GenEventClass::setOutfile(TString outfilepath){
  fOutfile = outfilepath;
}
