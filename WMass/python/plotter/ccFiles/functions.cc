#ifndef FUNCTIONS_H
#define FUNCTIONS_H

// #include <stdio.h>
// #include <stdlib.h>
#include <iostream>
#include <cstdlib> //as stdlib.h                 
#include <cstdio>
#include <map>
#include <string>
#include <cmath>
#include "TROOT.h"
#include "TVector2.h"
#include "Math/GenVector/LorentzVector.h"
#include "Math/GenVector/PtEtaPhiM4D.h"

#include "defines.h"

using namespace std;

typedef ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> > PtEtaPhiMVector;

//// UTILITY FUNCTIONS NOT IN TFORMULA ALREADY

// ROOT::RDF::RResultPtr<double> getRDFcolumnSum(const ROOT::RDataFrame& rdf, const std::string& column) {
//   return rdf::Sum<double>(column);
// }

double genWeightLargeClipped(const double& wgt, const double& max) {

  // max is already a positive number, no need to check or take absolute value here
  return static_cast<double>(std::copysign(1.0,wgt) * std::min<double>(std::abs(wgt),max));
  
}

double genWeightLargeRemoved(const double& wgt, const double& max) {

  // max is already a positive number, no need to check or take absolute value here
  return (std::abs(wgt) < max) ? wgt : 0.0;
  
}

Vec_i indices(const Vec_f& vec, const int& start = 0) {
    Vec_i res(vec.size(), 0);
    std::iota(std::begin(res), std::end(res), start);
    return res;
}

Vec_i indices(const int& size, const int& start = 0) {
    Vec_i res(size, 0);
    std::iota(std::begin(res), std::end(res), start);
    return res;
}


Vec_f scalarToRVec(const float& var, const int& size) {

  Vec_f res(size,var); // initialize to 0
  return res;
  
}

Vec_f scalarToRVec(const float var, const Vec_d& size) {
  Vec_f res(size.size(), var); // initialize to 0
  return res;
}

Vec_f scalarToRVec(const float& var, const Vec_f& size) {

  Vec_f res(size.size(),var); // initialize to 0
  return res;
  
}

Vec_f scalarToRVec(const float& var, const Vec_i& size) {

  Vec_f res(size.size(),var); // initialize to 0
  return res;
  
}

float deltaPhi(float phi1, float phi2) {
    float result = phi1 - phi2;
    while (result > float(M_PI)) result -= float(2*M_PI);
    while (result <= -float(M_PI)) result += float(2*M_PI);
    return result;
}

float if3(bool cond, float iftrue, float iffalse) {
    return cond ? iftrue : iffalse;
}

float deltaR2(float eta1, float phi1, float eta2, float phi2) {
    float deta = eta1-eta2;
    float dphi = deltaPhi(phi1,phi2);
    return deta*deta + dphi*dphi;
}
float deltaR(float eta1, float phi1, float eta2, float phi2) {
    return std::sqrt(deltaR2(eta1,phi1,eta2,phi2));
}

float Hypot(float x, float y) {
  return hypot(x,y);
}

float pt_2(float pt1, float phi1, float pt2, float phi2) {
    phi2 -= phi1;
    return hypot(pt1 + pt2 * std::cos(phi2), pt2*std::sin(phi2));
}

float rapidity_2(float pt1, float eta1, float phi1, float m1, float pt2, float eta2, float phi2, float m2) {
    //typedef ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> > PtEtaPhiMVector;
    PtEtaPhiMVector p41(pt1,eta1,phi1,m1);
    PtEtaPhiMVector p42(pt2,eta2,phi2,m2);
    return (p41+p42).Rapidity();
}

float mt_2(float pt1, float phi1, float pt2, float phi2) {
    return std::sqrt(2*pt1*pt2*(1-std::cos(phi1-phi2)));
}

float mass_2_ene(float ene1, float eta1, float phi1, float m1, float ene2, float eta2, float phi2, float m2) {
  //typedef ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> > PtEtaPhiMVector;
    PtEtaPhiMVector unitp41(1.0,eta1,phi1,m1);
    PtEtaPhiMVector unitp42(1.0,eta2,phi2,m2);
    double theta1 = unitp41.Theta();
    double theta2 = unitp42.Theta();
    double pt1 = ene1*fabs(sin(theta1));
    double pt2 = ene2*fabs(sin(theta2));
    PtEtaPhiMVector p41(pt1,eta1,phi1,m1);
    PtEtaPhiMVector p42(pt2,eta2,phi2,m2);
    return (p41+p42).M();
}

float mass_2(float pt1, float eta1, float phi1, float m1, float pt2, float eta2, float phi2, float m2) {
    PtEtaPhiMVector p41(pt1,eta1,phi1,m1);
    PtEtaPhiMVector p42(pt2,eta2,phi2,m2);
    return (p41+p42).M();
}

float invariantmass(const Vec_f& pt, const Vec_f& eta, const Vec_f& phi, const Vec_f& m) {
  PtEtaPhiMVector p1(pt[0], eta[0], phi[0], m[0]);
  PtEtaPhiMVector p2(pt[1], eta[1], phi[1], m[1]);
  return (p1 + p2).mass();
}

float rapidity(const Vec_f& pt, const Vec_f& eta, const Vec_f& phi, const Vec_f& m) {
  //typedef ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> > PtEtaPhiMVector;
  PtEtaPhiMVector p1(pt[0], eta[0], phi[0], m[0]);;
  PtEtaPhiMVector p2(pt[1], eta[1], phi[1], m[1]);
  return (p1 + p2).Rapidity();
}

float transversemomentum(const Vec_f& pt, const Vec_f& phi) {
  float phidiff = phi[1] - phi[0];
  return hypot(pt[0] + pt[1] * std::cos(phidiff), pt[1]*std::sin(phidiff));
}

Vec_b chargedParticleByEventParity(const ULong64_t& event, const Vec_b& plus, const Vec_b& minus) {

  if (isOddEvent(event)) return plus;
  else                   return minus;
  
}

Vec_b goodMuonTriggerCandidate(const Vec_i& TrigObj_id, const Vec_f& TrigObj_pt, const Vec_f& TrigObj_l1pt, const Vec_f& TrigObj_l2pt, const Vec_i& TrigObj_filterBits) {

   Vec_b res(TrigObj_id.size(),false); // initialize to 0   
   for (unsigned int i = 0; i < res.size(); ++i) {
       if (TrigObj_id[i]  != 13 ) continue;
       if (TrigObj_pt[i]   < 24.) continue;
       if (TrigObj_l1pt[i] < 22.) continue;
       if (! (( TrigObj_filterBits[i] & 8) || (TrigObj_l2pt[i] > 10. && (TrigObj_filterBits[i] & 2) )) ) continue;
       res[i] = true;
   }
   // res will be goodTrigObjs in RDF
   // e.g. RDF::Define("goodTrigObjs","goodMuonTriggerCandidate(TrigObj_id,TrigObj_pt,TrigObj_l1pt,TrigObj_l2pt,TrigObj_filterBits)")
   return res;
}

Vec_b hasTriggerMatch(const Vec_f& eta, const Vec_f& phi, const Vec_f& TrigObj_eta, const Vec_f& TrigObj_phi) {

   Vec_b res(eta.size(),false); // initialize to 0
   for (unsigned int i = 0; i < res.size(); ++i) {
      for (unsigned int jtrig = 0; jtrig < res.size(); ++jtrig) {
	  // use deltaR*deltaR < 0.3*0.3, to be faster 
          if (deltaR2(eta[i], phi[i], TrigObj_eta[jtrig], TrigObj_phi[jtrig]) < 0.09) {
              res[i] = true;
              break; // exit loop on trigger objects, and go to next muon
          }
      }
   }
   // res will be triggerMatchedMuons in RDF, like
   // RDF::Define("triggerMatchedMuons","hasTriggerMatch(Muon_eta,Muon_phi,TrigObj_eta[goodTrigObjs],TrigObj_phi[goodTrigObjs])")
   return res;

}

bool hasTriggerMatch(const float& eta, const float& phi, const Vec_f& TrigObj_eta, const Vec_f& TrigObj_phi) {

  for (unsigned int jtrig = 0; jtrig < TrigObj_eta.size(); ++jtrig) {
    if (deltaR2(eta, phi, TrigObj_eta[jtrig], TrigObj_phi[jtrig]) < 0.09) return true;
  }
  return false;
  
}

Vec_b cleanJetsFromMuons(const Vec_f& Jet_eta, const Vec_f& Jet_phi, const Vec_f& Muon_eta, const Vec_f& Muon_phi) {

   Vec_b res(Jet_eta.size(), true); // initialize to true and set to false whenever the jet overlaps with a muon
   for (unsigned int ij = 0; ij < res.size(); ++ij) {
     for (unsigned int im = 0; im < Muon_eta.size(); ++im) {
       if (deltaR2(Jet_eta[ij], Jet_phi[ij], Muon_eta[im], Muon_phi[im]) < 0.16) { // cone DR = 0.4
	 res[ij] = false;
	 break;
       }
     }
   }

   return res;
}

Vec_b cleanJetsFromLeptons(const Vec_f& Jet_eta, const Vec_f& Jet_phi, const Vec_f& Muon_eta, const Vec_f& Muon_phi, const Vec_f& Electron_eta, const Vec_f& Electron_phi) {

   Vec_b res(Jet_eta.size(), true); // initialize to true and set to false whenever the jet overlaps with a muon

   for (unsigned int ij = 0; ij < res.size(); ++ij) {

     for (unsigned int im = 0; im < Muon_eta.size(); ++im) {
       if (deltaR2(Jet_eta[ij], Jet_phi[ij], Muon_eta[im], Muon_phi[im]) < 0.16) { // cone DR = 0.4
	 res[ij] = false;
	 break;
       }
     }

     if (res[ij]) {
       for (unsigned int ie = 0; ie < Electron_eta.size(); ++ie) {
	 if (deltaR2(Jet_eta[ij], Jet_phi[ij], Electron_eta[ie], Electron_phi[ie]) < 0.16) { // cone DR = 0.4
	   res[ij] = false;
	   break;
	 }
       }
     }

   }

   return res;
}


Vec_f absoluteValue(const Vec_f& val) {

  Vec_f res(val.size(),0.0); // initialize to 0
  for (unsigned int i = 0; i < res.size(); ++i) {
    res[i] = std::abs(val[i]);
  }
  return res;

}

Vec_i absoluteValue(const Vec_i& val) {

  Vec_i res(val.size(),0.0); // initialize to 0
  for (unsigned int i = 0; i < res.size(); ++i) {
    res[i] = std::abs(val[i]);
  }
  return res;

}

#include "TRandom3.h"
TRandom3 *randy = NULL;

float mass_2_smeared(float pt1, float eta1, float phi1, float m1, float pt2, float eta2, float phi2, float m2, int isData) {
    float finalpt1;
    float finalpt2;
    if (isData){
        finalpt1 = pt1;
        finalpt2 = pt2;
    }
    else {
        if (!randy) randy = new TRandom3(42);
        finalpt1 = pt1*(1.+0.34/pt1/1.4142*randy->Gaus(0.,1.) ) - 0.06957/pt1;
        finalpt2 = pt2*(1.+0.34/pt2/1.4142*randy->Gaus(0.,1.) ) - 0.06957/pt2;
    }

//    std::cout << "initial pT1: " << pt1 << " corrected pT1: " << finalpt1 << std::endl;
//    std::cout << "initial pT2: " << pt2 << " corrected pT2: " << finalpt2 << std::endl;

    typedef ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> > PtEtaPhiMVector;
    PtEtaPhiMVector p41(finalpt1,eta1,phi1,m1);
    PtEtaPhiMVector p42(finalpt2,eta2,phi2,m2);

    PtEtaPhiMVector p43(pt1,eta1,phi1,m1);
    PtEtaPhiMVector p44(pt2,eta2,phi2,m2);

    float finalmll = (p41+p42).M();
    float initialm = (p43+p44).M();

    //std::cout << "initial mll " << initialm << " final mll " << finalmll << std::endl;


    return (p41+p42).M();
}

float eta_2(float pt1, float eta1, float phi1, float m1, float pt2, float eta2, float phi2, float m2) {
    typedef ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> > PtEtaPhiMVector;
    PtEtaPhiMVector p41(pt1,eta1,phi1,m1);
    PtEtaPhiMVector p42(pt2,eta2,phi2,m2);
    return (p41+p42).Eta();
}

float pt_3(float pt1, float phi1, float pt2, float phi2, float pt3, float phi3) {
    phi2 -= phi1;
    phi3 -= phi1;
    return hypot(pt1 + pt2 * std::cos(phi2) + pt3 * std::cos(phi3), pt2*std::sin(phi2) + pt3*std::sin(phi3));
}

float mass_3(float pt1, float eta1, float phi1, float m1, float pt2, float eta2, float phi2, float m2, float pt3, float eta3, float phi3, float m3) {
    typedef ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> > PtEtaPhiMVector;
    PtEtaPhiMVector p41(pt1,eta1,phi1,m1);
    PtEtaPhiMVector p42(pt2,eta2,phi2,m2);
    PtEtaPhiMVector p43(pt3,eta3,phi3,m3);
    return (p41+p42+p43).M();
}

float pt_4(float pt1, float phi1, float pt2, float phi2, float pt3, float phi3, float pt4, float phi4) {
    phi2 -= phi1;
    phi3 -= phi1;
    phi4 -= phi1;
    return hypot(pt1 + pt2 * std::cos(phi2) + pt3 * std::cos(phi3) + pt4 * std::cos(phi4), pt2*std::sin(phi2) + pt3*std::sin(phi3) + pt4*std::sin(phi4));
}
 
float mass_4(float pt1, float eta1, float phi1, float m1, float pt2, float eta2, float phi2, float m2, float pt3, float eta3, float phi3, float m3, float pt4, float eta4, float phi4, float m4) {
    typedef ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> > PtEtaPhiMVector;
    PtEtaPhiMVector p41(pt1,eta1,phi1,m1);
    PtEtaPhiMVector p42(pt2,eta2,phi2,m2);
    PtEtaPhiMVector p43(pt3,eta3,phi3,m3);
    PtEtaPhiMVector p44(pt4,eta4,phi4,m4);
    return (p41+p42+p43+p44).M();
}

float mt_llv(float ptl1, float phil1, float ptl2, float phil2, float ptv, float phiv) {
    float px = ptl1*std::cos(phil1) + ptl2*std::cos(phil2) + ptv*std::cos(phiv);
    float py = ptl1*std::sin(phil1) + ptl2*std::sin(phil2) + ptv*std::sin(phiv);
    float ht = ptl1+ptl2+ptv;
    return std::sqrt(std::max(0.f, ht*ht - px*px - py*py));
}

float mt_lllv(float ptl1, float phil1, float ptl2, float phil2, float ptl3, float phil3, float ptv, float phiv) {
    float px = ptl1*std::cos(phil1) + ptl2*std::cos(phil2) + ptl3*std::cos(phil3) + ptv*std::cos(phiv);
    float py = ptl1*std::sin(phil1) + ptl2*std::sin(phil2) + ptl3*std::sin(phil3) + ptv*std::sin(phiv);
    float ht = ptl1+ptl2+ptl3+ptv;
    return std::sqrt(std::max(0.f, ht*ht - px*px - py*py));
}


float mtw_wz3l(float pt1, float eta1, float phi1, float m1, float pt2, float eta2, float phi2, float m2, float pt3, float eta3, float phi3, float m3, float mZ1, float met, float metphi) 
{
    if (abs(mZ1 - mass_2(pt1,eta1,phi1,m1,pt2,eta2,phi2,m2)) < 0.01) return mt_2(pt3,phi3,met,metphi);
    if (abs(mZ1 - mass_2(pt1,eta1,phi1,m1,pt3,eta3,phi3,m3)) < 0.01) return mt_2(pt2,phi2,met,metphi);
    if (abs(mZ1 - mass_2(pt2,eta2,phi2,m2,pt3,eta3,phi3,m3)) < 0.01) return mt_2(pt1,phi1,met,metphi);
    return 0;
}

float mt_lu_cart(float lep_pt, float lep_phi, float u_x, float u_y)
{
    float lep_px = lep_pt*std::cos(lep_phi), lep_py = lep_pt*std::sin(lep_phi);
    float u = hypot(u_x,u_y);
    float uDotLep = u_x*lep_px + u_y*lep_py;
    return sqrt(2*lep_pt*sqrt(u*u+lep_pt*lep_pt+2*uDotLep) + 2*uDotLep + 2*lep_pt*lep_pt);
}

float u1_2(float met_pt, float met_phi, float ref_pt, float ref_phi) 
{
    float met_px = met_pt*std::cos(met_phi), met_py = met_pt*std::sin(met_phi);
    float ref_px = ref_pt*std::cos(ref_phi), ref_py = ref_pt*std::sin(ref_phi);
    float ux = - met_px + ref_px, uy = - met_py + ref_py;
    return (ux*ref_px + uy*ref_py)/ref_pt;
}
float u2_2(float met_pt, float met_phi, float ref_pt, float ref_phi)
{
    float met_px = met_pt*std::cos(met_phi), met_py = met_pt*std::sin(met_phi);
    float ref_px = ref_pt*std::cos(ref_phi), ref_py = ref_pt*std::sin(ref_phi);
    float ux = - met_px + ref_px, uy = - met_py + ref_py;
    return (ux*ref_py - uy*ref_px)/ref_pt;
}

float met_cal(float met_pt, float met_phi, float lep_pt, float lep_phi, float u_coeff, float u_syst)
{
    float met_px = met_pt*std::cos(met_phi), met_py = met_pt*std::sin(met_phi);
    float lep_px = lep_pt*std::cos(lep_phi), lep_py = lep_pt*std::sin(lep_phi);
    float ux = met_px + lep_px, uy = met_py + lep_py;
    float metcal_px = - u_coeff*ux*(1+u_syst) - lep_px, metcal_py = - u_coeff*uy*(1+u_syst) - lep_py;
    return hypot(metcal_px,metcal_py);
}

float _puw2016_nTrueInt_36fb[100] = {0.3505407355600995, 0.8996968628890968, 1.100322319466069, 0.9562526765089195, 1.0366251229154624, 1.0713954619016586, 0.7593488199769544, 0.47490309461978414, 0.7059895997695581, 0.8447022252423783, 0.9169159386164522, 1.0248924033173097, 1.0848877947714115, 1.1350984224561655, 1.1589888429954602, 1.169048420382294, 1.1650383018054549, 1.1507200023444994, 1.1152571438041776, 1.0739529436969637, 1.0458014000030829, 1.032500407707141, 1.0391236062781293, 1.041283620738903, 1.0412963370894526, 1.0558823002770783, 1.073481674823461, 1.0887053272606795, 1.1041701696801014, 1.123218903738397, 1.1157169321377927, 1.1052520327174429, 1.0697489590429388, 1.0144652740600584, 0.9402657069968621, 0.857142825520793, 0.7527112615290031, 0.6420618248685722, 0.5324755829715156, 0.4306470627563325, 0.33289171600176093, 0.24686361729094983, 0.17781595237914027, 0.12404411884835284, 0.08487088505600057, 0.056447805688061216, 0.03540829360547507, 0.022412461576677457, 0.013970541270658443, 0.008587896629717911, 0.004986410514292661, 0.00305102303701641, 0.001832072556146534, 0.0011570757619737708, 0.0008992999249003301, 0.0008241241729452477, 0.0008825716073180279, 0.001187003960081393, 0.0016454104270429153, 0.0022514113879764414, 0.003683196037880878, 0.005456695951503178, 0.006165248770884191, 0.007552675218762607, 0.008525338219226993, 0.008654690499815343, 0.006289068906974821, 0.00652551838513972, 0.005139581024893171, 0.005115751962934923, 0.004182527768384693, 0.004317593022028565, 0.0035749335962533355, 0.003773660372937113, 0.002618732319396435, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0};
float puw2016_nTrueInt_36fb(int nTrueInt) { if (nTrueInt<100) return _puw2016_nTrueInt_36fb[nTrueInt]; else return 0; }


//==================================================


bool valueInsideRange(float value, float low, float high) {

  if (value > low and value < high) return true;
  else                              return false;

}

//==================================================

bool valueOutsideRange(float value, float low, float high) {

  if (value < low or value > high) return true;
  else                             return false;

}

//==================================================

float varLepPlusFromPair(float var1, int pdgid1, float var2, int pdgid2) {
  
  // pdg ID > 0 for particles, i.e. negative leptons
  // check that two leptons have opposite charge, return dummy value if not
  if (pdgid1*pdgid2 > 0) return -9999.0;

  if (pdgid1 > 0) return var2;
  else            return var1;

}

//==================================================
float varLepMinusFromPair(float var1, int pdgid1, float var2, int pdgid2) {
  
  // pdg ID > 0 for particles, i.e. negative leptons
  // check that two leptons have opposite charge, return dummy value if not
  if (pdgid1*pdgid2 > 0) return -9999.0;

  if (pdgid1 > 0) return var1;
  else            return var2;

}

//==================================================
float varChargedLepFromPair(int requestedCharge, float var1, int pdgid1, float var2, int pdgid2) {
  
  // requestedCharge must be > 0 for positive charge and < 0 otherwise

  // pdg ID > 0 for particles, i.e. negative leptons
  // check that two leptons have opposite charge, return dummy value if not
  if (pdgid1*pdgid2 > 0)   return -9999.0;
  if (requestedCharge > 0) return (pdgid1 < 0) ? var1 : var2;
  else                     return (pdgid1 > 0) ? var1 : var2;

}

//==================================================
TRandom3 *randy_v2 = NULL;
double randomVarFromPair(float var1, float var2) {
  
  // pdg ID > 0 for particles, i.e. negative leptons
  // check that two leptons have opposite charge, return 0 if not
  if (!randy_v2) randy_v2 = new TRandom3(0);
  if (randy_v2->Rndm() > 0.5) return var1;
  else                        return var2;

}

double ptDiffCharge(float pt1, int charge1, float pt2, int charge2) {

  if (charge1 < 0) return pt1-pt2;
  else             return pt2-pt1;

}



void functions() {}

#endif
