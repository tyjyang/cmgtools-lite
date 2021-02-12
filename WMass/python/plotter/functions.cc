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
#include "TH2F.h"
#include "TVector2.h"
#include "Math/GenVector/LorentzVector.h"
#include "Math/GenVector/PtEtaPhiM4D.h"

#include <ROOT/RVec.hxx>

using namespace std;
using Vec_t = ROOT::VecOps::RVec<float>;
typedef ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> > PtEtaPhiMVector;

//// UTILITY FUNCTIONS NOT IN TFORMULA ALREADY

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
    float deta = std::abs(eta1-eta2);
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

float invariantmass(const Vec_t& pt, const Vec_t& eta, const Vec_t& phi, const Vec_t& m) {
  PtEtaPhiMVector p1(pt[0], eta[0], phi[0], m[0]);
  PtEtaPhiMVector p2(pt[1], eta[1], phi[1], m[1]);
  return (p1 + p2).mass();
}

float rapidity(const Vec_t& pt, const Vec_t& eta, const Vec_t& phi, const Vec_t& m) {
  //typedef ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> > PtEtaPhiMVector;
  PtEtaPhiMVector p1(pt[0], eta[0], phi[0], m[0]);;
  PtEtaPhiMVector p2(pt[1], eta[1], phi[1], m[1]);
  return (p1 + p2).Rapidity();
}

float transversemomentum(const Vec_t& pt, const Vec_t& phi) {
  float phidiff = phi[1]-phi[0];
  return hypot(pt[0] + pt[1] * std::cos(phidiff), pt[1]*std::sin(phidiff));
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

float _puw_2016UL[100] = {0.6566053285329763, 0.41144868907530213, 0.8554973866769888, 0.7816231832979451, 0.7620601163880255, 0.44294370235808755, 0.21337916188839076, 0.1861679518260925, 0.26982560763525065, 0.3390190538106666, 0.4661745977152918, 0.624071351258268, 0.7379568044102374, 0.8023578601527022, 0.8459955872327743, 0.8979482340401523, 0.9426709044822664, 0.9729954527622069, 0.9878074496383361, 0.9966455425761672, 1.0099691798677684, 1.0296749744202007, 1.049521661450095, 1.0631150862840977, 1.0718721179591304, 1.0804003995042164, 1.087603306855375, 1.0922484639576304, 1.1002285047846572, 1.1089492682437438, 1.114592478539555, 1.1174652378549064, 1.1225708382402302, 1.1279537608970822, 1.1302795330092479, 1.1367432184190533, 1.1485613627731475, 1.1640128752737666, 1.1833126030571468, 1.20703482909623, 1.2292592420842194, 1.245259922261253, 1.2699117167626013, 1.2899951098429536, 1.348435424711583, 1.3191011175054879, 1.3363059630911416, 1.2218986839950168, 1.1565777584939498, 1.015027055750588, 0.8105109874154754, 0.6490468375922713, 0.4510927481328359, 0.3160235084562487, 0.2899825531248701, 0.21472738198393557, 0.241741259020553, 0.22402553994211524, 0.3195568610430816, 0.5267473763363504, 0.5782504216992813, 0.6471873635674419, 1.075983188715943, 1.3688192568416286, 0.6767765558933262, 0.9204306145749046, 1.0, 1.0, 1.0, 1.0, 1.0, 1.4191352148340897, 1.0, 1.0, 1.0, 1.2691713536696225, 1.0, 1.0, 1.0, 1.8331043321564422, 0.3753243775602727, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0};
float puw_2016UL(int nTrueInt) { if (nTrueInt<100) return _puw_2016UL[nTrueInt]; else return 0.; }


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
