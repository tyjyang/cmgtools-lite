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

Vec_f scalarToRVec(const float& var, const int& size) {

  Vec_f res(size,var); // initialize to 0
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
    if (deltaR(eta, phi, TrigObj_eta[jtrig], TrigObj_phi[jtrig]) < 0.3) return true;
  }
  return false;
  
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

// PU weights obtained from w-mass-13TeV/makePUweight.py on 21/02/2021 with
// python w-mass-13TeV/makePUweight.py --isPreVFP --cropHighWeight 5.0 (weigths larger than 5.0 cropped to 5.0)
static double _pileupWeights_2016UL_preVFP[100] = {0.017891859791232846, 0.25164270158539037, 0.7328704571739738, 0.8181251564460522, 0.8649406472442838, 0.5020328254845017, 0.2658104554235502, 0.25291785415248746, 0.43067445087530876, 0.5807375723700408, 0.7389514962080089, 0.8638650142763253, 0.9396099928679109, 0.9868029206345122, 1.0196162822425623, 1.0588761189937739, 1.094898081939306, 1.127276742027867, 1.149414736534094, 1.161488702664545, 1.1573768695552893, 1.140958539701416, 1.1235625472598643, 1.1081928436148043, 1.0893527771484866, 1.06450924232797, 1.025793115894688, 0.9746156509680088, 0.9186728763168767, 0.8616085351610866, 0.8020484997740417, 0.7445213151217449, 0.6912429884481729, 0.6402308026900916, 0.5893401409704421, 0.5405190531190912, 0.4935088544982941, 0.4467075390612075, 0.400872887447364, 0.3570491413946734, 0.31357528032369686, 0.2736781052187231, 0.2375036546212086, 0.2058278731154476, 0.18003241229641784, 0.15153863589379638, 0.12968940637544868, 0.10682015753856883, 0.0894323307837256, 0.07452733858657992, 0.06777639972884231, 0.07298084541281942, 0.08950648766178804, 0.12260107330067226, 0.19384713077284463, 0.2615885627432867, 0.3660971270237179, 0.4535015485650514, 0.6361173877919268, 0.9375635292635457, 1.209396722848853, 1.4619395406726488, 1.6359665710892473, 2.447616386237363, 1.8923192890766891, 2.4418060236344252, 5.0, 4.512706593240671, 5.0, 1.0, 1.0, 3.529756976636144, 5.0, 1.0, 1.0, 4.924548597100737, 1.0, 1.0, 1.0, 5.0, 0.7766971440272932, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0};

double puw_2016UL_preVFP(const Float_t& nTrueInt) {
  if (nTrueInt < 100.0)
    return _pileupWeights_2016UL_preVFP[static_cast<int>(nTrueInt)];
  else
    return 1.0;
}

// PU weights obtained from w-mass-13TeV/makePUweight.py on 21/02/2021 with
// python w-mass-13TeV/makePUweight.py --cropHighWeight 5.0 (weigths larger than 5.0 cropped to 5.0)
static double _pileupWeights_2016UL_postVFP[100] = {1.4121640872041221, 0.5898511462462072, 1.0186345612199057, 0.730896932534784, 0.6625906412478361, 0.373128877493089, 0.15202831320083082, 0.10810730190416742, 0.07922625145984823, 0.05412342031906427, 0.14377404335405167, 0.34244424125458106, 0.5009390258406541, 0.5858261792960938, 0.6415575114065856, 0.7087553686702128, 0.7632803248061507, 0.7922900511232938, 0.79791290590669, 0.8040725702550436, 0.837776272330411, 0.8992865658393996, 0.9622224459286195, 1.0101762330880721, 1.0505710052018327, 1.1002458975776779, 1.160152637272387, 1.2303177838850345, 1.3121316684401438, 1.399730491999252, 1.4803639433958555, 1.555545145192651, 1.62951843812034, 1.7006707124728095, 1.766838609765344, 1.8382375612180932, 1.9197947496434073, 2.008033481757116, 2.1037955476243497, 2.2073425542550345, 2.299520727318126, 2.39195308241725, 2.480611355126296, 2.570845679574736, 2.6856231319335375, 2.6895659148729645, 2.7170021459898392, 2.5991989807089486, 2.4413565237114865, 2.1175342531673946, 1.7286310371345444, 1.328305769847878, 0.9074697528382656, 0.5759813346926844, 0.3791280909476406, 0.2009268885008322, 0.10690147967822879, 0.0494866019105493, 0.025804167809444428, 0.014219209393874692, 0.006970784343273052, 0.003290065597384393, 0.0014923751883421411, 0.0009483513777169922, 0.0003281848478328599, 0.00019967886362757994, 0.0005765101556252217, 9.283233427021984e-05, 7.6901988053659e-05, 1.0, 1.0, 5.55519172322179e-06, 7.2808080984987366e-06, 1.0, 1.0, 5.607454368413037e-07, 1.0, 1.0, 1.0, 5.631509273014954e-08, 2.7009051699005925e-09, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0};
double puw_2016UL_postVFP(const Float_t& nTrueInt) {
  if (nTrueInt < 100.0)
    return _pileupWeights_2016UL_postVFP[static_cast<int>(nTrueInt)];
  else
    return 1.0;
}

// PU weights obtained from w-mass-13TeV/makePUweight.py on 21/02/2021 with
// python w-mass-13TeV/makePUweight.py --doInclusiveData --cropHighWeight 5.0
static double _pileupWeights_2016UL_all[100] = {0.6589867265263256, 0.407153005452541, 0.8642665336352487, 0.7780170875306136, 0.7718988813841754, 0.44276200466885046, 0.21349273305130528, 0.18633307976609584, 0.26907642709007246, 0.3385971761172566, 0.4652852764274064, 0.6241125653491979, 0.737906409006983, 0.8024313714506003, 0.8457825555557247, 0.8978884643903559, 0.9424182163604483, 0.9732478201709418, 0.9877920527908223, 0.9971465869206139, 1.0104225673337546, 1.029836293392211, 1.0493773854755748, 1.0631242087901378, 1.0715206820747583, 1.0809411741877695, 1.087572442599078, 1.0921890499092428, 1.099587626637935, 1.109040290182502, 1.113942074345191, 1.1174350144952392, 1.1226677579743707, 1.1278275354344371, 1.1307611124930197, 1.137217931323903, 1.1493238007450803, 1.1646147220938796, 1.1838871040556067, 1.207825328150971, 1.2267250990262277, 1.2476738483134555, 1.268898244611795, 1.293277472733083, 1.3321182940320204, 1.318539045280375, 1.319351570903305, 1.2528311287666474, 1.1708614071082524, 1.0139143678128297, 0.8314474831080699, 0.6501868955062319, 0.4656109815851464, 0.3310680803476131, 0.2790404453468702, 0.23369595525977435, 0.2469173888590991, 0.26773301406225, 0.3554916500496387, 0.513004184209515, 0.6565139488475376, 0.7912436231532053, 0.8844253754740058, 1.3226235140320606, 1.0223702497271776, 1.3191405473774638, 5.0, 2.437779286510053, 3.8756227608092186, 1.0, 1.0, 1.906755763065583, 4.76272763315208, 1.0, 1.0, 2.6602114002279516, 1.0, 1.0, 1.0, 4.269141635962611, 0.419567065322298, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0};
double puw_2016UL_all(const Float_t& nTrueInt) {
  if (nTrueInt < 100.0)
    return _pileupWeights_2016UL_all[static_cast<int>(nTrueInt)];
  else
    return 1.0;
}

double puw_2016UL_era(const Float_t& nTrueInt, const int& era) {
  if (nTrueInt < 100.0) {
    if      (era == 1) return _pileupWeights_2016UL_preVFP[static_cast<int>(nTrueInt)];
    else if (era == 2) return _pileupWeights_2016UL_postVFP[static_cast<int>(nTrueInt)];
    else               return _pileupWeights_2016UL_all[static_cast<int>(nTrueInt)];
  } else {
    return 1.0;
  }
}


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
