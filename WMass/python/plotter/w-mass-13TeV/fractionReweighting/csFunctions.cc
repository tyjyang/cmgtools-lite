#include <iostream>
#include <cstdlib> //as stdlib.h                 
#include <cstdio>
#include <map>
#include <string>
#include <cmath>
#include "TH2F.h"
#include "TVector2.h"
#include "TVector3.h"
#include "TLorentzVector.h"
#include "Math/GenVector/LorentzVector.h"
#include "Math/GenVector/PtEtaPhiM4D.h"
#include "ROOT/RVec.hxx"

typedef ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> > PtEtaPhiMVector;

bool getStatusFlag(Int_t flags, int index) { return ((flags >> index) & 1); }

std::pair<int, int> getIndexLepI(int ilep, int ngen, ROOT::VecOps::RVec<Int_t> pdgid, ROOT::VecOps::RVec<Int_t> midx, ROOT::VecOps::RVec<Int_t> status, ROOT::VecOps::RVec<Int_t> flags){

    //int nlep = 0;

    int status746_1 = -99;
    int status746_2 = -99;
    int other_1     = -99;
    int other_2     = -99;

    for (int i=0; i<ngen; ++i){
        if (!(TMath::Abs(pdgid[i]) == 13 || TMath::Abs(pdgid[i]) == 14) ) continue;
        if ( status[i] == 746 ) {
            if (status746_1 < 0) status746_1 = i;
            else                 status746_2 = i;
        }
        else if (status[i] == 1 && getStatusFlag(flags[i],8) ) {
            if (other_1 < 0) other_1 = i;
            else             other_2 = i;
        }
    }

    std::pair<int, int> returnValue = std::make_pair(-99, -99);

    if (status746_1 >= 0 ) {
        if (status746_2 >= 0) 
            returnValue = std::make_pair(status746_1, status746_2);
        else 
            returnValue = std::make_pair(status746_1, other_1); //not sure that this is 100% bullet proof
    }
    else 
        returnValue = std::make_pair(other_1, other_2);

    return returnValue;

    // nlep++;
    // if(nlep == ilep) return i;
    // return -999;
}

std::pair<TVector3, TVector3> csBoostedProtons(TLorentzVector dilepton) {
    float protonMass = 0.938272;
    float energy = 6500;
    int zsign = dilepton.Z() > 0 ? 1 : -1;
    TLorentzVector proton1(0., 0., zsign*energy, hypot(energy, protonMass));
    TLorentzVector proton2(0., 0., -1*zsign*energy, hypot(energy, protonMass));
    proton1.Boost(-1*dilepton.BoostVector());
    proton2.Boost(-1*dilepton.BoostVector());
    return std::make_pair<TVector3, TVector3>(proton1.Vect(), proton2.Vect());
}

const TVector3 csframe(TLorentzVector dilepton) {
    std::pair<TVector3, TVector3> protons = csBoostedProtons(dilepton);
    TVector3 csAxis = (protons.first.Unit() - protons.second.Unit()).Unit();
    return csAxis;
}

const TVector3 csframeY(TLorentzVector dilepton) {
    std::pair<TVector3, TVector3> protons = csBoostedProtons(dilepton);
    TVector3 csYAxis = protons.first.Unit().Cross(protons.second.Unit());
    return csYAxis.Unit();
}

const TVector3 csframeX(TLorentzVector dilepton) {
    TVector3 csAxis = csframe(dilepton);
    TVector3 csYAxis = csframeY(dilepton);
    TVector3 csXAxis = csYAxis.Cross(csAxis);
    return csXAxis.Unit();
}

float cosThetaCS(PtEtaPhiMVector lplus, PtEtaPhiMVector lminus) {
    PtEtaPhiMVector dilepton = lplus + lminus;
    TLorentzVector dilep(dilepton.X(), dilepton.Y(), dilepton.Z(), dilepton.T());
    TLorentzVector boostedLep(lplus.X(), lplus.Y(), lplus.Z(), lplus.T());
    boostedLep.Boost(-1*dilep.BoostVector());
    const TVector3 csFrame = csframe(dilep);
    return cos(boostedLep.Angle(csFrame));
}

float phiCS(PtEtaPhiMVector lplus, PtEtaPhiMVector lminus) {
    PtEtaPhiMVector dilepton = lplus + lminus;
    TLorentzVector dilep(dilepton.X(), dilepton.Y(), dilepton.Z(), dilepton.T());
    TLorentzVector boostedLep(lplus.X(), lplus.Y(), lplus.Z(), lplus.T());
    boostedLep.Boost(-1*dilep.BoostVector());
    const TVector3 csFrameX = csframeX(dilep);
    const TVector3 csFrameY = csframeY(dilep);
    float phi = atan2(boostedLep.Vect()*csFrameY, boostedLep.Vect()*csFrameX);
    return phi >= 0 ? phi : phi + 2*M_PI;
}
