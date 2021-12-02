#ifndef systHistHelpers
#define systHistHelpers

#include "TH2.h"
#include "TH2D.h"
#include "TH3.h"
#include <vector>

// Fix this so it makes the copy here
template <typename T>
T mirrorHist(T& nominal, T& alternate, T& mirrored) {
    for (size_t b = 0; b <= nominal.GetNcells(); b++) {
        double y0 = nominal.GetBinContent(b);
        double yA = alternate.GetBinContent(b);
        double yM = (yA > 0 && y0 > 0) ? y0*y0/yA : 0;
        //std::cout << "Here yM is " << yM << " y0 is " << y0 << " ratio is " << y0/yA << std::endl;
        mirrored.SetBinContent(b, yM);
        mirrored.SetBinError(b, alternate.GetBinError(b)*y0/yA);
    }
    return mirrored;
}

TH2D projectTH2FromTH3(TH3& hist3D, const char* name, size_t binStart, size_t binEnd=0) {
    if (binEnd == 0)
        binEnd = binStart;
    hist3D.GetZaxis()->SetRange(binStart, binEnd);
    //Order yx matters to have consistent axes!
    TH2D hist2D = *static_cast<TH2D*>(hist3D.Project3D("yxe")); 
    hist2D.SetName(name);
    return hist2D;
}

std::vector<TH2D> allProjectedSystHists(TH3& hist3d, std::string name, int binStart=1, size_t binEnd=0) {
    if (binEnd == 0)
        binEnd = hist3d.GetNbinsZ()+1;

    if (binStart >= binEnd)
        throw std::range_error("Invalid range specified. binStart must be less than binEnd");

    std::vector<TH2D> allHists;
    allHists.reserve(binEnd-binStart);
    for (size_t i = binStart; i < binEnd; i++) {
        std::string newName = name+to_string(i);
        allHists.emplace_back(projectTH2FromTH3(hist3d, newName.c_str(), i, i));
    }

    return allHists;
}

template <typename T>
std::vector<T> mirroredHists(T& nominal, std::vector<T>& hists) {
    std::vector<T> allHists;
    allHists.reserve(hists.size());
    for (auto& h : hists) {
        int i = 0;
        std::string name = h.GetName()+std::string("_mirror");
        T newh = *static_cast<T*>(h.Clone((name+std::to_string(i++)).c_str()));
        allHists.emplace_back(mirrorHist(nominal, h, newh));
    }
    return allHists;
}

template <typename T>
std::array<T, 2> envelopHists(std::vector<T>& hists) {
    std::string name = hists.at(0).GetName();
    T hup = *static_cast<T*>(hists.at(0).Clone(name+"Up"));
    T hdown = *static_cast<T*>(hup.Clone(name+"Down"));
    for (size_t i = 0; i <= hup.GetNcells(); i++) {
        float minc = hup.GetBinContent(i);
        float mine = hup.GetBinError(i);
        float maxc = hup.GetBinContent(i);
        float maxe = hup.GetBinError(i);

        for (auto& h : hists) {
            float cont = h.GetBinContent(i);
            if (cont > maxc) {
                maxc = cont;
                maxe = h.GetBinError(i);
            }
            else if (cont < minc) {
                minc = cont;
                mine = h.GetBinError(i);
            }
        }
        hup.SetBinContent(maxc);
        hup.SetBinError(maxe);
        hdown.SetBinContent(minc);
        hdown.SetBinError(mine);
    }

    return {hup, hdown};
}

#endif

