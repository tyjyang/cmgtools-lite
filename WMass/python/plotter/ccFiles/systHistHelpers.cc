#ifndef systHistHelpers
#define systHistHelpers

#include "TH2.h"

TH2D mirrorShape(TH2D& nominal, TH2D& alternate, TH2D& mirrored) {
    for (size_t bx = 0; bx <= nominal.GetNbinsX()+1; bx++) {
        for (size_t by = 0; by <= nominal.GetNbinsY()+1; by++) {
            float y0 = nominal.GetBinContent(bx, by);
            float yA = alternate.GetBinContent(bx, by);
            float yM = (yA > 0 && y0 > 0) ? y0*y0/yA : 0;
            mirrored.SetBinContent(bx, by, yM);
            mirrored.SetBinError(bx, by, alternate.GetBinError(bx, by));
        }
    }
    return mirrored;
}

#endif

