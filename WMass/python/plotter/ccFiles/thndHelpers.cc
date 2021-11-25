#ifndef THnDHelpers
#define THnDHelpers

#include "THn.h"
//
// Would really rather create the TH(N+1)D here...
//THnD fillTHNplus1fromTHn(THnD& thn, int nbinLow=-1, int nbinHigh=-1, float scale=1.0) {
//
// THnD thnp1("luminosity_{k}", "luminosity_{k}", nNomiDim+1, array('i',nominalNbins+[2]), array('d',nominalBinMin+[0.5]), array('d',nominalBinMax+[2.5]))
// If binHigh is not specified the full range including under/overflow is filled
void fillTHNplus1fromTHn(THnD& thnp1, THnD& thn, int binLow=0, int binHigh=-1) {
    int ndim = thn.GetNdimensions();

    if (binHigh < 0)
        binHigh = thnp1.GetAxis(ndim)->GetNbins();

    if (binLow > binHigh)
        throw std::range_error("Invalid inputs! binLow must be less than binHigh");

    std::vector<int> binLookup(ndim+1, 0);
    for (size_t globalBin = 0; globalBin < thn.GetNbins(); globalBin++) {
        float binContent = thn.GetBinContent(globalBin, binLookup.data());
        float binError = thn.GetBinError(globalBin);
        for (size_t iNewDim = binLow; iNewDim <= binHigh+1; iNewDim++) {
            binLookup[ndim] = iNewDim;
            int globalBinTHnP1 = thnp1.GetBin(binLookup.data());
            thnp1.SetBinContent(globalBinTHnP1, binContent);
            thnp1.SetBinError(globalBinTHnP1, binError);
        }
    }

    //return thnp1;
}

#endif
