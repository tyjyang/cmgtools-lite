#ifndef THnDHelpers
#define THnDHelpers

#include "THn.h"
#include "TH1D.h"
#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <cstdlib> //as stdlib.h      
#include <cstdio>

//
// Would really rather create the TH(N+1)D here...
//THnD fillTHNplus1fromTHn(THnD& thn, int nbinLow=-1, int nbinHigh=-1, float scale=1.0) {
//
// THnD thnp1("luminosity_{k}", "luminosity_{k}", nNomiDim+1, array('i',nominalNbins+[2]), array('d',nominalBinMin+[0.5]), array('d',nominalBinMax+[2.5]))
// If binHigh is not specified the full range including under/overflow is filled
void fillTHNplus1fromTHn(THnD& thnp1, THnD& thn, int binLow=0, int binHigh=-1) {
    int ndim = thn.GetNdimensions();

    if (binHigh < 0)
        binHigh = thnp1.GetAxis(ndim)->GetNbins() + 1; // up to overflow

    if (binLow > binHigh)
        throw std::range_error("Invalid inputs! binLow must be less than binHigh");

    std::vector<int> binLookup(ndim+1, 0);
    for (Long64_t globalBin = 0; globalBin < thn.GetNbins(); globalBin++) {
        double binContent = thn.GetBinContent(globalBin, binLookup.data());
        double binError = thn.GetBinError(globalBin);
        for (int iNewDim = binLow; iNewDim <= binHigh; iNewDim++) {
            binLookup[ndim] = iNewDim;
            Long64_t globalBinTHnP1 = thnp1.GetBin(binLookup.data());
            thnp1.SetBinContent(globalBinTHnP1, binContent);
            thnp1.SetBinError(globalBinTHnP1, binError);
        }
    }

    //return thnp1;
}

template <typename T>
bool cropNegatives(T& h, size_t nbins, float cropValue, bool cropError) {
    bool hasCroppedBins = false;
    for (size_t globalBin = 0; globalBin <= nbins; globalBin++) {
        float binContent = h.GetBinContent(globalBin);
        if (binContent < 0.0) {
            hasCroppedBins = true;
            h.SetBinContent(globalBin, cropValue);
            if (cropError)
                h.SetBinError(globalBin, cropValue);
        }
    }
    return hasCroppedBins;
}

bool cropNegativeContent(THnD& h, bool silent = false, bool cropError = false, double cropValue = 0.0001) {

    Long64_t nbins = h.GetNbins();
    bool hasCroppedBins = cropNegatives<THnD>(h, nbins, cropValue, cropError);
    if (not silent and hasCroppedBins)
        std::cout << "Cropping negative bins for histogram " << h.GetName() << std::endl;

    return hasCroppedBins;
    
}


bool cropNegativeContent(TH1& h, bool silent = false, bool cropError = false, double cropValue = 0.0001) {

    Long64_t nbins = h.GetNcells();
    double integral = h.Integral();
    double integralNonNeg = 0.0;

    bool hasCroppedBins = cropNegatives<TH1>(h, nbins, cropValue, cropError);
    if (not silent and hasCroppedBins) {
        integralNonNeg = h.Integral();
        std::cout << h.GetName() << ": original integral = " << integral << " changed by " << integralNonNeg/integral << " (new/old)" << std::endl;
    }
    return hasCroppedBins;
    
}

#endif
