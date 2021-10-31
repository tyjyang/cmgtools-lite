#!/usr/bin/env python3

import re
import os, os.path
import logging
import argparse
import shutil

## safe batch mode
import sys
args = sys.argv[:]
sys.argv = ['-b']
import ROOT
sys.argv = args
ROOT.gROOT.SetBatch(True)
ROOT.PyConfig.IgnoreCommandLineOptions = True

from copy import *

def wptBinsScales(i):
    # 5% quantiles (to be redone on the new MC)
    wptbins = [0.0, 1.971, 2.949, 3.838, 4.733, 5.674, 6.684, 7.781, 8.979, 10.303, 11.777, 13.435, 15.332, 17.525, 20.115, 23.245, 27.173, 32.414, 40.151, 53.858, 13000.0]
    if len(wptbins)<2*i:
        print('You are asking too much from the wpt binning for decorrelation of scales')
    ptlo = wptbins[2*(i-1)]
    pthi = wptbins[2*i]
    return [ptlo, pthi]

# TODO: Use this for other systematics (it needs to be generalized a bit)
def write3DHist(label, pt_expr, eta_expr, nsyst, etapt_binning, xylabels, weight_axis, regex, outfile=None, systBinStart=-0.5, indexStart=0, addWeight=None, replaceWeight=None, nBinsZaxis=None):

    if nBinsZaxis == None:
        nBinsZaxis = nsyst
    syst_binning = "%d,%.1f,%.1f" % (nBinsZaxis, systBinStart, nBinsZaxis+systBinStart)
    expr_string = f"indices({nsyst},{indexStart})\:scalarToRVec({pt_expr},{nsyst})\:scalarToRVec({eta_expr},{nsyst})"
    weight_items = []
    if addWeight:
        weight_items.append(f"AddWeight='{addWeight}'")
    if replaceWeight:
        weight_items.append(f"ReplaceWeight='{replaceWeight}'")
    weight_str = ", ".join(weight_items) if len(weight_items) else ""

    line = f"{label}_: {expr_string} : {etapt_binning},{syst_binning};" \
        f" {xylabels}, ZTitle='{weight_axis}', {weight_str}, ProcessRegexp='{regex}'\n"
    print(line)
    if outfile:
        outfile.write(line+'\n')

parser = argparse.ArgumentParser()
#parser.add_argument('-n', '--name', dest='baseHistName', default='muon_eta_pt', type=str, help='Base name for histograms')
parser.add_argument('-x', '--xAxisName', default='Muon #eta', type=str, help='x axis name')
parser.add_argument('-y', '--yAxisName', default='Muon p_{T} (GeV)', type=str, help='y axis name')
parser.add_argument('-b', '--bins', dest="etaptBins", default='48,-2.4,2.4,29,26,55', type=str, help='Bins for eta-pt, passed as to TH2 (only supports uniform binning for now)')
parser.add_argument('--ptVar', default='Muon_pt[goodMuonsCharge][0]', type=str, help='Expression for variable on pt axis')
parser.add_argument('--etaVar', default='Muon_eta[goodMuonsCharge][0]', type=str, help='Expression for variable on eta axis')
parser.add_argument('-a', '--analysis', choices=["wlike","wmass"], default="wmass", help='Analysis type (some settings are customized accordingly)')
parser.add_argument('-o', '--output', dest="outputFile", default='', type=str, help='Output file to store lines (they are also printed on stdout anyway)')
parser.add_argument('--pdf-weights', dest='pdfWeights', choices=["nnpdf30","nnpdf31"], default="nnpdf31", help='PDF set to use')
args = parser.parse_args()

###################################
# SOME BASIC CONFIGS #
######################
#baseHistName = args.baseHistName
axisNames = "XTitle='{x}', YTitle='{y}'".format(x=args.xAxisName,y=args.yAxisName)
etaptBins = args.etaptBins
isWlike = args.analysis == "wlike"
####################################

printToFile = False
outf = None
if args.outputFile != "":
    outf = open(args.outputFile,"w")
    printToFile = True
    
#nominal
#line = "{n}: {y}\:Muon_eta[goodMuonsCharge][0]: {b}; {axis} \n".format(n=baseHistName,y=args.ptVar,b=etaptBins,axis=axisNames)
line = f"nominal_: {args.ptVar}\:{args.etaVar}: {etaptBins}; {axisNames} \n"

print(line)
if printToFile: outf.write(line+'\n')

# pdf + alphaS in NNPDF3.1, alphaS has 0.002 variations around 0.118, the 1 sigma variation should be 0.0015
if args.pdfWeights == "nnpdf31":
    write3DHist(label = "pdfNNPDF31",
                pt_expr = args.ptVar,
                eta_expr = args.etaVar,
                nsyst = 103, # for PDFs, len(LHEPdfWeight) == 103 because it has nominal + 102 weights (100 pdf + 2 alphaS)
                xylabels = axisNames,
                weight_axis = "PDF/alpha_{S} index",
                etapt_binning = etaptBins,
                regex = "W.*|Z.*",
                outfile = outf,
                systBinStart = 0.5,
                indexStart = 0,
                addWeight = "LHEPdfWeight",
                nBinsZaxis = 102 # we want 102,0.5,102.5 (could have been 103,-0.5,102.5, but then histogram bin 1 would not be pdf1 but nominal, histgram bin 2 would not be pdf2 but pdf1 and so on)
    )
    line = f"alphaS0117NNPDF31_: {args.ptVar}\:{args.etaVar}: {etaptBins}; {axisNames}, AddWeight='LHEPdfWeightAltSet5[0]', ProcessRegexp='W.*'\n"
    line += f"alphaS0119NNPDF31_: {args.ptVar}\:{args.etaVar}: {etaptBins}; {axisNames}, AddWeight='LHEPdfWeightAltSet6[0]', ProcessRegexp='W.*'\n"
    print(line)
    if printToFile: outf.write(line+'\n')

elif args.pdfWeights == "nnpdf30":
    # pdf in NNPDF3.0, alphaS is available in another branch
    write3DHist(label = "pdfNNPDF30",
                pt_expr = args.ptVar,
                eta_expr = args.etaVar,
                nsyst = 101, # 1 nominal + 100 symmetric hessians
                xylabels = axisNames,
                weight_axis = "NNPDF3.0 PDF index",
                etapt_binning = etaptBins,
                regex = "W.*", # FIXME: only W for now
                outfile = outf,
                systBinStart = 0.5,
                indexStart = 0,
                addWeight = "(1.0/LHEPdfWeightAltSet13[0])*LHEPdfWeightAltSet13",  # when running analysis with different PDFs as nominal, the nominal weight LHEPdfWeightAltSet13[0] is alredy applied to the process, so we need to factorize it out here
                nBinsZaxis = 100 # we want 100,0.5,100.5 (could have been 101,-0.5,100.5, but then histogram bin 1 would not be pdf1 but nominal, histgram bin 2 would not be pdf2 but pdf1 and so on)
    )

    line = f"alphaS0117NNPDF30_: {args.ptVar}\:{args.etaVar}: {etaptBins}; {axisNames}, AddWeight='LHEPdfWeightAltSet15[0]', ProcessRegexp='W.*'\n"
    line += f"alphaS0119NNPDF30_: {args.ptVar}\:{args.etaVar}: {etaptBins}; {axisNames}, AddWeight='LHEPdfWeightAltSet16[0]', ProcessRegexp='W.*'\n"
    print(line)
    if printToFile: outf.write(line+'\n')

###  END of PDF if statement

### QCD scales should be used according to central PDF set. For now will keep using NNPDF3.1 in all cases

# qcd scales (not Vpt binned, that one comes just afterward)
write3DHist(label = "qcdScale",
            pt_expr = args.ptVar,
            eta_expr = args.etaVar,
            nsyst = 18, #  was 18 on older samples, because had both NNPDF3.0 and NNPDF3.1 together
            xylabels = axisNames,
            weight_axis = "QCD scale index",
            etapt_binning = etaptBins,
            regex = "W.*" if isWlike else "Z.*",
            outfile = outf,
            systBinStart = -0.5,
            indexStart = 0,
            addWeight = "LHEScaleWeight" 
)

# qcd scales (Vpt binned)
NVTPBINS = 10
process_regexpr =  "Z.*" if isWlike else "W.*" # opposite with respect to unbinned QCD scales
for ipt in range(1,1+NVTPBINS):
    ptcut = wptBinsScales(ipt)
    write3DHist(label = "qcdScaleVptBin%d" % ipt,
                pt_expr = args.ptVar,
                eta_expr = args.etaVar,
                nsyst = 9, #  was 18 on older samples, because had both NNPDF3.0 and NNPDF3.1 together   
                xylabels = axisNames,
                weight_axis = "QCD scale index",
                etapt_binning = etaptBins,
                regex = "Z.*" if isWlike else "W.*", # opposite with respect to unbinned QCD scales 
                outfile = outf,
                systBinStart = -0.5,
                indexStart = 0,
                addWeight = f"qcdScaleWeight_VptBinned(LHEScaleWeight\,ptVgen\,{ptcut[0]}\,{ptcut[1]})" 
    )


## end of QCD scales
    
# eff. stat. nuisances, one nuisance per TnP bin, treated as uncorrelated
# function to use is _get_fullMuonSFvariation, which replace _get_fullMuonSF in the nominal weight, using ReplaceWeight
# NOTE: from September 2021 we have changed pt binning, now it is 15 pt bins from 24 to 65, the bins are: 24, 26, 28, 30, 32, 34, 36, 38, 40, 42, 44, 47, 50, 55, 60, 65
# for the analysis we may cut between 26 and 55, maybe even in a narrower window, so we would need a way to remove the unnecessary bins
# for now we can neglect the upper pt bins, for the lowest ones we need some tricks to be implemented, so for now we have 48(eta) * 13(pt) = 624 bins (up to 55)
write3DHist(label = "effStatTnP",
            pt_expr = args.ptVar,
            eta_expr = args.etaVar,
            nsyst = 624, # remember to edit weight below if changing this 
            xylabels = axisNames,
            weight_axis = "Eff. stat. nuisance index",
            etapt_binning = etaptBins,
            regex = "W.*|Z.*|Top|Diboson", # no fakes here yet
            outfile = outf,
            systBinStart = 0.5,
            indexStart = 1,
            replaceWeight = f"_get_fullMuonSF(->_get_fullMuonSFvariation(624\," 
)

# eff. syst. nuisance, for now one just nuisance correlating all bins (may decorrelate across eta but for now we keep it simple)
# so we only need to make a TH2
# function to use is _get_fullMuonSF_dataAltSig, which replace _get_fullMuonSF in the nominal weight, using ReplaceWeight
# this is basically the SF obtained with alternative data efficiency from fit with analytic function (nominal uses MC template for signal instead)
# this is just one variation, we would make it symmetric later
line = f"effSystTnP_: {args.ptVar}\:{args.etaVar}: {etaptBins}; {axisNames}, ReplaceWeight='_get_fullMuonSF(->_get_fullMuonSF_dataAltSig(', ProcessRegexp='W.*|Z.*|Top|Diboson'\n"
print(line)
if printToFile: outf.write(line+'\n')


# muon prefiring stat uncertainty, uncorrelated for each eta bin. Here the function is _get_newMuonPrefiringSFvariationStat, which replace _get_newMuonPrefiringSF in the nominal weight, using ReplaceWeight
write3DHist(label = "muonL1PrefireStat",
            pt_expr = args.ptVar,
            eta_expr = args.etaVar,
            nsyst = 11, # remember to edit weight below if changing this 
            xylabels = axisNames,
            weight_axis = "Muon L1 prefiring nuisance index",
            etapt_binning = etaptBins,
            regex = "W.*|Z.*|Top|Diboson", # no fakes here yet
            outfile = outf,
            systBinStart = 0.5,
            indexStart = 1,
            replaceWeight = f"_get_newMuonPrefiringSF(->_get_newMuonPrefiringSFvariationStat(11\," 
            #replaceWeight = f"_get_newMuonPrefiringSF(Muon_eta\,Muon_pt\,Muon_phi\,Muon_looseId\,eraVFP)->_get_newMuonPrefiringSFvariationStat(11\,Muon_eta\,Muon_pt\,Muon_phi\,Muon_looseId\,eraVFP)" 
)

# muon prefiring syst uncertainty, correlated across eta bins. Here the function is _get_newMuonPrefiringSFvariationSyst, which replace _get_newMuonPrefiringSF in the nominal weight, using ReplaceWeight
write3DHist(label = "muonL1PrefireSyst",
            pt_expr = args.ptVar,
            eta_expr = args.etaVar,
            nsyst = 3, # remember to edit weight below if changing this 
            xylabels = axisNames,
            weight_axis = "Muon L1 prefiring nuisance index",
            etapt_binning = etaptBins,
            regex = "W.*|Z.*|Top|Diboson", # no fakes here yet
            outfile = outf,
            systBinStart = 0.5,
            indexStart = 1,
            #replaceWeight = f"_get_newMuonPrefiringSF(Muon_eta\,Muon_pt\,Muon_phi\,Muon_looseId\,eraVFP)->_get_newMuonPrefiringSFvariationSyst(3\,Muon_eta\,Muon_pt\,Muon_phi\,Muon_looseId\,eraVFP)" 
            replaceWeight = f"_get_newMuonPrefiringSF(->_get_newMuonPrefiringSFvariationSyst(" 
)


# mass weights from LHEReweightingWeightCorrectMass
write3DHist(label = "massWeight",
            pt_expr = args.ptVar,
            eta_expr = args.etaVar,
            nsyst=23 if isWlike else 21,
            xylabels = axisNames,
            weight_axis = "Mass weight index",
            etapt_binning = etaptBins,
            regex = "Z.*" if isWlike else "W.*",
            outfile = outf,
            systBinStart = -0.5,
            indexStart = 0,
            addWeight = "MEParamWeight"# is LHEReweightingWeightCorrectMass on older samples
)


# for chan in ["Wm", "Wp", "Z"]:
#     process_regexpr = "Z.*" 
#     if "W" in chan:
#         process_regexpr = ".*".join(["W", "minus" if "m" in chan else "plus", ])
#     START = -0.5
#     NWEIGHTS = 45
#     write3DHist(label ="scetlibWeights", # no trailing '_', it is added inside directly (and mcPlots will add another one) 
#                 pt_expr =args.ptVar,
#                 eta_expr = args.etaVar,
#                 nsyst = NWEIGHTS,
#                 xylabels = axisNames,
#                 weight_axis = "SCETlib variation index",
#                 etapt_binning = etaptBins,
#                 regex = process_regexpr,
#                 outfile = outf,
#                 systBinStart = START,
#                 indexStart = 0,
#                 addWeight = f"scetlibWeights{chan}"
#     )

print('-'*30)
print("SUMMARY")
print('-'*30)
print("Analysis:       %s" % args.analysis)
#print("Base hist name: %s" % baseHistName)
print("Binning eta-pt: %s" % etaptBins)
print('-'*30)
if printToFile:
    outf.close()
    print("Output saved in file %s" % args.outputFile)
    print()
