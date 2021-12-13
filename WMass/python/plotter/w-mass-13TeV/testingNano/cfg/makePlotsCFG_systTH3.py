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
def write3DHist(label, pt_expr, eta_expr, nsyst, etapt_binning, xylabels, weight_axis, regex, outfile=None, systBinStart=-0.5, indexStart=0, addWeight=None, replaceWeight=None, nBinsZaxis=None, replaceCutByName=None, isPtScaleTest=False):

    if nBinsZaxis == None:
        nBinsZaxis = nsyst
    syst_binning = "%d,%.1f,%.1f" % (nBinsZaxis, systBinStart, nBinsZaxis+systBinStart)
    if isPtScaleTest:
        expr_string = f"indices({nsyst},{indexStart})\:{pt_expr}\:scalarToRVec({eta_expr},{nsyst})"
    else:
        expr_string = f"indices({nsyst},{indexStart})\:scalarToRVec({pt_expr},{nsyst})\:scalarToRVec({eta_expr},{nsyst})"
    weight_items = []
    if addWeight:
        weight_items.append(f"AddWeight='{addWeight}'")
    if replaceWeight:
        weight_items.append(f"ReplaceWeight='{replaceWeight}'")
    if replaceCutByName:
        weight_items.append(f"ReplaceCutByName='{replaceCutByName}'")
    weight_str = ", ".join(weight_items) if len(weight_items) else ""

    line = f"{label}_: {expr_string} : {etapt_binning},{syst_binning};" \
        f" {xylabels}, ZTitle='{weight_axis}', {weight_str}, ProcessRegexp='{regex}'\n"
    print(line)
    if outfile:
        outfile.write(line+'\n')

pdfMap = {
    "nnpdf31" : {
        "name" : "NNPDF31",
        "entries" : 103,
        "truncate" : False,
        "onlyW" : False,
        "weight" : "LHEPdfWeight",
        "alphas" : ["LHEPdfWeightAltSet5[0]", "LHEPdfWeightAltSet6[0]"],
    },
    "nnpdf30" : {
        "name" : "NNPDF30",
        "entries" : 101,
        "truncate" : False,
        "onlyW" : True,
        "weight" : "LHEPdfWeightAltSet13",
        "alphas" : ["LHEPdfWeightAltSet15[0]", "LHEPdfWeightAltSet16[0]"],
    },
    "ct18" : {
        "name" : "CT18",
        "entries" : 59,
        "onlyW" : True,
        "truncate" : True,
        "weight" : "LHEPdfWeightAltSet18",
        "alphas" : ["LHEPdfWeightAltSet18[59]", "LHEPdfWeightAltSet18[60]"],
    },
    "mmht" : {
        "name" : "MMHT",
        "entries" : 51,
        "onlyW" : True,
        "truncate" : False,
        "weight" : "LHEPdfWeightAltSet19",
        "alphas" : ["LHEPdfWeightAltSet20[1]", "LHEPdfWeightAltSet20[2]"],
    },
}



parser = argparse.ArgumentParser()
#parser.add_argument('-n', '--name', dest='baseHistName', default='muon_eta_pt', type=str, help='Base name for histograms')
parser.add_argument('-x', '--xAxisName', default='Muon #eta', type=str, help='x axis name')
parser.add_argument('-y', '--yAxisName', default='Muon p_{T} (GeV)', type=str, help='y axis name')
parser.add_argument('-b', '--bins', dest="etaptBins", default='48,-2.4,2.4,29,26,55', type=str, help='Bins for eta-pt, passed as to TH2 (only supports uniform binning for now)')
parser.add_argument('--ptVar', default='Muon_pt[goodMuonsCharge][0]', type=str, help='Expression for variable on pt axis')
parser.add_argument('--etaVar', default='Muon_eta[goodMuonsCharge][0]', type=str, help='Expression for variable on eta axis')
parser.add_argument('-a', '--analysis', choices=["wlike","wmass"], default="wmass", help='Analysis type (some settings are customized accordingly)')
parser.add_argument('-o', '--output', dest="outputFile", default='', type=str, help='Output file to store lines (they are also printed on stdout anyway)')
parser.add_argument('--pdf-weights', dest='pdfWeights', choices=pdfMap.keys(), default="nnpdf31", help='PDF set to use')
parser.add_argument('--ptVarScaleTest', default='customPtTest', type=str, help='Expression for variable on pt axis, specific for tests about pt scale (it doesn\'t get vectorized)')
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

pdfInfo = pdfMap[args.pdfWeights]
# Needed if you don't want to use all the PDF weights in the vector (mostly needed for CT18, which didn't get parsed correctly)
weight = pdfInfo["weight"]
entries=pdfInfo["entries"]
weightexpr = weight if "Alt" not in weight else f"1.0/{weight}[0]*" + (str(weight) if not pdfInfo["truncate"] else
    f" ROOT::VecOps::RVec<float>(std::begin({weight})\, std::begin({weight})+{entries})")
print(weightexpr)
#exit(1)

regex = "W.*|Z.*" if not pdfInfo["onlyW"] else "W.*"
write3DHist(label = "pdf"+pdfInfo["name"],
            pt_expr = args.ptVar,
            eta_expr = args.etaVar,
            nsyst = pdfInfo["entries"], 
            xylabels = axisNames,
            weight_axis = "PDF index",
            etapt_binning = etaptBins,
            regex = regex,
            outfile = outf,
            systBinStart = 0.5,
            indexStart = 0,
            addWeight = weightexpr,
            nBinsZaxis = pdfInfo["entries"]-1
)
line = f"alphaS0117%s_: {args.ptVar}\:{args.etaVar}: {etaptBins}; {axisNames}, AddWeight='%s', ProcessRegexp='W.*'\n" % (pdfInfo["name"], pdfInfo["alphas"][0])
line += f"alphaS0119%s_: {args.ptVar}\:{args.etaVar}: {etaptBins}; {axisNames}, AddWeight='%s', ProcessRegexp='W.*'\n" % (pdfInfo["name"], pdfInfo["alphas"][1])

print(line)
if printToFile: outf.write(line+'\n')

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

# muon momentum scale test
write3DHist(label = "muonPtScaleTest",
            pt_expr = args.ptVarScaleTest,
            eta_expr = args.etaVar,
            nsyst = 96, # 48 etaBins*2 as we do also up/down together, remember to edit weight below if changing this 
            xylabels = axisNames,
            weight_axis = "pT scale nuisance index (Up+Down)",
            etapt_binning = etaptBins,
            regex = "W.*|Z.*|Top|Diboson", # no fakes here yet
            outfile = outf,
            systBinStart = 0.5,
            indexStart = 1,
            replaceCutByName = "accept->valueInsideRange(customPtTest\,26\,55)",
            isPtScaleTest = True
)


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
