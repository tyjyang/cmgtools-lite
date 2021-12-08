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

def writeNDHist(label, varExpr, nsyst, binning, axisLabels, weightAxisLabel, procRegexp,
                outfile=None, systBinStart=-0.5, indexStart=0,
                addWeight=None, replaceWeight=None, nBinsSystAxis=None, replaceCutByName=None, isPtScaleTest=False):

    if nBinsSystAxis == None:
        nBinsSystAxis = nsyst
    syst_binning = "%d,%.1f,%.1f" % (nBinsSystAxis, systBinStart, nBinsSystAxis+systBinStart)
    syst_expr = f"indices({nsyst},{indexStart})"
    weight_items = []
    if addWeight:
        weight_items.append(f"AddWeight='{addWeight}'")
    if replaceWeight:
        weight_items.append(f"ReplaceWeight='{replaceWeight}'")
    if replaceCutByName:
        weight_items.append(f"ReplaceCutByName='{replaceCutByName}'")
    weight_str = ", ".join(weight_items) if len(weight_items) else ""

    fullAxisLabels = f"{axisLabels[:-1]}\,{weightAxisLabel}'"  # remove ' from axisLabels, append weightAxisLabel, and add ' back
    line = f"{label}_: {syst_expr}\:{varExpr} : {binning},{syst_binning}; {fullAxisLabels}, {weight_str}, ProcessRegexp='{procRegexp}'\n"
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
        "alphaRange" : "002",
    },
    "nnpdf30" : {
        "name" : "NNPDF30",
        "entries" : 101,
        "truncate" : False,
        "onlyW" : True,
        "weight" : "LHEPdfWeightAltSet13",
        "alphas" : ["LHEPdfWeightAltSet15[0]", "LHEPdfWeightAltSet16[0]"],
        "alphaRange" : "002",
    },
    "ct18" : {
        "name" : "CT18",
        "entries" : 59,
        "onlyW" : True,
        "truncate" : True,
        "weight" : "LHEPdfWeightAltSet18",
        "alphas" : ["LHEPdfWeightAltSet18[59]", "LHEPdfWeightAltSet18[60]"],
        "alphaRange" : "002",
    },
    "mmht" : {
        "name" : "MMHT",
        "entries" : 51,
        "onlyW" : True,
        "truncate" : False,
        "weight" : "LHEPdfWeightAltSet19",
        "alphas" : ["LHEPdfWeightAltSet20[1]", "LHEPdfWeightAltSet20[2]"],
        "alphaRange" : "002",
    },
}

parser = argparse.ArgumentParser()
parser.add_argument("-d", "--dimension", nargs="+", default=['Muon_eta[goodMuons][0];Muon #eta;48;-2.4;2.4',
                                                             'Muon_pt[goodMuons][0];Muon p_{T} (GeV);29;26;55',
                                                             'Muon_charge[goodMuons][0];Muon charge;2;-2;2',
                                                             'regionIsoMt(Muon_pfRelIso04_all[goodMuons][0]<0.15,transverseMass<40);isolation-mT bin;4;-0.5;3.5'],
                    help="Add dimension to fill new axis. Can pass multiple arguments as ';'-separated list of 5 elements (variable, title, number of bins, minimum, maximum)")
parser.add_argument('-a', '--analysis', choices=["wlike","wmass"], default="wmass",
                    help='Analysis type (some settings are customized accordingly)')
parser.add_argument('-o', '--output', dest="outputFile", default='', type=str,
                    help='Output file to store lines (they are also printed on stdout anyway)')
parser.add_argument('--pdf-weights', dest='pdfWeights', choices=pdfMap.keys(), default="nnpdf31", help='PDF set to use')
#parser.add_argument('--ptVarScaleTest', default='customPtTest', type=str,
#                    help='Expression for variable on pt axis, specific for tests about pt scale')
args = parser.parse_args()

###################################
# SOME BASIC CONFIGS #
######################
variables = []
titles    = []
bins   = []
for dim in args.dimension:
    tokens = dim.split(";")
    if len(tokens) != 5:
        print(f"Warning, bad formatted line in option --dimension: {dim}. It requires 5 fields separated by ';'.")
        quit()
    # some sanity check on binning
    if float(tokens[3]) > float(tokens[4]):
        print(f"Warning, bad formatted line in option --dimension: {dim}. Axis minimum is larger than maximum.")
        quit()
    variables.append(f"{tokens[0]}")
    titles.append(f"{tokens[1]}")
    bins.append(",".join(tokens[2:]))

expression = "\:".join(variables[::-1])
axisNames  = "NTitle='{a}'".format(a="\,".join(titles))
binning    = ",".join(bins)
print("-"*30)
print("Default arguments for nominal histogram")
print(f"Expression: {expression}")
print(f"Axis names: {axisNames}")
print(f"Binning   : {binning}")
print("-"*30)
#quit()

isWlike = args.analysis == "wlike"
####################################

printToFile = False
outf = None
if args.outputFile != "":
    outf = open(args.outputFile,"w")
    printToFile = True
    
#nominal
line = f"nominal_: {expression}: {binning}; {axisNames} \n"
print(line)
if printToFile: outf.write(line+'\n')

pdfInfo = pdfMap[args.pdfWeights]
# Needed if you don't want to use all the PDF weights in the vector (mostly needed for CT18, which didn't get parsed correctly)
weight = pdfInfo["weight"]
entries=pdfInfo["entries"]
weightexpr = weight if "Alt" not in weight else f"1.0/{weight}[0]*" + (str(weight) if not pdfInfo["truncate"] else
    f" ROOT::VecOps::RVec<float>(std::begin({weight})\, std::begin({weight})+{entries})")
print(weightexpr)

regex = "W.*|Z.*" if not pdfInfo["onlyW"] else "W.*"
writeNDHist(label = "pdf"+pdfInfo["name"],
            varExpr = expression,
            nsyst = pdfInfo["entries"], 
            axisLabels = axisNames,
            weightAxisLabel = "PDF/alpha_{S} index",
            binning = binning,
            procRegexp = regex,
            outfile = outf,
            systBinStart = 0.5,
            indexStart = 0,
            addWeight = weightexpr,
            nBinsSystAxis = pdfInfo["entries"]-1
)

writeNDHist(label="alphaS{arange}{name}".format(arange=pdfInfo["alphaRange"], name=pdfInfo["name"]),
    varExpr = expression,
    nsyst = 2,
    axisLabels = axisNames,
    weightAxisLabel = "alpha_s index",
    binning = binning,
    procRegexp = "W.*",
    outfile = outf,
    systBinStart = -0.5,
    indexStart = 1,
    addWeight = "ROOT::VecOps::RVec<float>{{{0}\,{1}}}".format(*pdfInfo["alphas"])
)

### QCD scales should be used according to central PDF set. For now will keep using NNPDF3.1 in all cases

# qcd scales (not Vpt binned, that one comes just afterward)
writeNDHist(label = "qcdScale",
            varExpr = expression,
            nsyst = 18, #  was 18 on older samples, because had both NNPDF3.0 and NNPDF3.1 together
            axisLabels = axisNames,
            weightAxisLabel = "QCD scale index",
            binning = binning,
            procRegexp = "W.*" if isWlike else "Z.*",
            outfile = outf,
            systBinStart = -0.5,
            indexStart = 0,
            addWeight = "LHEScaleWeight" 
)

# qcd scales (Vpt binned)
NVTPBINS = 10
for ipt in range(1,1+NVTPBINS):
    ptcut = wptBinsScales(ipt)
    writeNDHist(label = "qcdScaleVptBin%d" % ipt,
                varExpr = expression,
                nsyst = 9, #  was 18 on older samples, because had both NNPDF3.0 and NNPDF3.1 together   
                axisLabels = axisNames,
                weightAxisLabel = "QCD scale index (Vpt binned)",
                binning = binning,
                procRegexp = "Z.*" if isWlike else "W.*", # opposite with respect to unbinned QCD scales 
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
writeNDHist(label = "effStatTnP",
            varExpr = expression,
            nsyst = 624, # remember to edit weight below if changing this 
            axisLabels = axisNames,
            weightAxisLabel = "Eff. stat. nuisance index",
            binning = binning,
            procRegexp = "W.*|Z.*|Top|Diboson", # no fakes here yet
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
line = f"effSystTnP_: {expression}: {binning}; {axisNames}, ReplaceWeight='_get_fullMuonSF(->_get_fullMuonSF_dataAltSig(', ProcessRegexp='W.*|Z.*|Top|Diboson'\n"
print(line)
if printToFile: outf.write(line+'\n')


# muon prefiring stat uncertainty, uncorrelated for each eta bin. Here the function is _get_newMuonPrefiringSFvariationStat, which replace _get_newMuonPrefiringSF in the nominal weight, using ReplaceWeight
writeNDHist(label = "muonL1PrefireStat",
            varExpr = expression,
            nsyst = 11, # remember to edit weight below if changing this 
            axisLabels = axisNames,
            weightAxisLabel = "Muon L1 prefiring stat nuisance index",
            binning = binning,
            procRegexp = "W.*|Z.*|Top|Diboson", # no fakes here yet
            outfile = outf,
            systBinStart = 0.5,
            indexStart = 1,
            replaceWeight = f"_get_newMuonPrefiringSF(->_get_newMuonPrefiringSFvariationStat(11\," 
            #replaceWeight = f"_get_newMuonPrefiringSF(Muon_eta\,Muon_pt\,Muon_phi\,Muon_looseId\,eraVFP)->_get_newMuonPrefiringSFvariationStat(11\,Muon_eta\,Muon_pt\,Muon_phi\,Muon_looseId\,eraVFP)" 
)

# muon prefiring syst uncertainty, correlated across eta bins. Here the function is _get_newMuonPrefiringSFvariationSyst, which replace _get_newMuonPrefiringSF in the nominal weight, using ReplaceWeight
writeNDHist(label = "muonL1PrefireSyst",
            varExpr = expression,
            nsyst = 3, # remember to edit weight below if changing this 
            axisLabels = axisNames,
            weightAxisLabel = "Muon L1 prefiring syst nuisance index",
            binning = binning,
            procRegexp = "W.*|Z.*|Top|Diboson", # no fakes here yet
            outfile = outf,
            systBinStart = 0.5,
            indexStart = 1,
            #replaceWeight = f"_get_newMuonPrefiringSF(Muon_eta\,Muon_pt\,Muon_phi\,Muon_looseId\,eraVFP)->_get_newMuonPrefiringSFvariationSyst(3\,Muon_eta\,Muon_pt\,Muon_phi\,Muon_looseId\,eraVFP)" 
            replaceWeight = f"_get_newMuonPrefiringSF(->_get_newMuonPrefiringSFvariationSyst(" 
)


# mass weights from LHEReweightingWeightCorrectMass
writeNDHist(label = "massWeight",
            varExpr = expression,
            nsyst=23 if isWlike else 21,
            axisLabels = axisNames,
            weightAxisLabel = "Mass weight index",
            binning = binning,
            procRegexp = "Z.*" if isWlike else "W.*",
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
#     writeNDHist(label ="scetlibWeights", # no trailing '_', it is added inside directly (and mcPlots will add another one) 
#                 varExpr = expression,
#                 nsyst = NWEIGHTS,
#                 axisLabels = axisNames,
#                 weightAxisLabel = "SCETlib variation index",
#                 binning = binning,
#                 procRegexp = process_regexpr,
#                 outfile = outf,
#                 systBinStart = START,
#                 indexStart = 0,
#                 addWeight = f"scetlibWeights{chan}"
#     )

# muon momentum scale test
# writeNDHist(label = "muonPtScaleTest",
#             varExpr = expression,  # FIXME
#             nsyst = 96, # 48 etaBins*2 as we do also up/down together, remember to edit weight below if changing this 
#             axisLabels = axisNames,
#             weightAxisLabel = "pT scale nuisance index (Up+Down)",
#             binning = binning,
#             procRegexp = "W.*|Z.*|Top|Diboson", # no fakes here yet
#             outfile = outf,
#             systBinStart = 0.5,
#             indexStart = 1,
#             replaceCutByName = "accept->valueInsideRange(customPtTest\,26\,55)",
#             isPtScaleTest = True
# )


print('-'*30)
print("SUMMARY")
print('-'*30)
print(f"Analysis  : {args.analysis}")
#print("Base hist name: %s" % baseHistName)
print(f"Expression: {expression}")
print(f"Axis names: {axisNames}")
print(f"Binning   : {binning}")
print("-"*30)
if printToFile:
    outf.close()
    print(f"Output saved in file {args.outputFile}")
    print()
