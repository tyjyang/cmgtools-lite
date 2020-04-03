#!/usr/bin/env python
from shutil import copyfile
import re, sys, os, os.path, subprocess

from optparse import OptionParser
parser = OptionParser(usage="%prog testname ")
parser.add_option("--mu", dest="useMuon", default=False, action='store_true', help="Do fake rate for muons");
parser.add_option("--qcdmc", dest="addQCDMC", default=False, action='store_true', help="Add QCD MC in plots (but do not subtract from data)");
parser.add_option("--full2016data", dest="useFullData2016", default=False, action='store_true', help="Use all 2016 data (B to H, 35.9/fb). By default, only B to F are used (19.3/fb. Luminosity is automatically set depending on this choice");
parser.add_option("--charge", dest="charge", default="", type='string', help="Select charge: p for positive, n for negative");
parser.add_option("--lumi", dest="lumi", default=35.9 , type='float', help="Integrated luminosity");
parser.add_option("--test", dest="test", default="", type='string', help="pass the name of a folder (mandatory) to store test FR plots. It is created in plots/fake-rate/test/");
parser.add_option("--fqcd-ranges", dest="fqcd_ranges", default="0,40,50,120", type='string', help="Pass a list of 4 comma separated numbers that represents the ranges for the two mT regions to compute the fake rate");
parser.add_option("--pt", dest="ptvar", default="pt_granular", type='string', help="Select pT definition: pt_granular (default) or pt_coarse, or pt_finer (for muons)");
parser.add_option("--useSkim", dest="useSkim", default=False, action='store_true', help="Use skimmed sample for fake rates");
parser.add_option("--usePickle", dest="usePickle", default=False, action='store_true', help="Read SumWeigth from pickle file (in case the hisotgram with this number is not present in the ntuples)");
parser.add_option("--no-scaleFactors", dest="noScaleFactors", default=False, action='store_true', help="Don't use lepton scale factors (only PU weight)");
parser.add_option("--useSignedEta", dest="useSignedEta", default=False, action='store_true', help="Make fake rate for eta bins distinguishing eta sign");
parser.add_option("--makeTH3-eta-pt-passID", dest="makeTH3_eta_pt_passID", default=False, action='store_true', help="This option is special. It will make the following create only the TH3D histograms with |eta|, pt and passID. This will allow to compute the FR in any binning of pt and eta (using another macro). It overrides some other options");
parser.add_option("--addOpts", dest="addOpts", default="", type='string', help="Options to pass some other options from outside to build the command");
parser.add_option("--reweightZpt", dest="reweightZpt", default=False, action='store_true', help="Use W and Z with reweighted pT");
(options, args) = parser.parse_args()

useMuon = options.useMuon
addQCDMC = options.addQCDMC  # trying to add QCD MC to graphs to be compared
charge = str(options.charge)
testDir = str(options.test)
fqcd_ranges = str(options.fqcd_ranges)
ptvar = str(options.ptvar)
useFullData2016 = options.useFullData2016
useSkim = options.useSkim
usePickle = options.usePickle
useSignedEta = options.useSignedEta
addOpts = options.addOpts
luminosity = options.lumi
reweightZpt = options.reweightZpt

print "# useFullData2016 = " + str(useFullData2016)

if useSignedEta:
    fitvar = "etal1mu" if useMuon else "etal1"
    etaRange = [ '-2.4', '2.4'] if useMuon else [ '-2.5', '2.5']
else:
    fitvar = "absetal1mu" if useMuon else "absetal1"  
    etaRange = [ '0.0', '2.4'] if useMuon else [ '0.0', '2.5']  
    


if not useMuon and ptvar not in ["pt_coarse", "pt_granular"]:
    print "warning: unknown pt definition %s, use pt_coarse, pt_granular" % ptvar
    quit()

if useMuon:
    addQCDMC = True

plotterPath = str(os.environ.get('CMSSW_BASE'))
plotterPath = plotterPath + "/src/CMGTools/WMass/python/plotter/"

chargeSelection = ""
if charge != "":
    if charge == "p":
        chargeSelection = "-A onelep positive 'LepGood1_charge > 0'"
    elif charge == "n":
        chargeSelection = "-A onelep negative 'LepGood1_charge < 0'"
    else:
        print "%s is not a valid input for charge setting: use p or n" % charge
        quit()

T="/eos/cms/store/group/dpg_ecal/comm_ecal/localreco/TREES_1LEP_80X_V3/" 
if useSkim:
    #T="/eos/cms/store/group/dpg_ecal/comm_ecal/localreco/TREES_1LEP_80X_V3_FRELSKIM_V5/"
    T="/eos/cms/store/cmst3/group/wmass/mciprian/TREES_1LEP_80X_V3_FRELSKIM_V9/"
if useMuon:
    #T="/afs/cern.ch/work/e/emanuele/TREES/TREES_wlike_mu_V1/"
    T="/eos/cms/store/cmst3/group/wmass/mciprian/TREE_4_WMASS_skimMuonFR_21Jan2020/"
objName='tree' # name of TTree object in Root file, passed to option --obj in tree2yield.py
print "# used trees from: ",T

ptcorr = "ptMuFull(LepGood1_calPt,LepGood1_eta)" if useMuon else "ptElFull(LepGood1_calPt,LepGood1_eta)"
if useFullData2016:
    print "# Using full 2016 dataset"
    MCweightOption = ' -W "puw2016_nTrueInt_36fb(nTrueInt)*_get_electronSF_TriggerAndID(LepGood1_pdgId,LepGood1_calPt,LepGood1_eta)*LepGood1_SF2*eleSF_L1Eff(LepGood1_pt,LepGood1_eta)" ' # with L1 prefire 

if useMuon:
    #MCweightOption = ' -W "puw2016_nTrueInt_36fb(nTrueInt)*_get_muonSF_recoToSelection(LepGood1_pdgId,LepGood1_calPt,LepGood1_eta)*_get_muonSF_selectionToTrigger(LepGood1_pdgId,LepGood1_calPt,LepGood1_eta,LepGood1_charge)*prefireJetsWeight(LepGood1_eta)" ' 
    # test no sf, use them only for numerator
    #MCweightOption = ' -W "puw2016_nTrueInt_36fb(nTrueInt)*prefireJetsWeight(LepGood1_eta)*_get_muonSF_TriggerAndIDiso(LepGood1_pdgId,LepGood1_calPt,LepGood1_eta,LepGood1_charge,(LepGood1_relIso04 < 0.15))" ' 
    MCweightOption = ' -W "puw2016_nTrueInt_36fb(nTrueInt)*prefireJetsWeight(LepGood1_eta)*_get_muonSF_TriggerAndIDiso(LepGood1_pdgId,LepGood1_calPt,LepGood1_eta,LepGood1_charge)" ' 
    # _get_muonSF_recoToSelection(LepGood1_pdgId,LepGood1_calPt,LepGood1_eta)*_get_muonSF_selectionToTrigger(LepGood1_pdgId,LepGood1_calPt,LepGood1_eta,LepGood1_charge)

if options.noScaleFactors:
    print "# Warning: not using lepton scale factors: only PU weight"
    MCweightOption = ' -W "puw2016_nTrueInt_36fb(nTrueInt)*prefireJetsWeight(LepGood1_eta)" '

J=4

BASECONFIG = plotterPath + "w-mass-13TeV/wmass_e"
MCA = BASECONFIG+'/mca-80X_V5_FRskim.txt'
CUTFILE =BASECONFIG+'/qcd1l_SRtrees.txt'
XVAR=ptvar
FITVAR=fitvar
NUM = "fakeRateNumerator_el"

if useMuon:
    BASECONFIG=plotterPath + "w-mass-13TeV/wmass_mu"
    MCA=BASECONFIG+'/mca_forMuonFR.txt'
    CUTFILE=BASECONFIG+'/cuts_forMuonFR.txt'
    XVAR=ptvar
    FITVAR=fitvar
    NUM = "fakeRateNumerator_mu"

OPTIONS = MCA + " " + CUTFILE + " -f -P " + T + " --obj " + objName + " --s2v -j " + str(J) + " -l " + str(luminosity) + " " + str(addOpts)
 
# no friends for the moment
OPTIONS += ' -F Friends '+T+'/friends/tree_Friend_{cname}.root '
# OPTIONS += ' -F Friends '+T+'/friends/tree_FRFriend_{cname}.root '
# OPTIONS += ' --FMC Friends '+T+'/friends/tree_TrgFriend_{cname}.root '  # only for MC, they have trigger scale factors
OPTIONS += ' --fqcd-ranges %s' % fqcd_ranges.replace(","," ")
#OPTIONS += datasetOption

if usePickle:
    OPTIONS += ' --usePickle '

# event weight (NB: not needed for data, and efficiency sf not good for MC since here we have fake electrons)
# use PU reweighting for BF or BH
OPTIONS += MCweightOption

PBASE = plotterPath + "plots/fake-rate/el/"
if useMuon:
    PBASE = plotterPath + "plots/fake-rate/mu/"
if testDir != "":
    PBASE = PBASE.replace('plots/fake-rate/','plots/fake-rate/test/'+str(testDir)+'/')
    if useMuon:
        PBASE = PBASE.replace("plots/fake-rate/test/","plots/fake-rate/test_mu/")

if charge == "p":
    PBASE = PBASE + "pos/"
elif charge == "n":
    PBASE = PBASE + "neg/"
else:
    PBASE = PBASE + "comb/"

# EWKSPLIT="-p 'W_fake,W,Z,Top,DiBosons,data'"
EWKSPLIT="-p 'W_fake,W,Z,data'"
EWKEXCLUDE="--xp 'W_LO,Z_LO'"
if addQCDMC:
    EWKSPLIT="-p 'QCD,W,Z,Top,DiBosons,data'"
if reweightZpt: EWKSPLIT = EWKSPLIT.replace("W,Z,","Wpt,Zpt,")

MCEFF = "python " + plotterPath + "w-mass-13TeV/dataFakeRate.py " + OPTIONS + " " + EWKSPLIT + " " + EWKEXCLUDE +" --groupBy cut " + plotterPath + "w-mass-13TeV/make_fake_rates_sels.txt " + plotterPath + "w-mass-13TeV/make_fake_rates_xvars.txt  "
if addQCDMC:
    MCEFF += "--sp QCD "
else:
    MCEFF += "--sp W_fake "

MCEFF += " --sP " + NUM + " --sP " + XVAR + "  --sP " + FITVAR + " " + FITVAR + "  --ytitle 'Fake rate' "
MCEFF += " --fixRatioRange --maxRatioRange 0.7 1.29 "      # ratio for other plots
LEGEND=" --legend=TL --fontsize 0.05 --legendWidth 0.4"
RANGES=" --showRatio  --ratioRange 0.50 1.99 --yrange 0 1.0  --xcut 25 100 "

MCEFF += (LEGEND+RANGES)

if addQCDMC:
    MCGO=MCEFF + " --algo=fQCD --compare QCD_prefit,data_fqcd,data_prefit "
else:
    MCGO=MCEFF + " --algo=fQCD --compare W_fake_prefit,data_fqcd,data_prefit "

for i in range(0,len(etaRange)-1):
    thisRange = etaRange[i].replace(".","p").replace("-","m") + "_" + etaRange[i+1].replace(".","p").replace("-","m")
    if useSignedEta:
        print MCEFF + " -o " + PBASE + "/fr_sub_eta_" + thisRange + ".root --bare -A onelep eta 'LepGood1_eta>" + str(etaRange[i]) + " && LepGood1_eta<" + str(etaRange[i+1]) + "' " + str(chargeSelection) + "\n"
    else:
        print MCEFF + " -o " + PBASE + "/fr_sub_eta_" + thisRange + ".root --bare -A onelep eta 'abs(LepGood1_eta)>" + str(etaRange[i]) + " && abs(LepGood1_eta)<" + str(etaRange[i+1]) + "' " + str(chargeSelection) + "\n"

    print "\n\n"



