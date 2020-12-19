#!/usr/bin/env python
from shutil import copyfile
import re, sys, os, os.path, subprocess
import ROOT

from optparse import OptionParser
parser = OptionParser(usage="%prog ")
parser.add_option("--qcdmc", dest="addQCDMC", default=False, action='store_true', help="Add QCD MC in plots (but do not subtract from data)");
parser.add_option("-d", "--dry-run", dest="dryRun", default=False, action='store_true', help="Do not execute commands, just print them");
parser.add_option("--charge", dest="charge", default="", type='string', help="Select charge: p for positive, n for negative");
parser.add_option("--outdir", dest="outdir", default="./", type='string', help="Output folder");
parser.add_option("--tree-path", dest="treePath", default="", type='string', help="Path to trees on eos");
parser.add_option("--lumi", dest="lumi", default=35.9 , type='float', help="Integrated luminosity");
parser.add_option("--pt", dest="ptvar", default="pt", type='string', help="Select pT definition from make_fake_rates_xvars.txt");
parser.add_option("--no-scaleFactors", dest="noScaleFactors", default=False, action='store_true', help="Don't use lepton scale factors (only PU weight)");
parser.add_option("--useSignedEta", dest="useSignedEta", default=False, action='store_true', help="Make fake rate for eta bins distinguishing eta sign");
parser.add_option("--addOpts", dest="addOpts", default="", type='string', help="Options to pass some other options from outside to build the command");
parser.add_option("--reweightZpt", dest="reweightZpt", default=False, action='store_true', help="Use W and Z with reweighted pT");
(options, args) = parser.parse_args()

addQCDMC = options.addQCDMC  # trying to add QCD MC to graphs to be compared
charge = str(options.charge)
ptvar = str(options.ptvar)
useSignedEta = options.useSignedEta
addOpts = options.addOpts
luminosity = options.lumi
reweightZpt = options.reweightZpt

if useSignedEta:
    fitvar = "eta" 
    etaRange = [ '-2.4', '2.4']
else:
    fitvar = "abseta" 
    etaRange = [ '0.0', '2.4']

plotterPath = str(os.environ.get('CMSSW_BASE'))
plotterPath = plotterPath + "/src/CMGTools/WMass/python/plotter/"

chargeSelection = ""
if charge != "":
    if charge == "p":
        chargeSelection = "-A accept positive 'Muon_charge[0] > 0'"
    elif charge == "n":
        chargeSelection = "-A accept negative 'Muon_charge[0] < 0'"
    else:
        print "## %s is not a valid input for charge setting: use p or n" % charge
        quit()

T = options.treePath
objName = 'Events'
print "## used trees from: ",T

MCweightOption = ' -W "puWeight*PrefireWeight*_get_muonSF_fast_wmass(Muon_pt[0],Muon_eta[0],Muon_charge[0])" ' 
if options.noScaleFactors:
    print "## Warning: not using lepton scale factors: only PU weight"
    MCweightOption = ' -W "puWeight*PrefireWeight" '

J=4

BASECONFIG=plotterPath + "w-mass-13TeV/testingNano/cfg"
MCA=BASECONFIG+'/mca-fakeRate.txt'
CUTFILE=BASECONFIG+'/cuts_fakeRate.txt'
XVAR=ptvar
FITVAR=fitvar
NUM = "fakeRateNumerator"
IDFILE = plotterPath + "w-mass-13TeV/make_fake_rates_sels.txt "
XVARFILE =  plotterPath + "w-mass-13TeV/make_fake_rates_xvars.txt "

OPTIONS = MCA + " " + CUTFILE + " " + IDFILE + XVARFILE 
OPTIONS += " -f -P " + T + " --obj " + objName + " --s2v -j " + str(J) + " -l " + str(luminosity) + " " + str(addOpts) + " " + MCweightOption

if charge == "p":
    WPROC = " --pg 'Wmunu := Wmunu_plus' --pg 'Wtaunu := Wtaunu_plus' "
elif charge == "m":
    WPROC = " --pg 'Wmunu := Wmunu_minus' --pg 'Wtaunu := Wtaunu_minus' "
else:
    WPROC = " --pg 'Wmunu := Wmunu_plus,Wmunu_minus' --pg 'Wtaunu := Wtaunu_plus,Wtaunu_minus' "
DATAPROC = " --pg 'data := data_B,data_C,data_D,data_E,data_F,data_F_postVFP,data_G,data_H' "

EWKSPLIT = "-p 'Wmunu,Wtaunu,Zmumu,Ztautau,Top,DiBosons,data' " + DATAPROC + WPROC

if addQCDMC:
    EWKSPLIT += " -p QCD --sp QCD "
#if reweightZpt: EWKSPLIT = EWKSPLIT.replace("W,Z,","Wpt,Zpt,") # FIXME

MCEFF = "python " + plotterPath + "w-mass-13TeV/dataFakeRateLite.py " + OPTIONS + " " + EWKSPLIT

MCEFF += " --sP " + NUM + " --sP " + XVAR + "  --sP " + FITVAR + " --fitVar " + FITVAR

thisRange = etaRange[0].replace(".","p").replace("-","m") + "_" + etaRange[-1].replace(".","p").replace("-","m")
fout = "/fr_sub_eta_" + thisRange + ".root "
if useSignedEta:
    cut = " -A accept eta 'Muon_eta[0]>" + str(etaRange[0]) + " && Muon_eta[0]<" + str(etaRange[-1]) + "' " + str(chargeSelection)
else:
    cut = " -A accept eta 'abs(Muon_eta[0])>" + str(etaRange[0]) + " && abs(Muon_eta[0])<" + str(etaRange[-1]) + "' " + str(chargeSelection)

finalCommand = MCEFF + " -o " + options.outdir + fout + cut + "\n"

# copy some stuff for reference    
os.system("cp {f} {o}".format(f=MCA,o=options.outdir))
os.system("cp {f} {o}".format(f=CUTFILE,o=options.outdir))
try:
    with open(options.outdir+'/command.txt', 'w') as f:
        f.write(finalCommand+'\n')
except IOError as e:
    print("## Couldn't open or write to file %s (%s)." % (options.outdir+'/command.txt', e))

if options.dryRun:
    print finalCommand
else:
    print "## Executing the command\n"    
    print finalCommand
    print "\n\n"
    os.system(finalCommand)
