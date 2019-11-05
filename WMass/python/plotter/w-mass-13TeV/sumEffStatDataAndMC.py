import os, sys, math
import numpy as np
import ROOT
ROOT.gROOT.SetBatch(True)

sys.path.append(os.getcwd() + "/plotUtils/")
from utility import *


from optparse import OptionParser
parser = OptionParser(usage='%prog [options]')
parser.add_option('-f','--flavor', dest='flavor', default='', type='string', help='name of the channel (mu or el)')
parser.add_option('-c','--charge', dest='charge', default='', type='string', help='plus or minus if the efficiencies were derived separately for each charge (should be for muons only. If empty, assumes no charge splitting in the inputs')
parser.add_option(     '--eleEEeta0p1', dest='eleEEeta0p1',action="store_true", default=False, help='Only for electrons, to choose between two possible files. It True, use the one with EE made with 0.1 eta granularity')
(options, args) = parser.parse_args()

ROOT.TH1.SetDefaultSumw2()

isMu = True if options.flavor == "mu" else False
charge = options.charge
lepton = "muon" if isMu else "electron" 
eleEEeta0p1 = options.eleEEeta0p1

if isMu:
    inputpath = "plots/scaleFactors_Final/effSyst_fromRooFitResult_onlyStatUnc_mu{ch}/".format(ch="Minus" if charge == "minus" else "Plus")
    fdata = inputpath + "systEff_trgmu_%s_OnlyStatUnc_forData_mu.root" % charge
    fmc =   inputpath + "systEff_trgmu_%s_OnlyStatUnc_forMC_mu.root"   % charge
    fout = inputpath + "systEff_trgmu_%s_OnlyStatUnc_quadSum_mu.root"   % charge
else:
    if eleEEeta0p1:
        inputpath = "plots/scaleFactors_Final/effSyst_fromRooFitResult_onlyStatUnc_el_etaEE0p1/"
    else:
        inputpath = "plots/scaleFactors_Final/effSyst_fromRooFitResult_onlyStatUnc_el_etaEB0p1_etaEE0p2/"
    fdata = inputpath + "systEff_trgel_OnlyStatUnc_forData_el.root" 
    fmc =   inputpath + "systEff_trgel_OnlyStatUnc_forMC_el.root"   
    fout = inputpath + "systEff_trgel_OnlyStatUnc_quadSum_el.root"   

tfdata = ROOT.TFile.Open(fdata)
if not tfdata or not tfdata.IsOpen():
    print "Error opening file %s" % fdata
    quit()
tfmc = ROOT.TFile.Open(fmc)
if not tfmc or not tfmc.IsOpen():
    print "Error opening file %s" % fmc
    quit()

tfout = ROOT.TFile.Open(fout,"recreate")
if not tfout or not tfout.IsOpen():
    print "Error opening file %s" % fout
    quit()

hdata = {}
hmc = {}
hQuadSum = {}

for name in ["p0", "p1", "p2"]:
    hdata[name] = tfdata.Get(name)
    hdata[name].SetDirectory(0)
    hmc[name] = tfmc.Get(name)
    hmc[name].SetDirectory(0)
    hQuadSum[name] = hdata[name].Clone("quadSum_"+name)
    hQuadSum[name].SetTitle("Nuisance for Erf %s" % name)
    hQuadSum[name].Reset("ICESM")
    for ix in range(1,hQuadSum[name].GetNbinsX()+1):
        for iy in range(1,hQuadSum[name].GetNbinsY()+1):
            pdata = hdata[name].GetBinContent(ix,iy)
            pmc = hmc[name].GetBinContent(ix,iy)
            quadSum = math.sqrt(pdata * pdata + pmc * pmc)
            # use sign of data
            #pmax = pdata if abs(pdata) > abs(pmc) else pmc
            #sign = -1.0 if pmax < 0 else 1.0  # choose sign of the largest of the two in absolute value, not sure it makes sense, but sort of
            sign = -1.0 if pmc < 0 else 1.0
            hQuadSum[name].SetBinContent(ix, iy, sign * quadSum)

    drawCorrelationPlot(hQuadSum[name],"%s #eta" % lepton,"%s p_{T} [GeV]" % lepton,
                        "variation / nominal" + ("::-0.004,0.001" if name == "p0" else "::-0.01,0.01"),
                        hQuadSum[name].GetName(),"ForceTitle",inputpath,0,0,False,False,False,
                        1,palette=55)
    hQuadSum[name].Write(name)

tfout .Close()
tfdata.Close()
tfmc.Close()
