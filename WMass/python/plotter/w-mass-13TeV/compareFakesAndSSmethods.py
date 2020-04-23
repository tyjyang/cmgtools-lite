import ROOT, os, datetime, re, operator, math
from array import array
ROOT.gROOT.SetBatch(True)

#import utilities
#utilities = utilities.util()

import sys
sys.path.append(os.getcwd() + "/plotUtils/")
from utility import *

plotterPath = os.environ['CMSSW_BASE'] + "/src/CMGTools/WMass/python/plotter/"

fFR = plotterPath + "plots/Wlike/TREE_4_WLIKE_MU/comparisons/chPlus_oddEvts_withPrefire_trigSFonlyZ__1ifBothLepMatchTrig_elseOnPlusLepIfMatchElseOther_addTkMuTrigger_testWlikeAnalysis_mt45_compareFakesPredictions_addFRpol2/plots_zmm.root"
fSS = plotterPath + "plots/Wlike/TREE_4_WLIKE_MU/comparisons/chPlus_oddEvts_withPrefire_trigSFonlyZ__1ifBothLepMatchTrig_elseOnPlusLepIfMatchElseOther_addTkMuTrigger_testWlikeAnalysis_mt45_testSameSignLep_v3/plots_zmm.root"

outdir = plotterPath + "plots/Wlike/TREE_4_WLIKE_MU/comparisons/shapesFRandSSmethod/"


if __name__ == "__main__":

    from optparse import OptionParser
    parser = OptionParser(usage='%prog [options] ')
    #parser.add_option('-o','--outdir',     dest='outdir',     default='',   type='string', help='outdput directory save plots')
    (options, args) = parser.parse_args()

    # outdir = options.outdir
    if outdir:
        if not outdir.endswith("/"): outdir += "/"
        if not os.path.isdir(outdir):
            os.system('mkdir -p {od}'.format(od=outdir))
        os.system('cp {pf} {od}'.format(pf='/afs/cern.ch/user/g/gpetrucc/php/index.php',od=outdir))
    else:
        print "Warning: must pass output folder with option -o. Abort"
        quit()

    # no index: FR/PR from Marc
    # v2: with lepton SF, no Zpt, fit up to 55 for FR, fit PR with pol1
    # v3: no lepton SF, no Zpt, fit up to 55 for FR, fit PR with pol1
    # v4: with lepton SF, with Zpt, fit up to 60 for FR, fit both FR/PR with pol2
    # v5: with lepton SF, no Zpt, fit up to 60 for FR, fit both FR/PR with pol2
    procsFR = ["data_fakes", "data_fakes_v2", "data_fakes_v5"]  
    procsSS = ["SS"]

    plots = ["zmassLogy", "zmass", "mt_wlike", "pt", "eta"] # for pt and eta I need to adapt the name for each file
    xaxisLabels = ["mass of muon pair [GeV]",
                   "mass of muon pair [GeV]",
                   "W-like transverse mass [GeV]",
                   "muon p_{T} {GeV",
                   "muon #eta"]

    adjustSettings_CMS_lumi()
    c = ROOT.TCanvas("c","c",900,1000)

    # open the 2 root file
    tfFR = ROOT.TFile.Open(fFR,'read')
    tfSS = ROOT.TFile.Open(fSS,'read')

    for ip,pl in enumerate(plots):
        hfr = {}
        hss = {}
        ## file FR
        tfFR.cd()
        for name in procsFR:
            pname = pl
            if pl == "pt": 
                pname = "ptl1plus"
            elif pl == "eta":
                pname = "etal1plus"                
            fullPlotName =pname+"_"+name
            hfr[name] = tfFR.Get(fullPlotName)
            if not hfr[name]:
                print "Warning: could not get histogram " + fullPlotName
                quit()
            hfr[name].SetDirectory(0)
        ## file SS
        tfSS.cd()
        for name in procsSS:
            pname = pl
            if pl == "pt": 
                pname = "ptBothLep"
            elif pl == "eta":
                pname = "etaBothLep"         
            hsstmp = {}
            for tmp in ["data","signal"]:
                fullPlotName = pname+"_"+tmp
                hsstmp[tmp] = tfSS.Get(fullPlotName)
                if not hsstmp[tmp]:
                    print "Warning: could not get histogram " + fullPlotName
                    quit()                    
                hsstmp[tmp].SetDirectory(0)
            hss[name] = hsstmp["data"].Clone()
            hss[name].Add(hsstmp["signal"],-1.0)
            hss[name].SetDirectory(0)
            # rescale histograms filled with both lepton by dividing by 2
            #if pl in ["pt", "eta"]: 
            #    hss[name].Scale(0.5)        

        # now plotting, using ss at denominator
        hists = []
        legEntries = []
        for n in hss.keys():
            hists.append(hss[n])
            legEntries.append(n)
        for n in hfr.keys():
            hists.append(hfr[n])
            legEntries.append(n)
        drawNTH1(hists,legEntries,xaxisLabels[ip],"Events",pl,outdir,
                 draw_both0_noLog1_onlyLog2=1,
                 labelRatioTmp="FR/SS::0.5,1.5",
                 legendCoords="0.15,0.85,0.75,0.9;2",passCanvas=c)


    tfFR.Close()
    tfSS.Close()

    print ""
    print "THE END"
    print ""

