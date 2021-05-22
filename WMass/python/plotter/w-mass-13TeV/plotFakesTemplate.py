#!/usr/bin/env python3

# manage histograms made using runFakeRate.py
# it creates the QCD templates including all systematic variations, and also stores the histograms for other process
# in a root file which can be passed as input to makeHistogramsWMass.py
# this script also creates the variations for luminosity, which is a constant scaling for processes other than QCD
# while for QCD the systematics depends of the subtraction of other processes from data

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

sys.path.append(os.getcwd() + "/plotUtils/")
from utility import *
#from plotUtils.utility import *

sys.path.append(os.getcwd())
from cropNegativeTemplateBins import cropNegativeContent

logging.basicConfig(level=logging.INFO)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("rootfile", type=str, nargs=1, help="Input file with histograms")
    parser.add_argument("outdir",   type=str, nargs=1, help="Ouput folder for plots")
    parser.add_argument("-b", "--ptBins", required=True, type=str, help = "Binning for pt in a single region, formatted as in TH1 constructor, e.g. '29,26,55' (uniform binning expected for now)")
    parser.add_argument("-l", "--lumi", type=float, default="1.012", help="Luminosity uncertainty (1 + uncertainty, so e.g. 1.012 for 1.2% uncertainty )")
    parser.add_argument("--plot-fakes-syst", dest="plotFakesSyst", action="append", default=[], help="Plot this syst for fakes, passing it as 'systname=zbin1,zbin2,...'. Can specify multiple times")
    parser.add_argument("--pt-range-projection", dest="ptRangeProjection", default=(0,-1), type=float, nargs=2, help="Pt range to select bins to use for 1D projection (for upper range remember that upper bin edge belongs to next bin in ROOT)")
    args = parser.parse_args()

    if len(args.ptBins.split(',')) != 3:
        print("Error: the pt binning is expected to be of the form 'nPt,ptLow,ptHigh'. Abort")
        quit()
              
    fname = args.rootfile[0]
    outdir = args.outdir[0]
    if not outdir.endswith('/'):
        outdir += '/'
    createPlotDirAndCopyPhp(outdir)

    ROOT.TH1.SetDefaultSumw2()

    plotFakesSyst = {}
    outdirSystRatio = outdir + "systRatios_allRegions/"
    if len(args.plotFakesSyst):
        createPlotDirAndCopyPhp(outdirSystRatio)
        for x in args.plotFakesSyst:
            sname,bins = x.split('=')
            plotFakesSyst[sname] = [int(b) for b in bins.split(',')]

    
    lumi = float(args.lumi)
    #lumiUnc = 100.0*(lumi-1.0)
    #print(f"lumi = {lumi}  --> uncertainty on luminosity {lumiUnc}%")
    systNamesProcsHists = {}
    
    f = ROOT.TFile.Open(fname)
    nomihists = {}
    if not f or not f.IsOpen():
        raise RuntimeError(f"Error when opening file {fname}")
    for k in f.GetListOfKeys():
        name = k.GetName()
        syst,proc = name.split("__")
        if syst == "nominal":
            nomihists[proc] = f.Get(name)
            checkNullObj(nomihists[proc], objName=name)
            nomihists[proc].SetDirectory(0)
            wasCropped = cropNegativeContent(nomihists[proc], silent=True, cropError=False)
        else:
            if syst not in systNamesProcsHists:
                systNamesProcsHists[syst] = {} # this dict will have {proc : hist}
            systNamesProcsHists[syst][proc] = f.Get(name)
            checkNullObj(systNamesProcsHists[syst][proc], objName=name)
            systNamesProcsHists[syst][proc].SetDirectory(0)
            wasCropped = cropNegativeContent(systNamesProcsHists[syst][proc], silent=True, cropError=False)
    f.Close()

    nEtaBins = nomihists["data"].GetNbinsX()
    etaLow   = nomihists["data"].GetXaxis().GetBinLowEdge(1)
    etaHigh  = nomihists["data"].GetXaxis().GetBinLowEdge(1+nEtaBins)
    tokens = args.ptBins.split(',')
    nPtBins = int(tokens[0])
    ptLow   = float(tokens[1])
    ptHigh  = float(tokens[2])

    xAxisName = "Muon #eta"
    yAxisName = "Muon p_{T} (GeV)"

    systNamesProcsHists["luminosity"] = {}
    nPtBins4regions = nomihists["data"].GetNbinsY()
    ptLow4regions   = nomihists["data"].GetYaxis().GetBinLowEdge(1)
    ptHigh4regions  = nomihists["data"].GetYaxis().GetBinLowEdge(1+nPtBins4regions)

    # get data-MC for all regions
    hDataSubMC = nomihists["data"].Clone("dataSubMC")
    allProcnames = []
    for k in nomihists.keys():
        if k == "data": continue
        allProcnames.append(k)
        hDataSubMC.Add(nomihists[k],-1.0)
        # add luminosity Up and Down
        systNamesProcsHists["luminosity"][k] = ROOT.TH3D(f"luminosity__{k}", f"luminosity Up/Down by {args.lumi}%",
                                                         nEtaBins, etaLow, etaHigh,
                                                         nPtBins4regions, ptLow4regions, ptHigh4regions,
                                                         2, 0.5, 2.5)
        fillTH3binFromTH2(systNamesProcsHists["luminosity"][k], nomihists[k], 1, scaleFactor=lumi)
        fillTH3binFromTH2(systNamesProcsHists["luminosity"][k], nomihists[k], 2, scaleFactor=(1./lumi))
        systNamesProcsHists["luminosity"][k].GetXaxis().SetTitle(xAxisName)
        systNamesProcsHists["luminosity"][k].GetYaxis().SetTitle(yAxisName)

    print(f"PROCESSES: {allProcnames}")
        
    hDataSubMC_syst = {} # syst : TH3

    canvas = ROOT.TCanvas("canvas", "", 800, 700)
    canvas4regions = ROOT.TCanvas("canvas4regions", "", 800, 1200)
    adjustSettings_CMS_lumi()

    
    for sys in systNamesProcsHists.keys():
        # get numbr of Z bins from TH3
        anyKey = list(systNamesProcsHists[sys].keys())[0]
        hDataSubMC_syst[sys] = copy.deepcopy(systNamesProcsHists[sys][anyKey].Clone(f"{sys}_dataSubMC"))
        hDataSubMC_syst[sys].Reset("ICESM")
        nZbins = hDataSubMC_syst[sys].GetNbinsZ()
        procsWithSyst = list(systNamesProcsHists[sys].keys()) 
        for iz in range(1, 1+nZbins):
            tmpTH2dataSubMC = copy.deepcopy(nomihists["data"].Clone(f"tmpTH2dataSubMC_{sys}_{iz}"))
            for proc in allProcnames:
                if proc in procsWithSyst:
                    #tmpTH2mc = nomihists["data"].Clone("tmpTH2mc")
                    #fillTH2fromTH3zbin(tmpTH2mc, systNamesProcsHists[sys][proc], iz)
                    tmpTH2mc = getTH2fromTH3(systNamesProcsHists[sys][proc], f"tmpTH2mc_{sys}_{iz}_{proc}", iz)
                    tmpTH2dataSubMC.Add(tmpTH2mc, -1.0)
                    if sys in plotFakesSyst.keys() and iz in plotFakesSyst[sys]:
                        tmpTH2mc.Divide(nomihists[proc])
                        tmpTH2mc.SetTitle(f"{proc}: {sys}, bin {iz}")
                        zmin,zmax = getMinMaxHisto(tmpTH2mc, excludeEmpty=True, sumError=False,
                                                   excludeUnderflow=True, excludeOverflow=True,
                                                   excludeMin=0, excludeMax=None)
                        ratioAbsMax = max(abs(1.0 - zmin), abs(1.0 - zmax))
                        zmin = 1.0 - ratioAbsMax
                        zmax = 1.0 + ratioAbsMax
                                          
                        drawCorrelationPlot(tmpTH2mc,
                                            xAxisName,
                                            yAxisName,
                                            f"ratio: systematics / nominal::{zmin},{zmax}",
                                            f"ratio_{proc}_{sys}_{iz}",
                                            plotLabel="ForceTitle",
                                            outdir=outdirSystRatio,
                                            passCanvas=canvas4regions,
                                            drawOption="COLZ0")
                else:
                    tmpTH2dataSubMC.Add(nomihists[proc], -1.0)
            # not worth checking here, since the signal region (quadrant n.3) will always have negative bins
            #wasCropped = cropNegativeContent(tmpTH2dataSubMC, silent=False, cropError=False)
            #if wasCropped:
            #    print(f"Histogram bins cropped to 0 for syst {syst} and iz = {iz}")
            fillTH3binFromTH2(hDataSubMC_syst[sys], tmpTH2dataSubMC, iz)

            
    # now unpack the four regions within the histogram
    regionId = {0: "highIso_lowMt",
                1: "lowIso_lowMt",
                2: "highIso_highMt",
                3: "lowIso_highMt"}

    hFakes = {}
    
    outdirQCD = outdir + "templatesQCD/"
    createPlotDirAndCopyPhp(outdirQCD)
    
    # nominal is a TH2, systematics are TH3
    for k in regionId.keys():
        hFakes[k] = ROOT.TH2D(f"hFakes_{regionId[k]}", f"data-MC: {regionId[k]}",
                              nEtaBins, etaLow, etaHigh,
                              nPtBins,  ptLow,  ptHigh)

        ptBinOffset = k * nPtBins
        fillTH2fromTH2part(hFakes[k], hDataSubMC,
                           xbinLow=1, ybinLow=1,
                           xbinHigh=nEtaBins, ybinHigh=nPtBins,
                           xoffset=0, yoffset=ptBinOffset)
        # for ieta in range(1, 1 + nEtaBins):
        #     for ipt in range(1, 1 + nPtBins):
        #         content = hDataSubMC.GetBinContent(ieta, ipt + k * nPtBins)
        #         error   = hDataSubMC.GetBinError(  ieta, ipt + k * nPtBins)
        #         hFakes[k].SetBinContent(ieta, ipt, content)
        #         hFakes[k].SetBinError(  ieta, ipt, error)

        # check there are no negative bins
        wasCropped = cropNegativeContent(hFakes[k], silent=False, cropError=False)
        if wasCropped:
            print(f"Histogram cropped to 0 in {regionId[k]} region")
                
        drawCorrelationPlot(hFakes[k], nomihists["data"].GetXaxis().GetTitle(), nomihists["data"].GetYaxis().GetTitle(), "Events",
                            f"yields_dataSubMC_{regionId[k]}", plotLabel="ForceTitle", outdir=outdirQCD,
                            passCanvas=canvas, drawOption="COLZ0")


    templateQCD = copy.deepcopy(hFakes[1].Clone("templateQCD"))
    templateQCD.Divide(hFakes[0])
    templateQCD.Multiply(hFakes[2])
    templateQCD.SetTitle("QCD template")
    drawCorrelationPlot(templateQCD, nomihists["data"].GetXaxis().GetTitle(), nomihists["data"].GetYaxis().GetTitle(), "Events",
                        f"yields_templateQCD", plotLabel="ForceTitle", outdir=outdirQCD,
                        passCanvas=canvas, drawOption="COLZ0")
    drawCorrelationPlot(templateQCD,
                        nomihists["data"].GetXaxis().GetTitle(), nomihists["data"].GetYaxis().GetTitle(), "Absolute uncertainty",
                        f"absUncertainty_templateQCD", plotLabel="ForceTitle", outdir=outdirQCD,
                        passCanvas=canvas, plotError=True, drawOption="COLZ0")
    drawCorrelationPlot(templateQCD,
                        nomihists["data"].GetXaxis().GetTitle(), nomihists["data"].GetYaxis().GetTitle(), "Relative uncertainty",
                        f"relUncertainty_templateQCD", plotLabel="ForceTitle", outdir=outdirQCD,
                        passCanvas=canvas, plotRelativeError=True, drawOption="COLZ0")
    

    # open file to start saving shapes
    foutname = outdir + "wmass_shapes.root"
    fout = ROOT.TFile.Open(foutname, "RECREATE")
    if not fout or not fout.IsOpen():
        raise RuntimeError(f"Error when opening file {foutname}")
    fout.cd()
    templateQCD.Write("nominal__data_fakes")
    # hdata and hmc will be used later for plotting
    hdata = None
    hmc = []
    # get part of the histogram corresponding to signal region, which is number 3 (4th key)
    ptBinOffsetSigRegion = 3 * nPtBins
    for k in nomihists.keys():
        # nomihists[k].Write(nomihists[k].GetName()) ## no, this has the pt range for all the 4 regions, we don't want it
        hSigRegion = ROOT.TH2D(f"hSigRegion_{k}", "",
                               nEtaBins, etaLow, etaHigh,
                               nPtBins,  ptLow,  ptHigh)
        fillTH2fromTH2part(hSigRegion, nomihists[k],
                           xbinLow=1, ybinLow=1,
                           xbinHigh=nEtaBins, ybinHigh=nPtBins,
                           xoffset=0, yoffset=ptBinOffsetSigRegion)
        if k == "data":
            hdata = copy.deepcopy(hSigRegion.Clone("data"))
        else:
            hmc.append(copy.deepcopy(hSigRegion.Clone(f"{k}")))
        hSigRegion.Write(f"nominal__{k}")
    # add also data_fakes to array    
    hmc.append(copy.deepcopy(templateQCD.Clone(f"data_fakes")))

    # do not close file here, more histograms saved later
    
    print("Now dealing with systematic variations for fakes")
        
    # now for the other systs
    hFakes_syst = {}
    for sys in hDataSubMC_syst.keys():

        #print(f"Processing {sys} ...")
        hFakes_syst[sys] = {} # regionId : histogram
        # get numbr of Z bins from TH3
        nZbins = hDataSubMC_syst[sys].GetNbinsZ()
        zLow   = hDataSubMC_syst[sys].GetZaxis().GetBinLowEdge(1)
        zHigh  = hDataSubMC_syst[sys].GetZaxis().GetBinLowEdge(1 + nZbins)
        
        for k in regionId.keys():                
            hFakes_syst[sys][k] = ROOT.TH3D(f"hFakes_{sys}_{regionId[k]}", regionId[k],
                                            nEtaBins, etaLow, etaHigh,
                                            nPtBins,  ptLow,  ptHigh,
                                            nZbins,   zLow,   zHigh)
            ptBinOffset = k * nPtBins
            fillTH3fromTH3part(hFakes_syst[sys][k], hDataSubMC_syst[sys],
                               xbinLow=1, ybinLow=1, zbinLow=1,
                               xbinHigh=nEtaBins, ybinHigh=nPtBins, zbinHigh=nZbins,
                               xoffset=0, yoffset=ptBinOffset, zoffset=0)

            if sys in plotFakesSyst.keys() and k != 3:  # no need to plot signal region here, it is not used
                outdirQCDsys = f"{outdirQCD}systematics/{sys}/"
                createPlotDirAndCopyPhp(outdirQCDsys)
                for zbin in plotFakesSyst[sys]:
                    hfakeSys = getTH2fromTH3(hFakes_syst[sys][k], f"hFakes_{sys}_{regionId[k]}_{zbin}", int(zbin))
                    hfakeSys.SetTitle(f"{regionId[k]}")
                    drawCorrelationPlot(hfakeSys,
                                        hfakeSys.GetXaxis().GetTitle(),
                                        hfakeSys.GetYaxis().GetTitle(),
                                        f"Events ({sys} bin {zbin})",
                                        f"hFakes_{sys}_{regionId[k]}_{zbin}",
                                        plotLabel="ForceTitle",
                                        outdir=outdirQCDsys,
                                        passCanvas=canvas,
                                        drawOption="COLZ0")

            
            # check there are no negative bins
            wasCropped = cropNegativeContent(hFakes[k], silent=False, cropError=False)
            if wasCropped:
                print(f"Histogram cropped to 0 in {regionId[k]} region and syst {sys}")

        templateQCDsys = copy.deepcopy(hFakes_syst[sys][1].Clone(f"templateQCD_{sys}"))
        templateQCDsys.Divide(hFakes_syst[sys][0])
        templateQCDsys.Multiply(hFakes_syst[sys][2])
        templateQCDsys.SetTitle(f"QCD template: {sys}")
        templateQCDsys.Write(f"{sys}__data_fakes")
        print(f"Writing histogram {sys}__data_fakes")

    # get part of the histogram corresponding to signal region, which is number 3 (4th key)
    ptBinOffsetSigRegion = 3 * nPtBins # defined above actually
    for sys in systNamesProcsHists.keys():
        for proc in systNamesProcsHists[sys].keys():
            nZbins = systNamesProcsHists[sys][proc].GetNbinsZ()
            zLow = systNamesProcsHists[sys][proc].GetZaxis().GetBinLowEdge(1)
            zHigh = systNamesProcsHists[sys][proc].GetZaxis().GetBinLowEdge(1 + nZbins)
            histRegion3 = ROOT.TH3D(f"{sys}__{proc}_region3", f"{sys}__{proc}",
                                    nEtaBins, etaLow, etaHigh,
                                    nPtBins,  ptLow,  ptHigh,
                                    nZbins,   zLow,   zHigh)
            fillTH3fromTH3part(histRegion3, systNamesProcsHists[sys][proc],
                               xbinLow=1, ybinLow=1, zbinLow=1,
                               xbinHigh=nEtaBins, ybinHigh=nPtBins, zbinHigh=nZbins,
                               xoffset=0, yoffset=ptBinOffsetSigRegion, zoffset=0)
            wasCropped = cropNegativeContent(histRegion3, silent=False, cropError=False)
            if wasCropped:
                print(f"Histogram cropped to 0 in {regionId[3]} region for process {proc} and syst {sys}")
            histRegion3.Write(f"{sys}__{proc}")
            print(f"Writing histogram {sys}__{proc}")
            
    # now can close the file
    fout.Close()
    print(f"Shapes for all processes written in file {foutname}")
    print()

    
    outdirSR = outdir + "distributions_signalRegion/"
    createPlotDirAndCopyPhp(outdirSR)

    fShapesName = outdirSR + "plots.root"
    fShapes = ROOT.TFile.Open(fShapesName, "RECREATE")
    if not fShapes or fShapes.IsZombie():
        print(f"Error opening file {fShapesName}")
        quit()
        
    # now we may do some more plotting of data and MC in signal region
    hmc = sorted(hmc, key= lambda x: x.Integral()) # , reverse=True) 
    stack_eta = ROOT.THStack("stack_eta", "signal and backgrounds")
    stack_pt = ROOT.THStack("stack_pt", "signal and backgrounds")

    ratio2D = copy.deepcopy(hdata.Clone("dataOverMC2D"))
    den2D = copy.deepcopy(hdata.Clone("sigAndBkg2D"))
    den2D.Reset("ICESM")
    
    hdata.SetMarkerColor(ROOT.kBlack)
    hdata.SetLineColor(ROOT.kBlack)
    #hdata.SetLineWidth(2)
    hdata.SetMarkerStyle(20)
    hdata.SetMarkerSize(1)
    hdata.SetTitle("")

    # for projections along eta
    ptRange = ""
    if args.ptRangeProjection[0] < args.ptRangeProjection[1]:
        lowPtbin = max(1, hdata.GetYaxis().FindFixBin(args.ptRangeProjection[0]))
        highPtbin = min(hdata.GetNbinsY(), hdata.GetYaxis().FindFixBin(args.ptRangeProjection[1])) # hdata.GetNbinsY()
        ptRange = "_%gTo%g" % (hdata.GetYaxis().GetBinLowEdge(lowPtbin), hdata.GetYaxis().GetBinLowEdge(1+highPtbin))
        ptRange = ptRange.replace(".","p")
    else:
        lowPtbin = 0
        highPtbin = -1
    
    hdata_eta = hdata.ProjectionX("data_eta",lowPtbin,highPtbin,"e")
    hdata_pt  = hdata.ProjectionY("data_pt",0,-1,"e")

    legend = ROOT.TLegend(0.2,0.72,0.95,0.92)
    legend.SetFillColor(0)
    legend.SetFillStyle(0)
    legend.SetBorderSize(0)
    legend.SetNColumns(3)

    colors = {"Wmunu"      : ROOT.kRed+1,
              "Wmunu_plus" : ROOT.kRed+2,
              "Wmunu_minus": ROOT.kRed+1,
              "Zmumu"      : ROOT.kAzure+2,
              "Wtaunu"     : ROOT.kCyan+1,
              "Wtaunu_plus" : ROOT.kCyan+2,
              "Wtaunu_minus" : ROOT.kCyan+1,
              "Ztautau"    : ROOT.kSpring+9,
              "Top"        : ROOT.kGreen+2,
              "Diboson"    : ROOT.kViolet,
              "data_fakes" : ROOT.kGray}

    legEntries = {"Wmunu"      : "W#rightarrow#mu#nu",
                  "Wmunu_plus" : "W^{+}#rightarrow#mu#nu",
                  "Wmunu_minus": "W^{-}#rightarrow#mu#nu",
                  "Zmumu"      : "Z#rightarrow#mu#mu",
                  "Wtaunu"     : "W#rightarrow#tau#nu",
                  "Wtaunu_plus" : "W^{+}#rightarrow#tau#nu",
                  "Wtaunu_minus": "W^{-}#rightarrow#tau#nu",
                  "Ztautau"    : "Z#rightarrow#tau#tau",
                  "Top"        : "t quark",
                  "Diboson"    : "Diboson",
                  "data_fakes" : "QCD"}

    hdata.Write("data2D")
    hdata_eta.Write()
    hdata_pt.Write()
    
    legend.AddEntry(hdata, "Data", "EP")
    for h in hmc:
        h.SetTitle("")
        h.SetFillColor(colors[h.GetName()])
        h.SetLineColor(ROOT.kBlack)
        stack_eta.Add(h.ProjectionX(f"{h.GetName()}_eta",lowPtbin,highPtbin,"e"))
        stack_pt.Add( h.ProjectionY(f"{h.GetName()}_pt",0,-1,"e"))
        den2D.Add(h)
        h.Write()
    for i in range(len(hmc)-1, 0, -1):
        legend.AddEntry(hmc[i], legEntries[hmc[i].GetName()], "F")

    stack_eta.Write()
    stack_pt.Write()
    den2D.Write()

    ratio2D.Divide(den2D)
    ratio2D.Write()
    
    canvas1D = ROOT.TCanvas("canvas1D", "", 800, 900)
    
    drawTH1dataMCstack(hdata_eta, stack_eta, "Muon #eta", "Events", "muon_eta_signalRegion" + ptRange,
                       outdirSR, legend, ratioPadYaxisNameTmp="Data/MC::0.92,1.08", passCanvas=canvas1D,
                       drawLumiLatex=True, xcmsText=0.3, noLegenRatio=True
    )
    drawTH1dataMCstack(hdata_pt, stack_pt, "Muon p_{T} (GeV)", "Events", "muon_pt_signalRegion",
                       outdirSR, legend, ratioPadYaxisNameTmp="Data/MC::0.92,1.08", passCanvas=canvas1D,
                       drawLumiLatex=True, xcmsText=0.3, noLegenRatio=True
    )

    ratio2D.SetTitle("data / (signal + background)")
    drawCorrelationPlot(ratio2D, xAxisName, yAxisName, "data/MC",
                        f"muon_eta_pt_dataMCratio", plotLabel="ForceTitle", outdir=outdirSR,
                        palette=57, passCanvas=canvas, drawOption="COLZ0", skipLumi=True)
    drawCorrelationPlot(ratio2D, xAxisName, yAxisName, "data/MC",
                        f"muon_eta_pt_dataMCratio_absUncertainty", plotLabel="ForceTitle", outdir=outdirSR,
                        palette=57, passCanvas=canvas, drawOption="COLZ0", skipLumi=True, plotError=True)

    
    allHists = hmc + [hdata]
    for h in allHists:
        h.SetTitle(h.GetName())
        drawCorrelationPlot(h, xAxisName, yAxisName, "Events",
                            f"muon_eta_pt_{h.GetName()}", plotLabel="ForceTitle", outdir=outdirSR,
                            palette=57, passCanvas=canvas, drawOption="COLZ0", skipLumi=True)
