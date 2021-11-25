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

ROOT.gInterpreter.ProcessLine(".O3")
ROOT.gInterpreter.ProcessLine('#include "ccFiles/thndHelpers.cc"')

## this assumes the dimensions of the THn is eta-pt-charge-isoMt bin

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("rootfile", type=str, nargs=1, help="Input file with histograms")
    parser.add_argument("-l", "--lumi", type=float, default="1.012", help="Luminosity uncertainty (1 + uncertainty, so e.g. 1.012 for 1.2% uncertainty )")
    parser.add_argument("--plot-fakes-syst", dest="plotFakesSyst", action="append", default=[], help="Plot this syst for fakes, passing it as 'systname=zbin1,zbin2,...'. Can specify multiple times")
    parser.add_argument("--pt-range-projection", dest="ptRangeProjection", default=(0,-1), type=float, nargs=2, help="Pt range to select bins to use for 1D projection (for upper range remember that upper bin edge belongs to next bin in ROOT)")
    parser.add_argument(     '--crop-negative', dest='cropNegativeBin' , default=False , action='store_true',   help='Set negative bins to non negative value (it uses cropNegativeContent), but it might be slow')
    parser.add_argument('--makePlots', action='store_true', help="Output pretty plots")
    args = parser.parse_args()
           
    fname = args.rootfile[0]
    outdir = os.path.dirname(fname) + "/postprocessing/"
    createPlotDirAndCopyPhp(outdir)

    ROOT.TH1.SetDefaultSumw2()

    plotFakesSyst = {}
    if len(args.plotFakesSyst):
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
        print(f">>> Reading {name}")
        #if any(x in syst for x in ["pdf", "effStat", "qcdScale"]): continue # to test code in a faster way
        if syst == "nominal":
            nomihists[proc] = safeGetObject(f, name, detach=False)
            checkNullObj(nomihists[proc], objName=name)
            if args.cropNegativeBin:
                wasCropped = cropNegativeContent(nomihists[proc], silent=True, cropError=False)
        else:
            if syst not in systNamesProcsHists:
                systNamesProcsHists[syst] = {} # this dict will have {proc : hist}
            systNamesProcsHists[syst][proc] = safeGetObject(f, name, detach=False)
            checkNullObj(systNamesProcsHists[syst][proc], objName=name)
            if args.cropNegativeBin:
                wasCropped = cropNegativeContent(systNamesProcsHists[syst][proc], silent=True, cropError=False)
    f.Close()
    
    nNomiDim = nomihists["data"].GetNdimensions()
    # nEtaBins = nomihists["data"].GetAxis(0).GetNbins() # nbins with no underflow/overflow, usually 48 eta bins
    # etaLow   = nomihists["data"].GetAxis(0).GetBinLowEdge(1)
    # etaHigh  = nomihists["data"].GetAxis(0).GetBinLowEdge(1+nEtaBins)
    # nPtBins = nomihists["data"].GetAxis(1).GetNbins()
    # ptLow   = nomihists["data"].GetAxis(1).GetBinLowEdge(1)
    # ptHigh  = nomihists["data"].GetAxis(1).GetBinLowEdge(1+nPtBins)

    nominalNbins = [nomihists["data"].GetAxis(i).GetNbins() for i in range(nNomiDim)]
    nominalBinMin = [nomihists["data"].GetAxis(i).GetBinLowEdge(1) for i in range(nNomiDim)]
    nominalBinMax = [nomihists["data"].GetAxis(i).GetBinLowEdge(1+nomihists["data"].GetAxis(i).GetNbins()) for i in range(nNomiDim)]
    
    xAxisName = nomihists["data"].GetAxis(0).GetTitle() # "Muon #eta"
    yAxisName = nomihists["data"].GetAxis(1).GetTitle() # "Muon p_{T} (GeV)"

    ###
    # will use ROOT.fillTHNplus1fromTHn() to fill TH5 from TH4, for the opposite one can use projections
    ###
    
    # get data-MC for all regions
    # also use this loop to create the THn for luminosity, stacking both variations in the same histogram
    hDataSubMC = nomihists["data"].Clone("dataSubMC")
    allProcnames = []
    systNamesProcsHists["luminosity"] = {}
    for k in nomihists.keys():
        if k == "data": continue
        allProcnames.append(k)
        hDataSubMC.Add(nomihists[k], -1.0)
        # add luminosity variations in the list of systematics (with Up and Down variations in the same object)
        systNamesProcsHists["luminosity"][k] = ROOT.THnD("luminosity_{k}", "luminosity_{k}", nNomiDim+1, array('i',nominalNbins+[2]), array('d',nominalBinMin+[0.5]), array('d',nominalBinMax+[2.5]))
        ROOT.fillTHNplus1fromTHn(systNamesProcsHists["luminosity"][k], nomihists[k], 1, 1)
        ROOT.fillTHNplus1fromTHn(systNamesProcsHists["luminosity"][k], nomihists[k], 2, 2)
        
    print(f"PROCESSES: {allProcnames}")

    hDataSubMC_syst = {} # syst : TH5 (but may have TH4 for some, only eta-pt-charge-isoMt)
    
    canvas = ROOT.TCanvas("canvas", "", 800, 700)
    adjustSettings_CMS_lumi()

    print()
    print(">>> Now going to produce histograms for systematic variations of fakes")
    print()
    for sys in systNamesProcsHists.keys():
        print(f"   {sys}")
        # get numbr of Z bins from TH3
        procsWithSyst = list(systNamesProcsHists[sys].keys())
        anyKey = procsWithSyst[0]
        tmpSystHisto = systNamesProcsHists[sys][anyKey] # just to facilitate usage below
        systHasSameDimAsNomi = True if tmpSystHisto.GetNdimensions() == nNomiDim else False
        if systHasSameDimAsNomi:
            hDataSubMC_syst[sys] = nomihists["data"].Clone(f"{sys}_dataSubMC")
            for proc in allProcnames:
                if proc in procsWithSyst:
                    hDataSubMC_syst[sys].Add(systNamesProcsHists[sys][proc], -1.0)
                else:
                    hDataSubMC_syst[sys].Add(nomihists[proc], -1.0)
        else:
            nZbins = tmpSystHisto.GetAxis(nNomiDim).GetNbins() # number of bins on the axis for systematics
            nSystDim = nNomiDim + 1  # axis index for the systematic variations (usually 5th dimension)
            # the following are the number of bins, min, and max, for all axes, the name of the variable starts with syst just because these are the histograms embedding systematics
            systNbins  = [tmpSystHisto.GetAxis(i).GetNbins() for i in range(nSystDim)]
            systBinMin = [tmpSystHisto.GetAxis(i).GetBinLowEdge(1) for i in range(nSystDim)]
            systBinMax = [tmpSystHisto.GetAxis(i).GetBinLowEdge(1+nZbins) for i in range(nSystDim)]
            # create histogram for data with one more dimension and same bins as systematic histograms
            # same for other nominal histograms for processes not having this systematic uncertainty
            nomihistDimNp1 = {}
            print(f"      Augmenting data in N+1 dimensions")
            hDataSubMC_syst[sys] = ROOT.THnD(f"{sys}_dataSubMC", "", nSystDim, array('i',systNbins), array('d',systBinMin), array('d',systBinMax))
            ROOT.fillTHNplus1fromTHn(hDataSubMC_syst[sys], nomihists["data"], 0, nZbins+1)
            for proc in allProcnames:
                print(f"      Subtracting {proc}")
                if proc in procsWithSyst:
                    hDataSubMC_syst[sys].Add(systNamesProcsHists[sys][proc], -1.0)
                else:
                    nomihistDimNp1 = ROOT.THnD(f"nominal_{proc}_dimNp1", "", nSystDim, array('i',systNbins), array('d',systBinMin), array('d',systBinMax))
                    print(f"      Augmenting {proc} in N+1 dimensions")
                    ROOT.fillTHNplus1fromTHn(nomihistDimNp1, nomihists[proc], 0, nZbins+1)
                    hDataSubMC_syst[sys].Add(nomihistDimNp1, -1.0)

    
    print()
    print(">>> Now building histogram for fakes in signal region")
    print()
    # now unpack the four regions within the histogram

    ## should be careful about region used at denominator in the fake rate factor, if a bin is close to 0 or negative then cropping it to 0.0001 will result in a huge ratio if numerator is not 0
    # perhaps in those cases I should set the numerator to 0, or crop to 0 such that THn::Divide sets the ratio to 0 when denominator is 0

    # the id starts from 0 because it was associated to the bin center (axis defined with 4 bins from -0.5 to 3.5)
    # but now we have isoMt in a dedicated axis rather than by multiplying the pt axis by 4 as we were doing before
    # so we have to pick the bin number 4 to get the signal region (root axis has bin number starting from 1)
    regionId = {0: "highIso_lowMt",
                1: "lowIso_lowMt",
                2: "highIso_highMt",
                3: "lowIso_highMt"}
    # output will be TH3 with eta-pt-charge, not sure if we want to keep that as a TH3 of THn with dimension 3
    hFakes = {}
    
    outdirQCD = outdir + "templatesQCD/"
    createPlotDirAndCopyPhp(outdirQCD)
    
    axesForProjection = [i for i in range(nNomiDim-1)] # if nNomiDim = 4 we keep 3 axes, so range(3) = 0, 1, 2
    for k in regionId.keys():
        hDataSubMC.GetAxis(nNomiDim-1).SetRange(k+1, k+1) # select only one bin: note that region 0 is bin 1, region 1 is bin 2 and so on
        #hFakes[k] = hDataSubMC.ProjectionAny(len(axesForProjection), array("i", axesForProjection), False, "E") # protected method
        hFakes[k] = hDataSubMC.Projection(0, 1, 2, "E")
        hFakes[k].SetName(f"hFakes_{regionId[k]}")  # here all charges are in the same TH3
        hFakes[k].SetTitle(f"data-MC: {regionId[k]}")
        if args.cropNegativeBin:
            wasCropped = cropNegativeContent(hFakes[k], silent=False, cropError=False)
            if wasCropped:
                print(f"Histogram {hFakes[k].GetName()} cropped to 0 in {regionId[k]} region")

        for chbin in range(1,3):
            charge = "plus" if chbin == 2 else "minus"
            hFakes2D = getTH2fromTH3(hFakes[k], f"hFakes_{regionId[k]}_charge{charge}", chbin, chbin)
            hFakes2D.SetTitle(f"data-MC ({charge}): {regionId[k]}")
            if args.makePlots:
                drawCorrelationPlot(hFakes2D, xAxisName, yAxisName, "Events",
                                    f"yields_dataSubMC_{regionId[k]}_charge{charge}", plotLabel="ForceTitle", outdir=outdirQCD,
                                    passCanvas=canvas, drawOption="COLZ0")
        
    templateQCD = hFakes[1].Clone("templateQCD")
    templateQCD.Divide(hFakes[0])
    fakerateFactor = templateQCD.Clone("fakerateFactor")
    fakerateFactor.SetTitle("lowIso_lowMt / highIso_lowMt ")
    templateQCD.Multiply(hFakes[2])
    templateQCD.SetTitle("QCD template")

    if args.makePlots:
        for chbin in range(1,3):
            charge = "plus" if chbin == 2 else "minus"
            fakerateFactor2D = getTH2fromTH3(fakerateFactor, f"{fakerateFactor.GetName()}_charge{charge}", chbin, chbin)
            fakerateFactor2D.SetTitle(f"lowIso_lowMt / highIso_lowMt, charge {charge}")
            templateQCD2D = getTH2fromTH3(templateQCD, f"{templateQCD.GetName()}_charge{charge}", chbin, chbin)
            templateQCD2D.SetTitle(f"QCD template, charge {charge}")
            drawCorrelationPlot(fakerateFactor2D,
                                xAxisName, yAxisName, "Fakerate factor",
                                f"fakerateFactor_charge{charge}", plotLabel="ForceTitle", outdir=outdirQCD,
                                passCanvas=canvas, drawOption="COLZ0")
            drawCorrelationPlot(templateQCD2D, xAxisName, yAxisName, "Events",
                                f"yields_templateQCD_charge{charge}", plotLabel="ForceTitle", outdir=outdirQCD,
                                passCanvas=canvas, drawOption="COLZ0")
            drawCorrelationPlot(templateQCD2D,
                                xAxisName, yAxisName, "Absolute uncertainty",
                                f"absUncertainty_templateQCD_charge{charge}", plotLabel="ForceTitle", outdir=outdirQCD,
                                passCanvas=canvas, plotError=True, drawOption="COLZ0")
            drawCorrelationPlot(templateQCD2D,
                                xAxisName, yAxisName, "Relative uncertainty::0.0,1.0",
                                f"relUncertainty_templateQCD_charge{charge}", plotLabel="ForceTitle", outdir=outdirQCD,
                                passCanvas=canvas, plotRelativeError=True, drawOption="COLZ0")

        # open utility file to save QCD templates in all regions
        fqcdName = outdirQCD + "histogramsQCD.root"
        fqcd = safeOpenFile(fqcdName, mode="RECREATE")
        for k in regionId.keys():
            hFakesTmp = hFakes[k].Clone(f"tmp{k}")
            hFakesTmp.Write(f"yields_dataSubMC_{regionId[k]}")
        fakerateFactor.Write()
        templateQCDtoWrite = templateQCD.Clone(f"yields_templateQCD")
        templateQCDtoWrite.Write(f"yields_templateQCD")
        fqcd.Close()
    
    #########
    
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
    # get part of the histogram corresponding to signal region, which is number 4 (4th bin of the 4th axis of the TH4)
    for k in nomihists.keys():
        nomihists[k].GetAxis(nNomiDim-1).SetRange(4, 4)
        hSigRegion = nomihists[k].Projection(0, 1, 2, "E")
        hSigRegion.SetName(f"hSigRegion_{k}")
        hSigRegion.SetTitle("")
        hSigRegion.SetName(k)
        if k == "data":
            hdata = hSigRegion
        else:
            hmc.append(hSigRegion)
        hSigRegion.Write(f"nominal__{k}")
    # add also data_fakes to array
    templateQCD.SetName("data_fakes")
    hmc.append(templateQCD)

    # do not close file here, more histograms saved later

    print()
    print(">>> Now dealing with systematic variations for fakes")
    print()
        
    # now for the other systs
    hFakes_syst = {}
    for sys in hDataSubMC_syst.keys():

        #print(f"Processing {sys} ...")
        hFakes_syst[sys] = {} # regionId : histogram
        
        for k in regionId.keys():                
            hDataSubMC_syst[sys].GetAxis(nNomiDim-1).SetRange(k+1, k+1) # select only one bin: note that region 0 is bin 1, region 1 is bin 2 and so on (reminder: axis nNomiDim-1 is the isoMt axis, usually the 4th, so index 3)
            
            if hDataSubMC_syst[sys].GetNdimensions() == nNomiDim:     
                # this is the case where TH4 goes into a TH3, in general we should not assume that, namely we should always use ProjectionND to return a THn, and then in case convert to a TH3 when the dimension is correct
                hFakes_syst[sys][k] = hDataSubMC_syst[sys].Projection(0, 1, 2, "E") # returns a TH3
            else:
                # we have the syst axis too here, so need to pick axes 0, 1, 2, 4, while for the isoMt axis (axis index 3) we need to project in the desired isoMt bin, which is the third (by chance, nothing to do with the axis number) 
                axesIds = array("i", [i for i in range(nNomiDim+1) if i != (nNomiDim-1)]) # isoMt is the last axis of the nominal histograms
                hFakes_syst[sys][k] = hDataSubMC_syst[sys].ProjectionND(len(axesIds), axesIds, "E")

            hFakes_syst[sys][k].SetName(f"hFakes_{sys}_{regionId[k]}")  # here all charges are in the same TH3
            hFakes_syst[sys][k].SetTitle(regionId[k])

            ## keep it for now, but would not use it. It also needs to project TH3 into a TH2 for plotting (not yet implemented below, this was coming from the older setup)
            ##
            if sys in plotFakesSyst.keys() and k != 3:  # no need to plot signal region here, it is not used
                outdirQCDsys = f"{outdirQCD}systematics/{sys}/"
                createPlotDirAndCopyPhp(outdirQCDsys)
                isTHn = False
                if "THn" in hFakes_syst[sys][k].ClassName():
                    isTHn = True
                for chbin in range(1,3):
                    charge = "plus" if chbin == 2 else "minus"
                    for zbin in plotFakesSyst[sys]:
                        plotHistName = f"hFakes_{sys}_{regionId[k]}_charge{charge}_{zbin}"
                        hfakeSys_tmp = hFakes_syst[sys][k].Clone(f"{plotHistName}_TMP")
                        if isTHn:
                            hfakeSys_tmp.GetAxis(2).SetRange(chbin, chbin)
                            hfakeSys_tmp.GetAxis(3).SetRange(zbin, zbin)
                            hfakeSys = hfakeSys_tmp.Projection(1, 0, "E") # project into a TH2
                            hfakeSys.SetName(plotHistName)
                        else:
                            hfakeSys = getTH2fromTH3(hfakeSys_tmp, plotHistName, chbin) # in this case the histo is a TH3 eta-pt-charge, since there is only one sytematic index, so the 4th axis was not really present
                        hfakeSys.SetTitle(f"{regionId[k]} charge {charge}")
                        drawCorrelationPlot(hfakeSys,
                                            hfakeSys.GetXaxis().GetTitle(),
                                            hfakeSys.GetYaxis().GetTitle(),
                                            f"Events ({sys} bin {zbin})",
                                            plotHistName,
                                            plotLabel="ForceTitle",
                                            outdir=outdirQCDsys,
                                            passCanvas=canvas,
                                            drawOption="COLZ0")

            
            ## do I really need this?
            ## check there are no negative bins
            if args.cropNegativeBin:
                wasCropped = cropNegativeContent(hFakes_syst[sys][k], silent=False, cropError=False)
                if wasCropped:
                    print(f"Histogram cropped to 0 in {regionId[k]} region and syst {sys}")

        templateQCDsys = hFakes_syst[sys][1].Clone(f"templateQCD_{sys}")
        templateQCDsys.Divide(hFakes_syst[sys][0])
        templateQCDsys.Multiply(hFakes_syst[sys][2])
        templateQCDsys.SetTitle(f"QCD template: {sys}")
        templateQCDsys.Write(f"{sys}__data_fakes")
        print(f"Writing histogram {sys}__data_fakes")

    print()
    print(">>> Now extracting the signal region for all other histograms with systematics, from the proper iso/mT bin")
    print()
    # get part of the histogram corresponding to signal region, which is bin number 4 of the isoMt axis 
    for sys in systNamesProcsHists.keys():
        for proc in systNamesProcsHists[sys].keys():

            systNamesProcsHists[sys][proc].GetAxis(nNomiDim-1).SetRange(4, 4) # select only signal region bin from isoMt axis (reminder: axis nNomiDim-1 is the isoMt axis, usually the 4th)
            if systNamesProcsHists[sys][proc].GetNdimensions() == nNomiDim:     
                # this is the case where TH4 goes into a TH3, in general we should not assume that, and use ProjectionND to return a THn, and then in case convert to a TH3 when the dimension is correct
                histRegion3 = systNamesProcsHists[sys][proc].Projection(0, 1, 2, "E")
            else:
                # we have the syst axis too here, so need to pick axes 0, 1, 2, 4, while for the isoMt axis (axis index 3) we need to project in the desired isoMt bin, which is the third (by chance, nothing to do with the axis number) 
                axesIds = array("i", [i for i in range(nNomiDim+1) if i != (nNomiDim-1)]) # isoMt is the last axis of the nominal histograms
                histRegion3 = systNamesProcsHists[sys][proc].ProjectionND(len(axesIds), axesIds, "E")

            histRegion3.SetName(f"{sys}__{proc}_region3")  # name is just temporary
            histRegion3.SetTitle(f"{sys}__{proc}")

            if args.cropNegativeBin:
                wasCropped = cropNegativeContent(histRegion3, silent=False, cropError=False)
                if wasCropped:
                    print(f"Histogram cropped to 0 in {regionId[3]} region for process {proc} and syst {sys}")

            histRegion3.Write(f"{sys}__{proc}")
            print(f"Writing histogram {sys}__{proc}")
            
    # now can close the file
    fout.Close()
    print(f"Shapes for all processes written in file {foutname}")
    print()


    if not args.makePlots:
        quit()  # part below still needs to be updated
    
    outdirSRcommon = outdir + "distributions_signalRegion/"
    createPlotDirAndCopyPhp(outdirSRcommon)

    canvas1D = ROOT.TCanvas("canvas1D", "", 800, 900)

    for chbin in range(1,3):
        charge = "plus" if chbin == 2 else "minus"

        outdirSR = outdirSRcommon + f"/{charge}/"
        createPlotDirAndCopyPhp(outdirSR)
        
        fShapesName = outdirSR + "plots.root"
        fShapes = ROOT.TFile.Open(fShapesName, "RECREATE")
        if not fShapes or fShapes.IsZombie():
            print(f"Error opening file {fShapesName}")
            quit()

        hdata2D = getTH2fromTH3(hdata, f"{hdata.GetName()}", chbin, chbin)
        # now we may do some more plotting of data and MC in signal region
        hmc2D = [getTH2fromTH3(h, f"{h.GetName()}", chbin, chbin) for h in hmc]
        hmc2D = sorted(hmc2D, key= lambda x: x.Integral()) # , reverse=True) 
        stack_eta = ROOT.THStack("stack_eta", "signal and backgrounds")
        stack_pt = ROOT.THStack("stack_pt", "signal and backgrounds")

        ratio2D = copy.deepcopy(hdata2D.Clone("dataOverMC2D"))
        den2D = copy.deepcopy(hdata2D.Clone("sigAndBkg2D"))
        den2D.Reset("ICESM")
        den2Dnofakes = copy.deepcopy(den2D.Clone("sigAndBkgNoFakes2D"))

        hdata2D.SetMarkerColor(ROOT.kBlack)
        hdata2D.SetLineColor(ROOT.kBlack)
        #hdata2D.SetLineWidth(2)
        hdata2D.SetMarkerStyle(20)
        hdata2D.SetMarkerSize(1)
        hdata2D.SetTitle("")

        # for projections along eta
        ptRange = ""
        if args.ptRangeProjection[0] < args.ptRangeProjection[1]:
            lowPtbin = max(1, hdata2D.GetYaxis().FindFixBin(args.ptRangeProjection[0]))
            highPtbin = min(hdata2D.GetNbinsY(), hdata2D.GetYaxis().FindFixBin(args.ptRangeProjection[1])) # hdata2D.GetNbinsY()
            ptRange = "_%gTo%g" % (hdata2D.GetYaxis().GetBinLowEdge(lowPtbin), hdata2D.GetYaxis().GetBinLowEdge(1+highPtbin))
            ptRange = ptRange.replace(".","p")
        else:
            lowPtbin = 0
            highPtbin = -1

        hdata_eta = hdata2D.ProjectionX("data_eta",lowPtbin,highPtbin,"e")
        hdata_pt  = hdata2D.ProjectionY("data_pt",0,-1,"e")

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

        hdata2D.Write("data2D")
        hdata_eta.Write()
        hdata_pt.Write()

        legend.AddEntry(hdata2D, "Data", "EP")
        for h in hmc2D:
            h.SetTitle("")
            h.SetFillColor(colors[h.GetName().replace('_vpt','')])
            h.SetLineColor(ROOT.kBlack)
            stack_eta.Add(h.ProjectionX(f"{h.GetName()}_eta",lowPtbin,highPtbin,"e"))
            stack_pt.Add( h.ProjectionY(f"{h.GetName()}_pt",0,-1,"e"))
            den2D.Add(h)
            if h.GetName() != "data_fakes":
                den2Dnofakes.Add(h)
            h.Write()
        for i in range(len(hmc2D)-1, 0, -1):
            legend.AddEntry(hmc2D[i], legEntries[hmc2D[i].GetName().replace('_vpt','')], "F")

        stack_eta.Write()
        stack_pt.Write()
        den2D.Write()
        den2Dnofakes.Write()

        ratio2D.Divide(den2D)
        ratio2D.Write()

        drawTH1dataMCstack(hdata_eta, stack_eta, "Muon #eta", "Events", "muon_eta_signalRegion" + ptRange,
                           outdirSR, legend, ratioPadYaxisNameTmp="Data/MC::0.92,1.08", passCanvas=canvas1D,
                           drawLumiLatex=True, xcmsText=0.3, noLegendRatio=True
        )
        drawTH1dataMCstack(hdata_pt, stack_pt, "Muon p_{T} (GeV)", "Events", "muon_pt_signalRegion",
                           outdirSR, legend, ratioPadYaxisNameTmp="Data/MC::0.92,1.08", passCanvas=canvas1D,
                           drawLumiLatex=True, xcmsText=0.3, noLegendRatio=True
        )

        ratio2D.SetTitle("data / (signal + background)")
        drawCorrelationPlot(ratio2D, xAxisName, yAxisName, "data/MC",
                            f"muon_eta_pt_dataMCratio", plotLabel="ForceTitle", outdir=outdirSR,
                            palette=57, passCanvas=canvas, drawOption="COLZ0", skipLumi=True)
        drawCorrelationPlot(ratio2D, xAxisName, yAxisName, "data/MC",
                            f"muon_eta_pt_dataMCratio_absUncertainty", plotLabel="ForceTitle", outdir=outdirSR,
                            palette=57, passCanvas=canvas, drawOption="COLZ0", skipLumi=True, plotError=True)


        allHists = hmc2D + [hdata2D]
        for h in allHists:
            h.SetTitle(h.GetName())
            drawCorrelationPlot(h, xAxisName, yAxisName, "Events",
                                f"muon_eta_pt_{h.GetName()}", plotLabel="ForceTitle", outdir=outdirSR,
                                palette=57, passCanvas=canvas, drawOption="COLZ0", skipLumi=True)
