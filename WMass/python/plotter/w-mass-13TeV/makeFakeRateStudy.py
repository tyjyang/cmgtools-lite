#!/bin/env python

import ROOT, os, sys, re, array, math
import time

sys.path.append(os.getcwd() + "/plotUtils/")
from utility import *


# python w-helicity-13TeV/makeFakeRateStudy.py plots/distribution/muonPlots_newSF/SKIMS_muons_latest/FRclosure_MarcFR_plus/test_plots.root plots/distribution/muonPlots_newSF/SKIMS_muons_latest/FRclosure_MarcFR_minus/test_plots.root -o plots/test_fakesChargeAsym/muon_MarcFR_coarseBin/ -n dataSubEWK_ratioPlusOverMinus -z "(data - EWK) template ratio" -t "charge + / charge -" --variable ptl1_large__etal1 -p "data,data_fakes,W,Z,TauTopVV" -r 0.7 1.3 --h1Dbinning "25,0.7,1.3" --palette -1 --yRange 26 54 --rebin 6 4



if __name__ == "__main__":
            
    from optparse import OptionParser
    parser = OptionParser(usage='%prog [options] filePlus fileMinus')
    parser.add_option('-o','--outdir',      dest='outdir',      default='', type='string', help='output directory to save things')
    parser.add_option('-n','--outhistname', dest='outhistname', default='', type='string', help='name for output canvas')
    parser.add_option('-x','--xAxisTitle',  dest='xAxisTitle',  default='', type='string', help='X axis title. If not given, use the one from hist1')
    parser.add_option('-y','--yAxisTitle',  dest='yAxisTitle',  default='', type='string', help='Y axis title. If not given, use the one from hist1')
    parser.add_option('-z','--zAxisTitle',  dest='zAxisTitle',  default='', type='string', help='Z axis title. If not given, use the one from hist1')
    parser.add_option('-t','--histTitle',   dest='histTitle',   default='', type='string', help='Title to assign to output histogram. It is used as a label for the canvas')
    parser.add_option('-r','--ratioRange',  dest='ratioRange',  default=(0, 2),type="float", nargs=2, help="Min and max for the ratio in the plot")
    parser.add_option(     '--h1Dbinning',  dest='h1Dbinning',  default='50,0.9,1.1', type='string', help='Comma separated list of 3 numbers: nbins,min,max')
    parser.add_option('-v','--valBadRatio', dest='valBadRatio', default='0', type='float', help='Value to be used in case of bad ratio (division by 0). The 1D histogram is not filled in case of bad ratio')
    parser.add_option(     '--xRange'     , dest='xRange', default=(0,-1), type='float', nargs=2, help='Select range for X axis to plot. Also, bins outside this range are not considered in the 1D histogram. If min > max, the option is neglected')
    parser.add_option(     '--yRange'     , dest='yRange', default=(0,-1), type='float', nargs=2, help='Select range for Y axis to plot. Also, bins outside this range are not considered in the 1D histogram. If min > max, the option is neglected')
    parser.add_option(     '--palette'  , dest='palette',      default=55, type=int, help='Set palette: use a negative number to select a built-in one, otherwise the default is 55 (kRainbow)')
    parser.add_option(     '--variable',    dest='variable',    default='ptl1large__etal1', type='string', help='Variable to get histogram inside file')
    parser.add_option('-p','--processes',   dest='processes',   default='data,data_fakes,W,Z,TauTopVVFlips', type='string', help='Comma separated list of processes (to build histogram\'s name to be taken from input file)')
    parser.add_option(     '--rebin'     , dest='rebin', default=(-1,-1), type='int', nargs=2, help='Rebin bins along x and y axis by these factors. Keep one of the two ngative to rebin only the others')
    parser.add_option(      '--no-lumi-xsec',   dest='noLumiXsec',   default=False, action='store_true', help='Do not propagate luminosity and cross section uncertainty')
    (options, args) = parser.parse_args()

    if len(sys.argv) < 2:
        parser.print_usage()
        quit()

    f1 = args[0]
    f2 = args[1]
    ROOT.TH1.SetDefaultSumw2()
    ROOT.gStyle.SetPaintTextFormat('.3g')

    print ""

    if options.outdir:
        outname = options.outdir
        addStringToEnd(outname,"/",notAddIfEndswithMatch=True)
        createPlotDirAndCopyPhp(outname)
    else:
        print "Error: you should specify an output folder using option -o <name>. Exit"
        quit()
    if not options.outhistname:
        print "Error: you should specify an output histogram name using option -n <name>. "
        print "Exit"
        quit()


    procNames = options.processes.split(',')    
    fakesName = "MISSING FAKES"
    for name in procNames:
        if "data_fakes" in name: fakesName = name
    #fakesName = "data_fakes" if "data_fakes" in procNames else "data_fakes_test" if "data_fakes_test" in procNames else "MISSING_FAKES"
    if not all(x in procNames for x in ["data", fakesName]):
        print "Warning: This script requires processes named data and data_fakes. Exit"
        quit()
    hratio = 0

    procsPlus = {}
    procsMinus = {}
    procsPlus_noLumiXsec = {}
    procsMinus_noLumiXsec = {}

    # read histograms
    ######################
    # file +
    tf = ROOT.TFile.Open(f1)        
    for p in procNames:
        hname = options.variable + "_" + p
        procsPlus[p] = tf.Get(hname)
        if (procsPlus[p] == 0):
            print "Error: could not retrieve %s from input file %s. Exit" % (hname,f1)
            quit()
        else:
            procsPlus[p].SetDirectory(0)
    tf.Close()

    # ----------------------------------
    # file -
    tf = ROOT.TFile.Open(f2)        
    for p in procNames:
        hname = options.variable + "_" + p
        procsMinus[p] = tf.Get(hname)
        if (procsMinus[p] == 0):
            print "Error: could not retrieve %s from input file %s. Exit" % (hname,f2)
            quit()
        else:
            procsMinus[p].SetDirectory(0)
    tf.Close()
    ######################

    xsecUncProc = {"Wincl" : 0.038,
                   "Z" : 0.04,
                   "TauTopVVFlips" : 0.038,                   
                   "W" : 0.038,
                   "TauTopVV" : 0.038
                   }  # relative uncertainty

    (rebinx,rebiny) = options.rebin
    for p in procNames:
        if rebiny > 0:
            procsPlus[p].RebinY(rebiny)
            procsMinus[p].RebinY(rebiny)
        if rebinx > 0:
            procsPlus[p].RebinX(rebinx)
            procsMinus[p].RebinX(rebinx)

    for p in procNames:
        procsPlus_noLumiXsec[p]  = procsPlus[p].Clone(procsPlus[p].GetName()+"_noLumiXsec")
        procsMinus_noLumiXsec[p] = procsMinus[p].Clone(procsMinus[p].GetName()+"_noLumiXsec")

    # sum uncertainty for xsec (above) and 2.5% of luminosity (there would be a 1% of efficiency but it is negligible)
    # loop on one histogram, the binning is the same
    for ix in range(1,1+procsPlus["data"].GetNbinsX()):
        for iy in range(1,1+procsPlus["data"].GetNbinsY()):
            for p in procNames:
                if any(x == p for x in ["data", fakesName]): continue
                # charge plus
                newerr = procsPlus[p].GetBinError(ix, iy)
                bincontent = procsPlus[p].GetBinContent(ix, iy)
                newerr = pow(newerr,2) 
                if not options.noLumiXsec: newerr +=  pow(xsecUncProc[p]*bincontent,2) + pow(0.025*bincontent,2)
                procsPlus[p].SetBinError(ix, iy, math.sqrt(newerr))
                # charge minus
                newerr = procsMinus[p].GetBinError(ix, iy)
                bincontent = procsMinus[p].GetBinContent(ix, iy)
                newerr = pow(newerr,2)
                if not options.noLumiXsec: newerr += pow(xsecUncProc[p]*bincontent,2) + pow(0.025*bincontent,2)
                procsMinus[p].SetBinError(ix, iy, math.sqrt(newerr))
        
    data_subEWK_plus  = procsPlus["data"].Clone("data_subEWK_plus")
    data_subEWK_minus = procsMinus["data"].Clone("data_subEWK_minus")
    data_subEWK_plus_noLumiXsec  = procsPlus_noLumiXsec["data"].Clone("data_subEWK_plus_noLumiXsec")
    data_subEWK_minus_noLumiXsec = procsMinus_noLumiXsec["data"].Clone("data_subEWK_minus_noLumiXsec")

    for p in procNames:
        if any(x == p for x in ["data", fakesName]): continue
        data_subEWK_plus.Add(procsPlus[p], -1.)        
        data_subEWK_minus.Add(procsMinus[p], -1.)
        data_subEWK_plus_noLumiXsec.Add(procsPlus_noLumiXsec[p], -1.)        
        data_subEWK_minus_noLumiXsec.Add(procsMinus_noLumiXsec[p], -1.)

    xMin = options.xRange[0]
    xMax = options.xRange[1]
    yMin = options.yRange[0]
    yMax = options.yRange[1]

    # this ratio will only propagate statistical uncertainty on inputs
    #hratio = data_subEWK_plus_noLumiXsec.Clone("data_subEWK_PlusOverMinus")
    hratio = data_subEWK_plus.Clone("data_subEWK_PlusOverMinus")

    nbins,minx,maxx = options.h1Dbinning.split(',')
    hratioDistr = ROOT.TH1D(options.outhistname+"_1D","Distribution of ratio values",int(nbins),float(minx),float(maxx))

    nout = 0    
    for ix in range(1,1+hratio.GetNbinsX()):
        for iy in range(1,1+hratio.GetNbinsY()):
            xval = hratio.GetXaxis().GetBinCenter(ix) 
            yval = hratio.GetYaxis().GetBinCenter(iy) 
            hist2xbin = data_subEWK_minus.GetXaxis().FindFixBin(xval)
            hist2ybin = data_subEWK_minus.GetYaxis().FindFixBin(yval)
            if xMin < xMax:
                if xval < xMin or xval > xMax: continue
            if yMin < yMax:
                if yval < yMin or yval > yMax: continue            
            denval = data_subEWK_minus.GetBinContent(hist2xbin, hist2ybin)
            if denval != 0:                
                numval = hratio.GetBinContent(ix,iy)
                ratio = numval / denval
                #hratioDistr.Fill(ratio)
                #hratio.SetBinContent(ix,iy,ratio)
                if ratio < float(minx) or ratio > float(maxx): nout += 1
            else: 
                print "Warning: found division by 0 in one bin: setting ratio to " + str(options.valBadRatio)
                #hratio.SetBinContent(ix,iy,options.valBadRatio)

    hratio.Divide(data_subEWK_minus)
    #hratio.Divide(data_subEWK_minus_noLumiXsec)

    print "nout = " + str(nout)

    hratio.SetTitle(options.histTitle)

    if options.xAxisTitle: hratio.GetXaxis().SetTitle(options.xAxisTitle)    
    if options.yAxisTitle: hratio.GetYaxis().SetTitle(options.yAxisTitle)    
    if options.zAxisTitle: hratio.GetZaxis().SetTitle(options.zAxisTitle)    
    xAxisTitle = hratio.GetXaxis().GetTitle()
    yAxisTitle = hratio.GetYaxis().GetTitle()
    zAxisTitle = hratio.GetZaxis().GetTitle()

    if xMax > xMin and not "::" in xAxisTitle: 
        xAxisTitle = xAxisTitle + "::" + str(options.xRange[0]) + "," + str(options.xRange[1])
    if yMax > yMin and not "::" in yAxisTitle: 
        yAxisTitle = yAxisTitle + "::" + str(options.yRange[0]) + "," + str(options.yRange[1])

    adjustSettings_CMS_lumi()

    canvas2D = ROOT.TCanvas("canvas2D","",900,700)
    
    # the axis name can be used to set the range if it is in the format "name::min,maz"
    # if this is not already the case, use the selected range from the input option
    if not "::" in zAxisTitle:  
        zAxisTitle = zAxisTitle + "::" + str(options.ratioRange[0]) + "," + str(options.ratioRange[1])

    drawCorrelationPlot(hratio,xAxisTitle,yAxisTitle,zAxisTitle,
                        options.outhistname,"ForceTitle",outname,0,0,False,False,False,1,palette=options.palette,passCanvas=canvas2D, 
                        drawOption="colz0 text60 E")

    canvas = ROOT.TCanvas("canvas","",800,700)    

    dataSubEWK_over_fakes_plus = data_subEWK_plus.Clone("dataSubEWK_over_fakes_plus")
    dataSubEWK_over_fakes_plus.Divide(procsPlus[fakesName])
    dataSubEWK_over_fakes_minus = data_subEWK_minus.Clone("dataSubEWK_over_fakes_minus")
    dataSubEWK_over_fakes_minus.Divide(procsMinus[fakesName])
    dataSubEWK_over_fakes_plus.SetTitle("Charge +")
    dataSubEWK_over_fakes_minus.SetTitle("Charge -")

    dataSubEWK_over_fakes_plus_noLumiXsec = data_subEWK_plus_noLumiXsec.Clone("dataSubEWK_over_fakes_plus_noLumiXsec")
    dataSubEWK_over_fakes_plus_noLumiXsec.Divide(procsPlus_noLumiXsec[fakesName])
    dataSubEWK_over_fakes_minus_noLumiXsec = data_subEWK_minus_noLumiXsec.Clone("dataSubEWK_over_fakes_minus_noLumiXsec")
    dataSubEWK_over_fakes_minus_noLumiXsec.Divide(procsMinus_noLumiXsec[fakesName])


    zrange = "0.7,1.3"

    zAxisTitle = "(data - EWK) / fakes::" + zrange
    drawCorrelationPlot(dataSubEWK_over_fakes_plus,
                        xAxisTitle,yAxisTitle,zAxisTitle,
                        "dataSubEWK_over_fakes_plus","ForceTitle",outname,0,0,False,False,False,1,palette=options.palette,passCanvas=canvas2D,
                        drawOption="colz0 text60 E")

    zAxisTitle = "(data - EWK) / fakes::" + zrange
    drawCorrelationPlot(dataSubEWK_over_fakes_minus,
                        xAxisTitle,yAxisTitle,zAxisTitle,
                        "dataSubEWK_over_fakes_minus","ForceTitle",outname,0,0,False,False,False,1,palette=options.palette,passCanvas=canvas2D,
                        drawOption="colz0 text60 E")

    dataSubEWKoverFakes_ratioPlusOverMinus = dataSubEWK_over_fakes_plus_noLumiXsec.Clone("dataSubEWKoverFakes_ratioPlusOverMinus_noLumiXsec")
    dataSubEWKoverFakes_ratioPlusOverMinus.Divide(dataSubEWK_over_fakes_minus_noLumiXsec)
    #zrange = "0.7,1.3"
    dataSubEWKoverFakes_ratioPlusOverMinus.SetTitle("charge + / charge -")
    zAxisTitle = "ratio of (data-EWK)/fakes between charges::" + zrange
    drawCorrelationPlot(dataSubEWKoverFakes_ratioPlusOverMinus,
                        xAxisTitle,yAxisTitle,zAxisTitle,
                        "dataSubEWKoverFakes_ratioPlusOverMinus","ForceTitle",outname,0,0,False,False,False,1,palette=options.palette,passCanvas=canvas2D,
                        drawOption="colz0 text60 E")
    

    for i in range(1,dataSubEWKoverFakes_ratioPlusOverMinus.GetNbinsX()+1):
        for j in range(1,dataSubEWKoverFakes_ratioPlusOverMinus.GetNbinsY()+1):
            if xMin < xMax:
                xval = dataSubEWKoverFakes_ratioPlusOverMinus.GetXaxis().GetBinCenter(i)
                if xval < xMin or xval > xMax: continue
            if yMin < yMax:
                yval = dataSubEWKoverFakes_ratioPlusOverMinus.GetYaxis().GetBinCenter(j)
                if yval < yMin or yval > yMax: continue            

            hratioDistr.Fill(dataSubEWKoverFakes_ratioPlusOverMinus.GetBinContent(i,j))

    drawTH1(hratioDistr, 
            dataSubEWKoverFakes_ratioPlusOverMinus.GetZaxis().GetTitle(),
            "number of events",
            outname,
            "summary",
            "dataSubEWKoverFakes_ratioPlusOverMinus",
            passCanvas=canvas
            )


    fakesRatio = procsPlus[fakesName].Clone("fakesRatio")
    fakesRatio.Divide(procsMinus[fakesName])
    fakesRatio.SetTitle("charge + / charge -")
    zAxisTitle = "fakes template ratio::" + str(options.ratioRange[0]) + "," + str(options.ratioRange[1])
    drawCorrelationPlot(fakesRatio,
                        xAxisTitle,yAxisTitle,zAxisTitle,
                        "fakes_PlusOverMinus","ForceTitle",outname,0,0,False,False,False,1,palette=options.palette,passCanvas=canvas2D,
                        drawOption="colz0 text60 E")


