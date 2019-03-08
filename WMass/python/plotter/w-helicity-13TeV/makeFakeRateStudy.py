#!/bin/env python

import ROOT, os, sys, re, array, math
import time

sys.path.append(os.getcwd() + "/plotUtils/")
from utility import *

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
    parser.add_option('-p','--processes',   dest='processes',   default='data,data_fakes,Wincl,Z,TauTopVVFlips', type='string', help='Comma separated list of processes (to build histogram\'s name to be taken from input file)')
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
    hratio = 0

    procsPlus = {}
    procsMinus = {}

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
                   "TauTopVVFlips" : 0.038                   
                   }  # relative uncertainty

    for p in procNames:
        procsPlus[p].RebinY(3)
        procsPlus[p].RebinX(5)
        procsMinus[p].RebinY(3)
        procsMinus[p].RebinX(5)

    # patch
    #procsPlus["data_fakes"].Scale(1./1.1377)
    #procsMinus["data_fakes"].Scale(1./1.143)

    # sum uncertainty for xsec (above) and 2.5% of luminosity (there would be a 1% of efficiency but it is negligible)
    # loop on one histogram, the binning is the same
    for ix in range(1,1+procsPlus["data"].GetNbinsX()):
        for iy in range(1,1+procsPlus["data"].GetNbinsY()):
            for p in procNames:
                if any(x == p for x in ["data", "data_fakes"]): continue
                # charge plus
                newerr = procsPlus[p].GetBinError(ix, iy)
                bincontent = procsPlus[p].GetBinContent(ix, iy)
                newerr = pow(newerr,2) + pow(xsecUncProc[p]*bincontent,2) + pow(0.025*bincontent,2)
                procsPlus[p].SetBinError(ix, iy, math.sqrt(newerr))
                # charge minus
                newerr = procsMinus[p].GetBinError(ix, iy)
                bincontent = procsMinus[p].GetBinContent(ix, iy)
                newerr = pow(newerr,2) + pow(xsecUncProc[p]*bincontent,2) + pow(0.025*bincontent,2)
                procsMinus[p].SetBinError(ix, iy, math.sqrt(newerr))
        
    data_subEWK_plus  = procsPlus["data"].Clone("data_subEWK_plus")
    data_subEWK_minus = procsMinus["data"].Clone("data_subEWK_minus")
    for p in procNames:
        if any(x == p for x in ["data", "data_fakes"]): continue
        data_subEWK_plus.Add(procsPlus[p], -1.)        
        data_subEWK_minus.Add(procsMinus[p], -1.)

    xMin = options.xRange[0]
    xMax = options.xRange[1]
    yMin = options.yRange[0]
    yMax = options.yRange[1]

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
                hratioDistr.Fill(ratio)
                #hratio.SetBinContent(ix,iy,ratio)
                if ratio < float(minx) or ratio > float(maxx): nout += 1
            else: 
                print "Warning: found division by 0 in one bin: setting ratio to " + str(options.valBadRatio)
                #hratio.SetBinContent(ix,iy,options.valBadRatio)

    hratio.Divide(data_subEWK_minus)

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
    drawTH1(hratioDistr, 
            hratio.GetZaxis().GetTitle() if options.zAxisTitle else "ratio",
            "number of events",
            outname,
            "ratioDistribution",
            options.outhistname,
            passCanvas=canvas
            )
    

    dataSubEWK_over_fakes_plus = data_subEWK_plus.Clone("dataSubEWK_over_fakes_plus")
    dataSubEWK_over_fakes_plus.Divide(procsPlus["data_fakes"])
    dataSubEWK_over_fakes_minus = data_subEWK_minus.Clone("dataSubEWK_over_fakes_minus")
    dataSubEWK_over_fakes_minus.Divide(procsMinus["data_fakes"])
    dataSubEWK_over_fakes_plus.SetTitle("Charge +")
    dataSubEWK_over_fakes_minus.SetTitle("Charge -")

    zrange = "0.7,1.3"

    zAxisTitle = "(data - EWK) / QCD prediction::" + zrange
    drawCorrelationPlot(dataSubEWK_over_fakes_plus,
                        xAxisTitle,yAxisTitle,zAxisTitle,
                        "dataSubEWK_over_fakes_plus","ForceTitle",outname,0,0,False,False,False,1,palette=options.palette,passCanvas=canvas2D,
                        drawOption="colz0 text60 E")

    zAxisTitle = "(data - EWK) / QCD prediction::" + zrange
    drawCorrelationPlot(dataSubEWK_over_fakes_minus,
                        xAxisTitle,yAxisTitle,zAxisTitle,
                        "dataSubEWK_over_fakes_minus","ForceTitle",outname,0,0,False,False,False,1,palette=options.palette,passCanvas=canvas2D,
                        drawOption="colz0 text60 E")

    dataSubEWKoverFakes_ratioPlusOverMinus = dataSubEWK_over_fakes_plus.Clone("dataSubEWKoverFakes_ratioPlusOverMinus")
    dataSubEWKoverFakes_ratioPlusOverMinus.Divide(dataSubEWK_over_fakes_minus)
    zrange = "0.9,1.1"
    dataSubEWKoverFakes_ratioPlusOverMinus.SetTitle("charge + / charge -")
    zAxisTitle = "ratio of (data-EWK)/fakes between charges::" + zrange
    drawCorrelationPlot(dataSubEWKoverFakes_ratioPlusOverMinus,
                        xAxisTitle,yAxisTitle,zAxisTitle,
                        "dataSubEWKoverFakes_ratioPlusOverMinus","ForceTitle",outname,0,0,False,False,False,1,palette=options.palette,passCanvas=canvas2D,
                        drawOption="colz0 text60 E")
    

    fakesRatio = procsPlus["data_fakes"].Clone("fakesRatio")
    fakesRatio.Divide(procsMinus["data_fakes"])
    fakesRatio.SetTitle("charge + / charge -")
    zAxisTitle = "QCD prediction template ratio::" + str(options.ratioRange[0]) + "," + str(options.ratioRange[1])
    drawCorrelationPlot(fakesRatio,
                        xAxisTitle,yAxisTitle,zAxisTitle,
                        "fakes_PlusOverMinus","ForceTitle",outname,0,0,False,False,False,1,palette=options.palette,passCanvas=canvas2D,
                        drawOption="colz0 text60 E")
