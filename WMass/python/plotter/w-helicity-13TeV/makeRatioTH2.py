#!/bin/env python

# script to make ratio of two TH2 histograms
# the two objects can have different bin ranges, but in order for the ratio to make sense, 
# the binning of the smaller TH2 should be a subset of the other histogram's binning
# In any case, the ratio is made looping on the first histogram passed (which is cloned to get a new one with the same binning) and getting the bin center
# for each bin, then that value is used to look for the bin of the other histogram
# therefore, in case the two TH2 have different binning, make sure that the first one has the finer granularity

###################
# >>>>>>>>> WARNING <<<<<<<<<<<< 
###################
# Apparently it works in release CMSSW_10_2_0_pre4
# while in CMSSW_8_0_25 it produces a segmentatio fault from the setTDRStyle()
# Maybe it is just due to the root version

################################
# Exaples
################################

# python w-helicity-13TeV/makeRatioTH2.py ~/www/wmass/13TeV/scaleFactors/electron/fullID_extPt/smoothEfficiency_electrons_fullID.root scaleFactor ~/www/wmass/13TeV/scaleFactors/electron/fullID_noErfPlusLine/smoothEfficiency_electrons_fullID.root scaleFactor -o ~/www/wmass/13TeV/scaleFactors/ratio/electron/fullID_extPt__over__fullID/ -f ratio.root -n ratio2D -t "full ID scale factor ratio" -z ratio --ratioRange 0.9 1.1

# python w-helicity-13TeV/makeRatioTH2.py ../../data/fakerate/fakeRateSmoothed_el_mT40_35p9fb_signedEta_pt65_fullWMC_newTrigSF_fitpol2.root fr_pt_eta_ewk ~/www/wmass/13TeV/fake-rate/electron/FR_graphs/fakeRate_eta_pt_granular_mT40_35p9fb_signedEta_pt65_subtrAllMC_newTrigSF_fitpol2/fakeRateSmoothed_el_mT40_35p9fb_signedEta_pt65_subtrAlllMC_newTrigSF_fitpol2.root fr_pt_eta_ewk -o ~/www/wmass/13TeV/fake-rate/electron/FR_graphs/fakeRate_eta_pt_granular_mT40_35p9fb_signedEta_pt65_subtrAllMC_newTrigSF_fitpol2/ -f ratio_PR_WZ_allMC.root -n ratioPR_WZ_allMC -t "PR W,Z MC / all MC" -z "PR ratio" -r 0.98 1.02

################################
################################


import ROOT, os, sys, re, array, math

sys.path.append(os.getcwd() + "/plotUtils/")
from utility import *

ROOT.gROOT.SetBatch(True)
        
if __name__ == "__main__":
            
    from optparse import OptionParser
    parser = OptionParser(usage='%prog [options] file1 hist1 file2 hist2')
    parser.add_option('-o','--outdir',      dest='outdir',      default='', type='string', help='output directory to save things')
    parser.add_option('-f','--outfilename', dest='outfilename', default='', type='string', help='Name of output file to save results')
    parser.add_option('-n','--outhistname', dest='outhistname', default='', type='string', help='Name of output histogram saved in output file')
    parser.add_option('-x','--xAxisTitle',  dest='xAxisTitle',  default='', type='string', help='X axis title. If not given, use the one from hist1')
    parser.add_option('-y','--yAxisTitle',  dest='yAxisTitle',  default='', type='string', help='Y axis title. If not given, use the one from hist1')
    parser.add_option('-z','--zAxisTitle',  dest='zAxisTitle',  default='', type='string', help='Z axis title. If not given, use the one from hist1')
    parser.add_option('-t','--histTitle',   dest='histTitle',   default='', type='string', help='Title to assign to output histogram. It is used as a label for the canvas')
    parser.add_option('-r','--ratioRange',  dest='ratioRange',  default=(0, 2),type="float", nargs=2, help="Min and max for the ratio in the plot")
    parser.add_option(     '--h1Dbinning',  dest='h1Dbinning',  default='50,0.9,1.1', type='string', help='Comma separated list of 3 numbers: nbins,min,max')
    parser.add_option('-v','--valBadRatio', dest='valBadRatio', default='0', type='float', help='Value to be used in case of bad ratio (division by 0). The 1D histogram is not filled in case of bad ratio')
    (options, args) = parser.parse_args()

    if len(sys.argv) < 4:
        parser.print_usage()
        quit()

    f1 = args[0]
    h1 = args[1]
    f2 = args[2]
    h2 = args[3]

    ROOT.TH1.SetDefaultSumw2()

    print ""

    if options.outdir:
        outname = options.outdir
        addStringToEnd(outname,"/",notAddIfEndswithMatch=True)
        createPlotDirAndCopyPhp(outname)
    else:
        print "Error: you should specify an output folder using option -o <name>. Exit"
        quit()
    if not options.outfilename:
        print "Error: you should specify an output file name using option -f <name>. Exit"
        quit()
    if not options.outhistname:
        print "Error: you should specify an output histogram name using option -n <name>. Exit"
        quit()


    hratio = 0

    # file 1
    tf = ROOT.TFile.Open(f1)        
    hist1 =   tf.Get(h1)
    if (hist1 == 0):
        print "Error: could not retrieve %s from input file %s. Exit" % (h1,f1)
        quit()
    else:
        hist1.SetDirectory(0)
    tf.Close()

    # file2
    tf = ROOT.TFile.Open(f2)        
    hist2 =   tf.Get(h2)
    if (hist2 == 0):
        print "Error: could not retrieve %s from input file %s. Exit" % (h2,f2)
        quit()
    else:
        hist2.SetDirectory(0)
    tf.Close()

    hratio = hist1.Clone(options.outhistname)
    hratio.SetTitle(options.histTitle)

    nbins,minx,maxx = options.h1Dbinning.split(',')
    hratioDistr = ROOT.TH1D(options.outhistname+"_1D","Distribution of ratio values",int(nbins),float(minx),float(maxx))

    for ix in range(1,1+hratio.GetNbinsX()):
        for iy in range(1,1+hratio.GetNbinsY()):
            xval = hratio.GetXaxis().GetBinCenter(ix) 
            yval = hratio.GetYaxis().GetBinCenter(iy) 
            hist2xbin = hist2.GetXaxis().FindFixBin(xval)
            hist2ybin = hist2.GetYaxis().FindFixBin(yval)
            if hist2.GetBinContent(hist2xbin, hist2ybin) != 0:
                ratio = hratio.GetBinContent(ix,iy) / hist2.GetBinContent(hist2xbin, hist2ybin)
                hratioDistr.Fill(ratio)
                hratio.SetBinContent(ix,iy,ratio)
            else: 
                hratio.SetBinContent(ix,iy,options.valBadRatio)

    if options.xAxisTitle: hratio.GetXaxis().SetTitle(options.xAxisTitle)    
    if options.yAxisTitle: hratio.GetYaxis().SetTitle(options.yAxisTitle)    
    if options.zAxisTitle: hratio.GetZaxis().SetTitle(options.zAxisTitle)    
    xAxisTitle = hratio.GetXaxis().GetTitle()
    yAxisTitle = hratio.GetYaxis().GetTitle()
    zAxisTitle = hratio.GetZaxis().GetTitle()

    # print "xAxisTitle = " + xAxisTitle
    # print "yAxisTitle = " + yAxisTitle
    # print "zAxisTitle = " + zAxisTitle

    adjustSettings_CMS_lumi()
    # the axis name can be used to set the range if it is in the format "name::min,maz"
    # if this is not already the case, use the selected range from the input option
    if not "::" in zAxisTitle:  
        zAxisTitle = zAxisTitle + "::" + str(options.ratioRange[0]) + "," + str(options.ratioRange[1])
    drawCorrelationPlot(hratio,xAxisTitle,yAxisTitle,zAxisTitle,
                        options.outhistname,"ForceTitle",outname,0,0,False,False,False,1,palette=55)
    
    canvas = ROOT.TCanvas("canvas","",800,700)
    canvas.SetTickx(1)
    canvas.SetTicky(1)
    canvas.cd()
    canvas.SetLeftMargin(0.12)
    canvas.SetRightMargin(0.04)
    canvas.cd()

    hratioDistr.SetLineColor(ROOT.kBlack)
    hratioDistr.SetLineWidth(2)
    hratioDistr.GetXaxis().SetTitle(hratio.GetZaxis().GetTitle() if options.zAxisTitle else "ratio")
    hratioDistr.GetXaxis().SetTitleOffset(1.2)
    hratioDistr.GetXaxis().SetTitleSize(0.05)
    hratioDistr.GetXaxis().SetLabelSize(0.04)
    hratioDistr.GetYaxis().SetTitle("number of events")
    hratioDistr.GetYaxis().SetTitleOffset(1.15)
    hratioDistr.GetYaxis().SetTitleSize(0.05)
    hratioDistr.GetYaxis().SetLabelSize(0.04)
    hratioDistr.Draw("HIST")
    canvas.RedrawAxis("sameaxis")
    setTDRStyle()
    # force drawing stat box
    ROOT.gPad.Update()
    ROOT.gStyle.SetOptStat(1110)
    ROOT.gStyle.SetOptFit(1102)
    #
    for ext in ["png","pdf"]:
        canvas.SaveAs(outname + "ratioDistribution_{hname}.{ext}".format(hname=options.outhistname,ext=ext))
 
    ###########################
    # Now save things
    ###########################
    tf = ROOT.TFile.Open(outname+options.outfilename,'recreate')
    hratio.Write(options.outhistname)
    tf.Close()
    print ""
    print "Created file %s" % (outname+options.outfilename)
    print ""

                               
         
