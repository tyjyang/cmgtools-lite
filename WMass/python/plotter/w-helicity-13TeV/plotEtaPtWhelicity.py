import ROOT, os, datetime, re, operator, math
from array import array
ROOT.gROOT.SetBatch(True)

import sys
sys.path.append(os.getcwd() + "/plotUtils/")
from utility import *

if __name__ == "__main__":

    from optparse import OptionParser
    parser = OptionParser(usage="%prog inputfile [options]")
    parser.add_option('-o','--outdir', dest='outdir', default='.', type='string', help='output directory to save things')
    (options, args) = parser.parse_args()

    ROOT.TH1.SetDefaultSumw2()

    if len(args) < 1:
        parser.print_usage()
        quit()

    if not options.outdir:
        print "Error: please specify an output folder with option -o (a subfolder named as options.charge will be added inside)"
        quit()

    outname = options.outdir
    addStringToEnd(outname,"/",notAddIfEndswithMatch=True)
    createPlotDirAndCopyPhp(outname)

    Whists = {}
    infile = ROOT.TFile(args[0], 'read')
    if not infile:
        print "Error: I couldn't open file " + args[0]
        quit()
    for charge in ["plus", "minus"]:
        for pol in ["left", "right", "long"]:
            key = "ptl1__etal1_W{ch}_{pol}".format(ch=charge,pol=pol)
            print key
            hdummy = infile.Get(key)
            if not hdummy:
                print "Error, I couldn't get histogram " + key
                quit()
            Whists[key] = hdummy.Clone(key)
            Whists[key].SetDirectory(0)
    # wLp = infile.Get("ptl1__etal1_Wplus_left")
    # wRp = infile.Get("ptl1__etal1_Wplus_right")
    # w0p = infile.Get("ptl1__etal1_Wplus_long")
    # wLm = infile.Get("ptl1__etal1_Wminus_left")
    # wRm = infile.Get("ptl1__etal1_Wminus_right")
    # w0m = infile.Get("ptl1__etal1_Wminus_long")
    # wLp.SetDirectory(0)
    # wRp.SetDirectory(0)
    # w0p.SetDirectory(0)
    # wLm.SetDirectory(0)
    # wRm.SetDirectory(0)
    # w0m.SetDirectory(0)
    infile.Close()

    zmin =  99999.9
    zmax = -99999.9
    for key in Whists.keys():
        zminval,zmaxval = getMinMaxHisto(Whists[key],sumError=False)
        if zminval < zmin: zmin = zminval
        if zmaxval > zmax: zmax = zmaxval

    xaxisTitle = "reconstructed muon p_{T} [GeV]"
    yaxisTitle = "reconstructed muon #eta"
    zaxisTitle = "Events"

    adjustSettings_CMS_lumi()
    canvas = ROOT.TCanvas("canvas","",1000,800)

    for key in Whists.keys():
        Whists[key].SetTitle("")
        for izrange in ["sameZrange", "histZrange"]:
            
            real_outname = outname + izrange + "/"
            createPlotDirAndCopyPhp(real_outname)
            
            if izrange == "sameZrange": 
                real_zaxisTitle = zaxisTitle + "::{:.1f},{:.1f}".format(0.995*zmin,1.005*zmax)
            else:
                zminval,zmaxval = getMinMaxHisto(Whists[key],sumError=False)
                real_zaxisTitle = zaxisTitle + "::{:.1f},{:.1f}".format(0.995*zminval,1.005*zmaxval)
                
            canvas = drawCorrelationPlot(Whists[key], 
                                         xaxisTitle, yaxisTitle, real_zaxisTitle, 
                                         key, "ForceTitle", real_outname,1,1,False,False,False,1,
                                         passCanvas=canvas,palette=57, lumi="",
                                         writeSimulation=True)

            canvas.cd()
            text = ROOT.TLatex()
            text.SetNDC()
            text.SetTextSize(0.04)  # 0.03
            text.SetTextFont(42)
            text.DrawLatex(0.16, 0.95, "#bf{CMS} #it{Simulation Preliminary}")
            text.DrawLatex(0.63, 0.95, "35.9 fb^{-1} (13 TeV)")
            outname_test = real_outname # + "test/"
            createPlotDirAndCopyPhp(outname_test)
            canvas.RedrawAxis("sameaxis")
            for ext in [".pdf", ".png"]:
                canvas.SaveAs(outname_test + key + ext)
