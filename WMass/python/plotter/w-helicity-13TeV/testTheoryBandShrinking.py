#!/usr/bin/env python                                                       

# script to test theory band in plots: when doing the normalized cross section, the band shrink when the crossing point is located (crossover effect)

#from shutil import copyfile
import re, sys, os, os.path, ROOT, copy, math
from array import array
#import numpy as np
import root_numpy

ROOT.gROOT.SetBatch(True)

import utilities
utilities = utilities.util()

from make_diff_xsec_cards import getDiffXsecBinning
from make_diff_xsec_cards import templateBinning

sys.path.append(os.getcwd() + "/plotUtils/")
from utility import *


if __name__ == "__main__":

    from optparse import OptionParser
    parser = OptionParser(usage='%prog [options]')
    parser.add_option('-c','--charge', dest='charge', default='plus,minus', type='string', help='charges to run')
    parser.add_option('-b','--bin', dest='testbin', default=13, type='int', help='Bin number to test (vs eta)')
    parser.add_option('-o','--outdir', dest='outdir', default='', type='string', help='output directory to save things')
    (options, args) = parser.parse_args()

    ROOT.TH1.SetDefaultSumw2()
    ROOT.TH1.StatOverflows(True)

    charges = [x for x in options.charge.split(',')]
    for c in charges:
        if c not in ["plus", "minus"]:
            print "Error: unknown charge %s (select 'plus' or 'minus' or both separated by comma)" % c
            quit()
    hasBothCharges = True
    if len(charges) < 2:
        hasBothCharges = False

    if options.outdir:
        outdir = options.outdir
        addStringToEnd(outdir,"/",notAddIfEndswithMatch=True)
        outdir = outdir + "/"
        createPlotDirAndCopyPhp(outdir)
    else:
        print "Error: you should specify an output folder using option -o <name>. Exit"
        quit()


    print "Now retrieving all the theory variations, to make the theory bands later on"
    histDict = utilities.getTheoryHistDiffXsecFast(xsecWithWptWeights=False,ptmin=-1.0)
    print "DONE, all the histograms for the theory bands were taken"

    gr = ROOT.TGraph(60)  # nPDFs
    grnomi = ROOT.TGraph(1)  # nPDFs

    # NptBisec = 1000
    # grBisec = ROOT.TGraph(NptBisec)  # nPDFs
    # for i in range(NptBisec):
    #     grBisec.SetPoint(i,)

    testbin = options.testbin
    xsecVsEta = histDict["xsec"][1]
    xsec2D    = histDict["xsec"][0] # to get total

    nPts = 0
    for key in xsecVsEta:
        ch,nuis = str(key[0]),str(key[1])
        if ch not in charges: continue
        if any(x in nuis for x in ["pdf","nominal"]): 
            pass
        else:
            continue
        #print key

        hist = xsecVsEta[key]
        histTotXsec = xsec2D[key] 
        #for ibin in range(1,1+hist.GetNbinsX()):
        #    print "bin {n}: xsec = {xs:.2f}".format(n=ibin,xs=hist.GetBinContent(ibin))
        xsecToPb = 35900. #*60400.
        xsecbin = hist.GetBinContent(testbin)/xsecToPb
        totxsec = histTotXsec.Integral(1,18,1,18)/xsecToPb # inside acceptance |18|*18 eta-pt bins           
        if "nominal" in nuis:
            print "N bins vs |eta| = %d" % hist.GetNbinsX()
            print "N bins vs |eta|-pt = %d-%d" % (histTotXsec.GetNbinsX(),histTotXsec.GetNbinsY())     
        
        if "pdf" in nuis:
            ipdf = utilities.getNFromString(nuis)
            gr.SetPoint(ipdf-1, xsecbin, totxsec)        
            nPts += 1
        else:
            grnomi.SetPoint(0, xsecbin, totxsec)
        #if xsecbin < 1000.0:
        #    print "DEBUG: ipdf = {ip}   xsec = {xs}".format(ip=ipdf,xs=str(xsecbin))

    print "{n} points in graph".format(n=nPts)
    
    canvas = ROOT.TCanvas("canvas","",1000,900)
    canvas.SetTickx(1)
    canvas.SetTicky(1)
    canvas.cd()
    canvas.SetLeftMargin(0.16)
    canvas.SetRightMargin(0.1)
    canvas.SetFillColor(0)
    canvas.SetGrid()
    canvas.cd()
    
    gr.SetTitle("PDFs: |#eta| bin {n}".format(n=testbin))
    gr.SetMarkerStyle(ROOT.kFullTriangleDown)
    gr.SetMarkerSize(2)
    gr.SetMarkerColor(ROOT.kRed+1)
    gr.Draw("AP")
    grnomi.SetMarkerStyle(ROOT.kFullCircle)
    grnomi.SetMarkerSize(2)
    grnomi.SetMarkerColor(ROOT.kBlue)
    grnomi.Draw("P SAME")
    gr.GetXaxis().SetTitleSize(0.05)
    gr.GetXaxis().SetLabelSize(0.04)
    gr.GetYaxis().SetTitleOffset(1.5)
    gr.GetYaxis().SetTitleSize(0.05)
    gr.GetYaxis().SetLabelSize(0.04)
    gr.GetXaxis().SetTitle("cross section [pb]")    
    gr.GetYaxis().SetTitle("total cross section [pb]")
    gr.GetXaxis().SetMaxDigits(3)
    gr.GetYaxis().SetMaxDigits(4);
    leg = ROOT.TLegend(0.2,0.6,0.5,0.9)
    leg.SetFillStyle(0)
    leg.SetFillColor(0)
    leg.SetBorderSize(0)
    leg.SetTextFont(62)
    leg.AddEntry(gr,"{n} PDFs".format(n=nPts),"P")
    leg.AddEntry(grnomi,"nominal","P")
    leg.Draw("same")


    #gr.GetXaxis().SetRangeUser(6500.0, 7500.0)

    for ext in ["png", "pdf"]:
        canvas.SaveAs("{od}/testPDFband_xsecVsEta_bin{b}.{ext}".format(od=outdir,b=testbin,ext=ext))

