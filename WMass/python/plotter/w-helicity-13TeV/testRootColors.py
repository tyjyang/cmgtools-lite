import ROOT, os, datetime, re, operator, math
from array import array
ROOT.gROOT.SetBatch(True)

import sys
sys.path.append(os.getcwd() + "/plotUtils/")
from utility import *


if __name__ == "__main__":

    ROOT.gStyle.SetOptStat(0)

    from optparse import OptionParser
    parser = OptionParser(usage='%prog [options] ')
    parser.add_option('-o','--outdir',     dest='outdir',     default='./',   type='string', help='output folder')
    parser.add_option('-f','--fill-style',     dest='fillStyle',     default='1001',   type='int', help='fill style, default is solid')
    (options, args) = parser.parse_args()

    if options.outdir:
        if not os.path.isdir(options.outdir):
            os.system('mkdir -p {od}'.format(od=options.outdir))
        os.system('cp {pf} {od}'.format(pf='/afs/cern.ch/user/g/gpetrucc/php/index.php',od=options.outdir))

    colorsLong = [ROOT.kPink,
                  ROOT.kSpring, 
                  ROOT.kViolet, 
                  ROOT.kAzure, 
                  ROOT.kOrange, 
                  ROOT.kTeal]

    colorsLongName = ["pink", "spring", "violet", "azure", "orange", "teal"]

    colorsShort = [ROOT.kCyan,
                   ROOT.kMagenta,
                   ROOT.kYellow, 
                   ROOT.kGreen, 
                   ROOT.kRed, 
                   ROOT.kBlue]

    colorsShortName = ["cyan", "magenta", "yellow", "green", "red", "yellow"]

    stackL = ROOT.THStack("stackL","")
    stackS = ROOT.THStack("stackS","")


    histL = {}
    hL_start = -9.5 
    nhL = 20
    for ic,c in enumerate(colorsLong):  # ic needed if not stacking, to have larger bin content based on color
        # 20, -9.5,10.5
        for nc in range(nhL):
            icL = int(hL_start + float(nc) + 0.5)
            histL[(c,icL)] = ROOT.TH1F("hist_colorLong_%d_icL_%d" % (c,icL), "", nhL, hL_start, hL_start+float(nhL))
            histL[(c,icL)].Fill(icL,1) 
            histL[(c,icL)].SetFillColor(c+icL)
            histL[(c,icL)].SetFillStyle(options.fillStyle)
            histL[(c,icL)].SetLineColor(c+icL)
            stackL.Add(histL[(c,icL)])

    histS = {}
    hS_start = -10.5 
    nhS = 15
    for ic,c in enumerate(colorsShort):
        # 20, -9.5,10.5
        for nc in range(nhS):
            icS = int(hS_start + float(nc) + 0.5)
            histS[(c,icS)] = ROOT.TH1F("hist_colorShort_%d_icS_%d" % (c,icS), "", nhS, hS_start, hS_start+float(nhS))
            histS[(c,icS)].Fill(icS,1) 
            histS[(c,icS)].SetFillColor(c+icS)
            histS[(c,icS)].SetFillStyle(options.fillStyle)
            histS[(c,icS)].SetLineColor(c+icS)
            stackS.Add(histS[(c,icS)])

    canvas = ROOT.TCanvas("canvas","",1800,900)
    canvas.SetTickx(1)
    canvas.SetTicky(1)
    canvas.cd()
    canvas.SetLeftMargin(0.05)
    canvas.SetRightMargin(0.04)
    canvas.SetTopMargin(0.1)
    canvas.SetBottomMargin(0.12)
    canvas.cd()

    title = "colors (1-{n}): ".format(n=len(colorsLongName)) + ",  ".join("{c}".format(c=colorsLongName[i]) for i in range(len(colorsLongName))) + " (fill style {f})".format(f=options.fillStyle)
    hdummyL = ROOT.TH1F("hdummyL",title,1,hL_start,hL_start+float(nhL))
    hdummyL.SetBinContent(1,len(colorsLong))
    hdummyL.SetLineColor(ROOT.kBlack)
    hdummyL.SetLineWidth(2)
    hdummyL.GetXaxis().SetTitle("color index (based on ROOT color wheel)")
    hdummyL.GetXaxis().SetTitleOffset(1.2)
    hdummyL.GetXaxis().SetTitleSize(0.05)
    hdummyL.GetXaxis().SetLabelSize(0.04)
    hdummyL.GetYaxis().SetTitle("a.u.")
    hdummyL.GetYaxis().SetTitleOffset(0.5)
    hdummyL.GetYaxis().SetTitleSize(0.05)
    hdummyL.GetYaxis().SetLabelSize(0.04)
    hdummyL.GetXaxis().SetRangeUser(hL_start,hL_start+float(nhL))
    hdummyL.GetYaxis().SetRangeUser(0.0,6.0)
    # force drawing stat box
    hdummyL.SetStats(1)
    hdummyL.Draw("HIST")
    stackL.Draw("HIST SAME")
    canvas.RedrawAxis("sameaxis")
    setTDRStyle()
    ROOT.gStyle.SetOptTitle(1)
    for ext in ["png","pdf"]:
        canvas.SaveAs("{o}/colorsLong_fillStyle_{f}.{ext}".format(o=options.outdir,f=options.fillStyle,ext=ext))

    title = "colors (1-{n}): ".format(n=len(colorsShortName)) + ",  ".join("{c}".format(c=colorsShortName[i]) for i in range(len(colorsShortName))) + " (fill style {f})".format(f=options.fillStyle)
    hdummyS = ROOT.TH1F("hdummyS",title,1,hS_start,hS_start+float(nhS))
    hdummyS.SetBinContent(1,len(colorsShort))
    hdummyS.SetLineColor(ROOT.kBlack)
    hdummyS.SetLineWidth(2)
    hdummyS.GetXaxis().SetTitle("color index (based on ROOT color wheel)")
    hdummyS.GetXaxis().SetTitleOffset(1.2)
    hdummyS.GetXaxis().SetTitleSize(0.05)
    hdummyS.GetXaxis().SetLabelSize(0.04)
    hdummyS.GetYaxis().SetTitle("a.u.")
    hdummyS.GetYaxis().SetTitleOffset(0.5)
    hdummyS.GetYaxis().SetTitleSize(0.05)
    hdummyS.GetYaxis().SetLabelSize(0.04)
    hdummyS.GetXaxis().SetRangeUser(hS_start,hS_start+float(nhS))
    hdummyS.GetYaxis().SetRangeUser(0.0,6.0)
    # force drawing stat box
    hdummyS.SetStats(1)
    hdummyS.Draw("HIST")
    stackS.Draw("HIST SAME")
    canvas.RedrawAxis("sameaxis")
    setTDRStyle()
    ROOT.gStyle.SetOptTitle(1)
    for ext in ["png","pdf"]:
        canvas.SaveAs("{o}/colorsShort_fillStyle_{f}.{ext}".format(o=options.outdir,f=options.fillStyle,ext=ext))

