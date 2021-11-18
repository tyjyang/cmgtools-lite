#!/usr/bin/env python3

# plot ratios of histogram variation over nominal
# example
# python w-mass-13TeV/makeSystRatios.py cards/wmass_Wnnpdf30_noPDFonZ/Wmunu_plus_shapes_pseudodataWmunuFromNNPDF31.root plots/testNanoAOD/WmassPlots/Wnnpdf30_noPDFonZ/testFit_pseudodataWmunuFromNNPDF31_alphaS1sigma/makeSystRatios/ -s ".*mass.*100MeV.*" -p "Wmunu_plus,Wtaunu_plus" -u

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

def plotUnrolledHistogram(h, process, syst, outdir, canvas, hist2DforBins):

    canvas.cd()
    yTitleOffset = 0.7
    h.SetTitle(f"{process},      {syst}")
    h.GetXaxis().SetTitleSize(0.05)
    h.GetXaxis().SetTitleOffset(1.05)
    h.GetXaxis().SetLabelSize(0.04)
    h.GetXaxis().SetTitle("Unrolled muon #eta-p_{T} bin")
    h.GetYaxis().SetTitleSize(0.05)
    h.GetYaxis().SetLabelSize(0.04)
    h.GetYaxis().SetTitleOffset(yTitleOffset)
    h.GetYaxis().SetTitle(f"syst/nomi")
    h.GetYaxis().SetTickSize(0.01)

    ########
    h.SetLineColor(ROOT.kBlack)
    h.SetMarkerColor(ROOT.kBlack)
    #h.SetMarkerStyle(20)
    #h.SetMarkerSize(0.9)
    h.SetMarkerStyle(1)
    #h.SetFillColor(ROOT.kGray+3)
    #h.SetFillColorAlpha(ROOT.kGray+3, 0.4)
    
    miny, maxy = getMinMaxHisto(h)
    diff = maxy - miny
    #print(f"min,max = {miny}, {maxy}:   diff = {diff}")
    #print(f"bin min,max = {h.GetMinimumBin()}, {h.GetMaximumBin()}")
    miny -= diff * 0.1
    maxy += diff * 0.1
    #print(f"new min,max = {miny}, {maxy}")
    h.GetYaxis().SetRangeUser(miny, maxy)

    h.Draw("HIST")

    bintext = ROOT.TLatex()
    textSize = 0.04
    textAngle = 30
    bintext.SetTextSize(textSize)
    bintext.SetTextFont(42)
    bintext.SetTextAngle(textAngle)

    # to draw panels in the unrolled plots
    nptBins = int(hist2DforBins.GetNbinsY())
    etarange = float(hist2DforBins.GetNbinsX())

    vertline = ROOT.TLine(36, 0, 36, canvas.GetUymax())
    vertline.SetLineColor(ROOT.kBlack)
    vertline.SetLineStyle(3) # 2 larger hatches
    for i in range(1, nptBins): # do not need line at canvas borders
        vertline.DrawLine(etarange*i, miny, etarange*i, maxy)

    ptBinRanges = []
    for ipt in range(0, nptBins):
        ptBinRanges.append("#splitline{{[{ptmin},{ptmax}]}}{{GeV}}".format(ptmin=int(hist2DforBins.GetYaxis().GetBinLowEdge(ipt+1)),
                                                                           ptmax=int(hist2DforBins.GetYaxis().GetBinLowEdge(ipt+2))))
    offsetText = etarange / 4.0
    for i in range(0,len(ptBinRanges)): # we need nptBins texts
        bintext.DrawLatex(etarange*i + offsetText, maxy - 0.2*diff, ptBinRanges[i])        
        
    line = ROOT.TF1("horiz_line", "1", h.GetXaxis().GetBinLowEdge(1), h.GetXaxis().GetBinLowEdge(h.GetNbinsX()+1))
    line.SetLineWidth(1)
    line.SetLineColor(ROOT.kRed)
    line.Draw("Lsame")
        
    canvas.RedrawAxis("sameaxis")

    for ext in ["png","pdf"]:
        canvas.SaveAs(f"{outdir}unrolled_{process}_{syst}.{ext}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("rootfile", type=str, nargs=1, help="Input root file")
    parser.add_argument("outdir",   type=str, nargs=1, help="Folder for plots")
    parser.add_argument("-s", "--systematics",    type=str, default=".*pdf.*", help="Comma separated list of regular expressions to select systematics to make ratios with nominal")
    parser.add_argument("-p", "--processes",    type=str, default="Wmunu_plus", help="Comma separated list of processes to plot (full name please)")
    parser.add_argument(     '--nContours', dest='nContours',    default=51, type=int, help='Number of contours in palette. Default is 51 (let it be odd, so the central strip is white if not using --abs-value and the range is symmetric)')
    parser.add_argument(     '--palette'  , dest='palette',      default=55, type=int, help='Set palette: use 0 for a built-in one, 55 is kRainbow')
    parser.add_argument(     '--invertPalette', dest='invertePalette' , default=False , action='store_true',   help='Inverte color ordering in palette')
    parser.add_argument( "-u", '--unrolled', dest='doUnrolled' , default=False , action='store_true',   help='Add unrolled plot too')
    args = parser.parse_args()

    fname = args.rootfile[0]
    outdir = args.outdir[0] + "/"
    if outdir and not os.path.exists(outdir):
        os.makedirs(outdir)
        htmlpath = "./templates/index.php"
        shutil.copy(htmlpath, outdir)

    processes = args.processes.split(',')
    regexp_syst = re.compile(args.systematics.replace(',','|'))

    nominals = {p : None for p in processes}
    #print(nominals)
    ratios = {p : [] for p in processes}

    canvas = ROOT.TCanvas("canvas","",900,800) 

    canvas_unroll = ROOT.TCanvas("canvas_unroll","",3000,800) 
    leftMargin = 0.06
    rightMargin = 0.01
    bottomMargin = 0.12
    canvas_unroll.SetTickx(1)
    canvas_unroll.SetTicky(1)
    canvas_unroll.cd()
    canvas_unroll.SetLeftMargin(leftMargin)
    canvas_unroll.SetRightMargin(rightMargin)
    canvas_unroll.cd()
    canvas_unroll.SetBottomMargin(bottomMargin)

    f = ROOT.TFile.Open(fname)
    if not f or not f.IsOpen():
        raise RuntimeError(f"Error when opening file {fname}")
    # get nominals
    for p in processes:
        nominals[p] = f.Get(f"x_{p}")
        if not nominals[p]:
            print(f"Error getting nominal histogram for process {p}")
            quit()
        else:
            nominals[p].SetDirectory(0)
    # now make a full loop for systematics
    for k in f.GetListOfKeys():
        name = k.GetName()
        if not regexp_syst.match(name): continue
        tokens = name.split("_") # remove "x_" and name of nuisance
        #print(tokens)
        pname = "_".join(tokens[1:-1])
        if pname not in processes: continue
        sname = tokens[-1]
        ratio = f.Get(name)
        ratio.SetDirectory(0)
        #ratios[pname].append(htmp.Divide(nominals[pname]))
        ratio.Divide(nominals[pname])
        ratio.SetTitle(f"syst: {sname}")
        drawCorrelationPlot(ratio, "Muon #eta", "Muon p_{T} (GeV)", f"{pname}: syst / nominal",
                            name, plotLabel="ForceTitle", outdir=outdir,
                            smoothPlot=False, drawProfileX=False, scaleToUnitArea=False, draw_both0_noLog1_onlyLog2=1,
                            palette=args.palette, nContours=args.nContours, invertePalette=args.invertePalette,
                            passCanvas=canvas, drawOption="COLZ0")
        if args.doUnrolled:
            ratio_unrolled = unroll2Dto1D(ratio, newname=f"unrolled_{name}", cropNegativeBins=False)
            plotUnrolledHistogram(ratio_unrolled, pname, sname, outdir, canvas_unroll, ratio)
            
    print()
