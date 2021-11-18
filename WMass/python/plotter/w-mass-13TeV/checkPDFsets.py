#!/usr/bin/env python3

# script to check PDF bands and also alpha, using different PDF sets. The script currently relies on having produced a root file with all histograms, typically it would contain TH3/TH2 made with mcPlots.py as for the standard analysis before using makeHistogramsWMass.py 
# this should probably have different alphaS variations for different sets, we choose 0.001 variation but scale it to correspond to 0.0015 (1 sigma variation from global pdf fits)
# nominal alphaS should always be 0.118

# example
# python w-mass-13TeV/checkPDFsets.py plots/testNanoAOD/testNtuplesAltPDF/testCode_TH3_chargePlus/plots_wmass_testPDF.root -o plots/testNanoAOD/testNtuplesAltPDF/testCode_TH3_chargePlus/checkPDFsets/ -c plus -p Wmunu_plus_postVFP

import os, re
import argparse
from array import array

from rollingFunctions import roll1Dto2D
from postFitHistograms import dressed2DfromFit, unroll2Dto1D, chargeUnrolledBinShifts

from checkTheoryUncertaintyOnShapes import quadSumVariationHisto, setHistErrorFromHisto

import sys
args = sys.argv[:]
sys.argv = ['-b']
import ROOT
sys.argv = args
ROOT.gROOT.SetBatch(True)
ROOT.PyConfig.IgnoreCommandLineOptions = True

import utilities
utilities = utilities.util()

sys.path.append(os.getcwd() + "/plotUtils/")
from utility import *

ROOT.gStyle.SetOptStat(0)
#ROOT.gStyle.SetOptTitle(1)
ROOT.TH1.SetDefaultSumw2()

def getGraphAsymmErrorsFromHistograms(hcen, hyerrUp, hyerrDown, scaleXerr=1.0, name=None):
    # if scaleXerr in [0,1], use this fraction on half bin width along X as bin error.
    # 0 means seeting X error to 0
    nbins = hcen.GetNbinsX()
    gr = ROOT.TGraphAsymmErrors(nbins)
    for ib in range(nbins):
        ibPlus1 = ib + 1
        gr.SetPoint(ib, hcen.GetBinCenter(ibPlus1), hcen.GetBinContent(ibPlus1))  
        xerr = scaleXerr * hcen.GetBinWidth(ibPlus1) / 2.0
        gr.SetPointEXhigh(ib, xerr)
        gr.SetPointEXlow( ib, xerr)
        gr.SetPointEYhigh(ib, hyerrUp.GetBinContent(ibPlus1))
        gr.SetPointEYlow(ib,  hyerrDown.GetBinContent(ibPlus1))
    if name:
        gr.SetName(name)
    return gr

def makeRatioGraphWithHisto(gr, h1):
    for i in range(h1.GetNbinsX()):
        xval = h1.GetBinCenter(i+1)
        yval = h1.GetBinContent(i+1)
        if yval == 0.0:
            print(f"In makeRatioGraphWithHisto(): y_graph = {gr.Eval(xval)},  y_histo = {yval}")
            quit()
        gr.SetPoint(i, xval, gr.Eval(xval)/yval)
        gr.SetPointEYhigh(i, gr.GetErrorYhigh(i)/yval)
        gr.SetPointEYlow(i,  gr.GetErrorYlow(i)/yval)
    return gr

drawStylePDF = {"NNPDF3.1": {"fillcolor" : ROOT.kAzure+7,
                             "linecolor" : ROOT.kAzure+7,
                             "fillstyle" : 1001,
                             "fillcoloralpha" : 0.9},
                "NNPDF3.0": {"fillcolor" : ROOT.kYellow+1,
                             "linecolor" : ROOT.kYellow+2,
                             "fillstyle" : 1001,
                             "fillcoloralpha" : 0.8},
                "CT18"    : {"fillcolor" : ROOT.kGreen+2,
                             "linecolor" : ROOT.kGreen+3,
                             "fillstyle" : 1001,
                             "fillcoloralpha" : 0.5},
                "CT18Z"   : {"fillcolor" : ROOT.kOrange+7,
                             "linecolor" : ROOT.kOrange+2,
                             "fillstyle" : 1001,
                             "fillcoloralpha" : 0.5}
}


def setPDFplotStyle(h, pdfset):
    h.SetLineColor(drawStylePDF[pdfset]["linecolor"])
    h.SetFillColor(drawStylePDF[pdfset]["fillcolor"])
    h.SetFillStyle(drawStylePDF[pdfset]["fillstyle"])
    h.SetFillColorAlpha(drawStylePDF[pdfset]["fillcolor"], drawStylePDF[pdfset]["fillcoloralpha"])


def scaleAlphaHist(halpha, hnomi, scale):
    halpha.Add(hnomi, -1.0)
    halpha.Scale(scale)
    halpha.Add(hnomi)
    return halpha
    
    
def comparePDFandAlpha(hstat, hpdf, halphaUp, halphaDown,
                       process, outdir, canvas, hist2DforBins, postfix,
                       PDFset="NNPDF3.1", nomiPDFset="NNPDF3.1"):

    canvas.cd()
    yTitleOffset = 0.75
    hstat.SetTitle(f"{process},    prefit,    {PDFset} PDF and #alpha_{{S}} (1_{{ }}#sigma variation)")
    hstat.GetXaxis().SetLabelSize(0)
    hstat.GetXaxis().SetTitle("")  
    hstat.GetYaxis().SetTitleSize(0.05)
    hstat.GetYaxis().SetLabelSize(0.04)
    hstat.GetYaxis().SetTitleOffset(yTitleOffset)
    hstat.GetYaxis().SetTitle("Events")

    ########
    hstat.SetLineColor(ROOT.kBlack)
    hstat.SetMarkerColor(ROOT.kBlack)
    #hstat.SetMarkerStyle(20)
    #hstat.SetMarkerSize(0.9)
    hstat.SetMarkerStyle(1)
    hstat.SetFillColor(ROOT.kGray+3)
    hstat.SetFillColorAlpha(ROOT.kGray+3, 0.4)
    hstat.Draw("E2")
    # pdf NNPDF3.1
    setPDFplotStyle(hpdf, PDFset)
    #hpdf.Draw("E2 SAME") # skip drawing since invisible
    # pdf NNPDF3.0
    halphaUp.SetLineColor(ROOT.kPink-6)
    halphaUp.SetLineWidth(1)
    halphaUp.SetFillColor(0)
    halphaUp.SetMarkerStyle(0)
    #halphaUp.Draw("HIST SAME")
    halphaDown.SetLineColor(ROOT.kPink-4)
    halphaDown.SetLineWidth(1)
    halphaDown.SetFillColor(0)
    halphaDown.SetMarkerStyle(0)
    #halphaDown.Draw("HIST SAME")

    miny, maxy = getMinMaxMultiHisto([hstat, hpdf, halphaUp, halphaDown])
    diff = maxy - miny
    miny -= diff * 0.1
    maxy += diff * 0.25
    hstat.GetYaxis().SetRangeUser(miny, maxy)

    hstat.Draw("PL SAME")

    bintext = ROOT.TLatex()
    textSize = 0.025
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
        vertline.DrawLine(etarange*i, 0, etarange*i, maxy)

    ptBinRanges = []
    for ipt in range(0, nptBins):
        ptBinRanges.append("#splitline{{[{ptmin},{ptmax}]}}{{GeV}}".format(ptmin=int(hist2DforBins.GetYaxis().GetBinLowEdge(ipt+1)),
                                                                           ptmax=int(hist2DforBins.GetYaxis().GetBinLowEdge(ipt+2))))
    offsetText = etarange / 4.0
    for i in range(0,len(ptBinRanges)): # we need nptBins texts
        bintext.DrawLatex(etarange*i + offsetText, 0.8*maxy, ptBinRanges[i])        

    canvas.RedrawAxis("sameaxis")

    pad2 = ROOT.TPad("pad2","pad2", 0.0, 0.0, 1, 0.95)
    pad2.SetTickx(1)
    pad2.SetTicky(1)
    pad2.SetTopMargin(1-lowerPanelHeight)
    pad2.SetRightMargin(rightMargin)
    pad2.SetLeftMargin(leftMargin)
    pad2.SetBottomMargin(bottomMargin)
    pad2.SetFillColor(0)
    pad2.SetGridy(1)
    pad2.SetFillStyle(0)
    pad2.Draw()
    pad2.cd()
    ratio_pdf = copy.deepcopy(hpdf.Clone("ratio_pdf"))
    ratio_alphaUp = copy.deepcopy(halphaUp.Clone("ratio_alphaUp"))
    ratio_alphaDown = copy.deepcopy(halphaDown.Clone("ratio_alphaDown"))
    den_noerr = copy.deepcopy(hstat.Clone("den_noerr"))
    den = copy.deepcopy(hstat.Clone("den"))
    for iBin in range (1, den_noerr.GetNbinsX()+1):
        den_noerr.SetBinError(iBin, 0.)
    if ratio_pdf.InheritsFrom("TH1"):
        ratio_pdf.Divide(den_noerr)
    else:
        ratio_pdf = makeRatioGraphWithHisto(ratio_pdf, den_noerr)
    ratio_alphaUp.Divide(den_noerr)
    ratio_alphaDown.Divide(den_noerr)
    den.Divide(den_noerr)
    miny, maxy = getMinMaxMultiHisto([den, ratio_pdf, ratio_alphaUp, ratio_alphaDown])
    # print(f"ratio min and max = {miny} and {maxy}")
    den.SetTitle("")
    den.Draw("E2")
    den.GetXaxis().SetTitleSize(0.05)
    den.GetXaxis().SetTitleOffset(1.05)
    den.GetXaxis().SetLabelSize(0.04)
    den.GetXaxis().SetTitle("Unrolled muon #eta-p_{T} bin")
    den.GetYaxis().SetTitleSize(0.05)
    den.GetYaxis().SetLabelSize(0.04)
    den.GetYaxis().SetTitleOffset(yTitleOffset)
    den.GetYaxis().SetTitle(f"Ratio with {nomiPDFset}")

    diff = maxy - miny
    miny -= 0.1 * diff
    maxy += 0.4 * diff
    den.GetYaxis().SetRangeUser(miny, maxy)
    ratio_pdf.Draw("E2 SAME")
    ratio_alphaUp.Draw("HIST SAME")
    ratio_alphaDown.Draw("HIST SAME")
    den.Draw("E2 SAME")
    
    vertline2 = ROOT.TLine(36, 0, 36, pad2.GetUymax())
    vertline2.SetLineColor(ROOT.kBlack)
    vertline2.SetLineStyle(3) # 2 larger hatches
    for i in range(1, nptBins): # do not need line at canvas borders
        vertline2.DrawLine(etarange*i, miny, etarange*i, maxy)

    ylowLegend = 0.47# 0.40 for 2 columns
    leg = ROOT.TLegend(0.18, ylowLegend, 0.92, 0.54)
    #leg.SetFillStyle(0)
    leg.SetFillColor(0)
    #leg.SetFillColorAlpha(0,0.6)
    leg.SetBorderSize(0)
    #leg.SetTextFont(62)
    leg.AddEntry(den, f"stat ({nomiPDFset})", "PLF")
    leg.AddEntry(ratio_pdf, f"PDF ({PDFset})","LF")
    leg.AddEntry(ratio_alphaUp, f"#alpha_{{S}} Up ({PDFset})","L")
    leg.AddEntry(ratio_alphaDown, f"#alpha_{{S}} Dn ({PDFset})","L")
    leg.SetNColumns(4)
    leg.Draw("same")
    pad2.RedrawAxis("sameaxis")

    PDFsetNoDot = PDFset.replace(".","")
    for ext in ["png","pdf"]:
        canvas.SaveAs(f"{outdir}{process}_pdfAndAlphaUncertainty_{PDFsetNoDot}{postfix}.{ext}")


def comparePDFhessianAndNomi(hnomi, hpdfUp, hpdfDown,
                             process, outdir, canvas, hist2DforBins, postfix,
                             nomiPDFset="NNPDF3.1", hessian=1, xAxisName=None,
                             leftMargin=0.16, rightMargin=0.05, bottomMargin=0.12, lowerPanelHeight=0.3):

    canvas.cd()
    yTitleOffset = 0.75
    yTitleSize = 0.05
    yLabelSize = 0.04
    if not hist2DforBins:
        yTitleOffset = 1.4
        yTitleSize = 0.05
        yLabelSize = 0.04
        
    hnomi.SetTitle(f"{process},    prefit,    {nomiPDFset} PDFs, Hessian {hessian}")
    hnomi.GetXaxis().SetLabelSize(0)
    hnomi.GetXaxis().SetTitle("")  
    hnomi.GetYaxis().SetTitleSize(yTitleSize)
    hnomi.GetYaxis().SetLabelSize(yLabelSize)
    hnomi.GetYaxis().SetTitleOffset(yTitleOffset)
    hnomi.GetYaxis().SetTitle("Events")

    ########
    hnomi.SetLineColor(ROOT.kBlack)
    #hnomi.SetMarkerStyle(20)
    #hnomi.SetMarkerSize(0.9)
    if hist2DforBins:
        hnomi.SetMarkerColor(ROOT.kBlack)
        hnomi.SetMarkerStyle(1)
        hnomi.SetFillColor(ROOT.kGray+3)
        hnomi.SetFillColorAlpha(ROOT.kGray+3, 0.4)
        hnomi.Draw("E2")
    else:
        #hnomi.SetMarkerColor(ROOT.kBlack)
        #hnomi.SetMarkerStyle(20)
        hnomi.SetFillColor(ROOT.kWhite)
        #hnomi.SetFillColorAlpha(ROOT.kGray+3, 0.4)
        hnomi.Draw("HE")
    # pdf NNPDF3.1
    hpdfUp.SetLineColor(ROOT.kRed+2)
    hpdfUp.SetLineWidth(1)
    hpdfUp.SetFillColor(0)
    hpdfUp.SetMarkerStyle(0)
    #hpdfUp.Draw("HIST SAME")
    hpdfDown.SetLineColor(ROOT.kBlue)
    hpdfDown.SetLineWidth(1)
    hpdfDown.SetFillColor(0)
    hpdfDown.SetMarkerStyle(0)
    #hpdfDown.Draw("HIST SAME")
    if not hist2DforBins:
        hpdfUp.Draw("HIST SAME")
        hpdfDown.Draw("HIST SAME")

    miny, maxy = getMinMaxMultiHisto([hnomi, hpdfUp, hpdfDown])
    diff = maxy - miny
    miny -= diff * 0.1
    maxy += diff * 0.25
    if hist2DforBins:
        hnomi.GetYaxis().SetRangeUser(miny, maxy)
    else:
        hnomi.GetYaxis().SetRangeUser(0.8*miny, 1.2*maxy)
    hnomi.Draw("PL SAME")


    if hist2DforBins:
        bintext = ROOT.TLatex()
        textSize = 0.025
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
            vertline.DrawLine(etarange*i, 0, etarange*i, maxy)

        ptBinRanges = []
        for ipt in range(0, nptBins):
            ptBinRanges.append("#splitline{{[{ptmin},{ptmax}]}}{{GeV}}".format(ptmin=int(hist2DforBins.GetYaxis().GetBinLowEdge(ipt+1)),
                                                                               ptmax=int(hist2DforBins.GetYaxis().GetBinLowEdge(ipt+2))))
        offsetText = etarange / 4.0
        for i in range(0,len(ptBinRanges)): # we need nptBins texts
            bintext.DrawLatex(etarange*i + offsetText, 0.8*maxy, ptBinRanges[i])        

    canvas.RedrawAxis("sameaxis")

    pad2 = ROOT.TPad("pad2","pad2", 0.0, 0.0, 1, 0.95)
    pad2.SetTickx(1)
    pad2.SetTicky(1)
    pad2.SetTopMargin(1-lowerPanelHeight)
    pad2.SetRightMargin(rightMargin)
    pad2.SetLeftMargin(leftMargin)
    pad2.SetBottomMargin(bottomMargin) 
    pad2.SetFillColor(0)
    pad2.SetGridy(1)
    pad2.SetFillStyle(0)
    pad2.Draw()
    pad2.cd()
    ratio_pdfUp = copy.deepcopy(hpdfUp.Clone("ratio_pdfUp"))
    ratio_pdfDown = copy.deepcopy(hpdfDown.Clone("ratio_pdfDown"))
    den_noerr = copy.deepcopy(hnomi.Clone("den_noerr"))
    den = copy.deepcopy(hnomi.Clone("den"))
    for iBin in range (1, den_noerr.GetNbinsX()+1):
        den_noerr.SetBinError(iBin, 0.)
    ratio_pdfUp.Divide(den_noerr)
    ratio_pdfDown.Divide(den_noerr)
    den.Divide(den_noerr)
    den.SetFillColor(ROOT.kGray)
    den.SetMarkerColor(0)
    den.SetMarkerSize(0)
    miny, maxy = getMinMaxMultiHisto([den, ratio_pdfUp, ratio_pdfDown])
    den.SetTitle("")
    den.Draw("E2")
    den.GetXaxis().SetTitleSize(0.05)
    den.GetXaxis().SetTitleOffset(1.05)
    den.GetXaxis().SetLabelSize(0.04)
    if xAxisName:
        den.GetXaxis().SetTitle(xAxisName)
    else:
        den.GetXaxis().SetTitle("Unrolled muon #eta-p_{T} bin")
    den.GetYaxis().SetTitleSize(yTitleSize)
    den.GetYaxis().SetLabelSize(yLabelSize)
    den.GetYaxis().SetTitleOffset(yTitleOffset)
    den.GetYaxis().SetTitle(f"pdf/nomi")

    diff = maxy - miny
    miny -= 0.1 * diff
    maxy += 0.4 * diff
    den.GetYaxis().SetRangeUser(miny, maxy)
    ratio_pdfUp.Draw("HIST SAME")
    ratio_pdfDown.Draw("HIST SAME")

    if hist2DforBins:
        vertline2 = ROOT.TLine(36, 0, 36, pad2.GetUymax())
        vertline2.SetLineColor(ROOT.kBlack)
        vertline2.SetLineStyle(3) # 2 larger hatches
        for i in range(1, nptBins): # do not need line at canvas borders
            vertline2.DrawLine(etarange*i, miny, etarange*i, maxy)
    else:
        den.GetYaxis().SetNdivisions(5)
    
            
    ylowLegend = 0.47# 0.40 for 2 columns
    yhighLegend = 0.54
    if not hist2DforBins:
        ylowLegend = 0.83# 0.40 for 2 columns
        yhighLegend = 0.90
        
    leg = ROOT.TLegend(0.18, ylowLegend, 0.92, yhighLegend)
    #leg.SetFillStyle(0)
    leg.SetFillColor(0)
    #leg.SetFillColorAlpha(0,0.6)
    leg.SetBorderSize(0)
    #leg.SetTextFont(62)
    leg.AddEntry(den, f"stat ({nomiPDFset})", "PLF")
    leg.AddEntry(ratio_pdfUp, f"pdf {hessian} Up","L")
    leg.AddEntry(ratio_pdfDown, f"pdf {hessian} Dn","L")
    leg.SetNColumns(3)
    leg.Draw("same")
    pad2.RedrawAxis("sameaxis")

    PDFsetNoDot = nomiPDFset.replace(".","")
    savErrorLevel = ROOT.gErrorIgnoreLevel; ROOT.gErrorIgnoreLevel = ROOT.kError;
    for ext in ["png","pdf"]:
        canvas.SaveAs(f"{outdir}{process}_pdfHessian_{hessian}_{PDFsetNoDot}{postfix}.{ext}")
    ROOT.gErrorIgnoreLevel = savErrorLevel;


        

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("rootfile", type=str, nargs=1, help="Input file with prefit histograms")
    parser.add_argument('-o','--outdir', dest='outdir', default='.', type=str, help='output directory to save output')
    parser.add_argument('-c','--charge-channel', dest='chargeChannel', default='', choices=["plus", "minus"], type=str, help='Charge channel for which histograms were made')
    #parser.add_argument('-l','--lumi', default=36.3, type=float, help='Integrated luminosity to print in the plot')
    parser.add_argument('-p','--processes', default="Wmunu_plus_postVFP", type=str, help='Comma-separated list of processes for which the plot is requested')
    parser.add_argument('--postfix', default="", type=str, help='Postfix to be added to plots')
    # currently alpha is not added, the option is preliminary (it just edits the legend entry for now)
    parser.add_argument(     '--add-alpha', dest="addAlpha", action='store_true', help="Sum alpha in quadrature to PDF bands (default uses PDFs only)")
    # to add postfit shapes (makes sense for a fit to data)
    # check w-mass-13TeV/postFitHistograms.py to see how these works, it is a bit convoluted
    # these are currently unused, they were taken from w-mass-13TeV/checkTheoryUncertaintyOnShapes.py
    parser.add_argument('--add-postfit-shape', dest="addPostfitShape", default="", type=str, help='Add postfit shapes, taken from the root file passed here as argument')
    parser.add_argument('--fit-charges', dest='fitCharges', choices=['plus', 'minus', 'plus,minus'], default='plus,minus', type=str, help='Charges used in the fit output file')
    parser.add_argument('-m','--n-mask-chan', dest='nMaskedChannel', default=0, type=int, help='Number of masked channels in the fit for each charge (0 if not using masked channels because no signal POIs is used in the fit)')
    args = parser.parse_args()
    
    outdir = args.outdir
    addStringToEnd(outdir, "/", notAddIfEndswithMatch=True)
    createPlotDirAndCopyPhp(outdir)

    outdirPDFvsAlpha = outdir + "checkPDFvsAlpha/"
    createPlotDirAndCopyPhp(outdirPDFvsAlpha)

    
    procs = args.processes.split(",")
    if len(procs) == 0:
        print("Error: you must specify some processes using option -p")
        quit()

    leftMargin = 0.07
    rightMargin = 0.01
    bottomMargin = 0.12
    lowerPanelHeight = 0.55
    
    #adjustSettings_CMS_lumi()
    canvas = ROOT.TCanvas("canvas", "", 3000, 1200)
    #adjustSettings_CMS_lumi()
    canvas.SetTickx(1)
    canvas.SetTicky(1)
    canvas.cd()
    canvas.SetLeftMargin(leftMargin)
    canvas.SetRightMargin(rightMargin)
    #canvas.SetBottomMargin(0.15)
    canvas.SetBottomMargin(lowerPanelHeight)

    canvas1D = ROOT.TCanvas("canvas1D", "", 800, 800)
    #adjustSettings_CMS_lumi()
    canvas1D.SetTickx(1)
    canvas1D.SetTicky(1)
    canvas1D.SetLeftMargin(0.16)
    canvas1D.SetRightMargin(0.05)
    #canvas.SetBottomMargin(0.15)
    canvas1D.SetBottomMargin(0.3)

    canvas.cd()

    
    prefitCharge = args.chargeChannel
    postfixToAdd = "charge" + ("Plus" if prefitCharge == "plus" else "Minus")
    if args.postfix:
        args.postfix += f"_{postfixToAdd}"
    else:
        args.postfix = postfixToAdd
    postfix = f"_{args.postfix}" if args.postfix else ""

    infile = safeOpenFile(args.rootfile[0], mode="READ")

    for p in procs:
        
        h3_nnpdf31 = safeGetObject(infile, f"pdfNNPDF31__{p}")
        h2_nomi_nnpdf31_fromTH2 = safeGetObject(infile, f"nominal__{p}")
        h2_nomi_nnpdf31_statUnc = copy.deepcopy(h2_nomi_nnpdf31_fromTH2.Clone(f"etaPt_{p}_nomi_nnpdf31_statUnc"))
        h3_nnpdf30 = safeGetObject(infile, f"pdfNNPDF30__{p}")
        h3_CT18all = safeGetObject(infile, f"pdfCT18__{p}")
  
        # for NNPDF3.1 the nominal is already available as TH2, but let's extract it from the underflow bin (perhaps we can check it is filled correctly)
        # do not use 0,0 to pick underflow bin with TAxis::SetRange, but -1,0, see https://root.cern.ch/doc/master/classTAxis.html#aed523b084d6b3f24f6b1128d7810e199
        h2_nomi_nnpdf31 = getTH2fromTH3(h3_nnpdf31, f"etaPt_{p}_nomi_nnpdf31", -1, binEnd=0) # this is taking full integral of TH3 !! Strange, but for now I can use the following way to access the underflow
        h2_nomi_nnpdf31.Reset("ICESM")
        fillTH2fromTH3zbin(h2_nomi_nnpdf31, h3_nnpdf31, zbin=0)
        h2_pdfHessian_nnpdf31 = [getTH2fromTH3(h3_nnpdf31, f"etaPt_{p}_hessian{i}_nnpdf31", i, i) for i in range(1, 101)] # bin 102 and 103 are alphaS Down/Up with 0.002 variation, but we don't use them here
        h2_alphaS0117_nnpdf31 = safeGetObject(infile, f"alphaS0117NNPDF31__{p}")
        h2_alphaS0119_nnpdf31 = safeGetObject(infile, f"alphaS0119NNPDF31__{p}")
        
        # NNPDF3.0
        h2_nomi_nnpdf30 = getTH2fromTH3(h3_nnpdf30, f"etaPt_{p}_nomi_nnpdf30", -1, 0) # see comments above about underflow bin
        h2_nomi_nnpdf30.Reset("ICESM")
        fillTH2fromTH3zbin(h2_nomi_nnpdf30, h3_nnpdf30, zbin=0)
        h2_pdfHessian_nnpdf30 = [getTH2fromTH3(h3_nnpdf30, f"etaPt_{p}_hessian{i}_nnpdf30", i, i) for i in range(1, 101)] # bin 102 and 103 are alphaS Down/Up with 0.002 variation, but we don't use them here
        h2_alphaS0117_nnpdf30 = safeGetObject(infile, f"alphaS0117NNPDF30__{p}")
        h2_alphaS0119_nnpdf30 = safeGetObject(infile, f"alphaS0119NNPDF30__{p}")

        # for CT18 we have CT18 and CT18Z packed up in the same histogram with 122 bins
        # CT18: bin 1 is nominal, 2->59 hessian (58 in total), 60/61 alphaS Down/Up with 0.001 variation
        # CT18Z: as above starting from bin 62
        # note that hessians are not symmetric, so the logic is that 2 is Up, 3 is Down, 4 is Up, and so on (Up/Down are arbitrary in terms of being above or below the cross section)
        # also these are 90% uncertainty for PDFs: have to scale PDF variations down by Scale = 1/1.645 = 0.607903 for 90% --> 68% 
        #
        # CT18
        h2_nomi_CT18 = getTH2fromTH3(h3_CT18all, f"etaPt_{p}_nomi_CT18", 1, 1)
        h2_pdfHessian_CT18_Up   = [getTH2fromTH3(h3_CT18all, f"etaPt_{p}_hessian{i-1}_CT18", i, i) for i in range(2, 60, 2)]
        h2_pdfHessian_CT18_Down = [getTH2fromTH3(h3_CT18all, f"etaPt_{p}_hessian{i-1}_CT18", i, i) for i in range(3, 60, 2)]
        h2_alphaS0117_CT18 = getTH2fromTH3(h3_CT18all, f"etaPt_{p}_alphaS0117_CT18", 60, 60)
        h2_alphaS0119_CT18 = getTH2fromTH3(h3_CT18all, f"etaPt_{p}_alphaS0119_CT18", 61, 61)
        #
        # CT18Z
        h2_nomi_CT18Z = getTH2fromTH3(h3_CT18all, f"etaPt_{p}_nomi_CT18Z", 62, 62)
        h2_pdfHessian_CT18Z_Up   = [getTH2fromTH3(h3_CT18all, f"etaPt_{p}_hessian{i-62}_CT18Z", i, i) for i in range(63, 121, 2)]
        h2_pdfHessian_CT18Z_Down = [getTH2fromTH3(h3_CT18all, f"etaPt_{p}_hessian{i-62}_CT18Z", i, i) for i in range(64, 121, 2)]
        h2_alphaS0117_CT18Z = getTH2fromTH3(h3_CT18all, f"etaPt_{p}_alphaS0117_CT18Z", 121, 121)
        h2_alphaS0119_CT18Z = getTH2fromTH3(h3_CT18all, f"etaPt_{p}_alphaS0119_CT18Z", 122, 122)
        
        # scale diff by 1.5 since the 1 sigma uncertainty for alpha should be 0.0015, while the variation is 0.001 (nominal is 0.118)
        scaleAlphaHist(h2_alphaS0117_nnpdf31, h2_nomi_nnpdf31, 1.5) 
        scaleAlphaHist(h2_alphaS0119_nnpdf31, h2_nomi_nnpdf31, 1.5) 
        scaleAlphaHist(h2_alphaS0117_nnpdf30, h2_nomi_nnpdf30, 1.5) 
        scaleAlphaHist(h2_alphaS0119_nnpdf30, h2_nomi_nnpdf30, 1.5) 
        scaleAlphaHist(h2_alphaS0117_CT18, h2_nomi_CT18, 1.5) 
        scaleAlphaHist(h2_alphaS0119_CT18, h2_nomi_CT18, 1.5) 
        scaleAlphaHist(h2_alphaS0117_CT18Z, h2_nomi_CT18Z, 1.5) 
        scaleAlphaHist(h2_alphaS0119_CT18Z, h2_nomi_CT18Z, 1.5) 

        #print("Nominal NNPDF3.1 integral: %.2f vs %.2f" % (h2_nomi_nnpdf31_fromTH2.Integral(), h2_nomi_nnpdf31.Integral()))
        #print("Hessian 1 NNPDF3.1 integral: %.2f" % (h2_pdfHessian_nnpdf31[0].Integral()))
        h2_nomi_nnpdf31_fromTH2.Divide(h2_nomi_nnpdf31)
        canvas2D = ROOT.TCanvas("canvas2D", "canvas2D", 800, 800)
        drawCorrelationPlot(h2_nomi_nnpdf31_fromTH2,
                            "Muon #eta",
                            "Muon p_{T} (GeV)",
                            "ratio:  nominal/underflow",
                            f"testNominalRatio_{p}",
                            plotLabel="ForceTitle",
                            outdir=outdir,
                            passCanvas=canvas2D,
                            drawOption="COLZ0")

        canvas.cd()        
        
        # now make so to plot unrolled, and also pt and eta projections for easier check, and ratio to NNPDF3.1
        # for now neglect alpha, it may go in a different plot. Ratio is made wrt nominal (central value of NNPDF3.1)
        # we may decide to sum alpha in quadrature, but first one may want to see how large alpha is wrt quadrature sum of hessian variations
        # in the ratio we expect 4 bands for pdfs, and alpha would have additional lines
        
        pdfUnc_nnpdf31 = quadSumVariationHisto(h2_nomi_nnpdf31, h2_pdfHessian_nnpdf31, postfixName="pdfUnc_nnpdf31")
        setHistErrorFromHisto(h2_nomi_nnpdf31, pdfUnc_nnpdf31)
        h1_nomi_nnpdf31 = unroll2Dto1D(h2_nomi_nnpdf31, newname=f"unroll_{h2_nomi_nnpdf31.GetName()}", cropNegativeBins=False)
        h1_nomi_nnpdf31_statUnc = unroll2Dto1D(h2_nomi_nnpdf31_statUnc, newname=f"unroll_{h2_nomi_nnpdf31_statUnc.GetName()}", cropNegativeBins=False)
        h1_alphaS0117_nnpdf31 = unroll2Dto1D(h2_alphaS0117_nnpdf31, newname=f"unroll_{h2_alphaS0117_nnpdf31.GetName()}", cropNegativeBins=False)
        h1_alphaS0119_nnpdf31 = unroll2Dto1D(h2_alphaS0119_nnpdf31, newname=f"unroll_{h2_alphaS0119_nnpdf31.GetName()}", cropNegativeBins=False)
        
        pdfUnc_nnpdf30 = quadSumVariationHisto(h2_nomi_nnpdf30, h2_pdfHessian_nnpdf30, postfixName="pdfUnc_nnpdf30")
        setHistErrorFromHisto(h2_nomi_nnpdf30, pdfUnc_nnpdf30)
        h1_nomi_nnpdf30 = unroll2Dto1D(h2_nomi_nnpdf30, newname=f"unroll_{h2_nomi_nnpdf30.GetName()}", cropNegativeBins=False)
        h1_alphaS0117_nnpdf30 = unroll2Dto1D(h2_alphaS0117_nnpdf30, newname=f"unroll_{h2_alphaS0117_nnpdf30.GetName()}", cropNegativeBins=False)
        h1_alphaS0119_nnpdf30 = unroll2Dto1D(h2_alphaS0119_nnpdf30, newname=f"unroll_{h2_alphaS0119_nnpdf30.GetName()}", cropNegativeBins=False)
                
        # need a TGraph with Asymmetric errors here
        #
        # CT18
        pdfUnc_CT18_Up = quadSumVariationHisto(h2_nomi_CT18, h2_pdfHessian_CT18_Up, postfixName="pdfUnc_CT18_Up", scale=0.607903)
        pdfUnc_CT18_Down = quadSumVariationHisto(h2_nomi_CT18, h2_pdfHessian_CT18_Down, postfixName="pdfUnc_CT18_Down", scale=0.607903)
        unroll_pdfUnc_CT18_Up = unroll2Dto1D(pdfUnc_CT18_Up, newname=f"unroll_{pdfUnc_CT18_Up.GetName()}", cropNegativeBins=False)
        unroll_pdfUnc_CT18_Down = unroll2Dto1D(pdfUnc_CT18_Down, newname=f"unroll_{pdfUnc_CT18_Down.GetName()}", cropNegativeBins=False)
        h1_nomi_CT18 = unroll2Dto1D(h2_nomi_CT18, newname=f"unroll_{h2_nomi_CT18.GetName()}", cropNegativeBins=False)
        h1_alphaS0117_CT18 = unroll2Dto1D(h2_alphaS0117_CT18, newname=f"unroll_{h2_alphaS0117_CT18.GetName()}", cropNegativeBins=False)
        h1_alphaS0119_CT18 = unroll2Dto1D(h2_alphaS0119_CT18, newname=f"unroll_{h2_alphaS0119_CT18.GetName()}", cropNegativeBins=False)
        gr1_nomi_CT18 = getGraphAsymmErrorsFromHistograms(h1_nomi_CT18, unroll_pdfUnc_CT18_Up, unroll_pdfUnc_CT18_Down, name=f"unroll_graph_{h2_nomi_CT18.GetName()}")
        print()
        print()
        for ih,h in enumerate(h2_pdfHessian_CT18_Up):
            diffUp = h2_nomi_CT18.Integral() - h.Integral()
            diffDown = h2_nomi_CT18.Integral() - h2_pdfHessian_CT18_Down[ih].Integral()
            # check if variations are on the same side wrt nominal
            if diffUp * diffDown > 0:
                outdirHessian = outdir + "debugHessians_CT18/"
                createPlotDirAndCopyPhp(outdirHessian)
                where = "above" if diffUp < 0 else "below"
                print(f">>> CT18: warning with hessian {ih+1}: both variations {where} nominal")
                unroll_pdfHessian_CT18_Up = unroll2Dto1D(h2_pdfHessian_CT18_Up[ih], newname=f"unroll_pdfHessian{ih+1}_CT18_Up", cropNegativeBins=False)
                unroll_pdfHessian_CT18_Down = unroll2Dto1D(h2_pdfHessian_CT18_Down[ih], newname=f"unroll_pdfHessian{ih+1}_CT18_Down", cropNegativeBins=False)
                comparePDFhessianAndNomi(h1_nomi_CT18, unroll_pdfHessian_CT18_Up, unroll_pdfHessian_CT18_Down,
                                         p, outdirHessian, canvas, h2_nomi_nnpdf31_fromTH2, postfix,
                                         nomiPDFset="CT18", hessian=int(ih+1),
                                         leftMargin=leftMargin, rightMargin=rightMargin, bottomMargin=bottomMargin, lowerPanelHeight=lowerPanelHeight)
                h1eta_nomi_CT18 = h2_nomi_CT18.ProjectionX("px_eta",0,-1,"e")
                h1eta_pdfHessian_CT18_Up = h2_pdfHessian_CT18_Up[ih].ProjectionX("px_hessianUp_eta",0,-1,"e")
                h1eta_pdfHessian_CT18_Down = h2_pdfHessian_CT18_Down[ih].ProjectionX("px_hessianDown_eta",0,-1,"e")
                comparePDFhessianAndNomi(h1eta_nomi_CT18, h1eta_pdfHessian_CT18_Up, h1eta_pdfHessian_CT18_Down,
                                         p, outdirHessian, canvas1D, None, postfix+"_1Deta",
                                         nomiPDFset="CT18", hessian=int(ih+1), xAxisName="Muon #eta",
                                         leftMargin=0.16, rightMargin=0.05, bottomMargin=0.12, lowerPanelHeight=0.3)
                
        # CT18Z
        pdfUnc_CT18Z_Up = quadSumVariationHisto(h2_nomi_CT18Z, h2_pdfHessian_CT18Z_Up, postfixName="pdfUnc_CT18Z_Up", scale=0.607903)
        pdfUnc_CT18Z_Down = quadSumVariationHisto(h2_nomi_CT18Z, h2_pdfHessian_CT18Z_Down, postfixName="pdfUnc_CT18Z_Down", scale=0.607903)
        unroll_pdfUnc_CT18Z_Up = unroll2Dto1D(pdfUnc_CT18Z_Up, newname=f"unroll_{pdfUnc_CT18Z_Up.GetName()}", cropNegativeBins=False)
        unroll_pdfUnc_CT18Z_Down = unroll2Dto1D(pdfUnc_CT18Z_Down, newname=f"unroll_{pdfUnc_CT18Z_Down.GetName()}", cropNegativeBins=False)
        h1_nomi_CT18Z = unroll2Dto1D(h2_nomi_CT18Z, newname=f"unroll_{h2_nomi_CT18Z.GetName()}", cropNegativeBins=False)
        h1_alphaS0117_CT18Z = unroll2Dto1D(h2_alphaS0117_CT18Z, newname=f"unroll_{h2_alphaS0117_CT18Z.GetName()}", cropNegativeBins=False)
        h1_alphaS0119_CT18Z = unroll2Dto1D(h2_alphaS0119_CT18Z, newname=f"unroll_{h2_alphaS0119_CT18Z.GetName()}", cropNegativeBins=False)
        gr1_nomi_CT18Z = getGraphAsymmErrorsFromHistograms(h1_nomi_CT18Z, unroll_pdfUnc_CT18Z_Up, unroll_pdfUnc_CT18Z_Down, name=f"unroll_graph_{h2_nomi_CT18Z.GetName()}")
        # check for CT18 pdfs
        print()
        print()
        # check for CTZ18 pdfs
        for ih,h in enumerate(h2_pdfHessian_CT18Z_Up):
            diffUp = h2_nomi_CT18Z.Integral() - h.Integral()
            diffDown = h2_nomi_CT18Z.Integral() - h2_pdfHessian_CT18Z_Down[ih].Integral()
            # check if variations are on the same side wrt nominal
            if diffUp * diffDown > 0:
                outdirHessian = outdir + "debugHessians_CT18Z/"
                createPlotDirAndCopyPhp(outdirHessian)
                where = "above" if diffUp < 0 else "below"
                print(f">>> CT18Z: warning with hessian {ih+1}: both variations {where} nominal")
                unroll_pdfHessian_CT18Z_Up = unroll2Dto1D(h2_pdfHessian_CT18Z_Up[ih], newname=f"unroll_pdfHessian{ih+1}_CT18Z_Up", cropNegativeBins=False)
                unroll_pdfHessian_CT18Z_Down = unroll2Dto1D(h2_pdfHessian_CT18Z_Down[ih], newname=f"unroll_pdfHessian{ih+1}_CT18Z_Down", cropNegativeBins=False)
                comparePDFhessianAndNomi(h1_nomi_CT18Z, unroll_pdfHessian_CT18Z_Up, unroll_pdfHessian_CT18Z_Down,
                                         p, outdirHessian, canvas, h2_nomi_nnpdf31_fromTH2, postfix,
                                         nomiPDFset="CT18Z", hessian=int(ih+1),
                                         leftMargin=leftMargin, rightMargin=rightMargin, bottomMargin=bottomMargin, lowerPanelHeight=lowerPanelHeight)
                h1eta_nomi_CT18Z = h2_nomi_CT18Z.ProjectionX("px_eta",0,-1,"e")
                h1eta_pdfHessian_CT18Z_Up = h2_pdfHessian_CT18Z_Up[ih].ProjectionX("px_hessianUp_eta",0,-1,"e")
                h1eta_pdfHessian_CT18Z_Down = h2_pdfHessian_CT18Z_Down[ih].ProjectionX("px_hessianDown_eta",0,-1,"e")
                comparePDFhessianAndNomi(h1eta_nomi_CT18Z, h1eta_pdfHessian_CT18Z_Up, h1eta_pdfHessian_CT18Z_Down,
                                         p, outdirHessian, canvas1D, None, postfix+"_1Deta", 
                                         nomiPDFset="CT18Z", hessian=int(ih+1), xAxisName="Muon #eta",
                                         leftMargin=0.16, rightMargin=0.05, bottomMargin=0.12, lowerPanelHeight=0.3)
        print()
        print()



        
        canvas.cd()
        yTitleOffset = 0.75
        pdfText = "PDFs #oplus #alpha_{S}" if args.addAlpha else "PDFs only (no #alpha_{S})"
        h1_nomi_nnpdf31_statUnc.SetTitle(f"{p},    prefit,    {pdfText}")
        h1_nomi_nnpdf31_statUnc.GetXaxis().SetLabelSize(0)
        h1_nomi_nnpdf31_statUnc.GetXaxis().SetTitle("")  
        h1_nomi_nnpdf31_statUnc.GetYaxis().SetTitleSize(0.05)
        h1_nomi_nnpdf31_statUnc.GetYaxis().SetLabelSize(0.04)
        h1_nomi_nnpdf31_statUnc.GetYaxis().SetTitleOffset(yTitleOffset)
        h1_nomi_nnpdf31_statUnc.GetYaxis().SetTitle("Events")

        ########
        h1_nomi_nnpdf31_statUnc.Draw("E2")
        h1_nomi_nnpdf31_statUnc.SetLineColor(ROOT.kBlack)
        h1_nomi_nnpdf31_statUnc.SetLineWidth(1)
        h1_nomi_nnpdf31_statUnc.SetMarkerColor(ROOT.kBlack)
        #h1_nomi_nnpdf31_statUnc.SetMarkerStyle(20)
        #h1_nomi_nnpdf31_statUnc.SetMarkerSize(0.9)
        h1_nomi_nnpdf31_statUnc.SetMarkerStyle(1)
        h1_nomi_nnpdf31_statUnc.SetFillColor(ROOT.kGray+3)
        h1_nomi_nnpdf31_statUnc.SetFillColorAlpha(ROOT.kGray+3, 0.4)
        #h1_nomi_nnpdf31_statUnc.SetFillStyle(3002)
        # pdf NNPDF3.1
        setPDFplotStyle(h1_nomi_nnpdf31, "NNPDF3.1")
        #h1_nomi_nnpdf31.Draw("E2 SAME") # skip drawing since invisible
        # pdf NNPDF3.0
        setPDFplotStyle(h1_nomi_nnpdf30, "NNPDF3.0")
        #h1_nomi_nnpdf30.Draw("E2 SAME") # skip drawing since invisible
        # pdf CT18
        setPDFplotStyle(gr1_nomi_CT18, "CT18")
        #gr1_nomi_CT18.Draw("E2 SAME") # skip drawing since invisible
        # pdf CT18Z
        setPDFplotStyle(gr1_nomi_CT18Z, "CT18Z")
        #gr1_nomi_CT18Z.Draw("E2 SAME") # skip drawing since invisible

        miny, maxy = getMinMaxMultiHisto([h1_nomi_nnpdf31_statUnc, h1_nomi_nnpdf31, h1_nomi_nnpdf30, gr1_nomi_CT18, gr1_nomi_CT18Z])
        diff = maxy - miny
        miny -= diff * 0.1
        maxy += diff * 0.25
        h1_nomi_nnpdf31_statUnc.GetYaxis().SetRangeUser(miny, maxy)
        
        h1_nomi_nnpdf31_statUnc.Draw("PL SAME")
        
        bintext = ROOT.TLatex()
        textSize = 0.025
        textAngle = 30
        bintext.SetTextSize(textSize)
        bintext.SetTextFont(42)
        bintext.SetTextAngle(textAngle)

        # to draw panels in the unrolled plots
        nptBins = int(h2_nomi_nnpdf31_fromTH2.GetNbinsY())
        etarange = float(h2_nomi_nnpdf31_fromTH2.GetNbinsX())

        vertline = ROOT.TLine(36, 0, 36, canvas.GetUymax())
        vertline.SetLineColor(ROOT.kBlack)
        vertline.SetLineStyle(3) # 2 larger hatches
        for i in range(1, nptBins): # do not need line at canvas borders
            vertline.DrawLine(etarange*i, 0, etarange*i, maxy)

        ptBinRanges = []
        for ipt in range(0, nptBins):
            ptBinRanges.append("#splitline{{[{ptmin},{ptmax}]}}{{GeV}}".format(ptmin=int(h2_nomi_nnpdf31_fromTH2.GetYaxis().GetBinLowEdge(ipt+1)),
                                                                               ptmax=int(h2_nomi_nnpdf31_fromTH2.GetYaxis().GetBinLowEdge(ipt+2))))
        offsetText = etarange / 4.0
        for i in range(0,len(ptBinRanges)): # we need nptBins texts
            bintext.DrawLatex(etarange*i + offsetText, 0.8*maxy, ptBinRanges[i])        

        canvas.RedrawAxis("sameaxis")

        pad2 = ROOT.TPad("pad2","pad2", 0.0, 0.0, 1, 0.95)
        pad2.SetTickx(1)
        pad2.SetTicky(1)
        pad2.SetTopMargin(1-lowerPanelHeight)
        pad2.SetRightMargin(rightMargin)
        pad2.SetLeftMargin(leftMargin)
        pad2.SetBottomMargin(bottomMargin)
        pad2.SetFillColor(0)
        pad2.SetGridy(1)
        pad2.SetFillStyle(0)
        pad2.Draw()
        pad2.cd()
        ratio_nnpdf31 = copy.deepcopy(h1_nomi_nnpdf31.Clone("ratio_nnpdf31"))
        ratio_nnpdf30 = copy.deepcopy(h1_nomi_nnpdf30.Clone("ratio_nnpdf30"))
        ratio_CT18 = copy.deepcopy(gr1_nomi_CT18.Clone("ratio_CT18"))
        ratio_CT18Z = copy.deepcopy(gr1_nomi_CT18Z.Clone("ratio_CT18Z"))
        den_noerr = copy.deepcopy(h1_nomi_nnpdf31_statUnc.Clone("den_noerr"))
        den = copy.deepcopy(h1_nomi_nnpdf31_statUnc.Clone("den"))
        for iBin in range (1, den_noerr.GetNbinsX()+1):
            den_noerr.SetBinError(iBin, 0.)
        ratio_nnpdf31.Divide(den_noerr)
        ratio_nnpdf30.Divide(den_noerr)
        ratio_CT18 = makeRatioGraphWithHisto(ratio_CT18, den_noerr)
        ratio_CT18Z = makeRatioGraphWithHisto(ratio_CT18Z, den_noerr)
        den.Divide(den_noerr)
        miny, maxy = getMinMaxMultiHisto([den, ratio_nnpdf31, ratio_nnpdf30, ratio_CT18, ratio_CT18Z])
        # print(f"ratio min and max = {miny} and {maxy}")
        den.SetTitle("")
        den.Draw("E2")
        den.GetXaxis().SetTitleSize(0.05)
        den.GetXaxis().SetTitleOffset(1.05)
        den.GetXaxis().SetLabelSize(0.04)
        den.GetXaxis().SetTitle("Unrolled muon #eta-p_{T} bin")
        den.GetYaxis().SetTitleSize(0.05)
        den.GetYaxis().SetLabelSize(0.04)
        den.GetYaxis().SetTitleOffset(yTitleOffset)
        den.GetYaxis().SetTitle("Ratio with NNPDF3.1")

        diff = maxy - miny
        miny -= 0.1 * diff
        maxy += 0.4 * diff
        den.GetYaxis().SetRangeUser(miny, maxy)
        ratio_nnpdf31.Draw("E2 SAME")
        ratio_nnpdf30.Draw("E2 SAME")
        #ratio_CT18Z.Draw("LE2 SAME")
        ratio_CT18Z.Draw("L|| SAME")
        ratio_CT18.Draw("L|| SAME")
        #den.Draw("HE SAME")
        den.Draw("E2 SAME")
        vertline2 = ROOT.TLine(36, 0, 36, pad2.GetUymax())
        vertline2.SetLineColor(ROOT.kBlack)
        vertline2.SetLineStyle(3) # 2 larger hatches
        for i in range(1, nptBins): # do not need line at canvas borders
            vertline2.DrawLine(etarange*i, miny, etarange*i, maxy)
        
        ylowLegend = 0.47# 0.40 for 2 columns
        leg = ROOT.TLegend(0.18, ylowLegend, 0.92, 0.54)
        #leg.SetFillStyle(0)
        leg.SetFillColor(0)
        #leg.SetFillColorAlpha(0,0.6)
        #leg.SetBorderSize(0)
        #leg.SetTextFont(62)
        leg.AddEntry(den, "MC statistics", "PLEF")
        leg.AddEntry(ratio_nnpdf31, "NNPDF3.1","F")
        leg.AddEntry(ratio_nnpdf30, "NNPDF3.0","LF")
        leg.AddEntry(ratio_CT18, "CT18","LF")
        leg.AddEntry(ratio_CT18Z, "CT18Z","LF")
        leg.SetNColumns(5)
        leg.Draw("same")
        pad2.RedrawAxis("sameaxis")

        for ext in ["png","pdf"]:
            canvas.SaveAs(f"{outdir}{p}_pdfUncertainty{postfix}.{ext}")


        comparePDFandAlpha(h1_nomi_nnpdf31_statUnc, h1_nomi_nnpdf31, h1_alphaS0119_nnpdf31, h1_alphaS0117_nnpdf31, p, outdirPDFvsAlpha, canvas, h2_nomi_nnpdf31_fromTH2, postfix, PDFset="NNPDF3.1")
        comparePDFandAlpha(h1_nomi_nnpdf31_statUnc, h1_nomi_nnpdf30, h1_alphaS0119_nnpdf30, h1_alphaS0117_nnpdf30, p, outdirPDFvsAlpha, canvas, h2_nomi_nnpdf31_fromTH2, postfix, PDFset="NNPDF3.0")
        comparePDFandAlpha(h1_nomi_nnpdf31_statUnc, gr1_nomi_CT18, h1_alphaS0119_CT18, h1_alphaS0117_CT18, p, outdirPDFvsAlpha, canvas, h2_nomi_nnpdf31_fromTH2, postfix, PDFset="CT18")
        comparePDFandAlpha(h1_nomi_nnpdf31_statUnc, gr1_nomi_CT18Z, h1_alphaS0119_CT18Z, h1_alphaS0117_CT18Z, p, outdirPDFvsAlpha, canvas, h2_nomi_nnpdf31_fromTH2, postfix, PDFset="CT18Z")

        outfile = safeOpenFile(outdir+f"shapes_{p}.root", mode="RECREATE")
        outfile.cd()
        h2_nomi_nnpdf31.Write()
        h2_nomi_nnpdf30.Write()
        h2_nomi_CT18.Write()    
        h2_nomi_CT18Z.Write()
        outfile.Close()
        # restore previous file
        infile.cd()
        
    ########        
    infile.Close()

    
    

                        
                                 
        
