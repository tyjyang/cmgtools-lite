#!/usr/bin/env python3

# example
# python w-mass-13TeV/checkTheoryUncertaintyOnShapes.py cards/wmass_fixMassWeights_splitW/Wmunu_plus_shapes.root -o plots/testNanoAOD/WmassPlots_jetEta2p4_fixMassWeight_splitW/testUncertainty/theory/ -l 16.8 -p Wmunu_plus,data_fakes,Zmumu --add-postfit-shape cards/wmass_fixMassWeights_splitW/fit/hessian/fitresults_123456789_Asimov_clipSyst1p3_bbb1_cxs1.root --fit-charges "plus,minus" -m 0

import os, re
import argparse
from array import array

from rollingFunctions import roll1Dto2D
from postFitHistograms import dressed2DfromFit, unroll2Dto1D, chargeUnrolledBinShifts

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

def quadSumVariationHisto(hnomi, hvars, hname=None, postfixName="quadSumVar"):
    # hlist is supposed to contain a list of histograms for variations, for example pdf_i
    retname = hname if hname != None else f"{hnomi.GetName()}_{postfixName}"
    hret = copy.deepcopy(hnomi.Clone(retname))

    if hnomi.GetDimension() == 1:
        for ix in range(1, 1 + hnomi.GetNbinsX()):
            sumSquare = 0 
            for h in hvars:
                val = h.GetBinContent(ix) - hnomi.GetBinContent(ix)
                sumSquare += (val * val)
            hret.SetBinContent(ix, math.sqrt(sumSquare))

    elif hnomi.GetDimension() == 2:
        for ix in range(1, 1 + hnomi.GetNbinsX()):
            for iy in range(1, 1 + hnomi.GetNbinsY()):
                sumSquare = 0 
                for h in hvars:
                    val = h.GetBinContent(ix, iy) - hnomi.GetBinContent(ix, iy)
                    sumSquare += (val * val)
                hret.SetBinContent(ix, iy, math.sqrt(sumSquare))

    return hret


def setHistErrorFromHisto(hnomi, herr):
    if hnomi.GetDimension() == 1:
        for ix in range(1, 1 + hnomi.GetNbinsX()):
            hnomi.SetBinError(ix, herr.GetBinContent(ix))

    elif hnomi.GetDimension() == 2:
        for ix in range(1, 1 + hnomi.GetNbinsX()):
            for iy in range(1, 1 + homi.GetNbinsY()):
                hnomi.SetBinError(ix, iy, herr.GetBinContent(ix, iy))
    

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("rootfile", type=str, nargs=1, help="Input file with prefit histograms")
    parser.add_argument('-o','--outdir', dest='outdir', default='.', type=str, help='output directory to save output')
    parser.add_argument('-l','--lumi', default=36.3, type=float, help='Integrated luminosity to print in the plot')
    parser.add_argument('-p','--processes', default="", type=str, help='Comma-separated list of processes for which the plot is requested')
    parser.add_argument('--postfix', default="", type=str, help='Postfix to be added to plots')
    # to add postfit shapes (makes sense for a fit to data)
    # check w-mass-13TeV/postFitHistograms.py to see how these works, it is a bit convoluted
    parser.add_argument('--add-postfit-shape', dest="addPostfitShape", default="", type=str, help='Add postfit shapes, taken from the root file passed here as argument')
    parser.add_argument('--fit-charges', dest='fitCharges', choices=['plus', 'minus', 'plus,minus'], default='plus,minus', type=str, help='Charges used in the fit output file')
    parser.add_argument('-m','--n-mask-chan', dest='nMaskedChannel', default=0, type=int, help='Number of masked channels in the fit for each charge (0 if not using masked channels because no signal POIs is used in the fit)')
    args = parser.parse_args()
    
    outdir = args.outdir
    addStringToEnd(outdir, "/", notAddIfEndswithMatch=True)
    createPlotDirAndCopyPhp(outdir)

    procs = args.processes.split(",")
    if len(procs) == 0:
        print("Error: you must specify some processes using option -p")
        quit()

    postfitHists = {}

    hnomi = {}
    hpdf = {}
    halphaSUp = {}
    halphaSDown = {}
        
    #adjustSettings_CMS_lumi()
    canvas = ROOT.TCanvas("canvas", "", 800, 1200)
    #adjustSettings_CMS_lumi()
    canvas.SetTickx(1)
    canvas.SetTicky(1)
    canvas.cd()
    canvas.SetLeftMargin(0.15)
    canvas.SetRightMargin(0.04)
    canvas.cd()
    #canvas.SetBottomMargin(0.15)
    lowerPanelHeight = 0.55
    canvas.SetBottomMargin(lowerPanelHeight)

    prefitCharge = "plus" if "plus" in args.rootfile[0] else "minus"
    postfixToAdd = "charge" + ("Plus" if prefitCharge == "plus" else "Minus")
    if args.postfix:
        args.postfix += f"_{postfixToAdd}"
    else:
        args.postfix = postfixToAdd
    postfix = f"_{args.postfix}" if args.postfix else ""

    infile = safeOpenFile(args.rootfile[0], mode="READ")

    for p in procs:
        hnomi = {"2D" : None, "eta" : None, "pt": None}
        hpdf    = {"2D" : [], "eta" : [], "pt": []}
        halphaSUp = {"2D" : [], "eta" : [], "pt": []}
        halphaSDown = {"2D" : [], "eta" : [], "pt": []}
        hnomi["2D"] = safeGetObject(infile, f"x_{p}")
        hnomi["eta"] = hnomi["2D"].ProjectionX(hnomi["2D"].GetName()+"_eta", 0, -1, "e")
        hnomi["pt"] = hnomi["2D"].ProjectionY(hnomi["2D"].GetName()+"_pt", 0, -1, "e")
        for ipdf in range(1,101):
            hpdf["2D"].append(safeGetObject(infile, f"x_{p}_pdf{ipdf}Up"))
            hpdf["eta"].append(hpdf["2D"][-1].ProjectionX(hpdf["2D"][-1].GetName()+"_eta", 0, -1, "e"))
            hpdf["pt"].append(hpdf["2D"][-1].ProjectionY(hpdf["2D"][-1].GetName()+"_pt", 0, -1, "e"))
        
        halphaSUp["2D"] = safeGetObject(infile, f"x_{p}_alphaSUp")
        halphaSUp["eta"] = halphaSUp["2D"].ProjectionX(halphaSUp["2D"].GetName()+"_eta", 0, -1, "e")
        halphaSUp["pt"] = halphaSUp["2D"].ProjectionY(halphaSUp["2D"].GetName()+"_pt", 0, -1, "e")

        halphaSDown["2D"] = safeGetObject(infile, f"x_{p}_alphaSDown")
        halphaSDown["eta"] = halphaSDown["2D"].ProjectionX(halphaSDown["2D"].GetName()+"_eta", 0, -1, "e")
        halphaSDown["pt"] = halphaSDown["2D"].ProjectionY(halphaSDown["2D"].GetName()+"_pt", 0, -1, "e")

        hpfs = {"2D" : None, "eta" : None, "pt": None}
        
        if args.addPostfitShape:
            pfs_file = safeOpenFile(args.addPostfitShape, mode="READ")
            postfitHist = safeGetObject(pfs_file, f"expproc_{p}_postfit")
            charges = args.fitCharges.split(",")
            nCharges = len(charges)
            nMaskedChanPerCharge = args.nMaskedChannel # check if we actually have masked channels, we may not, default should be 0
            etabins = [round(hnomi["2D"].GetXaxis().GetBinLowEdge(i), 1) for i in range(1, 2 + hnomi["2D"].GetNbinsX())]
            ptbins = [round(hnomi["2D"].GetYaxis().GetBinLowEdge(i), 1) for i in range(1, 2 + hnomi["2D"].GetNbinsY())]
            recoBins = templateBinning(etabins, ptbins)
            nRecoBins = recoBins.NTotBins
            #following array is used to call function dressed2DfromFit()                                                                                                                           
            binning = [recoBins.Neta, recoBins.etaBins, recoBins.Npt, recoBins.ptBins]
            shifts = chargeUnrolledBinShifts(pfs_file, "mu", nCharges, nMaskedChanPerCharge)
            pfs_file.Close()
            #print(shifts)
            binshift = shifts[prefitCharge]
            hpfs["2D"] = dressed2DfromFit(postfitHist, binning, f"{p}_postfit", p, binshift, nMaskedCha=nMaskedChanPerCharge,
                                          isComb=False, # old option for helicity analysis, related to e-mu combination 
                                          isMuPlot=True,
                                          nRecoBins=nRecoBins)        
            hpfs["eta"] = hpfs["2D"].ProjectionX(hpfs["2D"].GetName()+"_eta", 0, -1, "e")
            hpfs["pt"] = hpfs["2D"].ProjectionY(hpfs["2D"].GetName()+"_pt", 0, -1, "e")

        pdfUnc = {"2D"  : quadSumVariationHisto(hnomi["2D"], hpdf["2D"], postfixName="pdf"),
                  "eta" : quadSumVariationHisto(hnomi["eta"], hpdf["eta"], postfixName="pdf"),
                  "pt"  : quadSumVariationHisto(hnomi["pt"], hpdf["pt"], postfixName="pdf")
        }

        for var in ["eta", "pt"]:

            hWithPdfUnc = copy.deepcopy(hnomi[var].Clone(f"{hnomi[var].GetName()}_withPDFUnc"))
            setHistErrorFromHisto(hWithPdfUnc, pdfUnc[var])

            canvas.cd()
            hnomi[var].SetTitle(f"{p},   prefit,   {args.lumi} fb^{{-1}}")
            #hnomi[var].GetXaxis().SetTitleSize(0.05)
            #hnomi[var].GetXaxis().SetLabelSize(0.04)
            hnomi[var].GetXaxis().SetLabelSize(0)
            #hnomi[var].GetXaxis().SetTitle("Muon " + ("#eta" if var == "eta" else "p_{T} [GeV]"))
            hnomi[var].GetXaxis().SetTitle("")  
            hnomi[var].GetYaxis().SetTitleSize(0.05)
            hnomi[var].GetYaxis().SetLabelSize(0.04)
            hnomi[var].GetYaxis().SetTitleOffset(1.45)
            hnomi[var].GetYaxis().SetTitle("Events")

            ########
            hnomi[var].Draw("HE")
            hnomi[var].SetLineColor(ROOT.kBlack)
            hnomi[var].SetMarkerColor(ROOT.kBlack)
            hnomi[var].SetMarkerStyle(20)
            hnomi[var].SetMarkerSize(1.2)
            hnomi[var].SetFillColor(0)
            # pdf
            hWithPdfUnc.SetFillColor(ROOT.kAzure+7)
            hWithPdfUnc.SetFillColorAlpha(ROOT.kAzure+7, 0.6)
            hWithPdfUnc.SetFillStyle(1001)
            hWithPdfUnc.Draw("E2 SAME")
            # alphaS
            halphaSUp[var].SetLineColor(ROOT.kPink-6)
            halphaSUp[var].SetLineWidth(2)
            halphaSUp[var].SetFillColor(0)
            halphaSUp[var].SetMarkerStyle(0)
            halphaSUp[var].Draw("HIST SAME")
            halphaSDown[var].SetLineColor(ROOT.kPink-4)
            halphaSDown[var].SetLineWidth(2)
            halphaSDown[var].SetFillColor(0)
            halphaSDown[var].SetMarkerStyle(0)
            halphaSDown[var].Draw("HIST SAME")
            # postfit yields in given
            if args.addPostfitShape:
                hpfs[var].SetMarkerStyle(26)
                hpfs[var].SetMarkerColor(ROOT.kBlue+1)
                hpfs[var].SetFillColor(0)
                hpfs[var].SetLineColor(ROOT.kBlue+1)
                hpfs[var].Draw("PE SAME")
                
            hnomi[var].Draw("HESAME")
            canvas.RedrawAxis("sameaxis")
            
            #postfix = f"_{args.postfix}" if args.postfix else ""
            #for ext in ["png","pdf"]:
            #    canvas.SaveAs(f"{outdir}{hnomi[var].GetName()}_theoryUnc{postfix}.{ext}")

            pad2 = ROOT.TPad("pad2","pad2", 0.0, 0.0, 1, 0.95)
            pad2.SetTickx(1)
            pad2.SetTicky(1)
            pad2.SetTopMargin(1-lowerPanelHeight)
            pad2.SetRightMargin(0.04)
            pad2.SetLeftMargin(0.15)
            pad2.SetBottomMargin(0.12)
            pad2.SetFillColor(0)
            pad2.SetGridy(1)
            pad2.SetFillStyle(0)
            pad2.Draw()
            pad2.cd()
            ratio = hWithPdfUnc.Clone("ratio")
            ratioAlphaSUp = halphaSUp[var].Clone("ratioAlphaSUp")
            ratioAlphaSDown = halphaSDown[var].Clone("ratioAlphaSDown")
            den_noerr = hnomi[var].Clone("den_noerr")
            den = hnomi[var].Clone("den")
            for iBin in range (1, den_noerr.GetNbinsX()+1):
                den_noerr.SetBinError(iBin, 0.)
            ratio.Divide(den_noerr)
            ratioAlphaSUp.Divide(den_noerr)
            ratioAlphaSDown.Divide(den_noerr)
            den.Divide(den_noerr)
            miny, maxy = getMinMaxMultiHisto([ratio, den, ratioAlphaSUp, ratioAlphaSDown])
            den.SetTitle("")
            den.Draw("HE")
            den.GetXaxis().SetTitleSize(0.05)
            den.GetXaxis().SetTitleOffset(1.05)
            den.GetXaxis().SetLabelSize(0.04)
            den.GetXaxis().SetTitle("Muon " + ("#eta" if var == "eta" else "p_{T} [GeV]"))
            den.GetYaxis().SetTitleSize(0.05)
            den.GetYaxis().SetLabelSize(0.04)
            den.GetYaxis().SetTitleOffset(1.45)
            den.GetYaxis().SetTitle("Ratio with nominal")
            
            diff = maxy - miny
            den.GetYaxis().SetRangeUser(miny - 0.1 * diff, maxy + 0.5 * diff)
            ratio.Draw("E2 SAME")
            ratioAlphaSUp.Draw("HIST SAME")
            ratioAlphaSDown.Draw("HIST SAME")
            den.Draw("HE SAME")

            ylowLegend = 0.40
            if args.addPostfitShape:
                ratioPostfit = hpfs[var].Clone("ratioPostfit")
                ratioPostfit.Divide(den_noerr)
                ratioPostfit.Draw("PE SAME")
                ylowLegend = 0.39
                
            leg = ROOT.TLegend(0.18, ylowLegend, 0.92, 0.54)
            leg.SetFillStyle(0)
            leg.SetFillColor(0)
            leg.SetFillColorAlpha(0,0.6)
            leg.SetBorderSize(0)
            #leg.SetTextFont(62)
            leg.SetNColumns(2)
            leg.AddEntry(den, "MC statistics", "PLE")
            leg.AddEntry(ratio, "PDFs","F")
            leg.AddEntry(ratioAlphaSUp, "#alpha_{S} up","L")
            leg.AddEntry(ratioAlphaSDown, "#alpha_{S} down","L")
            if args.addPostfitShape:
                leg.AddEntry(ratioPostfit, "postfit yields","PLE")
            leg.Draw("same")
            pad2.RedrawAxis("sameaxis")

            for ext in ["png","pdf"]:
                canvas.SaveAs(f"{outdir}{hnomi[var].GetName()}_theoryUnc{postfix}.{ext}")
            
        ########
        
        
    infile.Close()


                        
                                 
        
