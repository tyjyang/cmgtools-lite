#!/bin/env python

import os, sys, re, array, math

## safe batch mode                                 
args = sys.argv[:]
sys.argv = ['-b']
import ROOT
sys.argv = args
ROOT.gROOT.SetBatch(True)
ROOT.PyConfig.IgnoreCommandLineOptions = True

sys.path.append(os.getcwd() + "/plotUtils/")
from utility import *

import utilities
utilities = utilities.util()

def writeHistoIntoFile(h,f, name="", verbose=True):
    f.cd()
    h.Write(name if len(name) else h.GetName())
    if verbose: 
        print ">>> Saving histogram: {n}".format(n=name if len(name) else h.GetName())


def writeGraphIntoFile(g, f, name="", verbose=True):
    f.cd()
    g.Write(name if len(name) else g.GetName())
    if verbose: 
        print ">>> Saving graph: {n}".format(n=name if len(name) else g.GetName())


if __name__ == "__main__":


    from optparse import OptionParser
    parser = OptionParser(usage='%prog [options]')
    parser.add_option('-o','--outdir', dest='outdir', default='', type='string', help='output directory to save things')
    (options, args) = parser.parse_args()


    ROOT.TH1.SetDefaultSumw2()
    ROOT.TH1.StatOverflows(True)
    
    outname = options.outdir
    createPlotDirAndCopyPhp(outname)

    nBinsYw = 12
    print "Now retrieving all the theory variations, to make the theory bands later on"
    scale =  0.5 / 35900.0 # factor 0.5 for the double channel , and then dividing by lumi to go from yields to pb
    # there might be a factor 4 = 1/0.25 because dividing by bin width, but it will be done in plotYW directly (only relevant for aboslute xsec)
    histDict = utilities.getTheoryHistHelicityFast(xsecWithWptWeights=False,nBinsYw = nBinsYw,scale=scale)
    print "DONE, all the histograms for the theory bands were taken"

    # sanity check
    print "="*30
    print histDict["listkeys"]
    print "Found %d keys" % len(histDict["listkeys"])
    print "="*30

    pseudodata = ROOT.TH1D("pseudodata","",nBinsYw,0,0.25*nBinsYw) # 0.25 width, we only care about first 10 bins

    outfilename = "plotHelicityChargeAsymmetry.root"
    tfOut = ROOT.TFile.Open(outname+outfilename,'recreate')
    tfOut.cd()
    
    lepton = "lep"
    canvas1D = ROOT.TCanvas("canvas1D","",1000,1000)

    hAsymPDF = {}
    hAsymTotTheory = {}

    hA4PDF = {}
    hA4TotTheory = {}

    hA4PDFonly = {}
    hA4alphaSonly = {}
    hA4QCDonly = {}

    hA0PDF = {}
    hA0TotTheory = {}

    hXsecTotTheory = {}
    hXsecNormTotTheory = {}
    hXsecPDF = {}
    hXsecNormPDF = {}

    for pol in ["left", "right", "long", "unpolarized"]:
        hAsymPDF[pol] = ROOT.TGraphAsymmErrors(nBinsYw)
        hAsymPDF[pol].SetName("hAsymPDF_{pol}".format(pol=pol))
        hAsymPDF[pol].SetTitle("Charge asymmetry PDF: %s" % pol)
        #
        hAsymTotTheory[pol] = ROOT.TGraphAsymmErrors(nBinsYw)
        hAsymTotTheory[pol].SetName("hAsymTotTheory_{pol}".format(pol=pol))
        hAsymTotTheory[pol].SetTitle("Charge asymmetry TotTheory: %s" % pol)

        print "Make theory bands for asymmetry and pol %s" % pol
        utilities.getTheoryBandHelicity(hAsymTotTheory[pol], hAsymPDF[pol], histDict["listkeys"], histDict["asym"], nBinsYw, charge="all", pol=pol)
        print "DONE"
        writeGraphIntoFile(hAsymTotTheory[pol],tfOut)
        writeGraphIntoFile(hAsymPDF[pol],tfOut)

        # pseudodata = ROOT.TH1D("pseudodata","",10,0,2.5)
        # for ib in range(1,1+pseudodata.GetNbinsX()):
        #     xval = ROOT.Double(0)
        #     yval = ROOT.Double(0)
        #     hAsymPDF.GetPoint(ib-1,xval,yval)
        #     pseudodata.SetBinContent(ib,yval)
        #     pseudodata.SetBinError(ib,hAsymPDF.GetErrorY(ib-1))

        #     drawXsecAndTheoryband(pseudodata, hAsymTotTheory, "|y_{W}|", "charge asymmetry::-0.10,0.40","chargeAsym_unpol_yW_lep_dataAndExp",outname,labelRatioTmp="data-exp::-0.10,0.10",draw_both0_noLog1_onlyLog2=1,passCanvas=canvas1D,lumi="35.9",
        #                           lowerPanelHeight=0.35, moreTextLatex="unpolarized", 
        #                           invertRatio=False, useDifferenceInLowerPanel=True,
        #                           histMCpartialUnc=hAsymPDF, histMCpartialUncLegEntry="PDFs #oplus #alpha_{S}",skipPreliminary=True)

        for charge in ["plus", "minus"]:

            hXsecTotTheory[(charge,pol)] = ROOT.TGraphAsymmErrors(nBinsYw)
            hXsecTotTheory[(charge,pol)].SetName("hXsecTotTheory_{ch}_{p}".format(ch=charge,p=pol))
            hXsecTotTheory[(charge,pol)].SetTitle("abs cross section TotTheory: %s %s" % (charge, pol))
            #
            hXsecPDF[(charge,pol)] = ROOT.TGraphAsymmErrors(nBinsYw)
            hXsecPDF[(charge,pol)].SetName("hXsecPDF_{ch}_{p}".format(ch=charge,p=pol))
            hXsecPDF[(charge,pol)].SetTitle("abs cross section PDF: %s %s" % (charge, pol))
            #
            print "Make theory bands for absolute cross section and charge %s and pol %s" % (charge, pol)
            utilities.getTheoryBandHelicity(hXsecTotTheory[(charge,pol)], hXsecPDF[(charge,pol)], histDict["listkeys"], histDict["xsec"], nBinsYw, charge=charge, pol=pol)
            print "DONE"
            writeGraphIntoFile(hXsecTotTheory[(charge,pol)],tfOut)
            writeGraphIntoFile(hXsecPDF[(charge,pol)],tfOut)

            for ib in range(1,1+pseudodata.GetNbinsX()):
                xval = ROOT.Double(0)
                yval = ROOT.Double(0)
                hXsecPDF[(charge,pol)].GetPoint(ib-1,xval,yval)
                pseudodata.SetBinContent(ib,yval)
                pseudodata.SetBinError(ib,hXsecPDF[(charge,pol)].GetErrorY(ib-1))

            # range should be ::-200,3500
            drawXsecAndTheoryband(pseudodata, hXsecTotTheory[(charge,pol)], "|y_{W}|::0,2.5", "abs xsec","absXsec_{c}_{p}_yW_lep_dataAndExp".format(c=charge,p=pol),outname,labelRatioTmp="data/exp::0.8,1.2",draw_both0_noLog1_onlyLog2=1,passCanvas=canvas1D,lumi="35.9",
                                  lowerPanelHeight=0.35, moreTextLatex="{c}_{p}".format(c=charge,p=pol), 
                                  invertRatio=False, useDifferenceInLowerPanel=False,
                                  histMCpartialUnc=hXsecPDF[(charge,pol)], histMCpartialUncLegEntry="PDFs #oplus #alpha_{S}",skipPreliminary=True)


            hXsecNormTotTheory[(charge,pol)] = ROOT.TGraphAsymmErrors(nBinsYw)
            hXsecNormTotTheory[(charge,pol)].SetName("hXsecNormTotTheory_{ch}_{p}".format(ch=charge,p=pol))
            hXsecNormTotTheory[(charge,pol)].SetTitle("norm cross section TotTheory: %s %s" % (charge, pol))
            #
            hXsecNormPDF[(charge,pol)] = ROOT.TGraphAsymmErrors(nBinsYw)
            hXsecNormPDF[(charge,pol)].SetName("hXsecNormPDF_{ch}_{p}".format(ch=charge,p=pol))
            hXsecNormPDF[(charge,pol)].SetTitle("norm cross section PDF: %s %s" % (charge, pol))

            print "Make theory bands for normalized cross section and charge %s and pol %s" % (charge, pol)
            utilities.getTheoryBandHelicity(hXsecNormTotTheory[(charge,pol)], hXsecNormPDF[(charge,pol)], histDict["listkeys"],histDict["xsecnorm"], nBinsYw, charge=charge, pol=pol)            
            print "DONE"
            writeGraphIntoFile(hXsecNormTotTheory[(charge,pol)],tfOut)
            writeGraphIntoFile(hXsecNormPDF[(charge,pol)],tfOut)

            for ib in range(1,1+pseudodata.GetNbinsX()):
                xval = ROOT.Double(0)
                yval = ROOT.Double(0)
                hXsecNormPDF[(charge,pol)].GetPoint(ib-1,xval,yval)
                pseudodata.SetBinContent(ib,yval)
                pseudodata.SetBinError(ib,hXsecNormPDF[(charge,pol)].GetErrorY(ib-1))

            # range should be ::-200,3500
            drawXsecAndTheoryband(pseudodata, hXsecNormTotTheory[(charge,pol)], "|y_{W}|::0,2.5", "norm xsec","normXsec_{c}_{p}_yW_lep_dataAndExp".format(c=charge,p=pol),outname,labelRatioTmp="data/exp::0.8,1.2",draw_both0_noLog1_onlyLog2=1,passCanvas=canvas1D,lumi="35.9",
                                  lowerPanelHeight=0.35, moreTextLatex="{c}_{p}".format(c=charge,p=pol), 
                                  invertRatio=False, useDifferenceInLowerPanel=False,
                                  histMCpartialUnc=hXsecNormPDF[(charge,pol)], histMCpartialUncLegEntry="PDFs #oplus #alpha_{S}",skipPreliminary=True)

            # A4 and A0
            if pol == "unpolarized":
                # A4
                hA4PDF[charge] = ROOT.TGraphAsymmErrors(nBinsYw)
                hA4PDF[charge].SetName("hA4PDF_{charge}".format(charge=charge))
                hA4PDF[charge].SetTitle("A4 PDF: %s" % charge)
                #
                hA4TotTheory[charge] = ROOT.TGraphAsymmErrors(nBinsYw)
                hA4TotTheory[charge].SetName("hA4TotTheory_{charge}".format(charge=charge))
                hA4TotTheory[charge].SetTitle("A4 TotTheory: %s" % charge)

                hA4PDFonly[charge] = ROOT.TH1F("hA4PDFonly_{ch}".format(ch=charge),
                                               "A4 PDF only: {ch}".format(ch=charge),
                                                nBinsYw,0,0.25*nBinsYw)
                hA4alphaSonly[charge] = ROOT.TGraphAsymmErrors(nBinsYw)
                hA4alphaSonly[charge].SetName("hA4alphaSonly_{charge}".format(charge=charge))
                hA4alphaSonly[charge].SetTitle("A4 alphaS only: %s" % charge)

                hA4QCDonly[charge] = ROOT.TGraphAsymmErrors(nBinsYw)
                hA4QCDonly[charge].SetName("hA4QCDonly_{charge}".format(charge=charge))
                hA4QCDonly[charge].SetTitle("A4 QCD only: %s" % charge)

                print "Make theory bands for A4 and charge %s" % charge
                utilities.getTheoryBandHelicity(hA4TotTheory[charge], hA4PDF[charge], histDict["listkeys"], histDict["A4"], nBinsYw, charge=charge, pol=pol)
                print "DONE"
                writeGraphIntoFile(hA4TotTheory[charge],tfOut)
                writeGraphIntoFile(hA4PDF[charge],tfOut)

                utilities.checkTheoryBandHelicity(hA4PDFonly[charge], hA4alphaSonly[charge], hA4QCDonly[charge],histDict["listkeys"], histDict["A4"], nBinsYw, charge=charge,pol=pol)

                for ib in range(1,1+pseudodata.GetNbinsX()):
                    xval = ROOT.Double(0)
                    yval = ROOT.Double(0)
                    hA4PDF[charge].GetPoint(ib-1,xval,yval)
                    pseudodata.SetBinContent(ib,yval)
                    pseudodata.SetBinError(ib,hA4PDF[charge].GetErrorY(ib-1))

                # range should be ::-200,3500
                drawXsecAndTheoryband(pseudodata, hA4TotTheory[charge], "|y_{W}|::0,2.5", "A4::-1.0,2.0","A4_{c}_{p}_yW_lep_dataAndExp".format(c=charge,p=pol),outname,labelRatioTmp="data/exp::-0.1,0.1",draw_both0_noLog1_onlyLog2=1,passCanvas=canvas1D,lumi="35.9",lowerPanelHeight=0.35, moreTextLatex="{c}_{p}".format(c=charge,p=pol), 
                                      invertRatio=False, useDifferenceInLowerPanel=True,
                                      histMCpartialUnc=hA4PDF[charge], histMCpartialUncLegEntry="PDFs #oplus #alpha_{S}",skipPreliminary=True)
                
                drawCheckTheoryBand(hA4PDFonly[charge], hA4alphaSonly[charge], hA4QCDonly[charge],"|y_W|::0,2.5","rel. unc.::0.8,1.2", "theoryBands_A4_{ch}".format(ch=charge),outname,draw_both0_noLog1_onlyLog2=1,lumi="35.9",moreTextLatex="{c}_{p}".format(c=charge,p=pol))

                # A0
                hA0PDF[charge] = ROOT.TGraphAsymmErrors(nBinsYw)
                hA0PDF[charge].SetName("hA0PDF_{charge}".format(charge=charge))
                hA0PDF[charge].SetTitle("A0 PDF: %s" % charge)
                #
                hA0TotTheory[charge] = ROOT.TGraphAsymmErrors(nBinsYw)
                hA0TotTheory[charge].SetName("hA0TotTheory_{charge}".format(charge=charge))
                hA0TotTheory[charge].SetTitle("A0 TotTheory: %s" % charge)

                print "Make theory bands for A0 and charge %s" % charge
                utilities.getTheoryBandHelicity(hA0TotTheory[charge], hA0PDF[charge], histDict["listkeys"], histDict["A0"], nBinsYw, charge=charge, pol=pol)
                print "DONE"
                writeGraphIntoFile(hA0TotTheory[charge],tfOut)
                writeGraphIntoFile(hA0PDF[charge],tfOut)


                for ib in range(1,1+pseudodata.GetNbinsX()):
                    xval = ROOT.Double(0)
                    yval = ROOT.Double(0)
                    hA0PDF[charge].GetPoint(ib-1,xval,yval)
                    pseudodata.SetBinContent(ib,yval)
                    pseudodata.SetBinError(ib,hA0PDF[charge].GetErrorY(ib-1))

                # range should be ::-200,3500
                drawXsecAndTheoryband(pseudodata, hA0TotTheory[charge], "|y_{W}|::0,2.5", "A0::0.05,0.2","A0_{c}_{p}_yW_lep_dataAndExp".format(c=charge,p=pol),outname,labelRatioTmp="data/exp::-0.1,0.1",draw_both0_noLog1_onlyLog2=1,passCanvas=canvas1D,lumi="35.9",
                                      lowerPanelHeight=0.35, moreTextLatex="{c}_{p}".format(c=charge,p=pol), 
                                    invertRatio=False, useDifferenceInLowerPanel=True,
                                      histMCpartialUnc=hA0PDF[charge], histMCpartialUncLegEntry="PDFs #oplus #alpha_{S}",skipPreliminary=True)

    tfOut.Close()
