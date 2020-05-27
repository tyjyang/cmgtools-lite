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

    nBinsYw = 10
    print "Now retrieving all the theory variations, to make the theory bands later on"
    histDict = utilities.getTheoryHistHelicityFast(xsecWithWptWeights=False,nBinsYw = nBinsYw)
    print "DONE, all the histograms for the theory bands were taken"

    # sanity check
    print "="*30
    print histDict["listkeys"]
    print "Found %d keys" % len(histDict["listkeys"])
    print "="*30

    outfilename = "plotHelicityChargeAsymmetry.root"
    tfOut = ROOT.TFile.Open(outname+outfilename,'recreate')
    tfOut.cd()
    
    lepton = "lep"
    canvas1D = ROOT.TCanvas("canvas1D","",1000,1000)

    hAsymPDF = {}
    hAsymTotTheory = {}

    hA4PDF = {}
    hA4TotTheory = {}

    hA0PDF = {}
    hA0TotTheory = {}

    hXsecTotTheory = {}
    hXsecNormTotTheory = {}
    hXsecPDF = {}
    hXsecNormPDF = {}

    for pol in ["left", "right", "unpolarized"]:
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
            hXsecTotTheory[(charge,pol)].SetName("hXsec_{ch}_{p}_TotTheory".format(ch=charge,p=pol))
            hXsecTotTheory[(charge,pol)].SetTitle("abs cross section TotTheory: %s %s" % (charge, pol))
            #
            hXsecPDF[(charge,pol)] = ROOT.TGraphAsymmErrors(nBinsYw)
            hXsecPDF[(charge,pol)].SetName("hXsec_{ch}_{p}_PDF".format(ch=charge,p=pol))
            hXsecPDF[(charge,pol)].SetTitle("abs cross section PDF: %s %s" % (charge, pol))
            #
            print "Make theory bands for absolute cross section and charge %s and pol %s" % (charge, pol)
            utilities.getTheoryBandHelicity(hXsecTotTheory[(charge,pol)], hXsecPDF[(charge,pol)], histDict["listkeys"], histDict["xsec"], nBinsYw, charge=charge, pol=pol)
            print "DONE"
            writeGraphIntoFile(hXsecTotTheory[(charge,pol)],tfOut)
            writeGraphIntoFile(hXsecPDF[(charge,pol)],tfOut)

            hXsecNormTotTheory[(charge,pol)] = ROOT.TGraphAsymmErrors(nBinsYw)
            hXsecNormTotTheory[(charge,pol)].SetName("hXsecNorm_{ch}_{p}_TotTheory".format(ch=charge,p=pol))
            hXsecNormTotTheory[(charge,pol)].SetTitle("norm cross section TotTheory: %s %s" % (charge, pol))
            #
            hXsecNormPDF[(charge,pol)] = ROOT.TGraphAsymmErrors(nBinsYw)
            hXsecNormPDF[(charge,pol)].SetName("hXsecNorm_{ch}_{p}_PDF".format(ch=charge,p=pol))
            hXsecNormPDF[(charge,pol)].SetTitle("norm cross section PDF: %s %s" % (charge, pol))

            print "Make theory bands for normalized cross section and charge %s and pol %s" % (charge, pol)
            utilities.getTheoryBandHelicity(hXsecNormTotTheory[(charge,pol)], hXsecNormPDF[(charge,pol)], histDict["listkeys"],histDict["xsecnorm"], nBinsYw, charge=charge, pol=pol)            
            print "DONE"
            writeGraphIntoFile(hXsecNormTotTheory[(charge,pol)],tfOut)
            writeGraphIntoFile(hXsecNormPDF[(charge,pol)],tfOut)

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

                print "Make theory bands for A4 and charge %s" % charge
                utilities.getTheoryBandHelicity(hA4TotTheory[charge], hA4PDF[charge], histDict["listkeys"], histDict["A4"], nBinsYw, charge=charge, pol=pol)
                print "DONE"
                writeGraphIntoFile(hA4TotTheory[charge],tfOut)
                writeGraphIntoFile(hA4PDF[charge],tfOut)

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


    tfOut.Close()
