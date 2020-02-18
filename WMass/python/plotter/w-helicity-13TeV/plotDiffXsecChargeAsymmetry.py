#!/bin/env python

# python w-helicity-13TeV/plotDiffXsecChargeAsymmetry.py -i cards/diffXsec_mu_2018_07_11_group10_coarseBin/ -o plots/diffXsec/chargeAsymmetry/muon/ -c mu -t toys/diffXsec_mu_2018_07_11_group10_coarseBin/toys_comb_WchargeAsymmetry.root -n [--no-group-POI --hessian -s <suffix>]

import ROOT, os, sys, re, array, math

from make_diff_xsec_cards import getDiffXsecBinning
from make_diff_xsec_cards import templateBinning

sys.path.append(os.getcwd() + "/plotUtils/")
from utility import *

ROOT.gROOT.SetBatch(True)

import utilities
utilities = utilities.util()

def scaleTH2inRange(h2, nPtMinToScale=1, scale=0.5):
    nPt  = h2.GetNbinsY()
    nEta = h2.GetNbinsX()
    for ieta in range(1,nEta+1):
        for ipt in range(nPtMinToScale,nPt+1):
            h2.SetBinContent(ieta,ipt, scale*h2.GetBinContent(ieta,ipt))
            h2.SetBinError(ieta,ipt, scale*h2.GetBinError(ieta,ipt))


def writeHistoIntoFile(h,f, name="", verbose=True):
    f.cd()
    h.Write(name if len(name) else h.GetName())
    if verbose: 
        print ">>> Saving histogram: {n}".format(n=name if len(name) else h.GetName())


def normGraphByBinWidth1D(gr, h):

    for ib in range(h.GetNbinsX()):
        xval = h.GetBinCenter(ib+1)
        yval = gr.Eval(xval)
        binwidth =  h.GetBinWidth(ib+1)
        yerrhigh = gr.GetErrorYhigh(ib)
        yerrlow = gr.GetErrorYlow(ib)
        gr.SetPoint(ib, xval, yval/binwidth)
        gr.SetPointEYhigh(ib, yerrhigh/binwidth)
        gr.SetPointEYlow(ib, yerrlow/binwidth)

def normGraphByBinArea2D(gr, h):

    # the graph associated to the unrolled, to be made from the 2D (i.e. h), has bin wdth = 1, and centered on integers
    for ieta in range(h.GetNbinsX()):
        for ipt in range(h.GetNbinsY()):
            binarea = h.GetXaxis().GetBinWidth(ieta+1) * h.GetYaxis().GetBinWidth(ipt+1)
            ibin = ieta + ipt * h.GetNbinsX()
            xval = 1 + ibin
            yval = gr.Eval(xval)
            yerrhigh = gr.GetErrorYhigh(ibin)
            yerrlow = gr.GetErrorYlow(ibin)
            gr.SetPoint(ibin, xval, yval/binarea)
            gr.SetPointEYhigh(ibin, yerrhigh/binarea)
            gr.SetPointEYlow(ibin, yerrlow/binarea)

def scaleGraphByConstant1D(gr, c, scaleX=False):

    xvals = gr.GetX()
    yvals = gr.GetY()
    #print "scaleGraphByConstant1D: %s %s" % (gr.GetName(), ",".join(str(x) for x in  xvals))
    for ib in range(len(xvals)):
        xval = xvals[ib]
        yval = yvals[ib]
        if scaleX:
            xerrhigh = gr.GetErrorXhigh(ib)
            xerrlow = gr.GetErrorXlow(ib)
            gr.SetPoint(      ib, c * xval, yval)
            gr.SetPointEXhigh(ib, c * xerrhigh)
            gr.SetPointEXlow( ib, c * xerrlow)
        else:
            yerrhigh = gr.GetErrorYhigh(ib)
            yerrlow = gr.GetErrorYlow(ib)
            gr.SetPoint(ib, xval, c * yval)
            gr.SetPointEYhigh(ib, c * yerrhigh)
            gr.SetPointEYlow(ib,  c * yerrlow)


def getTH1fromTH2(h2D,h2Derr=None,unrollAlongX=True):  # unrollAlongX=True --> select rows, i.e. takes a stripe from x1 to xn at same y, then go to next stripe at next y
    nX = h2D.GetNbinsX()
    nY = h2D.GetNbinsY()
    nbins = nX * nY
    name = h2D.GetName() + "_unrollTo1D_" + ("eta" if unrollAlongX else "pt") 
    newh = ROOT.TH1D(name,h2D.GetTitle(),nbins,0.5,nbins+0.5)
    if 'TH2' not in h2D.ClassName(): raise RuntimeError, "Calling getTH1fromTH2 on something that is not TH2"
    for i in xrange(nX):
        for j in xrange(nY):
            if unrollAlongX:
                ibin = 1 + i + j * nX
            else:
                ibin = 1 + j + i * nY
            newh.SetBinContent(ibin,h2D.GetBinContent(i+1,j+1))
            if h2Derr: newh.SetBinError(ibin,h2Derr.GetBinContent(i+1,j+1))
            else:      newh.SetBinError(ibin,h2D.GetBinError(i+1,j+1))
    return newh


def getUnrolledGraph(gr, neta, npt, unrollAlongX=True):
    # basically return the same graph if need to unroll along eta (X), else rearrange the bins so to reproduce unrolled along pt
    name = gr.GetName() + "_unrollTo1D_" + ("eta" if unrollAlongX else "pt") 
    if unrollAlongX:         
        retgr = gr.Clone(name)
        return retgr
    else:
        retgr = gr.Clone(name)
        for ib in range(gr.GetN()):
            ieta = ib % neta
            ipt =  ib / neta
            ibnew = ipt + ieta * npt
            yval = gr.Eval(ibnew+1)
            yerrHigh = gr.GetErrorYhigh(ibnew)
            yerrLow = gr.GetErrorYlow(ibnew)
            retgr.SetPoint(ib, ib+1, yval)
            retgr.SetPointEYhigh(ib, yerrHigh)
            retgr.SetPointEYlow(ib, yerrLow)
        return retgr

if __name__ == "__main__":


    from optparse import OptionParser
    parser = OptionParser(usage='%prog [options]')
    parser.add_option('-i','--input', dest='inputdir', default='', type='string', help='input directory with all the cards inside. It is used to get other information')
    #parser.add_option(     '--no-group-POI', dest='noGroupPOI', default=False , action='store_true', help='Specify that _group_<N>_ is not present in name of POI')
    parser.add_option('-o','--outdir', dest='outdir', default='', type='string', help='output directory to save things')
    parser.add_option('-t','--toyfile', dest='toyfile', default='.', type='string', help='Root file with toys.')
    parser.add_option('-c','--channel', dest='channel', default='', type='string', help='name of the channel (mu or el). For combination, select lep')
    parser.add_option('-C','--charge', dest='charge', default='plus,minus', type='string', help='charges to run')
    parser.add_option('-s','--suffix', dest='suffix', default='', type='string', help='Suffix added to output dir (i.e, either to hessian or toys in the path)')
    parser.add_option('-f','--friend', dest='friend', default='', type='string', help='Root file with friend tree containing total xsec (it does not include the outliers). Tree name is assumed to be "toyFriend"')
    parser.add_option('-p','--palette', dest='palette', default='55', type='int', help='Palette for plots')
    parser.add_option('-l','--lumi-norm', dest='lumiNorm', default='-1', type='float', help='If > 0, divide cross section by this factor (lumi in 1/Pb)')
    parser.add_option(     '--lumiInt', dest='lumiInt', default='35.9', type='float', help='Integrated luminosity')
    parser.add_option('-n','--norm-width', dest='normWidth' , default=False , action='store_true',   help='Normalize cross section histograms dividing by bin width')
    parser.add_option(     '--hessian', dest='hessian' , default=False , action='store_true',   help='The file passed with -t is interpreted as hessian')
    parser.add_option(     '--fit-data', dest='fitData' , default=False , action='store_true',   help='If True, axis range in plots is customized for data')
    parser.add_option(     '--invert-ratio', dest='invertRatio' , default=False , action='store_true',   help='If True, make ratio as exp./obs. (only when plotting data)')
    parser.add_option('-e','--expected-toyfile', dest='exptoyfile', default='', type='string', help='Root file to get expected and make ratio with data (only work with option --fit-data). If SAME, use same file as data and take _gen variables to get the expected')
    parser.add_option(       '--pt-range-bkg', dest='pt_range_bkg', action="append", type="float", nargs=2, default=[], help='Bins with gen level pt in this range are treated as background in the datacard, so there is no POI for them. Takes two float as arguments (increasing order) and can specify multiple times. They should match bin edges and a bin is not considered as background if at least one edge is outside this range (therefore, one can also choose a slightly larger range to avoid floating precision issues).')
    parser.add_option(       '--eta-range-bkg', dest='eta_range_bkg', action="append", type="float", nargs=2, default=[], help='Bins with gen level pt in this range are treated as background in the datacard, so there is no POI for them. Takes two float as arguments (increasing order) and can specify multiple times. They should match bin edges and a bin is not considered as background if at least one edge is outside this range (therefore, one can also choose a slightly larger range to avoid floating precision issues).')
    parser.add_option(     '--pt-range', dest='ptRange', default='template', type='string', help='Comma separated pair of floats used to define the pt range. If "template" is passed, the template min and max values are used.')
    parser.add_option(     '--blind-data', dest='blindData' , default=False , action='store_true',   help='If True, data is substituded with hessian (it requires passing the file with option --expected-toyfile')
    parser.add_option(       '--use-xsec-wpt', dest='useXsecWithWptWeights' , default=False, action='store_true', help='Use xsec file made with W-pt weights')
    parser.add_option(       '--force-allptbins-theoryband', dest='forceAllPtBinsTheoryBand' , default=False, action='store_true', help='Use all pt bins for theory band, even if some of them were treated as background')
    parser.add_option(      '--combineElePt01asBkg', dest='combineElePt01asBkg', default=False, action='store_true', help='If True, flavour combination is made by keeping some electron bins as background (input to be managed slightly diffferently)')
    parser.add_option(      '--skipPreliminary', dest='skipPreliminary', default=False, action='store_true', help='Do not add "Preliminary" to text on top left of canvas')
    parser.add_option('--n-split-unrolled-bins', dest='nSplitUnrolledBins', default='1', type='int', help='In how many parts the unrolled plots should be split into (e.g. 2 means two plots)')
    (options, args) = parser.parse_args()

    ROOT.TH1.SetDefaultSumw2()
    ROOT.TH1.StatOverflows(True)
    
    channel = options.channel
    if channel not in ["el","mu","lep"]:
        print "Error: unknown channel %s (select 'el' or 'mu' or 'lep')" % channel
        quit()

    isMuElComb = False
    if channel == "lep": isMuElComb = True

    charges = [x for x in options.charge.split(',')]
    for c in charges:
        if c not in ["plus", "minus"]:
            print "Error: unknown charge %s (select 'plus' or 'minus' or both separated by comma)" % c
            quit()
    hasBothCharges = True
    if len(charges) < 2:
        hasBothCharges = False

    if options.blindData:
        if len(options.exptoyfile) == 0:
            print "WARNING: you used option --blind-data, but you didn't specify a hessian file with option --expected-toyfile. Exit"
            quit()
        if options.fitData:
            print "WARNING: options --blind-data and --fit-data are not compatible. The former will effectively make plots with Asimov. Exit"
            quit()            


    if not options.toyfile:
        print "Error: you should specify a file containing the toys using option -t <name>. Exit"
        quit()



    if options.inputdir:
        if not options.inputdir.endswith("/"): options.inputdir += "/"
        basedir = os.path.basename(os.path.normpath(options.inputdir))
        binfile = options.inputdir + "binningPtEta.txt"
    else:
        print "Error: you should specify an input folder containing all the cards using option -i <name>. Exit"
        quit()


    if options.outdir:
        outname = options.outdir
        addStringToEnd(outname,"/",notAddIfEndswithMatch=True)
        if options.hessian:     outname = outname + "/hessian"
        else              :     outname = outname + "/toys"
        if len(options.suffix): outname = outname + "_" + options.suffix
        if options.blindData:   outname = outname + "_blindData"
        outname = outname + "/"
        createPlotDirAndCopyPhp(outname)
    else:
        print "Error: you should specify an output folder using option -o <name>. Exit"
        quit()

    xsecChargefile = {}
    for charge in charges:
        xsecChargefile[charge] = options.inputdir + "W{lep}_{ch}_shapes_xsec.root".format(lep=channel, ch=charge)

    etaPtBinningVec = getDiffXsecBinning(binfile, "gen")
    genBins  = templateBinning(etaPtBinningVec[0],etaPtBinningVec[1])
    netabins = genBins.Neta
    nptbins  = genBins.Npt

    ptBinIsBackground = []  # will store a bool to assess whether the given ipt index is considered as background
    nPtBinsBkg = 0
    for bin in range(genBins.Npt):
        ptBinIsBackground.append(False)
    ptRangesBkg = options.pt_range_bkg
    if len(ptRangesBkg):
        hasPtRangeBkg = True
        print "Signal bins with gen pt in the following ranges were considered as background processes, so there is no POI for them"
        print options.pt_range_bkg            
        for index in range(genBins.Npt):
            for pair in ptRangesBkg:
        #print pair
                if genBins.ptBins[index] >= pair[0] and genBins.ptBins[index+1] <= pair[1]:
                    ptBinIsBackground[index] = True
                    nPtBinsBkg += 1
    else:
        hasPtRangeBkg = False

    firstPtSignalBin = 0
    for i,val in enumerate(ptBinIsBackground):
        if not val: 
            firstPtSignalBin = i 
            break

    etaBinIsBackground = []  # will store a bool to assess whether the given ieta index is considered as background
    nEtaBinsBkg = 0
    for bin in range(genBins.Neta):
        etaBinIsBackground.append(False)
    etaRangesBkg = options.eta_range_bkg
    if len(etaRangesBkg):
        hasEtaRangeBkg = True
        print "Signal bins with gen eta in the following ranges were considered as background processes, so there is no POI for them"
        print options.eta_range_bkg            
        for index in range(genBins.Neta):
            for pair in etaRangesBkg:
        #print pair
                if genBins.etaBins[index] >= pair[0] and genBins.etaBins[index+1] <= pair[1]:
                    etaBinIsBackground[index] = True
                    nEtaBinsBkg += 1
    else:
        hasEtaRangeBkg = False


    # this file will store all histograms for later usage
    outfilename = "plotDiffXsecChargeAsymmetry.root"
    tfOut = ROOT.TFile.Open(outname+outfilename,'recreate')
    tfOut.cd()

    lepton = "electron" if channel == "el" else "muon" if channel == "mu" else "lepton"
    Wchannel = "W #rightarrow %s#nu" % ("e" if channel == "el" else "#mu" if channel == "mu" else "l")
    ptRangeText = "p_{{T}}^{{{l}}} #in [{ptmin:2g}, {ptmax:3g}] GeV".format(l="e" if channel == "el" else "#mu" if channel == "mu" else "l", 
                                                                            ptmin=genBins.ptBins[firstPtSignalBin], ptmax=genBins.ptBins[-1])
    #etaRangeText = "|#eta|^{{{l}}} #in [{etamin:3g}, {etamax:3g}] GeV".format(l="e" if channel == "el" else "#mu" if channel == "mu" else "l", 
    #                                                                          etamin=genBins.etaBins[0], etamax=genBins.etaBins[-1])
    etaRangeText = "|#eta|^{{{l}}} < {etamax:3g}".format(l="e" if channel == "el" else "#mu" if channel == "mu" else "l", 
                                                                              etamin=genBins.etaBins[0], etamax=genBins.etaBins[-1])
    labelRatioDataExp = "exp./obs.::0.9,1.1" if options.invertRatio else "obs./exp.::0.9,1.1"

    hChAsymm = ROOT.TH2F("hChAsymm_{lep}".format(lep=lepton),"Charge asymmetry: {Wch}".format(Wch=Wchannel),
                         genBins.Neta, array('d',genBins.etaBins), genBins.Npt, array('d',genBins.ptBins))
    hChAsymmErr = ROOT.TH2F("hChAsymmErr_{lep}".format(lep=lepton),"Charge asymmetry uncertainty: {Wch}".format(Wch=Wchannel),
                            genBins.Neta, array('d',genBins.etaBins), genBins.Npt, array('d',genBins.ptBins))
    

    hChAsymm1Deta = ROOT.TH1F("hChAsymm1Deta_{lep}".format(lep=lepton),"Charge asymmetry: {Wch}".format(Wch=Wchannel),
                              genBins.Neta, array('d',genBins.etaBins))
    hChAsymm1Dpt = ROOT.TH1F("hChAsymm1Dpt_{lep}".format(lep=lepton),"Charge asymmetry: {Wch}".format(Wch=Wchannel),
                              genBins.Npt, array('d',genBins.ptBins))


    #hChargeAsymPDF = ROOT.TH2F("hChargeAsymPDF_{lep}".format(lep=lepton),"Charge asymmetry: {Wch}".format(Wch=Wchannel),
    #                           genBins.Neta, array('d',genBins.etaBins), genBins.Npt, array('d',genBins.ptBins))
    # hChargeAsymPDF_1Deta = ROOT.TH1F("hChargeAsymPDF_1Deta_{lep}".format(lep=lepton),"Charge asymmetry: {Wch}".format(Wch=Wchannel),
    #                                  genBins.Neta, array('d',genBins.etaBins))
    # hChargeAsymPDF_1Dpt = ROOT.TH1F("hChargeAsymPDF_1Dpt_{lep}".format(lep=lepton),"Charge asymmetry: {Wch}".format(Wch=Wchannel),
    #                                 genBins.Npt, array('d',genBins.ptBins))
    #hChargeAsymTotTheory = ROOT.TH2F("hChargeAsymTotTheory_{lep}".format(lep=lepton),"Charge asymmetry: {Wch}".format(Wch=Wchannel),
    #                           genBins.Neta, array('d',genBins.etaBins), genBins.Npt, array('d',genBins.ptBins))
    # hChargeAsymTotTheory_1Deta = ROOT.TH1F("hChargeAsymTotTheory_1Deta_{lep}".format(lep=lepton),"Charge asymmetry: {Wch}".format(Wch=Wchannel),
    #                                  genBins.Neta, array('d',genBins.etaBins))
    # hChargeAsymTotTheory_1Dpt = ROOT.TH1F("hChargeAsymTotTheory_1Dpt_{lep}".format(lep=lepton),"Charge asymmetry: {Wch}".format(Wch=Wchannel),
    #                                 genBins.Npt, array('d',genBins.ptBins))
    hChargeAsymPDF = ROOT.TGraphAsymmErrors(genBins.Neta*genBins.Npt)
    hChargeAsymPDF.SetName("hChargeAsymPDF_{lep}".format(lep=lepton))
    hChargeAsymPDF.SetTitle("Charge asymmetry: {Wch}".format(Wch=Wchannel))
    #
    hChargeAsymTotTheory = ROOT.TGraphAsymmErrors(genBins.Neta*genBins.Npt)
    hChargeAsymTotTheory.SetName("hChargeAsymTotTheory_{lep}".format(lep=lepton))
    hChargeAsymTotTheory.SetTitle("Charge asymmetry: {Wch}".format(Wch=Wchannel))
    #
    hChargeAsymPDF_1Deta = ROOT.TGraphAsymmErrors(genBins.Neta)
    hChargeAsymPDF_1Deta.SetName("hChargeAsymPDF_1Deta_{lep}".format(lep=lepton))
    hChargeAsymPDF_1Deta.SetTitle("Charge asymmetry: {Wch}".format(Wch=Wchannel))
    #
    hChargeAsymPDF_1Dpt = ROOT.TGraphAsymmErrors(genBins.Npt)
    hChargeAsymPDF_1Dpt.SetName("hChargeAsymPDF_1Dpt_{lep}".format(lep=lepton))
    hChargeAsymPDF_1Dpt.SetTitle("Charge asymmetry: {Wch}".format(Wch=Wchannel))
    #
    hChargeAsymTotTheory_1Deta = ROOT.TGraphAsymmErrors(genBins.Neta)
    hChargeAsymTotTheory_1Deta.SetName("hChargeAsymTotTheory_1Deta_{lep}".format(lep=lepton))
    hChargeAsymTotTheory_1Deta.SetTitle("Charge asymmetry: {Wch}".format(Wch=Wchannel))
    #
    hChargeAsymTotTheory_1Dpt = ROOT.TGraphAsymmErrors(genBins.Npt)
    hChargeAsymTotTheory_1Dpt.SetName("hChargeAsymTotTheory_1Dpt_{lep}".format(lep=lepton))
    hChargeAsymTotTheory_1Dpt.SetTitle("Charge asymmetry: {Wch}".format(Wch=Wchannel))

    hChAsymm.SetDirectory(tfOut)
    hChAsymmErr.SetDirectory(tfOut)
    hChAsymm1Deta.SetDirectory(tfOut)
    hChAsymm1Dpt.SetDirectory(tfOut)
    #hChargeAsymPDF.SetDirectory(tfOut)
    #hChargeAsymPDF_1Deta.SetDirectory(tfOut)
    #hChargeAsymPDF_1Dpt.SetDirectory(tfOut)
    #hChargeAsymTotTheory.SetDirectory(tfOut)
    #hChargeAsymTotTheory_1Deta.SetDirectory(tfOut)
    #hChargeAsymTotTheory_1Dpt.SetDirectory(tfOut)

    # filled later on
    hDiffXsecPDF = {}
    hDiffXsecPDF_1Deta = {}
    hDiffXsecPDF_1Dpt = {}
    hDiffXsecTotTheory = {}
    hDiffXsecTotTheory_1Deta = {}
    hDiffXsecTotTheory_1Dpt = {}
    #
    hDiffXsecNormPDF = {}
    hDiffXsecNormPDF_1Deta = {}
    hDiffXsecNormPDF_1Dpt = {}
    hDiffXsecNormTotTheory = {}
    hDiffXsecNormTotTheory_1Deta = {}
    hDiffXsecNormTotTheory_1Dpt = {}
    # for tests
    hDiffXsecPDFonly_1Deta = {}
    hDiffXsecPDFonly_1Dpt = {}
    hDiffXsecAlphaonly_1Deta = {}
    hDiffXsecAlphaonly_1Dpt = {}
    hDiffXsecQCDonly_1Deta = {}
    hDiffXsecQCDonly_1Dpt = {}

    hDiffXsecNormPDFonly_1Deta = {}
    hDiffXsecNormPDFonly_1Dpt = {}
    hDiffXsecNormAlphaonly_1Deta = {}
    hDiffXsecNormAlphaonly_1Dpt = {}
    hDiffXsecNormQCDonly_1Deta = {}
    hDiffXsecNormQCDonly_1Dpt = {}

    ptminTheoryHist = options.ptRange.split(",")[0] if options.ptRange != "template" else -1.0
    if options.forceAllPtBinsTheoryBand: ptminTheoryHist = -1.0
    print "Now retrieving all the theory variations, to make the theory bands later on"
    histDict = utilities.getTheoryHistDiffXsecFast(options.useXsecWithWptWeights, 
                                                   ptmin=ptminTheoryHist)
    print "DONE, all the histograms for the theory bands were taken"

    # sanity check
    print "="*30
    print histDict["listkeys"]
    print "Found %d keys" % len(histDict["listkeys"])
    print "="*30
    # print histDict["xsec"][1]
    # print "="*30
    # print histDict["xsecnorm"][1]
    # print "="*30
    # print histDict["asym"][1]
    # print "="*30
    # quit()
    # returned object is
    # ret = { "xsec"     : [htheory, htheory_1Deta, htheory_1Dpt],
    #         "xsecnorm" : [htheory_xsecnorm, htheory_1Deta_xsecnorm, htheory_1Dpt_xsecnorm],
    #         "asym"     : [htheory_asym, htheory_1Deta_asym, htheory_1Dpt_asym],
    #         "listkeys" : [key[1] for key in htheory]}
    #
    # where listkeys is pdf1, pdf2, alphaSUp, muRUp, muFDown, ...
    # the rest are the lists of all theory histograms (and nominal) for 2D, and 1D projections, for xsec, xsecnorm and asym
    # they are dictionary with key (charge,theory_nuis), where theory_nuis is returned by listkeys, and charge is all for asym and plus/minus for the xsecs

    if hasBothCharges:
        print "Make theory bands for asymmetry"
        utilities.getTheoryBandDiffXsec(hChargeAsymTotTheory, hChargeAsymPDF, histDict["listkeys"], histDict["asym"][0], 
                                        genBins.Neta, genBins.Npt, charge="all")
        utilities.getTheoryBandDiffXsec1Dproj(hChargeAsymTotTheory_1Deta, hChargeAsymPDF_1Deta, histDict["listkeys"], histDict["asym"][1],
                                              genBins.Neta, charge="all")
        utilities.getTheoryBandDiffXsec1Dproj(hChargeAsymTotTheory_1Dpt, hChargeAsymPDF_1Dpt, histDict["listkeys"], histDict["asym"][2], 
                                              genBins.Npt, charge="all")
        print "DONE"


    nbins = (genBins.Neta - nEtaBinsBkg) * (genBins.Npt - nPtBinsBkg)
    binCount = 1        

    mainfile = options.toyfile
    if options.blindData:        
        mainfile = options.exptoyfile

    f = ROOT.TFile(mainfile, 'read')
    tree = f.Get('fitresults')
    if options.friend != "":
        tree.AddFriend('toyFriend',options.friend)

    if hasBothCharges:

        print "Now filling histograms with charge asymmetry"

        for ieta in range(1,genBins.Neta+1):
            for ipt in range(1,genBins.Npt+1):
                if ptBinIsBackground[ipt-1]: continue
                if etaBinIsBackground[ieta-1]: continue
                #print "\rBin {num}/{tot}".format(num=binCount,tot=nbins),
                sys.stdout.write('Bin {num}/{tot}   \r'.format(num=binCount,tot=nbins))
                sys.stdout.flush()
                binCount += 1

                if options.hessian: 
                    #central = utilities.getDiffXsecAsymmetryFromHessian(ieta-1,ipt-1,options.toyfile)
                    central = utilities.getDiffXsecAsymmetryFromHessianFast(ieta-1,ipt-1,
                                                                            nHistBins=2000, minHist=0., maxHist=1.0, tree=tree)
                    error = utilities.getDiffXsecAsymmetryFromHessianFast(ieta-1,ipt-1,
                                                                          nHistBins=2000, minHist=0., maxHist=1.0, tree=tree, getErr=True)
                else:                
                    central,up,down = utilities.getDiffXsecAsymmetryFromToysFast(ieta-1,ipt-1,
                                                                                 nHistBins=2000, minHist=0., maxHist=1.0, tree=tree)
                    #central,up,down = utilities.getDiffXsecAsymmetryFromToys(ieta-1,ipt-1,options.toyfile)
                    error = up - central
                hChAsymm.SetBinError(ieta,ipt,error)         
                hChAsymmErr.SetBinContent(ieta,ipt,error)
                hChAsymm.SetBinContent(ieta,ipt,central)        

        for ieta in range(1,genBins.Neta+1):            
            if etaBinIsBackground[ieta-1]: continue
            central = utilities.getDiffXsecAsymmetry1DFromHessianFast(ieta-1,isIeta=True,nHistBins=1000, minHist=0., maxHist=1.0, tree=tree)
            error   = utilities.getDiffXsecAsymmetry1DFromHessianFast(ieta-1,isIeta=True,nHistBins=1000, minHist=0., maxHist=1.0, tree=tree, getErr=True)
            hChAsymm1Deta.SetBinContent(ieta,central)
            hChAsymm1Deta.SetBinError(ieta,error)

        for ipt in range(1,genBins.Npt+1):            
            if ptBinIsBackground[ipt-1]: continue
            central = utilities.getDiffXsecAsymmetry1DFromHessianFast(ipt-1,isIeta=False,nHistBins=1000, minHist=0., maxHist=1.0, tree=tree)
            error   = utilities.getDiffXsecAsymmetry1DFromHessianFast(ipt-1,isIeta=False,nHistBins=1000, minHist=0., maxHist=1.0, tree=tree, getErr=True)
            hChAsymm1Dpt.SetBinContent(ipt,central)
            hChAsymm1Dpt.SetBinError(ipt,error)

    setTDRStyle()

    canvas = ROOT.TCanvas("canvas","",800,700)
    canvas1D = ROOT.TCanvas("canvas1D","",1000,1000)

    xaxisTitle = 'dressed %s |#eta|' % lepton
    yaxisTitle = 'dressed %s p_{T} (GeV)' % lepton
    if options.ptRange != "template":
        ptmin,ptmax = options.ptRange.split(',')
        yaxisTitle = yaxisTitle + "::%s,%s" % (ptmin, ptmax)
    #zaxisTitle = "Asymmetry::%.3f,%.3f" % (hChAsymm.GetMinimum(),hChAsymm.GetMaximum())
    if options.fitData:
        zaxisTitle = "Asymmetry::0.05,0.35"
        if channel == "el": zaxisTitle = "Asymmetry::0.0,0.45"
        zaxisTitle = "Asymmetry::%.3f,%.3f" % (max(0,0.99*hChAsymm.GetBinContent(hChAsymm.GetMinimumBin())),
                                               min(0.30,hChAsymm.GetBinContent(hChAsymm.GetMaximumBin())))

    else:
        zaxisTitle = "Asymmetry::0.05,0.35"
        if channel == "el": zaxisTitle = "Asymmetry::0.05,0.35"
        zaxisTitle = "Asymmetry::%.3f,%.3f" % (0.99*hChAsymm.GetBinContent(hChAsymm.GetMinimumBin()),
                                               hChAsymm.GetBinContent(hChAsymm.GetMaximumBin()))

    if hasBothCharges:
        drawCorrelationPlot(hChAsymm,
                            xaxisTitle, yaxisTitle, zaxisTitle,
                            hChAsymm.GetName(),
                            "ForceTitle",outname,1,1,False,False,False,1, canvasSize="700,625",leftMargin=0.14,rightMargin=0.22,passCanvas=canvas,palette=options.palette)

        zaxisTitle = "Asymmetry uncertainty::%.3f,%.3f" % (max(0,0.9*hChAsymmErr.GetBinContent(hChAsymmErr.GetMinimumBin())),
                                                           min(0.3,hChAsymmErr.GetBinContent(hChAsymmErr.GetMaximumBin())))
        #zaxisTitle = "Asymmetry uncertainty::0,0.04" 
        drawCorrelationPlot(hChAsymmErr,
                            xaxisTitle, yaxisTitle, zaxisTitle,
                            hChAsymmErr.GetName(),
                            "ForceTitle",outname,1,1,False,False,False,1, canvasSize="700,625",leftMargin=0.14,rightMargin=0.22,passCanvas=canvas,palette=options.palette)

        hChAsymmRelErr = hChAsymmErr.Clone(hChAsymmErr.GetName().replace("AsymmErr","AsymmRelErr"))
        hChAsymmRelErr.Divide(hChAsymm)
        hChAsymmRelErr.SetTitle(hChAsymmErr.GetTitle().replace("uncertainty","rel. uncertainty"))
        hChAsymmRelErr.SetDirectory(tfOut)

        #zaxisTitle = "Asymmetry relative uncertainty::%.4f,%.4f" % (max(0,0.99*hChAsymmRelErr.GetBinContent(hChAsymmRelErr.GetMinimumBin())),
                                                                     # min(0.12,1.01*hChAsymmRelErr.GetBinContent(hChAsymmRelErr.GetMaximumBin()))
                                                                     #)
        #zaxisTitle = "Asymmetry relative uncertainty::0.01,0.3" 
        zaxisTitle = "Asymmetry relative uncertainty::%.4f,%.4f" % (max(0,hChAsymmRelErr.GetBinContent(hChAsymmRelErr.GetMinimumBin())),
                                                                    min(0.4,hChAsymmRelErr.GetBinContent(hChAsymmRelErr.GetMaximumBin()))) 
        drawCorrelationPlot(hChAsymmRelErr,
                            xaxisTitle, yaxisTitle, zaxisTitle,
                            hChAsymmRelErr.GetName(),
                            "ForceTitle",outname,1,1,False,False,False,1, canvasSize="700,625",leftMargin=0.14,rightMargin=0.22,passCanvas=canvas,palette=options.palette)


        #additionalText = "W #rightarrow {lep}#nu::0.2,0.5,0.4,0.6".format(lep="e" if channel == "el" else "#mu" if channel == "mu" else "l") # pass x1,y1,x2,y2
        #additionalText = "W #rightarrow {lep}#nu;{pttext}::0.2,0.6,0.5,0.7".format(lep="e" if channel == "el" else "#mu" if channel == "mu" else "l", pttext=ptRangeText) # pass x1,y1,x2,y2
        texCoord = "0.2,0.65"
        additionalText = "W #rightarrow {lep}#nu;{pttext}::{txc},0.08,0.04".format(lep="e" if channel == "el" else "#mu" if channel == "mu" else "l", 
                                                                                    pttext=ptRangeText,
                                                                                    txc=texCoord)
        legendCoords = "0.2,0.4,0.75,0.85"
        drawSingleTH1(hChAsymm1Deta,xaxisTitle,"charge asymmetry",
                      "chargeAsym1D_eta_{fl}".format(fl=channel),
                      outname,labelRatioTmp="Rel.Unc.::0.9,1.1",legendCoords=legendCoords, draw_both0_noLog1_onlyLog2=1, drawLineLowerPanel="",
                      passCanvas=canvas1D, lumi=options.lumiInt,
                      lowerPanelHeight=0.35, moreTextLatex=additionalText
                      )
        legendCoords = "0.2,0.4,0.45,0.55"
        texCoord = "0.6,0.85"
        additionalText = "W #rightarrow {lep}#nu;{etatext}::{txc},0.08,0.04".format(lep="e" if channel == "el" else "#mu" if channel == "mu" else "l", 
                                                                                   etatext=etaRangeText,
                                                                                   txc=texCoord)
        drawSingleTH1(hChAsymm1Dpt,yaxisTitle,"charge asymmetry",
                      "chargeAsym1D_pt_{fl}".format(fl=channel),
                      outname,labelRatioTmp="Rel.Unc.::0.9,1.1",legendCoords=legendCoords, draw_both0_noLog1_onlyLog2=1, drawLineLowerPanel="",
                      passCanvas=canvas1D, lumi=options.lumiInt,
                      lowerPanelHeight=0.35, moreTextLatex=additionalText
                      )

        writeHistoIntoFile(hChAsymm,tfOut)
        writeHistoIntoFile(hChAsymmErr,tfOut)
        writeHistoIntoFile(hChAsymmRelErr,tfOut)
        writeHistoIntoFile(hChAsymm1Deta,tfOut)
        writeHistoIntoFile(hChAsymm1Dpt,tfOut)

    icharge = 0

    # hMu = {}
    # hMuErr = {}
    # hDiffXsec = {}
    # hDiffXsecErr = {}
    # hDiffXsecNorm = {}
    # hDiffXsecNormErr = {}
    # hDiffXsec_1Deta = {}
    # hDiffXsecNorm_1Deta = {}


    for charge in charges:

        icharge += 1
        print ""
        xaxisTitle = 'dressed %s |#eta|' % lepton  # it will be modified later, so I have to restore it here


        chargeSign = "+" if charge == "plus" else "-"

        hMu= ROOT.TH2F("hMu_{lep}_{ch}".format(lep=lepton,ch=charge),
                              "signal strength: {Wch}".format(Wch=Wchannel.replace('W','W{chs}'.format(chs=chargeSign))),
                              genBins.Neta, array('d',genBins.etaBins), genBins.Npt, array('d',genBins.ptBins))
        hMuErr= ROOT.TH2F("hMuErr_{lep}_{ch}".format(lep=lepton,ch=charge),
                                 "signal strength uncertainty: {Wch}".format(Wch=Wchannel.replace('W','W{chs}'.format(chs=chargeSign))),
                                 genBins.Neta, array('d',genBins.etaBins), genBins.Npt, array('d',genBins.ptBins))

        hDiffXsec= ROOT.TH2F("hDiffXsec_{lep}_{ch}".format(lep=lepton,ch=charge),
                              "cross section: {Wch}".format(Wch=Wchannel.replace('W','W{chs}'.format(chs=chargeSign))),
                              genBins.Neta, array('d',genBins.etaBins), genBins.Npt, array('d',genBins.ptBins))
        hDiffXsecErr= ROOT.TH2F("hDiffXsecErr_{lep}_{ch}".format(lep=lepton,ch=charge),
                                 "cross section uncertainty: {Wch}".format(Wch=Wchannel.replace('W','W{chs}'.format(chs=chargeSign))),
                                 genBins.Neta, array('d',genBins.etaBins), genBins.Npt, array('d',genBins.ptBins))
        hDiffXsecNorm= ROOT.TH2F("hDiffXsecNorm_{lep}_{ch}".format(lep=lepton,ch=charge),
                                  "normalized cross section: {Wch}".format(Wch=Wchannel.replace('W','W{chs}'.format(chs=chargeSign))),
                                  genBins.Neta, array('d',genBins.etaBins), genBins.Npt, array('d',genBins.ptBins))
        hDiffXsecNormErr= ROOT.TH2F("hDiffXsecNormErr_{lep}_{ch}".format(lep=lepton,ch=charge),
                                     "normalized cross section uncertainty: {Wch}".format(Wch=Wchannel.replace('W','W{chs}'.format(chs=chargeSign))),
                                     genBins.Neta, array('d',genBins.etaBins), genBins.Npt, array('d',genBins.ptBins))

        hDiffXsec_1Deta= ROOT.TH1F("hDiffXsec_1Deta_{lep}_{ch}".format(lep=lepton,ch=charge),
                                            "cross section: {Wch}".format(Wch=Wchannel.replace('W','W{chs}'.format(chs=chargeSign))),
                                            genBins.Neta, array('d',genBins.etaBins))

        hDiffXsecNorm_1Deta= ROOT.TH1F("hDiffXsecNorm_1Deta_{lep}_{ch}".format(lep=lepton,ch=charge),
                                                "cross section: {Wch}".format(Wch=Wchannel.replace('W','W{chs}'.format(chs=chargeSign))),
                                                genBins.Neta, array('d',genBins.etaBins))

        hDiffXsec_1Dpt= ROOT.TH1F("hDiffXsec_1Dpt_{lep}_{ch}".format(lep=lepton,ch=charge),
                                            "cross section: {Wch}".format(Wch=Wchannel.replace('W','W{chs}'.format(chs=chargeSign))),
                                            genBins.Npt, array('d',genBins.ptBins))

        hDiffXsecNorm_1Dpt= ROOT.TH1F("hDiffXsecNorm_1Dpt_{lep}_{ch}".format(lep=lepton,ch=charge),
                                                "cross section: {Wch}".format(Wch=Wchannel.replace('W','W{chs}'.format(chs=chargeSign))),
                                                genBins.Npt, array('d',genBins.ptBins))


        # histogram to keep track of PDF/alphaS variations from the Asimov. 
        # It is like the nominal, but the uncertainty only includes the prefit pdf/alphaS variations
        #hDiffXsecPDF[charge] = ROOT.TH2F("hDiffXsec_{lep}_{ch}_PDF".format(lep=lepton,ch=charge),
        #                                 "cross section PDF: {Wch}".format(Wch=Wchannel.replace('W','W{chs}'.format(chs=chargeSign))),
        #                                 genBins.Neta, array('d',genBins.etaBins), genBins.Npt, array('d',genBins.ptBins))

        # hDiffXsecPDF_1Deta[charge] = ROOT.TH1F("hDiffXsec_{lep}_{ch}_PDF_1Deta".format(lep=lepton,ch=charge),
        #                                     "cross section PDF: {Wch}".format(Wch=Wchannel.replace('W','W{chs}'.format(chs=chargeSign))),
        #                                     genBins.Neta, array('d',genBins.etaBins))
        # hDiffXsecPDF_1Dpt[charge] = ROOT.TH1F("hDiffXsec_{lep}_{ch}_PDF_1Dpt".format(lep=lepton,ch=charge),
        #                                     "cross section PDF: {Wch}".format(Wch=Wchannel.replace('W','W{chs}'.format(chs=chargeSign))),
        #                                     genBins.Npt, array('d',genBins.ptBins))
        hDiffXsecPDF[charge] = ROOT.TGraphAsymmErrors(genBins.Neta * genBins.Npt)
        hDiffXsecPDF[charge].SetName("hDiffXsec_{lep}_{ch}_PDF".format(lep=lepton,ch=charge))
        hDiffXsecPDF[charge].SetTitle("cross section PDF: {Wch}".format(Wch=Wchannel.replace('W','W{chs}'.format(chs=chargeSign))))
        #
        hDiffXsecPDF_1Deta[charge] = ROOT.TGraphAsymmErrors(genBins.Neta)
        hDiffXsecPDF_1Deta[charge].SetName("hDiffXsec_{lep}_{ch}_PDF_1Deta".format(lep=lepton,ch=charge))
        hDiffXsecPDF_1Deta[charge].SetTitle("cross section PDF: {Wch}".format(Wch=Wchannel.replace('W','W{chs}'.format(chs=chargeSign))))
        #
        hDiffXsecPDF_1Dpt[charge] = ROOT.TGraphAsymmErrors(genBins.Npt)
        hDiffXsecPDF_1Dpt[charge].SetName("hDiffXsec_{lep}_{ch}_PDF_1Dpt".format(lep=lepton,ch=charge))
        hDiffXsecPDF_1Dpt[charge].SetTitle("cross section PDF: {Wch}".format(Wch=Wchannel.replace('W','W{chs}'.format(chs=chargeSign))))

        # histogram to keep track of total theory variations from the Asimov
        #hDiffXsecTotTheory[charge] = ROOT.TH2F("hDiffXsec_{lep}_{ch}_TotTheory".format(lep=lepton,ch=charge),
        #                                 "cross section TotTheory: {Wch}".format(Wch=Wchannel.replace('W','W{chs}'.format(chs=chargeSign))),
        #                                 genBins.Neta, array('d',genBins.etaBins), genBins.Npt, array('d',genBins.ptBins))
        # hDiffXsecTotTheory_1Deta[charge] = ROOT.TH1F("hDiffXsec_{lep}_{ch}_TotTheory_1Deta".format(lep=lepton,ch=charge),
        #                                     "cross section TotTheory: {Wch}".format(Wch=Wchannel.replace('W','W{chs}'.format(chs=chargeSign))),
        #                                     genBins.Neta, array('d',genBins.etaBins))
        # hDiffXsecTotTheory_1Dpt[charge] = ROOT.TH1F("hDiffXsec_{lep}_{ch}_TotTheory_1Dpt".format(lep=lepton,ch=charge),
        #                                     "cross section TotTheory: {Wch}".format(Wch=Wchannel.replace('W','W{chs}'.format(chs=chargeSign))),
        #                                     genBins.Npt, array('d',genBins.ptBins))
        hDiffXsecTotTheory[charge] = ROOT.TGraphAsymmErrors(genBins.Neta * genBins.Npt)
        hDiffXsecTotTheory[charge].SetName("hDiffXsec_{lep}_{ch}_TotTheory".format(lep=lepton,ch=charge))
        hDiffXsecTotTheory[charge].SetTitle("cross section TotTheory: {Wch}".format(Wch=Wchannel.replace('W','W{chs}'.format(chs=chargeSign))))
        #
        hDiffXsecTotTheory_1Deta[charge] = ROOT.TGraphAsymmErrors(genBins.Neta)
        hDiffXsecTotTheory_1Deta[charge].SetName("hDiffXsec_{lep}_{ch}_TotTheory_1Deta".format(lep=lepton,ch=charge))
        hDiffXsecTotTheory_1Deta[charge].SetTitle("cross section TotTheory: {Wch}".format(Wch=Wchannel.replace('W','W{chs}'.format(chs=chargeSign))))
        #
        hDiffXsecTotTheory_1Dpt[charge] = ROOT.TGraphAsymmErrors(genBins.Npt)
        hDiffXsecTotTheory_1Dpt[charge].SetName("hDiffXsec_{lep}_{ch}_TotTheory_1Dpt".format(lep=lepton,ch=charge))
        hDiffXsecTotTheory_1Dpt[charge].SetTitle("cross section TotTheory: {Wch}".format(Wch=Wchannel.replace('W','W{chs}'.format(chs=chargeSign))))

        ####
        # XSEC NORM
        #hDiffXsecNormPDF[charge] = ROOT.TH2F("hDiffXsecNorm_{lep}_{ch}_PDF".format(lep=lepton,ch=charge),
        #                                 "cross section PDF: {Wch}".format(Wch=Wchannel.replace('W','W{chs}'.format(chs=chargeSign))),
        #                                 genBins.Neta, array('d',genBins.etaBins), genBins.Npt, array('d',genBins.ptBins))

        # hDiffXsecNormPDF_1Deta[charge] = ROOT.TH1F("hDiffXsecNorm_{lep}_{ch}_PDF_1Deta".format(lep=lepton,ch=charge),
        #                                     "cross section PDF: {Wch}".format(Wch=Wchannel.replace('W','W{chs}'.format(chs=chargeSign))),
        #                                     genBins.Neta, array('d',genBins.etaBins))
        # hDiffXsecNormPDF_1Dpt[charge] = ROOT.TH1F("hDiffXsecNorm_{lep}_{ch}_PDF_1Dpt".format(lep=lepton,ch=charge),
        #                                     "cross section PDF: {Wch}".format(Wch=Wchannel.replace('W','W{chs}'.format(chs=chargeSign))),
        #                                     genBins.Npt, array('d',genBins.ptBins))
        hDiffXsecNormPDF[charge] = ROOT.TGraphAsymmErrors(genBins.Neta*genBins.Npt)
        hDiffXsecNormPDF[charge].SetName("hDiffXsecNorm_{lep}_{ch}_PDF".format(lep=lepton,ch=charge))
        hDiffXsecNormPDF[charge].SetTitle("cross section PDF: {Wch}".format(Wch=Wchannel.replace('W','W{chs}'.format(chs=chargeSign))))
        #
        hDiffXsecNormPDF_1Deta[charge] = ROOT.TGraphAsymmErrors(genBins.Neta)
        hDiffXsecNormPDF_1Deta[charge].SetName("hDiffXsecNorm_{lep}_{ch}_PDF_1Deta".format(lep=lepton,ch=charge))
        hDiffXsecNormPDF_1Deta[charge].SetTitle("cross section PDF: {Wch}".format(Wch=Wchannel.replace('W','W{chs}'.format(chs=chargeSign))))
        #
        hDiffXsecNormPDF_1Dpt[charge] = ROOT.TGraphAsymmErrors(genBins.Npt)
        hDiffXsecNormPDF_1Dpt[charge].SetName("hDiffXsecNorm_{lep}_{ch}_PDF_1Dpt".format(lep=lepton,ch=charge))
        hDiffXsecNormPDF_1Dpt[charge].SetTitle("cross section PDF: {Wch}".format(Wch=Wchannel.replace('W','W{chs}'.format(chs=chargeSign))))

        # histogram to keep track of total theory variations from the Asimov
        #hDiffXsecNormTotTheory[charge] = ROOT.TH2F("hDiffXsecNorm_{lep}_{ch}_TotTheory".format(lep=lepton,ch=charge),
        #                                 "cross section TotTheory: {Wch}".format(Wch=Wchannel.replace('W','W{chs}'.format(chs=chargeSign))),
        #                                 genBins.Neta, array('d',genBins.etaBins), genBins.Npt, array('d',genBins.ptBins))
        # hDiffXsecNormTotTheory_1Deta[charge] = ROOT.TH1F("hDiffXsecNorm_{lep}_{ch}_TotTheory_1Deta".format(lep=lepton,ch=charge),
        #                                     "cross section TotTheory: {Wch}".format(Wch=Wchannel.replace('W','W{chs}'.format(chs=chargeSign))),
        #                                     genBins.Neta, array('d',genBins.etaBins))
        # hDiffXsecNormTotTheory_1Dpt[charge] = ROOT.TH1F("hDiffXsecNorm_{lep}_{ch}_TotTheory_1Dpt".format(lep=lepton,ch=charge),
        #                                     "cross section TotTheory: {Wch}".format(Wch=Wchannel.replace('W','W{chs}'.format(chs=chargeSign))),
        #                                     genBins.Npt, array('d',genBins.ptBins))
        hDiffXsecNormTotTheory[charge] = ROOT.TGraphAsymmErrors(genBins.Neta*genBins.Npt)
        hDiffXsecNormTotTheory[charge].SetName("hDiffXsecNorm_{lep}_{ch}_TotTheory".format(lep=lepton,ch=charge))
        hDiffXsecNormTotTheory[charge].SetTitle("cross section TotTheory: {Wch}".format(Wch=Wchannel.replace('W','W{chs}'.format(chs=chargeSign))))
        #
        hDiffXsecNormTotTheory_1Deta[charge] = ROOT.TGraphAsymmErrors(genBins.Neta)
        hDiffXsecNormTotTheory_1Deta[charge].SetName("hDiffXsecNorm_{lep}_{ch}_TotTheory_1Deta".format(lep=lepton,ch=charge))
        hDiffXsecNormTotTheory_1Deta[charge].SetTitle("cross section TotTheory: {Wch}".format(Wch=Wchannel.replace('W','W{chs}'.format(chs=chargeSign))))
        #
        hDiffXsecNormTotTheory_1Dpt[charge] = ROOT.TGraphAsymmErrors(genBins.Npt)
        hDiffXsecNormTotTheory_1Dpt[charge].SetName("hDiffXsecNorm_{lep}_{ch}_TotTheory_1Dpt".format(lep=lepton,ch=charge))
        hDiffXsecNormTotTheory_1Dpt[charge].SetTitle("cross section TotTheory: {Wch}".format(Wch=Wchannel.replace('W','W{chs}'.format(chs=chargeSign))))



        # only for tests
        hDiffXsecPDFonly_1Deta[charge] = ROOT.TH1F("hDiffXsec_{lep}_{ch}_PDFonly_1Deta".format(lep=lepton,ch=charge),
                                                   "cross section PDFonly: {Wch}".format(Wch=Wchannel.replace('W','W{chs}'.format(chs=chargeSign))),
                                                   genBins.Neta, array('d',genBins.etaBins))
        hDiffXsecPDFonly_1Dpt[charge] = ROOT.TH1F("hDiffXsec_{lep}_{ch}_PDFonly_1Dpt".format(lep=lepton,ch=charge),
                                                  "cross section PDFonly: {Wch}".format(Wch=Wchannel.replace('W','W{chs}'.format(chs=chargeSign))),
                                                  genBins.Npt, array('d',genBins.ptBins))

        hDiffXsecNormPDFonly_1Deta[charge] = ROOT.TH1F("hDiffXsecNorm_{lep}_{ch}_PDFonly_1Deta".format(lep=lepton,ch=charge),
                                                   "cross section PDFonly: {Wch}".format(Wch=Wchannel.replace('W','W{chs}'.format(chs=chargeSign))),
                                                   genBins.Neta, array('d',genBins.etaBins))
        hDiffXsecNormPDFonly_1Dpt[charge] = ROOT.TH1F("hDiffXsecNorm_{lep}_{ch}_PDFonly_1Dpt".format(lep=lepton,ch=charge),
                                                  "cross section PDFonly: {Wch}".format(Wch=Wchannel.replace('W','W{chs}'.format(chs=chargeSign))),
                                                  genBins.Npt, array('d',genBins.ptBins))
        # hDiffXsecAlphaonly_1Deta[charge] = ROOT.TH1F("hDiffXsec_{lep}_{ch}_Alphaonly_1Deta".format(lep=lepton,ch=charge),
        #                                              "cross section Alphaonly: {Wch}".format(Wch=Wchannel.replace('W','W{chs}'.format(chs=chargeSign))),
        #                                              genBins.Neta, array('d',genBins.etaBins))
        # hDiffXsecAlphaonly_1Dpt[charge] = ROOT.TH1F("hDiffXsec_{lep}_{ch}_Alphaonly_1Dpt".format(lep=lepton,ch=charge),
        #                                             "cross section Alphaonly: {Wch}".format(Wch=Wchannel.replace('W','W{chs}'.format(chs=chargeSign))),
        #                                             genBins.Npt, array('d',genBins.ptBins))
        # hDiffXsecQCDonly_1Deta[charge] = ROOT.TH1F("hDiffXsec_{lep}_{ch}_QCDonly_1Deta".format(lep=lepton,ch=charge),
        #                                            "cross section QCDonly: {Wch}".format(Wch=Wchannel.replace('W','W{chs}'.format(chs=chargeSign))),
        #                                            genBins.Neta, array('d',genBins.etaBins))
        # hDiffXsecQCDonly_1Dpt[charge] = ROOT.TH1F("hDiffXsec_{lep}_{ch}_QCDonly_1Dpt".format(lep=lepton,ch=charge),
        #                                           "cross section QCDonly: {Wch}".format(Wch=Wchannel.replace('W','W{chs}'.format(chs=chargeSign))),
        #                                           genBins.Npt, array('d',genBins.ptBins))

        hDiffXsecAlphaonly_1Deta[charge] = ROOT.TGraphAsymmErrors(genBins.Neta)
        hDiffXsecAlphaonly_1Deta[charge].SetName("hDiffXsec_{lep}_{ch}_Alphaonly_1Deta".format(lep=lepton,ch=charge))
        hDiffXsecAlphaonly_1Deta[charge].SetTitle("cross section Alphaonly: {Wch}".format(Wch=Wchannel.replace('W','W{chs}'.format(chs=chargeSign))))

        hDiffXsecQCDonly_1Deta[charge] = ROOT.TGraphAsymmErrors(genBins.Neta)
        hDiffXsecQCDonly_1Deta[charge].SetName("hDiffXsec_{lep}_{ch}_QCDonly_1Deta".format(lep=lepton,ch=charge))
        hDiffXsecQCDonly_1Deta[charge].SetTitle("cross section QCDonly: {Wch}".format(Wch=Wchannel.replace('W','W{chs}'.format(chs=chargeSign))))

        hDiffXsecAlphaonly_1Dpt[charge] = ROOT.TGraphAsymmErrors(genBins.Npt)
        hDiffXsecAlphaonly_1Dpt[charge].SetName("hDiffXsec_{lep}_{ch}_Alphaonly_1Dpt".format(lep=lepton,ch=charge))
        hDiffXsecAlphaonly_1Dpt[charge].SetTitle("cross section Alphaonly: {Wch}".format(Wch=Wchannel.replace('W','W{chs}'.format(chs=chargeSign))))

        hDiffXsecQCDonly_1Dpt[charge] = ROOT.TGraphAsymmErrors(genBins.Npt)
        hDiffXsecQCDonly_1Dpt[charge].SetName("hDiffXsec_{lep}_{ch}_QCDonly_1Dpt".format(lep=lepton,ch=charge))
        hDiffXsecQCDonly_1Dpt[charge].SetTitle("cross section QCDonly: {Wch}".format(Wch=Wchannel.replace('W','W{chs}'.format(chs=chargeSign))))

        hDiffXsecNormAlphaonly_1Deta[charge] = ROOT.TGraphAsymmErrors(genBins.Neta)
        hDiffXsecNormAlphaonly_1Deta[charge].SetName("hDiffXsecNorm_{lep}_{ch}_Alphaonly_1Deta".format(lep=lepton,ch=charge))
        hDiffXsecNormAlphaonly_1Deta[charge].SetTitle("cross section Alphaonly: {Wch}".format(Wch=Wchannel.replace('W','W{chs}'.format(chs=chargeSign))))

        hDiffXsecNormQCDonly_1Deta[charge] = ROOT.TGraphAsymmErrors(genBins.Neta)
        hDiffXsecNormQCDonly_1Deta[charge].SetName("hDiffXsecNorm_{lep}_{ch}_QCDonly_1Deta".format(lep=lepton,ch=charge))
        hDiffXsecNormQCDonly_1Deta[charge].SetTitle("cross section QCDonly: {Wch}".format(Wch=Wchannel.replace('W','W{chs}'.format(chs=chargeSign))))

        hDiffXsecNormAlphaonly_1Dpt[charge] = ROOT.TGraphAsymmErrors(genBins.Npt)
        hDiffXsecNormAlphaonly_1Dpt[charge].SetName("hDiffXsecNorm_{lep}_{ch}_Alphaonly_1Dpt".format(lep=lepton,ch=charge))
        hDiffXsecNormAlphaonly_1Dpt[charge].SetTitle("cross section Alphaonly: {Wch}".format(Wch=Wchannel.replace('W','W{chs}'.format(chs=chargeSign))))

        hDiffXsecNormQCDonly_1Dpt[charge] = ROOT.TGraphAsymmErrors(genBins.Npt)
        hDiffXsecNormQCDonly_1Dpt[charge].SetName("hDiffXsecNorm_{lep}_{ch}_QCDonly_1Dpt".format(lep=lepton,ch=charge))
        hDiffXsecNormQCDonly_1Dpt[charge].SetTitle("cross section QCDonly: {Wch}".format(Wch=Wchannel.replace('W','W{chs}'.format(chs=chargeSign))))


        hMu.SetDirectory(tfOut)
        hMuErr.SetDirectory(tfOut)
        hDiffXsec.SetDirectory(tfOut)
        hDiffXsecErr.SetDirectory(tfOut)
        hDiffXsecNorm.SetDirectory(tfOut)
        hDiffXsecNormErr.SetDirectory(tfOut)
        hDiffXsec_1Deta.SetDirectory(tfOut)
        hDiffXsecNorm_1Deta.SetDirectory(tfOut)
        hDiffXsec_1Dpt.SetDirectory(tfOut)
        hDiffXsecNorm_1Dpt.SetDirectory(tfOut)
        #hDiffXsecPDF[charge].SetDirectory(tfOut)
        #hDiffXsecPDF_1Deta[charge].SetDirectory(tfOut)
        #hDiffXsecPDF_1Dpt[charge].SetDirectory(tfOut)
        #hDiffXsecTotTheory[charge].SetDirectory(tfOut)
        #hDiffXsecTotTheory_1Deta[charge].SetDirectory(tfOut)
        #hDiffXsecTotTheory_1Dpt[charge].SetDirectory(tfOut)
        #hDiffXsecNormPDF[charge].SetDirectory(tfOut)
        #hDiffXsecNormPDF_1Deta[charge].SetDirectory(tfOut)
        #hDiffXsecNormPDF_1Dpt[charge].SetDirectory(tfOut)
        #hDiffXsecNormTotTheory[charge].SetDirectory(tfOut)
        #hDiffXsecNormTotTheory_1Deta[charge].SetDirectory(tfOut)
        #hDiffXsecNormTotTheory_1Dpt[charge].SetDirectory(tfOut)


        # returned object is
        # ret = { "xsec"     : [htheory, htheory_1Deta, htheory_1Dpt],
        #         "xsecnorm" : [htheory_xsecnorm, htheory_1Deta_xsecnorm, htheory_1Dpt_xsecnorm],
        #         "asym"     : [htheory_asym, htheory_1Deta_asym, htheory_1Dpt_asym],
        #         "listkeys" : [key[1] for key in htheory]}
        #
        # where listkeys is pdf1, pdf2, alphaSUp, muRUp, muFDown, ...
        # the rest are the lists of all theory histograms (and nominal) for 2D, and 1D projections, for xsec, xsecnorm and asym
        # they are dictionary with key (charge,theory_nuis), where theory_nuis is returned by listkeys, and charge is all for asym and plus/minus for the xsecs
        
        utilities.checkTheoryBandDiffXsec1Dproj(hDiffXsecPDFonly_1Deta[charge], hDiffXsecAlphaonly_1Deta[charge], hDiffXsecQCDonly_1Deta[charge],
                                                histDict["listkeys"], histDict["xsec"][1], genBins.Neta, charge=charge)
        utilities.checkTheoryBandDiffXsec1Dproj(hDiffXsecPDFonly_1Dpt[charge], hDiffXsecAlphaonly_1Dpt[charge], hDiffXsecQCDonly_1Dpt[charge],
                                                histDict["listkeys"], histDict["xsec"][2], genBins.Npt, charge=charge)

        utilities.checkTheoryBandDiffXsec1Dproj(hDiffXsecNormPDFonly_1Deta[charge], hDiffXsecNormAlphaonly_1Deta[charge], hDiffXsecNormQCDonly_1Deta[charge],
                                                histDict["listkeys"], histDict["xsecnorm"][1], genBins.Neta, charge=charge)
        utilities.checkTheoryBandDiffXsec1Dproj(hDiffXsecNormPDFonly_1Dpt[charge], hDiffXsecNormAlphaonly_1Dpt[charge], hDiffXsecNormQCDonly_1Dpt[charge],
                                                histDict["listkeys"], histDict["xsecnorm"][2], genBins.Npt, charge=charge)
        
        print "Make theory bands for absolute cross section and charge %s" % charge
        utilities.getTheoryBandDiffXsec(hDiffXsecTotTheory[charge], hDiffXsecPDF[charge], histDict["listkeys"], histDict["xsec"][0], 
                                        genBins.Neta, genBins.Npt, charge=charge)
        utilities.getTheoryBandDiffXsec1Dproj(hDiffXsecTotTheory_1Deta[charge], hDiffXsecPDF_1Deta[charge], histDict["listkeys"], histDict["xsec"][1], 
                                              genBins.Neta, charge=charge)
        utilities.getTheoryBandDiffXsec1Dproj(hDiffXsecTotTheory_1Dpt[charge], hDiffXsecPDF_1Dpt[charge], histDict["listkeys"], histDict["xsec"][2], 
                                              genBins.Npt, charge=charge)
        print "DONE"

        print "Make theory bands for normalized cross section and charge %s" % charge
        utilities.getTheoryBandDiffXsec(hDiffXsecNormTotTheory[charge], hDiffXsecNormPDF[charge], histDict["listkeys"], 
                                        histDict["xsecnorm"][0],genBins.Neta, genBins.Npt, charge=charge)
        utilities.getTheoryBandDiffXsec1Dproj(hDiffXsecNormTotTheory_1Deta[charge], hDiffXsecNormPDF_1Deta[charge], histDict["listkeys"], 
                                              histDict["xsecnorm"][1], genBins.Neta, charge=charge)
        utilities.getTheoryBandDiffXsec1Dproj(hDiffXsecNormTotTheory_1Dpt[charge], hDiffXsecNormPDF_1Dpt[charge], histDict["listkeys"], 
                                              histDict["xsecnorm"][2],genBins.Npt, charge=charge)
        print "DONE"

        # utilities.getPDFbandFromXsec(hDiffXsecPDF[charge], charge, xsecChargefile[charge], 
        #                              genBins.Neta, genBins.Npt, firstPtBin=firstPtSignalBin, histoTotTheory=hDiffXsecTotTheory[charge])
        # if options.fitData or options.blindData: 
        #     utilities.getPDFbandFromXsec1D(hDiffXsecPDF_1Deta[charge], charge, xsecChargefile[charge], 
        #                                    genBins.Neta, genBins.Npt, firstOtherVarBin=firstPtSignalBin,isEta=True, 
        #                                    histoTotTheory=hDiffXsecTotTheory_1Deta[charge])
        #     utilities.getPDFbandFromXsec1D(hDiffXsecPDF_1Dpt[charge], charge, xsecChargefile[charge], 
        #                                    genBins.Neta, genBins.Npt, firstVarBin=firstPtSignalBin, isEta=False,
        #                                    histoTotTheory=hDiffXsecTotTheory_1Dpt[charge])


        #nbins = genBins.Neta * (genBins.Npt - nPtBinsBkg)
        binCount = 1        
        print "Now filling histograms with differential cross section (charge = %s)" % charge


        denExpression = ""
        
        if not options.hessian:
            if options.friend != "":
                denExpression = "totxsec_" + charge
            else:
                denExpression = utilities.getDenExpressionForNormDiffXsec(charge, genBins.Neta,genBins.Npt)

        central = 0
        error   = 0

        for ieta in range(1,genBins.Neta+1):
            for ipt in range(1,genBins.Npt+1):
                if ptBinIsBackground[ipt-1]: continue
                if etaBinIsBackground[ieta-1]: continue
                #print "\rBin {num}/{tot}".format(num=binCount,tot=nbins),
                sys.stdout.write('Bin {num}/{tot}   \r'.format(num=binCount,tot=nbins))
                sys.stdout.flush()
                binCount += 1

                # signal strength
                if options.hessian:                    
                    central = utilities.getSignalStrengthFromHessianFast(charge,ieta-1,ipt-1,
                                                                         nHistBins=100, minHist=0.5, maxHist=1.5, tree=tree)                    
                    error = utilities.getSignalStrengthFromHessianFast(charge,ieta-1,ipt-1,
                                                                       nHistBins=200, minHist=0.0, maxHist=1., tree=tree, getErr=True)                    
                else:
                    central,up,down = utilities.getSignalStrengthFromToysFast(charge,ieta-1,ipt-1, 
                                                                              nHistBins=200, minHist=0.5, maxHist=1.5, tree=tree)
                    error = up - central
                hMu.SetBinContent(ieta,ipt,central)        
                hMu.SetBinError(ieta,ipt,error)         
                hMuErr.SetBinContent(ieta,ipt,error)
                
                # normalized cross section
                if options.hessian:                    
                    central = utilities.getNormalizedDiffXsecFromHessianFast(charge,ieta-1,ipt-1,
                                                                             nHistBins=2000, minHist=0., maxHist=200., tree=tree)                    
                    error = utilities.getNormalizedDiffXsecFromHessianFast(charge,ieta-1,ipt-1,
                                                                           nHistBins=2000, minHist=0., maxHist=10., tree=tree, getErr=True)                    
                else:
                    central,up,down = utilities.getNormalizedDiffXsecFromToysFast(charge,ieta-1,ipt-1,
                                                                                  denExpression, nHistBins=1000, minHist=0., maxHist=0.1, tree=tree)
                    error = up - central
                hDiffXsecNorm.SetBinContent(ieta,ipt,central)        
                hDiffXsecNorm.SetBinError(ieta,ipt,error)         
                hDiffXsecNormErr.SetBinContent(ieta,ipt,error)
                
                # unnormalized cross section
                if options.hessian:
                    central = utilities.getDiffXsecFromHessianFast(charge,ieta-1,ipt-1,
                                                                   nHistBins=2000, minHist=0., maxHist=200., tree=tree)                    
                    error = utilities.getDiffXsecFromHessianFast(charge,ieta-1,ipt-1,
                                                                 nHistBins=2000, minHist=0., maxHist=200., tree=tree, getErr=True)                    
                else:                
                    central,up,down = utilities.getDiffXsecFromToysFast(charge,ieta-1,ipt-1,
                                                                        nHistBins=2000, minHist=0., maxHist=200., tree=tree)
                    error = up - central
                hDiffXsec.SetBinContent(ieta,ipt,central)        
                hDiffXsec.SetBinError(ieta,ipt,error)         
                hDiffXsecErr.SetBinContent(ieta,ipt,error)
                
        for ieta in range(1,genBins.Neta+1):
            if etaBinIsBackground[ieta-1]: continue
            central = utilities.getDiffXsec1DFromHessianFast(charge,ieta-1,isIeta=True,nHistBins=5000, minHist=0., maxHist=5000., tree=tree)
            error   = utilities.getDiffXsec1DFromHessianFast(charge,ieta-1,isIeta=True,nHistBins=5000, minHist=0., maxHist=5000., tree=tree, getErr=True)
            hDiffXsec_1Deta.SetBinContent(ieta, central)
            hDiffXsec_1Deta.SetBinError(ieta, error)

            central = utilities.getNormalizedDiffXsec1DFromHessianFast(charge,ieta-1,isIeta=True,nHistBins=1000, minHist=0., maxHist=1., tree=tree)
            error   = utilities.getNormalizedDiffXsec1DFromHessianFast(charge,ieta-1,isIeta=True,nHistBins=1000, minHist=0., maxHist=1., tree=tree, getErr=True)
            hDiffXsecNorm_1Deta.SetBinContent(ieta, central)
            hDiffXsecNorm_1Deta.SetBinError(ieta, error)

        for ipt in range(1,genBins.Npt+1):
            if ptBinIsBackground[ipt-1]: continue
            central = utilities.getDiffXsec1DFromHessianFast(charge,ipt-1,isIeta=False,nHistBins=5000, minHist=0., maxHist=5000., tree=tree)
            error   = utilities.getDiffXsec1DFromHessianFast(charge,ipt-1,isIeta=False,nHistBins=5000, minHist=0., maxHist=5000., tree=tree, getErr=True)
            hDiffXsec_1Dpt.SetBinContent(ipt, central)
            hDiffXsec_1Dpt.SetBinError(ipt, error)

            central = utilities.getNormalizedDiffXsec1DFromHessianFast(charge,ipt-1,isIeta=False,nHistBins=1000, minHist=0., maxHist=1., tree=tree)
            error   = utilities.getNormalizedDiffXsec1DFromHessianFast(charge,ipt-1,isIeta=False,nHistBins=1000, minHist=0., maxHist=1., tree=tree, getErr=True)
            hDiffXsecNorm_1Dpt.SetBinContent(ipt, central)
            hDiffXsecNorm_1Dpt.SetBinError(ipt, error)


        if options.normWidth:
            hDiffXsec.Scale(1.,"width")
            hDiffXsecErr.Scale(1.,"width")
            hDiffXsecNorm.Scale(1.,"width")
            hDiffXsecNormErr.Scale(1.,"width")
            hDiffXsec_1Deta.Scale(1.,"width")
            hDiffXsecNorm_1Deta.Scale(1.,"width")            
            hDiffXsec_1Dpt.Scale(1.,"width")
            hDiffXsecNorm_1Dpt.Scale(1.,"width")            
            normGraphByBinArea2D(hDiffXsecPDF[charge],hDiffXsec)
            normGraphByBinWidth1D(hDiffXsecPDF_1Deta[charge], hDiffXsec_1Deta)
            normGraphByBinWidth1D(hDiffXsecPDF_1Dpt[charge], hDiffXsec_1Dpt)
            #hDiffXsecPDF[charge].Scale(1.,"width")
            #hDiffXsecPDF_1Deta[charge].Scale(1.,"width")
            #hDiffXsecPDF_1Dpt[charge].Scale(1.,"width")
            normGraphByBinArea2D(hDiffXsecTotTheory[charge],hDiffXsec)
            normGraphByBinWidth1D(hDiffXsecTotTheory_1Deta[charge], hDiffXsec_1Deta)
            normGraphByBinWidth1D(hDiffXsecTotTheory_1Dpt[charge], hDiffXsec_1Dpt)            
            #hDiffXsecTotTheory[charge].Scale(1.,"width")
            #hDiffXsecTotTheory_1Deta[charge].Scale(1.,"width")
            #hDiffXsecTotTheory_1Dpt[charge].Scale(1.,"width")
            normGraphByBinArea2D(hDiffXsecNormPDF[charge],hDiffXsecNorm)
            normGraphByBinWidth1D(hDiffXsecNormPDF_1Deta[charge], hDiffXsecNorm_1Deta)
            normGraphByBinWidth1D(hDiffXsecNormPDF_1Dpt[charge], hDiffXsecNorm_1Dpt)
            #hDiffXsecNormPDF[charge].Scale(1.,"width")
            #hDiffXsecNormPDF_1Deta[charge].Scale(1.,"width")
            #hDiffXsecNormPDF_1Dpt[charge].Scale(1.,"width")
            normGraphByBinArea2D(hDiffXsecNormTotTheory[charge],hDiffXsecNorm)
            normGraphByBinWidth1D(hDiffXsecNormTotTheory_1Deta[charge], hDiffXsecNorm_1Deta)
            normGraphByBinWidth1D(hDiffXsecNormTotTheory_1Dpt[charge], hDiffXsecNorm_1Dpt)
            #hDiffXsecNormTotTheory[charge].Scale(1.,"width")
            #hDiffXsecNormTotTheory_1Deta[charge].Scale(1.,"width")
            #hDiffXsecNormTotTheory_1Dpt[charge].Scale(1.,"width")

        if options.lumiNorm > 0:
            scaleFactor = 1./options.lumiNorm
            hDiffXsec.Scale(scaleFactor)
            hDiffXsecErr.Scale(scaleFactor)
            # do not divide the normalized cross section, the scaling factor is already removed
            #hDiffXsecNorm.Scale(scaleFactor)
            #hDiffXsecNormErr.Scale(scaleFactor)
            hDiffXsec_1Deta.Scale(scaleFactor)        
            hDiffXsec_1Dpt.Scale(scaleFactor)        
            #hDiffXsecPDF[charge].Scale(scaleFactor)
            scaleGraphByConstant1D(hDiffXsecPDF[charge], scaleFactor)
            scaleGraphByConstant1D(hDiffXsecPDF_1Deta[charge], scaleFactor)
            scaleGraphByConstant1D(hDiffXsecPDF_1Dpt[charge], scaleFactor)
            #hDiffXsecPDF_1Deta[charge].Scale(scaleFactor)
            #hDiffXsecPDF_1Dpt[charge].Scale(scaleFactor)
            #hDiffXsecTotTheory[charge].Scale(scaleFactor)
            scaleGraphByConstant1D(hDiffXsecTotTheory[charge], scaleFactor)
            scaleGraphByConstant1D(hDiffXsecTotTheory_1Deta[charge], scaleFactor)
            scaleGraphByConstant1D(hDiffXsecTotTheory_1Dpt[charge], scaleFactor)
            #hDiffXsecTotTheory_1Deta[charge].Scale(scaleFactor)
            #hDiffXsecTotTheory_1Dpt[charge].Scale(scaleFactor)

        if isMuElComb:
            if options.combineElePt01asBkg:
                # 3 is the first bin to scale (skip first 2)
                scaleTH2inRange(hDiffXsec,nPtMinToScale=3,scale=0.5) 
                scaleTH2inRange(hDiffXsecErr,nPtMinToScale=3,scale=0.5) 
                scaleTH2inRange(hDiffXsecNorm,nPtMinToScale=3,scale=0.5)
                scaleTH2inRange(hDiffXsecNormErr,nPtMinToScale=3,scale=0.5)
                hDiffXsecNorm_1Deta.Scale(0.5) 
                hDiffXsecNorm_1Dpt.Scale(0.5) 
            else:
                hDiffXsec.Scale(0.5)
                hDiffXsecErr.Scale(0.5)                
            #hDiffXsecNorm.Scale(0.5)
            #hDiffXsecNormErr.Scale(0.5)
            hDiffXsec_1Deta.Scale(0.5)
            #hDiffXsecNorm_1Deta.Scale(0.5)
            hDiffXsec_1Dpt.Scale(0.5)
            #hDiffXsecNorm_1Dpt.Scale(0.5)
            ## no longer divide PFDs and theory, they are taken from a separate file now (which was used to make the shape_xsec.root ones)
            # hDiffXsecPDF[charge].Scale(0.5)
            # hDiffXsecPDF_1Deta[charge].Scale(0.5)
            # hDiffXsecPDF_1Dpt[charge].Scale(0.5)
            # hDiffXsecTotTheory[charge].Scale(0.5)
            # hDiffXsecTotTheory_1Deta[charge].Scale(0.5)
            # hDiffXsecTotTheory_1Dpt[charge].Scale(0.5)


        hDiffXsecRelErr = hDiffXsecErr.Clone(hDiffXsecErr.GetName().replace('XsecErr','XsecRelErr'))
        hDiffXsecRelErr.Divide(hDiffXsec)
        hDiffXsecRelErr.SetTitle(hDiffXsecErr.GetTitle().replace('uncertainty','rel.unc.'))
        hDiffXsecRelErr.SetDirectory(tfOut)

        hDiffXsecNormRelErr = hDiffXsecNormErr.Clone(hDiffXsecNormErr.GetName().replace('XsecNormErr','XsecNormRelErr'))
        hDiffXsecNormRelErr.Divide(hDiffXsecNorm)
        hDiffXsecNormRelErr.SetTitle(hDiffXsecNormErr.GetTitle().replace('uncertainty','rel.unc.'))
        hDiffXsecNormRelErr.SetDirectory(tfOut)

        # now starting to draw the cross sections

        zaxisTitle = "Signal strength #mu::0.98,1.02"
        if options.hessian: zaxisTitle = "Signal strength #mu::0.995,1.005"
        if options.fitData: zaxisTitle = "Signal strength #mu::0.75,1.25"
        drawCorrelationPlot(hMu,
                            xaxisTitle, yaxisTitle, zaxisTitle,
                            hMu.GetName(),
                            "ForceTitle",outname,1,1,False,False,False,1, canvasSize="700,625",leftMargin=0.14,rightMargin=0.22,passCanvas=canvas,palette=options.palette)
        writeHistoIntoFile(hMu,tfOut)

        zaxisTitle = "uncertainty on signal strength #mu::0.0,0.5"
        if options.hessian: zaxisTitle = "uncertainty on signal strength #mu::0.0,0.25"
        drawCorrelationPlot(hMuErr,
                            xaxisTitle, yaxisTitle, zaxisTitle,
                            hMuErr.GetName(),
                            "ForceTitle",outname,1,1,False,False,False,1, canvasSize="700,625",leftMargin=0.14,rightMargin=0.22,passCanvas=canvas,palette=options.palette)
        writeHistoIntoFile(hMuErr,tfOut)

        if options.fitData:
            if charge == "plus": zmin,zmax = 0,130
            else:                zmin,zmax = 0,130
        else:
            if charge == "plus": zmin,zmax = 10,120
            else:                zmin,zmax = 10,100

        zminHist = hDiffXsec.GetBinContent(hDiffXsec.GetMinimumBin())                                                   
        zaxisTitle = "d^{2}#sigma / d|#eta|dp_{T} (pb/GeV)::%.3f,%.3f" % (0.99*(zminHist if zminHist > 1. else zmin),
                                                                    min(zmax,hDiffXsec.GetBinContent(hDiffXsec.GetMaximumBin())))
        #zaxisTitle = "d^{2}#sigma / d|#eta|dp_{T} (pb/GeV)::%.3f,%.3f" % (zmin,zmax)
        drawCorrelationPlot(hDiffXsec,
                            xaxisTitle, yaxisTitle, zaxisTitle,
                            hDiffXsec.GetName(),
                            "ForceTitle",outname,1,1,False,False,False,1, canvasSize="700,625",leftMargin=0.14,rightMargin=0.22,passCanvas=canvas,palette=options.palette)
        writeHistoIntoFile(hDiffXsec,tfOut)

        if options.fitData:
            if charge == "plus": zmin,zmax = 0.5,4.5
            else:                zmin,zmax = 0.5,3.5
        else:
            if charge == "plus": zmin,zmax = 0.5,4.5
            else:                zmin,zmax = 0.5,3.5
        zaxisTitle = "uncertainty on d^{2}#sigma / d|#eta|dp_{T} (pb/GeV)::%.3f,%.3f" % (0.9*hDiffXsecErr.GetMinimum(),min(25,hDiffXsecErr.GetMaximum()))
        #zaxisTitle = "uncertainty on d^{2}#sigma / d|#eta|dp_{T} (pb/GeV)::%.3f,%.3f" % (zmin,zmax)
        drawCorrelationPlot(hDiffXsecErr,
                            xaxisTitle, yaxisTitle, zaxisTitle,
                            hDiffXsecErr.GetName(),
                            "ForceTitle",outname,1,1,False,False,False,1, canvasSize="700,625",leftMargin=0.14,rightMargin=0.22,passCanvas=canvas,palette=options.palette)
        writeHistoIntoFile(hDiffXsecErr,tfOut)

        #zaxisTitle = "rel. uncertainty on d^{2}#sigma / d|#eta|dp_{T}::%.3f,%.3f" % (0.9*hDiffXsecRelErr.GetMinimum(),hDiffXsecRelErr.GetMaximum())
        if options.fitData:
            zaxisTitle = "rel. uncertainty on d^{2}#sigma / d|#eta|dp_{T}::0.010,0.2"
        else:
            zaxisTitle = "rel. uncertainty on d^{2}#sigma / d|#eta|dp_{T}::0.025,0.2"
        drawCorrelationPlot(hDiffXsecRelErr,
                            xaxisTitle, yaxisTitle, zaxisTitle,
                            hDiffXsecRelErr.GetName(),
                            "ForceTitle",outname,1,1,False,False,False,1, canvasSize="700,625",leftMargin=0.14,rightMargin=0.22,passCanvas=canvas,palette=options.palette)
        writeHistoIntoFile(hDiffXsecRelErr,tfOut)

        if charge == "plus": zmin,zmax = 0.002,0.032
        else:                zmin,zmax = 0.002,0.032
        #zaxisTitle = "d^{2}#sigma / d|#eta|dp_{T} / #sigma_{tot}::%.3f,%.3f" % (0.9*hDiffXsecNorm.GetMinimum(),hDiffXsecNorm.GetMaximum())
        zaxisTitle = "d^{2}#sigma / d|#eta|dp_{T} / #sigma_{tot}::%.3f,%.3f" % (zmin,zmax)
        drawCorrelationPlot(hDiffXsecNorm,
                            xaxisTitle, yaxisTitle, zaxisTitle,
                            hDiffXsecNorm.GetName(),
                            "ForceTitle",outname,1,1,False,False,False,1, canvasSize="700,625",leftMargin=0.14,rightMargin=0.22,passCanvas=canvas,palette=options.palette)
        writeHistoIntoFile(hDiffXsecNorm,tfOut)

        if charge == "plus": zmin,zmax = 0.0,0.0015
        else:                zmin,zmax = 0.0,0.0015
        #zaxisTitle = "uncertainty on d^{2}#sigma / d|#eta|dp_{T} / #sigma_{tot}::%.3f,%.3f" % (0.9*hDiffXsecNormErr.GetMinimum(),hDiffXsecNormErr.GetMaximum())
        zaxisTitle = "uncertainty on d^{2}#sigma / d|#eta|dp_{T} / #sigma_{tot}::%.3f,%.3f" % (zmin,zmax)
        drawCorrelationPlot(hDiffXsecNormErr,
                            xaxisTitle, yaxisTitle, zaxisTitle,
                            hDiffXsecNormErr.GetName(),
                            "ForceTitle",outname,1,1,False,False,False,1, canvasSize="700,625",leftMargin=0.14,rightMargin=0.22,passCanvas=canvas,palette=options.palette)
        writeHistoIntoFile(hDiffXsecNormErr,tfOut)

        zaxisTitle = "rel. uncertainty on d^{2}#sigma / d|#eta|dp_{T} / #sigma_{tot}::%.4f,%.4f" % (min(0.01,0.9*hDiffXsecNormRelErr.GetMinimum()),min(0.2,hDiffXsecNormRelErr.GetBinContent(hDiffXsecNormRelErr.GetMaximumBin())))
        #zaxisTitle = "rel. uncertainty on d^{2}#sigma / d|#eta|dp_{T} / #sigma_{tot}::0,0.1"
        drawCorrelationPlot(hDiffXsecNormRelErr,
                            xaxisTitle, yaxisTitle, zaxisTitle,
                            hDiffXsecNormRelErr.GetName(),
                            "ForceTitle",outname,1,1,False,False,False,1, canvasSize="700,625",leftMargin=0.14,rightMargin=0.22,passCanvas=canvas,palette=options.palette)
        writeHistoIntoFile(hDiffXsecNormRelErr,tfOut)

        # with only W -> munu, coordinates can be 0.45,0.8,0.65,0.9 with TPaveText
        texCoord = "0.45,0.85" if charge == "plus" else "0.45,0.5"
        additionalText = "W^{{{chs}}} #rightarrow {lep}#nu;{pttext}::{txc},0.08,0.04".format(chs=" "+chargeSign,
                                                                                             lep="e" if channel == "el" else "#mu" if channel == "mu" else "l",
                                                                                             pttext=ptRangeText,
                                                                                             txc=texCoord) 
        legendCoords = "0.2,0.4,0.75,0.85" if charge == "plus" else "0.2,0.4,0.4,0.5"
        drawSingleTH1(hDiffXsec_1Deta,xaxisTitle,"d#sigma/d|#eta| [pb]",
                      "xsec_eta_abs_{ch}_{fl}".format(ch=charge,fl=channel),
                      outname,labelRatioTmp="Rel.Unc.::0.9,1.1",legendCoords=legendCoords, draw_both0_noLog1_onlyLog2=1,
                      passCanvas=canvas1D, lumi=options.lumiInt,
                      lowerPanelHeight=0.35, moreTextLatex=additionalText
                      )
        writeHistoIntoFile(hDiffXsec_1Deta,tfOut)

        drawSingleTH1(hDiffXsecNorm_1Deta,xaxisTitle,"d#sigma/d|#eta| / #sigma_{tot}",
                      "xsec_eta_norm_{ch}_{fl}".format(ch=charge,fl=channel),
                      outname,labelRatioTmp="Rel.Unc.::0.9,1.1",legendCoords=legendCoords, draw_both0_noLog1_onlyLog2=1,drawLineLowerPanel="",
                      passCanvas=canvas1D, lumi=options.lumiInt,
                      lowerPanelHeight=0.35, moreTextLatex=additionalText
                      )
        writeHistoIntoFile(hDiffXsecNorm_1Deta,tfOut)

        legendCoords = "0.65,0.85,0.75,0.85"
        texCoord = "0.2,0.5"
        additionalText = "W^{{{chs}}} #rightarrow {lep}#nu;{etatext}::{txc},0.08,0.04".format(chs=" "+chargeSign,
                                                                                             lep="e" if channel == "el" else "#mu" if channel == "mu" else "l",
                                                                                             etatext=etaRangeText,
                                                                                             txc=texCoord) 
        drawSingleTH1(hDiffXsec_1Dpt,yaxisTitle,"d#sigma/dp_{T} (pb/GeV)",
                      "xsec_pt_abs_{ch}_{fl}".format(ch=charge,fl=channel),
                      outname,labelRatioTmp="Rel.Unc.::0.9,1.1",legendCoords=legendCoords, draw_both0_noLog1_onlyLog2=1,
                      passCanvas=canvas1D, lumi=options.lumiInt,
                      lowerPanelHeight=0.35, moreTextLatex=additionalText
                      )
        writeHistoIntoFile(hDiffXsec_1Dpt,tfOut)

        drawSingleTH1(hDiffXsecNorm_1Dpt,yaxisTitle,"d#sigma/dp_{T} / #sigma_{tot} (1/GeV)",
                      "xsec_pt_norm_{ch}_{fl}".format(ch=charge,fl=channel),
                      outname,labelRatioTmp="Rel.Unc.::0.9,1.1",legendCoords=legendCoords, draw_both0_noLog1_onlyLog2=1,drawLineLowerPanel="",
                      passCanvas=canvas1D, lumi=options.lumiInt,
                      lowerPanelHeight=0.35, moreTextLatex=additionalText
                      )
        writeHistoIntoFile(hDiffXsecNorm_1Dpt,tfOut)

        writeHistoIntoFile(hDiffXsecPDF[charge],tfOut)
        writeHistoIntoFile(hDiffXsecPDF_1Deta[charge],tfOut)
        writeHistoIntoFile(hDiffXsecPDF_1Dpt[charge],tfOut)
        writeHistoIntoFile(hDiffXsecTotTheory[charge],tfOut)
        writeHistoIntoFile(hDiffXsecTotTheory_1Deta[charge],tfOut)
        writeHistoIntoFile(hDiffXsecTotTheory_1Dpt[charge],tfOut)
        writeHistoIntoFile(hDiffXsecNormPDF[charge],tfOut)
        writeHistoIntoFile(hDiffXsecNormPDF_1Deta[charge],tfOut)
        writeHistoIntoFile(hDiffXsecNormPDF_1Dpt[charge],tfOut)
        writeHistoIntoFile(hDiffXsecNormTotTheory[charge],tfOut)
        writeHistoIntoFile(hDiffXsecNormTotTheory_1Deta[charge],tfOut)
        writeHistoIntoFile(hDiffXsecNormTotTheory_1Dpt[charge],tfOut)
        if charge == "plus" and hasBothCharges:
            writeHistoIntoFile(hChargeAsymPDF,tfOut)
            writeHistoIntoFile(hChargeAsymPDF_1Deta,tfOut)
            writeHistoIntoFile(hChargeAsymPDF_1Dpt,tfOut)
            writeHistoIntoFile(hChargeAsymTotTheory,tfOut)
            writeHistoIntoFile(hChargeAsymTotTheory_1Deta,tfOut)
            writeHistoIntoFile(hChargeAsymTotTheory_1Dpt,tfOut)
            
        

######
        # now drawing a TH1 unrolling TH2
        chsize = 4500
        canvUnroll = ROOT.TCanvas("canvUnroll","",chsize,1500) # 3000,1500
        leftMargin = 0.12 if chsize < 3500 else 0.07
        rightMargin = 0.04 if chsize < 3500 else 0.01

        ratioYaxis = "Rel.Unc.::0.8,1.2"
        if channel == "el": ratioYaxis = "Rel.Unc.::0.8,1.2"

        xaxisTitle = "template global bin"
        vertLinesArg = ""
        # pass x1,y1,x2,y2
        additionalText = "W^{{{chs}}} #rightarrow {lep}#nu::0.78,0.84,0.88,0.9".format(chs=" "+chargeSign, 
                                                                                     lep="e" if channel == "el" else "#mu" if channel == "mu" else "l") 

        h1D_chargeAsym = {}
        h1D_pmaskedexp = {}
        h1D_pmaskedexp_norm = {}

        for unrollAlongEta in [True, False]:

            unrollVar = "eta" if unrollAlongEta else "pt"

            varBinRanges = []

            if unrollAlongEta:
                #xaxisTitle = xaxisTitle + " = 1 + ipt + ieta * %d; ipt in [%d,%d], ieta in [%d,%d]" % (nptbins-1,0,nptbins-1,0,netabins-1)
                xaxisTitle = "unrolled |#eta| bin: |#eta| #in [%.1f, %.1f]" % (genBins.etaBins[0], genBins.etaBins[-1])
                vertLinesArg = "{a},{b}".format(a=genBins.Npt,b=genBins.Neta)
                for ipt in range(0,genBins.Npt):
                    #varBinRanges.append("p_{{T}} #in [{ptmin:.2g}, {ptmax:.3g}]".format(ptmin=genBins.ptBins[ipt], ptmax=genBins.ptBins[ipt+1]))
                    varBinRanges.append("#splitline{{[{ptmin:.2g}, {ptmax:.3g}]}}{{GeV}}".format(ptmin=genBins.ptBins[ipt], ptmax=genBins.ptBins[ipt+1]))
            else:
                #xaxisTitle = xaxisTitle + " = 1 + ieta + ipt * %d; ipt in [%d,%d], ieta in [%d,%d]" % (netabins-1,0,nptbins-1,0,netabins-1)
                xaxisTitle = "unrolled p_{T} bin: p_{T} #in [%.3g, %.3g] GeV" % (genBins.ptBins[0], genBins.ptBins[-1])
                vertLinesArg = "{a},{b}".format(a=genBins.Neta,b=genBins.Npt)
                for ieta in range(0,genBins.Neta):
                    #varBinRanges.append("|#eta| #in [{etamin:.1f}, {etamax:.1f}]".format(etamin=genBins.etaBins[ieta], etamax=genBins.etaBins[ieta+1]))
                    varBinRanges.append("[{etamin:.1f}, {etamax:.1f}]".format(etamin=genBins.etaBins[ieta], etamax=genBins.etaBins[ieta+1]))

            tfOut.cd()
            h1D_pmaskedexp[(unrollAlongEta,charge)] = getTH1fromTH2(hDiffXsec, hDiffXsecErr, unrollAlongX=unrollAlongEta)        
            h1D_pmaskedexp[(unrollAlongEta,charge)].SetDirectory(tfOut)
            drawSingleTH1(h1D_pmaskedexp[(unrollAlongEta,charge)],xaxisTitle,"d^{2}#sigma/d|#eta|dp_{T} (pb/GeV)",
                          "unrolledXsec_{var}_abs_{ch}_{fl}".format(var=unrollVar,ch=charge,fl=channel),
                          outname,labelRatioTmp=ratioYaxis,draw_both0_noLog1_onlyLog2=1,passCanvas=canvUnroll, lumi=options.lumiInt,
                          drawVertLines=vertLinesArg, textForLines=varBinRanges, lowerPanelHeight=0.35, moreText=additionalText,
                          leftMargin=leftMargin, rightMargin=rightMargin)
            writeHistoIntoFile(h1D_pmaskedexp[(unrollAlongEta,charge)],tfOut)
            #h1D_pmaskedexp[unrollAlongEta].Write()

            h1D_pmaskedexp_norm[(unrollAlongEta,charge)] = getTH1fromTH2(hDiffXsecNorm, hDiffXsecNormErr, unrollAlongX=unrollAlongEta)        
            h1D_pmaskedexp_norm[(unrollAlongEta,charge)].SetDirectory(tfOut)
            drawSingleTH1(h1D_pmaskedexp_norm[(unrollAlongEta,charge)],xaxisTitle,"d^{2}#sigma/d|#eta|dp_{T} / #sigma_{tot} (1/GeV)",
                          "unrolledXsec_{var}_norm_{ch}_{fl}".format(var=unrollVar,ch=charge,fl=channel),
                          outname,labelRatioTmp=ratioYaxis,draw_both0_noLog1_onlyLog2=1,passCanvas=canvUnroll, lumi=options.lumiInt, drawLineLowerPanel="",
                          drawVertLines=vertLinesArg, textForLines=varBinRanges, lowerPanelHeight=0.35,moreText=additionalText,
                          leftMargin=leftMargin, rightMargin=rightMargin)
            writeHistoIntoFile(h1D_pmaskedexp_norm[(unrollAlongEta,charge)],tfOut)

            # do also charge asymmetry (only once if running on both charges)
            if icharge == 2:

                #xaxisTitle = xaxisTitle.replace("cross section", "charge asymmetry")
                additionalTextBackup = additionalText
                additionalText = "W #rightarrow {lep}#nu::0.78,0.84,0.88,0.9".format(lep="e" if channel == "el" else "#mu" if channel == "mu" else "l") # pass x1,y1,x2,y2

                h1D_chargeAsym[unrollAlongEta] = getTH1fromTH2(hChAsymm, h2Derr=None, unrollAlongX=unrollAlongEta)        
                h1D_chargeAsym[unrollAlongEta].SetDirectory(tfOut)
                drawSingleTH1(h1D_chargeAsym[unrollAlongEta], xaxisTitle,"charge asymmetry::0.0,1.0",
                              "unrolledChargeAsym_{var}_{fl}".format(var=unrollVar,fl=channel),
                              outname,labelRatioTmp=ratioYaxis,draw_both0_noLog1_onlyLog2=1,drawLineLowerPanel="", 
                              passCanvas=canvUnroll,lumi=options.lumiInt,
                              drawVertLines=vertLinesArg, textForLines=varBinRanges, lowerPanelHeight=0.35, moreText=additionalText,
                              leftMargin=leftMargin, rightMargin=rightMargin)
                writeHistoIntoFile(h1D_chargeAsym[unrollAlongEta],tfOut)

                additionalText = additionalTextBackup


# for data, add plot with ratio with expected
# might be used for toys to compare with asimov
        #if options.fitData:
        if True:
            if options.exptoyfile:
                
                getExpFromGen = False
                if options.exptoyfile == "SAME":
                    treeexp = tree
                    getExpFromGen = True
                else:
                    fexp = ROOT.TFile(options.exptoyfile, 'read')
                    treeexp = fexp.Get('fitresults')

                tfOut.cd() # to make sure the following histograms belong to it, let's see

                #hChargeAsym_exp = None
                # if icharge == 2:
                #     hChargeAsym_exp = ROOT.TH2F("hChargeAsym_{lep}_exp".format(lep=lepton),
                #                                 "charge asymmetry: {Wch}".format(Wch=Wchannel),
                #                                 genBins.Neta, array('d',genBins.etaBins), genBins.Npt, array('d',genBins.ptBins))

                #     hChAsymm1Deta_exp = ROOT.TH1F("hChAsymm1Deta_{lep}_exp".format(lep=lepton),"Charge asymmetry: {Wch}".format(Wch=Wchannel),
                #                                   genBins.Neta, array('d',genBins.etaBins))
                #     hChAsymm1Dpt_exp = ROOT.TH1F("hChAsymm1Dpt_{lep}_exp".format(lep=lepton),"Charge asymmetry: {Wch}".format(Wch=Wchannel),
                #                                   genBins.Npt, array('d',genBins.ptBins))

                #     hChargeAsym_exp.SetDirectory(tfOut)
                #     hChAsymm1Deta_exp.SetDirectory(tfOut)
                #     hChAsymm1Dpt_exp.SetDirectory(tfOut)


                # hDiffXsec_exp = ROOT.TH2F("hDiffXsec_{lep}_{ch}_exp".format(lep=lepton,ch=charge),
                #                           "cross section: {Wch}".format(Wch=Wchannel.replace('W','W{chs}'.format(chs=chargeSign))),
                #                           genBins.Neta, array('d',genBins.etaBins), genBins.Npt, array('d',genBins.ptBins))
                # hDiffXsecNorm_exp = ROOT.TH2F("hDiffXsecNorm_{lep}_{ch}_exp".format(lep=lepton,ch=charge),
                #                               "normalized cross section: {Wch}".format(Wch=Wchannel.replace('W','W{chs}'.format(chs=chargeSign))),
                #                               genBins.Neta, array('d',genBins.etaBins), genBins.Npt, array('d',genBins.ptBins))


                # hDiffXsec_1Deta_exp = ROOT.TH1F("hDiffXsec_1Deta_{lep}_{ch}_exp".format(lep=lepton,ch=charge),
                #                                 "cross section: {Wch}".format(Wch=Wchannel.replace('W','W{chs}'.format(chs=chargeSign))),
                #                                 genBins.Neta, array('d',genBins.etaBins))

                # hDiffXsecNorm_1Deta_exp = ROOT.TH1F("hDiffXsecNorm_1Deta_{lep}_{ch}_exp".format(lep=lepton,ch=charge),
                #                                     "cross section: {Wch}".format(Wch=Wchannel.replace('W','W{chs}'.format(chs=chargeSign))),
                #                                     genBins.Neta, array('d',genBins.etaBins))

                # hDiffXsec_1Dpt_exp = ROOT.TH1F("hDiffXsec_1Dpt_{lep}_{ch}_exp".format(lep=lepton,ch=charge),
                #                                 "cross section: {Wch}".format(Wch=Wchannel.replace('W','W{chs}'.format(chs=chargeSign))),
                #                                 genBins.Npt, array('d',genBins.ptBins))

                # hDiffXsecNorm_1Dpt_exp = ROOT.TH1F("hDiffXsecNorm_1Dpt_{lep}_{ch}_exp".format(lep=lepton,ch=charge),
                #                                     "cross section: {Wch}".format(Wch=Wchannel.replace('W','W{chs}'.format(chs=chargeSign))),
                #                                     genBins.Npt, array('d',genBins.ptBins))

                # hDiffXsec_exp.SetDirectory(tfOut)
                # hDiffXsecNorm_exp.SetDirectory(tfOut)
                # hDiffXsec_1Deta_exp.SetDirectory(tfOut)
                # hDiffXsecNorm_1Deta_exp.SetDirectory(tfOut)
                # hDiffXsec_1Dpt_exp.SetDirectory(tfOut)
                # hDiffXsecNorm_1Dpt_exp.SetDirectory(tfOut)

                # binCount = 0
                # print ""
                # print "Now reading Hessian to make ratio of data with expected"
                print "Now making ratio of data with expected"
                # for ieta in range(1,genBins.Neta+1):
                #     for ipt in range(1,genBins.Npt+1):
                #         if ptBinIsBackground[ipt-1]: continue
                #         if etaBinIsBackground[ieta-1]: continue
                #         #print "\rBin {num}/{tot}".format(num=binCount,tot=nbins),
                #         sys.stdout.write('Bin {num}/{tot}   \r'.format(num=binCount,tot=nbins))
                #         sys.stdout.flush()
                #         binCount += 1
                        
                #         if icharge == 2:
                #             # charge asymmetry
                #             central = utilities.getDiffXsecAsymmetryFromHessianFast(ieta-1,ipt-1,
                #                                                                    nHistBins=2000, minHist=0., maxHist=1., 
                #                                                                    tree=treeexp, getGen=getExpFromGen)         
                #             error = utilities.getDiffXsecAsymmetryFromHessianFast(ieta-1,ipt-1,
                #                                                                  nHistBins=2000, minHist=0., maxHist=1., 
                #                                                                  tree=treeexp, getErr=True)                    
                            
                #             hChargeAsym_exp.SetBinContent(ieta,ipt,central)        
                #             hChargeAsym_exp.SetBinError(ieta,ipt,error)         


                #         # normalized cross section
                #         central = utilities.getNormalizedDiffXsecFromHessianFast(charge,ieta-1,ipt-1,
                #                                                                  nHistBins=2000, minHist=0., maxHist=200., tree=treeexp, getGen=getExpFromGen)         
                #         error = utilities.getNormalizedDiffXsecFromHessianFast(charge,ieta-1,ipt-1,
                #                                                                nHistBins=2000, minHist=0., maxHist=200., tree=treeexp, getErr=True)                    

                #         hDiffXsecNorm_exp.SetBinContent(ieta,ipt,central)        
                #         hDiffXsecNorm_exp.SetBinError(ieta,ipt,error)         

                #         # unnormalized cross section
                #         central = utilities.getDiffXsecFromHessianFast(charge,ieta-1,ipt-1,
                #                                                        nHistBins=2000, minHist=0., maxHist=200., tree=treeexp, getGen=getExpFromGen)                    
                #         error = utilities.getDiffXsecFromHessianFast(charge,ieta-1,ipt-1,
                #                                                      nHistBins=2000, minHist=0., maxHist=200., tree=treeexp, getErr=True)                    
                #         hDiffXsec_exp.SetBinContent(ieta,ipt,central)        
                #         hDiffXsec_exp.SetBinError(ieta,ipt,error)         

                # # now 1D stuff
                # for ieta in range(1,genBins.Neta+1):

                #     if etaBinIsBackground[ieta-1]: continue

                #     if icharge == 2:
                #         central = utilities.getDiffXsecAsymmetry1DFromHessianFast(ieta-1,isIeta=True,nHistBins=1000, minHist=0., maxHist=1.0, tree=treeexp)
                #         error   = utilities.getDiffXsecAsymmetry1DFromHessianFast(ieta-1,isIeta=True,nHistBins=1000, minHist=0., maxHist=1.0, tree=treeexp, getErr=True)
                #         hChAsymm1Deta_exp.SetBinContent(ieta,central)
                #         hChAsymm1Deta_exp.SetBinError(ieta,error)
        
                #     central = utilities.getDiffXsec1DFromHessianFast(charge,ieta-1,isIeta=True,nHistBins=5000, minHist=0., maxHist=5000., tree=treeexp)
                #     error   = utilities.getDiffXsec1DFromHessianFast(charge,ieta-1,isIeta=True,nHistBins=5000, minHist=0., maxHist=5000., tree=treeexp, getErr=True)
                #     hDiffXsec_1Deta_exp.SetBinContent(ieta, central)
                #     hDiffXsec_1Deta_exp.SetBinError(ieta, error)

                #     central = utilities.getNormalizedDiffXsec1DFromHessianFast(charge,ieta-1,isIeta=True,nHistBins=1000, minHist=0., maxHist=1., tree=treeexp)
                #     error   = utilities.getNormalizedDiffXsec1DFromHessianFast(charge,ieta-1,isIeta=True,nHistBins=1000, minHist=0., maxHist=1., tree=treeexp, getErr=True)
                #     hDiffXsecNorm_1Deta_exp.SetBinContent(ieta, central)
                #     hDiffXsecNorm_1Deta_exp.SetBinError(ieta, error)

                # for ipt in range(1,genBins.Npt+1):

                #     if ptBinIsBackground[ipt-1]: continue

                #     if icharge == 2:
                #         central = utilities.getDiffXsecAsymmetry1DFromHessianFast(ipt-1,isIeta=False,nHistBins=1000, minHist=0., maxHist=1.0, tree=treeexp)
                #         error   = utilities.getDiffXsecAsymmetry1DFromHessianFast(ipt-1,isIeta=False,nHistBins=1000, minHist=0., maxHist=1.0, tree=treeexp, getErr=True)
                #         hChAsymm1Dpt_exp.SetBinContent(ipt,central)
                #         hChAsymm1Dpt_exp.SetBinError(ipt,error)
        
                #     central = utilities.getDiffXsec1DFromHessianFast(charge,ipt-1,isIeta=False,nHistBins=5000, minHist=0., maxHist=5000., tree=treeexp)
                #     error   = utilities.getDiffXsec1DFromHessianFast(charge,ipt-1,isIeta=False,nHistBins=5000, minHist=0., maxHist=5000., tree=treeexp, getErr=True)
                #     hDiffXsec_1Dpt_exp.SetBinContent(ipt, central)
                #     hDiffXsec_1Dpt_exp.SetBinError(ipt, error)

                #     central = utilities.getNormalizedDiffXsec1DFromHessianFast(charge,ipt-1,isIeta=False,nHistBins=1000, minHist=0., maxHist=1., tree=treeexp)
                #     error   = utilities.getNormalizedDiffXsec1DFromHessianFast(charge,ipt-1,isIeta=False,nHistBins=1000, minHist=0., maxHist=1., tree=treeexp, getErr=True)
                #     hDiffXsecNorm_1Dpt_exp.SetBinContent(ipt, central)
                #     hDiffXsecNorm_1Dpt_exp.SetBinError(ipt, error)


                # if options.normWidth:
                #     hDiffXsec_exp.Scale(1.,"width")
                #     hDiffXsecNorm_exp.Scale(1.,"width")
                #     hDiffXsec_1Deta_exp.Scale(1.,"width")
                #     hDiffXsecNorm_1Deta_exp.Scale(1.,"width")
                #     hDiffXsec_1Dpt_exp.Scale(1.,"width")
                #     hDiffXsecNorm_1Dpt_exp.Scale(1.,"width")
                    

                # if options.lumiNorm > 0:
                #     scaleFactor = 1./options.lumiNorm
                #     hDiffXsec_exp.Scale(scaleFactor)
                #     hDiffXsec_1Deta_exp.Scale(scaleFactor)
                #     hDiffXsec_1Dpt_exp.Scale(scaleFactor)

                # if isMuElComb:
                #     if options.combineElePt01asBkg:
                #         scaleTH2inRange(hDiffXsec_exp,nPtMinToScale=3,scale=0.5) 
                #         scaleTH2inRange(hDiffXsecNorm_exp,nPtMinToScale=3,scale=0.5) 
                #         hDiffXsecNorm_1Deta_exp.Scale(0.5)
                #         hDiffXsecNorm_1Dpt_exp.Scale(0.5)
                #     else:
                #         hDiffXsec_exp.Scale(0.5)
                #     hDiffXsec_1Deta_exp.Scale(0.5)
                #     hDiffXsec_1Dpt_exp.Scale(0.5)
                #     #hDiffXsecNorm_exp.Scale(0.5)
                #     #hDiffXsecNorm_1Deta_exp.Scale(0.5)
                #     #hDiffXsecNorm_1Dpt_exp.Scale(0.5)


                # # now drawing a TH1 unrolling TH2
                #h1D_pmaskedexp_exp = {}
                #h1D_pmaskedexp_norm_exp = {}
                #h1D_chargeAsym_exp = {}
                h1D_pmaskedexp_PDF = {}
                h1D_pmaskedexp_TotTheory = {}
                h1D_pmaskedexp_norm_PDF = {}
                h1D_pmaskedexp_norm_TotTheory = {}
                h1D_chargeAsym_PDF = {}
                h1D_chargeAsym_TotTheory = {}

                for unrollAlongEta in [True, False]:

                    unrollVar = "eta" if unrollAlongEta else "pt"

                    varBinRanges = []

                    if unrollAlongEta:
                        #xaxisTitle = xaxisTitle + " = 1 + ipt + ieta * %d; ipt in [%d,%d], ieta in [%d,%d]" % (nptbins-1,0,nptbins-1,0,netabins-1)
                        xaxisTitle = "unrolled |#eta| bin: |#eta| #in [%.1f, %.1f]" % (genBins.etaBins[0], genBins.etaBins[-1])
                        vertLinesArg = "{a},{b}".format(a=genBins.Npt,b=genBins.Neta)
                        for ipt in range(0,genBins.Npt):
                            #varBinRanges.append("p_{{T}} #in [{ptmin:.2g}, {ptmax:.3g}]".format(ptmin=genBins.ptBins[ipt], ptmax=genBins.ptBins[ipt+1]))
                            varBinRanges.append("#splitline{{[{ptmin:.2g}, {ptmax:.3g}]}}{{GeV}}".format(ptmin=genBins.ptBins[ipt], ptmax=genBins.ptBins[ipt+1]))
                    else:
                        #xaxisTitle = xaxisTitle + " = 1 + ieta + ipt * %d; ipt in [%d,%d], ieta in [%d,%d]" % (netabins-1,0,nptbins-1,0,netabins-1)
                        xaxisTitle = "unrolled p_{T} bin: p_{T} #in [%.3g, %.3g] GeV" % (genBins.ptBins[0], genBins.ptBins[-1])
                        vertLinesArg = "{a},{b}".format(a=genBins.Neta,b=genBins.Npt)
                        for ieta in range(0,genBins.Neta):
                            #varBinRanges.append("|#eta| #in [{etamin:.1f}, {etamax:.1f}]".format(etamin=genBins.etaBins[ieta], etamax=genBins.etaBins[ieta+1]))
                            varBinRanges.append("[{etamin:.1f}, {etamax:.1f}]".format(etamin=genBins.etaBins[ieta], etamax=genBins.etaBins[ieta+1]))

                    #h1D_pmaskedexp_PDF[(unrollAlongEta,charge)] = getTH1fromTH2(hDiffXsecPDF[charge], h2Derr=None, unrollAlongX=unrollAlongEta)   
                    #h1D_pmaskedexp_TotTheory[(unrollAlongEta,charge)] = getTH1fromTH2(hDiffXsecTotTheory[charge], h2Derr=None, unrollAlongX=unrollAlongEta) 
                    #h1D_pmaskedexp_PDF[(unrollAlongEta,charge)].SetDirectory(tfOut) 
                    #h1D_pmaskedexp_TotTheory[(unrollAlongEta,charge)].SetDirectory(tfOut) 
                    #h1D_pmaskedexp_exp[(unrollAlongEta,charge)] = getTH1fromTH2(hDiffXsec_exp, h2Derr=None, unrollAlongX=unrollAlongEta)        
                    #h1D_pmaskedexp_exp[(unrollAlongEta,charge)].SetDirectory(tfOut)
                    h1D_pmaskedexp_PDF[(unrollAlongEta,charge)] = getUnrolledGraph(hDiffXsecPDF[charge], genBins.Neta, genBins.Npt,
                                                                                   unrollAlongX=unrollAlongEta)
                    h1D_pmaskedexp_TotTheory[(unrollAlongEta,charge)] = getUnrolledGraph(hDiffXsecTotTheory[charge], genBins.Neta, genBins.Npt,
                                                                                         unrollAlongX=unrollAlongEta)
                    drawXsecAndTheoryband(h1D_pmaskedexp[(unrollAlongEta,charge)], h1D_pmaskedexp_TotTheory[(unrollAlongEta,charge)],
                                          xaxisTitle,"d^{2}#sigma/d|#eta|dp_{T} (pb/GeV)::%.3f,200" % (0 if unrollAlongEta else 0.0),
                                          "unrolledXsec_{var}_abs_{ch}_{fl}_dataAndExp".format(var=unrollVar,ch=charge,fl=channel),
                                          outname,labelRatioTmp=labelRatioDataExp,draw_both0_noLog1_onlyLog2=1,passCanvas=canvUnroll,lumi=options.lumiInt,
                                          drawVertLines=vertLinesArg, textForLines=varBinRanges, lowerPanelHeight=0.4, moreText=additionalText,
                                          leftMargin=leftMargin, rightMargin=rightMargin, invertRatio=options.invertRatio,
                                          histMCpartialUnc=h1D_pmaskedexp_PDF[(unrollAlongEta,charge)],histMCpartialUncLegEntry="PDFs #oplus #alpha_{S}",skipPreliminary=options.skipPreliminary, nSplitUnrolledBins=options.nSplitUnrolledBins)


                    #h1D_pmaskedexp_norm_PDF[(unrollAlongEta,charge)] = getTH1fromTH2(hDiffXsecNormPDF[charge], h2Derr=None, unrollAlongX=unrollAlongEta)  
                    #h1D_pmaskedexp_norm_TotTheory[(unrollAlongEta,charge)] = getTH1fromTH2(hDiffXsecNormTotTheory[charge], h2Derr=None, 
                    #                                                                       unrollAlongX=unrollAlongEta)  
                    #h1D_pmaskedexp_norm_PDF[(unrollAlongEta,charge)].SetDirectory(tfOut)  
                    #h1D_pmaskedexp_norm_TotTheory[(unrollAlongEta,charge)].SetDirectory(tfOut)
                    #h1D_pmaskedexp_norm_exp[(unrollAlongEta,charge)] = getTH1fromTH2(hDiffXsecNorm_exp, h2Derr=None, unrollAlongX=unrollAlongEta)        
                    #h1D_pmaskedexp_norm_exp[(unrollAlongEta,charge)].SetDirectory(tfOut)
                    h1D_pmaskedexp_norm_PDF[(unrollAlongEta,charge)] = getUnrolledGraph(hDiffXsecNormPDF[charge], genBins.Neta, genBins.Npt,
                                                                                        unrollAlongX=unrollAlongEta)
                    h1D_pmaskedexp_norm_TotTheory[(unrollAlongEta,charge)] = getUnrolledGraph(hDiffXsecNormTotTheory[charge], 
                                                                                              genBins.Neta, genBins.Npt, unrollAlongX=unrollAlongEta)
                    drawXsecAndTheoryband(h1D_pmaskedexp_norm[(unrollAlongEta,charge)], h1D_pmaskedexp_norm_TotTheory[(unrollAlongEta,charge)],
                                          xaxisTitle,"d^{2}#sigma/d|#eta|dp_{T} / #sigma_{tot} (1/GeV)::%s" % ("0.0,0.045" if unrollAlongEta else "0.0,0.045"),
                                          "unrolledXsec_{var}_norm_{ch}_{fl}_dataAndExp".format(var=unrollVar,ch=charge,fl=channel),
                                          outname,labelRatioTmp=labelRatioDataExp,draw_both0_noLog1_onlyLog2=1,passCanvas=canvUnroll, lumi=options.lumiInt,
                                          drawVertLines=vertLinesArg, textForLines=varBinRanges, lowerPanelHeight=0.4, moreText=additionalText,
                                          leftMargin=leftMargin, rightMargin=rightMargin, invertRatio=options.invertRatio,
                                          histMCpartialUnc=h1D_pmaskedexp_norm_PDF[(unrollAlongEta,charge)],histMCpartialUncLegEntry="PDFs #oplus #alpha_{S}",skipPreliminary=options.skipPreliminary,nSplitUnrolledBins=options.nSplitUnrolledBins)


                    # do also charge asymmetry (only once if running on both charges)
                    if icharge == 2:
                        
                        additionalTextBackup = additionalText
                        # pass x1,y1,x2,y2
                        additionalText = "W #rightarrow {lep}#nu::0.78,0.84,0.88,0.9".format(lep="e" if channel == "el" else "#mu" if channel == "mu" else "l") 
                        
                        labelRatioDataExp_asym = labelRatioDataExp
                        #labelRatioDataExp_asym = str(labelRatioDataExp.split("::")[0]) + "::0.8,1.5"
                        labelRatioDataExp_asym = str(labelRatioDataExp.split("::")[0]) + "::-0.1,0.1" # use difference
                        labelRatioDataExp_asym = labelRatioDataExp_asym.replace("/","-") # use difference
                        #h1D_chargeAsym_PDF[unrollAlongEta] = getTH1fromTH2(hChargeAsymPDF, h2Derr=None, unrollAlongX=unrollAlongEta)   
                        #h1D_chargeAsym_TotTheory[unrollAlongEta] = getTH1fromTH2(hChargeAsymTotTheory, h2Derr=None, unrollAlongX=unrollAlongEta)  
                        #h1D_chargeAsym_PDF[unrollAlongEta].SetDirectory(tfOut) 
                        #h1D_chargeAsym_TotTheory[unrollAlongEta].SetDirectory(tfOut) 
                        #h1D_chargeAsym_exp[unrollAlongEta] = getTH1fromTH2(hChargeAsym_exp, h2Derr=None, unrollAlongX=unrollAlongEta)        
                        #h1D_chargeAsym_exp[unrollAlongEta].SetDirectory(tfOut)
                        h1D_chargeAsym_PDF[unrollAlongEta] = getUnrolledGraph(hChargeAsymPDF,genBins.Neta, genBins.Npt, unrollAlongX=unrollAlongEta)
                        h1D_chargeAsym_TotTheory[unrollAlongEta] = getUnrolledGraph(hChargeAsymTotTheory,genBins.Neta, genBins.Npt,
                                                                                    unrollAlongX=unrollAlongEta) 
                        drawXsecAndTheoryband(h1D_chargeAsym[unrollAlongEta], h1D_chargeAsym_TotTheory[unrollAlongEta],
                                              xaxisTitle.replace("cross section", "charge asymmetry"),
                                              "charge asymmetry::0,0.65", "unrolledChargeAsym_{var}_{fl}_dataAndExp".format(var=unrollVar,fl=channel),
                                              outname,labelRatioTmp=labelRatioDataExp_asym,draw_both0_noLog1_onlyLog2=1,passCanvas=canvUnroll,
                                              lumi=options.lumiInt,drawVertLines=vertLinesArg, textForLines=varBinRanges, 
                                              lowerPanelHeight=0.4, moreText=additionalText,leftMargin=leftMargin, rightMargin=rightMargin, 
                                              invertRatio=options.invertRatio, useDifferenceInLowerPanel=True,
                                              histMCpartialUnc=h1D_chargeAsym_PDF[unrollAlongEta],histMCpartialUncLegEntry="PDFs #oplus #alpha_{S}",skipPreliminary=options.skipPreliminary,nSplitUnrolledBins=options.nSplitUnrolledBins)

                        additionalText = additionalTextBackup

                # now data/prediction for the 1D xsection and charge asymmetry
                ###############
                xaxisTitle = 'dressed %s |#eta|' % lepton
                additionalTextBackup = additionalText
                #additionalText = "W^{{{chs}}} #rightarrow {lep}#nu;{pttext}::0.45,0.8,0.75,0.9".format(chs=" "+chargeSign,lep="e" if channel == "el" else "#mu" if channel == "mu" else "l",
                #                                                                                       pttext=ptRangeText) # pass x1,y1,x2,y2
                texCoord = "0.6,0.85" if charge == "plus" else "0.6,0.85"
                additionalText = "W^{{{chs}}} #rightarrow {lep}#nu;{pttext}::{txc},0.08,0.045".format(chs=" "+chargeSign,
                                                                                                     lep="e" if channel == "el" else "#mu" if channel == "mu" else "l",
                                                                                                     pttext=ptRangeText,
                                                                                                     txc=texCoord) 
                legendCoords = "0.18,0.58,0.75,0.85" if charge == "plus" else "0.18,0.58,0.4,0.5" # does not include space for PDF line (created automatically inside)                
                drawXsecAndTheoryband(hDiffXsec_1Deta,hDiffXsecTotTheory_1Deta[charge],xaxisTitle,"d#sigma/d|#eta| [pb]",
                                      "xsec_eta_abs_{ch}_{fl}_dataAndExp".format(ch=charge,fl=channel),
                                      outname,labelRatioTmp=labelRatioDataExp,legendCoords=legendCoords, draw_both0_noLog1_onlyLog2=1,
                                      passCanvas=canvas1D, lumi=options.lumiInt,
                                      lowerPanelHeight=0.35, moreTextLatex=additionalText, 
                                      invertRatio=options.invertRatio, 
                                      histMCpartialUnc=hDiffXsecPDF_1Deta[charge], histMCpartialUncLegEntry="PDFs #oplus #alpha_{S}",skipPreliminary=options.skipPreliminary)
 
                drawXsecAndTheoryband(hDiffXsecNorm_1Deta,hDiffXsecNormTotTheory_1Deta[charge],xaxisTitle,"d#sigma/d|#eta| / #sigma_{tot}",
                                      "xsec_eta_norm_{ch}_{fl}_dataAndExp".format(ch=charge,fl=channel),
                                      outname,labelRatioTmp=labelRatioDataExp,legendCoords=legendCoords, draw_both0_noLog1_onlyLog2=1,
                                      passCanvas=canvas1D, lumi=options.lumiInt,
                                      lowerPanelHeight=0.35, moreTextLatex=additionalText, invertRatio=options.invertRatio,
                                      histMCpartialUnc=hDiffXsecNormPDF_1Deta[charge], histMCpartialUncLegEntry="PDFs #oplus #alpha_{S}",skipPreliminary=options.skipPreliminary)

                theoryText = "W^{{{chs}}} #rightarrow {lep}#nu;{etatext}::0.2,0.3,0.08,0.06".format(chs=" "+chargeSign,
                                                                                             lep="e" if channel == "el" else "#mu" if channel == "mu" else "l",
                                                                                             etatext=ptRangeText) 
                drawCheckTheoryBand(hDiffXsecPDFonly_1Deta[charge], hDiffXsecAlphaonly_1Deta[charge], hDiffXsecQCDonly_1Deta[charge],
                                    xaxisTitle,"rel. unc.::0.94,1.06", "theoryBands_eta_abs_{ch}_{fl}".format(ch=charge,fl=channel),
                                    outname,draw_both0_noLog1_onlyLog2=1,lumi=options.lumiInt,moreTextLatex=theoryText)
                drawCheckTheoryBand(hDiffXsecNormPDFonly_1Deta[charge], hDiffXsecNormAlphaonly_1Deta[charge], hDiffXsecNormQCDonly_1Deta[charge],
                                    xaxisTitle,"rel. unc.::0.94,1.06", "theoryBands_eta_norm_{ch}_{fl}".format(ch=charge,fl=channel),
                                    outname,draw_both0_noLog1_onlyLog2=1,lumi=options.lumiInt,moreTextLatex=theoryText)

                #legendCoords = "0.65,0.85,0.75,0.85" 
                legendCoords = "0.18,0.58,0.4,0.5" 
                texCoord = "0.65,0.85"
                additionalText = "W^{{{chs}}} #rightarrow {lep}#nu;{etatext}::{txc},0.08,0.045".format(chs=" "+chargeSign,
                                                                                                     lep="e" if channel == "el" else "#mu" if channel == "mu" else "l",
                                                                                                     etatext=etaRangeText,
                                                                                                     txc=texCoord) 
                drawXsecAndTheoryband(hDiffXsec_1Dpt,hDiffXsecTotTheory_1Dpt[charge],yaxisTitle,"d#sigma/dp_{T} (pb/GeV)",
                                      "xsec_pt_abs_{ch}_{fl}_dataAndExp".format(ch=charge,fl=channel),
                                      outname,labelRatioTmp=labelRatioDataExp,legendCoords=legendCoords, draw_both0_noLog1_onlyLog2=1,
                                      passCanvas=canvas1D, lumi=options.lumiInt,
                                      lowerPanelHeight=0.35, moreTextLatex=additionalText, 
                                      invertRatio=options.invertRatio, 
                                      histMCpartialUnc=hDiffXsecPDF_1Dpt[charge], histMCpartialUncLegEntry="PDFs #oplus #alpha_{S}",skipPreliminary=options.skipPreliminary)

                drawXsecAndTheoryband(hDiffXsecNorm_1Dpt,hDiffXsecNormTotTheory_1Dpt[charge],yaxisTitle,"d#sigma/dp_{T} / #sigma_{tot} (1/GeV)",
                                      "xsec_pt_norm_{ch}_{fl}_dataAndExp".format(ch=charge,fl=channel),
                                      outname,labelRatioTmp=labelRatioDataExp,legendCoords=legendCoords, draw_both0_noLog1_onlyLog2=1,
                                      passCanvas=canvas1D, lumi=options.lumiInt,
                                      lowerPanelHeight=0.35, moreTextLatex=additionalText, invertRatio=options.invertRatio,
                                      histMCpartialUnc=hDiffXsecNormPDF_1Dpt[charge], histMCpartialUncLegEntry="PDFs #oplus #alpha_{S}",skipPreliminary=options.skipPreliminary)

                theoryText = "W^{{{chs}}} #rightarrow {lep}#nu;{etatext}::0.2,0.3,0.08,0.06".format(chs=" "+chargeSign,
                                                                                             lep="e" if channel == "el" else "#mu" if channel == "mu" else "l",
                                                                                             etatext=etaRangeText) 

                drawCheckTheoryBand(hDiffXsecPDFonly_1Dpt[charge], hDiffXsecAlphaonly_1Dpt[charge], hDiffXsecQCDonly_1Dpt[charge],
                                    yaxisTitle,"rel. unc.::0.94,1.06", "theoryBands_pt_abs_{ch}_{fl}".format(ch=charge,fl=channel),
                                    outname,draw_both0_noLog1_onlyLog2=1,lumi=options.lumiInt,moreTextLatex=theoryText)

                drawCheckTheoryBand(hDiffXsecNormPDFonly_1Dpt[charge], hDiffXsecNormAlphaonly_1Dpt[charge], hDiffXsecNormQCDonly_1Dpt[charge],
                                    yaxisTitle,"rel. unc.::0.94,1.06", "theoryBands_pt_norm_{ch}_{fl}".format(ch=charge,fl=channel),
                                    outname,draw_both0_noLog1_onlyLog2=1,lumi=options.lumiInt,moreTextLatex=theoryText)

                additionalText = additionalTextBackup

                if icharge == 2:
                                        
                    additionalTextBackup = additionalText
                    #additionalText = "W #rightarrow {lep}#nu;{pttext}::0.2,0.6,0.5,0.7".format(lep="e" if channel == "el" else "#mu",pttext=ptRangeText) # pass x1,y1,x2,y2
                    texCoord = "0.2,0.68"
                    additionalText = "W #rightarrow {lep}#nu;{pttext}::{txc},0.08,0.045".format(lep="e" if channel == "el" else "#mu" if channel == "mu" else "l", 
                                                                                               pttext=ptRangeText,
                                                                                               txc=texCoord)
                    legendCoords = "0.2,0.6,0.75,0.85"
                    # now building PDF band for charge asymmetry (this is wrong, doesn't include correlations)
                    #hChAsymmPDF_1D = utilities.getChargeAsyFromTH1pair(hDiffXsecPDF_1Deta["plus"],hDiffXsecPDF_1Deta["minus"],name="asymPDF")
                    labelRatioDataExp_asym = labelRatioDataExp
                    labelRatioDataExp_asym = str(labelRatioDataExp.split("::")[0]) + "::-0.01,0.01"
                    labelRatioDataExp_asym = labelRatioDataExp_asym.replace("/","-") # use difference   
                    drawXsecAndTheoryband(hChAsymm1Deta, hChargeAsymTotTheory_1Deta, xaxisTitle, "charge asymmetry::0.06,0.25",
                                          "chargeAsym1D_eta_{fl}_dataAndExp".format(fl=channel),
                                          outname,labelRatioTmp=labelRatioDataExp_asym,draw_both0_noLog1_onlyLog2=1,
                                          passCanvas=canvas1D,lumi=options.lumiInt, legendCoords=legendCoords,
                                          lowerPanelHeight=0.35, moreTextLatex=additionalText, 
                                          invertRatio=options.invertRatio, useDifferenceInLowerPanel=True,
                                          histMCpartialUnc=hChargeAsymPDF_1Deta, histMCpartialUncLegEntry="PDFs #oplus #alpha_{S}",skipPreliminary=options.skipPreliminary)
                    
                    legendCoords = "0.48,0.88,0.75,0.85"
                    texCoord = "0.2,0.48"
                    additionalText = "W #rightarrow {lep}#nu;{pttext}::{txc},0.08,0.045".format(lep="e" if channel == "el" else "#mu" if channel == "mu" else "l", 
                                                                                               pttext=etaRangeText,
                                                                                               txc=texCoord)
                    labelRatioDataExp_asym = str(labelRatioDataExp.split("::")[0]) + "::-0.05,0.05"
                    labelRatioDataExp_asym = labelRatioDataExp_asym.replace("/","-") # use difference   
                    drawXsecAndTheoryband(hChAsymm1Dpt, hChargeAsymTotTheory_1Dpt, yaxisTitle, "charge asymmetry::0.0,0.25",
                                          "chargeAsym1D_pt_{fl}_dataAndExp".format(fl=channel),
                                          outname,labelRatioTmp=labelRatioDataExp_asym,draw_both0_noLog1_onlyLog2=1,
                                          passCanvas=canvas1D,lumi=options.lumiInt, legendCoords=legendCoords,
                                          lowerPanelHeight=0.35, moreTextLatex=additionalText, 
                                          invertRatio=options.invertRatio, useDifferenceInLowerPanel=True,
                                          histMCpartialUnc=hChargeAsymPDF_1Dpt, histMCpartialUncLegEntry="PDFs #oplus #alpha_{S}",skipPreliminary=options.skipPreliminary)

                    additionalText = additionalTextBackup
                
                
                if options.exptoyfile != "SAME": fexp.Close()

        print ""


    #tfOut.Write()
    tfOut.Close()
    print ""
    print "Created file %s" % (outname+outfilename)
    print ""
