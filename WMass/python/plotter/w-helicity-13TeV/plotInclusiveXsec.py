#!/bin/env python

import ROOT, os, sys, re, array, math

from make_diff_xsec_cards import getDiffXsecBinning
from make_diff_xsec_cards import templateBinning

from plotDiffXsecChargeAsymmetry import scaleGraphByConstant1D
from plotDiffXsecChargeAsymmetry import writeHistoIntoFile

sys.path.append(os.getcwd() + "/plotUtils/")
from utility import *

ROOT.gROOT.SetBatch(True)

import utilities
utilities = utilities.util()

if __name__ == "__main__":


    from optparse import OptionParser
    parser = OptionParser(usage='%prog [options]')
    parser.add_option('-i','--input', dest='inputdir', default='', type='string', help='input directory with all the cards inside. It is used to get other information')
    parser.add_option('-o','--outdir', dest='outdir', default='', type='string', help='output directory to save things')
    parser.add_option('-t','--toyfile', dest='toyfile', default='.', type='string', help='Root file with toys.')
    parser.add_option('-c','--channel', dest='channel', default='', type='string', help='name of the channel (mu or el). For combination, select lep')
    parser.add_option('-s','--suffix', dest='suffix', default='', type='string', help='Suffix added to output dir (i.e, either to hessian or toys in the path)')
    parser.add_option('-l','--lumi-norm', dest='lumiNorm', default='-1', type='float', help='If > 0, divide cross section by this factor (lumi in 1/Pb)')
    parser.add_option(     '--lumiInt', dest='lumiInt', default='35.9', type='float', help='Integrated luminosity')
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
    parser.add_option(      '--skipPreliminary', dest='skipPreliminary', default=False, action='store_true', help='Do not add "Preliminary" to text on top left of canvas')
    parser.add_option(     '--do-expected', dest='doExpected' , default=False , action='store_true',   help='If True, add "_expected" suffix to output name')
    parser.add_option(     '--no-BBB', dest='noBBB' , default=False , action='store_true',   help='If True, add "_noBBB" suffix to output name')
    (options, args) = parser.parse_args()

    ROOT.TH1.SetDefaultSumw2()
    ROOT.TH1.StatOverflows(True)
    
    channel = options.channel
    if channel not in ["el","mu","lep"]:
        print "Error: unknown channel %s (select 'el' or 'mu' or 'lep')" % channel
        quit()

    isMuElComb = False
    if channel == "lep": isMuElComb = True


    if options.blindData:
        if len(options.exptoyfile) == 0:
            print "WARNING: you used option --blind-data, but you didn't specify a hessian file with option --expected-toyfile. Exit"
            quit()
        if options.fitData:
            print "WARNING: options --blind-data and --fit-data are not compatible. The former will effectively make plots with Asimov. Exit"
            quit()            

    if not options.hessian:
        print "WARNING: options --hessian not set. Code currently supports only that option. Exit"
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


    pname = "inclusiveXsecAndAsymmetry"
    if options.doExpected: pname += "_expected"
    if options.noBBB: pname += "_noBBB"
    # this file will store all histograms for later usage
    outfilename = outname+pname+".root"  # "plotInclusiveXsec.root"
    tfOut = ROOT.TFile.Open(outfilename,'recreate')
    tfOut.cd()

    lepton = "electron" if channel == "el" else "muon" if channel == "mu" else "lepton"
    Wchannel = "W #rightarrow %s#nu" % ("e" if channel == "el" else "#mu" if channel == "mu" else "l")

    ptRangeText = "p_{{T}}^{{{l}}} #in [{ptmin:2g}, {ptmax:3g}] GeV".format(l="e" if channel == "el" else "#mu" if channel == "mu" else "l", 
                                                                            ptmin=genBins.ptBins[firstPtSignalBin], ptmax=genBins.ptBins[-1])
    #etaRangeText = "|#eta|^{{{l}}} #in [{etamin:3g}, {etamax:3g}] GeV".format(l="e" if channel == "el" else "#mu" if channel == "mu" else "l", 
    #                                                                          etamin=genBins.etaBins[0], etamax=genBins.etaBins[-1])
    etaRangeText = "|#eta^{{{l}}}| < {etamax:3g}".format(l="e" if channel == "el" else "#mu" if channel == "mu" else "l", 
                                                                              etamin=genBins.etaBins[0], etamax=genBins.etaBins[-1])

    # just one value
    hChAsymm = ROOT.TH1F("hChargeAsym_{lep}".format(lep=lepton),
                         "Charge asymmetry: {Wch}".format(Wch=Wchannel),
                         1,0,1)
    hChAsymm.SetDirectory(tfOut)
    #
    hChargeRatio = ROOT.TH1F("hChargeRatio_{lep}".format(lep=lepton),
                             "Charge ratio +/-: {Wch}".format(Wch=Wchannel),
                             1,0,1)
    hChargeRatio.SetDirectory(tfOut)
    #
    hXsecPlus = ROOT.TH1F("hXsec_{lep}_plus".format(lep=lepton),
                         "cross section: {Wch}".format(Wch=Wchannel.replace('W','W+')),
                         1,0,1)
    hXsecPlus.SetDirectory(tfOut)
    #
    hXsecMinus = ROOT.TH1F("hXsec_{lep}_minus".format(lep=lepton),
                           "cross section: {Wch}".format(Wch=Wchannel.replace('W','W-')),
                           1,0,1)
    hXsecMinus.SetDirectory(tfOut)
    #
    hXsec = ROOT.TH1F("hXsec_{lep}".format(lep=lepton),
                      "cross section: {Wch}".format(Wch=Wchannel),
                      1,0,1)
    hXsec.SetDirectory(tfOut)
    

    ##########################
    # graphs for theory uncertainty
    hChargeAsymPDF = ROOT.TGraphAsymmErrors(1)
    hChargeAsymPDF.SetName("hChargeAsym_{lep}_PDF".format(lep=lepton))
    hChargeAsymPDF.SetTitle("Charge asymmetry: {Wch}".format(Wch=Wchannel))
    #
    hChargeRatioPDF = ROOT.TGraphAsymmErrors(1)
    hChargeRatioPDF.SetName("hChargeRatio_{lep}_PDF".format(lep=lepton))
    hChargeRatioPDF.SetTitle("Charge ratio: {Wch}".format(Wch=Wchannel))
    #
    hXsecPlusPDF = ROOT.TGraphAsymmErrors(1)
    hXsecPlusPDF.SetName("hXsec_{lep}_plus_PDF".format(lep=lepton))
    hXsecPlusPDF.SetTitle("cross section PDF: {Wch}".format(Wch=Wchannel.replace('W','W+')))
    #
    hXsecMinusPDF = ROOT.TGraphAsymmErrors(1)
    hXsecMinusPDF.SetName("hXsec_{lep}_minus_PDF".format(lep=lepton))
    hXsecMinusPDF.SetTitle("cross section PDF: {Wch}".format(Wch=Wchannel.replace('W','W-')))
    #
    hXsecPDF = ROOT.TGraphAsymmErrors(1)
    hXsecPDF.SetName("hXsec_{lep}_minus_PDF".format(lep=lepton))
    hXsecPDF.SetTitle("cross section PDF: {Wch}".format(Wch=Wchannel))
    #
    hChargeAsymTotTheory = ROOT.TGraphAsymmErrors(1)
    hChargeAsymTotTheory.SetName("hChargeAsym_{lep}_TotTheory".format(lep=lepton))
    hChargeAsymTotTheory.SetTitle("Charge asymmetry: {Wch}".format(Wch=Wchannel))
    #
    hChargeRatioTotTheory = ROOT.TGraphAsymmErrors(1)
    hChargeRatioTotTheory.SetName("hChargeRatio_{lep}_TotTheory".format(lep=lepton))
    hChargeRatioTotTheory.SetTitle("Charge ratio: {Wch}".format(Wch=Wchannel))
    #
    hXsecPlusTotTheory = ROOT.TGraphAsymmErrors(1)
    hXsecPlusTotTheory.SetName("hXsec_{lep}_plus_TotTheory".format(lep=lepton))
    hXsecPlusTotTheory.SetTitle("cross section: {Wch}".format(Wch=Wchannel.replace('W','W+')))
    #
    hXsecMinusTotTheory = ROOT.TGraphAsymmErrors(1)
    hXsecMinusTotTheory.SetName("hXsec_{lep}_minus_TotTheory".format(lep=lepton))
    hXsecMinusTotTheory.SetTitle("cross section: {Wch}".format(Wch=Wchannel.replace('W','W-')))
    #
    hXsecTotTheory = ROOT.TGraphAsymmErrors(1)
    hXsecTotTheory.SetName("hXsec_{lep}_minus_TotTheory".format(lep=lepton))
    hXsecTotTheory.SetTitle("cross section: {Wch}".format(Wch=Wchannel))
    #
    ##########################


    ptminTheoryHist = options.ptRange.split(",")[0] if options.ptRange != "template" else -1.0
    if options.forceAllPtBinsTheoryBand: ptminTheoryHist = -1.0
    print "Now retrieving all the theory variations, to make the theory bands later on"
    histDict = utilities.getTheoryInclusiveXsecFast(options.useXsecWithWptWeights, 
                                                    ptmin=ptminTheoryHist)
    print "DONE, all the histograms for the theory bands were taken"
    # sanity check
    print "="*30
    print histDict["listkeys"]
    print "Found %d keys" % len(histDict["listkeys"])
    print "="*30


    print "Make theory bands for inclusive asymmetry"
    utilities.getTheoryBandInclusiveXsec(hChargeAsymTotTheory, hChargeAsymPDF, 
                                         histDict["listkeys"], histDict["asym"][0], 
                                         charge="all", valOnXaxis = True)
    print "Make theory bands for inclusive ratio"
    utilities.getTheoryBandInclusiveXsec(hChargeRatioTotTheory, hChargeRatioPDF, 
                                         histDict["listkeys"], histDict["ratio"][0], 
                                         charge="all", valOnXaxis = True)
    print "Make theory bands for absolute cross section and charge plus"
    utilities.getTheoryBandInclusiveXsec(hXsecPlusTotTheory, hXsecPlusPDF, 
                                         histDict["listkeys"], histDict["xsec"][0], 
                                         charge="plus", valOnXaxis = True)
    print "Make theory bands for absolute cross section and charge minus"
    utilities.getTheoryBandInclusiveXsec(hXsecMinusTotTheory, hXsecMinusPDF, 
                                         histDict["listkeys"], histDict["xsec"][0], 
                                         charge="minus", valOnXaxis = True)
    print "Make theory bands for absolute cross section summing both charges"
    utilities.getTheoryBandInclusiveXsec(hXsecTotTheory, hXsecPDF, 
                                         histDict["listkeys"], histDict["xsecIncl"][0], 
                                         charge="incl", valOnXaxis = True)

    # now reading file with measured value
    # note that we no longer need to read file with expected, the theory band already has the correct central value
    mainfile = options.toyfile
    if options.blindData:        
        mainfile = options.exptoyfile

    f = ROOT.TFile(mainfile, 'read')
    tree = f.Get('fitresults')
    
    central = utilities.getInclusiveAsymmetryFromHessianFast(nHistBins=2000, minHist=0., maxHist=1.0, tree=tree)
    error = utilities.getInclusiveAsymmetryFromHessianFast(nHistBins=2000, minHist=0., maxHist=1.0, tree=tree, getErr=True)
    hChAsymm.SetBinContent(1,central)    
    hChAsymm.SetBinError(1,error)    

    central = utilities.getInclusiveRatioFromHessianFast(nHistBins=2000, minHist=0., maxHist=1.0, tree=tree)
    error = utilities.getInclusiveRatioFromHessianFast(nHistBins=2000, minHist=0., maxHist=1.0, tree=tree, getErr=True)
    hChargeRatio.SetBinContent(1,central)    
    hChargeRatio.SetBinError(1,error)    


    central = utilities.getInclusiveXsecFromHessianFast("plus",nHistBins=10000, minHist=0., maxHist=50000., tree=tree)
    error = utilities.getInclusiveXsecFromHessianFast("plus",nHistBins=10000, minHist=0., maxHist=50000., tree=tree, getErr=True)                    
    hXsecPlus.SetBinContent(1,central)    
    hXsecPlus.SetBinError(1,error)    

    central = utilities.getInclusiveXsecFromHessianFast("minus",nHistBins=10000, minHist=0., maxHist=50000., tree=tree)
    error = utilities.getInclusiveXsecFromHessianFast("minus",nHistBins=10000, minHist=0., maxHist=50000., tree=tree, getErr=True)                    
    hXsecMinus.SetBinContent(1,central)    
    hXsecMinus.SetBinError(1,error)    

    central = utilities.getInclusiveXsecFromHessianFast("incl",nHistBins=10000, minHist=0., maxHist=50000., tree=tree)
    error = utilities.getInclusiveXsecFromHessianFast("incl",nHistBins=10000, minHist=0., maxHist=50000., tree=tree, getErr=True)                    
    hXsec.SetBinContent(1,central)    
    hXsec.SetBinError(1,error)    

    scaleFactor = 1.0
    if options.lumiNorm > 0:
        scaleFactor = 1./options.lumiNorm
        hXsecPlus.Scale(scaleFactor)
        hXsecMinus.Scale(scaleFactor)
        hXsec.Scale(scaleFactor)
        scaleGraphByConstant1D(hXsecPlusPDF, scaleFactor, scaleX=True)
        scaleGraphByConstant1D(hXsecMinusPDF, scaleFactor, scaleX=True)
        scaleGraphByConstant1D(hXsecPDF, scaleFactor, scaleX=True)
        scaleGraphByConstant1D(hXsecPlusTotTheory, scaleFactor, scaleX=True)
        scaleGraphByConstant1D(hXsecMinusTotTheory, scaleFactor, scaleX=True)
        scaleGraphByConstant1D(hXsecTotTheory, scaleFactor, scaleX=True)

    if isMuElComb:
        hXsecPlus.Scale(0.5)
        hXsecMinus.Scale(0.5)
        hXsec.Scale(0.5)

    writeHistoIntoFile(hChAsymm,tfOut)
    writeHistoIntoFile(hChargeRatio,tfOut)
    writeHistoIntoFile(hXsecPlus,tfOut)
    writeHistoIntoFile(hXsecMinus,tfOut)
    writeHistoIntoFile(hXsec,tfOut)
    writeHistoIntoFile(hChargeAsymPDF,tfOut)
    writeHistoIntoFile(hChargeRatioPDF,tfOut)
    writeHistoIntoFile(hXsecPlusPDF,tfOut)
    writeHistoIntoFile(hXsecMinusPDF,tfOut)
    writeHistoIntoFile(hXsecPDF,tfOut)
    writeHistoIntoFile(hChargeAsymTotTheory,tfOut)
    writeHistoIntoFile(hChargeRatioTotTheory,tfOut)
    writeHistoIntoFile(hXsecPlusTotTheory,tfOut)
    writeHistoIntoFile(hXsecMinusTotTheory,tfOut)
    writeHistoIntoFile(hXsecTotTheory,tfOut)

    ################
    # getting some impacts for inclusive cross section
    #key, target, label #(label is the name of the x axis label for that the TH2 for the proper target)
    impactDict = {"incl" :  ["sumpois",        "W_lep_sumxsec",            None],
                  "plus" :  ["sumpois",        "Wplus_lep_sumxsec",        None],
                  "minus" : ["sumpois",        "Wminus_lep_sumxsec",       None],
                  "ratio" : ["ratiometapois", "Wratio_lep_ratiometaratio", None],
                  "asym" :  ["chargemetapois", "Wasym_lep_chargemetaasym", None]
                }

    # could just have exp+theory (I have the proper impact already, it has everything except stat(s) and lumi)
    #impactToDisplay = ["luminosity", "allTheory", "allExp", "binByBinStat", "stat"]
    #impactTextLabel = ["lumi", "theo", "exp", "stat(MC)", "stat"]
    impactToDisplay = ["luminosity", "allTheoAndExp", "binByBinStat", "stat", "allExp", "allTheory"]
    impactTextLabel = ["lumi", "ExpAndTheoSyst", "stat(MC)", "stat", "expSyst", "theoSyst"]
    for key in impactDict.keys():
        th2name = 'nuisance_group_impact_{sfx}'.format(sfx=impactDict[key][0])
        impMat = f.Get(th2name)
        if impMat==None:
            print "ERROR: Cannot find the impact TH2 named ",th2name," in the input file. Maybe you didn't run --doImpacts?\nSkipping."
            sys.exit()
        for ib in range(1,impMat.GetNbinsX()+1):
            if impMat.GetXaxis().GetBinLabel(ib) == impactDict[key][1]: # this will be true only once
                impactDict[key][2] = ROOT.TH1D("impact_"+key,"",impMat.GetNbinsY(),0,impMat.GetNbinsY())
                for iy in range(1,impMat.GetNbinsY()+1):                    
                    impactDict[key][2].GetXaxis().SetBinLabel(iy,impMat.GetYaxis().GetBinLabel(iy))
                    impactVal = impMat.GetBinContent(ib,iy)
                    if key != "asym" and key != "ratio":
                        impactVal = impactVal * scaleFactor
                        if isMuElComb:
                            # xsec was twice the value for a single lepton
                            impactVal = impactVal * 0.5
                    impactDict[key][2].SetBinContent(iy,impactVal)
            else:
                continue

    ################


    dictHist = {"incl"  : [hXsec,      hXsecTotTheory,       hXsecPDF],
                "plus"  : [hXsecPlus,  hXsecPlusTotTheory,   hXsecPlusPDF],
                "minus" : [hXsecMinus, hXsecMinusTotTheory,  hXsecMinusPDF],
                "asym"  : [hChAsymm,   hChargeAsymTotTheory, hChargeAsymPDF],
                "ratio"  : [hChargeRatio,   hChargeRatioTotTheory, hChargeRatioPDF]
            }

    thislep = "e" if channel == "el" else "#mu" if channel == "mu" else "l"
    texts = {"incl"  : "W #rightarrow %s^{ }#nu" % thislep,
             "plus"  : "W^{ +} #rightarrow %s^{ + }#nu" % thislep,
             "minus" : "W^{ -} #rightarrow %s^{ - }#bar{#nu}" % thislep,
             "ratio"  : "W^{ + }/ W^{ - }",
             "asym"  : "(W^{ + }- W^{ - }) / (W^{ + }+ W^{ - })"
    }

    #sortedKeys = ["asym", "minus", "plus", "incl"]  # last one will appear on top
    sortedKeys = ["ratio", "minus", "plus", "incl"]  # last one will appear on top
    #sortedKeys = ["plus"]

    # at this point the histograms exist, need to plot them with partial uncertainties as well
    # for that, we need to access the impact plots
    # will plot value over theory central value, so that everythin is centered at 1
    # then I will write the actual numbers in the canvas.
    # for that, I need to revert the histograms and graphs
    adjustSettings_CMS_lumi()
    leftMargin = 0.04
    rightMargin = 0.04
    ROOT.gStyle.SetEndErrorSize(5)
    canvas = ROOT.TCanvas("canvas","",2500,1500)   
    canvas.SetTickx(0)
    canvas.SetTicky(1)
    canvas.cd()
    canvas.SetLeftMargin(leftMargin)
    canvas.SetRightMargin(rightMargin)
    canvas.SetTopMargin(0.08)
    canvas.SetBottomMargin(0.12)
    canvas.cd()
    xmin = 0.87
    xmax = 1.20
    ymin = 0
    ymax = 7

    # for theory band
    colorBandPart = {"line" : ROOT.kAzure+7,
                     "fill" : ROOT.kAzure+7} # try +6 and 50% transparent
    colorBandTot = {"line" : ROOT.kPink-6,
                    "fill" : ROOT.kPink-6} # try kRed+1 with 50% transparent #  + 6 # try -3, -9, +1, -4, -2, +6
    fcalphaPart = 0.6 # 0.5
    fcalphaTot = 0.66 # 0.9
    fsPart = 1001
    fsTot = 1001
    ########
    text = ROOT.TLatex()
    text.SetTextColor(ROOT.kBlack)
    text.SetTextSize(0.04)  # 0.03
    text.SetTextFont(42)

    xsectext = ROOT.TLatex()
    xsectext.SetTextColor(ROOT.kBlack)
    xsectext.SetTextSize(0.04)  # 0.03
    xsectext.SetTextFont(42)    
    
    cuttext = ROOT.TLatex()
    cuttext.SetTextColor(ROOT.kBlack)
    cuttext.SetTextSize(0.045)  # 0.03
    cuttext.SetTextFont(42)

    theoryTextColor = ROOT.kPink-6 # or kBlue+1
    
    ystart = 1
    nhists = len(sortedKeys)
    gr = {}
    htotTheory = {}
    hpdf = {}
    hpdfdummyWhite = {}

    xtext = 1.08
    latMeas = ROOT.TLatex()
    latMeas.SetTextFont(42)
    latMeas.SetTextSize(0.035)

    textfile = open(outname+pname+".txt","w+")
    textfile.write("#Measured quantity and uncertainty breakdown (note that we show theory and exp both separately and combined)\n")        
    textfile.write("#Uncertainties do not trivially sum in quadratue, as they are correlated.\n\n")        

    for ik,key in enumerate(sortedKeys):
        thisy = ystart + ik
        h = dictHist[key][0].Clone("dummy")

        # no, do not sum in qudrature the uncertainties, they do not sum in quadrature
        # totUncNoLumi = 0.0
        # uncLumi = 0.0
        iph = impactDict[key][2]
        # for ip in impactToDisplay:
        #     if ip == "luminosity": 
        #         uncLumi = iph.GetBinContent(iph.GetXaxis().FindFixBin(ip))
        #     else:
        #         unc = iph.GetBinContent(iph.GetXaxis().FindFixBin(ip))
        #         totUncNoLumi += unc*unc
        # totUncNoLumi = ROOT.TMath.Sqrt(totUncNoLumi)
        # if key == "asym" or key == "ratio":
        #     xsectextData = "%.4f #pm %.4f" % (h.GetBinContent(1), h.GetBinError(1))
        # else:
        #     xsectextData = "%.0f #pm %.0f (%.0f #oplus %.0f) pb" % (h.GetBinContent(1), h.GetBinError(1), uncLumi, totUncNoLumi)
        if key == "asym" or key == "ratio":
            xsectextData = "%.4f #pm %.4f" % (h.GetBinContent(1), h.GetBinError(1))
        else:
            xsectextData = "%.0f #pm %.0f pb" % (h.GetBinContent(1), h.GetBinError(1))

        fline = key
        if key == "asym" or key == "ratio":
            fline += ": {c:.4f} #pm {e:.4f}     ".format(c=h.GetBinContent(1),e=h.GetBinError(1))
            for i,ip in enumerate(impactToDisplay):
                fline += "   {v:.4f} ({itl})".format(v=iph.GetBinContent(iph.GetXaxis().FindFixBin(ip)),
                                                       itl=impactTextLabel[i])
        else:
            fline += ":    {c:.0f} #pm {e:.0f}       ".format(c=h.GetBinContent(1),e=h.GetBinError(1))
            for i,ip in enumerate(impactToDisplay):
                fline += "   {v:.0f} ({itl})".format(v=iph.GetBinContent(iph.GetXaxis().FindFixBin(ip)),
                                                       itl=impactTextLabel[i])
        fline += "\n\n"            
        textfile.write(fline)

        htotTheory[ik] = dictHist[key][1]
        hpdf[ik] = dictHist[key][2]
        centralTheoryVal = ROOT.Double(0)
        centralTheoryPointVal = ROOT.Double(0) # not needed anyway
        htotTheory[ik].GetPoint(0,centralTheoryVal,centralTheoryPointVal)
        if key == "asym" or key == "ratio":
            xsectextTheo = "%.4f^{ +%.4f}_{ -%.4f}" % (centralTheoryVal, 
                                                                      htotTheory[ik].GetErrorXhigh(0), 
                                                                      htotTheory[ik].GetErrorXlow(0))
        else:
            xsectextTheo = "%.0f^{ +%.0f}_{ -%.0f} pb" % (centralTheoryVal, 
                                                                   htotTheory[ik].GetErrorXhigh(0), 
                                                                   htotTheory[ik].GetErrorXlow(0))
        # if key == "asym":
        #     xsectextTheo = "%.4f + %.4f - %.4f" % (centralTheoryVal, 
        #                                                               htotTheory[ik].GetErrorXhigh(0), 
        #                                                               htotTheory[ik].GetErrorXlow(0))
        # else:
        #     xsectextTheo = "%.0f + %.0f - %.0f pb" % (centralTheoryVal, 
        #                                                            htotTheory[ik].GetErrorXhigh(0), 
        #                                                            htotTheory[ik].GetErrorXlow(0))

        scaleTheoryVal = 1./centralTheoryVal
        #print "x = %f   y = %f" % (centralTheoryVal, centralTheoryPointVal)
        scaleGraphByConstant1D(htotTheory[ik], scaleTheoryVal, scaleX=True)        
        scaleGraphByConstant1D(hpdf[ik], scaleTheoryVal, scaleX=True)        
        gr[ik] = ROOT.TGraphAsymmErrors(1)
        gr[ik].SetName(h.GetName()+ "_%d" % ik)
        gr[ik].SetTitle(h.GetTitle())
        gr[ik].SetPoint(0,h.GetBinContent(1)*scaleTheoryVal, thisy)
        gr[ik].SetPointError(0,h.GetBinError(1)*scaleTheoryVal, h.GetBinError(1)*scaleTheoryVal, 0,0)
        gr[ik].SetLineColor(ROOT.kBlack)
        gr[ik].SetMarkerColor(ROOT.kBlack)
        gr[ik].SetMarkerStyle(20)
        gr[ik].SetLineWidth(3)
        if key == "ratio":
            ROOT.gStyle.SetEndErrorSize(0)
            gr[ik].SetMarkerSize(2)
        else:
            ROOT.gStyle.SetEndErrorSize(5)
            gr[ik].SetMarkerSize(3)
        htotTheory[ik].SetPoint(0,1, thisy)                 
        hpdf[ik].SetPoint(0,1, thisy)                 

        if ik == 0:
            frame = ROOT.TH1D("frame","",100,0.5,1.5)        
            frame.SetLineColor(ROOT.kWhite)
            frame.SetMarkerColor(ROOT.kWhite)
            frame.SetMarkerSize(0)
            frame.GetXaxis().SetLabelSize(0.05)
            frame.GetXaxis().SetTitle("ratio (measured/prediction) of fiducial cross sections")
            frame.GetXaxis().SetTitleOffset(1.2)
            frame.GetXaxis().SetTitleSize(0.05)
            frame.GetXaxis().SetRangeUser(xmin,xmax)
            # remove Y axis title, labels and ticks
            frame.GetYaxis().SetTitle("")
            frame.GetYaxis().SetTitleOffset(999)
            frame.GetYaxis().SetTitleSize(0)
            frame.GetYaxis().SetLabelSize(0)
            frame.GetYaxis().SetRangeUser(ymin, ymax)    
            frame.GetYaxis().SetTickSize(0)
            frame.Draw("HIST")
            canvas.Update()
            ##        
        htotTheory[ik].SetFillStyle(fsTot)  # 1001 for solid, 3001 is better than 3002 for pdf, while 3002 is perfect for png
        htotTheory[ik].SetLineColor(colorBandTot["line"]) 
        htotTheory[ik].SetFillColor(colorBandTot["fill"]) 
        htotTheory[ik].SetFillColorAlpha(colorBandTot["fill"], fcalphaTot);
        htotTheory[ik].Draw("E2 SAME")
        ## draw a solid white as partial band, so that if the partial band has a non solid style, it won't be
        # spoilt by the total one that is drawn below it
        hpdfdummyWhite[ik] = hpdf[ik].Clone("hpdfdummyWhite_%d" % ik)
        hpdfdummyWhite[ik].SetFillColor(ROOT.kWhite)
        hpdfdummyWhite[ik].SetLineColor(ROOT.kWhite)
        hpdfdummyWhite[ik].SetFillStyle(1001)
        hpdfdummyWhite[ik].Draw("E2 SAME")
        hpdf[ik].SetFillStyle(fsPart)  # 1001 for solid, 3001 is better than 3002 for pdf, while 3002 is perfect for png
        hpdf[ik].SetLineColor(colorBandPart["line"]) 
        hpdf[ik].SetFillColor(colorBandPart["fill"]) 
        hpdf[ik].SetFillColorAlpha(colorBandPart["fill"], fcalphaPart);
        hpdf[ik].Draw("E2 SAME")
        gr[ik].Draw("E1P1 SAME") # redraw to have it on top        

        text.DrawLatex(xmin+0.01, thisy-0.1, texts[key])
        xsectext.SetTextColor(ROOT.kBlack)
        xsectext.DrawLatex(xtext, thisy+0.1, xsectextData)
        xsectext.SetTextColor(theoryTextColor)
        #xsectext.DrawLatex(xtext, thisy-0.4, xsectextTheo) # if no sub/super-scripts
        xsectext.DrawLatex(xtext, thisy-0.3, xsectextTheo) # with sub/super-scripts 
        if ik == (len(sortedKeys)-1):
            thisy = ystart + len(sortedKeys)
            text.DrawLatex(xmin+0.01, thisy-0.1, "#bf{Observable}")
            latMeas.SetTextColor(ROOT.kBlack)
            #latMeas.DrawLatex(xtext, thisy+0.1, "#bf{Measured #pm uncertainty (lumi #oplus syst)}")
            latMeas.DrawLatex(xtext, thisy+0.1, "#bf{Measured #pm total uncertainty}")
            latMeas.SetTextColor(theoryTextColor)
            latMeas.DrawLatex(xtext, thisy-0.4, "#bf{Theory prediction}")
    # end loop

    horizline = ROOT.TLine(xmin,0,xmax,0)
    horizline.SetLineColor(ROOT.kBlack)
    horizline.SetLineStyle(3) # 2 for denser hatches
    for iy in range(ystart, ystart + nhists + 1):
        yval = float(iy) - 0.5
        horizline.DrawLine(xmin,yval,xmax,yval)
    horizline.SetLineStyle(0)
    horizline.SetLineWidth(2)
    # horizline.DrawLine(xmin,yval,xmax,yval) # draw after legend

    leg = ROOT.TLegend(0.06,0.76,0.5,0.91)
    leg.SetFillColor(0)
    #leg.SetFillStyle(0)
    leg.SetBorderSize(0)
    leg.SetNColumns(1)
    leg.AddEntry(gr[0],"measured","LP")
    leg.AddEntry(htotTheory[ik],"MadGraph5_aMC@NLO","LF")
    leg.AddEntry(hpdf[ik],"PDFs #oplus #alpha_{S}","LF")
    leg.Draw("same")

    yval = float(ystart + nhists + 1) - 0.5
    vertline = ROOT.TLine(1,ymin, 1, yval)
    vertline.SetLineColor(ROOT.kRed+2)
    vertline.SetLineWidth(2) # 2 for denser hatches
    vertline.DrawLine(1,ymin, 1, yval)
    horizline.DrawLine(xmin,yval,xmax,yval) # draw after legend

    #cuttext.SetNDC();
    #cuttext.DrawLatex(0.65, 0.85, "Fiducial region:")
    #cuttext.DrawLatex(0.65, 0.8,  ptRangeText + ", " + etaRangeText)
    #cuttext.DrawLatex(0.65, 0.75, etaRangeText)
    cuttext.DrawLatex(xtext, ymin + 0.92 * (ymax - ymin), "#bf{Fiducial region:}")
    cuttext.DrawLatex(xtext, ymin + 0.84 * (ymax - ymin),  ptRangeText + ", " + etaRangeText)
    #latMeas.DrawLatex(xtext, ymin + 0.66 * (ymax - ymin), "#bf{Measured}")
    #latMeas.DrawLatex(xtext+0.12, ymin + 0.66 * (ymax - ymin), "#bf{Theory}")

    
    #setTDRStyle()
    latCMS = ROOT.TLatex()
    latCMS.SetNDC();
    latCMS.SetTextFont(42)
    latCMS.SetTextSize(0.045)
    latCMS.DrawLatex(leftMargin, 0.95, '#bf{CMS}' + (' #it{Preliminary}' if not options.skipPreliminary else ''))
    if options.lumiInt > 0 : 
        latCMS.DrawLatex(0.77+(0.04-rightMargin), 0.95, '%s fb^{-1} (13 TeV)' % options.lumiInt)
    else:   
        latCMS.DrawLatex(0.82+(0.04-rightMargin), 0.95, '(13 TeV)' % lumi)
    canvas.RedrawAxis("sameaxis")
    for ext in ["png", "pdf"]:
        canvas.SaveAs(outname + pname + "." + ext)

    textfile.close()

    tfOut.Close()
    print ""
    print "Created file %s" % (outname+outfilename)
    print ""
    f.Close()
