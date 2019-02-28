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


def getTH1fromTH2(h2D,h2Derr=None,unrollAlongX=True):  # unrollAlongX=True --> select rows, i.e. takes a stripe from x1 to xn at same y, then go to next stripe at next y
    nX = h2D.GetNbinsX()
    nY = h2D.GetNbinsY()
    nbins = nX * nY
    name = h2D.GetName() + "_unrollTo1D_" + "eta" if unrollAlongX else "pt" 
    newh = ROOT.TH1D(name,h2D.GetTitle(),nbins,0.5,nbins+0.5)
    if 'TH2' not in h2D.ClassName(): raise RuntimeError, "Calling getTH1fromTH2 on something that is not TH2"
    for i in xrange(nX):
        for j in xrange(nY):
            if unrollAlongX:
                bin = 1 + i + j * nX
            else:
                bin = 1 + j + i * nY
            newh.SetBinContent(bin,h2D.GetBinContent(i+1,j+1))
            if h2Derr: newh.SetBinError(bin,h2Derr.GetBinContent(i+1,j+1))
            else:      newh.SetBinError(bin,h2D.GetBinError(i+1,j+1))
    return newh



if __name__ == "__main__":


    from optparse import OptionParser
    parser = OptionParser(usage='%prog [options]')
    parser.add_option('-i','--input', dest='inputdir', default='', type='string', help='input directory with all the cards inside. It is used to get other information')
    #parser.add_option(     '--no-group-POI', dest='noGroupPOI', default=False , action='store_true', help='Specify that _group_<N>_ is not present in name of POI')
    parser.add_option('-o','--outdir', dest='outdir', default='', type='string', help='output directory to save things')
    parser.add_option('-t','--toyfile', dest='toyfile', default='.', type='string', help='Root file with toys.')
    parser.add_option('-c','--channel', dest='channel', default='', type='string', help='name of the channel (mu or el)')
    parser.add_option('-C','--charge', dest='charge', default='plus,minus', type='string', help='charges to run')
    parser.add_option('-s','--suffix', dest='suffix', default='', type='string', help='Suffix added to output dir (i.e, either to hessian or toys in the path)')
    parser.add_option('-f','--friend', dest='friend', default='', type='string', help='Root file with friend tree containing total xsec (it does not include the outliers). Tree name is assumed to be "toyFriend"')
    parser.add_option('-p','--palette', dest='palette', default='55', type='int', help='Palette for plots')
    parser.add_option('-l','--lumi-norm', dest='lumiNorm', default='-1', type='float', help='If > 0, divide cross section by this factor (lumi in 1/Pb)')
    parser.add_option(     '--lumiInt', dest='lumiInt', default='35.9', type='float', help='Integrated luminosity')
    parser.add_option('-n','--norm-width', dest='normWidth' , default=False , action='store_true',   help='Normalize cross section histograms dividing by bin width')
    parser.add_option(     '--hessian', dest='hessian' , default=False , action='store_true',   help='The file passed with -t is interpreted as hessian: will provide the central value of charge asymmetry but not the uncertainty')
    parser.add_option(     '--fit-data', dest='fitData' , default=False , action='store_true',   help='If True, axis range in plots is customized for data')
    # now doing both
    #parser.add_option(     '--unrollEta', dest='unrollEta' , default=False , action='store_true',   help='If True, make unroll along eta direction instead of pt')
    parser.add_option('-e','--expected-toyfile', dest='exptoyfile', default='', type='string', help='Root file to get expected and make ratio with data (only work with option --fit-data). If SAME, use same file as data and take _gen variables to get the expected')
    parser.add_option(       '--pt-range-bkg', dest='pt_range_bkg', action="append", type="float", nargs=2, default=[], help='Bins with gen level pt in this range are treated as background in the datacard, so there is no POI for them. Takes two float as arguments (increasing order) and can specify multiple times. They should match bin edges and a bin is not considered as background if at least one edge is outside this range (therefore, one can also choose a slightly larger range to avoid floating precision issues).')
    parser.add_option(       '--eta-range-bkg', dest='eta_range_bkg', action="append", type="float", nargs=2, default=[], help='Bins with gen level pt in this range are treated as background in the datacard, so there is no POI for them. Takes two float as arguments (increasing order) and can specify multiple times. They should match bin edges and a bin is not considered as background if at least one edge is outside this range (therefore, one can also choose a slightly larger range to avoid floating precision issues).')
    parser.add_option(     '--pt-range', dest='ptRange', default='template', type='string', help='Comma separated pair of floats used to define the pt range. If "template" is passed, the template min and max values are used.')
    (options, args) = parser.parse_args()

    ROOT.TH1.SetDefaultSumw2()
    ROOT.TH1.StatOverflows(True)
    
    channel = options.channel
    if channel not in ["el","mu"]:
        print "Error: unknown channel %s (select 'el' or 'mu')" % channel
        quit()

    charges = [x for x in options.charge.split(',')]
    for c in charges:
        if c not in ["plus", "minus"]:
            print "Error: unknown charge %s (select 'plus' or 'minus' or both separated by comma)" % c
            quit()



    if options.outdir:
        outname = options.outdir
        addStringToEnd(outname,"/",notAddIfEndswithMatch=True)
        if options.hessian: outname = outname + "/hessian"
        else              : outname = outname + "/toys"
        if len(options.suffix): outname = outname + "_" + options.suffix + "/"
        else                  : outname = outname + "/"
        createPlotDirAndCopyPhp(outname)
    else:
        print "Error: you should specify an output folder using option -o <name>. Exit"
        quit()

    if not options.toyfile:
        print "Error: you should specify a file containing the toys using option -t <name>. Exit"
        quit()


    if options.inputdir:
        if not options.inputdir.endswith("/"): options.inputdir += "/"
        basedir = os.path.basename(os.path.normpath(options.inputdir))
        binfile = options.inputdir + "binningPtEta.txt"
        if "_group" in basedir:
            tokens = basedir.split("_")
            for x in tokens:
                if "group" in x:
                    #ngroup = int("".join(str(s) for s in x if s.isdigit()))
                    ngroup = int(x.split("group")[1])
        else:
            ngroup = 1 
    else:
        print "Error: you should specify an input folder containing all the cards using option -i <name>. Exit"
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


    # if hasPtRangeBkg and not options.ptRange:
    #     for index in range(genBins.Npt):
    #         if ptBinIsBackground[index]:
    #             if firstBin: 

    lepton = "electron" if channel == "el" else "muon"
    Wchannel = "W #rightarrow %s#nu" % ("e" if channel == "el" else "#mu")
    ptRangeText = "p_{{T}}^{{{l}}} #in [{ptmin:3g}, {ptmax:3g}] GeV".format(l="e" if channel == "el" else "#mu", ptmin=genBins.ptBins[0], ptmax=genBins.ptBins[-1])


    hChAsymm = ROOT.TH2F("hChAsymm_{lep}".format(lep=lepton),"Charge asymmetry: {Wch}".format(Wch=Wchannel),
                         genBins.Neta, array('d',genBins.etaBins), genBins.Npt, array('d',genBins.ptBins))
    hChAsymmErr = ROOT.TH2F("hChAsymmErr_{lep}".format(lep=lepton),"Charge asymmetry uncertainty: {Wch}".format(Wch=Wchannel),
                            genBins.Neta, array('d',genBins.etaBins), genBins.Npt, array('d',genBins.ptBins))
    

    hChAsymm1Deta = ROOT.TH1F("hChAsymm1Deta_{lep}".format(lep=lepton),"Charge asymmetry: {Wch}".format(Wch=Wchannel),
                              genBins.Neta, array('d',genBins.etaBins))


    nbins = (genBins.Neta - nEtaBinsBkg) * (genBins.Npt - nPtBinsBkg)
    binCount = 1        
    print "Now filling histograms with charge asymmetry"

    f = ROOT.TFile(options.toyfile, 'read')
    tree = f.Get('fitresults')
    if options.friend != "":
        tree.AddFriend('toyFriend',options.friend)

    for ieta in range(1,genBins.Neta+1):
        for ipt in range(1,genBins.Npt+1):
            if ptBinIsBackground[ipt-1]: continue
            if etaBinIsBackground[ieta-1]: continue
            #print "\rBin {num}/{tot}".format(num=binCount,tot=nbins),
            sys.stdout.write('Bin {num}/{tot}   \r'.format(num=binCount,tot=nbins))
            sys.stdout.flush()
            binCount += 1

            if options.hessian: 
                #central = utilities.getDiffXsecAsymmetryFromHessian(channel,ieta-1,ipt-1,options.toyfile)
                central = utilities.getDiffXsecAsymmetryFromHessianFast(channel,ieta-1,ipt-1,
                                                                        nHistBins=2000, minHist=0., maxHist=1.0, tree=tree)
                error = utilities.getDiffXsecAsymmetryFromHessianFast(channel,ieta-1,ipt-1,
                                                                      nHistBins=2000, minHist=0., maxHist=1.0, tree=tree, getErr=True)
            else:                
                central,up,down = utilities.getDiffXsecAsymmetryFromToysFast(channel,ieta-1,ipt-1,
                                                                             nHistBins=2000, minHist=0., maxHist=1.0, tree=tree)
                #central,up,down = utilities.getDiffXsecAsymmetryFromToys(channel,ieta-1,ipt-1,options.toyfile)
                error = up - central
            hChAsymm.SetBinError(ieta,ipt,error)         
            hChAsymmErr.SetBinContent(ieta,ipt,error)
            hChAsymm.SetBinContent(ieta,ipt,central)        

    for ieta in range(1,genBins.Neta+1):            
        if etaBinIsBackground[ieta-1]: continue
        central = utilities.getDiffXsecAsymmetry1DFromHessianFast(channel,ieta-1,isIeta=True,nHistBins=1000, minHist=0., maxHist=1.0, tree=tree)
        error   = utilities.getDiffXsecAsymmetry1DFromHessianFast(channel,ieta-1,isIeta=True,nHistBins=1000, minHist=0., maxHist=1.0, tree=tree, getErr=True)
        hChAsymm1Deta.SetBinContent(ieta,central)
        hChAsymm1Deta.SetBinError(ieta,error)

    setTDRStyle()

    canvas = ROOT.TCanvas("canvas","",800,700)
    canvas1D = ROOT.TCanvas("canvas1D","",800,800)

    xaxisTitle = 'gen %s |#eta|' % lepton
    yaxisTitle = 'gen %s p_{T} [GeV]' % lepton
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


    #additionalText = "W #rightarrow {lep}#nu::0.2,0.5,0.4,0.6".format(lep="e" if channel == "el" else "#mu") # pass x1,y1,x2,y2
    #additionalText = "W #rightarrow {lep}#nu;{pttext}::0.2,0.6,0.5,0.7".format(lep="e" if channel == "el" else "#mu", pttext=ptRangeText) # pass x1,y1,x2,y2
    texCoord = "0.2,0.65"
    additionalText = "W #rightarrow {lep}#nu;{pttext}::{txc},0.08,0.04".format(lep="e" if channel == "el" else "#mu", 
                                                                               pttext=ptRangeText,
                                                                               txc=texCoord)
    legendCoords = "0.2,0.4,0.75,0.85"
    drawSingleTH1(hChAsymm1Deta,xaxisTitle,"charge asymmetry",
                  "chargeAsym1D_eta_{fl}".format(fl=channel),
                  outname,labelRatioTmp="Rel.Unc.::0.9,1.1",legendCoords=legendCoords, draw_both0_noLog1_onlyLog2=1, drawLineLowerPanel="",
                  passCanvas=canvas1D, lumi=options.lumiInt,
                  lowerPanelHeight=0.35, moreTextLatex=additionalText
                  )

    icharge = 0

    for charge in charges:

        icharge += 1
        print ""
        xaxisTitle = 'gen %s |#eta|' % lepton  # it will be modified later, so I have to restore it here


        chargeSign = "+" if charge == "plus" else "-"

        hMu = ROOT.TH2F("hMu_{lep}_{ch}".format(lep=lepton,ch=charge),
                              "signal strength: {Wch}".format(Wch=Wchannel.replace('W','W{chs}'.format(chs=chargeSign))),
                              genBins.Neta, array('d',genBins.etaBins), genBins.Npt, array('d',genBins.ptBins))
        hMuErr = ROOT.TH2F("hMuErr_{lep}_{ch}".format(lep=lepton,ch=charge),
                                 "signal strength uncertainty: {Wch}".format(Wch=Wchannel.replace('W','W{chs}'.format(chs=chargeSign))),
                                 genBins.Neta, array('d',genBins.etaBins), genBins.Npt, array('d',genBins.ptBins))

        hDiffXsec = ROOT.TH2F("hDiffXsec_{lep}_{ch}".format(lep=lepton,ch=charge),
                              "cross section: {Wch}".format(Wch=Wchannel.replace('W','W{chs}'.format(chs=chargeSign))),
                              genBins.Neta, array('d',genBins.etaBins), genBins.Npt, array('d',genBins.ptBins))
        hDiffXsecErr = ROOT.TH2F("hDiffXsecErr_{lep}_{ch}".format(lep=lepton,ch=charge),
                                 "cross section uncertainty: {Wch}".format(Wch=Wchannel.replace('W','W{chs}'.format(chs=chargeSign))),
                                 genBins.Neta, array('d',genBins.etaBins), genBins.Npt, array('d',genBins.ptBins))
        hDiffXsecNorm = ROOT.TH2F("hDiffXsecNorm_{lep}_{ch}".format(lep=lepton,ch=charge),
                                  "normalized cross section: {Wch}".format(Wch=Wchannel.replace('W','W{chs}'.format(chs=chargeSign))),
                                  genBins.Neta, array('d',genBins.etaBins), genBins.Npt, array('d',genBins.ptBins))
        hDiffXsecNormErr = ROOT.TH2F("hDiffXsecNormErr_{lep}_{ch}".format(lep=lepton,ch=charge),
                                     "normalized cross section uncertainty: {Wch}".format(Wch=Wchannel.replace('W','W{chs}'.format(chs=chargeSign))),
                                     genBins.Neta, array('d',genBins.etaBins), genBins.Npt, array('d',genBins.ptBins))

        hDiffXsec_1Deta = ROOT.TH1F("hDiffXsec_1Deta_{lep}_{ch}".format(lep=lepton,ch=charge),
                                    "cross section: {Wch}".format(Wch=Wchannel.replace('W','W{chs}'.format(chs=chargeSign))),
                                    genBins.Neta, array('d',genBins.etaBins))

        hDiffXsecNorm_1Deta = ROOT.TH1F("hDiffXsecNorm_1Deta_{lep}_{ch}".format(lep=lepton,ch=charge),
                                        "cross section: {Wch}".format(Wch=Wchannel.replace('W','W{chs}'.format(chs=chargeSign))),
                                        genBins.Neta, array('d',genBins.etaBins))


        #nbins = genBins.Neta * (genBins.Npt - nPtBinsBkg)
        binCount = 1        
        print "Now filling histograms with differential cross section (charge = %s)" % charge


        denExpression = ""
        
        if not options.hessian:
            if options.friend != "":
                denExpression = "totxsec_" + charge
            else:
                denExpression = utilities.getDenExpressionForNormDiffXsec(channel, charge, genBins.Neta,genBins.Npt)

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
                    central = utilities.getSignalStrengthFromHessianFast(channel,charge,ieta-1,ipt-1,
                                                                         nHistBins=100, minHist=0.5, maxHist=1.5, tree=tree)                    
                    error = utilities.getSignalStrengthFromHessianFast(channel,charge,ieta-1,ipt-1,
                                                                       nHistBins=200, minHist=0.0, maxHist=1., tree=tree, getErr=True)                    
                else:
                    central,up,down = utilities.getSignalStrengthFromToysFast(channel,charge,ieta-1,ipt-1, 
                                                                              nHistBins=200, minHist=0.5, maxHist=1.5, tree=tree)
                    error = up - central
                hMu.SetBinContent(ieta,ipt,central)        
                hMu.SetBinError(ieta,ipt,error)         
                hMuErr.SetBinContent(ieta,ipt,error)
                
                # normalized cross section
                if options.hessian:                    
                    central = utilities.getNormalizedDiffXsecFromHessianFast(channel,charge,ieta-1,ipt-1,
                                                                             nHistBins=2000, minHist=0., maxHist=200., tree=tree)                    
                    error = utilities.getNormalizedDiffXsecFromHessianFast(channel,charge,ieta-1,ipt-1,
                                                                           nHistBins=2000, minHist=0., maxHist=10., tree=tree, getErr=True)                    
                else:
                    central,up,down = utilities.getNormalizedDiffXsecFromToysFast(channel,charge,ieta-1,ipt-1,
                                                                                  denExpression, nHistBins=1000, minHist=0., maxHist=0.1, tree=tree)
                    error = up - central
                hDiffXsecNorm.SetBinContent(ieta,ipt,central)        
                hDiffXsecNorm.SetBinError(ieta,ipt,error)         
                hDiffXsecNormErr.SetBinContent(ieta,ipt,error)
                
                # unnormalized cross section
                if options.hessian:
                    central = utilities.getDiffXsecFromHessianFast(channel,charge,ieta-1,ipt-1,
                                                                   nHistBins=2000, minHist=0., maxHist=200., tree=tree)                    
                    error = utilities.getDiffXsecFromHessianFast(channel,charge,ieta-1,ipt-1,
                                                                 nHistBins=2000, minHist=0., maxHist=200., tree=tree, getErr=True)                    
                else:                
                    central,up,down = utilities.getDiffXsecFromToysFast(channel,charge,ieta-1,ipt-1,
                                                                        nHistBins=2000, minHist=0., maxHist=200., tree=tree)
                    error = up - central
                hDiffXsec.SetBinContent(ieta,ipt,central)        
                hDiffXsec.SetBinError(ieta,ipt,error)         
                hDiffXsecErr.SetBinContent(ieta,ipt,error)
                
        for ieta in range(1,genBins.Neta+1):
            if etaBinIsBackground[ieta-1]: continue
            central = utilities.getDiffXsec1DFromHessianFast(channel,charge,ieta-1,isIeta=True,nHistBins=5000, minHist=0., maxHist=5000., tree=tree)
            error   = utilities.getDiffXsec1DFromHessianFast(channel,charge,ieta-1,isIeta=True,nHistBins=5000, minHist=0., maxHist=5000., tree=tree, getErr=True)
            hDiffXsec_1Deta.SetBinContent(ieta, central)
            hDiffXsec_1Deta.SetBinError(ieta, error)

            central = utilities.getNormalizedDiffXsec1DFromHessianFast(channel,charge,ieta-1,isIeta=True,nHistBins=1000, minHist=0., maxHist=1., tree=tree)
            error   = utilities.getNormalizedDiffXsec1DFromHessianFast(channel,charge,ieta-1,isIeta=True,nHistBins=1000, minHist=0., maxHist=1., tree=tree, getErr=True)
            hDiffXsecNorm_1Deta.SetBinContent(ieta, central)
            hDiffXsecNorm_1Deta.SetBinError(ieta, error)


        if options.normWidth:
            hDiffXsec.Scale(1.,"width")
            hDiffXsecErr.Scale(1.,"width")
            hDiffXsecNorm.Scale(1.,"width")
            hDiffXsecNormErr.Scale(1.,"width")
            hDiffXsec_1Deta.Scale(1.,"width")
            hDiffXsecNorm_1Deta.Scale(1.,"width")

        if options.lumiNorm > 0:
            scaleFactor = 1./options.lumiNorm
            hDiffXsec.Scale(scaleFactor)
            hDiffXsecErr.Scale(scaleFactor)
            # do not divide the normalized cross section, the scaling factor is already removed
            #hDiffXsecNorm.Scale(scaleFactor)
            #hDiffXsecNormErr.Scale(scaleFactor)
            hDiffXsec_1Deta.Scale(scaleFactor)

        hDiffXsecRelErr = hDiffXsecErr.Clone(hDiffXsecErr.GetName().replace('XsecErr','XsecRelErr'))
        hDiffXsecRelErr.Divide(hDiffXsec)
        hDiffXsecRelErr.SetTitle(hDiffXsecErr.GetTitle().replace('uncertainty','rel.unc.'))

        hDiffXsecNormRelErr = hDiffXsecNormErr.Clone(hDiffXsecNormErr.GetName().replace('XsecNormErr','XsecNormRelErr'))
        hDiffXsecNormRelErr.Divide(hDiffXsecNorm)
        hDiffXsecNormRelErr.SetTitle(hDiffXsecNormErr.GetTitle().replace('uncertainty','rel.unc.'))

        # now starting to draw the cross sections

        zaxisTitle = "Signal strength #mu::0.98,1.02"
        if options.hessian: zaxisTitle = "Signal strength #mu::0.995,1.005"
        if options.fitData: zaxisTitle = "Signal strength #mu::0.5,1.5"
        drawCorrelationPlot(hMu,
                            xaxisTitle, yaxisTitle, zaxisTitle,
                            hMu.GetName(),
                            "ForceTitle",outname,1,1,False,False,False,1, canvasSize="700,625",leftMargin=0.14,rightMargin=0.22,passCanvas=canvas,palette=options.palette)

        zaxisTitle = "uncertainty on signal strength #mu::0.0,0.5"
        if options.hessian: zaxisTitle = "uncertainty on signal strength #mu::0.0,0.25"
        drawCorrelationPlot(hMuErr,
                            xaxisTitle, yaxisTitle, zaxisTitle,
                            hMuErr.GetName(),
                            "ForceTitle",outname,1,1,False,False,False,1, canvasSize="700,625",leftMargin=0.14,rightMargin=0.22,passCanvas=canvas,palette=options.palette)

        if options.fitData:
            if charge == "plus": zmin,zmax = 30,130
            else:                zmin,zmax = 10,110
        else:
            if charge == "plus": zmin,zmax = 30,120
            else:                zmin,zmax = 25,95

        zminHist = hDiffXsec.GetBinContent(hDiffXsec.GetMinimumBin())                                                   
        zaxisTitle = "d^{2}#sigma / d|#eta|dp_{T} [pb/GeV]::%.3f,%.3f" % (0.99*(zminHist if zminHist > 1. else zmin),
                                                                    min(zmax,hDiffXsec.GetBinContent(hDiffXsec.GetMaximumBin())))
        #zaxisTitle = "d^{2}#sigma / d|#eta|dp_{T} [pb/GeV]::%.3f,%.3f" % (zmin,zmax)
        drawCorrelationPlot(hDiffXsec,
                            xaxisTitle, yaxisTitle, zaxisTitle,
                            hDiffXsec.GetName(),
                            "ForceTitle",outname,1,1,False,False,False,1, canvasSize="700,625",leftMargin=0.14,rightMargin=0.22,passCanvas=canvas,palette=options.palette)

        if options.fitData:
            if charge == "plus": zmin,zmax = 0.5,4.5
            else:                zmin,zmax = 0.5,3.5
        else:
            if charge == "plus": zmin,zmax = 0.5,4.5
            else:                zmin,zmax = 0.5,3.5
        zaxisTitle = "uncertainty on d^{2}#sigma / d|#eta|dp_{T} [pb/GeV]::%.3f,%.3f" % (0.9*hDiffXsecErr.GetMinimum(),min(25,hDiffXsecErr.GetMaximum()))
        #zaxisTitle = "uncertainty on d^{2}#sigma / d|#eta|dp_{T} [pb/GeV]::%.3f,%.3f" % (zmin,zmax)
        drawCorrelationPlot(hDiffXsecErr,
                            xaxisTitle, yaxisTitle, zaxisTitle,
                            hDiffXsecErr.GetName(),
                            "ForceTitle",outname,1,1,False,False,False,1, canvasSize="700,625",leftMargin=0.14,rightMargin=0.22,passCanvas=canvas,palette=options.palette)

        #zaxisTitle = "rel. uncertainty on d^{2}#sigma / d|#eta|dp_{T}::%.3f,%.3f" % (0.9*hDiffXsecRelErr.GetMinimum(),hDiffXsecRelErr.GetMaximum())
        if options.fitData:
            zaxisTitle = "rel. uncertainty on d^{2}#sigma / d|#eta|dp_{T}::0.010,0.2"
        else:
            zaxisTitle = "rel. uncertainty on d^{2}#sigma / d|#eta|dp_{T}::0.025,0.2"
        drawCorrelationPlot(hDiffXsecRelErr,
                            xaxisTitle, yaxisTitle, zaxisTitle,
                            hDiffXsecRelErr.GetName(),
                            "ForceTitle",outname,1,1,False,False,False,1, canvasSize="700,625",leftMargin=0.14,rightMargin=0.22,passCanvas=canvas,palette=options.palette)

        if charge == "plus": zmin,zmax = 0.008,0.028
        else:                zmin,zmax = 0.008,0.029
        #zaxisTitle = "d^{2}#sigma / d|#eta|dp_{T} / #sigma_{tot}::%.3f,%.3f" % (0.9*hDiffXsecNorm.GetMinimum(),hDiffXsecNorm.GetMaximum())
        zaxisTitle = "d^{2}#sigma / d|#eta|dp_{T} / #sigma_{tot}::%.3f,%.3f" % (zmin,zmax)
        drawCorrelationPlot(hDiffXsecNorm,
                            xaxisTitle, yaxisTitle, zaxisTitle,
                            hDiffXsecNorm.GetName(),
                            "ForceTitle",outname,1,1,False,False,False,1, canvasSize="700,625",leftMargin=0.14,rightMargin=0.22,passCanvas=canvas,palette=options.palette)

        if charge == "plus": zmin,zmax = 0.0,0.0015
        else:                zmin,zmax = 0.0,0.0015
        #zaxisTitle = "uncertainty on d^{2}#sigma / d|#eta|dp_{T} / #sigma_{tot}::%.3f,%.3f" % (0.9*hDiffXsecNormErr.GetMinimum(),hDiffXsecNormErr.GetMaximum())
        zaxisTitle = "uncertainty on d^{2}#sigma / d|#eta|dp_{T} / #sigma_{tot}::%.3f,%.3f" % (zmin,zmax)
        drawCorrelationPlot(hDiffXsecNormErr,
                            xaxisTitle, yaxisTitle, zaxisTitle,
                            hDiffXsecNormErr.GetName(),
                            "ForceTitle",outname,1,1,False,False,False,1, canvasSize="700,625",leftMargin=0.14,rightMargin=0.22,passCanvas=canvas,palette=options.palette)

        zaxisTitle = "rel. uncertainty on d^{2}#sigma / d|#eta|dp_{T} / #sigma_{tot}::%.4f,%.4f" % (min(0.01,0.9*hDiffXsecNormRelErr.GetMinimum()),min(0.2,hDiffXsecNormRelErr.GetBinContent(hDiffXsecNormRelErr.GetMaximumBin())))
        #zaxisTitle = "rel. uncertainty on d^{2}#sigma / d|#eta|dp_{T} / #sigma_{tot}::0,0.1"
        drawCorrelationPlot(hDiffXsecNormRelErr,
                            xaxisTitle, yaxisTitle, zaxisTitle,
                            hDiffXsecNormRelErr.GetName(),
                            "ForceTitle",outname,1,1,False,False,False,1, canvasSize="700,625",leftMargin=0.14,rightMargin=0.22,passCanvas=canvas,palette=options.palette)


        # with only W -> munu, coordinates can be 0.45,0.8,0.65,0.9 with TPaveText
        texCoord = "0.45,0.85" if charge == "plus" else "0.45,0.5"
        additionalText = "W^{{{chs}}} #rightarrow {lep}#nu;{pttext}::{txc},0.08,0.04".format(chs=" "+chargeSign,lep="e" if channel == "el" else "#mu",
                                                                                             pttext=ptRangeText,
                                                                                             txc=texCoord) 
        legendCoords = "0.2,0.4,0.75,0.85" if charge == "plus" else "0.2,0.4,0.4,0.5"
        drawSingleTH1(hDiffXsec_1Deta,xaxisTitle,"d#sigma/d|#eta| [pb]",
                      "xsec_eta_abs_{ch}_{fl}".format(ch=charge,fl=channel),
                      outname,labelRatioTmp="Rel.Unc.::0.9,1.1",legendCoords=legendCoords, draw_both0_noLog1_onlyLog2=1,
                      passCanvas=canvas1D, lumi=options.lumiInt,
                      lowerPanelHeight=0.35, moreTextLatex=additionalText
                      )
        drawSingleTH1(hDiffXsecNorm_1Deta,xaxisTitle,"d#sigma/d|#eta| / #sigma_{tot}",
                      "xsec_eta_norm_{ch}_{fl}".format(ch=charge,fl=channel),
                      outname,labelRatioTmp="Rel.Unc.::0.9,1.1",legendCoords=legendCoords, draw_both0_noLog1_onlyLog2=1,
                      passCanvas=canvas1D, lumi=options.lumiInt,
                      lowerPanelHeight=0.35, moreTextLatex=additionalText
                      )



######
        # now drawing a TH1 unrolling TH2
        chsize = 4500
        canvUnroll = ROOT.TCanvas("canvUnroll","",chsize,1500) # 3000,1500
        leftMargin = 0.12 if chsize < 3500 else 0.06
        rightMargin = 0.04 if chsize < 3500 else 0.02

        ratioYaxis = "Rel.Unc.::0.8,1.2"
        if channel == "el": ratioYaxis = "Rel.Unc.::0.8,1.2"

        xaxisTitle = "template global bin"
        vertLinesArg = ""
        additionalText = "W^{{{chs}}} #rightarrow {lep}#nu::0.8,0.84,0.9,0.9".format(chs=" "+chargeSign, lep="e" if channel == "el" else "#mu") # pass x1,y1,x2,y2

        h1D_chargeAsym = {}
        h1D_pmaskedexp = {}
        h1D_pmaskedexp_norm = {}

        for unrollAlongEta in [True, False]:

            unrollVar = "eta" if unrollAlongEta else "pt"

            varBinRanges = []

            if unrollAlongEta:
                #xaxisTitle = xaxisTitle + " = 1 + ipt + ieta * %d; ipt in [%d,%d], ieta in [%d,%d]" % (nptbins-1,0,nptbins-1,0,netabins-1)
                xaxisTitle = "cross section unrolled along |#eta|: |#eta| #in [%.1f, %.1f]" % (genBins.etaBins[0], genBins.etaBins[-1])
                vertLinesArg = "{a},{b}".format(a=genBins.Npt,b=genBins.Neta)
                for ipt in range(0,genBins.Npt):
                    varBinRanges.append("p_{{T}} #in [{ptmin:3g}, {ptmax:.3g}]".format(ptmin=genBins.ptBins[ipt], ptmax=genBins.ptBins[ipt+1]))
            else:
                #xaxisTitle = xaxisTitle + " = 1 + ieta + ipt * %d; ipt in [%d,%d], ieta in [%d,%d]" % (netabins-1,0,nptbins-1,0,netabins-1)
                xaxisTitle = "cross section unrolled along p_{T}: p_{T} #in [%.3g, %.3g] GeV" % (genBins.ptBins[0], genBins.ptBins[-1])
                vertLinesArg = "{a},{b}".format(a=genBins.Neta,b=genBins.Npt)
                for ieta in range(0,genBins.Neta):
                    varBinRanges.append("|#eta| #in [{etamin:.1f}, {etamax:.1f}]".format(etamin=genBins.etaBins[ieta], etamax=genBins.etaBins[ieta+1]))

            h1D_pmaskedexp[unrollAlongEta] = getTH1fromTH2(hDiffXsec, hDiffXsecErr, unrollAlongX=unrollAlongEta)        
            drawSingleTH1(h1D_pmaskedexp[unrollAlongEta],xaxisTitle,"d^{2}#sigma/d|#eta|dp_{T} [pb/GeV]",
                          "unrolledXsec_{var}_abs_{ch}_{fl}".format(var=unrollVar,ch=charge,fl=channel),
                          outname,labelRatioTmp=ratioYaxis,draw_both0_noLog1_onlyLog2=1,passCanvas=canvUnroll, lumi=options.lumiInt,
                          drawVertLines=vertLinesArg, textForLines=varBinRanges, lowerPanelHeight=0.35, moreText=additionalText,
                          leftMargin=leftMargin, rightMargin=rightMargin)

            h1D_pmaskedexp_norm[unrollAlongEta] = getTH1fromTH2(hDiffXsecNorm, hDiffXsecNormErr, unrollAlongX=unrollAlongEta)        
            drawSingleTH1(h1D_pmaskedexp_norm[unrollAlongEta],xaxisTitle,"d^{2}#sigma/d|#eta|dp_{T} / #sigma_{tot} [1/GeV]",
                          "unrolledXsec_{var}_norm_{ch}_{fl}".format(var=unrollVar,ch=charge,fl=channel),
                          outname,labelRatioTmp=ratioYaxis,draw_both0_noLog1_onlyLog2=1,passCanvas=canvUnroll, lumi=options.lumiInt, drawLineLowerPanel="",
                          drawVertLines=vertLinesArg, textForLines=varBinRanges, lowerPanelHeight=0.35,moreText=additionalText,
                          leftMargin=leftMargin, rightMargin=rightMargin)

            # do also charge asymmetry (only once if running on both charges)
            if icharge == 1:

                xaxisTitle = xaxisTitle.replace("cross section", "charge asymmetry")
                additionalTextBackup = additionalText
                additionalText = "W #rightarrow {lep}#nu::0.8,0.84,0.9,0.9".format(lep="e" if channel == "el" else "#mu") # pass x1,y1,x2,y2

                h1D_chargeAsym[unrollAlongEta] = getTH1fromTH2(hChAsymm, h2Derr=None, unrollAlongX=unrollAlongEta)        
                drawSingleTH1(h1D_chargeAsym[unrollAlongEta], xaxisTitle,"charge asymmetry::0.0,1.0",
                              "unrolledChargeAsym_{var}_{fl}".format(var=unrollVar,fl=channel),
                              outname,labelRatioTmp=ratioYaxis,draw_both0_noLog1_onlyLog2=1,drawLineLowerPanel="", 
                              passCanvas=canvUnroll,lumi=options.lumiInt,
                              drawVertLines=vertLinesArg, textForLines=varBinRanges, lowerPanelHeight=0.35, moreText=additionalText,
                              leftMargin=leftMargin, rightMargin=rightMargin)

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

                hChargeAsym_exp = None
                if icharge == 1:
                    hChargeAsym_exp = ROOT.TH2F("hChargeAsym_{lep}_exp".format(lep=lepton),
                                                "charge asymmetry: {Wch}".format(Wch=Wchannel),
                                                genBins.Neta, array('d',genBins.etaBins), genBins.Npt, array('d',genBins.ptBins))

                    hChAsymm1Deta_exp = ROOT.TH1F("hChAsymm1Deta_{lep}_exp".format(lep=lepton),"Charge asymmetry: {Wch}".format(Wch=Wchannel),
                                                  genBins.Neta, array('d',genBins.etaBins))


                hDiffXsec_exp = ROOT.TH2F("hDiffXsec_{lep}_{ch}_exp".format(lep=lepton,ch=charge),
                                          "cross section: {Wch}".format(Wch=Wchannel.replace('W','W{chs}'.format(chs=chargeSign))),
                                          genBins.Neta, array('d',genBins.etaBins), genBins.Npt, array('d',genBins.ptBins))
                #hDiffXsecErr = ROOT.TH2F("hDiffXsecErr_{lep}_{ch}".format(lep=lepton,ch=charge),
                #                         "cross section uncertainty: {Wch}".format(Wch=Wchannel.replace('W','W{chs}'.format(chs=chargeSign))),
                #                         genBins.Neta, array('d',genBins.etaBins), genBins.Npt, array('d',genBins.ptBins))
                hDiffXsecNorm_exp = ROOT.TH2F("hDiffXsecNorm_{lep}_{ch}_exp".format(lep=lepton,ch=charge),
                                              "normalized cross section: {Wch}".format(Wch=Wchannel.replace('W','W{chs}'.format(chs=chargeSign))),
                                              genBins.Neta, array('d',genBins.etaBins), genBins.Npt, array('d',genBins.ptBins))
                #hDiffXsecNormErr = ROOT.TH2F("hDiffXsecNormErr_{lep}_{ch}".format(lep=lepton,ch=charge),
                #                             "normalized cross section uncertainty: {Wch}".format(Wch=Wchannel.replace('W','W{chs}'.format(chs=chargeSign))),
                #                             genBins.Neta, array('d',genBins.etaBins), genBins.Npt, array('d',genBins.ptBins))


                hDiffXsec_1Deta_exp = ROOT.TH1F("hDiffXsec_1Deta_{lep}_{ch}_exp".format(lep=lepton,ch=charge),
                                                "cross section: {Wch}".format(Wch=Wchannel.replace('W','W{chs}'.format(chs=chargeSign))),
                                                genBins.Neta, array('d',genBins.etaBins))

                hDiffXsecNorm_1Deta_exp = ROOT.TH1F("hDiffXsecNorm_1Deta_{lep}_{ch}_exp".format(lep=lepton,ch=charge),
                                                    "cross section: {Wch}".format(Wch=Wchannel.replace('W','W{chs}'.format(chs=chargeSign))),
                                                    genBins.Neta, array('d',genBins.etaBins))


                binCount = 0
                print ""
                print "Now reading Hessian to make ratio of data with expected"
                for ieta in range(1,genBins.Neta+1):
                    for ipt in range(1,genBins.Npt+1):
                        if ptBinIsBackground[ipt-1]: continue
                        if etaBinIsBackground[ieta-1]: continue
                        #print "\rBin {num}/{tot}".format(num=binCount,tot=nbins),
                        sys.stdout.write('Bin {num}/{tot}   \r'.format(num=binCount,tot=nbins))
                        sys.stdout.flush()
                        binCount += 1
                        
                        if icharge == 1:
                            # charge asymmetry
                            central = utilities.getDiffXsecAsymmetryFromHessianFast(channel,ieta-1,ipt-1,
                                                                                   nHistBins=2000, minHist=0., maxHist=1., 
                                                                                   tree=treeexp, getGen=getExpFromGen)         
                            error = utilities.getDiffXsecAsymmetryFromHessianFast(channel,ieta-1,ipt-1,
                                                                                 nHistBins=2000, minHist=0., maxHist=1., 
                                                                                 tree=treeexp, getErr=True)                    
                            
                            hChargeAsym_exp.SetBinContent(ieta,ipt,central)        
                            hChargeAsym_exp.SetBinError(ieta,ipt,error)         


                        # normalized cross section
                        central = utilities.getNormalizedDiffXsecFromHessianFast(channel,charge,ieta-1,ipt-1,
                                                                                 nHistBins=2000, minHist=0., maxHist=200., tree=treeexp, getGen=getExpFromGen)         
                        error = utilities.getNormalizedDiffXsecFromHessianFast(channel,charge,ieta-1,ipt-1,
                                                                               nHistBins=2000, minHist=0., maxHist=200., tree=treeexp, getErr=True)                    

                        hDiffXsecNorm_exp.SetBinContent(ieta,ipt,central)        
                        hDiffXsecNorm_exp.SetBinError(ieta,ipt,error)         

                        # unnormalized cross section
                        central = utilities.getDiffXsecFromHessianFast(channel,charge,ieta-1,ipt-1,
                                                                       nHistBins=2000, minHist=0., maxHist=200., tree=treeexp, getGen=getExpFromGen)                    
                        error = utilities.getDiffXsecFromHessianFast(channel,charge,ieta-1,ipt-1,
                                                                     nHistBins=2000, minHist=0., maxHist=200., tree=treeexp, getErr=True)                    
                        hDiffXsec_exp.SetBinContent(ieta,ipt,central)        
                        hDiffXsec_exp.SetBinError(ieta,ipt,error)         

                # now 1D stuff
                for ieta in range(1,genBins.Neta+1):

                    if etaBinIsBackground[ieta-1]: continue

                    if icharge == 1:
                        central = utilities.getDiffXsecAsymmetry1DFromHessianFast(channel,ieta-1,isIeta=True,nHistBins=1000, minHist=0., maxHist=1.0, tree=treeexp)
                        error   = utilities.getDiffXsecAsymmetry1DFromHessianFast(channel,ieta-1,isIeta=True,nHistBins=1000, minHist=0., maxHist=1.0, tree=treeexp, getErr=True)
                        hChAsymm1Deta_exp.SetBinContent(ieta,central)
                        hChAsymm1Deta_exp.SetBinError(ieta,error)
        
                    central = utilities.getDiffXsec1DFromHessianFast(channel,charge,ieta-1,isIeta=True,nHistBins=5000, minHist=0., maxHist=5000., tree=treeexp)
                    error   = utilities.getDiffXsec1DFromHessianFast(channel,charge,ieta-1,isIeta=True,nHistBins=5000, minHist=0., maxHist=5000., tree=treeexp, getErr=True)
                    hDiffXsec_1Deta_exp.SetBinContent(ieta, central)
                    hDiffXsec_1Deta_exp.SetBinError(ieta, error)

                    central = utilities.getNormalizedDiffXsec1DFromHessianFast(channel,charge,ieta-1,isIeta=True,nHistBins=1000, minHist=0., maxHist=1., tree=treeexp)
                    error   = utilities.getNormalizedDiffXsec1DFromHessianFast(channel,charge,ieta-1,isIeta=True,nHistBins=1000, minHist=0., maxHist=1., tree=treeexp, getErr=True)
                    hDiffXsecNorm_1Deta_exp.SetBinContent(ieta, central)
                    hDiffXsecNorm_1Deta_exp.SetBinError(ieta, error)


                if options.normWidth:
                    hDiffXsec_exp.Scale(1.,"width")
                    hDiffXsecNorm_exp.Scale(1.,"width")
                    hDiffXsec_1Deta_exp.Scale(1.,"width")
                    hDiffXsecNorm_1Deta_exp.Scale(1.,"width")
                    

                if options.lumiNorm > 0:
                    scaleFactor = 1./options.lumiNorm
                    hDiffXsec_exp.Scale(scaleFactor)
                    hDiffXsec_1Deta_exp.Scale(scaleFactor)

                # # now drawing a TH1 unrolling TH2

                for unrollAlongEta in [True, False]:

                    unrollVar = "eta" if unrollAlongEta else "pt"

                    varBinRanges = []

                    if unrollAlongEta:
                        #xaxisTitle = xaxisTitle + " = 1 + ipt + ieta * %d; ipt in [%d,%d], ieta in [%d,%d]" % (nptbins-1,0,nptbins-1,0,netabins-1)
                        xaxisTitle = "cross section unrolled along |#eta|: |#eta| #in [%.1f, %.1f]" % (genBins.etaBins[0], genBins.etaBins[-1])
                        vertLinesArg = "{a},{b}".format(a=genBins.Npt,b=genBins.Neta)
                        for ipt in range(0,genBins.Npt):
                            varBinRanges.append("p_{{T}} #in [{ptmin:3g}, {ptmax:.3g}]".format(ptmin=genBins.ptBins[ipt], ptmax=genBins.ptBins[ipt+1]))
                    else:
                        #xaxisTitle = xaxisTitle + " = 1 + ieta + ipt * %d; ipt in [%d,%d], ieta in [%d,%d]" % (netabins-1,0,nptbins-1,0,netabins-1)
                        xaxisTitle = "cross section unrolled along p_{T}: p_{T} #in [%.3g, %.3g] GeV" % (genBins.ptBins[0], genBins.ptBins[-1])
                        vertLinesArg = "{a},{b}".format(a=genBins.Neta,b=genBins.Npt)
                        for ieta in range(0,genBins.Neta):
                            varBinRanges.append("|#eta| #in [{etamin:.1f}, {etamax:.1f}]".format(etamin=genBins.etaBins[ieta], etamax=genBins.etaBins[ieta+1]))

                    h1D_pmaskedexp_exp = getTH1fromTH2(hDiffXsec_exp, h2Derr=None, unrollAlongX=unrollAlongEta)        
                    drawDataAndMC(h1D_pmaskedexp[unrollAlongEta], h1D_pmaskedexp_exp,xaxisTitle,"d^{2}#sigma/d|#eta|dp_{T} [pb/GeV]::0,220",
                                  "unrolledXsec_{var}_abs_{ch}_{fl}_dataAndExp".format(var=unrollVar,ch=charge,fl=channel),
                                  outname,labelRatioTmp="obs./exp.::0.8,1.2",draw_both0_noLog1_onlyLog2=1,passCanvas=canvUnroll,lumi=options.lumiInt,
                                  drawVertLines=vertLinesArg, textForLines=varBinRanges, lowerPanelHeight=0.35, moreText=additionalText,
                                  leftMargin=leftMargin, rightMargin=rightMargin)


                    h1D_pmaskedexp_norm_exp = getTH1fromTH2(hDiffXsecNorm_exp, h2Derr=None, unrollAlongX=unrollAlongEta)        
                    drawDataAndMC(h1D_pmaskedexp_norm[unrollAlongEta], h1D_pmaskedexp_norm_exp,
                                  xaxisTitle,"d^{2}#sigma/d|#eta|dp_{T} / #sigma_{tot} [1/GeV]",
                                  "unrolledXsec_{var}_norm_{ch}_{fl}_dataAndExp".format(var=unrollVar,ch=charge,fl=channel),
                                  outname,labelRatioTmp="obs./exp.::0.8,1.2",draw_both0_noLog1_onlyLog2=1,passCanvas=canvUnroll, lumi=options.lumiInt,
                                  drawVertLines=vertLinesArg, textForLines=varBinRanges, lowerPanelHeight=0.35, moreText=additionalText,
                                  leftMargin=leftMargin, rightMargin=rightMargin)


                    # do also charge asymmetry (only once if running on both charges)
                    if icharge == 1:
                        
                        additionalTextBackup = additionalText
                        additionalText = "W #rightarrow {lep}#nu::0.8,0.84,0.9,0.9".format(lep="e" if channel == "el" else "#mu") # pass x1,y1,x2,y2

                        h1D_chargeAsym_exp = getTH1fromTH2(hChargeAsym_exp, h2Derr=None, unrollAlongX=unrollAlongEta)        
                        drawDataAndMC(h1D_chargeAsym[unrollAlongEta], h1D_chargeAsym_exp,xaxisTitle.replace("cross section", "charge asymmetry"),
                                      "charge asymmetry::0,1.0", "unrolledChargeAsym_{var}_{fl}_dataAndExp".format(var=unrollVar,fl=channel),
                                      outname,labelRatioTmp="obs./exp.::0.8,1.2",draw_both0_noLog1_onlyLog2=1,passCanvas=canvUnroll,lumi=options.lumiInt,
                                      drawVertLines=vertLinesArg, textForLines=varBinRanges, lowerPanelHeight=0.35, moreText=additionalText,
                                      leftMargin=leftMargin, rightMargin=rightMargin)

                        additionalText = additionalTextBackup

                # now data/prediction for the 1D xsection and charge asymmetry
                ###############
                xaxisTitle = 'gen %s |#eta|' % lepton
                additionalTextBackup = additionalText
                #additionalText = "W^{{{chs}}} #rightarrow {lep}#nu;{pttext}::0.45,0.8,0.75,0.9".format(chs=" "+chargeSign,lep="e" if channel == "el" else "#mu",
                #                                                                                       pttext=ptRangeText) # pass x1,y1,x2,y2
                texCoord = "0.45,0.85" if charge == "plus" else "0.45,0.5"
                additionalText = "W^{{{chs}}} #rightarrow {lep}#nu;{pttext}::{txc},0.08,0.04".format(chs=" "+chargeSign,lep="e" if channel == "el" else "#mu",
                                                                                                     pttext=ptRangeText,
                                                                                                     txc=texCoord) 
                legendCoords = "0.2,0.4,0.75,0.85" if charge == "plus" else "0.2,0.4,0.4,0.5"
                drawDataAndMC(hDiffXsec_1Deta,hDiffXsec_1Deta_exp,xaxisTitle,"d#sigma/d|#eta| [pb]",
                              "xsec_eta_abs_{ch}_{fl}_dataAndExp".format(ch=charge,fl=channel),
                              outname,labelRatioTmp="obs./exp.::0.8,1.2",legendCoords=legendCoords, draw_both0_noLog1_onlyLog2=1,
                              passCanvas=canvas1D, lumi=options.lumiInt,
                              lowerPanelHeight=0.35, moreTextLatex=additionalText
                              )
                drawDataAndMC(hDiffXsecNorm_1Deta,hDiffXsecNorm_1Deta_exp,xaxisTitle,"d#sigma/d|#eta| / #sigma_{tot}",
                              "xsec_eta_norm_{ch}_{fl}_dataAndExp".format(ch=charge,fl=channel),
                              outname,labelRatioTmp="obs./exp.::0.8,1.2",legendCoords=legendCoords, draw_both0_noLog1_onlyLog2=1,
                              passCanvas=canvas1D, lumi=options.lumiInt,
                              lowerPanelHeight=0.35, moreTextLatex=additionalText
                              )
                additionalText = additionalTextBackup

                if icharge == 1:
                                        
                    additionalTextBackup = additionalText
                    #additionalText = "W #rightarrow {lep}#nu;{pttext}::0.2,0.6,0.5,0.7".format(lep="e" if channel == "el" else "#mu",pttext=ptRangeText) # pass x1,y1,x2,y2
                    texCoord = "0.2,0.65"
                    additionalText = "W #rightarrow {lep}#nu;{pttext}::{txc},0.08,0.04".format(lep="e" if channel == "el" else "#mu", 
                                                                                               pttext=ptRangeText,
                                                                                               txc=texCoord)
                    legendCoords = "0.2,0.4,0.75,0.85"
                    drawDataAndMC(hChAsymm1Deta, hChAsymm1Deta_exp, xaxisTitle, "charge asymmetry",
                                  "chargeAsym1D_eta_{fl}_dataAndExp".format(fl=channel),
                                  outname,labelRatioTmp="obs./exp.::0.8,1.2",draw_both0_noLog1_onlyLog2=1,
                                  passCanvas=canvas1D,lumi=options.lumiInt,
                                  lowerPanelHeight=0.35, moreTextLatex=additionalText
                                  )

                    additionalText = additionalTextBackup
                
                
                if options.exptoyfile != "SAME": fexp.Close()

        print ""
