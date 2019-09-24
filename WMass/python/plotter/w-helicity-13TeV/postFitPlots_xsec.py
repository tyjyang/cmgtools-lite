#!/usr/bin/env python
# USAGE: python postFitPlots.py wel_minus_floatPOI.root cards_el -o outputdir [--prefit]

import ROOT, os, re
from array import array
from rollingFunctions import roll1Dto2D, unroll2Dto1D

from make_diff_xsec_cards import getArrayParsingString
from make_diff_xsec_cards import getArrayBinNumberFromValue
from make_diff_xsec_cards import getDiffXsecBinning
from make_diff_xsec_cards import templateBinning
from make_diff_xsec_cards import get_ieta_ipt_from_process_name

ROOT.gStyle.SetOptStat(0)
ROOT.gROOT.SetBatch()
ROOT.TH1.SetDefaultSumw2()

import utilities
utilities = utilities.util()

import sys
sys.path.append(os.getcwd() + "/plotUtils/")
from utility import *

ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetPadRightMargin(0.13)

#h1d_shifted = ROOT.TH1D('shift','',1,0,1)

def normalizeTH2byBinWidth(h2):
    for ix in range(1,1+h2.GetNbinsX()):
        for iy in range(1,1+h2.GetNbinsY()):
            binWidth = h2.GetXaxis().GetBinWidth(ix) * h2.GetYaxis().GetBinWidth(iy)
            h2.SetBinContent(ix,iy, h2.GetBinContent(ix,iy)/binWidth)
            h2.SetBinError(ix,iy, h2.GetBinError(ix,iy)/binWidth)


def normalizeTH1unrolledSingleChargebyBinWidth(h1, h2, unrollAlongX=True):
    # h2 is just used to retrieve the binning, the content is not used
    # unrollAlongX should be true when X is eta and we are taking slices at costant pt
    for ix in range(1,1+h2.GetNbinsX()):
        for iy in range(1,1+h2.GetNbinsY()):
            ibin = 0
            if  unrollAlongX:
                ibin = ix + (iy-1) * h2.GetNbinsX()
            else:
                ibin = iy + (ix-1) * h2.GetNbinsY()
            binWidth = h2.GetXaxis().GetBinWidth(ix) * h2.GetYaxis().GetBinWidth(iy)
            h1.SetBinContent(ibin, h1.GetBinContent(ibin)/binWidth)
            h1.SetBinError(ibin,   h1.GetBinError(ibin)/binWidth)



def dressed2D(h1d,binning,name,title='',shift=0,nCharges=2,nMaskedCha=2, 
              isComb=False, isMuPlot=True, nRecoBins=0):
    if len(binning) == 4:
        n1 = binning[0]; bins1 = array('d', binning[1])
        n2 = binning[2]; bins2 = array('d', binning[3])
        h2_1 = ROOT.TH2F(name, title, n1, bins1, n2, bins2 )
    else:
        n1 = binning[0]; min1 = binning[1]; max1 = binning[2]
        n2 = binning[3]; min2 = binning[4]; max2 = binning[5]
        h2_1 = ROOT.TH2F(name, title, n1, min1, max1, n2, min2, max2)
    h1d_shifted = singleChargeUnrolled(h1d,shift,nCharges,nMaskedCha, 
                                       isComb=isComb, isMuPlot=isMuPlot, nRecoBins=nRecoBins)
    h2_backrolled_1 = roll1Dto2D(h1d_shifted, h2_1 )
    return h2_backrolled_1

# the expproc histograms have the unrolled histograms for both charges (if the charge-combined fit was made) and flavours (if the flavour combined fit was made)
# the order is charge plus, charge minus for charge combination
# if the flavour combinations was made, each charge sector is in turn divided in muons-electrons
# at the end, there are N bins where N is the number of charges in the fit times number of masked channels for each charge

def chargeUnrolledBinShifts(infile,channel,nCharges=2,nMaskedCha=2):
    # guess from a signal with charge defined name
    h1d = infile.Get('expproc_Z_postfit'.format(ch=channel))
    # shift the 1D to remove the empty bins of the other charge
    nbins = int((h1d.GetNbinsX()-nCharges*nMaskedCha)/2)
    ret = {}
    if h1d.Integral(0,nbins)==0:
        ret = {'plus': nbins, 'minus': 0}
    else:
        ret = {'plus': 0, 'minus': nbins}
    return ret

# isComb = isComb, isMuPlot=True if channelCombToPlot == "mu" else False, nRecoBins = nRecoBins
def singleChargeUnrolled(h1d,shift,nCharges=2,nMaskedCha=2, name="shift", isComb=False, isMuPlot=True,  nRecoBins=0):
    extrabins = 0 if 'obs' in h1d.GetName() else nCharges*nMaskedCha
    nbins = int((h1d.GetNbinsX()-extrabins)/2)
    shiftFlavour = 0
    if isComb:        
        if isMuPlot:
            nbins = nRecoBins        
            shiftFlavour = 0
        else:
            shiftFlavour = nbins - nRecoBins # here nRecoBins is the number for electrons, so shiftFlavour is the number of muon bins
            nbins = nRecoBins # only now we update nbins to the actual number of electro bins
    h1d_shifted = ROOT.TH1D(name,'',nbins,0,nbins)
    #h1d_shifted.SetBins(nbins,0.5,float(nbins)+0.5)
    for b in xrange(1,nbins+1):
        h1d_shifted.SetBinContent(b,h1d.GetBinContent(b+shift+shiftFlavour))
        h1d_shifted.SetBinError(b,h1d.GetBinError(b+shift+shiftFlavour))
    return h1d_shifted

def prepareLegend(legWidth=0.50,textSize=0.035,nColumns=3):
    (x1,y1,x2,y2) = (.75-legWidth, .73, .85, .90)
    leg = ROOT.TLegend(x1,y1,x2,y2)
    leg.SetNColumns(nColumns)
    leg.SetFillColor(0)
    leg.SetFillColorAlpha(0,0.6)
    leg.SetShadowColor(0)
    leg.SetLineColor(0)
    leg.SetBorderSize(0)
    leg.SetTextFont(42)
    leg.SetTextSize(textSize)
    return leg


def doRatioHists(data,total,maxRange,fixRange=False,ylabel="Data/pred.",yndiv=505,doWide=False,showStatTotLegend=False,textSize=0.035):
    ratio = data.Clone("data_div"); 
    ratio.Divide(total)

    unity = total.Clone("sim_div");
    
    rmin, rmax =  1,1
    for b in xrange(1,unity.GetNbinsX()+1):
        e,n = unity.GetBinError(b), unity.GetBinContent(b)
        unity.SetBinContent(b, 1 if n > 0 else 0)
        unity.SetBinError(b, e/n if n > 0 else 0)
        rmin = min([ rmin, 1-2*e/n if n > 0 else 1])
        rmax = max([ rmax, 1+2*e/n if n > 0 else 1])
    for b in xrange(1,unity.GetNbinsX()+1):
        if ratio.GetBinContent(b) == 0: continue
        rmin = min([ rmin, ratio.GetBinContent(b) - 2*ratio.GetBinError(b) ]) 
        rmax = max([ rmax, ratio.GetBinContent(b) + 2*ratio.GetBinError(b) ])  
    if rmin < maxRange[0] or fixRange: rmin = maxRange[0]; 
    if rmax > maxRange[1] or fixRange: rmax = maxRange[1];
    if (rmax > 3 and rmax <= 3.4): rmax = 3.4
    if (rmax > 2 and rmax <= 2.4): rmax = 2.4
    unity.SetFillStyle(1001);
    unity.SetFillColor(ROOT.kCyan);
    unity.SetMarkerStyle(1);
    unity.SetMarkerColor(ROOT.kCyan);
    ROOT.gStyle.SetErrorX(0.5);
    unity.Draw("E2");
    unity.GetYaxis().SetRangeUser(rmin,rmax);
    unity.GetXaxis().SetTitleFont(42)
    unity.GetXaxis().SetTitleSize(0.14)
    unity.GetXaxis().SetTitleOffset(0.9)
    unity.GetXaxis().SetLabelFont(42)
    unity.GetXaxis().SetLabelSize(0.1)
    unity.GetXaxis().SetLabelOffset(0.007)
    unity.GetYaxis().SetNdivisions(505)
    unity.GetYaxis().SetTitleFont(42)
    unity.GetYaxis().SetTitleSize(0.14)
    unity.GetYaxis().SetLabelFont(42)
    unity.GetYaxis().SetLabelSize(0.11)
    unity.GetYaxis().SetLabelOffset(0.007)
    unity.GetYaxis().SetDecimals(True) 
    unity.GetYaxis().SetTitle(ylabel)
    total.GetXaxis().SetLabelOffset(999) ## send them away
    total.GetXaxis().SetTitleOffset(999) ## in outer space
    total.GetYaxis().SetTitleSize(0.06)
    unity.GetYaxis().SetTitleOffset(0.20 if doWide else 1.48)
    total.GetYaxis().SetLabelSize(0.05)
    total.GetYaxis().SetLabelOffset(0.007)
    #$ROOT.gStyle.SetErrorX(0.0);
    line = ROOT.TLine(unity.GetXaxis().GetXmin(),1,unity.GetXaxis().GetXmax(),1)
    line.SetLineWidth(2);
    line.SetLineColor(58);
    line.Draw("L")
    ratio.Draw("E SAME" if ratio.ClassName() != "TGraphAsymmErrors" else "PZ SAME");
    leg1 = ROOT.TLegend(0.45, 0.8, 0.7, 0.9)
    leg1.SetFillColor(0)
    leg1.SetShadowColor(0)
    leg1.SetLineColor(0)
    leg1.SetTextFont(42)
    leg1.SetTextSize(0.035*0.7/0.3)
    leg1.AddEntry(unity, "total unc.", "F")
    leg1.Draw()

    return (ratio, unity, line)

def plotPostFitRatio(charge,channel,hratio,outdir,prefix,suffix,passCanvas=None,canvasSize="2400,600", 
                     drawVertLines="", textForLines=[]):

    ROOT.gStyle.SetPadLeftMargin(0.13)
    ROOT.gStyle.SetPadRightMargin(0.02)
    ROOT.gStyle.SetPadBottomMargin(0.3);
    
    plotformat = (int(canvasSize.split(",")[0]), int(canvasSize.split(",")[1]))
    c1 = passCanvas if passCanvas != None else ROOT.TCanvas("canvas","",plotformat[0],plotformat[1])
    c1 = ROOT.TCanvas("c1", "c1", plotformat[0], plotformat[1]); c1.Draw()
    c1.SetWindowSize(plotformat[0] + (plotformat[0] - c1.GetWw()), (plotformat[1] + (plotformat[1] - c1.GetWh())));

    ydiff = hratio.GetBinContent(hratio.GetMaximumBin()) - hratio.GetBinContent(hratio.GetMinimumBin())
    rmin = max(0.1, getMinimumTH(hratio,0.01))
    rmax = min(3., ydiff*0.2 + hratio.GetBinContent(hratio.GetMaximumBin()))
    ROOT.gStyle.SetErrorX(0.5);
    hratio.GetYaxis().SetRangeUser(rmin,rmax);
    hratio.GetXaxis().SetTitleFont(42)
    hratio.GetXaxis().SetTitleSize(0.14)
    hratio.GetXaxis().SetTitleOffset(0.9)
    hratio.GetXaxis().SetLabelFont(42)
    hratio.GetXaxis().SetLabelSize(0.1)
    hratio.GetXaxis().SetLabelOffset(0.007)
    hratio.GetYaxis().SetNdivisions(505)
    hratio.GetYaxis().SetTitleFont(42)
    hratio.GetYaxis().SetTitleSize(0.14)
    hratio.GetYaxis().SetLabelFont(42)
    hratio.GetYaxis().SetLabelSize(0.11)
    hratio.GetYaxis().SetLabelOffset(0.01)
    hratio.GetYaxis().SetDecimals(True) 
    hratio.GetYaxis().SetTitle('post-fit/pre-fit')
    hratio.GetXaxis().SetTitle('unrolled lepton (#eta,p_{T}) bin')
    hratio.GetYaxis().SetTitleOffset(0.40)
    hratio.SetLineColor(ROOT.kBlack)
    hratio.SetStats(0)
    hratio.Draw("HIST" if hratio.ClassName() != "TGraphAsymmErrors" else "PZ SAME");
    line = ROOT.TLine(hratio.GetXaxis().GetXmin(),1,hratio.GetXaxis().GetXmax(),1)
    line.SetLineWidth(2);
    line.SetLineColor(58);
    line.Draw("L")
    lat = ROOT.TLatex()
    lat.SetNDC(); lat.SetTextFont(42)
    lat.DrawLatex(0.15, 0.94, '#bf{CMS} #it{Preliminary}')
    lat.DrawLatex(0.85, 0.94, '35.9 fb^{-1} (13 TeV)')

    # draw vertical lines to facilitate reading of plot
    vertline = ROOT.TLine(36,0,36,c1.GetUymax())
    vertline.SetLineColor(ROOT.kBlack)
    vertline.SetLineStyle(2)
    bintext = ROOT.TLatex()
    #bintext.SetNDC()
    bintext.SetTextSize(0.025)  # 0.03
    bintext.SetTextFont(42)
    if len(textForLines): bintext.SetTextAngle(45 if "#eta" in textForLines[0] else 10)

    if len(drawVertLines):
        #print "drawVertLines = " + drawVertLines
        nptBins = int(drawVertLines.split(',')[0])
        etarange = float(drawVertLines.split(',')[1])        
        offsetXaxisHist = hratio.GetXaxis().GetBinLowEdge(0)
        sliceLabelOffset = 6. if "#eta" in textForLines[0] else 6.
        for i in range(1,nptBins): # do not need line at canvas borders
            vertline.DrawLine(etarange*i-offsetXaxisHist,rmin,etarange*i-offsetXaxisHist,rmax)
        if len(textForLines):
            for i in range(0,len(textForLines)): # we need nptBins texts
                #texoffset = 0.1 * (4 - (i%4))
                #ytext = (1. + texoffset)*ymax/2.  
                ytext = rmax - 0.1*(rmax - rmin)
                bintext.DrawLatex(etarange*i + etarange/sliceLabelOffset, ytext, textForLines[i])


    for ext in ['pdf', 'png']:
        c1.SaveAs('{odir}/{pfx}_{ch}_{flav}_diffXsec_{sfx}.{ext}'.format(odir=outdir,pfx=prefix,ch=charge,flav=channel,sfx=suffix,ext=ext))

def getNormSysts(channel):
    systs = {}
    sysfile = "{cmsswbase}/src/CMGTools/WMass/python/plotter/w-helicity-13TeV/wmass_{ch}/systsEnv.txt".format(cmsswbase=os.environ['CMSSW_BASE'],ch='e' if channel=='el' else 'mu')
    for line in open(sysfile, 'r'):
        if re.match("\s*#.*", line): continue
        line = re.sub("#.*","",line).strip()
        if len(line) == 0: continue
        field = [f.strip() for f in line.split(':')]
        if len(field) == 4 or field[4] == "lnN":
            (name, procmap, binmap, amount) = field[:4]
            if name not in systs: systs[name] = []
            systs[name].append((re.compile(procmap+"$"),amount))
    return systs


if __name__ == "__main__":
    from optparse import OptionParser
    parser = OptionParser(usage="%prog [options] fitresults.root cards_dir")
    parser.add_option('-o','--outdir', dest='outdir', default='.', type='string', help='output directory to save the matrix')
    parser.add_option('-m','--n-mask-chan', dest='nMaskedChannel', default=1, type='int', help='Number of masked channels in the fit for each charge')
    parser.add_option(     '--no2Dplot', dest="no2Dplot", default=False, action='store_true', help="Do not plot templates (but you can still save them in a root file with option -s)");
    parser.add_option(     '--no2Dplot-signal-bin', dest="no2DplotSignalBin", default=True, action='store_true', help="Do not plot templates for each signal bin");
    parser.add_option('-n','--norm-width', dest="normWidth", default=False, action='store_true', help="Normalize histograms by bin area");
    parser.add_option(     '--suffix', dest="suffix", default='', type='string', help="define suffix for each plot");
    parser.add_option('--plot-flavour-from-combination', dest="plotFlavourFromComb", default='', type='string', help="For flavour combination only, specify which flavour to plot: mu or el");
    #parser.add_option('--n-reco-bins', dest='nRecoBins', default=0, type='int', help='Number of reco bins for the flavour being plotted (only for combined fit)')
    #parser.add_option("--pt-bins-name-change", dest="ptBinsNameChange", type="string", default='0,1', help="Comma separated list of pt bins for which the process name should be changed");
    (options, args) = parser.parse_args()
    if len(args) < 2:
        parser.print_usage()
        quit()


    setTDRStyle()

    isComb = False
    channelCombToPlot = ""
    if options.plotFlavourFromComb:
        if options.plotFlavourFromComb not in ["mu","el"]:
            print "Error: option --plot-flavour-from-combination need 'mu' or 'el' as argument. Abort"
            quit()
        isComb = True
        channelCombToPlot = options.plotFlavourFromComb

    nCharges = 2; 
    nMaskedChanPerCharge = options.nMaskedChannel    

    binningfile = "binningPtEta.txt"
    if isComb:
        if channelCombToPlot == "mu": 
            binningfile = "binningPtEta_mu.txt"
        else:
            binningfile = "binningPtEta_el.txt"
    etaPtBinningFile = args[1]+"/" + binningfile
        
    print "binning file: %s " % etaPtBinningFile 
    # get eta-pt binning for both reco and gen
    etaPtBinningVec = getDiffXsecBinning(etaPtBinningFile, "reco")
    recoBins = templateBinning(etaPtBinningVec[0],etaPtBinningVec[1])
    genEtaPtBinningVec = getDiffXsecBinning(etaPtBinningFile, "gen")
    genBins  = templateBinning(genEtaPtBinningVec[0],genEtaPtBinningVec[1])

    nRecoBins = recoBins.NTotBins

    #following array is used to call function dressed2D()
    binning = [recoBins.Neta, recoBins.etaBins, recoBins.Npt, recoBins.ptBins]

    outname = options.outdir
    if isComb: 
        outname = outname + "/" + channelCombToPlot + "/"
    if not os.path.exists(outname):
        os.system("mkdir -p "+outname)
    if os.path.exists("/afs/cern.ch"): os.system("cp /afs/cern.ch/user/e/emanuele/public/index.php "+outname)
    outnamesub = outname+'/postfit_over_prefit/'
    if not os.path.exists(outnamesub):
        os.system("mkdir -p "+outnamesub)
    if os.path.exists("/afs/cern.ch"): os.system("cp /afs/cern.ch/user/e/emanuele/public/index.php "+outnamesub)

    savErrorLevel = ROOT.gErrorIgnoreLevel; ROOT.gErrorIgnoreLevel = ROOT.kError;
    
    infile = ROOT.TFile(args[0], 'read')
    #channel = 'mu' if any(['_mu_postfit' in k.GetName() for k in infile.GetListOfKeys()]) else 'el'
    channel = 'mu' if "diffXsec_mu_" in args[1] else 'el'
    if isComb:
        channel = channelCombToPlot
    print "It seems that you are plotting results for channel ",channel

    full_outfileName = '{odir}/plots_{sfx}.root'.format(odir=outname,sfx=options.suffix)
    outfile = ROOT.TFile(full_outfileName, 'recreate')
    print "Will save 2D templates in file --> " + full_outfileName

    shifts = chargeUnrolledBinShifts(infile,channel,nCharges,nMaskedChanPerCharge)

    colors = {'Wplus_lep'  : ROOT.kRed+2,
              'Wminus_lep' : ROOT.kRed+1,
              'Wplus_el'   : ROOT.kRed+2,
              'Wminus_el'  : ROOT.kRed+1,
              'Wplus_mu'   : ROOT.kRed+2,
              'Wminus_mu'  : ROOT.kRed+1,
              'outliers'   : ROOT.kOrange+2,
              'Top'        : ROOT.kGreen+2,  
              'DiBosons'   : ROOT.kViolet+2, 
              'TauDecaysW' : ROOT.kSpring+9,       
              'Z'          : ROOT.kAzure+2,  
              'Flips'      : ROOT.kCyan,   
              'data_fakes' : ROOT.kGray 
              }
         

    canvas2D = ROOT.TCanvas("canvas2D","",800,700)
    lep = "muon" if channel == "mu" else "electron"
    xaxisname2D = "{l} #eta".format(l=lep)
    yaxisname2D = "{l} p_{{T}} [GeV]".format(l=lep)

    verticalAxisName = "Events / bin [GeV^{-1 }]" if options.normWidth else "Events"
    verticalAxisNameProjX = "Events / bin" if options.normWidth else "Events"  # X is eta
    verticalAxisNameProjY = "Events / bin [GeV^{-1 }]" if options.normWidth else "Events"  # Y is pt


    cwide = ROOT.TCanvas("cwide","",2400,600)                      
    cnarrow = ROOT.TCanvas("cnarrow","",650,700)                      

    canvasRatio = ROOT.TCanvas("canvasRatio","",2400,600)

    # to draw panels in the unrolled plots
    ptBinRanges = []
    for ipt in range(0,recoBins.Npt):
        ptBinRanges.append("p_{{T}} #in [{ptmin:3g}, {ptmax:.3g}]".format(ptmin=recoBins.ptBins[ipt], ptmax=recoBins.ptBins[ipt+1]))

                      
    for charge in ['plus','minus']:
        binshift = shifts[charge]

        print "="*30
        print "charge: " + charge
        print "-"*30
        total_sig = {};     total_sig_unrolled = {}
        bkg_and_data = {};  bkg_and_data_unrolled = {}
        all_procs = {};     all_procs_unrolled = {};    

        single_sig_unrolled = {}; 
        ratios_unrolled = {}

        for prepost in ['postfit','prefit']:
            suffix = prepost+options.suffix
            # doing signal
            canv = ROOT.TCanvas()
            #chfl = '{ch}_{fl}'.format(ch=charge,fl=channel)
            chfl = '{ch}_lep'.format(ch=charge)
            #keyplot = 'W'+chfl+'_'+prepost
            keyplot = prepost
            sigkeyplot = 'W' + chfl + '_' + keyplot

            print "-"*30
            print "doing " + prepost
            print "-"*30
            nGenbins = genBins.Neta * genBins.Npt
            for genEtaPtBin in range(nGenbins): # nGenbins is for outliers (testing whether to use it below together with other backgrounds
                etabin = genEtaPtBin % genBins.Neta
                ptbin = int(genEtaPtBin / genBins.Neta)                
                genkeyplot = 'ieta_{ie}_ipt_{ip}_{k}'.format(ie=etabin,ip=ptbin,k=keyplot)
                genkey = 'ieta_{ie}_ipt_{ip}'.format(ie=etabin,ip=ptbin)

                # outliers managed with other background
                #if genEtaPtBin == nGenbins:
                #    genkeyplot = 'outliers_{k}'.format(k=keyplot)
                #    genkey = 'outliers'

                chs = '+' if charge == 'plus' else '-' 
                # note: for electrons, if some bins were treated as background (e.g. those below pt = 30 GeV), we will also have
                # expproc_W{ch}_el_* histograms, while expproc_W{ch}_lep_* will only contain muons, so we need to check for that
                h1_1 = None
                hname = "expproc_W{chfl}_{gk}".format(chfl=chfl, gk=genkeyplot)
                if channel == "el":
                    hname_testEl = hname.replace("_lep_","_{ch}_".format(ch=channel))
                    h1_1 = infile.Get(hname_testEl)
                    if not h1_1:
                        h1_1 = infile.Get(hname)
                else:
                    h1_1 = infile.Get(hname)
                if not h1_1:
                    print "Error: could not retrieve histogram '{hn}'. Exit".format(hn=hname)
                    quit()
                name2D = 'W{chfl}_{gk}'.format(chfl=chfl, gk=genkeyplot)
                # title2D = ""
                # if "outliers" in name2D:
                #     title2D = 'W{chs}: |#eta| > {etamax}#; p_{{T}} #notin [{ptmin:3g},{ptmax:.3g})'.format(etamax=genBins.etaBins[-1],
                #                                                                                            ptmin=genBins.ptBins[0],
                #                                                                                            ptmax=genBins.ptBins[-1],
                #                                                                                            chs=chs)
                # else:
                #     title2D = 'W{chs}: |#eta| #in [{etamin},{etamax})#; p_{{T}} #in [{ptmin:.3g},{ptmax:.3g})'.format(etamin=genBins.etaBins[etabin],
                #                                                                                                       etamax=genBins.etaBins[etabin+1],
                #                                                                                                       ptmin=genBins.ptBins[ptbin],
                #                                                                                                       ptmax=genBins.ptBins[ptbin+1],
                #                                                                                                       chs=chs)
                title2D = 'W{chs}: |#eta| #in [{etamin},{etamax})#; p_{{T}} #in [{ptmin:.3g},{ptmax:.3g})'.format(etamin=genBins.etaBins[etabin],
                                                                                                                  etamax=genBins.etaBins[etabin+1],
                                                                                                                  ptmin=genBins.ptBins[ptbin],
                                                                                                                  ptmax=genBins.ptBins[ptbin+1],
                                                                                                                  chs=chs)
                h2_backrolled_1 = dressed2D(h1_1,binning,name2D,title2D,binshift,nCharges=2,nMaskedCha=nMaskedChanPerCharge,
                                            isComb=isComb, isMuPlot=True if channelCombToPlot == "mu" else False, nRecoBins=nRecoBins)
                h1_unrolled = singleChargeUnrolled(h1_1,binshift, nMaskedCha=nMaskedChanPerCharge,name="unroll_"+name2D, 
                                                   isComb=isComb, isMuPlot=True if channelCombToPlot == "mu" else False, nRecoBins=nRecoBins)
                if options.normWidth:
                    #normalizeTH2byBinWidth(h2_backrolled_1)
                    normalizeTH1unrolledSingleChargebyBinWidth(h1_unrolled,h2_backrolled_1,True)
                single_sig_unrolled[genkeyplot] = h1_unrolled.Clone('W{chfl}_{gk}_unrolled'.format(chfl=chfl, gk=genkeyplot))
                single_sig_unrolled[genkeyplot].SetDirectory(None)
                if genEtaPtBin == 0:
                    total_sig[sigkeyplot] = h2_backrolled_1.Clone('{kp}'.format(kp=sigkeyplot))
                    total_sig_unrolled[sigkeyplot] = h1_unrolled.Clone('{kp}_unrolled'.format(kp=sigkeyplot))
                    total_sig[sigkeyplot].SetDirectory(None)
                    total_sig_unrolled[sigkeyplot].SetDirectory(None)
                else:
                    if genEtaPtBin < nGenbins:
                        total_sig[sigkeyplot].Add(h2_backrolled_1)
                        total_sig_unrolled[sigkeyplot].Add(h1_unrolled)
                if prepost=='prefit':
                    genbinpostfitkey = '{gkpost}'.format(gkpost=genkeyplot.replace('prefit','postfit')) 
                    if single_sig_unrolled[genbinpostfitkey]!=None:
                        ykeyratio = '{gknopf}_ratio'.format(gknopf=genkeyplot.replace('_prefit',''))
                        ratios_unrolled[ykeyratio] = single_sig_unrolled[genbinpostfitkey].Clone(ykeyratio)
                        ratios_unrolled[ykeyratio].SetDirectory(None)
                        ratios_unrolled[ykeyratio].Divide(single_sig_unrolled[genkeyplot])
                    else:
                        print "Error: something went wrong! Missing key " + genbinpostfitkey
                        quit()
                if not options.no2DplotSignalBin:
                    h2_backrolled_1.Write(name2D)
                    cname = "W{chfl}_{gk}_diffXsec_{sfx}".format(chfl=chfl,gk=genkey,sfx=suffix)
                    drawCorrelationPlot(h2_backrolled_1,xaxisname2D,yaxisname2D,verticalAxisName,cname, "", 
                                        outname, 0,0, False, False, False, 1, palette=57, passCanvas=canvas2D)

            if prepost=='prefit':
                postfitkey  = 'W'+chfl+'_postfit'
                if total_sig_unrolled[postfitkey]!=None:
                    keyratio = 'W'+chfl+'_ratio'
                    ratios_unrolled[keyratio] = total_sig_unrolled[postfitkey].Clone("W{chfl}_unrolled_ratio".format(chfl=chfl))
                    ratios_unrolled[keyratio].SetDirectory(None)
                    ratios_unrolled[keyratio].Divide(total_sig_unrolled[sigkeyplot])
                else:
                    print "Error: something went wrong! Missing key " + postfitkey
                    quit()

            if not options.no2Dplot:
                total_sig[sigkeyplot].Write(total_sig[sigkeyplot].GetName())
                cname = "W{chfl}_TOTALSIG_diffXsec_{sfx}".format(chfl=chfl,sfx=suffix)
                drawCorrelationPlot(total_sig[sigkeyplot],xaxisname2D,yaxisname2D,verticalAxisName,cname, "", 
                                    outname, 0,0, False, False, False, 1, palette=57, passCanvas=canvas2D)
     
     
            # do backgrounds now
            procs=["W{chfl}_outliers".format(chfl=chfl), "Flips","Z","Top","DiBosons","TauDecaysW","data_fakes","obs"]
            titleOut = 'W{chs}: |#eta| > {etamax}; p_{{T}} #notin [{ptmin:3g},{ptmax:.3g})'.format(etamax=genBins.etaBins[-1],
                                                                                                    ptmin=genBins.ptBins[0],
                                                                                                    ptmax=genBins.ptBins[-1],
                                                                                                    chs="+" if charge == "plus" else "-")            
            titles=["out", "charge flips","Drell-Yan","Top","di-bosons","W#rightarrow#tau#nu","QCD","data"]
            procsAndTitles = dict(zip(procs,titles))
            for i,p in enumerate(procs):
                if channel == "mu" and p == "Flips":
                    continue # muons don't have Flips components
                keyplot = p if 'obs' in p else p+'_'+prepost
                if chfl not in keyplot:
                    keyplot = chfl + "_" + keyplot
                hname = "expproc_{p}_{sfx}".format(p=p,sfx=prepost) if 'obs' not in p else "obs"
                h1_1 = infile.Get(hname)
                if not h1_1: 
                    print "Error: could not retrieve histogram '{hn}'. Exit".format(hn=hname)
                    quit()
                pname = p+"_"+prepost if 'obs' not in p else "obs"
                h1_unrolled =  singleChargeUnrolled(h1_1,binshift,nMaskedCha=nMaskedChanPerCharge,name="unroll_"+pname,
                                                    isComb=isComb, isMuPlot=True if channelCombToPlot == "mu" else False, nRecoBins=nRecoBins)
                h2_backrolled_1 = dressed2D(h1_1,binning,pname,titles[i],binshift,nMaskedCha=nMaskedChanPerCharge,
                                            isComb=isComb, isMuPlot=True if channelCombToPlot == "mu" else False, nRecoBins=nRecoBins)
                if options.normWidth:
                    #normalizeTH2byBinWidth(h2_backrolled_1)
                    normalizeTH1unrolledSingleChargebyBinWidth(h1_unrolled,h2_backrolled_1,True)
                bkg_and_data[keyplot] = h2_backrolled_1;  
                bkg_and_data_unrolled[keyplot] = h1_unrolled
                bkg_and_data[keyplot].SetDirectory(None); 
                bkg_and_data_unrolled[keyplot].SetDirectory(None)
                if prepost=='prefit' and p!='obs':
                    postfitkey =  p + '_postfit'
                    if "outliers" not in postfitkey: postfitkey = chfl + "_" + postfitkey
                    if bkg_and_data_unrolled[postfitkey]!=None:
                        keyratio = chfl + "_" + p + '_ratio'
                        ratios_unrolled[keyratio] = bkg_and_data_unrolled[postfitkey].Clone(keyratio)
                        ratios_unrolled[keyratio].SetDirectory(None)
                        ratios_unrolled[keyratio].Divide(bkg_and_data_unrolled[keyplot])
                    else:
                        print "Error: something went wrong! Missing key " + postfitkey
                        quit()

                if not options.no2Dplot:
                    h2_backrolled_1.Write(str(p))
                    cname = "{proc}_{chfl}_diffXsec_{sfx}".format(proc=p, chfl=chfl, sfx=suffix)
                    drawCorrelationPlot(h2_backrolled_1,xaxisname2D,yaxisname2D,verticalAxisName,cname, "", 
                                        outname, 0,0, False, False, False, 1, palette=57, passCanvas=canvas2D)
                    
            # now draw the 1D projections
            all_procs.update(total_sig);   
            all_procs_unrolled.update(total_sig_unrolled);
            all_procs.update(bkg_and_data); 
            all_procs_unrolled.update(bkg_and_data_unrolled); 
            chargesym={'plus':'+', 'minus':'-'}
            procsAndTitles['W{chfl}'.format(chfl=chfl)] = 'W^{{{sign}}}'.format(sign=chargesym[charge])
         
            for p,h in all_procs.iteritems():
                h.integral = h.Integral()
     
            # this has the uncertainty propagated with the full covariance matrix
            h1_expfull = infile.Get('expfull_{sfx}'.format(sfx=prepost))
            expfullName2D = 'expfull_{ch}_{sfx}'.format(ch=charge,sfx=prepost)
            h2_expfull_backrolled = dressed2D(h1_expfull,binning,expfullName2D,expfullName2D,binshift, nMaskedCha=nMaskedChanPerCharge,
                                              isComb=isComb, isMuPlot=True if channelCombToPlot == "mu" else False, nRecoBins=nRecoBins)
            h1_expfull_unrolled = singleChargeUnrolled(h1_expfull, binshift, nMaskedCha=nMaskedChanPerCharge, name="unroll_"+expfullName2D,
                                                       isComb=isComb, isMuPlot=True if channelCombToPlot == "mu" else False, nRecoBins=nRecoBins)
            # can normalize the unrolled, not the 2D
            if options.normWidth:
                #normalizeTH2byBinWidth(h2_expfull_backrolled)
                normalizeTH1unrolledSingleChargebyBinWidth(h1_expfull_unrolled,h2_expfull_backrolled,True)
     
            for projection in ['X','Y']:
                nbinsProj = all_procs[chfl + "_" + 'obs'].GetNbinsY() if projection=='X' else all_procs[chfl + "_" + 'obs'].GetNbinsX()
                hexpfull = h2_expfull_backrolled.ProjectionX("x_expfull_{sfx}_{ch}".format(sfx=prepost,ch=charge),1,nbinsProj,"e") if projection=='X' else h2_expfull_backrolled.ProjectionY("y_expfull_{sfx}_{ch}".format(sfx=prepost,ch=charge),1,nbinsProj,"e")
                hdata = all_procs[chfl + "_" + 'obs'].ProjectionX("x_data_{ch}".format(ch=charge),1,nbinsProj,"e") if projection=='X' else all_procs[chfl + "_" + 'obs'].ProjectionY("y_data_{ch}".format(ch=charge),1,nbinsProj,"e")
                if options.normWidth:
                    hdata.Scale(1., "width")
                    hexpfull.Scale(1., "width")

                htot = hdata.Clone('tot_{ch}'.format(ch=charge)); 
                htot.Reset("ICES"); 
                htot.Sumw2()
                stack = ROOT.THStack("stack_{sfx}_{ch}_proj{pj}".format(sfx=prepost,ch=charge,pj=projection),"")
                
                leg = prepareLegend()
            
                for key,histo in sorted(all_procs.iteritems(), key=lambda (k,v): (v.integral,k)):
                    keycolor = key.replace('_prefit','').replace('_postfit','').replace(chfl + "_", "" )
                    #print "key - keycolor:   %s   %s " % (key, keycolor)
                    if "outliers" in keycolor:
                        keycolor = "outliers"
                    if 'obs' in key or prepost not in key: continue
                    print "%s:  %s   %s " % (projection, key, histo.GetName())
                    proj1d = all_procs[key].ProjectionX(all_procs[key].GetName()+charge+"_px",0,-1,"e") if  projection=='X' else all_procs[key].ProjectionY(all_procs[key].GetName()+charge+"_py",0,-1,"e")
                    if options.normWidth:
                        proj1d.Scale(1., "width")
                    proj1d.SetFillColor(colors[keycolor])
                    stack.Add(proj1d)
                    htot.Add(proj1d) 
                    leg.AddEntry(proj1d,procsAndTitles["{k}".format(k=keycolor if "outliers" not in keycolor else "W{chfl}_outliers".format(chfl=chfl))],'F')
         
                leg.AddEntry(hdata,'data','PE')
                             
                xaxisProj = xaxisname2D if projection == "X" else yaxisname2D
                cnameProj = "projection{pj}_{chfl}_diffXsec_{sfx}".format(pj=projection, chfl=chfl.replace("_lep","_"+channel), sfx=suffix)
                ratioYlabel = "data/pred::" + ("0.9,1.1" if prepost == "prefit" else "0.99,1.01")
                verticalAxisNameProj = verticalAxisNameProjX if projection == "X" else verticalAxisNameProjY
                drawTH1dataMCstack(hdata,stack,xaxisProj, verticalAxisNameProj, cnameProj, outname, leg, ratioYlabel,
                                   1, passCanvas=cnarrow,hErrStack=hexpfull,lumi=35.9, )
            
            # now the unrolled, assign names to run monsterPull.py on them
            hdata_unrolled = singleChargeUnrolled(infile.Get('obs'),binshift,nCharges,nMaskedChanPerCharge, 
                                                  name='unrolled_{ch}_data'.format(ch=charge),
                                                  isComb=isComb, isMuPlot=True if channelCombToPlot == "mu" else False, nRecoBins=nRecoBins).Clone('unrolled_{ch}_data'.format(ch=charge))
            if options.normWidth:
                normalizeTH1unrolledSingleChargebyBinWidth(hdata_unrolled,h2_expfull_backrolled,True)            
            hdata_unrolled.SetDirectory(None)
            htot_unrolled  = hdata_unrolled.Clone('unrolled_{ch}_full'.format(ch=charge)); 
            htot_unrolled.Reset("ICES"); 
            htot_unrolled.Sumw2()
            htot_unrolled.SetDirectory(None)
            stack_unrolled = ROOT.THStack("stack_unrolled_{sfx}_{ch}".format(sfx=prepost,ch=charge),"") 
            leg_unrolled = prepareLegend(legWidth=0.10,nColumns=5)
            for key,histo in sorted(all_procs.iteritems(), key=lambda (k,v): (v.integral,k)):
                keycolor = key.replace('_prefit','').replace('_postfit','').replace(chfl + "_", "" )
                if "outliers" in keycolor: keycolor = "outliers"
                if key=='obs' or prepost not in key: continue
                #print "unrolled:  %s  %s %.3f" % (key, histo.GetName(), histo.Integral())
                print "unrolled {: >35} {: >35}   {y} ".format(key, histo.GetName(), y=str("%.3f" % histo.Integral()) )
                proc_unrolled = all_procs_unrolled[key]
                proc_unrolled.SetFillColor(colors[keycolor])
                stack_unrolled.Add(proc_unrolled)
                htot_unrolled.Add(proc_unrolled)
                #print ">>>  htot_unrolled.Integral() " + str(htot_unrolled.Integral())
                leg_unrolled.AddEntry(proc_unrolled,procsAndTitles[keycolor if "outliers" not in keycolor else "W{chfl}_outliers".format(chfl=chfl)],'F')
            print "Integral data  = " + str(hdata_unrolled.Integral())
            print "Integral stack = " + str(stack_unrolled.GetStack().Last().Integral())
            print "Integral htot  = " + str(htot_unrolled.Integral())
            leg_unrolled.AddEntry(hdata_unrolled,'data','PE')
            leg_unrolled.SetNColumns(5)
            #htot_unrolled.GetYaxis().SetRangeUser(0, 1.8*max(htot_unrolled.GetMaximum(), hdata_unrolled.GetMaximum()))
            #htot_unrolled.GetYaxis().SetTitle(verticalAxisName)
            #htot_unrolled.GetXaxis().SetTitle('unrolled lepton (#eta,p_{T}) bin')

            cnameUnroll = "unrolled_{chfl}_diffXsec_{sfx}".format(pj=projection, chfl=chfl.replace("_lep","_"+channel), sfx=suffix)
            XlabelUnroll = "unrolled template along #eta in %s channel:  #eta #in [%.1f, %.1f]" % (lep, recoBins.etaBins[0], recoBins.etaBins[-1])
            YlabelUnroll = verticalAxisName + "::%.2f,%.2f" % (0, 2.*hdata_unrolled.GetBinContent(hdata_unrolled.GetMaximumBin()))
            ratioYlabel = "data/pred::" + ("0.9,1.1" if prepost == "prefit" else "0.99,1.01")
            #ptBinRanges = []
            #for ipt in range(0,recoBins.Npt):
            #    ptBinRanges.append("p_{{T}} #in [{ptmin:3g},{ptmax:.3g}]".format(ptmin=recoBins.ptBins[ipt], ptmax=recoBins.ptBins[ipt+1]))
            drawTH1dataMCstack(hdata_unrolled,stack_unrolled, XlabelUnroll, YlabelUnroll, cnameUnroll, outname, leg, ratioYlabel,
                               1, passCanvas=cwide,hErrStack=h1_expfull_unrolled,lumi=35.9,wideCanvas=True, leftMargin=0.05,rightMargin=0.02, 
                               drawVertLines="{a},{b}".format(a=recoBins.Npt,b=recoBins.Neta), textForLines=ptBinRanges,etaptbinning=binning)
            hdata_unrolled.Write(); htot_unrolled.Write()

        # plot the postfit/prefit ratio
        print "NOW PLOTTING THE RATIOS..."
        for key,histo in ratios_unrolled.iteritems():
            matchToRemove = charge + "_lep_"
            newkey = key
            if matchToRemove in key: newkey = key.replace(matchToRemove,'')
            if options.no2DplotSignalBin and "ieta_" in newkey: continue
            print "Making unrolled ratio for ",newkey            
            outdir = outnamesub
            plotPostFitRatio(charge,channel,histo,outdir,'postfit2prefit_'+newkey,options.suffix, passCanvas=canvasRatio, 
                             drawVertLines="{a},{b}".format(a=recoBins.Npt,b=recoBins.Neta), textForLines=ptBinRanges)
                
    outfile.Close()

    ROOT.gErrorIgnoreLevel = savErrorLevel;

