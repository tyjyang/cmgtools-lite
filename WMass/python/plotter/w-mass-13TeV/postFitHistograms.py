#!/usr/bin/env python3

# example (only charge plus)
# python w-mass-13TeV/postFitHistograms.py cards/wmass_fixMassWeights_splitW/fit/hessian/fitresults_123456789_Asimov_clipSyst1p3_bbb1_cxs1.root -o plots/testNanoAOD/WmassPlots_jetEta2p4_fixMassWeight_splitW/afterFitPlots/postFitPlots/ --suffix postVFP -l 16.8 -c plus

import os, re
import argparse
from array import array
from rollingFunctions import roll1Dto2D

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
ROOT.TH1.SetDefaultSumw2()
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



def dressed2DfromFit(h1d, binning, name, title='', shift=0,
                     nCharges=2, nMaskedCha=2, 
                     isComb=False, isMuPlot=True, nRecoBins=0, invertXY=True):
    
    if len(binning) == 4:
        n1 = binning[0]; bins1 = array('d', binning[1])
        n2 = binning[2]; bins2 = array('d', binning[3])
        h2_1 = ROOT.TH2D(name, title, n1, bins1, n2, bins2 )
    else:
        n1 = binning[0]; min1 = binning[1]; max1 = binning[2]
        n2 = binning[3]; min2 = binning[4]; max2 = binning[5]
        h2_1 = ROOT.TH2D(name, title, n1, min1, max1, n2, min2, max2)
    h1d_shifted = singleChargeUnrolled(h1d, shift, nCharges, nMaskedCha, 
                                       isComb=isComb, isMuPlot=isMuPlot, nRecoBins=nRecoBins)
    h2_backrolled = roll1Dto2D(h1d_shifted, h2_1, invertXY)
    return h2_backrolled

# the expproc histograms have the unrolled histograms for both charges (if the charge-combined fit was made) and flavours (if the flavour combined fit was made)
# the order is charge plus, charge minus for charge combination
# if the flavour combinations was made, each charge sector is in turn divided in muons-electrons
# at the end, there are N bins where N is the number of charges in the fit times number of masked channels for each charge

def chargeUnrolledBinShifts(infile,channel,nCharges=2,nMaskedCha=2):
    # guess from a signal with charge defined name
    h1d = infile.Get('expproc_Zmumu_postfit')
    if not h1d:
        print("Error in chargeUnrolledBinShifts(): histogram expproc_Zmumu_postfit not found, please check")
        quit()
    # shift the 1D to remove the empty bins of the other charge
    nbins = int((h1d.GetNbinsX()-nCharges*nMaskedCha)/2)
    ret = {}
    if h1d.Integral(0,nbins)==0:
        ret = {'plus': nbins, 'minus': 0}
    else:
        ret = {'plus': 0, 'minus': nbins}
    return ret

# isComb = isComb, isMuPlot=True if channelCombToPlot == "mu" else False, nRecoBins = nRecoBins
def singleChargeUnrolled(h1d, shift, nCharges=2, nMaskedCha=2, name="shift", isComb=False, isMuPlot=True,  nRecoBins=0):
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
    for b in range(1, nbins+1):
        h1d_shifted.SetBinContent(b,h1d.GetBinContent(b+shift+shiftFlavour))
        h1d_shifted.SetBinError(b,h1d.GetBinError(b+shift+shiftFlavour))
    return h1d_shifted

def prepareLegend(legWidth=0.50, textSize=0.035, nColumns=3):
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


def plotPostFitRatio(charge, channel, hratio, outdir, prefix, suffix, passCanvas=None, canvasSize="2400,600", 
                     drawVertLines="", textForLines=[], lumi=36.3, yTitle='postfit/prefit yields'):
    
    plotformat = (int(canvasSize.split(",")[0]), int(canvasSize.split(",")[1]))
    c1 = None
    if passCanvas != None:
        c1 = passCanvas
        c1.SetLeftMargin(0.13)
        c1.SetRightMargin(0.02)
        c1.SetBottomMargin(0.3)
        c1.SetTopMargin(0.1)
    else:
        ROOT.gStyle.SetPadLeftMargin(0.13)
        ROOT.gStyle.SetPadRightMargin(0.02)
        ROOT.gStyle.SetPadBottomMargin(0.3)
        ROOT.gStyle.SetPadTopMargin(0.1)
        c1 = ROOT.TCanvas("c1", "c1", plotformat[0], plotformat[1])
        c1.Draw()
        c1.SetWindowSize(plotformat[0] + (plotformat[0] - c1.GetWw()), (plotformat[1] + (plotformat[1] - c1.GetWh())))

    
    #rmin = max(0.1, 0.9*getMinimumTH(hratio,0.01))
    #ydiff = hratio.GetBinContent(hratio.GetMaximumBin()) - rmin
    #rmax = min(3., ydiff*0.2 + hratio.GetBinContent(hratio.GetMaximumBin()))
    rmin,rmax = getMinMaxHisto(hratio)
    ydiff = rmax - rmin
    if ydiff > 0:
        rmax = min(3., ydiff*0.2 + rmax)
    else:
        rmax = rmax * 1.01
    ROOT.gStyle.SetErrorX(0.5);
    hratio.GetYaxis().SetRangeUser(rmin,rmax)
    hratio.GetXaxis().SetTitleFont(42)
    hratio.GetXaxis().SetTitleSize(0.14)
    hratio.GetXaxis().SetTitleOffset(0.9)
    hratio.GetXaxis().SetLabelFont(42)
    hratio.GetXaxis().SetLabelSize(0.1)
    hratio.GetXaxis().SetLabelOffset(0.007)
    hratio.GetYaxis().SetNdivisions(505)
    hratio.GetYaxis().SetTitleFont(42)
    hratio.GetYaxis().SetTitleSize(0.12)
    hratio.GetYaxis().SetLabelFont(42)
    hratio.GetYaxis().SetLabelSize(0.11)
    hratio.GetYaxis().SetLabelOffset(0.01)
    hratio.GetYaxis().SetDecimals(True) 
    hratio.GetYaxis().SetTitle(yTitle)
    hratio.GetXaxis().SetTitle('unrolled lepton (#eta, p_{T}) bin')
    hratio.GetYaxis().SetTitleOffset(0.45)
    hratio.SetLineColor(ROOT.kBlack)
    hratio.SetStats(0)
    if hratio.ClassName() != "TGraphAsymmErrors":
        hratio.Draw("E2")
        hLine = hratio.Clone(f"{hratio.GetName()}_lineOnly")
        hLine.SetLineColor(ROOT.kBlack)
        hLine.SetFillColor(0)
        hLine.Draw("HIST SAME")
    else:
        hratio.Draw("PZ SAME")
    line = ROOT.TLine(hratio.GetXaxis().GetXmin(),1,hratio.GetXaxis().GetXmax(),1)
    line.SetLineWidth(2);
    line.SetLineColor(58);
    line.Draw("L")
    lat = ROOT.TLatex()
    lat.SetNDC()
    lat.SetTextFont(42)
    lat.SetTextSize(0.05)
    lat.DrawLatex(0.15, 0.94, '#bf{CMS} #it{Preliminary}')
    lat.DrawLatex(0.85, 0.94, '%s fb^{-1} (13 TeV)' % lumi)

    # draw vertical lines to facilitate reading of plot
    vertline = ROOT.TLine(36,0,36,c1.GetUymax())
    vertline.SetLineColor(ROOT.kBlack)
    vertline.SetLineStyle(2)
    bintext = ROOT.TLatex()
    #bintext.SetNDC()
    bintext.SetTextSize(0.03)  # 0.03
    bintext.SetTextFont(42)
    if len(textForLines): bintext.SetTextAngle(0 if "GeV" in textForLines[0] else 45)

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
        c1.SaveAs('{odir}/{pfx}_{sfx}.{ext}'.format(odir=outdir,pfx=prefix,ch=charge,flav=channel,sfx=suffix,ext=ext))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("rootfile", type=str, nargs=1, help="Input file with histograms")
    #parser.add_argument("cardir",   type=str, nargs=1, help="Folder with cards inside")
    parser.add_argument('-o','--outdir', dest='outdir', default='.', type=str, help='output directory to save the matrix')
    parser.add_argument('-m','--n-mask-chan', dest='nMaskedChannel', default=0, type=int, help='Number of masked channels in the fit for each charge (0 if not using masked channels because no signal POIs is used in the fit)')
    parser.add_argument('-c','--charges', dest='charges', choices=['plus', 'minus', 'plus,minus'], default='plus,minus', type=str, help='Charges to process')
    parser.add_argument(     '--no2Dplot', dest="no2Dplot", action='store_true', help="Do not plot templates (but you can still save them in a root file with option -s)");
    parser.add_argument('-n','--norm-width', dest="normWidth", action='store_true', help="Normalize histograms by bin area (mainly if non uniform binning is used)");
    parser.add_argument(     '--suffix', dest="suffix", default='', type=str, help="define suffix for each plot");
    parser.add_argument('-l','--lumi', default=36.3, type=float, help='Integrated luminosity to print in the plot')
    args = parser.parse_args()

    setTDRStyle()

    charges = args.charges.split(",")
    nCharges = len(charges) 
    nMaskedChanPerCharge = args.nMaskedChannel # check if we actually have masked channels, we may not, default should be 0

    # hardcoded eta-pt reco binning for now
    etabins = [round(-2.4 + (0.1 * i), 1) for i in range(0,49)]
    ptbins = [round(26.0 + (1.0 * i), 1) for i in range(0,30)]
    recoBins = templateBinning(etabins, ptbins)
    nRecoBins = recoBins.NTotBins
    #following array is used to call function dressed2DfromFit()
    binning = [recoBins.Neta, recoBins.etaBins, recoBins.Npt, recoBins.ptBins]

    outname = args.outdir
    createPlotDirAndCopyPhp(outname)
    outnamesub = outname+'/postfit_over_prefit/'
    createPlotDirAndCopyPhp(outnamesub)

    savErrorLevel = ROOT.gErrorIgnoreLevel; ROOT.gErrorIgnoreLevel = ROOT.kError;
    
    infile = ROOT.TFile(args.rootfile[0], 'read')
    channel = 'mu' # for now harcoded, we only use muons
    lep = "muon" if channel == "mu" else "electron"
    ############
    # not really needed if only doing muons, but this script is adapted from helicity analysis, where electrons were also used
    isComb = False
    channelCombToPlot = channel
    ############
    xaxisname2D = "{l} #eta".format(l=lep)
    yaxisname2D = "{l} p_{{T}} [GeV]".format(l=lep)
    print(f"It seems that you are plotting results for channel {channel}")

    full_outfileName = '{odir}/plots{sfx}.root'.format(odir=outname,sfx=("_"+args.suffix) if args.suffix else "")
    outfile = ROOT.TFile(full_outfileName, 'recreate')
    print(f"Will save 2D templates in file --> {full_outfileName}")

    shifts = chargeUnrolledBinShifts(infile, channel, nCharges, nMaskedChanPerCharge)

    # FIXME: hardcoded for now, to be checked further
    process_features = {'Wmunu_plus'  : {"color" : ROOT.kRed+2,    "title" : "W^{+}#rightarrow#mu#nu"},
                        'Wmunu_minus' : {"color" : ROOT.kRed+1,    "title" : "W^{-}#rightarrow#mu#nu"},
                        'Top'         : {"color" : ROOT.kGreen+2,  "title" : "top"},
                        'Diboson'     : {"color" : ROOT.kViolet+2, "title" : "dibosons"},
                        'Wtaunu_plus' : {"color" : ROOT.kSpring+9, "title" : "W^{+}#rightarrow#tau#nu"},
                        'Wtaunu_minus': {"color" : ROOT.kSpring+8, "title" : "W^{-}#rightarrow#tau#nu"},
                        'Zmumu'       : {"color" : ROOT.kAzure+2,  "title" : "Z#rightarrow#mu#mu"},
                        'Ztautau'     : {"color" : ROOT.kCyan,     "title" : "Z#rightarrow#tau#tau"},
                        'data_fakes'  : {"color" : ROOT.kGray,     "title" : "QCD multijet"}
    }
         

    canvas2D = ROOT.TCanvas("canvas2D","",800,700)

    verticalAxisName = "Events / bin [GeV^{-1 }]" if args.normWidth else "Events"
    verticalAxisNameProjX = "Events / bin" if args.normWidth else "Events"  # X is eta
    verticalAxisNameProjY = "Events / bin [GeV^{-1 }]" if args.normWidth else "Events"  # Y is pt

    cwide = ROOT.TCanvas("cwide","",2400,600)                      
    cnarrow = ROOT.TCanvas("cnarrow","",650,700)                      

    canvasRatio = ROOT.TCanvas("canvasRatio","",2400,600)

    # to draw panels in the unrolled plots
    ptBinRanges = []
    for ipt in range(0,recoBins.Npt):
        #ptBinRanges.append("p_{{T}} #in [{ptmin:3g}, {ptmax:.3g}]".format(ptmin=recoBins.ptBins[ipt], ptmax=recoBins.ptBins[ipt+1]))
        ptBinRanges.append("#splitline{{[{ptmin},{ptmax}]}}{{GeV}}".format(ptmin=int(recoBins.ptBins[ipt]), ptmax=int(recoBins.ptBins[ipt+1])))

                      
    for charge in charges:
        binshift = shifts[charge]

        print("="*30)
        print(f"charge: {charge}")
        print("-"*30)
        all_procs = {}
        all_procs_unrolled = {}
        ratios_unrolled = {}
        ratios_unc_unrolled = {}
        unc_unrolled = {}
        
        for prepost in ['postfit','prefit']:

            suffix = prepost + args.suffix
            canv = ROOT.TCanvas()
            chfl = charge
            print("-"*30)
            print(f"Doing {prepost}")
            print("-"*30)

            procs = [str(key) for key in process_features.keys()]
            titles = [process_features[p]["title"] for p in procs]
            procs += ["obs"]
            titles += ["Data"]
            procsAndTitles = dict(zip(procs,titles))
            for i,p in enumerate(procs):
                if p == 'obs':
                    pname = p
                    hname = pname
                else:
                    pname = f"{p}_{prepost}"
                    hname = f"expproc_{p}_{prepost}"
                keyplot = f"{chfl}_{pname}"
                h1_1 = infile.Get(hname)
                if not h1_1: 
                    print(f"Error: could not retrieve histogram '{hname}'. Exit")
                    quit()
                # this may not be unrolled as we want, it depends on combinetf
                # in fact, the current version of combinetf unrolls vs pt (it starts from TH2D and concatenate pt shapes for each eta bin, so let's make the 2D first and then unroll as we want, which is versus eta
                # h1_unrolled = singleChargeUnrolled(h1_1, binshift, nMaskedCha=nMaskedChanPerCharge,
                #                                    name=f"unroll_{pname}",
                #                                    isComb=isComb, isMuPlot=True if channelCombToPlot == "mu" else False,
                #                                    nRecoBins=nRecoBins)
                h2_backrolled = dressed2DfromFit(h1_1, binning, pname, titles[i], binshift, nMaskedCha=nMaskedChanPerCharge,
                                                   isComb=isComb, isMuPlot=True if channelCombToPlot == "mu" else False,
                                                   nRecoBins=nRecoBins)
                h1_unrolled = unroll2Dto1D(h2_backrolled, newname=f"unroll_{pname}", cropNegativeBins=False)

                if args.normWidth:
                    #normalizeTH2byBinWidth(h2_backrolled)
                    normalizeTH1unrolledSingleChargebyBinWidth(h1_unrolled, h2_backrolled)
                all_procs[keyplot] = h2_backrolled;  
                all_procs_unrolled[keyplot] = h1_unrolled
                all_procs[keyplot].SetDirectory(0); 
                all_procs_unrolled[keyplot].SetDirectory(0)
                if prepost == 'prefit' and p != 'obs':
                    postfitkey = keyplot.replace("prefit", "postfit")
                    if all_procs_unrolled[postfitkey] != None and all_procs_unrolled[keyplot] != None:
                        keyratio = f"{chfl}_{p}_ratio"
                        # yields
                        ratios_unrolled[keyratio] = copy.deepcopy(all_procs_unrolled[postfitkey].Clone(keyratio))
                        ratios_unrolled[keyratio].SetDirectory(0)
                        ratios_unrolled[keyratio].Divide(all_procs_unrolled[keyplot])
                        ratios_unrolled[keyratio].SetFillColor(process_features[p]["color"])
                        for ibin in range(1,ratios_unrolled[keyratio].GetNbinsX()+1):
                            unc = 0.0 if all_procs_unrolled[keyplot].GetBinContent(ibin) == 0.0 else (all_procs_unrolled[postfitkey].GetBinError(ibin) / all_procs_unrolled[keyplot].GetBinContent(ibin))
                            ratios_unrolled[keyratio].SetBinError(ibin, unc)
                        # uncertainties
                        ratios_unc_unrolled[keyratio] = copy.deepcopy(all_procs_unrolled[postfitkey].Clone(keyratio+"_unc"))
                        ratios_unc_unrolled[keyratio].SetDirectory(0)
                        ratios_unc_unrolled[keyratio].SetFillColor(process_features[p]["color"])
                        unc_unrolled[postfitkey] = copy.deepcopy(all_procs_unrolled[postfitkey].Clone(postfitkey+"_unc"))
                        unc_unrolled[keyplot] = copy.deepcopy(all_procs_unrolled[keyplot].Clone(keyplot+"_unc"))
                        unc_unrolled[postfitkey].Reset("ICESM")
                        unc_unrolled[keyplot].Reset("ICESM")
                        for ibin in range(1,ratios_unc_unrolled[keyratio].GetNbinsX()+1):
                            val = 0.0 if all_procs_unrolled[keyplot].GetBinError(ibin) == 0.0 else (all_procs_unrolled[postfitkey].GetBinError(ibin) / all_procs_unrolled[keyplot].GetBinError(ibin))
                            ratios_unc_unrolled[keyratio].SetBinContent(ibin, val)
                            ratios_unc_unrolled[keyratio].SetBinError(ibin, 0.0)
                            unc_unrolled[postfitkey].SetBinContent(ibin, all_procs_unrolled[postfitkey].GetBinError(ibin))
                            unc_unrolled[keyplot].SetBinContent(ibin, all_procs_unrolled[keyplot].GetBinError(ibin))
                        # try plotting both
                        hists = [copy.deepcopy(all_procs_unrolled[postfitkey].Clone(f"postfit_{chfl}_{p}")),
                                 copy.deepcopy(all_procs_unrolled[keyplot].Clone(f"prefit_{chfl}_{p}"))]
                        vertLines="{a},{b}".format(a=recoBins.Npt,b=recoBins.Neta)
                        legs = ["postfit","prefit"]
                        # yields
                        drawNTH1(hists, legs, "unrolled lepton (#eta, p_{T}) bin", "Events", f"postfitAndprefit_yields_chan{chfl}_{p}", outnamesub, leftMargin=0.06, rightMargin=0.02, labelRatioTmp="postfit/prefit::0.8,1.2", legendCoords="0.45,0.8,0.92,1.0;2", passCanvas=cwide, drawLumiLatex=True, lumi=args.lumi, drawVertLines=vertLines, textForLines=ptBinRanges, yAxisExtendConstant=1.25, markerStyleFirstHistogram=1, fillStyleSecondHistogram=1001, colorVec=[ROOT.kGray], moreTextLatex=f"{process_features[p]['title']}::0.3,0.95,0.08,0.055")
                        # uncertainties
                        hists = [unc_unrolled[postfitkey], unc_unrolled[keyplot]]
                        drawNTH1(hists, legs, "unrolled lepton (#eta, p_{T}) bin", "Uncertainty", f"postfitAndprefit_uncertainty_chan{chfl}_{p}", outnamesub, leftMargin=0.06, rightMargin=0.02, labelRatioTmp="postfit/prefit", legendCoords="0.45,0.8,0.92,1.0;2", passCanvas=cwide, drawLumiLatex=True, lumi=args.lumi, drawVertLines=vertLines, textForLines=ptBinRanges, yAxisExtendConstant=1.25, markerStyleFirstHistogram=1, fillStyleSecondHistogram=1001, colorVec=[ROOT.kGray], moreTextLatex=f"{process_features[p]['title']}::0.3,0.95,0.08,0.055", setRatioRangeFromHisto=True, setOnlyLineRatio=True)

                        
                    else:
                        print("Error: something went wrong! Missing either {postfitkey} or {keyplot}")
                        quit()

                if not args.no2Dplot:
                    cname = f"{p}_{chfl}_{suffix}"
                    h2_backrolled.Write(cname)
                    drawCorrelationPlot(h2_backrolled, xaxisname2D, yaxisname2D, verticalAxisName, cname, "", 
                                        outname, 0,0, False, False, False, 1, palette=57, passCanvas=canvas2D)
                    
            # now draw the 1D projections         
            sortedKeys = list(sorted(all_procs.keys(), key = lambda k: all_procs[k].Integral()))
                
            # this has the uncertainty propagated with the full covariance matrix
            h1_expfull = infile.Get(f"expfull_{prepost}")
            expfullName2D = f"expfull_{charge}_{prepost}"
            h2_expfull_backrolled = dressed2DfromFit(h1_expfull, binning, expfullName2D, expfullName2D,
                                                     binshift, nMaskedCha=nMaskedChanPerCharge, isComb=isComb,
                                                     isMuPlot=True if channelCombToPlot == "mu" else False, nRecoBins=nRecoBins)
            h1_expfull_unrolled = unroll2Dto1D(h2_expfull_backrolled, newname=f"unroll_{expfullName2D}", cropNegativeBins=False)

            h2_expfull_backrolled.Write(f"expfull_{prepost}_{charge}")
            # can normalize the unrolled, not the 2D
            if args.normWidth:
                #normalizeTH2byBinWidth(h2_expfull_backrolled)
                normalizeTH1unrolledSingleChargebyBinWidth(h1_expfull_unrolled, h2_expfull_backrolled)
     
            for projection in ['X', 'Y']:
                if projection=='X':
                    nbinsProj = all_procs[f"{chfl}_obs"].GetNbinsY()
                    projName = f"expfull_{prepost}_{charge}_px"
                    projNameData = f"data_{prepost}_{charge}_px"
                    hexpfull = h2_expfull_backrolled.ProjectionX(projName,     1, nbinsProj, "e")
                    hdata = all_procs[f"{chfl}_obs"].ProjectionX(projNameData, 1, nbinsProj, "e")
                else:
                    nbinsProj = all_procs[f"{chfl}_obs"].GetNbinsX()
                    projName = f"expfull_{prepost}_{charge}_py"
                    projNameData = f"data_{prepost}_{charge}_py"
                    hexpfull = h2_expfull_backrolled.ProjectionY(projName, 1, nbinsProj, "e")
                    hdata = all_procs[f"{chfl}_obs"].ProjectionY(projNameData, 1, nbinsProj, "e")
                if args.normWidth:
                    hdata.Scale(1., "width")
                    hexpfull.Scale(1., "width")

                hdata.Write()
                hexpfull.Write()
                htot = hdata.Clone(f"tot_{charge}") 
                htot.Reset("ICES");
                htot.Sumw2()
                stack = ROOT.THStack(f"stack_{prepost}_{charge}_proj{projection}", "")
                
                leg = prepareLegend()

                listOfProj = []
                for key in sortedKeys:
                    histo = all_procs[key]
                    keycolor = key.replace('_prefit','').replace('_postfit','')
                    if keycolor.startswith(f"{charge}_"):
                        keycolor = "_".join(keycolor.split("_")[1:])
                    if 'obs' in key or prepost not in key: continue
                    print(f"{projection}:  {key}   {histo.GetName()}   {keycolor}")
                    if  projection == 'X':
                        proj1d = all_procs[key].ProjectionX(all_procs[key].GetName() + charge+"_px", 0, -1, "e")
                    else:
                        proj1d = all_procs[key].ProjectionY(all_procs[key].GetName() + charge+"_py", 0, -1, "e")
                    if args.normWidth:
                        proj1d.Scale(1., "width")
                    proj1d.SetFillColor(process_features[keycolor]["color"])
                    stack.Add(proj1d)
                    proj1d.Write()
                    htot.Add(proj1d) 
                    listOfProj.append([proj1d, procsAndTitles[keycolor]])
                    
                leg.AddEntry(hdata,'data','PE')
                for pair in reversed(listOfProj):
                    leg.AddEntry(pair[0], pair[1], 'F')
                                      
                xaxisProj = xaxisname2D if projection == "X" else yaxisname2D
                cnameProj = f"projection{projection}_{chfl}_{suffix}"
                ratioYlabel = "data/pred::" + ("0.9,1.1" if prepost == "prefit" else "0.99,1.01")
                verticalAxisNameProj = verticalAxisNameProjX if projection == "X" else verticalAxisNameProjY
                drawTH1dataMCstack(hdata, stack, xaxisProj, verticalAxisNameProj, cnameProj, outname, leg, ratioYlabel,
                                   1, passCanvas=cnarrow, hErrStack=hexpfull, lumi=args.lumi)
            
            # hdata_unrolled = singleChargeUnrolled(infile.Get('obs'), binshift, nCharges, nMaskedChanPerCharge, 
            #                                       name=f"unrolled_{charge}_data",
            #                                       isComb=isComb, isMuPlot=True if channelCombToPlot == "mu" else False,
            #                                       nRecoBins=nRecoBins).Clone('unrolled_{ch}_data'.format(ch=charge))
            dataName2D = f"data_{charge}"
            h2_data_backrolled = dressed2DfromFit(infile.Get('obs'), binning, dataName2D, dataName2D,
                                                     binshift, nMaskedCha=nMaskedChanPerCharge, isComb=isComb,
                                                     isMuPlot=True if channelCombToPlot == "mu" else False, nRecoBins=nRecoBins)
            hdata_unrolled = unroll2Dto1D(h2_data_backrolled, newname=f"unrolled_{charge}_data", cropNegativeBins=False)


            if args.normWidth:
                normalizeTH1unrolledSingleChargebyBinWidth(hdata_unrolled, h2_data_backrolled)
            hdata_unrolled.SetDirectory(0)
            htot_unrolled  = hdata_unrolled.Clone(f"unrolled_{charge}_full")
            htot_unrolled.Reset("ICES")
            htot_unrolled.Sumw2()
            htot_unrolled.SetDirectory(0)
            stack_unrolled = ROOT.THStack(f"stack_unrolled_{prepost}_{charge}", "")
            leg_unrolled = prepareLegend(legWidth=0.60, textSize=0.045, nColumns=5)
            listOfProj = []
            for key in sortedKeys:
                histo = all_procs[key]
                keycolor = key.replace('_prefit','').replace('_postfit','')
                if keycolor.startswith(f"{charge}_"):
                    keycolor = "_".join(keycolor.split("_")[1:])
                if key=='obs' or prepost not in key: continue
                print("unrolled {: >35} {: >35}   {y} ".format(key, histo.GetName(), y=str("%.3f" % histo.Integral())) )
                proc_unrolled = all_procs_unrolled[key]
                proc_unrolled.SetFillColor(process_features[keycolor]["color"])
                stack_unrolled.Add(proc_unrolled)
                htot_unrolled.Add(proc_unrolled)
                listOfProj.append([proc_unrolled, procsAndTitles[keycolor]])
                
            leg_unrolled.AddEntry(hdata_unrolled, 'data', 'PE')
            for pair in reversed(listOfProj):
                leg_unrolled.AddEntry(pair[0], pair[1], 'F')

            print("Integral data  = " + str(hdata_unrolled.Integral()))
            print("Integral stack = " + str(stack_unrolled.GetStack().Last().Integral()))
            print("Integral htot  = " + str(htot_unrolled.Integral()))
            #htot_unrolled.GetYaxis().SetRangeUser(0, 1.8*max(htot_unrolled.GetMaximum(), hdata_unrolled.GetMaximum()))
            #htot_unrolled.GetYaxis().SetTitle(verticalAxisName)
            #htot_unrolled.GetXaxis().SetTitle('unrolled lepton (#eta,p_{T}) bin')

            cnameUnroll = f"unrolled_{chfl}_{suffix}"
            XlabelUnroll = "unrolled template along #eta in %s channel:  #eta #in [%.1f, %.1f]" % (lep, recoBins.etaBins[0], recoBins.etaBins[-1])
            YlabelUnroll = verticalAxisName + "::%.2f,%.2f" % (0, 2.*hdata_unrolled.GetBinContent(hdata_unrolled.GetMaximumBin()))
            ratioYlabel = "data/pred::" + ("0.9,1.1" if prepost == "prefit" else "0.99,1.01")
            drawTH1dataMCstack(hdata_unrolled, stack_unrolled, XlabelUnroll, YlabelUnroll, cnameUnroll, outname,
                               leg_unrolled, ratioYlabel, 1, passCanvas=cwide, hErrStack=h1_expfull_unrolled, lumi=args.lumi,
                               wideCanvas=True, leftMargin=0.05,rightMargin=0.02, 
                               drawVertLines="{a},{b}".format(a=recoBins.Npt,b=recoBins.Neta),
                               textForLines=ptBinRanges, etaptbinning=binning, textSize=0.04, textAngle=0)
            #hdata_unrolled.Write()
            #htot_unrolled.Write()

        # plot the postfit/prefit ratio
        # these are probably no longer needed, we make similar distributions above, with the actual distributions and ratios
        # print("NOW PLOTTING THE RATIOS...")
        # for key,histo in ratios_unrolled.items():
        #     print(f"Making unrolled ratio of yields for {key}")
        #     plotPostFitRatio(charge, channel, histo, outnamesub, f"postfit2prefit_yields_chan{key}", args.suffix, 
        #                      drawVertLines="{a},{b}".format(a=recoBins.Npt,b=recoBins.Neta), textForLines=ptBinRanges,
        #                      lumi=args.lumi)
        # for key,histo in ratios_unc_unrolled.items():
        #     print(f"Making unrolled ratio of uncertainties for {key}")
        #     plotPostFitRatio(charge, channel, histo, outnamesub, f"postfit2prefit_uncertainty_chan{key}", args.suffix, 
        #                      drawVertLines="{a},{b}".format(a=recoBins.Npt,b=recoBins.Neta), textForLines=ptBinRanges,
        #                      lumi=args.lumi, yTitle="postfit/prefit #sigma")
                
    outfile.Close()

    ROOT.gErrorIgnoreLevel = savErrorLevel;

