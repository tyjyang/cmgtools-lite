#!/usr/bin/env python
# USAGE: python postFitPlots.py wel_minus_floatPOI.root cards_el -o outputdir [--prefit]
import ROOT, os, re, copy
from array import array
from rollingFunctions import roll1Dto2D, unroll2Dto1D, reverseUnroll2Dto1D

from make_diff_xsec_cards import getDiffXsecBinning
from make_diff_xsec_cards import templateBinning

import utilities
utilities = utilities.util()
print 'right pad margin', ROOT.gStyle.GetPadRightMargin()

def getbinning(splitline):
    bins = splitline[5]
    if '*' in bins:
        etabins = list( float(i) for i in bins.replace(' ','').split('*')[0].replace('\'[','').replace(']','').split(',') )
        ptbins  = list( float(i) for i in bins.replace(' ','').split('*')[1].replace('[','').replace(']\'','').split(',') )
        nbinseta = len(etabins)-1
        nbinspt  = len( ptbins)-1
        binning = [nbinseta, etabins, nbinspt, ptbins]
    else:
        bins = bins.replace('\'','')
        nbinseta = int( bins.split(',')[0] )
        nbinspt  = int( bins.split(',')[3] )
        etamin   = float( bins.split(',')[1] )
        etamax   = float( bins.split(',')[2] )
        ptmin    = float( bins.split(',')[4] )
        ptmax    = float( bins.split(',')[5] )
    return binning

ROOT.gStyle.SetOptStat(0)

marL = ROOT.gStyle.GetPadLeftMargin()
marR = ROOT.gStyle.GetPadRightMargin()
marT = ROOT.gStyle.GetPadTopMargin()
marB = ROOT.gStyle.GetPadBottomMargin()


def dressed2D(h1d,binning,name,title='',shift=0,nCharges=2,nMaskedCha=2):
    if len(binning) == 4:
        n1 = binning[0]; bins1 = array('d', binning[1])
        n2 = binning[2]; bins2 = array('d', binning[3])
        h2_1 = ROOT.TH2F(name, title, n1, bins1, n2, bins2 )
    else:
        n1 = binning[0]; min1 = binning[1]; max1 = binning[2]
        n2 = binning[3]; min2 = binning[4]; max2 = binning[5]
        h2_1 = ROOT.TH2F(name, title, n1, min1, max1, n2, min2, max2)
    #h1d_shifted = singleChargeUnrolled(h1d,shift,nCharges,nMaskedCha)
    h2_backrolled_1 = roll1Dto2D(h1d, h2_1 )
    h2_backrolled_1 .GetXaxis().SetTitle('lepton #eta')
    h2_backrolled_1 .GetYaxis().SetTitle('lepton p_{T} (GeV)')
    h2_backrolled_1 .GetZaxis().SetRangeUser(0.01*h2_backrolled_1.GetMaximum(),1.1*h2_backrolled_1.GetMaximum())
    return h2_backrolled_1

def chargeUnrolledBinShifts(infile,channel,nCharges=2,nMaskedCha=2):
    # guess from a signal with charge defined name
    h1d = infile.Get('expproc_Wplus_left_Ybin_0_postfit')
    # shift the 1D to remove the empty bins of the other charge
    nbins = int((h1d.GetNbinsX()-nCharges*nMaskedCha)/2)
    ret = {}
    if h1d.Integral(0,nbins)==0:
        ret = {'plus': nbins, 'minus': 0}
    else:
        ret = {'plus': 0, 'minus': nbins}
    return ret

def singleChargeUnrolled(h1d,shift,nCharges=2,nMaskedCha=2):
    extrabins = 0 if 'obs' in h1d.GetName() else nCharges*nMaskedCha
    nbins = int((h1d.GetNbinsX()-extrabins)/2)
    h1d_shifted = ROOT.TH1D('shift','',nbins,0,nbins) # from 0 to nbins, not from 1
    for b in xrange(nbins):
        h1d_shifted.SetBinContent(b+1,h1d.GetBinContent(b+shift+1))
        h1d_shifted.SetBinError(b+1,h1d.GetBinError(b+shift+1))
    return h1d_shifted

def prepareLegend(legWidth=0.50,textSize=0.035,xmin=0.75):
    (x1,y1,x2,y2) = (xmin-legWidth, .70, .85, .87)
    leg = ROOT.TLegend(x1,y1,x2,y2)
    leg.SetNColumns(3)
    leg.SetFillColor(0)
    leg.SetFillColorAlpha(0,0.6)
    leg.SetShadowColor(0)
    leg.SetLineColor(0)
    leg.SetBorderSize(0)
    leg.SetTextFont(42)
    leg.SetTextSize(textSize)
    return leg

def plotOne(charge,channel,stack,htot,hdata,legend,outdir,prefix,suffix,veryWide=False, 
            drawVertLines="", # "12,36": format --> N of sections (e.g: 12 pt bins), and N of bins in each section (e.g. 36 eta bins), assuming uniform bin width
            textForLines=[],  # text in each panel                       
            outtext=''
            ):

    ## Prepare split screen
    plotformat = (2400,600) if veryWide else (600,750)
    c1 = ROOT.TCanvas("c1", "c1", plotformat[0], plotformat[1]); c1.Draw()
    if not 'projection' in prefix:
        ROOT.gStyle.SetPadLeftMargin(0.07 if veryWide else 0.18)
        ROOT.gStyle.SetPadRightMargin(0.02 if veryWide else 0.13)
        ROOT.gStyle.SetPadTopMargin(marT)
        ROOT.gStyle.SetPadBottomMargin(marB)
        c1.SetWindowSize(plotformat[0] + (plotformat[0] - c1.GetWw()), (plotformat[1] + (plotformat[1] - c1.GetWh())));
        p1 = ROOT.TPad("pad1","pad1",0,0.29,1,0.99);
        p1.SetBottomMargin(0.03);
        p1.Draw();
        p2 = ROOT.TPad("pad2","pad2",0,0,1,0.31);
        p2.SetTopMargin(0);
        p2.SetBottomMargin(0.3);
        p2.SetFillStyle(0);
        p2.Draw();
        p1.cd();
    else:
        c1.GetPad(0).SetTopMargin   (0.09)
        c1.GetPad(0).SetBottomMargin(0.15)
        c1.GetPad(0).SetLeftMargin  (0.17)
        c1.GetPad(0).SetRightMargin (0.04)
        c1.GetPad(0).SetTickx(1)
        c1.GetPad(0).SetTicky(1)
        ##ROOT.gStyle.SetPadLeftMargin  (0.17)
        ##ROOT.gStyle.SetPadRightMargin (0.05)
        ##ROOT.gStyle.SetPadTopMargin   (0.12)
        ##ROOT.gStyle.SetPadBottomMargin(0.33)
        c1.cd()
    ## Draw absolute prediction in top frame
    offset = 0.45 if veryWide else 0.62
    htot.GetYaxis().SetTitleOffset(offset)
    htot.SetTitle('')
    htot.GetYaxis().SetTitleOffset(1.35)
    htot.GetYaxis().SetTitleSize(0.05)
    htot.GetXaxis().SetLabelSize(0.04)
    htot.GetYaxis().SetLabelSize(0.04)
    htot.GetXaxis().SetTitleSize(0.05)
    htot.Draw("HIST")
    #htot.SetLabelOffset(9999.0);
    #htot.SetTitleOffset(9999.0);
    stack.Draw("HIST F SAME")
    hdata.SetMarkerColor(ROOT.kBlack)
    hdata.SetMarkerStyle(ROOT.kFullCircle)
    hdata.Draw("E SAME")
    htot.Draw("AXIS SAME")
    totalError = utilities.doShadedUncertainty(htot)            
    legend.Draw()
    lat = ROOT.TLatex()
    lat.SetNDC(); lat.SetTextFont(42)
    if veryWide:
        x1 = 0.09; x2 = 0.85
    else:
        x1 = 0.26; x2 = 0.60
    lat.DrawLatex(x1, 0.92, '#bf{CMS}')# #it{Preliminary}')
    lat.DrawLatex(x2, 0.92, '35.9 fb^{-1} (13 TeV)')
    

    # set Y range a little above the current value
    ymaxBackup = htot.GetMaximum()
    vertline = ROOT.TLine(36,0,36,c1.GetUymax())
    vertline.SetLineColor(ROOT.kBlack)
    vertline.SetLineStyle(2)
    bintext = ROOT.TLatex()
    #bintext.SetNDC()
    bintext.SetTextSize(0.04)
    bintext.SetTextFont(42)

    if len(drawVertLines):
        #print "drawVertLines = " + drawVertLines
        nptBins = int(drawVertLines.split(',')[0])
        etarange = float(drawVertLines.split(',')[1])        
        for i in range(1,nptBins): # do not need line at canvas borders
            vertline.DrawLine(etarange*i,0,etarange*i,ymaxBackup)
        if len(textForLines):
            for i in range(0,len(textForLines)): # we need nptBins texts
                bintext.DrawLatex(etarange*i + etarange/6., 0.7 *ymaxBackup, textForLines[i])


    if not 'projection' in prefix:
        p2.cd()
        maxrange = [0.97,1.03] if prepost == 'postfit' else [0.90,1.10]
        if 'projection' in prefix and prepost == 'postfit':
            maxrange = [0.99,1.01]
        rdata,rnorm,rline = doRatioHists(hdata, htot, maxRange=maxrange, fixRange=True, doWide=veryWide)
    c1.cd()
    outf_basename = '{odir}/{pfx}_{ch}_{flav}_PFMT40_absY_{sfx}.'.format(odir=outdir,pfx=prefix,ch=charge,flav=channel,sfx=suffix)
    for ext in ['pdf', 'png']:
        c1.SaveAs(outf_basename+ext)
    if outtext:
        ofo = open(outf_basename+'txt', 'w')
        ofo.write(outtext)
        ofo.close()
    ROOT.gStyle.SetPadLeftMargin  (marL)
    ROOT.gStyle.SetPadRightMargin (marR)
    ROOT.gStyle.SetPadTopMargin   (marT)
    ROOT.gStyle.SetPadBottomMargin(marB)
        


def doRatioHists(data,total,maxRange,fixRange=False,ylabel="Data/pred.",yndiv=505,doWide=False,showStatTotLegend=False,textSize=0.035):
    ratio = data.Clone("data_div"); 
    #ratio.Divide(total)  # this combines the uncertainty on data with the one on total. 
    # Since the uncertainty on total is already in the band centered around 1, it is better to divide the data by a clone of total with 0 uncertainty
    total_noErr = total.Clone("total_noErr")
    for i in range (0,2+total_noErr.GetNbinsX()):
        total_noErr.SetBinError(i, 0)
    ratio.Divide(total_noErr)

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

def plotPostFitRatio(charge,channel,hratio,outdir,prefix,suffix, drawVertLines="", textForLines=[]):
    ROOT.gStyle.SetPadLeftMargin(0.13)
    ROOT.gStyle.SetPadRightMargin(0.02)
    ROOT.gStyle.SetPadBottomMargin(0.3);

    plotformat = (2400,600)
    c1 = ROOT.TCanvas("c1", "c1", plotformat[0], plotformat[1]); c1.Draw()
    c1.SetWindowSize(plotformat[0] + (plotformat[0] - c1.GetWw()), (plotformat[1] + (plotformat[1] - c1.GetWh())));

    ydiff = hratio.GetBinContent(hratio.GetMaximumBin()) - hratio.GetBinContent(hratio.GetMinimumBin())
    rmin = max(0.1, hratio.GetBinContent(hratio.GetMinimumBin())); rmax = min(10., ydiff*0.2 + hratio.GetBinContent(hratio.GetMaximumBin()))
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
    hratio.Draw("HIST" if hratio.ClassName() != "TGraphAsymmErrors" else "PZ SAME");
    line = ROOT.TLine(hratio.GetXaxis().GetXmin(),1,hratio.GetXaxis().GetXmax(),1)
    line.SetLineWidth(2);
    line.SetLineColor(58);
    line.Draw("L")
    lat = ROOT.TLatex()
    lat.SetNDC(); lat.SetTextFont(42)
    lat.DrawLatex(0.15, 0.94, '#bf{CMS}')# #it{Preliminary}')
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
        c1.SaveAs('{odir}/{pfx}_{ch}_{flav}_PFMT40_absY_{sfx}.{ext}'.format(odir=outdir,pfx=prefix,ch=charge,flav=channel,sfx=suffix,ext=ext))

    ROOT.gStyle.SetPadLeftMargin(marL)
    ROOT.gStyle.SetPadRightMargin(marR)
    ROOT.gStyle.SetPadBottomMargin(marB);

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

ROOT.gROOT.SetBatch()

if __name__ == "__main__":
    from optparse import OptionParser
    parser = OptionParser(usage="%prog [options] fitresults.root cards_dir")
    parser.add_option('-o','--outdir', dest='outdir', default='.', type='string', help='output directory to save the matrix')
    parser.add_option(     '--no2Dplot', dest="no2Dplot", default=False, action='store_true', help="Do not plot templates (but you can still save them in a root file with option -s)");
    parser.add_option('-m','--n-mask-chan', dest='nMaskedChannel', default=1, type='int', help='Number of masked channels in the fit for each charge')
    parser.add_option(     '--suffix', dest="suffix", default='', type='string', help="define suffix for each plot");
    parser.add_option(     '--reverseUnrolling', dest="reverseUnrolling", default=False, action='store_true', help="do the unrolling in the ohter way around (pt in eta chunks)");
    (options, args) = parser.parse_args()

    groupJobs=5 # used in make_helicity_cards.py
    nCharges = 2; nMaskedChanPerCharge = options.nMaskedChannel; 
    
    if len(args) < 1:
        parser.print_usage()
        quit()

    ## get the binning of the YW bins
    ybinfile = args[1]+'/binningYW.txt'
    ybinfile = open(ybinfile, 'r')
    ybins = eval(ybinfile.read())
    ybinfile.close()

    outname = options.outdir
    if not os.path.exists(outname):
        os.system("mkdir -p "+outname)
        if os.path.exists("/afs/cern.ch"): os.system("cp /afs/cern.ch/user/g/gpetrucc/php/index.php "+outname)
    outnamesub = outname+'/singleRapidities'
    if not os.path.exists(outnamesub):
        os.system("mkdir -p "+outnamesub)
        if os.path.exists("/afs/cern.ch"): os.system("cp /afs/cern.ch/user/g/gpetrucc/php/index.php "+outnamesub)

    savErrorLevel = ROOT.gErrorIgnoreLevel; ROOT.gErrorIgnoreLevel = ROOT.kError;
    
    infile = ROOT.TFile(args[0], 'read')
    channel = utilities.getChannelFromFitresults(infile)
    print "From the list of histograms it seems that you are plotting results for channel ",channel

    full_outfileName = '{odir}/plots_{sfx}.root'.format(odir=outname,sfx=options.suffix)
    outfile = ROOT.TFile(full_outfileName, 'recreate')
    print "Will save 2D templates in file --> " + full_outfileName

    ## get the pT and eta binning from the file in the directory
    etaPtBinningFile = args[1]+"/binningPtEta.txt"
    # get eta-pt binning for reco
    etaPtBinningVec = getDiffXsecBinning(etaPtBinningFile, "reco")
    recoBins = templateBinning(etaPtBinningVec[0],etaPtBinningVec[1])

    binning = [recoBins.Neta, recoBins.etaBins, recoBins.Npt, recoBins.ptBins]

    # to draw panels in the unrolled plots
    ptBinRanges = []
    for ipt in range(0,recoBins.Npt):
        ptBinRanges.append("p_{{T}} #in [{ptmin:3g}, {ptmax:.3g}]".format(ptmin=recoBins.ptBins[ipt], ptmax=recoBins.ptBins[ipt+1]))


    shifts = chargeUnrolledBinShifts(infile,channel,nCharges,nMaskedChanPerCharge)
    print 'THESE ARE THE BINSHIFTS', shifts

    reverse = 'reverse' if options.reverseUnrolling else ''

    for charge in ['plus','minus']:
        binshift = shifts[charge]
        
        total_sig = {};     total_sig_unrolled = {}
        bkg_and_data = {};  bkg_and_data_unrolled = {}
        all_procs = {};     all_procs_unrolled = {};    

        single_sig_unrolled = {}; 
        ratios_unrolled = {}

        for prepost in ['postfit','prefit']:
            suffix = prepost+options.suffix
            # doing signal
            canv = ROOT.TCanvas()
            unpolkeyplot = 'W'+charge+'_unpol_'+prepost
            for ipol,pol in enumerate(['right', 'left', 'long']):
                print "\tPOLARIZATION ",pol
                chpol = '{ch}_{pol}'.format(ch=charge,pol=pol)
                keyplot = 'W'+chpol+'_'+prepost
                nYbins=len(ybins[chpol])-1
                #for ybin in range(1):
                for ybin in range(nYbins): 
                    ykeyplot = 'Y{iy}_{k}'.format(iy=ybin,k=keyplot)
                    group = ybin/groupJobs + 1
                    jobind = min(groupJobs * group - 1,nYbins-1)
                    line = ybin%groupJobs
                    
                    ymin = ybins[chpol][ybin]
                    ymax = ybins[chpol][ybin+1]
         
                    chs = '+' if charge == 'plus' else '-' 
                    h1_1 = infile.Get('expproc_W{ch}_{pol}_Ybin_{ybin}_{sfx}'.format(ch=charge,pol=pol,ybin=ybin,sfx=prepost))
                    name2D = 'W{ch}_{pol}_Ybin_{ybin}'.format(ch=charge,pol=pol,ybin=ybin)
                    title2D = 'W{ch} {pol} : |Y_{{W}}| #in [{ymin},{ymax}]'.format(ymin=ymin,ymax=ymax,pol=pol,ybin=ybin,ch=chs)
                    h1_unrolled = singleChargeUnrolled(h1_1,binshift,nMaskedCha=options.nMaskedChannel)
                    h2_backrolled_1 = dressed2D(copy.deepcopy(h1_unrolled),binning,name2D,title=title2D,shift=binshift,nMaskedCha=options.nMaskedChannel)
                    if options.reverseUnrolling:
                        h1_unrolled = reverseUnroll2Dto1D(h2_backrolled_1,newname=h1_unrolled.GetName()+'_'+reverse)
                    single_sig_unrolled[ykeyplot] = h1_unrolled.Clone('Y{iy}_W{ch}_{pol}_{flav}_{rev}unrolled'.format(iy=ybin,ch=charge,pol=pol,flav=channel,rev=reverse))
                    single_sig_unrolled[ykeyplot].SetDirectory(None)
                    if ybin==0:
                        total_sig[keyplot] = h2_backrolled_1.Clone('W{ch}_{pol}_{flav}'.format(ch=charge,pol=pol,flav=channel))
                        total_sig_unrolled[keyplot] = h1_unrolled.Clone('W{ch}_{pol}_{flav}_{rev}unrolled'.format(ch=charge,pol=pol,flav=channel,rev=reverse))
                    else:
                        total_sig[keyplot].Add(h2_backrolled_1)
                        total_sig_unrolled[keyplot].Add(h1_unrolled)
                    total_sig[keyplot].SetDirectory(None)
                    total_sig_unrolled[keyplot].SetDirectory(None)
                    if prepost=='prefit':
                        ypostfitkey = 'Y{iy}_W{chpol}_postfit'.format(iy=ybin,chpol=chpol) 
                        if single_sig_unrolled[ypostfitkey]!=None:
                            ykeyratio = 'Y{iy}_W{chpol}_ratio'.format(iy=ybin,chpol=chpol)
                            ratios_unrolled[ykeyratio] = single_sig_unrolled[ypostfitkey].Clone(single_sig_unrolled[ypostfitkey].GetName()+'_ratio')
                            ratios_unrolled[ykeyratio].SetDirectory(None)
                            ratios_unrolled[ykeyratio].Divide(single_sig_unrolled[ykeyplot])
                    if not options.no2Dplot:
                        h2_backrolled_1.Draw('colz')
                        h2_backrolled_1.Write(name2D)
                        for ext in ['pdf', 'png']:
                            canv.SaveAs('{odir}/W{ch}_{pol}_W{ch}_{flav}_Ybin_{ybin}_PFMT40_absY_{sfx}.{ext}'.format(odir=outname,ch=charge,flav=channel,pol=pol,ybin=ybin,sfx=suffix,ext=ext))
                        total_sig[keyplot].Draw('colz')
                        total_sig[keyplot].Write(total_sig[keyplot].GetName())
                if ipol==0:
                    total_sig[unpolkeyplot] = total_sig['W'+chpol+'_'+prepost].Clone('W{ch}_unpol_{flav}'.format(ch=charge,flav=channel))
                    total_sig_unrolled[unpolkeyplot] = total_sig_unrolled['W'+chpol+'_'+prepost].Clone('W{ch}_unpol_{flav}_{rev}unrolled'.format(ch=charge,flav=channel,rev=reverse))
                else:
                    total_sig[unpolkeyplot].Add(total_sig['W'+chpol+'_'+prepost])
                    total_sig_unrolled[unpolkeyplot].Add(total_sig_unrolled['W'+chpol+'_'+prepost])
                if prepost=='prefit':
                    postfitkey  = 'W'+chpol+'_postfit'
                    if total_sig_unrolled[postfitkey]!=None:
                        keyratio = 'W'+chpol+'_ratio'
                        ratios_unrolled[keyratio] = total_sig_unrolled[postfitkey].Clone(total_sig_unrolled[postfitkey].GetName()+'_ratio')
                        ratios_unrolled[keyratio].SetDirectory(None)
                        ratios_unrolled[keyratio].Divide(total_sig_unrolled[keyplot])
         
                if not options.no2Dplot:
                    for ext in ['pdf', 'png']:
                        canv.SaveAs('{odir}/W{ch}_{pol}_{flav}_TOTALSIG_PFMT40_absY_{sfx}.{ext}'.format(odir=outname,ch=charge,flav=channel,pol=pol,sfx=suffix,ext=ext))

            if prepost=='prefit':
                unpolpostfitkey = 'W'+charge+'_unpol_postfit'
                if total_sig_unrolled[unpolpostfitkey]!=None:
                    keyratio = 'W'+charge+'_unpol_ratio'
                    ratios_unrolled[keyratio] = total_sig_unrolled[unpolpostfitkey].Clone(total_sig_unrolled[unpolpostfitkey].GetName()+'_ratio')
                    ratios_unrolled[keyratio].SetDirectory(None)
                    ratios_unrolled[keyratio].Divide(total_sig_unrolled[unpolkeyplot])
                                                          
            # do backgrounds now
            procs=["Flips","Z","Top","DiBosons","TauDecaysW","data_fakes","obs"]
            titles=["charge flips","Drell-Yan","Top","di-bosons","W#rightarrow#tau#nu","QCD","data"]
            procsAndTitles = dict(zip(procs,titles))
            for i,p in enumerate(procs):
                keyplot = p if 'obs' in p else p+'_'+prepost
                h1_1 = infile.Get('expproc_{p}_{sfx}'.format(p=p,sfx=prepost)) if 'obs' not in p else infile.Get('obs')
                if not h1_1: continue # muons don't have Flips components
                h1_unrolled =  singleChargeUnrolled(h1_1,binshift,nMaskedCha=options.nMaskedChannel)
                # march2_backrolled_1 = dressed2D(h1_1,binning,p,titles[i],binshift,nMaskedCha=options.nMaskedChannel)
                h2_backrolled_1 = dressed2D(copy.deepcopy(h1_unrolled),binning,p,titles[i],binshift,nMaskedCha=options.nMaskedChannel)
                if options.reverseUnrolling:
                    h1_unrolled = reverseUnroll2Dto1D(h2_backrolled_1,newname=keyplot+charge+'_'+reverse)
                bkg_and_data[keyplot] = h2_backrolled_1;  bkg_and_data_unrolled[keyplot] = h1_unrolled
                bkg_and_data[keyplot].SetDirectory(None); bkg_and_data_unrolled[keyplot].SetDirectory(None)
                if prepost=='prefit' and p!='obs':
                    postfitkey = p+'_postfit'
                    if bkg_and_data_unrolled[postfitkey]!=None:
                        keyratio = p+'_ratio'
                        ratios_unrolled[keyratio] = bkg_and_data_unrolled[postfitkey].Clone(bkg_and_data_unrolled[postfitkey].GetName()+'_ratio')
                        ratios_unrolled[keyratio].SetDirectory(None)
                        ratios_unrolled[keyratio].Divide(bkg_and_data_unrolled[keyplot])
                bkg_and_data_unrolled[keyplot].Write()
                if not options.no2Dplot:
                    h2_backrolled_1.Draw('colz')
                    h2_backrolled_1.Write(str(p))
                    for ext in ['pdf', 'png']:
                        canv.SaveAs('{odir}/{proc}_{ch}_{flav}_PFMT40_absY_{sfx}.{ext}'.format(odir=outname,proc=p,ch=charge,flav=channel,sfx=suffix,ext=ext))
                    
            # now draw the 1D projections
            all_procs.update(total_sig);   all_procs_unrolled.update(total_sig_unrolled);
            all_procs.update(bkg_and_data); all_procs_unrolled.update(bkg_and_data_unrolled); 
            polsym = {'left': 'L', 'right': 'R', 'long': '0'}; chargesym={'plus':'+', 'minus':'-'}
            for pol in ['right', 'left', 'long']:
                procsAndTitles['W{ch}_{pol}'.format(ch=charge,pol=pol)] = 'W^{{{sign}}}_{{{pol}}}'.format(sign=chargesym[charge],pol=polsym[pol])
            procsAndTitles['W{ch}_unpol'.format(ch=charge)] = 'W^{{{sign}}}'.format(sign=chargesym[charge])
         
            colors = {'Wplus_long' : ROOT.kGray+1,   'Wplus_left' : ROOT.kBlue-1,   'Wplus_right' : ROOT.kGreen+1,
                      'Wminus_long': ROOT.kYellow+1, 'Wminus_left': ROOT.kViolet-1, 'Wminus_right': ROOT.kAzure+1,
                      'Top'        : ROOT.kGreen+2,  'DiBosons'   : ROOT.kViolet+2, 'TauDecaysW'  : ROOT.kPink,       
                      'Z'          : ROOT.kAzure+2,  'Flips'      : ROOT.kGray+1,   'data_fakes'  : ROOT.kGray+2,
                      'Wplus_unpol': ROOT.kRed+2,    'Wminus_unpol': ROOT.kRed+2 }
         
            for p,h in all_procs.iteritems():
                h.integral = h.Integral()
     
            # this has the uncertainty propagated with the full covariance matrix
            h1_expfull = infile.Get('expfull_{sfx}'.format(sfx=prepost))
            print '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
            print 'NUMBER OF BINS FOR THE 1D HISTOGRAM', h1_expfull.GetNbinsX()
            h1_expfull_ch = singleChargeUnrolled(h1_expfull,binshift,nMaskedCha=options.nMaskedChannel)
            h2_expfull_backrolled = dressed2D(copy.deepcopy(h1_expfull_ch),binning,'expfull_{ch}_{sfx}'.format(ch=charge,sfx=prepost),'expfull_{ch}_{sfx}'.format(ch=charge,sfx=prepost),binshift,nMaskedCha=options.nMaskedChannel)
            if options.reverseUnrolling:
                h1_expfull_ch = reverseUnroll2Dto1D(h2_expfull_backrolled,newname=h1_expfull.GetName()+'_'+reverse)
            print '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
            print 'nbinsX of 2D histogram:', h2_expfull_backrolled.GetNbinsX()
            print 'nbinsY of 2D histogram:', h2_expfull_backrolled.GetNbinsY()
     
            for projection in ['X','Y']:
                nbinsProj = all_procs['obs'].GetNbinsY() if projection=='X' else all_procs['obs'].GetNbinsX()
                hexpfull = h2_expfull_backrolled.ProjectionX("x_expfull_{sfx}_{ch}".format(sfx=prepost,ch=charge),1,nbinsProj+1,"e") if projection=='X' else h2_expfull_backrolled.ProjectionY("y_expfull_{sfx}_{ch}".format(sfx=prepost,ch=charge),1,nbinsProj+1,"e")
                hdata = all_procs['obs'].ProjectionX("x_data_{ch}".format(ch=charge),1,nbinsProj+1,"e") if projection=='X' else all_procs['obs'].ProjectionY("y_data_{ch}".format(ch=charge),1,nbinsProj+1,"e")
                htot = hdata.Clone('tot_{ch}'.format(ch=charge)); htot.Reset(); htot.Sumw2()
                stack = ROOT.THStack("stack_{sfx}_{ch}".format(sfx=prepost,ch=charge),"")
                
                leg = prepareLegend()
            
                for key,histo in sorted(all_procs.iteritems(), key=lambda (k,v): (v.integral,k)):
                    keycolor = key.replace('_prefit','').replace('_postfit','')
                    if key=='obs' or 'unpol' in key or prepost not in key: continue
                    proj1d = all_procs[key].ProjectionX(all_procs[key].GetName()+charge+"_px",1,nbinsProj+1,"e") if  projection=='X' else all_procs[key].ProjectionY(all_procs[key].GetName()+charge+"_py",1,nbinsProj+1,"e")
                    proj1d.SetFillColor(colors[keycolor])
                    stack.Add(proj1d)
                    htot.Add(proj1d) 
                    leg.AddEntry(proj1d,procsAndTitles[keycolor],'F')
         
                leg.AddEntry(hdata,'data','PE')
                
                hexpfull.GetYaxis().SetRangeUser(0, 1.8*max(htot.GetMaximum(), hdata.GetMaximum()))
                hexpfull.GetYaxis().SetTitle('Events')
             
                plotOne(charge,channel,stack,hexpfull,hdata,leg,outname,'projection%s'%projection,suffix)
            
            # now the unrolled, assign names to run monsterPull.py on them
            hdata_unrolled = singleChargeUnrolled(infile.Get('obs'),binshift,nCharges,nMaskedChanPerCharge).Clone('unrolled_{ch}_data'.format(ch=charge))
            hdata_backrolled_1 = dressed2D(copy.deepcopy(hdata_unrolled),binning,'obs','data2D',binshift,nMaskedCha=options.nMaskedChannel)
            if options.reverseUnrolling:
                hdata_unrolled = reverseUnroll2Dto1D(hdata_backrolled_1,newname=hdata_unrolled.GetName()+'_'+reverse)
            hdata_unrolled.SetDirectory(None)
            htot_unrolled  = hdata_unrolled.Clone('unrolled_{sfx}_{ch}_full'.format(sfx=prepost,ch=charge)); htot_unrolled.Reset(); htot_unrolled.Sumw2()
            htot_unrolled.SetDirectory(None)
            stack_unrolled = ROOT.THStack("stack_unrolled_{sfx}_{ch}".format(sfx=prepost,ch=charge),"") 
            leg_unrolled = prepareLegend(legWidth=0.10)
            output_txt = ''
            print '++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++'
            print 'nbinsX of hdata_unrolled', hdata_unrolled.GetNbinsX()
            for key,histo in sorted(all_procs.iteritems(), key=lambda (k,v): (v.integral,k)):
                keycolor = key.replace('_prefit','').replace('_postfit','')
                if key=='obs' or 'unpol' in key or prepost not in key: continue
                proc_unrolled = all_procs_unrolled[key]
                proc_unrolled.SetFillColor(colors[keycolor])
                stack_unrolled.Add(proc_unrolled)
                htot_unrolled.Add(proc_unrolled)
                leg_unrolled.AddEntry(proc_unrolled,procsAndTitles[keycolor],'F')
                output_txt += '{pname:25s} {nevents:15.2f} +- {nerr:6.2f} \n'.format(pname=key, nevents=histo.Integral(), nerr=0.)
            output_txt += '------------------------\ndata: \t\t {nevents:.0f} \n'.format(nevents=hdata_unrolled.Integral())
            leg_unrolled.AddEntry(hdata_unrolled,'data','PE')
            htot_unrolled.GetYaxis().SetRangeUser(0, 1.8*max(htot_unrolled.GetMaximum(), hdata_unrolled.GetMaximum()))
            htot_unrolled.GetYaxis().SetTitle('Events')
            htot_unrolled.GetXaxis().SetTitle('unrolled lepton (#eta,p_{T}) bin')
            if options.reverseUnrolling:
                plotOne(charge,channel,stack_unrolled,htot_unrolled,hdata_unrolled,leg_unrolled,outname,'unrolled_reverse',suffix,True)
            else:
                plotOne(charge,channel,stack_unrolled,htot_unrolled,hdata_unrolled,leg_unrolled,outname,'unrolled',suffix,True, 
                        drawVertLines="{a},{b}".format(a=recoBins.Npt,b=recoBins.Neta), textForLines=ptBinRanges, outtext=output_txt)
            hdata_unrolled.Write(); htot_unrolled.Write()

        # plot the postfit/prefit ratio
        print "NOW PLOTTING THE RATIOS..."
        for key,histo in ratios_unrolled.iteritems():
            print "Making unrolled ratio for ",key
            outdir = outnamesub if key.startswith('Y') else outname
            if options.reverseUnrolling:
                plotPostFitRatio(charge,channel,histo,outdir,'postfit2prefit_reverse_'+key,options.suffix)
            else:
                plotPostFitRatio(charge,channel,histo,outdir,'postfit2prefit'+reverse+'_'+key,options.suffix, 
                                 drawVertLines="{a},{b}".format(a=recoBins.Npt,b=recoBins.Neta), textForLines=ptBinRanges)
                
    outfile.Close()

    ROOT.gErrorIgnoreLevel = savErrorLevel;

