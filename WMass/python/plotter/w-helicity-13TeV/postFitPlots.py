#!/usr/bin/env python
# USAGE: python postFitPlots.py wel_minus_floatPOI.root cards_el -o outputdir [--prefit]

import ROOT, os, re
from array import array
from CMGTools.WMass.plotter.mcPlots import doShadedUncertainty

import utilities
utilities = utilities.util()

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
ROOT.gStyle.SetPadRightMargin(0.13)

def roll1Dto2D(h1d, histo):  #,h2dname):#,plotfile,options):
    for i in xrange(1,h1d.GetNbinsX()+1):
        xbin = i % histo.GetNbinsX()
        if not xbin: xbin = xbin+histo.GetNbinsX()
        ybin = i / histo.GetNbinsX() + (1 if i%histo.GetNbinsX() else 0)
        val = h1d.GetBinContent(i)
        histo.SetBinContent(xbin,ybin,h1d.GetBinContent(i))
        histo.SetBinError(xbin,ybin,h1d.GetBinError(i))
    return histo

def unroll2Dto1D(h):
    nbins = h.GetNbinsX() * h.GetNbinsY()
    goodname = h.GetName()
    h.SetName(goodname+"_oldbinning")
    newh = ROOT.TH1D(goodname,h.GetTitle(),nbins,0.5,nbins+0.5)
    newh.Sumw2()
    if 'TH2' not in h.ClassName(): raise RuntimeError, "Calling rebin2Dto1D on something that is not TH2"
    for i in xrange(h.GetNbinsX()):
        for j in xrange(h.GetNbinsY()):
            bin = 1 + i + j*h.GetNbinsX()
            newh.SetBinContent(bin,h.GetBinContent(i+1,j+1))
            newh.SetBinError(bin,h.GetBinError(i+1,j+1))
    for bin in range(1,nbins+1):
        if newh.GetBinContent(bin)<0:
            print 'Warning: cropping to zero bin %d in %s (was %f)'%(bin,newh.GetName(),newh.GetBinContent(bin))
            newh.SetBinContent(bin,0)
    newh.SetLineWidth(h.GetLineWidth())
    newh.SetLineStyle(h.GetLineStyle())
    newh.SetLineColor(h.GetLineColor())
    return newh

def dressed2D(h1d,binning,name,title='',shift=0,nCharges=2,nMaskedCha=2):
    if len(binning) == 4:
        n1 = binning[0]; bins1 = array('d', binning[1])
        n2 = binning[2]; bins2 = array('d', binning[3])
        h2_1 = ROOT.TH2F(name, title, n1, bins1, n2, bins2 )
    else:
        n1 = binning[0]; min1 = binning[1]; max1 = binning[2]
        n2 = binning[3]; min2 = binning[4]; max2 = binning[5]
        h2_1 = ROOT.TH2F(name, title, n1, min1, max1, n2, min2, max2)
    h1d_shifted = singleChargeUnrolled(h1d,shift,nCharges,nMaskedCha)
    h2_backrolled_1 = roll1Dto2D(h1d_shifted, h2_1 )
    h2_backrolled_1 .GetXaxis().SetTitle('lepton #eta')
    h2_backrolled_1 .GetYaxis().SetTitle('lepton p_{T} (GeV)')
    h2_backrolled_1 .GetZaxis().SetRangeUser(0.01*h2_backrolled_1.GetMaximum(),1.1*h2_backrolled_1.GetMaximum())
    return h2_backrolled_1

def chargeUnrolledBinShifts(infile,channel,nCharges=2,nMaskedCha=2):
    # guess from a signal with charge defined name
    h1d = infile.Get('expproc_Wplus_left_Wplus_left_{ch}_Ybin_0_postfit'.format(ch=channel))
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
    h1d_shifted = ROOT.TH1D('shift','',nbins,1,nbins)
    for b in xrange(nbins):
        h1d_shifted.SetBinContent(b+1,h1d.GetBinContent(b+shift+1))
        h1d_shifted.SetBinError(b+1,h1d.GetBinError(b+shift+1))
    return h1d_shifted

def prepareLegend(legWidth=0.50,textSize=0.035):
    (x1,y1,x2,y2) = (.75-legWidth, .70, .85, .87)
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

def plotOne(charge,channel,stack,htot,hdata,legend,outdir,prefix,suffix,veryWide=False):
    ## Prepare split screen
    plotformat = (2400,600) if veryWide else (600,750)
    c1 = ROOT.TCanvas("c1", "c1", plotformat[0], plotformat[1]); c1.Draw()
    c1.SetWindowSize(plotformat[0] + (plotformat[0] - c1.GetWw()), (plotformat[1] + (plotformat[1] - c1.GetWh())));
    ROOT.gStyle.SetPadLeftMargin(0.07 if veryWide else 0.18)
    ROOT.gStyle.SetPadRightMargin(0.07 if veryWide else 0.13)
    p1 = ROOT.TPad("pad1","pad1",0,0.29,1,0.99);
    p1.SetBottomMargin(0.03);
    p1.Draw();
    p2 = ROOT.TPad("pad2","pad2",0,0,1,0.31);
    p2.SetTopMargin(0);
    p2.SetBottomMargin(0.3);
    p2.SetFillStyle(0);
    p2.Draw();
    p1.cd();
    ## Draw absolute prediction in top frame
    offset = 0.45 if veryWide else 0.62
    htot.GetYaxis().SetTitleOffset(offset)
    htot.Draw("HIST")
    #htot.SetLabelOffset(9999.0);
    #htot.SetTitleOffset(9999.0);
    stack.Draw("HIST F SAME")
    hdata.SetMarkerColor(ROOT.kBlack)
    hdata.SetMarkerStyle(ROOT.kFullCircle)
    hdata.Draw("E SAME")
    htot.Draw("AXIS SAME")
    totalError = doShadedUncertainty(htot)            
    legend.Draw()
    lat = ROOT.TLatex()
    lat.SetNDC(); lat.SetTextFont(42)
    if veryWide:
        x1 = 0.09; x2 = 0.85
    else:
        x1 = 0.16; x2 = 0.65
    lat.DrawLatex(x1, 0.92, '#bf{CMS} #it{Preliminary}')
    lat.DrawLatex(x2, 0.92, '36 fb^{-1} (13 TeV)')
    
    p2.cd()
    maxrange = [0.95,1.05] if prepost == 'postfit' else [0.90,1.10]
    rdata,rnorm,rline = doRatioHists(hdata, htot, maxRange=maxrange, fixRange=True, doWide=veryWide)
    c1.cd()
    for ext in ['pdf', 'png']:
        c1.SaveAs('{odir}/{pfx}_{ch}_{flav}_PFMT40_absY_{sfx}.{ext}'.format(odir=outdir,pfx=prefix,ch=charge,flav=channel,sfx=suffix,ext=ext))


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
    parser.add_option(     '--prefit', dest="prefit", default=False, action='store_true', help="");
    parser.add_option(     '--suffix', dest="suffix", default='', type='string', help="define suffix for each plot");
    groupJobs=5 # used in make_helicity_cards.py
    nCharges = 2; nMaskedChanPerCharge = 2; # edm hardcoded, but to be guessed 
    
    (options, args) = parser.parse_args()
    if len(args) < 1:
        parser.print_usage()
        quit()

    prepost = 'prefit' if options.prefit else 'postfit'
    suffix = prepost+options.suffix

    ## get the binning of the YW bins
    ybinfile = args[1]+'/binningYW.txt'
    ybinfile = open(ybinfile, 'r')
    ybins = eval(ybinfile.read())
    ybinfile.close()

    outname = options.outdir
    if not os.path.exists(outname):
        os.system("mkdir -p "+outname)
        if os.path.exists("/afs/cern.ch"): os.system("cp /afs/cern.ch/user/g/gpetrucc/php/index.php "+outname)

    infile = ROOT.TFile(args[0], 'read')
    channel = 'el' if any(['el_Ybin' in k.GetName() for k in infile.GetListOfKeys()]) else 'mu'
    print "From the list of histograms it seems that you are plotting results for channel ",channel

    full_outfileName = '{odir}/plots_{sfx}.root'.format(odir=outname,sfx=suffix)
    outfile = ROOT.TFile(full_outfileName, 'recreate')
    print "Will save 2D templates in file --> " + full_outfileName

    ## get the pT and eta binning from the file in the directory
    binninPtEtaFile = open(args[1]+'/binningPtEta.txt','r')
    bins = binninPtEtaFile.readlines()[1].split()[1]
    ## hack. easier
    etabins = list( float(i) for i in bins.replace(' ','').split('*')[0].replace('[','').replace(']','').split(',') )
    ptbins  = list( float(i) for i in bins.replace(' ','').split('*')[1].replace('[','').replace(']','').split(',') )
    nbinseta = len(etabins)-1
    nbinspt  = len( ptbins)-1
    binning = [nbinseta, etabins, nbinspt, ptbins]

    shifts = chargeUnrolledBinShifts(infile,channel,nCharges,nMaskedChanPerCharge)

    savErrorLevel = ROOT.gErrorIgnoreLevel; ROOT.gErrorIgnoreLevel = ROOT.kWarning;
    # doing signal
    for charge in ['plus','minus']:
        total_sig = {}
        total_sig_unrolled = {}
        canv = ROOT.TCanvas()
        binshift = shifts[charge]
        for pol in ['right', 'left', 'long']:
            print "\tPOLARIZATION ",pol
            chpol = '{ch}_{pol}'.format(ch=charge,pol=pol)
            nYbins=len(ybins[chpol])-1
            #for ybin in range(1):
            for ybin in range(nYbins): 
     
                group = ybin/groupJobs + 1
                jobind = min(groupJobs * group - 1,nYbins-1)
                line = ybin%groupJobs
                
                ymin = ybins[chpol][ybin]
                ymax = ybins[chpol][ybin+1]
     
                chs = '+' if charge == 'plus' else '-' 
                h1_1 = infile.Get('expproc_W{ch}_{pol}_W{ch}_{pol}_{flav}_Ybin_{ybin}_{sfx}'.format(ch=charge,pol=pol,flav=channel,ybin=ybin,sfx=prepost))
                name2D = 'W{ch}_{pol}_W{ch}_{pol}_{flav}_Ybin_{ybin}'.format(ch=charge,pol=pol,flav=channel,ybin=ybin)
                title2D = 'W{ch} {pol} : |Y_{{W}}| #in [{ymin},{ymax}]'.format(ymin=ymin,ymax=ymax,pol=pol,ybin=ybin,ch=chs)
                h2_backrolled_1 = dressed2D(h1_1,binning,name2D,title2D,binshift)
                h1_unrolled = singleChargeUnrolled(h1_1,binshift)
                if ybin==0:
                    total_sig["W"+chpol] = h2_backrolled_1.Clone('W{ch}_{pol}_{flav}'.format(ch=charge,pol=pol,flav=channel))
                    total_sig_unrolled["W"+chpol] = h1_unrolled.Clone('W{ch}_{pol}_{flav}_unrolled'.format(ch=charge,pol=pol,flav=channel))
                else:
                    total_sig["W"+chpol].Add(h2_backrolled_1)
                    total_sig_unrolled["W"+chpol].Add(h1_unrolled)
                total_sig["W"+chpol].SetDirectory(None)
                total_sig_unrolled["W"+chpol].SetDirectory(None)
                if not options.no2Dplot:
                    h2_backrolled_1.Draw('colz')
                    h2_backrolled_1.Write(name2D)
                    for ext in ['pdf', 'png']:
                        canv.SaveAs('{odir}/W{ch}_{pol}_W{ch}_{flav}_Ybin_{ybin}_PFMT40_absY_{sfx}.{ext}'.format(odir=outname,ch=charge,flav=channel,pol=pol,ybin=ybin,sfx=suffix,ext=ext))
                    total_sig["W"+chpol].Draw('colz')
                    total_sig["W"+chpol].Write(total_sig["W"+chpol].GetName())
     
            if not options.no2Dplot:
                for ext in ['pdf', 'png']:
                    canv.SaveAs('{odir}/W{ch}_{pol}_{flav}_TOTALSIG_PFMT40_absY_{sfx}.{ext}'.format(odir=outname,ch=charge,flav=channel,pol=pol,sfx=suffix,ext=ext))


        # do backgrounds now
        procs=["Flips","Z","Top","DiBosons","TauDecaysW","data_fakes","obs"]
        titles=["charge flips","Drell-Yan","Top","di-bosons","W#rightarrow#tau#nu","QCD","data"]
        procsAndTitles = dict(zip(procs,titles))
        bkg_and_data = {}
        bkg_and_data_unrolled = {}
        for i,p in enumerate(procs):
            h1_1 = infile.Get('expproc_{p}_{sfx}'.format(p=p,sfx=prepost)) if 'obs' not in p else infile.Get('obs')
            h1_unrolled =  singleChargeUnrolled(h1_1,binshift)
            if not h1_1: continue # muons don't have Flips components
            h2_backrolled_1 = dressed2D(h1_1,binning,p,titles[i],binshift)
            bkg_and_data[p] = h2_backrolled_1;  bkg_and_data_unrolled[p] = h1_unrolled
            bkg_and_data[p].SetDirectory(None); bkg_and_data_unrolled[p].SetDirectory(None)
            if not options.no2Dplot:
                h2_backrolled_1.Draw('colz')
                h2_backrolled_1.Write(str(p))
                for ext in ['pdf', 'png']:
                    canv.SaveAs('{odir}/{proc}_{ch}_{flav}_PFMT40_absY_{sfx}.{ext}'.format(odir=outname,proc=p,ch=charge,flav=channel,sfx=suffix,ext=ext))
                
        ROOT.gErrorIgnoreLevel = savErrorLevel;

        # now draw the 1D projections
        all_procs = total_sig.copy();   all_procs_unrolled = total_sig_unrolled.copy();
        all_procs.update(bkg_and_data); all_procs_unrolled.update(bkg_and_data_unrolled); 
        polsym = {'left': 'L', 'right': 'R', 'long': '0'}; chargesym={'plus':'+', 'minus':'-'}
        for pol in ['right', 'left', 'long']:
            procsAndTitles['W{ch}_{pol}'.format(ch=charge,pol=pol)] = 'W^{{{sign}}}_{{{pol}}}'.format(sign=chargesym[charge],pol=polsym[pol])
     
        colors = {'Wplus_long' : ROOT.kGray+1,   'Wplus_left' : ROOT.kBlue-1,   'Wplus_right' : ROOT.kGreen+1,
                  'Wminus_long': ROOT.kYellow+1, 'Wminus_left': ROOT.kViolet-1, 'Wminus_right': ROOT.kAzure+1,
                  'Top'        : ROOT.kGreen+2,  'DiBosons'   : ROOT.kViolet+2, 'TauDecaysW'  : ROOT.kPink,       
                  'Z'          : ROOT.kAzure+2,  'Flips'      : ROOT.kGray+1,   'data_fakes'  : ROOT.kGray+2 }
     
        for p,h in all_procs.iteritems():
            h.integral = h.Integral()

        # this has the uncertainty propagated with the full covariance matrix
        h1_expfull = infile.Get('expfull_{sfx}'.format(sfx=prepost))
        h2_expfull_backrolled = dressed2D(h1_expfull,binning,'expfull','expfull',binshift)

        for projection in ['X','Y']:
            nbinsProj = all_procs['obs'].GetNbinsY() if projection=='X' else all_procs['obs'].GetNbinsX()
            hexpfull = h2_expfull_backrolled.ProjectionX("x_expfull_{ch}".format(ch=charge),1,nbinsProj,"e") if projection=='X' else h2_expfull_backrolled.ProjectionY("y_expfull_{ch}".format(ch=charge),1,nbinsProj,"e")
            hdata = all_procs['obs'].ProjectionX("x_data_{ch}".format(ch=charge),1,nbinsProj,"e") if projection=='X' else all_procs['obs'].ProjectionY("y_data_{ch}".format(ch=charge),1,nbinsProj,"e")
            htot = hdata.Clone('tot'); htot.Reset(); htot.Sumw2()
            stack = ROOT.THStack("stack","")
            
            leg = prepareLegend()
        
            for key,histo in sorted(all_procs.iteritems(), key=lambda (k,v): (v.integral,k)):
                if key=='obs': continue
                proj1d = all_procs[key].ProjectionX(all_procs[key].GetName()+"_px",0,-1,"e") if  projection=='X' else all_procs[key].ProjectionY(all_procs[key].GetName()+"_py",0,-1,"e")
                proj1d.SetFillColor(colors[key])
                stack.Add(proj1d)
                htot.Add(proj1d) 
                leg.AddEntry(proj1d,procsAndTitles[key],'F')
     
            leg.AddEntry(hdata,'data','PE')
            
            hexpfull.GetYaxis().SetRangeUser(0, 1.8*max(htot.GetMaximum(), hdata.GetMaximum()))
            hexpfull.GetYaxis().SetTitle('Events')
         
            plotOne(charge,channel,stack,hexpfull,hdata,leg,outname,'projection%s'%projection,suffix)
        
        # now the unrolled, assign names to run monsterPull.py on them
        hdata_unrolled = singleChargeUnrolled(infile.Get('obs'),binshift,nCharges,nMaskedChanPerCharge).Clone('unrolled_{ch}_data'.format(ch=charge))
        hdata_unrolled.SetDirectory(None)
        htot_unrolled  = hdata_unrolled.Clone('unrolled_{ch}_full'.format(ch=charge)); htot_unrolled.Reset(); htot_unrolled.Sumw2()
        htot_unrolled.SetDirectory(None)
        stack_unrolled = ROOT.THStack("stack_unrolled","") 
        leg_unrolled = prepareLegend(legWidth=0.10)
        for key,histo in sorted(all_procs.iteritems(), key=lambda (k,v): (v.integral,k)):
            if key=='obs': continue
            proc_unrolled = all_procs_unrolled[key]
            proc_unrolled.SetFillColor(colors[key])
            stack_unrolled.Add(proc_unrolled)
            htot_unrolled.Add(proc_unrolled)
            leg_unrolled.AddEntry(proc_unrolled,procsAndTitles[key],'F')
        leg_unrolled.AddEntry(hdata_unrolled,'data','PE')
        htot_unrolled.GetYaxis().SetRangeUser(0, 1.8*max(htot_unrolled.GetMaximum(), hdata_unrolled.GetMaximum()))
        htot_unrolled.GetYaxis().SetTitle('Events')
        htot_unrolled.GetXaxis().SetTitle('unrolled lepton (#eta,p_{T}) bin')
        plotOne(charge,channel,stack_unrolled,htot_unrolled,hdata_unrolled,leg_unrolled,outname,'unrolled',suffix,True)
        hdata_unrolled.Write(); htot_unrolled.Write()

    outfile.Close()
