#!/usr/bin/env python
# USAGE: python postFitPlots.py wel_minus_floatPOI.root cards_el -o outputdir

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

def dressed2D(h1d,binning,name,title=''):
    if len(binning) == 4:
        n1 = binning[0]; bins1 = array('d', binning[1])
        n2 = binning[2]; bins2 = array('d', binning[3])
        h2_1 = ROOT.TH2F(name, title, n1, bins1, n2, bins2 )
    else:
        n1 = binning[0]; min1 = binning[1]; max1 = binning[2]
        n2 = binning[3]; min2 = binning[4]; max2 = binning[5]
        h2_1 = ROOT.TH2F(name, title, n1, min1, max1, n2, min2, max2)
    #h2_backrolled_1 = roll1Dto2D(h1_1, h2_1 )
    h2_backrolled_1 = roll1Dto2D(h1d, h2_1 )
    h2_backrolled_1 .GetXaxis().SetTitle('lepton #eta')
    h2_backrolled_1 .GetYaxis().SetTitle('lepton p_{T} (GeV)')
    h2_backrolled_1 .GetZaxis().SetRangeUser(0.01*h2_backrolled_1.GetMaximum(),1.1*h2_backrolled_1.GetMaximum())
    return h2_backrolled_1

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
    offset = 0.32 if doWide else 0.62
    unity.GetYaxis().SetTitleOffset(offset)
    unity.GetYaxis().SetLabelFont(42)
    unity.GetYaxis().SetLabelSize(0.11)
    unity.GetYaxis().SetLabelOffset(0.007)
    unity.GetYaxis().SetDecimals(True) 
    unity.GetYaxis().SetTitle(ylabel)
    total.GetXaxis().SetLabelOffset(999) ## send them away
    total.GetXaxis().SetTitleOffset(999) ## in outer space
    total.GetYaxis().SetTitleSize(0.06)
    total.GetYaxis().SetTitleOffset(0.75 if doWide else 1.48)
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
    parser.add_option(     '--no2Dplot', dest="no2Dplot", default=True, action='store_true', help="Do not plot templates (but you can still save them in a root file with option -s)");
    parser.add_option('-s','--save', dest='outfile_templates', default='templates_2D', type='string', help='pass name of output file to save 2D histograms (charge is automatically appended before extension). No need to specify extension, .root is automatically addedx')
    parser.add_option(     '--prefit', dest="prefit", default=False, action='store_true', help="");
    groupJobs=5 # used in make_helicity_cards.py

    (options, args) = parser.parse_args()
    if len(args) < 1:
        parser.print_usage()
        quit()

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
    charge = 'plus' if any(['Wplus' in k.GetName() for k in infile.GetListOfKeys()]) else 'minus'
    print "From the list of histograms it seems that you are plotting results for channel ",channel," and charge ",charge

    outfile_templates = options.outfile_templates
    if outfile_templates.endswith(".root"):
        outfile_templates = outfile_templates.replace(".root","_%s.root" % str(charge))
    else:
        outfile_templates = outfile_templates + "_" + str(charge) + ".root"

    full_outfileName = outname + "/" + outfile_templates
    outfile = ROOT.TFile(full_outfileName, 'recreate')
    print "Will save 2D templates in file --> " + full_outfileName

    # doing signal
    total_sig = {}
    canv = ROOT.TCanvas()
    for pol in ['right', 'left', 'long']:
        print "\tPOLARIZATION ",pol
        chpol = '{ch}_{pol}'.format(ch=charge,pol=pol)
        nYbins=len(ybins[chpol])-1
        #for ybin in range(1):
        for ybin in range(nYbins): 

            group = ybin/groupJobs + 1
            jobind = min(groupJobs * group - 1,nYbins-1)
            line = ybin%groupJobs
            

            jobsdir = args[1]+'/jobs/'
            jobfile_name = 'W{ch}_{pol}_{flav}_Ybin_{b}.sh'.format(ch=charge,pol=pol,flav=channel,b=jobind)
            tmp_jobfile = open(jobsdir+jobfile_name, 'r')
            lineno = (-groupJobs+line)*2 + 1 # white lines
            if ybin==nYbins-1: lineno = -1 # not general hack!!!
            tmp_line = tmp_jobfile.readlines()[lineno].split()
            ymin = list(i for i in tmp_line if '(genw_y)>' in i)[0].replace('\'','').split('>')[-1]
            ymax = list(i for i in tmp_line if '(genw_y)<' in i)[0].replace('\'','').split('<')[-1]
            
            binning = getbinning(tmp_line)

            chs = '+' if charge == 'plus' else '-' 
            h1_1 = infile.Get('expproc_W{ch}_{pol}_W{ch}_{pol}_{flav}_Ybin_{ybin}_postfit'.format(ch=charge,pol=pol,flav=channel,ybin=ybin))
            name2D = 'W{ch}_{pol}_W{ch}_{pol}_{flav}_Ybin_{ybin}'.format(ch=charge,pol=pol,flav=channel,ybin=ybin)
            title2D = 'W{ch} {pol} : |Y_{{W}}| #in [{ymin},{ymax}]'.format(ymin=ymin,ymax=ymax,pol=pol,ybin=ybin,ch=chs)
            h2_backrolled_1 = dressed2D(h1_1,binning,name2D,title2D)
            if ybin==0:
                total_sig["W"+chpol] = h2_backrolled_1.Clone('W{ch}_{pol}_{flav}'.format(ch=charge,pol=pol,flav=channel))
            else:
                total_sig["W"+chpol].Add(h2_backrolled_1)
            total_sig["W"+chpol].SetDirectory(None)
            h2_backrolled_1.Draw('colz')
            h2_backrolled_1.Write(name2D)
            if not options.no2Dplot:
                for ext in ['pdf', 'png']:
                    canv.SaveAs('{odir}/W{ch}_{pol}_W{ch}_{flav}_Ybin_{ybin}_PFMT40_absY.{ext}'.format(odir=outname,ch=charge,flav=channel,pol=pol,ybin=ybin,ext=ext))
        total_sig["W"+chpol].Draw('colz')
        total_sig["W"+chpol].Write(total_sig["W"+chpol].GetName())

        if not options.no2Dplot:
            for ext in ['pdf', 'png']:
                canv.SaveAs('{odir}/W{ch}_{pol}_{flav}_TOTALSIG_PFMT40_absY.{ext}'.format(odir=outname,ch=charge,flav=channel,pol=pol,ext=ext))


    # do backgrounds now
    procs=["Flips","Z","Top","DiBosons","TauDecaysW","data_fakes","obs"]
    titles=["charge flips","Drell-Yan","Top","di-bosons","W#rightarrow#tau#nu","QCD","toy data"]
    procsAndTitles = dict(zip(procs,titles))
    bkg_and_data = {}
    for i,p in enumerate(procs):
        h1_1 = infile.Get('expproc_{p}_postfit'.format(p=p)) if 'obs' not in p else infile.Get('obs')
        if not h1_1: continue # muons don't have Flips components
        h2_backrolled_1 = dressed2D(h1_1,binning,p,titles[i])
        bkg_and_data[p] = h2_backrolled_1
        bkg_and_data[p].SetDirectory(None)
        h2_backrolled_1.Draw('colz')
        h2_backrolled_1.Write(str(p))
        if not options.no2Dplot:
            for ext in ['pdf', 'png']:
                canv.SaveAs('{odir}/{proc}_{ch}_{flav}_PFMT40_absY.{ext}'.format(odir=outname,proc=p,ch=charge,flav=channel,ext=ext))
        
    outfile.Close()
    

    # now draw the 1D projections
    all_procs = total_sig.copy()
    all_procs.update(bkg_and_data)
    polsym = {'left': 'L', 'right': 'R', 'long': '0'}; chargesym={'plus':'+', 'minus':'-'}
    for pol in ['right', 'left', 'long']:
        procsAndTitles['W{ch}_{pol}'.format(ch=charge,pol=pol)] = 'W^{{{sign}}}_{{{pol}}}'.format(sign=chargesym[charge],pol=polsym[pol])

    colors = {'Wplus_long' : ROOT.kGray+1,   'Wplus_left' : ROOT.kBlue-1,   'Wplus_right' : ROOT.kGreen+1,
              'Wminus_long': ROOT.kYellow+1, 'Wminus_left': ROOT.kViolet-1, 'Wminus_right': ROOT.kAzure+1,
              'Top'        : ROOT.kGreen+2,  'DiBosons'   : ROOT.kViolet+2, 'TauDecaysW'  : ROOT.kPink,       
              'Z'          : ROOT.kAzure+2,  'Flips'      : ROOT.kGray+1,   'data_fakes'  : ROOT.kGray+2 }

    for p,h in all_procs.iteritems():
        h.integral = h.Integral()

    for projection in ['X','Y']:
        hdata = all_procs['obs'].ProjectionX("x_total",0,-1,"e") if projection=='X' else all_procs['obs'].ProjectionY("y_total",0,-1,"e")
        htot = hdata.Clone('tot'); htot.Reset(); htot.Sumw2()
        stack = ROOT.THStack("stack","")
        
        legWidth=0.50; textSize=0.035
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
    
        for key,histo in sorted(all_procs.iteritems(), key=lambda (k,v): (v.integral,k)):
            if key=='obs': continue
            proj1d = all_procs[key].ProjectionX(all_procs[key].GetName()+"_px",0,-1,"e") if  projection=='X' else all_procs[key].ProjectionY(all_procs[key].GetName()+"_py",0,-1,"e")
            proj1d.SetFillColor(colors[key])
            stack.Add(proj1d)
            htot.Add(proj1d)                    
            leg.AddEntry(proj1d,procsAndTitles[key],'F')

        leg.AddEntry(hdata,'toy data','PE')
        
        doRatio=True
        htot.GetYaxis().SetRangeUser(0, 1.8*max(htot.GetMaximum(), hdata.GetMaximum()))
        htot.GetYaxis().SetTitle('Events')
     
    
        ## Prepare split screen
        c1 = ROOT.TCanvas("c1", "c1", 600, 750); c1.Draw()
        c1.SetWindowSize(600 + (600 - c1.GetWw()), (750 + (750 - c1.GetWh())));
        ROOT.gStyle.SetPadLeftMargin(0.18)
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
        htot.Draw("HIST")
        #htot.SetLabelOffset(9999.0);
        #htot.SetTitleOffset(9999.0);
        stack.Draw("HIST F SAME")
        hdata.SetMarkerColor(ROOT.kBlack)
        hdata.SetMarkerStyle(ROOT.kFullCircle)
        hdata.Draw("E SAME")
        htot.Draw("AXIS SAME")
        totalError = doShadedUncertainty(htot)            
        leg.Draw()
        lat = ROOT.TLatex()
        lat.SetNDC(); lat.SetTextFont(42)
        lat.DrawLatex(0.16, 0.92, '#bf{CMS} #it{Preliminary}')
        lat.DrawLatex(0.65, 0.92, '36 fb^{-1} (13 TeV)')

        p2.cd()
        rdata,rnorm,rline = doRatioHists(hdata, htot, maxRange=[0.99,1.01], fixRange=True)
        c1.cd()
        for ext in ['pdf', 'png']:
            c1.SaveAs('{odir}/projection{proj}_{ch}_{flav}_PFMT40_absY.{ext}'.format(odir=outname,proj=projection,ch=charge,flav=channel,ext=ext))
        
