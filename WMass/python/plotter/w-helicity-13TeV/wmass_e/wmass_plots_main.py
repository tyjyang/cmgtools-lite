import optparse, subprocess, ROOT, datetime, math, array, copy, os, itertools, copy
ROOT.gROOT.SetBatch(True)

POLARIZATIONS = ['left','right','long']

def doLegend(histos,labels,styles,corner="TR",textSize=0.035,legWidth=0.18,legBorder=False,nColumns=1):
    nentries = len(histos)
    (x1,y1,x2,y2) = (.85-legWidth, .7 - textSize*max(nentries-3,0), .90, .91)
    if corner == "TR":
        (x1,y1,x2,y2) = (.85-legWidth, .7 - textSize*max(nentries-3,0), .90, .91)
    elif corner == "TC":
        (x1,y1,x2,y2) = (.5, .75 - textSize*max(nentries-3,0), .5+legWidth, .91)
    elif corner == "TL":
        (x1,y1,x2,y2) = (.2, .75 - textSize*max(nentries-3,0), .2+legWidth, .91)
    elif corner == "BR":
        (x1,y1,x2,y2) = (.85-legWidth, .33 + textSize*max(nentries-3,0), .90, .15)
    elif corner == "BC":
        (x1,y1,x2,y2) = (.5, .33 + textSize*max(nentries-3,0), .5+legWidth, .35)
    elif corner == "BL":
        (x1,y1,x2,y2) = (.2, .33 + textSize*max(nentries-3,0), .33+legWidth, .35)
    leg = ROOT.TLegend(x1,y1,x2,y2)
    leg.SetNColumns(nColumns)
    leg.SetFillColor(0)
    leg.SetFillColorAlpha(0,0.6)  # should make the legend semitransparent (second number is 0 for fully transparent, 1 for full opaque)
    #leg.SetFillStyle(0) # transparent legend, so it will not cover plots (markers of legend entries will cover it unless one changes the histogram FillStyle, but this has other effects on color, so better not touching the FillStyle)
    leg.SetShadowColor(0)
    if not legBorder:
        leg.SetLineColor(0)
        leg.SetBorderSize(0)  # remove border  (otherwise it is drawn with a white line, visible if it overlaps with plots
    leg.SetTextFont(42)
    leg.SetTextSize(textSize)
    for (plot,label,style) in zip(histos,labels,styles): leg.AddEntry(plot,label,style)
    leg.Draw()
    ## assign it to a global variable so it's not deleted
    global legend_
    legend_ = leg 
    return leg

def setRootStyle(th2=False):
    ROOT.gStyle.SetPadTopMargin(0.0)
    ROOT.gStyle.SetPadBottomMargin(0.0)
    ROOT.gStyle.SetPadLeftMargin(0.18)
    if th2:
        ROOT.gStyle.SetPadLeftMargin(0.16)                        
    ROOT.gStyle.SetTitleXOffset(1.5)
    ROOT.gStyle.SetTitleYOffset(1.1)

def kappaPtW(pt):
    offset = 0.044; slope = 0.005
    return offset + slope*math.sqrt(pt)

def wptPostFitRatios(options):
    maindir = '/eos/home-e/emanuele/www/Analysis/WMass/13TeV/plots/gen/wpt/'

    colors = {('plus','left'):ROOT.kBlue-1,('plus','right'):ROOT.kGreen+1,('plus','long'):ROOT.kGray+1,
              ('minus','left'):ROOT.kViolet-1,('minus','right'):ROOT.kAzure+1,('minus','long'):ROOT.kYellow+1,}
    
    setRootStyle()
    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetHatchesLineWidth(1) 
    ROOT.gStyle.SetHatchesSpacing(0.5)

    def makeRatio(hprefit,hpostfit,color):
        hpostfit.Scale(hprefit.Integral()/hpostfit.Integral())
        ratio = copy.deepcopy(hpostfit.Clone(hpostfit.GetName()+'_ratio'))
        ratio.Divide(hprefit)
        for b in xrange(ratio.GetNbinsX()):
            pt = ratio.GetXaxis().GetBinCenter(b+1)
            kappa = kappaPtW(pt)
            err_ratio = math.hypot(kappa,ratio.GetBinError(b+1))
            ratio.SetBinError(b+1,err_ratio)
            if ratio.GetBinContent(b+1)!=0:
                hpostfit.SetBinError(b+1,math.hypot(hpostfit.GetBinError(b+1),hpostfit.GetBinContent(b+1)*err_ratio/ratio.GetBinContent(b+1)))
                hprefit.SetBinError(b+1,math.hypot(hprefit.GetBinError(b+1),hprefit.GetBinContent(b+1)*kappa))
            else:
                hpostfit.SetBinError(b+1,0)
                hprefit.SetBinError(b+1,0)

        ratio.SetFillColorAlpha(color,0.20)
        ratio.SetLineColor(color)
        ratio.SetMarkerColor(color)
        ratio.GetYaxis().SetTitle('postfit / prefit')
        ratio.GetXaxis().SetTitleOffset(1.5)
        ratio.SetTitle('')
        ratio.GetYaxis().SetRangeUser(0.8,1.2)
        
        for h in [hpostfit,hprefit]:
            h.SetTitle('')
            h.SetMarkerSize(0)
            h.SetLineWidth(2)

        if 'wpt' in ratio.GetName():
            hprefit.GetXaxis().SetRangeUser(0,99)
            hpostfit.GetXaxis().SetRangeUser(0,99)
            ratio.GetXaxis().SetRangeUser(0,99)
            ratio.GetXaxis().SetTitle('p_{T}^{W} (GeV)')
        elif 'wy' in ratio.GetName():
            ratio.GetXaxis().SetRangeUser(-5,5)
            ratio.GetXaxis().SetTitle('Y_{W}')
        return ratio

    for var in ['wpt']: #['wpt','wy']:
        plots_unpol = {}
        c = ROOT.TCanvas('c','',1200,1200)
        charges = ['plus','minus']
        signs = {'plus':'+', 'minus':'-'}
        for charge in charges:
            if var=='wpt': c.SetLogx(1)
            else: c.SetLogx(0)
            rfiles = {}
            rfiles['prefit']  = ROOT.TFile.Open('{d}/wgen_nosel/w{charge}_{var}.root'.format(d=maindir,charge=charge,var=var))
            rfiles['postfit'] = ROOT.TFile.Open('{d}/wgen_nosel_qcdpostfit_barolo_fitel/w{charge}_{var}.root'.format(d=maindir,charge=charge,var=var))
            histos = []; labels = []; styles = []
            for i,pol in enumerate(POLARIZATIONS):
                prefit  = rfiles['prefit'] .Get('w{charge}_{var}_W{charge}_{pol}'.format(charge=charge,var=var,pol=pol))
                postfit = rfiles['postfit'].Get('w{charge}_{var}_W{charge}_{pol}'.format(charge=charge,var=var,pol=pol))
                if i==0:
                    prefit_unpol  = prefit .Clone(prefit .GetName()+'_unpol')
                    postfit_unpol = postfit.Clone(postfit.GetName()+'_unpol')
                else:
                    prefit_unpol .Add(prefit)
                    postfit_unpol.Add(postfit)
                ratio = makeRatio(prefit,postfit,colors[(charge,pol)])
                ratio.Draw('E3' if i==0 else 'E3 same')
                histos.append(ratio)
                labels.append('W{sign} {pol}'.format(sign=signs[charge],pol=pol))
                styles.append('pf')
            ## draw polarized
            legPol = doLegend(histos,labels,styles,corner="BL")
            legPol.Draw()
            for ext in ['pdf','png']:
                c.Print('{pdir}/{var}_postOverPrefit_{charge}.{ext}'.format(pdir=options.printDir,var=var,charge=charge,pol=pol,ext=ext))
            legPol.Clear()
            ## draw unpolarized
            ratio_unpol = makeRatio(prefit_unpol,postfit_unpol,ROOT.kRed+1)
            plots_unpol['{var}_w{charge}_prefit'.format(var=var,charge=charge)] = copy.deepcopy(prefit_unpol)
            plots_unpol['{var}_w{charge}_postfit'.format(var=var,charge=charge)] = copy.deepcopy(postfit_unpol)
            plots_unpol['{var}_w{charge}_ratio'.format(var=var,charge=charge)] = copy.deepcopy(ratio_unpol)
            ratio_unpol.Draw('E3')
            legUnpol = doLegend([ratio_unpol],['W{sign}'.format(sign=signs[charge])],['pf'],corner="BL")
            legUnpol.Draw()
            for ext in ['pdf','png']:
                c.Print('{pdir}/{var}_postOverPrefit_{charge}_unpol.{ext}'.format(pdir=options.printDir,var=var,charge=charge,pol=pol,ext=ext))        
            legUnpol.Clear()
        c.Clear()
        ## now plot the unpolarized W PT/Y on top and ratio in the bottom pad
        lMargin = 0.12
        rMargin = 0.05
        bMargin = 0.30
        tMargin = 0.07
        padTop = ROOT.TPad('padTop','',0.,0.4,1,0.98)
        padTop.SetLeftMargin(lMargin)
        padTop.SetRightMargin(rMargin)
        padTop.SetTopMargin(tMargin)
        padTop.SetBottomMargin(0)
        padTop.SetFrameBorderMode(0);
        padTop.SetBorderMode(0);
        padTop.SetBorderSize(0);
        padTop.Draw()

        padBottom = ROOT.TPad('padBottom','',0.,0.02,1,0.4)
        padBottom.SetLeftMargin(lMargin)
        padBottom.SetRightMargin(rMargin)
        padBottom.SetTopMargin(0)
        padBottom.SetBottomMargin(bMargin)
        padBottom.SetFrameBorderMode(0);
        padBottom.SetBorderMode(0);
        padBottom.SetBorderSize(0);
        padBottom.Draw()

        chargecol = {'plus': ROOT.kOrange+1, 'minus': ROOT.kAzure+1}
        # normalize the total plus + minus = 1
        for stage in ['prefit','postfit']:
            plus  = plots_unpol['{var}_wplus_{stage}'.format(var=var,stage=stage)]
            minus = plots_unpol['{var}_wminus_{stage}'.format(var=var,stage=stage)]
            chratio = minus.Integral()/plus.Integral()
            plus.Scale(1./(1+chratio)/plus.Integral()) ; minus.Scale(chratio/(1+chratio)/minus.Integral())
        graphs = {}
        ratios = []
        for i,charge in enumerate(charges):
            padTop.cd()
            prefit  = plots_unpol['{var}_w{charge}_prefit' .format(var=var,charge=charge)]
            postfit = plots_unpol['{var}_w{charge}_postfit'.format(var=var,charge=charge)]
            ratio   = plots_unpol['{var}_w{charge}_ratio'.format(var=var,charge=charge)]
            
            ## need a graph to make the legend filled with the band color. Crazy?
            prefit_gr  = ROOT.TGraphErrors(prefit.GetNbinsX()-1)
            postfit_gr = ROOT.TGraphErrors(postfit.GetNbinsX()-1)
            prefit_gr.SetName(prefit.GetName()+'_gr')
            postfit_gr.SetName(postfit.GetName()+'_gr')
            for b in range(prefit.GetNbinsX()-1): # exclude the last point with all overflows
                prefit_gr.SetPoint(b,prefit.GetBinCenter(b+1),prefit.GetBinContent(b+1))
                prefit_gr.SetPointError(b,prefit.GetBinCenter(b+1),prefit.GetBinError(b+1))
                postfit_gr.SetPoint(b,postfit.GetBinCenter(b+1),postfit.GetBinContent(b+1))
                postfit_gr.SetPointError(b,postfit.GetBinCenter(b+1),postfit.GetBinError(b+1))

            prefit_gr.GetYaxis().SetLabelFont(42)
            prefit_gr.GetYaxis().SetLabelSize(0.05)
            prefit_gr.GetYaxis().SetLabelOffset(0.01)
            prefit_gr.GetYaxis().SetTitleSize(0.05)
            prefit_gr.GetYaxis().SetTitleOffset(1.2)
            prefit_gr.GetYaxis().SetNdivisions(505)

            ratio.GetYaxis().SetLabelFont(42)
            ratio.GetYaxis().SetLabelSize(0.05*0.55/0.4)
            ratio.GetYaxis().SetLabelOffset(0.01)
            ratio.GetYaxis().SetTitleSize(0.05*0.55/0.4)
            ratio.GetYaxis().SetTitleOffset(0.5)
            ratio.GetYaxis().SetNdivisions(505)

            ratio.GetXaxis().SetLabelFont(42)
            ratio.GetXaxis().SetLabelSize(0.05*1./0.55)
            ratio.GetXaxis().SetTitleSize(0.05*1./0.55)
            ratio.GetXaxis().SetTitleOffset(1)

            prefit_gr.GetYaxis().SetTitle('Normalized events (a.u.)')
            if 'wpt' in var: prefit_gr.GetXaxis().SetRangeUser(0,99)
            prefit_gr.SetLineColor(chargecol[charge])
            prefit_gr.SetFillColorAlpha(chargecol[charge],0.8)
            prefit_gr.SetFillStyle(3335)
            prefit_gr.SetFillColor(ROOT.kBlack)
            postfit_gr.SetLineColor(chargecol[charge]+2)
            postfit_gr.SetFillColorAlpha(chargecol[charge]+2,0.8)

            graphs['{var}_w{charge}_postfit'.format(var=var,charge=charge)] = copy.deepcopy(postfit_gr)
            graphs['{var}_w{charge}_prefit'.format(var=var,charge=charge)] = copy.deepcopy(prefit_gr)
            graphs['{var}_w{charge}_postfit'.format(var=var,charge=charge)].SetTitle('')
            graphs['{var}_w{charge}_prefit'.format(var=var,charge=charge)].SetTitle('')

            graphs['{var}_w{charge}_prefit'.format(var=var,charge=charge)].Draw('A E3' if i==0 else 'E3') ##
            graphs['{var}_w{charge}_postfit'.format(var=var,charge=charge)].Draw('E3') ##

            padBottom.cd()
            ratio.SetMarkerColor(chargecol[charge])
            ratio.SetFillColorAlpha(chargecol[charge],0.20) 
            ratio.Draw('E3' if i==0 else 'E3 same')
            ratios.append(ratio)

        padBottom.cd()
        line = ROOT.TLine()
        line.DrawLine(ratio.GetXaxis().GetBinLowEdge(1), 1, ratio.GetXaxis().GetBinLowEdge(ratio.GetNbinsX()+1), 1)
        line.SetLineStyle(3)
        line.SetLineColor(ROOT.kBlack)

        one = copy.deepcopy(ratio.Clone('one'))
        for b in range(one.GetNbinsX()):
            one.SetBinContent(b+1,1)
        #one.SetFillColorAlpha(ROOT.kGray+2,0.2)
        one.SetFillColor(ROOT.kBlack)
        one.SetFillStyle(3335)
        one.SetMarkerSize(0)
        one.SetMarkerColor(0)
        one.Draw('E3 same')
        legendRatio = doLegend(ratios+[one],['W^{{{sign}}}'.format(sign=signs[ch]) for ch in charges]+['prefit'],['f' for i in range(3)],corner='BL',legWidth=0.70,textSize=0.070,nColumns=3)
        legendRatio.Draw()

        padTop.cd()
        legendPrePost = doLegend([graphs['{var}_w{charge}_{prepost}'.format(var=var,charge=ch,prepost=step)] for ch in charges for step in ['prefit','postfit']],
                                 ['W^{{{sign}}}_{{{prepost}}}'.format(sign=signs[ch],prepost=step) for ch in charges for step in ['prefit','postfit']], 
                                 ['f' for i in range(4)],corner='TR',legWidth=0.30,nColumns=2)
        legendPrePost.Draw()

        lat = ROOT.TLatex(); lat.SetNDC()
        lat.SetTextFont(42)
        lat.SetTextSize(0.05)
        lat.DrawLatex(lMargin, 1-tMargin+0.02, '#bf{CMS}')
        lat.DrawLatex(0.72,    1-tMargin+0.02, '35.9 fb^{-1} (13 TeV)')


        for ext in ['pdf','png','C']:
            c.Print('{pdir}/{var}_postOverPrefit_unpol.{ext}'.format(pdir=options.printDir,var=var,pol=pol,ext=ext))
        c.Close()

if __name__ == '__main__':
    parser = optparse.OptionParser(usage='usage: %prog [opts] ', version='%prog 1.0')
    parser.add_option("--pdir"      , "--print-dir"  , dest="printDir"     , type="string"       , default="plots", help="print out plots in this directory");
    parser.add_option('--pf'        , '--postfix'    , dest='postfix'      , type='string'       , default=''     , help='postfix for running each module')
    parser.add_option('--wpt'       , '--wptPostfit' , dest='wptPostfit'   , action='store_true' , default=False  , help='run wpt pre/post fit comparison')
    (opts, args) = parser.parse_args()
    
    os.system('mkdir -p {od}'.format(od=opts.printDir))
    os.system('cp /afs/cern.ch/user/m/mdunser/public/index.php {od}'.format(od=opts.printDir))

    if opts.wptPostfit:
        print 'make the wpt pre/post fit ratios'
        wptPostFitRatios(opts)
        
