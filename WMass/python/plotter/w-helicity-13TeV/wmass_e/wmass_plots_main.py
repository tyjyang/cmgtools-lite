import optparse, subprocess, ROOT, datetime, math, array, copy, os, itertools
ROOT.gROOT.SetBatch(True)

POLARIZATIONS = ['left','right','long']

def doLegend(histos,lables,styles,corner="TR",textSize=0.035,legWidth=0.18,legBorder=False,nColumns=1):
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
        (x1,y1,x2,y2) = (.5, .2 + textSize*max(nentries-3,0), .5+legWidth, .35)
    elif corner == "BL":
        (x1,y1,x2,y2) = (.2, .2 + textSize*max(nentries-3,0), .2+legWidth, .35)
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
    for (plot,label,style) in zip(histos,lables,styles): leg.AddEntry(plot,label,style)
    leg.Draw()
    ## assign it to a global variable so it's not deleted
    global legend_
    legend_ = leg 
    return leg

def setRootStyle(th2=False):
    ROOT.gStyle.SetPadTopMargin(0.10)
    ROOT.gStyle.SetPadBottomMargin(0.18)
    ROOT.gStyle.SetPadLeftMargin(0.18)
    if th2:
        ROOT.gStyle.SetPadLeftMargin(0.16)                        
    ROOT.gStyle.SetTitleXOffset(1.5)
    ROOT.gStyle.SetTitleYOffset(1.1)

def wptPostFitRatios(options):
    maindir = '/eos/home-e/emanuele/www/Analysis/WMass/13TeV/plots/gen/wpt/'

    colors = {('plus','left'):ROOT.kBlue-1,('plus','right'):ROOT.kGreen+1,('plus','long'):ROOT.kGray+1,
              ('minus','left'):ROOT.kViolet-1,('minus','right'):ROOT.kAzure+1,('minus','long'):ROOT.kYellow+1,}
    
    setRootStyle()
    ROOT.gStyle.SetOptStat(0)
    ## should propagate better the error instead of taking a flat envelope
    kappa = 0.03

    def makeRatio(hprefit,hpostfit,color):
        hpostfit.Scale(hprefit.Integral()/hpostfit.Integral())
        hpostfit.Divide(hprefit)
        for b in xrange(hpostfit.GetNbinsX()):
            hpostfit.SetBinError(b+1,math.hypot(kappa,hpostfit.GetBinError(b+1)))
        hpostfit.SetFillColorAlpha(color,0.20)
        hpostfit.SetLineColor(color)
        hpostfit.SetMarkerColor(color)
        hpostfit.GetYaxis().SetTitle('postfit / prefit')
        hpostfit.GetXaxis().SetTitleOffset(1.5)
        hpostfit.SetTitle('')
        hpostfit.GetYaxis().SetRangeUser(0.8,1.2)
        if 'wpt' in hpostfit.GetName():
            hpostfit.GetXaxis().SetRangeUser(1,100)
            hpostfit.GetXaxis().SetTitle('p_{T}^{W} (GeV)')
        elif 'wy' in hpostfit.GetName():
            hpostfit.GetXaxis().SetRangeUser(-5,5)
            hpostfit.GetXaxis().SetTitle('Y_{W}')

    for charge in ['plus','minus']:
        sign='+' if charge=='plus' else '-'
        for var in ['wpt','wy']:
            c = ROOT.TCanvas('c','',600,600)
            if var=='wpt': c.SetLogx(1)
            else: c.SetLogx(0)
            rfiles = {}
            rfiles['prefit']  = ROOT.TFile.Open('{d}/wgen_nosel/w{charge}_{var}.root'.format(d=maindir,charge=charge,var=var))
            rfiles['postfit'] = ROOT.TFile.Open('{d}/wgen_nosel_qcdpostfit/w{charge}_{var}.root'.format(d=maindir,charge=charge,var=var))
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
                makeRatio(prefit,postfit,colors[(charge,pol)])
                postfit.Draw('E3' if i==0 else 'E3 same')
                histos.append(postfit)
                labels.append('W{sign} {pol}'.format(sign=sign,pol=pol))
                styles.append('pf')
            ## draw polarized
            legPol = doLegend(histos,labels,styles,corner="BL")
            legPol.Draw()
            for ext in ['pdf','png']:
                c.Print('{pdir}/{var}_postOverPrefit_{charge}.{ext}'.format(pdir=options.printDir,var=var,charge=charge,pol=pol,ext=ext))
            legPol.Clear()
            ## draw unpolarized
            makeRatio(prefit_unpol,postfit_unpol,ROOT.kRed+1)
            postfit_unpol.Draw('E3')
            legUnpol = doLegend([postfit_unpol],['W{sign}'.format(sign=sign)],['pf'],corner="BL")
            legUnpol.Draw()
            for ext in ['pdf','png']:
                c.Print('{pdir}/{var}_postOverPrefit_{charge}_unpol.{ext}'.format(pdir=options.printDir,var=var,charge=charge,pol=pol,ext=ext))        
            legUnpol.Clear()
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
        
