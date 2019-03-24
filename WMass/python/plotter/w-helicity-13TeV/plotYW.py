#!/bin/env python
# USAGE: python plotYW.py --type toys --infile toys_wplus.root -y cards_el/binningYW.txt -C plus -o plots --xsecfiles Wel_plus_shapes_xsec.root [--normxsec]
#
# When run on the W+ W- combined fit it also makes charge asymmetry plot:
# python plotYW.py --type toys --infile toys_wboth.root -y cards_el/binningYW.txt -C plus,minus -o plots --xsecfiles Wel_plus_shapes_xsec.root,Wel_minus_shapes_xsec.root

import ROOT, datetime, array, os, math, re
ROOT.gROOT.SetBatch(True)
ROOT.PyConfig.IgnoreCommandLineOptions = True
ROOT.gSystem.Load("libHiggsAnalysisCombinedLimit")

import utilities
utilities = utilities.util()

REFMC = 'MC@NLO'

class valueClass:
    def __init__(self, name):
        self.name = name

        self.pol = 'left' if 'left' in self.name else 'right' if 'right' in self.name else 'long'
        if 'unpolarized' in self.name: self.pol = 'unpolarized'
        self.isleft  = self.pol == 'left'
        self.isright = self.pol == 'right'
        self.islong  = self.pol == 'long'
        self.isunpolarized = self.pol == 'unpolarized'

        self.charge = 'plus' if 'plus' in name else 'minus'
        self.ch     = '+' if 'plus' in name else '-'
        if 'asymmetry' in name:
            self.charge = self.ch = ''

        self.color  = ROOT.kBlue-9 if self.isleft else ROOT.kOrange-4 if self.isright else ROOT.kGray+1
        self.colorf = ROOT.kBlue-3 if self.isleft else ROOT.kOrange-3 if self.isright else ROOT.kGray+3
        if self.isunpolarized: 
            self.color = ROOT.kSpring-6
            self.colorf = ROOT.kSpring-7

        ## here all the arrays that will contain the values and errors etc.
        self.val       = array.array('f', []); self.ehi       = array.array('f', []); self.elo       = array.array('f', []);
        self.val_fit   = array.array('f', []); self.ehi_fit   = array.array('f', []); self.elo_fit   = array.array('f', []);
        self.relv      = array.array('f', []); self.relhi     = array.array('f', []); self.rello     = array.array('f', []);
        self.relv_fit  = array.array('f', []); self.relhi_fit = array.array('f', []); self.rello_fit = array.array('f', []);
        self.rap       = array.array('f', []); self.rlo       = array.array('f', []); self.rhi       = array.array('f', []);

    def makeGraphs(self):
        if len(self.val):
            self.graph = ROOT.TGraphAsymmErrors(len(self.val), self.rap, self.val, self.rlo, self.rhi, self.elo, self.ehi)
            self.graph.SetName('graph'+self.pol)
        if len(self.relv): 
            self.graph_rel= ROOT.TGraphAsymmErrors(len(self.relv), self.rap, self.relv, self.rlo, self.rhi, self.rello, self.relhi)
            self.graph_rel.SetName('graph'+self.pol+'_rel')
        if len(self.val_fit):
            self.graph_fit = ROOT.TGraphAsymmErrors(len(self.val_fit), self.rap, self.val_fit, self.rlo, self.rhi, self.elo_fit, self.ehi_fit)
            self.graph_fit.SetName('graph'+self.pol+'_fit')
        zeros = array.array('f',[0 for i in xrange(len(self.rlo))])
        if len(self.relv_fit):
            self.graph_fit_rel = ROOT.TGraphAsymmErrors(len(self.relv_fit), self.rap, self.relv_fit, zeros, zeros, self.rello_fit, self.relhi_fit)
            self.graph_fit_rel.SetName('graph'+self.pol+'_fit_rel')
        
        self.graphStyle()
        if len(self.relv) and len(self.relv_fit): self.makeMultiGraphRel()

    def makeMultiGraphRel(self):
        self.mg = ROOT.TMultiGraph()
        self.mg.SetTitle() ## no title 'W^{{{ch}}}: {p}'.format(ch=self.ch,p=self.pol))
        self.shiftPoints(self.graph_fit_rel)
        self.mg.Add(self.graph_rel,'P2')
        self.mg.Add(self.graph_fit_rel)

    def graphStyle(self):
        #fillstyles = {'left': 3017, 'right': 3018, 'long': 3016, 'unpolarized': 3020}
        fillstyles = {'left': 3002, 'right': 3144, 'long': 3016, 'unpolarized': 3020}
        if hasattr(self,'graph'):
            self.graph.SetLineColor(self.color)
            self.graph.SetFillColor(self.color)
            self.graph.SetFillStyle(fillstyles[self.pol])
        if hasattr(self,'graph_fit'):
            self.graph_fit.SetLineWidth(2)
            self.graph_fit.SetMarkerSize(1.0)
            self.graph_fit.SetMarkerStyle(ROOT.kFullCircle)
            self.graph_fit.SetMarkerColor(self.color)
            self.graph_fit.SetLineColor(self.color)
        if hasattr(self,'graph_rel'):
            self.graph_rel.SetLineWidth(5)
            self.graph_rel.SetLineColor(self.color)
            self.graph_rel.SetFillColor(self.color)
            self.graph_rel.SetFillStyle(fillstyles[self.pol])
        if hasattr(self,'graph_fit_rel'):
            self.graph_fit_rel.SetLineWidth(2)
            self.graph_fit_rel.SetMarkerSize(1.0)
            self.graph_fit_rel.SetMarkerStyle(ROOT.kFullCircle)
            self.graph_fit_rel.SetLineColor(self.color)
            self.graph_fit_rel.SetMarkerColor(self.color)

    def shiftPoints(self, graph):
        shifts = {'left': -0.01, 'right': 0.01, 'long': 0.0, 'unpolarized': 0.0}
        for p in xrange(graph.GetN()):
            x = ROOT.Double(0); y = ROOT.Double(0)
            graph.GetPoint(p,x,y)
            graph.SetPoint(p,x+shifts[self.pol],y)
            graph.SetPointEXhigh(p,0); graph.SetPointEXlow(p,0)

def plotValues(values,charge,channel,options):
        c2 = ROOT.TCanvas('foo','', 800, 800)
        c2.GetPad(0).SetTopMargin(0.09)
        c2.GetPad(0).SetBottomMargin(0.35)
        c2.GetPad(0).SetLeftMargin(0.17)
        c2.GetPad(0).SetRightMargin(0.04)
        c2.GetPad(0).SetTickx(1)
        c2.GetPad(0).SetTicky(1)

        ch = '#plus' if charge == 'plus' else '#minus'
        if charge == 'asymmetry': ch = ''
        date = datetime.date.today().isoformat()
        normstr = 'norm' if (options.normxsec and charge!='asymmetry') else ''

        lat = ROOT.TLatex()
        lat.SetNDC(); lat.SetTextFont(42)
        ## the four graphs exist now. now starting to draw them
        ## ===========================================================
        if sum(hasattr(values[pol],'graph') and hasattr(values[pol],'graph_fit') for pol in ['left','right','long'])==3:
            leg = ROOT.TLegend(0.40, 0.80, 0.90, 0.90)
            leg.SetFillStyle(0)
            leg.SetBorderSize(0)
            leg.AddEntry(values['left'] .graph     , 'W_{{L}} ({mc})'.format(mc=REFMC) , 'f')
            leg.AddEntry(values['left'] .graph_fit , 'W_{L} (fit)', 'pl')
            leg.AddEntry(values['right'].graph     , 'W_{{R}} ({mc})'.format(mc=REFMC) , 'f')
            leg.AddEntry(values['right'].graph_fit , 'W_{R} (fit)', 'pl') 
            leg.SetNColumns(2)
            if not options.nolong:
                leg.AddEntry(values['long'] .graph     , 'W_{{0}} ({mc})'.format(mc=REFMC) , 'f')
                leg.AddEntry(values['long'] .graph_fit , 'W_{0} (fit)', 'pl')
                leg.SetNColumns(3)

            values['left'].graph.SetTitle('W {ch}: Y_{{W}}'.format(ch=ch))
                
            mg = ROOT.TMultiGraph()
            mg.Add(values['left'] .graph,'P2')
            mg.Add(values['right'].graph,'P2')
            if not options.nolong: mg.Add(values['long'] .graph,'P2')
            mg.Add(values['left'] .graph_fit)
            mg.Add(values['right'].graph_fit)
            if not options.nolong: mg.Add(values['long'] .graph_fit)
     
            mg.Draw('Pa')
            mg.GetXaxis().SetRangeUser(0., options.maxRapidity) # max would be 6.
            mg.GetXaxis().SetTitle('|Y_{W}|')
            mg.GetXaxis().SetLabelSize(0)
            if charge=='asymmetry':
                mg.GetYaxis().SetTitle('Charge asymmetry')
                mg.GetYaxis().SetRangeUser(-0.1,0.4)
            else:
                if options.normxsec: 
                    mg.GetYaxis().SetTitle('#frac{d#sigma}{#sigma_{tot}^{fit}} / d|Y_{W}|')
                    mg.GetYaxis().SetRangeUser(-0.05,0.8 if options.maxRapidity > 2.9 else 0.4)
                else: 
                    mg.GetYaxis().SetTitle('d#sigma (pb) / d|Y_{W}|')
                    mg.GetYaxis().SetRangeUser(-200,3500)
            mg.GetYaxis().SetTitleSize(0.04)
            mg.GetYaxis().SetLabelSize(0.04)
            mg.GetYaxis().SetTitleOffset(2.0)
     
            leg.Draw('same')
     
        ## now make the relative error plot:
        ## ======================================
        if sum(hasattr(values[pol],'mg') for pol in ['left','right','long'])==3:

            pad2 = ROOT.TPad("pad2","pad2",0,0.,1,0.9)
            pad2.SetTopMargin(0.65)
            pad2.SetRightMargin(0.04)
            pad2.SetLeftMargin(0.17)
            pad2.SetFillColor(0)
            pad2.SetGridy(0)
            pad2.SetFillStyle(0)
            pad2.SetTicky(1)

            pad2.Draw()
            pad2.cd()
     
            line = ROOT.TF1("horiz_line","0" if charge=='asymmetry' else '1',0.0,3.0);
            line.SetLineColor(ROOT.kBlack);
            line.SetLineWidth(2);

            if charge=='asymmetry':
                yaxtitle = 'A_{Theory}-A_{Data}'
                yaxrange = (-0.1, 0.1)
            else:
                yaxtitle = '#sigma_{Theory}/#sigma_{Data}'
                yaxrange = (0.70, 1.30)
     

            helToPlot = ['left','right']
            if not options.nolong: helToPlot.append('long')
            for  ih,hel in enumerate(helToPlot):
     
                values[hel].mg.Draw('Pa' if ih==0 else 'P')
                if ih==0:
                    ## x axis fiddling
                    values[hel].mg.GetXaxis().SetRangeUser(0., options.maxRapidity)
                    values[hel].mg.GetXaxis().SetLabelSize(0.04)
                    ## y axis fiddling
                    values[hel].mg.GetYaxis().SetTitleOffset(1.8)
                    values[hel].mg.GetYaxis().SetTitleSize(0.04)
                    values[hel].mg.GetYaxis().SetLabelSize(0.04)
                    values[hel].mg.GetYaxis().SetTitle(yaxtitle)
                    values[hel].mg.GetYaxis().SetRangeUser(yaxrange[0],yaxrange[1])
                    values[hel].mg.GetYaxis().SetNdivisions(510)
                    values[hel].mg.GetYaxis().CenterTitle()
            line.Draw("Lsame");
            c2.cd()
            lat.DrawLatex(0.16, 0.92, '#bf{CMS} #it{Preliminary}')
            lat.DrawLatex(0.62, 0.92, '35.9 fb^{-1} (13 TeV)')
            lat.DrawLatex(0.20, 0.80,  'W^{{{ch}}} #rightarrow {lep}^{{{ch}}}{nu}'.format(ch=ch,lep="#mu" if channel == "mu" else "e",nu="#bar{#nu}" if charge=='minus' else "#nu"))
            lat.DrawLatex(0.90, 0.03, '|Y_{W}|')
        for ext in ['png', 'pdf']:
            c2.SaveAs('{od}/genAbsY{norm}_pdfs_{ch}{suffix}_{t}.{ext}'.format(od=options.outdir, norm=normstr, ch=charge, suffix=options.suffix, ext=ext,t=options.type))


def plotUnpolarizedValues(values,charge,channel,options):
        c2 = ROOT.TCanvas('foo','', 800, 800)
        c2.GetPad(0).SetTopMargin(0.09)
        c2.GetPad(0).SetBottomMargin(0.35)
        c2.GetPad(0).SetLeftMargin(0.15)
        c2.GetPad(0).SetRightMargin(0.04)
        c2.GetPad(0).SetTickx(1)
        c2.GetPad(0).SetTicky(1)

        ch = '#plus' if charge == 'plus' else '#minus'
        if charge == 'asymmetry': ch = ''
        date = datetime.date.today().isoformat()
        if 'values_a0' in values.name: normstr = 'A0'
        elif 'values_a4' in values.name: normstr = 'A4'
        elif 'values_sumxsec' in values.name: normstr = 'xsec'
        else: normstr = 'norm' if (options.normxsec and charge!='asymmetry') else ''

        lat = ROOT.TLatex()
        lat.SetNDC(); lat.SetTextFont(42)
        legx1, legx2, legy1, legy2 = 0.2, 0.5, 0.7, 0.85
        ## the graphs exist now. now starting to draw them
        ## ===========================================================
        if hasattr(values,'graph') and hasattr(values,'graph_fit'):
                
            mg = ROOT.TMultiGraph()
            mg.Add(values.graph,'P2')
            mg.Add(values.graph_fit)
            mg.Draw('Pa')
            mg.GetXaxis().SetRangeUser(0., options.maxRapidity) # max would be 6.
            mg.GetXaxis().SetTitle('|Y_{W}|')
            mg.GetXaxis().SetTitleOffset(5.5)
            mg.GetXaxis().SetLabelSize(0)
            if charge=='asymmetry':
                 mg.GetYaxis().SetTitle('Charge asymmetry')
                 mg.GetYaxis().SetRangeUser(-0.1,0.4)
            else:
                if normstr=='xsec':
                    if options.normxsec: 
                        mg.GetYaxis().SetTitle('#frac{d#sigma}{#sigma_{tot}^{fit}} / d|Y_{W}|')
                        mg.GetYaxis().SetRangeUser(-0.05,0.8 if options.maxRapidity > 2.7 else 0.4)
                    else: 
                        mg.GetYaxis().SetTitle('d#sigma (pb) / d|Y_{W}|')
                        mg.GetYaxis().SetRangeUser(1500,4500)
                else:
                    mg.GetYaxis().SetRangeUser(-0.05 if normstr=='A0' else -1,0.4 if normstr=='A0' else 2)
                    mg.GetYaxis().SetTitle('|A_{0}|' if normstr=='A0' else '|A_{4}|')
            mg.GetYaxis().SetTitleSize(0.06)
            mg.GetYaxis().SetLabelSize(0.04)
            mg.GetYaxis().SetTitleOffset(1.)
     
            leg = ROOT.TLegend(legx1, legy1, legx2, legy2)
            leg.SetFillStyle(0)
            leg.SetBorderSize(0)
            leg.AddEntry(values.graph_fit , 'data', 'pl')
            leg.AddEntry(values.graph     , REFMC, 'f')

            leg.Draw('same')
            lat.DrawLatex(0.16, 0.92, '#bf{CMS} #it{Preliminary}')
            lat.DrawLatex(0.62, 0.92, '35.9 fb^{-1} (13 TeV)')
            lat.DrawLatex(0.20, 0.50,  'W^{{{ch}}} #rightarrow {lep}^{{{ch}}}{nu}'.format(ch=ch,lep="#mu" if channel == "mu" else "e",nu="#bar{#nu}" if charge=='minus' else "#nu"))
     

        ## now make the relative error plot:
        ## ======================================
        if hasattr(values,'mg'):

            pad2 = ROOT.TPad("pad2","pad2",0,0.,1,0.9)
            pad2.SetTopMargin(0.65)
            pad2.SetRightMargin(0.04)
            pad2.SetLeftMargin(0.15)
            pad2.SetFillColor(0)
            pad2.SetGridy(0)
            pad2.SetFillStyle(0)
            pad2.SetTicky(1)

            pad2.Draw()
            pad2.cd()

            line = ROOT.TF1("horiz_line","0" if charge=='asymmetry' else '1',0.0,3.0);
            line.SetLineColor(ROOT.kBlack);
            line.SetLineWidth(2);
            yaxrange = (0,0)
            if charge=='asymmetry':
                yaxtitle = 'A_{Theory}-A_{Data}'
                yaxrange = (-0.2, 0.2)
            else:
                yaxtitle = '#sigma_{Theory}/#sigma_{Data}'
                yaxrange = (0.70, 1.30) if normstr=='xsec' else (-1.5, 3.5)

            values.mg.Draw('Pa')
            ## x axis fiddling
            values.mg.GetXaxis().SetTitle('|Y_{W}|')
            values.mg.GetXaxis().SetTitleOffset(1.)
            values.mg.GetXaxis().SetRangeUser(0., options.maxRapidity)
            values.mg.GetXaxis().SetTitleSize(0.1)
            values.mg.GetXaxis().SetLabelSize(0.04)
            ## y axis fiddling
            values.mg.GetYaxis().SetTitleOffset(1.2)
            values.mg.GetYaxis().SetTitleSize(0.06)
            values.mg.GetYaxis().SetLabelSize(0.04)
            values.mg.GetYaxis().SetTitle(yaxtitle)
            values.mg.GetYaxis().SetRangeUser(yaxrange[0],yaxrange[1])
            values.mg.GetYaxis().SetNdivisions(5)
            values.mg.GetYaxis().CenterTitle()

            line.Draw("Lsame");

        for ext in ['png', 'pdf']:
            c2.SaveAs('{od}/genAbsYUnpolarized{norm}_pdfs_{ch}{suffix}_{t}.{ext}'.format(od=options.outdir, norm=normstr, ch=charge, suffix=options.suffix, ext=ext,t=options.type))


NPDFs = 60
LUMINOSITY = 36000

if __name__ == "__main__":

    from optparse import OptionParser
    parser = OptionParser(usage='%prog ybinfile workspace.root toys.root [options] ')
    parser.add_option('-i', '--infile'      , dest='infile'   , default=''            , type='string', help='workspace converted from datacard')
    parser.add_option('-y', '--ybinfile'    , dest='ybinfile' , default=''            , type='string', help='file with the yw binning')

    parser.add_option('-t', '--type'        , dest='type'     , default='toys'        , type='string', help='run the plot from which postfit? toys/scans/hessian')
    parser.add_option(      '--toyfile'     , dest='toyfile'  , default=''            , type='string', help='file that has the toys')
    parser.add_option(      '--scandir'     , dest='scandir'  , default=''            , type='string', help='directory with all the scans')
    parser.add_option(      '--hessfile'    , dest='hessfile' , default=''            , type='string', help='file that contains the hessian errors in a dictionary')
    parser.add_option(      '--xsecfiles'    , dest='xsecfiles' , default=None          , type='string', help='files that contains the expected x sections with variations (one per charge,comma separated in the same order of the charges) ')
    parser.add_option('-C', '--charge'      , dest='charge'   , default='plus,minus'  , type='string', help='process given charge. default is both')
    parser.add_option('-c', '--channel'     , dest='channel'   , default='mu'  , type='string', help='Channel (mu|el)')
    parser.add_option('-o', '--outdir'      , dest='outdir'   , default='.'           , type='string', help='outdput directory to save the plots')
    parser.add_option(      '--suffix'      , dest='suffix'   , default=''            , type='string', help='suffix for the correlation matrix')
    parser.add_option('-n', '--normxsec'    , dest='normxsec' , default=False         , action='store_true',   help='if given, plot the differential xsecs normalized to the total xsec')
    parser.add_option(      '--nolong'      , dest='nolong'   , default=False         , action='store_true',   help='if given, do not plot longitudinal component')
    parser.add_option(      '--max-rap'     , dest='maxRapidity', default='2.75'       , type='float', help='Max value for rapidity range')
    (options, args) = parser.parse_args()


    if not os.path.isdir(options.outdir):
        os.system('mkdir {od}'.format(od=options.outdir))
        if os.path.exists("/afs/cern.ch"): os.system("cp /afs/cern.ch/user/g/gpetrucc/php/index.php {od}".format(od=options.outdir))

    if options.ybinfile:
        ybinfile = options.ybinfile
    else:
        ybinfile = os.path.dirname(os.path.abspath(options.infile))+'/binningYW.txt'


    ## get the central values and uncertainties depending on the type given:

    ## if --type=toys   , we expect a toyfile
    ## if --type=scans  , we expect a scan directory
    ## if --type=hessian, we expect a hessian file

    if   options.type == 'toys':
        valuesAndErrors = utilities.getFromToys(options.infile)
    elif options.type == 'scans':
        valuesAndErrors = utilities.getFromScans(options.infile)
    elif options.type == 'hessian':
        valuesAndErrors = utilities.getFromHessian(options.infile)
    else:
        print 'ERROR: none of your types is supported. specify either "toys", "scans", or "hessian"'
        sys.exit()

    #channel = 'mu' if any(re.match('.*_mu_Ybin_.*', param) for param in valuesAndErrors.keys()) else 'el'
    channel = options.channel
    print "From the list of parameters it seems that you are plotting results for channel ",channel

    ybinfile = open(ybinfile, 'r')
    ybins = eval(ybinfile.read())
    ybinfile.close()

    ## calculate the bin widths for the rapidity bins
    ybinwidths = {}
    for k,v in ybins.items():
        tmplist = list(abs(i - v[v.index(i)+1]) for i in v[:-1])
        ybinwidths[k] = [float('{n:.2f}'.format(n=i)) for i in tmplist]


    charges = options.charge.split(',')
    xsecfiles = options.xsecfiles.split(',')
    xsec_nominal_allCharges = {}; xsec_systematics_allCharges = {}
    polarizations = ['left','right', 'long']
    for ic,charge in enumerate(charges):

        sign = 1. if charge=='plus' else -1.

        ## this gets the pdf central variation binned in the correct format
        xsec_nominal = utilities.getXSecFromShapes(ybins,charge,xsecfiles[ic],channel, 0)
        xsec_nominal_allCharges[charge] = xsec_nominal

        value_syst = {}
        for pol in polarizations:
            histos = []
            values = []
            for ip in xrange(1,NPDFs+1):
                # print "Loading polarization %s, histograms for pdf %d" % (pol,ip)
                ## this gets the pdf variations after correctly rebinning the YW
                xsec_pdf = utilities.getXSecFromShapes(ybins,charge,xsecfiles[ic],channel,ip)
                values.append(xsec_pdf[pol])
            value_syst[pol] = values

        xsec_systematics = {}
        for pol in polarizations:
            #print "===> Running pol = ",pol
            xsec_systs=[]
            for iy,y in enumerate(ybinwidths['{ch}_{pol}'.format(ch=charge,pol=pol if not pol=='long' else 'right')]):
                xsec_nom = xsec_nominal[pol][iy]
                #print "\tBin iy={iy},y={y}. Nom = {nom} ".format(iy=iy,y=y,nom=nom)
                totUp=0; xsec_totUp=0
                for ip,pdf in enumerate(value_syst[pol]):
                    xsec_pdf = value_syst[pol][ip]
                    #print "\tip = {ip}  pdf = {pdf}".format(ip=ip,pdf=pdf[iy])
                    # debug
                    xsec_relsyst = abs(xsec_nom-xsec_pdf[iy])/xsec_nom if xsec_nom else 0.0
                    if xsec_relsyst>0.20:
                        print "SOMETHING WENT WRONG WITH THIS PDF: %d HAS RELATIVE SYST = %f. SKIPPING !" % (ip,relsyst)
                    else:
                        xsec_totUp += math.pow(xsec_relsyst*xsec_nom,2)
                xsec_totUp = math.sqrt(xsec_totUp)
                # print "Rel systematic for Y bin %d = +/-%.3f" % (iy,totUp/nom)
                # print "\tRel systematic on xsec for Y bin %d = +/-%.3f" % (iy,xsec_totUp/xsec_nom if xsec_nom else 0.)
                xsec_systs.append(xsec_totUp)
            xsec_systematics[pol]=xsec_systs
        xsec_systematics_allCharges[charge] = xsec_systematics

        angcoeff_nominal = {'sumxsec': [], 'a0': [], 'a4': []}
        angcoeff_systematics = {'sumxsec': [], 'a0': [], 'a4': []}
        for  iy,y in enumerate(ybinwidths['{ch}_{pol}'.format(ch=charge,pol=pol if not pol=='long' else 'right')]):
            xsec_unpolarized_nominal_iy = sum([xsec_nominal[pol][iy] for pol in polarizations])
            angcoeff_nominal['sumxsec'].append(xsec_unpolarized_nominal_iy)

            coeffs_val = utilities.getCoeffs(xsec_nominal['left'][iy],     xsec_nominal['right'][iy],     xsec_nominal['long'][iy],
                                             xsec_systematics['left'][iy], xsec_systematics['right'][iy], xsec_systematics['long'][iy])

            angcoeff_nominal['a0'].append(coeffs_val['a0'][0])
            angcoeff_nominal['a4'].append(sign*coeffs_val['a4'][0])

            xsec_unpolarized_iy = sum([xsec_systematics[pol][iy] for pol in polarizations])
            angcoeff_systematics['sumxsec'].append(xsec_unpolarized_iy)
            angcoeff_systematics['a0'].append(coeffs_val['a0'][1])
            angcoeff_systematics['a4'].append(coeffs_val['a4'][1])

        nOuterBinsToExclude = 1  ### EDM hardcoded: out of acceptance Y bins

        allValues = {}
        for pol in polarizations:
            cp = '{ch}_{pol}'.format(ch=charge,pol=pol)
            MAXYFORNORM = ybins[cp][-nOuterBinsToExclude-1] # exclude the outermost 2 bins which has huge error due to acceptance
            normsigmaIn = sum([xsec_nominal[allpol][iy] for allpol in polarizations for iy,y in enumerate(ybins[cp][:-1]) if abs(y)<MAXYFORNORM])
            normsigmaOut = sum([xsec_nominal[allpol][iy] for allpol in polarizations for iy,y in enumerate(ybins[cp][:-1]) if abs(y)>=MAXYFORNORM])
            normsigmaInFit = sum([valuesAndErrors['W{charge}_{pol}_W{charge}_{pol}_{ch}_Ybin_{iy}_pmaskedexp'.format(charge=charge,pol=allpol,ch=channel,iy=iy)][0]/LUMINOSITY for allpol in polarizations for iy,y in enumerate(ybins[cp][:-1]) if abs(y)<MAXYFORNORM])
            normsigmaOutFit = sum([valuesAndErrors['W{charge}_{pol}_W{charge}_{pol}_{ch}_Ybin_{iy}_pmaskedexp'.format(charge=charge,pol=allpol,ch=channel,iy=iy)][0]/LUMINOSITY for allpol in polarizations for iy,y in enumerate(ybins[cp][:-1]) if abs(y)>=MAXYFORNORM])

            print "total expected (fit) xsec up to |Y|<{maxy} = {sigma:.3f} ({fit:.3f}) pb".format(maxy=MAXYFORNORM,sigma=normsigmaIn,fit=normsigmaInFit)
            print "total expected (fit) xsec beyond |Y|>{maxy} = {sigma:.3f} ({fit:.3f}) pb".format(maxy=MAXYFORNORM,sigma=normsigmaOut,fit=normsigmaOutFit)

            tmp_val = valueClass('values_'+charge+'_'+pol)

            for iy,y in enumerate(ybinwidths['{ch}_{pol}'.format(ch=charge,pol=pol)]):
                normsigma = normsigmaInFit if abs(ybins[cp][iy])<MAXYFORNORM else normsigmaOutFit
                parname = 'W{charge}_{pol}_W{charge}_{pol}_{ch}_Ybin_{iy}'.format(charge=charge,pol=pol,ch=channel,iy=iy)

                scale = 1.
                if options.normxsec:
                    if   options.type == 'toys': 
                        xsec_fit = utilities.getNormalizedXsecFromToys(ybins,charge,pol,channel,iy,options.infile,MAXYFORNORM)
                    elif options.type == 'hessian':
                        xsec_fit = valuesAndErrors[parname+'_pmaskedexpnorm']
                    else:
                        print "--normxsec not implemented yet for scans."
                        sys.exit()
                else:
                    xsec_fit = valuesAndErrors[parname+'_pmaskedexp']
                    scale = LUMINOSITY

                if options.normxsec:
                    rfit     = xsec_nominal[pol][iy]/normsigma/xsec_fit[0]
                else:
                    rfit     = xsec_nominal[pol][iy]/xsec_fit[0]*scale
                rfit_err = rfit*abs(xsec_fit[0]-xsec_fit[1])/xsec_fit[0] 

                tmp_val.val.append(xsec_nominal[pol][iy]/ybinwidths[cp][iy])
                tmp_val.ehi.append(xsec_systematics[pol][iy]/ybinwidths[cp][iy])
                tmp_val.elo.append(xsec_systematics[pol][iy]/ybinwidths[cp][iy]) # symmetric for the expected
                if options.normxsec:
                    tmp_val.val[-1] = tmp_val.val[-1]/normsigma
                    tmp_val.ehi[-1] = tmp_val.ehi[-1]/normsigma
                    tmp_val.elo[-1] = tmp_val.elo[-1]/normsigma

                tmp_val.relv. append(rfit);
                tmp_val.rello.append(xsec_systematics[pol][iy]/xsec_nominal[pol][iy])
                tmp_val.relhi.append(xsec_systematics[pol][iy]/xsec_nominal[pol][iy]) # symmetric for the expected
                
                tmp_val.val_fit.append(xsec_fit[0]/ybinwidths[cp][iy]/scale)
                tmp_val.elo_fit.append(abs(xsec_fit[0]-xsec_fit[1])/ybinwidths[cp][iy]/scale)
                tmp_val.ehi_fit.append(abs(xsec_fit[0]-xsec_fit[2])/ybinwidths[cp][iy]/scale)

                units = '' if options.normxsec else '(pb)'
                print "par = {parname}, expected sigma = {sigma:.3f} {units}   fitted = {val:.3f} + {ehi:.3f} - {elo:.3f} {units}".format(parname=parname,
                                                                                                                                          sigma=tmp_val.val[-1],units=units,
                                                                                                                                          val=tmp_val.val_fit[-1],ehi=tmp_val.ehi_fit[-1],elo=tmp_val.elo_fit[-1])

                tmp_val.relv_fit .append(1.)
                tmp_val.rello_fit.append(rfit_err)
                tmp_val.relhi_fit.append(rfit_err)

                tmp_val.rap.append((ybins[cp][iy]+ybins[cp][iy+1])/2.)
                tmp_val.rlo.append(abs(ybins[cp][iy]-tmp_val.rap[-1]))
                tmp_val.rhi.append(abs(ybins[cp][iy]-tmp_val.rap[-1]))

            tmp_val.makeGraphs()

            allValues[pol] = tmp_val

        plotValues(allValues,charge,channel,options)

        if not options.normxsec: # this is only implemented for absolute xsecs
            # now do the unpolarized ones
            cp = 'plus_left' # this works if the binning for all the pol is the same
            xsec_params = ['sumxsec','a0','a4']
            for xs in xsec_params:
                tmp_val = valueClass('values_{xs}_{charge}_unpolarized'.format(xs=xs,charge=charge))
                for iy,y in enumerate(ybinwidths['{ch}_{pol}'.format(ch=charge,pol=pol)]):
                    parname = 'W{charge}_{ch}_Ybin_{iy}_{xs}'.format(charge=charge,ch=channel,iy=iy,xs=xs)
                    if xs=='sumxsec':
                        ybinwidth_scale = ybinwidths[cp][iy]
                        scale = LUMINOSITY
                    else:
                        ybinwidth_scale = 1.
                        scale = 1.

                    tmp_val.val.append(abs(angcoeff_nominal[xs][iy]/ybinwidth_scale))
                    experr = angcoeff_systematics[xs][iy]/ybinwidth_scale
                    tmp_val.ehi.append(experr)
                    tmp_val.elo.append(experr) # symmetric for the expected
         
                    xsec_fit = valuesAndErrors[parname]
                    scale = LUMINOSITY if xs=='sumxsec' else 1.
         
                    tmp_val.val_fit.append(xsec_fit[0]/ybinwidth_scale/scale)
                    tmp_val.elo_fit.append(abs(xsec_fit[0]-xsec_fit[1])/ybinwidth_scale/scale)
                    tmp_val.ehi_fit.append(abs(xsec_fit[0]-xsec_fit[2])/ybinwidth_scale/scale)

                    tmp_val.relv. append(tmp_val.val[-1]/tmp_val.val_fit[-1])
                    experrrel = angcoeff_systematics[xs][iy]/angcoeff_nominal[xs][iy]
                    tmp_val.rello.append(experrrel)
                    tmp_val.relhi.append(experrrel) # symmetric for the expected
                    
                    units = '(pb)' if xs=='sumxsec' else ''
                    print "par = {parname}, expected value = {sigma:.3f} {units}   fitted = {val:.3f} + {ehi:.3f} - {elo:.3f} {units}".format(parname=parname, sigma=tmp_val.val[-1],units=units,
                                                                                                                                              val=tmp_val.val_fit[-1],ehi=tmp_val.ehi_fit[-1],elo=tmp_val.elo_fit[-1])
                    tmp_val.relv_fit .append(1.)
                    tmp_val.rello_fit.append(tmp_val.elo_fit[-1]/tmp_val.val_fit[-1])
                    tmp_val.relhi_fit.append(tmp_val.ehi_fit[-1]/tmp_val.val_fit[-1])
         
                    tmp_val.rap.append((ybins[cp][iy]+ybins[cp][iy+1])/2.)
                    tmp_val.rlo.append(abs(ybins[cp][iy]-tmp_val.rap[-1]))
                    tmp_val.rhi.append(abs(ybins[cp][iy]-tmp_val.rap[-1]))
         
                tmp_val.makeGraphs()
                plotUnpolarizedValues(tmp_val,charge,channel,options)


    if len(charges)>1:
        print "Making charge asymmetry plots now..."
        asymmetryValues = {}
        
        for pol in polarizations:
            cp = 'plus_'+pol
            tmp_val = valueClass('asymmetry_'+pol)
            for iy,y in enumerate(ybinwidths[cp]):
                chasy_val = utilities.getChargeAsy(xsec_nominal_allCharges['plus'][pol][iy],     xsec_nominal_allCharges['minus'][pol][iy],
                                                   xsec_systematics_allCharges['plus'][pol][iy], xsec_systematics_allCharges['minus'][pol][iy])
                tmp_val.val .append(chasy_val['asy'][0])
                tmp_val.ehi.append(chasy_val['asy'][1])
                tmp_val.elo.append(chasy_val['asy'][1])

                if options.type == 'toys':
                    asy_fit = utilities.getAsymmetryFromToys(pol,channel,iy,options.infile)
                else:
                    asy_fit = valuesAndErrors['W_{pol}_{ch}_Ybin_{iy}_chargeasym'.format(pol=pol,ch=channel,iy=iy)]
                tmp_val.val_fit .append(asy_fit[0])
                tmp_val.elo_fit.append(abs(asy_fit[0]-asy_fit[1]))
                tmp_val.ehi_fit.append(abs(asy_fit[0]-asy_fit[2]))

                # on the charge asymmetry, which is A~0, better to show the difference 
                # Afit - Aexp wrt the ratio. The error on "exp" diff shows the error on Aexp, while the error bar the error on Afit
                tmp_val.relv. append(tmp_val.val[-1] - tmp_val.val_fit[-1])
                tmp_val.rello.append(chasy_val['asy'][1])
                tmp_val.relhi.append(chasy_val['asy'][1])

                tmp_val.relv_fit .append(0.)
                tmp_val.rello_fit.append(tmp_val.elo_fit[-1])
                tmp_val.relhi_fit.append(tmp_val.ehi_fit[-1])

                tmp_val.rap.append((ybins[cp][iy]+ybins[cp][iy+1])/2.)
                tmp_val.rlo.append(abs(ybins[cp][iy]-tmp_val.rap[-1]))
                tmp_val.rhi.append(abs(ybins[cp][iy]-tmp_val.rap[-1]))

            tmp_val.makeGraphs()
            asymmetryValues[pol] = tmp_val        
        plotValues(asymmetryValues,'asymmetry',channel,options)
            
        # now do the unpolarized ones
        tmp_val = valueClass('asymmetry_unpolarized')
        for iy,y in enumerate(ybinwidths['plus_left']): # this assumes that all the 3 polarizations have the same binning
            xval = {'plus': 0, 'minus': 0}; xerr = {'plus': 0, 'minus': 0}
            for charge in ['plus','minus']:
                for pol in polarizations:
                    xval[charge] += xsec_nominal_allCharges[charge][pol][iy]
                    xerr[charge] += pow(xsec_systematics_allCharges[charge][pol][iy],2)
                xerr[charge] = math.sqrt(xerr[charge])

            chasy_val = utilities.getChargeAsy(xval['plus'], xval['minus'],
                                               xerr['plus'], xerr['minus'])
            
            tmp_val.val .append(chasy_val['asy'][0])
            tmp_val.ehi.append(chasy_val['asy'][1])
            tmp_val.elo.append(chasy_val['asy'][1])

            if options.type == 'hessian': # should make the right expression from toys, if needed... 
                asy_fit = valuesAndErrors['W_{ch}_Ybin_{iy}_chargemetaasym'.format(ch=channel,iy=iy)]
            tmp_val.val_fit .append(asy_fit[0])
            tmp_val.elo_fit.append(abs(asy_fit[0]-asy_fit[1]))
            tmp_val.ehi_fit.append(abs(asy_fit[0]-asy_fit[2]))

            tmp_val.relv. append(tmp_val.val[-1] - tmp_val.val_fit[-1])
            tmp_val.rello.append(chasy_val['asy'][1])
            tmp_val.relhi.append(chasy_val['asy'][1])

            tmp_val.relv_fit .append(0.)
            tmp_val.rello_fit.append(tmp_val.elo_fit[-1])
            tmp_val.relhi_fit.append(tmp_val.ehi_fit[-1])

            tmp_val.rap.append((ybins[cp][iy]+ybins[cp][iy+1])/2.)
            tmp_val.rlo.append(abs(ybins[cp][iy]-tmp_val.rap[-1]))
            tmp_val.rhi.append(abs(ybins[cp][iy]-tmp_val.rap[-1]))

        tmp_val.makeGraphs()
        plotUnpolarizedValues(tmp_val,'asymmetry',channel,options)
            
