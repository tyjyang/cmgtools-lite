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
PRELIMINARY = '' # '#it{Preliminary}'

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

        # I tried the following two lines, but the next ones might be good as well
        #self.color  = ROOT.kBlue-7 if self.isleft else ROOT.kOrange+7 if self.isright else ROOT.kGray+2
        #self.colorf = ROOT.kBlue-4 if self.isleft else ROOT.kOrange+1 if self.isright else ROOT.kGray+3
        self.color  = ROOT.kBlue+2 if self.isleft else ROOT.kRed+1 if self.isright else ROOT.kGray+1
        self.colorf = ROOT.kAzure+1 if self.isleft else ROOT.kOrange+1 if self.isright else ROOT.kGray+3
        if self.isunpolarized: 
            self.color = ROOT.kSpring-6
            self.colorf = ROOT.kSpring-7

        ## here all the arrays that will contain the values and errors etc.
        self.val       = array.array('f', []); self.ehi       = array.array('f', []); self.elo       = array.array('f', []); self.ehi2   = array.array('f', []); self.elo2   = array.array('f', []);
        self.val_fit   = array.array('f', []); self.ehi_fit   = array.array('f', []); self.elo_fit   = array.array('f', []); 
        self.relv      = array.array('f', []); self.relhi     = array.array('f', []); self.rello     = array.array('f', []); self.relhi2 = array.array('f', []); self.rello2 = array.array('f', []);
        self.relv_fit  = array.array('f', []); self.relhi_fit = array.array('f', []); self.rello_fit = array.array('f', []);
        self.rap       = array.array('f', []); self.rlo       = array.array('f', []); self.rhi       = array.array('f', []);
        self.altval    = array.array('f', [])

    def makeGraphs(self):
        if len(self.val):
            self.graph = ROOT.TGraphAsymmErrors(len(self.val), self.rap, self.val, self.rlo, self.rhi, self.elo, self.ehi)
            self.graph.SetName('graph'+self.pol)
            self.graph.SetTitle('')
        if len(self.val) and len(self.ehi2) and len(self.elo2):
            self.graph2 = ROOT.TGraphAsymmErrors(len(self.val), self.rap, self.val, self.rlo, self.rhi, self.elo2, self.ehi2)
            self.graph2.SetName('graph2'+self.pol)
            self.graph2.SetTitle('')
        if len(self.relv):
            self.graph_rel= ROOT.TGraphAsymmErrors(len(self.relv), self.rap, self.relv, self.rlo, self.rhi, self.rello, self.relhi)
            self.graph_rel.SetName('graph'+self.pol+'_rel')
            self.graph_rel.SetTitle('')
        if len(self.relv) and len(self.relhi2) and len(self.rello2):
            self.graph2_rel= ROOT.TGraphAsymmErrors(len(self.relv), self.rap, self.relv, self.rlo, self.rhi, self.rello2, self.relhi2)
            self.graph2_rel.SetName('graph2'+self.pol+'_rel')
            self.graph2_rel.SetTitle('')
        if len(self.val_fit):
            self.graph_fit = ROOT.TGraphAsymmErrors(len(self.val_fit), self.rap, self.val_fit, self.rlo, self.rhi, self.elo_fit, self.ehi_fit)
            self.graph_fit.SetName('graph'+self.pol+'_fit')
            self.graph_fit.SetTitle('')
        zeros = array.array('f',[0 for i in xrange(len(self.rlo))])
        if len(self.relv_fit):
            self.graph_fit_rel = ROOT.TGraphAsymmErrors(len(self.relv_fit), self.rap, self.relv_fit, zeros, zeros, self.rello_fit, self.relhi_fit)
            self.graph_fit_rel.SetName('graph'+self.pol+'_fit_rel')
            self.graph_fit_rel.SetTitle('')
        if len(self.altval):
            self.altgraph = ROOT.TGraph(len(self.val), self.rap, self.altval)
            self.graph.SetName('altgraph'+self.pol)
            self.graph.SetTitle('')       

        self.graphStyle()
        if len(self.relv) and len(self.relv_fit): self.makeMultiGraphRel()

    def makeMultiGraphRel(self):
        self.mg = ROOT.TMultiGraph()
        self.mg.SetTitle() ## no title 'W^{{{ch}}}: {p}'.format(ch=self.ch,p=self.pol))
        self.shiftPoints(self.graph_fit_rel)
        self.mg.Add(self.graph_rel,'P2')
        self.mg.Add(self.graph_fit_rel)
        # if len(self.relhi2) and len(self.rello2):
        #     self.mg.Add(self.graph2_rel,'P2')

    def graphStyle(self):
        #fillstyles = {'left': 3244, 'right': 3001, 'long': 3144, 'unpolarized': 3001}
        #fillstyles_rel = {'left': 3244, 'right': 3001, 'long': 3144, 'unpolarized': 3001}
        #fillstyles = {'left': 3244, 'right': 3244, 'long': 3244, 'unpolarized': 3244}
        #fillstyles_rel = {'left': 3444, 'right': 3444, 'long': 3444, 'unpolarized': 3244}
        fillstyles     = {'left': 1001, 'right': 1001, 'long': 1001, 'unpolarized': 1001}
        fillstyles_rel = {'left': 1001, 'right': 1001, 'long': 1001, 'unpolarized': 1001}
        if hasattr(self,'graph'):
            self.graph.SetLineColor(self.color)
            self.graph.SetFillColorAlpha(self.colorf,0.30)
            self.graph.SetFillStyle(fillstyles[self.pol])
        if hasattr(self,'graph2'):
            self.graph2.SetLineColor(self.color)
            self.graph2.SetFillColor(self.colorf)
            self.graph2.SetFillStyle(fillstyles[self.pol])
        if hasattr(self,'graph_fit'):
            self.graph_fit.SetLineWidth(3)
            self.graph_fit.SetMarkerSize(1.0)
            self.graph_fit.SetMarkerStyle(ROOT.kFullCircle)
            self.graph_fit.SetMarkerColor(self.color)
            self.graph_fit.SetLineColor(self.color)
        if hasattr(self,'graph_rel'):
            self.graph_rel.SetLineWidth(5)
            self.graph_rel.SetLineColor(self.color)
            self.graph_rel.SetFillColorAlpha(self.colorf,0.30)
            self.graph_rel.SetFillStyle(fillstyles_rel[self.pol])
        if hasattr(self,'graph2_rel'):
            self.graph2_rel.SetLineWidth(5)
            self.graph2_rel.SetLineColor(self.color)
            self.graph2_rel.SetFillColor(self.colorf)
            self.graph2_rel.SetFillStyle(fillstyles_rel[self.pol])
        if hasattr(self,'graph_fit_rel'):
            self.graph_fit_rel.SetLineWidth(2)
            self.graph_fit_rel.SetMarkerSize(1.0)
            self.graph_fit_rel.SetMarkerStyle(ROOT.kFullCircle)
            self.graph_fit_rel.SetLineColor(self.color)
            self.graph_fit_rel.SetMarkerColor(self.color)
            self.graph_fit_rel.SetFillColor(ROOT.kGreen+3)
            self.graph_fit_rel.SetFillStyle(3001)
        if hasattr(self,'altgraph'):
            self.altgraph.SetLineColor(self.colorf+1)
            self.altgraph.SetLineWidth(2)
            self.altgraph.SetFillColor(self.colorf+1)
            self.altgraph.SetFillStyle(fillstyles[self.pol])

    def shiftPoints(self, graph):
        shift = 0.25/8.
        shifts = {'left': -shift, 'right': shift, 'long': 0.0, 'unpolarized': 0.0}
        for p in xrange(graph.GetN()):
            x = ROOT.Double(0); y = ROOT.Double(0)
            graph.GetPoint(p,x,y)
            graph.SetPoint(p,x+shifts[self.pol],y)
            graph.SetPointEXhigh(p,0); graph.SetPointEXlow(p,0)

def plotValues(values,charge,channel,options, polarizations=['left','right','long']):
        c2 = ROOT.TCanvas('foo','', 800, 800)
        c2.GetPad(0).SetTopMargin(0.09)
        c2.GetPad(0).SetBottomMargin(0.35)
        c2.GetPad(0).SetLeftMargin(0.17)
        c2.GetPad(0).SetRightMargin(0.04)
        c2.GetPad(0).SetTickx(1)
        c2.GetPad(0).SetTicky(1)

        skipLong = False
        if options.nolong: skipLong = True
        doAltExp = (len(values['left'].altval)>0)

        ch = '#plus' if charge == 'plus' else '#minus'
        if charge == 'asymmetry': ch = ''
        date = datetime.date.today().isoformat()
        normstr = 'norm' if (options.normxsec and charge!='asymmetry') else ''

        lat = ROOT.TLatex()
        lat.SetNDC(); lat.SetTextFont(42)
        ## the four graphs exist now. now starting to draw them
        ## ===========================================================
        if sum(hasattr(values[pol],'graph') and hasattr(values[pol],'graph_fit') for pol in polarizations)==len(polarizations):
            leg = ROOT.TLegend(0.2, 0.78 if skipLong else 0.75, 0.93, 0.88)
        #if sum(hasattr(values[pol],'graph') and hasattr(values[pol],'graph_fit') for pol in ['left','right','long'])==3:
        #    leg = ROOT.TLegend(0.40, 0.80, 0.90, 0.90)
            leg.SetFillStyle(0)
            leg.SetBorderSize(0)
            leg.AddEntry(values['left'] .graph_fit , 'W_{L} (fit)', 'pl')
            leg.AddEntry(values['left'] .graph     , 'W_{{L}} ({mc})'.format(mc=REFMC) , 'f')
            if doAltExp:
                leg.AddEntry(values['left'] .altgraph     , 'W_{{L}} ({mc}*)'.format(mc=REFMC) , 'l')
            leg.AddEntry(values['right'].graph_fit , 'W_{R} (fit)', 'pl') 
            leg.AddEntry(values['right'].graph     , 'W_{{R}} ({mc})'.format(mc=REFMC) , 'f')
            if doAltExp:
                leg.AddEntry(values['right'].altgraph     , 'W_{{R}} ({mc}*)'.format(mc=REFMC) , 'l')
            leg.SetNColumns(3 if doAltExp else 3)
            if not skipLong:
                leg.AddEntry(values['long'] .graph     , 'W_{{0}} ({mc})'.format(mc=REFMC) , 'f')
                leg.AddEntry(values['long'] .graph_fit , 'W_{0} (fit)', 'pl')

            values['left'].graph.SetTitle('W {ch}: y_{{W}}'.format(ch=ch))
                
            mg = ROOT.TMultiGraph()
            mg.Add(values['left'] .graph,'E2'); #mg.Add(values['left'] .graph2,'E2')
            mg.Add(values['right'].graph,'E2'); #mg.Add(values['right'].graph2,'E2')
            if not skipLong: mg.Add(values['long'] .graph,'E2'); #mg.Add(values['long'] .graph2,'E2')
            if doAltExp:
                mg.Add(values['left'] .altgraph,'L2')
                mg.Add(values['right'].altgraph,'L2')
            mg.Add(values['left'] .graph_fit)
            mg.Add(values['right'].graph_fit)
            if not skipLong: mg.Add(values['long'] .graph_fit)
     
            mg.Draw('Pa')
            mg.GetXaxis().SetRangeUser(0., options.maxRapidity) # max would be 6.
            mg.GetXaxis().SetTitle('')
            mg.GetXaxis().SetLabelSize(0)
            if charge=='asymmetry':
                mg.GetYaxis().SetTitle('Charge asymmetry')
                # -0.05,0.45 instead of -0.1,0.4 not to overlap with legend
                mg.GetYaxis().SetRangeUser(-0.05,0.45) 
            else:
                if options.normxsec: 
                    #mg.GetYaxis().SetTitle('#frac{d#sigma}{#sigma_{tot}^{fit}} / d|Y_{W}|')
                    mg.GetYaxis().SetTitle('d#sigma / d|y_{W}| / #sigma_{tot}')
                    mg.GetYaxis().SetRangeUser(0.,0.8 if options.maxRapidity > 2.9 else 0.25)
                else: 
                    mg.GetYaxis().SetTitle('d#sigma / d|y_{W}| (pb)')
                    mg.GetYaxis().SetRangeUser(-200,3500)
            mg.GetYaxis().SetTitleSize(0.04)
            mg.GetYaxis().SetLabelSize(0.04)
            mg.GetYaxis().SetTitleOffset(2.0)
            mg.GetYaxis().CenterTitle()
            mg.GetYaxis().SetDecimals()
     
            leg.Draw('same')

        # save for HEPData
        plotname = "{od}/genAbsY{norm}_pdfs_{ch}{suffix}_{t}".format(od=options.outdir, norm=normstr, ch=charge, suffix=options.suffix, t=options.type)
        rfHepName = '{n}.root'.format(n=plotname)
        rfHep = ROOT.TFile.Open(rfHepName,"recreate")
        if not rfHep:
            print "Error in plotYW.py: could not open root file %s" % rfHepName
            quit()
        rfHep.cd()
        polsToSave = ["left","right"] + ([] if skipLong else ["long"])
        for pts in polsToSave:
            values[pts].graph.Clone().Write("exp_{p}".format(p=pts))
            values[pts].graph_fit.Clone().Write("data_{p}".format(p=pts))
            if doAltExp:
                values[pts].altgraph.Clone().Write("expAlt_{p}".format(p=pts))
        rfHep.Close()
        # done with hepdata
     
        ## now make the relative error plot:
        ## ======================================
        if sum(hasattr(values[pol],'mg') for pol in polarizations)==len(polarizations):

            pad2 = ROOT.TPad("pad2","pad2",0,0.,1,0.9)
            pad2.SetTopMargin(0.65)
            pad2.SetRightMargin(0.04)
            pad2.SetLeftMargin(0.17)
            pad2.SetBottomMargin(0.14)
            pad2.SetFillColor(0)
            pad2.SetGridy(0)
            pad2.SetFillStyle(0)
            pad2.SetTicky(1)
            pad2.SetTickx(1)

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
                yaxrange = (0.80, 1.20)
     

            helToPlot = ['left','right']
            if not skipLong:  helToPlot.append('long')
            for  ih,hel in enumerate(helToPlot):
     
                values[hel].mg.Draw('Pa' if ih==0 else 'P')
                if ih==0:
                    ## x axis fiddling
                    values[hel].mg.GetXaxis().SetRangeUser(0., options.maxRapidity)
                    values[hel].mg.GetXaxis().SetLabelSize(0.04)
                    ## y axis fiddling
                    values[hel].mg.GetYaxis().SetTitleOffset(1.8)
                    values[hel].mg.GetYaxis().SetTitleSize(0.045)
                    values[hel].mg.GetYaxis().SetLabelSize(0.04)
                    values[hel].mg.GetYaxis().SetTitle(yaxtitle)
                    values[hel].mg.GetYaxis().SetRangeUser(yaxrange[0],yaxrange[1])
                    values[hel].mg.GetYaxis().SetNdivisions(4)
                    values[hel].mg.GetYaxis().CenterTitle()
                    values[hel].mg.GetYaxis().SetDecimals()
            line.Draw("Lsame");
            c2.cd()
            lat.DrawLatex(0.18, 0.94, '#bf{{CMS}} {prel}'.format(prel=PRELIMINARY))
            lat.DrawLatex(0.62, 0.94, '35.9 fb^{-1} (13 TeV)')
            flavor = "#mu" if channel == "mu" else "e" if channel=='el' else 'l'
            if charge == 'asymmetry':
                lat.DrawLatex(0.2, 0.63,  'W^{{{ch}}} #rightarrow {lep}^{{{ch} }}{nu}'.format(ch=ch,lep=flavor,nu="#bar{#nu}" if charge=='minus' else "#nu"))
            else:
                lat.DrawLatex(0.2, 0.43,  'W^{{ {ch}}} #rightarrow {lep}^{{ {ch} }}{nu}'.format(ch=ch,lep=flavor,nu="#bar{#nu}" if charge=='minus' else "#nu"))

            lat.DrawLatex(0.88, 0.03, '|y_{W}|')
        for ext in ['png', 'pdf']: #, 'root']:
            c2.SaveAs('{n}.{ext}'.format(n=plotname, ext=ext))


def plotUnpolarizedValues(values,charge,channel,options):
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
        doAltExp = (len(values.altval)>0)
        
        valkey = values.name.split('_')[1]

        lat = ROOT.TLatex()
        lat.SetNDC(); lat.SetTextFont(42)
        legx1, legx2, legy1, legy2 = 0.2, 0.5, 0.7, 0.85

        #of = ROOT.TFile.Open('{od}/genAbsYUnpolarized{norm}_pdfs_{ch}{suffix}_{t}.root'.format(od=options.outdir, norm=valkey, ch=charge, suffix=options.suffix, t=options.type),'recreate')

        # save for HEPData
        # will also save the TMultigraph used in plotYW_FEWZ.py
        plotname = "{od}/genAbsYUnpolarized{norm}_pdfs_{ch}{suffix}_{t}".format(od=options.outdir, norm=valkey, ch=charge, suffix=options.suffix, t=options.type)
        rfHepName = '{n}.root'.format(n=plotname)
        rfHep = ROOT.TFile.Open(rfHepName,"recreate")
        if not rfHep:
            print "Error in plotYW.py: could not open root file %s" % rfHepName
            quit()
        rfHep.cd()

        ## the graphs exist now. now starting to draw them
        ## ===========================================================
        if hasattr(values,'graph') and hasattr(values,'graph_fit'):
                
            mg = ROOT.TMultiGraph()
            mg.SetName('values')
            mg.Add(values.graph,'E2'); #mg.Add(values.graph2,'E2')
            mg.Add(values.graph_fit)
            if doAltExp:
                mg.Add(values.altgraph,'L2')
            mg.Draw('Pa')
            mg.GetXaxis().SetRangeUser(0., options.maxRapidity) # max would be 6.
            mg.GetXaxis().SetTitle('')
            mg.GetXaxis().SetTitleOffset(5.5)
            mg.GetXaxis().SetLabelSize(0)
            titles = {'asymmetry': 'Charge asymmetry',
                      'a0': 'A_{0}', 
                      'a4': '-A_{4}' if charge == 'plus' else 'A_{4}',
                      'sumxsec': 'd#sigma / d|y_{W}| (pb)',
                      'sumxsecnorm': 'd#sigma/d|y_{W}| / #sigma_{tot}'}
            ranges = {'asymmetry': (-0.1,0.4),
                      'a0': (0.07,0.2),
                      #'a4': (-1,2),
                      'a4': (-0.5,1.2),
                      'sumxsec': (1500,4500),
                      'sumxsecnorm': (0.1,0.5)}
            mg.GetYaxis().SetTitle(titles[valkey])
            mg.GetYaxis().SetRangeUser(ranges[valkey][0],ranges[valkey][1])

            mg.GetYaxis().SetTitleSize(0.04)
            mg.GetYaxis().SetLabelSize(0.04)
            mg.GetYaxis().SetTitleOffset(2.0)
            mg.GetYaxis().CenterTitle()
            mg.GetYaxis().SetDecimals()

            leg = ROOT.TLegend(legx1, legy1, legx2, legy2)
            leg.SetFillStyle(0)
            leg.SetBorderSize(0)
            leg.AddEntry(values.graph_fit , 'Measured', 'pl')
            leg.AddEntry(values.graph     , REFMC, 'f')
            if doAltExp:
                leg.AddEntry(values.graph     , REFMC+'*', 'l')

            leg.Draw('same')
            lat.DrawLatex(0.18, 0.94, '#bf{{CMS}} {prel}'.format(prel=PRELIMINARY))
            lat.DrawLatex(0.62, 0.94, '35.9 fb^{-1} (13 TeV)')
            flavor = "#mu" if channel == "mu" else "e" if channel=='el' else 'l'
            if charge == 'asymmetry':
                lat.DrawLatex(0.20, 0.40,  'W^{{{ch}}} #rightarrow {lep}^{{{ch} }}{nu}'.format(ch=ch,lep=flavor,nu="#bar{#nu}" if charge=='minus' else "#nu"))
            else:
                lat.DrawLatex(0.20, 0.40,  'W^{{ {ch}}} #rightarrow {lep}^{{ {ch} }}{nu}'.format(ch=ch,lep=flavor,nu="#bar{#nu}" if charge=='minus' else "#nu"))

            lat.DrawLatex(0.88, 0.03, '|y_{W}|')
            mg.Write() 

        # save for HepData, but for A4 in W+ I need to take graph changing sign to all points
        # because A4 is negative for W+, but the fit returns a positive number
        # for HepData I want to report the actual values as predicted from theory, even if
        # the plot has -A4 to make the comparison with minus charge easier
        # for the asymmetric errors, what is ErrorYhigh becomes ErrorYlow and viceversa
        if valkey == "a4" and charge == "plus":
            gA4_exp = values.graph.Clone("gA4_exp") 
            gA4_data = values.graph_fit.Clone("gA4_data") 
            if doAltExp:
                gA4_altexp = values.altgraph.Clone("gA4_altexp") 
            for ip in range(gA4_data.GetN()):
                xval = ROOT.Double(0)
                yval_data = ROOT.Double(0)  
                yval_exp = ROOT.Double(0)  
                gA4_data.GetPoint(ip,xval,yval_data)  
                gA4_data.SetPoint(ip,xval,-1.0*yval_data)
                tmperrhigh = gA4_data.GetErrorYhigh(ip) 
                gA4_data.SetPointEYhigh(ip, gA4_data.GetErrorYlow(ip))
                gA4_data.SetPointEYlow(ip, tmperrhigh)
                gA4_exp.GetPoint(ip,xval,yval_exp)  
                gA4_exp.SetPoint(ip,xval,-1.0*yval_exp)
                tmperrhigh = gA4_exp.GetErrorYhigh(ip) 
                gA4_exp.SetPointEYhigh(ip, gA4_exp.GetErrorYlow(ip))
                gA4_exp.SetPointEYlow(ip, tmperrhigh)
                if doAltExp:
                    yval_altexp = ROOT.Double(0)
                    gA4_altexp.GetPoint(ip,xval,yval_altexp)  
                    gA4_altexp.SetPoint(ip,xval,-1.0*yval_altexp)
            # now saves
            gA4_exp.Write("exp")
            gA4_data.Write("data")
            if doAltExp:
                gA4_altexp.Write("expAlt")
        else:
            values.graph.Clone().Write("exp")
            values.graph_fit.Clone().Write("data")
            if doAltExp:
                values.altgraph.Clone().Write("expAlt")

        ## now make the relative error plot:
        ## ======================================
        if hasattr(values,'mg'):

            pad2 = ROOT.TPad("pad2","pad2",0,0.,1,0.9)
            pad2.SetTopMargin(0.65)
            pad2.SetRightMargin(0.04)
            pad2.SetLeftMargin(0.17)
            pad2.SetBottomMargin(0.14)
            pad2.SetFillColor(0)
            pad2.SetGridy(0)
            pad2.SetFillStyle(0)
            pad2.SetTicky(1)
            pad2.SetTickx(1)

            pad2.Draw()
            pad2.cd()

            line = ROOT.TF1("horiz_line","0" if valkey in ['asymmetry','a4'] else '1',0.0,3.0);
            line.SetLineColor(ROOT.kBlack);
            line.SetLineWidth(2);
            yaxrange = (0,0)
            if valkey in ['asymmetry','a4']:
                yaxtitle = 'A_{Theory}-A_{Data}'
                if valkey == 'a4':
                    yaxrange = (-0.2, 0.2)
                else:
                    yaxrange = (-0.1, 0.1)
            else:
                yaxtitle = '#sigma_{Theory}/#sigma_{Data}'
                yaxrange = (0.80, 1.20) if 'xsec' in valkey else (0.8,1.2) if valkey=='a0' else (-0.5, 2.5)

            values.mg.SetName('ratios')
            values.mg.Draw('Pa')
            ## x axis fiddling
            values.mg.GetXaxis().SetTitle('')
            values.mg.GetXaxis().SetRangeUser(0., options.maxRapidity)
            values.mg.GetXaxis().SetTitleSize(0.14)
            values.mg.GetXaxis().SetLabelSize(0.04)
            ## y axis fiddling
            values.mg.GetYaxis().SetTitleOffset(1.8)
            values.mg.GetYaxis().SetTitleSize(0.045)
            values.mg.GetYaxis().SetLabelSize(0.04)
            values.mg.GetYaxis().SetTitle(yaxtitle)
            values.mg.GetYaxis().SetRangeUser(yaxrange[0],yaxrange[1])
            values.mg.GetYaxis().SetNdivisions(5)
            values.mg.GetYaxis().CenterTitle()
            values.mg.GetYaxis().SetDecimals()
            values.mg.Write()

            line.Draw("Lsame");

        rfHep.Close()
        # done with hepdata

        for ext in ['png', 'pdf']:
            c2.SaveAs('{n}.{ext}'.format(n=plotname, ext=ext))

NPDFs = 60
#LUMINOSITY = 35900
LUMINOSITY = 36000

def getGraphsTheoryXsecPrefit(wptReweighted=False, pdfOnly=False):
    # file with all theory bands prefit, without or with Wpt reweighting, normalized to native aMC@NLO xsec
    # made with w-helicity-13TeV/plotHelicityChargeAsymmetry.py if you need to redo it
    pdfpostfix = "_pdfOnly" if pdfOnly else ""
    wptpostfix = "_wptWeights" if wptReweighted else ""

    fname = "/afs/cern.ch/user/m/mciprian/www/wmass/13TeV/helicityAnalysis/checkTheoryBands{w}{p}/plotHelicityChargeAsymmetry.root".format(w=wptpostfix,p=pdfpostfix)

    f = ROOT.TFile.Open(fname,"READ")
    graphs = {}
    for k in f.GetListOfKeys():
        name = k.GetName()
        obj  = k.ReadObj()
        if not obj:
            raise RuntimeError('Unable to read object {n} in {fn}'.format(n=name,fn=fname))
        graphs[name] = obj
    f.Close()
    #
    # check the graphs are still in the dictionary, exit after first entry
    # yes it does, can comment
    # for ik,k in enumerate(graphs.keys()):
    #     if ik: 
    #         break
    #     graphs[k].Print()
        
    return graphs

if __name__ == "__main__":

    from optparse import OptionParser
    parser = OptionParser(usage='%prog ybinfile workspace.root toys.root [options] ')
    parser.add_option('-i', '--infile'      , dest='infile'   , default=''            , type='string', help='workspace converted from datacard')
    parser.add_option('-y', '--ybinfile'    , dest='ybinfile' , default=''            , type='string', help='file with the yw binning')

    parser.add_option('-t', '--type'        , dest='type'     , default='toys'        , type='string', help='run the plot from which postfit? toys/scans/hessian')
    parser.add_option(      '--toyfile'     , dest='toyfile'  , default=''            , type='string', help='file that has the toys')
    parser.add_option(      '--scandir'     , dest='scandir'  , default=''            , type='string', help='directory with all the scans')
    parser.add_option(      '--hessfile'    , dest='hessfile' , default=''            , type='string', help='file that contains the hessian errors in a dictionary')
    parser.add_option(      '--xsecfiles'   , dest='xsecfiles' , default=None         , type='string', help='files that contains the expected x sections with variations (one per charge,comma separated in the same order of the charges) ')
    parser.add_option(      '--altxsecfiles', dest='altxsecfiles' , default=None      , type='string', help='files that contains alternative xsecs as the nominal to be drawn as a line ')
    parser.add_option('-C', '--charge'      , dest='charge'   , default='plus,minus'  , type='string', help='process given charge. default is both')
    parser.add_option('-o', '--outdir'      , dest='outdir'   , default='.'           , type='string', help='outdput directory to save the plots')
    parser.add_option(      '--suffix'      , dest='suffix'   , default=''            , type='string', help='suffix for the correlation matrix')
    parser.add_option('-n', '--normxsec'    , dest='normxsec' , default=False         , action='store_true',   help='if given, plot the differential xsecs normalized to the total xsec')
    parser.add_option(      '--nolong'      , dest='nolong'   , default=False         , action='store_true',   help='if given, do not plot longitudinal component (but it assumes the POIs exist)')
    parser.add_option(      '--longBkg'     , dest='longBkg'  , default=False         , action='store_true',   help='if True, longitudinal component was treated as background, so the POIs are missing. Manage inputs accordingly')
    parser.add_option(     '--ybinsBkg', dest='ybinsBkg', type='string', default="10,11", help='Define which Y bins are to be considered as background. With format 14,15 ')
    parser.add_option(     '--ybinsOutAcc', dest='ybinsOutAcc', type='string', default="", help='Define which Y bins were put in OutAcc channel in the fit. With format 14,15 ')
    parser.add_option(      '--max-rap'     , dest='maxRapidity', default='2.5'       , type='float', help='Max value for rapidity range')
    parser.add_option(      '--hessianFromToy', default=0       , type=int, help='get entry from toy file to be hessian like')
    parser.add_option(      '--pdfonly'      , dest='pdfonly'   , default=False         , action='store_true',   help='if given, do not include alphaS and QCD scales in the theory bands')
    (options, args) = parser.parse_args()


    theoryBands = getGraphsTheoryXsecPrefit(wptReweighted=True,pdfOnly=options.pdfonly) 
    # graphs from 0 to 3.0 with 12 bins of 0.25 width, but the last two are out of acceptance and ignored (the last one should have been from 2.75 to 10 actually, but we don't care), the normalized xsec is made using all the 12 bins for all polarizations
    theoryBandsNoWpt = getGraphsTheoryXsecPrefit(wptReweighted=False,pdfOnly=options.pdfonly)
    nBinsTheoryGraph = 12
    hGrNameKey = "PDF" if options.pdfonly else "TotTheory"

    # should have used bands from Xsec without Wpt reweighting, but apparently the PDF uncertainty is smaller than expected, and this is probably due to the fact that the input cross sections were not made using appropriate helicity fractions for each pdf variation
    # instead, this is apparently done for the Wpt reweighted xsec (although I don't get exactly the same numbers as Josh, but very close)
    # so, the central value will be taken from the non-reweighted xsec, but the band from the reweighted one (it has slightly larger pdf uncertainty, and probably this is more correct)

    if not os.path.isdir(options.outdir):
        os.system('mkdir -p {od}'.format(od=options.outdir))
        if os.path.exists("/afs/cern.ch"): os.system("cp /afs/cern.ch/user/m/mciprian/public/index.php {od}".format(od=options.outdir))

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
        valuesAndErrors = utilities.getFromHessian(options.infile, takeEntry=options.hessianFromToy)
    else:
        print 'ERROR: none of your types is supported. specify either "toys", "scans", or "hessian"'
        sys.exit()

    ybinfile = open(ybinfile, 'r')
    ybins = eval(ybinfile.read())
    ybinfile.close()
    #print "YBINS: "
    #print ybins

    bkgYBins = []
    if options.ybinsBkg:
        bkgYBins = list(int(i) for i in options.ybinsBkg.split(','))        

    outAccYBins = []
    if options.ybinsOutAcc:
        outAccYBins = list(int(i) for i in options.ybinsOutAcc.split(','))

    if len(outAccYBins) and len(bkgYBins):
        print "Warning: I see there are rapidity bins treated as outAcc and/or background."
        print "Although this is not impossible, I suspect you are messing up with options --ybinsBkg and --ybinsOutAcc"
        print "Please check!"
        quit()

    ## calculate the bin widths for the rapidity bins
    ybinwidths = {}
    for k,v in ybins.items():
        tmplist = list(abs(i - v[v.index(i)+1]) for i in v[:-1])
        ybinwidths[k] = [float('{n:.2f}'.format(n=i)) for i in tmplist]
    #print "YBINWIDTHS: "
    #print ybinwidths

    charges = options.charge.split(',')
    xsecfiles = options.xsecfiles.split(',')
    doAltExp =  options.altxsecfiles
    xsec_nominal_allCharges = {}; 
    xsecnorm_nominal_allCharges = {}; 

    if doAltExp:
        alt_xsecfiles = options.altxsecfiles.split(',')
        alt_xsec_nominal_allCharges = {}
        alt_xsecnorm_nominal_allCharges = {}

    polarizations = ['left','right','long']
    signal_polarizations = ['left','right']
    if not options.longBkg:
        signal_polarizations.append('long')

    if 'lep' in os.path.basename(xsecfiles[0]):
        nChan = 2
        channel = 'lep'
    else:
        nChan = 1
        channel = 'mu' if 'mu' in os.path.basename(xsecfiles[0]) else 'el'
    print "From the xsec file names it seems that you are plotting results for channel ",channel

    for ic,charge in enumerate(charges):

        sign = 1. if charge=='plus' else -1.

        xsec_nominal = {}
        xsecnorm_nominal = {}
        alt_xsec_nominal = {}
        alt_xsecnorm_nominal = {}
        for pol in ["left","right","long"]:
            nameGraphTotTheory_xsec = "hXsec{n}_{c}_{p}".format(n=hGrNameKey,c=charge,p=pol)
            nameGraphTotTheory_xsecnorm = "hXsecNorm{n}_{c}_{p}".format(n=hGrNameKey,c=charge,p=pol)
            values = []
            valuesnorm = []
            for iy in range(nBinsTheoryGraph): # all bins of yW
                xval = ROOT.Double(0)
                yval = ROOT.Double(0)  
                yvalnorm = ROOT.Double(0)  
                theoryBandsNoWpt[nameGraphTotTheory_xsec].GetPoint(iy,xval,yval)  
                values.append(yval)
                theoryBandsNoWpt[nameGraphTotTheory_xsecnorm].GetPoint(iy,xval,yvalnorm)  
                valuesnorm.append(yvalnorm)
            xsec_nominal[pol] = values
            xsecnorm_nominal[pol] = valuesnorm
        xsec_nominal_allCharges[charge] = xsec_nominal        
        xsecnorm_nominal_allCharges[charge] = xsecnorm_nominal        

        if doAltExp: 
            for pol in ["left","right","long"]:
                altvalues = []
                altvaluesnorm = []
                nameGraphTotTheory_xsec = "hXsec{n}_{c}_{p}".format(n=hGrNameKey,c=charge,p=pol)
                nameGraphTotTheory_xsecnorm = "hXsecNorm{n}_{c}_{p}".format(n=hGrNameKey,c=charge,p=pol)
                for iy in range(nBinsTheoryGraph): # all bins of yW
                    xval = ROOT.Double(0)
                    altyval = ROOT.Double(0)  
                    altyvalnorm = ROOT.Double(0)  
                    theoryBands[nameGraphTotTheory_xsec].GetPoint(iy,xval,altyval)  
                    altvalues.append(altyval)
                    theoryBands[nameGraphTotTheory_xsecnorm].GetPoint(iy,xval,altyvalnorm)  
                    altvaluesnorm.append(altyvalnorm)
                alt_xsec_nominal[pol] = altvalues
                alt_xsecnorm_nominal[pol] = altvaluesnorm            
            alt_xsec_nominal_allCharges[charge] = alt_xsec_nominal
            alt_xsecnorm_nominal_allCharges[charge] = alt_xsecnorm_nominal
            # print "Charge %s " % charge
            # print "alt_xsec_nominal[left]:"
            # print alt_xsec_nominal["left"]
            # print "alt_xsec_nominal[right]:"
            # print alt_xsec_nominal["right"]

        angcoeff_nominal = {'sumxsec': [], 'a0': [], 'a4': []}; alt_angcoeff_nominal = {'sumxsec': [], 'a0': [], 'a4': []}

        nOuterBinsToExclude = 0  ### out of acceptance Y bins, or that were treated as background (not to be considered for the total xsec)
        if len(outAccYBins):
            nOuterBinsToExclude = len(outAccYBins)
        print "number of outer bins to exclude: " + str(nOuterBinsToExclude)

        cp = '{ch}_left'.format(ch=charge)
        MAXYFORNORM = ybins[cp][-nOuterBinsToExclude-1] # exclude the outermost 2 bins which has huge error due to acceptance
        print "MAXYFORNORM = " + str(MAXYFORNORM)
        normsigmaIn  = sum([xsec_nominal[allpol][iy] for allpol in polarizations for iy,y in enumerate(ybins[cp][:-1]) if abs(y)<MAXYFORNORM])
        print "NORMSIGMAIN: " + str(normsigmaIn)
        normsigmaOut = sum([xsec_nominal[allpol][iy] for allpol in polarizations for iy,y in enumerate(ybins[cp][:-1]) if abs(y)>=MAXYFORNORM])
        if doAltExp:
            alt_normsigmaIn  = sum([alt_xsec_nominal[allpol][iy] for allpol in polarizations for iy,y in enumerate(ybins[cp][:-1]) if abs(y)<MAXYFORNORM])
            alt_normsigmaOut = sum([alt_xsec_nominal[allpol][iy] for allpol in polarizations for iy,y in enumerate(ybins[cp][:-1]) if abs(y)>=MAXYFORNORM])

        allValues = {}
        for pol in signal_polarizations:
            print "total expected xsec up to |Y|<{maxy} = {sigma:.3f} pb".format(maxy=MAXYFORNORM,sigma=normsigmaIn)
            if len(outAccYBins):
                print "total expected xsec beyond |Y|>{maxy} = {sigma:.3f} pb".format(maxy=MAXYFORNORM,sigma=normsigmaOut)

            tmp_val = valueClass('values_'+charge+'_'+pol)

            for iy,y in enumerate(ybinwidths['{ch}_{pol}'.format(ch=charge,pol=pol)]):
                if iy in bkgYBins: continue

                normsigma = normsigmaIn if abs(ybins[cp][iy])<MAXYFORNORM else normsigmaOut
                if doAltExp: alt_normsigma = alt_normsigmaIn if abs(ybins[cp][iy])<MAXYFORNORM else alt_normsigmaOut
                parname = 'W{charge}_{pol}_Ybin_{iy}'.format(charge=charge,pol=pol,iy=iy)

                scale = 1.
                if options.longBkg:
                    suffix = 'pmaskedexp' if pol!='long' else 'sumxsec'
                else: 
                    suffix = 'pmaskedexp'
                if options.normxsec:
                    if   options.type == 'toys': 
                        xsec_fit = utilities.getNormalizedXsecFromToys(ybins,charge,pol,channel,iy,options.infile,MAXYFORNORM)
                    elif options.type == 'hessian':
                        xsec_fit = [x for x in valuesAndErrors['{par}_{sfx}norm'.format(par=parname,sfx=suffix)]]
                    else:
                        print "--normxsec not implemented yet for scans."
                        sys.exit()
                else:
                    xsec_fit = [x/float(nChan) for x in valuesAndErrors['{par}_{sfx}'.format(par=parname,sfx=suffix)]]
                    scale = LUMINOSITY

                if doAltExp: 
                    if options.normxsec:
                        #tmp_val.altval[-1] = tmp_val.altval[-1]/alt_normsigma
                        tmp_val.altval.append(alt_xsecnorm_nominal[pol][iy]/ybinwidths[cp][iy])
                        #print "iy = %d --> alt_val %s = %.3f" % (iy,pol,tmp_val.altval[-1])
                    else:
                        tmp_val.altval.append(alt_xsec_nominal[pol][iy]/ybinwidths[cp][iy])
                        #print "iy = %d --> alt_val %s = %.3f" % (iy,pol,tmp_val.altval[-1])

                # abs xsec
                nameGraphTotTheory_xsec = "hXsec{n}_{c}_{p}".format(n=hGrNameKey,c=charge,p=pol)
                xsec_ehi_totTheory = theoryBands[nameGraphTotTheory_xsec].GetErrorYhigh(iy)  
                xsec_elo_totTheory = theoryBands[nameGraphTotTheory_xsec].GetErrorYlow(iy)  
                nameGraphPDF_xsec = "hXsecPDF_{c}_{p}".format(c=charge,p=pol)
                xsec_ehi_pdf = theoryBands[nameGraphPDF_xsec].GetErrorYhigh(iy)  
                xsec_elo_pdf = theoryBands[nameGraphPDF_xsec].GetErrorYlow(iy)  
                # norm xsec
                nameGraphTotTheory_xsecnorm = "hXsecNorm{n}_{c}_{p}".format(n=hGrNameKey,c=charge,p=pol)
                xsecnorm_ehi_totTheory = theoryBands[nameGraphTotTheory_xsecnorm].GetErrorYhigh(iy)  
                xsecnorm_elo_totTheory = theoryBands[nameGraphTotTheory_xsecnorm].GetErrorYlow(iy)  
                nameGraphPDF_xsecnorm = "hXsecNormPDF_{c}_{p}".format(c=charge,p=pol)
                xsecnorm_ehi_pdf = theoryBands[nameGraphPDF_xsecnorm].GetErrorYhigh(iy)  
                xsecnorm_elo_pdf = theoryBands[nameGraphPDF_xsecnorm].GetErrorYlow(iy)  

                if options.normxsec:
                    tmp_val.val.append(xsecnorm_nominal[pol][iy]/ybinwidths[cp][iy])
                    tmp_val.ehi.append(xsecnorm_ehi_totTheory/ybinwidths[cp][iy])     
                    tmp_val.elo.append(xsecnorm_elo_totTheory/ybinwidths[cp][iy])
                    tmp_val.ehi2.append(xsecnorm_ehi_pdf/ybinwidths[cp][iy])
                    tmp_val.elo2.append(xsecnorm_elo_pdf/ybinwidths[cp][iy])
                    #rfit     = xsec_nominal[pol][iy]/normsigma/xsec_fit[0]
                    rfit     = xsecnorm_nominal[pol][iy]/xsec_fit[0]
                    tmp_val.rello.append(xsecnorm_elo_totTheory/xsecnorm_nominal[pol][iy])
                    tmp_val.relhi.append(xsecnorm_ehi_totTheory/xsecnorm_nominal[pol][iy])
                    tmp_val.rello2.append(xsecnorm_elo_pdf/xsecnorm_nominal[pol][iy])
                    tmp_val.relhi2.append(xsecnorm_ehi_pdf/xsecnorm_nominal[pol][iy])
                else:
                    tmp_val.val.append(xsec_nominal[pol][iy]/ybinwidths[cp][iy])
                    tmp_val.ehi.append(xsec_ehi_totTheory/ybinwidths[cp][iy])                
                    tmp_val.elo.append(xsec_elo_totTheory/ybinwidths[cp][iy])
                    tmp_val.ehi2.append(xsec_ehi_pdf/ybinwidths[cp][iy])
                    tmp_val.elo2.append(xsec_elo_pdf/ybinwidths[cp][iy])
                    rfit = xsec_nominal[pol][iy]/(xsec_fit[0]/scale) 
                    tmp_val.rello.append(xsec_elo_totTheory/xsec_nominal[pol][iy])
                    tmp_val.relhi.append(xsec_ehi_totTheory/xsec_nominal[pol][iy])
                    tmp_val.rello2.append(xsec_elo_pdf/xsec_nominal[pol][iy])
                    tmp_val.relhi2.append(xsec_ehi_pdf/xsec_nominal[pol][iy])
                
                tmp_val.val_fit.append(xsec_fit[0]/ybinwidths[cp][iy]/scale)
                tmp_val.elo_fit.append(abs(xsec_fit[0]-xsec_fit[1])/ybinwidths[cp][iy]/scale)
                tmp_val.ehi_fit.append(abs(xsec_fit[0]-xsec_fit[2])/ybinwidths[cp][iy]/scale)

                units = '' if options.normxsec else '(pb)'
                print "par = {parname}, expected sigma = {sigma:.3f} {units}   fitted = {val:.3f} + {ehi:.3f} - {elo:.3f} {units}".format(parname=parname,
                                                                                                                                          sigma=tmp_val.val[-1],units=units,
                                                                                                                                          val=tmp_val.val_fit[-1],ehi=tmp_val.ehi_fit[-1],elo=tmp_val.elo_fit[-1])

                tmp_val.relv. append(rfit);
                tmp_val.relv_fit .append(1.)
                #rfit_err = rfit*abs(xsec_fit[0]-xsec_fit[1])/xsec_fit[0]
                rfit_err = abs(xsec_fit[0]-xsec_fit[1])/xsec_fit[0] # if it is the error on the denominator of the ratio, which is the fit, then just take the relative uncertainty (rfit is exp/data)
                tmp_val.rello_fit.append(rfit_err)
                tmp_val.relhi_fit.append(rfit_err)

                tmp_val.rap.append((ybins[cp][iy]+ybins[cp][iy+1])/2.)
                tmp_val.rlo.append(abs(ybins[cp][iy]-tmp_val.rap[-1]))
                tmp_val.rhi.append(abs(ybins[cp][iy]-tmp_val.rap[-1]))

            tmp_val.makeGraphs()

            allValues[pol] = tmp_val

        plotValues(allValues,charge,channel,options, polarizations=signal_polarizations)

        # now do the unpolarized ones
        cp = 'plus_left' # this works if the binning for all the pol is the same

        ## first prepare the graphs to be used in the loop below
        # retrieve the graphs for unpolarized xsec
        nameGraphTotTheory_xsec = "hXsec{n}_{c}_unpolarized".format(n=hGrNameKey,c=charge)
        nameGraphTotTheory_xsecnorm = "hXsecNorm{n}_{c}_unpolarized".format(n=hGrNameKey,c=charge)
        nameGraphTotTheory_A4 = "hA4{n}_{c}".format(n=hGrNameKey,c=charge)
        nameGraphTotTheory_A0 = "hA0{n}_{c}".format(n=hGrNameKey,c=charge)
        nameGraphPDF_xsec = "hXsecPDF_{c}_unpolarized".format(c=charge)
        nameGraphPDF_xsecnorm = "hXsecNormPDF_{c}_unpolarized".format(c=charge)
        nameGraphPDF_A4 = "hA4PDF_{c}".format(c=charge)
        nameGraphPDF_A0 = "hA0PDF_{c}".format(c=charge)
        nameGraphs_TotTheory = {"sumxsec"     : nameGraphTotTheory_xsec,
                                "sumxsecnorm" : nameGraphTotTheory_xsecnorm,
                                "a0"          : nameGraphTotTheory_A0,
                                "a4"          : nameGraphTotTheory_A4
        }
        nameGraphs_PDF = {"sumxsec"     : nameGraphPDF_xsec,
                          "sumxsecnorm" : nameGraphPDF_xsecnorm,
                          "a0"          : nameGraphPDF_A0,
                          "a4"          : nameGraphPDF_A4
        }
        ## done with graphs
        
        xsec_params = ['sumxsecnorm','a0','a4'] if options.normxsec else ['sumxsec']
        for xs in xsec_params:
            tmp_val = valueClass('values_{xs}_{charge}_unpolarized'.format(xs=xs,charge=charge))
            print "total expected (fit) xsec up to |Y|<{maxy} = {sigma:.3f} pb".format(maxy=MAXYFORNORM,sigma=normsigmaIn)
            if len(outAccYBins):
                print "total expected (fit) xsec beyond |Y|>{maxy} = {sigma:.3f} pb".format(maxy=MAXYFORNORM,sigma=normsigmaOut)
            normsigma = normsigmaIn if abs(ybins[cp][iy])<MAXYFORNORM else normsigmaOut
            if doAltExp: alt_normsigma = alt_normsigmaIn if abs(ybins[cp][iy])<MAXYFORNORM else alt_normsigmaOut

            xs_nominal = []
            alt_xs_nominal = []

            for iy,y in enumerate(ybinwidths[cp]):
                # here will use the graphs for the bands, overriding the rest
                # will use band for totTheory, but it is the same
                xval = ROOT.Double(0)
                yval = ROOT.Double(0)  
                theoryBandsNoWpt[nameGraphs_TotTheory[xs]].GetPoint(iy,xval,yval)
                xs_nominal.append(yval)
                altyval = ROOT.Double(0)  
                theoryBands[nameGraphs_TotTheory[xs]].GetPoint(iy,xval,altyval)
                alt_xs_nominal.append(altyval)
            angcoeff_nominal[xs] = xs_nominal
            alt_angcoeff_nominal[xs] = alt_xs_nominal

            for iy,y in enumerate(ybinwidths[cp]):
                if iy in bkgYBins: continue
                parname = 'W{charge}_Ybin_{iy}_{xs}'.format(charge=charge,iy=iy,xs=xs)
                #print parname
                ybinwidth_scale = 1.
                scale = 1.
                if xs=='sumxsec':
                    ybinwidth_scale = ybinwidths[cp][iy]
                    scale = LUMINOSITY
                elif xs=='sumxsecnorm':
                    ybinwidth_scale = ybinwidths[cp][iy]
                #

                # can now use the actual key
                # note that the angular coefficients might come negative
                # although I think the fit returns positive because of 
                # the order of the L R 0 passed to polgroup.
                # Removing abs value here to see what the fit returns
                tmp_val.val.append(angcoeff_nominal[xs][iy]/ybinwidth_scale)
                if doAltExp:
                    tmp_val.altval.append(alt_angcoeff_nominal[xs][iy]/ybinwidth_scale)  
                experr_high = theoryBands[nameGraphs_TotTheory[xs]].GetErrorYhigh(iy)/ybinwidth_scale
                experr_low  = theoryBands[nameGraphs_TotTheory[xs]].GetErrorYlow(iy)/ybinwidth_scale
                tmp_val.ehi.append(experr_high)
                tmp_val.elo.append(experr_low)
                # pdfonly actually includes alpha, which is part of PDFs
                experr_pdfonly_high = theoryBands[nameGraphs_PDF[xs]].GetErrorYhigh(iy)/ybinwidth_scale
                experr_pdfonly_low  = theoryBands[nameGraphs_PDF[xs]].GetErrorYlow(iy)/ybinwidth_scale
                tmp_val.ehi2.append(experr_pdfonly_high)
                tmp_val.elo2.append(experr_pdfonly_low) 
                
                #print "NORMSIGMA = " + str(normsigma)
                
                if xs=='sumxsec':
                    scale *= float(nChan) 
                else:
                    pass # sumxsecnorm, a0, a4 are already normalized
        
                xsec_fit = valuesAndErrors[parname]
        
                val_fit_backup = xsec_fit[0]/ybinwidth_scale/scale
                tmp_val.val_fit.append(val_fit_backup)
                tmp_val.elo_fit.append(abs(xsec_fit[0]-xsec_fit[1])/ybinwidth_scale/scale)
                tmp_val.ehi_fit.append(abs(xsec_fit[0]-xsec_fit[2])/ybinwidth_scale/scale)

                if xs=='a4': # this is close to 0, bettter the difference wrt ratio
                    relv = tmp_val.val[-1] - val_fit_backup
                    experrrel_high = experr_high
                    experrrel_low = experr_low
                    experrrel_pdfonly_high = experr_pdfonly_high
                    experrrel_pdfonly_low  = experr_pdfonly_low
                else: 
                    relv = tmp_val.val[-1]/val_fit_backup
                    experrrel_high = experr_high/val_fit_backup
                    experrrel_low = experr_low/val_fit_backup
                    experrrel_pdfonly_high = experr_pdfonly_high/val_fit_backup
                    experrrel_pdfonly_low = experr_pdfonly_low/val_fit_backup

                tmp_val.relv. append(relv)
                tmp_val.rello.append(experrrel_low)
                tmp_val.relhi.append(experrrel_high) 
                tmp_val.rello2.append(experrrel_pdfonly_low)
                tmp_val.relhi2.append(experrrel_pdfonly_high)
                
                units = '(pb)' if xs=='sumxsec' else ''
                print "par = {parname}, expected value = {sigma:.3f} {units}   fitted = {val:.3f} + {ehi:.3f} - {elo:.3f} {units}".format(parname=parname, sigma=tmp_val.val[-1],units=units,val=tmp_val.val_fit[-1],ehi=tmp_val.ehi_fit[-1],elo=tmp_val.elo_fit[-1])
                if xs=='a4':
                    relv_fit = 0.
                    rello_fit = tmp_val.elo_fit[-1]
                    relhi_fit = tmp_val.ehi_fit[-1]
                else:
                    relv_fit = 1.
                    rello_fit = tmp_val.elo_fit[-1]/tmp_val.val_fit[-1]
                    relhi_fit = tmp_val.ehi_fit[-1]/tmp_val.val_fit[-1]
                tmp_val.relv_fit .append(relv_fit)
                tmp_val.rello_fit.append(rello_fit)
                tmp_val.relhi_fit.append(relhi_fit)
        
                tmp_val.rap.append((ybins[cp][iy]+ybins[cp][iy+1])/2.)
                tmp_val.rlo.append(abs(ybins[cp][iy]-tmp_val.rap[-1]))
                tmp_val.rhi.append(abs(ybins[cp][iy]-tmp_val.rap[-1]))
        
            tmp_val.makeGraphs()
            plotUnpolarizedValues(tmp_val,charge,channel,options)

                
    if len(charges)>1:
        print "Making charge asymmetry plots now..."
        asymmetryValues = {}
        
        for pol in signal_polarizations:
            cp = 'plus_'+pol
            tmp_val = valueClass('asymmetry_'+pol)

            # get graph for asymmetry to retrieve band
            # use the one with total theory (the QCD scales are basically cancelled)
            nameGraphs_TotTheory = "hAsym{n}_{p}".format(n=hGrNameKey,p=pol)
            nameGraphs_PDF = "hAsymPDF_{p}".format(p=pol)
        
            for iy,y in enumerate(ybinwidths[cp]):
                if iy in bkgYBins: continue

                xval = ROOT.Double(0)
                yval = ROOT.Double(0)                  
                altyval = ROOT.Double(0)                  
                theoryBandsNoWpt[nameGraphs_TotTheory].GetPoint(iy,xval,yval)
                theoryBands[nameGraphs_TotTheory].GetPoint(iy,xval,altyval)
                tmp_val.val .append(yval)
                if doAltExp:
                    tmp_val.altval .append(altyval)
                ehi = theoryBands[nameGraphs_TotTheory].GetErrorYhigh(iy)
                elo = theoryBands[nameGraphs_TotTheory].GetErrorYlow(iy)
                ehi_pdf = theoryBands[nameGraphs_PDF].GetErrorYhigh(iy)
                elo_pdf = theoryBands[nameGraphs_PDF].GetErrorYlow(iy)
                tmp_val.ehi.append(ehi)
                tmp_val.elo.append(elo)
                tmp_val.ehi2.append(ehi_pdf)
                tmp_val.elo2.append(elo_pdf)

                if options.type == 'toys':
                    asy_fit = utilities.getAsymmetryFromToys(pol,channel,iy,options.infile)
                else:
                    asy_fit = valuesAndErrors['W_{pol}_Ybin_{iy}_chargeasym'.format(pol=pol,iy=iy)]
                tmp_val.val_fit .append(asy_fit[0])
                tmp_val.elo_fit.append(abs(asy_fit[0]-asy_fit[1]))
                tmp_val.ehi_fit.append(abs(asy_fit[0]-asy_fit[2]))

                # on the charge asymmetry, which is A~0, better to show the difference 
                # Afit - Aexp wrt the ratio. The error on "exp" diff shows the error on Aexp, while the error bar the error on Afit
                tmp_val.relv. append(tmp_val.val[-1] - tmp_val.val_fit[-1])
                tmp_val.rello.append(elo)
                tmp_val.relhi.append(ehi)
                tmp_val.rello2.append(elo_pdf)
                tmp_val.relhi2.append(ehi_pdf)

                tmp_val.relv_fit .append(0.)
                tmp_val.rello_fit.append(tmp_val.elo_fit[-1])
                tmp_val.relhi_fit.append(tmp_val.ehi_fit[-1])

                tmp_val.rap.append((ybins[cp][iy]+ybins[cp][iy+1])/2.)
                tmp_val.rlo.append(abs(ybins[cp][iy]-tmp_val.rap[-1]))
                tmp_val.rhi.append(abs(ybins[cp][iy]-tmp_val.rap[-1]))

            tmp_val.makeGraphs()
            asymmetryValues[pol] = tmp_val        
        plotValues(asymmetryValues,'asymmetry',channel,options, polarizations=signal_polarizations)
            
        
        # get graph for asymmetry to retrieve band
        # use the one with total theory (the QCD scales are basically cancelled)
        nameGraphs_TotTheory = "hAsym{n}_unpolarized".format(n=hGrNameKey)
        nameGraphs_PDF = "hAsymPDF_unpolarized"

        # now do the unpolarized ones
        tmp_val = valueClass('values_asymmetry_unpolarized')
        for iy,y in enumerate(ybinwidths['plus_left']): # this assumes that all the 3 polarizations have the same binning

            if iy in bkgYBins: continue            
            xval = ROOT.Double(0)
            yval = ROOT.Double(0)                  
            altyval = ROOT.Double(0)                  
            theoryBandsNoWpt[nameGraphs_TotTheory].GetPoint(iy,xval,yval)
            theoryBands[nameGraphs_TotTheory].GetPoint(iy,xval,altyval)
            tmp_val.val .append(yval)
            if doAltExp:
                tmp_val.altval .append(altyval)
            ehi = theoryBands[nameGraphs_TotTheory].GetErrorYhigh(iy)
            elo = theoryBands[nameGraphs_TotTheory].GetErrorYlow(iy)
            ehi_pdf = theoryBands[nameGraphs_PDF].GetErrorYhigh(iy)
            elo_pdf = theoryBands[nameGraphs_PDF].GetErrorYlow(iy)
            tmp_val.ehi.append(ehi)
            tmp_val.elo.append(elo)
            tmp_val.ehi2.append(ehi_pdf)
            tmp_val.elo2.append(elo_pdf)

            if options.type == 'hessian': 
                asy_fit = valuesAndErrors['W_Ybin_{iy}_chargemetaasym'.format(iy=iy)]
            tmp_val.val_fit .append(asy_fit[0])
            tmp_val.elo_fit.append(abs(asy_fit[0]-asy_fit[1]))
            tmp_val.ehi_fit.append(abs(asy_fit[0]-asy_fit[2]))

            tmp_val.relv. append(tmp_val.val[-1] - tmp_val.val_fit[-1])
            tmp_val.rello.append(elo)
            tmp_val.relhi.append(ehi)
            tmp_val.rello2.append(elo_pdf)
            tmp_val.relhi2.append(ehi_pdf)

            tmp_val.relv_fit .append(0.)
            tmp_val.rello_fit.append(tmp_val.elo_fit[-1])
            tmp_val.relhi_fit.append(tmp_val.ehi_fit[-1])

            tmp_val.rap.append((ybins[cp][iy]+ybins[cp][iy+1])/2.)
            tmp_val.rlo.append(abs(ybins[cp][iy]-tmp_val.rap[-1]))
            tmp_val.rhi.append(abs(ybins[cp][iy]-tmp_val.rap[-1]))

        tmp_val.makeGraphs()
        #print "CHECK"
        plotUnpolarizedValues(tmp_val,'asymmetry',channel,options)
        #print "END CHECK"
            
