#!/bin/env python
# USAGE: python plotYW.py --type toys --infile toys_wplus.root -y cards_el/binningYW.txt -C plus -o plots --xsecfiles Wel_plus_shapes_xsec.root [--normxsec]
#
# When run on the W+ W- combined fit it also makes charge asymmetry plot:
# python plotYW.py --type toys --infile toys_wboth.root -y cards_el/binningYW.txt -C plus,minus -o plots --xsecfiles Wel_plus_shapes_xsec.root,Wel_minus_shapes_xsec.root

import ROOT, datetime, array, os, math, re
ROOT.gROOT.SetBatch(True)
ROOT.PyConfig.IgnoreCommandLineOptions = True
ROOT.gSystem.Load("libHiggsAnalysisCombinedLimit")
from plotYW import valueClass

import utilities
utilities = utilities.util()

REFMC = 'MC@NLO'

def makeCanvas():
    c2 = ROOT.TCanvas('foo','', 2000, 3000)
    
    c2.Range(0,0,1,1);
    c2.SetFillColor(0);
    c2.SetBorderMode(0);
    c2.SetBorderSize(2);
    c2.SetTickx(1);
    c2.SetTicky(1);
    c2.SetLeftMargin(0.05);
    c2.SetRightMargin(0.05);
    c2.SetTopMargin(0.05);
    c2.SetBottomMargin(0.20);
    c2.SetFrameFillStyle(0);
    c2.SetFrameBorderMode(0);
    return c2

def makePads(subpadsYLow,subpadsYUp):
    pads = []
    for ip in xrange(len(subpadsYLow)):
        pad = ROOT.TPad("pad{ip}".format(ip=ip),"",0.,subpadsYLow[ip],1.,subpadsYUp[ip])
        pad.SetFillColor(0);
        pad.SetBorderMode(0);
        pad.SetBorderSize(2);
        pad.SetTickx(1);
        pad.SetTicky(1);
        pad.SetLeftMargin(0.2);
        pad.SetRightMargin(0.05);
        pad.SetTopMargin(0.05);
        if ip==len(subpadsYLow)-1:
            pad.SetBottomMargin(0.15);
        else:
            pad.SetBottomMargin(0.05);
        pad.SetFrameFillStyle(0);
        pad.SetFrameBorderMode(0);
        pad.SetFrameFillStyle(0);
        pad.SetFrameBorderMode(0);
        pads.append(pad)
    return pads

def removeErrX(graph):
    for p in xrange(graph.GetN()):
        x = ROOT.Double(0); y = ROOT.Double(0)
        graph.GetPoint(p,x,y)
        graph.SetPointEXhigh(p,0); graph.SetPointEXlow(p,0)

def applyChannelStyles(values,polarizations):
    markerStyles = {'mu': ROOT.kOpenSquare, 'el': ROOT.kOpenCircle, 'lep': ROOT.kFullTriangleDown}
    markerSizes = {'mu': 2, 'el': 2, 'lep': 2}
    for flav in ['lep','mu','el']:
        for pol in polarizations:
            # remove the x error bars for the single channels
            if flav!='lep': removeErrX(values[(flav,pol)].graph_fit)
            values[(flav,pol)].graph.SetLineStyle(1)
            values[(flav,pol)].graph.SetLineWidth(2)
            values[(flav,pol)].graph.SetLineColor(values[(flav,pol)].color)
            values[(flav,pol)].graph_fit.SetMarkerStyle(markerStyles[flav])
            values[(flav,pol)].graph_fit.SetMarkerSize(markerSizes[flav])
            values[(flav,pol)].graph_fit_rel.SetMarkerStyle(markerStyles[flav])
            values[(flav,pol)].graph_fit_rel.SetMarkerSize(markerSizes[flav])
            values[(flav,pol)].graph_fit_rel.SetLineColor(ROOT.kBlack)
            values[(flav,pol)].graph_fit_rel.SetLineWidth(2)
            values[(flav,pol)].graph_fit_rel.SetFillColor(ROOT.kYellow-9)
            values[(flav,pol)].graph_fit_rel.SetFillStyle(1001)
            for p in xrange(values[(flav,pol)].graph_fit_rel.GetN()):
                values[(flav,pol)].graph_fit_rel.SetPointEXhigh(p,values[(flav,pol)].rhi[p])
                values[(flav,pol)].graph_fit_rel.SetPointEXlow(p,values[(flav,pol)].rlo[p])

def applyChannelStylesUnpol(values):
    markerStyles = {'mu': ROOT.kOpenSquare, 'el': ROOT.kOpenCircle, 'lep': ROOT.kFullTriangleDown}
    markerSizes = {'mu': 2, 'el': 2, 'lep': 2}
    tmpcolors = {'mu': ROOT.kAzure+7, 'el': ROOT.kRed-2, 'lep': ROOT.kGreen-2}
    for flav in ['lep','mu','el']:
        # remove the x error bars for the single channels
        values[flav].graph.SetLineStyle(1)
        values[flav].graph.SetLineWidth(2)
        values[flav].graph.SetLineColor(values[flav].color)
        values[flav].graph_fit.SetMarkerStyle(markerStyles[flav])
        values[flav].graph_fit.SetMarkerSize(markerSizes[flav])
        values[flav].graph_fit.SetMarkerColor(tmpcolors[flav])
        values[flav].graph_fit.SetLineColor(tmpcolors[flav])
        values[flav].graph_fit_rel.SetMarkerStyle(markerStyles[flav])
        values[flav].graph_fit_rel.SetMarkerSize(markerSizes[flav])
        values[flav].graph_fit_rel.SetLineColor(ROOT.kBlack)
        values[flav].graph_fit_rel.SetLineWidth(2)
        values[flav].graph_fit_rel.SetFillColor(ROOT.kYellow-9)
        values[flav].graph_fit_rel.SetFillStyle(1001)
        for p in xrange(values[flav].graph_fit_rel.GetN()):
            values[flav].graph_fit_rel.SetPointEXhigh(p,values[flav].rhi[p])
            values[flav].graph_fit_rel.SetPointEXlow(p,values[flav].rlo[p])

def makeUncorrRatio(gr1,gr2):
    ratio = ROOT.TGraphAsymmErrors(gr1.GetN())
    ratio.SetName('graph_ratio')
    for p in xrange(gr1.GetN()):
        x = ROOT.Double(0); y1 = ROOT.Double(0); y2 = ROOT.Double(0)
        gr1.GetPoint(p,x,y1); gr2.GetPoint(p,x,y2)
        r = y2/y1
        # consider el and mu mostly uncorrelated (since dominated by the MC stat)
        erhi = r*math.hypot(gr1.GetErrorYhigh(p)/y1,gr2.GetErrorYhigh(p)/y2); erlo=erhi
        ratio.SetPoint(p,x,r)
        ratio.SetPointError(p,gr1.GetErrorXlow(p),gr1.GetErrorXhigh(p),erlo,erhi)
    return ratio

def makeUncorrDiff(gr1,gr2):
    diff = ROOT.TGraphAsymmErrors(gr1.GetN())
    diff.SetName('graph_diff')
    for p in xrange(gr1.GetN()):
        x = ROOT.Double(0); y1 = ROOT.Double(0); y2 = ROOT.Double(0);
        gr2.GetPoint(p,x,y2); gr1.GetPoint(p,x,y1)
        d = y2-y1; 
        # consider el and mu mostly uncorrelated (since dominated by the MC stat)
        erhi = math.hypot(gr1.GetErrorYhigh(p),gr2.GetErrorYhigh(p)); erlo=erhi
        diff.SetPoint(p,x,d)
        diff.SetPointError(p,gr1.GetErrorXlow(p),gr1.GetErrorXhigh(p),erlo,erhi)
    return diff

def makeFullLegend(values):
    doAltExp = (len(values[('lep','left')].altval)>0)
    legs = []
    # expected values
    legExp = ROOT.TLegend(0.25, 0.75, 0.90, 0.85)
    legExp.SetFillStyle(0)
    legExp.SetBorderSize(0)
    legExp.AddEntry(values[('lep','left')] .graph     , 'W_{{L}} ({mc})'.format(mc=REFMC), 'f')
    legExp.AddEntry(values[('lep','right')].graph     , 'W_{{R}} ({mc})'.format(mc=REFMC), 'f')
    if doAltExp:
        legExp.SetNColumns(2)
        legExp.AddEntry(values[('lep','left')] .altgraph     , 'W_{{L}} ({mc}*)'.format(mc=REFMC) , 'l')
        legExp.AddEntry(values[('lep','right')].altgraph     , 'W_{{R}} ({mc}*)'.format(mc=REFMC) , 'l')
    legs.append(legExp)
    # data Left
    legLeft = ROOT.TLegend(0.25, 0.60, 0.90, 0.75)
    legLeft.SetFillStyle(0)
    legLeft.SetBorderSize(0)
    legLeft.AddEntry(values[('mu','left')] .graph_fit     , ' ', 'p')
    legLeft.AddEntry(values[('el','left')].graph_fit     , ' ', 'p')
    legLeft.AddEntry(values[('lep','left')].graph_fit     , ' ', 'p')
    legs.append(legLeft)
    # data Right
    legRight = ROOT.TLegend(0.35, 0.60, 0.90, 0.75)
    legRight.SetFillStyle(0)
    legRight.SetBorderSize(0)
    legRight.AddEntry(values[('mu','right')] .graph_fit    , 'data #mu' , 'p')
    legRight.AddEntry(values[('el','right')].graph_fit     , 'data e'   , 'p')
    legRight.AddEntry(values[('lep','right')].graph_fit    , 'data comb', 'p')
    legs.append(legRight)
    return legs

def makeFullLegendUnpol(values):
    doAltExp = (len(values['lep'].altval)>0)
    # expected values
    leg = ROOT.TLegend(0.25, 0.60, 0.50, 0.85)
    leg.SetFillStyle(0)
    leg.SetBorderSize(0)
    leg.AddEntry(values['lep'] .graph    , REFMC, 'f')
    if doAltExp:
        leg.AddEntry(values['lep'] .altgraph    , REFMC+'*', 'l')
    leg.AddEntry(values['mu'] .graph_fit , 'data #mu', 'pl')
    leg.AddEntry(values['el'] .graph_fit , 'data e', 'pl')
    leg.AddEntry(values['lep'].graph_fit , 'data comb ', 'pl')
    return leg

def plotValues(values,charge,channel,options):
    
    c2 = makeCanvas()

    ROOT.gStyle.SetHatchesLineWidth(2)
    subpadsYLow = [0.53, 0.40, 0.27, 0.16, 0.06]
    subpadsYUp =  [0.99, 0.50, 0.37, 0.27, 0.16]
    pads = makePads(subpadsYLow,subpadsYUp)

    skipLong = False
    if options.nolong or options.longBkg: skipLong = True
    doAltExp = (len(values[('lep','left')].altval)>0)

    flavors = ['mu','el','lep']; 
    polarizations = ['left','right']
    if not skipLong: polarizations.append('long')
    ch = '#plus' if charge == 'plus' else '#minus'
    if charge == 'asymmetry': ch = ''
    date = datetime.date.today().isoformat()
    normstr = 'norm' if (options.normxsec and charge!='asymmetry') else ''

    applyChannelStyles(values,polarizations)

    lat = ROOT.TLatex()
    lat.SetNDC(); lat.SetTextFont(42)
    ## the four graphs exist now. now starting to draw them
    ## ===========================================================
    if sum(hasattr(values[(flav,pol)],'graph') and hasattr(values[(flav,pol)],'graph_fit') for pol in polarizations for flav in flavors)==3*len(polarizations):
        c2.cd()
        pads[0].Draw(); pads[0].cd(); ROOT.gPad.SetBottomMargin(0);

        values[('lep','left')].graph.SetTitle('W {ch}: y_{{W}}'.format(ch=ch))
            
        mg = ROOT.TMultiGraph()
        mg.Add(values[('lep','left')] .graph,'P2')
        mg.Add(values[('lep','right')].graph,'P2')
        if not skipLong: mg.Add(values['lep','long'] .graph,'P2')
        for flav in flavors:
            mg.Add(values[(flav,'left')] .graph_fit)
            mg.Add(values[(flav,'right')].graph_fit)
            if not skipLong: mg.Add(values[(flav,'long')] .graph_fit)
            if doAltExp:
                mg.Add(values[(flav,'left')] .altgraph,'L2')
                mg.Add(values[(flav,'right')].altgraph,'L2')
    
        mg.Draw('Pa')
        mg.GetXaxis().SetRangeUser(0., options.maxRapidity) # max would be 6.
        mg.GetXaxis().SetTitle('|y_{W}|')
        mg.GetXaxis().SetLabelSize(0)
        if charge=='asymmetry':
            mg.GetYaxis().SetTitle('Charge asymmetry')
            mg.GetYaxis().SetRangeUser(-0.1,0.4)
        else:
            if options.normxsec: 
                mg.GetYaxis().SetTitle('d#sigma / #sigma_{tot}^{fit} / d|y_{W}|')
                mg.GetYaxis().SetRangeUser(-0.05,0.8 if options.maxRapidity > 2.9 else 0.4)
            else: 
                mg.GetYaxis().SetTitle('d#sigma (pb) / d|y_{W}|')
                mg.GetYaxis().SetRangeUser(-200,3500)
        mg.GetYaxis().SetTitleSize(0.05)
        mg.GetYaxis().SetLabelSize(0.04)
        mg.GetYaxis().SetTitleOffset(1.5)

        legs = makeFullLegend(values)     
        for leg in legs:
            leg.Draw("same")

    # save for HEPData
    plotname = '{od}/genAbsY{norm}_pdfs_{ch}{suffix}'.format(od=options.outdir, norm=normstr, ch=charge, suffix=options.suffix)
    rfHepName = '{n}.root'.format(n=plotname)
    rfHep = ROOT.TFile.Open(rfHepName,"recreate")
    if not rfHep:
        print "Error in plotYW.py: could not open root file %s" % rfHepName
        quit()
    rfHep.cd()
    polsToSave = ["left","right"] + ([] if skipLong else ["long"])
    for pts in polsToSave:
        for flav in flavors:
            values[(flav,pts)].graph_fit.Clone().Write("data_{f}_{p}".format(f=flav,p=pts))
            if flav == "lep":
                values[(flav,pts)].graph.Clone().Write("exp_{f}_{p}".format(f=flav,p=pts))
                if doAltExp:
                    values[(flav,pts)].altgraph.Clone().Write("expAlt_{f}_{p}".format(f=flav,p=pts))
    rfHep.Close()
    # done with hepdata

    lines = {}
    subpadLegends = {}
    if sum(hasattr(values[(flav,pol)],'mg') for pol in polarizations for flav in flavors)==3*len(polarizations):

        ## now make the relative error plot for combination:
        ## ======================================
        isubpad=1
        for ipol,pol in enumerate(['left','right']):
            c2.cd()
            pads[isubpad].Draw(); pads[isubpad].cd(); ROOT.gPad.SetBottomMargin(0); ROOT.gPad.SetTopMargin(0)

            if charge=='asymmetry':
                yaxtitle = 'Theory-Data'
                yaxrange = (-0.2, 0.2)
            else:
                yaxtitle = 'Theory/Data'
                yaxrange = (0.82, 1.18)
    
            values[('lep',pol)].graph_fit_rel.Draw("A2")
            values[('lep',pol)].graph_fit_rel.Draw("pe")
            values[('lep',pol)].graph_rel.Draw("2")
            ## x axis fiddling
            values[('lep',pol)].graph_fit_rel.GetXaxis().SetRangeUser(0., options.maxRapidity)
            values[('lep',pol)].graph_fit_rel.GetXaxis().SetLabelSize(0.04)
            values[('lep',pol)].graph_fit_rel.GetXaxis().SetLabelSize(0.0)
            values[('lep',pol)].graph_fit_rel.GetXaxis().SetTitleSize(0.0)
            ## y axis fiddling
            values[('lep',pol)].graph_fit_rel.GetYaxis().SetTitleOffset(0.4)
            values[('lep',pol)].graph_fit_rel.GetYaxis().SetTitleSize(0.18)
            values[('lep',pol)].graph_fit_rel.GetYaxis().SetLabelSize(0.12)
            values[('lep',pol)].graph_fit_rel.GetYaxis().SetTitle(yaxtitle)
            values[('lep',pol)].graph_fit_rel.GetYaxis().SetRangeUser(yaxrange[0],yaxrange[1])
            values[('lep',pol)].graph_fit_rel.GetYaxis().CenterTitle()
            values[('lep',pol)].graph_fit_rel.GetYaxis().SetDecimals()

            lines["horiz_line_comb"+pol] = ROOT.TF1("horiz_line_comb"+pol,"0" if charge=='asymmetry' else '1',0.0,3.0);
            lines["horiz_line_comb"+pol].SetLineColor(ROOT.kRed);
            lines["horiz_line_comb"+pol].SetLineWidth(2);
            lines["horiz_line_comb"+pol].SetLineStyle(ROOT.kDashed);
            lines["horiz_line_comb"+pol].Draw("Lsame");
            
            subpadLegends["leg_lep_"+pol] = ROOT.TLegend(0.25, 0.8, 0.5, 0.9); subpadLegends["leg_lep_"+pol].SetFillStyle(0); subpadLegends["leg_lep_"+pol].SetBorderSize(0)
            subpadLegends["leg_lep_"+pol].SetNColumns(2)
            subpadLegends["leg_lep_"+pol].AddEntry(values[('lep',pol)].graph_fit_rel,'Data','fpl')
            labels_pol = {'left':'L','right':'R','long':'0'} 
            subpadLegends["leg_lep_"+pol].AddEntry(values[('lep',pol)].graph_rel,'Theory (W_{{{pol}}})'.format(pol=labels_pol[pol]),'f')
            subpadLegends["leg_lep_"+pol].Draw('same')

            isubpad+=1

        # now the compatibility plots
        comps = {}
        for ipol,pol in enumerate(['left','right']):
            c2.cd()
            pads[isubpad].Draw(); pads[isubpad].cd(); ROOT.gPad.SetTopMargin(0)
            if isubpad<4: ROOT.gPad.SetBottomMargin(0)
            else: ROOT.gPad.SetBottomMargin(-0.3)

            if charge=='asymmetry':
                yaxtitle = 'A_{#mu}-A_{e}'
                yaxrange = (-0.3, 0.3)
                comps["comp_"+pol] = makeUncorrDiff(values[('el',pol)].graph_fit, values[('mu',pol)].graph_fit)
            else:
                yaxtitle = '#sigma_{#mu}/#sigma_{e}'
                yaxrange = (0.65, 1.35)
                comps["comp_"+pol] = makeUncorrRatio(values[('el',pol)].graph_fit, values[('mu',pol)].graph_fit)
            comps["comp_"+pol].SetName("comp_"+pol); comps["comp_"+pol].SetTitle("")
            comps["comp_"+pol].SetMarkerColor(values[('mu',pol)].color)
            comps["comp_"+pol].SetLineColor(values[('mu',pol)].color)
            comps["comp_"+pol].SetFillColor(values[('mu',pol)].color)
            comps["comp_"+pol].SetMarkerStyle(ROOT.kFullCircle)
            comps["comp_"+pol].SetMarkerSize(2)

            comps["comp_"+pol].Draw('Pa')
            ## x axis fiddling
            comps["comp_"+pol].GetXaxis().SetRangeUser(0., options.maxRapidity)
            comps["comp_"+pol].GetXaxis().SetLabelSize(0.1)
            ## y axis fiddling
            comps["comp_"+pol].GetYaxis().SetTitleOffset(0.4)
            comps["comp_"+pol].GetYaxis().SetTitleSize(0.18)
            comps["comp_"+pol].GetYaxis().SetLabelSize(0.12)
            comps["comp_"+pol].GetYaxis().SetTitle(yaxtitle)
            comps["comp_"+pol].GetYaxis().SetRangeUser(yaxrange[0],yaxrange[1])
            comps["comp_"+pol].GetYaxis().CenterTitle()
            comps["comp_"+pol].GetYaxis().SetDecimals()
            lines["horiz_line_comp_"+pol] = ROOT.TF1("horiz_line_comp_"+pol,"0" if charge=='asymmetry' else '1',0.0,3.0);
            lines["horiz_line_comp_"+pol].SetLineColor(ROOT.kRed);
            lines["horiz_line_comp_"+pol].SetLineWidth(2);
            lines["horiz_line_comp_"+pol].SetLineStyle(ROOT.kDashed);
            lines["horiz_line_comp_"+pol].Draw("Lsame");

            subpadLegends["comp_lep_"+pol] = ROOT.TLegend(0.25, 0.8, 0.5, 0.9); subpadLegends["comp_lep_"+pol].SetFillStyle(0); subpadLegends["comp_lep_"+pol].SetBorderSize(0)
            subpadLegends["comp_lep_"+pol].SetNColumns(2)
            labels_pol = {'left':'L','right':'R','long':'0'} 
            legtitle = "Data #mu - e" if charge=='asymmetry' else "Data #mu / e" 
            subpadLegends["comp_lep_"+pol].AddEntry(comps["comp_"+pol],'{legtitle} (W_{{{pol}}})'.format(legtitle=legtitle,pol=labels_pol[pol]),'pl')
            subpadLegends["comp_lep_"+pol].Draw('same')

            isubpad+=1

    c2.cd()
    # was 0.95 but for me it appears inside the frame
    lat.DrawLatex(0.2, 0.97, '#bf{CMS}') #it{Preliminary}')
    lat.DrawLatex(0.62, 0.97, '35.9 fb^{-1} (13 TeV)')
    lat.DrawLatex(0.25, 0.60,  'W^{{{ch}}} #rightarrow l^{{{ch}}}{nu}'.format(ch=ch,nu="#bar{#nu}" if charge=='minus' else "#nu"))
    lat.DrawLatex(0.85, 0.025, '|y_{W}|')
    for ext in ['png', 'pdf']:
        c2.SaveAs('{p}.{ext}'.format(p=plotname, ext=ext))


def plotUnpolarizedValues(values,charge,channel,options):

    c2 = makeCanvas()

    ROOT.gStyle.SetHatchesLineWidth(2)

    subpadsYLow = [0.52, 0.30, 0.08]
    subpadsYUp =  [0.98, 0.50, 0.28]
    pads = makePads(subpadsYLow,subpadsYUp)

    flavors = ['mu','el','lep'];
    polarizations = ['left','right']
    skipLong = False
    if options.nolong or options.longBkg: skipLong = True
    if not skipLong: 
        polarizations.append('long')
    ch = '#plus' if charge == 'plus' else '#minus'
    if charge == 'asymmetry': ch = ''
    date = datetime.date.today().isoformat()
    doAltExp = (len(values['lep'].altval)>0)

    valkey = values['lep'].name.split('_')[1]

    applyChannelStylesUnpol(values)

    lat = ROOT.TLatex()
    lat.SetNDC(); lat.SetTextFont(42)
    legx1, legx2, legy1, legy2 = 0.2, 0.5, 0.7, 0.85
    ## the graphs exist now. now starting to draw them
    ## ===========================================================
    if sum(hasattr(values[flav],'graph') and hasattr(values[flav],'graph_fit') for flav in flavors)==3:
            
        c2.cd()
        pads[0].Draw(); pads[0].cd(); ROOT.gPad.SetBottomMargin(0);

        mg = ROOT.TMultiGraph()
        mg.Add(values['lep'].graph,'P2')
        for flav in flavors:
            mg.Add(values[flav].graph_fit)
            if doAltExp:
                mg.Add(values[flav].altgraph,'L2')
        
        mg.Draw('Pa')
        mg.GetXaxis().SetRangeUser(0., options.maxRapidity) # max would be 6.
        mg.GetXaxis().SetTitle('|y_{W}|')
        mg.GetXaxis().SetTitleOffset(5.5)
        mg.GetXaxis().SetLabelSize(0)
        titles = {'asymmetry': 'Charge asymmetry',
                  'a0': 'A_{0}', 
                  'a4': 'A_{4}',
                  'sumxsec': 'd#sigma / d|y_{W}| (pb)',
                  'sumxsecnorm': 'd#sigma / #sigma_{tot}^{fit} / d|y_{W}|'}
        ranges = {'asymmetry': (-0.1,0.4),
                  'a0': (0.07,0.2),
                  'a4': (-1,2),
                  'sumxsec': (1700,4500),
                  'sumxsecnorm': (0.1,0.5)}
        mg.GetYaxis().SetTitle(titles[valkey])
        mg.GetYaxis().SetRangeUser(ranges[valkey][0],ranges[valkey][1])

        mg.GetYaxis().SetTitleSize(0.06)
        mg.GetYaxis().SetLabelSize(0.04)
        mg.GetYaxis().SetTitleOffset(1.4)
        mg.GetYaxis().CenterTitle()

        leg = makeFullLegendUnpol(values)
        leg.Draw('same')
    
    # save for HEPData
    plotname = '{od}/genAbsYUnpolarized{norm}_pdfs_{ch}{suffix}'.format(od=options.outdir, norm=valkey, ch=charge, suffix=options.suffix)
    rfHepName = '{n}.root'.format(n=plotname)
    rfHep = ROOT.TFile.Open(rfHepName,"recreate")
    if not rfHep:
        print "Error in plotYW.py: could not open root file %s" % rfHepName
        quit()
    rfHep.cd()
    for flav in flavors:
        values[flav].graph_fit.Clone().Write("data_{f}".format(f=flav))
        if flav == "lep":
            values[flav].graph.Clone().Write("exp_{f}".format(f=flav))
            if doAltExp:
                values[flav].altgraph.Clone().Write("expAlt_{f}".format(f=flav))
    rfHep.Close()
    # done with hepdata

    lines = {}
    subpadLegends = {}
    ## now make the relative error plot:
    ## ======================================
    if sum(hasattr(values[flav],'mg') for flav in flavors)==3:

        c2.cd()
        pads[1].Draw(); pads[1].cd(); # ROOT.gPad.SetBottomMargin(0); ROOT.gPad.SetTopMargin(0)

        if valkey=='asymmetry':
            yaxtitle = 'Theory-Data'
            yaxrange = (-0.2, 0.2)
            yaxtitlesize = 0.15
        else:
            yaxtitle = 'Theory/Data'
            yaxrange = (0.75, 1.25) if 'xsec' in valkey else (0.7,1.3) if valkey=='a0' else (-0.5, 2.5)
            yaxtitlesize = 0.15

        values['lep'].graph_fit_rel.Draw("A2")
        values['lep'].graph_fit_rel.Draw("pe")
        values['lep'].graph_rel.Draw("2")

        ## x axis fiddling
        values['lep'].graph_fit_rel.GetXaxis().SetRangeUser(0., options.maxRapidity)
        values['lep'].graph_fit_rel.GetXaxis().SetLabelSize(0)
        values['lep'].graph_fit_rel.GetXaxis().SetTitleSize(0)
        ## y axis fiddling
        values['lep'].graph_fit_rel.GetYaxis().SetTitleOffset(0.45)
        values['lep'].graph_fit_rel.GetYaxis().SetTitleSize(yaxtitlesize)
        values['lep'].graph_fit_rel.GetYaxis().SetLabelSize(0.12)
        values['lep'].graph_fit_rel.GetYaxis().SetTitle(yaxtitle)
        values['lep'].graph_fit_rel.GetYaxis().SetRangeUser(yaxrange[0],yaxrange[1])
        values['lep'].graph_fit_rel.GetYaxis().CenterTitle()
        values['lep'].graph_fit_rel.GetYaxis().SetDecimals()
        values['lep'].graph_fit_rel.GetYaxis().SetNdivisions(505,ROOT.kFALSE)

        lines["horiz_line_comb"] = ROOT.TF1("horiz_line_comb","0" if charge=='asymmetry' else '1',0.0,3.0);
        lines["horiz_line_comb"].SetLineColor(ROOT.kRed);
        lines["horiz_line_comb"].SetLineWidth(2);
        lines["horiz_line_comb"].SetLineStyle(ROOT.kDashed);
        lines["horiz_line_comb"].Draw("Lsame");
        
        subpadLegends["leg_lep"] = ROOT.TLegend(0.25, 0.8, 0.5, 0.9); subpadLegends["leg_lep"].SetFillStyle(0); subpadLegends["leg_lep"].SetBorderSize(0)
        subpadLegends["leg_lep"].SetNColumns(2)
        subpadLegends["leg_lep"].AddEntry(values['lep'].graph_fit_rel,'Data','fpl')
        subpadLegends["leg_lep"].AddEntry(values['lep'].graph_rel,'Theory','f')
        subpadLegends["leg_lep"].Draw('same')

        # now the compatibility plots
        c2.cd()
        pads[2].Draw(); pads[2].cd(); #ROOT.gPad.SetTopMargin(0)        ROOT.gPad.SetBottomMargin(-0.3)

        if charge=='asymmetry':
            yaxtitle = 'A_{#mu}-A_{e}'
            yaxrange = (-0.3, 0.3)
            comp = makeUncorrDiff(values['el'].graph_fit, values['mu'].graph_fit)
        else:
            yaxtitle = '#sigma_{#mu}/#sigma_{e}'
            yaxrange = (0.65, 1.35)
            comp = makeUncorrRatio(values['el'].graph_fit, values['mu'].graph_fit)
        comp.SetName("comp"); comp.SetTitle("")
        comp.SetMarkerColor(values['mu'].color)
        comp.SetLineColor(values['mu'].color)
        comp.SetFillColor(values['mu'].color)
        comp.SetMarkerStyle(ROOT.kFullCircle)
        comp.SetMarkerSize(2)

        comp.Draw('Pa')
        ## x axis fiddling
        comp.GetXaxis().SetRangeUser(0., options.maxRapidity)
        comp.GetXaxis().SetLabelSize(0.1)
        ## y axis fiddling
        comp.GetYaxis().SetTitleOffset(0.4)
        comp.GetYaxis().SetTitleSize(0.18)
        comp.GetYaxis().SetLabelSize(0.12)
        comp.GetYaxis().SetTitle(yaxtitle)
        comp.GetYaxis().SetRangeUser(yaxrange[0],yaxrange[1])
        comp.GetYaxis().CenterTitle()
        comp.GetYaxis().SetDecimals()
        lines["horiz_line_comp"] = ROOT.TF1("horiz_line_comp","0" if charge=='asymmetry' else '1',0.0,3.0);
        lines["horiz_line_comp"].SetLineColor(ROOT.kRed);
        lines["horiz_line_comp"].SetLineWidth(2);
        lines["horiz_line_comp"].SetLineStyle(ROOT.kDashed);
        lines["horiz_line_comp"].Draw("Lsame");

        subpadLegends["comp_lep"] = ROOT.TLegend(0.25, 0.8, 0.5, 0.9); subpadLegends["comp_lep"].SetFillStyle(0); subpadLegends["comp_lep"].SetBorderSize(0)
        subpadLegends["comp_lep"].SetNColumns(2)
        legtitle = "Data #mu - e" if charge=='asymmetry' else "Data #mu / e"
        subpadLegends["comp_lep"].AddEntry(comp,legtitle,'pl')
        subpadLegends["comp_lep"].Draw('same')

    c2.cd()
    # was 0.95 but for me it appears inside the frame
    lat.DrawLatex(0.2, 0.97, '#bf{CMS}') #it{Preliminary}')
    lat.DrawLatex(0.60, 0.97, '35.9 fb^{-1} (13 TeV)')
    lat.DrawLatex(0.25, 0.60,  'W^{{{ch}}} #rightarrow l^{{{ch}}}{nu}'.format(ch=ch,nu="#bar{#nu}" if charge=='minus' else "#nu"))
    lat.DrawLatex(0.85, 0.025, '|y_{W}|')
    for ext in ['png', 'pdf']: #, 'C', 'root']:
        c2.SaveAs('{n}.{ext}'.format(n=plotname, ext=ext))


NPDFs = 60
LUMINOSITY = 36000

if __name__ == "__main__":

    from optparse import OptionParser
    parser = OptionParser(usage='%prog ybinfile workspace.root toys.root [options] ')
    parser.add_option(      '--infile-mu'   , dest='infileMu' , default=''            , type='string', help='fitresults file with muon results')
    parser.add_option(      '--infile-el'   , dest='infileEl' , default=''            , type='string', help='fitresults file with electron results')
    parser.add_option(      '--infile-lep'  , dest='infileLep', default=''            , type='string', help='fitresults file with muon,electron combination results')
    parser.add_option('-y', '--ybinfile'    , dest='ybinfile' , default=''            , type='string', help='file with the yw binning')
    parser.add_option(      '--xsecfiles'   , dest='xsecfiles', default=None          , type='string', help='files that contains the expected x sections with variations (one per charge,comma separated in the same order of the charges) ')
    parser.add_option(      '--altxsecfiles', dest='altxsecfiles' , default=None      , type='string', help='files that contains alternative xsecs as the nominal to be drawn as a line ')
    parser.add_option('-C', '--charge'      , dest='charge'   , default='plus,minus'  , type='string', help='process given charge. default is both')
    parser.add_option('-o', '--outdir'      , dest='outdir'   , default='.'           , type='string', help='outdput directory to save the plots')
    parser.add_option(      '--suffix'      , dest='suffix'   , default=''            , type='string', help='suffix for the correlation matrix')
    parser.add_option('-n', '--normxsec'    , dest='normxsec' , default=False         , action='store_true',   help='if given, plot the differential xsecs normalized to the total xsec')
    parser.add_option(      '--nolong'      , dest='nolong'   , default=False         , action='store_true',   help='if given, do not plot longitudinal component')
    parser.add_option(      '--longBkg'     , dest='longBkg'  , default=False         , action='store_true',   help='if True, longitudinal component was treated as background, so the POIs are missing. Manage inputs accordingly')
    parser.add_option(     '--ybinsBkg', dest='ybinsBkg', type='string', default="10,11", help='Define which Y bins are to be considered as background. With format 14,15 ')
    parser.add_option(     '--ybinsOutAcc', dest='ybinsOutAcc', type='string', default="", help='Define which Y bins were put in OutAcc channel in the fit. With format 14,15 ')
    parser.add_option(      '--max-rap'     , dest='maxRapidity', default='2.5'       , type='float', help='Max value for rapidity range')
    (options, args) = parser.parse_args()


    if not os.path.isdir(options.outdir):
        os.system('mkdir {od}'.format(od=options.outdir))
        if os.path.exists("/afs/cern.ch"): os.system("cp /afs/cern.ch/user/m/mciprian/public/index.php {od}".format(od=options.outdir))

    if options.ybinfile:
        ybinfile = options.ybinfile
    else:
        ybinfile = os.path.dirname(os.path.abspath(options.infile))+'/binningYW.txt'


    ## get the central values and uncertainties (only Hessian implemented)
    valuesAndErrors = {}
    valuesAndErrors['mu'] = utilities.getFromHessian(options.infileMu)
    valuesAndErrors['el'] = utilities.getFromHessian(options.infileEl)
    valuesAndErrors['lep'] = utilities.getFromHessian(options.infileLep)

    ybinfile = open(ybinfile, 'r')
    ybins = eval(ybinfile.read())
    ybinfile.close()


    bkgYBins = []
    if options.ybinsBkg:
        bkgYBins = list(int(i) for i in options.ybinsBkg.split(','))        
    if options.longBkg:
        options.nolong = True

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

    charges = options.charge.split(',')
    xsecfiles = options.xsecfiles.split(',')
    doAltExp =  options.altxsecfiles
    xsec_nominal_allCharges = {}; xsec_systematics_allCharges = {}

    if doAltExp:
        alt_xsecfiles = options.altxsecfiles.split(',')
        alt_xsec_nominal_allCharges = {}

    polarizations = ['left','right','long']
    signal_polarizations = ['left','right']
    if not options.longBkg:
        signal_polarizations.append('long')
    flavors = ['mu','el','lep']

    nChan = 2
    channel = 'lep'

    for ic,charge in enumerate(charges):

        sign = 1. if charge=='plus' else -1.

        ## this gets the pdf central variation binned in the correct format
        xsec_nominal     = utilities.getXSecFromShapes(ybins,charge,xsecfiles[ic],0,nChan, polarizations=polarizations)
        xsec_qcdenvelope = utilities.getQCDScaleEnvelope(ybins,charge,xsecfiles[ic],nChan, polarizations=polarizations )
        xsec_alphas_syst = utilities.getQCDScaleEnvelope(ybins,charge,xsecfiles[ic],nChan, polarizations=polarizations, doAlphaS=True )
        xsec_nominal_allCharges[charge] = xsec_nominal

        if doAltExp: 
            alt_xsec_nominal = utilities.getXSecFromShapes(ybins,charge,alt_xsecfiles[ic],0,nChan, polarizations=polarizations )
            alt_xsec_nominal_allCharges[charge] = alt_xsec_nominal

        value_syst = {}
        for pol in polarizations:
            histos = []
            values = []
            for ip in xrange(1,NPDFs+1):
                # print "Loading polarization %s, histograms for pdf %d" % (pol,ip)
                ## this gets the pdf variations after correctly rebinning the YW
                xsec_pdf = utilities.getXSecFromShapes(ybins,charge,xsecfiles[ic],ip,nChan, polarizations=polarizations)
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
                        print "SOMETHING WENT WRONG WITH THIS PDF: %d HAS RELATIVE SYST = %f. SKIPPING !" % (ip,xsec_relsyst)
                    else:
                        xsec_totUp += math.pow(xsec_relsyst*xsec_nom,2)
                # now to the total PDF error, add the QCD scales envelope
                xsec_totUp += math.pow(xsec_qcdenvelope[pol][iy],2)
                # and add alphaS (which in principle is part of PDF variations)
                xsec_totUp += math.pow(xsec_alphas_syst[pol][iy],2)
                xsec_totUp = math.sqrt(xsec_totUp)
                # print "Rel systematic for Y bin %d = +/-%.3f" % (iy,totUp/nom)
                # print "\tRel systematic on xsec for Y bin %d = +/-%.3f" % (iy,xsec_totUp/xsec_nom if xsec_nom else 0.)
                xsec_systs.append(xsec_totUp)
            xsec_systematics[pol]=xsec_systs
        xsec_systematics_allCharges[charge] = xsec_systematics

        angcoeff_nominal = {'sumxsec': [], 'a0': [], 'a4': []}; alt_angcoeff_nominal = {'sumxsec': [], 'a0': [], 'a4': []}
        angcoeff_systematics = {'sumxsec': [], 'a0': [], 'a4': []}
        for  iy,y in enumerate(ybinwidths['{ch}_{pol}'.format(ch=charge,pol=pol if not pol=='long' else 'right')]):
            xsec_unpolarized_nominal_iy = sum([xsec_nominal[pol][iy] for pol in polarizations])
            angcoeff_nominal['sumxsec'].append(xsec_unpolarized_nominal_iy)
            xsec_unpolarized_iy = sum([xsec_systematics[pol][iy] for pol in polarizations])
            angcoeff_systematics['sumxsec'].append(xsec_unpolarized_iy)
            coeffs_val = utilities.getCoeffs(xsec_nominal['left'][iy],     xsec_nominal['right'][iy],     xsec_nominal['long'][iy],
                                             xsec_systematics['left'][iy], xsec_systematics['right'][iy], xsec_systematics['long'][iy])
            angcoeff_nominal['a0'].append(coeffs_val['a0'][0])
            angcoeff_nominal['a4'].append(sign*coeffs_val['a4'][0])
            angcoeff_systematics['a0'].append(coeffs_val['a0'][1])
            angcoeff_systematics['a4'].append(coeffs_val['a4'][1])
            if doAltExp:
                alt_xsec_unpolarized_nominal_iy = sum([alt_xsec_nominal[pol][iy] for pol in polarizations])
                alt_angcoeff_nominal['sumxsec'].append(alt_xsec_unpolarized_nominal_iy)
                alt_coeffs_val = utilities.getCoeffs(alt_xsec_nominal['left'][iy],     alt_xsec_nominal['right'][iy],     alt_xsec_nominal['long'][iy],
                                                     xsec_systematics['left'][iy],     xsec_systematics['right'][iy],     xsec_systematics['long'][iy])
                alt_angcoeff_nominal['a0'].append(alt_coeffs_val['a0'][0])
                alt_angcoeff_nominal['a4'].append(sign*alt_coeffs_val['a4'][0])

        nOuterBinsToExclude = 0  ### out of acceptance Y bins, or that were treated as background (not to be considered for the total xsec)
        if len(outAccYBins):
            nOuterBinsToExclude = len(outAccYBins)
        print "number of outer bins to exclude: " + str(nOuterBinsToExclude)

        cp = '{ch}_left'.format(ch=charge)
        MAXYFORNORM = ybins[cp][-nOuterBinsToExclude-1] # exclude the outermost 2 bins which has huge error due to acceptance
        print "MAXYFORNORM = " + str(MAXYFORNORM)
        normsigmaIn = sum([xsec_nominal[allpol][iy] for allpol in polarizations for iy,y in enumerate(ybins[cp][:-1]) if abs(y)<MAXYFORNORM])
        normsigmaOut = sum([xsec_nominal[allpol][iy] for allpol in polarizations for iy,y in enumerate(ybins[cp][:-1]) if abs(y)>=MAXYFORNORM])
        if doAltExp:
            alt_normsigmaIn  = sum([alt_xsec_nominal[allpol][iy] for allpol in polarizations for iy,y in enumerate(ybins[cp][:-1]) if abs(y)<MAXYFORNORM])
            alt_normsigmaOut = sum([alt_xsec_nominal[allpol][iy] for allpol in polarizations for iy,y in enumerate(ybins[cp][:-1]) if abs(y)>=MAXYFORNORM])

        allValues = {}
        normsigmaInFit = {}; normsigmaOutFit = {}
        for flav in flavors:
            nChan = 2 if flav=='lep' else 1
            for pol in signal_polarizations:
                print "total expected (fit for {flav}) xsec up to |Y|<{maxy} = {sigma:.3f} pb".format(flav=flav,maxy=MAXYFORNORM,sigma=normsigmaIn)
                if len(outAccYBins):
                    print "total expected (fit for {flav}) xsec beyond |Y|>{maxy} = {sigma:.3f} pb".format(flav=flav,maxy=MAXYFORNORM,sigma=normsigmaOut)
     
                tmp_val = valueClass('values_'+charge+'_'+pol+'_'+flav)
     
                for iy,y in enumerate(ybinwidths['{ch}_{pol}'.format(ch=charge,pol=pol)]):
                    if iy in bkgYBins: continue
                    # normsigma is used to normalize the expected: then, the sum should be the one on the expected
                    # which is also less sensitive to fluctuations, since in data some bins can be 0
                    normsigma = normsigmaIn if abs(ybins[cp][iy])<MAXYFORNORM else normsigmaOut
                    if doAltExp: alt_normsigma = alt_normsigmaIn if abs(ybins[cp][iy])<MAXYFORNORM else alt_normsigmaOut
                    parname = 'W{charge}_{pol}_Ybin_{iy}'.format(charge=charge,pol=pol,iy=iy)
     
                    scale = 1.
                    if options.normxsec:
                        xsec_fit = [x for x in valuesAndErrors[flav][parname+'_pmaskedexpnorm']]
                    else:
                        xsec_fit = [x/float(nChan) for x in valuesAndErrors[flav][parname+'_pmaskedexp']]
                        scale = LUMINOSITY
     
                    if options.normxsec:
                        rfit     = xsec_nominal[pol][iy]/normsigma/xsec_fit[0]
                    else:
                        rfit     = xsec_nominal[pol][iy]/xsec_fit[0]*scale
                    rfit_err = rfit*abs(xsec_fit[0]-xsec_fit[1])/xsec_fit[0] 
     
                    tmp_val.val.append(xsec_nominal[pol][iy]/ybinwidths[cp][iy])
                    if doAltExp: 
                        tmp_val.altval.append(alt_xsec_nominal[pol][iy]/ybinwidths[cp][iy])
                    tmp_val.ehi.append(xsec_systematics[pol][iy]/ybinwidths[cp][iy])
                    tmp_val.elo.append(xsec_systematics[pol][iy]/ybinwidths[cp][iy]) # symmetric for the expected
                    if options.normxsec:
                        tmp_val.val[-1] = tmp_val.val[-1]/normsigma
                        if doAltExp:
                            tmp_val.altval[-1] = tmp_val.altval[-1]/alt_normsigma
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
     
                allValues[(flav,pol)] = tmp_val

        plotValues(allValues,charge,channel,options)

        # now do the unpolarized ones
        cp = 'plus_left' # this works if the binning for all the pol is the same
        xsec_params = ['sumxsecnorm','a0','a4'] if options.normxsec else ['sumxsec']
        for xs in xsec_params:
            allValuesUnpol = {}
            for flav in flavors:
                nChan = 2 if flav=='lep' else 1
                tmp_val = valueClass('values_{xs}_{charge}_unpolarized'.format(xs=xs,charge=charge))
                normsigma = normsigmaIn if abs(ybins[cp][iy])<MAXYFORNORM else normsigmaOut
                if doAltExp: alt_normsigma = alt_normsigmaIn if abs(ybins[cp][iy])<MAXYFORNORM else alt_normsigmaOut
                for iy,y in enumerate(ybinwidths['{ch}_{pol}'.format(ch=charge,pol=pol)]):
                    if iy in bkgYBins: continue
                    parname = 'W{charge}_Ybin_{iy}_{xs}'.format(charge=charge,iy=iy,xs=xs)

                    ybinwidth_scale = 1.
                    scale = 1.
                    xskey=xs
                    if xs=='sumxsec':
                        ybinwidth_scale = ybinwidths[cp][iy]
                        scale = LUMINOSITY*float(nChan)
                    elif xs=='sumxsecnorm':
                        ybinwidth_scale = ybinwidths[cp][iy]
                        xskey = xs.replace('norm','') # to use the expected xsec values, which have to be normalized
        
                    tmp_val.val.append(abs(angcoeff_nominal[xskey][iy]/ybinwidth_scale))
                    if doAltExp:
                        tmp_val.altval.append(abs(alt_angcoeff_nominal[xskey][iy]/ybinwidth_scale))
                    experr = angcoeff_systematics[xskey][iy]/ybinwidth_scale
                    tmp_val.ehi.append(experr)
                    tmp_val.elo.append(experr) # symmetric for the expected

                    if xs=='sumxsecnorm':
                        tmp_val.val[-1] = tmp_val.val[-1]/normsigma
                        if doAltExp:
                            tmp_val.altval[-1] = tmp_val.altval[-1]/alt_normsigma
                        tmp_val.ehi[-1] = tmp_val.ehi[-1]/normsigma
                        tmp_val.elo[-1] = tmp_val.elo[-1]/normsigma
            
                    xsec_fit = valuesAndErrors[flav][parname]
            
                    tmp_val.val_fit.append(xsec_fit[0]/ybinwidth_scale/scale)
                    tmp_val.elo_fit.append(abs(xsec_fit[0]-xsec_fit[1])/ybinwidth_scale/scale)
                    tmp_val.ehi_fit.append(abs(xsec_fit[0]-xsec_fit[2])/ybinwidth_scale/scale)
        
                    tmp_val.relv. append(tmp_val.val[-1]/tmp_val.val_fit[-1])
                    experrrel = angcoeff_systematics[xskey][iy]/angcoeff_nominal[xskey][iy]

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
                allValuesUnpol[flav] = tmp_val
            plotUnpolarizedValues(allValuesUnpol,charge,channel,options)

    if len(charges)>1 and not options.normxsec:
        print "Making charge asymmetry plots now..."
        asymmetryValues = {}
        for flav in flavors:
            for pol in signal_polarizations:
                cp = 'plus_'+pol
                tmp_val = valueClass('asymmetry_'+pol)
                for iy,y in enumerate(ybinwidths[cp]):
                    if any(iy == x for x in bkgYBins): continue
                    chasy_val = utilities.getChargeAsy(xsec_nominal_allCharges['plus'][pol][iy],     xsec_nominal_allCharges['minus'][pol][iy],
                                                       xsec_systematics_allCharges['plus'][pol][iy], xsec_systematics_allCharges['minus'][pol][iy])
                    if doAltExp:
                        chasy_altval = utilities.getChargeAsy(alt_xsec_nominal_allCharges['plus'][pol][iy],     alt_xsec_nominal_allCharges['minus'][pol][iy],
                                                              xsec_systematics_allCharges['plus'][pol][iy],     xsec_systematics_allCharges['minus'][pol][iy]) # assumes the same error band for nominal and alternative, which is ~true for pt-unw and pt-w samples
                    tmp_val.val .append(chasy_val['asy'][0])
                    if doAltExp:
                        tmp_val.altval .append(chasy_altval['asy'][0])
                    tmp_val.ehi.append(chasy_val['asy'][1])
                    tmp_val.elo.append(chasy_val['asy'][1])
     
                    asy_fit = valuesAndErrors[flav]['W_{pol}_Ybin_{iy}_chargeasym'.format(pol=pol,iy=iy)]
                    # there is a crash when starting 'lep', don't know why
                    # so the comparison for asymmetry will be missing
                    #print "CHECK flav={flav} pol={pol} iy={iy}!!!".format(flav=flav,pol=pol,iy=iy)
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
                asymmetryValues[(flav,pol)] = tmp_val        
        plotValues(asymmetryValues,'asymmetry',channel,options)
                
        asymmetryValuesUnpol = {}
        for flav in flavors:
            # now do the unpolarized ones
            tmp_val = valueClass('values_asymmetry_unpolarized')
            for iy,y in enumerate(ybinwidths['plus_left']): # this assumes that all the 3 polarizations have the same binning
                if any(iy == x for x in bkgYBins): continue
                xval = {'plus': 0, 'minus': 0}; xerr = {'plus': 0, 'minus': 0}
                alt_xval = {'plus': 0, 'minus': 0}
                for charge in ['plus','minus']:
                    for pol in polarizations:
                        xval[charge] += xsec_nominal_allCharges[charge][pol][iy]
                        xerr[charge] += pow(xsec_systematics_allCharges[charge][pol][iy],2)
                        if doAltExp: alt_xval[charge] += alt_xsec_nominal_allCharges[charge][pol][iy]
                    xerr[charge] = math.sqrt(xerr[charge])
     
                chasy_val = utilities.getChargeAsy(xval['plus'], xval['minus'],
                                                   xerr['plus'], xerr['minus'])
                if doAltExp:
                    chasy_altval = utilities.getChargeAsy(alt_xval['plus'], alt_xval['minus'],
                                                          xerr['plus'], xerr['minus'])
                
                tmp_val.val .append(chasy_val['asy'][0])
                if doAltExp:
                    tmp_val.altval .append(chasy_altval['asy'][0])
                tmp_val.ehi.append(chasy_val['asy'][1])
                tmp_val.elo.append(chasy_val['asy'][1])
     
                asy_fit = valuesAndErrors[flav]['W_Ybin_{iy}_chargemetaasym'.format(iy=iy)]
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
            asymmetryValuesUnpol[flav] = tmp_val
        plotUnpolarizedValues(asymmetryValuesUnpol,'asymmetry',channel,options)
        
