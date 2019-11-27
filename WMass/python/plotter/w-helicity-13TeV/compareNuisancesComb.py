# usage: python plotExpObsPull.py --exp nuisances_pdf_fixedPOIs_hessian_bbb1_syst1_asimov.latex --obs nuisances_pdf_fixedPOIs_hessian_bbb1_syst1_data.latex 
import ROOT, random, os, sys, re
from array import array
from plotExpObsPull import getParams,niceNameQCDScales,niceNamePDFs

ROOT.gROOT.SetBatch()

if __name__ == "__main__":

    from optparse import OptionParser
    parser = OptionParser(usage='%prog [options] ')
    parser.add_option('--mu'    , dest='mu' , type="string", default='mu.latex' , help='Text file for the muon fit');
    parser.add_option('--el'    , dest='el' , type="string", default='el.latex' , help='Text file for the electron fit');
    parser.add_option('--lep'   , dest='lep', type="string", default='lep.latex', help='Text file for the mu+el fit');
    parser.add_option('--outdir', dest='outdir'  , type="string", default='results'       , help='Output directory');
    (options, args) = parser.parse_args()

    tokens = options.mu.split('_')
    name = '_'.join(tokens[1:3])

    mu_pulls  = getParams(options.mu)
    el_pulls  = getParams(options.el)
    lep_pulls = getParams(options.lep)
    npars = len(mu_pulls)
    
    x = array('f',[i+1 for i in range(npars)])

    zero   = array('f',[0 for i in range(npars)])

    yel  = array('f',[p[1] for p in el_pulls])
    eyel  = array('f',[p[2] for p in el_pulls])
    exel  = array('f', [0.4 for i in range(npars)])

    ymu  = array('f',[p[1] for p in mu_pulls])
    eymu  = array('f',[p[2] for p in mu_pulls])
    exmu  = array('f', [0.25 for i in range(npars)])

    ylep   = array('f',[p[1] for p in lep_pulls])
    eylep  = array('f',[p[2] for p in lep_pulls])

    c = ROOT.TCanvas('c','',1200,600)

    maxz = 5

    lat = ROOT.TLatex(); lat.SetNDC()
    lat.SetTextFont(42)
    lat.SetTextSize(0.05)

    marginR = 0.02
    marginL = 0.12
    c.SetRightMargin(marginR)
    c.SetLeftMargin(marginL)
    c.SetBottomMargin(0.21)
    c.SetTopMargin(0.10)
    
    leg = ROOT.TLegend(marginL, 0.02, 1.-marginR, 0.08)
    leg.SetNColumns(3)
    leg.SetColumnSeparation(0.1)
    #leg.SetFillStyle(0)
    leg.SetFillColor(0)
    leg.SetBorderSize(0)

    ROOT.gStyle.SetOptStat(0)
    dummyh = ROOT.TH1F('dummyh','',npars,x[0]-0.5,x[-1]+0.5)
    dummyh.GetYaxis().SetRangeUser(-maxz,maxz)
    dummyh.GetXaxis().SetRangeUser(-1.5,npars+1.5)
    dummyh.GetXaxis().LabelsOption('v')
    dummyh.GetYaxis().SetTitle('#theta - #theta^{0}')
    dummyh.GetYaxis().CenterTitle()
    dummyh.GetXaxis().SetTitleFont(42)
    dummyh.GetXaxis().SetTitleSize(0.05)
    dummyh.GetXaxis().SetLabelFont(42)
    dummyh.GetYaxis().SetTitleFont(42)
    dummyh.GetYaxis().SetTitleSize(0.07)
    dummyh.GetYaxis().SetTitleOffset(0.5)
    dummyh.GetYaxis().SetLabelFont(42)
    if 'pdf' in name:
        for i in range(npars):
            dummyh.GetXaxis().SetBinLabel(i+1,niceNamePDFs(mu_pulls[i][0]))
    elif 'mu' in name:
        for i in range(npars):
            dummyh.GetXaxis().SetBinLabel(i+1,niceNameQCDScales(mu_pulls[i][0]))
    else:
        for i in range(npars):
            dummyh.GetXaxis().SetBinLabel(i+1,mu_pulls[i][0])
    dummyh.Draw()

    # electrons
    gr_el = ROOT.TGraphErrors(npars,x,yel,exel,eyel)
    gr_el.GetYaxis().SetRangeUser(-maxz,maxz)
    gr_el.GetXaxis().SetRangeUser(-0.5,npars+0.5)
    gr_el.SetFillColorAlpha(ROOT.kAzure+6,0.8)
    gr_el.SetMarkerColor(ROOT.kAzure+6)
    gr_el.Draw('P2')
    leg.AddEntry(gr_el,'obs. el','f')
    
    gr_el.GetXaxis().SetNdivisions(npars+1,ROOT.kFALSE)
    gr_el.GetYaxis().SetNdivisions(2*int(maxz))

    # muons
    gr_mu = ROOT.TGraphErrors(npars,x,ymu,exmu,eymu)
    gr_mu.GetYaxis().SetRangeUser(-maxz,maxz)
    gr_mu.GetXaxis().SetRangeUser(-0.5,npars+0.5)
    gr_mu.SetFillColorAlpha(ROOT.kOrange+7,0.8)
    gr_mu.SetMarkerColor(ROOT.kOrange+7)
    gr_mu.Draw('P2')
    leg.AddEntry(gr_mu,'obs. #mu','f')

    # leptons
    gr_lep = ROOT.TGraphErrors(npars,x,ylep,zero,eymu)
    gr_lep.GetYaxis().SetRangeUser(-maxz,maxz)
    gr_lep.GetXaxis().SetRangeUser(-0.5,npars+0.5)
    gr_lep.SetMarkerStyle(ROOT.kFullCircle)
    gr_lep.SetMarkerSize(1)
    gr_lep.Draw('P EZ')
    leg.AddEntry(gr_lep,'obs. #mu + el')

    lat.DrawLatex(marginL, 0.92, '#bf{CMS}')
    lat.DrawLatex(0.80, 0.92, '35.9 fb^{-1} (13 TeV)')

    leg.Draw('same')

    line = ROOT.TLine()
    for z in range(-int(maxz),int(maxz)):
        line.DrawLine(0, z, npars+1.5, z)
        line.SetLineStyle(3)
        line.SetLineColor(ROOT.kBlack)
    
    for ext in ['pdf','png']:
        c.SaveAs('{outdir}/comb_{name}.{ext}'.format(outdir=options.outdir,name=name,ext=ext))
    

    
