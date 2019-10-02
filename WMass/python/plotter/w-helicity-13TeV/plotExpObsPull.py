# usage: python plotExpObsPull.py --exp nuisances_pdf_fixedPOIs_hessian_bbb1_syst1_asimov.latex --obs nuisances_pdf_fixedPOIs_hessian_bbb1_syst1_data.latex 
import ROOT, random, os, sys, re
from array import array

ROOT.gROOT.SetBatch()

def getParams(infile):
    f = open(infile,'r')
    params = []
    for l in f:
        if l.strip().startswith('\\') or  l.strip().startswith('&'): continue
        values = l.strip().split(' ')
        key = values[0]
        pull = [float(values[-3].rstrip(',')),float(values[-2].rstrip('}'))]
        params.append((key,pull[0],pull[1]))
    return params

def niceNameQCDScales(name):
    pol = 'L' if name.startswith('left') else 'R' if name.startswith('right') else '0'
    charge = '+' if name.endswith('plus') else '-'
    nuis = name.replace('left','').replace('right','').replace('long','').replace('plus','').replace('minus','')
    nuis = nuis.replace('R','_{R}').replace('F','_{F}')
    nuis = nuis.replace('mu','#mu')
    m = re.match('(\D+)(\d+)',nuis)
    nuis = m.group(1)
    index = m.group(2)
    niceName = 'W^{{{charge}}}_{{{pol}}} {nuis}^{{{index}}}'.format(charge=charge,pol=pol,nuis=nuis,index=index)
    return niceName
    
if __name__ == "__main__":

    from optparse import OptionParser
    parser = OptionParser(usage='%prog [options] ')
    parser.add_option('--exp'   , dest='expected', type="string", default='expected.latex', help='Text file for the expected');
    parser.add_option('--obs'   , dest='observed', type="string", default='observed.latex', help='Text file for the observed');
    parser.add_option('--outdir', dest='outdir'  , type="string", default='results'       , help='Output directory');
    (options, args) = parser.parse_args()

    tokens = options.expected.split('_')
    name = '_'.join(tokens[1:3])

    exp_pulls = getParams(options.expected)
    obs_pulls = getParams(options.observed)
    npars = len(exp_pulls)
    
    x = array('f',[i+1 for i in range(npars)])

    zero   = array('f',[0 for i in range(npars)])
    one    = array('f',[1 for i in range(npars)])
    exone  = array('f', [0.4 for i in range(npars)])

    eyexp  = array('f',[p[2] for p in exp_pulls])
    exexp  = array('f', [0.25 for i in range(npars)])

    yobs   = array('f',[p[1] for p in obs_pulls])
    eyobs  = array('f',[p[2] for p in obs_pulls])

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
    dummyh = ROOT.TH1F('dummyh','',npars,x[0],x[-1])
    dummyh.GetYaxis().SetRangeUser(-maxz,maxz)
    dummyh.GetXaxis().SetRangeUser(-0.5,npars+1.5)
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
        dummyh.GetXaxis().SetTitle('PDF Hessian index')
    else:
        for i in range(npars):
            dummyh.GetXaxis().SetBinLabel(i+1,niceNameQCDScales(exp_pulls[i][0]))
    dummyh.Draw()

    # prefit (trivial, all 1)
    gr_prefit = ROOT.TGraphErrors(npars,x,zero,exone,one)
    gr_prefit.SetTitle('')
    gr_prefit.SetFillColor(ROOT.kAzure+6)
    gr_prefit.SetMarkerColor(ROOT.kAzure+6)
    gr_prefit.Draw('P2')
    leg.AddEntry(gr_prefit,'pre-fit','f')
    
    gr_prefit.GetXaxis().SetNdivisions(npars)
    gr_prefit.GetYaxis().SetNdivisions(2*int(maxz))
    
    # expected 
    gr_expected = ROOT.TGraphErrors(npars,x,zero,exexp,eyexp)
    gr_expected.GetYaxis().SetRangeUser(-maxz,maxz)
    gr_expected.GetXaxis().SetRangeUser(-0.5,npars+0.5)
    gr_expected.SetFillColor(ROOT.kOrange+7)
    gr_expected.SetMarkerColor(ROOT.kOrange+7)
    gr_expected.Draw('P2')
    leg.AddEntry(gr_expected,'post-fit expected','f')
    
    # observed
    gr_observed = ROOT.TGraphErrors(npars,x,yobs,zero,eyexp)
    gr_observed.GetYaxis().SetRangeUser(-maxz,maxz)
    gr_observed.GetXaxis().SetRangeUser(-0.5,npars+0.5)
    gr_observed.SetMarkerStyle(ROOT.kFullCircle)
    gr_observed.SetMarkerSize(1)
    gr_observed.Draw('P EZ')
    leg.AddEntry(gr_observed,'post-fit observed')

    lat.DrawLatex(marginL, 0.92, '#bf{CMS}')
    lat.DrawLatex(0.80, 0.92, '35.9 fb^{-1} (13 TeV)')

    leg.Draw('same')

    line = ROOT.TLine()
    for z in range(-int(maxz),int(maxz)):
        line.DrawLine(0, z, npars+1.5, z)
        line.SetLineStyle(3)
        line.SetLineColor(ROOT.kBlack)
    
    for ext in ['pdf','png']:
        c.SaveAs('{outdir}/comp_{name}.{ext}'.format(outdir=options.outdir,name=name,ext=ext))
    

    
