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
    if 'alphaS' in params[0][0]:
        alphas = params[0]
        del params[0]
        params.append(alphas)
    return params


def getParamsFromRoot(infile):
    f = ROOT.TFile.Open(infile,"READ")
    if not f or not f.IsOpen():
        raise RuntimeError('Unable to open file {fn}'.format(fn=infile))
    h = f.Get("fit_s")
    if not h:
        raise RuntimeError('Unable to get histogram "fit_s" from file {fn}'.format(fn=infile))    
    params = []
    for ib in range(1,1+h.GetNbinsX()):
        key = h.GetXaxis().GetBinLabel(ib)
        pull = [h.GetBinContent(ib),h.GetBinError(ib)]
        params.append((key,pull[0],pull[1]))
    #if 'alphaS' in params[0][0]:
    #    alphas = params[0]
    #    del params[0]
    #    params.append(alphas)
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
    #niceName = 'W^{{{charge}}}_{{{pol}}} {nuis}^{{{index}}}'.format(charge=charge,pol=pol,nuis=nuis,index=index)
    # better to have index not as exponent, it overlaps with subscripts of previous label
    niceName = 'W^{{{charge}}}_{{{pol}}} {nuis} {index}'.format(charge=charge,pol=pol,nuis=nuis,index=index)

    return niceName

def niceNamePDFs(name):
    if 'pdf' in name:
        ihess = name.split('pdf')[-1]
        niceName = 'Hessian '+ihess
    elif 'alphaS' in name:
        niceName = '#alpha_{S}'
    return niceName
        

if __name__ == "__main__":

    from optparse import OptionParser
    parser = OptionParser(usage='%prog [options] ')
    parser.add_option('--exp'   , dest='expected', type="string", default='expected.latex', help='Text file for the expected, or root file');
    parser.add_option('--obs'   , dest='observed', type="string", default='observed.latex', help='Text file for the observed, or root file');
    parser.add_option('--outdir', dest='outdir'  , type="string", default='results'       , help='Output directory');
    parser.add_option('--name'   , dest='name', type="string", default='', help='name for output');
    (options, args) = parser.parse_args()

    if not options.name:
        print "Please pass a name for the output using option --name"
        quit()

    name = options.name
    if options.expected.endswith(".root"):
        print "Getting values from root file"
        exp_pulls = getParamsFromRoot(options.expected)
        obs_pulls = getParamsFromRoot(options.observed)
    else:
        exp_pulls = getParams(options.expected)
        obs_pulls = getParams(options.observed)
        
    npars = len(exp_pulls)
    
    x = array('d',[float(i+1) for i in range(npars)])
    
    zero   = array('d',[0.0 for i in range(npars)])
    one    = array('d',[1.0 for i in range(npars)])
    exone  = array('d', [0.4 for i in range(npars)])

    eyexp  = array('d',[float(p[2]) for p in exp_pulls])
    exexp  = array('d', [0.25 for i in range(npars)])

    yobs   = array('d',[float(p[1]) for p in obs_pulls])
    eyobs  = array('d',[float(p[2]) for p in obs_pulls])

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
    #dummyh.GetXaxis().LabelsOption('v') # gives warning: TAxis::Sort:0: RuntimeWarning: Cannot sort. No labels
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
            dummyh.GetXaxis().SetBinLabel(i+1,niceNamePDFs(exp_pulls[i][0]))
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
    leg.AddEntry(gr_prefit,'Pre-fit','f')
    
    gr_prefit.GetXaxis().SetNdivisions(npars+1,ROOT.kFALSE)
    gr_prefit.GetYaxis().SetNdivisions(2*int(maxz))
    
    # expected 
    gr_expected = ROOT.TGraphErrors(npars,x,zero,exexp,eyexp)
    gr_expected.GetYaxis().SetRangeUser(-maxz,maxz)
    gr_expected.GetXaxis().SetRangeUser(-0.5,npars+0.5)
    gr_expected.SetFillColor(ROOT.kOrange+7)
    gr_expected.SetMarkerColor(ROOT.kOrange+7)
    gr_expected.Draw('P2')
    leg.AddEntry(gr_expected,'Post-fit expected','f')
    
    # observed
    gr_observed = ROOT.TGraphErrors(npars,x,yobs,zero,eyobs)
    gr_observed.GetYaxis().SetRangeUser(-maxz,maxz)
    gr_observed.GetXaxis().SetRangeUser(-0.5,npars+0.5)
    gr_observed.SetMarkerStyle(ROOT.kFullCircle)
    gr_observed.SetMarkerSize(1)
    gr_observed.Draw('P EZ')
    leg.AddEntry(gr_observed,'Post-fit observed')

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
    
    # save graphs in root file to convert to hepdata later (no need to add prefit, trivial
    rfHepName = '{outdir}/comp_{name}.root'.format(outdir=options.outdir,name=name)
    rfHep = ROOT.TFile.Open(rfHepName,"recreate")
    if not rfHep:
        print "Error in plotExpObsPull.py: could not open root file %s" % rfHepName
        quit()
    rfHep.cd()
    dummyh.Write("xAxisLabels")
    gr_expected.Write("expected_postfit")
    gr_observed.Write("observed_postfit")
    rfHep.Close()

    # print "="*30
    # print "Printing obs_pulls for debugging"
    # for x in obs_pulls:
    #     print "%s    %.4f    %.4f" % (x[0], float(x[1]), float(x[2]))
    # gr_observed.Print()
    # print "="*30
    # print "Printing exp_pulls for debugging"
    # for x in exp_pulls:
    #     print "%s    %.4f    %.4f" % (x[0], float(x[1]), float(x[2]))
    # gr_expected.Print()
    # print "="*30
