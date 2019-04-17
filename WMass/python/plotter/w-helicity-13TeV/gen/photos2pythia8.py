import ROOT,os,re
ROOT.gROOT.SetBatch(True)

def makeCanvas():
    ROOT.gStyle.SetOptStat(0)
    c2 = ROOT.TCanvas('foo','', 600, 600)
    c2.Range(0,0,1,1);
    c2.SetFillColor(0);
    c2.SetBorderMode(0);
    c2.SetBorderSize(2);
    c2.SetTickx(1);
    c2.SetTicky(1);
    c2.SetLeftMargin(0.03);
    c2.SetRightMargin(0.05);
    c2.SetTopMargin(0.05);
    c2.SetBottomMargin(0.0);
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
        pad.SetFrameFillStyle(0);
        pad.SetFrameBorderMode(0);
        pad.SetFrameFillStyle(0);
        pad.SetFrameBorderMode(0);
        pads.append(pad)
    return pads


def makeRatios(args):
    
    infiles = args[0:2]
    outfile = args[2]

    ## take the histograms from the files
    histos = {}
    for filein in infiles:
        tf = ROOT.TFile(filein)
        for k in tf.GetListOfKeys():
            name = k.GetName()
            h = tf.Get(name)
            h.SetDirectory(None)
            if name in histos:
                histos[name].append(h)
            else:
                histos[name] = [h]
            
    ## make ratios of each variable
    for var,histoarr in histos.iteritems():
        for h in histoarr:
            h.Sumw2()
            if h.GetNbinsX()==1000:
                if var=='lptDressOverPreFSR': h.Rebin(4)
                else: h.Rebin(10)
            elif h.GetNbinsX()==5000: h.Rebin(50)
            else: h.Rebin(1)
            h.Scale(1./h.Integral())
            h.SetTitle('')
            if var=='lptDressOverPreFSR': h.GetXaxis().SetRangeUser(0.95,1.5)
            
        photos,pythia = histoarr
        photos.SetName(var+'_photos')
        pythia.SetName(var+'_pythia')
        ratio = photos.Clone(var+'_ratio')
        ratio.SetDirectory(None)
        ratio.Divide(pythia)
        histoarr.append(ratio)

    ## save everything in one file
    toutfile = ROOT.TFile(outfile,'recreate')
    for var,histoarr in histos.iteritems():
        for h in histoarr:
            h.Write()

    toutfile.Close()

    return histos



def plotRatios(histos,axislabels):
    
    subpadsYLow = [0.4,0.06]
    subpadsYUp = [0.95,0.4]

    for var,histoarr in histos.iteritems():
        canv = makeCanvas()
        pads = makePads(subpadsYLow,subpadsYUp)
        canv.cd()
        ## draw the single ones
        pads[0].Draw(); pads[0].cd(); ROOT.gPad.SetBottomMargin(0);
        if re.match('fsrpt.*',var) or var=='lptDressOverPreFSR':
            pads[0].SetLogy()
        photos,pythia,ratio = histoarr
        photos.SetMarkerColor(ROOT.kRed)
        photos.SetLineColor(ROOT.kRed)
        pythia.SetLineColor(ROOT.kBlack)
        photos.SetMaximum(1.1*max(photos.GetMaximum(),pythia.GetMaximum()))
        photos.Draw('pe')
        pythia.Draw('hist same')
        photos.GetYaxis().SetTitle('Normalized entries')
        photos.GetYaxis().SetTitleSize(0.07)
        photos.GetYaxis().SetDecimals()
        
        canv.cd()
        legend = ROOT.TLegend(0.25, 0.9, 0.90, 0.99);
        legend.SetFillStyle(0); legend.SetBorderSize(0)
        legend.SetNColumns(2)
        legend.AddEntry(photos,'PHOTOS','pl')
        legend.AddEntry(pythia,'PYTHIA 8','pl')
        legend.Draw()

        ## draw the ratio
        pads[1].Draw(); pads[1].cd(); ROOT.gPad.SetTopMargin(0); ROOT.gPad.SetBottomMargin(0.3); 
        ratio.SetMarkerColor(ROOT.kRed+1)
        ratio.SetLineColor(ROOT.kRed+1)
        ratio.GetYaxis().SetTitle('PHOTOS / PYTHIA 8')
        ratio.GetYaxis().SetTitleOffset(0.7)
        ratio.GetYaxis().SetTitleSize(0.08)
        ratio.GetYaxis().SetLabelSize(0.07)
        ratio.GetXaxis().SetLabelSize(0.09)
        ratio.GetYaxis().CenterTitle()
        ratio.GetYaxis().SetDecimals()
        ratio.Draw('pe')

        ratio.GetXaxis().SetTitleSize(0.12)
        ratio.GetXaxis().SetTitleOffset(0.9)
        if var in axislabels:
            ratio.GetXaxis().SetTitle(axislabels[var])
        else:
            ratio.GetXaxis().SetTitle(var)

        for ext in ['png','pdf']:
            canv.SaveAs("{var}.{ext}".format(var=var,ext=ext))
        
if __name__ == "__main__":

    from optparse import OptionParser
    parser = OptionParser(usage="%prog [options] photos.root pythia8.root photosOverPythia8.root ")
    (options, args) = parser.parse_args()
    if len(args)<3:
        print "Need photos.root pythia8.root photosOverPythia8.root"
        exit(0)
    
    histos = makeRatios(args)

    axis_labels = {'leta':'lepton #eta', 'lpt': 'lepton p_{T} [GeV]',
                   'wy':'W rapidity', 'wpt': 'W p_{T} [GeV]', 'wmass': 'W mass',
                   'fsrdr_close': '#Delta R(l,#gamma_{closest})', 'fsrpt_close': 'p_{T}(#gamma_{closest})',
                   'fsrdr_hard': '#Delta R(l,#gamma_{hardest})', 'fsrpt_hard': 'p_{T}(#gamma_{hardest})', 'fsrptfrac_hard': 'p_{T}(#gamma_{hardest})/p_{T}^{l}',
                   'nfsr': 'number of FSR #gamma',
                   'lptDressOverPreFSR': 'p_{T}^{l} / p_{T}^{pre-FSR l}'}

    plotRatios(histos,axis_labels)
