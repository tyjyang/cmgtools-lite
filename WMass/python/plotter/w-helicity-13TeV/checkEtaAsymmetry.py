import ROOT, os, array

## usage
## python checkEtaAsymmetry.py -i <inputRootFile>

if __name__ == "__main__":

    from optparse import OptionParser
    parser = OptionParser(usage='%prog [options] ')
    parser.add_option('-i', '--infile'         , dest='infile'     , default=''      , type='string', help='input file with the plot and ratio')
    parser.add_option('-v', '--var'            , dest='var'        , default='eta'   , type='string', help='name of the plotted variable')
    (options, args) = parser.parse_args()

    #var = os.path.basename(options.infile).split('.')[0]

    inf = ROOT.TFile(options.infile,'read')

    sig, bkg = 0, 0
    lok = inf.GetListOfKeys()
    for k in lok:
        if not options.var in k.GetName(): 
            continue
        if options.var in k.GetName() and '_signal' in k.GetName():
            sig  = inf.Get(options.var+'_signal')
        if options.var in k.GetName() and '_background' in k.GetName():
            bkg  = inf.Get(options.var+'_background')
        if options.var in k.GetName() and '_data' in k.GetName():
            data = inf.Get(options.var+'_data')
        
    if sig and bkg:
        bkg.Add(sig)
    elif bkg and not sig:
        bkg = bkg
    elif not bkg and sig:
        bkg = sig
    else:
        print 'something went massively wrong!'

    ratio = data.Clone('ratio')
    ratio.Sumw2()
    ratio.Divide(bkg)

    nhalf = ratio.GetXaxis().GetNbins()/2

    binedges = []
    values = []; errors = []

    for ip in range(nhalf*2):
        values.append( ratio.GetBinContent(ip+1) )
        errors.append( ratio.GetBinError  (ip+1) )
        binedges.append( ratio.GetXaxis().GetBinLowEdge(ip+nhalf+1) )

    binedges = binedges[:nhalf]

    histleft  = ROOT.TH1F('hist_left' , '', len(binedges)-1, array.array('d',binedges) )
    histright = ROOT.TH1F('hist_right', '', len(binedges)-1, array.array('d',binedges) )
    histleft.SetTitle('data-MC ratio for +/- eta')

    for ib in range(nhalf):
        histleft.SetBinContent(ib+1, values[nhalf-ib-1])
        histleft.SetBinError  (ib+1, errors[nhalf-ib-1])

        histright.SetBinContent(ib+1, values[nhalf+ib])
        histright.SetBinError  (ib+1, errors[nhalf+ib])

    histleft  .SetMarkerStyle(20); histleft  .SetMarkerColor(ROOT.kBlue-1); histleft  .SetLineColor(ROOT.kBlue-1)
    histright .SetMarkerStyle(21); histright .SetMarkerColor(ROOT.kRed +2); histright .SetLineColor(ROOT.kRed +2)

    canv = ROOT.TCanvas('foob', 'data-MC ratio for +/- #eta', 800, 600)
    ROOT.gStyle.SetOptStat(0)

    histleft.Draw('pe')
    histright.Draw('same pe')
    histleft.GetYaxis().SetRangeUser(0.9, 1.1)
    histleft.GetYaxis().SetTitle('data/MC ratio')
    histleft.GetXaxis().SetTitle('#eta_{lepton}')

    leg = ROOT.TLegend(0.2, 0.7, 0.5, 0.9)
    leg.SetFillStyle(0)
    leg.SetLineWidth(0)
    leg.AddEntry(histleft , 'negative eta', 'pl')
    leg.AddEntry(histright, 'positive eta', 'pl')
    leg.Draw('same')

    line = ROOT.TLine()
    line.SetLineStyle(3)
    line.SetLineWidth(2)
    line.SetLineColor(ROOT.kBlack)
    line.DrawLine(binedges[0], 1., binedges[-1], 1.)

    canv.SaveAs(os.path.dirname(options.infile)+'/'+options.var+'_asymmetry.pdf')
    canv.SaveAs(os.path.dirname(options.infile)+'/'+options.var+'_asymmetry.png')

