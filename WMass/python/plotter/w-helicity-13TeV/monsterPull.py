import ROOT, math, os

## usage python monsterPull.py -i <inputfile> -d <distribution>


if __name__ == '__main__':
    
    from optparse import OptionParser
    parser = OptionParser(usage='%prog workspace ntoys [prefix] [options] ')
    parser.add_option('-i', '--infile'      , dest='infile'      , type='string',                     help='input file name');
    parser.add_option('-d', '--distribution', dest='distribution', type='string', default='unrolled', help='name of distribution in file');
    (options, args) = parser.parse_args()

    infile = ROOT.TFile(options.infile, 'read')

    hist_bkg = infile.Get(options.distribution+'_background')
    hist_sig = infile.Get(options.distribution+'_signal')

    hist_dat = infile.Get(options.distribution+'_data')

    hist_bkg.Add(hist_sig)
    

    pull_min, pull_max = -5., 5.

    pulls = ROOT.TH1F('pulls', '', 100, pull_min, pull_max)
    pulls.GetXaxis().SetTitle('pull')
    pulls.GetYaxis().SetTitle('# entries')
    pulls.GetXaxis().SetLabelSize(0.04), pulls.GetXaxis().SetTitleSize(0.05), pulls.GetXaxis().SetTitleOffset(0.90)
    pulls.GetYaxis().SetLabelSize(0.04), pulls.GetYaxis().SetTitleSize(0.05), pulls.GetYaxis().SetTitleOffset(0.90)
    pulls.SetLineColor(ROOT.kBlack), pulls.SetLineWidth(2)
    
    for ib in range(1,hist_bkg.GetNbinsX()+1):
        if hist_bkg.GetBinError(ib) == 0. or hist_dat.GetBinError(ib) == 0.:
            continue

        pull_num = hist_dat.GetBinContent(ib) - hist_bkg.GetBinContent(ib)
        pull_den = math.sqrt(hist_dat.GetBinError(ib)**2 + hist_bkg.GetBinError(ib)**2)

        pull = pull_num/pull_den

        if pull < pull_min: pull = pull_min
        if pull > pull_max: pull = pull_max

        pulls.Fill(pull)

    ROOT.gStyle.SetOptStat(0)
    
    canv = ROOT.TCanvas()

    pulls.Fit('gaus')
    pulls.Draw()

    pulls.GetFunction('gaus').SetLineColor(ROOT.kGreen+2)
    pulls.GetFunction('gaus').SetLineWidth(2)
    chi2 = pulls.GetFunction('gaus').GetChisquare()
    ndf  = pulls.GetFunction('gaus').GetNDF()

    lat = ROOT.TLatex(); lat.SetNDC(); lat.SetTextFont(42)

    lat.DrawLatex(0.10, 0.92, '#bf{CMS} #it{Preliminary}')
    lat.DrawLatex(0.68, 0.92, '36 fb^{-1} (13 TeV)')

    lat.SetTextSize(0.03)
    ymax = 0.85
    lat.DrawLatex(0.15, ymax-0.00, 'mean: {a:.2f}'.format(a=pulls.GetMean()))
    lat.DrawLatex(0.15, ymax-0.05, '#sigma: {a:.2f}'.format(a=pulls.GetStdDev()))
    lat.DrawLatex(0.15, ymax-0.10, 'mean_{{fit}}: {a:.2f}'.format(a=pulls.GetFunction('gaus').GetParameter(1)))
    lat.DrawLatex(0.15, ymax-0.15, '#sigma_{{fit}}: {a:.2f}'.format(a=pulls.GetFunction('gaus').GetParameter(2)))
    lat.DrawLatex(0.15, ymax-0.20, '#chi^{{2}}/ndf: {a:.2f}/{b:.0f} = {c:.2f}'.format(a=chi2,b=ndf,c=chi2/ndf))



    for ext in ['png', 'pdf']:
        canv.SaveAs(os.path.dirname(options.infile)+'/pulls_'+options.distribution+'.'+ext)
