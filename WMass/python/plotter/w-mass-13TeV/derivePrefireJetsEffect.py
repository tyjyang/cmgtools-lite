import ROOT


infile = ROOT.TFile('/afs/cern.ch/user/m/mdunser/www/private/w-helicity-13TeV/prefiring/2019-03-01/etal1Wleft.root','read')

nominal = infile.Get('etal1Wleft_Wleft')
prefire = infile.Get('etal1Wleft_Wleftprefire')


ratio = prefire.Clone('ratio')
ratio.Divide(nominal)


for i in range(1,ratio.GetXaxis().GetNbins()+1):
    edge = ratio.GetXaxis().GetBinUpEdge(i)
    val  = ratio.GetBinContent(i)
    print 'below {n:.1f} the value is {m:.4f}'.format(n=edge,m=val)
