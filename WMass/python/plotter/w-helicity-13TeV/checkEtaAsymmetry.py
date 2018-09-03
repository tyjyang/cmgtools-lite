import ROOT, os, array

## usage
## python checkEtaAsymmetry.py -i <inputRootFile>

if __name__ == "__main__":

    from optparse import OptionParser
    parser = OptionParser(usage='%prog [options] ')
    parser.add_option('-i', '--infile'         , dest='infile'     , default=''      , type='string', help='input file with the plot and ratio')
    #parser.add_option('-v', '--variable'       , dest='variable'   , default='etalep', type='string', help='name of the plotted variable')
    (options, args) = parser.parse_args()

    var = os.path.basename(options.infile).split('.')[0]

    inf = ROOT.TFile(options.infile,'read')
    bkg  = inf.Get(var+'_Z')
    data = inf.Get(var+'_data')

    ratio = data.Clone('ratio')
    ratio.Divide(bkg)

    graph = ROOT.TGraph(ratio)

    nhalf = graph.GetN()/2

    xs = []; ys = []
    xlo= []; ylo= []

    for ip in range(nhalf*2):
        a = ROOT.Double(1.)
        b = ROOT.Double(1.)
        graph.GetPoint(ip, a, b)
        print 'getting point ', a, b
        xs.append(a)
        ys.append(b)
        xlo.append(graph.GetErrorXlow(ip+1))
        ylo.append(graph.GetErrorYlow(ip+1))

    graphleft  = ROOT.TGraph(nhalf, array.array('d', [i for i in range(nhalf)]), array.array('d', ys[:nhalf][::-1]) )
    graphright = ROOT.TGraph(nhalf, array.array('d', [i for i in range(nhalf)]), array.array('d', ys[nhalf:]      ) )

    graphleft  .SetMarkerStyle(20); graphleft  .SetMarkerColor(ROOT.kBlue-1); graphleft  .SetLineColor(ROOT.kBlue-1)
    graphright .SetMarkerStyle(21); graphright .SetMarkerColor(ROOT.kRed +2); graphright .SetLineColor(ROOT.kRed +2)

    mg = ROOT.TMultiGraph()
    mg.Add(graphleft)
    mg.Add(graphright)

    canv = ROOT.TCanvas('foob', 'data-MC ratio for +/- #eta', 800, 600)
    mg.Draw('ap')

    mg.GetYaxis().SetRangeUser(0.8, 1.2)
    mg.GetYaxis().SetTitle('data/MC ratio')
    mg.GetXaxis().SetTitle('bin number in eta')

    leg = ROOT.TLegend(0.2, 0.7, 0.5, 0.9)
    leg.SetFillStyle(0)
    leg.SetLineWidth(0)
    leg.AddEntry(graphleft , 'negative eta', 'pl')
    leg.AddEntry(graphright, 'positive eta', 'pl')
    leg.Draw('same')

    canv.SaveAs(os.path.dirname(options.infile)+'/'+var+'_asymmetry.pdf')
    canv.SaveAs(os.path.dirname(options.infile)+'/'+var+'_asymmetry.png')

