import ROOT, copy


## usage
## python smoothWithGraph.py -i <infile> -o <outfile> -h <histNameToSmooth>

if __name__ == "__main__":

    from optparse import OptionParser
    parser = OptionParser(usage='%prog [options] ')
    parser.add_option('-i', '--infile'         , dest='infile'     , default=''      , type='string', help='input file with the plot and ratio')
    parser.add_option('-o', '--outfile'        , dest='outfile'    , default=''      , type='string', help='output file name/path')
    parser.add_option(      '--hist'           , dest='hist'       , default='scaleFactor'      , type='string', help='name of histogram to smooth')
    (options, args) = parser.parse_args()

    infile = ROOT.TFile(options.infile, 'read')
    
    
    hist2d = infile.Get(options.hist)
    
    
    histcopy = copy.deepcopy(hist2d)
    histcopy.SetName(histcopy.GetName()+'_smoothedByGraph')
    graph2d = ROOT.TGraph2D(histcopy)
    xbinsize = 0.05
    graph2d.SetNpx ( int((graph2d.GetXmax() - graph2d.GetXmin())/xbinsize) )

    ## don't smooth along pT
    #; ybinsize = 5.
    #smoothed_2dg.SetNpy( int((smoothed_2dg.GetYmax() - smoothed_2dg.GetYmin())/ybinsize) )

    blablabla = graph2d.GetHistogram() ## have to call this, otherwise root will freak out
    
    outfile = ROOT.TFile(options.outfile, 'recreate')
    
    outfile.cd()
    blablabla.Write()
    outfile.Close()

