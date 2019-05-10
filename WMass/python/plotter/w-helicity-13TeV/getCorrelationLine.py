import ROOT, os, datetime, re, operator, math
from array import array
ROOT.gROOT.SetBatch(True)

# python w-helicity-13TeV/getCorrelationLine.py cards/diffXsec_mu_2019_04_09_newSystAndWtau_fixTriSF/fit/hessian/fitresults_123456789_Asimov_combinedLep_bbb1_cxs1.root -o plots/diffXsecAnalysis/muon/diffXsec_mu_2019_04_09_newSystAndWtau_fixTriSF/getCorrelationLine/ -p CMS_Wmu_sig_lepeff -m sumpoisnorm -n 50

from make_diff_xsec_cards import get_ieta_ipt_from_process_name
from subMatrix import niceName
from operator import itemgetter

import utilities
utilities = utilities.util()



if __name__ == "__main__":

    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetPaintTextFormat('.3f')

    date = datetime.date.today().isoformat()

    from optparse import OptionParser
    parser = OptionParser(usage='%prog toys.root [options] ')
    parser.add_option('-o','--outdir', dest='outdir',    default='', type='string', help='outdput directory to save the matrix')
    parser.add_option('-p','--param', dest='param',    default='', type='string', help='parameter for which you want to show the correlation matrix. Must be a single object')
    parser.add_option('-t','--type'  , dest='type'  ,    default='hessian', type='string', help='which type of input file: toys or hessian (default)')
    parser.add_option('-m','--matrix', dest='matrix',    default='', type='string', help='matrix to be used (name is correlation_matrix_channel<matrix>)')
    parser.add_option(     '--suffix', dest='suffix',    default='', type='string', help='suffix for the correlation matrix')
    parser.add_option(     '--vertical-labels-X', dest='verticalLabelsX',    default=False, action='store_true', help='Set labels on X axis vertically (sometimes they overlap if rotated)')
    parser.add_option(     '--title'  , dest='title',    default='', type='string', help='Title for matrix. Use 0 to remove title. By default, string passed to option -p is used')
    parser.add_option('-n','--show-N' , dest='showN',    default=10, type=int, help='Show the N nuisances more correlated (in absolute value) with the parameter given with --param.')
    (options, args) = parser.parse_args()

    ROOT.TColor.CreateGradientColorTable(3,
                                      array ("d", [0.00, 0.50, 1.00]),
                                      ##array ("d", [1.00, 1.00, 0.00]),
                                      ##array ("d", [0.70, 1.00, 0.34]),
                                      ##array ("d", [0.00, 1.00, 0.82]),
                                      array ("d", [0.00, 1.00, 1.00]),
                                      array ("d", [0.34, 1.00, 0.65]),
                                      array ("d", [0.82, 1.00, 0.00]),
                                      255,  0.95)


    if not options.type in ['toys', 'hessian']:
        print 'the given type needs to be either "toys", "scans", or "hessian"!!'
        sys.exit()

    if not options.matrix:
        print "Need to specify which matrix with option -m (e.g. -m 'channelpmaskedexpnorm')"
        quit()

    if not options.param:
        print "Need to specify a parameter with option -p"
        quit()

    if len(args) < 1:
        print 'You have to pass a root file with the fit result, either for toys or hessian'
        sys.exit()

    if options.outdir:
        if not os.path.exists(options.outdir):
            os.makedirs(options.outdir)
        os.system('cp {pf} {od}'.format(pf='/afs/cern.ch/user/g/gpetrucc/php/index.php',od=options.outdir))

    param = options.param
    print "Will do POI with the following name: ",param
    

    corr = {}
    sign = {}
    index = 0

    ### GET LIST OF PARAMETERS THAT MATCH THE SPECIFIED OPTION IN THE TOYFILE
    if options.type == 'toys':
        print "Toys not implemented, sorry. Only hessian for now"
        quit()
    elif options.type == 'hessian':
        hessfile = ROOT.TFile(args[0],'read')
        suffix = options.matrix
        corrmatrix = hessfile.Get('correlation_matrix_channel'+suffix)
        covmatrix  = hessfile.Get('covariance_matrix_channel'+suffix)
        for ib in range(1+corrmatrix.GetNbinsX()+1):
            if re.match(param, corrmatrix.GetXaxis().GetBinLabel(ib)):
                ## store mean and rms into the dictionaries from before
                ## also keep a list of the parameter names, for sorting
                index = ib
            
    ## construct the covariances and the correlations in one go.
    for bin in range(1+corrmatrix.GetNbinsX()+1):
        label = corrmatrix.GetXaxis().GetBinLabel(bin)
        if label == param: continue
        # save absolute value, but keep track of the sign
        bincontent = corrmatrix.GetBinContent(index,bin)
        corr[label] = abs(bincontent)
        sign[label] = -1 if bincontent < 0 else 1

    #sorted_keys = sorted(corr.items(), key=itemgetter(1), reverse=True)
    sorted_keys = sorted(corr.items(), key=itemgetter(1), reverse=True)
    inum = 1
    hist = ROOT.TH1D("hist","",options.showN,0,options.showN)
    for key, val in sorted_keys:
        print "%s   %s" % (key, val*sign[key])
        hist.GetXaxis().SetBinLabel(inum,niceName(key))
        hist.SetBinContent(inum,val*sign[key])
        inum += 1        
        if inum > options.showN: break

    c = ROOT.TCanvas("c","",1200,800)
    c.SetTickx(1)
    c.SetTicky(1)

    c.SetLeftMargin(0.12)
    c.SetRightMargin(0.05)
    c.SetBottomMargin(0.3)

    if options.verticalLabelsX: hist.LabelsOption("v","X")
    if hist.GetNbinsX() >= 20: hist.LabelsOption("v","X")

    hist.SetTitle("parameter: " + param + "    channel: " + options.matrix.replace("channel",""))
    if len(options.title): 
        if options.title == "0":
            hist.SetTitle("")
        else:
            hist.SetTitle(options.title)

    hist.GetYaxis().SetTitle("Correlation")
    hist.SetLineWidth(2)
    hist.SetLineColor(ROOT.kGreen+2)
    hist.SetFillColor(ROOT.kGreen+1)
    hist.SetFillColorAlpha(ROOT.kGreen+1, 0.35)
    #hist.SetFillStyle(3001)
    hist.Draw("B")
    miny = hist.GetBinContent(hist.GetMinimumBin())
    maxy = hist.GetBinContent(hist.GetMaximumBin())
    maxval = max(abs(miny),abs(maxy))
    maxval *= 1.1
    hist.GetYaxis().SetRangeUser(-maxval, maxval)
    c.SetGridx(1)
    c.SetGridy(1)
    c.RedrawAxis("sameaxis")

    if options.outdir:
        for i in ['pdf', 'png']:
            suff = '' if not options.suffix else '_'+options.suffix
            c.SaveAs(options.outdir+'/corrLine{suff}_{pn}_{ch}.{i}'.format(suff=suff,i=i,pn=param, ch=options.matrix.replace("channel","")))
        os.system('cp {pf} {od}'.format(pf='/afs/cern.ch/user/g/gpetrucc/php/index.php',od=options.outdir))


