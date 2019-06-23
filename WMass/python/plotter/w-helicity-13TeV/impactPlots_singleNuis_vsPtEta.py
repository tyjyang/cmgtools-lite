#!/bin/env python

# script to plot impacts versus pt-eta for a single nuisance parameter or a single group
# EXAMPLE
# python w-helicity-13TeV/impactPlots_singleNuis_vsPtEta.py cards/diffXsec_mu_2019_06_17_zptReweight//fit/data/fitresults_123456789_Data_zptReweight_bbb1_cxs1.root -o plots/diffXsecAnalysis_new/muon/diffXsec_mu_2019_06_17_zptReweight//impactPlots_singleNuis_vsPtEta/zptReweight_bbb1_cxs1/  --suffix Data --nContours 51 --margin '0.16,0.15,0.05,0.25' --canvasSize '1000,1000' --splitOutByTarget  --etaptbinfile cards/diffXsec_mu_2019_06_17_zptReweight//binningPtEta.txt  -c mu  --nuisgroups EffStat --target asym --palette 55

import ROOT, os, datetime, re, operator, math
from array import array
from subMatrix import niceName
ROOT.gROOT.SetBatch(True)

from make_diff_xsec_cards import get_ieta_ipt_from_process_name
from make_diff_xsec_cards import getDiffXsecBinning
from make_diff_xsec_cards import templateBinning

from postFitPlots import prepareLegend

from impactPlots import latexLabel
from impactPlots import prepareLegendV2
from impactPlots import niceSystName


import utilities
utilities = utilities.util()

if __name__ == "__main__":

    ROOT.gStyle.SetOptStat(0)

    from optparse import OptionParser
    parser = OptionParser(usage='%prog fitresults.root [options] ')
    parser.add_option('-o','--outdir',     dest='outdir',     default='',   type='string', help='outdput directory to save the matrix')
    parser.add_option(     '--nuis',       dest='nuis',       default='',   type='string', help='single nuis for which you want to show the impacts versus pt-eta')
    parser.add_option(     '--nuisgroups', dest='nuisgroups', default='',   type='string', help='nuis groups for which you want to show the correlation matrix. comma separated list of regexps')
    parser.add_option(     '--suffix',     dest='suffix',     default='',   type='string', help='suffix for the correlation matrix')
    parser.add_option(     '--target',     dest='target',     default='mu', type='string', help='target POI (can be mu,xsec,xsecnorm, or other things)')
    parser.add_option(     '--zrange',     dest='zrange',     default='', type='string', help='Pass range for z axis as "min,max". If not given, it is set automatically based on the extremes of the impacts')
    parser.add_option(     '--canvasSize', dest='canvasSize', default='1200,1000', type='string', help='Pass canvas dimensions as "width,height" ')
    parser.add_option(     '--draw-option', dest='drawOption', default='COLZ', type='string', help='Options for drawing TH2')
    parser.add_option(     '--margin',     dest='margin',     default='', type='string', help='Pass canvas margin as "left,right,top,bottom" ')
    parser.add_option(     '--nContours', dest='nContours',    default=51, type=int, help='Number of contours in palette. Default is 51 (let it be odd, so the central strip is white if not using --abs-value and the range is symmetric)')
    parser.add_option(     '--palette'  , dest='palette',      default=0, type=int, help='Set palette: default is a built-in one, 55 is kRainbow')
    parser.add_option(     '--invertPalette', dest='invertePalette' , default=False , action='store_true',   help='Inverte color ordering in palette')
    parser.add_option(     '--parNameCanvas', dest='parNameCanvas',    default='', type='string', help='The canvas name is built using the parameters selected with --nuis or nuisgroups, but one could use a completely different name')
    parser.add_option(     '--abs-value', dest='absValue' , default=False , action='store_true',   help='Use absolute values for impacts (groups are already positive)')
    parser.add_option('-a','--absolute',   dest='absolute',   default=False, action='store_true', help='absolute uncertainty (default is relative)')
    parser.add_option(     '--latex',      dest='latex',      default=False, action='store_true', help='target POI (can be mu,xsec,xsecnorm)')
    parser.add_option(     '--etaptbinfile',   dest='etaptbinfile',   default='',  type='string', help='use this file with the |eta|-pT binning')
    parser.add_option(     '--splitOutByTarget', dest='splitOutByTarget' , default=False , action='store_true',   help='Create a subfolder appending target to output directory, where plots will be saved')
    parser.add_option(     '--pt-min-signal'  , dest='ptMinSignal',  default=-1, type=float, help='specify the minimum pt for the bins considered as signal. If not given, take the very first pt value from the gen pt binning')
    parser.add_option('-c', '--channel',     dest='channel',     default='',   type='string', help='Specify channel (el|mu) to make legend. It is no longer guessed from the name of POIs. If empty, generic "l" will be used in legends instead of (#mu|e)')
    (options, args) = parser.parse_args()

    # palettes:
    # 69 + inverted, using TColor::InvertPalette(), kBeach
    # 70 + inverted, using TColor::InvertPalette(), kBlackBody
    # 107, kCool
    # 57 kBird
    # 55 kRainBow

    if len(options.nuis) and len(options.nuisgroups):
        print 'You can specify either single nuis (--nuis) or poi groups (--nuisgroups), not both!'
        sys.exit()
    if len(options.nuis)==0 and len(options.nuisgroups)==0:
        print "You must specify one nuisance parameter or one group. Please use the full name,no regular expression"
        quit()

    if len(options.channel) and options.channel not in ["el","mu"]:
        print "Warning: if you specify a channel with -c, it must be 'el' or 'mu'"
        quit()

    if len(options.nuisgroups): options.absValue = True

    ROOT.TColor.CreateGradientColorTable(3,
                                         array ("d", [0.00, 0.50, 1.00]),
                                         ##array ("d", [1.00, 1.00, 0.00]),
                                         ##array ("d", [0.70, 1.00, 0.34]),
                                         ##array ("d", [0.00, 1.00, 0.82]),
                                         array ("d", [0.00, 1.00, 1.00]),
                                         array ("d", [0.34, 1.00, 0.65]),
                                         array ("d", [0.82, 1.00, 0.00]),
                                         255,  0.95)

    if options.absValue:
        ROOT.TColor.CreateGradientColorTable(2,
                                             array ("d", [0.00, 1.00]),
                                             ##array ("d", [1.00, 1.00, 0.00]),
                                             ##array ("d", [0.70, 1.00, 0.34]),
                                             ##array ("d", [0.00, 1.00, 0.82]),
                                             array ("d", [1.00, 1.00]),
                                             array ("d", [1.00, 0.65]),
                                             array ("d", [1.00, 0.00]),
                                             255,  0.95)

    # flavour is used in legend, so it is a TLatex string
    # for generic lepton, the Latex \ell is not supported (once should use TMathText which allows to write in Latex style and has \ell, but
    # apparently it is not supported for PDFs (one solution is that one should save in .eps and convert the file to pdf format later)
    # see: https://root-forum.cern.ch/t/tlatex-vs-ell/7767/11
    flavour = "l" 
    if len(options.channel):
        if options.channel == "el": flavour = "e"
        else                      : flavour = "#mu"

    if options.outdir:
        # emanuele: why this?
        #options.outdir = options.outdir + "/" + options.target + "/"
        # mciprian: for better bookkeeping, because it was easier to look at many plots
        if options.splitOutByTarget:
            options.outdir = options.outdir + "/" + options.target + "/"
        ROOT.gROOT.SetBatch()
        if not os.path.isdir(options.outdir):
            os.system('mkdir -p {od}'.format(od=options.outdir))
        os.system('cp {pf} {od}'.format(pf='/afs/cern.ch/user/g/gpetrucc/php/index.php',od=options.outdir))

    nuis = options.nuis if len(options.nuis) else options.nuisgroups

    hessfile = ROOT.TFile(args[0],'read')
    valuesAndErrors = utilities.getFromHessian(args[0])

    group = 'group_' if len(options.nuisgroups) else ''
    if   options.target=='xsec':       target = 'pmaskedexp'
    elif options.target=='xsecnorm':   target = 'pmaskedexpnorm'
    elif options.target=='asym':  target = 'chargepois'  # stands for etaptasym
    else:                              target = 'mu'

    pois_regexps = ["W.*_ieta_.*_ipt_.*"]

    th2name = 'nuisance_{group}impact_{sfx}'.format(group=group,sfx=target)
    impMat = hessfile.Get(th2name)
    if impMat==None:
        print "ERROR: Cannot find the impact TH2 named ",th2name," in the input file. Maybe you didn't run --doImpacts?\nSkipping."
        sys.exit()
    
    # apparently, target chargemetapois for diff.xsec contains labels ending with chargemetaasym and chargemetatotalxsec, 
    # so one should filter one of them according to options.target
    # same for chargepois, that has both chargeasym and chargetotalxsec

    # use dictionary for pois_indices, so to keep track of poi's name and allow for some sorting (for example to allow plus and minus in same plot)
    pois = []; pois_indices = {}
    for ib in xrange(1,impMat.GetNbinsX()+1):
        for poi in pois_regexps:
            if re.match(poi, impMat.GetXaxis().GetBinLabel(ib)):
                if options.target=='etaasym':
                    if re.match(".*chargemetaasym", impMat.GetXaxis().GetBinLabel(ib)):
                        pois.append(impMat.GetXaxis().GetBinLabel(ib))
                        pois_indices[impMat.GetXaxis().GetBinLabel(ib)] = ib
                elif options.target=='etaptasym':
                    if re.match(".*chargeasym", impMat.GetXaxis().GetBinLabel(ib)):
                        pois.append(impMat.GetXaxis().GetBinLabel(ib))
                        pois_indices[impMat.GetXaxis().GetBinLabel(ib)] = ib
                else:
                    pois.append(impMat.GetXaxis().GetBinLabel(ib))
                    pois_indices[impMat.GetXaxis().GetBinLabel(ib)] = ib


    nuisance = ""; nuisance_index = ""
    for ib in xrange(1,impMat.GetNbinsY()+1):
        if re.match(nuis, impMat.GetYaxis().GetBinLabel(ib)):
            nuisance = impMat.GetYaxis().GetBinLabel(ib)
            nuisance_index = ib

    mat = {}

    pois = sorted(pois, key= lambda x: get_ieta_ipt_from_process_name(x) if ('_ieta_' in x and '_ipt_' in x) else 0)
    # sort by charge if needed
    pois = sorted(pois, key = lambda x: 0 if "plus" in x else 1 if "minus" in x else 2)

    for ipoi, poi in enumerate(pois):        
            mat[(poi,nuisance)] = impMat.GetBinContent(pois_indices[poi], nuisance_index)

    ch = 1200
    cw = 1200
    if options.canvasSize:
        cw = int(options.canvasSize.split(',')[0])
        ch = int(options.canvasSize.split(',')[1])
    c = ROOT.TCanvas("c","",cw,ch)
    c.SetGridx()
    c.SetGridy()
    if options.nContours: ROOT.gStyle.SetNumberContours(options.nContours)
    if options.palette:   ROOT.gStyle.SetPalette(options.palette)
    if options.invertePalette:    ROOT.TColor.InvertPalette()


    clm = 0.2
    crm = 0.2
    cbm = 0.2
    ctm = 0.1
    if options.margin:
        clm,crm,ctm,cbm = (float(x) for x in options.margin.split(','))
    c.SetLeftMargin(clm)
    c.SetRightMargin(crm)
    c.SetBottomMargin(cbm)
    c.SetTopMargin(ctm)
    c.SetTickx(1)
    c.SetTicky(1)

    ## make the TH2F (eta on X and pT on Y)
    etaPtBinningVec = getDiffXsecBinning(options.etaptbinfile, "gen")                                                                                                    
    genBins = templateBinning(etaPtBinningVec[0],etaPtBinningVec[1])                                                                                                     

    charges = ["allCharges"] if re.match("asym",options.target) else ["plus", "minus"]
    for charge in charges:

        th2_sub = ROOT.TH2F('sub_imp_matrix_%s' % charge, 
                            'nuisance: %s' % niceSystName(nuisance), 
                            genBins.Neta, array("d",genBins.etaBins), genBins.Npt, array("d",genBins.ptBins))
        #th2_sub.GetXaxis().SetTickLength(0.)
        #th2_sub.GetYaxis().SetTickLength(0.)
        th2_sub.GetXaxis().SetTitle("%s |#eta|" % "muon" if options.channel == "mu" else "electron")
        th2_sub.GetYaxis().SetTitle("%s p_{T} [GeV]" % "muon" if options.channel == "mu" else "electron")

        th2_sub.GetXaxis().SetTitleSize(0.05)
        th2_sub.GetXaxis().SetLabelSize(0.04) 
        th2_sub.GetXaxis().SetTitleOffset(0.95) 

        th2_sub.GetYaxis().SetTitleSize(0.05)
        th2_sub.GetYaxis().SetLabelSize(0.04)
        th2_sub.GetYaxis().SetTitleOffset(1.1) 
        
        poiName_target = {"mu":        "signal strength",
                          "xsec":      "cross section",
                          "xsecnorm":  "normalized cross section",
                          "asym":   "charge asymmetry",
                          }
        th2_sub.GetZaxis().SetTitle("impact on POI for {p} {units}".format(units='' if options.absolute else '(%)', p=poiName_target[options.target]))
        th2_sub.GetZaxis().SetTitleOffset(1.4)
        th2_sub.GetZaxis().SetTitleSize(0.05)
        th2_sub.GetZaxis().SetLabelSize(0.04)                                                                                                                               

        for p in pois:
            if charge != "allCharges" and not charge in p: continue
            if '_ieta_' in p and '_ipt_' in p:
                if options.absolute: 
                    val = mat[(p,nuisance)]
                else: 
                    val = 100*mat[(p,nuisance)]/valuesAndErrors[p][0] if valuesAndErrors[p][0] !=0 else 0.0
                ieta,ipt = get_ieta_ipt_from_process_name(p)
                th2_sub.SetBinContent(ieta+1, ipt+1, val)

        rmax = max(abs(th2_sub.GetMaximum()),abs(th2_sub.GetMinimum()))
        if not options.absolute: rmax = min(20.,rmax)
        if options.absValue:
            th2_sub.GetZaxis().SetRangeUser(0,rmax)
        else :
            th2_sub.GetZaxis().SetRangeUser(-rmax,rmax)
        if options.zrange:
            rmin,rmax = (float(x) for x in options.zrange.split(','))
            th2_sub.GetZaxis().SetRangeUser(rmin,rmax)
        th2_sub.GetXaxis().LabelsOption("v")
        if options.absolute: 
            ROOT.gStyle.SetPaintTextFormat('1.3f')
            if options.target != "mu":
                ROOT.gStyle.SetPaintTextFormat('.3g')
        else: 
            ROOT.gStyle.SetPaintTextFormat('1.1f')

        th2_sub.Draw(options.drawOption)    # colz0 text45

        if options.parNameCanvas: 
            cname = options.parNameCanvas
        else:
            nuisName = ('NuisGroup'+options.nuisgroups) if len(options.nuisgroups) else options.nuis
            cname = "{nn}_vs_EtaPt".format(nn=nuisName)

        suff = '' if not options.suffix else '_'+options.suffix

        if options.latex:
            txtfilename = 'smallImpacts{rel}{suff}_{target}_{nn}_On_{pn}.tex'.format(rel='Abs' if options.absolute else 'Rel',suff=suff,target=options.target,i=i,cn=cname)
            if options.outdir: txtfilename = options.outdir + '/' + txtfilename
            txtfile = open(txtfilename,'w')
            txtfile.write("\\begin{tabular}{l "+"   ".join(['r' for i in xrange(th2_sub.GetNbinsX())])+"} \\hline \n")
            txtfile.write('                    ' + " & ".join(['{poi}'.format(poi=latexLabel(th2_sub.GetXaxis().GetBinLabel(i+1))) for i in xrange(th2_sub.GetNbinsX())]) + "\\\\ \\hline \n")
            for j in xrange(th2_sub.GetNbinsY()):
                txtfile.write('{label:<20}& '.format(label=th2_sub.GetYaxis().GetBinLabel(j+1)) + " & ".join(['{syst:^.2f}'.format(syst=th2_sub.GetBinContent(i+1,j+1)) for i in xrange(th2_sub.GetNbinsX())]) + "\\\\ \n")
            txtfile.write("\\end{tabular}\n")
            txtfile.close()
            print "Saved the latex table in file ",txtfilename

        if options.outdir and not options.latex:
            for i in ['pdf', 'png']:
                suff = '' if not options.suffix else ('_'+options.suffix)
                c.SaveAs(options.outdir+'/smallImpacts{rel}{suff}_{target}_{cn}_{ch}.{i}'.format(rel='Abs' if options.absolute else 'Rel',suff=suff,target=options.target,i=i,cn=cname,ch=charge))

            os.system('cp {pf} {od}'.format(pf='/afs/cern.ch/user/g/gpetrucc/php/index.php',od=options.outdir))

