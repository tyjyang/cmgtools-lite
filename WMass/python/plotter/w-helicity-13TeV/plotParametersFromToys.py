#!/usr/bin/env python

## USAGE
## python plotParametersFromToys.py toys_wminus.root --pdir plots --suffix minus
## use option ' -a diffXsec ' for differential cross section analysis

import ROOT, os, re, datetime
from sys import argv, stdout, stderr, exit
ROOT.gROOT.SetBatch(True)

import utilities
utilities = utilities.util()

# import few functions from plotUtils/utility.h
from sys import path as syspath
syspath.append(os.getcwd() + "/plotUtils/")
from utility import addStringToEnd,createPlotDirAndCopyPhp

from make_diff_xsec_cards import get_ieta_ipt_from_process_name


lat = ROOT.TLatex(); lat.SetNDC()

def plotPars(inputFile, pois=None, selectString='', maxPullsPerPlot=30, plotdir='./', suffix='', analysis='helicity', excludeName=None):
 
    # warning: utilities.getHistosFromToys is working only for parameters centered around 0 (histogram is defined in -3,3)
    # this script will not work if using parameters with pmaskedexp
    if pois:
        if any (["pmaskedexp" in x for x in pois.split(",")]):
            print "Warning: current setup cannot work with parameters matching pmaskedexp: they are expected to be centered around 0. Exit"
            print "You might want to change behaviour of function utilities.getHistosFromToys to use those parameters as well" 
            exit(0)
    
    all_valuesAndErrors = utilities.getHistosFromToys(inputFile)
    
    matchString = ".*_mu_Ybin_.*" if analysis == "helicity" else ".*_mu_group_.*"
    channel = 'mu' if any(re.match(param,matchString) for param in all_valuesAndErrors.keys()) else 'el'
    print "From the list of parameters it seems that you are plotting results for channel ",channel
    
    print "pois = ", pois

    valuesAndErrors = {}

    excludedNames = excludeName.split(',') if excludeName else []

    if pois:
        poi_patts = pois.split(",")
        for ppatt in poi_patts:            
            for (k,v) in all_valuesAndErrors.iteritems():
                if any([re.match(x,k) for x in excludedNames]): continue
                if re.match(ppatt,k) and k not in valuesAndErrors: valuesAndErrors[k] = v
            #print valuesAndErrors.keys()

    params = valuesAndErrors.keys()
    print "----------------------------------"
    print "List of parameters:"
    print params
    print "----------------------------------"

    if any(re.match('pdf.*',x) for x in params):
        params = sorted(params, key = lambda x: int(x.split('pdf')[-1]), reverse=False)

    nPulls = 0
    pullLabelSize = 0.028
    pullSummaryMap = {}

    c = ROOT.TCanvas("c","",960,800)

    for name in params:

        print "Making pull for parameter ",name
        mean_p, sigma_p = (valuesAndErrors[name][0],valuesAndErrors[name][1]-valuesAndErrors[name][0])

        histo = valuesAndErrors[name][3]
        histo.Draw()
        histo.GetXaxis().SetTitle(histo.GetTitle());        histo.GetYaxis().SetTitle("no toys (%d total)" % histo.Integral());         histo.SetTitle("")
        histo.GetYaxis().SetTitleOffset(1.05);     histo.GetXaxis().SetTitleOffset(0.9);        histo.GetYaxis().SetTitleSize(0.05);        histo.GetXaxis().SetTitleSize(0.05);        histo.GetXaxis().SetTitle(name)
     
        fitPull = histo.Integral()>0
        if fitPull:
            histo.Fit("gaus","Q")
            fit = histo.GetFunction("gaus")
            fit.SetLineColor(4)
            lat.DrawLatex(0.12, 0.8, 'mean:     {me:.2f}'.format(me=fit.GetParameter(1)))
            lat.DrawLatex(0.12, 0.7, 'err :     {er:.2f}'.format(er=fit.GetParameter(2)))
            lat.DrawLatex(0.12, 0.6, 'chi2/ndf: {cn:.2f}'.format(cn=fit.GetChisquare()/fit.GetNDF() if fit.GetNDF()>0 else 999))
        
        os.system('cp /afs/cern.ch/user/g/gpetrucc/php/index.php '+plotdir)
        distrdir = plotdir+'/pulldistr'
        if not os.path.exists(distrdir): 
            os.system('mkdir '+distrdir)
            os.system('cp /afs/cern.ch/user/g/gpetrucc/php/index.php '+distrdir)
        for ext in ['png', 'pdf']:
            c.SaveAs("{pdir}/{name}_postfit_{suffix}_{channel}.{ext}".format(pdir=distrdir,name=name,suffix=suffix,channel=channel,ext=ext))

        if fitPull:
            tlatex = ROOT.TLatex(); tlatex.SetNDC(); 
            tlatex.SetTextSize(0.11)
            tlatex.SetTextColor(4);
            tlatex.DrawLatex(0.65,0.80,"Mean    : %.3f #pm %.3f" % (histo.GetFunction("gaus").GetParameter(1),histo.GetFunction("gaus").GetParError(1)))
            tlatex.DrawLatex(0.65,0.66,"Sigma   : %.3f #pm %.3f" % (histo.GetFunction("gaus").GetParameter(2),histo.GetFunction("gaus").GetParError(2)))
            
            tlatex.SetTextSize(0.11);
            tlatex.SetTextColor(1);
            tlatex.DrawLatex(0.65,0.33,"Post-fit #pm #sigma_{#theta}: %.3f #pm %.3f" % (mean_p, sigma_p))

            pullSummaryMap[name]=(histo.GetFunction("gaus").GetParameter(1),histo.GetFunction("gaus").GetParameter(2),
                                  histo.GetMean(),utilities.effSigma(histo))
            nPulls += 1
            
    if nPulls>0:
        print "Generating Pull Summaries...\n"
        nRemainingPulls = nPulls
        hc = ROOT.TCanvas("hc","",3000,2000); hc.SetGrid(0);
        pullPlots = 1;
        while nRemainingPulls > 0:
            nThisPulls = min(maxPullsPerPlot,nRemainingPulls)

            pull_rms      = ROOT.TH1F("pull_rms"     ,"", 3*nThisPulls+1,0,nThisPulls*3+1);
            pull_effsigma = ROOT.TH1F("pull_effsigma","", 3*nThisPulls+1,0,nThisPulls*3+1);
            pi=1
            sortedpulls = []
            if 'pdf' in pois:
                sortedpulls = sorted(pullSummaryMap.keys(), key=lambda k: int(k.split('pdf')[-1]))
            # following condition should not happen (parameter is not centered around 0)
            elif 'masked' in pois:
                keys = pullSummaryMap.keys()
                if analysis == "helicity":
                    keys_l = list(k for k in keys if 'left' in k)
                    keys_r = list(k for k in keys if 'right' in k)
                    norms_l = sorted(keys_l, key=lambda k: int(k.split('_')[-1]), reverse=False)
                    norms_r = sorted(keys_r, key=lambda k: int(k.split('_')[-1]), reverse=True)
                    sortedpulls = norms_r + norms_l
                else:
                    # differential xsection: get ieta, ipt index and use them as keys to sort
                    sortedpulls = sorted(keys, key = lambda x: get_ieta_ipt_from_process_name(x))
            else:
                sortedpulls = sorted(pullSummaryMap.keys(), key=lambda k: k[0])  # alphabetic order
            if len(sortedpulls)==0: break

            maxErr = 0;
            for name in sortedpulls:
                if pi>nThisPulls: break
                pull = pullSummaryMap[name]
                pull_rms     .SetBinContent(3*pi+0,pull[0]);  pull_rms     .SetBinError(3*pi+0,pull[1])
                pull_effsigma.SetBinContent(3*pi+1,pull[2]);  pull_effsigma.SetBinError(3*pi+1,pull[3])

                pull_rms.GetXaxis().SetBinLabel(3*pi,name)
                pull_rms.GetXaxis().SetTickSize(0.)
                
                max2=max(pull[1],pull[3])
                if max2>maxErr: maxErr=max2
                del pullSummaryMap[name]
                pi += 1
                nRemainingPulls -= 1

            #pull_rms .SetMarkerColor(46); pull_rms .SetLineColor(45); pull_rms .SetLineWidth(2); pull_rms .SetMarkerSize(1.0); pull_rms .SetMarkerStyle(21); 
            #pull_effsigma .SetMarkerColor(21); pull_effsigma .SetLineColor(24); pull_effsigma .SetLineWidth(2); pull_effsigma .SetMarkerSize(1.0); pull_effsigma .SetMarkerStyle(22)

            NmaxChar = 1
            for name in sortedpulls:
                NmaxChar = max(NmaxChar,len(name))
            if NmaxChar >= 8:
                hc.SetBottomMargin(0.2)  # hardcoded and not optimized I just looked at a random pull plot with the CMS_We_elescale and QCD scales parameters

            pull_rms .SetMarkerColor(ROOT.kRed+2); 
            pull_rms .SetLineColor(ROOT.kRed); pull_rms .SetLineWidth(3); 
            pull_rms .SetMarkerSize(1.0);      pull_rms .SetMarkerStyle(21); 

            pull_effsigma .SetMarkerColor(ROOT.kBlack); 
            pull_effsigma .SetLineColor(ROOT.kGray+2);   pull_effsigma .SetLineWidth(3); 
            pull_effsigma .SetMarkerSize(1.0);           pull_effsigma .SetMarkerStyle(22)

            pull_rms.SetLabelSize(pullLabelSize)
            pull_rms.LabelsOption("v")
            yrange = 3. 
            pull_rms.GetYaxis().SetRangeUser(-yrange,yrange)
            pull_rms.GetYaxis().SetTitle("pull summary (n#sigma)")
            pull_rms.Draw("p")
            pull_effsigma.Draw("p same")

            line0 = ROOT.TLine(pull_rms.GetXaxis().GetXmin(), 0., pull_rms.GetXaxis().GetXmax(), 0.); line0.SetLineStyle(7)
            line0.Draw('same')

            leg = ROOT.TLegend(0.60, 0.70, 0.85, 0.90)
            leg.SetFillStyle(0)
            leg.SetBorderSize(0)
            leg.AddEntry(pull_rms,"Gaussian #sigma")
            leg.AddEntry(pull_effsigma,"Effective #sigma")
            leg.Draw("same")
            param_group=pois.replace('.*','')
            for ext in ['png', 'pdf']:
                hc.SaveAs("{pdir}/pullSummaryToys_{params}_{igroup}_{suffix}_{c}.{ext}".format(pdir=plotdir,suffix=suffix,params=param_group,igroup=pullPlots,c=channel,ext=ext))
            pullPlots += 1


if __name__ == "__main__":

    ROOT.gROOT.SetBatch()
    ROOT.gStyle.SetOptStat(0)

    date = datetime.date.today().isoformat()

    from optparse import OptionParser
    parser = OptionParser(usage='%prog toys.root [options] ')
    parser.add_option(      '--parameters'  , dest='pois'     , default='pdf.*', type='string', help='comma separated list of regexp parameters to run. default is all parameters!')
    parser.add_option(      '--exclude-param'  , dest='excludeParam'     , default=None, type='string', help='Work as --parameters, but matches will be excluded from the list of parameters (e.g., can exclude a given pdf, pmaskednorm to match only pmasked and so on)')
    parser.add_option(      '--pdir'        , dest='plotdir'  , default='./'   , type='string', help='directory to save the likelihood scans')
    parser.add_option(      '--suffix'      , dest='suffix'   , default=''     , type='string', help='suffix to give to the plot files')
    parser.add_option('-a', '--analysis'    , dest='analysis' , default='helicity', type='string', help='Which analysis: helicity or diffXsec') 

    (options, args) = parser.parse_args()

    if len(args)<1: 
        print "need toyfile as argument. Exit."
        exit(0)

    if options.analysis not in ["helicity", "diffXsec"]:
        print "Warning: analysis not recognized, must be either \"helicity\" or \"diffXsec\""
        exit(0)
 
    outname = options.plotdir
    addStringToEnd(outname,"/",notAddIfEndswithMatch=True)
    createPlotDirAndCopyPhp(outname)

    toyfile = args[0]
    plotPars(toyfile,pois=options.pois,maxPullsPerPlot=30,plotdir=outname,suffix=options.suffix,analysis=options.analysis,excludeName=options.excludeParam)

