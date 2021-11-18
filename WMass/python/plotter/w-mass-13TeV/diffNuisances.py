#!/usr/bin/env python3

# example
# python w-mass-13TeV/diffNuisances.py --infile cards/wmass_fixMassWeights_splitW/fit/hessian/fitresults_123456789_Asimov_clipSyst1p3_bbb1_cxs1.root --outdir plots/testNanoAOD/WmassPlots_jetEta2p4_fixMassWeight_splitW/afterFitPlots/diffNuisances/ -a --format html --type hessian  --suffix Asimov --pois ".*" --uniqueString "allNuisances"
#
# python w-mass-13TeV/diffNuisances.py --infile cards/wmass_fixMassWeights_splitW/fit/hessian/fitresults_123456789_Asimov_clipSyst1p3_bbb1_cxs1.root --outdir plots/testNanoAOD/WmassPlots_jetEta2p4_fixMassWeight_splitW/afterFitPlots/diffNuisances/ -a --format html --type hessian  --suffix Asimov --pois "pdf.*|alphaS" --uniqueString "pdfsAndAlphaS"
#

import re, os
import datetime
from optparse import OptionParser

from subMatrix import niceName, niceNameHEPDATA

import argparse

import utilities
utilities = utilities.util()

## safe batch mode
import sys
args = sys.argv[:]
sys.argv = ['-b']
import ROOT
sys.argv = args
ROOT.gROOT.SetBatch(True)
ROOT.PyConfig.IgnoreCommandLineOptions = True

sys.path.append(os.getcwd() + "/plotUtils/")
from utility import *


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--infile'        , dest='infile'     , default=''        , type=str, help='file with the fitresult')
    parser.add_argument(      '--expected-infile'        , dest='expInfile'     , default=''        , type=str, help='file with the fitresult for expected, to plot together with observed (still to be implemented)')
    parser.add_argument('-o','--outdir', dest='outdir', default=None, type=str, help='If given, plot the pulls of the nuisances in this output directory')
    parser.add_argument("--vtol", "--val-tolerance", dest="vtol", default=0.30, type=float, help="Report nuisances whose value changes by more than this amount of sigmas")
    parser.add_argument("--vtol2", "--val-tolerance2", dest="vtol2", default=2.0, type=float, help="Report severely nuisances whose value changes by more than this amount of sigmas")
    parser.add_argument("-A", "--abs",      dest="absolute_values",    default=False,  action="store_true", help="Report also absolute values of nuisance values, not only the ones normalized to the input values")
    parser.add_argument("-a", "--all",      dest="show_all_parameters",    default=False,  action="store_true", help="Print all nuisances, even the ones which are unchanged w.r.t. pre-fit values.")
    parser.add_argument("-p", "--pois",      default=None,   type=str,  help="Name of the nuisances to be plotted (comma separated list of regexps)")
    parser.add_argument("-f", "--format", default="html", choices=['text', 'latex', 'html'], type=str,  help="Output format (in addition to pdf,png)")
    parser.add_argument(     '--suffix',default='', type=str, help='suffix for the correlation matrix')
    parser.add_argument(     '--uniqueString', default='', type=str, help='the output file names is nuisances_XXX_SUFF.EXT, with XXX obtained from the regular expression in input to select parameters, but can be replaced with this option (useful when using long regular expressions)')
    parser.add_argument('-t', '--type',  default='hessian', choices=["hessian", "scans", "toys"], type=str, help='run the plot from which postfit? toys/scans/hessian')
    parser.add_argument('--bm', '--bottom-margin' , dest='setBottomMargin'     , default=0.3        , type=float, help='Bottom margin for the canvas')
    parser.add_argument(     '--canvasSize',  default='', type=str, help='Pass canvas dimensions as "width,height". Default is 800,600, but it is automatically adjusted for large number of parameters')
    parser.add_argument(      "--y-title",      dest="ytitle",  default="S+B fit #theta",   type=str,  help="Title for Y axis")
    parser.add_argument(      "--y-setting",    dest="ysetting",  nargs=5, default=[-5.0,-3,0,3,5.0], type=float,  help="Settings to customize y axis: pass list of ymin,yhalfd,ycen,yhalfu,ymax, where horizontal lines are drawn")
    #parser.add_argument(      "--y-setting",    dest="ysetting",  type=lambda s: [float(item) for item in s.split(',')], default="-5.0,-3,0,3,5.0",  help="Settings to customize y axis: comma-separated list of ymin,yhalfd,ycen,yhalfu,ymax, where horizontal lines are drawn")
    parser.add_argument('-R', "--rank-nuisances-by", dest="rankNuisancesBy",  default=None, choices=["", "pull", "sigma"],  type=str,  help="Rank nuisances based on either sigma or pull. It is devised to work with --pois '.*', but of course you can further filter nuisances and/or pois")
    parser.add_argument('-N','--show-N' , dest='showN',    default=0, type=int, help='To be used with -R: it shows only the N nuisances ranked. If not positive, no limit is used')    
    parser.add_argument(     '--lower-limit-pull' , dest='lowerLimitPull', default=-1.0, type=float, help='To be used with -R. Take only nuisances with pull above this value (in absolute value). If negative, use no limit')    
    parser.add_argument(     '--upper-limit-sigma' , dest='upperLimitSigma', default=-1.0, type=float, help='To be used with -R. Take only nuisances with postfit sigma below this value . If negative, use no limit')    
    parser.add_argument(     '--use-hepdata-labels', dest='useHepdataLabels',    default=False, action='store_true', help='Write axis labels using latex for hepdata, with some name polishing')
    parser.add_argument('-c','--channel', default='mu', type=str, help='Channel (el|mu|lep), to force lepton flavor when writing some systematics in HEPDATA with --use-hepdata-labels')
    # following options are needed to prepare list of nuisances and POIs consistently with covariance matrices for hepdata
    parser.add_argument(     '--prepare-as-covariance-matrix', dest='prepareAsCovarianceMatrix',    default=False, action='store_true', help='Sort as in subMatrix.py to get better correspondance between matrix and list of POIs and nuisance parameters, when preparing material for hepdata')
    args = parser.parse_args()

    infile = args.infile
    infile_exp = args.expInfile
    plotObsWithExp = False
    if len(infile_exp):
        plotObsWithExp = True

    if args.outdir:
        createPlotDirAndCopyPhp(args.outdir)
    else:
        print("You must pass an output folder with option -o")
        quit()
        
    #valuesPrefit = dict((k,v) for k,v in valuesAndErrorsAll.items() if k.endswith('_gen'))
    pois_regexps = list(args.pois.split(','))

    if args.type == 'toys':
        valuesAndErrorsAll = utilities.getFromToys(infile,keepGen=True,params=pois_regexps)
        if plotObsWithExp:
            valuesAndErrorsAll_exp = utilities.getFromHessian(infile_exp,keepGen=True, params=pois_regexps)
    elif args.type == 'hessian':
        valuesAndErrorsAll = utilities.getFromHessian(infile,keepGen=True, params=pois_regexps)
        if plotObsWithExp:
            valuesAndErrorsAll_exp = utilities.getFromHessian(infile_exp,keepGen=True, params=pois_regexps)
    else:
        print("Sorry, scans method for option -t still has to be implemented")
        quit()
        
    print(f"Looking for regexp match {pois_regexps}")    
    valuesAndErrors = {}
    valuesErrors = {}
    valuesPrefit = {}

    for ppatt in pois_regexps:            
        for (k,v) in valuesAndErrorsAll.items():
            if re.match(ppatt,k):
                if k.endswith('_gen'):
                    valuesPrefit   [k]  = v #dict((k,v) for k,v in valuesPrefit.items() if re.match(ppatt.replace('$','_gen'),k))
                else:
                    valuesAndErrors[k]  = v #dict((k,v) for k,v in valuesAndErrors.items() if re.match(ppatt,k) and not k.endswith('_gen'))

    params = valuesAndErrors.keys()
    if len(params)==0:
        print("No parameters selected. Exiting.")
        quit()

    # if you are going to rank nuisances, just skip the sorting by names, to save some time
    if args.rankNuisancesBy:
        # this sorting will be managed later
        pass        
    else:
        params = sorted(params)
        if any(re.match('pdf.*',x) for x in params):
            # generally there will be alphaS along with pdfs
            # put alpha at last position
            params = sorted(params, key = lambda x: int(x.split('pdf')[-1]) if 'pdf' in x else 101 if 'alphaS' in x else 0, reverse=False)
        elif any(re.match('.*muR.*|.*muF.*',x) for x in params):
            params = sorted(params)
            params = sorted(params, key= lambda x: utilities.getNFromString(x))
            params = sorted(params, key= lambda x: int(re.sub('\D','',x)) if ('muRmuF' in x and x != "muRmuF")  else 0)
            params = sorted(params, key= lambda x: int(re.sub('\D','',x)) if ('muR' in x and x != "muR" and 'muRmuF' not in x) else 0)
            params = sorted(params, key= lambda x: int(re.sub('\D','',x)) if ('muF' in x and x != "muF" and 'muRmuF' not in x) else 0)
            params = sorted(params, key = lambda x: 0 if "Plus" in x else 1 if "Minus" in x else 2)
        elif any(re.match('FakesEtaUncorrelated.*',x) for x in params):
            params = sorted(params, key = lambda x: utilities.getNFromString(x), reverse=False)
        elif any(re.match('FakesPtUncorrelated.*',x) for x in params):
            params = sorted(params, key = lambda x: utilities.getNFromString(x), reverse=False)
        elif any(re.match('.*effStatTnP.*',x) for x in params):
            params = sorted(params, key = lambda x: utilities.getNFromString(x) , reverse=False)
        elif any(re.match('.*muonL1Prefire.*',x) for x in params):
            params = sorted(params, key = lambda x: utilities.getNFromString(x,chooseIndex=1) , reverse=False)

        print('='*30)
        print(f"Sorted params = {params}")


    nuis_p_i=0
    title = "#theta"


    hist_fit_s    = ROOT.TH1F("fit_s"   ,'',len(params),0,len(params))
    pmin, pmax = -3., 3.
    hist_fit_1d   = ROOT.TH1F("fit_1d " ,'',20,pmin,pmax)
    hist_fit_1d_e = ROOT.TH1F("fit_1d_e",'',59,pmin,pmax)
    hist_fit_1d   .SetLineColor(ROOT.kBlack); hist_fit_1d  .SetLineWidth(2)
    hist_fit_1d_e .SetLineColor(ROOT.kRed  ); hist_fit_1d_e.SetLineWidth(2)

    isFlagged = {}

    # maps from nuisance parameter name to the row to be printed in the table
    table = {}

    numberRankedNuisances = 0
    # loop over all fitted parameters
    for name in params:

        # keeps information to be printed about the nuisance parameter
        row = []
        flag = False

        ## this try catch catches the cases for which no gen (prefit) value exists. needed for the .* one
        try: 
            mean_p = valuesPrefit[name+'_gen'][0]
        except:
            # print(f"Exception caught: name = {name}")
            continue
        val_f,err_f = (valuesAndErrors[name][0],abs(valuesAndErrors[name][0]-valuesAndErrors[name][1]))

        if args.absolute_values:
            valShift = val_f
            sigShift = err_f
        else:
            valShift = val_f - mean_p
            sigShift = err_f

        if args.rankNuisancesBy:
            if args.lowerLimitPull > 0.0:
                if abs(valShift) <= args.lowerLimitPull: continue
            if args.upperLimitSigma > 0.0:
                if sigShift >= args.upperLimitSigma: continue
            numberRankedNuisances += 1

        if args.outdir: 
            nuis_p_i+=1
            hist_fit_s.SetBinContent(nuis_p_i,val_f)
            hist_fit_s.SetBinError(nuis_p_i,err_f)
            thisname = niceNameHEPDATA(name) if args.useHepdataLabels else niceName(name)
            hist_fit_s.GetXaxis().SetBinLabel(nuis_p_i, thisname)
            hist_fit_1d  .Fill(max(pmin,min(pmax-0.01,val_f)))
            hist_fit_1d_e.Fill(max(pmin,min(pmax-0.01,err_f-1.)))

        row += [" %+4.4f, %4.4f" % (valShift, sigShift)]

        if abs(val_f - mean_p) > args.vtol2*sigShift:

            # severely report this nuisance:
            # 
            # the best fit moved by more than 2.0 sigma or the uncertainty (sigma)
            # changed by more than 50% (default thresholds) w.r.t the prefit values
            isFlagged[name] = 2
            flag = True

        elif abs(val_f - mean_p) > args.vtol*sigShift:
            # report this nuisance:
            # 
            # the best fit moved by more than 0.3 sigma or the uncertainty (sigma)
            # changed by more than 10% (default thresholds) w.r.t the prefit values
            if args.show_all_parameters: isFlagged[name] = 1
            flag = True

        elif args.show_all_parameters:
            flag = True

        if flag or args.show_all_parameters: table[name] = row
    
    #end of loop over all fitted parameters

    #----------
    # print the results
    #----------
    if args.rankNuisancesBy:
        addSuffix = "rankBy{what}".format(what=args.rankNuisancesBy)
        if args.suffix:
            args.suffix += "_{add}".format(add=addSuffix)
        else:
            args.suffix = addSuffix
            
    uniquestr = ''.join(e for e in args.pois if e.isalnum()) ## removes all special characters from pois option
    outnameNoExt = "{od}/nuisances_{ps}_{suff}".format(od=args.outdir, ps=args.uniqueString if args.uniqueString != "" else uniquestr, suff=args.suffix)

    for ext in args.format.split(','):
        txtfilename = "{noext}.{ext}".format(noext=outnameNoExt, ext=ext)
        txtfile = open(txtfilename,'w')
     
        fmtstring = "%-40s     %15s"
        highlight = "*%s*"
        morelight = "!%s!"
        pmsub, sigsub = None, None
        if ext == 'text':
            txtfile.write(fmtstring % ('name', 's+b fit\n'))
        elif ext == 'latex':
            pmsub  = (r"(\S+) \+/- (\S+)", r"$\1 \\pm \2$")
            sigsub = ("sig", r"$\\sigma$")
            highlight = "\\textbf{%s}"
            morelight = "{{\\color{red}\\textbf{%s}}}"
            if args.absolute_values:
                fmtstring = "%-40s &  %15s & %30s \\\\"
                txtfile.write( "\\begin{tabular}{|l|r|r|r|} \\hline \n")
                txtfile.write( (fmtstring % ('name', 'pre fit', '$s+b$ fit')), " \\hline \n")
            else:
                fmtstring = "%-40s & %15s \\\\"
                txtfile.write( "\\begin{tabular}{|l|r|} \\hline \n")
                what = r"\Delta x/\sigma_{\text{in}}$, $\sigma_{\text{out}}/\sigma_{\text{in}}$"
                txtfile.write( fmtstring % ('', '$s+b$ fit\n'))
                txtfile.write((fmtstring % ('name', what)))
                txtfile.write(" \\hline \n")
        elif ext == 'html':
            pmsub  = (r"(\S+) \+/- (\S+)", r"\1 &plusmn; \2")
            sigsub = ("sig", r"&sigma;")
            highlight = "<b>%s</b>"
            morelight = "<strong>%s</strong>"
            txtfile.write("""
        <html><head><title>Comparison of nuisances</title>
        <style type="text/css">
            td, th { border-bottom: 1px solid black; padding: 1px 1em; }
            td { font-family: 'Consolas', 'Courier New', courier, monospace; }
            strong { color: red; font-weight: bolder; }
        </style>
        </head><body style="font-family: 'Verdana', sans-serif; font-size: 10pt;"><h1>Comparison of nuisances</h1>
        <table>
        """)
     
     
            if args.absolute_values:
                what = "x, &sigma;<sub>fit</sub>";
            else:
                what = "&Delta;x, &sigma;<sub>fit</sub>";
            txtfile.write("<tr><th>nuisance</th><th>signal fit<br/>%s</th>\n" % (what))
            fmtstring = "<tr><td><tt>%-40s</tt> </td><td> %-15s </td></tr>\n"
     
        names = table.keys()

        if args.rankNuisancesBy:    
            # rank by pull
            names = sorted(names, key = lambda x: abs(float((table[x][0].split(',')[0]).strip())), reverse=True)
            if args.rankNuisancesBy == "sigma":
                # sort by sigma, using the more constrained nuisances on top (i.e. lower sigma)
                names = sorted(names, key = lambda x: float((table[x][0].split(',')[1]).strip()))
        else:
            names = sorted(names) # first ordering
            #names = sorted(names, key = lambda x: int(x.replace('pdf','')) if 'pdf' in x else int(x.split('_')[-2]) if ('masked' in x or name.endswith('mu')) else -1)
            names = sorted(names, key= lambda x: utilities.getNFromString(x) if 'pdf' in x else 101 if 'alphaS' in x else 0)
            ## for mu* QCD scales, distinguish among muR and muRXX with XX in 1-10
            names = sorted(names, key= lambda x: int(re.sub('\D','',x)) if ('muRmuF' in x and x != "muRmuF")  else 0)
            names = sorted(names, key= lambda x: int(re.sub('\D','',x)) if ('muR' in x and x != "muR" and 'muRmuF' not in x) else 0)
            names = sorted(names, key= lambda x: int(re.sub('\D','',x)) if ('muF' in x and x != "muF" and 'muRmuF' not in x) else 0)
            names = sorted(names, key= lambda x: utilities.getNFromString(x) if 'effStatTnP' in x else 0)
            names = sorted(names, key = lambda x: utilities.getNFromString(x,chooseIndex=1) if 'muonL1Prefire' in x else 0)

            names = sorted(names, key= lambda x: utilities.getNFromString(x) if 'FakesEtaUncorrelated' in x else 0)
            names = sorted(names, key= lambda x: utilities.getNFromString(x) if 'FakesPtUncorrelated'  in x else 0)
     
        highlighters = { 1:highlight, 2:morelight };
        nameCounter = 0
        for n in names:
            if args.showN > 0 and nameCounter > args.showN: 
                break
            nameCounter += 1
            v = table[n]
            if ext == "latex": n = n.replace(r"_", r"\_")
            if pmsub  != None: v = [ re.sub(pmsub[0],  pmsub[1],  i) for i in v ]
            if sigsub != None: v = [ re.sub(sigsub[0], sigsub[1], i) for i in v ]
            if (n) in isFlagged: v[-1] = highlighters[isFlagged[n]] % v[-1]
            txtfile.write(fmtstring % (n, v[0]))
            txtfile.write('\n')
     
        if ext == "latex":
            txtfile.write(" \\hline\n\end{tabular}\n")
        elif ext == "html":
            txtfile.write("</table></body></html>\n")
        txtfile.close()
        print(f"Info: {ext} file {txtfilename} has been created")

    if args.outdir:
        line = ROOT.TLine()
        lat  = ROOT.TLatex(); lat.SetNDC(); lat.SetTextFont(42); lat.SetTextSize(0.04)
        ROOT.gStyle.SetOptStat(0)
        fout = ROOT.TFile.Open("{fn}.root".format(fn=outnameNoExt),"RECREATE")

        # need to shrink histogram, as some bins might have been removed when ranking. 
        # Also, use same sorting as table
        hist_fit_s_ranked = None
        nbins = len(params) # default
        if args.rankNuisancesBy:
            nbins = min(args.showN, len(names)) if args.showN > 0 else len(names)
            hist_fit_s_ranked    = ROOT.TH1F("fit_s_ranked"   ,'',nbins,0,nbins)
            index_name = 1
            for n in names:
                if args.showN > 0 and index_name > args.showN:
                    break 
                else:
                    thisname = niceNameHEPDATA(n) if args.useHepdataLabels else niceName(n)
                    binNumber = hist_fit_s.GetXaxis().FindFixBin(thisname)
                    if binNumber == -1:
                        print(f"Error when filling hist_fit_s_ranked. Could not find label {n} from hist_fit_s")
                        quit()
                    val = hist_fit_s.GetBinContent(binNumber)
                    err = hist_fit_s.GetBinError(binNumber)
                    hist_fit_s_ranked.SetBinContent(index_name,val)
                    hist_fit_s_ranked.SetBinError(index_name,err)
                    hist_fit_s_ranked.GetXaxis().SetBinLabel(index_name,n)
                    index_name += 1

        # customize canvas width a bit
        cw = 800
        ch = 600
        if nbins >= 100:
            cw = 2000            
        elif nbins >= 50:
            cw = 1200            
        if args.canvasSize:
            cw = int(args.canvasSize.split(',')[0])
            ch = int(args.canvasSize.split(',')[1])
        canvas_nuis = ROOT.TCanvas("nuisances", "nuisances", cw, ch)
        ## some style stuff
        #ymin,yhalfd,ycen,yhalfu,ymax = args.ysetting
        #ymin,yhalfd,ycen,yhalfu,ymax = map(float, list(args.ysetting.split(",")))
        ymin,yhalfd,ycen,yhalfu,ymax = args.ysetting
        if hist_fit_s_ranked != None:
            hist_fit_s = hist_fit_s_ranked

        canvas_nuis.SetTickx(1)
        canvas_nuis.SetTicky(1)
        hist_fit_s.GetYaxis().SetRangeUser(ymin-0.5,ymax+0.5)
        hist_fit_s.SetLineColor  (39)
        hist_fit_s.SetMarkerColor(ROOT.kGray+3)
        hist_fit_s.SetMarkerStyle(20)
        hist_fit_s.SetMarkerSize(1.0)
        hist_fit_s.SetLineWidth(2)
        hist_fit_s.Draw("PE1")
        hist_fit_s.GetXaxis().SetLabelSize(0.045 if nbins < 20 else 0.035)
        hist_fit_s.GetXaxis().LabelsOption("v")

        hist_fit_s.GetYaxis().SetTitle(args.ytitle)
        hist_fit_s.GetYaxis().SetTitleSize(0.05)
        hist_fit_s.GetYaxis().SetTitleOffset(0.50)

        clm = 0.1 if nbins < 100 else 0.05
        crm = 0.05 if nbins < 100 else 0.02
        cbm = args.setBottomMargin
        ctm = 0.1
        canvas_nuis.SetLeftMargin(clm)
        canvas_nuis.SetRightMargin(crm)
        canvas_nuis.SetBottomMargin(cbm)
        canvas_nuis.SetTopMargin(ctm)

        lat.DrawLatex(0.10, 0.92, '#bf{CMS} #it{Preliminary}')
        lat.DrawLatex(0.71 +(0.1-crm), 0.92, '36.3 fb^{-1} (13 TeV)')
        line.DrawLine(0., ycen, nbins, ycen)
        line.DrawLine(0., ymax, nbins, ymax)
        line.DrawLine(0., ymin, nbins, ymin)
        line.SetLineStyle(2); line.SetLineColor(ROOT.kRed)
        line.DrawLine(0., yhalfu, nbins, yhalfu)
        line.DrawLine(0., yhalfd, nbins, yhalfd)
        line.SetLineStyle(3)
        line.DrawLine(0., 1., nbins, 1.)
        line.DrawLine(0.,-1., nbins,-1.)
        hist_fit_s.Draw("PE1 same") ## draw again over the lines

        canvas_nuis.SetGridx()
        canvas_nuis.RedrawAxis()
        canvas_nuis.RedrawAxis('g')
        # leg=ROOT.TLegend(0.6,0.7,0.89,0.89)
        # leg.SetFillColor(0)
        # leg.SetTextFont(42)
        # leg.AddEntry(hist_fit_s,"S+B fit"   ,"EPL")
        # leg.Draw()
        #fout.WriteTObject(canvas_nuis)
        hist_fit_s.Write()
        fout.Close()

        for ext in ['png', 'pdf']:
            canvas_nuis.SaveAs("{noext}.{ext}".format(noext=outnameNoExt, ext=ext))
        if nbins >= 10:
            canv_constraints = ROOT.TCanvas('foobar', '', 800, 800)
            ## hist_fit_1d  .Scale(1./nbins)
            hist_fit_1d.GetXaxis().SetTitle('pulls and constraints')
            hist_fit_1d.GetYaxis().SetTitle('# of nuisance variations')
            #if hist_fit_1d.GetStdDev() > 0.01:
            hist_fit_1d  .SetMarkerSize(0.9)
            hist_fit_1d  .SetMarkerStyle(20)
            hist_fit_1d  .Draw("hist")
            hist_fit_1d.Fit('gaus')
            ff = hist_fit_1d.GetFunction('gaus')
            ff.SetLineColor(ROOT.kBlue-3)
            ff.SetLineWidth(2)
            ff.Draw('same')

            constraintColor = ROOT.kOrange+7
            canv_constraints.Update()
            rightmax = 1.5*hist_fit_1d_e.GetMaximum()
            scale = ROOT.gPad.GetUymax()/rightmax
            hist_fit_1d_e.SetLineColor(constraintColor)
            hist_fit_1d_e.SetLineWidth(2)
            hist_fit_1d_e.Scale(scale)
            hist_fit_1d_e.Draw("hist same")
            lat.DrawLatex(0.10, 0.92, '#bf{CMS}') #it{Preliminary}')
            lat.DrawLatex(0.61 , 0.92, '36.3 fb^{-1} (13 TeV)')
            ## draw an axis on the right side
            axis = ROOT.TGaxis(ROOT.gPad.GetUxmax(),ROOT.gPad.GetUymin(), ROOT.gPad.GetUxmax(), ROOT.gPad.GetUymax(),0,rightmax,510,"+L")
            axis.SetLineColor (constraintColor)
            axis.SetTextColor (constraintColor)
            axis.SetLabelColor(constraintColor)
            axis.Draw()

            leg1 = ROOT.TLegend(0.61, 0.6, 0.85, 0.85)
            leg1.SetLineWidth(0)
            leg1.SetFillStyle(0)
            leg1.SetTextSize(0.04)
            leg1.AddEntry(hist_fit_1d  , 'pulls', 'l')
            leg1.AddEntry(ff, '#splitline{{fit to pulls}}{{#scale[0.5]{{(#hat{{#mu}} = {mu:.2f}, #sigma = {si:.2f}) }} }}'.format(mu=ff.GetParameter(1), si=ff.GetParameter(2)), 'l')
            leg1.AddEntry(hist_fit_1d_e, '#splitline{{constraints}}{{ #scale[0.5]{{(mean = {a:.2f}) }} }}'.format(a=hist_fit_1d_e.GetMean()), 'pl')
            leg1.Draw('same')

            for ext in ['png', 'pdf']:
                canv_constraints.SaveAs("{noext}_pulls1D.{ext}".format(noext=outnameNoExt, ext=ext))

   

