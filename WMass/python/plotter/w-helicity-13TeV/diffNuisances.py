#!/usr/bin/env python
## USAGE:  python diffNuisances.py --infile fitresults.root --type hessian --pois ".*maskedexpnorm" --outdir plots -a --format html > xsecnorm.html
## add "--abs" to show the absolute value of the fitted nuisance and not the shift wrt the prefit

import re, os
from sys import argv, stdout, stderr, exit
import datetime
from optparse import OptionParser
#import HiggsAnalysis.CombinedLimit.calculate_pulls as CP 
from make_diff_xsec_cards import get_ieta_ipt_from_process_name
from make_diff_xsec_cards import get_ipt_from_process_name
from make_diff_xsec_cards import get_ieta_from_process_name
from subMatrix import niceName

import ROOT
ROOT.gROOT.SetBatch(True)
import utilities
utilities = utilities.util()


if __name__ == "__main__":

    parser = OptionParser(usage="usage: %prog [options] in.root  \nrun with --help to get list of options")
    parser.add_option("--vtol", "--val-tolerance", dest="vtol", default=0.30, type="float", help="Report nuisances whose value changes by more than this amount of sigmas")
    parser.add_option("--vtol2", "--val-tolerance2", dest="vtol2", default=2.0, type="float", help="Report severely nuisances whose value changes by more than this amount of sigmas")
    parser.add_option("-A", "--abs",      dest="absolute_values",    default=False,  action="store_true", help="Report also absolute values of nuisance values, not only the ones normalized to the input values")
    parser.add_option("-a", "--all",      dest="show_all_parameters",    default=False,  action="store_true", help="Print all nuisances, even the ones which are unchanged w.r.t. pre-fit values.")
    parser.add_option("-p", "--pois",      dest="pois",  default=None,   type="string",  help="Name of the nuisances to be fitted (comma separated list of regexps)")
    parser.add_option("-f", "--format",   dest="format", default="text", type="string",  help="Output format ('text', 'latex', 'html'")
    parser.add_option('-o','--outdir', dest='outdir', default=None, type='string', help='If given, plot the pulls of the nuisances in this output directory')
    parser.add_option(     '--suffix', dest='suffix', default='', type='string', help='suffix for the correlation matrix')
    parser.add_option('-t', '--type'        , dest='type'     , default='toys'        , type='string', help='run the plot from which postfit? toys/scans/hessian')
    parser.add_option('-i', '--infile'        , dest='infile'     , default=''        , type='string', help='file with the fitresult')
    parser.add_option('--bm', '--bottom-margin' , dest='setBottomMargin'     , default=0.3        , type='float', help='Bottom margin for the canvas')
    parser.add_option(     '--canvasSize', dest='canvasSize', default='', type='string', help='Pass canvas dimensions as "width,height". Default is 800,600, but it is automatically adjusted for large number of parameters')
    parser.add_option(      "--y-title",      dest="ytitle",  default="S+B fit #theta",   type="string",  help="Title for Y axis")

    (options, args) = parser.parse_args()
    infile = options.infile

    os.system('mkdir -p {od}'.format(od=options.outdir))
    os.system('cp ~emanuele/public/index.php {od}'.format(od=options.outdir))

    #valuesPrefit = dict((k,v) for k,v in valuesAndErrorsAll.iteritems() if k.endswith('_gen'))
    pois_regexps = list(options.pois.split(','))

    if   options.type == 'toys':
        valuesAndErrorsAll = utilities.getFromToys(infile,keepGen=True,params=pois_regexps)
    elif options.type == 'hessian':
        valuesAndErrorsAll = utilities.getFromHessian(infile,keepGen=True, params=pois_regexps)
    else:
        print 'ERROR: none of your types is supported. specify either "toys", "scans", or "hessian"'
        sys.exit()

    print 'looking for regexp match', pois_regexps    
    valuesAndErrors = {}
    valuesErrors = {}
    valuesPrefit = {}

    for ppatt in pois_regexps:            
        for (k,v) in valuesAndErrorsAll.iteritems():
            if re.match(ppatt,k):
                if k.endswith('_gen'):
                    valuesPrefit   [k]  = v #dict((k,v) for k,v in valuesPrefit.iteritems() if re.match(ppatt.replace('$','_gen'),k))
                else:
                    valuesAndErrors[k]  = v #dict((k,v) for k,v in valuesAndErrors.iteritems() if re.match(ppatt,k) and not k.endswith('_gen'))

    params = valuesAndErrors.keys()
    if len(params)==0:
        print "No parameters selected. Exiting."
        exit(1)

    params = sorted(params)
    if any(re.match('pdf.*',x) for x in params):
        # generally there will be alphaS along with pdfs
        params = sorted(params, key = lambda x: int(x.split('pdf')[-1]) if 'pdf' in x else 0, reverse=False)
    elif any(re.match('.*muR.*|.*muF.*',x) for x in params):
        params = sorted(params)
        # params = sorted(params, key= lambda x: int(re.sub('\D','',x)) if ('muRmuF' in x and x != "muRmuF")  else 0)
        # params = sorted(params, key= lambda x: int(re.sub('\D','',x)) if (''.join([j for j in x if not j.isdigit()]) == 'muR' and x != "muR") else 0)
        # params = sorted(params, key= lambda x: int(re.sub('\D','',x)) if (''.join([j for j in x if not j.isdigit()]) == 'muF' and x != "muF") else 0)
        params = sorted(params, key= lambda x: utilities.getNFromString(x))
        params = sorted(params, key = lambda x: 0 if "muRmuF" in x else 1 if "muR" in x else 2 if "muF" in x else 3)
        params = sorted(params, key = lambda x: 0 if "plus" in x else 1 if "minus" in x else 2)
    elif any(re.match('FakesEtaUncorrelated.*',x) for x in params):
        params = sorted(params, key = lambda x: utilities.getNFromString(x), reverse=False)
    elif any(re.match('FakesPtUncorrelated.*',x) for x in params):
        params = sorted(params, key = lambda x: utilities.getNFromString(x), reverse=False)
    elif any(re.match('.*EffStat.*',x) for x in params):
        params = sorted(params, key = lambda x: utilities.getNFromString(x) , reverse=False)
    elif any('masked' in x or x.endswith('mu') for x in params):
        if any('_ieta_' in x for x in params):
            # differential xsection: get ieta, ipt index and use them as keys to sort
            params = sorted(params, key = lambda x: -1 if "_ieta_" not in x else get_ieta_ipt_from_process_name(x), reverse=False )
        else:
            params = sorted(params, key = lambda x: int(x.split('_')[-2]) if ('masked' in x or x.endswith('_mu')) else -1, reverse=False)

    nuis_p_i=0
    title = "#theta"

    print "sorted params =",params

    #hist_fit_s  = ROOT.TH1F("fit_s"   ,"S+B fit Nuisances   ;;%s "%title,len(params),0,len(params))
    hist_fit_s    = ROOT.TH1F("fit_s"   ,'',len(params),0,len(params))
    pmin, pmax = -3., 3.
    hist_fit_1d   = ROOT.TH1F("fit_1d " ,'',20,pmin,pmax)
    hist_fit_1d_e = ROOT.TH1F("fit_1d_e",'',59,pmin,pmax)
    hist_fit_1d   .SetLineColor(ROOT.kBlack); hist_fit_1d  .SetLineWidth(2)
    hist_fit_1d_e .SetLineColor(ROOT.kRed  ); hist_fit_1d_e.SetLineWidth(2)

    isFlagged = {}

    # maps from nuisance parameter name to the row to be printed in the table
    table = {}

    # loop over all fitted parameters
    for name in params:

        # keeps information to be printed about the nuisance parameter
        row = []
        flag = False

        ## this try catch catches the cases for which no gen (prefit) value exists. needed for the .* one
        try: 
            mean_p = valuesPrefit[name+'_gen'][0]
        except:
            continue
        val_f,err_f = (valuesAndErrors[name][0],abs(valuesAndErrors[name][0]-valuesAndErrors[name][1]))
        row += [ "%+.4f +/- %.4f" % (val_f, err_f) ]        
        if options.outdir: 
            nuis_p_i+=1
            hist_fit_s.SetBinContent(nuis_p_i,val_f)
            hist_fit_s.SetBinError(nuis_p_i,err_f)
            hist_fit_s.GetXaxis().SetBinLabel(nuis_p_i,niceName(name.replace('CMS_','')))
            hist_fit_1d  .Fill(max(pmin,min(pmax-0.01,val_f)))
            hist_fit_1d_e.Fill(max(pmin,min(pmax-0.01,err_f-1.)))

        if options.absolute_values:
            valShift = val_f
            sigShift = err_f
        else:
            valShift = val_f - mean_p
            sigShift = err_f
        row[-1] = " %+4.4f, %4.4f" % (valShift, sigShift)

        if abs(val_f - mean_p) > options.vtol2*sigShift:

            # severely report this nuisance:
            # 
            # the best fit moved by more than 2.0 sigma or the uncertainty (sigma)
            # changed by more than 50% (default thresholds) w.r.t the prefit values

            isFlagged[name] = 2

            flag = True

        elif abs(val_f - mean_p) > options.vtol*sigShift:
            # report this nuisance:
            # 
            # the best fit moved by more than 0.3 sigma or the uncertainty (sigma)
            # changed by more than 10% (default thresholds) w.r.t the prefit values

            if options.show_all_parameters: isFlagged[name] = 1

            flag = True

        elif options.show_all_parameters:
            flag = True

        if flag or options.show_all_parameters: table[name] = row
    
    #end of loop over all fitted parameters

    #----------
    # print the results
    #----------
    for ext in options.format.split(','):
        uniquestr = ''.join(e for e in options.pois if e.isalnum()) ## removes all special characters from pois option
        txtfilename = "{od}/nuisances_{ps}_{suff}.{ext}".format(od=options.outdir, ps=uniquestr, suff=options.suffix, ext=ext)
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
            if options.absolute_values:
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
     
     
            if options.absolute_values:
                what = "x, &sigma;<sub>fit</sub>";
            else:
                what = "&Delta;x, &sigma;<sub>fit</sub>";
            txtfile.write("<tr><th>nuisance</th><th>signal fit<br/>%s</th>\n" % (what))
            fmtstring = "<tr><td><tt>%-40s</tt> </td><td> %-15s </td></tr>\n"
     
        names = table.keys()
        names = sorted(names) # first ordering
        #names = sorted(names, key = lambda x: int(x.replace('pdf','')) if 'pdf' in x else int(x.split('_')[-2]) if ('masked' in x or name.endswith('mu')) else -1)
        names = sorted(names, key= lambda x: utilities.getNFromString(x) if 'pdf' in x else 0)
        ## for mu* QCD scales, distinguish among muR and muRXX with XX in 1-10
        names = sorted(names, key= lambda x: -1 if "_ieta_" not in x else get_ieta_ipt_from_process_name(x) )
        names = sorted(names, key= lambda x: int(re.sub('\D','',x)) if ('muRmuF' in x and x != "muRmuF")  else 0)
        names = sorted(names, key= lambda x: int(re.sub('\D','',x)) if (''.join([j for j in x if not j.isdigit()]) == 'muR' and x != "muR") else 0)
        names = sorted(names, key= lambda x: int(re.sub('\D','',x)) if (''.join([j for j in x if not j.isdigit()]) == 'muF' and x != "muF") else 0)
        names = sorted(names, key= lambda x: utilities.getNFromString(x) if 'EffStat' in x else 0)
        names = sorted(names, key= lambda x: utilities.getNFromString(x) if 'FakesEtaUncorrelated' in x else 0)
        names = sorted(names, key= lambda x: utilities.getNFromString(x) if 'FakesPtUncorrelated'  in x else 0)

        if options.pois in ['.','.*']:
            names = sorted(names, key = lambda x: abs(float(table[x][0].split(',')[0])), reverse=True)
     
        highlighters = { 1:highlight, 2:morelight };
        for n in names:
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
        print "Info: ",ext," file ",txtfilename," has been created"

    if options.outdir:
        import ROOT
        line = ROOT.TLine()
        lat  = ROOT.TLatex(); lat.SetNDC(); lat.SetTextFont(42); lat.SetTextSize(0.04)
        ROOT.gROOT.SetBatch()
        ROOT.gStyle.SetOptStat(0)
        fout = ROOT.TFile("{od}/postfit_{suff}.root".format(od=options.outdir, suff=options.suffix),"RECREATE")

        # customize canvas width a bit
        cw = 800
        ch = 600
        if len(params) >= 100:
            cw = 2000            
        elif len(params) >= 50:
            cw = 1200            
        if options.canvasSize:
            cw = int(options.canvasSize.split(',')[0])
            ch = int(options.canvasSize.split(',')[1])
        canvas_nuis = ROOT.TCanvas("nuisances", "nuisances", cw, ch)
        ## some style stuff
        if '_mu' in options.pois:
            ymin,yhalfd,ycen,yhalfu,ymax = (0.,0.5,1,1.5,2.)
        else:
            ymin,yhalfd,ycen,yhalfu,ymax = (-5.,-3,0,3,5.)

        canvas_nuis.SetTickx(1)
        canvas_nuis.SetTicky(1)
        hist_fit_s.GetYaxis().SetRangeUser(ymin-0.5,ymax+0.5)
        hist_fit_s.SetLineColor  (39)
        hist_fit_s.SetMarkerColor(ROOT.kGray+3)
        hist_fit_s.SetMarkerStyle(20)
        hist_fit_s.SetMarkerSize(1.0)
        hist_fit_s.SetLineWidth(2)
        hist_fit_s.Draw("PE1")
        hist_fit_s.GetXaxis().SetLabelSize(0.045 if len(params) < 20 else 0.035)
        hist_fit_s.GetXaxis().LabelsOption("v")

        hist_fit_s.GetYaxis().SetTitle(options.ytitle)
        hist_fit_s.GetYaxis().SetTitleSize(0.05)
        hist_fit_s.GetYaxis().SetTitleOffset(0.80)

        clm = 0.1 if len(params) < 100 else 0.05
        crm = 0.05 if len(params) < 100 else 0.02
        cbm = options.setBottomMargin
        ctm = 0.1
        canvas_nuis.SetLeftMargin(clm)
        canvas_nuis.SetRightMargin(crm)
        canvas_nuis.SetBottomMargin(cbm)
        canvas_nuis.SetTopMargin(ctm)

        lat.DrawLatex(0.10, 0.92, '#bf{CMS} #it{Preliminary}')
        lat.DrawLatex(0.71 +(0.1-crm), 0.92, '35.9 fb^{-1} (13 TeV)')
        line.DrawLine(0., ycen, len(params), ycen)
        line.DrawLine(0., ymax, len(params), ymax)
        line.DrawLine(0., ymin, len(params), ymin)
        line.SetLineStyle(2); line.SetLineColor(ROOT.kRed)
        line.DrawLine(0., yhalfu, len(params), yhalfu)
        line.DrawLine(0., yhalfd, len(params), yhalfd)
        line.SetLineStyle(3)
        line.DrawLine(0., 1., len(params), 1.)
        line.DrawLine(0.,-1., len(params),-1.)
        hist_fit_s.Draw("PE1 same") ## draw again over the lines

        canvas_nuis.SetGridx()
        canvas_nuis.RedrawAxis()
        canvas_nuis.RedrawAxis('g')
        # leg=ROOT.TLegend(0.6,0.7,0.89,0.89)
        # leg.SetFillColor(0)
        # leg.SetTextFont(42)
        # leg.AddEntry(hist_fit_s,"S+B fit"   ,"EPL")
        # leg.Draw()
        fout.WriteTObject(canvas_nuis)

        for ext in ['png', 'pdf']:
            canvas_nuis.SaveAs("{od}/nuisances_{ps}_{suff}.{ext}".format(od=options.outdir, ps=uniquestr, suff=options.suffix, ext=ext))
        if len(params) >= 10:
            canv_constraints = ROOT.TCanvas('foobar', '', 800, 800)
            ## hist_fit_1d  .Scale(1./len(params))
            hist_fit_1d.GetXaxis().SetTitle('pulls and constraints')
            hist_fit_1d.GetYaxis().SetTitle('# of PDF variations')
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
            lat.DrawLatex(0.10, 0.92, '#bf{CMS} #it{Preliminary}')
            lat.DrawLatex(0.61 , 0.92, '35.9 fb^{-1} (13 TeV)')
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
            leg1.AddEntry(ff, '#splitline{{fit to pulls}}{{#scale[0.5]{{(#hat{{#mu}} = {mu:.2f}, #sigma = {si:.2f} }} }}'.format(mu=ff.GetParameter(1), si=ff.GetParameter(2)), 'l')
            leg1.AddEntry(hist_fit_1d_e, '#splitline{{constraints}}{{ #scale[0.5]{{(mean = {a:.2f} }} }}'.format(a=hist_fit_1d_e.GetMean()), 'pl')
            leg1.Draw('same')

            ##hist_fit_1d.GetYaxis().SetRangeUser(0., len(params)/2.)
            #hist_fit_1d_e.Draw(' same')
            for ext in ['png', 'pdf']:
                canv_constraints.SaveAs("{od}/nuisances_{ps}_{suff}_pulls1D.{ext}".format(od=options.outdir, ps=uniquestr, suff=options.suffix, ext=ext))

   

