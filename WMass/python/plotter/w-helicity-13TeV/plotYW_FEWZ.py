import ROOT, datetime, array, os, math, re
ROOT.gROOT.SetBatch(True)
ROOT.PyConfig.IgnoreCommandLineOptions = True
import numpy as np

from plotYWCompatibility import makeCanvas,makePads
from convertFEWZ import *

# hardcoded binwidth !!! Nice, eh?
bwidth = 0.25

def getRatios(mcvals,graph_data,rationame,ybincenters,asymmetry=False):
    ratios = []
    # this assumes the same binning
    for iy in range(graph_data.GetN()):
        x=ROOT.Double(0); yfit = ROOT.Double(0); 
        graph_data.GetPoint(iy,x,yfit)
        yfit_err = graph_data.GetErrorY(iy)
        #print "yfit = ",yfit," +/- ",yfit_err
        if not asymmetry:
            ratio = mcvals[iy][1]/yfit
            ## consider only the error on MC, since we show the error bar of data in the same ratio plot
            ratio_errUp = ratio*mcvals[iy][2]/mcvals[iy][1]
            ratio_errDn = ratio*mcvals[iy][3]/mcvals[iy][1]
        else:
            ratio = mcvals[iy][1] - yfit
            ratio_errUp = mcvals[iy][2]
            ratio_errDn = mcvals[iy][3]
        ratios.append((mcvals[iy][1], ratio, ratio_errUp, ratio_errDn))        
    print ratios

    ## make the TGraph out of them
    ratioGr = ROOT.TGraphAsymmErrors()
    ratioGr.SetName(rationame)
    for ir,r in enumerate(ratios):
        ratioGr.SetPoint(ir,ybincenters[ir],r[1])
        ratioGr.SetPointError(ir,0.5*bwidth,0.5*bwidth,r[2],r[3])
    return ratioGr

def makeFullLegend(graphs,titles,styles,ncols=1,xmin=0.25,ymin=0.7,xmax=0.9,ymax=0.9,textSize=0.05):
    leg = ROOT.TLegend(xmin,ymin,xmax,ymax)
    leg.SetNColumns(ncols)
    leg.SetFillStyle(0)
    leg.SetBorderSize(0)
    leg.SetTextFont(42)
    leg.SetTextSize(textSize)
    for i in range(len(graphs)):        
        # try to align item for 'data' in canvas, TLegend screws the position based on the text size and length
        if len(graphs)==2 and "FEWZ" in titles[i]:
            if "CT18" in titles[i]:
                leg.AddEntry(graphs[i], titles[i]+"            ", styles[i]) 
            else:
                leg.AddEntry(graphs[i], titles[i]+"       ", styles[i]) 
        else:
            leg.AddEntry(graphs[i], titles[i], styles[i])
    return leg

def plotOne(mcvals,mcname,color,fillstyle,ybincenters):
    graph = ROOT.TGraphAsymmErrors()
    graph.SetName(mcname)
    for point,yvals in enumerate(mcvals):
        if point==len(ybincenters): break
        y,xsec,errUp,errDn = yvals
        graph.SetPoint(point,ybincenters[point],xsec)
        graph.SetPointError(point,0.5*bwidth,0.5*bwidth,errDn,errUp)
        
    graph.SetLineColor(0)
    graph.SetFillStyle(fillstyle)
    #graph.SetFillColorAlpha(color,0.30)
    graph.SetFillColor(color)
    graph.Draw('E2')
    return graph

def plotOneRatio(ratioGraphs,mcname,color,fillstyle,isBottomPlot=False,asymmetry=False):

    for iratio,ratio in enumerate(ratioGraphs):
        if iratio==0:
            #ratio.SetFillColorAlpha(color,0.30)
            ratio.SetFillStyle(fillstyle)
            ratio.SetFillColor(color)
            ratio.SetLineColor(0)
            ratio.Draw('A2')
            if isBottomPlot:
                ratio.GetXaxis().SetLabelSize(0.15)
            else:
                ratio.GetXaxis().SetLabelSize(0)
            ratio.GetXaxis().SetRangeUser(0,2.5)
            if not asymmetry:
                ratio.GetYaxis().SetRangeUser(0.8,1.2)
                ratio.GetYaxis().SetTitle('Theory/Data')
            else:
                ratio.GetYaxis().SetRangeUser(-0.1,0.1)
                ratio.GetYaxis().SetTitle('Theory-Data')
            ratio.GetYaxis().SetNdivisions(505,ROOT.kFALSE)
            ratio.GetYaxis().SetLabelSize(0.12)
            ratio.GetYaxis().SetTitleSize(0.16)
            ratio.GetYaxis().SetTitleOffset(0.5)
            ratio.GetYaxis().CenterTitle()
            ratio.GetYaxis().SetDecimals()
        else:
            ratio.SetLineColor(ROOT.kBlack)
            ratio.SetMarkerColor(ROOT.kBlack)
            ratio.SetMarkerSize(2)
            ratio.Draw('pe')

    leg = makeFullLegend(ratioGraphs[:2],
                          [mcname,'data'],
                          ['f','pe'],ymin=0.70,xmax=0.7,ncols=2,textSize=0.09)
                                                     #0.75               #0.08
    leg.Draw()
    leg.SetName("leg"+mcname)

    line = ROOT.TF1("horiz_line"+mcname,"0" if asymmetry else '1',0.0,3.0);
    line.SetLineColor(ROOT.kRed);
    line.SetLineWidth(2);
    line.SetLineStyle(ROOT.kDashed);
    line.Draw("Lsame");

    ## assign it to a global variable so it's not deleted
    global legend_
    legend_ = leg
    global line_
    line_ = line
    return (leg,line)

def plotMCaNLO(values,ratios,fewzvals_nnpdf31,fewzvals_ct18,plotname,outdir):

    of = ROOT.TFile.Open('{od}/{plot}_fewz.root'.format(od=outdir, plot=plotname),'recreate')

    c2 = makeCanvas()
    ROOT.gStyle.SetHatchesLineWidth(1)
    ROOT.gStyle.SetHatchesSpacing(0.3)

    ch = '+' if 'plus' in plotname else '-' if 'minus' in plotname else ' '
    asymmetry = (ch==' ')

    ## 3 pads
    subpadsYLow = [0.52, 0.36, 0.20, 0.04]
    subpadsYUp =  [0.98, 0.52, 0.36, 0.20]
    #subpadsYLow = [0.52, 0.41, 0.30, 0.19, 0.08]
    #subpadsYUp =  [0.98, 0.52, 0.41, 0.30, 0.19]
    pads = makePads(subpadsYLow,subpadsYUp)
    ## these are needed to make them persistent
    legends = []
    lines = []

    ## make the values pad
    c2.cd()
    pads[0].Draw(); pads[0].cd()
    graphs = values.GetListOfGraphs()
    yranges = {'+': (2200,3600),
               '-': (1600,3500),
               ' ': (-0.05,0.30)
               }
    ybincenters = []
    for i,gr in enumerate(graphs): 
        if i==0:
            # for some reason, this info is lost
            gr.SetFillStyle(1001)
            #gr.SetFillColorAlpha(ROOT.kOrange-2,0.30)
            gr.SetFillColor(ROOT.kOrange-2)
            gr.SetLineColor(0)
            gr.Draw('AE2')
            gr.GetXaxis().SetRangeUser(0,2.5)
            gr.GetYaxis().SetRangeUser(yranges[ch][0],yranges[ch][1])
            gr.GetXaxis().SetLabelSize(0)
            gr.GetYaxis().SetTitleSize(0.06)
            gr.GetYaxis().SetLabelSize(0.06)
            gr.GetYaxis().SetTitleOffset(1.4)
            gr.GetYaxis().CenterTitle()
            for ip in range(gr.GetN()):
                x=ROOT.Double(0); dummy = ROOT.Double(0);
                gr.GetPoint(ip,x,dummy)
                ybincenters.append(x)
            if ch!=' ':
                gr.GetYaxis().SetTitle('d#sigma / d|y_{W}| (pb)')
            else:
                gr.GetYaxis().SetTitle('Charge asymmetry')
        elif i==1: 
            ## this is data, draw it on top
            pass
        elif i==2: gr.Draw('L2')
        else: print "BUH !"
        

    ## superimpose FEWZ + NNPDF3.1
    fewz_nnpdf31_Gr = plotOne(fewzvals_nnpdf31,'fewz_nnpdf31',ROOT.kMagenta+1,3385,ybincenters)

    ## superimpose FEWZ + CT18
    fewz_ct18_Gr = plotOne(fewzvals_ct18,'fewz_ct18',ROOT.kAzure+10,3358,ybincenters)

    ## superimpose the DATA on top
    graphs[1].SetLineColor(ROOT.kBlack)
    graphs[1].SetMarkerColor(ROOT.kBlack)
    graphs[1].SetMarkerSize(2)
    graphs[1].SetLineWidth(1)
    graphs[1].Draw('PE')

    leg0 = makeFullLegend([graphs[0],graphs[2],fewz_nnpdf31_Gr,fewz_ct18_Gr,graphs[1]],
                          ['MC@NLO NNPDF3.0','MC@NLO* NNPDF3.0','FEWZ NNPDF3.1','FEWZ CT18','data'],
                          ['f','l','f','f','pl'],ncols=2,textSize=0.03)
    leg0.Draw()

    ## save in a ROOT file
    of.cd()
    graphs[0].SetName('mcanlo')
    graphs[1].SetName('data')
    graphs[2].SetName('mcanlo_rwgt')
    for ig in range(3): 
        graphs[ig].Write()
    fewz_nnpdf31_Gr.Write()
    fewz_ct18_Gr.Write()
    of.Close()

    ## make the nominal ratio pad    
    c2.cd()
    pads[1].Draw(); pads[1].cd();
    mcnlo_ratios = ratios.GetListOfGraphs()
    leg1,line1 = plotOneRatio(mcnlo_ratios,'MC@NLO NNPDF3.0',ROOT.kOrange-2,1001,asymmetry=asymmetry)
    legends.append(leg1); lines.append(line1)

    # now the FEWZ + NNPDF3.1 ratio pad
    c2.cd()
    pads[2].Draw(); pads[2].cd();
    fewz_nnpdf31_ToData = getRatios(fewzvals_nnpdf31,graphs[1],'ratio_fewz_nnpdf31',ybincenters,asymmetry=asymmetry)
    leg2,line2 = plotOneRatio([fewz_nnpdf31_ToData,mcnlo_ratios[1]],'FEWZ NNPDF3.1',ROOT.kMagenta+1,3385,asymmetry=asymmetry)
    legends.append(leg2); lines.append(line2)

    # now the FEWZ + CT18 ratio pad 
    c2.cd()
    pads[3].Draw(); pads[3].cd();
    fewz_ct18_ToData = getRatios(fewzvals_ct18,graphs[1],'ratio_fewz_ct18',ybincenters,asymmetry=asymmetry)
    leg3,line3 = plotOneRatio([fewz_ct18_ToData,mcnlo_ratios[1]],'FEWZ CT18',ROOT.kAzure+10,3358,isBottomPlot=True,asymmetry=asymmetry)
    legends.append(leg3); lines.append(line3)

    c2.cd()
    lat = ROOT.TLatex()
    lat.SetNDC(); lat.SetTextFont(42)
    lat.DrawLatex(0.2, 0.97, '#bf{CMS}') #it{Preliminary}')
    lat.DrawLatex(0.60, 0.97, '35.9 fb^{-1} (13 TeV)')
    lat.DrawLatex(0.25, 0.57,  'W^{{{ch}}} #rightarrow l^{{{ch}}}{nu}'.format(ch=ch,nu="#bar{#nu}" if ch=='-' else "#nu"))
    lat2 = ROOT.TLatex()
    lat2.SetNDC(); lat2.SetTextFont(42);  lat2.SetTextSize(0.04);
    lat2.DrawLatex(0.90, 0.012, '|y_{W}|')


    for ext in ['png', 'pdf', 'C']:
        c2.SaveAs('{od}/{plot}_fewz.{ext}'.format(od=outdir, plot=plotname, ext=ext))



if __name__ == "__main__":

    from optparse import OptionParser
    parser = OptionParser(usage='%prog ybinfile workspace.root toys.root [options] ')
    parser.add_option('-d', '--dir'         , dest='dir'      , default='.'           , type='string', help='input dir where to take the fit graphs and  outdput directory to save the plots')
    parser.add_option(      '--pdfonly'     , dest='pdfonly'  , default=False         , action='store_true',   help='if given, do not include alphaS and QCD scales in the theory bands')
    (options, args) = parser.parse_args()

    if not os.path.isdir(options.dir):
        print "Need the input directory where to take the MC@NLO graphs"
        exit(1)
    else:
        theDir = options.dir

    plots = {'p':   'genAbsYUnpolarizedsumxsec_pdfs_plusfloatingPOIs_hessian_bbb1_syst1_data_lep_hessian',
             'm':   'genAbsYUnpolarizedsumxsec_pdfs_minusfloatingPOIs_hessian_bbb1_syst1_data_lep_hessian',
             'a': 'genAbsYUnpolarizedasymmetry_pdfs_asymmetryfloatingPOIs_hessian_bbb1_syst1_data_lep_hessian',
             }

    if options.pdfonly:
        print "===> WARNING! PDF ONLY UNCERTAINTIES (NOT QCD SCALES) <==="
        print "Are you taking the MC@NLO from the correct directory????? "

    print "now converting FEWZ..."
    fewz_nnpdf31 = getForPDF('nnpdf31', pdfsonly=options.pdfonly)
    fewz_ct18    = getForPDF('ct18'   , pdfsonly=options.pdfonly)

    print "and now it is THE time to plot XSECS!"
    for charge in ['p','m','a']:
        plot = plots[charge]
        infile = ROOT.TFile.Open('{d}/{p}.root'.format(d=theDir,p=plot))
        nominal_values = infile.Get('values')
        nominal_ratios = infile.Get('ratios')

        fewz_nnpdf31_sc = np.array(fewz_nnpdf31[charge])
        fewz_ct18_sc = np.array(fewz_ct18[charge])
        if charge!='a':
            # rescale by the Y bin width
            fewz_nnpdf31_sc[:,1:] *= 1./bwidth
            fewz_ct18_sc[:,1:] *= 1./bwidth

            ## get the nominal MC and data from our beloved graphs
        canv = plotMCaNLO(nominal_values,nominal_ratios,fewz_nnpdf31_sc,fewz_ct18_sc,plot,theDir)

    

    
        
