# usage: python make2DEffSystVariations.py electrons_triggerSF_covariance.root effsyst.root --pdir plots
import os, sys
import numpy as np
import ROOT
ROOT.gROOT.SetBatch(True)

sys.path.append(os.getcwd() + "/plotUtils/")
from utility import *


if __name__ == "__main__":
    from optparse import OptionParser
    parser = OptionParser(usage="%prog [options] covariancefile.root systfile.root ")
    parser.add_option("-c", "--covTH3", dest="covariance", type="string", default="hist_ErfCovMatrix_vs_eta_data", help="TH3D histogram containing eta as x axis and 3x3 cov matrix as other dimensions")
    parser.add_option("-p", "--parTH2", dest="parameters", type="string", default="hist_ErfParam_vs_eta_data", help="TH2D histogram containing eta as x-axis and 3 parameters of Erf function as y-axis")
    parser.add_option("-s", "--suffix", dest="suffix",     type="string", default=None, help="suffix for the ROOT file and plots")
    parser.add_option(      "--pdir",   dest="printDir",   type="string", default=None, help="dir where to put the ROOT file and plots")
    (options, args) = parser.parse_args()
    
    if len(args)<2:
        parser.print_help()
        exit(1)

    if "EtaPtCorrelatedEfficiency_cc.so" not in ROOT.gSystem.GetLibraries():
        print "Load C++ Worker"
        ROOT.gROOT.ProcessLine(".include /cvmfs/cms.cern.ch/slc6_amd64_gcc530/external/eigen/3.2.2/include")
        ROOT.gROOT.ProcessLine(".L %s/src/CMGTools/WMass/python/postprocessing/helpers/EtaPtCorrelatedEfficiency.cc+" % os.environ['CMSSW_BASE'])

    isEle = False
    if 'ele' in args[0]: isEle = True
    lepton = "electron" if isEle else "muon"

    outfilename = args[1]
    if options.suffix: outfilename = outfilename.replace('.root','_%s.root' % options.suffix)
    if isEle: outfilename = outfilename.replace('.root','_el.root')
    else    : outfilename = outfilename.replace('.root','_mu.root')
    if options.printDir: outfilename = '%s/%s' % (options.printDir,outfilename)

    # defined in plotUtils/utility.h
    createPlotDirAndCopyPhp(options.printDir)    
    
    outf = ROOT.TFile.Open(outfilename,'recreate')
    nbins_eta = 50; nbins_pt=200
    var0 = ROOT.TH2F('p0','',nbins_eta,-2.5,2.5,nbins_pt,25,45)
    var1 = var0.Clone('p1')
    var2 = var0.Clone('p2')
    systHistos = [var0,var1,var2]

    tf = ROOT.TFile.Open(args[0])
    covHisto = tf.Get(options.covariance)
    parHisto = tf.Get(options.parameters)

    systCalc = ROOT.EtaPtCorrelatedEfficiency(covHisto,parHisto)
    
    for ieta in xrange(nbins_eta):
        eta = var0.GetXaxis().GetBinCenter(ieta+1)
        for ipt in xrange(nbins_pt):
            pt = var0.GetYaxis().GetBinCenter(ipt+1)
            relSysts = np.array([0,0,0],dtype=float)
            systCalc.DoEffSyst(eta,pt,relSysts)
            for ivar in xrange(3):
                systHistos[ivar].SetBinContent(ieta+1,ipt+1,relSysts[ivar])
                #print "eta = %.2f, pt = %.2f, syst = %.3f" % (eta,pt,relSysts[ivar])

    outf.cd()
    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetPalette(ROOT.kRainBow)
    canv = ROOT.TCanvas("canv","",900,700)
    basename = os.path.basename(outfilename)
    basename = os.path.splitext(basename)[0]

    for ivar in xrange(3):

        # systHistos[ivar].SetTitle('nuisance for Erf p%d' % ivar)
        # systHistos[ivar].GetXaxis().SetTitle("lepton #eta")
        # systHistos[ivar].GetYaxis().SetTitle("lepton p_{T} (GeV)")
        # if 'ele' in args[0]:
        #     systHistos[ivar].GetYaxis().SetRangeUser(30,45)
        # systHistos[ivar].Draw('COLZ0')
        # if ivar==0:
        #     systHistos[ivar].GetZaxis().SetRangeUser(-0.002,0.002)
        # else:
        #     systHistos[ivar].GetZaxis().SetRangeUser(-0.01,0.01)
        # for ext in ['pdf', 'png']:
        #     canv.SaveAs('{odir}/{name}.{ext}'.format(odir='.' if not options.printDir else options.printDir,
        #                                              name='%s_p%d%s' % (basename,ivar,"_"+options.suffix if options.suffix else ''),
        #                                              ext=ext))

        systHistos[ivar].SetTitle('nuisance for Erf p%d' % ivar)
        xaxisTitle = "%s #eta" % lepton
        yaxisTitle = "%s p_{T} [GeV]" % lepton
        if isEle: yaxisTitle += "::30,45"
        zaxisTitle = "variation / nominal"
        if ivar==0: zaxisTitle += "::-0.004,0.001" # 0.002,0.002
        else      : zaxisTitle += "::-0.01,0.01"
        canvasName ='%s_p%d%s' % (basename,ivar,("_"+options.suffix) if options.suffix else '')

        # defined in plotUtils/utility.h
        drawCorrelationPlot(systHistos[ivar],
                            xaxisTitle, yaxisTitle, zaxisTitle,
                            canvasName,
                            "ForceTitle",
                            '.' if not options.printDir else options.printDir,
                            1,1,False,False,False,1,leftMargin=0.14,rightMargin=0.22,
                            palette=ROOT.kRainBow,passCanvas=canv)



        systHistos[ivar].Write()
    outf.Close()

