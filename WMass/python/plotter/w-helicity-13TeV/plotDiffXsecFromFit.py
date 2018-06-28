#!/usr/bin/env python

import ROOT, os
from array import array
from make_diff_xsec_cards import getArrayParsingString
from make_diff_xsec_cards import getArrayBinNumberFromValue

import sys
#sys.path.append(os.environ['CMSSW_BASE']+"/src/CMGTools/WMass/python/plotter/")
sys.path.append(os.getcwd() + "/plotUtils/")
from utility import *

ROOT.gStyle.SetOptStat(0)
#ROOT.gStyle.SetPadRightMargin(0.13)
ROOT.gROOT.SetBatch()

if __name__ == "__main__":

    from optparse import OptionParser
    parser = OptionParser(usage="%prog [options] fitresults.root")
    parser.add_option('-o','--outdir', dest='outdir', default='.', type='string', help='output directory to save things')
    parser.add_option('-c','--channel', dest='channel', default='el', type='string', help='Channel (el, mu)')
    parser.add_option('-C','--charge', dest='charge', default='plus,minus', type='string', help='Charges to consider')
    parser.add_option('-b','--etaPtbinning', dest='etaPtbinning', default='[-2.5,-1.566,-1.4442,0,1.4442,1.566,2.5]*[30,35,40,45]', type='string', help='eta-pt binning for templates (will have to implement reading it from file). Use -b file=<name> to read binning from file <name>')
    (options, args) = parser.parse_args()

    if len(args) < 1:
        parser.print_usage()
        quit()

    channel = options.channel
    if channel not in ["el","mu"]:
        print "Error: unknown channel %s (select 'el' or 'mu')" % channel
        quit()

    outname = options.outdir
    addStringToEnd(outname,"/",notAddIfEndswithMatch=True)
    createPlotDirAndCopyPhp(outname)

    if options.etaPtbinning.startswith("file="):
        etaPtbinningFile = options.etaPtbinning.replace("file=","")
        with open(etaPtbinningFile) as f:
            content = f.readlines()
        for x in content:
            tmpbinning = str(x).strip() #if not str(x).startswith("#")
        etabinning = tmpbinning.split('*')[0]    # this is like [a,b,c,...], and is of type string. We nedd to get an array  
        ptbinning  = tmpbinning.split('*')[1]
    else:
        etabinning = options.etaPtbinning.split('*')[0]    # this is like [a,b,c,...], and is of type string. We nedd to get an array  
        ptbinning  = options.etaPtbinning.split('*')[1]
    etabinning = getArrayParsingString(etabinning,makeFloat=True)
    ptbinning  = getArrayParsingString(ptbinning,makeFloat=True)
    binning = [len(etabinning)-1, etabinning, len(ptbinning)-1, ptbinning] 
    #print binning

    lepton = "electron" if channel == "el" else " muon"
    charges = options.charge.split(',')

    for charge in charges:

        infile = ROOT.TFile(args[0], 'read')
        print ""
        print "==> RUNNING FOR CHARGE ",charge
        chs = '+' if charge == 'plus' else '-' 
        adjustSettings_CMS_lumi()

        hmu = ROOT.TH2F("hmu",'#mu: W{chs}'.format(chs=chs),len(etabinning)-1, array('d',etabinning), len(ptbinning)-1, array('d',ptbinning))
        hmu_err = ROOT.TH2F("hmu_err",'#mu uncertainty: W{chs}'.format(chs=chs),len(etabinning)-1, array('d',etabinning), len(ptbinning)-1, array('d',ptbinning))

        h_pmaskedexpnorm_mu = ROOT.TH2F("h_pmaskedexpnorm_mu",'#mu pmaskedexpnorm: W{chs}'.format(chs=chs),
                                        len(etabinning)-1, array('d',etabinning), len(ptbinning)-1, array('d',ptbinning))
        h_pmaskedexpnorm_mu_err = ROOT.TH2F("h_pmaskedexpnorm_mu_err",'pmaskedexpnorm #mu uncertainty: W{chs}'.format(chs=chs),
                                            len(etabinning)-1, array('d',etabinning), len(ptbinning)-1, array('d',ptbinning))

        h_pmaskedexp_mu = ROOT.TH2F("h_pmaskedexp_mu",'#mu pmaskedexp: W{chs}'.format(chs=chs),
                                    len(etabinning)-1, array('d',etabinning), len(ptbinning)-1, array('d',ptbinning))
        h_pmaskedexp_mu_err = ROOT.TH2F("h_pmaskedexp_mu_err",'pmaskedexp #mu uncertainty: W{chs}'.format(chs=chs),
                                        len(etabinning)-1, array('d',etabinning), len(ptbinning)-1, array('d',ptbinning))

        tree = infile.Get("fitresults")
        entries = tree.GetEntries()
        leaves  = tree.GetListOfLeaves()
        for entry in range(entries):  # actually it is just 1 entry 
            tree.GetEntry(entry)
            for p in leaves:
                name = p.GetName()
                if not "W{ch}".format(ch=charge) in name: continue
                if any(x in name for x in ["_minosup", "_minosdown", "_gen"]): continue
                # name should be like Wplus_el_ieta_4_ipt_1_Wplus_el_group_6_mu
                tokens = name.split('_')
                for i,tkn in enumerate(tokens):
                    if tkn == "ieta": etabinIndex = int(tokens[i + 1])
                    if tkn == "ipt": ptbinIndex = int(tokens[i + 1])
                    # eta and pt index start from 0, the corresponding histogram bin number is bin+1

                print "%s = %f" % (name, p.GetValue())
                # could ask whether name.endswith("_mu"), but this would depend on the argument passed to --POIMode in text2tf.py
                if "pmaskedexpnorm" in name: 
                    if name.endswith("_err"): h_pmaskedexpnorm_mu_err.SetBinContent(etabinIndex+1, ptbinIndex+1,p.GetValue())
                    else:                     h_pmaskedexpnorm_mu    .SetBinContent(etabinIndex+1, ptbinIndex+1,p.GetValue())
                elif "pmaskedexp" in name:
                    if name.endswith("_err"): h_pmaskedexp_mu_err.SetBinContent(etabinIndex+1, ptbinIndex+1,p.GetValue())
                    else:                     h_pmaskedexp_mu    .SetBinContent(etabinIndex+1, ptbinIndex+1,p.GetValue())
                else:
                    if name.endswith("_err"): hmu_err.SetBinContent(etabinIndex+1, ptbinIndex+1,p.GetValue())
                    else:                     hmu    .SetBinContent(etabinIndex+1, ptbinIndex+1,p.GetValue())
                
        xaxisTitle = '%s #eta' % lepton
        yaxisTitle = '%s p_{T} [GeV]' % lepton
        #zaxisTitle = "#mu::%.3g,%.3g" % (hmu.GetMinimum(), hmu.GetMaximum())

        zaxisTitle = "#mu::0.99,1.01"
        drawCorrelationPlot(hmu, 
                            xaxisTitle, yaxisTitle, zaxisTitle, 
                            hmu.GetName(),
                            "ForceTitle",outname,1,1,False,False,False,1)

        zaxisTitle = "uncertainty on #mu::%.3g,%.3g" % (hmu_err.GetMinimum(), hmu_err.GetMaximum())
        drawCorrelationPlot(hmu_err, 
                            xaxisTitle, yaxisTitle, zaxisTitle, 
                            hmu_err.GetName(),
                            "ForceTitle",outname,1,1,False,False,False,1)

        zaxisTitle = "d#sigma/d#etadp_{T} / #sigma_{tot}::%.3g,%.3g" % (h_pmaskedexpnorm_mu.GetMinimum(), h_pmaskedexpnorm_mu.GetMaximum())
        drawCorrelationPlot(h_pmaskedexpnorm_mu, 
                            xaxisTitle, yaxisTitle, zaxisTitle, 
                            h_pmaskedexpnorm_mu.GetName(),
                            "ForceTitle",outname,1,1,False,False,False,1)

        zaxisTitle = "uncertainty on d#sigma/d#etadp_{T} / #sigma_{tot}::%.3g,%.3g" % (h_pmaskedexpnorm_mu_err.GetMinimum(), h_pmaskedexpnorm_mu_err.GetMaximum())
        drawCorrelationPlot(h_pmaskedexpnorm_mu_err, 
                            xaxisTitle, yaxisTitle, zaxisTitle, 
                            h_pmaskedexpnorm_mu_err.GetName(),
                            "ForceTitle",outname,1,1,False,False,False,1)

        zaxisTitle = "d#sigma/d#etadp_{T} [pb]::%.3g,%.3g" % (h_pmaskedexp_mu.GetMinimum(), h_pmaskedexp_mu.GetMaximum())
        drawCorrelationPlot(h_pmaskedexp_mu, 
                            xaxisTitle, yaxisTitle, zaxisTitle, 
                            h_pmaskedexp_mu.GetName(),
                            "ForceTitle",outname,1,1,False,False,False,1)

        zaxisTitle = "uncertainty on d#sigma/d#etadp_{T} [pb]::%.3g,%.3g" % (h_pmaskedexp_mu_err.GetMinimum(), h_pmaskedexp_mu_err.GetMaximum())
        drawCorrelationPlot(h_pmaskedexp_mu_err, 
                            xaxisTitle, yaxisTitle, zaxisTitle, 
                            h_pmaskedexp_mu_err.GetName(),
                            "ForceTitle",outname,1,1,False,False,False,1)

        h_pmaskedexp_mu_relErr = h_pmaskedexp_mu_err.Clone("h_pmaskedexp_mu_relErr")
        h_pmaskedexp_mu_relErr.Divide(h_pmaskedexp_mu)
        h_pmaskedexp_mu_relErr.SetTitle('W{chs}'.format(chs=chs))

        #zaxisTitle = "relative uncertainty on d#sigma/d#etadp_{T}::%.3g,%.3g" % (h_pmaskedexp_mu_relErr.GetMinimum(), h_pmaskedexp_mu_relErr.GetMaximum())
        zaxisTitle = "relative uncertainty on d#sigma/d#etadp_{T}::0.0,0.1" #% (h_pmaskedexp_mu_relErr.GetMinimum(), h_pmaskedexp_mu_relErr.GetMaximum())
        drawCorrelationPlot(h_pmaskedexp_mu_relErr, 
                            xaxisTitle, yaxisTitle, zaxisTitle, 
                            h_pmaskedexp_mu_relErr.GetName(),
                            "ForceTitle",outname,1,1,False,False,False,1)

        infile.Close()
