#!/usr/bin/env python                                                       

# to add smooth pt scales in a file that contains the standard version (used to fix the extreme pt bins)
 
#from shutil import copyfile
import re, sys, os, os.path, ROOT, copy, math
from array import array
#import numpy as np

from w_helicity_13TeV.utilities import util
utilities = util()

from w_helicity_13TeV.make_diff_xsec_cards import getXYBinsFromGlobalBin
from w_helicity_13TeV.make_diff_xsec_cards import getArrayParsingString
from w_helicity_13TeV.make_diff_xsec_cards import getDiffXsecBinning
from w_helicity_13TeV.make_diff_xsec_cards import templateBinning
from w_helicity_13TeV.make_diff_xsec_cards import get_ieta_ipt_from_process_name
from w_helicity_13TeV.make_diff_xsec_cards import get_ieta_from_process_name
from w_helicity_13TeV.make_diff_xsec_cards import get_ipt_from_process_name

from w_helicity_13TeV.templateRolling import roll1Dto2D, dressed2D
from w_helicity_13TeV.rollingFunctions import unroll2Dto1D

from w_helicity_13TeV.mergeCardComponentsAbsY import mirrorShape
from w_helicity_13TeV.mergeCardComponentsAbsY import putUncorrelatedFakes
from w_helicity_13TeV.mergeCardComponentsAbsY import putEffSystHistos # for electron L1-prefire uncertainty
from w_helicity_13TeV.mergeCardComponentsAbsY import addZOutOfAccPrefireSyst # for ele L1-prefire uncertainty on Z


# experimental function for momentum scales
def addSmoothLepScaleSyst(infile,regexp,charge, isMu, outdir=None, 
                          uncorrelateByCharge=False,
                          uncorrelateByEtaSide=False,
                          getExtremePtBinsFromOldScale=True):

    indir = os.path.dirname(os.path.abspath(infile))
    outdir = outdir if outdir != None else indir
    flav = 'mu' if isMu else 'el'

    # get eta-pt binning for reco 
    etaPtBinningVec = getDiffXsecBinning(indir+'/binningPtEta.txt', "reco")  # this get two vectors with eta and pt binning
    recoBins = templateBinning(etaPtBinningVec[0],etaPtBinningVec[1])        # this create a class to manage the binnings
    binning = [recoBins.Neta, recoBins.etaBins, recoBins.Npt, recoBins.ptBins]

    tmp_infile = ROOT.TFile(infile, 'read')
    outfile = ROOT.TFile(outdir+'/SmoothScaleSyst_{flav}_{ch}.root'.format(flav=flav,ch=charge), 'recreate')

    ## scale systematics as a function of eta
    etabins_mu = array('f',[0.0, 2.1, 2.4])

    scaleSyst_mu = ROOT.TH1F('scaleSyst_mu','',len(etabins_mu)-1,etabins_mu)
    scaleSyst_mu.SetBinContent(1,0.003)
    scaleSyst_mu.SetBinContent(2,0.010)
    
    etabins_el = array('f',[0.0, 1.0, 1.5, 2.1, 2.4])
    scaleSyst_el = ROOT.TH1F('scaleSyst_el','',len(etabins_el)-1,etabins_el)
    scaleSyst_el.SetBinContent(1,0.003)
    scaleSyst_el.SetBinContent(2,0.005)
    scaleSyst_el.SetBinContent(3,0.008)
    scaleSyst_el.SetBinContent(4,0.01)
        
    if isMu:
        scaleSyst = utilities.getExclusiveBinnedSyst(scaleSyst_mu)
        etabins = etabins_mu
    else:
        scaleSyst = utilities.getExclusiveBinnedSyst(scaleSyst_el)
        etabins = etabins_el

    
    for k in tmp_infile.GetListOfKeys():
        tmp_name = k.GetName()
        ## don't reweight any histos that don't match the regexp
        if not re.match(regexp, tmp_name): continue
        ## don't reweight any histos that are already variations of something else
        if 'Up' in tmp_name or 'Down' in tmp_name: continue
        process = tmp_name.split('_')[1]

        ## now should be left with only the ones we are interested in
        print 'reweighting smooth lepscale syst for process', tmp_name
        
        tmp_nominal = tmp_infile.Get(tmp_name)
        tmp_nominal_2d = dressed2D(tmp_nominal,binning, tmp_name+'backrolled')
        n_ptbins = tmp_nominal_2d.GetNbinsY()

        etasides = [-1, 1] if uncorrelateByEtaSide else [1]

        for side in etasides:
            sideTag = "" if not uncorrelateByEtaSide else "etaside{s}".format(s="P" if side > 0 else "M")
            for systBin in xrange(len(etabins)-1):
                etasyst = scaleSyst.GetXaxis().GetBinCenter(systBin+1)
                for shift_dir in ['Up','Down']:
                    outname_2d = tmp_nominal_2d.GetName().replace('backrolled','')+'_smooth{lep}scale{idx}{ch}{st}{shiftdir}'.format(lep=flav,idx=systBin,shiftdir=shift_dir,ch=charge if uncorrelateByCharge else "",st=sideTag)
                    tmp_scaledHisto = copy.deepcopy(tmp_nominal_2d.Clone(outname_2d))
                    tmp_scaledHisto.Reset()

                    outname_2d_old = ""
                    tmp_oldscale = None
                    if getExtremePtBinsFromOldScale: 
                        outname_2d_old = outname_2d.replace("smoothmuscale","CMS_Wmu_muscale") if isMu else outname_2d.replace("smoothelscale","CMS_We_elescale")
                        tmp_oldscale = tmp_infile.Get(outname_2d_old)
                        if not tmp_oldscale:
                            print "Error in addSmoothLepScaleSyst(): could not get histogram %s with old scales to support %s. Abort" % (outname_2d_old, outname_2d)
                            quit()


                    ## loop over all eta bins of the 2d histogram 
                    for ieta in range(1,tmp_nominal_2d.GetNbinsX()+1):
                        eta = tmp_nominal_2d.GetXaxis().GetBinCenter(ieta)
                        scale_syst = scaleSyst.GetBinContent(scaleSyst.GetXaxis().FindFixBin(etasyst))
                        for ipt in range(1,tmp_nominal_2d.GetNbinsY()+1):                            
                            if abs(eta)>etabins[systBin]:
                                if (uncorrelateByEtaSide and eta*float(side) < 0):
                                    # we want to uncorrelate, but this eta is the other side
                                    tmp_scaledHisto.SetBinContent(ieta, ipt, tmp_nominal_2d.GetBinContent(ieta,ipt))
                                else:
                                    ## assume uniform distribution within a bin
                                    if shift_dir == "Up":
                                        # for outliers take extreme bin from old scales and also 2nd based on shift
                                        if ipt == 1 or ("outliers" in tmp_name and any(ipt == int(x) for x in [2,tmp_nominal_2d.GetNbinsY()])):
                                            if getExtremePtBinsFromOldScale: 
                                                # can avoid unrolling, as TH1 are unrolled vs eta, so the first and last N=Neta bins of tmp_oldscale are already what we need
                                                #print "Getting pt bin %d from old scale" % ipt
                                                offsetEta = recoBins.Neta * (ipt - 1)
                                                tmp_scaledHisto.SetBinContent(ieta,ipt,tmp_oldscale.GetBinContent(ieta + offsetEta))
                                            else:
                                                continue
                                        else:
                                            # assuming 1% scale variations and a bin with pt in [25, 26] GeV 
                                            # (and previous bin from 24 to 25), only the events in the previous bin 
                                            # for which (1+0.01)*pt' > 25 will move to the bin in [25, 26], so this pt' 
                                            # is 25/1.01=24.7525 GeV, which corresponds to a fraction of events in 
                                            # previous bin of 24.75%
                                            # so, if ptl is the low-pt edge, and x is the scale variation (scaleSyst_mu/el)
                                            # we can write DeltaN_low = |ptl - ptl/(1+x)|/binWidth_low * binContent_low =
                                            # = x/(1+x) *ptl/binWidth_low * binContent_low
                                            pt      = tmp_nominal_2d.GetYaxis().GetBinUpEdge(ipt)
                                            pt_prev = tmp_nominal_2d.GetYaxis().GetBinLowEdge(ipt)
                                            nominal_val      = tmp_nominal_2d.GetBinContent(ieta,ipt)
                                            nominal_val_prev = tmp_nominal_2d.GetBinContent(ieta,ipt-1)
                                            pt_width      = tmp_nominal_2d.GetYaxis().GetBinWidth(ipt)
                                            pt_width_prev = tmp_nominal_2d.GetYaxis().GetBinWidth(ipt-1)
                                            factor = scale_syst / (1.0 + scale_syst)
                                            from_prev = factor * (pt_prev / pt_width_prev) * max(0,nominal_val_prev)
                                            to_right  = factor * (pt      / pt_width     ) * max(0,nominal_val)
                                            tmp_scaledHisto.SetBinContent(ieta,ipt,max(0,nominal_val + from_prev - to_right))
                                    else:
                                        # for outliers take extreme bin from old scales and also 2nd based on shift
                                        if ipt == tmp_nominal_2d.GetNbinsY() or ("outliers" in tmp_name and any(ipt == int(x) for x in [1,(tmp_nominal_2d.GetNbinsY()-1)])):
                                            if getExtremePtBinsFromOldScale: 
                                                # can avoid unrolling, as TH1 are unrolled vs eta, so the first and last N=Neta bins of tmp_oldscale are already what we need
                                                #print "Getting pt bin %d from old scale" % ipt
                                                offsetEta = recoBins.Neta * (ipt - 1)
                                                tmp_scaledHisto.SetBinContent(ieta,ipt,tmp_oldscale.GetBinContent(offsetEta + ieta))
                                            else:
                                                continue
                                        else:
                                            # same logic as before, but now (1-0.01)*pt' < 25 --> pt'< 25.2525 
                                            pt     = tmp_nominal_2d.GetYaxis().GetBinLowEdge(ipt)
                                            pt_seq = tmp_nominal_2d.GetYaxis().GetBinUpEdge(ipt)
                                            nominal_val     = tmp_nominal_2d.GetBinContent(ieta,ipt)
                                            nominal_val_seq = tmp_nominal_2d.GetBinContent(ieta,ipt+1)
                                            pt_width     = tmp_nominal_2d.GetYaxis().GetBinWidth(ipt)
                                            pt_width_seq = tmp_nominal_2d.GetYaxis().GetBinWidth(ipt+1)
                                            factor = scale_syst / (1.0 - scale_syst)
                                            from_seq = factor * (pt_seq / pt_width_seq) * max(0,nominal_val_seq)
                                            to_left  = factor * (pt     / pt_width    ) * max(0,nominal_val)
                                            tmp_scaledHisto.SetBinContent(ieta,ipt,max(0,nominal_val + from_seq - to_left))
                            else:
                                tmp_scaledHisto.SetBinContent(ieta, ipt, tmp_nominal_2d.GetBinContent(ieta,ipt))
                        ## since in the first bin we cannot foresee 2-neighbors migrations, let's assume same syst of bin i+1
                        if not getExtremePtBinsFromOldScale:
                            if shift_dir == "Up":
                                factor_bin2 = tmp_scaledHisto.GetBinContent(ieta,2)/tmp_nominal_2d.GetBinContent(ieta,2) if tmp_nominal_2d.GetBinContent(ieta,2) else 0
                                tmp_scaledHisto.SetBinContent(ieta,1,factor_bin2 *tmp_nominal_2d.GetBinContent(ieta,1))
                            else:
                                binLastMinus1 = tmp_nominal_2d.GetNbinsY() - 1
                                factor_binLastMinus1 = tmp_scaledHisto.GetBinContent(ieta,binLastMinus1)/tmp_nominal_2d.GetBinContent(ieta,binLastMinus1) if tmp_nominal_2d.GetBinContent(ieta,binLastMinus1) else 0
                                tmp_scaledHisto.SetBinContent(ieta,tmp_nominal_2d.GetNbinsY(),factor_binLastMinus1 *tmp_nominal_2d.GetBinContent(ieta,tmp_nominal_2d.GetNbinsY()))

                    ## re-roll the 2D to a 1D histo
                    tmp_scaledHisto_1d = unroll2Dto1D(tmp_scaledHisto, newname=tmp_scaledHisto.GetName().replace('2DROLLED',''), cropNegativeBins=False)

                    outfile.cd()
                    tmp_scaledHisto_1d.Write()
        # end loop on sides
    outfile.Close()
    print "done with the smooth lep scale variations"
    print "Output file: %s" % outfile.GetName()


# experimental function for muon momentum scales using Rochester corrections (to be tested)
# TO BE COMPLETED
def addSmoothMuonScaleSyst(infile,regexp,charge, isMu, outdir=None, 
                           uncorrelateByCharge=False,
                           uncorrelateByEtaSide=False,
                           getExtremePtBinsFromOldScale=True):

    indir = os.path.dirname(os.path.abspath(infile))
    outdir = outdir if outdir != None else indir
    flav = 'mu' if isMu else 'el'

    # get eta-pt binning for reco 
    etaPtBinningVec = getDiffXsecBinning(indir+'/binningPtEta.txt', "reco")  # this get two vectors with eta and pt binning
    recoBins = templateBinning(etaPtBinningVec[0],etaPtBinningVec[1])        # this create a class to manage the binnings
    binning = [recoBins.Neta, recoBins.etaBins, recoBins.Npt, recoBins.ptBins]

    tmp_infile = ROOT.TFile(infile, 'read')
    outfile = ROOT.TFile(outdir+'/SmoothScaleSyst_{flav}_{ch}.root'.format(flav=flav,ch=charge), 'recreate')

    ## scale systematics as a function of eta
    etabins_mu = array('f',[0.0, 2.1, 2.4])

    scaleSyst_mu = ROOT.TH1F('scaleSyst_mu','',len(etabins_mu)-1,etabins_mu)
    scaleSyst_mu.SetBinContent(1,0.003)
    scaleSyst_mu.SetBinContent(2,0.010)
    
    etabins_el = array('f',[0.0, 1.0, 1.5, 2.1, 2.4])
    scaleSyst_el = ROOT.TH1F('scaleSyst_el','',len(etabins_el)-1,etabins_el)
    scaleSyst_el.SetBinContent(1,0.003)
    scaleSyst_el.SetBinContent(2,0.005)
    scaleSyst_el.SetBinContent(3,0.008)
    scaleSyst_el.SetBinContent(4,0.01)
        
    if isMu:
        scaleSyst = utilities.getExclusiveBinnedSyst(scaleSyst_mu)
        etabins = etabins_mu
    else:
        scaleSyst = utilities.getExclusiveBinnedSyst(scaleSyst_el)
        etabins = etabins_el

    
    for k in tmp_infile.GetListOfKeys():
        tmp_name = k.GetName()
        ## don't reweight any histos that don't match the regexp
        if not re.match(regexp, tmp_name): continue
        ## don't reweight any histos that are already variations of something else
        if 'Up' in tmp_name or 'Down' in tmp_name: continue
        process = tmp_name.split('_')[1]

        ## now should be left with only the ones we are interested in
        print 'reweighting smooth lepscale syst for process', tmp_name
        
        tmp_nominal = tmp_infile.Get(tmp_name)
        tmp_nominal_2d = dressed2D(tmp_nominal,binning, tmp_name+'backrolled')
        n_ptbins = tmp_nominal_2d.GetNbinsY()

        etasides = [-1, 1] if uncorrelateByEtaSide else [1]

        for side in etasides:
            sideTag = "" if not uncorrelateByEtaSide else "etaside{s}".format(s="P" if side > 0 else "M")
            for systBin in xrange(len(etabins)-1):
                etasyst = scaleSyst.GetXaxis().GetBinCenter(systBin+1)
                for shift_dir in ['Up','Down']:
                    outname_2d = tmp_nominal_2d.GetName().replace('backrolled','')+'_smooth{lep}scale{idx}{ch}{st}{shiftdir}'.format(lep=flav,idx=systBin,shiftdir=shift_dir,ch=charge if uncorrelateByCharge else "",st=sideTag)
                    tmp_scaledHisto = copy.deepcopy(tmp_nominal_2d.Clone(outname_2d))
                    tmp_scaledHisto.Reset()

                    outname_2d_old = ""
                    tmp_oldscale = None
                    if getExtremePtBinsFromOldScale: 
                        outname_2d_old = outname_2d.replace("smoothmuscale","CMS_Wmu_muscale") if isMu else outname_2d.replace("smoothelscale","CMS_We_elescale")
                        tmp_oldscale = tmp_infile.Get(outname_2d_old)
                        if not tmp_oldscale:
                            print "Error in addSmoothLepScaleSyst(): could not get histogram %s with old scales to support %s. Abort" % (outname_2d_old, outname_2d)
                            quit()


                    ## loop over all eta bins of the 2d histogram 
                    for ieta in range(1,tmp_nominal_2d.GetNbinsX()+1):
                        eta = tmp_nominal_2d.GetXaxis().GetBinCenter(ieta)
                        scale_syst = scaleSyst.GetBinContent(scaleSyst.GetXaxis().FindFixBin(etasyst))
                        for ipt in range(1,tmp_nominal_2d.GetNbinsY()+1):                            
                            if abs(eta)>etabins[systBin]:
                                if (uncorrelateByEtaSide and eta*float(side) < 0):
                                    # we want to uncorrelate, but this eta is the other side
                                    tmp_scaledHisto.SetBinContent(ieta, ipt, tmp_nominal_2d.GetBinContent(ieta,ipt))
                                else:
                                    ## assume uniform distribution within a bin
                                    if shift_dir == "Up":
                                        # for outliers take extreme bin from old scales and also 2nd based on shift
                                        if ipt == 1 or ("outliers" in tmp_name and any(ipt == int(x) for x in [2,tmp_nominal_2d.GetNbinsY()])):
                                            if getExtremePtBinsFromOldScale: 
                                                # can avoid unrolling, as TH1 are unrolled vs eta, so the first and last N=Neta bins of tmp_oldscale are already what we need
                                                #print "Getting pt bin %d from old scale" % ipt
                                                offsetEta = recoBins.Neta * (ipt - 1)
                                                tmp_scaledHisto.SetBinContent(ieta,ipt,tmp_oldscale.GetBinContent(ieta + offsetEta))
                                            else:
                                                continue
                                        else:
                                            # assuming 1% scale variations and a bin with pt in [25, 26] GeV 
                                            # (and previous bin from 24 to 25), only the events in the previous bin 
                                            # for which (1+0.01)*pt' > 25 will move to the bin in [25, 26], so this pt' 
                                            # is 25/1.01=24.7525 GeV, which corresponds to a fraction of events in 
                                            # previous bin of 24.75%
                                            # so, if ptl is the low-pt edge, and x is the scale variation (scaleSyst_mu/el)
                                            # we can write DeltaN_low = |ptl - ptl/(1+x)|/binWidth_low * binContent_low =
                                            # = x/(1+x) *ptl/binWidth_low * binContent_low
                                            pt      = tmp_nominal_2d.GetYaxis().GetBinUpEdge(ipt)
                                            pt_prev = tmp_nominal_2d.GetYaxis().GetBinLowEdge(ipt)
                                            nominal_val      = tmp_nominal_2d.GetBinContent(ieta,ipt)
                                            nominal_val_prev = tmp_nominal_2d.GetBinContent(ieta,ipt-1)
                                            pt_width      = tmp_nominal_2d.GetYaxis().GetBinWidth(ipt)
                                            pt_width_prev = tmp_nominal_2d.GetYaxis().GetBinWidth(ipt-1)
                                            factor = scale_syst / (1.0 + scale_syst)
                                            from_prev = factor * (pt_prev / pt_width_prev) * max(0,nominal_val_prev)
                                            to_right  = factor * (pt      / pt_width     ) * max(0,nominal_val)
                                            tmp_scaledHisto.SetBinContent(ieta,ipt,max(0,nominal_val + from_prev - to_right))
                                    else:
                                        # for outliers take extreme bin from old scales and also 2nd based on shift
                                        if ipt == tmp_nominal_2d.GetNbinsY() or ("outliers" in tmp_name and any(ipt == int(x) for x in [1,(tmp_nominal_2d.GetNbinsY()-1)])):
                                            if getExtremePtBinsFromOldScale: 
                                                # can avoid unrolling, as TH1 are unrolled vs eta, so the first and last N=Neta bins of tmp_oldscale are already what we need
                                                #print "Getting pt bin %d from old scale" % ipt
                                                offsetEta = recoBins.Neta * (ipt - 1)
                                                tmp_scaledHisto.SetBinContent(ieta,ipt,tmp_oldscale.GetBinContent(offsetEta + ieta))
                                            else:
                                                continue
                                        else:
                                            # same logic as before, but now (1-0.01)*pt' < 25 --> pt'< 25.2525 
                                            pt     = tmp_nominal_2d.GetYaxis().GetBinLowEdge(ipt)
                                            pt_seq = tmp_nominal_2d.GetYaxis().GetBinUpEdge(ipt)
                                            nominal_val     = tmp_nominal_2d.GetBinContent(ieta,ipt)
                                            nominal_val_seq = tmp_nominal_2d.GetBinContent(ieta,ipt+1)
                                            pt_width     = tmp_nominal_2d.GetYaxis().GetBinWidth(ipt)
                                            pt_width_seq = tmp_nominal_2d.GetYaxis().GetBinWidth(ipt+1)
                                            factor = scale_syst / (1.0 - scale_syst)
                                            from_seq = factor * (pt_seq / pt_width_seq) * max(0,nominal_val_seq)
                                            to_left  = factor * (pt     / pt_width    ) * max(0,nominal_val)
                                            tmp_scaledHisto.SetBinContent(ieta,ipt,max(0,nominal_val + from_seq - to_left))
                            else:
                                tmp_scaledHisto.SetBinContent(ieta, ipt, tmp_nominal_2d.GetBinContent(ieta,ipt))
                        ## since in the first bin we cannot foresee 2-neighbors migrations, let's assume same syst of bin i+1
                        if not getExtremePtBinsFromOldScale:
                            if shift_dir == "Up":
                                factor_bin2 = tmp_scaledHisto.GetBinContent(ieta,2)/tmp_nominal_2d.GetBinContent(ieta,2) if tmp_nominal_2d.GetBinContent(ieta,2) else 0
                                tmp_scaledHisto.SetBinContent(ieta,1,factor_bin2 *tmp_nominal_2d.GetBinContent(ieta,1))
                            else:
                                binLastMinus1 = tmp_nominal_2d.GetNbinsY() - 1
                                factor_binLastMinus1 = tmp_scaledHisto.GetBinContent(ieta,binLastMinus1)/tmp_nominal_2d.GetBinContent(ieta,binLastMinus1) if tmp_nominal_2d.GetBinContent(ieta,binLastMinus1) else 0
                                tmp_scaledHisto.SetBinContent(ieta,tmp_nominal_2d.GetNbinsY(),factor_binLastMinus1 *tmp_nominal_2d.GetBinContent(ieta,tmp_nominal_2d.GetNbinsY()))

                    ## re-roll the 2D to a 1D histo
                    tmp_scaledHisto_1d = unroll2Dto1D(tmp_scaledHisto, newname=tmp_scaledHisto.GetName().replace('2DROLLED',''), cropNegativeBins=False)

                    outfile.cd()
                    tmp_scaledHisto_1d.Write()
        # end loop on sides
    outfile.Close()
    print "done with the smooth lep scale variations"
    print "Output file: %s" % outfile.GetName()



if __name__ == "__main__":

    from optparse import OptionParser
    parser = OptionParser(usage="%prog [options]")
    #parser.add_option("-i", "--indir",    dest="indir",  type="string", default="", help="Input folder");
    parser.add_option("-f", "--infile",   dest="infile", type="string", default="", help="Input file name");
    parser.add_option("-o", "--outdir",   dest="outdir", type="string", default="./", help="Output folder (current one as default)");
    parser.add_option("-p", "--processes",   dest="processes", type="string", default="x_Z|x_W.*|x_TauDecaysW", help="Regular expression for processes that need the pt scale systematics");
    #parser.add_option("-f", "--flavour", dest="flavour", type="string", default='el', help="Channel: either 'el' or 'mu'");
    # charge is only used for the file name, not for the match
    #parser.add_option("-c", "--charge",  dest="charge",  type="string", default='', help="Charge: either 'plus' or 'minus'");
    parser.add_option(       '--uncorrelate-ptscale-by-etaside', dest='uncorrelatePtScaleByEtaSide' , default=False, action='store_true', help='Uncorrelate momentum scales by eta sides')
    parser.add_option(       '--uncorrelate-ptscale-by-charge', dest='uncorrelatePtScaleByCharge' , default=False, action='store_true', help='Uncorrelate momentum scales by charge')
    (options, args) = parser.parse_args()

    # if not len(options.indir):
    #     print "Warning: you must specify a folder with signal root files with option --indir"
    #     quit()
    if not len(options.infile):
        print "Warning: you must specify an input file name with option --infile"
        quit()

    # manage output folder
    outdir = options.outdir
    if not outdir.endswith('/'): outdir += "/"
    if outdir != "./":
        if not os.path.exists(outdir):
            print "Creating folder", outdir
            os.system("mkdir -p " + outdir)

    #if options.flavour not in ["el", "mu"]:
    #    print "Warning: you must specify a lepton flavour with option -f el|mu"
    #    quit()
    #if options.charge not in ["plus", "minus"]:
    #    print "Warning: you must specify a charge with option -c plus|minus"
    #    quit()

    basenameFile = os.path.basename(options.infile)

    charge = ""
    flavour = ""
    if not any(ch in basenameFile for ch in ["plus", "minus"]):
        print "Warning: I could not understand charge plus|minus from file name. Abort"
        quit()
    else:
        charge = "plus" if "plus" in basenameFile else "minus"

    if not any(fl in basenameFile for fl in ["el_", "mu_"]):
        print "Warning: I could not understand flavour el|mu from file name. Abort"
        quit()
    else:
        flavour = "mu" if "mu_" in basenameFile else "el"

    isMu = True if flavour == "mu" else False

    #def addSmoothLepScaleSyst(infile,regexp,charge, isMu, outdir=None, 
    #                          uncorrelateByCharge=False,
    #                          uncorrelateByEtaSide=False,
    #                          getExtremePtBinsFromOldScale=True):

    addSmoothLepScaleSyst(options.infile, options.processes, charge, isMu,
                          outdir=outdir,
                          uncorrelateByCharge=options.uncorrelatePtScaleByCharge,
                          uncorrelateByEtaSide=options.uncorrelatePtScaleByEtaSide,
                          getExtremePtBinsFromOldScale=True)
