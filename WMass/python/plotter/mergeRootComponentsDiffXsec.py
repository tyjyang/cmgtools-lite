#!/usr/bin/env python                                                       
 
#from shutil import copyfile
import re, sys, os, os.path, ROOT, copy, math
from array import array
#import numpy as np
import root_numpy

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

from addSmoothPtScaleToShape import addSmoothLepScaleSyst
from addSmoothPtScaleToShape import addSmoothMuonScaleSyst # use this one for muons, with Rochester corrections
from addSmoothPtScaleToShape import addSmoothElectronScaleSyst
# experimental function for momentum scales
# def addSmoothLepScaleSyst(infile,regexp,charge, isMu, outdir=None, 
#                           uncorrelateByCharge=False,
#                           uncorrelateByEtaSide=False):

#     indir = outdir if outdir != None else options.inputdir
#     flav = 'mu' if isMu else 'el'

#     # get eta-pt binning for reco 
#     etaPtBinningVec = getDiffXsecBinning(indir+'/binningPtEta.txt', "reco")  # this get two vectors with eta and pt binning
#     recoBins = templateBinning(etaPtBinningVec[0],etaPtBinningVec[1])        # this create a class to manage the binnings
#     binning = [recoBins.Neta, recoBins.etaBins, recoBins.Npt, recoBins.ptBins]

#     tmp_infile = ROOT.TFile(infile, 'read')
#     outfile = ROOT.TFile(indir+'/SmoothScaleSyst_{flav}_{ch}.root'.format(flav=flav,ch=charge), 'recreate')

#     ## scale systematics as a function of eta
#     etabins_mu = array('f',[0.0, 2.1, 2.4])

#     scaleSyst_mu = ROOT.TH1F('scaleSyst_mu','',len(etabins_mu)-1,etabins_mu)
#     scaleSyst_mu.SetBinContent(1,0.003)
#     scaleSyst_mu.SetBinContent(2,0.010)
    
#     etabins_el = array('f',[0.0, 1.0, 1.5, 2.1, 2.4])
#     scaleSyst_el = ROOT.TH1F('scaleSyst_el','',len(etabins_el)-1,etabins_el)
#     scaleSyst_el.SetBinContent(1,0.003)
#     scaleSyst_el.SetBinContent(2,0.005)
#     scaleSyst_el.SetBinContent(3,0.008)
#     scaleSyst_el.SetBinContent(4,0.01)
        
#     if isMu:
#         scaleSyst = utilities.getExclusiveBinnedSyst(scaleSyst_mu)
#         etabins = etabins_mu
#     else:
#         scaleSyst = utilities.getExclusiveBinnedSyst(scaleSyst_el)
#         etabins = etabins_el

    
#     for k in tmp_infile.GetListOfKeys():
#         tmp_name = k.GetName()
#         ## don't reweight any histos that don't match the regexp
#         if not re.match(regexp, tmp_name): continue
#         ## don't reweight any histos that are already variations of something else
#         if 'Up' in tmp_name or 'Down' in tmp_name: continue
#         process = tmp_name.split('_')[1]

#         ## now should be left with only the ones we are interested in
#         print 'reweighting smooth lepscale syst for process', tmp_name
        
#         tmp_nominal = tmp_infile.Get(tmp_name)
#         tmp_nominal_2d = dressed2D(tmp_nominal,binning, tmp_name+'backrolled')
#         n_ptbins = tmp_nominal_2d.GetNbinsY()

#         etasides = [-1, 1] if uncorrelateByEtaSide else [1]

#         for side in etasides:
#             sideTag = "" if not uncorrelateByEtaSide else "etaside{s}".format(s="P" if side > 0 else "M")
#             for systBin in xrange(len(etabins)-1):
#                 etasyst = scaleSyst.GetXaxis().GetBinCenter(systBin+1)
#                 for shift_dir in ['Up','Down']:
#                     outname_2d = tmp_nominal_2d.GetName().replace('backrolled','')+'_smooth{lep}scale{idx}{ch}{st}{shiftdir}'.format(lep=flav,idx=systBin,shiftdir=shift_dir,ch=charge if uncorrelateByCharge else "",st=sideTag)
#                     tmp_scaledHisto = copy.deepcopy(tmp_nominal_2d.Clone(outname_2d))
#                     tmp_scaledHisto.Reset()

#                     ## loop over all eta bins of the 2d histogram 
#                     for ieta in range(1,tmp_nominal_2d.GetNbinsX()+1):
#                         eta = tmp_nominal_2d.GetXaxis().GetBinCenter(ieta)
#                         scale_syst = scaleSyst.GetBinContent(scaleSyst.GetXaxis().FindFixBin(etasyst))
#                         for ipt in range(1,tmp_nominal_2d.GetNbinsY()+1):                            
#                             if abs(eta)>etabins[systBin]:
#                                 if (uncorrelateByEtaSide and eta*float(side) < 0):
#                                     # we want to uncorrelate, but this eta is the other side
#                                     tmp_scaledHisto.SetBinContent(ieta, ipt, tmp_nominal_2d.GetBinContent(ieta,ipt))
#                                 else:
#                                     ## assume uniform distribution within a bin
#                                     if shift_dir == "Up":
#                                         if ipt == 1: continue
#                                         # assuming 1% scale variations and a bin with pt in [25, 26] GeV 
#                                         # (and previous bin from 24 to 25), only the events in the previous bin 
#                                         # for which (1+0.01)*pt' > 25 will move to the bin in [25, 26], so this pt' 
#                                         # is 25/1.01=24.7525 GeV, which corresponds to a fraction of events in 
#                                         # previous bin of 24.75%
#                                         # so, if ptl is the low-pt edge, and x is the scale variation (scaleSyst_mu/el)
#                                         # we can write DeltaN_low = |ptl - ptl/(1+x)|/binWidth_low * binContent_low =
#                                         # = x/(1+x) *ptl/binWidth_low * binContent_low
#                                         pt      = tmp_nominal_2d.GetYaxis().GetBinUpEdge(ipt)
#                                         pt_prev = tmp_nominal_2d.GetYaxis().GetBinLowEdge(ipt)
#                                         nominal_val      = tmp_nominal_2d.GetBinContent(ieta,ipt)
#                                         nominal_val_prev = tmp_nominal_2d.GetBinContent(ieta,ipt-1)
#                                         pt_width      = tmp_nominal_2d.GetYaxis().GetBinWidth(ipt)
#                                         pt_width_prev = tmp_nominal_2d.GetYaxis().GetBinWidth(ipt-1)
#                                         factor = scale_syst / (1.0 + scale_syst)
#                                         from_prev = factor * (pt_prev / pt_width_prev) * max(0,nominal_val_prev)
#                                         to_right  = factor * (pt      / pt_width     ) * max(0,nominal_val)
#                                         tmp_scaledHisto.SetBinContent(ieta,ipt,max(0,nominal_val + from_prev - to_right))
#                                     else:
#                                         if ipt == tmp_nominal_2d.GetNbinsY(): continue
#                                         # same logic as before, but now (1-0.01)*pt' < 25 --> pt'< 25.2525 
#                                         pt     = tmp_nominal_2d.GetYaxis().GetBinLowEdge(ipt)
#                                         pt_seq = tmp_nominal_2d.GetYaxis().GetBinUpEdge(ipt)
#                                         nominal_val     = tmp_nominal_2d.GetBinContent(ieta,ipt)
#                                         nominal_val_seq = tmp_nominal_2d.GetBinContent(ieta,ipt+1)
#                                         pt_width     = tmp_nominal_2d.GetYaxis().GetBinWidth(ipt)
#                                         pt_width_seq = tmp_nominal_2d.GetYaxis().GetBinWidth(ipt+1)
#                                         factor = scale_syst / (1.0 - scale_syst)
#                                         from_seq = factor * (pt_seq / pt_width_seq) * max(0,nominal_val_seq)
#                                         to_left  = factor * (pt     / pt_width    ) * max(0,nominal_val)
#                                         tmp_scaledHisto.SetBinContent(ieta,ipt,max(0,nominal_val + from_seq - to_left))
#                             else:
#                                 tmp_scaledHisto.SetBinContent(ieta, ipt, tmp_nominal_2d.GetBinContent(ieta,ipt))
#                         ## since in the first bin we cannot foresee 2-neighbors migrations, let's assume same syst of bin i+1
#                         if shift_dir == "Up":
#                             factor_bin2 = tmp_scaledHisto.GetBinContent(ieta,2)/tmp_nominal_2d.GetBinContent(ieta,2) if tmp_nominal_2d.GetBinContent(ieta,2) else 0
#                             tmp_scaledHisto.SetBinContent(ieta,1,factor_bin2 *tmp_nominal_2d.GetBinContent(ieta,1))
#                         else:
#                             binLastMinus1 = tmp_nominal_2d.GetNbinsY() - 1
#                             factor_binLastMinus1 = tmp_scaledHisto.GetBinContent(ieta,binLastMinus1)/tmp_nominal_2d.GetBinContent(ieta,binLastMinus1) if tmp_nominal_2d.GetBinContent(ieta,binLastMinus1) else 0
#                             tmp_scaledHisto.SetBinContent(ieta,tmp_nominal_2d.GetNbinsY(),factor_binLastMinus1 *tmp_nominal_2d.GetBinContent(ieta,tmp_nominal_2d.GetNbinsY()))

#                     ## re-roll the 2D to a 1D histo
#                     tmp_scaledHisto_1d = unroll2Dto1D(tmp_scaledHisto, newname=tmp_scaledHisto.GetName().replace('2DROLLED',''), cropNegativeBins=False)

#                     outfile.cd()
#                     tmp_scaledHisto_1d.Write()
#         # end loop on sides
#     outfile.Close()
#     print "done with the smooth lep scale variations"


def makeAlternateFromSymmetricRatio(alt, nomi, binning):
    # make ratio of alt and nomi, roll inot TH2 (pt-eta), symmetrize this ratio versus eta and then unroll again
    #
    # should I return the new alt or is it enough to modify the argument to really change it outside the function?
    if alt.GetNbinsX() != nomi.GetNbinsX():
        print "Error in makeAlternateFromSymmetricRatio: alt has %d bins, nomi has %d. Abort" % (alt.GetNbinsX(), nomi.GetNbinsX())
        quit()
    ratio1D = alt.Clone(alt.GetName()+"_ratio1D")
    ratio1D.Divide(nomi)
    name2D = alt.GetName() + "_ratio2D"
    ratio2D = dressed2D(ratio1D,binning,name2D,name2D)
    # now make the ratio symmetric versus central eta
    for ieta in range(1,ratio2D.GetNbinsX()+1):
        for ipt in range(1,ratio2D.GetNbinsY()+1):
            contentInThisEta = ratio2D.GetBinContent(ieta,ipt)  
            contentInSymmetricEta = ratio2D.GetBinContent(1+ratio2D.GetNbinsX()-ieta,ipt)  
            #ratio2D.SetBinContent(ieta,ipt, 0.5 * (contentInThisEta+contentInSymmetricEta))
            bin1D = (ipt-1) * ratio2D.GetNbinsX() + ieta
            ratio1D.SetBinContent(bin1D, 0.5 * (contentInThisEta+contentInSymmetricEta))
    # finally redefine the alt histogram, after unrolling again the ratio
    # ratio1D = unroll2Dto1D(ratio2D,newname=ratio1D.GetName(),cropNegativeBins=False) 
    #alt.Multiply(ratio1D,nomi)
    for ib in range(1, nomi.GetNbinsX()+1):
        alt.SetBinContent(ib, ratio1D.GetBinContent(ib) * nomi.GetBinContent(ib))

def makeAlternateFromSymmetricRatioV2(alt, nomi, binning):

    # make ratio of alt and nomi after symmetrizing a copy of num and den versus eta
    #
    # here we compute the ratio from the average of num and den in symmetric bins
    # i.e, we do ratio = (a+a')/(b+b'), where x and x' are bin contents at opposite eta
    # this is different than doing 0.5 * (a/b + a'/b'), i.e. the average of the ratios
    # if the nominal histogram is already very eta symmetric, then they are equivalent
    #print "Symmetrizing ratio of {a} and {n}".format(a=alt.GetName(), n=nomi.GetName())
    if alt.GetNbinsX() != nomi.GetNbinsX():
        print "Error in makeAlternateFromSymmetricRatioV2: alt has %d bins, nomi has %d. Abort" % (alt.GetNbinsX(), nomi.GetNbinsX())
        quit()
    altsym = alt.Clone(alt.GetName()+"_sym1D")
    nomisym = nomi.Clone(nomi.GetName()+"_sym1D")
    nEta = binning[0]
    nPt = binning[2]
    if nomi.GetNbinsX() != (nEta*nPt):
        print "Error in makeAlternateFromSymmetricRatioV2: nPt*nEta incompatible with number of bins of histogram. Abort"
        quit()

    # this assumes alt and nomi are unrolled so to have consecutive eta distributions at constant pt
    for ieta in range(1,nEta+1):
        for ipt in range(1,nPt+1):
            bin1D = (ipt-1) * nEta + ieta
            bin1Dsym = (ipt-1) * nEta + (1 + nEta - ieta) # i.e. for 48 eta bins take average of bin 1 and 48
            contentInThisEta = nomi.GetBinContent(bin1D)  
            contentInSymmetricEta = nomi.GetBinContent(bin1Dsym)  
            # skip multiplication by 0.5 to get average, it will cancel in the ratio
            nomisym.SetBinContent(bin1D, (contentInThisEta+contentInSymmetricEta))
            contentInThisEta = alt.GetBinContent(bin1D)  
            contentInSymmetricEta = alt.GetBinContent(bin1Dsym)  
            # skip multiplication by 0.5 to get average, it will cancel in the ratio
            altsym.SetBinContent(bin1D, (contentInThisEta+contentInSymmetricEta))

    ratio1D = altsym.Clone(altsym.GetName()+"_ratio1D")
    ratio1D.Divide(nomisym)
    # finally redefine the alt histogram, multiplying nomi by symmetrized ratio
    alt.Multiply(ratio1D,nomi)

def uncorrelateHistByEtaSide(hnomi, haltEtaPlus, haltEtaMinus, neta, npt):
    # to make it faster and avoid creating TH2, exploit the fact that hnomi and halt are unrolled along eta 
    #so eta distribution at each pt bin, also, neta is an even number
    netaHalf = neta/2
    #h2nomi = ROOT.TH2F("h2nomi_tmp","",neta,array('d',recoBins.etaBins),npt,array('d',recoBins.ptBins))
    #h2alt =  ROOT.TH2F("h2alt_tmp", "",neta,array('d',recoBins.etaBins),npt,array('d',recoBins.ptBins))
    for ieta in range(neta):
        for ipt in range(npt):
            globalbin = ieta + ipt * neta + 1
            if ieta < netaHalf: 
                haltEtaPlus.SetBinContent(globalbin, hnomi.GetBinContent(globalbin))
            else: 
                haltEtaMinus.SetBinContent(globalbin, hnomi.GetBinContent(globalbin))

def putTestEffSystHistosDiffXsec(infile,regexp,charge, outdir=None, isMu=True, suffix=""):

    # if not isMu:
    #     print "Electrons not implemented in putTestEffSystHistosDiffXsec(). Exit"
    #     quit()

    # this is mainly for tests, otherwise one should change the definition of weights when filling histograms
    # however, the efficiency systematics is almost flat in pt and varies only as a function of eta
    # this is not completely true for electrons, though

    indir = outdir if outdir != None else options.inputdir

    # get eta-pt binning for reco 
    etaPtBinningVec = getDiffXsecBinning(indir+'/binningPtEta.txt', "reco")  # this get two vectors with eta and pt binning
    recoBins = templateBinning(etaPtBinningVec[0],etaPtBinningVec[1])        # this create a class to manage the binnings
    binning = [recoBins.Neta, recoBins.etaBins, recoBins.Npt, recoBins.ptBins]
    etabins = recoBins.etaBins

    etaPtBinningVecGen = getDiffXsecBinning(indir+'/binningPtEta.txt', "gen")  # this get two vectors with eta and pt binning
    genBins = templateBinning(etaPtBinningVecGen[0],etaPtBinningVecGen[1])        # this create a class to manage the binnings


    tmp_infile = ROOT.TFile(infile, 'read')

    flavour = "mu" if isMu else "el"
    outfile = ROOT.TFile(indir+'/TestEffSyst_{fl}_{ch}{sfx}.root'.format(fl=flavour, ch=charge, sfx=suffix), 'recreate')

    ndone = 0
    for k in tmp_infile.GetListOfKeys():
        tmp_name = k.GetName()
        ## don't reweight any histos that don't match the regexp
        if not re.match(regexp, tmp_name): continue
        ## don't reweight any histos that are already variations of something else
        if 'Up' in tmp_name or 'Down' in tmp_name: continue

        ## now should be left with only the ones we are interested in
        # print 'reweighting testeffsyst nuisances for process', tmp_name
        
        tmp_nominal = tmp_infile.Get(tmp_name)
        tmp_nominal_2d = dressed2D(tmp_nominal,binning, tmp_name+'backrolled')

        effsystvals = {}
        if isMu:
            effsystvals[0.0] = 0.002
            effsystvals[1.0] = math.sqrt( math.pow(0.004,2) - math.pow(0.002,2))
            effsystvals[1.5] = math.sqrt( math.pow(0.014,2) - math.pow(0.004,2) - math.pow(0.002,2))
        else:            
            # this might get the wrong cmssw path if you typed cmsenv from a different release
            # L1EG_file = "{cmssw}/src/CMGTools/WMass/python/postprocessing/data/leptonSF/new2016_madeSummer2018/l1EG_eff.root".format(cmssw=os.environ["CMSSW_BASE"])
            # let's use relative path
            L1EG_file = "../postprocessing/data/leptonSF/new2016_madeSummer2018/l1EG_eff.root"
            L1EG_hist = "l1EG_eff"
            tf = ROOT.TFile(L1EG_file,'read')
            hL1 = tf.Get(L1EG_hist)
            if not hL1:
                print "Error in putTestEffSystHistosDiffXsec(): could not get histogram {h} in file {f}. Abort".format(h=L1EG_hist,f=L1EG_file)
                quit()
            # # FIX ME: here the prefiring syst is not added!
            #
            # numbers are deviced in such a way to obtain the commented value when summing in quadrature all terms before that
            effsystvals[0.0] =   0.006
            effsystvals[1.0] =   0.0053  # 0.008
            effsystvals[1.5] = 0.01    # 0.013
            effsystvals[2.1] =   0.0093  # 0.016
            # then we will no longer sum L1 prefiring term, that will be a new piece            

        for ik,key in enumerate(effsystvals):

            # for 2D xsec we assume that eta bins outside the corresponding gen bins are almost empty, so we will skip their reweigthing based on the gen eta
            if re.match("x_W.*_ieta_.*",tmp_name):
                ietagen,iptgen = get_ieta_ipt_from_process_name(tmp_name)
                #print "found process %s --> ietagen = %d" % (tmp_name, ietagen)
                # by uncommenting the following lines we skip histograms for bins which would not be affected
                # we might still want to have the histograms but just excluding the nuisance for some processes in the datacard
                # in this way we can at least define the template ratio when summing all histograms into an inclusive template in templateRolling.py
                # if ietagen >= 0:
                #     geneta = genBins.etaBins[ietagen]
                #     #print "geneta = %.3f" % geneta
                #     if geneta < key:   
                #         #print "skip reweighting of %s for TestEffSyst%d" % (tmp_name,ik) 
                #         continue                

            outname_2d = tmp_nominal_2d.GetName().replace('backrolled','')+'_{fl}TestEffSyst{ik}2DROLLED'.format(fl=flavour,ik=ik)

            tmp_scaledHisto_up = copy.deepcopy(tmp_nominal_2d.Clone(outname_2d+'Up'))
            tmp_scaledHisto_dn = copy.deepcopy(tmp_nominal_2d.Clone(outname_2d+'Down'))

            for ieta in range(1,tmp_scaledHisto_up.GetNbinsX()+1):
                for ipt in range(1,tmp_scaledHisto_up.GetNbinsY()+1):
                    etaval = tmp_scaledHisto_up.GetXaxis().GetBinCenter(ieta)
                    scaling = 0.0
                    if abs(etaval) >= key:
                        scaling += effsystvals[key]
                        # if isMu: 
                        #     scaling += effsystvals[key]                        
                        # else:
                        #     ptL1  = tmp_scaledHisto_up.GetYaxis().GetBinLowEdge(ipt)
                        #     etaL1 =  tmp_scaledHisto_up.GetXaxis().GetBinUpEdge(ieta) if (etaval < 0) else tmp_scaledHisto_up.GetXaxis().GetBinLowEdge(ieta)
                        #     L1RelativeUncertainty = 0
                        #     if (abs(etaL1) >= 1.479 and ptL1 >= 35.0):                                
                        #         ptval = tmp_scaledHisto_up.GetYaxis().GetBinCenter(ipt)
                        #         L1RelativeUncertainty = hL1.GetBinError(hL1.GetXaxis().FindFixBin(etaval),hL1.GetYaxis().FindFixBin(ptval))
                        #     scaling += math.sqrt(effsystvals[key]*effsystvals[key] + L1RelativeUncertainty*L1RelativeUncertainty)

                    tmp_bincontent = tmp_scaledHisto_up.GetBinContent(ieta, ipt)
                    tmp_bincontent_up = tmp_bincontent*(1.+scaling)
                    tmp_bincontent_dn = tmp_bincontent*(1.-scaling)
                    tmp_scaledHisto_up.SetBinContent(ieta, ipt, tmp_bincontent_up)
                    tmp_scaledHisto_dn.SetBinContent(ieta, ipt, tmp_bincontent_dn)

            ## re-roll the 2D to a 1D histo
            tmp_scaledHisto_up_1d = unroll2Dto1D(tmp_scaledHisto_up, newname=tmp_scaledHisto_up.GetName().replace('2DROLLED',''),cropNegativeBins=False)
            tmp_scaledHisto_dn_1d = unroll2Dto1D(tmp_scaledHisto_dn, newname=tmp_scaledHisto_dn.GetName().replace('2DROLLED',''),cropNegativeBins=False)

            outfile.cd()
            tmp_scaledHisto_up_1d.Write()
            tmp_scaledHisto_dn_1d.Write()

    outfile.Close()
    print 'done with the many reweightings for the testeffsyst'


def putEffStatHistosDiffXsec(infile,regexp,charge, outdir=None, isMu=True, suffix=""):

    # for differential cross section I don't use the same option for inputs, so I pass it from outside
    indir = outdir if outdir != None else options.inputdir

    # get eta-pt binning for reco 
    etaPtBinningVec = getDiffXsecBinning(indir+'/binningPtEta.txt', "reco")  # this get two vectors with eta and pt binning
    recoBins = templateBinning(etaPtBinningVec[0],etaPtBinningVec[1])        # this create a class to manage the binnings
    binning = [recoBins.Neta, recoBins.etaBins, recoBins.Npt, recoBins.ptBins]
    etabins = recoBins.etaBins

    basedir = '/afs/cern.ch/work/m/mdunser/public/cmssw/w-helicity-13TeV/CMSSW_8_0_25/src/CMGTools/WMass/python/postprocessing/data/'    
    if isMu:
        parfile_name = basedir+'/leptonSF/new2016_madeSummer2018/systEff_trgmu_{ch}_mu.root'.format(ch=charge)
    else:
        parfile_name = basedir+'/leptonSF/new2016_madeSummer2018/systEff_trgel.root'

    parfile = ROOT.TFile(parfile_name, 'read')
    

    tmp_infile = ROOT.TFile(infile, 'read')

    flavour = "mu" if isMu else "el"
    outfile = ROOT.TFile(indir+'/ErfParEffStat_{fl}_{ch}{sfx}.root'.format(fl=flavour, ch=charge, sfx=suffix), 'recreate')

    ndone = 0
    for k in tmp_infile.GetListOfKeys():
        tmp_name = k.GetName()
        ## don't reweight any histos that don't match the regexp
        if not re.match(regexp, tmp_name): continue
        ## don't reweight any histos that are already variations of something else
        if 'Up' in tmp_name or 'Down' in tmp_name: continue

        #if ndone: continue
        ndone += 1

        ## now should be left with only the ones we are interested in
        print 'reweighting erfpareffstat nuisances for process', tmp_name
        
        tmp_nominal = tmp_infile.Get(tmp_name)
        tmp_nominal_2d = dressed2D(tmp_nominal,binning, tmp_name+'backrolled')

        # the erfPar histogram has 50 bins for muons and electrons (for muons, the extremes are equal to the neighbours because eta goes up to 2.4 at most)
        # the binning in eta is always 0.1 wide
        # if the template has coarser granularity, that bin will have N variations, for the N 0.1-wide bins contained in it
        # this would require using the loop to reweight the histograms, because each on the N-th variation would only affect 1/N of the events in the wider bin
        # so, here we will reduce the scaling by 1/N as well.

        # get width of bins in multiples of 0.1
        # for electrons the region around the gap is odd
        # we don't need this value to be integer, otherwise keep in mind that 0.1 is not perfectly represented with float in python
        # which entails that 0.2/0.1 might not yield 2.0, but 1.999
        binwidths = []
        for ieta in range(1,tmp_nominal_2d.GetNbinsX()+1):
            binwidths.append(tmp_nominal_2d.GetXaxis().GetBinWidth(ieta)/0.1)  

        nEtaErfPar = 48 if isMu else 50

        ## loop over the three parameters
        for npar in range(3):
            parhist = parfile.Get('p'+str(npar))
            ## loop over all eta bins of the 2d histogram
            for ietaErf_tmp in range(1,nEtaErfPar+1):

                # FIXME, need more care for electron channel, template binning can be odd around the gap
                # identify template eta which contain that erfPar eta
                etabinOffset = 0
                # # the parhist histogram for muons is defined with 50 bins, where the two with |eta| in [2.4,2.5] are filled like the inner ones
                # # these bins for muons do not make sense: when looping on the bin number, we need to add an offset to skip the first bin
                # or we can just use both and have ErfParXEffStatY with Y from 2 to 49 instead than from 1 to 48 (1 and 50 won't be used)
                if isMu and parhist.GetNbinsX() == 50:
                    etabinOffset = 1
                ietaErf = ietaErf_tmp + etabinOffset

                tmp_parhist_val = parhist.GetXaxis().GetBinCenter(ietaErf)
                #ietaTemplate = tmp_nominal_2d.GetXaxis().FindFixBin( parhist.GetXaxis().GetBinCenter(ietaErf+etabinOffset) )
                ietaTemplate = tmp_nominal_2d.GetXaxis().FindFixBin( tmp_parhist_val )

                # if template has eta range narrower than parhist, continue
                if ietaTemplate == 0 or ietaTemplate == (1 + tmp_nominal_2d.GetNbinsX()):
                    continue

                if not isMu:
                    # if in the gap, skip this EffStat, it only creates troubles
                    if abs(tmp_parhist_val) > 1.4 and abs(tmp_parhist_val) < 1.5:
                        continue
                        #ietaTemplate = tmp_nominal_2d.GetXaxis().FindFixBin( 1.41 if tmp_parhist_val > 0 else -1.41)
                    elif abs(tmp_parhist_val) > 1.5 and abs(tmp_parhist_val) < 1.6:
                        continue
                        #ietaTemplate = tmp_nominal_2d.GetXaxis().FindFixBin( 1.57 if tmp_parhist_val > 0 else -1.57)

                # for electrons there is the gap:
                # eg, we can have 1.3, 1.4442, 1.5, 1.566, 1.7
                # and the bin centers of the erfPar are 1.35, 1.45, 1.55, ...
                # so with 1.35 we get the template bin from 1.3 to 1.4442
                # but with 1.45 we would fall in the empty bin, while we would want to reweight the template bin that would span between 1.4 to 1.5
                # in this case, the syst is ineffective. The same for 1.55, which select the other "side" of the empty bin around 1.5
                # one could device a fancy way to isolate the gap and assign a larger uncertainty around it but let's keep it simple

                #print "ErfPar{p}EffStat{ieta}".format(p=npar,ieta=ietaErf_tmp)
                outname_2d = tmp_nominal_2d.GetName().replace('backrolled','')+'_ErfPar{p}EffStat{ieta}{fl}{ch}2DROLLED'.format(p=npar,
                                                                                                                                ieta=ietaErf_tmp,
                                                                                                                                fl=flavour,
                                                                                                                                ch=charge)
            
                tmp_scaledHisto_up = copy.deepcopy(tmp_nominal_2d.Clone(outname_2d+'Up'))
                tmp_scaledHisto_dn = copy.deepcopy(tmp_nominal_2d.Clone(outname_2d+'Down'))
                
                ## loop over all pT bins in that bin of eta (which is ieta)
                for ipt in range(1,tmp_scaledHisto_up.GetNbinsY()+1):
                    tmp_bincontent = tmp_scaledHisto_up.GetBinContent(ietaTemplate, ipt)
                    ybincenter = tmp_scaledHisto_up.GetYaxis().GetBinCenter(ipt)
                    ## now get the content of the parameter variation!
                    phistybin = min(parhist.GetNbinsX(),parhist.GetYaxis().FindBin(ybincenter))
                    tmp_scale = parhist.GetBinContent(ietaErf, phistybin)
                    scaling = math.sqrt(2.)*tmp_scale
                    ## scale electrons by sqrt(2) due to the input file being charge inclusive
                    if flavour == 'el':
                        scaling *= math.sqrt(2)
                    # modify scaling if template bin has larger width than the ErfPar histogram
                    # assuming the two uncertainties have the same order of magnitude, we should scale down by sqrt(binWidth/0.1)
                    # however, we could just be conservative and not scale anything
                    #if binwidths[ietaTemplate-1] > 1.: scaling = scaling / binwidths[ietaTemplate-1]
                    #
                    # do scaling only for electrons for eta bins with width larger than or equal to 0.2, 
                    # as the effstat are probably overestimated (but let's try)
                    # should do for muons as well, but things are fine and we can stay a bit more conservative (note that this is only for DY and Tau)
                    if binwidths[ietaTemplate-1] > 1.8: scaling = scaling / math.sqrt( binwidths[ietaTemplate-1] )
                    ## scale up and down with what we got from the histo
                    tmp_bincontent_up = tmp_bincontent*(1.+scaling)
                    tmp_bincontent_dn = tmp_bincontent*(1.-scaling)
                    tmp_scaledHisto_up.SetBinContent(ietaTemplate, ipt, tmp_bincontent_up)
                    tmp_scaledHisto_dn.SetBinContent(ietaTemplate, ipt, tmp_bincontent_dn)

                ## re-roll the 2D to a 1D histo
                tmp_scaledHisto_up_1d = unroll2Dto1D(tmp_scaledHisto_up, newname=tmp_scaledHisto_up.GetName().replace('2DROLLED',''),cropNegativeBins=False)
                tmp_scaledHisto_dn_1d = unroll2Dto1D(tmp_scaledHisto_dn, newname=tmp_scaledHisto_dn.GetName().replace('2DROLLED',''),cropNegativeBins=False)

                outfile.cd()
                tmp_scaledHisto_up_1d.Write()
                tmp_scaledHisto_dn_1d.Write()
    outfile.Close()
    print 'done with the many reweightings for the erfpar effstat'


# for Z use this function to make EffStat variations instead of doing it with jobs. Signal is made with the loop, but could be done like this as well
# note that this function is based on an approximation, but since we apply the variation in an eta-pt dependent way, it is a very good approximation
# by doing it for signal as well, the loop and unpacking of TH3 to TH1 could be made much faster, maybe

# originally developed to make cards for diff xsec in lepton pt and |eta|
# merge signal shapes with Z and data+backgrounds
# first, merge shapes with Z and remove data from there (data is just the MC)
# then, merge signal, data and merged Z

# Z and data are assumed to be in folder given with --indir-bkg
# signal is assumed to be named W<flavour>_<charge>_shapes_signal.root

# test function
def putBinUncEffStatHistosDiffXsec(infile,regexp,charge, outdir=None, isMu=True, suffix="", useBinnedSF=True):

    if not isMu:
        print "putBinUncEffStatHistosDiffXsec() currently implemented for electrons only"
        quit()

    # for differential cross section I don't use the same option for inputs, so I pass it from outside
    indir = outdir if outdir != None else options.inputdir

    # get eta-pt binning for reco 
    etaPtBinningVec = getDiffXsecBinning(indir+'/binningPtEta.txt', "reco")  # this get two vectors with eta and pt binning
    recoBins = templateBinning(etaPtBinningVec[0],etaPtBinningVec[1])        # this create a class to manage the binnings
    binning = [recoBins.Neta, recoBins.etaBins, recoBins.Npt, recoBins.ptBins]
    etabins = recoBins.etaBins

    basedir = '/afs/cern.ch/work/m/mdunser/public/cmssw/w-helicity-13TeV/CMSSW_8_0_25/src/CMGTools/WMass/python/postprocessing/data/'    
    if isMu:
        parfile_name = basedir+'/leptonSF/new2016_madeSummer2018/smoothEfficiency_muons_{ch}_trigger.root'.format(ch=charge)
        staterrfile_name = "../postprocessing/data/leptonSF/new2016_madeSummer2018/TnPstuff/muon/triggerMuonEff{ch}_fromRooFitResult_onlyStatUnc.root".format(ch="Plus" if charge == "plus" else "Minus") # can't use CMSSW_BASE, because code runs from release for combinetf.pt, >= 10_3_X
    else:
        parfile_name = basedir+'/leptonSF/new2016_madeSummer2018/systEff_trgel.root'

    parfile = ROOT.TFile(parfile_name, 'read')
    staterrfile = ROOT.TFile(staterrfile_name, 'read')

    tmp_infile = ROOT.TFile(infile, 'read')

    flavour = "mu" if isMu else "el"
    outfile = ROOT.TFile(indir+'/BinUncEffStat_{fl}_{ch}{sfx}.root'.format(fl=flavour, ch=charge, sfx=suffix), 'recreate')

    ndone = 0
    for k in tmp_infile.GetListOfKeys():
        tmp_name = k.GetName()
        ## don't reweight any histos that don't match the regexp
        if not re.match(regexp, tmp_name): continue
        ## don't reweight any histos that are already variations of something else
        if 'Up' in tmp_name or 'Down' in tmp_name: continue

        #if ndone: continue
        ndone += 1

        ## now should be left with only the ones we are interested in
        print 'reweighting binunceffstat nuisances for process', tmp_name
        
        tmp_nominal = tmp_infile.Get(tmp_name)
        tmp_nominal_2d = dressed2D(tmp_nominal,binning, tmp_name+'backrolled')

        # the erfPar histogram has 50 bins for muons and electrons (for muons, the extremes are equal to the neighbours because eta goes up to 2.4 at most)
        # the binning in eta is always 0.1 wide
        # if the template has coarser granularity, that bin will have N variations, for the N 0.1-wide bins contained in it
        # this would require using the loop to reweight the histograms, because each on the N-th variation would only affect 1/N of the events in the wider bin
        # so, here we will reduce the scaling by 1/N as well.

        # get width of bins in multiples of 0.1
        # for electrons the region around the gap is odd
        # we don't need this value to be integer, otherwise keep in mind that 0.1 is not perfectly represented with float in python
        # which entails that 0.2/0.1 might not yield 2.0, but 1.999
        binwidths = []
        for ieta in range(1,tmp_nominal_2d.GetNbinsX()+1):
            binwidths.append(tmp_nominal_2d.GetXaxis().GetBinWidth(ieta)/0.1)  

        nEtaErfPar = 48 if isMu else 50

        ## loop over the three parameters
        for npar in range(1):
            parhist = parfile.Get('scaleFactor' + ("Original" if useBinnedSF else ""))
            staterrhist = staterrfile.Get('triggerSF_{ch}'.format(ch=charge))
            if not staterrhist:
                print "Error: histogram %s not found in %s" % ('triggerSF_{ch}'.format(ch=charge), staterrfile_name) 
                quit()
            ## loop over all eta bins of the 2d histogram
            for ietaErf_tmp in range(1,nEtaErfPar+1):

                # FIXME, need more care for electron channel, template binning can be odd around the gap
                # identify template eta which contain that erfPar eta

                etabinOffset = 0
                # # the parhist histogram for muons is defined with 50 bins, where the two with |eta| in [2.4,2.5] are filled like the inner ones
                # # these bins for muons do not make sense: when looping on the bin number, we need to add an offset to skip the first bin
                if isMu and parhist.GetNbinsX() == 50:
                     etabinOffset = 1
                ietaErf = ietaErf_tmp + etabinOffset

                tmp_parhist_val = parhist.GetXaxis().GetBinCenter(ietaErf)
                ietaTemplate = tmp_nominal_2d.GetXaxis().FindFixBin( tmp_parhist_val )

                # if template has eta range narrower than parhist, continue
                if ietaTemplate == 0 or ietaTemplate == (1 + tmp_nominal_2d.GetNbinsX()):
                    continue

                if not isMu:
                    # if in the gap, skip this EffStat, it only creates troubles 
                    if abs(tmp_parhist_val) > 1.4 and abs(tmp_parhist_val) < 1.5:
                        continue
                        #ietaTemplate = tmp_nominal_2d.GetXaxis().FindFixBin( 1.41 if tmp_parhist_val > 0 else -1.41) 
                    elif abs(tmp_parhist_val) > 1.5 and abs(tmp_parhist_val) < 1.6:
                        continue
                        #ietaTemplate = tmp_nominal_2d.GetXaxis().FindFixBin( 1.57 if tmp_parhist_val > 0 else -1.57) 
               

                # for electrons there is the gap:
                # eg, we can have 1.3, 1.4442, 1.5, 1.566, 1.7
                # and the bin centers of the erfPar are 1.35, 1.45, 1.55, ...
                # so with 1.35 we get the template bin from 1.3 to 1.4442
                # but with 1.45 we would fall in the empty bin, while we would want to reweight the template bin that would span between 1.4 to 1.5
                # in this case, the syst is ineffective. The same for 1.55, which select the other "side" of the empty bin around 1.5
                # one could device a fancy way to isolate the gap and assign a larger uncertainty around it but let's keep it simple

                #print "ErfPar{p}EffStat{ieta}".format(p=npar,ieta=ietaErf)
                outname_2d = tmp_nominal_2d.GetName().replace('backrolled','')+'_BinUncEffStat{ieta}{fl}{ch}2DROLLED'.format(p=npar,
                                                                                                                             ieta=ietaErf_tmp,
                                                                                                                             fl=flavour,
                                                                                                                             ch=charge)
            
                tmp_scaledHisto_up = copy.deepcopy(tmp_nominal_2d.Clone(outname_2d+'Up'))
                tmp_scaledHisto_dn = copy.deepcopy(tmp_nominal_2d.Clone(outname_2d+'Down'))
                
                ## loop over all pT bins in that bin of eta (which is ieta)
                for ipt in range(1,tmp_scaledHisto_up.GetNbinsY()+1):
                    tmp_bincontent = tmp_scaledHisto_up.GetBinContent(ietaTemplate, ipt)
                    ybincenter = tmp_scaledHisto_up.GetYaxis().GetBinCenter(ipt)
                    ## now get the content of the parameter variation!                    
                    phistybin = min(parhist.GetNbinsY(),parhist.GetYaxis().FindBin(ybincenter))
                    # if I need to reduce the variation as below, first get the variation (with proper sign), scale it, and then sum 1
                    binstatunc = staterrhist.GetBinError(staterrhist.GetXaxis().FindFixBin(tmp_parhist_val), 
                                                         staterrhist.GetYaxis().FindFixBin(ybincenter)) 
                    # inflate to account for other SF, here we are just using trigger
                    binstatunc *= (math.sqrt(2.) if isMu else 2.)  
                    tmp_scale_up = binstatunc / parhist.GetBinContent(ietaErf, phistybin)
                    tmp_scale_dn = -1.0 * tmp_scale_up
        
                    #tmp_scale_up =        parhist.GetBinError(ietaErf, phistybin) / parhist.GetBinContent(ietaErf, phistybin)
                    #tmp_scale_dn = -1.0 * parhist.GetBinError(ietaErf, phistybin) / parhist.GetBinContent(ietaErf, phistybin)
                    # modify scaling if template bin has larger width than the ErfPar histogram
                    # no longer scaling, trying to be conservative
                    if binwidths[ietaTemplate-1] > 1.8: 
                        tmp_scale_up = tmp_scale_up / math.sqrt( binwidths[ietaTemplate-1])
                        tmp_scale_dn = tmp_scale_dn / math.sqrt( binwidths[ietaTemplate-1])
                    tmp_scale_up = 1. +  tmp_scale_up
                    tmp_scale_dn = 1. +  tmp_scale_dn
                    ## scale up and down with what we got from the histo

                    tmp_bincontent_up = tmp_bincontent*tmp_scale_up
                    tmp_bincontent_dn = tmp_bincontent*tmp_scale_dn
                    tmp_scaledHisto_up.SetBinContent(ietaTemplate, ipt, tmp_bincontent_up)
                    tmp_scaledHisto_dn.SetBinContent(ietaTemplate, ipt, tmp_bincontent_dn)

                ## re-roll the 2D to a 1D histo
                tmp_scaledHisto_up_1d = unroll2Dto1D(tmp_scaledHisto_up, newname=tmp_scaledHisto_up.GetName().replace('2DROLLED',''),cropNegativeBins=False)
                tmp_scaledHisto_dn_1d = unroll2Dto1D(tmp_scaledHisto_dn, newname=tmp_scaledHisto_dn.GetName().replace('2DROLLED',''),cropNegativeBins=False)

                outfile.cd()
                tmp_scaledHisto_up_1d.Write()
                tmp_scaledHisto_dn_1d.Write()
    outfile.Close()
    print 'done with the many reweightings for the binunc effstat'


def putEtaPtBinUncEffStatHistosDiffXsec(infile,regexp,charge, outdir=None, isMu=True, suffix="", useBinnedSF=True):

    if not isMu:
        print "putEtaPtBinUncEffStatHistosDiffXsec() currently implemented for electrons only"
        quit()

    # for differential cross section I don't use the same option for inputs, so I pass it from outside
    indir = outdir if outdir != None else options.inputdir

    # get eta-pt binning for reco 
    etaPtBinningVec = getDiffXsecBinning(indir+'/binningPtEta.txt', "reco")  # this get two vectors with eta and pt binning
    recoBins = templateBinning(etaPtBinningVec[0],etaPtBinningVec[1])        # this create a class to manage the binnings
    binning = [recoBins.Neta, recoBins.etaBins, recoBins.Npt, recoBins.ptBins]
    etabins = recoBins.etaBins

    etaPtBinningVecGen = getDiffXsecBinning(indir+'/binningPtEta.txt', "gen")  # this get two vectors with eta and pt binning
    genBins = templateBinning(etaPtBinningVecGen[0],etaPtBinningVecGen[1])        # this create a class to manage the binnings

    basedir = '/afs/cern.ch/work/m/mdunser/public/cmssw/w-helicity-13TeV/CMSSW_8_0_25/src/CMGTools/WMass/python/postprocessing/data/'    
    if isMu:
        parfile_name = basedir+'/leptonSF/new2016_madeSummer2018/smoothEfficiency_muons_{ch}_trigger.root'.format(ch=charge)
        staterrfile_name = "../postprocessing/data/leptonSF/new2016_madeSummer2018/TnPstuff/muon/triggerMuonEff{ch}_fromRooFitResult_onlyStatUnc.root".format(ch="Plus" if charge == "plus" else "Minus") # can't use CMSSW_BASE, because code runs from release for combinetf.pt, >= 10_3_X
    else:
        parfile_name = basedir+'/leptonSF/new2016_madeSummer2018/systEff_trgel.root'

    parfile = ROOT.TFile(parfile_name, 'read')
    staterrfile = ROOT.TFile(staterrfile_name, 'read')

    tmp_infile = ROOT.TFile(infile, 'read')

    flavour = "mu" if isMu else "el"
    outfile = ROOT.TFile(indir+'/BinEtaPtUncorrUncEffStat_{fl}_{ch}{sfx}.root'.format(fl=flavour, ch=charge, sfx=suffix), 'recreate')

    ndone = 0
    for k in tmp_infile.GetListOfKeys():
        tmp_name = k.GetName()
        ## don't reweight any histos that don't match the regexp
        if not re.match(regexp, tmp_name): continue
        ## don't reweight any histos that are already variations of something else
        if 'Up' in tmp_name or 'Down' in tmp_name: continue

        #if ndone: continue
        ndone += 1

        ## now should be left with only the ones we are interested in
        print 'reweighting binetaptuncorrunceffstat nuisances for process', tmp_name
        
        tmp_nominal = tmp_infile.Get(tmp_name)
        tmp_nominal_2d = dressed2D(tmp_nominal,binning, tmp_name+'backrolled')

        # the erfPar histogram has 50 bins for muons and electrons (for muons, the extremes are equal to the neighbours because eta goes up to 2.4 at most)
        # the binning in eta is always 0.1 wide
        # if the template has coarser granularity, that bin will have N variations, for the N 0.1-wide bins contained in it
        # this would require using the loop to reweight the histograms, because each on the N-th variation would only affect 1/N of the events in the wider bin
        # so, here we will reduce the scaling by 1/N as well.

        # get width of bins in multiples of 0.1
        # for electrons the region around the gap is odd
        # we don't need this value to be integer, otherwise keep in mind that 0.1 is not perfectly represented with float in python
        # which entails that 0.2/0.1 might not yield 2.0, but 1.999
        binwidths = []
        for ieta in range(1,tmp_nominal_2d.GetNbinsX()+1):
            binwidths.append(tmp_nominal_2d.GetXaxis().GetBinWidth(ieta)/0.1)  

        nEtaErfPar = 48 if isMu else 50

        ## loop over the three parameters
        for npar in range(1):
            parhist = parfile.Get('scaleFactor' + ("Original" if useBinnedSF else ""))
            staterrhist = staterrfile.Get('triggerSF_{ch}'.format(ch=charge))
            if not staterrhist:
                print "Error: histogram %s not found in %s" % ('triggerSF_{ch}'.format(ch=charge), staterrfile_name) 
                quit()
            ## loop over all eta bins of the 2d histogram
            for ietaErf_tmp in range(1,nEtaErfPar+1):

                etabinOffset = 0
                # # the parhist histogram for muons is defined with 50 bins, where the two with |eta| in [2.4,2.5] are filled like the inner ones
                # # these bins for muons do not make sense: when looping on the bin number, we need to add an offset to skip the first bin
                if isMu and parhist.GetNbinsX() == 50:
                     etabinOffset = 1
                ietaErf = ietaErf_tmp + etabinOffset

                tmp_parhist_val = parhist.GetXaxis().GetBinCenter(ietaErf)
                ietaTemplate = tmp_nominal_2d.GetXaxis().FindFixBin( tmp_parhist_val )
                isSignalInAcc = False
                ietabin = 0.0
                iptbin = 0.0
                if re.match(".*_ieta_.*_ipt_.*",tmp_name):
                    ietabin,iptbin = get_ieta_ipt_from_process_name(tmp_name)
                    isSignalInAcc = True
                    if abs(tmp_parhist_val) < genBins.etaBins[ietabin] or abs(tmp_parhist_val) > genBins.etaBins[ietabin+1]:
                        continue
                
                # if template has eta range narrower than parhist, continue
                if ietaTemplate == 0 or ietaTemplate == (1 + tmp_nominal_2d.GetNbinsX()):
                    continue

                if not isMu:
                    # if in the gap, skip this EffStat, it only creates troubles 
                    if abs(tmp_parhist_val) > 1.4 and abs(tmp_parhist_val) < 1.5:
                        continue
                        #ietaTemplate = tmp_nominal_2d.GetXaxis().FindFixBin( 1.41 if tmp_parhist_val > 0 else -1.41) 
                    elif abs(tmp_parhist_val) > 1.5 and abs(tmp_parhist_val) < 1.6:
                        continue
                        #ietaTemplate = tmp_nominal_2d.GetXaxis().FindFixBin( 1.57 if tmp_parhist_val > 0 else -1.57) 
               

                # for electrons there is the gap:
                # eg, we can have 1.3, 1.4442, 1.5, 1.566, 1.7
                # and the bin centers of the erfPar are 1.35, 1.45, 1.55, ...
                # so with 1.35 we get the template bin from 1.3 to 1.4442
                # but with 1.45 we would fall in the empty bin, while we would want to reweight the template bin that would span between 1.4 to 1.5
                # in this case, the syst is ineffective. The same for 1.55, which select the other "side" of the empty bin around 1.5
                # one could device a fancy way to isolate the gap and assign a larger uncertainty around it but let's keep it simple

                #print "ErfPar{p}EffStat{ieta}".format(p=npar,ieta=ietaErf)
                
                ## loop over all pT bins in that bin of eta (which is ieta)
                for ipt in range(1, 1 + staterrhist.GetNbinsY()):

                    tmpptval = staterrhist.GetYaxis().GetBinCenter(ipt)
                    ptlow = staterrhist.GetYaxis().GetBinLowEdge(ipt)
                    pthigh = staterrhist.GetYaxis().GetBinUpEdge(ipt)
                    if pthigh <= recoBins.ptBins[0]: continue  # might happen only on electrons, neglect bins outside reco acceptance
                    if ptlow >= recoBins.ptBins[-1]: continue # if it is larger than last edge it is out of the template, so we skip it, but also if it is larger than the previous bin, because the algorithm would choose the next edge, which means we would still fall out of the template range
                    elif ptlow > recoBins.ptBins[-2]: continue # use > and not >= here, as the equality denotes a good range

                    # now we will reweight template's bins contained in ptlow, pthigh, choosing the edges just above if the binning does not coincide   
                    # e.g., if we have 25,27.5, the template bins are 26, 28
                    ptBinLow = -1
                    ptBinHigh = -1
                    for i in range(recoBins.Npt):
                        if ptBinLow < 0 and ptlow <= recoBins.ptBins[i]: 
                            ptBinLow = i
                            # to save time, check if upper edge is within acceptance, if it is out then set ptBinHigh to last bin and exit loop
                            if ptBinHigh < 0 and pthigh >= recoBins.ptBins[-2]:
                                ptBinHigh = recoBins.Npt
                                break
                        if ptBinHigh < 0 and pthigh <= recoBins.ptBins[i]: 
                            ptBinHigh = i
                            break # exit loop once the upper edge is matched
                    #print "pt range: SF [%.1f, %.1f], template [%.1f, %.1f] (%d bins)" % (ptlow,pthigh, 
                    #                                                                      recoBins.ptBins[ptBinLow], recoBins.ptBins[ptBinHigh],
                    #                                                                      1+ptBinHigh-ptBinLow)

                    # add 1 to use these indices as the bins of the histogram
                    ptBinLow += 1
                    ptBinHigh += 1

                    # at this point we know that the pt bin of the efficiency hist span the template bins from ptBinLow to ptBinHigh
                    # they could coincide, i.e. they can be a single bin
                    # we will only have to reweight those bins
                    # this algorithm is expected to lead to the reweighting of non-overlapping pt-segments of the template 
                    
                    # for signal make this pt variation only if we satisfy a given criterium
                    # it could be to do it only for bins which are close to the pt value of this bin
                    # I could choose only the central and neighboring bins, or all the bins if they have at least 1% or 0.5% of the integral
                    # former would be fine (and very simple to implement) if pt resolution were constant at all eta bins, but it is not the case
                    # if isSignalInAcc:
                    #     integralRatio = 0.005 # 0.5%
                    #     tmp_bincontent_integral = tmp_nominal_2d.Integral(ietaTemplate,ietaTemplate,  # single eta bin
                    #                                                       ptBinLow,ptBinHigh)         # potentially more pt bins
                    #     if tmp_bincontent_integral < (integralRatio * tmp_nominal_2d.Integral()):
                    #         continue

                    nametag= "BinEtaPtUncorrUncEffStat"
                    outname_2d = tmp_nominal_2d.GetName().replace('backrolled','')+'_{tag}{ieta}ipt{ipt}{fl}{ch}2DROLLED'.format(p=npar,
                                                                                                                                 tag=nametag,
                                                                                                                                 ieta=ietaErf_tmp,
                                                                                                                                 ipt=ipt,
                                                                                                                                 fl=flavour,
                                                                                                                                 ch=charge)
                                
                    tmp_scaledHisto_up = copy.deepcopy(tmp_nominal_2d.Clone(outname_2d+'Up'))
                    tmp_scaledHisto_dn = copy.deepcopy(tmp_nominal_2d.Clone(outname_2d+'Down'))
                    
                    # now reweight the corresponding pt bins (might be more than one)
                    for iptrange in range(ptBinLow,ptBinHigh+1):
                        tmp_bincontent = tmp_scaledHisto_up.GetBinContent(ietaTemplate, iptrange)
                        # following two lines would not be really needed, but we admit the possibility that the smoothed central values are being used
                        # in other words, the histograms with central values and uncertainties might not coincide (they would in the latest version of 
                        # the root file, i.e., the uncertainty stored in the SF histogram, either smoothed or original, should already be stat only)
                        ybincenter = tmp_scaledHisto_up.GetYaxis().GetBinCenter(iptrange)
                        phistybin = min(parhist.GetNbinsY(),parhist.GetYaxis().FindBin(ybincenter))
                        ## now get the content of the parameter variation!                    
                        # if I need to reduce the variation as below, first get the variation (with proper sign), scale it, and then sum 1
                        binstatunc = staterrhist.GetBinError(ietaErf,ipt) # ietaErf,ipt already define the bin of the histogram with binned uncertainties
                        # inflate to account for other SF, here we are just using trigger
                        binstatunc *= (math.sqrt(2.) if isMu else 2.)  
                        tmp_scale_up = binstatunc / parhist.GetBinContent(ietaErf, phistybin)
                        tmp_scale_dn = -1.0 * tmp_scale_up                    
                        tmp_scale_up = 1. +  tmp_scale_up
                        tmp_scale_dn = 1. +  tmp_scale_dn
                        ## scale up and down with what we got from the histo
                        tmp_bincontent_up = tmp_bincontent*tmp_scale_up
                        tmp_bincontent_dn = tmp_bincontent*tmp_scale_dn
                        tmp_scaledHisto_up.SetBinContent(ietaTemplate, iptrange, tmp_bincontent_up)
                        tmp_scaledHisto_dn.SetBinContent(ietaTemplate, iptrange, tmp_bincontent_dn)

                    ## re-roll the 2D to a 1D histo
                    tmp_scaledHisto_up_1d = unroll2Dto1D(tmp_scaledHisto_up, 
                                                         newname=tmp_scaledHisto_up.GetName().replace('2DROLLED',''),
                                                         cropNegativeBins=False)
                    tmp_scaledHisto_dn_1d = unroll2Dto1D(tmp_scaledHisto_dn, 
                                                         newname=tmp_scaledHisto_dn.GetName().replace('2DROLLED',''),
                                                         cropNegativeBins=False)

                    outfile.cd()
                    tmp_scaledHisto_up_1d.Write()
                    tmp_scaledHisto_dn_1d.Write()
    outfile.Close()
    print 'done with the many reweightings for the binunc effstat (eta-pt uncorrelated)'



## python mergeRootComponentsDiffXsec.py -f mu -c minus --indir-bkg  cards/diffXsec_mu_2018_11_24_group10_onlyBkg/part0/  --indir-sig cards/diffXsec_mu_2018_11_24_group10_onlyBkg/ -o cards/diffXsec_mu_2018_11_24_group10_onlyBkg/ -d

from optparse import OptionParser
parser = OptionParser(usage="%prog [options]")
parser.add_option("-d", "--dry-run", dest="dryrun",   action="store_true", default=False, help="Dry run: print commands but do not merge");
parser.add_option(      "--indir-bkg", dest="indirBkg", type="string", default="", help="Input folder with root files for Z, data and other backgrounds");
parser.add_option(      "--indir-sig", dest="indirSig", type="string", default="", help="Input folder with root files for signal");
parser.add_option("-o", "--outdir",    dest="outdir", type="string", default="./", help="Output folder (current one as default)");
parser.add_option("-n", "--name",      dest="name",   type="string", default="", help="Name for output root file (if not given, name is W<flavou>_<charge>_shapes.root ).");
parser.add_option("-s", "--suffix",    dest="suffix", type="string", default="", help="Suffix to add to output file before extension. Ineffective if using option -n < name>");
parser.add_option("-f", "--flavour",   dest="flavour", type="string", default='', help="Channel: either 'el' or 'mu'");
parser.add_option("-c", "--charge",    dest="charge", type="string", default='', help="Charge: either 'plus' or 'minus'");
parser.add_option(      "--etaBordersForFakesUncorr",    dest="etaBordersForFakesUncorr", type="string", default='0.5,1.0,1.5,2.0', help="Borders passed to function that creates the eta-uncorrelated normalization for fakes. Pass comma separated list (no 0 or outer edges, and only positive values");
parser.add_option(      "--no-qcdsyst-Z", dest="useQCDsystForZ", action="store_false", default=True, help="If False, do not store the muR,muF,muRmuF variations for Z (if they were present)");
#parser.add_option(      "--no-effstatsyst-Z", dest="useEffstatsystForZ", action="store_false", default=True, help="If False, do not create the Effstat variations for Z");
parser.add_option(       '--uncorrelate-fakes-by-charge', dest='uncorrelateFakesByCharge' , default=False, action='store_true', help='If True, nuisances for fakes are uncorrelated between charges (Eta, PtSlope, PtNorm)')
parser.add_option(       "--test-eff-syst", dest="testEffSyst",   action="store_true", default=False, help="Add some more nuisances to test efficiency systematics on signal and Z (possily other components as well)");
parser.add_option(       "--useBinUncEffStat", dest="useBinUncEffStat",   action="store_true", default=False, help="Add some more nuisances representing eff stat with SF shifted by their uncertainty (for tests)");
parser.add_option(       "--useBinEtaPtUncorrUncEffStat", dest="useBinEtaPtUncorrUncEffStat",   action="store_true", default=False, help="Add some more nuisances representing eff stat with SF shifted by their uncertainty (for tests)");
parser.add_option(       '--uncorrelate-QCDscales-by-charge', dest='uncorrelateQCDscalesByCharge' , default=False, action='store_true', help='Use charge-dependent QCD scales (on signal and tau)')
parser.add_option(       '--uncorrelate-ptscale-by-etaside', dest='uncorrelatePtScaleByEtaSide' , default=False, action='store_true', help='Uncorrelate momentum scales by eta sides')
parser.add_option(       '--uncorrelate-ptscale-by-charge', dest='uncorrelatePtScaleByCharge' , default=False, action='store_true', help='Uncorrelate momentum scales by charge')
parser.add_option(       '--use-smooth-ptscale', dest='useSmoothPtScales' , default=False, action='store_true', help='Use smooth pt scales instead of native ones')
parser.add_option(       '--uncorrelate-nuisances-by-charge', dest='uncorrelateNuisancesByCharge' , type="string", default='', help='Regular expression for nuisances to decorrelate between charges (note that some nuisances have dedicated options)')
parser.add_option(      "--symmetrize-syst-ratio",  dest="symSystRatio", type="string", default='', help="Regular expression matching systematics whose histograms should be obtained by making the ratio with nominal symmetric in eta");
parser.add_option(       '--distinguish-name-sig-as-bkg', dest='distinguishNameSigAsBkg' , default=False, action='store_true', help='Use different name (hardcoded for now) to identify name of signal processes that are treated as background. Mainly for first pt bins of electrons when combining with muons, so to keep treating them as background in the combination (requires another option to specify the bins)')
parser.add_option(       '--pt-range-bkg', dest='pt_range_bkg', action="append", type="float", nargs=2, default=[], help='Will treat signal templates with gen level pt in this range as background in the datacard. Takes two float as arguments (increasing order) and can specify multiple times. They should match bin edges and a bin is not considered as background if at least one edge is outside this range')
parser.add_option(       '--add-smooth-ptscale-extremePtFromStandard', dest='addSmoothPtScalesWithStandardOnExtremePtBins' , default=False, action='store_true', help='Add smooth pt scales along with native ones. The code behaves as if the old ones are used, but also add the smooth ones with the extreme pt bins being copied from the standard ones. It requires not using --use-smooth-ptscale')
parser.add_option(       '--use-native-MCatNLO-xsec'  , dest='useNativeMCatNLOxsecW', default=False, action='store_true', help='Use native MC@NLO xsec for W and tau (actually, just scale W, signal should already have it)')
parser.add_option(       '--use-analytic-smooth-pt-scales'  , dest='useAnalyticSmoothPtScales', default=False, action='store_true', help='Use smooth pt scales made with analytic method (using special histograms with larger pt acceptance for the first and last pt bins). For muons, this uses the actual uncertainties for the Rochester corrections')
#parser.add_option("--pt-bins-name-change", dest="ptBinsNameChange", type="string", default='0,1', help="Comma separated list of pt bins for which the process name should be changed");
(options, args) = parser.parse_args()
    
# manage output folder
outdir = options.outdir
if not outdir.endswith('/'): outdir += "/"
if outdir != "./":
    if not os.path.exists(outdir):
        print "Creating folder", outdir
        os.system("mkdir -p " + outdir)

if options.flavour not in ["el", "mu"]:
    print "Warning: you must specify a lepton flavour with option -f el|mu"
    quit()
if options.charge not in ["plus", "minus"]:
    print "Warning: you must specify a charge with option -c plus|minus"
    quit()

charge = options.charge
flavour = options.flavour

if not len(options.indirBkg):
    print "Warning: you must specify a folder with data and background root files with option --indir-bkg"
    quit()
if not len(options.indirSig):
    print "Warning: you must specify a folder with signal root files with option --indir-sig"
    quit()

if options.useBinEtaPtUncorrUncEffStat and options.useBinUncEffStat:
    print "Warning: useBinUncEffStat and useBinEtaPtUncorrUncEffStat are incompatible. Use only one of them"
    quit()

# define ultimate output file
shapename = ""
if options.name == "":
    shapename = "{od}W{fl}_{ch}_shapes.root".format(od=outdir, fl=flavour, ch=charge)
    if len(options.suffix): shapename = shapename.replace(".root","_{sf}.root".format(sf=options.suffix))
else:
    shapename = outdir + options.name

print ""
print "-"*20
print ""
#following array is used to call function dressed2D() 
etaPtBinningVec = getDiffXsecBinning(options.indirSig+'/binningPtEta.txt', "reco")  # this get two vectors with eta and pt binning
recoBins = templateBinning(etaPtBinningVec[0],etaPtBinningVec[1])        # this create a class to manage the binnings
recoBinning = [recoBins.Neta, recoBins.etaBins, recoBins.Npt, recoBins.ptBins]

etaPtBinningVecGen = getDiffXsecBinning(options.indirSig+'/binningPtEta.txt', "gen")  
genBins = templateBinning(etaPtBinningVecGen[0],etaPtBinningVecGen[1])      
# consider some signal bins as background                                        
ptBinIsBackground = []  # will store a bool to assess whether the given ipt index is considered as background   
hasPtRangeBkg = False
for bin in range(genBins.Npt):
    ptBinIsBackground.append(False)
ptRangesBkg = options.pt_range_bkg
if len(ptRangesBkg):
    hasPtRangeBkg = True
    print "Signal bins with gen pt in the following ranges will be considered as background processes"
    print options.pt_range_bkg
    for index in range(genBins.Npt):
        for pair in ptRangesBkg:
        #print pair                          
            if genBins.ptBins[index] >= pair[0] and genBins.ptBins[index+1] <= pair[1]:
                ptBinIsBackground[index] = True
else:
    hasPtRangeBkg = False

## prepare the relevant files. First merge Z with correct charge
zMatch = "^Z_{fl}_{ch}.*".format(fl=flavour,ch=charge)
zfiles = [os.path.join(dp, f) for dp, dn, fn in os.walk(options.indirBkg) for f in fn if (f.endswith('.input.root') and re.match(zMatch,f))]

#for f in zfiles:
#    print f

# merge Z files (will have to remove x_data_obs later)
tmpZfile="{od}Z_{fl}_{ch}_mergeTMP.root".format(od=outdir, fl=flavour, ch=charge)
cmdMerge = "hadd -f -k -O {tmp} {zfs}".format(tmp=tmpZfile, zfs=" ".join([f for f in zfiles]))
#print cmdMerge
Zfile="{od}Z_{fl}_{ch}_merge.root".format(od=outdir, fl=flavour, ch=charge)

print "Merging Z"
if not options.dryrun: os.system(cmdMerge)

print "Remove x_data_obs from Z, and replace 'x_Z_dy_' with 'x_Z_'"
print "Also changing Dn to Down"
if not options.useQCDsystForZ:
    print "Will reject the QCD scales variations on Z (muR, muF, muRmuF)"
nZcopied = 0
# open new Z file to remove data from input file
#----------------------------------
newname = ""
zpdf = []
znominal = None
nMirroredPDF = 0
if not options.dryrun:
    tf = ROOT.TFile.Open(tmpZfile,"READ")
    if not tf or not tf.IsOpen():
        raise RuntimeError('Unable to open file {fn}'.format(fn=tmpZfile))

    # open output file
    of = ROOT.TFile(Zfile,'recreate')
    if not of or not of.IsOpen():
        raise RuntimeError('Unable to open file {fn}'.format(fn=Zfile))

    for ikey,e in enumerate(tf.GetListOfKeys()):
        name = e.GetName()
        obj  = e.ReadObj()
        if not obj:
            raise RuntimeError('Unable to read object {n}'.format(n=name))
        if "data_obs" in name: continue
        if not options.useQCDsystForZ:
            if any(x in name for x in ["muR", "muF", "muRmuF"]): continue
        if "x_Z_dy_" in name: newname = name.replace("x_Z_dy_","x_Z_")
        else: newname = name
        if name == "x_Z": znominal = obj.Clone(name)
        if "_pdf" in name:
            zpdf.append(obj.Clone(newname))
            continue
        if newname.endswith("Dn"): newname = newname.replace("Dn","Down")
        newobj = obj.Clone(newname)
        newobj.Write(newname)
        nZcopied += 1

    for h in zpdf:
        nMirroredPDF += 1
        (alternate,mirror) = mirrorShape(znominal,h,h.GetName(),use2xNomiIfAltIsZero=True)
        for alt in [alternate,mirror]:
            alt.Write()

    of.Close()
    tf.Close()
#----------------------------------

print "Copied {n} histograms in {zf} (x_data_obs removed)".format(n=str(nZcopied),zf=Zfile)
print "Created {n} mirrored histograms for PDFs in {zf} ".format(n=str(nMirroredPDF),zf=Zfile)
print "Removing temporary file {tmp}".format(tmp=tmpZfile)
if not options.dryrun: os.system("rm {tmp}".format(tmp=tmpZfile))    
print "-"*20
print ""

fileZeffStat = ""
if options.useBinEtaPtUncorrUncEffStat or options.useBinUncEffStat:
    pass
else:
    fileZeffStat = "{od}ErfParEffStat_{fl}_{ch}.root".format(od=outdir, fl=flavour, ch=options.charge)  
    # this name is used inside putEffStatHistosDiffXsec (do not change it outside here)
    print "Now adding ErfParXEffStatYY systematics to x_Z process"
    print "Will create file --> {of}".format(of=fileZeffStat)
    if not options.dryrun: putEffStatHistosDiffXsec(Zfile, 'x_Z', charge, outdir, isMu=True if flavour=="mu" else False)
    print ""

    print ""
    print "-"*20
    print ""

fileZbinUncEffStat = ""  
if options.useBinUncEffStat:
    fileZbinUncEffStat = "{od}BinUncEffStat_{fl}_{ch}.root".format(od=outdir, fl=flavour, ch=options.charge)  
    # this name is used inside putBinUncEffStatHistosDiffXsec (do not change it outside here)
    print "Now adding BinUncEffStatEffStatYY systematics to x_Z process"
    print "Will create file --> {of}".format(of=fileZbinUncEffStat)
    if not options.dryrun: putBinUncEffStatHistosDiffXsec(Zfile, 'x_Z', charge, outdir, isMu=True if flavour=="mu" else False)
    print ""
    print ""
    print "-"*20
    print ""

fileZbinEtaPtUncorrUncEffStat = ""  
if options.useBinEtaPtUncorrUncEffStat:
    fileZbinEtaPtUncorrUncEffStat = "{od}BinEtaPtUncorrUncEffStat_{fl}_{ch}.root".format(od=outdir, fl=flavour, ch=options.charge)  
    # this name is used inside putBinUncEffStatHistosDiffXsec (do not change it outside here)
    print "Now adding BinEtaPtUncorrUncEffStatEffStatYYiptZZ systematics to x_Z process"
    print "Will create file --> {of}".format(of=fileZbinEtaPtUncorrUncEffStat)
    if not options.dryrun: putEtaPtBinUncEffStatHistosDiffXsec(Zfile, 'x_Z', charge, outdir, isMu=True if flavour=="mu" else False)
    print ""
    print ""
    print "-"*20
    print ""


## prepare the relevant files. First merge Tau with correct charge
tauMatch = "^TauDecaysW_{fl}_{ch}.*".format(fl=flavour,ch=charge)
taufiles = [os.path.join(dp, f) for dp, dn, fn in os.walk(options.indirBkg) for f in fn if (f.endswith('.input.root') and re.match(tauMatch,f))]

#for f in zfiles:
#    print f

# merge Z files (will have to remove x_data_obs later)
tmpTaufile="{od}TauDecaysW_{fl}_{ch}_mergeTMP.root".format(od=outdir, fl=flavour, ch=charge)
cmdMerge = "hadd -f -k -O {tmp} {zfs}".format(tmp=tmpTaufile, zfs=" ".join([f for f in taufiles]))
#print cmdMerge
Taufile="{od}TauDecaysW_{fl}_{ch}_merge.root".format(od=outdir, fl=flavour, ch=charge)

print "Merging Tau"
if not options.dryrun: os.system(cmdMerge)

print "Remove x_data_obs from Tau, and replace 'x_TauDecaysW_wtau_' with 'x_TauDecaysW_'"
print "Also changing Dn to Down"
#print "Also fixing (mu|ele)scale postfit to CMS_W(mu|e)_(mu|ele)scale"
#if not options.useQCDsystForZ:
#    print "Will reject the QCD scales variations on Z (muR, muF, muRmuF)"
nTaucopied = 0
# open new TAU file to remove data from input file
#----------------------------------
newname = ""
taupdf = []
taunominal = None
nMirroredPDF = 0
if not options.dryrun:
    tf = ROOT.TFile.Open(tmpTaufile,"READ")
    if not tf or not tf.IsOpen():
        raise RuntimeError('Unable to open file {fn}'.format(fn=tmpTaufile))

    # open output file
    of = ROOT.TFile(Taufile,'recreate')
    if not of or not of.IsOpen():
        raise RuntimeError('Unable to open file {fn}'.format(fn=Taufile))

    for ikey,e in enumerate(tf.GetListOfKeys()):
        name = e.GetName()
        obj  = e.ReadObj()
        if not obj:
            raise RuntimeError('Unable to read object {n}'.format(n=name))
        if "data_obs" in name: continue
        #if not options.useQCDsystForZ:
        #    if any(x in name for x in ["muR", "muF", "muRmuF"]): continue
        newname = name
        if "x_TauDecaysW_wtau_" in name: newname = newname.replace("x_TauDecaysW_wtau_","x_TauDecaysW_")        
        # patch for bad name due to syst missing in systEnv (but not mca) when running jobs
        # it happened in the muon channel only
        # should not be needed anymore now, but keep the lines for future reference
        # if not options.useSmoothPtScales:
        #     if re.match(".*(mu|ele)scale.*",newname) and not re.match(".*CMS_W(mu|e)_(mu|ele)scale.*",newname):
        #         if re.match(".*_(Up|Dn|Down)",newname): 
        #             newname = newname.replace("_Up","Up").replace("_Dn","Dn").replace("_Down","Down") 
        #         if flavour == "el": 
        #             newname = newname.replace("elescale", "CMS_We_elescale")
        #         else: 
        #             newname = newname.replace("muscale", "CMS_Wmu_muscale")
        if name == "x_TauDecaysW": taunominal = obj.Clone(name)
        if "_pdf" in name:
            taupdf.append(obj.Clone(newname))
            continue
        if newname.endswith("Dn"): newname = newname.replace("Dn","Down")
        newobj = obj.Clone(newname)
        if options.useNativeMCatNLOxsecW:
            # 60400./(3*20508.9)  first one is MC@NLO, second one is fewz3.1 
            # also add factor 1.014 for normalization to gen-xsec after Wpt reweighting (this is needed
            # because I copied the old files made with old norm factor)
            newobj.Scale(0.981688 * 1.014) 
        newobj.Write(newname)
        nTaucopied += 1

    for h in taupdf:
        nMirroredPDF += 1
        (alternate,mirror) = mirrorShape(taunominal,h,h.GetName(),use2xNomiIfAltIsZero=True)
        for alt in [alternate,mirror]:
            alt.Write()

    of.Close()
    tf.Close()
#----------------------------------

print "Copied {n} histograms in {zf} (x_data_obs removed)".format(n=str(nTaucopied),zf=Taufile)
print "Created {n} mirrored histograms for PDFs in {zf} ".format(n=str(nMirroredPDF),zf=Taufile)
print "Removing temporary file {tmp}".format(tmp=tmpTaufile)
if not options.dryrun: os.system("rm {tmp}".format(tmp=tmpTaufile))    
print "-"*20
print ""

fileTaueffStat = ""
if options.useBinEtaPtUncorrUncEffStat or options.useBinUncEffStat:
    pass
else:
    fileTaueffStat = "{od}ErfParEffStat_{fl}_{ch}_Tau.root".format(od=outdir, fl=flavour, ch=options.charge)  
    # this name is used inside putEffStatHistosDiffXsec (do not change it outside here)
    print "Now adding ErfParXEffStatYY systematics to x_TauDecaysW process"
    print "Will create file --> {of}".format(of=fileTaueffStat)
    if not options.dryrun: putEffStatHistosDiffXsec(Taufile, 'x_TauDecaysW', charge, outdir, isMu=True if flavour=="mu" else False, suffix="_Tau")
    print ""
    print ""
    print "-"*20
    print ""

fileTaubinUncEffStat = ""
if options.useBinUncEffStat:
    fileTaubinUncEffStat = "{od}BinUncEffStat_{fl}_{ch}_Tau.root".format(od=outdir, fl=flavour, ch=options.charge)  
    # this name is used inside putBinUncEffStatHistosDiffXsec (do not change it outside here)
    print "Now adding BinUncEffStatYY systematics to x_TauDecaysW process"
    print "Will create file --> {of}".format(of=fileTaubinUncEffStat)
    if not options.dryrun: putBinUncEffStatHistosDiffXsec(Taufile, 'x_TauDecaysW', charge, outdir, isMu=True if flavour=="mu" else False, suffix="_Tau")
    print ""
    print ""
    print "-"*20
    print ""

fileTaubinEtaPtUncorrUncEffStat = ""
if options.useBinEtaPtUncorrUncEffStat:
    fileTaubinEtaPtUncorrUncEffStat = "{od}BinEtaPtUncorrUncEffStat_{fl}_{ch}_Tau.root".format(od=outdir, fl=flavour, ch=options.charge)  
    # this name is used inside putBinUncEffStatHistosDiffXsec (do not change it outside here)
    print "Now adding BinEtaPtUncorrUncEffStatYYiptZZ systematics to x_TauDecaysW process"
    print "Will create file --> {of}".format(of=fileTaubinEtaPtUncorrUncEffStat)
    if not options.dryrun: putEtaPtBinUncEffStatHistosDiffXsec(Taufile, 'x_TauDecaysW', charge, outdir, isMu=True if flavour=="mu" else False, suffix="_Tau")
    print ""
    print ""
    print "-"*20
    print ""


dataAndBkgFile = "{obkg}bkg_and_data_{fl}_{ch}.input.root".format(obkg=options.indirBkg, fl=flavour, ch=charge)
dataAndBkgFileTmp = dataAndBkgFile.replace(".input.root","TMP.input.root")
print "Creating temporary file {bkg} to remove 'x_data' histogram".format(bkg=dataAndBkgFileTmp)
print "Also changing Dn to Down"
print "-"*30
print "WARNING: will also reject histograms whose name starts with x_Z, which are not supposed to stay in this file"
print "-"*30
# now remove x_data from bkg_and_data_*.input.root (we only need x_data_obs, this x_data should be removed)
#----------------------------------
if not options.dryrun:
    tf = ROOT.TFile.Open(dataAndBkgFile,"READ")
    if not tf or not tf.IsOpen():
        raise RuntimeError('Unable to open file {fn}'.format(fn=dataAndBkgFile))

    # open output file
    of = ROOT.TFile(dataAndBkgFileTmp,'recreate')
    if not of or not of.IsOpen():
        raise RuntimeError('Unable to open file {fn}'.format(fn=dataAndBkgFileTmp))
    
    nKeys = tf.GetNkeys()
    nCopiedKeys = 0
    for ikey,e in enumerate(tf.GetListOfKeys()):
        name = e.GetName()
        obj  = e.ReadObj()
        if not obj:
            raise RuntimeError('Unable to read object {n}'.format(n=name))
        if name == "x_data": continue
        if "x_Z_" in name: continue  # patch, before February 2019 the muscale and lepeff systs on Z appear in the background file by mistake
        newname = name
        if newname.endswith("Dn"): 
            newname = newname.replace("Dn","Down")        
        newobj = obj.Clone(newname)
        newobj.Write(newname)        
        nCopiedKeys += 1

    print "Copied {n}/{tot} from {bkg}".format(n=str(nCopiedKeys),tot=str(nKeys),bkg=dataAndBkgFile)
    of.Close()
    tf.Close()
#----------------------------------
print "-"*20
print ""

fileFakesEtaUncorr = "{od}FakesEtaUncorrelated_{fl}_{ch}.root".format(od=outdir, fl=flavour, ch=options.charge) 
fileFakesEtaChargeUncorr = "{od}FakesEtaChargeUncorrelated_{fl}_{ch}.root".format(od=outdir, fl=flavour, ch=options.charge) 
fileFakesPtSlopeUncorr  = "{od}FakesPtSlopeUncorrelated_{fl}_{ch}.root".format( od=outdir, fl=flavour, ch=options.charge) 
fileFakesPtNormUncorr  = "{od}FakesPtNormUncorrelated_{fl}_{ch}.root".format( od=outdir, fl=flavour, ch=options.charge) 
# these names are used inside putUncorrelatedFakes (do not change them outside here)
print "Now adding Fakes*Uncorrelated systematics to x_data_fakes process"
print "Will create file --> {of}".format(of=fileFakesEtaUncorr)
print "Will create file --> {of}".format(of=fileFakesEtaChargeUncorr)
print "Will create file --> {of}".format(of=fileFakesPtSlopeUncorr)
print "Will create file --> {of}".format(of=fileFakesPtNormUncorr)
etaBordersForFakes = [float(x) for x in options.etaBordersForFakesUncorr.split(',')]
if not options.dryrun: 
    putUncorrelatedFakes(dataAndBkgFileTmp, 'x_data_fakes', charge, outdir, isMu=True if flavour=="mu" else False, etaBordersTmp=etaBordersForFakes, 
                         doType='eta', uncorrelateCharges=options.uncorrelateFakesByCharge, isHelicity=False)
    putUncorrelatedFakes(dataAndBkgFileTmp, 'x_data_fakes', charge, outdir, isMu=True if flavour=="mu" else False, etaBordersTmp=etaBordersForFakes, 
                        doType='etacharge', uncorrelateCharges=True, isHelicity=False) # this is always true
    putUncorrelatedFakes(dataAndBkgFileTmp, 'x_data_fakes', charge, outdir, isMu=True if flavour=="mu" else False, etaBordersTmp=etaBordersForFakes, 
                         doType='ptslope', uncorrelateCharges=options.uncorrelateFakesByCharge, isHelicity=False)
    putUncorrelatedFakes(dataAndBkgFileTmp, 'x_data_fakes', charge, outdir, isMu=True if flavour=="mu" else False, etaBordersTmp=etaBordersForFakes, 
                         doType='ptnorm', uncorrelateCharges=options.uncorrelateFakesByCharge, isHelicity=False)

sigfile = "{osig}W{fl}_{ch}_shapes_signal.root".format(osig=options.indirSig, fl=flavour, ch=charge)
fileSignalbinEtaPtUncorrUncEffStat = ""  
if options.useBinEtaPtUncorrUncEffStat:
    fileSignalbinEtaPtUncorrUncEffStat = "{od}BinEtaPtUncorrUncEffStat_{fl}_{ch}_signal.root".format(od=outdir, fl=flavour, ch=options.charge)  
    # this name is used inside putBinUncEffStatHistosDiffXsec (do not change it outside here)
    print "Now adding BinEtaPtUncorrUncEffStatEffStatYYiptZZ systematics to signal process"
    print "Will create file --> {of}".format(of=fileSignalbinEtaPtUncorrUncEffStat)
    if not options.dryrun: putEtaPtBinUncEffStatHistosDiffXsec(sigfile, 'x_W{ch}_lep'.format(ch=charge), 
                                                               charge, outdir, isMu=True if flavour=="mu" else False, suffix="_signal")
    print ""

    print ""
    print "-"*20
    print ""


print "Now merging signal + Z + data + other backgrounds + Fakes*Uncorrelated + ZEffStat"

cmdFinalMerge="hadd -f -k -O {of} {sig} {zf} {tauf} {bkg} ".format(of=shapename,  
                                                                   sig=sigfile, 
                                                                   zf=Zfile, 
                                                                   tauf=Taufile, 
                                                                   bkg=dataAndBkgFileTmp)

cmdFinalMerge += " {fakesEta} {fakesEtaCharge} {fakesPtSlope} {fakesPtNorm} {zEffStat} {tauEffStat} ".format(fakesEta=fileFakesEtaUncorr,
                                                                                                             fakesEtaCharge=fileFakesEtaChargeUncorr,
                                                                                                             fakesPtSlope=fileFakesPtSlopeUncorr,
                                                                                                             fakesPtNorm=fileFakesPtNormUncorr,
                                                                                                             zEffStat=fileZeffStat,
                                                                                                             tauEffStat=fileTaueffStat)

if options.useBinUncEffStat:
    cmdFinalMerge += "{zBinUncEffStat} {tauBinUncEffStat} ".format(zBinUncEffStat=fileZbinUncEffStat,tauBinUncEffStat=fileTaubinUncEffStat)
elif options.useBinEtaPtUncorrUncEffStat:
    cmdFinalMerge += "{zBinUncEffStat} {tauBinUncEffStat} {sigBinUncEffStat}".format(zBinUncEffStat=fileZbinEtaPtUncorrUncEffStat,
                                                                                     tauBinUncEffStat=fileTaubinEtaPtUncorrUncEffStat,
                                                                                     sigBinUncEffStat=fileSignalbinEtaPtUncorrUncEffStat)



print "Final merging ..."
print cmdFinalMerge
if not options.dryrun: os.system(cmdFinalMerge)

print "-"*20
print ""
print "Removing temporary file {bkg}".format(bkg=dataAndBkgFileTmp)
if not options.dryrun: os.system("rm {bkg}".format(bkg=dataAndBkgFileTmp))

print "-"*20
print ""
print ""
print "Wrote root file in %s" % shapename
print ""

if options.testEffSyst:
    # test
    fileTestEffSyst = "{od}TestEffSyst_{fl}_{ch}.root".format(od=outdir, fl=flavour, ch=options.charge)  
    # this name is used inside putTestEffSystHistosDiffXsec (do not change it outside here)
    print "Now adding TestEffSystYY systematics to x_Z, x_TauDecaysW, and x_W.* process"
    print "Will create file --> {of}".format(of=fileTestEffSyst)
    shapenameNoTestEffSyst = shapename.replace("_shapes","_shapes_noTestEffSyst")
    print "Copying {od} into {odnew}".format(od=shapename,odnew=shapenameNoTestEffSyst)
    if not options.dryrun: 
        os.system("cp {of} {ofnew}".format(of=shapename,ofnew=shapenameNoTestEffSyst))
        putTestEffSystHistosDiffXsec(shapename, 'x_Z|x_W.*|x_TauDecaysW', charge, outdir, isMu=True if flavour=="mu" else False)
    cmdMergeWithTestEffSyst = "hadd -f -k -O {of} {shapeNoTestEffSyst} {onlyTestEffSyst}".format(of=shapename,
                                                                                                 shapeNoTestEffSyst=shapenameNoTestEffSyst,
                                                                                                 onlyTestEffSyst=fileTestEffSyst)
    print "Now finally merging TestEffSystYY systematics into {of} ...".format(of=shapename)
    print cmdMergeWithTestEffSyst
    rmcmd = "rm {shapeNoTestEffSyst}".format(shapeNoTestEffSyst=shapenameNoTestEffSyst)
    print rmcmd
    if not options.dryrun: 
        os.system(cmdMergeWithTestEffSyst)
        os.system(rmcmd)
    print ""
    print "Wrote again root file in %s" % shapename
    print ""

# now adding prefiring for electrons
# this is on top of testEffSyst, i.e. it is an independent source of uncertainty
if flavour == "el":
    fileL1EffSyst = "{od}L1PrefireEleEffSyst_{fl}.root".format(od=outdir, fl=flavour)
    fileL1ZOutOfAcc = "{od}ZOutOfAccPrefireSyst_{fl}.root".format(od=outdir, fl=flavour)
    shapenameNoL1EffSyst = shapename.replace("_shapes","_shapes_noL1EffSyst")
    print "Copying {od} into {odnew}".format(od=shapename,odnew=shapenameNoL1EffSyst)
    if not options.dryrun: 
        os.system("cp {of} {ofnew}".format(of=shapename,ofnew=shapenameNoL1EffSyst))
        putEffSystHistos(shapename, 'x_Z|x_W.*|x_TauDecaysW', doType='L1PrefireEle', 
                         outdir=outdir, isMu=True if flavour=="mu" else False, isHelicity=False)
        addZOutOfAccPrefireSyst(shapename,outdir)
    cmdMergeWithL1EffSyst = "hadd -f -k -O {of} {shapeNoL1EffSyst} {onlyL1EffSyst} {L1ZOutOfAcc}".format(of=shapename, shapeNoL1EffSyst=shapenameNoL1EffSyst, onlyL1EffSyst=fileL1EffSyst, L1ZOutOfAcc=fileL1ZOutOfAcc)

    print "Now finally merging L1EffSystYY and L1ZOutOfAcc systematics into {of} ...".format(of=shapename)
    print cmdMergeWithL1EffSyst
    rmcmd = "rm {shapeNoL1EffSyst}".format(shapeNoL1EffSyst=shapenameNoL1EffSyst)
    print rmcmd
    if not options.dryrun: 
        os.system(cmdMergeWithL1EffSyst)
        os.system(rmcmd)
    print ""
    print "Wrote again root file in %s" % shapename
    print ""
# end of prefiring
####################

# test: symmetrize ratio of some systs over nominal, and redefine the alternate as nominal times symmetrized ratio 
# symmetrizing the ratio of some alternate histograms with nominal versus eta
if options.symSystRatio:
    match_symSystRatio = re.compile(options.symSystRatio)
    print "Now I will symmetrize versus eta the ratio of some alternate histograms with nominal"
    shapenameOnlySymSyst = shapename.replace("_shapes","_shapes_OnlySymSyst")
    tfno = ROOT.TFile.Open(shapename,"READ")
    if not tfno or not tfno.IsOpen():
        raise RuntimeError('Unable to open file {fn}'.format(fn=shapename))
    if not options.dryrun:
        # open output file
        of = ROOT.TFile(shapenameOnlySymSyst,'recreate')
        if not of or not of.IsOpen():
            raise RuntimeError('Unable to open file {fn}'.format(fn=shapenameOnlySymSyst))
        nKeys = tfno.GetNkeys()
        nCopiedKeys = 0
        for ikey,e in enumerate(tfno.GetListOfKeys()):
            name = e.GetName()
            if match_symSystRatio.match(name) and any(x in name for x in ["Up","Down"]):
                obj  = e.ReadObj()
                if not obj:
                    raise RuntimeError('Unable to read object {n}'.format(n=name))
                thisAlternate = obj.Clone(name+"_TMP_SYMMETRIZED")
                #thisAlternate.SetDirectory(0)
                thisNominalName = ""
                # remove last token, with syst name (WARNING, might not work for some systs)
                if "_CMS_W" in name:
                    thisNominalName = name.split("_CMS_W")[0]
                else:
                    thisNominalName = "_".join(name.split("_")[:-1]) 
                thisNominal = tfno.Get(thisNominalName)
                if not thisNominal:
                    print "Error when symmetrizing nuisances versus eta. I couldn't find nominal histogram %s. Abort" % thisNominalName
                    quit()
                else:
                    makeAlternateFromSymmetricRatioV2(thisAlternate, thisNominal, binning=recoBinning)
                    thisAlternate.Write(name)                
                nCopiedKeys += 1
            #else:
            #    newobj = obj.Clone(name+"_CLONE")
            #    newobj.Write(name)
            # nCopiedKeys += 1
            #sys.stdout.write('Key {num}/{tot}   \r'.format(num=ikey+1,tot=nKeys))
            #sys.stdout.flush()    
            if (ikey+1) % 100 == 0:                                      
                sys.stdout.write('Key {0:.2%}     \r'.format(float(ikey+1)/nKeys))
                sys.stdout.flush()
        print "Copied {n}/{tot} from {fn}".format(n=str(nCopiedKeys),tot=str(nKeys),fn=shapename)
        of.Close()
        tfno.Close()
    #print "Moving {od} into {odnew}".format(od=shapenameOnlySymSyst,odnew=shapename)

    # copy shapename into a new version which will be updated
    print "Now copying shapes into a backup file before overwriting some keys with symmetrized histograms"
    shapenameWithSymSyst = shapename.replace("_shapes","_shapes_WithSymSyst")
    cpcmd = "cp {od} {odnew}".format(odnew=shapenameWithSymSyst,od=shapename)
    print cpcmd
    if not options.dryrun:
        os.system(cpcmd)
    # now open shapenameWithSymSyst and overwrite keys which we just redefined
    tfno = ROOT.TFile.Open(shapenameOnlySymSyst,"READ")
    if not tfno or not tfno.IsOpen():
        raise RuntimeError('Unable to open file {fn}'.format(fn=shapenameOnlySymSyst))
    if not options.dryrun:
        # open output file
        of = ROOT.TFile(shapenameWithSymSyst,'UPDATE')
        if not of or not of.IsOpen():
            raise RuntimeError('Unable to open file {fn}'.format(fn=shapenameWithSymSyst))
        nKeys = tfno.GetNkeys()
        nCopiedKeys = 0
        of.cd()
        for ikey,e in enumerate(tfno.GetListOfKeys()):
            name = e.GetName()
            obj  = e.ReadObj()
            if not obj:
                raise RuntimeError('Unable to read object {n}'.format(n=name))
            obj.Write(name,ROOT.TObject.kOverwrite)
            nCopiedKeys += 1
        print "Overwrote {n}/{tot} from {fn}".format(n=str(nCopiedKeys),tot=str(nKeys),fn=shapenameOnlySymSyst)
        of.Close()
        tfno.Close()

    print "Now overwriting original file with backup including also symmetrized keys"
    mvcmd = "mv {od} {odnew}".format(od=shapenameWithSymSyst,odnew=shapename)
    print mvcmd
    if not options.dryrun:
        os.system(mvcmd)

print ""


## uncorrelate QCD scales and other possible nuisances by charge (add plus/minus as postfix to histograms' name)
# use the loop to do other things as well
regexp_uncorrCharge = options.uncorrelateNuisancesByCharge
if options.uncorrelateQCDscalesByCharge:
    if len(regexp_uncorrCharge):
        regexp_uncorrCharge += "|.*mu(R|F)\d+"
    else:
        regexp_uncorrCharge = ".*mu(R|F)\d+"
# skip this part when using options.useAnalyticSmoothPtScales, because the old pt scales are no longer needed
if not options.useAnalyticSmoothPtScales and (options.uncorrelatePtScaleByCharge and not options.useSmoothPtScales):
    if len(regexp_uncorrCharge):
        regexp_uncorrCharge += "|.*CMS_W.*scale\d+"
    else:
        regexp_uncorrCharge = ".*CMS_W.*scale\d+"    
if len(regexp_uncorrCharge) or options.distinguishNameSigAsBkg:
    match_ipt = re.compile("x_W{ch}.*_ipt_".format(ch=charge))
    match_regexp_uncorrCharge = re.compile(regexp_uncorrCharge)
    # loop on final merged file, look for muRXX or muFXX and add charge    
    print "Now I will uncorrelate between charges nuisances matching '{regexp}'".format(regexp=regexp_uncorrCharge)
    print "I will also do other work with histograms according to options"
    shapenameQCDuncorr = shapename.replace("_shapes","_shapes_QCDscalesChargeUncorr")
    #print "I will then copy it back into the original after uncorrelating QCD scales between charges" 
    tfno = ROOT.TFile.Open(shapename,"READ")
    if not tfno or not tfno.IsOpen():
        raise RuntimeError('Unable to open file {fn}'.format(fn=shapename))
    if not options.dryrun:
        # open output file
        of = ROOT.TFile(shapenameQCDuncorr,'recreate')
        if not of or not of.IsOpen():
            raise RuntimeError('Unable to open file {fn}'.format(fn=shapenameQCDuncorr))
        nKeys = tfno.GetNkeys()
        nCopiedKeys = 0
        for ikey,e in enumerate(tfno.GetListOfKeys()):
            name = e.GetName()
            obj  = e.ReadObj()
            if not obj:
                raise RuntimeError('Unable to read object {n}'.format(n=name))
            newname = name
            if options.distinguishNameSigAsBkg and hasPtRangeBkg:                
                # here change name from _lep_ to _el_ if needed
                if match_ipt.match(newname):
                    ipt = get_ipt_from_process_name(newname)
                    if ptBinIsBackground[ipt]:
                        newname = newname.replace("_lep_","_{lep}_".format(lep=flavour))
            if len(regexp_uncorrCharge) and match_regexp_uncorrCharge.match(newname):
                if newname.endswith('Up'): newname = newname.replace('Up','{ch}Up'.format(ch=charge))
                elif newname.endswith('Down'): newname = newname.replace('Down','{ch}Down'.format(ch=charge))
            newobj = obj.Clone(newname)
            newobj.Write(newname)        
            nCopiedKeys += 1
            #sys.stdout.write('Key {num}/{tot}   \r'.format(num=ikey+1,tot=nKeys))
            #sys.stdout.flush()
            if (ikey+1) % 100 == 0:                                      
                sys.stdout.write('Key {0:.2%}     \r'.format(float(ikey+1)/nKeys))
                sys.stdout.flush()
        print "Copied {n}/{tot} from {fn}".format(n=str(nCopiedKeys),tot=str(nKeys),fn=shapename)
        of.Close()
        tfno.Close()
    #print "Moving {od} into {odnew}".format(od=shapenameQCDuncorr,odnew=shapename)
    mvcmd = "mv {od} {odnew}".format(od=shapenameQCDuncorr,odnew=shapename)
    print mvcmd
    if not options.dryrun:
        os.system(mvcmd)

print ""
regexp_ptscale = ".*CMS_W.*scale\d+.*"
splitTag = "_CMS_W"
#if options.useSmoothPtScales:
#    regexp_ptscale = ".*smooth.*scale\d+.*"
#    splitTag = "_smooth"
#
# skip this part when using options.useAnalyticSmoothPtScales, because it no longer needs the old ptscales
#
if not options.useAnalyticSmoothPtScales and (options.uncorrelatePtScaleByEtaSide and not options.useSmoothPtScales):
    # if options.useSmoothPtScales is True this part was already managed in addSmoothLepScaleSyst
    print "Now I will uncorrelated PtScales by eta side"
    shapenameOnlyPtScaleEtaSideuncorr = shapename.replace("_shapes","_shapes_OnlyPtScaleEtaSideUncorr")
    #print "I will then copy it back into the original after uncorrelating QCD scales between charges" 
    tfno = ROOT.TFile.Open(shapename,"READ")
    if not tfno or not tfno.IsOpen():
        raise RuntimeError('Unable to open file {fn}'.format(fn=shapename))
    if not options.dryrun:
        # open output file
        of = ROOT.TFile(shapenameOnlyPtScaleEtaSideuncorr,'recreate')
        if not of or not of.IsOpen():
            raise RuntimeError('Unable to open file {fn}'.format(fn=shapenameOnlyPtScaleEtaSideuncorr))
        nKeys = tfno.GetNkeys()
        nCopiedKeys = 0
        for ikey,e in enumerate(tfno.GetListOfKeys()):
            name = e.GetName()
            newname = name
            if re.match(regexp_ptscale,name):                
                obj  = e.ReadObj()
                if not obj:
                    raise RuntimeError('Unable to read object {n}'.format(n=name))
                # roll to 2D, fill other side as nominal, then clone and change name. Repeat for each side
                nominalName = newname.split(splitTag)[0]
                hnomi = tfno.Get(nominalName)
                if not hnomi:
                    raise RuntimeError('Unable to retrieve nominal histogram named {fn}'.format(fn=nominalName))
                if   newname.endswith('Up'):   newname = newname.replace('Up','etaSIDEUp')
                elif newname.endswith('Down'): newname = newname.replace('Down','etaSIDEDown')
                hetaPos = obj.Clone(newname.replace("SIDE","sideP"))
                hetaNeg = obj.Clone(newname.replace("SIDE","sideM"))
                uncorrelateHistByEtaSide(hnomi, hetaPos, hetaNeg, recoBins.Neta, recoBins.Npt)
                hetaPos.Write(hetaPos.GetName())
                hetaNeg.Write(hetaNeg.GetName())    
                nCopiedKeys += 2
            #else:
            #    newobj = obj.Clone(newname)
            #    newobj.Write(newname)        
            #    nCopiedKeys += 1
            #sys.stdout.write('Key {num}/{tot}   \r'.format(num=ikey+1,tot=nKeys))
            #sys.stdout.flush()
            if (ikey+1) % 100 == 0:                                      
                sys.stdout.write('Key {0:.2%}     \r'.format(float(ikey+1)/nKeys))
                sys.stdout.flush()
        print "Copied {n}/{tot} from {fn}".format(n=str(nCopiedKeys),tot=str(nKeys),fn=shapename)
        of.Close()
        tfno.Close()

    # copy shapename into a new version which will be updated
    print "Now copying shapes into a backup file before overwriting some keys with scale histograms"
    shapenameWithUncorrScale = shapename.replace("_shapes","_shapes_WithUncorrScale")
    cpcmd = "cp {od} {odnew}".format(odnew=shapenameWithUncorrScale,od=shapename)
    print cpcmd
    if not options.dryrun:
        os.system(cpcmd)
    # now open shapenameWithUncorrScale and overwrite keys which we just redefined
    tfno = ROOT.TFile.Open(shapenameOnlyPtScaleEtaSideuncorr,"READ")
    if not tfno or not tfno.IsOpen():
        raise RuntimeError('Unable to open file {fn}'.format(fn=shapenameOnlyPtScaleEtaSideuncorr))
    if not options.dryrun:
        # open output file
        of = ROOT.TFile(shapenameWithUncorrScale,'UPDATE')
        if not of or not of.IsOpen():
            raise RuntimeError('Unable to open file {fn}'.format(fn=shapenameWithUncorrScale))
        nKeys = tfno.GetNkeys()
        nCopiedKeys = 0
        of.cd()
        for ikey,e in enumerate(tfno.GetListOfKeys()):
            name = e.GetName()
            obj  = e.ReadObj()
            if not obj:
                raise RuntimeError('Unable to read object {n}'.format(n=name))
            obj.Write(name,ROOT.TObject.kOverwrite)
            nCopiedKeys += 1
        print "Overwrote {n}/{tot} from {fn}".format(n=str(nCopiedKeys),tot=str(nKeys),fn=shapenameOnlyPtScaleEtaSideuncorr)
        of.Close()
        tfno.Close()

    print "Now overwriting original file with backup including also new uncorrelated scales keys"
    #print "Moving {od} into {odnew}".format(od=shapenameWithUncorrScale,odnew=shapename)
    mvcmd = "mv {od} {odnew}".format(od=shapenameWithUncorrScale,odnew=shapename)
    print mvcmd
    if not options.dryrun:
        os.system(mvcmd)


print ""
print "-"*20
print ""


# now smooth momentum scales, obtained with addSmoothLepScaleSyst()
# nuis names smooth{lep}scale{idx}: lep=mu/el, idx goes from 0 to N-1, N=2(4) for mu(el)
# these are added on top of the unsmoothed ones (the standard ones).
# one can choose to define the extreme pt bins from the standard ones (this algorithm do not work for those bins)
# that's why we keep this step as the last one (we already need all the standard momentum scales)
if not options.useAnalyticSmoothPtScales and (options.useSmoothPtScales or options.addSmoothPtScalesWithStandardOnExtremePtBins):
    print "Now adding smooth pt scales"
    if options.uncorrelatePtScaleByCharge:
        print "pt scales will be uncorrelated between charges"
    if options.uncorrelatePtScaleByEtaSide:
        print "pt scales will be uncorrelated between eta sides"
    fileSmoothPtScales = "{od}SmoothScaleSyst_{fl}_{ch}.root".format(od=outdir, fl=flavour, ch=options.charge)
    shapenameNoSmoothPtScales = shapename.replace("_shapes","_shapes_noSmoothPtScales")
    print "Copying {od} into {odnew}".format(od=shapename,odnew=shapenameNoSmoothPtScales)
    if not options.dryrun: 
        os.system("cp {of} {ofnew}".format(of=shapename,ofnew=shapenameNoSmoothPtScales))
        addSmoothLepScaleSyst(shapename, 'x_Z|x_W.*|x_TauDecaysW', options.charge,
                              isMu=True if flavour=="mu" else False,  
                              outdir=outdir, 
                              uncorrelateByCharge=options.uncorrelatePtScaleByCharge,
                              uncorrelateByEtaSide=options.uncorrelatePtScaleByEtaSide,
                              getExtremePtBinsFromOldScale=True if options.addSmoothPtScalesWithStandardOnExtremePtBins else False)
    cmdMergeWithSmoothPtScales = "hadd -f -k -O {of} {shapeNoSmoothPtScales} {onlySmoothPtScales}".format(of=shapename, shapeNoSmoothPtScales=shapenameNoSmoothPtScales, onlySmoothPtScales=fileSmoothPtScales)

    print "Now finally merging smooth(mu|el)scaleYY systematics into {of} ...".format(of=shapename)
    print cmdMergeWithSmoothPtScales
    rmcmd = "rm {shapeNoSmoothPtScales}".format(shapeNoSmoothPtScales=shapenameNoSmoothPtScales)
    print rmcmd
    if not options.dryrun: 
        os.system(cmdMergeWithSmoothPtScales)
        os.system(rmcmd)
    print ""
    print "Wrote again root file in %s" % shapename
    print ""



# these are the ultimate momentum scales, made using only the analytic smoothing 
# they only need access to additional nominal histograms made with extended pt acceptance
# for muons, use Rochester corrections
if options.useAnalyticSmoothPtScales:
    print "Now adding smooth pt scales (made moving events based on analytic pt-dependent shift)"
    print "pt scales will be uncorrelated between charges"
    fileSmoothPtScales = ""
    if flavour=="mu":
        fileSmoothPtScales = "{od}SmoothScaleSystRochester_{fl}_{ch}.root".format(od=outdir, fl=flavour, ch=options.charge)  
    else:
        fileSmoothPtScales = "{od}SmoothScaleSystClosure_{fl}_{ch}.root".format(od=outdir, fl=flavour, ch=options.charge) 
    shapenameNoSmoothPtScales = shapename.replace("_shapes","_shapes_noSmoothPtScales")
    print "Copying {od} into {odnew}".format(od=shapename,odnew=shapenameNoSmoothPtScales)
    if not options.dryrun: 
        os.system("cp {of} {ofnew}".format(of=shapename,ofnew=shapenameNoSmoothPtScales))
        if flavour=="mu":
            addSmoothMuonScaleSyst(shapename, 'x_Z|x_W.*|x_TauDecaysW', options.charge, outdir=outdir)
        else:
            addSmoothElectronScaleSyst(shapename, 'x_Z|x_W.*|x_TauDecaysW', options.charge, outdir=outdir,
                                       uncorrelateByCharge=options.uncorrelatePtScaleByCharge,
                                       uncorrelateByEtaSide=options.uncorrelatePtScaleByEtaSide)

    cmdMergeWithSmoothPtScales = "hadd -f -k -O {of} {shapeNoSmoothPtScales} {onlySmoothPtScales}".format(of=shapename, shapeNoSmoothPtScales=shapenameNoSmoothPtScales, onlySmoothPtScales=fileSmoothPtScales)

    print "Now finally merging smooth(mu|el)scale(Stat|Syst)YY systematics into {of} ...".format(of=shapename)
    print cmdMergeWithSmoothPtScales
    rmcmd = "rm {shapeNoSmoothPtScales}".format(shapeNoSmoothPtScales=shapenameNoSmoothPtScales)
    print rmcmd
    if not options.dryrun: 
        os.system(cmdMergeWithSmoothPtScales)
        os.system(rmcmd)
    print ""
    print "Wrote again root file in %s" % shapename
    print ""




# can also use
# python changeNameSignalProcess.py -i cards/diffXsec_el_2019_07_20_latestScaleFactor_AllIn_IDwithMConlyStat/ -o cards/diffXsec_el_2019_07_20_latestScaleFactor_AllIn_IDwithMConlyStat/ --pt-bins-name-change "0,1" -r "_lep_" "_el_" -f el -c plus; python changeNameSignalProcess.py -i cards/diffXsec_el_2019_07_20_latestScaleFactor_AllIn_IDwithMConlyStat/ -o cards/diffXsec_el_2019_07_20_latestScaleFactor_AllIn_IDwithMConlyStat/ --pt-bins-name-change "0,1" -r "_lep_" "_el_" -f el -c minus
