#!/usr/bin/env python                                                       

# to add smooth pt scales in a file that contains the standard version (used to fix the extreme pt bins)
 
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
from w_helicity_13TeV.make_diff_xsec_cards import etaIsOutsideFilledRangeDiffXsec

from w_helicity_13TeV.templateRolling import roll1Dto2D, dressed2D
from w_helicity_13TeV.rollingFunctions import unroll2Dto1D

from w_helicity_13TeV.mergeCardComponentsAbsY import mirrorShape
from w_helicity_13TeV.mergeCardComponentsAbsY import putUncorrelatedFakes
from w_helicity_13TeV.mergeCardComponentsAbsY import putEffSystHistos # for electron L1-prefire uncertainty
from w_helicity_13TeV.mergeCardComponentsAbsY import addZOutOfAccPrefireSyst # for ele L1-prefire uncertainty on Z
from w_helicity_13TeV.mergeCardComponentsAbsY import cropHighSysts

from cropNegativeTemplateBins import cropNegativeContent

# copied from helicity merger, might not be ideal for 2D xsec
def addSmoothLeptonScaleSyst_algoHelicity(infile,regexp,charge,isMu,alternateShapeOnly=False,outdir=None):

    indir = os.path.dirname(os.path.abspath(infile))
    outdir = outdir if outdir != None else indir
    flav = 'mu' if isMu else 'el'

    # get eta-pt binning for reco 
    etaPtBinningVec = getDiffXsecBinning(indir+'/binningPtEta.txt', "reco")  # this get two vectors with eta and pt binning
    recoBins = templateBinning(etaPtBinningVec[0],etaPtBinningVec[1])        # this create a class to manage the binnings
    binning = [recoBins.Neta, recoBins.etaBins, recoBins.Npt, recoBins.ptBins]

    tmp_infile = ROOT.TFile(infile, 'read')
    outfile = ROOT.TFile(outdir+'/SmoothScaleSyst_{flav}_{ch}.root'.format(flav=flav,ch=charge), 'recreate')

    dirWithInputSyst = '/afs/cern.ch/work/e/emanuele/wmass/heppy/CMSSW_8_0_25/src/CMGTools/WMass/python/postprocessing/data/leptonScale'
    muscale_syst_f = ROOT.TFile.Open(dirWithInputSyst+'/mu/muscales.root') if isMu else ROOT.TFile.Open(dirWithInputSyst+'/el/elscales_extended.root')
    ## use one histogram to map the binning -> array
    binning_histo = muscale_syst_f.Get('stathis_eig_plus_0')
    ## convert histograms to arrays to get it faster
    maxstats = 99 if isMu else 97
    systrange = [2,5] if isMu else [0,1]
    stathists = [root_numpy.hist2array(muscale_syst_f.Get('stathis_eig_{ch}_{istat}'.format(ch=charge,istat=idx))) for idx in range(maxstats)]
    systhists = [root_numpy.hist2array(muscale_syst_f.Get('systhist_{ch}_{isyst}'.format(ch=charge,isyst=idx))) for idx in range(systrange[0],systrange[1]+1)]

    ## stat error from re-generated stat. replicas (diagonalized) 
    ## then 4 systematic uncertainties fully correlated in the eta/pt plane
    ## for muons, keep the original rochester corr naming Syst2-Syst6
    offset = 2 if isMu else 0
    allsysts = ['Stat{idx}'.format(idx=i) for i in range(len(stathists))] + ['Syst{idx}'.format(idx=i+offset) for i in range(len(systhists))]
    allhists = stathists + systhists
    systsAndHists = zip(allsysts,allhists)

    for k in tmp_infile.GetListOfKeys():
        tmp_name = k.GetName()
        ## don't reweight any histos that don't match the regexp
        if not re.match(regexp, tmp_name): continue
        ## don't reweight any histos that are already variations of something else
        if 'Up' in tmp_name or 'Down' in tmp_name: continue
        process = tmp_name.split('_')[1]

        ## now should be left with only the ones we are interested in
        print 'reweighting ',flav,' scale syst of type for process', tmp_name
        
        tmp_nominal = tmp_infile.Get(tmp_name)
        tmp_nominal_2d = dressed2D(tmp_nominal,binning, tmp_name+'backrolled')
        n_ptbins = tmp_nominal_2d.GetNbinsY()


        for syst,hist in systsAndHists:
            if re.match('Syst',syst):
                syst_ptbins = [('',0,1000)] if isMu else [('pt0',0,42),('pt1',42,1000)]
            else:
                syst_ptbins = [('',0,1000)]
            for systipt in syst_ptbins:
                for shift_dir in ['Up','Down']:
                    outname_2d = tmp_nominal_2d.GetName().replace('backrolled','')+'_smooth{lep}scale{idx}{systipt}{shiftdir}'.format(lep=flav,idx=syst,systipt=systipt[0],shiftdir=shift_dir)
                    tmp_scaledHisto = copy.deepcopy(tmp_nominal_2d.Clone(outname_2d))
     
                    ## loop over all eta bins of the 2d histogram 
                    for ieta in range(1,tmp_nominal_2d.GetNbinsX()+1):
                        eta = tmp_nominal_2d.GetXaxis().GetBinCenter(ieta)
                        for ipt in range(2,tmp_nominal_2d.GetNbinsY()+1):
                            pt      = tmp_nominal_2d.GetYaxis().GetBinCenter(ipt)
                            if systipt[1]<= pt < systipt[2]:
                                ## assume uniform distribution within a bin
                                etabin = max(1, min(binning_histo.GetNbinsX(), binning_histo.GetXaxis().FindFixBin(eta)))
                                ptbin  = max(1, min(binning_histo.GetNbinsY(), binning_histo.GetYaxis().FindFixBin(pt)))
         
                                pt_prev = tmp_nominal_2d.GetYaxis().GetBinCenter(ipt-1)
                                nominal_val      = tmp_nominal_2d.GetBinContent(ieta,ipt)
                                nominal_val_prev = tmp_nominal_2d.GetBinContent(ieta,ipt-1)
                                pt_width      = tmp_nominal_2d.GetYaxis().GetBinWidth(ipt)
                                pt_width_prev = tmp_nominal_2d.GetYaxis().GetBinWidth(ipt-1)
             
                                if isMu: scale_syst = 1-utilities.getRochesterUncertainty(charge,etabin-1,ptbin-1,hist,syst!='Syst3')
                                else:    scale_syst = 1-utilities.getRochesterUncertainty(charge,etabin-1,ptbin-1,hist,False)
                                if shift_dir=='Down': scale_syst = -1*scale_syst
         
                                from_prev = scale_syst * pt_prev / pt_width_prev * max(0,nominal_val_prev)
                                to_right  = scale_syst * pt      / pt_width      * max(0,nominal_val)
         
                                tmp_scaledHisto.SetBinContent(ieta,ipt,max(0,nominal_val + from_prev - to_right))
                        ## since in the first bin we cannot foresee 2-neighbors migrations, let's assume same syst of bin i+1
                        tmp_scaledHisto.SetBinContent(ieta,1,tmp_scaledHisto.GetBinContent(ieta,2)/tmp_nominal_2d.GetBinContent(ieta,2)*tmp_nominal_2d.GetBinContent(ieta,1) if tmp_nominal_2d.GetBinContent(ieta,2) else 0)
                    ## re-roll the 2D to a 1D histo
                    tmp_scaledHisto_1d = unroll2Dto1D(tmp_scaledHisto, newname=tmp_scaledHisto.GetName().replace('2DROLLED',''))
                    if alternateShapeOnly:
                        tmp_scaledHisto_1d.Scale(tmp_nominal.Integral()/tmp_scaledHisto_1d.Integral())
                    cropHighSysts(tmp_nominal,tmp_scaledHisto_1d,maxSyst=0.05)                
                    outfile.cd()
                    tmp_scaledHisto_1d.Write()
    outfile.Close()
    print "done with the smooth ",flav," scale variations"


# function for both muons and ele (muon part copied from the other function below
# use Rochester corrections for muons, and analogue systematics for electrons
def addSmoothLeptonScaleSyst(infile, regexp, charge, isMu,
                             outdir=None, 
                             alternateShapeOnly=False,
                             cropNegativeBin=False,  # set negative bin to 0 in both nominal and alternate
                             maxSyst=0.0):  # if > 0, set max/min ratio to 1+/- this value

    indir = os.path.dirname(os.path.abspath(infile))
    outdir = outdir if outdir != None else indir
    flav = 'mu' if isMu else 'el'
    lepton = 'muon' if isMu else 'electron'

    # get eta-pt binning for reco 
    etaPtBinningVec = getDiffXsecBinning(indir+'/binningPtEta.txt', "reco")  # this get two vectors with eta and pt binning
    recoBins = templateBinning(etaPtBinningVec[0],etaPtBinningVec[1])        # this create a class to manage the binnings
    binning = [recoBins.Neta, recoBins.etaBins, recoBins.Npt, recoBins.ptBins]
    etaPtBinningVecGen = getDiffXsecBinning(indir+'/binningPtEta.txt', "gen")  # this get two vectors with eta and pt binning
    genBins = templateBinning(etaPtBinningVecGen[0],etaPtBinningVecGen[1])        # this create a class to manage the binnings
    tmp_infile = ROOT.TFile(infile, 'read')
    fnameTag = "SmoothScaleSystRochester" if isMu else "SmoothScaleSystElectron"
    outfile = ROOT.TFile(outdir+'/{n}_{flav}_{ch}.root'.format(n=fnameTag,flav=flav,ch=charge), 'recreate')

    ##//////////////////////////
    # get nominal histograms with one more pt bin outside acceptance, to do smoothing
    #
    # define new binning for pt
    extendedPtBins = []
    if isMu:
        extendedPtBins = [25.0] + [float(x) for x in recoBins.ptBins] + [57.0]
    else:
        extendedPtBins = [29.0] + [float(x) for x in recoBins.ptBins] + [57.0]

    # now open file and get histograms
    # these are TH1 and should be rolled to TH2
    specialFolderTag = ""
    if "_1sigBin_4fixedPOI" in infile:
        specialFolderTag = "_1sigBin_4fixedPOI"
    
    file_nomi_recoPt25to57 = "cards/diffXsec_mu_2019_09_19_onlyZandTau_recoPt25to57{sft}/nominalShape_{ch}_recoPt25to57.root".format(ch=charge,sft=specialFolderTag)
    file_nomi_recoPt29to57 = "cards/diffXsec_el_2019_09_22_onlyZandTau_recoPt29to57{sft}/nominalShape_{ch}_recoPt29to57.root".format(ch=charge,sft=specialFolderTag)

    if "_1sigBin_4fixedPOI_ptMax45" in infile:
        file_nomi_recoPt25to57 = "cards/diffXsec_mu_2019_09_19_onlyZandTau_recoPt25to46p5_1sigBin_4fixedPOI/nominalShape_{ch}_recoPt25to46p5.root".format(ch=charge,sft=specialFolderTag)
        file_nomi_recoPt29to57 = "cards/diffXsec_el_2019_09_22_onlyZandTau_recoPt29to46p5_1sigBin_4fixedPOI/nominalShape_{ch}_recoPt29to46p5.root".format(ch=charge,sft=specialFolderTag)
        extendedPtBins[-1] = 46.5
        

    # it is 29-57 for electrons
    binning_recoPt2Xto57 = [recoBins.Neta, recoBins.etaBins, len(extendedPtBins)-1, extendedPtBins]

    hist_nomi_recoPt2Xto57 = ROOT.TFile.Open(file_nomi_recoPt25to57 if isMu else file_nomi_recoPt29to57)
    nominalHists = {}
    for ikey,e in enumerate(hist_nomi_recoPt2Xto57.GetListOfKeys()):
        name = e.GetName()
        obj  = e.ReadObj()
        if not obj:
            raise RuntimeError('Unable to read object {n}'.format(n=name))
        if cropNegativeBin: cropNegativeContent(obj,silent=False)
        nominalHists[name] = dressed2D(obj,binning_recoPt2Xto57, name+"_recoPt2Xto57_backrolled")
        nominalHists[name].SetDirectory(0)
    hist_nomi_recoPt2Xto57.Close()
    # these will be used later
    #
    ##//////////////////////////
    #muscale_syst_f = ROOT.TFile.Open('../postprocessing/data/leptonScale/mu/muscales.root')
    #muscale_syst_f = ROOT.TFile.Open('/afs/cern.ch/user/b/bendavid/cmspublic/muscales_extended.root')
    lepfile = "/afs/cern.ch/work/e/emanuele/wmass/heppy/CMSSW_8_0_25/src/CMGTools/WMass/python/postprocessing/data/leptonScale/mu/muscales.root" if isMu else "/afs/cern.ch/work/e/emanuele/wmass/heppy/CMSSW_8_0_25/src/CMGTools/WMass/python/postprocessing/data/leptonScale/el/elscales_extended.root"
    lepscale_syst_f = ROOT.TFile.Open(lepfile)
    ## use one histogram to map the binning -> array
    binning_histo = lepscale_syst_f.Get('stathis_eig_plus_0')
    ## convert histograms to arrays to get it faster
    maxstats = 99 if isMu else 97
    systrange = [2,5] if isMu else [0,1]
    stathists = [root_numpy.hist2array(lepscale_syst_f.Get('stathis_eig_{ch}_{istat}'.format(ch=charge,istat=idx))) for idx in range(maxstats)]
    systhists = [root_numpy.hist2array(lepscale_syst_f.Get('systhist_{ch}_{isyst}'.format(ch=charge,isyst=idx))) for idx in range(systrange[0],systrange[1]+1)]

    ## stat error from re-generated stat. replicas (diagonalized) 
    ## then 4 systematic uncertainties fully correlated in the eta/pt plane
    offset = 2 if isMu else 0
    allsysts = ['Stat{idx}'.format(idx=i) for i in range(len(stathists))] + ['Syst{idx}'.format(idx=i+offset) for i in range(len(systhists))]
    allhists = stathists + systhists
    systsAndHists = zip(allsysts,allhists)
    ##//////////////////////////

    for k in tmp_infile.GetListOfKeys():
        tmp_name = k.GetName()
        ## don't reweight any histos that don't match the regexp
        if not re.match(regexp, tmp_name): continue
        ## don't reweight any histos that are already variations of something else
        if 'Up' in tmp_name or 'Down' in tmp_name: continue
        process = tmp_name.split('_')[1] 

        ## now should be left with only the ones we are interested in
        print 'reweighting smooth %sscale syst for process %s' % (flav,tmp_name)
        
        tmp_nominal = tmp_infile.Get(tmp_name)
        if cropNegativeBin: cropNegativeContent(tmp_nominal,silent=False)
        tmp_nominal_2d = dressed2D(tmp_nominal,binning, tmp_name+'backrolled')
        n_ptbins = tmp_nominal_2d.GetNbinsY()

        for syst,hist in systsAndHists:
            if re.match('Syst',syst):                
                #syst_ptbins = [('',0,1000)] if isMu else [('pt0',0,42),('pt1',42,1000)]
                syst_ptbins = [('',0,1000)]
                if isMu:
                    syst_etabins = [('etaside1M',-2.4,-2.1),('etaside0M',-2.1,0),('etaside0P',0,2.1),('etaside1P',2.1,2.4)] 
                else:
                    # for 2D xsec the range goes up to 2.4, unlike for the helicity measurement
                    syst_etabins = [('etaside3M',-2.4,-2.1),('etaside2M',-2.1,-1.5),('etaside1M',-1.5,-1),('etaside0M',-1,0),('etaside0P',0,1),('etaside1P',1,1.5),('etaside2P',1.5,2.1),('etaside3P',2.1,2.4)]
            else:
                syst_ptbins = [('',0,1000)]
                syst_etabins = [('',-10,10)]
            for systipt in syst_ptbins:
                for systieta in syst_etabins:
                    # for syst, when syst_etabins has more than 1 item, check if signal bin matches
                    if len(syst_etabins) > 1 and "_ieta_" in tmp_name: 
                        #print "===== check ====="
                        gen_ieta = get_ieta_from_process_name(tmp_name)
                        maxEtaAbs = max(abs(systieta[1]),abs(systieta[2]))
                        minEtaAbs = min(abs(systieta[1]),abs(systieta[2]))
                        centerGenEta = 0.5 * (genBins.etaBins[gen_ieta] + genBins.etaBins[gen_ieta+1])
                        if centerGenEta < minEtaAbs or centerGenEta > maxEtaAbs:
                            continue

                    for shift_dir in ['Up','Down']:
                        chargesyst = '' if ((isMu and syst=='Syst3') or re.match('Stat',syst)) else charge
                        outname_2d = tmp_nominal_2d.GetName().replace('backrolled','')+'_smooth{lep}scale{idx}{systieta}{systipt}{ch}{shiftdir}'.format(lep=flav,idx=syst,systieta=systieta[0],systipt=systipt[0],ch=chargesyst,shiftdir=shift_dir)
                        tmp_scaledHisto = copy.deepcopy(tmp_nominal_2d.Clone(outname_2d))
                        tmp_scaledHisto.Reset()

                        ## loop over all eta bins of the 2d histogram 
                        ## for signal, can just do a single eta bin and its neighbors (the rest is empty)
                        for ieta in range(1,tmp_nominal_2d.GetNbinsX()+1):
                            eta = tmp_nominal_2d.GetXaxis().GetBinCenter(ieta)
                            if etaIsOutsideFilledRangeDiffXsec(tmp_name,eta,
                                                               NetaNeighbours=1,
                                                               genEtaBins=genBins.etaBins):
                                continue
                            # if re.match(".*W.*_ieta_.*",tmp_name):
                            #     genetabin = get_ieta_from_process_name(tmp_name)
                            #     # check range including first neighbours as well
                            #     # careful: genetabin goes from 0 to N(genBins)-1
                            #     if abs(eta) < genBins.etaBins[max(0,genetabin-1)] or abs(eta) >= genBins.etaBins[min(genBins.Neta,genetabin+2)]:
                            #         continue
                            # try to save some time for muons, most histograms only need eta, not pt (average is made)
                            etabin = max(1, min(binning_histo.GetNbinsX(), binning_histo.GetXaxis().FindFixBin(eta)))
                            tmp_syst_ptscale = 0
                            if isMu and syst!='Syst3': 
                                tmp_syst_ptscale = utilities.getRochesterUncertainty(charge,etabin-1,0,hist,True)  
                            for ipt in range(1,tmp_nominal_2d.GetNbinsY()+1):                            

                                ptcenter      = tmp_nominal_2d.GetYaxis().GetBinCenter(ipt)
                                ptbin  = max(1, min(binning_histo.GetNbinsY(), binning_histo.GetYaxis().FindFixBin(ptcenter)))
                                # scale_syst is something close to 1: if larger, we move events from left to right
                                # if lower, no event will move from left to right, but actually from right to left
                                syst_ptscale = 1 # it can be seen as 1 + dx with dx > or < 0
                                # the syst can be > or < 1, we define the variation based on the number in the histogram, and arbitrarily call Up the variation based on that.
                                # then we artificially define as Down the variation defined as 1 - (var-1) = 2 -var
                                # generally, the variation would be > 1 (which would mean scaling pt up
                                # but in some regions inside acceptance it could be < 1 ( but we still call it up when creating the alternate histogram)
                                #if (not isMu) or (isMu and syst == 'Syst3'): 
                                nominal_val = tmp_nominal_2d.GetBinContent(ieta,ipt)
                                pt_width    = tmp_nominal_2d.GetYaxis().GetBinWidth(ipt)

                                if (systipt[1]<= ptcenter < systipt[2]) and (systieta[1]<= eta < systieta[2]):
                                    if not (isMu and syst!='Syst3'):
                                        tmp_syst_ptscale = utilities.getRochesterUncertainty(charge,etabin-1,ptbin-1,hist,False)
                                    if shift_dir == "Up":
                                        syst_ptscale = tmp_syst_ptscale 
                                    else:
                                        syst_ptscale = 2. - tmp_syst_ptscale

                                    if syst_ptscale > 1:                            
                                        scale_syst = syst_ptscale - 1. 
                                        pt      = tmp_nominal_2d.GetYaxis().GetBinUpEdge(ipt)
                                        pt_prev = tmp_nominal_2d.GetYaxis().GetBinLowEdge(ipt)
                                        if ipt == 1:
                                            # take first pt bin from histo with extended pt range
                                            nameHisto = tmp_name if isMu else tmp_name.replace("_el_","_lep_")
                                            nominal_val_prev = nominalHists[nameHisto].GetBinContent(ieta,1)
                                            pt_width_prev    = nominalHists[nameHisto].GetYaxis().GetBinWidth(1)   
                                        else:
                                            nominal_val_prev = tmp_nominal_2d.GetBinContent(ieta,ipt-1)
                                            pt_width_prev = tmp_nominal_2d.GetYaxis().GetBinWidth(ipt-1)
                                        factor = scale_syst / (1.0 + scale_syst)  # > 0
                                        from_prev = factor * (pt_prev / pt_width_prev) * max(0,nominal_val_prev)
                                        to_next   = factor * (pt      / pt_width     ) * max(0,nominal_val)
                                        tmp_scaledHisto.SetBinContent(ieta,ipt,max(0,nominal_val + from_prev - to_next))
                                    else:
                                        scale_syst = 1 - syst_ptscale # so scale_syst > 0 again
                                        pt      = tmp_nominal_2d.GetYaxis().GetBinLowEdge(ipt)
                                        pt_next = tmp_nominal_2d.GetYaxis().GetBinUpEdge(ipt)
                                        if ipt == tmp_nominal_2d.GetNbinsY():
                                            # histograms with extended pt range have 2 more pt bins
                                            nameHisto = tmp_name if isMu else tmp_name.replace("_el_","_lep_")
                                            nominal_val_next = nominalHists[nameHisto].GetBinContent(ieta,ipt+2)
                                            pt_width_next    = nominalHists[nameHisto].GetYaxis().GetBinWidth(ipt+2)
                                        else:
                                            nominal_val_next = tmp_nominal_2d.GetBinContent(ieta,ipt+1)
                                            pt_width_next    = tmp_nominal_2d.GetYaxis().GetBinWidth(ipt+1)
                                        factor = scale_syst / (1.0 - scale_syst)
                                        from_next = factor * (pt_next / pt_width_next) * max(0,nominal_val_next)
                                        to_prev   = factor * (pt      / pt_width     ) * max(0,nominal_val)
                                        tmp_scaledHisto.SetBinContent(ieta,ipt,max(0,nominal_val + from_next - to_prev))
                                else:
                                    tmp_scaledHisto.SetBinContent(ieta,ipt,nominal_val)


                            # the following is not needed if cropping bins
                            # if maxSyst == 0.0 and "outliers" in tmp_name and shift_dir == "Up":
                            #     # ad hoc fix for second bin, which gets a very large ratio because the assumption
                            #     # of flat event density is maximally wrong (would be the same for first bin but let's neglect that)
                            #     nomiContent = tmp_nominal_2d.GetBinContent(ieta,2)
                            #     ratioPt3 = 1.0
                            #     if tmp_nominal_2d.GetBinContent(ieta,3) != 0:
                            #         ratioPt3 = tmp_scaledHisto.GetBinContent(ieta,3)/tmp_nominal_2d.GetBinContent(ieta,3)
                            #     tmp_scaledHisto.SetBinContent(ieta,2,ratioPt3*nomiContent)

                        ## re-roll the 2D to a 1D histo
                        tmp_scaledHisto_1d = unroll2Dto1D(tmp_scaledHisto, newname=tmp_scaledHisto.GetName().replace('2DROLLED',''), cropNegativeBins=cropNegativeBin)

                        if abs(maxSyst) > 0.0:
                            cropHighSysts(tmp_nominal,tmp_scaledHisto_1d,maxSyst=maxSyst)
                        outfile.cd()
                        tmp_scaledHisto_1d.Write()

    outfile.Close()
    print "done with the smooth %s scale variations" % lepton
    print "Output file: %s" % outfile.GetName()


#############


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


# experimental function for muon momentum scales using Rochester corrections (ok)
def addSmoothMuonScaleSyst(infile, regexp, charge, 
                           outdir=None, 
                           alternateShapeOnly=False):

    indir = os.path.dirname(os.path.abspath(infile))
    outdir = outdir if outdir != None else indir
    flav = 'mu'
    isMu = True

    # get eta-pt binning for reco 
    etaPtBinningVec = getDiffXsecBinning(indir+'/binningPtEta.txt', "reco")  # this get two vectors with eta and pt binning
    recoBins = templateBinning(etaPtBinningVec[0],etaPtBinningVec[1])        # this create a class to manage the binnings
    binning = [recoBins.Neta, recoBins.etaBins, recoBins.Npt, recoBins.ptBins]

    tmp_infile = ROOT.TFile(infile, 'read')
    outfile = ROOT.TFile(outdir+'/SmoothScaleSystRochester_{flav}_{ch}.root'.format(flav=flav,ch=charge), 'recreate')

    ##//////////////////////////
    # get nominal histograms with one more pt bin outside acceptance, to do smoothing
    #
    # define new binning for pt
    extendedPtBins = [25.0] + [float(x) for x in recoBins.ptBins] + [57.0]
    binning_recoPt25to57 = [recoBins.Neta, recoBins.etaBins, len(extendedPtBins)-1, extendedPtBins]
    # now open file and get histograms
    # these are TH1 and should be rolled to TH2
    file_nomi_recoPt25to57 = "cards/diffXsec_mu_2019_09_19_onlyZandTau_recoPt25to57/nominalShape_{ch}_recoPt25to57.root".format(ch=charge)
    hist_nomi_recoPt25to57 = ROOT.TFile.Open(file_nomi_recoPt25to57)
    nominalHists = {}
    for ikey,e in enumerate(hist_nomi_recoPt25to57.GetListOfKeys()):
        name = e.GetName()
        obj  = e.ReadObj()
        if not obj:
            raise RuntimeError('Unable to read object {n}'.format(n=name))
        nominalHists[name] = dressed2D(obj,binning_recoPt25to57, name+"_recoPt25to57_backrolled")
        nominalHists[name].SetDirectory(0)
    hist_nomi_recoPt25to57.Close()
    # these will be used later
    #
    ##//////////////////////////
    #muscale_syst_f = ROOT.TFile.Open('../postprocessing/data/leptonScale/mu/muscales.root')
    #muscale_syst_f = ROOT.TFile.Open('/afs/cern.ch/user/b/bendavid/cmspublic/muscales_extended.root')
    muscale_syst_f = ROOT.TFile.Open('/afs/cern.ch/work/e/emanuele/wmass/heppy/CMSSW_8_0_25/src/CMGTools/WMass/python/postprocessing/data/leptonScale/mu/muscales.root')
    ## use one histogram to map the binning -> array
    binning_histo = muscale_syst_f.Get('systhist_plus_2')
    ## convert histograms to arrays to get it faster
    stathists = [root_numpy.hist2array(muscale_syst_f.Get('stathis_eig_{ch}_{istat}'.format(ch=charge,istat=idx))) for idx in range(99)]
    systhists = [root_numpy.hist2array(muscale_syst_f.Get('systhist_{ch}_{isyst}'.format(ch=charge,isyst=idx))) for idx in range(2,6)]

    ## stat error from re-generated stat. replicas (diagonalized) 
    ## then 4 systematic uncertainties fully correlated in the eta/pt plane
    allsysts = ['Stat{idx}'.format(idx=i) for i in range(len(stathists))] + ['Syst{idx}'.format(idx=i+2) for i in range(len(systhists))]
    allhists = stathists + systhists
    systsAndHists = zip(allsysts,allhists)
    ##//////////////////////////

    
    for k in tmp_infile.GetListOfKeys():
        tmp_name = k.GetName()
        ## don't reweight any histos that don't match the regexp
        if not re.match(regexp, tmp_name): continue
        ## don't reweight any histos that are already variations of something else
        if 'Up' in tmp_name or 'Down' in tmp_name: continue
        process = tmp_name.split('_')[1] 

        ## now should be left with only the ones we are interested in
        print 'reweighting smooth muscale syst for process', tmp_name
        
        tmp_nominal = tmp_infile.Get(tmp_name)
        tmp_nominal_2d = dressed2D(tmp_nominal,binning, tmp_name+'backrolled')
        n_ptbins = tmp_nominal_2d.GetNbinsY()

        for syst,hist in systsAndHists:
            for shift_dir in ['Up','Down']:
                outname_2d = tmp_nominal_2d.GetName().replace('backrolled','')+'_smooth{lep}scale{idx}{shiftdir}'.format(lep=flav,idx=syst,shiftdir=shift_dir)
                tmp_scaledHisto = copy.deepcopy(tmp_nominal_2d.Clone(outname_2d))
                tmp_scaledHisto.Reset()

                ## loop over all eta bins of the 2d histogram 
                for ieta in range(1,tmp_nominal_2d.GetNbinsX()+1):
                    eta = tmp_nominal_2d.GetXaxis().GetBinCenter(ieta)
                    # try to save some time, most histograms only need eta, not pt (average is made)
                    etabin = max(1, min(binning_histo.GetNbinsX(), binning_histo.GetXaxis().FindFixBin(eta)))
                    tmp_syst_ptscale = 0
                    if syst!='Syst3': 
                        tmp_syst_ptscale = utilities.getRochesterUncertainty(charge,etabin-1,0,hist,True)
                    for ipt in range(1,tmp_nominal_2d.GetNbinsY()+1):                            

                        ptcenter      = tmp_nominal_2d.GetYaxis().GetBinCenter(ipt)
                        ptbin  = max(1, min(binning_histo.GetNbinsY(), binning_histo.GetYaxis().FindFixBin(ptcenter)))
                        # scale_syst is something close to 1: if larger, we move events from left to right
                        # if lower, no event will move from left to right, but actually from right to left
                        syst_ptscale = 1 # it can be seen as 1 + dx with dx > or < 0
                        # the syst can be > or < 1, we define the variation based on the number in the histogram, and arbitrarily call Up the variation based on that.
                        # then we artificially define as Down the variation defined as 1 - (var-1) = 2 -var
                        # generally, the variation would be > 1 (which would mean scaling pt up
                        # but in some regions inside acceptance it could be < 1 ( but we still call it up when creating the alternate histogram)
                        if syst == 'Syst3': 
                            tmp_syst_ptscale = utilities.getRochesterUncertainty(charge,etabin-1,ptbin-1,hist,False)
                        if shift_dir == "Up":
                            syst_ptscale = tmp_syst_ptscale 
                        else:
                            syst_ptscale = 2. - tmp_syst_ptscale
                    
                        nominal_val = tmp_nominal_2d.GetBinContent(ieta,ipt)
                        pt_width    = tmp_nominal_2d.GetYaxis().GetBinWidth(ipt)
                        if syst_ptscale > 1:                            
                            scale_syst = syst_ptscale - 1. 
                            pt      = tmp_nominal_2d.GetYaxis().GetBinUpEdge(ipt)
                            pt_prev = tmp_nominal_2d.GetYaxis().GetBinLowEdge(ipt)
                            if ipt == 1:
                                # take first pt bin from histo with extended pt range
                                nominal_val_prev = nominalHists[tmp_name].GetBinContent(ieta,1)
                                pt_width_prev    = nominalHists[tmp_name].GetYaxis().GetBinWidth(1)   
                            else:
                                nominal_val_prev = tmp_nominal_2d.GetBinContent(ieta,ipt-1)
                                pt_width_prev = tmp_nominal_2d.GetYaxis().GetBinWidth(ipt-1)
                            factor = scale_syst / (1.0 + scale_syst)  # > 0
                            from_prev = factor * (pt_prev / pt_width_prev) * max(0,nominal_val_prev)
                            to_next   = factor * (pt      / pt_width     ) * max(0,nominal_val)
                            tmp_scaledHisto.SetBinContent(ieta,ipt,max(0,nominal_val + from_prev - to_next))
                        else:
                            scale_syst = 1 - syst_ptscale # so scale_syst > 0 again
                            pt      = tmp_nominal_2d.GetYaxis().GetBinLowEdge(ipt)
                            pt_next = tmp_nominal_2d.GetYaxis().GetBinUpEdge(ipt)
                            if ipt == tmp_nominal_2d.GetNbinsY():
                                # histograms with extended pt range have 2 more pt bins
                                nominal_val_next = nominalHists[tmp_name].GetBinContent(ieta,ipt+2)
                                pt_width_next    = nominalHists[tmp_name].GetYaxis().GetBinWidth(ipt+2)
                            else:
                                nominal_val_next = tmp_nominal_2d.GetBinContent(ieta,ipt+1)
                                pt_width_next    = tmp_nominal_2d.GetYaxis().GetBinWidth(ipt+1)
                            factor = scale_syst / (1.0 - scale_syst)
                            from_next = factor * (pt_next / pt_width_next) * max(0,nominal_val_next)
                            to_prev   = factor * (pt      / pt_width     ) * max(0,nominal_val)
                            tmp_scaledHisto.SetBinContent(ieta,ipt,max(0,nominal_val + from_next - to_prev))

                    if "outliers" in tmp_name and shift_dir == "Up":
                        # ad hoc fix for second bin, which gets a very large ratio because the assumption
                        # of flat event density is maximally wrong (would be the same for first bin but let's neglect that)
                        nomiContent = tmp_nominal_2d.GetBinContent(ieta,2)
                        ratioPt3 = 1.0
                        if tmp_nominal_2d.GetBinContent(ieta,3) != 0:
                            ratioPt3 = tmp_scaledHisto.GetBinContent(ieta,3)/tmp_nominal_2d.GetBinContent(ieta,3)
                        tmp_scaledHisto.SetBinContent(ieta,2,ratioPt3*nomiContent)

                ## re-roll the 2D to a 1D histo
                tmp_scaledHisto_1d = unroll2Dto1D(tmp_scaledHisto, newname=tmp_scaledHisto.GetName().replace('2DROLLED',''), cropNegativeBins=False)

                outfile.cd()
                tmp_scaledHisto_1d.Write()

    outfile.Close()
    print "done with the smooth muon scale variations"
    print "Output file: %s" % outfile.GetName()


def addSmoothElectronScaleSyst(infile, regexp, charge, 
                               outdir=None, 
                               uncorrelateByCharge=False,
                               uncorrelateByEtaSide=False,
                               alternateShapeOnly=False):

    indir = os.path.dirname(os.path.abspath(infile))
    outdir = outdir if outdir != None else indir
    flav = 'el'
    isMu = False

    # get eta-pt binning for reco 
    etaPtBinningVec = getDiffXsecBinning(indir+'/binningPtEta.txt', "reco")  # this get two vectors with eta and pt binning
    recoBins = templateBinning(etaPtBinningVec[0],etaPtBinningVec[1])        # this create a class to manage the binnings
    binning = [recoBins.Neta, recoBins.etaBins, recoBins.Npt, recoBins.ptBins]

    tmp_infile = ROOT.TFile(infile, 'read')
    outfile = ROOT.TFile(outdir+'/SmoothScaleSystClosure_{flav}_{ch}.root'.format(flav=flav,ch=charge), 'recreate')

    ##//////////////////////////
    # get nominal histograms with one more pt bin outside acceptance, to do smoothing
    #
    # define new binning for pt
    extendedPtBins = [29.0] + [float(x) for x in recoBins.ptBins] + [57.0]
    binning_recoPt29to57 = [recoBins.Neta, recoBins.etaBins, len(extendedPtBins)-1, extendedPtBins]
    # now open file and get histograms
    # these are TH1 and should be rolled to TH2
    file_nomi_recoPt29to57 = "cards/diffXsec_el_2019_09_22_onlyZandTau_recoPt29to57/nominalShape_{ch}_recoPt29to57.root".format(ch=charge)
    hist_nomi_recoPt29to57 = ROOT.TFile.Open(file_nomi_recoPt29to57)
    nominalHists = {}
    for ikey,e in enumerate(hist_nomi_recoPt29to57.GetListOfKeys()):
        name = e.GetName()
        obj  = e.ReadObj()
        if not obj:
            raise RuntimeError('Unable to read object {n}'.format(n=name))
        nominalHists[name] = dressed2D(obj,binning_recoPt29to57, name+"_recoPt29to57_backrolled")
        nominalHists[name].SetDirectory(0)
    hist_nomi_recoPt29to57.Close()
    # these will be used later
    #
    ##//////////////////////////

    etabins_el = array('f',[0.0, 1.0, 1.5, 2.1, 2.4])
    scaleSyst_el = ROOT.TH1F('scaleSyst_el','',len(etabins_el)-1,etabins_el)
    scaleSyst_el.SetBinContent(1,0.003)
    scaleSyst_el.SetBinContent(2,0.005)
    scaleSyst_el.SetBinContent(3,0.008)
    scaleSyst_el.SetBinContent(4,0.01)
        
    scaleSyst = utilities.getExclusiveBinnedSyst(scaleSyst_el)
    etabins = etabins_el

    ###////////
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
                                    nominal_val = tmp_nominal_2d.GetBinContent(ieta,ipt)
                                    pt_width    = tmp_nominal_2d.GetYaxis().GetBinWidth(ipt)
                                    if shift_dir == "Up":                            
                                        pt      = tmp_nominal_2d.GetYaxis().GetBinUpEdge(ipt)
                                        pt_prev = tmp_nominal_2d.GetYaxis().GetBinLowEdge(ipt)
                                        if ipt == 1:
                                            # take first pt bin from histo with extended pt range
                                            nominal_val_prev = nominalHists[tmp_name.replace("_el_","_lep_")].GetBinContent(ieta,1)
                                            pt_width_prev    = nominalHists[tmp_name.replace("_el_","_lep_")].GetYaxis().GetBinWidth(1)   
                                        else:
                                            nominal_val_prev = tmp_nominal_2d.GetBinContent(ieta,ipt-1)
                                            pt_width_prev = tmp_nominal_2d.GetYaxis().GetBinWidth(ipt-1)
                                        factor = scale_syst / (1.0 + scale_syst)  # > 0
                                        from_prev = factor * (pt_prev / pt_width_prev) * max(0,nominal_val_prev)
                                        to_next   = factor * (pt      / pt_width     ) * max(0,nominal_val)
                                        tmp_scaledHisto.SetBinContent(ieta,ipt,max(0,nominal_val + from_prev - to_next))
                                    else:
                                        pt      = tmp_nominal_2d.GetYaxis().GetBinLowEdge(ipt)
                                        pt_next = tmp_nominal_2d.GetYaxis().GetBinUpEdge(ipt)
                                        if ipt == tmp_nominal_2d.GetNbinsY():
                                            # need last bin, and extended histo has two more bins than nominal
                                            nominal_val_next = nominalHists[tmp_name.replace("_el_","_lep_")].GetBinContent(ieta,ipt+2)
                                            pt_width_next    = nominalHists[tmp_name.replace("_el_","_lep_")].GetYaxis().GetBinWidth(ipt+2)
                                        else:
                                            # histograms with extended pt range have 2 more pt bins
                                            nominal_val_next = tmp_nominal_2d.GetBinContent(ieta,ipt+1)
                                            pt_width_next    = tmp_nominal_2d.GetYaxis().GetBinWidth(ipt+1)
                                        factor = scale_syst / (1.0 - scale_syst)
                                        from_next = factor * (pt_next / pt_width_next) * max(0,nominal_val_next)
                                        to_prev   = factor * (pt      / pt_width     ) * max(0,nominal_val)
                                        tmp_scaledHisto.SetBinContent(ieta,ipt,max(0,nominal_val + from_next - to_prev))

                        if "outliers" in tmp_name and shift_dir == "Up":
                            # ad hoc fix for second bin, which gets a very large ratio because the assumption
                            # of flat event density is maximally wrong (would be the same for first bin but let's neglect that)
                            nomiContent = tmp_nominal_2d.GetBinContent(ieta,2)
                            ratioPt3 = 1.0
                            if tmp_nominal_2d.GetBinContent(ieta,3) != 0:
                                ratioPt3 = tmp_scaledHisto.GetBinContent(ieta,3)/tmp_nominal_2d.GetBinContent(ieta,3)
                            tmp_scaledHisto.SetBinContent(ieta,2,ratioPt3*nomiContent)


                    ## re-roll the 2D to a 1D histo
                    tmp_scaledHisto_1d = unroll2Dto1D(tmp_scaledHisto, newname=tmp_scaledHisto.GetName().replace('2DROLLED',''), cropNegativeBins=False)

                    outfile.cd()
                    tmp_scaledHisto_1d.Write()

    outfile.Close()
    print "done with the smooth electron scale variations"
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
    parser.add_option("--function",   dest="function", type="string", default="addSmoothLeptonScaleSyst", help="Specify function to use");
    parser.add_option(       '--crop-negative-bin', dest='cropNegativeBin' , default=False, action='store_true', help='Crop negative bins in nominal and alternate')
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

    if options.function == "addSmoothLeptonScaleSyst":
        addSmoothLeptonScaleSyst(options.infile, options.processes, charge, isMu,
                                 outdir=outdir,
                                 alternateShapeOnly=False,
                                 cropNegativeBin=options.cropNegativeBin)        
    elif options.function == "addSmoothLeptonScaleSyst_algoHelicity":
        print "Using addSmoothLeptonScaleSyst_algoHelicity()"
        addSmoothLeptonScaleSyst_algoHelicity(options.infile, options.processes, charge, isMu,
                                              outdir=outdir)
    else:
        addSmoothLepScaleSyst(options.infile, options.processes, charge, isMu,
                              outdir=outdir,
                              uncorrelateByCharge=options.uncorrelatePtScaleByCharge,
                              uncorrelateByEtaSide=options.uncorrelatePtScaleByEtaSide,
                              getExtremePtBinsFromOldScale=True)
