#!/usr/bin/env python                                                       
 
#from shutil import copyfile
import re, sys, os, os.path, ROOT, copy, math
from array import array
#import numpy as np

from w_helicity_13TeV.make_diff_xsec_cards import getXYBinsFromGlobalBin
from w_helicity_13TeV.make_diff_xsec_cards import getArrayParsingString
from w_helicity_13TeV.make_diff_xsec_cards import getDiffXsecBinning
from w_helicity_13TeV.make_diff_xsec_cards import templateBinning
from w_helicity_13TeV.make_diff_xsec_cards import get_ieta_ipt_from_process_name
from w_helicity_13TeV.make_diff_xsec_cards import get_ieta_from_process_name

from w_helicity_13TeV.rollingFunctions import roll1Dto2D, dressed2D, unroll2Dto1D

from w_helicity_13TeV.mergeCardComponentsAbsY import mirrorShape
from w_helicity_13TeV.mergeCardComponentsAbsY import putUncorrelatedFakes
#from w_helicity_13TeV.mergeCardComponentsAbsY import putEffStatHistos   # use function below that can manage case with wider template bins along eta

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
        print 'reweighting testeffsyst nuisances for process', tmp_name
        
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
            effsystvals[1.479] = 0.01    # 0.013
            effsystvals[2.0] =   0.0093  # 0.016
            # then we will sum L1 prefiring term            

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
                        if isMu: 
                            scaling += effsystvals[key]                        
                        else:
                            ptL1  = tmp_scaledHisto_up.GetYaxis().GetBinLowEdge(ipt)
                            etaL1 =  tmp_scaledHisto_up.GetXaxis().GetBinUpEdge(ieta) if (etaval < 0) else tmp_scaledHisto_up.GetXaxis().GetBinLowEdge(ieta)
                            if (abs(etaL1) >= 1.479 and ptL1 >= 35.0):                                
                                ptval = tmp_scaledHisto_up.GetYaxis().GetBinCenter(ipt)
                                L1RelativeUncertainty = hL1.GetBinError(hL1.GetXaxis().FindFixBin(etaval),hL1.GetYaxis().FindFixBin(ptval))
                                scaling += math.sqrt(effsystvals[key]*effsystvals[key] + L1RelativeUncertainty*L1RelativeUncertainty)

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
            for ietaErf in range(1,nEtaErfPar+1):

                # FIXME, need more care for electron channel, template binning can be odd around the gap
                # identify template eta which contain that erfPar eta

                etabinOffset = 0
                # the parhist histogram for muons is defined with 50 bins, where the two with |eta| in [2.4,2.5] are filled like the inner ones
                # these bins for muons do not make sense: when looping on the bin number, we need to add an offset to skip the first bin
                if isMu and parhist.GetNbinsX() == 50:
                    etabinOffset = 1

                ietaTemplate = tmp_nominal_2d.GetXaxis().FindFixBin( parhist.GetXaxis().GetBinCenter(ietaErf+etabinOffset) )
                #ietaTemplate = tmp_nominal_2d.GetXaxis().FindFixBin( parhist.GetXaxis().GetBinCenter(ietaErf + parhistBinOffset) )

                # if template has eta range narrower than parhist, continue
                if ietaTemplate == 0 or ietaTemplate == (1 + tmp_nominal_2d.GetNbinsX()):
                    continue

                if not isMu:
                    if parhist.GetXaxis().GetBinCenter(ietaErf) > 1.4 and parhist.GetXaxis().GetBinCenter(ietaErf) < 1.5:
                        ietaTemplate = tmp_nominal_2d.GetXaxis().FindFixBin( 1.41 )

                # for electrons there is the gap:
                # eg, we can have 1.3, 1.4442, 1.5, 1.566, 1.7
                # and the bin centers of the erfPar are 1.35, 1.45, 1.55, ...
                # so with 1.35 we get the template bin from 1.3 to 1.4442
                # but with 1.45 we would fall in the empty bin, while we would want to reweight the template bin that would span between 1.4 to 1.5
                # in this case, the syst is ineffective. The same for 1.55, which select the other "side" of the empty bin around 1.5
                # one could device a fancy way to isolate the gap and assign a larger uncertainty around it but let's keep it simple

                #print "ErfPar{p}EffStat{ieta}".format(p=npar,ieta=ietaErf)
                outname_2d = tmp_nominal_2d.GetName().replace('backrolled','')+'_ErfPar{p}EffStat{ieta}{fl}{ch}2DROLLED'.format(p=npar,ieta=ietaErf,fl=flavour,ch=charge)
            
                tmp_scaledHisto_up = copy.deepcopy(tmp_nominal_2d.Clone(outname_2d+'Up'))
                tmp_scaledHisto_dn = copy.deepcopy(tmp_nominal_2d.Clone(outname_2d+'Down'))
                
                ## loop over all pT bins in that bin of eta (which is ieta)
                for ipt in range(1,tmp_scaledHisto_up.GetNbinsY()+1):
                    tmp_bincontent = tmp_scaledHisto_up.GetBinContent(ietaTemplate, ipt)
                    ybincenter = tmp_scaledHisto_up.GetYaxis().GetBinCenter(ipt)
                    ## now get the content of the parameter variation!
                    phistybin = parhist.GetYaxis().FindBin(ybincenter)
                    tmp_scale = parhist.GetBinContent(ietaErf, phistybin)
                    scaling = math.sqrt(2.)*tmp_scale
                    ## scale electrons by sqrt(2) due to the input file being charge inclusive
                    if flavour == 'el':
                        scaling *= math.sqrt(2)
                    # modify scaling if template bin has larger width than the ErfPar histogram
                    # assuming the two uncertainties have the same order of magnitude, we should scale down by sqrt(binWidth/0.1)
                    # however, we could just be conservative and not scale anything
                    #if binwidths[ietaTemplate-1] > 1.: scaling = scaling / binwidths[ietaTemplate-1]
                    #if binwidths[ietaTemplate-1] > 1.: scaling = scaling / math.sqrt( binwidths[ietaTemplate-1] )
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
def putBinUncEffStatHistosDiffXsec(infile,regexp,charge, outdir=None, isMu=True, suffix="", subtractSystPart = False):

    # subtractSystPart is needed to test only the real statistical component (stat+syst is saved in the inputs)
    # this is very hardcoded

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
    else:
        parfile_name = basedir+'/leptonSF/new2016_madeSummer2018/systEff_trgel.root'

    parfile = ROOT.TFile(parfile_name, 'read')
    

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
            parhist = parfile.Get('scaleFactor')
            ## loop over all eta bins of the 2d histogram
            for ietaErf in range(1,nEtaErfPar+1):

                # FIXME, need more care for electron channel, template binning can be odd around the gap
                # identify template eta which contain that erfPar eta

                etabinOffset = 0
                # the parhist histogram for muons is defined with 50 bins, where the two with |eta| in [2.4,2.5] are filled like the inner ones
                # these bins for muons do not make sense: when looping on the bin number, we need to add an offset to skip the first bin
                if isMu and parhist.GetNbinsX() == 50:
                    etabinOffset = 1

                ietaTemplate = tmp_nominal_2d.GetXaxis().FindFixBin( parhist.GetXaxis().GetBinCenter(ietaErf+etabinOffset) )
                #ietaTemplate = tmp_nominal_2d.GetXaxis().FindFixBin( parhist.GetXaxis().GetBinCenter(ietaErf + parhistBinOffset) )

                # if template has eta range narrower than parhist, continue
                if ietaTemplate == 0 or ietaTemplate == (1 + tmp_nominal_2d.GetNbinsX()):
                    continue

                if not isMu:
                    if parhist.GetXaxis().GetBinCenter(ietaErf) > 1.4 and parhist.GetXaxis().GetBinCenter(ietaErf) < 1.5:
                        ietaTemplate = tmp_nominal_2d.GetXaxis().FindFixBin( 1.41 )

                # for electrons there is the gap:
                # eg, we can have 1.3, 1.4442, 1.5, 1.566, 1.7
                # and the bin centers of the erfPar are 1.35, 1.45, 1.55, ...
                # so with 1.35 we get the template bin from 1.3 to 1.4442
                # but with 1.45 we would fall in the empty bin, while we would want to reweight the template bin that would span between 1.4 to 1.5
                # in this case, the syst is ineffective. The same for 1.55, which select the other "side" of the empty bin around 1.5
                # one could device a fancy way to isolate the gap and assign a larger uncertainty around it but let's keep it simple

                #print "ErfPar{p}EffStat{ieta}".format(p=npar,ieta=ietaErf)
                outname_2d = tmp_nominal_2d.GetName().replace('backrolled','')+'_BinUncEffStat{ieta}{fl}{ch}2DROLLED'.format(p=npar,ieta=ietaErf,fl=flavour,ch=charge)
            
                tmp_scaledHisto_up = copy.deepcopy(tmp_nominal_2d.Clone(outname_2d+'Up'))
                tmp_scaledHisto_dn = copy.deepcopy(tmp_nominal_2d.Clone(outname_2d+'Down'))
                
                maxyparhist = parhist.GetYaxis().GetBinUpEdge(parhist.GetNbinsY())

                ## loop over all pT bins in that bin of eta (which is ieta)
                for ipt in range(1,tmp_scaledHisto_up.GetNbinsY()+1):
                    tmp_bincontent = tmp_scaledHisto_up.GetBinContent(ietaTemplate, ipt)
                    ybincenter = tmp_scaledHisto_up.GetYaxis().GetBinCenter(ipt)
                    ## now get the content of the parameter variation!                    
                    phistybin = parhist.GetYaxis().FindBin(ybincenter)
                    if maxyparhist <= phistybin:
                        phistybin = parhist.GetNbinsY()
                    # if I need to reduce the variation as below, first get the variation (with proper sign), scale it, and then sum 1
                    tmp_scale_up =        parhist.GetBinError(ietaErf, phistybin) / parhist.GetBinContent(ietaErf, phistybin)
                    tmp_scale_dn = -1.0 * parhist.GetBinError(ietaErf, phistybin) / parhist.GetBinContent(ietaErf, phistybin)
                    # modify scaling if template bin has larger width than the ErfPar histogram
                    # no longer scaling, trying to be conservative
                    # if binwidths[ietaTemplate-1] > 1.: 
                    #     tmp_scale_up = tmp_scale_up / binwidths[ietaTemplate-1]
                    #     tmp_scale_dn = tmp_scale_dn / binwidths[ietaTemplate-1]
                    if sub:
                        if (abseta<1): syst = 0.002
        elif (abseta<1.5):   syst = 0.004;                                                                                        
                else                   syst = 0.014;                                                                                                                                     

                    tmp_scale_up = 1 +  tmp_scale_up
                    tmp_scale_dn = 1 +  tmp_scale_dn
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
        if "x_TauDecaysW_wtau_" in name: newname = name.replace("x_TauDecaysW_wtau_","x_TauDecaysW_")
        else: newname = name
        if name == "x_TauDecaysW": taunominal = obj.Clone(name)
        if "_pdf" in name:
            taupdf.append(obj.Clone(newname))
            continue
        if newname.endswith("Dn"): newname = newname.replace("Dn","Down")
        newobj = obj.Clone(newname)
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

fileTaueffStat = "{od}ErfParEffStat_{fl}_{ch}_Tau.root".format(od=outdir, fl=flavour, ch=options.charge)  
# this name is used inside putEffStatHistosDiffXsec (do not change it outside here)
print "Now adding ErfParXEffStatYY systematics to x_TauDecaysW process"
print "Will create file --> {of}".format(of=fileTaueffStat)
if not options.dryrun: putEffStatHistosDiffXsec(Taufile, 'x_TauDecaysW', charge, outdir, isMu=True if flavour=="mu" else False, suffix="_Tau")
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

print "Now merging signal + Z + data + other backgrounds + Fakes*Uncorrelated + ZEffStat"
sigfile = "{osig}W{fl}_{ch}_shapes_signal.root".format(osig=options.indirSig, fl=flavour, ch=charge)

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
    print "Now adding TestEffSystYY systematics to x_Z and x_W.* process"
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
    if not options.dryrun: os.system(cmdMergeWithTestEffSyst)
    print ""
    print "Wrote again root file in %s" % shapename
    print ""

print ""
print "-"*20
print ""
