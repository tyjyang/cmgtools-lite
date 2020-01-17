import ROOT, array, copy
import rollingFunctions as rf
from utilities import templateBinning, getBinning


def checkIfBinInRange(eta, pt, reta, rpt):
    inRange = -1
    for ir,r in enumerate(reta):
        if r[0] <= eta <= r[1] and rpt[ir][0] <= pt <= rpt[ir][1]:
            inRange = ir
    return inRange

def makeGenericMultipliers(binningFile, systName, pts, etas, sizes):


    if   type(pts) == str:
        ranges_pt  = [ [float(i.split(',')[0]),float(i.split(',')[1])] for i in pts .split(';') ]
        ranges_eta = [ [float(i.split(',')[0]),float(i.split(',')[1])] for i in etas.split(';') ]
        syst_sizes = [ float(i) for i in sizes.split(';') ]
    elif type(pts) == list:
        ranges_pt  = pts
        ranges_eta = etas
        syst_sizes = sizes
    else:
        print 'you need to give options to the generic systematic producer in either string or list format'
        print 'exititng...'
        exit(0)
    
    if not ( len(ranges_pt) == len(syst_sizes) == len(ranges_eta) ):
        print 'the option sizes has to be the same length of at least one of ranges_pt and ranges_eta'
        print 'exititng...'
        exit(0)

    # get eta-pt binning for both reco 
    etaPtBinningVec = getBinning(binningFile, "reco")  # this get two vectors with eta and pt binning
    recoBins = templateBinning(etaPtBinningVec[0],etaPtBinningVec[1])        # this create a class to manage the binnings
    binning = [recoBins.Neta, recoBins.etaBins, recoBins.Npt, recoBins.ptBins]                                                            

    h_syst_2d_up = ROOT.TH2D(systName+'Up', systName+'Up', 
                          recoBins.Neta, array.array('d', recoBins.etaBins), 
                          recoBins.Npt , array.array('d', recoBins.ptBins ))

    h_syst_2d_dn = h_syst_2d_up.Clone(h_syst_2d_up.GetName().replace('Up','Down'))
    h_syst_2d_dn.SetTitle(h_syst_2d_up.GetTitle().replace('Up','Down'))
    
    
    for ix in range(1,h_syst_2d_up.GetNbinsX()+1):
        for iy in range(1,h_syst_2d_up.GetNbinsY()+1):
            tmp_eta = h_syst_2d_up.GetXaxis().GetBinCenter(ix)
            tmp_pt  = h_syst_2d_up.GetYaxis().GetBinCenter(iy)
            isBinInRange = checkIfBinInRange(tmp_eta, tmp_pt, ranges_eta, ranges_pt)
            if isBinInRange > -1:
                h_syst_2d_up.SetBinContent(ix, iy, 1.+(syst_sizes[isBinInRange]-1.))
                h_syst_2d_up.SetBinError  (ix, iy, 0.)
                h_syst_2d_dn.SetBinContent(ix, iy, 1.-(syst_sizes[isBinInRange]-1.))
                h_syst_2d_dn.SetBinError  (ix, iy, 0.)
            else:
                h_syst_2d_up.SetBinContent(ix, iy, 1.)
                h_syst_2d_up.SetBinError  (ix, iy, 0.)
                h_syst_2d_dn.SetBinContent(ix, iy, 1.)
                h_syst_2d_dn.SetBinError  (ix, iy, 0.)

    h_syst_2d_up.GetZaxis().SetRangeUser(0.,2.)
    h_syst_2d_dn.GetZaxis().SetRangeUser(0.,2.)

    h_syst_1d_up = rf.unroll2Dto1D(h_syst_2d_up, h_syst_2d_up.GetName())
    h_syst_1d_dn = rf.unroll2Dto1D(h_syst_2d_dn, h_syst_2d_dn.GetName())

    return h_syst_2d_up, h_syst_2d_dn, h_syst_1d_up, h_syst_1d_dn
    
def makeEffStatMultipliers(binningFile, indir):

    allHistos = []
    
    dummy_central = makeGenericMultipliers(binningFile, 'dummy', [[0.,100.]], [[-2.4,2.4]], [1.])[0]
    
    for ch in ['plus', 'minus']:
        #print 'at charge', ch
        infile = ROOT.TFile(indir+'/systEff_trgmu_{ch}_mu.root'.format(ch=ch), 'READ')
        for par in ['p0', 'p1', 'p2']:
            #print 'at parameter', par
            npar = int(par[-1])
            tmp_varhisto_orig = infile.Get(par)

            tmp_varhisto_clone = tmp_varhisto_orig.Clone('Efficiencies_{ch}_{p}_original'.format(ch=ch,p=par))
            tmp_varhisto_clone.SetTitle(tmp_varhisto_clone.GetName())
            
            allHistos.append(copy.deepcopy(tmp_varhisto_clone))

            tmp_varhisto_goodbin = dummy_central.Clone(tmp_varhisto_clone.GetName().replace('_original',''))
            tmp_varhisto_goodbin.SetTitle(tmp_varhisto_goodbin.GetName())

            for ieta in range(1,tmp_varhisto_goodbin.GetNbinsX()+1):
                #print 'at eta', ieta
                eta = tmp_varhisto_goodbin.GetXaxis().GetBinCenter(ieta)

                outname_2d = 'EffStat{p}mu{ch}{eta}'.format(p=npar,eta=ieta,ch=ch)
            
                tmp_varhisto_up = copy.deepcopy(tmp_varhisto_goodbin.Clone(outname_2d+'Up'))  ; tmp_varhisto_up.SetTitle(tmp_varhisto_up.GetName())
                tmp_varhisto_dn = copy.deepcopy(tmp_varhisto_goodbin.Clone(outname_2d+'Down')); tmp_varhisto_dn.SetTitle(tmp_varhisto_dn.GetName())
                #rangemax = tmp_varhisto_orig.GetZaxis().GetXmax()
                tmp_varhisto_up.GetZaxis().SetRangeUser(0.99, 1.01)
                tmp_varhisto_dn.GetZaxis().SetRangeUser(0.99, 1.01)
                
                ## loop over all pT bins in that bin of eta (which is ieta)
                for ipt in range(1,tmp_varhisto_up.GetNbinsY()+1):
                    pt = tmp_varhisto_up.GetYaxis().GetBinCenter(ipt)
                    tmp_bincontent = tmp_varhisto_clone.GetBinContent(tmp_varhisto_clone.FindBin(eta, pt))
                    tmp_varhisto_up.SetBinContent(ieta, ipt, 1.+tmp_bincontent); tmp_varhisto_up.SetBinError  (ieta, ipt, 0.)
                    tmp_varhisto_dn.SetBinContent(ieta, ipt, 1.-tmp_bincontent); tmp_varhisto_dn.SetBinError  (ieta, ipt, 0.)

                ## re-roll the 2D to a 1D histo
                tmp_varhisto_up_1d = rf.unroll2Dto1D(tmp_varhisto_up, newname=tmp_varhisto_up.GetName())
                tmp_varhisto_dn_1d = rf.unroll2Dto1D(tmp_varhisto_dn, newname=tmp_varhisto_dn.GetName())

                allHistos += [tmp_varhisto_up, tmp_varhisto_dn]
                allHistos += [copy.deepcopy(tmp_varhisto_up_1d), copy.deepcopy(tmp_varhisto_dn_1d)]

    return allHistos
    
##    outfile = ROOT.TFile(outfile_name, 'RECREATE')
##    for h in allHistos:
##        h.Write()
##    outfile.Close()
        

