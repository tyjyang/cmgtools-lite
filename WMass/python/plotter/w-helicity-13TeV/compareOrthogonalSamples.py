import ROOT, copy, math
from rollingFunctions import roll1Dto2D, dressed2D, unroll2Dto1D
from make_diff_xsec_cards import getDiffXsecBinning
from make_diff_xsec_cards import templateBinning

ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetPadRightMargin(0.15)
ROOT.gROOT.SetBatch()

etaPtBinningVec = getDiffXsecBinning('../cards/helicity_2019_06_04_orthogonalHelicitySamples/binningPtEta.txt', "reco")  # this get two vectors with eta and pt binning
recoBins = templateBinning(etaPtBinningVec[0],etaPtBinningVec[1])        # this create a class to manage the binnings
binning = [recoBins.Neta, recoBins.etaBins, recoBins.Npt, recoBins.ptBins]



for iy,yhelbin in enumerate(['Wminus_left_mu_Ybin_0', 'Wplus_right_mu_Ybin_5', 'Wminus_right_mu_Ybin_8', 'Wplus_long_mu_Ybin_2']):

    canv = ROOT.TCanvas('foo'+yhelbin, '', 800, 800)
    canv.Divide(2,2)


    nom_f = ROOT.TFile('../cards/helicity_2019_05_15_rebinnigEta/part0/'+yhelbin+'.input.root')
    ort_f = ROOT.TFile('../cards/helicity_2019_06_04_orthogonalHelicitySamples/part0/'+yhelbin+'.input.root')
    
    nom_h = nom_f.Get('x_'+'_'.join(yhelbin.split('_')[:2])); nom_2d = dressed2D(nom_h, binning, yhelbin+'_2d')
    ort_h = ort_f.Get('x_'+'_'.join(yhelbin.split('_')[:2])); ort_2d = dressed2D(ort_h, binning, yhelbin+'_2d')

    ratio = nom_2d.Clone('ratio')

    ratio.Divide(ort_2d)
    
    ratio.GetZaxis().SetRangeUser(0.9,1.1)

    ratio.SetTitle('ratio for '+ ' '.join(yhelbin.split('_')).replace('minus','-').replace('plus','+') )

    pulls = ROOT.TH1F('pulls'+yhelbin,'pulls', 20, -2., 2.0)
    for ib in range(1,nom_h.GetNbinsX()+1):
        if math.sqrt(nom_h.GetBinError(ib)**2 + ort_h.GetBinError(ib)**2 ):
            pulls.Fill( (nom_h.GetBinContent(ib) - ort_h.GetBinContent(ib) ) / math.sqrt(nom_h.GetBinError(ib)**2 + ort_h.GetBinError(ib)**2 ) )

    
    canv.cd(1)
    nom_2d.Draw('colz')
    canv.cd(2)
    ort_2d.Draw('colz')
    canv.cd(3)
    ratio.Draw('colz')
    canv.cd(4)
    pulls.Draw('hist')

    canv.SaveAs('~/www/private/w-helicity-13TeV/orthogonalSamples/'+yhelbin+'.pdf')
    canv.SaveAs('~/www/private/w-helicity-13TeV/orthogonalSamples/'+yhelbin+'.png')
