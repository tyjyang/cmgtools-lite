import ROOT, os
from array import array

ROOT.gROOT.SetBatch()

#infile = ROOT.TFile('~mdunser/www/private/w-helicity-13TeV/templates/2018-09-29-unrolled/templates_2D_plus.root' ,'read')
infile = ROOT.TFile('/afs/cern.ch/user/m/mdunser/www/private/w-helicity-13TeV/templates/2019-09-30/templates_2D_plus.root' ,'read')

hists = []
names = []

h1 = infile.Get('Wplus_right_mu_Ybin_0') ; hists.append(h1); names.append('W^{+}_{right}: 0.0 < |y_{W}| < 0.25')
h2 = infile.Get('Wplus_right_mu_Ybin_2') ; hists.append(h2); names.append('W^{+}_{right}: 0.5 < |y_{W}| < 0.75')
h3 = infile.Get('Wplus_left_mu_Ybin_8')  ; hists.append(h3); names.append('W^{+}_{left}: 2.0 < |y_{W}| < 2.25')
#h4 = infile.Get('Wplus_left_Wplus_left_mu_Ybin_11') ; hists.append(h4)

ROOT.TColor.CreateGradientColorTable(3,
                                  array ("d", [0.00, 0.50, 1.00]),
                                  ##array ("d", [1.00, 1.00, 0.00]),
                                  ##array ("d", [0.70, 1.00, 0.34]),
                                  ##array ("d", [0.00, 1.00, 0.82]),
                                  array ("d", [1.00, 0.00, 1.00]),
                                  array ("d", [1.00, 0.34, 0.65]),
                                  array ("d", [1.00, 0.82, 0.00]),
                                  255,  0.95)

#ROOT.gStyle.SetPalette(56)


lat = ROOT.TLatex(); lat.SetNDC()
lat.SetTextFont(42)
lat.SetTextSize(0.05)

cutoff = 0.0015

colors = [ROOT.kAzure-7, ROOT.kRed+1, ROOT.kGreen+1]
styles = [3001, 3002, 3002, 3008]

ROOT.gStyle.SetHatchesSpacing(0.01)

c1 = ROOT.TCanvas('asdf','',1200,900)
c1.cd()

marginR = 0.02
marginL = 0.12
c1.SetRightMargin(marginR)
c1.SetLeftMargin(marginL)
c1.SetBottomMargin(0.21)
c1.SetTopMargin(0.10)
    
leg = ROOT.TLegend(marginL, 0.02, 1.-marginR, 0.08)
leg.SetNColumns(3)
leg.SetColumnSeparation(0.03)
#leg.SetFillStyle(0)
leg.SetFillColor(0)
leg.SetBorderSize(0)

for ih,hist in enumerate(hists):
    hist.SetTitle('')
    hist.Scale( 1./hist.Integral())
    
    for ix in range(1,hist.GetXaxis().GetNbins()+1):
        for iy in range(1,hist.GetYaxis().GetNbins()+1):
            if hist.GetBinContent(ix,iy) < cutoff:
                hist.SetBinContent(ix, iy, 0)
                
    hist.SetFillStyle(styles[ih])
    hist.SetFillColor(colors[ih])
    hist.SetLineColor(colors[ih])
    
    hist.GetYaxis().SetTitleFont(42)
    hist.GetYaxis().SetTitleSize(0.05)
    hist.GetYaxis().SetLabelSize(0.05)
    
    hist.GetXaxis().SetTitleFont(42)
    hist.GetXaxis().SetTitleSize(0.05)
    hist.GetXaxis().SetLabelSize(0.05)

    leg.AddEntry(hist, names[ih], 'f')
    
    hist.Draw('box same' if ih else 'box')
    
lat.DrawLatex(marginL, 0.92, '#bf{CMS} #it{Simulation}')
#lat.DrawLatex(0.73, 0.92, '36 fb^{-1} (13 TeV)')

leg.Draw('same')
    
c1.SaveAs('~mdunser/www/private/w-helicity-13TeV/paperPlots/template_example.pdf')
c1.SaveAs('~mdunser/www/private/w-helicity-13TeV/paperPlots/template_example.png')
