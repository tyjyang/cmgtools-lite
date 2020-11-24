## source /cvmfs/sft.cern.ch/lcg/app/releases/ROOT/6.22.02/x86_64-centos7-gcc48-opt/bin/thisroot.sh
import ROOT
from array import array

ROOT.gROOT.SetBatch()

def formatHisto(hist):                                                 
    hist.GetXaxis().SetTitleOffset(1.02)
    hist.GetXaxis().SetTitleSize(0.06)
    hist.GetXaxis().SetLabelSize(0.06)

    hist.GetYaxis().SetTitleOffset(1.02)
    hist.GetYaxis().SetTitleSize(0.06)
    hist.GetYaxis().SetLabelSize(0.06)

    hist.GetZaxis().SetLabelSize(0.06)

ROOT.gStyle.SetOptStat(0)

ROOT.ROOT.EnableImplicitMT()

df_read = ROOT.RDataFrame("Events", "root://eoscms.cern.ch//store/cmst3/group/wmass/w-mass-13TeV/NanoAOD/WplusJetsToMuNu_TuneCP5_13TeV-powhegMiNNLO-pythia8-photos/NanoAODv7/201022_231658/0000/SMP-RunIISummer16NanoAODv7-00336_*.root")#SMP-RunIISummer16NanoAODv7-00336_345.root")

print('done loading the RDF')


## load the macro with some functions for the fractions etc.
ROOT.gInterpreter.ProcessLine(".L csFunctions.cc+")


#tmpdf = df_read.Define("preFSRLeps", "preFSRLeptons(nGenPart, GenPart_pdgId, GenPart_genPartIdxMother, GenPart_status, GenPart_statusFlags, GenPart_pt, GenPart_eta, GenPart_phi, GenPart_mass)")
tmpdf = df_read.Define("lepindices", "getIndexLepI(1,nGenPart, GenPart_pdgId, GenPart_genPartIdxMother, GenPart_status, GenPart_statusFlags)")
tmpdf = tmpdf  .Define("lindex1", "lepindices.first")
tmpdf = tmpdf  .Define("lindex2", "lepindices.second")
#tmpdf = tmpdf  .Define("lindex1", "getIndexLepI(1,nGenPart, GenPart_pdgId, GenPart_genPartIdxMother, GenPart_status, GenPart_statusFlags)")
#tmpdf = tmpdf  .Define("lindex2", "getIndexLepI(2,nGenPart, GenPart_pdgId, GenPart_genPartIdxMother, GenPart_status, GenPart_statusFlags)")

tmpdf = tmpdf.Define("pmu1", "ROOT::Math::PtEtaPhiMVector(GenPart_pt[lindex1], GenPart_eta[lindex1], GenPart_phi[lindex1], GenPart_mass[lindex1])")
tmpdf = tmpdf.Define("pmu2", "ROOT::Math::PtEtaPhiMVector(GenPart_pt[lindex2], GenPart_eta[lindex2], GenPart_phi[lindex2], GenPart_mass[lindex2])")
#tmpdf = tmpdf.Define("pmu2", "ROOT::Math::PtEtaPhiMVector(LHEPart_pt[1], LHEPart_eta[1], LHEPart_phi[1], LHEPart_mass[1])")
newdf = tmpdf.Define("pV", "pmu1+pmu2")
dfv = newdf.Define("ptV", "pV.Pt()")
dfv = dfv.Define("yV", "pV.Rapidity()")
dfv = dfv.Define("absyV", "fabs(pV.Rapidity())")
dfv = dfv.Define("phiV", "pV.Phi()")
dfv = dfv.Define("mV", "pV.M()")
#dfv = dfv.Filter("mV < 101 && mZ > 81")


tmpdf2 = dfv.Define("pmup", "GenPart_pdgId[lindex1] < 0 ? pmu1 : pmu2")
tmpdf2 = tmpdf2.Define("pmum", "GenPart_pdgId[lindex1] > 0 ? pmu1 : pmu2")
tmpdf2 = tmpdf2.Define("pmumpt", "pmum.Pt()")

dfcs = tmpdf2.Define("costcs", "cosThetaCS(pmup, pmum)")
dfcs = dfcs.Define("thetacs", "TMath::ACos(costcs)")
dfcs = dfcs.Define("phics", "phiCS(pmup, pmum)")
dfcs = dfcs.Define("genWeightReg", "fabs(genWeight) > 4. ? 4.*genWeight/fabs(genWeight) : genWeight")

yqpt = array('d', [0.0, 2., 3., 4., 4.75, 5.5, 6.5, 8., 9., 10., 12., 14., 16., 18., 20., 23., 27., 32., 40., 55., 100.0])
nbins_pt = len(yqpt)-1
bins_pt = array('d', yqpt)
#bins_pt = [0.0, 2., 3., 4., 4.75, 5.5, 6.5, 8., 9., 10., 12., 14., 16., 18., 20., 23., 27., 32., 40., 55., 100.0]

yqyv = [0.0+i*0.25 for i in range(17)] + [10.] ## keep it the same as the fraction derivation
nbins_yv = len(yqyv)-1
bins_yv = array('d', yqyv)
#bins_yv = [0.0+i*0.25 for i in range(17)] + [10.]

dfcs = dfcs.Define('term_const' , '(1.+(costcs*costcs))')
dfcs = dfcs.Define('term_a0'    , '(1.-3.*(costcs*costcs))')
dfcs = dfcs.Define('term_a1'    , 'TMath::Sin(2.*TMath::ACos(costcs)) * TMath::Cos(phics)')
dfcs = dfcs.Define('term_a2'    , 'TMath::Sin(TMath::ACos(costcs))*TMath::Sin(TMath::ACos(costcs)) * TMath::Cos(2.*phics)')
dfcs = dfcs.Define('term_a3'    , 'TMath::Sin(TMath::ACos(costcs))*TMath::Cos(phics)')
dfcs = dfcs.Define('term_a4'    , 'costcs')
dfcs = dfcs.Define('term_a5'    , 'TMath::Sin(TMath::ACos(costcs))*TMath::Sin(TMath::ACos(costcs)) * TMath::Sin(2.*phics)')
dfcs = dfcs.Define('term_a6'    , 'TMath::Sin(2.*TMath::ACos(costcs)) * TMath::Sin(phics)')
dfcs = dfcs.Define('term_a7'    , 'TMath::Sin(TMath::ACos(costcs)) * TMath::Sin(phics)')

allterms  = [('const', 'term_const'), ('a0', 'term_a0'), ('a1', 'term_a1') , ('a2', 'term_a2'), ('a3', 'term_a3'), ('a4', 'term_a4') , ('a5', 'term_a5') , ('a6', 'term_a6') , ('a7', 'term_a7')]
prefactor = [('const', [1., 0.]), ('a0', [10./3., 2./3.]), ('a1', [5., 0.]), ('a2', [10., 0.]), ('a3', [4., 0.]), ('a4', [4., 0.]), ('a5', [5., 0.]), ('a6', [5., 0.]), ('a7', [4., 0.])]
nbins_arg =  50
bins_arg = array('d', [2*i/float(nbins_arg)-1 for i in range(nbins_arg+1)])
#bins_arg = [2*i/float(nbins_arg)-1 for i in range(nbins_arg+1)]


hists = []

print('now booking the 3d histos')
for term,arg in allterms:
    nx, x = nbins_arg, bins_arg
    ny, y = nbins_yv, bins_yv
    nz, z = nbins_pt, bins_pt
    if term == 'const': nx, x = nbins_arg, array('d', [2.+i/float(nbins_arg)-1 for i in range(nbins_arg+1)])
    if term == 'a0'   : nx, x = nbins_arg, array('d', [-2.+3.*i/float(nbins_arg)-0 for i in range(nbins_arg+1)])

    #if term == 'const': nx, x = nbins_arg, [2.+i/float(nbins_arg)-1 for i in range(nbins_arg+1)]
    #if term == 'a0'   : nx, x = nbins_arg, [-2.+3.*i/float(nbins_arg)-0 for i in range(nbins_arg+1)]
    #tmp_histo  = ROOT.TH3D('h3_{v}VsRapVsPt' .format(v=term),'h3_{v}VsRapVsPt' .format(v=term), nx, x, ny, y, nz, z); tmp_histo .Sumw2()
    #dfcs.Histo3D( tmp_histo, arg, 'yV', 'ptV', 'genWeightReg')

    for i in range(100):
        histo_3d_tmp = ROOT.TH3D('h3_{v}VsRapVsPt{i}' .format(v=term,i=i),'h3_{v}VsRapVsPt' .format(v=term), nx, x, ny, y, nz, z)
        histo_model_tmp = ROOT.RDF.TH3DModel(histo_3d_tmp)
        histo_tmp = dfcs.Histo3D( histo_model_tmp, arg, 'yV', 'ptV', 'genWeightReg')
        #histo_tmp = dfcs.Histo3D( ('h3_'+term+'VsRapVsPt','h3_'+term+'VsRapVsPt', 50, -3., 2., 20, 0., 10., 20, 0., 100.), arg, 'yV', 'ptV', 'genWeightReg')

        hists.append(histo_tmp)

#ROOT.RDF.RunGraphs(hists)

lat = ROOT.TLatex(); lat.SetNDC(); lat.SetTextSize(0.05)
c = ROOT.TCanvas()                    
c.SetLeftMargin  (0.15)
c.SetRightMargin (0.15)
c.SetBottomMargin(0.15)

plotsdir='~/www/private/w-mass-13TeV/fractions/W/2020-11-19/'
plotsdir='./'

print('now doing the coefficients')

outfile = ROOT.TFile(plotsdir+'/fractions.root', 'recreate')
outfile.cd()

for ih,hist in enumerate(hists):
    if not ih:
        hist.Draw()
    h2_rapVsPt  = hist.Project3D('zy')                           
    formatHisto(h2_rapVsPt)

    nbinsX = hist.GetXaxis().GetNbins()
    nbinsY = hist.GetYaxis().GetNbins()
    nbinsZ = hist.GetZaxis().GetNbins()
    for iy in range(1,nbinsY+1):
        for iz in range(1,nbinsZ+1):
            c.Clear()
            name = '{iy}_{iz}'.format(iy=iy-1,iz=iz-1)
            #if options.verbose: 
            #   print('at bin {n}'.format(n=name))
            c.SetName ('canv_'+name)
            c.SetTitle('canv_'+name)
            h_cm = hist.ProjectionX(name, iy if iy else 1, iy if iy else nbinsY, iz if iz else 1, iz if iz else nbinsZ)
            #h_cm.SetTitle('cosTheta_'+name)
            ylo = hist.GetYaxis().GetBinLowEdge(iy)
            yhi = hist.GetYaxis().GetBinUpEdge(iy)
            plo = hist.GetZaxis().GetBinLowEdge(iz)
            phi = hist.GetZaxis().GetBinUpEdge(iz)
            h_cm.SetTitle('{arg}: y_{{V}} #in [{ylo:.2f},{yhi:.2f}] , p_{{T}}^{{V}} #in [{plo:.2f},{phi:.2f}] {ts}'.format(ylo=ylo,yhi=yhi,plo=plo,phi=phi,arg=prefactor[ih][0],ts='preFSR') )
            h_cm.GetXaxis().SetTitle('cos #theta^{*}x #phi')
            
            #h_cm_norm = h_cm.Clone(h_cm.GetName()+'_norm')
            #h_cm_norm.Scale(1./h_cm_norm.Integral())
            
            parameter = prefactor[ih][1][0]*h_cm.GetMean() + prefactor[ih][1][1]
            
            h2_rapVsPt.SetBinContent(iy, iz, parameter)
            
            h_cm.Draw()
            lat.DrawLatex(0.45, 0.68, '{term}: {par:.4f}'.format(term=prefactor[ih][0],par=parameter))
            #if options.verbose:
            #c.SaveAs('{pd}/{term}_argument_{n}.pdf'.format(pd=plotsdir,n=name,term=prefactor[ih][0]))
            #c.SaveAs('{pd}/{term}_argument_{n}.png'.format(pd=plotsdir,n=name,term=prefactor[ih][0]))

    h2_rapVsPt.Draw('colz')
    #c.SaveAs('{pd}/{term}_distribution_plus.pdf'.format(pd=plotsdir,n=name,term=prefactor[ih][0]))
    #c.SaveAs('{pd}/{term}_distribution_plus.png'.format(pd=plotsdir,n=name,term=prefactor[ih][0]))
    h2_rapVsPt.SetName(prefactor[ih][0]+'_plus')
    h2_rapVsPt.Write()


## skip for now cols  = ROOT.RDFDetail.ColumnNames_t()
## skip for now for c in ["costcs", "phics", "yV", "ptV", "mV", "phiV", "genWeight", "genWeightReg"]:
## skip for now     cols.push_back(c)
## skip for now dfcs.Snapshot("Events", "/eos/cms/store/cmst3/group/wmass/w-mass-13TeV/newtest.root", cols)

