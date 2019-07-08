import ROOT, os, datetime, re, operator, math, sys
from array import array
ROOT.gROOT.SetBatch(True)

sys.path.append(os.getcwd() + "/plotUtils/")
from utility import *

ROOT.TH1.SetDefaultSumw2()


# create a TH2 with data and MC efficiency from txt files like [0]. The txt file has this content
# columns from 1 to 4 are eta and pt bin borders, then there are the efficiency in data and MC, each one followed by stat.only uncertainty, then for MC there is syst (bu it is -1, so not available), then there is again data made with alternate fits with stat and syst uncertainty (both are -1)

# [0] https://mdunser.web.cern.ch/mdunser/private/w-helicity-13TeV/tnpFits/muFullData_trigger_fineBin_noMu50_PLUS/triggerMu/egammaEffi.txt

charge = "plus"
#mypath = "/afs/cern.ch/work/m/mciprian/w_mass_analysis/heppy/CMSSW_8_0_25/src/CMGTools/WMass/python/plotter/"
mypath = "/afs/cern.ch/work/m/mciprian/w_mass_analysis/heppy/CMSSW_8_0_25/src/CMGTools/WMass/python/postprocessing/data/leptonSF/new2016_madeSummer2018/"
infileData = mypath + "mu{ch}_Run2016_noMu50_all_triggerMu.nominalFit.root".format(ch="Minus" if charge == "minus" else "Plus")
infileMC = mypath + "mu{ch}_DY_noMu50_triggerMu.altSigFit.root".format(ch="Minus" if charge == "minus" else "Plus")
outfile = mypath + 'triggerMuonEff{ch}_fromRooFitResult_onlyStatUnc.root'.format(ch="Minus" if charge == "minus" else "Plus")
plotoutdir = "/afs/cern.ch/user/m/mciprian/www/wmass/13TeV/test/muonTriggerSF_onlyStatUnc/"

# files used to get uncertainty when RooFitResults yields crazy numbers: these files below contain the uncertainties from the txt stored in TH2
infile_help = mypath + "triggerMuonEff{ch}_onlyStatUnc.root".format(ch="Minus" if charge == "minus" else "Plus")
hEffData_help = None
hEffMC_help = None
tf = ROOT.TFile(infile_help)
if not tf or not tf.IsOpen():
    print "Error opening file %s" % infile_help
    quit()
else:
    hEffData_help = tf.Get("effData_%s" % charge)
    hEffMC_help = tf.Get("effMC_%s" % charge)
    if not hEffData_help or not hEffMC_help:
        print "Error while getting histograms from file %s" % infile_help
        quit()
    else:
        hEffData_help.SetDirectory(0)
        hEffMC_help.SetDirectory(0)
    tf.Close()
        
etabins = [round(-2.40 + 0.10*x,1) for x in range(49)]
ptbins = [round(float(x),1) for x in [25, 27.5, 30, 33, 36, 39, 45, 52, 55]]
print etabins
print ptbins


hPassData = ROOT.TH1D("passDataEta_{ch}".format(ch=charge),"charge {ch}".format(ch=charge),len(etabins)-1, array("d",etabins))
hFailData = ROOT.TH1D("failDataEta_{ch}".format(ch=charge),"charge {ch}".format(ch=charge),len(etabins)-1, array("d",etabins))
hTotData = ROOT.TH1D("totDataEta_{ch}".format(ch=charge),"charge {ch}".format(ch=charge),len(etabins)-1, array("d",etabins))
heffData = ROOT.TH2D("effData_{ch}".format(ch=charge),"charge {ch}".format(ch=charge),len(etabins)-1, array("d",etabins), len(ptbins)-1, array("d",ptbins))

hPassMC = ROOT.TH1D("passMCEta_{ch}".format(ch=charge),"charge {ch}".format(ch=charge),len(etabins)-1, array("d",etabins))
hFailMC = ROOT.TH1D("failMCEta_{ch}".format(ch=charge),"charge {ch}".format(ch=charge),len(etabins)-1, array("d",etabins))
hTotMC = ROOT.TH1D("totMCEta_{ch}".format(ch=charge),"charge {ch}".format(ch=charge),len(etabins)-1, array("d",etabins))
heffMC = ROOT.TH2D("effMC_{ch}".format(ch=charge),"charge {ch}".format(ch=charge),len(etabins)-1, array("d",etabins), len(ptbins)-1, array("d",ptbins))
heffSF = ROOT.TH2D("triggerSF_{ch}".format(ch=charge),"charge {ch}".format(ch=charge),len(etabins)-1, array("d",etabins), len(ptbins)-1, array("d",ptbins))

heffUncRatio_DataOverMC = ROOT.TH2D("heffUncRatio_DataOverMC_{ch}".format(ch=charge),"charge {ch}".format(ch=charge),len(etabins)-1, array("d",etabins), len(ptbins)-1, array("d",ptbins))
htmp = ROOT.TH2D("htmp_{ch}".format(ch=charge),"charge {ch}".format(ch=charge),len(etabins)-1, array("d",etabins), len(ptbins)-1, array("d",ptbins))

#heffData_statUncOverTot = ROOT.TH2D("effData_statUncOverTot_{ch}".format(ch=charge),"charge {ch}".format(ch=charge),
#                                    len(etabins)-1, array("d",etabins), len(ptbins)-1, array("d",ptbins))

# start with data
tf = ROOT.TFile(infileData)
if not tf or not tf.IsOpen():
    print "Error opening file %s" % infileData
    quit()

for ipt in range(len(ptbins)-1):
    for ieta in range(len(etabins)-1):
        ibin = ipt * (len(etabins) - 1) + ieta
        if ibin < 10: ibin = "00" + str(ibin)
        elif ibin < 100: ibin  = "0" + str(ibin)
        # Name should be something like bin383_probe_lep_eta_2p30To2p40_probe_lep_pt_52p00To55p00_resP
        el = str(etabins[ieta]).replace('.','p').replace('-','m')+"0"
        eh = str(etabins[ieta+1]).replace('.','p').replace('-','m')+"0"
        pl = str(ptbins[ipt]).replace('.','p').replace('-','m')+"0"
        ph = str(ptbins[ipt+1]).replace('.','p').replace('-','m')+"0"
        frPname = "bin{n}_probe_lep_eta_{el}To{eh}_probe_lep_pt_{pl}To{ph}_resP".format(n=ibin,el=el,eh=eh,pl=pl,ph=ph)
        frFname = "bin{n}_probe_lep_eta_{el}To{eh}_probe_lep_pt_{pl}To{ph}_resF".format(n=ibin,el=el,eh=eh,pl=pl,ph=ph)
        frP = tf.Get(frPname)
        if not frP: 
            print "Error: object %s not found in file %s" % (frPname, infile)
            quit()
        frF = tf.Get(frFname)
        if not frF: 
            print "Error: object %s not found in file %s" % (frFname, infile)
            quit()
        # get RooRealVar with nSigP or nSigF
        rrv_npass = frP.floatParsFinal().find("nSigP")
        rrv_nfail = frF.floatParsFinal().find("nSigF")
        hPassData.SetBinContent(ieta+1,rrv_npass.getValV())        
        hPassData.SetBinError(ieta+1, rrv_npass.getError())        
        hFailData.SetBinContent(ieta+1,rrv_nfail.getValV())        
        hFailData.SetBinError(ieta+1, rrv_nfail.getError())        
        
    hTotData.Add(hPassData, hFailData)
    gr = ROOT.TGraphAsymmErrors(hPassData, hTotData, "cl=0.683 b(1,1) mode") 
    for i in range(gr.GetN()):
        heffData.SetBinContent(i+1, ipt+1, hPassData.GetBinContent(i+1)/hTotData.GetBinContent(i+1))
        uncertainty = gr.GetErrorY(i)  # or max(gr.GetErrorYhigh(i),gr.GetErrorYlow(i)) ?
        if uncertainty > 0.01: uncertainty = min(uncertainty, hEffData_help.GetBinError(i+1, ipt+1))
        heffData.SetBinError(i+1, ipt+1, uncertainty)
        heffUncRatio_DataOverMC.SetBinContent(i+1, ipt+1, uncertainty)

tf.Close()
print "Created file %s" % infileData

# ugly to repeat the same code, but faster to manage for now
tf = ROOT.TFile(infileMC)
if not tf or not tf.IsOpen():
    print "Error opening file %s" % infileMC
    quit()

for ipt in range(len(ptbins)-1):
    for ieta in range(len(etabins)-1):
        ibin = ipt * (len(etabins) - 1) + ieta
        if ibin < 10: ibin = "00" + str(ibin)
        elif ibin < 100: ibin  = "0" + str(ibin)
        # name should be something like bin383_probe_lep_eta_2p30To2p40_probe_lep_pt_52p00To55p00_resP
        el = str(etabins[ieta]).replace('.','p').replace('-','m')+"0"
        eh = str(etabins[ieta+1]).replace('.','p').replace('-','m')+"0"
        pl = str(ptbins[ipt]).replace('.','p').replace('-','m')+"0"
        ph = str(ptbins[ipt+1]).replace('.','p').replace('-','m')+"0"
        frPname = "bin{n}_probe_lep_eta_{el}To{eh}_probe_lep_pt_{pl}To{ph}_resP".format(n=ibin,el=el,eh=eh,pl=pl,ph=ph)
        frFname = "bin{n}_probe_lep_eta_{el}To{eh}_probe_lep_pt_{pl}To{ph}_resF".format(n=ibin,el=el,eh=eh,pl=pl,ph=ph)
        frP = tf.Get(frPname)
        if not frP: 
            print "Error: object %s not found in file %s" % (frPname, infile)
            quit()
        frF = tf.Get(frFname)
        if not frF: 
            print "Error: object %s not found in file %s" % (frFname, infile)
            quit()
        # get RooRealVar with nSigP or nSigF
        rrv_npass = frP.floatParsFinal().find("nSigP")
        rrv_nfail = frF.floatParsFinal().find("nSigF")
        hPassMC.SetBinContent(ieta+1,rrv_npass.getValV())        
        hPassMC.SetBinError(ieta+1, rrv_npass.getError())        
        hFailMC.SetBinContent(ieta+1,rrv_nfail.getValV())        
        hFailMC.SetBinError(ieta+1, rrv_nfail.getError())        
        
    hTotMC.Add(hPassMC, hFailMC)
    gr = ROOT.TGraphAsymmErrors(hPassMC, hTotMC, "cl=0.683 b(1,1) mode") 
    for i in range(gr.GetN()):
        heffMC.SetBinContent(i+1, ipt+1, hPassMC.GetBinContent(i+1)/hTotMC.GetBinContent(i+1))
        uncertainty = gr.GetErrorY(i)  # or max(gr.GetErrorYhigh(i),gr.GetErrorYlow(i)) ?
        if uncertainty > 0.01: uncertainty = min(uncertainty, hEffMC_help.GetBinError(i+1, ipt+1))
        heffMC.SetBinError(i+1, ipt+1, uncertainty)
        htmp.SetBinContent(i+1, ipt+1, uncertainty)

tf.Close()
print "Created file %s" % infileMC

heffSF.Divide(heffData,heffMC)  # uncertainties are uncorrelated, can use simple Divide to propagate uncertainty on SF
heffUncRatio_DataOverMC.Divide(htmp)

createPlotDirAndCopyPhp(plotoutdir)
drawCorrelationPlot(heffUncRatio_DataOverMC,"muon #eta","muon p_{T} [GeV]",
                    "#sigma(stat)_{data} / #sigma(stat)_{MC}: #mu trigger::0,1.0",
                    heffUncRatio_DataOverMC.GetName(),"ForceTitle",plotoutdir,0,0,False,False,False,
                    1,palette=55)

drawCorrelationPlot(heffData,"muon #eta","muon p_{T} [GeV]",
                    "#sigma(stat)_{data} for #mu trigger efficiency::0,0.02",
                    heffData.GetName()+"_statUncertainty","ForceTitle",plotoutdir,0,0,False,False,False,
                    1,palette=55,plotError=True)

drawCorrelationPlot(heffMC,"muon #eta","muon p_{T} [GeV]",
                    "#sigma(stat)_{MC} for #mu trigger efficiency::0,0.02",
                    heffMC.GetName()+"_statUncertainty","ForceTitle",plotoutdir,0,0,False,False,False,
                    1,palette=55,plotError=True)

tfout = ROOT.TFile.Open(outfile,'recreate')
heffData.Write(heffData.GetName())
heffMC.Write(heffMC.GetName())
heffSF.Write(heffSF.GetName())
#heffData_statUncOverTot.Write(heffData_statUncOverTot.GetName())
tfout.Close()
print "Created file %s" % outfile

