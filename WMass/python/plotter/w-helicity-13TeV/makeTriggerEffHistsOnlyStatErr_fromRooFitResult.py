import ROOT, os, datetime, re, operator, math, sys
from array import array
ROOT.gROOT.SetBatch(True)

sys.path.append(os.getcwd() + "/plotUtils/")
from utility import *

ROOT.TH1.SetDefaultSumw2()


def formatBinBoundary(num):
    ret = str(round(num,2))
    if ret.startswith("-"): 
        ret = ret.replace("-","m")
    if "." in ret:        
        ret = ret.replace(".","p")
        digits = str(ret.split("p")[1])
        if len(digits) == 1:
            ret = ret + "0"
        elif len(digits) > 2:
            ret = str(ret.split("p")[0]) + "p" + digits[:2]
    else:
        ret = ret + "p00"
    return ret


# create a TH2 with data and MC efficiency from txt files like [0]. The txt file has this content
# columns from 1 to 4 are eta and pt bin borders, then there are the efficiency in data and MC, each one followed by stat.only uncertainty, then for MC there is syst (bu it is -1, so not available), then there is again data made with alternate fits with stat and syst uncertainty (both are -1)

# [0] https://mdunser.web.cern.ch/mdunser/private/w-helicity-13TeV/tnpFits/muFullData_trigger_fineBin_noMu50_PLUS/triggerMu/egammaEffi.txt

isMu = False
charge = "plus" if isMu else "all"
#lepton = "muon" 
lepton = "muon" if isMu else "electron" 
subdir = "" if isMu else "trigger_etaEE0p1/"
#mypath = "/afs/cern.ch/work/m/mciprian/w_mass_analysis/heppy/CMSSW_8_0_25/src/CMGTools/WMass/python/plotter/"
mypath = "/afs/cern.ch/work/m/mciprian/w_mass_analysis/heppy/CMSSW_8_0_25/src/CMGTools/WMass/python/postprocessing/data/leptonSF/new2016_madeSummer2018/TnPstuff/{lep}/{sub}".format(lep=lepton,sub=subdir)

etabins = []
ptbins = []
if isMu:
    infileData = mypath + "mu{ch}_Run2016_noMu50_all_triggerMu.nominalFit.root".format(ch="Minus" if charge == "minus" else "Plus")
    #infileMC = mypath + "mu{ch}_DY_noMu50_triggerMu.altSigFit.root".format(ch="Minus" if charge == "minus" else "Plus")
    # following MC file has TH1 with mass distribution, for Passign and Failing probes, can just use integral
    infileMC = mypath + "mu{ch}_DY_noMu50_triggerMu.root".format(ch="Minus" if charge == "minus" else "Plus")
    outfile = mypath + 'triggerMuonEff{ch}_fromRooFitResult_onlyStatUnc.root'.format(ch="Minus" if charge == "minus" else "Plus")
    plotoutdir = "/afs/cern.ch/user/m/mciprian/www/wmass/13TeV/test/muonTriggerSF_onlyStatUnc/"
    infile_help = mypath + "triggerMuonEff{ch}_onlyStatUnc.root".format(ch="Minus" if charge == "minus" else "Plus")
    etabins = [round(-2.40 + 0.10*x,1) for x in range(49)]
    ptbins = [round(x,1) for x in [25, 27.5, 30, 33, 36, 39, 45, 52, 55]]
else:
    infileData = mypath + "el_Run2016B_triggerEl.nominalFit.root"
    # following MC file has TH1 with mass distribution, for Passign and Failing probes, can just use integral
    infileMC = mypath + "el_DY_triggerEl.root"
    outfile = mypath + 'triggerElectronEffAllCharge_fromRooFitResult_onlyStatUnc.root'
    plotoutdir = "/afs/cern.ch/user/m/mciprian/www/wmass/13TeV/test/electronTriggerSF_onlyStatUnc/{sub}".format(sub=subdir)
    infile_help = mypath + "triggerElectronEffAllCharge_onlyStatUnc.root".format(ch="Minus" if charge == "minus" else "Plus")
    if "etaEE0p2" in subdir:
        posetabins = [round(0.10*x,1) for x in range(1,14)] + [1.444, 1.566, 1.7, 1.9, 2.1, 2.3, 2.5]
        etabins = [-1.0*x for x in reversed(posetabins)] + [0.0] + posetabins
        ptbins = [round(x,1) for x in [25, 27.5, 30, 31.5, 33, 36, 39, 42, 45, 48, 52, 55]]
    else:
        posetabins = [round(0.10*x,1) for x in range(17,26)]
        etabins = [-1.0*x for x in reversed(posetabins)] + [-1.566, 1.566] + posetabins
        ptbins = [round(x,1) for x in [30, 33, 36, 40, 45]]

print etabins
print ptbins

# infile_help files used to get uncertainty when RooFitResults yields crazy numbers: they contain the efficiencies and uncertainties from the txt stored in TH2 (might need to create it with makeTriggerEffHistsOnlyStatErr.py 

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
nFits = 0
for k in tf.GetListOfKeys():
    obj=k.ReadObj()
    if obj.ClassName() == "RooFitResult": nFits += 1
nFits = nFits/2  # there are objects for pass and fail, but the bin index is the same
print "There are {n} RooFitResult objects in file".format(n=nFits)

for ipt in range(len(ptbins)-1):
    for ieta in range(len(etabins)-1):
        ibin = ipt * (len(etabins) - 1) + ieta
        if ibin < 10: 
            ibin = ("0" if nFits <= 100 else "00") + str(ibin)
        elif ibin < 100: 
            ibin  = ("" if nFits <= 100 else "0") + str(ibin)
        # Name should be something like bin383_probe_lep_eta_2p30To2p40_probe_lep_pt_52p00To55p00_resP
        el = formatBinBoundary(etabins[ieta])
        eh = formatBinBoundary(etabins[ieta+1])
        pl = formatBinBoundary(ptbins[ipt])
        ph = formatBinBoundary(ptbins[ipt+1])       
        frPname = "bin{n}_probe_lep_eta_{el}To{eh}_probe_lep_pt_{pl}To{ph}_resP".format(n=ibin,el=el,eh=eh,pl=pl,ph=ph)
        frFname = "bin{n}_probe_lep_eta_{el}To{eh}_probe_lep_pt_{pl}To{ph}_resF".format(n=ibin,el=el,eh=eh,pl=pl,ph=ph)
        frP = tf.Get(frPname)
        if not frP: 
            print "Error: object %s not found in file %s" % (frPname, infileData)
            quit()
        frF = tf.Get(frFname)
        if not frF: 
            print "Error: object %s not found in file %s" % (frFname, infileData)
            quit()
        # get RooRealVar with nSigP or nSigF
        rrv_npass = frP.floatParsFinal().find("nSigP")
        rrv_nfail = frF.floatParsFinal().find("nSigF")
        nPass = rrv_npass.getValV()
        nFail = rrv_nfail.getValV()
        errPass = rrv_npass.getError()
        errFail = rrv_nfail.getError()
        #print "nPass = %s +/ %s    nFail = %s +/ %s" % (rrv_npass.getValV(), errPass, rrv_nfail.getValV(), errFail)
        if nPass > 0:
            hPassData.SetBinContent(ieta+1,nPass)                
            hPassData.SetBinError(ieta+1, errPass)        
        else:
            hPassData.SetBinContent(ieta+1,0.0)                
            hPassData.SetBinError(ieta+1,  1.0)                    
        if nFail > 0:
            hFailData.SetBinContent(ieta+1,nFail)        
            hFailData.SetBinError(ieta+1, errFail)        
        else:
            hFailData.SetBinContent(ieta+1, 1.0)       
            hFailData.SetBinError(ieta+1,   1.0)

    #print "%s < eta < %s    %s < pt < %s " % (el, eh, pl, ph)
    hTotData.Add(hPassData, hFailData)
    gr = ROOT.TGraphAsymmErrors(hPassData, hTotData, "cl=0.683 b(1,1) mode") 
    for i in range(gr.GetN()):
        if hPassData.GetBinContent(i+1) > 0:
            heffData.SetBinContent(i+1, ipt+1, hPassData.GetBinContent(i+1)/hTotData.GetBinContent(i+1))
            uncertainty = gr.GetErrorY(i)  # or max(gr.GetErrorYhigh(i),gr.GetErrorYlow(i)) ?
        else:
            heffData.SetBinContent(i+1, ipt+1, 1.0)
            uncertainty = 1.0
        if isMu:
            if uncertainty > 0.01: 
                uncertainty = min(uncertainty, hEffData_help.GetBinError(i+1, ipt+1))
        else:
            if uncertainty > 0.02 and uncertainty < 1.0:
                uncertainty = min(uncertainty, hEffData_help.GetBinError(i+1, ipt+1))
        heffData.SetBinError(i+1, ipt+1, uncertainty)
        heffUncRatio_DataOverMC.SetBinContent(i+1, ipt+1, uncertainty)

tf.Close()
print "Created file %s" % infileData

# ugly to repeat the same code, but faster to manage for now
tf = ROOT.TFile(infileMC)
if not tf or not tf.IsOpen():
    print "Error opening file %s" % infileMC
    quit()
nFits = 0
for k in tf.GetListOfKeys():
    obj=k.ReadObj()
    if obj.ClassName() == "RooFitResult": nFits += 1
nFits = nFits/2

for ipt in range(len(ptbins)-1):
    for ieta in range(len(etabins)-1):
        ibin = ipt * (len(etabins) - 1) + ieta
        if ibin < 10: 
            ibin = ("0" if nFits <= 100 else "00") + str(ibin)
        elif ibin < 100: 
            ibin  = ("" if nFits <= 100 else "0") + str(ibin)
        # name should be something like bin383_probe_lep_eta_2p30To2p40_probe_lep_pt_52p00To55p00_Pass or Fail
        el = formatBinBoundary(etabins[ieta])
        eh = formatBinBoundary(etabins[ieta+1])
        pl = formatBinBoundary(ptbins[ipt])
        ph = formatBinBoundary(ptbins[ipt+1])
        frPname = "bin{n}_probe_lep_eta_{el}To{eh}_probe_lep_pt_{pl}To{ph}_Pass".format(n=ibin,el=el,eh=eh,pl=pl,ph=ph)
        frFname = "bin{n}_probe_lep_eta_{el}To{eh}_probe_lep_pt_{pl}To{ph}_Fail".format(n=ibin,el=el,eh=eh,pl=pl,ph=ph)
        frP = tf.Get(frPname)
        if not frP: 
            print "Error: object %s not found in file %s" % (frPname, infileMC)
            quit()
        frF = tf.Get(frFname)
        if not frF: 
            print "Error: object %s not found in file %s" % (frFname, infileMC)
            quit()
        # get histogram and retrieve integral and uncertainty (sqrt(N))
        npass = frP.Integral()
        nfail = frF.Integral()
        hPassMC.SetBinContent(ieta+1,npass if npass > 0 else 1.0)        
        hPassMC.SetBinError(ieta+1, math.sqrt(npass) if npass > 0 else 1.0)        
        hFailMC.SetBinContent(ieta+1,nfail if nfail > 0 else 1.0)        
        hFailMC.SetBinError(ieta+1, math.sqrt(nfail) if nfail > 0 else 1.0)        
        
    hTotMC.Add(hPassMC, hFailMC)
    gr = ROOT.TGraphAsymmErrors(hPassMC, hTotMC, "cl=0.683 b(1,1) mode") 
    for i in range(gr.GetN()):
        if hTotMC.GetBinContent(i+1) > 0:
            heffMC.SetBinContent(i+1, ipt+1, hPassMC.GetBinContent(i+1)/hTotMC.GetBinContent(i+1))
        else:
            heffMC.SetBinContent(i+1, ipt+1, 1.0)
        uncertainty = gr.GetErrorY(i)  # or max(gr.GetErrorYhigh(i),gr.GetErrorYlow(i)) ?
        #if uncertainty > 0.01: uncertainty = min(uncertainty, hEffMC_help.GetBinError(i+1, ipt+1))
        heffMC.SetBinError(i+1, ipt+1, uncertainty)
        htmp.SetBinContent(i+1, ipt+1, uncertainty)

tf.Close()
print "Created file %s" % infileMC

heffSF.Divide(heffData,heffMC)  # uncertainties are uncorrelated, can use simple Divide to propagate uncertainty on SF
heffUncRatio_DataOverMC.Divide(htmp)

createPlotDirAndCopyPhp(plotoutdir)
xaxisname = "{lep} #eta".format(lep=lepton)
yaxisname = "{lep} p_{{T}} [GeV]".format(lep=lepton)

drawCorrelationPlot(heffUncRatio_DataOverMC,xaxisname,yaxisname,
                    "#sigma(stat)_{data} / #sigma(stat)_{MC}: %s trigger::0,1.0" % ("#mu" if isMu else "e"),
                    heffUncRatio_DataOverMC.GetName(),"ForceTitle",plotoutdir,0,0,False,False,False,
                    1,palette=55)

drawCorrelationPlot(heffData,xaxisname,yaxisname,
                    "#sigma(stat)_{data} for %s trigger efficiency::0,0.02" % ("#mu" if isMu else "e"),
                    heffData.GetName()+"_statUncertainty","ForceTitle",plotoutdir,0,0,False,False,False,
                    1,palette=55,plotError=True)

drawCorrelationPlot(heffMC,xaxisname,yaxisname,
                    "#sigma(stat)_{MC} for %s trigger efficiency::0,0.02" % ("#mu" if isMu else "e"),
                    heffMC.GetName()+"_statUncertainty","ForceTitle",plotoutdir,0,0,False,False,False,
                    1,palette=55,plotError=True)

tfout = ROOT.TFile.Open(outfile,'recreate')
heffData.Write(heffData.GetName())
heffMC.Write(heffMC.GetName())
heffSF.Write(heffSF.GetName())
#heffData_statUncOverTot.Write(heffData_statUncOverTot.GetName())
tfout.Close()
print "Created file %s" % outfile

