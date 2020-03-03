import ROOT, os, datetime, re, operator, math, sys
from array import array
ROOT.gROOT.SetBatch(True)

sys.path.append(os.getcwd() + "/plotUtils/")
from utility import *

ROOT.TH1.SetDefaultSumw2()


# create a TH2 with data and MC efficiency from txt files like [0]. The txt file has this content
# columns from 1 to 4 are eta and pt bin borders, then there are the efficiency in data and MC, each one followed by stat.only uncertainty, then for MC there is syst (bu it is -1, so not available), then there is again data made with alternate fits with stat and syst uncertainty (both are -1)

# [0] https://mdunser.web.cern.ch/mdunser/private/w-helicity-13TeV/tnpFits/muFullData_trigger_fineBin_noMu50_PLUS/triggerMu/egammaEffi.txt

isMu = True
charge = "minus" if isMu else "all"
#lepton = "muon" 
lepton = "muon" if isMu else "electron" 
subdir = "" if isMu else "trigger_etaEE0p1/"
#mypath = "/afs/cern.ch/work/m/mciprian/w_mass_analysis/heppy/CMSSW_8_0_25/src/CMGTools/WMass/python/plotter/"
mypath = "/afs/cern.ch/work/m/mciprian/w_mass_analysis/heppy/CMSSW_8_0_25/src/CMGTools/WMass/python/postprocessing/data/leptonSF/new2016_madeSummer2018/TnPstuff/{lep}/{sub}".format(lep=lepton,sub=subdir)
infile = mypath + "trigger{lep}Eff{ch}.txt".format(lep="Muon" if isMu else "Electron", ch="Minus" if charge == "minus" else "Plus" if charge == "plus" else "AllCharge")
outfile = infile.replace('.txt','_onlyStatUnc.root')
plotoutdir = "/afs/cern.ch/user/m/mciprian/www/wmass/13TeV/test/{lep}/trigger_{ch}/{sub}".format(lep=lepton,sub=subdir,ch=charge)

noAlternativeFitData = False
if subdir == "trigger_etaEE0p1/":
    noAlternativeFitData = True
    
etabins = []
ptbins = []
if isMu:
    etabins = [round(-2.40 + 0.10*x,1) for x in range(49)]
    ptbins = [25, 27.5, 30, 33, 36, 39, 45, 52, 55]
else:
    if "etaEE0p2" in subdir:
        posetabins = [round(0.10*x,1) for x in range(1,14)] + [1.444, 1.566, 1.7, 1.9, 2.1, 2.3, 2.5]
        etabins = [-1.0*x for x in reversed(posetabins)] + [0.0] + posetabins
        ptbins = [25, 27.5, 30, 31.5, 33, 36, 39, 42, 45, 48, 52, 55] 
    else:
        posetabins = [round(0.10*x,1) for x in range(17,26)]
        etabins = [-1.0*x for x in reversed(posetabins)] + [-1.566, 1.566] + posetabins
        ptbins = [30, 33, 36, 40, 45]

print etabins
print ptbins

heffData = ROOT.TH2D("effData_{ch}".format(ch=charge),"",len(etabins)-1, array("d",etabins), len(ptbins)-1, array("d",ptbins))
heffMC = ROOT.TH2D("effMC_{ch}".format(ch=charge),"",len(etabins)-1, array("d",etabins), len(ptbins)-1, array("d",ptbins))
heffSF = ROOT.TH2D("triggerSF_{ch}".format(ch=charge),"",len(etabins)-1, array("d",etabins), len(ptbins)-1, array("d",ptbins))

heffData_statUncOverTot = ROOT.TH2D("effData_statUncOverTot_{ch}".format(ch=charge),"charge {ch}".format(ch=charge),
                                    len(etabins)-1, array("d",etabins), len(ptbins)-1, array("d",ptbins))
heffData_systUnc = ROOT.TH2D("effData_systUnc_{ch}".format(ch=charge),"charge {ch}".format(ch=charge),
                             len(etabins)-1, array("d",etabins), len(ptbins)-1, array("d",ptbins))


with open(infile) as f:
    lines = f.readlines()
    for line in lines:
        if line.startswith('###'): continue
        columns = line.split()
        etabincenter = (float(columns[0]) + float(columns[1]))/2.
        ptbincenter = (float(columns[2]) + float(columns[3]))/2.
        ieta = heffData.GetXaxis().FindFixBin(etabincenter)
        ipt  = heffData.GetYaxis().FindFixBin(ptbincenter)
        heffData.SetBinContent(ieta,ipt,float(columns[4]))
        heffData.SetBinError(ieta,ipt,float(columns[5]))
        heffMC.SetBinContent(ieta,ipt,float(columns[6]))
        heffMC.SetBinError(ieta,ipt,float(columns[7]))
        systUncData = abs(float(columns[9]) - float(columns[4]))
        if noAlternativeFitData: 
            systUncData = 0.0
        totUncData = math.sqrt( float(columns[5]) * float(columns[5]) + systUncData * systUncData)
        heffData_statUncOverTot.SetBinContent(ieta,ipt,float(columns[5])/totUncData)
        heffData_systUnc.SetBinContent(ieta,ipt,systUncData)

heffSF.Divide(heffData,heffMC)  # uncertainties are uncorrelated, can use simple Divide to propagate uncertainty on SF

drawCorrelationPlot(heffData_statUncOverTot,"%s #eta" % lepton, "%s p_{T} [GeV]" % lepton,
                    "#sigma(stat) / #sigma(Tot): {l} trigger data efficiency".format(l="#mu" if isMu else "e"),
                    heffData_statUncOverTot.GetName(),"ForceTitle",plotoutdir,0,0,False,False,False,
                    1,palette=55)
drawCorrelationPlot(heffData_systUnc,"%s #eta" % lepton, "%s p_{T} [GeV]" % lepton,
                    "#sigma(syst): {l} trigger data efficiency".format(l="#mu" if isMu else "e"),
                    heffData_systUnc.GetName(),"ForceTitle",plotoutdir,0,0,False,False,False,
                    1,palette=55)

tf = ROOT.TFile.Open(outfile,'recreate')
heffData.Write(heffData.GetName())
heffMC.Write(heffMC.GetName())
heffSF.Write(heffSF.GetName())
heffData_statUncOverTot.Write(heffData_statUncOverTot.GetName())
heffData_systUnc.Write(heffData_systUnc.GetName())
tf.Close()
print "Created file %s " % outfile
