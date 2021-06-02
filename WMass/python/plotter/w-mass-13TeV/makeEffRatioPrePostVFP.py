#!/usr/bin/env python3

# make SF for data/data and MC/MC from efficiency ratios

import re
import os, os.path
import logging
import argparse
import shutil

## safe batch mode
import sys
args = sys.argv[:]
sys.argv = ['-b']
import ROOT
sys.argv = args
ROOT.gROOT.SetBatch(True)
ROOT.PyConfig.IgnoreCommandLineOptions = True

from copy import *

sys.path.append(os.getcwd() + "/plotUtils/")
from utility import *

logging.basicConfig(level=logging.INFO)

def getAntiisoEfficiency(hiso, name):
    hantiiso = copy.deepcopy(hiso.Clone(name))
    hantiiso.Reset("ICESM")
    for ix in range(1, 1 + hiso.GetNbinsX()):
        for iy in range(1, 1 + hiso.GetNbinsX()):
            val = hiso.GetBinContent(ix, iy)
            err = hiso.GetBinError(  ix, iy)
            relerr = 1.0 if err == 0.0 else val/err
            hantiiso.SetBinContent(ix, iy, 1.0 - val)
            # keep same relative uncertainty
            hantiiso.SetBinError(  ix, iy, relerr * hantiiso.GetBinContent(ix, iy))
    return hantiiso

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("rootfile", type=str, nargs=1)
    parser.add_argument("outdir",   type=str, nargs=1)
    
    args = parser.parse_args()

    fname = args.rootfile[0]
    outdir = args.outdir[0] + "/"
    if outdir and not os.path.exists(outdir):
        os.makedirs(outdir)
        htmlpath = "./templates/index.php"
        shutil.copy(htmlpath, outdir)

    outname = outdir + os.path.basename(fname.replace(".root","_NEW.root"))
    
    histsEffMC =   {"pre" : {}, "post" : {} }
    histsEffData = {"pre" : {}, "post" : {}}

    adjustSettings_CMS_lumi()
    canvas = ROOT.TCanvas("canvas","",800,800)

    f = ROOT.TFile.Open(fname)
    if not f or not f.IsOpen():
        raise RuntimeError(f"Error when opening file {fname}")
    for k in f.GetListOfKeys():
        name = k.GetName()
        if not name.startswith("eff"): continue
        if "_BtoH_" in name: continue
        if not any(x in name for x in ["_BtoF_", "_GtoH"]): continue
        if "SF2D_Data_" in name or "SF2D_MC_" in name: continue # these will be created here, if already present let's skip them
        htmp = f.Get(name)
        htmp.SetTitle(name)
        drawCorrelationPlot(htmp,"Muon #eta","Muon p_{T} (GeV)","efficiency",
                            name,plotLabel="ForceTitle",outdir=outdir,
                            smoothPlot=False, drawProfileX=False, scaleToUnitArea=False,
                            draw_both0_noLog1_onlyLog2=1,passCanvas=canvas)
        era = ""
        if   "_BtoF_" in name: era = "pre"
        elif "_GtoH_" in name: era = "post"
        newname = name.replace("BtoF_","").replace("GtoH_","")
        if "effData" in name:
            newname = newname.replace("effData_","")
            histsEffData[era][newname] = copy.deepcopy(htmp.Clone("data_"+newname))
            histsEffData[era][newname].SetDirectory(0)
            #print(name + "   "+histsEffData[era][newname].GetName())
        elif "effMC" in name:
            newname = newname.replace("effMC_","")
            histsEffMC[era][newname] = copy.deepcopy(htmp.Clone("mc_"+newname))
            histsEffMC[era][newname].SetDirectory(0)
            #print(name + "   "+histsEffMC[era][newname].GetName())
    f.Close()
    
    print("----------------")

    # make efficiency ratios between eras for data and MC
    effKeys = list(histsEffMC["pre"].keys())

    # for antiiso the file only has SF, but not the efficiencies. Let's make them here
    # note that they are not saved for now
    
    #newEfficienciesToCopy = []
    for era in ["pre", "post"]:

        if not any(x == "antiiso_both" for x in effKeys):
            histsEffData[era]["antiiso_both"] = getAntiisoEfficiency(histsEffData[era]["iso_both"], "data_antiiso_both")
            histsEffMC[era]["antiiso_both"] = getAntiisoEfficiency(histsEffMC[era]["iso_both"], "mc_antiiso_both")
            #newEfficienciesToCopy.append(histsEffData[era]["antiiso_both"])
            #newEfficienciesToCopy.append(histsEffMC[era]["antiiso_both"])

        if not any(x == "antiisonotrig_both" for x in effKeys):
            histsEffData[era]["antiisonotrig_both"] = getAntiisoEfficiency(histsEffData[era]["isonotrig_both"], "data_antiisonotrig_both")
            histsEffMC[era]["antiisonotrig_both"] = getAntiisoEfficiency(histsEffMC[era]["isonotrig_both"], "mc_antiisonotrig_both")
            #newEfficienciesToCopy.append(histsEffData[era]["antiisonotrig_both"])
            #newEfficienciesToCopy.append(histsEffMC[era]["antiisonotrig_both"])
            
    # get keys again
    effKeys = list(histsEffMC["pre"].keys())
    
    hRatioData = {}
    hRatioMC = {}    

    for k in effKeys:
        #print(k)
        hRatioData[k] = copy.deepcopy(histsEffData["pre"][k].Clone("SF2D_Data_preOverPost_"+k))
        hRatioData[k].Divide(histsEffData["post"][k])
        hRatioMC[k] = copy.deepcopy(histsEffMC["pre"][k].Clone("SF2D_MC_preOverPost_"+k))
        hRatioMC[k].Divide(histsEffMC["post"][k])

        drawCorrelationPlot(hRatioData[k],"Muon #eta","Muon p_{T} (GeV)","Scale factor pre/post data::0.9,1.1",
                            f"effRatioData_PreOverPost_{k}",plotLabel="",outdir=outdir,
                            smoothPlot=False, drawProfileX=False, scaleToUnitArea=False,
                            draw_both0_noLog1_onlyLog2=1,passCanvas=canvas)
        drawCorrelationPlot(hRatioMC[k],"Muon #eta","Muon p_{T} (GeV)","Scale factor pre/post MC::0.9,1.1",
                            f"effRatioMC_PreOverPost_{k}",plotLabel="",outdir=outdir,
                            smoothPlot=False, drawProfileX=False, scaleToUnitArea=False,
                            draw_both0_noLog1_onlyLog2=1,passCanvas=canvas)

    originalFile = ROOT.TFile.Open(fname,"READ")
    if not originalFile or not originalFile.IsOpen():
        raise RuntimeError(f"Error when opening file {fname}")
        
    f = ROOT.TFile.Open(outname,"RECREATE")
    # copy all keys in original file, and add the new ones too
    if not f or not f.IsOpen():
        raise RuntimeError(f"Error when opening file {outname}")
    f.cd()
    for k in effKeys:
        hRatioData[k].Write()
        hRatioMC[k].Write()
    for k in originalFile.GetListOfKeys():
        name = k.GetName()
        obj  = k.ReadObj()
        if "SF2D_Data_" in name or "SF2D_MC_" in name: continue # do not copy the older version of the data/data or MC/MC scale factors if present
        if not obj:
            raise RuntimeError('Unable to read object {n}'.format(n=name))
        obj.Write(name)
    f.Close()
    originalFile.Close()
    print(effKeys)
    print(f"Everything was written in file\n{outname}\n")
