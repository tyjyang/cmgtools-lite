#!/usr/bin/env python3

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
        htmp = f.Get(name)
        if not name.startswith("eff"): continue
        if "_BtoH_" in name: continue
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
        if not obj:
            raise RuntimeError('Unable to read object {n}'.format(n=name))
        obj.Write(name)
    f.Close()
    originalFile.Close()
    print(f"Everything was written in file\n{outname}\n")
