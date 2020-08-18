import ROOT, os, datetime, re, operator, math, sys
from array import array
ROOT.gROOT.SetBatch(True)

from make_diff_xsec_cards import get_ieta_ipt_from_process_name
from make_diff_xsec_cards import get_ieta_from_process_name
from make_diff_xsec_cards import get_ipt_from_process_name
from make_diff_xsec_cards import getDiffXsecBinning
from make_diff_xsec_cards import templateBinning

import utilities
utilities = utilities.util()

## ===================================================================
## USAGE:
## needs as infile a toys.root with limit tree from toys
## takes a comma separated list of regular expressions as input via --params
## if no output directory is given, it will just plot the smaller correlation matrix
## if output directory is given, it will save it there as pdf and png

## example:
## python w-helicity-13TeV/subMatrix.py toys.root --params alph,muR,muF,.*Ybin.*2,pdf12,pdf56,pdf42 --outdir <output_directory> --type toys/hessian
## examples of common regexps:
## NORMs: 'norm_.*'
## PDFs: range [1,20]: '^pdf([1-9]|1[0-9]|20)$'
## ===================================================================
#def SetCorrMatrixPalette():
#ROOT.gStyle.SetPalette()

def getBinAreaFromParamName(pname,genBins,isHelicity=False):
    if isHelicity:
        genYwBins = genBins
        dyw = 1.0
        # use left pol to retrieve binning if none of the keys below is found (e.g. it happens for unpolarized)
        # similarly for charge, use minus if no match
        # this works because the binning is the same
        charge = "plus" if "plus" in pname else "minus" 
        pol = "left" if 'left' in pname else "right" if "right" in pname else "long" if "long" in pname else "left"
        charge_pol = "{ch}_{pol}".format(ch=charge,pol=pol)  
        if "_Ybin_" in pname:
            iyw = int((pname.split("_Ybin_")[1]).split("_")[0])
            dyw = genYwBins[charge_pol][iyw+1] - genYwBins[charge_pol][iyw]
        return dyw
    else:
        genEtaPtBins = genBins
        dpt = 1.0
        deta = 1.0
        if "_ieta_" in pname:        
            ieta = int((pname.split("_ieta_")[1]).split("_")[0])
            deta = genEtaPtBins.etaBins[ieta+1] - genEtaPtBins.etaBins[ieta]
        if "_ipt_" in pname:        
            ipt = int((pname.split("_ipt_")[1]).split("_")[0])
            dpt = genEtaPtBins.ptBins[ipt+1] - genEtaPtBins.ptBins[ipt]
        return deta * dpt

def lepInFakeSystForSort(name):
    if re.match("Uncorrelated\d+mu",name): 
        return 1
    else:
        return 2


def niceName(name,genBins="",forceLep="",drawRangeInLabel=False):

    if re.match("W.*sumxsec",name) and all(x not in name for x in ["_Ybin", "_ieta_", "_ipt_"]):
        # inclusive xsec poi
        nn  = '#mu' if '_mu_' in name else 'el' if '_el_' in name else 'l'
        if forceLep: nn = '#mu' if forceLep == "mu" else 'e' if forceLep == "el" else 'l'
        chs = "+" if "plus" in name else "-" if "minus" in name else ""
        return "W{ch}({l}#nu) fiducial".format(ch=chs,l=nn)

    if re.match("Wasym.*chargemetaasym",name) and all(x not in name for x in ["_Ybin", "_ieta_", "_ipt_"]):
        # inclusive xsec poi
        nn  = '#mu' if '_mu_' in name else 'el' if '_el_' in name else 'l'
        if forceLep: nn = '#mu' if forceLep == "mu" else 'e' if forceLep == "el" else 'l'
        return "W({l}#nu) (asymmetry) fiducial".format(l=nn)

    if re.match("Wratio.*ratiometaratio",name) and all(x not in name for x in ["_Ybin", "_ieta_", "_ipt_"]):
        # inclusive xsec poi
        nn  = '#mu' if '_mu_' in name else 'el' if '_el_' in name else 'l'
        if forceLep: nn = '#mu' if forceLep == "mu" else 'e' if forceLep == "el" else 'l'
        return "W+/W- fiducial".format(l=nn)


    if '_Ybin' in name:
        genYwBins = genBins
        if drawRangeInLabel:
            charge = "plus" if "plus" in name else "minus" if "minus" in name else ""
            wch = "W+" if "plus" in name else "W-" if "minus" in name else "W"
            lep = '#mu' if '_mu_' in name else 'el' if '_el_' in name else 'l'
            if forceLep: lep = '#mu' if forceLep == "mu" else 'e' if forceLep == "el" else 'l'
            pol = "left" if 'left' in name else "right" if "right" in name else "long" if "long" in name else "unpolarized"
            ywl,ywh = 0.0,0.0
            if genYwBins:
                iyw = int((name.split("_Ybin_")[1]).split("_")[0])
                # get key to read binning (same for all charges and polarizations, but let's stay general
                # in case charge was not set, take key for "plus" charge and cross fingers
                # similarly, use left for binning when pol is unpolarized
                charge_pol = "{ch}_{pol}".format(ch=charge if charge != "" else "plus",pol=pol if pol!="unpolarized" else "left")  
                ywl = genYwBins[charge_pol][iyw]
                ywh = genYwBins[charge_pol][iyw+1]
            #nn = "{wch}#rightarrow{lep}#nu {pol}: |y_{{W}}| #in [{ywl:1.2f},{ywh:1.2f}]".format(wch=wch,lep=lep,pol=pol,ywl=ywl,ywh=ywh)
            nn = "{wch} {pol}: |y_{{W}}| #in [{ywl:1.2f},{ywh:1.2f}]".format(wch=wch,pol=pol,ywl=ywl,ywh=ywh)
            
            if name.startswith("norm_"):  # normalization systematics for bins not fitted as pois
                nn = "norm.syst. " + nn
        else:
            nn  = '#mu: ' if '_mu_' in name else 'el: ' if '_el_' in name else ''
            if 'plus' in name: nn += 'W+ '
            elif 'minus' in name: nn += 'W- '
            else: nn += 'W '
            if 'left' in name: nn += 'left '
            elif 'right' in name: nn += 'right '
            elif 'long' in name: nn += 'long '
            else: nn += 'unpolarized '
            idx = -2 if (name.endswith('mu') or any([x in name for x in ['masked','sumxsec','charge','a0','a4']])) else -1
            nn += name.split('_')[idx]
            if 'eff_unc' in name:
                nn = '#epsilon_{unc}^{'+nn+'}'
        return nn

    elif '_ieta_' in name and '_ipt_' in name:
        genEtaPtBins = genBins
        ieta,ipt = get_ieta_ipt_from_process_name(name)
        if drawRangeInLabel:
            etal,etah = 0.0,0.0
            ptl,pth = 0.0,0.0
            if genEtaPtBins:
                etal = genEtaPtBins.etaBins[ieta]
                etah = genEtaPtBins.etaBins[ieta+1]
                ptl = genEtaPtBins.ptBins[ipt]
                pth = genEtaPtBins.ptBins[ipt+1]
            wch = 'W+' if 'plus' in name else 'W-' if 'minus' in name else 'W'
            lep = '#mu' if '_mu_' in name else 'el' if '_el_' in name else 'l'
            if forceLep: lep = '#mu' if forceLep == "mu" else 'e' if forceLep == "el" else 'l'
            #nn = "{wch}#rightarrow{lep}#nu: |#eta|-p_{{T}} #in [{etal:1.1f},{etah:1.1f}]-[{ptl:3g},{pth:3g}]".format(wch=wch,lep=lep,etal=etal,etah=etah,ptl=ptl,pth=pth)
            nn = "{wch} |#eta|-p_{{T}} #in [{etal:1.1f},{etah:1.1f}]-[{ptl:3g},{pth:3g}]".format(wch=wch,etal=etal,etah=etah,ptl=ptl,pth=pth)
        else:
            nn  = '#mu: ' if '_mu_' in name else 'el: ' if '_el_' in name else ''
            nn += 'W+ ' if 'plus' in name else 'W- ' if 'minus' in name else 'W '
            nn += "i#eta, ip_{{T}} = {neta}, {npt} ".format(neta=ieta,npt=ipt)
        return nn

    elif '_ieta_' in name:
        genEtaPtBins = genBins
        # nn  = '#mu: ' if '_mu_' in name else 'el: ' if '_el_' in name else ''
        # nn += 'W+ ' if 'plus' in name else 'W- ' if 'minus' in name else 'W '
        ieta = int((name.split("_ieta_")[1]).split("_")[0])
        # nn += "i#eta = {neta}".format(neta=ieta)
        nn = ""
        if drawRangeInLabel:
            etal,etah = 0.0,0.0
            if genEtaPtBins:
                etal = genEtaPtBins.etaBins[ieta]
                etah = genEtaPtBins.etaBins[ieta+1]
            wch = 'W+ ' if 'plus' in name else 'W- ' if 'minus' in name else 'W '
            lep = '#mu' if '_mu_' in name else 'el' if '_el_' in name else 'l'
            if forceLep: lep = '#mu' if forceLep == "mu" else 'e' if forceLep == "el" else 'l'
            nn = "{wch} |#eta^{{{lep}}}| #in [{etal:1.1f},{etah:1.1f}]".format(wch=wch,lep=lep,etal=etal,etah=etah)
        else:
            nn  = '#mu: ' if '_mu_' in name else 'el: ' if '_el_' in name else ''
            nn += 'W+ ' if 'plus' in name else 'W- ' if 'minus' in name else 'W '            
            nn += "i#eta = {neta}".format(neta=ieta)
            
        return nn

    elif '_ipt_' in name:
        genEtaPtBins = genBins
        # nn  = '#mu: ' if '_mu_' in name else 'el: ' if '_el_' in name else ''
        # nn += 'W+ ' if 'plus' in name else 'W- ' if 'minus' in name else 'W '
        ipt = int((name.split("_ipt_")[1]).split("_")[0])
        # nn += "ip_{{T}} = {npt}".format(npt=ipt)
        nn = ""
        if drawRangeInLabel:
            ptl,pth = 0.0,0.0
            if genEtaPtBins:
                ptl = genEtaPtBins.ptBins[ipt]
                pth = genEtaPtBins.ptBins[ipt+1]
            wch = 'W+ ' if 'plus' in name else 'W- ' if 'minus' in name else 'W '
            lep = '#mu' if '_mu_' in name else 'el' if '_el_' in name else 'l'
            if forceLep: lep = '#mu' if forceLep == "mu" else 'e' if forceLep == "el" else 'l'
            nn = "{wch} p_{{T}}^{{{lep}}} #in [{ptl:3g},{pth:3g}]".format(wch=wch,lep=lep,ptl=ptl,pth=pth)
        else:
            nn  = '#mu: ' if '_mu_' in name else 'el: ' if '_el_' in name else ''
            nn += 'W+ ' if 'plus' in name else 'W- ' if 'minus' in name else 'W '
            nn += "ip_{{T}} = {npt}".format(npt=ipt)
        return nn

    elif "CMS_" in name:
        # keep Wmu or We now that we do combination, they are different sources
        # if "CMS_Wmu" in name:
        #     return name.replace("CMS_Wmu_","")
        # elif "CMS_We" in name:
        #     return name.replace("CMS_We_","")        
        #else:
        #    return name
        #return name.replace("CMS_","")
        return name

    elif re.match( "Fakes(Eta|EtaCharge|PtNorm|PtSlope)Uncorrelated.*",name):
        num = re.findall(r'\d+', name) # get number
        pfx = name.split(num[0])[1]    # split on number and read what's on the right
        leptonCharge = ""
        if len(pfx):
            leptonCharge = "{lep}{chs}".format(lep="#mu" if "mu" in pfx else "e", chs = "+" if "plus" in pfx else "-" if "minus" in pfx else "")
        tmpvar = ""
        if "FakesEtaCharge" in name: tmpvar = "#eta-ch"
        elif "FakesEta" in name: tmpvar = "#eta"
        elif "FakesPtNorm" in name: tmpvar = "p_{T}-norm"
        elif "FakesPtSlope" in name: tmpvar = "p_{T}-shape"    
        return "QCD bkg {var}-{n} {lepCh}".format(var=tmpvar, n=num[0], lepCh=leptonCharge)

    elif re.match(".*EffStat\d+.*",name):
        num = re.findall(r'\d+', name) # get number (there will be two of them, need the second)
        pfx = name.split("EffStat"+str(num[1]))[1]    # split on second number and read what's on the right
        leptonCharge = ""
        if len(pfx):
            leptonCharge = "{lep}{chs}".format(lep="#mu" if "mu" in pfx else "e", chs = "+" if "plus" in pfx else "-" if "minus" in pfx else "")
        return "Eff.stat. {n1}-{n2} {lepCh}".format(n1=num[0],n2=num[1],lepCh=leptonCharge)


    else:  
        return name

def niceNameHEPDATA(name,genBins="",forceLep="",drawRangeInLabel=False, isHelicity=False):

    # try a simple but understandable form
    # use plain latex and not tlatex for math symbols: i.e. $\mu$ instead of #mu

    if re.match("W.*sumxsec",name) and all(x not in name for x in ["_Ybin", "_ieta_", "_ipt_"]):
        # inclusive xsec poi
        nn  = '\mu' if '_mu_' in name else 'el' if '_el_' in name else 'l'
        if forceLep: nn = '\mu' if forceLep == "mu" else 'e' if forceLep == "el" else 'l'
        chs = "+" if "plus" in name else "-" if "minus" in name else ""
        if chs == "":
            return "$W\\rightarrow {l}\\nu$ fiducial".format(l=nn)
        else:
            return "$W^{{{ch}}}\\rightarrow {l}\\nu$ fiducial".format(ch=chs,l=nn)

    if re.match("Wasym.*chargemetaasym",name) and all(x not in name for x in ["_Ybin", "_ieta_", "_ipt_"]):
        # inclusive xsec poi
        nn  = '\mu' if '_mu_' in name else 'el' if '_el_' in name else 'l'
        if forceLep: nn = '\mu' if forceLep == "mu" else 'e' if forceLep == "el" else 'l'
        return "$W\\rightarrow {l}\\nu$ (asymmetry) fiducial".format(l=nn)

    if re.match("Wratio.*ratiometaratio",name) and all(x not in name for x in ["_Ybin", "_ieta_", "_ipt_"]):
        # inclusive xsec poi
        #nn  = '\mu' if '_mu_' in name else 'el' if '_el_' in name else 'l'
        #if forceLep: nn = '\mu' if forceLep == "mu" else 'e' if forceLep == "el" else 'l'
        return "$W^{+}$/$W^{-}$ fiducial"


    if '_Ybin' in name:
        genYwBins = genBins
        if drawRangeInLabel:
            charge = "plus" if "plus" in name else "minus" if "minus" in name else ""
            wch = "W^{+}" if "plus" in name else "W^{-}" if "minus" in name else "W"
            lep = '\mu' if '_mu_' in name else 'el' if '_el_' in name else 'l'
            if forceLep: lep = '\mu' if forceLep == "mu" else 'e' if forceLep == "el" else 'l'
            pol = "left" if 'left' in name else "right" if "right" in name else "long" if "long" in name else "unpolarized"
            pol_sub = "L" if 'left' in name else "R" if "right" in name else "0" if "long" in name else ""
            ywl,ywh = 0.0,0.0
            if genYwBins:
                iyw = int((name.split("_Ybin_")[1]).split("_")[0])
                # get key to read binning (same for all charges and polarizations, but let's stay general
                # in case charge was not set, take key for "plus" charge and cross fingers
                # similarly, use left for binning when pol is unpolarized
                charge_pol = "{ch}_{pol}".format(ch=charge if charge != "" else "plus",pol=pol if pol!="unpolarized" else "left")  
                ywl = genYwBins[charge_pol][iyw]
                ywh = genYwBins[charge_pol][iyw+1]
                #nn = "{wch}#rightarrow{lep}#nu {pol}: |y_{{W}}| #in [{ywl:1.2f},{ywh:1.2f}]".format(wch=wch,lep=lep,pol=pol,ywl=ywl,ywh=ywh)
            nn = "${wch}_{{{pol}}}$, $|y_{{W}}| \in$ [{ywl:1.2f},{ywh:1.2f}]".format(wch=wch,pol=pol_sub,ywl=ywl,ywh=ywh)
            
            if name.startswith("norm_"):  # normalization systematics for bins not fitted as pois
                nn = "norm.syst. " + nn
        else:
            nn  = '#mu: ' if '_mu_' in name else 'el: ' if '_el_' in name else ''
            if 'plus' in name: nn += 'W+ '
            elif 'minus' in name: nn += 'W- '
            else: nn += 'W '
            if 'left' in name: nn += 'left '
            elif 'right' in name: nn += 'right '
            elif 'long' in name: nn += 'long '
            else: nn += 'unpolarized '
            idx = -2 if (name.endswith('mu') or any([x in name for x in ['masked','sumxsec','charge','a0','a4']])) else -1
            nn += name.split('_')[idx]
            if 'eff_unc' in name:
                nn = '#epsilon_{unc}^{'+nn+'}'
        return nn

    elif '_ieta_' in name and '_ipt_' in name:
        genEtaPtBins = genBins
        ieta,ipt = get_ieta_ipt_from_process_name(name)
        if drawRangeInLabel:
            etal,etah = 0.0,0.0
            ptl,pth = 0.0,0.0
            if genEtaPtBins:
                etal = genEtaPtBins.etaBins[ieta]
                etah = genEtaPtBins.etaBins[ieta+1]
                ptl = genEtaPtBins.ptBins[ipt]
                pth = genEtaPtBins.ptBins[ipt+1]
            wch = 'W^{+}' if 'plus' in name else 'W^{-}' if 'minus' in name else 'W'
            lep = '\mu' if '_mu_' in name else 'el' if '_el_' in name else 'l'
            if forceLep: lep = '\mu' if forceLep == "mu" else 'e' if forceLep == "el" else 'l'
            #nn = "{wch}#rightarrow{lep}#nu: |#eta|-p_{{T}} #in [{etal:1.1f},{etah:1.1f}]-[{ptl:3g},{pth:3g}]".format(wch=wch,lep=lep,etal=etal,etah=etah,ptl=ptl,pth=pth)
            nn = "${wch}$, $|\eta|$-$p_{{T}} \in$ [{etal:1.1f},{etah:1.1f}]-[{ptl:3g},{pth:3g}]".format(wch=wch,etal=etal,etah=etah,ptl=ptl,pth=pth)
        else:
            nn  = '#mu: ' if '_mu_' in name else 'el: ' if '_el_' in name else ''
            nn += 'W+ ' if 'plus' in name else 'W- ' if 'minus' in name else 'W '
            nn += "i#eta, ip_{{T}} = {neta}, {npt} ".format(neta=ieta,npt=ipt)
        return nn

    elif '_ieta_' in name:
        genEtaPtBins = genBins
        # nn  = '#mu: ' if '_mu_' in name else 'el: ' if '_el_' in name else ''
        # nn += 'W+ ' if 'plus' in name else 'W- ' if 'minus' in name else 'W '
        ieta = int((name.split("_ieta_")[1]).split("_")[0])
        # nn += "i#eta = {neta}".format(neta=ieta)
        nn = ""
        if drawRangeInLabel:
            etal,etah = 0.0,0.0
            if genEtaPtBins:
                etal = genEtaPtBins.etaBins[ieta]
                etah = genEtaPtBins.etaBins[ieta+1]
            wch = 'W^{+}' if 'plus' in name else 'W^{-}' if 'minus' in name else 'W'
            lep = '\mu' if '_mu_' in name else 'el' if '_el_' in name else 'l'
            if forceLep: lep = '\mu' if forceLep == "mu" else 'e' if forceLep == "el" else 'l'
            nn = "${wch}$, $|\eta^{{{lep}}}| \in$ [{etal:1.1f},{etah:1.1f}]".format(wch=wch,lep=lep,etal=etal,etah=etah)
        else:
            nn  = '#mu: ' if '_mu_' in name else 'el: ' if '_el_' in name else ''
            nn += 'W+ ' if 'plus' in name else 'W- ' if 'minus' in name else 'W '            
            nn += "i#eta = {neta}".format(neta=ieta)
            
        return nn

    elif '_ipt_' in name:
        genEtaPtBins = genBins
        # nn  = '#mu: ' if '_mu_' in name else 'el: ' if '_el_' in name else ''
        # nn += 'W+ ' if 'plus' in name else 'W- ' if 'minus' in name else 'W '
        ipt = int((name.split("_ipt_")[1]).split("_")[0])
        # nn += "ip_{{T}} = {npt}".format(npt=ipt)
        nn = ""
        if drawRangeInLabel:
            ptl,pth = 0.0,0.0
            if genEtaPtBins:
                ptl = genEtaPtBins.ptBins[ipt]
                pth = genEtaPtBins.ptBins[ipt+1]
            wch = 'W^{+}' if 'plus' in name else 'W^{-}' if 'minus' in name else 'W'
            lep = '\mu' if '_mu_' in name else 'el' if '_el_' in name else 'l'
            if forceLep: lep = '\mu' if forceLep == "mu" else 'e' if forceLep == "el" else 'l'
            nn = "${wch}$, $p_{{T}}^{{{lep}}} \in$ [{ptl:3g},{pth:3g}]".format(wch=wch,lep=lep,ptl=ptl,pth=pth)
        else:
            nn  = '#mu: ' if '_mu_' in name else 'el: ' if '_el_' in name else ''
            nn += 'W+ ' if 'plus' in name else 'W- ' if 'minus' in name else 'W '
            nn += "ip_{{T}} = {npt}".format(npt=ipt)
        return nn

    elif re.match( "Fakes(Eta|EtaCharge|PtNorm|PtSlope)Uncorrelated.*",name):

        num = re.findall(r'\d+', name) # get number
        pfx = name.split(num[0])[1]    # split on number and read what's on the right
        leptonCharge = ""
        if len(pfx):
            if "plus" in pfx or "minus" in pfx:
                leptonCharge = "{lep}^{chs}".format(lep="\mu" if "mu" in pfx else "e", chs = "+" if "plus" in pfx else "-")
            else:
                leptonCharge = "{lep}".format(lep="\mu" if "mu" in pfx else "e")
        tmpvar = ""
        FakesBins = []
        pt_or_eta = ""
        # eta binning for fakes systematics (PtNorm uses another one chosen below)
        if isHelicity:
            if "mu" in pfx:
                FakesBins = [-2.4] + [round(-2.0+0.5*i,1) for i in range(9)] + [2.4] 
            else:
                FakesBins = [-2.5] + [round(-2.4+0.2*i,1) for i in range(25)] + [2.5]
        else:
            FakesBins = [-2.4,-2.1,-1.5,-1.0,-0.5,0.0,0.5,1.0,1.5,2.1,2.4]

        if "FakesEtaCharge" in name: 
            tmpvar = "$\eta$-norm-chargeUncorr"
            pt_or_eta = "\eta"
        elif "FakesEta" in name: 
            tmpvar = "$\eta$-norm"
            pt_or_eta = "\eta"
        elif "FakesPtSlope" in name: 
            tmpvar = "$p_{T}$-shape"    
            pt_or_eta = "\eta"
        elif "FakesPtNorm" in name: 
            tmpvar = "$p_{T}$-norm"
            pt_or_eta = "p_{T}"
            # for pt norm has to bin on pt
            if isHelicity:
                FakesBins = [26, 29, 32, 35, 38, 41, 45] if "mu" in pfx else [30, 35, 40, 45]
            else:
                FakesBins = [26, 33, 36, 40.5, 45, 50, 56] if "mu" in pfx else [30, 36, 40.5, 45, 50, 56]
            
        #print "{n} FakeBins[{m}]".format(n=name,m=num)
        n = int(num[0])
        vlow = FakesBins[n-1]
        vhigh = FakesBins[n]
        return "QCD bkg {var}, ${l}<{v}<{h}$, ${lepCh}$".format(var=tmpvar, l=vlow, v=pt_or_eta, h=vhigh, lepCh=leptonCharge)

    elif re.match(".*EffStat\d+.*",name):

        ## utility arrays for the effstat binning
        _etaBinsEffStatMu = [round(-2.4 + 0.1* x,1) for x in range(0,49)]
        _etaBinsEffStatEl_diffXsec = [round(-2.5 + 0.1* x,1) for x in range(0,51)] # then first and last bins were not used, so the index in the nuis name starts from 2 and ends at 49
        _etaBinsEffStatEl_helicity = [-2.5,-2.4,-2.3,-2.2,-2.1,-2.0,-1.9,-1.8,-1.7,-1.6,-1.566,-1.4442,-1.4,-1.3,-1.2,-1.1,-1.0,-0.9,-0.8,-0.7,-0.6,-0.5,-0.4,-0.3,-0.2,-0.1,0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.4442,1.566,1.6,1.7,1.8,1.9,2.0,2.1,2.2,2.3,2.4,2.5] # here the index in nuis name goes from 1 to 52, because the nuisances were defined based on the reco-eta binning rather than on the actual systs, so the gap between EB and EE induces two more bins than it should have beeen normally

        num = re.findall(r'\d+', name) # get number (there will be two of them, need the second)
        pfx = name.split("EffStat"+str(num[1]))[1]    # split on second number and read what's on the right
        leptonCharge = ""
        if len(pfx):
            leptonCharge = "{lep}^{chs}".format(lep="\mu" if "mu" in pfx else "e", chs = "+" if "plus" in pfx else "-" if "minus" in pfx else "")
        n2 = int(num[1])
        etaBinsEffStat = _etaBinsEffStatMu if "mu" in pfx else _etaBinsEffStatEl_helicity if isHelicity else _etaBinsEffStatEl_diffXsec
        etalow  = str(etaBinsEffStat[n2-1])
        etahigh = str(etaBinsEffStat[n2])
        return "Eff.stat. param.{n1}, ${l}<\eta<{h}$, ${lepCh}$".format(n1=num[0],l=etalow,h=etahigh,lepCh=leptonCharge)

    elif "L1PrefireEleEffSyst" in name:
        num = re.findall(r'\d+', name) # get number
        #print "{n} L1PrefireEleEffSys[{m}]".format(n=name,m=num)
        n = int(num[1]) # use second number, the first is '1' in 'L1Prefire...'
        etabinsPrefire = []
        if isHelicity:
            etabinsPrefire = [-2.5,  -2.35, -2.2, -2.0, -1.5, 1.5, 2.0, 2.2, 2.35, 2.5]
            # n will go from 0 to 7 (skipping bin in barrel, so only 8 bins)
            if n > 3: n = n + 1 # have to skip bin in [-1.5, 1.5]
        else:
            etabinsPrefire = [-2.4, -2.1, -1.9, -1.5, 1.5, 1.9, 2.1, 2.4]
            # n will go from 0 to 5 (skipping bin in barrel, so only 6 bins)
            if n > 2: n = n + 1 # have to skip bin in [-1.5, 1.5]
        low = etabinsPrefire[n]
        high = etabinsPrefire[n+1]
        return "L1-trigger electron eff.syst., ${l}<\eta<{h}$".format(l=low,h=high)

    elif "OutOfAccPrefireSyst" in name:
        sign = "<" if "0" in name else ">"
        return "L1-trigger electron eff.syst. for Drell-Yan, $\eta{sign}0$".format(sign=sign)

    elif re.match( "smooth(el|mu)scale.*",name):
        num = re.findall(r'\d+', name) # get number 
        n = 0
        n2 = 0
        n3 = 0
        lep = "\mu" if "smoothmu" in name else "e"
        if "scaleStat" in name: 
            n = int(num[0])
            return "$p_{{T}}^{{{lep}}}$ scale stat.{n}".format(lep=lep,n=n)
        else:
            n = int(num[0])
            n2 = int(num[1])
            if isHelicity:
                etabinsPtSyst = [-2.4, -2.1, 0.0, 2.1, 2.4] if "smoothmu" in name else [-2.5, -2.1, -1.5, -1.0, 0.0, 1.0, 1.5, 2.1, 2.5]
                low = str(etabinsPtSyst[n2])
                high = str(etabinsPtSyst[n2+1])
                ptText = ""
                if lep == "e": # electron syst is divided in two pt bins at 42 GeV
                    n2 = int(num[2])
                    ptText = ", $p_{{T}}^{{{lep}}}{sign}42$".format(lep=lep,sign="<" if n2==0 else ">")
                return "$p_{{T}}^{{{lep}}}$ scale syst.{n}, ${l}<\eta<{h}${ch}{t}".format(lep=lep,n=n,l=low,h=high,ch=", charge +" if "plus" in name else ", charge -" if "minus" in name else "",t=ptText)
            else:
                etabinsPtSyst = [0.0, 2.1, 2.4] if "smoothmu" in name else [0.0, 1.0, 1.5, 2.1, 2.4]
                # match the 'P' to select positive eta side
                if re.match(".*etaside\d+P(plus|minus)*",name): 
                    low = str(etabinsPtSyst[n2])
                    high = str(etabinsPtSyst[n2+1])
                else:
                    low = "-"+str(etabinsPtSyst[n2+1])
                    high = "-"+str(etabinsPtSyst[n2])
                return "$p_{{T}}^{{{lep}}}$ scale syst.{n}, ${l}<\eta<{h}${ch}".format(lep=lep,n=n,l=low,h=high,ch=", charge +" if "plus" in name else ", charge -" if "minus" in name else "")

    elif "TnPEffSyst" in name or "TestEffSyst" in name:
        num = re.findall(r'\d+', name) # get number
        #print "{n} EffSyst[{m}]".format(n=name,m=num)
        n = int(num[0])
        if isHelicity:
            etabinsEffSyst = [0,1.0,1.5,2.4] if "mu" in name else [0,1.0,1.5,2.0,2.5]
        else:
            etabinsEffSyst = [0,1.0,1.5,2.4] if "mu" in name else [0,1.0,1.5,1.9,2.1,2.4]
        low = etabinsEffSyst[n]
        high = etabinsEffSyst[n+1]
        # use double '\' for abs, because '\a' is a special character
        return "eff.syst., ${l}<|\eta|<{h}$, ${lep}$".format(lep="\mu" if "mu" in name else "e",l=low,h=high)

    elif "fsr" in name:
        return "{lep} QED final state radiation".format(lep="muon" if "Mu" in name else "electron")

    elif name == "mW":
        return "$m_W$"

    elif any(x in name for x in ["muR", "muF", "muRmuF"]):
        scale = "$\mu_{R}\mu_{F}$" if "muRmuF" in name else "$\mu_{R}$" if "muR" in name else "$\mu_{F}$"
        charge = "+" if "plus" in name else "-" if "minus" in name else ""
        pol = "L" if "left" in name else "R" if "right" in name else "0" if "long" in name else ""
        if any(x == name for x in ["muR", "muF", "muRmuF"]):
            return "{s} (Drell-Yan bkg)".format(s=scale)
        else:
            num = re.findall(r'\d+', name) # get number
            n = int(num[0]) # goes from 1 to 10
            ptWbins = [0.0, 2.9, 4.7, 6.7, 9.0, 11.8, 15.0, 20.1, 27.2, 40.2, 13000.0]
            low = ptWbins[n-1]
            high = ptWbins[n]
            if isHelicity:
                if pol != "":
                    wch = "$W^{{{c}}}_{{{p}}}$".format(c=charge,p=pol)
                    txt = "(W signal)"
                else:
                    wch = "$W^{{{c}}}$".format(c=charge)
                    txt = "($W\\rightarrow\\tau\\nu$ bkg)"
            else:
                wch = "$W^{{{c}}}$".format(c=charge)
                txt = "(W signal and $W\\rightarrow\\tau\\nu$ bkg)"
            if n == 10:
                return "{s} {w}, $p_{{T}}^{{W}}>{l}$ {t}".format(w=wch,s=scale,l=low,t=txt)
            else:
                return "{s} {w}, ${l}<p_{{T}}^{{W}}<{h}$ {t}".format(w=wch,s=scale,l=low,h=high,t=txt)

    elif "alpha" in name:
        return "$\\alpha_{S}$"

    elif "pdf" in name:
        num = re.findall(r'\d+', name) # get number                                               
        n = int(num[0]) # goes from 1 to 10            return "$\alpha_{S}$"
        return "Hessian {i}".format(i=n)

    elif name.startswith("CMS_"):
        if "lumi" in name:
            return "luminosity"
        elif "Tau" in name:
            return "$W\\rightarrow\\tau\\nu$ bkg norm."
        elif "VV" in name:
            return "Diboson bkg norm."
        elif "Top" in name:
            return "t quark bkg norm."
        elif "flips" in name:
            return "Charge flips bkg norm."
        elif "bkg_lepeff" in name:
            return "eff.syst. bkg, ${l}$".format(l="e" if "We" in name else "\mu")
        elif "lepVeto" in name:
            return "second lepton veto, ${l}$".format(l="e" if "We" in name else "\mu")
        else:
            return name

    else:  
        return name
        

if __name__ == "__main__":

    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetPaintTextFormat('.3f')

    date = datetime.date.today().isoformat()

    from optparse import OptionParser
    parser = OptionParser(usage='%prog toys.root [options] ')
    parser.add_option('-o','--outdir', dest='outdir',    default='', type='string', help='output directory to save the matrix')
    parser.add_option('-p','--params', dest='params',    default='', type='string', help='parameters for which you want to show the correlation matrix. comma separated list of regexps')
    parser.add_option('-t','--type'  , dest='type'  ,    default='toys', type='string', help='which type of input file: toys(default),scans, or hessian')
    parser.add_option(     '--suffix', dest='suffix',    default='', type='string', help='suffix for the correlation matrix')
    parser.add_option(     '--parNameCanvas', dest='parNameCanvas',    default='', type='string', help='The canvas name is built using the parameters selected with --params. If they are many, better to pass a name, like QCDscales or PDF for example')
    parser.add_option(     '--nContours', dest='nContours',    default=51, type=int, help='Number of contours in palette. Default is 51 (keep it odd: no correlation is white with our palette)')
    parser.add_option(     '--palette'  , dest='palette',      default=0, type=int, help='Set palette: default is a built-in one, 55 is kRainbow')
    parser.add_option(     '--vertical-labels-X', dest='verticalLabelsX',    default=False, action='store_true', help='Set labels on X axis vertically (sometimes they overlap if rotated)')
    parser.add_option(     '--title'  , dest='title',    default='', type='string', help='Title for matrix ')
    parser.add_option(     '--show-more-correlated' , dest='showMoreCorrelated',    default=0, type=int, help='Show the N nuisances more correlated (in absolute value) with the parameters given with --params. If 0, do not do this part')
    parser.add_option('-m','--matrix-type', dest='matrixType',    default='channelpmaskedexpnorm', type='string', help='Select which matrix to read from file')
    parser.add_option(     '--margin',     dest='margin',     default='', type='string', help='Pass canvas margin as "left,right,top,bottom" ')
    parser.add_option(     '--canvasSize', dest='canvasSize', default='', type='string', help='Pass canvas dimensions as "width,height" ')
    parser.add_option(     '--etaptbinfile',   dest='etaptbinfile',   default='',  type='string', help='eta-pt binning used for labels with 2D xsec')
    parser.add_option(     '--ywbinfile',   dest='ywbinfile',   default='',  type='string', help='Yw binning used for labels with helicity')
    parser.add_option('-c','--channel',     dest='channel',     default='', type='string', help='Channel (el|mu|lep), if not given it is inferred from the inputs, but might be wrong depending on naming conventions')
    parser.add_option(     '--show-all-nuisances', dest='showAllNuisances',    default=False, action='store_true', help='Show all nuisances in the matrix (e.g. to prepare HEPdata entries): this implies that option --params is only used to filter POIs')
    parser.add_option('--which-matrix',  dest='whichMatrix',     default='both', type='string', help='Which matrix: covariance|correlation|both')
    parser.add_option(     '--divide-covariance-by-bin-area', dest='divideCovarianceBybinArea',    default=False, action='store_true', help='Divide covariance by bin area (2D xsec) or bin width (1D xsec)')
    parser.add_option(     '--skipLatexOnTop', dest='skipLatexOnTop',    default=False, action='store_true', help='Do not write "CMS blabla" on top (mainly useful when a title is needed)')
    parser.add_option(     '--use-hepdata-labels', dest='useHepdataLabels',    default=False, action='store_true', help='Write axis labels using latex for hepdata, with some name polishing')
    parser.add_option(     '--old-paper-style', dest='oldPaperStyle',    default=False, action='store_true', help='This option is meant to reproduce the plots in the paper, where the style was slightly different than what is implemented now (e.g. for axis labels)')
    (options, args) = parser.parse_args()

    if options.divideCovarianceBybinArea:
        if not options.whichMatrix == "covariance":
            print "Warning: you are not using a covariance matrix, but correlation one"
            print "Then, option --divide-covariance-by-bin-area will be ignored"
            options.divideCovarianceBybinArea = False
        if not options.etaptbinfile and not options.ywbinfile:
            print "Error: option --divide-covariance-by-bin-area require a file with binning."
            print "Please specify one using either --etaptbinfile or --ywbinfile"
            quit()
        
    if options.oldPaperStyle and options.useHepdataLabels:
        print "Error: options --old-paper-style and --use-hepdata-labels are not compatible."
        print "Please select only one of them and try again"
        quit()

    ROOT.TColor.CreateGradientColorTable(3,
                                      array ("d", [0.00, 0.50, 1.00]),
                                      ##array ("d", [1.00, 1.00, 0.00]),
                                      ##array ("d", [0.70, 1.00, 0.34]),
                                      ##array ("d", [0.00, 1.00, 0.82]),
                                      array ("d", [0.00, 1.00, 1.00]),
                                      array ("d", [0.34, 1.00, 0.65]),
                                      array ("d", [0.82, 1.00, 0.00]),
                                      255,  0.95)


    if not options.type in ['toys', 'scans', 'hessian']:
        print 'the given type needs to be either "toys", "scans", or "hessian"!!'
        sys.exit()

    if len(args) < 1:
        print 'You have to pass a root file with the fit result, either for toys or hessian'
        sys.exit()

    if options.outdir:
        ROOT.gROOT.SetBatch()
        if not os.path.isdir(options.outdir):
            os.system('mkdir -p {od}'.format(od=options.outdir))
        os.system('cp {pf} {od}'.format(pf='/afs/cern.ch/user/m/mciprian/public/index.php',od=options.outdir))

    pois_regexps = list(options.params.split(','))
    print "Filtering POIs with the following regex: ",pois_regexps

    if options.etaptbinfile and options.ywbinfile:
        print "Error: options --etaptbinfile and --ywbinfile are incompatible. Choose only one based on the fit"
        quit()

    if options.etaptbinfile:
        print "HERE for 2D xsec"
        etaPtBinningVec = getDiffXsecBinning(options.etaptbinfile, "gen")
        genBins = templateBinning(etaPtBinningVec[0],etaPtBinningVec[1])
    elif options.ywbinfile:
        print "HERE for helicity"
        ybinfile = open(options.ywbinfile, 'r')
        genBins = eval(ybinfile.read())
        ybinfile.close()
    else:
        genBins = ""
    
    print "="*30
    print "Gen binning:"
    print "-"*30
    if options.etaptbinfile:
        genBins.printBinAll() 
    else:
        print genBins
    print "="*30


    params = []; indices = []

    ## directly store the mean and RMS into a dictionary
    fitvals = {}; fiterrs = {}

    cov = {}; corr = {}

    nNuisances = 0
    nPois = 0

    ### GET LIST OF PARAMETERS THAT MATCH THE SPECIFIED OPTION IN THE TOYFILE
    if options.type == 'toys':
        toyfile = ROOT.TFile(args[0], 'read')
        _tree = toyfile.Get('fitresults')
        lol = _tree.GetListOfLeaves()

        for l in lol:
            ## skip a bunch of those we don't want
            if '_err'   in l.GetName(): continue
            if '_minos' in l.GetName(): continue
            if '_gen'   in l.GetName(): continue
            if '_In'    in l.GetName(): continue
            for poi in pois_regexps:
                if re.match(poi, l.GetName()):
                    ## draw the parameter into a histogram
                    _tree.Draw(l.GetName()+'>>h_'+l.GetName())
                    ## find that histogram and clone it
                    h = ROOT.gROOT.FindObject('h_'+l.GetName()).Clone()
                    ## store mean and rms into the dictionaries from before
                    fitvals[l.GetName()] = h.GetMean()
                    fiterrs[l.GetName()] = h.GetRMS()
                    ## also keep a list of the parameter names, for sorting
                    params.append(l.GetName())

    elif options.type == 'hessian':
        hessfile = ROOT.TFile(args[0],'read')
        suffix = options.matrixType
        corrmatrix = hessfile.Get('correlation_matrix_'+suffix)
        covmatrix  = hessfile.Get('covariance_matrix_'+suffix)
        for ib in range(1,corrmatrix.GetNbinsX()+1):
            #if nNuisances > 5 and nPois > 10: break # only for tests
            for poi in pois_regexps:
                if re.match(poi, corrmatrix.GetXaxis().GetBinLabel(ib)):
                    ## store mean and rms into the dictionaries from before
                    ## also keep a list of the parameter names, for sorting
                    params .append(corrmatrix.GetXaxis().GetBinLabel(ib))
                    indices.append(ib)
                    nPois += 1
                elif options.showAllNuisances:
                    # exclude pois, hopefully no nuisance has the name starting with W!
                    if not corrmatrix.GetXaxis().GetBinLabel(ib).startswith("W"):
                        params .append(corrmatrix.GetXaxis().GetBinLabel(ib))
                        indices.append(ib)
                        nNuisances += 1

    if options.showAllNuisances:
        print "nPois = %d" % nPois
        print "nNuisances = %d" % nNuisances

    filter_matrixType_poiPostfix = {"channelpmaskedexp"     : "pmaskedexp",
                                    "channelpmaskedexpnorm" : "pmaskedexpnorm",
                                    "channelsumpois"        : "sumxsec",
                                    "channelsumpoisnorm"    : "sumxsecnorm",
                                    "channelchargepois"     : "chargeasym",
                                    "channelchargemetapois" : "chargemetaasym",
                                    "channelratiometapois"  : "ratiometaratio",
                                    "channelpolpois"        : "a4",  #  check if it is correct
                                    #"channelpolpois"        : "unpolarizedxsec", 
                                    "channelnone"           : "pmaskedexp" # dummy, there are no POIs in this case
    }
    poiPostfix = filter_matrixType_poiPostfix[options.matrixType]
    poiHasToBeScaled = True if poiPostfix in ["pmaskedexp", "sumxsec"] else False
    ## construct the covariances and the correlations in one go.
    ## for the covariance, need to scale xsec content to port it from number of events to pb (and for combination need to divide by 2), but only for absolute things
    # also normalize by bin area the  numbers for absolute and normalized cross section (as it is done on the plots)
    # the yields/pb scale factor is 35900 for 2D xsec, and 36000 for helicity
    # this is because the yields for helicity were made using 36/fb
    scalefactor_poi = 1./(36000.0 if options.ywbinfile else 35900.0)
    if options.channel == "lep":
        scalefactor_poi /= 2.            

    #######
    # NOTE 
    #######
    # should change signs for A4 in the W+, because the fit returns a positive number, but theory
    # predicts A4<0 for W+, so we manually change sign in the plot
    # this must also reflect in the covariance matrix
    # i.e. if the covariance with a nuisance and a yW bin for A4 in W+ channel is negative, it means scaling the former up moves the latter down, but since A4 should have opposite sign the direction of the shift is swapped as well


    print "Preparing list of parameters to build the matrix ..."

    for ip1, p1 in enumerate(params):
        for ip2, p2 in enumerate(params):
            if options.type == 'toys':
                var = '({x}-{x0})*({y}-{y0})'.format(x=p1,x0=fitvals[p1],y=p2,y0=fitvals[p2])
                _tree.Draw('{var}>>h_{x}_{y}'.format(var=var,x=p1,y=p2))
                h = ROOT.gROOT.FindObject('h_{x}_{y}'.format(x=p1,y=p2)).Clone()
                cov [(p1,p2)] = h.GetMean()
                corr[(p1,p2)] = cov[(p1,p2)]/(fiterrs[p1]*fiterrs[p2])
                # this part is obsolete and might not have all the proper scaling factors used for hessian below
            elif options.type == 'hessian':
                scalefactor = 1.0
                # normalization by bin area/width only for 2D/1D xsec, but might be added for rapidity as well
                scalefactorSignA4wplus = 1.0
                if options.divideCovarianceBybinArea: 
                    if any(x in p1 for x in ["_ieta_","_ipt_"]):
                        scalefactor /= getBinAreaFromParamName(p1,genBins)
                    if any(x in p2 for x in ["_ieta_","_ipt_"]):
                        scalefactor /= getBinAreaFromParamName(p2,genBins)
                    if "_Ybin_" in p1:
                        scalefactor /= getBinAreaFromParamName(p1,genBins,isHelicity=True)
                    if "_Ybin_" in p2:
                        scalefactor /= getBinAreaFromParamName(p2,genBins,isHelicity=True)
                if poiHasToBeScaled:
                    if p1.endswith(poiPostfix):     
                        scalefactor *= scalefactor_poi
                        if poiPostfix == "a4" and "Wplus_Ybin_" in p1:
                            scalefactorSignA4wplus *= -1.0 # see previous comment on A4 for W+
                if p2.endswith(poiPostfix):     
                        scalefactor *= scalefactor_poi
                        if poiPostfix == "a4" and "Wplus_Ybin_" in p2:
                            scalefactorSignA4wplus *= -1.0 # see previous comment on A4 for W+
                cov [(p1,p2)] = scalefactorSignA4wplus * scalefactor * covmatrix .GetBinContent(indices[ip1],indices[ip2])
                corr[(p1,p2)] = scalefactorSignA4wplus * corrmatrix.GetBinContent(indices[ip1],indices[ip2])
        

    print "===> Build covariance matrix from this set of params: ", params

    p_tmp = set(params)
    params = list(p_tmp)

    # to help sorting with helicity
    # if using more helicity and Y bins, sort by hel,Ybin
    helSorted = { "left" : 1, "right" : 2, "long" : 3}
    chargeSorted = { "Wplus" : 1, "Wminus" : 2}
    lepSorted = { "mu" : 1, "el" : 2}

    ## sort the floatParams. alphabetically, except for pdfs, which are sorted by number
    ## for mu* QCD scales, distinguish among muR and muRXX with XX in 1-10
    
    # why is this commented? Isn't it nice that different charge and polarizations are grouped together?
    # see old example here: 
    # http://mciprian.web.cern.ch/mciprian/wmass/13TeV/helicityAnalysis/electron/fromEmanuele/13-12-18/fitresults_poim1_exp1_bbb1/subMatrix/smallCorrelation_Wplusminus_leftrightYbin_2-5.png
    #params = sorted(params, key= lambda x: (int(chargeSorted[x.split('_')[0]]),int(helSorted[x.split('_')[1]]),int(x.split('_')[-1])) if '_Ybin_' in x else 0)
    # one might want to invert the order of charge and helicity for the sorting

    params = sorted(params, key = lambda x: 0 if "plus" in x else 1 if "minus" in x else 2)
    if options.ywbinfile:
        params = sorted(params, key= lambda x: utilities.getNFromString(x) if '_Ybin_' in x else 0)
        params = sorted(params, key= lambda x: (1 if "left" in x else 2 if "right" in x else "3") if '_Ybin_' in x else 0)
        params = sorted(params, key= lambda x: (int(chargeSorted[x.split('_')[0]]),int(helSorted[x.split('_')[1]]),int(x.split('_')[-2])) if '_Ybin_' in x else 0)
    else:
        params = sorted(params, key= lambda x: get_ieta_from_process_name(x) if ('_ieta_' in x) else 0)
        params = sorted(params, key= lambda x: get_ipt_from_process_name(x) if ('_ipt_' in x) else 0)
        params = sorted(params, key= lambda x: get_ieta_ipt_from_process_name(x) if ('_ieta_' in x and '_ipt_' in x) else 0)
    params = sorted(params, key= lambda x: 1 if x.endswith(poiPostfix) else 0)
 
    # sort if not using all params (otherwise the order should be already ok, as it is taken from the original matrix)
    if not options.showAllNuisances:

        params = sorted(params, key= lambda x: int(x.replace('pdf','')) if 'pdf' in x else 100 if 'alphaS' in x else 0)
        params = sorted(params, key= lambda x: int(re.sub('\D','',x)) if ('muRmuF' in x and x != "muRmuF")  else 0)
        params = sorted(params, key= lambda x: int(re.sub('\D','',x)) if (''.join([j for j in x if not j.isdigit()]) == 'muR' and x != "muR") else 0)
        params = sorted(params, key= lambda x: int(re.sub('\D','',x)) if (''.join([j for j in x if not j.isdigit()]) == 'muF' and x != "muF") else 0)
        params = sorted(params, key= lambda x: utilities.getNFromString(x) if 'EffStat' in x else 0)            
        params = sorted(params, key= lambda x: lepInFakeSystForSort(x) if 'Fakes' in x else 0)   
        params = sorted(params, key= lambda x: utilities.getNFromString(x) if 'FakesEtaUncorrelated' in x else 0)    
        params = sorted(params, key= lambda x: utilities.getNFromString(x) if 'FakesPtNormUncorrelated' in x else 0)
        params = sorted(params, key= lambda x: utilities.getNFromString(x) if 'FakesEtaChargeUncorrelated' in x else 0) 
        params = sorted(params, key= lambda x: utilities.getNFromString(x) if 'FakesPtSlopeUncorrelated' in x else 0)   
        # sort by charge if needed     
        # I think it is useful that different charges are separated in the matrix, it is easier to read it 
        #params = sorted(params, key = lambda x: 0 if "right" in x else -1 if "left" in x else -1)

    print "sorted params = ", params

    print "="*30
    print "Now going to create the matrix"
    print "="*30
    ch = 1200
    cw = 1200
    if options.canvasSize:
        cw = int(options.canvasSize.split(',')[0])
        ch = int(options.canvasSize.split(',')[1])
    c = ROOT.TCanvas("c","",cw,ch,)
    c.SetGridx()
    c.SetGridy()
    #ROOT.gStyle.SetPalette(55)
    #ROOT.gStyle.SetNumberContours(200); # default is 20 (values on palette go from -1 to 1)
    if options.nContours: ROOT.gStyle.SetNumberContours(options.nContours)
    if options.palette:   ROOT.gStyle.SetPalette(options.palette)


    clm = 0.15
    crm = 0.15
    cbm = 0.15
    ctm = 0.07
    if options.margin:
        clm,crm,ctm,cbm = (float(x) for x in options.margin.split(','))
    c.SetLeftMargin(clm)
    c.SetRightMargin(crm)
    c.SetBottomMargin(cbm)
    c.SetTopMargin(ctm)

    ## make the new, smaller TH2D correlation matrix
    nbins = len(params)
    th2_sub = ROOT.TH2D('sub_corr_matrix', 'correlation matrix', nbins, 0., nbins, nbins, 0., nbins)
    th2_cov = ROOT.TH2D('sub_cov_matrix',  'covariance matrix', nbins, 0., nbins, nbins, 0., nbins)

    if 'Wplus' in options.params and 'pmaskedexpnorm' in options.params:
        pass
        #th2_sub.SetTitle('correlations of W^{+} processes')
    if 'Wminus' in options.params and 'pmaskedexpnorm' in options.params:
        pass
        #th2_sub.SetTitle('correlations of W^{-} processes')
    if 'pdf' in options.params:
        #th2_sub.SetTitle('correlations of PDF nuisance parameters')
        th2_sub.GetXaxis().SetLabelSize(0.025)
        th2_sub.GetYaxis().SetLabelSize(0.025)
        #th2_cov.SetTitle('covariance of PDF nuisance parameters')
        th2_cov.GetXaxis().SetLabelSize(0.025)
        th2_cov.GetYaxis().SetLabelSize(0.025)

    th2_sub.GetXaxis().SetTickLength(0.)
    th2_sub.GetYaxis().SetTickLength(0.)
    th2_cov.GetXaxis().SetTickLength(0.)
    th2_cov.GetYaxis().SetTickLength(0.)
    
    ## pretty nested loop. enumerate the tuples
    nParams = len(params)
    # set axis labels
    print "Setting Labels"
    for i,x in enumerate(params):
        sys.stdout.write('Row {num}/{tot}   \r'.format(num=i,tot=nParams))
        sys.stdout.flush()
        if options.useHepdataLabels:
            new_x = niceNameHEPDATA(x,genBins=genBins,forceLep=options.channel,drawRangeInLabel=True,isHelicity=True if options.ywbinfile else False)
        else:
            new_x = niceName(x,genBins=genBins,forceLep=options.channel,drawRangeInLabel=False if options.oldPaperStyle else True)
        th2_sub.GetXaxis().SetBinLabel(i+1, new_x)
        th2_sub.GetYaxis().SetBinLabel(i+1, new_x)
        th2_cov.GetXaxis().SetBinLabel(i+1, new_x)
        th2_cov.GetYaxis().SetBinLabel(i+1, new_x)
         
    print "Setting Values"
    for i,x in enumerate(params):
        for j,y in enumerate(params):
            if j>i: break
            sys.stdout.write('Row {num}/{tot}   Column {col}/{tot}   \r'.format(num=i,tot=nParams, col=j))
            sys.stdout.flush()
            ## note that the matrix is symmetric
            if not options.whichMatrix == "covariance":
                th2_sub.SetBinContent(i+1, j+1, corr[(x,y)])
                th2_sub.SetBinContent(j+1, i+1, corr[(x,y)])
            if not options.whichMatrix == "correlation":
                th2_cov.SetBinContent(i+1, j+1, cov [(x,y)])
                th2_cov.SetBinContent(j+1, i+1, cov [(x,y)])

    th2_sub.GetZaxis().SetRangeUser(-1, 1)
    
    covMax = max(abs(th2_cov.GetMaximum()), abs(th2_cov.GetMinimum()))
    th2_cov.GetZaxis().SetRangeUser(-1.*covMax, covMax)

    print "="*30
    print "Now finally drawing the matrix"
    print "="*30
    matricesToPlot = []
    if options.whichMatrix == "both": 
        matricesToPlot = [th2_sub, th2_cov]
    elif options.whichMatrix == "covariance": 
        matricesToPlot = [th2_cov]
    else:
        matricesToPlot = [th2_sub]

    for im,tmp_mat in enumerate(matricesToPlot):

        if options.whichMatrix == "both":
            corcov = 'Correlation' if not im else 'Covariance'
        else:
            corcov = 'Covariance' if options.whichMatrix == "covariance" else "Correlation"

        tmp_mat.SetTitle("")
        if options.title: 
            tmp_mat.SetTitle(options.title)

        if options.outdir:
            ROOT.gStyle.SetPaintTextFormat('1.2f')
            if len(params)<30: tmp_mat.Draw('colz text45')
            else: tmp_mat.Draw('colz')

            lat = ROOT.TLatex()
            lat.SetNDC(); lat.SetTextFont(42)
            if not options.skipLatexOnTop:
                offsetLatex = c.GetLeftMargin()-0.15
                rightTextOffset = 0.56 if options.oldPaperStyle else 0.58
                lat.DrawLatex(0.15+offsetLatex, 0.95, '#bf{CMS}') #it{Preliminary}')
                lat.DrawLatex(rightTextOffset+offsetLatex, 0.95, '35.9 fb^{-1} (13 TeV)')

            if options.verticalLabelsX: tmp_mat.LabelsOption("v","X")
            if nbins >= 20: tmp_mat.LabelsOption("v","X")

            if options.parNameCanvas: 
                paramsName = options.parNameCanvas
            else : 
                paramsName = options.params.replace(',','AND')
                for x in ['.', '*', '$', '^', '|', '[', ']', '(', ')']:
                    paramsName = paramsName.replace(x,'')

            suff = '' if not options.suffix else '_'+options.suffix
            outfname = options.outdir+'/small{corcov}{suff}_{pn}'.format(suff=suff,pn=paramsName,corcov=corcov)
            for i in ['pdf', 'png']:
                c.SaveAs('{ofn}.{i}'.format(ofn=outfname,i=i))
            # save matrix in root file
            matRootFile = ROOT.TFile.Open("{ofn}.root".format(ofn=outfname),"recreate")
            matRootFile.cd()
            tmp_mat.Write()
            matRootFile.Close("matrix{corcov}".format(corcov=corcov))
            os.system('cp {pf} {od}'.format(pf='/afs/cern.ch/user/m/mciprian/public/index.php',od=options.outdir))


    if options.showMoreCorrelated:
        pass
