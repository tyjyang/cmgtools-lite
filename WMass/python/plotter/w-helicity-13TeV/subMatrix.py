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
        return "W({l}#nu) fiducial".format(l=nn)

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
            nn = "{wch}#rightarrow{lep}#nu  {pol}: |y_{{W}}| #in [{ywl:1.2f},{ywh:1.2f}]".format(wch=wch,lep=lep,pol=pol,ywl=ywl,ywh=ywh)
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
            nn = "{wch}#rightarrow{lep}#nu: |#eta|-p_{{T}} #in [{etal:1.1f},{etah:1.1f}]-[{ptl:3g},{pth:3g}]".format(wch=wch,lep=lep,etal=etal,etah=etah,ptl=ptl,pth=pth)
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
    (options, args) = parser.parse_args()

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
        etaPtBinningVec = getDiffXsecBinning(options.etaptbinfile, "gen")
        genBins = templateBinning(etaPtBinningVec[0],etaPtBinningVec[1])
    elif options.ywbinfile:
        ybinfile = open(options.ywbinfile, 'r')
        genBins = eval(ybinfile.read())
        ybinfile.close()
    else:
        genBins = ""
    
    print "="*30
    print "Gen binning:"
    print "-"*30
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
    scalefactor_poi = 1./35900.
    if options.channel == "lep":
        scalefactor_poi /= 2.

    print "Preparing list of parameters to build the matrix ..."

    for ip1, p1 in enumerate(params):
        for ip2, p2 in enumerate(params):
            if options.type == 'toys':
                var = '({x}-{x0})*({y}-{y0})'.format(x=p1,x0=fitvals[p1],y=p2,y0=fitvals[p2])
                _tree.Draw('{var}>>h_{x}_{y}'.format(var=var,x=p1,y=p2))
                h = ROOT.gROOT.FindObject('h_{x}_{y}'.format(x=p1,y=p2)).Clone()
                cov [(p1,p2)] = h.GetMean()
                corr[(p1,p2)] = cov[(p1,p2)]/(fiterrs[p1]*fiterrs[p2])
            elif options.type == 'hessian':
                scalefactor = 1.0
                # normalization by bin area/width only for 2D/1D xsec, but might be added for rapidity as well
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
                    if p2.endswith(poiPostfix):     
                        scalefactor *= scalefactor_poi
                cov [(p1,p2)] = scalefactor * covmatrix .GetBinContent(indices[ip1],indices[ip2])
                corr[(p1,p2)] = corrmatrix.GetBinContent(indices[ip1],indices[ip2])
        

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
        params = sorted(params, key = lambda x: 0 if "right" in x else -1 if "left" in x else -1)

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
        new_x = niceName(x,genBins=genBins,forceLep=options.channel,drawRangeInLabel=True)
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
                lat.DrawLatex(0.15, 0.95, '#bf{CMS}') #it{Preliminary}')
                lat.DrawLatex(0.56, 0.95, '35.9 fb^{-1} (13 TeV)')

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
