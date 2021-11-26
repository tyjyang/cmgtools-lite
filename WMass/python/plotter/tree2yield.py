#!/usr/bin/env python
from math import *
import re
import os, os.path
from array import array
import logging
from collections import OrderedDict
import numbaDefines
import json

## safe batch mode
import sys
args = sys.argv[:]
sys.argv = ['-b']
import ROOT
sys.argv = args
ROOT.gROOT.SetBatch(True)
ROOT.PyConfig.IgnoreCommandLineOptions = True

ROOT.gInterpreter.ProcessLine(".O3")

from copy import *

from cutsFile import *
from fakeRate import *
#from mcCorrections import *

if "/functions_cc.so" not in ROOT.gSystem.GetLibraries(): 
    compileMacro("ccFiles/functions.cc")

if "/pileupWeights_cc.so" not in ROOT.gSystem.GetLibraries(): 
    compileMacro("ccFiles/pileupWeights.cc")

if "/jsonManager_cc.so" not in ROOT.gSystem.GetLibraries(): 
    compileMacro("ccFiles/jsonManager.cc")

if "/functionsWMass_cc.so" not in ROOT.gSystem.GetLibraries(): 
    compileMacro("ccFiles/functionsWMass.cc")

if "/RoccoR_cc.so" not in ROOT.gSystem.GetLibraries(): 
    compileMacro("ccFiles/RoccoR.cc")

if "/roccorManager_cc.so" not in ROOT.gSystem.GetLibraries(): 
    compileMacro("ccFiles/roccorManager.cc")

def setLogging(verbosity):
    verboseLevel = [logging.CRITICAL, logging.ERROR, logging.WARNING, logging.INFO, logging.DEBUG]
    logging.basicConfig(format='%(levelname)s: %(message)s', level=verboseLevel[min(4,verbosity)])

class PlotSpec:
    def __init__(self,name,expr,bins,opts):
        self.name = name
        self.expr = expr
        self.bins = bins
        self.opts = opts
        self.logs = {}
    def hasOption(self,name):
        return (name in self.opts)
    def getOption(self,name,default=None):
        return self.opts[name] if (name in self.opts) else default
    def setOption(self,name,value):
        self.opts[name] = value
    def setLog(self,name,value):
        self.logs[name] = value
    def allLogs(self):
        return self.logs.items()

def stylePlot(plot,spec,getOption):
        if "THn" in plot.ClassName():
            # to be further implemented
            nDim = plot.GetNdimensions()
            defaultTitle = ",".join(f"Dimension {idim}" for idim in range(nDim))
            titles = spec.getOption('NTitle',defaultTitle).split(",")                
            logging.debug(f"{plot.GetName()}: {titles}")
            for idim in range(nDim):
                if idim < len(titles):
                    plot.GetAxis(idim).SetTitle(f"{titles[idim]}")
                else:
                    plot.GetAxis(idim).SetTitle(f"Dimension {idim}")
        else:
            ## Sample specific-options, from self
            if getOption('FillColor',None) != None:
                plot.SetFillColor(getOption('FillColor',0))
                plot.SetFillStyle(getOption('FillStyle',1001))
            else:
                plot.SetFillStyle(0)
                plot.SetLineWidth(getOption('LineWidth',1))
            plot.SetLineColor(getOption('LineColor',1))
            plot.SetLineStyle(getOption('LineStyle',1))
            plot.SetMarkerColor(getOption('MarkerColor',1))
            plot.SetMarkerStyle(getOption('MarkerStyle',20))
            plot.SetMarkerSize(getOption('MarkerSize',1.1))
            ## Plot specific-options, from spec
            plot.GetYaxis().SetTitle(spec.getOption('YTitle',"Events"))
            plot.GetXaxis().SetTitle(spec.getOption('XTitle',spec.name))
            if "TH3" not in plot.ClassName():
                plot.GetXaxis().SetNdivisions(spec.getOption('XNDiv',510))
                plot.GetXaxis().SetMoreLogLabels(True)
            else:
                plot.GetZaxis().SetTitle(spec.getOption('ZTitle',"Events"))

def makeHistFromBinsAndSpec(name,expr,bins,plotspec):
        profile1D      = plotspec.getOption('Profile1D',False) if plotspec != None else False
        profile2D      = plotspec.getOption('Profile2D',False) if plotspec != None else False
        nvars = expr.replace("::","--").count(":")+1
        if nvars == 1 or (nvars == 2 and profile1D):
            if bins[0] == "[":
                edges = [ float(f) for f in bins[1:-1].split(",") ]
                if profile1D: 
                    histo = ROOT.TProfile(name,name,len(edges)-1,array('f',edges))
                else:
                    histo = ROOT.TH1D(name,name,len(edges)-1,array('f',edges))
            else:
                (nb,xmin,xmax) = bins.split(",")
                if profile1D:
                    histo = ROOT.TProfile(name,name,int(nb),float(xmin),float(xmax))
                else:
                    histo = ROOT.TH1D(name,name,int(nb),float(xmin),float(xmax))
        elif nvars == 2 or (nvars == 3 and profile2D):
            if bins[0] == "[":
                xbins, ybins = bins.split("*")
                xedges = [ float(f) for f in xbins[1:-1].split(",") ]
                yedges = [ float(f) for f in ybins[1:-1].split(",") ]
                if profile2D:
                    histo = ROOT.TProfile2D(name,name,len(xedges)-1,array('d',xedges),len(yedges)-1,array('d',yedges))
                else:
                    histo = ROOT.TH2D(name,name,len(xedges)-1,array('f',xedges),len(yedges)-1,array('f',yedges))
            else:
                (nbx,xmin,xmax,nby,ymin,ymax) = bins.split(",")
                if profile2D:
                    histo = ROOT.TProfile2D(name,name,int(nbx),float(xmin),float(xmax),int(nby),float(ymin),float(ymax))
                else:
                    histo = ROOT.TH2D(name,name,int(nbx),float(xmin),float(xmax),int(nby),float(ymin),float(ymax))
        elif nvars == 3:
            if bins[0] == "[":
                xbins, ybins, zbins = bins.split("*")
                xedges = [ float(f) for f in xbins[1:-1].split(",") ]
                yedges = [ float(f) for f in ybins[1:-1].split(",") ]
                zedges = [ float(f) for f in zbins[1:-1].split(",") ]
                histo = ROOT.TH3D(name,name,len(xedges)-1,array('f',xedges),len(yedges)-1,array('f',yedges),len(zedges)-1,array('f',zedges))
        elif nvars >= 4:
            if bins[0] == "[":
                logging.error("In makeHistFromBinsAndSpec(): only uniform binning currently supported for histograms with dimension > 3.")
                quit()
            else:
                #(nbx,xmin,xmax,nby,ymin,ymax,nbz,zmin,zmax,nbd4,d4min,d4max) = bins.split(",")
                #nbins = [nbx, nby, nbz, nbd4]
                #vmin = [xmin,ymin, zmin, d4min]
                #vmax = [xmax,ymax, zmax, d4max]
                #histo = ROOT.THnD(name,name, 4, array('i',nbins), array('f',vmin), array('f',vmax))
                binRanges = bins.split(",")
                nbins = []
                vmin = []
                vmax = []
                for iv in range(nvars):
                    trueIndex = 3*iv
                    nbins.append(int(binRanges[trueIndex]))
                    vmin.append(float(binRanges[trueIndex+1]))
                    vmax.append(float(binRanges[trueIndex+2]))
                #print(f"Defining {nvars}-dimensional histograms")    
                #print(nbins)
                #print(vmin)
                #print(vmax)
                histo = ROOT.THnD(name, name, nvars, array('i',nbins), array('d',vmin), array('d',vmax))
                #histo = ROOT.THnD(name, name, nvars, nbins, vmin, vmax)
                    
        # elif nvars == 5:
        #     if bins[0] == "[":
        #         logging.error("In makeHistFromBinsAndSpec(): only uniform binning supported for histograms with 5 dimensions.")
        #         quit()
        #     else:
        #         (nbx,xmin,xmax,nby,ymin,ymax,nbz,zmin,zmax,nbd4,d4min,d4max,nbd5,d5min,d5max) = bins.split(",")
        #         nbins = [nbx, nby, nbz, nbd4, nbd5]
        #         vmin = [xmin,ymin, zmin, d4min, d5min]
        #         vmax = [xmax,ymax, zmax, d4max, d5max]
        #         histo = ROOT.THnD(name,name, 5, array('i',nbins), array('f',vmin), array('f',vmax))
        else:
            raise RuntimeError("Can't make a plot with %d dimensions" % nvars)
        histo.Sumw2()
        return histo

def cropNegativeBins(histo):
            if "TH1" in histo.ClassName():
                for b in range(0,histo.GetNbinsX()+2):
                    if histo.GetBinContent(b) < 0: histo.SetBinContent(b, 0.0)
            elif "TH2" in histo.ClassName():
                for bx in range(0,histo.GetNbinsX()+2):
                    for by in range(0,histo.GetNbinsY()+2):
                        if histo.GetBinContent(bx,by) < 0: histo.SetBinContent(bx,by, 0.0)
            elif "TH3" in histo.ClassName():
                for bx in range(0,histo.GetNbinsX()+2):
                    for by in range(0,histo.GetNbinsY()+2):
                        for bz in range(0,histo.GetNbinsZ()+2):
                            if histo.GetBinContent(bx,by,bz) < 0: histo.SetBinContent(bx,by,bz, 0.0)


class TreeToYield:
    def __init__(self,root,options,scaleFactor=1.0,name=None,cname=None,settings={},objname=None, frienddir=None):
        self._name  = name  if name != None else root
        self._cname = cname if cname != None else self._name
        #print('in treetoyield init this is root and type(root)', root, type(root))
        self._fname = root
        self._sumGenWeights = None # filled with RDF later
        self._isInit = False
        self._options = options
        self._objname = objname if objname else options.obj
        self._frienddir = frienddir if frienddir else ""
        self._weight  = (options.weight and 'data' not in self._name)
        self._isdata = 'data' in self._name
        self._eraVFP = "preVFP" if "preVFP" in self._cname else "postVFP" if "postVFP" in self._cname else "all"
        self._lumiWeight = self._options.lumiWeight if self._options.lumiWeight != None else self._options.lumiDict[self._eraVFP] # defined also for data, but scaling done only for MC
        if self._isdata or not options.weightString:
            self._weightString = "1"
        else:
            self._weightString = "*".join("("+str(w)+")" for w in options.weightString)
        self._scaleFactor = scaleFactor
        self._settings = settings
        # for fake rate we may not need the older method, now that we have RDF    
        if 'FakeRate' in settings:
            self._FR = FakeRate(settings['FakeRate'],self._lumiWeight)
            ## add additional weight correction.
            ## note that the weight receives the other mcCorrections, but not itself
            frweight = self.adaptExpr(self._FR.weight(), cut=True)
            ## modify cuts to get to control region. order is important
            #self._mcCorrs = self._mcCorrs + self._FR.cutMods()  + self._FR.mods()
            self._weightString = self.adaptExpr(self._weightString, cut=True) + "* (" + frweight + ")"
            self._weight = True
        else:
            self._weightString = self.adaptExpr(self._weightString, cut=True)
        if self._options.forceunweight: self._weight = False
        for macro in self._options.loadMacro:
            libname = macro.replace(".cc","_cc.so").replace(".cxx","_cxx.so")
            if libname not in ROOT.gSystem.GetLibraries():
                ROOT.gROOT.ProcessLine(".L %s+" % macro);
        self._cutReport = None
        self._rdfDefs = OrderedDict()
        self._rdfAlias = {}
        if self._options.rdfDefineFile:
            self._rdfDefs = self.getRdfDefinitions(self._options.rdfDefineFile)
        if self._options.rdfDefine:
            self._rdfDefs.update(self.getRdfDefinitions(self._options.rdfDefine))

        if self._options.rdfAlias:
            self._rdfAlias = self.getRdfDefinitions(self._options.rdfAlias)
        self._rdfDefsWeightColumns = {}

            
    def getRdfDefinitions(self, data):
        lines = data
        if type(data) == str and os.path.isfile(data):
            with open(data) as f:
                lines = [x.strip() for x in f if not x.startswith("#") and len(x) > 0]

        ret = OrderedDict()
        for l in lines:                    
            # This is why you don't define your own config format
            l = l.replace("::", "--")
            tokens = [x.strip().replace("--", "::") for x in l.split(":")]
            if len(tokens) < 2:
                continue
            name = tokens[0]
            define = tokens[1]
            procRegexp = tokens[2] if len(tokens) > 2 else ".*"
            if len(tokens) == 4:
                altExpr = tokens[3]
            else:
                altExpr = None
                
            ret[name] = {"define" : define,
                         "procRegexp" : procRegexp,
                         "altExpr" : altExpr}
        return ret

    def printRdfDefinitions(self):
        logging.info("-"*40)
        logging.info("Have these RDF definitions")
        logging.info("-"*40)
        for key,value in self._rdfDefs.iteritems():
            logging.info("{define}) : {args} (for process {procRegexp})".format(name=key, **value))
        logging.info("-"*40)

    def getCutReport(self):
        return self._cutReport
    
    def setScaleFactor(self,scaleFactor):
        if (not self._options.forceunweight) and scaleFactor != 1: self._weight = True
        self._scaleFactor = scaleFactor
    def getScaleFactor(self):
        return self._scaleFactor
    def setRDF(self, rdf):
        self._fname = rdf
    def name(self):
        return self._name
    def cname(self):
        return self._cname
    def hasOption(self,name):
        return (name in self._settings)
    def getOption(self,name,default=None):
        if name in self._settings: return self._settings[name]
        return default
    def setOption(self,name,value):
        self._settings[name] = value
    def adaptDataMCExpr(self,expr):
        ret = expr
        if self._isdata:
            ret = re.sub(r'\$MC\{.*?\}', '', re.sub(r'\$DATA\{(.*?)\}', r'\1', expr));
        else:
            ret = re.sub(r'\$DATA\{.*?\}', '', re.sub(r'\$MC\{(.*?)\}', r'\1', expr));
        return ret
    def adaptExpr(self,expr,cut=False):
        ret = self.adaptDataMCExpr(expr)
        return ret
    def _init(self):
        self._tree  = self._fname ## just set the tree to the rdataframe
        #print 'done with init'
        self._isInit = True
        
    def getTree(self):
        if not self._isInit: self._init()
        return self._tree

    def prettyPrint(self,report):
        # maximum length of the cut descriptions
        clen = max([len(cut) for cut,yields in report]) + 3
        cfmt = "%%-%ds" % clen;

        fmtlen = 12
        nfmtL = "    %8d"
        nfmtS = "    %8.3f" if self._weight else nfmtL

        if self._options.errors:
            nfmtS+=u"%8.3f"
            nfmtL+=u"%8.3f"
            fmtlen+=8
        if self._options.fractions:
            nfmtS+="%7.1f%%"
            nfmtL+="%7.1f%%"
            fmtlen+=8

        logging.info("cut".center(clen),"yield".center(fmtlen))
        logging.info("-"*((fmtlen+1)+clen))
        for i,(cut,(nev,err)) in enumerate(report):
            logging.info(cfmt % cut,)
            den = report[i-1][1][0] if i>0 else 0
            fraction = nev/float(den) if den > 0 else 1
            toPrint = (nev,)
            if self._options.errors:    toPrint+=(err,)
            if self._options.fractions: toPrint+=(fraction*100,)
            if self._weight and nev < 1000: logging.info(nfmtS % toPrint,)
            else                          : logging.info(nfmtL % toPrint,)
            logging.info("")
    def _stylePlot(self,plot,spec):
        return stylePlot(plot,spec,self.getOption)
    def getManyPlots(self,plotspecs,cut,fsplit=None,closeTreeAfter=False):
        # getManyPlotsRaw should not trigger the loop yet, but refineManyPlots does, as it acts on the histograms
        # so, when using runGraphs, these should be called separately after collecting all histograms for all processes
        rets = self.getManyPlotsRaw(cut, plotspecs, fsplit=fsplit, closeTreeAfter=closeTreeAfter)
        # now other operations on histograms for this cname
        # need to suspend it if I want to use runGraph
        rets = self.refineManyPlots(rets, plotspecs)
        return rets
    def refineManyPlots(self, rets, plotspecs):

        # operations on histograms, which trigger the loop, should be included inside here if possible
        
        # normalize histograms by sum of gen weights if appropriate
        # crop bins with negative content if desired
        # and change directory for them        
        if not self._isdata and self._options.weight:
            message = f"Process {self._name} ({self._cname}) -> lumi weight = {self._lumiWeight}"
            sumGenWeights = self._sumGenWeights.GetValue() if self._sumGenWeights else 1.0
            message += f"    sum gewWeights = {sumGenWeights}"
            print(message)
            for histo in rets:
                # this is going to trigger the loop, but now all histograms exist so it is ok
                histo.Scale(self._lumiWeight/sumGenWeights)
                self.negativeCheck(histo)
                if not histo.ClassName() == 'THnT<double>':
                    histo.SetDirectory(0)
        else:
            for histo in rets:
                if not histo.ClassName() == 'THnT<double>':
                    histo.SetDirectory(0)
                    
        # not all histograms are made for all plots, so need to prune the list of pspec to be used below
        plotspecsPruned = [p for p in plotspecs if re.match(p.getOption('ProcessRegexp', ".*"), self._name)]
                    
        # fold overflow
        for iret,ret in enumerate(rets):
            if ret.ClassName() in [ "TH1F", "TH1D" ] :
                n = ret.GetNbinsX()
                if plotspecsPruned[iret].getOption('IncludeOverflows',True) and ("TProfile" not in ret.ClassName()):
                    ret.SetBinContent(1,ret.GetBinContent(0)+ret.GetBinContent(1))
                    ret.SetBinContent(n,ret.GetBinContent(n+1)+ret.GetBinContent(n))
                    ret.SetBinError(1,hypot(ret.GetBinError(0),ret.GetBinError(1)))
                    ret.SetBinError(n,hypot(ret.GetBinError(n+1),ret.GetBinError(n)))
                    ret.SetBinContent(0,0)
                    ret.SetBinContent(n+1,0)
                    ret.SetBinContent(0,0)
                    ret.SetBinContent(n+1,0)
                if plotspecsPruned[iret].getOption('IncludeOverflow',False) and ("TProfile" not in ret.ClassName()):
                    ret.SetBinContent(n,ret.GetBinContent(n+1)+ret.GetBinContent(n))
                    ret.SetBinError(n,hypot(ret.GetBinError(n+1),ret.GetBinError(n)))
                    ret.SetBinContent(n+1,0)
                    ret.SetBinContent(n+1,0)
                if plotspecsPruned[iret].getOption('IncludeUnderflow',False) and ("TProfile" not in ret.ClassName()):
                    ret.SetBinContent(1,ret.GetBinContent(0)+ret.GetBinContent(1))
                    ret.SetBinError(1,hypot(ret.GetBinError(0),ret.GetBinError(1)))
                    ret.SetBinContent(0,0)
                    ret.SetBinContent(0,0)
                rebin = plotspecsPruned[iret].getOption('rebinFactor',0)
                if plotspecsPruned[iret].bins[0] != "[" and rebin > 1 and n > 5:
                    while n % rebin != 0: rebin -= 1
                    if rebin != 1: ret.Rebin(rebin)
                if plotspecsPruned[iret].getOption('Density',False):
                    for b in range(1,n+1):
                        ret.SetBinContent( b, ret.GetBinContent(b) / ret.GetXaxis().GetBinWidth(b) )
                        ret.SetBinError(   b, ret.GetBinError(b) / ret.GetXaxis().GetBinWidth(b) )
            self._stylePlot(ret,plotspecsPruned[iret])
            ret._cname = self._cname
        return rets
    # apparently not used anywhere
    def getWeightForCut(self,cut):
        if self._weight:
            cut = "(%s)*(%s)*(%s)" % (self._weightString, self._scaleFactor, self.adaptExpr(cut,cut=True))
        else:
            cut = self.adaptExpr(cut,cut=True)
        return cut
    def defineColumnFromExpression(self, colExpr, colName, expr):
        if expr in list(colExpr.keys()):
            logging.debug(f"Expression {expr} already exists, I will reuse it without defining new column")
            colName = colExpr[expr]
        else:                    
            # if we are here the expression changed, but it could be the column name is the same as an existing ones (e.g. for the event weight string, for which the column is defined for each plot but the column name does not depend on it, so need to rename the column
            if colName in list(colExpr.values()):
                matchedCnames = list(filter(lambda x : re.match(colName+"_\d+",x), list(colExpr.values())))
                cnindex = len(matchedCnames)
                oldname = colName
                colName = oldname + "_%d" % (cnindex)
                logging.debug("Column %s already defined, name changed to %s" % (oldname, colName))
            else:
                logging.debug("Defining new column %s" % colName)
            colExpr[expr] = colName
            self._tree = self._tree.Define(colName, expr)
        return (colName, colExpr)
    def getManyPlotsRaw(self,cut,plotspecs,fsplit=None,closeTreeAfter=False):
        if not self._isInit: self._init()
        retlist = []

        ## define the Sum of genWeights before the filter
        if not self._isdata and self._options.weight:
            self._sumGenWeights = self._tree.Sum("genWeightClip") # uses floats
            # self._sumGenWeights = ROOT.getRDFcolumnSum(self._tree,"myweight") # uses doubles, but to be fixed
            ## now filter the rdf with just the cut

        # define isData to handle data and MC dynamically
        # needed before filters, as we use this for the json selection
        if self._cname.startswith("data"):
            self._tree = self._tree.Define("isData", "Data")
        else:
            self._tree = self._tree.Define("isData", "MC")

        # preselection, if any, before defines
        if self._options.printYieldsRDF:
            logging.debug('='*30)
            logging.debug(" Preselection before defines")
            logging.debug('-'*30)
            for (cutname,cutexpr,extra) in cut.cuts():
                if "BEFOREDEFINE" in extra and extra["BEFOREDEFINE"] == True:
                    logging.debug(f" {cutname}: {cutexpr}")
                    self._tree = self._tree.Filter(cutexpr,cutname)
            logging.debug('='*30)
                    
        # this column must be defined according to the process name, which is expected to include pre or postVFP
        # the eraVFP column will be used in calls to some functions to specify the era (e.g. for scale factors)
        if "preVFP" in self._cname:
            self._tree = self._tree.Define("eraVFP", "BToF")
        elif "postVFP" in self._cname:
            self._tree = self._tree.Define("eraVFP", "GToH")
        else:
            self._tree = self._tree.Define("eraVFP", "BToH")

        # defines
        for name, entry in self._rdfDefs.items():
            if re.match(entry["procRegexp"], self._cname):
                logging.debug("Defining %s as %s" % (name, entry["define"]))
                self._tree = self._tree.Define(name, entry["define"])
            elif entry["altExpr"]:
                self._tree = self._tree.Define(name, entry["altExpr"])
        # aliases
        for name,entry in self._rdfAlias.items():
            if re.match(entry["procRegexp"], self._cname):
                self._tree = self._tree.Alias(name, entry["define"])

        # filters
        cutsAsWeight = {}
        if self._options.printYieldsRDF:
            logging.debug('='*30)
            logging.debug(" Selection after defines")
            logging.debug('-'*30)
            for (cutname,cutexpr,extra) in cut.cuts():
                if "BEFOREDEFINE" not in extra:
                    logging.debug(f" {cutname}: {cutexpr}; {extra}")
                    # in case one cut has to be skipped or overridden, it still need to be added as dummy in the list of filters, or RDF may complain when doing all processes in the same loop
                    if "EXCLUDEPROCESS" in extra and re.match(extra["EXCLUDEPROCESS"],self._cname):
                        if "ALTEXPR" in extra:
                            self._tree = self._tree.Filter(extra["ALTEXPR"],cutname)
                        else:
                            self._tree = self._tree.Filter("1>0",cutname)
                    else:
                        if "ASWEIGHT" in extra and extra["ASWEIGHT"] == True:
                            # if it is data, set weight flag to true, it might have been false unless using some scaling from MCA file
                            # usually one could use ASWEIGHT_MCONLY to enable it only for MC and not data
                            if self._cname.startswith("data"):
                                self._weight = True
                            self._tree = self._tree.Filter("1>0",cutname)
                            cutsAsWeight[cutname] = cutexpr
                        elif "ASWEIGHT_MCONLY" in extra and extra["ASWEIGHT_MCONLY"] == True:
                            if self._cname.startswith("data"):
                                self._tree = self._tree.Filter(cutexpr,cutname)
                            else:
                                self._tree = self._tree.Filter("1>0",cutname)
                                cutsAsWeight[cutname] = cutexpr
                        else:
                            self._tree = self._tree.Filter(cutexpr,cutname)

            logging.debug('='*30)
        else:
            # usually we never use this piece, it is faster to have multiple filters rather than a single one, it runs much faster
            fullcut = cut
            self._tree = self._tree.Filter(fullcut)
            
        # do not call it, or it will trigger the loop now
        #print "sumGenWeights"
        #print sumGenWeights.GetValue()
        histos = []

        # keep list of defined columns used as variables to fill histograms, to avoid defining the same
        # expression multiple times, which may affect performances
        column_expr = {}

        ## weight has a common part for each plot (lumiWeight applied at the end as a scaling factor)
        ### Note:
        ### data is usually unweighted, but one can add a scaling using the 3rd field of the process line in the MCA file, or using AddWeight (which adds it to the 3rd field, creating it if only two fields were present)
        ### when this happens, the weight is set in mcAnalysis.py as a scalefactor, so it will be included in self._scaleFactor. This also activate self._weight for data as well, but self._weightString
        ### (which usually is passed as option -W) is set to "1.0" for data, so wgtCommon == (1.0)*self._scaleFactor
        if self._weight:
            wgtCommon = "(%s)*(%s)" % (self._weightString, self._scaleFactor)
            if "ReplaceWeight" in self._settings:
                oldw,neww = self._settings["ReplaceWeight"].split("->")
                if oldw not in wgtCommon:
                    logging.warning(f"'{oldw}' not found in wgtCommon for process {self._name}. Please check that the weight string is defined as you expect.")
                    quit()                
                wgtCommon = wgtCommon.replace(oldw,neww) # this is for all plots, so we update the common weight directly
        else:
            wgtCommon = '1.'
        for cutname in cutsAsWeight.keys():
            wgtCommon += f"*({cutsAsWeight[cutname]})"
        #print(f"Process {self._name}, wgtCommon = {wgtCommon}")
            
            
        for plotspec in plotspecs:

            # check if histogram has to be done for this process
            # if the flag is not set it is assumed it has to be done
            if plotspec.getOption('ProcessRegexp', None):            
                regexp = plotspec.getOption('ProcessRegexp', ".*")
                if re.match(regexp, self._name):
                    logging.debug("Going to plot %s for process %s" % (plotspec.name,self._name))
                else:
                    logging.debug("Skipping plot %s for process %s" % (plotspec.name,self._name))
                    continue
            else:
                logging.debug("Going to plot %s for process %s (all accepted)" % (plotspec.name,self._name))
            tmp_histo = makeHistFromBinsAndSpec(self._cname+'_'+plotspec.name,plotspec.expr,plotspec.bins,plotspec)

            tmp_weight = self._cname+'_weight'
            #print "In getManyPlotsRaw"
            #print(tmp_weight)

            if plotspec.getOption('AddWeight', None):
                wgt = wgtCommon + "*(%s)" % plotspec.getOption('AddWeight','1')
            else:
                wgt = wgtCommon
            if plotspec.getOption('ReplaceWeight', None):
                # modify wgt defined just above, even though usually AddWeight and ReplaceWeight are not used together, but it can happen that they are
                oldw,neww = plotspec.getOption('ReplaceWeight', "1->1").split("->")
                if oldw not in wgt:
                    logging.warning(f"'{oldw}' not found in wgt for plot {plotspec.name} and process {self._name}. Please check that the weight string is defined as you expect.")
                    quit()
                wgt = wgt.replace(oldw,neww) 
            if plotspec.getOption('ReplaceCutByName', None):
                cutname,neww = plotspec.getOption('ReplaceCutByName', "alwaystrue->1").split("->")
                if cutname in cutsAsWeight.keys():
                    if cutsAsWeight[cutname] not in wgt:
                        logging.warning(f"'{cutsAsWeight[cutname]}' not found in wgt for plot {plotspec.name} and process {self._name}. Please check that all flags are correctly set in the cut file.")
                        quit()
                    wgt = wgt.replace(cutsAsWeight[cutname],neww)
                    #print(f"Replacing weight per plot ({plotspec.name}): wgt = {wgt}")
                    #print()
            (tmp_weight, self._rdfDefsWeightColumns) = self.defineColumnFromExpression(self._rdfDefsWeightColumns, tmp_weight, wgt)
            
            logging.debug("%s weight string = %s " % (tmp_weight, wgt))

            # TODO
            # using THxDModel is not necessary, will change it at some point
            if tmp_histo.ClassName() == 'TH1D':
                tmp_histo_model = ROOT.RDF.TH1DModel(tmp_histo)
                expr_x = plotspec.expr
                tmp_varx    = self._cname+'_'+plotspec.name+'_varx'
                (tmp_varx, column_expr) = self.defineColumnFromExpression(column_expr, tmp_varx, expr_x)
                tmp_histo  = self._tree.Histo1D(tmp_histo_model, tmp_varx, tmp_weight)
            elif tmp_histo.ClassName() == 'TH2D':
                tmp_histo_model = ROOT.RDF.TH2DModel(tmp_histo)
                expr_y,expr_x = plotspec.expr.split(':')
                tmp_varx    = self._cname+'_'+plotspec.name+'_varx'
                tmp_vary    = self._cname+'_'+plotspec.name+'_vary'
                (tmp_varx, column_expr) = self.defineColumnFromExpression(column_expr, tmp_varx, expr_x)
                (tmp_vary, column_expr) = self.defineColumnFromExpression(column_expr, tmp_vary, expr_y)
                tmp_histo  = self._tree.Histo2D(tmp_histo_model, tmp_varx, tmp_vary, tmp_weight)
            elif tmp_histo.ClassName() == 'TH3D':
                tmp_histo_model = ROOT.RDF.TH3DModel(tmp_histo)
                expr_z,expr_y,expr_x = plotspec.expr.split(':')
                tmp_varx    = self._cname+'_'+plotspec.name+'_varx'
                tmp_vary    = self._cname+'_'+plotspec.name+'_vary'
                tmp_varz    = self._cname+'_'+plotspec.name+'_varz'
                (tmp_varx, column_expr) = self.defineColumnFromExpression(column_expr, tmp_varx, expr_x)
                (tmp_vary, column_expr) = self.defineColumnFromExpression(column_expr, tmp_vary, expr_y)
                (tmp_varz, column_expr) = self.defineColumnFromExpression(column_expr, tmp_varz, expr_z)
                tmp_histo  = self._tree.Histo3D(tmp_histo_model, tmp_varx, tmp_vary, tmp_varz, tmp_weight)
            elif tmp_histo.ClassName() == 'THnT<double>':
                tmp_histo_model = ROOT.RDF.THnDModel(tmp_histo)                
                nDim = tmp_histo.GetNdimensions()
                #print(f"plotspec.expr = {plotspec.expr}")
                expr_n = plotspec.expr.split(':')[::-1]
                #print(f"expr_n = {expr_n}")
                tmp_varn = []
                for idim in range(nDim):
                    tmp_varn.append(f"{self._cname}_{plotspec.name}_var{idim}d")
                    (tmp_varn[idim], column_expr) = self.defineColumnFromExpression(column_expr, tmp_varn[idim], expr_n[idim])
                # append weight to list of columns to use, since this is how the constructor works    
                tmp_varn.append(tmp_weight)
                tmp_histo  = self._tree.HistoND(tmp_histo_model, tmp_varn)
                #print(f"Creating RDF object with type {type(tmp_histo)}")
                
            histos.append(tmp_histo)

        if self._options.printYieldsRDF:
            #print("="*30)
            #print("Preparing yields for %s" % self._cname)
            #print("-"*30)
            self._cutReport = self._tree.Report()
            #cutReport.Print()
            #print("="*30)
        if closeTreeAfter: self._tfile.Close()
        return histos

    def negativeCheck(self,histo):
        if not self._options.allowNegative and not self._name in self._options.negAllowed: 
            cropNegativeBins(histo)
    def __str__(self):
        mystr = ""
        mystr += str(self._fname) + '\n'
        #mystr += str(self._tfile) + '\n'
        mystr += str(self._weight) + '\n'
        mystr += str(self._scaleFactor)
        return mystr
def _copyPlotStyle(self,plotfrom,plotto):
        plotto.SetFillStyle(plotfrom.GetFillStyle())
        plotto.SetFillColor(plotfrom.GetFillColor())
        plotto.SetMarkerStyle(plotfrom.GetMarkerStyle())
        plotto.SetMarkerColor(plotfrom.GetMarkerColor())
        plotto.SetMarkerSize(plotfrom.GetMarkerSize())
        plotto.SetLineStyle(plotfrom.GetLineStyle())
        plotto.SetLineColor(plotfrom.GetLineColor())
        plotto.SetLineWidth(plotfrom.GetLineWidth())
        plotto.GetXaxis().SetTitle(plotfrom.GetXaxis().GetTitle())
        plotto.GetYaxis().SetTitle(plotfrom.GetYaxis().GetTitle())
        plotto.GetZaxis().SetTitle(plotfrom.GetZaxis().GetTitle())
        plotto.GetXaxis().SetNdivisions(plotfrom.GetXaxis().GetNdivisions())
        plotto.GetYaxis().SetNdivisions(plotfrom.GetYaxis().GetNdivisions())
        plotto.GetZaxis().SetNdivisions(plotfrom.GetZaxis().GetNdivisions())

def addTreeToYieldOptions(parser):
    parser.add_argument('--lumi-dict', dest='lumiDict', type=json.loads, default='{"all":36.3,"preVFP":19.5,"postVFP":16.8}', help="default luminosity values for inclusive samples or pre|postVFP, chosen according to component name in MCA text file")
    parser.add_argument('--lumi-weight', dest='lumiWeight', type=float, default=None, help="Use this value as lumi weight (for the label in plot one needs option -l defined in mcPlots.py)")
    parser.add_argument("-u", "--unweight", dest="weight", action="store_false", help="Don't use weights (in MC events), note weights are still used if a fake rate file is given");
    parser.add_argument("--uf", "--unweight-forced", dest="forceunweight", action="store_true", help="Do not use weight even if a fake rate file is given.");
    parser.add_argument("-W", "--weightString", dest="weightString", action="append", default=[], help="Use weight (in MC events), can specify multiple times");
    parser.add_argument("-f", "--final", action="store_true", help="Just compute final yield after all cuts (no longer active with RDF)");
    parser.add_argument("-e", "--errors", action="store_true", help="Include uncertainties in the reports");
    parser.add_argument("--tf", "--text-format", dest="txtfmt", type=str, default="txt", choices=["txt","tsv","csv","dsv","ssv"], help="Output format: txt,tsv,csv,dsv,ssv");
    parser.add_argument("-S", "--start-at-cut", dest="startCut", type=str, help="Run selection starting at the cut matched by this regexp, included.") 
    parser.add_argument("-U", "--up-to-cut", dest="upToCut", type=str, help="Run selection only up to the cut matched by this regexp, included.") 
    parser.add_argument("-X", "--exclude-cut", dest="cutsToExclude", action="append", default=[], help="Cuts to exclude (regexp matching cut name), can specify multiple times.") 
    parser.add_argument("-E", "--enable-cut", dest="cutsToEnable", action="append", default=[], help="Cuts to enable if they were disabled in the cut file (regexp matching cut name), can specify multiple times.") 
    parser.add_argument("-I", "--invert-cut", dest="cutsToInvert", action="append", default=[], help="Cuts to invert (regexp matching cut name), can specify multiple times.") 
    parser.add_argument("-R", "--replace-cut", dest="cutsToReplace", action="append", default=[], nargs=3, help="Cuts to invert (regexp of old cut name, new name, new cut); can specify multiple times.") 
    parser.add_argument("-A", "--add-cut", dest="cutsToAdd", action="append", default=[], nargs=3, help="Cuts to insert (regexp of cut name after which this cut should go, new name, new cut); can specify multiple times.") 
    parser.add_argument("--RV", "--replace-var", dest="varsToReplace", action="append", default=[], nargs=2, help="Variable to replace(old var name, new var name); can specify multiple times.")
    parser.add_argument("--obj", "--objname", dest="obj", default='Events', help="Pattern for the name of the TTree inside the file");
    parser.add_argument("-G", "--no-fractions", dest="fractions", action="store_false", default=True, help="Don't print the fractions");
    parser.add_argument("--neg", "--allow-negative-results", dest="allowNegative", action="store_true", default=False, help="If the total yield is negative, keep it so rather than truncating it to zero") 
    parser.add_argument("--neglist", dest="negAllowed", action="append", default=[], help="Give process names where negative values are allowed")
    parser.add_argument("--max-entries", dest="maxEntries", default=1000000000000, type=int, help="Max entries to process in each tree") 
    parser.add_argument("--rdf-range", dest="rdfRange", default=[], nargs=3, type=int, help="Arguments to use RDF::Range to sets entries to process with RDataFrame, passing begin, end, stride (it sets multithreading to use 1 thread)") 
    parser.add_argument("--max-entries-not-data", dest="maxEntriesNotData", action="store_true", help="When --max-entries is used, make if effective only for non data processes (needed for some tests because MC is rescaled to luminosity, data cannot)") 
    parser.add_argument("-L", "--load-macro", dest="loadMacro", type=str, action="append", default=[], help="Load the following macro, with .L <file>+");
    parser.add_argument("--rdf-define-file", dest="rdfDefineFile", type=str, default="", help="Load file with some definitions to be used in RDataFrame");
    parser.add_argument("--rdf-define", dest="rdfDefine", type=str, action="append", default=[], help="Add some definitions to be used in RDataFrame on the fly, formatted as 'name:definition:regularEpressionForCname'. If the regexp is '.*' it can be omitted. Note that these can override those passed with --rdf-define-file");
    parser.add_argument("--rdf-alias", dest="rdfAlias", type=str, action="append", default=[], help="Use Alias to rename columns in RDataFrame on the fly, formatted as 'newname:oldname:regularEpressionForCname'. If the regexp is '.*' it can be omitted. It is evaluated after the Define commands");

def mergeReports(reports):
    import copy
    one = copy.deepcopy(reports[0])
    for i,(c,x) in enumerate(one):
        one[i][1][1] = pow(one[i][1][1], 2)
    for two in reports[1:]:
        for i,(c,x) in enumerate(two):
            one[i][1][0] += x[0]
            one[i][1][1] += pow(x[1],2)
            one[i][1][2] += x[2]
    for i,(c,x) in enumerate(one):
        one[i][1][1] = sqrt(one[i][1][1])
    return one

def mergePlots(name,plots):
    one = plots[0].Clone(name)
    if "TGraph" in one.ClassName():
        others = ROOT.TList()
        for two in plots[1:]: 
            twoc = two.Clone(name+'foobar')
            others.Add(twoc)
        one.Merge(others)
    else:         
        for two in plots[1:]: 
            twoc = two.Clone(name+'foobar')
            one.Add(twoc)
    return one

