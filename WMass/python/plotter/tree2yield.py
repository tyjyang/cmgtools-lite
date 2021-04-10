#!/usr/bin/env python
from math import *
import re
import os, os.path
from array import array
import logging
from collections import OrderedDict
import numbaDefines

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
from mcCorrections import *

if "/functions_cc.so" not in ROOT.gSystem.GetLibraries(): 
    compileMacro("ccFiles/functions.cc")

if "/jsonManager_cc.so" not in ROOT.gSystem.GetLibraries(): 
    compileMacro("ccFiles/jsonManager.cc")

if "/w-mass-13TeV/functionsWMass_cc.so" not in ROOT.gSystem.GetLibraries(): 
    compileMacro("ccFiles/functionsWMass.cc")

def setLogging(verbosity):
    verboseLevel = [logging.CRITICAL, logging.ERROR, logging.WARNING, logging.INFO, logging.DEBUG]
    logging.basicConfig(format='%(levelname)s: %(message)s', level=verboseLevel[min(4,verbosity)])

# not sure we still need it for RDF
def scalarToVector(x):
    x0 = x
    x = re.sub(r"(LepGood|LepCorr|JetFwd|Jet|JetClean|Jet_Clean|GenTop|SV|PhoGood|TauGood|Tau)(\d)_(\w+)", lambda m : "%s_%s[%d]" % (m.group(1),m.group(3),int(m.group(2))-1), x)
    x = re.sub(r"\bmet\b", "met_pt", x)
    return x

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
        if "TH3" not in plot.ClassName():
            plot.GetYaxis().SetTitle(spec.getOption('YTitle',"Events"))
            plot.GetXaxis().SetTitle(spec.getOption('XTitle',spec.name))
            plot.GetXaxis().SetNdivisions(spec.getOption('XNDiv',510))
            plot.GetXaxis().SetMoreLogLabels(True)

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
            ez,ey,ex = [ e.replace("--","::") for e in expr.replace("::","--").split(":") ]
            if bins[0] == "[":
                xbins, ybins, zbins = bins.split("*")
                xedges = [ float(f) for f in xbins[1:-1].split(",") ]
                yedges = [ float(f) for f in ybins[1:-1].split(",") ]
                zedges = [ float(f) for f in zbins[1:-1].split(",") ]
                histo = ROOT.TH3D(name,name,len(xedges)-1,array('f',xedges),len(yedges)-1,array('f',yedges),len(zedges)-1,array('f',zedges))
            else:
                (nbx,xmin,xmax,nby,ymin,ymax,nbz,zmin,zmax) = bins.split(",")
                histo = ROOT.TH3D(name,name,int(nbx),float(xmin),float(xmax),int(nby),float(ymin),float(ymax),int(nbz),float(zmin),float(zmax))
            histo.GetXaxis().SetTitle(ex)
            histo.GetYaxis().SetTitle(ey)
            histo.GetZaxis().SetTitle(ez)
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
        self._sumGenWeights = None # filled with RDF later, and used if options.sumGenWeightFromHisto is False
        self._isInit = False
        self._options = options
        self._objname = objname if objname else options.obj
        self._frienddir = frienddir if frienddir else ""
        self._weight  = (options.weight and 'data' not in self._name)
        self._isdata = 'data' in self._name
        if self._isdata or not options.weightString:
            self._weightString = "1"
        else:
            self._weightString = "*".join("("+str(w)+")" for w in options.weightString)
        self._scaleFactor = scaleFactor
        self._fullYield = 0 # yield of the full sample, as if it passed the full skim and all cuts
        self._fullNevt = 0 # number of events of the full sample, as if it passed the full skim and all cuts
        self._settings = settings
        loadMCCorrections(options)            ## make sure this is loaded
        self._mcCorrs = globalMCCorrections() ##  get defaults
        if 'SkipDefaultMCCorrections' in settings: ## unless requested to 
            self._mcCorrs = []                     ##  skip them
        if self._isdata: 
# bug: does not work
#            self._mcCorrs = [c for c in self._mcCorrs if c.alsoData] ## most don't apply to data, some do 
            newcorrs=[]
            for corr in self._mcCorrs:
                newcorr = copy(corr)
                newlist = copy(newcorr._corrections)
                newlist = [icorr for icorr in newlist if icorr.alsoData]
                newcorr._corrections = newlist
                newcorrs.append(newcorr)
            self._mcCorrs=newcorrs
        if 'MCCorrections' in settings:
            self._mcCorrs = self._mcCorrs[:] # make copy
            for cfile in settings['MCCorrections'].split(','): 
                self._mcCorrs.append( MCCorrections(cfile) )
        if self._mcCorrs and self._scaleFactor and self._scaleFactor != 1.0:
            # apply MC corrections to the scale factor
            self._scaleFactor = self.adaptExpr(self._scaleFactor, cut=True)
        if 'FakeRate' in settings:
            self._FR = FakeRate(settings['FakeRate'],self._options.lumi)
            ## add additional weight correction.
            ## note that the weight receives the other mcCorrections, but not itself
            frweight = self.adaptExpr(self._FR.weight(), cut=True)
            ## modify cuts to get to control region. order is important
            self._mcCorrs = self._mcCorrs + self._FR.cutMods()  + self._FR.mods()
            self._weightString = self.adaptExpr(self._weightString, cut=True) + "* (" + frweight + ")"
            self._weight = True
        else:
            self._weightString = self.adaptExpr(self._weightString, cut=True)
        if self._options.forceunweight: self._weight = False
        for macro in self._options.loadMacro:
            libname = macro.replace(".cc","_cc.so").replace(".cxx","_cxx.so")
            if libname not in ROOT.gSystem.GetLibraries():
                ROOT.gROOT.ProcessLine(".L %s+" % macro);
        self._appliedCut = None
        self._cutReport = None
        self._elist = None
        self._entries = None
        self._rdfDefs = {}
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
            name = tokens[0]
            define = tokens[1]
            procRegexp = tokens[2] if len(tokens) > 2 else ".*"

            ret[name] = {"define" : define, "procRegexp" : procRegexp}
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
    
        #print "Done creation  %s for task %s in pid %d " % (self._fname, self._name, os.getpid())
    def setScaleFactor(self,scaleFactor,mcCorrs=True):
        if (not self._options.forceunweight) and scaleFactor != 1: self._weight = True
        if mcCorrs and self._mcCorrs and scaleFactor and scaleFactor != 1.0:
            # apply MC corrections to the scale factor
            self._scaleFactor = self.adaptExpr(scaleFactor, cut=True)
        else:
            self._scaleFactor = scaleFactor
    def getScaleFactor(self):
        return self._scaleFactor
    def setRDF(self, rdf):
        self._fname = rdf
    def setFullYield(self,fullYield):
        self._fullYield = fullYield
    def setFullNevt(self,fullNevt):
        self._fullNevt = fullNevt
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
        for mcc in self._mcCorrs:
            ret = mcc(ret,self._name,self._cname,cut)
        return ret
    def _init(self):
        self._tree  = self._fname ## just set the tree to the rdataframe

        # all the part below for friends might be removed
        self._friends = []
        friendOpts = self._options.friendTrees[:]
        friendOpts += [ ('sf/t', d+"/evVarFriend_{cname}.root") for d in self._options.friendTreesSimple]
        friendOpts += (self._options.friendTreesData if self._isdata else self._options.friendTreesMC)
        friendOpts += [ ('sf/t', d+"/evVarFriend_{cname}.root") for d in (self._options.friendTreesDataSimple if self._isdata else self._options.friendTreesMCSimple) ]
        if 'Friends' in self._settings: friendOpts += self._settings['Friends']
        if 'FriendsSimple' in self._settings: friendOpts += [ ('sf/t', d+"/evVarFriend_{cname}.root") for d in self._settings['FriendsSimple'] ]
        for tf_tree,tf_file in friendOpts:
            # print 'Adding friend',tf_tree,tf_file
            basepath = None
            if self._options.nanoaodTree:
                friend_cname = os.path.basename(self._cname)
            else:
                for treepath in getattr(self._options, 'path', []):
                    if self._cname in os.listdir(treepath):
                        basepath = treepath
                        break
                if not basepath:
                    raise RuntimeError("%s -- ERROR: %s process not found in paths (%s)" % (__name__, cname, repr(options.path)))
                friend_cname = self._cname
            tf_filename = tf_file.format(name=self._name, cname=friend_cname, P=basepath, fd=self._frienddir)
            if os.path.exists(tf_filename+".url"):
                tf_filename = open(tf_filename+".url","r").readline().strip()
            # print 'Adding friend',tf_tree,tf_filename
            tf = self._tree.AddFriend(tf_tree, tf_filename),
            self._friends.append(tf)
        #print 'done with init'
        self._isInit = True
        
    def getTree(self):
        if not self._isInit: self._init()
        return self._tree

    def getEntries(self,useEList=True,closeFileAfterwards=True):
        if useEList and self._elist: 
            return self._elist.GetN()
        if self._entries is None:
            if closeFileAfterwards and (not self._isInit):
                if "root://" in self._fname: ROOT.gEnv.SetValue("XNet.Debug", -1); # suppress output about opening connections
                tfile = ROOT.TFile.Open(self._fname)
                if not tfile: raise RuntimeError("Cannot open %s\n" % self._fname)
                t = tfile.Get(self._objname)
                if not t: raise RuntimeError("Cannot find tree %s in file %s\n" % (self._objname, self._fname))
                self._entries = t.GetEntries()
            else:
                self._entries = self.getTree().GetEntries()
        return self._entries
    def hasEntries(self,useEList=True):
        if useEList and self._elist:
            return True
        return (self._entries != None)
    def setEntries(self,entries):
        self._entries = entries
    def getYields(self,cuts,noEntryLine=False,fsplit=None):
        if not self._isInit: self._init()
        report = []; cut = ""
        cutseq = [ ['entry point','1'] ]
        if noEntryLine: cutseq = []
        sequential = False
        if self._options.nMinusOne or self._options.nMinusOneInverted: 
            if self._options.nMinusOneSelection:
                cutseq = cuts.nMinusOneSelectedCuts(self._options.nMinusOneSelection,inverted=self._options.nMinusOneInverted)
            else:
                cutseq = cuts.nMinusOneCuts(inverted=self._options.nMinusOneInverted)
            cutseq += [ ['all',cuts.allCuts()] ]
            sequential = False
        elif self._options.final:
            cutseq = [ ['all', cuts.allCuts()] ]
        else:
            cutseq += cuts.cuts();
            sequential = True
        for cn,cv in cutseq:
            if sequential:
                if cut: cut += " && "
                cut += "(%s)" % cv
            else:
                cut = cv
            report.append((cn,self._getYield(self._tree,cut,fsplit=fsplit)))
        return report
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
            if self._options.nMinusOne or self._options.nMinusOneInverted: 
                fraction = report[-1][1][0]/nev if nev > 0 else 1
            toPrint = (nev,)
            if self._options.errors:    toPrint+=(err,)
            if self._options.fractions: toPrint+=(fraction*100,)
            if self._weight and nev < 1000: logging.info(nfmtS % toPrint,)
            else                          : logging.info(nfmtL % toPrint,)
            logging.info("")
    def _getCut(self,cut,noweight=False):
        if self._weight and not noweight:
            if self._isdata: cut = "(%s)     *(%s)*(%s)" % (self._weightString,                    self._scaleFactor, self.adaptExpr(cut,cut=True))
            else:            cut = "(%s)*(%s)*(%s)*(%s)" % (self._weightString,self._options.lumi, self._scaleFactor, self.adaptExpr(cut,cut=True))
        else: 
            cut = self.adaptExpr(cut,cut=True)
        if self._options.doS2V:
            cut  = scalarToVector(cut)
        return cut
    def _getYield(self,tree,cut,fsplit=None,cutNeedsPreprocessing=True):
        cut = self._getCut(cut) if cutNeedsPreprocessing else cut
        if self._weight:
#            print cut
            ROOT.gROOT.cd()
            if ROOT.gROOT.FindObject("dummy") != None: ROOT.gROOT.FindObject("dummy").Delete()
            histo = ROOT.TH1D("dummy","dummy",1,0.0,1.0); histo.Sumw2()
            (firstEntry, maxEntries) = self._rangeToProcess(fsplit)
            nev = tree.Draw("0.5>>dummy", cut, "goff", maxEntries, firstEntry)
            self.negativeCheck(histo)
            return [ histo.GetBinContent(1), histo.GetBinError(1), nev ]
        else:
            print("HERE FOR _getYield()")
            quit()
            if self._options.doS2V:
                cut  = scalarToVector(cut)
            (firstEntry, maxEntries) = self._rangeToProcess(fsplit)
            npass = tree.Draw("1",self.adaptExpr(cut,cut=True),"goff", maxEntries, firstEntry);
            return [ npass, sqrt(npass), npass ]
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
        if not self._isdata and self._options.weight and not self._options.sumGenWeightFromHisto and len(self._options.maxGenWeightProc) > 0:
            sumGenWeights = self._sumGenWeights.GetValue()
            print(f"Process {self._name} ({self._cname}) -> sum gewWeights = {sumGenWeights}")
            for histo in rets:
                # this is going to trigger the loop, but now all histograms exist so it is ok
                histo.Scale(1./sumGenWeights)
                self.negativeCheck(histo)
                histo.SetDirectory(0)

                
        # fold overflow
        for iret,ret in enumerate(rets):
            if ret.ClassName() in [ "TH1F", "TH1D" ] :
                n = ret.GetNbinsX()
                if plotspecs[iret].getOption('IncludeOverflows',True) and ("TProfile" not in ret.ClassName()):
                    ret.SetBinContent(1,ret.GetBinContent(0)+ret.GetBinContent(1))
                    ret.SetBinContent(n,ret.GetBinContent(n+1)+ret.GetBinContent(n))
                    ret.SetBinError(1,hypot(ret.GetBinError(0),ret.GetBinError(1)))
                    ret.SetBinError(n,hypot(ret.GetBinError(n+1),ret.GetBinError(n)))
                    ret.SetBinContent(0,0)
                    ret.SetBinContent(n+1,0)
                    ret.SetBinContent(0,0)
                    ret.SetBinContent(n+1,0)
                if plotspecs[iret].getOption('IncludeOverflow',False) and ("TProfile" not in ret.ClassName()):
                    ret.SetBinContent(n,ret.GetBinContent(n+1)+ret.GetBinContent(n))
                    ret.SetBinError(n,hypot(ret.GetBinError(n+1),ret.GetBinError(n)))
                    ret.SetBinContent(n+1,0)
                    ret.SetBinContent(n+1,0)
                if plotspecs[iret].getOption('IncludeUnderflow',False) and ("TProfile" not in ret.ClassName()):
                    ret.SetBinContent(1,ret.GetBinContent(0)+ret.GetBinContent(1))
                    ret.SetBinError(1,hypot(ret.GetBinError(0),ret.GetBinError(1)))
                    ret.SetBinContent(0,0)
                    ret.SetBinContent(0,0)
                rebin = plotspecs[iret].getOption('rebinFactor',0)
                if plotspecs[iret].bins[0] != "[" and rebin > 1 and n > 5:
                    while n % rebin != 0: rebin -= 1
                    if rebin != 1: ret.Rebin(rebin)
                if plotspecs[iret].getOption('Density',False):
                    for b in range(1,n+1):
                        ret.SetBinContent( b, ret.GetBinContent(b) / ret.GetXaxis().GetBinWidth(b) )
                        ret.SetBinError(   b, ret.GetBinError(b) / ret.GetXaxis().GetBinWidth(b) )
            self._stylePlot(ret,plotspecs[iret])
            ret._cname = self._cname
        return rets
    # apparently not used anywhere
    def getWeightForCut(self,cut):
        if self._weight:
            if self._isdata: cut = "(%s)     *(%s)*(%s)" % (self._weightString,                    self._scaleFactor, self.adaptExpr(cut,cut=True))
            else:            cut = "(%s)*(%s)*(%s)*(%s)" % (self._weightString,self._options.lumi, self._scaleFactor, self.adaptExpr(cut,cut=True))
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
        if not self._isdata and self._options.weight and not self._options.sumGenWeightFromHisto and len(self._options.maxGenWeightProc) > 0:
            self._sumGenWeights = self._tree.Sum("myweight") # uses floats
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
                    
        # defines
        for name, entry in self._rdfDefs.items():
            if re.match(entry["procRegexp"], self._cname):
                logging.debug("Defining %s as %s" % (name, entry["define"]))
                self._tree = self._tree.Define(name, entry["define"])
        # aliases
        for name,entry in self._rdfAlias.items():
            if re.match(entry["procRegexp"], self._cname):
                self._tree = self._tree.Alias(name, entry["define"])
            
        if self._options.printYieldsRDF:
            logging.debug('='*30)
            logging.debug(" Selection after defines")
            logging.debug('-'*30)
            for (cutname,cutexpr,extra) in cut.cuts():
                if "BEFOREDEFINE" not in extra:
                    logging.debug(f" {cutname}: {cutexpr}")
                    self._tree = self._tree.Filter(cutexpr,cutname)
            logging.debug('='*30)
        else:
            fullcut = cut
            self._tree = self._tree.Filter(fullcut)

        # this column must be defined according to the process name, which is expected to include pre or postVFP
        # the eraVFP column will be used in calls to some functions to specify the era (e.g. for scale factors)
        if "preVFP" in self._cname:
            self._tree = self._tree.Define("eraVFP", "BToF")
        elif "postVFP" in self._cname:
            self._tree = self._tree.Define("eraVFP", "GToH")
        else:
            self._tree = self._tree.Define("eraVFP", "BToH")
            
        # do not call it, or it will trigger the loop now
        #print "sumGenWeights"
        #print sumGenWeights.GetValue()
        histos = []

        # keep list of defined columns used as variables to fill histograms, to avoid defining the same
        # expression multiple times, which may affect performances
        column_expr = {}

        ## weight has a common part for each plot
        if self._weight:
            if self._isdata: wgtCommon = "(%s)     *(%s)" % (self._weightString,                    self._scaleFactor)           
            else:            wgtCommon = "(%s)*(%s)*(%s)" % (self._weightString,self._options.lumi, self._scaleFactor)
        else:
            wgtCommon = '1.' ## wtf is that marc

        
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
                
            (tmp_weight, self._rdfDefsWeightColumns) = self.defineColumnFromExpression(self._rdfDefsWeightColumns, tmp_weight, wgt)
            
            logging.debug("%s weight string = %s " % (tmp_weight, wgt))
            
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
    def processEvents(self,eventLoop,cut):
        logging.info("In processEvents in tree2yield.py")
        if not self._isInit: self._init()
        cut = self.adaptExpr(cut,cut=True)
        if self._options.doS2V:
            cut  = scalarToVector(cut)
            self._tree.vectorTree = True 
        eventLoop.beginComponent(self)
        eventLoop.loop(self._tree, getattr(self._options, 'maxEvents', -1), cut=cut)
        eventLoop.endComponent(self)
    def applyCutAndElist(self,cut,elist):
        if self._appliedCut != None and self._appliedCut != cut: 
            logging.warning("changing applied cut from %s to %s\n" % (self._appliedCut, cut))
        self._appliedCut = cut
        self._elist = elist
    # no longer needed with RDF
    # def cutToElist(self,cut,fsplit=None):
    #     logging.info('beginning of cuttoelist')
    #     if not self._isInit: self._init()
    #     logging.info('afetr init')
    #     ##marcif self._weight:
    #     ##marc    if self._isdata: cut = "(%s)     *(%s)*(%s)" % (self._weightString,                    self._scaleFactor, self.adaptExpr(cut,cut=True))
    #     ##marc    else:            cut = "(%s)*(%s)*(%s)*(%s)" % (self._weightString,self._options.lumi, self._scaleFactor, self.adaptExpr(cut,cut=True))
    #     ##marcelse: cut = self.adaptExpr(cut,cut=True)
    #     cut = self.adaptExpr(cut,cut=True)
    #     if self._options.doS2V: cut  = scalarToVector(cut)
    #     (firstEntry, maxEntries) = self._rangeToProcess(fsplit)
    #     logging.debug('now gonna perform the draw command 4')
    #     ## marc self._tree.Draw('>>elist', cut, 'entrylist', maxEntries, firstEntry)
    #     ## apply a cut, rdf style
    #     logging.debug(' i am filtering with', cut)
    #     self._tree = self._tree.Filter(cut) #Draw('>>elist', cut, 'entrylist', maxEntries, firstEntry)
    #     ## marc elist = ROOT.gDirectory.Get('elist')
    #     ## marc if self._tree.GetEntries()==0 and elist==None: elist = ROOT.TEntryList("elist",cut) # empty list if tree is empty, elist would be a ROOT.nullptr TObject otherwise
    #     ## marc return elist
    def _rangeToProcess(self,fsplit):
        if fsplit != None and fsplit != (0,1):
            if self._options.maxEntriesNotData and self._isdata:
                allEntries = self.getEntries()
            else:
                allEntries = min(self.getEntries(), self._options.maxEntries)
            chunkSize = int(ceil(allEntries/float(fsplit[1])))
            firstEntry = chunkSize * fsplit[0]
            maxEntries = chunkSize # the last chunk may go beyond the end of the tree, but ROOT stops anyway so we don't care
        else:
            firstEntry = 0
            if (self._options.maxEntriesNotData and self._isdata):
                maxEntries = self.getEntries() 
            else:
                maxEntries = self._options.maxEntries 
        return (firstEntry, maxEntries)
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
    parser.add_argument("-l", "--lumi", type=float, default="1", help="Luminosity (in 1/fb)");
    parser.add_argument("-u", "--unweight", dest="weight", action="store_false", help="Don't use weights (in MC events), note weights are still used if a fake rate file is given");
    parser.add_argument("--uf", "--unweight-forced", dest="forceunweight", action="store_true", help="Do not use weight even if a fake rate file is given.");
    parser.add_argument("-W", "--weightString", dest="weightString", action="append", default=[], help="Use weight (in MC events), can specify multiple times");
    parser.add_argument("-f", "--final", action="store_true", help="Just compute final yield after all cuts");
    parser.add_argument("-e", "--errors", action="store_true", help="Include uncertainties in the reports");
    parser.add_argument("--tf", "--text-format", dest="txtfmt", type=str, default="txt", choices=["txt","tsv","csv","dsv","ssv"], help="Output format: txt,tsv,csv,dsv,ssv");
    parser.add_argument("-S", "--start-at-cut", dest="startCut", type=str, help="Run selection starting at the cut matched by this regexp, included.") 
    parser.add_argument("-U", "--up-to-cut", dest="upToCut", type=str, help="Run selection only up to the cut matched by this regexp, included.") 
    parser.add_argument("-X", "--exclude-cut", dest="cutsToExclude", action="append", default=[], help="Cuts to exclude (regexp matching cut name), can specify multiple times.") 
    parser.add_argument("-E", "--enable-cut", dest="cutsToEnable", action="append", default=[], help="Cuts to enable if they were disabled in the cut file (regexp matching cut name), can specify multiple times.") 
    parser.add_argument("-I", "--invert-cut", dest="cutsToInvert", action="append", default=[], help="Cuts to invert (regexp matching cut name), can specify multiple times.") 
    parser.add_argument("-R", "--replace-cut", dest="cutsToReplace", action="append", default=[], nargs=3, help="Cuts to invert (regexp of old cut name, new name, new cut); can specify multiple times.") 
    parser.add_argument("-A", "--add-cut", dest="cutsToAdd", action="append", default=[], nargs=3, help="Cuts to insert (regexp of cut name after which this cut should go, new name, new cut); can specify multiple times.") 
    parser.add_argument("-N", "--n-minus-one", dest="nMinusOne", action="store_true", help="Compute n-minus-one yields and plots")
    parser.add_argument("--select-n-minus-one", dest="nMinusOneSelection", type=str, help="Select which cuts to do N-1 for (comma separated list of regexps)")
    parser.add_argument("--NI", "--inv-n-minus-one", dest="nMinusOneInverted", action="store_true", help="Compute n-minus-one yields and plots")
    parser.add_argument("--obj", "--objname", dest="obj", default='Events', help="Pattern for the name of the TTree inside the file");
    parser.add_argument("--RV", "--replace-var", dest="varsToReplace", action="append", default=[], nargs=2, help="Variable to replace(old var name, new var name); can specify multiple times.")
    parser.add_argument("-G", "--no-fractions", dest="fractions", action="store_false", default=True, help="Don't print the fractions");
    parser.add_argument("-F", "--add-friend", dest="friendTrees", action="append", default=[], nargs=2, help="Add a friend tree (treename, filename). Can use {name}, {cname} patterns in the treename") 
    parser.add_argument("--Fs", "--add-friend-simple", dest="friendTreesSimple", action="append", default=[], nargs=1, help="Add friends in a directory. The rootfile must be called evVarFriend_{cname}.root and tree must be called 't' in a subdir 'sf' inside the rootfile.") 
    parser.add_argument("--FMC", "--add-friend-mc", dest="friendTreesMC", action="append", default=[], nargs=2, help="Add a friend tree (treename, filename) to MC only. Can use {name}, {cname} patterns in the treename") 
    parser.add_argument("--FD", "--add-friend-data", dest="friendTreesData", action="append", default=[], nargs=2, help="Add a friend tree (treename, filename) to data trees only. Can use {name}, {cname} patterns in the treename") 
    parser.add_argument("--FMCs", "--add-friend-mc-simple", dest="friendTreesMCSimple", action="append", default=[], nargs=1, help="Add friends in a directory to MC only. The rootfile must be called evVarFriend_{cname}.root and tree must be called 't' in a subdir 'sf' inside the rootfile.") 
    parser.add_argument("--FDs", "--add-friend-data-simple", dest="friendTreesDataSimple", action="append", default=[], nargs=1, help="Add friends in a directory to data only. The rootfile must be called evVarFriend_{cname}.root and tree must be called 't' in a subdir 'sf' inside the rootfile.") 
    parser.add_argument("--mcc", "--mc-corrections", dest="mcCorrs", action="append", default=[], nargs=1, help="Load the following file of mc to data corrections") 
    parser.add_argument("--s2v", "--scalar2vector", dest="doS2V", action="store_true", help="Do scalar to vector conversion") 
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

