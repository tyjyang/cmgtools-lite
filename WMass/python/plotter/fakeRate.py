import re
import os
import ROOT
import logging

def compileMacro(x,basedir=os.environ['PWD']):
    #ROOT.gROOT.ProcessLine(".L %s/%s+" % (os.environ['CMSSW_BASE'],x));
    #success = ROOT.gSystem.CompileMacro("%s/%s" % (os.environ['PWD'],x),"k")
    success = ROOT.gSystem.CompileMacro("%s" % (x),"k")
    if not success:
        logging.error("Loading and compiling %s failed! Exit" % x)
        quit()

if "/fakeRate_cc.so" not in ROOT.gSystem.GetLibraries(): 
   compileMacro("ccFiles/fakeRate.cc")
   #compileMacro("fakeRate.cc")

##marc from CMGTools.TTHAnalysis.plotter.mcCorrections import SimpleCorrection
from mcCorrections import SimpleCorrection
#from CMGTools.WMass.plotter.mcCorrections import SimpleCorrection
##marc from CMGTools.WMass.plotter.cutsFile import *
from cutsFile import *


_loads = {}
class FakeRate:
    def __init__(self,filestring,lumi=None):
        files = filestring.split(",")
        self._weight = None
        self._mods = []
        self._cutMods = []
        for file in files:
            stream = open(file,'r')
        for line in stream:
            if len(line.strip()) == 0 or line.strip()[0] == '#': continue
            while line.strip()[-1] == "\\":
                line = line.strip()[:-1] + stream.next()
            fields = [x.strip() for x in line.split(":")]
            if fields[0] == "weight":
                if self._weight is not None: raise RuntimeError("Duplicate weight definition in fake rate file "+file)
                self._weight = fields[1]
            elif fields[0] == "change": 
                self._mods.append( SimpleCorrection(fields[1],fields[2]) )
            elif fields[0] == "cut-change": 
                self._cutMods.append( SimpleCorrection(fields[1],fields[2],onlyForCuts=True) )
            elif fields[0].startswith("load-histo"):
                ## marc data = "%s/src/CMGTools/WMass/data/" % os.environ['CMSSW_BASE'];
                data = "../../data/"
                fname = fields[2].replace("$DATA",data)
                hname = fields[3] if len(fields) >= 4 else fields[1]
                if fields[1] in _loads:
                    if _loads[fields[1]] != (fname,hname):
                        #print "Conflicting load for %s: %s, %s vs older %s " % (fields[1],fields[2],hname, _loads[fields[1]])
                        pass
                    else:
                        #print "Duplicate load for %s: %s, %s" % (fields[1],fields[2],hname)
                        pass
                else:
                    _loads[fields[1]] = (fname,hname)
                if not '3D' in fields[0]:
                    ROOT.loadFRHisto(fields[1],fname,hname)
                else:
                    ROOT.loadFRHisto3D(fields[1],fname,hname)
            elif fields[0] == 'norm-lumi-override':
                if self._weight is None: raise RuntimeError("norm-lumi-override must follow weight declaration in fake rate file "+file)
                if not lumi: raise RuntimeError("lumi not set in options, cannot apply norm-lumi-override")
                logging.warning("normalization overridden from %s/fb to %s/fb in fake rate file %s" % (lumi,fields[1],file))
                self._weight = '((%s)*(%s)/(%s))' % (self._weight,fields[1],lumi)
            elif fields[0] == 'cut-file':
                if self._weight is None: raise RuntimeError("cut-file must follow weight declaration in fake rate file "+file)
                addcuts = CutsFile(fields[1],options=None,ignoreEmptyOptionsEnforcement=True)
                self._weight = '((%s)*(%s))' % (self._weight,addcuts.allCuts(doProduct=True))
#                    print "WARNING: cuts loaded from fake rate file "+file
            else:
                raise RuntimeError("Unknown directive "+fields[0])
        if file==files[0]:
            if self._weight is None: raise RuntimeError("Missing weight definition in fake rate file "+file)
#        if len(self._cutMods) == 0: print "WARNING: no directives to change cuts in fake rate file "+filestring
    def weight(self): 
        return self._weight
    def mods(self): 
        return self._mods
    def cutMods(self): 
        return self._cutMods
