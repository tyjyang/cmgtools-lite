#!/usr/bin/env python
from math import *
import re
import os, os.path
from array import array
from copy import *
import logging

class CutsFile:
    def __init__(self,txtfileOrCuts,options=None,ignoreEmptyOptionsEnforcement=False):
        if not ignoreEmptyOptionsEnforcement:
            if not options: raise RuntimeError('options is empty when passed to CutsFile constructor')
        if type(txtfileOrCuts) == list:
            self._cuts = deepcopy(txtfileOrCuts[:])
        elif isinstance(txtfileOrCuts,CutsFile):
            self._cuts = deepcopy(txtfileOrCuts.cuts())
        else:
            self._cuts = []
            file = open(txtfileOrCuts, "r")
            if not file: raise RuntimeError("Cannot open "+txtfileOrCuts+"\n")
            if options:
                for cr,cn,cv in options.cutsToAdd:
                    if re.match(cr,"entry point"): self._cuts.append((cn,cv,{}))
            for line in file:
              try:
                line = line.strip()
                if len(line) == 0 or line[0] == '#': continue
                line = re.sub(r"(?<!\\)#.*","",line)  ## regexp black magic: match a # only if not preceded by a \!
                line = line.replace(r"\#","#")        ## and now we just unescape the remaining #'s
                while line[-1] == "\\":
                    line = line[:-1] + " " + file.next().strip()
                    line = re.sub(r"(?<!\\)#.*","",line)  ## regexp black magic: match a # only if not preceded by a \!
                    line = line.replace(r"\#","#")        ## and now we just unescape the remaining #'s
                extra = {}
                if ";" in line:
                    (line, extra) = self.getExtraFromLine(line)
                (name,cut) = [x.strip().replace(";",":") for x in line.replace("\:",";").split(":")]
                if name == "entry point" and cut == "1": continue
                if options:
                    for cn in options.cutsToEnable:
                        if re.search(cn,name): extra["Disable"] = False 
                if "Disable" in extra and extra["Disable"] == True: continue
                if options:
                    if options.startCut and not re.search(options.startCut,name): continue
                    if options.startCut and re.search(options.startCut,name): options.startCut = None
                self._cuts.append((name,cut,extra))
                if options:
                    for cr,cn,cv in options.cutsToAdd:
                        if re.match(cr,name):
                            if ";" in cv:
                                (actualNewCut, newExtra) = self.getExtraFromLine(cv)
                            else:
                                (actualNewCut, newExtra) = (cv, {})
                            self._cuts.append((cn,actualNewCut,newExtra))
                    if options.upToCut and re.search(options.upToCut,name):
                        break
              except ValueError(e):
                logging.error("Failed to parse cut line [%s]" % line.strip())
                raise 
            if options:
                for ci in options.cutsToInvert:  self.invert(ci)
                for ci in options.cutsToExclude: self.remove(ci)
                for cr,cn,cv in options.cutsToReplace: self.replace(cr,cn,cv)
                for vr,vn in options.varsToReplace: self.replaceCutVariable(vr,vn)
    def __str__(self):
        newstring = ""
        for cut in self._cuts:
            newstring += "{0} : {1};{2}\n".format(cut[0],cut[1],", ".join(f"{k}={cut[2][k]}" for k in list(cut[2].keys())))
        return newstring[:-1]

    def getExtraFromLine(self, line):
        retExtra = {}
        if ";" not in line:
            return (line, retExtra)
        (retline,more) = line.split(";")[:2]
        for setting in [f.replace(';',',').strip() for f in more.replace('\\,',';').split(',')]:
            if "=" in setting:
                (key,val) = [f.strip() for f in setting.split("=")]
                retExtra[key] = eval(val)# if val in ["True", "False"] else str(val) 
            else: retExtra[setting] = True
        return (retline, retExtra)
    
    def remove(self,cut):
        self._cuts = [(cn,cv,extra) for (cn,cv,extra) in self._cuts if not re.search(cut,cn)]
        return self
    def invert(self,cut):
        for i,(cn,cv,extra) in enumerate(self._cuts[:]):
            if re.search(cut,cn):
                if cn.startswith("not ") and re.match(r"!\(.*\)", cv):
                    self._cuts[i] = (cn[4:], cv[2:-1])
                else:
                    self._cuts[i] = ("not "+cn, "!("+cv+")")
        return self
    def replace(self,cut,newname,newcut):       
        for i,(cn,cv,extra) in enumerate(self._cuts[:]):
            if re.search(cut,cn):
                if ";" in newcut:
                    (actualNewCut, newExtra) = self.getExtraFromLine(newcut)
                else:
                    (actualNewCut, newExtra) = (newcut, {})
                self._cuts[i] = (newname, actualNewCut, newExtra)
        return self
    def cuts(self):
        return self._cuts[:]
    def sequentialCuts(self):
        # for simplicity this just keep the extra of the line being added, which might not be what one wants
        # but probably for sequential cuts those extra would not even be present
        if len(self._cuts) == 0: return []
        ret = [ (self._cuts[0][0], "(%s)" % self._cuts[0][1], self._cuts[0][2]) ]
        for (cn,cv,extra) in self._cuts[1:]:
            ret.append( ( cn, "%s && (%s)" % (ret[-1][1], cv), extra ) )
        return ret
    def nMinusOne(self,inverted=False):
        return CutsFile(self.nMinusOneCuts(inverted=inverted))
    def nMinusOneCuts(self,inverted=False):
        ret = []
        for cn,cv,extra in self._cuts[1:]:
            nm1 = " && ".join("(%s)" % cv1 for cn1,cv1,extra in self._cuts if cn1 != cn)
            if inverted:
                ret.append(("fail only "+cn, "%s && !(%s)" % (nm1,cv,extra)))
            else:
                ret.append(("all but "+cn, nm1,extra))
        return ret
    def nMinusOneSelectedCuts(self,selection,inverted=False):
        ret = []
        cutmatches = [ re.compile(c.strip()) for c in selection.split(",") if c.strip() ]
        for cn,cv in self._cuts[1:]:
            if any(m for m in cutmatches if re.match(m,cn)):
                nm1 = " && ".join("(%s)" % cv1 for cn1,cv1 in self._cuts if cn1 != cn)
                ret.append(("all but "+cn, nm1))
        return ret
    def allCuts(self,n=-1,doProduct=False):
        jstring = " * " if doProduct else " && "
        return jstring.join("(%s)" % x[1] for x in (self._cuts[0:n+1] if n != -1 and n+1 < len(self._cuts) else self._cuts))
    def addAfter(self,cut,newname,newcut):
        for i,(cn,cv,extra) in enumerate(self._cuts[:]):
            if re.search(cut,cn):
                if ";" in newcut:
                    (actualNewCut, newExtra) = self.getExtraFromLine(newcut)
                else:
                    (actualNewCut, newExtra) = (newcut, {})
                self._cuts.insert(i+1,(newname, actualNewCut, newExtra))
                break
        return self
    def insert(self,index,newname,newcut, extra={}):
        self._cuts.insert(index,(newname, newcut, extra))
        return self
    def add(self,newname,newcut, extra={}):
        self._cuts.append((newname,newcut,extra))
        return self
    def setParams(self,paramMap):
        self._cuts = [ (cn.format(**paramMap), cv.format(**paramMap), extra.format(**paramMap)) for (cn,cv,extra) in self._cuts ]
    def cartesianProduct(self,other):
        return CutsFile( [ ("%s && %s" % (cn1,cn2), "(%s) && (%s)" % (cv1,cv2)) for (cn1,cv1) in self._cuts for (cn2,cv2) in other.cuts() ] )
    def replaceCutVariable(self,oldvar,newvar):
        for i,(cn,cv,extra) in enumerate(self._cuts[:]):
            if re.search(oldvar,cv):
                if ";" in newvar:
                    (actualNewVar, newExtra) = self.getExtraFromLine(newvar)
                else:
                    (actualNewVar, newExtra) = (newcut, {})
                newcut = cv.replace(oldvar,actualNewVar)
                self._cuts[i] = (cn, newcut, newExtra)

