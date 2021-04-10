#!/usr/bin/env python

from tree2yield import *
from projections import *
from figuresOfMerit import FOM_BY_NAME
import pickle, re, random, time, glob, math
import argparse
import json

ROOT.ROOT.EnableImplicitMT()
#ROOT.ROOT.TTreeProcessorMT.SetMaxTasksPerFilePerWorker(1)

#_T0 = long(ROOT.gSystem.Now())

## These must be defined as standalone functions, to allow running them in parallel
## these might actually be no longer needed with RDF
def _runYields(args):
    key,tty,cuts,noEntryLine,fsplit = args
    return (key, tty.getYields(cuts,noEntryLine=noEntryLine,fsplit=fsplit))

def _runPlot(args):
    key,tty,plotspec,cut,fsplit,closeTree = args
    timer = ROOT.TStopwatch()
    #print "Starting plot %s for %s, %s" % (plotspec.name,key,tty._cname)
    ret = (key,tty.getPlot(plotspec,cut,fsplit=fsplit,closeTreeAfter=closeTree))
    #print "Done plot %s for %s, %s, fsplit %s in %s s, at %.2f; entries = %d, time/entry = %.3f ms" % (plotspec.name,key,tty._cname,fsplit,timer.RealTime(), 0.001*(long(ROOT.gSystem.Now()) - _T0), ret[1].GetEntries(), (long(ROOT.gSystem.Now()) - _T0)/float(ret[1].GetEntries()))
    return ret

# def _runApplyCut(args):
#     key,tty,cut,fsplit = args
#     return (key, tty.cutToElist(cut,fsplit=fsplit))

def _runGetEntries(args):
    key,tty = args
    return (key, tty.getEntries())

def initializeJson(jsonMap, jsonFile):
    if not os.path.isfile(jsonFile):
         raise ValueError("Could not find json file %s" % jsonFile)
    info = json.load(open(jsonFile))

    for k,v in info.items():
        vec = ROOT.std.vector["std::pair<unsigned int, unsigned int>"]()
        for combo in v:
            pair = ROOT.std.pair["unsigned int", "unsigned int"](*[int(c) for c in combo])
            vec.push_back(pair)
        jsonMap[int(k)] = vec

class MCAnalysis:
    def __init__(self,samples,options):
        self._options = options
        self._allData     = {}
        self._data        = []
        self._signals     = []
        self._backgrounds = [] 
        self._isSignal    = {}
        self._rank        = {} ## keep ranks as in the input text file
        self._projection  = Projections(options.project, options) if options.project != None else None
        self._premap = []
        self._optionsOnlyProcesses = {}
        self.init_defaults = {}
        for premap in options.premap:
            to,fro = premap.split("=")
            if to[-1] == ":": to = to[:-1]
            to = to.strip()
            for k in fro.split(","):
                self._premap.append((re.compile(k.strip()+"$"), to))
        if hasattr(ROOT, "initializeScaleFactors"):
            logging.info("Initializing histograms with scale factors")
            ROOT.initializeScaleFactors(self._options.scaleFactorFile)
        if hasattr(ROOT, "jsonMap_all"):
            initializeJson(ROOT.jsonMap_all, self._options.json)
            logging.info(f"Initialized json files for data: {self._options.json}")
        self.allRDF = {}    # dictionary with rdf per process (the key is the process cname)
        self.allHistos = {} # dictionary with a key for process name and a list of histograms (or dictionary to have a name for each histogram)
        self.allCutReports = {}
        self.readMca(samples,options)
     
    def getSumGenWeightMCfromHisto(self, pname, rootfile, verbose=False):

        maxGenWgt = None
        sumGenWeights = 1.0
        nUnweightedEvents = 1.0

        ROOT.gEnv.SetValue("TFile.AsyncReading", 1);
        #tmp_rootfile = ROOT.TXNetFile(rootfile+"?readaheadsz=65535")
        tmp_rootfile = ROOT.TFile.Open(rootfile)
 
        histo_sumgenweight = tmp_rootfile.Get('hGenWeights')
        if not histo_sumgenweight:
            raise RuntimeError("Can't get histogram hGenWeights from %s.\nMake sure it is available when using option --max-genWeight-procs" % rootfile)

        if self._options.weight and len(self._options.maxGenWeightProc):
            # get sum of weights from histograms, filtering some events with large weights
            for procRegExp,tmp_maxGenWgt in self._options.maxGenWeightProc:
                if re.match(procRegExp,pname):
                    #tmp_maxGenWgt is a string, convert to float
                    maxGenWgt = fabs(float(tmp_maxGenWgt))
                    log10_maxGenWgt = math.log10(maxGenWgt)
                    minBin = histo_sumgenweight.GetXaxis().FindFixBin(-log10_maxGenWgt)
                    maxBin = histo_sumgenweight.GetXaxis().FindFixBin(log10_maxGenWgt)
                    sumGenWeights = histo_sumgenweight.Integral(minBin,maxBin)                  
                    if self._options.clipGenWeightToMax:
                        # then the other histogram, whose integral is the number of events before any cut
                        histo_numweight = tmp_rootfile.Get('hNumWeights')
                        if not histo_numweight:
                            raise RuntimeError("Can't get histogram hNumWeights from %s.\nMake sure it is available when using option --max-genWeight-procs and --clip-genWeight-toMax" % rootfile)
                        # get actual upper threshold based on the bin edge
                        maxGenWgt = math.pow(10.0,histo_sumgenweight.GetXaxis().GetBinUpEdge(maxBin))
                        if verbose:
                            logging.info("Process %s -> actual genWeight threshold set to %s" % (pname,str(maxGenWgt)))
                        nLargeWeight = -1*histo_numweight.Integral(0,max(0,minBin-1)) + histo_numweight.Integral(min(histo_numweight.GetNbinsX()+1,maxBin+1), histo_numweight.GetNbinsX()+1)
                        sumGenWeights += (nLargeWeight*maxGenWgt) 
                        nUnweightedEvents = histo_numweight.Integral(0,histo_numweight.GetNbinsX()+1)
        else:
            sumGenWeights = histo_sumgenweight.Integral(0,1+histo_sumgenweight.GetNbinsX())
            nUnweightedEvents = histo_sumgenweight.GetEntries()            
            
        tmp_rootfile.Close()                            
        return (maxGenWgt, sumGenWeights, nUnweightedEvents)
 

    def readMca(self,samples,options,addExtras={},field0_addExtras=""):

        # when called with not empty addExtras, issue a warning in case you are overwriting settings
        if "ALLOW_OVERWRITE_SETTINGS" in addExtras:
            if addExtras["ALLOW_OVERWRITE_SETTINGS"] == True:
                logging.warning("Found setting ALLOW_OVERWRITE_SETTINGS=True in %s. " % str(field0_addExtras))
                loggin.warning("Will use new settings overwriting those in mca included file %s." % str(samples))
                
        for line in open(samples,'r'):
            if re.match("\s*#.*", line): continue
            line = re.sub(r"(?<!\\)#.*","",line)  ## regexp black magic: match a # only if not preceded by a \!
            line = line.replace(r"\#","#")        ## and now we just unescape the remaining #'s
            extra = {}
            if ";" in line:
                (line,more) = line.split(";")[:2]
                for setting in [f.replace(';',',').strip() for f in more.replace('\\,',';').split(',')]:
                    if setting == "": continue
                    if "=" in setting: 
                        (key,val) = [f.strip() for f in setting.split("=",1)]
                        extra[key] = eval(val)
                    else: extra[setting] = True
            for k,v in addExtras.items():
                if k in extra: 
                    if "ALLOW_OVERWRITE_SETTINGS" in addExtras and addExtras["ALLOW_OVERWRITE_SETTINGS"] == True:
                        pass
                    else:
                        raise RuntimeError('You are trying to overwrite an extra option (%s - %s ) already set (did you forget ALLOW_OVERWRITE_SETTINGS=True ?)' % (k, v))
                extra[k] = v
            field = [f.strip() for f in line.split(':')]
            if len(field) == 1 and field[0] == "*":
                if len(self._allData): raise RuntimeError("MCA defaults ('*') can be specified only before all processes")
                #print "Setting the following defaults for all samples: "
                for k,v in extra.items():
                    #print "\t%s: %r" % (k,v)
                    self.init_defaults[k] = v
                continue
            else:
                for k,v in self.init_defaults.items():
                    if k not in extra: extra[k] = v
            if len(field) <= 1: continue
            if "SkipMe" in extra and extra["SkipMe"] == True and not options.allProcesses: continue
            field0noPostFix = field[0]
            if 'PostFix' in extra:
                hasPlus = (field[0][-1]=='+')
                if hasPlus: field[0] = field[0][:-1]
                field[0] += extra['PostFix']
                if hasPlus: field[0]+='+'
            signal = False
            pname = field[0]
            if pname[-1] == "+": 
                signal = True
                pname = pname[:-1]
            ## if we remap process names, do it
            for x,newname in self._premap:
                if re.match(x,pname):
                    pname = newname
           ## If we have a user-defined list of processes as signal
            if len(options.processesAsSignal):
                signal = False
                for p0 in options.processesAsSignal:
                    for p in p0.split(","):
                        if re.match(p+"$", pname): signal = True
            ## Options only processes
            if field[1] == "-": 
                self._optionsOnlyProcesses[field[0]] = extra
                self._isSignal[field[0]] = signal
                continue
            if field[1] == "+": # include an mca into another one, usage:   otherprocesses : + ; IncludeMca="path/to/other/mca.txt"
                if 'IncludeMca' not in extra: raise RuntimeError('You have declared a component with IncludeMca format, but not included this option')
                extra_to_pass = copy(extra)
                del extra_to_pass['IncludeMca']
                self.readMca(extra['IncludeMca'],options,addExtras=extra_to_pass,field0_addExtras=field0noPostFix) # call readMca recursively on included mca files
                continue
            # Customize with additional weight if requested
            if 'AddWeight' in extra:
                if len(field)<2: raise RuntimeError('You are trying to set an additional weight, but there is no weight initially defined for this component')
                elif len(field)==2: field.append(extra['AddWeight'])
                else: field[2] = '(%s)*(%s)'%(field[2],extra['AddWeight'])
            ## If we have a selection of process names, apply it
            skipMe = (len(options.processes) > 0)
            for p0 in options.processes:
                for p in p0.split(","):
                    if re.match(p+"$", pname): skipMe = False
            for p0 in options.processesToExclude:
                for p in p0.split(","):
                    if re.match(p+"$", pname): skipMe = True
            for p0 in options.filesToExclude:
                for p in p0.split(","):
                    if re.match(p+"$", field[1]): skipMe = True
            if skipMe: continue
            cnamesWithPath = []
            
            #print ">>> printing line"
            #print line
            #print ">>>"
            if options.nanoaodTree:    
                field1_regexp = re.compile(field[1])
                if len(options.filterProcessFiles):
                    for procRegexp,filterRegExp in options.filterProcessFiles:
                        if re.match(procRegexp,field[0]): 
                            field1_regexp = re.compile(filterRegExp)
                                
                subpath = extra['SubPath'] if 'SubPath' in extra else ".*"
                subPath_regexp = re.compile(subpath)
                pathsToSearch = options.path
                if 'TreePath' in extra:
                    if '[XXX]' in extra['TreePath']:                        
                        if 'XXX' in extra:
                            XXX = extra['XXX'].split(',')
                            pathsToSearch = [extra['TreePath'].replace("[XXX]", x) for x in XXX]
                        else:
                            logging.error("Process %s '[XXX]' found in TreePath, but XXX not defined in MCA file" % (pname))
                            quit()
                    else:                        
                        pathsToSearch = [extra['TreePath']]
                else:
                    if not options.path:
                        logging.warning("You didn't specify a path to ntuples with option -P, but process %s has no 'TreePath' key in the MCA file. Please specify a valid path." % pname)
                        quit()
                logging.info("Process %s -> searching for files in %s/%s" % (pname,repr(pathsToSearch),subpath))
                for treepath in pathsToSearch:
                    for dirpath, dirnames, filenames in os.walk(treepath):
                        # either use regexp or simply that subpath is in the path (so one doesn't need 
                        # to start subpath with .* for a proper match if a word without wildcards is used)
                        if subPath_regexp.match(dirpath+"/") or subpath in (dirpath+"/"):
                            filtered_fnames = [f for f in filenames if f.endswith(".root") and field1_regexp.match(f)]
                            for filename in filtered_fnames:
                                cnamesWithPath.append(os.path.join(dirpath, filename))
                logging.info("Process %s -> %d files selected" % (pname,len(cnamesWithPath)))
            else:
                match=re.match("(\S+)\*",field[1]) 
                if match:
                    for treepath in options.path:
                        if 'SubPath' in extra:
                            cnamesWithPath += (list(glob.iglob("%s/%s/%s*" % (treepath,extra['SubPath'],match.group(0)))))
                        else:
                            cnamesWithPath += (list(glob.iglob("%s/%s*" % (treepath,match.group(0)))))
                    cnames = [os.path.basename(cname) for cname in cnamesWithPath]
                else: cnames = [ x.strip() for x in field[1].split("+") ]

            total_w = 0.; to_norm = False; ttys = [];
            is_w = -1
            pname0 = pname

            if options.weight and len(options.maxGenWeightProc):
                for procRegExp,tmp_maxGenWgt in options.maxGenWeightProc:
                    if re.match(procRegExp,pname):
                        if options.clipGenWeightToMax:
                            logging.info("Process %s -> clipping genWeight to |x| < %s" % (pname,tmp_maxGenWgt))
                        else:
                            logging.info("Process %s -> rejecting events with genWeight > %s" % (pname,tmp_maxGenWgt))
            
            ## cname = os.path.basename(cnameWithPath)
            if options.useCnames: pname = pname0+"."+cname
            for (ffrom, fto) in options.filesToSwap:
                if cname == ffrom: cname = fto
            treename = extra["TreeName"] if "TreeName" in extra else options.tree 
            objname  = extra["ObjName"]  if "ObjName"  in extra else options.obj
            subpath  = extra["SubPath"]  if "SubPath"  in extra else ""
            friendDir = extra['FriendDir'] if 'FriendDir' in extra else ""
            if '[XXX]' in friendDir:
                if "XXX" in extra:
                    XXX = extra['XXX'].split(',')
                    friendDir_tmp = friendDir
                    for x in XXX:                            
                        if friendDir_tmp.replace("[XXX]",x) in cnameWithPath:
                            friendDir = friendDir.replace("[XXX]",x)
                            # print "FriendDir = ",friendDir
                else:
                    logging.error("Process %s '[XXX]' found in FriendDir, but XXX not defined in MCA file" % (pname))
                    quit()

            if subpath != "" and not subpath.endswith("/"):
                subpath += "/"
            #print "subpath = %s" % subpath

            
            ### CHECKPOINT
            nAllFiles = len(cnamesWithPath)
            tmp_names = ROOT.std.vector('string')()
            for n in cnamesWithPath: 
                tmp_names.push_back(n.replace('/eos/cms/','root://eoscms.cern.ch//'))
                #if '/eos/cms/' in n and not n.startswith("root://"):
                #    tmp_names.push_back('root://eoscms.cern.ch//' + n)
            if field[0] in list(self.allRDF.keys()):
                # case with multiple lines in MCA with same common name and different processes attached to it
                matchedCnames = list(filter(lambda x : re.match(field[0]+"__ext\d+$",x),self.allRDF.keys())) # start with cname followed by __ ext and at least one digit (and ending like that
                cnindex = len(matchedCnames) 
                field[0] = f"{field[0]}__ext{cnindex}" # add number to component to distinguish it from other ones
            if self._options.rdfRange:
                ROOT.ROOT.DisableImplicitMT() # to be called before creating the RDataFrame
            self.allRDF[field[0]] = ROOT.RDataFrame(objname, tmp_names)
            #tmp_rdf = ROOT.RDataFrame(objname, tmp_names)
            if self._options.rdfRange:
                begin,end,stride = self._options.rdfRange
                self.allRDF[field[0]] = self.allRDF[field[0]].Range(begin,end,stride)

            tmp_rdf = self.allRDF[field[0]]            
            tty = TreeToYield(tmp_rdf, options, settings=extra, name=pname, cname=field[0], objname=objname, frienddir=friendDir)
            ttys.append(tty)
            if signal: 
                self._signals.append(tty)
                self._isSignal[pname] = True
            elif pname == "data":
                self._data.append(tty)
            else:
                self._isSignal[pname] = False
                self._backgrounds.append(tty)
            if pname in self._allData: self._allData[pname].append(tty)
            else                     : self._allData[pname] =     [tty]
            if "data" not in pname:
                if options.noHeppyTree:
                    if options.weight:
                        if (is_w==0): raise RuntimeError("Can't put together a weighted and an unweighted component (%s)" % cnames)
                        is_w = 1
                        total_w = 1.0  # not really used for general trees at the moment
                        scale = "(%s)" % field[2]
                    else:
                        scale = "1"
                        total_w = 1 # not really used for general trees at the moment
                        if (is_w==1): raise RuntimeError("Can't put together a weighted and an unweighted component (%s)" % cnames)
                        is_w = 0
                elif options.nanoaodTree:

                    maxGenWgt = None
                    sumGenWeights = 1.0
                    nUnweightedEvents = 1.0
                    if options.sumGenWeightFromHisto:
                        # useful for skims when rejecting events with fancy weights
                        # if not doing such removal, it is still better to use the actual
                        # information in the Runs tree, as done below
                        for iFile,f in enumerate(cnamesWithPath):
                            (maxGenWgt, tmp_sumGenWeights, tmp_nUnweightedEvents) = self.getSumGenWeightMCfromHisto(pname, f, verbose=(iFile==0))
                            sys.stdout.write('INFO >>> preparing files {0:.2%}  \r'.format(float(iFile+1)/nAllFiles))
                            sys.stdout.flush()
                            sumGenWeights += tmp_sumGenWeights
                            nUnweightedEvents += tmp_nUnweightedEvents
                    else:
                        if options.weight and len(options.maxGenWeightProc):
                            # get sum of weights from Events tree, filtering some events with large weights
                            # this assumes the trees are unskimmed to correctly compute the sum!!
                            for procRegExp,tmp_maxGenWgt in options.maxGenWeightProc:
                                if re.match(procRegExp,pname):
                                    #tmp_maxGenWgt is a string, convert to float
                                    maxGenWgt = abs(float(tmp_maxGenWgt))                        
                                    if options.clipGenWeightToMax:
                                        # set weights > max to max
                                        #tmp_rdf = tmp_rdf.Define('myweight', 'std::copysign(1.0,genWeight) * std::min<float>(std::abs(genWeight),{s})'.format(s=maxGenWgt))
                                        tmp_rdf = tmp_rdf.Define('myweight', 'genWeightLargeClipped(genWeight,{s})'.format(s=maxGenWgt))
                                    else:
                                        # reject event with weight > max
                                        #tmp_rdf = tmp_rdf.Define('myweight', 'genWeight*(abs(genWeight) < {s})'.format(s=str(maxGenWgt)))                                      
                                        tmp_rdf = tmp_rdf.Define('myweight', 'genWeightLargeRemoved(genWeight,{s})'.format(s=maxGenWgt))   
                                    # set again tmp_rdf
                                    tty.setRDF(tmp_rdf)
                        elif options.weight:
                            # we immediately trigger the loop here, but each Runs tree has one event
                            # so it should be quite fast to get the sum like this
                            tmp_rdf_runs = ROOT.RDataFrame("Runs", tmp_names)
                            _sumGenWeights = tmp_rdf_runs.Sum('genEventSumw')
                            _nUnweightedEvents = tmp_rdf_runs.Sum('genEventCount')
                            sumGenWeights = _sumGenWeights.GetValue()
                            nUnweightedEvents = _nUnweightedEvents.GetValue()

                    if options.weight and True: # True for now, later on this could explicitly require using the actual genWeights as opposed to using sum of unweighted events for MC (see the case for cmgtools below)
                        if (is_w==0): raise RuntimeError("Can't put together a weighted and an unweighted component (%s)" % cnames)
                        is_w = 1; 
                        if options.sumGenWeightFromHisto or (len(options.maxGenWeightProc) == 0):
                            total_w = sumGenWeights
                        else:
                            # in this case division by sum of gen weights is done at the end by scaling histograms
                            # if we tried to compute the sum now, this would trigger the loop on rdf
                            # but it would be a waste of time
                            total_w = 1.0 
                        if maxGenWgt != None:
                            scale = '(%s)*myweight' % field[2]                            
                        else:
                            scale = "genWeight*(%s)" % field[2]
                    elif not options.weight:
                        scale = "1"
                        total_w = 1
                        if (is_w==1): raise RuntimeError("Can't put together a weighted and an unweighted component (%s)" % cnames)
                        is_w = 0
                    else: # case for unweighted events for MC
                        if (is_w==1): raise RuntimeError("Can't put together a weighted and an unweighted component (%s)" % cnames)
                        is_w = 0;
                        total_w = nUnweightedEvents
                        scale = "(%s)" % field[2]
                if len(field) == 4: scale += "*("+field[3]+")"
                for p0,s in options.processesToScale:
                    for p in p0.split(","):
                        if re.match(p+"$", pname): scale += "*("+s+")"
                to_norm = True
            elif len(field) == 2:
                pass
            elif len(field) == 3:
                tty.setScaleFactor(field[2])
            else:
                logging.warning("Poorly formatted line: ", field)
                raise RuntimeError                    
            # Adjust free-float and fixed from command line
            for p0 in options.processesToFloat:
                for p in p0.split(","):
                    if re.match(p+"$", pname): tty.setOption('FreeFloat', True)
            for p0 in options.processesToFix:
                for p in p0.split(","):
                    if re.match(p+"$", pname): tty.setOption('FreeFloat', False)
            for p0, p1 in options.processesToPeg:
                for p in p0.split(","):
                    if re.match(p+"$", pname): tty.setOption('PegNormToProcess', p1)
            for p0, p1 in options.processesToSetNormSystematic:
                for p in p0.split(","):
                    if re.match(p+"$", pname): tty.setOption('NormSystematic', float(p1))
            if pname not in self._rank: self._rank[pname] = len(self._rank)

            if to_norm: 
                if options.sumGenWeightFromHisto:
                    ">>> Total sumgenweights from histograms = %s (process = %s)" % (str(total_w),pname)
                for tty in ttys: 
                    if options.weight: 
                        tty.setScaleFactor("%s*%g" % (scale, 1000.0/total_w))
                    else: 
                        tty.setScaleFactor(1)
    def listProcesses(self,allProcs=False):
        ret = self.listSignals(allProcs=allProcs) + self.listBackgrounds(allProcs=allProcs)
        if 'data' in self._allData.keys(): ret.append('data')
        return ret
    def listOptionsOnlyProcesses(self):
        return self._optionsOnlyProcesses.keys()
    def isBackground(self,process):
        return process != 'data' and not self._isSignal[process]
    def isSignal(self,process):
        return self._isSignal[process]
    def listSignals(self,allProcs=False):
        ret = [ p for p in self._allData.keys() if p != 'data' and self._isSignal[p] and (self.getProcessOption(p, 'SkipMe') != True or allProcs) ]
        ret.sort(key = lambda n : self._rank[n])
        return ret
    def listBackgrounds(self,allProcs=False):
        ret = [ p for p in self._allData.keys() if p != 'data' and not self._isSignal[p] and (self.getProcessOption(p, 'SkipMe') != True or allProcs) ]
        ret.sort(key = lambda n : self._rank[n])
        return ret
    def hasProcess(self,process):
        return process in self._allData
    def scaleProcess(self,process,scaleFactor):
        for tty in self._allData[process]: tty.setScaleFactor(scaleFactor)
    def scaleUpProcess(self,process,scaleFactor):
        for tty in self._allData[process]: 
            tty.setScaleFactor( "((%s) * (%s))" % (tty.getScaleFactor(),scaleFactor) )
    def getProcessOption(self,process,name,default=None,noThrow=False):
        if process in self._allData:
            return self._allData[process][0].getOption(name,default=default)
        elif process in self._optionsOnlyProcesses:
            options = self._optionsOnlyProcesses[process]
            return options[name] if name in options else default
        elif noThrow:
            return default
        else: raise RuntimeError("Can't get option %s for undefined process %s" % (name,process))
    def setProcessOption(self,process,name,value):
        if process in self._allData:
            return self._allData[process][0].setOption(name,value)
        elif process in self._optionsOnlyProcesses:
            self._optionsOnlyProcesses[process][name] = value
        else: raise RuntimeError("Can't set option %s for undefined process %s" % (name,process))
    def getScales(self,process):
        return [ tty.getScaleFactor() for tty in self._allData[process] ] 
    def setScales(self,process,scales):
        for (tty,factor) in zip(self._allData[process],scales): tty.setScaleFactor(factor,mcCorrs=False)
    def getYields(self,cuts,process=None,nodata=False,makeSummary=False,noEntryLine=False):
        ## first figure out what we want to do
        tasks = []
        for key,ttys in self._allData.items():
            if key == 'data' and nodata: continue
            if process != None and key != process: continue
            for tty in ttys:
                tasks.append((key,tty,cuts,noEntryLine,None))
        ## then do the work
        if self._options.splitFactor > 1 or  self._options.splitFactor == -1:
            tasks = self._splitTasks(tasks)
        retlist = self._processTasks(_runYields, tasks,name="yields")
        ## then gather results with the same process
        mergemap = {}
        for (k,v) in retlist: 
            if k not in mergemap: mergemap[k] = []
            mergemap[k].append(v)
        ## and finally merge them
        ret = dict([ (k,mergeReports(v)) for k,v in mergemap.items() ])

        rescales = []
        self.compilePlotScaleMap(self._options.plotscalemap,rescales)
        for p,v in ret.items():
            for regexp in rescales:
                if re.match(regexp[0],p): ret[p]=[v[0], [x*regexp[1] for x in v[1]]]

        regroups = [] # [(compiled regexp,target)]
        self.compilePlotMergeMap(self._options.plotmergemap,regroups)
        for regexp in regroups: ret = self.regroupReports(ret,regexp)

        # if necessary project to different lumi, energy,
        if self._projection:
            self._projection.scaleReport(ret)
        # and comute totals
        if makeSummary:
            allSig = []; allBg = []
            for (key,val) in ret.items():
                if key != 'data':
                    if self._isSignal[key]: allSig.append(ret[key])
                    else: allBg.append(ret[key])
            if self._signals and 'signal' not in ret and len(allSig) > 0:
                ret['signal'] = mergeReports(allSig)
            if self._backgrounds and 'background' not in ret and len(allBg) > 0:
                ret['background'] = mergeReports(allBg)
        return ret
    def getPlotsRaw(self,name,expr,bins,cut,process=None,nodata=False,makeSummary=False,closeTreeAfter=False):
        return self.getPlots(PlotSpec(name,expr,bins,{}),cut,process,nodata,makeSummary,closeTreeAfter)
    def getPlots(self,plotspecs,cut,process=None,nodata=False,makeSummary=False,closeTreeAfter=False):
        ret = { }
        allSig = []; allBg = []
        retlist = []

        print("Collecting plots for all processes")
        for key,ttys in self._allData.items():
            if key == 'data' and nodata: continue
            if process != None and key != process: continue
            for tty in ttys:
                logging.info(f"Processing {key} ({tty.cname()})...")
                if self._options.useRunGraphs:
                    self.allHistos[tty.cname()] = tty.getManyPlotsRaw(cut, plotspecs, None, False)
                    # retlist will be created afterwards, when also manipulating histograms
                else:
                    retlist.append( (key, tty.getManyPlots(plotspecs, cut, None, False)) )   
                if self._options.printYieldsRDF:
                    self.allCutReports[tty.cname()] = tty.getCutReport()
                    
        if self._options.useRunGraphs:        
            # assign all histograms to the graph, to allow for simultaneous loop
            logging.info("Have these processes and components")
            logging.info(list(self.allRDF.keys()))       
            allHistogramsList = []
            for key in list(self.allRDF.keys()):
                allHistogramsList.extend(hist for hist in self.allHistos[key])
            #logging.info("Have %d histograms in total" % len(allHistogramsList))
            logging.info("Creating ROOT.RDF.RunGraphs and triggering loop")
            ROOT.RDF.RunGraphs(allHistogramsList)
            # now other operations that were temporarily skipped to avoid triggering the loop
            logging.info("Done with loop, now some refinements on histograms")
            # ugly to repeat it, but for now it works
            for key,ttys in self._allData.items():
                if key == 'data' and nodata: continue
                if process != None and key != process: continue
                for tty in ttys:
                    retlist.append( (key, tty.refineManyPlots(self.allHistos[tty.cname()], plotspecs)) )
            logging.info("Done :)")        

        # manage yields (a repetition of the part below for plots, but for now let's keep them separate for debugging
        if self._options.printYieldsRDF:

            yieldsPerProcess = {}
            for key,ttys in self._allData.items():
                if key == 'data' and nodata: continue
                if process != None and key != process: continue
                if key not in yieldsPerProcess:
                    # to reuse mergeReports() function in tree2yield.py, the format is an array of array of ntuples
                    # ntuple is (cutname, yields), where yields is an array [weighted yield, error, unweighted yields]
                    # for data weighted and unweighted are just the unweighted, and error is sqrt(yield)
                    # the array of ntuples contain the cutflow for a given part of a process (e.g preVFP of Zmumu)
                    # the array of array has an element (being the array of ntuples) for each part
                    # the third element should not be necessary, but should be kept for now to use the old workflow 
                    yieldsPerProcess[key] = []
                for tty in ttys:
                    tmpvec = []
                    for cutInfo in self.allCutReports[tty.cname()]: # will it work?
                        yields = cutInfo.GetPass()
                        tmpvec.append( (cutInfo.GetName(), [yields, math.sqrt(yields), yields]) )
                        yieldsPerProcess[key].append(tmpvec)

            ## and finally merge them
            retyields = dict([ (k,mergeReports(v)) for k,v in yieldsPerProcess.items() ])
            rescales = []
            self.compilePlotScaleMap(self._options.plotscalemap,rescales)
            for p,v in retyields.items():
                for regexp in rescales:
                    if re.match(regexp[0],p): retyields[p]=[v[0], [x*regexp[1] for x in v[1]]]
            regroups = [] # [(compiled regexp,target)]
            self.compilePlotMergeMap(self._options.plotmergemap,regroups)
            for regexp in regroups:
                retyields = self.regroupReports(retyields,regexp)
            # if necessary project to different lumi, energy,
            if self._projection:
                self._projection.scaleReport(retyields)
            # and compute totals
            # this is actually already done in prettyPrint, so can skip here (or one has to change pretty print to skip them if already present)
            if makeSummary:
                allSig = []; allBg = []
                for (key,val) in retyields.items():
                    if key != 'data':
                        if self._isSignal[key]:
                            logging.debug(f">>> {key} is signal")
                            allSig.append(retyields[key])
                        else:
                            logging.debug(f">>> {key} is background")
                            allBg.append(retyields[key])
                if self._signals and 'signal' not in retyields and len(allSig) > 0:
                    retyields['signal'] = mergeReports(allSig)
                if self._backgrounds and 'background' not in retyields and len(allBg) > 0:
                    retyields['background'] = mergeReports(allBg)
            #return retyields
            self.prettyPrint(retyields)
        # done with yields
            
        rets = []
        for ip,plotspec in enumerate(plotspecs):
            mergemap = {}
            #logging.debug(">>>>>> plotspec %s" % plotspec.name)
            for (k,v) in retlist:
                #logging.debug("------- process %s" % k)
                if not k in mergemap: mergemap[k] = []
                if plotspec.getOption('ProcessRegexp', None):
                    regexp = plotspec.getOption('ProcessRegexp', ".*")
                    if re.match(regexp, k):
                        mergemap[k].append(v[ip])
                    else:
                        mergemap[k].append(None)
                        v.insert(ip,None)
                else:
                    mergemap[k].append(v[ip])
                #logging.debug(" ".join(x.GetName() if x != None else "None" for x in v))

            ret = dict([ (k,mergePlots(plotspec.name+"_"+k,v)) for k,v in mergemap.items() if all(x != None for x in v)])

            rescales = []
            self.compilePlotScaleMap(self._options.plotscalemap,rescales)
            for p,v in ret.items():
                for regexp in rescales:
                    if re.match(regexp[0],p): v.Scale(regexp[1])

            regroups = [] # [(compiled regexp,target)]
            self.compilePlotMergeMap(self._options.plotmergemap,regroups)
            for regexp in regroups: ret = self.regroupPlots(ret,regexp,plotspec)

            # if necessary project to different lumi, energy,
            if self._projection:
                self._projection.scalePlots(ret)
            if makeSummary:
                allSig = [v for k,v in ret.items() if k != 'data'and self._isSignal[k] == True  ]
                allBg  = [v for k,v in ret.items() if k != 'data'and self._isSignal[k] == False ]
                if self._signals and not 'signal' not in ret and len(allSig) > 0:
                    ret['signal'] = mergePlots(plotspec.name+"_signal", allSig)
                    ret['signal'].summary = True
                if self._backgrounds and 'background' not in ret and len(allBg) > 0:
                    ret['background'] = mergePlots(plotspec.name+"_background",allBg)
                    ret['background'].summary = True

            rets.append(ret)
        return rets
    ## CANDIDATE FOR REMOVAL
    def prepareForSplit(self):
        ttymap = {}
        for key,ttys in self._allData.items():
            for tty in ttys:
                if not tty.hasEntries(): 
                    #print "For tty %s/%s, I don't have the number of entries" % (tty._name, tty._cname)
                    ttymap[id(tty)] = tty
        if len(ttymap):
            retlist = self._processTasks(_runGetEntries, ttymap.items(), name="GetEntries")
            for ttid, entries in retlist:
                ttymap[ttid].setEntries(entries)
    def prettyPrint(self,reports,makeSummary=True):
        allSig = []; allBg = []
        for key in reports:
            #print(f"XXXXXX  {key}")
            if all(key != x for x in ['data', 'background', 'signal']):
                if self._isSignal[key]:
                    allSig.append((key,reports[key]))
                else:
                    allBg.append((key,reports[key]))
        allSig.sort(key = lambda x: self._rank[x[0]])
        allBg.sort( key = lambda x: self._rank[x[0]])
        table = allSig + allBg
        if makeSummary:
            if len(allSig)>1:
                if "signal" in reports:
                    table.append(('ALL SIG',reports['signal']))
                else:
                    table.append(('ALL SIG',mergeReports([v for n,v in allSig])))
            if len(allBg)>1:
                if "background" in reports:
                    table.append(('ALL BKG',reports['background']))
                else:
                    table.append(('ALL BKG',mergeReports([v for n,v in allBg])))
        if "data" in reports: table += [ ('DATA', reports['data']) ]
        for fomname in self._options.figureOfMerit:
            fom = FOM_BY_NAME[fomname]
            nrows = len(table[0][1])
            table += [ (fomname, [ (None, [fom(self, reports, row), 0, 0]) for row in xrange(nrows) ] ) ]

        # maximum length of the cut descriptions
        clen = max([len(cut) for cut,yields in table[0][1]]) + 3
        cfmt = "%%-%ds" % clen;

        fmtlen = 10
        nfmtL = "  %8d"
        nfmtS = "  %8.2f" if self._options.weight else nfmtL
        nfmtX = "  %8.4f" if self._options.weight else nfmtL

        if self._options.errors:
            nfmtS+=u" %7.2f"
            nfmtX+=u" %7.4f"
            nfmtL+=u" %7.2f"
            fmtlen+=9
        if self._options.fractions:
            nfmtS+=" %7.1f%%"
            nfmtX+=" %7.1f%%"
            nfmtL+=" %7.1f%%"
            fmtlen+=8

        stdoutBackup = sys.stdout
        sys.stdout = open(self._options.yieldsOutfile, "w")
        
        if self._options.txtfmt == "txt":
            print("CUT".center(clen), end='')
            for h,r in table: 
                if len("   "+h) <= fmtlen:
                    print(("   "+h).center(fmtlen), end='')
                elif len(h) <= fmtlen:
                    print(h.center(fmtlen), end='')
                else:
                    print(h[:fmtlen], end='')
            print("")
            print("-"*((fmtlen+1)*len(table)+clen))
            for i,(cut,dummy) in enumerate(table[0][1]):
                print(cfmt % cut, end='')
                for name,report in table:
                    (nev,err,nev_run_upon) = report[i][1]
                    den = report[i-1][1][0] if i>0 else 0
                    fraction = nev/float(den) if den > 0 else 1
                    if self._options.nMinusOne: 
                        fraction = report[-1][1][0]/float(nev) if nev > 0 else 1
                    elif self._options.nMinusOneInverted: 
                        fraction = float(nev)/report[-1][1][0] if report[-1][1][0] > 0 else 1
                    toPrint = (nev,)
                    if self._options.errors:    toPrint+=(err,)
                    if self._options.fractions: toPrint+=(fraction*100,)
                    if self._options.weight and nev < 1000: print(( nfmtS if nev > 0.2 else nfmtX) % toPrint, end='')
                    else                                  : print(nfmtL % toPrint, end='')
                print("")
        elif self._options.txtfmt in ("tsv","csv","dsv","ssv"):
            sep = { 'tsv':"\t", 'csv':",", 'dsv':';', 'ssv':' ' }[self._options.txtfmt]
            if len(table[0][1]) == 1:
                for k,r in table:
                    if sep in k:
                        if self._options.txtfmt in ("tsv","ssv"):
                            k = k.replace(sep,"_")
                        else:
                            k = '"'+k.replace('"','""')+'"'
                    (nev,err,fraction) = r[0][1][0], r[0][1][1], 1.0
                    toPrint = (nev,)
                    if self._options.errors:    toPrint+=(err,)
                    if self._options.fractions: toPrint+=(fraction*100,)
                    if self._options.weight and nev < 1000: ytxt = ( nfmtS if nev > 0.2 else nfmtX) % toPrint
                    else                                  : ytxt = nfmtL % toPrint
                    print("%s%s%s" % (k,sep,sep.join(ytxt.split())))
                print("")

        sys.stdout.close()
        sys.stdout = stdoutBackup
        logging.info(f"Yields saved in file {self._options.yieldsOutfile}")
                
    # apparently used nowhere
    #def _getYields(self,ttylist,cuts):
    #    return mergeReports([tty.getYields(cuts) for tty in ttylist])
    def __str__(self):
        mystr = ""
        for a in self._allData:
            mystr += str(a) + '\n' 
        for a in self._data:
            mystr += str(a) + '\n' 
        for a in self._signals:
            mystr += str(a) + '\n' 
        for a in self._backgrounds:
            mystr += str(a) + '\n'
        return mystr[:-1]
    def processEvents(self,eventLoop,cut):
        for p in self.listProcesses():
            for tty in self._allData[p]:
                tty.processEvents(eventLoop,cut)
    def compilePlotMergeMap(self,inlist,relist):
        for m in inlist:
            to,fro = m.split("=")
            if to[-1] == "+": to = to[:-1]
            else: raise RuntimeError('Incorrect plotmergemap format: %s'%m)
            to = to.strip()
            for k in fro.split(","):
                relist.append((re.compile(k.strip()+"$"), to))
    def compilePlotScaleMap(self,inlist,relist):
        for m in inlist:
            dset,scale = m.split("=")
            if dset[-1] == "*": dset = dset[:-1]
            else: raise RuntimeError('Incorrect plotscalemap format: %s'%m)
            relist.append((re.compile(dset.strip()+"$"),float(scale)))
    def regroupReports(self,pmap,regexp):
        patt, to = regexp
        mergemap={}
        for (k,v) in pmap.items():
            k2 = k
            if k2 != to and re.match(patt,k2): k2 = to
            if k2 not in mergemap: mergemap[k2]=[]
            mergemap[k2].append(v)
        return dict([ (k,mergeReports(v)) for k,v in mergemap.items() ])
    def regroupPlots(self,pmap,regexp,pspec):
        patt, to = regexp
        mergemap={}
        for (k,v) in pmap.items():
            k2 = k
            if k2 != to and re.match(patt,k2): k2 = to
            if k2 not in mergemap: mergemap[k2]=[]
            mergemap[k2].append(v)
        return dict([ (k,mergePlots(pspec.name+"_"+k,v)) for k,v in mergemap.items() ])
    def stylePlot(self,process,plot,pspec,mayBeMissing=False):
        if process in self._allData:
            for tty in self._allData[process]: 
                tty._stylePlot(plot,pspec)
                break
        elif process in self._optionsOnlyProcesses:
            opts = self._optionsOnlyProcesses[process]
            stylePlot(plot, pspec, lambda x: opts[x[0]] if key in opts else default)
        elif not mayBeMissing:
            raise KeyError("Process %r not found" % process)
    def _processTasks(self,func,tasks,name=None):
        #timer = ROOT.TStopwatch()
        #print "Starting job %s with %d tasks, %d threads" % (name,len(tasks),self._options.jobs)
        logging.info('self._options.jobs', self._options.jobs)
        if self._options.jobs == 0: 
            retlist = map(func, tasks)
        else:
            from multiprocessing import Pool
            pool = Pool(self._options.jobs)
            retlist = pool.map(func, tasks, 1)
            pool.close()
            pool.join()
        #print "Done %s in %s s at %.2f " % (name,timer.RealTime(),0.001*(long(ROOT.gSystem.Now()) - _T0))
        return retlist
    def _splitTasks(self,tasks):
        nsplit = self._options.splitFactor
        if nsplit == -1: nsplit = self._options.jobs
        if nsplit <= 1: return tasks
        newtasks = []
        if not self._options.splitDynamic:
            for task in tasks:
                for fsplit in [ (i,nsplit) for i in xrange(nsplit) ]:
                    newtasks.append( tuple( (list(task)[:-1]) + [fsplit] ) )
        else:
            self.prepareForSplit() 
            #print "Original task list has %d entries; split factor %d." % (len(tasks), nsplit)
            maxent = max( task[1].getEntries() for task in tasks )
            grain  = maxent / nsplit / 1 # factor 2 may be optimized
            #print "Largest task has %d entries. Will use %d as grain " % (maxent, grain)
            if grain < 10: return tasks # sanity check
            newtasks_wsize = []
            for task in tasks:
                tty = task[1]; 
                entries = tty.getEntries()
                chunks  = min(max(1, int(round(entries/grain))), nsplit)
                fsplits = [ (i,chunks) for i in xrange(chunks) ]
                #print "    task %s/%s has %d entries. N/g = %.1f, chunks = %d" % (tty._name, tty._cname, entries, entries/float(grain), chunks)
                for fsplit in fsplits:
                    newtasks_wsize.append( (entries/float(chunks), tuple( (list(task)[:-1]) + [fsplit] ) ) )
            if self._options.splitSort:
                newtasks_wsize.sort(key = lambda x : x[0], reverse = True)
                #for s,t in newtasks_wsize: print "\t%9d %s/%s %s" % (s,t[1]._name, t[1]._cname, t[-1])
            newtasks = [ task for (size,task) in newtasks_wsize ]
        #print "New task list has %d entries; actual split factor %.2f" % (len(newtasks), len(newtasks)/float(len(tasks)))
        return newtasks


def addMCAnalysisOptions(parser,addTreeToYieldOnesToo=True):
    if addTreeToYieldOnesToo: addTreeToYieldOptions(parser)
    parser.add_argument("-j", "--jobs", type=int, default=0, help="Use N threads");
    parser.add_argument("--split-factor", dest="splitFactor", type=int, default=0, help="Use N chunks per sample (-1 means to use the same as what passed to -j, which appears to work well in the average case)");
    parser.add_argument("--split-static", dest="splitDynamic", action="store_false", help="Make the splitting dynamic (reduce the chunks for small samples)");
    parser.add_argument("--split-nosort", dest="splitSort", action="store_false", help="Make the splitting dynamic (reduce the chunks for small samples)");
    parser.add_argument("-P", "--path", action="append", type=str, default=[], help="Path to directory with input trees and pickle files. Can supply multiple paths which will be searched in order. (default: ./") 
    parser.add_argument("--RP", "--remote-path", dest="remotePath", type=str, help="path to remote directory with trees, but not other metadata (default: same as path)") 
    parser.add_argument("-p", "--process", dest="processes", type=str, default=[], action="append", help="Processes to print (comma-separated list of regexp, can specify multiple ones)");
    parser.add_argument("--pg", "--pgroup", dest="premap", type=str, default=[], action="append", help="Group proceses into one. Syntax is '<newname> := (comma-separated list of regexp)', can specify multiple times. Note that it is applied _before_ -p, --sp and --xp");
    parser.add_argument("--xf", "--exclude-files", dest="filesToExclude", type=str, default=[], action="append", help="Files to exclude (comma-separated list of regexp, can specify multiple ones)");
    parser.add_argument("--xp", "--exclude-process", dest="processesToExclude", type=str, default=[], action="append", help="Processes to exclude (comma-separated list of regexp, can specify multiple ones)");
    parser.add_argument("--sf", "--swap-files", dest="filesToSwap", type=str, default=[], nargs=2, action="append", help="--swap-files X Y uses file Y instead of X in the MCA");
    parser.add_argument("--sp", "--signal-process", dest="processesAsSignal", type=str, default=[], action="append", help="Processes to set as signal (overriding the '+' in the text file)");
    parser.add_argument("--float-process", "--flp", dest="processesToFloat", type=str, default=[], action="append", help="Processes to set as freely floating (overriding the 'FreeFloat' in the text file; affects e.g. mcPlots with --fitData)");
    parser.add_argument("--fix-process", "--fxp", dest="processesToFix", type=str, default=[], action="append", help="Processes to set as not freely floating (overriding the 'FreeFloat' in the text file; affects e.g. mcPlots with --fitData)");
    parser.add_argument("--peg-process", dest="processesToPeg", type=str, default=[], nargs=2, action="append", help="--peg-process X Y make X scale as Y (equivalent to set PegNormToProcess=Y in the mca.txt)");
    parser.add_argument("--scale-process", dest="processesToScale", type=str, default=[], nargs=2, action="append", help="--scale-process X Y make X scale by Y (equivalent to add it in the mca.txt)");
    parser.add_argument("--process-norm-syst", dest="processesToSetNormSystematic", type=str, default=[], nargs=2, action="append", help="--process-norm-syst X Y sets the NormSystematic of X to be Y (for plots, etc. Overrides mca.txt)");
    parser.add_argument("--AP", "--all-processes", dest="allProcesses", action="store_true", help="Include also processes that are marked with SkipMe=True in the MCA.txt")
    parser.add_argument("--use-cnames",  dest="useCnames", action="store_true", help="Use component names instead of process names (for debugging)")
    parser.add_argument("--project", dest="project", type=str, help="Project to a scenario (e.g 14TeV_300fb_scenario2)")
    parser.add_argument("--plotgroup", dest="plotmergemap", type=str, default=[], action="append", help="Group plots into one. Syntax is '<newname> := (comma-separated list of regexp)', can specify multiple times. Note it is applied after plotting.")
    parser.add_argument("--scaleplot", dest="plotscalemap", type=str, default=[], action="append", help="Scale plots by this factor (before grouping). Syntax is '<newname> := (comma-separated list of regexp)', can specify multiple times.")
    parser.add_argument("-t", "--tree", default='treeProducerWMass', help="Pattern for tree name");
    parser.add_argument("--fom", "--figure-of-merit", dest="figureOfMerit", type=str, default=[], action="append", help="Add this figure of merit to the output table (S/B, S/sqrB, S/sqrSB)")
    parser.add_argument("--max-genWeight-procs", dest="maxGenWeightProc", type=str, nargs=2, action="append", default=[], help="maximum genWeight to be used for a given MC process (first value is a regular expression for the process, second is the max weight). This option effectively applies a cut on genWeight, and also modifies the sum of genweights. Can be specified more than once for different processes. This will cut away events with larger weights");
    parser.add_argument("--clip-genWeight-toMax", dest="clipGenWeightToMax", action="store_true", help="It only works with --nanoaod-tree when using --max-genWeight-procs, setting large weights to the max instead of rejecting the event");
    parser.add_argument("--no-heppy-tree", dest="noHeppyTree", action="store_true", help="Set to true to read root files when they were not made with Heppy (different convention for path names, might need to be adapted)");
    parser.add_argument("--nanoaod-tree", dest="nanoaodTree", action="store_true", help="Set to true to read root files from nanoAOD");
    parser.add_argument("--filter-proc-files", dest="filterProcessFiles", type=str, nargs=2, action="append", default=[], help="Can use this option to override second field on each process line in MCA file, so to select few files without modifying the MCA file (e.g. for tests). E.g. --filter-proc-files 'W.*' '.*_12_.*' to only use files with _12_ in their name. Only works with option --nanoaod-tree");
    parser.add_argument("--sum-genWeight-fromHisto", dest="sumGenWeightFromHisto", action="store_true", help="If True, compute sum of gen weights from histogram (when using --nanoaod-tree)");
    parser.add_argument("--no-rdf-runGraphs", dest="useRunGraphs", action="store_false", help="If True, use RDF::RunGraphs to make all histograms for all processes at once");
    parser.add_argument("-v", "--verbose", type=int, default=3, choices=[0,1,2,3,4], help="Set verbosity level with logging, the larger the more verbose");
    parser.add_argument("--json", type=str, default="pileupStuff/Cert_271036-284044_13TeV_Legacy2016_Collisions16_JSON.txt", help="json file to filter data. Default is the one for UL2016") 
    parser.add_argument("--scale-factor-file", dest="scaleFactorFile", type=str, default="./testMuonSF/scaleFactorProduct_31Mar2021.root", help="File to be used for scale factors")
    parser.add_argument("sampleFile", type=str, help="Text file with sample definitions");
    parser.add_argument("cutFile", type=str, help="Text file with cut definitions");

    
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    addMCAnalysisOptions(parser)
    args = parser.parse_args()
    options = args
    setLogging(args.verbose)
    if not options.path and not options.nanoaodTree: options.path = ['./']
    tty = TreeToYield(args.sampleFile,options) if ".root" in args.sampleFile else MCAnalysis(args.sampleFile,options)
    cf  = CutsFile(args.cutFile,options)
    # following part assumes that one can pass more cut files, but we don't use it for now
    for cutFile in args[2:]:
        temp = CutsFile(cutFile,options)
        for cut in temp.cuts():
            cf.add(cut[0],cut[1])
    report = tty.getYields(cf)#, process=options.process)
    tty.prettyPrint(report)
