#!/usr/bin/env python
#from tree2yield import *
from CMGTools.WMass.plotter.tree2yield import *
from CMGTools.WMass.plotter.projections import *
from CMGTools.WMass.plotter.figuresOfMerit import FOM_BY_NAME
import pickle, re, random, time, glob, math

#_T0 = long(ROOT.gSystem.Now())

## These must be defined as standalone functions, to allow runing them in parallel
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

def _runApplyCut(args):
    key,tty,cut,fsplit = args
    return (key, tty.cutToElist(cut,fsplit=fsplit))

def _runGetEntries(args):
    key,tty = args
    return (key, tty.getEntries())


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
        self.readMca(samples,options)


    def getSumGenWeightMCfromHisto(self, pname, rootfile, verbose=False):

        maxGenWgt = None
        sumGenWeights = 1.0
        nUnweightedEvents = 1.0

        ROOT.gEnv.SetValue("TFile.AsyncReading", 1);
        #tmp_rootfile = ROOT.TXNetFile(rootfile+"?readaheadsz=65535")
        tmp_rootfile = ROOT.TFile(rootfile+"?readaheadsz=65535")
 
        if self._options.weight and len(self._options.maxGenWeightProc):
            # get sum of weights from histograms, filtering some events with large weights
            for procRegExp,tmp_maxGenWgt in self._options.maxGenWeightProc:
                if re.match(procRegExp,pname):
                    #tmp_maxGenWgt is a string, convert to float
                    maxGenWgt = fabs(float(tmp_maxGenWgt))
                    log10_maxGenWgt = math.log10(maxGenWgt)
                    histo_sumgenweight = tmp_rootfile.Get('hGenWeights')
                    if not histo_sumgenweight:
                        raise RuntimeError, "Can't get histogram hGenWeights from %s.\nMake sure it is available when using option --max-genWeight-procs" % rootfile
                    minBin = histo_sumgenweight.GetXaxis().FindFixBin(-log10_maxGenWgt)
                    maxBin = histo_sumgenweight.GetXaxis().FindFixBin(log10_maxGenWgt)
                    sumGenWeights = histo_sumgenweight.Integral(minBin,maxBin)                  
                    if self._options.clipGenWeightToMax:
                        # then the other histogram, whose integral is the number of events before any cut
                        histo_numweight = tmp_rootfile.Get('hNumWeights')
                        if not histo_numweight:
                            raise RuntimeError, "Can't get histogram hNumWeights from %s.\nMake sure it is available when using option --max-genWeight-procs and --clip-genWeight-toMax" % rootfile      
                        # get actual upper threshold based on the bin edge
                        maxGenWgt = math.pow(10.0,histo_sumgenweight.GetXaxis().GetBinUpEdge(maxBin))
                        if verbose:
                            print "INFO >>> Process %s -> actual genWeight threshold set to %s" % (pname,str(maxGenWgt))
                        nLargeWeight = histo_numweight.Integral(0,max(0,minBin-1)) + histo_numweight.Integral(min(histo_numweight.GetNbinsX()+1,maxBin+1), histo_numweight.GetNbinsX()+1)
                        sumGenWeights += (nLargeWeight*maxGenWgt) 
                        nUnweightedEvents = histo_numweight.Integral(0,histo_numweight.GetNbinsX()+1)
        else:
            # get sum of weights from Runs tree (faster, but only works if not cutting away large genWeight)
            tmp_tree = tmp_rootfile.Get('Runs')
            if not tmp_tree or tmp_tree == None:
                raise RuntimeError, "Can't get tree Runs from %s.\n" % rootfile
            tmp_hist = ROOT.TH1D("sumweights","",1,0,10)
            tmp_hist2 = ROOT.TH1D("sumcount","",1,0,10)
            tmp_tree.Draw("1>>sumweights", "genEventSumw")
            tmp_hist = ROOT.gROOT.FindObject("sumweights")
            sumGenWeights = tmp_hist.Integral()
            tmp_tree.Draw("1>>sumcount", "genEventCount")
            tmp_hist2 = ROOT.gROOT.FindObject("sumcount")
            nUnweightedEvents = tmp_hist2.Integral()

        tmp_rootfile.Close()                            
        return (maxGenWgt, sumGenWeights, nUnweightedEvents)


    def getSumGenWeightMC(self, pname, rootfile):

        isTChain = False
        chain = None
        if not isinstance(rootfile,str):
            # then rootfile is a list of files
            print "Building chain to compute sum of gen weights ..."
            isTChain = True
            chain = ROOT.TChain("Events","chainToComputeSumWeight_Events")
            chainRuns = ROOT.TChain("Runs","chainToComputeSumWeight_Runs")
            for f in rootfile:
                chain.Add(f)
                chainRuns.Add(f)
            print "Done building chain!"


        maxGenWgt = None
        sumGenWeights = 1.0
        nUnweightedEvents = 1.0
        if self._options.weight and len(self._options.maxGenWeightProc):
            # get sum of weights from Events tree, filtering some events with large weights
            # this assumes the trees are unskimmed to correctly compute the sum!!
            for procRegExp,tmp_maxGenWgt in self._options.maxGenWeightProc:
                if re.match(procRegExp,pname):
                    #tmp_maxGenWgt is a string, convert to float
                    maxGenWgt = fabs(float(tmp_maxGenWgt))
                    if isTChain:
                        tmp_tree = chain
                        if not tmp_tree or tmp_tree == None:
                            raise RuntimeError, "Can't get chain with tree Events in function getSumGenWeightMC().\n"
                        nEvents = tmp_tree.GetEntries()
                        tmp_tree_runs = chainRuns
                        if not tmp_tree_runs or tmp_tree_runs == None:
                            raise RuntimeError, "Can't get chain with tree Runs in function getSumGenWeightMC().\n"
                        nGenEvents = 0
                        for i,event in enumerate(tmp_tree_runs): # this is few events, equal to number of files
                            nGenEvents += event.genEventCount
        
                    else:
                        ROOT.gEnv.SetValue("TFile.AsyncReading", 1);
                        #tmp_rootfile = ROOT.TXNetFile(rootfile+"?readaheadsz=65535")
                        tmp_rootfile = ROOT.TFile(rootfile+"?readaheadsz=65535")
                        tmp_tree = tmp_rootfile.Get('Events')
                        if not tmp_tree or tmp_tree == None:
                            raise RuntimeError, "Can't get tree Events from %s.\n" % rootfile
                    # check the tree was not skimmed                                    
                        nEvents = tmp_tree.GetEntries()
                        tmp_tree_runs = tmp_rootfile.Get('Runs')
                        if not tmp_tree_runs or tmp_tree_runs == None:
                            raise RuntimeError, "Can't get tree Runs from %s.\n" % rootfile
                        for i,event in enumerate(tmp_tree_runs): # should be one event only
                            nGenEvents = event.genEventCount
                            if i: break

                    if nEvents != nGenEvents:
                        print "nEvents = %d    nGenEvents = %f" % (nEvents,nGenEvents)
                        raise RuntimeError, "You are trying to remove or clip large gen weights, so I am recomputing the sum of gen weights excluding them.\nHowever, it seems the file\n%s\n you are using contains a skimmed tree.\nIn this way the sum of gen weights will be wrong" % rootfile
                        ## check was ok

                    tmp_hist = ROOT.TH1D("sumweights","",1,0,10)
                    if self._options.clipGenWeightToMax:
                        # set weight to max
                        # TMath::Sign(a,b) returns a with the sign of b
                        nUnweightedEvents = tmp_tree.Draw("1>>sumweights", "TMath::Sign(TMath::Min(abs(genWeight),%s),genWeight)" % str(maxGenWgt))           
                    else:
                        # reject event with weight > max
                        nUnweightedEvents = tmp_tree.Draw("1>>sumweights", "genWeight*(abs(genWeight) < %s)" % str(maxGenWgt))           
                    #tmp_hist = ROOT.gROOT.FindObject("sumweights")
                    sumGenWeights = tmp_hist.Integral()
                    if not isTChain:
                        tmp_rootfile.Close()

        else:
            if isTChain:
                tmp_tree = chain
                if not tmp_tree or tmp_tree == None:
                    raise RuntimeError, "Can't get chain with tree Runs\n"
            else:
                # get sum of weights from Runs tree (faster, but only works if not cutting away large genWeight)
                ROOT.gEnv.SetValue("TFile.AsyncReading", 1);
                #tmp_rootfile = ROOT.TXNetFile(rootfile+"?readaheadsz=65535")
                tmp_rootfile = ROOT.TFile(rootfile+"?readaheadsz=65535")
                tmp_tree = tmp_rootfile.Get('Runs')
                if not tmp_tree or tmp_tree == None:
                    raise RuntimeError, "Can't get tree Runs from %s.\n" % rootfile

            tmp_hist = ROOT.TH1D("sumweights","",1,0,10)
            tmp_hist2 = ROOT.TH1D("sumcount","",1,0,10)
            tmp_tree.Draw("1>>sumweights", "genEventSumw")
            tmp_hist = ROOT.gROOT.FindObject("sumweights")
            sumGenWeights = tmp_hist.Integral()
            tmp_tree.Draw("1>>sumcount", "genEventCount")
            tmp_hist2 = ROOT.gROOT.FindObject("sumcount")
            nUnweightedEvents = tmp_hist2.Integral()
            if not isTChain:
                tmp_rootfile.Close()                            

        if isTChain:
            print "Done computing sum of gen weights!"
            #print ">>> pname = %s   sumGenWeights = %f " % (pname,sumGenWeights)
        return (maxGenWgt, sumGenWeights, nUnweightedEvents)


    def readMca(self,samples,options,addExtras={},field0_addExtras=""):

        # when called with not empty addExtras, issue a warning in case you are overwriting settings
        if "ALLOW_OVERWRITE_SETTINGS" in addExtras:
            if addExtras["ALLOW_OVERWRITE_SETTINGS"] == True:
                print "### WARNING: found setting ALLOW_OVERWRITE_SETTINGS=True in %s. " % str(field0_addExtras)
                print "### Will use new settings overwriting those in mca included file %s." % str(samples)                    

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
            for k,v in addExtras.iteritems():
                if k in extra: 
                    if "ALLOW_OVERWRITE_SETTINGS" in addExtras and addExtras["ALLOW_OVERWRITE_SETTINGS"] == True:
                        pass
                    else:
                        raise RuntimeError, 'You are trying to overwrite an extra option (%s - %s ) already set (did you forget ALLOW_OVERWRITE_SETTINGS=True ?)' % (k, v)
                extra[k] = v
            field = [f.strip() for f in line.split(':')]
            if len(field) == 1 and field[0] == "*":
                if len(self._allData): raise RuntimeError, "MCA defaults ('*') can be specified only before all processes"
                #print "Setting the following defaults for all samples: "
                for k,v in extra.iteritems():
                    #print "\t%s: %r" % (k,v)
                    self.init_defaults[k] = v
                continue
            else:
                for k,v in self.init_defaults.iteritems():
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
                if 'IncludeMca' not in extra: raise RuntimeError, 'You have declared a component with IncludeMca format, but not included this option'
                extra_to_pass = copy(extra)
                del extra_to_pass['IncludeMca']
                self.readMca(extra['IncludeMca'],options,addExtras=extra_to_pass,field0_addExtras=field0noPostFix) # call readMca recursively on included mca files
                continue
            # Customize with additional weight if requested
            if 'AddWeight' in extra:
                if len(field)<2: raise RuntimeError, 'You are trying to set an additional weight, but there is no weight initially defined for this component'
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
                            print "Error >>> Process %s '[XXX]' found in TreePath, but XXX not defined in MCA file" % (pname)
                            quit()
                    else:                        
                        pathsToSearch = [extra['TreePath']]
                else:
                    if not options.path:
                        print "Warning: you didn't specify a path to ntuples with option -P, but process %s has no 'TreePath' key in the MCA file. Please specify a valid path." % pname
                        quit()
                print "INFO >>> Process %s -> searching files in %s" % (pname,repr(pathsToSearch))
                for treepath in pathsToSearch:
                    for dirpath, dirnames, filenames in os.walk(treepath):
                        # either use regexp or simply that subpath is in the path (so one doesn't need 
                        # to start subpath with .* for a proper match if a word without wildcards is used)
                        if subPath_regexp.match(dirpath+"/") or subpath in (dirpath+"/"):
                            filtered_fnames = [f for f in filenames if f.endswith(".root") and field1_regexp.match(f)]
                            for filename in filtered_fnames:
                                cnamesWithPath.append(os.path.join(dirpath, filename))
                print "INFO >>> Process %s -> %d files selected" % (pname,len(cnamesWithPath))
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
                            print "INFO >>> Process %s -> clipping genWeight to |x| < %s" % (pname,tmp_maxGenWgt)
                        else:
                            print "INFO >>> Process %s -> rejecting events with genWeight > %s" % (pname,tmp_maxGenWgt)

            ### CHECKPOINT
            tch = None # TChain, yet to be defined, done below
            list_rootfiles = [] # files that will be used in the TChain, if applicable (should be equal to allFiles)
            allFiles = cnamesWithPath if len(cnamesWithPath) else cnames
            nAllFiles = len(allFiles)

            iFile = 0
            for cnameWithPath in allFiles:
                cname = os.path.basename(cnameWithPath)
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
                        print "Error >>> Process %s '[XXX]' found in FriendDir, but XXX not defined in MCA file" % (pname)
                        quit()

                if subpath != "" and not subpath.endswith("/"):
                    subpath += "/"
                #print "subpath = %s" % subpath

                basepath = None
                if options.nanoaodTree and os.path.exists(os.path.dirname(cnameWithPath)):
                    basepath = cnameWithPath
                    if not basepath:
                        raise RuntimeError("%s -- ERROR: %s process not found in paths (%s)" % (__name__, cname, pathsToSearch))
                else:
                    for treepath in options.path:
                        if os.path.exists(treepath+"/"+subpath+cname):
                            basepath = treepath
                            break
                    if not basepath:
                        raise RuntimeError("%s -- ERROR: %s process not found in paths (%s)" % (__name__, cname, repr(options.path)))

                rootfile = "%s/%s/%s/%s_tree.root" % (basepath, cname, treename, treename)
                if options.noHeppyTree:
                    # under development, to run on any kind of root file
                    #print "cname = %s" % cname
                    if "TreeName" in extra:
                        rootfile = "%s/%s.root" % (basepath, treename)                    
                    else:
                        rootfile = "%s/%s" % (basepath, cname)                    
                    #print "rootfile: %s" % rootfile
                    #print "objname : %s" % objname
                elif options.nanoaodTree:
                    # under development, to run on nanoaod (similar to case above, but that 
                    # is for another purpose, so keep it separate for now)
                    rootfile = cnameWithPath
                    #print "cname = %s" % cname
                    #print "rootfile: %s" % rootfile
                    #print "objname : %s" % objname
                    #
                    prepath = ''
                    if not 'root:/' in rootfile:
                        if '/eos/user/' in rootfile: 
                            prepath = 'root://eosuser.cern.ch//'
                        elif '/eos/cms/store/' in rootfile: 
                            prepath = 'root://eoscms.cern.ch//'
                    rootfile = prepath+rootfile
                else:
                    if options.remotePath:
                        rootfile = "root:%s/%s/%s_tree.root" % (options.remotePath, cname, treename)
                    elif os.path.exists(rootfile+".url"): #(not os.path.exists(rootfile)) and :
                        rootfile = open(rootfile+".url","r").readline().strip()
                    elif (not os.path.exists(rootfile)) and os.path.exists("%s/%s/%s/tree.root" % (basepath, cname, treename)):
                        # Heppy calls the tree just 'tree.root'
                        rootfile = "%s/%s/%s/tree.root" % (basepath, cname, treename)
                        prepath = ''
                        if not 'root:/' in rootfile:
                            if   '/eos/user/'      in rootfile: prepath = 'root://eosuser.cern.ch//'
                            elif '/eos/cms/store/' in rootfile: prepath = 'root://eoscms.cern.ch//'
                        rootfile = prepath+rootfile
                    elif (not os.path.exists(rootfile)) and os.path.exists("%s/%s/%s/tree.root.url" % (basepath, cname, treename)):
                        # Heppy calls the tree just 'tree.root'
                        rootfile = "%s/%s/%s/tree.root" % (basepath, cname, treename)
                        rootfile = open(rootfile+".url","r").readline().strip()

                ## needed temporarily
                pckfile = basepath+"/%s/skimAnalyzerCount/SkimReport.pck" % cname

                ## MIDDLE CHECKPOINT
                if options.useTChain:
                    if tch == None:
                        tch = ROOT.TChain(objname,pname)
                    #print rootfile
                    nCurrentFileInChain = len(list_rootfiles)
                    if nCurrentFileInChain%2 == 0:
                        sys.stdout.write('>>> Creating chain {0:.2%}  \r'.format(float(nCurrentFileInChain+1)/nAllFiles))
                        sys.stdout.flush()
                    tch.Add(rootfile)
                    list_rootfiles.append(rootfile)
                    if len(list_rootfiles) == nAllFiles:
                        tty = TreeToYield(tch, options, settings=extra, name=pname, cname=cname, objname=objname, frienddir=friendDir); 
                        ttys.append(tty)
                    else:
                        # basically accumulate all files for the chain, and then go to the rest of the code, 
                        # which would be the same for each file belonging to a given process, so we should not be
                        # missing anythin, hopefully
                        continue
                else:
                    tty = TreeToYield(rootfile, options, settings=extra, name=pname, cname=cname, objname=objname, frienddir=friendDir); 
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
                            if (is_w==0): raise RuntimeError, "Can't put together a weighted and an unweighted component (%s)" % cnames
                            is_w = 1
                            total_w = 1.0  # not really used for general trees at the moment
                            scale = "(%s)" % field[2]
                        else:
                            scale = "1"
                            total_w = 1 # not really used for general trees at the moment
                            if (is_w==1): raise RuntimeError, "Can't put together a weighted and an unweighted component (%s)" % cnames
                            is_w = 0
                    elif options.nanoaodTree:

                        #(maxGenWgt, sumGenWeights, nUnweightedEvents) = self.getSumGenWeightMC(pname, list_rootfiles if options.useTChain else rootfile)
                        (maxGenWgt, sumGenWeights, nUnweightedEvents) = self.getSumGenWeightMCfromHisto(pname, rootfile, verbose=(iFile==0))
                        sys.stdout.write('INFO >>> preparing files {0:.2%}  \r'.format(float(iFile+1)/nAllFiles))
                        sys.stdout.flush()



                        if options.weight and True: # True for now, later on this could explicitly require using the actual genWeights as opposed to using sum of unweighted events for MC (see the case for cmgtools below)
                            if (is_w==0): raise RuntimeError, "Can't put together a weighted and an unweighted component (%s)" % cnames
                            is_w = 1; 
                            total_w += sumGenWeights
                            if maxGenWgt != None:
                                if options.clipGenWeightToMax:
                                    scale = "(%s)*TMath::Sign(TMath::Min(abs(genWeight),%s),genWeight)" % (field[2],str(maxGenWgt))
                                else:
                                    scale = "genWeight*(%s)*(abs(genWeight)<%s)" % (field[2],str(maxGenWgt))
                            else:
                                scale = "genWeight*(%s)" % field[2]
                        elif not options.weight:
                            scale = "1"
                            total_w = 1
                            if (is_w==1): raise RuntimeError, "Can't put together a weighted and an unweighted component (%s)" % cnames
                            is_w = 0
                        else: # case for unweighted events for MC
                            if (is_w==1): raise RuntimeError, "Can't put together a weighted and an unweighted component (%s)" % cnames
                            is_w = 0;
                            # total_w += n_count ## ROOT histo version
                            total_w += nUnweightedEvents
                            scale = "(%s)" % field[2]
                    else:
                        useHistoForWeight = False
                        maxGenWgt = None
                        if options.weight and len(options.maxGenWeightProc):
                            for procRegExp,tmp_maxGenWgt in options.maxGenWeightProc:
                                if re.match(procRegExp,pname):
                                    #tmp_maxGenWgt is a string, convert to float
                                    maxGenWgt = abs(float(tmp_maxGenWgt))
                                    log10_maxGenWgt = math.log10(maxGenWgt)
                                    useHistoForWeight = True
                                    ROOT.gEnv.SetValue("TFile.AsyncReading", 1);
                                    tmp_rootfile = ROOT.TXNetFile(rootfile+"?readaheadsz=65535")
                                    histo_sumgenweight = tmp_rootfile.Get('distrGenWeights')
                                    if not histo_sumgenweight:
                                        raise RuntimeError, "Can't get histogram distrGenWeights from %s.\nMake sure it is available when using option --max-genWeight-procs" % rootfile
                                    minBin = histo_sumgenweight.GetXaxis().FindFixBin(-log10_maxGenWgt)
                                    maxBin = histo_sumgenweight.GetXaxis().FindFixBin(log10_maxGenWgt)
                                    n_sumgenweight = histo_sumgenweight.Integral(minBin,maxBin)
                                    tmp_rootfile.Close()

                        ## get the counts from the histograms instead of pickle file (smart, but extra load for ROOT from EOS it seems)
                        # ROOT.gEnv.SetValue("TFile.AsyncReading", 1);
                        # tmp_rootfile = ROOT.TXNetFile(rootfile+"?readaheadsz=65535")
                        # histo_count        = tmp_rootfile.Get('Count')
                        # histo_sumgenweight = tmp_rootfile.Get('SumGenWeights')
                        # n_count        = histo_count       .GetBinContent(1)
                        # n_sumgenweight = (histo_sumgenweight.GetBinContent(1) if histo_sumgenweight else n_count)
                        # tmp_rootfile.Close()

                        #do always read this file as well, not to modify the following 'if' statement
                        pckobj  = pickle.load(open(pckfile,'r'))
                        counters = dict(pckobj)                            
                        # if ( n_count != n_sumgenweight ) and options.weight: ## ROOT histo version
                        ## this needed for now...
                        if ('Sum Weights' in counters) and options.weight:
                            if (is_w==0): raise RuntimeError, "Can't put together a weighted and an unweighted component (%s)" % cnames
                            is_w = 1; 
                            if useHistoForWeight:
                                total_w += n_sumgenweight ## ROOT histo version
                                scale = "genWeight*(%s)*(abs(genWeight)<%s)" % (field[2], str(maxGenWgt))
                            else:                            
                                total_w += counters['Sum Weights']
                                scale = "genWeight*(%s)" % field[2]
                        elif not options.weight:
                            scale = "1"
                            total_w = 1
                            if (is_w==1): raise RuntimeError, "Can't put together a weighted and an unweighted component (%s)" % cnames
                            is_w = 0
                        else:
                            if (is_w==1): raise RuntimeError, "Can't put together a weighted and an unweighted component (%s)" % cnames
                            is_w = 0;
                            # total_w += n_count ## ROOT histo version
                            total_w += counters['All Events']
                            scale = "(%s)" % field[2]
                    # closes if options.noHeppyTree
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
                    print "Poorly formatted line: ", field
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
                iFile += 1
            ### END CHECKPOINT
            if to_norm: 
                print ">>> Total sumgenweights = %s (process = %s)" % (str(total_w),pname)
                for tty in ttys: 
                    if options.weight: 
                        tty.setScaleFactor("%s*%g" % (scale, 1000.0/total_w))
                    else: 
                        tty.setScaleFactor(1)
        #if len(self._signals) == 0: raise RuntimeError, "No signals!"
        #if len(self._backgrounds) == 0: raise RuntimeError, "No backgrounds!"
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
        else: raise RuntimeError, "Can't get option %s for undefined process %s" % (name,process)
    def setProcessOption(self,process,name,value):
        if process in self._allData:
            return self._allData[process][0].setOption(name,value)
        elif process in self._optionsOnlyProcesses:
            self._optionsOnlyProcesses[process][name] = value
        else: raise RuntimeError, "Can't set option %s for undefined process %s" % (name,process)
    def getScales(self,process):
        return [ tty.getScaleFactor() for tty in self._allData[process] ] 
    def setScales(self,process,scales):
        for (tty,factor) in zip(self._allData[process],scales): tty.setScaleFactor(factor,mcCorrs=False)
    def getYields(self,cuts,process=None,nodata=False,makeSummary=False,noEntryLine=False):
        ## first figure out what we want to do
        tasks = []
        for key,ttys in self._allData.iteritems():
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
        ret = dict([ (k,mergeReports(v)) for k,v in mergemap.iteritems() ])

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
            for (key,val) in ret.iteritems():
                if key != 'data':
                    if self._isSignal[key]: allSig.append(ret[key])
                    else: allBg.append(ret[key])
            if self._signals and not ret.has_key('signal') and len(allSig) > 0:
                ret['signal'] = mergeReports(allSig)
            if self._backgrounds and not ret.has_key('background') and len(allBg) > 0:
                ret['background'] = mergeReports(allBg)
        return ret
    def getPlotsRaw(self,name,expr,bins,cut,process=None,nodata=False,makeSummary=False,closeTreeAfter=False):
        return self.getPlots(PlotSpec(name,expr,bins,{}),cut,process,nodata,makeSummary,closeTreeAfter)
    def getPlots(self,plotspec,cut,process=None,nodata=False,makeSummary=False,closeTreeAfter=False):
        ret = { }
        allSig = []; allBg = []
        tasks = []
        for key,ttys in self._allData.iteritems():
            if key == 'data' and nodata: continue
            if process != None and key != process: continue
            for tty in ttys:
                tasks.append((key,tty,plotspec,cut,None,closeTreeAfter))
        if self._options.splitFactor > 1 or  self._options.splitFactor == -1:
            tasks = self._splitTasks(tasks)
        retlist = self._processTasks(_runPlot, tasks, name="plot "+plotspec.name)
        ## then gather results with the same process
        mergemap = {}
        for (k,v) in retlist: 
            if k not in mergemap: mergemap[k] = []
            mergemap[k].append(v)
        ## and finally merge them
        ret = dict([ (k,mergePlots(plotspec.name+"_"+k,v)) for k,v in mergemap.iteritems() ])

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
            allSig = [v for k,v in ret.iteritems() if k != 'data'and self._isSignal[k] == True  ]
            allBg  = [v for k,v in ret.iteritems() if k != 'data'and self._isSignal[k] == False ]
            if self._signals and not ret.has_key('signal') and len(allSig) > 0:
                ret['signal'] = mergePlots(plotspec.name+"_signal", allSig)
                ret['signal'].summary = True
            if self._backgrounds and not ret.has_key('background') and len(allBg) > 0:
                ret['background'] = mergePlots(plotspec.name+"_background",allBg)
                ret['background'].summary = True

        #print "DONE getPlots at %.2f" % (0.001*(long(ROOT.gSystem.Now()) - _T0))
        return ret
    def prepareForSplit(self):
        ttymap = {}
        for key,ttys in self._allData.iteritems():
            for tty in ttys:
                if not tty.hasEntries(): 
                    #print "For tty %s/%s, I don't have the number of entries" % (tty._name, tty._cname)
                    ttymap[id(tty)] = tty
        if len(ttymap):
            retlist = self._processTasks(_runGetEntries, ttymap.items(), name="GetEntries")
            for ttid, entries in retlist:
                ttymap[ttid].setEntries(entries)
    def applyCut(self,cut):
        tasks = []; revmap = {}
        for key,ttys in self._allData.iteritems():
            for tty in ttys:
                revmap[id(tty)] = tty
                tasks.append( (id(tty), tty, cut, None) )
        if self._options.splitFactor > 1 or self._options.splitFactor == -1:
            tasks = self._splitTasks(tasks)
        retlist = self._processTasks(_runApplyCut, tasks, name="apply cut "+cut)
        if self._options.splitFactor > 1 or self._options.splitFactor == -1:
            aggregated = {}
            for ttid, elist in retlist:
                if ttid not in aggregated: aggregated[ttid] = elist
                else:                      aggregated[ttid].Add(elist)
            retlist = aggregated.items()
        for ttid, elist in retlist:
            tty = revmap[ttid]
            tty.applyCutAndElist(cut, elist)
    def clearCut(self):
        for key,ttys in self._allData.iteritems():
            for tty in ttys:
                tty.clearCut() 
    def prettyPrint(self,reports,makeSummary=True):
        allSig = []; allBg = []
        for key in reports:
            if key != 'data':
                if self._isSignal[key]: allSig.append((key,reports[key]))
                else: allBg.append((key,reports[key]))
        allSig.sort(key = lambda (n,v): self._rank[n])
        allBg.sort( key = lambda (n,v): self._rank[n])
        table = allSig + allBg
        if makeSummary:
            if len(allSig)>1:
                table.append(('ALL SIG',mergeReports([v for n,v in allSig])))
            if len(allBg)>1:
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

        if self._options.txtfmt == "text":
            print "CUT".center(clen),
            for h,r in table: 
                if len("   "+h) <= fmtlen:
                    print ("   "+h).center(fmtlen),
                elif len(h) <= fmtlen:
                    print h.center(fmtlen),
                else:
                    print h[:fmtlen],
            print ""
            print "-"*((fmtlen+1)*len(table)+clen)
            for i,(cut,dummy) in enumerate(table[0][1]):
                print cfmt % cut,
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
                    if self._options.weight and nev < 1000: print ( nfmtS if nev > 0.2 else nfmtX) % toPrint,
                    else                                  : print nfmtL % toPrint,
                print ""
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
                    print "%s%s%s" % (k,sep,sep.join(ytxt.split()))
                print ""

    def _getYields(self,ttylist,cuts):
        return mergeReports([tty.getYields(cuts) for tty in ttylist])
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
            else: raise RuntimeError, 'Incorrect plotmergemap format: %s'%m
            to = to.strip()
            for k in fro.split(","):
                relist.append((re.compile(k.strip()+"$"), to))
    def compilePlotScaleMap(self,inlist,relist):
        for m in inlist:
            dset,scale = m.split("=")
            if dset[-1] == "*": dset = dset[:-1]
            else: raise RuntimeError, 'Incorrect plotscalemap format: %s'%m
            relist.append((re.compile(dset.strip()+"$"),float(scale)))
    def regroupReports(self,pmap,regexp):
        patt, to = regexp
        mergemap={}
        for (k,v) in pmap.items():
            k2 = k
            if k2 != to and re.match(patt,k2): k2 = to
            if k2 not in mergemap: mergemap[k2]=[]
            mergemap[k2].append(v)
        return dict([ (k,mergeReports(v)) for k,v in mergemap.iteritems() ])
    def regroupPlots(self,pmap,regexp,pspec):
        patt, to = regexp
        mergemap={}
        for (k,v) in pmap.items():
            k2 = k
            if k2 != to and re.match(patt,k2): k2 = to
            if k2 not in mergemap: mergemap[k2]=[]
            mergemap[k2].append(v)
        return dict([ (k,mergePlots(pspec.name+"_"+k,v)) for k,v in mergemap.iteritems() ])
    def stylePlot(self,process,plot,pspec,mayBeMissing=False):
        if process in self._allData:
            for tty in self._allData[process]: 
                tty._stylePlot(plot,pspec)
                break
        elif process in self._optionsOnlyProcesses:
            opts = self._optionsOnlyProcesses[process]
            stylePlot(plot, pspec, lambda key,default : opts[key] if key in opts else default)
        elif not mayBeMissing:
            raise KeyError, "Process %r not found" % process
    def _processTasks(self,func,tasks,name=None):
        #timer = ROOT.TStopwatch()
        #print "Starting job %s with %d tasks, %d threads" % (name,len(tasks),self._options.jobs)
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
                newtasks_wsize.sort(key = lambda (size,task) : size, reverse = True)
                #for s,t in newtasks_wsize: print "\t%9d %s/%s %s" % (s,t[1]._name, t[1]._cname, t[-1])
            newtasks = [ task for (size,task) in newtasks_wsize ]
        #print "New task list has %d entries; actual split factor %.2f" % (len(newtasks), len(newtasks)/float(len(tasks)))
        return newtasks


def addMCAnalysisOptions(parser,addTreeToYieldOnesToo=True):
    if addTreeToYieldOnesToo: addTreeToYieldOptions(parser)
    parser.add_option("-j", "--jobs",           dest="jobs", type="int", default=0, help="Use N threads");
    parser.add_option("--split-factor",         dest="splitFactor", type="int", default=0, help="Use N chunks per sample (-1 means to use the same as what passed to -j, which appears to work well in the average case)");
    #parser.add_option("--split-dynamic",         dest="splitDynamic", action="store_true", default=True, help="Make the splitting dynamic (reduce the chunks for small samples)");
    parser.add_option("--split-static",         dest="splitDynamic", action="store_false", default=True, help="Make the splitting dynamic (reduce the chunks for small samples)");
    #parser.add_option("--split-sort",         dest="splitSort", action="store_true", default=True, help="Make the splitting dynamic (reduce the chunks for small samples)");
    parser.add_option("--split-nosort",         dest="splitSort", action="store_false", default=True, help="Make the splitting dynamic (reduce the chunks for small samples)");
    parser.add_option("-P", "--path", dest="path", action="append", type="string", default=[], help="Path to directory with input trees and pickle files. Can supply multiple paths which will be searched in order. (default: ./") 
    parser.add_option("--RP", "--remote-path",   dest="remotePath",  type="string", default=None,      help="path to remote directory with trees, but not other metadata (default: same as path)") 
    parser.add_option("-p", "--process", dest="processes", type="string", default=[], action="append", help="Processes to print (comma-separated list of regexp, can specify multiple ones)");
    parser.add_option("--pg", "--pgroup", dest="premap", type="string", default=[], action="append", help="Group proceses into one. Syntax is '<newname> := (comma-separated list of regexp)', can specify multiple times. Note that it is applied _before_ -p, --sp and --xp");
    parser.add_option("--xf", "--exclude-files", dest="filesToExclude", type="string", default=[], action="append", help="Files to exclude (comma-separated list of regexp, can specify multiple ones)");
    parser.add_option("--xp", "--exclude-process", dest="processesToExclude", type="string", default=[], action="append", help="Processes to exclude (comma-separated list of regexp, can specify multiple ones)");
    parser.add_option("--sf", "--swap-files", dest="filesToSwap", type="string", default=[], nargs=2, action="append", help="--swap-files X Y uses file Y instead of X in the MCA");
    parser.add_option("--sp", "--signal-process", dest="processesAsSignal", type="string", default=[], action="append", help="Processes to set as signal (overriding the '+' in the text file)");
    parser.add_option("--float-process", "--flp", dest="processesToFloat", type="string", default=[], action="append", help="Processes to set as freely floating (overriding the 'FreeFloat' in the text file; affects e.g. mcPlots with --fitData)");
    parser.add_option("--fix-process", "--fxp", dest="processesToFix", type="string", default=[], action="append", help="Processes to set as not freely floating (overriding the 'FreeFloat' in the text file; affects e.g. mcPlots with --fitData)");
    parser.add_option("--peg-process", dest="processesToPeg", type="string", default=[], nargs=2, action="append", help="--peg-process X Y make X scale as Y (equivalent to set PegNormToProcess=Y in the mca.txt)");
    parser.add_option("--scale-process", dest="processesToScale", type="string", default=[], nargs=2, action="append", help="--scale-process X Y make X scale by Y (equivalent to add it in the mca.txt)");
    parser.add_option("--process-norm-syst", dest="processesToSetNormSystematic", type="string", default=[], nargs=2, action="append", help="--process-norm-syst X Y sets the NormSystematic of X to be Y (for plots, etc. Overrides mca.txt)");
    parser.add_option("--AP", "--all-processes", dest="allProcesses", action="store_true", help="Include also processes that are marked with SkipMe=True in the MCA.txt")
    parser.add_option("--use-cnames",  dest="useCnames", action="store_true", help="Use component names instead of process names (for debugging)")
    parser.add_option("--project", dest="project", type="string", help="Project to a scenario (e.g 14TeV_300fb_scenario2)")
    parser.add_option("--plotgroup", dest="plotmergemap", type="string", default=[], action="append", help="Group plots into one. Syntax is '<newname> := (comma-separated list of regexp)', can specify multiple times. Note it is applied after plotting.")
    parser.add_option("--scaleplot", dest="plotscalemap", type="string", default=[], action="append", help="Scale plots by this factor (before grouping). Syntax is '<newname> := (comma-separated list of regexp)', can specify multiple times.")
    parser.add_option("-t", "--tree",          dest="tree", default='treeProducerWMass', help="Pattern for tree name");
    parser.add_option("--fom", "--figure-of-merit", dest="figureOfMerit", type="string", default=[], action="append", help="Add this figure of merit to the output table (S/B, S/sqrB, S/sqrSB)")
    parser.add_option("--max-genWeight-procs", dest="maxGenWeightProc", type="string", nargs=2, action="append", default=[], help="maximum genWeight to be used for a given MC process (first value is a regular expression for the process, second is the max weight). This option effectively applies a cut on genWeight, and also modifies the sum of genweights. Can be specified more than once for different processes. This will cut away events with larger weights");
    parser.add_option("--clip-genWeight-toMax",         dest="clipGenWeightToMax", action="store_true", default=False, help="It only works with --nanoaod-tree when using --max-genWeight-procs, setting large weights to the max instead of rejecting the event");
    parser.add_option("--no-heppy-tree",         dest="noHeppyTree", action="store_true", default=False, help="Set to true to read root files when they were not made with Heppy (different convention for path names, might need to be adapted)");
    parser.add_option("--nanoaod-tree",         dest="nanoaodTree", action="store_true", default=False, help="Set to true to read root files from nanoAOD");
    parser.add_option("--filter-proc-files", dest="filterProcessFiles", type="string", nargs=2, action="append", default=[], help="Can use this option to override second field on each process line in MCA file, so to select few files without modifying the MCA file (e.g. for tests). E.g. --filter-proc-files 'W.*' '.*_12_.*' to only use files with _12_ in their name. Only works with option --nanoaod-tree");
    parser.add_option("--use-tchain",         dest="useTChain", action="store_true", default=False, help="Use TChain to manage all files per process (should be faster). Only implemented when using --nanoaod-tree, and doesn' work with --split-factor N and N > 1");

if __name__ == "__main__":
    from optparse import OptionParser
    parser = OptionParser(usage="%prog [options] tree.root cuts.txt")
    addMCAnalysisOptions(parser)
    (options, args) = parser.parse_args()
    if not options.path and not options.nanoaodTree: options.path = ['./']
    tty = TreeToYield(args[0],options) if ".root" in args[0] else MCAnalysis(args[0],options)
    cf  = CutsFile(args[1],options)
    for cutFile in args[2:]:
        temp = CutsFile(cutFile,options)
        for cut in temp.cuts():
            cf.add(cut[0],cut[1])
    report = tty.getYields(cf)#, process=options.process)
    tty.prettyPrint(report)
