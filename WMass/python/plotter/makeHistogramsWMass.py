#!/usr/bin/env python

# developed from makeShapeCards.py

from CMGTools.WMass.plotter.mcAnalysis import *
import re, sys, os, os.path
import math

def mirrorShape(nominal,alternate,mirror):
    # assumes any regularization (e.g. cropping negative bin content to make it 0) already  happened outside
    # same for normalization

    # FIXME: need to add possible conditions when either y0 or yA are negative
    # while this choice is only relevant for some bins with low stat in MC, it can induce non trivial 
    # effects on constraints of some particular nuisance parameters

    for b in xrange(1,nominal.GetNbinsX()+1):
        y0 = nominal.GetBinContent(b)
        yA = alternate.GetBinContent(b)
        yM = max(0,2*y0-yA)
        mirror.SetBinContent(b, yM)
    return mirror


def rebin2Dto1D(h,funcstring):
    nbins,fname = funcstring.split(':',1)
    func = getattr(ROOT,fname)
    nbins = int(nbins)
    goodname = h.GetName()
    h.SetName(goodname+"_oldbinning")
    newh = ROOT.TH1D(goodname,h.GetTitle(),nbins,0.5,nbins+0.5)
    x = h.GetXaxis()
    y = h.GetYaxis()
    allowed = range(1,nbins+1)
    if 'TH2' not in h.ClassName(): raise RuntimeError, "Calling rebin2Dto1D on something that is not TH2"
    for i in xrange(x.GetNbins()):
        for j in xrange(y.GetNbins()):
            bin = int(func(x.GetBinCenter(i+1),y.GetBinCenter(j+1)))
            if bin not in allowed: raise RuntimeError, "Binning function gives not admissible result"
            newh.SetBinContent(bin,newh.GetBinContent(bin)+h.GetBinContent(i+1,j+1))
            newh.SetBinError(bin,math.hypot(newh.GetBinError(bin),h.GetBinError(i+1,j+1)))
    for bin in range(1,nbins+1):
        if newh.GetBinContent(bin)<0:
            print 'Warning: cropping to zero bin %d in %s (was %f)'%(bin,newh.GetName(),newh.GetBinContent(bin))
            newh.SetBinContent(bin,0)
    newh.SetLineWidth(h.GetLineWidth())
    newh.SetLineStyle(h.GetLineStyle())
    newh.SetLineColor(h.GetLineColor())
    return newh

def unroll2Dto1D(h):
    nbins = h.GetNbinsX() * h.GetNbinsY()
    goodname = h.GetName()
    h.SetName(goodname+"_oldbinning")
    newh = ROOT.TH1D(goodname,h.GetTitle(),nbins,0.5,nbins+0.5)
    if 'TH2' not in h.ClassName(): raise RuntimeError, "Calling unroll2Dto1D on something that is not TH2"
    for i in xrange(h.GetNbinsX()):
        for j in xrange(h.GetNbinsY()):
            bin = 1 + i + j*h.GetNbinsX()
            newh.SetBinContent(bin,h.GetBinContent(i+1,j+1))
            newh.SetBinError(bin,h.GetBinError(i+1,j+1))
    for bin in range(1,nbins+1):
        if newh.GetBinContent(bin)<0:
            print 'Warning: found bin with negative event weight! will set to it though!'
            #print 'Warning: cropping to zero bin %d in %s (was %f)'%(bin,newh.GetName(),newh.GetBinContent(bin))
            #newh.SetBinContent(bin,0)
    return newh

from optparse import OptionParser
parser = OptionParser(usage="%prog [options] mc.txt cuts.txt var bins systs.txt ")
addMCAnalysisOptions(parser)
parser.add_option("-o",   "--out",    dest="outname", type="string", default=None, help="output name") 
parser.add_option("--od", "--outdir", dest="outdir", type="string", default=None, help="output name") 
parser.add_option("-v", "--verbose",  dest="verbose",  default=0,  type="int",    help="Verbosity level (0 = quiet, 1 = verbose, 2+ = more)")
parser.add_option("--asimov", dest="asimov", action="store_true", default=False, help="Asimov")
parser.add_option("--2d-binning-function",dest="binfunction", type="string", default=None, help="Function used to bin the 2D histogram: for now can be None or unroll2Dto1D")
parser.add_option("--infile",dest="infile", type="string", default=None, help="File to read histos from (to reuse the one made with --savefile)")
parser.add_option("--savefile",dest="savefile", type="string", default=None, help="File to save histos to (this has only those produced by getPlotsRaw() )")

(options, args) = parser.parse_args()

# can be deleted I think
#options.weight = True
#options.final  = True
#options.allProcesses  = True

if "/functions_cc.so" not in ROOT.gSystem.GetLibraries():
    compileMacro("src/CMGTools/WMass/python/plotter/functions.cc")

mca  = MCAnalysis(args[0],options)
cuts = CutsFile(args[1],options)

binname = os.path.basename(args[1]).replace(".txt","") if options.outname == None else options.outname
outdir  = options.outdir+"/" if options.outdir else "./"


myout = outdir;
if not os.path.exists(myout): os.mkdir(myout)

report={}

if options.infile!=None:
    infile = ROOT.TFile(myout+binname+".bare.root","read")
    for p in mca.listSignals(True)+mca.listBackgrounds(True)+['data']:
        h = infile.Get(p)
        if h: report[p] = h
else:
    report = mca.getPlotsRaw("x", args[2], args[3], cuts.allCuts(), nodata=options.asimov)

if options.savefile!=None:
    savefile = ROOT.TFile(myout+binname+".bare.root","recreate")
    for n,h in report.iteritems(): savefile.WriteTObject(h,n)
    savefile.Close()

if options.asimov:
    tomerge = []
    for p in mca.listSignals() + mca.listBackgrounds():
        if p in report: tomerge.append(report[p])
    if len(tomerge): report['data_obs'] = mergePlots("x_data_obs", tomerge) 
else:
    report['data_obs'] = report['data'].Clone("x_data_obs") 

allyields = dict([(p,h.Integral()) for p,h in report.iteritems()])
procs = []
for i,s in enumerate(mca.listSignals()):
    if allyields[s] == 0: continue
    procs.append(s)
for i,b in enumerate(mca.listBackgrounds()):
    if allyields[b] == 0: continue
    procs.append(b)

systs = {} # not needed if not creating dummy cards, but keep for now (it was ised for lnN nuisances, which do not need a new template)
systsEnv = {}
for sysfile in args[4:]:
    for line in open(sysfile, 'r'):
        if re.match("\s*#.*", line): continue
        line = re.sub("#.*","",line).strip()
        if len(line) == 0: continue
        field = [f.strip() for f in line.split(':')]
        if len(field) < 4:
            raise RuntimeError, "Malformed line %s in file %s"%(line.strip(),sysfile)
        elif len(field) == 4 or field[4] == "lnN":
            (name, procmap, binmap, amount) = field[:4]
            if re.match(binmap+"$",binname) == None: continue
            if name not in systs: systs[name] = []
            systs[name].append((re.compile(procmap+"$"),amount))
        elif field[4] in ["envelop","shapeOnly","templates","templatesShapeOnly","alternateShape","alternateShapeOnly"] or '2D' in field[4]:
            (name, procmap, binmap, amount) = field[:4]
            if re.match(binmap+"$",binname) == None: continue
            if name not in systs: systsEnv[name] = []
            systsEnv[name].append((re.compile(procmap+"$"),amount,field[4]))
        elif field[4] in ["stat_foreach_shape_bins"]:
            (name, procmap, binmap, amount) = field[:4]
            if re.match(binmap+"$",binname) == None: continue
            if name not in systsEnv: systsEnv[name] = []
            systsEnv[name].append((re.compile(procmap+"$"),amount,field[4],field[5].split(',')))
        else:
            raise RuntimeError, "Unknown systematic type %s" % field[4]
    if options.verbose:
        print "Loaded %d systematics" % len(systs)
        print "Loaded %d envelop systematics" % len(systsEnv)


if options.binfunction:
    newhistos={}
    _to_be_rebinned={}
    for n,h in report.iteritems(): _to_be_rebinned[h.GetName()]=h
    for n,h in _to_be_rebinned.iteritems():
        thisname = h.GetName()
        if(options.binfunction=='unroll2Dto1D'): newhistos[thisname]=unroll2Dto1D(h)
        else: newhistos[thisname]=rebin2Dto1D(h,options.binfunction)
    for n,h in report.iteritems(): report[n] = newhistos[h.GetName().replace('_oldbinning','')]
    allyields = dict([(p,h.Integral()) for p,h in report.iteritems()])
    procs = []
    for i,s in enumerate(mca.listSignals()):
        if allyields[s] == 0: continue
        procs.append(s)
    for i,b in enumerate(mca.listBackgrounds()):
        if allyields[b] == 0: continue
        procs.append(b)

systsEnv2={}
for name in systsEnv.keys():
    modes = [entry[2] for entry in systsEnv[name]]
    for _m in modes:
        if _m!=modes[0]: raise RuntimeError, "Not supported"
    if (any([re.match(x+'.*',modes[0]) for x in ["envelop","shapeOnly"]])): continue # do only this before rebinning
    # we plan to have only templates* or alternateShape* nuisance parameters
    effmap0  = {}
    effmap12 = {}
    mode = ""
    for p in procs:
        effect = "-"
        effect0  = "-"
        effect12 = "-"
        mode = ""
        for entry in systsEnv[name]:
            procmap,amount,mode = entry[:3]
            if re.match(procmap, p):
                effect = float(amount) if mode not in ["templates","templatesShapeOnly","alternateShape", "alternateShapeOnly"] else amount
        if mca._projection != None and effect not in ["-","0","1",1.0,0.0] and type(effect) == type(1.0):
            effect = mca._projection.scaleSyst(name, effect)
        if effect == "-" or effect == "0": 
            effmap0[p]  = "-" 
            effmap12[p] = "-" 
            continue
        
        if mode in ["templates","templatesShapeOnly"]:
            nominal = report[p]
            p0Up = report["%s_%s_Up" % (p, effect)]
            p0Dn = report["%s_%s_Dn" % (p, effect)]
            if not p0Up or not p0Dn: 
                raise RuntimeError, "Missing templates %s_%s_(Up,Dn) for %s" % (p,effect,name)
            p0Up.SetName("%s_%sUp"   % (nominal.GetName(),name))
            p0Dn.SetName("%s_%sDown" % (nominal.GetName(),name))
            if p0Up.Integral()<=0 or p0Dn.Integral()<=0:
                if p0Up.Integral()<=0 and p0Dn.Integral()<=0: raise RuntimeError, 'ERROR: both template variations have negative or zero integral: %s, Nominal %f, Up %f, Down %f'%(p,nominal.Integral(),p0Up.Integral(),p0Dn.Integral())
                print 'Warning: I am going to fix a template prediction that would have negative or zero integral: %s, Nominal %f, Up %f, Down %f'%(p,nominal.Integral(),p0Up.Integral(),p0Dn.Integral())
                for b in xrange(1,nominal.GetNbinsX()+1):
                    y0 = nominal.GetBinContent(b)
                    yA = p0Up.GetBinContent(b) if p0Up.Integral()>0 else p0Dn.GetBinContent(b)
                    yM = y0
                    if (y0 > 0 and yA > 0):
                        yM = y0*y0/yA
                    elif yA == 0:
                        yM = 2*y0
                    if p0Up.Integral()>0: p0Dn.SetBinContent(b, yM)
                    else: p0Up.SetBinContent(b, yM)
                print 'The integral is now: %s, Nominal %f, Up %f, Down %f'%(p,nominal.Integral(),p0Up.Integral(),p0Dn.Integral())
            if mode == 'templatesShapeOnly':
                p0Up.Scale(nominal.Integral()/p0Up.Integral())
                p0Dn.Scale(nominal.Integral()/p0Dn.Integral())
            report[str(p0Up.GetName())[2:]] = p0Up
            report[str(p0Dn.GetName())[2:]] = p0Dn
            effect0  = "1"
            effect12 = "-"
            if mca._projection != None:
                mca._projection.scaleSystTemplate(name,nominal,p0Up)
                mca._projection.scaleSystTemplate(name,nominal,p0Dn)
        elif mode in ["alternateShape", "alternateShapeOnly"]:
            # for k in report:
            #     print "%s : %s" % (k, report[k])
            if options.verbose:
                print "CHECKPOINT: %s   %s" % (p,mode)
            nominal = report[p]
            alternate = report[effect]
            if mca._projection != None:
                mca._projection.scaleSystTemplate(name,nominal,alternate)
            alternate.SetName("%s_%sUp" % (nominal.GetName(),name))
            if mode == "alternateShapeOnly":
                alternate.Scale(nominal.Integral()/alternate.Integral())
            mirror = nominal.Clone("%s_%sDown" % (nominal.GetName(),name))
            mirror = mirrorShape(nominal,alternate,mirror) # this does not alter normalization
            if mode == "alternateShapeOnly":
                # keep same normalization
                mirror.Scale(nominal.Integral()/mirror.Integral())
            if options.verbose:
                print "CHECKPOINT: int(nomi) ", str(nominal.Integral())
                print "CHECKPOINT: int(alte) ", str(alternate.Integral())
                print "CHECKPOINT: int(mirr) ", str(mirror.Integral())
            report[alternate.GetName()] = alternate
            report[mirror.GetName()] = mirror
            effect0  = "1"
            effect12 = "-"
        effmap0[p]  = effect0 
        effmap12[p] = effect12 
    if mode not in ["stat_foreach_shape_bins"]: systsEnv2[name] = (effmap0,effmap12,mode)

systsEnv = {}
systsEnv.update(systsEnv2)

workspace = ROOT.TFile.Open(myout+binname+".input.root", "RECREATE")
for n,h in report.iteritems():
    if options.verbose: print "\t%s (%8.3f events)" % (h.GetName(),h.Integral())
    workspace.WriteTObject(h,h.GetName())
workspace.Close()

print "Wrote to ",myout+binname+".input.root"
# now check goodness of file, if bad returns an exit code different from 0
f = ROOT.TFile.Open(myout+binname+".input.root", "READ")
if f.IsZombie():    
    print 'file is probably corrupted'
    sys.exit(100)
if f.TestBit(ROOT.TFile.kRecovered):
    print 'file is in fishy state, was recovered'
    sys.exit(101)
if f.GetListOfKeys().GetSize()==0:
    print 'file is bad, has no keys'
    sys.exit(102)
print "File is in good state :)"
