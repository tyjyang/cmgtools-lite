#!/usr/bin/env python
#from mcAnalysis import *
#from CMGTools.WMass.plotter.mcEfficiencies import *
from CMGTools.WMass.plotter.mcPlots import *
import itertools

if "/fakeRate_cc.so" not in ROOT.gSystem.GetLibraries():
    compileMacro("src/CMGTools/WMass/python/plotter/fakeRate.cc")

def makeDataSub(report,mca):
    data_sub      = report['data'].Clone(report['data'].GetName()+'_sub')
    data_sub_syst = report['data'].Clone(report['data'].GetName()+'_sub_syst')
    for p in mca.listBackgrounds():
        if p not in report: continue
        b = report[p]
        data_sub.Add(b, -1.0)
        data_sub_syst.Add(b, -1.0)
        syst = mca.getProcessOption(p,'NormSystematic',0.) 
        #print "subtracting background %s from data with systematic %r" % (p,syst)
        if syst <= 0: continue
        if "TH1" in b.ClassName():
            for bx in xrange(1,b.GetNbinsX()+1):
                data_sub_syst.SetBinError(bx, hypot(data_sub_syst.GetBinError(bx), syst * b.GetBinContent(bx)))
        elif "TH2" in b.ClassName():
            for (bx,by) in itertools.product(range(1,b.GetNbinsX()+1), range(1,b.GetNbinsY()+1)):
                data_sub_syst.SetBinError(bx, by, hypot(data_sub_syst.GetBinError(bx, by), syst * b.GetBinContent(bx, by)))
        elif "TH3" in b.ClassName():
            for (bx,by,bz) in itertools.product(range(1,b.GetNbinsX()+1), range(1,b.GetNbinsY()+1), range(1,b.GetNbinsZ()+1)):
                data_sub_syst.SetBinError(bx, by, bz, hypot(data_sub_syst.GetBinError(bx, by, bz), syst * b.GetBinContent(bx, by, bz)))
    report['data_sub']      = data_sub
    report['data_sub_syst'] = data_sub_syst

def makeEff(mca,cut,idplot,xvarplot,mainOptions=None):
    import copy
    is2D = (":" in xvarplot.expr.replace("::","--"))
    options = copy.copy(idplot.opts)
    options.update(xvarplot.opts)
    mybins = copy.copy(xvarplot.bins)
    if xvarplot.bins[0] == "[":
        mybins += "*[-0.5,0.5,1.5]"
    else:
        mybins += ",2,-0.5,1.5"
    print "making eff for idplot = ",idplot.name,"  ", idplot.expr, " vs ",xvarplot.expr
    pspec = PlotSpec("%s_vs_%s"  % (idplot.name, xvarplot.name), 
                     "%s:%s" % (idplot.expr,xvarplot.expr),
                     mybins,
                     options) 
    report = mca.getPlots(pspec,cut,makeSummary=True)
    if 'signal' in report and 'background' in report:
        report['total'] = mergePlots(pspec.name+"_total", [ report[s] for s in ('signal','background') ] )
    if 'data' in report and 'background' in report:
        makeDataSub(report, mca)
    return report


if __name__ == "__main__":
    from optparse import OptionParser
    parser = OptionParser(usage="%prog [options] mc.txt cuts.txt ids.txt xvars.txt")
    parser.add_option("--fitVar",  dest="fitVar",  default="",  type="string", help="Variable to fit (used to be mT, now it is usually eta)")
    addPlotMakerOptions(parser,  addAlsoMCAnalysis=True)
    (options, args) = parser.parse_args()

    mca  = MCAnalysis(args[0],options)
    procs = mca.listProcesses()
    cut = CutsFile(args[1],options).allCuts()
    ids   = PlotFile(args[2],options).plots()
    xvars = PlotFile(args[3],options).plots()
    fitvarname = options.fitVar
    
    print "Processes = ",procs

    outname  = options.out if options.out else (args[2].replace(".txt","")+".root")
    if os.path.dirname(outname) != "":
        dirname = os.path.dirname(outname)
        options.printDir = dirname
        if not os.path.exists(dirname):
            os.system("mkdir -p "+dirname)
            if os.path.exists("/afs/cern.ch"): os.system("cp /afs/cern.ch/user/m/mciprian/public/index.php "+dirname)
        # copy mca and cut file to output folder
        os.system("cp %s %s " % (args[0], dirname))
        os.system("cp %s %s " % (args[1], dirname))

            
    outfile  = ROOT.TFile(outname,"RECREATE")
    plotter = PlotMaker(outfile, options)

    vars2d = []
    print "xvars: %s" % " ".join(x.name for x in xvars)
    for xspec in xvars:
        if xspec.name == fitvarname: continue
        for fspec in xvars:
            if fspec.name != fitvarname: continue
            if xspec.bins[0] == "[":
                if fspec.bins[0] == "[":        
                    fbins = fspec.bins
                else:
                    (nbins,fmin,fmax) = map(float, fspec.bins.split(','))
                    fbins = "[" + ",".join(map(str, [fmin+i*(fmax-fmin)/nbins for i in xrange(0,int(nbins+1))]))  + "]"
                bins2d = xspec.bins + "*" + fbins
            elif fspec.bins[0] == "[":
                (nbins,xmin,xmax) = map(float, xspec.bins.split(','))
                xbins = "[" + ",".join(map(str, [xmin+i*(xmax-xmin)/nbins for i in xrange(0,int(nbins+1))])) + "]"
                bins2d = xbins + "*" + fspec.bins
            else:
                bins2d = xspec.bins + "," + fspec.bins
            print "Doing %s_%s: %s"  % (fspec.name, xspec.name, bins2d)
            pspec = PlotSpec("%s_%s"  % (fspec.name, xspec.name), "%s:%s" % (fspec.expr, xspec.expr), bins2d, xspec.opts) 
            pspec.xvar = xspec
            pspec.fvar = fspec
            vars2d.append(pspec) 

    backup = options.globalRebin;
    options.globalRebin = 1
    hists = [ (y,x2.xvar,x2.fvar,x2,makeEff(mca,cut,y,x2,mainOptions=options)) for y in ids for x2 in vars2d ]
    print "N(hists) = %d" % len(hists)
    options.globalRebin = backup

    for (yspec,xspec,fspec,x2d,report) in hists:
        for k,v in report.iteritems(): 
            print "{n}:  {i}".format(n=v.GetName(),i=v.Integral())
            outfile.WriteTObject(v)
    outfile.ls()
    outfile.Close()
    print "Output saved to %s. exiting." % outname
    exit()


