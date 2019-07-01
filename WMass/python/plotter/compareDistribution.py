#!/usr/bin/env python                                                                                                                                                        
#from shutil import copyfile
import re, sys, os, os.path, ROOT
from array import array
#import numpy as np

from optparse import OptionParser
parser = OptionParser(usage="%prog [options] f1.root h1 f2.root h2")
parser.add_option("-o", "--outdir",    dest="outdir", type="string", default="./", help="Output folder (current one as default)");
parser.add_option("-n", "--name",      dest="name",   type="string", default="", help="Name for output root file");
parser.add_option("-x", "--xname",    dest="xname",   type="string", default="", help="Name for x axis (h1 default if empty)");
parser.add_option("-y", "--yname",    dest="yname",   type="string", default="", help="Name for y axis (h2 default if empty)");
parser.add_option(      "--xrange",    dest="xrange",   type="string", default="", help="x axis range (comma separated pair of floats)");
parser.add_option(      "--yrange",    dest="yrange",   type="string", default="", help="y axis range (comma separated pair of floats)");
parser.add_option(      "--allow-different-nbinsX",    dest="allowDifferentNbinsX", default=False, action='store_true', help="Allow different number of bins on second histogram. Will select the first N bins");
parser.add_option(      "--dummyError",    dest="dummyError",   type="string", default="", help="Set this dummy number as error bars");
(options, args) = parser.parse_args()

ROOT.TH1.SetDefaultSumw2()

if len(sys.argv) < 4:
    parser.print_usage()
    quit()

f1 = args[0]
name1 = args[1]
f2 = args[2]
name2 = args[3]
    
# manage output folder
outdir = options.outdir
if not outdir.endswith('/'): outdir += "/"
if outdir != "./":
    if not os.path.exists(outdir):
        print "Creating folder", outdir
        os.system("mkdir -p " + outdir)

if os.path.exists("/afs/cern.ch"): 
    os.system("cp /afs/cern.ch/user/g/gpetrucc/php/index.php "+ outdir)

tf1 = ROOT.TFile.Open(f1,"READ")
if not tf1 or not tf1.IsOpen():
    raise RuntimeError('Unable to open file {fn}'.format(fn=f1))
h1 = tf1.Get(name1)
if not h1:
    raise RuntimeError('Unable to read histogram {n1} from file 1'.format(n1=name1))
h1.SetDirectory(0)
tf1.Close()

# open output file
tf2 = ROOT.TFile(f2,'READ')
if not tf2 or not tf2.IsOpen():
    raise RuntimeError('Unable to open file {fn}'.format(fn=f2))
h2 = tf2.Get(name2)
if not h2:
    raise RuntimeError('Unable to read histogram {n2} from file 2'.format(n2=name2))
h2.SetDirectory(0)
tf2.Close()

integralRatio = h1.Integral() / h2.Integral()

h1.Scale(1./h1.Integral())
h2.Scale(1./h2.Integral())
h2tmp = 0
if h1.GetNbinsX() != h1.GetNbinsY() and options.allowDifferentNbinsX:
    xbinedges = [float(h1.GetXaxis().GetBinLowEdge(i)) for i in range(1,h1.GetNbinsX()+2)]
    h2tmp = ROOT.TH1D(h2.GetName()+"_tmp",h2.GetTitle(),len(xbinedges)-1, array('d',xbinedges))
    for i in range(1,h2tmp.GetNbinsX()+1):
        h2tmp.SetBinContent(i, h2.GetBinContent(i))
        h2tmp.SetBinError(i, h2.GetBinError(i))

hratio = h1.Clone("ratio")
divideRes = True
if options.allowDifferentNbinsX:
    divideRes = hratio.Divide(h2tmp)
else:
    divideRes = hratio.Divide(h2)
if not divideRes:
    raise RuntimeError('Error when dividing histograms')
hratio.GetXaxis().SetTitle(options.xname if len(options.xname) else h1.GetXaxis().GetTitle() )
hratio.GetYaxis().SetTitle(options.yname if len(options.yname) else h1.GetYaxis().GetTitle() )
if len(options.xrange): 
    xmin,xmax = options.xrange.split(',')
    hratio.GetXaxis().SetRangeUser(float(xmin),float(xmax))
if len(options.yrange): 
    ymin,ymax = options.yrange.split(',')
    hratio.GetYaxis().SetRangeUser(float(ymin),float(ymax))

if len(options.dummyError):
    err = float(options.dummyError)
    for i in range(1,1+hratio.GetNbinsX()):
        hratio.SetBinError(i,err)

canvasName =  options.name if  options.name else ("comparison_" + h1.GetName() + "__" + h2.GetName())
canvas = ROOT.TCanvas("canvas","",700,600)
canvas.SetTickx(1)
canvas.SetTicky(1)
canvas.SetGridy(1)
canvas.cd()
canvas.SetLeftMargin(0.16)
canvas.SetRightMargin(0.06)
canvas.SetBottomMargin(0.15)
canvas.cd()

hratio.SetStats(0)
hratio.SetTitle("")
hratio.SetFillColor(0)
hratio.SetLineColor(ROOT.kBlack)
hratio.SetLineWidth(2)
hratio.Draw("HE")
canvas.RedrawAxis("sameaxis")
canvas.SaveAs(outdir + canvasName + ".png")


print ""
print "Normalization ratio: h1/h2 = ", str(integralRatio)
print ""

for i in range(1,hratio.GetNbinsX()):
    print "bin %d: ratio = %.5f +/- %.5f" % (i, hratio.GetBinContent(i), hratio.GetBinError(i))
print ""

