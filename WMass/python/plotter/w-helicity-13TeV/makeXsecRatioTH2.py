#!/bin/env python

# script to make ratio of two TH2 histograms for cross section
# variation/nominal. If both up and down exists, only Up is used

###################
# >>>>>>>>> WARNING <<<<<<<<<<<< 
###################
# Apparently it works in release CMSSW_10_2_0_pre4
# while in CMSSW_8_0_25 it produces a segmentatio fault from the setTDRStyle()
# Maybe it is just due to the root version

################################
# Exaples
################################


################################
################################


import ROOT, os, sys, re, array, math
import time

sys.path.append(os.getcwd() + "/plotUtils/")
from utility import *

ROOT.gROOT.SetBatch(True)
        
if __name__ == "__main__":
            
    from optparse import OptionParser
    parser = OptionParser(usage='%prog [options] file1')
    parser.add_option('-o','--outdir',      dest='outdir',      default='', type='string', help='output directory to save things')
    #parser.add_option('-f','--outfilename', dest='outfilename', default='', type='string', help='Name of output file to save results')
    #parser.add_option('-n','--outhistname', dest='outhistname', default='', type='string', help='Name of output histogram saved in output file. If FILE is used, take same name as the output file, removing the extension')
    parser.add_option(     '--hname'     ,  dest='hname'     ,  default='gen_ptl1_absetal1_preFSR', type='string', help='Prefix for histogram name')
    parser.add_option('-c','--charge'    ,  dest='charge'    ,  default='plus', type='string', help='Charge to select histogram name')
    parser.add_option('-f','--flavor'    ,  dest='flavor'    ,  default='mu'  , type='string', help='Select channel (mu|el) to select histogram name')
    parser.add_option('-x','--xAxisTitle',  dest='xAxisTitle',  default='', type='string', help='X axis title. If not given, use the one from hist1')
    parser.add_option('-y','--yAxisTitle',  dest='yAxisTitle',  default='', type='string', help='Y axis title. If not given, use the one from hist1')
    parser.add_option('-z','--zAxisTitle',  dest='zAxisTitle',  default='', type='string', help='Z axis title. If not given, use the one from hist1')
    parser.add_option('-t','--histTitle',   dest='histTitle',   default='', type='string', help='Title to assign to output histogram. It is used as a label for the canvas')
    parser.add_option('-r','--ratioRange',  dest='ratioRange',  default=(0.99, 1.01),type="float", nargs=2, help="Min and max for the ratio in the plot")
    parser.add_option(     '--xRange'     , dest='xRange', default=(0,-1), type='float', nargs=2, help='Select range for X axis to plot. If min > max, the option is neglected')
    parser.add_option(     '--yRange'     , dest='yRange', default=(0,-1), type='float', nargs=2, help='Select range for Y axis to plot. If min > max, the option is neglected')
    (options, args) = parser.parse_args()

    if len(sys.argv) < 1:
        parser.print_usage()
        quit()

    cmssw_version = os.environ['CMSSW_VERSION']
    if int(cmssw_version.split('_')[1]) < 10:
        print ""
        print ">>>> WARNING:" 
        print "current version of CMSSW ({version}) might lead to segmentation fault when running this script!".format(version=cmssw_version) 
        print "Or, the Z axis range might not be set properly"
        print "If you experience any of these issues, try running the script from a release CMSSW_X_Y_Z with X > 8 (10_2 suggested)"
        print "Versions with Y > 2 seems to work, but lead to some warnings and print garbage on stdout"
        print "Continuing in 3 seconds ..."
        print ""
        time.sleep(3)

    f1 = args[0]

    ROOT.TH1.SetDefaultSumw2()

    print ""

    if options.charge not in ["plus", "minus"]:
        print "Error: charge must be plus or minus. Exit"
        quit()
    if options.flavor not in ["el", "mu"]:
        print "Error: flavor must be el or mu. Exit"
        quit()

    if options.outdir:
        outname = options.outdir
        addStringToEnd(outname,"/",notAddIfEndswithMatch=True)
        outname = outname + options.charge + "/"
        createPlotDirAndCopyPhp(outname)
    else:
        print "Error: you should specify an output folder using option -o <name>. Exit"
        quit()
    # if not options.outfilename:
    #     print "Error: you should specify an output file name using option -f <name>. Exit"
    #     quit()
    # if not options.outhistname:
    #     print "Error: you should specify an output histogram name using option -n <name>. "
    #     print "If FILE is used, take same name as the output file, removing the extension"
    #     print "Exit"
    #     quit()


    #if options.outhistname == "FILE":
    #    options.outhistname = options.outfilename.split('.')[0]

    hratio = 0
    hnomi = 0
    cname = "{pf}_W{c}_{f}_central".format(pf=options.hname, c=options.charge, f=options.flavor)

    # file 1
    tf = ROOT.TFile.Open(f1)        
    hnomi =   tf.Get(cname)
    if (hnomi == 0):
        print "Error: could not retrieve %s from input file %s. Exit" % (hname,f1)
        quit()
    else:
        hnomi.SetDirectory(0)

    xMin = options.xRange[0]
    xMax = options.xRange[1]
    yMin = options.yRange[0]
    yMax = options.yRange[1]

    hnames = []
    for i in range(1,61):
        hnames.append("{pf}_W{c}_{f}_pdf{ipdf}".format(pf=options.hname, c=options.charge, f=options.flavor, ipdf=i))
    for i in range(1,11):
        hnames.append("{pf}_W{c}_{f}_muR{iqcd}Up".format(pf=options.hname, c=options.charge, f=options.flavor, iqcd=i))
        hnames.append("{pf}_W{c}_{f}_muRmuF{iqcd}Up".format(pf=options.hname, c=options.charge, f=options.flavor, iqcd=i))
        hnames.append("{pf}_W{c}_{f}_muRmuF{iqcd}Up".format(pf=options.hname, c=options.charge, f=options.flavor, iqcd=i))
    hnames.append("{pf}_W{c}_{f}_mW_80470".format(pf=options.hname, c=options.charge, f=options.flavor))
    hnames.append("{pf}_W{c}_{f}_alphaSUp".format(pf=options.hname, c=options.charge, f=options.flavor))
    
    canvas2D = ROOT.TCanvas("canvas","",700,700)

    for n in hnames:

        hnum = tf.Get(n)
        if (hnum == 0):
            print "Error: could not retrieve %s from input file %s. Exit" % (n,f1)
            quit()

        hratio = hnum.Clone("ratio_" + n.replace("{pf}_".format(pf=options.hname),""))
        hratio.SetTitle(options.histTitle)
        hratio.Divide(hnomi)
    
        if options.xAxisTitle: hratio.GetXaxis().SetTitle(options.xAxisTitle)    
        if options.yAxisTitle: hratio.GetYaxis().SetTitle(options.yAxisTitle)    
        if options.zAxisTitle: hratio.GetZaxis().SetTitle(options.zAxisTitle)    
        xAxisTitle = hratio.GetXaxis().GetTitle()
        yAxisTitle = hratio.GetYaxis().GetTitle()
        zAxisTitle = hratio.GetZaxis().GetTitle()

        if xMax > xMin and not "::" in xAxisTitle: 
            xAxisTitle = xAxisTitle + "::" + str(options.xRange[0]) + "," + str(options.xRange[1])
        if yMax > yMin and not "::" in yAxisTitle: 
            yAxisTitle = yAxisTitle + "::" + str(options.yRange[0]) + "," + str(options.yRange[1])
                
    # print "xAxisTitle = " + xAxisTitle
    # print "yAxisTitle = " + yAxisTitle
    # print "zAxisTitle = " + zAxisTitle


        adjustSettings_CMS_lumi()
    # the axis name can be used to set the range if it is in the format "name::min,maz"
    # if this is not already the case, use the selected range from the input option
        if not "::" in zAxisTitle:  
            zAxisTitle = zAxisTitle + "::" + str(options.ratioRange[0]) + "," + str(options.ratioRange[1])
        drawCorrelationPlot(hratio,xAxisTitle,yAxisTitle,zAxisTitle,
                            hratio.GetName(),"ForceTitle",outname,0,0,False,False,False,1,palette=55,passCanvas=canvas2D)
    

    tf.Close()
 
    print ""
    print "DONE" 
    print ""

                               
         
