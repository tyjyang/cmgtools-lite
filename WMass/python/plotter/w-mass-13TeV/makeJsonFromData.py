#!/bin/env python


################################
# Examples
################################

# python w-mass-13TeV/makeJsonFromData.py -i /data/shared/originalNANO/NanoV8Data/Run2016B/ -o lumiStuff/ -n NanoAOD_json_Run2016B.txt --max-entries 100 -v

# Once the json is made, one can compute the integrated luminosity with this command:
# brilcalc lumi --normtag /cvmfs/cms-bril.cern.ch/cms-lumi-pog/Normtags/normtag_PHYSICS.json -u /fb -i lumiStuff/NanoAOD_json_Run2016B.txt --without-checkjson

# --normtag adds a filter based on latest luminosity calibration: see https://twiki.cern.ch/twiki/bin/view/CMS/TWikiLUM#CurRec
# Running brilcalc requires having brilcalc installed on lxplus: see https://cms-service-lumi.web.cern.ch/cms-service-lumi/brilwsdoc.html#brilcalc

################################
################################


import os, re, array, math
import argparse
## safe batch mode
import sys
args = sys.argv[:]
sys.argv = ['-b']
import ROOT
sys.argv = args
ROOT.gROOT.SetBatch(True)
ROOT.PyConfig.IgnoreCommandLineOptions = True

sys.path.append(os.getcwd() + "/plotUtils/")
from utility import *

if __name__ == "__main__":
            
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('-i','--input', dest='inputdir', default='', type=str, help='input directory with ntuples.')
    parser.add_argument('-o','--outdir', dest='outdir', default='', type=str, help='output directory to save things (should end with /)')
    parser.add_argument('-n','--outfilename', dest='outfilename', default='', type=str, help='Name of output file to save json')
    parser.add_argument('-m','--match',     dest='match',     default='', type=str, help='Filter data files matching this regular expression')
    parser.add_argument(     '--max-entries',     dest='maxentries',     default=0, type=int, help='How many events to process (for tests). Default is 0, which is all')
    parser.add_argument('-v',"--verbose",   dest="verbose",    action="store_true", default=False, help="Print debugging stuff");
    args = parser.parse_args()
    

    ROOT.TH1.SetDefaultSumw2()

    if args.outdir:
        outname = args.outdir
        if not outname.endswith("/"): outname += "/"
        os.system("mkdir -p {outdir}".format(outdir=outname))
    else:
        print("Error: you should specify an output folder using option -o <name>. Exit")
        quit()

    if not args.outfilename:
        print("Error: you should specify an output file name using option -n <name>. Exit")
        quit()

    outfilename = args.outfilename
    outfileFullName = outname + outfilename

    files = []
    for dirpath, dirnames, filenames in os.walk(args.inputdir):
        filtered_fnames = [f for f in filenames if f.endswith(".root")]
        for filename in filtered_fnames:
            files.append(os.path.join(dirpath, filename))
            if args.match != "":
                files = [f for f in files if any(re.match(x, f) for x in args.match.split(','))]

    if args.verbose:
        printLine("-",20)
        print("Reading these files (%d files)" % len(files))
        for i,f in enumerate(files):
            print("%d) %s" % (i,f))

    chain = ROOT.TChain("LuminosityBlocks")
    
    for f in files:
        chain.Add(f)

    entries = chain.GetEntries()
    print("Chain has %d entries" % entries)
    leaves  = chain.GetListOfLeaves()

    # for entry in range(10):
    #     chain.GetEntry(entry)            
    #     print("run = %d,    LS = %d " % (chain.run.GetValue(),chain.lumi.GetValue())
    
    i = 0
    runList = []
    map_run_LS = {}
    currentlumi = 0
    currentrun = 0
    for event in chain:
        if not args.verbose:  # if using option -v, the carriage return is spoiled by the other print statement in this loop
            if i%1000 == 0:
                sys.stdout.write('Reading chain: {frac:.2f}%   \r'.format(frac=100.*float(i)/entries))
                sys.stdout.flush()
        run = event.run
        lumi = event.luminosityBlock
        if run != currentrun:
            currentrun = run
            currentlumi = 0  # when new run starts, reset current lumi, to avoid not updating it correctly
        if lumi != currentlumi: 
            currentlumi = lumi
            if args.verbose:
                print("%d) run = %d,    LS = %d " % (i,run,lumi))
        i += 1
        if run not in map_run_LS: 
            runList.append(run)
            map_run_LS[run] = []
            map_run_LS[run].append(lumi)
        else:
            if lumi not in map_run_LS[run]: map_run_LS[run].append(lumi)
        if i == args.maxentries: break


    print("")
    leaves = 0
    chain = 0
    

    ##########################
    # start preparing the json
    ##########################

    outjson = open(outfileFullName,'w')
    outjson.write('{\n')

    # first sort the run number (should not be needed, but let's be sure)
    runList = sorted(runList)
    
    printLine("-",20)
    print("There are %d runs" % len(runList))
    if args.verbose:
        print("Printing dictionary with runs and lumisections")
        printLine("-",20)
        for run in runList:
            print("{run} : {lumis}".format(run=run, lumis=sorted(map_run_LS[run])))
    printLine("-",20)
    print("")

    # now, for each run, we have to make LS blocks: 1,2,4,5,6,8,11,12,13... --> [1,2],[4,6],[8,8],[11,...]
    for run in runList:
        lumiList = sorted(map_run_LS[run])

        # already printed before
        #if args.verbose:
        #    print("{run} : [{lumis}]".format(run=run, lumis=",".join(str(x) for x in lumiList)) 

        istart = 0
        lumiblocks = []
        endOfList = False
        if args.verbose:
            printLine("-", 20)
            print(">>> Run {run}".format(run=run))

        while not endOfList:
            lumi_low = lumiList[istart]

            # range() function skips last value: use range up to len - 1 because we need to evaluate  j+1, which must be a valid index (<= len-1)
            if istart == (len(lumiList) - 1):
                lumi_high = lumiList[-1]
                endOfList = True
            else:
                for j in range(istart, len(lumiList)-1):
                    if int(lumiList[j+1] - lumiList[j]) > 1:
                        lumi_high = lumiList[j]
                        istart = j+1
                        break
                    else:
                        # in this case we might have reached the end of list 
                        # for example, we are evaluating the second last number and the last one would belong to the same block
                        # or generally the list ends with a bunch of consecutive LS
                        if int(j+1) == int(len(lumiList)-1):
                            lumi_high = lumiList[-1]
                            endOfList = True

            #print("low,high = %d,%d" % (lumi_low, lumi_high)

            lumiblocks.append("[{low},{high}]".format(low=lumi_low,high=lumi_high))
            if args.verbose:
                print("Appending block --> [{low},{high}]".format(low=lumi_low,high=lumi_high))
                  
        line = '"{run}": [ {blocks} ]{comma}\n'.format(run=run,blocks=", ".join(str(x) for x in lumiblocks),comma="," if run != runList[-1] else "")
        if args.verbose:
            print("Writing in file --> {line}".format(line=line))

        outjson.write(line)            
        lumiBlocks = []  # reset array to save space, not sure it does what I want, but I got killed once before the end of the script
        lumiList = []    # reset array to save space, not sure it does what I want, but I got killed once before the end of the script



    outjson.write('}\n')
    outjson.write('\n')
    outjson.close()
    ###########################
    # Now save things
    ###########################
    print("")
    print("Created file %s" % outfileFullName)
    print("")
                               
         
