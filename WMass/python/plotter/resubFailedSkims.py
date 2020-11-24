import array
import json
import ROOT
import re 
import os

##
## script to check failed jobs when running skims on CMGTools ntuples
##
# WARNING: the script currently checks for the presence of the folder, it doesn't check whether the root file is in a good state. This feature can be added in the future
##
##


if __name__ == "__main__":
    from optparse import OptionParser
    parser = OptionParser(usage='%prog [options] ')
    parser.add_option('-i','--indir',     dest='indir',     default='',   type='string', help='Folder with the condor file (the one with all the logs)')
    #parser.add_option('-o','--outdir',     dest='outdir',     default='SAME',   type='string', help='outdput directory to write the new condor files and logs (can be the same as input). If SAME, use same as input folder')
    parser.add_option('--ntuples-dir-in',     dest='ntuplesDirIn',     default='',   type='string', help='Path to input unskimmed ntuples')
    parser.add_option('--ntuples-dir-out',     dest='ntuplesDirOut',     default='',   type='string', help='Path to output skimmed ntuples')
    parser.add_option('-d', '--dry-run', dest='dryRun' , default=False , action='store_true',   help='Print commands, but do not run them')
    parser.add_option('-f', '--friends', dest='friends' , default=False , action='store_true',   help='Specify if doing friend trees, things need to be handled differently (paths in --ntuples-dir-in and --ntuples-dir-out do no need to end with "/friends/"')
    (options, args) = parser.parse_args()

    
    if not options.indir:
        print "WARNING: you must specify an input folder with option -i. Abort"
        quit()

    outdir = os.path.dirname(options.indir)
    if "_resub" in outdir:
        outdirNoNum = outdir.split("_resub")[-1]
        outdir.replace("resub"+outdirNoNum, "resub"+str(int(outdirNoNum)+1))
        outdir += "/"
    else:
        outdir += "_resub1/"
    if not os.path.isdir(outdir):
            os.system('mkdir -p {od}'.format(od=outdir))

    if not options.ntuplesDirIn:
        print "WARNING: you must specify the path to input ntuples with option --ntuples-dir-in. Abort"
        quit()

    if not options.ntuplesDirOut:
        print "WARNING: you must specify the path to skimmed ntuples with option --ntuples-dir-out. Abort"
        quit()

    if options.friends:        
        if "/friends" not in options.ntuplesDirIn:  options.ntuplesDirIn  += "/friends/"
        if "/friends" not in options.ntuplesDirOut: options.ntuplesDirOut += "/friends/"

    ##### check missing pieces in the skim
    input_samples = [f for f in os.listdir(options.ntuplesDirIn) if not f.endswith(".txt") and not "friends" in f]
    output_samples = [f for f in os.listdir(options.ntuplesDirOut) if not f.endswith(".txt") and not "friends" in f]
    missing_samples = []
    for s in input_samples:
        if s not in output_samples:
            if options.friends:
                missing_samples.append(s.replace("tree_Friend_","").replace(".root",""))
            else:
                missing_samples.append(s)
    if len(missing_samples):
        print"="*30
        print "There are %d missing samples%s:" % (len(missing_samples),
                                                   " for friend trees" if options.friends else "" )    
        print"-"*30
        for i,x in enumerate(sorted(missing_samples)): print "%d) %s" % (i+1,x)
        print"="*30
    #####

    condor_files = [os.path.join(options.indir, f) for f in os.listdir(options.indir) if (f.endswith(".condor") and os.path.isfile(os.path.join(options.indir, f)) and "resub" not in f)]
    #print condor_files
    file_has_missing_samples = {}
    for thisf in condor_files:
        #print "-"*30
        #print thisf
        outf = outdir + os.path.basename(thisf)
        #print outf
        file_has_missing_samples[outf] = False # initialize with False
        of = open(outf,"w")
        f  = open(thisf,"r")
        for line in f:
            if any(line.startswith(x) for x in ["arguments", "queue", "Error", "  ", "\n"]):
                if line.startswith("arguments"):
                    # check if line has to be copied, then also add queue line
                    component = line.split("--component")[1].lstrip().split()[0]
                    if component in missing_samples:
                        file_has_missing_samples[outf] = True
                        of.write("\n"+line)
                        if options.friends:
                            of.write("Error = " + outdir + component + ".error\n")
                        of.write("queue 1\n")
            else:
                #print line.rstrip().lstrip()
                of.write(line)
        f.close()
        of.close()    

    print "-"*30
    print ""

    if any(file_has_missing_samples[x] for x in file_has_missing_samples.keys()):
        print "There are some jobs to resubmit :-("
        for key in file_has_missing_samples:
            if file_has_missing_samples[key]:
                cmd = "condor_submit %s" % key
                print cmd
                if not options.dryRun: 
                    os.system(cmd)
    else:
        print "There are no jobs to resubmit :-)"

    print "-"*30
    print ""

        
