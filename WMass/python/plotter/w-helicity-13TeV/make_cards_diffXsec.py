import os
from datetime import datetime

# python w-helicity-13TeV/make_cards_diffXsec.py -C el -q cmscaf1nd --groupSignalBy 10 --syst -d

from optparse import OptionParser
parser = OptionParser(usage="%prog [options]")
parser.add_option("-d", "--dry-run",    dest="dryRun",   action="store_true", default=False, help="Do not run the job, only print the command");
parser.add_option("-p", "--print-only", dest="printOnly",   action="store_true", default=False, help="Just print commands, do not execute anything");
parser.add_option("-f", "--force",      dest="force",   action="store_true", default=False, help="Force running without question below (useful only when using PDF systematics)");
parser.add_option("-o", "--outdir", dest="outdir", type="string", default='', help="Name of output folder. It will ignore option -s and the automatic date in the name. This option is particulaly useful to produce output in an already existing folder (e.g. when you do signal and background in different days, which would create a new folder). Warning: with condor it might overwrite things!!!");
parser.add_option("-C", "--channel", dest="channel", type="string", default='el', help="Channel. either 'el' or 'mu'");
parser.add_option("-s", "--suffix", dest="suffix", type="string", default=None, help="Append a suffix to the default outputdir (diffXsec_<date>)");
parser.add_option("-q", "--queue", dest="queue", type="string", default="cmscaf1nd", help="Select the queue to use");
parser.add_option("-r", "--run", dest="run", type="string", default="sb", help="Which components to run: s for signal, b for backgrounds or sb for both");
parser.add_option("--syst", dest="addSyst", action="store_true", default=False, help="Add PDF and QCD scale systematics to the signal (need incl_sig directive in the MCA file)");
parser.add_option("--skip-syst-jobs", dest="skipSystJobs", action="store_true", default=False, help="Skip running jobs for syst (useful to retrieve only nominal one for tau and Z and run locally). Cannot just skip option --syst, because otherwise tau and Z are made with other backgrounds");
#### options for differential xsec
parser.add_option(      "--xsec-sigcard-binned", dest="xsec_sigcard_binned",   action="store_true", default=False, help="When doing differential cross-section, will make 1 signal card for each 2D template bin (default is False because the number of cards easily gets huge)");
parser.add_option(      "--groupSignalBy",       dest="groupSignalBy", type="int", default='0', help="Group signal bins in bunches of N (pass N as argument). Default is 0, meaning not using this option. This option will reduce the number of chunk datacards for signal,but jobs will last for longer");
parser.add_option("-l", "--leptonDef", dest="leptonDef", type="string", default='dressed', help="Lepton definition 'preFSR' or 'dressed'");
parser.add_option("--useLSF", action='store_true', default=False, help="force use LSF. default is using condor");
parser.add_option("--usePickle", dest="usePickle", action="store_true", default=False, help="Read Sum Weights from Pickle file (needed only if using old samples that did not have the histogram inside). By default, the histogram is used. Check in mcAnalysis whether it is really implemented");
(options, args) = parser.parse_args()


if not options.xsec_sigcard_binned and not options.groupSignalBy:
    print ""
    print "You haven't selected any of options --groupSignalBy and --xsec-sigcard-binned."
    print "This implies you are trying to make a single signal template without gen eta/pt categories"
    print "Maybe this was not your intention, are you sure you want to proceed? [y/N]\n"
    if raw_input()!='y':
        print 'Aborting'
        exit()

if options.xsec_sigcard_binned and options.addSyst and not options.force:
    print ""
    print "You are trying to run the differential cross-section measurement making a signal template/card for each bin of the 2D templates."
    print "In addition, you are trying to add PDF and QCD scale systematics (60 variations for PDF)"
    print "This is a huge number of outputs, are you sure you want to proceed? [y/N]\n"
    if raw_input()!='y':
        print 'Aborting'
        exit()
    

channel = options.channel
if channel not in ["mu", "el"]:
    print "Error in make_cards_diffXsec.py: unknown lepton flavour, use -c el|mu. Exit"
    exit(0)

if options.leptonDef == "preFSR":
    genLepVar = "GenLepPreFSR"
    genwVar = "prefsrw"
elif options.leptonDef == "dressed":
    genLepVar = "GenLepDressed"
    genwVar = "genw"
else:
    print "Error in make_cards_diffXsec.py: unknown lepton definition, use -l preFSR|dressed. Exit"
    exit(0)

PROG="w-helicity-13TeV/make_diff_xsec_cards.py"
outdirbase = "diffXsec"
QUEUE=str(options.queue)

if channel == "el":

    BASECONFIG="w-helicity-13TeV/wmass_e"
    MCA=BASECONFIG+'/mca-80X-wenu-xsec.txt'
    # binning for the reco template, eta on x axis, pt on y axis. 
    #BINNING="\"[-2.5,-2.3,-2.1,-1.9,-1.7,-1.566,-1.4442,-1.2,-1.0,-0.8,-0.6,-0.4,0.2,0,0.2,0.4,0.6,0.8,1.0,1.2,1.4442,1.566,1.7,1.9,2.1,2.3,2.5]*[30,33,36,39,42,45]\""
    #BINNING="\"[-2.5,-2.3,-2.1,-1.9,-1.7,-1.566,-1.4442,-1.2,-1.0,-0.8,-0.6,-0.4,-0.2,0,0.2,0.4,0.6,0.8,1.0,1.2,1.4442,1.566,1.7,1.9,2.1,2.3,2.5]*[30,33,36,39,42,45]\""
    #BINNING="\"[-2.5,-1.566,-1.4442,0,1.4442,1.566,2.5]*[30,35,40,45]\""
    #BINNING="\"[-2.5,-2.3,-2.1,-1.9,-1.7,-1.566,-1.4442,-1.3,-1.2,-1.1,-1.0,-0.9,-0.8,-0.7,-0.6,-0.5,-0.4,-0.3,-0.2,-0.1,0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4442,1.566,1.7,1.9,2.1,2.3,2.5]*[30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45]\""
    #BINNING="\"[-2.5,-2.3,-2.1,-1.9,-1.7,-1.566,-1.4442,-1.2,-1.0,-0.8,-0.6,-0.4,-0.2,0,0.2,0.4,0.6,0.8,1.0,1.2,1.4442,1.566,1.7,1.9,2.1,2.3,2.5]*[30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45]\""
    # gen binning is used to defined the signal categories (it is the equivalent of Yw for the rapidity/helicity
    # we use |eta| here for the gen
    #GENBINNING="\"[0,0.2,0.4,0.6,0.8,1.0,1.2,1.4442,1.566,1.7,1.9,2.1,2.3,2.5]*[30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45]\""
    #etabinning='[-2.5,-2.4,-2.3,-2.15,-2.0,-1.85,-1.7,-1.566,-1.4442,-1.3,-1.15,-1.0,-0.9,-0.8,-0.7,-0.6,-0.5,-0.4,-0.3,-0.2,-0.1,0.,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.15,1.3,1.45,1.6,1.75,1.9,2.05,2.2,2.4,2.5]'

    #GENBINNING="\"[0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.15,1.3,1.45,1.6,1.75,1.9,2.05,2.2,2.35,2.5]*[26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45.46.47.48.49.50]\""

    # eta for gen and reco can be different. The gen should be consistent with muons

    #binsEta = list(i*0.1 for i in range(1,14)) + list(x for x in [1.4442,1.566,1.65]) + list(1.65+0.15*i for i in range(1,6))  # 0.1 from 0 to 1.3; 0.15 from 1.65 to 2.4, then 2.4 and 2.5, the gap region is 1.3, 1.4442, 1.566, 1.65
    #binsEta = list(i*0.1 for i in range(1,25))  # 0.1, also around the gap, from 0 to 2.4 (not 2.5, bad data/MC)
    binsEta = list(i*0.1 for i in range(1,14)) + list((1.3 + i*0.2) for i in range(1,5)) + [2.4]
    negstring = '[' + ','.join(reversed(  list(str(-1.*i) for i in binsEta) ))
    posstring = ','.join( list(str(i) for i in binsEta) ) + ']'
    etabinning= negstring+',0.,'+posstring

    ## variable binning in pt
    #ptbinning = '['+','.join(str(i) for i in range(26,56))+']'  # 1 from 26 to 55
    #ptbinning = '['+','.join(str(i) for i in range(30,44,2))+',45]'  # 1 from 26 to 45
    ptbinning = "[26,28,30,31.5,33,34.5,36,37.5,39,40.5,42,43.5,45,46.5,48,50,52,54,56]"  #47.5,50,52.5,55]"
    #ptbinning = "[26,28,30,33,36,39,42,45,48,51,54,56]"  #47.5,50,52.5,55]"
    #ptbinning = "[26,28,30,32,34,36,38,40,42,44,46,48,50,52,54,56]"  #47.5,50,52.5,55]"

    GENBINNING = "'[0.,"+posstring+"*"+ptbinning+"'"
    #recoetabinning = "[-2.5,-2.4,-2.3,-2.2,-2.1,-2.0,-1.9,-1.8,-1.7,-1.6,-1.566,-1.5,-1.4442,-1.4,-1.3,-1.2,-1.1,-1.0,-0.9,-0.8,-0.7,-0.6,-0.5,-0.4,-0.3,-0.2,-0.1,0.,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.4442,1.5,1.566,1.6,1.7,1.8,1.9,2.0,2.1,2.2,2.3,2.4,2.5]"
    #BINNING    = "'"+recoetabinning+"*[30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45]'"
    #recoetabinning = "[-2.4,-2.2,-2.0,-1.8,-1.6,-1.566,-1.4442,-1.4,-1.2,-1.1,-1.0,-0.9,-0.8,-0.7,-0.6,-0.5,-0.4,-0.3,-0.2,-0.1,0.,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.4,1.4442,1.566,1.6,1.8,2.0,2.2,2.4]"
    #recoetabinning = "[-2.4,-2.3,-2.2,-2.1,-2.0,-1.9,-1.8,-1.7,-1.6,-1.566,-1.5,-1.4442,-1.4,-1.3,-1.2,-1.1,-1.0,-0.9,-0.8,-0.7,-0.6,-0.5,-0.4,-0.3,-0.2,-0.1,0.,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.4442,1.5,1.566,1.6,1.7,1.8,1.9,2.0,2.1,2.2,2.3,2.4]"
    recoetabinning = "[-2.4,-2.1,-1.9,-1.7,-1.566,-1.5,-1.4442,-1.3,-1.2,-1.1,-1.0,-0.9,-0.8,-0.7,-0.6,-0.5,-0.4,-0.3,-0.2,-0.1,0.,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4442,1.5,1.566,1.7,1.9,2.1,2.4]"
    #    BINNING    = "'"+recoetabinning+"*[30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45]'"
    #recoptbinning = '['+','.join(str(i) for i in range(30,44,2))+',45]'  # 1 from 26 to 45
    recoptbinning = "[30,31.5,33,34.5,36,37.5,39,40.5,42,43.5,45,46.5,48,50,52,54,56]"  #47.5,50,52.5,55]"
    #recoptbinning = "[30,33,36,39,42,45,48,51,54,56]"  #47.5,50,52.5,55]"
    #recoptbinning = "[30,32,34,36,38,40,42,44,46,48,50,52,54,56]"
    #recoptbinning = "[30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56]" #,47.5,50,52.5,55]"
    BINNING    = "'"+recoetabinning+"*" + recoptbinning + "'"

    # for quick tests
    #BINNING    = "'[-2.4,-2.1,-1.6,-1.4,-1.0,0,1.0,1.4,1.6,2.1,2.4]*[26,30,35,40,45]'"
    #GENBINNING = "'[0,1.0,1.4,1.6,2.1,2.4]*[26,30,35,40,45]'"
    #GENBINNING = "'[0,0.6,1.2,1.8,2.4]*[30,35,40,45]'"


    CUTFILE=BASECONFIG+'/wenu_80X_xsec.txt'
    SYSTFILE=BASECONFIG+'/systsEnv.txt'
    #TREEPATH="/eos/cms/store/group/dpg_ecal/comm_ecal/localreco/TREES_1LEP_80X_V3_WENUSKIM_V5_TINY"
    #TREEPATH="/eos/cms/store/cmst3/group/wmass/w-helicity-13TeV/trees/TREES_electrons_1l_V6_TINY"
    #TREEPATH="/eos/cms/store/cmst3/group/wmass/mciprian/TREES_1LEP_80X_V3_WSKIM_NEW"
    TREEPATH="/afs/cern.ch/work/m/mciprian/w_mass_analysis/heppy/CMSSW_8_0_25/src/CMGTools/WMass/python/plotter/TREE_4_XSEC_AFS"
    VAR="\"ptElFull(LepGood1_calPt,LepGood1_eta):LepGood1_eta\""
    #WEIGHTSTRING=" \'puw2016_nTrueInt_36fb(nTrueInt)*trgSF_We(LepGood1_pdgId,LepGood1_pt,LepGood1_eta,2)*leptonSF_We(LepGood1_pdgId,LepGood1_pt,LepGood1_eta)\' "
    #WEIGHTSTRING=" \'puw2016_nTrueInt_36fb(nTrueInt)*lepSF(LepGood1_pdgId,LepGood1_pt,LepGood1_eta,LepGood1_SF1,LepGood1_SF2,LepGood1_SF3)\' "
    WEIGHTSTRING=" \'puw2016_nTrueInt_36fb(nTrueInt)*_get_electronSF_TriggerAndID(LepGood1_pdgId,LepGood1_calPt,LepGood1_eta)*LepGood1_SF2 * eleSF_L1Eff(LepGood1_pt,LepGood1_eta)\' "
    LUMI=35.9
    binnedSignalMCA = "mca-80X-wenu-sigInclCharge_binned_eta_pt.txt"

else:

    BASECONFIG   = 'w-helicity-13TeV/wmass_mu'
    MCA          = BASECONFIG+'/mca-wmu-xsec.txt'
    CUTFILE      = BASECONFIG+'/cuts_wmu_xsec.txt'
    SYSTFILE     = BASECONFIG+'/systsEnv.txt'
    #TREEPATH     = '/eos/user/m/mdunser/w-helicity-13TeV/trees/TREES_latest_new_1muskim/'
    TREEPATH     = '/afs/cern.ch/work/m/mdunser/public/wmassTrees/SKIMS_muons_latest/'
    VAR          = '\'ptMuFull(LepGood1_calPt,LepGood1_eta):LepGood1_eta\''

    ## variable binning in eta
    #binsEta = list(i*0.1 for i in range(1,21)) + list((2.0 + i*0.2) for i in range(1,3))
    #binsEta = list(i*0.1 for i in range(1,13)) + list((1.2 + i*0.2) for i in range(1,7))
    #binsEta = list(i*0.2 for i in range(1,13)) 
    binsEta = list(i*0.1 for i in range(1,14)) + list((1.3 + i*0.2) for i in range(1,5)) + [2.4]
    #binsEta = list(i*0.1 for i in range(1,25))
    negstring = '['+','.join(reversed(  list(str(-1.*i) for i in binsEta) ))
    posstring = ','.join( list(str(i) for i in binsEta) )+']'
    etabinning= negstring+',0.,'+posstring
    #etabinning= '[0.,'+posstring
    #etabinning= negstring+', 0.0]'

    ## variable binning in pt
    #ptbinning = '['+','.join(str(i) for i in range(26,48,2))+']'
    #
    #ptbinning = "[26,28,30,31.5,33,34.5,36,37.5,39,40.5,42,43.5,45,46.5,48,50,52,54,56]" #,47.5,50,52.5,55]"
    ptbinning = "[26,28,30,31.5,33,34.5,36,37.5,39,40.5,42,43.5,45]" #,47.5,50,52.5,55]"
    #
    #ptbinning = "[26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56]" #,47.5,50,52.5,55]"

    # different gen binning
    binsEta = list(i*0.1 for i in range(1,14)) + list((1.3 + i*0.2) for i in range(1,5)) + [2.4]
    #posstring = ','.join( list(str(i) for i in binsEta) )+']'
    posstring = '2.4]'

    BINNING      = '\''+etabinning+'*'+ptbinning+'\''
    #ptbinning = "[26,28,30,32,34,36,38,40,42,44,46,48,50,52,54,56]" #,47.5,50,52.5,55]"
    #ptbinning = "[26,28,30,31.5,33,34.5,36,37.5,39,40.5,42,43.5,45,46.5,48,50,52,54,56]"
    ptbinning = "[26,45]"
    GENBINNING = "'[0.,"+posstring+"*"+ptbinning+"'"
    #WEIGHTSTRING=  ' \'puw2016_nTrueInt_36fb(nTrueInt)*LepGood_SF1[0]*LepGood_SF2[0]\' '  # ok, SF4 is for e only, and SF3 is basically 1
    #WEIGHTSTRING=  ' \'puw2016_nTrueInt_36fb(nTrueInt)*_get_muonSF_recoToSelection(LepGood1_pdgId,LepGood1_calPt,LepGood1_eta)*_get_muonSF_selectionToTrigger(LepGood1_pdgId,LepGood1_calPt,LepGood1_eta,LepGood1_charge)*prefireJetsWeight(LepGood1_eta)\' '  # this has charge dependent non-eta-smoothed SF for trigger, eta-smoothed for reco2sel
    #
    WEIGHTSTRING=  ' \'puw2016_nTrueInt_36fb(nTrueInt)*_get_muonSF_selectionToTrigger(LepGood1_pdgId,LepGood1_calPt,LepGood1_eta,LepGood1_charge)*_get_muonSF_recoToSelection(LepGood1_pdgId,LepGood1_calPt,LepGood1_eta)*prefireJetsWeight(LepGood1_eta)\' '  # this has charge dependent non-eta-smoothed SF for trigger, eta-smoothed for reco2sel
    # recoToSelection could be LepGood_SF2[0]
    #WEIGHTSTRING=  ' \'puw2016_nTrueInt_36fb(nTrueInt)*_get_muonSF_selectionToTrigger(LepGood1_pdgId,LepGood1_calPt,LepGood1_eta,LepGood1_charge,0.0,0,1)*_get_muonSF_recoToSelection(LepGood1_pdgId,LepGood1_calPt,LepGood1_eta,1)*prefireJetsWeight(LepGood1_eta)\' '  # this has charge dependent non-eta-smoothed SF for trigger, eta-smoothed for reco2sel
    LUMI = 35.9    
    binnedSignalMCA = "mca-80X-wmunu-sigInclCharge_binned_eta_pt.txt"

    #BINNING    = "'[-2.4,-1.6,-1.0,0,1.0,1.6,2.4]*[26,30,35,40,45]'"
    #GENBINNING = "'[0,1.0,1.6,2.4]*[26,30,35,40,45]'"
    #GENBINNING = "'[0,0.6,1.2,1.8,2.4]*[26,30,35,40,45]'"


print ""
print "================================="
print "GEN BINNING (eta - pt)"
print "---------------------------------"
print GENBINNING
print "================================="
print "RECO BINNING (eta - pt)"
print "---------------------------------"
print BINNING
#print etabinning
#print ptbinning
print "================================="
print ""

OUTDIR="%s_%s_%s" % (outdirbase,channel,datetime.now().strftime("%Y_%m_%d"))
if options.groupSignalBy: OUTDIR += ("_group%s" % str(options.groupSignalBy))
if options.suffix: OUTDIR += ("_%s" % options.suffix)

if options.outdir != "": OUTDIR = options.outdir

components=""
if "s" in options.run:
    components += " -s "
if "b" in options.run:
    components += " -b "

###############
# create the mca for signal bins: it is needed only if you want to group signal bins, because you have to use option -p of mcAnalysis to select them as different processes
#
print "Creating signal MCA file: {base}/mca-includes/{mca}".format(base=BASECONFIG,mca=binnedSignalMCA)
makeMCAcommand="python w-helicity-13TeV/printBinnnedSignalMCA.py -o {base}/mca-includes/ -n {mca}  -b {genbin} -x 'abs({genLepVar}_eta[0])' -y '{genLepVar}_pt[0]' -c {ch} -l {genLepDef}".format(base=BASECONFIG,mca=binnedSignalMCA,genbin=GENBINNING,genLepVar=genLepVar,ch=channel, genLepDef=options.leptonDef)
if not options.printOnly:
    os.system(makeMCAcommand)
else:
    print ""
    print makeMCAcommand
    print ""
print "================================="
print ""

###############

cmd="python " + " ".join([PROG,MCA,CUTFILE,VAR,BINNING,SYSTFILE,OUTDIR]) + \
    (" -W %s " % WEIGHTSTRING) + (" -P %s " % TREEPATH) + (" -q %s " % QUEUE) + components + \
    (" -C %s " % channel) + (" -l %f " % LUMI) + (" --gen-binning %s " % GENBINNING) + (" -w %s " % genwVar)  
if options.xsec_sigcard_binned: cmd += '  --xsec-sigcard-binned '
if options.groupSignalBy: cmd += '  --groupSignalBy %d ' % options.groupSignalBy
if options.dryRun: cmd += '  --dry-run '
if options.addSyst: cmd += '  --pdf-syst --qcd-syst '
if options.useLSF: cmd += '  --useLSF '
if options.skipSystJobs: cmd += '  --skip-syst-jobs '
if options.usePickle: cmd += ' --usePickle ' 
if not options.printOnly:
    os.system(cmd)
else:
    print cmd
    print ""

print ""
