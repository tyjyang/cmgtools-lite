import ROOT,os,re
ROOT.gROOT.SetBatch(True)

## USAGE:
## python photos.py -c mu --indir /eos/cms/store/user/arapyan/mc/wp_munu_pythia8/GEN wp_pythia8.root
## python photos.py -c mu --indir /eos/cms/store/user/arapyan/mc/wp_munu_photos/GEN wp_photos.root

if __name__ == "__main__":
    from optparse import OptionParser
    parser = OptionParser(usage="%prog [options] outputfilename ")
    parser.add_option("-c", "--channel"   ,  dest="channel", type='string', default='mu',  help="run tnp ntuples for muons/electrons. default mu")
    parser.add_option("-i", "--indir"     ,  dest="indir"  , type='string', default=''  ,  help="directory with the trees")
    parser.add_option("-m", "--maxentries",  type='int'    , default=-1  ,  help="run only these many entries.")
    (options, args) = parser.parse_args()

    if len(args)<1:
        print "Need at least the name of the output ROOT file. Exiting."
        exit(0)

    ## use xrootd 
    xrdindir  = options.indir
    if '/eos/cms/store/' in xrdindir and not 'eoscms' in xrdindir:
        xrdindir = 'root://eoscms.cern.ch/'+xrdindir
    if '/eos/user/' in xrdindir and not 'eosuser' in xrdindir:
        xrdindir = 'root://eosuser.cern.ch/'+xrdindir

    ## skip strangely bad files
    if 'mu' in options.channel:
        skipfiles = ['wp_munu_pythia8/GEN/GEN_1520000','wp_munu_pythia8/GEN/GEN_5380000','wp_munu_pythia8/GEN/GEN_5620000','wp_munu_pythia8/GEN/GEN_640000','wp_munu_pythia8/GEN/GEN_840000',
                     'wm_munu_pythia8/GEN/GEN_1060000','wm_munu_pythia8/GEN/GEN_3460000','wm_munu_pythia8/GEN/GEN_4660000','wm_munu_pythia8/GEN/GEN_5420000',
                     'wp_munu_photos/GEN/GEN_3380000','wp_munu_photos/GEN/GEN_4100000','wp_munu_photos/GEN/GEN_4600000','wp_munu_photos/GEN/GEN_6980000',
                     'wm_munu_photos/GEN/GEN_4960000']
    else: 
        skipfiles = ['wp_enu_photos/GEN/GEN_1560000','wp_enu_photos/GEN/GEN_2120000','wp_enu_photos/GEN/GEN_500000',
                     'wm_enu_pythia8/GEN/GEN_320000',
                     'wm_enu_photos/GEN/GEN_3120000','wm_enu_photos/GEN/GEN_3320000','wm_enu_photos/GEN/GEN_6840000']

    # check the goodness of files
    treefiles = []
    for rf in os.listdir(options.indir):
        treefile = xrdindir+'/'+rf
        tf = ROOT.TFile(treefile)
        good = True
        if not tf: good=False
        else:
            if tf.IsZombie() or tf.TestBit(ROOT.TFile.kRecovered): good=False
        if any(x in treefile for x in skipfiles): good=False
        if good:
            treefiles.append(treefile)
        else:
            print 'file ',treefile,' is BAD.'

#    treefiles = ['/eos/cms/store/user/arapyan/mc/wm_enu_pythia8/GEN/GEN_6020000.root']

    chain = ROOT.TChain('Events')
    for fname in treefiles:
        print 'adding file {file} to the chain '.format(file=fname)
        chain.Add(fname)

    print 'compiling the ROOT macro'
    ROOT.gROOT.ProcessLine(".L GenEventClass.C+")
    worker = ROOT.GenEventClass(chain);
    worker.setOutfile(args[0])
    if 'mu' in options.channel:
        worker.setFlavor(13)
    else: 
        worker.setFlavor(11)

    worker.Loop(options.maxentries)
