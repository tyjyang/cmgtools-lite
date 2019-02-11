import ROOT, os
ROOT.gROOT.SetBatch(True)

## USAGE:
## python leptonTnP.py -i /eos/user/m/mdunser/w-helicity-13TeV/trees/TREES_recoLeptons_and_triggerMatch_latest/ -c mu -o indir

## option -o indir : saves the output trees in the input directory. if none given it saves it where it's run
## option -c       : has to be mu or el
## option -a       : either "trigger" or "selection"


if __name__ == "__main__":
    from optparse import OptionParser
    parser = OptionParser(usage="%prog [options] mc.txt cuts.txt treeDir outputDirSkims ")
    parser.add_option("-c", "--channel"   ,  dest="channel", type='string', default='mu',  help="run tnp ntuples for muons/electrons. default mu")
    parser.add_option("-i", "--indir"     ,  dest="indir"  , type='string', default=''  ,  help="directory with the trees")
    parser.add_option("-m", "--maxentries",  type='int'    , default=-1  ,  help="run only these many entries.")
    parser.add_option("-o", "--outdir"    ,  type='string'    , default='./'  ,  help="save the files in this outdir")
    parser.add_option("-a", "--analyzer"  ,  dest="analyzer", type='string', default='trigger',  help="run tnp ntuples for trigger/selection. default trigger")
    (options, args) = parser.parse_args() 


    if options.outdir == 'indir':
        options.outdir = options.indir

    xrdindir  = options.indir
    xrdoutdir = options.outdir
    ## use xrootd
    if '/eos/cms/store/' in xrdindir and not 'eoscms' in xrdindir:
        xrdindir = 'root://eoscms.cern.ch/'+xrdindir
    if '/eos/user/' in xrdindir and not 'eosuser' in xrdindir:
        xrdindir = 'root://eosuser.cern.ch/'+xrdindir

    ## use xrootd also for output directory and file
    if '/eos/cms/store/' in xrdoutdir and not 'eoscms' in xrdoutdir:
        xrdoutdir = 'root://eoscms.cern.ch/'+xrdoutdir
    if '/eos/user/' in xrdoutdir and not 'eosuser' in xrdoutdir:
        xrdoutdir = 'root://eosuser.cern.ch/'+xrdoutdir

    for sd in os.listdir(options.indir):
        treefile = '/'.join([options.indir,sd,'treeProducerWMass','tree.root'])

        if not os.path.isfile(treefile): continue
        print treefile

        if (options.channel == 'el' and 'SingleMu' in sd): continue
        if (options.channel == 'mu' and 'SingleEl' in sd): continue
        ##if options.analyzer == 'selection' and 'DYJets' in sd: continue
        if 'DYJets' in sd: continue

        print 'runing on file', treefile

        tmp_file = ROOT.TFile(treefile, 'read')
        tmp_tree = tmp_file.Get('tree')
        
        ## make the instance of the worker
        ##if 'el' in options.channel or 'mu' in options.channel:
        ## friends now also for muons!
        ffile = xrdindir+'/friendsNEWSCALE/tree_Friend_'+sd+'.root'
        tmp_friend_file = ROOT.TFile(ffile)
        tmp_friend_tree = tmp_friend_file.Get('Friends')
        ##else: 
        ##    tmp_friend_tree = None
    
        if options.analyzer=='trigger'     : 
            print 'compiling trigger ntuple producer'
            ROOT.gROOT.ProcessLine(".L TnPNtuplesTriggerEfficiency.C+" )
            tnp_worker = ROOT.TnPNtuplesTriggerEfficiency(tmp_tree,tmp_friend_tree)
        elif options.analyzer=='selection' : 
            print 'compiling selection ntuple producer'
            ROOT.gROOT.ProcessLine(".L TnPNtuplesSelectionEfficiency.C+" )
            tnp_worker = ROOT.TnPNtuplesSelectionEfficiency(tmp_tree,tmp_friend_tree)
        else: 
            print "analyzer can be either trigger/selection."
            exit(1)

        if 'mu' in options.channel:
            tnp_worker.setFlavor(13)
        else: 
            tnp_worker.setFlavor(11)

        tnp_worker.setOutfile(xrdoutdir+'/'+options.analyzer+'TnP_'+os.path.basename(sd)+'_{ch}.root'.format(ch=options.channel))
        tnp_worker.Loop(options.maxentries)
    
        del tnp_worker ## this segfaults otherwise

