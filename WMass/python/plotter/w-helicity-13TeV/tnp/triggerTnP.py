import ROOT, os

## USAGE:
## python triggerTnP.py -i /eos/user/m/mdunser/w-helicity-13TeV/trees/TREES_2018-07-06-2l_triggerMatch_MUONS/ -c mu -o indir

## option -o indir : saves the output trees in the input directory. if none given it saves it where it's run
## option -c       : has to be mu or el


ROOT.gROOT.ProcessLine(".L TnPNtuplesTriggerEfficiency.C+" )

if __name__ == "__main__":
    from optparse import OptionParser
    parser = OptionParser(usage="%prog [options] mc.txt cuts.txt treeDir outputDirSkims ")
    parser.add_option("-c", "--channel"   ,  dest="channel", type='string', default='mu',  help="run tnp ntuples for muons/electrons. default mu")
    parser.add_option("-i", "--indir"     ,  dest="indir"  , type='string', default=''  ,  help="directory with the trees")
    parser.add_option("-m", "--maxentries",  type='int'    , default=-1  ,  help="run only these many entries.")
    parser.add_option("-o", "--outdir"    ,  type='string'    , default='./'  ,  help="save the files in this outdir")
    (options, args) = parser.parse_args() 


    if options.outdir == 'indir':
        options.outdir = options.indir

    for sd in os.listdir(options.indir):
        treefile = '/'.join([options.indir,sd,'treeProducerWMass','tree.root'])
        if not os.path.isfile(treefile): continue

        print 'runing on file', treefile

        tmp_file = ROOT.TFile(treefile)
        tmp_tree = tmp_file.Get('tree')
        
        ## make the instance of the worker
        if 'mu' in options.channel:
            tnp_worker = ROOT.TnPNtuplesTriggerEfficiency(tmp_tree)
            tnp_worker.setFlavor(13)
        else: 
            ffile = options.indir+'/friends/+'+os.path.basename(treefile) ## FIX THIS
            tnp_worker = ROOT.TnPNtuplesTriggerEfficiency(tmp_tree)
            tnp_worker.setFlavor(11)

        tnp_worker.setOutfile(options.outdir+'/'+'triggerTnP_'+os.path.basename(sd)+'_{ch}.root'.format(ch=options.channel))
        tnp_worker.Loop(options.maxentries)
    
        del tnp_worker ## this segfaults otherwise

