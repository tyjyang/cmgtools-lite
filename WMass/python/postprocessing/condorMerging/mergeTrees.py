import ROOT, commands, os, sys, optparse, datetime


## USAGE:
## ====================================
## python mergeTrees.py -d <inputDirectoryWithStructure> -o <directoryOnEos>

## the input directory must have all the XXXX_ChunkXYZ directories in it
## the output directory on eos does not necessarily have to exist

def getAverageFileSize(p, thing):
    avg = thing[0]
    n   = thing[1]
    size = os.path.getsize(p)/1e9 ## in GB
    newavg = (n*avg + size)/(n+1)
    retval = [newavg, n+1]
    return retval

def checkIntegrity(url):
    """check if file is not corrupted"""
    f=ROOT.TFile.Open(url)
    try:
        if f.IsZombie():    
            raise UserWarning(url,'is probably corrupted')
        if f.TestBit(ROOT.TFile.kRecovered):
            raise UserWarning(url,'is in fishy state, was recovered')
        if f.GetListOfKeys().GetSize()==0:
            raise UserWarning(url,'has no keys...fishy...')
    except:
        raise UserWarning(url,'is unusable at all')
    f.Close()

if __name__ == '__main__':

    parser = optparse.OptionParser(usage='usage: %prog [opts] ', version='%prog 1.0')
    parser.add_option('-d', '--directory', type=str, default='',    help='directory with the output directories [%default]')
    parser.add_option(      '--strict',              default=False, help='perform a strict check on the integrity of the root files [%default]', action='store_true')
    parser.add_option('--skipAllButTreeAndSkimReport', dest='skipAllButTreeAndSkimReport',
                      default=False,action='store_true',
                      help='When hadding stuff, skip hadd of some pck files that might be missing but not necessary (you only need the tree and the report from skimAnalyzerCount for number of events before any selection')
    parser.add_option(      '--dryRun',              default=False, help='dry run (do not submit to condor) [%default]', action='store_true')
    parser.add_option('-m', '--maxsize'  , type=float, default=6. ,    help='maximum size of the output parts. [%default] gb')
    parser.add_option('-o', '--targetdir', type=str, default='',    help='target directory for the output [%default]')
    (options, args) = parser.parse_args()

    if not options.targetdir:
        print 'no target directory given... exiting'
        sys.exit(0)
    if not options.directory:
        print 'no source directory given... exiting'
        sys.exit(0)
    if options.strict:
        print 'Will perform a strict integrity check on every single file'
        print 'Submitted jobs won\'t perform check'
        
    cmsenv = os.environ['CMSSW_BASE']

    datasets = set()
    dss = {}

    badChunksList=[]
    for isd,sd in enumerate(os.listdir(options.directory)):
        if not 'Chunk' in sd: continue
        if not os.path.isdir(options.directory+'/'+sd): continue
        dsname = '_'.join(sd.split('_')[:-1])
        datasets.add(dsname)
        if dsname not in dss.keys():
            dss[dsname] = {}
            dss[dsname]['files']   = []
            dss[dsname]['chunks']  = []
            dss[dsname]['avgsize'] = [0., 0]
        try:
            ## this is basically already a check for processed chunks
            urlfile = options.directory+'/'+sd+'/treeProducerWMass/tree.root.url'
            tmp_f = open(urlfile,'r')
            tmp_root = tmp_f.readlines()[0].replace('\n','')
            if options.strict:
                checkIntegrity(tmp_root)
            dss[dsname]['files'  ] .append(tmp_root)
            dss[dsname]['chunks' ] .append(os.path.abspath(options.directory+'/'+sd))
            if dss[dsname]['avgsize'][1] < 300: ## only calculate avg on the first 300. this is sparta!
                newavgsize = getAverageFileSize('/'.join(tmp_root.split('/')[3:]), dss[dsname]['avgsize'])
                dss[dsname]['avgsize'] = newavgsize
            tmp_f.close()
        except UserWarning as uw:
            badChunksList.append(options.directory+'/'+sd)
            continue
        except Exception:
            badChunksList.append(options.directory+'/'+sd)
            continue

    date = datetime.date.today().isoformat()
    condorSubmitCommands = []
    for ds in dss.keys():
        n_chunksPerPart = int(options.maxsize/dss[ds]['avgsize'][0])
        print 'for dataset {d} i will merge {n} files per chunk for roughly {nc} chunks'.format(d=ds,n=n_chunksPerPart,nc=len(dss[ds]['files'])/n_chunksPerPart+1)
        tmp_condor_filename = 'condor_merge_{ds}_{d}.condor'.format(ds=ds,d=date)
        tmp_condor = open(tmp_condor_filename,'w')
        tmp_condor.write('''Executable = mergeScript.sh
use_x509userproxy = true
Log        = log_merge_{ds}_$(ProcId).log
Output     = log_merge_{ds}_$(ProcId).out
Error      = log_merge_{ds}_$(ProcId).error
getenv      = True

transfer_input_files  = fullMergeTrees.py

environment = "LS_SUBCWD={here}"
request_memory = 2000
+MaxRuntime = 7200\n\n'''.format(ds=ds, here=os.environ['PWD']))
        if os.environ['USER'] in ['mciprian']:
            tmp_condor.write('+AccountingGroup = "group_u_CMS.CAF.ALCA"\n\n')

        n_part = 1
        while len(dss[ds]['chunks']):
            chunks = ','.join(dss[ds]['chunks'][:n_chunksPerPart])
            tmp_condor.write('arguments = {cmssw} {n} {chunks} {td} {strict} {skipSome}\n'.format(cmssw=cmsenv,n=n_part,chunks=chunks,td=options.targetdir,strict=options.strict,skipSome=options.skipAllButTreeAndSkimReport))
            tmp_condor.write('queue 1\n\n')
            ## now remove the chunks from the list that are in this job:
            dss[ds]['chunks'] = dss[ds]['chunks'][n_chunksPerPart:]
            n_part += 1
        tmp_condor.close()
        condorSubmitCommands.append('condor_submit {cf}'.format(cf=tmp_condor_filename))

    if options.dryRun:
        print 'this was just a dry run...'
    else:
        print 'submitting to condor...'
        for sc in condorSubmitCommands:
            print sc
            os.system(sc)

    #list files to resubmit
    print 'Listing files to resubmit below'
    for f in badChunksList:
        print 'cmgResubChunk -q HTCondor -t 4000 %s'%f


    print '...done'

        

