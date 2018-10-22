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

if __name__ == '__main__':

    parser = optparse.OptionParser(usage='usage: %prog [opts] ', version='%prog 1.0')
    parser.add_option('-d', '--directory', type=str, default='', help='directory with the output directories')
    parser.add_option('-m', '--maxsize'  , type=int, default=6 , help='maximum size of the output parts. default 6 gb')
    parser.add_option('-o', '--targetdir', type=str, default='', help='target directory for the output')
    (options, args) = parser.parse_args()

    if not options.targetdir:
        print 'no target directory given... exiting'
        sys.exit(0)
    if not options.directory:
        print 'no source directory given... exiting'
        sys.exit(0)

    cmsenv = os.environ['CMSSW_BASE']

    datasets = set()
    dss = {}

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
            tmp_f = open(options.directory+'/'+sd+'/treeProducerWMass/tree.root.url','r')
            tmp_root = tmp_f.readlines()[0].replace('\n','')
            dss[dsname]['files'  ] .append(tmp_root)
            dss[dsname]['chunks' ] .append(options.directory+'/'+sd)
            if dss[dsname]['avgsize'][1] < 100: ## only calculate avg on the first 100
                newavgsize = getAverageFileSize('/'.join(tmp_root.split('/')[3:]), dss[dsname]['avgsize'])
                dss[dsname]['avgsize'] = newavgsize
            tmp_f.close()
        except:
            continue

    date = datetime.date.today().isoformat()
    condorSubmitCommands = []
    for ds in dss.keys():
        n_chunksPerPart = int(options.maxsize/dss[ds]['avgsize'][0])
        print 'for dataset {d} i will merge {n} files per chunk for roughly {nc} chunks'.format(d=ds,n=n_chunksPerPart,nc=len(dss[ds]['files'])/n_chunksPerPart+1)
        tmp_condor_filename = 'condor_merge_{ds}_{d}.condor'.format(ds=ds,d=date)
        tmp_condor = open(tmp_condor_filename,'w')
        tmp_condor.write('''Executable = mergeScript.sh
use_x509userproxy = $ENV(X509_USER_PROXY)
Log        = merge_{ds}.log
Output     = merge_{ds}.out
Error      = merge_{ds}.error
getenv      = True

transfer_input_files  = fullMergeTrees.py

environment = "LS_SUBCWD={here}"
request_memory = 4000
+MaxRuntime = 7200\n\n'''.format(ds=ds, here=os.environ['PWD']))
        n_part = 1
        while len(dss[ds]['chunks']):
            chunks = ','.join(dss[ds]['chunks'][:n_chunksPerPart])
            tmp_condor.write('arguments = {cmssw} {n} {chunks} {td}\n'.format(cmssw=cmsenv,n=n_part,chunks=chunks,td=options.targetdir))
            tmp_condor.write('queue 1\n\n')
            ## now remove the chunks from the list that are in this job:
            dss[ds]['chunks'] = dss[ds]['chunks'][n_chunksPerPart:]
            n_part += 1
        tmp_condor.close()
        condorSubmitCommands.append('condor_submit {cf}'.format(cf=tmp_condor_filename))

    print 'submitting to condor...'
    for sc in condorSubmitCommands:
        print sc
        os.system(sc)
    print '...done'

        

