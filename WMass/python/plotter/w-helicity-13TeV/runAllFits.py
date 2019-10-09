# USAGE: python runAllFits.py cards_el el
import os,sys,datetime,random
from submitToys import makeCondorFile,jobstring_tf


from optparse import OptionParser
parser = OptionParser(usage="%prog [options] cardsDir flavor ")
parser.add_option("-q", "--queue",      default=False, dest="queue",  action="store_true",   help="Run jobs on condor instead of locally");
parser.add_option('-r', '--runtime',    default=24, type=int,   help='New runtime for condor resubmission in hours. default: 24h (combined fits may be long)');
parser.add_option(      '--regularize', default=False,  dest="regularize", action="store_true",   help='do regularization with optimized taus');
parser.add_option(      '--useSciPy',   default=False, action='store_true',  help='use the slower SciPy minimizer (should be not necessary any more)');
parser.add_option('-s', '--suffix',     default='', type="string", help='use suffix for the output dir');
parser.add_option(      '--toys',       default=0,  type=int, help='run that number of toys, randomly generated');
parser.add_option('-t', '--threads',    default=None,  type=int, dest='nThreads', help='use nThreads in the fit (suggested 2 for single charge, 1 for combination)')
(options, args) = parser.parse_args()

cardsdir = os.path.abspath(args[0])
channel = args[1]
if channel not in ['mu','el','lep']:
    print "Channel must be either mu or el or lep (el-mu combination). Exiting."
    sys.exit()

pois = [('poim1',''),('poim0',' --POIMode none ')]
#pois = [('poim0',' --POIMode none ')]
expected = [('exp1',' -t -1 '),('exp0',' -t 0 ')]
#BBBs = [('bbb1',' --binByBinStat --correlateXsecStat '),('bbb0','')]
BBBs = [('bbb1',' --binByBinStat --correlateXsecStat ')]

date = datetime.date.today().isoformat()
absopath  = os.path.abspath('fitsout_'+date)
if options.suffix: absopath += '_'+options.suffix

if options.queue:
    jobdir = absopath+'/jobs/'
    if not os.path.isdir(jobdir):
        os.system('mkdir -p {od}'.format(od=jobdir))
    logdir = absopath+'/logs/'
    if not os.path.isdir(logdir):
        os.system('mkdir -p {od}'.format(od=logdir))
    errdir = absopath+'/errs/'
    if not os.path.isdir(errdir):
        os.system('mkdir -p {od}'.format(od=errdir))
    outdirCondor = absopath+'/outs/'
    if not os.path.isdir(outdirCondor):
        os.system('mkdir -p {od}'.format(od=outdirCondor))

srcfiles = []

# values optimized for el/mu. For lep use muon value
# taureg = 600 if channel=='el' else 932.5 ## values for useExp
taureg = 600 if channel=='el' else 950

threads_opt = ' --nThreads {nt}'.format(nt=options.nThreads) if options.nThreads else ''

for itoy in range(options.toys) if options.toys else range(1):
    if options.toys:
        irand = random.randint(1,1e9)
        expected = [('toy'+str(irand),' -t 1 ')]
        pois = [('poim1','')]
        BBBs = [('bbb1',' --binByBinStat --correlateXsecStat ')]
        seed = irand

    for ipm,POImode in pois:
        card = cardsdir+"/W{chan}_card_withXsecMask_noRochesterStat.hdf5".format(chan=channel) if ipm=='poim1' else cardsdir+'/W{chan}_card_noRochesterStat.hdf5'.format(chan=channel)
        doImpacts = ' --doImpacts ' if ipm=='poim1' else ''
        smoothnessTest = ' --doSmoothnessTest --doh5Output --allowNegativePOI ' if ipm=='poim1' else ''
        regularize = ' --doRegularization --regularizationUseLog --regularizationTau {tau} '.format(tau=taureg) if ipm=='poim1' and options.regularize else ''
        for iexp,exp in expected:
            saveHist = ' --saveHists --computeHistErrors '
            for ibbb,bbb in BBBs:
                pfx = '--postfix {ipm}_{iexp}_{ibbb}'.format(ipm=ipm, iexp=iexp, ibbb=ibbb)
                cmd = 'combinetf.py {poimode} {exp} {bbb} {saveh} {imp} {pfx} {card} {reg} {st} {to} --fitverbose 9'.format(poimode=POImode, exp=exp, bbb=bbb, saveh=saveHist, imp=doImpacts, pfx=pfx, card=card, reg=regularize, st=smoothnessTest, to=threads_opt)
                if options.useSciPy: 
                    cmd += ' --useSciPyMinimizer '
                if options.toys:
                    cmd += ' --seed {irand} '.format(irand=irand)
                if options.queue:
                    job_file_name = jobdir+'jobfit_{ipm}_{iexp}_{ibbb}.sh'.format(ipm=ipm, iexp=iexp, ibbb=ibbb)
                    log_file_name = job_file_name.replace('.sh','.log')
                    tmp_file = open(job_file_name, 'w')
                    tmp_filecont = jobstring_tf
                    tmp_filecont = tmp_filecont.replace('COMBINESTRING', cmd)
                    tmp_filecont = tmp_filecont.replace('CMSSWBASE', os.environ['CMSSW_BASE']+'/src/')
                    tmp_filecont = tmp_filecont.replace('OUTDIR', absopath)
                    tmp_file.write(tmp_filecont)
                    tmp_file.close()
                    srcfiles.append(job_file_name)
                else:
                    print "running ",cmd
                    os.system(cmd)
                    os.system('mv fitresults_123456789.root fitresults_{ipm}_{iexp}_{ibbb}.root'.format(ipm=ipm, iexp=iexp, ibbb=ibbb))
                    print "fit done. Moving to the next one."
                
if options.queue:
    cf = makeCondorFile(jobdir,srcfiles,options, logdir, errdir, outdirCondor)
    subcmd = 'condor_submit {rf} '.format(rf = cf)
    print subcmd

