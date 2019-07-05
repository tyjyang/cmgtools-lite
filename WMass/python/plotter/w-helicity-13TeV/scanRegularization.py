import os
from submitToys import makeCondorFile,jobstring_tf

from optparse import OptionParser
parser = OptionParser(usage="%prog [options] ")
parser.add_option('-r', '--runtime',    default=1, type=int,   help='New runtime for condor resubmission in hours. default: 24h (combined fits may be long)')
parser.add_option("-q", "--queue",      default=False, dest="queue",  action="store_true",   help="Run jobs on condor instead of locally");
parser.add_option("-c", "--channel",    default='mu', dest="channel", type="string",   help="Can be mu or el");
parser.add_option(      "--useExp",     default=False, dest="useExp", action="store_true",   help="The range of the scan changes drastically if regularizationUseExpected");
(options, args) = parser.parse_args()

cwd = os.getcwd()

if options.useExp:
    points = [5*i for i in xrange(1,201)]
    expStr = ' --regularizationUseExpected '
else:
    points = [2*i for i in xrange(1,601)]
    expStr = ''

srcfiles = []
for  tau in points:
    cmd = 'combinetf.py -t -1 --binByBinStat --correlateXsecStat --postfix poim1_exp1_bbb1_tau{reg:.2f} W{flav}_card_withXsecMask_syst0.hdf5 --doRegularization --regularizationUseLog {expected} --regularizationTau {reg} --fitverbose 9'.format(reg=tau,flav=options.channel,expected=expStr)
    print "==> now testing tau = ",tau
    print cmd
    if options.queue:
        job_file_name = 'jobfit_tau{reg:.2f}.sh'.format(reg=tau)
        log_file_name = job_file_name.replace('.sh','.log')
        tmp_file = open(job_file_name, 'w')
        tmp_filecont = jobstring_tf
        tmp_filecont = tmp_filecont.replace('COMBINESTRING', cmd)
        tmp_filecont = tmp_filecont.replace('CMSSWBASE', os.environ['CMSSW_BASE']+'/src/')
        tmp_filecont = tmp_filecont.replace('OUTDIR', cwd)
        tmp_file.write(tmp_filecont)
        tmp_file.close()
        srcfiles.append(job_file_name)
    else:
        os.system(cmd)
    
if options.queue:
    cf = makeCondorFile(cwd,srcfiles,options,cwd,cwd,cwd)
    subcmd = 'condor_submit {rf} '.format(rf = cf)
    print subcmd

