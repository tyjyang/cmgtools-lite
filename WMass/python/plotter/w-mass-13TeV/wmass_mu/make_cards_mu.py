import os
from datetime import datetime

PROG         = 'w-mass-13TeV/make_wmass_cards.py'
BASECONFIG   = 'w-mass-13TeV/testingNano/cfg'
#MCA          = BASECONFIG+'/mca-test.txt'
MCA          = BASECONFIG+'/mca-testHistForCard.txt'
#CUTFILE      = BASECONFIG+'/cuts_test.txt'
CUTFILE      = BASECONFIG+'/cuts_testHistForCard.txt'
SYSTFILE     = BASECONFIG+'/systsEnv.txt'

QUEUE        = '1nd'
VAR          = '\'Muon_pt[0]:Muon_eta[0]\''

TREEPATH     = '/eos/cms/store/cmst3/group/wmass/w-mass-13TeV/postNANO/dec2020/'

binningeta = [-2.4 + i*0.1 for i in range(49) ]
binningeta = [float('{a:.3f}'.format(a=i)) for i in binningeta]

etabinning = '['+','.join('{a:.1f}'.format(a=i) for i in binningeta)+']'

## variable binning in pt
ptbinning = '['+','.join(str(i) for i in range(26,46))+']'

BINNING      = '\''+etabinning+'*'+ptbinning+'\''
WEIGHTSTRING = ' \'puw2016_nTrueInt_36fb(Pileup_nTrueInt)*PrefireWeight\' ' # to be updated
OUTDIR       = 'wmass_%s' % datetime.now().strftime('%Y_%m_%d')
    
if __name__ == '__main__':
    from optparse import OptionParser
    parser = OptionParser(usage='%prog [options]')
    parser.add_option('-s', '--suffix' , dest='suffix' , type='string'      , default=None , help='Append a suffix to the default outputdir (wmass_<date>)');
    parser.add_option('-d', '--dry-run', dest='dryRun' , action='store_true', default=False, help='Do not run the job, only print the command');
    parser.add_option("--syst"         , dest="addSyst", action="store_true", default=False, help="Add PDF systematics to the signal (need incl_sig directive in the MCA file)");
    parser.add_option("--genw"                         , action="store_true", default=False, help="use genw (dressed leptons) instead of prefsrw.");
    parser.add_option("-r", "--run", dest="run", type="string", default="sb", help="Which components to run: s for signal, b for backgrounds or sb for both");
    parser.add_option("--max-genWeight", dest="maxGenWeight", type="string", default="50000.0", help="Maximum gen weight to be used for Z and W samples (with any decay). Weights larger than this value will be set to it.");
    (options, args) = parser.parse_args()
    
    if options.suffix: OUTDIR += ('_%s' % options.suffix)

    components = []
    if "s" in options.run:
        components.append(" -s ")
    if "b" in options.run:
        components.append(" -b ")


    for c in components:
        cmd='python ' + ' '.join([PROG,MCA,CUTFILE,VAR,BINNING,SYSTFILE,OUTDIR,'-C mu']) + \
            (' -W %s ' % WEIGHTSTRING) + (' -P %s ' % TREEPATH) + (' -q %s ' % QUEUE) + c
        if options.dryRun: cmd += '  --dry-run '
        if options.addSyst: cmd += '  --pdf-syst --qcd-syst '
        if not options.genw: cmd += ' --wvar prefsrw '
        cmd += ' --nanoaod-tree --max-genWeight-procs "W.*|Z.*" "{m}" --clip-genWeight-toMax '.format(m=options.maxGenWeight)
        cmd += ' -g 5 '
        cmd += ' --decorrelateSignalScales '
        cmd += ' --vpt-weight Z '#--vpt-weight W --vpt-weight TauDecaysW '
        os.system(cmd)
