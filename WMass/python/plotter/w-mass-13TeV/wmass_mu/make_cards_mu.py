import os
from datetime import datetime

PROG         = 'w-mass-13TeV/make_wmass_cards.py'
BASECONFIG   = 'w-mass-13TeV/testingNano/cfg'
MCA          = BASECONFIG+'/mca-testHistForCard.txt'
CUTFILE      = BASECONFIG+'/cuts_testHistForCard.txt'
CUTFILE_WLIKE = BASECONFIG+'/cuts_testHistForCard_wlike.txt'
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
WEIGHTSTRING = ' \'puWeight*_get_muonSF_fast_wmass(Muon_pt[0],Muon_eta[0],Muon_charge[0])*PrefireWeight\' '  # _get_muonSF_fast_wmass applies trigger and selection SF all at once
WEIGHTSTRING_WLIKE = ' \'puWeight*_get_muonSF_fast_wlike(Muon_pt[0],Muon_eta[0],Muon_charge[0],Muon_pt[1],Muon_eta[1],Muon_charge[1],event)*PrefireWeight\' ' 
OUTDIR       = 'wmass_%s' % datetime.now().strftime('%Y_%m_%d')
    
if __name__ == '__main__':
    from optparse import OptionParser
    parser = OptionParser(usage='%prog [options]')
    parser.add_option('-s', '--suffix' , dest='suffix' , type='string'      , default=None , help='Append a suffix to the default outputdir (wmass_<date>)');
    parser.add_option("-o", "--outdir", dest="outdir", type="string", default='', help="Name of output folder. It will ignore option -s and the automatic date in the name. This option is particulaly useful to produce output in an already existing folder (e.g. when you do signal and background in different days, which would create a new folder). Warning: with condor it might overwrite things!!!");
    parser.add_option('-d', '--dry-run', dest='dryRun' , action='store_true', default=False, help='Do not run the job, only print the command');
    parser.add_option("-p", "--print-only", dest="printOnly",   action="store_true", default=False, help="Just print commands, do not execute anything in this script");
    parser.add_option("--syst"         , dest="addSyst", action="store_true", default=False, help="Add PDF systematics to the signal (need incl_sig directive in the MCA file)");
    parser.add_option("-r", "--run", dest="run", type="string", default="sb", help="Which components to run: s for signal, b for backgrounds or sb for both");
    parser.add_option("--max-genWeight", dest="maxGenWeight", type="string", default="50000.0", help="Maximum gen weight to be used for Z and W samples (with any decay). Weights larger than this value will be set to it.");
    parser.add_option('--wlike', dest='wlike', action="store_true", default=False, help="Make cards for the wlike analysis. Default is wmass");
    (options, args) = parser.parse_args()
    
    if options.suffix: 
        OUTDIR += ('_%s' % options.suffix)
    if options.outdir != "": 
        OUTDIR = options.outdir

    components = []
    if "s" in options.run:
        components.append(" -s ")
    if "b" in options.run:
        components.append(" -b ")

    if options.wlike:
        CUTFILE = CUTFILE_WLIKE
        WEIGHTSTRING = WEIGHTSTRING_WLIKE

    for c in components:
        cmd='python ' + ' '.join([PROG,MCA,CUTFILE,VAR,BINNING,SYSTFILE,OUTDIR,'-C mu']) + \
            (' -W %s ' % WEIGHTSTRING) + (' -P %s ' % TREEPATH) + (' -q %s ' % QUEUE) + c
        if options.dryRun: cmd += '  --dry-run '
        if options.addSyst: cmd += '  --pdf-syst --qcd-syst '        
        if options.wlike: cmd += ' --wlike '
        if 'b' in c:
            cmd += ' -n bkg '
        elif 's' in c:
            cmd += ' -n sig '
        optsForHisto = ' --nanoaod-tree --max-genWeight-procs \'W.*|Z.*\' \'{m}\' --clip-genWeight-toMax '.format(m=options.maxGenWeight)
        cmd += ' --add-option "{opt}"'.format(opt=optsForHisto)
        cmd += ' -g 5 '
        cmd += ' --decorrelateSignalScales '
        cmd += ' --vpt-weight Zmumu --vpt-weight Ztautau' #--vpt-weight Wmunu --vpt-weight Wtaunu ' # check if easier to use regular expressions to catch cases inside PROG        
        if options.printOnly:
            print cmd
            print ""
        else:
            os.system(cmd)
