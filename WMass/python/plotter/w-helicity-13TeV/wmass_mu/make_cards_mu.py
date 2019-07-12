import os
from datetime import datetime

PROG         = 'w-helicity-13TeV/make_helicity_cards.py'
BASECONFIG   = 'w-helicity-13TeV/wmass_mu'
MCA          = BASECONFIG+'/mca-wmu-helicity.txt'
CUTFILE      = BASECONFIG+'/cuts_wmu.txt'
SYSTFILE     = BASECONFIG+'/systsEnv.txt'

QUEUE        = '2nd'
VAR          = '\'ptMuFull(LepGood1_calPt,LepGood1_eta):LepGood1_eta\''

TREEPATH     = '/afs/cern.ch/work/m/mdunser/public/wmassTrees/SKIMS_muons_latest/' ## this makes cremedelacreme
# TREEPATH     = '/afs/cern.ch/work/m/mdunser/public/wmassTrees/SKIMS_muons_withFSR/' ## this is marco's skim

binningeta = [-2.4 + i*0.1 for i in range(49) ]
#binningeta = [-2.2 + i*0.2 for i in range(23) ]
binningeta = [float('{a:.3f}'.format(a=i)) for i in binningeta]

etabinning = '['+','.join('{a:.1f}'.format(a=i) for i in binningeta)+']'

## variable binning in pt
ptbinning = '['+','.join(str(i) for i in range(26,46))+']'
#ptbinning = '['+','.join(str(i) for i in range(30,46))+']'
#ptbinning = '['+','.join(str(26+1.5*i) for i in range(0,14))+']'

BINNING      = '\''+etabinning+'*'+ptbinning+'\''
## do with histogram. sick of friends !! WEIGHTSTRING = ' \'puw2016_nTrueInt_36fb(nTrueInt)*LepGood_SF1[0]*LepGood_SF2[0]\' '
WEIGHTSTRING = ' \'puw2016_nTrueInt_36fb(nTrueInt)*_get_muonSF_selectionToTrigger(LepGood_pdgId[0],LepGood_calPt[0],LepGood_eta[0],LepGood_charge[0])*LepGood_SF2[0]*prefireJetsWeight(LepGood_eta[0])\' '
OUTDIR       = 'helicity_%s' % datetime.now().strftime('%Y_%m_%d')

#components=[' -b ']#, ' -s ']
components=[' -s ', ' -b ']
    

if __name__ == '__main__':
    from optparse import OptionParser
    parser = OptionParser(usage='%prog [options]')
    parser.add_option('-d', '--dry-run', dest='dryRun',   action='store_true', default=False, help='Do not run the job, only print the command');
    parser.add_option('-s', '--suffix', dest='suffix', type='string', default=None, help='Append a suffix to the default outputdir (helicity_<date>)');
    parser.add_option("--syst", dest="addSyst", action="store_true", default=False, help="Add PDF systematics to the signal (need incl_sig directive in the MCA file)");
    parser.add_option("--long-bkg", dest="longBkg",    action="store_true", default=False, help="Treat the longitudinal polarization as one background template.");
    parser.add_option("--genw", action="store_true", default=False, help="use genw (dressed leptons) instead of prefsrw.");
    (options, args) = parser.parse_args()
    
    if options.suffix: OUTDIR += ('_%s' % options.suffix)

    for c in components:
        cmd='python ' + ' '.join([PROG,MCA,CUTFILE,VAR,BINNING,SYSTFILE,OUTDIR,'-C mu']) + \
            (' -W %s ' % WEIGHTSTRING) + (' -P %s ' % TREEPATH) + (' -q %s ' % QUEUE) + c
        if options.dryRun: cmd += '  --dry-run '
        if options.addSyst: cmd += '  --pdf-syst --qcd-syst '
        if options.longBkg: cmd += ' --long-bkg '
        if not options.genw: cmd += ' --wvar prefsrw '
        cmd += ' -g 5 '
        cmd += ' --decorrelateSignalScales '
        cmd += ' --vpt-weight Z '#--vpt-weight W --vpt-weight TauDecaysW '
        os.system(cmd)
