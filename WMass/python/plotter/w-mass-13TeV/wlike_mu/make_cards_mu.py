import os
from datetime import datetime

PROG         = 'w-mass-13TeV/make_wmass_cards.py'
BASECONFIG   = 'w-mass-13TeV/wlike_mu'
MCA          = BASECONFIG+'/mca-wlike-mass.txt'
CUTFILE      = BASECONFIG+'/cuts_wmu.txt'
SYSTFILE     = BASECONFIG+'/systsEnv.txt'

QUEUE        = '1nd'
VAR          = '\'ptMuFull(returnChargeVal(LepGood1_calPt,LepGood1_charge,LepGood2_calPt,LepGood2_charge,evt),returnChargeVal(LepGood1_eta,LepGood1_charge,LepGood2_eta,LepGood2_charge,evt)):returnChargeVal(LepGood1_eta,LepGood1_charge,LepGood2_eta,LepGood2_charge,evt)\''

TREEPATH     = '/afs/cern.ch/work/e/emanuele/TREES/SKIM_2LEP_wlike_mu_V1/'

binningeta = [-2.4 + i*0.1 for i in range(49) ]
binningeta = [float('{a:.3f}'.format(a=i)) for i in binningeta]

etabinning = '['+','.join('{a:.1f}'.format(a=i) for i in binningeta)+']'

## variable binning in pt
ptbinning = '['+','.join(str(i) for i in range(26,46))+']'

BINNING      = '\''+etabinning+'*'+ptbinning+'\''
## do with histogram. sick of friends !! WEIGHTSTRING = ' \'puw2016_nTrueInt_36fb(nTrueInt)*LepGood_SF1[0]*LepGood_SF2[0]\' '
WEIGHTSTRING = ' \'puw2016_nTrueInt_36fb(nTrueInt)*_get_muonSF_recoToSelection(LepGood_pdgId[0],LepGood_calPt[0],LepGood_eta[0])*_get_muonSF_recoToSelection(LepGood_pdgId[1],LepGood_calPt[1],LepGood_eta[1])*prefireJetsWeight(LepGood_eta[0])*prefireJetsWeight(LepGood_eta[1])\' '

OUTDIR       = 'wlike_%s' % datetime.now().strftime('%Y_%m_%d')

components=[' -s ', ' -b ']
    

if __name__ == '__main__':
    from optparse import OptionParser
    parser = OptionParser(usage='%prog [options]')
    parser.add_option('-s', '--suffix' , dest='suffix' , type='string'      , default=None , help='Append a suffix to the default outputdir (helicity_<date>)');
    parser.add_option('-d', '--dry-run', dest='dryRun' , action='store_true', default=False, help='Do not run the job, only print the command');
    parser.add_option("--syst"         , dest="addSyst", action="store_true", default=False, help="Add PDF systematics to the signal (need incl_sig directive in the MCA file)");
    parser.add_option("--genw"                         , action="store_true", default=False, help="use genw (dressed leptons) instead of prefsrw.");
    (options, args) = parser.parse_args()
    
    if options.suffix: OUTDIR += ('_%s' % options.suffix)

    for c in components:
        cmd='python ' + ' '.join([PROG,MCA,CUTFILE,VAR,BINNING,SYSTFILE,OUTDIR,'-C mu']) + \
            (' -W %s ' % WEIGHTSTRING) + (' -P %s ' % TREEPATH) + (' -q %s ' % QUEUE) + c
        if options.dryRun: cmd += '  --dry-run '
        if options.addSyst: cmd += '  --pdf-syst --qcd-syst --qed-syst '
        if not options.genw: cmd += ' --wvar prefsrw '
        cmd += ' -g 30 '
        cmd += ' --decorrelateSignalScales '
        cmd += ' --vpt-weight Z --vpt-weight W --vpt-weight TauDecaysW '
        cmd += ' --wlike '
        os.system(cmd)
