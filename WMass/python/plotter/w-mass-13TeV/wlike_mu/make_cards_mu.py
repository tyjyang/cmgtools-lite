import os
from datetime import datetime

PROG         = 'w-mass-13TeV/make_wmass_cards.py'
BASECONFIG   = 'w-mass-13TeV/wlike_mu'
MCA          = BASECONFIG+'/mca-wlike-mass.txt'
CUTFILE      = BASECONFIG+'/cuts_wmu.txt'
SYSTFILE     = BASECONFIG+'/systsEnv.txt'

QUEUE        = '1nd'
#VAR          = '\'ptMuFull(returnChargeVal(LepGood1_kalPt,LepGood1_charge,LepGood2_kalPt,LepGood2_charge,evt),returnChargeVal(LepGood1_eta,LepGood1_charge,LepGood2_eta,LepGood2_charge,evt)):returnChargeVal(LepGood1_eta,LepGood1_charge,LepGood2_eta,LepGood2_charge,evt)\''

VAR          = '\'returnChargeVal(LepGood1_pt,LepGood1_charge,LepGood2_pt,LepGood2_charge,evt):returnChargeVal(LepGood1_eta,LepGood1_charge,LepGood2_eta,LepGood2_charge,evt)\''

TREEPATH     = '/afs/cern.ch/work/e/emanuele/TREES/SKIM_2LEP_wlike_mu_V2/'

binningeta = [-2.4 + i*0.1 for i in range(49) ]
binningeta = [float('{a:.3f}'.format(a=i)) for i in binningeta]

etabinning = '['+','.join('{a:.1f}'.format(a=i) for i in binningeta)+']'

## variable binning in pt
ptbinning = '['+','.join(str(i) for i in range(26,54))+']'  # for Wlike 46 is too low

BINNING      = '\''+etabinning+'*'+ptbinning+'\''
## change _pt => _kalPt if the SFs are finally done binning in the KaMuCa muon pt
WEIGHTSTRING = ' \'puw2016_nTrueInt_36fb(nTrueInt)*_get_muonSF_recoToSelection(LepGood1_pdgId,LepGood1_pt,LepGood1_eta)*_get_muonSF_recoToSelection(LepGood2_pdgId,LepGood2_pt,LepGood2_eta)*prefireJetsWeight(LepGood1_eta)*prefireJetsWeight(LepGood2_eta)\' '

OUTDIR       = 'wlike_%s' % datetime.now().strftime('%Y_%m_%d')
    

if __name__ == '__main__':
    from optparse import OptionParser
    parser = OptionParser(usage='%prog [options]')
    parser.add_option('-s', '--suffix' , dest='suffix' , type='string'      , default=None , help='Append a suffix to the default outputdir (helicity_<date>)');
    parser.add_option('-d', '--dry-run', dest='dryRun' , action='store_true', default=False, help='Do not run the job, only print the command');
    parser.add_option("--syst"         , dest="addSyst", action="store_true", default=False, help="Add PDF systematics to the signal (need incl_sig directive in the MCA file)");
    parser.add_option("--genw"                         , action="store_true", default=False, help="use genw (dressed leptons) instead of prefsrw.");
    parser.add_option("-r", "--run", dest="run", type="string", default="sb", help="Which components to run: s for signal, b for backgrounds or sb for both");
    parser.add_option("-p", "--print-only", dest="printOnly",   action="store_true", default=False, help="Just print commands of this script, do not execute anything");
    (options, args) = parser.parse_args()
    
    if options.suffix: OUTDIR += ('_%s' % options.suffix)

    components=[]
    if "s" in options.run:
        components.append(" -s ")
    if "b" in options.run:
        components.append(" -b ")

    print ""
    print "================================="
    print "RECO BINNING (eta - pt)"
    print "---------------------------------"
    print BINNING
    print 'eta:', etabinning
    print 'pt :', ptbinning
    print "================================="
    print ""

    for c in components:
        cmd='python ' + ' '.join([PROG,MCA,CUTFILE,VAR,BINNING,SYSTFILE,OUTDIR,'-C mu']) + \
            (' -W %s ' % WEIGHTSTRING) + (' -P %s ' % TREEPATH) + (' -q %s ' % QUEUE) + c
        if options.dryRun: cmd += '  --dry-run '
        if options.addSyst: cmd += '  --pdf-syst --qcd-syst --qed-syst --kamuca-syst '
        if not options.genw: cmd += ' --wvar prefsrw '
        cmd += ' -g 30 '
        cmd += ' --decorrelateSignalScales '
        cmd += ' --vpt-weight Z --vpt-weight W --vpt-weight TauDecaysW '
        cmd += ' --wlike '
        cmd += ' --inclusive '
        if not options.printOnly:
            os.system(cmd)
        else:
            print cmd
            print ""
        
        print ""

