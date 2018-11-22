import ROOT, optparse, os

## usage:
## python resultPlots.py --config results_config_mu.txt --outdir <outputdirectory> (--runtoys)

if __name__ == '__main__':
    parser = optparse.OptionParser(usage='usage: %prog [opts] ', version='%prog 1.0')
    parser.add_option('', '--config' , type='string'       , default=''    , help='config file with location of fitresults files')
    parser.add_option('', '--outdir' , type='string'       , default=''    , help='output directory with directory structure and plots')
    parser.add_option('', '--runtoys', action='store_true' , default=False , help='run also toys, not only hessian. takes longer')
    parser.add_option('', '--make'   , type='string'       , default='all' , help='run all (default) or only parts (nuis,rap,corr,syst)')
    (options, args) = parser.parse_args()


    ## this loads the path for all the results in a config
    results = {}
    execfile(options.config, results)

    ## make the output directory first
    os.system('mkdir -p {od}'.format(od=options.outdir))

    ## get the channel from the filename of the config.txt
    muEl = 'mu' if 'mu' in options.config else 'el'

    ## check if running also toys
    toysHessian = ['hessian']
    if options.runtoys:
        toysHessian += ['toys']

    ## plot syst ratios
    ## ================================
    if options.make in ['all', 'syst']:
        print 'running systRatios for a few things'
        tmp_outdir = options.outdir+'/systRatios/'
        os.system('mkdir -p {od}'.format(od=tmp_outdir))
        os.system('cp ~mdunser/public/index.php {od}'.format(od=tmp_outdir))
        systs  = ['pdf', 'alphaS']
        systs += ['muR'   +str(i) for i in range(1,11)]
        systs += ['muF'   +str(i) for i in range(1,11)]
        systs += ['muRmuF'+str(i) for i in range(1,11)]
        for nuis in systs:
            os.system('python w-helicity-13TeV/systRatios.py --outdir {od} -s {p} {d} {ch}'.format(od=tmp_outdir, p=nuis, d=results['cardsdir'], ch=muEl))

    ## plot correlation matrices
    ## ================================
    if options.make in ['all', 'corr']:
        print 'making correlation matrices'
        tmp_outdir = options.outdir+'/correlationMatrices/'
        os.system('mkdir -p {od}'.format(od=tmp_outdir))
        os.system('cp ~mdunser/public/index.php {od}'.format(od=tmp_outdir))
        for t in toysHessian:
            for nuis in ['pdf', 'muR,muF,muRmuF,alphaS,wpt', 'CMS_']:
                basecmd = 'python w-helicity-13TeV/subMatrix.py '
                os.system(basecmd+' {inf} --outdir {od} --params {p} --type {t} --suffix floatingPOIs_{t} '.format(od=tmp_outdir, t=t, p=nuis, inf=results['both_floatingPOIs_'+t]))
                os.system(basecmd+' {inf} --outdir {od} --params {p} --type {t} --suffix fixedPOIs_{t} '   .format(od=tmp_outdir, t=t, p=nuis, inf=results['both_fixedPOIs_'+t]))

    ## plot rapidity spectra
    ## ================================
    if options.make in ['all', 'rap']:
        print 'plotting rapidity spectra'
        tmp_outdir = options.outdir+'/rapiditySpectra/'
        os.system('mkdir -p {od}'.format(od=tmp_outdir))
        os.system('cp ~mdunser/public/index.php {od}'.format(od=tmp_outdir))
        for t in toysHessian:
            cmd  = 'python w-helicity-13TeV/plotYW.py '
            cmd += ' -C plus,minus --xsecfiles {xp},{xm} '.format(xp=results['xsecs_plus'],xm=results['xsecs_minus'])
            cmd += ' -y {cd}/binningYW.txt --infile {inf} --outdir {od} --type {t} --suffix floatingPOIs_{t} '.format(cd=results['cardsdir'],od=tmp_outdir, t=t, inf=results['both_floatingPOIs_'+t])
            os.system(cmd)
            print "===> plotting normalized xsecs..."
            cmd += ' --normxsec '
            os.system(cmd)

    ## do this at the end, it takes the longest
    ## diff nuisances
    ## ================================
    if options.make in ['all', 'nuis']:
        print 'running diffNuisances'
        tmp_outdir = options.outdir+'/diffNuisances/'
        os.system('mkdir -p {od}'.format(od=tmp_outdir))
        os.system('cp ~mdunser/public/index.php {od}'.format(od=tmp_outdir))
        for t in toysHessian:
            for nuis in ['pdf', 'muR,muF,wpt', 'CMS_', 'ErfPar0EffStat', 'ErfPar1EffStat', 'ErfPar2EffStat']:
                diffNuisances_cmd = 'python w-helicity-13TeV/diffNuisances.py --outdir {od} --pois {p}'.format(od=tmp_outdir, p=nuis)
                os.system('{cmd} --infile {inf} --suffix floatingPOIs_{t} --type {t} '.format(cmd=diffNuisances_cmd, inf=results['both_floatingPOIs_'+t], t=t))
                os.system('{cmd} --infile {inf} --suffix fixedPOIs_{t}    --type {t} '.format(cmd=diffNuisances_cmd, inf=results['both_fixedPOIs_'+t]   , t=t))

