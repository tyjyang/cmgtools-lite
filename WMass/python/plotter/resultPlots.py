import ROOT, optparse, os, re, itertools

## usage:
## python resultPlots.py --config results_config_mu.txt --outdir <outputdirectory> (--runtoys)

if __name__ == '__main__':
    parser = optparse.OptionParser(usage='usage: %prog [opts] ', version='%prog 1.0')
    parser.add_option('', '--config' , type='string'       , default=''    , help='config file with location of fitresults files')
    parser.add_option('', '--outdir' , type='string'       , default=''    , help='output directory with directory structure and plots')
    parser.add_option('', '--runtoys', action='store_true' , default=False , help='run also toys, not only hessian. takes longer')
    parser.add_option('', '--make'   , type='string'       , default='all' , help='run all (default) or only parts (nuis,rap,corr,syst,post,imp)')
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

    charges = ['plus','minus']
    ## plot syst ratios
    ## ================================
    if options.make in ['all', 'syst']:
        print 'running systRatios for a few things'
        tmp_outdir = options.outdir+'/systRatios/'
        os.system('mkdir -p {od}'.format(od=tmp_outdir))
        os.system('cp /afs/cern.ch/user/m/mdunser/public/index.php {od}'.format(od=tmp_outdir))
        systs = []
        systs += ['pdf', 'alphaS', 'mW', 'fsr']
        systs += ['muR'   +str(i) for i in range(1,11)]
        systs += ['muF'   +str(i) for i in range(1,11)]
        systs += ['muRmuF'+str(i) for i in range(1,11)]
        systs += ['smoothelscale{idx}'.format(idx=i) for i in range(0,4)]
        systs += ['smoothmuscale{idx}'.format(idx=i) for i in range(0,2)]
        systs += ['CMS_Wmu_FR_norm']
        systs += ['CMS_Wmu_FRmu_slope']
        systs += ['']
        nTnPUnc = 3 if muEl=='mu' else 4
        systs += ['TnPEffSyst{idx}{flav}'.format(idx=i,flav=muEl) for i in range(0,nTnPUnc)]
        systs += ['L1PrefireEleEffSyst{idx}el'.format(idx=i) for i in range(0,9)]
        systs += ['OutOfAccPrefireSyst{idx}el'.format(idx=i) for i in range(0,2)]
        nEtaUnc = 10 if muEl=='mu' else 26
        systs += [','.join(['FakesEtaUncorrelated{idx}{flav}'.format(idx=i,flav=muEl) for i in xrange(1,nEtaUnc+1)])]
        systs += [','.join(['FakesPtNormUncorrelated{idx}{flav}'.format(idx=i,flav=muEl) for i in xrange(1,nEtaUnc+1)])]
        systs += [','.join(['FakesPtSlopeUncorrelated{idx}{flav}'.format(idx=i,flav=muEl) for i in xrange(1,nEtaUnc+1)])]
        systs += [','.join(['FakesEtaChargeUncorrelated{idx}{flav}{charge}'.format(idx=i,flav=muEl,charge=charge) for i in xrange(1,nEtaUnc+1) for charge in charges])]
        tmp_systs = []
        for p,s,n in itertools.product(['long', 'left', 'right'], ['muR', 'muF', 'muRmuF'], range(1,11)):
            tmp_systs.append('{p}{s}{n}'.format(p=p,s=s,n=n))
        systs += [','.join(tmp_systs)]

        for nuis in systs:
            cmd = 'python w-helicity-13TeV/systRatios.py --unrolled --outdir {od} -s {p} {d} {ch}'.format(od=tmp_outdir, p=nuis, d=results['cardsdir'], ch=muEl)
            print "Running: ",cmd
            os.system(cmd)

    ## plot correlation matrices
    ## ================================
    if options.make in ['all', 'corr']:
        print 'making correlation matrices'
        tmp_outdir = options.outdir+'/correlationMatrices/'
        os.system('mkdir -p {od}'.format(od=tmp_outdir))
        os.system('cp /afs/cern.ch/user/m/mdunser/public/index.php {od}'.format(od=tmp_outdir))
        for t in toysHessian:
            for tmp_file in [i for i in results.keys() if 'both_floatingPOIs_'+t in i]:
                tmp_suffix = '_'.join(tmp_file.split('_')[1:])
                nuisancesAndPOIs = ['CMS_,W.*long', 'pdf', 'muR,muF,muRmuF,alphaS,wpt', 'CMS_', 'ErfPar', 'CMS_,W.*left',  'CMS_,W.*right', 'CMS_,mW,fsr']
                if 'floatingPOIs' in tmp_file:
                    #nuisancesAndPOIs += ['W{charge}_{pol}.*pmaskedexpnorm'.format(charge=charge,pol=pol) for charge in ['plus','minus'] for pol in ['left','right','long'] ]
                    nuisancesAndPOIs += ['W{charge}_right.*pmaskedexpnorm,W{charge}_left.*pmaskedexpnorm'.format(charge=charge) for charge in ['plus','minus'] ]
                for nuis in nuisancesAndPOIs:
                    cmd = 'python w-helicity-13TeV/subMatrix.py {inf} --outdir {od} --params {p} --type {t} --suffix {suf} '.format(od=tmp_outdir, t=t, p=nuis, inf=results[tmp_file], suf=tmp_suffix)
                    print cmd
                    os.system(cmd)

    ## plot rapidity spectra
    ## ==============================================================
    if options.make in ['all', 'rap', 'rapcomp']:
        print 'plotting rapidity spectra'
        tmp_outdir = options.outdir+'/rapiditySpectra/'
        os.system('mkdir -p {od}'.format(od=tmp_outdir))
        os.system('cp /afs/cern.ch/user/m/mdunser/public/index.php {od}'.format(od=tmp_outdir))
        for t in toysHessian:
            for tmp_file in [i for i in results.keys() if 'both_floatingPOIs_'+t in i and not i.endswith('_el') and not i.endswith('_mu')]:
                fitflavor = os.path.splitext(os.path.basename(options.config.split('_')[-1]))[0]
                tmp_suffix = '_'.join(tmp_file.split('_')[1:])
                normstr = [' ', ' --normxsec ']
                if options.make in ['rapcomp','all'] and sum(tmp_file+'_'+lepflav in results.keys() for lepflav in ['el','mu'])==2:
                    print 'NOW plotting combined YW...'
                    cmd  = 'python w-helicity-13TeV/plotYWCompatibility.py '
                    cmd += ' -C plus,minus --xsecfiles {xp},{xm} -y {cd}/binningYW.txt '.format(xp=results['xsecs_plus'],xm=results['xsecs_minus'],cd=results['cardsdir'])
                    cmd += ' --infile-lep {infl} --infile-mu {infmu} --infile-el {infel} --outdir {od} --suffix {suf}  --longBkg '.format(od=tmp_outdir, t=t, suf=tmp_suffix+'_'+fitflavor+'_comp', infl=results[tmp_file],infmu=results[tmp_file+'_mu'],infel=results[tmp_file+'_el'])
                    if 'alt_xsecs_plus' in results.keys() and 'alt_xsecs_minus' in results.keys(): 
                        cmd += ' --altxsecfiles {xp},{xm} '.format(xp=results['alt_xsecs_plus'],xm=results['alt_xsecs_minus'])
                    for norm in normstr:
                        print cmd+norm
                        os.system(cmd+norm)
                # do the single plots for both the single channels and the combined (w/o the el-mu comparison)
                cmd  = 'python w-helicity-13TeV/plotYW.py '
                cmd += ' -C plus,minus --xsecfiles {xp},{xm} -y {cd}/binningYW.txt '.format(xp=results['xsecs_plus'],xm=results['xsecs_minus'],cd=results['cardsdir'])
                cmd += ' --infile {inf} --outdir {od} --type {t} --suffix {suf}  --longBkg '.format(od=tmp_outdir, t=t, suf=tmp_suffix+'_'+fitflavor, inf=results[tmp_file])
                if 'alt_xsecs_plus' in results.keys() and 'alt_xsecs_minus' in results.keys(): 
                    cmd += ' --altxsecfiles {xp},{xm} '.format(xp=results['alt_xsecs_plus'],xm=results['alt_xsecs_minus'])
                for norm in normstr:
                    print cmd+norm
                    os.system(cmd+norm)

    ## plot rapidity spectra for many toys, do not do this by default
    ## ==============================================================
    if options.make in ['manyToys']:
        print 'plotting rapidity spectra'
        tmp_outdir = options.outdir+'/rapiditySpectraManyToys/'
        os.system('mkdir -p {od}'.format(od=tmp_outdir))
        os.system('cp /afs/cern.ch/user/m/mdunser/public/index.php {od}'.format(od=tmp_outdir))
        for t in ['hessian']:
            for tmp_file in os.listdir(results['dirWithToyFits']):
                if not '.root' in tmp_file or not 'fitresults' in tmp_file or not 'toy' in tmp_file: continue
                fitflavor = 'mu'
                tmp_suffix = '_'.join(tmp_file.split('_')[2:]).replace('.root','')
                normstr = [' ', ' --normxsec ']
                cmd  = 'python w-helicity-13TeV/plotYW.py '
                cmd += ' -C plus,minus --xsecfiles {xp},{xm} -y {cd}/binningYW.txt '.format(xp=results['xsecs_plus'],xm=results['xsecs_minus'],cd=results['cardsdir'])
                cmd += ' --infile {inf} --outdir {od} --type {t} --suffix {suf} --nolong '.format(od=tmp_outdir, t=t, suf=tmp_suffix+'_'+fitflavor, inf=results['dirWithToyFits']+'/'+tmp_file)
                for norm in normstr:
                    print cmd+norm
                    os.system(cmd+norm)
                    
    ## plot postfit plots
    ## ================================
    if options.make in ['all', 'post']:
        print 'making postfit plots'
        for tmp_file in [i for i in results.keys() if 'postfit_' in i]:
            tmp_outdir = options.outdir+'/postFitPlots/'
            os.system('mkdir -p {od}'.format(od=tmp_outdir))
            os.system('cp /afs/cern.ch/user/m/mdunser/public/index.php {od}'.format(od=tmp_outdir))
            for reverse in ['',' --reverseUnrolling ']:
                print "NOW RUNNING WITH REVERSE UNROLLING OPTION: ",reverse
                cmd  = 'python w-helicity-13TeV/postFitPlots.py --no2Dplot '
                cmd += ' {inf} {cd} --outdir {od} --suffix {suf} {rev} '.format(inf=results[tmp_file], cd=results['cardsdir'], od=tmp_outdir, suf=tmp_file.replace('postfit',''), rev=reverse)
                os.system(cmd)
            for charge in ['plus','minus']:
                cmdpull = 'python w-helicity-13TeV/monsterPull.py -i {od}/plots_{suf}.root -d unrolled_{ch} --suffix {suf}'.format(od=tmp_outdir,suf=tmp_file.replace('postfit',''),ch=charge)
                os.system(cmdpull)

    ## plot impacts
    ## ================================
    if options.make in ['all','imp','imp1D']:
        print 'making impact plots'
        tmp_outdir = options.outdir+'/impactPlots/'
        poisGroups = ['W{ch}.*{pol}.*'.format(ch=charge,pol=pol) for charge in ['plus','minus'] for pol in ['left','right','long']]
        singleNuisGroups = ['CMS.*','pdf.*','mu.*','ErfPar0.*','ErfPar1.*','ErfPar2.*']
        #targets = ['xsec','xsecnorm','asym','unpolasym','unpolxsec', 'A0', 'A4']
        targets = ['xsec','xsecnorm','asym','unpolasym','unpolxsec','unpolxsecnorm','A0','A4']
        POIsForSummary = {'xsec': '.*', 'xsecnorm': '.*', 'asym': '.*chargeasym', 'unpolasym': '.*chargemetaasym', 'unpolxsec': '.*sumxsec', 'unpolxsecnorm': '.*sumxsecnorm', 'A0': '.*a0', 'A4': '.*a4'}
        for t in toysHessian:
            for tmp_file in [i for i in results.keys() if re.match('both_floatingPOIs_{toyhess}'.format(toyhess=t),i)]:
                tmp_suffix = '_'.join(tmp_file.split('_')[1:])
                if False: #if options.make != 'imp1D':
                    for target in targets[:2]:
                        for poig in poisGroups:
                            print "RUNNING SINGLE NUISANCE IMPACTS..."
                            for nuis in singleNuisGroups:
                                cmd = 'python w-helicity-13TeV/impactPlots.py {fr} -o {od} --nuis {nuis} --pois {pois} --target {tg} --suffix {sfx} --longBkg '.format(fr=results[tmp_file], od=tmp_outdir, nuis=nuis, pois=poig, tg=target, sfx=tmp_suffix)
                                print "running ",cmd
                                os.system(cmd)
                            # now do the groups of systematics
                            print "RUNNING GROUPED NUISANCE IMPACTS..."
                            cmd = 'python w-helicity-13TeV/impactPlots.py {fr} -o {od} --nuisgroups .* --pois {pois} --target {tg} --suffix {sfx} --longBkg '.format(fr=results[tmp_file], od=tmp_outdir, pois=poig, tg=target, sfx=tmp_suffix)
                            print "running ",cmd
                            os.system(cmd)
                            # now make the latex tables for the nuisance groups
                            print "RUNNING TABLES FOR GROUPED NUISANCE IMPACTS..."
                            for charge in ['plus','minus']:
                                cmd = 'python w-helicity-13TeV/impactPlots.py {fr} -o {od} --latex --nuisgroups .* --pois "W{charge}.*(left|right).*(bin_0|bin_4|bin_7|bin_9)" --target {tg} --suffix {sfx} --longBkg '.format(fr=results[tmp_file], od=tmp_outdir, pois=poig, tg=target, sfx=tmp_suffix, charge=charge) 
                                os.system(cmd)
                # now do the 1D summaries
                for target in targets:
                    print "RUNNING 1D SUMMARIES OF SYSTEMATICS..."
                    cmd = 'python w-helicity-13TeV/impactPlots.py {fr} -o {od} --nuisgroups .* --pois {pois} -y {cd}/binningYW.txt --target {tg} --suffix summary_{sfx} --longBkg '.format(fr=results[tmp_file], od=tmp_outdir, pois=POIsForSummary[target], cd=results['cardsdir'], tg=target, sfx=tmp_suffix)
                    if re.match('.*asym|A\d',target): cmd += ' --absolute '
                    print cmd
                    os.system(cmd)

    ## do this at the end, it takes the longest
    ## diff nuisances
    ## ================================
    if options.make in ['all', 'nuis']:
        print 'running diffNuisances'
        tmp_outdir = options.outdir+'/diffNuisances/'
        os.system('mkdir -p {od}'.format(od=tmp_outdir))
        os.system('cp ~emanuele/public/index.php {od}'.format(od=tmp_outdir))
        for t in toysHessian:
            for tmp_file in [i for i in results.keys() if re.match('both_(floating|fixed)POIs_{toyhess}'.format(toyhess=t),i)]:
                tmp_suffix = '_'.join(tmp_file.split('_')[1:])
                nuisancesAndPOIs = ['.*', 'FakesPtNormUncorrelated', 'FakesPtSlopeUncorrelated', 'FakesEtaUncorrelated', 'ZEtaChargeUncorrelated', 'pdf', 'muR,muF,muRmuF,alphaS,wpt,mW', 'CMS_,fsr', 'ErfPar']
                nuisancesAndPOIs += ['longmu', 'leftmu', 'rightmu']
                if 'floatingPOIs' in results[tmp_file]: nuisancesAndPOIs += ['W{charge}_{pol}.*_mu'.format(charge=charge,pol=pol) for charge in ['plus','minus'] for pol in ['left','right','long'] ]
                for nuis in nuisancesAndPOIs:
                    diffNuisances_cmd = 'python w-helicity-13TeV/diffNuisances.py --all --format "html,latex" --outdir {od} --pois {p}'.format(od=tmp_outdir, p=nuis)
                    os.system('{cmd} --infile {inf} --suffix {suf} --type {t} '.format(cmd=diffNuisances_cmd, inf=results[tmp_file], suf=tmp_suffix, t=t))
                    print '{cmd} --infile {inf} --suffix {suf} --type {t} '.format(cmd=diffNuisances_cmd, inf=results[tmp_file], suf=tmp_suffix, t=t)
