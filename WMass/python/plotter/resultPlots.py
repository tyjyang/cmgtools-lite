import ROOT, optparse, os, re, itertools

## usage:
## python resultPlots.py --config results_config_mu.txt --outdir <outputdirectory> (--runtoys)

if __name__ == '__main__':
    parser = optparse.OptionParser(usage='usage: %prog [opts] ', version='%prog 1.0')
    parser.add_option('', '--config' , type='string'       , default=''    , help='config file with location of fitresults files')
    parser.add_option('', '--outdir' , type='string'       , default=''    , help='output directory with directory structure and plots')
    parser.add_option('', '--runtoys', action='store_true' , default=False , help='run also toys, not only hessian. takes longer')
    parser.add_option('', '--make'   , type='string'       , default='all' , help='run all (default) or only parts (nuis,rap,corr,syst,post,imp,templ)')
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

    if options.make in ['all','templ']:
        print 'make the nominal templates'
        tmp_outdir = options.outdir+'/templates2D/'
        os.system('mkdir -p {od}'.format(od=tmp_outdir))
        os.system('cp /afs/cern.ch/user/m/mdunser/public/index.php {od}'.format(od=tmp_outdir))
        cmd = 'python w-helicity-13TeV/make2DTemplates.py --outdir {od} {d} {ch}'.format(od=tmp_outdir,d=results['cardsdir'], ch=muEl) 
        os.system(cmd)

    charges = ['plus','minus']
    ## plot syst ratios
    ## ================================
    if options.make in ['all', 'syst']:
        print 'running systRatios for a few things'
        tmp_outdir = options.outdir+'/systRatios/'
        os.system('mkdir -p {od}'.format(od=tmp_outdir))
        os.system('cp /afs/cern.ch/user/m/mdunser/public/index.php {od}'.format(od=tmp_outdir))
        systs = []
        systs += ['pdf', 'alphaS']
        systs += ['mW', 'smoothfsr']
        systs += ['{pol}muR{i}{charge}'   .format(pol=pol,i=i,charge=charge) for i in range(1,11) for pol in ['left','right','long'] for charge in ['plus','minus']]
        systs += ['{pol}muF{i}{charge}'   .format(pol=pol,i=i,charge=charge) for i in range(1,11) for pol in ['left','right','long'] for charge in ['plus','minus']]
        systs += ['{pol}muRmuF{i}{charge}'.format(pol=pol,i=i,charge=charge) for i in range(1,11) for pol in ['left','right','long'] for charge in ['plus','minus']]
        systs += ['smoothelscaleStat{idx}'.format(idx=i) for i in range(0,98)] + ['smoothelscaleSyst{idx}eta{ieta}pt{ipt}{ch}'.format(idx=i,ieta=ieta,ipt=ipt,ch=charge) for i in range(2) for ieta in range(8) for ipt in range(2) for charge in ['plus','minus']] 
        systs += ['smoothmuscaleStat{idx}'.format(idx=i) for i in range(0,99)] + ['smoothmuscaleSyst{idx}eta{ieta}{ch}'.format(idx=i,ieta=ieta,ch=charge) for i in range(2,6) for ieta in range(4) for charge in (['plus','minus'] if i!=3 else [''])]
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
        for p,s,n,c in itertools.product(['long', 'left', 'right'], ['muR', 'muF', 'muRmuF'], range(1,11), ['plus','minus']):
            tmp_systs.append('{p}{s}{n}{c}'.format(p=p,s=s,n=n,c=c))
        systs += [','.join(tmp_systs)]
        for nuis in systs:
            cmd = 'python w-helicity-13TeV/systRatios.py --unrolled -a --outdir {od} -s {p} {d} {ch}'.format(od=tmp_outdir, p=nuis, d=results['cardsdir'], ch=muEl)
            print "Running: ",cmd
            os.system(cmd)
            # print 'running systRatios for a few single rapidity  bins (0,8,9)'
            # os.system('mkdir -p {od}/singleRapidity'.format(od=tmp_outdir))
            # os.system('cp /afs/cern.ch/user/m/mdunser/public/index.php {od}/singleRapidity'.format(od=tmp_outdir))
            # for rapBin in [0,8,9]:
            #     cmd = 'python w-helicity-13TeV/systRatios.py --unrolled --outdir {od}/singleRapidity -s {p} {d} {ch} --singleRap {iy}'.format(od=tmp_outdir, p=nuis, iy=rapBin, d=results['cardsdir'], ch=muEl)
            #     print "Now running: ",cmd
            #     os.system(cmd)

    ## plot correlation matrices
    ## ================================
    if options.make in ['all', 'corr']:
        print 'making correlation matrices'
        tmp_outdir = options.outdir+'/correlationMatrices/'
        os.system('mkdir -p {od}'.format(od=tmp_outdir))
        os.system('cp /afs/cern.ch/user/m/mdunser/public/index.php {od}'.format(od=tmp_outdir))
        for t in toysHessian:
            for tmp_file in [i for i in results.keys() if 'both_floatingPOIs_'+t in i or 'both_fixedPOIs_'+t in i]:
                tmp_suffix = '_'.join(tmp_file.split('_')[1:])
                nuisancesAndPOIs = ['pdf,alphaS', 'CMS_,W.*long', 'muR,muF,muRmuF,wpt', 'CMS_', 'ErfPar', 'CMS_,W.*left',  'CMS_,W.*right', 'CMS_,mW,fsr']
                if 'floatingPOIs' in tmp_file:
                    #nuisancesAndPOIs += ['W{charge}_{pol}.*pmaskedexpnorm'.format(charge=charge,pol=pol) for charge in ['plus','minus'] for pol in ['left','right','long'] ]
                    nuisancesAndPOIs += ['W{charge}_right.*pmaskedexpnorm,W{charge}_left.*pmaskedexpnorm'.format(charge=charge) for charge in ['plus','minus'] ]
                for nuis in nuisancesAndPOIs:
                    cmd = 'python w-helicity-13TeV/subMatrix.py {inf} --outdir {od} --params {p} --type {t} --suffix {suf} '.format(od=tmp_outdir, t=t, p=nuis, inf=results[tmp_file], suf=tmp_suffix)
                    if 'fixedPOIs' in tmp_file:
                        cmd += ' -m channelnone '
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
                nuisancesAndPOIs = []
                nuisancesAndPOIs += ['.*', 'FakesPtNormUncorrelated', 'FakesPtSlopeUncorrelated', 'FakesEtaUncorrelated', 'ZEtaChargeUncorrelated', 'pdf,alphaS', 'mW,smoothfsr,CMS_lumi_13TeV', 'ErfPar']
                nuisancesAndPOIs += ['longmu', 'leftmu', 'rightmu']
                nuisancesAndPOIs += ['norm_W']
                if 'floatingPOIs' in results[tmp_file]: nuisancesAndPOIs += ['W{charge}_{pol}.*_mu'.format(charge=charge,pol=pol) for charge in ['plus','minus'] for pol in ['left','right','long'] ]
                for nuis in nuisancesAndPOIs:
                    diffNuisances_cmd = 'python w-helicity-13TeV/diffNuisances.py --all --format "html,latex" --outdir {od} --pois {p}'.format(od=tmp_outdir, p=nuis)
                    os.system('{cmd} --infile {inf} --suffix {suf} --type {t} '.format(cmd=diffNuisances_cmd, inf=results[tmp_file], suf=tmp_suffix, t=t))
                    print '{cmd} --infile {inf} --suffix {suf} --type {t} '.format(cmd=diffNuisances_cmd, inf=results[tmp_file], suf=tmp_suffix, t=t)
        ## these are the refinement plots for the paper
        os.system('python w-helicity-13TeV/plotExpObsPull.py --exp {od}/nuisances_pdfalphaS_fixedPOIs_hessian_bbb1_syst1_asimov.latex --obs {od}/nuisances_pdfalphaS_fixedPOIs_hessian_bbb1_syst1_data.latex --outdir {od}'.format(od=tmp_outdir))
        for pol in ['left','right']:
            os.system('python w-helicity-13TeV/plotExpObsPull.py --exp {od}/nuisances_{pol}mu_floatingPOIs_hessian_bbb1_syst1_asimov.latex --obs {od}/nuisances_{pol}mu_floatingPOIs_hessian_bbb1_syst1_data.latex --outdir {od}'.format(pol=pol,od=tmp_outdir))
        ## these are to check the compatibility of common uncertainties between the channels
        for group in ['normW','mWsmoothfsrCMSlumi13TeV','pdfalphaS','mWsmoothfsrCMSlumi13TeV'] + ['{pol}mu'.format(pol=pol) for pol in 'left','right','long']:
            os.system('python w-helicity-13TeV/compareNuisancesComb.py --el {od}/nuisances_{gr}_floatingPOIs_hessian_bbb1_syst1_data_el.latex --mu {od}/nuisances_{gr}_floatingPOIs_hessian_bbb1_syst1_data_mu.latex --lep {od}/nuisances_{gr}_floatingPOIs_hessian_bbb1_syst1_data.latex --outdir {od}'.format(gr=group,od=tmp_outdir))
