#!/usr/bin/env python
# may use as: cat w-helicity-13TeV/wmass_e/zee_catlist.txt | xargs -i python w-helicity-13TeV/wmass_e/wmass_plots.py plots/testZskim {} > runplots.sh
# may use as: cat w-helicity-13TeV/wmass_e/wgen_catlist.txt | xargs -i python w-helicity-13TeV/wmass_e/wmass_plots.py plots/gen {} > runplots.sh

# to submit a subset of the plots in the plots.txt on condor may do:

import sys,re,os,datetime
sys.path.insert(0,os.path.abspath(os.getcwd()+"/w-helicity-13TeV"))
from submitToys import makeCondorFile,jobstring_tf

ODIR=sys.argv[1]

FASTTEST=''
#FASTTEST='--max-entries 1000 '

dowhat = "plots" 
#dowhat = "dumps" 
#dowhat = "yields" 

TREES = "-F Friends '{P}/friends/tree_Friend_{cname}.root' "
TREESONLYSKIMW = "-P /afs/cern.ch/work/e/emanuele/TREES/TREES_electrons_1l_V6_TINY"
TREESONLYSKIMZ = "-P /data1/emanuele/wmass/TREES_2018-07-17-recoLeptons" # trees with the Trigger Match object
TREESONLYFULL  = "-P /eos/cms/store/group/dpg_ecal/comm_ecal/localreco/TREES_1LEP_80X_V3"

def base(selection,useSkim=True):

    if 'wenu' in selection: TREESONLYSKIM=TREESONLYSKIMW
    elif 'wgen' in selection: TREESONLYSKIM=TREESONLYSKIMW
    elif 'zee' in selection: TREESONLYSKIM=TREESONLYSKIMZ
    else:
        raise RuntimeError, 'Unknown selection %s' % selection

    CORE=' '.join([TREES,TREESONLYSKIM if useSkim else TREESONLYFULL])
    if 'cmsphys06' in os.environ['HOSTNAME']: CORE = CORE.replace('/eos/cms/store/group/dpg_ecal/comm_ecal/localreco/','/data1/emanuele/wmass/')

    CORE+=" -f -j 8 -l 35.9 --s2v "+FASTTEST
    if dowhat == "plots": 
        CORE+=" --lspam '#bf{CMS} #it{Preliminary}' --legendWidth 0.20 --legendFontSize 0.035 "
        if selection != "wgen":
            CORE+=" --showRatio --maxRatioRange 0.90 1.10 --fixRatioRange "

    if selection=='wenu':
        GO="%s w-helicity-13TeV/wmass_e/mca-80X-wenu-simpleplots.txt w-helicity-13TeV/wmass_e/wenu_80X.txt "%CORE
        GO="%s -W 'puw2016_nTrueInt_36fb(nTrueInt)*lepSF(LepGood1_pdgId,LepGood1_pt,LepGood1_eta,LepGood1_SF1,LepGood1_SF2,LepGood1_SF3)*prefireJetsWeight(LepGood_eta[0])'"%GO
        if dowhat in ["plots","ntuple"]: 
            plotfile = "wenu_plots.txt "
            GO+=" w-helicity-13TeV/wmass_e/"+wenu_plots.txt
    elif selection=='wgen':
        GO="%s -W 1 "%CORE
        GO="%s w-helicity-13TeV/wmass_e/mca-80X-wenu-helicity.txt w-helicity-13TeV/wmass_e/wenu_80X.txt "%GO
        GO="%s -p Wplus_long,Wplus_left,Wplus_right --plotmode=nostack "%GO
        #GO="%s --sP wplus_mtwtk,wplus_etal1,wplus_ptl1,wplus_etal1gen,wplus_ptl1gen,wplus_wy,wplus_wpt "%GO
        if dowhat in ["plots","ntuple"]: GO+=" w-helicity-13TeV/wmass_e/wenu_plots.txt "        
    elif selection=='zee':
        GO="%s w-helicity-13TeV/wmass_e/mca-80X-zee.txt w-helicity-13TeV/wmass_e/zee.txt --fitData --flp Z "%CORE
        GO="%s -W 'puw2016_nTrueInt_36fb(nTrueInt)*eleSF_HLT_2lfriends(LepGood1_matchedTrgObjElePt,LepGood1_SF1,LepGood2_matchedTrgObjElePt,LepGood2_SF1)*LepGood1_SF2*LepGood1_SF3*LepGood2_SF2*LepGood2_SF3*eleSF_L1Eff(LepGood1_pt,LepGood1_eta)*eleSF_L1Eff(LepGood2_pt,LepGood2_eta)' --sp 'Z' "%GO
        if dowhat in ["plots","ntuple"]: GO+=" w-helicity-13TeV/wmass_e/zee_plots.txt "
    else:
        raise RuntimeError, 'Unknown selection'

    return GO

def procs(GO,mylist):
    return GO+' '+" ".join([ '-p %s'%l for l in mylist ])
def sigprocs(GO,mylist):
    return procs(GO,mylist)+' --showIndivSigs --noStackSig'
def runIt(GO,name,plots=[],noplots=[]):
    if   dowhat == "plots":  print 'python mcPlots.py',"--pdir %s/%s"%(ODIR,name),GO,' '.join(['--sP \'%s\''%p for p in plots]),' '.join(['--xP \'%s\''%p for p in noplots])
    elif dowhat == "yields": print 'echo %s; python mcAnalysis.py'%name,GO
    elif dowhat == "dumps":  print 'echo %s; python mcDump.py'%name,GO
    elif dowhat == "ntuple": print 'echo %s; python mcNtuple.py'%name,GO
def submitIt(GO,name,plots=[],noplots=[],opts=None):

    date = datetime.date.today().isoformat()
    jobdir= os.path.abspath("./jobsplots_{name}_{date}/".format(name=name,date=date))
    if not os.path.isdir(jobdir): os.mkdir(jobdir)

    srcfiles = []
    for i,pl in enumerate(plots):
        cmd = "python mcPlots.py --pdir {pdir} --sP {plot} -o {pdir}/{plot}.root {GO}".format(pdir=os.path.abspath(ODIR+"/"+name),plot=pl,GO=GO)
        if 'plus' in pl: cmd = swapcharge(cmd,'plus')
        elif 'minus' in pl: cmd = swapcharge(cmd,'minus')

        job_file_name = jobdir+'/plot_{i}.sh'.format(i=i)
        log_file_name = job_file_name.replace('.sh','.log')
        tmp_file = open(job_file_name, 'w')
        tmp_filecont = jobstring_tf
        tmp_filecont = tmp_filecont.replace('COMBINESTRING', cmd)
        tmp_filecont = tmp_filecont.replace('CMSSWBASE', os.environ['CMSSW_BASE']+'/src/')
        tmp_filecont = tmp_filecont.replace('OUTDIR', os.path.abspath(os.getcwd()))
        tmp_file.write(tmp_filecont)
        tmp_file.close()
        srcfiles.append(job_file_name)
    cf = makeCondorFile(jobdir,srcfiles,opts,jobdir,jobdir,jobdir)
    print "condor file is in ",jobdir+"/condor_submit.condor"
def add(GO,opt):
    return '%s %s'%(GO,opt)
def remove(GO,opt):
    return GO.replace(opt,'')
def swapcharge(GO,charge):
    return GO.replace('plus',charge)
def setwide(x):
    x2 = add(x,'--wide')
    x2 = x2.replace('--legendWidth 0.35','--legendWidth 0.20')
    return x2
def fulltrees(x,selection):
    if 'wenu' in selection: TREESONLYSKIM=TREESONLYSKIMW
    elif 'zee' in selection: TREESONLYSKIM=TREESONLYSKIMZ
    else:
        raise RuntimeError, 'Unknown selection'
    return x.replace(TREESONLYSKIM,TREESONLYFULL)

allow_unblinding = True

if __name__ == '__main__':

    from optparse import OptionParser
    parser = OptionParser(usage="%prog outdir what [options] ")
    parser.add_option("--dry-run", dest="dryRun",    action="store_true", default=False, help="Do not run the job, only print the command");
    parser.add_option('-r', '--runtime',    default=0, type=int,   help='New runtime for condor resubmission in hours. default: 0 (it means interactively)');
    (options, args) = parser.parse_args()

    torun = args[1]

    if (not allow_unblinding) and '_data' in torun and (not any([re.match(x.strip()+'$',torun) for x in ['.*_appl.*','cr_.*']])): raise RuntimeError, 'You are trying to unblind!'

    x=""
    if 'zee' in torun:
        x = base('zee')
        if '_incl' in torun: x = add(x," --sP 'etalep,ptlep,zmass_egm,zmass_corr' ")
        if '_ebeb' in torun: x = add(x,"-A alwaystrue ebeb 'max(abs(LepGood1_eta),abs(LepGood2_eta))<1.44'  --sP 'zmass_egm,zmass_corr' ")
        if '_notebeb' in torun: x = add(x,"-A alwaystrue notebeb 'max(abs(LepGood1_eta),abs(LepGood2_eta))>1.57' --sP 'zmass_egm,zmass_corr' ")
        if '_gg' in torun: x = add(x,"-A alwaystrue goldgold 'min(LepGood1_r9,LepGood2_r9)>0.94' --sP 'zmass_egm,zmass_corr' ")
        if '_notgg' in torun: x = add(x,"-A alwaystrue notgoldgold 'min(LepGood1_r9,LepGood2_r9)<0.94' --sP 'zmass_egm,zmass_corr' ")
        #if '_w_reweight' in torun and dowhat=="plots": x = add(x,"--sP 'z_mll,pt1,pt2,ptZ,scaledptZ,costheta_cs,phi_cs,sumAiPi,y_vs_ctheta,y_vs_phi,y_vs_sumAiPi' ")
        #if '_genpt' in torun: x = add(x,"--sP 'gen_ptv,gen_scaledptv' --xp 'data' -p 'Z' ")
    elif 'wenu' in torun:
        x = base('wenu')
        x = x = add(x," --sP 'etalep,ptlep' ")
    elif 'wgen' in torun:
        x = base('wgen')
        if 'nosel' in torun:
            x = add(x," -U 'alwaystrue ' ")
        if 'fullsel' in torun:
            x = add(x," -W 'puw2016_nTrueInt_36fb(nTrueInt)*trgSF_We(LepGood1_pdgId,LepGood1_pt,LepGood1_eta,2)*leptonSF_We(LepGood1_pdgId,LepGood1_pt,LepGood1_eta)' ")
        if torun.endswith('pdfs'):
            for i in xrange(60):
                wgt = 'hessWgt[%d]/genWeight' % i
                output = '%s/%s_pdf%i.root' % (ODIR,torun,i)
                xw = x
                xw = add(xw," -W '{wgt}' -o {output}".format(wgt=wgt,output=output))
                submitIt( xw,'%s_%s' % (torun,i), opts=options )
        if 'qcdpostfit' in torun:
            #x = x.replace('mca-80X-wenu-helicity.txt','mca-80X-wenu-helicity-postfit.txt')
            poldic = {'long':0,'left':1,'right':2}
            chargedic = {'plus':1,'minus':-1}
            scalings = ''
            for charge,sign in chargedic.iteritems():
                for pol,ipol in poldic.iteritems():
                    scalings += '--scale-process W{charge}_{pol} "postfitQCDWeight(pt_2(GenLepDressed_pt[0],GenLepDressed_phi[0],GenPromptNu_pt[0],GenPromptNu_phi[0]),{ipol},{sign})" '.format(charge=charge,pol=pol,ipol=ipol,sign=sign)
            x = add(x,scalings)

    plotlist = args[2] # if empty, to all the ones of the txt file
    if not torun.endswith('pdfs'): 
        plotlistfile = open(plotlist,'r')
        plots = [p.strip() for p in plotlistfile.readlines()]
        print "I will select these plots from the plotfile:"
        print plots
        if options.runtime>0: submitIt(x,'%s'%torun,plots,opts=options)
        else: runIt(x,'%s'%torun,plots)
