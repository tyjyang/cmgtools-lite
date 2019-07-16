import os,sys,re

from optparse import OptionParser
parser = OptionParser(usage="%prog [options] cardsdir")
parser.add_option("--cel", dest="cardsEl",  type="string", default="cards_el", help="Directory with electron cards");
parser.add_option("--cmu", dest="cardsMu",  type="string", default="cards_mu", help="Directory with muon cards");
(options, args) = parser.parse_args()

outdir = args[0]

suff = ['','_withXsecMask']
#maskedChannels = ['InAcc','OutAcc']
maskedChannels = ['InAcc'] # if the Y bins 10,11 are treated as bkg, then should not be masked channels
inputdirs = {'el': options.cardsEl, 'mu': options.cardsMu}

# gymnastics:
# it needs to combine the el and mu cards with xsecmasks which are the only ones that contain the sum and polgroups
# then it needs not to duplicate the xsec masked channels, so include the muon only one
# then, needs to replace the xsec file with one containing twice the events, in order --binByBinStat --correlateXsecStat to work properly

print "Starting the cards combination..."
for s in suff:
    freezePOI = True if len(s)==0 else False
    print "## Start the combination for freezePOI = ",freezePOI
    datacards=[]
    for flav in ['el','mu']:
        datacards.append('{inputdir}/W{flav}_card{suffix}.txt'.format(inputdir=os.path.abspath(inputdirs[flav]),flav=flav,suffix=s))
    combinedCard = outdir+"/Wlep_card{suffix}.txt".format(suffix=s)
    os.system('rm -f '+combinedCard)
    if freezePOI:
        ccCmd = 'combineCards.py --noDirPrefix '+' '.join(['{dcfile}'.format(dcfile=dc) for dc in datacards])+' > '+combinedCard
        os.system(ccCmd)
        txt2hdf5Cmd = 'text2hdf5.py '+combinedCard
    else:
        ccCmd = 'combineCards.py --noDirPrefix '+' '.join(['{dcfile}'.format(dcfile=dc) for dc in datacards])+' --xc .*_xsec_.* xsecs_plus_InAcc={cmu}/Wmu_plus_xsec_InAcc_card.txt xsecs_minus_InAcc={cmu}/Wmu_minus_xsec_InAcc_card.txt > '.format(cmu=options.cardsMu) + combinedCard
        ## here running the combine cards command first
        print ccCmd+'\n'
        os.system(ccCmd)
        ## now replace the single flav xsec with twice it
        cardlines = [line.rstrip('\n') for line in open(combinedCard,'r')]
        for l in cardlines:
            match = re.match('.*xsecs_(\S+)_InAcc.*',l)
            if l.startswith("shapes") and match:
                charge = match.group(1)
                rfile = l.split()[3]
                newpath = os.path.abspath('{clep}/Wlep_{ch}_shapes_xsec.root'.format(clep=outdir,ch=charge))
                haddCmd = 'hadd -f {newpath}'.format(newpath=newpath)+' '+' '.join([rfile for i in xrange(2)])
                os.system(haddCmd)
                ## this is only for plotting Y spectra
                haddCmd = haddCmd.replace('shapes_xsec.root','shapes_xsec_baremc.root')
                os.system(haddCmd)
                os.system('sed -i "s|{origpath}|{newpath}|g" {card} '.format(origpath=rfile,newpath=newpath,card=combinedCard))
        maskchan = [' --maskedChan xsecs_{ch}_InAcc'.format(ch=ch) for ch in ['plus','minus']]
        txt2hdf5Cmd = 'text2hdf5.py {maskch} --X-allow-no-background {cf}'.format(maskch=' '.join(maskchan),cf=combinedCard)
    ## here making the TF meta file
    print txt2hdf5Cmd+'\n\n'
    os.system(txt2hdf5Cmd)
