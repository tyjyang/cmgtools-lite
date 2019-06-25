import os,sys

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

print "Starting the cards combination..."
for s in suff:
    freezePOI = True if len(s)==0 else False
    print "## Start the combination for freezePOI = ",freezePOI
    datacards=[]
    for flav in ['el','mu']:
        datacards.append('{inputdir}/W{flav}_card{suffix}.txt'.format(inputdir=os.path.abspath(inputdirs[flav]),flav=flav,suffix=s))
    combinedCard = outdir+"/Wlep_card{suffix}.txt".format(suffix=s)
    ccCmd = 'combineCards.py --noDirPrefix '+' '.join(['{dcfile}'.format(dcfile=dc) for dc in datacards])+' > '+combinedCard
    if freezePOI:
        txt2hdf5Cmd = 'text2hdf5.py '+combinedCard
    else:
        maskchan = [' --maskedChan ch{icha}_W{flav}_{charge}_xsec_{mc}'.format(icha=iflav+1,flav=flav,charge=ch,mc=mc) for iflav,flav in enumerate(['el','mu']) for ch in ['plus','minus'] for mc in maskedChannels]
        txt2hdf5Cmd = 'text2hdf5.py {maskch} --X-allow-no-background {cf}'.format(maskch=' '.join(maskchan),cf=combinedCard)
    ## here running the combine cards command first
    print ccCmd+'\n'
    os.system(ccCmd)
    ## here making the TF meta file
    print txt2hdf5Cmd+'\n\n'
    os.system(txt2hdf5Cmd)
