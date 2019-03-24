# USAGE: python runAllFits.py cards_el el
import os,sys

cardsdir = sys.argv[1]
channel = sys.argv[2]
if channel not in ['mu','el','lep']:
    print "Channel must be either mu or el or lep (el-mu combination). Exiting."
    sys.exit()

pois = [('poim1',''),('poim0',' --POIMode none ')]
expected = [('exp1',' -t -1 '),('exp0',' -t 0 ')]
BBBs = [('bbb1',' --binByBinStat --correlateXsecStat '),('bbb0','')]

for ipm,POImode in pois:
    card = cardsdir+"/W{chan}_card_withXsecMask.hdf5".format(chan=channel) if ipm=='poim1' else cardsdir+'/W{chan}_card.hdf5'.format(chan=channel)
    doImpacts = ' --doImpacts ' if ipm=='poim1' else ''
    for iexp,exp in expected:
        saveHist = ' --saveHists --computeHistErrors '
        for ibbb,bbb in BBBs:
            cmd = 'combinetf.py {poimode} {exp} {bbb} {saveh} {imp} {card} --fitverbose 9'.format(poimode=POImode, exp=exp, bbb=bbb, saveh=saveHist, imp=doImpacts, card=card)
            print "running ",cmd
            os.system(cmd)
            os.system('mv fitresults_123456789.root fitresults_{ipm}_{iexp}_{ibbb}.root'.format(ipm=ipm, iexp=iexp, ibbb=ibbb))
            print "fit done. Moving to the next one."

            
