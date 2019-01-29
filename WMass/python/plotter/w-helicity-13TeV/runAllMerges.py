# USAGE: python runAllMerges.py cards_el el
import os,sys

cardsdir = sys.argv[1]
channel = sys.argv[2]
if channel not in ['mu','el']:
    print "Channel must be either mu or el. Exiting."
    sys.exit()

fps = ['',' --fp ']
charges = [' -C plus,minus', '--comb' ]

for ch in charges:
    for fp in fps:
        cmd = './mergeCardComponentsAbsY.py -b W{channel} {ch} -i {cdir} {fp}'.format(channel=channel,cdir=cardsdir,ch=ch,fp=fp)
        print cmd
        os.system(cmd)
