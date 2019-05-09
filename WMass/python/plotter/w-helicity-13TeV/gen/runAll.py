import os
maindir = '/eos/cms/store/user/arapyan/mc'

channel = 'e' # change this to "e" or "mu"

for ch in ['p','m']:
    for gen in ['pythia8','photos']:
        cmd = 'python photos.py -c {channel} --indir {main}/w{charge}_{channel}nu_{generator}/GEN w{charge}_{channel}nu_{generator}.root >& w{charge}_{channel}nu_{generator}.log &'.format(main=maindir,charge=ch,generator=gen,channel=channel)
        print cmd
