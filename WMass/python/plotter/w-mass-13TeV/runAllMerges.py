import os,sys

if __name__ == "__main__":
    from optparse import OptionParser
    parser = OptionParser(usage='%prog -c <channel> [-m]')
    parser.add_option("-c", "--channel",   dest="channel",  type="string",  default="el",        help="channel (el,mu). Assumes cards dir = cards_<channel>")
    parser.add_option('-m','--merge-root', dest='mergeRoot', default=False, action='store_true', help='Merge the root files with the inputs also')    
    (options, args) = parser.parse_args()

    mergeroot = ' -m ' if options.mergeRoot else ''
    basecmd = 'python mergeCardComponentsAbsY.py -b W{channel} -C plus,minus -i cards_{channel}'.format(channel=options.channel)
    
    cmds = [basecmd + mergeroot,
            basecmd + ' --comb ',
            basecmd + ' --fp ',
            basecmd + ' --fp --comb '
            ]
    
    for c in cmds:
        print "Executing now ",c," ..."
        os.system(c)

