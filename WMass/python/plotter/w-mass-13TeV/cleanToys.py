import os
import ROOT

if __name__ == "__main__":

    from optparse import OptionParser
    parser = OptionParser(usage='%prog -d toysdir')
    parser.add_option('-d'  , '--dir'  , dest='dirToys'           , type='string'           , default='./' , help='Directory with toys')
    (options, args) = parser.parse_args()

    dirWithToys = options.dirToys
    absDirWithToys = os.path.abspath(dirWithToys)

    cleanDir = 'clean_'+dirWithToys
    srcfile = open('clean.sh','w')
    srcfile.write("echo Creating new dir with links in "+cleanDir+"\n")
    srcfile.write('mkdir '+cleanDir+'\n')
    srcfile.write('cd '+cleanDir+'\n')
    srcfile.write('linkChunks.sh '+absDirWithToys+'\n')
    srcfile.write('echo now checking bad toys...\n')
    checkScript = os.environ['CMSSW_BASE']+'/src/CMGTools/WMass/scripts/checkRootFiles.py'
    badtoys = 'badtoys.txt'
    srcfile.write('python {script} -d {toydir} -f {outtxt}\n'.format(script=checkScript,toydir=absDirWithToys,outtxt=badtoys))
    srcfile.write('echo "List of bad toys in '+badtoys+'"\n')
    srcfile.write('source '+badtoys+'\n')
    srcfile.write('rm '+badtoys+'\n')
    srcfile.write('cd - \n')
    srcfile.write('echo "The dir '+cleanDir+' has only the good ROOT files of '+dirWithToys+'"\n')
    srcfile.write("echo ENJOY.\n")
    srcfile.close()
    print "Now source clean.sh"

    
              
