import os,sys,glob
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import ROOT

from optparse import OptionParser
parser = OptionParser(usage="%prog [options] ")
parser.add_option("-c", "--channel",    default='mu', dest="channel", type="string",   help="Can be mu or el");
parser.add_option("-s", "--suffix",     default='',   dest="suffix",  type="string",   help="Add asuffix to the output file");
(options, args) = parser.parse_args()

histo = 'correlation_matrix_channelpmaskedexp'
# this corresponds to:
# 4 = W+L, Ybin = 1
# 5 = W+L, Ybin = 2
# chosen to have the most of the stat and to have spills from 2 adjacent bins
covelement = (4,5)


fitresults = glob.glob('fitresults_123456789*.root')

x = []; y = []
for i,f in enumerate(fitresults):
    rf = ROOT.TFile.Open(f)
    fitres = rf.Get('fitresults')
    br = fitres.GetBranch("taureg")
    br.GetEntry(0)
    tau = fitres.taureg
    corrmat = rf.Get(histo)
    # check that the bins are always what we think they are
    lab1 = corrmat.GetXaxis().GetBinLabel(covelement[0])
    lab2 = corrmat.GetXaxis().GetBinLabel(covelement[1])
    if lab1!='Wplus_left_Ybin_1_pmaskedexp' or lab2!='Wplus_left_Ybin_2_pmaskedexp':
        print "ERROR! They are not the bins you think!"
        sys.exit(1)
    corr = corrmat.GetBinContent(covelement[0],covelement[1])
    print "point i = ",i," has tau = ",tau," with corr = ",corr
    x.append(tau)
    y.append(corr)

plt.figure(figsize=(8, 8))
plt.plot(x,y,'bs')
plt.xlabel(r'regularization $\tau$',fontsize=16)
plt.ylabel(r'correlation $W^+_L$ Ybin=1 vs Ybin=2',fontsize=16)
plt.grid(True)
if options.suffix:  sfx = '_'+options.suffix
for e in ['png','pdf']: 
    plt.savefig("tauscan_{flav}{suffix}.{ext}".format(ext=e,flav=options.channel,suffix=sfx))

