import ROOT, sys, os, itertools, array

## usage:
## python make3D.py /afs/cern.ch/user/m/mdunser/www/private/w-helicity-13TeV/cosThetaFits/2019-07-02-to10/

args = sys.argv

indir = args[1]

files = [f for f in os.listdir(indir) if 'fractions' in f and '.root' in f and not '3D' in f]
files = sorted(files, key = lambda x: int(x.split('_')[-1].replace('.root','').replace('pdf','')) if 'pdf' in x.split('_')[-1] else 0 if 'nominal' in x else -1) ## should sort also alphaS. but maybe not
files = sorted(files, key = lambda x: -1 if not 'alphaSUp' in x else 0)
files = sorted(files, key = lambda x: -1 if not 'alphaSDown' in x else 0)


for i in files:
    print i

tmp_infile = ROOT.TFile(indir+'/'+files[0], 'read')

tmp_hist = tmp_infile.Get('fractionR_plus_0')
binsx = array.array('f', tmp_hist.GetXaxis().GetXbins())
binsy = array.array('f', tmp_hist.GetYaxis().GetXbins())
binsz = array.array('f', [i - 0.5 for i in range(63)] )

tmp_infile.Close()

hists_3d = []

for i,(pol,ch) in enumerate(itertools.product(['R','L','0'],['plus', 'minus'])):

    tmp_hist_3d = ROOT.TH3F('fraction{p}_{ch}'.format(p=pol,ch=ch), 
                            'fraction{p}_{ch}'.format(p=pol,ch=ch), 
                            len(binsx)-1, binsx,
                            len(binsy)-1, binsy,
                            len(binsz)-1, binsz )

    for ip,p in enumerate(files):
        #if ip == 61: continue
        ipdf = 'nominal' if not ip else 'alphaSUp' if ip == 61 else 'alphaSDown' if ip == 62 else 'pdf'+str(ip)

        print 'from file', p, 'trying to get histogram', 'fraction{p}_{ch}_{i}'.format(p=pol,ch=ch,i=ip)
        tmp_f = ROOT.TFile(indir+'/'+p, 'read')
        tmp_h = tmp_f.Get('fraction{p}_{ch}_{i}'.format(p=pol,ch=ch,i=ip))

        for ix, iy in itertools.product(range(1,tmp_h.GetXaxis().GetNbins()+1), range(1,tmp_h.GetYaxis().GetNbins()+1)):
            tmp_hist_3d.SetBinContent(ix, iy, ip+1, tmp_h.GetBinContent(ix, iy))
            tmp_hist_3d.SetBinError  (ix, iy, ip+1, tmp_h.GetBinError  (ix, iy))

        tmp_f.Close()

    hists_3d.append(tmp_hist_3d)

outfile = ROOT.TFile(indir+'/fractions_pdfVariations_3Dhist.root','recreate')
for i in hists_3d:
    i.Write()
outfile.Close()

print 'wrote file', outfile.GetName()
