import ROOT, os, array, re, math


class s:
    def __init__(self, inf):
        self.inf = inf
        self.setValues()

    def setValues(self):
        tmp_f = open(self.inf, 'r')
        tmp_l = tmp_f.readlines()
        for l in tmp_l:
            if 'Signal processes scaled' in l:
                tmp_scale = float(l.split()[-1])
            if 'BACKGROUND' in l:
                tmp_b  = float(l.split()[1])
                tmp_be = float(l.split()[3])
            if 'SIGNAL' in l:
                tmp_s  = float(l.split()[1])
                tmp_se = float(l.split()[3])
            if 'DATA' in l:
                tmp_d  = float(l.split()[1])

        self.sig   = tmp_s
        self.sig_e = tmp_se
        self.bkg   = tmp_b
        self.bkg_e = tmp_be

        self.data  = tmp_d

        self.scale   = tmp_scale
        self.origsig = self.sig/self.scale
        self.newscale = (self.data-self.bkg-self.bkg_e)/self.origsig

        self.scale_e = self.scale-self.newscale
        
        

ROOT.TColor.CreateGradientColorTable(3, array.array ("d", [0.00, 0.50, 1.00]),
                                  ##array ("d", [1.00, 1.00, 0.00]),
                                  ##array ("d", [0.70, 1.00, 0.34]),
                                  ##array ("d", [0.00, 1.00, 0.82]),
                                  array.array ("d", [0.00, 1.00, 1.00]),
                                  array.array ("d", [0.34, 1.00, 0.65]),
                                  array.array ("d", [0.82, 1.00, 0.00]),
                                  255,  0.95)

def getnumbers(s):
    a = re.findall(r'\d+', s)
    a = [float(i) for i in a]
    return a


indir = '/afs/cern.ch/user/m/mdunser/www/private/w-helicity-13TeV/fakeClosureKinematic/'

scales = {}

for sd in os.listdir(indir):
    if not os.path.isdir(indir+'/'+sd): continue
    if not '2019-03-05' in sd: continue
    tmp_fn = indir+'/'+sd+'/etal1.txt'
    vals = s(tmp_fn)
    ## tmp_f =  open(indir+'/'+sd+'/etal1.txt','r')
    ## tmp_l = tmp_f.readlines()
    ## for l in tmp_l:
    ##     if 'Signal processes scaled' in l:
    ##         tmp_s = float(l.split()[-1])
    
    scales[sd.split('-')[-1]] = vals

arr_eta = array.array('d', [0., 0.5, 1.5, 2.4])
arr_pt  = array.array('d', [26, 29, 32, 35, 38, 41, 45])
    
th2_p = ROOT.TH2F('p','scale factors for +', len(arr_eta)-1, arr_eta, len(arr_pt)-1, arr_pt)
th2_m = ROOT.TH2F('m','scale factors for -', len(arr_eta)-1, arr_eta, len(arr_pt)-1, arr_pt)
th2_p .GetXaxis().SetTitle('|#eta|')
th2_m .GetXaxis().SetTitle('|#eta|')
th2_p .GetYaxis().SetTitle('p_{T}')
th2_m .GetYaxis().SetTitle('p_{T}')

for k,v in scales.items():
    pt  = sum(getnumbers(k.split('_')[0]))/2.
    eta = sum(getnumbers(k.split('_')[1].replace('p','')))/20.

    b = th2_p.FindBin(eta, pt)
    if k.split('_')[-1] == 'plus':
        th2_p.SetBinContent(b, v.scale)
        th2_p.SetBinError  (b, v.scale_e)
    else:
        th2_m.SetBinContent(b, v.scale)
        th2_m.SetBinError  (b, v.scale_e)

th2_p.GetZaxis().SetRangeUser(0.65, 1.35)
th2_m.GetZaxis().SetRangeUser(0.65, 1.35)

ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetPaintTextFormat('5.3g')

c = ROOT.TCanvas('c','', 1200, 600)
c.Divide(2,1)
c.cd(1)

th2_p.Draw('colz text e')
c.cd(2)
th2_m.Draw('colz text e')
