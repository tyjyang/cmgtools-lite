import ROOT, itertools, copy, os, sys, math
from array import array

ROOT.gROOT.SetBatch()
ROOT.gStyle.SetOptStat(0)

ROOT.TColor.CreateGradientColorTable(3,
                                  array ("d", [0.00, 0.50, 1.00]),
                                  array ("d", [0.00, 1.00, 1.00]),
                                  array ("d", [0.34, 1.00, 0.65]),
                                  array ("d", [0.82, 1.00, 0.00]),
                                  255,  0.95)


args = sys.argv

basedir = str(args[1]) #'/afs/cern.ch/user/m/mdunser/www/private/w-helicity-13TeV/cosThetaFits/2019-06-18-pdfFractions-to5-multiCore/' #multiCore/'
os.system('cp ~/public/index.php '+basedir)

files = {}

for i in os.listdir(basedir):
    if not '.root' in i: continue
    pdfvar = i.split('_')[-1].replace('.root','')
    files[pdfvar] = basedir+'/'+i
    

f_nom = ROOT.TFile(files['nominal'], 'read')
nominals = {}

for pol, ch in itertools.product(['R','L','0'],['plus', 'minus']):
    tmp_name = 'fraction{p}_{c}_0'  .format(p=pol, c=ch)
    tmp_hist = f_nom.Get(tmp_name)
    nominals[tmp_name] = copy.deepcopy(tmp_hist.Clone('nom_'+tmp_name))


rangex = range(1,nominals['fractionL_plus_0'].GetXaxis().GetNbins()+1)
rangey = range(1,nominals['fractionL_plus_0'].GetYaxis().GetNbins()+1)

sumsq = {}
for pol, ch in itertools.product(['R','L','0'],['plus', 'minus']):
    sumsq[pol+ch] = []
    for ix,iy in itertools.product(rangex,rangey):
        sumsq[pol+ch].append(0.)

##f_nom.Close()

print 'this is nominals', nominals


for ipdf in range(1,61):
    print 'at ipdf', ipdf
    tmp_filename = files['pdf'+str(ipdf)]

    ##if not os.path.isfile(tmp_filename):
    ##    print 'file not found... check it!'
    ##    print '       ', tmp_filename
    ##    continue
    f_pdf = ROOT.TFile(tmp_filename, 'read')

    for pol, ch in itertools.product(['R','L','0'],['plus', 'minus']):
        tmp_nom = copy.deepcopy(nominals['fraction{p}_{c}_0'  .format(p=pol, c=ch)] )

        tmp_pdf = f_pdf.Get('fraction{p}_{c}_{i}'.format(p=pol, c=ch, i=ipdf))

        tmp_ratio = tmp_pdf.Clone('ratio_pdf_'+str(ipdf))
        tmp_ratio.SetTitle('ratio for f_{{{p}}}^{{{c}}} and pdf {i}'.format(p=pol, c=ch, i=ipdf))

        tmp_ratio.Reset()

        for ib,(ix,iy) in enumerate(itertools.product(rangex,rangey)):
            try:
                rat = tmp_pdf.GetBinContent(ix, iy) / tmp_nom.GetBinContent(ix, iy)
            except:
                rat = 1.
            sumsq[pol+ch][ib] = sumsq[pol+ch][ib] + (1.-rat)**2
            tmp_ratio.SetBinContent(ix, iy, rat)

        ## print 'bin content before', tmp_ratio.GetBinContent(1,1)
        ## tmp_ratio.Divide(tmp_nom)
        ## print 'bin content after ', tmp_ratio.GetBinContent(1,1)
        canv = ROOT.TCanvas('pdf'+str(ipdf), 'pdf'+str(ipdf), 1200, 800)
        canv.cd()
        canv.SetRightMargin (0.15)
        canv.SetLeftMargin  (0.15)
        canv.SetBottomMargin(0.15)
        tmp_ratio.GetZaxis().SetRangeUser(0.95, 1.05)

        tmp_ratio.Draw('colz')
        tmp_ratio.GetZaxis().SetRangeUser(0.95, 1.05)

        canv.Update()
        canv.SaveAs(basedir+'/ratioToNominal_fraction{p}_{c}_pdf{i}.pdf'.format(p=pol, c=ch, i=ipdf))
        canv.SaveAs(basedir+'/ratioToNominal_fraction{p}_{c}_pdf{i}.png'.format(p=pol, c=ch, i=ipdf))

print '====================================='
#print sumsq
print '====================================='

    ##f_pdf.Close()
for pol, ch in itertools.product(['R','L','0'],['plus', 'minus']):
    tmp_nom = copy.deepcopy(nominals['fraction{p}_{c}_0'  .format(p=pol, c=ch)] )

    tmp_ratio = tmp_nom.Clone('ratio_pdf_'+str(ipdf))
    tmp_ratio.SetTitle('#sum_{ipdf} ( #frac{f_{ipdf} }{f_{nom} } )^{2} for f_{'+pol+'}^{'+('+' if ch == 'plus' else '-')+'}')
    tmp_ratio.Reset()

    for ib,(ix,iy) in enumerate(itertools.product(rangex,rangey)):
        tmp_ratio.SetBinContent(ix, iy, 1.+math.sqrt(sumsq[pol+ch][ib]))

    canv = ROOT.TCanvas('sumsqratio'+pol+ch, 'sumsqratio'+pol+ch, 1200, 800)
    canv.cd()
    canv.SetRightMargin (0.15)
    canv.SetLeftMargin  (0.15)
    canv.SetBottomMargin(0.15)
    tmp_ratio.GetZaxis().SetRangeUser(0.95, 1.05)

    tmp_ratio.Draw('colz')

    canv.Update()
    canv.SaveAs(basedir+'/ratioToNominal_fraction{p}_{c}_sumSquared.pdf'.format(p=pol, c=ch))
    canv.SaveAs(basedir+'/ratioToNominal_fraction{p}_{c}_sumSquared.png'.format(p=pol, c=ch))
