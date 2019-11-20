import ROOT, copy, math, os
from array import array
from rollingFunctions import roll1Dto2D,dressed2D

ROOT.gROOT.SetBatch()

colors = [2,3,4,6,7,8,9,46,ROOT.kOrange,ROOT.kViolet]

def histoStyle(histo,ind):
    histo.SetMarkerStyle(20)
    histo.SetMarkerSize(0.9)
    histo.SetMarkerColor(ROOT.kBlack)
    histo.SetLineWidth (2)
    histo.SetLineColor (colors[ind])
    histo.GetXaxis().SetLabelSize(0.045)
    histo.GetYaxis().SetLabelSize(0.045)
    histo.GetXaxis().SetTitleSize(0.050)
    histo.GetYaxis().SetTitleSize(0.045)


c = ROOT.TCanvas('foobar','',1200,900)
c.SetTopMargin    (0.07)
c.SetLeftMargin   (0.10)
c.SetRightMargin  (0.14)
c.SetBottomMargin (0.14)
#pad = ROOT.TPad('foo', '', 0.0,0.0,0.9,0.9)
ROOT.gStyle.SetOptStat(0)
ROOT.gROOT.SetBatch()
#c.Draw()

ROOT.gStyle.SetPalette(ROOT.kThermometer)#LightTemperature)

lat = ROOT.TLatex()
lat.SetNDC()
lat.SetTextSize(0.045)
lat.SetTextFont(42)
#lat.SetFillColor(ROOT.kWhite)


if __name__ == "__main__":
    from optparse import OptionParser
    parser = OptionParser(usage="%prog [options] shapesdir channel")
    parser.add_option('-i','--indir' , default='.', type='string', help='directory that has the shapes files')
    parser.add_option('-o','--outdir', default='.', type='string', help='outdput directory to save the plots')
    (options, args) = parser.parse_args()

    if not os.path.isdir(options.outdir):
        os.system('mkdir -p '+options.outdir)
        os.system('cp ~mdunser/public/index.php '+options.outdir)

    for charge in ['plus', 'minus']:

        shapefile = ROOT.TFile('{ind}/Wmu_{ch}_shapes.root'.format(ind=options.indir,ch=charge),'read')
        
        pols = ['a{p}'.format(p=i) for i in ['c']+range(8)]
        for ip,pol in enumerate(pols):
            
            testhistoname = 'x_W{ch}_{pol}'.format(ch=charge,pol=pol)
            
            testhisto = shapefile.Get(testhistoname)

            integral = testhisto.Integral()
            niceint = '{a:.2f} M events'.format(a=integral/1e6)
            
            binninPtEtaFile = open('{ind}/binningPtEta.txt'.format(ind=options.indir),'r')
            bins = binninPtEtaFile.readlines()[1].split()[1]
            etabins = list( float(i) for i in bins.replace(' ','').split('*')[0].replace('[','').replace(']','').split(',') )
            ptbins  = list( float(i) for i in bins.replace(' ','').split('*')[1].replace('[','').replace(']','').split(',') )
            nbinseta = len(etabins)-1
            nbinspt  = len( ptbins)-1
            binning = [nbinseta, etabins, nbinspt, ptbins]
            
            testhisto_unrolled =  dressed2D(testhisto,binning,'backrolled_testhisto')
            c.cd()


            histoStyle(testhisto_unrolled,0)
            testhisto_unrolled.Draw('colz')

            text = 'W {ch}: {p} - {i}'.format(ch=charge,p=pol if not 'ac' in pol else 'unpolarized',i=niceint)

            lat.SetTextAlign(11)
            lat.DrawLatex(0.10, 0.94, text)

            testhisto_unrolled.GetZaxis().SetRangeUser(0., 1.1*testhisto_unrolled.GetMaximum())
            c.SaveAs(options.outdir+'/simpleTemplate_'+testhistoname+'.pdf')
            c.SaveAs(options.outdir+'/simpleTemplate_'+testhistoname+'.png')

            projX = testhisto_unrolled.ProjectionX(testhistoname+'_projETA')
            histoStyle(projX,ip)
            projX.GetYaxis().SetRangeUser(0., 1.1*projX.GetMaximum())
            projX.Draw('pe')
            lat.SetTextAlign(21)
            lat.DrawLatex(0.50, 0.94, text)

            c.SaveAs(options.outdir+'/simpleTemplate_'+projX.GetName()+'.pdf')
            c.SaveAs(options.outdir+'/simpleTemplate_'+projX.GetName()+'.png')

            projY = testhisto_unrolled.ProjectionY(testhistoname+'_projPT')
            histoStyle(projY,ip)
            projY.GetYaxis().SetRangeUser(0., 1.1*projY.GetMaximum())
            projY.Draw('pe')
            lat.DrawLatex(0.50, 0.94, text)

            c.SaveAs(options.outdir+'/simpleTemplate_'+projY.GetName()+'.pdf')
            c.SaveAs(options.outdir+'/simpleTemplate_'+projY.GetName()+'.png')
