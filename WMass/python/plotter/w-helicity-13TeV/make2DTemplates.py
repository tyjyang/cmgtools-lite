import ROOT, copy, math
from array import array
from rollingFunctions import roll1Dto2D,dressed2D

ROOT.gROOT.SetBatch()


c = ROOT.TCanvas()
ROOT.gStyle.SetOptStat(0)
ROOT.gROOT.SetBatch()

if __name__ == "__main__":
    from optparse import OptionParser
    parser = OptionParser(usage="%prog [options] shapesdir channel")
    parser.add_option('-o','--outdir', dest='outdir', default='.', type='string', help='outdput directory to save the plots')
    (options, args) = parser.parse_args()

    for charge in ['plus', 'minus']:

        shapefile = ROOT.TFile('/afs/cern.ch/work/e/emanuele/wmass/heppy/CMSSW_8_0_25/src/CMGTools/WMass/python/plotter/w-helicity-13TeV/cards_mu/Wmu_{ch}_shapes.root'.format(ch=charge),'read')
        shapesfile = ROOT.TFile("{indir}/W{flav}_{ch}_shapes.root".format(indir=args[0],flav=args[1],ch=charge))
        
        for pol in ['right', 'left', 'long']:
            for iybin,ybin in enumerate(range(12)):
            
                testhistoname = 'x_W{ch}_{pol}_Ybin_{ybin}'.format(ch=charge,pol=pol,ybin=ybin)
            
                testhisto = shapefile.Get(testhistoname)
            
                if not iybin:
                    binninPtEtaFile = open('/afs/cern.ch/work/e/emanuele/wmass/heppy/CMSSW_8_0_25/src/CMGTools/WMass/python/plotter/w-helicity-13TeV/cards_mu/binningPtEta.txt','r')
                    bins = binninPtEtaFile.readlines()[1].split()[1]
                    ## hack. easier
                    etabins = list( float(i) for i in bins.replace(' ','').split('*')[0].replace('[','').replace(']','').split(',') )
                    ptbins  = list( float(i) for i in bins.replace(' ','').split('*')[1].replace('[','').replace(']','').split(',') )
                    nbinseta = len(etabins)-1
                    nbinspt  = len( ptbins)-1
                    binning = [nbinseta, etabins, nbinspt, ptbins]
            
                testhisto_unrolled =  dressed2D(testhisto,binning,'backrolled_testhisto')
                c.cd()
                #histo.SetTitle(testhistoname)
                testhisto_unrolled.Draw('colz')
                c.SaveAs(options.outdir+'/simpleTemplate_'+testhistoname+'.pdf')
                c.SaveAs(options.outdir+'/simpleTemplate_'+testhistoname+'.png')
