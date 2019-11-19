import ROOT, copy, math, os
from array import array
from rollingFunctions import roll1Dto2D,dressed2D

ROOT.gROOT.SetBatch()


c = ROOT.TCanvas()
ROOT.gStyle.SetOptStat(0)
ROOT.gROOT.SetBatch()

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
        for pol in pols:
            
            testhistoname = 'x_W{ch}_{pol}'.format(ch=charge,pol=pol)
            
            testhisto = shapefile.Get(testhistoname)
            
            binninPtEtaFile = open('{ind}/binningPtEta.txt'.format(ind=options.indir),'r')
            bins = binninPtEtaFile.readlines()[1].split()[1]
            etabins = list( float(i) for i in bins.replace(' ','').split('*')[0].replace('[','').replace(']','').split(',') )
            ptbins  = list( float(i) for i in bins.replace(' ','').split('*')[1].replace('[','').replace(']','').split(',') )
            nbinseta = len(etabins)-1
            nbinspt  = len( ptbins)-1
            binning = [nbinseta, etabins, nbinspt, ptbins]
            # get eta-pt binning for both reco 

            ## etaPtBinningVec = getDiffXsecBinning(options.indir+'/binningPtEta.txt', "reco")  # this get two vectors with eta and pt binning
            ## recoBins = templateBinning(etaPtBinningVec[0],etaPtBinningVec[1])        # this create a class to manage the binnings
            ## binning = [recoBins.Neta, recoBins.etaBins, recoBins.Npt, recoBins.ptBins]
            
            testhisto_unrolled =  dressed2D(testhisto,binning,'backrolled_testhisto')
            c.cd()
            #histo.SetTitle(testhistoname)
            testhisto_unrolled.Draw('colz')

            testhisto_unrolled.GetZaxis().SetRangeUser(0., 1.1*testhisto_unrolled.GetMaximum())
            c.SaveAs(options.outdir+'/simpleTemplate_'+testhistoname+'.pdf')
            c.SaveAs(options.outdir+'/simpleTemplate_'+testhistoname+'.png')

            projX = testhisto_unrolled.ProjectionX(testhistoname+'_projETA')
            projX.Draw()
            c.SaveAs(options.outdir+'/simpleTemplate_'+projX.GetName()+'.pdf')
            c.SaveAs(options.outdir+'/simpleTemplate_'+projX.GetName()+'.png')

            projY = testhisto_unrolled.ProjectionY(testhistoname+'_projPT')
            projY.Draw()
            c.SaveAs(options.outdir+'/simpleTemplate_'+projY.GetName()+'.pdf')
            c.SaveAs(options.outdir+'/simpleTemplate_'+projY.GetName()+'.png')
