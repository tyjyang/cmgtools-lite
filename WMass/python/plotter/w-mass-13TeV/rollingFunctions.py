from array import array

## safe batch mode
import sys
args = sys.argv[:]
sys.argv = ['-b']
import ROOT
sys.argv = args
ROOT.gROOT.SetBatch(True)
ROOT.PyConfig.IgnoreCommandLineOptions = True

def roll1Dto2D(h1d, histo, invertXY=False):
    # need to know whether the unrolled histograms have consecutive pt shapes (usually y axis of 2D plots), in which case invertXY must be True, or eta shapes (usually x axis of 2D plots)
    # we used to have eta shapes in unrolled TH1D passed to combinetf, but now we feed combinetf directly with 2D histograms, and it internally unrolls them in the other axis
    nBinsXaxis2D = histo.GetNbinsY() if invertXY else histo.GetNbinsX()
    for i in range(1,h1d.GetNbinsX()+1):
        # histogram bin is numbered starting from 1, so add 1
        xbin = (i - 1) % nBinsXaxis2D + 1  
        ybin = (i - 1) / nBinsXaxis2D + 1
        xbin = int(xbin)
        ybin = int(ybin)
        if invertXY:
            # our 2D plots will always have eta on x axis, so we must also swap xbin and ybin defined above
            xbin,ybin = ybin,xbin
        histo.SetBinContent(xbin, ybin, h1d.GetBinContent(i))
        histo.SetBinError(xbin, ybin, h1d.GetBinError(i)) 

    return histo

def dressed2D(h1d, binning, name, title='', invertXY=False):
    if len(binning) == 4:
        n1 = binning[0]; bins1 = array('d', binning[1])
        n2 = binning[2]; bins2 = array('d', binning[3])
        h2_1 = ROOT.TH2D(name, title, n1, bins1, n2, bins2 )
    else:
        n1 = binning[0]; min1 = binning[1]; max1 = binning[2]
        n2 = binning[3]; min2 = binning[4]; max2 = binning[5]
        h2_1 = ROOT.TH2D(name, title, n1, min1, max1, n2, min2, max2)
    h2_backrolled_1 = roll1Dto2D(h1d, h2_1, invertXY)
    return h2_backrolled_1


def unroll2Dto1D(h, newname='', cropNegativeBins=True, silent=False):
    nbins = h.GetNbinsX() * h.GetNbinsY()
    goodname = h.GetName()
    h.SetName(goodname+"_2d")
    newh = ROOT.TH1D(goodname if not newname else newname,h.GetTitle(),nbins,0.5,nbins+0.5)
    newh.Sumw2()
    if 'TH2' not in h.ClassName(): raise RuntimeError("Calling rebin2Dto1D on something that is not TH2")
    for i in range(h.GetNbinsX()):
        for j in range(h.GetNbinsY()):
            bin = 1 + i + j*h.GetNbinsX()
            newh.SetBinContent(bin,h.GetBinContent(i+1,j+1))
            newh.SetBinError(bin,h.GetBinError(i+1,j+1))
    if cropNegativeBins:
        for bin in range(1,nbins+1):
            if newh.GetBinContent(bin)<0:
                if not silent:
                    print('Warning: cropping to zero bin %d in %s (was %f)'%(bin,newh.GetName(),newh.GetBinContent(bin)))
                newh.SetBinContent(bin,0)
    newh.SetLineWidth(h.GetLineWidth())
    newh.SetLineStyle(h.GetLineStyle())
    newh.SetLineColor(h.GetLineColor())
    return newh

def reverseUnroll2Dto1D(h, newname='', cropNegativeBins=True):
    h2swap = ROOT.TH2D(h.GetName()+'_swapped','',
                       h.GetNbinsY(),h.GetYaxis().GetBinLowEdge(1),h.GetYaxis().GetBinUpEdge(h.GetNbinsY()),
                       h.GetNbinsX(),h.GetXaxis().GetBinLowEdge(1),h.GetXaxis().GetBinUpEdge(h.GetNbinsX()))
    for ieta in range(1,h.GetNbinsX()+1):
        for ipt in range(1,h.GetNbinsY()+1):
            h2swap.SetBinContent(ipt,ieta,h.GetBinContent(ieta,ipt))
            h2swap.SetBinError(ipt,ieta,h.GetBinError(ieta,ipt))
    return unroll2Dto1D(h2swap, newname, cropNegativeBins)
