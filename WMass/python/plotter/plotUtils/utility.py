#!/usr/bin/env python

from shutil import copyfile
import re, sys, os, os.path, subprocess, json, ROOT, copy, math
import numpy as np
from array import array

from CMS_lumi import *

_canvas_pull = ROOT.TCanvas("_canvas_pull","",800,800)
    

#########################################################################

def addStringToEnd(name, matchToAdd, notAddIfEndswithMatch=False):
    if notAddIfEndswithMatch and name.endswith(matchToAdd):
        return name
    elif not name.endswith(matchToAdd):
        return name + matchToAdd

#########################################################################

def getZaxisReasonableExtremesTH2(h,nSigma=3,minZtoUse=None,maxZtoUse=None):

    htmp = ROOT.TH1D("htmp","",1000,h.GetMinimum(),h.GetMaximum())
    nbins = h.GetNbinsX() * h.GetNbinsY()    
    for ibin in range (1,nbins+1):
        val = h.GetBinContent(ibin)
        canFill = True
        if minZtoUse != None:
            if val < minZtoUse: canFill = False
        if maxZtoUse != None:
            if val > maxZtoUse: canFill = False
        if canFill: htmp.Fill(val)

    mean = htmp.GetMean()
    stddev = htmp.GetStdDev()
    retmin = max(h.GetMinimum(),mean - nSigma*stddev)
    retmax = min(h.GetMaximum(),mean + nSigma*stddev)
    return retmin,retmax


#########################################################################

def getMinMaxHisto(h, excludeEmpty=True, sumError=True, 
                   excludeUnderflow=True, excludeOverflow=True,
                   excludeMin=None, excludeMax=None):
    
    # Warning, fix this function, GetBinContent with TH2 is not that simple, there are the underflow and overflow in each row and column
    # must check whether bin is underflow or overflow
    # therefore, the global bin is obtained as the number of bins +2, multiplied for each axis

    # excludeEmpty = True exclude bins with content 0.0. Useful when a histogram is filled with values in, for example, [1,2] but hassome empty bins
    # excludeMin/Max are used to select a range in which to look for maximum and minimum, useful to reject outliers, crazy or empty bins and so on
    # for histograms with non-negative values, excludeEmpty=True is equal to excludeMin==0.0

    # sumError is used to add or subtract error when looking for min/max (to have full error band in range)
    # when using excludeMin/Max, the errors are still ignored when evaluating the range

    # the better combination of options depends on dimension: for a TH1 is useful to visualize the error band in the plot range, while for a TH2 
    # only the bin content is interesting in the plot (the error is not reported with TH2::Draw, unless plotting it in a 3D space

    # one might exploit excludeMin/Max to select a rage depending on the distribution on the histogram bin content
    # for example, one can pass excludeMin=h.GetMean()-2*h.GetStdDev() and excludeMax=h.GetMean()+2*h.GetStdDev() so to 
    # select a range of 2 sigma around the mean

    dim = h.GetDimension()
    nbins = 0
    if   dim == 1: nbins = h.GetNbinsX() + 2
    elif dim == 2: nbins = (h.GetNbinsX() + 2) * (h.GetNbinsY() + 2)
    elif dim == 3: nbins = (h.GetNbinsX() + 2) * (h.GetNbinsY() + 2) * (h.GetNbinsZ() + 2)
    else:
        print "Error in getMaxHisto(): dim = %d is not supported. Exit" % dim
        quit()

    maxval = 0
    minval = 0
    firstValidBin = -1
    for ibin in range (1,nbins+1):
        if excludeUnderflow and h.IsBinUnderflow(ibin): continue
        if excludeOverflow and h.IsBinOverflow(ibin): continue
        tmpmax = h.GetBinContent(ibin)
        tmpmin = h.GetBinContent(ibin)
        if excludeEmpty and tmpmin == 0.0: continue
        if excludeMin != None and tmpmin <= excludeMin: continue
        if excludeMax != None and tmpmax >= excludeMax: continue
        if firstValidBin < 0: 
            #print "ibin %d:   tmpmin,tmpmax = %.2f, %.2f" % (ibin,tmpmin,tmpmax)
            firstValidBin = ibin
        if sumError:
            tmpmin -= h.GetBinError(ibin)
            tmpmax += h.GetBinError(ibin)
        if firstValidBin > 0 and ibin == firstValidBin:
            #the first time we pick a non empty bin, we set min and max to the histogram content in that bin
            minval = tmpmin
            maxval = tmpmax
            #print "#### ibin %d:   min,max = %.2f, %.2f" % (ibin,minval,maxval)
        else:
            minval = min(minval,tmpmin)
            maxval = max(maxval,tmpmax)
        #print "ibin %d:   min,max = %.2f, %.2f" % (ibin,minval,maxval)
    
    return minval,maxval

#########################################################################

def getMinimumTH(h, excludeMin=None):
    # get minimum excluding some values. For example, if an histogram has an empty bin, one might want to get the minimum such that it is > 0
    # underflow are not considered
    
    dim = h.GetDimension()
    retmin = sys.float_info.max

    if dim == 1:
        for ix in range(1,h.GetNbinsX()+1):
            if retmin > h.GetBinContent(ix):
                if excludeMin != None:
                    if h.GetBinContent(ix) > excludeMin: retmin = h.GetBinContent(ix)
                else:
                    retmin = h.GetBinContent(ix)

    elif dim == 2:
        for ix in range(1,h.GetNbinsX()+1):
            for iy in range(1,h.GetNbinsY()+1):
                if retmin > h.GetBinContent(ix,iy):
                    if excludeMin != None:
                        if h.GetBinContent(ix,iy) > excludeMin: retmin = h.GetBinContent(ix,iy)
                    else:
                        retmin = h.GetBinContent(ix,iy)

    elif dim == 3:
        for ix in range(1,h.GetNbinsX()+1):
            for iy in range(1,h.GetNbinsY()+1):
                for iz in range(1,h.GetNbinsZ()+1):
                    if retmin > h.GetBinContent(ix,iy,iz):
                        if excludeMin != None:
                            if h.GetBinContent(ix,iy,iz) > excludeMin: retmin = h.GetBinContent(ix,iy,iz)
                        else:
                            retmin = h.GetBinContent(ix,iy,iz)
                            

    else:
        raise RuntimeError, "Error in getMinimumTH(): unsupported histogram's dimension (%d)" % dim

    return retmin

#########################################################################

def getMaximumTH(h, excludeMax=None):
    # get maximum excluding some values. For example, if an histogram has a crazy bin, one might want to get the maximum value that is lower than that
    # overflow are not considered
    
    dim = h.GetDimension()
    retmax = sys.float_info.min

    if dim == 1:
        for ix in range(1,h.GetNbinsX()+1):
            if retmax < h.GetBinContent(ix):
                if excludeMax != None:
                    if h.GetBinContent(ix) < excludeMax: retmax = h.GetBinContent(ix)
                else:
                    retmax = h.GetBinContent(ix)

    elif dim == 2:
        for ix in range(1,h.GetNbinsX()+1):
            for iy in range(1,h.GetNbinsY()+1):
                if retmax < h.GetBinContent(ix,iy):
                    if excludeMax != None:
                        if h.GetBinContent(ix,iy) < excludeMax: retmax = h.GetBinContent(ix,iy)                        
                    else:
                        retmax = h.GetBinContent(ix,iy)

    elif dim == 3:
        for ix in range(1,h.GetNbinsX()+1):
            for iy in range(1,h.GetNbinsY()+1):
                for iz in range(1,h.GetNbinsZ()+1):
                    if retmax < h.GetBinContent(ix,iy,iz):
                        if excludeMax != None:
                            if h.GetBinContent(ix,iy,iz) < excludeMax: retmax = h.GetBinContent(ix,iy,iz)                            
                        else:
                            retmax = h.GetBinContent(ix,iy,iz)

    else:
        raise RuntimeError, "Error in getMaximumTH(): unsupported histogram's dimension (%d)" % dim

    return retmax


#########################################################################


def createPlotDirAndCopyPhp(outdir):
    if outdir != "./":
        if not os.path.exists(outdir):
            os.system("mkdir -p "+outdir)
            if os.path.exists("/afs/cern.ch"): os.system("cp /afs/cern.ch/user/e/emanuele/public/index.php "+outdir)
    

#########################################################################

def getAxisRangeFromUser(axisNameTmp="", 
                         separator="::", 
                         rangeSeparator=","
                         ):
  
    setXAxisRangeFromUser = False;
    fields = axisNameTmp.split(separator)
    axisName = fields[0]
    
    if len(fields) > 1:
        setXAxisRangeFromUser = True;
        xmin = float(fields[1].split(rangeSeparator)[0])
        xmax = float(fields[1].split(rangeSeparator)[1])
    else:
        xmin = 0
        xmax = 0
        
    return axisName,setXAxisRangeFromUser,xmin,xmax


#########################################################################

def adjustSettings_CMS_lumi():

    ## dummy function to be called before using any other fucntion calling CMS_lumi
    ## for some reason, the settings of the very first plot are screwed up.
    ## To fix this issue, it is enough to call it to a dummy plot
    dummy = ROOT.TH1D("dummy","",10,0,10)
    for i in range(1,1+dummy.GetNbinsX()):
        dummy.SetBinContent(i,i)
    dummy.GetXaxis().SetTitle("x axis")
    dummy.GetYaxis().SetTitle("y axis")
    cdummy = ROOT.TCanvas("cdummy","",600,600)
    dummy.Draw("HE")
    CMS_lumi(cdummy,"",True,False)
    setTDRStyle()        
    ## no need to save the canvas    


#########################################################################

def drawTH1(htmp,
            labelXtmp="xaxis",
            labelYtmp="Events",
            outdir= "./",
            prefix = "distribution",
            outhistname = "",
            canvasSize="700,625",
            passCanvas=None,
            moreTextLatex=""
            ):


    addStringToEnd(outdir,"/",notAddIfEndswithMatch=True)
    createPlotDirAndCopyPhp(outdir)

    labelX,setXAxisRangeFromUser,xmin,xmax = getAxisRangeFromUser(labelXtmp)
    labelY,setYAxisRangeFromUser,ymin,ymax = getAxisRangeFromUser(labelYtmp)

    h = htmp.Clone("htmp")

    cw,ch = canvasSize.split(',')
    canvas = passCanvas if passCanvas != None else ROOT.TCanvas("canvas","",int(cw),int(ch))
    canvas.SetTickx(1)
    canvas.SetTicky(1)
    canvas.cd()
    canvas.SetLeftMargin(0.12)
    canvas.SetRightMargin(0.04)
    canvas.cd()

    h.SetLineColor(ROOT.kBlack)
    h.SetLineWidth(2)
    h.GetXaxis().SetTitle(labelX)
    h.GetXaxis().SetTitleOffset(1.2)
    h.GetXaxis().SetTitleSize(0.05)
    h.GetXaxis().SetLabelSize(0.04)
    h.GetYaxis().SetTitle(labelY)
    h.GetYaxis().SetTitleOffset(1.15)
    h.GetYaxis().SetTitleSize(0.05)
    h.GetYaxis().SetLabelSize(0.04)
    if (setXAxisRangeFromUser): h.GetXaxis().SetRangeUser(xmin,xmax)
    if (setYAxisRangeFromUser): h.GetYaxis().SetRangeUser(ymin,ymax)
    # force drawing stat box
    h.SetStats(1)
    h.Draw("HIST")
    canvas.RedrawAxis("sameaxis")
    setTDRStyle()
    # force drawing stat box
    ROOT.gStyle.SetOptStat(111110)
    ROOT.gStyle.SetOptFit(1102)
    #
    if len(moreTextLatex):
        realtext = moreTextLatex.split("::")[0]
        x1,y1,ypass,textsize = 0.75,0.8,0.08,0.035
        if "::" in moreTextLatex:
            x1,y1,ypass,textsize = (float(x) for x in (moreTextLatex.split("::")[1]).split(","))            
        lat = ROOT.TLatex()
        lat.SetNDC();
        lat.SetTextFont(42)        
        lat.SetTextSize(textsize)
        for itx,tx in enumerate(realtext.split(";")):
            lat.DrawLatex(x1,y1-itx*ypass,tx)

    for ext in ["png","pdf"]:
        canvas.SaveAs(outdir + "{pfx}{hname}.{ext}".format(pfx=prefix,hname=("_"+outhistname) if len(outhistname) else "",ext=ext))



#########################################################################




# function to draw 2D histograms, can also plot profile along X on top
def drawCorrelationPlot(h2D_tmp,
                        labelXtmp="xaxis", labelYtmp="yaxis", labelZtmp="zaxis",
                        canvasName="default", plotLabel="", outdir="./",
                        rebinFactorY=0,
                        rebinFactorX=0,
                        smoothPlot=True,
                        drawProfileX=True,
                        scaleToUnitArea=True,
                        draw_both0_noLog1_onlyLog2=0,
                        leftMargin=0.16,
                        rightMargin=0.20,
                        nContours=51,
                        palette=55,
                        canvasSize="700,625",
                        passCanvas=None,
                        bottomMargin=0.1,
                        plotError=False,
                        lumi=None,
                        drawOption = "colz",
                        writeSimulation=False):


    # if h2D.GetName() == "scaleFactor_origBinPt":
    #     print "="*20
    #     print "Check: hist %s: Z axis title = %s" % (h2D.GetName(),labelZtmp)
    #     print "="*20

    ROOT.TH1.SetDefaultSumw2()
    adjustSettings_CMS_lumi()

    if (rebinFactorX): 
        if isinstance(rebinFactorX, int): h2D_tmp.RebinY(rebinFactorX)
        else:                             h2D_tmp.RebinY(len(rebinFactorX)-1,"",array('d',rebinFactorX)) # case in which rebinFactorX is a list of bin edges

    if (rebinFactorY): 
        if isinstance(rebinFactorY, int): h2D_tmp.RebinY(rebinFactorY)
        else:                             h2D_tmp.RebinY(len(rebinFactorY)-1,"",array('d',rebinFactorY)) # case in which rebinFactorX is a list of bin edges

    if plotError:
        herr = h2D_tmp.Clone(h2D_tmp.GetName()+"_err")
        herr.Reset("ICESM")
        for i in range(1,herr.GetNbinsX()+1):
            for j in range(1,herr.GetNbinsY()+1):
                herr.SetBinContent(i,j,h2D_tmp.GetBinError(i,j))
        h2D = herr
    else:
        h2D = h2D_tmp

    ROOT.TColor.CreateGradientColorTable(3,
                                         array ("d", [0.00, 0.50, 1.00]),
                                         ##array ("d", [1.00, 1.00, 0.00]),        
                                         ##array ("d", [0.70, 1.00, 0.34]),        
                                         ##array ("d", [0.00, 1.00, 0.82]),        
                                         array ("d", [0.00, 1.00, 1.00]),
                                         array ("d", [0.34, 1.00, 0.65]),
                                         array ("d", [0.82, 1.00, 0.00]),
                                         255,  0.95)

    if palette > 0: ROOT.gStyle.SetPalette(palette)  # 55:raibow palette ; 57: kBird (blue to yellow, default) ; 107 kVisibleSpectrum ; 77 kDarkRainBow 
    ROOT.gStyle.SetNumberContours(nContours) # default is 20 

    labelX,setXAxisRangeFromUser,xmin,xmax = getAxisRangeFromUser(labelXtmp)
    labelY,setYAxisRangeFromUser,ymin,ymax = getAxisRangeFromUser(labelYtmp)
    labelZ,setZAxisRangeFromUser,zmin,zmax = getAxisRangeFromUser(labelZtmp)
    
    cw,ch = canvasSize.split(',')
    #canvas = ROOT.TCanvas("canvas",h2D.GetTitle() if plotLabel == "ForceTitle" else "",700,625)    
    canvas = passCanvas if passCanvas != None else ROOT.TCanvas("canvas","",int(cw),int(ch))
    canvas.SetTickx(1)
    canvas.SetTicky(1)
    canvas.SetLeftMargin(leftMargin)
    canvas.SetRightMargin(rightMargin)
    canvas.SetBottomMargin(bottomMargin)
    canvas.cd()

    addStringToEnd(outdir,"/",notAddIfEndswithMatch=True)
    createPlotDirAndCopyPhp(outdir)
    # normalize to 1
    if (scaleToUnitArea): h2D.Scale(1./h2D.Integral())

    h2DGraph = 0

    h2DPlot = 0
    if (not smoothPlot): h2DPlot = h2D
    else:
        h2DGraph = ROOT.TGraph2D()
        h2DGraph.SetNpx(300)
        h2DGraph.SetNpy(300)
        nPoint = 0
        for iBinX in range (1,1+h2D.GetNbinsX()):
            for iBinY in range(1,1+h2D.GetNbinsY()):
                h2DGraph.SetPoint(nPoint,h2D.GetXaxis().GetBinCenter(iBinX),h2D.GetYaxis().GetBinCenter(iBinY),h2D.GetBinContent(iBinX,iBinY))
                nPoint += 1
            

        h2DPlot = h2DGraph.GetHistogram()

    if plotLabel == "ForceTitle":
        h2DPlot.SetTitle(h2D_tmp.GetTitle())
  
    h2DPlot.GetXaxis().SetTitle(labelX)
    h2DPlot.GetYaxis().SetTitle(labelY)
    h2DPlot.GetXaxis().SetTitleSize(0.05)
    h2DPlot.GetXaxis().SetLabelSize(0.04)
    h2DPlot.GetXaxis().SetTitleOffset(0.95) # 1.1 goes outside sometimes, maybe depends on root version or canvas width
    h2DPlot.GetYaxis().SetTitleSize(0.05)
    h2DPlot.GetYaxis().SetLabelSize(0.04)
    h2DPlot.GetYaxis().SetTitleOffset(1.1)
    h2DPlot.GetZaxis().SetTitleSize(0.05)
    h2DPlot.GetZaxis().SetLabelSize(0.04)
    h2DPlot.GetZaxis().SetTitleOffset(1.2)

    h2DPlot.GetZaxis().SetTitle(labelZ) 
    h2DPlot.Draw(drawOption)

    if (setXAxisRangeFromUser): h2DPlot.GetXaxis().SetRangeUser(xmin,xmax)
    if (setYAxisRangeFromUser): h2DPlot.GetYaxis().SetRangeUser(ymin,ymax)
    if (setZAxisRangeFromUser): h2DPlot.GetZaxis().SetRangeUser(zmin,zmax)


    # if h2D.GetName() == "scaleFactor_origBinPt":
    #     print "="*20
    #     print "Check: hist %s: Z axis title = %s" % (h2DPlot.GetName(),h2DPlot.GetZaxis().GetTitle())
    #     print "="*20

    # # attempt to make Z axis title farther depending on how many digits are printed
    # maxZaxisVal = h2DPlot.GetBinContent(h2DPlot.GetMaximumBin())
    # if (setZAxisRangeFromUser): maxZaxisVal = zmax

    # if maxZaxisVal >= 1.0:
    #     rootYear = int(str(ROOT.gROOT.GetVersionDate())[:4])        
    #     if (rootYear > 2016):
    #         h2DPlot.GetZaxis().SetMaxDigits(3)
    #     else:
    #         print "Warning in drawCorrelationPlot: TAxis::SetMaxDigits() not implemented for ROOT versions before 2017 (rough estimate)"
    #         print "Will not exit, but instruction will be neglected"
    #     if maxZaxisVal > 9999.:
    #         h2DPlot.GetZaxis().SetTitleOffset(h2DPlot.GetZaxis().GetTitleOffset()+0.15)
    #         print "Changing title offset by 0.15"
    # else:
    #     i = 1
    #     tryNext = True
    #     while (tryNext and i < 6):
    #         tmpVal = maxZaxisVal * pow(10,i)
    #         if tmpVal >= 1.0: tryNext = False 
    #         else: i += 1
    #     if i > 1:            
    #         print "Max Z axis < 1, will try to adjust distance of Z axis title to Z axis"
    #         print "i = %d: will move Z axis offset by 0.45" % i
    #         # for numbers like 0.025 or with more 0 after ., make increase distance between Z axis title and the Z axis
    #         h2DPlot.GetZaxis().SetTitleOffset(h2DPlot.GetZaxis().GetTitleOffset()+0.45)

    h2DPlot.GetZaxis().SetTitleOffset(h2DPlot.GetZaxis().GetTitleOffset()+0.4)


    h2DProfile = 0
    if drawProfileX:
        h2DProfile = h2D.ProfileX("%s_pfx" %h2D.GetName())
        h2DProfile.SetMarkerColor(ROOT.kBlack)
        h2DProfile.SetMarkerStyle(20)
        h2DProfile.SetMarkerSize(1)
        h2DProfile.Draw("EPsame")
        
    # not yet implemented
    setTDRStyle()
    if not plotLabel == "ForceTitle": 
        if lumi != None: CMS_lumi(canvas,lumi,True,False,isSimulation=writeSimulation)
        else:            CMS_lumi(canvas,"",True,False,isSimulation=writeSimulation)
    #setTDRStyle()
    #print ">>>>>>>>>>>>>> check <<<<<<<<<<<<<<<<<<<"

    if plotLabel == "ForceTitle":
        ROOT.gStyle.SetOptTitle(1)        

    #h2DPlot.GetZaxis().SetMaxDigits(1)  #for N>99, should use scientific notation, I'd like to make it work only with negative exponential but haven't succeeded yet
    # canvas.Modified()
    # canvas.Update()

    leg = ROOT.TLegend(0.39,0.75,0.89,0.95)
    leg.SetFillStyle(0)
    leg.SetFillColor(0)
    leg.SetBorderSize(0)
    leg.SetTextFont(62)
    if plotLabel not in ["", "ForceTitle"]: leg.AddEntry(0,plotLabel,"")
    if drawProfileX: leg.AddEntry(0,"Correlation = %.2f" % h2DPlot.GetCorrelationFactor(),"")
    leg.Draw("same")

    if (draw_both0_noLog1_onlyLog2 == 0 or draw_both0_noLog1_onlyLog2 == 1):
        for ext in ['png', 'pdf']:
            canvas.SaveAs('{od}/{cn}.{ext}'.format(od=outdir, cn=canvasName, ext=ext))
        
    if (draw_both0_noLog1_onlyLog2 == 0 or draw_both0_noLog1_onlyLog2 == 2):
        canvas.SetLogz()
        for ext in ['png', 'pdf']:
            canvas.SaveAs('{od}/{cn}_logZ.{ext}'.format(od=outdir, cn=canvasName, ext=ext))
        canvas.SetLogz(0)
        
    return canvas

##########################################################


def drawSingleTH1(h1,
                  labelXtmp="xaxis", labelYtmp="yaxis",
                  canvasName="default", outdir="./",
                  rebinFactorX=0,
                  draw_both0_noLog1_onlyLog2=0,                  
                  leftMargin=0.15,
                  rightMargin=0.04,
                  labelRatioTmp="Rel.Unc.::0.5,1.5",
                  drawStatBox=False,
                  legendCoords="0.15,0.35,0.8,0.9",  # x1,x2,y1,y2
                  canvasSize="600,700",  # use X,Y to pass X and Y size     
                  lowerPanelHeight = 0.3,  # number from 0 to 1, 0.3 means 30% of space taken by lower panel. 0 means do not draw lower panel with relative error
                  drawLineLowerPanel="luminosity uncertainty::0.025", # if not empty, draw band at 1+ number after ::, and add legend with title
                  passCanvas=None,
                  lumi=None,
                  drawVertLines="", # "12,36": format --> N of sections (e.g: 12 pt bins), and N of bins in each section (e.g. 36 eta bins), assuming uniform bin width
                  textForLines=[],                       
                  moreText="",
                  moreTextLatex=""
                  ):

    # moreText is used to pass some text to write somewhere (TPaveText is used)
    # e.g.  "stuff::x1,y1,x2,y2"  where xi and yi are the coordinates for the text
    # one can add more lines using the ";" key. FOr example, "stuff1;stuff2::x1,y1,x2,y2"
    # the coordinates should be defined taking into account how many lines will be drawn
    # if the coordinates are not passed (no "::"), then default ones are used, but this might not be satisfactory

    # moreTextLatex is similar, but used TLatex, and the four coordinates are x1,y1,ypass,textsize
    # where x1 and y1 are the coordinates the first line, and ypass is how much below y1 the second line is (and so on for following lines)

    if (rebinFactorX): 
        if isinstance(rebinFactorX, int): h1.Rebin(rebinFactorX)
        # case in which rebinFactorX is a list of bin edges
        else:                             h1.Rebin(len(rebinFactorX)-1,"",array('d',rebinFactorX)) 

    xAxisName,setXAxisRangeFromUser,xmin,xmax = getAxisRangeFromUser(labelXtmp)
    yAxisName,setYAxisRangeFromUser,ymin,ymax = getAxisRangeFromUser(labelYtmp)
    yRatioAxisName,setRatioYAxisRangeFromUser,yminRatio,ymaxRatio = getAxisRangeFromUser(labelRatioTmp)

    yAxisTitleOffset = 1.45 if leftMargin > 0.1 else 0.6

    addStringToEnd(outdir,"/",notAddIfEndswithMatch=True)
    createPlotDirAndCopyPhp(outdir)

    cw,ch = canvasSize.split(',')
    #canvas = ROOT.TCanvas("canvas",h2D.GetTitle() if plotLabel == "ForceTitle" else "",700,625)
    canvas = passCanvas if passCanvas != None else ROOT.TCanvas("canvas","",int(cw),int(ch))
    canvas.SetTickx(1)
    canvas.SetTicky(1)
    canvas.cd()
    canvas.SetLeftMargin(leftMargin)
    canvas.SetRightMargin(rightMargin)
    canvas.cd()

    pad2 = 0
    if lowerPanelHeight: 
        canvas.SetBottomMargin(lowerPanelHeight)
        pad2 = ROOT.TPad("pad2","pad2",0,0.,1,0.9)
        pad2.SetTopMargin(1-lowerPanelHeight)
        pad2.SetRightMargin(rightMargin)
        pad2.SetLeftMargin(leftMargin)
        pad2.SetFillColor(0)
        pad2.SetGridy(1)
        pad2.SetFillStyle(0)


    frame = h1.Clone("frame")
    frame.GetXaxis().SetLabelSize(0.04)
    frame.SetStats(0)

    h1.SetLineColor(ROOT.kBlack)
    h1.SetMarkerColor(ROOT.kBlack)
    h1.SetMarkerStyle(20)
    h1.SetMarkerSize(0)

    #ymax = max(ymax, max(h1.GetBinContent(i)+h1.GetBinError(i) for i in range(1,h1.GetNbinsX()+1)))
    # if min and max were not set, set them based on histogram content
    if ymin == ymax == 0.0:
        ymin,ymax = getMinMaxHisto(h1,excludeEmpty=True,sumError=True)            
        ymin *= 0.9
        ymax *= (1.1 if leftMargin > 0.1 else 2.0)
        if ymin < 0: ymin = 0
        #print "drawSingleTH1() >>> Histo: %s     minY,maxY = %.2f, %.2f" % (h1.GetName(),ymin,ymax)

    # print "#### WARNING ####"
    # print "Hardcoding ymin = 0 in function drawSingleTH1(): change it if it is not what you need"
    # print "#################"
    # ymin = 0 # hardcoded

    if lowerPanelHeight:
        h1.GetXaxis().SetLabelSize(0)
        h1.GetXaxis().SetTitle("")  
    else:
        h1.GetXaxis().SetTitle(xAxisName)
        h1.GetXaxis().SetTitleOffset(1.2)
        h1.GetXaxis().SetTitleSize(0.05)
        h1.GetXaxis().SetLabelSize(0.04)
    h1.GetYaxis().SetTitle(yAxisName)
    h1.GetYaxis().SetTitleOffset(yAxisTitleOffset) 
    h1.GetYaxis().SetTitleSize(0.05)
    h1.GetYaxis().SetLabelSize(0.04)
    h1.GetYaxis().SetRangeUser(ymin, ymax)
    h1.GetYaxis().SetTickSize(0.01)
    if setXAxisRangeFromUser: h1.GetXaxis().SetRangeUser(xmin,xmax)
    h1.Draw("HIST")
    h1err = h1.Clone("h1err")
    h1err.SetFillColor(ROOT.kRed+2)
    h1err.SetFillStyle(3001)  # 3001 is better than 3002 for pdf, while 3002 is perfect for png
    #h1err.SetFillStyle(3002)
    #h1err.SetFillStyle(3005)
    h1err.Draw("E2same")
    #h1.Draw("HIST same")

    nColumnsLeg = 1
    if ";" in legendCoords: 
        nColumnsLeg = int(legendCoords.split(";")[1])
    legcoords = [float(x) for x in (legendCoords.split(";")[0]).split(',')]
    lx1,lx2,ly1,ly2 = legcoords[0],legcoords[1],legcoords[2],legcoords[3]
    leg = ROOT.TLegend(lx1,ly1,lx2,ly2)
    #leg.SetFillColor(0)
    #leg.SetFillStyle(0)
    #leg.SetBorderSize(0)
    leg.SetNColumns(nColumnsLeg)
    leg.AddEntry(h1,"Value","L")
    leg.AddEntry(h1err,"Uncertainty","F")
    leg.Draw("same")
    canvas.RedrawAxis("sameaxis")

    if drawStatBox:
        ROOT.gPad.Update()
        ROOT.gStyle.SetOptStat(1110)
        ROOT.gStyle.SetOptFit(1102)
    else:
        h1.SetStats(0)

    vertline = ROOT.TLine(36,0,36,canvas.GetUymax())
    vertline.SetLineColor(ROOT.kBlack)
    vertline.SetLineStyle(3)
    bintext = ROOT.TLatex()
    #bintext.SetNDC()
    bintext.SetTextSize(0.025)  # 0.03
    bintext.SetTextFont(42)
    if len(textForLines): bintext.SetTextAngle(45 if "#eta" in textForLines[0] else 30)

    if len(drawVertLines):
        #print "drawVertLines = " + drawVertLines
        nptBins = int(drawVertLines.split(',')[0])
        etarange = float(drawVertLines.split(',')[1])        
        offsetXaxisHist = h1.GetXaxis().GetBinLowEdge(0)
        sliceLabelOffset = 6. 
        for i in range(1,nptBins): # do not need line at canvas borders
            #vertline.DrawLine(offsetXaxisHist+etarange*i,0,offsetXaxisHist+etarange*i,canvas.GetUymax())
            vertline.DrawLine(etarange*i-offsetXaxisHist,0,etarange*i-offsetXaxisHist,ymax)
        if len(textForLines):
            for i in range(0,len(textForLines)): # we need nptBins texts
                #texoffset = 0.1 * (4 - (i%4))
                #ytext = (1. + texoffset)*ymax/2.  
                ytext = (1.1)*ymax/2.  
                bintext.DrawLatex(etarange*i + etarange/sliceLabelOffset, ytext, textForLines[i])

    # redraw legend, or vertical lines appear on top of it
    leg.Draw("same")

    if len(moreText):
        realtext = moreText.split("::")[0]
        x1,y1,x2,y2 = 0.7,0.8,0.9,0.9
        if "::" in moreText:
            x1,y1,x2,y2 = (float(x) for x in (moreText.split("::")[1]).split(","))
        pavetext = ROOT.TPaveText(x1,y1,x2,y2,"NB NDC")
        for tx in realtext.split(";"):
            pavetext.AddText(tx)
        pavetext.SetFillColor(0)
        pavetext.SetFillStyle(0)
        pavetext.SetBorderSize(0)
        pavetext.SetLineColor(0)
        pavetext.Draw("same")

    if len(moreTextLatex):
        realtext = moreTextLatex.split("::")[0]
        x1,y1,ypass,textsize = 0.75,0.8,0.08,0.035
        if "::" in moreTextLatex:
            x1,y1,ypass,textsize = (float(x) for x in (moreTextLatex.split("::")[1]).split(","))            
        lat = ROOT.TLatex()
        lat.SetNDC();
        lat.SetTextFont(42)        
        lat.SetTextSize(textsize)
        for itx,tx in enumerate(realtext.split(";")):
            lat.DrawLatex(x1,y1-itx*ypass,tx)


  # TPaveText *pvtxt = NULL;
  # if (yAxisName == "a.u.") {
  #   pvtxt = new TPaveText(0.6,0.6,0.95,0.7, "BR NDC")
  #   pvtxt.SetFillColor(0)
  #   pvtxt.SetFillStyle(0)
  #   pvtxt.SetBorderSize(0)
  #   pvtxt.AddText(Form("norm num/den = %.2f +/- %.2f",IntegralRatio,ratioError))
  #   pvtxt.Draw()
  # }

    setTDRStyle()
    if leftMargin > 0.1:
        if lumi != None: CMS_lumi(canvas,lumi,True,False)
        else:            CMS_lumi(canvas,"",True,False)
    else:
        latCMS = ROOT.TLatex()
        latCMS.SetNDC();
        latCMS.SetTextFont(42)
        latCMS.SetTextSize(0.045)
        latCMS.DrawLatex(0.1, 0.95, '#bf{CMS} #it{Preliminary}')
        if lumi != None: latCMS.DrawLatex(0.85, 0.95, '%s fb^{-1} (13 TeV)' % lumi)
        else:            latCMS.DrawLatex(0.90, 0.95, '(13 TeV)' % lumi)

    if lowerPanelHeight:
        pad2.Draw()
        pad2.cd()

        frame.Reset("ICES")
        if setRatioYAxisRangeFromUser: frame.GetYaxis().SetRangeUser(yminRatio,ymaxRatio)
        #else:                          
        #frame.GetYaxis().SetRangeUser(0.5,1.5)
        frame.GetYaxis().SetNdivisions(5)
        frame.GetYaxis().SetTitle(yRatioAxisName)
        frame.GetYaxis().SetTitleOffset(yAxisTitleOffset)
        frame.GetYaxis().SetTitleSize(0.05)
        frame.GetYaxis().SetLabelSize(0.04)
        frame.GetYaxis().CenterTitle()
        frame.GetXaxis().SetTitle(xAxisName)
        if setXAxisRangeFromUser: frame.GetXaxis().SetRangeUser(xmin,xmax)
        frame.GetXaxis().SetTitleOffset(1.2)
        frame.GetXaxis().SetTitleSize(0.05)

        ratio = h1.Clone("ratio")
        den_noerr = h1.Clone("den_noerr")
        for iBin in range (1,den_noerr.GetNbinsX()+1):
            den_noerr.SetBinError(iBin,0.)

        ratio.Divide(den_noerr)
        ratio.SetFillColor(ROOT.kGray+1)
        #den_noerr.SetFillColor(ROOT.kGray)
        frame.Draw()
        ratio.SetMarkerSize(0)
        ratio.SetMarkerStyle(0) # important to remove dots at y = 1
        ratio.Draw("E2same")

        # print values
        # if "unrolledXsec_eta_abs_" in canvasName or "xsec_eta_abs_" in canvasName:
        #     print "eta-bin   Rel.Unc (%)"
        #     for i in range (1,ratio.GetNbinsX()+1):
        #         if "unrolledXsec_eta_abs_" in canvasName:
        #             if i%etarange == 1:
        #                 print textForLines[int((i-1)/etarange)]
        #             print "%s  %s" % (str(int(i%etarange)),str(100*ratio.GetBinError(i)))
        #         else:
        #             print "%s  %s" % (str(i),str(100*ratio.GetBinError(i)))
   
        line = ROOT.TF1("horiz_line","1",ratio.GetXaxis().GetBinLowEdge(1),ratio.GetXaxis().GetBinLowEdge(ratio.GetNbinsX()+1))
        line.SetLineColor(ROOT.kRed)
        line.SetLineWidth(1)
        line.Draw("Lsame")

        if drawLineLowerPanel:
            legEntry,yline = drawLineLowerPanel.split('::')
            line2 = ROOT.TF1("horiz_line_2",str(1+float(yline)),ratio.GetXaxis().GetBinLowEdge(1),ratio.GetXaxis().GetBinLowEdge(ratio.GetNbinsX()+1))
            line3 = ROOT.TF1("horiz_line_3",str(1-float(yline)),ratio.GetXaxis().GetBinLowEdge(1),ratio.GetXaxis().GetBinLowEdge(ratio.GetNbinsX()+1))
            line2.SetLineColor(ROOT.kBlue)
            line2.SetLineWidth(1)
            line2.Draw("Lsame")
            line3.SetLineColor(ROOT.kBlue)
            line3.SetLineWidth(1)
            line3.Draw("Lsame")
            x1leg2 = 0.2 if leftMargin > 0.1 else 0.07
            x2leg2 = 0.5 if leftMargin > 0.1 else 0.27
            y1leg2 = 0.25 if leftMargin > 0.1 else 0.3
            y2leg2 = 0.35 if leftMargin > 0.1 else 0.35
            leg2 = ROOT.TLegend(x1leg2, y1leg2, x2leg2, y2leg2)
            leg2.SetFillColor(0)
            leg2.SetFillStyle(0)
            leg2.SetBorderSize(0)
            leg2.AddEntry(line2,legEntry,"L")
            leg2.Draw("same")

        
        pad2.RedrawAxis("sameaxis")


    if draw_both0_noLog1_onlyLog2 != 2:
        canvas.SaveAs(outdir + canvasName + ".png")
        canvas.SaveAs(outdir + canvasName + ".pdf")

    if draw_both0_noLog1_onlyLog2 != 1:        
        if yAxisName == "a.u.": 
            h1.GetYaxis().SetRangeUser(max(0.0001,h1.GetMinimum()*0.8),h1.GetMaximum()*100)
        else:
            h1.GetYaxis().SetRangeUser(max(0.001,h1.GetMinimum()*0.8),h1.GetMaximum()*100)
        canvas.SetLogy()
        canvas.SaveAs(outdir + canvasName + "_logY.png")
        canvas.SaveAs(outdir + canvasName + "_logY.pdf")
        canvas.SetLogy(0)


################################################################

def drawDataAndMC(h1, h2,
                  labelXtmp="xaxis", labelYtmp="yaxis",
                  canvasName="default", outdir="./",
                  #rebinFactorX=0,
                  draw_both0_noLog1_onlyLog2=0,                  
                  leftMargin=0.15,
                  rightMargin=0.04,
                  labelRatioTmp="Data/pred.::0.5,1.5",
                  drawStatBox=False,
                  #legendCoords="0.15,0.35,0.8,0.9",  # x1,x2,y1,y2
                  legendCoords="0.15,0.65,0.85,0.9;3",  # for unrolled use 1 line # x1,x2,y1,y2
                  canvasSize="600,700",  # use X,Y to pass X and Y size     
                  lowerPanelHeight = 0.3,  # number from 0 to 1, 0.3 means 30% of space taken by lower panel. 0 means do not draw lower panel with relative error
                  #drawLineLowerPanel="lumi. uncertainty::0.025" # if not empty, draw band at 1+ number after ::, and add legend with title
                  #drawLineLowerPanel="", # if not empty, draw band at 1+ number after ::, and add legend with title
                  passCanvas=None,
                  lumi=None,
                  drawVertLines="", # "12,36": format --> N of sections (e.g: 12 pt bins), and N of bins in each section (e.g. 36 eta bins), assuming uniform bin width
                  textForLines=[],                       
                  moreText="",
                  moreTextLatex="",
                  invertRatio = False,  # make expected over observed if True
                  histMCpartialUnc = None,
                  histMCpartialUncLegEntry = "",
                  useDifferenceInLowerPanel = False,
                  noLegendLowerPanel = False,
                  legendEntries = []
                  ):

    # moreText is used to pass some text to write somewhere (TPaveText is used)
    # e.g.  "stuff::x1,y1,x2,y2"  where xi and yi are the coordinates for the text
    # one can add more lines using the ";" key. FOr example, "stuff1;stuff2::x1,y1,x2,y2"
    # the coordinates should be defined taking into account how many lines will be drawn
    # if the coordinates are not passed (no "::"), then default ones are used, but this might not be satisfactory

    # moreTextLatex is similar, but used TLatex, and the four coordinates are x1,y1,ypass,textsize
    # where x1 and y1 are the coordinates the first line, and ypass is how much below y1 the second line is (and so on for following lines)

    # h1 is data, h2 in MC

    #if (rebinFactorX): 
    #    if isinstance(rebinFactorX, int): h1.Rebin(rebinFactorX)
    #    # case in which rebinFactorX is a list of bin edges
    #    else:                             h1.Rebin(len(rebinFactorX)-1,"",array('d',rebinFactorX)) 

    xAxisName,setXAxisRangeFromUser,xmin,xmax = getAxisRangeFromUser(labelXtmp)
    yAxisName,setYAxisRangeFromUser,ymin,ymax = getAxisRangeFromUser(labelYtmp)
    yRatioAxisName,setRatioYAxisRangeFromUser,yminRatio,ymaxRatio = getAxisRangeFromUser(labelRatioTmp)

    yAxisTitleOffset = 1.45 if leftMargin > 0.1 else 0.6

    addStringToEnd(outdir,"/",notAddIfEndswithMatch=True)
    createPlotDirAndCopyPhp(outdir)

    cw,ch = canvasSize.split(',')
    #canvas = ROOT.TCanvas("canvas",h2D.GetTitle() if plotLabel == "ForceTitle" else "",700,625)
    canvas = passCanvas if passCanvas != None else ROOT.TCanvas("canvas","",int(cw),int(ch))
    canvas.SetTickx(1)
    canvas.SetTicky(1)
    canvas.cd()
    canvas.SetLeftMargin(leftMargin)
    canvas.SetRightMargin(rightMargin)
    canvas.cd()

    pad2 = 0
    if lowerPanelHeight: 
        canvas.SetBottomMargin(lowerPanelHeight)
        pad2 = ROOT.TPad("pad2","pad2",0,0.,1,0.9)
        pad2.SetTopMargin(1-lowerPanelHeight)
        pad2.SetRightMargin(rightMargin)
        pad2.SetLeftMargin(leftMargin)
        pad2.SetFillColor(0)
        pad2.SetGridy(1)
        pad2.SetFillStyle(0)


    frame = h1.Clone("frame")
    frame.GetXaxis().SetLabelSize(0.04)
    frame.SetStats(0)

    h1.SetLineColor(ROOT.kBlack)
    h1.SetMarkerColor(ROOT.kBlack)
    h1.SetMarkerStyle(20)
    h1.SetMarkerSize(1)

    #ymax = max(ymax, max(h1.GetBinContent(i)+h1.GetBinError(i) for i in range(1,h1.GetNbinsX()+1)))
    # if min and max were not set, set them based on histogram content
    if ymin == ymax == 0.0:
        ymin,ymax = getMinMaxHisto(h1,excludeEmpty=True,sumError=True)            
        ymin *= 0.9
        ymax *= (1.1 if leftMargin > 0.1 else 2.0)
        if ymin < 0: ymin = 0
        #print "drawSingleTH1() >>> Histo: %s     minY,maxY = %.2f, %.2f" % (h1.GetName(),ymin,ymax)

    # print "#### WARNING ####"
    # print "Hardcoding ymin = 0 in function drawDataAndMC(): change it if it is not what you need"
    # print "#################"
    # ymin = 0 # hardcoded

    if lowerPanelHeight:
        h1.GetXaxis().SetLabelSize(0)
        h1.GetXaxis().SetTitle("")  
    else:
        h1.GetXaxis().SetTitle(xAxisName)
        h1.GetXaxis().SetTitleOffset(1.2)
        h1.GetXaxis().SetTitleSize(0.05)
        h1.GetXaxis().SetLabelSize(0.04)
    h1.GetYaxis().SetTitle(yAxisName)
    h1.GetYaxis().SetTitleOffset(yAxisTitleOffset)
    h1.GetYaxis().SetTitleSize(0.05)
    h1.GetYaxis().SetLabelSize(0.04)
    h1.GetYaxis().SetRangeUser(ymin, ymax)    
    h1.GetYaxis().SetTickSize(0.01)
    if setXAxisRangeFromUser: h1.GetXaxis().SetRangeUser(xmin,xmax)
    h1.Draw("EP")
    #h1err = h1.Clone("h1err")
    #h1err.SetFillColor(ROOT.kRed+2)
    h2.SetFillStyle(1001)  # 3001 is better than 3002 for pdf, while 3002 is perfect for png
    h2.SetLineColor(ROOT.kRed+2)  # kGreen+2
    h2.SetFillColor(ROOT.kRed+1)    # kGreen
    h2.SetLineWidth(1)
    h2.Draw("E2 SAME")
    h2line = h2.Clone("h2line")
    h2line.SetFillColor(0)
    h3 = None
    if histMCpartialUnc != None:
        h3 = histMCpartialUnc.Clone("histMCpartialUnc")
        h3.SetFillColor(ROOT.kGreen)
        h3.SetFillStyle(1001)  # 3001, 3144 , 3244, 3003
        #h3.SetFillStyle(3244)  # 3144 , 3244, 3003
        h3.Draw("E2 SAME")
        #for i in range(1,1+h3.GetNbinsX()):
        #    print "PDF band: bin %d  val +/- error = %.3f +/- %.3f" % (i, h3.GetBinContent(i),h3.GetBinError(i))
    h2line.Draw("HIST SAME")
    h1.Draw("EP SAME")

    nColumnsLeg = 1
    if ";" in legendCoords: 
        nColumnsLeg = int(legendCoords.split(";")[1])
        #if histMCpartialUnc != None: nColumnsLeg = 3
    legcoords = [float(x) for x in (legendCoords.split(";")[0]).split(',')]
    lx1,lx2,ly1,ly2 = legcoords[0],legcoords[1],legcoords[2],legcoords[3]
    if histMCpartialUnc != None:
        ly2 = ly2 + 0.5 * (ly2 - ly1) # add one more row
    leg = ROOT.TLegend(lx1,ly1,lx2,ly2)
    #leg.SetFillColor(0)
    #leg.SetFillStyle(0)
    #leg.SetBorderSize(0)
    leg.SetNColumns(nColumnsLeg)
    if len(legendEntries):
        leg.AddEntry(h1,str(legendEntries[0]),"LPE")
        leg.AddEntry(h2,str(legendEntries[1]),"LF")        
    else:
        leg.AddEntry(h1,"measured","LPE")
        leg.AddEntry(h2,"MadGraph5_aMC@NLO","LF")
    if histMCpartialUnc != None:
        leg.AddEntry(h3,histMCpartialUncLegEntry,"LF")
    #leg.AddEntry(h1err,"Uncertainty","LF")
    leg.Draw("same")
    canvas.RedrawAxis("sameaxis")

    if drawStatBox:
        ROOT.gPad.Update()
        ROOT.gStyle.SetOptStat(1110)
        ROOT.gStyle.SetOptFit(1102)
    else:
        h1.SetStats(0)

    vertline = ROOT.TLine(36,0,36,canvas.GetUymax())
    vertline.SetLineColor(ROOT.kBlack)
    vertline.SetLineStyle(3) # 2 for denser hatches
    bintext = ROOT.TLatex()
    #bintext.SetNDC()
    bintext.SetTextSize(0.025)  # 0.03
    bintext.SetTextFont(42)
    if len(textForLines):
        bintext.SetTextAngle(45 if "#eta" in textForLines[0] else 30)

    if len(drawVertLines):
        #print "drawVertLines = " + drawVertLines
        nptBins = int(drawVertLines.split(',')[0])
        etarange = float(drawVertLines.split(',')[1])        
        offsetXaxisHist = h1.GetXaxis().GetBinLowEdge(0)
        sliceLabelOffset = 6. 
        for i in range(1,nptBins): # do not need line at canvas borders
            #vertline.DrawLine(offsetXaxisHist+etarange*i,0,offsetXaxisHist+etarange*i,canvas.GetUymax())
            vertline.DrawLine(etarange*i-offsetXaxisHist,0,etarange*i-offsetXaxisHist,ymax)
        if len(textForLines):
            for i in range(0,len(textForLines)): # we need nptBins texts
                #texoffset = 0.1 * (4 - (i%4))
                #ytext = (1. + texoffset)*ymax/2.  
                ytext = 0.6*ymax                  
                bintext.DrawLatex(etarange*i + etarange/sliceLabelOffset, ytext, textForLines[i])

    # redraw legend, or vertical lines appear on top of it
    leg.Draw("same")

    if len(moreText):
        realtext = moreText.split("::")[0]
        x1,y1,x2,y2 = 0.7,0.8,0.9,0.9
        if "::" in moreText:
            x1,y1,x2,y2 = (float(x) for x in (moreText.split("::")[1]).split(","))
        pavetext = ROOT.TPaveText(x1,y1,x2,y2,"NB NDC")
        for tx in realtext.split(";"):
            pavetext.AddText(tx)
        pavetext.SetFillColor(0)
        pavetext.SetFillStyle(0)
        pavetext.SetBorderSize(0)
        pavetext.SetLineColor(0)
        pavetext.Draw("same")

    if len(moreTextLatex):
        realtext = moreTextLatex.split("::")[0]
        x1,y1,ypass,textsize = 0.75,0.8,0.08,0.035
        if "::" in moreTextLatex:
            x1,y1,ypass,textsize = (float(x) for x in (moreTextLatex.split("::")[1]).split(","))            
        lat = ROOT.TLatex()
        lat.SetNDC();
        lat.SetTextFont(42)        
        lat.SetTextSize(textsize)
        for itx,tx in enumerate(realtext.split(";")):
            lat.DrawLatex(x1,y1-itx*ypass,tx)

  # TPaveText *pvtxt = NULL;
  # if (yAxisName == "a.u.") {
  #   pvtxt = new TPaveText(0.6,0.6,0.95,0.7, "BR NDC")
  #   pvtxt.SetFillColor(0)
  #   pvtxt.SetFillStyle(0)
  #   pvtxt.SetBorderSize(0)
  #   pvtxt.AddText(Form("norm num/den = %.2f +/- %.2f",IntegralRatio,ratioError))
  #   pvtxt.Draw()
  # }

    setTDRStyle()
    if leftMargin > 0.1:
        if lumi != None: CMS_lumi(canvas,lumi,True,False)
        else:            CMS_lumi(canvas,"",True,False)
    else:
        latCMS = ROOT.TLatex()
        latCMS.SetNDC();
        latCMS.SetTextFont(42)
        latCMS.SetTextSize(0.045)
        latCMS.DrawLatex(0.1, 0.95, '#bf{CMS} #it{Preliminary}')
        if lumi != None: latCMS.DrawLatex(0.85, 0.95, '%s fb^{-1} (13 TeV)' % lumi)
        else:            latCMS.DrawLatex(0.90, 0.95, '(13 TeV)' % lumi)

    if lowerPanelHeight:
        pad2.Draw()
        pad2.cd()

        frame.Reset("ICES")
        if setRatioYAxisRangeFromUser: frame.GetYaxis().SetRangeUser(yminRatio,ymaxRatio)
        #else:                          
        #frame.GetYaxis().SetRangeUser(0.5,1.5)
        frame.GetYaxis().SetNdivisions(5)
        frame.GetYaxis().SetTitle(yRatioAxisName)
        frame.GetYaxis().SetTitleOffset(yAxisTitleOffset)
        frame.GetYaxis().SetTitleSize(0.05)
        frame.GetYaxis().SetLabelSize(0.04)
        frame.GetYaxis().CenterTitle()
        frame.GetXaxis().SetTitle(xAxisName)
        if setXAxisRangeFromUser: frame.GetXaxis().SetRangeUser(xmin,xmax)
        frame.GetXaxis().SetTitleOffset(1.2)
        frame.GetXaxis().SetTitleSize(0.05)

        #ratio = copy.deepcopy(h1.Clone("ratio"))
        #den_noerr = copy.deepcopy(h2.Clone("den_noerr"))
        ratio = None
        den_noerr = None
        den = None
        if invertRatio:
            ratio = h2.Clone("ratio")
            den_noerr = h1.Clone("den_noerr")
            den = h1.Clone("den")
        else:
            ratio = h1.Clone("ratio")
            den_noerr = h2.Clone("den_noerr")
            den = h2.Clone("den")

        for iBin in range (1,den_noerr.GetNbinsX()+1):
            den_noerr.SetBinError(iBin,0.)

        if useDifferenceInLowerPanel:
            ratio.Add(den_noerr,-1.0)
            den.Add(den_noerr,-1.0)
        else:
            ratio.Divide(den_noerr)
            den.Divide(den_noerr)

        if invertRatio:
            if histMCpartialUnc == None:
                ratio.SetFillColor(ROOT.kGray+1)
                ratio.SetFillStyle(1001)  # make it solid again
                ratio.SetLineColor(ROOT.kGray+3)
            else:
                ratio.SetFillColor(ROOT.kRed+1) # kGreen+1
                ratio.SetFillStyle(1001)  # 1001 to make it solid again
                ratio.SetLineColor(ROOT.kRed+2) # kGreen+2                       
                ratio.SetLineWidth(1) # make it smaller when it is drawn on top of something
        else:
            if histMCpartialUnc == None:
                den.SetFillColor(ROOT.kGray+1)
                den.SetFillStyle(1001)  # make it solid again
                den.SetLineColor(ROOT.kRed)        
        if histMCpartialUnc != None:
            h3ratio = h3.Clone("h3ratio")
            if useDifferenceInLowerPanel:
                h3ratio.Add(den_noerr,-1.0)
            else:
                h3ratio.Divide(den_noerr)
            #h3ratio.SetFillStyle(3144) # 1001 for solid, 3144 instead of 3244, to have more dense hatches
            h3ratio.SetFillStyle(1001) # 
            h3ratio.SetFillColor(ROOT.kGreen+1)  # kRed-4
            h3ratio.SetLineColor(ROOT.kGreen+2)  # kRed-4

        frame.Draw()        
        if invertRatio:
            den.SetMarkerSize(0.85)
            den.SetMarkerStyle(20) 
            #den.Draw("EPsame")    # draw after red line, not now
            ratio.Draw("E2same")
        else:    
            ratio.SetMarkerSize(0.85)
            ratio.SetMarkerStyle(20) 
            den.Draw("E2same")
            #ratio.Draw("EPsame") # draw after red line, not now

        if histMCpartialUnc != None:
            h3ratio.Draw("E2 same")
        
        line = ROOT.TF1("horiz_line","0" if useDifferenceInLowerPanel else "1",
                        ratio.GetXaxis().GetBinLowEdge(1),ratio.GetXaxis().GetBinLowEdge(ratio.GetNbinsX()+1))
        line.SetLineWidth(1)
        line.Draw("Lsame")
        if invertRatio:
            ratioline = ratio.Clone("ratioline")
            ratioline.SetFillColor(0)
            if histMCpartialUnc != None and len(drawVertLines):
                # unrolled, here the red line would be too large and given the fluctuations is could be seen as an additional band
                ratioline.SetLineColor(histMCpartialUnc.GetLineColor())
            ratioline.SetFillStyle(0)
            if histMCpartialUnc != None: ratioline.Draw("HIST same") # to draw the line inside the band for the expected
            den.Draw("EPsame")
        else: 
            ratio.Draw("EPsame")

        x1leg2 = 0.15 if leftMargin > 0.1 else 0.07
        x2leg2 = 0.55 if leftMargin > 0.1 else 0.27
        y1leg2 = 0.28 if leftMargin > 0.1 else 0.3
        y2leg2 = 0.33 if leftMargin > 0.1 else 0.35
        leg2 = ROOT.TLegend(x1leg2, y1leg2, x2leg2, y2leg2)
        #leg2 = ROOT.TLegend(0.07,0.30,0.27,0.35)
        leg2.SetFillColor(0)
        leg2.SetFillStyle(0)
        leg2.SetBorderSize(0)
        if invertRatio:
            leg2.AddEntry(ratio,"total theory uncertainty","LF")
        else:
            leg2.AddEntry(den,"total theory uncertainty","LF")
        if not noLegendLowerPanel: leg2.Draw("same")
        if histMCpartialUnc != None:
            leg2bis = ROOT.TLegend(x1leg2 + 0.45, y1leg2, x2leg2 + 0.45, y2leg2)
            leg2bis.SetFillColor(0)
            leg2bis.SetFillStyle(0)
            leg2bis.SetBorderSize(0)
            leg2bis.AddEntry(h3ratio,histMCpartialUncLegEntry,"LF")
            if not noLegendLowerPanel: leg2bis.Draw("same")

        pad2.RedrawAxis("sameaxis")


    if draw_both0_noLog1_onlyLog2 != 2:
        canvas.SaveAs(outdir + canvasName + ".png")
        canvas.SaveAs(outdir + canvasName + ".pdf")

    if draw_both0_noLog1_onlyLog2 != 1:        
        if yAxisName == "a.u.": 
            h1.GetYaxis().SetRangeUser(max(0.0001,h1.GetMinimum()*0.8),h1.GetMaximum()*100)
        else:
            h1.GetYaxis().SetRangeUser(max(0.001,h1.GetMinimum()*0.8),h1.GetMaximum()*100)
        canvas.SetLogy()
        canvas.SaveAs(outdir + canvasName + "_logY.png")
        canvas.SaveAs(outdir + canvasName + "_logY.pdf")
        canvas.SetLogy(0)
        
          
#########################################################################

def drawTH1dataMCstack(h1, thestack, 
                       labelXtmp="xaxis", labelYtmp="yaxis",
                       canvasName="default", 
                       outdir="./",
                       legend=None,
                       ratioPadYaxisNameTmp="data/MC::0.5,1.5", 
                       draw_both0_noLog1_onlyLog2=0,
                       #minFractionToBeInLegend=0.001,
                       fillStyle=3001,
                       leftMargin=0.16,
                       rightMargin=0.05,
                       nContours=50,
                       palette=55,
                       canvasSize="700,625",
                       passCanvas=None,
                       normalizeMCToData=False,
                       hErrStack=None,   # might need to define an error on the stack in a special way
                       lumi=None,
                       yRangeScaleFactor=1.5, # if range of y axis is not explicitely passed, use (max-min) times this value
                       wideCanvas=False,
                       drawVertLines="", # "12,36": format --> N of sections (e.g: 12 pt bins), and N of bins in each section (e.g. 36 eta bins), assuming uniform bin width
                       textForLines=[], 
                       etaptbinning=[]
                       ):

    # if normalizing stack to same area as data, we need to modify the stack
    # however, the stack might be used outside the function. In order to avoid any changes in the stack, it is copied here just for the plot

    ROOT.TH1.SetDefaultSumw2()
    adjustSettings_CMS_lumi()
    
    ROOT.gStyle.SetPalette(palette)  # 55:raibow palette ; 57: kBird (blue to yellow, default) ; 107 kVisibleSpectrum ; 77 kDarkRainBow 
    ROOT.gStyle.SetNumberContours(nContours) # default is 20 

    labelX,setXAxisRangeFromUser,xmin,xmax = getAxisRangeFromUser(labelXtmp)
    labelY,setYAxisRangeFromUser,ymin,ymax = getAxisRangeFromUser(labelYtmp)
    labelRatioY,setRatioYAxisRangeFromUser,yminRatio,ymaxRatio = getAxisRangeFromUser(ratioPadYaxisNameTmp)
    
    cw,ch = canvasSize.split(',')
    #canvas = ROOT.TCanvas("canvas",h2D.GetTitle() if plotLabel == "ForceTitle" else "",700,625)    
    canvas = passCanvas if passCanvas != None else ROOT.TCanvas("canvas","",int(cw),int(ch))
    canvas.SetTickx(1)
    canvas.SetTicky(1)
    canvas.SetLeftMargin(leftMargin)
    canvas.SetRightMargin(rightMargin)
    canvas.SetBottomMargin(0.3)
    canvas.cd()

    addStringToEnd(outdir,"/",notAddIfEndswithMatch=True)
    createPlotDirAndCopyPhp(outdir)

    dataNorm = h1.Integral()
    stackNorm = 0.0

    #dummystack = thestack
    dummystack = ROOT.THStack("dummy_{sn}".format(sn=thestack.GetName()),"")
    for hist in thestack.GetHists():        
        stackNorm += hist.Integral()
    for hist in thestack.GetHists():        
        hnew = copy.deepcopy(hist.Clone("dummy_{hn}".format(hn=hist.GetName())))
        if normalizeMCToData:
            hnew.Scale(dataNorm/stackNorm)
        dummystack.Add(hnew)    
        
    stackCopy = dummystack.GetStack().Last() # used to make ratioplot without affecting the plot and setting maximum
    # the error of the last should be the sum in quadrature of the errors of single components, as the Last is the sum of them
    # however, better to recreate it
    stackErr = stackCopy
    if hErrStack != None:
        stackErr = copy.deepcopy(hErrStack.Clone("stackErr"))

    print "drawTH1dataMCstack():  integral(data):  " + str(h1.Integral()) 
    print "drawTH1dataMCstack():  integral(stack): " + str(stackCopy.Integral()) 
    print "drawTH1dataMCstack():  integral(herr):  " + str(stackErr.Integral()) 

    h1.SetStats(0)
    titleBackup = h1.GetTitle()
    h1.SetTitle("")

    pad2 = ROOT.TPad("pad2","pad2",0,0.,1,0.9)
    pad2.SetTopMargin(0.7)
    pad2.SetRightMargin(rightMargin)
    pad2.SetLeftMargin(leftMargin)
    pad2.SetFillColor(0)
    pad2.SetGridy(1)
    pad2.SetFillStyle(0)
    
    frame = h1.Clone("frame")
    frame.GetXaxis().SetLabelSize(0.04)
    frame.SetStats(0)

    h1.SetLineColor(ROOT.kBlack)
    h1.SetMarkerColor(ROOT.kBlack)
    h1.SetMarkerStyle(20)
    h1.SetMarkerSize(1)

    h1.GetXaxis().SetLabelSize(0)
    h1.GetXaxis().SetTitle("")
    h1.GetYaxis().SetTitle(labelY)
    h1.GetYaxis().SetTitleOffset(0.5 if wideCanvas else 1.5)
    h1.GetYaxis().SetTitleSize(0.05)
    h1.GetYaxis().SetLabelSize(0.04)
    h1.GetYaxis().SetTickSize(0.01)
    ymaxBackup = 0
    if setYAxisRangeFromUser: 
        ymaxBackup = ymax
        h1.GetYaxis().SetRangeUser(ymin,ymax)
    else:
        ymaxBackup = max(h1.GetBinContent(h1.GetMaximumBin()),stackCopy.GetBinContent(stackCopy.GetMaximumBin())) * yRangeScaleFactor
        h1.GetYaxis().SetRangeUser(0.0, ymaxBackup)
    if setXAxisRangeFromUser: h1.GetXaxis().SetRangeUser(xmin,xmax)
    h1.Draw("EP")
    dummystack.Draw("HIST SAME")
    h1.Draw("EP SAME")

    vertline = ROOT.TLine(36,0,36,canvas.GetUymax())
    vertline.SetLineColor(ROOT.kBlack)
    vertline.SetLineStyle(3) # 2 larger hatches
    bintext = ROOT.TLatex()
    #bintext.SetNDC()
    bintext.SetTextSize(0.03)
    bintext.SetTextFont(42)
    if len(textForLines) > 15:
        bintext.SetTextAngle(10)

    if len(drawVertLines):
        #print "drawVertLines = " + drawVertLines
        nptBins = int(drawVertLines.split(',')[0])
        etarange = float(drawVertLines.split(',')[1])        
        for i in range(1,nptBins): # do not need line at canvas borders
            #vertline.DrawLine(etarange*i,0,etarange*i,canvas.GetUymax())
            vertline.DrawLine(etarange*i,0,etarange*i,ymaxBackup)
        if len(textForLines):
            offsetText = etarange / (6. if len(textForLines) > 15 else 4.)
            for i in range(0,len(textForLines)): # we need nptBins texts
                bintext.DrawLatex(etarange*i + offsetText, 1.1*ymaxBackup/2., textForLines[i])

    # legend.SetFillColor(0)
    # legend.SetFillStyle(0)
    # legend.SetBorderSize(0)
    legend.Draw("same")
    canvas.RedrawAxis("sameaxis")

    reduceSize = False
    offset = 0
    # check whether the Y axis will have exponential notatio
    if h1.GetBinContent(h1.GetMaximumBin()) > 1000000:
        reduceSize = True
        offset = 0.1
    if wideCanvas: 
        offset = 0.1
        latCMS = ROOT.TLatex()
        latCMS.SetNDC();
        latCMS.SetTextFont(42)
        latCMS.SetTextSize(0.045)
        latCMS.DrawLatex(0.1, 0.95, '#bf{CMS} #it{Preliminary}')
        if lumi != None: latCMS.DrawLatex(0.85, 0.95, '%s fb^{-1} (13 TeV)' % lumi)
        else:            latCMS.DrawLatex(0.90, 0.95, '(13 TeV)' % lumi)
    else:    
        if lumi != None: CMS_lumi(canvas,lumi,True,False, reduceSize, offset)
        else:            CMS_lumi(canvas,"",True,False)    

    setTDRStyle()

    pad2.Draw();
    pad2.cd();

    frame.Reset("ICES")
    if setRatioYAxisRangeFromUser: frame.GetYaxis().SetRangeUser(yminRatio,ymaxRatio)
    #else:                          
    #frame.GetYaxis().SetRangeUser(0.5,1.5)
    frame.GetYaxis().SetNdivisions(5)
    frame.GetYaxis().SetTitle(labelRatioY)
    frame.GetYaxis().SetTitleOffset(0.5 if wideCanvas else 1.5)
    frame.GetYaxis().SetTitleSize(0.05)
    frame.GetYaxis().SetLabelSize(0.04)
    frame.GetYaxis().CenterTitle()
    frame.GetXaxis().SetTitle(labelX)
    if setXAxisRangeFromUser: frame.GetXaxis().SetRangeUser(xmin,xmax)
    frame.GetXaxis().SetTitleOffset(1.2)
    frame.GetXaxis().SetTitleSize(0.05)

    #ratio = copy.deepcopy(h1.Clone("ratio"))
    #den_noerr = copy.deepcopy(stackErr.Clone("den_noerr"))
    ratio = h1.Clone("ratio")
    den_noerr = stackErr.Clone("den_noerr")
    den = stackErr.Clone("den")
    for iBin in range (1,den_noerr.GetNbinsX()+1):
        den_noerr.SetBinError(iBin,0.)

    ratio.Divide(den_noerr)
    den.Divide(den_noerr)
    den.SetFillColor(ROOT.kCyan)
    den.SetFillStyle(1001)  # make it solid again
    #den.SetLineColor(ROOT.kRed)
    frame.Draw()        
    ratio.SetMarkerSize(0.85)
    ratio.SetMarkerStyle(20) 
    den.Draw("E2same")
    ratio.Draw("EPsame")

    # if not "unrolled_" in canvasName:
    #     for i in range(1,1+ratio.GetNbinsX()):
    #         print "Error data bin {bin}: {val}".format(bin=i,val=ratio.GetBinError(i))

    line = ROOT.TF1("horiz_line","1",ratio.GetXaxis().GetBinLowEdge(1),ratio.GetXaxis().GetBinLowEdge(ratio.GetNbinsX()+1))
    line.SetLineColor(ROOT.kRed)
    line.SetLineWidth(1)  # 1, not 2, which is too wide for canvas with large width
    line.Draw("Lsame")

    leg2 = ROOT.TLegend(0.2,0.25,0.4,0.30)
    leg2.SetFillColor(0)
    leg2.SetFillStyle(0)
    leg2.SetBorderSize(0)
    leg2.AddEntry(den,"tot. unc. exp.","LF")
    leg2.Draw("same")

    pad2.RedrawAxis("sameaxis")


    if draw_both0_noLog1_onlyLog2 != 2:
        canvas.SaveAs(outdir + canvasName + ".png")
        canvas.SaveAs(outdir + canvasName + ".pdf")

    if draw_both0_noLog1_onlyLog2 != 1:        
        if yAxisName == "a.u.": 
            h1.GetYaxis().SetRangeUser(max(0.0001,h1.GetMinimum()*0.8),h1.GetMaximum()*100)
        else:
            h1.GetYaxis().SetRangeUser(max(0.001,h1.GetMinimum()*0.8),h1.GetMaximum()*100)
            canvas.SetLogy()
            canvas.SaveAs(outdir + canvasName + "_logY.png")
            canvas.SaveAs(outdir + canvasName + "_logY.pdf")
            canvas.SetLogy(0)
            

    h1.SetTitle(titleBackup)
  
    if "unrolled" in canvasName:

        _canvas_pull.SetTickx(1)
        _canvas_pull.SetTicky(1)
        _canvas_pull.SetGridx(1)
        _canvas_pull.SetGridy(1)
        _canvas_pull.SetTopMargin(0.1)
        _canvas_pull.SetBottomMargin(0.12)
        _canvas_pull.SetLeftMargin(0.12)
        _canvas_pull.SetRightMargin(0.04)
        # make pulls
        pulltitle = "unrolled W^{{{ch}}} {pf}".format(ch="+" if "plus" in canvasName else "-", pf="postfit" if "postfit" in canvasName else "prefit")
        hpull = ROOT.TH1D("hpull_"+canvasName,pulltitle,51,-5,5)    
        hpull.SetStats(1)
        _canvas_pull.cd()
        for i in range (1,ratio.GetNbinsX()+1):
            errTotDen = ratio.GetBinError(i)*ratio.GetBinError(i) + den.GetBinError(i)*den.GetBinError(i)            
            if errTotDen > 0.0:
                hpull.Fill((ratio.GetBinContent(i)-1)/math.sqrt(errTotDen))
        hpull.Draw("HIST")
        hpull.GetXaxis().SetTitle("pull")
        hpull.GetYaxis().SetTitle("Events")
        hpull.SetLineColor(ROOT.kBlack)
        hpull.SetLineWidth(2)
        ROOT.gStyle.SetOptTitle(1)                
        ROOT.gStyle.SetOptStat(111110)
        ROOT.gStyle.SetOptFit(1102)
        _canvas_pull.RedrawAxis("sameaxis")
        _canvas_pull.SaveAs(outdir + "pull_" + canvasName + ".png")    
        _canvas_pull.SaveAs(outdir + "pull_" + canvasName + ".pdf")

        if len(etaptbinning):
            _canvas_pull.SetGridx(0)
            _canvas_pull.SetGridy(0)            
            _canvas_pull.SetRightMargin(0.16)            
            h2pull = ROOT.TH2D("h2pull_"+canvasName, pulltitle.replace("unrolled","rolled") ,
                               etaptbinning[0], array('d', etaptbinning[1]), etaptbinning[2], array('d', etaptbinning[3]))
            hpull.Reset("ICESM")  # will use again for pulls in EE only
            hpullEEp = hpull.Clone("hpullEEp")
            hpullEEm = hpull.Clone("hpullEEm")
            for i in range (1,ratio.GetNbinsX()+1):
                etabin = (i-1)%etaptbinning[0] + 1
                ptbin = (i-1)/etaptbinning[0] + 1
                errTotDen = ratio.GetBinError(i)*ratio.GetBinError(i) + den.GetBinError(i)*den.GetBinError(i)            
                if errTotDen > 0.0:
                    pullVal = (ratio.GetBinContent(i)-1)/math.sqrt(errTotDen)
                    h2pull.SetBinContent(etabin,ptbin, pullVal)
                    if abs(etaptbinning[1][etabin]) >= 1.499: 
                        hpull.Fill(pullVal)
                        if etaptbinning[1][etabin] > 0: hpullEEp.Fill(pullVal)
                        else:                           hpullEEm.Fill(pullVal)
            h2pull.GetXaxis().SetTitle("%s #eta" % "muon" if "muon" in labelX else "electron")
            h2pull.GetYaxis().SetTitle("%s p_{T}" % "muon" if "muon" in labelX else "electron")
            h2pull.GetZaxis().SetTitle("pull")
            h2pull.SetStats(0)
            h2pull.GetZaxis().SetRangeUser(-3,3)
            h2pull.Draw("COLZ")
            _canvas_pull.RedrawAxis("sameaxis")
            _canvas_pull.SaveAs(outdir + "pull2D_" + canvasName + ".png")
            _canvas_pull.SaveAs(outdir + "pull2D_" + canvasName + ".pdf")

            # add pulls for EE only
            _canvas_pull.SetTickx(1)
            _canvas_pull.SetTicky(1)
            _canvas_pull.SetGridx(1)
            _canvas_pull.SetGridy(1)
            _canvas_pull.SetTopMargin(0.1)
            _canvas_pull.SetBottomMargin(0.12)
            _canvas_pull.SetLeftMargin(0.12)
            _canvas_pull.SetRightMargin(0.04)
            hpull.Draw("HIST")
            hpull.GetXaxis().SetTitle("pull (only |#eta| >= 1.5)")
            hpull.GetYaxis().SetTitle("Events")
            hpullEEp.SetLineWidth(2)
            hpullEEp.SetLineColor(ROOT.kOrange+2)
            hpullEEp.SetFillColor(ROOT.kOrange+1)
            hpullEEp.SetFillStyle(3001)
            hpullEEm.SetLineWidth(2)
            hpullEEm.SetLineColor(ROOT.kBlue+2)
            hpullEEm.SetFillColor(ROOT.kAzure+1)
            hpullEEm.SetFillStyle(3244)    
            hpullEEp.Draw("HIST SAME")
            hpullEEm.Draw("HIST SAME")
            hpull.Draw("HIST SAME")
            legEE = ROOT.TLegend(0.15,0.5,0.45,0.8)
            legEE.SetFillStyle(0)
            legEE.SetFillColor(0)
            legEE.SetBorderSize(0)
            legEE.AddEntry(hpull,"|#eta| >= 1.5","L")
            legEE.AddEntry(hpullEEp,"#eta >= 1.5","LF")
            legEE.AddEntry(hpullEEm,"#eta<= -1.5","LF")
            legEE.Draw("same")
            # pEEtext = ROOT.TLatex()
            # pEEtext.SetTextSize(0.1)
            # pEEtext.SetTextFont(42)
            # pEEtext.SetTextColor(ROOT.kRed+2)
            # pEEtext.DrawLatex(0.1,0.6,"#eta >=  1.5")
            # pEEtext2 = ROOT.TLatex()
            # pEEtext2.SetTextSize(0.1)
            # pEEtext2.SetTextFont(42)
            # pEEtext2.SetTextColor(ROOT.kBlue+2)
            # pEEtext2.DrawLatex(0.1,0.3,"#eta <= -1.5")
            _canvas_pull.RedrawAxis("sameaxis")
            _canvas_pull.SaveAs(outdir + "pull_onlyEE_" + canvasName + ".png")    
            _canvas_pull.SaveAs(outdir + "pull_onlyEE_" + canvasName + ".pdf")
            
#########################################################################

def drawMuElComparison(hlep, hmu, hel,
                       labelXtmp="xaxis", labelYtmp="yaxis",
                       canvasName="default", outdir="./",
                       #rebinFactorX=0,
                       draw_both0_noLog1_onlyLog2=0,                  
                       leftMargin=0.15,
                       rightMargin=0.04,
                       labelRatioTmp="Data/pred.::0.5,1.5",
                       drawStatBox=False,
                       legendCoords="0.15,0.35,0.8,0.9",  # x1,x2,y1,y2
                       canvasSize="600,700",  # use X,Y to pass X and Y size     
                       lowerPanelHeight = 0.3,  # number from 0 to 1, 0.3 means 30% of space taken by lower panel. 0 means do not draw lower panel with relative error
                       #drawLineLowerPanel="lumi. uncertainty::0.025" # if not empty, draw band at 1+ number after ::, and add legend with title
                       #drawLineLowerPanel="", # if not empty, draw band at 1+ number after ::, and add legend with title
                       passCanvas=None,
                       lumi=None,
                       drawVertLines="", # "12,36": format --> N of sections (e.g: 12 pt bins), and N of bins in each section (e.g. 36 eta bins), assuming uniform bin width
                       textForLines=[],                       
                       moreText="",
                       moreTextLatex="",
                       skipPreliminary=True,
                       saveRootFileForHepdata = True                       
                   ):

    # moreText is used to pass some text to write somewhere (TPaveText is used)
    # e.g.  "stuff::x1,y1,x2,y2"  where xi and yi are the coordinates for the text
    # one can add more lines using the ";" key. FOr example, "stuff1;stuff2::x1,y1,x2,y2"
    # the coordinates should be defined taking into account how many lines will be drawn
    # if the coordinates are not passed (no "::"), then default ones are used, but this might not be satisfactory

    # moreTextLatex is similar, but used TLatex, and the four coordinates are x1,y1,ypass,textsize
    # where x1 and y1 are the coordinates the first line, and ypass is how much below y1 the second line is (and so on for following lines)

    # h1 is data, h2 in MC
    
    combColor = ROOT.kGray # ROOT.kYellow

    #if (rebinFactorX): 
    #    if isinstance(rebinFactorX, int): h1.Rebin(rebinFactorX)
    #    # case in which rebinFactorX is a list of bin edges
    #    else:                             h1.Rebin(len(rebinFactorX)-1,"",array('d',rebinFactorX)) 

    xAxisName,setXAxisRangeFromUser,xmin,xmax = getAxisRangeFromUser(labelXtmp)
    yAxisName,setYAxisRangeFromUser,ymin,ymax = getAxisRangeFromUser(labelYtmp)
    yRatioAxisName,setRatioYAxisRangeFromUser,yminRatio,ymaxRatio = getAxisRangeFromUser(labelRatioTmp)

    yAxisTitleOffset = 1.45 if leftMargin > 0.1 else 0.6

    addStringToEnd(outdir,"/",notAddIfEndswithMatch=True)
    createPlotDirAndCopyPhp(outdir)

    cw,ch = canvasSize.split(',')
    #canvas = ROOT.TCanvas("canvas",h2D.GetTitle() if plotLabel == "ForceTitle" else "",700,625)
    canvas = passCanvas if passCanvas != None else ROOT.TCanvas("canvas","",int(cw),int(ch))
    canvas.SetTickx(1)
    canvas.SetTicky(1)
    canvas.cd()
    canvas.SetLeftMargin(leftMargin)
    canvas.SetRightMargin(rightMargin)
    canvas.cd()

    pad2 = 0
    if lowerPanelHeight: 
        canvas.SetBottomMargin(lowerPanelHeight)
        pad2 = ROOT.TPad("pad2","pad2",0,0.,1,0.9)
        pad2.SetTopMargin(1-lowerPanelHeight)
        pad2.SetRightMargin(rightMargin)
        pad2.SetLeftMargin(leftMargin)
        pad2.SetFillColor(0)
        pad2.SetGridy(1)
        pad2.SetFillStyle(0)


    frame = hlep.Clone("frame")
    frame.GetXaxis().SetLabelSize(0.04 if leftMargin > 0.1 else 0.05)
    frame.SetStats(0)

    #hlep.SetLineColor(ROOT.kGreen+2) # kBlack
    #hlep.SetMarkerColor(ROOT.kGreen+2)
    hlep.SetLineColor(combColor+3) # kBlack
    hlep.SetMarkerColor(combColor+3)
    hlep.SetFillColor(combColor)
    hlep.SetMarkerStyle(20)
    hlep.SetMarkerSize(1.3)

    #ymax = max(ymax, max(hlep.GetBinContent(i)+hlep.GetBinError(i) for i in range(1,hlep.GetNbinsX()+1)))
    # if min and max were not set, set them based on histogram content
    if ymin == ymax == 0.0:
        ymin,ymax = getMinMaxHisto(hlep,excludeEmpty=True,sumError=True)            
        ymin *= 0.9
        ymax *= (1.1 if leftMargin > 0.1 else 2.0)
        if ymin < 0: ymin = 0
        #print "drawSingleTH1() >>> Histo: %s     minY,maxY = %.2f, %.2f" % (hlep.GetName(),ymin,ymax)

    # print "#### WARNING ####"
    # print "Hardcoding ymin = 0 in function drawDataAndMC(): change it if it is not what you need"
    # print "#################"
    # ymin = 0 # hardcoded

    if lowerPanelHeight:
        hlep.GetXaxis().SetLabelSize(0)
        hlep.GetXaxis().SetTitle("")  
    else:
        hlep.GetXaxis().SetTitle(xAxisName)
        hlep.GetXaxis().SetTitleOffset(1.2)
        hlep.GetXaxis().SetTitleSize(0.05)
        hlep.GetXaxis().SetLabelSize(0.04 if leftMargin > 0.1 else 0.05)
    hlep.GetYaxis().SetTitle(yAxisName)
    hlep.GetYaxis().SetTitleOffset(yAxisTitleOffset)
    hlep.GetYaxis().SetTitleSize(0.05)
    hlep.GetYaxis().SetLabelSize(0.04)
    hlep.GetYaxis().SetDecimals()
    hlep.GetYaxis().SetRangeUser(ymin, ymax)    
    hlep.GetYaxis().SetTickSize(0.01)
    if setXAxisRangeFromUser: hlep.GetXaxis().SetRangeUser(xmin,xmax)
    hlep.Draw("E2P0")
    #hleperr = hlep.Clone("hleperr")
    #hleperr.SetFillColor(ROOT.kRed+2)

    # hmu.SetFillStyle(3001)  # 3001 is better than 3002 for pdf, while 3002 is perfect for png
    # hmu.SetLineColor(ROOT.kGreen+2)  #kRed+2
    # hmu.SetFillColor(ROOT.kGreen)
    # hmu.SetLineWidth(1)
    # hmu.Draw("E2 SAME")
    # hmuline = hmu.Clone("hmuline")
    # hmuline.SetFillColor(0)
    # h3 = None
    # hmuline.Draw("HIST SAME")
    # hlep.Draw("EP SAME")

    hmu.SetLineColor(ROOT.kBlue) # kgreen+2
    hmu.SetMarkerColor(ROOT.kBlue+1)
    hmu.SetMarkerStyle(21)
    hmu.SetMarkerSize(1.3)
    hmu.Draw("E0X0P0 SAME")

    hel.SetLineColor(ROOT.kRed+1)  # kRed+2
    hel.SetMarkerColor(ROOT.kRed+2)
    hel.SetMarkerStyle(21)
    hel.SetMarkerSize(1.3)
    hel.Draw("E0X0P0 SAME")

    nColumnsLeg = 1
    if ";" in legendCoords: 
        nColumnsLeg = int(legendCoords.split(";")[1])
    legcoords = [float(x) for x in (legendCoords.split(";")[0]).split(',')]
    lx1,lx2,ly1,ly2 = legcoords[0],legcoords[1],legcoords[2],legcoords[3]
    if leftMargin < 0.1:
        if len(moreText):
            realtext = moreText.split("::")[0]
            lx2 = lx2 + 0.333 * (lx2 - lx1) # increase length of legend
            nColumnsLeg += 1

    leg = ROOT.TLegend(lx1,ly1,lx2,ly2)
    #leg.SetFillColor(0)
    #leg.SetFillStyle(0)
    #leg.SetBorderSize(0)
    leg.SetNColumns(nColumnsLeg)
    leg.AddEntry(hlep,"combination","LFPE")
    leg.AddEntry(hmu,"muon","LPE")
    leg.AddEntry(hel,"electron","LPE")
    if len(moreText):
        leg.AddEntry(None,realtext,"")
    #leg.AddEntry(hleperr,"Uncertainty","LF")
    leg.Draw("same")
    canvas.RedrawAxis("sameaxis")

    if drawStatBox:
        ROOT.gPad.Update()
        ROOT.gStyle.SetOptStat(1110)
        ROOT.gStyle.SetOptFit(1102)
    else:
        hlep.SetStats(0)

    vertline = ROOT.TLine(36,0,36,canvas.GetUymax())
    vertline.SetLineColor(ROOT.kBlack)
    vertline.SetLineStyle(3) # 2 for denser hatches
    bintext = ROOT.TLatex()
    #bintext.SetNDC()
    bintext.SetTextSize(0.038)  # 0.03
    bintext.SetTextFont(42)
    if len(textForLines):
        bintext.SetTextAngle(18)

    if len(drawVertLines):
        #print "drawVertLines = " + drawVertLines
        nptBins = int(drawVertLines.split(',')[0])
        etarange = float(drawVertLines.split(',')[1])        
        offsetXaxisHist = hlep.GetXaxis().GetBinLowEdge(0)
        sliceLabelOffset = 6. 
        for i in range(1,nptBins): # do not need line at canvas borders
            #vertline.DrawLine(offsetXaxisHist+etarange*i,0,offsetXaxisHist+etarange*i,canvas.GetUymax())
            vertline.DrawLine(etarange*i-offsetXaxisHist,0,etarange*i-offsetXaxisHist,ymax)
        if len(textForLines):
            for i in range(0,len(textForLines)): # we need nptBins texts
                #texoffset = 0.1 * (4 - (i%4))
                #ytext = (1. + texoffset)*ymax/2.  
                factor = 0.62 if "#eta" in xAxisName else 0.67
                ytext = factor*ymax                  
                bintext.DrawLatex(etarange*i + etarange/sliceLabelOffset, ytext, textForLines[i])

    # redraw legend, or vertical lines appear on top of it
    leg.Draw("same")

    if len(moreText) and leftMargin > 0.1:
        realtext = moreText.split("::")[0]
        x1,y1,x2,y2 = 0.7,0.8,0.9,0.9
        if "::" in moreText:
            x1,y1,x2,y2 = (float(x) for x in (moreText.split("::")[1]).split(","))
        pavetext = ROOT.TPaveText(x1,y1,x2,y2,"NB NDC")
        for tx in realtext.split(";"):
            pavetext.AddText(tx)
        pavetext.SetFillColor(0)
        pavetext.SetFillStyle(0)
        pavetext.SetBorderSize(0)
        pavetext.SetLineColor(0)
        pavetext.Draw("same")

    if len(moreTextLatex):
        realtext = moreTextLatex.split("::")[0]
        x1,y1,ypass,textsize = 0.75,0.8,0.08,0.035
        if "::" in moreTextLatex:
            x1,y1,ypass,textsize = (float(x) for x in (moreTextLatex.split("::")[1]).split(","))            
        lat = ROOT.TLatex()
        lat.SetNDC();
        lat.SetTextFont(42)        
        lat.SetTextSize(textsize)
        for itx,tx in enumerate(realtext.split(";")):
            lat.DrawLatex(x1,y1-itx*ypass,tx)

  # TPaveText *pvtxt = NULL;
  # if (yAxisName == "a.u.") {
  #   pvtxt = new TPaveText(0.6,0.6,0.95,0.7, "BR NDC")
  #   pvtxt.SetFillColor(0)
  #   pvtxt.SetFillStyle(0)
  #   pvtxt.SetBorderSize(0)
  #   pvtxt.AddText(Form("norm num/den = %.2f +/- %.2f",IntegralRatio,ratioError))
  #   pvtxt.Draw()
  # }

    setTDRStyle()
    if leftMargin > 0.1:
        if lumi != None: CMS_lumi(canvas,lumi,True,skipPreliminary, offsetLumi=0.02)
        else:            CMS_lumi(canvas,"",True,skipPreliminary)
    else:
        latCMS = ROOT.TLatex()
        latCMS.SetNDC();
        latCMS.SetTextFont(42)
        latCMS.SetTextSize(0.045)
        latCMS.DrawLatex(leftMargin, 0.95, '#bf{CMS}' + (' #it{Preliminary}' if not skipPreliminary else ''))
        if lumi != None: latCMS.DrawLatex(0.85+(0.04-rightMargin), 0.95, '%s fb^{-1} (13 TeV)' % lumi)
        else:            latCMS.DrawLatex(0.90+(0.04-rightMargin), 0.95, '(13 TeV)' % lumi)

    if lowerPanelHeight:
        pad2.Draw()
        pad2.cd()

        frame.Reset("ICES")
        if setRatioYAxisRangeFromUser: frame.GetYaxis().SetRangeUser(yminRatio,ymaxRatio)
        #else:                          
        #frame.GetYaxis().SetRangeUser(0.5,1.5)
        frame.GetYaxis().SetDecimals()
        frame.GetYaxis().SetNdivisions(5)
        frame.GetYaxis().SetTitle(yRatioAxisName)
        frame.GetYaxis().SetTitleOffset(yAxisTitleOffset)
        frame.GetYaxis().SetTitleSize(0.05)
        frame.GetYaxis().SetLabelSize(0.04)
        frame.GetYaxis().CenterTitle()
        frame.GetXaxis().SetTitle(xAxisName)
        if setXAxisRangeFromUser: frame.GetXaxis().SetRangeUser(xmin,xmax)
        frame.GetXaxis().SetTitleOffset(1.2)
        frame.GetXaxis().SetTitleSize(0.05)

        #ratio = copy.deepcopy(hlep.Clone("ratio"))
        #den_noerr = copy.deepcopy(h2.Clone("den_noerr"))
        ratiomu = None
        ratioel = None
        den_noerr = None
        den = None
        ratiomu = hmu.Clone("ratiomu")
        ratioel = hel.Clone("ratioel")
        den_noerr = hlep.Clone("den_noerr")
        den = hlep.Clone("den")
       
        for iBin in range (1,den_noerr.GetNbinsX()+1):
            den_noerr.SetBinError(iBin,0.)

        ratiomu.Divide(den_noerr)
        ratioel.Divide(den_noerr)
        den.Divide(den_noerr)
        #den.SetFillColor(ROOT.kGray)
        #den.SetLineColor(ROOT.kBlack)        
        den.SetFillColor(combColor)
        den.SetLineColor(ROOT.kBlack)        
        den.SetFillStyle(1001)  # make it solid again
        den.SetLineWidth(2)        

        frame.Draw()        
        ratiomu.SetMarkerSize(1.6)
        ratioel.SetMarkerSize(1.6)
        den.SetMarkerSize(0) 
        den.Draw("E2same")
        ratiomu.Draw("E0X0P0 same")
        ratioel.Draw("E0X0P0 same")
        
 
        line = ROOT.TF1("horiz_line","1",den.GetXaxis().GetBinLowEdge(1),den.GetXaxis().GetBinLowEdge(den.GetNbinsX()+1))
        line.SetLineColor(ROOT.kBlack)
        line.SetLineWidth(1)
        line.Draw("Lsame")
        # if invertRatio:
        #     ratioline = ratio.Clone("ratioline")
        #     ratioline.SetFillColor(0)
        #     ratioline.SetFillStyle(0)
        #     den.Draw("EPsame")
        # else: 
        #     ratio.Draw("EPsame")

        x1leg2 = 0.15 if leftMargin > 0.1 else 0.07
        x2leg2 = 0.55 if leftMargin > 0.1 else 0.32
        y1leg2 = 0.33 if leftMargin > 0.1 else 0.33
        y2leg2 = 0.38 if leftMargin > 0.1 else 0.4
        leg2 = ROOT.TLegend(x1leg2, y1leg2, x2leg2, y2leg2)
        #leg2 = ROOT.TLegend(0.07,0.30,0.27,0.35)
        leg2.SetFillColor(0)
        leg2.SetFillStyle(0)
        leg2.SetBorderSize(0)
        leg2.AddEntry(den,"combination","LF")
        #leg2.Draw("same")  # not sure I need to add this legend
        pad2.RedrawAxis("sameaxis")


    if draw_both0_noLog1_onlyLog2 != 2:
        canvas.SaveAs(outdir + canvasName + ".png")
        canvas.SaveAs(outdir + canvasName + ".pdf")

    if draw_both0_noLog1_onlyLog2 != 1:        
        if yAxisName == "a.u.": 
            hlep.GetYaxis().SetRangeUser(max(0.0001,hlep.GetMinimum()*0.8),hlep.GetMaximum()*100)
        else:
            hlep.GetYaxis().SetRangeUser(max(0.001,hlep.GetMinimum()*0.8),hlep.GetMaximum()*100)
        canvas.SetLogy()
        canvas.SaveAs(outdir + canvasName + "_logY.png")
        canvas.SaveAs(outdir + canvasName + "_logY.pdf")
        canvas.SetLogy(0)

    if saveRootFileForHepdata:
        rootfilename = outdir + canvasName + ".root"
        rfHep = ROOT.TFile.Open(rootfilename,"recreate")
        if not rfHep:
            print "Error in drawXsecAndTheoryband(): could not open root file %s" % rootfilename
            quit()
        rfHep.cd()
        hlep.Write("combination")
        hmu.Write("muon")
        hel.Write("electron")            
        rfHep.Close()
        
          
################################################################

def drawCheckTheoryBand(h1, h2, h3,
                        labelXtmp="xaxis", labelYtmp="yaxis",
                        canvasName="default", outdir="./",
                        #rebinFactorX=0,
                        draw_both0_noLog1_onlyLog2=0,                  
                        leftMargin=0.15,
                        rightMargin=0.04,
                        labelRatioTmp="rel. unc..::0.95,1.05",
                        drawStatBox=False,
                        #legendCoords="0.15,0.35,0.8,0.9",  # x1,x2,y1,y2
                        legendCoords="0.15,0.85,0.85,0.9;3",  # for unrolled use 1 line # x1,x2,y1,y2
                        canvasSize="1000,700",  # use X,Y to pass X and Y size     
                        lowerPanelHeight = 0.0,  # number from 0 to 1, 0.3 means 30% of space taken by lower panel. 0 means do not draw lower panel with relative error
                        passCanvas=None,
                        lumi=None,
                        drawVertLines="", # "12,36": format --> N of sections (e.g: 12 pt bins), and N of bins in each section (e.g. 36 eta bins), assuming uniform bin width
                        textForLines=[],                       
                        moreText="",
                        moreTextLatex="",
                        invertRatio = False,  # make expected over observed if True
                        useDifferenceInLowerPanel = False,
                        noLegendLowerPanel = False,
                        legendEntries = []
                    ):

    # moreText is used to pass some text to write somewhere (TPaveText is used)
    # e.g.  "stuff::x1,y1,x2,y2"  where xi and yi are the coordinates for the text
    # one can add more lines using the ";" key. FOr example, "stuff1;stuff2::x1,y1,x2,y2"
    # the coordinates should be defined taking into account how many lines will be drawn
    # if the coordinates are not passed (no "::"), then default ones are used, but this might not be satisfactory

    # moreTextLatex is similar, but used TLatex, and the four coordinates are x1,y1,ypass,textsize
    # where x1 and y1 are the coordinates the first line, and ypass is how much below y1 the second line is (and so on for following lines)

    # h1 is data, h2 in MC

    #if (rebinFactorX): 
    #    if isinstance(rebinFactorX, int): h1.Rebin(rebinFactorX)
    #    # case in which rebinFactorX is a list of bin edges
    #    else:                             h1.Rebin(len(rebinFactorX)-1,"",array('d',rebinFactorX)) 

    xAxisName,setXAxisRangeFromUser,xmin,xmax = getAxisRangeFromUser(labelXtmp)
    yAxisName,setYAxisRangeFromUser,ymin,ymax = getAxisRangeFromUser(labelYtmp)
    yRatioAxisName,setRatioYAxisRangeFromUser,yminRatio,ymaxRatio = getAxisRangeFromUser(labelRatioTmp)

    yAxisTitleOffset = 1.45 if leftMargin > 0.1 else 0.6

    addStringToEnd(outdir,"/",notAddIfEndswithMatch=True)
    createPlotDirAndCopyPhp(outdir)

    cw,ch = canvasSize.split(',')
    #canvas = ROOT.TCanvas("canvas",h2D.GetTitle() if plotLabel == "ForceTitle" else "",700,625)
    canvas = passCanvas if passCanvas != None else ROOT.TCanvas("canvas","",int(cw),int(ch))
    canvas.SetTickx(1)
    canvas.SetTicky(1)
    canvas.cd()
    canvas.SetLeftMargin(leftMargin)
    canvas.SetRightMargin(rightMargin)
    canvas.cd()

    pad2 = 0
    if lowerPanelHeight: 
        canvas.SetBottomMargin(lowerPanelHeight)
        pad2 = ROOT.TPad("pad2","pad2",0,0.,1,0.9)
        pad2.SetTopMargin(1-lowerPanelHeight)
        pad2.SetRightMargin(rightMargin)
        pad2.SetLeftMargin(leftMargin)
        pad2.SetFillColor(0)
        pad2.SetGridy(1)
        pad2.SetFillStyle(0)


    frame = h1.Clone("frame")
    frame.GetXaxis().SetLabelSize(0.04)
    frame.SetStats(0)

    #ymax = max(ymax, max(h1.GetBinContent(i)+h1.GetBinError(i) for i in range(1,h1.GetNbinsX()+1)))
    # if min and max were not set, set them based on histogram content
    if ymin == ymax == 0.0:
        ymin,ymax = getMinMaxHisto(h1,excludeEmpty=True,sumError=True)            
        ymin *= 0.9
        ymax *= (1.1 if leftMargin > 0.1 else 2.0)
        if ymin < 0: ymin = 0
        #print "drawSingleTH1() >>> Histo: %s     minY,maxY = %.2f, %.2f" % (h1.GetName(),ymin,ymax)

    # print "#### WARNING ####"
    # print "Hardcoding ymin = 0 in function drawDataAndMC(): change it if it is not what you need"
    # print "#################"
    # ymin = 0 # hardcoded

    if lowerPanelHeight:
        h1.GetXaxis().SetLabelSize(0)
        h1.GetXaxis().SetTitle("")  
    else:
        h1.GetXaxis().SetTitle(xAxisName)
        h1.GetXaxis().SetTitleOffset(1.2)
        h1.GetXaxis().SetTitleSize(0.05)
        h1.GetXaxis().SetLabelSize(0.04)
    h1.GetYaxis().SetTitle(yAxisName)
    h1.GetYaxis().SetTitleOffset(yAxisTitleOffset)
    h1.GetYaxis().SetTitleSize(0.05)
    h1.GetYaxis().SetLabelSize(0.04)
    h1.GetYaxis().SetRangeUser(ymin, ymax)    
    h1.GetYaxis().SetTickSize(0.01)
    if setXAxisRangeFromUser: h1.GetXaxis().SetRangeUser(xmin,xmax)


    ratio = None
    den_noerr = None
    den = None
    ratio1 = h1.Clone("ratio1")
    ratio2 = h2.Clone("ratio2")
    ratio3 = h3.Clone("ratio3")
    den_noerr = h1.Clone("den_noerr")

    for iBin in range (1,den_noerr.GetNbinsX()+1):
        den_noerr.SetBinError(iBin,0.)

    h2IsGraph = False
    h3IsGraph = False
    ratio1.Divide(den_noerr)

    if ratio2.InheritsFrom("TH1"):
        ratio2.Divide(den_noerr)
    else:
        h2IsGraph = True
        for i in range(den_noerr.GetNbinsX()):
            ratio2.SetPoint(i, den_noerr.GetBinCenter(i+1), ratio2.Eval(den_noerr.GetBinCenter(i+1)) / h1.GetBinContent(i+1) )
            ratio2.SetPointEYhigh(i, ratio2.GetErrorYhigh(i) / h1.GetBinContent(i+1))
            ratio2.SetPointEYlow( i, ratio2.GetErrorYlow(i)  / h1.GetBinContent(i+1))

    if ratio3.InheritsFrom("TH1"):
        ratio3.Divide(den_noerr)
    else:
        h3IsGraph = True
        for i in range(den_noerr.GetNbinsX()):
            ratio3.SetPoint(i, den_noerr.GetBinCenter(i+1), ratio3.Eval(den_noerr.GetBinCenter(i+1)) / h1.GetBinContent(i+1) )
            ratio3.SetPointEYhigh(i, ratio3.GetErrorYhigh(i) / h1.GetBinContent(i+1))
            ratio3.SetPointEYlow( i, ratio3.GetErrorYlow(i)  / h1.GetBinContent(i+1))


    ratio1.SetFillColor(ROOT.kGreen) # kGreen+1
    ratio1.SetFillStyle(3001)  # 1001 to make it solid again
    ratio1.SetLineColor(ROOT.kGreen+1) # kGreen+2                       
    ratio1.SetLineWidth(1) # make it smaller when it is drawn on top of something
    ratio2.SetFillColor(ROOT.kRed+1) # kGreen+1
    ratio2.SetFillStyle(3244)  # 1001 to make it solid again
    ratio2.SetLineColor(ROOT.kRed+2) # kGreen+2                       
    ratio2.SetLineWidth(1) # make it smaller when it is drawn on top of something
    ratio3.SetLineColor(ROOT.kBlack)
    ratio3.SetMarkerColor(ROOT.kBlack)
    ratio3.SetMarkerStyle(20)
    ratio3.SetMarkerSize(1)
    
    ratio1.Draw("E2")
    ratio2.Draw("F2 SAME" if h2IsGraph else "E2 SAME")
    ratio3.Draw("EP SAME")
    line = ROOT.TF1("horiz_line","1",
                    ratio1.GetXaxis().GetBinLowEdge(1),ratio1.GetXaxis().GetBinLowEdge(ratio1.GetNbinsX()+1))
    line.SetLineColor(ROOT.kRed+3)
    line.SetLineWidth(1)
    line.Draw("Lsame")

    nColumnsLeg = 1
    if ";" in legendCoords: 
        nColumnsLeg = int(legendCoords.split(";")[1])
    legcoords = [float(x) for x in (legendCoords.split(";")[0]).split(',')]
    lx1,lx2,ly1,ly2 = legcoords[0],legcoords[1],legcoords[2],legcoords[3]
    leg = ROOT.TLegend(lx1,ly1,lx2,ly2)
    #leg.SetFillColor(0)
    #leg.SetFillStyle(0)
    #leg.SetBorderSize(0)
    leg.SetNColumns(3)
    leg.AddEntry(ratio1,"PDFs","LF")
    leg.AddEntry(ratio2,"#alpha_{S}","LF")        
    leg.AddEntry(ratio3,"QCD scales","LPE")        
    leg.Draw("same")
    canvas.RedrawAxis("sameaxis")

    if drawStatBox:
        ROOT.gPad.Update()
        ROOT.gStyle.SetOptStat(1110)
        ROOT.gStyle.SetOptFit(1102)
    else:
        h1.SetStats(0)

    vertline = ROOT.TLine(36,0,36,canvas.GetUymax())
    vertline.SetLineColor(ROOT.kBlack)
    vertline.SetLineStyle(3) # 2 for denser hatches
    bintext = ROOT.TLatex()
    #bintext.SetNDC()
    bintext.SetTextSize(0.025)  # 0.03
    bintext.SetTextFont(42)
    if len(textForLines):
        bintext.SetTextAngle(45 if "#eta" in textForLines[0] else 30)

    if len(drawVertLines):
        #print "drawVertLines = " + drawVertLines
        nptBins = int(drawVertLines.split(',')[0])
        etarange = float(drawVertLines.split(',')[1])        
        offsetXaxisHist = h1.GetXaxis().GetBinLowEdge(0)
        sliceLabelOffset = 6. 
        for i in range(1,nptBins): # do not need line at canvas borders
            #vertline.DrawLine(offsetXaxisHist+etarange*i,0,offsetXaxisHist+etarange*i,canvas.GetUymax())
            vertline.DrawLine(etarange*i-offsetXaxisHist,0,etarange*i-offsetXaxisHist,ymax)
        if len(textForLines):
            for i in range(0,len(textForLines)): # we need nptBins texts
                #texoffset = 0.1 * (4 - (i%4))
                #ytext = (1. + texoffset)*ymax/2.  
                ytext = 0.6*ymax                  
                bintext.DrawLatex(etarange*i + etarange/sliceLabelOffset, ytext, textForLines[i])

    # redraw legend, or vertical lines appear on top of it
    leg.Draw("same")

    if len(moreText):
        realtext = moreText.split("::")[0]
        x1,y1,x2,y2 = 0.7,0.8,0.9,0.9
        if "::" in moreText:
            x1,y1,x2,y2 = (float(x) for x in (moreText.split("::")[1]).split(","))
        pavetext = ROOT.TPaveText(x1,y1,x2,y2,"NB NDC")
        for tx in realtext.split(";"):
            pavetext.AddText(tx)
        pavetext.SetFillColor(0)
        pavetext.SetFillStyle(0)
        pavetext.SetBorderSize(0)
        pavetext.SetLineColor(0)
        pavetext.Draw("same")

    if len(moreTextLatex):
        realtext = moreTextLatex.split("::")[0]
        x1,y1,ypass,textsize = 0.75,0.8,0.08,0.035
        if "::" in moreTextLatex:
            x1,y1,ypass,textsize = (float(x) for x in (moreTextLatex.split("::")[1]).split(","))            
        lat = ROOT.TLatex()
        lat.SetNDC();
        lat.SetTextFont(42)        
        lat.SetTextSize(textsize)
        for itx,tx in enumerate(realtext.split(";")):
            lat.DrawLatex(x1,y1-itx*ypass,tx)

  # TPaveText *pvtxt = NULL;
  # if (yAxisName == "a.u.") {
  #   pvtxt = new TPaveText(0.6,0.6,0.95,0.7, "BR NDC")
  #   pvtxt.SetFillColor(0)
  #   pvtxt.SetFillStyle(0)
  #   pvtxt.SetBorderSize(0)
  #   pvtxt.AddText(Form("norm num/den = %.2f +/- %.2f",IntegralRatio,ratioError))
  #   pvtxt.Draw()
  # }

    setTDRStyle()
    if leftMargin > 0.1:
        if lumi != None: CMS_lumi(canvas,lumi,True,False)
        else:            CMS_lumi(canvas,"",True,False)
    else:
        latCMS = ROOT.TLatex()
        latCMS.SetNDC();
        latCMS.SetTextFont(42)
        latCMS.SetTextSize(0.045)
        latCMS.DrawLatex(0.1, 0.95, '#bf{CMS} #it{Preliminary}')
        if lumi != None: latCMS.DrawLatex(0.85, 0.95, '%s fb^{-1} (13 TeV)' % lumi)
        else:            latCMS.DrawLatex(0.90, 0.95, '(13 TeV)' % lumi)

    if lowerPanelHeight:
        pad2.Draw()
        pad2.cd()

        frame.Reset("ICES")
        if setRatioYAxisRangeFromUser: frame.GetYaxis().SetRangeUser(yminRatio,ymaxRatio)
        #else:                          
        #frame.GetYaxis().SetRangeUser(0.5,1.5)
        frame.GetYaxis().SetNdivisions(5)
        frame.GetYaxis().SetTitle(yRatioAxisName)
        frame.GetYaxis().SetTitleOffset(yAxisTitleOffset)
        frame.GetYaxis().SetTitleSize(0.05)
        frame.GetYaxis().SetLabelSize(0.04)
        frame.GetYaxis().CenterTitle()
        frame.GetXaxis().SetTitle(xAxisName)
        if setXAxisRangeFromUser: frame.GetXaxis().SetRangeUser(xmin,xmax)
        frame.GetXaxis().SetTitleOffset(1.2)
        frame.GetXaxis().SetTitleSize(0.05)

        #ratio = copy.deepcopy(h1.Clone("ratio"))
        #den_noerr = copy.deepcopy(h2.Clone("den_noerr"))
        ratio = None
        den_noerr = None
        den = None
        if invertRatio:
            ratio = h2.Clone("ratio")
            den_noerr = h1.Clone("den_noerr")
            den = h1.Clone("den")
        else:
            ratio = h1.Clone("ratio")
            den_noerr = h2.Clone("den_noerr")
            den = h2.Clone("den")

        for iBin in range (1,den_noerr.GetNbinsX()+1):
            den_noerr.SetBinError(iBin,0.)

        if useDifferenceInLowerPanel:
            ratio.Add(den_noerr,-1.0)
            den.Add(den_noerr,-1.0)
        else:
            ratio.Divide(den_noerr)
            den.Divide(den_noerr)

        if invertRatio:
            if histMCpartialUnc == None:
                ratio.SetFillColor(ROOT.kGray+1)
                ratio.SetFillStyle(1001)  # make it solid again
                ratio.SetLineColor(ROOT.kGray+3)
            else:
                ratio.SetFillColor(ROOT.kRed+1) # kGreen+1
                ratio.SetFillStyle(1001)  # 1001 to make it solid again
                ratio.SetLineColor(ROOT.kRed+2) # kGreen+2                       
                ratio.SetLineWidth(1) # make it smaller when it is drawn on top of something
        else:
            if histMCpartialUnc == None:
                den.SetFillColor(ROOT.kGray+1)
                den.SetFillStyle(1001)  # make it solid again
                den.SetLineColor(ROOT.kRed)        
        if histMCpartialUnc != None:
            h3ratio = h3.Clone("h3ratio")
            if useDifferenceInLowerPanel:
                h3ratio.Add(den_noerr,-1.0)
            else:
                h3ratio.Divide(den_noerr)
            #h3ratio.SetFillStyle(3144) # 1001 for solid, 3144 instead of 3244, to have more dense hatches
            h3ratio.SetFillStyle(1001) # 
            h3ratio.SetFillColor(ROOT.kGreen+1)  # kRed-4
            h3ratio.SetLineColor(ROOT.kGreen+2)  # kRed-4

        frame.Draw()        
        if invertRatio:
            den.SetMarkerSize(0.85)
            den.SetMarkerStyle(20) 
            #den.Draw("EPsame")    # draw after red line, not now
            ratio.Draw("E2same")
        else:    
            ratio.SetMarkerSize(0.85)
            ratio.SetMarkerStyle(20) 
            den.Draw("E2same")
            #ratio.Draw("EPsame") # draw after red line, not now

        if histMCpartialUnc != None:
            h3ratio.Draw("E2 same")
        
        line = ROOT.TF1("horiz_line","0" if useDifferenceInLowerPanel else "1",
                        ratio.GetXaxis().GetBinLowEdge(1),ratio.GetXaxis().GetBinLowEdge(ratio.GetNbinsX()+1))
        line.SetLineWidth(1)
        line.Draw("Lsame")
        if invertRatio:
            ratioline = ratio.Clone("ratioline")
            ratioline.SetFillColor(0)
            if histMCpartialUnc != None and len(drawVertLines):
                # unrolled, here the red line would be too large and given the fluctuations is could be seen as an additional band
                ratioline.SetLineColor(histMCpartialUnc.GetLineColor())
            ratioline.SetFillStyle(0)
            if histMCpartialUnc != None: ratioline.Draw("HIST same") # to draw the line inside the band for the expected
            den.Draw("EPsame")
        else: 
            ratio.Draw("EPsame")

        x1leg2 = 0.15 if leftMargin > 0.1 else 0.07
        x2leg2 = 0.55 if leftMargin > 0.1 else 0.27
        y1leg2 = 0.28 if leftMargin > 0.1 else 0.3
        y2leg2 = 0.33 if leftMargin > 0.1 else 0.35
        leg2 = ROOT.TLegend(x1leg2, y1leg2, x2leg2, y2leg2)
        #leg2 = ROOT.TLegend(0.07,0.30,0.27,0.35)
        leg2.SetFillColor(0)
        leg2.SetFillStyle(0)
        leg2.SetBorderSize(0)
        if invertRatio:
            leg2.AddEntry(ratio,"total theory uncertainty","LF")
        else:
            leg2.AddEntry(den,"total theory uncertainty","LF")
        if not noLegendLowerPanel: leg2.Draw("same")
        if histMCpartialUnc != None:
            leg2bis = ROOT.TLegend(x1leg2 + 0.45, y1leg2, x2leg2 + 0.45, y2leg2)
            leg2bis.SetFillColor(0)
            leg2bis.SetFillStyle(0)
            leg2bis.SetBorderSize(0)
            leg2bis.AddEntry(h3ratio,histMCpartialUncLegEntry,"LF")
            if not noLegendLowerPanel: leg2bis.Draw("same")

        pad2.RedrawAxis("sameaxis")


    if draw_both0_noLog1_onlyLog2 != 2:
        canvas.SaveAs(outdir + canvasName + ".png")
        canvas.SaveAs(outdir + canvasName + ".pdf")

    if draw_both0_noLog1_onlyLog2 != 1:        
        if yAxisName == "a.u.": 
            h1.GetYaxis().SetRangeUser(max(0.0001,h1.GetMinimum()*0.8),h1.GetMaximum()*100)
        else:
            h1.GetYaxis().SetRangeUser(max(0.001,h1.GetMinimum()*0.8),h1.GetMaximum()*100)
        canvas.SetLogy()
        canvas.SaveAs(outdir + canvasName + "_logY.png")
        canvas.SaveAs(outdir + canvasName + "_logY.pdf")
        canvas.SetLogy(0)

#########################################################################

# warning: h1 is histogram, h2 and histMCpartialUnc is a graph, there was a reason for that, but I do not remember

def drawXsecAndTheoryband(h1, h2,  # h1 is data, h2 is total uncertainty band
                          labelXtmp="xaxis", labelYtmp="yaxis",
                          canvasName="default", outdir="./",
                          #rebinFactorX=0,
                          draw_both0_noLog1_onlyLog2=0,                  
                          leftMargin=0.15,
                          rightMargin=0.04,
                          labelRatioTmp="Data/pred.::0.5,1.5",
                          drawStatBox=False,
                          #legendCoords="0.15,0.35,0.8,0.9",  # x1,x2,y1,y2
                          legendCoords="0.20,0.65,0.86,0.93;3",  # for unrolled use 1 line # x1,x2,y1,y2
                          canvasSize="600,700",  # use X,Y to pass X and Y size     
                          lowerPanelHeight = 0.3,  # number from 0 to 1, 0.3 means 30% of space taken by lower panel. 0 means do not draw lower panel with relative error
                          passCanvas=None,
                          lumi=None,
                          drawVertLines="", # "12,36": format --> N of sections (e.g: 12 pt bins), and N of bins in each section (e.g. 36 eta bins), assuming uniform bin width
                          textForLines=[],                       
                          moreText="",
                          moreTextLatex="",
                          invertRatio = False,  # make expected over observed if True
                          histMCpartialUnc = None,
                          histMCpartialUncLegEntry = "",
                          useDifferenceInLowerPanel = False,
                          noLegendLowerPanel = True,
                          legendEntries = [],
                          skipPreliminary = True,
                          nSplitUnrolledBins = 1,  # integer telling in how many parts the unrolled must be split
                          saveRootFileForHepdata = True
                      ):

    if h1.GetNbinsX()%nSplitUnrolledBins:
        print "="*30
        print "Error in drawXsecAndTheoryband(): number of histogram bins is not a multiple of the splitting number."
        print "ABORT"
        print "="*30
        quit()
    
    # moreText is used to pass some text to write somewhere (TPaveText is used)
    # e.g.  "stuff::x1,y1,x2,y2"  where xi and yi are the coordinates for the text
    # one can add more lines using the ";" key. FOr example, "stuff1;stuff2::x1,y1,x2,y2"
    # the coordinates should be defined taking into account how many lines will be drawn
    # if the coordinates are not passed (no "::"), then default ones are used, but this might not be satisfactory

    # moreTextLatex is similar, but used TLatex, and the four coordinates are x1,y1,ypass,textsize
    # where x1 and y1 are the coordinates the first line, and ypass is how much below y1 the second line is (and so on for following lines)

    # h1 is data, h2 in MC

    #if (rebinFactorX): 
    #    if isinstance(rebinFactorX, int): h1.Rebin(rebinFactorX)
    #    # case in which rebinFactorX is a list of bin edges
    #    else:                             h1.Rebin(len(rebinFactorX)-1,"",array('d',rebinFactorX)) 

    # colorBandpart = {"line" : ROOT.kRed,
    #                  "fill" : ROOT.kRed-9}
    # colorBandTot = {"line" : ROOT.kGreen,
    #                 "fill" : ROOT.kGreen-9}
    # colorBandPart = {"line" : ROOT.kCyan+2,
    #                  "fill" : ROOT.kCyan}
    # colorBandTot = {"line" : ROOT.kOrange+2,
    #                 "fill" : ROOT.kOrange}
    # colorBandPart = {"line" : ROOT.kCyan+2,
    #                  "fill" : ROOT.kCyan-9} # was -7
    # colorBandTot = {"line" : ROOT.kOrange+7,
    #                 "fill" : ROOT.kOrange+6} # try -3, -9, +1, -4, -2, +6
    # colorBandPart = {"line" : ROOT.kAzure+7,
    #                  "fill" : ROOT.kAzure+6} # try -9
    # colorBandTot = {"line" : ROOT.kOrange+10,
    #                 "fill" : ROOT.kOrange-9} # + 6 # try -3, -9, +1, -4, -2, +6

    # colorBandPart = {"line" : ROOT.kAzure+7,
    #                  "fill" : ROOT.kAzure+7} # try +6 and 50% transparent
    # colorBandTot = {"line" : ROOT.kOrange+10,
    #                 "fill" : ROOT.kOrange+7} # try kRed+1 with 50% transparent #  + 6 # try -3, -9, +1, -4, -2, +6
    # fcalphaPart = 0.5
    # fcalphaTot = 0.9
    # fsPart = 1001
    # fsTot = 1001

    colorBandPart = {"line" : ROOT.kBlue,
                     "fill" : ROOT.kAzure+7} # try +6 and 50% transparent
    colorBandTot = {"line" : ROOT.kPink-6,
                    "fill" : ROOT.kPink-6} # try kRed+1 with 50% transparent #  + 6 # try -3, -9, +1, -4, -2, +6
    centralLineColor = colorBandPart["line"]
    linewidth = 2 if leftMargin > 0.1 else 1
    fcalphaPart = 0.6 # 0.5
    fcalphaTot = 0.66 # 0.9
    fsPart = 1001
    fsTot = 1001

    yAxisTitleOffset = 1.45 if leftMargin > 0.1 else 0.63

    addStringToEnd(outdir,"/",notAddIfEndswithMatch=True)
    createPlotDirAndCopyPhp(outdir)

    for isplit in range(nSplitUnrolledBins):

        xAxisName,setXAxisRangeFromUser,xmin,xmax = getAxisRangeFromUser(labelXtmp)
        yAxisName,setYAxisRangeFromUser,ymin,ymax = getAxisRangeFromUser(labelYtmp)
        yRatioAxisName,setRatioYAxisRangeFromUser,yminRatio,ymaxRatio = getAxisRangeFromUser(labelRatioTmp)

        cw,ch = canvasSize.split(',')
        #canvas = ROOT.TCanvas("canvas",h2D.GetTitle() if plotLabel == "ForceTitle" else "",700,625)
        canvas = passCanvas if passCanvas != None else ROOT.TCanvas("canvas","",int(cw),int(ch))
        canvas.SetTickx(1)
        canvas.SetTicky(1)
        canvas.cd()
        canvas.SetLeftMargin(leftMargin)
        canvas.SetRightMargin(rightMargin)
        canvas.cd()

        pad2 = 0
        if lowerPanelHeight: 
            canvas.SetBottomMargin(lowerPanelHeight)
            pad2 = ROOT.TPad("pad2","pad2",0,0.,1,0.9)
            pad2.SetTopMargin(1-lowerPanelHeight)
            pad2.SetRightMargin(rightMargin)
            pad2.SetLeftMargin(leftMargin)
            pad2.SetFillColor(0)
            pad2.SetGridy(1)
            pad2.SetFillStyle(0)

        h1_true = None
        h2_true = None
        histMCpartialUnc_true = None
        textForLines_true = []
        canvasName_true = ""
        nbinsSplit = h1.GetNbinsX()
        if nSplitUnrolledBins == 1:
            h1_true = h1
            h2_true = h2
            histMCpartialUnc_true = histMCpartialUnc
            textForLines_true = textForLines
            canvasName_true = canvasName
        else:
            # start to split things
            ntxtLines = len(textForLines)/nSplitUnrolledBins
            for islice in range(ntxtLines):
                textForLines_true.append(textForLines[isplit*ntxtLines+islice])
            canvasName_true = canvasName + "_part%d"%(isplit+1)
            nbinsSplit = h1.GetNbinsX()/nSplitUnrolledBins
            h1_true = ROOT.TH1F(h1.GetName()+"_part%d"%(isplit+1),h1.GetTitle(),
                                nbinsSplit,0.5+isplit*nbinsSplit,(isplit+1)*nbinsSplit+0.5)
            for ibin in range(1,1+nbinsSplit):                
                h1bin = ibin+isplit*nbinsSplit
                h1_true.SetBinContent(ibin, h1.GetBinContent(h1bin))
                h1_true.SetBinError(  ibin, h1.GetBinError(h1bin))
            ##
            ## h2_true = ROOT.TH1F(h2.GetName()+"_part%d"%(isplit+1),h2.GetTitle(),
            ##                     nbinsSplit,0.5,nbinsSplit+0.5)
            ## histMCpartialUnc_true = ROOT.TH1F(histMCpartialUnc.GetName()+"_part%d"%(isplit+1),
            ##                                   histMCpartialUnc.GetTitle(),
            ##                                   nbinsSplit,0.5,nbinsSplit+0.5)

            # the following are graphs
            h2_true = ROOT.TGraphAsymmErrors(nbinsSplit)
            h2_true.SetName(h2.GetName()+"_part%d"%(isplit+1))
            h2_true.SetTitle(h2.GetTitle())
            for ip in range(nbinsSplit):
                ibin = ip + isplit * nbinsSplit
                hx = ROOT.Double(0.0)
                hy = ROOT.Double(0.0)
                h2.GetPoint(ibin,hx,hy)
                h2_true.SetPoint(ip,hx,hy)
                h2_true.SetPointError(ip,
                                      h2.GetErrorXlow(ibin),h2.GetErrorXhigh(ibin),
                                      h2.GetErrorYlow(ibin),h2.GetErrorYhigh(ibin))
            #
            histMCpartialUnc_true = ROOT.TGraphAsymmErrors(nbinsSplit)
            histMCpartialUnc_true.SetName(histMCpartialUnc.GetName()+"_part%d"%(isplit+1))
            histMCpartialUnc_true.SetTitle(histMCpartialUnc.GetTitle())
            for ip in range(nbinsSplit):
                ibin = ip + isplit * nbinsSplit
                hx = ROOT.Double(0.0)
                hy = ROOT.Double(0.0)
                histMCpartialUnc.GetPoint(ibin,hx,hy)
                histMCpartialUnc_true.SetPoint(ip,hx,hy)
                histMCpartialUnc_true.SetPointError(ip,
                                                    histMCpartialUnc.GetErrorXlow(ibin),histMCpartialUnc.GetErrorXhigh(ibin),
                                                    histMCpartialUnc.GetErrorYlow(ibin),histMCpartialUnc.GetErrorYhigh(ibin))


        frame = h1_true.Clone("frame")
        frame.GetXaxis().SetLabelSize(0.04 if leftMargin > 0.1 else 0.05)
        frame.SetStats(0)

        h1_true.SetLineColor(ROOT.kBlack)
        h1_true.SetMarkerColor(ROOT.kBlack)
        h1_true.SetMarkerStyle(20)
        h1_true.SetMarkerSize(1.2)

        #ymax = max(ymax, max(h1_true.GetBinContent(i)+h1_true.GetBinError(i) for i in range(1,h1_true.GetNbinsX()+1)))
        # if min and max were not set, set them based on histogram content
        if ymin == ymax == 0.0:
            ymin,ymax = getMinMaxHisto(h1_true,excludeEmpty=True,sumError=True)            
            ymin *= 0.9
            ymax *= (1.1 if leftMargin > 0.1 else 2.0)
            if ymin < 0: ymin = 0
            #print "drawSingleTH1() >>> Histo: %s     minY,maxY = %.2f, %.2f" % (h1.GetName(),ymin,ymax)

        # print "#### WARNING ####"
        # print "Hardcoding ymin = 0 in function drawDataAndMC(): change it if it is not what you need"
        # print "#################"
        # ymin = 0 # hardcoded

        if lowerPanelHeight:
            h1_true.GetXaxis().SetLabelSize(0)
            h1_true.GetXaxis().SetTitle("")  
        else:
            h1_true.GetXaxis().SetTitle(xAxisName)
            h1_true.GetXaxis().SetTitleOffset(1.2)
            h1_true.GetXaxis().SetTitleSize(0.05)
            h1_true.GetXaxis().SetLabelSize(0.04)
            h1_true.GetXaxis().SetLabelSize(0.05)
        h1_true.GetYaxis().SetTitle(yAxisName)
        h1_true.GetYaxis().SetTitleOffset(yAxisTitleOffset)
        h1_true.GetYaxis().SetTitleSize(0.05)
        h1_true.GetYaxis().SetLabelSize(0.04)
        h1_true.GetYaxis().SetDecimals()
        h1_true.GetYaxis().SetRangeUser(ymin, ymax)    
        h1_true.GetYaxis().SetTickSize(0.01)
        if setXAxisRangeFromUser: h1_true.GetXaxis().SetRangeUser(xmin,xmax)
        h1_true.Draw("EP0")
        #h1err = h1.Clone("h1err")
        #h1err.SetFillColor(ROOT.kRed+2)
        h2_true.SetFillStyle(fsTot)  # 1001 for solid, 3001 is better than 3002 for pdf, while 3002 is perfect for png
        h2_true.SetLineColor(colorBandTot["line"])  # kGreen+2
        h2_true.SetFillColor(colorBandTot["fill"])    # kGreen
        h2_true.SetFillColorAlpha(colorBandTot["fill"], fcalphaTot);
        h2_true.SetLineWidth(linewidth)
        h2_true.Draw("E2 SAME")
        h3 = None
        if histMCpartialUnc_true != None:
            h3 = histMCpartialUnc_true.Clone("histMCpartialUnc_true")
            h3dummyWhite = histMCpartialUnc_true.Clone("h3dummyWhite")
            h3dummyWhite.SetFillColor(ROOT.kWhite)
            h3dummyWhite.SetLineColor(ROOT.kWhite)
            h3dummyWhite.SetFillStyle(1001)
            h3dummyWhite.Draw("E2 SAME")
            h3.SetFillColor(colorBandPart["fill"])
            h3.SetLineColor(colorBandPart["line"])
            h3.SetFillColorAlpha(colorBandPart["fill"], fcalphaPart);
            h3.SetLineWidth(linewidth)
            h3.SetFillStyle(fsPart)  # 1001, 3001, 3144 , 3244, 3003
            #h3.SetFillStyle(3244)  # 3144 , 3244, 3003
            h3.Draw("E2 SAME")
            #for i in range(1,1+h3.GetNbinsX()):
            #    print "PDF band: bin %d  val +/- error = %.3f +/- %.3f" % (i, h3.GetBinContent(i),h3.GetBinError(i))
        h2line = None
        h2line = h1_true.Clone("h2line")
        for i in range(1,h2line.GetNbinsX()+1):
            xval = h2line.GetBinCenter(i)
            yval = h2_true.Eval(xval)
            h2line.SetBinContent(i,yval)
        h2line.SetFillColor(0)
        #h2line.SetLineColor(h2.GetLineColor())
        h2line.SetLineColor(centralLineColor)
        h2line.SetLineWidth(linewidth)
        h2line.Draw("HIST SAME")
        h1_true.Draw("EP0 SAME")

        nColumnsLeg = 1
        if ";" in legendCoords: 
            nColumnsLeg = int(legendCoords.split(";")[1])
            #if histMCpartialUnc_true != None: nColumnsLeg = 3
        legcoords = [float(x) for x in (legendCoords.split(";")[0]).split(',')]
        lx1,lx2,ly1,ly2 = legcoords[0],legcoords[1],legcoords[2],legcoords[3]
        if leftMargin > 0.1:
            if histMCpartialUnc_true != None:
                ly2 = ly2 + 0.5 * (ly2 - ly1) # add one more row (for unrolled everything is on the same line
        else:
            if len(moreText):
                realtext = moreText.split("::")[0]
                lx2 = lx2 + 0.333 * (lx2 - lx1) # increase length of legend
                nColumnsLeg += 1

        leg = ROOT.TLegend(lx1,ly1,lx2,ly2)
        #leg.SetFillColor(0)
        #leg.SetFillStyle(0)
        #leg.SetBorderSize(0)
        leg.SetNColumns(nColumnsLeg)
        if len(legendEntries):
            leg.AddEntry(h1_true,str(legendEntries[0]),"LPE")
            leg.AddEntry(h2_true,str(legendEntries[1]),"LF")        
        else:
            leg.AddEntry(h1_true,"measured","LPE")
            leg.AddEntry(h2_true,"MadGraph5_aMC@NLO","LF")
        if histMCpartialUnc_true != None:
            leg.AddEntry(h3,histMCpartialUncLegEntry,"LF")
        if len(moreText):
            leg.AddEntry(None,realtext+ "  ","")
        #leg.AddEntry(h1err,"Uncertainty","LF")
        leg.Draw("same")
        canvas.RedrawAxis("sameaxis")

        if drawStatBox:
            ROOT.gPad.Update()
            ROOT.gStyle.SetOptStat(1110)
            ROOT.gStyle.SetOptFit(1102)
        else:
            h1_true.SetStats(0)

        vertline = ROOT.TLine(36,0,36,canvas.GetUymax())
        vertline.SetLineColor(ROOT.kBlack)
        vertline.SetLineStyle(3) # 2 for denser hatches
        bintext = ROOT.TLatex()
        #bintext.SetNDC()
        textsize = {1: 0.038,
                    2: 0.04,
                    3: 0.05}
        bintext.SetTextSize(textsize[nSplitUnrolledBins])  # 0.03
        bintext.SetTextFont(42)
        if len(textForLines_true):
            if nSplitUnrolledBins > 2:
                pass
                #bintext.SetTextAngle(10)
            elif nSplitUnrolledBins > 1:
                bintext.SetTextAngle(15)
            else:
                #bintext.SetTextAngle(45 if "#eta" in textForLines_true[0] else 30)
                bintext.SetTextAngle(18)

        if len(drawVertLines):
            #print "drawVertLines = " + drawVertLines
            nptBins = int(drawVertLines.split(',')[0])/nSplitUnrolledBins
            etarange = float(drawVertLines.split(',')[1])        
            offsetXaxisHist = h1_true.GetXaxis().GetBinLowEdge(0)
            sliceLabelOffset = 6. 
            globalbinOffset = isplit * nbinsSplit
            for i in range(1,nptBins): # do not need line at canvas borders
                #vertline.DrawLine(offsetXaxisHist+etarange*i,0,offsetXaxisHist+etarange*i,canvas.GetUymax())
                vertline.DrawLine(etarange*i-offsetXaxisHist, 0, 
                                  etarange*i-offsetXaxisHist, ymax)
            if len(textForLines_true):
                for i in range(0,len(textForLines_true)): # we need nptBins texts
                    #texoffset = 0.1 * (4 - (i%4))
                    #ytext = (1. + texoffset)*ymax/2.  
                    ytext = 0.67*ymax  # was 0.6                
                    bintext.DrawLatex(globalbinOffset + etarange*i + etarange/sliceLabelOffset, 
                                      ytext, 
                                      textForLines_true[i])

        # redraw legend, or vertical lines appear on top of it
        leg.Draw("same")

        if len(moreText) and leftMargin > 0.1:
            realtext = moreText.split("::")[0]
            x1,y1,x2,y2 = 0.7,0.8,0.9,0.9
            if "::" in moreText:
                x1,y1,x2,y2 = (float(x) for x in (moreText.split("::")[1]).split(","))
            pavetext = ROOT.TPaveText(x1,y1,x2,y2,"NB NDC")
            for tx in realtext.split(";"):
                pavetext.AddText(tx)
            pavetext.SetFillColor(ROOT.kWhite)
            #pavetext.SetFillStyle(0)
            pavetext.SetBorderSize(0)
            pavetext.SetLineColor(0)
            pavetext.Draw("same")

        if len(moreTextLatex):
            realtext = moreTextLatex.split("::")[0]
            x1,y1,ypass,textsize = 0.75,0.8,0.08,0.035
            if "::" in moreTextLatex:
                x1,y1,ypass,textsize = (float(x) for x in (moreTextLatex.split("::")[1]).split(","))            
            lat = ROOT.TLatex()
            lat.SetNDC();
            lat.SetTextFont(42)        
            lat.SetTextSize(textsize)
            for itx,tx in enumerate(realtext.split(";")):
                lat.DrawLatex(x1,y1-itx*ypass,tx)

      # TPaveText *pvtxt = NULL;
      # if (yAxisName == "a.u.") {
      #   pvtxt = new TPaveText(0.6,0.6,0.95,0.7, "BR NDC")
      #   pvtxt.SetFillColor(0)
      #   pvtxt.SetFillStyle(0)
      #   pvtxt.SetBorderSize(0)
      #   pvtxt.AddText(Form("norm num/den = %.2f +/- %.2f",IntegralRatio,ratioError))
      #   pvtxt.Draw()
      # }

        setTDRStyle()
        if leftMargin > 0.1:
            if lumi != None: CMS_lumi(canvas,lumi,True,skipPreliminary, offsetLumi=0.02)
            else:            CMS_lumi(canvas,"",True,skipPreliminary)
        else:
            latCMS = ROOT.TLatex()
            latCMS.SetNDC();
            latCMS.SetTextFont(42)
            latCMS.SetTextSize(0.045)
            latCMS.DrawLatex(leftMargin, 0.95, '#bf{CMS}' + (' #it{Preliminary}' if not skipPreliminary else ''))
            if lumi != None: latCMS.DrawLatex(0.85+(0.04-rightMargin), 0.95, '%s fb^{-1} (13 TeV)' % lumi)
            else:            latCMS.DrawLatex(0.90+(0.04-rightMargin), 0.95, '(13 TeV)' % lumi)

        if lowerPanelHeight:
            pad2.Draw()
            pad2.cd()

            frame.Reset("ICES")
            if setRatioYAxisRangeFromUser: frame.GetYaxis().SetRangeUser(yminRatio,ymaxRatio)
            #else:                          
            #frame.GetYaxis().SetRangeUser(0.5,1.5)
            frame.GetYaxis().SetNdivisions(5)
            frame.GetYaxis().SetTitle(yRatioAxisName)
            frame.GetYaxis().SetTitleOffset(yAxisTitleOffset)
            frame.GetYaxis().SetTitleSize(0.05)
            frame.GetYaxis().SetLabelSize(0.04)
            frame.GetYaxis().SetDecimals()
            frame.GetYaxis().CenterTitle()
            frame.GetXaxis().SetTitle(xAxisName)
            if setXAxisRangeFromUser: frame.GetXaxis().SetRangeUser(xmin,xmax)
            frame.GetXaxis().SetTitleOffset(1.2)
            frame.GetXaxis().SetTitleSize(0.05)

            #ratio = copy.deepcopy(h1.Clone("ratio"))
            #den_noerr = copy.deepcopy(h2.Clone("den_noerr"))
            ratio = None
            den_noerr = None
            den = None
            den_noerr_TH1 = None
            if invertRatio:
                ratio = h2_true.Clone("ratio")
                den_noerr = h1_true.Clone("den_noerr")
                den = h1_true.Clone("den")
            else:
                ratio = h1_true.Clone("ratio")
                den_noerr = h2_true.Clone("den_noerr")
                den = h2_true.Clone("den")

            if den_noerr.InheritsFrom("TH1"):
                for iBin in range (1,den_noerr.GetNbinsX()+1):
                    den_noerr.SetBinError(iBin,0.)
                if useDifferenceInLowerPanel:
                    ratio.Add(den_noerr,-1.0)
                    den.Add(den_noerr,-1.0)
                else:
                    ratio.Divide(den_noerr)
                    den.Divide(den_noerr)
            else:
                den_noerr_TH1 = h1_true.Clone("den_noerr_TH1")
                den_noerr_TH1.Reset("ICESM")                        
                for i in range(den_noerr_TH1.GetNbinsX()):
                    yval = den_noerr.Eval( den_noerr_TH1.GetBinCenter(i+1) )
                    den_noerr_TH1.SetBinContent(i+1, yval)
                    den_noerr_TH1.SetBinError(i+1, 0)
                    den_noerr.SetPointEYhigh(i, 0)
                    den_noerr.SetPointEYlow( i, 0)

                if useDifferenceInLowerPanel:
                    ratio.Add(den_noerr_TH1,-1.0)
                    for i in range(den_noerr_TH1.GetNbinsX()):
                        xval = den_noerr_TH1.GetBinCenter(i+1)
                        yval = den_noerr_TH1.GetBinContent(i+1)
                        den.SetPoint(i, xval, 0.0)
                else:
                    ratio.Divide(den_noerr_TH1)
                    for i in range(den_noerr_TH1.GetNbinsX()):
                        xval = den_noerr_TH1.GetBinCenter(i+1)
                        yval = den_noerr_TH1.GetBinContent(i+1)
                        den.SetPoint(i, xval, 1.0)
                        den.SetPointEYhigh(i, den.GetErrorYhigh(i)/yval)
                        den.SetPointEYlow(i, den.GetErrorYlow(i)/yval)

            if invertRatio:
                if histMCpartialUnc_true == None:
                    ratio.SetFillColor(ROOT.kGray+1)
                    ratio.SetFillStyle(fsPart)  # 1001 make it solid again
                    ratio.SetLineColor(ROOT.kGray+3)
                else:                
                    ratio.SetFillColor(colorBandPart["fill"]) # kGreen+1
                    ratio.SetFillColorAlpha(colorBandPart["fill"], fcalphaPart);
                    ratio.SetFillStyle(fsPart)  # 1001 to make it solid again
                    #ratio.SetLineColor(colorBandTot["line"]) # kGreen+2                       
                    ratio.SetLineColor(centralLineColor) # kGreen+2                       
                    ratio.SetLineWidth(linewidth) # make it smaller when it is drawn on top of something
            else:
                if histMCpartialUnc_true == None:
                    den.SetFillColor(ROOT.kGray+1)
                    den.SetFillStyle(1001)  # make it solid again
                    den.SetLineColor(ROOT.kRed)        

            if histMCpartialUnc_true != None:
                h3ratio = h3.Clone("h3ratio")
                h3ratio.SetFillStyle(fsPart) # 
                h3ratio.SetFillColor(colorBandPart["fill"])  # kRed-4
                h3ratio.SetLineColor(colorBandPart["line"])  # kRed-4
                h3ratio.SetFillColorAlpha(colorBandPart["fill"], fcalphaPart);
                h3ratiodummyWhite = h3.Clone("h3ratiodummyWhite")
                h3ratiodummyWhite.SetFillColor(ROOT.kWhite)
                h3ratiodummyWhite.SetLineColor(ROOT.kWhite)
                h3ratiodummyWhite.SetFillStyle(1001)
                if h3ratio.InheritsFrom("TH1"):
                    if useDifferenceInLowerPanel:
                        h3ratio.Add(den_noerr,-1.0)

                    else:
                        h3ratio.Divide(den_noerr)
                    #h3ratio.SetFillStyle(3144) # 1001 for solid, 3144 instead of 3244, to have more dense hatches
                else:
                    if useDifferenceInLowerPanel:
                        for i in range(den_noerr_TH1.GetNbinsX()):
                            xval = den_noerr_TH1.GetBinCenter(i+1)
                            yval = den_noerr_TH1.GetBinContent(i+1)
                            h3ratio.SetPoint(i, xval, 0.0)
                    else:
                        for i in range(den_noerr_TH1.GetNbinsX()):
                            xval = den_noerr_TH1.GetBinCenter(i+1)
                            yval = den_noerr_TH1.GetBinContent(i+1)
                            h3ratio.SetPoint(i, xval, 1.0)
                            h3ratio.SetPointEYhigh(i, h3ratio.GetErrorYhigh(i)/yval)
                            h3ratio.SetPointEYlow(i,  h3ratio.GetErrorYlow(i)/yval)

                h3ratiodummyWhite = h3ratio.Clone("h3ratiodummyWhite")
                h3ratiodummyWhite.SetFillColor(ROOT.kWhite)
                h3ratiodummyWhite.SetLineColor(ROOT.kWhite)
                h3ratiodummyWhite.SetFillStyle(1001)

            frame.Draw()        
            if invertRatio:
                den.SetMarkerSize(0.95)
                den.SetMarkerStyle(20) 
                #den.Draw("EPsame")    # draw after red line, not now            
                ratio.Draw("F2same")
            else:    
                ratio.SetMarkerSize(0.95)
                ratio.SetMarkerStyle(20) 
                den.Draw("F2same")
                #ratio.Draw("EPsame") # draw after red line, not now

            if histMCpartialUnc_true != None:
                h3ratiodummyWhite.Draw("F2 same")
                h3ratio.Draw("F2 same")

            line = ROOT.TF1("horiz_line","0" if useDifferenceInLowerPanel else "1",
                            ratio.GetXaxis().GetBinLowEdge(1),ratio.GetXaxis().GetBinLowEdge(ratio.GetNbinsX()+1))
            line.SetLineWidth(linewidth)
            line.SetLineColor(centralLineColor)
            line.Draw("Lsame")
            # if invertRatio:
            #     ratioline = ratio.Clone("ratioline")
            #     ratioline.SetFillColor(0)
            #     if histMCpartialUnc_true != None and len(drawVertLines):
            #         # unrolled, here the red line would be too large and given the fluctuations is could be seen as an additional band
            #         ratioline.SetLineColor(histMCpartialUnc_true.GetLineColor())
            #     ratioline.SetFillStyle(0)
            #     if histMCpartialUnc_true != None: 
            #         if ratioline.InheritsFrom("TH1"):
            #             ratioline.Draw("HIST same") # to draw the line inside the band for the expected
            #         else:
            #             ratioline.Draw("L same")
            #     den.Draw("EPsame")
            #else: 
            #    ratio.Draw("EPsame")
            if invertRatio:
                den.Draw("E0P0same0") 
            else:
                ratio.Draw("E0P0same0")

            x1leg2 = 0.15 if leftMargin > 0.1 else 0.07
            x2leg2 = 0.55 if leftMargin > 0.1 else 0.27
            y1leg2 = 0.28 if leftMargin > 0.1 else 0.3
            y2leg2 = 0.33 if leftMargin > 0.1 else 0.35
            leg2 = ROOT.TLegend(x1leg2, y1leg2, x2leg2, y2leg2)
            #leg2 = ROOT.TLegend(0.07,0.30,0.27,0.35)
            leg2.SetFillColor(0)
            leg2.SetFillStyle(0)
            leg2.SetBorderSize(0)
            if invertRatio:
                leg2.AddEntry(ratio,"total theory uncertainty","LF")
            else:
                leg2.AddEntry(den,"total theory uncertainty","LF")
            if not noLegendLowerPanel: leg2.Draw("same")
            if histMCpartialUnc_true != None:
                leg2bis = ROOT.TLegend(x1leg2 + 0.45, y1leg2, x2leg2 + 0.45, y2leg2)
                leg2bis.SetFillColor(0)
                leg2bis.SetFillStyle(0)
                leg2bis.SetBorderSize(0)
                leg2bis.AddEntry(h3ratio,histMCpartialUncLegEntry,"LF")
                if not noLegendLowerPanel: leg2bis.Draw("same")

            pad2.RedrawAxis("sameaxis")
            ROOT.gPad.RedrawAxis()
        
        if draw_both0_noLog1_onlyLog2 != 2:
            canvas.SaveAs(outdir + canvasName_true + ".png")
            canvas.SaveAs(outdir + canvasName_true + ".pdf")

        if draw_both0_noLog1_onlyLog2 != 1:        
            if yAxisName == "a.u.": 
                h1_true.GetYaxis().SetRangeUser(max(0.0001,h1_true.GetMinimum()*0.8),h1_true.GetMaximum()*100)
            else:
                h1_true.GetYaxis().SetRangeUser(max(0.001,h1_true.GetMinimum()*0.8),h1_true.GetMaximum()*100)
            canvas.SetLogy()
            canvas.SaveAs(outdir + canvasName_true + "_logY.png")
            canvas.SaveAs(outdir + canvasName_true + "_logY.pdf")
            canvas.SetLogy(0)

        if saveRootFileForHepdata:
            rootfilename = outdir + canvasName_true + ".root"
            rfHep = ROOT.TFile.Open(rootfilename,"recreate")
            if not rfHep:
                print "Error in drawXsecAndTheoryband(): could not open root file %s" % rootfilename
                quit()
            rfHep.cd()
            h1_true.Write("data")
            h2_true.Write("theoryTotalPrediction")
            h3.Write("theoryPDFsAndAlphaS")            
            rfHep.Close()
