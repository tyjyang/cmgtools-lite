#!/usr/bin/env python                                                       
 
#from shutil import copyfile
import re, sys, os, os.path, ROOT, copy, math
from array import array
#import numpy as np

from w_helicity_13TeV.make_diff_xsec_cards import getXYBinsFromGlobalBin
from w_helicity_13TeV.make_diff_xsec_cards import getArrayParsingString
from w_helicity_13TeV.make_diff_xsec_cards import getDiffXsecBinning
from w_helicity_13TeV.make_diff_xsec_cards import templateBinning

from w_helicity_13TeV.rollingFunctions import roll1Dto2D, dressed2D, unroll2Dto1D

from w_helicity_13TeV.mergeCardComponentsAbsY import mirrorShape
from w_helicity_13TeV.mergeCardComponentsAbsY import putUncorrelatedFakes
#from w_helicity_13TeV.mergeCardComponentsAbsY import putEffStatHistos   # use function below that can manage case with wider template bins along eta

def putEffStatHistosDiffXsec(infile,regexp,charge, outdir=None, isMu=True):

    # for differential cross section I don't use the same option for inputs, so I pass it from outside
    indir = outdir if outdir != None else options.inputdir

    # this doesn't work for differential cross section because I have line 2 with the comment for gen binning
    # in all my scripts I developed some specific functions to read the binning, either reco or gen: here we only need reco
    # binninPtEtaFile = open(indir+'/binningPtEta.txt','r')
    # bins = binninPtEtaFile.readlines()[1].split()[1]
    # etabins = list( float(i) for i in bins.replace(' ','').split('*')[0].replace('[','').replace(']','').split(',') )
    # ptbins  = list( float(i) for i in bins.replace(' ','').split('*')[1].replace('[','').replace(']','').split(',') )
    # nbinseta = len(etabins)-1
    # nbinspt  = len( ptbins)-1
    # binning = [nbinseta, etabins, nbinspt, ptbins]

    # get eta-pt binning for reco 
    etaPtBinningVec = getDiffXsecBinning(indir+'/binningPtEta.txt', "reco")  # this get two vectors with eta and pt binning
    recoBins = templateBinning(etaPtBinningVec[0],etaPtBinningVec[1])        # this create a class to manage the binnings
    binning = [recoBins.Neta, recoBins.etaBins, recoBins.Npt, recoBins.ptBins]
    etabins = recoBins.etaBins

    # not needed
    #gen_etaPtBinningVec = getDiffXsecBinning(indir+'/binningPtEta.txt', "gen")  # this get two vectors with eta and pt binning
    #genBins = templateBinning(gen_etaPtBinningVec[0],gen_etaPtBinningVec[1])        # this create a class to manage the binnings


    basedir = '/afs/cern.ch/work/m/mdunser/public/cmssw/w-helicity-13TeV/CMSSW_8_0_25/src/CMGTools/WMass/python/postprocessing/data/'    
    if isMu:
        parfile_name = basedir+'/leptonSF/new2016_madeSummer2018/systEff_trgmu_{ch}_mu.root'.format(ch=charge)
    else:
        parfile_name = basedir+'/leptonSF/new2016_madeSummer2018/systEff_trgel.root'

    parfile = ROOT.TFile(parfile_name, 'read')
    

    tmp_infile = ROOT.TFile(infile, 'read')

    flavour = "mu" if isMu else "el"
    outfile = ROOT.TFile(indir+'/ErfParEffStat_{fl}_{ch}.root'.format(fl=flavour, ch=charge), 'recreate')

    ndone = 0
    for k in tmp_infile.GetListOfKeys():
        tmp_name = k.GetName()
        ## don't reweight any histos that don't match the regexp
        if not re.match(regexp, tmp_name): continue
        ## don't reweight any histos that are already variations of something else
        if 'Up' in tmp_name or 'Down' in tmp_name: continue

        #if ndone: continue
        ndone += 1

        ## now should be left with only the ones we are interested in
        print 'reweighting erfpareffstat nuisances for process', tmp_name
        
        tmp_nominal = tmp_infile.Get(tmp_name)
        tmp_nominal_2d = dressed2D(tmp_nominal,binning, tmp_name+'backrolled')

        # the erfPar histogram has 50 bins for muons and electrons (for muons, the extremes are equal to the neighbours because eta goes up to 2.4 at most)
        # the binning in eta is always 0.1 wide
        # if the template has coarser granularity, that bin will have N variations, for the N 0.1-wide bins contained in it
        # this would require using the loop to reweight the histograms, because each on the N-th variation would only affect 1/N of the events in the wider bin
        # so, here we will reduce the scaling by 1/N as well.

        # get width of bins in multiples of 0.1
        # for electrons the region around the gap is odd
        # we don't need this value to be integer, otherwise keep in mind that 0.1 is not perfectly represented with float in python
        # which entails that 0.2/0.1 might not yield 2.0, but 1.999
        binwidths = []
        for ieta in range(1,tmp_nominal_2d.GetNbinsX()+1):
            binwidths.append(tmp_nominal_2d.GetXaxis().GetBinWidth(ieta)/0.1)  

        nEtaErfPar = 48 if isMu else 50

        ## loop over the three parameters
        for npar in range(3):
            parhist = parfile.Get('p'+str(npar))
            ## loop over all eta bins of the 2d histogram
            for ietaErf in range(1,nEtaErfPar+1):

                # FIXME, need more care for electron channel, template binning can be odd around the gap
                # identify template eta which contain that erfPar eta

                etabinOffset = 0
                # the parhist histogram for muons is defined with 50 bins, where the two with |eta| in [2.4,2.5] are filled like the inner ones
                # these bins for muons do not make sense: when looping on the bin number, we need to add an offset to skip the first bin
                if isMu and parhist.GetNbinsX() == 50:
                    etabinOffset = 1

                ietaTemplate = tmp_nominal_2d.GetXaxis().FindFixBin( parhist.GetXaxis().GetBinCenter(ietaErf+etabinOffset) )
                #ietaTemplate = tmp_nominal_2d.GetXaxis().FindFixBin( parhist.GetXaxis().GetBinCenter(ietaErf + parhistBinOffset) )

                # if template has eta range narrower than parhist, continue
                if ietaTemplate == 0 or ietaTemplate == (1 + tmp_nominal_2d.GetNbinsX()):
                    continue

                if not isMu:
                    if parhist.GetXaxis().GetBinCenter(ietaErf) > 1.4 and parhist.GetXaxis().GetBinCenter(ietaErf) < 1.5:
                        ietaTemplate = tmp_nominal_2d.GetXaxis().FindFixBin( 1.41 )

                # for electrons there is the gap:
                # eg, we can have 1.3, 1.4442, 1.5, 1.566, 1.7
                # and the bin centers of the erfPar are 1.35, 1.45, 1.55, ...
                # so with 1.35 we get the template bin from 1.3 to 1.4442
                # but with 1.45 we would fall in the empty bin, while we would want to reweight the template bin that would span between 1.4 to 1.5
                # in this case, the syst is ineffective. The same for 1.55, which select the other "side" of the empty bin around 1.5
                # one could device a fancy way to isolate the gap and assign a larger uncertainty around it but let's keep it simple

                #print "ErfPar{p}EffStat{ieta}".format(p=npar,ieta=ietaErf)
                outname_2d = tmp_nominal_2d.GetName().replace('backrolled','')+'_ErfPar{p}EffStat{ieta}{fl}{ch}2DROLLED'.format(p=npar,ieta=ietaErf,fl=flavour,ch=charge)
            
                tmp_scaledHisto_up = copy.deepcopy(tmp_nominal_2d.Clone(outname_2d+'Up'))
                tmp_scaledHisto_dn = copy.deepcopy(tmp_nominal_2d.Clone(outname_2d+'Down'))
                
                ## loop over all pT bins in that bin of eta (which is ieta)
                for ipt in range(1,tmp_scaledHisto_up.GetNbinsY()+1):
                    tmp_bincontent = tmp_scaledHisto_up.GetBinContent(ietaTemplate, ipt)
                    ybincenter = tmp_scaledHisto_up.GetYaxis().GetBinCenter(ipt)
                    ## now get the content of the parameter variation!
                    phistybin = parhist.GetYaxis().FindBin(ybincenter)
                    tmp_scale = parhist.GetBinContent(ietaErf, phistybin)
                    scaling = math.sqrt(2.)*tmp_scale
                    ## scale electrons by sqrt(2) due to the input file being charge inclusive
                    if flavour == 'el':
                        scaling *= math.sqrt(2)
                    # modify scaling if template bin has larger width than the ErfPar histogram
                    if binwidths[ietaTemplate-1] > 1.: scaling = scaling / binwidths[ietaTemplate-1]
                    ## scale up and down with what we got from the histo

                    tmp_bincontent_up = tmp_bincontent*(1.+scaling)
                    tmp_bincontent_dn = tmp_bincontent*(1.-scaling)
                    tmp_scaledHisto_up.SetBinContent(ietaTemplate, ipt, tmp_bincontent_up)
                    tmp_scaledHisto_dn.SetBinContent(ietaTemplate, ipt, tmp_bincontent_dn)

                ## re-roll the 2D to a 1D histo
                tmp_scaledHisto_up_1d = unroll2Dto1D(tmp_scaledHisto_up, newname=tmp_scaledHisto_up.GetName().replace('2DROLLED',''))
                tmp_scaledHisto_dn_1d = unroll2Dto1D(tmp_scaledHisto_dn, newname=tmp_scaledHisto_dn.GetName().replace('2DROLLED',''))

                outfile.cd()
                tmp_scaledHisto_up_1d.Write()
                tmp_scaledHisto_dn_1d.Write()
    outfile.Close()
    print 'done with the many reweightings for the erfpar effstat'


# for Z use this function to make EffStat variations instead of doing it with jobs. Signal is made with the loop, but could be done like this as well
# note that this function is based on an approximation, but since we apply the variation in an eta-pt dependent way, it is a very good approximation
# by doing it for signal as well, the loop and unpacking of TH3 to TH1 could be made much faster, maybe

# originally developed to make cards for diff xsec in lepton pt and |eta|
# merge signal shapes with Z and data+backgrounds
# first, merge shapes with Z and remove data from there (data is just the MC)
# then, merge signal, data and merged Z

# Z and data are assumed to be in folder given with --indir-bkg
# signal is assumed to be named W<flavour>_<charge>_shapes_signal.root

## python mergeRootComponentsDiffXsec.py -f mu -c minus --indir-bkg  cards/diffXsec_mu_2018_11_24_group10_onlyBkg/part0/  --indir-sig cards/diffXsec_mu_2018_11_24_group10_onlyBkg/ -o cards/diffXsec_mu_2018_11_24_group10_onlyBkg/ -d

from optparse import OptionParser
parser = OptionParser(usage="%prog [options]")
parser.add_option("-d", "--dry-run", dest="dryrun",   action="store_true", default=False, help="Dry run: print commands but do not merge");
parser.add_option(      "--indir-bkg", dest="indirBkg", type="string", default="", help="Input folder with root files for Z, data and other backgrounds");
parser.add_option(      "--indir-sig", dest="indirSig", type="string", default="", help="Input folder with root files for signal");
parser.add_option("-o", "--outdir",    dest="outdir", type="string", default="./", help="Output folder (current one as default)");
parser.add_option("-n", "--name",      dest="name",   type="string", default="", help="Name for output root file (if not given, name is W<flavou>_<charge>_shapes.root ).");
parser.add_option("-s", "--suffix",    dest="suffix", type="string", default="", help="Suffix to add to output file before extension. Ineffective if using option -n < name>");
parser.add_option("-f", "--flavour",   dest="flavour", type="string", default='', help="Channel: either 'el' or 'mu'");
parser.add_option("-c", "--charge",    dest="charge", type="string", default='', help="Charge: either 'plus' or 'minus'");
parser.add_option(      "--etaBordersForFakesUncorr",    dest="etaBordersForFakesUncorr", type="string", default='0.5,1.0,1.5,2.0', help="Borders passed to function that creates the eta-uncorrelated normalization for fakes. Pass comma separated list (no 0 or outer edges, and only positive values");
parser.add_option(      "--no-qcdsyst-Z", dest="useQCDsystForZ", action="store_false", default=True, help="If False, do not store the muR,muF,muRmuF variations for Z (if they were present)");
#parser.add_option(      "--no-effstatsyst-Z", dest="useEffstatsystForZ", action="store_false", default=True, help="If False, do not create the Effstat variations for Z");
parser.add_option(       '--uncorrelate-fakes-by-charge', dest='uncorrelateFakesByCharge' , default=False, action='store_true', help='If True, nuisances for fakes are uncorrelated between charges (Eta, PtSlope, PtNorm)')
(options, args) = parser.parse_args()
    
# manage output folder
outdir = options.outdir
if not outdir.endswith('/'): outdir += "/"
if outdir != "./":
    if not os.path.exists(outdir):
        print "Creating folder", outdir
        os.system("mkdir -p " + outdir)

if options.flavour not in ["el", "mu"]:
    print "Warning: you must specify a lepton flavour with option -f el|mu"
    quit()
if options.charge not in ["plus", "minus"]:
    print "Warning: you must specify a charge with option -c plus|minus"
    quit()

charge = options.charge
flavour = options.flavour
#binname = options.bin if len(options.bin) else "W%s" % flavour

if not len(options.indirBkg):
    print "Warning: you must specify a folder with data and background root files with option --indir-bkg"
    quit()
if not len(options.indirSig):
    print "Warning: you must specify a folder with signal root files with option --indir-sig"
    quit()

# define ultimate output file
shapename = ""
if options.name == "":
    shapename = "{od}W{fl}_{ch}_shapes.root".format(od=outdir, fl=flavour, ch=charge)
    if len(options.suffix): shapename = shapename.replace(".root","_{sf}.root".format(sf=options.suffix))
else:
    shapename = outdir + options.name

print ""
print "-"*20
print ""

## prepare the relevant files. First merge Z with correct charge
zMatch = "^Z_{fl}_{ch}.*".format(fl=flavour,ch=charge)
zfiles = [os.path.join(dp, f) for dp, dn, fn in os.walk(options.indirBkg) for f in fn if (f.endswith('.input.root') and re.match(zMatch,f))]

#for f in zfiles:
#    print f

# merge Z files (will have to remove x_data_obs later)
tmpZfile="{od}Z_{fl}_{ch}_mergeTMP.root".format(od=outdir, fl=flavour, ch=charge)
cmdMerge = "hadd -f -k -O {tmp} {zfs}".format(tmp=tmpZfile, zfs=" ".join([f for f in zfiles]))
#print cmdMerge
Zfile="{od}Z_{fl}_{ch}_merge.root".format(od=outdir, fl=flavour, ch=charge)

print "Merging Z"
if not options.dryrun: os.system(cmdMerge)

print "Remove x_data_obs from Z, and replace 'x_Z_dy_' with 'x_Z_'"
print "Also changing Dn to Down"
if not options.useQCDsystForZ:
    print "Will reject the QCD scales variations on Z (muR, muF, muRmuF)"
nZcopied = 0
# open new Z file to remove data from input file
#----------------------------------
newname = ""
zpdf = []
znominal = None
nMirroredPDF = 0
if not options.dryrun:
    tf = ROOT.TFile.Open(tmpZfile,"READ")
    if not tf or not tf.IsOpen():
        raise RuntimeError('Unable to open file {fn}'.format(fn=tmpZfile))

    # open output file
    of = ROOT.TFile(Zfile,'recreate')
    if not of or not of.IsOpen():
        raise RuntimeError('Unable to open file {fn}'.format(fn=Zfile))

    for ikey,e in enumerate(tf.GetListOfKeys()):
        name = e.GetName()
        obj  = e.ReadObj()
        if not obj:
            raise RuntimeError('Unable to read object {n}'.format(n=name))
        if "data_obs" in name: continue
        if not options.useQCDsystForZ:
            if any(x in name for x in ["muR", "muF", "muRmuF"]): continue
        if "x_Z_dy_" in name: newname = name.replace("x_Z_dy_","x_Z_")
        else: newname = name
        if name == "x_Z": znominal = obj.Clone(name)
        if "_pdf" in name:
            zpdf.append(obj.Clone(newname))
            continue
        if newname.endswith("Dn"): newname = newname.replace("Dn","Down")
        newobj = obj.Clone(newname)
        newobj.Write(newname)
        nZcopied += 1

    for h in zpdf:
        nMirroredPDF += 1
        (alternate,mirror) = mirrorShape(znominal,h,h.GetName(),use2xNomiIfAltIsZero=True)
        for alt in [alternate,mirror]:
            alt.Write()

    of.Close()
    tf.Close()
#----------------------------------

print "Copied {n} histograms in {zf} (x_data_obs removed)".format(n=str(nZcopied),zf=Zfile)
print "Created {n} mirrored histograms for PDFs in {zf} ".format(n=str(nMirroredPDF),zf=Zfile)
print "Removing temporary file {tmp}".format(tmp=tmpZfile)
if not options.dryrun: os.system("rm {tmp}".format(tmp=tmpZfile))    
print "-"*20
print ""

fileZeffStat = "{od}ErfParEffStat_{fl}_{ch}.root".format(od=outdir, fl=flavour, ch=options.charge)  
# this name is used inside putEffStatHistosDiffXsec (do not change it outside here)
print "Now adding ErfParXEffStatYY systematics to x_Z process"
print "Will create file --> {of}".format(of=fileZeffStat)
if not options.dryrun: putEffStatHistosDiffXsec(Zfile, 'x_Z', charge, outdir, isMu=True if flavour=="mu" else False)
print ""

dataAndBkgFile = "{obkg}bkg_and_data_{fl}_{ch}.input.root".format(obkg=options.indirBkg, fl=flavour, ch=charge)
dataAndBkgFileTmp = dataAndBkgFile.replace(".input.root","TMP.input.root")
print "Creating temporary file {bkg} to remove 'x_data' histogram".format(bkg=dataAndBkgFileTmp)
print "Also changing Dn to Down"
print "-"*30
print "WARNING: will also reject histograms whose name starts with x_Z, which are not supposed to stay in this file"
print "-"*30
# now remove x_data from bkg_and_data_*.input.root (we only need x_data_obs, this x_data should be removed)
#----------------------------------
if not options.dryrun:
    tf = ROOT.TFile.Open(dataAndBkgFile,"READ")
    if not tf or not tf.IsOpen():
        raise RuntimeError('Unable to open file {fn}'.format(fn=dataAndBkgFile))

    # open output file
    of = ROOT.TFile(dataAndBkgFileTmp,'recreate')
    if not of or not of.IsOpen():
        raise RuntimeError('Unable to open file {fn}'.format(fn=dataAndBkgFileTmp))
    
    nKeys = tf.GetNkeys()
    nCopiedKeys = 0
    for ikey,e in enumerate(tf.GetListOfKeys()):
        name = e.GetName()
        obj  = e.ReadObj()
        if not obj:
            raise RuntimeError('Unable to read object {n}'.format(n=name))
        if name == "x_data": continue
        if "x_Z_" in name: continue  # patch, before February 2019 the muscale and lepeff systs on Z appear in the background file by mistake
        newname = name
        if newname.endswith("Dn"): 
            newname = newname.replace("Dn","Down")        
        newobj = obj.Clone(newname)
        newobj.Write(newname)        
        nCopiedKeys += 1

    print "Copied {n}/{tot} from {bkg}".format(n=str(nCopiedKeys),tot=str(nKeys),bkg=dataAndBkgFile)
    of.Close()
    tf.Close()
#----------------------------------
print "-"*20
print ""

fileFakesEtaUncorr = "{od}FakesEtaUncorrelated_{fl}_{ch}.root".format(od=outdir, fl=flavour, ch=options.charge) 
fileFakesPtSlopeUncorr  = "{od}FakesPtSlopeUncorrelated_{fl}_{ch}.root".format( od=outdir, fl=flavour, ch=options.charge) 
fileFakesPtNormUncorr  = "{od}FakesPtNormUncorrelated_{fl}_{ch}.root".format( od=outdir, fl=flavour, ch=options.charge) 
# these names are used inside putUncorrelatedFakes (do not change them outside here)
print "Now adding FakesEtaUncorrelated and FakesPtSlopeUncorrelated and FakesPtNormUncorrelated systematics to x_data_fakes process"
print "Will create file --> {of}".format(of=fileFakesEtaUncorr)
print "Will create file --> {of}".format(of=fileFakesPtSlopeUncorr)
print "Will create file --> {of}".format(of=fileFakesPtNormUncorr)
etaBordersForFakes = [float(x) for x in options.etaBordersForFakesUncorr.split(',')]
if not options.dryrun: 
    putUncorrelatedFakes(dataAndBkgFileTmp, 'x_data_fakes', charge, outdir, isMu=True if flavour=="mu" else False, etaBordersTmp=etaBordersForFakes, 
                         doType='eta', uncorrelateCharges=options.uncorrelateFakesByCharge)
    putUncorrelatedFakes(dataAndBkgFileTmp, 'x_data_fakes', charge, outdir, isMu=True if flavour=="mu" else False, etaBordersTmp=etaBordersForFakes, 
                         doType='ptslope', uncorrelateCharges=options.uncorrelateFakesByCharge)
    putUncorrelatedFakes(dataAndBkgFileTmp, 'x_data_fakes', charge, outdir, isMu=True if flavour=="mu" else False, etaBordersTmp=etaBordersForFakes, 
                         doType='ptnorm', uncorrelateCharges=options.uncorrelateFakesByCharge)

print "Now merging signal + Z + data + other backgrounds + Fakes*Uncorrelated + ZEffStat"
sigfile = "{osig}W{fl}_{ch}_shapes_signal.root".format(osig=options.indirSig, fl=flavour, ch=charge)

cmdFinalMerge="hadd -f -k -O {of} {sig} {zf} {bkg} {fakesEta} {fakesPtSlope} {fakesPtNorm} {zEffStat}".format(of=shapename, 
                                                                                                              sig=sigfile, 
                                                                                                              zf=Zfile, 
                                                                                                              bkg=dataAndBkgFileTmp, 
                                                                                                              fakesEta=fileFakesEtaUncorr,
                                                                                                              fakesPtSlope=fileFakesPtSlopeUncorr,
                                                                                                              fakesPtNorm=fileFakesPtNormUncorr,
                                                                                                              zEffStat=fileZeffStat)
print "Final merging ..."
print cmdFinalMerge
if not options.dryrun: os.system(cmdFinalMerge)

print "-"*20
print ""
print "Removing temporary file {bkg}".format(bkg=dataAndBkgFileTmp)
if not options.dryrun: os.system("rm {bkg}".format(bkg=dataAndBkgFileTmp))

print "-"*20
print ""
print ""
print "Wrote root file in %s" % shapename
print ""
