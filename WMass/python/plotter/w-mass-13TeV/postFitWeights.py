#!/usr/bin/env python
import ROOT,re,math
import numpy as np
from scipy.linalg import cholesky
from scipy.stats import norm
from array import array
from make_helicity_cards import wptBinsScales

def qcdScalesNames():
    return [pol+qcdscale+str(ptbin)+charge for pol in ['left','right','long'] for qcdscale in ['muR','muF','muRmuF'] for ptbin in xrange(1,11) for charge in ['plus','minus']]

def pdfNames():
    return ['pdf{i}'.format(i=i) for i in xrange(1,61)]+['alphaS']

def fitValues(frfile,params):
    tf = ROOT.TFile(frfile)
    fitres = tf.Get("fitresults")
    scales = []
    fitres.GetEntry(0)
    for name in qcdScalesNames():
        scales.append( (name,getattr(fitres,name),getattr(fitres,name+'_err')) )
    tf.Close()
    return scales

def getSubMatrix(frfile,pars):
    tf = ROOT.TFile(frfile)
    covmatrix = tf.Get('covariance_matrix_channelmu')
    params = []; indices = []
    for p in pars:
        for ib in range(1+covmatrix.GetNbinsX()+1):
            if re.match("^"+p[0]+"$", covmatrix.GetXaxis().GetBinLabel(ib)):
                ## store mean and rms into the dictionaries from before
                ## also keep a list of the parameter names, for sorting
                params .append(covmatrix.GetXaxis().GetBinLabel(ib))
                indices.append(ib)

    cov = np.zeros(shape=(len(params),len(params)))
    ## construct the covariances and the correlations in one go.
    for ip1, p1 in enumerate(params):
        for ip2, p2 in enumerate(params):
            cov [ip1,ip2] = covmatrix .GetBinContent(indices[ip1],indices[ip2])
    tf.Close()
    return cov

def getCorrelatedParameters(params,covmat):
    c = cholesky(covmat,lower=True)
    # generate samples from N independent normally distributed random variables
    # with mean the measured bias and std. dev. 1
    num_samples = int(1e+6)
    par_samples = []
    print "generating ",num_samples," Toy MC samples to build the correlated samples..."
    for par in params:
        x = norm.rvs(loc=par[1],size=(num_samples))
        par_samples.append(x)
    biases = np.stack(par_samples)
    # Convert the fit biases to correlated random variables with the covmat
    y = np.dot(c,biases)
    # this should be ~the same as the original covariance matrix 
    #print "cov corr = ",np.cov(y)
    corr_pars = {}
    for ip,par in enumerate(params):
        corr_pars[par[0]] = (np.mean(y[ip]),np.std(y[ip]))
    return corr_pars

def getUnpolSyst(syst):
    unpolsyst = syst.replace('left','').replace('right','').replace('long','')
    unpolsyst = unpolsyst.replace('plus','').replace('minus','')
    return unpolsyst
    
def loadKappas(xsecfile):
    tf = ROOT.TFile(xsecfile)
    print "Now loading the arrays of kappas from ",xsecfile,". May take time..."
    kappas = {}
    edges = []
    for charge in ['plus','minus']:
        for pol in ['left','right','long']:
            for syst in qcdScalesNames():
                unpolsyst = getUnpolSyst(syst)
                # don't waste time to load copies of the same syst. Just rename them later
                if (charge,pol,unpolsyst) in kappas: continue
                hcen = 'w{charge}_wpt_central_W{charge}_{pol}'.format(charge=charge,pol=pol) 
                hup  = 'w{charge}_wpt_{syst}Up_W{charge}_{pol}'.format(charge=charge,syst=unpolsyst,pol=pol)
                hdn  = hup.replace('Up','Dn')
                thcen = tf.Get(hcen)
                thup  = tf.Get(hup)
                thdn  = tf.Get(hdn)
                nbins = thup.GetNbinsX()
                kappasup = np.zeros(shape=nbins)
                kappasdn = np.zeros(shape=nbins)
                for b in xrange(0,nbins):
                    cen,up,dn = thcen.GetBinContent(b+1), thup.GetBinContent(b+1), thdn.GetBinContent(b+1)
                    kappasup[b] = math.log(max(up/cen,dn/cen)) if min(min(up,dn),cen)>0 else 0.
                    kappasdn[b] = math.log(min(up/cen,dn/cen)) if min(min(up,dn),cen)>0 else 0.
                    if len(edges)<nbins: edges.append(thdn.GetBinLowEdge(b+1))
                kappas[(charge,pol,unpolsyst)] = (kappasup,kappasdn)
    edges.append(thdn.GetBinLowEdge(nbins+1)) # they are all the same. Take the last bin upper edge from the last used one
    tf.Close()
    return edges,kappas

def kappa(kappas,charge,pol,syst,kbin):
    unpolsyst = getUnpolSyst(syst)
    kappa_up = kappas[(charge,pol,unpolsyst)][0]
    kappa_dn = kappas[(charge,pol,unpolsyst)][1]
    return (kappa_up[kbin],kappa_dn[kbin])

def getQCDWeightProd(kappas,thetas,pol,charge,nuisbin,kbin):
    thetak = {} # not really necessary, useful for debugging
    for qcdscale in ['muR','muF','muRmuF']:
        nuis = pol+qcdscale+str(nuisbin)+charge
        theta = thetas[nuis][0]
        ## should inteprolate, but let's just separate Up/Dn
        side = 0 if theta>0 else 1
        k = kappa(kappas,charge,pol,nuis,kbin)[side]
        thetak[nuis] = (abs(theta),k)
    #print "Calc wgt for ",pol,"  ",charge
    #print "thetak = ",thetak
    wgt = np.prod([math.exp(nuis[0]*nuis[1]) for k,nuis in thetak.iteritems()])
    return wgt

if __name__ == "__main__":

    from optparse import OptionParser
    parser = OptionParser(usage="%prog [options] fitresults.root")
    (options, args) = parser.parse_args()
    fitresfile = args[0]

    polarizations = ['left','right','long']
    charges = ['plus','minus']

    ## load the kappas
    ## theta=1  corresponds to doubling the scale
    ## theta=-1 corresponds to halfing it
    ## Translating this into a cross section requires knowing by how much the xsec changed from the X2 or x1/2 variations in the first place, Ie the k_i
    xsecfile = "/afs/cern.ch/user/m/mdunser/public/wpTspectraFixed.root"
    kedges,kappas = loadKappas(xsecfile)

    postfit_weights_file = ROOT.TFile('postfit_wgts.root','recreate')

    ## QCD scales
    scales = fitValues(fitresfile,qcdScalesNames())
    covsubmat = getSubMatrix(fitresfile,scales)
    thetas = getCorrelatedParameters(scales,covsubmat)
    
    # nuisances have a coarse pt binning
    nuis_edges = [wptBinsScales(i)[0] for i in xrange(1,11)] + [13000.0]
    hnuisscales = ROOT.TH1F('hnuisscales','',len(nuis_edges)-1,array('f',nuis_edges))

    # kappas have a more granular binning (e.g. 1 GeV in wpt)
    hscale = ROOT.TH1F('scales_postfit_wgts','',len(kedges)-1,array('f',kedges))
    hwgts = []
    for pol in polarizations:
        for charge in charges:
            hwgt = hscale.Clone('weights_{pol}{charge}'.format(pol=pol,charge=charge))
            hwgt.SetDirectory(None)
            for ipt in xrange(len(kedges)-1):
                iptcenter = hwgt.GetBinCenter(ipt+1)
                nuisbin = hnuisscales.FindBin(iptcenter) 
                wgt = getQCDWeightProd(kappas,thetas,pol,charge,nuisbin,ipt)
                hwgt.SetBinContent(ipt,wgt)
            hwgts.append(hwgt)
    postfit_weights_file.cd()
    for h in hwgts:
        h.Write()

    postfit_weights_file.Close()
    print "DONE. Weights are stored in ",postfit_weights_file

                        
    
