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
    num_samples = 100000
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
    corr_pars = []
    for ip,par in enumerate(params):
        corr_pars.append( (par[0],np.mean(y[ip]),np.std(y[ip])) )
    return corr_pars

def getWeightProd(allthetas,pattern):
    seltheta = filter(lambda x: re.match(pattern,x[0]),allthetas)
    ## weight is e^(k1theta1)*e^(k2theta2). The thetas here are already normalized, so k is already there
    return np.prod([math.exp(x[1]) for x in seltheta])

if __name__ == "__main__":

    from optparse import OptionParser
    parser = OptionParser(usage="%prog [options] fitresults.root")
    (options, args) = parser.parse_args()
    fitresfile = args[0]

    postfit_weights_file = ROOT.TFile('postfit_wgts.root','recreate')

    ## QCD scales
    scales = fitValues(fitresfile,qcdScalesNames())
    covsubmat = getSubMatrix(fitresfile,scales)
    theta = getCorrelatedParameters(scales,covsubmat)
    edges   = [float(wptBinsScales(ipt)[0]) for ipt in xrange(1,11)] + [float(13000.0)]
    hscales = ROOT.TH1F('scales_postfit_wgts','',len(edges)-1,array('f',edges))
    for ipt in xrange(1,11):
        hscales.SetBinContent( ipt, getWeightProd(theta,'.*(muR|muF){ipt}(plus|minus)'.format(ipt=ipt)) )
    postfit_weights_file.cd()
    hscales.Write()

    ## PDFs
    pdfs = fitValues(fitresfile,pdfNames())
    covsubmat = getSubMatrix(fitresfile,pdfs)
    theta = getCorrelatedParameters(pdfs,covsubmat)
    edges   = [float(0.5+i) for i in xrange(61)]
    hpdfs = ROOT.TH1F('pdfs_postfit_wgts','',len(edges)-1,array('f',edges))
    for ipdf in xrange(1,61):
        hpdfs.SetBinContent( ipdf, getWeightProd(theta,'pdf{ipdf}'.format(ipdf=ipdf)) )
    postfit_weights_file.cd()
    hpdfs.Write()

    postfit_weights_file.Close()
    print "DONE. Weights are stored in ",postfit_weights_file

                        
    
