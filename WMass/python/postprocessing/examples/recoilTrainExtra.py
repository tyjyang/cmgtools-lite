import ROOT
import os
import numpy as np
ROOT.PyConfig.IgnoreCommandLineOptions = True
from math import *

from CMGTools.WMass.postprocessing.framework.datamodel import Collection,Object
from CMGTools.WMass.postprocessing.framework.eventloop import Module
from PhysicsTools.HeppyCore.utils.deltar import deltaR,deltaPhi
     

class RecoilTrainExtraProducer(Module):
    '''
    Analyzes the event recoil and derives the corrections in pt and phi to bring the estimators to the true vector boson recoil.
    '''
    def __init__(self):
        self.llMassWindow=(91.,15.)

    def beginJob(self):
        """actions taken start of the job"""
        pass

    def endJob(self):
        """actions taken end of the job"""
        pass

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        """actions taken at the start of the file"""

        #start by calling init readers
        self.initReaders(inputTree)

        #define output
        self.out = wrappedOutputTree    
        self.out.branch('isGood','F')
        self.out.branch('isW','F')
        self.out.branch('isZ','F')


    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        """actions taken at the end of the file"""
        pass

    def initReaders(self,tree): # this function gets the pointers to Value and ArrayReaders and sets them in the C++ worker class                                 
        try:
            self.nGenPart = tree.valueReader("nGenPart")
            for var in ["pt","eta","phi","mass","pdgId","status"] : 
                setattr(self,"GenPart_"+var, tree.arrayReader("GenPart_"+var))
        except:
            print '[eventRecoilAnalyzer][Warning] Unable to attach to generator-level particles, only reco info will be made available'
        self._ttreereaderversion = tree._ttreereaderversion

    def getMCTruth(self,event):
        """computes the MC truth for this event"""

        #dummy values
        visibleV,V=VisibleVectorBoson(selLeptons=[]),ROOT.TLorentzVector(0,0,0,0)

        #nothing to do for data :)
        if event.isData : 
            return visibleV,V

        #get the neutrinos (if available, otherwise pass it)
        nuSum=ROOT.TLorentzVector(0,0,0,0)
        try:
            ngenNu = self.out._branches["nGenPromptNu"].buff[0]
            for i in xrange(0,ngenNu):
                p4=ROOT.TLorentzVector(0,0,0,0)
                p4.SetPtEtaPhiM( self.out._branches["GenPromptNu_pt"].buff[i],
                                 self.out._branches["GenPromptNu_eta"].buff[i],
                                 self.out._branches["GenPromptNu_phi"].buff[i],
                                 self.out._branches["GenPromptNu_mass"].buff[i] )
                nuSum+=p4
        except:
            pass

        #construct the visible boson
        try:
            
            dressedLepIds=[]
            dressedLeps=[]            
            for i in xrange(0,self.out._branches["nGenLepDressed"].buff[0]):
                dressedLeps.append( ROOT.TLorentzVector(0,0,0,0) )
                dressedLeps[-1].SetPtEtaPhiM( self.out._branches["GenLepDressed_pt"].buff[i],
                                              self.out._branches["GenLepDressed_eta"].buff[i],
                                              self.out._branches["GenLepDressed_phi"].buff[i],
                                              self.out._branches["GenLepDressed_mass"].buff[i] )
                dressedLepIds.append( self.out._branches["GenLepDressed_pdgId"].buff[i] )

            nLep=len(dressedLeps)
            visibleV=VisibleVectorBoson(selLeptons=[dressedLeps[0]])
            V=visibleV.p4()+nuSum
            if len(dressLeps)>1:
                ll=dressedLeps[0]+dressedLeps[1]
                if abs(dressedLepIds[0])==abs(dressedLepIds[1]):
                    if abs(ll.M()-self.llMassWindow[0])<self.llMassWindow[1] :
                        visibleV=VisibleVectorBoson(selLeptons=[dressedLeps[0],dressedLeps[1]])
                        V=visibleV.p4
        except Exception,e:
            print e
            pass

        return visibleV,V

    def getVisibleV(self,event):
        """
        tries to reconstruct a Z from the selected leptons
        for a W only the leading lepton is returned
        """
        
        lepColl = Collection(event, "LepGood")
        nLep = len(lepColl)

        visibleV=None
        if nLep==0: return visibleV

        visibleV=VisibleVectorBoson(selLeptons=[lepColl[0].p4()])

        #check if the two leading leptons build a Z
        if nLep>1: 
            if abs(lepColl[0].pdgId)==abs(lepColl[1].pdgId):
                ll=lepColl[0].p4()+lepColl[1].p4()
                if abs(ll.M()-self.llMassWindow[0])<self.llMassWindow[1]:
                    visibleV=VisibleVectorBoson(selLeptons=[lepColl[0].p4(),lepColl[1].p4()])

        return visibleV

    def summarizeMetEstimator(self,event, metObjects=['met'], weights=[] ):
        """gets the summary out of a met estimator"""
        p4,sumEt,count=None,0,0
        for i in xrange(0,len(metObjects)):
            wgt=weights[i]
            imet=Object(event,metObjects[i])
            if p4 is None: p4 = imet.p4()*wgt
            else : p4 += imet.p4()*wgt
            sumEt += imet.sumEt*wgt
            try:
                count += getattr(event,'%s_Count'%metObjects[i])*wgt
            except:
                pass
        return p4,sumEt,count

    def getRecoRecoil(self,event,metP4,sumEt,visibleV):
        """computes the MC truth for this event"""

        #charged recoil estimators
        h=ROOT.TVector3(-metP4.Px(),-metP4.Py(),0.)
        ht=sumEt
        if visibleV:
            for l in visibleV.legs:
                ht -= l.Pt()
                h  -= l.Vect()

        return h,max(ht,0.),metP4.M()
    
    def analyze(self, event):
        """process event, return True (go to next module) or False (fail, go to next event)"""

        if event._tree._ttreereaderversion > self._ttreereaderversion:
            self.initReaders(event._tree)


        self.out.fillBranch('isGood',1.)
        self.out.fillBranch('isW',1.)
        self.out.fillBranch('isZ',1.)


        #MC truth
        gen_visibleV,gen_V=self.getMCTruth(event)

        #selected leptons at reco level
        visibleV=self.getVisibleV(event)

        """
        #leading PF candidates
        self.out.fillBranch('leadch_pt',    event.leadCharged_pt)
        self.out.fillBranch('leadch_phi',   event.leadCharged_phi)
        self.out.fillBranch('leadneut_pt',  event.leadNeutral_pt)
        self.out.fillBranch('leadneut_phi', event.leadNeutral_phi)

        #met estimators
        metEstimators={
            'met'               : self.summarizeMetEstimator(event,
                                                             ['met'], 
                                                             [1]),
            'puppimet'          : self.summarizeMetEstimator(event,
                                                             ['puppimet'], 
                                                             [1]),
            'ntmet'             : self.summarizeMetEstimator(event,
                                                             ['ntMet'],       
                                                             [1]),
            'ntcentralmet'      : self.summarizeMetEstimator(event,
                                                             ['ntCentralMet'], 
                                                             [1]),             
            'tkmet'             : self.summarizeMetEstimator(event,
                                                             ['tkMetPVLoose'], 
                                                             [1]),  
            'chsmet'            : self.summarizeMetEstimator(event,
                                                             ['tkMetPVchs'], 
                                                             [1]),  
            'npvmet'            : self.summarizeMetEstimator(event,
                                                             ['tkMetPUPVLoose'],                               
                                                             [1]),  
            'ntnpv'             : self.summarizeMetEstimator(event,
                                                             ['ntMet','tkMetPUPVLoose'],                       
                                                             [1,1]),  
            'centralntnpv'      : self.summarizeMetEstimator(event,
                                                             ['ntCentralMet','tkMetPUPVLoose'],                
                                                             [1,1]),  
            'centralmetdbeta'   : self.summarizeMetEstimator(event,
                                                             ['tkMetPVLoose','ntCentralMet','tkMetPUPVLoose'], 
                                                             [1,1,0.5]),  
            }
       
        metEstimatorsList=metEstimators.keys()
        if not event.isData: metEstimatorsList+=["truth","gen"]
        
        #recoil estimators
        for metType in metEstimatorsList:

            if not metType in ["truth","gen","puppimet", 'ntmet','ntcentralmet', 'tkmet', 'npvmet', 'ntnpv']: continue

            #some may need to be specified
            if metType=="truth":
                vis   = gen_visibleV.VectorT()
                metP4 = gen_V
                h     = ROOT.TVector3(-gen_V.Px(),-gen_V.Py(),0)                
                ht    = gen_V.Pt()
                m     = gen_V.M()
                count = 0
            elif metType=='gen':
                vis   = gen_visibleV.VectorT()
                metP4=ROOT.TLorentzVector(0,0,0,0)
                metP4.SetPtEtaPhiM(event.tkGenMet_pt,0,event.tkGenMet_phi,0.)
                h     = gen_h
                ht    = gen_ht
                m     = 0
                count = 0
            else:
                vis = visibleV.VectorT()
                metP4, sumEt, count = metEstimators[metType]
                h, ht, m = self.getRecoRecoil(event=event,metP4=metP4,sumEt=sumEt,visibleV=visibleV)

            #save information to tree
            pt=h.Pt()
            phi=h.Phi()
            metphi=metP4.Phi()
            sphericity=pt/ht if ht>0 else -1               
            #e1=gen_V.Pt()/h.Pt() if h.Pt()>0 else -1
            #e2=deltaPhi(gen_V.Phi()+np.pi,h.Phi())
            #mt2= 2*vis.Pt()*((vis+h).Pt())+vis.Pt()**2+vis.Dot(h) 
            #mt=np.sqrt(mt2) if mt2>=0. else -np.sqrt(-mt2)

            self.out.fillBranch('%s_recoil_pt'%metType,            pt)
            self.out.fillBranch('%s_recoil_phi'%metType,           phi)
            self.out.fillBranch('%s_m'%metType,                    m)
            self.out.fillBranch('%s_recoil_sphericity'%metType,    sphericity)
            self.out.fillBranch('%s_n'%metType,                    count)
            #self.out.fillBranch('%s_recoil_e1'%metType,            e1)
            #self.out.fillBranch('%s_recoil_e2'%metType,            e2)
            #self.out.fillBranch('%s_mt'%metType,                   mt)
            self.out.fillBranch('%s_dphi2met'%metType,             deltaPhi(metEstimators['met'][0].Phi(),             metphi) )
            self.out.fillBranch('%s_dphi2puppimet'%metType,        deltaPhi(metEstimators['puppimet'][0].Phi(),        metphi) )
            self.out.fillBranch('%s_dphi2ntnpv'%metType,           deltaPhi(metEstimators['ntnpv'][0].Phi(),           metphi) )
            self.out.fillBranch('%s_dphi2centralntnpv'%metType,    deltaPhi(metEstimators['centralntnpv'][0].Phi(),    metphi) )
            self.out.fillBranch('%s_dphi2centralmetdbeta'%metType, deltaPhi(metEstimators['centralmetdbeta'][0].Phi(), metphi) )
            self.out.fillBranch('%s_dphi2leadch'%metType,          deltaPhi(event.leadCharged_phi,                     metphi) )
            self.out.fillBranch('%s_dphi2leadneut'%metType,        deltaPhi(event.leadNeutral_phi,                     metphi) )

            """

        return True


# define modules using the syntax 'name = lambda : constructor' to avoid having them loaded when not needed
recoilTrainExtra = lambda : RecoilTrainExtraProducer()
