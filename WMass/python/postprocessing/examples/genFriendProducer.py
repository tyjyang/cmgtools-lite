import ROOT
import os, array, copy
import numpy as np
ROOT.PyConfig.IgnoreCommandLineOptions = True

from CMGTools.WMass.postprocessing.framework.datamodel import Collection
from CMGTools.WMass.postprocessing.framework.eventloop import Module
from PhysicsTools.HeppyCore.utils.deltar import deltaPhi
from math import *

class SimpleVBoson:
    def __init__(self,legs,pdgId=24):
        self.pdgId=pdgId
        self.legs = legs
        if len(legs)<2:
            print "ERROR: making a VBoson w/ < 2 legs!"
        self.pt1 = legs[0].Pt()
        self.pt2 = legs[1].Pt()
        self.dphi = self.legs[0].DeltaPhi(self.legs[1])        
        self.deta = self.legs[0].Eta()-self.legs[1].Eta()
        self.px1 = legs[0].Px(); self.py1 = legs[0].Py();
        self.px2 = legs[1].Px(); self.py2 = legs[1].Py();
    def pt(self):
        return (self.legs[0]+self.legs[1]).Pt()
    def phi(self):
        return (self.legs[0]+self.legs[1]).Phi()
    def y(self):
        return (self.legs[0]+self.legs[1]).Rapidity()
    def mt(self):
        return sqrt(2*self.pt1*self.pt2*(1-cos(self.dphi)))
    def ux(self):
        return (-self.px1-self.px2)
    def uy(self):
        return (-self.py1-self.py2)    
    def mll(self):
        return sqrt(2*self.pt1*self.pt2*(cosh(self.deta)-cos(self.dphi)))

class KinematicVars:
    def __init__(self,beamE=6500):
        self.beamE = beamE
    def CSFrame(self,dilepton):
        pMass = 0.938272
        sign = np.sign(dilepton.Z())
        proton1 = ROOT.TLorentzVector(0.,0.,sign*self.beamE,hypot(self.beamE,pMass))  
        proton2 = ROOT.TLorentzVector(0.,0.,-sign*self.beamE,hypot(self.beamE,pMass))
        proton1.Boost(-dilepton.BoostVector()) 
        proton2.Boost(-dilepton.BoostVector())
        CSAxis = (proton1.Vect().Unit()-proton2.Vect().Unit()).Unit()
        yAxis = (proton1.Vect().Unit()).Cross((proton2.Vect().Unit()));
        yAxis = yAxis.Unit();
        xAxis = yAxis.Cross(CSAxis);
        xAxis = xAxis.Unit();
        return (xAxis,yAxis,CSAxis)
    def cosThetaCS(self,lplus,lminus):
        dilep = lplus + lminus
        boostedLep = ROOT.TLorentzVector(lminus)
        boostedLep.Boost(-dilep.BoostVector())
        csframe = self.CSFrame(dilep)
        return cos(boostedLep.Angle(csframe[2]))
    def cosTheta2D(self,w,lep):
        boostedLep = ROOT.TLorentzVector(lep)
        neww = ROOT.TLorentzVector()
        neww.SetPxPyPzE(w.Px(),w.Py(),0.,w.E())
        boostedLep.Boost(-neww.BoostVector())
        cost2d = (boostedLep.Px()*w.Px() + boostedLep.Py()*w.Py()) / (boostedLep.Pt()*w.Pt())
        return cost2d
    def cosThetaCM(self,lplus,lminus):
        dilep = lplus + lminus
        boostedLep = ROOT.TLorentzVector(lminus)
        boostedLep.Boost(-dilep.BoostVector())
        modw = sqrt(dilep.X()*dilep.X() + dilep.Y()*dilep.Y() + dilep.Z()*dilep.Z())
        modm = sqrt(boostedLep.X()*boostedLep.X() + boostedLep.Y()*boostedLep.Y() + boostedLep.Z()*boostedLep.Z())
        cos = (dilep.X()*boostedLep.X() + dilep.Y()*boostedLep.Y() + dilep.Z()*boostedLep.Z())/modw/modm
        return cos
    def phiCS(self,lplus,lminus):
        dilep = lplus + lminus
        boostedLep = ROOT.TLorentzVector(lminus)
        boostedLep.Boost(-dilep.BoostVector())
        csframe = self.CSFrame(dilep)
        phi = atan2((boostedLep.Vect()*csframe[1]),(boostedLep.Vect()*csframe[0]))
        if(phi<0): return phi + 2*ROOT.TMath.Pi()
        else: return phi

class GenQEDJetProducer(Module):
    def __init__(self,deltaR,beamEn=7000.,saveVertexZ=False,noSavePreFSR=False):
        self.beamEn=beamEn
        self.saveVertexZ=saveVertexZ
        self.noSavePreFSR=noSavePreFSR
        self.deltaR = deltaR
        self.vars = ("pt","eta","phi","mass","pdgId","fsrDR","fsrPt")
        self.genwvars = ("charge","pt","mass","y","costcs","phics","decayId")
        ## if "genQEDJetHelper_cc.so" not in ROOT.gSystem.GetLibraries():
        ##     print "Load C++ Worker"
        ##     ROOT.gROOT.ProcessLine(".L %s/src/CMGTools/WMass/python/postprocessing/helpers/genQEDJetHelper.cc+" % os.environ['CMSSW_BASE'])
        ## else:
        ##     print "genQEDJetHelper_cc.so found in ROOT libraries"
        ## self._worker = ROOT.GenQEDJetHelper(deltaR)
        self.pdfWeightOffset = 9 #index of first mc replica weight (careful, this should not be the nominal weight, which is repeated in some mc samples).  The majority of run2 LO madgraph_aMC@NLO samples with 5fs matrix element and pdf would use index 10, corresponding to pdf set 263001, the first alternate mc replica for the nominal pdf set 263000 used for these samples
        self.nMCReplicasWeights = 100 #number of input weights (100 for NNPDF 3.0)
        self.nHessianWeights    = 100 #number of hessian weights in the MC
        self.massWeights = range(-20,21) ## order them by integer range(80300, 80505, 5) #masses in MeV
        if "PDFWeightsHelper_cc.so" not in ROOT.gSystem.GetLibraries():
            ROOT.gROOT.ProcessLine(".include /cvmfs/cms.cern.ch/slc6_amd64_gcc530/external/eigen/3.2.2/include")
            ROOT.gROOT.ProcessLine(".L %s/src/CMGTools/WMass/python/postprocessing/helpers/PDFWeightsHelper.cc+" % os.environ['CMSSW_BASE'])
        mc2hessianCSV = "%s/src/CMGTools/WMass/python/postprocessing/data/gen/NNPDF30_nlo_as_0118_hessian_60.csv" % os.environ['CMSSW_BASE']
        self._pdfHelper = ROOT.PDFWeightsHelper()
        self._pdfHelper.Init(self.nMCReplicasWeights,self.nHessianWeights,mc2hessianCSV)
    def beginJob(self):
        pass
    def endJob(self):
        pass
    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.initReaders(inputTree) # initReaders must be called in beginFile
        self.out = wrappedOutputTree
        self.out.branch("pdfCentralWgt"  , "F")
        self.out.branch("partonId1"      , "I")
        self.out.branch("partonId2"      , "I")
        self.out.branch("nGenLepDressed" , "I")
        self.out.branch("nGenLepBare"    , "I")
        self.out.branch("nGenPromptNu"   , "I")
        for V in self.vars:
            self.out.branch("GenLepDressed_"+V, "F", lenVar="nGenLepDressed")
            if not self.noSavePreFSR: self.out.branch("GenLepPreFSR_" +V, "F", lenVar="nGenLepPreFSR")
            self.out.branch("GenLepBare_"   +V, "F", lenVar="nGenLepBare")
            self.out.branch("GenPromptNu_"  +V, "F", lenVar="nGenPromptNu")
        if self.saveVertexZ:
            self.out.branch("GenLepDressed_vertexZ", "F", lenVar="nGenLepDressed")
            self.out.branch("GenLepBare_vertexZ", "F", lenVar="nGenLepDressed")

        for V in self.genwvars:
            self.out.branch("genw_"   +V, "F")
            if not self.noSavePreFSR: self.out.branch("prefsrw_"+V, "F")
        for N in range(1,self.nHessianWeights+1):
            self.out.branch("hessWgt"+str(N), "H")
        for scale in ['muR','muF',"muRmuF","alphaS"]:
            for idir in ['Up','Dn']:
                self.out.branch("qcd_{scale}{idir}".format(scale=scale,idir=idir), "H")
        for imass in self.massWeights:
            masssign = 'm' if imass < 0 else 'p' if imass > 0 else ''
            self.out.branch("mass_{s}{mass}".format(s=masssign,mass=abs(imass)), "F")
    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass

    def initReaders(self,tree): # this function gets the pointers to Value and ArrayReaders and sets them in the C++ worker class
        ## marc try:
        ## marc     self.nGenPart = tree.valueReader("nGenPart")
        ## marc     for B in ("pt","eta","phi","mass","pdgId","isPromptHard","motherId","status") : setattr(self,"GenPart_"+B, tree.arrayReader("GenPart_"+B))
        ## marc     self._worker.setGenParticles(self.nGenPart,self.GenPart_pt,self.GenPart_eta,self.GenPart_phi,self.GenPart_mass,self.GenPart_pdgId,self.GenPart_isPromptHard,self.GenPart_motherId,self.GenPart_status)
        ## marc     #self.nLHEweight = tree.valueReader("nLHEweight")
        ## marc     #self.LHEweight_wgt = tree.arrayReader("LHEweight_wgt")
        ## marc except:
        ## marc     print '[genFriendProducer][Warning] Unable to attach to generator-level particles (data only?). No info will be produced'
        self._ttreereaderversion = tree._ttreereaderversion # self._ttreereaderversion must be set AFTER all calls to tree.valueReader or tree.arrayReader

    def mcRep2Hess(self,lheweights):
        nomlhew = lheweights[0]
        inPdfW = np.array(lheweights[self.pdfWeightOffset:self.pdfWeightOffset+self.nMCReplicasWeights],dtype=float)
        outPdfW = np.array([0 for i in xrange(self.nHessianWeights)],dtype=float)
        self._pdfHelper.DoMC2Hessian(nomlhew,inPdfW,outPdfW)
        pdfeigweights = [wgt/nomlhew for wgt in outPdfW.tolist()]
        return pdfeigweights

    def qcdScaleWgtIdx(self,mur="0",muf="0"):
        # mapping from https://indico.cern.ch/event/459797/contribution/2/attachments/1181555/1800214/mcaod-Feb15-2016.pdf
        idx_map={}
        idx_map[("0" ,"0" )] = 0
        idx_map[("0" ,"Up")] = 1
        idx_map[("0" ,"Dn")] = 2
        idx_map[("Up","0" )] = 3
        idx_map[("Up","Up")] = 4
        idx_map[("Up","Dn")] = 5
        idx_map[("Dn","0" )] = 6
        idx_map[("Dn","Up")] = 7
        idx_map[("Dn","Dn")] = 8
        if (mur,muf) not in idx_map: raise Exception('Scale variation muR={mur},muF={muf}'.format(mur=mur,muf=muf))
        return idx_map[(mur,muf)]

    def bwWeight(self,genMass,imass, isW):
        # mass and width from MiNNLO samples, using PDG inputs + width-dependent scheme https://indico.cern.ch/event/908015/contributions/3820269
        (m0,gamma) = (80351.972,2084.299) if isW else (91153.481, 2494.266) # MeV 
        newmass = m0 + imass*5.
        s_hat = pow(genMass,2)
        return (pow(s_hat - m0*m0,2) + pow(gamma*m0,2)) / (pow(s_hat - newmass*newmass,2) + pow(gamma*newmass,2))

    def getNeutrino(self):
        nus = []
        for p in self.genParts:
            if abs(p.pdgId) in [12, 14, 16] : ##and p.isPromptFinalState > 0: ## requiring isPromptFinalState results in 12% of events without neutrino
                nus.append(p)

        ##nus = sorted(nus, key = lambda x: x.pt)
        nus.sort(key = lambda x: x.pt, reverse=True)

        if len(nus) < 1:
            return 0, 0

        nu = nus[0]
        nuvec = ROOT.TLorentzVector()
        nuvec.SetPtEtaPhiM(nu.pt, nu.eta, nu.phi, nu.mass if nu.mass > 0. else 0.)

        return nuvec, nu.pdgId

    def getBareLeptons(self, strictlyPrompt=True):
        leptons = []
        for p in self.genParts:
            if not abs(p.pdgId) in [11,13,15] : continue
            if strictlyPrompt:
                if abs(p.pdgId) in [11, 13] and not p.isPromptFinalState : continue
                if abs(p.pdgId) == 15       and not p.isPromptDecayed : continue
            lepton = ROOT.TLorentzVector()
            lepton.SetPtEtaPhiM(p.pt, p.eta, p.phi, p.mass if p.mass >= 0. else 0.)
            if self.saveVertexZ:
                leptons.append( [lepton, p.pdgId, p.vertexZ] ) # only works if vertexZ was saved in genParts
            else:
                leptons.append( [lepton, p.pdgId, 0] ) 
            

        #leptons = sorted(leptons, key = lambda x: x[0].Pt() )
        leptons.sort(key = lambda x: x[0].Pt(), reverse=True )
        return leptons

    def getDressedLeptons(self, bareLeptons, strictlyPrompt=True, cone=0.1):

        if len(bareLeptons) ==0: 
            return bareLeptons
            ##print 'ERROR: DID NOT FIND A LEPTON FOR DRESSING !!!!!'

        leptons = copy.deepcopy(bareLeptons)
        for p in self.genParts:
            if not abs(p.pdgId)==22: continue
            if strictlyPrompt and not p.isPromptFinalState: continue
            tmp_photon = ROOT.TLorentzVector()
            tmp_photon.SetPtEtaPhiM(p.pt, p.eta, p.phi, p.mass if p.mass >= 0. else 0.)

            #first match in cone
            for il in xrange(0,len(leptons)):
                if leptons[il][0].DeltaR(tmp_photon) > cone: continue
                leptons[il][0] = leptons[il][0]+tmp_photon
                break

        ##print 'this is dressed leptons', leptons

        return leptons
    
    def getPreFSRLepton(self):
        
        lepInds = []
        for ip,p in enumerate(self.genParts):
            if abs(p.pdgId) in [11, 13] and p.fromHardProcessFinalState:
                lepInds.append( (p, ip) )
            if abs(p.pdgId) == 15       and p.fromHardProcessDecayed:
                lepInds.append( (p, ip) )
            ##if abs(p.pdgId) == 22:
            ##    if p.motherIndex >=0:
            ##        print 'photon with mother', self.genParts[int(p.motherIndex)].pdgId, 'and status', p.status

        if len(lepInds) < 1:
            return []

        retLeps = []

        for (lep, ind) in lepInds:
            finallep = lep

            ## see if there's a mother with same pdgId and lastBeforeFSR
            while (self.genParts[int(lep.motherIndex)].pdgId == lep.pdgId):
                if self.genParts[int(lep.motherIndex)].lastBeforeFSR:
                    finallep = self.genParts[int(lep.motherIndex)]
                    break
                ind = int(lep.motherIndex)
                lep = self.genParts[ind]
                if lep.motherIndex < 0:
                    break

            lepton = ROOT.TLorentzVector()
            lepton.SetPtEtaPhiM(finallep.pt, finallep.eta, finallep.phi, finallep.mass)

            retLeps.append( (lepton, finallep.pdgId) )

        ## if len(retLeps) == 4:
        ##     print '============'
        ##     for (lep, ind) in lepInds:
        ##         print 'this is pt {pt:4.2f} pdg id {pdg} motherpdgid {mid}'.format(pt=lep.pt,pdg= lep.pdgId, mid= self.genParts[int(lep.motherIndex)].pdgId)

        ##     print 'below the original leptons'
        ##     for ip,p in enumerate(self.genParts):
        ##         if abs(p.pdgId) in [11, 13]:
        ##             print 'pt {pt:4.2f} eta {eta:4.2f} phi {phi:4.2f} status {st}'.format(pt= p.pt, eta=p.eta, phi=p.phi, st=p.status)
        ##             lepInds.append( (p, ip) )
        ##     print '================='

        ##     print 'invariant mass of the first two leptons:', (retLeps[0][0]+retLeps[1][0]).M()
        ##     print 'invariant mass of the last  two leptons:', (retLeps[2][0]+retLeps[3][0]).M()
            

        return retLeps

    def getFSRPhotons(self, lepton, leptonPdgId, cone=10):

        photons = []
        for p in self.genParts:
            if abs(p.pdgId) != 22 : continue
            if int(p.motherIndex) < 0 : continue
            if self.genParts[int(p.motherIndex)].pdgId != leptonPdgId: continue 
            if not p.isPromptFinalState : continue
            photon = ROOT.TLorentzVector()
            photon.SetPtEtaPhiM(p.pt, p.eta, p.phi, p.mass if p.mass >= 0. else 0.)
            if photon.DeltaR(lepton) < cone: 
                photons.append(photon)

        photons.sort(key = lambda x: x.Pt(), reverse=True )

        return photons

    def analyze(self, event):
        """process event, return True (go to next module) or False (fail, go to next event)"""
        #nothing to do if this is data
        if event.isData: return True

        try:
            lhe_wgts      = Collection(event, 'LHEweight')
        except:
            lhe_wgts = 0
        self.genParts = Collection(event, 'GenPart')

        #this is a stupid check if variables are present. do this on the first gen particle
        makeBosons = False
        for ip,p in enumerate(self.genParts):
            if ip:
                break
            if hasattr(p, 'isPromptFinalState'):
                makeBosons = True
            

        if event._tree._ttreereaderversion > self._ttreereaderversion: # do this check at every event, as other modules might have read further branches
            self.initReaders(event._tree)

        ## check if the proper Ws can be made
        bareLeptonCollection                = self.getBareLeptons(strictlyPrompt=makeBosons)
        dressedLeptonCollection             = self.getDressedLeptons(bareLeptonCollection,strictlyPrompt=makeBosons)
        dressedLepton, dressedLeptonPdgId, dressedLeptonVertexZ   = dressedLeptonCollection[0] if len(dressedLeptonCollection) else (0,0,0)
        bareLepton,    bareLeptonPdgId, bareLeptonVertexZ      = bareLeptonCollection[0]    if len(bareLeptonCollection)    else (0,0,0)
        #(preFSRLepton , preFSRLeptonPdgId ) = self.getPreFSRLepton()     if makeBosons else (0,0)
        (neutrino     , neutrinoPdgId     ) = self.getNeutrino()         if makeBosons else (0,0)
        ## new to work also for Z
        if not self.noSavePreFSR:
            preFSRLeptonsAndPdgIds = self.getPreFSRLepton() if makeBosons else []
            # ## take the two hardest preFSR leptons if there are 4, which happens with photos
            if len(preFSRLeptonsAndPdgIds) == 4:
                preFSRLeptonsAndPdgIds = preFSRLeptonsAndPdgIds[:2]

        ## keep commented to save all particles
        #
        # if len(dressedLeptonCollection) == 4:
        #     dressedLeptonCollection = dressedLeptonCollection[:2]
        # ## =================================


        #always produce the dressed lepton collection
        if len(dressedLeptonCollection):
            self.out.fillBranch("nGenLepDressed", len(dressedLeptonCollection))
            retL={}
            #leptonsToTake = range(0,1 if makeBosons else len(dressedLeptonCollection))
            retL["pt"]    = [leppy[0].Pt()  for leppy in dressedLeptonCollection ]
            retL["eta"]   = [leppy[0].Eta() for leppy in dressedLeptonCollection ]
            retL["phi"]   = [leppy[0].Phi() for leppy in dressedLeptonCollection ]
            retL["pdgId"] = [leppy[1]       for leppy in dressedLeptonCollection ]
            retL["mass"]  = [leppy[0].M() if leppy[0].M() >= 0. else 0. ]
            retL["vertexZ"] = [leppy[2]       for leppy in dressedLeptonCollection ]

            retL["fsrDR"] = []
            retL["fsrPt"] = []
            for  leppy in dressedLeptonCollection:
                photonsFSR = self.getFSRPhotons(leppy[0],leppy[1])
                retL["fsrDR"].append(photonsFSR[0].DeltaR(leppy[0]) if len(photonsFSR)>0 else -1)
                retL["fsrPt"].append(photonsFSR[0].Pt() if len(photonsFSR)>0 else -1)

            for V in self.vars:
                self.out.fillBranch("GenLepDressed_"+V, retL[V])
            if self.saveVertexZ:
                self.out.fillBranch("GenLepDressed_vertexZ", retL["vertexZ"])

        #always produce the bare lepton collection
        if len(bareLeptonCollection)>0:
            retB={}
            #leptonsToTake = range(0,1 if makeBosons else len(bareLeptonCollection))
            leptonsToTake = range(0,len(bareLeptonCollection))
            retB["pt"]    = [bareLeptonCollection[i][0].Pt() for i in leptonsToTake ]
            retB["eta"]   = [bareLeptonCollection[i][0].Eta() for i in leptonsToTake ]
            retB["phi"]   = [bareLeptonCollection[i][0].Phi() for i in leptonsToTake ]
            retB["mass"]  = [bareLeptonCollection[i][0].M() if dressedLeptonCollection[i][0].M() >= 0. else 0. for i in leptonsToTake]
            retB["pdgId"] = [bareLeptonCollection[i][1] for i in leptonsToTake ]
            #retB["vertexZ"] = [bareLeptonCollection[i][2] for i in leptonsToTake ]
            retB["fsrDR"] = []
            retB["fsrPt"] = []
            for  i in leptonsToTake:
                photonsFSR = self.getFSRPhotons(bareLeptonCollection[i][0],bareLeptonCollection[i][1])
                retB["fsrDR"].append(photonsFSR[0].DeltaR(bareLeptonCollection[i][0]) if len(photonsFSR)>0 else -1)
                retB["fsrPt"].append(photonsFSR[0].Pt() if len(photonsFSR)>0 else -1)

            self.out.fillBranch("nGenLepBare", leptonsToTake[-1])
            for V in self.vars:
                self.out.fillBranch("GenLepBare_"+V, retB[V])
            if self.saveVertexZ:
                self.out.fillBranch("GenLepBare_vertexZ", retL["vertexZ"])

        #W-specific
        if neutrino:
            retN={}
            retN["pt"]    = [neutrino.Pt()  ]
            retN["eta"]   = [neutrino.Eta() ]
            retN["phi"]   = [neutrino.Phi() ]
            retN["mass"]  = [neutrino.M() if neutrino.M() >= 0. else 0.  ] ## weird protection... also already checked. but better safe than sorry
            retN["pdgId"] = [neutrinoPdgId  ]
            retN["fsrDR"] = [-1 ]
            retN["fsrPt"] = [-1 ]
            self.out.fillBranch("nGenPromptNu", 1)
            for V in self.vars:
                self.out.fillBranch("GenPromptNu_"+V, retN[V])
            #self.out.fillBranch("GenPromptNu_pdgId", [pdgId for pdgId in nuPdgIds])            

        if not self.noSavePreFSR:
            if len(preFSRLeptonsAndPdgIds):
                self.out.fillBranch("nGenLepPreFSR", len(preFSRLeptonsAndPdgIds))
                retP={}
                retP["pt"]    = [leppy[0].Pt()  for leppy in preFSRLeptonsAndPdgIds]
                retP["eta"]   = [leppy[0].Eta() for leppy in preFSRLeptonsAndPdgIds]
                retP["phi"]   = [leppy[0].Phi() for leppy in preFSRLeptonsAndPdgIds]
                retP["pdgId"] = [leppy[1]       for leppy in preFSRLeptonsAndPdgIds]
                retP["mass"]  = [leppy[0].M()   if leppy[0].M() >= 0. else 0. for leppy in preFSRLeptonsAndPdgIds]
                retP["fsrDR"] = []
                retP["fsrPt"] = []
                for leppy in preFSRLeptonsAndPdgIds:
                    photonsFSR = self.getFSRPhotons(leppy[0],leppy[1])
                    retP["fsrDR"].append(photonsFSR[0].DeltaR(leppy[0]) if len(photonsFSR)>0 else -1 )
                    retP["fsrPt"].append(photonsFSR[0].Pt() if len(photonsFSR)>0 else -1 )

                for V in self.vars:
                    self.out.fillBranch("GenLepPreFSR_"+V, retP[V])


            # this part is to be revisited if I change the number of saved GenLepDressed or bare
            # there should be 2 or 4 (for photos), when there are 4 the first 2 are saved
            # but I might want to save all
            #if len(preFSRLeptonsAndPdgIds) == 2 or (len(preFSRLeptonsAndPdgIds) == 1 and neutrino):
            if len(preFSRLeptonsAndPdgIds):


                if not neutrino: ## these are Zs
                    isWBoson = False
                    prefsrw = preFSRLeptonsAndPdgIds[0][0] + preFSRLeptonsAndPdgIds[1][0]
                    self.out.fillBranch("prefsrw_charge" , 0)
                    l1, l1pdg, l2 = preFSRLeptonsAndPdgIds[0][0], preFSRLeptonsAndPdgIds[0][1], preFSRLeptonsAndPdgIds[1][0]
                    # convention for phiCS: use l- direction for W-, use neutrino for W+
                    (lplus,lminus) = (l2,l1) if l1pdg<0 else (l1,l2)
                    self.out.fillBranch("prefsrw_decayId", abs(l1pdg))

                else: ## these are Ws
                    isWBoson = True
                    prefsrw = preFSRLeptonsAndPdgIds[0][0] + neutrino
                    self.out.fillBranch("prefsrw_charge" , float(-1*np.sign(preFSRLeptonsAndPdgIds[0][1])))
                    l1, l1pdg = preFSRLeptonsAndPdgIds[0][0], preFSRLeptonsAndPdgIds[0][1]
                    # convention for phiCS: use l- direction for W-, use neutrino for W+
                    (lplus,lminus) = (neutrino,l1) if l1pdg<0 else (l1,neutrino)
                    self.out.fillBranch("prefsrw_decayId", abs(neutrinoPdgId))

                self.out.fillBranch("prefsrw_pt"     , prefsrw.Pt())
                self.out.fillBranch("prefsrw_y"      , prefsrw.Rapidity())
                self.out.fillBranch("prefsrw_mass"   , prefsrw.M())
                kv = KinematicVars(self.beamEn)
                self.out.fillBranch("prefsrw_costcs" , kv.cosThetaCS(lplus, lminus))
                self.out.fillBranch("prefsrw_phics"  , kv.phiCS     (lplus, lminus))
                for imass in self.massWeights:
                    masssign = 'm' if imass < 0 else 'p' if imass > 0 else ''
                    self.out.fillBranch("mass_{s}{mass}".format(s=masssign,mass=abs(imass)), self.bwWeight(genMass=prefsrw.M()*1000,imass=imass,isW=isWBoson))

        #if len(dressedLeptonCollection) == 2 or (len(dressedLeptonCollection) == 1 and neutrino):
        if len(dressedLeptonCollection):
            if not neutrino: ## these are Zs
                isWBoson = False
                genw = dressedLeptonCollection[0][0] + dressedLeptonCollection[1][0]
                self.out.fillBranch("genw_charge" , 0)
                l1, l1pdg, l2 = dressedLeptonCollection[0][0], dressedLeptonCollection[0][1], dressedLeptonCollection[1][0]
                # convention for phiCS: use l- direction for W-, use neutrino for W+
                (lplus,lminus) = (l2,l1) if l1pdg<0 else (l1,l2)
                self.out.fillBranch("genw_decayId", abs(l1pdg))

            else: ## these are Ws
                isWBoson = True
                genw = dressedLeptonCollection[0][0] + neutrino
                self.out.fillBranch("genw_charge" , float(-1*np.sign(dressedLeptonCollection[0][1])))
                l1, l1pdg = dressedLeptonCollection[0][0], dressedLeptonCollection[0][1]
                # convention for phiCS: use l- direction for W-, use neutrino for W+
                (lplus,lminus) = (neutrino,l1) if l1pdg<0 else (l1,neutrino)
                self.out.fillBranch("genw_decayId", abs(neutrinoPdgId))

            self.out.fillBranch("genw_pt"     , genw.Pt())
            self.out.fillBranch("genw_y"      , genw.Rapidity())
            self.out.fillBranch("genw_mass"   , genw.M())
            kv = KinematicVars(self.beamEn)
            self.out.fillBranch("genw_costcs" , kv.cosThetaCS(lplus, lminus))
            self.out.fillBranch("genw_phics"  , kv.phiCS     (lplus, lminus))
            for imass in self.massWeights:
                masssign = 'm' if imass < 0 else 'p' if imass > 0 else ''
                self.out.fillBranch("mass_{s}{mass}".format(s=masssign,mass=abs(imass)), self.bwWeight(genMass=genw.M()*1000,imass=imass,isW=isWBoson))

            ## remove these variables to save space
            #self.out.fillBranch("genw_mt"   , sqrt(2*lplus.Pt()*lminus.Pt()*(1.-cos(deltaPhi(lplus.Phi(),lminus.Phi())) )))
            #self.out.fillBranch("genw_eta",genw.Eta())
            #self.out.fillBranch("genw_phi",genw.Phi())
            #self.out.fillBranch("genw_costcm",kv.cosThetaCM(lplus,lminus))
            #self.out.fillBranch("genw_cost2d",kv.cosTheta2D(genw,dressedLeptons[0]))
        else: ## no neutrino found!
            ##if not len(dressedLeptons): 
            ##    print '================================'
            ##    print 'no dressed leptons found!'
            ##    print 'no dressed leptons fround :lumi:evt: {a}:{b}:{c}'.format(a=getattr(event, "run"),b=getattr(event, "lumi"),c=getattr(event, "evt"))
            ##else:
            ##    print '================================'
            ##    print 'no neutrinos found, in run:lumi:evt: {a}:{b}:{c}'.format(a=getattr(event, "run"),b=getattr(event, "lumi"),c=getattr(event, "evt"))
            for V in self.genwvars:
                self.out.fillBranch("genw_"+V   , -999)
                if not self.noSavePreFSR: self.out.fillBranch("prefsrw_"+V, -999)
            for imass in self.massWeights:
                masssign = 'm' if imass < 0 else 'p' if imass > 0 else ''
                self.out.fillBranch("mass_{s}{mass}".format(s=masssign,mass=abs(imass)), 1.)

        if hasattr(event,"genWeight"):
            self.out.fillBranch("partonId1" , getattr(event, "id1") if hasattr(event, 'id1') else -999)
            self.out.fillBranch("partonId2" , getattr(event, "id2") if hasattr(event, 'id2') else -999)
        else:
            self.out.fillBranch("partonId1" , -999 )
            self.out.fillBranch("partonId2" , -999 )

        centralPDFWeight = 1.
        if lhe_wgts:
            lheweights = [w.wgt for w in lhe_wgts]

            ## set here the nominal weight for the PDF set we want to use:
            ## the index should be:
            ##  18 for NNPDF31_nnlo_hessian_pdfas
            ## 129 for NNPDF31_nnlo_as_0118_CMSW1_hessian_100, where: CMSW1 => No CMS W data
            ## 230 for NNPDF31_nnlo_as_0118_CMSW2_hessian_100, where: CMSW2 => No collider W data
            ## 331 for NNPDF31_nnlo_as_0118_CMSW3_hessian_100, where: CMSW3 => No CMS W,Z data
            ## 432 for NNPDF31_nnlo_as_0118_CMSW4_hessian_100, where: CMSW4 => No collider W,Z data
            ## 533 for CT14nnlo             ATTENTION, number of hessianWeights is: 29
            ## 592 for MMHT2014nnlo68cl     ATTENTION, number of hessianWeights is: 51, alphaS separate
            ## 646 for ABMP16_5_nnlo        ATTENTION, number of hessianWeights is: 30
            ## 676 for HERAPDF20_NNLO_EIG   ATTENTION, number of hessianWeights is: 29
            ## 705 for HERAPDF20_NNLO_VAR   ATTENTION, number of hessianWeights is: 14, alphaS separate

            centralPDFIndex  = 432 ## for CMSW4
            centralPDFWeight = lheweights[centralPDFIndex] ## 423 for CMSW4

            ## the PDF variations are saved as ratios w/r/t the nominal weight of the selected sample
            ## i.e. running a template with pdf6 up will be pdfCentralWgt * hessWgt6,
            ## while the nominal template will be just pdfCentralWgt
            for N in range(1,self.nHessianWeights+1):
                self.out.fillBranch("hessWgt"+str(N), lheweights[centralPDFIndex+N]/centralPDFWeight if centralPDFWeight else 0.)

            qcd0Wgt=lheweights[self.qcdScaleWgtIdx()]
            for ii,idir in enumerate(['Up','Dn']):
                self.out.fillBranch("qcd_muR{idir}"   .format(idir=idir), lheweights[self.qcdScaleWgtIdx(mur=idir)]/qcd0Wgt)
                self.out.fillBranch("qcd_muF{idir}"   .format(idir=idir), lheweights[self.qcdScaleWgtIdx(muf=idir)]/qcd0Wgt)
                self.out.fillBranch("qcd_muRmuF{idir}".format(idir=idir), lheweights[self.qcdScaleWgtIdx(mur=idir,muf=idir)]/qcd0Wgt) # only correlated variations are physical
                
                ## alphaS now different in new MC. CMSW4 does not have a alphaS variation, so taking the one of the nominal PDF set
                self.out.fillBranch("qcd_alphaS{idir}".format(idir=idir), lheweights[110+ii] ) ## 110 is alphaS-up of nominal pdf set, 111 is alphaS-down of nominal pdf set

        ## save the central PDF weight that is either 1. or set in the previous block of code
        self.out.fillBranch("pdfCentralWgt", centralPDFWeight )

        return True

# define modules using the syntax 'name = lambda : constructor' to avoid having them loaded when not needed

genQEDJets14TeV = lambda : GenQEDJetProducer(deltaR=0.1,beamEn=7000.)
genQEDJets = lambda : GenQEDJetProducer(deltaR=0.1,beamEn=6500.,noSavePreFSR=True)
genQEDJetsWithVertex = lambda : GenQEDJetProducer(deltaR=0.1,beamEn=6500.,saveVertexZ=True, noSavePreFSR=True)

