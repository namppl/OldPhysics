import FWCore.ParameterSet.Config as cms

isMC = True

process = cms.Process("NtupleProduction")

## MessageLogger
process.load("FWCore.MessageLogger.MessageLogger_cfi")

## Options and Output Report
process.options   = cms.untracked.PSet( 
  wantSummary = cms.untracked.bool(True) 
)
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

## Source
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
       #'root://eoscms//eos/cms/store/mc/RunIISpring15DR74/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/Asympt25ns_MCRUN2_74_V9-v3/10000/FEF074BD-9714-E511-8FE8-0025907DC9C4.root'
       'file:DY_MINIAOD.root'
    )
)

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100) )

## Geometry and Detector Conditions (needed for a few patTuple production steps)
process.load("Configuration.Geometry.GeometryRecoDB_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
from Configuration.AlCa.autoCond import autoCond
if isMC:
   process.GlobalTag.globaltag = autoCond['startup']
else:
   process.GlobalTag.globaltag = autoCond['com10']

process.load("Configuration.StandardSequences.MagneticField_cff")

## Output Module Configuration (expects a path 'p')
from PhysicsTools.PatAlgos.patEventContent_cff import patEventContent
process.out = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string('patTuple_skim.root'),
    splitLevel = cms.untracked.int32(0),
    # save only events passing the full path
    #SelectEvents   = cms.untracked.PSet( SelectEvents = cms.vstring('p') ),
    # save PAT Layer 1 output; you need a '*' to
    # unpack the list of commands 'patEventContent'
    outputCommands = cms.untracked.vstring('keep *_*_*_*')
)

# load the PAT config
process.load("PhysicsTools.PatAlgos.patSequences_cff")
#process.patJets.addTagInfos     = cms.bool(False)

process.TFileService = cms.Service("TFileService",
  fileName = cms.string('ntuple_skim.root')
)

process.goodOfflinePrimaryVertices = cms.EDFilter("VertexSelector",
   src = cms.InputTag("offlineSlimmedPrimaryVertices"),
   cut = cms.string("!isFake && ndof > 4 && abs(z) < 24 && position.Rho < 2"), # tracksSize() > 3 for the older cut
   filter = cms.bool(True),   # otherwise it won't filter the events, just produce an empty vertex collection.
)

process.noscraping = cms.EDFilter("FilterOutScraping",
   applyfilter = cms.untracked.bool(True),
   debugOn = cms.untracked.bool(False),
   numtrack = cms.untracked.uint32(10),
   thresh = cms.untracked.double(0.25)
)

from PhysicsTools.SelectorUtils.tools.vid_id_tools import *
# turn on VID producer, indicate data format  to be
# DataFormat.AOD or DataFormat.MiniAOD, as appropriate

useAOD = False

if useAOD == True :
    dataFormat = DataFormat.AOD
else :
    dataFormat = DataFormat.MiniAOD

switchOnVIDPhotonIdProducer(process, dataFormat)

# define which IDs we want to produce
my_id_modules = ['RecoEgamma.PhotonIdentification.Identification.cutBasedPhotonID_Spring15_25ns_V1_cff']

#add them to the VID producer
for idmod in my_id_modules:
    setupAllVIDIdsInModule(process,idmod,setupVIDPhotonSelection)

process.tree = cms.EDAnalyzer("NtupleMaker",
  isMC = cms.untracked.bool(True),
  TriggerResults = cms.InputTag("TriggerResults","","HLT"),
  TriggerObjects = cms.InputTag("selectedPatTrigger"),
  TriggerSummary = cms.InputTag("hltTriggerSummaryAOD"),
  isMINIAOD = cms.untracked.bool(True),
  Vertex = cms.InputTag("offlineSlimmedPrimaryVertices"),
  BeamSpot = cms.InputTag("offlineBeamSpot"),
  Muons = cms.InputTag("slimmedMuons"),
  Electrons = cms.InputTag("slimmedElectrons"),
  Conversions = cms.InputTag("reducedEgamma:reducedConversions"),
  GenParticles = cms.InputTag("prunedGenParticles"),
  Generator = cms.InputTag("generator"),
  Photons = cms.InputTag("slimmedPhotons"),
  rho = cms.InputTag("fixedGridRhoFastjetAll"),
  full5x5SigmaIEtaIEtaMap = cms.InputTag("photonIDValueMapProducer:phoFull5x5SigmaIEtaIEta"),
  phoChargedIsolation = cms.InputTag("photonIDValueMapProducer:phoChargedIsolation"),
  phoNeutralHadronIsolation = cms.InputTag("photonIDValueMapProducer:phoNeutralHadronIsolation"),
  phoPhotonIsolation = cms.InputTag("photonIDValueMapProducer:phoPhotonIsolation"),
  Jets = cms.InputTag("slimmedJets"),
  MET = cms.InputTag("slimmedMETs")
)

#process.recoTree.ElectronCalib = "calibratedPatElectrons"
#process.recoTree.Electron = "selectedPatElectrons"
#process.recoTree.Jet = "selectedPatJets"
#process.recoTree.MET = "patMETs"
#process.recoTree.FilterType = 1
#process.recoTree.IntLuminosity = 1.
#process.recoTree.CrossSection = 1666.
#process.recoTree.FilterEfficiency = 1.
#process.recoTree.TotalNevents = 1963296.

#process.load('PhysicsTools.PatAlgos.patSequences_cff')
#from SUSYBSMAnalysis.Zprime2muAnalysis.PATTools import pruneMCLeptons
#pruneMCLeptons(process, use_sim=False) # need to decide whether to move AODOnly() call in here, if so use_sim should just be set False

#from SUSYBSMAnalysis.Zprime2muAnalysis.PATTools import AODOnly
#AODOnly(process)

#process.load("RecoTauTag/Configuration/RecoPFTauTag_cff")
#from PhysicsTools.PatAlgos.tools.coreTools import *
#removeSpecificPATObjects( process, ['Taus'] )
#process.patDefaultSequence.remove( process.patTaus )

if isMC==False:
   removeMCMatching(process, ['All'])

#from RecoJets.JetProducers.kt4PFJets_cfi import kt4PFJets
#process.kt6PFJetsForIsolation = kt4PFJets.clone( rParam = 0.6, doRhoFastjet = True)
#process.kt6PFJetsForIsolation.Rho_EtaMax = cms.double(2.5)

# MET correction
#process.load("JetMETCorrections.Type1MET.pfMETCorrections_cff") # Type I
#process.load("JetMETCorrections.Type1MET.pfMETCorrections_cff")
#process.load("JetMETCorrections.Type1MET.pfMETCorrectionType0_cfi")
#process.pfType1CorrectedMet.applyType0Corrections = cms.bool(False)
#process.pfType1CorrectedMet.srcType1Corrections = cms.VInputTag(
#  cms.InputTag('pfMETcorrType0'),
#  cms.InputTag('pfJetMETcorr', 'type1')        
#)

# energy scale correction
#process.RandomNumberGeneratorService = cms.Service("RandomNumberGeneratorService",
#  calibratedPatElectrons = cms.PSet(
#    initialSeed = cms.untracked.uint32(1),
#    engineName = cms.untracked.string('TRandom3')
#  ),
#)

#process.load("EgammaAnalysis.ElectronTools.calibratedPatElectrons_cfi")
# dataset to correct
#process.calibratedPatElectrons.isMC = cms.bool(True)
#process.calibratedPatElectrons.inputDataset = cms.string("Summer12_LegacyPaper")
#process.calibratedPatElectrons.updateEnergyError = cms.bool(True)
#process.calibratedPatElectrons.correctionsType = cms.int32(2) 
#process.calibratedPatElectrons.combinationType = cms.int32(3)
#process.calibratedPatElectrons.verbose = cms.bool(False)
#process.calibratedPatElectrons.synchronization = cms.bool(False)

#process.load('EgammaAnalysis.ElectronTools.electronRegressionEnergyProducer_cfi')
#process.eleRegressionEnergy.energyRegressionType = cms.uint32(2)
#process.eleRegressionEnergy.inputElectronsTag = cms.InputTag('selectedPatElectrons')
#process.eleRegressionEnergy.produceValueMaps = cms.bool(True)

#process.eleIsoSequence = setupPFElectronIso(process, 'gsfElectrons')
#process.phoIsoSequence = setupPFPhotonIso(process, 'photons')

# boosted Z->ee isolation
#process.load("RecoLocalCalo.EcalRecAlgos.EcalSeverityLevelESProducer_cfi")
#process.heepIdNoIso = cms.EDProducer("HEEPIdValueMapProducer",
 #    eleLabel = cms.InputTag("gsfElectrons"),
  #   barrelCuts = cms.PSet(heepBarrelCutsZpSUSY),
   #  endcapCuts = cms.PSet(heepEndcapCutsZpSUSY),
    # eleIsolEffectiveAreas = cms.PSet(heepEffectiveAreas),
#     eleRhoCorrLabel = cms.InputTag("kt6PFJets", "rho"),
 #    verticesLabel = cms.InputTag("offlinePrimaryVertices"),
  #   applyRhoCorrToEleIsol = cms.bool(True),
   #  writeIdAsInt = cms.bool(True)
#)
#process.heepIdNoIso.barrelCuts.cuts=cms.string("et:detEta:ecalDriven:dEtaIn:dPhiIn:hadem:e2x5Over5x5:nrMissHits:dxy")
#process.heepIdNoIso.endcapCuts.cuts=cms.string("et:detEta:ecalDriven:dEtaIn:dPhiIn:hadem:sigmaIEtaIEta:nrMissHits:dxy")
#process.heepIdNoIsoEles = cms.EDProducer("tsw::HEEPGsfProducer", cutValueMap = cms.InputTag("heepIdNoIso"),
 #                                        inputGsfEles = cms.InputTag("gsfElectrons")  )

#from PhysicsTools.PatAlgos.patEventContent_cff import *
#process.out.outputCommands += patTriggerEventContent
#process.out.outputCommands += patExtraAodEventContent
#process.out.outputCommands += patEventContentNoCleaning
#process.out.outputCommands.extend(cms.untracked.vstring(
#  'keep *_*_*_*',
#))


#process.heepIdNoIsoOrig = cms.EDProducer("HEEPIdValueMapProducer",
   #  eleLabel = cms.InputTag("gsfElectrons"),
   #  barrelCuts = cms.PSet(heepBarrelCuts),
   #  endcapCuts = cms.PSet(heepEndcapCuts),
   #  eleIsolEffectiveAreas = cms.PSet(heepEffectiveAreas),
    # eleRhoCorrLabel = cms.InputTag("kt6PFJets", "rho"),
   #  verticesLabel = cms.InputTag("offlinePrimaryVertices"),
  #   applyRhoCorrToEleIsol = cms.bool(True),
 #    writeIdAsInt = cms.bool(True)
#)
#process.heepIdNoIsoOrig.barrelCuts.cuts=cms.string("et:detEta:ecalDriven:dEtaIn:dPhiIn:hadem:e2x5Over5x5:nrMissHits:dxy")
#process.heepIdNoIsoOrig.endcapCuts.cuts=cms.string("et:detEta:ecalDriven:dEtaIn:dPhiIn:hadem:sigmaIEtaIEta:nrMissHits:dxy")
#process.heepIdNoIsoOrigEles = cms.EDProducer("tsw::HEEPGsfProducer", cutValueMap = cms.InputTag("heepIdNoIsoOrig"),
 #                                        inputGsfEles = cms.InputTag("gsfElectrons")  )

# Boosted Z ModEleIso: 1b) Calculating the modified iso. values using BstdZeeTools EDProducer
#from TSWilliams.BstdZeeTools.bstdzeemodisolproducer_cff import *
#process.modElectronIso = cms.EDProducer("BstdZeeModIsolProducer", 
 #     bstdZeeModIsolParams, vetoGsfEles = cms.InputTag("heepIdNoIsoEles") )
#process.modElectronIsoOrig = cms.EDProducer("BstdZeeModIsolProducer", 
    #  bstdZeeModIsolParams, vetoGsfEles = cms.InputTag("heepIdNoIsoOrigEles") )

#process.patElectrons.userIsolation.user = cms.VPSet(
 #   cms.PSet(src = cms.InputTag("modElectronIso","track")),
  #  cms.PSet(src = cms.InputTag("modElectronIso","ecal")),
   # cms.PSet(src = cms.InputTag("modElectronIso","hcalDepth1")),
    #cms.PSet(src = cms.InputTag("modElectronIsoOrig","track")),
    #cms.PSet(src = cms.InputTag("modElectronIsoOrig","ecal")),
    #cms.PSet(src = cms.InputTag("modElectronIsoOrig","hcalDepth1")),
#)

# Let it run
process.p = cms.Path(
    #process.noscraping*
    process.goodOfflinePrimaryVertices*
    process.egmPhotonIDSequence*
    #process.PFTau*
    #process.kt6PFJetsForIsolation*
    #process.type0PFMEtCorrection*
    #process.producePFMETCorrections*
    #(process.heepIdNoIso+process.heepIdNoIsoEles+process.modElectronIso)*
    #(process.heepIdNoIsoOrig+process.heepIdNoIsoOrigEles+process.modElectronIsoOrig)*
    #process.eleRegressionEnergy*
    #process.calibratedPatElectrons*
    #(process.pfParticleSelectionSequence + process.eleIsoSequence + process.phoIsoSequence)*
    #process.patDefaultSequence*
    process.tree
)

#process.outpath = cms.EndPath(process.out)

