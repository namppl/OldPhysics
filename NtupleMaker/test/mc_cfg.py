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
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff")
process.GlobalTag.globaltag = cms.string("76X_mcRun2_asymptotic_v12") 
process.load("TrackingTools/TransientTrack/TransientTrackBuilder_cfi")

# load the PAT config
#process.load("PhysicsTools.PatAlgos.patSequences_cff")

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

useAOD = False

if useAOD == True :
    dataFormat = DataFormat.AOD
else :
    dataFormat = DataFormat.MiniAOD

switchOnVIDPhotonIdProducer(process, dataFormat)

my_id_modules = ['RecoEgamma.PhotonIdentification.Identification.cutBasedPhotonID_Spring15_25ns_V1_cff']

for idmod in my_id_modules:
    setupAllVIDIdsInModule(process,idmod,setupVIDPhotonSelection)

process.tree = cms.EDAnalyzer("NtupleMaker",
    isMC = cms.untracked.bool(True),
    TriggerResults = cms.InputTag("TriggerResults","","HLT"),
    TriggerObjects = cms.InputTag("selectedPatTrigger"),
    TriggerSummary = cms.InputTag("hltTriggerSummaryAOD"),
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

process.p = cms.Path(
    #process.noscraping*
    process.goodOfflinePrimaryVertices*
    process.egmPhotonIDSequence*
    process.tree
)