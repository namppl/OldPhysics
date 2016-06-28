import FWCore.ParameterSet.Config as cms


process = cms.Process("NtupleProduction")
isMC = True

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
       ''
    )
)

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

## Geometry and Detector Conditions (needed for a few patTuple production steps)
process.load("Configuration.Geometry.GeometryRecoDB_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff")
process.GlobalTag.globaltag = cms.string("80X_mcRun2_asymptotic_2016_v3")
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
switchOnVIDElectronIdProducer(process, dataFormat)

my_phoId_modules = ['RecoEgamma.PhotonIdentification.Identification.mvaPhotonID_Spring15_25ns_nonTrig_V2_cff']
my_eleId_modules = ['RecoEgamma.ElectronIdentification.Identification.cutBasedElectronID_Spring15_25ns_V1_cff']
#my_eleId_modules = ['RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Spring15_25ns_Trig_V1_cff']

for idmod in my_phoId_modules:
    setupAllVIDIdsInModule(process,idmod,setupVIDPhotonSelection)
for idmod in my_eleId_modules:
    setupAllVIDIdsInModule(process,idmod,setupVIDElectronSelection)

process.tree = cms.EDAnalyzer("NtupleMaker",
    isMC = cms.untracked.bool(True),
	putMuons = cms.untracked.bool(True),
	putElectrons = cms.untracked.bool(True),
	putDileptons = cms.untracked.bool(False),
	putPhotons = cms.untracked.bool(True),
	putJets = cms.untracked.bool(False),
	putMETs = cms.untracked.bool(False),
    TriggerResults = cms.InputTag("TriggerResults","","HLT"),
    TriggerObjects = cms.InputTag("selectedPatTrigger"),
    TriggerSummary = cms.InputTag("hltTriggerSummaryAOD"),
    Vertex = cms.InputTag("offlineSlimmedPrimaryVertices"),
    BeamSpot = cms.InputTag("offlineBeamSpot"),
    PFCandidates = cms.InputTag("packedPFCandidates"),
    Muons = cms.InputTag("slimmedMuons"),
    Electrons = cms.InputTag("slimmedElectrons"),
    Conversions = cms.InputTag("reducedEgamma:reducedConversions"),
    GenParticles = cms.InputTag("prunedGenParticles"),
    Generator = cms.InputTag("generator"),
    Photons = cms.InputTag("slimmedPhotons"),
    rho = cms.InputTag("fixedGridRhoFastjetAll"),
    miniRho = cms.InputTag("fixedGridRhoFastjetCentralNeutral"),
    phoMediumIdBoolMap = cms.InputTag("egmPhotonIDs:mvaPhoID-Spring15-25ns-nonTrig-V2-wp90"),
    phoMediumIdFullInfoMap = cms.InputTag("egmPhotonIDs:mvaPhoID-Spring15-25ns-nonTrig-V2-wp90"),
    mvaValuesMap     = cms.InputTag("photonMVAValueMapProducer:PhotonMVAEstimatorRun2Spring15NonTrig25nsV2Values"),
    mvaCategoriesMap = cms.InputTag("photonMVAValueMapProducer:PhotonMVAEstimatorRun2Spring15NonTrig25nsV2Categories"),
    full5x5SigmaIEtaIEtaMap = cms.InputTag("photonIDValueMapProducer:phoFull5x5SigmaIEtaIEta"),
    phoChargedIsolation = cms.InputTag("photonIDValueMapProducer:phoChargedIsolation"),
    phoNeutralHadronIsolation = cms.InputTag("photonIDValueMapProducer:phoNeutralHadronIsolation"),
    phoPhotonIsolation = cms.InputTag("photonIDValueMapProducer:phoPhotonIsolation"),
    Jets = cms.InputTag("slimmedJets"),
    MET = cms.InputTag("slimmedMETs"),
    eleMediumIdMap = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Spring15-25ns-V1-standalone-medium"),
    eleTightIdMap = cms.InputTag("egmGsfElectronIDs:cutBasedElectronID-Spring15-25ns-V1-standalone-tight")
    #eleMediumIdMap = cms.InputTag("egmGsfElectronIDs:mvaEleID-Spring15-25ns-Trig-V1-wp90"),
    #eleTightIdMap = cms.InputTag("egmGsfElectronIDs:mvaEleID-Spring15-25ns-Trig-V1-wp80"),
    #eleMvaValuesMap     = cms.InputTag("electronMVAValueMapProducer:ElectronMVAEstimatorRun2Spring15Trig25nsV1Values"),
    #eleMvaCategoriesMap = cms.InputTag("electronMVAValueMapProducer:ElectronMVAEstimatorRun2Spring15Trig25nsV1Categories")
)

process.p = cms.Path(
    #process.noscraping*
    process.goodOfflinePrimaryVertices*
    process.egmPhotonIDSequence*
    process.egmGsfElectronIDSequence*
    process.tree
)
