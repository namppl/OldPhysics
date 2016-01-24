// -*- C++ -*-
//
// Package:    Physics/NtupleMaker
// Class:      NtupleMaker
// 
/**\class NtupleMaker NtupleMaker.cc Physics/NtupleMaker/plugins/NtupleMaker.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Kyungwook Nam
//         Created:  Sat, 23 Jan 2016 15:57:13 GMT
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h" 
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "DataFormats/EgammaCandidates/interface/Conversion.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
#include "DataFormats/EgammaCandidates/interface/ConversionFwd.h"
#include "DataFormats/EgammaCandidates/interface/Conversion.h"
#include "DataFormats/EgammaCandidates/interface/Photon.h"
#include "DataFormats/EgammaCandidates/interface/PhotonFwd.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/PatCandidates/interface/Tau.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/METReco/interface/MET.h"
#include "DataFormats/METReco/interface/HcalNoiseSummary.h"
#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/JetReco/interface/PFJetCollection.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "JetMETCorrections/Objects/interface/JetCorrectionsRecord.h"
#include "JetMETCorrections/Objects/interface/JetCorrector.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "PhysicsTools/RecoUtils/interface/CandCommonVertexFitter.h"
#include "RecoEgamma/EgammaTools/interface/ConversionTools.h"
#include "RecoEgamma/EgammaTools/interface/EffectiveAreas.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"
#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"
#include "RecoVertex/VertexTools/interface/InvariantMassFromVertex.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/TransientTrack/interface/GsfTransientTrack.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"

#include <TTree.h>
#include <TROOT.h>
#include <TSystem.h>
#include <TFile.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TFrame.h>
#include <TMath.h>
#include <TF1.h>
#include <TLorentzVector.h>
#include <TVector3.h>
#include <boost/foreach.hpp>
#include <iostream>
#include <map>
#include <string>
#include <iomanip>
#include <cassert>
     
#include "NtupleMaker.h"

using namespace std;
using namespace reco;
using namespace edm;
using namespace pat;
using namespace isodeposit;


//
// class declaration
//

class NtupleMaker : public edm::one::EDAnalyzer<edm::one::SharedResources>  {

public:
    explicit NtupleMaker(const edm::ParameterSet&);
    ~NtupleMaker();

    static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


private:
    virtual void beginJob() override;
    virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
    virtual void endJob() override;

    virtual void fillTriggers(const edm::Event& iEvent);
    virtual void fillMuons(const edm::Event& iEvent);
    virtual void fillElectrons(const edm::Event& iEvent);
    virtual void fillDileptons(const edm::Event& iEvent);
    virtual void fillPhotons(const edm::Event& iEvent);
    virtual void fillJets(const edm::Event& iEvent);
    virtual void fillMETs(const edm::Event& iEvent);
    virtual void fillGenParticles(const edm::Event& iEvent);


    // ----------member data ---------------------------
    bool isMC;
    edm::EDGetTokenT<edm::TriggerResults> triggerToken;
    edm::EDGetTokenT<std::vector<pat::TriggerObjectStandAlone>> triggerobjectsToken;
    edm::EDGetTokenT<trigger::TriggerEvent> triggersummaryToken;
    edm::EDGetTokenT<reco::BeamSpot> beamspotToken;
    edm::EDGetTokenT<reco::VertexCollection> vertexToken;
    edm::EDGetTokenT<std::vector<pat::Muon>> muonsToken;  
    edm::EDGetTokenT<std::vector<pat::Electron>> electronsToken;
    edm::EDGetTokenT<std::vector<reco::Conversion>> conversionsToken;
    edm::EDGetTokenT<edm::View<pat::Photon>> photonToken;
    edm::EDGetTokenT<double> rhoToken;
    edm::EDGetTokenT<edm::ValueMap<float>> full5x5SigmaIEtaIEtaMapToken;
    edm::EDGetTokenT<edm::ValueMap<float>> phoChargedIsolationToken;
    edm::EDGetTokenT<edm::ValueMap<float>> phoNeutralHadronIsolationToken;
    edm::EDGetTokenT<edm::ValueMap<float>> phoPhotonIsolationToken;
    edm::FileInPath effAreaChHadronsFile;
    edm::FileInPath effAreaNeuHadronsFile;
    edm::FileInPath effAreaPhotonsFile;
    edm::EDGetTokenT<std::vector<pat::Jet>> jetToken;
    edm::EDGetTokenT<std::vector<pat::MET>> metToken;
    edm::EDGetTokenT<std::vector<reco::GenParticle>> genparticlesToken;
    edm::EDGetTokenT<GenEventInfoProduct> generatorToken;


    ESHandle<MagneticField> B;
    edm::ESHandle<TransientTrackBuilder> theTTBuilder;    
    edm::Handle<reco::BeamSpot> beamSpotHandle;
    reco::BeamSpot beamSpot;
    edm::Handle<std::vector<reco::Vertex>> pvHandle;
    reco::Vertex vtx;

    std::map<std::string,TTree*> tree_;
    NtupleEvent event_;   

};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
NtupleMaker::NtupleMaker(const edm::ParameterSet& iConfig):
    isMC                           (iConfig.getUntrackedParameter<bool>("isMC")),    
    triggerToken                   (consumes<edm::TriggerResults>                       (iConfig.getParameter<edm::InputTag>("TriggerResults"))),
    triggerobjectsToken            (consumes<std::vector<pat::TriggerObjectStandAlone>> (iConfig.getParameter<edm::InputTag>("TriggerObjects"))),
    triggersummaryToken            (consumes<trigger::TriggerEvent>                     (iConfig.getParameter<edm::InputTag>("TriggerSummary"))),
    beamspotToken                  (consumes<reco::BeamSpot>                            (iConfig.getParameter<edm::InputTag>("BeamSpot"))),
    vertexToken                    (consumes<std::vector<reco::Vertex>>                 (iConfig.getParameter<edm::InputTag>("Vertex"))),
    muonsToken                     (consumes<std::vector<pat::Muon>>                    (iConfig.getParameter<edm::InputTag>("Muons"))),
    electronsToken                 (consumes<std::vector<pat::Electron>>                (iConfig.getParameter<edm::InputTag>("Electrons"))),
    conversionsToken               (consumes<std::vector<reco::Conversion>>             (iConfig.getParameter<edm::InputTag>("Conversions"))),
    photonToken                    (consumes<edm::View<pat::Photon>>                    (iConfig.getParameter<edm::InputTag>("Photons"))),
    rhoToken                       (consumes<double>                                    (iConfig.getParameter<edm::InputTag>("rho"))),
    full5x5SigmaIEtaIEtaMapToken   (consumes<edm::ValueMap<float>>                      (iConfig.getParameter<edm::InputTag>("full5x5SigmaIEtaIEtaMap"))),
    phoChargedIsolationToken       (consumes<edm::ValueMap<float>>                      (iConfig.getParameter<edm::InputTag>("phoChargedIsolation"))),
    phoNeutralHadronIsolationToken (consumes<edm::ValueMap<float>>                      (iConfig.getParameter<edm::InputTag>("phoNeutralHadronIsolation"))),
    phoPhotonIsolationToken        (consumes<edm::ValueMap<float>>                      (iConfig.getParameter<edm::InputTag>("phoPhotonIsolation"))),
    effAreaChHadronsFile           (iConfig.getUntrackedParameter<edm::FileInPath>( "effAreaChHadFile", edm::FileInPath("RecoEgamma/PhotonIdentification/data/PHYS14/effAreaPhotons_cone03_pfChargedHadrons_V2.txt") )),
    effAreaNeuHadronsFile          (iConfig.getUntrackedParameter<edm::FileInPath>( "effAreaNeuHadFile", edm::FileInPath("RecoEgamma/PhotonIdentification/data/PHYS14/effAreaPhotons_cone03_pfNeutralHadrons_V2.txt") )),
    effAreaPhotonsFile             (iConfig.getUntrackedParameter<edm::FileInPath>( "effAreaPhoFile", edm::FileInPath("RecoEgamma/PhotonIdentification/data/PHYS14/effAreaPhotons_cone03_pfPhotons_V2.txt") )),
    jetToken                       (consumes<std::vector<pat::Jet>>                     (iConfig.getParameter<edm::InputTag>("Jets"))),
    metToken                       (consumes<std::vector<pat::MET>>                     (iConfig.getParameter<edm::InputTag>("MET"))),
    genparticlesToken              (consumes<std::vector<reco::GenParticle>>            (iConfig.getParameter<edm::InputTag>("GenParticles"))),
    generatorToken                 (consumes<GenEventInfoProduct>                       (iConfig.getParameter<edm::InputTag>("Generator")))
{}


NtupleMaker::~NtupleMaker() {

   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)
}


//
// member functions
//

// ------------ method called for each event  ------------
void NtupleMaker::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {

    event_.triggers.clear();
    event_.triggerobjects.clear();
    event_.muons.clear();
    event_.electrons.clear();
    event_.dimuons.clear();
    event_.dielectrons.clear();
    event_.emus.clear();
    event_.photons.clear();
    event_.jets.clear();
    event_.MET.clear();
    event_.genparticles.clear();

    iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",theTTBuilder);
    iSetup.get<IdealMagneticFieldRecord>().get(B);

    iEvent.getByToken(beamspotToken, beamSpotHandle);
    beamSpot = (*beamSpotHandle);

    iEvent.getByToken(vertexToken, pvHandle);
    vtx = pvHandle->front();

    event_.run       = iEvent.id().run();
    event_.event     = iEvent.id().event();
    event_.lumi      = iEvent.id().luminosityBlock();
    event_.nVertices = pvHandle->size();

    fillTriggers(iEvent);
    fillMuons(iEvent);
    fillElectrons(iEvent);
    fillDileptons(iEvent);
    fillPhotons(iEvent);
    fillJets(iEvent);
    fillMETs(iEvent);
    if(isMC) fillGenParticles(iEvent);
    tree_["physicsTree"]->Fill();

}//


// ------------ method called once each job just before starting event loop  ------------
void NtupleMaker::beginJob() {

    TH1::SetDefaultSumw2() ;
    edm::Service<TFileService> fs;
    tree_["physicsTree"] = fs->make<TTree>("physicsTree","physicsTree");
    tree_["physicsTree"]->Branch("event" ,&event_);

}//

// ------------ method called once each job just after ending the event loop  ------------
void NtupleMaker::endJob() {

}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void NtupleMaker::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
   //The following says we do not know what parameters are allowed so do no validation
   // Please change this to state exactly what you do use, even if it is no parameters
   edm::ParameterSetDescription desc;
   desc.setUnknown();
   descriptions.addDefault(desc);
}

void NtupleMaker::fillTriggers(const edm::Event &iEvent) {

    NtupleTrigger trigger_;
    NtupleTriggerObject triggerobject_;

    edm::Handle< TriggerResults > ResultHandle;
    iEvent.getByToken(triggerToken,ResultHandle);

    edm::Handle<std::vector<pat::TriggerObjectStandAlone>> trigCol;
    iEvent.getByToken(triggerobjectsToken, trigCol);

    string trigs[] = {
        "HLT_IsoMu17_eta2p1_v*",
        "HLT_IsoMu20_v*",
        "HLT_IsoTkMu20_v*", 
        "HLT_IsoTkMu20_eta2p1_v*", 
        "HLT_Mu50_v*",
        "HLT_Mu45_eta2p1_v*",       
        "HLT_Mu17_TrkIsoVVL_v*",
        "HLT_DoubleIsoMu17_eta2p1_v*",                
        "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v*",
        "HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v*",
        "HLT_Mu27_TkMu8_v*",
        "HLT_Ele22_eta2p1_WP75_Gsf_v",
        "HLT_Ele22_eta2p1_WPLoose_Gsf_v",
        "HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ"   
    };

    cout<<"HLT"<<endl;
    const unsigned ntrigs = sizeof(trigs)/sizeof(*trigs);
    const unsigned nTrigg = ResultHandle.product()->size();
    std::vector<std::pair<std::string, int>> indices;
    if(nTrigg==0) cout<<"No trigger result"<<endl;
    edm::TriggerNames trigList = iEvent.triggerNames(*ResultHandle);

    for( unsigned i = 0; i != nTrigg; i++ ) {
        string _trigName = trigList.triggerName(i);
        for( unsigned j = 0; j != ntrigs; j++ ) {
            if( _trigName.find(trigs[j].substr(0,trigs[j].find("*"))) != std::string::npos ) {
                trigger_.name = trigs[j];
                if( ResultHandle.product()->accept(i) ) {
                    cout<<"Fired: "<<trigs[j]<<endl;
                    trigger_.isFired = true;
                }
                else {
                    trigger_.isFired = false;
                }
                event_.triggers.push_back(trigger_);
            }
        }
    }

    for(unsigned i=0; i<trigCol->size(); i++) {

        pat::TriggerObjectStandAlone trigObj = trigCol->at(i);
        trigObj.unpackPathNames(trigList);
        const std::vector< std::string > pathNames = trigObj.pathNames();

        for(std::vector<std::string>::const_iterator name = pathNames.begin(); name != pathNames.end(); ++name){
            for(unsigned int j=0; j<ntrigs; j++) {
                if( name->find(trigs[j].substr(0,trigs[j].find("*"))) != std::string::npos && trigObj.hasPathName(*name,true,true) ) {
                    std::cout<<trigs[j]<<std::endl;
                    triggerobject_.eta = trigObj.eta();
                    triggerobject_.phi = trigObj.phi();
                    triggerobject_.name = trigs[j];
                    event_.triggerobjects.push_back(triggerobject_);
                }
            }
        }
    }

}//

// ------------ method fills ntuple with the muon information
void NtupleMaker::fillMuons(const edm::Event& iEvent) {

    NtupleMuon mu_;

    edm::Handle<std::vector<pat::Muon>> muons;
    iEvent.getByToken(muonsToken,muons);

    event_.nMuons = muons->size();
    if( muons->size() == 0 ) return;

    for( std::vector<pat::Muon>::const_iterator mu = muons->begin(); mu != muons->end(); ++mu ) {    

        mu_.isStandAloneMuon = mu->isStandAloneMuon();
        mu_.isGlobalMuon     = mu->isGlobalMuon();
        mu_.isTrackerMuon    = mu->isTrackerMuon();
        mu_.isPFMuon         = mu->isPFMuon();
        mu_.px               = mu->px();
        mu_.py               = mu->py();
        mu_.pz               = mu->pz();
        mu_.pt               = mu->pt();
        mu_.eta              = mu->eta();
        mu_.phi              = mu->phi();
        mu_.charge           = mu->charge();
        mu_.nChambers        = mu->numberOfChambers();
        mu_.stationMask      = mu->stationMask();
        mu_.nMatchedStations = mu->numberOfMatchedStations();

        if( mu->isGlobalMuon() ) { // Global Muon

            reco::TrackRef glbTrack  =  mu->globalTrack();

            if( glbTrack.isNonnull() ) {

                const reco::HitPattern & glbhit = glbTrack->hitPattern();

                mu_.normalizedChi2  =  glbTrack->normalizedChi2();
                mu_.nValidHits      =  glbTrack->numberOfValidHits();         
                mu_.nValidMuonHits  = glbhit.numberOfValidMuonHits();
                mu_.qoverp          = glbTrack->qoverp();
                mu_.theta           = glbTrack->theta();
                mu_.lambda          = glbTrack->lambda();
                mu_.dxy             = glbTrack->dxy();
                mu_.d0              = glbTrack->d0();
                mu_.dsz             = glbTrack->dsz();
                mu_.dz              = glbTrack->dz();
                mu_.dxyBS           = glbTrack->dxy(beamSpot.position());
                mu_.dszBS           = glbTrack->dsz(beamSpot.position());
                mu_.dzBS            = glbTrack->dz(beamSpot.position());
                mu_.vx              = glbTrack->vx();
                mu_.vy              = glbTrack->vy();
                mu_.vz              = glbTrack->vz();

            }

            reco::TrackRef trackerTrack = mu->innerTrack();

            if( trackerTrack.isNonnull() ) {

                const reco::HitPattern & inhit = trackerTrack->hitPattern();

                mu_.nValidTrackerHits = inhit.numberOfValidTrackerHits();
                mu_.nValidPixelHits   = inhit.numberOfValidPixelHits();
                mu_.nTrackerLayers    = inhit.trackerLayersWithMeasurement();

            }

        } // Global Muon 
        else if( mu->isStandAloneMuon() ) { // STA Muon

            reco::TrackRef muonTrack  =  mu->outerTrack();

            if( muonTrack.isNonnull() ) {

                const reco::HitPattern & muonhit  =  muonTrack->hitPattern();

                mu_.nValidMuonHits  =  muonhit.numberOfValidMuonHits();

            }

        } // STA Muon
        else if(mu->isTrackerMuon() ) { // Tracker Muon

            reco::TrackRef trackerTrack = mu->innerTrack();

            if( trackerTrack.isNonnull() ) {

                const reco::HitPattern & inhit = trackerTrack->hitPattern();

                mu_.normalizedChi2    = trackerTrack->normalizedChi2();
                mu_.nValidHits        = trackerTrack->numberOfValidHits();                    
                mu_.nValidTrackerHits = inhit.numberOfValidTrackerHits();
                mu_.nValidPixelHits   = inhit.numberOfValidPixelHits();
                mu_.nTrackerLayers    = inhit.trackerLayersWithMeasurement();
                mu_.qoverp            = trackerTrack->qoverp();
                mu_.theta             = trackerTrack->theta();
                mu_.lambda            = trackerTrack->lambda();
                mu_.dxy               = trackerTrack->dxy();
                mu_.d0                = trackerTrack->d0();
                mu_.dsz               = trackerTrack->dsz();
                mu_.dz                = trackerTrack->dz();
                mu_.dxyBS             = trackerTrack->dxy(beamSpot.position());
                mu_.dszBS             = trackerTrack->dsz(beamSpot.position());
                mu_.dzBS              = trackerTrack->dz(beamSpot.position());
                mu_.vx                = trackerTrack->vx();
                mu_.vy                = trackerTrack->vy();
                mu_.vz                = trackerTrack->vz();

            }
        } // Tracker Muon

        // MuonBestTrack
        if( mu->muonBestTrack().isNonnull() ) {

            mu_.muonBestTrack_nTrackerLayers = mu->muonBestTrack()->hitPattern().trackerLayersWithMeasurement();
            mu_.muonBestTrack_px             = mu->muonBestTrack()->px();
            mu_.muonBestTrack_py             = mu->muonBestTrack()->py();
            mu_.muonBestTrack_pz             = mu->muonBestTrack()->pz();
            mu_.muonBestTrack_pt             = mu->muonBestTrack()->pt();
            mu_.muonBestTrack_ptError        = mu->muonBestTrack()->ptError();

            if( !pvHandle->empty() && !pvHandle->front().isFake() ) {

                mu_.dxyVTX = mu->muonBestTrack()->dxy(vtx.position());
                mu_.dszVTX = mu->muonBestTrack()->dsz(vtx.position());
                mu_.dzVTX  = mu->muonBestTrack()->dz(vtx.position());
            }

        }

        // InnerTrack
        if( mu->innerTrack().isNonnull() ) {

            mu_.innerTrack_nTrackerLayers = mu->innerTrack()->hitPattern().trackerLayersWithMeasurement();
            mu_.innerTrack_px             = mu->innerTrack()->px();
            mu_.innerTrack_py             = mu->innerTrack()->py();
            mu_.innerTrack_pz             = mu->innerTrack()->pz();
            mu_.innerTrack_pt             = mu->innerTrack()->pt();
            mu_.innerTrack_ptError        = mu->innerTrack()->ptError();

            if( !pvHandle->empty() && !pvHandle->front().isFake() ) {

                mu_.innerTrack_dxyVTX = mu->innerTrack()->dxy(vtx.position());
                mu_.innerTrack_dszVTX = mu->innerTrack()->dsz(vtx.position());
                mu_.innerTrack_dzVTX  = mu->innerTrack()->dz(vtx.position());

            }

        }

        // TunePMuonBestTrack
        if( mu->tunePMuonBestTrack().isNonnull() ) {

            mu_.tunePMuonBestTrack_nTrackerLayers = mu->tunePMuonBestTrack()->hitPattern().trackerLayersWithMeasurement();
            mu_.tunePMuonBestTrack_px             = mu->tunePMuonBestTrack()->px();
            mu_.tunePMuonBestTrack_py             = mu->tunePMuonBestTrack()->py();
            mu_.tunePMuonBestTrack_pz             = mu->tunePMuonBestTrack()->pz();
            mu_.tunePMuonBestTrack_pt             = mu->tunePMuonBestTrack()->pt();
            mu_.tunePMuonBestTrack_ptError        = mu->tunePMuonBestTrack()->ptError();

            if( !pvHandle->empty() && !pvHandle->front().isFake() ) {

                mu_.tunePMuonBestTrack_dxyVTX = mu->tunePMuonBestTrack()->dxy(vtx.position());
                mu_.tunePMuonBestTrack_dszVTX = mu->tunePMuonBestTrack()->dsz(vtx.position());
                mu_.tunePMuonBestTrack_dzVTX  = mu->tunePMuonBestTrack()->dz(vtx.position());

            }

        }   

        // Isolation
        mu_.isolationR03_sumpt    = mu->isolationR03().sumPt;
        mu_.isolationR03_hadEt    = mu->isolationR03().hadEt;
        mu_.isolationR03_emEt     = mu->isolationR03().emEt;
        mu_.isolationR05_sumpt    = mu->isolationR05().sumPt;
        mu_.isolationR05_hadEt    = mu->isolationR05().hadEt;
        mu_.isolationR05_emEt     = mu->isolationR05().emEt; 
        mu_.PfChargedHadronIsoR04 = mu->pfIsolationR04().sumChargedHadronPt;
        mu_.PfNeutralHadronIsoR04 = mu->pfIsolationR04().sumNeutralHadronEt;
        mu_.PfGammaIsoR04         = mu->pfIsolationR04().sumPhotonEt;
        mu_.PfChargedHadronIsoR03 = mu->pfIsolationR03().sumChargedHadronPt;
        mu_.PfNeutralHadronIsoR03 = mu->pfIsolationR03().sumNeutralHadronEt;
        mu_.PfGammaIsoR03         = mu->pfIsolationR03().sumPhotonEt;

        event_.muons.push_back(mu_);
    }

}//

void NtupleMaker::fillElectrons(const edm::Event& iEvent) {

    NtupleElectron el_;

    edm::Handle< double > rhoH;
    iEvent.getByToken(rhoToken,rhoH);
    event_.rho = *rhoH;

    edm::Handle<std::vector<reco::Conversion>> conversions;
    iEvent.getByToken(conversionsToken, conversions);

    edm::Handle<std::vector<pat::Electron>> electrons;
    iEvent.getByToken(electronsToken, electrons);

    event_.nElectrons = electrons->size();
    if( electrons->size() == 0 ) return;

    for( std::vector<pat::Electron>::const_iterator el = electrons->begin(); el != electrons->end(); ++el ) {

        el_.pt     = el->pt();
        el_.eta    = el->eta();
        el_.rap    = el->rapidity();
        el_.phi    = el->phi();
        el_.E      = el->energy();
        el_.charge = el->charge();

        double R = sqrt(el->superCluster()->x()*el->superCluster()->x() + el->superCluster()->y()*el->superCluster()->y() +el->superCluster()->z()*el->superCluster()->z());
        double Rt = sqrt(el->superCluster()->x()*el->superCluster()->x() + el->superCluster()->y()*el->superCluster()->y());

        el_.enSC = el->superCluster()->energy();
        el_.preEnSC  = el->superCluster()->preshowerEnergy();
        el_.rawEnSC = el->superCluster()->rawEnergy();
        el_.etSC =  (el->superCluster()->energy())*(Rt/R);
        el_.etaSC = el->superCluster()->eta();
        el_.phiSC = el->superCluster()->phi();

        // ECAL
        el_.sigmaIetaIeta = el->full5x5_sigmaIetaIeta();
        el_.E1x5          = el->e1x5();
        el_.E2x5          = el->e2x5Max();
        el_.E5x5          = el->e5x5();
        el_.hOverE        = el->hcalOverEcal();
        el_.etaScWidth    = el->superCluster()->etaWidth();
        el_.phiScWidth    = el->superCluster()->phiWidth();
        el_.r9            = el->r9();

        // ECAL + Track
        el_.dEtaIn = el->deltaEtaSuperClusterTrackAtVtx();
        el_.dPhiIn = el->deltaPhiSuperClusterTrackAtVtx();
        // |1/E-1/p| = |1/E - EoverPinner/E| is computed below. The if protects against ecalEnergy == inf or zero
        // (always the case for miniAOD for electrons <5 GeV)
        if( el->ecalEnergy() == 0 ){
          //printf("Electron energy is zero!\n");
          el_.ooEmooP = 1e30;
        }else if( !std::isfinite(el->ecalEnergy())){
          //printf("Electron energy is not finite!\n");
          el_.ooEmooP = 1e30;
        }else{
          el_.ooEmooP = fabs(1.0/el->ecalEnergy() - el->eSuperClusterOverP()/el->ecalEnergy() );
        }

        // Isolation
        GsfElectron::PflowIsolationVariables pfIso = el->pfIsolationVariables();

        // Compute isolation with delta beta correction for PU
        el_.isoChargedHadrons = pfIso.sumChargedHadronPt;
        el_.isoNeutralHadrons = pfIso.sumNeutralHadronEt;
        el_.isoPhotons        = pfIso.sumPhotonEt;
        el_.isoChargedFromPU  = pfIso.sumPUPt;

        float abseta = fabs(el->superCluster()->eta());

        // The effective areas constants file in the local release or default CMSSW, whichever is found
        edm::FileInPath eaConstantsFile("RecoEgamma/ElectronIdentification/data/PHYS14/effAreaElectrons_cone03_pfNeuHadronsAndPhotons.txt");
        EffectiveAreas effectiveAreas(eaConstantsFile.fullPath());
        float eA = effectiveAreas.getEffectiveArea(abseta);

        el_.isoDeltaBeta = (pfIso.sumChargedHadronPt + max<float>( 0.0, pfIso.sumNeutralHadronEt + pfIso.sumPhotonEt - 0.5 * pfIso.sumPUPt))/(el->pt());
        el_.isoRho       = (pfIso.sumChargedHadronPt + max<float>( 0.0, pfIso.sumNeutralHadronEt + pfIso.sumPhotonEt - event_.rho * eA))/(el->pt());

        // Track - Impact Parameter, Conversion rejection, Converted
        //cout<<"gsfTrack1"<<endl;
        reco::GsfTrackRef theTrack = el->gsfTrack();
        el_.d0 = (-1) * theTrack->dxy(vtx.position());
        el_.dz = theTrack->dz(vtx.position());

        el_.expectedMissingInnerHits =el->gsfTrack()->hitPattern().numberOfHits(reco::HitPattern::MISSING_INNER_HITS);

        bool passConvVeto = !ConversionTools::hasMatchedConversion(*el, conversions, beamSpotHandle->position());
        el_.passConversionVeto = passConvVeto;
        el_.brem =el->fbrem();
        //cout<<"gsfTrack1"<<endl;
        //mvaValue = (*mvaValues)[el] );
        //mvaCategory = (*mvaCategories)[el] );

        el_.eleInBarrel =el->isEB();
        el_.eleInEndcap =el->isEE();

        // ECAL driven
        el_.eleEcalDrivenSeed =el->ecalDrivenSeed();
        //cout<<"ECAL driven: "<<el->ecalDrivenSeed()<<endl;

        event_.electrons.push_back(el_);
    }

}//

void NtupleMaker::fillDileptons(const edm::Event& iEvent) {

    NtupleDimuon dimu_;
    NtupleDielectron diel_;
    NtupleEmu emu_;

    edm::Handle<std::vector<pat::Muon>> muons;
    iEvent.getByToken(muonsToken,muons);

    edm::Handle<std::vector<pat::Electron>> electrons;
    iEvent.getByToken(electronsToken, electrons);

    for( unsigned i = 0; i != muons->size(); i++ ) {
        const pat::Muon mu = muons->at(i);
        for( unsigned j = 0; j != muons->size(); j++ ) {
            const pat::Muon mu2 = muons->at(j);
            if( i >= j ) continue;

            if((mu.isGlobalMuon()||mu.isTrackerMuon())&&(mu2.isGlobalMuon()||mu2.isTrackerMuon())) {

                reco::TrackRef tkTrack = mu.innerTrack();
                reco::TrackRef tkTrack2 = mu2.innerTrack();
                reco::TransientTrack muTransient1(tkTrack, B.product());
                reco::TransientTrack muTransient2(tkTrack2, B.product());

                vector<reco::TransientTrack> DimuonTracksTrk;
                DimuonTracksTrk.push_back(muTransient1);
                DimuonTracksTrk.push_back(muTransient2);
                KalmanVertexFitter KalmanFitterTrk(true);
                CachingVertex<5> vertexTrk;
                TransientVertex vtxtmptrk;
                bool isVertexTrk = true;
                try {
                    vertexTrk = KalmanFitterTrk.vertex(DimuonTracksTrk);
                    vtxtmptrk = KalmanFitterTrk.vertex(DimuonTracksTrk);
                }
                catch( exception & err ) {
                    isVertexTrk = false;
                }
                if( isVertexTrk && vertexTrk.isValid() ) {
                // inv. mass refit using the Dimuon vtx
                    InvariantMassFromVertex imfvTrk;
                    static const double muon_mass = 0.1056583;
                    const CachingVertex<5>& vtxTrk = vertexTrk;
                    //InvariantMassFromVertex::LorentzVector new_p4Trk = imfvTrk.p4(vtxTrk, muon_mass);
                    Measurement1D new_massTrk = imfvTrk.invariantMass(vtxTrk, muon_mass);

                    dimu_.vertexFitChi2     = vtxTrk.totalChiSquared();
                    dimu_.vertexFitNdof     = vtxTrk.degreesOfFreedom();
                    dimu_.vertexFitChi2Ndof = vtxTrk.totalChiSquared()/vtxTrk.degreesOfFreedom();
                    dimu_.vertexFitProb     = TMath::Prob(vtxTrk.totalChiSquared(),(int)vtxTrk.degreesOfFreedom());
                    dimu_.X = i;
                    dimu_.Y = j;
                }

                // cosmic variable
                double cosine = acos(-tkTrack->momentum().Dot(tkTrack2->momentum()/tkTrack->p()/tkTrack2->p()));
                dimu_.openingAngle = cosine;
                //CosAngle.push_back(cosine);
                event_.dimuons.push_back(dimu_);
            }
        }
    }
    for( unsigned int i = 0; i != electrons->size(); i++ ) {
        const pat::Electron el = electrons->at(i);
        for( unsigned int j = 0; j != electrons->size(); j++ ) {
            const pat::Electron el2 = electrons->at(j);
            if( i >= j ) continue;

            //cout<<"gsfTrack2"<<endl;
            reco::GsfTrackRef elecTrk = el.gsfTrack();
            reco::GsfTrackRef elecTrk2 = el2.gsfTrack();    
            //cout<<"gsfTrack2"<<endl;

            if( elecTrk.isNonnull() && elecTrk2.isNonnull() ) {
                //reco::GsfTransientTrack eTransient1(elecTrk, B.product());
                //reco::GsfTransientTrack eTransient2(elecTrk2, B.product());

                vector<reco::TransientTrack> dielecTracksTrk;
                dielecTracksTrk.push_back(theTTBuilder->build(elecTrk));
                dielecTracksTrk.push_back(theTTBuilder->build(elecTrk2));
                KalmanVertexFitter KalmanFitterTrk(true);
                CachingVertex<5> vertexTrk;
                TransientVertex vtxtmptrk;
                bool isVertexTrk = true;
                try {
                    vertexTrk = KalmanFitterTrk.vertex(dielecTracksTrk);
                    vtxtmptrk = KalmanFitterTrk.vertex(dielecTracksTrk);
                }
                catch( exception & err ) {
                    isVertexTrk = false;
                }
                if( isVertexTrk && vertexTrk.isValid() ) {
                    // inv. mass refit using the dielec vtx
                    InvariantMassFromVertex imfvTrk;
                    static const double elec_mass = 0.000511;
                    const CachingVertex<5>& vtxTrk = vertexTrk;
                    //InvariantMassFromVertex::LorentzVector new_p4Trk = imfvTrk.p4(vtxTrk, elec_mass);
                    Measurement1D new_massTrk = imfvTrk.invariantMass(vtxTrk, elec_mass);

                    diel_.vertexFitChi2     = vtxTrk.totalChiSquared();
                    diel_.vertexFitNdof     = vtxTrk.degreesOfFreedom();
                    diel_.vertexFitChi2Ndof = vtxTrk.totalChiSquared()/vtxTrk.degreesOfFreedom();
                    diel_.vertexFitProb     = TMath::Prob(vtxTrk.totalChiSquared(),(int)vtxTrk.degreesOfFreedom());
                    diel_.X = i;
                    diel_.Y = j;
                }

                double cosine = acos(-elecTrk->momentum().Dot(elecTrk2->momentum()/elecTrk->p()/elecTrk2->p()));
                diel_.openingAngle = cosine;

                event_.dielectrons.push_back(diel_);
            }
        }
    }
    for( unsigned i = 0; i != muons->size(); i++ ) {
        const pat::Muon mu = muons->at(i);
        if( !mu.isTrackerMuon() || !mu.isGlobalMuon() ) continue;
        reco::TrackRef tkTrack = mu.innerTrack(); 

        for( unsigned j = 0; j != electrons->size(); j++ ) {
            const pat::Electron ielec2 = electrons->at(j);
            //cout<<"gsfTrack3"<<endl;
            reco::GsfTrackRef elecTrk2 = ielec2.gsfTrack();
            //cout<<"gsfTrack3"<<endl;

            if( tkTrack.isNonnull() && elecTrk2.isNonnull() ) {
                reco::TransientTrack muTransient1(tkTrack, B.product());

                vector<reco::TransientTrack> emuTracksTrk;
                emuTracksTrk.push_back(muTransient1);
                emuTracksTrk.push_back(theTTBuilder->build(elecTrk2));

                KalmanVertexFitter KalmanFitterTrk(true);
                CachingVertex<5> vertexTrk;
                TransientVertex vtxtmptrk;
                bool isVertexTrk = true;
                try {
                    vertexTrk = KalmanFitterTrk.vertex(emuTracksTrk);
                    vtxtmptrk = KalmanFitterTrk.vertex(emuTracksTrk);
                }
                catch( exception & err ) {
                    isVertexTrk = false;
                }
                if( isVertexTrk && vertexTrk.isValid() ) {
                    // inv. mass refit using the emu vtx
                    const CachingVertex<5>& vtxTrk = vertexTrk;
                    emu_.vertexFitChi2     = vtxTrk.totalChiSquared();
                    emu_.vertexFitNdof     = vtxTrk.degreesOfFreedom();
                    emu_.vertexFitChi2Ndof = vtxTrk.totalChiSquared()/vtxTrk.degreesOfFreedom();
                    emu_.vertexFitProb     = TMath::Prob(vtxTrk.totalChiSquared(),(int)vtxTrk.degreesOfFreedom());
                    emu_.X = i;
                    emu_.Y = j;
                }

                double cosine = acos(-tkTrack->momentum().Dot(elecTrk2->momentum()/tkTrack->p()/elecTrk2->p()));
                emu_.openingAngle = cosine;

                event_.emus.push_back(emu_);
            }
        }
    }

}//

void NtupleMaker::fillPhotons(const edm::Event& iEvent) {

    NtuplePhoton pho_;

    edm::Handle<edm::View<pat::Photon>> photons;
    iEvent.getByToken(photonToken, photons);

    edm::Handle< double > rhoHH;
    iEvent.getByToken(rhoToken,rhoHH);
    float rho_ = *rhoHH;

    edm::Handle<edm::ValueMap<float> > full5x5SigmaIEtaIEtaMap;
    iEvent.getByToken(full5x5SigmaIEtaIEtaMapToken, full5x5SigmaIEtaIEtaMap);

    edm::Handle<edm::ValueMap<float> > phoChargedIsolationMap;
    iEvent.getByToken(phoChargedIsolationToken, phoChargedIsolationMap);

    edm::Handle<edm::ValueMap<float> > phoNeutralHadronIsolationMap;
    iEvent.getByToken(phoNeutralHadronIsolationToken, phoNeutralHadronIsolationMap);

    edm::Handle<edm::ValueMap<float> > phoPhotonIsolationMap;
    iEvent.getByToken(phoPhotonIsolationToken, phoPhotonIsolationMap);

    EffectiveAreas effAreaChHadrons_( effAreaChHadronsFile.fullPath() );
    EffectiveAreas effAreaNeuHadrons_( effAreaNeuHadronsFile.fullPath() );
    EffectiveAreas effAreaPhotons_( effAreaPhotonsFile.fullPath() );

    event_.nPhotons = photons->size();
    if( event_.nPhotons == 0 ) return;

    for(unsigned i=0; i< photons->size(); ++i) {
        const auto pho = photons->ptrAt(i);

        pho_.pt                    = pho->pt();
        pho_.eta                   = pho->eta();
        pho_.phi                   = pho->phi();
        pho_.etaSC                 = pho->superCluster()->eta();
        pho_.phiSC                 = pho->superCluster()->phi();
        pho_.HoverE                = pho->hadTowOverEm();
        pho_.hasPixelSeed          = (Int_t)pho->hasPixelSeed();
        pho_.Full5x5_SigmaIEtaIEta = (*full5x5SigmaIEtaIEtaMap)[ pho ];

        float chIso = (double)(*phoChargedIsolationMap)[pho];
        float nhIso = (double)(*phoNeutralHadronIsolationMap)[pho];
        float phIso = (double)(*phoPhotonIsolationMap)[pho];
        pho_.ChIso = chIso;
        pho_.NhIso = nhIso;
        pho_.PhIso = phIso;

        float abseta = fabs( pho->superCluster()->eta());
        pho_.ChIsoWithEA = std::max( (float)0.0, chIso - rho_*effAreaChHadrons_.getEffectiveArea(abseta) );
        pho_.NhIsoWithEA = std::max( (float)0.0, nhIso - rho_*effAreaNeuHadrons_.getEffectiveArea(abseta) );
        pho_.PhIsoWithEA = std::max( (float)0.0, phIso - rho_*effAreaPhotons_.getEffectiveArea(abseta) );

        event_.photons.push_back(pho_);
    }

}//

void NtupleMaker::fillJets(const edm::Event& iEvent) {

    NtupleJet jet_;

    edm::Handle<std::vector<pat::Jet>> jets;
    iEvent.getByToken(jetToken,jets);

    event_.nJets = jets->size();
    if( event_.nJets == 0 ) return;

    for( vector<pat::Jet>::const_iterator jet = jets->begin(); jet != jets->end(); ++jet ) {
        jet_.pt       = jet->pt();
        jet_.eta      = jet->eta();
        jet_.phi      = jet->phi();
        jet_.charge   = jet->jetCharge();
        jet_.flavor   = jet->partonFlavour();
        jet_.bTag     = jet->bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags");
        jet_.CHfrac   = jet->chargedHadronEnergyFraction();
        jet_.NHfrac   = jet->neutralHadronEnergyFraction();
        jet_.NHEMfrac = jet->neutralEmEnergyFraction();
        jet_.CHEMfrac = jet->chargedEmEnergyFraction();
        jet_.CHmulti  = jet->chargedMultiplicity();
        jet_.NHmulti  = jet->neutralMultiplicity();

        event_.jets.push_back(jet_);
    }

}//

void NtupleMaker::fillMETs(const edm::Event& iEvent) {

    NtupleMET MET_;

    edm::Handle<std::vector<pat::MET>> MET;
    iEvent.getByToken(metToken,MET);

    MET_.pt          = MET->front().uncorPt();
    MET_.phi         = MET->front().uncorPhi();
    MET_.px          = MET->front().uncorPx();
    MET_.py          = MET->front().uncorPy();
    MET_.sumEt       = MET->front().uncorSumEt();
    MET_.Type1_pt    = MET->front().pt();
    MET_.Type1_phi   = MET->front().phi();
    MET_.Type1_px    = MET->front().px();
    MET_.Type1_py    = MET->front().py();
    MET_.Type1_sumEt = MET->front().sumEt();

    event_.MET.push_back(MET_);

}//

void NtupleMaker::fillGenParticles(const edm::Event& iEvent) {

    NtupleGenParticle particle_;

    edm::Handle<std::vector<reco::GenParticle>> particles;
    iEvent.getByToken(genparticlesToken, particles);

    edm::Handle<GenEventInfoProduct> genInfo;
    iEvent.getByToken(generatorToken, genInfo);

    if( genInfo->weight() > 0 ) event_.weight = 1.0;
    else event_.weight = -1.0;

    event_.nGenParticles = particles->size();
    if( particles->size() == 0 ) return;

    int _GennParticles = 0;
    for( std::vector<reco::GenParticle>::const_iterator particle = particles->begin(); particle != particles->end(); ++particle ) {

        if( abs(particle->pdgId()) == 13 || abs(particle->pdgId()) == 11 || abs(particle->pdgId()) == 15 ) {

            particle_.id     = particle->pdgId(); 
            particle_.pt     = particle->pt(); 
            particle_.px     = particle->px();
            particle_.py     = particle->py();
            particle_.pz     = particle->pz();
            particle_.eta    = particle->eta();
            particle_.phi    = particle->phi();
            particle_.charge = particle->charge();
            particle_.status = particle->status();
            particle_.mother = particle->mother(0)->pdgId();

            particle_.isHardProcess                                = particle->isHardProcess();
            particle_.isPromptFinalState                           = particle->isPromptFinalState();
            particle_.fromHardProcessFinalState                    = particle->fromHardProcessFinalState();
            particle_.fromHardProcessBeforeFSR                     = particle->fromHardProcessBeforeFSR();
            particle_.isDirectHardProcessTauDecayProductFinalState = particle->isDirectHardProcessTauDecayProductFinalState();
            particle_.isDirectPromptTauDecayProductFinalState      = particle->isDirectPromptTauDecayProductFinalState();

            _GennParticles++;
            event_.genparticles.push_back(particle_);
        }
    }
    event_.nGenParticles = _GennParticles;

}//

//define this as a plug-in
DEFINE_FWK_MODULE(NtupleMaker);
