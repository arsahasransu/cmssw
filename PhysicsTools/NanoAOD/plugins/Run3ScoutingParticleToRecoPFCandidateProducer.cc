// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/METReco/interface/PFMET.h"
#include "DataFormats/METReco/interface/PFMETCollection.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/BTauReco/interface/JetTag.h"
#include "DataFormats/Scouting/interface/Run3ScoutingParticle.h"
#include "DataFormats/Candidate/interface/LeafCandidate.h"
#include "DataFormats/Common/interface/OrphanHandle.h"

#include "DataFormats/Math/interface/deltaR.h"

#include "DataFormats/Math/interface/libminifloat.h"

#include "SimGeneral/HepPDTRecord/interface/ParticleDataTable.h"
#include "fastjet/contrib/SoftKiller.hh"

class Run3ScoutingParticleToRecoPFCandidateProducer : public edm::stream::EDProducer<> {
public:
  explicit Run3ScoutingParticleToRecoPFCandidateProducer(const edm::ParameterSet &);
  ~Run3ScoutingParticleToRecoPFCandidateProducer() override;

  static void fillDescriptions(edm::ConfigurationDescriptions &descriptions);
  void beginStream(edm::StreamID) override {}
  void produce(edm::Event & iEvent, edm::EventSetup const & setup) override;
  void endStream() override {}

  void createPFCandidates(edm::Handle<std::vector<Run3ScoutingParticle>> scoutingparticleHandle, std::unique_ptr<reco::PFCandidateCollection> & pfcands);
  void createPFCandidatesSK(edm::Handle<std::vector<Run3ScoutingParticle>> scoutingparticleHandle, std::unique_ptr<reco::PFCandidateCollection> & pfcands);
  reco::PFCandidate createPFCand(Run3ScoutingParticle scoutingparticle);
  void clearVars();

private:
  const edm::EDGetTokenT<std::vector<Run3ScoutingParticle>> input_scoutingparticle_token_;
  const edm::ESGetToken<HepPDT::ParticleDataTable, edm::DefaultRecord> particletable_token_;
  bool use_softKiller_;
  bool use_CHS_;
  const HepPDT::ParticleDataTable* pdTable_;

  std::vector<int8_t> vertexIndex;
  std::vector<float> normchi2;
  std::vector<float> dz;
  std::vector<float> dxy;
  std::vector<float> dzsig;
  std::vector<float> dxysig;
  std::vector<int> lostInnerHits;
  std::vector<int> quality;
  std::vector<float> trkPt;
  std::vector<float> trkEta;
  std::vector<float> trkPhi;
};

//
// constructors and destructor
//
Run3ScoutingParticleToRecoPFCandidateProducer::Run3ScoutingParticleToRecoPFCandidateProducer(const edm::ParameterSet &iConfig)
    : input_scoutingparticle_token_(consumes(iConfig.getParameter<edm::InputTag>("scoutingparticle"))),
      particletable_token_(esConsumes<HepPDT::ParticleDataTable, edm::DefaultRecord>()),
      use_softKiller_(iConfig.getParameter<bool>("softKiller")),
      use_CHS_(iConfig.getParameter<bool>("CHS")) {

  //register products
  produces<reco::PFCandidateCollection>();
  produces<edm::ValueMap<int>>("vertexIndex");
  produces<edm::ValueMap<float>>("normchi2");
  produces<edm::ValueMap<float>>("dz");
  produces<edm::ValueMap<float>>("dxy");
  produces<edm::ValueMap<float>>("dzsig");
  produces<edm::ValueMap<float>>("dxysig");
  produces<edm::ValueMap<int>>("lostInnerHits");
  produces<edm::ValueMap<int>>("quality");
  produces<edm::ValueMap<float>>("trkPt");
  produces<edm::ValueMap<float>>("trkEta");
  produces<edm::ValueMap<float>>("trkPhi");
}

Run3ScoutingParticleToRecoPFCandidateProducer::~Run3ScoutingParticleToRecoPFCandidateProducer() = default;

reco::PFCandidate Run3ScoutingParticleToRecoPFCandidateProducer::createPFCand(Run3ScoutingParticle scoutingparticle) {
  auto m = pdTable_->particle(HepPDT::ParticleID(scoutingparticle.pdgId())) != nullptr
                    ? pdTable_->particle(HepPDT::ParticleID(scoutingparticle.pdgId()))->mass()
                    : -99.f;
  auto q = pdTable_->particle(HepPDT::ParticleID(scoutingparticle.pdgId())) != nullptr
                    ? pdTable_->particle(HepPDT::ParticleID(scoutingparticle.pdgId()))->charge()
                    : -99.f;
  //if (m < -90 ||  q < -90) continue;
  
  float px = scoutingparticle.pt() * cos(scoutingparticle.phi());
  float py = scoutingparticle.pt() * sin(scoutingparticle.phi());
  float pz = scoutingparticle.pt() * sinh(scoutingparticle.eta());
  float p = scoutingparticle.pt() * cosh(scoutingparticle.eta());
  float energy = std::sqrt(p*p + m*m);
  reco::Particle::LorentzVector p4(px, py, pz, energy); 

  static const reco::PFCandidate dummy; 
  auto pfcand = reco::PFCandidate(q, p4, dummy.translatePdgIdToType(scoutingparticle.pdgId()));

  bool relativeTrackVars = scoutingparticle.relative_trk_vars();
  vertexIndex.push_back(scoutingparticle.vertex());
  normchi2.push_back(scoutingparticle.normchi2());
  dz.push_back(scoutingparticle.dz());
  dxy.push_back(scoutingparticle.dxy());
  dzsig.push_back(scoutingparticle.dzsig());
  dxysig.push_back(scoutingparticle.dxysig());
  lostInnerHits.push_back(scoutingparticle.lostInnerHits());
  quality.push_back(scoutingparticle.quality());
  trkPt.push_back(relativeTrackVars ? scoutingparticle.trk_pt() + scoutingparticle.pt() : scoutingparticle.trk_pt());
  trkEta.push_back(relativeTrackVars ? scoutingparticle.trk_eta() + scoutingparticle.eta() : scoutingparticle.trk_eta());
  trkPhi.push_back(relativeTrackVars ? scoutingparticle.trk_phi() + scoutingparticle.phi() : scoutingparticle.trk_phi());

  return pfcand;  
}

void Run3ScoutingParticleToRecoPFCandidateProducer::createPFCandidates(edm::Handle<std::vector<Run3ScoutingParticle>> scoutingparticleHandle, std::unique_ptr<reco::PFCandidateCollection> & pfcands) {
  for (unsigned int icand = 0; icand < scoutingparticleHandle->size(); ++icand) {

   auto& scoutingparticle = (*scoutingparticleHandle)[icand];

   if (use_CHS_ and scoutingparticle.vertex() != 0) continue;

   auto pfcand = createPFCand(scoutingparticle);
   pfcands->push_back(pfcand);
  }
}

void Run3ScoutingParticleToRecoPFCandidateProducer::createPFCandidatesSK(edm::Handle<std::vector<Run3ScoutingParticle>> scoutingparticleHandle, std::unique_ptr<reco::PFCandidateCollection> & pfcands) {
  std::vector<fastjet::PseudoJet> fj;

  for (auto iter = scoutingparticleHandle->begin(), ibegin = scoutingparticleHandle->begin(), iend = scoutingparticleHandle->end(); iter != iend; ++iter) {
    auto m = pdTable_->particle(HepPDT::ParticleID(iter->pdgId())) != nullptr
                        ? pdTable_->particle(HepPDT::ParticleID(iter->pdgId()))->mass()
                        : -99.f;
    if (m < -90) continue;
    math::PtEtaPhiMLorentzVector p4(iter->pt(), iter->eta(), iter->phi(), m);
    fj.push_back(fastjet::PseudoJet(p4.px(), p4.py(), p4.pz(), p4.energy()));
    fj.back().set_user_index(iter - ibegin);
  } 

  fastjet::contrib::SoftKiller soft_killer(5, 0.4);
  std::vector<fastjet::PseudoJet> soft_killed_particles = soft_killer(fj);

  for (auto &particle : soft_killed_particles) {
    const Run3ScoutingParticle scoutingparticle = scoutingparticleHandle->at(particle.user_index());
    auto pfcand = createPFCand(scoutingparticle);
    pfcands->push_back(pfcand);
  }
}

// ------------ method called to produce the data  ------------
void Run3ScoutingParticleToRecoPFCandidateProducer::produce(edm::Event & iEvent, edm::EventSetup const & setup) {
  using namespace edm;

  auto pdt = setup.getHandle(particletable_token_);
  pdTable_ = pdt.product();

  Handle<std::vector<Run3ScoutingParticle>> scoutingparticleHandle;
  iEvent.getByToken(input_scoutingparticle_token_, scoutingparticleHandle);

  auto pfcands = std::make_unique<reco::PFCandidateCollection>();

  if (use_softKiller_) {
     createPFCandidatesSK(scoutingparticleHandle, pfcands);
  } else {
     createPFCandidates(scoutingparticleHandle, pfcands);
  }
  
  edm::OrphanHandle<reco::PFCandidateCollection> oh = iEvent.put(std::move(pfcands));
  
  std::unique_ptr<edm::ValueMap<int>> vertexIndex_VM(new edm::ValueMap<int>());
  edm::ValueMap<int>::Filler filler_vertexIndex(*vertexIndex_VM);
  filler_vertexIndex.insert(oh, vertexIndex.begin(), vertexIndex.end());
  filler_vertexIndex.fill();
  iEvent.put(std::move(vertexIndex_VM), "vertexIndex");
  
  std::unique_ptr<edm::ValueMap<float>> normchi2_VM(new edm::ValueMap<float>());
  edm::ValueMap<float>::Filler filler_normchi2(*normchi2_VM);
  filler_normchi2.insert(oh, normchi2.begin(), normchi2.end());
  filler_normchi2.fill();
  iEvent.put(std::move(normchi2_VM), "normchi2");

  std::unique_ptr<edm::ValueMap<float>> dz_VM(new edm::ValueMap<float>());
  edm::ValueMap<float>::Filler filler_dz(*dz_VM);
  filler_dz.insert(oh, dz.begin(), dz.end());
  filler_dz.fill();
  iEvent.put(std::move(dz_VM), "dz");

  std::unique_ptr<edm::ValueMap<float>> dxy_VM(new edm::ValueMap<float>());
  edm::ValueMap<float>::Filler filler_dxy(*dxy_VM);
  filler_dxy.insert(oh, dxy.begin(), dxy.end());
  filler_dxy.fill();
  iEvent.put(std::move(dxy_VM), "dxy");

  std::unique_ptr<edm::ValueMap<float>> dzsig_VM(new edm::ValueMap<float>());
  edm::ValueMap<float>::Filler filler_dzsig(*dzsig_VM);
  filler_dzsig.insert(oh, dzsig.begin(), dzsig.end());
  filler_dzsig.fill();
  iEvent.put(std::move(dzsig_VM), "dzsig");

  std::unique_ptr<edm::ValueMap<float>> dxysig_VM(new edm::ValueMap<float>());
  edm::ValueMap<float>::Filler filler_dxysig(*dxysig_VM);
  filler_dxysig.insert(oh, dxysig.begin(), dxysig.end());
  filler_dxysig.fill();
  iEvent.put(std::move(dxysig_VM), "dxysig");

  std::unique_ptr<edm::ValueMap<int>> lostInnerHits_VM(new edm::ValueMap<int>());
  edm::ValueMap<int>::Filler filler_lostInnerHits(*lostInnerHits_VM);
  filler_lostInnerHits.insert(oh, lostInnerHits.begin(), lostInnerHits.end());
  filler_lostInnerHits.fill();
  iEvent.put(std::move(lostInnerHits_VM), "lostInnerHits");

  std::unique_ptr<edm::ValueMap<int>> quality_VM(new edm::ValueMap<int>());
  edm::ValueMap<int>::Filler filler_quality(*quality_VM);
  filler_quality.insert(oh, quality.begin(), quality.end());
  filler_quality.fill();
  iEvent.put(std::move(quality_VM), "quality");

  std::unique_ptr<edm::ValueMap<float>> trkPt_VM(new edm::ValueMap<float>());
  edm::ValueMap<float>::Filler filler_trkPt(*trkPt_VM);
  filler_trkPt.insert(oh, trkPt.begin(), trkPt.end());
  filler_trkPt.fill();
  iEvent.put(std::move(trkPt_VM), "trkPt");

  std::unique_ptr<edm::ValueMap<float>> trkEta_VM(new edm::ValueMap<float>());
  edm::ValueMap<float>::Filler filler_trkEta(*trkEta_VM);
  filler_trkEta.insert(oh, trkEta.begin(), trkEta.end());
  filler_trkEta.fill();
  iEvent.put(std::move(trkEta_VM), "trkEta");

  std::unique_ptr<edm::ValueMap<float>> trkPhi_VM(new edm::ValueMap<float>());
  edm::ValueMap<float>::Filler filler_trkPhi(*trkPhi_VM);
  filler_trkPhi.insert(oh, trkPhi.begin(), trkPhi.end());
  filler_trkPhi.fill();
  iEvent.put(std::move(trkPhi_VM), "trkPhi");

  clearVars();
}

void Run3ScoutingParticleToRecoPFCandidateProducer::clearVars() { 
  vertexIndex.clear();
  normchi2.clear();
  dz.clear();
  dxy.clear();
  dzsig.clear();
  dxysig.clear();
  lostInnerHits.clear();
  quality.clear();
  trkPt.clear();
  trkEta.clear();
  trkPhi.clear();
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void Run3ScoutingParticleToRecoPFCandidateProducer::fillDescriptions(edm::ConfigurationDescriptions &descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("scoutingparticle", edm::InputTag("hltScoutingPFPacker"));
  desc.add<bool>("softKiller", false);
  desc.add<bool>("CHS", false);
  descriptions.addWithDefaultLabel(desc);
}

// declare this class as a framework plugin
DEFINE_FWK_MODULE(Run3ScoutingParticleToRecoPFCandidateProducer);
