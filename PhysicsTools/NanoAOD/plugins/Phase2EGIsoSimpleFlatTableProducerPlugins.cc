#include "PhysicsTools/NanoAOD/interface/SimpleFlatTableProducer.h"

#include "DataFormats/L1Trigger/interface/VertexWord.h"
typedef SimpleFlatTableProducer<l1t::VertexWord> SimpleL1TVertexWordFlatTableProducer;

#include "DataFormats/L1TCorrelator/interface/TkElectron.h"
typedef SimpleFlatTableProducer<l1t::TkElectron> SimpleL1TTkElectronFlatTableProducer;

#include "DataFormats/L1TCorrelator/interface/TkEm.h"
typedef SimpleFlatTableProducer<l1t::TkEm> SimpleL1TTkEmFlatTableProducer;

#include "DataFormats/L1TParticleFlow/interface/PFCandidate.h"
typedef SimpleFlatTableProducer<l1t::PFCandidate> SimpleL1TPFCandidateFlatTableProducer;

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(SimpleL1TVertexWordFlatTableProducer);
DEFINE_FWK_MODULE(SimpleL1TTkElectronFlatTableProducer);
DEFINE_FWK_MODULE(SimpleL1TTkEmFlatTableProducer);
DEFINE_FWK_MODULE(SimpleL1TPFCandidateFlatTableProducer);