#ifndef DataFormatsL1TCorrelator_TkEm_h
#define DataFormatsL1TCorrelator_TkEm_h

// -*- C++ -*-
//
// Package:     L1Trigger
// Class  :     TkEm
//

#include "DataFormats/L1Trigger/interface/L1Candidate.h"
#include "DataFormats/Common/interface/Ptr.h"

#include "DataFormats/L1Trigger/interface/EGamma.h"

#include "DataFormats/L1TrackTrigger/interface/TTTypes.h"
#include "DataFormats/L1TParticleFlow/interface/gt_datatypes.h"
#include "DataFormats/L1TParticleFlow/interface/egamma.h"
#include "FWCore/Utilities/interface/Exception.h"

#include <ap_int.h>

namespace l1t {

  class TkEm : public L1Candidate {
  public:
    TkEm();

    TkEm(const LorentzVector& p4, float tkisol = -999., float tkisolPV = -999);

    TkEm(const LorentzVector& p4, const edm::Ptr<L1Candidate>& egCaloPtr, float tkisol = -999., float tkisolPV = -999);

    enum class HWEncoding { None, CT, GT };

    // ---------- const member functions ---------------------

    const edm::Ptr<L1Candidate>& egCaloPtr() const { return egCaloPtr_; }

    float trkIsol() const { return trkIsol_; }          // not constrained to the PV, just track ptSum
    float trkIsolPV() const { return trkIsolPV_; }      // constrained to the PV by DZ
    float pfIsol() const { return pfIsol_; }            // not constrained to the PV, just track ptSum
    float pfIsolPV() const { return pfIsolPV_; }        // constrained to the PV by DZ
    float puppiIsol() const { return puppiIsol_; }      // not constrained to the PV, just track ptSum
    float puppiIsolPV() const { return puppiIsolPV_; }  // constrained to the PV by DZ

    // ---------- member functions ---------------------------

    void setTrkIsol(float TrkIsol) { trkIsol_ = TrkIsol; }
    void setTrkIsolPV(float TrkIsolPV) { trkIsolPV_ = TrkIsolPV; }
    void setPFIsol(float pfIsol) { pfIsol_ = pfIsol; }
    void setPFIsolPV(float pfIsolPV) { pfIsolPV_ = pfIsolPV; }
    void setPuppiIsol(float puppiIsol) { puppiIsol_ = puppiIsol; }
    void setPuppiIsolPV(float puppiIsolPV) { puppiIsolPV_ = puppiIsolPV; }
    void setEgCaloPtr(const edm::Ptr<L1Candidate>& egPtr) { egCaloPtr_ = egPtr; }

    void setHwCalo(const l1ct::EGIsoObj& eg) {
      hwCaloEta_ = eg.intEta();
      hwCaloPhi_ = eg.intPhi();
    }

    template <int N>
    void setEgBinaryWord(ap_uint<N> word, HWEncoding encoding) {
      egBinaryWord0_ = word;
      egBinaryWord1_ = (word >> 32);
      egBinaryWord2_ = (word >> 64);
      encoding_ = encoding;
    }

    l1gt::Photon hwObj() const {
      if (encoding() != HWEncoding::GT) {
        throw cms::Exception("RuntimeError") << "TkEm::hwObj : encoding is not in GT format!" << std::endl;
      }
      return l1gt::Photon::unpack_ap(egBinaryWord<l1gt::Photon::BITWIDTH>());
    }

    // Override the legacy L1Candidate int accessors: in the Phase2 flow those int
    // members are never filled (they always read out 0). Instead, decode the hardware
    // (integer fixed-point) kinematics from the packed EG binary word, which is set by
    // the Layer1/Layer2 producers. Returns 0 if the word has not been set (None).
    int hwPt() const {
      if (encoding() == HWEncoding::CT) {
        return l1ct::EGIsoObj::unpack(egBinaryWord<l1ct::EGIsoObj::BITWIDTH>()).intPt();
      }
      if (encoding() == HWEncoding::GT) {
        const l1gt::pt_t pt = hwObj().v3.pt;
        return ap_uint<l1gt::pt_t::width>(pt.range()).to_int();
      }
      return 0;
    }
    int hwEta() const {  // eta at the calo face (CT: global beta) / (GT: global eta)
      if (encoding() == HWEncoding::CT) {
        return l1ct::EGIsoObj::unpack(egBinaryWord<l1ct::EGIsoObj::BITWIDTH>()).intEta();
      }
      if (encoding() == HWEncoding::GT) {
        return hwObj().v3.eta.to_int();
      }
      return 0;
    }
    int hwPhi() const {  // phi at the calo face (CT: global phi) / (GT: global phi)
      if (encoding() == HWEncoding::CT) {
        return l1ct::EGIsoObj::unpack(egBinaryWord<l1ct::EGIsoObj::BITWIDTH>()).intPhi();
      }
      if (encoding() == HWEncoding::GT) {
        return hwObj().v3.phi.to_int();
      }
      return 0;
    }

    // Digitized eta/phi of the matched EM cluster at the calo face (global coords,
    // LSB = pi/720). Set from the l1ct::EGIsoObj at object creation (photons & electrons).
    int hwCaloEta() const { return hwCaloEta_; }
    int hwCaloPhi() const { return hwCaloPhi_; }

    template <int N>
    ap_uint<N> egBinaryWord() const {
      return ap_uint<N>(egBinaryWord0_) | (ap_uint<N>(egBinaryWord1_) << 32) | (ap_uint<N>(egBinaryWord2_) << 64);
    }

    HWEncoding encoding() const { return encoding_; }

  private:
    edm::Ptr<L1Candidate> egCaloPtr_;
    float trkIsol_;
    float trkIsolPV_;
    float pfIsol_;
    float pfIsolPV_;
    float puppiIsol_;
    float puppiIsolPV_;
    uint32_t egBinaryWord0_;
    uint32_t egBinaryWord1_;
    uint32_t egBinaryWord2_;
    HWEncoding encoding_;
    int hwCaloEta_{0};
    int hwCaloPhi_{0};
  };
}  // namespace l1t

#endif
