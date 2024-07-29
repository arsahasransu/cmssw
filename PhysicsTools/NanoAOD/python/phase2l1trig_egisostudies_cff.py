import FWCore.ParameterSet.Config as cms
from PhysicsTools.NanoAOD.nano_eras_cff import *
from PhysicsTools.NanoAOD.common_cff import *
from PhysicsTools.NanoAOD.simpleCandidateFlatTableProducer_cfi import simpleCandidateFlatTableProducer

l1_float_precision_ = 20

P4Vars = cms.PSet(
    energy = Var("energy", float, precision=l1_float_precision_),
    pt = Var("pt", float, precision=l1_float_precision_),
    eta = Var("eta", float, precision=l1_float_precision_),
    phi = Var("phi", float, precision=l1_float_precision_)
)

V3Vars = cms.PSet(
    vx = Var("vx", float, precision=l1_float_precision_),
    vy = Var("vy", float, precision=l1_float_precision_),
    vz = Var("vz", float, precision=l1_float_precision_)
)


p2l1egisoFinalGenParticles = cms.EDProducer("GenParticlePruner",
    src = cms.InputTag("genParticles::HLT"),
    select = cms.vstring(
	"drop *",
        "keep++ abs(pdgId) == 15 & (pt > 15 ||  isPromptDecayed() )",#  keep full tau decay chain for some taus
	#"drop status==1 & pt < 1", #drop soft stable particle in tau decay
        "keep+ abs(pdgId) == 15 ",  #  keep first gen decay product for all tau
        "+keep pdgId == 22 && status == 1 && (pt > 10 || isPromptFinalState())", # keep gamma above 10 GeV (or all prompt) and its first parent
	"+keep abs(pdgId) == 11 || abs(pdgId) == 13 || abs(pdgId) == 15", #keep leptons, with at most one mother back in the history
	"drop abs(pdgId)= 2212 && abs(pz) > 1000", #drop LHC protons accidentally added by previous keeps
        "keep (400 < abs(pdgId) < 600) || (4000 < abs(pdgId) < 6000)", #keep all B and C hadrons
        "keep abs(pdgId) == 12 || abs(pdgId) == 14 || abs(pdgId) == 16",   # keep neutrinos
	"keep status == 3 || (status > 20 && status < 30)", #keep matrix element summary
        "keep isHardProcess() ||  fromHardProcessDecayed()  || fromHardProcessFinalState() || (statusFlags().fromHardProcess() && statusFlags().isLastCopy())",  #keep event summary based on status flags
	"keep  (status > 70 && status < 80 && pt > 15) ", # keep high pt partons right before hadronization
        "keep abs(pdgId) == 23 || abs(pdgId) == 24 || abs(pdgId) == 25 || abs(pdgId) == 37 ",   # keep VIP(articles)s
        #"keep abs(pdgId) == 310 && abs(eta) < 2.5 && pt > 1 ",                                                     # keep K0
        "keep (1000001 <= abs(pdgId) <= 1000039 ) || ( 2000001 <= abs(pdgId) <= 2000015)", #keep SUSY fiction particles
   )
)

p2l1egisoGenParticleTable = simpleCandidateFlatTableProducer.clone(
    src = cms.InputTag("p2l1egisoFinalGenParticles"),
    name = cms.string("GenPart"),
    doc = cms.string("all gen particles "),
    variables = cms.PSet(
         P4Vars, V3Vars,
         mass = Var("?!((abs(pdgId)>=1 && abs(pdgId)<=5) || (abs(pdgId)>=11 && abs(pdgId)<=16) || pdgId==21 || pdgId==111 || abs(pdgId)==211 || abs(pdgId)==421 || abs(pdgId)==411 || (pdgId==22 && mass<1))?mass:0", float,precision="?((abs(pdgId)==6 || abs(pdgId)>1000000) && statusFlags().isLastCopy())?20:8",doc="Mass stored for all particles with the exception of quarks (except top), leptons/neutrinos, photons with mass < 1 GeV, gluons, pi0(111), pi+(211), D0(421), and D+(411). For these particles, you can lookup the value from PDG."),
         pdgId  = Var("pdgId", int, doc="PDG id"),
         status  = Var("status", int, doc="Particle status. 1=stable"),
         charge = Var("charge", int, doc="charge of the particle"),
         genPartIdxMother = Var("?numberOfMothers>0?motherRef(0).key():-1", "int16", doc="index of the mother particle"),
         statusFlags = (Var(
            "statusFlags().isLastCopyBeforeFSR()                  * 16384 +"
            "statusFlags().isLastCopy()                           * 8192  +"
            "statusFlags().isFirstCopy()                          * 4096  +"
            "statusFlags().fromHardProcessBeforeFSR()             * 2048  +"
            "statusFlags().isDirectHardProcessTauDecayProduct()   * 1024  +"
            "statusFlags().isHardProcessTauDecayProduct()         * 512   +"
            "statusFlags().fromHardProcess()                      * 256   +"
            "statusFlags().isHardProcess()                        * 128   +"
            "statusFlags().isDirectHadronDecayProduct()           * 64    +"
            "statusFlags().isDirectPromptTauDecayProduct()        * 32    +"
            "statusFlags().isDirectTauDecayProduct()              * 16    +"
            "statusFlags().isPromptTauDecayProduct()              * 8     +"
            "statusFlags().isTauDecayProduct()                    * 4     +"
            "statusFlags().isDecayedLeptonHadron()                * 2     +"
            "statusFlags().isPrompt()                             * 1      ",
            "uint16", doc=("gen status flags stored bitwise, bits are: "
                "0 : isPrompt, "
                "1 : isDecayedLeptonHadron, "
                "2 : isTauDecayProduct, "
                "3 : isPromptTauDecayProduct, "
                "4 : isDirectTauDecayProduct, "
                "5 : isDirectPromptTauDecayProduct, "
                "6 : isDirectHadronDecayProduct, "
                "7 : isHardProcess, "
                "8 : fromHardProcess, "
                "9 : isHardProcessTauDecayProduct, "
                "10 : isDirectHardProcessTauDecayProduct, "
                "11 : fromHardProcessBeforeFSR, "
                "12 : isFirstCopy, "
                "13 : isLastCopy, "
                "14 : isLastCopyBeforeFSR, ")
            )),
    )
)

p2l1egisoL1EmuVtxTable = cms.EDProducer("SimpleL1TVertexWordFlatTableProducer",
    src = cms.InputTag("l1tVertexFinderEmulator", "L1VerticesEmulation:L1P2GT"),
    name = cms.string("L1EmuVtx"),
    doc = cms.string("L1 Emulated Vertices"),
    variables = cms.PSet(
        zL1EmuVtx = Var("z0", float),
        sumPtL1EmuVtx = Var("pt", float),
        qualityL1EmuVtx = Var("quality", int),
    ),
)

p2l1egisoL1EGTable = cms.EDProducer("SimpleTriggerL1EGFlatTableProducer",
    src = cms.InputTag("l1tEGammaClusterEmuProducer::L1P2GT"),
    minBX = cms.int32(-2),
    maxBX = cms.int32(2),                           
    cut = cms.string(""), 
    name= cms.string("L1EGClusEmu"),
    doc = cms.string(""),
    extension = cms.bool(False), 
    variables = cms.PSet(P4Vars)
)

p2l1egisoLayer1EGEETable = cms.EDProducer("SimpleTriggerL1EGFlatTableProducer",
    src = cms.InputTag("l1tLayer1EG", "L1EgEE:L1P2GT"),
    minBX = cms.int32(-2),
    maxBX = cms.int32(2),                           
    cut = cms.string(""), 
    name= cms.string("L1Layer1EGEE"),
    doc = cms.string(""),
    extension = cms.bool(False), 
    variables = cms.PSet(P4Vars)
)

IsolVars = cms.PSet(
    trkIsol = Var("trkIsol", float),
    trkIsolPV = Var("trkIsolPV", float),
    pfIsol = Var("pfIsol", float),
    pfIsolPV = Var("pfIsolPV", float),
    puppiIsol = Var("puppiIsol", float),
    puppiIsolPV = Var("puppiIsolPV", float)
)

p2l1egisoL1TkEleEBTable = cms.EDProducer("SimpleL1TTkElectronFlatTableProducer",
    src = cms.InputTag("l1tLayer1EG", "L1TkEleEB:L1P2GT"),
    name= cms.string("L1TkEleEB"),
    doc = cms.string("L1 Tk Electron in EB"),
    variables = cms.PSet(P4Vars, IsolVars,
                         trkzVtx = Var("trkzVtx", float)
                     )
)

p2l1egisoL1TkEleEETable = cms.EDProducer("SimpleL1TTkElectronFlatTableProducer",
    src = cms.InputTag("l1tLayer1EG", "L1TkEleEE:L1P2GT"),
    name= cms.string("L1TkEleEE"),
    doc = cms.string("L1 Tk Electron in EE"),
    variables = cms.PSet(P4Vars, IsolVars,
                         trkzVtx = Var("trkzVtx", float)
                     )
)

p2l1egisoL1TkEmEBTable = cms.EDProducer("SimpleL1TTkEmFlatTableProducer",
    src = cms.InputTag("l1tLayer1EG", "L1TkEmEB:L1P2GT"),
    name= cms.string("L1TkEmEB"),
    doc = cms.string("L1 Tk EM in EB"),
    variables = cms.PSet(P4Vars, IsolVars)
)

p2l1egisoL1TkEmEETable = cms.EDProducer("SimpleL1TTkEmFlatTableProducer",
    src = cms.InputTag("l1tLayer1EG", "L1TkEmEE:L1P2GT"),
    name= cms.string("L1TkEmEE"),
    doc = cms.string("L1 Tk EM in EE"),
    variables = cms.PSet(P4Vars, IsolVars)
)

L1PFCandVars = cms.PSet(
    charge = Var("charge", int, doc="charge of the particle"),
    dxy = Var("dxy", float, doc="transverse impact parameter"),
    z0 = Var("z0", float, doc="longitudnal impact parameter"),
    puppiWeight = Var("puppiWeight", float, doc="Puppi weight, -1 if unavailable"),
    pdgId = Var("pdgId", int, doc="PDG id"),
    clusPt = Var("?pfCluster.isNull()?-9999.0:pfCluster().pt()", float, doc="pt of the associated cluster"),
    clusEta = Var("?pfCluster.isNull()?-9999.0:pfCluster().eta()", float, doc="eta of the associated cluster"),
    clusPhi = Var("?pfCluster.isNull()?-9999.0:pfCluster().phi()", float, doc="phi of the associated cluster"),
)

p2l1egisoL1PFTable = cms.EDProducer("SimpleL1TPFCandidateFlatTableProducer",
    src = cms.InputTag("l1tLayer1", "PF:L1P2GT"),
    name= cms.string("L1PFCand"),
    doc = cms.string("L1 PF Objects"),
    variables = cms.PSet(P4Vars, V3Vars, L1PFCandVars)
)

p2l1egisoL1PuppiTable = cms.EDProducer("SimpleL1TPFCandidateFlatTableProducer",
    src = cms.InputTag("l1tLayer1", "Puppi:L1P2GT"),
    name= cms.string("L1PuppiCand"),
    doc = cms.string("L1 Puppi Objects"),
    variables = cms.PSet(P4Vars, V3Vars, L1PFCandVars)
)

phase2L1EGIsoTablesTask = cms.Task(p2l1egisoFinalGenParticles, p2l1egisoGenParticleTable, 
                                   p2l1egisoL1EmuVtxTable,
                                   p2l1egisoL1EGTable, p2l1egisoLayer1EGEETable,
                                   p2l1egisoL1TkEleEBTable, p2l1egisoL1TkEleEETable,
                                   p2l1egisoL1TkEmEBTable, p2l1egisoL1TkEmEETable,
                                   p2l1egisoL1PFTable, p2l1egisoL1PuppiTable)