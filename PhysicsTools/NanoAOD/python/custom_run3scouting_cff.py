import FWCore.ParameterSet.Config as cms
from  PhysicsTools.NanoAOD.common_cff import *
from PhysicsTools.NanoAOD.simpleCandidateFlatTableProducer_cfi import simpleCandidateFlatTableProducer

def AddRun3Scouting(process):
  """
  Add all Run 3 scouting objects except jets and PFCandidates
  """
  process.photonScoutingTable = cms.EDProducer("SimpleRun3ScoutingPhotonFlatTableProducer",
       src = cms.InputTag("hltScoutingEgammaPacker"),
       cut = cms.string(""),
       name = cms.string("ScoutingPhoton"),
       doc  = cms.string("Photon Scouting Informations"),
       singleton = cms.bool(False), # the number of entries is variable
       extension = cms.bool(False), # this is the main table for the muons
       variables = cms.PSet(
           pt = Var('pt', 'float', precision=14, doc='super-cluster (SC) pt'),
           eta = Var('eta', 'float', precision=14, doc='SC eta'),
           phi = Var('phi', 'float', precision=14, doc='SC phi'),
           m = Var('m', 'float', precision=14, doc='SC mass'),
           sigmaIetaIeta = Var('sigmaIetaIeta', 'float', precision=14, doc='sigmaIetaIeta of the SC, calculated with full 5x5 region, noise cleaned'),
           hOverE = Var('hOverE', 'float', precision=14, doc='Energy in HCAL / Energy in ECAL'),
           ecalIso = Var('ecalIso', 'float', precision=14, doc='Isolation of SC in the ECAL'),
           hcalIso = Var('hcalIso', 'float', precision=14, doc='Isolation of SC in the HCAL'),
           r9 = Var('r9', 'float', precision=14, doc='Photon SC r9 as defined in https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideEgammaShowerShape'),
           sMin = Var('sMin', 'float', precision=14, doc='minor moment of the SC shower shape'),
           sMaj = Var('sMaj', 'float', precision=14, doc='major moment of the SC shower shape'),
           seedId = Var('seedId', 'int', doc='ECAL ID of the SC seed'),
       )
  )

  process.electronScoutingTable = cms.EDProducer("SimpleRun3ScoutingElectronFlatTableProducer",
       src = cms.InputTag("hltScoutingEgammaPacker"),
       cut = cms.string(""),
       name = cms.string("ScoutingElectron"),
       doc  = cms.string("Electron Scouting Informations"),
       singleton = cms.bool(False), # the number of entries is variable
       extension = cms.bool(False), # this is the main table for the muons
       variables = cms.PSet(
           pt = Var('pt', 'float', precision=14, doc='super-cluster (SC) pt'),
           eta = Var('eta', 'float', precision=14, doc='SC eta'),
           phi = Var('phi', 'float', precision=14, doc='SC phi'),
           m = Var('m', 'float', precision=14, doc='SC mass'),
           d0 = Var('d0', 'float', precision=14, doc='track d0'),
           dz = Var('dz', 'float', precision=14, doc='track dz'),
           dEtaIn = Var('dEtaIn', 'float', precision=14, doc='#Delta#eta(SC seed, track pixel seed)'),
           dPhiIn = Var('dPhiIn', 'float', precision=14, doc='#Delta#phi(SC seed, track pixel seed)'),
           sigmaIetaIeta = Var('sigmaIetaIeta', 'float', precision=14, doc='sigmaIetaIeta of the SC, calculated with full 5x5 region, noise cleaned'),
           hOverE = Var('hOverE', 'float', precision=14, doc='Energy in HCAL / Energy in ECAL'),
           ooEMOop = Var('ooEMOop', 'float', precision=14, doc='1/E(SC) - 1/p(track momentum)'),
           missingHits = Var('missingHits', 'int', doc='missing hits in the tracker'),
           charge = Var('charge', 'int', doc='track charge'),
           ecalIso = Var('ecalIso', 'float', precision=14, doc='Isolation of SC in the ECAL'),
           hcalIso = Var('hcalIso', 'float', precision=14, doc='Isolation of SC in the HCAL'),
           trackIso = Var('trackIso', 'float', precision=14, doc='Isolation of electron track in the tracker'),
           r9 = Var('r9', 'float', precision=14, doc='ELectron SC r9 as defined in https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideEgammaShowerShape'),
           sMin = Var('sMin', 'float', precision=14, doc='minor moment of the SC shower shape'),
           sMaj = Var('sMaj', 'float', precision=14, doc='major moment of the SC shower shape'),
           seedId = Var('seedId', 'int', doc='ECAL ID of the SC seed'),
       )
  )

  process.muonScoutingTable = cms.EDProducer("SimpleRun3ScoutingMuonFlatTableProducer",
       src = cms.InputTag("hltScoutingMuonPacker"),
       cut = cms.string(""),
       name = cms.string("ScoutingMuon"),
       doc  = cms.string("Muon Scouting Informations"),
       singleton = cms.bool(False), # the number of entries is variable
       extension = cms.bool(False), # this is the main table for the muons
       variables = cms.PSet(
           pt = Var('pt', 'float', precision=14, doc='pt'),
           eta = Var('eta', 'float', precision=14, doc='eta'),
           phi = Var('phi', 'float', precision=14, doc='phi'),
           m = Var('m', 'float', precision=14, doc='m'),
           type = Var('type', 'int', doc='type'),
           charge = Var('charge', 'int', doc='charge'),
           normchi2 = Var('normalizedChi2', 'float', precision=14, doc='normalizedChi2'),
           ecalIso = Var('ecalIso', 'float', precision=14, doc='ecalIso'),
           hcalIso = Var('hcalIso', 'float', precision=14, doc='hcalIso'),
           nValidStandAloneMuonHits = Var('nValidStandAloneMuonHits', 'int', doc='nValidStandAloneMuonHits'),
           nStandAloneMuonMatchedStations = Var('nStandAloneMuonMatchedStations', 'int', doc='nStandAloneMuonMatchedStations'),
           nValidRecoMuonHits = Var('nValidRecoMuonHits', 'int', doc='nValidRecoMuonHits'),
           nRecoMuonChambers = Var('nRecoMuonChambers', 'int', doc='nRecoMuonChambers'),
           nRecoMuonChambersCSCorDT = Var('nRecoMuonChambersCSCorDT', 'int', doc='nRecoMuonChambersCSCorDT'),
           nRecoMuonMatches = Var('nRecoMuonMatches', 'int', doc='nRecoMuonMatches'),
           nRecoMuonMatchedStations = Var('nRecoMuonMatchedStations', 'int', doc='nRecoMuonMatchedStations'),
           nRecoMuonExpectedMatchedStations = Var('nRecoMuonExpectedMatchedStations', 'int', doc='nRecoMuonExpectedMatchedStations'),
           recoMuonStationMask = Var('recoMuonStationMask', 'int', doc='recoMuonStationMask'),
           nRecoMuonMatchedRPCLayers = Var('nRecoMuonMatchedRPCLayers', 'int', doc='nRecoMuonMatchedRPCLayers'),
           recoMuonRPClayerMask = Var('recoMuonRPClayerMask', 'int', doc='recoMuonRPClayerMask'),
           nValidPixelHits = Var('nValidPixelHits', 'int', doc='nValidPixelHits'),
           nValidStripHits = Var('nValidStripHits', 'int', doc='nValidStripHits'),
           nPixelLayersWithMeasurement = Var('nPixelLayersWithMeasurement', 'int', doc='nPixelLayersWithMeasurement'),
           nTrackerLayersWithMeasurement = Var('nTrackerLayersWithMeasurement', 'int', doc='nTrackerLayersWithMeasurement'),
           trk_chi2 = Var('trk_chi2', 'float', precision=14, doc='trk_chi2'),
           trk_ndof = Var('trk_ndof', 'float', precision=14, doc='trk_ndof'),
           trk_dxy = Var('trk_dxy', 'float', precision=14, doc='trk_dxy'),
           trk_dz = Var('trk_dz', 'float', precision=14, doc='trk_dz'),
           trk_qoverp = Var('trk_qoverp', 'float', precision=14, doc='trk_qoverp'),
           trk_lambda = Var('trk_lambda', 'float', precision=14, doc='trk_lambda'),
           trk_pt = Var('trk_pt', 'float', precision=14, doc='trk_pt'),
           trk_phi = Var('trk_phi', 'float', precision=14, doc='trk_phi'),
           trk_eta = Var('trk_eta', 'float', precision=14, doc='trk_eta'),
           trk_dxyError = Var('trk_dxyError', 'float', precision=14, doc='trk_dxyError'),
           trk_dzError = Var('trk_dzError', 'float', precision=14, doc='trk_dzError'),
           trk_qoverpError = Var('trk_qoverpError', 'float', precision=14, doc='trk_qoverpError'),
           trk_lambdaError = Var('trk_lambdaError', 'float', precision=14, doc='trk_lambdaError'),
           trk_phiError = Var('trk_phiError', 'float', precision=14, doc='trk_phiError'),
           trk_dsz = Var('trk_dsz', 'float', precision=14, doc='trk_dsz'),
           trk_dszError = Var('trk_dszError', 'float', precision=14, doc='trk_dszError'),
           trk_qoverp_lambda_cov = Var('trk_qoverp_lambda_cov', 'float', precision=14, doc='trk_qoverp_lambda_cov'),
           trk_qoverp_phi_cov = Var('trk_qoverp_phi_cov', 'float', precision=14, doc='trk_qoverp_phi_cov'),
           trk_qoverp_dxy_cov = Var('trk_qoverp_dxy_cov', 'float', precision=14, doc='trk_qoverp_dxy_cov'),
           trk_qoverp_dsz_cov = Var('trk_qoverp_dsz_cov', 'float', precision=14, doc='trk_qoverp_dsz_cov'),
           trk_lambda_phi_cov = Var('trk_lambda_phi_cov', 'float', precision=14, doc='trk_lambda_phi_cov'),
           trk_lambda_dxy_cov = Var('trk_lambda_dxy_cov', 'float', precision=14, doc='trk_lambda_dxy_cov'),
           trk_lambda_dsz_cov = Var('trk_lambda_dsz_cov', 'float', precision=14, doc='trk_lambda_dsz_cov'),
           trk_phi_dxy_cov = Var('trk_phi_dxy_cov', 'float', precision=14, doc='trk_phi_dxy_cov'),
           trk_phi_dsz_cov = Var('trk_phi_dsz_cov', 'float', precision=14, doc='trk_phi_dsz_cov'),
           trk_dxy_dsz_cov = Var('trk_dxy_dsz_cov', 'float', precision=14, doc='trk_dxy_dsz_cov'),
           trk_vx = Var('trk_vx', 'float', precision=14, doc='trk_vx'),
           trk_vy = Var('trk_vy', 'float', precision=14, doc='trk_vy'),
           trk_vz = Var('trk_vz', 'float', precision=14, doc='trk_vz'),
       )
  )

  process.trackScoutingTable = cms.EDProducer("SimpleRun3ScoutingTrackFlatTableProducer",
       src = cms.InputTag("hltScoutingTrackPacker"),
       cut = cms.string(""),
       name = cms.string("ScoutingTrack"),
       doc  = cms.string("Track Scouting Informations"),
       singleton = cms.bool(False), # the number of entries is variable
       extension = cms.bool(False), # this is the main table for the muons
       variables = cms.PSet(
           pt = Var('tk_pt', 'float', precision=14, doc='pt'),
           eta = Var('tk_eta', 'float', precision=14, doc='eta'),
           phi = Var('tk_phi', 'float', precision=14, doc='phi'),
           chi2 = Var('tk_chi2', 'float', precision=14, doc='chi2'),
           ndof = Var('tk_ndof', 'float', precision=14, doc='ndof'),
           charge = Var('tk_charge', 'int', doc='charge'),
           dxy = Var('tk_dxy', 'float', precision=14, doc='dxy'),
           dz = Var('tk_dz', 'float', precision=14, doc='dz'),
           nValidPixelHits = Var('tk_nValidPixelHits', 'int', doc='nValidPixelHits'),
           nValidStripHits = Var('tk_nValidStripHits', 'int', doc='nValidStripHits'),
           nTrackerLayersWithMeasurement = Var('tk_nTrackerLayersWithMeasurement', 'int', doc='nTrackerLayersWithMeasurement'),
           qoverp = Var('tk_qoverp', 'float', precision=14, doc='qoverp'),
           _lambda = Var('tk_lambda', 'float', precision=14, doc='lambda'),
           dxyError = Var('tk_dxy_Error', 'float', precision=14, doc='dxyError'),
           dzError = Var('tk_dz_Error', 'float', precision=14, doc='dzError'),
           qoverpError = Var('tk_qoverp_Error', 'float', precision=14, doc='qoverpError'),
           lambdaError = Var('tk_lambda_Error', 'float', precision=14, doc='lambdaError'),
           phiError = Var('tk_phi_Error', 'float', precision=14, doc='phiError'),
           dsz = Var('tk_dsz', 'float', precision=14, doc='dsz'),
           dszError = Var('tk_dsz_Error', 'float', precision=14, doc='dszError'),
           qoverp_lambda_cov = Var('tk_qoverp_lambda_cov', 'float', precision=14, doc='qoverp_lambda_cov'),
           qoverp_phi_cov = Var('tk_qoverp_phi_cov', 'float', precision=14, doc='qoverp_phi_cov'),
           qoverp_dxy_cov = Var('tk_qoverp_dxy_cov', 'float', precision=14, doc='qoverp_dxy_cov'),
           qoverp_dsz_cov = Var('tk_qoverp_dsz_cov', 'float', precision=14, doc='qoverp_dsz_cov'),
           lambda_phi_cov = Var('tk_lambda_phi_cov', 'float', precision=14, doc='lambda_phi_cov'),
           lambda_dxy_cov = Var('tk_lambda_dxy_cov', 'float', precision=14, doc='lambda_dxy_cov'),
           lambda_dsz_cov = Var('tk_lambda_dsz_cov', 'float', precision=14, doc='lambda_dsz_cov'),
           phi_dxy_cov = Var('tk_phi_dxy_cov', 'float', precision=14, doc='phi_dxy_cov'),
           phi_dsz_cov = Var('tk_phi_dsz_cov', 'float', precision=14, doc='phi_dsz_cov'),
           dxy_dsz_cov = Var('tk_dxy_dsz_cov', 'float', precision=14, doc='dxy_dsz_cov'),
           vtxInd = Var('tk_vtxInd', 'int', doc='vtxInd'),
           vx = Var('tk_vx', 'float', precision=14, doc='vx'),
           vy = Var('tk_vy', 'float', precision=14, doc='vy'),
           vz = Var('tk_vz', 'float', precision=14, doc='vz'),
       )
  )

  process.primaryvertexScoutingTable = cms.EDProducer("SimpleRun3ScoutingVertexFlatTableProducer",
       src = cms.InputTag("hltScoutingPrimaryVertexPacker", "primaryVtx"),
       cut = cms.string(""),
       name = cms.string("ScoutingPrimaryVertex"),
       doc  = cms.string("PrimaryVertex Scouting Informations"),
       singleton = cms.bool(False), # the number of entries is variable
       extension = cms.bool(False), # this is the main table for the muons
       variables = cms.PSet(
           x = Var('x', 'float', precision=14, doc='x'),
           y = Var('y', 'float', precision=14, doc='y'),
           z = Var('z', 'float', precision=14, doc='z'),
           xError = Var('xError', 'float', precision=14, doc='xError'),
           yError = Var('yError', 'float', precision=14, doc='yError'),
           zError = Var('zError', 'float', precision=14, doc='zError'),
           tracksSize = Var('tracksSize', 'int', doc='tracksSize'),
           chi2 = Var('chi2', 'float', precision=14, doc='chi2'),
           ndof = Var('ndof', 'int', doc='ndof'),
           isValidVtx = Var('isValidVtx', 'bool', doc='isValidVtx'),
       )
  )

  process.displacedvertexScoutingTable = cms.EDProducer("SimpleRun3ScoutingVertexFlatTableProducer",
       src = cms.InputTag("hltScoutingMuonPacker","displacedVtx"),
       cut = cms.string(""),
       name = cms.string("ScoutingDisplacedVertex"),
       doc  = cms.string("DisplacedVertex Scouting Informations"),
       singleton = cms.bool(False), # the number of entries is variable
       extension = cms.bool(False), # this is the main table for the muons
       variables = cms.PSet(
           x = Var('x', 'float', precision=14, doc='x'),
           y = Var('y', 'float', precision=14, doc='y'),
           z = Var('z', 'float', precision=14, doc='z'),
           xError = Var('xError', 'float', precision=14, doc='xError'),
           yError = Var('yError', 'float', precision=14, doc='yError'),
           zError = Var('zError', 'float', precision=14, doc='zError'),
           tracksSize = Var('tracksSize', 'int', doc='tracksSize'),
           chi2 = Var('chi2', 'float', precision=14, doc='chi2'),
           ndof = Var('ndof', 'int', doc='ndof'),
           isValidVtx = Var('isValidVtx', 'bool', doc='isValidVtx'),
       )
  )

  process.rhoScoutingTable = cms.EDProducer("GlobalVariablesTableProducer",
      name = cms.string(""),
      variables = cms.PSet(
          ScoutingRho = ExtVar( cms.InputTag("hltScoutingPFPacker", "rho"), "double", doc = "rho from all PF Candidates, no foreground removal (for isolation of prompt photons)" ),
      )
  )

  process.metScoutingTable = cms.EDProducer("GlobalVariablesTableProducer",
      name = cms.string("ScoutingMET"),
      variables = cms.PSet(
          pt = ExtVar( cms.InputTag("hltScoutingPFPacker", "pfMetPt"), "double", doc = "Scouting MET pt"),
          phi = ExtVar( cms.InputTag("hltScoutingPFPacker", "pfMetPhi"), "double", doc = "Scouting MET phi"),
      )
  )

  run3ScoutingTaskName = "run3ScoutingTask"
  setattr(process, run3ScoutingTaskName, cms.Task(
      getattr(process,"photonScoutingTable"),
      getattr(process,"muonScoutingTable"),
      getattr(process,"electronScoutingTable"),
      getattr(process,"trackScoutingTable"),
      getattr(process,"primaryvertexScoutingTable"),
      getattr(process,"displacedvertexScoutingTable"),
      getattr(process,"rhoScoutingTable"),
      getattr(process,"metScoutingTable"),
    )
  )
  process.customRun3ScoutingNanoAODTask.add(getattr(process,run3ScoutingTaskName))
  return process

def ConvertScoutingToReco(process):
  """
  Convert Run 3 scouting particles to recoPFCandidates
  """
  process.pfcands = cms.EDProducer(
     "Run3ScoutingParticleToRecoPFCandidateProducer",
     scoutingparticle=cms.InputTag("hltScoutingPFPacker"),
   )

  process.customRun3ScoutingNanoAODTask.add(cms.Task(getattr(process,"pfcands")))
  return process

def AddParticles(process):
  """
  Add particles
  """
  process.particleTable = cms.EDProducer("SimpleCandidateFlatTableProducer",
    src = cms.InputTag("pfcands"),
    name = cms.string("ScoutingParticle"),
    cut = cms.string(""),
    doc = cms.string("ScoutingParticle"),
    singleton = cms.bool(False),
    extension = cms.bool(False), # this is the main table
    externalVariables = cms.PSet(
       vertexIndex = ExtVar(cms.InputTag("pfcands", "vertexIndex"), int, doc="vertex index"),
       trkNormchi2 = ExtVar(cms.InputTag("pfcands", "normchi2"), float, doc="normchi2 of best track", precision=6),
       trkDz = ExtVar(cms.InputTag("pfcands", "dz"), float, doc="dz of best track", precision=6),
       trkDxy = ExtVar(cms.InputTag("pfcands", "dxy"), float, doc="dxy of best track", precision=6),
       trkDzsig = ExtVar(cms.InputTag("pfcands", "dzsig"), float, doc="dzsig of best track", precision=6),
       trkDxysig = ExtVar(cms.InputTag("pfcands", "dxysig"), float, doc="dxysig of best track", precision=6),
       trkLostInnerHits = ExtVar(cms.InputTag("pfcands", "lostInnerHits"), int, doc="lostInnerHits of best track"),
       trkQuality = ExtVar(cms.InputTag("pfcands", "quality"), int, doc="quality of best track"),
       trkPt = ExtVar(cms.InputTag("pfcands", "trkPt"), float, doc="pt of best track", precision=6),
       trkEta = ExtVar(cms.InputTag("pfcands", "trkEta"), float, doc="eta of best track", precision=6),
       trkPhi = ExtVar(cms.InputTag("pfcands", "trkPhi"), float, doc="phi of best track", precision=6),
    ),
    variables = cms.PSet(
       CandVars,
    ),
  )
  process.customRun3ScoutingNanoAODTask.add(cms.Task(getattr(process,"particleTable")))
  return process

def AddAK4PFJets(process):
  """
  Add AK4 PF jets
  """
  from RecoJets.JetProducers.ak4PFJets_cfi import ak4PFJets
  process.ak4Jets = ak4PFJets.clone(
     src = ("pfcands"),
  )

  process.ak4ParticleNetJetTagInfos = cms.EDProducer("DeepBoostedJetTagInfoProducer",
      jet_radius = cms.double( 0.4 ),
      min_jet_pt = cms.double( 5.0 ),
      max_jet_eta = cms.double( 2.5 ),
      min_pt_for_track_properties = cms.double( 0.95 ),
      min_pt_for_pfcandidates = cms.double( 0.1 ),
      use_puppiP4 = cms.bool( False ),
      include_neutrals = cms.bool( True ),
      sort_by_sip2dsig = cms.bool( False ),
      min_puppi_wgt = cms.double( -1.0 ),
      flip_ip_sign = cms.bool( False ),
      sip3dSigMax = cms.double( -1.0 ),
      use_hlt_features = cms.bool( False ),
      pf_candidates = cms.InputTag( "pfcands" ),
      jets = cms.InputTag( "ak4Jets" ),
      puppi_value_map = cms.InputTag( "" ),
      use_scouting_features = cms.bool( True ),
      normchi2_value_map = cms.InputTag("pfcands", "normchi2"),
      dz_value_map = cms.InputTag("pfcands", "dz"),
      dxy_value_map = cms.InputTag("pfcands", "dxy"),
      dzsig_value_map = cms.InputTag("pfcands", "dzsig"),
      dxysig_value_map = cms.InputTag("pfcands", "dxysig"),
      lostInnerHits_value_map = cms.InputTag("pfcands", "lostInnerHits"),
      quality_value_map = cms.InputTag("pfcands", "quality"),
      trkPt_value_map = cms.InputTag("pfcands", "trkPt"),
      trkEta_value_map = cms.InputTag("pfcands", "trkEta"),
      trkPhi_value_map = cms.InputTag("pfcands", "trkPhi"),
  )

  from RecoBTag.ONNXRuntime.boostedJetONNXJetTagsProducer_cfi import boostedJetONNXJetTagsProducer
  process.ak4ParticleNetJetTags = cms.EDProducer("BoostedJetONNXValueMapProducer",
      jets = cms.InputTag("ak4Jets"),
      src = cms.InputTag("ak4ParticleNetJetTagInfos"),
      preprocess_json = cms.string("Models/preprocess_flavourtag.json"),
      model_path = cms.FileInPath("Models/flavourtag.onnx"),
      flav_names = cms.vstring(["probb", "probbb","probc", "probcc", "probuds", "probg", "probundef"]),
      debugMode = cms.untracked.bool(False),
  )

  process.ak4JetTable = cms.EDProducer("SimpleCandidateFlatTableProducer",
      src = cms.InputTag("ak4Jets"),
      name = cms.string("ScoutingJet"),
      cut = cms.string(""),
      doc = cms.string("ScoutingJet"),
      singleton = cms.bool(False),
      extension = cms.bool(False), # this is the main table
      externalVariables = cms.PSet(
         particleNet_prob_b = ExtVar(cms.InputTag('ak4ParticleNetJetTags:probb'), float, doc="ParticleNet prob b", precision=10),
         particleNet_prob_bb = ExtVar(cms.InputTag('ak4ParticleNetJetTags:probbb'), float, doc="ParticleNet prob b", precision=10),
         particleNet_prob_c = ExtVar(cms.InputTag('ak4ParticleNetJetTags:probc'), float, doc="ParticleNet prob c", precision=10),
         particleNet_prob_cc = ExtVar(cms.InputTag('ak4ParticleNetJetTags:probcc'), float, doc="ParticleNet prob cc", precision=10),
         particlenet_prob_uds = ExtVar(cms.InputTag('ak4ParticleNetJetTags:probuds'), float, doc="particlenet prob uds", precision=10),
         particleNet_prob_g = ExtVar(cms.InputTag('ak4ParticleNetJetTags:probg'), float, doc="ParticleNet prob g", precision=10),
         particleNet_prob_undef = ExtVar(cms.InputTag('ak4ParticleNetJetTags:probundef'), float, doc="ParticleNet prob undef", precision=10),
      ),
      variables = cms.PSet(
         P4Vars,
         area = Var("jetArea()", float, doc="jet catchment area, for JECs",precision=10),
         chHEF = Var("chargedHadronEnergy()/(chargedHadronEnergy()+neutralHadronEnergy()+photonEnergy()+electronEnergy()+muonEnergy())", float, doc="charged Hadron Energy Fraction", precision= 6),
         neHEF = Var("neutralHadronEnergy()/(chargedHadronEnergy()+neutralHadronEnergy()+photonEnergy()+electronEnergy()+muonEnergy())", float, doc="neutral Hadron Energy Fraction", precision= 6),
         chEmEF = Var("(electronEnergy()+muonEnergy())/(chargedHadronEnergy()+neutralHadronEnergy()+photonEnergy()+electronEnergy()+muonEnergy())", float, doc="charged Electromagnetic Energy Fraction", precision= 6),
         neEmEF = Var("(photonEnergy())/(chargedHadronEnergy()+neutralHadronEnergy()+photonEnergy()+electronEnergy()+muonEnergy())", float, doc="neutral Electromagnetic Energy Fraction", precision= 6),
         muEmEF = Var("(muonEnergy())/(chargedHadronEnergy()+neutralHadronEnergy()+photonEnergy()+electronEnergy()+muonEnergy())", float, doc="muon Energy Fraction", precision= 6),
         nCh = Var("chargedHadronMultiplicity()", int, doc="number of charged hadrons in the jet"),
         nNh = Var("neutralHadronMultiplicity()", int, doc="number of neutral hadrons in the jet"),
         nMuons = Var("muonMultiplicity()", int, doc="number of muons in the jet"),
         nElectrons = Var("electronMultiplicity()", int, doc="number of electrons in the jet"),
         nPhotons = Var("photonMultiplicity()", int, doc="number of photons in the jet"),
         nConstituents = Var("numberOfDaughters()", "uint8", doc="Number of particles in the jet")
      ),
  )

  ak4PFJetTaskName = "ak4PFJetTask"
  setattr(process, ak4PFJetTaskName, cms.Task(
      getattr(process,"ak4Jets"),
      getattr(process,"ak4ParticleNetJetTagInfos"),
      getattr(process,"ak4ParticleNetJetTags"),
      getattr(process,"ak4JetTable"),
    )
  )
  process.customRun3ScoutingNanoAODTask.add(getattr(process,ak4PFJetTaskName))
  return process

def AddAK8PFJets(process):
  """
  Add AK8 PF jets
  """
  from RecoJets.JetProducers.ak4PFJets_cfi import ak4PFJets
  process.ak8Jets = ak4PFJets.clone(
     src = ("pfcands"),
     rParam   = 0.8,
     jetPtMin = 50.0,
  )

  process.ak8JetsSoftDrop = ak4PFJets.clone(
     src = ("pfcands"),
     rParam   = 0.8,
     jetPtMin = 50.0,
     useSoftDrop = cms.bool(True),
     zcut = cms.double(0.1),
     beta = cms.double(0.0),
     R0   = cms.double(0.8),
     useExplicitGhosts = cms.bool(True),
     writeCompound = cms.bool(True),
     jetCollInstanceName=cms.string("SubJets"),
  )

  process.ak8JetsSoftDropMass = cms.EDProducer("RecoJetDeltaRValueMapProducer",
     src = cms.InputTag("ak8Jets"),
     matched = cms.InputTag("ak8JetsSoftDrop"),                                         
     distMax = cms.double(0.8),
     value = cms.string('mass')  
  )

  from RecoJets.JetProducers.ECF_cff import ecfNbeta1
  process.ecfNbeta1 = ecfNbeta1.clone(src = cms.InputTag("ak8Jets"), srcWeights="")

  from RecoJets.JetProducers.nJettinessAdder_cfi import Njettiness
  process.Njettiness = Njettiness.clone(src = cms.InputTag("ak8Jets"), srcWeights="")

  process.ak8ParticleNetJetTagInfos = cms.EDProducer("DeepBoostedJetTagInfoProducer",
      jet_radius = cms.double( 0.8 ),
      min_jet_pt = cms.double( 5.0 ),
      max_jet_eta = cms.double( 2.5 ),
      min_pt_for_track_properties = cms.double( 0.95 ),
      min_pt_for_pfcandidates = cms.double( 0.1 ),
      use_puppiP4 = cms.bool( False ),
      include_neutrals = cms.bool( True ),
      sort_by_sip2dsig = cms.bool( False ),
      min_puppi_wgt = cms.double( -1.0 ),
      flip_ip_sign = cms.bool( False ),
      sip3dSigMax = cms.double( -1.0 ),
      use_hlt_features = cms.bool( False ),
      pf_candidates = cms.InputTag( "pfcands" ),
      jets = cms.InputTag( "ak8Jets" ),
      puppi_value_map = cms.InputTag( "" ),
      use_scouting_features = cms.bool( True ),
      normchi2_value_map = cms.InputTag("pfcands", "normchi2"),
      dz_value_map = cms.InputTag("pfcands", "dz"),
      dxy_value_map = cms.InputTag("pfcands", "dxy"),
      dzsig_value_map = cms.InputTag("pfcands", "dzsig"),
      dxysig_value_map = cms.InputTag("pfcands", "dxysig"),
      lostInnerHits_value_map = cms.InputTag("pfcands", "lostInnerHits"),
      quality_value_map = cms.InputTag("pfcands", "quality"),
      trkPt_value_map = cms.InputTag("pfcands", "trkPt"),
      trkEta_value_map = cms.InputTag("pfcands", "trkEta"),
      trkPhi_value_map = cms.InputTag("pfcands", "trkPhi"),
  )

  from RecoBTag.ONNXRuntime.boostedJetONNXJetTagsProducer_cfi import boostedJetONNXJetTagsProducer
  process.ak8ParticleNetJetTags = cms.EDProducer("BoostedJetONNXValueMapProducer",
      jets = cms.InputTag("ak8Jets"),
      src = cms.InputTag("ak8ParticleNetJetTagInfos"),
      preprocess_json = cms.string("Models/preprocess_doublebtag.json"),
      model_path = cms.FileInPath("Models/doublebtag.onnx"),
      flav_names = cms.vstring(["probHbb", "probHcc","probHqq", "probQCDall"]),
      debugMode = cms.untracked.bool(False),
  )

  process.ak8ParticleNetMassRegressionJetTags = cms.EDProducer("BoostedJetONNXValueMapProducer",
      jets = cms.InputTag("ak8Jets"), 
      src = cms.InputTag("ak8ParticleNetJetTagInfos"),
      preprocess_json = cms.string("Models/preprocess_massreg.json"),
      model_path = cms.FileInPath("Models/massreg.onnx"),
      flav_names = cms.vstring(["mass"]),
      debugMode = cms.untracked.bool(False),
  )

  process.ak8JetTable = cms.EDProducer("SimpleCandidateFlatTableProducer",
      src = cms.InputTag("ak8Jets"),
      name = cms.string("ScoutingFatJet"),
      cut = cms.string(""),
      doc = cms.string("ScoutingFatJet"),
      singleton = cms.bool(False),
      extension = cms.bool(False), # this is the main table
      externalVariables = cms.PSet(
         msoftdrop = ExtVar(cms.InputTag('ak8JetsSoftDropMass'), float, doc="Softdrop mass", precision=10),
         n2b1 = ExtVar(cms.InputTag('ecfNbeta1:ecfN2'), float, doc="N2 with beta=1", precision=10),
         n3b1 = ExtVar(cms.InputTag('ecfNbeta1:ecfN3'), float, doc="N3 with beta=1", precision=10),
         tau1 = ExtVar(cms.InputTag('Njettiness:tau1'), float, doc="Nsubjettiness (1 axis)", precision=10),
         tau2 = ExtVar(cms.InputTag('Njettiness:tau2'), float, doc="Nsubjettiness (2 axis)", precision=10),
         tau3 = ExtVar(cms.InputTag('Njettiness:tau3'), float, doc="Nsubjettiness (3 axis)", precision=10),
         tau4 = ExtVar(cms.InputTag('Njettiness:tau4'), float, doc="Nsubjettiness (4 axis)", precision=10),
         particleNet_mass = ExtVar(cms.InputTag('ak8ParticleNetMassRegressionJetTags:mass'), float, doc="ParticleNet regress mass", precision=10),
         particleNet_prob_Hbb = ExtVar(cms.InputTag('ak8ParticleNetJetTags:probHbb'), float, doc="ParticleNet prob Hbb", precision=10),
         particleNet_prob_Hcc = ExtVar(cms.InputTag('ak8ParticleNetJetTags:probHcc'), float, doc="ParticleNet prob Hcc", precision=10),
         particleNet_prob_Hqq = ExtVar(cms.InputTag('ak8ParticleNetJetTags:probHqq'), float, doc="ParticleNet prob Hqq", precision=10),
         particleNet_prob_QCD = ExtVar(cms.InputTag('ak8ParticleNetJetTags:probQCDall'), float, doc="ParticleNet probbQCD", precision=10),
      ),
      variables = cms.PSet(
         P4Vars,
         area = Var("jetArea()", float, doc="jet catchment area, for JECs",precision=10),
         chHEF = Var("chargedHadronEnergy()/(chargedHadronEnergy()+neutralHadronEnergy()+photonEnergy()+electronEnergy()+muonEnergy())", float, doc="charged Hadron Energy Fraction", precision= 6),
         neHEF = Var("neutralHadronEnergy()/(chargedHadronEnergy()+neutralHadronEnergy()+photonEnergy()+electronEnergy()+muonEnergy())", float, doc="neutral Hadron Energy Fraction", precision= 6),
         chEmEF = Var("(electronEnergy()+muonEnergy())/(chargedHadronEnergy()+neutralHadronEnergy()+photonEnergy()+electronEnergy()+muonEnergy())", float, doc="charged Electromagnetic Energy Fraction", precision= 6),
         neEmEF = Var("(photonEnergy())/(chargedHadronEnergy()+neutralHadronEnergy()+photonEnergy()+electronEnergy()+muonEnergy())", float, doc="neutral Electromagnetic Energy Fraction", precision= 6),
         muEmEF = Var("(muonEnergy())/(chargedHadronEnergy()+neutralHadronEnergy()+photonEnergy()+electronEnergy()+muonEnergy())", float, doc="muon Energy Fraction", precision= 6),
         nCh = Var("chargedHadronMultiplicity()", int, doc="number of charged hadrons in the jet"),
         nNh = Var("neutralHadronMultiplicity()", int, doc="number of neutral hadrons in the jet"),
         nMuons = Var("muonMultiplicity()", int, doc="number of muons in the jet"),
         nElectrons = Var("electronMultiplicity()", int, doc="number of electrons in the jet"),
         nPhotons = Var("photonMultiplicity()", int, doc="number of photons in the jet"),
         nConstituents = Var("numberOfDaughters()", "uint8", doc="Number of particles in the jet")
      ),
  )

  ak8PFJetTaskName = "ak8PFJetTask"
  setattr(process, ak8PFJetTaskName, cms.Task(
      getattr(process,"ak8Jets"),
      getattr(process,"ak8JetsSoftDrop"),
      getattr(process,"ak8JetsSoftDropMass"),
      getattr(process,"ecfNbeta1"),
      getattr(process,"Njettiness"),
      getattr(process,"ak8ParticleNetJetTagInfos"),
      getattr(process,"ak8ParticleNetJetTags"),
      getattr(process,"ak8ParticleNetMassRegressionJetTags"),
      getattr(process,"ak8JetTable"),
    )
  )
  process.customRun3ScoutingNanoAODTask.add(getattr(process,ak8PFJetTaskName))

  return process

def AddAK4MatchToGen(process):
  process.ak4MatchGen = cms.EDProducer("RecoJetToGenJetDeltaRValueMapProducer",
      src = cms.InputTag("ak4Jets"),
      matched = cms.InputTag("slimmedGenJets"),
      distMax = cms.double(0.4),
      value = cms.string("index"),
  )
  ak4MatchGenTaskName = "ak4MatchGenTask"
  setattr(process, ak4MatchGenTaskName, cms.Task(process.ak4MatchGen))
  externalVariables = getattr(process.ak4JetTable, 'externalVariables', cms.PSet())
  externalVariables.genJetIdx = ExtVar(cms.InputTag("ak4MatchGen"), int, doc="gen jet idx")
  process.ak4JetTable.externalVariables = externalVariables
  process.customRun3ScoutingNanoAODTask.add(getattr(process,ak4MatchGenTaskName))
  
  return process

def AddAK8MatchToGen(process):
  process.ak8MatchGen = cms.EDProducer("RecoJetToGenJetDeltaRValueMapProducer",
      src = cms.InputTag("ak8Jets"),
      matched = cms.InputTag("slimmedGenJetsAK8"),
      distMax = cms.double(0.8),
      value = cms.string("index"),
  )
  ak8MatchGenTaskName = "ak8MatchGenTask"
  setattr(process, ak8MatchGenTaskName, cms.Task(process.ak8MatchGen))
  externalVariables = getattr(process.ak8JetTable, 'externalVariables', cms.PSet())
  externalVariables.genJetAK8Idx = ExtVar(cms.InputTag("ak8MatchGen"), int, doc="gen jet idx")
  process.ak8JetTable.externalVariables = externalVariables
  process.customRun3ScoutingNanoAODTask.add(getattr(process,ak8MatchGenTaskName))
  
  return process

#===========================================================================
#
# CUSTOMIZATION function
#
#===========================================================================
def PrepRun3ScoutingCustomNanoAOD(process,runOnMC):
  process.customRun3ScoutingNanoAODTask = cms.Task()
  process = ConvertScoutingToReco(process)
  #process = AddRun3Scouting(process)
  #process = AddParticles(process)
  #process = AddAK4PFJets(process)
  #process = AddAK8PFJets(process)

  if (runOnMC):
    process = AddAK4MatchToGen(process)
    process = AddAK8MatchToGen(process)

  process.schedule.associate(process.customRun3ScoutingNanoAODTask)

  return process

def PrepRun3ScoutingCustomNanoAOD_MC(process):
  process = PrepRun3ScoutingCustomNanoAOD(process,runOnMC=True)

  return process

def PrepRun3ScoutingCustomNanoAOD_Data(process):
  process = PrepRun3ScoutingCustomNanoAOD(process,runOnMC=False)

  return process
