import FWCore.ParameterSet.Config as cms
from  PhysicsTools.NanoAOD.common_cff import *
from PhysicsTools.NanoAOD.simpleCandidateFlatTableProducer_cfi import simpleCandidateFlatTableProducer

process.run3ScoutingTask = cms.Task()

def AddRun3Scouting(proc):
  """
  Add all Run 3 scouting objects except jets and PFCandidates
  """

  setattr(proc,"photonScoutingTable", simpleCandidateFlatTableProducer.clone(
        src = cms.InputTag("hltScoutingEgammaPacker"),
        cut = cms.string(""),
        name = cms.string("ScoutingPhoton"),
        doc  = cms.string("Photon Scouting Informations"),
        singleton = cms.bool(False), # the number of entries is variable
        extension = cms.bool(False), # this is the main table for the muons
        variables = cms.PSet(
            pt = Var('pt', 'float', precision=14, doc='pt'),
            eta = Var('eta', 'float', precision=14, doc='eta'),
            phi = Var('phi', 'float', precision=14, doc='phi'),
            m = Var('m', 'float', precision=14, doc='m'),
            sigmaIetaIeta = Var('sigmaIetaIeta', 'float', precision=14, doc='sigmaIetaIeta'),
            hOverE = Var('hOverE', 'float', precision=14, doc='hOverE'),
            ecalIso = Var('ecalIso', 'float', precision=14, doc='ecalIso'),
            hcalIso = Var('hcalIso', 'float', precision=14, doc='hcalIso'),
            trackIso = Var('trkIso', 'float', precision=14, doc='trackIso'),
            r9 = Var('r9', 'float', precision=14, doc='r9'),
            sMin = Var('sMin', 'float', precision=14, doc='sMin'),
            sMaj = Var('sMaj', 'float', precision=14, doc='sMaj'),
            seedId = Var('seedId', 'int', doc='seedId'),
        )
     )
  )

  process.run3ScoutingTask.add(cms.Task(getattr(proc,"photonScoutingTable")))

#===========================================================================
#
# CUSTOMIZATION function
#
#===========================================================================
def PrepRun3ScoutingCustomNanoAOD(process,runOnMC):

  process = AddRun3Scouting(process)
  process.schedule.associate(process.run3ScoutingTask)

def PrepRun3ScoutingCustomNanoAOD_MC(process):
  process = PrepRun3ScoutingCustomNanoAOD(process,runOnMC=True)

  return process

def PrepRun3ScoutingCustomNanoAOD_Data(process):
  process = PrepRun3ScoutingCustomNanoAOD(process,runOnMC=False)
  return process
