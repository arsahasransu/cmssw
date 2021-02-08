# Auto generated configuration file
# using: 
# Revision: 1.19 
# Source: /local/reps/CMSSW/CMSSW/Configuration/Applications/python/ConfigBuilder.py,v 
# with command line options: HLT2 --step=HLT:DisplacedLeptonMuMuHT_Trial5 --era=Run2_2018 --data --conditions auto:run3_hlt_GRun --filein root://xrootd-cms.infn.it///store/data/Run2018D/EphemeralHLTPhysics1/RAW/v1/000/320/617/00000/0E0A4B9E-D994-E811-AA44-02163E017CDA.root --processName=HLT2 -n 100 --no_exec
import FWCore.ParameterSet.Config as cms

from Configuration.Eras.Era_Run2_2018_cff import Run2_2018

process = cms.Process('HLT2',Run2_2018)

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff')
process.load('HLTrigger.Configuration.HLT_DisplacedLeptonMuMuHT_Trial5_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1),
    output = cms.optional.untracked.allowed(cms.int32,cms.PSet)
)
  









# Input source
process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring('file:/eos/cms/store/data/Run2018D/EphemeralHLTPhysics8/RAW/v1/000/323/775/00000/043E7677-75F6-2F40-824F-3E189D052B7F.root',
                                                              'file:/eos/cms/store/data/Run2018D/EphemeralHLTPhysics8/RAW/v1/000/323/775/00000/0827FBAA-6E09-2549-94CD-A061A5313127.root',
                                                              'file:/eos/cms/store/data/Run2018D/EphemeralHLTPhysics8/RAW/v1/000/323/775/00000/54429C76-5510-BF48-842C-6C4B7A7CE2D7.root',
                                                              'file:/eos/cms/store/data/Run2018D/EphemeralHLTPhysics8/RAW/v1/000/323/775/00000/911883B1-A4BD-3240-A162-CA9D6B3B8EA6.root',
                                                              'file:/eos/cms/store/data/Run2018D/EphemeralHLTPhysics8/RAW/v1/000/323/775/00000/E213E146-4910-C744-A7A3-81C49E3DD2BF.root',
                                                              'file:/eos/cms/store/data/Run2018D/EphemeralHLTPhysics8/RAW/v1/000/323/775/00000/0E695C4C-7D5C-C641-AE82-3E5901DF846F.root',
                                                              'file:/eos/cms/store/data/Run2018D/EphemeralHLTPhysics8/RAW/v1/000/323/775/00000/65F52FA8-D922-8D4A-9959-4205515C5F13.root',
                                                              'file:/eos/cms/store/data/Run2018D/EphemeralHLTPhysics8/RAW/v1/000/323/775/00000/9A7F5493-40D6-604C-9971-EBC19326FF62.root',
                                                              'file:/eos/cms/store/data/Run2018D/EphemeralHLTPhysics8/RAW/v1/000/323/775/00000/E7792316-8C19-6E4F-9BCB-6FAB57A5B603.root',
                                                              'file:/eos/cms/store/data/Run2018D/EphemeralHLTPhysics8/RAW/v1/000/323/775/00000/0EC7851C-BBE9-8C46-89C8-A08514945D3A.root',
                                                              'file:/eos/cms/store/data/Run2018D/EphemeralHLTPhysics8/RAW/v1/000/323/775/00000/67EFF9D9-7789-5840-991C-94440CAF803C.root'),
    secondaryFileNames = cms.untracked.vstring()
)

process.options = cms.untracked.PSet(
    FailPath = cms.untracked.vstring(),
    IgnoreCompletely = cms.untracked.vstring(),
    Rethrow = cms.untracked.vstring(),
    SkipEvent = cms.untracked.vstring(),
    allowUnscheduled = cms.obsolete.untracked.bool,
    canDeleteEarly = cms.untracked.vstring(),
    emptyRunLumiMode = cms.obsolete.untracked.string,
    eventSetup = cms.untracked.PSet(
        forceNumberOfConcurrentIOVs = cms.untracked.PSet(

        ),
        numberOfConcurrentIOVs = cms.untracked.uint32(1)
    ),
    fileMode = cms.untracked.string('FULLMERGE'),
    forceEventSetupCacheClearOnNewRun = cms.untracked.bool(False),
    makeTriggerResults = cms.obsolete.untracked.bool,
    numberOfConcurrentLuminosityBlocks = cms.untracked.uint32(1),
    numberOfConcurrentRuns = cms.untracked.uint32(1),
    numberOfStreams = cms.untracked.uint32(0),
    numberOfThreads = cms.untracked.uint32(1),
    printDependencies = cms.untracked.bool(False),
    sizeOfStackForThreadsInKB = cms.optional.untracked.uint32,
    throwIfIllegalParameter = cms.untracked.bool(True),
    wantSummary = cms.untracked.bool(False)
)

# Production Info
process.configurationMetadata = cms.untracked.PSet(
    annotation = cms.untracked.string('HLT2 nevts:100'),
    name = cms.untracked.string('Applications'),
    version = cms.untracked.string('$Revision: 1.19 $')
)

# Output definition

process.RECOSIMoutput = cms.OutputModule("PoolOutputModule",
    dataset = cms.untracked.PSet(
        dataTier = cms.untracked.string(''),
        filterName = cms.untracked.string('')
    ),
    fileName = cms.untracked.string('HLT2_HLT.root'),
    outputCommands = process.RECOSIMEventContent.outputCommands,
    splitLevel = cms.untracked.int32(0)
)
process.RECOSIMoutput.outputCommands.append('keep *_hltL3NoFiltersNoVtxMuonCandidates_*_*')  
process.RECOSIMoutput.outputCommands.append('keep *_hltL3NoFiltersNoVtxMuons_*_*')
process.RECOSIMoutput.outputCommands.append('keep *_*BeamSpot_*_*')

# Additional output definition

# Other statements
from HLTrigger.Configuration.CustomConfigs import ProcessName
process = ProcessName(process)

from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run3_hlt_GRun', '')

# Path and EndPath definitions
process.endjob_step = cms.EndPath(process.endOfProcess)
process.RECOSIMoutput_step = cms.EndPath(process.RECOSIMoutput)

#process.demo = cms.EDAnalyzer('MyTriggerAnalyzerRAW')
#process.TFileService = cms.Service("TFileService",
#                                   fileName = cms.string( "DisplacedLepton_MuMuHT.root" )
#                               )
#process.demo_step = cms.EndPath(process.demo)

# Schedule definition
process.schedule = cms.Schedule()
process.schedule.extend(process.HLTSchedule)
process.schedule.extend([process.endjob_step,process.RECOSIMoutput_step])
from PhysicsTools.PatAlgos.tools.helpers import associatePatAlgosToolsTask
associatePatAlgosToolsTask(process)


# Customisation from command line

# Add early deletion of temporary data products to reduce peak memory need
from Configuration.StandardSequences.earlyDeleteSettings_cff import customiseEarlyDelete
process = customiseEarlyDelete(process)
# End adding early deletion
