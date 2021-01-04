from CRABClient.UserUtilities import config
config = config()

config.General.requestName = 'SoftDisplacedLeptonHLT_crabRunSkim2'
config.General.workArea = 'crab_SoftDisplacedLeptonHLT'
config.General.transferOutputs = True
config.General.instance = 'prod'

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'MyHLT_HLT.py'
config.JobType.inputFiles = ['json_2018D_Ephemeral_20181022.json']

config.Data.inputDataset = '/EphemeralHLTPhysics1/Run2018D-v1/RAW'
config.Data.inputDBS = 'global'
config.Data.splitting = 'LumiBased'
config.Data.unitsPerJob = 20
config.Data.lumiMask = 'json_2018D_Ephemeral_20181022.json'
config.Data.publication = True
config.Data.outputDatasetTag = 'DisplacedLeptonHLT_EphemeralHLTSkim2_asahasra'
config.Data.ignoreLocality = True

config.Site.whitelist = ['T2_US*','T2_CH*']
config.Site.storageSite = 'T2_BE_IIHE'
