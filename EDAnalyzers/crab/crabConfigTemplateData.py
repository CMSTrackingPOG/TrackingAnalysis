from WMCore.Configuration import Configuration
config = Configuration()

config.section_('General')
config.General.requestName = 'REQUESTNAME'
config.General.transferLogs = True
config.section_('JobType')
config.JobType.psetName = '../test/residuals.py'
config.JobType.pluginName = 'Analysis'
config.JobType.pyCfgParams = ['isData=1', 'withBS=0', 'is2017=1']
config.JobType.allowUndistributedCMSSW = True
#config.JobType.maxMemoryMB = 4000
##config.JobType.maxJobRuntimeMin = 2749 # min

config.section_('Data')
config.Data.splitting='LumiBased'
config.Data.totalUnits = -1
config.Data.unitsPerJob = 50

#config.Data.allowNonValidInputDataset = True
#config.Data.lumiMask = 'JSON/Cert_271036-284044_13TeV_Legacy2016_Collisions16_JSON.txt'
config.Data.lumiMask = 'JSON/Cert_294927-306462_13TeV_UL2017_Collisions17_GoldenJSON.txt'
#config.Data.lumiMask = 'JSON/Cert_314472-325175_13TeV_Legacy2018_Collisions18_JSON.txt'
config.Data.publication = False
config.Data.inputDataset = 'INPUTDATASET'
config.Data.outputDatasetTag = 'PUBLISHDATANAME'
config.Data.publishDBS = 'https://cmsweb.cern.ch/dbs/prod/phys03/DBSWriter'
config.Data.outLFNDirBase = 'OUTLFN'

config.section_('User')
config.User.voGroup = 'becms'
config.section_('Site')
config.Site.storageSite = 'T2_BE_IIHE'
#config.Site.whitelist = ['T2_BE_IIHE']
