from WMCore.Configuration import Configuration
config = Configuration()

config.section_('General')
config.General.requestName = 'REQUESTNAME'
config.General.transferLogs = True
config.section_('JobType')
config.JobType.psetName = '../test/residuals.py'
config.JobType.pluginName = 'Analysis'
config.JobType.pyCfgParams = ['withBS=1']
#config.JobType.allowUndistributedCMSSW = True
config.JobType.maxMemoryMB = 4000

config.section_('Data')
config.Data.splitting='FileBased'
#config.Data.splitting='EventAwareLumiBased'
config.Data.totalUnits = -1
config.Data.unitsPerJob = 1
#config.Data.unitsPerJob = 5000

#config.Data.allowNonValidInputDataset = True
config.Data.publication = False
config.Data.inputDataset = 'INPUTDATASET'
config.Data.outputDatasetTag = 'PUBLISHDATANAME'
config.Data.publishDBS = 'https://cmsweb.cern.ch/dbs/prod/phys03/DBSWriter'
config.Data.outLFNDirBase = 'OUTLFN'

config.section_('User')
config.User.voGroup = 'becms'
config.section_('Site')
config.Site.storageSite = 'T2_BE_IIHE'
config.Site.whitelist = ['T2_BE_IIHE']
