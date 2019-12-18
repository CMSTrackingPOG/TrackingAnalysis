import FWCore.ParameterSet.Config as cms

readFiles = cms.untracked.vstring()
secFiles = cms.untracked.vstring() 
source = cms.Source ("PoolSource",fileNames = readFiles, secondaryFileNames = secFiles)
readFiles.extend( [
#'file:/afs/cern.ch/work/k/kskovpen/GEN-SIM-DIGI-RAW-HLTDEBUG.root'
'/store/relval/CMSSW_10_6_0/RelValQCD_FlatPt_15_3000HS_13UP17/FEVTDEBUGHLT/PUpmx25ns_106X_mc2017_realistic_v3_ulhlt17hs_pmx-v1/10000/ECF91A89-E170-1E4A-827B-D7E06212D1A6.root'
#'/store/relval/CMSSW_10_6_0/ZeroBias/FEVTDEBUGHLT/103X_dataRun2_HLT_relval_v8_RelVal_2016B-v1/10000/F6017F7C-761F-1044-ABA7-80A710D7CFAD.root'
]);

# /RelValTTbar_13/CMSSW_11_0_0_pre10-PU25ns_110X_mcRun2_asymptotic_v2-v1/GEN-SIM-RECO
# /RelValTTbar_13/CMSSW_11_0_0_pre10-110X_mcRun2_asymptotic_v2_FastSim-v1/GEN-SIM-DIGI-RECO
# /MinBias_TuneCUETP8M1_13TeV-pythia8/RunIIFall15DR76-PU25nsData2015v1Recodebug_76X_mcRun2_asymptotic_RunIIFall15DR76_v1-v1/GEN-SIM-RECODEBUG
# /ZeroBias/Run2015D-PromptReco-v4/RECO

# /RelValQCD_FlatPt_15_3000HS_13UP17/CMSSW_10_6_0-PUpmx25ns_106X_mc2017_realistic_v3_ulhlt17hs_pmx-v1/FEVTDEBUGHLT
# /ZeroBias/CMSSW_10_6_0-103X_dataRun2_HLT_relval_v8_RelVal_2016B-v1/FEVTDEBUGHLT

# go to https://twiki.cern.ch/twiki/bin/view/CMS/PdmVRelVals2019
# run runTheMatrix.py command to get the info on the cmsDriver.py
# runTheMatrix.py -n
# runTheMatrix.py --what pileup -l 25202 --dryRun

secFiles.extend( [ ]);

process = cms.Process("IpResolutions")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 10000
# import of standard configurations
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load("Configuration.StandardSequences.MagneticField_cff")
#process.load('Configuration.StandardSequences.RawToDigi_cff')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')

#from SimGeneral.MixingModule.pixelDigitizer_cfi import *

process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.GlobalTag.globaltag = '110X_mcRun2_asymptotic_v2'

##process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
##process.GlobalTag.globaltag = 'SET_GLOBALTAG'

#process.load("SimTracker.Configuration.SimTracker_cff")

##process.load("SimTracker.TrackAssociation.trackingParticleRecoTrackAsssociation_cfi")
process.load("SimGeneral.TrackingAnalysis.simHitTPAssociation_cfi")
##process.simHitTPAssocProducer.simHitSrc = cms.VInputTag(cms.InputTag("fastSimProducer","TrackerHits"), cms.InputTag("MuonSimHits","MuonCSCHits"), cms.InputTag("MuonSimHits","MuonDTHits"), cms.InputTag("MuonSimHits","MuonRPCHits"))
##process.load("Validation.RecoTrack.MultiTrackValidator_cff")
##process.load("Validation.RecoTrack.TrackValidation_cff")

#process.RandomNumberGeneratorService = cms.Service("RandomNumberGeneratorService",
#simSiPixelDigis = cms.PSet(initialSeed = cms.untracked.uint32(1234567),
#                           engineName = cms.untracked.string('TRandom3')),
#SiPixelDigitizer = cms.PSet(initialSeed = cms.untracked.uint32(1234567),
#                           engineName = cms.untracked.string('TRandom3')),
#mix = cms.PSet(initialSeed = cms.untracked.uint32(1234568),
#               engineName = cms.untracked.string('TRandom3'))                                     
#)

#process.RandomNumberGeneratorService = cms.Service("RandomNumberGeneratorService",
#  g4SimHits = cms.PSet(
#  initialSeed = cms.untracked.uint32(123456789),
#  engineName = cms.untracked.string('TRandom3')
#  ),
#  LHCTransport = cms.PSet(
#  initialSeed = cms.untracked.uint32(321456789),
#  engineName = cms.untracked.string('TRandom3')
#  )
#)

process.load("SimGeneral.MixingModule.mixNoPU_cfi")
process.load("SimGeneral.MixingModule.trackingTruthProducerSelection_cfi")
process.trackingParticles.simHitCollections = cms.PSet( )
process.mix.playback = cms.untracked.bool(False)
process.mix.digitizers = cms.PSet(
  mergedtruth = cms.PSet(process.trackingParticles)
)
for a in process.aliases: delattr(process, a)

process.load("SimTracker.TrackAssociation.LhcParametersDefinerForTP_cfi")

#process.load("SimTracker.TrackerHitAssociation.tpClusterProducer_cfi")
#process.tpClusterProducer.pixelSimLinkSrc = cms.InputTag("mix", "simSiPixelDigis")

process.load("SimTracker.TrackAssociatorProducers.trackAssociatorByHits_cfi")
#process.load("SimTracker.TrackAssociation.trackingParticleRecoTrackAsssociation_cfi")

#process.load("SimTracker.TrackAssociatorProducers.quickTrackAssociatorByHits_cfi")
#process.quickTrackAssociatorByHits.SimToRecoDenominator = cms.string('reco')
#process.quickTrackAssociatorByHits.useClusterTPAssociation = cms.bool(True)

import DQMServices.Core.DQMStore_cfi as DQM
DQM.DQMStore.collateHistograms = False

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100) )

process.source = source

process.load('IpResoStudies.EDAnalyzers.vertexResponsesAndTrueResolutions_cfi')

process.TFileService = cms.Service("TFileService", 
      fileName = cms.string("output.root"),
      closeFileFast = cms.untracked.bool(True)
)

process.HLTMinBias = cms.EDFilter("HLTHighLevel",
                                  TriggerResultsTag = cms.InputTag("TriggerResults","","HLT"),
                                  HLTPaths = cms.vstring('HLT_AK4PFJet30_v*'),
#                                  HLTPaths = cms.vstring('HLT_ZeroBias_v*'),
                                  eventSetupPathsKey = cms.string(''), # not empty => use read paths from AlCaRecoTriggerBitsRcd via this key
                                  andOr = cms.bool(True),              # how to deal with multiple triggers: True (OR) accept if ANY is true, False (AND) accept if ALL are true
                                  throw = cms.bool(True)               # throw exception on unknown path names
                                  )

process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(False) )

process.p = cms.Path(
#process.HLTMinBias*
#process.mix*
process.simHitTPAssocProducer*
#process.tpClusterProducer*
#process.quickTrackAssociatorByHits*
process.trackAssociatorByHits*
process.vertexResponsesAndTrueResolutions
)
