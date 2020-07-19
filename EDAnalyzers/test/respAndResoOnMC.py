import FWCore.ParameterSet.Config as cms

readFiles = cms.untracked.vstring()
secFiles = cms.untracked.vstring() 
source = cms.Source ("PoolSource",fileNames = readFiles, secondaryFileNames = secFiles)
readFiles.extend( [
'/store/mc/RunIISummer19UL17RECO/SingleNeutrino/AODSIM/FEVTDEBUG_106X_mc2017_realistic_v6-v2/20000/880C9813-73EA-574A-9208-A0DC4B177C13.root'
]);

secFiles.extend( [
'/store/mc/RunIISummer19UL17DIGI/SingleNeutrino/GEN-SIM-DIGI-RAW/FEVTDEBUG_106X_mc2017_realistic_v6-v2/210002/19281869-E1AD-9C41-986F-7D87A5A94602.root'
] );

process = cms.Process("IpResolutions")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 10000

process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')

process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.GlobalTag.globaltag = '106X_mc2017_realistic_v6'

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

process.load('TrackingAnalysis.EDAnalyzers.vertexResponsesAndTrueResolutions_cfi')

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
