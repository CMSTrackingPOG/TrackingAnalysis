import FWCore.ParameterSet.Config as cms

readFiles = cms.untracked.vstring()
secFiles = cms.untracked.vstring() 
source = cms.Source ("PoolSource",fileNames = readFiles, secondaryFileNames = secFiles)
readFiles.extend( [
'/store/data/Run2016B/ZeroBias/AOD/07Aug17_ver1-v1/50000/D2F70B66-BD8B-E711-8761-E0071B7BC1A1.root'
]);

secFiles.extend([ ]);

process = cms.Process("IpResiduals")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 5000

# import of standard configurations
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.RawToDigi_cff')
process.load('Configuration.StandardSequences.Reconstruction_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.EventContent.EventContent_cff')

#process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
#from Configuration.AlCa.autoCond import autoCond
#process.GlobalTag.globaltag = autoCond['com10_2011']

process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.GlobalTag.globaltag = 'FT_R_53_LV5::All'

##process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1000) ) # Data 
##process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(10000) ) # MC noPU

process.source = source

process.HLTMinBias = cms.EDFilter("HLTHighLevel",
                                  TriggerResultsTag = cms.InputTag("TriggerResults","","HLT"),
                                  #HLTPaths = cms.vstring('HLT_L1_BscMinBiasOR_BptxPlusORMinus'),
                                  HLTPaths = cms.vstring('HLT_*'),
                                  eventSetupPathsKey = cms.string(''), # not empty => use read paths from AlCaRecoTriggerBitsRcd via this key
                                  andOr = cms.bool(True),              # how to deal with multiple triggers: True (OR) accept if ANY is true, False (AND) accept if ALL are true
                                  throw = cms.bool(True)               # throw exception on unknown path names
                                  )


# process.load("RecoTracker.TrackProducer.TrackRefitters_cff")
# # remove the following lines if you run on RECO files
# process.TrackRefitter.src = 'generalTracks'
# process.TrackRefitter.NavigationSchool = ''

## PV refit
# process.load("TrackingTools.TransientTrack.TransientTrackBuilder_cfi")
# from RecoVertex.PrimaryVertexProducer.OfflinePrimaryVertices_cfi import offlinePrimaryVertices
# process.offlinePrimaryVerticesFromRefittedTrks  = offlinePrimaryVertices.clone()
# process.offlinePrimaryVerticesFromRefittedTrks.TrackLabel                                       = cms.InputTag("TrackRefitter")
# process.offlinePrimaryVerticesFromRefittedTrks.vertexCollections.maxDistanceToBeam              = 1
# process.offlinePrimaryVerticesFromRefittedTrks.TkFilterParameters.maxNormalizedChi2             = 20
# process.offlinePrimaryVerticesFromRefittedTrks.TkFilterParameters.minSiliconLayersWithHits      = 5
# process.offlinePrimaryVerticesFromRefittedTrks.TkFilterParameters.maxD0Significance             = 5.0
# process.offlinePrimaryVerticesFromRefittedTrks.TkFilterParameters.minPixelLayersWithHits        = 2

process.load('TrackingAnalysis.EDAnalyzers.residuals_cfi')
#process.residuals.TrackLabel = cms.InputTag("TrackRefitter")
#process.residuals.VertexLabel = cms.InputTag("offlinePrimaryVerticesFromRefittedTrks")

process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string("output.root"),
                                   closeFileFast = cms.untracked.bool(True)
                                   )

process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(False) )

process.p = cms.Path(#process.HLTMinBias*
    # process.offlineBeamSpot                        +
    # process.TrackRefitter                          +
    # process.offlinePrimaryVerticesFromRefittedTrks +
    process.residuals
)
