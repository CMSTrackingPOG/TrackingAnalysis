import FWCore.ParameterSet.Config as cms

readFiles = cms.untracked.vstring()
secFiles = cms.untracked.vstring() 
source = cms.Source ("PoolSource",fileNames = readFiles, secondaryFileNames = secFiles)
readFiles.extend( [  #before run introducing prescaling for jetTriggers (135445)
#'/store/relval/CMSSW_10_6_0/RelValQCD_FlatPt_15_3000HS_13UP17/FEVTDEBUGHLT/PUpmx25ns_106X_mc2017_realistic_v3_ulhlt17hs_pmx-v1/10000/ECF91A89-E170-1E4A-827B-D7E06212D1A6.root'
#'/store/relval/CMSSW_10_6_0/RelValMinBias_13/GEN-SIM-RECO/106X_mcRun2_asymptotic_v3-v1/10000/DD15408C-2088-4F47-822F-F1D5AED95767.root'
#'/store/relval/CMSSW_10_6_0/ZeroBias/RECO/106X_dataRun2_relval_v9_RelVal_2016B-v1/10000/0BCB036A-1778-3C4E-8683-62185F9C8E95.root'
'/store/data/Run2016B/ZeroBias/AOD/07Aug17_ver1-v1/50000/D2F70B66-BD8B-E711-8761-E0071B7BC1A1.root'
]);

# /ZeroBias/CMSSW_10_6_0-106X_dataRun2_relval_v9_RelVal_2016B-v1/RECO
# /RelValMinBias_13/CMSSW_10_6_0-106X_mcRun2_asymptotic_v3-v1/GEN-SIM-RECO

# /RelValTTbar_13/CMSSW_11_0_0_pre10-110X_mcRun2_asymptotic_v2-v1/GEN-SIM-RECO 
# /RelValTTbar_13/CMSSW_11_0_0_pre10-PU25ns_110X_mcRun2_asymptotic_v2-v1/GEN-SIM-RECO

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

process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
##process.GlobalTag.globaltag = '110X_mcRun2_asymptotic_v2'
process.GlobalTag.globaltag = '80X_dataRun2_2016LegacyRepro_v4'

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

process.load('IpResoStudies.EDAnalyzers.residuals_cfi')
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
