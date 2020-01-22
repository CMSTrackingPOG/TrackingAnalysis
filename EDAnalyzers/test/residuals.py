import FWCore.ParameterSet.Config as cms
from FWCore.ParameterSet.VarParsing import VarParsing

options = VarParsing('analysis')
options.register('withBS',False,VarParsing.multiplicity.singleton,VarParsing.varType.bool,'Primary vertex reconstruction with BS constraint')
options.parseArguments()

readFiles = cms.untracked.vstring()
secFiles = cms.untracked.vstring() 
source = cms.Source ("PoolSource",fileNames = readFiles, secondaryFileNames = secFiles)
readFiles.extend( [
#'/store/data/Run2011A/MinimumBias/AOD/12Oct2013-v1/20001/6A9F67BD-8746-E311-AC78-0025901D5CDC.root'
'root://maite.iihe.ac.be/pnfs/iihe/cms/ph/sc4/store/data/Run2012A/MinimumBias/AOD/22Jan2013-v1/20000/1097F002-8D67-E211-8536-003048FFCC18.root'
#'/store/data/Run2012A/MinimumBias/AOD/22Jan2013-v1/20003/F286D9D5-8667-E211-949A-003048679188.root'
#'/store/mc/Fall11/MinBias_Tune4C_7TeV-pythia8/AODSIM/PU_S6_START44_V9B-v1/0001/FEAD5287-9736-E111-AF86-0030487D8121.root'
#'/store/mc/Summer12_DR53X/MinBias_Tune4C_8TeV-pythia8/AODSIM/PU_S10_START53_V7A-v1/0000/1EB5119C-52DA-E111-9679-0018F3D09642.root'
]);

secFiles.extend([ ]);

process = cms.Process("IpResiduals")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport = cms.untracked.PSet( reportEvery = cms.untracked.int32(10) )
#process.MessageLogger.cerr.FwkReport = cms.untracked.PSet( reportEvery = cms.untracked.int32(1) )
                
#process.MessageLogger.cout.placeholder = cms.untracked.bool(False)
#process.MessageLogger.cout.threshold = cms.untracked.string('INFO')
#process.MessageLogger.cout.default = cms.untracked.PSet( limit = cms.untracked.int32(10000000) )
#process.MessageLogger.cout.FwkReport = cms.untracked.PSet( reportEvery = cms.untracked.int32(1000) )
#process.MessageLogger.debugModules = cms.untracked.vstring('*')

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
#process.GlobalTag.globaltag = 'FT_R_53_LV5::All'
#FT_R_53_V18::All # 2012A-C
#FT_R_53_V21::All # 2012D
#process.GlobalTag.globaltag = 'START44_V9B::All'
#process.GlobalTag.globaltag = 'START53_V7A::All'
process.GlobalTag.globaltag = 'FT_R_53_V21::All'

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
#process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100) )

process.source = source

process.HLT = cms.EDFilter("HLTHighLevel",
                           TriggerResultsTag = cms.InputTag("TriggerResults","","HLT"),
                           HLTPaths = cms.vstring('HLT_ZeroBiasPixel_DoubleTrack_v*','HLT_ZeroBias_v*'),
#                           HLTPaths = cms.vstring('HLT_*'),
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

process.load('RecoVertex.PrimaryVertexProducer.OfflinePrimaryVertices_cfi')
#process.load('RecoVertex.PrimaryVertexProducer.OfflinePrimaryVerticesWithBS_cfi')

primVtx = process.offlinePrimaryVertices
PVSelParameters = cms.PSet( maxDistanceToBeam = primVtx.vertexCollections[0].maxDistanceToBeam )
process.offlinePrimaryVerticesRerun = primVtx.clone( PVSelParameters = PVSelParameters,
                                                     useBeamConstraint = primVtx.vertexCollections[0].useBeamConstraint,
                                                     verbose = False,
                                                     algorithm = primVtx.vertexCollections[0].algorithm,
                                                     minNdof = primVtx.vertexCollections[0].minNdof
)

if options.withBS:
    PVSelParameters = cms.PSet( maxDistanceToBeam = primVtx.vertexCollections[1].maxDistanceToBeam )
    process.offlinePrimaryVerticesRerun = primVtx.clone( PVSelParameters = PVSelParameters,
                                                         useBeamConstraint = primVtx.vertexCollections[1].useBeamConstraint,
                                                         verbose = False,
                                                         algorithm = primVtx.vertexCollections[1].algorithm,
                                                         minNdof = primVtx.vertexCollections[1].minNdof
    )

process.offlinePrimaryVerticesRerun.TkClusParameters.algorithm = cms.string("DA")
#print process.offlinePrimaryVerticesRerun.dumpPython()

process.load('TrackingAnalysis.EDAnalyzers.residuals_cfi')

if options.withBS:
    process.residuals.VertexPrimaryLabel = cms.InputTag('offlinePrimaryVerticesWithBS')

#process.residuals.TrackLabel = cms.InputTag("TrackRefitter")
#process.residuals.VertexLabel = cms.InputTag("offlinePrimaryVerticesFromRefittedTrks")

process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string("output.root"),
                                   closeFileFast = cms.untracked.bool(True)
                                   )

process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(False) )

process.p = cms.Path(process.HLT*
    # process.offlineBeamSpot                        +
    # process.TrackRefitter                          +
    # process.offlinePrimaryVerticesFromRefittedTrks +
    process.offlinePrimaryVerticesRerun*
    process.residuals
)
