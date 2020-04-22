import FWCore.ParameterSet.Config as cms
from FWCore.ParameterSet.VarParsing import VarParsing

options = VarParsing('analysis')
options.register('withBS',False,VarParsing.multiplicity.singleton,VarParsing.varType.bool,'Primary vertex reconstruction with BS constraint')
options.register('isData',False,VarParsing.multiplicity.singleton,VarParsing.varType.bool,'Run on data')
options.parseArguments()

readFiles = cms.untracked.vstring()
secFiles = cms.untracked.vstring() 
source = cms.Source ("PoolSource",fileNames = readFiles, secondaryFileNames = secFiles)
readFiles.extend( [
'/store/mc/RunIISummer19UL17RECO/SingleNeutrino/AODSIM/106X_mc2017_realistic_v6-v2/30000/01989700-04AA-5F49-BDFE-C79CEA5EECE3.root'
#'/store/data/Run2017F/ZeroBias/AOD/09Aug2019_UL2017-v1/00000/E543BEF4-7BDC-F141-87EC-734B88DC42D6.root'
#'/store/mc/RunIIAutumn18DRPremix/SingleNeutrino/AODSIM/forRECO_102X_upgrade2018_realistic_v15_ext1-v1/100000/03EA430F-EE65-EE48-8435-A9F4A1E08D43.root'
#'/store/data/Run2018C/ZeroBias/AOD/17Sep2018-v1/110000/0773B928-0421-A44D-A3F9-B09C054FBE63.root'
]);

secFiles.extend([ ]);

process = cms.Process("IpResiduals")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport = cms.untracked.PSet( reportEvery = cms.untracked.int32(1) )

process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
#process.GlobalTag.globaltag = '106X_dataRun2_v20'
process.GlobalTag.globaltag = '106X_mc2017_realistic_v6'
#process.GlobalTag.globaltag = '102X_dataRun2_v12'
#process.GlobalTag.globaltag = '102X_upgrade2018_realistic_v20'

process.load("CondCore.CondDB.CondDB_cfi")
process.load('Configuration.Geometry.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
#process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(10) )

process.source = source

#process.HLT = cms.EDFilter("HLTHighLevel",
#                           TriggerResultsTag = cms.InputTag("TriggerResults","","HLT"),
#                           HLTPaths = cms.vstring('HLT_ZeroBiasPixel_DoubleTrack_v*','HLT_ZeroBias_v*'),
#                           eventSetupPathsKey = cms.string(''), # not empty => use read paths from AlCaRecoTriggerBitsRcd via this key
#                           andOr = cms.bool(True),              # how to deal with multiple triggers: True (OR) accept if ANY is true, False (AND) accept if ALL are true
#                           throw = cms.bool(True)               # throw exception on unknown path names
#)

process.load("TrackingTools.TransientTrack.TransientTrackBuilder_cfi")
process.load('RecoVertex.PrimaryVertexProducer.OfflinePrimaryVertices_cfi')

primVtx = process.offlinePrimaryVertices
process.offlinePrimaryVerticesRerun = primVtx.clone( maxDistanceToBeam = primVtx.vertexCollections[0].maxDistanceToBeam,
                                                     useBeamConstraint = primVtx.vertexCollections[0].useBeamConstraint,
                                                     chi2cutoff = primVtx.vertexCollections[0].chi2cutoff,
                                                     verbose = False,
                                                     algorithm = primVtx.vertexCollections[0].algorithm,
                                                     minNdof = primVtx.vertexCollections[0].minNdof,
                                                     label = primVtx.vertexCollections[0].label
)

if options.withBS:
    process.offlinePrimaryVerticesRerun = primVtx.clone( maxDistanceToBeam = primVtx.vertexCollections[1].maxDistanceToBeam,
                                                         useBeamConstraint = primVtx.vertexCollections[1].useBeamConstraint,
                                                         chi2cutoff = primVtx.vertexCollections[1].chi2cutoff,
                                                         verbose = False,
                                                         algorithm = primVtx.vertexCollections[1].algorithm,
                                                         minNdof = primVtx.vertexCollections[1].minNdof,
                                                         label = primVtx.vertexCollections[1].label
    )

#print process.offlinePrimaryVerticesRerun.dumpPython()

process.load('TrackingAnalysis.EDAnalyzers.residuals_cfi')
process.residuals.BeamSpotConfig = ''

if options.withBS:
    process.residuals.VertexPrimaryLabel = cms.InputTag('offlinePrimaryVerticesWithBS')
    process.residuals.BeamSpotConfig = 'WithBS'
    
process.residuals.RunOnData = False
if options.isData:
    process.residuals.RunOnData = True

process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string("output.root"),
                                   closeFileFast = cms.untracked.bool(True)
                                   )

process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(False) )

process.p = cms.Path(
    process.offlinePrimaryVerticesRerun*
    process.residuals
)
