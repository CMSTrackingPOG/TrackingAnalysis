import FWCore.ParameterSet.Config as cms
from FWCore.ParameterSet.VarParsing import VarParsing

options = VarParsing('analysis')
options.register('withBS', False, VarParsing.multiplicity.singleton, VarParsing.varType.bool, 'Primary vertex reconstruction with BS constraint')
options.register('isData', False, VarParsing.multiplicity.singleton, VarParsing.varType.bool, 'Run on data')
options.register('doTruth', False, VarParsing.multiplicity.singleton, VarParsing.varType.bool, 'Include MC truth information')
options.register('addLost', False, VarParsing.multiplicity.singleton, VarParsing.varType.bool, 'Add lost tracks')
options.register('is2016preVFP', False, VarParsing.multiplicity.singleton, VarParsing.varType.bool, 'Run on 2016 preVFP MC')
options.register('is2017', False, VarParsing.multiplicity.singleton, VarParsing.varType.bool, 'Run on 2017 MC')
options.register('is2018', False, VarParsing.multiplicity.singleton, VarParsing.varType.bool, 'Run on 2018 MC')
options.parseArguments()

readFiles = cms.untracked.vstring()
secFiles = cms.untracked.vstring()

source = cms.Source("PoolSource",fileNames = readFiles, secondaryFileNames = secFiles)

readFiles.extend( [
#'/store/data/Run2016B/JetHT/MINIAOD/ver2_HIPM_UL2016_MiniAODv2-v2/230000/45CF386B-286D-7545-B097-238839251127.root'
'/store/mc/RunIISummer20UL16MiniAODAPVv2/QCD_Pt-15to7000_TuneCP5_Flat2018_13TeV_pythia8/MINIAODSIM/106X_mcRun2_asymptotic_preVFP_v11-v1/120000/D8B934DD-D7FA-7E4C-9513-FA2300288486.root'
]);

process = cms.Process("IpResiduals")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport = cms.untracked.PSet( reportEvery = cms.untracked.int32(1) )

process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

if options.isData: process.GlobalTag.globaltag = '106X_dataRun2_v35'
elif options.is2018: process.GlobalTag.globaltag = '106X_upgrade2018_realistic_v15_L1v1'
elif options.is2017: process.GlobalTag.globaltag = '106X_mc2017_realistic_v8'
elif options.is2016preVFP: process.GlobalTag.globaltag = '106X_mcRun2_asymptotic_preVFP_v11'
else: process.GlobalTag.globaltag = '106X_mcRun2_asymptotic_v17'

process.load("CondCore.CondDB.CondDB_cfi")
process.load('Configuration.Geometry.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')

#process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(50) )

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

#primVtx = process.offlinePrimaryVertices
#process.offlinePrimaryVerticesRerun = primVtx.clone( TrackLabel = cms.InputTag("lostTracks"),
#                                                     maxDistanceToBeam = primVtx.vertexCollections[0].maxDistanceToBeam,
#                                                     useBeamConstraint = primVtx.vertexCollections[0].useBeamConstraint,
#                                                     chi2cutoff = primVtx.vertexCollections[0].chi2cutoff,
#                                                     verbose = False,
#                                                     algorithm = primVtx.vertexCollections[0].algorithm,
#                                                     minNdof = primVtx.vertexCollections[0].minNdof,
#                                                     label = primVtx.vertexCollections[0].label
#)

#if options.withBS:
#    process.offlinePrimaryVerticesRerun = primVtx.clone( TrackLabel = cms.InputTag("lostTracks"),
#                                                         maxDistanceToBeam = primVtx.vertexCollections[1].maxDistanceToBeam,
#                                                         useBeamConstraint = primVtx.vertexCollections[1].useBeamConstraint,
#                                                         chi2cutoff = primVtx.vertexCollections[1].chi2cutoff,
#                                                         verbose = False,
#                                                         algorithm = primVtx.vertexCollections[1].algorithm,
#                                                         minNdof = primVtx.vertexCollections[1].minNdof,
#                                                         label = primVtx.vertexCollections[1].label
#    )

#print process.offlinePrimaryVerticesRerun.dumpPython()

process.load('TrackingAnalysis.EDAnalyzers.residuals_cfi')
process.residuals.BeamSpotConfig = ''
if options.withBS:
    process.residuals.BeamSpotConfig = 'WithBS'

process.residuals.AddLostTracks = False
if options.addLost:
    process.residuals.AddLostTracks = True
    
process.residuals.RunOnData = False
if options.isData:
    process.residuals.RunOnData = True

process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string("output.root"),
                                   closeFileFast = cms.untracked.bool(True)
                                   )

process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(False) )

process.p = cms.Path(
#    process.offlinePrimaryVerticesRerun*
    process.residuals
)
