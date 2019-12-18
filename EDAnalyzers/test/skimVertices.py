import FWCore.ParameterSet.Config as cms

readFiles = cms.untracked.vstring()
secFiles = cms.untracked.vstring() 
source = cms.Source ("PoolSource",fileNames = readFiles, secondaryFileNames = secFiles)
readFiles.extend( [  #before run introducing prescaling for jetTriggers (135445)
'/store/data/Commissioning10/MinimumBias/RAW-RECO/GOODCOLL-May27thSkim_v5/0020/FE49901C-9E72-DF11-8384-0018F3D0962C.root',
'/store/data/Commissioning10/MinimumBias/RAW-RECO/GOODCOLL-May27thSkim_v5/0020/F2BFED72-9E72-DF11-A02B-0018F3D096E8.root',
'/store/data/Commissioning10/MinimumBias/RAW-RECO/GOODCOLL-May27thSkim_v5/0020/EA91F5D8-9D72-DF11-A90C-002618943894.root',
'/store/data/Commissioning10/MinimumBias/RAW-RECO/GOODCOLL-May27thSkim_v5/0020/E86BF28A-9E72-DF11-845E-0018F3D096E8.root',
'/store/data/Commissioning10/MinimumBias/RAW-RECO/GOODCOLL-May27thSkim_v5/0020/E673C3FA-9D72-DF11-97F1-0018F3D09600.root',
'/store/data/Commissioning10/MinimumBias/RAW-RECO/GOODCOLL-May27thSkim_v5/0020/E2EBC24B-9D72-DF11-96B4-0018F3D0962C.root',
'/store/data/Commissioning10/MinimumBias/RAW-RECO/GOODCOLL-May27thSkim_v5/0020/E06FF96A-9D72-DF11-A9C3-002618943894.root',
'/store/data/Commissioning10/MinimumBias/RAW-RECO/GOODCOLL-May27thSkim_v5/0020/E05C5938-9C72-DF11-83A6-0018F3D09600.root',
'/store/data/Commissioning10/MinimumBias/RAW-RECO/GOODCOLL-May27thSkim_v5/0020/DAAFA989-9E72-DF11-B55E-0018F3D0961A.root',
'/store/data/Commissioning10/MinimumBias/RAW-RECO/GOODCOLL-May27thSkim_v5/0020/D8311F50-9C72-DF11-9CF1-0018F3D09600.root',
'/store/data/Commissioning10/MinimumBias/RAW-RECO/GOODCOLL-May27thSkim_v5/0020/D4D9AA9D-9C72-DF11-B547-001A92811732.root',
'/store/data/Commissioning10/MinimumBias/RAW-RECO/GOODCOLL-May27thSkim_v5/0020/D2079EE2-9C72-DF11-849B-0018F3D096E8.root',
'/store/data/Commissioning10/MinimumBias/RAW-RECO/GOODCOLL-May27thSkim_v5/0020/D055CE02-9C72-DF11-BC06-0018F3D096F8.root',
'/store/data/Commissioning10/MinimumBias/RAW-RECO/GOODCOLL-May27thSkim_v5/0020/C8E21D60-9C72-DF11-8544-0018F3D0962C.root',
'/store/data/Commissioning10/MinimumBias/RAW-RECO/GOODCOLL-May27thSkim_v5/0020/C2D15BE8-9C72-DF11-A391-002618943894.root',
'/store/data/Commissioning10/MinimumBias/RAW-RECO/GOODCOLL-May27thSkim_v5/0020/BED87A53-9C72-DF11-B2E9-003048678FA6.root',
'/store/data/Commissioning10/MinimumBias/RAW-RECO/GOODCOLL-May27thSkim_v5/0020/B4CE00C4-9B72-DF11-A7B1-0018F3D096E8.root',
'/store/data/Commissioning10/MinimumBias/RAW-RECO/GOODCOLL-May27thSkim_v5/0020/B481581A-9F72-DF11-82A0-0018F3D0962C.root',
'/store/data/Commissioning10/MinimumBias/RAW-RECO/GOODCOLL-May27thSkim_v5/0020/AAC54980-9F72-DF11-B3EE-0018F3D096E8.root',
'/store/data/Commissioning10/MinimumBias/RAW-RECO/GOODCOLL-May27thSkim_v5/0020/9EFBB04F-9E72-DF11-AE61-002618943894.root',
'/store/data/Commissioning10/MinimumBias/RAW-RECO/GOODCOLL-May27thSkim_v5/0020/9EA59E66-9F72-DF11-871E-002618943858.root',
'/store/data/Commissioning10/MinimumBias/RAW-RECO/GOODCOLL-May27thSkim_v5/0020/9A67196C-9D72-DF11-91B1-001A92811732.root',
'/store/data/Commissioning10/MinimumBias/RAW-RECO/GOODCOLL-May27thSkim_v5/0020/96727029-A072-DF11-B8F4-001A92811732.root',
'/store/data/Commissioning10/MinimumBias/RAW-RECO/GOODCOLL-May27thSkim_v5/0020/9247EF6C-9E72-DF11-B86A-001A92811732.root',
'/store/data/Commissioning10/MinimumBias/RAW-RECO/GOODCOLL-May27thSkim_v5/0020/90D67B9B-9C72-DF11-9F48-001A92971B68.root',
'/store/data/Commissioning10/MinimumBias/RAW-RECO/GOODCOLL-May27thSkim_v5/0020/8AC1B396-9D72-DF11-9CCB-0018F3D096E8.root',
'/store/data/Commissioning10/MinimumBias/RAW-RECO/GOODCOLL-May27thSkim_v5/0020/886B2677-9D72-DF11-AB57-0018F3D09600.root',
'/store/data/Commissioning10/MinimumBias/RAW-RECO/GOODCOLL-May27thSkim_v5/0020/884284F6-9D72-DF11-9445-0018F3D0962C.root',
'/store/data/Commissioning10/MinimumBias/RAW-RECO/GOODCOLL-May27thSkim_v5/0020/7AAC36EE-9C72-DF11-B629-0018F3D0962C.root',
'/store/data/Commissioning10/MinimumBias/RAW-RECO/GOODCOLL-May27thSkim_v5/0020/7807CC04-9D72-DF11-B0FF-0018F3D09600.root',
'/store/data/Commissioning10/MinimumBias/RAW-RECO/GOODCOLL-May27thSkim_v5/0020/70F4EE84-9C72-DF11-9584-002618943858.root',
'/store/data/Commissioning10/MinimumBias/RAW-RECO/GOODCOLL-May27thSkim_v5/0020/70B6997A-9E72-DF11-8CFA-001A92971B68.root',
'/store/data/Commissioning10/MinimumBias/RAW-RECO/GOODCOLL-May27thSkim_v5/0020/6E7BEC7B-9B72-DF11-945B-001A92971B88.root',
'/store/data/Commissioning10/MinimumBias/RAW-RECO/GOODCOLL-May27thSkim_v5/0020/6C820C9B-9D72-DF11-8521-0018F3D0962C.root',
'/store/data/Commissioning10/MinimumBias/RAW-RECO/GOODCOLL-May27thSkim_v5/0020/682896E9-9B72-DF11-BA70-002618943894.root',
'/store/data/Commissioning10/MinimumBias/RAW-RECO/GOODCOLL-May27thSkim_v5/0020/60A7403D-9F72-DF11-85C2-0018F3D0961A.root',


]);



secFiles.extend( [ ]);


process = cms.Process("IpResiduals")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1000
# import of standard configurations
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.StandardSequences.GeometryExtended_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.RawToDigi_cff')
process.load('Configuration.StandardSequences.Reconstruction_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.EventContent.EventContent_cff')

process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
#process.GlobalTag.globaltag = 'START3X_V26A::All'
process.GlobalTag.globaltag = 'GR10_P_V4::All'

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100000) )

process.source = source
#======================================
# L1 (from GOODCOLL)
#======================================
process.load('L1TriggerConfig.L1GtConfigProducers.L1GtTriggerMaskTechTrigConfig_cff')
from HLTrigger.HLTfilters.hltLevel1GTSeed_cfi import hltLevel1GTSeed
process.L1MinBias = hltLevel1GTSeed.clone(
    L1TechTriggerSeeding = cms.bool(True),
    L1SeedsLogicalExpression = cms.string('(40 OR 41) AND NOT (36 OR 37 OR 38 OR 39) AND NOT ((42 AND NOT 43) OR (43 AND NOT 42))')
    )

#======================================
# HLT
#======================================
process.HLTMinBias = cms.EDFilter("HLTHighLevel",
                               TriggerResultsTag = cms.InputTag("TriggerResults","","HLT"),
                               HLTPaths = cms.vstring('HLT_L1_BscMinBiasOR_BptxPlusORMinus'), 
#                               HLTPaths = cms.vstring('HLT_Jet15U'), 
                               eventSetupPathsKey = cms.string(''), # not empty => use read paths from AlCaRecoTriggerBitsRcd via this key
                               andOr = cms.bool(True),             # how to deal with multiple triggers: True (OR) accept if ANY is true, False (AND) accept if ALL are true
                               throw = cms.bool(True)    # throw exception on unknown path names
                               )

#====================================
# vertex selector
#====================================
process.goodVertices = cms.EDFilter("VertexSelector",
   src = cms.InputTag("offlinePrimaryVertices"),
   cut = cms.string(" xError > 0.0012 &&" +
                    " xError < 0.0035 &&" +
                    " yError > 0.0012 &&" +
                    " yError < 0.0035 &&" +
                    " zError > 0.0017 &&" +
                    " zError < 0.0040 &&" +
                    " !isFake" ),
   filter = cms.bool(True),   # otherwise it won't filter the events, just produce an empty vertex collection.
)


process.out = cms.OutputModule("PoolOutputModule",
#    fileName = cms.untracked.string('vtx.DATA.Jet15U.root'),
    fileName = cms.untracked.string('vtx.DATA.MinBias.root'),
    outputCommands = cms.untracked.vstring(
    'drop *',
    'keep recoVertexs_offlinePrimaryVertices__RECO'
        ),
    SelectEvents = cms.untracked.PSet(
         SelectEvents = cms.vstring('selection')
     )
)

process.options = cms.untracked.PSet(
     wantSummary = cms.untracked.bool(True)
)

process.selection = cms.Path(#process.L1MinBias
                             process.HLTMinBias
                             #*process.goodVertices
                             )
process.e = cms.EndPath(process.out)

#process.schedule = cms.Schedule(process.p,process.e)

