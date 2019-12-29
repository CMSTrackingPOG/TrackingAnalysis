import FWCore.ParameterSet.Config as cms

residuals = cms.EDAnalyzer("Residuals",

                           # Beam spot
                           BeamSpotLabel = cms.InputTag("offlineBeamSpot"),

                           # Rho
                           RhoLabel = cms.InputTag("fixedGridRhoFastjetAll"),

                           # Trigger results
                           TriggerResultsLabel = cms.InputTag("TriggerResults","","HLT"),
                           
                           # Selection of Tracks
                           TrackLabel = cms.InputTag("generalTracks"),
                           TkMinPt = cms.double(0.3),
                           TkMinXLayers = cms.int32(7),
                           TkMaxMissedOuterLayers = cms.int32(4),
                           TkMaxMissedInnerLayers = cms.int32(0),

                           # Selection of Vertices
                           VertexLabel = cms.InputTag("offlinePrimaryVerticesRerun"),
                           VtxTracksSizeMin = cms.int32(2),
                           VtxTracksSizeMax = cms.int32(300),

                           # Vertex selection for Jet6U trigger
                           VtxErrorXMin = cms.double(0.0015),
                           VtxErrorXMax = cms.double(0.0037),
                           VtxErrorYMin = cms.double(0.0015),
                           VtxErrorYMax = cms.double(0.0037),
                           VtxErrorZMin = cms.double(0.0020),
                           VtxErrorZMax = cms.double(0.0036),

                           # Vertex selection for MinBias trigger
                           ##VtxErrorXMin = cms.double(0.0020),
                           ##VtxErrorXMax = cms.double(0.0055),
                           ##VtxErrorYMin = cms.double(0.0020),
                           ##VtxErrorYMax = cms.double(0.0055),
                           ##VtxErrorZMin = cms.double(0.0025),
                           ##VtxErrorZMax = cms.double(0.0060),
)
