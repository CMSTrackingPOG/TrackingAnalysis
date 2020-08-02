import FWCore.ParameterSet.Config as cms

vertexResponsesAndTrueResolutions = cms.EDAnalyzer("VertexResponsesAndTrueResolutions",

  TrackLabel = cms.InputTag("generalTracks"),     
  TkMinPt = cms.double(0.3),
  TkMinXLayers = cms.int32(7),
  TkMaxMissedOuterLayers = cms.int32(4),
  TkMaxMissedInnerLayers = cms.int32(0),

  BeamSpotLabel = cms.InputTag("offlineBeamSpot"),
  BeamSpotConfig = cms.string(""),
  
  VertexLabel = cms.InputTag("offlinePrimaryVertices"),       
  VtxTracksSizeMin = cms.int32(2),
  VtxTracksSizeMax = cms.int32(300),

  VtxErrorXMin = cms.double(0.0015),
  VtxErrorXMax = cms.double(0.0037),
  VtxErrorYMin = cms.double(0.0015),
  VtxErrorYMax = cms.double(0.0037),
  VtxErrorZMin = cms.double(0.0020),
  VtxErrorZMax = cms.double(0.0036),
)
