import FWCore.ParameterSet.Config as cms

residuals = cms.EDAnalyzer("Residuals",

                           # Run on data
                           RunOnData = cms.bool(False),

                           # Read MC Truth information
                           DoTruth = cms.bool(False),
                           
                           # Beam spot
                           BeamSpotLabel = cms.InputTag("offlineBeamSpot"),
                           BeamSpotConfig = cms.string(""),

                           # Rho
                           RhoLabel = cms.InputTag("fixedGridRhoFastjetAll"),

                           # Trigger results
                           TriggerResultsLabel = cms.InputTag("TriggerResults", "", "HLT"),

                           # Track Associator and Tracking Particles
                           TrackingParticleLabel = cms.InputTag("mix", "MergedTrackTruth"),
                           TrackingVertexLabel = cms.InputTag("mix", "MergedTrackTruth"),
#                           TrackAssociatorLabel = cms.InputTag("quickTrackAssociatorByHits"),
#                           TrackAssociatorLabel = cms.InputTag("trackAssociatorByChi2"),
                           VertexAssociatorLabel = cms.InputTag("vertexAssociator"),
                           TrackAssociatorLabel = cms.InputTag("trackAssociatorByPull"),
                           
                           # For Track&VertexHistory
                           vertexProducer = cms.untracked.InputTag("offlinePrimaryVerticesRerun"),
                           trackingTruth = cms.untracked.InputTag("mix", "MergedTrackTruth"),
                           vertexAssociator = cms.untracked.InputTag('vertexAssociator'),
                           bestMatchByMaxValue = cms.untracked.bool(True),
                           enableRecoToSim = cms.untracked.bool(True),
                           enableSimToReco = cms.untracked.bool(False),
                           
                           hepMC = cms.untracked.InputTag("generatorSmeared"),
                           beamSpot = cms.untracked.InputTag("offlineBeamSpot"),
                           badPull = cms.untracked.double(3.0),
                           longLivedDecayLength = cms.untracked.double(1e-14), # > 10^-14 sec of lifetime
                           vertexClusteringDistance = cms.untracked.double(0.0001), # 1mum
                           numberOfInnerLayers = cms.untracked.uint32(2),
                           minTrackerSimHits = cms.untracked.uint32(3),
                           
                           trackProducer = cms.untracked.InputTag('generalTracks'),
                           trackAssociator = cms.untracked.InputTag("trackAssociatorByPull"),
                           hitAssociator = cms.PSet(associatePixel = cms.bool(False), # False when running on AOD
                                                    associateStrip = cms.bool(False), # False when running on AOD
                                                    associateRecoTracks = cms.bool(True)
                           ),
                           
                           # Pileup
                           puInfoLabel = cms.InputTag("addPileupInfo"),

                           # Selection of Tracks
                           TrackLabel = cms.InputTag("generalTracks"),
                           TkMinPt = cms.double(0.0),
                           TkMinXLayers = cms.int32(7),
                           TkMaxMissedOuterLayers = cms.int32(4),
                           TkMaxMissedInnerLayers = cms.int32(0),

                           # Selection of Vertices
                           VertexLabel = cms.InputTag("offlinePrimaryVerticesRerun"),
#                           VertexPrimaryLabel = cms.InputTag("offlinePrimaryVertices"),
#                           VertexPrimaryLabel = cms.InputTag("offlinePrimaryVertices","WithBS"),
                           VtxTracksSizeMin = cms.int32(2),
                           VtxTracksSizeMax = cms.int32(1000),

                           # Track jets
                           TrackJetsLabel = cms.InputTag("ak4TrackJets","","RECO"),

                           # PF jets
                           PFJetsLabel = cms.InputTag("ak4PFJets","","RECO"),
                           
                           # Vertex selection for Jet6U trigger
#                           VtxErrorXMin = cms.double(0.0015),
#                           VtxErrorXMax = cms.double(0.0037),
#                           VtxErrorYMin = cms.double(0.0015),
#                           VtxErrorYMax = cms.double(0.0037),
#                           VtxErrorZMin = cms.double(0.0020),
#                           VtxErrorZMax = cms.double(0.0036),

                           # Event filter
                           EventScale = cms.int32(1),
#                           TrackScale = cms.int32(100),
                           TrackScale = cms.int32(1),
                           
                           # Vertex selection for MinBias trigger
                           ##VtxErrorXMin = cms.double(0.0020),
                           ##VtxErrorXMax = cms.double(0.0055),
                           ##VtxErrorYMin = cms.double(0.0020),
                           ##VtxErrorYMax = cms.double(0.0055),
                           ##VtxErrorZMin = cms.double(0.0025),
                           ##VtxErrorZMax = cms.double(0.0060),
)
