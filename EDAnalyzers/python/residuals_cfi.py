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
                           
                           trackProducer = cms.untracked.InputTag('lostTracks'),
                           trackAssociator = cms.untracked.InputTag("trackAssociatorByPull"),
                           hitAssociator = cms.PSet(associatePixel = cms.bool(False), # False when running on AOD
                                                    associateStrip = cms.bool(False), # False when running on AOD
                                                    associateRecoTracks = cms.bool(True)
                           ),
                           
                           TkFilterParameters = cms.PSet(algorithm=cms.string('filter'),
                                                         maxNormalizedChi2 = cms.double(20.0),
                                                         minPixelLayersWithHits = cms.int32(2),
                                                         minSiliconLayersWithHits = cms.int32(5),
                                                         maxD0Significance = cms.double(5.0),
                                                         minPt = cms.double(0.0),
                                                         maxEta = cms.double(5.0),
                                                         trackQuality = cms.string("any")
                           ),                                
                           
                           TkClusParameters = cms.PSet(algorithm = cms.string("DA_vect"),
                                                       TkDAClusParameters = cms.PSet(coolingFactor = cms.double(0.6),  #  moderate annealing speed
                                                                                     Tmin = cms.double(4.),            #  end of annealing
                                                                                     Tpurge = cms.double(2.0),         # cleaning
                                                                                     Tstop = cms.double(0.5),          # end of annealing
                                                                                     uniquetrkweight = cms.double(0.8), # require at least two tracks with this weight at T=Tpurge
                                                                                     zmerge = cms.double(1e-2),        # merge intermediat clusters separated by less than zmerge
                                                                                     vertexSize = cms.double(0.01),    #  ~ resolution / sqrt(Tmin)
                                                                                     d0CutOff = cms.double(3.),        # downweight high IP tracks 
                                                                                     dzCutOff = cms.double(4.)         # outlier rejection after freeze-out (T<Tmin)
                                                       )
                           ),

                           VxFitterParameters = cms.PSet(algorithm=cms.string('AVF'),
                                                         minNdof = cms.double(0.0),
                                                         maxDistanceToBeam = cms.double(1.0)
                           ),

                           VxFitterBSParameters = cms.PSet(algorithm=cms.string('AVFBS'),
                                                           minNdof = cms.double(2.0),
                                                           maxDistanceToBeam = cms.double(1.0)
                           ),
                           
                           # Pileup
                           puInfoLabel = cms.InputTag("slimmedAddPileupInfo"),

                           # Selection of Tracks
                           TrackLabel = cms.InputTag("lostTracks"),
                           TkMinPt = cms.double(0.0),
                           TkMinXLayers = cms.int32(7),
                           TkMaxMissedOuterLayers = cms.int32(4),
                           TkMaxMissedInnerLayers = cms.int32(0),

                           # Selection of Vertices
                           VertexLabel = cms.InputTag("offlineSlimmedPrimaryVertices"),
#                           VertexPrimaryLabel = cms.InputTag("offlinePrimaryVertices"),
#                           VertexPrimaryLabel = cms.InputTag("offlinePrimaryVertices","WithBS"),
                           VtxTracksSizeMin = cms.int32(2),
                           VtxTracksSizeMax = cms.int32(1000),

                           # Track jets
                           TrackJetsLabel = cms.InputTag("ak4TrackJets","","RECO"),

                           # PF jets
                           PFJetsLabel = cms.InputTag("slimmedJets"),
                           
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
