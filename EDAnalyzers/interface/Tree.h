#ifndef TREE_H
#define TREE_H

#include <TTree.h>
#include <string>
#include <iostream>
#include <vector>

#define null -777

class ResTree
{ 
   
 public:
   
   ResTree(TTree *_tree);
   
   TTree *tree;
   
   void Init();
   void CreateBranches(int buff, bool runOnData);

   int ev_run;
   Int_t ev_id;
   int ev_lumi;
   int ev_bunchCrossing;
   int ev_orbitNumber;
   float ev_rho;
   int ev_nPV;

   bool trig_ZeroBias_pass;
   
   bool trig_ZeroBias_part0_pass;
   bool trig_ZeroBias_part1_pass;
   bool trig_ZeroBias_part2_pass;
   bool trig_ZeroBias_part3_pass;
   bool trig_ZeroBias_part4_pass;
   bool trig_ZeroBias_part5_pass;
   bool trig_ZeroBias_part6_pass;
   bool trig_ZeroBias_part7_pass;

   bool trig_PFJet40_pass;
   bool trig_PFJet60_pass;
   bool trig_PFJet80_pass;
   bool trig_PFJet140_pass;
   bool trig_PFJet200_pass;
   bool trig_PFJet260_pass;
   bool trig_PFJet320_pass;
   bool trig_PFJet400_pass;
   bool trig_PFJet450_pass;
   bool trig_PFJet500_pass;
   bool trig_PFJet550_pass;

   bool trig_AK4PFJet30_pass;
   bool trig_AK4PFJet50_pass;
   bool trig_AK4PFJet80_pass;
   bool trig_AK4PFJet100_pass;
   bool trig_AK4PFJet120_pass;
   
   bool trig_PFHT180_pass;
   bool trig_PFHT250_pass;
   bool trig_PFHT370_pass;
   bool trig_PFHT430_pass;
   bool trig_PFHT510_pass;
   bool trig_PFHT590_pass;
   bool trig_PFHT680_pass;
   bool trig_PFHT780_pass;
   bool trig_PFHT890_pass;
   bool trig_PFHT1050_pass;
   bool trig_PFHT350_pass;
   
   int mc_pu_intime_NumInt;
   int mc_pu_trueNumInt;
   int mc_pu_before_npu;
   int mc_pu_after_npu;
   
   int mc_pu_Npvi;
   std::vector<int> mc_pu_Nzpositions;
   std::vector<int> mc_pu_BunchCrossing;
   std::vector<std::vector<float> > mc_pu_zpositions;
   std::vector<std::vector<float> > mc_pu_sumpT_lowpT;
   std::vector<std::vector<float> > mc_pu_sumpT_highpT;
   std::vector<std::vector<int> > mc_pu_ntrks_lowpT;
   std::vector<std::vector<int> > mc_pu_ntrks_highpT;
   
   int bs_type;
   float bs_x0;
   float bs_y0;
   float bs_z0;
   float bs_x_zpv;
   float bs_y_zpv;
   float bs_sigmaZ;
   float bs_dxdz;
   float bs_dydz;
   float bs_BeamWidthX;
   float bs_BeamWidthY;
   float bs_x0Error;
   float bs_y0Error;
   float bs_z0Error;
   float bs_sigmaZ0Error;
   float bs_dxdzError;
   float bs_dydzError;
   float bs_BeamWidthXError;
   float bs_BeamWidthYError;
   float bs_emittanceX;
   float bs_emittanceY;
   float bs_betaStar;
   
   std::vector<bool> pv_IsValid;
   std::vector<bool> pv_IsFake;
   std::vector<int> pv_NTracks;
   std::vector<float> pv_SumTrackPt;
   std::vector<float> pv_SumTrackPt2;
   std::vector<float> pv_fracHighPurity;
   std::vector<float> pv_chi2;
   std::vector<int> pv_ndof;
   std::vector<float> pv_x;
   std::vector<float> pv_y;
   std::vector<float> pv_z;
   std::vector<float> pv_xError;
   std::vector<float> pv_yError;
   std::vector<float> pv_zError;
   
   std::vector<std::vector<float> > pv_trk_weight;
   std::vector<std::vector<bool> > pv_trk_isHighPurity;
   std::vector<std::vector<int> > pv_trk_algo;
   std::vector<std::vector<int> > pv_trk_originalAlgo;
   
   std::vector<std::vector<int> > pv_trk_idx;
   
   std::vector<std::vector<float> > pv_trk_pt;
   std::vector<std::vector<float> > pv_trk_px;
   std::vector<std::vector<float> > pv_trk_py;
   std::vector<std::vector<float> > pv_trk_pz;
   std::vector<std::vector<float> > pv_trk_p;
   std::vector<std::vector<float> > pv_trk_eta;
   std::vector<std::vector<float> > pv_trk_phi;
   
   std::vector<std::vector<int> > pv_trk_nTrackerLayers;
   std::vector<std::vector<int> > pv_trk_nPixelBarrelLayers;
   std::vector<std::vector<int> > pv_trk_nPixelEndcapLayers;
   std::vector<std::vector<int> > pv_trk_nStripLayers;
   
   std::vector<std::vector<int> > pv_trk_nValid;
   std::vector<std::vector<float> > pv_trk_fValid;
   std::vector<std::vector<int> > pv_trk_nValidTracker;
   std::vector<std::vector<int> > pv_trk_nValidPixelBarrel;
   std::vector<std::vector<int> > pv_trk_nValidPixelEndcap;
   std::vector<std::vector<int> > pv_trk_nValidStrip;
   
   std::vector<std::vector<int> > pv_trk_nMissed;
   std::vector<std::vector<int> > pv_trk_nMissedOut;
   std::vector<std::vector<int> > pv_trk_nMissedIn;
   std::vector<std::vector<int> > pv_trk_nMissedTrackerOut;
   std::vector<std::vector<int> > pv_trk_nMissedTrackerIn;
   std::vector<std::vector<int> > pv_trk_nMissedPixelBarrelOut;
   std::vector<std::vector<int> > pv_trk_nMissedPixelBarrelIn;
   std::vector<std::vector<int> > pv_trk_nMissedPixelEndcapOut;
   std::vector<std::vector<int> > pv_trk_nMissedPixelEndcapIn;
   
   std::vector<std::vector<bool> > pv_trk_hasPixelBarrelLayer1;
   std::vector<std::vector<bool> > pv_trk_hasPixelEndcapLayer1;
   std::vector<std::vector<bool> > pv_trk_hasPixelBarrelLayer2;
   std::vector<std::vector<bool> > pv_trk_hasPixelEndcapLayer2;
   std::vector<std::vector<bool> > pv_trk_hasPixelBarrelLayer3;
   std::vector<std::vector<bool> > pv_trk_hasPixelEndcapLayer3;
   std::vector<std::vector<bool> > pv_trk_hasPixelBarrelLayer4;
   std::vector<std::vector<bool> > pv_trk_hasPixelEndcapLayer4;
   
   std::vector<std::vector<int> > pv_trk_quality;
   std::vector<std::vector<float> > pv_trk_normalizedChi2;
   std::vector<std::vector<int> > pv_trk_ndof;
   std::vector<std::vector<int> > pv_trk_charge;
   std::vector<std::vector<float> > pv_trk_qoverp;
   std::vector<std::vector<float> > pv_trk_qoverpError;
   std::vector<std::vector<float> > pv_trk_theta;
   std::vector<std::vector<float> > pv_trk_thetaError;
   std::vector<std::vector<float> > pv_trk_lambda;
   std::vector<std::vector<float> > pv_trk_lambdaError;
   std::vector<std::vector<float> > pv_trk_ptError;
   std::vector<std::vector<float> > pv_trk_etaError;
   std::vector<std::vector<float> > pv_trk_phiError;
   
   std::vector<std::vector<float> > pv_trk_d0;
   std::vector<std::vector<float> > pv_trk_dz;
   std::vector<std::vector<float> > pv_trk_d0_pv;
   std::vector<std::vector<float> > pv_trk_dz_pv;
   std::vector<std::vector<float> > pv_trk_d0_bs;
   std::vector<std::vector<float> > pv_trk_d0_bs_zpca;
   std::vector<std::vector<float> > pv_trk_d0_bs_zpv;
   std::vector<std::vector<float> > pv_trk_dz_bs;
   std::vector<std::vector<float> > pv_trk_d0Err;
   std::vector<std::vector<float> > pv_trk_dzErr;
   
   std::vector<bool> pv_mc_hasMatch;
   std::vector<std::vector<float> >  pv_mc_matchQuality;
   std::vector<std::vector<bool> > pv_mc_isFake;
   std::vector<std::vector<bool> > pv_mc_isPrimaryVertex;
   std::vector<std::vector<bool> > pv_mc_isSecondaryVertex;
   std::vector<std::vector<bool> > pv_mc_isTertiaryVertex;
   std::vector<std::vector<bool> > pv_mc_isSignalEvent;
   std::vector<std::vector<bool> > pv_mc_isBWeakDecay;
   std::vector<std::vector<bool> > pv_mc_isCWeakDecay;
   std::vector<std::vector<bool> > pv_mc_isTauDecay;
   std::vector<std::vector<bool> > pv_mc_isKsDecay;
   std::vector<std::vector<bool> > pv_mc_isLambdaDecay;
   std::vector<std::vector<bool> > pv_mc_isJpsiDecay;
   std::vector<std::vector<bool> > pv_mc_isXiDecay;
   std::vector<std::vector<bool> > pv_mc_isOmegaDecay;
   std::vector<std::vector<bool> > pv_mc_isSigmaPlusDecay;
   std::vector<std::vector<bool> > pv_mc_isSigmaMinusDecay;
   std::vector<std::vector<bool> > pv_mc_isLongLivedDecay;
   
   std::vector<std::vector<bool> > pv_mc_isKnownProcess;
   std::vector<std::vector<bool> > pv_mc_isUndefinedProcess;
   std::vector<std::vector<bool> > pv_mc_isUnknownProcess;
   std::vector<std::vector<bool> > pv_mc_isPrimaryProcess;
   std::vector<std::vector<bool> > pv_mc_isHadronicProcess;
   std::vector<std::vector<bool> > pv_mc_isDecayProcess;
   std::vector<std::vector<bool> > pv_mc_isComptonProcess;
   std::vector<std::vector<bool> > pv_mc_isAnnihilationProcess;
   std::vector<std::vector<bool> > pv_mc_isEIoniProcess;
   std::vector<std::vector<bool> > pv_mc_isHIoniProcess;
   std::vector<std::vector<bool> > pv_mc_isMuIoniProcess;
   std::vector<std::vector<bool> > pv_mc_isPhotonProcess;
   std::vector<std::vector<bool> > pv_mc_isMuPairProdProcess;
   std::vector<std::vector<bool> > pv_mc_isConversionsProcess;
   std::vector<std::vector<bool> > pv_mc_isEBremProcess;
   std::vector<std::vector<bool> > pv_mc_isSynchrotronRadiationProcess;
   std::vector<std::vector<bool> > pv_mc_isMuBremProcess;
   std::vector<std::vector<bool> > pv_mc_isMuNuclProcess;
   std::vector<std::vector<bool> > pv_mc_isUnknown;

   std::vector<std::vector<bool> > pv_mc_inVolume;
   std::vector<std::vector<float> > pv_mc_x;
   std::vector<std::vector<float> > pv_mc_y;
   std::vector<std::vector<float> > pv_mc_z;
   std::vector<std::vector<float> > pv_mc_t;
   std::vector<std::vector<int> > pv_mc_nGenVtx;
   std::vector<std::vector<int> > pv_mc_nSimVtx;
   std::vector<std::vector<int> > pv_mc_nDaughterTracks;
   std::vector<std::vector<int> > pv_mc_nSourceTracks;
   
   std::vector<bool> pv_IsValid_p1;
   std::vector<bool> pv_IsFake_p1;
   std::vector<int> pv_NTracks_p1;
   std::vector<float> pv_SumTrackPt_p1;
   std::vector<float> pv_SumTrackPt2_p1;
   std::vector<float> pv_fracHighPurity_p1;
   std::vector<std::vector<int> > pv_vtxTkIdx_p1;
   std::vector<float> pv_chi2_p1;
   std::vector<int> pv_ndof_p1;
   std::vector<float> pv_x_p1;
   std::vector<float> pv_y_p1;
   std::vector<float> pv_z_p1;
   std::vector<float> pv_xError_p1;
   std::vector<float> pv_yError_p1;
   std::vector<float> pv_zError_p1;
   
   std::vector<bool> pv_IsValid_p2;
   std::vector<bool> pv_IsFake_p2;
   std::vector<int> pv_NTracks_p2;
   std::vector<float> pv_SumTrackPt_p2;
   std::vector<float> pv_SumTrackPt2_p2;
   std::vector<float> pv_fracHighPurity_p2;
   std::vector<std::vector<int> > pv_vtxTkIdx_p2;
   std::vector<float> pv_chi2_p2;
   std::vector<int> pv_ndof_p2;
   std::vector<float> pv_x_p2;
   std::vector<float> pv_y_p2;
   std::vector<float> pv_z_p2;
   std::vector<float> pv_xError_p2;
   std::vector<float> pv_yError_p2;
   std::vector<float> pv_zError_p2;
   
   int pfjet_n;
   std::vector<float> pfjet_pt;
   std::vector<float> pfjet_eta;
   std::vector<float> pfjet_phi;
   std::vector<float> pfjet_E;

   std::vector<bool> trk_mc_hasMatch;
   std::vector<std::vector<float> > trk_mc_matchQuality;
   
   std::vector<std::vector<int> > trk_mc_pdgId;   
   std::vector<std::vector<int> > trk_mc_origin;
   std::vector<std::vector<int> > trk_mc_status;
   
   std::vector<std::vector<float> > trk_mc_pt;
   std::vector<std::vector<float> > trk_mc_px;
   std::vector<std::vector<float> > trk_mc_py;
   std::vector<std::vector<float> > trk_mc_pz;
   std::vector<std::vector<float> > trk_mc_E;
   std::vector<std::vector<float> > trk_mc_p;
   std::vector<std::vector<float> > trk_mc_eta;
   std::vector<std::vector<float> > trk_mc_phi;
   
   std::vector<std::vector<int> > trk_mc_numberOfHits;
   std::vector<std::vector<int> > trk_mc_numberOfTrackerHits;
   std::vector<std::vector<int> > trk_mc_numberOfTrackerLayers;
   
   std::vector<std::vector<float> > trk_mc_dxy_center;
   std::vector<std::vector<float> > trk_mc_dz_center;
   std::vector<std::vector<float> > trk_mc_dxy_pv;
   std::vector<std::vector<float> > trk_mc_dz_pv;
   std::vector<std::vector<float> > trk_mc_dxy_bs;
   std::vector<std::vector<float> > trk_mc_dz_bs;
   
   std::vector<std::vector<float> > trk_mc_dxy_tp_center;
   std::vector<std::vector<float> > trk_mc_dz_tp_center;
   std::vector<std::vector<float> > trk_mc_dxy_tp_pv;
   std::vector<std::vector<float> > trk_mc_dz_tp_pv;
   std::vector<std::vector<float> > trk_mc_dxy_tp_bs;
   std::vector<std::vector<float> > trk_mc_dz_tp_bs;
   
   std::vector<std::vector<float> > trk_mc_vtx_x;
   std::vector<std::vector<float> > trk_mc_vtx_y;
   std::vector<std::vector<float> > trk_mc_vtx_z;
   std::vector<std::vector<float> > trk_mc_vtx_pca_x;
   std::vector<std::vector<float> > trk_mc_vtx_pca_y;
   std::vector<std::vector<float> > trk_mc_vtx_pca_z;
   
   std::vector<std::vector<bool> > trk_mc_isFake;
   std::vector<std::vector<bool> > trk_mc_isBad;
   std::vector<std::vector<bool> > trk_mc_isBadInnerHits;
   std::vector<std::vector<bool> > trk_mc_isSharedInnerHits;
   std::vector<std::vector<bool> > trk_mc_isSignalEvent;
   std::vector<std::vector<bool> > trk_mc_isTrackerSimHits;
   std::vector<std::vector<bool> > trk_mc_isBottom;
   std::vector<std::vector<bool> > trk_mc_isCharm;
   std::vector<std::vector<bool> > trk_mc_isLight;
   std::vector<std::vector<bool> > trk_mc_isMuon;
   
   std::vector<std::vector<bool> > trk_mc_isBWeakDecay;
   std::vector<std::vector<bool> > trk_mc_isCWeakDecay;
   std::vector<std::vector<bool> > trk_mc_isChargePionDecay;
   std::vector<std::vector<bool> > trk_mc_isChargeKaonDecay;
   std::vector<std::vector<bool> > trk_mc_isTauDecay;
   std::vector<std::vector<bool> > trk_mc_isKsDecay;
   std::vector<std::vector<bool> > trk_mc_isLambdaDecay;
   std::vector<std::vector<bool> > trk_mc_isJpsiDecay;
   std::vector<std::vector<bool> > trk_mc_isXiDecay;
   std::vector<std::vector<bool> > trk_mc_isOmegaDecay;
   std::vector<std::vector<bool> > trk_mc_isSigmaPlusDecay;
   std::vector<std::vector<bool> > trk_mc_isSigmaMinusDecay;
   std::vector<std::vector<bool> > trk_mc_isLongLivedDecay;
	     
   std::vector<std::vector<bool> > trk_mc_isKnownProcess;
   std::vector<std::vector<bool> > trk_mc_isUndefinedProcess;
   std::vector<std::vector<bool> > trk_mc_isUnknownProcess;
   std::vector<std::vector<bool> > trk_mc_isPrimaryProcess;
   std::vector<std::vector<bool> > trk_mc_isHadronicProcess;
   std::vector<std::vector<bool> > trk_mc_isDecayProcess;
   std::vector<std::vector<bool> > trk_mc_isComptonProcess;
   std::vector<std::vector<bool> > trk_mc_isAnnihilationProcess;
   std::vector<std::vector<bool> > trk_mc_isEIoniProcess;
   std::vector<std::vector<bool> > trk_mc_isHIoniProcess;
   std::vector<std::vector<bool> > trk_mc_isMuIoniProcess;
   std::vector<std::vector<bool> > trk_mc_isPhotonProcess;
   std::vector<std::vector<bool> > trk_mc_isMuPairProdProcess;
   std::vector<std::vector<bool> > trk_mc_isConversionsProcess;
   std::vector<std::vector<bool> > trk_mc_isEBremProcess;
   std::vector<std::vector<bool> > trk_mc_isSynchrotronRadiationProcess;
   std::vector<std::vector<bool> > trk_mc_isMuBremProcess;
   std::vector<std::vector<bool> > trk_mc_isMuNuclProcess;
   
   std::vector<std::vector<bool> > trk_mc_isFromBWeakDecayMuon;
   std::vector<std::vector<bool> > trk_mc_isFromCWeakDecayMuon;
   std::vector<std::vector<bool> > trk_mc_isDecayOnFlightMuon;
   std::vector<std::vector<bool> > trk_mc_isFromChargePionMuon;
   std::vector<std::vector<bool> > trk_mc_isFromChargeKaonMuon;
   
   std::vector<std::vector<bool> > trk_mc_isPrimaryVertex;
   std::vector<std::vector<bool> > trk_mc_isSecondaryVertex;
   std::vector<std::vector<bool> > trk_mc_isTertiaryVertex;
   
   std::vector<std::vector<bool> > trk_mc_isUnknown;

   std::vector<float> trk_pt;
   std::vector<float> trk_px;
   std::vector<float> trk_py;
   std::vector<float> trk_pz;
   std::vector<float> trk_p;
   std::vector<float> trk_eta;
   std::vector<float> trk_phi;
   
   std::vector<int> trk_idx;
   
   std::vector<int> trk_nTrackerLayers;
   std::vector<int> trk_nPixelBarrelLayers;
   std::vector<int> trk_nPixelEndcapLayers;
   std::vector<int> trk_nStripLayers;
   
   std::vector<int> trk_nValid;
   std::vector<float> trk_fValid;
   std::vector<int> trk_nValidTracker;
   std::vector<int> trk_nValidPixelBarrel;
   std::vector<int> trk_nValidPixelEndcap;
   std::vector<int> trk_nValidStrip;

   std::vector<int> trk_nMissed;
   std::vector<int> trk_nMissedOut;
   std::vector<int> trk_nMissedIn;
   std::vector<int> trk_nMissedTrackerOut;
   std::vector<int> trk_nMissedTrackerIn;
   std::vector<int> trk_nMissedPixelBarrelOut;
   std::vector<int> trk_nMissedPixelBarrelIn;
   std::vector<int> trk_nMissedPixelEndcapOut;
   std::vector<int> trk_nMissedPixelEndcapIn;
   
   std::vector<bool> trk_hasPixelBarrelLayer1;
   std::vector<bool> trk_hasPixelEndcapLayer1;
   std::vector<bool> trk_hasPixelBarrelLayer2;
   std::vector<bool> trk_hasPixelEndcapLayer2;
   std::vector<bool> trk_hasPixelBarrelLayer3;
   std::vector<bool> trk_hasPixelEndcapLayer3;
   std::vector<bool> trk_hasPixelBarrelLayer4;
   std::vector<bool> trk_hasPixelEndcapLayer4;
   
   std::vector<int> trk_quality;
   std::vector<bool> trk_isHighPurity;
   std::vector<float> trk_normalizedChi2;
   std::vector<int> trk_ndof;
   std::vector<int> trk_charge;
   std::vector<float> trk_qoverp;
   std::vector<float> trk_qoverpError;
   std::vector<float> trk_theta;
   std::vector<float> trk_thetaError;
   std::vector<float> trk_lambda;
   std::vector<float> trk_lambdaError;
   std::vector<float> trk_ptError;
   std::vector<float> trk_etaError;
   std::vector<float> trk_phiError;
   
   std::vector<float> trk_d0;
   std::vector<float> trk_dz;
   std::vector<float> trk_d0_pv;
   std::vector<float> trk_dz_pv;
   std::vector<float> trk_d0_bs;
   std::vector<float> trk_d0_bs_zpca;
   std::vector<float> trk_d0_bs_zpv;
   std::vector<float> trk_dz_bs;
   std::vector<float> trk_d0Err;
   std::vector<float> trk_dzErr;
   std::vector<float> trk_d0_pv_NoRefit;
   std::vector<float> trk_dz_pv_NoRefit;

   // Tracks from TrackJets
   
   std::vector<bool> trk_jet_found;
   
   std::vector<float> trk_jet_pt;
   std::vector<float> trk_jet_eta;
   std::vector<float> trk_jet_phi;
   std::vector<int> trk_jet_nTracks;
   
   std::vector<float> trk_jet_pv_x;
   std::vector<float> trk_jet_pv_y;
   std::vector<float> trk_jet_pv_z;

   std::vector<bool> trk_jetTrk_found;
   
   std::vector<float> trk_jetTrk_deltaR;
   
   std::vector<float> trk_jetTrk_pt;
   std::vector<float> trk_jetTrk_px;
   std::vector<float> trk_jetTrk_py;
   std::vector<float> trk_jetTrk_pz;
   std::vector<float> trk_jetTrk_p;
   std::vector<float> trk_jetTrk_eta;
   std::vector<float> trk_jetTrk_phi;
   
   std::vector<int> trk_jetTrk_nTrackerLayers;
   std::vector<int> trk_jetTrk_nPixelBarrelLayers;
   std::vector<int> trk_jetTrk_nPixelEndcapLayers;
   std::vector<int> trk_jetTrk_nStripLayers;
   
   std::vector<int> trk_jetTrk_nValid;
   std::vector<float> trk_jetTrk_fValid;
   std::vector<int> trk_jetTrk_nValidTracker;
   std::vector<int> trk_jetTrk_nValidPixelBarrel;
   std::vector<int> trk_jetTrk_nValidPixelEndcap;
   std::vector<int> trk_jetTrk_nValidStrip;
   
   std::vector<int> trk_jetTrk_nMissed;
   std::vector<int> trk_jetTrk_nMissedOut;
   std::vector<int> trk_jetTrk_nMissedIn;
   std::vector<int> trk_jetTrk_nMissedTrackerOut;
   std::vector<int> trk_jetTrk_nMissedTrackerIn;
   std::vector<int> trk_jetTrk_nMissedPixelBarrelOut;
   std::vector<int> trk_jetTrk_nMissedPixelBarrelIn;
   std::vector<int> trk_jetTrk_nMissedPixelEndcapOut;
   std::vector<int> trk_jetTrk_nMissedPixelEndcapIn;
   
   std::vector<bool> trk_jetTrk_hasPixelBarrelLayer1;
   std::vector<bool> trk_jetTrk_hasPixelEndcapLayer1;
   std::vector<bool> trk_jetTrk_hasPixelBarrelLayer2;
   std::vector<bool> trk_jetTrk_hasPixelEndcapLayer2;
   std::vector<bool> trk_jetTrk_hasPixelBarrelLayer3;
   std::vector<bool> trk_jetTrk_hasPixelEndcapLayer3;
   std::vector<bool> trk_jetTrk_hasPixelBarrelLayer4;
   std::vector<bool> trk_jetTrk_hasPixelEndcapLayer4;
   
   std::vector<int> trk_jetTrk_quality;
   std::vector<float> trk_jetTrk_normalizedChi2;
   std::vector<int> trk_jetTrk_ndof;
   std::vector<int> trk_jetTrk_charge;
   std::vector<float> trk_jetTrk_qoverp;
   std::vector<float> trk_jetTrk_qoverpError;
   std::vector<float> trk_jetTrk_theta;
   std::vector<float> trk_jetTrk_thetaError;
   std::vector<float> trk_jetTrk_lambda;
   std::vector<float> trk_jetTrk_lambdaError;
   std::vector<float> trk_jetTrk_ptError;
   std::vector<float> trk_jetTrk_etaError;
   std::vector<float> trk_jetTrk_phiError;
   
   std::vector<float> trk_jetTrk_d0;
   std::vector<float> trk_jetTrk_dz;
   std::vector<float> trk_jetTrk_d0_pv;
   std::vector<float> trk_jetTrk_dz_pv;
   std::vector<float> trk_jetTrk_d0_bs;
   std::vector<float> trk_jetTrk_d0_bs_zpca;
   std::vector<float> trk_jetTrk_d0_bs_zpv;
   std::vector<float> trk_jetTrk_dz_bs;
   std::vector<float> trk_jetTrk_d0Err;
   std::vector<float> trk_jetTrk_dzErr;
   std::vector<float> trk_jetTrk_d0_pv_NoRefit;
   std::vector<float> trk_jetTrk_dz_pv_NoRefit;

   // Tracks from PFJets
   
   std::vector<bool> trk_pfjet_found;
   
   std::vector<float> trk_pfjet_pt;
   std::vector<float> trk_pfjet_eta;
   std::vector<float> trk_pfjet_phi;
   std::vector<int> trk_pfjet_nTracks;

   std::vector<bool> trk_pfjetTrk_found;
   
   std::vector<float> trk_pfjetTrk_deltaR;
   
   std::vector<float> trk_pfjetTrk_pt;
   std::vector<float> trk_pfjetTrk_px;
   std::vector<float> trk_pfjetTrk_py;
   std::vector<float> trk_pfjetTrk_pz;
   std::vector<float> trk_pfjetTrk_p;
   std::vector<float> trk_pfjetTrk_eta;
   std::vector<float> trk_pfjetTrk_phi;
   
   std::vector<int> trk_pfjetTrk_nTrackerLayers;
   std::vector<int> trk_pfjetTrk_nPixelBarrelLayers;
   std::vector<int> trk_pfjetTrk_nPixelEndcapLayers;
   std::vector<int> trk_pfjetTrk_nStripLayers;
   
   std::vector<int> trk_pfjetTrk_nValid;
   std::vector<float> trk_pfjetTrk_fValid;
   std::vector<int> trk_pfjetTrk_nValidTracker;
   std::vector<int> trk_pfjetTrk_nValidPixelBarrel;
   std::vector<int> trk_pfjetTrk_nValidPixelEndcap;
   std::vector<int> trk_pfjetTrk_nValidStrip;
   
   std::vector<int> trk_pfjetTrk_nMissed;
   std::vector<int> trk_pfjetTrk_nMissedOut;
   std::vector<int> trk_pfjetTrk_nMissedIn;
   std::vector<int> trk_pfjetTrk_nMissedTrackerOut;
   std::vector<int> trk_pfjetTrk_nMissedTrackerIn;
   std::vector<int> trk_pfjetTrk_nMissedPixelBarrelOut;
   std::vector<int> trk_pfjetTrk_nMissedPixelBarrelIn;
   std::vector<int> trk_pfjetTrk_nMissedPixelEndcapOut;
   std::vector<int> trk_pfjetTrk_nMissedPixelEndcapIn;
   
   std::vector<bool> trk_pfjetTrk_hasPixelBarrelLayer1;
   std::vector<bool> trk_pfjetTrk_hasPixelEndcapLayer1;
   std::vector<bool> trk_pfjetTrk_hasPixelBarrelLayer2;
   std::vector<bool> trk_pfjetTrk_hasPixelEndcapLayer2;
   std::vector<bool> trk_pfjetTrk_hasPixelBarrelLayer3;
   std::vector<bool> trk_pfjetTrk_hasPixelEndcapLayer3;
   std::vector<bool> trk_pfjetTrk_hasPixelBarrelLayer4;
   std::vector<bool> trk_pfjetTrk_hasPixelEndcapLayer4;
   
   std::vector<int> trk_pfjetTrk_quality;
   std::vector<float> trk_pfjetTrk_normalizedChi2;
   std::vector<int> trk_pfjetTrk_ndof;
   std::vector<int> trk_pfjetTrk_charge;
   std::vector<float> trk_pfjetTrk_qoverp;
   std::vector<float> trk_pfjetTrk_qoverpError;
   std::vector<float> trk_pfjetTrk_theta;
   std::vector<float> trk_pfjetTrk_thetaError;
   std::vector<float> trk_pfjetTrk_lambda;
   std::vector<float> trk_pfjetTrk_lambdaError;
   std::vector<float> trk_pfjetTrk_ptError;
   std::vector<float> trk_pfjetTrk_etaError;
   std::vector<float> trk_pfjetTrk_phiError;
   
   std::vector<float> trk_pfjetTrk_d0;
   std::vector<float> trk_pfjetTrk_dz;
   std::vector<float> trk_pfjetTrk_d0_pv;
   std::vector<float> trk_pfjetTrk_dz_pv;
   std::vector<float> trk_pfjetTrk_d0_bs;
   std::vector<float> trk_pfjetTrk_d0_bs_zpca;
   std::vector<float> trk_pfjetTrk_d0_bs_zpv;
   std::vector<float> trk_pfjetTrk_dz_bs;
   std::vector<float> trk_pfjetTrk_d0Err;
   std::vector<float> trk_pfjetTrk_dzErr;
   std::vector<float> trk_pfjetTrk_d0_pv_NoRefit;
   std::vector<float> trk_pfjetTrk_dz_pv_NoRefit;
};

#endif
