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
   
   bool pv_IsValid;
   bool pv_IsFake;
   int pv_NTracks;
   float pv_SumTrackPt;
   float pv_SumTrackPt2;
   float pv_chi2;
   int pv_ndof;
   float pv_x;
   float pv_y;
   float pv_z;
   float pv_xError;
   float pv_yError;
   float pv_zError;

   bool pv_IsValid_p1;
   bool pv_IsFake_p1;
   int pv_NTracks_p1;
   float pv_SumTrackPt_p1;
   float pv_SumTrackPt2_p1;
   float pv_chi2_p1;
   int pv_ndof_p1;
   float pv_x_p1;
   float pv_y_p1;
   float pv_z_p1;
   float pv_xError_p1;
   float pv_yError_p1;
   float pv_zError_p1;

   bool pv_IsValid_p2;
   bool pv_IsFake_p2;
   int pv_NTracks_p2;
   float pv_SumTrackPt_p2;
   float pv_SumTrackPt2_p2;
   float pv_chi2_p2;
   int pv_ndof_p2;
   float pv_x_p2;
   float pv_y_p2;
   float pv_z_p2;
   float pv_xError_p2;
   float pv_yError_p2;
   float pv_zError_p2;
   
   std::vector<float> trk_pt;
   std::vector<float> trk_px;
   std::vector<float> trk_py;
   std::vector<float> trk_pz;
   std::vector<float> trk_p;
   std::vector<float> trk_eta;
   std::vector<float> trk_phi;
   
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

};

#endif
