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
   void CreateBranches(int buff);

   int ev_run;
   Int_t ev_id;
   int ev_lumi;
   float ev_rho;

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
   float pv_chi2_p2;
   int pv_ndof_p2;
   float pv_x_p2;
   float pv_y_p2;
   float pv_z_p2;
   float pv_xError_p2;
   float pv_yError_p2;
   float pv_zError_p2;

   std::vector<float> trk_pt;
   std::vector<float> trk_p;
   std::vector<float> trk_eta;
   std::vector<float> trk_phi;
   std::vector<int> trk_nXLayers;
   std::vector<int> trk_nMissedOut;
   std::vector<int> trk_nMissedIn;
   std::vector<int> trk_hasPXL;
   std::vector<int> trk_quality;
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
   
};

#endif
