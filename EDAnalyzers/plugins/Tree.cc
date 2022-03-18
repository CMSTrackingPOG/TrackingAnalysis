#include "TrackingAnalysis/EDAnalyzers/interface/Tree.h"

ResTree::ResTree(TTree* _tree)
{
   tree = _tree;
}

void ResTree::Init()
{
   ev_run = null;
   ev_id = null;
   ev_lumi = null;
   ev_bunchCrossing = null;
   ev_orbitNumber = null;
   ev_time = null;
   ev_rho = null;
   ev_nPV = null;
   
   trig_ZeroBias_pass = 0;
   
   trig_ZeroBias_part0_pass = 0;
   trig_ZeroBias_part1_pass = 0;
   trig_ZeroBias_part2_pass = 0;
   trig_ZeroBias_part3_pass = 0;
   trig_ZeroBias_part4_pass = 0;
   trig_ZeroBias_part5_pass = 0;
   trig_ZeroBias_part6_pass = 0;
   trig_ZeroBias_part7_pass = 0;

   trig_PFJet40_pass = 0;
   trig_PFJet60_pass = 0;
   trig_PFJet80_pass = 0;
   trig_PFJet140_pass = 0;
   trig_PFJet200_pass = 0;
   trig_PFJet260_pass = 0;
   trig_PFJet320_pass = 0;
   trig_PFJet400_pass = 0;
   trig_PFJet450_pass = 0;
   trig_PFJet500_pass = 0;
   trig_PFJet550_pass = 0;

   trig_AK4PFJet30_pass = 0;
   trig_AK4PFJet50_pass = 0;
   trig_AK4PFJet80_pass = 0;
   trig_AK4PFJet100_pass = 0;
   trig_AK4PFJet120_pass = 0;

   trig_PFHT180_pass = 0;
   trig_PFHT250_pass = 0;
   trig_PFHT370_pass = 0;
   trig_PFHT430_pass = 0;
   trig_PFHT510_pass = 0;
   trig_PFHT590_pass = 0;
   trig_PFHT680_pass = 0;
   trig_PFHT780_pass = 0;
   trig_PFHT890_pass = 0;
   trig_PFHT1050_pass = 0;
   trig_PFHT350_pass = 0;
   
   mc_pu_intime_NumInt = null;
   mc_pu_trueNumInt = null;
   mc_pu_before_npu = null;
   mc_pu_after_npu = null;
   
   mc_pu_Npvi = null;
   mc_pu_Nzpositions.clear();
   mc_pu_BunchCrossing.clear();
   for(unsigned int i=0;i<mc_pu_zpositions.size();i++) mc_pu_zpositions[i].clear(); 
   mc_pu_zpositions.clear();
   for(unsigned int i=0;i<mc_pu_sumpT_lowpT.size();i++) mc_pu_sumpT_lowpT[i].clear(); 
   mc_pu_sumpT_lowpT.clear();
   for(unsigned int i=0;i<mc_pu_sumpT_highpT.size();i++) mc_pu_sumpT_highpT[i].clear(); 
   mc_pu_sumpT_highpT.clear();
   for(unsigned int i=0;i<mc_pu_ntrks_lowpT.size();i++) mc_pu_ntrks_lowpT[i].clear(); 
   mc_pu_ntrks_lowpT.clear();
   for(unsigned int i=0;i<mc_pu_ntrks_highpT.size();i++) mc_pu_ntrks_highpT[i].clear(); 
   mc_pu_ntrks_highpT.clear();
     
   bs_type = null;
   bs_x0 = null;
   bs_y0 = null;
   bs_z0 = null;
   bs_x_zpv = null;
   bs_y_zpv = null;
   bs_sigmaZ = null;
   bs_dxdz = null;
   bs_dydz = null;
   bs_BeamWidthX = null;
   bs_BeamWidthY = null;
   bs_x0Error = null;
   bs_y0Error = null;
   bs_z0Error = null;
   bs_sigmaZ0Error = null;
   bs_dxdzError = null;
   bs_dydzError = null;
   bs_BeamWidthXError = null;
   bs_BeamWidthYError = null;
   bs_emittanceX = null;
   bs_emittanceY = null;
   bs_betaStar = null;

   int pv_ntrk = pv_NTracks.size();

   pv_IsValid.clear();
   pv_IsFake.clear();
   pv_NTracks.clear();
   pv_SumTrackPt.clear();
   pv_SumTrackPt2.clear();
   pv_fracHighPurity.clear();
   pv_chi2.clear();
   pv_ndof.clear();
   pv_x.clear();
   pv_y.clear();
   pv_z.clear();
   pv_xError.clear();
   pv_yError.clear();
   pv_zError.clear();
   
   for(int i=0;i<pv_ntrk;i++)
     {	
	pv_trk_weight[i].clear();
	pv_trk_isHighPurity[i].clear();
	pv_trk_algo[i].clear();
	pv_trk_originalAlgo[i].clear();
   
	pv_trk_idx[i].clear();
	
	pv_trk_pvN[i].clear();
	pv_trk_pv1N[i].clear();
	pv_trk_pv2N[i].clear();
	
	pv_trk_pvunbiased_IsValid[i].clear();
	pv_trk_pvunbiased_IsFake[i].clear();
	pv_trk_pvunbiased_NTracks[i].clear();
	pv_trk_pvunbiased_SumTrackPt[i].clear();
	pv_trk_pvunbiased_SumTrackPt2[i].clear();
	pv_trk_pvunbiased_fracHighPurity[i].clear();
	pv_trk_pvunbiased_chi2[i].clear();
	pv_trk_pvunbiased_ndof[i].clear();
	pv_trk_pvunbiased_x[i].clear();
	pv_trk_pvunbiased_y[i].clear();
	pv_trk_pvunbiased_z[i].clear();
	pv_trk_pvunbiased_xError[i].clear();
	pv_trk_pvunbiased_yError[i].clear();
	pv_trk_pvunbiased_zError[i].clear();
	
	pv_trk_d0_pvunbiased[i].clear();
	pv_trk_dz_pvunbiased[i].clear();
	pv_trk_d0_bs_zpvunbiased[i].clear();

	pv_trk_pvunbiased_IsValid_p1[i].clear();
	pv_trk_pvunbiased_IsFake_p1[i].clear();
	pv_trk_pvunbiased_NTracks_p1[i].clear();
	pv_trk_pvunbiased_SumTrackPt_p1[i].clear();
	pv_trk_pvunbiased_SumTrackPt2_p1[i].clear();
	pv_trk_pvunbiased_fracHighPurity_p1[i].clear();
	pv_trk_pvunbiased_chi2_p1[i].clear();
	pv_trk_pvunbiased_ndof_p1[i].clear();
	pv_trk_pvunbiased_x_p1[i].clear();
	pv_trk_pvunbiased_y_p1[i].clear();
	pv_trk_pvunbiased_z_p1[i].clear();
	pv_trk_pvunbiased_xError_p1[i].clear();
	pv_trk_pvunbiased_yError_p1[i].clear();
	pv_trk_pvunbiased_zError_p1[i].clear();

	pv_trk_d0_pvunbiased_p1[i].clear();
	pv_trk_dz_pvunbiased_p1[i].clear();
	pv_trk_d0_bs_zpvunbiased_p1[i].clear();

	pv_trk_pvunbiased_IsValid_p2[i].clear();
	pv_trk_pvunbiased_IsFake_p2[i].clear();
	pv_trk_pvunbiased_NTracks_p2[i].clear();
	pv_trk_pvunbiased_SumTrackPt_p2[i].clear();
	pv_trk_pvunbiased_SumTrackPt2_p2[i].clear();
	pv_trk_pvunbiased_fracHighPurity_p2[i].clear();
	pv_trk_pvunbiased_chi2_p2[i].clear();
	pv_trk_pvunbiased_ndof_p2[i].clear();
	pv_trk_pvunbiased_x_p2[i].clear();
	pv_trk_pvunbiased_y_p2[i].clear();
	pv_trk_pvunbiased_z_p2[i].clear();
	pv_trk_pvunbiased_xError_p2[i].clear();
	pv_trk_pvunbiased_yError_p2[i].clear();
	pv_trk_pvunbiased_zError_p2[i].clear();
	
	pv_trk_d0_pvunbiased_p2[i].clear();
	pv_trk_dz_pvunbiased_p2[i].clear();
	pv_trk_d0_bs_zpvunbiased_p2[i].clear();

	pv_trk_pt[i].clear();
	pv_trk_px[i].clear();
	pv_trk_py[i].clear();
	pv_trk_pz[i].clear();
	pv_trk_p[i].clear();
	pv_trk_eta[i].clear();
	pv_trk_phi[i].clear();

	pv_trk_nTrackerLayers[i].clear();
	pv_trk_nPixelBarrelLayers[i].clear();
	pv_trk_nPixelEndcapLayers[i].clear();
	pv_trk_nStripLayers[i].clear();
   
	pv_trk_nValid[i].clear();
	pv_trk_fValid[i].clear();
	pv_trk_nValidTracker[i].clear();
	pv_trk_nValidPixelBarrel[i].clear();
	pv_trk_nValidPixelEndcap[i].clear();
	pv_trk_nValidStrip[i].clear();
   
	pv_trk_nMissed[i].clear();
	pv_trk_nMissedOut[i].clear();
	pv_trk_nMissedIn[i].clear();
	pv_trk_nMissedTrackerOut[i].clear();
	pv_trk_nMissedTrackerIn[i].clear();
	pv_trk_nMissedPixelBarrelOut[i].clear();
	pv_trk_nMissedPixelBarrelIn[i].clear();
	pv_trk_nMissedPixelEndcapOut[i].clear();
	pv_trk_nMissedPixelEndcapIn[i].clear();
   
	pv_trk_hasPixelBarrelLayer1[i].clear();
	pv_trk_hasPixelEndcapLayer1[i].clear();
	pv_trk_hasPixelBarrelLayer2[i].clear();
	pv_trk_hasPixelEndcapLayer2[i].clear();
	pv_trk_hasPixelBarrelLayer3[i].clear();
	pv_trk_hasPixelEndcapLayer3[i].clear();
	pv_trk_hasPixelBarrelLayer4[i].clear();
	pv_trk_hasPixelEndcapLayer4[i].clear();

	pv_trk_quality[i].clear();
	pv_trk_normalizedChi2[i].clear();
	pv_trk_ndof[i].clear();
	pv_trk_charge[i].clear();
	pv_trk_qoverp[i].clear();
	pv_trk_qoverpError[i].clear();
	pv_trk_theta[i].clear();
	pv_trk_thetaError[i].clear();
	pv_trk_lambda[i].clear();
	pv_trk_lambdaError[i].clear();
	pv_trk_ptError[i].clear();
	pv_trk_etaError[i].clear();
	pv_trk_phiError[i].clear();
	
	pv_trk_d0[i].clear();
	pv_trk_dz[i].clear();
	pv_trk_d0_pv[i].clear();
	pv_trk_dz_pv[i].clear();
	pv_trk_d0_bs[i].clear();
	pv_trk_d0_bs_zpca[i].clear();
	pv_trk_d0_bs_zpv[i].clear();
	pv_trk_dz_bs[i].clear();
	pv_trk_d0Err[i].clear();
	pv_trk_dzErr[i].clear();
     }   

   pv_trk_weight.clear();
   pv_trk_isHighPurity.clear();
   pv_trk_algo.clear();
   pv_trk_originalAlgo.clear();
   
   pv_trk_idx.clear();
   
   pv_trk_pvN.clear();
   pv_trk_pv1N.clear();
   pv_trk_pv2N.clear();
   
   pv_trk_pvunbiased_IsValid.clear();
   pv_trk_pvunbiased_IsFake.clear();
   pv_trk_pvunbiased_NTracks.clear();
   pv_trk_pvunbiased_SumTrackPt.clear();
   pv_trk_pvunbiased_SumTrackPt2.clear();
   pv_trk_pvunbiased_fracHighPurity.clear();
   pv_trk_pvunbiased_chi2.clear();
   pv_trk_pvunbiased_ndof.clear();
   pv_trk_pvunbiased_x.clear();
   pv_trk_pvunbiased_y.clear();
   pv_trk_pvunbiased_z.clear();
   pv_trk_pvunbiased_xError.clear();
   pv_trk_pvunbiased_yError.clear();
   pv_trk_pvunbiased_zError.clear();
   
   pv_trk_d0_pvunbiased.clear();
   pv_trk_dz_pvunbiased.clear();
   pv_trk_d0_bs_zpvunbiased.clear();
   
   pv_trk_pvunbiased_IsValid_p1.clear();
   pv_trk_pvunbiased_IsFake_p1.clear();
   pv_trk_pvunbiased_NTracks_p1.clear();
   pv_trk_pvunbiased_SumTrackPt_p1.clear();
   pv_trk_pvunbiased_SumTrackPt2_p1.clear();
   pv_trk_pvunbiased_fracHighPurity_p1.clear();
   pv_trk_pvunbiased_chi2_p1.clear();
   pv_trk_pvunbiased_ndof_p1.clear();
   pv_trk_pvunbiased_x_p1.clear();
   pv_trk_pvunbiased_y_p1.clear();
   pv_trk_pvunbiased_z_p1.clear();
   pv_trk_pvunbiased_xError_p1.clear();
   pv_trk_pvunbiased_yError_p1.clear();
   pv_trk_pvunbiased_zError_p1.clear();
   
   pv_trk_d0_pvunbiased_p1.clear();
   pv_trk_dz_pvunbiased_p1.clear();
   pv_trk_d0_bs_zpvunbiased_p1.clear();

   pv_trk_pvunbiased_IsValid_p2.clear();
   pv_trk_pvunbiased_IsFake_p2.clear();
   pv_trk_pvunbiased_NTracks_p2.clear();
   pv_trk_pvunbiased_SumTrackPt_p2.clear();
   pv_trk_pvunbiased_SumTrackPt2_p2.clear();
   pv_trk_pvunbiased_fracHighPurity_p2.clear();
   pv_trk_pvunbiased_chi2_p2.clear();
   pv_trk_pvunbiased_ndof_p2.clear();
   pv_trk_pvunbiased_x_p2.clear();
   pv_trk_pvunbiased_y_p2.clear();
   pv_trk_pvunbiased_z_p2.clear();
   pv_trk_pvunbiased_xError_p2.clear();
   pv_trk_pvunbiased_yError_p2.clear();
   pv_trk_pvunbiased_zError_p2.clear();
   
   pv_trk_d0_pvunbiased_p2.clear();
   pv_trk_dz_pvunbiased_p2.clear();
   pv_trk_d0_bs_zpvunbiased_p2.clear();
   
   pv_trk_pt.clear();
   pv_trk_px.clear();
   pv_trk_py.clear();
   pv_trk_pz.clear();
   pv_trk_p.clear();
   pv_trk_eta.clear();
   pv_trk_phi.clear();
   
   pv_trk_nTrackerLayers.clear();
   pv_trk_nPixelBarrelLayers.clear();
   pv_trk_nPixelEndcapLayers.clear();
   pv_trk_nStripLayers.clear();
   
   pv_trk_nValid.clear();
   pv_trk_fValid.clear();
   pv_trk_nValidTracker.clear();
   pv_trk_nValidPixelBarrel.clear();
   pv_trk_nValidPixelEndcap.clear();
   pv_trk_nValidStrip.clear();
   
   pv_trk_nMissed.clear();
   pv_trk_nMissedOut.clear();
   pv_trk_nMissedIn.clear();
   pv_trk_nMissedTrackerOut.clear();
   pv_trk_nMissedTrackerIn.clear();
   pv_trk_nMissedPixelBarrelOut.clear();
   pv_trk_nMissedPixelBarrelIn.clear();
   pv_trk_nMissedPixelEndcapOut.clear();
   pv_trk_nMissedPixelEndcapIn.clear();
   
   pv_trk_hasPixelBarrelLayer1.clear();
   pv_trk_hasPixelEndcapLayer1.clear();
   pv_trk_hasPixelBarrelLayer2.clear();
   pv_trk_hasPixelEndcapLayer2.clear();
   pv_trk_hasPixelBarrelLayer3.clear();
   pv_trk_hasPixelEndcapLayer3.clear();
   pv_trk_hasPixelBarrelLayer4.clear();
   pv_trk_hasPixelEndcapLayer4.clear();
   
   pv_trk_quality.clear();
   pv_trk_normalizedChi2.clear();
   pv_trk_ndof.clear();
   pv_trk_charge.clear();
   pv_trk_qoverp.clear();
   pv_trk_qoverpError.clear();
   pv_trk_theta.clear();
   pv_trk_thetaError.clear();
   pv_trk_lambda.clear();
   pv_trk_lambdaError.clear();
   pv_trk_ptError.clear();
   pv_trk_etaError.clear();
   pv_trk_phiError.clear();
   
   pv_trk_d0.clear();
   pv_trk_dz.clear();
   pv_trk_d0_pv.clear();
   pv_trk_dz_pv.clear();
   pv_trk_d0_bs.clear();
   pv_trk_d0_bs_zpca.clear();
   pv_trk_d0_bs_zpv.clear();
   pv_trk_dz_bs.clear();
   pv_trk_d0Err.clear();
   pv_trk_dzErr.clear();
   
   pv_IsValid_p1.clear();
   pv_IsFake_p1.clear();
   pv_NTracks_p1.clear();
   pv_SumTrackPt_p1.clear();
   pv_SumTrackPt2_p1.clear();
   pv_fracHighPurity_p1.clear();
   for(unsigned int i=0;i<pv_vtxTkIdx_p1.size();i++) pv_vtxTkIdx_p1[i].clear();
   pv_vtxTkIdx_p1.clear();
   pv_chi2_p1.clear();
   pv_ndof_p1.clear();
   pv_x_p1.clear();
   pv_y_p1.clear();
   pv_z_p1.clear();
   pv_xError_p1.clear();
   pv_yError_p1.clear();
   pv_zError_p1.clear();

   pv_IsValid_p2.clear();
   pv_IsFake_p2.clear();
   pv_NTracks_p2.clear();
   pv_SumTrackPt_p2.clear();
   pv_SumTrackPt2_p2.clear();
   pv_fracHighPurity_p2.clear();
   for(unsigned int i=0;i<pv_vtxTkIdx_p2.size();i++) pv_vtxTkIdx_p2[i].clear();
   pv_vtxTkIdx_p2.clear();
   pv_chi2_p2.clear();
   pv_ndof_p2.clear();
   pv_x_p2.clear();
   pv_y_p2.clear();
   pv_z_p2.clear();
   pv_xError_p2.clear();
   pv_yError_p2.clear();
   pv_zError_p2.clear();
   
   pfjet_n = null;
   pfjet_pt.clear();
   pfjet_eta.clear();
   pfjet_phi.clear();
   pfjet_E.clear();
   
   // Tracks from PFJets
   
/*   trk_pfjet_found.clear();
   
   trk_pfjet_pt.clear();
   trk_pfjet_eta.clear();
   trk_pfjet_phi.clear();
   trk_pfjet_nTracks.clear();

   trk_pfjetTrk_found.clear();
   
   trk_pfjetTrk_deltaR.clear();
   
   trk_pfjetTrk_pt.clear();
   trk_pfjetTrk_px.clear();
   trk_pfjetTrk_py.clear();
   trk_pfjetTrk_pz.clear();
   trk_pfjetTrk_p.clear();
   trk_pfjetTrk_eta.clear();
   trk_pfjetTrk_phi.clear();
   
   trk_pfjetTrk_nTrackerLayers.clear();
   trk_pfjetTrk_nPixelBarrelLayers.clear();
   trk_pfjetTrk_nPixelEndcapLayers.clear();
   trk_pfjetTrk_nStripLayers.clear();
   
   trk_pfjetTrk_nValid.clear();
   trk_pfjetTrk_fValid.clear();
   trk_pfjetTrk_nValidTracker.clear();
   trk_pfjetTrk_nValidPixelBarrel.clear();
   trk_pfjetTrk_nValidPixelEndcap.clear();
   trk_pfjetTrk_nValidStrip.clear();
   
   trk_pfjetTrk_nMissed.clear();
   trk_pfjetTrk_nMissedOut.clear();
   trk_pfjetTrk_nMissedIn.clear();
   trk_pfjetTrk_nMissedTrackerOut.clear();
   trk_pfjetTrk_nMissedTrackerIn.clear();
   trk_pfjetTrk_nMissedPixelBarrelOut.clear();
   trk_pfjetTrk_nMissedPixelBarrelIn.clear();
   trk_pfjetTrk_nMissedPixelEndcapOut.clear();
   trk_pfjetTrk_nMissedPixelEndcapIn.clear();
   
   trk_pfjetTrk_hasPixelBarrelLayer1.clear();
   trk_pfjetTrk_hasPixelEndcapLayer1.clear();
   trk_pfjetTrk_hasPixelBarrelLayer2.clear();
   trk_pfjetTrk_hasPixelEndcapLayer2.clear();
   trk_pfjetTrk_hasPixelBarrelLayer3.clear();
   trk_pfjetTrk_hasPixelEndcapLayer3.clear();
   trk_pfjetTrk_hasPixelBarrelLayer4.clear();
   trk_pfjetTrk_hasPixelEndcapLayer4.clear();
   
   trk_pfjetTrk_quality.clear();
   trk_pfjetTrk_normalizedChi2.clear();
   trk_pfjetTrk_ndof.clear();
   trk_pfjetTrk_charge.clear();
   trk_pfjetTrk_qoverp.clear();
   trk_pfjetTrk_qoverpError.clear();
   trk_pfjetTrk_theta.clear();
   trk_pfjetTrk_thetaError.clear();
   trk_pfjetTrk_lambda.clear();
   trk_pfjetTrk_lambdaError.clear();
   trk_pfjetTrk_ptError.clear();
   trk_pfjetTrk_etaError.clear();
   trk_pfjetTrk_phiError.clear();
   
   trk_pfjetTrk_d0.clear();
   trk_pfjetTrk_dz.clear();
   trk_pfjetTrk_d0_pv.clear();
   trk_pfjetTrk_dz_pv.clear();
   trk_pfjetTrk_d0_bs.clear();
   trk_pfjetTrk_d0_bs_zpca.clear();
   trk_pfjetTrk_d0_bs_zpv.clear();
   trk_pfjetTrk_dz_bs.clear();
   trk_pfjetTrk_d0Err.clear();
   trk_pfjetTrk_dzErr.clear();
   trk_pfjetTrk_d0_pv_NoRefit.clear();
   trk_pfjetTrk_dz_pv_NoRefit.clear();*/
}

void ResTree::CreateBranches(int buff = 32000, bool runOnData = false)
{
   tree->Branch("ev_run", &ev_run, "ev_run/I", buff);
   tree->Branch("ev_id", &ev_id, "ev_id/I", buff);
   tree->Branch("ev_lumi", &ev_lumi, "ev_lumi/I", buff);
   tree->Branch("ev_bunchCrossing", &ev_bunchCrossing, "ev_bunchCrossing/I", buff);
   tree->Branch("ev_orbitNumber", &ev_orbitNumber, "ev_orbitNumber/I", buff);
   tree->Branch("ev_time", &ev_time, "ev_time/I", buff);
   tree->Branch("ev_rho", &ev_rho, "ev_rho/F", buff);
   tree->Branch("ev_nPV", &ev_nPV, "ev_nPV/I", buff);
   
   tree->Branch("trig_ZeroBias_pass", &trig_ZeroBias_pass, "trig_ZeroBias_pass/O", buff);
   
   tree->Branch("trig_ZeroBias_part0_pass", &trig_ZeroBias_part0_pass, "trig_ZeroBias_part0_pass/O", buff);
   tree->Branch("trig_ZeroBias_part1_pass", &trig_ZeroBias_part1_pass, "trig_ZeroBias_part1_pass/O", buff);
   tree->Branch("trig_ZeroBias_part2_pass", &trig_ZeroBias_part2_pass, "trig_ZeroBias_part2_pass/O", buff);
   tree->Branch("trig_ZeroBias_part3_pass", &trig_ZeroBias_part3_pass, "trig_ZeroBias_part3_pass/O", buff);
   tree->Branch("trig_ZeroBias_part4_pass", &trig_ZeroBias_part4_pass, "trig_ZeroBias_part4_pass/O", buff);
   tree->Branch("trig_ZeroBias_part5_pass", &trig_ZeroBias_part5_pass, "trig_ZeroBias_part5_pass/O", buff);
   tree->Branch("trig_ZeroBias_part6_pass", &trig_ZeroBias_part6_pass, "trig_ZeroBias_part6_pass/O", buff);
   tree->Branch("trig_ZeroBias_part7_pass", &trig_ZeroBias_part7_pass, "trig_ZeroBias_part7_pass/O", buff);
   
   tree->Branch("trig_PFJet40_pass", &trig_PFJet40_pass, "trig_PFJet40_pass/O", buff);
   tree->Branch("trig_PFJet60_pass", &trig_PFJet60_pass, "trig_PFJet60_pass/O", buff);
   tree->Branch("trig_PFJet80_pass", &trig_PFJet80_pass, "trig_PFJet80_pass/O", buff);
   tree->Branch("trig_PFJet140_pass", &trig_PFJet140_pass, "trig_PFJet140_pass/O", buff);
   tree->Branch("trig_PFJet200_pass", &trig_PFJet200_pass, "trig_PFJet200_pass/O", buff);
   tree->Branch("trig_PFJet260_pass", &trig_PFJet260_pass, "trig_PFJet260_pass/O", buff);
   tree->Branch("trig_PFJet320_pass", &trig_PFJet320_pass, "trig_PFJet320_pass/O", buff);
   tree->Branch("trig_PFJet400_pass", &trig_PFJet400_pass, "trig_PFJet400_pass/O", buff);
   tree->Branch("trig_PFJet450_pass", &trig_PFJet450_pass, "trig_PFJet450_pass/O", buff);
   tree->Branch("trig_PFJet500_pass", &trig_PFJet500_pass, "trig_PFJet500_pass/O", buff);
   tree->Branch("trig_PFJet550_pass", &trig_PFJet550_pass, "trig_PFJet550_pass/O", buff);
   
   tree->Branch("trig_AK4PFJet30_pass", &trig_AK4PFJet30_pass, "trig_AK4PFJet30_pass/O", buff);
   tree->Branch("trig_AK4PFJet50_pass", &trig_AK4PFJet50_pass, "trig_AK4PFJet50_pass/O", buff);
   tree->Branch("trig_AK4PFJet80_pass", &trig_AK4PFJet80_pass, "trig_AK4PFJet80_pass/O", buff);
   tree->Branch("trig_AK4PFJet100_pass", &trig_AK4PFJet100_pass, "trig_AK4PFJet100_pass/O", buff);
   tree->Branch("trig_AK4PFJet120_pass", &trig_AK4PFJet120_pass, "trig_AK4PFJet120_pass/O", buff);
   
   tree->Branch("trig_PFHT180_pass", &trig_PFHT180_pass, "trig_PFHT180_pass/O", buff);
   tree->Branch("trig_PFHT250_pass", &trig_PFHT250_pass, "trig_PFHT250_pass/O", buff);
   tree->Branch("trig_PFHT370_pass", &trig_PFHT370_pass, "trig_PFHT370_pass/O", buff);
   tree->Branch("trig_PFHT430_pass", &trig_PFHT430_pass, "trig_PFHT430_pass/O", buff);
   tree->Branch("trig_PFHT510_pass", &trig_PFHT510_pass, "trig_PFHT510_pass/O", buff);
   tree->Branch("trig_PFHT590_pass", &trig_PFHT590_pass, "trig_PFHT590_pass/O", buff);
   tree->Branch("trig_PFHT680_pass", &trig_PFHT680_pass, "trig_PFHT680_pass/O", buff);
   tree->Branch("trig_PFHT780_pass", &trig_PFHT780_pass, "trig_PFHT780_pass/O", buff);
   tree->Branch("trig_PFHT890_pass", &trig_PFHT890_pass, "trig_PFHT890_pass/O", buff);
   tree->Branch("trig_PFHT1050_pass", &trig_PFHT1050_pass, "trig_PFHT1050_pass/O", buff);
   tree->Branch("trig_PFHT350_pass", &trig_PFHT350_pass, "trig_PFHT350_pass/O", buff);
   
   if( !runOnData ) 
     {
	tree->Branch("mc_pu_intime_NumInt", &mc_pu_intime_NumInt, "mc_pu_intime_NumInt/I", buff);
	tree->Branch("mc_pu_trueNumInt", &mc_pu_trueNumInt, "mc_pu_trueNumInt/I", buff);
	tree->Branch("mc_pu_before_npu", &mc_pu_before_npu, "mc_pu_before_npu/I", buff);
	tree->Branch("mc_pu_after_npu", &mc_pu_after_npu, "mc_pu_after_npu/I", buff);
	
	tree->Branch("mc_pu_Npvi", &mc_pu_Npvi, "mc_pu_Npvi/I", buff);
	tree->Branch("mc_pu_Nzpositions", "std::vector<int>", &mc_pu_Nzpositions, buff);
	tree->Branch("mc_pu_BunchCrossing", "std::vector<int>", &mc_pu_BunchCrossing, buff);
	tree->Branch("mc_pu_zpositions", "std::vector<std::vector<float> >", &mc_pu_zpositions, buff);
	tree->Branch("mc_pu_sumpT_lowpT", "std::vector<std::vector<float> >", &mc_pu_sumpT_lowpT, buff);
	tree->Branch("mc_pu_sumpT_highpT", "std::vector<std::vector<float> >", &mc_pu_sumpT_highpT, buff);
	tree->Branch("mc_pu_ntrks_lowpT", "std::vector<std::vector<int> >", &mc_pu_ntrks_lowpT, buff);
	tree->Branch("mc_pu_ntrks_highpT", "std::vector<std::vector<int> >", &mc_pu_ntrks_highpT, buff);
     }   
   
   tree->Branch("bs_type", &bs_type, "bs_type/I", buff);
   tree->Branch("bs_x0", &bs_x0, "bs_x0/F", buff);
   tree->Branch("bs_y0", &bs_y0, "bs_y0/F", buff);
   tree->Branch("bs_z0", &bs_z0, "bs_z0/F", buff);
   tree->Branch("bs_x_zpv", &bs_x_zpv, "bs_x_zpv/F", buff);
   tree->Branch("bs_y_zpv", &bs_y_zpv, "bs_y_zpv/F", buff);
   tree->Branch("bs_sigmaZ", &bs_sigmaZ, "bs_sigmaZ/F", buff);
   tree->Branch("bs_dxdz", &bs_dxdz, "bs_dxdz/F", buff);
   tree->Branch("bs_dydz", &bs_dydz, "bs_dydz/F", buff);
   tree->Branch("bs_BeamWidthX", &bs_BeamWidthX, "bs_BeamWidthX/F", buff);
   tree->Branch("bs_BeamWidthY", &bs_BeamWidthY, "bs_BeamWidthY/F", buff);
   tree->Branch("bs_x0Error", &bs_x0Error, "bs_x0Error/F", buff);
   tree->Branch("bs_y0Error", &bs_y0Error, "bs_y0Error/F", buff);
   tree->Branch("bs_z0Error", &bs_z0Error, "bs_z0Error/F", buff);
   tree->Branch("bs_sigmaZ0Error", &bs_sigmaZ0Error, "bs_sigmaZ0Error/F", buff);
   tree->Branch("bs_dxdzError", &bs_dxdzError, "bs_dxdzError/F", buff);
   tree->Branch("bs_dydzError", &bs_dydzError, "bs_dydzError/F", buff);
   tree->Branch("bs_BeamWidthXError", &bs_BeamWidthXError, "bs_BeamWidthXError/F", buff);
   tree->Branch("bs_BeamWidthYError", &bs_BeamWidthYError, "bs_BeamWidthYError/F", buff);
   tree->Branch("bs_emittanceX", &bs_emittanceX, "bs_emittanceX/F", buff);
   tree->Branch("bs_emittanceY", &bs_emittanceY, "bs_emittanceY/F", buff);
   tree->Branch("bs_betaStar", &bs_betaStar, "bs_betaStar/F", buff);

   tree->Branch("pv_IsValid", "std::vector<bool>", &pv_IsValid, buff);
   tree->Branch("pv_IsFake", "std::vector<bool>", &pv_IsFake, buff);
   tree->Branch("pv_NTracks", "std::vector<int>", &pv_NTracks, buff);
   tree->Branch("pv_SumTrackPt", "std::vector<float>", &pv_SumTrackPt, buff);
   tree->Branch("pv_SumTrackPt2", "std::vector<float>", &pv_SumTrackPt2, buff);
   tree->Branch("pv_fracHighPurity", "std::vector<float>", &pv_fracHighPurity, buff);
   tree->Branch("pv_chi2", "std::vector<float>", &pv_chi2, buff);
   tree->Branch("pv_ndof", "std::vector<int>", &pv_ndof, buff);
   tree->Branch("pv_x", "std::vector<float>", &pv_x, buff);
   tree->Branch("pv_y", "std::vector<float>", &pv_y, buff);
   tree->Branch("pv_z", "std::vector<float>", &pv_z, buff);
   tree->Branch("pv_xError", "std::vector<float>", &pv_xError, buff);
   tree->Branch("pv_yError", "std::vector<float>", &pv_yError, buff);
   tree->Branch("pv_zError", "std::vector<float>", &pv_zError, buff);
   
   tree->Branch("pv_trk_weight", "std::vector<std::vector<float> >", &pv_trk_weight, buff);
   tree->Branch("pv_trk_isHighPurity", "std::vector<std::vector<bool> >", &pv_trk_isHighPurity, buff);
   tree->Branch("pv_trk_algo", "std::vector<std::vector<int> >", &pv_trk_algo, buff);
   tree->Branch("pv_trk_originalAlgo", "std::vector<std::vector<int> >", &pv_trk_originalAlgo, buff);
   
   tree->Branch("pv_trk_idx", "std::vector<std::vector<int> >", &pv_trk_idx, buff);
   
   tree->Branch("pv_trk_pvN", "std::vector<std::vector<int> >", &pv_trk_pvN, buff);
   tree->Branch("pv_trk_pv1N", "std::vector<std::vector<int> >", &pv_trk_pv1N, buff);
   tree->Branch("pv_trk_pv2N", "std::vector<std::vector<int> >", &pv_trk_pv2N, buff);

   tree->Branch("pv_trk_pvunbiased_IsValid", "std::vector<std::vector<bool> >", &pv_trk_pvunbiased_IsValid, buff);
   tree->Branch("pv_trk_pvunbiased_IsFake", "std::vector<std::vector<bool> >", &pv_trk_pvunbiased_IsFake, buff);
   tree->Branch("pv_trk_pvunbiased_NTracks", "std::vector<std::vector<int> >", &pv_trk_pvunbiased_NTracks, buff);
   tree->Branch("pv_trk_pvunbiased_SumTrackPt", "std::vector<std::vector<float> >", &pv_trk_pvunbiased_SumTrackPt, buff);
   tree->Branch("pv_trk_pvunbiased_SumTrackPt2", "std::vector<std::vector<float> >", &pv_trk_pvunbiased_SumTrackPt2, buff);
   tree->Branch("pv_trk_pvunbiased_fracHighPurity", "std::vector<std::vector<float> >", &pv_trk_pvunbiased_fracHighPurity, buff);
   tree->Branch("pv_trk_pvunbiased_chi2", "std::vector<std::vector<float> >", &pv_trk_pvunbiased_chi2, buff);
   tree->Branch("pv_trk_pvunbiased_ndof", "std::vector<std::vector<int> >", &pv_trk_pvunbiased_ndof, buff);
   tree->Branch("pv_trk_pvunbiased_x", "std::vector<std::vector<float> >", &pv_trk_pvunbiased_x, buff);
   tree->Branch("pv_trk_pvunbiased_y", "std::vector<std::vector<float> >", &pv_trk_pvunbiased_y, buff);
   tree->Branch("pv_trk_pvunbiased_z", "std::vector<std::vector<float> >", &pv_trk_pvunbiased_z, buff);
   tree->Branch("pv_trk_pvunbiased_xError", "std::vector<std::vector<float> >", &pv_trk_pvunbiased_xError, buff);
   tree->Branch("pv_trk_pvunbiased_yError", "std::vector<std::vector<float> >", &pv_trk_pvunbiased_yError, buff);
   tree->Branch("pv_trk_pvunbiased_zError", "std::vector<std::vector<float> >", &pv_trk_pvunbiased_zError, buff);

   tree->Branch("pv_trk_d0_pvunbiased", "std::vector<std::vector<float> >", &pv_trk_d0_pvunbiased, buff);
   tree->Branch("pv_trk_dz_pvunbiased", "std::vector<std::vector<float> >", &pv_trk_dz_pvunbiased, buff);
   tree->Branch("pv_trk_d0_bs_zpvunbiased", "std::vector<std::vector<float> >", &pv_trk_d0_bs_zpvunbiased, buff);
   
   tree->Branch("pv_trk_pvunbiased_IsValid_p1", "std::vector<std::vector<bool> >", &pv_trk_pvunbiased_IsValid_p1, buff);
   tree->Branch("pv_trk_pvunbiased_IsFake_p1", "std::vector<std::vector<bool> >", &pv_trk_pvunbiased_IsFake_p1, buff);
   tree->Branch("pv_trk_pvunbiased_NTracks_p1", "std::vector<std::vector<int> >", &pv_trk_pvunbiased_NTracks_p1, buff);
   tree->Branch("pv_trk_pvunbiased_SumTrackPt_p1", "std::vector<std::vector<float> >", &pv_trk_pvunbiased_SumTrackPt_p1, buff);
   tree->Branch("pv_trk_pvunbiased_SumTrackPt2_p1", "std::vector<std::vector<float> >", &pv_trk_pvunbiased_SumTrackPt2_p1, buff);
   tree->Branch("pv_trk_pvunbiased_fracHighPurity_p1", "std::vector<std::vector<float> >", &pv_trk_pvunbiased_fracHighPurity_p1, buff);
   tree->Branch("pv_trk_pvunbiased_chi2_p1", "std::vector<std::vector<float> >", &pv_trk_pvunbiased_chi2_p1, buff);
   tree->Branch("pv_trk_pvunbiased_ndof_p1", "std::vector<std::vector<int> >", &pv_trk_pvunbiased_ndof_p1, buff);
   tree->Branch("pv_trk_pvunbiased_x_p1", "std::vector<std::vector<float> >", &pv_trk_pvunbiased_x_p1, buff);
   tree->Branch("pv_trk_pvunbiased_y_p1", "std::vector<std::vector<float> >", &pv_trk_pvunbiased_y_p1, buff);
   tree->Branch("pv_trk_pvunbiased_z_p1", "std::vector<std::vector<float> >", &pv_trk_pvunbiased_z_p1, buff);
   tree->Branch("pv_trk_pvunbiased_xError_p1", "std::vector<std::vector<float> >", &pv_trk_pvunbiased_xError_p1, buff);
   tree->Branch("pv_trk_pvunbiased_yError_p1", "std::vector<std::vector<float> >", &pv_trk_pvunbiased_yError_p1, buff);
   tree->Branch("pv_trk_pvunbiased_zError_p1", "std::vector<std::vector<float> >", &pv_trk_pvunbiased_zError_p1, buff);

   tree->Branch("pv_trk_d0_pvunbiased_p1", "std::vector<std::vector<float> >", &pv_trk_d0_pvunbiased_p1, buff);
   tree->Branch("pv_trk_dz_pvunbiased_p1", "std::vector<std::vector<float> >", &pv_trk_dz_pvunbiased_p1, buff);
   tree->Branch("pv_trk_d0_bs_zpvunbiased_p1", "std::vector<std::vector<float> >", &pv_trk_d0_bs_zpvunbiased_p1, buff);

   tree->Branch("pv_trk_pvunbiased_IsValid_p2", "std::vector<std::vector<bool> >", &pv_trk_pvunbiased_IsValid_p2, buff);
   tree->Branch("pv_trk_pvunbiased_IsFake_p2", "std::vector<std::vector<bool> >", &pv_trk_pvunbiased_IsFake_p2, buff);
   tree->Branch("pv_trk_pvunbiased_NTracks_p2", "std::vector<std::vector<int> >", &pv_trk_pvunbiased_NTracks_p2, buff);
   tree->Branch("pv_trk_pvunbiased_SumTrackPt_p2", "std::vector<std::vector<float> >", &pv_trk_pvunbiased_SumTrackPt_p2, buff);
   tree->Branch("pv_trk_pvunbiased_SumTrackPt2_p2", "std::vector<std::vector<float> >", &pv_trk_pvunbiased_SumTrackPt2_p2, buff);
   tree->Branch("pv_trk_pvunbiased_fracHighPurity_p2", "std::vector<std::vector<float> >", &pv_trk_pvunbiased_fracHighPurity_p2, buff);
   tree->Branch("pv_trk_pvunbiased_chi2_p2", "std::vector<std::vector<float> >", &pv_trk_pvunbiased_chi2_p2, buff);
   tree->Branch("pv_trk_pvunbiased_ndof_p2", "std::vector<std::vector<int> >", &pv_trk_pvunbiased_ndof_p2, buff);
   tree->Branch("pv_trk_pvunbiased_x_p2", "std::vector<std::vector<float> >", &pv_trk_pvunbiased_x_p2, buff);
   tree->Branch("pv_trk_pvunbiased_y_p2", "std::vector<std::vector<float> >", &pv_trk_pvunbiased_y_p2, buff);
   tree->Branch("pv_trk_pvunbiased_z_p2", "std::vector<std::vector<float> >", &pv_trk_pvunbiased_z_p2, buff);
   tree->Branch("pv_trk_pvunbiased_xError_p2", "std::vector<std::vector<float> >", &pv_trk_pvunbiased_xError_p2, buff);
   tree->Branch("pv_trk_pvunbiased_yError_p2", "std::vector<std::vector<float> >", &pv_trk_pvunbiased_yError_p2, buff);
   tree->Branch("pv_trk_pvunbiased_zError_p2", "std::vector<std::vector<float> >", &pv_trk_pvunbiased_zError_p2, buff);

   tree->Branch("pv_trk_d0_pvunbiased_p2", "std::vector<std::vector<float> >", &pv_trk_d0_pvunbiased_p2, buff);
   tree->Branch("pv_trk_dz_pvunbiased_p2", "std::vector<std::vector<float> >", &pv_trk_dz_pvunbiased_p2, buff);
   tree->Branch("pv_trk_d0_bs_zpvunbiased_p2", "std::vector<std::vector<float> >", &pv_trk_d0_bs_zpvunbiased_p2, buff);
   
   tree->Branch("pv_trk_pt", "std::vector<std::vector<float> >", &pv_trk_pt, buff);
   tree->Branch("pv_trk_px", "std::vector<std::vector<float> >", &pv_trk_px, buff);
   tree->Branch("pv_trk_py", "std::vector<std::vector<float> >", &pv_trk_py, buff);
   tree->Branch("pv_trk_pz", "std::vector<std::vector<float> >", &pv_trk_pz, buff);
   tree->Branch("pv_trk_p", "std::vector<std::vector<float> >", &pv_trk_p, buff);
   tree->Branch("pv_trk_eta", "std::vector<std::vector<float> >", &pv_trk_eta, buff);
   tree->Branch("pv_trk_phi", "std::vector<std::vector<float> >", &pv_trk_phi, buff);
   
   tree->Branch("pv_trk_nTrackerLayers", "std::vector<std::vector<int> >", &pv_trk_nTrackerLayers, buff);
   tree->Branch("pv_trk_nPixelBarrelLayers", "std::vector<std::vector<int> >", &pv_trk_nPixelBarrelLayers, buff);
   tree->Branch("pv_trk_nPixelEndcapLayers", "std::vector<std::vector<int> >", &pv_trk_nPixelEndcapLayers, buff);
   tree->Branch("pv_trk_nStripLayers", "std::vector<std::vector<int> >", &pv_trk_nStripLayers, buff);
   
   tree->Branch("pv_trk_nValid", "std::vector<std::vector<int> >", &pv_trk_nValid, buff);
   tree->Branch("pv_trk_fValid", "std::vector<std::vector<float> >", &pv_trk_fValid, buff);
   tree->Branch("pv_trk_nValidTracker", "std::vector<std::vector<int> >", &pv_trk_nValidTracker, buff);
   tree->Branch("pv_trk_nValidPixelBarrel", "std::vector<std::vector<int> >", &pv_trk_nValidPixelBarrel, buff);
   tree->Branch("pv_trk_nValidPixelEndcap", "std::vector<std::vector<int> >", &pv_trk_nValidPixelEndcap, buff);
   tree->Branch("pv_trk_nValidStrip", "std::vector<std::vector<int> >", &pv_trk_nValidStrip, buff);
   
   tree->Branch("pv_trk_nMissed", "std::vector<std::vector<int> >", &pv_trk_nMissed, buff);
   tree->Branch("pv_trk_nMissedOut", "std::vector<std::vector<int> >", &pv_trk_nMissedOut, buff);
   tree->Branch("pv_trk_nMissedIn", "std::vector<std::vector<int> >", &pv_trk_nMissedIn, buff);
   tree->Branch("pv_trk_nMissedTrackerOut", "std::vector<std::vector<int> >", &pv_trk_nMissedTrackerOut, buff);
   tree->Branch("pv_trk_nMissedTrackerIn", "std::vector<std::vector<int> >", &pv_trk_nMissedTrackerIn, buff);
   tree->Branch("pv_trk_nMissedPixelBarrelOut", "std::vector<std::vector<int> >", &pv_trk_nMissedPixelBarrelOut, buff);
   tree->Branch("pv_trk_nMissedPixelBarrelIn", "std::vector<std::vector<int> >", &pv_trk_nMissedPixelBarrelIn, buff);
   tree->Branch("pv_trk_nMissedPixelEndcapOut", "std::vector<std::vector<int> >", &pv_trk_nMissedPixelEndcapOut, buff);
   tree->Branch("pv_trk_nMissedPixelEndcapIn", "std::vector<std::vector<int> >", &pv_trk_nMissedPixelEndcapIn, buff);
   
   tree->Branch("pv_trk_hasPixelBarrelLayer1", "std::vector<std::vector<bool> >", &pv_trk_hasPixelBarrelLayer1, buff);
   tree->Branch("pv_trk_hasPixelEndcapLayer1", "std::vector<std::vector<bool> >", &pv_trk_hasPixelEndcapLayer1, buff);
   tree->Branch("pv_trk_hasPixelBarrelLayer2", "std::vector<std::vector<bool> >", &pv_trk_hasPixelBarrelLayer2, buff);
   tree->Branch("pv_trk_hasPixelEndcapLayer2", "std::vector<std::vector<bool> >", &pv_trk_hasPixelEndcapLayer2, buff);
   tree->Branch("pv_trk_hasPixelBarrelLayer3", "std::vector<std::vector<bool> >", &pv_trk_hasPixelBarrelLayer3, buff);
   tree->Branch("pv_trk_hasPixelEndcapLayer3", "std::vector<std::vector<bool> >", &pv_trk_hasPixelEndcapLayer3, buff);
   tree->Branch("pv_trk_hasPixelBarrelLayer4", "std::vector<std::vector<bool> >", &pv_trk_hasPixelBarrelLayer4, buff);
   tree->Branch("pv_trk_hasPixelEndcapLayer4", "std::vector<std::vector<bool> >", &pv_trk_hasPixelEndcapLayer4, buff);
   
   tree->Branch("pv_trk_quality", "std::vector<std::vector<int> >", &pv_trk_quality, buff);
   tree->Branch("pv_trk_normalizedChi2", "std::vector<std::vector<float> >", &pv_trk_normalizedChi2, buff);
   tree->Branch("pv_trk_ndof", "std::vector<std::vector<int> >", &pv_trk_ndof, buff);
   tree->Branch("pv_trk_charge", "std::vector<std::vector<int> >", &pv_trk_charge, buff);
   tree->Branch("pv_trk_qoverp", "std::vector<std::vector<float> >", &pv_trk_qoverp, buff);
   tree->Branch("pv_trk_qoverpError", "std::vector<std::vector<float> >", &pv_trk_qoverpError, buff);
   tree->Branch("pv_trk_theta", "std::vector<std::vector<float> >", &pv_trk_theta, buff);
   tree->Branch("pv_trk_thetaError", "std::vector<std::vector<float> >", &pv_trk_thetaError, buff);
   tree->Branch("pv_trk_lambda", "std::vector<std::vector<float> >", &pv_trk_lambda, buff);
   tree->Branch("pv_trk_lambdaError", "std::vector<std::vector<float> >", &pv_trk_lambdaError, buff);
   tree->Branch("pv_trk_ptError", "std::vector<std::vector<float> >", &pv_trk_ptError, buff);
   tree->Branch("pv_trk_etaError", "std::vector<std::vector<float> >", &pv_trk_etaError, buff);
   tree->Branch("pv_trk_phiError", "std::vector<std::vector<float> >", &pv_trk_phiError, buff);
   
   tree->Branch("pv_trk_d0", "std::vector<std::vector<float> >", &pv_trk_d0, buff);
   tree->Branch("pv_trk_dz", "std::vector<std::vector<float> >", &pv_trk_dz, buff);
   tree->Branch("pv_trk_d0_pv", "std::vector<std::vector<float> >", &pv_trk_d0_pv, buff);
   tree->Branch("pv_trk_dz_pv", "std::vector<std::vector<float> >", &pv_trk_dz_pv, buff);
   tree->Branch("pv_trk_d0_bs", "std::vector<std::vector<float> >", &pv_trk_d0_bs, buff);
   tree->Branch("pv_trk_d0_bs_zpca", "std::vector<std::vector<float> >", &pv_trk_d0_bs_zpca, buff);
   tree->Branch("pv_trk_d0_bs_zpv", "std::vector<std::vector<float> >", &pv_trk_d0_bs_zpv, buff);
   tree->Branch("pv_trk_dz_bs", "std::vector<std::vector<float> >", &pv_trk_dz_bs, buff);
   tree->Branch("pv_trk_d0Err", "std::vector<std::vector<float> >", &pv_trk_d0Err, buff);
   tree->Branch("pv_trk_dzErr", "std::vector<std::vector<float> >", &pv_trk_dzErr, buff);
   
   tree->Branch("pv_IsValid_p1", "std::vector<bool>", &pv_IsValid_p1, buff);
   tree->Branch("pv_IsFake_p1", "std::vector<bool>", &pv_IsFake_p1, buff);
   tree->Branch("pv_NTracks_p1", "std::vector<int>", &pv_NTracks_p1, buff);
   tree->Branch("pv_SumTrackPt_p1", "std::vector<float>", &pv_SumTrackPt_p1, buff);
   tree->Branch("pv_SumTrackPt2_p1", "std::vector<float>", &pv_SumTrackPt2_p1, buff);
   tree->Branch("pv_fracHighPurity_p1", "std::vector<float>", &pv_fracHighPurity_p1, buff);
   tree->Branch("pv_vtxTkIdx_p1", "std::vector<std::vector<int> >", &pv_vtxTkIdx_p1, buff);
   tree->Branch("pv_chi2_p1", "std::vector<float>", &pv_chi2_p1, buff);
   tree->Branch("pv_ndof_p1", "std::vector<int>", &pv_ndof_p1, buff);
   tree->Branch("pv_x_p1", "std::vector<float>", &pv_x_p1, buff);
   tree->Branch("pv_y_p1", "std::vector<float>", &pv_y_p1, buff);
   tree->Branch("pv_z_p1", "std::vector<float>", &pv_z_p1, buff);
   tree->Branch("pv_xError_p1", "std::vector<float>", &pv_xError_p1, buff);
   tree->Branch("pv_yError_p1", "std::vector<float>", &pv_yError_p1, buff);
   tree->Branch("pv_zError_p1", "std::vector<float>", &pv_zError_p1, buff);

   tree->Branch("pv_IsValid_p2", "std::vector<bool>", &pv_IsValid_p2, buff);
   tree->Branch("pv_IsFake_p2", "std::vector<bool>", &pv_IsFake_p2, buff);
   tree->Branch("pv_NTracks_p2", "std::vector<int>", &pv_NTracks_p2, buff);
   tree->Branch("pv_SumTrackPt_p2", "std::vector<float>", &pv_SumTrackPt_p2, buff);
   tree->Branch("pv_SumTrackPt2_p2", "std::vector<float>", &pv_SumTrackPt2_p2, buff);
   tree->Branch("pv_fracHighPurity_p2", "std::vector<float>", &pv_fracHighPurity_p2, buff);
   tree->Branch("pv_vtxTkIdx_p2", "std::vector<std::vector<int> >", &pv_vtxTkIdx_p2, buff);
   tree->Branch("pv_chi2_p2", "std::vector<float>", &pv_chi2_p2, buff);
   tree->Branch("pv_ndof_p2", "std::vector<int>", &pv_ndof_p2, buff);
   tree->Branch("pv_x_p2", "std::vector<float>", &pv_x_p2, buff);
   tree->Branch("pv_y_p2", "std::vector<float>", &pv_y_p2, buff);
   tree->Branch("pv_z_p2", "std::vector<float>", &pv_z_p2, buff);
   tree->Branch("pv_xError_p2", "std::vector<float>", &pv_xError_p2, buff);
   tree->Branch("pv_yError_p2", "std::vector<float>", &pv_yError_p2, buff);
   tree->Branch("pv_zError_p2", "std::vector<float>", &pv_zError_p2, buff);   
   
   tree->Branch("pfjet_n", &pfjet_n, "pfjet_n/I", buff);
   tree->Branch("pfjet_pt", "std::vector<float>", &pfjet_pt, buff);
   tree->Branch("pfjet_eta", "std::vector<float>", &pfjet_eta, buff);
   tree->Branch("pfjet_phi", "std::vector<float>", &pfjet_phi, buff);
   tree->Branch("pfjet_E", "std::vector<float>", &pfjet_E, buff);
   
   // Tracks from PFJets
   
/*   tree->Branch("trk_pfjet_found", "std::vector<bool>", &trk_pfjet_found, buff);
   
   tree->Branch("trk_pfjet_pt", "std::vector<float>", &trk_pfjet_pt, buff);
   tree->Branch("trk_pfjet_eta", "std::vector<float>", &trk_pfjet_eta, buff);
   tree->Branch("trk_pfjet_phi", "std::vector<float>", &trk_pfjet_phi, buff);
   tree->Branch("trk_pfjet_nTracks", "std::vector<int>", &trk_pfjet_nTracks, buff);

   tree->Branch("trk_pfjetTrk_found", "std::vector<bool>", &trk_pfjetTrk_found, buff);
   
   tree->Branch("trk_pfjetTrk_deltaR", "std::vector<float>", &trk_pfjetTrk_deltaR, buff);
   
   tree->Branch("trk_pfjetTrk_pt", "std::vector<float>", &trk_pfjetTrk_pt, buff);
   tree->Branch("trk_pfjetTrk_px", "std::vector<float>", &trk_pfjetTrk_px, buff);
   tree->Branch("trk_pfjetTrk_py", "std::vector<float>", &trk_pfjetTrk_py, buff);
   tree->Branch("trk_pfjetTrk_pz", "std::vector<float>", &trk_pfjetTrk_pz, buff);
   tree->Branch("trk_pfjetTrk_p", "std::vector<float>", &trk_pfjetTrk_p, buff);
   tree->Branch("trk_pfjetTrk_eta", "std::vector<float>", &trk_pfjetTrk_eta, buff);
   tree->Branch("trk_pfjetTrk_phi", "std::vector<float>", &trk_pfjetTrk_phi, buff);
   
   tree->Branch("trk_pfjetTrk_nTrackerLayers", "std::vector<int>", &trk_pfjetTrk_nTrackerLayers, buff);
   tree->Branch("trk_pfjetTrk_nPixelBarrelLayers", "std::vector<int>", &trk_pfjetTrk_nPixelBarrelLayers, buff);
   tree->Branch("trk_pfjetTrk_nPixelEndcapLayers", "std::vector<int>", &trk_pfjetTrk_nPixelEndcapLayers, buff);
   tree->Branch("trk_pfjetTrk_nStripLayers", "std::vector<int>", &trk_pfjetTrk_nStripLayers, buff);
   
   tree->Branch("trk_pfjetTrk_nValid", "std::vector<int>", &trk_pfjetTrk_nValid, buff);
   tree->Branch("trk_pfjetTrk_fValid", "std::vector<float>", &trk_pfjetTrk_fValid, buff);
   tree->Branch("trk_pfjetTrk_nValidTracker", "std::vector<int>", &trk_pfjetTrk_nValidTracker, buff);
   tree->Branch("trk_pfjetTrk_nValidPixelBarrel", "std::vector<int>", &trk_pfjetTrk_nValidPixelBarrel, buff);
   tree->Branch("trk_pfjetTrk_nValidPixelEndcap", "std::vector<int>", &trk_pfjetTrk_nValidPixelEndcap, buff);
   tree->Branch("trk_pfjetTrk_nValidStrip", "std::vector<int>", &trk_pfjetTrk_nValidStrip, buff);
   
   tree->Branch("trk_pfjetTrk_nMissed", "std::vector<int>", &trk_pfjetTrk_nMissed, buff);
   tree->Branch("trk_pfjetTrk_nMissedOut", "std::vector<int>", &trk_pfjetTrk_nMissedOut, buff);
   tree->Branch("trk_pfjetTrk_nMissedIn", "std::vector<int>", &trk_pfjetTrk_nMissedIn, buff);
   tree->Branch("trk_pfjetTrk_nMissedTrackerOut", "std::vector<int>", &trk_pfjetTrk_nMissedTrackerOut, buff);
   tree->Branch("trk_pfjetTrk_nMissedTrackerIn", "std::vector<int>", &trk_pfjetTrk_nMissedTrackerIn, buff);
   tree->Branch("trk_pfjetTrknMissedPixelBarrelOut", "std::vector<int>", &trk_pfjetTrk_nMissedPixelBarrelOut, buff);
   tree->Branch("trk_pfjetTrknMissedPixelBarrelIn", "std::vector<int>", &trk_pfjetTrk_nMissedPixelBarrelIn, buff);
   tree->Branch("trk_pfjetTrknMissedPixelEndcapOut", "std::vector<int>", &trk_pfjetTrk_nMissedPixelEndcapOut, buff);
   tree->Branch("trk_pfjetTrknMissedPixelEndcapIn", "std::vector<int>", &trk_pfjetTrk_nMissedPixelEndcapIn, buff);
   
   tree->Branch("trk_pfjetTrk_hasPixelBarrelLayer1", "std::vector<bool>", &trk_pfjetTrk_hasPixelBarrelLayer1, buff);
   tree->Branch("trk_pfjetTrk_hasPixelEndcapLayer1", "std::vector<bool>", &trk_pfjetTrk_hasPixelEndcapLayer1, buff);
   tree->Branch("trk_pfjetTrk_hasPixelBarrelLayer2", "std::vector<bool>", &trk_pfjetTrk_hasPixelBarrelLayer2, buff);
   tree->Branch("trk_pfjetTrk_hasPixelEndcapLayer2", "std::vector<bool>", &trk_pfjetTrk_hasPixelEndcapLayer2, buff);
   tree->Branch("trk_pfjetTrk_hasPixelBarrelLayer3", "std::vector<bool>", &trk_pfjetTrk_hasPixelBarrelLayer3, buff);
   tree->Branch("trk_pfjetTrk_hasPixelEndcapLayer3", "std::vector<bool>", &trk_pfjetTrk_hasPixelEndcapLayer3, buff);
   tree->Branch("trk_pfjetTrk_hasPixelBarrelLayer4", "std::vector<bool>", &trk_pfjetTrk_hasPixelBarrelLayer4, buff);
   tree->Branch("trk_pfjetTrk_hasPixelEndcapLayer4", "std::vector<bool>", &trk_pfjetTrk_hasPixelEndcapLayer4, buff);
   
   tree->Branch("trk_pfjetTrk_quality", "std::vector<int>", &trk_pfjetTrk_quality, buff);
   tree->Branch("trk_pfjetTrk_normalizedChi2", "std::vector<float>", &trk_pfjetTrk_normalizedChi2, buff);
   tree->Branch("trk_pfjetTrk_ndof", "std::vector<int>", &trk_pfjetTrk_ndof, buff);
   tree->Branch("trk_pfjetTrk_charge", "std::vector<int>", &trk_pfjetTrk_charge, buff);
   tree->Branch("trk_pfjetTrk_qoverp", "std::vector<float>", &trk_pfjetTrk_qoverp, buff);
   tree->Branch("trk_pfjetTrk_qoverpError", "std::vector<float>", &trk_pfjetTrk_qoverpError, buff);
   tree->Branch("trk_pfjetTrk_theta", "std::vector<float>", &trk_pfjetTrk_theta, buff);
   tree->Branch("trk_pfjetTrk_thetaError", "std::vector<float>", &trk_pfjetTrk_thetaError, buff);
   tree->Branch("trk_pfjetTrk_lambda", "std::vector<float>", &trk_pfjetTrk_lambda, buff);
   tree->Branch("trk_pfjetTrk_lambdaError", "std::vector<float>", &trk_pfjetTrk_lambdaError, buff);
   tree->Branch("trk_pfjetTrk_ptError", "std::vector<float>", &trk_pfjetTrk_ptError, buff);
   tree->Branch("trk_pfjetTrk_etaError", "std::vector<float>", &trk_pfjetTrk_etaError, buff);
   tree->Branch("trk_pfjetTrk_phiError", "std::vector<float>", &trk_pfjetTrk_phiError, buff);
   
   tree->Branch("trk_pfjetTrk_d0", "std::vector<float>", &trk_pfjetTrk_d0, buff);
   tree->Branch("trk_pfjetTrk_dz", "std::vector<float>", &trk_pfjetTrk_dz, buff);
   tree->Branch("trk_pfjetTrk_d0_pv", "std::vector<float>", &trk_pfjetTrk_d0_pv, buff);
   tree->Branch("trk_pfjetTrk_dz_pv", "std::vector<float>", &trk_pfjetTrk_dz_pv, buff);
   tree->Branch("trk_pfjetTrk_d0_bs", "std::vector<float>", &trk_pfjetTrk_d0_bs, buff);
   tree->Branch("trk_pfjetTrk_d0_bs_zpca", "std::vector<float>", &trk_pfjetTrk_d0_bs_zpca, buff);
   tree->Branch("trk_pfjetTrk_d0_bs_zpv", "std::vector<float>", &trk_pfjetTrk_d0_bs_zpv, buff);
   tree->Branch("trk_pfjetTrk_dz_bs", "std::vector<float>", &trk_pfjetTrk_dz_bs, buff);
   tree->Branch("trk_pfjetTrk_d0Err", "std::vector<float>", &trk_pfjetTrk_d0Err, buff);
   tree->Branch("trk_pfjetTrk_dzErr", "std::vector<float>", &trk_pfjetTrk_dzErr, buff);
   tree->Branch("trk_pfjetTrk_d0_pv_NoRefit", "std::vector<float>", &trk_pfjetTrk_d0_pv_NoRefit, buff);
   tree->Branch("trk_pfjetTrk_dz_pv_NoRefit", "std::vector<float>", &trk_pfjetTrk_dz_pv_NoRefit, buff);*/
}
