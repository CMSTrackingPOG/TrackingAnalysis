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

	pv_trk_mc_dxy_pvunbiased[i].clear();
	pv_trk_mc_dz_pvunbiased[i].clear();
	
	pv_trk_mc_dxy_tp_pvunbiased[i].clear();
	pv_trk_mc_dz_tp_pvunbiased[i].clear();
	
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
	
	pv_trk_d0_tv[i].clear();
	pv_trk_dz_tv[i].clear();
     }   

   pv_trk_weight.clear();
   pv_trk_isHighPurity.clear();
   pv_trk_algo.clear();
   pv_trk_originalAlgo.clear();
   
   pv_trk_idx.clear();
   
   pv_trk_pvN.clear();
   
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
   
   pv_trk_mc_dxy_pvunbiased.clear();
   pv_trk_mc_dz_pvunbiased.clear();

   pv_trk_mc_dxy_tp_pvunbiased.clear();
   pv_trk_mc_dz_tp_pvunbiased.clear();
   
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
   
   pv_trk_d0_tv.clear();
   pv_trk_dz_tv.clear();
   
   pv_mc_hasMatch.clear();   
   int pv_mc_n = pv_mc_matchQuality.size();
   for( int i=0;i<pv_mc_n;i++ )
     {	
	pv_mc_matchQuality[i].clear();
	pv_mc_isFake[i].clear();
	pv_mc_isPrimaryVertex[i].clear();
	pv_mc_isSecondaryVertex[i].clear();
	pv_mc_isTertiaryVertex[i].clear();
	pv_mc_isSignalEvent[i].clear();
	pv_mc_isBWeakDecay[i].clear();
	pv_mc_isCWeakDecay[i].clear();
	pv_mc_isTauDecay[i].clear();
	pv_mc_isKsDecay[i].clear();
	pv_mc_isLambdaDecay[i].clear();
	pv_mc_isJpsiDecay[i].clear();
	pv_mc_isXiDecay[i].clear();
	pv_mc_isOmegaDecay[i].clear();
	pv_mc_isSigmaPlusDecay[i].clear();
	pv_mc_isSigmaMinusDecay[i].clear();
	pv_mc_isLongLivedDecay[i].clear();
	
	pv_mc_isKnownProcess[i].clear();
	pv_mc_isUndefinedProcess[i].clear();
	pv_mc_isUnknownProcess[i].clear();
	pv_mc_isPrimaryProcess[i].clear();
	pv_mc_isHadronicProcess[i].clear();
	pv_mc_isDecayProcess[i].clear();
	pv_mc_isComptonProcess[i].clear();
	pv_mc_isAnnihilationProcess[i].clear();
	pv_mc_isEIoniProcess[i].clear();
	pv_mc_isHIoniProcess[i].clear();
	pv_mc_isMuIoniProcess[i].clear();
	pv_mc_isPhotonProcess[i].clear();
	pv_mc_isMuPairProdProcess[i].clear();
	pv_mc_isConversionsProcess[i].clear();
	pv_mc_isEBremProcess[i].clear();
	pv_mc_isSynchrotronRadiationProcess[i].clear();
	pv_mc_isMuBremProcess[i].clear();
	pv_mc_isMuNuclProcess[i].clear();
	pv_mc_isUnknown[i].clear();

	pv_mc_inVolume[i].clear();
	pv_mc_x[i].clear();
	pv_mc_y[i].clear();
	pv_mc_z[i].clear();
	pv_mc_t[i].clear();
	pv_mc_nGenVtx[i].clear();
	pv_mc_nSimVtx[i].clear();
	pv_mc_nDaughterTracks[i].clear();
	pv_mc_nSourceTracks[i].clear();
     }

   pv_mc_matchQuality.clear();
   pv_mc_isFake.clear();
   pv_mc_isPrimaryVertex.clear();
   pv_mc_isSecondaryVertex.clear();
   pv_mc_isTertiaryVertex.clear();
   pv_mc_isSignalEvent.clear();
   pv_mc_isBWeakDecay.clear();
   pv_mc_isCWeakDecay.clear();
   pv_mc_isTauDecay.clear();
   pv_mc_isKsDecay.clear();
   pv_mc_isLambdaDecay.clear();
   pv_mc_isJpsiDecay.clear();
   pv_mc_isXiDecay.clear();
   pv_mc_isOmegaDecay.clear();
   pv_mc_isSigmaPlusDecay.clear();
   pv_mc_isSigmaMinusDecay.clear();
   pv_mc_isLongLivedDecay.clear();
   
   pv_mc_isKnownProcess.clear();
   pv_mc_isUndefinedProcess.clear();
   pv_mc_isUnknownProcess.clear();
   pv_mc_isPrimaryProcess.clear();
   pv_mc_isHadronicProcess.clear();
   pv_mc_isDecayProcess.clear();
   pv_mc_isComptonProcess.clear();
   pv_mc_isAnnihilationProcess.clear();
   pv_mc_isEIoniProcess.clear();
   pv_mc_isHIoniProcess.clear();
   pv_mc_isMuIoniProcess.clear();
   pv_mc_isPhotonProcess.clear();
   pv_mc_isMuPairProdProcess.clear();
   pv_mc_isConversionsProcess.clear();
   pv_mc_isEBremProcess.clear();
   pv_mc_isSynchrotronRadiationProcess.clear();
   pv_mc_isMuBremProcess.clear();
   pv_mc_isMuNuclProcess.clear();
   pv_mc_isUnknown.clear();

   pv_mc_inVolume.clear();
   pv_mc_x.clear();
   pv_mc_y.clear();
   pv_mc_z.clear();
   pv_mc_t.clear();
   pv_mc_nGenVtx.clear();
   pv_mc_nSimVtx.clear();
   pv_mc_nDaughterTracks.clear();
   pv_mc_nSourceTracks.clear();
   
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

   trk_mc_hasMatch.clear();
   int trk_mc_n = trk_mc_matchQuality.size();
   for( int i=0;i<trk_mc_n;i++ )
     {	
	trk_mc_matchQuality[i].clear();
	
	trk_mc_pdgId[i].clear();
	trk_mc_origin[i].clear();
	trk_mc_status[i].clear();
	
	trk_mc_pt[i].clear();
	trk_mc_px[i].clear();
	trk_mc_py[i].clear();
	trk_mc_pz[i].clear();
	trk_mc_E[i].clear();
	trk_mc_p[i].clear();
	trk_mc_eta[i].clear();
	trk_mc_phi[i].clear();
	
	trk_mc_numberOfHits[i].clear();
	trk_mc_numberOfTrackerHits[i].clear();
	trk_mc_numberOfTrackerLayers[i].clear();

	trk_mc_dxy_center[i].clear();
	trk_mc_dz_center[i].clear();
	trk_mc_dxy_pv[i].clear();
	trk_mc_dz_pv[i].clear();
	trk_mc_dxy_bs[i].clear();
	trk_mc_dz_bs[i].clear();
	
	trk_mc_dxy_tp_center[i].clear();
	trk_mc_dz_tp_center[i].clear();
	trk_mc_dxy_tp_pv[i].clear();
	trk_mc_dz_tp_pv[i].clear();
	trk_mc_dxy_tp_bs[i].clear();
	trk_mc_dz_tp_bs[i].clear();
	
	trk_mc_vtx_x[i].clear();
	trk_mc_vtx_y[i].clear();
	trk_mc_vtx_z[i].clear();
	trk_mc_vtx_pca_x[i].clear();
	trk_mc_vtx_pca_y[i].clear();
	trk_mc_vtx_pca_z[i].clear();
	
	trk_mc_isFake[i].clear();
	trk_mc_isBad[i].clear();
	trk_mc_isBadInnerHits[i].clear();
	trk_mc_isSharedInnerHits[i].clear();
	trk_mc_isSignalEvent[i].clear();
	trk_mc_isTrackerSimHits[i].clear();
	trk_mc_isBottom[i].clear();
	trk_mc_isCharm[i].clear();
	trk_mc_isLight[i].clear();
	trk_mc_isMuon[i].clear();
	
	trk_mc_isBWeakDecay[i].clear();
	trk_mc_isCWeakDecay[i].clear();
	trk_mc_isChargePionDecay[i].clear();
	trk_mc_isChargeKaonDecay[i].clear();
	trk_mc_isTauDecay[i].clear();
	trk_mc_isKsDecay[i].clear();
	trk_mc_isLambdaDecay[i].clear();
	trk_mc_isJpsiDecay[i].clear();
	trk_mc_isXiDecay[i].clear();
	trk_mc_isOmegaDecay[i].clear();
	trk_mc_isSigmaPlusDecay[i].clear();
	trk_mc_isSigmaMinusDecay[i].clear();
	trk_mc_isLongLivedDecay[i].clear();
	
	trk_mc_isKnownProcess[i].clear();
	trk_mc_isUndefinedProcess[i].clear();
	trk_mc_isUnknownProcess[i].clear();
	trk_mc_isPrimaryProcess[i].clear();
	trk_mc_isHadronicProcess[i].clear();
	trk_mc_isDecayProcess[i].clear();
	trk_mc_isComptonProcess[i].clear();
	trk_mc_isAnnihilationProcess[i].clear();
	trk_mc_isEIoniProcess[i].clear();
	trk_mc_isHIoniProcess[i].clear();
	trk_mc_isMuIoniProcess[i].clear();
	trk_mc_isPhotonProcess[i].clear();
	trk_mc_isMuPairProdProcess[i].clear();
	trk_mc_isConversionsProcess[i].clear();
	trk_mc_isEBremProcess[i].clear();
	trk_mc_isSynchrotronRadiationProcess[i].clear();
	trk_mc_isMuBremProcess[i].clear();
	trk_mc_isMuNuclProcess[i].clear();
	
	trk_mc_isFromBWeakDecayMuon[i].clear();
	trk_mc_isFromCWeakDecayMuon[i].clear();
	trk_mc_isDecayOnFlightMuon[i].clear();
	trk_mc_isFromChargePionMuon[i].clear();
	trk_mc_isFromChargeKaonMuon[i].clear();
	
	trk_mc_isPrimaryVertex[i].clear();
	trk_mc_isSecondaryVertex[i].clear();
	trk_mc_isTertiaryVertex[i].clear();
	
	trk_mc_isUnknown[i].clear();
     }   
   
   trk_mc_matchQuality.clear();
   
   trk_mc_pdgId.clear();
   trk_mc_origin.clear();
   trk_mc_status.clear();
   
   trk_mc_pt.clear();
   trk_mc_px.clear();
   trk_mc_py.clear();
   trk_mc_pz.clear();
   trk_mc_E.clear();
   trk_mc_p.clear();
   trk_mc_eta.clear();
   trk_mc_phi.clear();

   trk_mc_numberOfHits.clear();
   trk_mc_numberOfTrackerHits.clear();
   trk_mc_numberOfTrackerLayers.clear();

   trk_mc_dxy_center.clear();
   trk_mc_dz_center.clear();
   trk_mc_dxy_pv.clear();
   trk_mc_dz_pv.clear();
   trk_mc_dxy_bs.clear();
   trk_mc_dz_bs.clear();
   
   trk_mc_dxy_tp_center.clear();
   trk_mc_dz_tp_center.clear();
   trk_mc_dxy_tp_pv.clear();
   trk_mc_dz_tp_pv.clear();
   trk_mc_dxy_tp_bs.clear();
   trk_mc_dz_tp_bs.clear();
   
   trk_mc_vtx_x.clear();
   trk_mc_vtx_y.clear();
   trk_mc_vtx_z.clear();
   trk_mc_vtx_pca_x.clear();
   trk_mc_vtx_pca_y.clear();
   trk_mc_vtx_pca_z.clear();
   
   trk_mc_isFake.clear();
   trk_mc_isBad.clear();
   trk_mc_isBadInnerHits.clear();
   trk_mc_isSharedInnerHits.clear();
   trk_mc_isSignalEvent.clear();
   trk_mc_isTrackerSimHits.clear();
   trk_mc_isBottom.clear();
   trk_mc_isCharm.clear();
   trk_mc_isLight.clear();
   trk_mc_isMuon.clear();
   
   trk_mc_isBWeakDecay.clear();
   trk_mc_isCWeakDecay.clear();
   trk_mc_isChargePionDecay.clear();
   trk_mc_isChargeKaonDecay.clear();
   trk_mc_isTauDecay.clear();
   trk_mc_isKsDecay.clear();
   trk_mc_isLambdaDecay.clear();
   trk_mc_isJpsiDecay.clear();
   trk_mc_isXiDecay.clear();
   trk_mc_isOmegaDecay.clear();
   trk_mc_isSigmaPlusDecay.clear();
   trk_mc_isSigmaMinusDecay.clear();
   trk_mc_isLongLivedDecay.clear();
   
   trk_mc_isKnownProcess.clear();
   trk_mc_isUndefinedProcess.clear();
   trk_mc_isUnknownProcess.clear();
   trk_mc_isPrimaryProcess.clear();
   trk_mc_isHadronicProcess.clear();
   trk_mc_isDecayProcess.clear();
   trk_mc_isComptonProcess.clear();
   trk_mc_isAnnihilationProcess.clear();
   trk_mc_isEIoniProcess.clear();
   trk_mc_isHIoniProcess.clear();
   trk_mc_isMuIoniProcess.clear();
   trk_mc_isPhotonProcess.clear();
   trk_mc_isMuPairProdProcess.clear();
   trk_mc_isConversionsProcess.clear();
   trk_mc_isEBremProcess.clear();
   trk_mc_isSynchrotronRadiationProcess.clear();
   trk_mc_isMuBremProcess.clear();
   trk_mc_isMuNuclProcess.clear();
   
   trk_mc_isFromBWeakDecayMuon.clear();
   trk_mc_isFromCWeakDecayMuon.clear();
   trk_mc_isDecayOnFlightMuon.clear();
   trk_mc_isFromChargePionMuon.clear();
   trk_mc_isFromChargeKaonMuon.clear();
   
   trk_mc_isPrimaryVertex.clear();
   trk_mc_isSecondaryVertex.clear();
   trk_mc_isTertiaryVertex.clear();
   
   trk_mc_isUnknown.clear();

   trk_pt.clear();
   trk_px.clear();
   trk_py.clear();
   trk_pz.clear();
   trk_p.clear();
   trk_eta.clear();
   trk_phi.clear();
   
   trk_idx.clear();
   
   trk_nTrackerLayers.clear();
   trk_nPixelBarrelLayers.clear();
   trk_nPixelEndcapLayers.clear();
   trk_nStripLayers.clear();
   
   trk_nValid.clear();
   trk_fValid.clear();
   trk_nValidTracker.clear();
   trk_nValidPixelBarrel.clear();
   trk_nValidPixelEndcap.clear();
   trk_nValidStrip.clear();

   trk_nMissed.clear();
   trk_nMissedOut.clear();
   trk_nMissedIn.clear();
   trk_nMissedTrackerOut.clear();
   trk_nMissedTrackerIn.clear();
   trk_nMissedPixelBarrelOut.clear();
   trk_nMissedPixelBarrelIn.clear();
   trk_nMissedPixelEndcapOut.clear();
   trk_nMissedPixelEndcapIn.clear();
   
   trk_hasPixelBarrelLayer1.clear();
   trk_hasPixelEndcapLayer1.clear();
   trk_hasPixelBarrelLayer2.clear();
   trk_hasPixelEndcapLayer2.clear();
   trk_hasPixelBarrelLayer3.clear();
   trk_hasPixelEndcapLayer3.clear();
   trk_hasPixelBarrelLayer4.clear();
   trk_hasPixelEndcapLayer4.clear();
   
   trk_quality.clear();
   trk_isHighPurity.clear();
   trk_normalizedChi2.clear();
   trk_ndof.clear();
   trk_charge.clear();
   trk_qoverp.clear();
   trk_qoverpError.clear();
   trk_theta.clear();
   trk_thetaError.clear();
   trk_lambda.clear();
   trk_lambdaError.clear();
   trk_ptError.clear();
   trk_etaError.clear();
   trk_phiError.clear();
   
   trk_d0.clear();
   trk_dz.clear();
   trk_d0_pv.clear();
   trk_dz_pv.clear();
   trk_d0_bs.clear();
   trk_d0_bs_zpca.clear();
   trk_d0_bs_zpv.clear();
   trk_dz_bs.clear();
   trk_d0Err.clear();
   trk_dzErr.clear();
   trk_d0_pv_NoRefit.clear();
   trk_dz_pv_NoRefit.clear();
   
   trk_d0_tv.clear();
   trk_dz_tv.clear();

   // Tracks from TrackJets
   
   trk_jet_found.clear();
   
   trk_jet_pt.clear();
   trk_jet_eta.clear();
   trk_jet_phi.clear();
   trk_jet_nTracks.clear();
   
   trk_jet_pv_x.clear();
   trk_jet_pv_y.clear();
   trk_jet_pv_z.clear();

   trk_jetTrk_found.clear();
   
   trk_jetTrk_deltaR.clear();
   
   trk_jetTrk_pt.clear();
   trk_jetTrk_px.clear();
   trk_jetTrk_py.clear();
   trk_jetTrk_pz.clear();
   trk_jetTrk_p.clear();
   trk_jetTrk_eta.clear();
   trk_jetTrk_phi.clear();
   
   trk_jetTrk_nTrackerLayers.clear();
   trk_jetTrk_nPixelBarrelLayers.clear();
   trk_jetTrk_nPixelEndcapLayers.clear();
   trk_jetTrk_nStripLayers.clear();
   
   trk_jetTrk_nValid.clear();
   trk_jetTrk_fValid.clear();
   trk_jetTrk_nValidTracker.clear();
   trk_jetTrk_nValidPixelBarrel.clear();
   trk_jetTrk_nValidPixelEndcap.clear();
   trk_jetTrk_nValidStrip.clear();
   
   trk_jetTrk_nMissed.clear();
   trk_jetTrk_nMissedOut.clear();
   trk_jetTrk_nMissedIn.clear();
   trk_jetTrk_nMissedTrackerOut.clear();
   trk_jetTrk_nMissedTrackerIn.clear();
   trk_jetTrk_nMissedPixelBarrelOut.clear();
   trk_jetTrk_nMissedPixelBarrelIn.clear();
   trk_jetTrk_nMissedPixelEndcapOut.clear();
   trk_jetTrk_nMissedPixelEndcapIn.clear();
   
   trk_jetTrk_hasPixelBarrelLayer1.clear();
   trk_jetTrk_hasPixelEndcapLayer1.clear();
   trk_jetTrk_hasPixelBarrelLayer2.clear();
   trk_jetTrk_hasPixelEndcapLayer2.clear();
   trk_jetTrk_hasPixelBarrelLayer3.clear();
   trk_jetTrk_hasPixelEndcapLayer3.clear();
   trk_jetTrk_hasPixelBarrelLayer4.clear();
   trk_jetTrk_hasPixelEndcapLayer4.clear();
   
   trk_jetTrk_quality.clear();
   trk_jetTrk_normalizedChi2.clear();
   trk_jetTrk_ndof.clear();
   trk_jetTrk_charge.clear();
   trk_jetTrk_qoverp.clear();
   trk_jetTrk_qoverpError.clear();
   trk_jetTrk_theta.clear();
   trk_jetTrk_thetaError.clear();
   trk_jetTrk_lambda.clear();
   trk_jetTrk_lambdaError.clear();
   trk_jetTrk_ptError.clear();
   trk_jetTrk_etaError.clear();
   trk_jetTrk_phiError.clear();
   
   trk_jetTrk_d0.clear();
   trk_jetTrk_dz.clear();
   trk_jetTrk_d0_pv.clear();
   trk_jetTrk_dz_pv.clear();
   trk_jetTrk_d0_bs.clear();
   trk_jetTrk_d0_bs_zpca.clear();
   trk_jetTrk_d0_bs_zpv.clear();
   trk_jetTrk_dz_bs.clear();
   trk_jetTrk_d0Err.clear();
   trk_jetTrk_dzErr.clear();
   trk_jetTrk_d0_pv_NoRefit.clear();
   trk_jetTrk_dz_pv_NoRefit.clear();

   // Tracks from PFJets
   
   trk_pfjet_found.clear();
   
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
   trk_pfjetTrk_dz_pv_NoRefit.clear();   
}

void ResTree::CreateBranches(int buff = 32000, bool runOnData = false)
{
   tree->Branch("ev_run", &ev_run, "ev_run/I", buff);
   tree->Branch("ev_id", &ev_id, "ev_id/I", buff);
   tree->Branch("ev_lumi", &ev_lumi, "ev_lumi/I", buff);
   tree->Branch("ev_bunchCrossing", &ev_bunchCrossing, "ev_bunchCrossing/I", buff);
   tree->Branch("ev_orbitNumber", &ev_orbitNumber, "ev_orbitNumber/I", buff);
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
   
   tree->Branch("pv_trk_mc_dxy_pvunbiased", "std::vector<std::vector<float> >", &pv_trk_mc_dxy_pvunbiased, buff);
   tree->Branch("pv_trk_mc_dz_pvunbiased", "std::vector<std::vector<float> >", &pv_trk_mc_dz_pvunbiased, buff);
   
   tree->Branch("pv_trk_mc_dxy_tp_pvunbiased", "std::vector<std::vector<float> >", &pv_trk_mc_dxy_tp_pvunbiased, buff);
   tree->Branch("pv_trk_mc_dz_tp_pvunbiased", "std::vector<std::vector<float> >", &pv_trk_mc_dz_tp_pvunbiased, buff);
   
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
   
   tree->Branch("pv_trk_d0_tv", "std::vector<std::vector<float> >", &pv_trk_d0_tv, buff);
   tree->Branch("pv_trk_dz_tv", "std::vector<std::vector<float> >", &pv_trk_dz_tv, buff);
   
   tree->Branch("pv_mc_hasMatch", "std::vector<bool>", &pv_mc_hasMatch, buff);
   tree->Branch("pv_mc_matchQuality", "std::vector<std::vector<float> >", &pv_mc_matchQuality, buff);
   
   tree->Branch("pv_mc_isFake", "std::vector<std::vector<bool> >", &pv_mc_isFake, buff);
   tree->Branch("pv_mc_isPrimaryVertex", "std::vector<std::vector<bool> >", &pv_mc_isPrimaryVertex, buff);
   tree->Branch("pv_mc_isSecondaryVertex", "std::vector<std::vector<bool> >", &pv_mc_isSecondaryVertex, buff);
   tree->Branch("pv_mc_isTertiaryVertex", "std::vector<std::vector<bool> >", &pv_mc_isTertiaryVertex, buff);
   tree->Branch("pv_mc_isSignalEvent", "std::vector<std::vector<bool> >", &pv_mc_isSignalEvent, buff);
   tree->Branch("pv_mc_isBWeakDecay", "std::vector<std::vector<bool> >", &pv_mc_isBWeakDecay, buff);
   tree->Branch("pv_mc_isCWeakDecay", "std::vector<std::vector<bool> >", &pv_mc_isCWeakDecay, buff);
   tree->Branch("pv_mc_isTauDecay", "std::vector<std::vector<bool> >", &pv_mc_isTauDecay, buff);
   tree->Branch("pv_mc_isKsDecay", "std::vector<std::vector<bool> >", &pv_mc_isKsDecay, buff);
   tree->Branch("pv_mc_isLambdaDecay", "std::vector<std::vector<bool> >", &pv_mc_isLambdaDecay, buff);
   tree->Branch("pv_mc_isJpsiDecay", "std::vector<std::vector<bool> >", &pv_mc_isJpsiDecay, buff);
   tree->Branch("pv_mc_isXiDecay", "std::vector<std::vector<bool> >", &pv_mc_isXiDecay, buff);
   tree->Branch("pv_mc_isOmegaDecay", "std::vector<std::vector<bool> >", &pv_mc_isOmegaDecay, buff);
   tree->Branch("pv_mc_isSigmaPlusDecay", "std::vector<std::vector<bool> >", &pv_mc_isSigmaPlusDecay, buff);
   tree->Branch("pv_mc_isSigmaMinusDecay", "std::vector<std::vector<bool> >", &pv_mc_isSigmaMinusDecay, buff);
   tree->Branch("pv_mc_isLongLivedDecay", "std::vector<std::vector<bool> >", &pv_mc_isLongLivedDecay, buff);
   
   tree->Branch("pv_mc_isKnownProcess", "std::vector<std::vector<bool> >", &pv_mc_isKnownProcess, buff);
   tree->Branch("pv_mc_isUndefinedProcess", "std::vector<std::vector<bool> >", &pv_mc_isUndefinedProcess, buff);
   tree->Branch("pv_mc_isUnknownProcess", "std::vector<std::vector<bool> >", &pv_mc_isUnknownProcess, buff);
   tree->Branch("pv_mc_isPrimaryProcess", "std::vector<std::vector<bool> >", &pv_mc_isPrimaryProcess, buff);
   tree->Branch("pv_mc_isHadronicProcess", "std::vector<std::vector<bool> >", &pv_mc_isHadronicProcess, buff);
   tree->Branch("pv_mc_isDecayProcess", "std::vector<std::vector<bool> >", &pv_mc_isDecayProcess, buff);
   tree->Branch("pv_mc_isComptonProcess", "std::vector<std::vector<bool> >", &pv_mc_isComptonProcess, buff);
   tree->Branch("pv_mc_isAnnihilationProcess", "std::vector<std::vector<bool> >", &pv_mc_isAnnihilationProcess, buff);
   tree->Branch("pv_mc_isEIoniProcess", "std::vector<std::vector<bool> >", &pv_mc_isEIoniProcess, buff);
   tree->Branch("pv_mc_isHIoniProcess", "std::vector<std::vector<bool> >", &pv_mc_isHIoniProcess, buff);
   tree->Branch("pv_mc_isMuIoniProcess", "std::vector<std::vector<bool> >", &pv_mc_isMuIoniProcess, buff);
   tree->Branch("pv_mc_isPhotonProcess", "std::vector<std::vector<bool> >", &pv_mc_isPhotonProcess, buff);
   tree->Branch("pv_mc_isMuPairProdProcess", "std::vector<std::vector<bool> >", &pv_mc_isMuPairProdProcess, buff);
   tree->Branch("pv_mc_isConversionsProcess", "std::vector<std::vector<bool> >", &pv_mc_isConversionsProcess, buff);
   tree->Branch("pv_mc_isEBremProcess", "std::vector<std::vector<bool> >", &pv_mc_isEBremProcess, buff);
   tree->Branch("pv_mc_isSynchrotronRadiationProcess", "std::vector<std::vector<bool> >", &pv_mc_isSynchrotronRadiationProcess, buff);
   tree->Branch("pv_mc_isMuBremProcess", "std::vector<std::vector<bool> >", &pv_mc_isMuBremProcess, buff);
   tree->Branch("pv_mc_isMuNuclProcess", "std::vector<std::vector<bool> >", &pv_mc_isMuNuclProcess, buff);
   tree->Branch("pv_mc_isUnknown", "std::vector<std::vector<bool> >", &pv_mc_isUnknown, buff);
   
   tree->Branch("pv_mc_inVolume", "std::vector<std::vector<bool> >", &pv_mc_inVolume, buff);
   tree->Branch("pv_mc_x", "std::vector<std::vector<float> >", &pv_mc_x, buff);
   tree->Branch("pv_mc_y", "std::vector<std::vector<float> >", &pv_mc_y, buff);
   tree->Branch("pv_mc_z", "std::vector<std::vector<float> >", &pv_mc_z, buff);
   tree->Branch("pv_mc_t", "std::vector<std::vector<float> >", &pv_mc_t, buff);
   tree->Branch("pv_mc_nGenVtx", "std::vector<std::vector<int> >", &pv_mc_nGenVtx, buff);
   tree->Branch("pv_mc_nSimVtx", "std::vector<std::vector<int> >", &pv_mc_nSimVtx, buff);
   tree->Branch("pv_mc_nDaughterTracks", "std::vector<std::vector<int> >", &pv_mc_nDaughterTracks, buff);
   tree->Branch("pv_mc_nSourceTracks", "std::vector<std::vector<int> >", &pv_mc_nSourceTracks, buff);
   
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

   tree->Branch("trk_mc_hasMatch", "std::vector<bool>", &trk_mc_hasMatch, buff);
   tree->Branch("trk_mc_matchQuality", "std::vector<std::vector<float> >", &trk_mc_matchQuality, buff);
   
   tree->Branch("trk_mc_pdgId", "std::vector<std::vector<int> >", &trk_mc_pdgId, buff);
   tree->Branch("trk_mc_origin", "std::vector<std::vector<int> >", &trk_mc_origin, buff);
   tree->Branch("trk_mc_status", "std::vector<std::vector<int> >", &trk_mc_status, buff);
   
   tree->Branch("trk_mc_pt", "std::vector<std::vector<float> >", &trk_mc_pt, buff);
   tree->Branch("trk_mc_px", "std::vector<std::vector<float> >", &trk_mc_px, buff);
   tree->Branch("trk_mc_py", "std::vector<std::vector<float> >", &trk_mc_py, buff);
   tree->Branch("trk_mc_pz", "std::vector<std::vector<float> >", &trk_mc_pz, buff);
   tree->Branch("trk_mc_E", "std::vector<std::vector<float> >", &trk_mc_E, buff);
   tree->Branch("trk_mc_p", "std::vector<std::vector<float> >", &trk_mc_p, buff);
   tree->Branch("trk_mc_eta", "std::vector<std::vector<float> >", &trk_mc_eta, buff);
   tree->Branch("trk_mc_phi", "std::vector<std::vector<float> >", &trk_mc_phi, buff);
   
   tree->Branch("trk_mc_numberOfHits", "std::vector<std::vector<int> >", &trk_mc_numberOfHits, buff);
   tree->Branch("trk_mc_numberOfTrackerHits", "std::vector<std::vector<int> >", &trk_mc_numberOfTrackerHits, buff);
   tree->Branch("trk_mc_numberOfTrackerLayers", "std::vector<std::vector<int> >", &trk_mc_numberOfTrackerLayers, buff);
   
   tree->Branch("trk_mc_dxy_center", "std::vector<std::vector<float> >", &trk_mc_dxy_center, buff);
   tree->Branch("trk_mc_dz_center", "std::vector<std::vector<float> >", &trk_mc_dz_center, buff);
   tree->Branch("trk_mc_dxy_pv", "std::vector<std::vector<float> >", &trk_mc_dxy_pv, buff);
   tree->Branch("trk_mc_dz_pv", "std::vector<std::vector<float> >", &trk_mc_dz_pv, buff);
   tree->Branch("trk_mc_dxy_bs", "std::vector<std::vector<float> >", &trk_mc_dxy_bs, buff);
   tree->Branch("trk_mc_dz_bs", "std::vector<std::vector<float> >", &trk_mc_dz_bs, buff);

   tree->Branch("trk_mc_dxy_tp_center", "std::vector<std::vector<float> >", &trk_mc_dxy_tp_center, buff);
   tree->Branch("trk_mc_dz_tp_center", "std::vector<std::vector<float> >", &trk_mc_dz_tp_center, buff);
   tree->Branch("trk_mc_dxy_tp_pv", "std::vector<std::vector<float> >", &trk_mc_dxy_tp_pv, buff);
   tree->Branch("trk_mc_dz_tp_pv", "std::vector<std::vector<float> >", &trk_mc_dz_tp_pv, buff);
   tree->Branch("trk_mc_dxy_tp_bs", "std::vector<std::vector<float> >", &trk_mc_dxy_tp_bs, buff);
   tree->Branch("trk_mc_dz_tp_bs", "std::vector<std::vector<float> >", &trk_mc_dz_tp_bs, buff);

   tree->Branch("trk_mc_vtx_x", "std::vector<std::vector<float> >", &trk_mc_vtx_x, buff);
   tree->Branch("trk_mc_vtx_y", "std::vector<std::vector<float> >", &trk_mc_vtx_y, buff);
   tree->Branch("trk_mc_vtx_z", "std::vector<std::vector<float> >", &trk_mc_vtx_z, buff);
   tree->Branch("trk_mc_vtx_pca_x", "std::vector<std::vector<float> >", &trk_mc_vtx_pca_x, buff);
   tree->Branch("trk_mc_vtx_pca_y", "std::vector<std::vector<float> >", &trk_mc_vtx_pca_y, buff);
   tree->Branch("trk_mc_vtx_pca_z", "std::vector<std::vector<float> >", &trk_mc_vtx_pca_z, buff);
   
   tree->Branch("trk_mc_isFake", "std::vector<std::vector<bool> >", &trk_mc_isFake, buff);
   tree->Branch("trk_mc_isBad", "std::vector<std::vector<bool> >", &trk_mc_isBad, buff);
   tree->Branch("trk_mc_isBadInnerHits", "std::vector<std::vector<bool> >", &trk_mc_isBadInnerHits, buff);
   tree->Branch("trk_mc_isSharedInnerHits", "std::vector<std::vector<bool> >", &trk_mc_isSharedInnerHits, buff);
   tree->Branch("trk_mc_isSignalEvent", "std::vector<std::vector<bool> >", &trk_mc_isSignalEvent, buff);
   tree->Branch("trk_mc_isTrackerSimHits", "std::vector<std::vector<bool> >", &trk_mc_isTrackerSimHits, buff);
   tree->Branch("trk_mc_isBottom", "std::vector<std::vector<bool> >", &trk_mc_isBottom, buff);
   tree->Branch("trk_mc_isCharm", "std::vector<std::vector<bool> >", &trk_mc_isCharm, buff);
   tree->Branch("trk_mc_isLight", "std::vector<std::vector<bool> >", &trk_mc_isLight, buff);
   tree->Branch("trk_mc_isMuon", "std::vector<std::vector<bool> >", &trk_mc_isMuon, buff);
   
   tree->Branch("trk_mc_isBWeakDecay", "std::vector<std::vector<bool> >", &trk_mc_isBWeakDecay, buff);
   tree->Branch("trk_mc_isCWeakDecay", "std::vector<std::vector<bool> >", &trk_mc_isCWeakDecay, buff);
   tree->Branch("trk_mc_isChargePionDecay", "std::vector<std::vector<bool> >", &trk_mc_isChargePionDecay, buff);
   tree->Branch("trk_mc_isChargeKaonDecay", "std::vector<std::vector<bool> >", &trk_mc_isChargeKaonDecay, buff);
   tree->Branch("trk_mc_isTauDecay", "std::vector<std::vector<bool> >", &trk_mc_isTauDecay, buff);
   tree->Branch("trk_mc_isKsDecay", "std::vector<std::vector<bool> >", &trk_mc_isKsDecay, buff);
   tree->Branch("trk_mc_isLambdaDecay", "std::vector<std::vector<bool> >", &trk_mc_isLambdaDecay, buff);
   tree->Branch("trk_mc_isJpsiDecay", "std::vector<std::vector<bool> >", &trk_mc_isJpsiDecay, buff);
   tree->Branch("trk_mc_isXiDecay", "std::vector<std::vector<bool> >", &trk_mc_isXiDecay, buff);
   tree->Branch("trk_mc_isOmegaDecay", "std::vector<std::vector<bool> >", &trk_mc_isOmegaDecay, buff);
   tree->Branch("trk_mc_isSigmaPlusDecay", "std::vector<std::vector<bool> >", &trk_mc_isSigmaPlusDecay, buff);
   tree->Branch("trk_mc_isSigmaMinusDecay", "std::vector<std::vector<bool> >", &trk_mc_isSigmaMinusDecay, buff);
   tree->Branch("trk_mc_isLongLivedDecay", "std::vector<std::vector<bool> >", &trk_mc_isLongLivedDecay, buff);
   
   tree->Branch("trk_mc_isKnownProcess", "std::vector<std::vector<bool> >", &trk_mc_isKnownProcess, buff);
   tree->Branch("trk_mc_isUndefinedProcess", "std::vector<std::vector<bool> >", &trk_mc_isUndefinedProcess, buff);
   tree->Branch("trk_mc_isUnknownProcess", "std::vector<std::vector<bool> >", &trk_mc_isUnknownProcess, buff);
   tree->Branch("trk_mc_isPrimaryProcess", "std::vector<std::vector<bool> >", &trk_mc_isPrimaryProcess, buff);
   tree->Branch("trk_mc_isHadronicProcess", "std::vector<std::vector<bool> >", &trk_mc_isHadronicProcess, buff);
   tree->Branch("trk_mc_isDecayProcess", "std::vector<std::vector<bool> >", &trk_mc_isDecayProcess, buff);
   tree->Branch("trk_mc_isComptonProcess", "std::vector<std::vector<bool> >", &trk_mc_isComptonProcess, buff);
   tree->Branch("trk_mc_isAnnihilationProcess", "std::vector<std::vector<bool> >", &trk_mc_isAnnihilationProcess, buff);
   tree->Branch("trk_mc_isEIoniProcess", "std::vector<std::vector<bool> >", &trk_mc_isEIoniProcess, buff);
   tree->Branch("trk_mc_isHIoniProcess", "std::vector<std::vector<bool> >", &trk_mc_isHIoniProcess, buff);
   tree->Branch("trk_mc_isMuIoniProcess", "std::vector<std::vector<bool> >", &trk_mc_isMuIoniProcess, buff);
   tree->Branch("trk_mc_isPhotonProcess", "std::vector<std::vector<bool> >", &trk_mc_isPhotonProcess, buff);
   tree->Branch("trk_mc_isMuPairProdProcess", "std::vector<std::vector<bool> >", &trk_mc_isMuPairProdProcess, buff);
   tree->Branch("trk_mc_isConversionsProcess", "std::vector<std::vector<bool> >", &trk_mc_isConversionsProcess, buff);
   tree->Branch("trk_mc_isEBremProcess", "std::vector<std::vector<bool> >", &trk_mc_isEBremProcess, buff);
   tree->Branch("trk_mc_isSynchrotronRadiationProcess", "std::vector<std::vector<bool> >", &trk_mc_isSynchrotronRadiationProcess, buff);
   tree->Branch("trk_mc_isMuBremProcess", "std::vector<std::vector<bool> >", &trk_mc_isMuBremProcess, buff);
   tree->Branch("trk_mc_isMuNuclProcess", "std::vector<std::vector<bool> >", &trk_mc_isMuNuclProcess, buff);
   
   tree->Branch("trk_mc_isFromBWeakDecayMuon", "std::vector<std::vector<bool> >", &trk_mc_isFromBWeakDecayMuon, buff);
   tree->Branch("trk_mc_isFromCWeakDecayMuon", "std::vector<std::vector<bool> >", &trk_mc_isFromCWeakDecayMuon, buff);
   tree->Branch("trk_mc_isDecayOnFlightMuon", "std::vector<std::vector<bool> >", &trk_mc_isDecayOnFlightMuon, buff);
   tree->Branch("trk_mc_isFromChargePionMuon", "std::vector<std::vector<bool> >", &trk_mc_isFromChargePionMuon, buff);
   tree->Branch("trk_mc_isFromChargeKaonMuon", "std::vector<std::vector<bool> >", &trk_mc_isFromChargeKaonMuon, buff);
   
   tree->Branch("trk_mc_isPrimaryVertex", "std::vector<std::vector<bool> >", &trk_mc_isPrimaryVertex, buff);
   tree->Branch("trk_mc_isSecondaryVertex", "std::vector<std::vector<bool> >", &trk_mc_isSecondaryVertex, buff);
   tree->Branch("trk_mc_isTertiaryVertex", "std::vector<std::vector<bool> >", &trk_mc_isTertiaryVertex, buff);
   
   tree->Branch("trk_mc_isUnknown", "std::vector<std::vector<bool> >", &trk_mc_isUnknown, buff);
   
   tree->Branch("trk_pt", "std::vector<float>", &trk_pt, buff);
   tree->Branch("trk_px", "std::vector<float>", &trk_px, buff);
   tree->Branch("trk_py", "std::vector<float>", &trk_py, buff);
   tree->Branch("trk_pz", "std::vector<float>", &trk_pz, buff);
   tree->Branch("trk_p", "std::vector<float>", &trk_p, buff);
   tree->Branch("trk_eta", "std::vector<float>", &trk_eta, buff);
   tree->Branch("trk_phi", "std::vector<float>", &trk_phi, buff);
   
   tree->Branch("trk_idx", "std::vector<int>", &trk_idx, buff);

   tree->Branch("trk_nTrackerLayers", "std::vector<int>", &trk_nTrackerLayers, buff);
   tree->Branch("trk_nPixelBarrelLayers", "std::vector<int>", &trk_nPixelBarrelLayers, buff);
   tree->Branch("trk_nPixelEndcapLayers", "std::vector<int>", &trk_nPixelEndcapLayers, buff);
   tree->Branch("trk_nStripLayers", "std::vector<int>", &trk_nStripLayers, buff);
   
   tree->Branch("trk_nValid", "std::vector<int>", &trk_nValid, buff);
   tree->Branch("trk_fValid", "std::vector<float>", &trk_fValid, buff);
   tree->Branch("trk_nValidTracker", "std::vector<int>", &trk_nValidTracker, buff);
   tree->Branch("trk_nValidPixelBarrel", "std::vector<int>", &trk_nValidPixelBarrel, buff);
   tree->Branch("trk_nValidPixelEndcap", "std::vector<int>", &trk_nValidPixelEndcap, buff);
   tree->Branch("trk_nValidStrip", "std::vector<int>", &trk_nValidStrip, buff);
   
   tree->Branch("trk_nMissed", "std::vector<int>", &trk_nMissed, buff);
   tree->Branch("trk_nMissedOut", "std::vector<int>", &trk_nMissedOut, buff);
   tree->Branch("trk_nMissedIn", "std::vector<int>", &trk_nMissedIn, buff);
   tree->Branch("trk_nMissedTrackerOut", "std::vector<int>", &trk_nMissedTrackerOut, buff);
   tree->Branch("trk_nMissedTrackerIn", "std::vector<int>", &trk_nMissedTrackerIn, buff);
   tree->Branch("trk_nMissedPixelBarrelOut", "std::vector<int>", &trk_nMissedPixelBarrelOut, buff);
   tree->Branch("trk_nMissedPixelBarrelIn", "std::vector<int>", &trk_nMissedPixelBarrelIn, buff);
   tree->Branch("trk_nMissedPixelEndcapOut", "std::vector<int>", &trk_nMissedPixelEndcapOut, buff);
   tree->Branch("trk_nMissedPixelEndcapIn", "std::vector<int>", &trk_nMissedPixelEndcapIn, buff);

   tree->Branch("trk_hasPixelBarrelLayer1", "std::vector<bool>", &trk_hasPixelBarrelLayer1, buff);
   tree->Branch("trk_hasPixelEndcapLayer1", "std::vector<bool>", &trk_hasPixelEndcapLayer1, buff);
   tree->Branch("trk_hasPixelBarrelLayer2", "std::vector<bool>", &trk_hasPixelBarrelLayer2, buff);
   tree->Branch("trk_hasPixelEndcapLayer2", "std::vector<bool>", &trk_hasPixelEndcapLayer2, buff);
   tree->Branch("trk_hasPixelBarrelLayer3", "std::vector<bool>", &trk_hasPixelBarrelLayer3, buff);
   tree->Branch("trk_hasPixelEndcapLayer3", "std::vector<bool>", &trk_hasPixelEndcapLayer3, buff);
   tree->Branch("trk_hasPixelBarrelLayer4", "std::vector<bool>", &trk_hasPixelBarrelLayer4, buff);
   tree->Branch("trk_hasPixelEndcapLayer4", "std::vector<bool>", &trk_hasPixelEndcapLayer4, buff);
   
   tree->Branch("trk_quality", "std::vector<int>", &trk_quality, buff);
   tree->Branch("trk_isHighPurity", "std::vector<bool>", &trk_isHighPurity, buff);
   tree->Branch("trk_normalizedChi2", "std::vector<float>", &trk_normalizedChi2, buff);
   tree->Branch("trk_ndof", "std::vector<int>", &trk_ndof, buff);
   tree->Branch("trk_charge", "std::vector<int>", &trk_charge, buff);
   tree->Branch("trk_qoverp", "std::vector<float>", &trk_qoverp, buff);
   tree->Branch("trk_qoverpError", "std::vector<float>", &trk_qoverpError, buff);
   tree->Branch("trk_theta", "std::vector<float>", &trk_theta, buff);
   tree->Branch("trk_thetaError", "std::vector<float>", &trk_thetaError, buff);
   tree->Branch("trk_lambda", "std::vector<float>", &trk_lambda, buff);
   tree->Branch("trk_lambdaError", "std::vector<float>", &trk_lambdaError, buff);
   tree->Branch("trk_ptError", "std::vector<float>", &trk_ptError, buff);
   tree->Branch("trk_etaError", "std::vector<float>", &trk_etaError, buff);
   tree->Branch("trk_phiError", "std::vector<float>", &trk_phiError, buff);
   
   tree->Branch("trk_d0", "std::vector<float>", &trk_d0, buff);
   tree->Branch("trk_dz", "std::vector<float>", &trk_dz, buff);
   tree->Branch("trk_d0_pv", "std::vector<float>", &trk_d0_pv, buff);
   tree->Branch("trk_dz_pv", "std::vector<float>", &trk_dz_pv, buff);
   tree->Branch("trk_d0_bs", "std::vector<float>", &trk_d0_bs, buff);
   tree->Branch("trk_d0_bs_zpca", "std::vector<float>", &trk_d0_bs_zpca, buff);
   tree->Branch("trk_d0_bs_zpv", "std::vector<float>", &trk_d0_bs_zpv, buff);
   tree->Branch("trk_dz_bs", "std::vector<float>", &trk_dz_bs, buff);
   tree->Branch("trk_d0Err", "std::vector<float>", &trk_d0Err, buff);
   tree->Branch("trk_dzErr", "std::vector<float>", &trk_dzErr, buff);
   tree->Branch("trk_d0_pv_NoRefit", "std::vector<float>", &trk_d0_pv_NoRefit, buff);
   tree->Branch("trk_dz_pv_NoRefit", "std::vector<float>", &trk_dz_pv_NoRefit, buff);
   
   tree->Branch("trk_d0_tv", "std::vector<float>", &trk_d0_tv, buff);
   tree->Branch("trk_dz_tv", "std::vector<float>", &trk_dz_tv, buff);

   // Tracks from TrackJets
   
   tree->Branch("trk_jet_found", "std::vector<bool>", &trk_jet_found, buff);
   
   tree->Branch("trk_jet_pt", "std::vector<float>", &trk_jet_pt, buff);
   tree->Branch("trk_jet_eta", "std::vector<float>", &trk_jet_eta, buff);
   tree->Branch("trk_jet_phi", "std::vector<float>", &trk_jet_phi, buff);
   tree->Branch("trk_jet_nTracks", "std::vector<int>", &trk_jet_nTracks, buff);
   
   tree->Branch("trk_jet_pv_x", "std::vector<float>", &trk_jet_pv_x, buff);
   tree->Branch("trk_jet_pv_y", "std::vector<float>", &trk_jet_pv_y, buff);
   tree->Branch("trk_jet_pv_z", "std::vector<float>", &trk_jet_pv_z, buff);
   
   tree->Branch("trk_jetTrk_found", "std::vector<bool>", &trk_jetTrk_found, buff);
   
   tree->Branch("trk_jetTrk_deltaR", "std::vector<float>", &trk_jetTrk_deltaR, buff);
   
   tree->Branch("trk_jetTrk_pt", "std::vector<float>", &trk_jetTrk_pt, buff);
   tree->Branch("trk_jetTrk_px", "std::vector<float>", &trk_jetTrk_px, buff);
   tree->Branch("trk_jetTrk_py", "std::vector<float>", &trk_jetTrk_py, buff);
   tree->Branch("trk_jetTrk_pz", "std::vector<float>", &trk_jetTrk_pz, buff);
   tree->Branch("trk_jetTrk_p", "std::vector<float>", &trk_jetTrk_p, buff);
   tree->Branch("trk_jetTrk_eta", "std::vector<float>", &trk_jetTrk_eta, buff);
   tree->Branch("trk_jetTrk_phi", "std::vector<float>", &trk_jetTrk_phi, buff);
   
   tree->Branch("trk_jetTrk_nTrackerLayers", "std::vector<int>", &trk_jetTrk_nTrackerLayers, buff);
   tree->Branch("trk_jetTrk_nPixelBarrelLayers", "std::vector<int>", &trk_jetTrk_nPixelBarrelLayers, buff);
   tree->Branch("trk_jetTrk_nPixelEndcapLayers", "std::vector<int>", &trk_jetTrk_nPixelEndcapLayers, buff);
   tree->Branch("trk_jetTrk_nStripLayers", "std::vector<int>", &trk_jetTrk_nStripLayers, buff);
   
   tree->Branch("trk_jetTrk_nValid", "std::vector<int>", &trk_jetTrk_nValid, buff);
   tree->Branch("trk_jetTrk_fValid", "std::vector<float>", &trk_jetTrk_fValid, buff);
   tree->Branch("trk_jetTrk_nValidTracker", "std::vector<int>", &trk_jetTrk_nValidTracker, buff);
   tree->Branch("trk_jetTrk_nValidPixelBarrel", "std::vector<int>", &trk_jetTrk_nValidPixelBarrel, buff);
   tree->Branch("trk_jetTrk_nValidPixelEndcap", "std::vector<int>", &trk_jetTrk_nValidPixelEndcap, buff);
   tree->Branch("trk_jetTrk_nValidStrip", "std::vector<int>", &trk_jetTrk_nValidStrip, buff);
   
   tree->Branch("trk_jetTrk_nMissed", "std::vector<int>", &trk_jetTrk_nMissed, buff);
   tree->Branch("trk_jetTrk_nMissedOut", "std::vector<int>", &trk_jetTrk_nMissedOut, buff);
   tree->Branch("trk_jetTrk_nMissedIn", "std::vector<int>", &trk_jetTrk_nMissedIn, buff);
   tree->Branch("trk_jetTrk_nMissedTrackerOut", "std::vector<int>", &trk_jetTrk_nMissedTrackerOut, buff);
   tree->Branch("trk_jetTrk_nMissedTrackerIn", "std::vector<int>", &trk_jetTrk_nMissedTrackerIn, buff);
   tree->Branch("trk_jetTrknMissedPixelBarrelOut", "std::vector<int>", &trk_jetTrk_nMissedPixelBarrelOut, buff);
   tree->Branch("trk_jetTrknMissedPixelBarrelIn", "std::vector<int>", &trk_jetTrk_nMissedPixelBarrelIn, buff);
   tree->Branch("trk_jetTrknMissedPixelEndcapOut", "std::vector<int>", &trk_jetTrk_nMissedPixelEndcapOut, buff);
   tree->Branch("trk_jetTrknMissedPixelEndcapIn", "std::vector<int>", &trk_jetTrk_nMissedPixelEndcapIn, buff);
   
   tree->Branch("trk_jetTrk_hasPixelBarrelLayer1", "std::vector<bool>", &trk_jetTrk_hasPixelBarrelLayer1, buff);
   tree->Branch("trk_jetTrk_hasPixelEndcapLayer1", "std::vector<bool>", &trk_jetTrk_hasPixelEndcapLayer1, buff);
   tree->Branch("trk_jetTrk_hasPixelBarrelLayer2", "std::vector<bool>", &trk_jetTrk_hasPixelBarrelLayer2, buff);
   tree->Branch("trk_jetTrk_hasPixelEndcapLayer2", "std::vector<bool>", &trk_jetTrk_hasPixelEndcapLayer2, buff);
   tree->Branch("trk_jetTrk_hasPixelBarrelLayer3", "std::vector<bool>", &trk_jetTrk_hasPixelBarrelLayer3, buff);
   tree->Branch("trk_jetTrk_hasPixelEndcapLayer3", "std::vector<bool>", &trk_jetTrk_hasPixelEndcapLayer3, buff);
   tree->Branch("trk_jetTrk_hasPixelBarrelLayer4", "std::vector<bool>", &trk_jetTrk_hasPixelBarrelLayer4, buff);
   tree->Branch("trk_jetTrk_hasPixelEndcapLayer4", "std::vector<bool>", &trk_jetTrk_hasPixelEndcapLayer4, buff);
   
   tree->Branch("trk_jetTrk_quality", "std::vector<int>", &trk_jetTrk_quality, buff);
   tree->Branch("trk_jetTrk_normalizedChi2", "std::vector<float>", &trk_jetTrk_normalizedChi2, buff);
   tree->Branch("trk_jetTrk_ndof", "std::vector<int>", &trk_jetTrk_ndof, buff);
   tree->Branch("trk_jetTrk_charge", "std::vector<int>", &trk_jetTrk_charge, buff);
   tree->Branch("trk_jetTrk_qoverp", "std::vector<float>", &trk_jetTrk_qoverp, buff);
   tree->Branch("trk_jetTrk_qoverpError", "std::vector<float>", &trk_jetTrk_qoverpError, buff);
   tree->Branch("trk_jetTrk_theta", "std::vector<float>", &trk_jetTrk_theta, buff);
   tree->Branch("trk_jetTrk_thetaError", "std::vector<float>", &trk_jetTrk_thetaError, buff);
   tree->Branch("trk_jetTrk_lambda", "std::vector<float>", &trk_jetTrk_lambda, buff);
   tree->Branch("trk_jetTrk_lambdaError", "std::vector<float>", &trk_jetTrk_lambdaError, buff);
   tree->Branch("trk_jetTrk_ptError", "std::vector<float>", &trk_jetTrk_ptError, buff);
   tree->Branch("trk_jetTrk_etaError", "std::vector<float>", &trk_jetTrk_etaError, buff);
   tree->Branch("trk_jetTrk_phiError", "std::vector<float>", &trk_jetTrk_phiError, buff);
   
   tree->Branch("trk_jetTrk_d0", "std::vector<float>", &trk_jetTrk_d0, buff);
   tree->Branch("trk_jetTrk_dz", "std::vector<float>", &trk_jetTrk_dz, buff);
   tree->Branch("trk_jetTrk_d0_pv", "std::vector<float>", &trk_jetTrk_d0_pv, buff);
   tree->Branch("trk_jetTrk_dz_pv", "std::vector<float>", &trk_jetTrk_dz_pv, buff);
   tree->Branch("trk_jetTrk_d0_bs", "std::vector<float>", &trk_jetTrk_d0_bs, buff);
   tree->Branch("trk_jetTrk_d0_bs_zpca", "std::vector<float>", &trk_jetTrk_d0_bs_zpca, buff);
   tree->Branch("trk_jetTrk_d0_bs_zpv", "std::vector<float>", &trk_jetTrk_d0_bs_zpv, buff);
   tree->Branch("trk_jetTrk_dz_bs", "std::vector<float>", &trk_jetTrk_dz_bs, buff);
   tree->Branch("trk_jetTrk_d0Err", "std::vector<float>", &trk_jetTrk_d0Err, buff);
   tree->Branch("trk_jetTrk_dzErr", "std::vector<float>", &trk_jetTrk_dzErr, buff);
   tree->Branch("trk_jetTrk_d0_pv_NoRefit", "std::vector<float>", &trk_jetTrk_d0_pv_NoRefit, buff);
   tree->Branch("trk_jetTrk_dz_pv_NoRefit", "std::vector<float>", &trk_jetTrk_dz_pv_NoRefit, buff);

   // Tracks from PFJets
   
   tree->Branch("trk_pfjet_found", "std::vector<bool>", &trk_pfjet_found, buff);
   
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
   tree->Branch("trk_pfjetTrk_dz_pv_NoRefit", "std::vector<float>", &trk_pfjetTrk_dz_pv_NoRefit, buff);   
}
