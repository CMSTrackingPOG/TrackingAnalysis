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
   for(unsigned int i=0;i<mc_pu_zpositions.size();i++) 
     mc_pu_zpositions[i].clear();
   mc_pu_zpositions.clear();
   for(unsigned int i=0;i<mc_pu_sumpT_lowpT.size();i++) 
     mc_pu_sumpT_lowpT[i].clear();
   mc_pu_sumpT_lowpT.clear();
   for(unsigned int i=0;i<mc_pu_sumpT_highpT.size();i++) 
     mc_pu_sumpT_highpT[i].clear();
   mc_pu_sumpT_highpT.clear();
   for(unsigned int i=0;i<mc_pu_ntrks_lowpT.size();i++) 
     mc_pu_ntrks_lowpT[i].clear();
   mc_pu_ntrks_lowpT.clear();
   for(unsigned int i=0;i<mc_pu_ntrks_highpT.size();i++) 
     mc_pu_ntrks_highpT[i].clear();
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
   
   pv_IsValid = 0;
   pv_IsFake = 0;
   pv_NTracks = null;
   pv_SumTrackPt = null;
   pv_SumTrackPt2 = null;
   pv_chi2 = null;
   pv_ndof = null;
   pv_x = null;
   pv_y = null;
   pv_z = null;
   pv_xError = null;
   pv_yError = null;
   pv_zError = null;

   pv_IsValid_p1 = 0;
   pv_IsFake_p1 = 0;
   pv_NTracks_p1 = null;
   pv_SumTrackPt_p1 = null;
   pv_SumTrackPt2_p1 = null;
   pv_chi2_p1 = null;
   pv_ndof_p1 = null;
   pv_x_p1 = null;
   pv_y_p1 = null;
   pv_z_p1 = null;
   pv_xError_p1 = null;
   pv_yError_p1 = null;
   pv_zError_p1 = null;

   pv_IsValid_p2 = 0;
   pv_IsFake_p2 = 0;
   pv_NTracks_p2 = null;
   pv_SumTrackPt_p2 = null;
   pv_SumTrackPt2_p2 = null;
   pv_chi2_p2 = null;
   pv_ndof_p2 = null;
   pv_x_p2 = null;
   pv_y_p2 = null;
   pv_z_p2 = null;
   pv_xError_p2 = null;
   pv_yError_p2 = null;
   pv_zError_p2 = null;

   trk_pt.clear();
   trk_px.clear();
   trk_py.clear();
   trk_pz.clear();
   trk_p.clear();
   trk_eta.clear();
   trk_phi.clear();
   
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
   
   tree->Branch("pv_IsValid", &pv_IsValid, "pv_IsValid/O", buff);
   tree->Branch("pv_IsFake", &pv_IsFake, "pv_IsFake/O", buff);
   tree->Branch("pv_NTracks", &pv_NTracks, "pv_NTracks/I", buff);
   tree->Branch("pv_SumTrackPt", &pv_SumTrackPt, "pv_SumTrackPt/F", buff);
   tree->Branch("pv_SumTrackPt2", &pv_SumTrackPt2, "pv_SumTrackPt2/F", buff);
   tree->Branch("pv_chi2", &pv_chi2, "pv_chi2/F", buff);
   tree->Branch("pv_ndof", &pv_ndof, "pv_ndof/I", buff);
   tree->Branch("pv_x", &pv_x, "pv_x/F", buff);
   tree->Branch("pv_y", &pv_y, "pv_y/F", buff);
   tree->Branch("pv_z", &pv_z, "pv_z/F", buff);
   tree->Branch("pv_xError", &pv_xError, "pv_xError/F", buff);
   tree->Branch("pv_yError", &pv_yError, "pv_yError/F", buff);
   tree->Branch("pv_zError", &pv_zError, "pv_zError/F", buff);

   tree->Branch("pv_IsValid_p1", &pv_IsValid_p1, "pv_IsValid_p1/O", buff);
   tree->Branch("pv_IsFake_p1", &pv_IsFake_p1, "pv_IsFake_p1/O", buff);
   tree->Branch("pv_NTracks_p1", &pv_NTracks_p1, "pv_NTracks_p1/I", buff);
   tree->Branch("pv_SumTrackPt_p1", &pv_SumTrackPt_p1, "pv_SumTrackPt_p1/F", buff);
   tree->Branch("pv_SumTrackPt2_p1", &pv_SumTrackPt2_p1, "pv_SumTrackPt2_p1/F", buff);
   tree->Branch("pv_chi2_p1", &pv_chi2_p1, "pv_chi2_p1/F", buff);
   tree->Branch("pv_ndof_p1", &pv_ndof_p1, "pv_ndof_p1/I", buff);
   tree->Branch("pv_x_p1", &pv_x_p1, "pv_x_p1/F", buff);
   tree->Branch("pv_y_p1", &pv_y_p1, "pv_y_p1/F", buff);
   tree->Branch("pv_z_p1", &pv_z_p1, "pv_z_p1/F", buff);
   tree->Branch("pv_xError_p1", &pv_xError_p1, "pv_xError_p1/F", buff);
   tree->Branch("pv_yError_p1", &pv_yError_p1, "pv_yError_p1/F", buff);
   tree->Branch("pv_zError_p1", &pv_zError_p1, "pv_zError_p1/F", buff);

   tree->Branch("pv_IsValid_p2", &pv_IsValid_p2, "pv_IsValid_p2/O", buff);
   tree->Branch("pv_IsFake_p2", &pv_IsFake_p2, "pv_IsFake_p2/O", buff);
   tree->Branch("pv_NTracks_p2", &pv_NTracks_p2, "pv_NTracks_p2/I", buff);
   tree->Branch("pv_SumTrackPt_p2", &pv_SumTrackPt_p2, "pv_SumTrackPt_p2/F", buff);
   tree->Branch("pv_SumTrackPt2_p2", &pv_SumTrackPt2_p2, "pv_SumTrackPt2_p2/F", buff);
   tree->Branch("pv_chi2_p2", &pv_chi2_p2, "pv_chi2_p2/F", buff);
   tree->Branch("pv_ndof_p2", &pv_ndof_p2, "pv_ndof_p2/I", buff);
   tree->Branch("pv_x_p2", &pv_x_p2, "pv_x_p2/F", buff);
   tree->Branch("pv_y_p2", &pv_y_p2, "pv_y_p2/F", buff);
   tree->Branch("pv_z_p2", &pv_z_p2, "pv_z_p2/F", buff);
   tree->Branch("pv_xError_p2", &pv_xError_p2, "pv_xError_p2/F", buff);
   tree->Branch("pv_yError_p2", &pv_yError_p2, "pv_yError_p2/F", buff);
   tree->Branch("pv_zError_p2", &pv_zError_p2, "pv_zError_p2/F", buff);
   
   tree->Branch("trk_pt", "std::vector<float>", &trk_pt, buff);
   tree->Branch("trk_px", "std::vector<float>", &trk_px, buff);
   tree->Branch("trk_py", "std::vector<float>", &trk_py, buff);
   tree->Branch("trk_pz", "std::vector<float>", &trk_pz, buff);
   tree->Branch("trk_p", "std::vector<float>", &trk_p, buff);
   tree->Branch("trk_eta", "std::vector<float>", &trk_eta, buff);
   tree->Branch("trk_phi", "std::vector<float>", &trk_phi, buff);
   
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

}

