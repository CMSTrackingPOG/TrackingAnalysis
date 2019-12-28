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
   ev_rho = null;

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
   pv_chi2_p2 = null;
   pv_ndof_p2 = null;
   pv_x_p2 = null;
   pv_y_p2 = null;
   pv_z_p2 = null;
   pv_xError_p2 = null;
   pv_yError_p2 = null;
   pv_zError_p2 = null;

   trk_pt.clear();
   trk_p.clear();
   trk_eta.clear();
   trk_phi.clear();
   trk_nXLayers.clear();
   trk_nMissedOut.clear();
   trk_nMissedIn.clear();
   trk_hasPXL.clear();
   trk_quality.clear();
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
}

void ResTree::CreateBranches(int buff = 32000)
{
   tree->Branch("ev_run", &ev_run, "ev_run/I", buff);
   tree->Branch("ev_id", &ev_id, "ev_id/I", buff);
   tree->Branch("ev_lumi", &ev_lumi, "ev_lumi/I", buff);
   tree->Branch("ev_rho", &ev_rho, "ev_rho/F", buff);

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
   tree->Branch("pv_chi2_p2", &pv_chi2_p2, "pv_chi2_p2/F", buff);
   tree->Branch("pv_ndof_p2", &pv_ndof_p2, "pv_ndof_p2/I", buff);
   tree->Branch("pv_x_p2", &pv_x_p2, "pv_x_p2/F", buff);
   tree->Branch("pv_y_p2", &pv_y_p2, "pv_y_p2/F", buff);
   tree->Branch("pv_z_p2", &pv_z_p2, "pv_z_p2/F", buff);
   tree->Branch("pv_xError_p2", &pv_xError_p2, "pv_xError_p2/F", buff);
   tree->Branch("pv_yError_p2", &pv_yError_p2, "pv_yError_p2/F", buff);
   tree->Branch("pv_zError_p2", &pv_zError_p2, "pv_zError_p2/F", buff);

   tree->Branch("trk_pt", "std::vector<float>", &trk_pt, buff);
   tree->Branch("trk_p", "std::vector<float>", &trk_p, buff);
   tree->Branch("trk_eta", "std::vector<float>", &trk_eta, buff);
   tree->Branch("trk_phi", "std::vector<float>", &trk_phi, buff);
   tree->Branch("trk_nXLayers", "std::vector<int>", &trk_nXLayers, buff);
   tree->Branch("trk_nMissedOut", "std::vector<int>", &trk_nMissedOut, buff);
   tree->Branch("trk_nMissedIn", "std::vector<int>", &trk_nMissedIn, buff);
   tree->Branch("trk_hasPXL", "std::vector<int>", &trk_hasPXL, buff);
   tree->Branch("trk_quality", "std::vector<int>", &trk_quality, buff);
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
}

