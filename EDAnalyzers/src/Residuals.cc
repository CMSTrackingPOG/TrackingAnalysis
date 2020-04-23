// add offline PFHT, double check what G uses for pileup reweighting
// try to access tracks from PFJets

// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "MagneticField/Engine/interface/MagneticField.h" 
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h" 

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/Common/interface/RefToPtr.h"

#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
#include "HLTrigger/HLTcore/interface/HLTPrescaleProvider.h"

#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"

#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"
 
#include "FWCore/Framework/interface/ESHandle.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/TrackReco/interface/HitPattern.h"

#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"

#include "DataFormats/JetReco/interface/TrackJet.h"
#include "DataFormats/JetReco/interface/PFJet.h"

#include <TrackingTools/TrajectoryState/interface/PerigeeConversions.h>
#include <TrackingTools/TrajectoryState/interface/TrajectoryStateClosestToPoint.h>
#include <TrackingTools/PatternTools/interface/TSCPBuilderNoMaterial.h>

#include "TrackingAnalysis/EDAnalyzers/interface/VertexReProducer.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include <CommonTools/UtilAlgos/interface/TFileService.h>

#include "TROOT.h"
#include "TH1F.h"
#include "TTree.h"
#include "TRandom3.h"
#include "Compression.h"

#include "TrackingAnalysis/EDAnalyzers/interface/Tree.h"

namespace
{
   bool sortPt(const reco::TransientTrack & t1,
	       const reco::TransientTrack & t2) 
     {
	return t1.track().pt() > t2.track().pt();
     }
}

class Residuals : public edm::EDAnalyzer 
{  
   
 public:
   explicit Residuals(const edm::ParameterSet& pset);
   ~Residuals();
    
 private:
   virtual void beginRun(const edm::Run&, const edm::EventSetup&);
   virtual void analyze(const edm::Event&, const edm::EventSetup&);
   virtual void endRun();
   
   bool trackSelection(const reco::Track& track) const;
   bool vertexSelection(const reco::Vertex& vertex) const;
   
   float getDeltaR(float eta1, float phi1, float eta2, float phi2);
   
   class TrackEqual 
     {	
      public:
	TrackEqual( const edm::Ptr<reco::Track> & t) : track_( t ) {}
	
	bool operator()( const edm::Ptr<reco::Track> & t ) const 
	  {
	     return t->pt()==track_->pt();
	  }
	
      private:
	const edm::Ptr<reco::Track> & track_;
     };

   class TrackEqualRef
     {	
      public:
	TrackEqualRef( const reco::TrackRef & t) : track_( t ) {}

	bool operator()( const reco::TrackRef & t ) const 
	  {
	     return t->pt()==track_->pt();
	  }
	
      private:
	const reco::TrackRef & track_;
     };
   
   class VertexEqual
     {	
      public:
	VertexEqual( const reco::Vertex::Point & p) : p_( p ) {}
	
	bool operator()( const reco::Vertex::Point & p ) const 
	  {
	     return (p.x()==p_.x() && p.y()==p_.y() && p.z()==p_.z());
	  }
	
      private:
	const reco::Vertex::Point & p_;
     };
   
   // ----------member data ---------------------------
   edm::EDGetTokenT<reco::VertexCollection> thePVToken_;
//   edm::EDGetTokenT<reco::VertexCollection> thePVPrimaryToken_;
   edm::EDGetTokenT<reco::TrackCollection> theTracksToken_;
   edm::EDGetTokenT<reco::BeamSpot> theBeamspotToken_;
   edm::EDGetTokenT<double> theRhoToken_;
   edm::EDGetTokenT< vector<reco::TrackJet> > theTrackJetsToken_;
   edm::EDGetTokenT< vector<reco::PFJet> > thePFJetsToken_;
   edm::EDGetTokenT<edm::TriggerResults> theTriggerBitsToken_;
   edm::EDGetTokenT<std::vector<PileupSummaryInfo> > puInfoToken_;

   // --- track selection variables
   double tkMinPt;
   int tkMinXLayers,tkMaxMissedOuterLayers,tkMaxMissedInnerLayers;

   // --- vertex selection variables
   unsigned int vtxTracksSizeMin;  
   unsigned int vtxTracksSizeMax;  
//   double vtxErrorXMin,vtxErrorXMax;
//   double vtxErrorYMin,vtxErrorYMax;
//   double vtxErrorZMin,vtxErrorZMax;

   std::string beamSpotConfig;
   
   bool runOnData;

   int eventScale;
   int trackScale;
   
   TRandom3 *rnd;
   
   HLTConfigProvider hltConfig_;
   HLTPrescaleProvider hltPrescale_;
   
   const edm::Service<TFileService> fs;
   ResTree* ftree;
   
   int ncount;
};

Residuals::Residuals(const edm::ParameterSet& pset):
   hltPrescale_(pset,consumesCollector(),*this)
{
   edm::InputTag TrackCollectionTag_ = pset.getParameter<edm::InputTag>("TrackLabel");
   theTracksToken_= consumes<reco::TrackCollection>(TrackCollectionTag_);
   
   edm::InputTag VertexCollectionTag_ = pset.getParameter<edm::InputTag>("VertexLabel");
   thePVToken_ = consumes<reco::VertexCollection>(VertexCollectionTag_);

//   edm::InputTag VertexPrimaryCollectionTag_ = pset.getParameter<edm::InputTag>("VertexPrimaryLabel");
//   thePVPrimaryToken_ = consumes<reco::VertexCollection>(VertexPrimaryCollectionTag_);
   
   edm::InputTag BeamspotTag_ = edm::InputTag("offlineBeamSpot");
   theBeamspotToken_ = consumes<reco::BeamSpot>(BeamspotTag_);

   edm::InputTag RhoTag_ = pset.getParameter<edm::InputTag>("RhoLabel");
   theRhoToken_ = consumes<double>(RhoTag_);

   edm::InputTag TrackJetsTag_ = pset.getParameter<edm::InputTag>("TrackJetsLabel");
   theTrackJetsToken_= consumes< vector<reco::TrackJet> >(TrackJetsTag_);

   edm::InputTag PFJetsTag_ = pset.getParameter<edm::InputTag>("PFJetsLabel");
   thePFJetsToken_= consumes< vector<reco::PFJet> >(PFJetsTag_);
   
   edm::InputTag TriggerBitsTag_ = pset.getParameter<edm::InputTag>("TriggerResultsLabel");
   theTriggerBitsToken_ = consumes<edm::TriggerResults>(TriggerBitsTag_);
   
   edm::InputTag PUInfoTag_ = pset.getParameter<edm::InputTag>("puInfoLabel");
   puInfoToken_ = consumes<std::vector<PileupSummaryInfo> >(PUInfoTag_);
   
   beamSpotConfig = pset.getParameter<std::string>("BeamSpotConfig");
   
   tkMinPt = pset.getParameter<double>("TkMinPt");
   tkMinXLayers = pset.getParameter<int>("TkMinXLayers");
   tkMaxMissedOuterLayers = pset.getParameter<int>("TkMaxMissedOuterLayers");
   tkMaxMissedInnerLayers = pset.getParameter<int>("TkMaxMissedInnerLayers");
   
   vtxTracksSizeMin = pset.getParameter<int>("VtxTracksSizeMin");
   vtxTracksSizeMax = pset.getParameter<int>("VtxTracksSizeMax");
//   vtxErrorXMin     = pset.getParameter<double>("VtxErrorXMin");
//   vtxErrorXMax     = pset.getParameter<double>("VtxErrorXMax");
//   vtxErrorYMin     = pset.getParameter<double>("VtxErrorYMin");
//   vtxErrorYMax     = pset.getParameter<double>("VtxErrorYMax");
//   vtxErrorZMin     = pset.getParameter<double>("VtxErrorZMin");
//   vtxErrorZMax     = pset.getParameter<double>("VtxErrorZMax");

   eventScale = pset.getParameter<int>("EventScale");
   trackScale = pset.getParameter<int>("TrackScale");
   
   runOnData = pset.getParameter<bool>("RunOnData");
   
   rnd = new TRandom3();

   TFile& f = fs->file();
   f.SetCompressionAlgorithm(ROOT::kZLIB);
   f.SetCompressionLevel(9);
   ftree = new ResTree(fs->make<TTree>("tree","tree"));   
   ftree->CreateBranches(32000,runOnData);
   
   ncount = 0;
}

Residuals::~Residuals()
{
   delete rnd;
}

void Residuals::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   ncount++;
   if( (ncount-1) % eventScale != 0 && eventScale > 0 ) return;
   
   using namespace edm;
   using namespace reco;
   using namespace std;

   ftree->Init();

   Handle<TrackCollection> tracks;
   iEvent.getByToken(theTracksToken_, tracks);
   
   Handle<VertexCollection> vtxH;
   iEvent.getByToken(thePVToken_, vtxH);

//   Handle<VertexCollection> vtxP;
//   iEvent.getByToken(thePVPrimaryToken_, vtxP);

   Handle<BeamSpot> pvbeamspot;
   iEvent.getByToken(theBeamspotToken_, pvbeamspot);
   
   if( !vtxH.isValid() ) return;

   if( vtxH->size() == 0 ) return;

   Handle< vector<reco::TrackJet> > trackJets;
   iEvent.getByToken(theTrackJetsToken_, trackJets);

   Handle< vector<reco::PFJet> > pfJets;
   iEvent.getByToken(thePFJetsToken_, pfJets);
   
   VertexReProducer revertex(vtxH, iEvent, beamSpotConfig);
   
   TrackCollection tracksAll;
   tracksAll.assign(tracks->begin(), tracks->end());

   // refit primary vertices
   vector<TransientVertex> pvs = revertex.makeVertices(tracksAll, *pvbeamspot, iSetup);   
   
   if( pvs.empty() ) return;     
   
   TransientVertex vtxTrans = pvs.front();
   reco::Vertex vtx = reco::Vertex(vtxTrans);
   
   if( !vertexSelection(vtx) ) return;

   ESHandle<MagneticField> theMF;
   iSetup.get<IdealMagneticFieldRecord>().get(theMF);

   Handle<double> rhoPtr;
   iEvent.getByToken(theRhoToken_, rhoPtr);
   
   ftree->ev_run = iEvent.id().run();
   ftree->ev_id = iEvent.id().event();
   ftree->ev_lumi = iEvent.id().luminosityBlock();
   ftree->ev_bunchCrossing = iEvent.bunchCrossing();
   ftree->ev_orbitNumber = iEvent.orbitNumber();
   ftree->ev_rho = *rhoPtr;
   ftree->ev_nPV =  pvs.size();

   Handle<TriggerResults> triggerBits;
   iEvent.getByToken(theTriggerBitsToken_, triggerBits);
   const TriggerNames &names = iEvent.triggerNames(*triggerBits);
   for( unsigned int i=0;i<names.size();i++ )
     {
	TString trigName = TString(names.triggerName(i));

//	std::cout << i << " " << trigName << std::endl;
	
	bool pass = (triggerBits->accept(i) ? true : false);
	
	if( trigName.Contains("HLT_ZeroBias_v") ) ftree->trig_ZeroBias_pass = pass;
	
	else if( trigName.Contains("HLT_ZeroBias_part0_v") ) ftree->trig_ZeroBias_part0_pass = pass;
	else if( trigName.Contains("HLT_ZeroBias_part1_v") ) ftree->trig_ZeroBias_part1_pass = pass;
	else if( trigName.Contains("HLT_ZeroBias_part2_v") ) ftree->trig_ZeroBias_part2_pass = pass;
	else if( trigName.Contains("HLT_ZeroBias_part3_v") ) ftree->trig_ZeroBias_part3_pass = pass;
	else if( trigName.Contains("HLT_ZeroBias_part4_v") ) ftree->trig_ZeroBias_part4_pass = pass;
	else if( trigName.Contains("HLT_ZeroBias_part5_v") ) ftree->trig_ZeroBias_part5_pass = pass;
	else if( trigName.Contains("HLT_ZeroBias_part6_v") ) ftree->trig_ZeroBias_part6_pass = pass;
	else if( trigName.Contains("HLT_ZeroBias_part7_v") ) ftree->trig_ZeroBias_part7_pass = pass;
	
	else if( trigName.Contains("HLT_PFJet40_v") ) ftree->trig_PFJet40_pass = pass;
	else if( trigName.Contains("HLT_PFJet60_v") ) ftree->trig_PFJet60_pass = pass;
	else if( trigName.Contains("HLT_PFJet80_v") ) ftree->trig_PFJet80_pass = pass;
	else if( trigName.Contains("HLT_PFJet140_v") ) ftree->trig_PFJet140_pass = pass;
	else if( trigName.Contains("HLT_PFJet200_v") ) ftree->trig_PFJet200_pass = pass;
	else if( trigName.Contains("HLT_PFJet260_v") ) ftree->trig_PFJet260_pass = pass;
	else if( trigName.Contains("HLT_PFJet320_v") ) ftree->trig_PFJet320_pass = pass;
	else if( trigName.Contains("HLT_PFJet400_v") ) ftree->trig_PFJet400_pass = pass;
	else if( trigName.Contains("HLT_PFJet450_v") ) ftree->trig_PFJet450_pass = pass;
	else if( trigName.Contains("HLT_PFJet500_v") ) ftree->trig_PFJet500_pass = pass;
	else if( trigName.Contains("HLT_PFJet550_v") ) ftree->trig_PFJet550_pass = pass;
	
	else if( trigName.Contains("HLT_AK4PFJet30_v") ) ftree->trig_AK4PFJet30_pass = pass;
	else if( trigName.Contains("HLT_AK4PFJet50_v") ) ftree->trig_AK4PFJet50_pass = pass;
	else if( trigName.Contains("HLT_AK4PFJet80_v") ) ftree->trig_AK4PFJet80_pass = pass;
	else if( trigName.Contains("HLT_AK4PFJet100_v") ) ftree->trig_AK4PFJet100_pass = pass;
	else if( trigName.Contains("HLT_AK4PFJet120_v") ) ftree->trig_AK4PFJet120_pass = pass;
	
	else if( trigName.Contains("HLT_PFHT180_v") ) ftree->trig_PFHT180_pass = pass;
	else if( trigName.Contains("HLT_PFHT250_v") ) ftree->trig_PFHT250_pass = pass;	
	else if( trigName.Contains("HLT_PFHT370_v") ) ftree->trig_PFHT370_pass = pass;
	else if( trigName.Contains("HLT_PFHT430_v") ) ftree->trig_PFHT430_pass = pass;
	else if( trigName.Contains("HLT_PFHT510_v") ) ftree->trig_PFHT510_pass = pass;
	else if( trigName.Contains("HLT_PFHT590_v") ) ftree->trig_PFHT590_pass = pass;
	else if( trigName.Contains("HLT_PFHT680_v") ) ftree->trig_PFHT680_pass = pass;
	else if( trigName.Contains("HLT_PFHT780_v") ) ftree->trig_PFHT780_pass = pass;
	else if( trigName.Contains("HLT_PFHT890_v") ) ftree->trig_PFHT890_pass = pass;
	else if( trigName.Contains("HLT_PFHT1050_v") ) ftree->trig_PFHT1050_pass = pass;
	else if( trigName.Contains("HLT_PFHT350_v") ) ftree->trig_PFHT350_pass = pass;
	  
//	     std::cout << i << " " << trigName << " L1=" << L1prescale << " HLT=" << HLTprescale << std::endl;
     }

   if( !runOnData )
     {	
	edm::Handle<std::vector< PileupSummaryInfo> > pileupInfo;
	iEvent.getByToken(puInfoToken_,pileupInfo);

	ftree->mc_pu_Npvi = pileupInfo->size();
	for(std::vector<PileupSummaryInfo>::const_iterator pvi=pileupInfo->begin();pvi!=pileupInfo->end();pvi++)
	  {
	     signed int n_bc = pvi->getBunchCrossing();
	     ftree->mc_pu_BunchCrossing.push_back(n_bc);
	     if( n_bc == 0 )
	       {		  
		  ftree->mc_pu_intime_NumInt = pvi->getPU_NumInteractions();
		  ftree->mc_pu_trueNumInt = pvi->getTrueNumInteractions();
	       }	     
	     else if( n_bc == -1 ) ftree->mc_pu_before_npu = pvi->getPU_NumInteractions();
	     else if( n_bc == +1 ) ftree->mc_pu_after_npu  = pvi->getPU_NumInteractions();
	     
	     std::vector<float> mc_pu_zpositions;
	     std::vector<float> mc_pu_sumpT_lowpT;
	     std::vector<float> mc_pu_sumpT_highpT;
	     std::vector<int> mc_pu_ntrks_lowpT;
	     std::vector<int> mc_pu_ntrks_highpT;
	     
	     ftree->mc_pu_Nzpositions.push_back(pvi->getPU_zpositions().size());
	     for( unsigned int ipu=0;ipu<pvi->getPU_zpositions().size();ipu++ )
	       {		  
		  mc_pu_zpositions.push_back((pvi->getPU_zpositions())[ipu]);
		  mc_pu_sumpT_lowpT.push_back((pvi->getPU_sumpT_lowpT())[ipu]);
		  mc_pu_sumpT_highpT.push_back((pvi->getPU_sumpT_highpT())[ipu]);
		  mc_pu_ntrks_lowpT.push_back((pvi->getPU_ntrks_lowpT())[ipu]);
		  mc_pu_ntrks_highpT.push_back((pvi->getPU_ntrks_highpT())[ipu]);
	       }	     
	     
	     ftree->mc_pu_zpositions.push_back(mc_pu_zpositions);
	     ftree->mc_pu_sumpT_lowpT.push_back(mc_pu_sumpT_lowpT);
	     ftree->mc_pu_sumpT_highpT.push_back(mc_pu_sumpT_highpT);
	     ftree->mc_pu_ntrks_lowpT.push_back(mc_pu_ntrks_lowpT);
	     ftree->mc_pu_ntrks_highpT.push_back(mc_pu_ntrks_highpT);
	  }	
     }   
   
   double micron = 10000;

   // Beam spot   
   ftree->bs_type = pvbeamspot->type();
   ftree->bs_x0 = pvbeamspot->x0();
   ftree->bs_y0 = pvbeamspot->y0();
   ftree->bs_z0 = pvbeamspot->z0();
   ftree->bs_x_zpv = pvbeamspot->x(vtx.z());
   ftree->bs_y_zpv = pvbeamspot->y(vtx.z());
   ftree->bs_sigmaZ = pvbeamspot->sigmaZ();
   ftree->bs_dxdz = pvbeamspot->dxdz();
   ftree->bs_dydz = pvbeamspot->dydz();
   ftree->bs_BeamWidthX = pvbeamspot->BeamWidthX();
   ftree->bs_BeamWidthY = pvbeamspot->BeamWidthY();
   ftree->bs_x0Error = pvbeamspot->x0Error();
   ftree->bs_y0Error = pvbeamspot->y0Error();
   ftree->bs_z0Error = pvbeamspot->z0Error();
   ftree->bs_sigmaZ0Error = pvbeamspot->sigmaZ0Error();
   ftree->bs_dxdzError = pvbeamspot->dxdzError();
   ftree->bs_dydzError = pvbeamspot->dydzError();
   ftree->bs_BeamWidthXError = pvbeamspot->BeamWidthXError();
   ftree->bs_BeamWidthYError = pvbeamspot->BeamWidthYError();
   ftree->bs_emittanceX = pvbeamspot->emittanceX();
   ftree->bs_emittanceY = pvbeamspot->emittanceY();
   ftree->bs_betaStar = pvbeamspot->betaStar();

   // Primary vertex
   double trackSumPt = 0;
   double trackSumPt2 = 0;
   
   std::vector<reco::TransientTrack> vtxTracks = vtxTrans.originalTracks();   
   stable_sort(vtxTracks.begin(), vtxTracks.end(), sortPt);
   
   for( std::vector<reco::TransientTrack>::const_iterator it = vtxTracks.begin(); it != vtxTracks.end(); it++ )
     {
	reco::Track trk = (*it).track();
	trackSumPt += trk.pt();
	trackSumPt2 += pow(trk.pt(),2);
     }
   
   ftree->pv_IsValid = vtx.isValid();
   ftree->pv_IsFake = vtx.isFake();
   ftree->pv_NTracks = vtx.tracksSize();
   ftree->pv_SumTrackPt = trackSumPt;
   ftree->pv_SumTrackPt2 = trackSumPt2;
   ftree->pv_chi2 = vtx.chi2();	
   ftree->pv_ndof = vtx.ndof();   
   ftree->pv_x = vtx.x()*micron;
   ftree->pv_y = vtx.y()*micron;
   ftree->pv_z = vtx.z()*micron;
   ftree->pv_xError = vtx.xError()*micron;
   ftree->pv_yError = vtx.yError()*micron;
   ftree->pv_zError = vtx.zError()*micron;

   // Vertex split method
   reco::TrackCollection vtxTkCollection1;
   reco::TrackCollection vtxTkCollection2;

   double trackSumPt_p1 = 0;
   double trackSumPt2_p1 = 0;

   double trackSumPt_p2 = 0;
   double trackSumPt2_p2 = 0;
   
   for( std::vector<reco::TransientTrack>::const_iterator it = vtxTracks.begin(); it != vtxTracks.end(); it++ )
     {
	reco::Track trk = (*it).track();

	if( rnd->Rndm() > 0.5 )
	  {
	     vtxTkCollection1.push_back(trk);
	     trackSumPt_p1 += trk.pt();
	     trackSumPt2_p1 += pow(trk.pt(),2);
	  }	
	else
	  {	     
	     vtxTkCollection2.push_back(trk);
	     trackSumPt_p2 += trk.pt();
	     trackSumPt2_p2 += pow(trk.pt(),2);
	  }
     }

   vector<TransientVertex> pvs1 = revertex.makeVertices(vtxTkCollection1, *pvbeamspot, iSetup);
   vector<TransientVertex> pvs2 = revertex.makeVertices(vtxTkCollection2, *pvbeamspot, iSetup);

   if( !pvs1.empty() && !pvs2.empty() )
     {
	reco::Vertex vtx1 = reco::Vertex(pvs1.front());
	reco::Vertex vtx2 = reco::Vertex(pvs2.front());

	ftree->pv_IsValid_p1 = vtx1.isValid();
	ftree->pv_IsValid_p2 = vtx2.isValid();
	
	ftree->pv_IsFake_p1 = vtx1.isFake();
	ftree->pv_IsFake_p2 = vtx2.isFake();
	
	ftree->pv_NTracks_p1 = vtxTkCollection1.size();
	ftree->pv_NTracks_p2 = vtxTkCollection2.size();
	
	ftree->pv_SumTrackPt_p1 = trackSumPt_p1;
	ftree->pv_SumTrackPt_p2 = trackSumPt_p2;

	ftree->pv_SumTrackPt2_p1 = trackSumPt2_p1;
	ftree->pv_SumTrackPt2_p2 = trackSumPt2_p2;
	
	ftree->pv_chi2_p1 = vtx1.chi2();
	ftree->pv_chi2_p2 = vtx2.chi2();
	
	ftree->pv_ndof_p1 = vtx1.ndof();
        ftree->pv_ndof_p2 = vtx2.ndof();
	
	ftree->pv_x_p1 = vtx1.x()*micron;
	ftree->pv_y_p1 = vtx1.y()*micron;
        ftree->pv_z_p1 = vtx1.z()*micron;
	ftree->pv_xError_p1 = vtx1.xError()*micron;
	ftree->pv_yError_p1 = vtx1.yError()*micron;
	ftree->pv_zError_p1 = vtx1.zError()*micron;
	
	ftree->pv_x_p2 = vtx2.x()*micron;
	ftree->pv_y_p2 = vtx2.y()*micron;
	ftree->pv_z_p2 = vtx2.z()*micron;
	ftree->pv_xError_p2 = vtx2.xError()*micron;
	ftree->pv_yError_p2 = vtx2.yError()*micron;
	ftree->pv_zError_p2 = vtx2.zError()*micron;
     }

   // PFJets
   
   unsigned int nFPJets = pfJets->size();
   ftree->pfjet_n = nFPJets;

   for( unsigned int ij=0;ij<nFPJets;ij++ )
     {
	reco::PFJet jet = pfJets->at(ij);

	ftree->pfjet_pt.push_back( jet.pt() );
	ftree->pfjet_eta.push_back( jet.eta() );
	ftree->pfjet_phi.push_back( jet.phi() );
	ftree->pfjet_E.push_back( jet.energy() );
     }   
   
   // Tracks
   float trackProb = 1./float(trackScale);
   int nTracks = tracks->size();
   
   std::cout << "Tracks = " << nTracks << ", use = " << int(float(nTracks)/float(trackScale)) << std::endl;
   
   for( TrackCollection::const_iterator itk = tracks->begin(); itk != tracks->end(); ++itk )
     {
	if( rnd->Rndm() > trackProb && trackScale > 0 ) continue;
	
	// --- track selection ---
	if( ! trackSelection(*itk) ) continue;
	// ---
     
	TrackCollection newTkCollection;
	newTkCollection.assign(tracks->begin(), itk);
	newTkCollection.insert(newTkCollection.end(),itk+1,tracks->end());

	//newTkCollection.insert(newTkCollection.end(),itk,tracks->end()); // only for debugging purpose

	//cout << "tracks before,after size: " << tracks->size() << " , " << newTkCollection.size() << endl;

	// Refit the primary vertex
	vector<TransientVertex> pvs = revertex.makeVertices(newTkCollection, *pvbeamspot, iSetup);
	//cout << "vertices before,after: " << vtxH->size() << " , " << pvs.size() << endl;
	
	if( pvs.empty() ) continue;

	reco::Vertex newPV = reco::Vertex(pvs.front());
	Track::Point vtxPosition = Track::Point(newPV.position().x(),
						newPV.position().y(),
						newPV.position().z());
	
	if( ! vertexSelection(newPV) ) continue;
	
	ftree->trk_pt.push_back( itk->pt() );
	ftree->trk_px.push_back( itk->px() );
	ftree->trk_py.push_back( itk->py() );
	ftree->trk_pz.push_back( itk->pz() );
	ftree->trk_p.push_back( itk->p() );
	ftree->trk_eta.push_back( itk->eta() );
	ftree->trk_phi.push_back( itk->phi() );
	
	ftree->trk_nTrackerLayers.push_back( itk->hitPattern().trackerLayersWithMeasurement() );
	ftree->trk_nPixelBarrelLayers.push_back( itk->hitPattern().pixelBarrelLayersWithMeasurement() );
	ftree->trk_nPixelEndcapLayers.push_back( itk->hitPattern().pixelEndcapLayersWithMeasurement() );
	ftree->trk_nStripLayers.push_back( itk->hitPattern().stripLayersWithMeasurement() );
	
	ftree->trk_nValid.push_back( itk->numberOfValidHits() );
	ftree->trk_fValid.push_back( itk->validFraction() );
	ftree->trk_nValidTracker.push_back( itk->hitPattern().numberOfValidTrackerHits() );
	ftree->trk_nValidPixelBarrel.push_back( itk->hitPattern().numberOfValidPixelBarrelHits() );
	ftree->trk_nValidPixelEndcap.push_back( itk->hitPattern().numberOfValidPixelEndcapHits() );
	ftree->trk_nValidStrip.push_back( itk->hitPattern().numberOfValidStripHits() );
	
	ftree->trk_nMissed.push_back( itk->numberOfLostHits() );
	ftree->trk_nMissedOut.push_back( itk->hitPattern().numberOfLostHits(HitPattern::MISSING_OUTER_HITS) );
	ftree->trk_nMissedIn.push_back( itk->hitPattern().numberOfLostHits(HitPattern::MISSING_INNER_HITS) );
	ftree->trk_nMissedTrackerOut.push_back( itk->hitPattern().numberOfLostTrackerHits(HitPattern::MISSING_OUTER_HITS) );
	ftree->trk_nMissedTrackerIn.push_back( itk->hitPattern().numberOfLostTrackerHits(HitPattern::MISSING_INNER_HITS) );
	ftree->trk_nMissedPixelBarrelOut.push_back( itk->hitPattern().numberOfLostPixelBarrelHits(HitPattern::MISSING_OUTER_HITS) );
	ftree->trk_nMissedPixelBarrelIn.push_back( itk->hitPattern().numberOfLostPixelBarrelHits(HitPattern::MISSING_INNER_HITS) );
	ftree->trk_nMissedPixelEndcapOut.push_back( itk->hitPattern().numberOfLostPixelEndcapHits(HitPattern::MISSING_OUTER_HITS) );
	ftree->trk_nMissedPixelEndcapIn.push_back( itk->hitPattern().numberOfLostPixelEndcapHits(HitPattern::MISSING_INNER_HITS) );
	
	ftree->trk_hasPixelBarrelLayer1.push_back( itk->hitPattern().hasValidHitInPixelLayer(PixelSubdetector::SubDetector::PixelBarrel, 1) );
	ftree->trk_hasPixelEndcapLayer1.push_back( itk->hitPattern().hasValidHitInPixelLayer(PixelSubdetector::SubDetector::PixelEndcap, 1) );
	ftree->trk_hasPixelBarrelLayer2.push_back( itk->hitPattern().hasValidHitInPixelLayer(PixelSubdetector::SubDetector::PixelBarrel, 2) );
	ftree->trk_hasPixelEndcapLayer2.push_back( itk->hitPattern().hasValidHitInPixelLayer(PixelSubdetector::SubDetector::PixelEndcap, 2) );
	ftree->trk_hasPixelBarrelLayer3.push_back( itk->hitPattern().hasValidHitInPixelLayer(PixelSubdetector::SubDetector::PixelBarrel, 3) );
	ftree->trk_hasPixelEndcapLayer3.push_back( itk->hitPattern().hasValidHitInPixelLayer(PixelSubdetector::SubDetector::PixelEndcap, 3) );
	ftree->trk_hasPixelBarrelLayer4.push_back( itk->hitPattern().hasValidHitInPixelLayer(PixelSubdetector::SubDetector::PixelBarrel, 4) );
	ftree->trk_hasPixelEndcapLayer4.push_back( itk->hitPattern().hasValidHitInPixelLayer(PixelSubdetector::SubDetector::PixelEndcap, 4) );

	ftree->trk_quality.push_back( itk->qualityMask() );
	ftree->trk_normalizedChi2.push_back( itk->normalizedChi2() );
	ftree->trk_ndof.push_back( itk->ndof() );
	ftree->trk_charge.push_back( itk->charge() );
	ftree->trk_qoverp.push_back( itk->qoverp() );
	ftree->trk_qoverpError.push_back( itk->qoverpError() );
	ftree->trk_theta.push_back( itk->theta() );
	ftree->trk_thetaError.push_back( itk->thetaError() );
	ftree->trk_lambda.push_back( itk->lambda() );
	ftree->trk_lambdaError.push_back( itk->lambdaError() );
	ftree->trk_ptError.push_back( itk->ptError() );
	ftree->trk_etaError.push_back( itk->etaError() );
	ftree->trk_phiError.push_back( itk->phiError() );
	
	ftree->trk_d0.push_back( itk->dxy() * micron );
	ftree->trk_dz.push_back( itk->dz() * micron );
	ftree->trk_d0_pv.push_back( itk->dxy(vtxPosition) * micron );
	ftree->trk_dz_pv.push_back( itk->dz(vtxPosition) * micron );
	ftree->trk_d0_bs.push_back( itk->dxy(pvbeamspot->position()) * micron );
	ftree->trk_d0_bs_zpca.push_back( itk->dxy(*pvbeamspot) * micron );
	ftree->trk_d0_bs_zpv.push_back( itk->dxy(pvbeamspot->position(vtx.z())) * micron );
	ftree->trk_dz_bs.push_back( itk->dz(pvbeamspot->position()) * micron );
	ftree->trk_d0Err.push_back( itk->d0Error() * micron );
	ftree->trk_dzErr.push_back( itk->dzError() * micron );
	ftree->trk_d0_pv_NoRefit.push_back( itk->dxy(vtxH->front().position()) * micron );
	ftree->trk_dz_pv_NoRefit.push_back( itk->dz(vtxH->front().position()) * micron );

	// Tracks from TrackJets
	
	float drMin = 10E+10;
	int iJet = -1;
	int iTrackMin = -1;

	for( unsigned int ij=0;ij<trackJets->size();ij++ )
	  {
	     std::vector<edm::Ptr<reco::Track> > trks = trackJets->at(ij).tracks();

	     const reco::TrackRef trkRef = reco::TrackRef(tracks, itk - tracks->begin());
	     
	     std::vector<edm::Ptr<reco::Track> >::const_iterator itt = find_if(trks.begin(), trks.end(), TrackEqual(edm::refToPtr(trkRef)));
	     if( itt != trks.end() )
	       {		  		  
		  size_t pos = itt - trks.begin();
		  
		  iJet = ij;
		  
		  for( unsigned int it=0;it<trackJets->at(ij).numberOfTracks();it++ )
		    {
		       if( it == pos ) continue;
		       
		       float dr = getDeltaR(trks[pos]->eta(), trks[pos]->phi(), trks[it]->eta(), trks[it]->phi());
		       if( dr < drMin )
			 {
			    drMin = dr;
			    iTrackMin = it;
			 }
		    }
		  
		  break;
	       }
	  }
	
	if( iJet >= 0 )
	  {
	     reco::TrackJet jet = trackJets->at(iJet);
	     
	     const reco::VertexRef vtxj = jet.primaryVertex();
	     
	     ftree->trk_jet_found.push_back( true );
	     
	     ftree->trk_jet_pt.push_back( jet.pt() );
	     ftree->trk_jet_eta.push_back( jet.eta() );
	     ftree->trk_jet_phi.push_back( jet.phi() );
	     ftree->trk_jet_nTracks.push_back( jet.numberOfTracks() );
	     
	     ftree->trk_jet_pv_x.push_back( vtxj->x() );
	     ftree->trk_jet_pv_y.push_back( vtxj->y() );
	     ftree->trk_jet_pv_z.push_back( vtxj->z() );
	  }
	else
	  {
	     ftree->trk_jet_found.push_back( false );
	     
	     ftree->trk_jet_pt.push_back( null );
	     ftree->trk_jet_eta.push_back( null );
	     ftree->trk_jet_phi.push_back( null );
	     ftree->trk_jet_nTracks.push_back( null );
	     
	     ftree->trk_jet_pv_x.push_back( null );
	     ftree->trk_jet_pv_y.push_back( null );
	     ftree->trk_jet_pv_z.push_back( null );
	  }	
	     
	if( iTrackMin >= 0 )
	  {		  
	     edm::Ptr<reco::Track> trkCl = trackJets->at(iJet).track(iTrackMin);
	     
	     ftree->trk_jetTrk_found.push_back( true );
	     
	     ftree->trk_jetTrk_deltaR.push_back( drMin );
	     
	     ftree->trk_jetTrk_pt.push_back( trkCl->pt() );
	     ftree->trk_jetTrk_px.push_back( trkCl->px() );
	     ftree->trk_jetTrk_py.push_back( trkCl->py() );
	     ftree->trk_jetTrk_pz.push_back( trkCl->pz() );
	     ftree->trk_jetTrk_p.push_back( trkCl->p() );
	     ftree->trk_jetTrk_eta.push_back( trkCl->eta() );
	     ftree->trk_jetTrk_phi.push_back( trkCl->phi() );
	     
	     ftree->trk_jetTrk_nTrackerLayers.push_back( trkCl->hitPattern().trackerLayersWithMeasurement() );
	     ftree->trk_jetTrk_nPixelBarrelLayers.push_back( trkCl->hitPattern().pixelBarrelLayersWithMeasurement() );
	     ftree->trk_jetTrk_nPixelEndcapLayers.push_back( trkCl->hitPattern().pixelEndcapLayersWithMeasurement() );
	     ftree->trk_jetTrk_nStripLayers.push_back( trkCl->hitPattern().stripLayersWithMeasurement() );
	     
	     ftree->trk_jetTrk_nValid.push_back( trkCl->numberOfValidHits() );
	     ftree->trk_jetTrk_fValid.push_back( trkCl->validFraction() );
	     ftree->trk_jetTrk_nValidTracker.push_back( trkCl->hitPattern().numberOfValidTrackerHits() );
	     ftree->trk_jetTrk_nValidPixelBarrel.push_back( trkCl->hitPattern().numberOfValidPixelBarrelHits() );
	     ftree->trk_jetTrk_nValidPixelEndcap.push_back( trkCl->hitPattern().numberOfValidPixelEndcapHits() );
	     ftree->trk_jetTrk_nValidStrip.push_back( trkCl->hitPattern().numberOfValidStripHits() );
	     
	     ftree->trk_jetTrk_nMissed.push_back( trkCl->numberOfLostHits() );
	     ftree->trk_jetTrk_nMissedOut.push_back( trkCl->hitPattern().numberOfLostHits(HitPattern::MISSING_OUTER_HITS) );
	     ftree->trk_jetTrk_nMissedIn.push_back( trkCl->hitPattern().numberOfLostHits(HitPattern::MISSING_INNER_HITS) );
	     ftree->trk_jetTrk_nMissedTrackerOut.push_back( trkCl->hitPattern().numberOfLostTrackerHits(HitPattern::MISSING_OUTER_HITS) );
	     ftree->trk_jetTrk_nMissedTrackerIn.push_back( trkCl->hitPattern().numberOfLostTrackerHits(HitPattern::MISSING_INNER_HITS) );
	     ftree->trk_jetTrk_nMissedPixelBarrelOut.push_back( trkCl->hitPattern().numberOfLostPixelBarrelHits(HitPattern::MISSING_OUTER_HITS) );
	     ftree->trk_jetTrk_nMissedPixelBarrelIn.push_back( trkCl->hitPattern().numberOfLostPixelBarrelHits(HitPattern::MISSING_INNER_HITS) );
	     ftree->trk_jetTrk_nMissedPixelEndcapOut.push_back( trkCl->hitPattern().numberOfLostPixelEndcapHits(HitPattern::MISSING_OUTER_HITS) );
	     ftree->trk_jetTrk_nMissedPixelEndcapIn.push_back( trkCl->hitPattern().numberOfLostPixelEndcapHits(HitPattern::MISSING_INNER_HITS) );
	     
	     ftree->trk_jetTrk_hasPixelBarrelLayer1.push_back( trkCl->hitPattern().hasValidHitInPixelLayer(PixelSubdetector::SubDetector::PixelBarrel, 1) );
	     ftree->trk_jetTrk_hasPixelEndcapLayer1.push_back( trkCl->hitPattern().hasValidHitInPixelLayer(PixelSubdetector::SubDetector::PixelEndcap, 1) );
	     ftree->trk_jetTrk_hasPixelBarrelLayer2.push_back( trkCl->hitPattern().hasValidHitInPixelLayer(PixelSubdetector::SubDetector::PixelBarrel, 2) );
	     ftree->trk_jetTrk_hasPixelEndcapLayer2.push_back( trkCl->hitPattern().hasValidHitInPixelLayer(PixelSubdetector::SubDetector::PixelEndcap, 2) );
	     ftree->trk_jetTrk_hasPixelBarrelLayer3.push_back( trkCl->hitPattern().hasValidHitInPixelLayer(PixelSubdetector::SubDetector::PixelBarrel, 3) );
	     ftree->trk_jetTrk_hasPixelEndcapLayer3.push_back( trkCl->hitPattern().hasValidHitInPixelLayer(PixelSubdetector::SubDetector::PixelEndcap, 3) );
	     ftree->trk_jetTrk_hasPixelBarrelLayer4.push_back( trkCl->hitPattern().hasValidHitInPixelLayer(PixelSubdetector::SubDetector::PixelBarrel, 4) );
	     ftree->trk_jetTrk_hasPixelEndcapLayer4.push_back( trkCl->hitPattern().hasValidHitInPixelLayer(PixelSubdetector::SubDetector::PixelEndcap, 4) );
	     
	     ftree->trk_jetTrk_quality.push_back( trkCl->qualityMask() );
	     ftree->trk_jetTrk_normalizedChi2.push_back( trkCl->normalizedChi2() );
	     ftree->trk_jetTrk_ndof.push_back( trkCl->ndof() );
	     ftree->trk_jetTrk_charge.push_back( trkCl->charge() );
	     ftree->trk_jetTrk_qoverp.push_back( trkCl->qoverp() );
	     ftree->trk_jetTrk_qoverpError.push_back( trkCl->qoverpError() );
	     ftree->trk_jetTrk_theta.push_back( trkCl->theta() );
	     ftree->trk_jetTrk_thetaError.push_back( trkCl->thetaError() );
	     ftree->trk_jetTrk_lambda.push_back( trkCl->lambda() );
	     ftree->trk_jetTrk_lambdaError.push_back( trkCl->lambdaError() );
	     ftree->trk_jetTrk_ptError.push_back( trkCl->ptError() );
	     ftree->trk_jetTrk_etaError.push_back( trkCl->etaError() );
	     ftree->trk_jetTrk_phiError.push_back( trkCl->phiError() );
	     
	     ftree->trk_jetTrk_d0.push_back( trkCl->dxy() * micron );
	     ftree->trk_jetTrk_dz.push_back( trkCl->dz() * micron );
	     ftree->trk_jetTrk_d0_pv.push_back( trkCl->dxy(vtxPosition) * micron );
	     ftree->trk_jetTrk_dz_pv.push_back( trkCl->dz(vtxPosition) * micron );
	     ftree->trk_jetTrk_d0_bs.push_back( trkCl->dxy(pvbeamspot->position()) * micron );
	     ftree->trk_jetTrk_d0_bs_zpca.push_back( trkCl->dxy(*pvbeamspot) * micron );
	     ftree->trk_jetTrk_d0_bs_zpv.push_back( trkCl->dxy(pvbeamspot->position(vtx.z())) * micron );
	     ftree->trk_jetTrk_dz_bs.push_back( trkCl->dz(pvbeamspot->position()) * micron );
	     ftree->trk_jetTrk_d0Err.push_back( trkCl->d0Error() * micron );
	     ftree->trk_jetTrk_dzErr.push_back( trkCl->dzError() * micron );
	     ftree->trk_jetTrk_d0_pv_NoRefit.push_back( trkCl->dxy(vtxH->front().position()) * micron );
	     ftree->trk_jetTrk_dz_pv_NoRefit.push_back( trkCl->dz(vtxH->front().position()) * micron );
	  }
	else
	  {
	     ftree->trk_jetTrk_found.push_back( false );
	     
	     ftree->trk_jetTrk_deltaR.push_back( null );
	     
	     ftree->trk_jetTrk_pt.push_back( null );
	     ftree->trk_jetTrk_px.push_back( null );
	     ftree->trk_jetTrk_py.push_back( null );
	     ftree->trk_jetTrk_pz.push_back( null );
	     ftree->trk_jetTrk_p.push_back( null );
	     ftree->trk_jetTrk_eta.push_back( null );
	     ftree->trk_jetTrk_phi.push_back( null );
	     
	     ftree->trk_jetTrk_nTrackerLayers.push_back( null );
	     ftree->trk_jetTrk_nPixelBarrelLayers.push_back( null );
	     ftree->trk_jetTrk_nPixelEndcapLayers.push_back( null );
	     ftree->trk_jetTrk_nStripLayers.push_back( null );
	     
	     ftree->trk_jetTrk_nValid.push_back( null );
	     ftree->trk_jetTrk_fValid.push_back( null );
	     ftree->trk_jetTrk_nValidTracker.push_back( null );
	     ftree->trk_jetTrk_nValidPixelBarrel.push_back( null );
	     ftree->trk_jetTrk_nValidPixelEndcap.push_back( null );
	     ftree->trk_jetTrk_nValidStrip.push_back( null );
	     
	     ftree->trk_jetTrk_nMissed.push_back( null );
	     ftree->trk_jetTrk_nMissedOut.push_back( null );
	     ftree->trk_jetTrk_nMissedIn.push_back( null );
	     ftree->trk_jetTrk_nMissedTrackerOut.push_back( null );
	     ftree->trk_jetTrk_nMissedTrackerIn.push_back( null );
	     ftree->trk_jetTrk_nMissedPixelBarrelOut.push_back( null );
	     ftree->trk_jetTrk_nMissedPixelBarrelIn.push_back( null );
	     ftree->trk_jetTrk_nMissedPixelEndcapOut.push_back( null );
	     ftree->trk_jetTrk_nMissedPixelEndcapIn.push_back( null );
	     
	     ftree->trk_jetTrk_hasPixelBarrelLayer1.push_back( false );
	     ftree->trk_jetTrk_hasPixelEndcapLayer1.push_back( false );
	     ftree->trk_jetTrk_hasPixelBarrelLayer2.push_back( false );
	     ftree->trk_jetTrk_hasPixelEndcapLayer2.push_back( false );
	     ftree->trk_jetTrk_hasPixelBarrelLayer3.push_back( false );
	     ftree->trk_jetTrk_hasPixelEndcapLayer3.push_back( false );
	     ftree->trk_jetTrk_hasPixelBarrelLayer4.push_back( false );
	     ftree->trk_jetTrk_hasPixelEndcapLayer4.push_back( false );
	     
	     ftree->trk_jetTrk_quality.push_back( null );
	     ftree->trk_jetTrk_normalizedChi2.push_back( null );
	     ftree->trk_jetTrk_ndof.push_back( null );
	     ftree->trk_jetTrk_charge.push_back( null );
	     ftree->trk_jetTrk_qoverp.push_back( null );
	     ftree->trk_jetTrk_qoverpError.push_back( null );
	     ftree->trk_jetTrk_theta.push_back( null );
	     ftree->trk_jetTrk_thetaError.push_back( null );
	     ftree->trk_jetTrk_lambda.push_back( null );
	     ftree->trk_jetTrk_lambdaError.push_back( null );
	     ftree->trk_jetTrk_ptError.push_back( null );
	     ftree->trk_jetTrk_etaError.push_back( null );
	     ftree->trk_jetTrk_phiError.push_back( null );
	     
	     ftree->trk_jetTrk_d0.push_back( null );
	     ftree->trk_jetTrk_dz.push_back( null );
	     ftree->trk_jetTrk_d0_pv.push_back( null );
	     ftree->trk_jetTrk_dz_pv.push_back( null );
	     ftree->trk_jetTrk_d0_bs.push_back( null );
	     ftree->trk_jetTrk_d0_bs_zpca.push_back( null );
	     ftree->trk_jetTrk_d0_bs_zpv.push_back( null );
	     ftree->trk_jetTrk_dz_bs.push_back( null );
	     ftree->trk_jetTrk_d0Err.push_back( null );
	     ftree->trk_jetTrk_dzErr.push_back( null );
	     ftree->trk_jetTrk_d0_pv_NoRefit.push_back( null );
	     ftree->trk_jetTrk_dz_pv_NoRefit.push_back( null );
	  }

	// Tracks from PFJets

	float drMinPF = 10E+10;
	int iJetPF = -1;
	int iTrackMinPF = -1;
	
	for( unsigned int ij=0;ij<pfJets->size();ij++ )
	  {
	     reco::PFJet jet = pfJets->at(ij);
	     reco::TrackRefVector trks = jet.getTrackRefs();
	     const reco::TrackRef trkRef = reco::TrackRef(tracks, itk - tracks->begin());
	     edm::RefVector<TrackCollection>::const_iterator itt = find_if(trks.begin(), trks.end(), TrackEqualRef(trkRef));
	     if( itt != trks.end() )
	       {
		  size_t pos = itt - trks.begin();
		  
		  iJetPF = ij;
		  
		  for( unsigned int it=0;it<trks.size();it++ )
		    {
		       if( it == pos ) continue;
		       
		       float dr = getDeltaR(trks[pos]->eta(), trks[pos]->phi(), trks[it]->eta(), trks[it]->phi());
		       if( dr < drMin )
			 {
			    drMinPF = dr;
			    iTrackMinPF = it;
			 }
		    }
		  
		  break;
	       }
	  }

	if( iJetPF >= 0 )
	  {
	     reco::PFJet jet = pfJets->at(iJetPF);
	     
	     ftree->trk_pfjet_found.push_back( true );
	     
	     ftree->trk_pfjet_pt.push_back( jet.pt() );
	     ftree->trk_pfjet_eta.push_back( jet.eta() );
	     ftree->trk_pfjet_phi.push_back( jet.phi() );
	     ftree->trk_pfjet_nTracks.push_back( jet.getTrackRefs().size() );
	  }
	else
	  {
	     ftree->trk_pfjet_found.push_back( false );
	     
	     ftree->trk_pfjet_pt.push_back( null );
	     ftree->trk_pfjet_eta.push_back( null );
	     ftree->trk_pfjet_phi.push_back( null );
	     ftree->trk_pfjet_nTracks.push_back( null );
	  }	

	if( iTrackMinPF >= 0 )
	  {		  
	     reco::TrackRef trkCl = pfJets->at(iJetPF).getTrackRefs().at(iTrackMinPF);
	     
	     ftree->trk_pfjetTrk_found.push_back( true );
	     
	     ftree->trk_pfjetTrk_deltaR.push_back( drMinPF );
	     
	     ftree->trk_pfjetTrk_pt.push_back( trkCl->pt() );
	     ftree->trk_pfjetTrk_px.push_back( trkCl->px() );
	     ftree->trk_pfjetTrk_py.push_back( trkCl->py() );
	     ftree->trk_pfjetTrk_pz.push_back( trkCl->pz() );
	     ftree->trk_pfjetTrk_p.push_back( trkCl->p() );
	     ftree->trk_pfjetTrk_eta.push_back( trkCl->eta() );
	     ftree->trk_pfjetTrk_phi.push_back( trkCl->phi() );
	     
	     ftree->trk_pfjetTrk_nTrackerLayers.push_back( trkCl->hitPattern().trackerLayersWithMeasurement() );
	     ftree->trk_pfjetTrk_nPixelBarrelLayers.push_back( trkCl->hitPattern().pixelBarrelLayersWithMeasurement() );
	     ftree->trk_pfjetTrk_nPixelEndcapLayers.push_back( trkCl->hitPattern().pixelEndcapLayersWithMeasurement() );
	     ftree->trk_pfjetTrk_nStripLayers.push_back( trkCl->hitPattern().stripLayersWithMeasurement() );
	     
	     ftree->trk_pfjetTrk_nValid.push_back( trkCl->numberOfValidHits() );
	     ftree->trk_pfjetTrk_fValid.push_back( trkCl->validFraction() );
	     ftree->trk_pfjetTrk_nValidTracker.push_back( trkCl->hitPattern().numberOfValidTrackerHits() );
	     ftree->trk_pfjetTrk_nValidPixelBarrel.push_back( trkCl->hitPattern().numberOfValidPixelBarrelHits() );
	     ftree->trk_pfjetTrk_nValidPixelEndcap.push_back( trkCl->hitPattern().numberOfValidPixelEndcapHits() );
	     ftree->trk_pfjetTrk_nValidStrip.push_back( trkCl->hitPattern().numberOfValidStripHits() );
	     
	     ftree->trk_pfjetTrk_nMissed.push_back( trkCl->numberOfLostHits() );
	     ftree->trk_pfjetTrk_nMissedOut.push_back( trkCl->hitPattern().numberOfLostHits(HitPattern::MISSING_OUTER_HITS) );
	     ftree->trk_pfjetTrk_nMissedIn.push_back( trkCl->hitPattern().numberOfLostHits(HitPattern::MISSING_INNER_HITS) );
	     ftree->trk_pfjetTrk_nMissedTrackerOut.push_back( trkCl->hitPattern().numberOfLostTrackerHits(HitPattern::MISSING_OUTER_HITS) );
	     ftree->trk_pfjetTrk_nMissedTrackerIn.push_back( trkCl->hitPattern().numberOfLostTrackerHits(HitPattern::MISSING_INNER_HITS) );
	     ftree->trk_pfjetTrk_nMissedPixelBarrelOut.push_back( trkCl->hitPattern().numberOfLostPixelBarrelHits(HitPattern::MISSING_OUTER_HITS) );
	     ftree->trk_pfjetTrk_nMissedPixelBarrelIn.push_back( trkCl->hitPattern().numberOfLostPixelBarrelHits(HitPattern::MISSING_INNER_HITS) );
	     ftree->trk_pfjetTrk_nMissedPixelEndcapOut.push_back( trkCl->hitPattern().numberOfLostPixelEndcapHits(HitPattern::MISSING_OUTER_HITS) );
	     ftree->trk_pfjetTrk_nMissedPixelEndcapIn.push_back( trkCl->hitPattern().numberOfLostPixelEndcapHits(HitPattern::MISSING_INNER_HITS) );
	     
	     ftree->trk_pfjetTrk_hasPixelBarrelLayer1.push_back( trkCl->hitPattern().hasValidHitInPixelLayer(PixelSubdetector::SubDetector::PixelBarrel, 1) );
	     ftree->trk_pfjetTrk_hasPixelEndcapLayer1.push_back( trkCl->hitPattern().hasValidHitInPixelLayer(PixelSubdetector::SubDetector::PixelEndcap, 1) );
	     ftree->trk_pfjetTrk_hasPixelBarrelLayer2.push_back( trkCl->hitPattern().hasValidHitInPixelLayer(PixelSubdetector::SubDetector::PixelBarrel, 2) );
	     ftree->trk_pfjetTrk_hasPixelEndcapLayer2.push_back( trkCl->hitPattern().hasValidHitInPixelLayer(PixelSubdetector::SubDetector::PixelEndcap, 2) );
	     ftree->trk_pfjetTrk_hasPixelBarrelLayer3.push_back( trkCl->hitPattern().hasValidHitInPixelLayer(PixelSubdetector::SubDetector::PixelBarrel, 3) );
	     ftree->trk_pfjetTrk_hasPixelEndcapLayer3.push_back( trkCl->hitPattern().hasValidHitInPixelLayer(PixelSubdetector::SubDetector::PixelEndcap, 3) );
	     ftree->trk_pfjetTrk_hasPixelBarrelLayer4.push_back( trkCl->hitPattern().hasValidHitInPixelLayer(PixelSubdetector::SubDetector::PixelBarrel, 4) );
	     ftree->trk_pfjetTrk_hasPixelEndcapLayer4.push_back( trkCl->hitPattern().hasValidHitInPixelLayer(PixelSubdetector::SubDetector::PixelEndcap, 4) );
	     
	     ftree->trk_pfjetTrk_quality.push_back( trkCl->qualityMask() );
	     ftree->trk_pfjetTrk_normalizedChi2.push_back( trkCl->normalizedChi2() );
	     ftree->trk_pfjetTrk_ndof.push_back( trkCl->ndof() );
	     ftree->trk_pfjetTrk_charge.push_back( trkCl->charge() );
	     ftree->trk_pfjetTrk_qoverp.push_back( trkCl->qoverp() );
	     ftree->trk_pfjetTrk_qoverpError.push_back( trkCl->qoverpError() );
	     ftree->trk_pfjetTrk_theta.push_back( trkCl->theta() );
	     ftree->trk_pfjetTrk_thetaError.push_back( trkCl->thetaError() );
	     ftree->trk_pfjetTrk_lambda.push_back( trkCl->lambda() );
	     ftree->trk_pfjetTrk_lambdaError.push_back( trkCl->lambdaError() );
	     ftree->trk_pfjetTrk_ptError.push_back( trkCl->ptError() );
	     ftree->trk_pfjetTrk_etaError.push_back( trkCl->etaError() );
	     ftree->trk_pfjetTrk_phiError.push_back( trkCl->phiError() );
	     
	     ftree->trk_pfjetTrk_d0.push_back( trkCl->dxy() * micron );
	     ftree->trk_pfjetTrk_dz.push_back( trkCl->dz() * micron );
	     ftree->trk_pfjetTrk_d0_pv.push_back( trkCl->dxy(vtxPosition) * micron );
	     ftree->trk_pfjetTrk_dz_pv.push_back( trkCl->dz(vtxPosition) * micron );
	     ftree->trk_pfjetTrk_d0_bs.push_back( trkCl->dxy(pvbeamspot->position()) * micron );
	     ftree->trk_pfjetTrk_d0_bs_zpca.push_back( trkCl->dxy(*pvbeamspot) * micron );
	     ftree->trk_pfjetTrk_d0_bs_zpv.push_back( trkCl->dxy(pvbeamspot->position(vtx.z())) * micron );
	     ftree->trk_pfjetTrk_dz_bs.push_back( trkCl->dz(pvbeamspot->position()) * micron );
	     ftree->trk_pfjetTrk_d0Err.push_back( trkCl->d0Error() * micron );
	     ftree->trk_pfjetTrk_dzErr.push_back( trkCl->dzError() * micron );
	     ftree->trk_pfjetTrk_d0_pv_NoRefit.push_back( trkCl->dxy(vtxH->front().position()) * micron );
	     ftree->trk_pfjetTrk_dz_pv_NoRefit.push_back( trkCl->dz(vtxH->front().position()) * micron );
	  }
	else
	  {
	     ftree->trk_pfjetTrk_found.push_back( false );
	     
	     ftree->trk_pfjetTrk_deltaR.push_back( null );
	     
	     ftree->trk_pfjetTrk_pt.push_back( null );
	     ftree->trk_pfjetTrk_px.push_back( null );
	     ftree->trk_pfjetTrk_py.push_back( null );
	     ftree->trk_pfjetTrk_pz.push_back( null );
	     ftree->trk_pfjetTrk_p.push_back( null );
	     ftree->trk_pfjetTrk_eta.push_back( null );
	     ftree->trk_pfjetTrk_phi.push_back( null );
	     
	     ftree->trk_pfjetTrk_nTrackerLayers.push_back( null );
	     ftree->trk_pfjetTrk_nPixelBarrelLayers.push_back( null );
	     ftree->trk_pfjetTrk_nPixelEndcapLayers.push_back( null );
	     ftree->trk_pfjetTrk_nStripLayers.push_back( null );
	     
	     ftree->trk_pfjetTrk_nValid.push_back( null );
	     ftree->trk_pfjetTrk_fValid.push_back( null );
	     ftree->trk_pfjetTrk_nValidTracker.push_back( null );
	     ftree->trk_pfjetTrk_nValidPixelBarrel.push_back( null );
	     ftree->trk_pfjetTrk_nValidPixelEndcap.push_back( null );
	     ftree->trk_pfjetTrk_nValidStrip.push_back( null );
	     
	     ftree->trk_pfjetTrk_nMissed.push_back( null );
	     ftree->trk_pfjetTrk_nMissedOut.push_back( null );
	     ftree->trk_pfjetTrk_nMissedIn.push_back( null );
	     ftree->trk_pfjetTrk_nMissedTrackerOut.push_back( null );
	     ftree->trk_pfjetTrk_nMissedTrackerIn.push_back( null );
	     ftree->trk_pfjetTrk_nMissedPixelBarrelOut.push_back( null );
	     ftree->trk_pfjetTrk_nMissedPixelBarrelIn.push_back( null );
	     ftree->trk_pfjetTrk_nMissedPixelEndcapOut.push_back( null );
	     ftree->trk_pfjetTrk_nMissedPixelEndcapIn.push_back( null );
	     
	     ftree->trk_pfjetTrk_hasPixelBarrelLayer1.push_back( false );
	     ftree->trk_pfjetTrk_hasPixelEndcapLayer1.push_back( false );
	     ftree->trk_pfjetTrk_hasPixelBarrelLayer2.push_back( false );
	     ftree->trk_pfjetTrk_hasPixelEndcapLayer2.push_back( false );
	     ftree->trk_pfjetTrk_hasPixelBarrelLayer3.push_back( false );
	     ftree->trk_pfjetTrk_hasPixelEndcapLayer3.push_back( false );
	     ftree->trk_pfjetTrk_hasPixelBarrelLayer4.push_back( false );
	     ftree->trk_pfjetTrk_hasPixelEndcapLayer4.push_back( false );
	     
	     ftree->trk_pfjetTrk_quality.push_back( null );
	     ftree->trk_pfjetTrk_normalizedChi2.push_back( null );
	     ftree->trk_pfjetTrk_ndof.push_back( null );
	     ftree->trk_pfjetTrk_charge.push_back( null );
	     ftree->trk_pfjetTrk_qoverp.push_back( null );
	     ftree->trk_pfjetTrk_qoverpError.push_back( null );
	     ftree->trk_pfjetTrk_theta.push_back( null );
	     ftree->trk_pfjetTrk_thetaError.push_back( null );
	     ftree->trk_pfjetTrk_lambda.push_back( null );
	     ftree->trk_pfjetTrk_lambdaError.push_back( null );
	     ftree->trk_pfjetTrk_ptError.push_back( null );
	     ftree->trk_pfjetTrk_etaError.push_back( null );
	     ftree->trk_pfjetTrk_phiError.push_back( null );
	     
	     ftree->trk_pfjetTrk_d0.push_back( null );
	     ftree->trk_pfjetTrk_dz.push_back( null );
	     ftree->trk_pfjetTrk_d0_pv.push_back( null );
	     ftree->trk_pfjetTrk_dz_pv.push_back( null );
	     ftree->trk_pfjetTrk_d0_bs.push_back( null );
	     ftree->trk_pfjetTrk_d0_bs_zpca.push_back( null );
	     ftree->trk_pfjetTrk_d0_bs_zpv.push_back( null );
	     ftree->trk_pfjetTrk_dz_bs.push_back( null );
	     ftree->trk_pfjetTrk_d0Err.push_back( null );
	     ftree->trk_pfjetTrk_dzErr.push_back( null );
	     ftree->trk_pfjetTrk_d0_pv_NoRefit.push_back( null );
	     ftree->trk_pfjetTrk_dz_pv_NoRefit.push_back( null );
	  }	
     }

   ftree->tree->Fill();
}

void Residuals::beginRun(const edm::Run& iRun, const edm::EventSetup& iSetup)
{
   bool changed = true;
   if( !hltConfig_.init(iRun, iSetup, "HLT", changed) )
     std::cout << "Warning, didn't find HLTConfigProvider with label "
     << "HLT" << " in run " << iRun.run() << std::endl;

   if( !hltPrescale_.init(iRun, iSetup, "HLT", changed) )
     std::cout << "Warning, didn't find HLTPrescaleProvider with label "
     << "HLT" << " in run " << iRun.run() << std::endl;
}

void Residuals::endRun() 
{
}

bool Residuals::trackSelection(const reco::Track& track) const 
{   
   using namespace reco;
   
   if( track.pt() < tkMinPt ) return false;
//   if( track.hitPattern().trackerLayersWithMeasurement() < tkMinXLayers ) return false;
//   if( track.trackerExpectedHitsOuter().numberOfLostHits() > tkMaxMissedOuterLayers ) return false;
//   if( track.trackerExpectedHitsInner().numberOfLostHits() > tkMaxMissedInnerLayers ) return false;   
   if( ! track.quality(reco::TrackBase::highPurity) ) return false;
//   if( ! (track.hitPattern().hasValidHitInPixelLayer(PixelSubdetector::PixelBarrel, 1) ||
//	  track.hitPattern().hasValidHitInPixelLayer(PixelSubdetector::PixelEndcap, 1)) ) return false;
   
   return true;
}

bool Residuals::vertexSelection(const reco::Vertex& vertex) const
{
   if( vertex.tracksSize()>vtxTracksSizeMax || vertex.tracksSize()<vtxTracksSizeMin ) return false;
//   if( vertex.xError() < vtxErrorXMin || vertex.xError() > vtxErrorXMax ) return false;
//   if( vertex.yError() < vtxErrorYMin || vertex.yError() > vtxErrorYMax ) return false;
//   if( vertex.zError() < vtxErrorZMin || vertex.zError() > vtxErrorZMax ) return false;
   
   return true;
}

float Residuals::getDeltaR(float eta1, float phi1, float eta2, float phi2)
{      
   float DeltaPhi = fabs(phi2 - phi1);
   if (DeltaPhi > 3.141593 ) DeltaPhi -= 2.*3.141593;
   return sqrt( (eta2-eta1)*(eta2-eta1) + DeltaPhi*DeltaPhi );
}

//define this as a plug-in
DEFINE_FWK_MODULE(Residuals);
