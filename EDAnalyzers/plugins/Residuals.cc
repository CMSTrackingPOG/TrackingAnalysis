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

#include "SimTracker/Records/interface/TrackAssociatorRecord.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"
#include "SimTracker/TrackAssociation/interface/TrackingParticleIP.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingVertex.h"
#include "SimDataFormats/Associations/interface/TrackToTrackingParticleAssociator.h"
#include "SimTracker/TrackAssociation/plugins/ParametersDefinerForTPESProducer.h"
#include "SimTracker/TrackAssociation/interface/ParametersDefinerForTP.h"
#include "SimDataFormats/Associations/interface/VertexToTrackingVertexAssociator.h"
#include "SimGeneral/HepPDTRecord/interface/ParticleDataTable.h"
#include "SimTracker/TrackHistory/interface/VertexClassifier.h"
#include "SimTracker/TrackHistory/interface/TrackClassifier.h"

#include "DataFormats/JetReco/interface/TrackJet.h"
#include "DataFormats/JetReco/interface/PFJet.h"

#include <TrackingTools/TrajectoryState/interface/PerigeeConversions.h>
#include <TrackingTools/TrajectoryState/interface/TrajectoryStateClosestToPoint.h>
#include <TrackingTools/PatternTools/interface/TSCPBuilderNoMaterial.h>

#include "TrackingTools/IPTools/interface/IPTools.h"

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
	
	TrackEqual( const edm::Ptr<reco::Track> & t ) : track_( t ) {}
	
	bool operator()( const edm::Ptr<reco::Track> & t ) const {
	   return t->pt()==track_->pt();
	}
	
      private:
	
	const edm::Ptr<reco::Track> & track_;
     };

   class TrackEqualReco
     {	
      public:
	
	TrackEqualReco( const reco::Track & t ) : track_( t ) {}
	
	bool operator()( const reco::Track & t ) const {
	   return t.pt()==track_.pt();
	}
	
      private:
	
	const reco::Track & track_;
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
   edm::EDGetTokenT<reco::TrackCollection> theTracksToken_;
   edm::EDGetTokenT<edm::View<reco::Track> > theTrackViewsToken_;
   edm::EDGetTokenT<reco::BeamSpot> theBeamspotToken_;
   edm::EDGetTokenT<double> theRhoToken_;
   edm::EDGetTokenT< vector<reco::TrackJet> > theTrackJetsToken_;
   edm::EDGetTokenT< vector<reco::PFJet> > thePFJetsToken_;
   edm::EDGetTokenT<edm::TriggerResults> theTriggerBitsToken_;
   edm::EDGetTokenT<std::vector<PileupSummaryInfo> > puInfoToken_;
   edm::EDGetTokenT<reco::TrackToTrackingParticleAssociator> theTrackAssociatorToken_;
   edm::EDGetTokenT<TrackingParticleCollection> theTrackingParticleToken_;
   edm::EDGetTokenT<TrackingVertexCollection> theTrackingVertexToken_;
   edm::EDGetTokenT<reco::VertexToTrackingVertexAssociator> theVertexAssociatorToken_;

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
   bool doTruth;

   int eventScale;
   int trackScale;
   
   int PVFitNTracksMin;
   int PVFitNMax;
   
   TRandom3 *rnd;
   
   HLTConfigProvider hltConfig_;
   HLTPrescaleProvider hltPrescale_;
   
   VertexClassifier vtxClassifier_;
   TrackClassifier trkClassifier_;
   
   const edm::Service<TFileService> fs;
   ResTree* ftree;
   
   int ncount;
};

Residuals::Residuals(const edm::ParameterSet& pset):
   hltPrescale_(pset, consumesCollector(), *this),
   vtxClassifier_(pset, consumesCollector()),
   trkClassifier_(pset, consumesCollector())
{
   edm::InputTag TrackCollectionTag_ = pset.getParameter<edm::InputTag>("TrackLabel");
   theTracksToken_= consumes<reco::TrackCollection>(TrackCollectionTag_);

   theTrackViewsToken_= consumes<edm::View<reco::Track> >(TrackCollectionTag_);
   
   edm::InputTag VertexCollectionTag_ = pset.getParameter<edm::InputTag>("VertexLabel");
   thePVToken_ = consumes<reco::VertexCollection>(VertexCollectionTag_);
   
   edm::InputTag BeamspotTag_ = edm::InputTag("offlineBeamSpot");
   theBeamspotToken_ = consumes<reco::BeamSpot>(BeamspotTag_);

   edm::InputTag RhoTag_ = pset.getParameter<edm::InputTag>("RhoLabel");
   theRhoToken_ = consumes<double>(RhoTag_);

   edm::InputTag TrackJetsTag_ = pset.getParameter<edm::InputTag>("TrackJetsLabel");
   theTrackJetsToken_ = consumes< vector<reco::TrackJet> >(TrackJetsTag_);

   edm::InputTag PFJetsTag_ = pset.getParameter<edm::InputTag>("PFJetsLabel");
   thePFJetsToken_ = consumes< vector<reco::PFJet> >(PFJetsTag_);
   
   edm::InputTag TriggerBitsTag_ = pset.getParameter<edm::InputTag>("TriggerResultsLabel");
   theTriggerBitsToken_ = consumes<edm::TriggerResults>(TriggerBitsTag_);
   
   edm::InputTag PUInfoTag_ = pset.getParameter<edm::InputTag>("puInfoLabel");
   puInfoToken_ = consumes<std::vector<PileupSummaryInfo> >(PUInfoTag_);

   edm::InputTag TrackingParticleTag_ = pset.getParameter<edm::InputTag>("TrackingParticleLabel");
   theTrackingParticleToken_ = consumes<TrackingParticleCollection>(TrackingParticleTag_);

   edm::InputTag TrackingVertexTag_ = pset.getParameter<edm::InputTag>("TrackingVertexLabel");
   theTrackingVertexToken_ = consumes<TrackingVertexCollection>(TrackingVertexTag_);
   
   edm::InputTag TrackAssociatorTag_ = pset.getParameter<edm::InputTag>("TrackAssociatorLabel");
   theTrackAssociatorToken_ = consumes<reco::TrackToTrackingParticleAssociator>(TrackAssociatorTag_);
   
   edm::InputTag VertexAssociatorTag_ = pset.getParameter<edm::InputTag>("VertexAssociatorLabel");
   theVertexAssociatorToken_ = consumes<reco::VertexToTrackingVertexAssociator>(VertexAssociatorTag_);
   
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
   
   PVFitNTracksMin = pset.getParameter<int>("PVFitNTracksMin");
   PVFitNMax = pset.getParameter<int>("PVFitNMax");
   
   runOnData = pset.getParameter<bool>("RunOnData");
   doTruth = pset.getParameter<bool>("DoTruth");
   
   rnd = new TRandom3();

   TFile& f = fs->file();
   f.SetCompressionAlgorithm(ROOT::kZLIB);
   f.SetCompressionLevel(9);
   ftree = new ResTree(fs->make<TTree>("tree", "tree"));
   ftree->CreateBranches(32000, runOnData);
   
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
   
   Handle<View<Track> >  trackViews;
   iEvent.getByToken(theTrackViewsToken_, trackViews);
   
   Handle<VertexCollection> vtxH;
   iEvent.getByToken(thePVToken_, vtxH);

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
   
   // refit primary vertices and put it in a new handle (note: TrackBaseRefs are different)
   vector<TransientVertex> pvs = revertex.makeVertices(tracksAll, *pvbeamspot, iSetup);
   if( pvs.empty() ) return;
   
   reco::VertexCollection pvr;
   for(unsigned int ipv=0;ipv<pvs.size();ipv++) 
     pvr.push_back(reco::Vertex(pvs[ipv]));
   
   edm::Handle<reco::VertexCollection> pvrh = edm::Handle(const_cast<reco::VertexCollection*>(&pvr), vtxH.provenance());
  
   reco::Vertex vtx = pvr.front();
   TransientVertex vtxTrans = pvs.front();
   
   if( !vertexSelection(vtx) ) return;

   ESHandle<MagneticField> theMF;
   iSetup.get<IdealMagneticFieldRecord>().get(theMF);

   Handle<double> rhoPtr;
   iEvent.getByToken(theRhoToken_, rhoPtr);

   // General event info
   ftree->ev_run = iEvent.id().run();
   ftree->ev_id = iEvent.id().event();
   ftree->ev_lumi = iEvent.id().luminosityBlock();
   ftree->ev_bunchCrossing = iEvent.bunchCrossing();
   ftree->ev_orbitNumber = iEvent.orbitNumber();
   ftree->ev_rho = *rhoPtr;
   ftree->ev_nPV =  pvr.size();

   // Trigger
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

   // Pileup info
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
   
   Handle<TrackingVertexCollection> trackingVertex;
   ESHandle<ParametersDefinerForTP> parametersDefinerTP;
   
   reco::RecoToSimCollection recSimCollTracks;

   if( doTruth && !runOnData )
     {
	vtxClassifier_.newEvent(iEvent, iSetup);
	trkClassifier_.newEvent(iEvent, iSetup);
	
	Handle<TrackingParticleCollection> trackingParticle;
	Handle<reco::TrackToTrackingParticleAssociator> trackAssociator;
	Handle<reco::VertexToTrackingVertexAssociator> vertexAssociator;

	reco::VertexRecoToSimCollection recSimCollVtx;
	
	iEvent.getByToken(theTrackingParticleToken_, trackingParticle);
	iEvent.getByToken(theTrackingVertexToken_, trackingVertex);
	iEvent.getByToken(theTrackAssociatorToken_, trackAssociator);
	iEvent.getByToken(theVertexAssociatorToken_, vertexAssociator);
	
	if( !trackingParticle.isValid() ) {
	   std::cout << "Can not access Tracking Particles" << std::endl; 
	   exit(1);
	}	
	if( !trackingVertex.isValid() ) {
	   std::cout << "Can not access Tracking Vertices" << std::endl;
	   exit(1);
	}	
	if( !trackAssociator.isValid() ) {
	   std::cout << "Can not access Tracking Associator" << std::endl;
	   exit(1);
	}
	if( !vertexAssociator.isValid() ) {
	   std::cout << "Can not access Vertex Associator" << std::endl;
	   exit(1);	   
	}
	
	std::cout << "TrackingParticles=" << trackingParticle->size() << std::endl;
	std::cout << "TrackingVertices=" << trackingVertex->size() << std::endl;
	std::cout << "Tracks=" << trackViews->size() << std::endl;
	std::cout << "Vertices=" << pvr.size() << std::endl;
	
	iSetup.get<TrackAssociatorRecord>().get("LhcParametersDefinerForTP", parametersDefinerTP);
	
	// match reco-sim tracks
	recSimCollTracks = trackAssociator->associateRecoToSim(trackViews, trackingParticle);

	// match reco-sim vertices
	recSimCollVtx = vertexAssociator->associateRecoToSim(vtxH, trackingVertex, recSimCollTracks);

	for( unsigned int ipv=0;ipv<pvr.size();ipv++ )
	  {	     	
	     reco::VertexRef recoVertex(pvrh, ipv);
	     auto matched = recSimCollVtx.find(recoVertex);

	     bool pv_mc_hasMatch = 0;
	     std::vector<float> pv_mc_matchQuality;
	     std::vector<bool> pv_mc_isFake;
	     std::vector<bool> pv_mc_isPrimaryVertex;
	     std::vector<bool> pv_mc_isSecondaryVertex;
	     std::vector<bool> pv_mc_isTertiaryVertex;
	     std::vector<bool> pv_mc_isSignalEvent;
	     std::vector<bool> pv_mc_isBWeakDecay;
	     std::vector<bool> pv_mc_isCWeakDecay;
	     std::vector<bool> pv_mc_isTauDecay;
	     std::vector<bool> pv_mc_isKsDecay;
	     std::vector<bool> pv_mc_isLambdaDecay;
	     std::vector<bool> pv_mc_isJpsiDecay;
	     std::vector<bool> pv_mc_isXiDecay;
	     std::vector<bool> pv_mc_isOmegaDecay;
	     std::vector<bool> pv_mc_isSigmaPlusDecay;
	     std::vector<bool> pv_mc_isSigmaMinusDecay;
	     std::vector<bool> pv_mc_isLongLivedDecay;
	     
	     std::vector<bool> pv_mc_isKnownProcess;
	     std::vector<bool> pv_mc_isUndefinedProcess;
	     std::vector<bool> pv_mc_isUnknownProcess;
	     std::vector<bool> pv_mc_isPrimaryProcess;
	     std::vector<bool> pv_mc_isHadronicProcess;
	     std::vector<bool> pv_mc_isDecayProcess;
	     std::vector<bool> pv_mc_isComptonProcess;
	     std::vector<bool> pv_mc_isAnnihilationProcess;
	     std::vector<bool> pv_mc_isEIoniProcess;
	     std::vector<bool> pv_mc_isHIoniProcess;
	     std::vector<bool> pv_mc_isMuIoniProcess;
	     std::vector<bool> pv_mc_isPhotonProcess;
	     std::vector<bool> pv_mc_isMuPairProdProcess;
	     std::vector<bool> pv_mc_isConversionsProcess;
	     std::vector<bool> pv_mc_isEBremProcess;
	     std::vector<bool> pv_mc_isSynchrotronRadiationProcess;
	     std::vector<bool> pv_mc_isMuBremProcess;
	     std::vector<bool> pv_mc_isMuNuclProcess;
	     std::vector<bool> pv_mc_isUnknown;
	     
	     std::vector<bool> pv_mc_inVolume;
	     std::vector<float> pv_mc_x;
	     std::vector<float> pv_mc_y;
	     std::vector<float> pv_mc_z;
	     std::vector<float> pv_mc_t;
	     std::vector<int> pv_mc_nGenVtx;
	     std::vector<int> pv_mc_nSimVtx;
	     std::vector<int> pv_mc_nDaughterTracks;
	     std::vector<int> pv_mc_nSourceTracks;
	
	     if( matched != recSimCollVtx.end() )
	       {
		  pv_mc_hasMatch = 1;

		  for(const auto vertexRefQuality: matched->val) 
		    {
		       const TrackingVertexRef* tvPtr = &(vertexRefQuality.first);
		       const TrackingVertexRef& tv = *tvPtr;
		       pv_mc_matchQuality.push_back( vertexRefQuality.second );
		       
		       vtxClassifier_.evaluate(*tvPtr);
		       
		       pv_mc_isFake.push_back( vtxClassifier_.is(VertexCategories::Fake) ); // no match to any simulated vertex
		       pv_mc_isPrimaryVertex.push_back( vtxClassifier_.is(VertexCategories::PrimaryVertex) ); // no other vertex within vertexClusteringDistance
		       pv_mc_isSecondaryVertex.push_back( vtxClassifier_.is(VertexCategories::SecondaryVertex) ); // one vertex found within vertexClusteringDistance
		       pv_mc_isTertiaryVertex.push_back( vtxClassifier_.is(VertexCategories::TertiaryVertex) ); // two or more vertices found within vertexClusteringDistance
		       pv_mc_isSignalEvent.push_back( vtxClassifier_.is(VertexCategories::SignalEvent) ); // produced by the signal part of the crossing frame
		       pv_mc_isBWeakDecay.push_back( vtxClassifier_.is(VertexCategories::BWeakDecay) );
		       pv_mc_isCWeakDecay.push_back( vtxClassifier_.is(VertexCategories::CWeakDecay) );
		       pv_mc_isTauDecay.push_back( vtxClassifier_.is(VertexCategories::TauDecay) );
		       pv_mc_isKsDecay.push_back( vtxClassifier_.is(VertexCategories::KsDecay) );
		       pv_mc_isLambdaDecay.push_back( vtxClassifier_.is(VertexCategories::LambdaDecay) );
		       pv_mc_isJpsiDecay.push_back( vtxClassifier_.is(VertexCategories::JpsiDecay) );
		       pv_mc_isXiDecay.push_back( vtxClassifier_.is(VertexCategories::XiDecay) );
		       pv_mc_isOmegaDecay.push_back( vtxClassifier_.is(VertexCategories::OmegaDecay) );
		       pv_mc_isSigmaPlusDecay.push_back( vtxClassifier_.is(VertexCategories::SigmaPlusDecay) );
		       pv_mc_isSigmaMinusDecay.push_back( vtxClassifier_.is(VertexCategories::SigmaMinusDecay) );
		       pv_mc_isLongLivedDecay.push_back( vtxClassifier_.is(VertexCategories::LongLivedDecay) );
		       
		       pv_mc_isKnownProcess.push_back( vtxClassifier_.is(VertexCategories::KnownProcess) );
		       pv_mc_isUndefinedProcess.push_back( vtxClassifier_.is(VertexCategories::UndefinedProcess) );
		       pv_mc_isUnknownProcess.push_back( vtxClassifier_.is(VertexCategories::UnknownProcess) );
		       pv_mc_isPrimaryProcess.push_back( vtxClassifier_.is(VertexCategories::PrimaryProcess) );
		       pv_mc_isHadronicProcess.push_back( vtxClassifier_.is(VertexCategories::HadronicProcess) );
		       pv_mc_isDecayProcess.push_back( vtxClassifier_.is(VertexCategories::DecayProcess) );
		       pv_mc_isComptonProcess.push_back( vtxClassifier_.is(VertexCategories::ComptonProcess) );
		       pv_mc_isAnnihilationProcess.push_back( vtxClassifier_.is(VertexCategories::AnnihilationProcess) );
		       pv_mc_isEIoniProcess.push_back( vtxClassifier_.is(VertexCategories::EIoniProcess) );
		       pv_mc_isHIoniProcess.push_back( vtxClassifier_.is(VertexCategories::HIoniProcess) );
		       pv_mc_isMuIoniProcess.push_back( vtxClassifier_.is(VertexCategories::MuIoniProcess) );
		       pv_mc_isPhotonProcess.push_back( vtxClassifier_.is(VertexCategories::PhotonProcess) );
		       pv_mc_isMuPairProdProcess.push_back( vtxClassifier_.is(VertexCategories::MuPairProdProcess) );
		       pv_mc_isConversionsProcess.push_back( vtxClassifier_.is(VertexCategories::ConversionsProcess) );
		       pv_mc_isEBremProcess.push_back( vtxClassifier_.is(VertexCategories::EBremProcess) );
		       pv_mc_isSynchrotronRadiationProcess.push_back( vtxClassifier_.is(VertexCategories::SynchrotronRadiationProcess) );
		       pv_mc_isMuBremProcess.push_back( vtxClassifier_.is(VertexCategories::MuBremProcess) );
		       pv_mc_isMuNuclProcess.push_back( vtxClassifier_.is(VertexCategories::MuNuclProcess) );
		       pv_mc_isUnknown.push_back( vtxClassifier_.is(VertexCategories::Unknown) );
		       
		       pv_mc_inVolume.push_back( tv->inVolume() );
		       pv_mc_x.push_back( tv->position().x() * micron );
		       pv_mc_y.push_back( tv->position().y() * micron );
		       pv_mc_z.push_back( tv->position().z() * micron );
		       pv_mc_t.push_back( tv->position().t() );
		       pv_mc_nGenVtx.push_back( tv->nGenVertices() );
		       pv_mc_nSimVtx.push_back( tv->nG4Vertices() );
		       pv_mc_nDaughterTracks.push_back( tv->nDaughterTracks() );
		       pv_mc_nSourceTracks.push_back( tv->nSourceTracks() );
		    }
	       }
	     
	     ftree->pv_mc_hasMatch.push_back( pv_mc_hasMatch );
	     ftree->pv_mc_matchQuality.push_back( pv_mc_matchQuality );
	     ftree->pv_mc_isFake.push_back( pv_mc_isFake );
	     ftree->pv_mc_isPrimaryVertex.push_back( pv_mc_isPrimaryVertex );
	     ftree->pv_mc_isSecondaryVertex.push_back( pv_mc_isSecondaryVertex );
	     ftree->pv_mc_isTertiaryVertex.push_back( pv_mc_isTertiaryVertex );
	     ftree->pv_mc_isSignalEvent.push_back( pv_mc_isSignalEvent );
	     ftree->pv_mc_isBWeakDecay.push_back( pv_mc_isBWeakDecay );
	     ftree->pv_mc_isCWeakDecay.push_back( pv_mc_isCWeakDecay );
	     ftree->pv_mc_isTauDecay.push_back( pv_mc_isTauDecay );
	     ftree->pv_mc_isKsDecay.push_back( pv_mc_isKsDecay );
	     ftree->pv_mc_isLambdaDecay.push_back( pv_mc_isLambdaDecay );
	     ftree->pv_mc_isJpsiDecay.push_back( pv_mc_isJpsiDecay );
	     ftree->pv_mc_isXiDecay.push_back( pv_mc_isXiDecay );
	     ftree->pv_mc_isOmegaDecay.push_back( pv_mc_isOmegaDecay );
	     ftree->pv_mc_isSigmaPlusDecay.push_back( pv_mc_isSigmaPlusDecay );
	     ftree->pv_mc_isSigmaMinusDecay.push_back( pv_mc_isSigmaMinusDecay );
	     ftree->pv_mc_isLongLivedDecay.push_back( pv_mc_isLongLivedDecay );
	     
	     ftree->pv_mc_isKnownProcess.push_back( pv_mc_isKnownProcess );
	     ftree->pv_mc_isUndefinedProcess.push_back( pv_mc_isUndefinedProcess );
	     ftree->pv_mc_isUnknownProcess.push_back( pv_mc_isUnknownProcess );
	     ftree->pv_mc_isPrimaryProcess.push_back( pv_mc_isPrimaryProcess );
	     ftree->pv_mc_isHadronicProcess.push_back( pv_mc_isHadronicProcess );
	     ftree->pv_mc_isDecayProcess.push_back( pv_mc_isDecayProcess );
	     ftree->pv_mc_isComptonProcess.push_back( pv_mc_isComptonProcess );
	     ftree->pv_mc_isAnnihilationProcess.push_back( pv_mc_isAnnihilationProcess );
	     ftree->pv_mc_isEIoniProcess.push_back( pv_mc_isEIoniProcess );
	     ftree->pv_mc_isHIoniProcess.push_back( pv_mc_isHIoniProcess );
	     ftree->pv_mc_isMuIoniProcess.push_back( pv_mc_isMuIoniProcess );
	     ftree->pv_mc_isPhotonProcess.push_back( pv_mc_isPhotonProcess );
	     ftree->pv_mc_isMuPairProdProcess.push_back( pv_mc_isMuPairProdProcess );
	     ftree->pv_mc_isConversionsProcess.push_back( pv_mc_isConversionsProcess );
	     ftree->pv_mc_isEBremProcess.push_back( pv_mc_isEBremProcess );
	     ftree->pv_mc_isSynchrotronRadiationProcess.push_back( pv_mc_isSynchrotronRadiationProcess );
	     ftree->pv_mc_isMuBremProcess.push_back( pv_mc_isMuBremProcess );
	     ftree->pv_mc_isMuNuclProcess.push_back( pv_mc_isMuNuclProcess );
	     ftree->pv_mc_isUnknown.push_back( pv_mc_isUnknown );

	     ftree->pv_mc_inVolume.push_back( pv_mc_inVolume );
	     ftree->pv_mc_x.push_back( pv_mc_x );
	     ftree->pv_mc_y.push_back( pv_mc_y );
	     ftree->pv_mc_z.push_back( pv_mc_z );
	     ftree->pv_mc_t.push_back( pv_mc_t );
	     ftree->pv_mc_nGenVtx.push_back( pv_mc_nGenVtx );
	     ftree->pv_mc_nSimVtx.push_back( pv_mc_nSimVtx );
	     ftree->pv_mc_nDaughterTracks.push_back( pv_mc_nDaughterTracks );
	     ftree->pv_mc_nSourceTracks.push_back( pv_mc_nSourceTracks );
	  }
     }   

   // Primary vertex   
   for( unsigned int ipv=0;ipv<pvr.size();ipv++ )
     {	     	   
	std::vector<reco::TransientTrack> vtxTracks = pvs[ipv].originalTracks();
	stable_sort(vtxTracks.begin(), vtxTracks.end(), sortPt);
	
	int nTracks = pvr[ipv].tracksSize();
	
	Track::Point vtxPosition = Track::Point(pvr[ipv].position().x(),
						pvr[ipv].position().y(),
						pvr[ipv].position().z());	

	float pv_SumTrackPt = 0.;
	float pv_SumTrackPt2 = 0.;
	float pv_fracHighPurity = 0.;
	
	std::vector<float> pv_trackWeight;
	std::vector<bool> pv_trk_isHighPurity;
	std::vector<int> pv_trk_algo;
	std::vector<int> pv_trk_originalAlgo;
		
	std::vector<int> pv_trk_idx;
	
	std::vector<int> pv_trk_pvN;
	std::vector<int> pv_trk_pv1N;
	std::vector<int> pv_trk_pv2N;

	std::vector<bool> pv_trk_pvunbiased_IsValid;
	std::vector<bool> pv_trk_pvunbiased_IsFake;
	std::vector<int> pv_trk_pvunbiased_NTracks;
	std::vector<float> pv_trk_pvunbiased_SumTrackPt;
	std::vector<float> pv_trk_pvunbiased_SumTrackPt2;
	std::vector<float> pv_trk_pvunbiased_fracHighPurity;
	std::vector<float> pv_trk_pvunbiased_chi2;
	std::vector<int> pv_trk_pvunbiased_ndof;
	std::vector<float> pv_trk_pvunbiased_x;
	std::vector<float> pv_trk_pvunbiased_y;
	std::vector<float> pv_trk_pvunbiased_z;
	std::vector<float> pv_trk_pvunbiased_xError;
	std::vector<float> pv_trk_pvunbiased_yError;
	std::vector<float> pv_trk_pvunbiased_zError;
		  
	std::vector<float> pv_trk_d0_pvunbiased;
	std::vector<float> pv_trk_dz_pvunbiased;
	std::vector<float> pv_trk_d0_bs_zpvunbiased;

	std::vector<bool> pv_trk_pvunbiased_IsValid_p1;
	std::vector<bool> pv_trk_pvunbiased_IsFake_p1;
	std::vector<int> pv_trk_pvunbiased_NTracks_p1;
	std::vector<float> pv_trk_pvunbiased_SumTrackPt_p1;
	std::vector<float> pv_trk_pvunbiased_SumTrackPt2_p1;
	std::vector<float> pv_trk_pvunbiased_fracHighPurity_p1;
	std::vector<float> pv_trk_pvunbiased_chi2_p1;
	std::vector<int> pv_trk_pvunbiased_ndof_p1;
	std::vector<float> pv_trk_pvunbiased_x_p1;
	std::vector<float> pv_trk_pvunbiased_y_p1;
	std::vector<float> pv_trk_pvunbiased_z_p1;
	std::vector<float> pv_trk_pvunbiased_xError_p1;
	std::vector<float> pv_trk_pvunbiased_yError_p1;
	std::vector<float> pv_trk_pvunbiased_zError_p1;
		  
	std::vector<float> pv_trk_d0_pvunbiased_p1;
	std::vector<float> pv_trk_dz_pvunbiased_p1;
	std::vector<float> pv_trk_d0_bs_zpvunbiased_p1;
	
	std::vector<bool> pv_trk_pvunbiased_IsValid_p2;
	std::vector<bool> pv_trk_pvunbiased_IsFake_p2;
	std::vector<int> pv_trk_pvunbiased_NTracks_p2;
	std::vector<float> pv_trk_pvunbiased_SumTrackPt_p2;
	std::vector<float> pv_trk_pvunbiased_SumTrackPt2_p2;
	std::vector<float> pv_trk_pvunbiased_fracHighPurity_p2;
	std::vector<float> pv_trk_pvunbiased_chi2_p2;
	std::vector<int> pv_trk_pvunbiased_ndof_p2;
	std::vector<float> pv_trk_pvunbiased_x_p2;
	std::vector<float> pv_trk_pvunbiased_y_p2;
	std::vector<float> pv_trk_pvunbiased_z_p2;
	std::vector<float> pv_trk_pvunbiased_xError_p2;
	std::vector<float> pv_trk_pvunbiased_yError_p2;
	std::vector<float> pv_trk_pvunbiased_zError_p2;
		  
	std::vector<float> pv_trk_d0_pvunbiased_p2;
	std::vector<float> pv_trk_dz_pvunbiased_p2;
	std::vector<float> pv_trk_d0_bs_zpvunbiased_p2;
	
	std::vector<float> pv_trk_mc_dxy_pvunbiased;
	std::vector<float> pv_trk_mc_dz_pvunbiased;	
	std::vector<float> pv_trk_mc_dxy_tp_pvunbiased;
	std::vector<float> pv_trk_mc_dz_tp_pvunbiased;

	std::vector<float> pv_trk_mc_dxy_pvunbiased_p1;
	std::vector<float> pv_trk_mc_dz_pvunbiased_p1;	
	std::vector<float> pv_trk_mc_dxy_tp_pvunbiased_p1;
	std::vector<float> pv_trk_mc_dz_tp_pvunbiased_p1;

	std::vector<float> pv_trk_mc_dxy_pvunbiased_p2;
	std::vector<float> pv_trk_mc_dz_pvunbiased_p2;	
	std::vector<float> pv_trk_mc_dxy_tp_pvunbiased_p2;
	std::vector<float> pv_trk_mc_dz_tp_pvunbiased_p2;
	
	std::vector<float> pv_trk_pt;
	std::vector<float> pv_trk_px;
	std::vector<float> pv_trk_py;
	std::vector<float> pv_trk_pz;
	std::vector<float> pv_trk_p;
	std::vector<float> pv_trk_eta;
	std::vector<float> pv_trk_phi;

	std::vector<int> pv_trk_nTrackerLayers;
	std::vector<int> pv_trk_nPixelBarrelLayers;
	std::vector<int> pv_trk_nPixelEndcapLayers;
	std::vector<int> pv_trk_nStripLayers;
	     
	std::vector<int> pv_trk_nValid;
	std::vector<float> pv_trk_fValid;
	std::vector<int> pv_trk_nValidTracker;
	std::vector<int> pv_trk_nValidPixelBarrel;
	std::vector<int> pv_trk_nValidPixelEndcap;
	std::vector<int> pv_trk_nValidStrip;

	std::vector<int> pv_trk_nMissed;
	std::vector<int> pv_trk_nMissedOut;
	std::vector<int> pv_trk_nMissedIn;
	std::vector<int> pv_trk_nMissedTrackerOut;
	std::vector<int> pv_trk_nMissedTrackerIn;
	std::vector<int> pv_trk_nMissedPixelBarrelOut;
	std::vector<int> pv_trk_nMissedPixelBarrelIn;
	std::vector<int> pv_trk_nMissedPixelEndcapOut;
	std::vector<int> pv_trk_nMissedPixelEndcapIn;
	     
	std::vector<bool> pv_trk_hasPixelBarrelLayer1;
	std::vector<bool> pv_trk_hasPixelEndcapLayer1;
	std::vector<bool> pv_trk_hasPixelBarrelLayer2;
	std::vector<bool> pv_trk_hasPixelEndcapLayer2;
	std::vector<bool> pv_trk_hasPixelBarrelLayer3;
	std::vector<bool> pv_trk_hasPixelEndcapLayer3;
	std::vector<bool> pv_trk_hasPixelBarrelLayer4;
	std::vector<bool> pv_trk_hasPixelEndcapLayer4;
	     
	std::vector<int> pv_trk_quality;
	std::vector<float> pv_trk_normalizedChi2;
	std::vector<int> pv_trk_ndof;
	std::vector<int> pv_trk_charge;
	std::vector<float> pv_trk_qoverp;
	std::vector<float> pv_trk_qoverpError;
	std::vector<float> pv_trk_theta;
	std::vector<float> pv_trk_thetaError;
	std::vector<float> pv_trk_lambda;
	std::vector<float> pv_trk_lambdaError;
	std::vector<float> pv_trk_ptError;
	std::vector<float> pv_trk_etaError;
	std::vector<float> pv_trk_phiError;
	
	std::vector<float> pv_trk_d0;
	std::vector<float> pv_trk_dz;
	std::vector<float> pv_trk_d0_pv;
	std::vector<float> pv_trk_dz_pv;
	std::vector<float> pv_trk_d0_bs;
	std::vector<float> pv_trk_d0_bs_zpca;
	std::vector<float> pv_trk_d0_bs_zpv;
	std::vector<float> pv_trk_dz_bs;
	std::vector<float> pv_trk_d0Err;
	std::vector<float> pv_trk_dzErr;
	
	std::vector<float> pv_trk_d0_tv;
	std::vector<float> pv_trk_dz_tv;
	
	TrackCollection initPVTkCollection;
	for( std::vector<reco::TransientTrack>::const_iterator it = vtxTracks.begin(); it != vtxTracks.end(); it++ )
	  {
	     reco::Track trk = (*it).track();
	     initPVTkCollection.push_back(trk);
	  }	
	
	int iTrk = 0;
	for( std::vector<reco::TransientTrack>::const_iterator it = vtxTracks.begin(); it != vtxTracks.end(); it++ )
	  {	     
	     reco::Track trk = (*it).track();

	     // Since TrackBaseRefs are not preserved after vertex refitting, do a pt-based track matching
	     TrackCollection::const_iterator itt = find_if(tracks->begin(), tracks->end(), TrackEqualReco(trk));

	     if( itt != tracks->end() ) pv_trk_idx.push_back( itt - tracks->begin() );
	     else pv_trk_idx.push_back( -1 );
	     
	     // Remove the track from the PV track collection
	     TrackCollection newPVTkCollection;
	     newPVTkCollection.assign(initPVTkCollection.begin(), initPVTkCollection.begin()+iTrk);
	     newPVTkCollection.insert(newPVTkCollection.end(), initPVTkCollection.begin()+iTrk+1, initPVTkCollection.end());

	     // Split vertex tracks in two groups (unbiased)
	     reco::TrackCollection vtxTkCollection1;
	     reco::TrackCollection vtxTkCollection2;
	     
	     float SumTrackPt_p1 = 0;
	     float SumTrackPt2_p1 = 0;
	     
	     float SumTrackPt_p2 = 0;
	     float SumTrackPt2_p2 = 0;
	     
	     float pv_fracHighPurity_p1 = 0;
	     float pv_fracHighPurity_p2 = 0;
	     
	     for( std::vector<reco::Track>::const_iterator itv = newPVTkCollection.begin(); itv != newPVTkCollection.end(); itv++ )
	       {
		  reco::Track pvtrk = (*itv);
		  
		  if( rnd->Rndm() > 0.5 )
		    {
		       vtxTkCollection1.push_back(pvtrk);
		       SumTrackPt_p1 += pvtrk.pt();
		       SumTrackPt2_p1 += pvtrk.pt()*pvtrk.pt();
		       pv_fracHighPurity_p1 += pvtrk.quality(reco::TrackBase::highPurity);
		    }	
		  else
		    {	     
		       vtxTkCollection2.push_back(pvtrk);
		       SumTrackPt_p2 += pvtrk.pt();
		       SumTrackPt2_p2 += pvtrk.pt()*pvtrk.pt();
		       pv_fracHighPurity_p2 += pvtrk.quality(reco::TrackBase::highPurity);
		    }
	       }	     
	     
	     // Remake vertices
	     vector<TransientVertex> pvst = revertex.makeVertices(newPVTkCollection, *pvbeamspot, iSetup);
	     
	     vector<TransientVertex> pvst1 = revertex.makeVertices(vtxTkCollection1, *pvbeamspot, iSetup);
	     vector<TransientVertex> pvst2 = revertex.makeVertices(vtxTkCollection2, *pvbeamspot, iSetup);
	     
	     pv_trk_pvN.push_back( pvst.size() );
	     pv_trk_pv1N.push_back( pvst1.size() );
	     pv_trk_pv2N.push_back( pvst2.size() );
	     
	     if( !pvst1.empty() && !pvst2.empty() )
	       {
		  reco::Vertex vtx1 = reco::Vertex(pvst1.front());
		  reco::Vertex vtx2 = reco::Vertex(pvst2.front());
		  
		  pv_trk_pvunbiased_IsValid_p1.push_back( vtx1.isValid() );
		  pv_trk_pvunbiased_IsValid_p2.push_back( vtx2.isValid() );
		  
		  pv_trk_pvunbiased_IsFake_p1.push_back( vtx1.isFake() );
		  pv_trk_pvunbiased_IsFake_p2.push_back( vtx2.isFake() );
		  
		  pv_trk_pvunbiased_NTracks_p1.push_back( vtxTkCollection1.size() );
		  pv_trk_pvunbiased_NTracks_p2.push_back( vtxTkCollection2.size() );
		  
		  pv_trk_pvunbiased_SumTrackPt_p1.push_back( SumTrackPt_p1 );
		  pv_trk_pvunbiased_SumTrackPt_p2.push_back( SumTrackPt_p2 );
		  
		  pv_trk_pvunbiased_SumTrackPt2_p1.push_back( SumTrackPt2_p1 );	     
		  pv_trk_pvunbiased_SumTrackPt2_p2.push_back( SumTrackPt2_p2 );
		  
		  pv_trk_pvunbiased_fracHighPurity_p1.push_back( pv_fracHighPurity_p1 );
		  pv_trk_pvunbiased_fracHighPurity_p2.push_back( pv_fracHighPurity_p2 );
		  
		  pv_trk_pvunbiased_chi2_p1.push_back( vtx1.chi2() );
		  pv_trk_pvunbiased_chi2_p2.push_back( vtx2.chi2() );
		  
		  pv_trk_pvunbiased_ndof_p1.push_back( vtx1.ndof() );
		  pv_trk_pvunbiased_ndof_p2.push_back( vtx2.ndof() );
		  
		  pv_trk_pvunbiased_x_p1.push_back( vtx1.x()*micron );
		  pv_trk_pvunbiased_y_p1.push_back( vtx1.y()*micron );
		  pv_trk_pvunbiased_z_p1.push_back( vtx1.z()*micron );
		  pv_trk_pvunbiased_xError_p1.push_back( vtx1.xError()*micron );
		  pv_trk_pvunbiased_yError_p1.push_back( vtx1.yError()*micron );
		  pv_trk_pvunbiased_zError_p1.push_back( vtx1.zError()*micron );
		  
		  pv_trk_pvunbiased_x_p2.push_back( vtx2.x()*micron );
		  pv_trk_pvunbiased_y_p2.push_back( vtx2.y()*micron );
		  pv_trk_pvunbiased_z_p2.push_back( vtx2.z()*micron );
		  pv_trk_pvunbiased_xError_p2.push_back( vtx2.xError()*micron );
		  pv_trk_pvunbiased_yError_p2.push_back( vtx2.yError()*micron );
		  pv_trk_pvunbiased_zError_p2.push_back( vtx2.zError()*micron );

		  Track::Point vtxPositionUnbiased1 = Track::Point(vtx1.position().x(), vtx1.position().y(), vtx1.position().z());
		  Track::Point vtxPositionUnbiased2 = Track::Point(vtx2.position().x(), vtx2.position().y(), vtx2.position().z());
		  
		  pv_trk_d0_pvunbiased_p1.push_back( trk.dxy(vtxPositionUnbiased1) * micron );
		  pv_trk_dz_pvunbiased_p1.push_back( trk.dz(vtxPositionUnbiased1) * micron );
		  pv_trk_d0_bs_zpvunbiased_p1.push_back( trk.dxy(pvbeamspot->position(vtxPositionUnbiased1.z())) * micron );

		  pv_trk_d0_pvunbiased_p2.push_back( trk.dxy(vtxPositionUnbiased2) * micron );
		  pv_trk_dz_pvunbiased_p2.push_back( trk.dz(vtxPositionUnbiased2) * micron );
		  pv_trk_d0_bs_zpvunbiased_p2.push_back( trk.dxy(pvbeamspot->position(vtxPositionUnbiased2.z())) * micron );

		  if( doTruth && !runOnData )
		    {
		       RefToBase<Track> trkRef(trackViews, pv_trk_idx.back());
		       
		       auto matched = recSimCollTracks.find(trkRef);
		       
		       if( matched != recSimCollTracks.end() )
			 {
			    const TrackingParticleRef* tpPtr = &((matched->val)[0].first);
			    const TrackingParticleRef& tp = *tpPtr;

			    TrackingParticle::Point vertex = tp->vertex();
			    TrackingParticle::Vector momentum = tp->momentum();

			    TrackingParticle::Point vertexTP = parametersDefinerTP->vertex(iEvent, iSetup, tp);
			    TrackingParticle::Vector momentumTP = parametersDefinerTP->momentum(iEvent, iSetup, tp);
			    
			    GlobalPoint gp1 = GlobalPoint(vtx1.position().x(), vtx1.position().y(), vtx1.position().z());
			    GlobalPoint gp2 = GlobalPoint(vtx2.position().x(), vtx2.position().y(), vtx2.position().z());

			    pv_trk_mc_dxy_pvunbiased_p1.push_back( TrackingParticleIP::dxy(vertex, momentum, gp1) * micron );
			    pv_trk_mc_dz_pvunbiased_p1.push_back( TrackingParticleIP::dz(vertex, momentum, gp1) * micron );
			    pv_trk_mc_dxy_tp_pvunbiased_p1.push_back( TrackingParticleIP::dxy(vertexTP, momentumTP, gp1) * micron );
			    pv_trk_mc_dz_tp_pvunbiased_p1.push_back( TrackingParticleIP::dz(vertexTP, momentumTP, gp1) * micron );

			    pv_trk_mc_dxy_pvunbiased_p2.push_back( TrackingParticleIP::dxy(vertex, momentum, gp2) * micron );
			    pv_trk_mc_dz_pvunbiased_p2.push_back( TrackingParticleIP::dz(vertex, momentum, gp2) * micron );
			    pv_trk_mc_dxy_tp_pvunbiased_p2.push_back( TrackingParticleIP::dxy(vertexTP, momentumTP, gp2) * micron );
			    pv_trk_mc_dz_tp_pvunbiased_p2.push_back( TrackingParticleIP::dz(vertexTP, momentumTP, gp2) * micron );
			 }
		       else
			 {
			    pv_trk_mc_dxy_pvunbiased_p1.push_back( null );
			    pv_trk_mc_dz_pvunbiased_p1.push_back( null );
			    pv_trk_mc_dxy_tp_pvunbiased_p1.push_back( null );
			    pv_trk_mc_dz_tp_pvunbiased_p1.push_back( null );

			    pv_trk_mc_dxy_pvunbiased_p2.push_back( null );
			    pv_trk_mc_dz_pvunbiased_p2.push_back( null );
			    pv_trk_mc_dxy_tp_pvunbiased_p2.push_back( null );
			    pv_trk_mc_dz_tp_pvunbiased_p2.push_back( null );
			 }
		    }		  
		  else
		    {		       
		       pv_trk_mc_dxy_pvunbiased_p1.push_back( null );
		       pv_trk_mc_dz_pvunbiased_p1.push_back( null );
		       pv_trk_mc_dxy_tp_pvunbiased_p1.push_back( null );
		       pv_trk_mc_dz_tp_pvunbiased_p1.push_back( null );

		       pv_trk_mc_dxy_pvunbiased_p2.push_back( null );
		       pv_trk_mc_dz_pvunbiased_p2.push_back( null );
		       pv_trk_mc_dxy_tp_pvunbiased_p2.push_back( null );
		       pv_trk_mc_dz_tp_pvunbiased_p2.push_back( null );
		    }
	       }	     
	     else
	       {
		  pv_trk_pvunbiased_IsValid_p1.push_back( 0 );
		  pv_trk_pvunbiased_IsValid_p2.push_back( 0 );
		  
		  pv_trk_pvunbiased_IsFake_p1.push_back( 0 );
		  pv_trk_pvunbiased_IsFake_p2.push_back( 0 );
		  
		  pv_trk_pvunbiased_NTracks_p1.push_back( 0 );
		  pv_trk_pvunbiased_NTracks_p2.push_back( 0 );
		  
		  pv_trk_pvunbiased_SumTrackPt_p1.push_back( null );
		  pv_trk_pvunbiased_SumTrackPt_p2.push_back( null );
		  
		  pv_trk_pvunbiased_SumTrackPt2_p1.push_back( null );
		  pv_trk_pvunbiased_SumTrackPt2_p2.push_back( null );
		  
		  pv_trk_pvunbiased_fracHighPurity_p1.push_back( null );
		  pv_trk_pvunbiased_fracHighPurity_p2.push_back( null );
		  
		  pv_trk_pvunbiased_chi2_p1.push_back( null );
		  pv_trk_pvunbiased_chi2_p2.push_back( null );
		  
		  pv_trk_pvunbiased_ndof_p1.push_back( null );
		  pv_trk_pvunbiased_ndof_p2.push_back( null );
		  
		  pv_trk_pvunbiased_x_p1.push_back( null );
		  pv_trk_pvunbiased_y_p1.push_back( null );
		  pv_trk_pvunbiased_z_p1.push_back( null );
		  pv_trk_pvunbiased_xError_p1.push_back( null );
		  pv_trk_pvunbiased_yError_p1.push_back( null );
		  pv_trk_pvunbiased_zError_p1.push_back( null );
		  
		  pv_trk_pvunbiased_x_p2.push_back( null );
		  pv_trk_pvunbiased_y_p2.push_back( null );
		  pv_trk_pvunbiased_z_p2.push_back( null );
		  pv_trk_pvunbiased_xError_p2.push_back( null );
		  pv_trk_pvunbiased_yError_p2.push_back( null );
		  pv_trk_pvunbiased_zError_p2.push_back( null );

		  pv_trk_d0_pvunbiased_p1.push_back( null );
		  pv_trk_dz_pvunbiased_p1.push_back( null );
		  pv_trk_d0_bs_zpvunbiased_p1.push_back( null );

		  pv_trk_d0_pvunbiased_p2.push_back( null );
		  pv_trk_dz_pvunbiased_p2.push_back( null );
		  pv_trk_d0_bs_zpvunbiased_p2.push_back( null );
	       }
	     
	     if( !pvst.empty() )
	       {
		  reco::Vertex vtxt = reco::Vertex(pvst.front());

		  float unbiasedSumTrackPt = 0.;
		  float unbiasedSumTrackPt2 = 0.;
		  float unbiasedFracHighPurity = 0.;
		  
		  for( TrackCollection::const_iterator itt = newPVTkCollection.begin(); itt != newPVTkCollection.end(); itt++ )
		    {
		       unbiasedSumTrackPt += (*itt).pt();
		       unbiasedSumTrackPt2 += (*itt).pt()*(*itt).pt();
		       unbiasedFracHighPurity += (*itt).quality(reco::TrackBase::highPurity);
		    }		  
		  int nTracksUnbiased = vtxt.tracksSize();
		  if( nTracksUnbiased ) unbiasedFracHighPurity /= float(nTracksUnbiased);
		  
		  pv_trk_pvunbiased_IsValid.push_back( vtxt.isValid() );
		  pv_trk_pvunbiased_IsFake.push_back( vtxt.isFake() );		  
		  pv_trk_pvunbiased_NTracks.push_back( newPVTkCollection.size() );
		  pv_trk_pvunbiased_SumTrackPt.push_back( unbiasedSumTrackPt );
		  pv_trk_pvunbiased_SumTrackPt2.push_back( unbiasedSumTrackPt2 );
		  pv_trk_pvunbiased_fracHighPurity.push_back( unbiasedFracHighPurity );
		  pv_trk_pvunbiased_chi2.push_back( vtxt.chi2() );		  
		  pv_trk_pvunbiased_ndof.push_back( vtxt.ndof() );		  
		  pv_trk_pvunbiased_x.push_back( vtxt.x()*micron );
		  pv_trk_pvunbiased_y.push_back( vtxt.y()*micron );
		  pv_trk_pvunbiased_z.push_back( vtxt.z()*micron );
		  pv_trk_pvunbiased_xError.push_back( vtxt.xError()*micron );
		  pv_trk_pvunbiased_yError.push_back( vtxt.yError()*micron );
		  pv_trk_pvunbiased_zError.push_back( vtxt.zError()*micron );

		  Track::Point vtxPositionUnbiased = Track::Point(vtxt.position().x(), vtxt.position().y(), vtxt.position().z());
		  
		  pv_trk_d0_pvunbiased.push_back( trk.dxy(vtxPositionUnbiased) * micron );
		  pv_trk_dz_pvunbiased.push_back( trk.dz(vtxPositionUnbiased) * micron );
		  pv_trk_d0_bs_zpvunbiased.push_back( trk.dxy(pvbeamspot->position(vtxPositionUnbiased.z())) * micron );
		  
		  if( doTruth && !runOnData )
		    {
		       RefToBase<Track> trkRef(trackViews, pv_trk_idx.back());
		       
		       auto matched = recSimCollTracks.find(trkRef);
		       
		       if( matched != recSimCollTracks.end() )
			 {
			    const TrackingParticleRef* tpPtr = &((matched->val)[0].first);
			    const TrackingParticleRef& tp = *tpPtr;

			    TrackingParticle::Point vertex = tp->vertex();
			    TrackingParticle::Vector momentum = tp->momentum();

			    TrackingParticle::Point vertexTP = parametersDefinerTP->vertex(iEvent, iSetup, tp);
			    TrackingParticle::Vector momentumTP = parametersDefinerTP->momentum(iEvent, iSetup, tp);
			    
			    GlobalPoint gp = GlobalPoint(vtxt.position().x(), vtxt.position().y(), vtxt.position().z());
			    
			    pv_trk_mc_dxy_pvunbiased.push_back( TrackingParticleIP::dxy(vertex, momentum, gp) * micron );
			    pv_trk_mc_dz_pvunbiased.push_back( TrackingParticleIP::dz(vertex, momentum, gp) * micron );			    
			    pv_trk_mc_dxy_tp_pvunbiased.push_back( TrackingParticleIP::dxy(vertexTP, momentumTP, gp) * micron );
			    pv_trk_mc_dz_tp_pvunbiased.push_back( TrackingParticleIP::dz(vertexTP, momentumTP, gp) * micron );
			 }
		       else
			 {
			    pv_trk_mc_dxy_pvunbiased.push_back( null );
			    pv_trk_mc_dz_pvunbiased.push_back( null );			    
			    pv_trk_mc_dxy_tp_pvunbiased.push_back( null );
			    pv_trk_mc_dz_tp_pvunbiased.push_back( null );			    
			 }
		    }
		  else
		    {		       
		       pv_trk_mc_dxy_pvunbiased.push_back( null );
		       pv_trk_mc_dz_pvunbiased.push_back( null );
		       pv_trk_mc_dxy_tp_pvunbiased.push_back( null );
		       pv_trk_mc_dz_tp_pvunbiased.push_back( null );
		    }		  
	       }
	     else
	       {
		  pv_trk_pvunbiased_IsValid.push_back( 0 );
		  pv_trk_pvunbiased_IsFake.push_back( 0 );
		  pv_trk_pvunbiased_NTracks.push_back( null );
		  pv_trk_pvunbiased_SumTrackPt.push_back( null );
		  pv_trk_pvunbiased_SumTrackPt2.push_back( null );
		  pv_trk_pvunbiased_fracHighPurity.push_back( null );
		  pv_trk_pvunbiased_chi2.push_back( null );		  
		  pv_trk_pvunbiased_ndof.push_back( null );
		  pv_trk_pvunbiased_x.push_back( null );
		  pv_trk_pvunbiased_y.push_back( null );
		  pv_trk_pvunbiased_z.push_back( null );
		  pv_trk_pvunbiased_xError.push_back( null );
		  pv_trk_pvunbiased_yError.push_back( null );
		  pv_trk_pvunbiased_zError.push_back( null );
		  
		  pv_trk_d0_pvunbiased.push_back( null );
		  pv_trk_dz_pvunbiased.push_back( null );
		  pv_trk_d0_bs_zpvunbiased.push_back( null );
	       }	     
	     
	     pv_SumTrackPt += trk.pt();
	     pv_SumTrackPt2 += trk.pt()*trk.pt();
	     pv_fracHighPurity += trk.quality(reco::TrackBase::highPurity);
	     
	     if( pvs[ipv].hasTrackWeight() ) pv_trackWeight.push_back( pvs[ipv].trackWeight(*it) );
	     else pv_trackWeight.push_back( null );
	     
	     pv_trk_isHighPurity.push_back( trk.quality(reco::TrackBase::highPurity) );
	     pv_trk_algo.push_back( trk.algo() );
	     pv_trk_originalAlgo.push_back( trk.originalAlgo() );
	     
	     pv_trk_pt.push_back( trk.pt() );
	     pv_trk_px.push_back( trk.px() );
	     pv_trk_py.push_back( trk.py() );
	     pv_trk_pz.push_back( trk.pz() );
	     pv_trk_p.push_back( trk.p() );
	     pv_trk_eta.push_back( trk.eta() );
	     pv_trk_phi.push_back( trk.phi() );
	     
	     pv_trk_nTrackerLayers.push_back( trk.hitPattern().trackerLayersWithMeasurement() );
	     pv_trk_nPixelBarrelLayers.push_back( trk.hitPattern().pixelBarrelLayersWithMeasurement() );
	     pv_trk_nPixelEndcapLayers.push_back( trk.hitPattern().pixelEndcapLayersWithMeasurement() );
	     pv_trk_nStripLayers.push_back( trk.hitPattern().stripLayersWithMeasurement() );
	     
	     pv_trk_nValid.push_back( trk.numberOfValidHits() );
	     pv_trk_fValid.push_back( trk.validFraction() );
	     pv_trk_nValidTracker.push_back( trk.hitPattern().numberOfValidTrackerHits() );
	     pv_trk_nValidPixelBarrel.push_back( trk.hitPattern().numberOfValidPixelBarrelHits() );
	     pv_trk_nValidPixelEndcap.push_back( trk.hitPattern().numberOfValidPixelEndcapHits() );
	     pv_trk_nValidStrip.push_back( trk.hitPattern().numberOfValidStripHits() );
	     
	     pv_trk_nMissed.push_back( trk.numberOfLostHits() );
	     pv_trk_nMissedOut.push_back( trk.hitPattern().numberOfLostHits(HitPattern::MISSING_OUTER_HITS) );
	     pv_trk_nMissedIn.push_back( trk.hitPattern().numberOfLostHits(HitPattern::MISSING_INNER_HITS) );
	     pv_trk_nMissedTrackerOut.push_back( trk.hitPattern().numberOfLostTrackerHits(HitPattern::MISSING_OUTER_HITS) );
	     pv_trk_nMissedTrackerIn.push_back( trk.hitPattern().numberOfLostTrackerHits(HitPattern::MISSING_INNER_HITS) );
	     pv_trk_nMissedPixelBarrelOut.push_back( trk.hitPattern().numberOfLostPixelBarrelHits(HitPattern::MISSING_OUTER_HITS) );
	     pv_trk_nMissedPixelBarrelIn.push_back( trk.hitPattern().numberOfLostPixelBarrelHits(HitPattern::MISSING_INNER_HITS) );
	     pv_trk_nMissedPixelEndcapOut.push_back( trk.hitPattern().numberOfLostPixelEndcapHits(HitPattern::MISSING_OUTER_HITS) );
	     pv_trk_nMissedPixelEndcapIn.push_back( trk.hitPattern().numberOfLostPixelEndcapHits(HitPattern::MISSING_INNER_HITS) );
	     
	     pv_trk_hasPixelBarrelLayer1.push_back( trk.hitPattern().hasValidHitInPixelLayer(PixelSubdetector::SubDetector::PixelBarrel, 1) );
	     pv_trk_hasPixelEndcapLayer1.push_back( trk.hitPattern().hasValidHitInPixelLayer(PixelSubdetector::SubDetector::PixelEndcap, 1) );
	     pv_trk_hasPixelBarrelLayer2.push_back( trk.hitPattern().hasValidHitInPixelLayer(PixelSubdetector::SubDetector::PixelBarrel, 2) );
	     pv_trk_hasPixelEndcapLayer2.push_back( trk.hitPattern().hasValidHitInPixelLayer(PixelSubdetector::SubDetector::PixelEndcap, 2) );
	     pv_trk_hasPixelBarrelLayer3.push_back( trk.hitPattern().hasValidHitInPixelLayer(PixelSubdetector::SubDetector::PixelBarrel, 3) );
	     pv_trk_hasPixelEndcapLayer3.push_back( trk.hitPattern().hasValidHitInPixelLayer(PixelSubdetector::SubDetector::PixelEndcap, 3) );
	     pv_trk_hasPixelBarrelLayer4.push_back( trk.hitPattern().hasValidHitInPixelLayer(PixelSubdetector::SubDetector::PixelBarrel, 4) );
	     pv_trk_hasPixelEndcapLayer4.push_back( trk.hitPattern().hasValidHitInPixelLayer(PixelSubdetector::SubDetector::PixelEndcap, 4) );
	     
	     pv_trk_quality.push_back( trk.qualityMask() );
	     pv_trk_normalizedChi2.push_back( trk.normalizedChi2() );
	     pv_trk_ndof.push_back( trk.ndof() );
	     pv_trk_charge.push_back( trk.charge() );
	     pv_trk_qoverp.push_back( trk.qoverp() );
	     pv_trk_qoverpError.push_back( trk.qoverpError() );
	     pv_trk_theta.push_back( trk.theta() );
	     pv_trk_thetaError.push_back( trk.thetaError() );
	     pv_trk_lambda.push_back( trk.lambda() );
	     pv_trk_lambdaError.push_back( trk.lambdaError() );
	     pv_trk_ptError.push_back( trk.ptError() );
	     pv_trk_etaError.push_back( trk.etaError() );
	     pv_trk_phiError.push_back( trk.phiError() );
	     
	     pv_trk_d0.push_back( trk.dxy() * micron );
	     pv_trk_dz.push_back( trk.dz() * micron );
	     pv_trk_d0_pv.push_back( trk.dxy(vtxPosition) * micron );
	     pv_trk_dz_pv.push_back( trk.dz(vtxPosition) * micron );
	     pv_trk_d0_bs.push_back( trk.dxy(pvbeamspot->position()) * micron );
	     pv_trk_d0_bs_zpca.push_back( trk.dxy(*pvbeamspot) * micron );
	     pv_trk_d0_bs_zpv.push_back( trk.dxy(pvbeamspot->position(vtxPosition.z())) * micron );
	     pv_trk_dz_bs.push_back( trk.dz(pvbeamspot->position()) * micron );
	     pv_trk_d0Err.push_back( trk.d0Error() * micron );
	     pv_trk_dzErr.push_back( trk.dzError() * micron );
	     
	     if( doTruth && !runOnData && ftree->pv_mc_hasMatch[ipv] )
	       {
		  Track::Point tvPosition = Track::Point((ftree->pv_mc_x[ipv][0])/micron, (ftree->pv_mc_y[ipv][0])/micron, (ftree->pv_mc_z[ipv][0])/micron);
		  
		  pv_trk_d0_tv.push_back( trk.dxy(tvPosition) * micron );
		  pv_trk_dz_tv.push_back( trk.dz(tvPosition) * micron );
	       }
	     else
	       {
		  pv_trk_d0_tv.push_back( null );
		  pv_trk_dz_tv.push_back( null );
	       }
	     
	     iTrk++;
	  }
	if( nTracks ) pv_fracHighPurity /= float(nTracks);
   
	ftree->pv_IsValid.push_back( pvr[ipv].isValid() );
	ftree->pv_IsFake.push_back( pvr[ipv].isFake() );
	ftree->pv_NTracks.push_back( nTracks );
	ftree->pv_SumTrackPt.push_back( pv_SumTrackPt );
	ftree->pv_SumTrackPt2.push_back( pv_SumTrackPt2 );
	ftree->pv_fracHighPurity.push_back( pv_fracHighPurity );
	ftree->pv_chi2.push_back( pvr[ipv].chi2() );
	ftree->pv_ndof.push_back( pvr[ipv].ndof() );
	ftree->pv_x.push_back( pvr[ipv].x()*micron );
	ftree->pv_y.push_back( pvr[ipv].y()*micron );
	ftree->pv_z.push_back( pvr[ipv].z()*micron );
	ftree->pv_xError.push_back( pvr[ipv].xError()*micron );
	ftree->pv_yError.push_back( pvr[ipv].yError()*micron );
	ftree->pv_zError.push_back( pvr[ipv].zError()*micron );
	
	ftree->pv_trk_weight.push_back( pv_trackWeight );
	ftree->pv_trk_isHighPurity.push_back( pv_trk_isHighPurity );
	ftree->pv_trk_algo.push_back( pv_trk_algo );
	ftree->pv_trk_originalAlgo.push_back( pv_trk_originalAlgo );
	
	ftree->pv_trk_idx.push_back( pv_trk_idx );

	ftree->pv_trk_pvN.push_back( pv_trk_pvN );
	ftree->pv_trk_pv1N.push_back( pv_trk_pv1N );
	ftree->pv_trk_pv2N.push_back( pv_trk_pv2N );

	ftree->pv_trk_pvunbiased_IsValid.push_back( pv_trk_pvunbiased_IsValid );
	ftree->pv_trk_pvunbiased_IsFake.push_back( pv_trk_pvunbiased_IsFake );
	ftree->pv_trk_pvunbiased_NTracks.push_back( pv_trk_pvunbiased_NTracks );
	ftree->pv_trk_pvunbiased_SumTrackPt.push_back( pv_trk_pvunbiased_SumTrackPt );
	ftree->pv_trk_pvunbiased_SumTrackPt2.push_back( pv_trk_pvunbiased_SumTrackPt2 );
	ftree->pv_trk_pvunbiased_fracHighPurity.push_back( pv_trk_pvunbiased_fracHighPurity );
	ftree->pv_trk_pvunbiased_chi2.push_back( pv_trk_pvunbiased_chi2 );
	ftree->pv_trk_pvunbiased_ndof.push_back( pv_trk_pvunbiased_ndof );
	ftree->pv_trk_pvunbiased_x.push_back( pv_trk_pvunbiased_x );
	ftree->pv_trk_pvunbiased_y.push_back( pv_trk_pvunbiased_y );
	ftree->pv_trk_pvunbiased_z.push_back( pv_trk_pvunbiased_z );
	ftree->pv_trk_pvunbiased_xError.push_back( pv_trk_pvunbiased_xError );
	ftree->pv_trk_pvunbiased_yError.push_back( pv_trk_pvunbiased_yError );
	ftree->pv_trk_pvunbiased_zError.push_back( pv_trk_pvunbiased_zError );
		  
	ftree->pv_trk_d0_pvunbiased.push_back( pv_trk_d0_pvunbiased );
	ftree->pv_trk_dz_pvunbiased.push_back( pv_trk_dz_pvunbiased );
	ftree->pv_trk_d0_bs_zpvunbiased.push_back( pv_trk_d0_bs_zpvunbiased );

	ftree->pv_trk_mc_dxy_pvunbiased.push_back( pv_trk_mc_dxy_pvunbiased );
	ftree->pv_trk_mc_dz_pvunbiased.push_back( pv_trk_mc_dz_pvunbiased );	
	ftree->pv_trk_mc_dxy_tp_pvunbiased.push_back( pv_trk_mc_dxy_tp_pvunbiased );
	ftree->pv_trk_mc_dz_tp_pvunbiased.push_back( pv_trk_mc_dz_tp_pvunbiased );

	ftree->pv_trk_mc_dxy_pvunbiased_p1.push_back( pv_trk_mc_dxy_pvunbiased_p1 );
	ftree->pv_trk_mc_dz_pvunbiased_p1.push_back( pv_trk_mc_dz_pvunbiased_p1 );	
	ftree->pv_trk_mc_dxy_tp_pvunbiased_p1.push_back( pv_trk_mc_dxy_tp_pvunbiased_p1 );
	ftree->pv_trk_mc_dz_tp_pvunbiased_p1.push_back( pv_trk_mc_dz_tp_pvunbiased_p1 );

	ftree->pv_trk_mc_dxy_pvunbiased_p2.push_back( pv_trk_mc_dxy_pvunbiased_p2 );
	ftree->pv_trk_mc_dz_pvunbiased_p2.push_back( pv_trk_mc_dz_pvunbiased_p2 );	
	ftree->pv_trk_mc_dxy_tp_pvunbiased_p2.push_back( pv_trk_mc_dxy_tp_pvunbiased_p2 );
	ftree->pv_trk_mc_dz_tp_pvunbiased_p2.push_back( pv_trk_mc_dz_tp_pvunbiased_p2 );
	
	ftree->pv_trk_pvunbiased_IsValid_p1.push_back( pv_trk_pvunbiased_IsValid_p1 );
	ftree->pv_trk_pvunbiased_IsFake_p1.push_back( pv_trk_pvunbiased_IsFake_p1 );
	ftree->pv_trk_pvunbiased_NTracks_p1.push_back( pv_trk_pvunbiased_NTracks_p1 );
	ftree->pv_trk_pvunbiased_SumTrackPt_p1.push_back( pv_trk_pvunbiased_SumTrackPt_p1 );
	ftree->pv_trk_pvunbiased_SumTrackPt2_p1.push_back( pv_trk_pvunbiased_SumTrackPt2_p1 );
	ftree->pv_trk_pvunbiased_fracHighPurity_p1.push_back( pv_trk_pvunbiased_fracHighPurity_p1 );
	ftree->pv_trk_pvunbiased_chi2_p1.push_back( pv_trk_pvunbiased_chi2_p1 );
	ftree->pv_trk_pvunbiased_ndof_p1.push_back( pv_trk_pvunbiased_ndof_p1 );
	ftree->pv_trk_pvunbiased_x_p1.push_back( pv_trk_pvunbiased_x_p1 );
	ftree->pv_trk_pvunbiased_y_p1.push_back( pv_trk_pvunbiased_y_p1 );
	ftree->pv_trk_pvunbiased_z_p1.push_back( pv_trk_pvunbiased_z_p1 );
	ftree->pv_trk_pvunbiased_xError_p1.push_back( pv_trk_pvunbiased_xError_p1 );
	ftree->pv_trk_pvunbiased_yError_p1.push_back( pv_trk_pvunbiased_yError_p1 );
	ftree->pv_trk_pvunbiased_zError_p1.push_back( pv_trk_pvunbiased_zError_p1 );
		  
	ftree->pv_trk_d0_pvunbiased_p1.push_back( pv_trk_d0_pvunbiased_p1 );
	ftree->pv_trk_dz_pvunbiased_p1.push_back( pv_trk_dz_pvunbiased_p1 );
	ftree->pv_trk_d0_bs_zpvunbiased_p1.push_back( pv_trk_d0_bs_zpvunbiased_p1 );

	ftree->pv_trk_pvunbiased_IsValid_p2.push_back( pv_trk_pvunbiased_IsValid_p2 );
	ftree->pv_trk_pvunbiased_IsFake_p2.push_back( pv_trk_pvunbiased_IsFake_p2 );
	ftree->pv_trk_pvunbiased_NTracks_p2.push_back( pv_trk_pvunbiased_NTracks_p2 );
	ftree->pv_trk_pvunbiased_SumTrackPt_p2.push_back( pv_trk_pvunbiased_SumTrackPt_p2 );
	ftree->pv_trk_pvunbiased_SumTrackPt2_p2.push_back( pv_trk_pvunbiased_SumTrackPt2_p2 );
	ftree->pv_trk_pvunbiased_fracHighPurity_p2.push_back( pv_trk_pvunbiased_fracHighPurity_p2 );
	ftree->pv_trk_pvunbiased_chi2_p2.push_back( pv_trk_pvunbiased_chi2_p2 );
	ftree->pv_trk_pvunbiased_ndof_p2.push_back( pv_trk_pvunbiased_ndof_p2 );
	ftree->pv_trk_pvunbiased_x_p2.push_back( pv_trk_pvunbiased_x_p2 );
	ftree->pv_trk_pvunbiased_y_p2.push_back( pv_trk_pvunbiased_y_p2 );
	ftree->pv_trk_pvunbiased_z_p2.push_back( pv_trk_pvunbiased_z_p2 );
	ftree->pv_trk_pvunbiased_xError_p2.push_back( pv_trk_pvunbiased_xError_p2 );
	ftree->pv_trk_pvunbiased_yError_p2.push_back( pv_trk_pvunbiased_yError_p2 );
	ftree->pv_trk_pvunbiased_zError_p2.push_back( pv_trk_pvunbiased_zError_p2 );
		  
	ftree->pv_trk_d0_pvunbiased_p2.push_back( pv_trk_d0_pvunbiased_p2 );
	ftree->pv_trk_dz_pvunbiased_p2.push_back( pv_trk_dz_pvunbiased_p2 );
	ftree->pv_trk_d0_bs_zpvunbiased_p2.push_back( pv_trk_d0_bs_zpvunbiased_p2 );
	
	ftree->pv_trk_pt.push_back( pv_trk_pt );
	ftree->pv_trk_px.push_back( pv_trk_px );
	ftree->pv_trk_py.push_back( pv_trk_py );
	ftree->pv_trk_pz.push_back( pv_trk_pz );
	ftree->pv_trk_p.push_back( pv_trk_p );
	ftree->pv_trk_eta.push_back( pv_trk_eta );
	ftree->pv_trk_phi.push_back( pv_trk_phi );

	ftree->pv_trk_nTrackerLayers.push_back( pv_trk_nTrackerLayers );
	ftree->pv_trk_nPixelBarrelLayers.push_back( pv_trk_nPixelBarrelLayers );
	ftree->pv_trk_nPixelEndcapLayers.push_back( pv_trk_nPixelEndcapLayers );
	ftree->pv_trk_nStripLayers.push_back( pv_trk_nStripLayers );
	     
	ftree->pv_trk_nValid.push_back( pv_trk_nValid );
	ftree->pv_trk_fValid.push_back( pv_trk_fValid );
	ftree->pv_trk_nValidTracker.push_back( pv_trk_nValidTracker );
	ftree->pv_trk_nValidPixelBarrel.push_back( pv_trk_nValidPixelBarrel );
	ftree->pv_trk_nValidPixelEndcap.push_back( pv_trk_nValidPixelEndcap );
	ftree->pv_trk_nValidStrip.push_back( pv_trk_nValidStrip );
	     
	ftree->pv_trk_nMissed.push_back( pv_trk_nMissed );
	ftree->pv_trk_nMissedOut.push_back( pv_trk_nMissedOut );
	ftree->pv_trk_nMissedIn.push_back( pv_trk_nMissedIn );
	ftree->pv_trk_nMissedTrackerOut.push_back( pv_trk_nMissedTrackerOut );
	ftree->pv_trk_nMissedTrackerIn.push_back( pv_trk_nMissedTrackerIn );
	ftree->pv_trk_nMissedPixelBarrelOut.push_back( pv_trk_nMissedPixelBarrelOut );
	ftree->pv_trk_nMissedPixelBarrelIn.push_back( pv_trk_nMissedPixelBarrelIn );
	ftree->pv_trk_nMissedPixelEndcapOut.push_back( pv_trk_nMissedPixelEndcapOut );
	ftree->pv_trk_nMissedPixelEndcapIn.push_back( pv_trk_nMissedPixelEndcapIn );
	     
	ftree->pv_trk_hasPixelBarrelLayer1.push_back( pv_trk_hasPixelBarrelLayer1 );
	ftree->pv_trk_hasPixelEndcapLayer1.push_back( pv_trk_hasPixelEndcapLayer1 );
	ftree->pv_trk_hasPixelBarrelLayer2.push_back( pv_trk_hasPixelBarrelLayer2 );
	ftree->pv_trk_hasPixelEndcapLayer2.push_back( pv_trk_hasPixelEndcapLayer2 );
	ftree->pv_trk_hasPixelBarrelLayer3.push_back( pv_trk_hasPixelBarrelLayer3 );
	ftree->pv_trk_hasPixelEndcapLayer3.push_back( pv_trk_hasPixelEndcapLayer3 );
	ftree->pv_trk_hasPixelBarrelLayer4.push_back( pv_trk_hasPixelBarrelLayer4 );
	ftree->pv_trk_hasPixelEndcapLayer4.push_back( pv_trk_hasPixelEndcapLayer4 );
	     
	ftree->pv_trk_quality.push_back( pv_trk_quality );
	ftree->pv_trk_normalizedChi2.push_back( pv_trk_normalizedChi2 );
	ftree->pv_trk_ndof.push_back( pv_trk_ndof );
	ftree->pv_trk_charge.push_back( pv_trk_charge );
	ftree->pv_trk_qoverp.push_back( pv_trk_qoverp );
	ftree->pv_trk_qoverpError.push_back( pv_trk_qoverpError );
	ftree->pv_trk_theta.push_back( pv_trk_theta );
	ftree->pv_trk_thetaError.push_back( pv_trk_thetaError );
	ftree->pv_trk_lambda.push_back( pv_trk_lambda );
	ftree->pv_trk_lambdaError.push_back( pv_trk_lambdaError );
	ftree->pv_trk_ptError.push_back( pv_trk_ptError );
	ftree->pv_trk_etaError.push_back( pv_trk_etaError );
	ftree->pv_trk_phiError.push_back( pv_trk_phiError );
	
	ftree->pv_trk_d0.push_back( pv_trk_d0 );
	ftree->pv_trk_dz.push_back( pv_trk_dz );
	ftree->pv_trk_d0_pv.push_back( pv_trk_d0_pv );
	ftree->pv_trk_dz_pv.push_back( pv_trk_dz_pv );
	ftree->pv_trk_d0_bs.push_back( pv_trk_d0_bs );
	ftree->pv_trk_d0_bs_zpca.push_back( pv_trk_d0_bs_zpca );
	ftree->pv_trk_d0_bs_zpv.push_back( pv_trk_d0_bs_zpv );
	ftree->pv_trk_dz_bs.push_back( pv_trk_dz_bs );
	ftree->pv_trk_d0Err.push_back( pv_trk_d0Err );
	ftree->pv_trk_dzErr.push_back( pv_trk_dzErr );
	
	ftree->pv_trk_d0_tv.push_back( pv_trk_d0_tv );
	ftree->pv_trk_dz_tv.push_back( pv_trk_dz_tv );
     }   

   // Vertex split method
   for( unsigned int ipv=0;ipv<pvr.size();ipv++ )
     {	     	   
	std::vector<reco::TransientTrack> vtxTracks = pvs[ipv].originalTracks();
	stable_sort(vtxTracks.begin(), vtxTracks.end(), sortPt);
	
	int nTracks = pvr[ipv].tracksSize();
	
	if( nTracks < PVFitNTracksMin || int(ipv) > PVFitNMax ) break;

	reco::TrackCollection vtxTkCollection1;
	reco::TrackCollection vtxTkCollection2;
	
	std::vector<int> vtxTkIdx1;
	std::vector<int> vtxTkIdx2;

	float SumTrackPt_p1 = 0;
	float SumTrackPt2_p1 = 0;
	
	float SumTrackPt_p2 = 0;
	float SumTrackPt2_p2 = 0;
	
	float pv_fracHighPurity_p1 = 0;
	float pv_fracHighPurity_p2 = 0;
	
	int iTrk = 0;
	
	for( std::vector<reco::TransientTrack>::const_iterator it = vtxTracks.begin(); it != vtxTracks.end(); it++ )
	  {
	     reco::Track trk = (*it).track();
	     
	     if( rnd->Rndm() > 0.5 )
	       {
		  vtxTkCollection1.push_back(trk);
		  SumTrackPt_p1 += trk.pt();
		  SumTrackPt2_p1 += trk.pt()*trk.pt();
		  pv_fracHighPurity_p1 += trk.quality(reco::TrackBase::highPurity);
		  vtxTkIdx1.push_back(iTrk);
	       }	
	     else
	       {	     
		  vtxTkCollection2.push_back(trk);
		  SumTrackPt_p2 += trk.pt();
		  SumTrackPt2_p2 += trk.pt()*trk.pt();
		  pv_fracHighPurity_p2 += trk.quality(reco::TrackBase::highPurity);
		  vtxTkIdx2.push_back(iTrk);
	       }
	     
	     iTrk++;
	  }
	if( nTracks ) 
	  {
	     pv_fracHighPurity_p1 /= float(nTracks);
	     pv_fracHighPurity_p2 /= float(nTracks);
	  }	

	vector<TransientVertex> pvs1 = revertex.makeVertices(vtxTkCollection1, *pvbeamspot, iSetup);
	vector<TransientVertex> pvs2 = revertex.makeVertices(vtxTkCollection2, *pvbeamspot, iSetup);
	
	if( !pvs1.empty() && !pvs2.empty() )
	  {
	     reco::Vertex vtx1 = reco::Vertex(pvs1.front());
	     reco::Vertex vtx2 = reco::Vertex(pvs2.front());
	     
	     ftree->pv_IsValid_p1.push_back( vtx1.isValid() );
	     ftree->pv_IsValid_p2.push_back( vtx2.isValid() );
	
	     ftree->pv_IsFake_p1.push_back( vtx1.isFake() );
	     ftree->pv_IsFake_p2.push_back( vtx2.isFake() );
	     
	     ftree->pv_NTracks_p1.push_back( vtxTkCollection1.size() );
	     ftree->pv_NTracks_p2.push_back( vtxTkCollection2.size() );
	     
	     ftree->pv_SumTrackPt_p1.push_back( SumTrackPt_p1 );
	     ftree->pv_SumTrackPt_p2.push_back( SumTrackPt_p2 );

	     ftree->pv_SumTrackPt2_p1.push_back( SumTrackPt2_p1 );	     
	     ftree->pv_SumTrackPt2_p2.push_back( SumTrackPt2_p2 );
	     
	     ftree->pv_fracHighPurity_p1.push_back( pv_fracHighPurity_p1 );
	     ftree->pv_fracHighPurity_p2.push_back( pv_fracHighPurity_p2 );

	     ftree->pv_vtxTkIdx_p1.push_back( vtxTkIdx1 );
	     ftree->pv_vtxTkIdx_p2.push_back( vtxTkIdx2 );
	     
	     ftree->pv_chi2_p1.push_back( vtx1.chi2() );
	     ftree->pv_chi2_p2.push_back( vtx2.chi2() );
	     
	     ftree->pv_ndof_p1.push_back( vtx1.ndof() );
	     ftree->pv_ndof_p2.push_back( vtx2.ndof() );
	     
	     ftree->pv_x_p1.push_back( vtx1.x()*micron );
	     ftree->pv_y_p1.push_back( vtx1.y()*micron );
	     ftree->pv_z_p1.push_back( vtx1.z()*micron );
	     ftree->pv_xError_p1.push_back( vtx1.xError()*micron );
	     ftree->pv_yError_p1.push_back( vtx1.yError()*micron );
	     ftree->pv_zError_p1.push_back( vtx1.zError()*micron );
	     
	     ftree->pv_x_p2.push_back( vtx2.x()*micron );
	     ftree->pv_y_p2.push_back( vtx2.y()*micron );
	     ftree->pv_z_p2.push_back( vtx2.z()*micron );
	     ftree->pv_xError_p2.push_back( vtx2.xError()*micron );
	     ftree->pv_yError_p2.push_back( vtx2.yError()*micron );
	     ftree->pv_zError_p2.push_back( vtx2.zError()*micron );
	  }
	else
	  {
	     ftree->pv_IsValid_p1.push_back( 0 );
	     ftree->pv_IsValid_p2.push_back( 0 );
	
	     ftree->pv_IsFake_p1.push_back( 0 );
	     ftree->pv_IsFake_p2.push_back( 0 );
	     
	     ftree->pv_NTracks_p1.push_back( 0 );
	     ftree->pv_NTracks_p2.push_back( 0 );
	     
	     ftree->pv_SumTrackPt_p1.push_back( null );
	     ftree->pv_SumTrackPt_p2.push_back( null );

	     ftree->pv_SumTrackPt2_p1.push_back( null );
	     ftree->pv_SumTrackPt2_p2.push_back( null );

	     ftree->pv_fracHighPurity_p1.push_back( null );
	     ftree->pv_fracHighPurity_p2.push_back( null );
	     
	     ftree->pv_vtxTkIdx_p1.push_back( vtxTkIdx1 );
	     ftree->pv_vtxTkIdx_p2.push_back( vtxTkIdx2 );
	     
	     ftree->pv_chi2_p1.push_back( null );
	     ftree->pv_chi2_p2.push_back( null );
	     
	     ftree->pv_ndof_p1.push_back( null );
	     ftree->pv_ndof_p2.push_back( null );
	     
	     ftree->pv_x_p1.push_back( null );
	     ftree->pv_y_p1.push_back( null );
	     ftree->pv_z_p1.push_back( null );
	     ftree->pv_xError_p1.push_back( null );
	     ftree->pv_yError_p1.push_back( null );
	     ftree->pv_zError_p1.push_back( null );
	     
	     ftree->pv_x_p2.push_back( null );
	     ftree->pv_y_p2.push_back( null );
	     ftree->pv_z_p2.push_back( null );
	     ftree->pv_xError_p2.push_back( null );
	     ftree->pv_yError_p2.push_back( null );
	     ftree->pv_zError_p2.push_back( null );
	  }	
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
   
   edm::ESHandle<TransientTrackBuilder> theB;
   iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder", theB);
   
   std::cout << "Tracks = " << nTracks << ", use = " << int(float(nTracks)/float(trackScale)) << std::endl;
   
   for( TrackCollection::const_iterator itk = tracks->begin(); itk != tracks->end(); ++itk )
     {
	if( rnd->Rndm() > trackProb && trackScale > 0 ) continue;
	
	// --- track selection ---
//	if( ! trackSelection(*itk) ) continue;
	// ---
	
	TrackCollection newTkCollection;
	newTkCollection.assign(tracks->begin(), itk);
	newTkCollection.insert(newTkCollection.end(), itk+1, tracks->end());

	//newTkCollection.insert(newTkCollection.end(),itk,tracks->end()); // only for debugging purpose

	//cout << "tracks before,after size: " << tracks->size() << " , " << newTkCollection.size() << endl;

	// Refit the primary vertex
	vector<TransientVertex> pvs = revertex.makeVertices(newTkCollection, *pvbeamspot, iSetup);
	//cout << "vertices before,after: " << vtxH->size() << " , " << pvs.size() << endl;
	
	if( pvs.empty() ) continue;

	reco::Vertex newPV = reco::Vertex(pvs.front());
	Track::Point vtxPosition = Track::Point(newPV.position().x(), newPV.position().y(), newPV.position().z());
	
	if( ! vertexSelection(newPV) ) continue;

	if( doTruth && !runOnData )
	  {	     
	     RefToBase<Track> trkRef(trackViews, itk - tracks->begin());

	     auto matched = recSimCollTracks.find(trkRef);
	     
	     bool trk_mc_hasMatch = 0;
	     std::vector<float> trk_mc_matchQuality;
	     
	     std::vector<int> trk_mc_pdgId;
	     std::vector<int> trk_mc_origin;
	     std::vector<int> trk_mc_status;

	     std::vector<float> trk_mc_pt;
	     std::vector<float> trk_mc_px;
	     std::vector<float> trk_mc_py;
	     std::vector<float> trk_mc_pz;
	     std::vector<float> trk_mc_E;
	     std::vector<float> trk_mc_p;
	     std::vector<float> trk_mc_eta;
	     std::vector<float> trk_mc_phi;

	     std::vector<int> trk_mc_numberOfHits;
	     std::vector<int> trk_mc_numberOfTrackerHits;
	     std::vector<int> trk_mc_numberOfTrackerLayers;

	     std::vector<float> trk_mc_dxy_center;
	     std::vector<float> trk_mc_dz_center;
	     std::vector<float> trk_mc_dxy_pv;
	     std::vector<float> trk_mc_dz_pv;
	     std::vector<float> trk_mc_dxy_bs;
	     std::vector<float> trk_mc_dz_bs;
	     
	     std::vector<float> trk_mc_dxy_tp_center;
	     std::vector<float> trk_mc_dz_tp_center;
	     std::vector<float> trk_mc_dxy_tp_pv;
	     std::vector<float> trk_mc_dz_tp_pv;
	     std::vector<float> trk_mc_dxy_tp_bs;
	     std::vector<float> trk_mc_dz_tp_bs;

	     std::vector<float> trk_mc_vtx_x;
	     std::vector<float> trk_mc_vtx_y;
	     std::vector<float> trk_mc_vtx_z;
	     std::vector<float> trk_mc_vtx_pca_x;
	     std::vector<float> trk_mc_vtx_pca_y;
	     std::vector<float> trk_mc_vtx_pca_z;
	     
	     std::vector<bool> trk_mc_isFake;
	     std::vector<bool> trk_mc_isBad;
	     std::vector<bool> trk_mc_isBadInnerHits;
	     std::vector<bool> trk_mc_isSharedInnerHits;
	     std::vector<bool> trk_mc_isSignalEvent;
	     std::vector<bool> trk_mc_isTrackerSimHits;
	     std::vector<bool> trk_mc_isBottom;
	     std::vector<bool> trk_mc_isCharm;
	     std::vector<bool> trk_mc_isLight;
	     std::vector<bool> trk_mc_isMuon;
		       
	     std::vector<bool> trk_mc_isBWeakDecay;
	     std::vector<bool> trk_mc_isCWeakDecay;
	     std::vector<bool> trk_mc_isChargePionDecay;
	     std::vector<bool> trk_mc_isChargeKaonDecay;
	     std::vector<bool> trk_mc_isTauDecay;
	     std::vector<bool> trk_mc_isKsDecay;
	     std::vector<bool> trk_mc_isLambdaDecay;
	     std::vector<bool> trk_mc_isJpsiDecay;
	     std::vector<bool> trk_mc_isXiDecay;
	     std::vector<bool> trk_mc_isOmegaDecay;
	     std::vector<bool> trk_mc_isSigmaPlusDecay;
	     std::vector<bool> trk_mc_isSigmaMinusDecay;
	     std::vector<bool> trk_mc_isLongLivedDecay;
	     
	     std::vector<bool> trk_mc_isKnownProcess;
	     std::vector<bool> trk_mc_isUndefinedProcess;
	     std::vector<bool> trk_mc_isUnknownProcess;
	     std::vector<bool> trk_mc_isPrimaryProcess;
	     std::vector<bool> trk_mc_isHadronicProcess;
	     std::vector<bool> trk_mc_isDecayProcess;
	     std::vector<bool> trk_mc_isComptonProcess;
	     std::vector<bool> trk_mc_isAnnihilationProcess;
	     std::vector<bool> trk_mc_isEIoniProcess;
	     std::vector<bool> trk_mc_isHIoniProcess;
	     std::vector<bool> trk_mc_isMuIoniProcess;
	     std::vector<bool> trk_mc_isPhotonProcess;
	     std::vector<bool> trk_mc_isMuPairProdProcess;
	     std::vector<bool> trk_mc_isConversionsProcess;
	     std::vector<bool> trk_mc_isEBremProcess;
	     std::vector<bool> trk_mc_isSynchrotronRadiationProcess;
	     std::vector<bool> trk_mc_isMuBremProcess;
	     std::vector<bool> trk_mc_isMuNuclProcess;
	     
	     std::vector<bool> trk_mc_isFromBWeakDecayMuon;
	     std::vector<bool> trk_mc_isFromCWeakDecayMuon;
	     std::vector<bool> trk_mc_isDecayOnFlightMuon;
	     std::vector<bool> trk_mc_isFromChargePionMuon;
	     std::vector<bool> trk_mc_isFromChargeKaonMuon;
	     
	     std::vector<bool> trk_mc_isPrimaryVertex;
	     std::vector<bool> trk_mc_isSecondaryVertex;
	     std::vector<bool> trk_mc_isTertiaryVertex;
	     
	     std::vector<bool> trk_mc_isUnknown;
	     
	     if( matched != recSimCollTracks.end() )
	       {
		  trk_mc_hasMatch = 1;
		 
		  for(const auto trkRefQuality: matched->val)
		    {		       
		       const TrackingParticleRef* tpPtr = &(trkRefQuality.first);
		       const TrackingParticleRef& tp = *tpPtr;
		       trk_mc_matchQuality.push_back( trkRefQuality.second );
		       
		       trk_mc_pt.push_back( tp->pt() );
		       trk_mc_px.push_back( tp->px() );
		       trk_mc_py.push_back( tp->py() );
		       trk_mc_pz.push_back( tp->pz() );
		       trk_mc_E.push_back( tp->energy() );
		       trk_mc_p.push_back( tp->p() );
		       trk_mc_eta.push_back( tp->eta() );
		       trk_mc_phi.push_back( tp->phi() );

		       trk_mc_numberOfHits.push_back( tp->numberOfHits() );
		       trk_mc_numberOfTrackerHits.push_back( tp->numberOfTrackerHits() );
		       trk_mc_numberOfTrackerLayers.push_back( tp->numberOfTrackerLayers() );
		       
		       trk_mc_pdgId.push_back( tp->pdgId() ); // pdgId of gen particle, otherwise of the first sim track
		       trk_mc_status.push_back( tp->status() ); // status code of gen particle, otherwise returns -99
		       
		       // TP properties at production
		       TrackingParticle::Point vertex = tp->vertex();
		       trk_mc_vtx_x.push_back( vertex.x() * micron );
		       trk_mc_vtx_y.push_back( vertex.y() * micron );
		       trk_mc_vtx_z.push_back( vertex.z() * micron );

		       TrackingParticle::Vector momentum = tp->momentum();
		       
		       trk_mc_dxy_center.push_back( TrackingParticleIP::dxy(vertex, momentum, Track::Point(0., 0., 0.)) * micron );
		       trk_mc_dz_center.push_back( TrackingParticleIP::dz(vertex, momentum, Track::Point(0., 0., 0.)) * micron );

		       GlobalPoint newPVgp = GlobalPoint(newPV.position().x(), newPV.position().y(), newPV.position().z());
		       
		       trk_mc_dxy_pv.push_back( TrackingParticleIP::dxy(vertex, momentum, newPVgp) * micron );
		       trk_mc_dz_pv.push_back( TrackingParticleIP::dz(vertex, momentum, newPVgp) * micron );

		       GlobalPoint BSgp = GlobalPoint(pvbeamspot->position().x(), pvbeamspot->position().y(), pvbeamspot->position().z());
		       
		       trk_mc_dxy_bs.push_back( TrackingParticleIP::dxy(vertex, momentum, BSgp) * micron );
		       trk_mc_dz_bs.push_back( TrackingParticleIP::dz(vertex, momentum, BSgp) * micron );
		       
		       // TP properties at PCA wrt BS
		       TrackingParticle::Point vertexTP = parametersDefinerTP->vertex(iEvent, iSetup, tp);
		       trk_mc_vtx_pca_x.push_back( vertexTP.x() * micron );
		       trk_mc_vtx_pca_y.push_back( vertexTP.y() * micron );
		       trk_mc_vtx_pca_z.push_back( vertexTP.z() * micron );
		       
		       TrackingParticle::Vector momentumTP = parametersDefinerTP->momentum(iEvent, iSetup, tp);

		       trk_mc_dxy_tp_center.push_back( TrackingParticleIP::dxy(vertexTP, momentumTP, Track::Point(0., 0., 0.)) * micron );
		       trk_mc_dz_tp_center.push_back( TrackingParticleIP::dz(vertexTP, momentumTP, Track::Point(0., 0., 0.)) * micron );
		       
		       trk_mc_dxy_tp_pv.push_back( TrackingParticleIP::dxy(vertexTP, momentumTP, newPVgp) * micron );
		       trk_mc_dz_tp_pv.push_back( TrackingParticleIP::dz(vertexTP, momentumTP, newPVgp) * micron );
		       
		       trk_mc_dxy_tp_bs.push_back( TrackingParticleIP::dxy(vertexTP, momentumTP, BSgp) * micron );
		       trk_mc_dz_tp_bs.push_back( TrackingParticleIP::dz(vertexTP, momentumTP, BSgp) * micron );
		       
		       TrackingVertexRef tv(trackingVertex, 0); // only consider the highest sum-pT^2 vertex
		       if( tp->parentVertex().get() != tv.get() )
			 {
			    if( tp->genParticles().size() ) trk_mc_origin.push_back(0); // Non prompt long lived pythia particle
			    else trk_mc_origin.push_back(1); // Geant4 particle
			 }
		       else trk_mc_origin.push_back(2); // Genuinely prompt
		       
		       trkClassifier_.evaluate(*tpPtr);
		       
		       trk_mc_isFake.push_back( trkClassifier_.is(TrackCategories::Fake) ); // no match no any sim track
		       trk_mc_isBad.push_back( trkClassifier_.is(TrackCategories::Bad) ); // has a large d0 pull
		       trk_mc_isBadInnerHits.push_back( trkClassifier_.is(TrackCategories::BadInnerHits) );
		       trk_mc_isSharedInnerHits.push_back( trkClassifier_.is(TrackCategories::SharedInnerHits) );
		       trk_mc_isSignalEvent.push_back( trkClassifier_.is(TrackCategories::SignalEvent) );
		       trk_mc_isTrackerSimHits.push_back( trkClassifier_.is(TrackCategories::TrackerSimHits) );
		       trk_mc_isBottom.push_back( trkClassifier_.is(TrackCategories::Bottom) );
		       trk_mc_isCharm.push_back( trkClassifier_.is(TrackCategories::Charm) );
		       trk_mc_isLight.push_back( trkClassifier_.is(TrackCategories::Light) );
		       trk_mc_isMuon.push_back( trkClassifier_.is(TrackCategories::Muon) );
		       
		       trk_mc_isBWeakDecay.push_back( trkClassifier_.is(TrackCategories::BWeakDecay) );
		       trk_mc_isCWeakDecay.push_back( trkClassifier_.is(TrackCategories::CWeakDecay) );
		       trk_mc_isChargePionDecay.push_back( trkClassifier_.is(TrackCategories::ChargePionDecay) );
		       trk_mc_isChargeKaonDecay.push_back( trkClassifier_.is(TrackCategories::ChargeKaonDecay) );
		       trk_mc_isTauDecay.push_back( trkClassifier_.is(TrackCategories::TauDecay) );
		       trk_mc_isKsDecay.push_back( trkClassifier_.is(TrackCategories::KsDecay) );
		       trk_mc_isLambdaDecay.push_back( trkClassifier_.is(TrackCategories::LambdaDecay) );
		       trk_mc_isJpsiDecay.push_back( trkClassifier_.is(TrackCategories::JpsiDecay) );
		       trk_mc_isXiDecay.push_back( trkClassifier_.is(TrackCategories::XiDecay) );
		       trk_mc_isOmegaDecay.push_back( trkClassifier_.is(TrackCategories::OmegaDecay) );
		       trk_mc_isSigmaPlusDecay.push_back( trkClassifier_.is(TrackCategories::SigmaPlusDecay) );
		       trk_mc_isSigmaMinusDecay.push_back( trkClassifier_.is(TrackCategories::SigmaMinusDecay) );
		       trk_mc_isLongLivedDecay.push_back( trkClassifier_.is(TrackCategories::LongLivedDecay) );
		       
		       trk_mc_isKnownProcess.push_back( trkClassifier_.is(TrackCategories::KnownProcess) );
		       trk_mc_isUndefinedProcess.push_back( trkClassifier_.is(TrackCategories::UndefinedProcess) );
		       trk_mc_isUnknownProcess.push_back( trkClassifier_.is(TrackCategories::UnknownProcess) );
		       trk_mc_isPrimaryProcess.push_back( trkClassifier_.is(TrackCategories::PrimaryProcess) );
		       trk_mc_isHadronicProcess.push_back( trkClassifier_.is(TrackCategories::HadronicProcess) );
		       trk_mc_isDecayProcess.push_back( trkClassifier_.is(TrackCategories::DecayProcess) );
		       trk_mc_isComptonProcess.push_back( trkClassifier_.is(TrackCategories::ComptonProcess) );
		       trk_mc_isAnnihilationProcess.push_back( trkClassifier_.is(TrackCategories::AnnihilationProcess) );
		       trk_mc_isEIoniProcess.push_back( trkClassifier_.is(TrackCategories::EIoniProcess) );
		       trk_mc_isHIoniProcess.push_back( trkClassifier_.is(TrackCategories::HIoniProcess) );
		       trk_mc_isMuIoniProcess.push_back( trkClassifier_.is(TrackCategories::MuIoniProcess) );
		       trk_mc_isPhotonProcess.push_back( trkClassifier_.is(TrackCategories::PhotonProcess) );
		       trk_mc_isMuPairProdProcess.push_back( trkClassifier_.is(TrackCategories::MuPairProdProcess) );
		       trk_mc_isConversionsProcess.push_back( trkClassifier_.is(TrackCategories::ConversionsProcess) );
		       trk_mc_isEBremProcess.push_back( trkClassifier_.is(TrackCategories::EBremProcess) );
		       trk_mc_isSynchrotronRadiationProcess.push_back( trkClassifier_.is(TrackCategories::SynchrotronRadiationProcess) );
		       trk_mc_isMuBremProcess.push_back( trkClassifier_.is(TrackCategories::MuBremProcess) );
		       trk_mc_isMuNuclProcess.push_back( trkClassifier_.is(TrackCategories::MuNuclProcess) );
		       
		       trk_mc_isFromBWeakDecayMuon.push_back( trkClassifier_.is(TrackCategories::FromBWeakDecayMuon) );
		       trk_mc_isFromCWeakDecayMuon.push_back( trkClassifier_.is(TrackCategories::FromCWeakDecayMuon) );
		       trk_mc_isDecayOnFlightMuon.push_back( trkClassifier_.is(TrackCategories::DecayOnFlightMuon) );
		       trk_mc_isFromChargePionMuon.push_back( trkClassifier_.is(TrackCategories::FromChargePionMuon) );
		       trk_mc_isFromChargeKaonMuon.push_back( trkClassifier_.is(TrackCategories::FromChargeKaonMuon) );
		       
		       trk_mc_isPrimaryVertex.push_back( trkClassifier_.is(TrackCategories::PrimaryVertex) );
		       trk_mc_isSecondaryVertex.push_back( trkClassifier_.is(TrackCategories::SecondaryVertex) );
		       trk_mc_isTertiaryVertex.push_back( trkClassifier_.is(TrackCategories::TertiaryVertex) );
		       
		       trk_mc_isUnknown.push_back( trkClassifier_.is(TrackCategories::Unknown) );
		    }
	       }
	     
	     ftree->trk_mc_hasMatch.push_back( trk_mc_hasMatch );
	     ftree->trk_mc_matchQuality.push_back( trk_mc_matchQuality );
	     
	     ftree->trk_mc_pdgId.push_back( trk_mc_pdgId );
	     ftree->trk_mc_origin.push_back( trk_mc_origin );
	     ftree->trk_mc_status.push_back( trk_mc_status );

	     ftree->trk_mc_numberOfHits.push_back( trk_mc_numberOfHits );
	     ftree->trk_mc_numberOfTrackerHits.push_back( trk_mc_numberOfTrackerHits );
	     ftree->trk_mc_numberOfTrackerLayers.push_back( trk_mc_numberOfTrackerLayers );
	     
	     ftree->trk_mc_pt.push_back( trk_mc_pt );
	     ftree->trk_mc_px.push_back( trk_mc_px );
	     ftree->trk_mc_py.push_back( trk_mc_py );
	     ftree->trk_mc_pz.push_back( trk_mc_pz );
	     ftree->trk_mc_E.push_back( trk_mc_E );
	     ftree->trk_mc_p.push_back( trk_mc_p );
	     ftree->trk_mc_eta.push_back( trk_mc_eta );
	     ftree->trk_mc_phi.push_back( trk_mc_phi );

	     ftree->trk_mc_dxy_center.push_back( trk_mc_dxy_center );
	     ftree->trk_mc_dz_center.push_back( trk_mc_dz_center );
	     ftree->trk_mc_dxy_pv.push_back( trk_mc_dxy_pv );
	     ftree->trk_mc_dz_pv.push_back( trk_mc_dz_pv );
	     ftree->trk_mc_dxy_bs.push_back( trk_mc_dxy_bs );
	     ftree->trk_mc_dz_bs.push_back( trk_mc_dz_bs );
	     
	     ftree->trk_mc_dxy_tp_center.push_back( trk_mc_dxy_tp_center );
	     ftree->trk_mc_dz_tp_center.push_back( trk_mc_dz_tp_center );
	     ftree->trk_mc_dxy_tp_pv.push_back( trk_mc_dxy_tp_pv );
	     ftree->trk_mc_dz_tp_pv.push_back( trk_mc_dz_tp_pv );
	     ftree->trk_mc_dxy_tp_bs.push_back( trk_mc_dxy_tp_bs );
	     ftree->trk_mc_dz_tp_bs.push_back( trk_mc_dz_tp_bs );
	     
	     ftree->trk_mc_vtx_x.push_back( trk_mc_vtx_x );
	     ftree->trk_mc_vtx_y.push_back( trk_mc_vtx_y );
	     ftree->trk_mc_vtx_z.push_back( trk_mc_vtx_z );
	     ftree->trk_mc_vtx_pca_x.push_back( trk_mc_vtx_pca_x );
	     ftree->trk_mc_vtx_pca_y.push_back( trk_mc_vtx_pca_y );
	     ftree->trk_mc_vtx_pca_z.push_back( trk_mc_vtx_pca_z );
	     
	     ftree->trk_mc_isFake.push_back( trk_mc_isFake );
	     ftree->trk_mc_isBad.push_back( trk_mc_isBad );
	     ftree->trk_mc_isBadInnerHits.push_back( trk_mc_isBadInnerHits );
	     ftree->trk_mc_isSharedInnerHits.push_back( trk_mc_isSharedInnerHits );
	     ftree->trk_mc_isSignalEvent.push_back( trk_mc_isSignalEvent );
	     ftree->trk_mc_isTrackerSimHits.push_back( trk_mc_isTrackerSimHits );
	     ftree->trk_mc_isBottom.push_back( trk_mc_isBottom );
	     ftree->trk_mc_isCharm.push_back( trk_mc_isCharm );
	     ftree->trk_mc_isLight.push_back( trk_mc_isLight );
	     ftree->trk_mc_isMuon.push_back( trk_mc_isMuon );
	     
	     ftree->trk_mc_isBWeakDecay.push_back( trk_mc_isBWeakDecay );
	     ftree->trk_mc_isCWeakDecay.push_back( trk_mc_isCWeakDecay );
	     ftree->trk_mc_isChargePionDecay.push_back( trk_mc_isChargePionDecay );
	     ftree->trk_mc_isChargeKaonDecay.push_back( trk_mc_isChargeKaonDecay );
	     ftree->trk_mc_isTauDecay.push_back( trk_mc_isTauDecay );
	     ftree->trk_mc_isKsDecay.push_back( trk_mc_isKsDecay );
	     ftree->trk_mc_isLambdaDecay.push_back( trk_mc_isLambdaDecay );
	     ftree->trk_mc_isJpsiDecay.push_back( trk_mc_isJpsiDecay );
	     ftree->trk_mc_isXiDecay.push_back( trk_mc_isXiDecay );
	     ftree->trk_mc_isOmegaDecay.push_back( trk_mc_isOmegaDecay );
	     ftree->trk_mc_isSigmaPlusDecay.push_back( trk_mc_isSigmaPlusDecay );
	     ftree->trk_mc_isSigmaMinusDecay.push_back( trk_mc_isSigmaMinusDecay );
	     ftree->trk_mc_isLongLivedDecay.push_back( trk_mc_isLongLivedDecay );
	     
	     ftree->trk_mc_isKnownProcess.push_back( trk_mc_isKnownProcess );
	     ftree->trk_mc_isUndefinedProcess.push_back( trk_mc_isUndefinedProcess );
	     ftree->trk_mc_isUnknownProcess.push_back( trk_mc_isUnknownProcess );
	     ftree->trk_mc_isPrimaryProcess.push_back( trk_mc_isPrimaryProcess );
	     ftree->trk_mc_isHadronicProcess.push_back( trk_mc_isHadronicProcess );
	     ftree->trk_mc_isDecayProcess.push_back( trk_mc_isDecayProcess );
	     ftree->trk_mc_isComptonProcess.push_back( trk_mc_isComptonProcess );
	     ftree->trk_mc_isAnnihilationProcess.push_back( trk_mc_isAnnihilationProcess );
	     ftree->trk_mc_isEIoniProcess.push_back( trk_mc_isEIoniProcess );
	     ftree->trk_mc_isHIoniProcess.push_back( trk_mc_isHIoniProcess );
	     ftree->trk_mc_isMuIoniProcess.push_back( trk_mc_isMuIoniProcess );
	     ftree->trk_mc_isPhotonProcess.push_back( trk_mc_isPhotonProcess );
	     ftree->trk_mc_isMuPairProdProcess.push_back( trk_mc_isMuPairProdProcess );
	     ftree->trk_mc_isConversionsProcess.push_back( trk_mc_isConversionsProcess );
	     ftree->trk_mc_isEBremProcess.push_back( trk_mc_isEBremProcess );
	     ftree->trk_mc_isSynchrotronRadiationProcess.push_back( trk_mc_isSynchrotronRadiationProcess );
	     ftree->trk_mc_isMuBremProcess.push_back( trk_mc_isMuBremProcess );
	     ftree->trk_mc_isMuNuclProcess.push_back( trk_mc_isMuNuclProcess );
	     
	     ftree->trk_mc_isFromBWeakDecayMuon.push_back( trk_mc_isFromBWeakDecayMuon );
	     ftree->trk_mc_isFromCWeakDecayMuon.push_back( trk_mc_isFromCWeakDecayMuon );
	     ftree->trk_mc_isDecayOnFlightMuon.push_back( trk_mc_isDecayOnFlightMuon );
	     ftree->trk_mc_isFromChargePionMuon.push_back( trk_mc_isFromChargePionMuon );
	     ftree->trk_mc_isFromChargeKaonMuon.push_back( trk_mc_isFromChargeKaonMuon );
	     
	     ftree->trk_mc_isPrimaryVertex.push_back( trk_mc_isPrimaryVertex );
	     ftree->trk_mc_isSecondaryVertex.push_back( trk_mc_isSecondaryVertex );
	     ftree->trk_mc_isTertiaryVertex.push_back( trk_mc_isTertiaryVertex );
	     
	     ftree->trk_mc_isUnknown.push_back( trk_mc_isUnknown );
	  }
	
	// Reco track
	ftree->trk_pt.push_back( itk->pt() );
	ftree->trk_px.push_back( itk->px() );
	ftree->trk_py.push_back( itk->py() );
	ftree->trk_pz.push_back( itk->pz() );
	ftree->trk_p.push_back( itk->p() );
	ftree->trk_eta.push_back( itk->eta() );
	ftree->trk_phi.push_back( itk->phi() );
	
	ftree->trk_idx.push_back( itk - tracks->begin() );
	
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
	ftree->trk_isHighPurity.push_back( itk->quality(reco::TrackBase::highPurity) );
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
	
	// Impact parameters are given at the innermost (reference) point on track (TrackBase)
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

	if( doTruth && !runOnData && ftree->ev_nPV > 0 && ftree->pv_mc_hasMatch[0] )
	  {
	     Track::Point tvPosition = Track::Point((ftree->pv_mc_x[0][0])/micron, (ftree->pv_mc_y[0][0])/micron, (ftree->pv_mc_z[0][0])/micron);
	     
	     ftree->trk_d0_tv.push_back( itk->dxy(tvPosition) * micron );
	     ftree->trk_dz_tv.push_back( itk->dz(tvPosition) * micron );
	  }
	else
	  {
	     ftree->trk_d0_tv.push_back( null );
	     ftree->trk_dz_tv.push_back( null );
	  }
	
/*	reco::TransientTrack tranitk = (*theB).build(*itk);

	GlobalPoint pvPos(vtx.position().x(), vtx.position().y(), vtx.position().z());
	GlobalPoint bsPos(pvbeamspot->position().x(), pvbeamspot->position().y(), pvbeamspot->position().z());
	GlobalVector dir(itk->px(), itk->py(), itk->pz());
	
	reco::Vertex bspv(pvbeamspot->position(), pvbeamspot->rotatedCovariance3D());
	
	float d0_perigee_pv = -tranitk.trajectoryStateClosestToPoint(pvPos).perigeeParameters().transverseImpactParameter() * micron;
	float d0_perigee_bs = -tranitk.trajectoryStateClosestToPoint(bsPos).perigeeParameters().transverseImpactParameter() * micron;
	float dz_perigee_pv = tranitk.trajectoryStateClosestToPoint(pvPos).perigeeParameters().longitudinalImpactParameter() * micron;
	float dz_perigee_bs = tranitk.trajectoryStateClosestToPoint(bsPos).perigeeParameters().longitudinalImpactParameter() * micron;
	
	float d0_IPTools_pv = IPTools::signedTransverseImpactParameter(tranitk,dir,vtx).second.value() * micron;
	float d0Err_IPTools_pv = IPTools::signedTransverseImpactParameter(tranitk,dir,vtx).second.error() * micron;
	float d0Sign_IPTools_pv = IPTools::signedTransverseImpactParameter(tranitk,dir,vtx).second.significance();
	float d0_IPTools_bs = IPTools::signedTransverseImpactParameter(tranitk,dir,bspv).second.value() * micron;
	float d0Err_IPTools_bs = IPTools::signedTransverseImpactParameter(tranitk,dir,bspv).second.error() * micron;
	float d0Sign_IPTools_bs = IPTools::signedTransverseImpactParameter(tranitk,dir,bspv).second.significance();
	
	float ip3d_IPTools_pv = IPTools::signedImpactParameter3D(tranitk,dir,vtx).second.value() * micron;
	float ip3dErr_IPTools_pv = IPTools::signedImpactParameter3D(tranitk,dir,vtx).second.error() * micron;
	float ip3dSign_IPTools_pv = IPTools::signedImpactParameter3D(tranitk,dir,vtx).second.significance();
	float ip3d_IPTools_bs = IPTools::signedImpactParameter3D(tranitk,dir,bspv).second.value() * micron;
	float ip3dErr_IPTools_bs = IPTools::signedImpactParameter3D(tranitk,dir,bspv).second.error() * micron;
	float ip3dSign_IPTools_bs = IPTools::signedImpactParameter3D(tranitk,dir,bspv).second.significance();
	
	float ip3dl_IPTools_pv = IPTools::linearizedSignedImpactParameter3D(tranitk,dir,vtx).second.value() * micron;
	float ip3dlErr_IPTools_pv = IPTools::linearizedSignedImpactParameter3D(tranitk,dir,vtx).second.error() * micron;
	float ip3dlSign_IPTools_pv = IPTools::linearizedSignedImpactParameter3D(tranitk,dir,vtx).second.significance();
	float ip3dl_IPTools_bs = IPTools::linearizedSignedImpactParameter3D(tranitk,dir,bspv).second.value() * micron;
	float ip3dlErr_IPTools_bs = IPTools::linearizedSignedImpactParameter3D(tranitk,dir,bspv).second.error() * micron;
	float ip3dlSign_IPTools_bs = IPTools::linearizedSignedImpactParameter3D(tranitk,dir,bspv).second.significance();*/
	
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
