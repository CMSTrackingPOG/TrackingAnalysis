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

#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"

#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"
 
#include "FWCore/Framework/interface/ESHandle.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/TrackReco/interface/HitPattern.h"

#include <TrackingTools/TrajectoryState/interface/PerigeeConversions.h>
#include <TrackingTools/TrajectoryState/interface/TrajectoryStateClosestToPoint.h>
#include <TrackingTools/PatternTools/interface/TSCPBuilderNoMaterial.h>

#include "IpResoStudies/EDAnalyzers/interface/VertexReProducer.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include <CommonTools/UtilAlgos/interface/TFileService.h>

#include "TH1F.h"
#include "TTree.h"
#include "TRandom3.h"

struct treeVarIP
{
   double pt;
   double p;
   double eta;
   double phi;
   int nXLayers;
   int nMissedOut;
   int nMissedIn;
   int hasPXL;
   int quality;
   double d0;
   double dz;
   double d0Err;
   double dzErr;
   double d0NoRefit;
   double dzNoRefit;
   double d0ErrNoRefit;
   double dzErrNoRefit;

   bool pvIsValid;
   bool pvIsFake;
   int pvNTracks;
   double pvSumTrackPt;
   double pvchi2;
   int pvndof;
   double pvx;
   double pvy;
   double pvz;
   double pvxError;
   double pvyError;
   double pvzError;
};

struct treeVarPV
{
   bool pvIsValid;
   bool pvIsFake;
   int pvNTracks;
   double pvSumTrackPt;
   double pvchi2;
   int pvndof;
   double pvx;
   double pvy;
   double pvz;
   double pvxError;
   double pvyError;
   double pvzError;
   
   bool pv1IsValid;
   bool pv1IsFake;
   int pv1NTracks;
   double pv1SumTrackPt;
   double pv1chi2;
   int pv1ndof;
   double pv1x;
   double pv1y;
   double pv1z;
   double pv1xError;
   double pv1yError;
   double pv1zError;
   
   bool pv2IsValid;
   bool pv2IsFake;
   int pv2NTracks;
   double pv2SumTrackPt;
   double pv2chi2;
   int pv2ndof;
   double pv2x;
   double pv2y;
   double pv2z;
   double pv2xError;
   double pv2yError;
   double pv2zError;
};

class Residuals : public edm::EDAnalyzer 
{  
   
 public:
   explicit Residuals(const edm::ParameterSet& pset);
   ~Residuals();
    
 private:
   virtual void beginJob() ;
   virtual void analyze(const edm::Event&, const edm::EventSetup&);
   virtual void endJob() ;

   bool trackSelection(const reco::Track& track) const;
   bool vertexSelection(const reco::Vertex& vertex) const;
  
   // ----------member data ---------------------------
   edm::EDGetTokenT<reco::VertexCollection> thePVToken_;
   edm::EDGetTokenT<reco::TrackCollection> theTracksToken_;
   edm::EDGetTokenT<reco::BeamSpot> theBeamspotToken_;

   // --- track selection variables
   double tkMinPt;
   int tkMinXLayers,tkMaxMissedOuterLayers,tkMaxMissedInnerLayers;

   // --- vertex selection variables
   unsigned int vtxTracksSizeMin;  
   unsigned int vtxTracksSizeMax;  
   double vtxErrorXMin,vtxErrorXMax;
   double vtxErrorYMin,vtxErrorYMax;
   double vtxErrorZMin,vtxErrorZMax;

   TH1F *h_d0;
   TTree *treeIP;
   TTree *treePV;
   treeVarIP trip;
   treeVarPV trpv;
   TRandom3 *rnd;
};

Residuals::Residuals(const edm::ParameterSet& pset)
{
   edm::InputTag TrackCollectionTag_ = pset.getParameter<edm::InputTag>("TrackLabel");
   theTracksToken_= consumes<reco::TrackCollection>(TrackCollectionTag_);
   
   edm::InputTag VertexCollectionTag_ = pset.getParameter<edm::InputTag>("VertexLabel");
   thePVToken_ = consumes<reco::VertexCollection>(VertexCollectionTag_);
   
   edm::InputTag BeamspotTag_ = edm::InputTag("offlineBeamSpot");
   theBeamspotToken_ = consumes<reco::BeamSpot>(BeamspotTag_);
   
   tkMinPt = pset.getParameter<double>("TkMinPt");    
   tkMinXLayers = pset.getParameter<int>("TkMinXLayers");
   tkMaxMissedOuterLayers = pset.getParameter<int>("TkMaxMissedOuterLayers");
   tkMaxMissedInnerLayers = pset.getParameter<int>("TkMaxMissedInnerLayers");
   
   vtxTracksSizeMin = pset.getParameter<int>("VtxTracksSizeMin");
   vtxTracksSizeMax = pset.getParameter<int>("VtxTracksSizeMax");
   vtxErrorXMin     = pset.getParameter<double>("VtxErrorXMin");
   vtxErrorXMax     = pset.getParameter<double>("VtxErrorXMax");
   vtxErrorYMin     = pset.getParameter<double>("VtxErrorYMin");
   vtxErrorYMax     = pset.getParameter<double>("VtxErrorYMax");
   vtxErrorZMin     = pset.getParameter<double>("VtxErrorZMin");
   vtxErrorZMax     = pset.getParameter<double>("VtxErrorZMax");

   rnd = new TRandom3();
   
   edm::Service<TFileService> fs;
   
   treeIP = fs->make<TTree>("treeIP", "recoTrack IP residuals");
   
   treeIP->Branch("pt", &trip.pt, "pt/D");
   treeIP->Branch("p", &trip.p, "p/D");
   treeIP->Branch("eta", &trip.eta, "eta/D");
   treeIP->Branch("phi", &trip.phi, "phi/D");
   treeIP->Branch("nXLayers", &trip.nXLayers, "nXLayers/I");
   treeIP->Branch("nMissedOut", &trip.nMissedOut, "nMissedOut/I");
   treeIP->Branch("nMissedIn", &trip.nMissedIn, "nMissedIn/I");
   treeIP->Branch("hasPXL", &trip.hasPXL, "hasPXL/I");
   treeIP->Branch("quality", &trip.quality, "quality/I");
   treeIP->Branch("d0", &trip.d0, "d0/D");
   treeIP->Branch("dz", &trip.dz, "dz/D");
   treeIP->Branch("d0Err", &trip.d0Err, "d0Err/D");
   treeIP->Branch("dzErr", &trip.dzErr, "dzErr/D");
   treeIP->Branch("d0NoRefit", &trip.d0NoRefit, "d0NoRefit/D");
   treeIP->Branch("dzNoRefit", &trip.dzNoRefit, "dzNoRefit/D");
   treeIP->Branch("d0ErrNoRefit", &trip.d0ErrNoRefit, "d0ErrNoRefit/D");
   treeIP->Branch("dzErrNoRefit", &trip.dzErrNoRefit, "dzErrNoRefit/D");

   treeIP->Branch("pvIsValid", &trip.pvIsValid, "pvIsValid/O");
   treeIP->Branch("pvIsFake", &trip.pvIsFake, "pvIsFake/O");
   treeIP->Branch("pvNTracks", &trip.pvNTracks, "pvNTracks/I");
   treeIP->Branch("pvSumTrackPt", &trip.pvSumTrackPt, "pvSumTrackPt/D");
   treeIP->Branch("pvchi2", &trip.pvchi2, "pvchi2/D");
   treeIP->Branch("pvndof", &trip.pvndof, "pvndof/I");
   treeIP->Branch("pvx", &trip.pvx, "pvx/D");
   treeIP->Branch("pvy", &trip.pvy, "pvy/D");
   treeIP->Branch("pvz", &trip.pvz, "pvz/D");
   treeIP->Branch("pvxError", &trip.pvxError, "pvxError/D");
   treeIP->Branch("pvyError", &trip.pvyError, "pvyError/D");
   treeIP->Branch("pvzError", &trip.pvzError, "pvzError/D");
   
   treePV = fs->make<TTree>("treePV", "PV resolution study");
   
   treePV->Branch("pvIsValid", &trpv.pvIsValid, "pvIsValid/O");
   treePV->Branch("pvIsFake", &trpv.pvIsFake, "pvIsFake/O");
   treePV->Branch("pvNTracks", &trpv.pvNTracks, "pvNTracks/I");
   treePV->Branch("pvSumTrackPt", &trpv.pvSumTrackPt, "pvSumTrackPt/D");
   treePV->Branch("pvchi2", &trpv.pvchi2, "pvchi2/D");
   treePV->Branch("pvndof", &trpv.pvndof, "pvndof/I");
   treePV->Branch("pvx", &trpv.pvx, "pvx/D");
   treePV->Branch("pvy", &trpv.pvy, "pvy/D");
   treePV->Branch("pvz", &trpv.pvz, "pvz/D");
   treePV->Branch("pvxError", &trpv.pvxError, "pvxError/D");
   treePV->Branch("pvyError", &trpv.pvyError, "pvyError/D");
   treePV->Branch("pvzError", &trpv.pvzError, "pvzError/D");
   
   treePV->Branch("pv1IsValid", &trpv.pv1IsValid, "pv1IsValid/O");
   treePV->Branch("pv1IsFake", &trpv.pv1IsFake, "pv1IsFake/O");
   treePV->Branch("pv1NTracks", &trpv.pv1NTracks, "pv1NTracks/I");
   treePV->Branch("pv1SumTrackPt", &trpv.pv1SumTrackPt, "pv1SumTrackPt/D");
   treePV->Branch("pv1chi2", &trpv.pv1chi2, "pv1chi2/D");
   treePV->Branch("pv1ndof", &trpv.pv1ndof, "pv1ndof/I");
   treePV->Branch("pv1x", &trpv.pv1x, "pv1x/D");
   treePV->Branch("pv1y", &trpv.pv1y, "pv1y/D");
   treePV->Branch("pv1z", &trpv.pv1z, "pv1z/D");
   treePV->Branch("pv1xError", &trpv.pv1xError, "pv1xError/D");
   treePV->Branch("pv1yError", &trpv.pv1yError, "pv1yError/D");
   treePV->Branch("pv1zError", &trpv.pv1zError, "pv1zError/D");
   
   treePV->Branch("pv2IsValid", &trpv.pv2IsValid, "pv2IsValid/O");
   treePV->Branch("pv2IsFake", &trpv.pv2IsFake, "pv2IsFake/O");
   treePV->Branch("pv2NTracks", &trpv.pv2NTracks, "pv2NTracks/I");
   treePV->Branch("pv2SumTrackPt", &trpv.pv2SumTrackPt, "pv2SumTrackPt/D");
   treePV->Branch("pv2chi2", &trpv.pv2chi2, "pv2chi2/D");
   treePV->Branch("pv2ndof", &trpv.pv2ndof, "pv2ndof/I");
   treePV->Branch("pv2x", &trpv.pv2x, "pv2x/D");
   treePV->Branch("pv2y", &trpv.pv2y, "pv2y/D");
   treePV->Branch("pv2z", &trpv.pv2z, "pv2z/D");
   treePV->Branch("pv2xError", &trpv.pv2xError, "pv2xError/D");
   treePV->Branch("pv2yError", &trpv.pv2yError, "pv2yError/D");
   treePV->Branch("pv2zError", &trpv.pv2zError, "pv2zError/D");
}

Residuals::~Residuals()
{
   delete rnd;
}

void Residuals::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;
   using namespace reco;
   using namespace std;

   std::cout << "point1" << std::endl;
   Handle<TrackCollection> tracks;
   iEvent.getByToken(theTracksToken_,tracks);
   std::cout << "point2" << std::endl;
   
//   ESHandle<TransientTrackBuilder> transTrackBuilder;
//   iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder", transTrackBuilder);
   
   //std::cout << "size of track collection: "<< tracks->size() << std::endl;
   
   std::cout << "point3" << std::endl;
   Handle<VertexCollection> vtxH;
   iEvent.getByToken(thePVToken_, vtxH);
   std::cout << "point4 " << thePVToken_ < 

   if (!vtxH.isValid()) return;

   //std::cout << "size of vtx collection: "<< vtxH->size() << std::endl;

   ESHandle<MagneticField> theMF;
   iSetup.get<IdealMagneticFieldRecord>().get(theMF);

   VertexReProducer revertex(vtxH, iEvent);

   //Handle<TrackCollection> pvtracks;
   //iEvent.getByLabel(revertex.inputTracks(),   pvtracks);
   //Handle<BeamSpot>        pvbeamspot;
   //iEvent.getByLabel(revertex.inputBeamSpot(), pvbeamspot);

   Handle<BeamSpot> pvbeamspot;
   iEvent.getByToken(theBeamspotToken_, pvbeamspot);

   //if(tracks.id() != pvtracks.id())
   // cout << "WARNING: the tracks originally used for PV are not the same used in this analyzer."
   //	  << "Is this really what you want?" << endl;

   //if (pvbeamspot.id() != theBeamSpot.id()) 
   //  edm::LogWarning("Inconsistency") << "The BeamSpot used for PV reco is not the same used in this analyzer.";

   double micron = 10000;
   
   // ------ check if the vertex is good enough -------
   if( vtxH->size() == 0 ) return;
   if( ! vertexSelection(vtxH->front()) ) return;
   // -------------------------------------------------

   /*
   cout << "original vtx x,y,z: " 
	<< vtxH->front().position().x() << " , "
	<< vtxH->front().position().y() << " , "
	<< vtxH->front().position().z() << endl;
   */

   const reco::Vertex vtx = vtxH->front();

   double trackSumPt = 0;   
   for( std::vector<reco::TrackBaseRef>::const_iterator it = vtx.tracks_begin(); it != vtx.tracks_end(); it++ )
     trackSumPt += (*it)->pt();
   
   trip.pvIsValid = vtx.isValid();
   trip.pvIsFake = vtx.isFake();	
   trip.pvNTracks = vtx.tracksSize();	
   trip.pvSumTrackPt = trackSumPt;	
   trip.pvchi2 = vtx.chi2();	
   trip.pvndof = vtx.ndof();   
   trip.pvx = vtx.x()*micron;
   trip.pvy = vtx.y()*micron;
   trip.pvz = vtx.z()*micron;
   trip.pvxError = vtx.xError()*micron;
   trip.pvyError = vtx.yError()*micron;
   trip.pvzError = vtx.zError()*micron;
   
   for( TrackCollection::const_iterator itk = tracks->begin(); itk != tracks->end(); ++itk )
     {
	// --- track selection ---
	if( ! trackSelection(*itk) ) continue;
	// ---
     
	TrackCollection newTkCollection;
	newTkCollection.assign(tracks->begin(), itk);
	newTkCollection.insert(newTkCollection.end(),itk+1,tracks->end());

	//newTkCollection.insert(newTkCollection.end(),itk,tracks->end()); // only for debugging purpose

	//cout << "tracks before,after size: " << tracks->size() << " , " << newTkCollection.size() << endl;
  
	// --- from Giovanni to refit the prim.vertex     
	vector<TransientVertex> pvs = revertex.makeVertices(newTkCollection, *pvbeamspot, 1, iSetup);
	//cout << "vertices before,after: " << vtxH->size() << " , " << pvs.size() << endl;

	if( pvs.empty() ) continue;

	reco::Vertex newPV = reco::Vertex(pvs.front());
	Track::Point vtxPosition = Track::Point(newPV.position().x(),
						newPV.position().y(),
						newPV.position().z());
	// ---
	if( ! vertexSelection(newPV) ) continue;

	double d0 = itk->dxy(vtxPosition);
	double dz = itk->dz(vtxPosition);

	double d0NoRefit = itk->dxy(vtxH->front().position());
	double dzNoRefit = itk->dz(vtxH->front().position());
	
	//cout << "d0:" << d0 << " dz:" << dz << std::endl;

	//Filling the tree
	trip.pt  = itk->pt();
	trip.p   = itk->p();
	trip.eta = itk->eta();
	trip.phi = itk->phi();
	trip.nXLayers   = itk->hitPattern().trackerLayersWithMeasurement();
	trip.nMissedOut = itk->hitPattern().numberOfLostHits(HitPattern::MISSING_OUTER_HITS);
	trip.nMissedIn  = itk->hitPattern().numberOfLostHits(HitPattern::MISSING_INNER_HITS);
	trip.hasPXL     = (itk->hitPattern().hasValidHitInPixelLayer(PixelSubdetector::SubDetector::PixelBarrel, 1) || 
			   itk->hitPattern().hasValidHitInPixelLayer(PixelSubdetector::SubDetector::PixelEndcap, 1));
	trip.quality = itk->qualityMask();
	trip.d0 = d0*micron;
	trip.dz = dz*micron;
	trip.d0Err = itk->d0Error()*micron;
	trip.dzErr = itk->dzError()*micron;
	trip.d0NoRefit = d0NoRefit*micron;
	trip.dzNoRefit = dzNoRefit*micron;
	trip.d0ErrNoRefit = itk->d0Error()*micron;
	trip.dzErrNoRefit = itk->dzError()*micron;
	
	treeIP->Fill();
     }

   std::vector<reco::TrackBaseRef> vtxTkCollection1;
   std::vector<reco::TrackBaseRef> vtxTkCollection2;

   double trackSumPt1 = 0;
   double trackSumPt2 = 0;
   
   for( std::vector<reco::TrackBaseRef>::const_iterator it = vtx.tracks_begin(); it != vtx.tracks_end(); it++ )
     {
	if( rnd->Rndm() > 0.5 )
	  {
	     vtxTkCollection1.push_back(*it);
	     trackSumPt1 += (*it)->pt();
	  }	
	else
	  {	     
	     vtxTkCollection2.push_back(*it);
	     trackSumPt2 += (*it)->pt();
	  }	
     }

   vector<TransientVertex> pvs1 = revertex.makeVertices(vtxTkCollection1, *pvbeamspot, 1, iSetup);
   vector<TransientVertex> pvs2 = revertex.makeVertices(vtxTkCollection2, *pvbeamspot, 1, iSetup);
//   vector<TransientVertex> pvs1 = revertex.makeVertices(vtxTkCollection1, *pvbeamspot, 0, iSetup);
//   vector<TransientVertex> pvs2 = revertex.makeVertices(vtxTkCollection2, *pvbeamspot, 0, iSetup);

   if( !pvs1.empty() && !pvs2.empty() )
     {	
	reco::Vertex vtx1 = reco::Vertex(pvs1.front());
	reco::Vertex vtx2 = reco::Vertex(pvs2.front());
	
	trpv.pvIsValid = vtx.isValid();
	trpv.pv1IsValid = vtx1.isValid();
	trpv.pv2IsValid = vtx2.isValid();
	
	trpv.pvIsFake = vtx.isFake();
	trpv.pv1IsFake = vtx1.isFake();
	trpv.pv2IsFake = vtx2.isFake();
	
	trpv.pvNTracks = vtx.tracksSize();
	trpv.pv1NTracks = vtxTkCollection1.size();
	trpv.pv2NTracks = vtxTkCollection2.size();
	
	trpv.pvSumTrackPt = trackSumPt;
	trpv.pv1SumTrackPt = trackSumPt1;
	trpv.pv2SumTrackPt = trackSumPt2;
	
	trpv.pvchi2 = vtx.chi2();
	trpv.pv1chi2 = vtx1.chi2();
	trpv.pv2chi2 = vtx2.chi2();
	
	trpv.pvndof = vtx.ndof();
	trpv.pv1ndof = vtx1.ndof();
	trpv.pv2ndof = vtx2.ndof();
	
	trpv.pvx = vtx.x()*micron;
	trpv.pvy = vtx.y()*micron;
	trpv.pvz = vtx.z()*micron;
	trpv.pvxError = vtx.xError()*micron;
	trpv.pvyError = vtx.yError()*micron;
	trpv.pvzError = vtx.zError()*micron;
	
	trpv.pv1x = vtx1.x()*micron;
	trpv.pv1y = vtx1.y()*micron;
	trpv.pv1z = vtx1.z()*micron;
	trpv.pv1xError = vtx1.xError()*micron;
	trpv.pv1yError = vtx1.yError()*micron;
	trpv.pv1zError = vtx1.zError()*micron;
	
	trpv.pv2x = vtx2.x()*micron;
	trpv.pv2y = vtx2.y()*micron;
	trpv.pv2z = vtx2.z()*micron;
	trpv.pv2xError = vtx2.xError()*micron;
	trpv.pv2yError = vtx2.yError()*micron;
	trpv.pv2zError = vtx2.zError()*micron;
	
	treePV->Fill();
     }   
}

void Residuals::beginJob()
{
}

void Residuals::endJob() 
{
}

bool Residuals::trackSelection(const reco::Track& track) const 
{   
   using namespace reco;
   
   if( track.pt() < tkMinPt ) return false;
   if( track.hitPattern().trackerLayersWithMeasurement() < tkMinXLayers ) return false;
   if( track.hitPattern().numberOfLostHits(HitPattern::MISSING_OUTER_HITS) > tkMaxMissedOuterLayers ) return false;
   if( track.hitPattern().numberOfLostHits(HitPattern::MISSING_INNER_HITS) > tkMaxMissedInnerLayers ) return false;   
   if( ! track.quality(reco::TrackBase::highPurity) ) return false;
//   if( ! (track.hitPattern().hasValidHitInPixelLayer(PixelSubdetector::PixelBarrel, 1) ||
//	  track.hitPattern().hasValidHitInPixelLayer(PixelSubdetector::PixelEndcap, 1)) ) return false;
   
   return true;
}

bool Residuals::vertexSelection(const reco::Vertex& vertex) const
{
   if( vertex.tracksSize()>vtxTracksSizeMax || vertex.tracksSize()<vtxTracksSizeMin ) return false;
   if( vertex.xError() < vtxErrorXMin || vertex.xError() > vtxErrorXMax ) return false;
   if( vertex.yError() < vtxErrorYMin || vertex.yError() > vtxErrorYMax ) return false;
   if( vertex.zError() < vtxErrorZMin || vertex.zError() > vtxErrorZMax ) return false;
   
   return true;
}

//define this as a plug-in
DEFINE_FWK_MODULE(Residuals);
