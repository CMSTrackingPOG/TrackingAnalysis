// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/Common/interface/TriggerResults.h"

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

#include <TrackingTools/TrajectoryState/interface/PerigeeConversions.h>
#include <TrackingTools/TrajectoryState/interface/TrajectoryStateClosestToPoint.h>
#include <TrackingTools/PatternTools/interface/TSCPBuilderNoMaterial.h>
#include "SimDataFormats/Associations/interface/TrackToTrackingParticleAssociator.h"
#include "SimTracker/TrackerHitAssociation/interface/TrackerHitAssociator.h"
#include "SimTracker/Records/interface/TrackAssociatorRecord.h"

#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingVertex.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingVertexContainer.h"
#include "SimTracker/TrackAssociation/plugins/ParametersDefinerForTPESProducer.h"

#include "SimDataFormats/Track/interface/SimTrackContainer.h"

#include "IpResoStudies/EDAnalyzers/interface/VertexReProducer.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include <CommonTools/UtilAlgos/interface/TFileService.h>

#include "TH1F.h"
#include "TTree.h"

class TrackAssociatorByHits;

struct treeReso
{
   double pt;
   double p;
   double eta;
   double phi;
   int nXLayers;
   int nMissedOut;
   int nMissedIn;
   int hasPXL;
   int type;
   double dxyReso;
   double dzReso;
};

struct treeResp
{
   double pt;
   double p;
   double eta;
   double phi;
   int nXLayers;
   int nMissedOut;
   int nMissedIn;
   int hasPXL;
   int type;
   double dxyResp;
   double dzResp;
};
  
class VertexResponsesAndTrueResolutions : public edm::EDAnalyzer 
{   
 public:
   
   explicit VertexResponsesAndTrueResolutions(const edm::ParameterSet& pset);
   ~VertexResponsesAndTrueResolutions();
    
 private:
   
   virtual void beginJob();
   virtual void analyze(const edm::Event&, const edm::EventSetup&);
   virtual void endJob();
    
   bool trackSelection(const reco::Track& track) const;
   bool vertexSelection(const reco::Vertex& vertex) const;
  
   edm::EDGetTokenT<reco::BeamSpot> theBeamspotToken_;
   edm::EDGetTokenT<reco::TrackCollection> theTracksToken_;
   edm::EDGetTokenT<edm::View<reco::Track> > theTracksViewToken_;
   edm::EDGetTokenT<reco::VertexCollection> thePVToken_;
   
   edm::EDGetTokenT<edm::SimTrackContainer> theSimTrackToken_;
   
   edm::EDGetTokenT<reco::TrackToTrackingParticleAssociator> theAssociatorToken_;
   edm::EDGetTokenT<TrackingParticleCollection> theTPCollectionHToken_;
   edm::EDGetTokenT<TrackingVertexCollection> theTVCollectionHToken_;
   
   edm::EDGetTokenT<edm::TriggerResults> triggerBitsToken_;
   
   edm::InputTag trackLabel;
   edm::InputTag vertexLabel;

   double tkMinPt;
   int tkMinXLayers,tkMaxMissedOuterLayers,tkMaxMissedInnerLayers;

   unsigned int vtxTracksSizeMin;  
   unsigned int vtxTracksSizeMax;  
   double vtxErrorXMin,vtxErrorXMax;
   double vtxErrorYMin,vtxErrorYMax;
   double vtxErrorZMin,vtxErrorZMax;

   TH1I *h_trackTypes;
   TTree *tReso;
   TTree *tResp;
   treeReso reso;
   treeResp resp;
};

VertexResponsesAndTrueResolutions::VertexResponsesAndTrueResolutions(const edm::ParameterSet& pset)
{
   trackLabel  = pset.getParameter<edm::InputTag>("TrackLabel");    
   vertexLabel = pset.getParameter<edm::InputTag>("VertexLabel");    
   
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

   edm::InputTag BeamspotTag_ = edm::InputTag("offlineBeamSpot");
   theBeamspotToken_ = consumes<reco::BeamSpot>(BeamspotTag_);
   
   edm::InputTag TrackCollectionTag_ = pset.getParameter<edm::InputTag>("TrackLabel");
   theTracksToken_= consumes<reco::TrackCollection>(TrackCollectionTag_);
   theTracksViewToken_= consumes<edm::View<reco::Track> >(TrackCollectionTag_);
   
   edm::InputTag VertexCollectionTag_ = pset.getParameter<edm::InputTag>("VertexLabel");
   thePVToken_ = consumes<reco::VertexCollection>(VertexCollectionTag_);
   
   edm::InputTag AssociatorTag_ = edm::InputTag("trackAssociatorByHits");
//   edm::InputTag AssociatorTag_ = edm::InputTag("quickTrackAssociatorByHits");
   theAssociatorToken_ = consumes<reco::TrackToTrackingParticleAssociator>(AssociatorTag_);

   edm::InputTag TPCollectionHTag_ = edm::InputTag("mix", "MergedTrackTruth");
   theTPCollectionHToken_ = consumes<TrackingParticleCollection>(TPCollectionHTag_);

   edm::InputTag TVCollectionHTag_ = edm::InputTag("mix", "MergedTrackTruth");
   theTVCollectionHToken_ = consumes<TrackingVertexCollection>(TVCollectionHTag_);

   edm::InputTag TriggerResultsTag_ = edm::InputTag("TriggerResults", "", "HLT");
   triggerBitsToken_ = consumes<edm::TriggerResults>(TriggerResultsTag_);
   
//   edm::InputTag SimTrackTag_ = edm::InputTag("g4SimHits");
//   theSimTrackToken_ = consumes<edm::SimTrackContainer>(SimTrackTag_);
   
   edm::Service<TFileService> fs;
   tReso = fs->make<TTree>("tReso", "resolutions");
   tResp = fs->make<TTree>("tResp", "vertex smearing");
   tReso->Branch("reso",&reso.pt,"pt/D:p/D:eta/D:phi/D:nXLayers/I:nMissedOut/I:nMissedIn/I:hasPXL/I:type/I:dxyReso/D:dzReso/D");
   tResp->Branch("resp",&resp.pt,"pt/D:p/D:eta/D:phi/D:nXLayers/I:nMissedOut/I:nMissedIn/I:hasPXL/I:type/I:dxyResp/D:dzResp/D");
   h_trackTypes = fs->make<TH1I>("trackTypes", "track types", 7, 0, 7 );
}

VertexResponsesAndTrueResolutions::~VertexResponsesAndTrueResolutions()
{
}

void VertexResponsesAndTrueResolutions::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;
   using namespace reco;
   using namespace std;
   
   edm::Handle<reco::TrackToTrackingParticleAssociator> theAssociator;
   iEvent.getByToken(theAssociatorToken_, theAssociator);
   
   Handle<TrackingParticleCollection>  TPCollectionH;
   iEvent.getByToken(theTPCollectionHToken_, TPCollectionH);
   
   Handle<TrackingVertexCollection> TVCollectionH;
   iEvent.getByToken(theTVCollectionHToken_, TVCollectionH);

   ESHandle<ParametersDefinerForTP> parametersDefinerTP; 
   iSetup.get<TrackAssociatorRecord>().get("LhcParametersDefinerForTP",parametersDefinerTP); 

   ESHandle<MagneticField> theMF;   
   iSetup.get<IdealMagneticFieldRecord>().get(theMF);

   Handle<VertexCollection> vtxH;
   iEvent.getByToken(thePVToken_, vtxH);

//   Handle<edm::SimTrackContainer> theSimTrack;
//   iEvent.getByToken(theSimTrackToken_, theSimTrack);
   
   VertexReProducer revertex(vtxH, iEvent);
   
//   Handle<TrackCollection> pvtracks;   
//   iEvent.getByLabel(revertex.inputTracks(),   pvtracks);

//   Handle<BeamSpot> pvbeamspot; 
//   iEvent.getByLabel(revertex.inputBeamSpot(), pvbeamspot);

   Handle<BeamSpot> pvbeamspot;
   iEvent.getByToken(theBeamspotToken_, pvbeamspot);
   
//   Handle<TrackCollection> tracks;  
//   iEvent.getByLabel(trackLabel, tracks);

//   if( tracks.id() != pvtracks.id() )
//     cout << "WARNING: the tracks originally used for PV are not the same used in this analyzer." 
//     << "Is this really what you want?" << endl;
  
   Handle<TrackCollection> tracks;
   iEvent.getByToken(theTracksToken_, tracks);
   
   Handle<View<Track> > trackViews;   
   iEvent.getByToken(theTracksViewToken_, trackViews);

//   edm::Handle<edm::TriggerResults> triggerBits;
//   iEvent.getByToken(triggerBitsToken_,triggerBits);
//   const edm::TriggerNames &names = iEvent.triggerNames(*triggerBits);
//   for(unsigned int i=0;i<names.size();i++)
//     {
//	std::cout << i << " " << names.triggerName(i) << std::endl;
//     }   
   
   //if (pvbeamspot.id() != theBeamSpot.id()) 
   //edm::LogWarning("Inconsistency") << "The BeamSpot used for PV reco is not the same used in this analyzer.";

   reco::RecoToSimCollection recSimColl = theAssociator->associateRecoToSim(trackViews, TPCollectionH);

   unsigned int counter;
   TrackCollection::const_iterator itk;
   for(itk = tracks->begin(), counter = 0; itk != tracks->end(); ++itk, ++counter)
     {
	RefToBase<Track> refTk(trackViews, counter); 
	h_trackTypes->Fill(0.); //fill bin for all tracks

	// --- track selection ---
	if( ! trackSelection(*itk) ) continue;
	h_trackTypes->Fill(1.); //fill bin for all tracks passing the track selection
	// ---
     
	// reco-sim association
	if( recSimColl.find(refTk) == recSimColl.end() ) continue;
	h_trackTypes->Fill(2.); //fill bin for all tracks, passing the track selection, which are matched to TP

	std::vector<std::pair<TrackingParticleRef, double> > tp = recSimColl[refTk];
	TrackingParticleRef tpr = tp.begin()->first;

	// ========== evaluating True MC-derived resolutions ================
	auto momentumTP_bs = parametersDefinerTP->momentum(iEvent,iSetup,tpr);
	auto vertexTP_bs = parametersDefinerTP->vertex(iEvent,iSetup,tpr);
	
	//double qoverpSim = tpr->charge()/sqrt(momentumTP.x()*momentumTP.x()+
	//     momentumTP.y()*momentumTP.y()+momentumTP.z()*momentumTP.z());
	//double lambdaSim = M_PI/2-momentumTP.theta();
	//double phiSim    = momentumTP.phi();
	double dxySim = (-vertexTP_bs.x()*sin(momentumTP_bs.phi())+vertexTP_bs.y()*cos(momentumTP_bs.phi()));
	double dzSim = vertexTP_bs.z() - (vertexTP_bs.x()*momentumTP_bs.x()+vertexTP_bs.y()*momentumTP_bs.y())/
	  sqrt(momentumTP_bs.perp2()) * momentumTP_bs.z()/sqrt(momentumTP_bs.perp2());

	double dxyRec = itk->dxy(pvbeamspot->position());
	double dzRec = itk->dz(pvbeamspot->position());
	
	double dxyRes = dxyRec - dxySim;
	double dzRes = dzRec - dzSim;
     
	reso.pt  = tpr->pt();
	reso.p = tpr->p();
	reso.eta = tpr->eta();
	reso.phi = tpr->phi();
	reso.nXLayers = itk->hitPattern().trackerLayersWithMeasurement();
	reso.nMissedOut = itk->hitPattern().numberOfLostHits(HitPattern::MISSING_OUTER_HITS);
	reso.nMissedIn = itk->hitPattern().numberOfLostHits(HitPattern::MISSING_INNER_HITS);
	reso.hasPXL =  (itk->hitPattern().hasValidHitInPixelLayer(PixelSubdetector::SubDetector::PixelBarrel, 1) || 
			itk->hitPattern().hasValidHitInPixelLayer(PixelSubdetector::SubDetector::PixelEndcap, 1));
	reso.dxyReso = dxyRes*10000.;
	reso.dzReso  = dzRes*10000.;
	// ============================ DONE ==============================

	// ===== check if the associated TP is really a particle with genuine zero IP ====
	int entryType;
	TrackingVertexRef tv(TVCollectionH,0);
	if( tpr->parentVertex().get() != tv.get() )
	  {
	     //cout << "matched tp is not a prompt particle. it has status,charge: "
	     //   << tpr->status() << " , " << tpr->charge() << endl;
	     if( tpr->genParticles().size() )
	       {
		  //cout << "It is pythia particle from decay of long living particle" << endl;
		  entryType = 5;
	       }
	     else 
	       {
		  //cout << "no genP mother. It is geant particle" << endl;
		  entryType = 6;
	       }       
	  }
	else
	  {
	     //fill bin for all tracks that pass the track selection, are matched to TP AND are genuinely prompt
	     entryType = 3;
	  } 
	h_trackTypes->Fill(entryType);
	reso.type = entryType;
	resp.type = entryType;
	// ============================

	tReso->Fill();
            
	// ========== evaluating vertex "smearing" impact on IP ================
	// ------ check if the vertex is good enough -------
	if( vtxH->size() == 0 ) return;
	if( ! vertexSelection(vtxH->front()) ) return;
	// -------------------------------------------------

//	 cout << "original vtx x,y,z: " 
//	 << vtxH->front().position().x() << " , "
//	 << vtxH->front().position().y() << " , "
//	 << vtxH->front().position().z() << endl;

	auto momentumTP = tpr->momentum(); 
	auto vertexTP = tpr->vertex(); 

	const FreeTrajectoryState ftsAtProduction(GlobalPoint(vertexTP.x(),vertexTP.y(),vertexTP.z()),
						  GlobalVector(momentumTP.x(),momentumTP.y(),momentumTP.z()),
						  TrackCharge(tpr->charge()),
						  theMF.product());

//	 cout << "tp bs vertex: "
//	 << vertexTP_bs.x() << " , "
//	 << vertexTP_bs.y() << " , "
//	 << vertexTP_bs.z() << endl;
	
//	 cout << "tp vertex: " 
//	 << vertexTP.x() << " , "
//	 << vertexTP.y() << " , "
//	 << vertexTP.z() << endl;
     
	GlobalPoint vtxGPosition(vtxH->front().position().x(),
				 vtxH->front().position().y(),
				 vtxH->front().position().z());
	
	TrajectoryStateClosestToPoint closestToPointVTX = TSCPBuilderNoMaterial()(ftsAtProduction,vtxGPosition);

	if( !closestToPointVTX.isValid() ) continue;

	GlobalPoint closestStatePointVTX = closestToPointVTX.position();
	GlobalVector closestStateVectorVTX = closestToPointVTX.momentum();

	/*
	 cout << "tp closestPoint: " 
	 << closestStatePointVTX.x() << " , "
	 << closestStatePointVTX.y() << " , "
	 << closestStatePointVTX.z() << endl;
	 */

	GlobalPoint v1(closestStatePointVTX.x()-vtxGPosition.x(),
		       closestStatePointVTX.y()-vtxGPosition.y(),
		       closestStatePointVTX.z()-vtxGPosition.z());
     
	//double qoverpSim = tsAtClosestApproach.trackStateAtPCA().charge()/closestStateVector.mag();
	//double lambdaSim = M_PI/2-closestStateVector.theta();
	//double phiSim    = closestStateVector.phi();
	//cout << "dxySim: " << dxySim << endl;
	double dxyResp = (-v1.x()*sin(closestStateVectorVTX.phi())+v1.y()*cos(closestStateVectorVTX.phi()));
	double dzResp = v1.z() - (v1.x()*closestStateVectorVTX.x()+v1.y()*closestStateVectorVTX.y())/
	  closestStateVectorVTX.perp() * closestStateVectorVTX.z()/closestStateVectorVTX.perp();

	resp.pt = tpr->pt();     
	resp.p = tpr->p();
	resp.eta = tpr->eta();
	resp.phi = tpr->phi();
	resp.nXLayers = itk->hitPattern().trackerLayersWithMeasurement();
	resp.nMissedOut = itk->hitPattern().numberOfLostHits(HitPattern::MISSING_OUTER_HITS);
	resp.nMissedIn = itk->hitPattern().numberOfLostHits(HitPattern::MISSING_INNER_HITS);
	resp.hasPXL = (itk->hitPattern().hasValidHitInPixelLayer(PixelSubdetector::SubDetector::PixelBarrel, 1) || 
		       itk->hitPattern().hasValidHitInPixelLayer(PixelSubdetector::SubDetector::PixelEndcap, 1));
	resp.dxyResp = dxyResp*10000.;
	resp.dzResp = dzResp*10000.;
	tResp->Fill();     
	// ============================ DONE ==============================
     }
}

void VertexResponsesAndTrueResolutions::beginJob()
{
}

void VertexResponsesAndTrueResolutions::endJob() 
{
}

bool VertexResponsesAndTrueResolutions::trackSelection(const reco::Track& track) const 
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

bool VertexResponsesAndTrueResolutions::vertexSelection(const reco::Vertex& vertex) const
{
   if( vertex.tracksSize() > vtxTracksSizeMax || vertex.tracksSize() < vtxTracksSizeMin ) return false;
   if( vertex.xError() < vtxErrorXMin || vertex.xError() > vtxErrorXMax ) return false;
   if( vertex.yError() < vtxErrorYMin || vertex.yError() > vtxErrorYMax ) return false;
   if( vertex.zError() < vtxErrorZMin || vertex.zError() > vtxErrorZMax ) return false;
   
   return true;
}

DEFINE_FWK_MODULE(VertexResponsesAndTrueResolutions);
