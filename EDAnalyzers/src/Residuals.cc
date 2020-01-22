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

#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"

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

#include "TrackingAnalysis/EDAnalyzers/interface/VertexReProducer.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include <CommonTools/UtilAlgos/interface/TFileService.h>

#include "TROOT.h"
#include "TH1F.h"
#include "TTree.h"
#include "TRandom3.h"
#include "Compression.h"

#include "TrackingAnalysis/EDAnalyzers/interface/Tree.h"

class Residuals : public edm::EDAnalyzer 
{  
   
 public:
   explicit Residuals(const edm::ParameterSet& pset);
   ~Residuals();
    
 private:
   virtual void beginRun(const edm::Run&, const edm::EventSetup&);
   virtual void analyze(const edm::Event&, const edm::EventSetup&);
   virtual void endRun() ;

   bool trackSelection(const reco::Track& track) const;
   bool vertexSelection(const reco::Vertex& vertex) const;
  
   // ----------member data ---------------------------
   edm::InputTag thePVLabel_;
   edm::InputTag thePVPrimaryLabel_;
   edm::InputTag theTrackLabel_;
   edm::InputTag theBeamspotLabel_;
   edm::InputTag theRhoLabel_;
   edm::InputTag theTriggerBitsLabel_;

   // --- track selection variables
   double tkMinPt;
   int tkMinXLayers,tkMaxMissedOuterLayers,tkMaxMissedInnerLayers;

   // --- vertex selection variables
   unsigned int vtxTracksSizeMin;  
   unsigned int vtxTracksSizeMax;  
   double vtxErrorXMin,vtxErrorXMax;
   double vtxErrorYMin,vtxErrorYMax;
   double vtxErrorZMin,vtxErrorZMax;

   int eventScale;
   
   TRandom3 *rnd;
   
   HLTConfigProvider hltConfig_;
   
   const edm::Service<TFileService> fs;
   ResTree* ftree;
   
   int ncount;
};

Residuals::Residuals(const edm::ParameterSet& pset)
{
   thePVLabel_  = pset.getParameter<edm::InputTag>("VertexLabel");
   thePVPrimaryLabel_  = pset.getParameter<edm::InputTag>("VertexPrimaryLabel");
   theTrackLabel_  = pset.getParameter<edm::InputTag>("TrackLabel");
   theBeamspotLabel_  = pset.getParameter<edm::InputTag>("BeamSpotLabel");
   theRhoLabel_  = pset.getParameter<edm::InputTag>("RhoLabel");
   theTriggerBitsLabel_  = pset.getParameter<edm::InputTag>("TriggerResultsLabel");
   
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

   eventScale = pset.getParameter<int>("EventScale");
   
   rnd = new TRandom3();

   TFile& f = fs->file();
   f.SetCompressionAlgorithm(ROOT::kZLIB);
   f.SetCompressionLevel(9);
   ftree = new ResTree(fs->make<TTree>("tree","tree"));   
   ftree->CreateBranches(32000);
   
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

   Handle<BeamSpot> pvbeamspot;
   iEvent.getByLabel(theBeamspotLabel_, pvbeamspot);
   
   Handle<TrackCollection> tracks;
   iEvent.getByLabel(theTrackLabel_,tracks);
   
//   ESHandle<TransientTrackBuilder> transTrackBuilder;
//   iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder", transTrackBuilder);
   
//   std::cout << "size of track collection: "<< tracks->size() << std::endl;

   Handle<VertexCollection> vtxH;
   iEvent.getByLabel(thePVLabel_, vtxH);

   Handle<VertexCollection> vtxP;
   iEvent.getByLabel(thePVPrimaryLabel_, vtxP);
   
   if( !vtxH.isValid() || !vtxP.isValid() ) return;

   if( vtxP->size() == 0 || vtxH->size() == 0 ) return;
   
   VertexReProducer revertex(vtxH, iEvent);
   
   TrackCollection tracksAll;
   tracksAll.assign(tracks->begin(), tracks->end());

   // refit primary vertices
   vector<TransientVertex> pvs0 = revertex.makeVertices(tracksAll, *pvbeamspot, iSetup);
   
   if( pvs0.empty() ) return;
   
   reco::Vertex vtx0 = reco::Vertex(pvs0.front());
   reco::Vertex vtx = vtxP->front();
   
   // refitted vertices are the same as the original ones
   if( fabs(vtx.x()-vtx0.x()) > 10E-10 || fabs(vtx.y()-vtx0.y()) > 10E-10 || fabs(vtx.z()-vtx0.z()) > 10E-10 ) return;
   
   if( !vertexSelection(vtx) ) return;

   ESHandle<MagneticField> theMF;
   iSetup.get<IdealMagneticFieldRecord>().get(theMF);

   //Handle<TrackCollection> pvtracks;
   //iEvent.getByLabel(revertex.inputTracks(),   pvtracks);
   //Handle<BeamSpot>        pvbeamspot;
   //iEvent.getByLabel(revertex.inputBeamSpot(), pvbeamspot);
   
   Handle<double> rhoPtr;
   iEvent.getByLabel(theRhoLabel_, rhoPtr);   
   
   ftree->ev_run = iEvent.id().run();
   ftree->ev_id = iEvent.id().event();
   ftree->ev_lumi = iEvent.id().luminosityBlock();
   ftree->ev_rho = *rhoPtr;
   
   Handle<TriggerResults> triggerBits;
   iEvent.getByLabel(theTriggerBitsLabel_, triggerBits);
   const TriggerNames &names = iEvent.triggerNames(*triggerBits);
   for( unsigned int i=0;i<names.size();i++ )
     {
	std::string trigName = names.triggerName(i);
//	std::cout << i << " " << trigName << std::endl;
	
	if( TString(trigName).Contains("HLT_Physics_v") ||
	    TString(trigName).Contains("HLT_PixelTracks_Multiplicity70_v") ||
	    TString(trigName).Contains("HLT_PixelTracks_Multiplicity80_v") ||
	    TString(trigName).Contains("HLT_PixelTracks_Multiplicity90_v") ||
	    TString(trigName).Contains("HLT_Random_v") ||
	    TString(trigName).Contains("HLT_ZeroBiasPixel_DoubleTrack_v") ||
	    TString(trigName).Contains("HLT_ZeroBias_v") )
	  {	    
	     std::pair<int,int> prescales =
	       hltConfig_.prescaleValues(iEvent, iSetup, trigName);
	
	     int L1ps = prescales.first;
	     int HLTps = prescales.second;
	     bool pass = (triggerBits->accept(i) ? true : false); 
	     	     
	     if( TString(trigName).Contains("HLT_ZeroBiasPixel_DoubleTrack_v") )
	       {
		  ftree->trig_ZeroBiasPixel_DoubleTrack_pass = pass;
		  ftree->trig_ZeroBiasPixel_DoubleTrack_L1ps = L1ps;
		  ftree->trig_ZeroBiasPixel_DoubleTrack_HLTps = HLTps;
	       }	     
	     else if( TString(trigName).Contains("HLT_ZeroBias_v") )
	       {
		  ftree->trig_ZeroBias_pass = pass;
		  ftree->trig_ZeroBias_L1ps = L1ps;
		  ftree->trig_ZeroBias_HLTps = HLTps;
	       }	     
	     
//	     std::cout << i << " " << trigName << " L1=" << L1prescale << " HLT=" << HLTprescale << std::endl;
	  }	
     }
   
   //if(tracks.id() != pvtracks.id())
   // cout << "WARNING: the tracks originally used for PV are not the same used in this analyzer."
   //	  << "Is this really what you want?" << endl;

   //if (pvbeamspot.id() != theBeamSpot.id()) 
   //  edm::LogWarning("Inconsistency") << "The BeamSpot used for PV reco is not the same used in this analyzer.";

   double micron = 10000;

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

   double trackSumPt = 0;
   for( std::vector<reco::TrackBaseRef>::const_iterator it = vtx.tracks_begin(); it != vtx.tracks_end(); it++ )
     trackSumPt += (*it)->pt();
   
   ftree->pv_IsValid = vtx.isValid();
   ftree->pv_IsFake = vtx.isFake();	
   ftree->pv_NTracks = vtx.tracksSize();
   ftree->pv_SumTrackPt = trackSumPt;	
   ftree->pv_chi2 = vtx.chi2();	
   ftree->pv_ndof = vtx.ndof();   
   ftree->pv_x = vtx.x()*micron;
   ftree->pv_y = vtx.y()*micron;
   ftree->pv_z = vtx.z()*micron;
   ftree->pv_xError = vtx.xError()*micron;
   ftree->pv_yError = vtx.yError()*micron;
   ftree->pv_zError = vtx.zError()*micron;

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
	
	ftree->pv_SumTrackPt_p1 = trackSumPt1;
	ftree->pv_SumTrackPt_p2 = trackSumPt2;
	
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
  
	// Refit the primary vertex
	vector<TransientVertex> pvs = revertex.makeVertices(newTkCollection, *pvbeamspot, iSetup);
	//cout << "vertices before,after: " << vtxH->size() << " , " << pvs.size() << endl;
	
	if( pvs.empty() ) continue;

	reco::Vertex newPV = reco::Vertex(pvs.front());
	Track::Point vtxPosition = Track::Point(newPV.position().x(),
						newPV.position().y(),
						newPV.position().z());
	
	// ---
	if( ! vertexSelection(newPV) ) continue;

	double pt = itk->pt();
	double p = itk->p();
	double eta = itk->eta();
	double phi = itk->phi();
	
	int nXLayers = itk->hitPattern().trackerLayersWithMeasurement();
	int nMissedOut = itk->trackerExpectedHitsOuter().numberOfLostHits();
	int nMissedIn = itk->trackerExpectedHitsInner().numberOfLostHits();
	int hasPXL = (itk->hitPattern().hasValidHitInFirstPixelBarrel() || itk->hitPattern().hasValidHitInFirstPixelEndcap());
	double quality = itk->qualityMask();

	double d0 = itk->dxy();
	double dz = itk->dz();
	
	double d0_pv = itk->dxy(vtxPosition);
	double dz_pv = itk->dz(vtxPosition);
	double d0NoRefit_pv = itk->dxy(vtxH->front().position());
	double dzNoRefit_pv = itk->dz(vtxH->front().position());

	double d0_bs = itk->dxy(pvbeamspot->position());
	double d0_bs_zpca = itk->dxy(*pvbeamspot);
	double dz_bs = itk->dz(pvbeamspot->position());
	
	double d0_bs_zpv = itk->dxy(pvbeamspot->position(vtx.z()));
	
	double d0Error = itk->d0Error();
	double dzError = itk->dzError();
	
	//cout << "d0:" << d0 << " dz:" << dz << std::endl;

	ftree->trk_pt.push_back( pt );
	ftree->trk_p.push_back( p );
	ftree->trk_eta.push_back( eta );
	ftree->trk_phi.push_back( phi );
	ftree->trk_nXLayers.push_back( nXLayers );
	ftree->trk_nMissedOut.push_back( nMissedOut );
	ftree->trk_nMissedIn.push_back( nMissedIn );
	ftree->trk_hasPXL.push_back( hasPXL );
	ftree->trk_quality.push_back( quality );
	ftree->trk_d0.push_back( d0*micron );
	ftree->trk_dz.push_back( dz*micron );
	ftree->trk_d0_pv.push_back( d0_pv*micron );
	ftree->trk_dz_pv.push_back( dz_pv*micron );
	ftree->trk_d0_bs.push_back( d0_bs*micron );
	ftree->trk_d0_bs_zpca.push_back( d0_bs_zpca*micron );
	ftree->trk_d0_bs_zpv.push_back( d0_bs_zpv*micron );
	ftree->trk_dz_bs.push_back( dz_bs*micron );
	ftree->trk_d0Err.push_back( d0Error*micron );
	ftree->trk_dzErr.push_back( dzError*micron );
	ftree->trk_d0_pv_NoRefit.push_back( d0NoRefit_pv*micron );
	ftree->trk_dz_pv_NoRefit.push_back( dzNoRefit_pv*micron );
     }

   ftree->tree->Fill();
}

void Residuals::beginRun(const edm::Run& iRun, const edm::EventSetup& iSetup)
{
   bool changed = true;
   if( !hltConfig_.init(iRun, iSetup, "HLT", changed) )
     std::cout << "Warning, didn't find HLTConfigProvider with label "
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

//define this as a plug-in
DEFINE_FWK_MODULE(Residuals);
