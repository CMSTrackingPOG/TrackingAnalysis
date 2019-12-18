#include "TrackingAnalysis/EDAnalyzers/interface/VertexReProducer.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "FWCore/ParameterSet/interface/Registry.h"
#include "FWCore/Common/interface/Provenance.h"
#include "DataFormats/Provenance/interface/ProductProvenance.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

VertexReProducer::VertexReProducer(const edm::Handle<reco::VertexCollection> &handle, const edm::Event &iEvent) 
{
   const edm::Provenance *prov = handle.provenance();
   if( prov == nullptr ) throw cms::Exception("CorruptData") << "Vertex handle doesn't have provenance.";
   edm::ParameterSet psetFromProvenance = edm::parameterSet(*prov);
   
   bool is_primary_available = false;
   const edm::Provenance *parent_prov = prov;
   if( edm::moduleName(*prov) != "PrimaryVertexProducer" )
     {
	std::vector<edm::BranchID> parents = prov->productProvenance()->parentage().parents();
	for( std::vector<edm::BranchID>::const_iterator it = parents.begin(), ed = parents.end(); it != ed; ++it )
	  {
	     edm::Provenance parprov = iEvent.getProvenance(*it);
	     if( parprov.friendlyClassName() == "recoVertexs" )
	       {  
		  // for AOD actually this the parent we should look for
		  parent_prov = &parprov;
		  psetFromProvenance = edm::parameterSet(parprov);
		  is_primary_available = true;
		  break;
	       }
	  }
     } else is_primary_available = true;
   
   if( is_primary_available )
     prov = parent_prov;
   else
     throw cms::Exception("Configuration") << "Vertices to re-produce don't come from a PrimaryVertexProducer \n";
   
   configure(psetFromProvenance);

   // Now we also dig out the ProcessName used for the reco::Tracks and reco::Vertices
/*kskovpen   std::vector<edm::BranchID> parents = prov->productProvenance()->parentage().parents();
   bool foundTracks = false;
   bool foundBeamSpot = false;
   for( std::vector<edm::BranchID>::const_iterator it = parents.begin(), ed = parents.end(); it != ed; ++it )
     {
	edm::Provenance parprov = iEvent.getProvenance(*it);
	if( parprov.friendlyClassName() == "recoTracks" )
	  {
	     tracksTag_ = edm::InputTag(parprov.moduleLabel(), parprov.productInstanceName(), parprov.processName());
	     foundTracks = true;
	     if( parprov.moduleLabel() != "generalTracks" )
	       foundTracks = false;  // this is necessary since we are asking for that in onia2mumu
	  } 
	else if( parprov.friendlyClassName() == "recoBeamSpot" )
	  {
	     beamSpotTag_ = edm::InputTag(parprov.moduleLabel(), parprov.productInstanceName(), parprov.processName());
	     foundBeamSpot = true;
	     if( parprov.moduleLabel() != "offlineBeamSpot" )
	       foundBeamSpot = false;  // this is necessary since we are asking for that in onia2mumu
	  }
     }
   if( !foundTracks || !foundBeamSpot )
     {
	//edm::LogWarning("VertexReProducer_MissingParentage") <<
	throw cms::Exception("Configuration")
	  << "Can't find correct parentage info for vertex collection inputs: " << (foundTracks ? "" : "generalTracks ")
	    << (foundBeamSpot ? "" : "offlineBeamSpot") << "\n";
     }*/
}

void VertexReProducer::configure(const edm::ParameterSet &iConfig) 
{
   config_ = iConfig;
   tracksTag_ = iConfig.getParameter<edm::InputTag>("TrackLabel");
   beamSpotTag_ = iConfig.getParameter<edm::InputTag>("beamSpotLabel");
   algo_.reset(new PrimaryVertexProducerAlgorithm(iConfig));

//   std::cout<<"configuration"<< iConfig.dump();
}

std::vector<TransientVertex> VertexReProducer::makeVertices(const reco::TrackCollection &tracks,
							    const reco::BeamSpot &bs,
							    bool useBeamSpot,
							    const edm::EventSetup &iSetup) const
{
   edm::ESHandle<TransientTrackBuilder> theB;
   iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder", theB);
   
   std::vector<reco::TransientTrack> t_tks;
   t_tks.reserve(tracks.size());
   
   for( reco::TrackCollection::const_iterator it = tracks.begin(), ed = tracks.end(); it != ed; ++it )
     {
	t_tks.push_back((*theB).build(*it));
	if( useBeamSpot ) t_tks.back().setBeamSpot(bs);
     }
   
   return algo_->vertices(t_tks, bs, "");
}

std::vector<TransientVertex> VertexReProducer::makeVertices(const std::vector<reco::TrackBaseRef> &tracks,
							    const reco::BeamSpot &bs,
							    bool useBeamSpot,
							    const edm::EventSetup &iSetup) const
{
   edm::ESHandle<TransientTrackBuilder> theB;
   iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder", theB);
   
   std::vector<reco::TransientTrack> t_tks;
   t_tks.reserve(tracks.size());
   
   for( std::vector<reco::TrackBaseRef>::const_iterator it = tracks.begin(), ed = tracks.end(); it != ed; ++it )
     {
	reco::TrackBaseRef tk = *it;
	t_tks.push_back((*theB).build(*tk));
	if( useBeamSpot ) t_tks.back().setBeamSpot(bs);
     }
   
   return algo_->vertices(t_tks, bs, "");
}
