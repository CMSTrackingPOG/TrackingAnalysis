#include "TrackingAnalysis/EDAnalyzers/interface/VertexReProducer.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "FWCore/ParameterSet/interface/Registry.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/getProducerParameterSet.h"
#include "DataFormats/Provenance/interface/Provenance.h"
#include "DataFormats/Provenance/interface/ProductProvenance.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "DataFormats/Provenance/interface/BranchDescription.h"

VertexReProducer::VertexReProducer(const edm::Handle<reco::VertexCollection> &handle, const edm::Event &iEvent) 
{
   const edm::Provenance *prov = handle.provenance();
   if( prov == nullptr ) throw cms::Exception("CorruptData") << "Vertex handle doesn't have provenance.";

   if( prov->moduleName() == "PrimaryVertexProducer" )
     {
	const edm::ParameterSet *psetFromProvenance = getProducerParameterSet(*prov);

	configure(*psetFromProvenance);

	std::vector<edm::BranchID> parents = prov->productProvenance()->parentage().parents();
	
	bool foundTracks = false;
	bool foundBeamSpot = false;
	for( std::vector<edm::BranchID>::const_iterator it = parents.begin(), ed = parents.end(); it != ed; ++it )
	  {
	     edm::Provenance parprov = iEvent.getProvenance(*it);
	     if( parprov.friendlyClassName() == "recoTracks" )
	       {
		  tracksTag_ = edm::InputTag(parprov.moduleLabel(), parprov.productInstanceName(), parprov.processName());
		  foundTracks = true;
	       }
	     else if( parprov.friendlyClassName() == "recoBeamSpot" )
	       {
		  beamSpotTag_ = edm::InputTag(parprov.moduleLabel(), parprov.productInstanceName(), parprov.processName());
		  foundBeamSpot = true;
	       }	
	  }
	
	if( !foundTracks || !foundBeamSpot )
	  {
	     edm::LogWarning("VertexReProducer_MissingParentage") << 
	       "Can't find parentage info for vertex collection inputs: " << 
	       (foundTracks ? "" : "tracks ") << (foundBeamSpot ? "" : "beamSpot") << "\n";
	  }
     }   
   else throw cms::Exception("Configuration") << "Vertices to re-produce don't come from a PrimaryVertexProducer \n";
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
							    const edm::EventSetup &iSetup) const
{
   edm::ESHandle<TransientTrackBuilder> theB;
   iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder", theB);
   
   std::vector<reco::TransientTrack> t_tks;
   t_tks.reserve(tracks.size());

   for( reco::TrackCollection::const_iterator it = tracks.begin(), ed = tracks.end(); it != ed; ++it )
     {
	reco::TransientTrack transientTrack = (*theB).build(*it);
	t_tks.push_back(transientTrack);
	t_tks.back().setBeamSpot(bs);
     }

   return algo_->vertices(t_tks, bs);
}

std::vector<TransientVertex> VertexReProducer::makeVertices(const std::vector<reco::TrackBaseRef> &tracks,
							    const reco::BeamSpot &bs,
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
	t_tks.back().setBeamSpot(bs);
     }
   
   return algo_->vertices(t_tks, bs);
}
