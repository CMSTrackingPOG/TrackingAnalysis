#include "TrackingAnalysis/EDAnalyzers/interface/VertexReProducer.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/ParameterSet/interface/Registry.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/getProducerParameterSet.h"
#include "DataFormats/Provenance/interface/Provenance.h"
#include "DataFormats/Provenance/interface/ProductProvenance.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "RecoVertex/AdaptiveVertexFit/interface/AdaptiveVertexFitter.h"

#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"

#include "DataFormats/Provenance/interface/BranchDescription.h"

VertexReProducer::VertexReProducer(const edm::ParameterSet &iConfig)
{
   beamSpotConfig_ = iConfig.getParameter<std::string>("BeamSpotConfig");
   
   std::string trackSelectionAlgorithm = iConfig.getParameter<edm::ParameterSet>("TkFilterParameters").getParameter<std::string>("algorithm");
   if( trackSelectionAlgorithm == "filter" )
     {
	theTrackFilter_ = new TrackFilterForPVFinding(iConfig.getParameter<edm::ParameterSet>("TkFilterParameters"));
     }
   theTrackClusterizer_ = new DAClusterizerInZ_vect(iConfig.getParameter<edm::ParameterSet>("TkClusParameters").getParameter<edm::ParameterSet>("TkDAClusParameters"));
   
   vertexSelector_ = new VertexCompatibleWithBeam(VertexDistanceXY(), iConfig.getParameter<edm::ParameterSet>("VxFitterParameters").getParameter<double>("maxDistanceToBeam"));
   
   minNdof_ = iConfig.getParameter<edm::ParameterSet>("VxFitterParameters").getParameter<double>("minNdof");
}

std::vector<TransientVertex> VertexReProducer::makeVertices(const reco::TrackCollection &tracks,
							    const reco::BeamSpot &bs,
							    const edm::EventSetup &iSetup) const
{
   bool validBS = true;
   VertexState beamVertexState(bs);
   if ( (beamVertexState.error().cxx() <= 0.) ||
	(beamVertexState.error().cyy() <= 0.) ||
	(beamVertexState.error().czz() <= 0.) ) 
     {
	validBS = false;
	edm::LogError("UnusableBeamSpot") << "Beamspot with invalid errors " << beamVertexState.error().matrix();
     }   
   
   edm::ESHandle<TransientTrackBuilder> theB;
   iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder", theB);
   
   std::vector<reco::TransientTrack> t_tks;
   t_tks.reserve(tracks.size());

   for( size_t i=0;i<tracks.size();i++ )
     {
	const reco::Track &pseudoTrack = tracks[i];
	
	reco::TransientTrack transTrack = (*theB).build(pseudoTrack);
	t_tks.push_back(transTrack);
	t_tks.back().setBeamSpot(bs);
     }
   
   std::vector<reco::TransientTrack> seltks = theTrackFilter_->select(t_tks);
   
   std::vector<std::vector<reco::TransientTrack> > clusters = theTrackClusterizer_->clusterize(seltks);
   
   AdaptiveVertexFitter fit;
   
   std::vector<TransientVertex> pvs;
   
   for( std::vector< std::vector<reco::TransientTrack> >::const_iterator iclus=clusters.begin();iclus!=clusters.end();iclus++ )
     {
	TransientVertex v;
	
	if( beamSpotConfig_ == "WithBS" && validBS && ((*iclus).size()>1) )
	  {	     	
	     v = fit.vertex(*iclus, bs);
	  }
	else if( beamSpotConfig_ != "WithBS" && ((*iclus).size()>1) )
	  {
	     v = fit.vertex(*iclus);
	  }	
	
	if( v.isValid() && (v.degreesOfFreedom() >= minNdof_) &&
	    (!validBS || (*(vertexSelector_))(v, beamVertexState)) ) pvs.push_back(v);
     }
   
   if( pvs.size() > 1 )
     {
	sort(pvs.begin(), pvs.end(), VertexHigherPtSquared());
     }   
   
   return pvs;
}

VertexReProducer::~VertexReProducer()
{
   if( theTrackFilter_ ) delete theTrackFilter_;
   if( theTrackClusterizer_ ) delete theTrackClusterizer_;
   if( vertexSelector_ ) delete vertexSelector_;
}
