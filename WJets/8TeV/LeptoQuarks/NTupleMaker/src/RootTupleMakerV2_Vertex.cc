#include "Leptoquarks/RootTupleMakerV2/interface/RootTupleMakerV2_Vertex.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"


RootTupleMakerV2_Vertex::RootTupleMakerV2_Vertex(const edm::ParameterSet& iConfig) :
    inputTag(iConfig.getParameter<edm::InputTag>("InputTag")),
    prefix  (iConfig.getParameter<std::string>  ("Prefix")),
    suffix  (iConfig.getParameter<std::string>  ("Suffix"))
{
  produces <std::vector<float> > ( prefix + "X" + suffix );
  produces <std::vector<float> > ( prefix + "Y" + suffix );
  produces <std::vector<float> > ( prefix + "Z" + suffix );
  produces <std::vector<float> > ( prefix + "XErr" + suffix );
  produces <std::vector<float> > ( prefix + "YErr" + suffix );
  produces <std::vector<float> > ( prefix + "ZErr" + suffix );
  produces <std::vector<float> > ( prefix + "Rho" + suffix );
  produces <std::vector<float> > ( prefix + "Chi2" + suffix );
  produces <std::vector<float> > ( prefix + "NDF" + suffix );
  produces <std::vector<int> >    ( prefix + "NTracks" + suffix );
  produces <std::vector<int> >    ( prefix + "NTracksW05" + suffix );
  produces <std::vector<bool> >   ( prefix + "IsFake" + suffix );
}

void RootTupleMakerV2_Vertex::
produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {

  std::auto_ptr<std::vector<float> >  x  ( new std::vector<float>()  );
  std::auto_ptr<std::vector<float> >  y  ( new std::vector<float>()  );
  std::auto_ptr<std::vector<float> >  z  ( new std::vector<float>()  );
  std::auto_ptr<std::vector<float> >  xErr  ( new std::vector<float>()  );
  std::auto_ptr<std::vector<float> >  yErr  ( new std::vector<float>()  );
  std::auto_ptr<std::vector<float> >  zErr  ( new std::vector<float>()  );
  std::auto_ptr<std::vector<float> >  rho  ( new std::vector<float>()  );
  std::auto_ptr<std::vector<float> >  chi2  ( new std::vector<float>()  );
  std::auto_ptr<std::vector<float> >  ndf  ( new std::vector<float>()  );
  std::auto_ptr<std::vector<int> >     ntracks  ( new std::vector<int>()  );
  std::auto_ptr<std::vector<int> >     ntracksw05  ( new std::vector<int>()  );
  std::auto_ptr<std::vector<bool> >    isfake  ( new std::vector<bool>()  );

  //-----------------------------------------------------------------
  edm::Handle<reco::VertexCollection> primaryVertices;
  iEvent.getByLabel(inputTag,primaryVertices);

  if(primaryVertices.isValid()) {
    edm::LogInfo("RootTupleMakerV2_VertexInfo") << "Total # Primary Vertices: " << primaryVertices->size();

    for( reco::VertexCollection::const_iterator it=primaryVertices->begin() ; it!=primaryVertices->end() ; ++it ) {
      x->push_back( it->x() );
      y->push_back( it->y() );
      z->push_back( it->z() );
      xErr->push_back( it->xError() );
      yErr->push_back( it->yError() );
      zErr->push_back( it->zError() );
      rho->push_back( it->position().rho() );
      chi2->push_back( it->chi2() );
      ndf->push_back( it->ndof() );
      ntracks->push_back( int(it->tracksSize()) );
      ntracksw05->push_back( it->nTracks(0.5) ); // number of tracks in the vertex with weight above 0.5
      isfake->push_back( it->isFake() );
    }
  } else {
    edm::LogError("RootTupleMakerV2_VertexError") << "Error! Can't get the product " << inputTag;
  }

  //-----------------------------------------------------------------
  // put vectors in the event
  iEvent.put( x, prefix + "X" + suffix );
  iEvent.put( y, prefix + "Y" + suffix );
  iEvent.put( z, prefix + "Z" + suffix );
  iEvent.put( xErr, prefix + "XErr" + suffix );
  iEvent.put( yErr, prefix + "YErr" + suffix );
  iEvent.put( zErr, prefix + "ZErr" + suffix );
  iEvent.put( rho, prefix + "Rho" + suffix );
  iEvent.put( chi2, prefix + "Chi2" + suffix );
  iEvent.put( ndf, prefix + "NDF" + suffix );
  iEvent.put( ntracks, prefix + "NTracks" + suffix );
  iEvent.put( ntracksw05, prefix + "NTracksW05" + suffix );
  iEvent.put( isfake, prefix + "IsFake" + suffix );
}
