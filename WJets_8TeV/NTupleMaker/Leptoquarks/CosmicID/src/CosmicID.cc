// -*- C++ -*-
//
// Package:    CosmicID
// Class:      CosmicID
// 
/**\class CosmicID CosmicID.cc Leptoquarks/CosmicID/src/CosmicID.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Francesco Santanastasio,8 R-019,+41227675765,
//         Created:  Fri May 13 14:49:14 CEST 2011
// $Id: CosmicID.cc,v 1.1 2011/05/13 14:02:25 santanas Exp $
//
//


// system include files
#include <memory>
#include <iostream>
#include <string>
#include <vector>


// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/Common/interface/ValueMap.h"

#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonCosmicCompatibility.h"

//
// class declaration
//

class CosmicID : public edm::EDProducer {
   public:
      explicit CosmicID(const edm::ParameterSet&);
      ~CosmicID();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

   private:
      virtual void beginJob() ;
      virtual void produce(edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;
      
      virtual void beginRun(edm::Run&, edm::EventSetup const&);
      virtual void endRun(edm::Run&, edm::EventSetup const&);
      virtual void beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);
      virtual void endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);

  // ----------member data ---------------------------
  edm::InputTag src_;
  std::string result_;

};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
CosmicID::CosmicID(const edm::ParameterSet& iConfig)
{
  //register your products
  /* Examples
     produces<ExampleData2>();
     
     //if do put with a label
     produces<ExampleData2>("label");
     
     //if you want to put into the Run
     produces<ExampleData2,InRun>();
  */
  produces<edm::ValueMap<float> >().setBranchAlias("CosmicDiscriminators");
  
  //now do what ever other initialization is needed
  src_= iConfig.getParameter<edm::InputTag>("src");
  result_ = iConfig.getParameter<std::string>("result");
}


CosmicID::~CosmicID()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to produce the data  ------------
void
CosmicID::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  using namespace edm;
  using namespace reco;

  Handle<edm::ValueMap<reco::MuonCosmicCompatibility> > CosmicMap;
  iEvent.getByLabel( src_, CosmicMap );
  edm::Handle<reco::MuonCollection> muons;
  iEvent.getByLabel("muons",muons);
  std::vector<float> values;
  values.reserve(muons->size());

  unsigned int muonIdx = 0;
  for(reco::MuonCollection::const_iterator muon = muons->begin();
      muon != muons->end(); ++muon) {
    reco::MuonRef muonRef(muons, muonIdx);
    reco::MuonCosmicCompatibility muonCosmicCompatibility = (*CosmicMap)[muonRef];

    if(result_ == "cosmicCompatibility") values.push_back(muonCosmicCompatibility.cosmicCompatibility);
    if(result_ == "timeCompatibility") values.push_back(muonCosmicCompatibility.timeCompatibility);
    if(result_ == "backToBackCompatibility") values.push_back(muonCosmicCompatibility.backToBackCompatibility);
    if(result_ == "overlapCompatibility") values.push_back(muonCosmicCompatibility.overlapCompatibility);
    ++muonIdx;
  }

  std::auto_ptr<edm::ValueMap<float> > out(new edm::ValueMap<float>());
  edm::ValueMap<float>::Filler filler(*out);
  filler.insert(muons, values.begin(), values.end());
  filler.fill();

  // put value map into event
  iEvent.put(out);

  /* This is an event example
  //Read 'ExampleData' from the Event
  Handle<ExampleData> pIn;
  iEvent.getByLabel("example",pIn);

  //Use the ExampleData to create an ExampleData2 which 
  // is put into the Event
  std::auto_ptr<ExampleData2> pOut(new ExampleData2(*pIn));
  iEvent.put(pOut);
  */

  /* this is an EventSetup example
  //Read SetupData from the SetupRecord in the EventSetup
  ESHandle<SetupData> pSetup;
  iSetup.get<SetupRecord>().get(pSetup);
  */
}

// ------------ method called once each job just before starting event loop  ------------
void 
CosmicID::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
CosmicID::endJob() {
}

// ------------ method called when starting to processes a run  ------------
void 
CosmicID::beginRun(edm::Run&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void 
CosmicID::endRun(edm::Run&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void 
CosmicID::beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
CosmicID::endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
CosmicID::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(CosmicID);
