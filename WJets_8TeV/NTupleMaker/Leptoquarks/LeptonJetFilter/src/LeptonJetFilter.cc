// -*- C++ -*-
//
// Package:    LeptonJetFilter
// Class:      LeptonJetFilter
//
/**\class LeptonJetFilter LeptonJetFilter.cc AnalysisCode/LeptonJetFilter/src/LeptonJetFilter.cc

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Jeff Temple
//         Created:  Mon Aug  3 13:02:30 CEST 2009
// $Id: LeptonJetFilter.cc,v 1.9 2011/08/31 20:49:23 santanas Exp $
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDFilter.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/EgammaReco/interface/SuperCluster.h"
#include "DataFormats/EgammaReco/interface/SuperClusterFwd.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Tau.h"
#include "RecoEgamma/EgammaTools/interface/HoECalculator.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/EgammaCandidates/interface/Photon.h"
#include "DataFormats/EgammaCandidates/interface/PhotonFwd.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"


#include "TH1.h"

using namespace edm;
using namespace std;
using namespace reco;

//
// class declaration
//

class LeptonJetFilter : public edm::EDFilter {
   public:
      explicit LeptonJetFilter(const edm::ParameterSet&);
      ~LeptonJetFilter();

   private:
      virtual void beginJob() ;
      virtual bool filter(edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;

      // ----------member data ---------------------------
      int photMin_, jetsMin_, muonsMin_, electronsMin_, tausMin_;
      int photMax_, jetsMax_, muonsMax_, electronsMax_, tausMax_;
      double  photET_, photEta_, photHoE_, jetPT_, jetEta_, muPT_, muEta_, elecPT_, elecEta_, tauPT_, tauEta_;
      bool useMuID_, useElecID_, usePhotID_, useTauID_;
      string muID_, elecID_, photID_, tauID_;
      bool debug_, allEventsPassFilter_;
      bool counteitherleptontype_;
      bool customfilterGenRecoMuJet_;
      edm::InputTag photLabel_, jetLabel_, tauLabel_, muLabel_, elecLabel_,genjetLabel_, genparticleLabel_;

      int TotalCount;
      int PassedCount;

      edm::Service<TFileService> histServ;
      TH1I* hCount;
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
LeptonJetFilter::LeptonJetFilter(const edm::ParameterSet& iConfig)
{
   //now do what ever initialization is needed

  // specify labels
  photLabel_  = iConfig.getUntrackedParameter<edm::InputTag>("photLabel");
  jetLabel_  = iConfig.getParameter<edm::InputTag>("jetLabel");
  tauLabel_  = iConfig.getParameter<edm::InputTag>("tauLabel");
  muLabel_   = iConfig.getParameter<edm::InputTag>("muLabel");
  elecLabel_ = iConfig.getParameter<edm::InputTag>("elecLabel");

  genjetLabel_  = iConfig.getParameter<edm::InputTag>("genjetLabel");
  genparticleLabel_   = iConfig.getParameter<edm::InputTag>("genparticleLabel");

  photMin_   = iConfig.getParameter<int>("photMin");
  photMax_   = iConfig.getParameter<int>("photMax");
  photET_    = iConfig.getParameter<double>("photET");
  photEta_   = iConfig.getParameter<double>("photEta");
  photHoE_   = iConfig.getParameter<double>("photHoE");
  usePhotID_ = iConfig.getParameter<bool>("usePhotID");
  photID_    = iConfig.getParameter<string>("photID");

  jetsMin_ = iConfig.getParameter<int>("jetsMin");
  jetsMax_ = iConfig.getParameter<int>("jetsMax");
  jetPT_   = iConfig.getParameter<double>("jetPT");
  jetEta_  = iConfig.getParameter<double>("jetEta");

  muonsMin_ = iConfig.getParameter<int>("muonsMin");
  muonsMax_ = iConfig.getParameter<int>("muonsMax");
  muPT_     = iConfig.getParameter<double>("muPT");
  muEta_    = iConfig.getParameter<double>("muEta");
  useMuID_  = iConfig.getParameter<bool>("useMuID");
  muID_     = iConfig.getParameter<string>("muID");

  electronsMin_ = iConfig.getParameter<int>("electronsMin");
  electronsMax_ = iConfig.getParameter<int>("electronsMax");
  elecPT_       = iConfig.getParameter<double>("elecPT");
  elecEta_      = iConfig.getParameter<double>("elecEta");
  useElecID_    = iConfig.getParameter<bool>("useElecID");
  elecID_       = iConfig.getParameter<string>("elecID");

  tausMin_  = iConfig.getParameter<int>("tausMin");
  tausMax_  = iConfig.getParameter<int>("tausMax");
  tauPT_    = iConfig.getParameter<double>("tauPT" );
  tauEta_   = iConfig.getParameter<double>("tauEta") ;
  useTauID_ = iConfig.getParameter<bool> ("useTauID");
  tauID_    = iConfig.getParameter<string>("tauID") ;

  allEventsPassFilter_ = iConfig.getUntrackedParameter<bool>("allEventsPassFilter");
  debug_   = iConfig.getUntrackedParameter<bool>("debug");
  counteitherleptontype_ = iConfig.getParameter<bool>("counteitherleptontype");
  customfilterGenRecoMuJet_ = iConfig.getParameter<bool>("customfilterGenRecoMuJet");

  TotalCount=0;
  PassedCount=0;
  TFileDirectory COUNT=histServ->mkdir("EventCount");
  hCount=COUNT.make<TH1I>("EventCounter","Event Counter",2,-0.5,1.5);
  hCount->GetXaxis()->SetBinLabel(1,"all events");
  hCount->GetXaxis()->SetBinLabel(2,"passed");
}


LeptonJetFilter::~LeptonJetFilter()
{

   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called on each new Event  ------------
bool
LeptonJetFilter::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  ++TotalCount;

  if( allEventsPassFilter_ )
    {
      ++PassedCount;
      if (debug_) cout <<"PASSED!"<<endl;
      return true;
    }

  // get photons
  Handle<PhotonCollection> photons;
  iEvent.getByLabel(photLabel_, photons);

  // get jets
  edm::Handle<edm::View<reco::Candidate> > jets;
  iEvent.getByLabel(jetLabel_, jets);

  // get leptons
  edm::Handle<edm::View<reco::Candidate> > electrons;
  iEvent.getByLabel(elecLabel_, electrons);

  edm::Handle<edm::View<reco::Candidate> > muons;
  iEvent.getByLabel(muLabel_, muons);

  edm::Handle<edm::View<reco::Candidate> > taus;
  iEvent.getByLabel(tauLabel_, taus);

  float MaxGenMuon = 0.0;
  float MaxRecMuon = 0.0;
  float MaxGenJet = 0.0;
  float MaxRecJet = 0.0;

  // Step 0: count photons
  int nphotons = 0;
  // std::cout<<" --------------- SKIM --------------------"<<std::endl;

  for (PhotonCollection::const_iterator it = photons->begin(); it != photons->end();++it)
    {
      if (debug_) cout << "Photons:" << endl;
      if (debug_) cout << "pT: " << it->pt() << " eta: " <<  it->eta() << " phi: " <<  it->phi() << " HoE: "<<it->hadronicOverEm()<<endl;
      bool passID = true;

      if (usePhotID_)
        {
          const pat::Photon *photon = dynamic_cast<const pat::Photon *>(&*it);
          if (!(photon->photonID(photID_)>0.)) passID = false;
        }

      if (it->pt()>photET_ && fabs(it->eta())<photEta_ && passID && it->hadronicOverEm()<photHoE_)
        {
          ++nphotons;
        }
    }

  if (debug_) cout <<"# Photons = "<<nphotons<<endl;
  // If not enough photons found (or too many found), return false
  // if (photMin_>-1 && nphotons<photMin_) return false;
  // if (photMax_>-1 && nphotons>photMax_) return false;


  // Step 1:  Count jets
  int njets=0;

  for (edm::View<reco::Candidate>::const_iterator it = jets->begin(); it != jets->end(); ++it)
    {
      if (debug_) cout << "Jet:" << endl;
      if (debug_) cout << "pT: " << it->pt() << " eta: " <<  it->eta() << " phi: " <<  it->phi() << endl;
      if (it->pt()>jetPT_ && fabs(it->eta())<jetEta_)
        {
          ++njets;
        }
      if (it->pt() > MaxRecJet) MaxRecJet = it->pt();
      // std::cout<<"Jet: "<<it->pt()<<std::endl;
    }

  if (debug_) cout <<"# Jets = "<<njets<<endl;
  // If not enough jets found (or too many found), return false
  // if (jetsMin_>-1 && njets<jetsMin_) return false;
  // if (jetsMax_>-1 && njets>jetsMax_) return false;

  // Step 2:  Count leptons
  int nmuons=0;
  int nelectrons=0;
  int ntaus=0;
  // count muons
  for (edm::View<reco::Candidate>::const_iterator it = muons->begin(); it != muons->end();++it)
    {
      if (debug_) cout << "Muon:" << endl;
      if (debug_) cout << "pT: " << it->pt() << " eta: " <<  it->eta() << " phi: " <<  it->phi() << endl;
      bool passID = true;

      if (useMuID_)
        {
          const pat::Muon *muon = dynamic_cast<const pat::Muon *>(&*it);
          if (!(muon->muonID(muID_))) passID = false;
        }

      if (it->pt()>muPT_ && fabs(it->eta())<muEta_ && passID)
        {
          ++nmuons;
        }
      if (it->pt() > MaxRecMuon) MaxRecMuon = it->pt();
      // std::cout<<"Muon: "<<it->pt()<<std::endl;
    }

  if (debug_) cout <<"# Muons = "<<nmuons<<endl;






  if( !iEvent.isRealData() ) {
    Handle<GenParticleCollection> genParticles;
    iEvent.getByLabel(genparticleLabel_, genParticles);

    for (size_t i = 0; i < genParticles->size(); ++i) {
      const Candidate & p = (*genParticles)[i];

        // Only gen muons
        if (abs(p.pdgId())!=13) continue;
        if (p.pt() > MaxGenMuon) MaxGenMuon = p.pt();

        // std::cout<<"GenMuon: "<<p.pt()<<std::endl;
    }

  }



  if( !iEvent.isRealData() ) {

    edm::Handle<reco::GenJetCollection> genJets;
    iEvent.getByLabel(genjetLabel_, genJets);

    for ( reco::GenJetCollection::const_iterator genJet = genJets->begin(); genJet != genJets->end(); ++genJet ) 
    {
    if (genJet->pt() > MaxGenJet) MaxGenJet = genJet->pt();
    // std::cout<<"GenJet: "<< genJet->pt()<<std::endl;

    }
  }






  // count electrons
  for (edm::View<reco::Candidate>::const_iterator it = electrons->begin(); it != electrons->end();++it)
    {
      if (debug_) cout << "Electron:" << endl;
      if (debug_) cout << "pT: " << it->pt() << " eta: " <<  it->eta() << " phi: " <<  it->phi() << endl;
      bool passID = true;

      if (useElecID_)
        {
          const pat::Electron *electron = dynamic_cast<const pat::Electron *>(&*it);
          if (!(electron->electronID(elecID_)>0.)) passID = false;
        }

      if (it->pt()>elecPT_ && fabs(it->eta())<elecEta_ && passID)
        {
          ++nelectrons;
        }
    }

  if (debug_) cout <<"# Electrons = "<<nelectrons<<endl;

  // count taus
  for (edm::View<reco::Candidate>::const_iterator it = taus->begin(); it != taus->end();++it){
   
    if (debug_) cout << "Tau:" << endl;
    if (debug_) cout << "pT: " << it->pt() << " eta: " <<  it->eta() << " phi: " <<  it->phi() << endl;
    bool passID = true;

    if ( useTauID_ ) { 
      const pat::Tau * tau = dynamic_cast<const pat::Tau*>(&*it);
      if ( tau->tauID(tauID_) < 0.5 ) passID = false;
    }

    if (it->pt()>tauPT_ && fabs(it->eta())<tauEta_ && passID) {
      ++ntaus;
    }   
    
  }

  if (customfilterGenRecoMuJet_==true)
  {
    // std::cout<<"CUSTOM FILTER ACTIVATED"<<std::endl;
    // std::cout<<"MaxMuon: "<<MaxRecMuon<<std::endl;
    // std::cout<<"MaxGenMuon: "<<MaxGenMuon<<std::endl;
    // std::cout<<"MaxJet: "<<MaxRecJet<<std::endl;
    // std::cout<<"MaxGenJet: "<<MaxGenJet<<std::endl;

    bool _hasmuon = (MaxGenMuon > 22.)||(MaxRecMuon > 22.);
    bool _hasjet = (MaxGenJet > 25.)||(MaxRecJet > 25.);
    bool _passcustomfilter = (_hasmuon && _hasjet);    
    // std::cout<<_hasmuon<<"  "<<_hasjet<<"  "<<_passcustomfilter<<std::endl;
    return _passcustomfilter;

  }

  // if (debug_) cout <<"# Tuas = "<<ntaus<<endl;

  // If we require both electron and muon condition to be met, check electron here
  // if (counteitherleptontype_==false)
  //   {
  //     if (electronsMin_>-1 && nelectrons<electronsMin_) return false;
  //     if (electronsMax_>-1 && nelectrons>electronsMax_) return false;
  //   }

  // // If we require both electron and muon condition to be met, check muon here
  // if (counteitherleptontype_==false)
  //   {
  //     if (muonsMin_>-1 && nmuons<muonsMin_) return false;
  //     if (muonsMax_>-1 && nmuons>muonsMax_) return false;
  //   }

  // if (counteitherleptontype_==false)
  //   {
  //     if (tausMin_>-1 && ntaus<tausMin_) return false;
  //     if (tausMax_>-1 && ntaus>tausMax_) return false; 
  //   }

  // // Otherwise, only fail cut if neither electron nor muon meet expectations
  // if (counteitherleptontype_==true)
  //   {
  //     if ((muonsMin_>-1 && nmuons<muonsMin_) && (electronsMin_>-1 && nelectrons<electronsMin_) && (tausMin_>-1 && ntaus<tausMin_)) return false;
  //     if ((muonsMax_>-1 && nmuons>muonsMax_) && (electronsMax_>-1 && nelectrons>electronsMax_) && (tausMax_>-1 && ntaus>tausMax_)) return false;
  //   }

  // ++PassedCount;
  // if (debug_) cout <<"PASSED!"<<endl;
   return true;
}

// ------------ method called once each job just before starting event loop  ------------
void
LeptonJetFilter::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void
LeptonJetFilter::endJob() {
  hCount->SetBinContent(1,TotalCount);
  hCount->SetBinContent(2,PassedCount);
  cout <<"Total events = "<<TotalCount<<"  Events Passed = "<<PassedCount<<endl;
}

//define this as a plug-in
DEFINE_FWK_MODULE(LeptonJetFilter);
