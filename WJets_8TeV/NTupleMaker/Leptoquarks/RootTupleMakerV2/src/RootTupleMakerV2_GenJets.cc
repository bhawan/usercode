#include "Leptoquarks/RootTupleMakerV2/interface/RootTupleMakerV2_GenJets.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "DataFormats/JetReco/interface/GenJetCollection.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "TLorentzVector.h"


RootTupleMakerV2_GenJets::RootTupleMakerV2_GenJets(const edm::ParameterSet& iConfig) :
    inputTag(iConfig.getParameter<edm::InputTag>("InputTag")),
    //inputTagP(iConfig.getParameter<edm::InputTag>("InputTagP")),
    prefix  (iConfig.getParameter<std::string>  ("Prefix")),
    suffix  (iConfig.getParameter<std::string>  ("Suffix")),
    maxSize (iConfig.getParameter<unsigned int> ("MaxSize"))
{
  produces <std::vector<float> > ( prefix + "Eta" + suffix );
  produces <std::vector<float> > ( prefix + "Phi" + suffix );
  produces <std::vector<float> > ( prefix + "P" + suffix );
  produces <std::vector<float> > ( prefix + "Pt" + suffix );
  produces <std::vector<float> > ( prefix + "Energy" + suffix );
  produces <std::vector<float> > ( prefix + "EMF" + suffix );
  produces <std::vector<float> > ( prefix + "HADF" + suffix );
  produces <std::vector<float> > ( prefix + "NUF" + suffix );
  produces <std::vector<float> > ( prefix + "LEPF" + suffix );

}

void RootTupleMakerV2_GenJets::
produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
  // using namespace std;
  std::auto_ptr<std::vector<float> >  eta  ( new std::vector<float>()  );
  std::auto_ptr<std::vector<float> >  phi  ( new std::vector<float>()  );
  std::auto_ptr<std::vector<float> >  p  ( new std::vector<float>()  );
  std::auto_ptr<std::vector<float> >  pt  ( new std::vector<float>()  );
  std::auto_ptr<std::vector<float> >  energy  ( new std::vector<float>()  );
  std::auto_ptr<std::vector<float> >  emf  ( new std::vector<float>()  );
  std::auto_ptr<std::vector<float> >  hadf  ( new std::vector<float>()  );
  std::auto_ptr<std::vector<float> >  nuf  ( new std::vector<float>()  );
  std::auto_ptr<std::vector<float> >  lepf  ( new std::vector<float>()  );

  //-----------------------------------------------------------------
  if( !iEvent.isRealData() ) {
   // edm::Handle<reco::GenParticleCollection> genParticlesP;
 //   iEvent.getByLabel(inputTagP, genParticlesP);

    edm::Handle<reco::GenJetCollection> genJets;
    iEvent.getByLabel(inputTag, genJets);

    if( genJets.isValid() ) {
      edm::LogInfo("RootTupleMakerV2_GenJetsInfo") << "Total # GenJets: " << genJets->size();
      // std::cout<<" --------------------------------------------------- "<<std::endl;
      for( reco::GenJetCollection::const_iterator it = genJets->begin(); it != genJets->end(); ++it ) {
        // exit from loop when you reach the required number of GenJets
        if(eta->size() >= maxSize)
          break;

        float nufrac = 0.0;
        float lepfrac = 0.0;
        TLorentzVector j;
        j.SetPtEtaPhiM(it->pt(),it->eta(),it->phi(),0.0);

        // Get status 1 photons
        // for( reco::GenParticleCollection::const_iterator itp = genParticlesP->begin(); itp != genParticlesP->end(); ++itp )
        // {
        //   TLorentzVector _p;
        //   _p.SetPtEtaPhiM(itp->pt(), itp->eta(), itp->phi(),0.0);

        //   std::cout<<"   -- "<<_p.DeltaR(j)<<"    "<<itp->pdgId()<<std::endl;
        // std::cout<<"   -- "<<(it->getGenConstituents())->size()<<std::endl;
        // }
        std::vector<const reco::GenParticle*> theJetConstituents  = it->getGenConstituents();


        for( std::vector<const reco::GenParticle*>::const_iterator aCandidate = theJetConstituents.begin();
              aCandidate != theJetConstituents.end();
              ++aCandidate)
          {
             int pdgId = std::abs((*aCandidate)->pdgId());
             float cpt = (*aCandidate)->energy()/it->energy();
             // std::cout<<  "  --> "<<cpt<<"  "<<pdgId<<std::endl;
             if ((pdgId==12)||(pdgId == 14)||(pdgId == 16)||(pdgId == 18))
             {
              // std::cout<<"                         NEUTRINO!!"<<std::endl;
              nufrac += cpt;

              TLorentzVector _p;
              _p.SetPtEtaPhiM((*aCandidate)->pt(), (*aCandidate)->eta(), (*aCandidate)->phi(),0.0);

              // if ( (_p.DeltaR(j) < 0.05) &&  ( fabs(_p.Pt() - j.Pt())/j.Pt() < 0.05 )) 
              //   std::cout<<"                             --> NU FLAG: "<<_p.DeltaR(j)<<"  "<<fabs(_p.Pt() - j.Pt())/j.Pt() <<std::endl;

              }
             if ((pdgId==11)||(pdgId == 13)||(pdgId == 15)||(pdgId == 17))
             {
              // std::cout<<"                         CHARGED LEPTON !!"<<std::endl;
              lepfrac += cpt;
              }
          }

        // fill in all the vectors
        eta->push_back( it->eta() );
        phi->push_back( it->phi() );
        p->push_back( it->p() );
        pt->push_back( it->pt() );
        energy->push_back( it->energy() );
        emf->push_back( it->emEnergy()/it->energy() );
        hadf->push_back( it->hadEnergy()/it->energy() );
        nuf->push_back(nufrac);
        lepf->push_back(lepfrac);

        // std::cout<<(1.0-  it->emEnergy()/it->energy()  -it->hadEnergy()/it->energy()  )<<std::endl;
        // std::cout<<nufrac<<std::endl;
        // std::cout<<lepfrac<<std::endl;

      }
    } else {
      edm::LogError("RootTupleMakerV2_GenJetsError") << "Error! Can't get the product " << inputTag;
    }
  }

  //-----------------------------------------------------------------
  // put vectors in the event
  iEvent.put( eta, prefix + "Eta" + suffix );
  iEvent.put( phi, prefix + "Phi" + suffix );
  iEvent.put( p, prefix + "P" + suffix );
  iEvent.put( pt, prefix + "Pt" + suffix );
  iEvent.put( energy, prefix + "Energy" + suffix );
  iEvent.put( emf, prefix + "EMF" + suffix );
  iEvent.put( hadf, prefix + "HADF" + suffix );
  iEvent.put( nuf, prefix + "NUF" + suffix );
  iEvent.put( lepf, prefix + "LEPF" + suffix );

}
