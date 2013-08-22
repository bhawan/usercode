#include "Leptoquarks/RootTupleMakerV2/interface/RootTupleMakerV2_GenParticles.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "TLorentzVector.h"


RootTupleMakerV2_GenParticles::RootTupleMakerV2_GenParticles(const edm::ParameterSet& iConfig) :
    inputTag(iConfig.getParameter<edm::InputTag>("InputTag")),
    prefix  (iConfig.getParameter<std::string>  ("Prefix")),
    suffix  (iConfig.getParameter<std::string>  ("Suffix")),
    maxSize (iConfig.getParameter<unsigned int> ("MaxSize"))
{
  produces <std::vector<float> > ( prefix + "Eta" + suffix );
  produces <std::vector<float> > ( prefix + "Phi" + suffix );
  produces <std::vector<float> > ( prefix + "P" + suffix );
  produces <std::vector<float> > ( prefix + "Px" + suffix );
  produces <std::vector<float> > ( prefix + "Py" + suffix );
  produces <std::vector<float> > ( prefix + "Pz" + suffix );
  produces <std::vector<float> > ( prefix + "Pt" + suffix );
  produces <std::vector<float> > ( prefix + "Energy" + suffix );
  produces <std::vector<int> >    ( prefix + "PdgId" + suffix );
  produces <std::vector<float> > ( prefix + "VX" + suffix );
  produces <std::vector<float> > ( prefix + "VY" + suffix );
  produces <std::vector<float> > ( prefix + "VZ" + suffix );
  produces <std::vector<int> >    ( prefix + "NumDaught" + suffix );
  produces <std::vector<int> >    ( prefix + "Status" + suffix );
  produces <std::vector<int> >    ( prefix + "MotherIndex" + suffix );
}

void RootTupleMakerV2_GenParticles::
produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {

  std::auto_ptr<std::vector<float> >  eta  ( new std::vector<float>()  );
  std::auto_ptr<std::vector<float> >  phi  ( new std::vector<float>()  );
  std::auto_ptr<std::vector<float> >  p  ( new std::vector<float>()  );
  std::auto_ptr<std::vector<float> >  px  ( new std::vector<float>()  );
  std::auto_ptr<std::vector<float> >  py  ( new std::vector<float>()  );
  std::auto_ptr<std::vector<float> >  pz  ( new std::vector<float>()  );
  std::auto_ptr<std::vector<float> >  pt  ( new std::vector<float>()  );
  std::auto_ptr<std::vector<float> >  energy  ( new std::vector<float>()  );
  std::auto_ptr<std::vector<int> >     pdgId ( new std::vector<int>()  );
  std::auto_ptr<std::vector<float> >  vx  ( new std::vector<float>()  );
  std::auto_ptr<std::vector<float> >  vy  ( new std::vector<float>()  );
  std::auto_ptr<std::vector<float> >  vz  ( new std::vector<float>()  );
  std::auto_ptr<std::vector<int> >     numDaught  ( new std::vector<int>()  );
  std::auto_ptr<std::vector<int> >     status  ( new std::vector<int>()  );
  std::auto_ptr<std::vector<int> >     motherIndex  ( new std::vector<int>()  );
  
  //-----------------------------------------------------------------

  std::vector<TLorentzVector> genmuons;
  std::vector<TLorentzVector> _muonsfromW;
  std::vector<TLorentzVector> _status1muonsfromW;
  std::vector<int> _muonsfromW_ids;
  std::vector<int> _status1muonsfromW_ids;

  std::vector<TLorentzVector> genelectrons;
  std::vector<TLorentzVector> _electronsfromW;
  std::vector<TLorentzVector> _status1electronsfromW;
  std::vector<int> _electronsfromW_ids;
  std::vector<int> _status1electronsfromW_ids;

  std::vector<int> pdgids;
  std::vector<TLorentzVector> genphotons;


  // std::cout<<" ------------- event ------------ "<<std::endl;
  if( !iEvent.isRealData() ) {
    edm::Handle<reco::GenParticleCollection> genParticles;
    iEvent.getByLabel(inputTag, genParticles);

    if( genParticles.isValid() ) {
      edm::LogInfo("RootTupleMakerV2_GenParticlesInfo") << "Total # GenParticles: " << genParticles->size();


      // Get status 1 photons
      for( reco::GenParticleCollection::const_iterator it = genParticles->begin(); it != genParticles->end(); ++it )
      {

                 // Only photons
        if (abs(it->pdgId())!=22) continue;
                 // Level - 1
        if (abs(it->status())!=1) continue;

        // continue; //sherpa

        TLorentzVector _photon;
        _photon.SetPtEtaPhiM(it->pt(), it->eta(), it->phi(),0.0);
        genphotons.push_back(_photon);
      }

      // Get muons & electrons which come from the W (status 3)
      for( reco::GenParticleCollection::const_iterator it = genParticles->begin(); it != genParticles->end(); ++it )
      {
        // Only muons from W
        if (abs(it->pdgId())!=13) continue;

        // Only electrons from W
        if (abs(it->pdgId())!=11) continue;

        // std::cout<<"Muon Pt/ID:  "<<it->pt()<<"  "<<it->pdgId()<<"  Status: "<<(it->status())<<"  MOTHER: "<<it->mother()->pdgId()<<std::endl;
        if (abs(it->mother()->pdgId()) != 24) continue; //madgraph
        // if (abs(it->status())!=3) continue; //sherpa

        TLorentzVector _muonfromW;
        TLorentzVector _electronfromW;
    
        _muonfromW.SetPtEtaPhiM(it->pt(), it->eta(), it->phi(),0.0);
        _muonsfromW.push_back(_muonfromW);
        _muonsfromW_ids.push_back(it->pdgId());
   
        _electronfromW.SetPtEtaPhiM(it->pt(), it->eta(), it->phi(),0.0);
        _electronsfromW.push_back(_electronfromW);
        _electronsfromW_ids.push_back(it->pdgId());


 }

      // Match status 1 muons to muons from W

      for( unsigned int im = 0 ; im!=_muonsfromW.size(); ++im )
      {
        TLorentzVector _closestmatchingmuon;
        float _closestdr = 9999999.99;
        int _closestmatchingmuon_id;
        for( reco::GenParticleCollection::const_iterator it = genParticles->begin(); it != genParticles->end(); ++it ) 
        {
          // Only status 1 muons
          if (abs(it->pdgId())!=13) continue;
          if (abs(it->status())!=1) continue; //madgraph
          // if (abs(it->status())!=3) continue; //sherpa

          TLorentzVector _status1muonfromW;

          _status1muonfromW.SetPtEtaPhiM(it->pt(), it->eta(), it->phi(),0.0);

          float _thisdr = _status1muonfromW.DeltaR(_muonsfromW[im]);
          if (_thisdr < _closestdr)
          {
            _closestdr = _thisdr;
            _closestmatchingmuon = _status1muonfromW;
            _closestmatchingmuon_id = it->pdgId();
          }
        }

        if (_closestdr < 0.3) 
        {
            _status1muonsfromW.push_back(_closestmatchingmuon);
            _status1muonsfromW_ids.push_back(_closestmatchingmuon_id);
        }
      }


      //Match status 1 electrons to electrons from W

for( unsigned int ie = 0 ; ie!=_electronsfromW.size(); ++ie )
      {
        TLorentzVector _closestmatchingelectron;
        float _closestdr = 9999999.99;
        int _closestmatchingelectron_id;
        for( reco::GenParticleCollection::const_iterator it = genParticles->begin(); it != genParticles->end(); ++it )
        {
          // Only status 1 electrons
          if (abs(it->pdgId())!=11) continue;
          if (abs(it->status())!=1) continue; //madgraph
          // if (abs(it->status())!=3) continue; //sherpa

          TLorentzVector _status1electronfromW;

          _status1electronfromW.SetPtEtaPhiM(it->pt(), it->eta(), it->phi(),0.0);

          float _thisdr = _status1electronfromW.DeltaR(_electronsfromW[ie]);
          if (_thisdr < _closestdr)
          {
            _closestdr = _thisdr;
            _closestmatchingelectron = _status1electronfromW;
            _closestmatchingelectron_id = it->pdgId();
          }
        }

        if (_closestdr < 0.3)
        {
            _status1electronsfromW.push_back(_closestmatchingelectron);
            _status1electronsfromW_ids.push_back(_closestmatchingelectron_id);
        }
      }


      for( unsigned int im = 0 ; im!=_status1muonsfromW.size(); ++im )
      {
        TLorentzVector _muon;
        int _muon_id = _status1muonsfromW_ids[im];
        _muon = _status1muonsfromW[im];
        std::vector<TLorentzVector> matchedphotons;

        // eta->push_back( _muon.Eta() );
        // phi->push_back( _muon.Phi() );
        // pt->push_back( _muon.Pt() );

        // std::cout<<"  Status 1 Muon:"<<_muon.Pt()<<"  "<<_muon.Eta()<<"  "<<_muon.Phi()<<std::endl;

        for (unsigned int ig = 0 ; ig!=genphotons.size(); ++ig)
        {
          if (fabs(_muon.DeltaR(genphotons[ig])) < 0.1)
          {
            matchedphotons.push_back(genphotons[ig]);
            // std::cout<<"   * "<<_muon.DeltaR(genphotons[ig])<<std::endl;
          }
        }

        for (unsigned int ig = 0 ; ig!=matchedphotons.size(); ++ig)
        {
          _muon = _muon+matchedphotons[ig];
        }
        // std::cout<<"   Dressed Muon:"<<_muon.Pt()<<"  "<<_muon.Eta()<<"  "<<_muon.Phi()<<std::endl;

        genmuons.push_back(_muon);
        int _newid = 0;
        if (_muon_id>0)  _newid = (700+_muon_id);
        if (_muon_id<0)  _newid = (-700+_muon_id);
        pdgids.push_back(_newid);
        // std::cout<<"  new ID: "<<_newid<<std::endl;
        // std::cout<<_muon.Pt()<<"  "<<700+it->pdgId()<<std::endl;

      }
   
for( unsigned int ie = 0 ; ie!=_status1electronsfromW.size(); ++ie )
      {
        TLorentzVector _electron;
        int _electron_id = _status1electronsfromW_ids[ie];
        _electron = _status1electronsfromW[ie];
        std::vector<TLorentzVector> matchedphotons;

        // std::cout<<"  Status 1 Muon:"<<_muon.Pt()<<"  "<<_muon.Eta()<<"  "<<_muon.Phi()<<std::endl;

        for (unsigned int ig = 0 ; ig!=genphotons.size(); ++ig)
        {
          if (fabs(_electron.DeltaR(genphotons[ig])) < 0.1)
          {
            matchedphotons.push_back(genphotons[ig]);
            // std::cout<<"   * "<<_muon.DeltaR(genphotons[ig])<<std::endl;
          }
        }

        for (unsigned int ig = 0 ; ig!=matchedphotons.size(); ++ig)
        {
          _electron = _electron+matchedphotons[ig];
        }
        // std::cout<<"   Dressed Muon:"<<_muon.Pt()<<"  "<<_muon.Eta()<<"  "<<_muon.Phi()<<std::endl;

        genelectrons.push_back(_electron);
        int _newid = 0;
        if (_electron_id>0)  _newid = (700+_electron_id);
        if (_electron_id<0)  _newid = (-700+_electron_id);
        pdgids.push_back(_newid);
        // std::cout<<"  new ID: "<<_newid<<std::endl;
        // std::cout<<_muon.Pt()<<"  "<<700+it->pdgId()<<std::endl;

      }


///////
 } 
    else 
    {
      edm::LogError("RootTupleMakerV2_GenParticlesError") << "Error! Can't get the product " << inputTag;
    }


  for (unsigned int imu = 0; imu!=genmuons.size(); ++imu)
  {
    eta->push_back( genmuons[imu].Eta() );
    phi->push_back( genmuons[imu].Phi() );
    pt->push_back( genmuons[imu].Pt() );
    pdgId->push_back( pdgids[imu] );
    motherIndex->push_back(799);
    status->push_back(799);
  }

 for (unsigned int iee = 0; iee!=genelectrons.size(); ++iee)
  {
    eta->push_back( genelectrons[iee].Eta() );
    phi->push_back( genelectrons[iee].Phi() );
    pt->push_back( genelectrons[iee].Pt() );
    pdgId->push_back( pdgids[iee] );
    motherIndex->push_back(799);
    status->push_back(799);
  }




 //-----------------------------------------------------------------
  if( !iEvent.isRealData() ) {
    edm::Handle<reco::GenParticleCollection> genParticles;
    iEvent.getByLabel(inputTag, genParticles);

    if( genParticles.isValid() ) {
      edm::LogInfo("RootTupleMakerV2_GenParticlesInfo") << "Total # GenParticles: " << genParticles->size();

      for( reco::GenParticleCollection::const_iterator it = genParticles->begin(); it != genParticles->end(); ++it ) {
        // exit from loop when you reach the required number of GenParticles
        if(eta->size() >= maxSize)
          break;

        int abid = abs(it->pdgId());
        bool _keep = false;

        if ( ( abid >= 11 ) &&  (abid <= 18 ) ) _keep = true;
        if (abid == 24) _keep = true;

        if (_keep == false) continue;

        // fill in all the vectors
        eta->push_back( it->eta() );
        phi->push_back( it->phi() );
        pt->push_back( it->pt() );
        pdgId->push_back( it->pdgId() );
        status->push_back( it->status() );
        motherIndex->push_back( it->mother()->pdgId() );
      }
    } else {
      edm::LogError("RootTupleMakerV2_GenParticlesError") << "Error! Can't get the product " << inputTag;
    }
  }

    
  }
  //-----------------------------------------------------------------
  // put vectors in the event
  iEvent.put( eta, prefix + "Eta" + suffix );
  iEvent.put( phi, prefix + "Phi" + suffix );
  iEvent.put( p, prefix + "P" + suffix );
  iEvent.put( px, prefix + "Px" + suffix );
  iEvent.put( py, prefix + "Py" + suffix );
  iEvent.put( pz, prefix + "Pz" + suffix );
  iEvent.put( pt, prefix + "Pt" + suffix );
  iEvent.put( energy, prefix + "Energy" + suffix );
  iEvent.put( pdgId, prefix + "PdgId" + suffix );
  iEvent.put( vx, prefix + "VX" + suffix );
  iEvent.put( vy, prefix + "VY" + suffix );
  iEvent.put( vz, prefix + "VZ" + suffix );
  iEvent.put( numDaught, prefix + "NumDaught" + suffix );
  iEvent.put( status, prefix + "Status" + suffix );
  iEvent.put( motherIndex, prefix + "MotherIndex" + suffix );
}
