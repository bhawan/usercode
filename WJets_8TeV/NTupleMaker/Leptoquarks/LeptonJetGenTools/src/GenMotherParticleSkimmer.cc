// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "RecoTauTag/TauTagTools/interface/GeneratorTau.h"

class GenMotherParticleSkimmer : public edm::EDProducer {
public:
  explicit GenMotherParticleSkimmer(const edm::ParameterSet&);
  ~GenMotherParticleSkimmer();
  
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
  
private:
  virtual void beginJob();
  virtual void produce(edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;
  
  virtual void beginRun(edm::Run&, edm::EventSetup const&);
  virtual void endRun(edm::Run&, edm::EventSetup const&);
  virtual void beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);
  virtual void endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);

  // user-defined information

  const edm::InputTag m_inputTag;
  const std::vector<int> m_motherPDGIds;
  const std::vector<int> m_motherPDGIdsVetoed;
  const std::vector<int> m_daughterPDGIds;
  const std::vector<int> m_daughterPDGStatuses;
  const bool m_storeFinalStateOnly;
  

};

GenMotherParticleSkimmer::GenMotherParticleSkimmer(const edm::ParameterSet& iConfig):
  m_inputTag            ( iConfig.getParameter<edm::InputTag>    ("InputTag")),
  m_motherPDGIds        ( iConfig.getParameter<std::vector<int> >("MotherPDGIds")),
  m_motherPDGIdsVetoed  ( iConfig.getParameter<std::vector<int> >("MotherPDGIdsVetoed")),
  m_daughterPDGIds      ( iConfig.getParameter<std::vector<int> >("DaughterPDGIds")),
  m_daughterPDGStatuses ( iConfig.getParameter<std::vector<int> >("DaughterPDGStatuses")),
  m_storeFinalStateOnly ( iConfig.getParameter<bool>             ("StoreFinalStateOnly"))
{
  produces<reco::GenParticleCollection>();
}

GenMotherParticleSkimmer::~GenMotherParticleSkimmer() {}

void GenMotherParticleSkimmer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
  
  std::auto_ptr<reco::GenParticleCollection> outputGenParticles ( new reco::GenParticleCollection());
  
  if( !iEvent.isRealData() ) {
    edm::Handle<reco::GenParticleCollection> inputGenParticles;
    iEvent.getByLabel(m_inputTag, inputGenParticles);
    
    if( inputGenParticles.isValid() ) {
      edm::LogInfo("GenMotherParticleSkimmer") << "Total # GenParticles: " << inputGenParticles->size();

      reco::GenParticleCollection::const_iterator it  = inputGenParticles->begin();
      reco::GenParticleCollection::const_iterator end = inputGenParticles->end();

      for (; it != end; ++it) { 
	
	//------------------------------------------------------------------------
	// Critical variables
	//------------------------------------------------------------------------

	int status = it -> status () ;
	int pdgid  = it -> pdgId  () ; 


	//------------------------------------------------------------------------
	// Is this one of the particle IDs that we want?
	//------------------------------------------------------------------------
	
	std::vector<int>::const_iterator interesting_id_iterator = std::find (m_daughterPDGIds.begin(), 
									      m_daughterPDGIds.end(), 
									      abs(pdgid) );
	

	if ( interesting_id_iterator == m_daughterPDGIds.end() ) continue;

	//------------------------------------------------------------------------
	// Is this one of the particle statuses that we want?
	//------------------------------------------------------------------------
	
	std::vector<int>::const_iterator interesting_status_iterator = std::find (m_daughterPDGStatuses.begin(), 
										  m_daughterPDGStatuses.end(), 
										  status );

	
	if ( interesting_status_iterator == m_daughterPDGStatuses.end() ) continue;

	//------------------------------------------------------------------------
	// Let's only store the final-state particles
	//------------------------------------------------------------------------

	if ( m_storeFinalStateOnly ) { 

	  reco::GenParticleRefVector::iterator iDaughter     = it -> daughterRefVector().begin();
	  reco::GenParticleRefVector::iterator iDaughter_end = it -> daughterRefVector().end();

	  bool final_state_particle = true;

	  for (; iDaughter != iDaughter_end; ++iDaughter) 
	    if ( (*iDaughter)->pdgId() == pdgid ) final_state_particle = false;
	  
	  if (!final_state_particle) continue;

	}

	//------------------------------------------------------------------------
	// Mother particle selection code
	//------------------------------------------------------------------------

	// variables of interest

	std::vector<int> motherChain;
	int numberOfMothers = 1;
	int maxNumberOfTries = 15;
	int numberOfTries = 0;
	bool foundInterestingMother = false;
	bool foundVetoedMother = false;

	// First: check that there is exactly one mother:
	
	if ( it -> numberOfMothers() != 1 ) continue;

	// Now check the rest of the decay chain
	
	const reco::Candidate* mother = it -> mother(0) ;
	
	while ( numberOfMothers == 1 && numberOfTries < maxNumberOfTries ) { 
	  
	  // Useful information about the mother particle

	  int mother_pdgid    = mother -> pdgId();
	  numberOfMothers = mother -> numberOfMothers();

	  // Is the mother of this particle in the list of mothers we want to save?
	  
	  std::vector<int>::const_iterator interesting_mother_iterator = std::find ( m_motherPDGIds.begin(),
	  									     m_motherPDGIds.end(),
	  									     abs(mother_pdgid) );
	  
	  foundInterestingMother = interesting_mother_iterator != m_motherPDGIds.end();
	  

	  //break the while loop for 'tries' if m_motherPDGId.size()==1 and it is already found.

          if( m_motherPDGIds.size()==1 && foundInterestingMother ) numberOfTries=maxNumberOfTries;


	  // Is the mother of this particle in the list of mothers we want to veto?
	  
	  std::vector<int>::const_iterator vetoed_mother_iterator = std::find ( m_motherPDGIdsVetoed.begin(),
										m_motherPDGIdsVetoed.end(),
										abs(mother_pdgid) );
	  
	  if(!foundVetoedMother) foundVetoedMother = vetoed_mother_iterator != m_motherPDGIdsVetoed.end();


	  // add mother PDG to chain of mothers

	  motherChain.push_back ( abs(mother_pdgid) );


	  // Look at the next mother, and iterate the number of tries
	  mother = mother -> mother(0);
	  numberOfTries++;
	  
	}
	
	//------------------------------------------------------------------------
	// If you never found an interesting mother, or if you found a vetoed mother: bail
	//------------------------------------------------------------------------

	if (  foundVetoedMother      ) continue;
	if ( !foundInterestingMother ) continue;
	
	//------------------------------------------------------------------------
	// Did we get ALL of the required mothers?
	//------------------------------------------------------------------------

	int nMotherPDGIds = m_motherPDGIds.size();
	bool found_all_requested_mother_pdgs = true;
	std::vector<int>::iterator motherChainIterator;


	for (int iMotherPDG = 0; iMotherPDG < nMotherPDGIds; ++iMotherPDG){
	  int requested_mother_pdgid = m_motherPDGIds[iMotherPDG];
	  motherChainIterator = std::find ( motherChain.begin(), motherChain.end(), abs(requested_mother_pdgid) ) ;
	  found_all_requested_mother_pdgs = found_all_requested_mother_pdgs && ( motherChainIterator != motherChain.end() ) ;
	}
	
	if ( ! found_all_requested_mother_pdgs ) continue;
	
	//------------------------------------------------------------------------
	// Success!  This is a good particle.  Add it to the collection!
	//------------------------------------------------------------------------

	outputGenParticles -> push_back ( *it );

	//------------------------------------------------------------------------
	// DEBUG: print mother chain
	//------------------------------------------------------------------------
	
	/*
	//int nMothersInChain = motherChain.size();
	if ( nMothersInChain != 0 ) { 
	std::cout << "DEBUG2:: Particle ID = " << pdgid << ", mother chain = ";
	for (int iMother = 0; iMother < nMothersInChain; ++iMother) { 
	std::cout << motherChain[iMother] << " ";
	}
	std::cout << std::endl;
	}
	*/
	  
      }
    }
  }

  iEvent.put ( outputGenParticles );
  
}

void GenMotherParticleSkimmer::beginJob(){}

void GenMotherParticleSkimmer::endJob() {}

void GenMotherParticleSkimmer::beginRun(edm::Run&, edm::EventSetup const&) {}

void GenMotherParticleSkimmer::endRun(edm::Run&, edm::EventSetup const&) {}

void GenMotherParticleSkimmer::beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&) {}

void GenMotherParticleSkimmer::endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&) {}

void GenMotherParticleSkimmer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(GenMotherParticleSkimmer);
