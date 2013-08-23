#include "Leptoquarks/RootTupleMakerV2/interface/RootTupleMakerV2_MET.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "DataFormats/PatCandidates/interface/MET.h"

RootTupleMakerV2_MET::RootTupleMakerV2_MET(const edm::ParameterSet& iConfig) :
    inputTag(iConfig.getParameter<edm::InputTag>("InputTag")),
    prefix  (iConfig.getParameter<std::string>  ("Prefix")),
    suffix  (iConfig.getParameter<std::string>  ("Suffix")),
    store_uncorrected_MET (iConfig.getParameter<bool>  ("StoreUncorrectedMET")),
    store_MET_significance (iConfig.getParameter<bool>  ("StoreMETSignificance"))
{
  produces <std::vector<float> > ( prefix + "MET" + suffix );
  produces <std::vector<float> > ( prefix + "METPhi" + suffix );
  produces <std::vector<float> > ( prefix + "SumET" + suffix );
  if ( store_uncorrected_MET ) {
    produces <std::vector<float> > ( prefix + "METUncorr" + suffix );
    produces <std::vector<float> > ( prefix + "METPhiUncorr" + suffix );
    produces <std::vector<float> > ( prefix + "SumETUncorr" + suffix );
  }
  if ( store_MET_significance ) {
    produces <std::vector<float> > ( prefix + "METSig" + suffix );
    produces <std::vector<float> > ( prefix + "METSigMatrixDXX" + suffix );
    produces <std::vector<float> > ( prefix + "METSigMatrixDXY" + suffix );
    produces <std::vector<float> > ( prefix + "METSigMatrixDYX" + suffix );
    produces <std::vector<float> > ( prefix + "METSigMatrixDYY" + suffix );
  }
}

void RootTupleMakerV2_MET::
produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {

  std::auto_ptr<std::vector<float> >  met  ( new std::vector<float>()  );
  std::auto_ptr<std::vector<float> >  metphi  ( new std::vector<float>()  );
  std::auto_ptr<std::vector<float> >  sumet  ( new std::vector<float>()  );
  std::auto_ptr<std::vector<float> >  metuncorr  ( new std::vector<float>()  );
  std::auto_ptr<std::vector<float> >  metphiuncorr  ( new std::vector<float>()  );
  std::auto_ptr<std::vector<float> >  sumetuncorr  ( new std::vector<float>()  );
  std::auto_ptr<std::vector<float> >  metsig  ( new std::vector<float>()  );  
  std::auto_ptr<std::vector<float> >  metsigmatrixdxx  ( new std::vector<float>()  );  
  std::auto_ptr<std::vector<float> >  metsigmatrixdxy  ( new std::vector<float>()  );  
  std::auto_ptr<std::vector<float> >  metsigmatrixdyx  ( new std::vector<float>()  );  
  std::auto_ptr<std::vector<float> >  metsigmatrixdyy  ( new std::vector<float>()  );  

  //-----------------------------------------------------------------
  edm::Handle<std::vector<pat::MET> > mets;
  iEvent.getByLabel(inputTag, mets);

  if(mets.isValid()) {
    edm::LogInfo("RootTupleMakerV2_METInfo") << "Total # METs: " << mets->size();

    for( std::vector<pat::MET>::const_iterator it = mets->begin(); it != mets->end(); ++it ) {

      // fill in all the vectors
      met->push_back( it->pt() );
      metphi->push_back( it->phi() );
      sumet->push_back( it->sumEt() );
      
      if ( store_uncorrected_MET ) {
	metuncorr->push_back( it->uncorrectedPt(pat::MET::uncorrALL) );
	metphiuncorr->push_back( it->uncorrectedPhi(pat::MET::uncorrALL) );
	sumetuncorr->push_back( it->sumEt() - it->corSumEt(pat::MET::uncorrALL) );
      }

      if ( store_MET_significance ) {
        //-- See https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideMETSignificance#Known_Issues
	float sigmaX2= it->getSignificanceMatrix()(0,0);
	float sigmaY2= it->getSignificanceMatrix()(1,1);
	float significance = -1;
	if(sigmaX2<1.e10 && sigmaY2<1.e10) significance = it->significance();
	//--

	metsig->push_back( significance );
	metsigmatrixdxx->push_back( it->getSignificanceMatrix()(0,0) );
	metsigmatrixdxy->push_back( it->getSignificanceMatrix()(0,1) );
	metsigmatrixdyx->push_back( it->getSignificanceMatrix()(1,0) );
	metsigmatrixdyy->push_back( it->getSignificanceMatrix()(1,1) );
	//See DataFormats/METReco/src/MET.cc
      }

    }
  } else {
    edm::LogError("RootTupleMakerV2_METError") << "Error! Can't get the product " << inputTag;
  }

  //-----------------------------------------------------------------
  // put vectors in the event
  iEvent.put( met, prefix + "MET" + suffix );
  iEvent.put( metphi, prefix + "METPhi" + suffix );
  iEvent.put( sumet, prefix + "SumET" + suffix );
  if ( store_uncorrected_MET ) {
    iEvent.put( metuncorr, prefix + "METUncorr" + suffix );
    iEvent.put( metphiuncorr, prefix + "METPhiUncorr" + suffix );
    iEvent.put( sumetuncorr, prefix + "SumETUncorr" + suffix );
  }
  if ( store_MET_significance ) {
    iEvent.put( metsig, prefix + "METSig" + suffix );
    iEvent.put( metsigmatrixdxx, prefix + "METSigMatrixDXX" + suffix );
    iEvent.put( metsigmatrixdxy, prefix + "METSigMatrixDXY" + suffix );
    iEvent.put( metsigmatrixdyx, prefix + "METSigMatrixDYX" + suffix );
    iEvent.put( metsigmatrixdyy, prefix + "METSigMatrixDYY" + suffix );
  }
}
