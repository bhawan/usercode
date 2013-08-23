#include "Leptoquarks/RootTupleMakerV2/interface/RootTupleMakerV2_Taus.h"
#include "Leptoquarks/RootTupleMakerV2/interface/PatUtilities.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Tau.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/Scalers/interface/DcsStatus.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "RecoEgamma/EgammaTools/interface/ConversionFinder.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include <DataFormats/Common/interface/Ref.h>
#include "DataFormats/TauReco/interface/PFTauDiscriminator.h"
#include "DataFormats/TauReco/interface/PFTau.h"

RootTupleMakerV2_Taus::RootTupleMakerV2_Taus(const edm::ParameterSet& iConfig) :
  inputTag    (iConfig.getParameter<edm::InputTag>("InputTag")),
  vtxInputTag (iConfig.getParameter<edm::InputTag>("VertexInputTag")),
  prefix      (iConfig.getParameter<std::string>("Prefix")),
  suffix      (iConfig.getParameter<std::string>("Suffix")),
  maxSize     (iConfig.getParameter<unsigned int>("MaxSize")),
  isSCTau     (iConfig.getParameter<bool>("isSCTau")),
  isHPSTau    (iConfig.getParameter<bool>("isHPSTau"))
{
/*  produces <std::vector<float> > ( prefix + "Eta"                              + suffix );
  produces <std::vector<float> > ( prefix + "Phi"                              + suffix );
  produces <std::vector<float> > ( prefix + "Pt"                               + suffix );
  produces <std::vector<float> > ( prefix + "Et"                               + suffix );
  produces <std::vector<int> >    ( prefix + "Charge"                           + suffix );
  produces <std::vector<int> >    ( prefix + "IsPFTau"                          + suffix );
  produces <std::vector<int> >    ( prefix + "IsCaloTau"                        + suffix );
  produces <std::vector<int> >    ( prefix + "DecayMode"                        + suffix );
  produces <std::vector<float> > ( prefix + "EmFraction"                       + suffix );
  produces <std::vector<float> > ( prefix + "Hcal3x3OverPLead"                 + suffix );
  produces <std::vector<float> > ( prefix + "HcalMaxOverPLead"                 + suffix );
  produces <std::vector<float> > ( prefix + "HcalTotOverPLead"                 + suffix );
  produces <std::vector<float> > ( prefix + "IsolationPFChargedHadrCandsPtSum" + suffix );
  produces <std::vector<float> > ( prefix + "IsolationPFGammaCandsEtSum"       + suffix );
  produces <std::vector<float> > ( prefix + "LeadPFChargedHadrCandsignedSipt"  + suffix );
  produces <std::vector<float> > ( prefix + "EtaLeadCharged"                   + suffix );
  produces <std::vector<float> > ( prefix + "PhiLeadCharged"                   + suffix );
  produces <std::vector<float> > ( prefix + "PtLeadCharged"                    + suffix );
  produces <std::vector<float> > ( prefix + "PhiphiMoment"                     + suffix );
  produces <std::vector<float> > ( prefix + "EtaetaMoment"                     + suffix );
  produces <std::vector<float> > ( prefix + "EtaphiMoment"                     + suffix );
  produces <std::vector<float> > ( prefix + "EcalStripSumEOverPLead"           + suffix ); 
  produces <std::vector<float> > ( prefix + "BremsRecoveryEOverPLead"          + suffix );
  produces <std::vector<float> > ( prefix + "MaximumHCALPFClusterEt"           + suffix ); 
  produces <std::vector<float> > ( prefix + "MatchedGenParticlePt"             + suffix ); 
  produces <std::vector<float> > ( prefix + "MatchedGenParticleEta"            + suffix ); 
  produces <std::vector<float> > ( prefix + "MatchedGenParticlePhi"            + suffix ); 
  produces <std::vector<float> > ( prefix + "MatchedGenJetPt"                  + suffix ); 
  produces <std::vector<float> > ( prefix + "MatchedGenJetEta"                 + suffix ); 
  produces <std::vector<float> > ( prefix + "MatchedGenJetPhi"                 + suffix ); 
*/ 
 //
  // ShrinkingCone PFTau Specific
  if(isSCTau){
    //shrinkingCone PFTau Discriminators (SCTau)
  /*  produces <std::vector<float> >    ( prefix + "LeadingTrackFindingDiscr"            + suffix );
    produces <std::vector<float> >    ( prefix + "LeadingTrackPtCutDiscr"              + suffix );
    produces <std::vector<float> >    ( prefix + "LeadingPionPtCutDiscr"               + suffix );
    produces <std::vector<float> >    ( prefix + "IsolationDiscr"                      + suffix );
    produces <std::vector<float> >    ( prefix + "TrackIsolationDiscr"                 + suffix );
    produces <std::vector<float> >    ( prefix + "EcalIsolationDiscr"                  + suffix );
    produces <std::vector<float> >    ( prefix + "IsolationUsingLeadingPionDiscr"      + suffix );
    produces <std::vector<float> >    ( prefix + "TrackIsolationUsingLeadingPionDiscr" + suffix );
    produces <std::vector<float> >    ( prefix + "EcalIsolationUsingLeadingPionDiscr"  + suffix );
    produces <std::vector<float> >    ( prefix + "AgainstElectronDiscr"                + suffix );
    produces <std::vector<float> >    ( prefix + "AgainstMuonDiscr"                    + suffix );
    produces <std::vector<float> >    ( prefix + "TaNCDiscr"                           + suffix );
    produces <std::vector<float> >    ( prefix + "TaNCfrOnePercentDiscr"               + suffix );
    produces <std::vector<float> >    ( prefix + "TaNCfrHalfPercentDiscr"              + suffix );
    produces <std::vector<float> >    ( prefix + "TaNCfrQuarterPercentDiscr"           + suffix );
    produces <std::vector<float> >    ( prefix + "TaNCfrTenthPercentDiscr"             + suffix );
*/
  }
  //
  // HPS PFTau Specific
  if(isHPSTau){
    //hps PFTau Discriminators (HPSTau)
/*   produces <std::vector<float> >    ( prefix + "DecayModeFindingDiscr"  + suffix );
    //
    produces <std::vector<float> >    ( prefix + "AgainstElectronLooseDiscr"        + suffix );
    produces <std::vector<float> >    ( prefix + "AgainstElectronMediumDiscr"       + suffix );
    produces <std::vector<float> >    ( prefix + "AgainstElectronTightDiscr"        + suffix );
    produces <std::vector<float> >    ( prefix + "AgainstElectronMVADiscr"          + suffix );
    produces <std::vector<float> >    ( prefix + "AgainstElectronMVA2rawDiscr"      + suffix );
    produces <std::vector<float> >    ( prefix + "AgainstElectronMVA2categoryDiscr" + suffix );
    produces <std::vector<float> >    ( prefix + "AgainstElectronVLooseMVA2Discr"   + suffix );
    produces <std::vector<float> >    ( prefix + "AgainstElectronLooseMVA2Discr"    + suffix );
    produces <std::vector<float> >    ( prefix + "AgainstElectronMediumMVA2Discr"   + suffix );
    produces <std::vector<float> >    ( prefix + "AgainstElectronTightMVA2Discr"    + suffix );
    produces <std::vector<float> >    ( prefix + "AgainstElectronMVA3rawDiscr"      + suffix );
    produces <std::vector<float> >    ( prefix + "AgainstElectronMVA3categoryDiscr" + suffix );
    produces <std::vector<float> >    ( prefix + "AgainstElectronLooseMVA3Discr"    + suffix );
    produces <std::vector<float> >    ( prefix + "AgainstElectronMediumMVA3Discr"   + suffix );
    produces <std::vector<float> >    ( prefix + "AgainstElectronTightMVA3Discr"    + suffix );
    produces <std::vector<float> >    ( prefix + "AgainstElectronVTightMVA3Discr"   + suffix );
    produces <std::vector<float> >    ( prefix + "AgainstElectronDeadECALDiscr"     + suffix );
    //
    produces <std::vector<float> >    ( prefix + "AgainstMuonLooseDiscr"   + suffix );
    produces <std::vector<float> >    ( prefix + "AgainstMuonMediumDiscr"  + suffix );
    produces <std::vector<float> >    ( prefix + "AgainstMuonTightDiscr"   + suffix );
    produces <std::vector<float> >    ( prefix + "AgainstMuonLoose2Discr"  + suffix );
    produces <std::vector<float> >    ( prefix + "AgainstMuonMedium2Discr" + suffix );
    produces <std::vector<float> >    ( prefix + "AgainstMuonTight2Discr"  + suffix );
    //
    produces <std::vector<float> >    ( prefix + "VLooseIsolationDiscr"                           + suffix );
    produces <std::vector<float> >    ( prefix + "LooseIsolationDiscr"                            + suffix );
    produces <std::vector<float> >    ( prefix + "MediumIsolationDiscr"                           + suffix );
    produces <std::vector<float> >    ( prefix + "TightIsolationDiscr"                            + suffix );
    produces <std::vector<float> >    ( prefix + "VLooseIsolationDeltaBetaCorrDiscr"              + suffix );
    produces <std::vector<float> >    ( prefix + "LooseIsolationDeltaBetaCorrDiscr"               + suffix );
    produces <std::vector<float> >    ( prefix + "MediumIsolationDeltaBetaCorrDiscr"              + suffix );
    produces <std::vector<float> >    ( prefix + "TightIsolationDeltaBetaCorrDiscr"               + suffix );
    produces <std::vector<float> >    ( prefix + "VLooseCombinedIsolationDeltaBetaCorrDiscr"      + suffix );
    produces <std::vector<float> >    ( prefix + "LooseCombinedIsolationDeltaBetaCorrDiscr"       + suffix );
    produces <std::vector<float> >    ( prefix + "MediumCombinedIsolationDeltaBetaCorrDiscr"      + suffix );
    produces <std::vector<float> >    ( prefix + "TightCombinedIsolationDeltaBetaCorrDiscr"       + suffix );
    produces <std::vector<float> >    ( prefix + "CombinedIsolationDeltaBetaCorr3HitsDiscr"       + suffix );
    produces <std::vector<float> >    ( prefix + "LooseCombinedIsolationDeltaBetaCorr3HitsDiscr"  + suffix );
    produces <std::vector<float> >    ( prefix + "MediumCombinedIsolationDeltaBetaCorr3HitsDiscr" + suffix );
    produces <std::vector<float> >    ( prefix + "TightCombinedIsolationDeltaBetaCorr3HitsDiscr"  + suffix );
    produces <std::vector<float> >    ( prefix + "IsolationMVArawDiscr"                           + suffix );
    produces <std::vector<float> >    ( prefix + "LooseIsolationMVADiscr"                         + suffix );
    produces <std::vector<float> >    ( prefix + "MediumIsolationMVADiscr"                        + suffix );
    produces <std::vector<float> >    ( prefix + "TightIsolationMVADiscr"                         + suffix );
    produces <std::vector<float> >    ( prefix + "LooseIsolationMVA2Discr"                        + suffix );
    produces <std::vector<float> >    ( prefix + "MediumIsolationMVA2Discr"                       + suffix );
    produces <std::vector<float> >    ( prefix + "TightIsolationMVA2Discr"                        + suffix );
    //
    // HPSTau Signal PFCandidates Info
    produces <std::vector<float> >    ( prefix + "SignalPFChargedHadrCandsPt"    + suffix );
    produces <std::vector<float> >    ( prefix + "SignalPFChargedHadrCandsEta"   + suffix );
    produces <std::vector<float> >    ( prefix + "SignalPFChargedHadrCandsPhi"   + suffix );
    produces <std::vector<float> >    ( prefix + "SignalPFChargedHadrCandsCount" + suffix );
    produces <std::vector<float> >    ( prefix + "SignalPFNeutrHadrCandsPt"      + suffix );
    produces <std::vector<float> >    ( prefix + "SignalPFNeutrHadrCandsEta"     + suffix );
    produces <std::vector<float> >    ( prefix + "SignalPFNeutrHadrCandsPhi"     + suffix );
    produces <std::vector<float> >    ( prefix + "SignalPFNeutrHadrCandsCount"   + suffix );
    produces <std::vector<float> >    ( prefix + "SignalPFGammaCandsPt"          + suffix );
    produces <std::vector<float> >    ( prefix + "SignalPFGammaCandsEta"         + suffix );
    produces <std::vector<float> >    ( prefix + "SignalPFGammaCandsPhi"         + suffix );
    produces <std::vector<float> >    ( prefix + "SignalPFGammaCandsCount"       + suffix );
    //
    // HPSTau Vertex Info
    produces <std::vector<int> >    ( prefix + "VtxIndex"      + suffix );
    produces <std::vector<float> > ( prefix + "VtxDistXY"     + suffix );
    produces <std::vector<float> > ( prefix + "VtxDistZ"      + suffix );
    produces <std::vector<float> > ( prefix + "LeadVtxDistXY" + suffix );
    produces <std::vector<float> > ( prefix + "LeadVtxDistZ"  + suffix );
    //
    // --------------------------------------------------------------------------------------- //
    // HPS Tau Optional Isolation information
    //produces <std::vector<float> >    ( prefix + "IsolationPFChargedHadrCandsPt"    + suffix );
    //produces <std::vector<float> >    ( prefix + "IsolationPFChargedHadrCandsEta"   + suffix );
    //produces <std::vector<float> >    ( prefix + "IsolationPFChargedHadrCandsPhi"   + suffix );
    //produces <std::vector<float> >    ( prefix + "IsolationPFChargedHadrCandsCount" + suffix );
    //produces <std::vector<float> >    ( prefix + "IsolationPFNeutrHadrCandsPt"      + suffix );
    //produces <std::vector<float> >    ( prefix + "IsolationPFNeutrHadrCandsEta"     + suffix );
    //produces <std::vector<float> >    ( prefix + "IsolationPFNeutrHadrCandsPhi"     + suffix );
    //produces <std::vector<float> >    ( prefix + "IsolationPFNeutrHadrCandsCount"   + suffix );
    //produces <std::vector<float> >    ( prefix + "IsolationPFGammaCandsPt"          + suffix );
    //produces <std::vector<float> >    ( prefix + "IsolationPFGammaCandsEta"         + suffix );
    //produces <std::vector<float> >    ( prefix + "IsolationPFGammaCandsPhi"         + suffix );
    //produces <std::vector<float> >    ( prefix + "IsolationPFGammaCandsCount"       + suffix );
    // --------------------------------------------------------------------------------------- //
*/    //
  }
  //
}

void RootTupleMakerV2_Taus::
produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
  //
  std::auto_ptr<std::vector<float> >  eta  ( new std::vector<float>()  );
  std::auto_ptr<std::vector<float> >  phi  ( new std::vector<float>()  );
  std::auto_ptr<std::vector<float> >  pt  ( new std::vector<float>()   );
  std::auto_ptr<std::vector<float> >  et  ( new std::vector<float>()   );
  std::auto_ptr<std::vector<int> >     charge  ( new std::vector<int>()   );
  std::auto_ptr<std::vector<int> >     ispftau  ( new std::vector<int>()   );
  std::auto_ptr<std::vector<int> >     iscalotau  ( new std::vector<int>()   );
  std::auto_ptr<std::vector<int> >     decaymode  ( new std::vector<int>()   );
  std::auto_ptr<std::vector<float> >  emfraction  ( new std::vector<float>()   );
  std::auto_ptr<std::vector<float> >  hcal3x3overplead  ( new std::vector<float>()   );
  std::auto_ptr<std::vector<float> >  hcalmaxoverplead  ( new std::vector<float>()   );
  std::auto_ptr<std::vector<float> >  hcaltotoverplead  ( new std::vector<float>()   );
  std::auto_ptr<std::vector<float> >  isolationpfchargedhadrcandsptsum  ( new std::vector<float>()   );
  std::auto_ptr<std::vector<float> >  isolationpfgammacandsetsum  ( new std::vector<float>()   );
  std::auto_ptr<std::vector<float> >  leadpfchargedhadrcandsignedsipt ( new std::vector<float>()   );
  std::auto_ptr<std::vector<float> >  etaleadcharged  ( new std::vector<float>()   );
  std::auto_ptr<std::vector<float> >  phileadcharged  ( new std::vector<float>()   );
  std::auto_ptr<std::vector<float> >  ptleadcharged  ( new std::vector<float>()   );
  std::auto_ptr<std::vector<float> >  phiphimoment  ( new std::vector<float>()   );
  std::auto_ptr<std::vector<float> >  etaetamoment  ( new std::vector<float>()   );
  std::auto_ptr<std::vector<float> >  etaphimoment  ( new std::vector<float>()   );
  std::auto_ptr<std::vector<float> >  ecalstripsumeoverplead  ( new std::vector<float>()   );
  std::auto_ptr<std::vector<float> >  bremsrecoveryeoverplead  ( new std::vector<float>()   );
  std::auto_ptr<std::vector<float> >  maximumhcalpfclusteret  ( new std::vector<float>()   );
  std::auto_ptr<std::vector<float> >  matchedgenparticlept ( new std::vector<float>()   );
  std::auto_ptr<std::vector<float> >  matchedgenparticleeta ( new std::vector<float>()   );
  std::auto_ptr<std::vector<float> >  matchedgenparticlephi ( new std::vector<float>()   );
  std::auto_ptr<std::vector<float> >  matchedgenjetpt ( new std::vector<float>()   );
  std::auto_ptr<std::vector<float> >  matchedgenjeteta ( new std::vector<float>()   );
  std::auto_ptr<std::vector<float> >  matchedgenjetphi ( new std::vector<float>()   );
  //
  //shrinkingCone PFTau Discriminators (SCTau)
  std::auto_ptr<std::vector<float> >  leadingtrackfindingdiscr  ( new std::vector<float>()   );
  std::auto_ptr<std::vector<float> >  leadingtrackptcutdiscr  ( new std::vector<float>()   );
  std::auto_ptr<std::vector<float> >  leadingpionptcutdiscr  ( new std::vector<float>()   );
  std::auto_ptr<std::vector<float> >  isolationdiscr  ( new std::vector<float>()   );
  std::auto_ptr<std::vector<float> >  trackisolationdiscr  ( new std::vector<float>()   );
  std::auto_ptr<std::vector<float> >  ecalisolationdiscr  ( new std::vector<float>()   );
  std::auto_ptr<std::vector<float> >  isolationusingleadingpiondiscr  ( new std::vector<float>()   );
  std::auto_ptr<std::vector<float> >  trackisolationusingleadingpiondiscr  ( new std::vector<float>()   );
  std::auto_ptr<std::vector<float> >  ecalisolationusingleadingpiondiscr  ( new std::vector<float>()   );
  std::auto_ptr<std::vector<float> >  againstelectrondiscr  ( new std::vector<float>()   );
  std::auto_ptr<std::vector<float> >  againstmuondiscr  ( new std::vector<float>()   );
  std::auto_ptr<std::vector<float> >  tancdiscr  ( new std::vector<float>()   );
  std::auto_ptr<std::vector<float> >  tancfronepercentdiscr  ( new std::vector<float>()   );
  std::auto_ptr<std::vector<float> >  tancfrhalfpercentdiscr  ( new std::vector<float>()   );
  std::auto_ptr<std::vector<float> >  tancfrquarterpercentdiscr  ( new std::vector<float>()   );
  std::auto_ptr<std::vector<float> >  tancfrtenthpercentdiscr  ( new std::vector<float>()   );
  //
  //hps PFTau Discriminators (HPSTau)
  std::auto_ptr<std::vector<float> >  decaymodefindingdiscr  ( new std::vector<float>()   );
  //
  std::auto_ptr<std::vector<float> >  againstelectronloosediscr        ( new std::vector<float>()   );
  std::auto_ptr<std::vector<float> >  againstelectronmediumdiscr       ( new std::vector<float>()   );
  std::auto_ptr<std::vector<float> >  againstelectrontightdiscr        ( new std::vector<float>()   );
  std::auto_ptr<std::vector<float> >  againstelectronmvadiscr          ( new std::vector<float>()   );
  std::auto_ptr<std::vector<float> >  againstelectronmva2rawdiscr      ( new std::vector<float>()   );
  std::auto_ptr<std::vector<float> >  againstelectronmva2categorydiscr ( new std::vector<float>()   );
  std::auto_ptr<std::vector<float> >  againstelectronvloosemva2discr   ( new std::vector<float>()   );
  std::auto_ptr<std::vector<float> >  againstelectronloosemva2discr    ( new std::vector<float>()   );
  std::auto_ptr<std::vector<float> >  againstelectronmediummva2discr   ( new std::vector<float>()   );
  std::auto_ptr<std::vector<float> >  againstelectrontightmva2discr    ( new std::vector<float>()   );
  std::auto_ptr<std::vector<float> >  againstelectronmva3rawdiscr      ( new std::vector<float>()   );
  std::auto_ptr<std::vector<float> >  againstelectronmva3categorydiscr ( new std::vector<float>()   );
  std::auto_ptr<std::vector<float> >  againstelectronloosemva3discr    ( new std::vector<float>()   );
  std::auto_ptr<std::vector<float> >  againstelectronmediummva3discr   ( new std::vector<float>()   );
  std::auto_ptr<std::vector<float> >  againstelectrontightmva3discr    ( new std::vector<float>()   );
  std::auto_ptr<std::vector<float> >  againstelectronvtightmva3discr   ( new std::vector<float>()   );
  std::auto_ptr<std::vector<float> >  againstelectrondeadecaldiscr     ( new std::vector<float>()   );
  //
  std::auto_ptr<std::vector<float> >  againstmuonloosediscr   ( new std::vector<float>()   );
  std::auto_ptr<std::vector<float> >  againstmuonmediumdiscr  ( new std::vector<float>()   );
  std::auto_ptr<std::vector<float> >  againstmuontightdiscr   ( new std::vector<float>()   );
  std::auto_ptr<std::vector<float> >  againstmuonloose2discr  ( new std::vector<float>()   );
  std::auto_ptr<std::vector<float> >  againstmuonmedium2discr ( new std::vector<float>()   );
  std::auto_ptr<std::vector<float> >  againstmuontight2discr  ( new std::vector<float>()   );
  //
  std::auto_ptr<std::vector<float> >  vlooseisolationdiscr                           ( new std::vector<float>()   );
  std::auto_ptr<std::vector<float> >  looseisolationdiscr                            ( new std::vector<float>()   );
  std::auto_ptr<std::vector<float> >  mediumisolationdiscr                           ( new std::vector<float>()   );
  std::auto_ptr<std::vector<float> >  tightisolationdiscr                            ( new std::vector<float>()   );
  std::auto_ptr<std::vector<float> >  vlooseisolationdeltabetacorrdiscr              ( new std::vector<float>()   );
  std::auto_ptr<std::vector<float> >  looseisolationdeltabetacorrdiscr               ( new std::vector<float>()   );
  std::auto_ptr<std::vector<float> >  mediumisolationdeltabetacorrdiscr              ( new std::vector<float>()   );
  std::auto_ptr<std::vector<float> >  tightisolationdeltabetacorrdiscr               ( new std::vector<float>()   );
  std::auto_ptr<std::vector<float> >  vloosecombinedisolationdeltabetacorrdiscr      ( new std::vector<float>()   );
  std::auto_ptr<std::vector<float> >  loosecombinedisolationdeltabetacorrdiscr       ( new std::vector<float>()   );
  std::auto_ptr<std::vector<float> >  mediumcombinedisolationdeltabetacorrdiscr      ( new std::vector<float>()   );
  std::auto_ptr<std::vector<float> >  tightcombinedisolationdeltabetacorrdiscr       ( new std::vector<float>()   );
  std::auto_ptr<std::vector<float> >  combinedisolationdeltabetacorr3hitsdiscr       ( new std::vector<float>()   );
  std::auto_ptr<std::vector<float> >  loosecombinedisolationdeltabetacorr3hitsdiscr  ( new std::vector<float>()   );
  std::auto_ptr<std::vector<float> >  mediumcombinedisolationdeltabetacorr3hitsdiscr ( new std::vector<float>()   );
  std::auto_ptr<std::vector<float> >  tightcombinedisolationdeltabetacorr3hitsdiscr  ( new std::vector<float>()   );
  std::auto_ptr<std::vector<float> >  isolationmvarawdiscr                           ( new std::vector<float>()   );
  std::auto_ptr<std::vector<float> >  looseisolationmvadiscr                         ( new std::vector<float>()   );
  std::auto_ptr<std::vector<float> >  mediumisolationmvadiscr                        ( new std::vector<float>()   );
  std::auto_ptr<std::vector<float> >  tightisolationmvadiscr                         ( new std::vector<float>()   );
  std::auto_ptr<std::vector<float> >  looseisolationmva2discr                        ( new std::vector<float>()   );
  std::auto_ptr<std::vector<float> >  mediumisolationmva2discr                       ( new std::vector<float>()   );
  std::auto_ptr<std::vector<float> >  tightisolationmva2discr                        ( new std::vector<float>()   );
  //
  //Signal Particles (HPSTau)  
  std::auto_ptr<std::vector<float> >  signalpfchargedhadrcandspt    ( new std::vector<float>()   );
  std::auto_ptr<std::vector<float> >  signalpfchargedhadrcandseta   ( new std::vector<float>()   );
  std::auto_ptr<std::vector<float> >  signalpfchargedhadrcandsphi   ( new std::vector<float>()   );
  std::auto_ptr<std::vector<float> >  signalpfchargedhadrcandscount ( new std::vector<float>()   );
  std::auto_ptr<std::vector<float> >  signalpfneutrhadrcandspt      ( new std::vector<float>()   );
  std::auto_ptr<std::vector<float> >  signalpfneutrhadrcandseta     ( new std::vector<float>()   );
  std::auto_ptr<std::vector<float> >  signalpfneutrhadrcandsphi     ( new std::vector<float>()   );
  std::auto_ptr<std::vector<float> >  signalpfneutrhadrcandscount   ( new std::vector<float>()   );
  std::auto_ptr<std::vector<float> >  signalpfgammacandspt          ( new std::vector<float>()   );
  std::auto_ptr<std::vector<float> >  signalpfgammacandseta         ( new std::vector<float>()   );
  std::auto_ptr<std::vector<float> >  signalpfgammacandsphi         ( new std::vector<float>()   );
  std::auto_ptr<std::vector<float> >  signalpfgammacandscount       ( new std::vector<float>()   );
  //                                                                                                                                                                          
  // HPSTau Vertex Info 
  std::auto_ptr<std::vector<int> >     vtxIndex   ( new std::vector<int>   ()  );
  std::auto_ptr<std::vector<float> >  vtxDistXY  ( new std::vector<float>()  );
  std::auto_ptr<std::vector<float> >  vtxDistZ   ( new std::vector<float>()  );
  std::auto_ptr<std::vector<float> >  vtx0DistXY ( new std::vector<float>()  );
  std::auto_ptr<std::vector<float> >  vtx0DistZ  ( new std::vector<float>()  );
  //
  // --------------------------------------------------------------------------------------- //
  // HPS Tau Optional Isolation information
  //std::auto_ptr<std::vector<float> >  isolationpfchargedhadrcandspt    ( new std::vector<float>()   );
  //std::auto_ptr<std::vector<float> >  isolationpfchargedhadrcandseta   ( new std::vector<float>()   );
  //std::auto_ptr<std::vector<float> >  isolationpfchargedhadrcandsphi   ( new std::vector<float>()   );
  //std::auto_ptr<std::vector<float> >  isolationpfchargedhadrcandscount ( new std::vector<float>()   );
  //std::auto_ptr<std::vector<float> >  isolationpfneutrhadrcandspt      ( new std::vector<float>()   );
  //std::auto_ptr<std::vector<float> >  isolationpfneutrhadrcandseta     ( new std::vector<float>()   );
  //std::auto_ptr<std::vector<float> >  isolationpfneutrhadrcandsphi     ( new std::vector<float>()   );
  //std::auto_ptr<std::vector<float> >  isolationpfneutrhadrcandscount   ( new std::vector<float>()   );
  //std::auto_ptr<std::vector<float> >  isolationpfgammacandspt          ( new std::vector<float>()   );
  //std::auto_ptr<std::vector<float> >  isolationpfgammacandseta         ( new std::vector<float>()   );
  //std::auto_ptr<std::vector<float> >  isolationpfgammacandsphi         ( new std::vector<float>()   );
  //std::auto_ptr<std::vector<float> >  isolationpfgammacandscount       ( new std::vector<float>()   );
  // --------------------------------------------------------------------------------------- //
  //  
  edm::Handle<reco::VertexCollection> primaryVertices;
  iEvent.getByLabel(vtxInputTag,primaryVertices);
  //
  edm::Handle<std::vector<pat::Tau> > taus;
  iEvent.getByLabel(inputTag, taus);
  //
  if(taus.isValid()) {
    edm::LogInfo("RootTupleMakerV2_TausInfo") << "Total # Taus: " << taus->size();

    std::vector<pat::Tau>::const_iterator it     = taus -> begin();
    std::vector<pat::Tau>::const_iterator it_end = taus -> end();
    //
    for (; it != it_end; ++it ) { 
      if ( eta->size() > maxSize ) break;
      //
      // SCTau Discriminators are at: (cvs up -r 1.53 PhysicsTools/PatAlgos/python/tools/tauTools.py)
      // http://cmssw.cvs.cern.ch/cgi-bin/cmssw.cgi/CMSSW/PhysicsTools/PatAlgos/python/tools/tauTools.py?revision=1.53&view=markup
      // Out of the Box, Clean SC Taus are given by switchToPFTau:
      // tauID("leadingTrackFinding") > 0.5 
      // tauID("leadingPionPtCut") > 0.5 
      // tauID("byIsolationUsingLeadingPion") > 0.5
      // tauID("againstMuon") > 0.5 
      // tauID("againstElectron") > 0.5'
      // (signalPFChargedHadrCands.size() = 1 | signalPFChargedHadrCands.size() = 3)'
      //
      if(isSCTau){
	againstelectrondiscr                -> push_back( it->tauID("againstElectron" )               );//applied in switchToPFTau
	againstmuondiscr                    -> push_back( it->tauID("againstMuon")                    );//applied in switchToPFTau 
	ecalisolationdiscr                  -> push_back( it->tauID("ecalIsolation")                  );
	ecalisolationusingleadingpiondiscr  -> push_back( it->tauID("byIsolationUsingLeadingPion")    );
	isolationdiscr                      -> push_back( it->tauID("byIsolation")                    );
	isolationusingleadingpiondiscr      -> push_back( it->tauID("byIsolationUsingLeadingPion")    );//applied in switchToPFTau 
	leadingtrackfindingdiscr            -> push_back( it->tauID("leadingTrackFinding")            );//applied in switchToPFTau
	leadingtrackptcutdiscr              -> push_back( it->tauID("leadingTrackPtCut")              );
	leadingpionptcutdiscr               -> push_back( it->tauID("leadingPionPtCut")               );//applied in switchToPFTau 
	trackisolationdiscr                 -> push_back( it->tauID("trackIsolation")                 );
	trackisolationusingleadingpiondiscr -> push_back( it->tauID("trackIsolationUsingLeadingPion") );
	//
	tancdiscr                 -> push_back( it->tauID("byTaNC")                 );
	tancfronepercentdiscr     -> push_back( it->tauID("byTaNCfrOnePercent")     );
	tancfrhalfpercentdiscr    -> push_back( it->tauID("byTaNCfrHalfPercent")    );
	tancfrquarterpercentdiscr -> push_back( it->tauID("byTaNCfrQuarterPercent") );
	tancfrtenthpercentdiscr   -> push_back( it->tauID("byTaNCfrTenthPercent")   );
      }
      //
      // --  PAT V08-09-51 -- https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuidePATReleaseNotes52X#V08_09_51  --
      //
      // HPSTau Discriminators are at: 
      // http://cmssw.cvs.cern.ch/cgi-bin/cmssw.cgi/CMSSW/PhysicsTools/PatAlgos/python/producersLayer1/tauProducer_cfi.py?revision=1.31.6.4&view=markup
      //
      // Clean HPS Taus are defined at:
      // http://cmssw.cvs.cern.ch/cgi-bin/cmssw.cgi/CMSSW/PhysicsTools/PatAlgos/python/cleaningLayer1/tauCleaner_cfi.py?revision=1.10&view=markup
      //
      // Out of the Box, Clean HPS Taus are given by cleanPatTaus: 
      // tauID("decayModeFinding") > 0.5                            ------ KEPT
      // tauID("byLooseCombinedIsolationDeltaBetaCorr") > 0.5       ------ REMOVED
      // tauID("againstMuonMedium") > 0.5                           ------ REMOVED
      // tauID("againstElectronMedium") > 0.5                       ------ REMOVED
      // pt > 20.0                                                  ------ CHANGED to pt>15.0
      // abs(eta) < 2.3                                             ------ CHANGED to abs(eta)<2.5
      //
      if(isHPSTau){
	decaymodefindingdiscr      -> push_back( it->tauID("decayModeFinding")      );//applied in cleanPatTaus, KEPT in cmssw_cfg.py
	//
	againstelectronloosediscr        -> push_back( it->tauID("againstElectronLoose")        );
	againstelectronmediumdiscr       -> push_back( it->tauID("againstElectronMedium")       );//applied in cleanPatTaus, REMOVED in cmssw_cfg.py
	againstelectrontightdiscr        -> push_back( it->tauID("againstElectronTight")        );
	againstelectronmvadiscr          -> push_back( it->tauID("againstElectronMVA")          );
	againstelectronmva2rawdiscr      -> push_back( it->tauID("againstElectronMVA2raw")      );
	againstelectronmva2categorydiscr -> push_back( it->tauID("againstElectronMVA2category") );
	againstelectronvloosemva2discr   -> push_back( it->tauID("againstElectronVLooseMVA2")   );
	againstelectronloosemva2discr    -> push_back( it->tauID("againstElectronLooseMVA2")    );
	againstelectronmediummva2discr   -> push_back( it->tauID("againstElectronMediumMVA2")   );
	againstelectrontightmva2discr    -> push_back( it->tauID("againstElectronTightMVA2")    );
	againstelectronmva3rawdiscr      -> push_back( it->tauID("againstElectronMVA3raw")      );
	againstelectronmva3categorydiscr -> push_back( it->tauID("againstElectronMVA3category") );
	againstelectronloosemva3discr    -> push_back( it->tauID("againstElectronLooseMVA3")    );
	againstelectronmediummva3discr   -> push_back( it->tauID("againstElectronMediumMVA3")   );
	againstelectrontightmva3discr    -> push_back( it->tauID("againstElectronTightMVA3")    );
	againstelectronvtightmva3discr   -> push_back( it->tauID("againstElectronVTightMVA3")   );
	againstelectrondeadecaldiscr     -> push_back( it->tauID("againstElectronDeadECAL")     );
	//
	againstmuonloosediscr   -> push_back( it->tauID("againstMuonLoose")   );
	againstmuonmediumdiscr  -> push_back( it->tauID("againstMuonMedium")  );//applied in cleanPatTaus, REMOVED in cmssw_cfg.py
	againstmuontightdiscr   -> push_back( it->tauID("againstMuonTight")   );
	againstmuonloose2discr  -> push_back( it->tauID("againstMuonLoose2")  );
	againstmuonmedium2discr -> push_back( it->tauID("againstMuonMedium2") );
	againstmuontight2discr  -> push_back( it->tauID("againstMuonTight2")  );
	//
	vlooseisolationdiscr              -> push_back( it->tauID("byVLooseIsolation") );
	looseisolationdiscr               -> push_back( it->tauID("byLooseIsolation")  );
	mediumisolationdiscr              -> push_back( it->tauID("byMediumIsolation") );
	tightisolationdiscr               -> push_back( it->tauID("byTightIsolation" ) );
	vlooseisolationdeltabetacorrdiscr -> push_back( it->tauID("byVLooseIsolationDeltaBetaCorr") );
	looseisolationdeltabetacorrdiscr  -> push_back( it->tauID("byLooseIsolationDeltaBetaCorr")  );
	mediumisolationdeltabetacorrdiscr -> push_back( it->tauID("byMediumIsolationDeltaBetaCorr") );
	tightisolationdeltabetacorrdiscr  -> push_back( it->tauID("byTightIsolationDeltaBetaCorr")  );
	vloosecombinedisolationdeltabetacorrdiscr -> push_back( it->tauID("byVLooseCombinedIsolationDeltaBetaCorr") );
	loosecombinedisolationdeltabetacorrdiscr  -> push_back( it->tauID("byLooseCombinedIsolationDeltaBetaCorr" ) );//applied in cleanPatTaus, REMOVED in cmssw_cfg.py
	mediumcombinedisolationdeltabetacorrdiscr -> push_back( it->tauID("byMediumCombinedIsolationDeltaBetaCorr") );
	tightcombinedisolationdeltabetacorrdiscr  -> push_back( it->tauID("byTightCombinedIsolationDeltaBetaCorr" ) );
	combinedisolationdeltabetacorr3hitsdiscr       -> push_back( it->tauID("byCombinedIsolationDeltaBetaCorrRaw3Hits" )   );
	loosecombinedisolationdeltabetacorr3hitsdiscr  -> push_back( it->tauID("byLooseCombinedIsolationDeltaBetaCorr3Hits")  );
	mediumcombinedisolationdeltabetacorr3hitsdiscr -> push_back( it->tauID("byMediumCombinedIsolationDeltaBetaCorr3Hits") );
	tightcombinedisolationdeltabetacorr3hitsdiscr  -> push_back( it->tauID("byTightCombinedIsolationDeltaBetaCorr3Hits")  );
	isolationmvarawdiscr     -> push_back( it->tauID("byIsolationMVAraw")     );
	looseisolationmvadiscr   -> push_back( it->tauID("byLooseIsolationMVA")   );
	mediumisolationmvadiscr  -> push_back( it->tauID("byMediumIsolationMVA")  );
	tightisolationmvadiscr   -> push_back( it->tauID("byTightIsolationMVA")   );
	looseisolationmva2discr  -> push_back( it->tauID("byLooseIsolationMVA2")  );
	mediumisolationmva2discr -> push_back( it->tauID("byMediumIsolationMVA2") );
	tightisolationmva2discr  -> push_back( it->tauID("byTightIsolationMVA2")  );
	//
      }
      //
      //
      eta -> push_back ( (float)(it -> eta()) ) ;
      phi -> push_back ( (float)(it -> phi()) ) ;
      pt  -> push_back ( (float)(it -> pt() ) ) ;
      et  -> push_back ( (float)(it -> et() ) ) ;
      charge  -> push_back ( (int)(it -> charge() ) ) ;
      if(  it ->isPFTau()  ){ispftau   -> push_back ( 1 ) ;}// all should be PF Tau
      if( !it ->isPFTau()  ){ispftau   -> push_back ( 0 ) ;}// all should be PF Tau
      if(  it ->isCaloTau()){iscalotau -> push_back ( 1 ) ;}// all should be PF Tau
      if( !it ->isCaloTau()){iscalotau -> push_back ( 0 ) ;}// all should be PF Tau
      decaymode                        -> push_back( (float)(it->decayMode())  );
      emfraction                       -> push_back( (float)(it->emFraction()) );
      hcal3x3overplead                 -> push_back( (float)(it->hcal3x3OverPLead()) );
      hcalmaxoverplead                 -> push_back( (float)(it->hcalMaxOverPLead()) );
      hcaltotoverplead                 -> push_back( (float)(it->hcalTotOverPLead()) );
      isolationpfchargedhadrcandsptsum -> push_back( (float)(it->isolationPFChargedHadrCandsPtSum()) );
      isolationpfgammacandsetsum       -> push_back( (float)(it->isolationPFGammaCandsEtSum())       );
      leadpfchargedhadrcandsignedsipt  -> push_back( (float)(it->leadPFChargedHadrCandsignedSipt())  );
      reco::PFCandidateRef leadPFChargedHadrCand_Ref = it->leadPFChargedHadrCand();
      if(leadPFChargedHadrCand_Ref.isNonnull()){// this check is needed in case hpsTau fails decayModeFinding.
	etaleadcharged                   -> push_back( (float)(leadPFChargedHadrCand_Ref->eta()) );
	phileadcharged                   -> push_back( (float)(leadPFChargedHadrCand_Ref->phi()) );
	ptleadcharged                    -> push_back( (float)(leadPFChargedHadrCand_Ref->pt())  );
      }
      phiphimoment                     -> push_back( (float)(it->phiphiMoment()) );
      etaetamoment                     -> push_back( (float)(it->etaetaMoment()) );
      etaphimoment                     -> push_back( (float)(it->etaphiMoment()) );
      ecalstripsumeoverplead           -> push_back( (float)(it->ecalStripSumEOverPLead()) );
      bremsrecoveryeoverplead          -> push_back( (float)(it->bremsRecoveryEOverPLead()) );
      maximumhcalpfclusteret           -> push_back( (float)(it->maximumHCALPFClusterEt()) );
      //
      // ----------- Vertex association ----------- //
      float minVtxDist3D = 9999.;
      int    vtxIndex_    = -1;
      float vtxDistXY_   = -9999.;
      float vtxDistZ_    = -9999.;      
      float vtx0DistXY_  = -9999.;
      float vtx0DistZ_   = -9999.;

      if( primaryVertices.isValid() ) {
	edm::LogInfo("RootTupleMakerV2_TausInfo") << "Total # Primary Vertices: " << primaryVertices->size();

	int i_vertex = 0;

        for( reco::VertexCollection::const_iterator v_it=primaryVertices->begin() ; v_it!=primaryVertices->end() ; ++v_it ) {
	  float distX  = (it->vertex()).x()-(v_it->position()).x();
	  float distY  = (it->vertex()).y()-(v_it->position()).y();
	  float distZ  = (it->vertex()).z()-(v_it->position()).z();
	  float distXY = sqrt(pow(distX,2) + pow(distY,2));
	  float dist3D = sqrt(pow(distXY,2) + pow(distZ,2));
	  
	  if ( i_vertex == 0 ) {  //leading vertex, by default sorted by sum(pt^2)
	    vtx0DistXY_ = distXY;
	    vtx0DistZ_  = distZ ;
	  }

          if( dist3D<minVtxDist3D ) {
            minVtxDist3D = dist3D;
            vtxIndex_    = int(std::distance(primaryVertices->begin(),v_it));
            vtxDistXY_   = distXY;
            vtxDistZ_    = distZ;
          }

	  i_vertex++;
        }
      } else {
	edm::LogError("RootTupleMakerV2_TausError") << "Error! Can't get the product " << vtxInputTag;
      }
      vtxIndex                 -> push_back( vtxIndex_  );
      vtxDistXY                -> push_back( vtxDistXY_ );
      vtxDistZ                 -> push_back( vtxDistZ_  );
      vtx0DistXY               -> push_back( vtx0DistXY_);
      vtx0DistZ                -> push_back( vtx0DistZ_ );
      //
      // ----------- Gen-Reco Matching ----------- //
      float genparPt =-999.; float genjetPt =-999.;
      float genparEta=-999.; float genjetEta=-999.;
      float genparPhi=-999.; float genjetPhi=-999.;
      if ( !iEvent.isRealData() ) {
	for(uint igen = 0 ; igen < it->genParticleRefs().size() ; ++igen ){//genParticleRefs().size() is either 0 or 1
	  genparPt=it->genParticle(igen)->pt();
	  genparEta=it->genParticle(igen)->eta();
	  genparPhi=it->genParticle(igen)->phi();
	}
	if( it->genJet() ){
	  genjetPt=it->genJet()->pt();
	  genjetEta=it->genJet()->eta();
	  genjetPhi=it->genJet()->phi();
	}
      }
      matchedgenparticlept   -> push_back ( (float)(genparPt)  );
      matchedgenparticleeta  -> push_back ( (float)(genparEta) );
      matchedgenparticlephi  -> push_back ( (float)(genparPhi) );
      matchedgenjetpt        -> push_back ( (float)(genjetPt)  );
      matchedgenjeteta       -> push_back ( (float)(genjetEta) );
      matchedgenjetphi       -> push_back ( (float)(genjetPhi) );
      //
      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      //  --  User Isolation and isoDeposit Methods -- Feb 2012
      // 
      // User Isolation definitions and "keys" are at:
      //  http://cmssw.cvs.cern.ch/cgi-bin/cmssw.cgi/CMSSW/DataFormats/PatCandidates/interface/Isolation.h?revision=1.7&view=markup
      //  http://cmssw.cvs.cern.ch/cgi-bin/cmssw.cgi/CMSSW/PhysicsTools/PatAlgos/python/producersLayer1/tauProducer_cfi.py?revision=1.27&view=markup
      //   it -> userIsolation(pat::PfChargedHadronIso);
      //   it -> userIsolation(pat::PfNeutralHadronIso);
      //   it -> userIsolation(pat::PfGammaIso);
      //
      // IsoDeposit is defined at:
      //  http://cmssdt.cern.ch/SDT/doxygen/CMSSW_4_2_5/doc/html/d4/d0a/classreco_1_1IsoDeposit.html
      //   (it -> isoDeposit(pat::PfChargedHadronIso)) -> candEnergy();
      //   (it -> isoDeposit(pat::PfChargedHadronIso)) -> depositWithin(0.5);
      //
      //
      // In comparison:
      // "pat::tau -> isolationPFCands()" is taken directly from the reco::PFTau object.
      // The parameters (cone size, pt threshold) are defined at:
      // http://cmssw.cvs.cern.ch/cgi-bin/cmssw.cgi/CMSSW/RecoTauTag/RecoTau/python/PFRecoTauQualityCuts_cfi.py?revision=1.4&view=markup
      // http://cmssw.cvs.cern.ch/cgi-bin/cmssw.cgi/CMSSW/RecoTauTag/RecoTau/python/PFRecoTauProducer_cfi.py?revision=1.23&view=markup
      // http://cmssw.cvs.cern.ch/cgi-bin/cmssw.cgi/CMSSW/PhysicsTools/IsolationAlgos/plugins/PFTauExtractor.cc?revision=1.1&view=markup
      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      //
      if(isHPSTau){
	reco::PFCandidateRefVector signalPFChargedHadrCands_RefVector = it->signalPFChargedHadrCands();
	float signalPFChargedHadrCands_RefVector_size=0;
	if(signalPFChargedHadrCands_RefVector.isNonnull()){
	  signalPFChargedHadrCands_RefVector_size=(float)(signalPFChargedHadrCands_RefVector.size());
	  reco::PFCandidateRefVector::iterator itSignalChargedHad     = signalPFChargedHadrCands_RefVector.begin();
	  reco::PFCandidateRefVector::iterator itSignalChargedHad_end = signalPFChargedHadrCands_RefVector.end();
	  for (; itSignalChargedHad != itSignalChargedHad_end; ++itSignalChargedHad ) {
	    signalpfchargedhadrcandspt  -> push_back((float)(*itSignalChargedHad)->pt());
	    signalpfchargedhadrcandseta -> push_back((float)(*itSignalChargedHad)->eta());
	    signalpfchargedhadrcandsphi -> push_back((float)(*itSignalChargedHad)->phi());
	  }
	}
	signalpfchargedhadrcandscount -> push_back(signalPFChargedHadrCands_RefVector_size);
	//
	reco::PFCandidateRefVector signalPFNeutrHadrCands_RefVector = it->signalPFNeutrHadrCands();
	float signalPFNeutrHadrCands_RefVector_size=0;
	if(signalPFNeutrHadrCands_RefVector.isNonnull()){
  	  signalPFNeutrHadrCands_RefVector_size=(float)(signalPFNeutrHadrCands_RefVector.size());
	  reco::PFCandidateRefVector::iterator itSignalNeutrHadr     = signalPFNeutrHadrCands_RefVector.begin();
	  reco::PFCandidateRefVector::iterator itSignalNeutrHadr_end = signalPFNeutrHadrCands_RefVector.end();
	  for (; itSignalNeutrHadr != itSignalNeutrHadr_end; ++itSignalNeutrHadr ) {
	    signalpfneutrhadrcandspt  -> push_back((float)(*itSignalNeutrHadr)->pt());	
	    signalpfneutrhadrcandseta -> push_back((float)(*itSignalNeutrHadr)->eta());	
	    signalpfneutrhadrcandsphi -> push_back((float)(*itSignalNeutrHadr)->phi());	
	  }
	}
	signalpfneutrhadrcandscount -> push_back(signalPFNeutrHadrCands_RefVector_size);
	//
	reco::PFCandidateRefVector signalPFGammaCands_RefVector = it->signalPFGammaCands();
	float signalPFGammaCands_RefVector_size=0;
	if(signalPFGammaCands_RefVector.isNonnull()){
	  signalPFGammaCands_RefVector_size=(float)(signalPFGammaCands_RefVector.size());
	  reco::PFCandidateRefVector::iterator itSignalGamma     = signalPFGammaCands_RefVector.begin();
	  reco::PFCandidateRefVector::iterator itSignalGamma_end = signalPFGammaCands_RefVector.end();
	  for (; itSignalGamma != itSignalGamma_end; ++itSignalGamma ) {
	    signalpfgammacandspt  -> push_back((float)(*itSignalGamma)->pt());
	    signalpfgammacandseta -> push_back((float)(*itSignalGamma)->eta());
	    signalpfgammacandsphi -> push_back((float)(*itSignalGamma)->phi());
	  }
	}
	signalpfgammacandscount -> push_back(signalPFGammaCands_RefVector_size);
	//	
	// --------------------------------------------------------------------------------------- //
	// HPS Tau Optional Isolation information 
	//reco::PFCandidateRefVector isoPFChargedHadrCands_RefVector = it->isolationPFChargedHadrCands();
	//isolationpfchargedhadrcandscount ->push_back((float)(isoPFChargedHadrCands_RefVector.size()));
	//if(isoPFChargedHadrCands_RefVector.isNonnull()){
	//  reco::PFCandidateRefVector::iterator itIsoChargedHad     = isoPFChargedHadrCands_RefVector.begin();
	//  reco::PFCandidateRefVector::iterator itIsoChargedHad_end = isoPFChargedHadrCands_RefVector.end();
	//  for (; itIsoChargedHad != itIsoChargedHad_end; ++itIsoChargedHad ) {
	//    isolationpfchargedhadrcandspt  -> push_back((float)(*itIsoChargedHad)->pt());
	//    isolationpfchargedhadrcandseta -> push_back((float)(*itIsoChargedHad)->eta());
	//    isolationpfchargedhadrcandsphi -> push_back((float)(*itIsoChargedHad)->phi());
	//  }
	//}
	//reco::PFCandidateRefVector isoPFNeutrHadrCands_RefVector = it->isolationPFNeutrHadrCands();
	//isolationpfneutrhadrcandscount -> push_back((float)(isoPFNeutrHadrCands_RefVector.size()));
	//if(isoPFNeutrHadrCands_RefVector.isNonnull()){
	//  reco::PFCandidateRefVector::iterator itIsoNeutrHadr     = isoPFNeutrHadrCands_RefVector.begin();
	//  reco::PFCandidateRefVector::iterator itIsoNeutrHadr_end = isoPFNeutrHadrCands_RefVector.end();
	//  for (; itIsoNeutrHadr != itIsoNeutrHadr_end; ++itIsoNeutrHadr ) {
	//    isolationpfneutrhadrcandspt  -> push_back((float)(*itIsoNeutrHadr)->pt());
	//    isolationpfneutrhadrcandseta -> push_back((float)(*itIsoNeutrHadr)->eta());
	//    isolationpfneutrhadrcandsphi -> push_back((float)(*itIsoNeutrHadr)->phi());
	//  }
	//}
	//reco::PFCandidateRefVector isoPFGammaCands_RefVector = it->isolationPFGammaCands();
	//isolationpfgammacandscount -> push_back((float)(isoPFGammaCands_RefVector.size()));
	//if(isoPFGammaCands_RefVector.isNonnull()){
	//  reco::PFCandidateRefVector::iterator itIsoGamma     = isoPFGammaCands_RefVector.begin();
	//  reco::PFCandidateRefVector::iterator itIsoGamma_end = isoPFGammaCands_RefVector.end();
	//  for (; itIsoGamma != itIsoGamma_end; ++itIsoGamma ) {
	//    isolationpfgammacandspt  -> push_back((float)(*itIsoGamma)->pt());
	//    isolationpfgammacandseta -> push_back((float)(*itIsoGamma)->eta());
	//    isolationpfgammacandsphi -> push_back((float)(*itIsoGamma)->phi());
	//  }
	//}
	// --------------------------------------------------------------------------------------- //
	//
      }
      //
    }
  } else {
    edm::LogError("RootTupleMakerV2_TausError") << "Error! Can't get the product " << inputTag;
  }
 /* iEvent.put( eta,                              prefix + "Eta"                               + suffix );
  iEvent.put( phi,                              prefix + "Phi"                               + suffix );
  iEvent.put( pt,                               prefix + "Pt"                                + suffix );
  iEvent.put( et,                               prefix + "Et"                                + suffix );
  iEvent.put( charge,                           prefix + "Charge"                            + suffix );
  iEvent.put( ispftau,                          prefix + "IsPFTau"                           + suffix );
  iEvent.put( iscalotau,                        prefix + "IsCaloTau"                         + suffix );
  iEvent.put( decaymode,                        prefix + "DecayMode"                         + suffix );
  iEvent.put( emfraction,                       prefix + "EmFraction"                        + suffix );
  iEvent.put( hcal3x3overplead,                 prefix + "Hcal3x3OverPLead"                  + suffix );
  iEvent.put( hcalmaxoverplead,                 prefix + "HcalMaxOverPLead"                  + suffix );
  iEvent.put( hcaltotoverplead,                 prefix + "HcalTotOverPLead"                  + suffix );
  iEvent.put( isolationpfchargedhadrcandsptsum, prefix + "IsolationPFChargedHadrCandsPtSum"  + suffix );
  iEvent.put( isolationpfgammacandsetsum,       prefix + "IsolationPFGammaCandsEtSum"        + suffix );
  iEvent.put( leadpfchargedhadrcandsignedsipt,  prefix + "LeadPFChargedHadrCandsignedSipt"   + suffix );
  iEvent.put( etaleadcharged,                   prefix + "EtaLeadCharged"                    + suffix );
  iEvent.put( phileadcharged,                   prefix + "PhiLeadCharged"                    + suffix );
  iEvent.put( ptleadcharged,                    prefix + "PtLeadCharged"                     + suffix );
  iEvent.put( phiphimoment,                     prefix + "PhiphiMoment"                      + suffix );
  iEvent.put( etaetamoment,                     prefix + "EtaetaMoment"                      + suffix );
  iEvent.put( etaphimoment,                     prefix + "EtaphiMoment"                      + suffix );
  iEvent.put( ecalstripsumeoverplead,           prefix + "EcalStripSumEOverPLead"            + suffix );
  iEvent.put( bremsrecoveryeoverplead,          prefix + "BremsRecoveryEOverPLead"           + suffix );
  iEvent.put( maximumhcalpfclusteret,           prefix + "MaximumHCALPFClusterEt"            + suffix );
  iEvent.put( matchedgenparticlept,             prefix + "MatchedGenParticlePt"              + suffix );
  iEvent.put( matchedgenparticleeta,            prefix + "MatchedGenParticleEta"             + suffix );
  iEvent.put( matchedgenparticlephi,            prefix + "MatchedGenParticlePhi"             + suffix );
  iEvent.put( matchedgenjetpt,                  prefix + "MatchedGenJetPt"                   + suffix );
  iEvent.put( matchedgenjeteta,                 prefix + "MatchedGenJetEta"                  + suffix );
  iEvent.put( matchedgenjetphi,                 prefix + "MatchedGenJetPhi"                  + suffix );
*/
  //
  if(isSCTau){
/*    iEvent.put( leadingtrackfindingdiscr,               prefix + "LeadingTrackFindingDiscr"            + suffix );
    iEvent.put( leadingtrackptcutdiscr,                 prefix + "LeadingTrackPtCutDiscr"              + suffix );
    iEvent.put( leadingpionptcutdiscr,                  prefix + "LeadingPionPtCutDiscr"               + suffix );
    iEvent.put( isolationdiscr,                         prefix + "IsolationDiscr"                      + suffix );
    iEvent.put( trackisolationdiscr,                    prefix + "TrackIsolationDiscr"                 + suffix );
    iEvent.put( ecalisolationdiscr,                     prefix + "EcalIsolationDiscr"                  + suffix );
    iEvent.put( isolationusingleadingpiondiscr,         prefix + "IsolationUsingLeadingPionDiscr"      + suffix );
    iEvent.put( trackisolationusingleadingpiondiscr,    prefix + "TrackIsolationUsingLeadingPionDiscr" + suffix );
    iEvent.put( ecalisolationusingleadingpiondiscr,     prefix + "EcalIsolationUsingLeadingPionDiscr"  + suffix );
    iEvent.put( againstelectrondiscr,                   prefix + "AgainstElectronDiscr"                + suffix );
    iEvent.put( againstmuondiscr,                       prefix + "AgainstMuonDiscr"                    + suffix );
    iEvent.put( tancdiscr,                              prefix + "TaNCDiscr"                           + suffix );
    iEvent.put( tancfronepercentdiscr,                  prefix + "TaNCfrOnePercentDiscr"               + suffix );
    iEvent.put( tancfrhalfpercentdiscr,                 prefix + "TaNCfrHalfPercentDiscr"              + suffix );
    iEvent.put( tancfrquarterpercentdiscr,              prefix + "TaNCfrQuarterPercentDiscr"           + suffix );
*/    iEvent.put( tancfrtenthpercentdiscr,                prefix + "TaNCfrTenthPercentDiscr"             + suffix );
 
 }
  //
  if(isHPSTau){
/*    iEvent.put( decaymodefindingdiscr,                     prefix + "DecayModeFindingDiscr"                       + suffix );
    //
    iEvent.put( againstelectronloosediscr,                 prefix + "AgainstElectronLooseDiscr"                   + suffix );
    iEvent.put( againstelectronmediumdiscr,                prefix + "AgainstElectronMediumDiscr"                  + suffix );
    iEvent.put( againstelectrontightdiscr,                 prefix + "AgainstElectronTightDiscr"                   + suffix );
    iEvent.put( againstelectronmvadiscr,                   prefix + "AgainstElectronMVADiscr"                     + suffix );
    iEvent.put( againstelectronmva2rawdiscr,               prefix + "AgainstElectronMVA2rawDiscr"                 + suffix );
    iEvent.put( againstelectronmva2categorydiscr,          prefix + "AgainstElectronMVA2categoryDiscr"            + suffix );
    iEvent.put( againstelectronvloosemva2discr,            prefix + "AgainstElectronVLooseMVA2Discr"              + suffix );
    iEvent.put( againstelectronloosemva2discr,             prefix + "AgainstElectronLooseMVA2Discr"               + suffix );
    iEvent.put( againstelectronmediummva2discr,            prefix + "AgainstElectronMediumMVA2Discr"              + suffix );
    iEvent.put( againstelectrontightmva2discr,             prefix + "AgainstElectronTightMVA2Discr"               + suffix );
    iEvent.put( againstelectronmva3rawdiscr,               prefix + "AgainstElectronMVA3rawDiscr"                 + suffix );
    iEvent.put( againstelectronmva3categorydiscr,          prefix + "AgainstElectronMVA3categoryDiscr"            + suffix );
    iEvent.put( againstelectronloosemva3discr,             prefix + "AgainstElectronLooseMVA3Discr"               + suffix );
    iEvent.put( againstelectronmediummva3discr,            prefix + "AgainstElectronMediumMVA3Discr"              + suffix );
    iEvent.put( againstelectrontightmva3discr,             prefix + "AgainstElectronTightMVA3Discr"               + suffix );
    iEvent.put( againstelectronvtightmva3discr,            prefix + "AgainstElectronVTightMVA3Discr"              + suffix );
    iEvent.put( againstelectrondeadecaldiscr,              prefix + "AgainstElectronDeadECALDiscr"                + suffix );
    //
    iEvent.put( againstmuonloosediscr,                     prefix + "AgainstMuonLooseDiscr"                       + suffix );
    iEvent.put( againstmuonmediumdiscr,                    prefix + "AgainstMuonMediumDiscr"                      + suffix );
    iEvent.put( againstmuontightdiscr,                     prefix + "AgainstMuonTightDiscr"                       + suffix );
    iEvent.put( againstmuonloose2discr,                    prefix + "AgainstMuonLoose2Discr"                      + suffix );
    iEvent.put( againstmuonmedium2discr,                   prefix + "AgainstMuonMedium2Discr"                     + suffix );
    iEvent.put( againstmuontight2discr,                    prefix + "AgainstMuonTight2Discr"                      + suffix );

    iEvent.put( vlooseisolationdiscr,                           prefix + "VLooseIsolationDiscr"                           + suffix );
    iEvent.put( looseisolationdiscr,                            prefix + "LooseIsolationDiscr"                            + suffix );
    iEvent.put( mediumisolationdiscr,                           prefix + "MediumIsolationDiscr"                           + suffix );
    iEvent.put( tightisolationdiscr,                            prefix + "TightIsolationDiscr"                            + suffix );
    iEvent.put( vlooseisolationdeltabetacorrdiscr,              prefix + "VLooseIsolationDeltaBetaCorrDiscr"              + suffix );
    iEvent.put( looseisolationdeltabetacorrdiscr,               prefix + "LooseIsolationDeltaBetaCorrDiscr"               + suffix );
    iEvent.put( mediumisolationdeltabetacorrdiscr,              prefix + "MediumIsolationDeltaBetaCorrDiscr"              + suffix );
    iEvent.put( tightisolationdeltabetacorrdiscr,               prefix + "TightIsolationDeltaBetaCorrDiscr"               + suffix );
    iEvent.put( vloosecombinedisolationdeltabetacorrdiscr,      prefix + "VLooseCombinedIsolationDeltaBetaCorrDiscr"      + suffix );
    iEvent.put( loosecombinedisolationdeltabetacorrdiscr,       prefix + "LooseCombinedIsolationDeltaBetaCorrDiscr"       + suffix );
    iEvent.put( mediumcombinedisolationdeltabetacorrdiscr,      prefix + "MediumCombinedIsolationDeltaBetaCorrDiscr"      + suffix );
    iEvent.put( tightcombinedisolationdeltabetacorrdiscr,       prefix + "TightCombinedIsolationDeltaBetaCorrDiscr"       + suffix );
    iEvent.put( combinedisolationdeltabetacorr3hitsdiscr,       prefix + "CombinedIsolationDeltaBetaCorr3HitsDiscr"       + suffix );
    iEvent.put( loosecombinedisolationdeltabetacorr3hitsdiscr,  prefix + "LooseCombinedIsolationDeltaBetaCorr3HitsDiscr"  + suffix );
    iEvent.put( mediumcombinedisolationdeltabetacorr3hitsdiscr, prefix + "MediumCombinedIsolationDeltaBetaCorr3HitsDiscr" + suffix );
    iEvent.put( tightcombinedisolationdeltabetacorr3hitsdiscr,  prefix + "TightCombinedIsolationDeltaBetaCorr3HitsDiscr"  + suffix );
    iEvent.put( isolationmvarawdiscr,                           prefix + "IsolationMVArawDiscr"                           + suffix );
    iEvent.put( looseisolationmvadiscr,                         prefix + "LooseIsolationMVADiscr"                         + suffix );
    iEvent.put( mediumisolationmvadiscr,                        prefix + "MediumIsolationMVADiscr"                        + suffix );
    iEvent.put( tightisolationmvadiscr,                         prefix + "TightIsolationMVADiscr"                         + suffix );
    iEvent.put( looseisolationmva2discr,                        prefix + "LooseIsolationMVA2Discr"                        + suffix );
    iEvent.put( mediumisolationmva2discr,                       prefix + "MediumIsolationMVA2Discr"                       + suffix );
    iEvent.put( tightisolationmva2discr,                        prefix + "TightIsolationMVA2Discr"                        + suffix );
    //
    iEvent.put( vtxIndex,    prefix + "VtxIndex"       + suffix );
    iEvent.put( vtxDistXY,   prefix + "VtxDistXY"      + suffix );
    iEvent.put( vtxDistZ,    prefix + "VtxDistZ"       + suffix );
    iEvent.put( vtx0DistXY,  prefix + "LeadVtxDistXY"  + suffix );
    iEvent.put( vtx0DistZ,   prefix + "LeadVtxDistZ"   + suffix );
    //
    iEvent.put( signalpfchargedhadrcandspt,   prefix +  "SignalPFChargedHadrCandsPt"   + suffix );
    iEvent.put( signalpfchargedhadrcandseta,  prefix +  "SignalPFChargedHadrCandsEta"  + suffix );
    iEvent.put( signalpfchargedhadrcandsphi,  prefix +  "SignalPFChargedHadrCandsPhi"  + suffix );
    iEvent.put( signalpfchargedhadrcandscount,prefix +  "SignalPFChargedHadrCandsCount"+ suffix );
    iEvent.put( signalpfneutrhadrcandspt,     prefix +  "SignalPFNeutrHadrCandsPt"     + suffix );
    iEvent.put( signalpfneutrhadrcandseta,    prefix +  "SignalPFNeutrHadrCandsEta"    + suffix );
    iEvent.put( signalpfneutrhadrcandsphi,    prefix +  "SignalPFNeutrHadrCandsPhi"    + suffix );
    iEvent.put( signalpfneutrhadrcandscount,  prefix +  "SignalPFNeutrHadrCandsCount"  + suffix );
    iEvent.put( signalpfgammacandspt,         prefix +  "SignalPFGammaCandsPt"         + suffix );
    iEvent.put( signalpfgammacandseta,        prefix +  "SignalPFGammaCandsEta"        + suffix );
    iEvent.put( signalpfgammacandsphi,        prefix +  "SignalPFGammaCandsPhi"        + suffix );
    iEvent.put( signalpfgammacandscount,      prefix +  "SignalPFGammaCandsCount"      + suffix );
    //
    // --------------------------------------------------------------------------------------- //
    // HPS Tau Optional Isolation information 
    //iEvent.put( isolationpfchargedhadrcandspt,   prefix +  "IsolationPFChargedHadrCandsPt"   + suffix );
    //iEvent.put( isolationpfchargedhadrcandseta,  prefix +  "IsolationPFChargedHadrCandsEta"  + suffix );
    //iEvent.put( isolationpfchargedhadrcandsphi,  prefix +  "IsolationPFChargedHadrCandsPhi"  + suffix );
    //iEvent.put( isolationpfchargedhadrcandscount,prefix +  "IsolationPFChargedHadrCandsCount"+ suffix );
    //iEvent.put( isolationpfneutrhadrcandspt,     prefix +  "IsolationPFNeutrHadrCandsPt"     + suffix );
    //iEvent.put( isolationpfneutrhadrcandseta,    prefix +  "IsolationPFNeutrHadrCandsEta"    + suffix );
    //iEvent.put( isolationpfneutrhadrcandsphi,    prefix +  "IsolationPFNeutrHadrCandsPhi"    + suffix );
    //iEvent.put( isolationpfneutrhadrcandscount,  prefix +  "IsolationPFNeutrHadrCandsCount"  + suffix );
    //iEvent.put( isolationpfgammacandspt,         prefix +  "IsolationPFGammaCandsPt"         + suffix );
    //iEvent.put( isolationpfgammacandseta,        prefix +  "IsolationPFGammaCandsEta"        + suffix );
    //iEvent.put( isolationpfgammacandsphi,        prefix +  "IsolationPFGammaCandsPhi"        + suffix );
    //iEvent.put( isolationpfgammacandscount,      prefix +  "IsolationPFGammaCandsCount"      + suffix );
*/  // --------------------------------------------------------------------------------------- //
  }
}
