import FWCore.ParameterSet.Config as cms

genTausFromLQs = cms.EDProducer('GenMotherParticleSkimmer',
  InputTag            = cms.InputTag('genParticles'),
  MotherPDGIds        = cms.vint32( 42 ),
  MotherPDGIdsVetoed  = cms.vint32( 6, 24 ),                               
  DaughterPDGIds      = cms.vint32( 15  ),                            
  DaughterPDGStatuses = cms.vint32( 2 ),
  StoreFinalStateOnly = cms.bool ( True )                                 
)

genTausFromLQTops = cms.EDProducer('GenMotherParticleSkimmer',
  InputTag            = cms.InputTag('genParticles'),
  MotherPDGIds        = cms.vint32( 42, 6, 24 ),
  MotherPDGIdsVetoed  = cms.vint32( ), 
  DaughterPDGIds      = cms.vint32( 15  ),                            
  DaughterPDGStatuses = cms.vint32( 2 ),
  StoreFinalStateOnly = cms.bool ( True )                                 
)

