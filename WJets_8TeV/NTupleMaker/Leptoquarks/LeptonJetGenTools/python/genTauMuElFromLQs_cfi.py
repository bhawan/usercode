import FWCore.ParameterSet.Config as cms

#LQ3ToTTau
genTausFromLQTaus = cms.EDProducer('GenMotherParticleSkimmer',
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

genMuonsFromLQTaus = cms.EDProducer('GenMotherParticleSkimmer',
  InputTag            = cms.InputTag('genParticles'),
  MotherPDGIds        = cms.vint32( 42, 15 ),
  MotherPDGIdsVetoed  = cms.vint32( 6, 111, 213 ),
  DaughterPDGIds      = cms.vint32( 13 ),    
  DaughterPDGStatuses = cms.vint32( 1 ),
  StoreFinalStateOnly = cms.bool ( False )                                 
)
 
genMuonsFromLQTops = cms.EDProducer('GenMotherParticleSkimmer',
  InputTag            = cms.InputTag('genParticles'),
  MotherPDGIds        = cms.vint32( 42, 6, 24 ),
  MotherPDGIdsVetoed  = cms.vint32( 111, 213 ),                                    
  DaughterPDGIds      = cms.vint32( 13  ),    
  DaughterPDGStatuses = cms.vint32( 1 ),
  StoreFinalStateOnly = cms.bool ( False )                                 
)

genElectronsFromLQTaus = cms.EDProducer('GenMotherParticleSkimmer',
  InputTag            = cms.InputTag('genParticles'),
  MotherPDGIds        = cms.vint32( 42, 15 ),
  MotherPDGIdsVetoed  = cms.vint32( 6, 111, 213 ),                                        
  DaughterPDGIds      = cms.vint32( 11 ),    
  DaughterPDGStatuses = cms.vint32( 1 ),
  StoreFinalStateOnly = cms.bool ( False )                                 
)

genElectronsFromLQTops = cms.EDProducer('GenMotherParticleSkimmer',
  InputTag            = cms.InputTag('genParticles'),
  MotherPDGIds        = cms.vint32( 42, 6, 24 ),
  MotherPDGIdsVetoed  = cms.vint32( 111, 213 ),                                        
  DaughterPDGIds      = cms.vint32( 11  ),    
  DaughterPDGStatuses = cms.vint32( 1 ),
  StoreFinalStateOnly = cms.bool ( False )                                 
)
