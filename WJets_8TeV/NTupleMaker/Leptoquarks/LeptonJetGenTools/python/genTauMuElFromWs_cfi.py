import FWCore.ParameterSet.Config as cms

#WToLNu: L=e/mu/tau
genTausFromWs = cms.EDProducer('GenMotherParticleSkimmer',
  InputTag            = cms.InputTag('genParticles'),
  MotherPDGIds        = cms.vint32( 24 ),
  MotherPDGIdsVetoed  = cms.vint32( ), 
  DaughterPDGIds      = cms.vint32( 15  ),                            
  DaughterPDGStatuses = cms.vint32( 2 ),
  StoreFinalStateOnly = cms.bool ( True )                                 
)
 
genMuonsFromWs = cms.EDProducer('GenMotherParticleSkimmer',
  InputTag            = cms.InputTag('genParticles'),
  MotherPDGIds        = cms.vint32( 24 ),
  MotherPDGIdsVetoed  = cms.vint32( 111, 213 ),                                    
  DaughterPDGIds      = cms.vint32( 13  ),    
  DaughterPDGStatuses = cms.vint32( 1 ),
  StoreFinalStateOnly = cms.bool ( False )                                 
)

genElectronsFromWs = cms.EDProducer('GenMotherParticleSkimmer',
  InputTag            = cms.InputTag('genParticles'),
  MotherPDGIds        = cms.vint32( 24 ),
  MotherPDGIdsVetoed  = cms.vint32( 111, 213 ),                                        
  DaughterPDGIds      = cms.vint32( 11  ),    
  DaughterPDGStatuses = cms.vint32( 1 ),
  StoreFinalStateOnly = cms.bool ( False )                                 
)
