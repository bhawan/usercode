import FWCore.ParameterSet.Config as cms

#DYJetsToLL L=e/mu/tau 
genTausFromZs= cms.EDProducer('GenMotherParticleSkimmer',
  InputTag            = cms.InputTag('genParticles'),
  MotherPDGIds        = cms.vint32( 23 ),
  MotherPDGIdsVetoed  = cms.vint32( ), 
  DaughterPDGIds      = cms.vint32( 15  ),    
  DaughterPDGStatuses = cms.vint32( 2 ),
  StoreFinalStateOnly = cms.bool ( True )                                 
)

genMuonsFromZs= cms.EDProducer('GenMotherParticleSkimmer',
  InputTag            = cms.InputTag('genParticles'),
  MotherPDGIds        = cms.vint32( 23 ),
  MotherPDGIdsVetoed  = cms.vint32( 111, 213 ), 
  DaughterPDGIds      = cms.vint32( 13  ),    
  DaughterPDGStatuses = cms.vint32( 1 ),
  StoreFinalStateOnly = cms.bool ( False )                                 
)

genElectronsFromZs= cms.EDProducer('GenMotherParticleSkimmer',
  InputTag            = cms.InputTag('genParticles'),
  MotherPDGIds        = cms.vint32( 23 ),
  MotherPDGIdsVetoed  = cms.vint32( 111, 213  ), 
  DaughterPDGIds      = cms.vint32( 11  ),
  DaughterPDGStatuses = cms.vint32( 1 ),
  StoreFinalStateOnly = cms.bool ( False )                                 
)
