import FWCore.ParameterSet.Config as cms

cosmicIDproducer = cms.EDProducer('CosmicID',
                                  src = cms.InputTag("cosmicsVeto"),
                                  result = cms.string("cosmicCompatibility")
                                  #result = cms.string("timeCompatibility")
                                  #result = cms.string("backToBackCompatibility")
                                  #result = cms.string("overlapCompatibility")
                                  )
