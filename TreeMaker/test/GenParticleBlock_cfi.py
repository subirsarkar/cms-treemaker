import FWCore.ParameterSet.Config as cms

genParticleBlock = cms.EDAnalyzer("GenParticleBlock",
  verbosity = cms.untracked.int32(0),
  genParticleSrc = cms.InputTag('prunedGenParticles')
)
