import FWCore.ParameterSet.Config as cms

genJetBlock = cms.EDAnalyzer('GenJetBlock',
  verbosity = cms.untracked.int32(0),
  genJetSrc = cms.InputTag('slimmedGenJets')
)
