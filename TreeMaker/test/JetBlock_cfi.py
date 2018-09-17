import FWCore.ParameterSet.Config as cms

jetBlock = cms.EDAnalyzer("JetBlock",
  verbosity = cms.untracked.int32(0),
  jetSrc = cms.InputTag('skimmedJetswqg')
)
