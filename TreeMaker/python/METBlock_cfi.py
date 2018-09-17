import FWCore.ParameterSet.Config as cms

metBlock = cms.EDAnalyzer('METBlock',
  verbosity = cms.untracked.int32(0),
  metSrc = cms.InputTag('slimmedMETs'),
  corrmetSrc = cms.InputTag('slimmedMETsNoHF'), #noHF met temporarily
  puppimetSrc = cms.InputTag('slimmedMETsPuppi')
)
