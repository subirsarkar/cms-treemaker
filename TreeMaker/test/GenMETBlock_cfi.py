import FWCore.ParameterSet.Config as cms

genMETBlock = cms.EDAnalyzer('GenMETBlock',
  verbosity = cms.untracked.int32(0),
  genMETSrc = cms.InputTag('genMetTrue')
)
