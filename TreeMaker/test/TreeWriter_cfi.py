import FWCore.ParameterSet.Config as cms

treeWriter = cms.EDAnalyzer("TreeMakerModule",
  verbosity = cms.untracked.int32(1),
  createTree = cms.bool(False)
)

