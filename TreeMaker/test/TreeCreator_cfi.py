import FWCore.ParameterSet.Config as cms

treeCreator = cms.EDAnalyzer("TreeMakerModule",
  verbosity = cms.untracked.int32(1),
  createTree = cms.bool(True)
)
