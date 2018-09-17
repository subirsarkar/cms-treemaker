import FWCore.ParameterSet.Config as cms

vertexBlock = cms.EDAnalyzer("VertexBlock",
  verbosity = cms.untracked.int32(0),
  vertexSrc = cms.InputTag('selectedPrimaryVertices')
)
