import FWCore.ParameterSet.Config as cms

packedPFCandidateBlock = cms.EDAnalyzer("PackedPFCandidateBlock",
  verbosity = cms.untracked.int32(0),
  skimCandidate = cms.untracked.bool(True),
  pdgTosave = cms.vint32(22), # 11,13,211
  minCandPt = cms.untracked.double(2.0), 
  pfCands   = cms.InputTag('packedPFCandidates'),
  vertexSrc = cms.InputTag('selectedPrimaryVertices')
)
