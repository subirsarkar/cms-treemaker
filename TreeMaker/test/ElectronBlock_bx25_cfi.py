import FWCore.ParameterSet.Config as cms

electronBlock = cms.EDAnalyzer("ElectronBlock",
  verbosity = cms.untracked.int32(0),
  beamSpotCorr = cms.untracked.bool(True),
  useTrigMode = cms.untracked.bool(False),
  offlineBeamSpot = cms.InputTag('offlineBeamSpot'),
  vertexSrc = cms.InputTag('selectedPrimaryVertices'),
  electronSrc = cms.InputTag('eleMVAproducer'),
  pfCands = cms.InputTag("packedPFCandidates")
)
