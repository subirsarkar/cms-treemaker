import FWCore.ParameterSet.Config as cms

muonBlock = cms.EDAnalyzer("MuonBlock",
  verbosity = cms.untracked.int32(0),
  beamSpotCorr = cms.untracked.bool(True),
  muonID = cms.untracked.string('GlobalMuonPromptTight'),
  muonSrc = cms.InputTag('gcleanMuons'),
  vertexSrc = cms.InputTag('selectedPrimaryVertices'),
  offlineBeamSpot = cms.InputTag('offlineBeamSpot'),
  pfCands = cms.InputTag("packedPFCandidates")
)
