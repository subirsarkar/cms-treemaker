import FWCore.ParameterSet.Config as cms

tauBlock = cms.EDAnalyzer("TauBlock",
  verbosity = cms.untracked.int32(0),
  beamSpotCorr = cms.untracked.bool(True),
  patTauSrc = cms.InputTag('slimmedTaus'),
  vertexSrc = cms.InputTag('selectedPrimaryVertices'),
  offlineBeamSpot = cms.InputTag('offlineBeamSpot')
)
