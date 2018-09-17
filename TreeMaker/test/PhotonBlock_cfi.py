import FWCore.ParameterSet.Config as cms

photonBlock = cms.EDAnalyzer("PhotonBlock",
  verbosity = cms.untracked.int32(0),
  photonSrc = cms.untracked.InputTag('slimmedPhotons'),
  pfCands   = cms.InputTag("packedPFCandidates") 
)
