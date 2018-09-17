import FWCore.ParameterSet.Config as cms

eventBlock = cms.EDAnalyzer("EventBlock",
  verbosity = cms.untracked.int32(0),
  vertexMinimumNDOF = cms.untracked.uint32(4),
  vertexMaxAbsZ = cms.untracked.double(24.),
  vertexMaxd0 = cms.untracked.double(2.),
  hpTrackThreshold = cms.untracked.double(0.25),
  l1InputTag = cms.InputTag('gtDigis'),
  vertexTag = cms.InputTag('offlineSlimmedPrimaryVertices'),
  puSummaryTag = cms.InputTag('slimmedAddPileupInfo'),
  selectedVtxTag = cms.InputTag('selectedPrimaryVertices'),
  rhoTag = cms.InputTag('kt6PFJets','rho'),                         
  rhoNeutralTag = cms.InputTag('kt6PFNeutralJetsForVtxMultReweighting', 'rho'),
  fixedGridRhoAllTag = cms.InputTag("fixedGridRhoAll"),
  fixedGridRhoFastjetAllTag = cms.InputTag("fixedGridRhoFastjetAll"),
  fixedGridRhoFastjetAllCaloTag = cms.InputTag("fixedGridRhoFastjetAllCalo"),
  fixedGridRhoFastjetCentralCaloTag = cms.InputTag("fixedGridRhoFastjetCentralCalo"),
  fixedGridRhoFastjetCentralChargedPileUpTag = cms.InputTag("fixedGridRhoFastjetCentralChargedPileUp"),
  fixedGridRhoFastjetCentralNeutralTag = cms.InputTag("fixedGridRhoFastjetCentralNeutral")
)
