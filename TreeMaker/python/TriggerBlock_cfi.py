import FWCore.ParameterSet.Config as cms

triggerBlock = cms.EDAnalyzer("TriggerBlock",
  verbosity = cms.untracked.int32(0),
  hltPathsOfInterest = cms.vstring ( 
        # Single Lepton:
        'HLT_Ele25_eta2p1_WPTight_Gsf_v',
        'HLT_Ele27_WPTight_Gsf_v',
        'HLT_Ele27_eta2p1_WPLoose_Gsf_v',
        'HLT_Ele32_eta2p1_WPTight_Gsf_v',
        'HLT_IsoMu20_v',
        'HLT_IsoTkMu20_v',
        'HLT_IsoMu22_v',
        'HLT_IsoTkMu22_v',
        'HLT_IsoMu24_v',
        'HLT_IsoTkMu24_v',
        # Dilepton
        'HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v',
        'HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v',
        'HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v',
        'HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v',
        'HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL_v',
        'HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v',
        'HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v',
        'HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v',
        'HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_v',
        'HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_v',
        # TriLepton
        'HLT_TripleMu_12_10_5_v',
        'HLT_Ele16_Ele12_Ele8_CaloIdL_TrackIdL_v',
        'HLT_DiMu9_Ele9_CaloIdL_TrackIdL_v',
        'HLT_Mu8_DiEle12_CaloIdL_TrackIdL_v'
  ),
  l1Tag = cms.InputTag('gtDigis'),
  hltTag = cms.InputTag('TriggerResults','','HLT')
)

triggerBlockZTnP = cms.EDAnalyzer("TriggerBlock",
  verbosity = cms.untracked.int32(0),
  hltPathsOfInterest = cms.vstring ( 
 	'HLT_IsoMu20_v*',
  	'HLT_IsoTkMu20_v*',
        'HLT_IsoMu22_v*',
        'HLT_IsoTkMu22_v*'
  ),
  l1Tag = cms.InputTag('gtDigis'),
  hltTag = cms.InputTag('TriggerResults','','HLT')
)
