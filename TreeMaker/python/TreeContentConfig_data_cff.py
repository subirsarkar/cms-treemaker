import FWCore.ParameterSet.Config as cms

from AnalysisSpace.TreeMaker.EventBlock_cfi import eventBlock
from AnalysisSpace.TreeMaker.TriggerBlock_cfi import triggerBlock,triggerBlockZTnP
from AnalysisSpace.TreeMaker.TriggerObjectBlock_cfi import triggerObjectBlock
from AnalysisSpace.TreeMaker.VertexBlock_cfi import vertexBlock
from AnalysisSpace.TreeMaker.ElectronBlock_cfi import electronBlock
from AnalysisSpace.TreeMaker.MuonBlock_cfi import muonBlock
from AnalysisSpace.TreeMaker.TauBlock_cfi import tauBlock
from AnalysisSpace.TreeMaker.METBlock_cfi import metBlock
from AnalysisSpace.TreeMaker.JetBlock_cfi import jetBlock
from AnalysisSpace.TreeMaker.PackedPFCandidateBlock_cfi import packedPFCandidateBlock

treeContentSequence = cms.Sequence(
   eventBlock
 + triggerBlock
 + triggerObjectBlock
 + vertexBlock
 + electronBlock
 + muonBlock
 + tauBlock
 + metBlock
 + jetBlock
 + packedPFCandidateBlock
)

treeContentSequenceZTnPData = cms.Sequence(
  eventBlock
  + triggerBlockZTnP 
  + triggerObjectBlock
  + vertexBlock
  + electronBlock
  + muonBlock
  + packedPFCandidateBlock
)
