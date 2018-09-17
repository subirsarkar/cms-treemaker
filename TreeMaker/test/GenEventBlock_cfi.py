import FWCore.ParameterSet.Config as cms

genEventBlock = cms.EDAnalyzer("GenEventBlock",
         verbosity = cms.untracked.int32(0),
               tag = cms.untracked.string("initrwgt"),
  isStandardSample = cms.untracked.bool(True), # set to True for centrally produced Hgg samples, set to False for HH2B2G samples
            pdfset = cms.untracked.string("NNPDF30_lo_as_0130_nf_4.LHgrid"), # for HH2B2G samples set here the pdf set name
   doAlphasWeights = cms.untracked.bool(True), # set to True for centrally produced Hgg samples, set to False for HH2B2G samples 
   doScaleWeights  = cms.untracked.bool(True),
    nPdfEigWeights = cms.uint32(60),
     mc2hessianCSV = cms.FileInPath('PhysicsTools/HepMCCandAlgos/data/NNPDF30_lo_as_0130_hessian_60.csv'),
   StorePDFWeights = cms.untracked.bool(False),
   GenEventInfoTag = cms.InputTag('generator'),
       LHEEventTag = cms.InputTag('externalLHEProducer', '', 'LHE'),
     PDFWeightsTag = cms.InputTag('pdfWeights', 'cteq66')
)
