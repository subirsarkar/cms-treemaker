#ifndef __AnalysisSpace_TreeMaker_GenEventBlock_h
#define __AnalysisSpace_TreeMaker_GenEventBlock_h

#include <string>
#include <vector>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "DataFormats/Common/interface/Ref.h"

#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Provenance/interface/EventID.h"
#include "PhysicsTools/HepMCCandAlgos/interface/PDFWeightsHelper.h"

namespace vhtm {
  class GenEvent;
}
class GenEventBlock : public edm::one::EDAnalyzer<edm::one::WatchRuns>
{
private:
  void beginJob() override;
  void beginRun(edm::Run const& iRun, edm::EventSetup const& iSetup) override;
  void analyze(edm::Event const& iEvent, edm::EventSetup const& iSetup) override;
  void endRun(edm::Run const& iEvent, edm::EventSetup const&) override {}
  void endJob() override {}

public:
  explicit GenEventBlock(const edm::ParameterSet& iConfig);
  virtual ~GenEventBlock();

private:
  std::vector<vhtm::GenEvent>* list_;

  const int verbosity_;
  const bool storePDFWeights_;

  // shamelessly taken from flashgg::PDFWeightObjectProducer
  const std::string tag_;
  const bool isStandardSample_;
  const bool doAlphasWeights_;
  const bool doScaleWeights_;
  const unsigned int nPdfEigWeights_;
  const std::string pdfset_;
  const edm::FileInPath mc2hessianCSV_;

  std::string pdfid_1;
  std::string pdfid_2;
  std::string alphas_id_1;
  std::string alphas_id_2;

  std::vector<int> pdf_indices;
  std::vector<int> alpha_indices;
  std::vector<int> scale_indices;
  PDFWeightsHelper pdfWeightsHelper_; // tool from HepMCCandAlgos/interface/PDFWeightsHelper
  std::string generatorType_;
  std::vector<double> inpdfweights;

  const edm::InputTag genEventTag_;
  const edm::InputTag lheEventTag_;
  const edm::InputTag pdfWeightsTag_;
  const edm::EDGetTokenT<GenEventInfoProduct> genEventToken_;
  const edm::EDGetTokenT<LHEEventProduct> lheEventToken_;
  const edm::EDGetTokenT<LHERunInfoProduct> lheRunToken_;
  const edm::EDGetTokenT<std::vector<double>> pdfWeightsToken_;
};
#endif
