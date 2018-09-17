#ifndef __AnalysisSpace_TreeMaker_EventBlock_h
#define __AnalysisSpace_TreeMaker_EventBlock_h

#include <string>
#include <vector>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "DataFormats/VertexReco/interface/VertexFwd.h"

namespace vhtm {
  class Event;
}
class EventBlock : public edm::one::EDAnalyzer<edm::one::WatchRuns>
{
private:
  void beginJob() override;
  void beginRun(edm::Run const& iRun, edm::EventSetup const& iSetup) override {}
  void endRun(edm::Run const& iRun, edm::EventSetup const& iSetup) override {}
  void analyze(edm::Event const& iEvent, edm::EventSetup const& iSetup) override;
  void endJob() override {}

public:
  explicit EventBlock(const edm::ParameterSet& iConfig);
  virtual ~EventBlock();
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  std::vector<vhtm::Event>* list_;

  const int verbosity_;
  const unsigned int vtxMinNDOF_;
  const double vtxMaxAbsZ_;
  const double vtxMaxd0_;

  const edm::InputTag vertexTag_;
  const edm::InputTag selectedVertexTag_;
  const edm::InputTag puSummaryTag_;
  const edm::InputTag rhoTag_;
  const edm::InputTag rhoNeutralTag_;
  const edm::InputTag fixedGridRhoAllTag_;
  const edm::InputTag fixedGridRhoFastjetAllTag_;
  const edm::InputTag fixedGridRhoFastjetAllCaloTag_;
  const edm::InputTag fixedGridRhoFastjetCentralCaloTag_;
  const edm::InputTag fixedGridRhoFastjetCentralChargedPileUpTag_;
  const edm::InputTag fixedGridRhoFastjetCentralNeutralTag_;
  
  const edm::EDGetTokenT<reco::VertexCollection> vertexToken_;
  const edm::EDGetTokenT<reco::VertexCollection> selectedVertexToken_;
  const edm::EDGetTokenT<std::vector<PileupSummaryInfo> > puSummaryToken_;
  const edm::EDGetTokenT<double> rhoToken_;
  const edm::EDGetTokenT<double> rhoNeutralToken_;
  const edm::EDGetTokenT<double> fixedGridRhoAllToken_;
  const edm::EDGetTokenT<double> fixedGridRhoFastjetAllToken_;
  const edm::EDGetTokenT<double> fixedGridRhoFastjetAllCaloToken_;
  const edm::EDGetTokenT<double> fixedGridRhoFastjetCentralCaloToken_;
  const edm::EDGetTokenT<double> fixedGridRhoFastjetCentralChargedPileUpToken_;
  const edm::EDGetTokenT<double> fixedGridRhoFastjetCentralNeutralToken_;
};
#endif
