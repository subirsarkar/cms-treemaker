#ifndef __AnalysisSpace_TreeMaker_VertexBlock_h
#define __AnalysisSpace_TreeMaker_VertexBlock_h

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
  class Vertex;
}
class VertexBlock : public edm::one::EDAnalyzer<edm::one::WatchRuns>
{
private:
  void beginJob();
  void beginRun(edm::Run const& iRun, edm::EventSetup const& iSetup) override {};
  void endRun(edm::Run const& iRun, edm::EventSetup const& iSetup) override {}
  void analyze(edm::Event const& iEvent, edm::EventSetup const& iSetup) override;
  void endJob() override {}

public:
  explicit VertexBlock(const edm::ParameterSet& iConfig);
  virtual ~VertexBlock();
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
  enum {
    kMaxVertex_ = 150
  };

private:
  std::vector<vhtm::Vertex>* list_;
  int fnVertex_;
  int verbosity_;
  edm::InputTag vertexTag_;
  edm::EDGetTokenT<reco::VertexCollection> vertexToken_;
};
#endif
