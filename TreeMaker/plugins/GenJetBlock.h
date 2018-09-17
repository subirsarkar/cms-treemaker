#ifndef __AnalysisSpace_TreeMaker_GenJetBlock_h
#define __AnalysisSpace_TreeMaker_GenJetBlock_h

#include <string>
#include <vector>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

namespace vhtm {
  class GenJet;
}

class GenJetBlock : public edm::one::EDAnalyzer<edm::one::WatchRuns>
{
private:
  void beginJob() override;
  void beginRun(edm::Run const& iRun, edm::EventSetup const& iSetup) override {}
  void endRun(edm::Run const& iRun, edm::EventSetup const& iSetup) override {}
  void analyze(edm::Event const& iEvent, edm::EventSetup const& iSetup) override;
  void endJob() override {}

public:
  explicit GenJetBlock(const edm::ParameterSet& iConfig);
  virtual ~GenJetBlock() {}

  enum {
    kMaxGenJet = 100
  };

private:
  std::vector<vhtm::GenJet>* list_;
  int fnGenJet_;

  const int verbosity_;
  const edm::InputTag genJetTag_;
  const edm::EDGetTokenT<reco::GenJetCollection> genJetToken_;
};
#endif
