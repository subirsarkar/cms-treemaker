#ifndef __AnalysisSpace_TreeMaker_GenMETBlock_h
#define __AnalysisSpace_TreeMaker_GenMETBlock_h

#include <string>
#include <vector>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "AnalysisSpace/TreeMaker/interface/PhysicsObjects.h"

namespace vhtm {
  class GenMET;
}

class GenMETBlock : public edm::one::EDAnalyzer<edm::one::WatchRuns>
{
private:
  void beginJob() override;
  void beginRun(edm::Run const& iRun, edm::EventSetup const& iSetup) override {}
  void endRun(edm::Run const& iRun, edm::EventSetup const& iSetup) override {}
  void analyze(edm::Event const& iEvent, edm::EventSetup const& iSetup) override;
  void endJob() override {}

public:
  explicit GenMETBlock(const edm::ParameterSet& iConfig);
  virtual ~GenMETBlock() {}

  enum {
    kMaxGenMET_ = 5
  };

private:
  std::vector<vhtm::GenMET>* list_;
  int fnGenMET_;

  const int verbosity_;
  const edm::InputTag genMETTag_;
  const edm::EDGetTokenT<reco::GenMETCollection> genMETToken_;
};
#endif
