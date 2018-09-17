#ifndef __AnalysisSpace_TreeMaker_METBlock_h
#define __AnalysisSpace_TreeMaker_METBlock_h

#include <string>
#include <vector>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "DataFormats/PatCandidates/interface/MET.h"

namespace vhtm {
  class MET;
}
class METBlock : public edm::one::EDAnalyzer<edm::one::WatchRuns>
{
private:
  void beginJob() override;
  void beginRun(edm::Run const& iRun, edm::EventSetup const& iSetup) override {}
  void endRun(edm::Run const& iRun, edm::EventSetup const& iSetup) override {}
  void analyze(edm::Event const& iEvent, edm::EventSetup const& iSetup) override;
  void endJob() override {}

public:
  explicit METBlock(const edm::ParameterSet& iConfig);
  virtual ~METBlock() {}

  enum {
    kMaxMET_ = 5
  };

  void fillMET(const edm::Event& iEvent,
               const edm::EventSetup& iSetup,
               std::vector<vhtm::MET>* list,
               int& fnMET,
               const edm::InputTag& iTag,  
               const edm::EDGetTokenT<pat::METCollection>& token);

private:
  std::vector<vhtm::MET>* pfList_;
  int fnPFMET_;

  std::vector<vhtm::MET>* corrList_;
  int fnCorrMET_;

  std::vector<vhtm::MET>* puppiList_;
  int fnPuppiMET_;

  const int verbosity_;
  const edm::InputTag pfMETTag_;
  const edm::InputTag corrMETTag_;
  const edm::InputTag puppiMETTag_;

  const edm::EDGetTokenT<pat::METCollection> pfMETToken_;
  const edm::EDGetTokenT<pat::METCollection> corrMETToken_;
  const edm::EDGetTokenT<pat::METCollection> puppiMETToken_;
};
#endif
