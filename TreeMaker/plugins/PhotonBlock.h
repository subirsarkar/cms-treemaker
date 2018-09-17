#ifndef __AnalysisSpace_TreeMaker_PhotonBlock_h
#define __AnalysisSpace_TreeMaker_PhotonBlock_h

#include <string>
#include <vector>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "DataFormats/Common/interface/Ref.h"

#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Provenance/interface/EventID.h"
#include "FWCore/ParameterSet/interface/ProcessDesc.h"

#include "AnalysisSpace/TreeMaker/interface/PhysicsObjects.h"

namespace vhtm {
  class Photon;
}

class PhotonBlock : public edm::one::EDAnalyzer<edm::one::WatchRuns> 
{
 private:
  void beginJob() override;
  void beginRun(edm::Run const& iRun, edm::EventSetup const& iSetup) override {}
  void endRun(edm::Run const& iRun, edm::EventSetup const& iSetup) override {}
  void analyze(edm::Event const& iEvent, edm::EventSetup const& iSetup) override;
  void endJob() override {}

  void calcIsoFromPF(double cone, edm::Handle<pat::PackedCandidateCollection>& pfs, const pat::Photon& v, std::vector<double>& iso);

 public:
  explicit PhotonBlock(const edm::ParameterSet& iConfig);
  virtual ~PhotonBlock() {}

  enum {
    kMaxPhoton = 100
  };
 private:
  std::vector<vhtm::Photon>* list_;
  int fnPhoton_;

  int verbosity_;
  const edm::InputTag photonTag_;
  const edm::InputTag pfcandTag_;
  const edm::EDGetTokenT<pat::PhotonCollection> photonToken_;
  const edm::EDGetTokenT<pat::PackedCandidateCollection> pfToken_;
};
#endif
