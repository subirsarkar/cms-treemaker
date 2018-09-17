#ifndef __AnalysisSpace_TreeMaker_PackedPFCandidateBlock_h
#define __AnalysisSpace_TreeMaker_PackedPFCandidateBlock_h

#include <string>
#include <vector>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "DataFormats/Common/interface/Ref.h"

#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Provenance/interface/EventID.h"
#include "FWCore/ParameterSet/interface/ProcessDesc.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "AnalysisSpace/TreeMaker/interface/PhysicsObjects.h"

class PackedPFCandidateBlock: public edm::one::EDAnalyzer<edm::one::WatchRuns> 
{
private:
  void beginJob() override;
  void beginRun(edm::Run const& iRun, edm::EventSetup const& iSetup) override {}
  void endRun(edm::Run const& iRun, edm::EventSetup const& iSetup) override {}
  void analyze(edm::Event const& iEvent, edm::EventSetup const& iSetup) override;
  void endJob() override {}

public:
  explicit PackedPFCandidateBlock(const edm::ParameterSet& iConfig);
  virtual ~PackedPFCandidateBlock();
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
  void calcIsoFromPF(const pat::PackedCandidate& v,
		     edm::Handle<pat::PackedCandidateCollection>& pfs,
		     double cone, std::vector<double>& iso);
  enum {
    kMaxPackedPFCandidate = 1000
  };

private:
  std::vector<vhtm::PackedPFCandidate>* list_;
  int fnPackedPFCandidate_;

  const int verbosity_;
  bool skimCand_;
  std::vector<int> pdgToSaveList_;
  double minCandPt_;

  const edm::InputTag pfCandTag_;
  const edm::InputTag vertexTag_;

  const edm::EDGetTokenT<pat::PackedCandidateCollection> pfToken_;
  const edm::EDGetTokenT<reco::VertexCollection> vertexToken_;
};
#endif
