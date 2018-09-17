#ifndef __AnalysisSpace_TreeMaker_MuonBlock_h
#define __AnalysisSpace_TreeMaker_MuonBlock_h

#include <string>
#include <vector>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/Utils/interface/StringCutObjectSelector.h"

namespace vhtm {
  class Muon;
}
class MuonBlock : public edm::one::EDAnalyzer<edm::one::WatchRuns>
{
 private:
  void beginJob() override;
  void beginRun(edm::Run const& iRun, edm::EventSetup const& iSetup) override {}
  void endRun(edm::Run const& iRun, edm::EventSetup const& iSetup) override {}
  void analyze(edm::Event const& iEvent, edm::EventSetup const& iSetup) override;
  void endJob() override {}

  void calcIsoFromPF(double cone, edm::Handle<pat::PackedCandidateCollection>& pfs, const pat::Muon& v, std::vector<double>& iso);     
  bool isTrackerHighPt(const pat::Muon & mu, const reco::Vertex & primaryVertex);

 public:
  explicit MuonBlock(const edm::ParameterSet& iConfig);
  virtual ~MuonBlock();
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
  enum {
    kMaxMuon_ = 100
  };

 private:
  std::vector<vhtm::Muon>* list_;
  int fnMuon_;

  const int verbosity_;
  const bool bsCorr_;
  const std::string muonID_;

  const edm::InputTag muonTag_;
  const edm::InputTag vertexTag_;
  const edm::InputTag bsTag_;
  const edm::InputTag pfcandTag_;
  
  const edm::EDGetTokenT<pat::MuonCollection> muonToken_;
  const edm::EDGetTokenT<reco::VertexCollection> vertexToken_;
  const edm::EDGetTokenT<reco::BeamSpot> bsToken_;
  const edm::EDGetTokenT<pat::PackedCandidateCollection> pfToken_;

  // Use default criteria to choose the best muon
  bool defaultBestMuon_;

  // Cut on the pair of objects together
  typedef std::pair<const reco::Muon *, const reco::Muon *> MuonPointerPair;
  StringCutObjectSelector<MuonPointerPair, true> bestMuonSelector_;

  bool isSameMuon(const pat::Muon &mu1, const pat::Muon &mu2) const;
  bool isBetterMuon(const pat::Muon &mu1, const pat::Muon &mu2) const;
};
#endif
