#ifndef __AnalysisSpace_TreeMaker_ElectronBlock_h
#define __AnalysisSpace_TreeMaker_ElectronBlock_h

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
#include "FWCore/ParameterSet/interface/ProcessDesc.h"
#include "DataFormats/Common/interface/ValueMap.h"

#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"

namespace vhtm {
  class Electron;
}
class ElectronBlock : public edm::one::EDAnalyzer<edm::one::WatchRuns>
{
private:
  void beginJob() override;
  void beginRun(edm::Run const& iRun, edm::EventSetup const& iSetup) override {}
  void endRun(edm::Run const& iRun, edm::EventSetup const& iSetup) override {}
  void analyze(edm::Event const& iEvent, edm::EventSetup const& iSetup) override;
  void endJob() override {}

public:
  explicit ElectronBlock(const edm::ParameterSet& iConfig);
  virtual ~ElectronBlock();
  void calcIsoFromPF(double cone, edm::Handle<pat::PackedCandidateCollection>& pfs, const pat::Electron& v, std::vector<double>& iso);     
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

  enum {
    kMaxElectron_ = 100
  };

private:
  std::vector<vhtm::Electron>* list_;
  int fnElectron_;

  int verbosity_;
  bool bsCorr_;
  bool trigMode_;

  const edm::InputTag bsTag_;
  const edm::InputTag vertexTag_;
  const edm::InputTag electronTag_;
  const edm::InputTag pfcandTag_;
  
  const edm::EDGetTokenT<reco::BeamSpot> bsToken_;
  const edm::EDGetTokenT<reco::VertexCollection> vertexToken_;
  const edm::EDGetTokenT<pat::ElectronCollection> electronToken_;
  const edm::EDGetTokenT<pat::PackedCandidateCollection> pfToken_;

  // ID decisions objects
  edm::EDGetTokenT<edm::ValueMap<bool>> eleMediumIdMapToken_;
  edm::EDGetTokenT<edm::ValueMap<bool>> eleTightIdMapToken_;
  edm::EDGetTokenT<edm::ValueMap<bool>> eleHZZIdMapToken_;
  
  const edm::EDGetToken gsfelectronTokenMVAId_;
};
#endif
