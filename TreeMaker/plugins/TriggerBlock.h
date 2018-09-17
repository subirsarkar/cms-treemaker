#ifndef __AnalysisSpace_TreeMaker_TriggerBlock_h
#define __AnalysisSpace_TreeMaker_TriggerBlock_h

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
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
#include "HLTrigger/HLTcore/interface/HLTPrescaleProvider.h"
#include "DataFormats/Common/interface/TriggerResults.h"


class TriggerBlock : public edm::one::EDAnalyzer<edm::one::WatchRuns>
		     //class TriggerBlock: public edm::one::EDAnalyzer<edm::one::WatchLuminosityBlocks, edm::one::WatchRuns> 
{
private:
  void beginJob() override;
  void beginRun(edm::Run const& iRun, edm::EventSetup const& iSetup) override;
  void analyze(edm::Event const& iEvent, edm::EventSetup const& iSetup) override;

  //  void beginLuminosityBlock(edm::LuminosityBlock const&,  edm::EventSetup const&) override {}
  //void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override {}
  void endRun(edm::Run const& iEvent, edm::EventSetup const&) override {}
  void endJob() override {}

public:
  explicit TriggerBlock(const edm::ParameterSet& iConfig);
  virtual ~TriggerBlock();
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:

  const int verbosity_;
  const std::vector<std::string> hltPathsOfInterest_;

  const edm::InputTag l1Tag_;
  const edm::InputTag hltTag_;

  std::vector<std::string> matchedPathList_;

  std::vector<int>* l1physbits_;
  std::vector<int>* l1techbits_;
  std::vector<std::string>* hltpaths_;
  std::vector<int>* hltresults_;
  std::vector<int>* hltprescales_;

  const edm::EDGetTokenT<L1GlobalTriggerReadoutRecord> l1Token_;
  const edm::EDGetTokenT<edm::TriggerResults> hltToken_;
  HLTPrescaleProvider hltPrescaleProvider_;
};
#endif
