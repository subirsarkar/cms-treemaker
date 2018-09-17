#ifndef __AnalysisSpace_TreeMaker_TreeMakerModule_h
#define __AnalysisSpace_TreeMaker_TreeMakerModule_h

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

class TreeMakerModule : public edm::one::EDAnalyzer<edm::one::WatchRuns>
{
private:
  void beginJob();
  void beginRun(edm::Run const& iRun, edm::EventSetup const& iSetup) override {};
  void endRun(edm::Run const& iRun, edm::EventSetup const& iSetup) override {}
  void analyze(edm::Event const& iEvent, edm::EventSetup const& iSetup) override;
  void endJob() override;

public:
  explicit TreeMakerModule(const edm::ParameterSet& iConfig);
  virtual ~TreeMakerModule() {}
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
private:
  const int verbosity_;
  const bool createTree_;
};
#endif
