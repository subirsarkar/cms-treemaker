#include <iostream>
#include <string>
#include <vector>
#include <cassert>
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "TTree.h"

#include "AnalysisSpace/TreeMaker/plugins/TreeMakerModule.h"
#include "AnalysisSpace/TreeMaker/interface/PhysicsObjects.h"
#include "AnalysisSpace/TreeMaker/interface/Utility.h"

TreeMakerModule::TreeMakerModule(const edm::ParameterSet& iConfig) :
  verbosity_(iConfig.getUntrackedParameter<int>("verbosity", 0)),
  createTree_(iConfig.getParameter<bool>("createTree"))
{
}
void TreeMakerModule::beginJob()
{
  if (!createTree_) return;
  edm::Service<TFileService> fs;
  fs->file().cd("/");
  TTree* tree = fs->make<TTree>("vhtree", "Physics Analysis Level TTree");
  assert(tree);
  fs->file().ls();
}
void TreeMakerModule::analyze(edm::Event const& iEvent, edm::EventSetup const& iSetup) {
  // Get TTree pointer
  if (createTree_) return;
  TTree* tree = vhtm::Utility::getTree("vhtree");
  if (tree) tree->Fill();
}
void TreeMakerModule::endJob() {
}
// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
TreeMakerModule::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}
#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(TreeMakerModule);
