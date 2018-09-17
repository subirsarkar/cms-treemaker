#include "TTree.h"

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "AnalysisSpace/TreeMaker/interface/Utility.h"

TTree* vhtm::Utility::getTree(const std::string& tree_name) {
  edm::Service<TFileService> fs;
  TTree *tree = dynamic_cast<TTree*>(fs->file().FindObjectAny(tree_name.c_str()));
  assert(tree);
  return tree;
}
