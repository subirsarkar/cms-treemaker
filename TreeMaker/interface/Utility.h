#ifndef __AnalysisSpace_TreeMaker_Utility_h
#define __AnalysisSpace_TreeMaker_Utility_h

#include <string>
#include "TTree.h"
class TTree;

namespace vhtm {
  class Utility {
  public:
    static TTree* getTree(const std::string& tree_name="vhtree");
  };
}
#endif
