#include "TTree.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "DataFormats/METReco/interface/GenMET.h"
#include "DataFormats/METReco/interface/GenMETFwd.h"

#include "AnalysisSpace/TreeMaker/interface/PhysicsObjects.h"
#include "AnalysisSpace/TreeMaker/plugins/GenMETBlock.h"
#include "AnalysisSpace/TreeMaker/interface/Utility.h"

GenMETBlock::GenMETBlock(const edm::ParameterSet& iConfig) :
  verbosity_(iConfig.getUntrackedParameter<int>("verbosity", 0)),
  genMETTag_(iConfig.getParameter<edm::InputTag>("genMETSrc")),
  genMETToken_(consumes<reco::GenMETCollection>(genMETTag_))
{
}
void GenMETBlock::beginJob()
{
  // Get TTree pointer
  TTree* tree = vhtm::Utility::getTree("vhtree");
  list_ = new std::vector<vhtm::GenMET>();
  tree->Branch("GenMET", "std::vector<vhtm::GenMET>", &list_, 32000, -1);
  tree->Branch("nGenMET", &fnGenMET_, "fnGenMET_/I");
}
void GenMETBlock::analyze(edm::Event const& iEvent, edm::EventSetup const& iSetup) {
  // Reset the vector and nObj variables
  list_->clear();
  fnGenMET_ = 0;

  if (!iEvent.isRealData()) {
    edm::Handle<reco::GenMETCollection> metColl;
    bool found = iEvent.getByToken(genMETToken_, metColl);
    if (found && metColl.isValid()) {
      edm::LogInfo("GenMETBlock") << "Total # GenMETs: " << metColl->size();
      for (const auto& v: *metColl) {
        if (list_->size() == kMaxGenMET_) {
          edm::LogInfo("GenMETBlock") << "Too many GenMET, fnGenMET = " << list_->size();
          break;
        }
	vhtm::GenMET genMet;

        // fill in all the vectors
        genMet.met    = v.pt();
        genMet.metphi = v.phi();
        genMet.sumet  = v.sumEt();

        list_->push_back(genMet);
      }
      fnGenMET_ = list_->size(); 
    }
    else {
      edm::LogError("GenMETBlock") << "Error >> Failed to get GenMETCollection for label: "
                                   << genMETTag_;
    }
  }
}
#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(GenMETBlock);
