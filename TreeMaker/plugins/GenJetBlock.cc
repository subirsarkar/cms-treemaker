#include <iostream>

#include "TTree.h"
#include "DataFormats/JetReco/interface/GenJetCollection.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "AnalysisSpace/TreeMaker/plugins/GenJetBlock.h"
#include "AnalysisSpace/TreeMaker/interface/PhysicsObjects.h"
#include "AnalysisSpace/TreeMaker/interface/Utility.h"

GenJetBlock::GenJetBlock(const edm::ParameterSet& iConfig) :
  verbosity_(iConfig.getUntrackedParameter<int>("verbosity", 0)),
  genJetTag_(iConfig.getParameter<edm::InputTag>("genJetSrc")),
  genJetToken_(consumes<reco::GenJetCollection>(genJetTag_))   
{
}
void GenJetBlock::beginJob() {
  // Get TTree pointer
  TTree* tree = vhtm::Utility::getTree("vhtree");
  list_ = new std::vector<vhtm::GenJet>();
  tree->Branch("GenJet", "std::vector<vhtm::GenJet>", &list_, 32000, -1);
  tree->Branch("nGenJet", &fnGenJet_, "fnGenJet_/I");
}
void GenJetBlock::analyze(edm::Event const& iEvent, edm::EventSetup const& iSetup) {
  // Reset the vector and the nObj variables
  list_->clear();
  fnGenJet_ = 0;

  if (!iEvent.isRealData()) {
    edm::Handle<reco::GenJetCollection> genJets;
    bool found = iEvent.getByToken(genJetToken_, genJets);

    if (found && genJets.isValid()) {
      edm::LogInfo("GenJetBlock") << "Total # GenJets: " << genJets->size();
      for (const auto& v: *genJets) {
        if (list_->size() == kMaxGenJet) {
	  edm::LogInfo("GenJetBlock") << "Too many Gen Jets, fnGenJet = " << list_->size();
	  break;
        }
	vhtm::GenJet jet;

        // fill in all the vectors
        jet.eta = v.eta();
        jet.phi = v.phi();
        jet.p   = v.p();
        jet.pt  = v.pt();
        double energy = v.energy();
        jet.energy = energy;
        jet.emf    = (energy > 0) ? v.emEnergy()/energy : 0;
        jet.hadf   = (energy > 0) ? v.hadEnergy()/energy : 0;
 
        list_->push_back(jet);
      }
      fnGenJet_ = list_->size();
    }
    else {
      edm::LogError("GenJetBlock") << "Error >> Failed to get GenJetCollection for label: "
                                   << genJetTag_;
    }
  }
}
#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(GenJetBlock);
