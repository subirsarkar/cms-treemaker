#include <iostream>
#include <iomanip>
#include <string>
#include <sstream>

#include "TTree.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "AnalysisSpace/TreeMaker/interface/PhysicsObjects.h"
#include "AnalysisSpace/TreeMaker/plugins/GenParticleBlock.h"
#include "AnalysisSpace/TreeMaker/interface/Utility.h"

using namespace std;

GenParticleBlock::GenParticleBlock(const edm::ParameterSet& iConfig) :
  verbosity_(iConfig.getUntrackedParameter<int>("verbosity", 0)),
  genParticleTag_(iConfig.getParameter<edm::InputTag>("genParticleSrc")),
  genParticleToken_(consumes<reco::GenParticleCollection>(genParticleTag_))
{
}
GenParticleBlock::~GenParticleBlock() {
  delete list_;
}
void GenParticleBlock::beginJob() {
  // Get TTree pointer
  TTree* tree = vhtm::Utility::getTree("vhtree");
  list_ = new std::vector<vhtm::GenParticle>();
  tree->Branch("GenParticle", "std::vector<vhtm::GenParticle>", &list_, 32000, -1);
  tree->Branch("nGenParticle", &fnGenParticle_, "fnGenParticle_/I");
}
void GenParticleBlock::analyze(edm::Event const& iEvent, edm::EventSetup const& iSetup) {
  // Reset the vector and the nObj variables
  list_->clear();
  fnGenParticle_ = 0;

  if (!iEvent.isRealData()) {
    edm::Handle<reco::GenParticleCollection> genParticles;
    bool found = iEvent.getByToken(genParticleToken_, genParticles);
    if (found && genParticles.isValid()) {
      edm::LogInfo("GenParticleBlock") << "Total # GenParticles: " << genParticles->size();
      if (verbosity_ > 0) {
	std::cout << std::setprecision(2);
        std::cout << "indx    status    pdgId  charge     eta      phi      pt     energy             mID                             dID"
		  << std::endl;
      }
      //      for (const reco:;GenParticle& v:auto it = genParticles->begin(); it != genParticles->end(); ++it) {
      for (const auto& v: *genParticles) {
        if (list_->size() == kMaxGenParticle_) {
           edm::LogInfo("GenParticleBlock") << "Too many GenParticles, fnGenParticle = "
                                            << list_->size();
           break;
        }
        int pdgid = v.pdgId();
        double pt = v.pt();

        vhtm::GenParticle gp;

        // fill in all the vectors
        gp.eta       = v.eta();
        gp.phi       = v.phi();
        gp.p         = v.p();
        gp.px        = v.px();
        gp.py        = v.py();
        gp.pz        = v.pz();
        gp.pt        = pt;
        gp.energy    = v.energy();
        gp.pdgId     = pdgid;
        gp.vx        = v.vx();
        gp.vy        = v.vy();
        gp.vz        = v.vz();
        gp.status    = v.status();
        gp.charge    = v.charge();
        gp.numDaught = v.numberOfDaughters();
        gp.numMother = v.numberOfMothers();

        // First mother
        int idx = -1;
        for (auto mit = genParticles->begin(); mit != genParticles->end(); ++mit) {
          if (v.mother() == &(*mit)) {
            idx = std::distance(genParticles->begin(), mit);
            break;
          }
        }
        gp.motherIndex = idx;

        // Mothers
        ostringstream mID;
        for (size_t j = 0; j < v.numberOfMothers(); ++j) {
          const reco::Candidate* m = v.mother(j);
          for (auto mit = genParticles->begin(); mit != genParticles->end(); ++mit) {
            if (m == &(*mit) ) {
              int idx = std::distance(genParticles->begin(), mit);
              gp.motherIndices.push_back(idx);
 	      mID << " " << idx; 
              break;
            }
          }
        }

        // Daughters
        ostringstream dID;
        for (size_t j = 0; j < v.numberOfDaughters(); ++j) {
          const reco::Candidate* d = v.daughter(j);
          for (auto mit = genParticles->begin(); mit != genParticles->end(); ++mit) {
            if (d == &(*mit) ) {
              int idx = std::distance(genParticles->begin(), mit);
              gp.daughtIndices.push_back(idx);
	      dID << " " << idx; 
              break;
            }
          }
        }

        // Dump
        if (verbosity_ > 0) {
          string ms = mID.str();
	  if (!ms.length()) ms = " -";

          string ds = dID.str();
	  if (!ds.length()) ds = " -";

	  std::cout << setw(4)  << fnGenParticle_++
		    << setw(8)  << v.status()
		    << setw(10) << v.pdgId()
                    << setw(8)  << v.charge()
		    << setw(10) << v.eta()
		    << setw(9)  << v.phi()
		    << setw(9)  << v.pt()
		    << setw(9)  << v.energy()
		    << setw(16) << ms 
		    << ds
		    << std::endl;
        }
        // add particle to the list
        list_->push_back(gp);
      }
      fnGenParticle_ = list_->size();
    }
    else {
      edm::LogError("GenParticleBlock") << "Error >> Failed to get GenParticleCollection for label: "
                                        << genParticleTag_;
    }
  }
}
// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
GenParticleBlock::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}
DEFINE_FWK_MODULE(GenParticleBlock);
