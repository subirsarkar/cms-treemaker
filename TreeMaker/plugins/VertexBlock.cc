#include <iostream>
#include <iomanip>
#include <sstream>

#include "TTree.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "DataFormats/VertexReco/interface/Vertex.h"

#include "AnalysisSpace/TreeMaker/interface/PhysicsObjects.h"
#include "AnalysisSpace/TreeMaker/plugins/VertexBlock.h"
#include "AnalysisSpace/TreeMaker/interface/Utility.h"

VertexBlock::VertexBlock(const edm::ParameterSet& iConfig) :
  verbosity_(iConfig.getUntrackedParameter<int>("verbosity", 0)),
  vertexTag_(iConfig.getParameter<edm::InputTag>("vertexSrc")),
  vertexToken_(consumes<reco::VertexCollection>(vertexTag_))
{}
VertexBlock::~VertexBlock()
{
  delete list_;
}
void VertexBlock::beginJob() {
  // Get TTree pointer
  TTree* tree = vhtm::Utility::getTree("vhtree");
  list_ = new std::vector<vhtm::Vertex>();
  tree->Branch("Vertex", "std::vector<vhtm::Vertex>", &list_, 32000, -1);
  tree->Branch("nVertex", &fnVertex_, "fnVertex_/I");
}
void VertexBlock::analyze(edm::Event const& iEvent, edm::EventSetup const& iSetup) {
  // Reset the vector and the nObj variables
  list_->clear();
  fnVertex_ = 0;

  edm::Handle<reco::VertexCollection> primaryVertices;
  bool found = iEvent.getByToken(vertexToken_, primaryVertices);

  if (found && primaryVertices.isValid()) {
    edm::LogInfo("VertexBlock") << "Total # of Primary Vertices: " 
                                << primaryVertices->size();
    if (verbosity_) {
      std::cout << std::setprecision(2);
      std::cout << "   indx      x      y      z    rho     chi2      ndf" << std::endl;
    }
    for (const auto& v: *primaryVertices) {
      if (list_->size() == kMaxVertex_) {
        edm::LogInfo("VertexBlock") << "Too many Vertex, fnVertex = " 
                                    << list_->size();
        break;
      }
      vhtm::Vertex vertex;
      vertex.x          = v.x();
      vertex.y          = v.y();
      vertex.z          = v.z();
      vertex.xErr       = v.xError();
      vertex.yErr       = v.yError();
      vertex.zErr       = v.zError();
      vertex.rho        = v.position().rho();
      vertex.chi2       = v.chi2();
      vertex.ndf        = v.ndof();
      vertex.isfake     = v.isFake();
      vertex.isvalid    = v.isValid();
      //vertex.sumPt      = v.p4().pt();

      if (verbosity_)
	std::cout << std::setw(7) << list_->size()
                  << std::setw(7) << vertex.x 
                  << std::setw(7) << vertex.y 
                  << std::setw(7) << vertex.y
                  << std::setw(7) << vertex.rho
		  << std::setw(9) << vertex.chi2
		  << std::setw(9) << vertex.ndf
                  << std::endl;
      list_->push_back(vertex);
    }
    fnVertex_ = list_->size();
  }
  else {
    edm::LogError("VertexBlock") << "Error! Failed to get VertexCollection for label: "
                                 << vertexTag_;
  }
}
// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
VertexBlock::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}
#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(VertexBlock);
