#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutRecord.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GtFdlWord.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"

#include "TTree.h"

#include "AnalysisSpace/TreeMaker/interface/PhysicsObjects.h"
#include "AnalysisSpace/TreeMaker/plugins/EventBlock.h"
#include "AnalysisSpace/TreeMaker/interface/Utility.h"

EventBlock::EventBlock(const edm::ParameterSet& iConfig):
  verbosity_(iConfig.getUntrackedParameter<int>("verbosity", 0)),
  vtxMinNDOF_(iConfig.getUntrackedParameter<unsigned int>("vertexMinimumNDOF", 4)),
  vtxMaxAbsZ_(iConfig.getUntrackedParameter<double>("vertexMaxAbsZ", 24.)),
  vtxMaxd0_(iConfig.getUntrackedParameter<double>("vertexMaxd0", 2.0)),
  vertexTag_(iConfig.getParameter<edm::InputTag>("vertexTag")),
  selectedVertexTag_(iConfig.getParameter<edm::InputTag>("selectedVtxTag")),
  puSummaryTag_(iConfig.getParameter<edm::InputTag>("puSummaryTag")),
  rhoTag_(iConfig.getParameter<edm::InputTag>("rhoTag")),
  rhoNeutralTag_(iConfig.getParameter<edm::InputTag>("rhoNeutralTag")),
  fixedGridRhoAllTag_(iConfig.getParameter<edm::InputTag>("fixedGridRhoAllTag")),
  fixedGridRhoFastjetAllTag_(iConfig.getParameter<edm::InputTag>("fixedGridRhoFastjetAllTag")),
  fixedGridRhoFastjetAllCaloTag_(iConfig.getParameter<edm::InputTag>("fixedGridRhoFastjetAllCaloTag")),
  fixedGridRhoFastjetCentralCaloTag_(iConfig.getParameter<edm::InputTag>("fixedGridRhoFastjetCentralCaloTag")),
  fixedGridRhoFastjetCentralChargedPileUpTag_(iConfig.getParameter<edm::InputTag>("fixedGridRhoFastjetCentralChargedPileUpTag")),
  fixedGridRhoFastjetCentralNeutralTag_(iConfig.getParameter<edm::InputTag>("fixedGridRhoFastjetCentralNeutralTag")),
  vertexToken_(consumes<reco::VertexCollection>(vertexTag_)),
  selectedVertexToken_(consumes<reco::VertexCollection>(selectedVertexTag_)),
  puSummaryToken_(consumes<std::vector<PileupSummaryInfo> >(puSummaryTag_)),
  rhoToken_(consumes<double>(rhoTag_)),
  rhoNeutralToken_(consumes<double>(rhoNeutralTag_)),
  fixedGridRhoAllToken_(consumes<double>(fixedGridRhoAllTag_)),
  fixedGridRhoFastjetAllToken_(consumes<double>(fixedGridRhoFastjetAllTag_)),
  fixedGridRhoFastjetAllCaloToken_(consumes<double>(fixedGridRhoFastjetAllCaloTag_)),
  fixedGridRhoFastjetCentralCaloToken_(consumes<double>(fixedGridRhoFastjetCentralCaloTag_)),
  fixedGridRhoFastjetCentralChargedPileUpToken_(consumes<double>(fixedGridRhoFastjetCentralChargedPileUpTag_)),
  fixedGridRhoFastjetCentralNeutralToken_(consumes<double>(fixedGridRhoFastjetCentralNeutralTag_))  
{
}
EventBlock::~EventBlock() {
  delete list_;
}
void EventBlock::beginJob() {
  // Get TTree pointer
  TTree* tree = vhtm::Utility::getTree("vhtree");

  list_ = new std::vector<vhtm::Event>();
  tree->Branch("Event", "std::vector<vhtm::Event>", &list_, 32000, -1);
}
void EventBlock::analyze(edm::Event const& iEvent, edm::EventSetup const& iSetup) {
  // Reset the vector
  list_->clear();

  // Create Event Object
  vhtm::Event ev;
  ev.run   = iEvent.id().run();
  ev.event = iEvent.id().event();
  ev.lumis = iEvent.id().luminosityBlock();
  ev.bunch = iEvent.bunchCrossing();
  ev.orbit = iEvent.orbitNumber();

  double sec = iEvent.time().value() >> 32 ;
  double usec = 0xFFFFFFFF & iEvent.time().value();
  double conv = 1e6;
  ev.time   = sec+usec/conv;
  ev.isdata = iEvent.isRealData();
 
  // Rho Part
  edm::Handle<double> fixedGridRhoAll;
  iEvent.getByToken(fixedGridRhoAllToken_,fixedGridRhoAll);
  ev.fGridRhoAll = *fixedGridRhoAll;

  edm::Handle<double> fixedGridRhoFastjetAll;
  iEvent.getByToken(fixedGridRhoFastjetAllToken_,fixedGridRhoFastjetAll);
  ev.fGridRhoFastjetAll = *fixedGridRhoFastjetAll;

  edm::Handle<double> fixedGridRhoFastjetAllCalo;
  iEvent.getByToken(fixedGridRhoFastjetAllCaloToken_,fixedGridRhoFastjetAllCalo);
  ev.fGridRhoFastjetAllCalo = *fixedGridRhoFastjetAllCalo;

  edm::Handle<double> fixedGridRhoFastjetCentralCalo;
  iEvent.getByToken(fixedGridRhoFastjetCentralCaloToken_,fixedGridRhoFastjetCentralCalo);
  ev.fGridRhoFastjetCentralCalo = *fixedGridRhoFastjetCentralCalo;

  edm::Handle<double> fixedGridRhoFastjetCentralChargedPileUp;
  iEvent.getByToken(fixedGridRhoFastjetCentralChargedPileUpToken_,fixedGridRhoFastjetCentralChargedPileUp);
  ev.fGridRhoFastjetCentralChargedPileUp = *fixedGridRhoFastjetCentralChargedPileUp;

  edm::Handle<double> fixedGridRhoFastjetCentralNeutral;
  iEvent.getByToken(fixedGridRhoFastjetCentralNeutralToken_,fixedGridRhoFastjetCentralNeutral);
  ev.fGridRhoFastjetCentralNeutral = *fixedGridRhoFastjetCentralNeutral;

  // Good Primary Vertex Part
  edm::Handle<reco::VertexCollection> primaryVertices;
  bool found = iEvent.getByToken(vertexToken_, primaryVertices);

  if (found && primaryVertices.isValid()) {
    edm::LogInfo("EventBlock") << "Total # Primary Vertices: " << primaryVertices->size();
    for (const auto& v: *primaryVertices) {
      if (!v.isFake() &&
           v.ndof() > vtxMinNDOF_ &&
	   std::abs(v.z()) <= vtxMaxAbsZ_ &&
	   std::abs(v.position().rho()) <= vtxMaxd0_)
      {
        ev.hasPrimaryVertex = true;
        break;
      }
    }
  }
  else {
    edm::LogError("EventBlock") << "Error! Failed to get VertexCollection for label: "
                                << vertexTag_;
  }

  // Access PU information
  if (!iEvent.isRealData()) {
    edm::Handle<std::vector<PileupSummaryInfo> > puInfo;
    found = iEvent.getByToken(puSummaryToken_, puInfo);
    if (found && puInfo.isValid()) {
      for (const auto& v: *puInfo) {
	ev.bunchCrossing.push_back(v.getBunchCrossing()); // in-time and out-of-time indices
	ev.nPU.push_back(v.getPU_NumInteractions());
	ev.trueNInt.push_back(v.getTrueNumInteractions());
      }
    }
    // More info about PU is here:
    // https://twiki.cern.ch/twiki/bin/viewauth/CMS/PileupInformation#Accessing_PileupSummaryInfo_in_r
  }
  // Access rhoNeutral
  edm::Handle<double> rhoNeutral;
  found = iEvent.getByToken(rhoNeutralToken_, rhoNeutral);
  if (found) ev.rhoNeutral = *rhoNeutral;
  
  // Vertex Container
  edm::Handle<reco::VertexCollection> spVertices;
  found = iEvent.getByToken(selectedVertexToken_, spVertices);
  if (found) ev.nvtx = spVertices->size();

  list_->push_back(ev);
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
EventBlock::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}
#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(EventBlock);
