#include "TTree.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "AnalysisSpace/TreeMaker/interface/PhysicsObjects.h"
#include "AnalysisSpace/TreeMaker/plugins/METBlock.h"
#include "AnalysisSpace/TreeMaker/interface/Utility.h"

METBlock::METBlock(const edm::ParameterSet& iConfig) :
  verbosity_(iConfig.getUntrackedParameter<int>("verbosity", 0)),
  pfMETTag_(iConfig.getParameter<edm::InputTag>("metSrc")),
  corrMETTag_(iConfig.getParameter<edm::InputTag>("corrmetSrc")),
  puppiMETTag_(iConfig.getParameter<edm::InputTag>("puppimetSrc")),
  pfMETToken_(consumes<pat::METCollection>(pfMETTag_)),
  corrMETToken_(consumes<pat::METCollection>(corrMETTag_)),
  puppiMETToken_(consumes<pat::METCollection>(puppiMETTag_))
{
}
void METBlock::beginJob()
{
  // Get TTree pointer
  TTree* tree = vhtm::Utility::getTree("vhtree");

  pfList_ = new std::vector<vhtm::MET>();
  tree->Branch("MET", "std::vector<vhtm::MET>", &pfList_, 32000, -1);
  tree->Branch("nMET", &fnPFMET_, "fnPFMET_/I");

  corrList_ = new std::vector<vhtm::MET>();
  tree->Branch("corrMET", "std::vector<vhtm::MET>", &corrList_, 32000, -1);
  tree->Branch("corrnMET", &fnCorrMET_, "fnCorrMET_/I");

  puppiList_ = new std::vector<vhtm::MET>();
  tree->Branch("puppiMET", "std::vector<vhtm::MET>", &puppiList_, 32000, -1);
  tree->Branch("puppinMET", &fnPuppiMET_, "fnPuppiMET_/I");
}
void METBlock::analyze(edm::Event const& iEvent, edm::EventSetup const& iSetup) {
  fillMET(iEvent, iSetup, pfList_, fnPFMET_, pfMETTag_, pfMETToken_);
  fillMET(iEvent, iSetup, corrList_, fnCorrMET_, corrMETTag_, corrMETToken_);
  fillMET(iEvent, iSetup, puppiList_, fnPuppiMET_, puppiMETTag_, puppiMETToken_);
}
void METBlock::fillMET(const edm::Event& iEvent,
                       const edm::EventSetup& iSetup,
                       std::vector<vhtm::MET>* list,
                       int& nMET,
                       const edm::InputTag& iTag,
                       const edm::EDGetTokenT<pat::METCollection>& token)
{
  // Reset the TClonesArray and the nObj variables
  list->clear();
  nMET = 0;

  edm::Handle<pat::METCollection> metColl;
  bool found = iEvent.getByToken(token, metColl);

  if (found && metColl.isValid()) {
    edm::LogInfo("METBlock") << "Total # PAT METs: " << metColl->size();
    for (const auto& v: *metColl) {
      if (list->size() == kMaxMET_) {
        edm::LogInfo("METBlock") << "Too many PAT MET, nMET = " << list->size()
				 << ", label: " << iTag;
        break;
      }
      vhtm::MET mobj;
      // fill in all the vectors
      mobj.met          = v.pt();
      mobj.metphi       = v.phi();
      mobj.sumet        = v.sumEt();
#if 0
      mobj.metuncorr    = v.uncorrectedPt(pat::MET::uncorrALL);
      mobj.metphiuncorr = v.uncorrectedPhi(pat::MET::uncorrALL);
      mobj.sumetuncorr  = v.sumEt() - v.corSumEt(pat::MET::uncorrALL);
#endif
      mobj.metJESUp     = v.shiftedPt(pat::MET::JetEnUp);
      mobj.metJESDn     = v.shiftedPt(pat::MET::JetEnDown);

      mobj.metuncorr    = v.uncorPt();
      mobj.metphiuncorr = v.uncorPhi();
      mobj.sumetuncorr  = v.uncorSumEt();

      list->push_back(mobj);
    }
    nMET = list->size();      
  }
  else {
    edm::LogError("METBlock") << "Error! Failed to get pat::MET collection for label: "
                              << iTag;
  }
}
#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(METBlock);
