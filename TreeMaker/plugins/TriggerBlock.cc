#include <iostream>
#include <algorithm>

#include "TTree.h"
#include "TPRegexp.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutRecord.h"

#include "AnalysisSpace/TreeMaker/plugins/TriggerBlock.h"
#include "AnalysisSpace/TreeMaker/interface/Utility.h"

static const unsigned int NmaxL1AlgoBit = 128;
static const unsigned int NmaxL1TechBit = 64;

// Constructor
TriggerBlock::TriggerBlock(const edm::ParameterSet& iConfig) :
  verbosity_(iConfig.getUntrackedParameter<int>("verbosity", 0)),
  hltPathsOfInterest_(iConfig.getParameter<std::vector<std::string>>("hltPathsOfInterest")),
  l1Tag_(iConfig.getParameter<edm::InputTag>("l1Tag")),
  hltTag_(iConfig.getParameter<edm::InputTag>("hltTag")),
  l1Token_(consumes<L1GlobalTriggerReadoutRecord>(l1Tag_)),
  hltToken_(consumes<edm::TriggerResults>(hltTag_)),
  hltPrescaleProvider_(iConfig, consumesCollector(), *this)
{
}
TriggerBlock::~TriggerBlock() {
  delete l1physbits_;
  delete l1techbits_;
  delete hltpaths_;
  delete hltresults_;
  delete hltprescales_;
}
void TriggerBlock::beginJob()
{
  std::string tree_name = "vhtree";
  TTree* tree = vhtm::Utility::getTree(tree_name);

  l1physbits_ = new std::vector<int>();
  tree->Branch("l1physbits", "vector<int>", &l1physbits_);

  l1techbits_ = new std::vector<int>();
  tree->Branch("l1techbits", "vector<int>", &l1techbits_);

  hltpaths_ = new std::vector<std::string>();
  tree->Branch("hltpaths", "vector<string>", &hltpaths_);

  hltresults_ = new std::vector<int>();
  tree->Branch("hltresults", "vector<int>", &hltresults_);

  hltprescales_ = new std::vector<int>();
  tree->Branch("hltprescales", "vector<int>", &hltprescales_);
}
void TriggerBlock::beginRun(edm::Run const& iRun, edm::EventSetup const& iSetup) {
  bool changed(true);
  if (hltPrescaleProvider_.init(iRun, iSetup, hltTag_.process(), changed)) {
    HLTConfigProvider const& hltConfig = hltPrescaleProvider_.hltConfigProvider();
    if (changed) {
      // if init returns TRUE, initialisation has succeeded!
      if (verbosity_)
      std::cout << "TriggerBlock: " << "HLT config with process name "
		<< hltTag_.process() 
		<< " successfully extracted"
		<< std::endl;
      matchedPathList_.clear();
      const auto& pathList = hltConfig.triggerNames();
      for (const auto& path: pathList) {
	if (hltPathsOfInterest_.size()) {
	  int nmatch = 0;
	  for (const auto& kt: hltPathsOfInterest_)
	    nmatch += TPRegexp(kt).Match(path);
	  if (!nmatch) continue;
	}
	matchedPathList_.push_back(path);
      }
    }
  }
  else {
    // if init returns FALSE, initialisation has NOT succeeded, which indicates a problem
    // with the file and/or code and needs to be investigated!
    edm::LogError("TriggerBlock") << "Error! HLT config extraction with process name "
				  << hltTag_.process() << " failed";
    // In this case, all access methods will return empty values!
  }
  //bool isConfigChanged = false;
  //hltPrescaleProvider_.init(iRun, iSetup, "HLT", isConfigChanged);
}
void TriggerBlock::analyze(edm::Event const& iEvent, edm::EventSetup const& iSetup) { 
  // Reset the vectors
  l1physbits_->clear();
  l1techbits_->clear();
  hltpaths_->clear();
  hltresults_->clear();
  hltprescales_->clear();

#if 0
  edm::Handle<L1GlobalTriggerReadoutRecord> l1GtReadoutRecord;
  bool found = iEvent.getByToken(l1Token_, l1GtReadoutRecord);
  if (found && l1GtReadoutRecord.isValid()) {
    edm::LogInfo("TriggerBlock") << "Successfully obtained L1GlobalTriggerReadoutRecord for label: "
                                 << l1Tag_;

    for (unsigned int i = 0; i < NmaxL1AlgoBit; ++i) 
      l1physbits_->push_back(l1GtReadoutRecord->decisionWord()[i] ? 1 : 0);

    for (unsigned int i = 0; i < NmaxL1TechBit; ++i) 
      l1techbits_->push_back(l1GtReadoutRecord->technicalTriggerWord()[i] ? 1 : 0 );
  }
  else 
    edm::LogError("TriggerBlock") << "Error >> Failed to get L1GlobalTriggerReadoutRecord for label: "
                                  << l1Tag_;
#endif

  edm::Handle<edm::TriggerResults> triggerResults;
  bool hltfound = iEvent.getByToken(hltToken_, triggerResults);
  if (hltfound && triggerResults.isValid()) {
    if (verbosity_) 
      std::cout << "TriggerBlock: Successfully obtained edm::TriggerResults for tag: " << hltTag_
		<< std::endl;

    HLTConfigProvider const& hltConfig = hltPrescaleProvider_.hltConfigProvider();
    for (const auto& path: matchedPathList_) {
      hltpaths_->push_back(path);

      int fired = -1;
      unsigned int index = hltConfig.triggerIndex(path);
      if (index < triggerResults->size())
        fired = (triggerResults->accept(index)) ? 1 : 0;
      else
	std::cerr << "TriggerBlock: Requested HLT path \"" << path << "\" does not exist";
      
      hltresults_->push_back(fired);

      int prescale = -1;
      const int prescaleSet = hltPrescaleProvider_.prescaleSet(iEvent, iSetup);
      if (prescaleSet < 0) {
        if (verbosity_)
	std::cout << "TriggerBlock: The prescale set index number could not be obtained for HLT path: " << path
		  << std::endl;
      }
      else {
        prescale = hltConfig.prescaleValue(prescaleSet, path);
      }
      hltprescales_->push_back(prescale);

      if (verbosity_) {
        auto a = hltPrescaleProvider_.prescaleValues(iEvent, iSetup, path);

	std::cout << "TriggerBlock: Path: " << path
		  << ", prescale: " << prescale 
		  << ", fired: " << fired
		  << ", PrescaleValues L1: " << a.first 
		  << ", PrescaleValues HLT: " << a.second
		  << std::endl;
        const auto& d = hltConfig.moduleLabels(path);
        for (const auto& v: d)
	  std::cout << "TriggerBlock:\tModule Labels: " << v << std::endl;
        
      }    
    } 
    if (verbosity_) {
      const auto& b = hltConfig.prescaleLabels();
      for (const auto& v: b)
	std::cout << "TriggerBlock:\tPrescale Labels: " << v
		  << std::endl;

      const auto& c = hltConfig.prescaleTable();
      for (const auto& ptr: c) {
	std::cout << "TriggerBlock: Key  : " << ptr.first << ": " << std::endl;
        for (const auto& v: ptr.second)
	  std::cout << "TriggerBlock value: " << v << std::endl;
      }
    } 
  } 
  else {
    std::cerr << "TriggerBlock: Failed to get TriggerResults for label: " << hltTag_ << std::endl;
  }
}
// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
TriggerBlock::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}
DEFINE_FWK_MODULE(TriggerBlock);
