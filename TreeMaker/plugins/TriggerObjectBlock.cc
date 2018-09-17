#include <iostream>
#include <algorithm>

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "PhysicsTools/PatUtils/interface/TriggerHelper.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "FWCore/Common/interface/TriggerNames.h"

#include "AnalysisSpace/TreeMaker/plugins/TriggerObjectBlock.h"
#include "AnalysisSpace/TreeMaker/interface/Utility.h"
#include "AnalysisSpace/TreeMaker/interface/PhysicsObjects.h"

#include "TMath.h"
#include "TTree.h"
#include "TPRegexp.h"

// Constructor
TriggerObjectBlock::TriggerObjectBlock(const edm::ParameterSet& iConfig) :
  verbosity_(iConfig.getUntrackedParameter<int>("verbosity", 0)),
  minTrigObjPt_(iConfig.getUntrackedParameter<double>("minTrigObjPt", 5.0)),
  hltTag_(iConfig.getParameter<edm::InputTag>("hltTag")),
  objectTag_(iConfig.getParameter<edm::InputTag>("triggerObjectTag")),
  hltToken_(consumes<edm::TriggerResults>(hltTag_)),
  objectToken_(consumes<pat::TriggerObjectStandAloneCollection>(objectTag_))
{
}
TriggerObjectBlock::~TriggerObjectBlock() {
  delete list_;
}
void TriggerObjectBlock::beginJob()
{
  std::string tree_name = "vhtree";
  TTree* tree = vhtm::Utility::getTree(tree_name);

  list_ = new std::vector<vhtm::TriggerObject>();
  tree->Branch("TriggerObject", "std::vector<vhtm::TriggerObject>", &list_, 32000, -1);
  tree->Branch("nTriggerObject", &fnTriggerObject_, "fnTriggerObject_/I");
}
void TriggerObjectBlock::beginRun(edm::Run const& iRun, edm::EventSetup const& iSetup) {
  bool changed = true;
  if (hltConfig_.init(iRun, iSetup, hltTag_.process(), changed)) {
    // if init returns TRUE, initialisation has succeeded!
    edm::LogInfo("TriggerObjectBlock") << "HLT config with process name "
				       << hltTag_.process()
				       << " successfully extracted";
  }
  else {
    // if init returns FALSE, initialisation has NOT succeeded, which indicates a problem
    // with the file and/or code and needs to be investigated!
    edm::LogError("TriggerObjectBlock") << "Error! HLT config extraction with process name "
					<< hltTag_.process() << " failed";
    // In this case, all access methods will return empty values!
  }
}
void TriggerObjectBlock::analyze(edm::Event const& iEvent, edm::EventSetup const& iSetup) {
  // Reset the vector and the nObj variables
  list_->clear();
  fnTriggerObject_ = 0;

  if (verbosity_ > 1) {
    edm::LogInfo("TriggerObjectBlock") << setiosflags(std::ios::fixed);
    edm::LogInfo("TriggerObjectBlock") << "Indx Eta Phi Pt Energy =Trigger path list=";
  }
 
  edm::Handle<edm::TriggerResults> triggerBits;
  bool found = iEvent.getByToken(hltToken_, triggerBits);
  if (found && triggerBits.isValid()) {
    edm::Handle<pat::TriggerObjectStandAloneCollection> triggerObjects;
    found = iEvent.getByToken(objectToken_, triggerObjects);
    if (found && triggerObjects.isValid()) {
      // Find the triggerNames
      const auto& names = iEvent.triggerNames(*triggerBits);
      for (auto obj: *triggerObjects) {
	if (list_->size() == kMaxTriggerObject_) {
	  edm::LogInfo("TriggerObjectBlock") << "Too many Trigger Objects (HLT), fnTriggerObject = "
					     << list_->size();
	  break;
	}
	if (obj.pt() < minTrigObjPt_) continue;
	obj.unpackPathNames(names);
	
        if (verbosity_) TriggerObjectBlock::printObjectInfo(obj);

	std::map<std::string, unsigned int> pathInfoMap;
	for (const auto& v: obj.pathNames(false)) {
	  int val = -1;
	  if      (obj.hasPathName(v, true,  true))  val = 3; 
	  else if (obj.hasPathName(v, false, true))  val = 2; 
	  else if (obj.hasPathName(v, true,  false)) val = 1; 
	  else if (obj.hasPathName(v, false, false)) val = 0; 
	  if (val > -1) pathInfoMap.insert({v, val});
	}
	
	vhtm::TriggerObject _tobj;
	_tobj.eta      = obj.eta();
	_tobj.phi      = obj.phi();
	_tobj.pt       = obj.pt();
	_tobj.energy   = obj.energy();
	_tobj.pathList = pathInfoMap;
	list_->push_back(_tobj);
	
	if (verbosity_ > 1) {
	  edm::LogInfo("TriggerObjectBlock") << std::setprecision(2);
	  edm::LogInfo("TriggerObjectBlock") << std::setw(4) << list_->size()
					     << std::setw(8) << _tobj.eta
					     << std::setw(8) << _tobj.phi
					     << std::setw(8) << _tobj.pt
					     << std::setw(8) << _tobj.energy;
	  for (const auto& jt: _tobj.pathList)
	    edm::LogInfo("TriggerObjectBlock") << "\t\t\t\t\t" << jt.first << " " << jt.second;;
	}
      }
      fnTriggerObject_ = list_->size();
    }
    else
      edm::LogError("TriggerObjectBlock") << "Failed to get TriggerObjects for label: "
					  << objectTag_;
  }
  else
    edm::LogError("TriggerObjectBlock") << "Failed to get TriggerResults for label: "
					<< hltTag_;
}
void TriggerObjectBlock::printObjectInfo(const pat::TriggerObjectStandAlone& obj) 
{
  edm::LogInfo("TriggerObjectBlock") << "\tTrigger object:  pt " << obj.pt() 
				     << ", eta " << obj.eta() 
				     << ", phi " << obj.phi();
  // Print trigger object collection and type
  edm::LogInfo("TriggerObjectBlock") << "\t   Collection: " << obj.collection() << "\n"
				     << "\t   Type IDs:   ";
#if 0
  for (unsigned h = 0; h < obj.filterIds().size(); ++h) 
    edm::LogInfo("TriggerObjectBlock") << " " << obj.filterIds()[h] ;
#endif
  for (const auto& v: obj.filterIds()) 
    edm::LogInfo("TriggerObjectBlock") << " " << v;

  // Print associated trigger filters
  edm::LogInfo("TriggerObjectBlock") << "\t   Filters:    ";
#if 0
  for (unsigned h = 0; h < obj.filterLabels().size(); ++h) 
    edm::LogInfo("TriggerObjectBlock") << " " << obj.filterLabels()[h];
#endif
  for (const auto& v: obj.filterLabels()) 
    edm::LogInfo("TriggerObjectBlock") << " " << v;

  const auto& pathNamesAll  = obj.pathNames(false);
  const auto& pathNamesLast = obj.pathNames(true);

  // Print all trigger paths, for each one record also if the object is associated to a 'l3' filter (always true for the
  // definition used in the PAT trigger producer) and if it's associated to the last filter of a successfull path (which
  // means that this object did cause this trigger to succeed; however, it doesn't work on some multi-object triggers)
  edm::LogInfo("TriggerObjectBlock") << "\t   Paths (" 
				     << pathNamesAll.size()
				     << "/"
				     << pathNamesLast.size()
				     << "):    ";
  for (const auto& p: pathNamesAll) {
    bool isBoth = obj.hasPathName(p, true,  true); 
    bool isL3   = obj.hasPathName(p, false, true); 
    bool isLF   = obj.hasPathName(p, true,  false); 
    bool isNone = obj.hasPathName(p, false, false); 

    edm::LogInfo("TriggerObjectBlock") << "   " << p;
    if (isBoth)                              edm::LogInfo("TriggerObjectBlock") << "(L,3)";
    if (isL3 && !isBoth)                     edm::LogInfo("TriggerObjectBlock") << "(*,3)";
    if (isLF && !isBoth)                     edm::LogInfo("TriggerObjectBlock") << "(L,*)";
    if (isNone && !isBoth && !isL3 && !isLF) edm::LogInfo("TriggerObjectBlock") << "(*,*)";
  }
}
#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(TriggerObjectBlock);
