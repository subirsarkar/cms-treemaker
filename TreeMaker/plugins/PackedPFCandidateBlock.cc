#include <iostream>
#include <algorithm>

#include "TTree.h"
#include "TPRegexp.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/Framework/interface/ESHandle.h"

#include "DataFormats/GeometryVector/interface/GlobalVector.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "DataFormats/Math/interface/deltaR.h"

#include "AnalysisSpace/TreeMaker/plugins/PackedPFCandidateBlock.h"
#include "AnalysisSpace/TreeMaker/interface/Utility.h"

using std::setprecision;
using std::setw;
using std::ios;
using std::setiosflags;
using std::cout;
using std::endl;

// Constructor
PackedPFCandidateBlock::PackedPFCandidateBlock(const edm::ParameterSet& iConfig) :
  verbosity_(iConfig.getUntrackedParameter<int>("verbosity", 0)),
  skimCand_(iConfig.getUntrackedParameter<bool>("skimCandidate", false)),
  pdgToSaveList_(iConfig.getParameter<std::vector<int>>("pdgTosave")),
  minCandPt_(iConfig.getUntrackedParameter<double>("minCandPt", 2.0)),
  pfCandTag_(iConfig.getParameter<edm::InputTag>("pfCands")), 
  vertexTag_(iConfig.getParameter<edm::InputTag>("vertexSrc")),
  pfToken_(consumes<pat::PackedCandidateCollection>(pfCandTag_)),
  vertexToken_(consumes<reco::VertexCollection>(vertexTag_))
{
}
PackedPFCandidateBlock::~PackedPFCandidateBlock() {
  delete list_;
}
void PackedPFCandidateBlock::beginJob()
{
  // Get TTree pointer
  std::string tree_name = "vhtree";
  TTree* tree = vhtm::Utility::getTree(tree_name);
  list_ = new std::vector<vhtm::PackedPFCandidate>();
  tree->Branch("PackedPFCandidate", "std::vector<vhtm::PackedPFCandidate>", &list_, 32000, -1);
  tree->Branch("nPackedPFCandidate", &fnPackedPFCandidate_, "fnPackedPFCandidate_/I");
}
void PackedPFCandidateBlock::analyze(edm::Event const& iEvent, edm::EventSetup const& iSetup) 
{
  // Reset the vector and the nObj variables
  list_->clear();
  fnPackedPFCandidate_ = 0;

  edm::Handle<pat::PackedCandidateCollection> pfs;
  bool found = iEvent.getByToken(pfToken_, pfs);
  
  if (found && pfs.isValid()) {
    edm::LogInfo("PackedPFCandidateBlock") << "Total # PackedPFCandidate: " << pfs->size();

    edm::Handle<reco::VertexCollection> primaryVertices;
    iEvent.getByToken(vertexToken_, primaryVertices);

    if (verbosity_)
      cout << "   indx  pdgId     pt    eta    phi  energy   mass charge fromPV hPurity" << endl;
    for (auto const& v: *pfs) {
      if (list_->size() == kMaxPackedPFCandidate) {
        edm::LogInfo("PackedPFCandidateBlock") << "Too many PackedPFCandidates, fnPackedPFCandidate = " 
                                               << fnPackedPFCandidate_; 
        break;
      }
      
      int pdgId = std::abs(v.pdgId());
      if ( v.pt() <= minCandPt_ ||
	   (skimCand_ && std::find(pdgToSaveList_.begin(), pdgToSaveList_.end(), pdgId) == pdgToSaveList_.end()) ) continue;
      
      vhtm::PackedPFCandidate pfCand;
      
      pfCand.pt = v.pt();
      pfCand.eta = v.eta();
      pfCand.phi = v.phi();
      pfCand.phiAtVtx = v.phiAtVtx();
      pfCand.energy = v.energy();
      pfCand.trackHighPurity = v.trackHighPurity();
      
      pfCand.pdgId = pdgId;
      pfCand.charge = v.charge();
      
      pfCand.vx = v.vx();
      pfCand.vz = v.vz();
      pfCand.vz = v.vz();
      
      pfCand.fromPV = v.fromPV();
      pfCand.dxy = v.dxy();
      pfCand.dz = v.dz();
      pfCand.dxyError = v.dxyError();  
      pfCand.dzError = v.dzError();   
      
      // dxy & dz wrt the Event vertex
      const reco::Vertex& vit = primaryVertices->front(); // Highest sumPt vertex
      pfCand.dxyEV = v.dxy(vit.position());
      pfCand.dzEV  = v.dz(vit.position());
      pfCand.dzAssociatedPV = v.dzAssociatedPV();
      
      pfCand.numberOfHits = v.numberOfHits();
      pfCand.numberOfPixelHits = v.numberOfPixelHits();
      pfCand.pixelLayersWithMeasurement = v.pixelLayersWithMeasurement();
      pfCand.stripLayersWithMeasurement = v.stripLayersWithMeasurement();
      pfCand.lostInnerHits = v.lostInnerHits();

      std::vector<double> isotemp;   
      calcIsoFromPF(v, pfs, 0.3, isotemp);
      pfCand.isolationMap["c30"] = isotemp;
      
      isotemp.clear();   
      calcIsoFromPF(v, pfs, 0.35, isotemp);
      pfCand.isolationMap["c35"] = isotemp;

      isotemp.clear();   
      calcIsoFromPF(v, pfs, 0.4, isotemp);
      pfCand.isolationMap["c40"] = isotemp;

      if (verbosity_) {
	cout.precision(2);
        cout << std::fixed;
	cout << setw(7) << list_->size()
	     << setw(7) << pdgId
	     << setw(7) << v.pt() 
	     << setw(7) << v.eta()
	     << setw(7) << v.phi()
	     << setw(8) << v.energy()
             << setprecision(3)
	     << setw(7) << v.mass()
             << setprecision(2)
	     << setw(7) << v.charge()
	     << setw(7) << v.fromPV()
	     << setw(8) << v.trackHighPurity()
	     << endl;
      }
      list_->push_back(pfCand);
    }
    fnPackedPFCandidate_ = list_->size();
  }
  else {
    edm::LogError("PackedPFCandidateBlock") << "Error >> Failed to get pat::PackedPFCandidate for label: " 
                                            << pfCandTag_;
  }
}
void PackedPFCandidateBlock::calcIsoFromPF(const pat::PackedCandidate& v, 
                                           edm::Handle<pat::PackedCandidateCollection>& pfs, 
                                           double cone, std::vector<double>& iso)
{
  // initialize sums
  double chargedHadSum = 0., 
    chargedParticleSum = 0., 
    neutralSum = 0., 
    photonSum = 0., 
    pileupSum  = 0;

  // now get a list of the PF candidates used to build this lepton, so to exclude them
  std::vector<reco::CandidatePtr> footprint;
  for (unsigned int i = 0; i < v.numberOfSourceCandidatePtrs(); ++i)
    footprint.push_back(v.sourceCandidatePtr(i));
 
  // now loop on pf candidates
  for (unsigned int i = 0; i < pfs->size(); ++i) {
    const pat::PackedCandidate& pf = (*pfs)[i];
    double dRcone = deltaR(v, pf);
    if (dRcone < cone) {
      // pfcandidate-based footprint removal
      if (std::find(footprint.begin(), footprint.end(), reco::CandidatePtr(pfs, i)) != footprint.end()) continue;

      double pt = pf.pt(); 
      int pdg = std::abs(pf.pdgId());
      if (pf.charge() == 0) {
        if (pt > 0.5 && dRcone > 0.01) {
          if (pdg == 22)
            photonSum += pt;
          else 
            neutralSum += pt;
        }
      } 
      else {
        if (pt > 0.2 && dRcone > 0.0001) { 
          if (pf.vertexRef().isNonnull() &&  pf.fromPV() >= 2) {
            chargedParticleSum += pt;
            if (pdg != 13 && pdg != 11) chargedHadSum += pt;
          } 
          else 
            pileupSum += pt;
        }
      }
    }
  }
  if (verbosity_ > 1) cout << "isoValues: (" << chargedHadSum << "," 
			   << neutralSum << "," << photonSum << "," 
			   << pileupSum << ")" 
			   << endl;
  iso.push_back(chargedHadSum);
  iso.push_back(chargedParticleSum);
  iso.push_back(neutralSum);
  iso.push_back(photonSum);
  iso.push_back(pileupSum);
}
// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
PackedPFCandidateBlock::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  // The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}
DEFINE_FWK_MODULE(PackedPFCandidateBlock);
