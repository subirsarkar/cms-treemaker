#include <iostream>
#include<iomanip>
#include "TTree.h"
#include "TClass.h"
using std::setw;

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/TrackReco/interface/HitPattern.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "Math/GenVector/VectorUtil.h"
#include "DataFormats/Math/interface/deltaR.h"

#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "DataFormats/GeometryVector/interface/VectorUtil.h"

#include "DataFormats/MuonReco/interface/MuonIsolation.h"
#include "DataFormats/MuonReco/interface/MuonPFIsolation.h"

#include "DataFormats/MuonReco/interface/MuonSelectors.h"

#include "AnalysisSpace/TreeMaker/interface/PhysicsObjects.h"
#include "AnalysisSpace/TreeMaker/plugins/MuonBlock.h"
#include "AnalysisSpace/TreeMaker/interface/Utility.h"

MuonBlock::MuonBlock(const edm::ParameterSet& iConfig):
  verbosity_(iConfig.getUntrackedParameter<int>("verbosity", 0)),
  bsCorr_(iConfig.getUntrackedParameter<bool>("beamSpotCorr", true)),
  muonID_(iConfig.getUntrackedParameter<std::string>("muonID", "GlobalMuonPromptTight")),
  muonTag_(iConfig.getParameter<edm::InputTag>("muonSrc")),
  vertexTag_(iConfig.getParameter<edm::InputTag>("vertexSrc")),
  bsTag_(iConfig.getParameter<edm::InputTag>("offlineBeamSpot")),
  pfcandTag_(iConfig.getParameter<edm::InputTag>("pfCands")),
  muonToken_(consumes<pat::MuonCollection>(muonTag_)),
  vertexToken_(consumes<reco::VertexCollection>(vertexTag_)),
  bsToken_(consumes<reco::BeamSpot>(bsTag_)),
  pfToken_(consumes<pat::PackedCandidateCollection>(pfcandTag_)),
  defaultBestMuon_(!iConfig.existsAs<std::string>("customArbitration")),
  bestMuonSelector_(defaultBestMuon_ ? std::string("") : iConfig.getParameter<std::string>("customArbitration"))
{
}
MuonBlock::~MuonBlock() {
}
void MuonBlock::beginJob()
{
  // Get TTree pointer
  TTree* tree = vhtm::Utility::getTree("vhtree");
  list_ = new std::vector<vhtm::Muon>();
  tree->Branch("Muon", "std::vector<vhtm::Muon>", &list_, 32000, -1);
  tree->Branch("nMuon", &fnMuon_, "fnMuon_/I");
}
void MuonBlock::analyze(edm::Event const& iEvent, edm::EventSetup const& iSetup) {
  // Reset the vector and the nObj variables
  list_->clear();
  fnMuon_ = 0;

  edm::Handle<pat::MuonCollection> muons;
  bool found = iEvent.getByToken(muonToken_, muons);

  edm::Handle<pat::PackedCandidateCollection> pfs;
  iEvent.getByToken(pfToken_, pfs);

  if (found && muons.isValid()) {
    edm::Handle<reco::VertexCollection> primaryVertices;
    iEvent.getByToken(vertexToken_, primaryVertices);

    edm::Handle<reco::BeamSpot> beamSpot;
    iEvent.getByToken(bsToken_, beamSpot);

    edm::LogInfo("MuonBlock") << "Total # of Muons: " << muons->size();

    // Muon Ghost cleaning
    unsigned int nMu = muons->size();
    std::vector<bool> good(nMu, true);
    for (unsigned int i = 0; i < nMu; ++i) {
      const pat::Muon& mu1 = muons->at(i);
      if (!mu1.track().isNonnull()) good[i] = false; 
      if (!good[i]) continue;
      int nSegments1 = mu1.numberOfMatches(reco::Muon::SegmentArbitration);
      for (unsigned int j = i+1; j < nMu; ++j) {
        const pat::Muon& mu2 = muons->at(j);
        if (isSameMuon(mu1, mu2)) continue;
        if (!good[j] || !mu2.track().isNonnull()) continue;
        int nSegments2 = mu2.numberOfMatches(reco::Muon::SegmentArbitration);
        if (nSegments2 == 0 || nSegments1 == 0) continue;
        double sf = muon::sharedSegments(mu1, mu2)/std::min<double>(nSegments1, nSegments2);
	if (sf > 0.499) {
	  if (isBetterMuon(mu1, mu2)) {
	    good[j] = false;
	  } else {
	    good[i] = false;
	  }
	}
      }
    }

    for ( unsigned int i = 0; i < nMu; ++i ) {
      if (list_->size() == kMaxMuon_) {
	edm::LogInfo("MuonBlock") << "Too many PAT Muons, fnMuon = " << list_->size();
	break;
      }
      const pat::Muon& v = muons->at(i);
      bool ghostCleaned = good[i] || (v.isGlobalMuon() && v.numberOfMatches() >= 2);
      vhtm::Muon muon;
      muon.isGlobalMuon  = v.isGlobalMuon() ? true : false;
      muon.isTrackerMuon = v.isTrackerMuon() ? true : false;
      muon.isPFMuon      = v.isPFMuon();
    
      muon.isghostCleaned = ghostCleaned;

      muon.eta     = v.eta();
      muon.phi     = v.phi();
      muon.pt      = v.pt();
      muon.p       = v.p();
      muon.energy  = v.energy();
      muon.charge  = v.charge();
      
      reco::TrackRef tk = v.muonBestTrack();
      bool hasTkinRef = tk.isNonnull();
       
      if (hasTkinRef) {
        muon.tkNChi2 = tk->normalizedChi2(); 
        double trkd0 = tk->d0();
        double trkdz = tk->dz();
        muon.muonBestTrackType = v.muonBestTrackType();
 
        if (bsCorr_) {
          if (beamSpot.isValid()) {
            trkd0 = -(tk->dxy(beamSpot->position()));
            trkdz = tk->dz(beamSpot->position());
          }
          else
            edm::LogError("MuonsBlock") << "Error >> Failed to get reco::BeamSpot for label: "
                                        << bsTag_;
	}
	muon.trkD0 = trkd0;
	muon.trkDz = trkdz;
      }  

      double globalChi2 = v.isGlobalMuon() ? v.normChi2() : 9999.;
      muon.globalChi2 = globalChi2;
      muon.passID     = v.muonID(muonID_) ? true : false;

      muon.chi2LocalPosition = v.combinedQuality().chi2LocalPosition;
      muon.trkKink = v.combinedQuality().trkKink;
      muon.isLooseMuon  = v.isLooseMuon();
      muon.isMediumMuon = v.isMediumMuon();

      bool goodGlb = v.isGlobalMuon() &&
                     v.normChi2() < 3 &&
                     v.combinedQuality().chi2LocalPosition < 12 &&
                     v.combinedQuality().trkKink < 20;
      bool goodMedMuon = v.innerTrack()->validFraction() >= 0.8 &&
	                 v.segmentCompatibility() >= (goodGlb ? 0.303 : 0.451);
      muon.isGoodMedMuon = goodMedMuon;

      double dxyWrtPV = -99.;
      double dzWrtPV = -99.;
      if (primaryVertices.isValid()) {
        edm::LogInfo("MuonBlock") << "Total # Primary Vertices: " << primaryVertices->size();

        const reco::Vertex& vit = primaryVertices->front(); // Highest sumPt vertex
        muon.passTrackerhighPtid = isTrackerHighPt(v, vit);
        muon.isTightMuon = v.isTightMuon(vit);

        dxyWrtPV = tk->dxy(vit.position());
        dzWrtPV  = tk->dz(vit.position());
     
        muon.dxyPV = dxyWrtPV;
        muon.dzPV  = dzWrtPV;

        // Vertex association
        double minVtxDist3D = 9999.;
        int indexVtx = -1;
        double vertexDistZ = 9999.;
        for (auto vit = primaryVertices->begin(); vit != primaryVertices->end(); ++vit) {
          double dxy = 9999.;
          double dz = 9999.;
          if (hasTkinRef) {
            dxy = tk->dxy(vit->position());
            dz = tk->dz(vit->position());
            double dist3D = std::sqrt(pow(dxy,2) + pow(dz,2));
            if (dist3D < minVtxDist3D) {
              minVtxDist3D = dist3D;
              indexVtx = static_cast<int>(std::distance(primaryVertices->begin(), vit));
              vertexDistZ = dz;
            }
          }
        }

        muon.vtxDist3D = minVtxDist3D;
        muon.vtxIndex = indexVtx;
        muon.vtxDistZ = vertexDistZ;
      }
      else {
        edm::LogError("MuonBlock") << "Error >> Failed to get reco::VertexCollection for label: "
                                   << vertexTag_;
      }

      // Hit pattern
      if (hasTkinRef) {
	const reco::HitPattern& hitp = tk->hitPattern(); // or global track
	muon.pixHits = hitp.numberOfValidPixelHits();
	muon.trkHits = hitp.numberOfValidTrackerHits();
	muon.muoHits = hitp.numberOfValidMuonHits();
	muon.matches = v.numberOfMatches();
	muon.trackerLayersWithMeasurement = hitp.trackerLayersWithMeasurement();
      }
      
      int numMuonStations = 0;
      unsigned int stationMask = static_cast<unsigned int>(v.stationMask(reco::Muon::SegmentAndTrackArbitration));
      for (int i = 0; i < 8; ++i)  // eight stations, eight bits
        if (stationMask & (1<<i)) ++numMuonStations;

      // PF Isolation
      const reco::MuonPFIsolation& pfIso03 = v.pfIsolationR03();
      muon.pfChargedIsoR03 = pfIso03.sumChargedParticlePt;
      muon.pfChargedHadIsoR03 = pfIso03.sumChargedHadronPt;
      muon.pfNeutralHadIsoR03 = pfIso03.sumNeutralHadronEt;
      muon.pfPhotonIso03 = pfIso03.sumPhotonEt;
      muon.sumPUPt03 = pfIso03.sumPUPt;

      float absiso = pfIso03.sumChargedHadronPt + std::max(0.0, pfIso03.sumNeutralHadronEt + pfIso03.sumPhotonEt - 0.5 * pfIso03.sumPUPt);
      float iso = absiso/(v.p4().pt());
      muon.pfRelIso03 = iso;

      const reco::MuonPFIsolation& pfIso04 = v.pfIsolationR04();
      muon.sumChargedParticlePt = pfIso04.sumChargedParticlePt;
      muon.sumChargedHadronPt = pfIso04.sumChargedHadronPt;
      muon.sumNeutralHadronEt = pfIso04.sumNeutralHadronEt;
      muon.sumPhotonEt =  pfIso04.sumPhotonEt;
      muon.sumPUPt = pfIso04.sumPUPt;
      absiso = pfIso04.sumChargedHadronPt + std::max(0.0, pfIso04.sumNeutralHadronEt + pfIso04.sumPhotonEt - 0.5 * pfIso04.sumPUPt);
      iso = absiso/(v.p4().pt());
      muon.pfRelIso04 = iso;

      // IP information
      muon.dB = v.dB(pat::Muon::PV2D);
      muon.edB = v.edB(pat::Muon::PV2D);

      muon.dB3D = v.dB(pat::Muon::PV3D);
      muon.edB3D = v.edB(pat::Muon::PV3D);

      // UW recommendation
      muon.isGlobalMuonPromptTight = muon::isGoodMuon(v, muon::GlobalMuonPromptTight);
      muon.isAllArbitrated         = muon::isGoodMuon(v, muon::AllArbitrated);
      muon.nChambers               = v.numberOfChambers();
      muon.nMatches                = v.numberOfMatches();
      muon.nMatchedStations        = v.numberOfMatchedStations();
      muon.stationMask             = v.stationMask();
      muon.stationGapMaskDistance  = v.stationGapMaskDistance();
      muon.stationGapMaskPull      = v.stationGapMaskPull();

      double ptError = (hasTkinRef) ? tk->ptError()/tk->pt() : 0.05;
      bool muonID = muon.isGlobalMuon && 
	            muon.isTrackerMuon && 
	            muon.isGlobalMuonPromptTight && 
	            muon.isAllArbitrated && 
	            std::fabs(dxyWrtPV) < 0.02 && 
	            std::fabs(dzWrtPV) < 0.2 && 
	            globalChi2 < 10 && 
	            ptError < 0.1 && 
	            muon.trkHits >= 10 && 
	            muon.pixHits >= 1 && 
	            numMuonStations >= 2 && 
	            muon.nMatches >= 1;
      muon.muonID = muonID;

      // Vertex information
      const reco::Candidate::Point& vertex = v.vertex();
      muon.vx = vertex.x();
      muon.vy = vertex.y();
      muon.vz = vertex.z();

      // Isolation from packed PF candidates 
      std::vector<double> isotemp;
      calcIsoFromPF(0.15, pfs, v, isotemp);
      muon.isolationMap["c15"] = isotemp;

      isotemp.clear();
      calcIsoFromPF(0.20, pfs, v, isotemp);
      muon.isolationMap["c20"] = isotemp;

      isotemp.clear();
      calcIsoFromPF(0.25, pfs, v, isotemp);
      muon.isolationMap["c25"] = isotemp;

      isotemp.clear();
      calcIsoFromPF(0.30, pfs, v, isotemp);
      muon.isolationMap["c30"] = isotemp;

      isotemp.clear();
      calcIsoFromPF(0.35, pfs, v, isotemp);
      muon.isolationMap["c35"] = isotemp;

      isotemp.clear();
      calcIsoFromPF(0.40, pfs, v, isotemp);
      muon.isolationMap["c40"] = isotemp;

      isotemp.clear();
      calcIsoFromPF(0.45, pfs, v, isotemp);
      muon.isolationMap["c45"] = isotemp;

      muon.nSegments = v.numberOfMatches(reco::Muon::SegmentArbitration);
      list_->push_back(muon);
    }
    fnMuon_ = list_->size();
    if (verbosity_) {
      int indx = 0;
      
      std::cout << "\n Printing Muon Collection Information for event \n ";
      
      std::cout << std::fixed << std::setprecision(3);
      std::cout << "  indx      pT     eta     phi  charge     dxy      dz  global tracker      PF  nMatch  Type       SIP ghostCleaned  reliso"
		<< std::endl;
      for (const auto& muon: *list_ ) {
	double iso = muon.sumChargedHadronPt + std::max(0., muon.sumNeutralHadronEt + muon.sumPhotonEt - 0.5 * muon.sumPUPt);
	std::cout << setw(6) << ++indx 
		  << setw(8) << muon.pt 
		  << setw(8) << muon.eta
		  << setw(8) << muon.phi
		  << setw(8) << muon.charge
		  << setw(8) << muon.dxyPV
		  << setw(8) << muon.dzPV
		  << setw(8) << muon.isGlobalMuon
		  << setw(8) << muon.isTrackerMuon
		  << setw(8) << muon.isPFMuon
		  << setw(8) << muon.matches
		  << setw(6) << muon.muonBestTrackType
		  << setw(10) << muon.dB3D/muon.edB3D
		  << setw(13) << muon.isghostCleaned
		  << setw(8) << iso/muon.pt
		  << std::endl;
      }
    }
  }
  else {
    edm::LogError("MuonBlock") << "Error >> Failed to get pat::Muon collection for label: "
                               << muonTag_;
  }
}

bool MuonBlock::isTrackerHighPt(const pat::Muon & mu, const reco::Vertex & primaryVertex) 
{
  //reco::TrackRef
  const auto& bestrkRef = mu.muonBestTrack();
  const auto& intrkRef = mu.innerTrack();
  if (!bestrkRef.isNonnull() || !intrkRef.isNonnull() )      return false;		
  return ( mu.numberOfMatchedStations() > 1 
           && (bestrkRef->ptError()/bestrkRef->pt()) < 0.3 
           && std::abs(bestrkRef->dxy(primaryVertex.position())) < 0.2 
           && std::abs(bestrkRef->dz(primaryVertex.position())) < 0.5 
           && intrkRef->hitPattern().numberOfValidPixelHits() > 0 
           && intrkRef->hitPattern().trackerLayersWithMeasurement() > 5 );
  
}

void MuonBlock::calcIsoFromPF(double cone, edm::Handle<pat::PackedCandidateCollection>& pfs, const pat::Muon& v, std::vector<double>& iso)
{
  // initialize sums
  double chargedhad = 0., chargedSum = 0., neutral = 0., photon = 0., pileup  = 0;
  // now get a list of the PF candidates used to build this lepton, so to exclude them
  std::vector<reco::CandidatePtr> footprint;
  std::vector< std::pair<int,double> >  chargedhadPdgPt,chargedSumPdgPt,neutralPdgPt,photonPdgPt,pileupPdgPt;
  for (unsigned int i = 0, n = v.numberOfSourceCandidatePtrs(); i < n; ++i) {
    footprint.push_back(v.sourceCandidatePtr(i));
  }
  // now loop on pf candidates
  for (unsigned int i = 0, n = pfs->size(); i < n; ++i) {
    const pat::PackedCandidate &pf = (*pfs)[i];
    if (deltaR(pf,v) < cone) {
      //pfcandidate-based footprint removal
      if (std::find(footprint.begin(), footprint.end(), reco::CandidatePtr(pfs,i)) != footprint.end()) {
        continue;
      }
      if (pf.charge() == 0) {
	if (pf.pt() > 0.5) {
          if( pf.pdgId() == 22 ) {
            photon += pf.pt();
            photonPdgPt.push_back(std::make_pair( pf.pdgId(),pf.pt() ));
          }
          else {
            neutral += pf.pt();
            neutralPdgPt.push_back(std::make_pair( pf.pdgId(),pf.pt() ));
          }
        }
      } else if (pf.vertexRef().isNonnull() &&  pf.fromPV() >= 2) {
        int pdg = std::abs(pf.pdgId());
        if( pdg!=13 && pdg!=11  ) {
	  chargedhad += pf.pt();
          chargedSum += pf.pt();
          chargedhadPdgPt.push_back(std::make_pair( pf.pdgId(),pf.pt() ));
          chargedSumPdgPt.push_back(std::make_pair( pf.pdgId(),pf.pt() ));
        } else {
	  chargedSum += pf.pt();
	  chargedSumPdgPt.push_back(std::make_pair( pf.pdgId(),pf.pt() ));
        }
      } else {
	if (pf.pt() > 0.5) {
          pileup += pf.pt();
          pileupPdgPt.push_back(std::make_pair( pf.pdgId(),pf.pt() ));
        }
      }
    }
  }
  iso.push_back(chargedhad);
  iso.push_back(chargedSum);
  iso.push_back(neutral);
  iso.push_back(photon);
  iso.push_back(pileup);

  if (verbosity_) {
    if (cone == 0.3) {
      std::cout << "Infomation of PackedPFCandidate Isolation Calculation" << std::endl;
      std::cout << "Number of Charged Hadrons within cone " << cone << chargedhadPdgPt.size() << std::endl;
      for (auto& v: chargedhadPdgPt) 
	std::cout << std::setw(4) << v.first << std::setw(6) << v.second << std::endl; 
      
      std::cout << "Number of Charged Hadrons + e/muwithin cone " << cone << chargedSumPdgPt.size() << std::endl;
      for (auto& v: chargedSumPdgPt) 
	std::cout<< std::setw(4) << v.first << std::setw(6) << v.second << std::endl; 
      
      std::cout << "Number of Neutral  Hadrons within cone " << cone << neutralPdgPt.size() << std::endl;
      for (auto& v: neutralPdgPt) 
	std::cout << std::setw(4) << v.first << std::setw(6) << v.second << std::endl; 
      
      std::cout << "Number of Photons within cone " << cone << photonPdgPt.size() << std::endl;
      for(auto& v: photonPdgPt) 
	std::cout << std::setw(4) << v.first << std::setw(6) << v.second << std::endl; 
      
      std::cout << "Number of PU contirbution within cone " << cone << pileupPdgPt.size() << std::endl;
      for(auto& v: pileupPdgPt) 
	std::cout << std::setw(4) << v.first << std::setw(6) << v.second << std::endl; 
    }
  }
}
bool MuonBlock::isSameMuon(const pat::Muon &mu1, const pat::Muon &mu2) const {
  return (& mu1 == & mu2)  ||
    (mu1.originalObjectRef() == mu2.originalObjectRef()) ||
    (mu1.reco::Muon::innerTrack().isNonnull() ?
     mu1.reco::Muon::innerTrack() == mu2.reco::Muon::innerTrack() :
     mu1.reco::Muon::outerTrack() == mu2.reco::Muon::outerTrack());
}
bool MuonBlock::isBetterMuon(const pat::Muon &mu1, const pat::Muon &mu2) const {
  if (!defaultBestMuon_) {
    MuonPointerPair pair = { &mu1, &mu2 };
    return bestMuonSelector_(pair);
  }
  if (mu2.track().isNull()) return true;
  if (mu1.track().isNull()) return false;
  if (mu1.isPFMuon() != mu2.isPFMuon()) return mu1.isPFMuon();
  if (mu1.isGlobalMuon() != mu2.isGlobalMuon()) return mu1.isGlobalMuon();
  if (mu1.charge() == mu2.charge() && deltaR2(mu1,mu2) < 0.0009) {
    return mu1.track()->ptError()/mu1.track()->pt() < mu2.track()->ptError()/mu2.track()->pt();
  } 
  else {
    int nm1 = mu1.numberOfMatches(reco::Muon::SegmentArbitration);
    int nm2 = mu2.numberOfMatches(reco::Muon::SegmentArbitration);     
    return (nm1 != nm2 ? nm1 > nm2 : mu1.pt() > mu2.pt());
  }
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
MuonBlock::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}
#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(MuonBlock);
