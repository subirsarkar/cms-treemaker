#include <iostream>
#include <algorithm>
#include "TFile.h"
#include "TTree.h"
#include "TROOT.h"
#include "TVector3.h"
#include "TLorentzVector.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "DataFormats/TrackReco/interface/HitPattern.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"

#include "AnalysisSpace/TreeMaker/interface/PhysicsObjects.h"
#include "AnalysisSpace/TreeMaker/plugins/ElectronBlock.h"
#include "AnalysisSpace/TreeMaker/interface/Utility.h"

// Constructor
ElectronBlock::ElectronBlock(const edm::ParameterSet& iConfig):
  verbosity_(iConfig.getUntrackedParameter<int>("verbosity", 0)),
  bsCorr_(iConfig.getUntrackedParameter<bool>("beamSpotCorr", false)),
  trigMode_(iConfig.getUntrackedParameter<bool>("useTrigMode", false)),
  bsTag_(iConfig.getParameter<edm::InputTag>("offlineBeamSpot")),
  vertexTag_(iConfig.getParameter<edm::InputTag>("vertexSrc")),
  electronTag_(iConfig.getParameter<edm::InputTag>("electronSrc")),
  pfcandTag_(iConfig.getParameter<edm::InputTag>("pfCands")),
  bsToken_(consumes<reco::BeamSpot>(bsTag_)),
  vertexToken_(consumes<reco::VertexCollection>(vertexTag_)),
  electronToken_(consumes<pat::ElectronCollection>(electronTag_)),
  pfToken_(consumes<pat::PackedCandidateCollection>(pfcandTag_))
{
}
ElectronBlock::~ElectronBlock() {
  delete list_;
}
void ElectronBlock::beginJob()
{
  // Get TTree pointer
  std::string tree_name = "vhtree";
  TTree* tree = vhtm::Utility::getTree(tree_name);
  list_ = new std::vector<vhtm::Electron>();
  tree->Branch("Electron", "std::vector<vhtm::Electron>", &list_, 32000, -1);
  tree->Branch("nElectron", &fnElectron_, "fnElectron_/I");
}
void ElectronBlock::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  // Reset the vector and the nObj variables
  list_->clear();
  fnElectron_ = 0;

  edm::Handle<pat::ElectronCollection> electrons;
  bool found = iEvent.getByToken(electronToken_, electrons);

  edm::Handle<pat::PackedCandidateCollection> pfs;
  iEvent.getByToken(pfToken_, pfs);

  if (found && electrons.isValid()) {
    edm::Handle<reco::BeamSpot> beamSpot;
    if (bsCorr_) iEvent.getByToken(bsToken_, beamSpot);

    edm::Handle<reco::VertexCollection> primaryVertices;
    iEvent.getByToken(vertexToken_, primaryVertices);

    edm::LogInfo("ElectronBlock") << "Total # PAT Electrons: " << electrons->size();
    
    for (const auto& v: *electrons) {
      if (list_->size() == kMaxElectron_) {
        edm::LogInfo("ElectronBlock") << "Too many PAT Electrons, fnElectron = "
                                      << list_->size();
        break;
      }
      bool hasGsfTrack = v.gsfTrack().isNonnull() ? true : false;

      vhtm::Electron electron;
      electron.ecalDriven      = v.ecalDrivenSeed();
      electron.eta             = v.eta();
      electron.phi             = v.phi();
      electron.pt              = v.pt();
      electron.hasGsfTrack     = hasGsfTrack;
      electron.energy          = v.energy();
      electron.caloEnergy      = v.ecalEnergy();
      electron.charge          = v.charge();
  
      float nMissingHits = 0;
      double dxyWrtPV = -99.;
      double dzWrtPV = -99.;
      
      electron.isGap = v.isGap();

      electron.BDT = v.userFloat("electronBDT"); 
      electron.BDTpreComp = v.userFloat("ElectronMVAEstimatorRun2Spring15NonTrig25nsV1Values");
      if (hasGsfTrack) {
        reco::GsfTrackRef tk = v.gsfTrack();
        electron.trackPt = tk->pt();

        // Hit pattern
        const reco::HitPattern& hitp = tk->hitPattern();
        electron.pixHits = hitp.numberOfValidPixelHits();
        electron.trkHits = hitp.numberOfValidTrackerHits();

        electron.nValidHits  = tk->numberOfValidHits();
        //electron.missingHits = hitp.numberOfHits(reco::HitPattern::MISSING_INNER_HITS);
        electron.missingHits = hitp.numberOfAllHits(reco::HitPattern::MISSING_INNER_HITS);

        double trkd0 = tk->d0();
        double trkdz = tk->dz();
        if (bsCorr_) {
          if (beamSpot.isValid()) {
            trkd0 = -(tk->dxy(beamSpot->position()));
            trkdz = tk->dz(beamSpot->position());
          }
          else
            edm::LogError("ElectronBlock") << "Error >> Failed to get BeamSpot for label: "
                                           << bsTag_;
        }
        electron.trkD0 = trkd0;
        electron.trkDz = trkdz;

        if (primaryVertices.isValid()) {
          const reco::Vertex& vit = primaryVertices->front(); // Highest sumPt vertex
          dxyWrtPV = tk->dxy(vit.position());
          dzWrtPV  = tk->dz(vit.position());
          electron.dxyPV = dxyWrtPV;
          electron.dzPV  = dzWrtPV;

          // Vertex association
          double minVtxDist3D = 9999.;
          int indexVtx = -1;
          double vertexDistZ = 9999.;
          edm::LogInfo("ElectronBlock") << "Total # Primary Vertices: " << primaryVertices->size();
          for (auto vit = primaryVertices->begin(); vit != primaryVertices->end(); ++vit) {
            double dxy = tk->dxy(vit->position());
            double dz  = tk->dz(vit->position());
            double dist3D = std::sqrt(pow(dxy, 2) + pow(dz, 2));
            if (dist3D < minVtxDist3D) {
              minVtxDist3D = dist3D;
              indexVtx = int(std::distance(primaryVertices->begin(), vit));
              vertexDistZ = dz;
            }
          }
          electron.vtxDist3D = minVtxDist3D;
          electron.vtxIndex = indexVtx;
          electron.vtxDistZ = vertexDistZ;
        }
        else {
          edm::LogError("ElectronBlock") << "Error >> Failed to get VertexCollection for label: "
                                         << vertexTag_;
        }
      }
      // ID variables
      float dPhi  = v.deltaPhiSuperClusterTrackAtVtx();
      float dEta  = v.deltaEtaSuperClusterTrackAtVtx();
      float sihih = v.sigmaIetaIeta();
      float hoe   = v.hadronicOverEm(); // v.hcalOverEcal();

      electron.hoe           = hoe;
      electron.eop           = v.eSuperClusterOverP();
      electron.sigmaEtaEta   = v.sigmaEtaEta();
      electron.sigmaIEtaIEta = sihih;
      electron.deltaPhiTrkSC = dPhi;
      electron.deltaEtaTrkSC = dEta;
      electron.classif       = v.classification();

      // SC associated with electron
      electron.scEn  = v.superCluster()->energy();
      electron.scEta = v.superCluster()->eta();
      electron.scPhi = v.superCluster()->phi();
      electron.scET  = v.superCluster()->energy()/cosh(v.superCluster()->eta());
      electron.scRawEnergy = v.superCluster()->rawEnergy();
      electron.relIso = (v.trackIso() + v.ecalIso() + v.hcalIso())/v.pt();

      // PF Isolation
      reco::GsfElectron::PflowIsolationVariables pfIso = v.pfIsolationVariables();
      electron.sumChargedHadronPt = pfIso.sumChargedHadronPt;
      electron.sumPUPt = pfIso.sumPUPt;
      electron.sumNeutralHadronEt = pfIso.sumNeutralHadronEt ;
      electron.sumPhotonEt = pfIso.sumPhotonEt;
      float absiso = pfIso.sumChargedHadronPt + std::max(0.0, pfIso.sumNeutralHadronEt + pfIso.sumPhotonEt - 0.5 * pfIso.sumPUPt);
      float iso = absiso/(v.p4().pt());
      electron.pfRelIso = iso;

      // isolation information
      electron.chargedHadronIso = v.chargedHadronIso();
      electron.neutralHadronIso = v.neutralHadronIso();
      electron.photonIso        = v.photonIso();
  
      // IP information
      electron.dB  = v.dB(pat::Electron::PV2D);
      electron.edB = v.edB(pat::Electron::PV2D);

      electron.dB3D = v.dB(pat::Electron::PV3D);
      electron.edB3D = v.edB(pat::Electron::PV3D);

      // Bremstrahlung information
      electron.nBrems = v.numberOfBrems();
      electron.fbrem  = v.fbrem();

      // other useful quantities
      bool myTrigPresel = false;
      double pt = v.pt();
      if (std::abs(v.superCluster()->eta()) < 1.485) {
	if (sihih < 0.014 &&
	    hoe < 0.15 &&
	    v.dr03TkSumPt()/pt < 0.2 &&
	    v.dr03EcalRecHitSumEt()/pt < 0.2 &&
	    v.dr03HcalTowerSumEt()/pt < 0.2 &&
	    nMissingHits == 0)
	  myTrigPresel = true;
      }
      else {
	if (sihih < 0.035 &&
	    hoe < 0.10 &&
	    v.dr03TkSumPt()/pt < 0.2 &&
	    v.dr03EcalRecHitSumEt()/pt < 0.2 &&
	    v.dr03HcalTowerSumEt()/pt < 0.2 &&
	    nMissingHits == 0)
	  myTrigPresel = true;
      }
      electron.mvaPreselection = myTrigPresel;
      
      bool flagEB = (v.isEB() ? true : false);
      bool flagEE = (v.isEE() ? true : false); 
      bool passconv = v.passConversionVeto();
      electron.hasMatchedConv = (passconv ? false : true);
  
      bool mvaPreselection = passconv && nMissingHits <= 0 && dxyWrtPV < 0.02 && dzWrtPV < 0.1 &&
	((flagEB 
          && sihih < 0.01
	  && std::fabs(dEta) < 0.007
	  && std::fabs(dPhi) < 0.15
          && hoe < 0.12
	  && v.dr03TkSumPt()/pt < 0.20
	  && (std::max(v.dr03EcalRecHitSumEt() - 1.0, 0.0))/pt < 0.20
	  && v.dr03HcalTowerSumEt()/pt < 0.20
	  ) ||
	 (flagEE 
          && sihih < 0.03
	  && std::fabs(dEta) < 0.009
	  && std::fabs(dPhi) < 0.10
          && hoe < 0.10
	  && v.dr03TkSumPt()/pt < 0.20
	  && (std::max(v.dr03EcalRecHitSumEt() - 1.0, 0.0))/pt < 0.20
	  && v.dr03HcalTowerSumEt()/pt < 0.20
	  ));
      electron.isTriggerElectron = mvaPreselection;

      // Fiducial flag
      int fidFlag = 0;
      if (v.isEB())        fidFlag |= (1 << 0);
      if (v.isEE())        fidFlag |= (1 << 1);
      if (v.isEBEtaGap())  fidFlag |= (1 << 2);
      if (v.isEBPhiGap())  fidFlag |= (1 << 3);
      if (v.isEERingGap()) fidFlag |= (1 << 4);
      if (v.isEEDeeGap())  fidFlag |= (1 << 5);
      if (v.isEBEEGap())   fidFlag |= (1 << 6);
      electron.fidFlag = fidFlag;

      // Vertex information
      const reco::Candidate::Point& vertex = v.vertex();
      electron.vx = vertex.x();
      electron.vy = vertex.y();
      electron.vz = vertex.z();

      for (const auto& pa: v.electronIDs())
	electron.idmap[pa.first] = pa.second;

#if 0
      electron.isPassMedium_MVAIDWP90 = v.userInt("isPassMedium_MVAIDWP90"); 
      electron.isPassTight_MVAIDWP80 = v.userInt("isPassTight_MVAIDWP80"); 
      electron.mvaId = v.userFloat("mvaId");
#endif
      // Isolation from packed PF candidates 
      std::vector<double> isotemp;

      electron.ElectronMVAEstimatorRun2Spring15NonTrig25nsV1Values = v.userFloat("ElectronMVAEstimatorRun2Spring15NonTrig25nsV1Values");
      calcIsoFromPF(0.15, pfs, v, isotemp);
      electron.isolationMap["c15"] = isotemp;

      isotemp.clear();
      calcIsoFromPF(0.20, pfs, v, isotemp);
      electron.isolationMap["c20"] = isotemp;

      isotemp.clear();
      calcIsoFromPF(0.25, pfs, v, isotemp);
      electron.isolationMap["c25"] = isotemp;

      isotemp.clear();
      calcIsoFromPF(0.30, pfs, v, isotemp);
      electron.isolationMap["c30"] = isotemp;

      isotemp.clear();
      calcIsoFromPF(0.35, pfs, v, isotemp);
      electron.isolationMap["c35"] = isotemp;

      isotemp.clear();
      calcIsoFromPF(0.40, pfs, v, isotemp);
      electron.isolationMap["c40"] = isotemp;

      isotemp.clear();
      calcIsoFromPF(0.45, pfs, v, isotemp);
      electron.isolationMap["c45"] = isotemp;

      // associatedPackedPFCandidatesP4 part
      const auto& pfref = v.associatedPackedPFCandidates();
      //for (edm::RefVector<pat::PackedCandidateCollection>::const_iterator it = pfref.begin();  it != pfref.end(); ++it) {
      for (const auto& el: pfref) {
        TLorentzVector p4ref;
        p4ref.SetPtEtaPhiE(el->pt(), el->eta(), el->phi(), el->energy());
        electron.associatedPackedPFCandidatesP4.push_back(p4ref);
      }

      list_->push_back(electron);
    }
    fnElectron_ = list_->size();
  }
  else {
    edm::LogError("ElectronBlock") << "Error >> Failed to get pat::Electron Collection for label: "
                                   << electronTag_;
  }
}
void ElectronBlock::calcIsoFromPF(double cone, edm::Handle<pat::PackedCandidateCollection>& pfs, const pat::Electron& v, std::vector<double>& iso)
{
  // initialize sums
  double chargedHadSum = 0., 
    chargedSum = 0., 
    neutralSum = 0., 
    photonSum = 0., 
    pileupSum  = 0;

  // now get a list of the PF candidates used to build this lepton, so to exclude them
  std::vector<reco::CandidatePtr> footprint;
  for (unsigned int i = 0, n = v.numberOfSourceCandidatePtrs(); i < n; ++i)
    footprint.push_back(v.sourceCandidatePtr(i));

  // now loop on pf candidates
  for (unsigned int i = 0, n = pfs->size(); i < n; ++i) {
    const pat::PackedCandidate& pf = (*pfs)[i];
    if (deltaR(pf,v) < cone) {

      // pfcandidate-based footprint removal
      if (std::find(footprint.begin(), footprint.end(), reco::CandidatePtr(pfs, i)) != footprint.end())
        continue;

      if (pf.charge() == 0) {
        if (pf.pt() > 0.5) {
          if (pf.vertexRef().isNonnull() && pf.pdgId() == 22)
	    photonSum += pf.pt();
          else 
            neutralSum += pf.pt();
        }
      } 
      else if (pf.fromPV() >= 2) {
        int pdg = std::abs(pf.pdgId());
        if (pdg != 13 && pdg != 11) chargedHadSum += pf.pt();
	chargedSum += pf.pt();
      } 
      else {
        if (pf.pt() > 0.5) pileupSum += pf.pt();
      }
    }
  }
  iso.push_back(chargedHadSum);
  iso.push_back(chargedSum);
  iso.push_back(neutralSum);
  iso.push_back(photonSum);
  iso.push_back(pileupSum);
}
// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
ElectronBlock::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

DEFINE_FWK_MODULE(ElectronBlock);
