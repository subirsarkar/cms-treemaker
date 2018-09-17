#include <iostream>
#include <algorithm>

#include "TTree.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/GeometryVector/interface/GlobalVector.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/EgammaReco/interface/SuperCluster.h"
#include "DataFormats/EgammaCandidates/interface/Conversion.h"
#include "DataFormats/EgammaCandidates/interface/ConversionFwd.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"

#include "AnalysisSpace/TreeMaker/plugins/PhotonBlock.h"
#include "AnalysisSpace/TreeMaker/interface/Utility.h"

// Constructor
PhotonBlock::PhotonBlock(const edm::ParameterSet& iConfig) :
  verbosity_(iConfig.getUntrackedParameter<int>("verbosity", 0)),
  photonTag_(iConfig.getParameter<edm::InputTag>("photonSrc")),
  pfcandTag_(iConfig.getParameter<edm::InputTag>("pfCands")), 
  photonToken_(consumes<pat::PhotonCollection>(photonTag_)),
  pfToken_(consumes<pat::PackedCandidateCollection>(pfcandTag_))
{}
void PhotonBlock::beginJob() 
{
  // Get TTree pointer
  std::string tree_name = "vhtree";
  TTree* tree = vhtm::Utility::getTree(tree_name);
  list_ = new std::vector<vhtm::Photon>();
  tree->Branch("Photon", "std::vector<vhtm::Photon>", &list_, 32000, -1);
  tree->Branch("nPhoton", &fnPhoton_, "fnPhoton_/I");
}
void PhotonBlock::analyze(edm::Event const& iEvent, edm::EventSetup const& iSetup) {
  // Reset the TClonesArray and the nObj variables
  list_->clear();
  fnPhoton_ = 0;

  edm::Handle<pat::PhotonCollection> photons;
  bool found = iEvent.getByToken(photonToken_, photons);

  edm::Handle<pat::PackedCandidateCollection> pfs;
  iEvent.getByToken(pfToken_, pfs);
  
  if (found && photons.isValid()) {
    edm::LogInfo("PhotonBlock") << "Total # PAT Photons: " << photons->size();
    for (auto const& v: *photons) {
      if (list_->size() == kMaxPhoton) {
	edm::LogInfo("PhotonBlock") << "Too many PAT Photon, fnPhoton = " 
                                    << fnPhoton_; 
	break;
      }
      vhtm::Photon photon;
      photon.et     = v.et();
      photon.eta    = v.eta();
      photon.clusterEta = v.caloPosition().eta();
      photon.phi    = v.phi();
      photon.clusterPhi = v.caloPosition().phi();
      photon.energy = v.energy();
      photon.theta  = v.theta();
      photon.vx     = v.vx();
      photon.vy     = v.vy();
      photon.vz     = v.vz();

      const reco::SuperClusterRef sCluster = v.superCluster(); 
      //if( !sCluster ) {
      photon.scEnergy    = sCluster->energy();
      photon.scEta       = sCluster->eta();
      photon.scPhi       = sCluster->phi();
      photon.scSize      = sCluster->clustersSize();
      photon.scEtaWidth  = sCluster->etaWidth();
      photon.scPhiWidth  = sCluster->phiWidth();
      photon.scEt        = sCluster->energy()/cosh(sCluster->eta());
      photon.scRawEnergy = sCluster->rawEnergy();
      photon.scx         = sCluster->x();
      photon.scy         = sCluster->y();
      photon.scz         = sCluster->z();
     // }

      photon.isoEcalRecHit03    = v.ecalRecHitSumEtConeDR03();
      photon.isoHcalRecHit03    = v.hcalTowerSumEtConeDR03();
      photon.isoSolidTrkCone03  = v.trkSumPtSolidConeDR03();
      photon.isoHollowTrkCone03 = v.trkSumPtHollowConeDR03();
      photon.nTrkSolidCone03    = v.nTrkSolidConeDR03();
      photon.nTrkHollowCone03   = v.nTrkHollowConeDR03();

      photon.isoEcalRecHit04    = v.ecalRecHitSumEtConeDR04();
      photon.isoHcalRecHit04    = v.hcalTowerSumEtConeDR04();
      photon.isoSolidTrkCone04  = v.trkSumPtSolidConeDR04();
      photon.isoHollowTrkCone04 = v.trkSumPtHollowConeDR04();
      photon.nTrkSolidCone04    = v.nTrkSolidConeDR04();
      photon.nTrkHollowCone04   = v.nTrkHollowConeDR04();

      photon.hasPixelSeed       = v.hasPixelSeed(); 
      photon.ecalIso            = v.ecalIso();
      photon.hcalIso            = v.hcalIso();
      photon.trackIso           = v.trackIso();

      photon.chargedHadIso      = v.chargedHadronIso();
      photon.neutralHadIso      = v.neutralHadronIso();
      photon.photonIso          = v.photonIso();

      int fidFlag = 0;
      if (v.isEB())        fidFlag |= (1 << 0);
      if (v.isEE())        fidFlag |= (1 << 1);
      if (v.isEBEtaGap())  fidFlag |= (1 << 2);
      if (v.isEBPhiGap())  fidFlag |= (1 << 3);
      if (v.isEERingGap()) fidFlag |= (1 << 4);
      if (v.isEEDeeGap())  fidFlag |= (1 << 5);
      if (v.isEBEEGap())   fidFlag |= (1 << 6);
      photon.fidFlag = fidFlag;

      photon.isEB               = v.isEB() ? true : false;
      photon.isEE               = v.isEE() ? true : false;
      photon.isEBGap            = v.isEBGap() ? true : false ;
      photon.isEEGap            = v.isEEGap() ? true : false;
      photon.isEBEEGap          = v.isEBEEGap() ? true : false;

      photon.r9                 = v.r9();
      photon.hoe                = v.hadronicOverEm();
      photon.sigmaEtaEta        = v.sigmaEtaEta();
      photon.sigmaIEtaIEta      = v.sigmaIetaIeta();
      photon.e1x5               = v.e1x5();
      photon.e2x5               = v.e2x5();
      photon.e3x3               = v.e3x3();
      photon.e5x5               = v.e5x5();
      photon.r1x5               = v.r1x5();
      photon.r2x5               = v.r2x5();
      photon.maxEnergyXtal      = v.maxEnergyXtal();

      photon.hasConversionTracks = v.hasConversionTracks();      
      photon.passElectronVeto = v.passElectronVeto();
/*
      if (v.hasConversionTracks()) {
        const reco::ConversionRefVector conversions = v.conversions();
        if( !conversions.empty()  ) { 
        for (edm::RefVector<reco::ConversionCollection>::const_iterator jt  = conversions.begin();
                                                                        jt != conversions.end(); 
                                                                      ++jt) 
	{
          const reco::Conversion& obj = (**jt);
	  if (obj.nTracks() < 2 or
              !obj.conversionVertex().isValid()) continue;
          photon.nTracks = obj.nTracks();
          photon.isConverted = obj.isConverted();
          photon.pairInvMass = obj.pairInvariantMass();
          photon.pairCotThetaSeparation
     	        = obj.pairCotThetaSeparation();

	  math::XYZVectorF  mom = obj.pairMomentum();
          photon.pairPx = mom.x();
          photon.pairPy = mom.y();
          photon.pairPz = mom.z();

          const reco::Vertex &vertex = obj.conversionVertex();
          photon.conv_vx = vertex.x();
          photon.conv_vy = vertex.y();
          photon.conv_vz = vertex.z();

	  photon.eovp              = obj.EoverP();
	  photon.zpv               = obj.zOfPrimaryVertexFromTracks();
	  photon.distOfMinApproach = obj.distOfMinimumApproach();
	  photon.dPhiTracksAtVtx   = obj.dPhiTracksAtVtx();
	  photon.dPhiTracksAtEcal  = obj.dPhiTracksAtEcal();
	  photon.dEtaTracksAtEcal  = obj.dEtaTracksAtEcal();
        }    
        }
      }
*/
      std::vector<double> isotemp;
      calcIsoFromPF(0.15, pfs, v, isotemp);
      photon.isolationMap["c15"] = isotemp;
      isotemp.clear();
      calcIsoFromPF(0.20, pfs, v, isotemp);
      photon.isolationMap["c20"] = isotemp;
      isotemp.clear();
      calcIsoFromPF(0.25, pfs, v, isotemp);
      photon.isolationMap["c25"] = isotemp;
      isotemp.clear();
      calcIsoFromPF(0.30, pfs, v, isotemp);
      photon.isolationMap["c30"] = isotemp;
      isotemp.clear();
      calcIsoFromPF(0.35, pfs, v, isotemp);
      photon.isolationMap["c35"] = isotemp;
      isotemp.clear();
      calcIsoFromPF(0.40, pfs, v, isotemp);
      photon.isolationMap["c40"] = isotemp;
      isotemp.clear();
      calcIsoFromPF(0.45, pfs, v, isotemp);
      photon.isolationMap["c45"] = isotemp;

      list_->push_back(photon);
    }
    fnPhoton_ = list_->size();
  }
  else {
    edm::LogError("PhotonBlock") << "Error >> Failed to get pat::Photon for label: " 
                                 << photonTag_;
  }
}
void PhotonBlock::calcIsoFromPF(double cone, edm::Handle<pat::PackedCandidateCollection>& pfs, const pat::Photon& v, std::vector<double>& iso)
{
  // initialize sums
  double charged = 0, 
    neutral = 0, 
    pileup = 0;
  // now get a list of the PF candidates used to build this lepton, so to exclude them
  std::vector<reco::CandidatePtr> footprint;
  for (unsigned int i = 0, n = v.numberOfSourceCandidatePtrs(); i < n; ++i)
    footprint.push_back(v.sourceCandidatePtr(i));

  // now loop on pf candidates
  for (unsigned int i = 0, n = pfs->size(); i < n; ++i) {
    const pat::PackedCandidate &pf = (*pfs)[i];
    if (deltaR(pf,v) < cone) {
      // pfcandidate-based footprint removal
      if (std::find(footprint.begin(), footprint.end(), reco::CandidatePtr(pfs,i)) != footprint.end()) 
	continue;
    
      if (pf.charge() == 0) {
	if (pf.pt() > 0.5) neutral += pf.pt();
      } 
      else if (pf.fromPV() >= 2) {
	if (pf.pt() > 0.2) charged += pf.pt();
      } 
      else {
	if (pf.pt() > 0.5) pileup += pf.pt();  
      }
    }
  }
  iso.push_back(charged); 
  iso.push_back(neutral);
  iso.push_back(pileup);
}
#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(PhotonBlock);
