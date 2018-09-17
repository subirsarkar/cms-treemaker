#include <iostream>

#include "TTree.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/Ref.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/PatCandidates/interface/TauPFEssential.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "FWCore/ParameterSet/interface/FileInPath.h"
#include "AnalysisSpace/TreeMaker/interface/PhysicsObjects.h"
#include "AnalysisSpace/TreeMaker/plugins/TauBlock.h"
#include "AnalysisSpace/TreeMaker/interface/Utility.h"

namespace tb
{
  template<typename T>
  bool isValidRef(const edm::Ref<T>& ref) {
    return ( (ref.isAvailable() || ref.isTransient()) && ref.isNonnull() );
  }
}
TauBlock::TauBlock(const edm::ParameterSet& iConfig) :
  verbosity_(iConfig.getUntrackedParameter<int>("verbosity", 0)),
  bsCorr_(iConfig.getUntrackedParameter<bool>("beamSpotCorr", true)),
  tauTag_(iConfig.getParameter<edm::InputTag>("patTauSrc")),
  vertexTag_(iConfig.getParameter<edm::InputTag>("vertexSrc")),
  bsTag_(iConfig.getParameter<edm::InputTag>("offlineBeamSpot")),
  tauToken_(consumes<pat::TauCollection>(tauTag_)),
  vertexToken_(consumes<reco::VertexCollection>(vertexTag_)),
  bsToken_(consumes<reco::BeamSpot>(bsTag_))
{
}
TauBlock::~TauBlock() {}
void TauBlock::beginJob()
{
  // Get TTree pointer
  TTree* tree = vhtm::Utility::getTree("vhtree");
  list_ = new std::vector<vhtm::Tau>();
  tree->Branch("Tau", "std::vector<vhtm::Tau>", &list_, 32000, -1);
  tree->Branch("nTau", &fnTau_, "fnTau_/I");
}
void TauBlock::analyze(edm::Event const& iEvent, edm::EventSetup const& iSetup) {
  // Reset the vector and the nObj variables
  list_->clear();
  fnTau_ = 0;

  edm::Handle<pat::TauCollection> taus;
  bool found = iEvent.getByToken(tauToken_, taus);
  
  if (found && taus.isValid()) {
    edm::Handle<reco::VertexCollection> primaryVertices;
    iEvent.getByToken(vertexToken_, primaryVertices);

    edm::Handle<reco::BeamSpot> beamSpot;
    if (bsCorr_)
      iEvent.getByToken(bsToken_, beamSpot);

    edm::LogInfo("TauBlock") << "Total # PAT Taus: " << taus->size();
    for (const auto& v: *taus) {
      if (list_->size() == kMaxTau_) {
        edm::LogInfo("TauBlock") << "Too many PAT Taus, fnTau = " << list_->size();
        break;
      }
      vhtm::Tau tau;

      // Store Tau variables
      tau.eta    = v.eta();
      tau.phi    = v.phi();
      tau.pt     = v.pt();
      tau.energy = v.energy();
      tau.charge = v.charge();
      if (v.leadChargedHadrCand().isNonnull()) {
        // We know that it returns a PackedCandidate
        //const pat::PackedCandidate& trk = dynamic_cast<const pat::PackedCandidate&>(*v.leadChargedHadrCand());
        const auto& trk = dynamic_cast<const pat::PackedCandidate&>(*v.leadChargedHadrCand());

        if (primaryVertices.isValid()) {
          edm::LogInfo("TauBlock") << "Total # Primary Vertices: " 
                                   << primaryVertices->size();

          // IP of leadChargedHadrCand wrt event PV
          const auto& vit = primaryVertices->front();
          tau.dxyPV = trk.dxy(vit.position());
          tau.dzPV  = trk.dz(vit.position());
          tau.dxyPVError = trk.dxyError();
          tau.dzPVError = trk.dzError();
       
          // IP of leadChargedHadrCand wrt closest PV
          // Vertex association
          double minVtxDist3D = 9999.;
          int indexVtx = -1;
          double vertexDz = 9999.;
          double vertexDxy = 9999.;
          for (auto vit = primaryVertices->begin(); vit != primaryVertices->end(); ++vit) {
            double dxy = trk.dxy(vit->position());
            double dz  = trk.dz(vit->position());
            double dist3D = std::sqrt(pow(dxy,2) + pow(dz,2));
            if (dist3D < minVtxDist3D) {
              minVtxDist3D = dist3D;
              indexVtx = int(std::distance(primaryVertices->begin(), vit));
              vertexDxy = dxy;
              vertexDz = dz;
            }
          }
          tau.vtxIndex = indexVtx;
          tau.vtxDxy   = vertexDxy;
          tau.vtxDz    = vertexDz;
        }
        else {
          edm::LogError("TauBlock") << "Error >> Failed to get VertexCollection for label: "
                                    << vertexTag_;
        }
      }
      // Leading particle pT
      tau.leadChargedParticlePt = 
           v.leadChargedHadrCand().isNonnull() ? v.leadChargedHadrCand()->pt(): 0.;
      tau.leadNeutralParticlePt = 
           v.leadNeutralCand().isNonnull() ? v.leadNeutralCand()->et(): 0.;
      tau.leadParticlePt = 
           v.leadCand().isNonnull() ? v.leadCand()->et(): 0.;

      // Signal Constituents
      // Charged hadrons
      for (const auto& iCand: v.signalChargedHadrCands()) {
        const reco::Candidate& cand = (*iCand);
	vhtm::Candidate c(cand.pt(), cand.eta(), cand.phi());
        tau.sigChHadList.push_back(c);
      }
      // Neutral hadrons
      for (const auto& iCand: v.signalNeutrHadrCands()) {
        const reco::Candidate& cand = (*iCand);
	vhtm::Candidate c(cand.pt(), cand.eta(), cand.phi());
        tau.sigNeHadList.push_back(c);
      }
      // Photons
      double sigPtSumGamma = 0;
      for (const auto& iCand: v.signalGammaCands()) {
        const reco::Candidate& cand = (*iCand);
	vhtm::Candidate c(cand.pt(), cand.eta(), cand.phi());
        tau.sigGammaList.push_back(c);
        sigPtSumGamma += cand.pt();
      }
      // Isolation Constituents
      // Charged hadrons
      double isoPtSumChHad = 0;
      for (const auto& iCand: v.isolationChargedHadrCands()) {
        const reco::Candidate& cand = (*iCand);
	vhtm::Candidate c(cand.pt(), cand.eta(), cand.phi());
        tau.isoChHadList.push_back(c);
        isoPtSumChHad += cand.pt();
      }
      // Neutral hadrons
      double isoPtSumNeHad = 0;
      for (const auto& iCand: v.isolationNeutrHadrCands()) {
        const reco::Candidate& cand = (*iCand);
	vhtm::Candidate c(cand.pt(), cand.eta(), cand.phi());
        tau.isoNeHadList.push_back(c);
        isoPtSumNeHad += cand.pt();
      }
      // Photons
      double isoPtSumGamma = 0;
      for (const auto& iCand: v.isolationGammaCands()) {
        const reco::Candidate& cand = (*iCand);
	vhtm::Candidate c(cand.pt(), cand.eta(), cand.phi());
        tau.isoGammaList.push_back(c);
        isoPtSumGamma += cand.pt();
      }
      // PtSum
      tau.ptSumChargedHadronsIsoCone = isoPtSumChHad;
      tau.ptSumNeutralHadronsIsoCone = isoPtSumNeHad;
      tau.ptSumPhotonsIsoCone        = isoPtSumGamma;

      // tau id. discriminators
      tau.decayModeFinding = v.tauID("decayModeFinding");
      tau.decayModeFindingNewDMs = v.tauID("decayModeFindingNewDMs");
      tau.decayMode = v.decayMode();

      // discriminators against muons
      tau.againstMuonLoose3 = v.tauID("againstMuonLoose3");
      tau.againstMuonTight3 = v.tauID("againstMuonTight3");

      // discriminators against electrons
      tau.againstElectronVLooseMVA = v.tauID("againstElectronVLooseMVA6");
      tau.againstElectronLooseMVA  = v.tauID("againstElectronLooseMVA6");
      tau.againstElectronMediumMVA = v.tauID("againstElectronMediumMVA6");
      tau.againstElectronTightMVA  = v.tauID("againstElectronTightMVA6");
      tau.againstElectronVTightMVA = v.tauID("againstElectronVTightMVA6");

      // DB Corrected Isolation
      tau.byLooseCombinedIsolationDeltaBetaCorr3Hits  = v.tauID("byLooseCombinedIsolationDeltaBetaCorr3Hits");
      tau.byMediumCombinedIsolationDeltaBetaCorr3Hits = v.tauID("byMediumCombinedIsolationDeltaBetaCorr3Hits");
      tau.byTightCombinedIsolationDeltaBetaCorr3Hits  = v.tauID("byTightCombinedIsolationDeltaBetaCorr3Hits");
      tau.byCombinedIsolationDeltaBetaCorrRaw3Hits    = v.tauID("byCombinedIsolationDeltaBetaCorrRaw3Hits");

      tau.byVLooseIsolationMVArun2v1DBoldDMwLT = v.tauID("byVLooseIsolationMVArun2v1DBoldDMwLT");
      tau.byLooseIsolationMVArun2v1DBoldDMwLT  = v.tauID("byLooseIsolationMVArun2v1DBoldDMwLT");
      tau.byMediumIsolationMVArun2v1DBoldDMwLT = v.tauID("byMediumIsolationMVArun2v1DBoldDMwLT");
      tau.byTightIsolationMVArun2v1DBoldDMwLT  = v.tauID("byTightIsolationMVArun2v1DBoldDMwLT");
      tau.byVTightIsolationMVArun2v1DBoldDMwLT = v.tauID("byVTightIsolationMVArun2v1DBoldDMwLT");

      // Isolation variables suggested by Christian
      tau.chargedIsoPtSum = v.tauID("chargedIsoPtSum");
      tau.neutralIsoPtSum = v.tauID("neutralIsoPtSum");
      tau.puCorrPtSum     = v.tauID("puCorrPtSum");
#if 0
      // who wants to store all the tauIDs?
      std::cout << ">>> tauID(label) = value" << std::endl;
      for (const auto& pa: v.tauIDs())
	std::cout << pa.first << "=" << pa.second << std::endl;
#endif
      // kinematic variables for PFJet associated to PFTau
      tau.jetPt  = v.pfEssential().p4Jet_.Pt();
      tau.jetEta = v.pfEssential().p4Jet_.Eta();
      tau.jetPhi = v.pfEssential().p4Jet_.Phi();

      // The following has to be computed from the PackedCandidates of Tau within the signal cone(?)
      tau.emFraction = (sigPtSumGamma > 0) ? sigPtSumGamma/v.pt() : 0;

      // Vertex information
      const auto& vertex = v.vertex();
      tau.vx = vertex.x();
      tau.vy = vertex.y();
      tau.vz = vertex.z();

      tau.zvertex = v.vz(); // distance from the primary vertex
      tau.mass    = v.p4().M();

      tau.dxySig = v.dxy_Sig();

      // add particle to the list
      list_->push_back(tau);
    }
    fnTau_ = list_->size();
  }
  else {
    edm::LogError("TauBlock") << "Error! Failed to get pat::Tau collection for label: "
                              << tauTag_;
  }
}
#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(TauBlock);
