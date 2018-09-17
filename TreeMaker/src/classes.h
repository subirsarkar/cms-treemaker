#include "AnalysisSpace/TreeMaker/interface/PhysicsObjects.h"

#include <vector>
#include<map>
#include<string>
namespace {
  struct dictionary {
    vhtm::Event rv1;
    vhtm::GenEvent rv2;
    vhtm::Electron rv3;
    vhtm::GenParticle rv4;
    vhtm::GenJet rv5;
    vhtm::GenMET rv6;
    vhtm::MET rv7;
    vhtm::Tau rv8;
    vhtm::Muon rv9;
    vhtm::Jet rva;
    vhtm::Vertex rvb;
    vhtm::TriggerObject rvd;
    vhtm::Candidate rve;
    vhtm::Photon rvf;
    vhtm::PackedPFCandidate rvg;

    std::vector<vhtm::Electron> vrv1;
    std::vector<vhtm::GenParticle> vrv2;
    std::vector<vhtm::GenJet> vrv3;
    std::vector<vhtm::GenMET> vrv4;
    std::vector<vhtm::Tau> vrv5;
    std::vector<vhtm::Muon> vrv6;
    std::vector<vhtm::Jet> vrv7;
    std::vector<vhtm::Vertex> vrv8;
    std::vector<vhtm::TriggerObject> vrva;
    std::vector<vhtm::MET> vrvb;
    std::vector<vhtm::Event> vrvc;
    std::vector<vhtm::GenEvent> vrvd;
    std::vector<vhtm::Candidate> vrve;
    std::vector<vhtm::Photon> vrvf;
    std::vector<vhtm::PackedPFCandidate> vrvg;
    std::vector<double> vrvl;
    std::map<double, std::vector<double> > vrvm;
    std::map< std::string, std::vector<double> > vrvn;
  };
}
