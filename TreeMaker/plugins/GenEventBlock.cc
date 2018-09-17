#include <algorithm>
#include <iostream>
#include <memory>
#include <iterator>
#include <cctype>

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include <boost/foreach.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string.hpp>

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/LHERunInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"
#include "FWCore/Framework/interface/Run.h"

#include "TTree.h"
#include "AnalysisSpace/TreeMaker/interface/PhysicsObjects.h"
#include "AnalysisSpace/TreeMaker/plugins/GenEventBlock.h"
#include "AnalysisSpace/TreeMaker/interface/Utility.h"

using std::cout;
using std::cerr;
using std::endl;
using std::vector;
using std::string;
using std::stoi;
using std::stringstream;

GenEventBlock::GenEventBlock(const edm::ParameterSet& iConfig) :
         verbosity_(iConfig.getUntrackedParameter<int>("verbosity", 0)),
   storePDFWeights_(iConfig.getUntrackedParameter<bool>("StorePDFWeights", false)),
               tag_(iConfig.getUntrackedParameter<string>("tag", "initrwgt")),
  isStandardSample_(iConfig.getUntrackedParameter<bool>("isStandardSample", true)),
   doAlphasWeights_(iConfig.getUntrackedParameter<bool>("doAlphasWeights", true)),
    doScaleWeights_(iConfig.getUntrackedParameter<bool>("doScaleWeights", true)),
    nPdfEigWeights_(iConfig.getParameter<unsigned int>("nPdfEigWeights")),
            pdfset_(iConfig.getUntrackedParameter<string>("pdfset", "NNPDF30_lo_as_0130.LHgrid")),
     mc2hessianCSV_(iConfig.getParameter<edm::FileInPath>("mc2hessianCSV")),
       genEventTag_(iConfig.getParameter<edm::InputTag>("GenEventInfoTag")),
       lheEventTag_(iConfig.getParameter<edm::InputTag>("LHEEventTag")),
     pdfWeightsTag_(iConfig.getParameter<edm::InputTag>("PDFWeightsTag")),
     genEventToken_(consumes<GenEventInfoProduct>(genEventTag_)),
     lheEventToken_(consumes<LHEEventProduct>(lheEventTag_)),
       lheRunToken_(consumes<LHERunInfoProduct,edm::InRun>(lheEventTag_)),
   pdfWeightsToken_(consumes<vector<double>>(pdfWeightsTag_))
{
}
GenEventBlock::~GenEventBlock() {
  delete list_;
}
void GenEventBlock::beginJob()
{
  // Get TTree pointer
  TTree* tree = vhtm::Utility::getTree("vhtree");

  list_ = new vector<vhtm::GenEvent>();
  tree->Branch("GenEvent", "vector<vhtm::GenEvent>", &list_, 32000, -1);
}
void GenEventBlock::beginRun(edm::Run const& iRun, edm::EventSetup const& iSetup) {
  pdf_indices.clear();
  scale_indices.clear();
  alpha_indices.clear();

  edm::Handle<LHERunInfoProduct> lheRunHandle; 
  iRun.getByLabel(lheEventTag_, lheRunHandle);
  //bool found = iRun.getByToken(lheRunToken_, lheRunHandle);
  
  if (lheRunHandle.isValid()) { // found &&
    const LHERunInfoProduct& lheRunInfo = *(lheRunHandle.product());
    //--- get info from LHERunInfoProduct
    vector<string> weight_lines;
    for (vector<LHERunInfoProduct::Header>::const_iterator iter  = lheRunInfo.headers_begin(); 
                                                           iter != lheRunInfo.headers_end();
   	                                                 ++iter) 
    {
      vector<string> lines = iter->lines();
      if ((iter->tag()).compare(tag_) == 0)
	weight_lines = iter->lines();

      for (const auto& l: lines) {
	if (l.find("powheg") != string::npos) {
	  generatorType_ = "powheg";
	  break;
	}
      }
    }
    // --- Get the PDF ids -
    // See Josh's slides 13-15: https://indico.cern.ch/event/459797/contribution/2/attachments/1181555/1800214/mcaod-Feb15-2016.pdf
    int pdfidx = lheRunInfo.heprup().PDFSUP.first;
    if (pdfidx == -1 && generatorType_ == "powheg") pdfidx = 260000;    
    cout << "This sample was generated with the following PDFs: " << pdfidx << endl;
    
    // --- Get min and max pdf index for 100 replicas
    if (!isStandardSample_) {pdfidx = 0;} 
    pdfid_1 = boost::lexical_cast<string>(pdfidx + 1);
    pdfid_2 = boost::lexical_cast<string>(pdfidx + 100);
    cout << "PDFs min and max id for MC replicas: " << pdfid_1 << "   " << pdfid_2 <<endl;

    // --- Get alphas id
    if (doAlphasWeights_) {
      switch (pdfidx) {
      case 292200:
	alphas_id_1 = "292301";
	alphas_id_2 = "292302";
	break;
      case 292000:
	alphas_id_1 = "292101";
	alphas_id_2 = "292102";
	break;
      case 260000:
	alphas_id_1 = "265000";
	alphas_id_2 = "266000";
	break;
      case 260400:
	alphas_id_1 = "265400";
	alphas_id_2 = "266400";
	break;
      default:
	cout << "Info: pdfidx: " << pdfidx << endl;
	break;
      }
      cout << "alpha_s min and max id: " << alphas_id_1 << "   " << alphas_id_2 << endl;
    }

    if (!weight_lines.size()) {
      cerr << "==> lheRunInfo header may not contain required information! " << endl;
      return;
    }
    // --- for amcatnlo samples, remove last line
    if (weight_lines.back() == "  <")
      weight_lines.pop_back();
                
    // --- Convert weight_lines in stringstream object and populate boost::property_tree
    stringstream ss;
    std::copy(weight_lines.begin(), weight_lines.end(), std::ostream_iterator<string>(ss,""));
    cout << ss.str() << endl;
    boost::property_tree::ptree pt;
    read_xml(ss, pt);

    // --- Name of the weightgroup
    string scalevar = "scale_variation";
    string pdfvar   = "PDF_variation";
    if (!isStandardSample_) {
      pdfvar = pdfset_;
      scalevar = "Central scale variation";
    }
    // --- Loop over elements to get PDF, alpha_s and scale weights
    BOOST_FOREACH(boost::property_tree::ptree::value_type const& v, pt.get_child("")) {
      if (v.first == "weightgroup") {
	cout << "First data: " << v.first.data() << endl;
                
	boost::property_tree::ptree subtree = (boost::property_tree::ptree) v.second ;
	boost::optional<string> weightgroupname1 = v.second.get_optional<string>("<xmlattr>.name");
	boost::optional<string> weightgroupname2 = v.second.get_optional<string>("<xmlattr>.type");
        
	if (weightgroupname1) cout << weightgroupname1.get() << endl;
	if (weightgroupname2) cout << weightgroupname2.get() << endl;
	
	// -- PDFs + alpha_s weights
	if ((weightgroupname1 && weightgroupname1.get() == pdfvar) || (weightgroupname2 && weightgroupname2.get() == pdfvar)) {               
	  BOOST_FOREACH(boost::property_tree::ptree::value_type &vs,subtree)
	    if (vs.first == "weight") {
	      cout << vs.first <<  "   " << vs.second.get<string>("<xmlattr>.id")  << "  " << vs.second.data() << endl;
	      
	      string strwid = vs.second.get<string>("<xmlattr>.id");
	      string strw   = vs.second.data();
	      
	      int id = stoi(strwid);
	      vector<string> strs;
	      if (isStandardSample_) {
		boost::split(strs, strw, boost::is_any_of("="));
	      }
	      else {
		boost::split(strs, strw, boost::is_any_of(" "));
	      }
	      
	      int pdfindex = stoi(strs.back());
	      cout << "pdfindex " << pdfindex << endl;
	      
	      if (pdfindex >= stoi(pdfid_1) && pdfindex <= stoi(pdfid_2)) pdf_indices.push_back(id);
	      if (doAlphasWeights_) {
		if (pdfindex == stoi(alphas_id_1) || pdfindex == stoi(alphas_id_2)) alpha_indices.push_back(id);
	      }
	    }
	}
	// -- Scale weights
	if ((weightgroupname1 && weightgroupname1.get() == scalevar) || (weightgroupname2 && weightgroupname2.get() == scalevar)) {
	  BOOST_FOREACH(boost::property_tree::ptree::value_type &vs, subtree)
	    if (vs.first == "weight") {
	      string strwid = vs.second.get<string>("<xmlattr>.id");
	      int id = stoi(strwid);
	      scale_indices.push_back(id);
	    }
	}
      }
    }
  }
  else {
    cerr << "Failed to retrieve LHERunInfoProduct with tag: " << lheEventTag_ << endl;
  }
}
void GenEventBlock::analyze(edm::Event const& iEvent, edm::EventSetup const& iSetup) {
  // Reset the vector
  list_->clear();

  if (!iEvent.isRealData()) {
    // Create Event Object
    vhtm::GenEvent ev;

    // GenEventInfo Part
    edm::Handle<GenEventInfoProduct> genEvtHandle;
    bool found = iEvent.getByToken(genEventToken_, genEvtHandle);
    if (found && genEvtHandle.isValid()) {
      edm::LogInfo("GenEventBlock") << "Success. Obtained GenEventInfoProduct for label: "
                                    << genEventTag_;
      ev.processID = genEvtHandle->signalProcessID();
      ev.ptHat     = genEvtHandle->hasBinningValues()
                         ? genEvtHandle->binningValues()[0] : 0.;
      ev.nMEPartons = genEvtHandle->nMEPartons();
      ev.qScale   = genEvtHandle->qScale();
      ev.alphaQCD = genEvtHandle->alphaQCD();
      ev.alphaQED = genEvtHandle->alphaQED();
      if (genEvtHandle->pdf()) {
	ev.qScalePdf = genEvtHandle->pdf()->scalePDF;
	ev.x1Pdf     = genEvtHandle->pdf()->x.first;
	ev.x2Pdf     = genEvtHandle->pdf()->x.second;
	ev.id1Pdf    = genEvtHandle->pdf()->id.first;
	ev.id2Pdf    = genEvtHandle->pdf()->id.second;
      }

      double evtWeight  = genEvtHandle->weight();
      ev.evtWeight = evtWeight;

      edm::Handle<LHEEventProduct> lheEvtHandle;
      found = iEvent.getByToken(lheEventToken_, lheEvtHandle);
      if (found && lheEvtHandle.isValid()) {
        inpdfweights.clear(); 
        for (size_t i = 0; i < lheEvtHandle->weights().size(); ++i) {
	  int id_i = stoi(lheEvtHandle->weights()[i].id);
	  
	  // --- Get PDF weights
	  for (size_t j = 0; j < pdf_indices.size(); ++j) {
	    int id_j = pdf_indices[j];
	    if (id_i == id_j) {
	      double weight = lheEvtHandle->weights()[i].wgt;
	      inpdfweights.push_back(weight);
	    }
	  }
            
	  // --- Get alpha_s weights
	  if (doAlphasWeights_) {
	    for (size_t k = 0; k < alpha_indices.size(); ++k) {
	      int id_k = alpha_indices[k];
	      if (id_i == id_k) {
		double wt = lheEvtHandle->weights()[i].wgt;
		ev.alphasWeightList.push_back(wt);
	      }
	    }
	  }
            
	  // --- Get qcd scale weights
	  if (doScaleWeights_) { 
	    for (size_t k = 0 ; k < scale_indices.size(); ++k) {
	      int id_k = scale_indices[k];
	      if (id_i == id_k) {
		double wt = lheEvtHandle->weights()[i].wgt;
		ev.qcdScaleWeightList.push_back(wt);
	      }
	    }
	  }
	}
        if (pdf_indices.size() > 0) {
	  // --- Get MCtoHessian PDF weights
	  pdfWeightsHelper_.Init(pdf_indices.size(), nPdfEigWeights_, mc2hessianCSV_);
	  vector<double> outpdfweights(nPdfEigWeights_);
        
	  double nomlheweight = lheEvtHandle->weights()[0].wgt;
	  pdfWeightsHelper_.DoMC2Hessian(nomlheweight, inpdfweights.data(), outpdfweights.data());
	  for (size_t i = 0; i < nPdfEigWeights_; ++i) {
	    double wgtval = outpdfweights[i];
	    float real_weight = wgtval * evtWeight/nomlheweight;
            
	    // the is the weight to be used for evaluating uncertainties with hessian weights
	    ev.pdfWeightList.push_back(real_weight);
	  }    
	}
        // https://github.com/KappaAnalysis/Kappa
	// Get generator level HT
	double lheHt = 0.;
	int lheNOutPartons = 0;
	const lhef::HEPEUP& lheEvent = lheEvtHandle->hepeup();
	const std::vector<lhef::HEPEUP::FiveVector>& lheParticles = lheEvent.PUP;
	for (size_t i = 0; i < lheParticles.size(); ++i) {
	  int id = std::abs(lheEvent.IDUP[i]);
	  int status = lheEvent.ISTUP[i];
	  if (status == 1 && ((id >= 1 && id <= 6) || id == 21)) { // quarks and gluons
	    lheHt += std::sqrt(std::pow(lheParticles[i][0], 2.) + std::pow(lheParticles[i][1], 2.));
	    ++lheNOutPartons;
	  }
	}
	ev.lheHt = lheHt;
	ev.lheNOutPartons = lheNOutPartons;
      }
      else {
        if (verbosity_ > 0)
        edm::LogError("GenEventBlock") << "Error! Failed to get LHEEventProduct for label: "
                                     << lheEventTag_;
      }
    }
    else {
      if (verbosity_ > 0)
      edm::LogError("GenEventBlock") << "Error! Failed to get GenEventInfoProduct for label: "
                                     << genEventTag_;
    }
    // PDF Weights Part
    if (storePDFWeights_) {
      edm::Handle<vector<double>> pdfWeightsHandle;
      found = iEvent.getByToken(pdfWeightsToken_, pdfWeightsHandle);

      if (found && pdfWeightsHandle.isValid()) {
        if (verbosity_ > 0)
        edm::LogInfo("GenEventBlock") << "Success. Obtained PDF handle for label: "
                                      << pdfWeightsTag_;
	copy(pdfWeightsHandle->begin(), pdfWeightsHandle->end(), ev.pdfWeights.begin());
      }
      else {
        if (verbosity_ > 0)
        edm::LogError("GenEventBlock") << "Error! Failed to get PDF handle for label: "
                                       << pdfWeightsTag_;
      }
    }
    list_->push_back(ev);
  }
}
#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(GenEventBlock);
