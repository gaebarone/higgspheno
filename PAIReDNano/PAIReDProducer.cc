#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/StreamID.h"

#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/Candidate/interface/CandidateFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include "RecoBTag/FeatureTools/interface/deep_helpers.h"
#include "DataFormats/Candidate/interface/VertexCompositePtrCandidate.h"
using namespace btagbtvdeep;

#include "CommonTools/Utils/interface/StringCutObjectSelector.h"
#include "DataFormats/NanoAOD/interface/FlatTable.h"

#include "inEllipse.h"

template<typename T>
class PAIReDJetProducer : public edm::stream::EDProducer<> {
public:
  explicit PAIReDJetProducer(const edm::ParameterSet &);
  ~PAIReDJetProducer() override;

  static void fillDescriptions(edm::ConfigurationDescriptions &descriptions);

private:
  void produce(edm::Event &, const edm::EventSetup &) override;

  //const std::string name_;
  const std::string name_;

  edm::EDGetTokenT<edm::View<T>> jet_token_;
  edm::EDGetTokenT<reco::CandidateView> cand_token_;
  
};

//
// constructors and destructor
//
template< typename T>
PAIReDJetProducer<T>::PAIReDJetProducer(const edm::ParameterSet &iConfig)
    : name_(iConfig.getParameter<std::string>("name")),
      jet_token_(consumes<edm::View<T>>(iConfig.getParameter<edm::InputTag>("jets"))),
      cand_token_(consumes<reco::CandidateView>(iConfig.getParameter<edm::InputTag>("candidates"))){
  produces<nanoaod::FlatTable>(name_);
}

template< typename T>
PAIReDJetProducer<T>::~PAIReDJetProducer() {}

template< typename T>
void PAIReDJetProducer<T>::produce(edm::Event &iEvent, const edm::EventSetup &iSetup) {
  // elements in all these collections must have the same order!
  std::vector<int> idx_jet1, idx_jet2, true_composition;
  std::vector<float> bb_score, cc_score, other_score;
  std::vector<float> regression_pt, regression_mass, regression_phi, regression_eta;

  auto jets = iEvent.getHandle(jet_token_);
  auto cands = iEvent.getHandle(cand_token_);

  for (unsigned i_jet = 0; i_jet < jets->size(); ++i_jet) {
    const auto &jet1 = jets->at(i_jet);
    idx_jet1.push_back(i_jet);
    for (unsigned j_jet = i_jet+1; j_jet < jets->size(); ++j_jet) {
      auto paired_jet = std::make_unique<std::vector<reco::CandidatePtr>>();
      const auto &jet2 = jets->at(j_jet);
      idx_jet2.push_back(j_jet);
      for (unsigned i_cand = 0; i_cand < cands->size(); ++i_cand) {
        const auto &cand = cands->at(i_cand);
        if inEllipse(jet1, jet2, cand){
          paired_jet.push_back(cand);
        }
      }
      // get truth

      // get training scores and regression information
      bb_score.push_back(/*insert score*/);
      cc_score.push_back(/*insert score*/);
      other_score.push_back(/*insert score*/);
      regression_pt.push_back(/*insert result*/);
      regression_mass.push_back(/*insert result*/);
      regression_phi.push_back(/*insert result*/);
      regression_eta.push_back(/*insert result*/);
    }
  }

  auto pjTable = std::make_unique<nanoaod::FlatTable>(idx_jet1->size(), name_, false);
  // We fill from here only stuff that cannot be created with the SimpleFlatTableProducer
  pjTable->addColumn<int>("idx_jet1", idx_jet1, "Index of constituent jet 1", nanoaod::FlatTable::IntColumn);
  pjTable->addColumn<int>("idx_jet2", idx_jet2, "Index of constituent jet 2", nanoaod::FlatTable::IntColumn);
  pjTable->addColumn<int>("true_composition", true_composition, "Truth value for jet (0 other, 1 bb, 2 cc)", nanoaod::FlatTable::IntColumn);
  pjTable->addColumn<float>("regression_pt", regression_pt, "pt of paired jet", nanoaod::FlatTable::FloatColumn, 10);
  pjTable->addColumn<float>("regression_mass", regression_mass, "mass of paired jet", nanoaod::FlatTable::FloatColumn, 10);
  pjTable->addColumn<float>("regression_eta", regression_eta, "eta of paired jet", nanoaod::FlatTable::FloatColumn, 10);
  pjTable->addColumn<float>("regression_phi", regression_phi, "phi of paired jet", nanoaod::FlatTable::FloatColumn, 10);
  pjTable->addColumn<float>("bb_score", bb_score, "Model score for bb jet", nanoaod::FlatTable::FloatColumn, 10);
  pjTable->addColumn<float>("cc_score", cc_score, "Model score for cc jet", nanoaod::FlatTable::FloatColumn, 10);
  pjTable->addColumn<float>("other_score", other_score, "Model score for other jet", nanoaod::FlatTable::FloatColumn, 10);

  iEvent.put(std::move(pjTable), name_);

}


template< typename T>
void PAIReDProducer<T>::fillDescriptions(edm::ConfigurationDescriptions &descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<std::string>("name", "PAIReDJets");
  desc.add<edm::InputTag>("jets", edm::InputTag("slimmedJetsAK4"));
  desc.add<edm::InputTag>("candidates", edm::InputTag("packedPFCandidates"));
  descriptions.addWithDefaultLabel(desc);
}

typedef JetConstituentTableProducer<pat::Jet> PAIReDJetProducer;
typedef JetConstituentTableProducer<reco::GenJet> PAIReDGenJetProducer;

DEFINE_FWK_MODULE(PAIReDJetProducer);
DEFINE_FWK_MODULE(PAIReDGenJetProducer);