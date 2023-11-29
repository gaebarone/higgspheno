// RecoBTag/FeatureTools/plugins/PAIReDTagInfoProducer.cc
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/StreamID.h"
#include "FWCore/Utilities/interface/ESGetToken.h"

#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"

#include "DataFormats/BTauReco/interface/ShallowTagInfo.h"

#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"

#include "RecoBTag/FeatureTools/interface/JetConverter.h"
#include "RecoBTag/FeatureTools/interface/ShallowTagInfoConverter.h"
#include "RecoBTag/FeatureTools/interface/SecondaryVertexConverter.h"
#include "RecoBTag/FeatureTools/interface/NeutralCandidateConverter.h"
#include "RecoBTag/FeatureTools/interface/ChargedCandidateConverter.h"

#include "RecoBTag/FeatureTools/interface/TrackInfoBuilder.h"
#include "RecoBTag/FeatureTools/interface/sorting_modules.h"

#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/Candidate/interface/VertexCompositePtrCandidate.h"

#include "RecoBTag/FeatureTools/interface/deep_helpers.h"

#include "FWCore/ParameterSet/interface/Registry.h"
#include "FWCore/Common/interface/Provenance.h"
#include "DataFormats/Provenance/interface/ProductProvenance.h"

#include "DataFormats/BTauReco/interface/SeedingTrackFeatures.h"
#include "DataFormats/BTauReco/interface/TrackPairFeatures.h"
#include "RecoBTag/FeatureTools/interface/TrackPairInfoBuilder.h"
#include "RecoBTag/FeatureTools/interface/SeedingTrackInfoBuilder.h"
#include "RecoBTag/FeatureTools/interface/SeedingTracksConverter.h"

#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "RecoBTag/TrackProbability/interface/HistogramProbabilityEstimator.h"
class HistogramProbabilityEstimator;
#include <memory>

#include <typeinfo>
#include "CondFormats/BTauObjects/interface/TrackProbabilityCalibration.h"
#include "CondFormats/DataRecord/interface/BTagTrackProbability2DRcd.h"
#include "CondFormats/DataRecord/interface/BTagTrackProbability3DRcd.h"
#include "FWCore/Framework/interface/EventSetupRecord.h"
#include "FWCore/Framework/interface/EventSetupRecordImplementation.h"
#include "FWCore/Framework/interface/EventSetupRecordKey.h"

#include "DataFormats/BTauReco/interface/SecondaryVertexFeatures.h"

#include "DataFormats/JetReco/interface/Jet.h"
#include "DataFormats/VertexReco/interface/Vertex.h"

#include "inEllipse.h"

#include "FWCore/Framework/interface/makeRefToBaseProdFrom.h"

#include "DataFormats/BTauReco/interface/JetTag.h"

#include "PhysicsTools/ONNXRuntime/interface/ONNXRuntime.h"

#include "DataFormats/NanoAOD/interface/FlatTable.h"

#include "testPAIReD.h"

using namespace cms::Ort;
using namespace btagbtvdeep;

class PAIReDONNXProducer : public edm::stream::EDProducer<edm::GlobalCache<ONNXRuntime>> {
public:
  explicit PAIReDONNXProducer(const edm::ParameterSet&, const ONNXRuntime*);
  ~PAIReDONNXProducer() override = default;

  static void fillDescriptions(edm::ConfigurationDescriptions&);

  static std::unique_ptr<ONNXRuntime> initializeGlobalCache(const edm::ParameterSet&);
  static void globalEndJob(const ONNXRuntime*);

private:
  typedef reco::JetTagCollection JetTagCollection;
  typedef reco::VertexCompositePtrCandidateCollection SVCollection;
  typedef reco::VertexCollection VertexCollection;

  const std::string name_;
  const edm::EDGetTokenT<edm::View<reco::Jet>> jet_token_;
  edm::EDGetTokenT<reco::CandidateView> cand_token_;
  edm::EDGetTokenT<VertexCollection> vtx_token_;
  const edm::EDGetTokenT<SVCollection> sv_token_;
  const edm::ESGetToken<TransientTrackBuilder, TransientTrackRecord> track_builder_token_;
  //edm::EDGetTokenT<edm::ValueMap<float>> puppi_value_map_token_;
  std::vector<std::string> input_names_;
  std::vector<std::string> output_names_;

  void beginStream(edm::StreamID) override {}
  void produce(edm::Event &, const edm::EventSetup &) override;
  void make_inputs(edm::Event& iEvent, const edm::EventSetup& iSetup, edm::Handle<edm::View<reco::Jet>> jets, unsigned i_jet, unsigned j_jet, edm::Handle<SVCollection> svs, edm::Handle<edm::View<reco::Candidate>> cands, edm::Handle<VertexCollection> vtxs);
  void get_input_sizes(edm::Event& iEvent, const edm::EventSetup& iSetup, edm::Handle<edm::View<reco::Jet>> jets, unsigned i_jet, unsigned j_jet, edm::Handle<SVCollection> svs, edm::Handle<edm::View<reco::Candidate>> cands, edm::Handle<VertexCollection> vtxs);
  void endStream() override {}

  enum InputIndexes {
    kChargedCandidates = 0,
    kNeutralCandidates = 1,
    kVertices = 2,
    kChargedCandidates4Vec = 3,
    kNeutralCandidates4Vec = 4,
    kVertices4Vec = 5
  };
  unsigned n_cpf_;
  constexpr static unsigned n_features_cpf_ = 27;
  constexpr static unsigned n_pairwise_features_cpf_ = 4;
  unsigned n_npf_;
  constexpr static unsigned n_features_npf_ = 12;
  constexpr static unsigned n_pairwise_features_npf_ = 4;
  unsigned n_sv_;
  constexpr static unsigned n_features_sv_ = 18;
  constexpr static unsigned n_pairwise_features_sv_ = 4;
  std::vector<unsigned> input_sizes_;
  std::vector<std::vector<int64_t>> input_shapes_;  // shapes of each input group (-1 for dynamic axis)

  // hold the input data
  FloatArrays data_;
  // TEST STUFF
  std::vector<std::vector<float>> test_data;
  std::vector<std::string> test_data_names = {"pf_cand_e", "pf_cand_px", "pf_cand_py", "pf_cand_pz"};
  std::vector<float> pf_cand_es;
  std::vector<float> pf_cand_pxs;
  std::vector<float> pf_cand_pys;
  std::vector<float> pf_cand_pzs;
  bool fill_test = true;
  // END OF TEST STUFF
};

PAIReDONNXProducer::PAIReDONNXProducer(const edm::ParameterSet& iConfig, const ONNXRuntime* cache)
    : name_(iConfig.getParameter<std::string>("name")),
      jet_token_(consumes<edm::View<reco::Jet>>(iConfig.getParameter<edm::InputTag>("jets"))),
      cand_token_(consumes<reco::CandidateView>(iConfig.getParameter<edm::InputTag>("candidates"))),
      vtx_token_(consumes<VertexCollection>(iConfig.getParameter<edm::InputTag>("vertices"))),
      sv_token_(consumes<SVCollection>(iConfig.getParameter<edm::InputTag>("secondary_vertices"))), 
      track_builder_token_(esConsumes<TransientTrackBuilder, TransientTrackRecord>(edm::ESInputTag("", "TransientTrackBuilder"))){
      //puppi_value_map_token_(consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("puppi_value_map"))){
  produces<nanoaod::FlatTable>(name_);
}

void PAIReDONNXProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<std::string>("name", "PAIReDJets");
  desc.add<edm::InputTag>("jets", edm::InputTag("ak4PFJetsCHS"));
  desc.add<edm::InputTag>("candidates", edm::InputTag("packedPFCandidates"));
  desc.add<edm::InputTag>("vertices", edm::InputTag("offlineSlimmedPrimaryVertices"));
  desc.add<edm::InputTag>("secondary_vertices", edm::InputTag("slimmedSecondaryVertices"));
  //desc.add<edm::InputTag>("puppi_value_map", edm::InputTag("puppi"));
  //desc.add<std::vector<std::string>>("input_names", {"charged_candidate_features", "neutral_candidate_features", "secondary_vertex_features", "charged_candidate_4vec", "neutral_candidate_4vec", "secondary_vertex_4vec"});
  //desc.add<std::vector<std::string>>("output_names", {"probbb", "probbcc", "probudsg", "reg_mass", "reg_pt", "reg_eta", "reg_phi"});
  desc.add<std::vector<std::string>>("input_names", {"pf_points", "pf_features", "pf_vectors", "pf_mask"});
  desc.add<std::vector<std::string>>("output_names", {"label_ll", "label_cc", "label_bb"});
  desc.add<edm::FileInPath>("model_path", edm::FileInPath("ProdTutorial/PAIReDONNXProducer/ONNX/PAIReDEllipse_ParT_20230520134246.onnx"));
  descriptions.add("PAIReDJetTable", desc);
}

std::unique_ptr<ONNXRuntime> PAIReDONNXProducer::initializeGlobalCache(
    const edm::ParameterSet& iConfig) {
  return std::make_unique<ONNXRuntime>(iConfig.getParameter<edm::FileInPath>("model_path").fullPath());
}

void PAIReDONNXProducer::globalEndJob(const ONNXRuntime* cache) {}

void PAIReDONNXProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
  // initialize output variables
  std::vector<int> idx_jet1, idx_jet2, true_composition;
  std::vector<float> bb_score, cc_score, other_score, ll_score;
  std::vector<float> regression_pt, regression_mass, regression_phi, regression_eta;
  edm::Handle<reco::VertexCollection> vtxs;
  iEvent.getByToken(vtx_token_, vtxs);
  edm::Handle<edm::View<reco::Jet>> jets;
  iEvent.getByToken(jet_token_, jets);
  edm::Handle<edm::View<reco::Candidate>> cands;
  iEvent.getByToken(cand_token_, cands);

  edm::Handle<reco::VertexCompositePtrCandidateCollection> svs;
  iEvent.getByToken(sv_token_, svs);

  for (unsigned i_jet = 0; i_jet < jets->size(); ++i_jet) {
    for (unsigned j_jet = i_jet+1; j_jet < jets->size(); ++j_jet) {
      std::vector<float> outputs(output_names_.size(), -1.0);
      // store constituent jets and find truth value
      idx_jet1.emplace_back(i_jet);
      idx_jet2.emplace_back(j_jet);
      get_input_sizes(iEvent, iSetup, jets, i_jet, j_jet, svs, cands, vtxs);
      // run prediction with dynamic batch size per event
      /*input_shapes_ = {{(int64_t)1, (int64_t)n_cpf_, (int64_t)n_features_cpf_},
                       {(int64_t)1, (int64_t)n_npf_, (int64_t)n_features_npf_},
                       {(int64_t)1, (int64_t)n_sv_, (int64_t)n_features_sv_},
                       {(int64_t)1, (int64_t)n_cpf_, (int64_t)n_pairwise_features_cpf_},
                       {(int64_t)1, (int64_t)n_npf_, (int64_t)n_pairwise_features_npf_},
                       {(int64_t)1, (int64_t)n_sv_, (int64_t)n_pairwise_features_sv_}};*/
      input_shapes_ = {{(int64_t)1, (int64_t)n_npf_, (int64_t)4},
                       {(int64_t)1, (int64_t)n_npf_, (int64_t)16},
                       {(int64_t)1, (int64_t)n_npf_, (int64_t)4},
                       {(int64_t)1, (int64_t)n_npf_, (int64_t)1}};
      outputs = globalCache()->run(input_names_, data_, input_shapes_)[0];
      assert(outputs.size() == output_names_.size());
      // store data
      ll_score.emplace_back(outputs[0]);
      cc_score.emplace_back(outputs[1]);
      bb_score.emplace_back(outputs[2]);
      /*bb_score.emplace_back(outputs[0]);
      cc_score.emplace_back(outputs[1]);
      other_score.emplace_back(outputs[2]);
      regression_mass.emplace_back(outputs[3]);
      regression_pt.emplace_back(outputs[4]);
      regression_eta.emplace_back(outputs[5]);
      regression_phi.emplace_back(outputs[6]);*/
      fill_test = false; // TEST STUFF
    }
  }
  saveParams(test_data, test_data_names); // TEST STUFF
  auto pjTable = std::make_unique<nanoaod::FlatTable>(idx_jet1.size(), name_, false);
  /*pjTable->addColumn<int>("idx_jet1", idx_jet1, "Index of constituent jet 1");
  pjTable->addColumn<int>("idx_jet2", idx_jet2, "Index of constituent jet 2");
  pjTable->addColumn<float>("regression_mass", regression_mass, "mass of paired jet", 10);
  pjTable->addColumn<float>("regression_pt", regression_pt, "pt of paired jet", 10);
  pjTable->addColumn<float>("regression_eta", regression_eta, "eta of paired jet", 10);
  pjTable->addColumn<float>("regression_phi", regression_phi, "phi of paired jet", 10);
  pjTable->addColumn<float>("bb_score", bb_score, "Model score for bb jet", 10);
  pjTable->addColumn<float>("cc_score", cc_score, "Model score for cc jet", 10);
  pjTable->addColumn<float>("other_score", other_score, "Model score for other jet", 10);*/

  pjTable->addColumn<float>("ll_score", ll_score, "Model score for ll jet", 10);
  pjTable->addColumn<float>("cc_score", cc_score, "Model score for cc jet", 10);
  pjTable->addColumn<float>("bb_score", bb_score, "Model score for bb jet", 10);

  iEvent.put(std::move(pjTable), name_);
}

void PAIReDONNXProducer::get_input_sizes(edm::Event& iEvent, 
                                                const edm::EventSetup& iSetup,
                                                edm::Handle<edm::View<reco::Jet>> jets, 
                                                unsigned i_jet, 
                                                unsigned j_jet, 
                                                edm::Handle<SVCollection> svs, 
                                                edm::Handle<edm::View<reco::Candidate>> cands, 
                                                edm::Handle<VertexCollection> vtxs) {
  unsigned int num_svs = 0;
  unsigned int num_cpfs = 0;
  unsigned int num_npfs = 0;
  const auto &jet1 = jets->at(i_jet);
  const auto &jet2 = jets->at(j_jet);
  // loop through secondary vertices and count those in the ellipse
  for (const auto& sv : *svs) {
    if (inEllipse(jet1.eta(), jet1.phi(), jet2.eta(), jet2.phi(), sv.eta(), sv.phi())) num_svs += 1;
  }
  // loop through particle flow candidates and count those in the ellipse
  for (unsigned i_cand = 0; i_cand < cands->size(); ++i_cand) {
    const reco::Candidate *cand = &(cands->at(i_cand));
    if (!inEllipse(jet1.eta(), jet1.phi(), jet2.eta(), jet2.phi(), (*cand).eta(), (*cand).phi())) continue;
    if (cand->charge() != 0) num_cpfs += 1;
    else num_npfs += 1;
  }
  // adjust counting as needed and make data the right size
  n_cpf_ = std::max((unsigned int)1, num_cpfs);
  n_npf_ = std::max((unsigned int)1, num_npfs);
  n_sv_ = std::max((unsigned int)1, num_svs);
  //n_cpf_ = std::min((unsigned int)25, n_cpf_);
  //n_npf_ = std::min((unsigned int)25, n_npf_);
  //n_sv_ = std::min((unsigned int)5, n_sv_);
  /*input_sizes_ = {
      n_cpf_ * n_features_cpf_,
      n_npf_ * n_features_npf_,
      n_sv_ * n_features_sv_,
      n_cpf_ * n_pairwise_features_cpf_,
      n_npf_ * n_pairwise_features_npf_,
      n_sv_ * n_pairwise_features_sv_,
  };*/

  input_sizes_ = {
      n_npf_ * 4,
      n_npf_ * 16,
      n_npf_ * 4,
      n_npf_ * 1,
  };
  // init data storage
  data_.clear();
  for (const auto& len : input_sizes_) {
    data_.emplace_back(1 * len, 0);
  }
  // generate the inputs
  make_inputs(iEvent, iSetup, jets, i_jet, j_jet, svs, cands, vtxs);
}

void PAIReDONNXProducer::make_inputs(edm::Event& iEvent, 
                                           const edm::EventSetup& iSetup,
                                           edm::Handle<edm::View<reco::Jet>> jets, 
                                           unsigned i_jet, 
                                           unsigned j_jet, 
                                           edm::Handle<SVCollection> svs, 
                                           edm::Handle<edm::View<reco::Candidate>> cands, 
                                           edm::Handle<VertexCollection> vtxs) {
  // initialize input facilitation
  float* ptr = nullptr;
  //float* start = nullptr;
  unsigned offset = 0;  
  // now begin to process and fill data
  const auto& pv = vtxs->at(0);
  double jet_radius_ = 0.4;
  // get primary jet data
  const auto &jet1 = jets->at(i_jet);
  edm::RefToBase<reco::Jet> jet_ref1(jets, i_jet);
  math::XYZVector jet_dir1 = jet1.momentum().Unit();
  GlobalVector jet_ref_track_dir1(jet1.px(), jet1.py(), jet1.pz());
  // get secondary jet data
  const auto &jet2 = jets->at(j_jet);
  edm::RefToBase<reco::Jet> jet_ref2(jets, j_jet);
  math::XYZVector jet_dir2 = jet2.momentum().Unit();
  GlobalVector jet_ref_track_dir2(jet2.px(), jet2.py(), jet2.pz());
  // get puppi weight and track builder
  //edm::Handle<edm::ValueMap<float>> puppi_value_map;
  //iEvent.getByToken(puppi_value_map_token_, puppi_value_map);
  edm::ESHandle<TransientTrackBuilder> track_builder = iSetup.getHandle(track_builder_token_);
  // sort secondary vertices by dxy
  auto svs_sorted = *svs;
  std::sort(svs_sorted.begin(), svs_sorted.end(), [&pv](const auto& sva, const auto& svb) {
    return btagbtvdeep::sv_vertex_comparator(sva, svb, pv);
  });
  const auto& svs_unsorted = *svs;
  // loop through secondary vertices
  /*
  std::size_t sv_n = 0;
  for (const auto& sv : svs_sorted) {
    if (!inEllipse(jet1.eta(), jet1.phi(), jet2.eta(), jet2.phi(), sv.eta(), sv.phi())) continue;
    // save secondary vertex variables
    ptr = &data_[kVertices][offset + sv_n * n_features_sv_];
    start = ptr;
    *ptr = catch_infs_and_bound(std::fabs(reco::deltaR(sv, jet_dir1)) - 0.5, 0, -2, 0); // deltaR 1
    *(++ptr) = catch_infs_and_bound(std::fabs(reco::deltaR(sv, jet_dir2)) - 0.5, 0, -2, 0); // deltaR 2
    *(++ptr) = catch_infs_and_bound(std::fabs(sv.eta() - jet1.eta()) - 0.5, 0, -2, 0); // etarel 1
    *(++ptr) = catch_infs_and_bound(std::fabs(sv.eta() - jet2.eta()) - 0.5, 0, -2, 0); // etarel 2
    *(++ptr) = catch_infs_and_bound(std::fabs(reco::deltaPhi(sv.phi(), jet1.phi())) - 0.5, 0, -2, 0); // phirel 1
    *(++ptr) = catch_infs_and_bound(std::fabs(reco::deltaPhi(sv.phi(), jet2.phi())) - 0.5, 0, -2, 0); // phirel 2
    *(++ptr) = sv.energy() / jet1.energy(); // en ratio 1
    *(++ptr) = sv.energy() / jet2.energy(); // en ratio 2
    *(++ptr) = sv.pt(); // pt
    *(++ptr) = sv.numberOfDaughters(); // ntracks
    *(++ptr) = sv.vertexChi2(); // chi2
    *(++ptr) = catch_infs_and_bound(sv.vertexChi2() / sv.vertexNdof(), 1000, -1000, 1000); // norm chi2
    *(++ptr) = sv.mass(); // mass
    const auto& dxy_meas = vertexDxy(sv, pv);
    *(++ptr) = dxy_meas.value(); // dxy
    *(++ptr) = catch_infs_and_bound(dxy_meas.value() / dxy_meas.error(), 0, -1, 800); // dxy sig
    const auto& d3d_meas = vertexD3d(sv, pv);
    *(++ptr) = d3d_meas.value(); // d3d
    *(++ptr) = catch_infs_and_bound(d3d_meas.value() / d3d_meas.error(), 0, -1, 800); // d3d sig
    *(++ptr) = vertexDdotP(sv, pv); // cos theta sv pv
    assert(start + n_features_sv_ - 1 == ptr);
    // save secondary vertex four vector
    ptr = &data_[kVertices4Vec][offset + sv_n * n_pairwise_features_sv_];
    start = ptr;
    *ptr = sv.px();
    *(++ptr) = sv.py();
    *(++ptr) = sv.pz();
    *(++ptr) = sv.energy();
    assert(start + n_pairwise_features_sv_ - 1 == ptr);
    sv_n += 1;
  }
  */
  // sort pf candidates
  std::vector<btagbtvdeep::SortingClass<size_t>> c_sorted_j1, n_sorted_j1, c_sorted_j2, n_sorted_j2;
  std::map<unsigned int, btagbtvdeep::TrackInfoBuilder> trackinfos_j1, trackinfos_j2;
  for (unsigned i_cand = 0; i_cand < cands->size(); ++i_cand) {
    const reco::Candidate *cand = &(cands->at(i_cand));
    // note that here we already ensure the candidate is within the ellipse
    if (!inEllipse(jet1.eta(), jet1.phi(), jet2.eta(), jet2.phi(), (*cand).eta(), (*cand).phi())) continue;
    if (cand->charge() != 0) {
      // build track info and keep charged candidates
      auto& track_info_j1 = trackinfos_j1.emplace(i_cand, track_builder).first->second;
      auto& track_info_j2 = trackinfos_j2.emplace(i_cand, track_builder).first->second;
      track_info_j1.buildTrackInfo(cand, jet_dir1, jet_ref_track_dir1, pv);
      track_info_j2.buildTrackInfo(cand, jet_dir2, jet_ref_track_dir2, pv);
      c_sorted_j1.emplace_back(i_cand, track_info_j1.getTrackSip2dSig(), -btagbtvdeep::mindrsvpfcand(svs_unsorted, cand), cand->pt() / jet1.pt());
      c_sorted_j2.emplace_back(i_cand, track_info_j2.getTrackSip2dSig(), -btagbtvdeep::mindrsvpfcand(svs_unsorted, cand), cand->pt() / jet2.pt());
    } else {
      // add neutral candidates
      n_sorted_j1.emplace_back(i_cand, -1, -btagbtvdeep::mindrsvpfcand(svs_unsorted, cand), cand->pt() / jet1.pt());
      n_sorted_j2.emplace_back(i_cand, -1, -btagbtvdeep::mindrsvpfcand(svs_unsorted, cand), cand->pt() / jet2.pt());
    }
  }
  // sort by relation to jet 1
  std::sort(c_sorted_j1.begin(), c_sorted_j1.end(), btagbtvdeep::SortingClass<std::size_t>::compareByABCInv);
  std::sort(n_sorted_j1.begin(), n_sorted_j1.end(), btagbtvdeep::SortingClass<std::size_t>::compareByABCInv);
  std::vector<size_t> c_sortedidx, n_sortedidx;
  // this puts 0 everywhere and the right position in ind
  c_sortedidx = btagbtvdeep::invertSortingVector(c_sorted_j1);
  n_sortedidx = btagbtvdeep::invertSortingVector(n_sorted_j1);
  // loop through charged pf candidates
  /*
  for (size_t i_cand = 0; i_cand < c_sortedidx.size(); ++i_cand) {
    auto entry = c_sortedidx.at(i_cand);
    const reco::Candidate *cand = &(cands->at(entry));
    auto packed_cand = dynamic_cast<const pat::PackedCandidate*>(cand);
    auto cand_ptr = dynamic_cast<reco::CandidatePtr>(cand);
    float puppiw = (*puppi_value_map)[cand_ptr];
    float drminpfcandsv = btagbtvdeep::mindrsvpfcand(svs_unsorted, cand);
    auto& track_info_j1 = trackinfos_j1.at(entry);
    auto& track_info_j2 = trackinfos_j2.at(entry);
    const reco::Track& PseudoTrack = packed_cand->pseudoTrack();
    reco::TransientTrack transientTrack;
    transientTrack = track_builder->build(PseudoTrack);
    // fill charged candidate variables
    ptr = &data_[kChargedCandidates][offset + i_cand * n_features_cpf_];
    start = ptr;
    *ptr = catch_infs_and_bound(track_info_j1.getTrackEtaRel(), 0, -5, 15); // track etarel 1
    *(++ptr) = catch_infs_and_bound(track_info_j2.getTrackEtaRel(), 0, -5, 15); // track etarel 2
    *(++ptr) = catch_infs_and_bound(track_info_j1.getTrackPtRel(), 0, -1, 4); // track ptrel 1
    *(++ptr) = catch_infs_and_bound(track_info_j2.getTrackPtRel(), 0, -1, 4); // track ptrel 2
    *(++ptr) = catch_infs_and_bound(track_info_j1.getTrackPPar(), 0, -1e5, 1e5); // track ppar 1
    *(++ptr) = catch_infs_and_bound(track_info_j2.getTrackPPar(), 0, -1e5, 1e5); // track ppar 2
    *(++ptr) = catch_infs_and_bound(track_info_j1.getTrackDeltaR(), 0, -5, 5); // track deltaR 1
    *(++ptr) = catch_infs_and_bound(track_info_j2.getTrackDeltaR(), 0, -5, 5); // track deltaR 2
    *(++ptr) = catch_infs_and_bound(track_info_j1.getTrackPParRatio(), 0, -10, 100); // track ppar ratio 1
    *(++ptr) = catch_infs_and_bound(track_info_j2.getTrackPParRatio(), 0, -10, 100); // track ppar ratio 2
    *(++ptr) = catch_infs_and_bound(track_info_j1.getTrackSip2dVal(), 0, -1, 70); // track sip2d val 1
    *(++ptr) = catch_infs_and_bound(track_info_j2.getTrackSip2dVal(), 0, -1, 70); // track sip2d val 2
    *(++ptr) = catch_infs_and_bound(track_info_j1.getTrackSip2dSig(), 0, -1, 4e4); // track sip2d sig 1
    *(++ptr) = catch_infs_and_bound(track_info_j2.getTrackSip2dSig(), 0, -1, 4e4); // track sip2d sig 2
    *(++ptr) = catch_infs_and_bound(track_info_j1.getTrackSip3dVal(), 0, -1, 1e5); // track sip3d val 1
    *(++ptr) = catch_infs_and_bound(track_info_j2.getTrackSip3dVal(), 0, -1, 1e5); // track sip3d val 2
    *(++ptr) = catch_infs_and_bound(track_info_j1.getTrackSip3dSig(), 0, -1, 4e4); // track sip3d sig 1
    *(++ptr) = catch_infs_and_bound(track_info_j2.getTrackSip3dSig(), 0, -1, 4e4); // track sip3d sig 2
    *(++ptr) = catch_infs_and_bound(track_info_j1.getTrackJetDistVal(), 0, -20, 1); // track jet dist val 1
    *(++ptr) = catch_infs_and_bound(track_info_j2.getTrackJetDistVal(), 0, -20, 1); // track jet dist val 2
    *(++ptr) = catch_infs_and_bound((packed_cand->pt() * puppiw) / jet1.pt(), 0, -1, 0, -1); // ptrel 1
    *(++ptr) = catch_infs_and_bound((packed_cand->pt() * puppiw) / jet2.pt(), 0, -1, 0, -1); // ptrel 2
    *(++ptr) = catch_infs_and_bound(drminpfcandsv, 0, -1. * jet_radius_, 0, -1. * jet_radius_); // drminsv
    *(++ptr) = packed_cand->pvAssociationQuality(); // vtx_ass
    *(++ptr) = puppiw; // puppiw
    *(++ptr) = catch_infs_and_bound(PseudoTrack.normalizedChi2(), 300, -1, 300); // chi2
    *(++ptr) = PseudoTrack.qualityMask(); // quality
    assert(start + n_features_cpf_ - 1 == ptr);
    // fill charged candidate 4-vectors
    ptr = &data_[kChargedCandidates4Vec][offset + i_cand * n_pairwise_features_cpf_];
    start = ptr;
    *ptr = packed_cand->px();
    *(++ptr) = packed_cand->py();
    *(++ptr) = packed_cand->pz();
    *(++ptr) = packed_cand->energy();
    assert(start + n_pairwise_features_cpf_ - 1 == ptr);
  }*/
  // loop through neutral candidates
  for (size_t i_cand = 0; i_cand < n_sortedidx.size(); ++i_cand) {
    auto entry = n_sortedidx.at(i_cand);
    const reco::Candidate *cand = &(cands->at(i_cand));
    auto packed_cand = dynamic_cast<const pat::PackedCandidate*>(cand);
    //auto cand_ptr = dynamic_cast<reco::CandidatePtr>(cand);
    auto puppiw = packed_cand->puppiWeight();//(*puppi_value_map)[cand_ptr];
    float drminpfcandsv = btagbtvdeep::mindrsvpfcand(svs_unsorted, cand);
    // fill neutral candidate variables
    ptr = &data_[0][offset + i_cand * 4];
    *ptr = catch_infs_and_bound(std::fabs(packed_cand->eta() - jet1.eta()), 0, -2, 0, -0.5); // etarel 1
    *(++ptr) = catch_infs_and_bound(std::fabs(reco::deltaPhi(packed_cand->phi(), jet1.phi())), 0, -2, 0, -0.5); // phirel 1
    *(++ptr) = catch_infs_and_bound(std::fabs(packed_cand->eta() - jet2.eta()), 0, -2, 0, -0.5); // etarel 2
    *(++ptr) = catch_infs_and_bound(std::fabs(reco::deltaPhi(packed_cand->phi(), jet2.phi())), 0, -2, 0, -0.5); // phirel 2
    ptr = &data_[1][offset + i_cand * 16];
    *ptr = catch_infs_and_bound((packed_cand->pt() * puppiw) / jet1.pt(), 0, -1, 0, -1); // ptrel 1
    *(++ptr) = catch_infs_and_bound((packed_cand->pt() * puppiw) / jet2.pt(), 0, -1, 0, -1); // ptrel 2
    *(++ptr) = catch_infs_and_bound((packed_cand->pt() * puppiw) / jet1.pt(), 0, -1, 0, -1); // ptrel 1
    *(++ptr) = catch_infs_and_bound((packed_cand->pt() * puppiw) / jet2.pt(), 0, -1, 0, -1); // ptrel 2
    *(++ptr) = catch_infs_and_bound(reco::deltaR(*packed_cand, jet1), 0, -0.6, 0, -0.6); // deltaR 1
    *(++ptr) = catch_infs_and_bound(reco::deltaR(*packed_cand, jet2), 0, -0.6, 0, -0.6); // deltaR 2
    *(++ptr) = 0;
    *(++ptr) = 0;
    *(++ptr) = 0;
    *(++ptr) = 0;
    *(++ptr) = 0;
    *(++ptr) = catch_infs_and_bound(std::fabs(packed_cand->eta() - jet1.eta()), 0, -2, 0, -0.5); // etarel 1
    *(++ptr) = catch_infs_and_bound(std::fabs(reco::deltaPhi(packed_cand->phi(), jet1.phi())), 0, -2, 0, -0.5); // phirel 1
    *(++ptr) = catch_infs_and_bound(std::fabs(packed_cand->eta() - jet2.eta()), 0, -2, 0, -0.5); // etarel 2
    *(++ptr) = catch_infs_and_bound(std::fabs(reco::deltaPhi(packed_cand->phi(), jet2.phi())), 0, -2, 0, -0.5); // phirel 2
    *(++ptr) = 0;
    ptr = &data_[2][offset + i_cand * 4];
    *ptr = packed_cand->px();
    *(++ptr) = packed_cand->py();
    *(++ptr) = packed_cand->pz();
    *(++ptr) = packed_cand->energy();
    ptr = &data_[3][offset + i_cand * 1];
    *ptr = 0;
    // TEST STUFF
    if (fill_test) {
      pf_cand_es.emplace_back(packed_cand->energy());
      pf_cand_pxs.emplace_back(packed_cand->px());
      pf_cand_pys.emplace_back(packed_cand->py());
      pf_cand_pzs.emplace_back(packed_cand->pz());
    }
    // END OF TEST STUFF
    /*
    ptr = &data_[0][offset + i_cand * n_features_npf_];
    start = ptr;
    *ptr = catch_infs_and_bound((packed_cand->pt() * puppiw) / jet1.pt(), 0, -1, 0, -1); // ptrel 1
    *(++ptr) = catch_infs_and_bound((packed_cand->pt() * puppiw) / jet2.pt(), 0, -1, 0, -1); // ptrel 2
    *(++ptr) = catch_infs_and_bound(std::fabs(packed_cand->eta() - jet1.eta()), 0, -2, 0, -0.5); // etarel 1
    *(++ptr) = catch_infs_and_bound(std::fabs(packed_cand->eta() - jet2.eta()), 0, -2, 0, -0.5); // etarel 2
    *(++ptr) = catch_infs_and_bound(std::fabs(reco::deltaPhi(packed_cand->phi(), jet1.phi())), 0, -2, 0, -0.5); // phirel 1
    *(++ptr) = catch_infs_and_bound(std::fabs(reco::deltaPhi(packed_cand->phi(), jet2.phi())), 0, -2, 0, -0.5); // phirel 2
    *(++ptr) = catch_infs_and_bound(reco::deltaR(*packed_cand, jet1), 0, -0.6, 0, -0.6); // deltaR 1
    *(++ptr) = catch_infs_and_bound(reco::deltaR(*packed_cand, jet2), 0, -0.6, 0, -0.6); // deltaR 2
    if (std::abs(packed_cand->pdgId()) == 22) { // isGamma
      *(++ptr) = 1;
    } else {
      *(++ptr) = 0;
    }
    *(++ptr) = packed_cand->hcalFraction(); // hadFrac
    *(++ptr) = catch_infs_and_bound(drminpfcandsv, 0, -1. * jet_radius_, 0, -1. * jet_radius_); // drminsv
    *(++ptr) = puppiw; // puppiw
    assert(start + n_features_npf_ - 1 == ptr);
    // fill neutral candidate 4-vectors
    ptr = &data_[kNeutralCandidates4Vec][offset + i_cand * n_pairwise_features_npf_];
    start = ptr;
    *ptr = packed_cand->px();
    *(++ptr) = packed_cand->py();
    *(++ptr) = packed_cand->pz();
    *(++ptr) = packed_cand->energy();
    assert(start + n_pairwise_features_npf_ - 1 == ptr);*/
  }
}

//define this as a plug-in

//typedef PAIReDONNXProducer<1> PAIReDTableProducer;
DEFINE_FWK_MODULE(PAIReDONNXProducer);