#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/Framework/interface/makeRefToBaseProdFrom.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/StreamID.h"

#include "DataFormats/BTauReco/interface/JetTag.h"

#include "DataFormats/BTauReco/interface/PAIReDTagInfo.h"
#include "DataFormats/BTauReco/interface/PAIReDFeatures.h"

#include "PhysicsTools/ONNXRuntime/interface/ONNXRuntime.h"

using namespace cms::Ort;

class PAIReDONNXJetTagsProducer : public edm::stream::EDProducer<edm::GlobalCache<ONNXRuntime>> {
public:
  explicit PAIReDONNXJetTagsProducer(const edm::ParameterSet&, const ONNXRuntime*);
  ~PAIReDONNXJetTagsProducer() override = default;

  static void fillDescriptions(edm::ConfigurationDescriptions&);

  static std::unique_ptr<ONNXRuntime> initializeGlobalCache(const edm::ParameterSet&);
  static void globalEndJob(const ONNXRuntime*);

private:
  typedef std::vector<reco::PAIReDTagInfo> TagInfoCollection;
  typedef reco::JetTagCollection JetTagCollection;

  void produce(edm::Event&, const edm::EventSetup&) override;

  void make_inputs(btagbtvdeep::PAIReDFeatures features);
  void get_input_sizes(const reco::FeaturesTagInfo<btagbtvdeep::PAIReDFeatures> taginfo);

  const edm::EDGetTokenT<TagInfoCollection> src_;
  std::vector<std::string> flav_names_;
  std::vector<std::string> input_names_;
  std::vector<std::string> output_names_;

  enum InputIndexes {
    kChargedCandidates_j1 = 0,
    kChargedCandidates_j2 = 1,
    kNeutralCandidates_j1 = 2,
    kNeutralCandidates_j2 = 3,
    kVertices = 4,
    kChargedCandidates4Vec = 5,
    kNeutralCandidates4Vec = 6,
    kVertices4Vec = 7
  };
  unsigned n_cpf_;
  constexpr static unsigned n_features_cpf_ = 16;
  constexpr static unsigned n_pairwise_features_cpf_ = 4;
  unsigned n_npf_;
  constexpr static unsigned n_features_npf_ = 8;
  constexpr static unsigned n_pairwise_features_npf_ = 4;
  unsigned n_sv_;
  constexpr static unsigned n_features_sv_ = 14;
  constexpr static unsigned n_pairwise_features_sv_ = 4;
  std::vector<unsigned> input_sizes_;
  std::vector<std::vector<int64_t>> input_shapes_;  // shapes of each input group (-1 for dynamic axis)

  // hold the input data
  FloatArrays data_;
};

PAIReDJetTagsProducer::PAIReDONNXJetTagsProducer(const edm::ParameterSet& iConfig,
                                                                                     const ONNXRuntime* cache)
    : src_(consumes<TagInfoCollection>(iConfig.getParameter<edm::InputTag>("src"))),
      flav_names_(iConfig.getParameter<std::vector<std::string>>("flav_names")),
      input_names_(iConfig.getParameter<std::vector<std::string>>("input_names")),
      output_names_(iConfig.getParameter<std::vector<std::string>>("output_names")) {
  // get output names from flav_names
  for (const auto& flav_name : flav_names_) {
    produces<JetTagCollection>(flav_name);
  }
}

void PAIReDONNXJetTagsProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  // pfParticleTransformerAK4JetTags
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("src", edm::InputTag("pfPAIReDTagInfos"));
  desc.add<std::vector<std::string>>("input_names", {"input_1", "input_2", "input_3", "input_4", "input_5", "input_6", "input_7", "input_8"});
  desc.add<edm::FileInPath>("model_path",
                            // DONT FORGET TO EDIT ME |
                            // DONT FORGET TO EDIT ME |
                            // DONT FORGET TO EDIT ME \/
                            edm::FileInPath("PATHTOPAIREDTAGGERMODEL.onnx"));
                          // DONT FORGET TO EDIT ME /\
                          // DONT FORGET TO EDIT ME |
                          // DONT FORGET TO EDIT ME |
  desc.add<std::vector<std::string>>("output_names", {"softmax"});
  desc.add<std::vector<std::string>>(
      "flav_names", std::vector<std::string>{"probbb", "probbcc", "probudsg"});

  descriptions.add("pfPAIReDJetTags", desc);
}

std::unique_ptr<ONNXRuntime> PAIReDONNXJetTagsProducer::initializeGlobalCache(
    const edm::ParameterSet& iConfig) {
  return std::make_unique<ONNXRuntime>(iConfig.getParameter<edm::FileInPath>("model_path").fullPath());
}

void PAIReDONNXJetTagsProducer::globalEndJob(const ONNXRuntime* cache) {}

void PAIReDONNXJetTagsProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
  edm::Handle<> tag_infos;
  iEvent.getByToken(src_, tag_infos);

  // initialize output collection
  std::vector<std::unique_ptr<JetTagCollection>> output_tags;
  if (!tag_infos->empty()) {
    auto jet_ref = tag_infos->begin()->jet();
    auto ref2prod = edm::makeRefToBaseProdFrom(jet_ref, iEvent);
    for (std::size_t i = 0; i < flav_names_.size(); i++) {
      output_tags.emplace_back(std::make_unique<JetTagCollection>(ref2prod));
    }
  } else {
    for (std::size_t i = 0; i < flav_names_.size(); i++) {
      output_tags.emplace_back(std::make_unique<JetTagCollection>());
    }
  }

  for (unsigned jet_n = 0; jet_n < tag_infos->size(); ++jet_n) {
    const auto& taginfo = (*tag_infos)[jet_n];
    std::vector<float> outputs(flav_names_.size(), -1.0);
    if (taginfo.features().is_filled) {
      get_input_sizes(taginfo);

      // run prediction with dynamic batch size per event
      input_shapes_ = {{(int64_t)1, (int64_t)n_cpf_, (int64_t)n_features_cpf_},
                       {(int64_t)1, (int64_t)n_cpf_, (int64_t)n_features_cpf_},
                       {(int64_t)1, (int64_t)n_npf_, (int64_t)n_features_npf_},
                       {(int64_t)1, (int64_t)n_npf_, (int64_t)n_features_npf_},
                       {(int64_t)1, (int64_t)n_sv_, (int64_t)n_features_sv_},
                       {(int64_t)1, (int64_t)n_cpf_, (int64_t)n_pairwise_features_cpf_},
                       {(int64_t)1, (int64_t)n_npf_, (int64_t)n_pairwise_features_npf_},
                       {(int64_t)1, (int64_t)n_sv_, (int64_t)n_pairwise_features_sv_}};

      outputs = globalCache()->run(input_names_, data_, input_shapes_, output_names_, 1)[0];
      assert(outputs.size() == flav_names_.size());
    }

    const auto& jet_ref = tag_infos->at(jet_n).jet();
    for (std::size_t flav_n = 0; flav_n < flav_names_.size(); flav_n++) {
      (*(output_tags[flav_n]))[jet_ref] = outputs[flav_n];
    }
  }

  // put into the event
  for (std::size_t flav_n = 0; flav_n < flav_names_.size(); ++flav_n) {
    iEvent.put(std::move(output_tags[flav_n]), flav_names_[flav_n]);
  }
}

void PAIReDONNXJetTagsProducer::get_input_sizes(
    const reco::PAIReDFeaturesTagInfo<btagbtvdeep::PAIReDFeatures> taginfo) {
  const auto& features = taginfo.features();

  unsigned int n_cpf = features.c_pf_features_j1.size();
  unsigned int n_npf = features.n_pf_features_j1.size();
  unsigned int n_vtx = features.sv_features.size();

  n_cpf_ = std::max((unsigned int)1, n_cpf);
  n_npf_ = std::max((unsigned int)1, n_npf);
  n_sv_ = std::max((unsigned int)1, n_vtx);

  n_cpf_ = std::min((unsigned int)25, n_cpf_);
  n_npf_ = std::min((unsigned int)25, n_npf_);
  n_sv_ = std::min((unsigned int)5, n_sv_);
  input_sizes_ = {
      n_cpf_ * n_features_cpf_,
      n_cpf_ * n_features_cpf_,
      n_npf_ * n_features_npf_,
      n_npf_ * n_features_npf_,
      n_sv_ * n_features_sv_,
      n_cpf_ * n_pairwise_features_cpf_,
      n_npf_ * n_pairwise_features_npf_,
      n_sv_ * n_pairwise_features_sv_,
  };
  // init data storage
  data_.clear();
  for (const auto& len : input_sizes_) {
    data_.emplace_back(1 * len, 0);
  }

  make_inputs(features);
}

void PAIReDONNXJetTagsProducer::make_inputs(btagbtvdeep::PAIReDFeatures features) {
  float* ptr = nullptr;
  const float* start = nullptr;
  unsigned offset = 0;

  // c_pf candidates for jet 1
  auto max_c_pf_n = std::min(features.c_pf_features_j1.size(), (std::size_t)n_cpf_);
  for (std::size_t c_pf_n = 0; c_pf_n < max_c_pf_n; c_pf_n++) {
    const auto& c_pf_features_j1 = features.c_pf_features_j1.at(c_pf_n);
    ptr = &data_[kChargedCandidates_j1][offset + c_pf_n * n_features_cpf_];
    start = ptr;
    *ptr = c_pf_features_j1.btagPf_trackEtaRel;
    *(++ptr) = c_pf_features_j1.btagPf_trackPtRel;
    *(++ptr) = c_pf_features_j2.btagPf_trackPtRel;
    *(++ptr) = c_pf_features_j1.btagPf_trackPPar;
    *(++ptr) = c_pf_features_j2.btagPf_trackPPar;
    *(++ptr) = c_pf_features_j1.btagPf_trackDeltaR;
    *(++ptr) = c_pf_features_j2.btagPf_trackDeltaR;
    *(++ptr) = c_pf_features_j1.btagPf_trackPParRatio;
    *(++ptr) = c_pf_features_j2.btagPf_trackPParRatio;
    *(++ptr) = c_pf_features_j1.btagPf_trackSip2dVal;
    *(++ptr) = c_pf_features_j2.btagPf_trackSip2dVal;
    *(++ptr) = c_pf_features_j1.btagPf_trackSip2dSig;
    *(++ptr) = c_pf_features_j2.btagPf_trackSip2dSig;
    *(++ptr) = c_pf_features_j1.btagPf_trackSip3dVal;
    *(++ptr) = c_pf_features_j2.btagPf_trackSip3dVal;
    *(++ptr) = c_pf_features_j1.btagPf_trackSip3dSig;
    *(++ptr) = c_pf_features_j2.btagPf_trackSip3dSig;
    *(++ptr) = c_pf_features_j1.btagPf_trackJetDistVal;
    *(++ptr) = c_pf_features_j2.btagPf_trackJetDistVal;
    *(++ptr) = c_pf_features_j1.ptrel;
    *(++ptr) = c_pf_features_j2.ptrel;
    *(++ptr) = c_pf_features_j1.drminsv;
    *(++ptr) = c_pf_features_j1.vtx_ass;
    *(++ptr) = c_pf_features_j1.puppiw;
    *(++ptr) = c_pf_features_j1.chi2;
    *(++ptr) = c_pf_features_j1.quality;
    assert(start + n_features_cpf_ - 1 == ptr);
  }

  // n_pf candidates
  auto max_n_pf_n = std::min(features.n_pf_features_j1.size(), (std::size_t)n_npf_);
  for (std::size_t n_pf_n = 0; n_pf_n < max_n_pf_n; n_pf_n++) {
    const auto& n_pf_features_j1 = features.n_pf_features_j1.at(n_pf_n);
    ptr = &data_[kNeutralCandidates_j1][offset + n_pf_n * n_features_npf_];
    start = ptr;
    *ptr = n_pf_features_j1.ptrel;
    *(++ptr) = n_pf_features_j2.ptrel;
    *(++ptr) = n_pf_features_j1.etarel;
    *(++ptr) = n_pf_features_j2.etarel;
    *(++ptr) = n_pf_features_j1.phirel;
    *(++ptr) = n_pf_features_j1.deltaR;
    *(++ptr) = n_pf_features_j2.deltaR;
    *(++ptr) = n_pf_features_j1.isGamma; // CHECK
    *(++ptr) = n_pf_features_j1.hadFrac; // CHECK
    *(++ptr) = n_pf_features_j1.drminsv;
    *(++ptr) = n_pf_features_j1.puppiw;
    assert(start + n_features_npf_ - 1 == ptr);
  }

  // sv candidates
  auto max_sv_n = std::min(features.sv_features.size(), (std::size_t)n_sv_);
  for (std::size_t sv_n = 0; sv_n < max_sv_n; sv_n++) {
    const auto& sv_features = features.sv_features.at(sv_n);
    ptr = &data_[kVertices][offset + sv_n * n_features_sv_];
    start = ptr;
    *ptr = sv_features.pt;
    *(++ptr) = sv_features.deltaR;
    *(++ptr) = sv_features.mass;
    *(++ptr) = sv_features.etarel;
    *(++ptr) = sv_features.phirel;
    *(++ptr) = sv_features.ntracks;
    *(++ptr) = sv_features.chi2;
    *(++ptr) = sv_features.normchi2;
    *(++ptr) = sv_features.dxy;
    *(++ptr) = sv_features.dxysig;
    *(++ptr) = sv_features.d3d;
    *(++ptr) = sv_features.d3dsig;
    *(++ptr) = sv_features.costhetasvpv;
    *(++ptr) = sv_features.enratio;
    assert(start + n_features_sv_ - 1 == ptr);
  }

  // cpf pairwise features (4-vectors)
  auto max_cpf_n = std::min(features.c_pf_features_j1.size(), (std::size_t)n_cpf_);
  for (std::size_t cpf_n = 0; cpf_n < max_cpf_n; cpf_n++) {
    const auto& cpf_pairwise_features = features.c_pf_features_j1.at(cpf_n);
    ptr = &data_[kChargedCandidates4Vec][offset + cpf_n * n_pairwise_features_cpf_];
    start = ptr;
    *ptr = cpf_pairwise_features.px;
    *(++ptr) = cpf_pairwise_features.py;
    *(++ptr) = cpf_pairwise_features.pz;
    *(++ptr) = cpf_pairwise_features.e;

    assert(start + n_pairwise_features_cpf_ - 1 == ptr);
  }

  // npf pairwise features (4-vectors)
  auto max_npf_n = std::min(features.n_pf_features_j1.size(), (std::size_t)n_npf_);
  for (std::size_t npf_n = 0; npf_n < max_npf_n; npf_n++) {
    const auto& npf_pairwise_features = features.n_pf_features_j1.at(npf_n);
    ptr = &data_[kNeutralCandidates4Vec][offset + npf_n * n_pairwise_features_npf_];
    start = ptr;
    *ptr = npf_pairwise_features.px;
    *(++ptr) = npf_pairwise_features.py;
    *(++ptr) = npf_pairwise_features.pz;
    *(++ptr) = npf_pairwise_features.e;

    assert(start + n_pairwise_features_npf_ - 1 == ptr);
  }

  // sv pairwise features (4-vectors)
  auto max_sv_N = std::min(features.sv_features.size(), (std::size_t)n_sv_);
  for (std::size_t sv_N = 0; sv_N < max_sv_N; sv_N++) {
    const auto& sv_pairwise_features = features.sv_features.at(sv_N);
    ptr = &data_[kVertices4Vec][offset + sv_N * n_pairwise_features_sv_];
    start = ptr;
    *ptr = sv_pairwise_features.px;
    *(++ptr) = sv_pairwise_features.py;
    *(++ptr) = sv_pairwise_features.pz;
    *(++ptr) = sv_pairwise_features.e;

    assert(start + n_pairwise_features_sv_ - 1 == ptr);
  }
}

//define this as a plug-in
DEFINE_FWK_MODULE(PAIReDONNXJetTagsProducer);
