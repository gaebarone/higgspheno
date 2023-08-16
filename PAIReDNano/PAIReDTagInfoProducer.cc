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

#include "DataFormats/BTauReco/interface/PAIReDTagInfo.h"
#include "DataFormats/BTauReco/interface/PAIReDFeatures.h"

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

#include "inEllipse.h"

template<typename T>
class PAIReDTagInfoProducer : public edm::stream::EDProducer<> {
public:
  explicit PAIReDTagInfoProducer(const edm::ParameterSet &);
  ~PAIReDTagInfoProducer() override;

  static void fillDescriptions(edm::ConfigurationDescriptions &descriptions);

private:
  typedef std::vector<reco::PAIReDTagInfo> PAIReDTagInfoCollection;
  typedef reco::VertexCompositePtrCandidateCollection SVCollection;
  typedef reco::VertexCollection VertexCollection;
  void produce(edm::Event &, const edm::EventSetup &) override;

  //const std::string name_;
  const std::string name_;

  const edm::EDGetTokenT<edm::View<T>> jet_token_;
  const edm::EDGetTokenT<VertexCollection> vtx_token_;
  const edm::EDGetTokenT<SVCollection> sv_token_;
  const edm::EDGetTokenT<reco::CandidateView> cand_token_;
  
};

//
// constructors and destructor
//
template< typename T>
PAIReDTagInfoProducer<T>::PAIReDTagInfoProducer(const edm::ParameterSet &iConfig)
    : jet_token_(consumes<edm::View<T>>(iConfig.getParameter<edm::InputTag>("jets"))),
      cand_token_(consumes<reco::CandidateView>(iConfig.getParameter<edm::InputTag>("candidates"))),
      vtx_token_(consumes<VertexCollection>(iConfig.getParameter<edm::InputTag>("vertices"))),
      sv_token_(consumes<SVCollection>(iConfig.getParameter<edm::InputTag>("secondary_vertices"))),
      track_builder_token_(
          esConsumes<TransientTrackBuilder, TransientTrackRecord>(edm::ESInputTag("", "TransientTrackBuilder"))),
      puppi_value_map_token_(consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("puppi_value_map"))){
  produces<PAIReDTagInfoCollection>;
}

template< typename T>
PAIReDTagInfoProducer<T>::~PAIReDTagInfoProducer() {}

template< typename T>
void PAIReDTagInfoProducer<T>::fillDescriptions(edm::ConfigurationDescriptions &descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("jets", edm::InputTag("slimmedJetsAK4"));
  desc.add<edm::InputTag>("vertices", edm::InputTag("offlinePrimaryVertices"));
  desc.add<edm::InputTag>("secondary_vertices", edm::InputTag("inclusiveCandidateSecondaryVertices"));
  desc.add<edm::InputTag>("candidates", edm::InputTag("packedPFCandidates"));
  descriptions.addWithDefaultLabel(desc);
}

template< typename T>
void PAIReDTagInfoProducer<T>::produce(edm::Event &iEvent, const edm::EventSetup &iSetup) {
  auto output_tag_infos = std::make_unique<PAIReDTagInfoCollection>();
  // produce empty TagInfos in case no primary vertex
  edm::Handle<VertexCollection> vtxs;
  iEvent.getByToken(vtx_token_, vtxs);
  if (vtxs->empty()) {
    iEvent.put(std::move(output_tag_infos));
    return;
  }
  // get jets, pf candidates, secondary vertices, puppi weights, track builder
  edm::Handle<edm::View<reco::Jet>> jets;
  iEvent.getByToken(jet_token_, jets);
  double jet_radius_ = 0.4;
  bool is_weighted_jet_ = true;
  bool flip_ = false;
  edm::Handle<edm::View<reco::Candidate>> cands;
  iEvent.getByToken(candidateToken_, cands);
  edm::Handle<SVCollection> svs;
  iEvent.getByToken(sv_token_, svs);
  edm::Handle<edm::ValueMap<float>> puppi_value_map;
  iEvent.getByToken(puppi_value_map_token_, puppi_value_map);
  edm::ESHandle<TransientTrackBuilder> track_builder = iSetup.getHandle(track_builder_token_);

  //loop through primary jets
  for (unsigned i_jet = 0; i_jet < jets->size(); ++i_jet) {
    const auto &jet1 = jets->at(i_jet);
    edm::RefToBase<reco::Jet> jet_ref1(jets, i_jet);
    const auto* pat_jet1 = dynamic_cast<const pat::Jet*>(&jet1);
    math::XYZVector jet_dir1 = jet1.momentum().Unit();
    GlobalVector jet_ref_track_dir1(jet1.px(), jet1.py(), jet1.pz());
    puppi_weight1 = (*puppi_value_map)[jet1->getPFConstituent(i_jet)];
    //loop through secondary jets
    for (unsigned j_jet = i_jet+1; j_jet < jets->size(); ++j_jet) {
      const auto &jet2 = jets->at(j_jet);
      edm::RefToBase<reco::Jet> jet_ref2(jets, j_jet);
      const auto* pat_jet2 = dynamic_cast<const pat::Jet*>(&jet2);
      math::XYZVector jet_dir2 = jet2.momentum().Unit();
      GlobalVector jet_ref_track_dir2(jet2.px(), jet2.py(), jet2.pz());
      puppi_weight2 = (*puppi_value_map)[jet2->getPFConstituent(j_jet)];
      btagbtvdeep::PAIReDFeatures features;
      // sort secondary vertices by dxy
      auto svs_sorted = *svs;
      std::sort(svs_sorted.begin(), svs_sorted.end(), [&pv](const auto& sva, const auto& svb) {
        return btagbtvdeep::sv_vertex_comparator(sva, svb, pv);
      });
      // add in secondary vertex information (wrt both jets) if in the ellipse
      for (const auto& sv : svs_sorted) {
        if (!inEllipse(jet1, jet2, sv)) continue;
        features.sv_features_j1.emplace_back();
        features.sv_features_j2.emplace_back();
        auto& sv_features_j1 = features.sv_features_j1.back();
        auto& sv_features_j2 = features.sv_features_j2.back();
        btagbtvdeep::svToFeatures(sv, pv, jet1, sv_features_j1);
        btagbtvdeep::svToFeatures(sv, pv, jet2, sv_features_j2);
      }
      const auto& svs_unsorted = *svs;
      // go through pf candidates and save them to add in order
      std::vector<btagbtvdeep::SortingClass<size_t>> c_sorted_j1, n_sorted_j1, c_sorted_j2, n_sorted_j2;
      std::map<unsigned int, btagbtvdeep::TrackInfoBuilder> trackinfos_j1, trackinfos_j2;
      for (unsigned i_cand = 0; i_cand < cands->size(); ++i_cand) {
        const auto &cand = cands->at(i_cand);
        if !inEllipse(jet1, jet2, cand) continue;
        if (cand->charge() != 0) {
          // build track info and keep charged candidates
          auto& trackinfo_j1 = trackinfos_j1.emplace(i_cand, track_builder).first->second;
          auto& trackinfo_j2 = trackinfos_j2.emplace(i_cand, track_builder).first->second;
          trackinfo_j1.buildTrackInfo(cand, jet_dir1, jet_ref_track_dir1, pv);
          trackinfo_j2.buildTrackInfo(cand, jet_dir2, jet_ref_track_dir2, pv);
          c_sorted_j1.emplace_back(i_cand, trackinfo_j1.getTrackSip2dSig(), -btagbtvdeep::mindrsvpfcand(svs_unsorted, cand), cand->pt() / jet1.pt());
          c_sorted_j2.emplace_back(i_cand, trackinfo_j2.getTrackSip2dSig(), -btagbtvdeep::mindrsvpfcand(svs_unsorted, cand), cand->pt() / jet2.pt());
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
      // set right size to vectors
      features.c_pf_features_j1.clear();
      features.c_pf_features_j1.resize(c_sorted.size());
      features.c_pf_features_j2.clear();
      features.c_pf_features_j2.resize(c_sorted.size());
      features.n_pf_features_j1.clear();
      features.n_pf_features_j1.resize(n_sorted.size());
      features.n_pf_features_j2.clear();
      features.n_pf_features_j2.resize(n_sorted.size());
      // now go back through and place candidates in order
      for (unsigned i_cand = 0; i_cand < cands->size(); ++i_cand) {
        const auto &cand = cands->at(i_cand);
        // check if in ellipse and get relevant info
        if !inEllipse(jet1, jet2, cand) continue;
        auto packed_cand = dynamic_cast<const pat::PackedCandidate*>(cand);
        auto cand_ptr = dynamic_cast<reco::CandidatePtr>cand
        float puppiw = (*puppi_value_map)[cand_ptr];
        float drminpfcandsv = btagbtvdeep::mindrsvpfcand(svs_unsorted, cand);
        float distminpfcandsv = 0;
        // deal with charged candidates
        if (cand->charge() != 0) {
          auto entry = c_sortedidx.at(i_cand);
          auto& c_pf_features_j1 = features.c_pf_features_j1.at(entry);
          auto& c_pf_features_j2 = features.c_pf_features_j2.at(entry);
          //build track info
          auto& trackinfo_j1 = trackinfos_j1.at(i_cand);
          auto& trackinfo_j2 = trackinfos_j2.at(i_cand);
          const reco::Track& PseudoTrack = packed_cand->pseudoTrack();
          reco::TransientTrack transientTrack;
          transientTrack = track_builder->build(PseudoTrack);
          distminpfcandsv = btagbtvdeep::mindistsvpfcand(svs_unsorted, transientTrack);
          btagbtvdeep::packedCandidateToFeatures(packed_cand, jet1, trackinfo_j1, is_weighted_jet_, drminpfcandsv, static_cast<float>(jet_radius_), puppiw, c_pf_features_j1, flip_, distminpfcandsv);
          btagbtvdeep::packedCandidateToFeatures(packed_cand, jet2, trackinfo_j2, is_weighted_jet_, drminpfcandsv, static_cast<float>(jet_radius_), puppiw, c_pf_features_j2, flip_, distminpfcandsv);
        } else { // deal with neutral candidates
          auto entry = n_sortedidx.at(i);
          auto& n_pf_features_j1 = features.n_pf_features_j1.at(entry);
          auto& n_pf_features_j2 = features.n_pf_features_j2.at(entry);
          btagbtvdeep::packedCandidateToFeatures(packed_cand, jet1, is_weighted_jet_, drminpfcandsv, static_cast<float>(jet_radius_), puppiw, n_pf_features_j1);
          btagbtvdeep::packedCandidateToFeatures(packed_cand, jet2, is_weighted_jet_, drminpfcandsv, static_cast<float>(jet_radius_), puppiw, n_pf_features_j2);
        }
      }
      // add features to output
      output_tag_infos->emplace_back(features, jet_ref1, jet_ref2);
    }
  }
  iEvent.put(std::move(output_tag_infos));
}

DEFINE_FWK_MODULE(PAIReDTagInfoProducer);