#ifndef HADRONIC_H
#define HADRONIC_H

#include "../includes/kinematics_include.h"
#include "../includes/helperfunctions_include.h"
#include "../../common_includes/ghost_tagging.h"
#include "../../common_includes/make_paired.h"


//------------------------------------------------------------------------------------------------------------------------------------------------------------
// pseudo b tagging
//------------------------------------------------------------------------------------------------------------------------------------------------------------


bool ghost_btagging(TClonesArray *branchGenParticle, Jet *jet, double jet_radius = 0.4) {
                                                                                                                                        
  for(int i=0; i<((TClonesArray*)branchGenParticle)->GetEntries(); i++){

    GenParticle *particle=(GenParticle*) branchGenParticle->At(i);

    if (jet->P4().DeltaR(particle->P4()) > jet_radius) continue;                                                                                                                                                        
    if (Rivet::PID::hasBottom(particle->PID)) return true;
                                                                                                                                                                                                  
  }           

  return false;

}

bool is_my_b_tag(Jet *jet, TClonesArray *branchGenParticle = nullptr, int seed = 0,  double jet_radius = 0.4, double eff = 0.9, double fake_eff = 0.01) {                                                                                                                                                                                                      

  bool is_b = ghost_btagging(branchGenParticle, jet, jet_radius);
  TRandom3 random_eff, random_fake_eff;

  random_eff.SetSeed(seed);
  random_fake_eff.SetSeed(seed);

  bool pass_eff = random_eff.Binomial( 1, eff ) > 0;
  bool pass_fake_eff = random_fake_eff.Binomial( 1, fake_eff ) > 0;

  if( is_b ) return pass_eff;
  else return pass_fake_eff;

}

//------------------------------------------------------------------------------------------------------------------------------------------------------------
// jets
//------------------------------------------------------------------------------------------------------------------------------------------------------------


vector <int> get_all_jets( string analysis_type = "reco", TClonesArray *branchJet = nullptr, double j_pT_min = 20 ) {

    vector <int> all_jet_indices;

    for( int i = 0; i<(int)branchJet->GetEntries(); i++ ) {

        Jet *jet = (Jet*) branchJet->At(i);
        if( jet->PT > j_pT_min )  all_jet_indices.push_back(i); 

    }

    return all_jet_indices;

}


//------------------------------------------------------------------------------------------------------------------------------------------------------------
// higgs - paired
//------------------------------------------------------------------------------------------------------------------------------------------------------------


// returns higgs
std::vector<std::pair< std::map<TString, float>, std::map<TString, std::vector<float>>>> get_higgs( string analysis_type = "reco", TClonesArray *branchJet = nullptr, TClonesArray *branchGenParticle = nullptr, TClonesArray *branchPFCand = nullptr ){

    std::vector<std::pair< std::map<TString, float>, std::map<TString, std::vector<float>>>> paired_b_jets;
    std::vector<std::pair< std::map<TString, float>, std::map<TString, std::vector<float>>>>  paired_jets = paired::PAIReDjointEvent( branchGenParticle, branchPFCand, branchJet, 0.4, false, false, true, 1.0, false );

    for(int i=0; i<(int)paired_jets.size(); i++){

      std::pair< std::map<TString, float>, std::map<TString, std::vector<float>>> paired_jet = paired_jets.at(i);

      if( paired_jet.first["isbtagged"] > 0 ) paired_b_jets.push_back(paired_jet);

    }

    return paired_b_jets;

}

pair< TLorentzVector, TLorentzVector> make_higgs( string analysis_type = "reco", std::vector<std::pair< std::map<TString, float>, std::map<TString, std::vector<float>>>> paired_b_jets = {} ) {

    TLorentzVector b1, b2;

    std::map<TString, float> b_jet;
    vector <int> b_jets;

    b_jet = paired_b_jets.at(0).first;

    b_jets.push_back(b_jet["jet1_index"]);
    b_jets.push_back(b_jet["jet2_index"]);

    b1.SetPtEtaPhiM(b_jet["jet1_pt"],b_jet["jet1_eta"],b_jet["jet1_phi"],b_jet["jet1_mass"]);
    b2.SetPtEtaPhiM(b_jet["jet2_pt"],b_jet["jet2_eta"],b_jet["jet2_phi"],b_jet["jet2_mass"]);
        
    return make_pair( b1, b2 );

}

/*
//------------------------------------------------------------------------------------------------------------------------------------------------------------
// higgs - pseudo b tagging
//------------------------------------------------------------------------------------------------------------------------------------------------------------

// returns vec containing b jets
vector <int> get_higgs( string analysis_type = "reco", TClonesArray *branchJet = nullptr, TClonesArray *branchGenParticle = nullptr ){
    
    double j_pT_min = 20.0;

    double eff = 1; //0.85;
    double fake_eff = 0; //0.02;

    vector <int> b_jets;
    vector <int> not_b_jets;

    for( int i = 0; i<(int)branchJet->GetEntries(); i++ ) {

        Jet *jet = (Jet*)branchJet->At(i);
        if( jet->PT < j_pT_min) continue;

        if( is_my_b_tag( jet, branchGenParticle, 0, 0.4, eff, fake_eff) ) b_jets.push_back(i); 
        else not_b_jets.push_back(i);
        
    }

    sort_by_pT( analysis_type, b_jets, branchJet, "jet");

    return b_jets;

}

pair< TLorentzVector, TLorentzVector> make_higgs( string analysis_type = "reco",  vector <int> b_jets = {-99}, TClonesArray *branchJet = nullptr ) {

    Jet *b1jet = (Jet*) branchJet->At( b_jets[0] );
    Jet *b2jet = (Jet*) branchJet->At( b_jets[1] );
    TLorentzVector b1 = b1jet -> P4();
    TLorentzVector b2 = b2jet -> P4();
    
    return make_pair(b1, b2);

}
*/

//------------------------------------------------------------------------------------------------------------------------------------------------------------
// vbf jets
//------------------------------------------------------------------------------------------------------------------------------------------------------------


// returns vec containing vbf jet candidates
vector <int> get_vbfjets( string analysis_type = "reco", TClonesArray *branchJet = nullptr, TClonesArray *branchGenParticle = nullptr ){
    
    double j_pT_min = 20.0;

    double eff = 1; //0.85;
    double fake_eff = 0; //0.02;

    vector <int> b_jets;
    vector <int> not_b_jets;

    // ADD ETA SEP REQ

    for( int i = 0; i<(int)branchJet->GetEntries(); i++ ) {

        Jet *jet = (Jet*)branchJet->At(i);
        if( jet->PT < j_pT_min) continue;

        if( is_my_b_tag( jet, branchGenParticle, 0, 0.4, eff, fake_eff) ) b_jets.push_back(i); 
        else not_b_jets.push_back(i);
        
    }

    sort_by_pT( analysis_type, not_b_jets, branchJet, "jet");

    return not_b_jets;

}

pair< TLorentzVector, TLorentzVector> make_vbfjets( string analysis_type = "reco", vector <int> not_b_jets = {-99}, TClonesArray *branchJet = nullptr) {

    Jet *jet1 = (Jet*) branchJet->At( not_b_jets[0] );
    Jet *jet2 = (Jet*) branchJet->At( not_b_jets[1] );
    TLorentzVector vbf_jet1 = jet1 -> P4();
    TLorentzVector vbf_jet2 = jet2 -> P4();
    
    return make_pair(vbf_jet1, vbf_jet2);

}


#endif