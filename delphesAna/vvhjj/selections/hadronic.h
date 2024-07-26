#ifndef HADRONIC_H
#define HADRONIC_H

#include "../includes/kinematics_include.h"
#include "../includes/helperfunctions_include.h"
// #include "../../common_includes/ghost_tagging.h"
// #include "../../common_includes/make_paired.h"
#include "../includes/cutflow_include.h"

#include "../includes/paired_include.h"
#include "../includes/ghosttag_include.h"


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
// all jets
//------------------------------------------------------------------------------------------------------------------------------------------------------------


vector <int> get_all_jets( string analysis_type = "reco", vector<int> true_leps = {0}, TClonesArray *branchJet = nullptr, TClonesArray *branchGenParticle = nullptr , TClonesArray *branchElectron = nullptr , TClonesArray *branchMuon = nullptr, double j_pT_min = 20 , bool debug_bool = false ) {

    vector <int> all_jet_indices, skimmed_jet_indices;

    if ( analysis_type == "reco" || analysis_type == "particle" ) {

        for( int i = 0; i<(int)branchJet->GetEntries(); i++ ) {

            Jet *genjet = (Jet*) branchJet->At(i);

            if( genjet->PT > j_pT_min ) all_jet_indices.push_back(i);

        }

        skimmed_jet_indices = rm_jetlep_overlap( analysis_type, all_jet_indices, true_leps, branchJet, branchGenParticle );
    
    } 

    if ( debug_bool && skimmed_jet_indices.size()>0 ) {

        for ( int i=0; i<skimmed_jet_indices.size(); i++ ) {

            Jet *jet = (Jet*) branchJet->At( skimmed_jet_indices[i] );

            debug_print( debug_bool, analysis_type + " j : pt = " + to_string(jet->PT) + " , eta = " + to_string(jet->Eta) + " , phi = " + to_string(jet->Phi) );

        }

    }

    return skimmed_jet_indices;

}

/*
vector <int> get_all_genjets( string analysis_type = "reco",  vector<int> leps_particle = {0}, TClonesArray *branchGenParticle = nullptr, TClonesArray *branchJet = nullptr, double j_pT_min = 20 , bool debug_bool = false ) {

    vector <int> all_genjet_indices;

    if ( analysis_type == "particle" ) {

        for( int i = 0; i<(int)branchJet->GetEntries(); i++ ) {

            Jet *genjet = (Jet*) branchJet->At(i);

            if( genjet->PT > j_pT_min ) all_genjet_indices.push_back(i);

        }

        rm_gen_jetlep_overlap( all_genjet_indices, leps_particle, branchJet, branchGenParticle );

        for( int i = 0; i<all_genjet_indices.size(); i++ ) {

            Jet *genjet_ = (Jet*) branchJet->At( all_genjet_indices[i] );
            debug_print( debug_bool, analysis_type + " jet " + " : pt = " + to_string(genjet_->PT) + " , eta = " + to_string(genjet_->Eta) + " , phi = " + to_string(genjet_->Phi));
        
        }
    
    }

    return all_genjet_indices;

}
*/

//------------------------------------------------------------------------------------------------------------------------------------------------------------
// higgs - paired
//------------------------------------------------------------------------------------------------------------------------------------------------------------

// returns higgs
std::vector<std::pair< std::map<TString, float>, std::map<TString, std::vector<float>>>> get_bb( string analysis_type = "reco", pair<int, int> vbf_jets = {0,0}, TClonesArray *branchJet = nullptr, TClonesArray *branchGenParticle = nullptr, TClonesArray *branchPFCand = nullptr ){

    std::vector<std::pair< std::map<TString, float>, std::map<TString, std::vector<float>>>> paired_b_jets;
    // std::vector<std::pair< std::map<TString, float>, std::map<TString, std::vector<float>>>> skimmed_paired_b_jets;

    if ( analysis_type == "reco" || analysis_type == "particle" ) { // redundant

        std::vector<std::pair< std::map<TString, float>, std::map<TString, std::vector<float>>>>  paired_jets = paired::PAIReDjointEvent( branchGenParticle, branchPFCand, branchJet, 0.4, false, false, true, 1.0, false );

        for(int i=0; i<(int)paired_jets.size(); i++) {

            std::pair< std::map<TString, float>, std::map<TString, std::vector<float>>> paired_jet = paired_jets.at(i);

         // cout << " isbtagged: " << paired_jet.first["isbtagged"] << " label_bb " << paired_jet.first["label_bb"] << endl;

            if( paired_jet.first["isbtagged"] > 0 ) {
                
                if ( vbf_jets.first == paired_jet.first["jet1_index"] || vbf_jets.second == paired_jet.first["jet1_index"] || vbf_jets.first ==  paired_jet.first["jet2_index"] || vbf_jets.second ==  paired_jet.first["jet2_index"]) continue; // check if already taken as vbf jet
                else paired_b_jets.push_back(paired_jet);

            }

        }
    
    }

    return paired_b_jets;

}


pair< TLorentzVector, TLorentzVector> make_bb( string analysis_type = "reco", std::vector<std::pair< std::map<TString, float>, std::map<TString, std::vector<float>>>> paired_b_jets = {} ) {

    TLorentzVector b1, b2;

    if ( analysis_type == "reco" || analysis_type == "particle" ) {

        vector < int > b_jets;

        std::map<TString, float> b_jet = paired_b_jets.at(0).first;

        b_jets.push_back(b_jet["jet1_index"]);
        b_jets.push_back(b_jet["jet2_index"]);

        b1.SetPtEtaPhiM(b_jet["jet1_pt"],b_jet["jet1_eta"],b_jet["jet1_phi"],b_jet["jet1_mass"]);
        b2.SetPtEtaPhiM(b_jet["jet2_pt"],b_jet["jet2_eta"],b_jet["jet2_phi"],b_jet["jet2_mass"]);

    } 

    return make_pair( b1, b2 );

}


//------------------------------------------------------------------------------------------------------------------------------------------------------------
// higgs - pseudo b tagging
//------------------------------------------------------------------------------------------------------------------------------------------------------------

/* 

// returns vec containing b jets
vector <int> get_bb( string analysis_type = "reco", TClonesArray *branchJet = nullptr, TClonesArray *branchGenParticle = nullptr ){
    
    double j_pT_min = 20.0;

    double eff = 1; //0.85;
    double fake_eff = 0; //0.02;

    vector <int> b_jets;
    vector <int> not_b_jets;

    for( int i = 0; i<(int)branchJet->GetEntries(); i++ ) {

        Jet *jet = (Jet*)branchJet->At(i);
        if( jet->PT < j_pT_min) continue;

        // if( is_my_b_tag( jet, branchGenParticle, 0, 0.4, eff, fake_eff) ) b_jets.push_back(i); 
        if( jet->BTag == 1) b_jets.push_back(i);
        else not_b_jets.push_back(i);
        
    }

    sort_by_pT( analysis_type, b_jets, branchJet, "jet");

    return b_jets;

}

pair< TLorentzVector, TLorentzVector> make_bb( string analysis_type = "reco",  vector <int> b_jets = {-99}, TClonesArray *branchJet = nullptr ) {

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
pair<int, int> get_jj( string analysis_type = "reco", vector<int> jets = {0}, TClonesArray *branchJet = nullptr, TClonesArray *branchGenParticle = nullptr ) {

    double deta_max = 0.0;

    pair<int, int> vbf_pair = make_pair( -1, -1 );

    if ( jets.size() > 1) {

        for (int i = 0; i < jets.size(); i++) {

            Jet* j1 = (Jet*)branchJet->At(jets[i]);

            for (int j = i + 1; j < jets.size(); j++) {

                Jet* j2 = (Jet*)branchJet->At(jets[j]);

                double deta = abs(j1->Eta - j2->Eta);

                if ( deta < 2.5 ) continue;
                if ( deta > deta_max ) {

                    deta_max = deta;

                    if (j1->PT > j2->PT) vbf_pair = make_pair(jets[i], jets[j]);
                    else vbf_pair = make_pair(jets[j], jets[i]);

                }

            }

        }

        // to suppress ttHbb
        if ( vbf_pair.first != -1 && vbf_pair.second != -1 ) {

            Jet* j1 = (Jet*)branchJet->At( vbf_pair.first ); Jet* j2 = (Jet*)branchJet->At( vbf_pair.second );
            
            if ( is_my_b_tag( j1, branchGenParticle, 0, 0.4, 0.87, 0.01) && is_my_b_tag( j2, branchGenParticle, 0, 0.4, 0.87, 0.01) ) vbf_pair = make_pair( -1, -1 ) ;

        }

    }

    return vbf_pair;

}


pair< TLorentzVector, TLorentzVector> make_jj( string analysis_type = "reco", pair<int, int> vbf_pair = {-99, -99}, TClonesArray *branchJet = nullptr) {

    Jet *jet1 = (Jet*) branchJet->At( vbf_pair.first );
    Jet *jet2 = (Jet*) branchJet->At( vbf_pair.second );
    TLorentzVector vbf_jet1 = jet1 -> P4();
    TLorentzVector vbf_jet2 = jet2 -> P4();
    
    return make_pair(vbf_jet1, vbf_jet2);

}


//------------------------------------------------------------------------------------------------------------------------------------------------------------
// parton ( both hbb jj)
//------------------------------------------------------------------------------------------------------------------------------------------------------------


pair < int, pair < int, int >> get_hbb_parton( TClonesArray *branchGenParticle = nullptr, bool debug_bool = false ) {

    bool found_h = false; bool found_hbb = false;

    pair<int, int> bb = make_pair(-99, -99);
    pair<int, pair<int, int>> hbb = make_pair(-99, bb);

    GenParticle *d1 = nullptr; GenParticle *d2 = nullptr; 

    for(int i=0; i<(int)branchGenParticle->GetEntries(); i++) {

        GenParticle *particle = (GenParticle*) branchGenParticle->At(i);

        if ( particle->D1 != -1 &&  particle->D2 != -1) {
        
            GenParticle *d1 = (GenParticle*) branchGenParticle->At( particle->D1 );
            GenParticle *d2 = (GenParticle*) branchGenParticle->At( particle->D2 );

            if ( particle->PID == 25 && d1->PID != 25 && d2->PID != 25 ) {

                found_h = true;

                if ( abs(d1->PID) == 5 && abs(d2->PID) == 5 ) {

                    found_hbb = true;

                    bb = make_pair( particle->D1, particle->D2 );
                    hbb = make_pair( i, bb );

                    debug_print( debug_bool, " " );
                    debug_print( debug_bool, "parton hbb m: " + to_string( particle->Mass ) + " hbb pT: " +  to_string( particle->PT ) );
                    debug_print( debug_bool, "parton b1 pT: " + to_string( d1->PT ) + " b1 eta: " + to_string( d1->Eta ) + " b1 phi: " + to_string( d1->Phi ) );
                    debug_print( debug_bool, "parton b2 pT: " + to_string( d2->PT ) + " b2 eta: " + to_string( d2->Eta ) + " b2 phi: " + to_string( d2->Phi ) );
                    debug_print( debug_bool, " " );

                }

            } 
        
    }
        
    }

    if ( found_h && found_hbb ) return hbb;
    else return make_pair(-99, make_pair(-99, -99));

}

// this might need some adjusting
pair < int, int > get_jj_parton( TClonesArray *branchGenParticle = nullptr, bool debug_bool = false  ) {

    vector < int > j_i; pair < int, int > jj = make_pair(-99,-99);

    for(int i=0; i<15; i++) {

        GenParticle* particle = (GenParticle*)branchGenParticle->At(i);

        if ( abs(particle->PID) == 1 ||  abs(particle->PID) == 2 ||  abs(particle->PID) == 3 ||  abs(particle->PID) == 4 ||  abs(particle->PID) == 5 ||  abs(particle->PID) == 6 ) {

            if ( j_i.size() == 2 ) break;
            if ( (particle->Status == 23 || particle->Status == 41 || particle->Status == 42) && abs(particle->Eta) < 100 ) {
                
                j_i.push_back( i );

            }
            
        }

    }

    if ( j_i.size() == 2) {

        GenParticle *j1 = (GenParticle*) branchGenParticle->At( j_i[0] );
        GenParticle *j2 = (GenParticle*) branchGenParticle->At( j_i[1] );

        if ( abs(j1->PID) == 5 && abs(j2->PID) == 5 ) return jj;

        debug_print( debug_bool, "parton j1 pT: " + to_string( j1->PT ) + " j1 eta: " + to_string( j1->Eta ) + " j1 phi: " + to_string( j1->Phi ) );
        debug_print( debug_bool, "parton j2 pT: " + to_string( j2->PT ) + " j2 eta: " + to_string( j2->Eta ) + " j2 phi: " + to_string( j2->Phi ) );
    
        jj = make_pair( j_i[0], j_i[1] );
    }

    return jj;

}

pair< TLorentzVector, pair< TLorentzVector,TLorentzVector >> make_parton_bb( string analysis_type = "parton",pair < int, pair < int, int >> hbb_parton = {-99,{-99, -99}}, TClonesArray *branchGenParticle = nullptr ) {

    TLorentzVector h = ((GenParticle*) branchGenParticle->At( hbb_parton.first )) -> P4();
    TLorentzVector b1 = ((GenParticle*) branchGenParticle->At( hbb_parton.second.first )) -> P4();
    TLorentzVector b2 = ((GenParticle*) branchGenParticle->At( hbb_parton.second.second )) -> P4();
    
    pair< TLorentzVector,TLorentzVector > bb = make_pair( b1, b2 );

    return make_pair( h, bb );

}


pair< TLorentzVector,TLorentzVector > make_parton_jj( string analysis_type = "parton",  pair < int, int > jj_parton = {-99, -99}, TClonesArray *branchGenParticle = nullptr ) {

    TLorentzVector j1 = ((GenParticle*) branchGenParticle->At( jj_parton.first )) -> P4();
    TLorentzVector j2 = ((GenParticle*) branchGenParticle->At( jj_parton.second )) -> P4();
    
    return make_pair( j1, j2 );

}



#endif