#ifndef LEPTONIC_H
#define LEPTONIC_H

#include "../includes/kinematics_include.h"
#include "../includes/helperfunctions_include.h"

//------------------------------------------------------------------------------------------------------------------------------------------------------------
// leptons
//------------------------------------------------------------------------------------------------------------------------------------------------------------

// returns a vector sorted by pT of all leptons in event
vector<int> get_all_leptons( string analysis_type = "reco", TClonesArray *branchElectron = nullptr, TClonesArray *branchMuon = nullptr, double pT_min = 15, double eta_max = 2.5 ) {

    vector <int> all_lep_indices;

    // select leptons

    if( branchElectron == nullptr ) return all_lep_indices;
    if( branchMuon == nullptr ) return all_lep_indices;

    for( int i=0; i<(int)branchElectron->GetEntries(); i++ ) {

        if ( analysis_type == "reco" ) { 

                Electron *e = (Electron *)branchElectron->At(i);
                Muon *mu = (Muon *)branchElectron->At(i);

                if ( e -> PT > pT_min && fabs( e -> Eta ) < eta_max ) all_lep_indices.push_back(i);
                if ( mu -> PT > pT_min && fabs( mu -> Eta ) < eta_max ) all_lep_indices.push_back(i);

        } else if ( analysis_type == "particle" ) {

            GenParticle *lep = (GenParticle *)branchElectron->At(i);

            if( lep->Status !=1 ) continue;

            if( abs( lep->PID ) != 11 ) continue;
            if( abs( lep->PID ) != 13 ) continue;


            if ( lep -> PT > pT_min && fabs( lep -> Eta ) < eta_max ) all_lep_indices.push_back(i);
              
        }

    }

    return all_lep_indices; 
  
}

// returns a vector sorted by pT of all passing leptons of a certain flavour and charge
vector<int> get_leptons( string analysis_type = "reco", TClonesArray *branchElectron = nullptr, double pT_min = 15, double eta_max = 2.5, string lep_flavor = "electron", int charge = -1 ) {

    vector <int> lep_indices;

    // select leptons

    if( branchElectron == nullptr ) return lep_indices;

    for( int i=0; i<(int)branchElectron->GetEntries(); i++ ) {

        if ( analysis_type == "reco" ){ 

            if ( lep_flavor == "electron" ) {

                Electron *lep = (Electron *)branchElectron->At(i);
                if ( lep -> PT > pT_min && fabs( lep -> Eta ) < eta_max && lep -> Charge == charge ) lep_indices.push_back(i);

            } else if ( lep_flavor == "muon" ) {

                Muon *lep = (Muon *)branchElectron->At(i);
                if ( lep -> PT > pT_min && fabs( lep -> Eta ) < eta_max && lep -> Charge == charge ) lep_indices.push_back(i);

            }
        }

        else if ( analysis_type == "particle" ){

            GenParticle *lep = (GenParticle *)branchElectron->At(i);

            if( lep->Status !=1 ) continue;

            if ( lep_flavor == "electron" ) {
                
                if ( charge == 1 ) {
                    if ( lep -> PT > pT_min && fabs( lep -> Eta ) < eta_max && lep->PID == -11) lep_indices.push_back(i);
                } else if ( charge == -1 ) {
                    if ( lep -> PT > pT_min && fabs( lep -> Eta ) < eta_max && lep->PID == 11) lep_indices.push_back(i);
                }

            } else if ( lep_flavor == "muon" ) {

                if ( charge == 1 ) {
                    if ( lep -> PT > pT_min && fabs( lep -> Eta ) < eta_max && lep->PID == -13) lep_indices.push_back(i);
                } else if ( charge == -1 ) {
                    if ( lep -> PT > pT_min && fabs( lep -> Eta ) < eta_max && lep->PID == 13) lep_indices.push_back(i);
                }
            }

        }

    }
  
    // sort the indices by pT
    sort_by_pT( analysis_type, lep_indices, branchElectron, lep_flavor );

    return lep_indices; 
  
}




//------------------------------------------------------------------------------------------------------------------------------------------------------------
// w
//------------------------------------------------------------------------------------------------------------------------------------------------------------


// returns w = < et, < w leps >>
pair< int, vector <int>> get_w_leptonic( string analysis_type = "reco", TClonesArray *branchElectron = nullptr, TClonesArray *branchMuon = nullptr ) {

    int event_type = -1;

    TLorentzVector l1, l2;

    // define w selection criteria
    double e_pT_min = 15.0;
    double e_eta_max = 2.5;
    double mu_pT_min = 15.0;
    double mu_eta_max = 2.5;

    vector <int> e_min_indices;
    vector <int> e_plus_indices;
    vector <int> mu_min_indices;
    vector <int> mu_plus_indices;
    
    // get e+ e- mu+ mu-
    e_min_indices = get_leptons( analysis_type, branchElectron, e_pT_min, e_eta_max, "electron", -1 );
    e_plus_indices = get_leptons( analysis_type, branchElectron, e_pT_min, e_eta_max, "electron", 1 );
    mu_min_indices = get_leptons( analysis_type, branchMuon, mu_pT_min, mu_eta_max, "muon", -1 );
    mu_plus_indices = get_leptons( analysis_type, branchMuon, mu_pT_min, mu_eta_max, "muon", 1 );

    vector <int> w_leps;
   
    if(  ( e_min_indices.size() + e_plus_indices.size() + mu_min_indices.size() + mu_plus_indices.size() ) >= 2 ) {

        if ( mu_min_indices.size() > 0 && mu_plus_indices.size() > 0 ) { // case mu- mu+

            event_type = 0;

            w_leps.push_back( mu_min_indices[0] );
            w_leps.push_back( mu_plus_indices[0] );

        } else if ( e_min_indices.size() > 0 && e_plus_indices.size() > 0 ) {  // case e- e+ 

            event_type = 1; 

            w_leps.push_back( e_min_indices[0] );
            w_leps.push_back( e_plus_indices[0] );

        } else if ( mu_min_indices.size() > 0 && e_plus_indices.size() > 0 ) { // case mu- e+

            event_type = 2;

            w_leps.push_back( mu_min_indices[0] );
            w_leps.push_back( e_plus_indices[0] );

        } else if ( e_min_indices.size() > 0 && mu_plus_indices.size() > 0 ) {  // case e- mu+

            event_type = 3; 

            w_leps.push_back( e_min_indices[0] );
            w_leps.push_back( mu_plus_indices[0] );

        }

    }

    if ( event_type == 0 ) get_leading_subleading( analysis_type, w_leps, branchMuon, branchMuon, "muon", "muon" );
    else if ( event_type == 1 ) get_leading_subleading( analysis_type, w_leps, branchElectron, branchElectron, "electron", "electron" );
    else if ( event_type == 2 ) get_leading_subleading( analysis_type, w_leps, branchMuon, branchElectron, "muon", "electron" );
    else if ( event_type == 3 ) get_leading_subleading( analysis_type, w_leps, branchElectron, branchMuon, "electron", "muon" );

    return make_pair( event_type, w_leps );

}


pair<TLorentzVector, TLorentzVector> make_w_leptonic( string analysis_type = "reco", pair< int, vector <int>> w_leps = {-99, {-99}}, TClonesArray *branchElectron = nullptr, TClonesArray *branchMuon = nullptr ,  TClonesArray *branchMissingET = nullptr ) {

    int event_type = w_leps.first;

    TLorentzVector w1, w2;

    TLorentzVector met = ( (MissingET*)branchMissingET->At( 0 ) ) -> P4();

    if ( analysis_type == "reco") {

        if ( event_type == 0 ) {

            TLorentzVector lep1 = ( (Muon*)branchMuon->At( w_leps.second[0] ) ) -> P4();
            TLorentzVector lep2 = ( (Muon*)branchMuon->At( w_leps.second[1] ) ) -> P4();
            w1 = lep1 + met;
            w2 = lep2 + met;

        } else if ( event_type == 1 ) {

            TLorentzVector lep1 = ( (Electron*)branchElectron->At( w_leps.second[0] ) ) -> P4();
            TLorentzVector lep2 = ( (Electron*)branchElectron->At( w_leps.second[1] ) ) -> P4();
            w1 = lep1 + met;
            w2 = lep2 + met;

        } else if ( event_type == 2 ) {

            TLorentzVector lep1 = ( (Muon*)branchMuon->At( w_leps.second[0] ) ) -> P4();
            TLorentzVector lep2 = ( (Electron*)branchElectron->At( w_leps.second[1] ) ) -> P4();
            w1 = lep1 + met;
            w2 = lep2 + met;

        } else if ( event_type == 3 ) {

            TLorentzVector lep1 = ( (Electron*)branchElectron->At( w_leps.second[0] ) ) -> P4();
            TLorentzVector lep2 = ( (Muon*)branchMuon->At( w_leps.second[1] ) ) -> P4();
            w1 = lep1 + met;
            w2 = lep2 + met;

        }

    } if ( analysis_type == "particle") {

        if ( event_type == 0 ) {

            TLorentzVector lep1 = ( (GenParticle*)branchMuon->At( w_leps.second[0] ) ) -> P4();
            TLorentzVector lep2 = ( (GenParticle*)branchMuon->At( w_leps.second[1] ) ) -> P4();
            w1 = lep1 + met;
            w2 = lep2 + met;

        } else if ( event_type == 1 ) {

            TLorentzVector lep1 = ( (GenParticle*)branchElectron->At( w_leps.second[0] ) ) -> P4();
            TLorentzVector lep2 = ( (GenParticle*)branchElectron->At( w_leps.second[1] ) ) -> P4();
            w1 = lep1 + met;
            w2 = lep2 + met;

        } else if ( event_type == 2 ) {

            TLorentzVector lep1 = ( (GenParticle*)branchMuon->At( w_leps.second[0] ) ) -> P4();
            TLorentzVector lep2 = ( (GenParticle*)branchElectron->At( w_leps.second[1] ) ) -> P4();
            w1 = lep1 + met;
            w2 = lep2 + met;

        } else if ( event_type == 3 ) {

            TLorentzVector lep1 = ( (GenParticle*)branchElectron->At( w_leps.second[0] ) ) -> P4();
            TLorentzVector lep2 = ( (GenParticle*)branchMuon->At( w_leps.second[1] ) ) -> P4();
            w1 = lep1 + met;
            w2 = lep2 + met;

        }

    }

    return make_pair( w1, w2 );

}

//------------------------------------------------------------------------------------------------------------------------------------------------------------
// z
//------------------------------------------------------------------------------------------------------------------------------------------------------------


// returns z = < et, < z1 leps, z2 leps>>
pair< int, pair< vector <int>, vector <int>>> get_z_leptonic( string analysis_type = "reco", TClonesArray *branchElectron = nullptr, TClonesArray *branchMuon = nullptr ) {

    int event_type = -1;

    TLorentzVector l1, l2, l3, l4;

    vector <int> e_min_indices;
    vector <int> e_plus_indices;
    vector <int> mu_min_indices;
    vector <int> mu_plus_indices;

    // define z selection criteria
    double e_pT_min = 15.0;
    double e_eta_max = 2.5;
    double mu_pT_min = 15.0;
    double mu_eta_max = 2.5;

    // get e+ e- mu+ mu-
    e_min_indices = get_leptons( analysis_type, branchElectron, e_pT_min, e_eta_max, "electron", -1 );
    e_plus_indices = get_leptons( analysis_type, branchElectron, e_pT_min, e_eta_max, "electron", 1 );
    mu_min_indices = get_leptons( analysis_type, branchMuon, mu_pT_min, mu_eta_max, "muon", -1 );
    mu_plus_indices = get_leptons( analysis_type, branchMuon, mu_pT_min, mu_eta_max, "muon", 1 );

    vector <int> e_indices;
    vector <int> mu_indices;
    concatenate_indices( e_min_indices, e_indices );
    concatenate_indices( e_plus_indices, e_indices );
    concatenate_indices( mu_min_indices, mu_indices );
    concatenate_indices( mu_plus_indices, mu_indices );

    sort_by_pT( analysis_type, e_indices, branchElectron, "electron" );
    sort_by_pT( analysis_type, mu_indices, branchMuon, "muon" );
 
    vector <int> z1_leps;
    vector <int> z2_leps;

    if(  ( e_min_indices.size() + e_min_indices.size() + mu_min_indices.size() + mu_plus_indices.size() ) >= 4 ) {

        if ( mu_min_indices.size() >= 2 && mu_plus_indices.size() >= 2 ) { // case mu- mu+ , mu- mu+

            event_type = 0;

            z1_leps.push_back( mu_min_indices[0] );
            z1_leps.push_back( mu_plus_indices[0] );
            z2_leps.push_back( mu_min_indices[1] );
            z2_leps.push_back( mu_plus_indices[1] );

        } else if ( e_min_indices.size() >= 2 && e_plus_indices.size() >= 2 ) {  // case e- e+ , e- e+

            event_type = 1; 

            z1_leps.push_back( e_min_indices[0] );
            z1_leps.push_back( e_plus_indices[0] );
            z2_leps.push_back( e_min_indices[1] );
            z2_leps.push_back( e_plus_indices[1] );

        } else if ( e_indices.size() >= 2 && mu_indices.size() >= 2 ) { // case mu- mu+ , e- e+ or case e- e+ , mu- mu+

            TLorentzVector e1 = ( (Electron *) branchElectron->At( e_indices[0] ) ) -> P4();
            TLorentzVector m1 = ( (Muon *) branchMuon->At( mu_indices[0] ) ) -> P4();

            if ( m1.Pt() > e1.Pt()) {

                event_type = 2;

                z1_leps.push_back( mu_min_indices[0] );
                z1_leps.push_back( mu_plus_indices[0] );
                z2_leps.push_back( e_min_indices[0] );
                z2_leps.push_back( e_plus_indices[0] );

            } else if ( e1.Pt() > m1.Pt()) {

                event_type = 3;

                z1_leps.push_back( e_min_indices[0] );
                z1_leps.push_back( e_plus_indices[0] );
                z2_leps.push_back( mu_min_indices[0] );
                z2_leps.push_back( mu_plus_indices[0] );

            }

        }  

    }

    return make_pair(event_type, make_pair(z1_leps, z2_leps));

}

pair<TLorentzVector, TLorentzVector> make_z_leptonic( string analysis_type = "reco", pair< int, pair< vector <int>, vector <int>>> z_leps = {-99, {{-99}, {-99}}}, TClonesArray *branchElectron = nullptr, TClonesArray *branchMuon = nullptr ) {

    int event_type = z_leps.first;
    vector <int> z1_leps = z_leps.second.first;
    vector <int> z2_leps = z_leps.second.first;

    TLorentzVector z1, z2;

    if ( analysis_type == "reco" ) {

        if ( event_type == 0 ) {

            TLorentzVector lep1 = ( (Muon*)branchMuon->At( z1_leps[0] ) ) -> P4();
            TLorentzVector lep2 = ( (Muon*)branchMuon->At( z1_leps[1] ) ) -> P4();
            TLorentzVector lep3 = ( (Muon*)branchMuon->At( z2_leps[0] ) ) -> P4();
            TLorentzVector lep4 = ( (Muon*)branchMuon->At( z2_leps[1] ) ) -> P4();
            z1 = lep1 + lep2;
            z2 = lep3 + lep4;

        } else if ( event_type == 1 ) {

            TLorentzVector lep1 = ( (Electron*)branchElectron->At( z1_leps[0] ) ) -> P4();
            TLorentzVector lep2 = ( (Electron*)branchElectron->At( z1_leps[1] ) ) -> P4();
            TLorentzVector lep3 = ( (Electron*)branchElectron->At( z2_leps[0] ) ) -> P4();
            TLorentzVector lep4 = ( (Electron*)branchElectron->At( z2_leps[1] ) ) -> P4();
            z1 = lep1 + lep2;
            z2 = lep3 + lep4;

        } else if ( event_type == 2 ) {

            TLorentzVector lep1 = ( (Muon*)branchMuon->At( z1_leps[0] ) ) -> P4();
            TLorentzVector lep2 = ( (Muon*)branchMuon->At( z1_leps[1] ) ) -> P4();
            TLorentzVector lep3 = ( (Electron*)branchElectron->At( z2_leps[0] ) ) -> P4();
            TLorentzVector lep4 = ( (Electron*)branchElectron->At( z2_leps[1] ) ) -> P4();
            z1 = lep1 + lep2;
            z2 = lep3 + lep4;

        } else if ( event_type == 3 ) {

            TLorentzVector lep1 = ( (Electron*)branchElectron->At( z1_leps[0] ) ) -> P4();
            TLorentzVector lep2 = ( (Electron*)branchElectron->At( z1_leps[1] ) ) -> P4();
            TLorentzVector lep3 = ( (Muon*)branchMuon->At( z2_leps[0] ) ) -> P4();
            TLorentzVector lep4 = ( (Muon*)branchMuon->At( z2_leps[1] ) ) -> P4();
            z1 = lep1 + lep2;
            z2 = lep3 + lep4;

        }

    } else if ( analysis_type == "particle" ) {

       if ( event_type == 0 ) {

            TLorentzVector lep1 = ( (GenParticle*)branchMuon->At( z1_leps[0] ) ) -> P4();
            TLorentzVector lep2 = ( (GenParticle*)branchMuon->At( z1_leps[1] ) ) -> P4();
            TLorentzVector lep3 = ( (GenParticle*)branchMuon->At( z2_leps[0] ) ) -> P4();
            TLorentzVector lep4 = ( (GenParticle*)branchMuon->At( z2_leps[1] ) ) -> P4();
            z1 = lep1 + lep2;
            z2 = lep3 + lep4;

        } else if ( event_type == 1 ) {

            TLorentzVector lep1 = ( (GenParticle*)branchElectron->At( z1_leps[0] ) ) -> P4();
            TLorentzVector lep2 = ( (GenParticle*)branchElectron->At( z1_leps[1] ) ) -> P4();
            TLorentzVector lep3 = ( (GenParticle*)branchElectron->At( z2_leps[0] ) ) -> P4();
            TLorentzVector lep4 = ( (GenParticle*)branchElectron->At( z2_leps[1] ) ) -> P4();
            z1 = lep1 + lep2;
            z2 = lep3 + lep4;

        } else if ( event_type == 2 ) {

            TLorentzVector lep1 = ( (GenParticle*)branchMuon->At( z1_leps[0] ) ) -> P4();
            TLorentzVector lep2 = ( (GenParticle*)branchMuon->At( z1_leps[1] ) ) -> P4();
            TLorentzVector lep3 = ( (GenParticle*)branchElectron->At( z2_leps[0] ) ) -> P4();
            TLorentzVector lep4 = ( (GenParticle*)branchElectron->At( z2_leps[1] ) ) -> P4();
            z1 = lep1 + lep2;
            z2 = lep3 + lep4;

        } else if ( event_type == 3 ) {

            TLorentzVector lep1 = ( (GenParticle*)branchElectron->At( z1_leps[0] ) ) -> P4();
            TLorentzVector lep2 = ( (GenParticle*)branchElectron->At( z1_leps[1] ) ) -> P4();
            TLorentzVector lep3 = ( (GenParticle*)branchMuon->At( z2_leps[0] ) ) -> P4();
            TLorentzVector lep4 = ( (GenParticle*)branchMuon->At( z2_leps[1] ) ) -> P4();
            z1 = lep1 + lep2;
            z2 = lep3 + lep4;

        }
    
    }

    return make_pair( z1, z2 );

}

#endif