#ifndef LEPTONIC_H
#define LEPTONIC_H

#include "../includes/kinematics_include.h"
#include "../includes/helperfunctions_include.h"
#include "../includes/cutflow_include.h"

//------------------------------------------------------------------------------------------------------------------------------------------------------------
// leptons
//------------------------------------------------------------------------------------------------------------------------------------------------------------


// returns a vector of all passing leptons of a certain flavour and charge
vector<int> get_leptons( string analysis_type = "reco", TClonesArray *branchElectron = nullptr, TClonesArray *branchGenParticle = nullptr, TClonesArray *branchJet = nullptr, double pT_min = 15, double eta_max = 2.5, string lep_flavor = "electron", int charge = -1, bool debug_bool = false ) {

    vector <int> lep_indices; 

    // select leptons
    for( int i=0; i<(int)branchElectron->GetEntries(); i++ ) {

        if ( analysis_type == "reco" ){ 

            if ( lep_flavor == "electron" ) {

                Electron *lep = (Electron *)branchElectron->At(i);

                if ( (lep -> PT > pT_min) && (fabs( lep -> Eta ) < eta_max) && (lep -> Charge == charge) ) {

                    lep_indices.push_back(i);
                    if (lep->Charge == -1) debug_print( debug_bool, analysis_type + " e- : pt = " + to_string(lep->PT) + " , eta = " + to_string(lep->Eta) + " , phi = " + to_string(lep->Phi) );
                    if (lep->Charge == 1) debug_print( debug_bool, analysis_type + " e+ : pt = " + to_string(lep->PT) + " , eta = " + to_string(lep->Eta) + " , phi = " + to_string(lep->Phi) );

                }

            } else if ( lep_flavor == "muon" ) {

                Muon *lep = (Muon *)branchElectron->At(i);

                if ( (lep -> PT > pT_min) && (fabs( lep -> Eta ) < eta_max) && (lep -> Charge == charge) ) {
    
                    lep_indices.push_back(i);
                    if (lep->Charge == -1) debug_print( debug_bool, analysis_type + " mu- : pt = " + to_string(lep->PT) + " , eta = " + to_string(lep->Eta) + " , phi = " + to_string(lep->Phi) );
                    if (lep->Charge == 1) debug_print( debug_bool, analysis_type + " mu+ : pt = " + to_string(lep->PT) + " , eta = " + to_string(lep->Eta) + " , phi = " + to_string(lep->Phi) );

                }

            }

        } else if ( analysis_type == "particle" ){

            GenParticle *lep = (GenParticle *)branchElectron->At(i);

            if ( analysis_type == "particle" ) {
                if( lep->Status !=1 ) continue;
                if( !is_prompt( i, branchElectron )) continue;
            }

            if ( lep_flavor == "electron" ) {
                
                if ( charge == -1 ) {

                    if ( (lep->PT > pT_min) && (fabs( lep->Eta ) < eta_max) && (lep->PID == 11)) {
                        
                        lep_indices.push_back(i);
                        debug_print( debug_bool, analysis_type + " e- : pt = " + to_string(lep->PT) + " , eta = " + to_string(lep->Eta) + " , phi = " + to_string(lep->Phi) );

                    }

                } else if ( charge == 1 ) {

                    if ( (lep->PT > pT_min) && (fabs( lep->Eta ) < eta_max) && (lep->PID == -11)) {
                        
                        lep_indices.push_back(i);
                        debug_print( debug_bool, analysis_type + " e+ : pt = " + to_string(lep->PT) + " , eta = " + to_string(lep->Eta) + " , phi = " + to_string(lep->Phi) );

                    }
                }

            } else if ( lep_flavor == "muon" ) {

                if ( charge == -1 ) {

                    if ( (lep->PT > pT_min) && (fabs( lep->Eta ) < eta_max) && (lep->PID == 13)) {
                        
                        lep_indices.push_back(i);
                        debug_print( debug_bool, analysis_type + " mu- : pt = " + to_string(lep->PT) + " , eta = " + to_string(lep->Eta) + " , phi = " + to_string(lep->Phi) );
                    
                    }

                } else if ( charge == 1 ) {

                    if ( (lep->PT > pT_min) && (fabs( lep->Eta ) < eta_max) && (lep->PID == -13)) {
                        
                        lep_indices.push_back(i);
                        debug_print( debug_bool, analysis_type + " mu+ : pt = " + to_string(lep->PT) + " , eta = " + to_string(lep->Eta) + " , phi = " + to_string(lep->Phi) );
                    
                    }

                }

            }

        } 

    }

    return lep_indices; 
  
}



vector<int> get_jet_leptons( string analysis_type = "reco", TClonesArray *branchElectron = nullptr, TClonesArray *branchGenParticle = nullptr, TClonesArray *branchJet = nullptr, vector <int> leps_parton = {-99}, vector <int> true_leps = {-99}, double pT_min = 15, double eta_max = 2.5, string lep_flavor = "electron", int charge = -1, bool debug_bool = false ) {

    vector < int > jet_leps;
    vector < int > skimmed_jet_leps;

    for( int i=0; i<(int)branchJet->GetEntries(); i++ ) {

        Jet* jet = (Jet*) branchJet->At( i );

        // if ( jet->Charge != charge ) continue;
        if ( jet->PT < pT_min || jet->Eta > eta_max ) continue;

        for (int j = 0; j < leps_parton.size(); j++) {

            GenParticle* genlep = (GenParticle*) branchGenParticle->At( leps_parton[j] );

            if ( sqrt( pow( jet->Eta - genlep->Eta, 2 ) + pow( jet->Phi - genlep->Phi, 2 ) ) > 0.4 ) continue; // matching ( > )

            if ( lep_flavor == "electron" && abs(genlep->PID) == 11 ) jet_leps.push_back(i);
            else if ( lep_flavor == "muon" && abs(genlep->PID) == 13 ) jet_leps.push_back(i);

        }

    }

    skimmed_jet_leps = rm_jetlep_overlap( analysis_type , jet_leps, true_leps, branchJet, branchGenParticle );

    for ( int i=0; i<skimmed_jet_leps.size(); i++ ){

            Jet* lep = (Jet*) branchJet->At( skimmed_jet_leps[i] );

            if ( lep_flavor == "electron" && lep->Charge == -1 ) debug_print( debug_bool, analysis_type + " (j) e- : pt = " + to_string(lep->PT) + " , eta = " + to_string(lep->Eta) + " , phi = " + to_string(lep->Phi) );
            if ( lep_flavor == "electron" && lep->Charge == 1 ) debug_print( debug_bool, analysis_type + " (j) e+ : pt = " + to_string(lep->PT) + " , eta = " + to_string(lep->Eta) + " , phi = " + to_string(lep->Phi) );
            if ( lep_flavor == "muon" && lep->Charge == -1 ) debug_print( debug_bool, analysis_type + " (j) mu- : pt = " + to_string(lep->PT) + " , eta = " + to_string(lep->Eta) + " , phi = " + to_string(lep->Phi) );
            if ( lep_flavor == "muon" && lep->Charge == 1 ) debug_print( debug_bool, analysis_type + " (j) mu+ : pt = " + to_string(lep->PT) + " , eta = " + to_string(lep->Eta) + " , phi = " + to_string(lep->Phi) );

        }

    return skimmed_jet_leps;

}



//------------------------------------------------------------------------------------------------------------------------------------------------------------
// w
//------------------------------------------------------------------------------------------------------------------------------------------------------------


// returns w = < et, < w leps >>
/* pair< int, vector <int>> get_ww_leptonic( string analysis_type = "reco", vector < int > eminus = {-99}, vector < int > eplus = {-99} , vector < int > muminus = {-99} , vector < int > muplus = {-99} ) {

    int event_type = -1;

    vector <int> w_leps;
   
    if(  ( eminus.size() + eplus.size() + muminus.size() + muplus.size() ) >= 2 ) {

        if ( muminus.size() > 0 && muplus.size() > 0 ) { // case mu- mu+

            event_type = 0;

            w_leps.push_back( muminus[0] );
            w_leps.push_back( muplus[0] );

        } else if ( eminus.size() > 0 && eplus.size() > 0 ) {  // case e- e+ 

            event_type = 1; 

            w_leps.push_back( eminus[0] );
            w_leps.push_back( eplus[0] );

        } else if ( muminus.size() > 0 && eplus.size() > 0 ) { // case mu- e+

            event_type = 2;

            w_leps.push_back( muminus[0] );
            w_leps.push_back( eplus[0] );

        } else if ( eminus.size() > 0 && muplus.size() > 0 ) {  // case e- mu+

            event_type = 3; 

            w_leps.push_back( eminus[0] );
            w_leps.push_back( muplus[0] );

        }

    }

    return make_pair( event_type, w_leps );

}
*/


pair< int, pair < vector < TLorentzVector > , vector < TLorentzVector > > > get_vv_leptonic( string analysis_type = "reco", string v = "w", vector < int > electrons = {-99}, vector < int > muons = {-99}, vector < int > jet_e = {-99}, vector < int > jet_mu = {-99}, TClonesArray *branchElectron = nullptr, TClonesArray *branchMuon = nullptr, TClonesArray *branchJet = nullptr, TClonesArray *branchMissingET = nullptr) {

    int event_type = -1;

    double target_mass; if ( v == "w" ) target_mass = 80; if ( v == "z" ) target_mass = 90;

    vector <pair < int, TLorentzVector >> leps;

    if ( analysis_type == "reco" ) {
           
        if( electrons.size() > 0 ) { // e = 1

            for ( int i=0; i<electrons.size(); i++ ) leps.push_back( make_pair( 1, ((Electron*) branchElectron->At( electrons[i] ))->P4() ));

        }

        if( muons.size() > 0 ) { // mu = 0

            for ( int i=0; i<muons.size(); i++ ) leps.push_back( make_pair( 0, ((Muon*) branchMuon->At( muons[i] ))->P4() )); 

        }

    } else if ( analysis_type == "particle" ) {

        if( electrons.size() > 0 ) { // e = 1

            for ( int i=0; i<electrons.size(); i++ ) leps.push_back( make_pair( 1, ((GenParticle*) branchElectron->At( electrons[i] ))->P4() ));
        }

        if( muons.size() > 0 ) { // mu = 0

            for ( int i=0; i<muons.size(); i++ ) leps.push_back( make_pair( 0, ((GenParticle*) branchMuon->At( muons[i] ))->P4() ));

        }

    }

    if( jet_e.size() > 0 ) { // e = 1

        for ( int i=0; i<jet_e.size(); i++ ) leps.push_back( make_pair( 1, ((Jet*) branchJet->At( jet_e[i] ))->P4() )); 

    }

    if( jet_mu.size() > 0 ) { // mu = 0

        for ( int i=0; i<jet_mu.size(); i++ ) leps.push_back( make_pair( 0, ((Jet*) branchJet->At( jet_mu[i] ))->P4() )); 

    }

    vector < TLorentzVector > selected_ls;  vector < TLorentzVector > selected_vs;


    if ( v == "w" && leps.size() > 0 ) {
  
        TLorentzVector met;

        // met option i
        // met = ( (MissingET*)branchMissingET->At( 0 ) ) -> P4();

        // met option ii
        MissingET *missinget = (MissingET*)branchMissingET->At( 0 );
        met.SetPtEtaPhiM((missinget->MET)/2,0,missinget->Phi,0);  

        pair < double, TLorentzVector > leading_lep;  pair < double, TLorentzVector > subleading_lep;
        double smallest_leading_mass_diff = 999.99; double smallest_subleading_mass_diff = 999.99;

        if (  leps.size() == 1 ) leading_lep = make_pair( 999, leps[0].second );

        if (  leps.size() > 1 ) {

            for (size_t i = 0; i < leps.size(); ++i) {

                TLorentzVector w = leps[i].second + met;
                double mass_diff = abs(w.M() - target_mass); 

                if (mass_diff < smallest_leading_mass_diff) {

                    smallest_subleading_mass_diff = smallest_leading_mass_diff;
                    subleading_lep = leading_lep;

                    smallest_leading_mass_diff = mass_diff;
                    leading_lep = make_pair( leps[i].first, leps[i].second);

                } else if (mass_diff < smallest_subleading_mass_diff) {

                    smallest_subleading_mass_diff = mass_diff;
                    subleading_lep = make_pair( leps[i].first, leps[i].second);

                }

            }

        }

        if ( (leading_lep.second + subleading_lep.second).M() > 10  ) {

            if ( leading_lep.first == 0 && subleading_lep.first == 0) event_type = 0;
            if ( leading_lep.first == 1 && subleading_lep.first == 1) event_type = 1;
            if ( leading_lep.first == 0 && subleading_lep.first == 1) event_type = 2;
            if ( leading_lep.first == 1 && subleading_lep.first == 0) event_type = 3;

            selected_ls.push_back( leading_lep.second ); selected_ls.push_back( subleading_lep.second );
            selected_vs.push_back( leading_lep.second + met ); selected_vs.push_back( subleading_lep.second + met );

        }

    }


    if (v == "z" && leps.size() > 1) {

        pair<double, TLorentzVector> z1_l1, z1_l2, z2_l1, z2_l2;

        double smallest_mass_diff = 999.99; double second_smallest_mass_diff = 999.99;

        for (size_t i = 0; i < leps.size(); ++i) {

            for (size_t j = i + 1; j < leps.size(); ++j) {

                if (leps[i].first == leps[j].first) {

                    TLorentzVector z_candidate = leps[i].second + leps[j].second;
                    double mass_diff = abs(z_candidate.M() - 91.1876); // Mass of Z boson in GeV/c^2

                    if (mass_diff < smallest_mass_diff) {

                        second_smallest_mass_diff = smallest_mass_diff;
                        z2_l1 = z1_l1;
                        z2_l2 = z1_l2;

                        smallest_mass_diff = mass_diff;

                        z1_l1 = make_pair(leps[i].first, leps[i].second);
                        z1_l2 = make_pair(leps[j].first, leps[j].second);

                    } else if (mass_diff < second_smallest_mass_diff) {

                        second_smallest_mass_diff = mass_diff;
                        z2_l1 = make_pair(leps[i].first, leps[i].second);
                        z2_l2 = make_pair(leps[j].first, leps[j].second);

                    }

                }

            }

        }

        if ( (z1_l1.second + z1_l2.second).M() > 10 && (z2_l1.second + z2_l2.second).M() > 10 ) {

            if ( ( z1_l1.first == 0 && z1_l2.first == 0 ) && ( z2_l1.first == 0 && z2_l2.first == 0 ) ) event_type = 0;
            if ( ( z1_l1.first == 1 && z1_l2.first == 1 ) && ( z2_l1.first == 1 && z2_l2.first == 1 ) ) event_type = 1;
            if ( ( z1_l1.first == 0 && z1_l2.first == 0 ) && ( z2_l1.first == 1 && z2_l2.first == 1 ) ) event_type = 2;
            if ( ( z1_l1.first == 1 && z1_l2.first == 1 ) && ( z2_l1.first == 0 && z2_l2.first == 0 ) ) event_type = 3;

            selected_ls.push_back( z1_l1.second ); selected_ls.push_back( z1_l2.second ); selected_ls.push_back( z2_l1.second ); selected_ls.push_back( z2_l2.second );
            selected_vs.push_back( z1_l1.second + z1_l2.second ); selected_vs.push_back( z2_l1.second + z2_l2.second );

        }

    }

    return make_pair( event_type, make_pair ( selected_vs, selected_ls ));

}

/*
pair< vector<TLorentzVector> , vector<TLorentzVector>> make_vv_leptonic( string analysis_type = "reco", string v = "w", pair< int, vector <int>> w_leps = {-99, {-99}}, TClonesArray *branchElectron = nullptr, TClonesArray *branchMuon = nullptr, TClonesArray *branchMissingET = nullptr ) {

    int event_type = w_leps.first;

    TLorentzVector met, lep1, lep2, w1, w2;
    vector < TLorentzVector > leps, ws;

    // option i
    // met = ( (MissingET*)branchMissingET->At( 0 ) ) -> P4();

    // option ii
    MissingET *missinget = (MissingET*)branchMissingET->At( 0 );
    met.SetPtEtaPhiM((missinget->MET)/2,0,missinget->Phi,0);


    if ( analysis_type == "reco") {

        if ( event_type == 0 ) {

            lep1 = ( (Muon*)branchMuon->At( w_leps.second[0] ) ) -> P4();
            lep2 = ( (Muon*)branchMuon->At( w_leps.second[1] ) ) -> P4();

        } else if ( event_type == 1 ) {

            lep1 = ( (Electron*)branchElectron->At( w_leps.second[0] ) ) -> P4();
            lep2 = ( (Electron*)branchElectron->At( w_leps.second[1] ) ) -> P4();

        } else if ( event_type == 2 ) {

            lep1 = ( (Muon*)branchMuon->At( w_leps.second[0] ) ) -> P4();
            lep2 = ( (Electron*)branchElectron->At( w_leps.second[1] ) ) -> P4();

        } else if ( event_type == 3 ) {

            lep1 = ( (Electron*)branchElectron->At( w_leps.second[0] ) ) -> P4();
            lep2 = ( (Muon*)branchMuon->At( w_leps.second[1] ) ) -> P4();

        }

    } if ( analysis_type == "particle") {

        if ( event_type == 0 ) {

            lep1 = ( (GenParticle*)branchMuon->At( w_leps.second[0] ) ) -> P4();
            lep2 = ( (GenParticle*)branchMuon->At( w_leps.second[1] ) ) -> P4();

        } else if ( event_type == 1 ) {

            lep1 = ( (GenParticle*)branchElectron->At( w_leps.second[0] ) ) -> P4();
            lep2 = ( (GenParticle*)branchElectron->At( w_leps.second[1] ) ) -> P4();

        } else if ( event_type == 2 ) {

            lep1 = ( (GenParticle*)branchMuon->At( w_leps.second[0] ) ) -> P4();
            lep2 = ( (GenParticle*)branchElectron->At( w_leps.second[1] ) ) -> P4();

        } else if ( event_type == 3 ) {

            lep1 = ( (GenParticle*)branchElectron->At( w_leps.second[0] ) ) -> P4();
            lep2 = ( (GenParticle*)branchMuon->At( w_leps.second[1] ) ) -> P4();

        }

    }

    if ( lep1.Pt() > lep2.Pt() ) {

        w1 = lep1 + met; w2 = lep2 + met;
        leps.push_back(lep1); leps.push_back(lep2);
    
    } else if ( lep1.Pt() < lep2.Pt() ) {
        
        w1 = lep2 + met; w2 = lep1 + met;
        leps.push_back(lep2); leps.push_back(lep1);
    
    }

    ws.push_back(w1); ws.push_back(w2);

    return make_pair( leps, ws );

}
*/

//------------------------------------------------------------------------------------------------------------------------------------------------------------
// z
//------------------------------------------------------------------------------------------------------------------------------------------------------------

/*
// returns z = < et, < z1 leps, z2 leps>>
pair< int, pair< vector <int>, vector <int>>> get_zz_leptonic( string analysis_type = "reco", vector < int > eminus = {0} , vector < int > eplus = {0} , vector < int > muminus = {0} , vector < int > muplus = {0}, TClonesArray *branchElectron = nullptr, TClonesArray *branchMuon = nullptr , TClonesArray *branchJet = nullptr ) {

    int event_type = -1;

    TLorentzVector l1, l2, l3, l4;

    vector <int> e;
    vector <int> mu;
    concatenate_indices( eminus, e );
    concatenate_indices( eplus, e );
    concatenate_indices( muminus, mu );
    concatenate_indices( muplus, mu );

    // sort_by_pT( analysis_type, e, branchElectron, "electron" );
    // sort_by_pT( analysis_type, mu, branchMuon, "muon" );
 
    vector <int> z1_leps;
    vector <int> z2_leps;

    if(  ( eminus.size() + eplus.size() + muminus.size() + muplus.size() ) >= 4 ) {

        if ( muminus.size() >= 2 && muplus.size() >= 2 ) { // case mu- mu+ , mu- mu+

            event_type = 0;

            z1_leps.push_back( muminus[0] );
            z1_leps.push_back( muplus[0] );
            z2_leps.push_back( muminus[1] );
            z2_leps.push_back( muplus[1] );

        } else if ( eminus.size() >= 2 && eplus.size() >= 2 ) {  // case e- e+ , e- e+

            event_type = 1; 

            z1_leps.push_back( eminus[0] );
            z1_leps.push_back( eplus[0] );
            z2_leps.push_back( eminus[1] );
            z2_leps.push_back( eplus[1] );

        } else if ( eminus.size() >= 1 && eplus.size() >= 1 && muminus.size() >= 1 && muplus.size() >= 1 ) { // case mu- mu+ , e- e+ or case e- e+ , mu- mu+

            TLorentzVector e_m = ( (Electron *) branchElectron->At( eminus[0] ) ) -> P4();
            TLorentzVector e_p = ( (Electron *) branchElectron->At( eplus[0] ) ) -> P4();
            TLorentzVector mu_m = ( (Muon *) branchMuon->At( muminus[0] ) ) -> P4();
            TLorentzVector mu_p = ( (Muon *) branchMuon->At( muplus[0] ) ) -> P4();

            if ( ( mu_m + mu_p ).Pt() > ( e_m + e_p ).Pt() ) {

                event_type = 2;

                z1_leps.push_back( muminus[0] );
                z1_leps.push_back( muplus[0] );
                z2_leps.push_back( eminus[0] );
                z2_leps.push_back( eplus[0] );

            } else if ( ( mu_m + mu_p ).Pt() < ( e_m + e_p ).Pt() ) {

                event_type = 3;

                z1_leps.push_back( eminus[0] );
                z1_leps.push_back( eplus[0] );
                z2_leps.push_back( muminus[0] );
                z2_leps.push_back( muplus[0] );

            }

        }  

    }

    return make_pair(event_type, make_pair(z1_leps, z2_leps));

}

pair< vector<TLorentzVector> , vector<TLorentzVector>> make_zz_leptonic( string analysis_type = "reco", pair< int, pair< vector <int>, vector <int>>> z_leps = {-99, {{-99}, {-99}}}, TClonesArray *branchElectron = nullptr, TClonesArray *branchMuon = nullptr ) {

    int event_type = z_leps.first;
    vector <int> z1_leps = z_leps.second.first;
    vector <int> z2_leps = z_leps.second.first;

    vector <TLorentzVector> leps, zs;
    TLorentzVector lep1, lep2, lep3, lep4, z1, z2;

    if ( analysis_type == "reco" ) {

        if ( event_type == 0 ) {

            lep1 = ( (Muon*)branchMuon->At( z1_leps[0] ) ) -> P4();
            lep2 = ( (Muon*)branchMuon->At( z1_leps[1] ) ) -> P4();
            lep3 = ( (Muon*)branchMuon->At( z2_leps[0] ) ) -> P4();
            lep4 = ( (Muon*)branchMuon->At( z2_leps[1] ) ) -> P4();
            z1 = lep1 + lep2;
            z2 = lep3 + lep4;

        } else if ( event_type == 1 ) {

            lep1 = ( (Electron*)branchElectron->At( z1_leps[0] ) ) -> P4();
            lep2 = ( (Electron*)branchElectron->At( z1_leps[1] ) ) -> P4();
            lep3 = ( (Electron*)branchElectron->At( z2_leps[0] ) ) -> P4();
            lep4 = ( (Electron*)branchElectron->At( z2_leps[1] ) ) -> P4();
            z1 = lep1 + lep2;
            z2 = lep3 + lep4;

        } else if ( event_type == 2 ) {

            lep1 = ( (Muon*)branchMuon->At( z1_leps[0] ) ) -> P4();
            lep2 = ( (Muon*)branchMuon->At( z1_leps[1] ) ) -> P4();
            lep3 = ( (Electron*)branchElectron->At( z2_leps[0] ) ) -> P4();
            lep4 = ( (Electron*)branchElectron->At( z2_leps[1] ) ) -> P4();
            z1 = lep1 + lep2;
            z2 = lep3 + lep4;

        } else if ( event_type == 3 ) {

            lep1 = ( (Electron*)branchElectron->At( z1_leps[0] ) ) -> P4();
            lep2 = ( (Electron*)branchElectron->At( z1_leps[1] ) ) -> P4();
            lep3 = ( (Muon*)branchMuon->At( z2_leps[0] ) ) -> P4();
            lep4 = ( (Muon*)branchMuon->At( z2_leps[1] ) ) -> P4();
            z1 = lep1 + lep2;
            z2 = lep3 + lep4;

        }

    } else if ( analysis_type == "particle" ) {

       if ( event_type == 0 ) {

            lep1 = ( (GenParticle*)branchMuon->At( z1_leps[0] ) ) -> P4();
            lep2 = ( (GenParticle*)branchMuon->At( z1_leps[1] ) ) -> P4();
            lep3 = ( (GenParticle*)branchMuon->At( z2_leps[0] ) ) -> P4();
            lep4 = ( (GenParticle*)branchMuon->At( z2_leps[1] ) ) -> P4();
            z1 = lep1 + lep2;
            z2 = lep3 + lep4;

        } else if ( event_type == 1 ) {

            lep1 = ( (GenParticle*)branchElectron->At( z1_leps[0] ) ) -> P4();
            lep2 = ( (GenParticle*)branchElectron->At( z1_leps[1] ) ) -> P4();
            lep3 = ( (GenParticle*)branchElectron->At( z2_leps[0] ) ) -> P4();
            lep4 = ( (GenParticle*)branchElectron->At( z2_leps[1] ) ) -> P4();
            z1 = lep1 + lep2;
            z2 = lep3 + lep4;

        } else if ( event_type == 2 ) {

            lep1 = ( (GenParticle*)branchMuon->At( z1_leps[0] ) ) -> P4();
            lep2 = ( (GenParticle*)branchMuon->At( z1_leps[1] ) ) -> P4();
            lep3 = ( (GenParticle*)branchElectron->At( z2_leps[0] ) ) -> P4();
            lep4 = ( (GenParticle*)branchElectron->At( z2_leps[1] ) ) -> P4();
            z1 = lep1 + lep2;
            z2 = lep3 + lep4;

        } else if ( event_type == 3 ) {

            lep1 = ( (GenParticle*)branchElectron->At( z1_leps[0] ) ) -> P4();
            lep2 = ( (GenParticle*)branchElectron->At( z1_leps[1] ) ) -> P4();
            lep3 = ( (GenParticle*)branchMuon->At( z2_leps[0] ) ) -> P4();
            lep4 = ( (GenParticle*)branchMuon->At( z2_leps[1] ) ) -> P4();
            z1 = lep1 + lep2;
            z2 = lep3 + lep4;

        }
    
    }

    leps.push_back(lep1); leps.push_back(lep2); leps.push_back(lep3); leps.push_back(lep4);
    zs.push_back(z1); zs.push_back(z2);

    return make_pair( leps, zs );

}
*/

//------------------------------------------------------------------------------------------------------------------------------------------------------------
// parton ( both w and z )
//------------------------------------------------------------------------------------------------------------------------------------------------------------


pair < int, pair < vector < int > , vector < int >>> get_vv_parton( TClonesArray *branchGenParticle = nullptr, bool debug_bool = false  ) { // < et , < < vs >, < ls > >

    vector <int> vs; vector <int> ls;
    pair < vector <int> , vector <int> > vvleps;
       
    int et = -1;

    GenParticle *l1; GenParticle *l2; GenParticle *l3; GenParticle *l4;

    bool is_w = false; bool is_z = false;

    for(int i=0; i<(int)branchGenParticle->GetEntries(); i++){

        GenParticle *particle=(GenParticle*)branchGenParticle->At(i);

        if (particle->D1 != -1 && particle->D2 != -1) {

            GenParticle *daughter1 = (GenParticle*) branchGenParticle->At(particle->D1);
            GenParticle *daughter2 = (GenParticle*) branchGenParticle->At(particle->D2);

            if ( abs(daughter1->PID) == 11 || abs(daughter1->PID) == 13 || abs(daughter2->PID) == 11 || abs(daughter2->PID) == 13 ) {

                if( abs(particle->PID) == 24 && abs(daughter1->PID) != 24 && abs(daughter2->PID) != 24 ) {
                    
                    is_w = true; vs.push_back(i);
                
                }

                if( abs(particle->PID) == 23 && abs(daughter1->PID) != 23 && abs(daughter2->PID) != 23 ) {
                    
                is_z = true; vs.push_back(i);
                
                }

            }

        }

    }

    if( vs.size() >= 2 ) {
            
        sort_by_pT( "parton", vs, branchGenParticle, "");

        // get w daughters + pid
        GenParticle *v1 = (GenParticle*) branchGenParticle->At( vs[0] ); GenParticle *v2 = (GenParticle*) branchGenParticle->At( vs[1] );
        GenParticle *v1d1 = (GenParticle*) branchGenParticle->At( v1->D1 ); GenParticle *v2d1 = (GenParticle*) branchGenParticle->At( v2->D1 );
        GenParticle *v1d2 = (GenParticle*) branchGenParticle->At( v1->D2 ); GenParticle *v2d2 = (GenParticle*) branchGenParticle->At( v2->D2 );

        // get leps
        if( abs( v1d1 -> PID ) == 11 || abs( v1d1 -> PID ) == 13 ) ls.push_back( v1->D1 );
        if( abs( v1d2 -> PID ) == 11 || abs( v1d2 -> PID ) == 13 ) ls.push_back( v1->D2 );
        if( abs( v2d1 -> PID ) == 11 || abs( v2d1 -> PID ) == 13 ) ls.push_back( v2->D1 );
        if( abs( v2d2 -> PID ) == 11 || abs( v2d2 -> PID ) == 13 ) ls.push_back( v2->D2 );

        l1 = (GenParticle*) branchGenParticle->At( ls[0] ); l2 = (GenParticle*) branchGenParticle->At( ls[1] );

        if( is_z ) {
            l3 = (GenParticle*) branchGenParticle->At( ls[2] ); l4 = (GenParticle*) branchGenParticle->At( ls[3] );
        }

        // event type

        if ( is_w ) {
            if ( abs( l1->PID ) == 13 && abs( l2->PID ) == 13 ) et = 0;
            else if ( abs( l1->PID ) == 11 && abs( l2->PID ) == 11 ) et = 1;
            else if ( abs( l1->PID ) == 13 && abs( l2->PID ) == 11 ) et = 2;
            else if ( abs( l1->PID ) == 11 && abs( l2->PID ) == 13 ) et = 3;
        }

        if  ( is_z ) {
            if ( abs( l1->PID ) == 13 && abs( l2->PID ) == 13 && abs( l3->PID ) == 13 && abs( l4->PID ) == 13) et = 0;
            else if ( abs( l1->PID ) == 11 && abs( l2->PID ) == 11 && abs( l3->PID ) == 11 && abs( l4->PID ) == 11 ) et = 1;
            else if ( abs( l1->PID ) == 13 && abs( l2->PID ) == 13 && abs( l3->PID ) == 11 && abs( l4->PID ) == 11 ) et = 2;
            else if ( abs( l1->PID ) == 11 && abs( l2->PID ) == 11 && abs( l3->PID ) == 13 && abs( l4->PID ) == 13 ) et = 3;
        }

    }

    if ( debug_bool ) {

        debug_print( debug_bool, " " );

        for ( int i=0; i<vs.size(); i++){

            GenParticle *v = (GenParticle*) branchGenParticle->At( vs[i] );
            debug_print( debug_bool, "parton v" + to_string(i+1) + " pid: " + to_string( v->PID ) + " v m: " +  to_string( v->Mass ) + " v pT: " +  to_string( v->PT ) );

        }

        for ( int i=0; i<ls.size(); i++){

            GenParticle *l = (GenParticle*) branchGenParticle->At( ls[i] );
            debug_print( debug_bool, "parton l" + to_string(i+1) + " pid: " + to_string( l->PID ) + " l pT: " +  to_string( l->PT ) + " l eta: " +  to_string( l->Eta ) + " l phi: " +  to_string( l->Phi ) );

        }

    }

    vvleps = make_pair( vs, ls );

    return make_pair( et, vvleps );

}


pair< TLorentzVector,TLorentzVector > make_parton_vv( string analysis_type = "parton", pair < int, pair < vector < int > , vector < int >>> vvleps = {-99,{{-99},{-99}}}, TClonesArray *branchGenParticle = nullptr ) {

    TLorentzVector v1 = ((GenParticle*) branchGenParticle->At( vvleps.second.first[0] )) -> P4();
    TLorentzVector v2 = ((GenParticle*) branchGenParticle->At( vvleps.second.first[1] )) -> P4();
    
    return make_pair( v1, v2 );

}

pair< TLorentzVector,TLorentzVector > make_parton_2l( string analysis_type = "parton", pair < int, pair < vector < int > , vector < int >>> vvleps = {-99,{{-99},{-99}}}, TClonesArray *branchGenParticle = nullptr ) {

    TLorentzVector l1 = ((GenParticle*) branchGenParticle->At( vvleps.second.second[0] )) -> P4();
    TLorentzVector l2 = ((GenParticle*) branchGenParticle->At( vvleps.second.second[1] )) -> P4();
    
    return make_pair( l1, l2 );

}

pair< pair< TLorentzVector,TLorentzVector >, pair< TLorentzVector,TLorentzVector > > make_parton_4l( string analysis_type = "parton", pair < int, pair < vector < int > , vector < int >>> vvleps = {-99,{{-99},{-99}}}, TClonesArray *branchGenParticle = nullptr ) {

    TLorentzVector l1 = ((GenParticle*) branchGenParticle->At( vvleps.second.second[0] )) -> P4();
    TLorentzVector l2 = ((GenParticle*) branchGenParticle->At( vvleps.second.second[1] )) -> P4();
    TLorentzVector l3 = ((GenParticle*) branchGenParticle->At( vvleps.second.second[2] )) -> P4();
    TLorentzVector l4 = ((GenParticle*) branchGenParticle->At( vvleps.second.second[3] )) -> P4();
    
    pair< TLorentzVector,TLorentzVector > l1l2 = make_pair( l1, l2 );
    pair< TLorentzVector,TLorentzVector > l3l4 = make_pair( l3, l4 );

    return make_pair( l1l2, l3l4 );

}

#endif