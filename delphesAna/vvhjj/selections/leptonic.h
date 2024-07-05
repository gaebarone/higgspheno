#ifndef LEPTONIC_H
#define LEPTONIC_H

#include "../includes/kinematics_include.h"
#include "../includes/helperfunctions_include.h"

//------------------------------------------------------------------------------------------------------------------------------------------------------------
// leptons
//------------------------------------------------------------------------------------------------------------------------------------------------------------

// returns a vector sorted by pT of all leptons in event

/*
vector<int> get_all_leptons( string analysis_type = "reco", TClonesArray *branchElectron = nullptr, TClonesArray *branchMuon = nullptr, double pT_min = 15, double eta_max = 2.5, bool debug_bool = false ) {

    vector <int> all_lep_indices;

    // select leptons

    if( branchElectron == nullptr ) return all_lep_indices;
    if( branchMuon == nullptr ) return all_lep_indices;

    if ( analysis_type == "reco" ) {

        for( int i=0; i<(int)branchElectron->GetEntries(); i++ ) {

            Electron *e = (Electron *)branchElectron->At(i);

            // debug_print( debug_bool, analysis_type + " e  " + " : pt = " + to_string(e->PT) + " , eta = " + to_string(e->Eta) + " , phi = " + to_string(e->Phi)); 
            if ( e -> PT > pT_min && fabs( e -> Eta ) < eta_max && e->Charge==-1 ) debug_print( debug_bool, analysis_type + " e-  " + " : pt = " + to_string(e->PT) + " , eta = " + to_string(e->Eta) + " , phi = " + to_string(e->Phi)); 
            if ( e -> PT > pT_min && fabs( e -> Eta ) < eta_max && e->Charge==1 ) debug_print( debug_bool, analysis_type + " e+  " + " : pt = " + to_string(e->PT) + " , eta = " + to_string(e->Eta) + " , phi = " + to_string(e->Phi)); 

            if ( e -> PT > pT_min && fabs( e -> Eta ) < eta_max ) all_lep_indices.push_back(i);

        } 

        for( int i=0; i<(int)branchMuon->GetEntries(); i++ ) {

            Muon *mu = (Muon *)branchMuon->At(i);

            // debug_print( debug_bool, analysis_type + " mu " + " : pt = " + to_string(mu->PT) + " , eta = " + to_string(mu->Eta) + " , phi = " + to_string(mu->Phi) ); 
            if ( mu -> PT > pT_min && fabs( mu -> Eta ) < eta_max && mu->Charge==-1 ) debug_print( debug_bool, analysis_type + " mu- " + " : pt = " + to_string(mu->PT) + " , eta = " + to_string(mu->Eta) + " , phi = " + to_string(mu->Phi) ); 
            if ( mu -> PT > pT_min && fabs( mu -> Eta ) < eta_max && mu->Charge==1 ) debug_print( debug_bool, analysis_type + " mu+ " + " : pt = " + to_string(mu->PT) + " , eta = " + to_string(mu->Eta) + " , phi = " + to_string(mu->Phi) ); 

            if ( mu -> PT > pT_min && fabs( mu -> Eta ) < eta_max ) all_lep_indices.push_back(i);

        } 
        
    } else if ( analysis_type == "particle" ) {

        for( int i=0; i<(int)branchElectron->GetEntries(); i++ ) {

            GenParticle *lep = (GenParticle *)branchElectron->At(i);

            // cout << "index" << i << "pt: " << lep->PT << "eta: " << lep->Eta << "phi: " << lep->Phi << endl;

            if( lep->Status !=1 ) continue;

            // if( abs( lep->PID ) == 11 ) debug_print( debug_bool, analysis_type + " e " + " : pt = " + to_string(lep->PT ) + " , eta = " + to_string(lep->Eta) + " , phi = " + to_string(lep->Phi) + " , status = " + to_string(lep->Status) );
            // if( abs( lep->PID ) == 13 ) debug_print( debug_bool, analysis_type + " mu " + " : pt = " + to_string(lep->PT ) + " , eta = " + to_string(lep->Eta) + " , phi = " + to_string(lep->Phi) + " , status = " + to_string(lep->Status) );
            if ( lep->PID == 11 && (lep -> PT > pT_min) && (fabs( lep -> Eta ) < eta_max) ) debug_print( debug_bool, analysis_type + " e-  " + " : pt = " + to_string(lep->PT ) + " , eta = " + to_string(lep->Eta) + " , phi = " + to_string(lep->Phi) + " , status = " + to_string(lep->Status) );
            if ( lep->PID == -11 && (lep -> PT > pT_min) && (fabs( lep -> Eta ) < eta_max) ) debug_print( debug_bool, analysis_type + " e+  " + " : pt = " + to_string(lep->PT ) + " , eta = " + to_string(lep->Eta) + " , phi = " + to_string(lep->Phi) + " , status = " + to_string(lep->Status) );
            if ( lep->PID == 13 && (lep -> PT > pT_min) && (fabs( lep -> Eta ) < eta_max) ) debug_print( debug_bool, analysis_type + " mu- " + " : pt = " + to_string(lep->PT ) + " , eta = " + to_string(lep->Eta) + " , phi = " + to_string(lep->Phi) + " , status = " + to_string(lep->Status) );
            if ( lep->PID == -13 && (lep -> PT > pT_min) && (fabs( lep -> Eta ) < eta_max) ) debug_print( debug_bool, analysis_type + " mu+ " + " : pt = " + to_string(lep->PT ) + " , eta = " + to_string(lep->Eta) + " , phi = " + to_string(lep->Phi) + " , status = " + to_string(lep->Status) );

            if ( (abs( lep->PID ) == 11 || abs( lep->PID ) == 13) && (lep -> PT > pT_min) && (fabs( lep -> Eta ) < eta_max) ) all_lep_indices.push_back(i);
              
        }

    }

    return all_lep_indices; 
  
}
*/

// returns a vector sorted by pT of all passing leptons of a certain flavour and charge
vector<int> get_leptons( string analysis_type = "reco", TClonesArray *branchElectron = nullptr, double pT_min = 15, double eta_max = 2.5, string lep_flavor = "electron", int charge = -1, bool debug_bool = false ) {

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

            if( lep->Status !=1 ) continue;
            if( !is_prompt( i, branchElectron )) continue;

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

    // sort the indices by pT
    sort_by_pT( analysis_type, lep_indices, branchElectron, lep_flavor );

    return lep_indices; 
  
}




//------------------------------------------------------------------------------------------------------------------------------------------------------------
// w
//------------------------------------------------------------------------------------------------------------------------------------------------------------


// returns w = < et, < w leps >>
pair< int, vector <int>> get_ww_leptonic( string analysis_type = "reco", vector < int > eminus = {0}, vector < int > eplus = {0} , vector < int > muminus = {0} , vector < int > muplus = {0} , TClonesArray *branchElectron = nullptr, TClonesArray *branchMuon = nullptr, TClonesArray *branchMissingET = nullptr ) {

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

/* SORTING BY MASS - STILL WORKING ON THIS

pair< int, vector <int>> get_ww_leptonic( string analysis_type = "reco", TClonesArray *branchGenParticle = nullptr, TClonesArray *branchMissingET = nullptr, string sign = "same") {

    // int event_type = -1;
    int event_type = 99;

    // define w selection criteria
    double pT_min = 15.0;
    double eta_max = 2.5;

    vector <int> w_leps;

    vector<int> leptons = get_all_leptons( "particle", branchGenParticle, branchGenParticle, pT_min, eta_max, false );

    debug_print( debug_bool, " ------ " + analysis_type + " ww selection ------ " );

    if(  leptons.size() >= 2 ) {

        sort_by_w_mass(analysis_type, leptons, branchGenParticle, branchMissingET);

        if ( analysis_type == "reco" ) apply_efficiency( leptons, 0.95 );
        
        w_leps.push_back( leptons[0] );
        w_leps.push_back( leptons[1] );

        TLorentzVector l1 = ( (GenParticle*)branchGenParticle->At( w_leps[0] ) )->P4();
        TLorentzVector l2 = ( (GenParticle*)branchGenParticle->At( w_leps[1] ) )->P4();
        TLorentzVector nu = ( (MissingET*)branchMissingET->At( 0 ) )->P4();


        debug_print( debug_bool, " " );
        debug_print( debug_bool, analysis_type + " event type: " + to_string(event_type) );
        debug_print( debug_bool, analysis_type + " w1: pt " + to_string( l1.Pt() + nu.Pt()/2 ) + " , m = " + to_string( l1.M() + nu.M()/2 ) );
        debug_print( debug_bool, analysis_type + " w1: pt " + to_string( l2.Pt() + nu.Pt()/2 ) + " , m = " + to_string( l2.M() + nu.M()/2 ) );
        debug_print( debug_bool, analysis_type + " l1: pt " + to_string( l1.Pt() ) + " , eta = " + to_string( l1.Eta() ) );
        debug_print( debug_bool, analysis_type + " l2: pt " + to_string( l2.Pt() ) + " , eta = " + to_string( l2.Eta() ) );
        debug_print( debug_bool, " " );

    } else {

        debug_print( debug_bool, " " );
        debug_print( debug_bool, "no ww pair" );
        debug_print( debug_bool, " " );

    }

    return make_pair( event_type, w_leps );

}
*/

pair< vector<TLorentzVector> , vector<TLorentzVector>> make_ww_leptonic( string analysis_type = "reco", pair< int, vector <int>> w_leps = {-99, {-99}}, TClonesArray *branchElectron = nullptr, TClonesArray *branchMuon = nullptr ,  TClonesArray *branchMissingET = nullptr ) {

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

//------------------------------------------------------------------------------------------------------------------------------------------------------------
// z
//------------------------------------------------------------------------------------------------------------------------------------------------------------


// returns z = < et, < z1 leps, z2 leps>>
pair< int, pair< vector <int>, vector <int>>> get_zz_leptonic( string analysis_type = "reco", vector < int > eminus = {0} , vector < int > eplus = {0} , vector < int > muminus = {0} , vector < int > muplus = {0}, TClonesArray *branchElectron = nullptr, TClonesArray *branchMuon = nullptr ) {

    int event_type = -1;

    TLorentzVector l1, l2, l3, l4;

    vector <int> e;
    vector <int> mu;
    concatenate_indices( eminus, e );
    concatenate_indices( eplus, e );
    concatenate_indices( muminus, mu );
    concatenate_indices( muplus, mu );

    sort_by_pT( analysis_type, e, branchElectron, "electron" );
    sort_by_pT( analysis_type, mu, branchMuon, "muon" );
 
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

#endif