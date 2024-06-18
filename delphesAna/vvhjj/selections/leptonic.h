#ifndef LEPTONIC_H
#define LEPTONIC_H

#include includes/kinematics_include.h


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

            if( lep_particle->Status !=1 ) continue;

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

        if ( analysis_type == "reco"  && lep_flavor == "electron" ) sort_by_pT( "reco", lep_indices, branchElectron, "electron" );
        if ( analysis_type == "reco"  && lep_flavor == "muon" ) sort_by_pT( "reco", lep_indices, branchElectron, "muon" );
        if ( analysis_type == "particle" && lep_flavor == "electron" ) sort_by_pT( "particle", lep_indices, branchElectron, "electron" );
        if ( analysis_type == "particle" && ep_flavor == "muon" ) sort_by_pT( "particle", lep_indices, branchElectron, "muon" );

    return lep_indices; 
  
}




//------------------------------------------------------------------------------------------------------------------------------------------------------------
// w
//------------------------------------------------------------------------------------------------------------------------------------------------------------


// returns w = < et, <w1, w2>>
pair< int, pair< TLorentzVector, TLorentzVector >> get_w_leptonic( string analysis_type = "reco", TClonesArray *branchElectron = nullptr, TClonesArray *branchMuon = nullptr ) {

    int event_type_reco = -1;

    TLorentzVector l1, l2;

    // define w selection criteria
    double e_pT_min = 15.0;
    double e_eta_max = 2.5;
    double mu_pT_min = 15.0;
    double mu_eta_max = 2.5;

    // get e+ e- mu+ mu-
    vector <int> e_min_indices_reco = get_leptons( "reco", branchElectron, e_pT_min, e_eta_max, "electron", -1 );
    vector <int> e_plus_indices_reco = get_leptons( "reco", branchElectron, e_pT_min, e_eta_max, "electron", 1 );
    vector <int> mu_min_indices_reco = get_leptons( "reco", branchMuon, mu_pT_min, mu_eta_max, "muon", -1 );
    vector <int> mu_plus_indices_reco = get_leptons( "reco", branchMuon, mu_pT_min, mu_eta_max, "muon", 1 );
 
    vector <int> w_leps;

    if(  ( e_min_indices_reco.size() + e_min_indices_reco.size() + mu_min_indices_reco.size() + mu_plus_indices_reco.size() ) >= 2 ) {

        if ( mu_min_indices_reco.size() > 0 && mu_plus_indices_reco.size() > 0 ) { // case mu- mu+

            event_type_reco = 0;

            w_leps.push_back( mu_min_indices_reco[0] );
            w_leps.push_back( mu_plus_indices_reco[0] );

            get_leading_subleading(w_leps, branchMuon, branchMuon, "muon", "muon");

            l1 = ( (Muon *) branchMuon->At( w_leps[0] ) ) -> P4();
            l2 = ( (Muon *) branchMuon->At( w_leps[1] ) ) -> P4();

        } else if ( e_min_indices_reco.size() > 0 && e_plus_indices_reco.size() > 0 ) {  // case e- e+ 

            event_type_reco = 1; 

            w_leps.push_back( e_min_indices_reco[0] );
            w_leps.push_back( e_plus_indices_reco[0] );

            get_leading_subleading(w_leps, branchElectron, branchElectron, "electron", "electron");

            l1 = ( (Electron *) branchElectron->At( w_leps[0] ) ) -> P4();
            l2 = ( (Electron *) branchElectron->At( w_leps[1] ) ) -> P4();

        } else if ( mu_min_indices_reco.size() > 0 && e_plus_indices_reco.size() > 0 ) { // case mu- e+

            event_type_reco = 2;

            w_leps.push_back( mu_min_indices_reco[0] );
            w_leps.push_back( e_plus_indices_reco[0] );

            get_leading_subleading(w_leps, branchMuon, branchElectron, "muon", "electron");

            l1 = ( (Muon *) branchMuon->At( w_leps[0] ) ) -> P4();
            l2 = ( (Electron *) branchElectron->At( w_leps[1] ) ) -> P4();

        } else if ( e_min_indices_reco.size() > 0 && mu_plus_indices_reco.size() > 0 ) {  // case e- mu+

            event_type_reco = 3; 

            w_leps.push_back( e_min_indices_reco[0] );
            w_leps.push_back( mu_plus_indices_reco[0] );

            sort2_by_pT(w_leps, branchElectron, branchMuon, "Electron", "Muon");

            l1 = ( (Electron *) branchElectron->At( w_leps[0] ) ) -> P4();
            l2 = ( (Muon *) branchMuon->At( w_leps[1] ) ) -> P4();

        }

    }

    met = ( (MissingET*) branchMissingET -> At(0) ) -> P4();

    TLorentzVector w1, w2;

    w1 = l1 + met;
    w2 = l2 + met;

    return < event_type_reco, < w1, w2 > >

}


//------------------------------------------------------------------------------------------------------------------------------------------------------------
// z
//------------------------------------------------------------------------------------------------------------------------------------------------------------


// returns z = < et, <z1, z2>>
pair< int, pair< TLorentzVector, TLorentzVector >> get_z_leptonic( string analysis_type = "reco", TClonesArray *branchElectron = nullptr, TClonesArray *branchMuon = nullptr ) {

    int event_type_reco = -1;

    TLorentzVector l1, l2, l3, l4;

    // define z selection criteria
    double e_pT_min = 15.0;
    double e_eta_max = 2.5;
    double mu_pT_min = 15.0;
    double mu_eta_max = 2.5;

    // get e+ e- mu+ mu-
    vector <int> e_min_indices_reco = get_leptons( "reco", branchElectron, e_pT_min, e_eta_max, "electron", -1 );
    vector <int> e_plus_indices_reco = get_leptons( "reco", branchElectron, e_pT_min, e_eta_max, "electron", 1 );
    vector <int> mu_min_indices_reco = get_leptons( "reco", branchMuon, mu_pT_min, mu_eta_max, "muon", -1 );
    vector <int> mu_plus_indices_reco = get_leptons( "reco", branchMuon, mu_pT_min, mu_eta_max, "muon", 1 );

    vector <int> e_indices_reco;
    vector <int> mu_indices_reco;
    concatenate_indices( e_min_indices_reco, e_indices_reco );
    concatenate_indices( e_plus_indices_reco, e_indices_reco );
    concatenate_indices( mu_min_indices_reco, mu_indices_reco );
    concatenate_indices( mu_plus_indices_reco, mu_indices_reco );
    sort_by_pT( "reco", e_indices_reco, branchElectron, "electron" );
    sort_by_pT( "reco", mu_indices_reco, branchMuon, "muon" );

 
    vector <int> z1_leps;
    vector <int> z2_leps;

    if(  ( e_min_indices_reco.size() + e_min_indices_reco.size() + mu_min_indices_reco.size() + mu_plus_indices_reco.size() ) >= 4 ) {

        if ( mu_min_indices_reco.size() >= 2 && mu_plus_indices_reco.size() >= 2 ) { // case mu- mu+ , mu- mu+

            event_type_reco = 0;

            z1_leps.push_back( mu_min_indices_reco[0] );
            z1_leps.push_back( mu_plus_indices_reco[0] );
            z2_leps.push_back( mu_min_indices_reco[1] );
            z2_leps.push_back( mu_plus_indices_reco[1] );

            l1 = ( (Muon *) branchMuon->At( z1_leps[0] ) ) -> P4();
            l2 = ( (Muon *) branchMuon->At( z1_leps[1] ) ) -> P4();
            l3 = ( (Muon *) branchMuon->At( z2_leps[0] ) ) -> P4();
            l4 = ( (Muon *) branchMuon->At( z2_leps[1] ) ) -> P4();

        } else if ( e_min_indices_reco.size() >= 2 && e_plus_indices_reco.size() >= 2 ) {  // case e- e+ , e- e+

            event_type_reco = 1; 

            z1_leps.push_back( e_min_indices_reco[0] );
            z1_leps.push_back( e_plus_indices_reco[0] );
            z2_leps.push_back( e_min_indices_reco[1] );
            z2_leps.push_back( e_plus_indices_reco[1] );

            l1 = ( (Electron *) branchElectron->At( z1_leps[0] ) ) -> P4();
            l2 = ( (Electron *) branchElectron->At( z1_leps[1] ) ) -> P4();
            l3 = ( (Electron *) branchElectron->At( z2_leps[0] ) ) -> P4();
            l4 = ( (Electron *) branchElectron->At( z2_leps[1] ) ) -> P4();

        } else if ( e_indices_reco.size() >= 2 && mu_indices_reco.size() >= 2 ) { // case mu- mu+ , e- e+ or case e- e+ , mu- mu+

            e1 = ( (Electron *) branchElectron->At( e_indices_reco[0] ) ) -> P4();
            m1 = ( (Muon *) branchMuon->At( mu_indices_reco[0] ) ) -> P4();

            if ( m1.Pt() > e1.Pt()) {

                event_type_reco = 2;

                z1_leps.push_back( mu_min_indices_reco[0] );
                z1_leps.push_back( mu_plus_indices_reco[0] );
                z2_leps.push_back( e_min_indices_reco[0] );
                z2_leps.push_back( e_plus_indices_reco[0] );

                l1 = ( (Muon *) branchMuon->At( z1_leps[0] ) ) -> P4();
                l2 = ( (Muon *) branchMuon->At( z1_leps[1] ) ) -> P4();
                l3 = ( (Electron *) branchElectron->At( z2_leps[0] ) ) -> P4();
                l4 = ( (Electron *) branchElectron->At( z2_leps[1] ) ) -> P4();

            } else if ( e1.Pt() > m1.Pt()) {

                event_type_reco = 3;

                z1_leps.push_back( e_min_indices_reco[0] );
                z1_leps.push_back( e_plus_indices_reco[0] );
                z2_leps.push_back( mu_min_indices_reco[0] );
                z2_leps.push_back( mu_plus_indices_reco[0] );

                l1 = ( (Electron *) branchElectron->At( z1_leps[0] ) ) -> P4();
                l2 = ( (Electron *) branchElectron->At( z1_leps[1] ) ) -> P4();
                l3 = ( (Muon *) branchMuon->At( z2_leps[0] ) ) -> P4();
                l4 = ( (Muon *) branchMuon->At( z2_leps[1] ) ) -> P4();

            }

        }  

    }

    z1 = l1 + l2;
    z2 = l3 + l4;

    return < event_type_reco, < z1, z2 > >

}


#endif
