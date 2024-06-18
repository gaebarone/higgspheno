#ifndef HADRONIC_H
#define HADRONIC_H

#include includes/kinematics_include.h

//------------------------------------------------------------------------------------------------------------------------------------------------------------
// jets
//------------------------------------------------------------------------------------------------------------------------------------------------------------


vector <int> get_all_jets( TClonesArray *branchJet = nullptr, double j_pT_min, double j_eta_min ) {

    vector <int> all_jet_indices;

    for( int i = 0; i<(int)branchJet->GetEntries(); i++ ) {

        Jet *jet = (Jet*) branchJet->At(i);

        if( jet->PT > j_pT_min && jet->Eta > j_eta_min)  all_jet_indices.push_back(i); 

    }

    return all_jet_indices;

}


//------------------------------------------------------------------------------------------------------------------------------------------------------------
// higgs
//------------------------------------------------------------------------------------------------------------------------------------------------------------

// returns higgs
pair< TLorentzVector, TLorentzVector> get_higgs( TClonesArray *branchJet = nullptr, TClonesArray *branchGenParticle = nullptr, TClonesArray *branchPFCand = nullptr ){

    TLorentzVector b1, b2;

    std::vector<std::pair< std::map<TString, float>, std::map<TString, std::vector<float>>>>  paired_jets = paired::PAIReDjointEvent( branchGenParticle, branchPFCand, branchJet, 0.4, false, false, true, 1.0, false );
    std::vector<std::pair< std::map<TString, float>, std::map<TString, std::vector<float>>>> paired_b_jets;

    for(int i=0; i<(int)paired_jets.size(); i++){

      std::pair< std::map<TString, float>, std::map<TString, std::vector<float>>> paired_jet = paired_jets.at(i);

      if( paired_jet.first["isbtagged"] > 0 ) paired_b_jets.push_back(paired_jet);

    }

    vector <int> b_jets;
    std::map<TString, float> b_jet;

    b_jet = paired_b_jets.at(0).first;

    b_jets.push_back(b_jet["jet1_index"]);
    b_jets.push_back(b_jet["jet2_index"]);

    b1.SetPtEtaPhiM(b_jet["jet1_pt"],b_jet["jet1_eta"],b_jet["jet1_phi"],b_jet["jet1_mass"]);
    b2.SetPtEtaPhiM(b_jet["jet2_pt"],b_jet["jet2_eta"],b_jet["jet2_phi"],b_jet["jet2_mass"]);

    pair< TLorentzVector, TLorentzVector> bb = < b1, b2 >;

    return bb;
}


//------------------------------------------------------------------------------------------------------------------------------------------------------------
// vbf jets
//------------------------------------------------------------------------------------------------------------------------------------------------------------


// returns pair containing vbf jets
pair< TLorentzVector, TLorentzVector> get_vbf_jets( TClonesArray *branchJet=nullptr){
    
    double j_pT_min = 20.0;
    double j_eta_min = 0.0;

    vector <int> vbf_jet_indices;

    for( int i = 0; i<(int)branchJet->GetEntries(); i++ ) {

        Jet *jet = (Jet*) branchJet->At(i);

        if( jet->PT < j_pT_min && jet->Eta < j_eta_min) continue;

        // if( jet->BTag>0) { // delphes b tagging
        if( isMyBTag( jet, branchGenParticle, 0, 0.4, btagEff, fakeEff) && abs(jet->Eta) < 2.5 ) {  // pseudo b tagging
        
            b_jets.push_back(i); 

        } else {

            not_b_jets.push_back(i);

      }

    }

    Jet *jet1 = (Jet*) branchJet->At( not_b_jets[0] );
    Jet *jet2 = (Jet*) branchJet->At( not_b_jets[1] );
    vbf_jet1 = jet1 -> P4();
    vbf_jet2 = jet1  -> P4();

    if ( deltaEta(vbf_jet1, vbf_jet2) =< 2.5 ) continue;
    else {

        pair< TLorentzVector, TLorentzVector> vbf_jets = < vbf_jet1, vbf_jet2 >;

    }

    return vbf_jets; 

}