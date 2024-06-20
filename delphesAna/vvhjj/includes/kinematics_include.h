#ifndef KINEMATICS_INCLUDE_H
#define KINEMATICS_INCLUDE_H


#include <cmath>

//------------------------------------------------------------------------------------------------------------------------------------------------------------
// kinematics
//------------------------------------------------------------------------------------------------------------------------------------------------------------


double deltaPhi(TLorentzVector &lorentzvector1, TLorentzVector &lorentzvector2) {

    // return (lorentzvector1.Phi() > lorentzvector2.Phi() ? -1:+1)*TMath::Abs(lorentzvector2.Phi() - lorentzvector1.Phi());
    return (lorentzvector1.Phi() - lorentzvector2.Phi());
}

double deltaEta(TLorentzVector &lorentzvector1, TLorentzVector &lorentzvector2) {

    return (lorentzvector1.Eta() - lorentzvector2.Eta());

}

double deltaR(TLorentzVector &lorentzvector1, TLorentzVector &lorentzvector2) {

    double dPhi = lorentzvector1.Phi() - lorentzvector2.Phi();
    double dEta = lorentzvector1.Eta() - lorentzvector2.Eta();
    return std::sqrt(dPhi*dPhi + dEta*dEta);
    
}

double massTransverse(TLorentzVector &lorentzvector1, TLorentzVector &lorentzvector2) {


    double pT_lepton = lorentzvector1.Pt();
    double pT_miss = lorentzvector2.Pt();
    double dPhi = deltaPhi(lorentzvector1, lorentzvector2);
    return std::sqrt(2 * pT_lepton * pT_miss * (1 - std::cos(dPhi)));

}

//------------------------------------------------------------------------------------------------------------------------------------------------------------
// sorting
//------------------------------------------------------------------------------------------------------------------------------------------------------------

// sort a vector of same type
void sort_by_pT( string analysis_type = "reco", vector<int> indices = vector<int>(), TClonesArray *branchElectron = nullptr, string object = "electron"){

    if ( analysis_type == "reco" && object == "electron" ) {

        sort( indices.begin(), indices.end(), [ branchElectron ] ( const int& lhs, const int& rhs ) {
        return ( (Electron*) branchElectron -> At( lhs )) -> PT > ( (Electron*) branchElectron -> At( rhs )) -> PT;
        });

    } else if ( analysis_type == "reco" && object == "muon" ) {

        sort( indices.begin(), indices.end(), [ branchElectron ] ( const int& lhs, const int& rhs ) {
        return ( (Muon*) branchElectron -> At( lhs )) -> PT > ( (Muon*) branchElectron -> At( rhs )) -> PT;
        });

    } else if ( analysis_type == "reco" && object == "jet" ) {

        sort( indices.begin(), indices.end(), [ branchElectron ] ( const int& lhs, const int& rhs ) {
        return ( (Jet*) branchElectron -> At( lhs )) -> PT > ( (Jet*) branchElectron -> At( rhs )) -> PT;
        });

    } else if ( analysis_type == "particle" ) {

        sort( indices.begin(), indices.end(), [ branchElectron ] ( const int& lhs, const int& rhs ) {
        return ( (GenParticle*) branchElectron -> At( lhs )) -> PT > ( (GenParticle*) branchElectron -> At( rhs )) -> PT;
        });

    } else throw std::runtime_error( "unknown branch" );

}

// get leading and subleading given any two leps
void get_leading_subleading(string analysis_type = "reco", vector<int> leps = {-99}, TClonesArray *branchElectron = nullptr, TClonesArray *branchMuon = nullptr, string branch1 = "electron" , string branch2 = "muon"){

  double lep1_pT = 0;
  double lep2_pT = 0;

    if ( analysis_type == "reco" ) {

        if (branch1 == "electron") {
            lep1_pT = ((Electron*) branchElectron->At(leps[0]))->PT;
        } else if (branch1 == "muon") {
            lep1_pT = ((Muon*) branchMuon->At(leps[0]))->PT;
        }

        if (branch2 == "electron") {
            lep2_pT = ((Electron*) branchElectron->At(leps[1]))->PT;
        } else if (branch2 == "muon") {
            lep2_pT = ((Muon*) branchMuon->At(leps[1]))->PT;
        }

    } else if ( analysis_type == "particle" ) {

        if (branch1 == "electron") {
            lep1_pT = ((GenParticle*) branchElectron->At(leps[0]))->PT;
        } else if (branch1 == "muon") {
            lep1_pT = ((GenParticle*) branchMuon->At(leps[0]))->PT;
        }

        if (branch2 == "electron") {
            lep2_pT = ((GenParticle*) branchElectron->At(leps[1]))->PT;
        } else if (branch2 == "muon") {
            lep2_pT = ((GenParticle*) branchMuon->At(leps[1]))->PT;
        }

    }


    if (lep1_pT > lep2_pT) return;
    else swap(leps[0], leps[1]);

}



#endif