#ifndef HELPERFUNCTIONS_H
#define HELPERFUNCTIONS_H

#include <iostream>
#include <vector>
#include <cstdlib>
#include <ctime>
#include <algorithm>
#include <random>
#include <cmath>


void event_dump ( string analysis_type = "particle", TClonesArray *branchGenParticle = nullptr, TClonesArray *branchElectron = nullptr, TClonesArray *branchMuon = nullptr, bool dumpbool = false ) { 
  
    if ( dumpbool ) {

        for(int i=0; i<(int)branchGenParticle->GetEntries(); i++) {

            GenParticle *particle = (GenParticle*) branchGenParticle->At( i );

            // if ( particle->Status != 23 ) continue;

            cout << " particle index " << i << "  - pid " << particle->PID << " , status " << particle->Status << " , pT " << particle->PT << " , eta " << particle->Eta << endl;

            if ( particle->D1 != -1 && particle->D2 != -1 ) {

            GenParticle *d1 = (GenParticle*) branchGenParticle->At( particle->D1 );
            GenParticle *d2 = (GenParticle*) branchGenParticle->At( particle->D2 );

            cout << "       d1 index " << particle->D1 << "  - pid " << d1->PID << " , status " << d1->Status << " , pT " << d1->PT << " , eta " << d1->Eta << endl;
            cout << "       d2 index " << particle->D2 << "  - pid " << d2->PID << " , status " << d2->Status << " , pT " << d2->PT << " , eta " << d2->Eta << endl;

            }

        }

        cout << " " << endl;
    
    }

}


void debug_print(bool debug_bool, const std::string& message) {

  if (debug_bool) {

      std::cout << message << std::endl;

  }

}


template <typename T>
void concatenate_indices(std::vector<T>& out, const std::vector<T>& in) {
    out.insert(out.end(), in.begin(), in.end());
}

template <typename T>
void max(T a, T b) {
    return (a > b) ? a : b;
}

bool is_int_in_vector(const std::vector<int>& vec, int value) {
    return std::find(vec.begin(), vec.end(), value) != vec.end();
}

bool is_prompt( int particle_index = 0, TClonesArray *branchGenParticle = nullptr ) {

    GenParticle *particle  = (GenParticle*) branchGenParticle->At( particle_index );
    GenParticle *mother    = (GenParticle*) branchGenParticle->At( particle->M1 );

    if( mother->PID == particle->PID ) return is_prompt( particle->M1, branchGenParticle );

    if( abs(mother->PID) == 23 || abs(mother->PID) == 24 ) return true; // w or z
    else return false;

}


vector<int> rm_jetlep_overlap( string analysis_type = "reco", vector<int> jets = {-99}, vector<int> leps = {-99}, TClonesArray *branchJet = nullptr, TClonesArray *branchGenParticle = nullptr) {

    // preference given to leptons

    vector<int> skimmed_jets;

    for (int j = 0; j < jets.size(); ++j) {

        Jet* jet = (Jet*) branchJet->At(jets[j]);

        bool overlap = false;

        for (int i = 0; i < leps.size(); ++i) {
            
            GenParticle* genlep = (GenParticle*) branchGenParticle->At(leps[i]);

            if (sqrt(pow(jet->Eta - genlep->Eta, 2) + pow(jet->Phi - genlep->Phi, 2)) < 0.4) {
                overlap = true;
                break;
            }

        }

        if (!overlap) {
            skimmed_jets.push_back(jets[j]);
        }

    }

    return skimmed_jets;

}



double smear_pT( int index = 0, TClonesArray *branchGenParticle = nullptr ) {

    GenParticle* lep = (GenParticle*) branchGenParticle->At( index );

    double pt = lep->PT; double eta = abs(lep->Eta);  double resolution = 0.0;

    if ( abs( lep->PID) == 11 ) {

        if (eta <= 0.5) {
            resolution = std::sqrt(0.03 * 0.03 + pt * pt * 1.3e-3 * 1.3e-3);
        } else if (0.5 < eta && eta <= 1.5) {
            resolution = std::sqrt(0.05 * 0.05 + pt * pt * 1.7e-3 * 1.7e-3);
        } else if (1.5 < eta && eta <= 2.5) {
            resolution = std::sqrt(0.15 * 0.15 + pt * pt * 3.1e-3 * 3.1e-3);
        }

    } else if ( abs( lep->PID) == 13 ) {

        if (eta <= 0.5) {
            resolution = std::sqrt(0.01 * 0.01 + pt * pt * 1.0e-4 * 1.0e-4);
        } else if (0.5 < eta && eta <= 1.5) {
            resolution = std::sqrt(0.015 * 0.015 + pt * pt * 1.5e-4 * 1.5e-4);
        } else if (1.5 < eta && eta <= 2.5) {
            resolution = std::sqrt(0.025 * 0.025 + pt * pt * 3.5e-4 * 3.5e-4);
        }

    }
    
    std::default_random_engine generator(std::random_device{}());
    double sample = std::normal_distribution<>(0, 1)(generator);
    
    return pt + sample * resolution;

}


vector <int> apply_efficiency( vector <int> leps_particle = {0}, TClonesArray *branchGenParticle = nullptr, bool debug_bool = false ) {

    vector<int> leps_reco;

    for (int i = 0; i < leps_particle.size(); ++i) {

        GenParticle* lep = (GenParticle*) branchGenParticle->At( leps_particle[i] );

        double pt = lep->PT; double eta = abs(lep->Eta); double efficiency = 0.0;

        if ( abs( lep->PID) == 11 ) {

            if (pt <= 0.1) {
                efficiency = 0.00;
            } else if (eta <= 1.5) {
                if (pt > 0.1 && pt <= 1.0) {
                    efficiency = 0.73;
                } else if (pt > 1.0 && pt <= 100.0) {
                    efficiency = 0.95;
                } else if (pt > 100.0) {
                    efficiency = 0.99;
                }
            } else if (1.5 <= eta && eta <= 2.5) {
                if (pt > 0.1 && pt <= 1.0) {
                    efficiency = 0.50;
                } else if (pt > 1.0 && pt <= 100.0) {
                    efficiency = 0.83;
                } else if (pt > 100.0) {
                    efficiency = 0.90;
                }
            } else {
                efficiency = 0.00;
            }

        } else if ( abs( lep->PID) == 13 ) {

            if (pt <= 0.1) {
                efficiency = 0.00;
            } else if (eta <= 1.5) {
                if (pt > 0.1 && pt <= 1.0) {
                    efficiency = 0.75;
                } else if (pt > 1.0 && pt <= 1000.0) {
                    efficiency = 0.99;
                } else if (pt > 1000.0) {
                    efficiency = 0.99 * std::exp(0.5 - pt * 5.0e-4);
                }
            } else if (1.5 <= eta && eta <= 2.5) {
                if (pt > 0.1 && pt <= 1.0) {
                    efficiency = 0.70;
                } else if (pt > 1.0 && pt <= 1000.0) {
                    efficiency = 0.98;
                } else if (pt > 1000.0) {
                    efficiency = 0.98 * std::exp(0.5 - pt * 5.0e-4);
                }
            } else {
                efficiency = 0.00;
            }

        }

        double random_value = static_cast<double>(std::rand()) / RAND_MAX;

        if (random_value < efficiency) {
            
            leps_reco.push_back(leps_particle[i]);
            if ( abs(lep->PID) == 11 && lep->Charge == -1) debug_print( debug_bool, "reco e- : pt = " + to_string( smear_pT( leps_particle[i], branchGenParticle ) ) + " , eta = " + to_string(lep->Eta) + " , phi = " + to_string(lep->Phi) );
            if ( abs(lep->PID) == 11 && lep->Charge == 1) debug_print( debug_bool, "reco e+ : pt = " + to_string( smear_pT( leps_particle[i], branchGenParticle ) ) + " , eta = " + to_string(lep->Eta) + " , phi = " + to_string(lep->Phi) );
            if ( abs(lep->PID) == 13 && lep->Charge == -1) debug_print( debug_bool, "reco mu- : pt = " + to_string( smear_pT( leps_particle[i], branchGenParticle ) ) + " , eta = " + to_string(lep->Eta) + " , phi = " + to_string(lep->Phi) );
            if ( abs(lep->PID) == 13 && lep->Charge == 1) debug_print( debug_bool, "reco mu+ : pt = " + to_string( smear_pT( leps_particle[i], branchGenParticle ) ) + " , eta = " + to_string(lep->Eta) + " , phi = " + to_string(lep->Phi) );

        }

    }

    return leps_reco;
}

#endif