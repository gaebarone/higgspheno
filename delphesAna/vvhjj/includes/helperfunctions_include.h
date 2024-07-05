#ifndef HELPERFUNCTIONS_H
#define HELPERFUNCTIONS_H

#include <iostream>
#include <vector>
#include <cstdlib>
#include <ctime>
#include <algorithm>
#include <random>
#include <cmath>


void debug_print(bool debug_bool, const std::string& message) {

  if (debug_bool) {

      std::cout << message << std::endl;

  }

}



template <typename T>
void concatenate_indices(std::vector<T>& out, const std::vector<T>& in) {
    out.insert(out.end(), in.begin(), in.end());
}


bool is_prompt( int particle_index = 0, TClonesArray *branchGenParticle = nullptr ) {

    GenParticle *particle  = (GenParticle*) branchGenParticle->At( particle_index );
    GenParticle *mother    = (GenParticle*) branchGenParticle->At( particle->M1 );

    if( mother->PID == particle->PID ) return is_prompt( particle->M1, branchGenParticle );

    if( abs(mother->PID) == 23 || abs(mother->PID) == 24 ) return true; // w or z
    else return false;

}



void rm_jetlep_overlap( vector<int>& jets_particle, vector<int>& leps_particle, TClonesArray *branchJet = nullptr, TClonesArray *branchElectron = nullptr, TClonesArray *branchMuon = nullptr) {

    // preference give to leptons

    vector<int> skimmed_jets_particle;

    for (int i = 0; i < jets_particle.size(); ++i) {
        bool found_match = false;
        Jet* genjet = (Jet*) branchJet->At( jets_particle[i] );

        for (int j = 0; j < leps_particle.size(); ++j) {

            GenParticle* genlep = (GenParticle*) branchElectron->At( leps_particle[j] );

            if ( sqrt( pow( genjet->Eta - genlep->Eta, 2 ) + pow( genjet->Phi - genlep->Phi, 2 ) ) < 0.4 ) {
                
                found_match = true; 
                break;

            }

        }

        if ( !found_match ) skimmed_jets_particle.push_back(i);

    }

    jets_particle = skimmed_jets_particle;

}   




void apply_efficiency (std::vector<int>& vec, double efficiency) {
    if (efficiency < 0.0 || efficiency > 1.0) {
        std::cerr << "efficiency must be between 0 and 1." << std::endl;
        return;
    }

    int totalItems = vec.size();
    double itemsToDeleteDouble = totalItems * (1 - efficiency);
    int itemsToDelete = std::floor(itemsToDeleteDouble);

    // Introduce randomness to sometimes delete 1 more item
    double fractionalPart = itemsToDeleteDouble - itemsToDelete;
    std::random_device rd;
    std::mt19937 g(rd());
    std::uniform_real_distribution<> dis(0.0, 1.0);
    if (dis(g) < fractionalPart) {
        itemsToDelete++;
    }

    std::vector<int> indices(totalItems);
    for (int i = 0; i < totalItems; ++i) {
        indices[i] = i;
    }

    std::shuffle(indices.begin(), indices.end(), g);

    for (int i = 0; i < itemsToDelete; ++i) {
        int indexToDelete = indices[i];
        vec.erase(vec.begin() + indexToDelete);
        for (int j = i + 1; j < itemsToDelete; ++j) {
            if (indices[j] > indexToDelete) {
                indices[j]--;
            }
        }
    }
}

#endif