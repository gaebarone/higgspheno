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

bool is_parent( int particle_index, TClonesArray *branchGenParticle, int parent_PID) {

    GenParticle *particle  = (GenParticle*) branchGenParticle->At( particle_index );
    GenParticle *mother    = (GenParticle*) branchGenParticle->At( particle->M1 );

    if( mother->PID == particle->PID ) return is_parent( particle->M1, branchGenParticle, parent_PID );

    if( abs(mother->PID) == parent_PID ) return true;
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


double smear_pT_e( double pt, double eta ) {


    double resolution = 0.0;

    if (eta <= 0.5) {
        resolution = std::sqrt(0.03 * 0.03 + pt * pt * 1.3e-3 * 1.3e-3);
    } else if (0.5 < eta && eta <= 1.5) {
        resolution = std::sqrt(0.05 * 0.05 + pt * pt * 1.7e-3 * 1.7e-3);
    } else if (1.5 < eta && eta <= 2.5) {
        resolution = std::sqrt(0.15 * 0.15 + pt * pt * 3.1e-3 * 3.1e-3);
    }
    
    std::default_random_engine generator(std::random_device{}());
    double sample = std::normal_distribution<>(0, 1)(generator);
    
    return pt + sample * resolution;

}

double smear_pT_mu( double pt, double eta ) {


    double resolution = 0.0;

    if (eta <= 0.5) {
        resolution = std::sqrt(0.01 * 0.01 + pt * pt * 1.0e-4 * 1.0e-4);
    } else if (0.5 < eta && eta <= 1.5) {
        resolution = std::sqrt(0.015 * 0.015 + pt * pt * 1.5e-4 * 1.5e-4);
    } else if (1.5 < eta && eta <= 2.5) {
        resolution = std::sqrt(0.025 * 0.025 + pt * pt * 3.5e-4 * 3.5e-4);
    }
    
    std::default_random_engine generator(std::random_device{}());
    double sample = std::normal_distribution<>(0, 1)(generator);
    
    return pt + sample * resolution;

}


double smear_pT( double pt, double eta, int et, int lep_number, string boson ) { // this is for TLorentzVectors

    double pT; int pid = 0;

    if ( boson == "w" ) {

        if ( et == 0 ) pid = 13;

        if ( et == 1 ) pid = 11;

        if ( et == 2 && lep_number == 1 ) pid = 13; 
        if ( et == 2 && lep_number == 2 ) pid = 11;

        if ( et == 3 && lep_number == 1 ) pid = 11; 
        if ( et == 3 && lep_number == 2 ) pid = 13;

    }

    if ( boson == "z" ) {

        if ( et == 0 ) pid = 13;

        if ( et == 1 ) pid = 11;

        if ( et == 2 && lep_number == 1 ) pid = 13; 
        if ( et == 2 && lep_number == 2 ) pid = 13;
        if ( et == 2 && lep_number == 3 ) pid = 11; 
        if ( et == 2 && lep_number == 4 ) pid = 11;

        if ( et == 3 && lep_number == 1 ) pid = 11; 
        if ( et == 3 && lep_number == 2 ) pid = 11;
        if ( et == 3 && lep_number == 3 ) pid = 13; 
        if ( et == 3 && lep_number == 4 ) pid = 13;

    }

    if ( pid == 11 ) pT = smear_pT_e( pt, eta );
    if ( pid == 13 ) pT = smear_pT_mu( pt, eta );

    return pT;

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
            if ( abs(lep->PID) == 11 && lep->Charge == -1) debug_print( debug_bool, "reco e- : pt = " + to_string( smear_pT_e( lep->PT, lep->Eta ) ) + " , eta = " + to_string(lep->Eta) + " , phi = " + to_string(lep->Phi) );
            if ( abs(lep->PID) == 11 && lep->Charge == 1) debug_print( debug_bool, "reco e+ : pt = " + to_string( smear_pT_e( lep->PT, lep->Eta ) ) + " , eta = " + to_string(lep->Eta) + " , phi = " + to_string(lep->Phi) );
            if ( abs(lep->PID) == 13 && lep->Charge == -1) debug_print( debug_bool, "reco mu- : pt = " + to_string( smear_pT_mu( lep->PT, lep->Eta ) ) + " , eta = " + to_string(lep->Eta) + " , phi = " + to_string(lep->Phi) );
            if ( abs(lep->PID) == 13 && lep->Charge == 1) debug_print( debug_bool, "reco mu+ : pt = " + to_string( smear_pT_mu( lep->PT, lep->Eta ) ) + " , eta = " + to_string(lep->Eta) + " , phi = " + to_string(lep->Phi) );

        }

    }

    return leps_reco;
}



std::vector<double> get_probabilities(std::vector<double> &log_probabilities) {
    double sum = 0;

    std::vector<double> probabilities (log_probabilities.size(), 0);

    for (std::vector<double>::const_iterator i = log_probabilities.begin(); i != log_probabilities.end(); ++i) {
        sum += exp(*i);
    }

    for (unsigned int i = 0; i < log_probabilities.size(); i++) {
        probabilities[i] = exp(log_probabilities[i]) / sum;
    }

    return probabilities;
}



void make_csv( vector < string > column_names ) {

    std::ofstream csvFile("bdt/vvhjj.csv");

    if (csvFile.is_open()) {

        for (size_t i = 0; i < column_names.size(); ++i) {

            csvFile << column_names[i];
            if (i < column_names.size() - 1) csvFile << ",";
        
        }
 
        csvFile << "\n";
        csvFile.close();

        std::cout << "CSV file created successfully " << std::endl;

    } else std::cerr << "Error: Could not create the file " << std::endl;

}


void fill_csv( string process_name, string selection, double cross_section, double weight, TLorentzVector b1_reco, TLorentzVector b2_reco, TLorentzVector j1_reco, TLorentzVector j2_reco, TLorentzVector l1_reco, TLorentzVector l2_reco, TLorentzVector l3_reco, TLorentzVector l4_reco, TLorentzVector w1_reco, TLorentzVector w2_reco, TLorentzVector z1_reco, TLorentzVector z2_reco, TLorentzVector met) {
       
    std::ofstream csv_file("bdt/vvhjj.csv", std::ios::app);
    
    if (!csv_file.is_open()) {
        std::cerr << "Error: Unable to open the CSV file!" << std::endl;
        return;
    }

   csv_file 


        << process_name << "," 
        << selection << ","
        << cross_section << "," 
        << weight << "," 


        << b1_reco.Pt() << "," 
        << b1_reco.Eta() << "," 
        << b1_reco.Phi() << "," 

        << b2_reco.Pt() << "," 
        << b2_reco.Eta() << "," 
        << b2_reco.Phi() << "," 

        << (b1_reco+b2_reco).M() << "," 
        << (b1_reco+b2_reco).Pt() << "," 
        << delta_phi(b1_reco,b2_reco) << ","
        << delta_eta(b1_reco,b2_reco) << ","
        << delta_r(b1_reco,b2_reco) << ","


        << j1_reco.Pt() << "," 
        << j1_reco.Eta() << "," 
        << j1_reco.Phi() << "," 

        << j2_reco.Pt() << "," 
        << j2_reco.Eta() << "," 
        << j2_reco.Phi() << ","

        << (j1_reco+j2_reco).M() << "," 
        << (j1_reco+j2_reco).Pt() << "," 
        << delta_phi(j1_reco,j2_reco) << ","
        << delta_eta(j1_reco,j2_reco) << ","
        << delta_r(j1_reco,j2_reco) << ","


        << l1_reco.Pt() << "," 
        << l1_reco.Eta() << "," 
        << l1_reco.Phi() << "," 

        << l2_reco.Pt() << "," 
        << l2_reco.Eta() << "," 
        << l2_reco.Phi() << ","

        << (l1_reco+l2_reco).M() << "," 
        << (l1_reco+l2_reco).Pt() << "," 
        << delta_phi(l1_reco,l2_reco) << ","
        << delta_eta(l1_reco,l2_reco) << ","
        << delta_r(l1_reco,l2_reco) << ","

        << l3_reco.Pt() << "," 
        << l3_reco.Eta() << "," 
        << l3_reco.Phi() << "," 

        << l4_reco.Pt() << "," 
        << l4_reco.Eta() << "," 
        << l4_reco.Phi() << ","

        << (l3_reco+l4_reco).M() << "," 
        << (l3_reco+l4_reco).Pt() << "," 
        << delta_phi(l3_reco,l4_reco) << ","
        << delta_eta(l3_reco,l4_reco) << ","
        << delta_r(l3_reco,l4_reco) << ","


        << w1_reco.Pt() << "," 
        << w1_reco.Eta() << "," 
        << w1_reco.Phi() << "," 

        << w2_reco.Pt() << "," 
        << w2_reco.Eta() << "," 
        << w2_reco.Phi() << ","

        << (w1_reco+w2_reco).M() << "," 
        << (w1_reco+w2_reco).Pt() << "," 
        << delta_phi(w1_reco,w2_reco) << ","
        << delta_eta(w1_reco,w2_reco) << ","
        << delta_r(w1_reco,w2_reco) << ","


        << z1_reco.Pt() << "," 
        << z1_reco.Eta() << "," 
        << z1_reco.Phi() << "," 

        << z2_reco.Pt() << "," 
        << z2_reco.Eta() << "," 
        << z2_reco.Phi() << ","

        << (z1_reco+z2_reco).M() << "," 
        << (z1_reco+z2_reco).Pt() << "," 
        << delta_phi(z1_reco,z2_reco) << ","
        << delta_eta(z1_reco,z2_reco) << ","
        << delta_r(z1_reco,z2_reco) << ","


        << met.Pt() << "\n";


    csv_file.close();

}


#endif