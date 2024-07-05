// not needed 

//#ifdef __CLING__
//R__LOAD_LIBRARY("libDelphes")
//#endif

//#define ONNXRUN
#ifdef ONNXRUN
//#include <onnxruntime/core/session/onnxruntime_cxx_api.h>
//#include "core/session/onnxruntime_cxx_api.h"
#include <onnxruntime_cxx_api.h>
#endif


// #define MDEBUG

#define MSEED 1234 

bool debug_bool = false;

#include "../common_includes/trasnform_inputs.h"
#include <unordered_map>
#include "HepMC/GenParticle.h"
#include "classes/DelphesClasses.h"
#include "classes/DelphesLHEFReader.h"
#include "external/ExRootAnalysis/ExRootTreeReader.h"
#include "../common_includes/ghost_tagging.h"
#include "../common_includes/combinations.h"
//#include "../common_includes/get_cross_section.h"
#include "../common_includes/make_paired.h"
#include <iostream>
#include <fstream>
#include <string>
#include <map>
#include "TH1F.h"
#include "TH2F.h"
#include "TClonesArray.h"
#include "TTree.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <TROOT.h>
#include <TLorentzVector.h>
#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>
#include <TMath.h>
#include <Rtypes.h>
#include <TString.h>
#include <TRandom.h>
#include <TRandom3.h>
#include "TParticle.h"
#include "TCanvas.h"
#include "TProfile.h"
#include <vector>
//#include "selections.h"
//#include "parton_selections.h"
#include <iomanip>
#include  <string.h>
#include <cmath>

#include "selections/leptonic.h"
#include "selections/hadronic.h"

#include "includes/cutflow_include.h"
#include "includes/crossx_include.h"
#include "includes/hist_include.h"
#include "includes/weights_include.h"
#include "includes/selections_include.h"
#include "includes/kinematics_include.h"
#include "includes/helperfunctions_include.h"
 
#include "lepAnalyzer.h"

using namespace std;

#ifdef ONNXRUN
using namespace ::Ort;
#endif


//------------------------------------------------------------------------------------------------------------------------------------------------------------
// MISC
//------------------------------------------------------------------------------------------------------------------------------------------------------------


template <typename T>
	T VectorProduct(const std::vector<T>& v){
		return std::accumulate(v.begin(), v.end(), 1, std::multiplies<T>());
	};

std::string print_shape(const std::vector<int64_t>& v){
  std::stringstream ss("");
  for (size_t i = 0; i < v.size() - 1; i++)
    ss << v[i] << "x";
  ss << v[v.size() - 1];
  return ss.str();
}

int calculate_product(const std::vector<int64_t>& v) {
  int total = 1;
  for (auto& i : v) total *= i;
  return total;
}


//------------------------------------------------------------------------------------------------------------------------------------------------------------
// Z ANALYZER
//------------------------------------------------------------------------------------------------------------------------------------------------------------

void printMap(const std::map<std::string, std::pair<int, double>>& map) {
    for (const auto& pair : map) {
        std::cout << " CUT: " << pair.first << " = " << pair.second.first << std::endl;
    }
}

void zAnalyzer(const char *input_file, const char *output_file, const char *process_name, bool debug_bool){

  #ifdef __CLING__
    gSystem->Load("libDelphes");
  #endif

// cutflow table

  vector <TH1F*> listOfTH1;
  vector <TH2F*> listOfTH2;
  vector <TProfile*> listOfTProfiles;

  DefineSelections();
  
  cut_list_reco = cut_sel_process_reco["all"];
  cut_list_particle = cut_sel_process_particle["all"];
  cutList_parton = cut_sel_process_parton["all"];

  std::map< string, bool > enable_cut_reco;  
  std::map< string, bool > enable_cut_particle; 

  std::map<string, bool> enableCutParton;
 
  for(vector<string>::iterator it=cut_sel_process_reco["all"].begin(); it!=cut_sel_process_reco["all"].end(); it++){
    enable_cut_reco[ (*it)] = hasCut(cut_list_reco, (*it));
  }
  for(vector<string>::iterator it=cut_sel_process_particle["all"].begin(); it!=cut_sel_process_particle["all"].end(); it++){
    enable_cut_particle[ (*it)] = hasCut(cut_list_particle, (*it));
  }
  for(vector<string>::iterator it=cut_sel_process_parton["all"].begin(); it!=cut_sel_process_parton["all"].end(); it++){
    enableCutParton[ (*it)] = hasCut(cutList_parton, (*it));
  }

  // new cft

  int cutVal_reco = 0;
  double cutValW_reco = 0;

  std::map<string, std::pair<int,double>> cft_total_reco;
  std::vector<std::map<std::string, std::pair<int, double>>> cft_allevents_reco;

  for(int i=0; i<(int) cut_list_reco.size(); i++) { 
    cft_total_reco[cut_list_reco.at(i)] = make_pair(0,0.0); 
  }
 
  int cutVal_particle = 0;
  double cutValW_particle = 0;

  std::map<string, std::pair<int,double>> cft_total_particle;
  std::vector<std::map<std::string, std::pair<int, double>>> cft_allevents_particle;

  for(int i=0; i<(int) cut_list_particle.size(); i++) { 
    cft_total_particle[cut_list_particle.at(i)] = make_pair(0,0.0); 
  }

  // new cft
  
  int cutVal_parton = 0;
  double cutValW_parton = 0;

  std::map<string, std::pair<int,double>> cutFlowMap_parton;
  for(int i=0; i<(int) cutList_parton.size(); i++) { 
    cutFlowMap_parton[cutList_parton.at(i)] = make_pair(0,0.0); 
  }
 
  vector <string> selType={"reco","particle","parton"};
  std::map<string, vector<string>> cutFlowMByType;
  cutFlowMByType["reco"] = cut_list_reco;
  cutFlowMByType["particle"] = cut_list_particle;
  cutFlowMByType["parton"]=cutList_parton;
 
  std::map<string,TH1F*> cutFlowHists;
  std::map<string,TProfile*> cutFlowEffs;

  typedef std::map<std::string, std::pair<int,double>> cutFlowMapDef;
  std::map<string, cutFlowMapDef* > cutFlowMapAll;
    cutFlowMapAll["reco"] =  & cft_total_reco;
    cutFlowMapAll["particle"] = & cft_total_particle;
    cutFlowMapAll["parton"] = & cutFlowMap_parton;
  
  for(std::vector<string>::iterator it=selType.begin(); it!=selType.end(); it++){
    cutFlowHists[(*it)]=new TH1F(Form("hSel_%s",(*it).c_str()),"",cutFlowMByType[(*it)].size(),0,cutFlowMByType[(*it)].size()+1);
    cutFlowEffs[(*it)]=new TProfile(Form("hEff_%s",(*it).c_str()),"",cutFlowMByType[(*it)].size(),0,cutFlowMByType[(*it)].size()+1);
    
    listOfTH1.push_back(cutFlowHists[(*it)]);
    listOfTProfiles.push_back((cutFlowEffs[(*it)]));
  }
  
// delphes

  TChain chain("Delphes");
  chain.Add(input_file);

  ExRootTreeReader *treeReader = new ExRootTreeReader(&chain);
  Long64_t numberOfEntries = treeReader->GetEntries();
  Long64_t numEntries = get_total_events(process_name);  
  if(numEntries==-1) numEntries=numberOfEntries;
  cout<<"NUMBER OF ENTRIES: "<<numEntries<<endl;
  double cross_section = get_cross_section(process_name);
  cout<<"CROSS SECTION: "<<cross_section<<endl;
  Float_t totalWeight = 0.0;

  TClonesArray *branchJet = treeReader->UseBranch("Jet");
  TClonesArray *branchElectron = treeReader->UseBranch("Electron");
  TClonesArray *branchMuon = treeReader->UseBranch("Muon");
  TClonesArray *branchEvent = treeReader->UseBranch("Event");
  TClonesArray *branchGenParticle = treeReader->UseBranch("Particle");
  TClonesArray *branchGenJet = treeReader->UseBranch("GenJet");
  TClonesArray *branchMissingET = treeReader->UseBranch("MissingET");
  TClonesArray *branchGenMissingET = treeReader->UseBranch("GenMissingET");
  TClonesArray *branchWeight  = treeReader->UseBranch("Weight");

  TClonesArray *branchPFCand = nullptr;

  TBranch *branch = nullptr;

  for (int i = 0; i < chain.GetListOfBranches()->GetEntries(); ++i) {

    branch = dynamic_cast<TBranch*>(chain.GetListOfBranches()->At(i));

    if (strcmp(branch->GetName(), "ParticleFlowCandidate") == 0){

      branchPFCand = treeReader->UseBranch("ParticleFlowCandidate");

    }
  }

  TH1F *hWeight = new TH1F("weights", "weight", 50, 0.0, 1.0);
  listOfTH1.push_back(hWeight);

  TFile *hists= new TFile(output_file,"recreate");


//------------------------------------------------------------------------------------------------------------------------------------------------------------
// initialize histograms
//------------------------------------------------------------------------------------------------------------------------------------------------------------


vector<string> selections = { "HWWJJ", "WWJJ" }; // , "HZZJJ", "(HWW)JJ", "(HZZ)JJ", "WWJJ", "ZZJJ" }; 

map< string, histograms > sel_hist_map;

for (auto selection = selections.begin(); selection != selections.end(); ++selection) {
  
  sel_hist_map[ *selection ] = histograms( *selection );
  sel_hist_map[ *selection ].initialize_bb();
  sel_hist_map[ *selection ].initialize_jj();
  sel_hist_map[ *selection ].initialize_ww();
  sel_hist_map[ *selection ].initialize_zz();

}


//------------------------------------------------------------------------------------------------------------------------------------------------------------
// DEF QUANTITIES
//------------------------------------------------------------------------------------------------------------------------------------------------------------


  double nPassed=0;
  double totWeightedEntries=0;
  int nPassedRaw=0;

  double Lumi=1000*300; // 1000 to convert from pb to fb, 300 for end of run 3
  cout << "LUMI: "<< Lumi << endl;


//------------------------------------------------------------------------------------------------------------------------------------------------------------
// WEIGHTS
//------------------------------------------------------------------------------------------------------------------------------------------------------------


  double sumOfWeights=0;

  TH1F *hClosure = new TH1F("hClosure","hClosure",1,0,1);
  listOfTH1.push_back(hClosure);

  for(Int_t entry = 0; entry < numberOfEntries; ++entry){

    // load branches with data from specified event
    treeReader->ReadEntry(entry);
    HepMCEvent *event = (HepMCEvent*) branchEvent -> At(0);
    totalWeight += event->Weight;

  }

  cout << "TOTAL WEIGHT: "<< totalWeight << endl;
    

//------------------------------------------------------------------------------------------------------------------------------------------------------------
// EVENT LOOP
//------------------------------------------------------------------------------------------------------------------------------------------------------------


  if (debug_bool) numberOfEntries = 25;

  for(Int_t entry = 0; entry < numberOfEntries; ++entry) {

    debug_print(debug_bool, "----------------------------------------------------------------------------------------------------");
    debug_print(debug_bool, "----------------------------------------------------------------------------------------------------");
    debug_print(debug_bool, " ");
    debug_print(debug_bool, " ****** EVENT: " + to_string(entry) + " ****** ");
    debug_print(debug_bool, " ");

    treeReader->ReadEntry(entry);
    HepMCEvent *event = (HepMCEvent*) branchEvent -> At(0);
    Float_t weight = event->Weight*Lumi*cross_section*numberOfEntries/(numEntries*totalWeight);
    Float_t test_weight = event->Weight*cross_section*numberOfEntries/(numEntries*totalWeight);
    hWeight -> Fill(event->Weight, test_weight);


  //------------------------------------------------------------------------------------------------------------------------------------------------------------
  // inititalize event cft
  //------------------------------------------------------------------------------------------------------------------------------------------------------------
   
    std::map<string, std::pair<int,double>> cft_event_reco;
    std::map<string, std::pair<int,double>> cft_event_particle;

    for(int i=0; i<(int) cut_list_reco.size(); i++) { 
      cft_event_reco[cut_list_reco.at(i)] = make_pair(0,0.0); 
    }

    for(int i=0; i<(int) cut_list_particle.size(); i++) { 
      cft_event_particle[cut_list_particle.at(i)] = make_pair(0,0.0); 
    }

    update_cft(cft_event_reco, cft_total_reco, "initial - reco", weight);
    update_cft(cft_event_particle, cft_total_particle, "initial - particle", weight);


  //------------------------------------------------------------------------------------------------------------------------------------------------------------
  // inititalize physics objects
  //------------------------------------------------------------------------------------------------------------------------------------------------------------


  string analysis_type;

  // reco
    bool found_bb_reco = false;
    bool found_jj_reco = false;
    bool found_ww_reco = false;
    bool found_zz_reco = false;

    vector < int > eminus_reco, eplus_reco, muminus_reco, muplus_reco, leps_reco, jets_reco;
    int nleps_reco, njets_reco;

    std::vector<std::pair< std::map<TString, float>, std::map<TString, std::vector<float>>>> bb_reco_i;         // paired b jets (higgs)
    // vector < int > bb_reco_i;                                                                                // pseudo b tagged jets (higgs)
    pair<int, int> jj_reco_i;                                                                                   // pair of vbf jets
    pair< int, vector <int>> ww_reco_i;                                                                         // event type and pair of leading subleading w index 
    pair< int, pair< vector <int>, vector <int>>> zz_reco_i;                                                    // event type and pair of leading subleading z index

    pair< TLorentzVector, TLorentzVector> bb_reco;                    // < b1, b2 >
    pair< TLorentzVector, TLorentzVector> jj_reco;                    // < j1, j2 >
    pair< vector<TLorentzVector>, vector<TLorentzVector>> ww_reco;    // < <leps>, <ws> >
    pair< vector<TLorentzVector>, vector<TLorentzVector>> zz_reco;    // < <leps>, <zs> >

  // particle
    bool found_bb_particle = false;
    bool found_jj_particle = false;
    bool found_ww_particle = false;
    bool found_zz_particle = false;

    vector < int > eminus_particle, eplus_particle, muminus_particle, muplus_particle, leps_particle, jets_particle;
    int nleps_particle, njets_particle;

    std::vector<std::pair< std::map<TString, float>, std::map<TString, std::vector<float>>>> bb_particle_i;      // paired b jets (higgs)
    // vector < int > bb_particle_i;                                                                             // pseudo b tagged jets (higgs)
    pair<int, int> jj_particle_i;                                                                                // pair of vbf jets
    pair< int, vector <int>> ww_particle_i;                                                                      // event type and pair of leading subleading w index
    pair< int, pair< vector <int>, vector <int>>> zz_particle_i;                                                 // event type and pair of leading subleading z index

    pair< TLorentzVector, TLorentzVector> bb_particle;                    // < b1, b2 >
    pair< TLorentzVector, TLorentzVector> jj_particle;                    // < j1, j2 >
    pair< vector<TLorentzVector>, vector<TLorentzVector>> ww_particle;    // < <leps>, <ws> >
    pair< vector<TLorentzVector>, vector<TLorentzVector>> zz_particle;    // < <leps>, <zs> >


  //------------------------------------------------------------------------------------------------------------------------------------------------------------
  // reco selection
  //------------------------------------------------------------------------------------------------------------------------------------------------------------

  analysis_type = "reco";

  debug_print( debug_bool, " " );
  debug_print( debug_bool, "------------------------------------------------------------------------" );
  debug_print( debug_bool, " reco leps and jets - passing pT and eta cut " );
  debug_print( debug_bool, "------------------------------------------------------------------------" );
  debug_print( debug_bool, " " );

    // all leps and jets
    eminus_reco = get_leptons( analysis_type, branchElectron, 15.0 , 2.5, "electron", -1, debug_bool );
    eplus_reco = get_leptons( analysis_type, branchElectron, 15.0 , 2.5, "electron", 1, debug_bool );
    muminus_reco = get_leptons( analysis_type, branchMuon, 15.0 , 2.5, "muon", -1, debug_bool );
    muplus_reco = get_leptons( analysis_type, branchMuon, 15.0 , 2.5, "muon", 1, debug_bool );
    concatenate_indices(leps_reco, eminus_reco); concatenate_indices(leps_reco, eplus_reco); concatenate_indices(leps_reco, muminus_reco); concatenate_indices(leps_reco, muplus_reco);

    jets_reco = get_all_jets( "reco", leps_reco, branchJet, branchElectron, branchMuon, 20 , debug_bool );

    nleps_reco = eminus_reco.size() + eplus_reco.size() + muminus_reco.size() + muplus_reco.size();
    njets_reco = jets_reco.size();

  debug_print( debug_bool, " " );
  debug_print( debug_bool, analysis_type + " num leps: " + to_string( nleps_reco ) + " ( # e- " + to_string( eminus_reco.size() ) + " # e+ " + to_string( eplus_reco.size() ) + " # mu- " + to_string( muminus_reco.size() ) + " # mu+ " + to_string( muplus_reco.size() ) + " ) " );
  debug_print( debug_bool, "number of reco jets: " + to_string(njets_reco) );
  debug_print( debug_bool, " " );

  // get
    jj_reco_i = get_jj( "reco", jets_reco, branchJet );
    if ( njets_reco >= 4 )   bb_reco_i = get_bb( "reco", jj_reco_i, branchJet, branchGenParticle, branchPFCand );

    if ( nleps_reco < 4 )         ww_reco_i = get_ww_leptonic( "reco", eminus_reco, eplus_reco, muminus_reco, muplus_reco, branchElectron, branchMuon, branchMissingET );
    else if ( nleps_reco >= 4 )   zz_reco_i = get_zz_leptonic( "reco", eminus_reco, eplus_reco, muminus_reco, muplus_reco, branchElectron, branchMuon );

  // set
    if ( bb_reco_i.size() >=2 ) found_bb_reco = true;
    if ( jj_reco_i.first != -1 && jj_reco_i.second != -1 ) found_jj_reco = true;
    if ( ww_reco_i.first != -1 && ww_reco_i.second.size() >= 2) found_ww_reco = true;
    if ( zz_reco_i.first != -1 && zz_reco_i.second.first.size() >= 2 && zz_reco_i.second.second.size() >= 2) found_zz_reco = true;


  //------------------------------------------------------------------------------------------------------------------------------------------------------------
  // particle selection
  //------------------------------------------------------------------------------------------------------------------------------------------------------------

    analysis_type = "particle";

  debug_print( debug_bool, " " );
  debug_print( debug_bool, "------------------------------------------------------------------------" );
  debug_print( debug_bool, " particle leps and jets - passing pT and eta cut " );
  debug_print( debug_bool, "------------------------------------------------------------------------" );
  debug_print( debug_bool, " " );

    // all leps and jets
    eminus_particle = get_leptons( analysis_type, branchGenParticle, 15.0 , 2.5, "electron", -1, debug_bool );
    eplus_particle = get_leptons( analysis_type, branchGenParticle, 15.0 , 2.5, "electron", 1, debug_bool );
    muminus_particle = get_leptons( analysis_type, branchGenParticle, 15.0 , 2.5, "muon", -1, debug_bool );
    muplus_particle = get_leptons( analysis_type, branchGenParticle, 15.0 , 2.5, "muon", 1, debug_bool );
    concatenate_indices(leps_particle, eminus_particle); concatenate_indices(leps_particle, eplus_particle); concatenate_indices(leps_particle, muminus_particle); concatenate_indices(leps_particle, muplus_particle);

    jets_particle = get_all_jets( "particle", leps_particle, branchGenJet, branchGenParticle, branchGenParticle, 20 , debug_bool );

    nleps_particle = eminus_particle.size() + eplus_particle.size() + muminus_particle.size() + muplus_particle.size();
    njets_particle = jets_particle.size();

  debug_print( debug_bool, " " );
  debug_print( debug_bool, analysis_type + " num leps: " + to_string( nleps_particle ) + " ( # e- " + to_string( eminus_particle.size() ) + " # e+ " + to_string( eplus_particle.size() ) + " # mu- " + to_string( muminus_particle.size() ) + " # mu+ " + to_string( muplus_particle.size() ) + " ) " );
  debug_print( debug_bool, "number of particle jets: " + to_string(njets_particle) );
  debug_print( debug_bool, " " );

  // get
    jj_particle_i = get_jj( "particle", jets_particle, branchGenJet );
    if ( njets_particle >= 4 ) bb_particle_i = get_bb( "particle", jj_particle_i, branchGenJet, branchGenParticle, branchPFCand );
    
    if ( nleps_particle < 4 )         ww_particle_i = get_ww_leptonic( "particle", eminus_particle, eplus_particle, muminus_particle, muplus_particle, branchGenParticle, branchGenParticle, branchGenMissingET );
    else if ( nleps_particle >= 4 )   zz_particle_i = get_zz_leptonic( "particle", eminus_particle, eplus_particle, muminus_particle, muplus_particle, branchGenParticle, branchGenParticle );

  // set 
    if ( jj_particle_i.first != -1 && jj_particle_i.second != -1 )  found_jj_particle = true;
    if ( bb_particle_i.size() >= 2 )  found_bb_particle = true;
    if ( ww_particle_i.first != -1 && ww_particle_i.second.size() >= 2) found_ww_particle = true;
    if ( zz_particle_i.first != -1 && zz_particle_i.second.first.size() >= 2 && zz_particle_i.second.second.size() >= 2)  found_zz_particle = true;


  //------------------------------------------------------------------------------------------------------------------------------------------------------------
  // parton selection
  //------------------------------------------------------------------------------------------------------------------------------------------------------------
/*
  // higgs
    int switchVal_parton = 0;
  
    if(enableCutParton["initial parton"]){
      increaseCount(cutFlowMap_parton,"initial parton",weight);
    }

    bool HiggsRecord=FillHiggsTruthRecord(branchGenParticle,h_parton,b1_parton,b2_parton,j1_parton,j2_parton);

    if(enableCutParton["Higgs Candidate"]){
      if(switchVal_parton == 0 && HiggsRecord) increaseCount(cutFlowMap_parton,"Higgs Candidate",weight);
      else  switchVal_parton = 1;
    }

    if(HiggsRecord){
      double bbdeltaPhiparton = deltaPhi(b1_parton, b2_parton);
      double bbdeltaEtarparton = deltaEta(b1_parton, b2_parton);
      double bbdeltaRrparton = deltaR(b1_parton, b2_parton);
    }

    // leps
    int thisPartonEventType=-1;

    vector <int> goodE_parton_indices  = GoodElectronPartonIndices(branchGenParticle, analysis);
    vector <int> goodMu_parton_indices = GoodMuonPartonIndices(branchGenParticle, analysis);
    // goodE_size_parton->Fill(goodE_parton_indices.size(),weight);
    // goodMu_size_parton->Fill(goodMu_parton_indices.size(),weight);

    vector <int> ZPartonIndices;
    vector <int> WPartonIndices;
    bool foundZZ = false;
    bool foundWW = false;

    if(analysis == "HZZJJ"){

      ZPartonIndices = GetZPartonIndices(branchGenParticle, analysis);
      
      if(ZPartonIndices.size() > 1) foundZZ = true;

      if(enableCutParton["ZZ parton"]){
        if(switchVal_parton == 0 && foundZZ) increaseCount(cutFlowMap_parton,"ZZ parton",weight);
        else switchVal_parton = 1;
      }
  
      if(switchVal_parton == 0) getPartonZLeps(thisPartonEventType, ZPartonIndices, branchGenParticle, z1_parton, z2_parton, l1_parton, l2_parton, l3_parton, l4_parton, q1_parton, q2_parton, q3_parton, q4_parton);

      if(foundZZ){
        double zzdeltaPhiparton = deltaPhi(z1_parton, z2_parton);
        double zzdeltaEtarparton = deltaEta(z1_parton, z2_parton);
        double zzdeltaRrparton = deltaR(z1_parton, z2_parton);
      }


    } if(analysis == "HWWJJ") { 

      WPartonIndices = GetWPartonIndices(branchGenParticle, analysis);

      if(WPartonIndices.size() > 1) foundWW = true;

      if(enableCutParton["WW parton"]){
        if(switchVal_parton == 0 && foundWW) increaseCount(cutFlowMap_parton,"WW parton",weight);
        else switchVal_parton = 1;
      }

      if(switchVal_parton == 0) getPartonWLeps(thisPartonEventType, WPartonIndices, branchGenParticle, w1_parton, w2_parton, l1_parton, l2_parton, q1_parton, q2_parton);

      // partonET->Fill(thisPartonEventType,weight);

      if(foundWW) {
	
        double wwdeltaPhiparton = deltaPhi(w1_parton, w2_parton);
        double wwdeltaEtarparton = deltaEta(w1_parton, w2_parton);
        double wwdeltaRrparton = deltaR(w1_parton, w2_parton);

        // lepPT_partonV.at(0)->Fill(l1_parton.Pt());
        // lepPT_partonV.at(1)->Fill(l2_parton.Pt());
	
      }
    }
*/

  //------------------------------------------------------------------------------------------------------------------------------------------------------------
  // cuts
  //------------------------------------------------------------------------------------------------------------------------------------------------------------


  if ( found_bb_reco ) update_cft(cft_event_reco, cft_total_reco, "found bb - reco", weight);
  if ( found_jj_reco ) update_cft(cft_event_reco, cft_total_reco, "found jj - reco", weight);
  if ( found_ww_reco ) update_cft(cft_event_reco, cft_total_reco, "found ww - reco", weight);
  if ( found_zz_reco ) update_cft(cft_event_reco, cft_total_reco, "found zz - reco", weight);

  if ( found_bb_particle ) update_cft(cft_event_particle, cft_total_particle, "found bb - particle", weight);
  if ( found_jj_particle ) update_cft(cft_event_particle, cft_total_particle, "found jj - particle", weight);
  if ( found_ww_particle ) update_cft(cft_event_particle, cft_total_particle, "found ww - particle", weight);
  if ( found_zz_particle ) update_cft(cft_event_particle, cft_total_particle, "found zz - particle", weight);


  //------------------------------------------------------------------------------------------------------------------------------------------------------------
  // get tlorentz vectors
  //------------------------------------------------------------------------------------------------------------------------------------------------------------


  TLorentzVector b1_reco, b2_reco, b1_particle, b2_particle;
  TLorentzVector j1_reco, j2_reco, j1_particle, j2_particle;

  TLorentzVector l1_reco, l2_reco, l3_reco, l4_reco, l1_particle, l2_particle, l3_particle, l4_particle;
  TLorentzVector w1_reco, w2_reco, w1_particle, w2_particle;
  TLorentzVector z1_reco, z2_reco, z1_particle, z2_particle;

  // reco 

    if ( found_bb_reco ) bb_reco = make_bb( "reco", bb_reco_i ); b1_reco = bb_reco.first; b2_reco = bb_reco.second;
    if ( found_jj_reco ) jj_reco = make_jj( "reco", jj_reco_i, branchJet ); j1_reco = jj_reco.first; j2_reco = jj_reco.second;
    if ( found_ww_reco ) {
      ww_reco = make_ww_leptonic( "reco", ww_reco_i, branchElectron, branchMuon, branchMissingET ); 
      l1_reco = ww_reco.first[0]; l2_reco = ww_reco.first[1]; w1_reco = ww_reco.second[0]; w2_reco = ww_reco.second[1];
    }
    if ( found_zz_reco ) {
      zz_reco = make_zz_leptonic( "reco", zz_reco_i, branchElectron, branchMuon ); 
      l1_reco = zz_reco.first[0]; l2_reco = zz_reco.first[1]; l3_reco = zz_reco.first[2]; l4_reco = zz_reco.first[3]; z1_reco = zz_reco.second[0]; z2_reco = zz_reco.second[1];
    }

  // particle

    if ( found_bb_particle ) bb_particle = make_bb( "particle", bb_particle_i ); b1_particle = bb_particle.first; b2_particle = bb_particle.second;
    if ( found_jj_particle ) jj_particle = make_jj( "particle", jj_particle_i, branchGenJet ); j1_particle = jj_particle.first; j2_particle = jj_particle.second;
    if ( found_ww_particle ) {
      ww_particle = make_ww_leptonic( "particle", ww_particle_i, branchGenParticle, branchGenParticle, branchGenMissingET ); 
      l1_particle = ww_particle.first[0]; l2_particle = ww_particle.first[1]; w1_particle = ww_particle.second[0]; w2_particle = ww_particle.second[1];
    }
    if ( found_zz_particle ) {
      zz_particle = make_zz_leptonic( "particle", zz_particle_i, branchGenParticle, branchGenParticle ); 
      l1_particle = zz_particle.first[0]; l2_particle = zz_particle.first[1]; l3_particle = zz_particle.first[2]; l4_particle = zz_particle.first[3]; z1_particle = zz_particle.second[0]; z2_particle = zz_particle.second[1];
    }

  //------------------------------------------------------------------------------------------------------------------------------------------------------------
  // dump selection info
  //------------------------------------------------------------------------------------------------------------------------------------------------------------


  if ( debug_bool ) {

    debug_print( debug_bool, " " );
    debug_print( debug_bool, "------------------------------------------------------------------------" );
    debug_print( debug_bool, "reco selection" );
    debug_print( debug_bool, "------------------------------------------------------------------------" );
    debug_print( debug_bool, " " );

    if (found_bb_reco ) {
      debug_print( debug_bool, " " );
      debug_print( debug_bool, " reco hbb m: " + to_string( (b1_reco+b2_reco).M() ) + " hbb pT: " +  to_string( (b1_reco+b2_reco).Pt() ) );
      debug_print( debug_bool, " reco b1 pT: " + to_string( b1_reco.Pt() ) + " b1 eta: " + to_string( b1_reco.Eta() ) + " b1 phi: " + to_string( b1_reco.Phi() ) );
      debug_print( debug_bool, " reco b2 pT: " + to_string( b2_reco.Pt() ) + " b2 eta: " + to_string( b2_reco.Eta() ) + " b2 phi: " + to_string( b2_reco.Phi() ) );
    }
    if (!found_bb_reco) debug_print( debug_bool, " no reco hbb found" );
    if (found_jj_reco ) {
      debug_print( debug_bool, " " );
      debug_print( debug_bool, " reco jj m: " + to_string( (j1_reco+j2_reco).M() ) + " jj pT: " +  to_string( (j1_reco+j2_reco).Pt() ) );
      debug_print( debug_bool, " reco j1 pT: " + to_string( j1_reco.Pt() ) + " j1 eta: " + to_string( j1_reco.Eta() ) + " j1 phi: " + to_string( j1_reco.Phi() ) );
      debug_print( debug_bool, " reco j2 pT: " + to_string( j2_reco.Pt() ) + " j2 eta: " + to_string( j2_reco.Eta() ) + " j2 phi: " + to_string( j2_reco.Phi() ) );
    }
    if (!found_jj_reco) debug_print( debug_bool, " no reco jj found" );
    if (found_ww_reco ) {
      debug_print( debug_bool, " " );
      debug_print( debug_bool, " reco l1 pT: " + to_string( l1_reco.Pt() ) + " reco l1 eta: " + to_string( l1_reco.Eta() ) + " reco l1 phi: " + to_string( l1_reco.Phi() ) );
      debug_print( debug_bool, " reco l2 pT: " + to_string( l2_reco.Pt() ) + " reco l2 eta: " + to_string( l2_reco.Eta() ) + " reco l2 phi: " + to_string( l2_reco.Phi() ) );
      debug_print( debug_bool, " reco w1 pT: " + to_string( w1_reco.M() ) + " reco w1 mT: " + to_string( w1_reco.Mt() ) );
      debug_print( debug_bool, " reco w2 pT: " + to_string( w2_reco.M() ) + " reco w2 mT: " + to_string( w2_reco.Mt() ) );
      debug_print( debug_bool, " " );
    }
    if (!found_ww_reco) debug_print( debug_bool, " no reco ww found" );
    if (found_zz_reco ) {
      debug_print( debug_bool, " " );
      debug_print( debug_bool, " reco l1 pT: " + to_string( l1_reco.Pt() ) + " reco l1 eta: " + to_string( l1_reco.Eta() ) + " reco l1 phi: " + to_string( l1_reco.Phi() ) );
      debug_print( debug_bool, " reco l2 pT: " + to_string( l2_reco.Pt() ) + " reco l2 eta: " + to_string( l2_reco.Eta() ) + " reco l2 phi: " + to_string( l2_reco.Phi() ) );
      debug_print( debug_bool, " reco l3 pT: " + to_string( l3_reco.Pt() ) + " reco l3 eta: " + to_string( l3_reco.Eta() ) + " reco l3 phi: " + to_string( l3_reco.Phi() ) );
      debug_print( debug_bool, " reco l4 pT: " + to_string( l4_reco.Pt() ) + " reco l4 eta: " + to_string( l4_reco.Eta() ) + " reco l4 phi: " + to_string( l4_reco.Phi() ) );
      debug_print( debug_bool, " reco z1 pT: " + to_string( z1_reco.Pt() ) + " reco z1 m: " + to_string( z1_reco.M() ) );
      debug_print( debug_bool, " reco z2 pT: " + to_string( z2_reco.Pt() ) + " reco z2 m: " + to_string( z2_reco.M() ) );
      debug_print( debug_bool, " " );
    }
    if (!found_zz_reco) debug_print( debug_bool, " no reco zz found" );

    debug_print( debug_bool, " " );
    debug_print( debug_bool, "------------------------------------------------------------------------" );
    debug_print( debug_bool, "particle selection" );
    debug_print( debug_bool, "------------------------------------------------------------------------" );
    debug_print( debug_bool, " " );

    if (found_bb_particle ) {
      debug_print( debug_bool, " " );
      debug_print( debug_bool, " particle hbb m: " + to_string( (b1_particle+b2_particle).M() ) + " hbb pT: " +  to_string( (b1_particle+b2_particle).Pt() ) );
      debug_print( debug_bool, " particle b1 pT: " + to_string( b1_particle.Pt() ) + " b1 eta: " + to_string( b1_particle.Eta() ) + " b1 phi: " + to_string( b1_particle.Phi() ) );
      debug_print( debug_bool, " particle b2 pT: " + to_string( b2_particle.Pt() ) + " b2 eta: " + to_string( b2_particle.Eta() ) + " b2 phi: " + to_string( b2_particle.Phi() ) );
    }
    if (!found_bb_particle) debug_print( debug_bool, " no particle hbb found" );
    if (found_jj_particle ) {
      debug_print( debug_bool, " " );
      debug_print( debug_bool, " particle jj m: " + to_string( (j1_particle+j2_particle).M() ) + " jj pT: " +  to_string( (j1_particle+j2_particle).Pt() ) );
      debug_print( debug_bool, " particle j1 pT: " + to_string( j1_particle.Pt() ) + " j1 eta: " + to_string( j1_particle.Eta() ) + " j1 phi: " + to_string( j1_particle.Phi() ) );
      debug_print( debug_bool, " particle j2 pT: " + to_string( j2_particle.Pt() ) + " j2 eta: " + to_string( j2_particle.Eta() ) + " j2 phi: " + to_string( j2_particle.Phi() ) );
    }
    if (!found_jj_particle) debug_print( debug_bool, " no particle jj found" );
    if (found_ww_particle ) {
      debug_print( debug_bool, " " );
      debug_print( debug_bool, " particle l1 pT: " + to_string( l1_particle.Pt() ) + " particle l1 eta: " + to_string( l1_particle.Eta() ) + " particle l1 phi: " + to_string( l1_particle.Phi() ) );
      debug_print( debug_bool, " particle l2 pT: " + to_string( l2_particle.Pt() ) + " particle l2 eta: " + to_string( l2_particle.Eta() ) + " particle l2 phi: " + to_string( l2_particle.Phi() ) );
      debug_print( debug_bool, " particle w1 pT: " + to_string( w1_particle.Pt() ) + " particle w1 mT: " + to_string( w1_particle.Mt() ) );
      debug_print( debug_bool, " particle w2 pT: " + to_string( w2_particle.Pt() ) + " particle w2 mT: " + to_string( w2_particle.Mt() ) );
      debug_print( debug_bool, " " );
    }
    if (!found_ww_particle) debug_print( debug_bool, " no particle ww found" );
    if (found_zz_particle ) {
      debug_print( debug_bool, " " );
      debug_print( debug_bool, " particle l1 pT: " + to_string( l1_particle.Pt() ) + " particle l1 eta: " + to_string( l1_particle.Eta() ) + " particle l1 phi: " + to_string( l1_particle.Phi() ) );
      debug_print( debug_bool, " particle l2 pT: " + to_string( l2_particle.Pt() ) + " particle l2 eta: " + to_string( l2_particle.Eta() ) + " particle l2 phi: " + to_string( l2_particle.Phi() ) );
      debug_print( debug_bool, " particle l3 pT: " + to_string( l3_particle.Pt() ) + " particle l3 eta: " + to_string( l3_particle.Eta() ) + " particle l3 phi: " + to_string( l3_particle.Phi() ) );
      debug_print( debug_bool, " particle l4 pT: " + to_string( l4_particle.Pt() ) + " particle l4 eta: " + to_string( l4_particle.Eta() ) + " particle l4 phi: " + to_string( l4_particle.Phi() ) );
      debug_print( debug_bool, " particle z1 pT: " + to_string( z1_particle.Pt() ) + " particle z1 m: " + to_string( z1_particle.M() ) );
      debug_print( debug_bool, " particle z2 pT: " + to_string( z2_particle.Pt() ) + " particle z2 m: " + to_string( z2_particle.M() ) );
      debug_print( debug_bool, " " );
    }
    if (!found_zz_particle) debug_print( debug_bool, " no particle zz found" );

    debug_print( debug_bool, " " );
    debug_print( debug_bool, "------------------------------------------------------------------------" );
    debug_print( debug_bool, " " );



  }


  //------------------------------------------------------------------------------------------------------------------------------------------------------------
  // fill hists
  //------------------------------------------------------------------------------------------------------------------------------------------------------------


    for (const auto& selection : selections) {

      enable_cut_reco[ selection ] = true;
      enable_cut_particle[ selection ] = true;

      for (const auto& cut : cut_sel_process_reco[selection]) {
        if ( cft_event_reco[cut].first == 0 ) enable_cut_reco[ selection ] = false;
        if ( cft_event_particle[cut].first == 0 ) enable_cut_particle[ selection ] = false;
      }

      if ( enable_cut_reco[ selection ] ) sel_hist_map[ selection ].fill_bb( "reco", (b1_reco+b2_reco).Pt(), (b1_reco+b2_reco).M() , delta_eta(b1_reco,b2_reco) , delta_phi(b1_reco,b2_reco) , delta_r(b1_reco,b2_reco) );
      if ( enable_cut_reco[ selection ] ) sel_hist_map[ selection ].fill_jj( "reco", (j1_reco+j2_reco).Pt(), (j1_reco+j2_reco).M() , delta_eta(j1_reco,j2_reco) , delta_phi(j1_reco,j2_reco) , delta_r(j1_reco,j2_reco) );
      if ( enable_cut_reco[ selection ] ) sel_hist_map[ selection ].fill_ww( "reco", w1_reco.Pt(), w2_reco.Pt(), w1_reco.Mt(), w2_reco.Mt(), delta_eta(w1_reco,w2_reco) , delta_phi(w1_reco,w2_reco) , delta_r(w1_reco,w2_reco) );
      if ( enable_cut_reco[ selection ] ) sel_hist_map[ selection ].fill_zz( "reco", z1_reco.Pt(), z2_reco.Pt(), z1_reco.M(), z1_reco.Mt(), delta_eta(z1_reco,z2_reco) , delta_phi(z1_reco,z2_reco) , delta_r(z1_reco,z2_reco) );

      if ( enable_cut_particle[ selection ] ) sel_hist_map[ selection ].fill_bb( "particle", (b1_particle+b2_particle).Pt(), (b1_particle+b2_particle).M() , delta_eta(b1_particle,b2_particle) , delta_phi(b1_particle,b2_particle) , delta_r(b1_particle,b2_particle) );
      if ( enable_cut_particle[ selection ] ) sel_hist_map[ selection ].fill_jj( "particle", (j1_particle+j2_particle).Pt(), (j1_particle+j2_particle).M() , delta_eta(j1_particle,j2_particle) , delta_phi(j1_particle,j2_particle) , delta_r(j1_particle,j2_particle) );
      if ( enable_cut_particle[ selection ] ) sel_hist_map[ selection ].fill_ww( "particle", w1_particle.Pt(), w2_particle.Pt(), w1_particle.Mt(), w2_particle.Mt(), delta_eta(w1_particle,w2_particle) , delta_phi(w1_particle,w2_particle) , delta_r(w1_particle,w2_particle) );
      if ( enable_cut_particle[ selection ] ) sel_hist_map[ selection ].fill_zz( "particle", z1_particle.Pt(), z2_particle.Pt(), z1_particle.M(), z2_particle.M(), delta_eta(z1_particle,z2_particle) , delta_phi(z1_particle,z2_particle) , delta_r(z1_particle,z2_particle) );


    }


  //------------------------------------------------------------------------------------------------------------------------------------------------------------
  // end of event loop
  //------------------------------------------------------------------------------------------------------------------------------------------------------------

    if( entry % 1000 == 0 ){
      cout<<"Processed "<<entry<< " / " <<numberOfEntries <<" "<< entry/ numberOfEntries *100 <<" %"<<endl;
      PrintCutFlow( cft_total_reco, cut_list_reco, "Reco");
      PrintCutFlow( cft_total_particle, cut_list_particle, "Particle");
      PrintCutFlow( cutFlowMap_parton, cutList_parton, "Parton");
    }

    nPassed+=weight;
    nPassedRaw++;

    cutVal_reco++; cutValW_reco+=weight;
    cutVal_particle++; cutValW_particle+=weight;
    cutVal_parton++; cutValW_parton+=weight;

    // cout << "------------------------------" << endl;
    // cout << "Event: " << entry << endl;
    // printMap(cft_event_reco);
    // cout << "------------------------------" << endl;

    cft_allevents_reco.push_back(cft_event_reco);
    cft_allevents_particle.push_back(cft_event_particle);

  }


//------------------------------------------------------------------------------------------------------------------------------------------------------------
// fill and print cft
//------------------------------------------------------------------------------------------------------------------------------------------------------------

  cout << " " << endl;
  cout << "RECO CUT FLOW" << endl;
    PrintCutFlow(cft_total_reco, cut_list_reco, "Reco");
  cout << "PARTICLE CUT FLOW" << endl;
    PrintCutFlow(cft_total_particle, cut_list_particle, "Particle");
  cout << "PARTON CUT FLOW" << endl;
    PrintCutFlow(cutFlowMap_parton,cutList_parton, "Parton");
  
  for(std::vector<string>::iterator it=selType.begin(); it!=selType.end(); it++){
    FillCutFlow(cutFlowHists[(*it)],cutFlowEffs[(*it)],*cutFlowMapAll[(*it)],cutFlowMByType[(*it)], (*it));
  }
  
//------------------------------------------------------------------------------------------------------------------------------------------------------------
// WRITE HISTOGRAMS
//------------------------------------------------------------------------------------------------------------------------------------------------------------

  hists->cd();

  for (map <string, histograms>:: iterator hist = sel_hist_map.begin(); hist != sel_hist_map.end(); ++hist) {
  
    (*hist).second.write_all_hist( );

  }

  hists->Close();

}

//------------------------------------------------------------------------------------------------------------------------------------------------------------
// MAIN
//------------------------------------------------------------------------------------------------------------------------------------------------------------


int main(int argc, char* argv[]){

  const char *input_file = argv[1];
  const char *output_file = argv[2];
  const char *process_name = argv[3];
  const char *debug = argv[4];

  if (argc > 4 && strcmp(argv[4], "--debug") == 0) {

    debug_bool = true;
    cout << " running in debug mode " << endl;

  }
  
  zAnalyzer(input_file, output_file, process_name, debug);

  return 0;
}














//------------------------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------------------------
// STUFF I AM AFRAID TO DELETE PERMANENTLY
//------------------------------------------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------------------------------------------


  //------------------------------------------------------------------------------------------------------------------------------------------------------------
  // reco lep debug stuf
  //------------------------------------------------------------------------------------------------------------------------------------------------------------

/*
    int thisRecoEventType=-1;

    double e_pT_min = 15.0;
    double e_eta_max = 2.5;
    double mu_pT_min = 15.0;
    double mu_eta_max = 2.5;
      
    // get e+e- mu+mu-
    vector <int> e_min_indices_reco = get_leptons( "reco", branchElectron, e_pT_min, e_eta_max, "electron", -1 );
    vector <int> e_plus_indices_reco = get_leptons( "reco", branchElectron, e_pT_min, e_eta_max, "electron", 1 );
    vector <int> mu_min_indices_reco = get_leptons( "reco", branchMuon, mu_pT_min, mu_eta_max, "muon", -1 );
    vector <int> mu_plus_indices_reco = get_leptons( "reco", branchMuon, mu_pT_min, mu_eta_max, "muon", 1 );

  #ifdef MDEBUG
    cout << " ------------------------ " << endl;
    cout << "goodE_min_reco_indices: " << goodE_min_reco_indices.size() << endl;
    cout << "goodE_plus_reco_indices: " << goodE_plus_reco_indices.size() << endl;
    cout << "goodMu_min_reco_indices: " << goodMu_min_reco_indices.size() << endl;
    cout << "goodMu_plus_reco_indices: " << goodMu_plus_reco_indices.size() << endl;

    TLorentzVector emin1_reco, eplus1_reco, mmin1_reco, mplus1_reco;
    cout << " ------------------------ " << endl;
    
    if (goodE_min_reco_indices.size()>0) {
      emin1_reco = ((Electron *) branchElectron->At(goodE_min_reco_indices[0]))->P4();
      cout << "emin1_reco Pt: " << emin1_reco.Pt() << " emin1_reco Eta: " << emin1_reco.Eta() << endl;
    }
    if (goodE_plus_reco_indices.size()>0) {
      eplus1_reco = ((Electron *) branchElectron->At(goodE_plus_reco_indices[0]))->P4();
      cout << "eplus1_reco Pt: " << eplus1_reco.Pt() << " eplus1_reco Eta: " << eplus1_reco.Eta() << endl;
    }
    if (goodMu_min_reco_indices.size()>0) {
      mmin1_reco = ((Muon *) branchMuon->At(goodMu_min_reco_indices[0]))->P4();
      cout << "mmin1_reco Pt: " << mmin1_reco.Pt() << " mmin1_reco Eta: " << mmin1_reco.Eta() << endl;
    }
    if (goodMu_plus_reco_indices.size()>0) {
      mplus1_reco = ((Muon *) branchMuon->At(goodMu_plus_reco_indices[0]))->P4();
      cout << "mplus1_reco Pt: " << mplus1_reco.Pt() << " mplus1_reco Eta: " << mplus1_reco.Eta() << endl;
    } 
  #endif

    // e = e+ & e- , mu = mu+ & mu- 
    vector<int> goodE_reco_indices;
    vector<int> goodMu_reco_indices;
    concatenate_indices(goodE_min_reco_indices, goodE_reco_indices);
    concatenate_indices(goodE_plus_reco_indices, goodE_reco_indices);
    concatenate_indices(goodMu_min_reco_indices, goodMu_reco_indices);
    concatenate_indices(goodMu_plus_reco_indices, goodMu_reco_indices);
*/

  //------------------------------------------------------------------------------------------------------------------------------------------------------------
  // PARTICLE - HIGGS
  //------------------------------------------------------------------------------------------------------------------------------------------------------------
/*
    int switchVal_particle = 0;
    bool foundHiggs_particle = false;

    if(enableCutParticle["initial - particle"]){
      increaseCount(cutFlowMap_particle,"initial - particle",weight);
    }
 
    vector <int> goodJetIndexParticle=GoodJetIndices(branchGenJet);

    if(enableCutParticle["jet pT > 20 - particle"]) {
      if(switchVal_particle == 0 && goodJetIndexParticle.size() > 0) increaseCount(cutFlowMap_particle,"jet pT > 20 - particle", weight);
      else switchVal_particle = 1;
    }

    std::vector<std::pair< std::map<TString, float>, std::map<TString, std::vector<float>>>>  pairedJetParticle=paired::PAIReDjointEvent(branchGenParticle,branchGenParticle,branchGenJet,0.4,false,false,true,1.0,false);
    //cout<<"PAIRED lables bb "<<pairedJet.first["label_bb"]<<" cc "<<pairedJet.first["label_cc"]<<" ll "<<pairedJet.first["label_ll"]<<" indices 1: "<<pairedJet.first["jet1_index"]<<" 2: "<<pairedJet.first["jet1_index"]<<endl;
    std::vector<std::pair< std::map<TString, float>, std::map<TString, std::vector<float>>>>  pairedJetBParticle;

    if(enableCutParticle["1 PAIReD jet - particle"]) {
      if(switchVal_particle == 0 && pairedJetParticle.size() > 0) increaseCount(cutFlowMap_particle,"1 PAIReD jet - particle",weight);
      else  switchVal_particle = 1;
    }

    for(int i=0; i<(int)pairedJetParticle.size(); i++){
      std::pair< std::map<TString, float>, std::map<TString, std::vector<float>>> thisPairedParticle=pairedJetParticle.at(i);
      if( thisPairedParticle.first["isbtagged"] > 0) {
        pairedJetBParticle.push_back(thisPairedParticle);
      }
    }

    vector <int> btagIndexParticle;
    int pairedJetSize_particle = pairedJetParticle.size();
    int pairedBJetSize_particle = pairedJetBParticle.size();

    if(enableCutParticle["1 bb PAIReD jet - particle"]) {
      if(switchVal_particle == 0 && pairedJetBParticle.size()>0){
        increaseCount(cutFlowMap_particle,"1 bb PAIReD jet - particle",weight);
        foundHiggs_particle = true;
      } else switchVal_particle = 1;
    }

    std::map<TString, float> paired_jet_particle;

    if (switchVal_particle == 0 && foundHiggs_particle){

      paired_jet_particle = pairedJetBParticle.at(0).first;

      btagIndexParticle.push_back(paired_jet_particle["jet1_index"]);
      btagIndexParticle.push_back(paired_jet_particle["jet2_index"]);

      b1_particle.SetPtEtaPhiM(paired_jet_particle["jet1_pt"],paired_jet_particle["jet1_eta"],paired_jet_particle["jet1_phi"],paired_jet_particle["jet1_mass"]);
      b2_particle.SetPtEtaPhiM(paired_jet_particle["jet2_pt"],paired_jet_particle["jet2_eta"],paired_jet_particle["jet2_phi"],paired_jet_particle["jet2_mass"]);

      h_particle = b1_particle + b2_particle; // dijet

      double bbdeltaPhiparticle = deltaPhi(b1_particle, b2_particle);
      double bbdeltaEtarparticle = deltaEta(b1_particle, b2_particle);
      double bbdeltaRrparticle = deltaR(b1_particle, b2_particle);
 
    }
 
    //------------------------------------------------------------------------------------------------------------------------------------------------------------
    // PARTICLE - VBF JETS
    //------------------------------------------------------------------------------------------------------------------------------------------------------------


    vector <int> nonHiggsJetParticle;
    vector<pair<int,int>> vbfJetIndexParticle;
    vector<vector <int>> vbfJetIndexParticleComb;
    vector<pair<int,int>> vbfJetIndexParticle_dEta;
    int vbfJetIndexParticleCandidate = -1;

    // check that they do not belong to higgs + sort by pT
    for(int i=0; i<(int)goodJetIndexParticle.size(); i++) {
      if( goodJetIndexParticle[i] == paired_jet_particle["jet1_index"] || goodJetIndexParticle[i] == paired_jet_particle["jet2_index"] ) continue;
    nonHiggsJetParticle.push_back(goodJetIndexParticle[i]);
    SortByPtIndices(nonHiggsJetParticle,branchGenJet);
    }

    // check that there are at least two + make combinations
    if(enableCutParticle["2 VBF jet - particle"]) {
      if(switchVal_particle==0 && nonHiggsJetParticle.size() > 1 ) {
	      increaseCount(cutFlowMap_particle,"2 VBF jet - particle",weight);
	      vbfJetIndexParticleComb=combinationsNoRepetitionAndOrderDoesNotMatter(2,nonHiggsJetParticle);
      }
      else switchVal_particle=1;
    }

    // convert vector of vector of ints to vector of pairs of ints
    vbfJetIndexParticle=GetvbfJetIndex(vbfJetIndexParticleComb);

    // loop and take those w dEta > 2.5
    for (int i=0; i<(int)vbfJetIndexParticle.size(); i++) {
      if( fabs((((Jet*)branchGenJet->At(vbfJetIndexParticle[i].first))->Eta - ((Jet*)branchGenJet->At(vbfJetIndexParticle[i].second))->Eta)) <= 2.5 ) {
	      continue; 
      } else { 
	      vbfJetIndexParticle_dEta.push_back(vbfJetIndexParticle.at(i));
	      vbfJetIndexParticleCandidate = i;
      }
    }

    if(enableCutParticle["2.5 deltaEta VBF jet - particle"]) {
      if(switchVal_particle==0 && vbfJetIndexParticle_dEta.size()>0) increaseCount(cutFlowMap_particle,"2.5 deltaEta VBF jet - particle",weight);
      else switchVal_particle=1;
    }

    // sort them again by eta for leading/subleading
    SortByEtaIndices(vbfJetIndexParticle_dEta,branchGenJet);

    Jet *jet1_particle =nullptr;
    Jet *jet2_particle =nullptr;
  
    if(switchVal_particle==0 && vbfJetIndexParticle.size()>0) {

      jet1_particle = (Jet*) branchGenJet->At(vbfJetIndexParticle[0].first);
      jet2_particle = (Jet*) branchGenJet->At(vbfJetIndexParticle[0].second);
      j1_particle=jet1_particle->P4();
      j2_particle=jet2_particle->P4();

      double jjdeltaPhiparticle = deltaPhi(j1_particle, j2_particle);
      double jjdeltaEtarparticle = deltaEta(j1_particle, j2_particle);
      double jjdeltaRrparticle = deltaR(j1_particle, j2_particle);
  
    }


  //------------------------------------------------------------------------------------------------------------------------------------------------------------
  // PARTICLE - LEPTONS 
  //------------------------------------------------------------------------------------------------------------------------------------------------------------

  int thisParticleEventType=-1;

  double lepPTMinFid = 15.0;
  double lepEtaMaxFid = 2.5;

  // get e+e- mu+mu-
  vector <int> goodE_min_particle_indices = get_good_particle_lepton_indices(branchGenParticle, lepPTMinFid, lepEtaMaxFid, analysis, 11);
  vector <int> goodE_plus_particle_indices = get_good_particle_lepton_indices(branchGenParticle,lepPTMinFid, lepEtaMaxFid, analysis, -11);
  vector <int> goodMu_min_particle_indices = get_good_particle_lepton_indices(branchGenParticle, lepPTMinFid, lepEtaMaxFid, analysis, 13);
  vector <int> goodMu_plus_particle_indices = get_good_particle_lepton_indices(branchGenParticle, lepPTMinFid, lepEtaMaxFid, analysis, -13);

#ifdef MDEBUG
  cout << " ------------------------ " << endl;
  cout << "goodE_min_particle_indices: " << goodE_min_particle_indices.size() << endl;
  cout << "goodE_plus_particle_indices: " << goodE_plus_particle_indices.size() << endl;
  cout << "goodMu_min_particle_indices: " << goodMu_min_particle_indices.size() << endl;
  cout << "goodMu_plus_particle_indices: " << goodMu_plus_particle_indices.size() << endl;

  TLorentzVector emin1_particle, eplus1_particle, mmin1_particle, mplus1_particle;
  cout << " ------------------------ " << endl;

  if (goodE_min_particle_indices.size()>0) {
    emin1_particle = ((GenParticle *) branchGenParticle->At(goodE_min_particle_indices[0]))->P4();
    cout << "emin1_particle Pt: " << emin1_particle.Pt() << " emin1_particle Eta: " << emin1_particle.Eta() << endl;
  }
  if (goodE_plus_particle_indices.size()>0) {
    eplus1_particle = ((GenParticle *) branchGenParticle->At(goodE_plus_particle_indices[0]))->P4();
    cout << "eplus1_particle Pt: " << eplus1_particle.Pt() << " eplus1_particle Eta: " << eplus1_particle.Eta() << endl;
  }
  if (goodMu_min_particle_indices.size()>0) {
    mmin1_particle = ((GenParticle *) branchGenParticle->At(goodMu_min_particle_indices[0]))->P4();
    cout << "mmin1_particle Pt: " << mmin1_particle.Pt() << " mmin1_particle Eta: " << mmin1_particle.Eta() << endl;
  }
  if (goodMu_plus_particle_indices.size()>0) {
    mplus1_particle = ((GenParticle *) branchGenParticle->At(goodMu_plus_particle_indices[0]))->P4();
    cout << "mplus1_particle Pt: " << mplus1_particle.Pt() << " mplus1_particle Eta: " << mplus1_particle.Eta() << endl;
  } 
#endif 

  // e = e+ & e- , mu = mu+ & mu- 
  vector<int> goodE_particle_indices;
  vector<int> goodMu_particle_indices;
  concatenate_indices( goodE_min_particle_indices,goodE_particle_indices);
  concatenate_indices( goodE_plus_particle_indices,goodE_particle_indices);
  concatenate_indices( goodMu_min_particle_indices,goodMu_particle_indices);
  concatenate_indices( goodMu_plus_particle_indices,goodMu_particle_indices);

  vector<int> goodLep_particle_indices;
  concatenate_indices( goodE_particle_indices,goodLep_particle_indices);
  concatenate_indices( goodMu_particle_indices,goodLep_particle_indices);

#ifdef MDEBUG
  // Debug 
  cout<<" After selection  particle:"<<endl;
  for(std::vector<int>::iterator it=goodLep_particle_indices.begin(); it!=goodLep_particle_indices.end(); it++){
    GenParticle *particle=(GenParticle*) branchGenParticle->At(*it); 
    dumpParticle(particle,*it);
  }
  // lep details
  cout << "E/Mu: " << goodE_particle_indices.size() << "/" << goodMu_particle_indices.size() << endl;
  //end debug 
#endif 
  
  goodE_size_particle->Fill(goodE_particle_indices.size(),weight);
  goodMu_size_particle->Fill(goodMu_particle_indices.size(),weight);

  if (goodE_particle_indices.size()>0) e1_particle = ((GenParticle *) branchGenParticle->At(goodE_particle_indices[0]))->P4();
  if (goodE_particle_indices.size()>1) e2_particle = ((GenParticle *) branchGenParticle->At(goodE_particle_indices[1]))->P4();
  if (goodMu_particle_indices.size()>0) m1_particle = ((GenParticle *) branchGenParticle->At(goodMu_particle_indices[0]))->P4();
  if (goodMu_particle_indices.size()>1) m2_particle = ((GenParticle *) branchGenParticle->At(goodMu_particle_indices[1]))->P4();

  lead_e_pt_particle->Fill(e1_particle.Pt(),weight);
  lead_e_eta_particle->Fill(e1_particle.Eta(),weight);
  lead_e_phi_particle->Fill(e1_particle.Phi(),weight);
  sublead_e_pt_particle->Fill(e2_particle.Pt(),weight);
  sublead_e_eta_particle->Fill(e2_particle.Eta(),weight);
  sublead_e_phi_particle->Fill(e2_particle.Phi(),weight);

  lead_mu_pt_particle->Fill(m1_particle.Pt(),weight);
  lead_mu_eta_particle->Fill(m1_particle.Eta(),weight);
  lead_mu_phi_particle->Fill(m1_particle.Phi(),weight);
  sublead_mu_pt_particle->Fill(m2_particle.Pt(),weight);
  sublead_mu_eta_particle->Fill(m2_particle.Eta(),weight);
  sublead_mu_phi_particle->Fill(m2_particle.Phi(),weight);
  
// V details

  vector<pair<int,pair<int,int>>> ZParticlePairIndices;
  vector<pair<pair<int,int>,int>> WParticlePairIndices;

  if(enableCutParticle["lep pT & eta cut - particle"]){
    if (switchVal_particle==0) {
      if (goodE_particle_indices.size() > 0 || goodMu_particle_indices.size() > 0) increaseCount(cutFlowMap_particle,"lep pT & eta cut - particle",weight);
    } else  switchVal_particle=1;
  }

  if(analysis == "HZZJJ"){

    // form pairs for each flavour
    vector< pair<int,int>> elecZParticlePairIndices=GetelecParticlePairIndices(branchGenParticle,goodE_particle_indices); 
    vector< pair<int,int>> muZParticlePairIndices=GetmuParticlePairIndices(branchGenParticle,goodMu_particle_indices);
    
    ZParticlePairIndices=GetParticlePairIndices(elecZParticlePairIndices,muZParticlePairIndices,branchGenParticle); // 0 for electron 1 for muon

    if(enableCutParticle["OSSF - particle"]){
      if (switchVal_particle ==0 && ZParticlePairIndices.size()>=2 )  increaseCount(cutFlowMap_particle,"OSSF - particle",weight);
      else switchVal_particle=1;
    }
         
    if(switchVal_particle==0 && ZParticlePairIndices.size()>=2){
      if( ZParticlePairIndices[0].first == 1 && ZParticlePairIndices[1].first == 1) thisParticleEventType=0;
      else if( ZParticlePairIndices[0].first == 0 && ZParticlePairIndices[1].first == 0) thisParticleEventType=1;
      else if( ZParticlePairIndices[0].first == 1 && ZParticlePairIndices[1].first == 0) thisParticleEventType=2;
      else if( ZParticlePairIndices[0].first == 0 && ZParticlePairIndices[1].first == 1) thisParticleEventType=3;
    }

    particleET->Fill(thisParticleEventType, weight);
    
    getParticleZLeps(thisParticleEventType, ZParticlePairIndices, branchGenParticle, l1_particle, l2_particle, l3_particle, l4_particle, q1_particle, q2_particle, q3_particle, q4_particle);

    if( switchVal_particle == 0 && thisParticleEventType != -1 && ZParticlePairIndices.size() >= 2 ){
      z1_particle=l1_particle + l2_particle;
      z2_particle=l3_particle + l4_particle;
      double zzdeltaPhiparticle = deltaPhi(z1_particle, z2_particle);
      double zzdeltaEtarparticle = deltaEta(z1_particle, z2_particle);
      double zzdeltaRrparticle = deltaR(z1_particle, z2_particle);
    }

    // WW 

   } else if(analysis == "HWWJJ") {

    WParticlePairIndices = GetWParticlePairIndices(goodE_particle_indices, goodMu_particle_indices, branchGenParticle, branchMissingET);

#ifdef MDEBUG
cout<<" switch val "<<switchVal_particle<<endl;
#endif     
    
    // FOR OFOS SWITCH mu mu / e e EVENT TYPE TO -1
    if((goodE_particle_indices.size() + goodMu_particle_indices.size()) >= 2 ){

#ifdef MDEBUG
cout<<" Mumin "	<<goodMu_min_particle_indices.size()<<" Muplus "<<goodMu_plus_particle_indices.size()<<" Emin "<<goodE_min_particle_indices.size()<<" Eplus "<<goodE_plus_particle_indices.size()<<endl;
#endif

    if (goodMu_min_particle_indices.size() > 0 && goodMu_plus_particle_indices.size() > 0) thisParticleEventType = 0;
    if (goodE_min_particle_indices.size() > 0 && goodE_plus_particle_indices.size() > 0) thisParticleEventType = 1;
    if (goodMu_min_particle_indices.size() > 0 && goodE_plus_particle_indices.size() > 0) thisParticleEventType = 2;
    if (goodE_min_particle_indices.size() > 0 && goodMu_plus_particle_indices.size() > 0) thisParticleEventType = 3;
    }

#ifdef MDEBUG
cout << thisParticleEventType << endl;
#endif

      particleET->Fill(thisParticleEventType);

      if(enableCutParticle["OSOF - particle"]){
        if (switchVal_particle==0 && thisParticleEventType != -1) increaseCount(cutFlowMap_particle,"OSOF - particle",weight);
        else switchVal_particle=1;
      }

      getWParticle(thisParticleEventType, WParticlePairIndices, branchGenParticle,branchMissingET, l1_particle, l2_particle, q1_particle, q2_particle, met1, met2);

      hllm_0_15_particle->Fill((l1_particle+l2_particle).M(),weight);

      if(enableCutParticle["mll > 10 - particle"]) {
        if(switchVal_particle == 0 &&  (l1_particle+l2_particle).M() >= 10) increaseCount(cutFlowMap_particle,"mll > 10 - particle",weight);
        else switchVal_particle = 1;
      }

      if( switchVal_particle == 0 ) {

        w1_particle=l1_particle + met1;
        w2_particle=l2_particle + met2;
        double wwdeltaPhiparticle = deltaPhi(w1_particle, w2_particle);
        double wwdeltaEtarparticle = deltaEta(w1_particle, w2_particle);
        double wwdeltaRrparticle = deltaR(w1_particle, w2_particle);

      }

    }
*/




  //------------------------------------------------------------------------------------------------------------------------------------------------------------
  // HISTOGRAMS - RECO
  //------------------------------------------------------------------------------------------------------------------------------------------------------------

  // 1D
/*
  if () recoET->Fill(ww.first);
  if () recoET->Fill(zz.first);

    // higgs
    if(switchVal_reco==0){
      if(foundHiggs_reco){
        hHpTreco -> Fill(h_reco.Pt(),weight);
        hHmreco -> Fill(h_reco.M(),weight); 
	      hbbdeltaPhireco -> Fill(deltaPhi(b1_reco, b2_reco),weight);
	      hbbdeltaEtareco -> Fill(deltaEta(b1_reco, b2_reco),weight);
	      hbbdeltaRreco -> Fill(deltaR(b1_reco, b2_reco),weight);
      }
    }

    if(switchVal_particle==0){
      if(foundHiggs_particle){
        hHpTparticle -> Fill(h_particle.Pt(),weight);
        hHmparticle -> Fill(h_particle.M(),weight); 
	      hbbdeltaPhiparticle -> Fill(deltaPhi(b1_particle, b2_particle),weight);
	      hbbdeltaEtaparticle -> Fill(deltaEta(b1_particle, b2_particle),weight);
	      hbbdeltaRparticle -> Fill(deltaR(b1_particle, b2_particle),weight);
      }
    }

    if(switchVal_parton==0){
      if(HiggsRecord){
        hHpTparton -> Fill(h_parton.Pt(),weight);
        hHmparton -> Fill(h_parton.M(),weight); 
	      hbbdeltaPhiparton -> Fill(deltaPhi(b1_parton, b2_parton),weight);
	      hbbdeltaEtaparton -> Fill(deltaEta(b1_parton, b2_parton),weight);
	      hbbdeltaRparton -> Fill(deltaR(b1_parton, b2_parton),weight);
      }
    }

    // vbfj
    if(switchVal_reco==0){
      if(vbfJetIndex.size()>0){
        hjjpTreco->Fill(j1_reco.Pt()+j2_reco.Pt(),weight);
        hjjdeltaPhireco->Fill(deltaPhi(j1_reco, j2_reco),weight); 
        hjjdeltaEtareco->Fill(deltaEta(j1_reco, j2_reco), weight);
        hjjdeltaRreco -> Fill(deltaR(j1_reco, j2_reco),weight);
      }
    }

    if(switchVal_particle==0){
      if(vbfJetIndexParticle.size()>0){
        hjjpTparticle->Fill(j1_particle.Pt()+j2_particle.Pt(),weight);
        hjjdeltaPhiparticle->Fill(deltaPhi(j1_particle, j2_particle),weight); 
        hjjdeltaEtaparticle->Fill(deltaEta(j1_particle, j2_particle), weight);
        hjjdeltaRparticle -> Fill(deltaR(j1_particle, j2_particle),weight);
      }
    }

    // z
    if(switchVal_reco==0){
      if(thisRecoEventType!=-1 && ZRecoPairIndices.size()>=2){
        hZ1pTreco->Fill(z1_reco.Pt(),weight);
        hZ2pTreco->Fill(z2_reco.Pt(),weight);
        hZ1mreco->Fill(z1_reco.M(),weight);
        hZ2mreco->Fill(z2_reco.M(),weight);
        hZZdeltaPhireco->Fill(zzdeltaPhireco,weight); 
        hZZdeltaEtareco->Fill(zzdeltaEtareco, weight);
        hZZdeltaRreco -> Fill(zzdeltaRreco,weight);
      }
    }

    if(switchVal_particle==0){
      if(thisParticleEventType!=-1 && ZParticlePairIndices.size()>=2){
        hZ1pTparticle->Fill(z1_particle.Pt(),weight);
        hZ2pTparticle->Fill(z2_particle.Pt(),weight);
        hZ1mparticle->Fill(z1_particle.M(),weight);
        hZ2mparticle->Fill(z2_particle.M(),weight);
        hZZdeltaPhiparticle->Fill(zzdeltaPhiparticle,weight); 
        hZZdeltaEtaparticle->Fill(zzdeltaEtaparticle, weight);
        hZZdeltaRparticle -> Fill(zzdeltaRparticle,weight);
      }
    }

    if(switchVal_parton==0 ){
      if(foundZZ){
        hZ1pTparton->Fill(z1_parton.Pt(),weight);
        hZ2pTparton->Fill(z2_parton.Pt(),weight);
        hZ1mparton->Fill(z1_parton.M(),weight);
        hZ2mparton->Fill(z2_parton.M(),weight);
        hZZdeltaPhiparton->Fill(zzdeltaPhiparton,weight); 
        hZZdeltaEtaparton->Fill(zzdeltaEtaparton, weight);
        hZZdeltaRparton -> Fill(zzdeltaRparton,weight);
      }
    }

    // w 
    if(switchVal_reco==0){
      if(thisRecoEventType!=-1 && wleps.size()>=2){
        hW1pTreco->Fill(w1_reco.Pt(),weight);
        hW1mreco->Fill(massTransverse(l1_reco, met),weight);
        hW2pTreco->Fill(w2_reco.Pt(),weight);
        hW2mreco->Fill(massTransverse(l2_reco, met),weight);

        //hWWpTreco->Fill((w1_reco + w2_reco).Pt(),weight);
        //hWWmreco->Fill((w1_reco + w2_reco).M(),weight);
        //hllpTreco->Fill((l1_reco+l2_reco).Pt(),weight);
        //hllmreco->Fill((l1_reco+l2_reco).M(),weight);

        hWWdeltaPhireco->Fill(wwdeltaPhireco,weight); 
        hWWdeltaEtareco->Fill(wwdeltaEtareco, weight);
        hWWdeltaRreco -> Fill(wwdeltaRreco,weight);
      }
    }

    if(switchVal_particle==0){
      if(thisParticleEventType!=-1 && WParticlePairIndices.size()>=2){
        hW1pTparticle->Fill(w1_particle.Pt(),weight);
        hW1mparticle->Fill(massTransverse(l1_particle, met),weight);
        hW2pTparticle->Fill(w2_particle.Pt(),weight);
        hW2mparticle->Fill(massTransverse(l2_particle, met),weight);

        //hllpTparticle->Fill((l1_particle+l2_particle).Pt(),weight);
        //hllmparticle->Fill((l1_particle+l2_particle).M(),weight);
        //hWWpTparticle->Fill((w1_particle + w2_sparticle).Pt(),weight);
        //hWWmparticle->Fill((w1_particle + w2_particle).M(),weight);

        hWWdeltaPhiparticle->Fill(wwdeltaPhiparticle,weight); 
        hWWdeltaEtaparticle->Fill(wwdeltaEtaparticle, weight);
        hWWdeltaRparticle -> Fill(wwdeltaRparticle,weight);
      }
    }

    if(switchVal_parton==0 ){
      if(foundWW){
        hW1pTparton->Fill(w1_parton.Pt(),weight);
        hW2pTparton->Fill(w2_parton.Pt(),weight);
        hW1mparton->Fill(w1_parton.M(),weight);
        hW2mparton->Fill(w2_parton.M(),weight);
        hWWdeltaPhiparton->Fill(wwdeltaPhiparton,weight); 
        hWWdeltaEtaparton->Fill(wwdeltaEtaparton, weight);
        hWWdeltaRparton -> Fill(wwdeltaRparton,weight);
      }
    }

  //------------------------------------------------------------------------------------------------------------------------------------------------------------
  // FILL HISTOGRAMS - 2D
  //------------------------------------------------------------------------------------------------------------------------------------------------------------


  // 2D - parton(1) particle(2) reco(3)

    if(switchVal_parton==0 && switchVal_particle==0 ){
      hHpT12Comp -> Fill(h_parton.Pt(), h_particle.Pt(), weight);
      hHm12Comp -> Fill(h_parton.M(), h_particle.M(), weight);
      hbbdeltaPhi12Comp -> Fill(bbdeltaPhiparton, bbdeltaPhiparticle, weight);
      hbbdeltaEta12Comp -> Fill(bbdeltaEtaparton, bbdeltaEtaparticle, weight);
      hjjpT12Comp -> Fill(j1_parton.Pt()+j2_parton.Pt(),j1_particle.Pt()+j2_particle.Pt(), weight);
      hjjdeltaPhi12Comp -> Fill(jjdeltaPhiparton, jjdeltaPhiparticle, weight);

      hz1pT12Comp->Fill(z1_parton.Pt(), z1_particle.Pt(), weight);
      hz1m12Comp->Fill(z1_parton.M(), z1_particle.M(), weight);
      hz2pT12Comp->Fill(z2_parton.Pt(), z2_particle.Pt(), weight);
      hz2m12Comp->Fill(z2_parton.M(), z2_particle.M(), weight);

      hw1pT12Comp->Fill(w1_parton.Pt(), w1_particle.Pt(), weight);
      hw1m12Comp->Fill(w1_parton.M(), w1_particle.M(), weight);
      hw2pT12Comp->Fill(w2_parton.Pt(), w2_particle.Pt(), weight);
      hw2m12Comp->Fill(w2_parton.M(), w2_particle.M(), weight);

      hl1pT12Comp->Fill(l1_parton.Pt(), l1_particle.Pt(), weight);
      hl2pT12Comp->Fill(l2_parton.Pt(), l2_particle.Pt(), weight);
    }
    if(switchVal_particle==0  && switchVal_reco==0 ){
      hHpT23Comp -> Fill(h_particle.Pt(), h_reco.Pt(), weight);
      hHm23Comp -> Fill(h_particle.M(), h_reco.M(), weight);
      hbbdeltaPhi23Comp -> Fill(bbdeltaPhiparticle, bbdeltaPhireco, weight);
      hbbdeltaEta23Comp -> Fill(bbdeltaEtaparticle, bbdeltaEtareco, weight);
      hjjpT23Comp -> Fill(j1_particle.Pt()+j2_particle.Pt(),j1_reco.Pt()+j2_reco.Pt(), weight);
      hjjdeltaPhi23Comp -> Fill(jjdeltaPhiparticle, jjdeltaPhireco, weight);

      hz1pT23Comp->Fill(z1_particle.Pt(), z1_reco.Pt(), weight);
      hz1m23Comp->Fill(z1_particle.M(), z1_reco.M(), weight);
      hz2pT23Comp->Fill(z2_particle.Pt(), z2_reco.Pt(), weight);
      hz2m23Comp->Fill(z2_particle.M(), z2_reco.M(), weight);

      hw1pT23Comp->Fill(w1_particle.Pt(), w1_reco.Pt(), weight);
      hw1m23Comp->Fill(w1_particle.M(), w1_reco.M(), weight);
      hw2pT23Comp->Fill(w2_particle.Pt(), w2_reco.Pt(), weight);
      hw2m23Comp->Fill(w2_particle.M(), w2_reco.M(), weight);

      hl1pT23Comp->Fill(l1_particle.Pt(), l1_reco.Pt(), weight);
      hl2pT23Comp->Fill(l2_particle.Pt(), l2_reco.Pt(), weight);

      //leadbscore23->Fill(leadbscore_particle,leadbscore_reco,weight);
      //subleadbscore23->Fill(subleadbscore_particle,subleadbscore_reco,weight);
    }
    if(switchVal_parton==0  && switchVal_reco==0 ){
      hHpT13Comp -> Fill(h_parton.Pt(), h_reco.Pt(), weight);
      hHm13Comp -> Fill(h_parton.M(), h_reco.M(), weight);
      hbbdeltaPhi13Comp -> Fill(bbdeltaPhiparton, bbdeltaPhireco, weight);
      hbbdeltaEta13Comp -> Fill(bbdeltaEtaparton, bbdeltaEtareco, weight);
      hjjpT13Comp -> Fill(j1_parton.Pt()+j2_parton.Pt(),j1_reco.Pt()+j2_reco.Pt(), weight);
      hjjdeltaPhi13Comp -> Fill(jjdeltaPhiparton, jjdeltaPhireco, weight);

      hz1pT13Comp->Fill(z1_parton.Pt(), z1_reco.Pt(), weight);
      hz1m13Comp->Fill(z1_parton.M(), z1_reco.M(), weight);
      hz2pT13Comp->Fill(z2_parton.Pt(), z2_reco.Pt(), weight);
      hz2m13Comp->Fill(z2_parton.M(), z2_reco.M(), weight);

      hw1pT13Comp->Fill(w1_parton.Pt(), w1_reco.Pt(), weight);
      hw1m13Comp->Fill(w1_parton.M(), w1_reco.M(), weight);
      hw2pT13Comp->Fill(w2_parton.Pt(), w2_reco.Pt(), weight);
      hw2m13Comp->Fill(w2_parton.M(), w2_reco.M(), weight);

      hl1pT13Comp->Fill(l1_parton.Pt(), l1_reco.Pt(), weight);
      hl2pT13Comp->Fill(l2_parton.Pt(), l2_reco.Pt(), weight);
    }
   
    if(switchVal_reco==0){
      hbbjjdeltaPhicompreco->Fill(bbdeltaPhireco,jjdeltaPhireco,weight);
      hbbdeltaEtajjdeltaPhicompreco->Fill(bbdeltaEtareco,jjdeltaPhireco,weight);
      hHpTl1l2deltaPhicompreco->Fill(h_reco.Pt(), l1l2deltaPhireco, weight);
      hHpTl3l4deltaPhicompreco->Fill(h_reco.Pt(), l3l4deltaPhireco, weight);
      hHz1pTcompreco->Fill(h_reco.Pt(), z1_reco.Pt(), weight);
      hHz2pTcompreco->Fill(h_reco.Pt(), z2_reco.Pt(), weight);
      hHzzpTcompreco->Fill(h_reco.Pt(), z1_reco.Pt() + z2_reco.Pt(), weight);
      hHpTzzdeltaPhicompreco->Fill(h_reco.Pt(), zzdeltaPhireco, weight);
      hHpTzzdeltaEtacompreco->Fill(h_reco.Pt(), zzdeltaEtareco, weight);
      hbbzzdeltaPhicompreco->Fill(bbdeltaPhireco, zzdeltaPhireco, weight);
      hbbzzdeltaEtacompreco->Fill(bbdeltaEtareco, zzdeltaEtareco, weight);   
    }

    if(switchVal_particle==0){
      hbbjjdeltaPhicompparticle->Fill(bbdeltaPhiparticle,jjdeltaPhiparticle,weight);
      hbbdeltaEtajjdeltaPhicompparticle->Fill(bbdeltaEtaparticle,jjdeltaPhiparticle,weight);
      hHpTl1l2deltaPhicompparticle->Fill(h_particle.Pt(), l1l2deltaPhiparticle, weight);
      hHpTl3l4deltaPhicompparticle->Fill(h_particle.Pt(), l3l4deltaPhiparticle, weight);
      hHz1pTcompparticle->Fill(h_particle.Pt(), z1_particle.Pt(), weight);
      hHz2pTcompparticle->Fill(h_particle.Pt(), z2_particle.Pt(), weight);
      hHzzpTcompparticle->Fill(h_particle.Pt(), z1_particle.Pt() + z2_particle.Pt(), weight);
      hHpTzzdeltaPhicompparticle->Fill(h_particle.Pt(), zzdeltaPhiparticle, weight);
      hHpTzzdeltaEtacompparticle->Fill(h_particle.Pt(), zzdeltaEtaparticle, weight);                                       
      hbbzzdeltaPhicompparticle->Fill(bbdeltaPhiparticle, zzdeltaPhiparticle, weight);
      hbbzzdeltaEtacompparticle->Fill(bbdeltaEtaparticle, zzdeltaEtaparticle, weight);
    }

    if(switchVal_parton){
      hbbjjdeltaPhicompparton->Fill(bbdeltaPhiparton,jjdeltaPhiparton,weight);
      hbbdeltaEtajjdeltaPhicompparton->Fill(bbdeltaEtaparton,jjdeltaPhiparton,weight);
      hHpTl1l2deltaPhicompparton->Fill(h_parton.Pt(), l1l2deltaPhiparton, weight);
      hHpTl3l4deltaPhicompparton->Fill(h_parton.Pt(), l3l4deltaPhiparton, weight);
      hHz1pTcompparton->Fill(h_parton.Pt(), z1_parton.Pt(), weight);
      hHz2pTcompparton->Fill(h_parton.Pt(), z2_parton.Pt(), weight);
      hHzzpTcompparton->Fill(h_parton.Pt(), z1_parton.Pt() + z2_parton.Pt(), weight);
      hHpTzzdeltaPhicompparton->Fill(h_parton.Pt(), zzdeltaPhiparton, weight);
      hHpTzzdeltaEtacompparton->Fill(h_parton.Pt(), zzdeltaEtaparton, weight);                                         
      hbbzzdeltaPhicompparton->Fill(bbdeltaPhiparton, zzdeltaPhiparton, weight);
      hbbzzdeltaEtacompparton->Fill(bbdeltaEtaparton, zzdeltaEtaparton, weight);
    }

//------------------------------------------------------------------------------------------------------------------------------------------------------------
// FILL HISTOGRAMS - MISC
//------------------------------------------------------------------------------------------------------------------------------------------------------------


    // w ET- reco 
    if(switchVal_reco==0){
      if(thisRecoEventType==0 && wleps.size()>=2){
        hllpTET0reco->Fill((l1_reco+l2_reco).Pt(),weight);
        hllmET0reco->Fill((l1_reco+l2_reco).M(),weight);
      } else if(thisRecoEventType==1 && wleps.size()>=2){
        hllpTET1reco->Fill((l1_reco+l2_reco).Pt(),weight);
        hllmET1reco->Fill((l1_reco+l2_reco).M(),weight);
      } else if(thisRecoEventType==2 && wleps.size()>=2){
        hllpTET2reco->Fill((l1_reco+l2_reco).Pt(),weight);
        hllmET2reco->Fill((l1_reco+l2_reco).M(),weight);
      } else if(thisRecoEventType==3 && wleps.size()>=2){
        hllpTET3reco->Fill((l1_reco+l2_reco).Pt(),weight);
        hllmET3reco->Fill((l1_reco+l2_reco).M(),weight);
      }
    }


    // w ET- particle 
    if(switchVal_particle==0){
      if(thisParticleEventType==0 && WParticlePairIndices.size()>=2){
        hllpTET0particle->Fill((l1_particle+l2_particle).Pt(),weight);
        hllmET0particle->Fill((l1_particle+l2_particle).M(),weight);
      } else if(thisParticleEventType==1 && WParticlePairIndices.size()>=2){
        hllpTET1particle->Fill((l1_particle+l2_particle).Pt(),weight);
        hllmET1particle->Fill((l1_particle+l2_particle).M(),weight);
      } else if(thisParticleEventType==2 && WParticlePairIndices.size()>=2){
        hllpTET2particle->Fill((l1_particle+l2_particle).Pt(),weight);
        hllmET2particle->Fill((l1_particle+l2_particle).M(),weight);
      } else if(thisParticleEventType==3 && WParticlePairIndices.size()>=2){
        hllpTET3particle->Fill((l1_particle+l2_particle).Pt(),weight);
        hllmET3particle->Fill((l1_particle+l2_particle).M(),weight);
      }
    }
*/



/*
  //------------------------------------------------------------------------------------------------------------------------------------------------------------
  // ONNX - WORK IN PROGRESS
  //------------------------------------------------------------------------------------------------------------------------------------------------------------

    bool doONNX=true;

#ifdef ONNXRUN

    if(doONNX){
      // onnxruntime setup
    //string model_file="/Users/gaetano/Documents/universita/SnowMass2020/Analysis/brown-cern/higgsandmore/delphesAna/vvhjj/delphesModel.onnx";
    string model_file="/Users/gaetano/Documents/universita/SnowMass2020/Analysis/brown-cern/higgsandmore/delphesAna/vvhjj/delphesModel_changed.onnx";
    //string model_file="/Users/gaetano/Documents/universita/SnowMass2020/Analysis/brown-cern/higgsandmore/delphesAna/vvhjj/examples_sv_Jan.onnx";
    

    auto providers = Ort::GetAvailableProviders();
    for (auto provider : providers) {
      std::cout << provider << std::endl;
    }
    // cout<<endl;
    
    Ort::Env env = Ort::Env(OrtLoggingLevel::ORT_LOGGING_LEVEL_VERBOSE, "Default");
    Ort::SessionOptions sessionOptions;
    sessionOptions.SetIntraOpNumThreads(1);
    sessionOptions.SetGraphOptimizationLevel(GraphOptimizationLevel::ORT_DISABLE_ALL);
    // Optimization will take time and memory during startup
    //sessionOptions.SetGraphOptimizationLevel(GraphOptimizationLevel::ORT_DISABLE_ALL);
    
    Ort::Session session = Ort::Session(env, model_file.c_str(), sessionOptions);
    Ort::AllocatorWithDefaultOptions allocator;
    
    // Demonstration of getting input node info by code
    size_t num_input_nodes = 0;
    std::vector<const char*>* input_node_names = nullptr; // Input node names
    //std::vector<const char*>* output_node_names = new std::vector<const char*>();
    std::vector<const char*> output_node_names;
    std::vector<std::vector<int64_t>> input_node_dims;    // Input node dimension.
    ONNXTensorElementDataType type;                       // Used to print input info
    Ort::TypeInfo* type_info;
    
    num_input_nodes = session.GetInputCount();
    input_node_names = new std::vector<const char*>;
    for (int i = 0; i < num_input_nodes; i++) {
      
      char* tempstring = new char[strlen(session.GetInputNameAllocated(i, allocator).get()) + 1];
      snprintf(tempstring, strlen(session.GetInputNameAllocated(i, allocator).get()) + 1, session.GetInputNameAllocated(i, allocator).get());
      input_node_names->push_back(tempstring);
      type_info = new Ort::TypeInfo(session.GetInputTypeInfo(i));
      auto tensor_info = type_info->GetTensorTypeAndShapeInfo();
      cout<<"tensor info "<< tensor_info<<endl;
      cout<<"tensor size "<<VectorProduct(tensor_info.GetShape())<<endl;
      cout<<"tensor shape "<<tensor_info.GetShape()<<endl;
      ///cout<<"tensor element count "<<tensor_info.GetElementCount()<<endl;
      //cout<<"tensor dimensions count "<<tensor_info.GetDimensionsCount()<<endl;
      //star:vector<int64_t> input_dims
      
      
      
      type = tensor_info.GetElementType();
      input_node_dims.push_back(tensor_info.GetShape());

      //for (int j = 0; j < input_node_dims.size(); j++) {
      //if (input_node_dims[j] == -1)
      //{
      //  input_node_dims[j] = 1;
      //}
      //printf("Input %d : dim %d=%jd\n", i, j, input_node_dims[j]);
      //}

      // print input shapes/dims
      printf("Input %d : name=%s\n", i, input_node_names->back());
      printf("Input %d : num_dims=%zu\n", i, input_node_dims.back().size());
      for (int j = 0; j < input_node_dims.back().size(); j++)
	printf("Input %d : dim %d=%jd\n", i, j, input_node_dims.back()[j]);
      printf("Input %d : type=%d\n", i, type);
      
      delete(type_info);
    }
    
    
    
    // Set output node name explicitly
    output_node_names.push_back("output");

    cout<<"------"<<endl;
    size_t inputCount = session.GetInputCount();
    for (int i = 0; i < inputCount; ++i) {
        auto name = session.GetInputNameAllocated(i, allocator);
        auto shape = session.GetInputTypeInfo(i).GetTensorTypeAndShapeInfo().GetShape();

        std::cout << "Input Number: " << i << std::endl;
        std::cout << " Input Name: " << name.get() << std::endl;
        std::cout << " Input Shape: " << shape << std::endl;
    }

    size_t outputCount = session.GetOutputCount();
    for (int i = 0; i < outputCount; ++i) {
        auto name = session.GetOutputNameAllocated(i, allocator);
        auto shape = session.GetOutputTypeInfo(i).GetTensorTypeAndShapeInfo().GetShape();

        std::cout << "Output Number: " << i << std::endl;
        std::cout << " Output Name: " << name.get() << std::endl;
        std::cout << " Output Shape: " << shape << std::endl;
    }

    


    /*
    std::vector<float>* input_tensor_values;    // Raw input
    std::vector<Ort::Value> inputTensor;        // Onnxruntime allowed input
    
    // this will make the input into 1,3,640,640
    cv::Mat blob = cv::dnn::blobFromImage(image, 1 / 255.0, cv::Size(640, 640), (0, 0, 0), false, false);
    size_t input_tensor_size = blob.total();
    input_tensor_values = new std::vector<float>((float*)blob.data, (float*)blob.data + input_tensor_size);
    
    try {
      inputTensor.emplace_back(Ort::Value::CreateTensor<float>(memory_info, input_tensor_values->data(), input_tensor_size, input_node_dims[0].data(), input_node_dims[0].size()));
    }
    catch (Ort::Exception oe) {
      std::cout << "ONNX exception caught: " << oe.what() << ". Code: " << oe.GetOrtErrorCode() << ".\n";
      return -1;
      }
    
    
    
    cout<<"HERE"<<endl;
  
   
    cout << "KABOOM "<<endl;
    
    auto memoryInfo = Ort::MemoryInfo::CreateCpu(OrtAllocatorType::OrtDeviceAllocator, OrtMemType::OrtMemTypeCPUOutput);
    

    
      std::cout << "Start warming up" << endl;
     
     
     
      std::vector<Ort::Value> input_tensors;
      std::vector<Ort::Value> output_tensors;
      std::cout << "################### befor run:##############" << endl;
      //std::cout << "input node name:" << inputNodeNames[0] << endl;
      //std::cout << "output0 node name:" << outputNodeNames[0] << endl;
      for (int i = 0; i < num_input_nodes; i++) {

	auto name = session.GetInputNameAllocated(i, allocator);
        auto shape = session.GetInputTypeInfo(i).GetTensorTypeAndShapeInfo().GetShape();
	size_t input_tensor_length = VectorProduct(shape);
	cout<<"Tensor size "<<input_tensor_length<<endl;
	float temp[input_tensor_length];
	
	type_info = new Ort::TypeInfo(session.GetInputTypeInfo(i));
	auto tensor_info = type_info->GetTensorTypeAndShapeInfo();
	cout<<"tensor info "<< tensor_info<<endl;
	cout<<"tensor size "<<VectorProduct(tensor_info.GetShape())<<endl;
	cout<<"tensor shape "<<tensor_info.GetShape()<<endl;

	input_tensors.push_back(Ort::Value::CreateTensor<float>(
							      memoryInfo, temp, input_tensor_length, tensor_info.GetShape().data(),
							      tensor_info.GetShape().size()));
      	
    }

      
      
      //input_tensors.push_back(Ort::Value::CreateTensor<float>(
      //						      memoryInfo, temp, input_tensor_length, input_tensor_info.GetShape().data(),
      //						      input_tensor_info.GetShape().size()));

      cout<<" Loop "<<endl;

      //const int64_t shape=3; //inputTensorShape.data()
      //input_tensors.push_back(Ort::Value::CreateTensor<float>(
      //						       memoryInfo, temp, 1,&shape,
      //						      1));
      
      //for (int i = 0; i < 1; i++) {
      //output_tensors = session.Run(Ort::RunOptions{ nullptr },
      //			     inputNodeNames.data(),
      //			     input_tensors.data(),
      //			     inputNodeNames.size(),
      //			     outputNodeNames.data(),
      //			     outputNodeNames.size());
      //}
      //std::cout << "################### after run:##############" << endl;
      //std::cout << "input node name:" << inputNodeNames[0] << endl;
      //std::cout << "output0 node name:" << outputNodeNames[0] << endl;
      //std::cout << "output1 node name:" << outputNodeNames[1] << endl;
    
    
      std::cout << "*********************************** test onnx ok  ***************************************" << endl;
    }

#endif
*/