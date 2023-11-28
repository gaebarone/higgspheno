#ifdef __CLING__
R__LOAD_LIBRARY("libDelphes")
#endif

#include "classes/DelphesClasses.h"
#include "classes/DelphesLHEFReader.h"
#include "external/ExRootAnalysis/ExRootTreeReader.h"
#include "../common_includes/ghost_tagging.h"
#include "../common_includes/get_cross_section.h"
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
#include "../common_includes/combinations.h" 
#include <iomanip>

typedef std::map<std::string, std::pair<int,double>> cutFlowMapDef;

using namespace std;

double btagEff=0.85;
double fakeEff=0.01;

//------------------------------------------------------------------------------------------------------------------------------------------------------------
// initialize combinations
//------------------------------------------------------------------------------------------------------------------------------------------------------------

vector <string> types={"4mu","4e", "2mu2e", "2e2mu"};
std::vector <string> cutList_reco={"initial reco", "1 btag reco", "2 good j reco", "2 b-like jet pairs reco", "found bb reco", "2 vbfj reco", "vbfj pairs","2.5 deltaEta vbf reco","OSFL"};
std::vector <string> cutList_particle={"initial particle", "1 btag particle", "2 good j particle", "2 b-like jet pairs part", "found bb particle", "2 vbfj particle", "comb vbf part","2.5 deltaEta vbf particle","OSFL"};
std::vector <string> cutList_parton={"initial parton", "Higgs Candidate", "ZZ parton"};


void FillCutFlow(TH1F* hSel, TProfile *hEff,std::map<string, std::pair<int,double>> cutFlowMap, std::vector <string> cutList, string label){
  //  cout<<" Filling cutflow "<<label<<endl;
    
  for(int i=0; i<(int) cutList.size(); i++) {

    const std::string cutName = cutList[i];
    double passed_reco =  cutFlowMap[cutName].second;
    double efficiency_reco = 100.00 * cutFlowMap[cutName].second / cutFlowMap[cutList[0]].second;
    //cout<<"Bin "<<i+1<<" passed "<<passed_reco<<" eff "<<efficiency_reco<<endl;
    hSel->GetXaxis()->SetBinLabel(i+1,cutName.c_str());
    hEff->GetXaxis()->SetBinLabel(i+1,cutName.c_str());
  
    
    hSel->SetBinContent(i+1,passed_reco);
    hEff->Fill(i+1.0,efficiency_reco);
  } 
} 

void PrintCutFlow(std::map<string, std::pair<int,double>> cutFlowMap, std::vector <string> cutList, string label){

std::cout<<std::left<<std::setw(25)<<label<<" Cut"<<std::setw(10)<<label<<" Passed"<<std::setw(15)<<" Rel Eff "<< std::setw(15)<<label <<" Efficiency"<< std::endl;
    for(int i=0; i<(int) cutList.size(); i++) {
        const std::string cutName = cutList[i];
        double passed_reco =  cutFlowMap[cutName].first; // switch to second if one wants to use weighted events!
        double efficiency_reco = 100.00 * cutFlowMap[cutName].second / cutFlowMap[cutList[0]].second;

	double relEff= 100;;
	if(i>0)
	  relEff = 100.00 * cutFlowMap[cutName].second / cutFlowMap[cutList.at(i-1)].second; 
	std::cout<<std::left<<std::setw(25)<<cutName<<std::setw(25)<<passed_reco<< std::setw(25)<<relEff <<std::setw(25)<<efficiency_reco<< std::endl;
    }

    cout<<endl;
}


void increaseCount(std::map<string, std::pair<int,double>> & cutFlowMap, string cutName, double weight){
  cutFlowMap[cutName]=make_pair(cutFlowMap[cutName].first+1,cutFlowMap[cutName].second+weight);
}

void remove_overlaps(vector< pair<int,int>> muPairIndices){
  for( vector< pair<int,int>>::iterator it=muPairIndices.begin(); it!=muPairIndices.end(); it++){
    pair<int,int> one=(*it);
    for(vector< pair<int,int>>::iterator it2=it+1; it2!=muPairIndices.end(); it2++){
      pair<int,int> two=(*it2);
      if( one.first == two.first || one.first == two.second || one.second == two.second || one.second == two.second) muPairIndices.erase(it2--);
    }
  }
}

bool FillHiggsTruthRecord(TClonesArray *branchGenParticle, TLorentzVector & h_parton, TLorentzVector & b1_parton, TLorentzVector & b2_parton, TLorentzVector & j1_parton, TLorentzVector &j2_parton){
    
    bool foundHiggsOrinTree=false;
    bool foundHiggsToBbb=false;
    bool foundHiggsToBbbJets=false;
    GenParticle *daughter1=nullptr;
    GenParticle *daughter2=nullptr;

    for(int i=0; i<(int)branchGenParticle->GetEntries(); i++) {
      GenParticle *particle=(GenParticle*) branchGenParticle->At(i);
      int d1_pid = 9999;
      if (particle->D1 != -1) {
        daughter1 = (GenParticle*) branchGenParticle->At(particle->D1);
        d1_pid = daughter1 -> PID;
      }
      int d2_pid = 9999;
      if (particle->D2 != -1) {
        daughter2 = (GenParticle*) branchGenParticle->At(particle->D2);
        d2_pid = daughter2 -> PID;
      } 
      // higgs parton
      if (particle->PID == 25 && d1_pid != 25 && d2_pid != 25) {
        h_parton = particle->P4();
	foundHiggsOrinTree=true;
	
        // check for b parton children
        if (abs(d1_pid) == 5 && abs(d2_pid) == 5) {
	  foundHiggsToBbb=true; 
          if (daughter1 -> PT > daughter2 -> PT) {
            b1_parton = daughter1 -> P4();
            b2_parton = daughter2 -> P4();
          } else {
            b1_parton = daughter2 -> P4();
            b2_parton = daughter1 -> P4();
          }
        }
      } else if ((particle->PID == 1 || particle->PID == 2 || particle->PID == 3 || particle->PID == 4 || particle->PID == 6) && d1_pid != particle->PID && d2_pid != particle->PID){
	foundHiggsToBbbJets=true;
	if (daughter1 -> PT > daughter2 -> PT) {
	  j1_parton = daughter1 -> P4();
	  j2_parton = daughter2 -> P4();
	} else {
	  j1_parton = daughter2 -> P4();
	  j2_parton = daughter1 -> P4();
	}
      }
    }
    return (foundHiggsOrinTree && foundHiggsToBbb && foundHiggsToBbbJets); 
  }


std::vector<std::vector<int>> combinationsNoRepetitionAndOrderDoesNotMatter (int subsetSize, std::vector<int> setOfNumbers){
  std::vector<std::vector<int> > subsets{};
  subsets.reserve (count_each_combination (setOfNumbers.begin (), setOfNumbers.begin () + subsetSize, setOfNumbers.end ()));
  for_each_combination (setOfNumbers.begin (), setOfNumbers.begin () + subsetSize, setOfNumbers.end (), [&subsets] (auto first, auto last) {
    subsets.push_back (std::vector<int>{ first, last });
    return false;
  });
  return subsets;
}

std::vector <std::vector <int>> getAllCombinations(vector<int> inputVector, int k){
  std::vector<vector<int>> combinations; 
  std::vector<int> selector(inputVector.size());
    std::fill(selector.begin(), selector.begin() + k, 1);

    do {
        std::vector<int> selectedIds;
        std::vector<int> selectedVectorElements;
        for (int i = 0; i < inputVector.size(); i++) {
            if (selector[i]) {
                selectedIds.push_back(i);
            }
        }
        for (auto& id : selectedIds) {
            selectedVectorElements.push_back(inputVector[id]);
        }
        combinations.push_back(selectedVectorElements);
    } while (std::prev_permutation(selector.begin(), selector.end()));

    return combinations;
}

//------------------------------------------------------------------------------------------------------------------------------------------------------------
// histogram settings
//------------------------------------------------------------------------------------------------------------------------------------------------------------

void PrintCanvas(TCanvas *c=nullptr,string name="default"){
  std::vector <string> types={"jpg"}; 
  for(std::vector<string>::iterator it=types.begin(); it!=types.end(); it++) {
    c->Print(Form("jpg/%s.%s",name.c_str(),(*it).c_str()),(*it).c_str());
  }
}

void draw_hist(TH1 *histo, const char *name, const char *title, const char *axistitle) {
  TCanvas *c = new TCanvas(name, title, 1500,1200);
  c->cd();
  histo->GetXaxis()->SetTitle(axistitle);
  histo->SetMinimum(0.0);
  histo->Draw("hist e");
  PrintCanvas(c, name);
}

void draw_hist2(TH2F*histo, const char *name, const char *title, const char *xaxistitle, const char *yaxistitle) {
  TCanvas *c = new TCanvas(name, title, 1500, 1200);
  histo->GetXaxis()->SetTitle(xaxistitle);
  histo->GetYaxis()->SetTitle(yaxistitle);
  //histo->SetStatX(0.875);
  //histo->SetStatY(0.875);
  histo->Draw("COLZ");
  PrintCanvas(c, name);
}


GenParticle* find_parent(TClonesArray *branchGenParticle, GenParticle *particle, int target) {
  
  // check particle itself
  if (particle -> PID != target){
    //cout<<endl;
    //if( particle != nullptr) cout<<"Done first Parent is "<<particle->PID<<endl;
    return particle;
  }
  
  
  // safely access daughters
  int d1_pid = 9999;
  int d2_pid = 9999;
  GenParticle *parent1;
  GenParticle *parent2;
  if (particle->M1 != -1) {
    parent1 = (GenParticle*) branchGenParticle->At(particle->M1);
    d1_pid = parent1 -> PID;
  }
  if (particle->M2 != -1) {
    parent2 = (GenParticle*) branchGenParticle->At(particle->M2);
    d2_pid = parent2 -> PID;
  }
  //cout<<"Particle "<<particle->PID<<" parent 1 "<<parent1->PID<<" 2 "<<parent2->PID<<endl;

  if (d1_pid == target ) {
    //cout<<"-> same particle  for parent 1 "<<endl;
    return find_parent(branchGenParticle, parent1,target);
  }
  else
    return parent1;
  /*
  // recursive call on parents or return null
  if ( d1_pid != target && d2_pid==target){
    //cout<<"Done Parent 1 is "<<parent1->PID<<endl;
    return parent1; 
  }
  else if( d2_pid !=target && d1_pid==target){
    //cout<<"Done Parent 2 is "<<parent2->PID<<endl;
      return parent2;
  }
  else if (d1_pid == target ) {
    //cout<<"-> same particle  for parent 1 "<<endl;
    return find_parent(branchGenParticle, parent1,target);
  }
  else if (d2_pid == target) {
    //cout<<"-> same particle  for parent 2 "<<endl;
    return find_parent(branchGenParticle, parent2, target);
  }
  else{
    //cout<<"Done Parent is "<<particle->PID<<endl;
    //cout<<endl;
    return particle;
  }
  */
}



//------------------------------------------------------------------------------------------------------------------------------------------------------------
// define find_status1_child
//------------------------------------------------------------------------------------------------------------------------------------------------------------

GenParticle* find_status1_child(TClonesArray *branchGenParticle, GenParticle *particle, int target) {
  // check particle itself
  if (abs(particle -> PID) == target && particle -> Status == 1) {
    return particle;
  }
  // safely access daughters
  int d1_pid = 9999;
  int d2_pid = 9999;
  GenParticle *daughter1;
  GenParticle *daughter2;
  if (particle->D1 != -1) {
    daughter1 = (GenParticle*) branchGenParticle->At(particle->D1);
    d1_pid = daughter1 -> PID;
  }
  if (particle->D2 != -1) {
    daughter2 = (GenParticle*) branchGenParticle->At(particle->D2);
    d2_pid = daughter2 -> PID;
  }
  // recursive call on daughters or return null
  if (abs(d1_pid) == target) {
    return find_status1_child(branchGenParticle, daughter1, target);
  } else if (abs(d2_pid) == target) {
    return find_status1_child(branchGenParticle, daughter2, target);
  } else {
    return nullptr;
  }
}

//------------------------------------------------------------------------------------------------------------------------------------------------------------
// num_entries
//------------------------------------------------------------------------------------------------------------------------------------------------------------

Long64_t get_num_entries(const char *inputName) {
  TChain chain("Delphes");
  chain.Add(inputName);
  ExRootTreeReader *treeReader = new ExRootTreeReader(&chain);
  return treeReader->GetEntries();
}

Long64_t get_total_num_entries(const char *process_name) {
  std::string inputFileName = std::string(process_name) + "_inputs.txt";
  std::ifstream inputFile(inputFileName.c_str());
  std::string line;
  TChain chain("Delphes");
  Long64_t total = 0;
  while (std::getline(inputFile, line)) {
    total += get_num_entries(line.c_str());
  }
  return total;
}



//------------------------------------------------------------------------------------------------------------------------------------------------------------
// z analyzer
//------------------------------------------------------------------------------------------------------------------------------------------------------------

// void zAnalyzer(const char *inputFile,const char *outputFile, int kappaVal = 8) {
void zAnalyzer(const char *inputFile,const char *outputFile, int analysisType=0) {

#ifdef __CLING__
  gSystem->Load("libDelphes");
#endif
 // const double cross_section = get_cross_section(process_name);
  TChain chain("Delphes");
  chain.Add(inputFile);

  ExRootTreeReader *treeReader = new ExRootTreeReader(&chain);
  Long64_t numberOfEntries = treeReader->GetEntries();

  TClonesArray *branchJet = treeReader->UseBranch("Jet");
  TClonesArray *branchElectron = treeReader->UseBranch("Electron");
  TClonesArray *branchMuon = treeReader->UseBranch("Muon");
  TClonesArray *branchEvent = treeReader->UseBranch("Event");
  TClonesArray *branchGenParticle = treeReader->UseBranch("Particle");
  TClonesArray *branchGenJet = treeReader->UseBranch("GenJet");
  TClonesArray *branchWeight   = treeReader->UseBranch("Weight");

  bool fill_1D = true;
  bool fill_2D = true;

  const double mBins = 20;
  const double pTBins = 20;
  const double phiBins = 10;
  const double etaBins = 10;
  const double RBins = 10;
  const double cosBins = 10;

  const double hpTmin = 0;
  const double jpTmin = 0;
  const double zpTmin = 0;
  const double lpTmin = 0;
  const double hpTmax = 500;
  const double jpTmax = 750;
  const double zpTmax = 500;
  const double lpTmax = 750;

  const double hmmin = 0;
  const double jmmin = 0;
  const double zmmin = 0;
  const double hmmax = 200;
  const double jmmax = 200;
  const double zmmax = 200;

  const double hetamin = 0;
  const double jetamin = 0;
  const double zetamin = 0;
  const double letamin = 0;
  const double hetamax = 5;
  const double jetamax = 5;
  const double zetamax = 5;
  const double letamax = 5;

  const double hRmin = 0;
  const double jRmin = 0;
  const double zRmin = 0;
  const double lRmin = 0;
  const double hRmax = 5;
  const double jRmax = 5;
  const double zRmax = 5;
  const double lRmax = 5;

  vector <TH1F*> listOfTH1;
  vector <TH2F*> listOfTH2;
  vector <TProfile*> listOfTPorifles;
  // cutflows
  int cutVal_reco = 0;
double cutValW_reco = 0;
 
//std::vector <string> cutList_reco={"initial reco"};
std::map<string, std::pair<int,double>> cutFlowMap_reco;
for(int i=0; i<(int) cutList_reco.size(); i++) { 
  cutFlowMap_reco[cutList_reco.at(i)] = make_pair(0,0.0); 
 }
 
 int cutVal_particle = 0;
 double cutValW_particle = 0;
 
 //std::vector <string> cutList_particle={"initial particle"};
 std::map<string, std::pair<int,double>> cutFlowMap_particle;
 for(int i=0; i<(int) cutList_particle.size(); i++) { 
   cutFlowMap_particle[cutList_particle.at(i)] = make_pair(0,0.0); 
 }
 
 int cutVal_parton = 0;
 double cutValW_parton = 0;
 std::map<string, std::pair<int,double>> cutFlowMap_parton;
 for(int i=0; i<(int) cutList_parton.size(); i++) { 
   cutFlowMap_parton[cutList_parton.at(i)] = make_pair(0,0.0); 
 }

 
 vector <string> selType={"reco","particle","parton"};
 std::map<string, vector<string>> cutFlowMByType;
 cutFlowMByType["reco"]=cutList_reco;
 cutFlowMByType["particle"]=cutList_particle;
 cutFlowMByType["parton"]=cutList_parton;
 
  
  std::map<string,TH1F*> cutFlowHists;
  std::map<string,TProfile*> cutFlowEffs;

  std::map<string, cutFlowMapDef* > cutFlowMapAll;
  cutFlowMapAll["reco"] =  & cutFlowMap_reco;
  cutFlowMapAll["particle"] = & cutFlowMap_particle;
  cutFlowMapAll["parton"] = & cutFlowMap_parton;
  
  
  
  for(std::vector<string>::iterator it=selType.begin(); it!=selType.end(); it++){
    cutFlowHists[(*it)]=new TH1F(Form("hSel_%s",(*it).c_str()),"",cutFlowMByType[(*it)].size(),0,cutFlowMByType[(*it)].size()+1);
    cutFlowEffs[(*it)]=new TProfile(Form("hEff_%s",(*it).c_str()),"",cutFlowMByType[(*it)].size(),0,cutFlowMByType[(*it)].size()+1);
    
    listOfTH1.push_back(cutFlowHists[(*it)]);
    listOfTPorifles.push_back((cutFlowEffs[(*it)]));
  }
   
  
  
  // higgs
  // 1D
  // reco 
  // pT + m
  TH1F *hHpTreco = new TH1F("hbb_pT_reco", "p^{T}_{hbb}_reco", pTBins, hpTmin, hpTmax); 	listOfTH1.push_back(hHpTreco);
  TH1F *hHmreco = new TH1F("hbb_m_reco", "m_{hbb}_reco", mBins, hmmin, hmmax); listOfTH1.push_back(hHmreco);
  TH1F *hb1pTreco = new TH1F("b1_pT_reco", "p^{T}_{b1}_reco", pTBins, hpTmin, hpTmax); listOfTH1.push_back(hb1pTreco);
  TH1F *hb1mreco = new TH1F("b1_m_reco", "m_{b1}_reco", mBins, hmmin, hmmax);listOfTH1.push_back(hb1mreco);
  TH1F *hb2pTreco = new TH1F("b2_pT_reco", "p^{T}_{b2}_reco", pTBins, hpTmin, hpTmax); listOfTH1.push_back(hb2pTreco);
  TH1F *hb2mreco = new TH1F("b2_m_reco", "m_{b2}_reco", mBins, hmmin, hmmax); listOfTH1.push_back(hb2mreco);
  // phi
  TH1F *hHphireco = new TH1F("hbb_#phi_reco", "#phi_{hbb}_reco", phiBins, -TMath::Pi(), +TMath::Pi()); listOfTH1.push_back(hHphireco);
  TH1F *hb1phireco = new TH1F("b1_#phi_reco", "#phi_{b1}_reco", phiBins, -TMath::Pi(), +TMath::Pi()); listOfTH1.push_back(hb1phireco);
  TH1F *hb2phireco = new TH1F("b2_#phi_reco", "#phi_{b2}_reco", phiBins, -TMath::Pi(), +TMath::Pi()); listOfTH1.push_back(hb2phireco);
  TH1F *hbbdeltaPhireco = new TH1F("bb_#Delta#phi_reco", "#Delta#phi_{bb}_reco", phiBins, -TMath::Pi(), +TMath::Pi()); listOfTH1.push_back(hbbdeltaPhireco);
  // eta
  TH1F *hHetareco = new TH1F("hbb_#eta_reco", "#eta_{hbb}_reco", etaBins, hetamin, hetamax); listOfTH1.push_back(hHetareco);
  TH1F *hb1etareco = new TH1F("b1_#eta_reco", "#eta_{b1}_reco", etaBins, hetamin, hetamax); listOfTH1.push_back(hb1etareco);
  TH1F *hb2etareco = new TH1F("b2_#eta_reco", "#eta_{b2}_reco", etaBins, hetamin, hetamax); listOfTH1.push_back(hb2etareco);
  TH1F *hbbdeltaEtareco = new TH1F("bb_#Delta#eta_reco", "#Delta#eta_{bb}_reco", etaBins, hetamin, hetamax); listOfTH1.push_back(hbbdeltaEtareco);
  // R
  TH1F *hHRreco = new TH1F("hbb_R_reco", "R_{hbb}_reco", RBins, hRmin, hRmax); listOfTH1.push_back(hHRreco);
  TH1F *hb1Rreco = new TH1F("b1_R_reco", "R_{b1}_reco", RBins, hRmin, hRmax); listOfTH1.push_back(hb1Rreco);
  TH1F *hb2Rreco = new TH1F("b2_R_reco", "R_{b2}_reco", RBins, hRmin, hRmax); listOfTH1.push_back(hb2Rreco);
  TH1F *hbbdeltaRreco = new TH1F("bb_#DeltaR_reco", "#DeltaR_{bb}_reco", RBins, hRmin, hRmax); listOfTH1.push_back(hbbdeltaRreco);
  
  
  // particle
  // pT + m
  TH1F *hHpTparticle = new TH1F("hbb_pT_particle", "p^{T}_{hbb}_particle", pTBins, hpTmin, hpTmax); listOfTH1.push_back(hHpTparticle);
  TH1F *hHmparticle = new TH1F("hbb_m_particle", "m_{hbb}_particle", mBins, hmmin, hmmax); listOfTH1.push_back(hHmparticle);
  TH1F *hb1pTparticle = new TH1F("b1_pT_particle", "p^{T}_{b1}_particle", pTBins, hpTmin, hpTmax); listOfTH1.push_back(hb1pTparticle);
  TH1F *hb1mparticle = new TH1F("b1_m_particle", "m_{b1}_particle", mBins, hmmin, hmmax); listOfTH1.push_back(hb1mparticle);
  TH1F *hb2pTparticle = new TH1F("b2_pT_particle", "p^{T}_{b2}_particle", pTBins, hpTmin, hpTmax); listOfTH1.push_back(hb2pTparticle);
  TH1F *hb2mparticle = new TH1F("b2_m_particle", "m_{b2}_particle", mBins, hmmin, hmmax); listOfTH1.push_back(hb2mparticle);
  // phi
  TH1F *hHphiparticle = new TH1F("hbb_#phi_particle", "#phi_{hbb}_particle", phiBins, -TMath::Pi(), +TMath::Pi()); listOfTH1.push_back(hHphiparticle);
  TH1F *hb1phiparticle = new TH1F("b1_#phi_particle", "#phi_{b1}_particle", phiBins, -TMath::Pi(), +TMath::Pi()); listOfTH1.push_back(hb1phiparticle);
  TH1F *hb2phiparticle = new TH1F("b2_#phi_particle", "#phi_{b2}_particle", phiBins, -TMath::Pi(), +TMath::Pi()); listOfTH1.push_back(hb2phiparticle);
  TH1F *hbbdeltaPhiparticle = new TH1F("bb_#Delta#phi_particle", "#Delta#phi_{bb}_particle", phiBins, -TMath::Pi(), +TMath::Pi()); listOfTH1.push_back(hbbdeltaPhiparticle);
  // eta
  TH1F *hHetaparticle = new TH1F("hbb_#eta_particle", "#eta_{hbb}_particle", etaBins, hetamin, hetamax); listOfTH1.push_back(hHetaparticle);
  TH1F *hb1etaparticle = new TH1F("b1_#eta_particle", "#eta_{b1}_particle", etaBins, hetamin, hetamax); listOfTH1.push_back(hb1etaparticle);
  TH1F *hb2etaparticle = new TH1F("b2_#eta_particle", "#eta_{b2}_particle", etaBins, hetamin, hetamax); listOfTH1.push_back(hb2etaparticle);
  TH1F *hbbdeltaEtaparticle = new TH1F("bb_#Delta#eta_particle", "#Delta#eta_{bb}_particle", etaBins, hetamin, hetamax); listOfTH1.push_back(hbbdeltaEtaparticle);
        // R
  TH1F *hHRparticle = new TH1F("hbb_R_particle", "R_{hbb}_particle", RBins, hRmin, hRmax); listOfTH1.push_back(hHRparticle);
  TH1F *hb1Rparticle = new TH1F("b1_R_particle", "R_{b1}_particle", RBins, hRmin, hRmax); listOfTH1.push_back(hb1Rparticle);
  TH1F *hb2Rparticle = new TH1F("b2_R_particle", "R_{b2}_particle", RBins, hRmin, hRmax); listOfTH1.push_back(hb2Rparticle);
  TH1F *hbbdeltaRparticle = new TH1F("bb_#DeltaR_particle", "#DeltaR_{bb}_particle", RBins, hRmin, hRmax); listOfTH1.push_back(hbbdeltaRparticle);
  
  // parton
  // pT + m
  TH1F *hHpTparton = new TH1F("hbb_pT_parton", "p^{T}_{hbb}_parton", pTBins, hpTmin, hpTmax); listOfTH1.push_back(hHpTparton);
  TH1F *hHmparton = new TH1F("hbb_m_parton", "m_{hbb}_parton", mBins, hmmin,  hmmax); listOfTH1.push_back(hHmparton);
  TH1F *hb1pTparton = new TH1F("b1_pT_parton", "p^{T}_{b1}_parton", pTBins, hpTmin, hpTmax); listOfTH1.push_back(hb1pTparton);
  TH1F *hb1mparton = new TH1F("b1_m_parton", "m_{b1}_parton", mBins, hmmin, hmmax); listOfTH1.push_back(hb1mparton);
  TH1F *hb2pTparton = new TH1F("b2_pT_parton", "p^{T}_{b2}_parton", pTBins, hpTmin, hpTmax); listOfTH1.push_back(hb2pTparton);
  TH1F *hb2mparton = new TH1F("b2_m_parton", "m_{b2}_parton", mBins, hmmin, hmmax); listOfTH1.push_back(hb2mparton);
  // phi
  TH1F *hHphiparton = new TH1F("hbb_#phi_parton", "#phi_{hbb}_parton", phiBins, -TMath::Pi(), +TMath::Pi()); listOfTH1.push_back(hHphiparton);
  TH1F *hb1phiparton = new TH1F("b1_#phi_parton", "#phi_{b1}_parton", phiBins, -TMath::Pi(), +TMath::Pi()); listOfTH1.push_back(hb1phiparton);
  TH1F *hb2phiparton = new TH1F("b2_#phi_parton", "#phi_{b2}_parton", phiBins, -TMath::Pi(), +TMath::Pi()); listOfTH1.push_back(hb2phiparton);
  TH1F *hbbdeltaPhiparton = new TH1F("bb_#Delta#phi_parton", "#Delta#phi_{bb}_parton", phiBins, -TMath::Pi(), +TMath::Pi()); listOfTH1.push_back(hbbdeltaPhiparton);
  // eta
  TH1F *hHetaparton = new TH1F("hbb_#eta_parton", "#eta_{hbb}_parton", etaBins, hetamin, hetamax); listOfTH1.push_back(hHetaparton);
  TH1F *hb1etaparton = new TH1F("b1_#eta_parton", "#eta_{b1}_parton", etaBins, hetamin, hetamax); listOfTH1.push_back(hb1etaparton);
  TH1F *hb2etaparton = new TH1F("b2_#eta_parton", "#eta_{b2}_parton", etaBins, hetamin, hetamax); listOfTH1.push_back(hb2etaparton);
  TH1F *hbbdeltaEtaparton = new TH1F("bb_#Delta#eta_parton", "#Delta#eta_{bb}_parton", etaBins, hetamin, hetamax); listOfTH1.push_back(hbbdeltaEtaparton);
  // R
  TH1F *hHRparton = new TH1F("hbb_R_parton", "R_{hbb}_parton", RBins, hRmin, hRmax); listOfTH1.push_back(hHRparton);
  TH1F *hb1Rparton = new TH1F("b1_R_parton", "R_{b1}_parton", RBins, hRmin, hRmax); listOfTH1.push_back(hb1Rparton);
  TH1F *hb2Rparton = new TH1F("b2_R_parton", "R_{b2}_parton", RBins, hRmin, hRmax); listOfTH1.push_back(hb2Rparton);
  TH1F *hbbdeltaRparton = new TH1F("bb_#DeltaR_parton", "#DeltaR_{bb}_parton", RBins, hRmin, hRmax); listOfTH1.push_back(hbbdeltaRparton);
  
  // jets
  // reco 
  // pT + m
  TH1F *hjjpTreco = new TH1F("jj_pT_reco", "p^{T}_{jj}_reco", pTBins, jpTmin, jpTmax); listOfTH1.push_back(hjjpTreco);
  TH1F *hj1pTreco = new TH1F("j1_pT_reco", "p^{T}_{j1}_reco", pTBins, jpTmin, jpTmax); listOfTH1.push_back(hj1pTreco);
  TH1F *hj2pTreco = new TH1F("j2_pT_reco", "p^{T}_{j2}_reco", pTBins, jpTmin, jpTmax); listOfTH1.push_back(hj2pTreco);
  // phi
  TH1F *hj1phireco = new TH1F("j1_#phi_reco", "#phi_{j1}_reco", phiBins, -TMath::Pi(), +TMath::Pi()); listOfTH1.push_back(hj1phireco);
  TH1F *hj2phireco = new TH1F("j2_#phi_reco", "#phi_{j2}_reco", phiBins, -TMath::Pi(), +TMath::Pi()); listOfTH1.push_back(hj2phireco);
  TH1F *hjjdeltaPhireco = new TH1F("jj_#Delta#phi_reco", "#Delta#phi_{jj}_reco", phiBins, -TMath::Pi(), +TMath::Pi()); listOfTH1.push_back(hjjdeltaPhireco);
  // eta
  TH1F *hj1etareco = new TH1F("j1_#eta_reco", "#eta_{j1}_reco", etaBins, jetamin, jetamax); listOfTH1.push_back(hj1etareco);
  TH1F *hj2etareco = new TH1F("j2_#eta_reco", "#eta_{j2}_reco", etaBins, jetamin, jetamax); listOfTH1.push_back(hj2etareco);
  TH1F *hjjdeltaEtareco = new TH1F("jj_#Delta#eta_reco", "#Delta#eta_{jj}_reco", etaBins, jetamin, jetamax); listOfTH1.push_back(hjjdeltaEtareco);
  // R
  TH1F *hj1Rreco = new TH1F("j1_R_reco", "R_{j1}_reco", RBins, jRmin, jRmax); listOfTH1.push_back(hj1Rreco);
  TH1F *hj2Rreco = new TH1F("j2_R_reco", "R_{j2}_reco", RBins, jRmin, jRmax); listOfTH1.push_back(hj2Rreco);
  TH1F *hjjdeltaRreco = new TH1F("jj_#DeltaR_reco", "#DeltaR_{jj}_reco", RBins, jRmin, jRmax); listOfTH1.push_back(hjjdeltaRreco);
  
  // particle
  // pT + m
  TH1F *hjjpTparticle = new TH1F("jj_pT_particle", "p^{T}_{jj}_particle", pTBins, jpTmin, jpTmax); listOfTH1.push_back(hjjpTparticle);
  TH1F *hj1pTparticle = new TH1F("j1_pT_particle", "p^{T}_{j1}_particle", pTBins, jpTmin, jpTmax); listOfTH1.push_back(hj1pTparticle);
  TH1F *hj2pTparticle = new TH1F("j2_pT_particle", "p^{T}_{j2}_particle", pTBins, jpTmin, jpTmax); listOfTH1.push_back(hj2pTparticle);
  // phi
  TH1F *hj1phiparticle = new TH1F("j1_#phi_particle", "#phi_{j1}_particle", phiBins, -TMath::Pi(), +TMath::Pi()); listOfTH1.push_back(hj1phiparticle);
  TH1F *hj2phiparticle = new TH1F("j2_#phi_particle", "#phi_{j2}_particle", phiBins, -TMath::Pi(), +TMath::Pi()); listOfTH1.push_back(hj2phiparticle);
  TH1F *hjjdeltaPhiparticle = new TH1F("jj_#Delta#phi_particle", "#Delta#phi_{jj}_particle", phiBins, -TMath::Pi(), +TMath::Pi()); listOfTH1.push_back(hjjdeltaPhiparticle);
  // eta
  TH1F *hj1etaparticle = new TH1F("j1_#eta_particle", "#eta_{j1}_particle", etaBins, jetamin, jetamax); listOfTH1.push_back(hj1etaparticle);
        TH1F *hj2etaparticle = new TH1F("j2_#eta_particle", "#eta_{j2}_particle", etaBins, jetamin, jetamax); listOfTH1.push_back(hj2etaparticle);
        TH1F *hjjdeltaEtaparticle = new TH1F("jj_#Delta#eta_particle", "#Delta#eta_{jj}_particle", etaBins, jetamin, jetamax); listOfTH1.push_back(hjjdeltaEtaparticle);
        // R
        TH1F *hj1Rparticle = new TH1F("j1_R_particle", "R_{j1}_particle", RBins, jRmin, jRmax); listOfTH1.push_back(hj1Rparticle);
        TH1F *hj2Rparticle = new TH1F("j2_R_particle", "R_{j2}_particle", RBins, jRmin, jRmax); listOfTH1.push_back(hj2Rparticle);
        TH1F *hjjdeltaRparticle = new TH1F("jj_#DeltaR_particle", "#DeltaR_{jj}_particle", RBins, jRmin, jRmax); listOfTH1.push_back(hjjdeltaRparticle);

// z 
    // reco
        // pT + m
        TH1F *hz1pTreco = new TH1F("z1_pT_reco", "p^{T}_{z1}_reco", pTBins, zpTmin, zpTmax); listOfTH1.push_back(hz1pTreco);
        TH1F *hz2pTreco = new TH1F("z2_pT_reco", "p^{T}_{z2}_reco", pTBins, zpTmin, zpTmax); listOfTH1.push_back(hz2pTreco);
        TH1F *hz1mreco = new TH1F("z1_m_reco", "m_{z1}_reco", mBins, zmmin, zmmax); listOfTH1.push_back(hz1mreco);
        TH1F *hz2mreco = new TH1F("z2_m_reco", "m_{z2}_reco", mBins, zmmin, zmmax); listOfTH1.push_back(hz2mreco);
        // phi
        TH1F *hz1phireco = new TH1F("z1_#phi_reco", "#phi_{z1}_reco", phiBins, -TMath::Pi(), +TMath::Pi()); listOfTH1.push_back(hz1phireco);
        TH1F *hz2phireco = new TH1F("z2_#phi_reco", "#phi_{z2}_reco", phiBins, -TMath::Pi(), +TMath::Pi()); listOfTH1.push_back(hz2phireco);
        TH1F *hzzdeltaPhireco = new TH1F("zz_#Delta#phi_reco", "#Delta#phi_{zz}_reco", phiBins, -TMath::Pi(), +TMath::Pi()); listOfTH1.push_back(hzzdeltaPhireco);
        // eta
        TH1F *hz1etareco = new TH1F("z1_#eta_reco", "#eta_{z1}_reco", etaBins, zetamin, zetamax);listOfTH1.push_back(hz1etareco);
        TH1F *hz2etareco = new TH1F("z2_#eta_reco", "#eta_{z2}_reco", etaBins, zetamin, zetamax);listOfTH1.push_back(hz2etareco);
        TH1F *hzzdeltaEtareco = new TH1F("zz_#Delta#eta_reco", "#Delta#eta_{zz}_reco", etaBins, zetamin, zetamax);listOfTH1.push_back(hzzdeltaEtareco);
        // R
        TH1F *hz1Rreco = new TH1F("z1_R_reco", "R_{z1}_reco", RBins, zRmin, zRmax); listOfTH1.push_back(hz1Rreco);
        TH1F *hz2Rreco = new TH1F("z2_R_reco", "R_{z2}_reco", RBins, zRmin, zRmax); listOfTH1.push_back(hz2Rreco);
        TH1F *hzzdeltaRreco = new TH1F("zz_#DeltaR_reco", "#DeltaR_{zz}_reco", RBins, zRmin, zRmax); listOfTH1.push_back(hzzdeltaRreco);
        // cos 
        TH1F *hz1cosThetareco = new TH1F("z1_cos#theta_reco", "cos#theta_{z1}_reco", cosBins, -1, 1); listOfTH1.push_back(hz1cosThetareco);
        TH1F *hz2cosThetareco = new TH1F("z2_cos#theta_reco", "cos#theta_{z2}_reco", cosBins, -1, 1); listOfTH1.push_back(hz2cosThetareco);
    // particle
        // pT + m
        TH1F *hz1pTparticle = new TH1F("z1_pT_particle", "p^{T}_{z1}_particle", pTBins, zpTmin, zpTmax);listOfTH1.push_back(hz1pTparticle);
        TH1F *hz2pTparticle = new TH1F("z2_pT_particle", "p^{T}_{z2}_particle", pTBins, zpTmin, zpTmax);listOfTH1.push_back(hz2pTparticle);
        TH1F *hz1mparticle = new TH1F("z1_m_particle", "m_{z1}_particle", mBins, zmmin, zmmax);listOfTH1.push_back(hz1mparticle);
        TH1F *hz2mparticle = new TH1F("z2_m_particle", "m_{z2}_particle", mBins, zmmin, zmmax);listOfTH1.push_back(hz2mparticle);
        // phi
        TH1F *hz1phiparticle = new TH1F("z1_#phi_particle", "#phi_{z1}_particle", phiBins, -TMath::Pi(), +TMath::Pi());listOfTH1.push_back(hz1phiparticle);
        TH1F *hz2phiparticle = new TH1F("z2_#phi_particle", "#phi_{z2}_particle", phiBins, -TMath::Pi(), +TMath::Pi());listOfTH1.push_back(hz2phiparticle);
        TH1F *hzzdeltaPhiparticle = new TH1F("zz_#Delta#phi_particle", "#Delta#phi_{zz}_particle", phiBins, -TMath::Pi(), +TMath::Pi()); listOfTH1.push_back(hzzdeltaPhiparticle);
        // eta
        TH1F *hz1etaparticle = new TH1F("z1_#eta_particle", "#eta_{z1}_particle", etaBins, zetamin, zetamax); listOfTH1.push_back(hz1etaparticle);
        TH1F *hz2etaparticle = new TH1F("z2_#eta_particle", "#eta_{z2}_particle", etaBins, zetamin, zetamax); listOfTH1.push_back(hz2etaparticle);
        TH1F *hzzdeltaEtaparticle = new TH1F("zz_#Delta#eta_particle", "#Delta#eta_{zz}_particle", etaBins, zetamin, zetamax); listOfTH1.push_back(hzzdeltaEtaparticle);
        // R
        TH1F *hz1Rparticle = new TH1F("z1_R_particle", "R_{z1}_particle", RBins, zRmin, zRmax); listOfTH1.push_back(hz1Rparticle);
        TH1F *hz2Rparticle = new TH1F("z2_R_particle", "R_{z2}_particle", RBins, zRmin, zRmax); listOfTH1.push_back(hz2Rparticle);
        TH1F *hzzdeltaRparticle = new TH1F("zz_#DeltaR_particle", "#DeltaR_{zz}_particle", RBins, zRmin, zRmax); listOfTH1.push_back(hzzdeltaRparticle);
    // parton
        // pT + m
        TH1F *hz1pTparton = new TH1F("z1_pT_parton", "p^{T}_{z1}_parton", pTBins, zpTmin, zpTmax); listOfTH1.push_back(hz1pTparton);
        TH1F *hz2pTparton = new TH1F("z2_pT_parton", "p^{T}_{z2}_parton", pTBins, zpTmin, zpTmax); listOfTH1.push_back(hz2pTparton);
        TH1F *hz1mparton = new TH1F("z1_m_parton", "m_{z1}_parton", mBins, zmmin, zmmax); listOfTH1.push_back(hz1mparton);
        TH1F *hz2mparton = new TH1F("z2_m_parton", "m_{z2}_parton", mBins, zmmin, zmmax); listOfTH1.push_back(hz2mparton);
        // phi
        TH1F *hz1phiparton = new TH1F("z1_#phi_parton", "#phi_{z1}_parton", phiBins, -TMath::Pi(), +TMath::Pi()); listOfTH1.push_back(hz1phiparton);
        TH1F *hz2phiparton = new TH1F("z2_#phi_parton", "#phi_{z2}_parton", phiBins, -TMath::Pi(), +TMath::Pi()); listOfTH1.push_back(hz2phiparton);
        TH1F *hzzdeltaPhiparton = new TH1F("zz_#Delta#phi_parton", "#Delta#phi_{zz}_parton", phiBins, -TMath::Pi(), +TMath::Pi()); listOfTH1.push_back(hzzdeltaPhiparton);
        // eta
        TH1F *hz1etaparton = new TH1F("z1_#eta_parton", "#eta_{z1}_parton", etaBins, zetamin, zetamax); listOfTH1.push_back(hz1etaparton);
        TH1F *hz2etaparton = new TH1F("z2_#eta_parton", "#eta_{z2}_parton", etaBins, zetamin, zetamax); listOfTH1.push_back(hz2etaparton);
        TH1F *hzzdeltaEtaparton = new TH1F("zz_#Delta#eta_parton", "#Delta#eta_{zz}_parton", etaBins, zetamin, zetamax); listOfTH1.push_back(hzzdeltaEtaparton);
        // R
        TH1F *hz1Rparton = new TH1F("z1_R_parton", "R_{z1}_parton", RBins, zRmin, zRmax); listOfTH1.push_back(hz1Rparton);
        TH1F *hz2Rparton = new TH1F("z2_R_parton", "R_{z2}_parton", RBins, zRmin, zRmax); listOfTH1.push_back(hz2Rparton);
        TH1F *hzzdeltaRparton = new TH1F("zz_#DeltaR_parton", "#DeltaR_{zz}_parton", RBins, zRmin, zRmax); listOfTH1.push_back(hzzdeltaRparton);

// leptons
    // reco
        // pT
        TH1F *hl1pTreco = new TH1F("l1_pT_reco", "p^{T}_{l1}_reco", pTBins, lpTmin, lpTmax); listOfTH1.push_back(hl1pTreco);
        TH1F *hl2pTreco = new TH1F("l2_pT_reco", "p^{T}_{l2}_reco", pTBins, lpTmin, lpTmax); listOfTH1.push_back(hl2pTreco);
        TH1F *hl3pTreco = new TH1F("l3_pT_reco", "p^{T}_{l3}_reco", pTBins, lpTmin, lpTmax); listOfTH1.push_back(hl3pTreco);
        TH1F *hl4pTreco = new TH1F("l4_pT_reco", "p^{T}_{l4}_reco", pTBins, lpTmin, lpTmax);  listOfTH1.push_back(hl4pTreco);
        // phi
        TH1F *hl1phireco = new TH1F("l1_#phi_reco", "#phi_{l1}_reco", phiBins, -TMath::Pi(), +TMath::Pi()); listOfTH1.push_back(hl1phireco);
        TH1F *hl2phireco = new TH1F("l2_#phi_reco", "#phi_{l2}_reco", phiBins, -TMath::Pi(), +TMath::Pi()); listOfTH1.push_back(hl2phireco);
        TH1F *hl3phireco = new TH1F("l3_#phi_reco", "#phi_{l3}_reco", phiBins, -TMath::Pi(), +TMath::Pi()); listOfTH1.push_back(hl3phireco);
        TH1F *hl4phireco = new TH1F("l4_#phi_reco", "#phi_{l4}_reco", phiBins, -TMath::Pi(), +TMath::Pi()); listOfTH1.push_back(hl4phireco);
        TH1F *hl1l2deltaPhireco = new TH1F("l1l2_#Delta#phi_reco", "#Delta#phi_{l1l2}_reco", phiBins, -TMath::Pi(), +TMath::Pi()); listOfTH1.push_back(hl1l2deltaPhireco);
        TH1F *hl3l4deltaPhireco = new TH1F("l3l4_#Delta#phi_reco", "#Delta#phi_{l3l4}_reco", phiBins, -TMath::Pi(), +TMath::Pi()); listOfTH1.push_back(hl3l4deltaPhireco);
        TH1F *hl1l2deltaPhiBoostreco = new TH1F("l1l2_#Delta#phi_Boost_reco", "#Delta#phi_{l1l2}_reco", phiBins, -TMath::Pi(), +TMath::Pi()); listOfTH1.push_back(hl1l2deltaPhiBoostreco);
        TH1F *hl3l4deltaPhiBoostreco = new TH1F("l3l4_#Delta#phi_Boost_reco", "#Delta#phi_{l3l4}_reco", phiBins, -TMath::Pi(), +TMath::Pi()); listOfTH1.push_back(hl3l4deltaPhiBoostreco);
        // eta
        TH1F *hl1etareco = new TH1F("l1_#eta_reco", "#eta_{l1}_reco", etaBins, letamin, letamax);listOfTH1.push_back(hl1etareco);
        TH1F *hl2etareco = new TH1F("l2_#eta_reco", "#eta_{l2}_reco", etaBins, letamin, letamax); listOfTH1.push_back(hl2etareco);
        TH1F *hl3etareco = new TH1F("l3_#eta_reco", "#eta_{l3}_reco", etaBins, letamin, letamax); listOfTH1.push_back(hl3etareco);
        TH1F *hl4etareco = new TH1F("l4_#eta_reco", "#eta_{l4}_reco", etaBins, letamin, letamax); listOfTH1.push_back(hl4etareco);
        TH1F *hl1l2deltaEtareco = new TH1F("l1l2_#Delta#eta_reco", "#Delta#eta_{l1l2}_reco", etaBins, letamin, letamax); listOfTH1.push_back(hl1l2deltaEtareco);
        TH1F *hl3l4deltaEtareco = new TH1F("l3l4_#Delta#eta_reco", "#Delta#eta_{l3l4}_reco", etaBins, letamin, letamax); listOfTH1.push_back(hl3l4deltaEtareco);
        TH1F *hl1l2deltaEtaBoostreco = new TH1F("l1l2_#Delta#eta_Boost_reco", "#Delta#eta_{l1l2}_Boost_reco", etaBins, letamin, letamax); listOfTH1.push_back(hl1l2deltaEtaBoostreco);
        TH1F *hl3l4deltaEtaBoostreco = new TH1F("l3l4_#Delta#eta_Boost_reco", "#Delta#eta_{l3l4}_Boost_reco", etaBins, letamin, letamax); listOfTH1.push_back(hl3l4deltaEtaBoostreco);
        // R
        TH1F *hl1Rreco = new TH1F("l1_R_reco", "R_{l1}_reco", RBins, lRmin, lRmax); listOfTH1.push_back(hl1Rreco);
        TH1F *hl2Rreco = new TH1F("l2_R_reco", "R_{l2}_reco", RBins, lRmin, lRmax); listOfTH1.push_back(hl2Rreco);
        TH1F *hl3Rreco = new TH1F("l3_R_reco", "R_{l3}_reco", RBins, lRmin, lRmax); listOfTH1.push_back(hl3Rreco);
        TH1F *hl4Rreco = new TH1F("l4_R_reco", "R_{l4}_reco", RBins, lRmin, lRmax); listOfTH1.push_back(hl4Rreco);
        TH1F *hl1l2deltaRreco = new TH1F("l1l2_#Delta R_reco", "#Delta R_{l1l2}_reco", RBins, lRmin, lRmax); listOfTH1.push_back(hl1l2deltaRreco);
        TH1F *hl3l4deltaRreco = new TH1F("l3l4_#Delta R_reco", "#Delta R_{l3l4}_reco", RBins, lRmin, lRmax); listOfTH1.push_back(hl3l4deltaRreco);
        // cos
        TH1F *hl1cosThetareco = new TH1F("l1_cos#theta_reco", "cos#theta_{l1}_reco", cosBins, -1, 1); listOfTH1.push_back(hl1cosThetareco);
        TH1F *hl2cosThetareco = new TH1F("l2_cos#theta_reco", "cos#theta_{l2}_reco", cosBins, -1, 1); listOfTH1.push_back(hl2cosThetareco);
        TH1F *hl3cosThetareco = new TH1F("l3_cos#theta_reco", "cos#theta_{l3}_reco", cosBins, -1, 1); listOfTH1.push_back(hl3cosThetareco);
        TH1F *hl4cosThetareco = new TH1F("l4_cos#theta_reco", "cos#theta_{l4}_reco", cosBins, -1, 1); listOfTH1.push_back(hl4cosThetareco);
        TH1F *hfourlcosThetareco = new TH1F("fourl_cos#theta_reco", "cos#theta_{fourl}_reco", cosBins, -1, 1); listOfTH1.push_back(hfourlcosThetareco);
        TH1F *hl1cosThetaBoostreco = new TH1F("l1_cos#theta_Boost_reco", "cos#theta_{l1}_Boost_reco", cosBins, -1, 1); listOfTH1.push_back(hl1cosThetaBoostreco);
        TH1F *hl2cosThetaBoostreco = new TH1F("l2_cos#theta_Boost_reco", "cos#theta_{l2}_Boost_reco", cosBins, -1, 1); listOfTH1.push_back(hl2cosThetaBoostreco);
        TH1F *hl3cosThetaBoostreco = new TH1F("l3_cos#theta_Boost_reco", "cos#theta_{l3}_Boost_reco", cosBins, -1, 1); listOfTH1.push_back(hl3cosThetaBoostreco);
        TH1F *hl4cosThetaBoostreco = new TH1F("l4_cos#theta_Boost_reco", "cos#theta_{l4}_Boost_reco", cosBins, -1, 1); listOfTH1.push_back(hl4cosThetaBoostreco);
        TH1F *hfourlcosThetaBoostreco = new TH1F("fourl_cos#theta_Boost_reco", "cos#theta_{fourl}_Boost_reco", cosBins, -1, 1); listOfTH1.push_back(hfourlcosThetaBoostreco);
        TH1F *hl1l2CScosThetareco = new TH1F("l1l2_cos#theta_{CS}_reco", "cos#theta_{CSl1l2}_reco", cosBins, -1, 1); listOfTH1.push_back(hl1l2CScosThetareco);
        TH1F *hl3l4CScosThetareco = new TH1F("l3l4_cos#theta_{CS}_reco", "cos#theta_{CSl3l4}_reco", cosBins, -1, 1); listOfTH1.push_back(hl3l4CScosThetareco);

    // particle
        // pT
        TH1F *hl1pTparticle = new TH1F("l1_pT_particle", "p^{T}_{l1}_particle", pTBins, lpTmin, lpTmax); listOfTH1.push_back(hl1pTparticle);
        TH1F *hl2pTparticle = new TH1F("l2_pT_particle", "p^{T}_{l2}_particle", pTBins, lpTmin, lpTmax); listOfTH1.push_back(hl2pTparticle);
        TH1F *hl3pTparticle = new TH1F("l3_pT_particle", "p^{T}_{l3}_particle", pTBins, lpTmin, lpTmax); listOfTH1.push_back(hl3pTparticle);
        TH1F *hl4pTparticle = new TH1F("l4_pT_particle", "p^{T}_{l4}_particle", pTBins, lpTmin, lpTmax); listOfTH1.push_back(hl4pTparticle);
        // phi
        TH1F *hl1phiparticle = new TH1F("l1_#phi_particle", "#phi_{l1}_particle", phiBins, -TMath::Pi(), +TMath::Pi()); listOfTH1.push_back(hl1phiparticle);
        TH1F *hl2phiparticle = new TH1F("l2_#phi_particle", "#phi_{l2}_particle", phiBins, -TMath::Pi(), +TMath::Pi()); listOfTH1.push_back(hl2phiparticle);
        TH1F *hl3phiparticle = new TH1F("l3_#phi_particle", "#phi_{l3}_particle", phiBins, -TMath::Pi(), +TMath::Pi()); listOfTH1.push_back(hl3phiparticle);
        TH1F *hl4phiparticle = new TH1F("l4_#phi_particle", "#phi_{l4}_particle", phiBins, -TMath::Pi(), +TMath::Pi()); listOfTH1.push_back(hl4phiparticle);
        TH1F *hl1l2deltaPhiparticle = new TH1F("l1l2_#Delta#phi_particle", "#Delta#phi_{l1l2}_particle", phiBins, -TMath::Pi(), +TMath::Pi()); listOfTH1.push_back(hl1l2deltaPhiparticle);
        TH1F *hl3l4deltaPhiparticle = new TH1F("l3l4_#Delta#phi_particle", "#Delta#phi_{l3l4}_particle", phiBins, -TMath::Pi(), +TMath::Pi()); listOfTH1.push_back(hl3l4deltaPhiparticle);
        TH1F *hl1l2deltaPhiBoostparticle = new TH1F("l1l2_#Delta#phi_Boost_particle", "#Delta#phi_{l1l2}_particle", phiBins, -TMath::Pi(), +TMath::Pi()); listOfTH1.push_back(hl1l2deltaPhiBoostparticle);
        TH1F *hl3l4deltaPhiBoostparticle = new TH1F("l3l4_#Delta#phi_Boost_particle", "#Delta#phi_{l3l4}_particle", phiBins, -TMath::Pi(), +TMath::Pi()); listOfTH1.push_back(hl3l4deltaPhiBoostparticle);
        // eta
        TH1F *hl1etaparticle = new TH1F("l1_#eta_particle", "#eta_{l1}_particle", etaBins, letamin, letamax); listOfTH1.push_back(hl1etaparticle);
        TH1F *hl2etaparticle = new TH1F("l2_#eta_particle", "#eta_{l2}_particle", etaBins, letamin, letamax); listOfTH1.push_back(hl2etaparticle);
        TH1F *hl3etaparticle = new TH1F("l3_#eta_particle", "#eta_{l3}_particle", etaBins, letamin, letamax); listOfTH1.push_back(hl3etaparticle);
        TH1F *hl4etaparticle = new TH1F("l4_#eta_particle", "#eta_{l4}_particle", etaBins, letamin, letamax); listOfTH1.push_back(hl4etaparticle);
        TH1F *hl1l2deltaEtaparticle = new TH1F("l1l2_#Delta#eta_particle", "#Delta#eta_{l1l2}_particle", etaBins, letamin, letamax); listOfTH1.push_back(hl1l2deltaEtaparticle);
        TH1F *hl3l4deltaEtaparticle = new TH1F("l3l4_#Delta#eta_particle", "#Delta#eta_{l3l4}_particle", etaBins, letamin, letamax); listOfTH1.push_back(hl3l4deltaEtaparticle);
        TH1F *hl1l2deltaEtaBoostparticle = new TH1F("l1l2_#Delta#eta_Boost_particle", "#Delta#eta_{l1l2}_Boost_particle", etaBins, letamin, letamax); listOfTH1.push_back(hl1l2deltaEtaBoostparticle);
        TH1F *hl3l4deltaEtaBoostparticle = new TH1F("l3l4_#Delta#eta_Boost_particle", "#Delta#eta_{l3l4}_Boost_particle", etaBins, letamin, letamax);listOfTH1.push_back(hl3l4deltaEtaBoostparticle);
        // R
        TH1F *hl1Rparticle = new TH1F("l1_R_particle", "R_{l1}_particle", RBins, lRmin, lRmax); listOfTH1.push_back(hl1Rparticle);
        TH1F *hl2Rparticle = new TH1F("l2_R_particle", "R_{l2}_particle", RBins, lRmin, lRmax); listOfTH1.push_back(hl2Rparticle);
        TH1F *hl3Rparticle = new TH1F("l3_R_particle", "R_{l3}_particle", RBins, lRmin, lRmax); listOfTH1.push_back(hl3Rparticle);
        TH1F *hl4Rparticle = new TH1F("l4_R_particle", "R_{l4}_particle", RBins, lRmin, lRmax); listOfTH1.push_back(hl4Rparticle);
        TH1F *hl1l2deltaRparticle = new TH1F("l1l2_#Delta R_particle", "#Delta R_{l1l2}_particle", RBins, lRmin, lRmax); listOfTH1.push_back(hl1l2deltaRparticle);
        TH1F *hl3l4deltaRparticle = new TH1F("l3l4_#Delta R_particle", "#Delta R_{l3l4}_particle", RBins, lRmin, lRmax); listOfTH1.push_back(hl3l4deltaRparticle);
        // cos
        TH1F *hl1cosThetaparticle = new TH1F("l1_cos#theta_particle", "cos#theta_{l1}_particle", cosBins, -1, 1); listOfTH1.push_back(hl1cosThetaparticle);
        TH1F *hl2cosThetaparticle = new TH1F("l2_cos#theta_particle", "cos#theta_{l2}_particle", cosBins, -1, 1); listOfTH1.push_back(hl2cosThetaparticle);
        TH1F *hl3cosThetaparticle = new TH1F("l3_cos#theta_particle", "cos#theta_{l3}_particle", cosBins, -1, 1); listOfTH1.push_back(hl3cosThetaparticle);
        TH1F *hl4cosThetaparticle = new TH1F("l4_cos#theta_particle", "cos#theta_{l4}_particle", cosBins, -1, 1); listOfTH1.push_back(hl4cosThetaparticle);
        TH1F *hfourlcosThetaparticle = new TH1F("fourl_cos#theta_particle", "cos#theta_{fourl}_particle", cosBins, -1, 1); listOfTH1.push_back(hfourlcosThetaparticle);
        TH1F *hl1cosThetaBoostparticle = new TH1F("l1_cos#theta_Boost_particle", "cos#theta_{l1}_Boost_particle", cosBins, -1, 1); listOfTH1.push_back(hl1cosThetaBoostparticle);
        TH1F *hl2cosThetaBoostparticle = new TH1F("l2_cos#theta_Boost_particle", "cos#theta_{l2}_Boost_particle", cosBins, -1, 1); listOfTH1.push_back(hl2cosThetaBoostparticle);
        TH1F *hl3cosThetaBoostparticle = new TH1F("l3_cos#theta_Boost_particle", "cos#theta_{l3}_Boost_particle", cosBins, -1, 1); listOfTH1.push_back(hl3cosThetaBoostparticle);
        TH1F *hl4cosThetaBoostparticle = new TH1F("l4_cos#theta_Boost_particle", "cos#theta_{l4}_Boost_particle", cosBins, -1, 1); listOfTH1.push_back(hl4cosThetaBoostparticle);
        TH1F *hfourlcosThetaBoostparticle = new TH1F("fourl_cos#theta_Boost_particle", "cos#theta_{fourl}_Boost_particle", cosBins, -1, 1); listOfTH1.push_back(hfourlcosThetaBoostparticle);
        TH1F *hl1l2CScosThetaparticle = new TH1F("l1l2_cos#theta_{CS}_particle", "cos#theta_{CSl1l2}_particle", cosBins, -1, 1); listOfTH1.push_back(hl1l2CScosThetaparticle);
        TH1F *hl3l4CScosThetaparticle = new TH1F("l3l4_cos#theta_{CS}_particle", "cos#theta_{CSl3l4}_particle", cosBins, -1, 1); listOfTH1.push_back(hl3l4CScosThetaparticle);

    // parton
        // pT
        TH1F *hl1pTparton = new TH1F("l1_pT_parton", "p^{T}_{l1}_parton", pTBins, lpTmin, lpTmax); listOfTH1.push_back(hl1pTparton);
        TH1F *hl2pTparton = new TH1F("l2_pT_parton", "p^{T}_{l2}_parton", pTBins, lpTmin, lpTmax); listOfTH1.push_back(hl2pTparton);
        TH1F *hl3pTparton = new TH1F("l3_pT_parton", "p^{T}_{l3}_parton", pTBins, lpTmin, lpTmax);  listOfTH1.push_back(hl3pTparton);
        TH1F *hl4pTparton = new TH1F("l4_pT_parton", "p^{T}_{l4}_parton", pTBins, lpTmin, lpTmax); listOfTH1.push_back(hl4pTparton);
        // phi
        TH1F *hl1phiparton = new TH1F("l1_#phi_parton", "#phi_{l1}_parton", phiBins, -TMath::Pi(), +TMath::Pi()); listOfTH1.push_back(hl1phiparton);
        TH1F *hl2phiparton = new TH1F("l2_#phi_parton", "#phi_{l2}_parton", phiBins, -TMath::Pi(), +TMath::Pi()); listOfTH1.push_back(hl2phiparton);
        TH1F *hl3phiparton = new TH1F("l3_#phi_parton", "#phi_{l3}_parton", phiBins, -TMath::Pi(), +TMath::Pi()); listOfTH1.push_back(hl3phiparton);
        TH1F *hl4phiparton = new TH1F("l4_#phi_parton", "#phi_{l4}_parton", phiBins, -TMath::Pi(), +TMath::Pi()); listOfTH1.push_back(hl4phiparton);
        TH1F *hl1l2deltaPhiparton = new TH1F("l1l2_#Delta#phi_parton", "#Delta#phi_{l1l2}_parton", phiBins, -TMath::Pi(), +TMath::Pi()); listOfTH1.push_back(hl1l2deltaPhiparton);
        TH1F *hl3l4deltaPhiparton = new TH1F("l3l4_#Delta#phi_parton", "#Delta#phi_{l3l4}_parton", phiBins, -TMath::Pi(), +TMath::Pi()); listOfTH1.push_back(hl3l4deltaPhiparton);
        TH1F *hl1l2deltaPhiBoostparton = new TH1F("l1l2_#Delta#phi_Boost_parton", "#Delta#phi_{l1l2}_parton", phiBins, -TMath::Pi(), +TMath::Pi()); listOfTH1.push_back(hl1l2deltaPhiBoostparton);
        TH1F *hl3l4deltaPhiBoostparton = new TH1F("l3l4_#Delta#phi_Boost_parton", "#Delta#phi_{l3l4}_parton", phiBins, -TMath::Pi(), +TMath::Pi()); listOfTH1.push_back(hl3l4deltaPhiBoostparton);
        // eta
        TH1F *hl1etaparton = new TH1F("l1_#eta_parton", "#eta_{l1}_parton", etaBins, letamin, letamax); listOfTH1.push_back(hl1etaparton);
        TH1F *hl2etaparton = new TH1F("l2_#eta_parton", "#eta_{l2}_parton", etaBins, letamin, letamax); listOfTH1.push_back(hl2etaparton);
        TH1F *hl3etaparton = new TH1F("l3_#eta_parton", "#eta_{l3}_parton", etaBins, letamin, letamax); listOfTH1.push_back(hl3etaparton);
        TH1F *hl4etaparton = new TH1F("l4_#eta_parton", "#eta_{l4}_parton", etaBins, letamin, letamax); listOfTH1.push_back(hl4etaparton);
        TH1F *hl1l2deltaEtaparton = new TH1F("l1l2_#Delta#eta_parton", "#Delta#eta_{l1l2}_parton", etaBins, letamin, letamax); listOfTH1.push_back(hl1l2deltaEtaparton);
        TH1F *hl3l4deltaEtaparton = new TH1F("l3l4_#Delta#eta_parton", "#Delta#eta_{l3l4}_parton", etaBins, letamin, letamax); listOfTH1.push_back(hl3l4deltaEtaparton);
        TH1F *hl1l2deltaEtaBoostparton = new TH1F("l1l2_#Delta#eta_Boost_parton", "#Delta#eta_{l1l2}_Boost_parton", etaBins, letamin, letamax); listOfTH1.push_back(hl1l2deltaEtaBoostparton);
        TH1F *hl3l4deltaEtaBoostparton = new TH1F("l3l4_#Delta#eta_Boost_parton", "#Delta#eta_{l3l4}_Boost_parton", etaBins, letamin, letamax); listOfTH1.push_back(hl3l4deltaEtaBoostparton);
        // R
        TH1F *hl1Rparton = new TH1F("l1_R_parton", "R_{l1}_parton", RBins, lRmin, lRmax); listOfTH1.push_back(hl1Rparton);
        TH1F *hl2Rparton = new TH1F("l2_R_parton", "R_{l2}_parton", RBins, lRmin, lRmax); listOfTH1.push_back(hl2Rparton);
        TH1F *hl3Rparton = new TH1F("l3_R_parton", "R_{l3}_parton", RBins, lRmin, lRmax); listOfTH1.push_back(hl3Rparton);
        TH1F *hl4Rparton = new TH1F("l4_R_parton", "R_{l4}_parton", RBins, lRmin, lRmax); listOfTH1.push_back(hl4Rparton);
        TH1F *hl1l2deltaRparton = new TH1F("l1l2_#Delta R_parton", "#Delta R_{l1l2}_parton", RBins, lRmin, lRmax); listOfTH1.push_back(hl1l2deltaRparton);
        TH1F *hl3l4deltaRparton = new TH1F("l3l4_#Delta R_parton", "#Delta R_{l3l4}_parton", RBins, lRmin, lRmax); listOfTH1.push_back(hl3l4deltaRparton);
        // cos
        TH1F *hl1cosThetaparton = new TH1F("l1_cos#theta_parton", "cos#theta_{l1}_parton", cosBins, -1, 1); listOfTH1.push_back(hl1cosThetaparton);
        TH1F *hl2cosThetaparton = new TH1F("l2_cos#theta_parton", "cos#theta_{l2}_parton", cosBins, -1, 1); listOfTH1.push_back(hl2cosThetaparton);
        TH1F *hl3cosThetaparton = new TH1F("l3_cos#theta_parton", "cos#theta_{l3}_parton", cosBins, -1, 1); listOfTH1.push_back(hl3cosThetaparton);
        TH1F *hl4cosThetaparton = new TH1F("l4_cos#theta_parton", "cos#theta_{l4}_parton", cosBins, -1, 1); listOfTH1.push_back(hl4cosThetaparton);
        TH1F *hfourlcosThetaparton = new TH1F("fourl_cos#theta_parton", "cos#theta_{fourl}_parton", cosBins, -1, 1); listOfTH1.push_back(hfourlcosThetaparton);
        TH1F *hl1cosThetaBoostparton = new TH1F("l1_cos#theta_Boost_parton", "cos#theta_{l1}_Boost_parton", cosBins, -1, 1); listOfTH1.push_back(hl1cosThetaBoostparton);
        TH1F *hl2cosThetaBoostparton = new TH1F("l2_cos#theta_Boost_parton", "cos#theta_{l2}_Boost_parton", cosBins, -1, 1); listOfTH1.push_back(hl2cosThetaBoostparton);
        TH1F *hl3cosThetaBoostparton = new TH1F("l3_cos#theta_Boost_parton", "cos#theta_{l3}_Boost_parton", cosBins, -1, 1); listOfTH1.push_back(hl3cosThetaBoostparton);
        TH1F *hl4cosThetaBoostparton = new TH1F("l4_cos#theta_Boost_parton", "cos#theta_{l4}_Boost_parton", cosBins, -1, 1); listOfTH1.push_back(hl4cosThetaBoostparton);
        TH1F *hfourlcosThetaBoostparton = new TH1F("fourl_cos#theta_Boost_parton", "cos#theta_{fourl}_Boost_parton", cosBins, -1, 1); listOfTH1.push_back(hfourlcosThetaBoostparton);
        TH1F *hl1l2CScosThetaparton = new TH1F("l1l2_cos#theta_{CS}_parton", "cos#theta_{CSl1l2}_parton", cosBins, -1, 1); listOfTH1.push_back(hl1l2CScosThetaparton);
        TH1F *hl3l4CScosThetaparton = new TH1F("l3l4_cos#theta_{CS}_parton", "cos#theta_{CSl3l4}_parton", cosBins, -1, 1); listOfTH1.push_back(hl3l4CScosThetaparton); 

// comp
// 2D - parton(1) particle(2) reco(3)

    // higgs
      // pT + m
        TH2F *hHpT12Comp = new TH2F("H_pT_comp_12", "p_{T}^{hbb}", pTBins, hpTmin, hpTmax, pTBins, hpTmin, hpTmax); listOfTH2.push_back(hHpT12Comp);
        TH2F *hHm12Comp = new TH2F("H_m_comp_12", "m_{hbb}", mBins, hmmin, hmmax, mBins, hmmin, hmmax); listOfTH2.push_back(hHm12Comp);
        TH2F *hHpT23Comp = new TH2F("H_pT_comp_23", "p_{T}^{hbb}", pTBins, hpTmin, hpTmax, pTBins, hpTmin, hpTmax);  listOfTH2.push_back(hHpT23Comp);
        TH2F *hHm23Comp = new TH2F("H_m_comp_23", "m_{hbb}", mBins, hmmin, hmmax, mBins, hmmin, hmmax); listOfTH2.push_back(hHm23Comp);
        TH2F *hHpT13Comp = new TH2F("H_pT_comp_13", "p_{T}^{hbb}", pTBins, hpTmin, hpTmax, pTBins, hpTmin, hpTmax); listOfTH2.push_back(hHpT13Comp);
        TH2F *hHm13Comp = new TH2F("H_m_comp_13", "m_{hbb}", mBins, hmmin, hmmax, mBins, hmmin, hmmax); listOfTH2.push_back(hHm13Comp);
        TH2F*hbbdeltaPhi12Comp = new TH2F("bb_#Delta#phi_comp_12", "#Delta#phi_{bb}", phiBins, -TMath::Pi(),+TMath::Pi(), phiBins, -TMath::Pi(),+TMath::Pi()); listOfTH2.push_back(hbbdeltaPhi12Comp);
        TH2F*hbbdeltaPhi23Comp = new TH2F("bb_#Delta#phi_comp_23", "#Delta#phi_{bb}", phiBins, -TMath::Pi(),+TMath::Pi(), phiBins, -TMath::Pi(),+TMath::Pi()); listOfTH2.push_back(hbbdeltaPhi23Comp);
        TH2F*hbbdeltaPhi13Comp = new TH2F("bb_#Delta#phi_comp_13", "#Delta#phi_{bb}", phiBins, -TMath::Pi(),+TMath::Pi(), phiBins, -TMath::Pi(),+TMath::Pi()); listOfTH2.push_back(hbbdeltaPhi13Comp);
        TH2F*hbbdeltaEta12Comp = new TH2F("bb_#Delta#eta_comp_12", "#Delta#eta_{bb}", etaBins, hetamin, hetamax, etaBins, hetamin, hetamax); listOfTH2.push_back(hbbdeltaEta12Comp);
        TH2F*hbbdeltaEta23Comp = new TH2F("bb_#Delta#eta_comp_23", "#Delta#eta_{bb}", etaBins, hetamin, hetamax, etaBins, hetamin, hetamax); listOfTH2.push_back(hbbdeltaEta23Comp);
        TH2F*hbbdeltaEta13Comp = new TH2F("bb_#Delta#eta_comp_13", "#Delta#eta_{bb}", etaBins, hetamin, hetamax, etaBins, hetamin, hetamax); listOfTH2.push_back(hbbdeltaEta13Comp);
    // jets
        TH2F*hjjpT12Comp = new TH2F("jj_pT_comp_12", "p_{T}^{jj}", pTBins, jpTmin, jpTmax, pTBins, jpTmin, jpTmax); listOfTH2.push_back(hjjpT12Comp);
        TH2F*hjjpT23Comp = new TH2F("jj_pT_comp_23", "p_{T}^{jj}", pTBins, jpTmin, jpTmax, pTBins, jpTmin, jpTmax); listOfTH2.push_back(hjjpT23Comp);
        TH2F*hjjpT13Comp = new TH2F("jj_pT_comp_13", "p_{T}^{jj}", pTBins, jpTmin, jpTmax, pTBins, jpTmin, jpTmax); listOfTH2.push_back(hjjpT13Comp);
        TH2F*hj1Phi23Comp = new TH2F("j1_#phi_comp_23", "#phi_{j1}", phiBins, -TMath::Pi(),+TMath::Pi(), phiBins, -TMath::Pi(),+TMath::Pi()); listOfTH2.push_back(hj1Phi23Comp);
        TH2F*hj2Phi23Comp = new TH2F("j2_#phi_comp_23", "#phi_{j2}", phiBins, -TMath::Pi(),+TMath::Pi(), phiBins, -TMath::Pi(),+TMath::Pi()); listOfTH2.push_back(hj2Phi23Comp);
        TH2F*hjjdeltaPhi12Comp = new TH2F("jj_#Delta#phi_comp_12", "#Delta#phi_{jj}", phiBins, -TMath::Pi(),+TMath::Pi(), phiBins, -TMath::Pi(),+TMath::Pi()); listOfTH2.push_back(hjjdeltaPhi12Comp);
        TH2F*hjjdeltaPhi23Comp = new TH2F("jj_#Delta#phi_comp_23", "#Delta#phi_{jj}", phiBins, -TMath::Pi(),+TMath::Pi(), phiBins, -TMath::Pi(),+TMath::Pi()); listOfTH2.push_back(hjjdeltaPhi23Comp);
        TH2F*hjjdeltaPhi13Comp = new TH2F("jj_#Delta#phi_comp_13", "#Delta#phi_{jj}", phiBins, -TMath::Pi(),+TMath::Pi(), phiBins, -TMath::Pi(),+TMath::Pi()); listOfTH2.push_back(hjjdeltaPhi13Comp);
    //z
      // pT + m
        TH2F*hz1pT12Comp = new TH2F("z1_pT_comp_12", "p^{T}_{z1}", pTBins, zpTmin, zpTmax, pTBins, zpTmin, zpTmax); listOfTH2.push_back(hz1pT12Comp);
        TH2F*hz1m12Comp = new TH2F("z1_m_comp_12", "m_{z1}", mBins, zmmin, zmmax, mBins, zmmin, zmmax); listOfTH2.push_back(hz1m12Comp);
        TH2F*hz2pT12Comp = new TH2F("z2_pT_comp_12", "p^{T}_{z2}", pTBins, zpTmin, zpTmax, pTBins, zpTmin, zpTmax); listOfTH2.push_back(hz2pT12Comp);
        TH2F*hz2m12Comp = new TH2F("z2_m_comp_12", "m_{z2}", mBins, zmmin, zmmax, mBins, zmmin, zmmax); listOfTH2.push_back(hz2m12Comp);
        TH2F*hz1pT23Comp = new TH2F("z1_pT_comp_23", "p^{T}_{z1}", pTBins, zpTmin, zpTmax, pTBins, zpTmin, zpTmax); listOfTH2.push_back(hz1pT23Comp);
        TH2F*hz1m23Comp = new TH2F("z1_m_comp_23", "m_{z1}", mBins, zmmin, zmmax, mBins, zmmin, zmmax); listOfTH2.push_back(hz1m23Comp);
        TH2F*hz2pT23Comp = new TH2F("z2_pT_comp_23", "p^{T}_{z2}", pTBins, zpTmin, zpTmax, pTBins, zpTmin, zpTmax); listOfTH2.push_back(hz2pT23Comp);
        TH2F*hz2m23Comp = new TH2F("z2_m_comp_23", "m_{z2}", mBins, zmmin, zmmax, mBins, zmmin, zmmax); listOfTH2.push_back(hz2m23Comp);
        TH2F*hz1pT13Comp = new TH2F("z1_pT_comp_13", "p^{T}_{z1}", pTBins, zpTmin, zpTmax, pTBins, zpTmin, zpTmax); listOfTH2.push_back(hz1pT13Comp);
        TH2F*hz1m13Comp = new TH2F("z1_m_comp_13", "m_{z1}", mBins, zmmin, zmmax, mBins, zmmin, zmmax); listOfTH2.push_back(hz1m13Comp);
        TH2F*hz2pT13Comp = new TH2F("z2_pT_comp_13", "p^{T}_{z2}", pTBins, zpTmin, zpTmax, pTBins, zpTmin, zpTmax); listOfTH2.push_back(hz2pT13Comp);
        TH2F*hz2m13Comp = new TH2F("z2_m_comp_13", "m_{z2}", mBins, zmmin, zmmax, mBins, zmmin, zmmax); listOfTH2.push_back(hz2m13Comp);

        TH2F*hl1l2deltaPhiHpTcompreco = new TH2F("l1l2_delta#phi_H_pT_comp_reco", "l1l2_delta#phi_H_pT_comp_reco", phiBins, -TMath::Pi(),+TMath::Pi(), pTBins, hpTmin, hpTmax); listOfTH2.push_back(hl1l2deltaPhiHpTcompreco);
         TH2F*hl3l4deltaPhiHpTcompreco = new TH2F("l3l4_delta#phi_H_pT_comp_reco", "l1l2_delta#phi_H_pT_comp_reco", phiBins, -TMath::Pi(),+TMath::Pi(), pTBins, hpTmin, hpTmax); listOfTH2.push_back(hl3l4deltaPhiHpTcompreco);

        TH2F*hHz1pTcompreco = new TH2F("h_z1_pT_comp_reco", "h_z1_pT_comp_reco", pTBins, hpTmin, hpTmax, pTBins, zpTmin, zpTmax); listOfTH2.push_back(hHz1pTcompreco);
        TH2F*hHz2pTcompreco = new TH2F("h_z2_pT_comp_reco", "h_z2_pT_comp_reco", pTBins, hpTmin, hpTmax, pTBins, zpTmin, zpTmax); listOfTH2.push_back(hHz2pTcompreco);
        TH2F*hHz1pTcompparticle = new TH2F("h_z1_pT_comp_particle", "h_z1_pT_comp_particle", pTBins, hpTmin, hpTmax, pTBins, zpTmin, zpTmax); listOfTH2.push_back(hHz1pTcompparticle);
        TH2F*hHz2pTcompparticle = new TH2F("h_z2_pT_comp_particle", "h_z2_pT_comp_particle", pTBins, hpTmin, hpTmax, pTBins, zpTmin, zpTmax); listOfTH2.push_back(hHz2pTcompparticle);
        TH2F*hHz1pTcompparton = new TH2F("h_z1_pT_comp_parton", "h_z1_pT_comp_parton", pTBins, hpTmin, hpTmax, pTBins, zpTmin, zpTmax); listOfTH2.push_back(hHz1pTcompparton);
        TH2F*hHz2pTcompparton = new TH2F("h_z2_pT_comp_parton", "h_z2_pT_comp_parton", pTBins, hpTmin, hpTmax, pTBins, zpTmin, zpTmax); listOfTH2.push_back(hHz2pTcompparton);

        TH2F*hHpTz1etacompreco = new TH2F("h_pT_z1_eta_comp_reco", "h_pT_z1_eta_comp_reco", pTBins, hpTmin, hpTmax, etaBins, zetamin, zetamax); listOfTH2.push_back(hHpTz1etacompreco);
        TH2F*hHpTz2etacompreco = new TH2F("h_pT_z2_eta_comp_reco", "h_pT_z2_eta_comp_reco", pTBins, hpTmin, hpTmax, etaBins, zetamin, zetamax); listOfTH2.push_back(hHpTz2etacompreco);
        TH2F*hHpTz1etacompparticle = new TH2F("h_pT_z1_eta_comp_particle", "h_pT_z1_eta_comp_particle", pTBins, hpTmin, hpTmax, etaBins, zetamin, zetamax); listOfTH2.push_back(hHpTz1etacompparticle);
        TH2F*hHpTz2etacompparticle = new TH2F("h_pT_z2_eta_comp_particle", "h_pT_z2_eta_comp_particle", pTBins, hpTmin, hpTmax, etaBins, zetamin, zetamax); listOfTH2.push_back(hHpTz2etacompparticle);
        TH2F*hHpTz1etacompparton = new TH2F("h_pT_z1_eta_comp_parton", "h_pT_z1_eta_comp_parton", pTBins, hpTmin, hpTmax, etaBins, zetamin, zetamax); listOfTH2.push_back(hHpTz1etacompparton);
        TH2F*hHpTz2etacompparton= new TH2F("h_pT_z2_eta_comp_parton", "h_pT_z2_eta_comp_parton", pTBins, hpTmin, hpTmax, etaBins, zetamin, zetamax); listOfTH2.push_back(hHpTz2etacompparton);

        TH2F*hHzzpTcompreco = new TH2F("h_zz_pT_comp_reco", "h_zz_pT_comp_reco", pTBins, hpTmin, 2*hpTmax, pTBins, zpTmin, 2*zpTmax); listOfTH2.push_back(hHzzpTcompreco);
        TH2F*hHzzpTcompparticle = new TH2F("h_zz_pT_comp_particle", "h_zz_pT_comp_particle", pTBins, hpTmin, 2*hpTmax, pTBins, zpTmin, 2*zpTmax); listOfTH2.push_back(hHzzpTcompparticle);
        TH2F*hHzzpTcompparton = new TH2F("h_zz_pT_comp_parton", "h_zz_pT_comp_parton", pTBins, hpTmin, 2*hpTmax, pTBins, zpTmin, 2*zpTmax); listOfTH2.push_back(hHzzpTcompparton);

	TProfile *kappaLambda = new TProfile("kappaLambda", "kappaLambda", 40, -20, 20);
	kappaLambda -> GetXaxis() -> SetTitle("#kappa_{#lambda}");
	
  double  nPassed=0;
  double Lumi=1;//3e3;
  double totWeightedEntries=0;
  int  nPassedRaw=0;
  int nQuads=0;

  
  TLorentzVector j1_reco, j1_particle,  j1_parton;
  TLorentzVector j2_reco, j2_particle,  j2_parton;

  TLorentzVector b1_reco, b1_particle,  b1_parton;
  TLorentzVector b2_reco, b2_particle,  b2_parton;

  TLorentzVector h_reco, h_parton, h_particle;

  TLorentzVector e1_reco, e1_particle,  e1_parton;
  TLorentzVector e2_reco, e2_particle,  e2_parton;
  TLorentzVector e3_reco, e3_particle,  e3_parton;
  TLorentzVector e4_reco, e4_particle,  e4_parton;

  TLorentzVector m1_reco, m1_particle,  m1_parton;
  TLorentzVector m2_reco, m2_particle,  m2_parton;
  TLorentzVector m3_reco, m3_particle,  m3_parton;
  TLorentzVector m4_reco, m4_particle,  m4_parton;

  TLorentzVector l1_reco, l1_particle,  l1_parton;
  TLorentzVector l2_reco, l2_particle,  l2_parton;
  TLorentzVector l3_reco, l3_particle,  l3_parton;
  TLorentzVector l4_reco, l4_particle,  l4_parton;

  TLorentzVector z1_reco, z1_particle,  z1_parton;
  TLorentzVector z2_reco, z2_particle,  z2_parton;

  TLorentzVector fourl_reco, fourl_particle, fourl_parton; 

  int q1_reco=0;
  int q2_reco=0;
  int q3_reco=0;
  int q4_reco=0;

  int q1_particle=0;
  int q2_particle=0;
  int q3_particle=0;
  int q4_particle=0;

// GB: Here it will do strange thing in case it does not pass the cut...
  double l1l2deltaPhireco= -99999;
  double l3l4deltaPhireco=-99999;
  
  double l1l2deltaEtareco= -99999;
  double l3l4deltaEtareco=-99999; 
  
  double l1l2deltaRreco=-99999;
  double l3l4deltaRreco=-99999;
  
  double l1cosThetareco=-99999;
  double l2cosThetareco=-99999;
  double l3cosThetareco=-99999;
  double l4cosThetareco=-99999;
  double fourlcosThetareco=-99999;
   double l1cosThetaBoostreco=-99999;
  double l2cosThetaBoostreco=-99999;
  double l3cosThetaBoostreco=-99999;
  double l4cosThetaBoostreco=-99999;
  double fourlcosThetaBoostreco=-99999;
  double l1l2deltaPhiBoostreco=-99999;
  double l3l4deltaPhiBoostreco=-99999;
  double l1l2deltaEtaBoostreco=-99999;
  double l3l4deltaEtaBoostreco=-99999;
  // collins soper frame
  double l1l2CScosThetareco=-99999;
  double l3l4CScosThetareco=-99999;
  
  double zzdeltaPhireco= -99999;
  double zzdeltaEtareco= -99999;
  double zzdeltaRreco=-99999;

  double bbdeltaPhireco = 9999;
  double bbdeltaEtareco = 9999;
  double bbdeltaRreco = -9999;



  double l1cosThetaBoostparton=-9999;
  double l2cosThetaBoostparton=-9999;
  double l3cosThetaBoostparton=-9999;
  double l4cosThetaBoostparton=-9999;
  double fourlcosThetaBoostparton=-9999;

  double l1l2deltaPhiparton=-9999;
    double l3l4deltaPhiparton=-9999;
    double l1l2deltaEtaparton=-9999;
    double l3l4deltaEtaparton=-9999;
    
    double l1l2deltaRparton=-9999;
    double l3l4deltaRparton=-9999;
    
    double l1cosThetaparton=-9999;
    double l2cosThetaparton=-9999;
    double l3cosThetaparton=-9999;
    double l4cosThetaparton=-9999;
    double fourlcosThetaparton=-9999;
    
    int q1_parton = -9999;
    int q2_parton = -9999;
    int q3_parton = -9999;
    int q4_parton = -9999;

    double bbdeltaRparton=-9999;
    double jjdeltaPhiparton = -9999;
    double jjdeltaEtaparton = -9999;
    double bbdeltaPhiparton = -9999;
    double bbdeltaEtaparton = -9999;

    double l1l2CScosThetaparton=-9999;
    double l3l4CScosThetaparton=-9999;
    
    double zzdeltaPhiparton=-9999;
    double zzdeltaEtaparton=-9999;
    double zzdeltaRparton=-9999;
  
  
/*
  DelphesLHEFReader *reader = new DelphesLHEFReader;
  std::map<int,double> mapKappaLambda;
    mapKappaLambda[1]=+10.0;
    mapKappaLambda[2]=-10.0;
    mapKappaLambda[3]=-5.0;
    mapKappaLambda[4]=-2.0;
    mapKappaLambda[5]=-1.0;
    mapKappaLambda[6]=-.1;
    mapKappaLambda[7]=+.1;
    mapKappaLambda[8]=+1; //problematic?
    mapKappaLambda[9]=+2.0;
    mapKappaLambda[10]=+5.0;
    mapKappaLambda[11]=+100;
*/

  double sumOfWeights=0;
  TH1F *hClosure = new TH1F("hClosure","hClosure",1,0,1);
  listOfTH1.push_back(hClosure);
  

//------------------------------------------------------------------------------------------------------------------------------------------------------------
// loop
//------------------------------------------------------------------------------------------------------------------------------------------------------------
  for(Int_t entry = 0; entry < numberOfEntries; ++entry) {
    //if( entry > 100) break;
    treeReader->ReadEntry(entry);
    std::map<int,double> kappaLambdaWeights;
    HepMCEvent *event = (HepMCEvent*) branchEvent -> At(0);
    Float_t weight = event->Weight;///numberOfEntries*Lumi;
    sumOfWeights+=weight;
    weight = weight * Lumi / numberOfEntries; 
    hClosure->Fill(0.5,weight);
    totWeightedEntries+=weight;
    

/*
 // cout << "Reading Event: " << entry << ", Weight:" << event->Weight << endl;
    Float_t weight = event->Weight/numberOfEntries*Lumi;
    totWeightedEntries+=weight;
    double nomXsec=6.9801096e-09;//1.5428498234628936e-06;
    //LHEFEvent event = (LHEFEvent) branchEvent -> At(0);
    reader->AnalyzeWeight((ExRootTreeBranch*) branchWeight);
   // cout<<"Internal event weights"<<" nominal weight "<<weight<<endl;
    //for(int i=145; i<(int)branchWeight->GetEntries(); i++){
    if(kappaVal > 0 && kappaVal < 12){
      int startVal=154;
      for(int i=startVal; i<startVal+10; i++){
        Weight *WeveWeight=(Weight*)branchWeight->At(i);
        LHEFWeight *eveWeight=(LHEFWeight*)branchWeight->At(i);
     // cout<<" Event "<<entry<<" weight entry i "<<i<<" "<<eveWeight->Weight<<" id is " <<eveWeight->GetUniqueID()<<" ID "<<eveWeight->ID<<endl;
     // cout<<" Event "<<entry<<" weight entry i "<<i<<" "<<WeveWeight->Weight<<" id is " <<WeveWeight->GetUniqueID()<<endl;
        kappaLambdaWeights[i-startVal+1]=WeveWeight->Weight;
     //  cout<<"Kappa lambda value "<<mapKappaLambda[i-startVal+1]<<endl;
        kappaLambda->Fill(mapKappaLambda[i-startVal+1],WeveWeight->Weight/nomXsec); 
      }
      // cout<<endl;
    double kappa10=1;
    kappa10=kappaLambdaWeights[kappaVal];
    weight*=kappa10/kappaLambdaWeights[8];
    }
*/

//------------------------------------------------------------------------------------------------------------------------------------------------------------
// reco
//------------------------------------------------------------------------------------------------------------------------------------------------------------

// higgs

int switchVal_reco = 0;

// cout<<"Start reco analysis "<<endl;

    vector <int> btagIndex;
    vector <int> noBtag;
    vector <int> goodJetIndex;

    for(int i=0; i<(int)branchJet->GetEntries(); i++){
      Jet *jet=(Jet*) branchJet->At(i);
    //if( jet->PT < 20) continue;
    //if (fabs(jet->Eta) > 4.4) continue; 
    //  if( jet->BTag>0) btagIndex.push_back(i);
      if( isMyBTag(jet, branchGenParticle,0,0.4,btagEff,fakeEff) && abs(jet->Eta) < 2.5 ) btagIndex.push_back(i); 

      else noBtag.push_back(i);
           goodJetIndex.push_back(i);
    }

    increaseCount(cutFlowMap_reco,"initial reco",weight);
    //cutFlowMap_reco["initial reco"] = {cutVal_reco,cutValW_reco};

   sort(btagIndex.begin(), btagIndex.end(), [branchJet](const int& lhs, const int& rhs) {
       return ((Jet*)branchJet->At(lhs))->PT > ((Jet*)branchJet->At(rhs))->PT;
     });
   sort(noBtag.begin(), noBtag.end(), [branchJet](const int& lhs, const int& rhs) {
       return ((Jet*)branchJet->At(lhs))->PT > ((Jet*)branchJet->At(rhs))->PT;
     });
   sort(goodJetIndex.begin(), goodJetIndex.end(), [branchJet](const int& lhs, const int& rhs) {
       return ((Jet*)branchJet->At(lhs))->PT > ((Jet*)branchJet->At(rhs))->PT;
     });

   //cout<<"Cut 1 btag reco";
   if(switchVal_reco == 0 && btagIndex.size() > 1) { // at least one b tag 
     increaseCount(cutFlowMap_reco,"1 btag reco",weight);
     //cutFlowMap_reco["1 btag reco"] = {cutVal_reco,cutValW_reco};
     //cout<<" passed "<<endl;
   } else{
     switchVal_reco = 1;
   //cout<<" failed "<<endl;
   }
   
  
   if(switchVal_reco == 0 && goodJetIndex.size() > 2) { // at least two jets 
     increaseCount(cutFlowMap_reco,"2 good j reco",weight);
     //cutFlowMap_reco["2 good j reco"] = {cutVal_reco,cutValW_reco};
   } else  switchVal_reco = 1;
   
   
   Jet *b1=nullptr;
   Jet *b2=nullptr;
   
   vector<pair<int,int>> bJetPairs;
   vector<vector <int>> bJetPairsComb;

   if(switchVal_reco == 0 && goodJetIndex.size()>1 )
     bJetPairsComb= combinationsNoRepetitionAndOrderDoesNotMatter(2,goodJetIndex);
   //  vector<vector <int>> bJetPairsComb=combinationsNoRepetitionAndOrderDoesNotMatter(2,btagIndex);
   
   // This is a cut !!
   if( switchVal_reco == 0  && bJetPairsComb.size() >= 1) { ////continue; // need at least two good jets;x
     increaseCount(cutFlowMap_reco,"2 b-like jet pairs reco",weight);
     //cutFlowMap_reco["2 b-like jet pairs reco"] = {cutVal_reco,cutValW_reco};
   }
   else switchVal_reco = 1;
   
   for(int i=0; i<(int)bJetPairsComb.size(); i++)
     bJetPairs.push_back(make_pair(bJetPairsComb[i][0],bJetPairsComb[i][1]));
   
   if( bJetPairs.size() > 1) 
     sort(bJetPairs.begin(), bJetPairs.end(), [branchJet](const pair<int,int> lhs, const pair<int,int> rhs) {
	 return fabs(((((Jet*)branchJet->At(lhs.first))->P4() + ((Jet*)branchJet->At(lhs.second))->P4())).M() - 125 ) <
	   fabs( ((((Jet*)branchJet->At(rhs.first))->P4() + ((Jet*)branchJet->At(rhs.second))->P4()).M()) - 125 ) ; 
       });
   
   pair <int,int> higgsbbcandidate;
   bool foundBjet=false; 
   for(int i=0; i<(int) bJetPairs.size(); i++){
     b1=(Jet*)branchJet->At(bJetPairs[i].first);
     b2=(Jet*)branchJet->At(bJetPairs[i].second);
     //if( b1->BTag>0 && b2->BTag>0) {
     // Attention 
     if( isMyBTag(b1, branchGenParticle,0,0.4,btagEff,fakeEff) && abs(b1->Eta) < 2.5 || (isMyBTag(b2, branchGenParticle,0,0.4,btagEff,fakeEff) && abs(b2->Eta)<2.5) ) {
       higgsbbcandidate=bJetPairs[i];
       foundBjet=true;
       break;
     }
   }

   if(switchVal_reco == 0 && foundBjet) { // b pair
     increaseCount(cutFlowMap_reco,"found bb reco",weight);
   } else  switchVal_reco = 1;
   
   if(switchVal_reco == 0 && b1 !=nullptr && b2 !=nullptr && foundBjet){
     b1_reco = b1->P4();
     b2_reco = b2->P4();
     h_reco = b1_reco + b2_reco;

     if (b1_reco.Eta() > b2_reco.Eta()) {
       bbdeltaPhireco = remainder( b1_reco.Phi() - b2_reco.Phi(), 2*M_PI );
       bbdeltaEtareco = b1_reco.Eta() - b2_reco.Eta();
     }
     else{
       bbdeltaPhireco = remainder( b2_reco.Phi() - b1_reco.Phi(), 2*M_PI );
       bbdeltaEtareco = b2_reco.Eta() - b1_reco.Eta();
     }
     
     // double bbdeltaPhireco=(b1_reco.Phi() > b2_reco.Phi() ? -1:+1)*TMath::Abs(b2_reco.Phi() - b1_reco.Phi());
     // double bbdeltaEtareco=(b1_reco.Eta() > b2_reco.Eta() ? -1:+1)*TMath::Abs(b2_reco.Eta() - b1_reco.Eta());
     // double bbdeltaPhireco= (b1_reco.Phi() - b2_reco.Phi());
     // double bbdeltaEtareco= (b1_reco.Eta() - b2_reco.Eta());
     bbdeltaRreco=sqrt((bbdeltaPhireco*bbdeltaPhireco)+(bbdeltaEtareco*bbdeltaEtareco));
   }
   
   
   
    

// VBF jets
    vector <int> nonHiggsJet;

    for(int i=0; i<(int)goodJetIndex.size(); i++){
      if( goodJetIndex[i] == higgsbbcandidate.first  || goodJetIndex[i] == higgsbbcandidate.second) continue;
      nonHiggsJet.push_back(i);
    }

    if(nonHiggsJet.size() > 1) {
      sort(nonHiggsJet.begin(), nonHiggsJet.end(), [branchJet](const int& lhs, const int& rhs) {
	  return ((Jet*)branchJet->At(lhs))->PT > ((Jet*)branchJet->At(rhs))->PT;
	});}
    
    if(switchVal_reco == 0 && nonHiggsJet.size() > 1) { // at least 2 vbf jets
    increaseCount(cutFlowMap_reco,"2 vbfj reco",weight);
    //cutFlowMap_reco["2 vbfj reco"] = {cutVal_reco,cutValW_reco};
  
  } else  switchVal_reco = 1;
  
  vector<pair<int,int>> vbfJetIndex;
  vector<vector <int>> vbfJetIndexComb;
  if(switchVal_reco==0 && nonHiggsJet.size() > 1 ) {
    increaseCount(cutFlowMap_reco,"vbfj pairs",weight);
    vbfJetIndexComb=combinationsNoRepetitionAndOrderDoesNotMatter(2,nonHiggsJet);
  }
  else switchVal_reco=1;

  //if(switchVal_reco==0)    
  for(int i=0; i<(int)vbfJetIndexComb.size(); i++)
    vbfJetIndex.push_back(make_pair(vbfJetIndexComb[i][0],vbfJetIndexComb[i][1]));
  
  // GB: What is this ? 
  //  if( branchMuon->GetEntries() + branchElectron->GetEntries() < 4) continue;
  
  //if( switchVal_reco==0 && vbfJetIndex.size() > 1) 
  //sort(vbfJetIndex.begin(), vbfJetIndex.end(), [branchJet](const pair<int,int> lhs, const pair<int,int> rhs) {
  //	return fabs((((Jet*)branchJet->At(lhs.first))->Eta - ((Jet*)branchJet->At(lhs.second))->Eta) ) >
  //	  fabs((((Jet*)branchJet->At(rhs.first))->Eta - ((Jet*)branchJet->At(rhs.second))->Eta) ) ; 
  //  });

  // sort by pT; 
  //if( switchVal_reco==0 && vbfJetIndex.size() > 1) 
  //sort(vbfJetIndex.begin(), vbfJetIndex.end(), [branchJet](const pair<int,int> lhs, const pair<int,int> rhs) {
  //	return (((Jet*)branchJet->At(lhs.first))->PT) < (((Jet*)branchJet->At(lhs.second))->PT);
  //  });
  
  
    int vbfJetsIndexCandidate = -1;
    vector<pair<int,int>> vbfJetIndex2;
    // loop and take first w eta > 2.5
    for (int i=0; i<(int)vbfJetIndex.size(); i++) {
      if( fabs((((Jet*)branchJet->At(vbfJetIndex[i].first))->Eta - ((Jet*)branchJet->At(vbfJetIndex[i].second))->Eta)) <= 2.5 ) {
	continue; 
      }
      else {
	vbfJetIndex2.push_back(vbfJetIndex.at(i));
	vbfJetsIndexCandidate = i;
      }
    }

    vbfJetIndex=vbfJetIndex2;
    if(switchVal_reco==0 && vbfJetIndex.size()>0)
      increaseCount(cutFlowMap_reco,"2.5 deltaEta vbf reco",weight);
    else
      switchVal_reco=1; 
      
    
    //Jet *jet1 = (Jet*) branchJet->At(vbfJetIndex[vbfJetsIndexCandidate].first);
    //Jet *jet2 = (Jet*) branchJet->At(vbfJetIndex[vbfJetsIndexCandidate].second);
    
    // GB, well what happens here if we don't have two jets in the event?
    //Jet *jet1 = (Jet*) branchJet->At(vbfJetIndex[0].first);
    //Jet *jet2 = (Jet*) branchJet->At(vbfJetIndex[0].second);
    
      
    // re sort by eta 
    if( switchVal_reco==0 && vbfJetIndex.size() > 1) 
      sort(vbfJetIndex.begin(), vbfJetIndex.end(), [branchJet](const pair<int,int> lhs, const pair<int,int> rhs) {
	  return fabs((((Jet*)branchJet->At(lhs.first))->Eta - ((Jet*)branchJet->At(lhs.second))->Eta) ) >
	    fabs((((Jet*)branchJet->At(rhs.first))->Eta - ((Jet*)branchJet->At(rhs.second))->Eta) ) ; 
	});
    
    Jet *jet1 =nullptr;
    Jet *jet2 =nullptr;
  double jjdeltaPhireco =  -9999;
  double jjdeltaEtareco= -9999;
  double jjdeltaRreco = -9999; 
    
  
  if(switchVal_reco==0){
    jet1 = (Jet*) branchJet->At(vbfJetIndex[0].first);
    jet2 = (Jet*) branchJet->At(vbfJetIndex[0].second);
    j1_reco=jet1->P4();
    j2_reco=jet2->P4();
    
    if (j1_reco.Eta() > j2_reco.Eta()) {
      jjdeltaPhireco = remainder( j1_reco.Phi() - j2_reco.Phi(), 2*M_PI );
    }
    else{
      jjdeltaPhireco = remainder( j2_reco.Phi() - j1_reco.Phi(), 2*M_PI );
    }
    // double jjdeltaPhireco=(j1_reco.Phi() > j2_reco.Phi() ? -1:+1)*TMath::Abs(j2_reco.Phi() - j1_reco.Phi());
    // double jjdeltaEtareco=(j1_reco.Eta() > j2_reco.Eta() ? -1:+1)*TMath::Abs(j2_reco.Eta() - j1_reco.Eta());
    // double jjdeltaPhireco= (j1_reco.Phi() - j2_reco.Phi());
    jjdeltaEtareco= (j1_reco.Eta() - j2_reco.Eta());
    jjdeltaRreco=sqrt((jjdeltaPhireco*jjdeltaPhireco)+(jjdeltaEtareco*jjdeltaEtareco));
  }

  
  //---------------------------------------------------------------------------------------------
  //---------------------------------------------------------------------------------------------
  
// leptons + z

    vector <int> goodE_reco_indices; 
    for(int i=0; i<(int)branchElectron->GetEntries(); i++){
      Electron *el_reco = (Electron *) branchElectron->At(i);

      // pT and eta cuts 
      if( el_reco->PT > 5.0 && fabs(el_reco->Eta) < 2.5) goodE_reco_indices.push_back(i);
      
    }

    // sort the indices by pT ;
    sort(goodE_reco_indices.begin(), goodE_reco_indices.end(), [branchElectron](const int& lhs, const int& rhs) {
	return ((Electron*)branchElectron->At(lhs))->PT > ((Electron*)branchElectron->At(rhs))->PT;
      });

    vector <int> goodMu_reco_indices; 
    for(int i=0; i<(int)branchMuon->GetEntries(); i++){
      Muon *mu_reco = (Muon *) branchMuon->At(i);
    
    // pT and eta cuts 
      if( mu_reco->PT > 5.0 && fabs(mu_reco->Eta) < 2.5) goodMu_reco_indices.push_back(i);

    }

    // sort the indices by pT ;
    sort(goodMu_reco_indices.begin(), goodMu_reco_indices.end(), [branchMuon](const int& lhs, const int& rhs) {
	return ((Muon*)branchMuon->At(lhs))->PT > ((Muon*)branchMuon->At(rhs))->PT;
      });

    //cout<<endl;
    //cout<<"Electrons "<<endl;
    //for(int i=0;i<(int)goodE_reco_indices.size(); i++)
    //cout<<" e Pt: "<<((Electron*)branchElectron->At(goodE_reco_indices[i]))->PT<<" eta "<<((Electron*)branchElectron->At(goodE_reco_indices[i]))->Eta<<" i "<<i<<endl;
    //cout<<"Muons "<<endl;
    //for(int i=0;i<(int)goodMu_reco_indices.size(); i++)
    //cout<<" mu Pt: "<<((Muon*)branchMuon->At(goodMu_reco_indices[i]))->PT<<" eta "<<((Muon*)branchMuon->At(goodMu_reco_indices[i]))->Eta<<" i "<<i<<endl;
    
    // form pairs for each flavour

    // electrons 
    vector< pair<int,int>> elecRecoPairIndices;
    vector< pair<int,int>> elecRecoPairIndicesIn;
    vector <vector<int>> elecRecoPairIndices_;
    if( goodE_reco_indices.size() > 1 )
    elecRecoPairIndices_=combinationsNoRepetitionAndOrderDoesNotMatter(2,goodE_reco_indices);

    // remove all combinations not satisfying criteria 
    for(int i=0;i<(int)elecRecoPairIndices_.size(); i++){
      int elecRecoIndex=elecRecoPairIndices_[i].at(0);
      int elecRecoIndex2=elecRecoPairIndices_[i].at(1);
      
      Electron *el1_reco=(Electron*) branchElectron->At(elecRecoIndex);
      Electron *el2_reco=(Electron*) branchElectron->At(elecRecoIndex2);

      if( el1_reco->Charge == el2_reco->Charge ) continue;
      
      elecRecoPairIndicesIn.push_back(make_pair(elecRecoIndex,elecRecoIndex2));
    }

    elecRecoPairIndices=elecRecoPairIndicesIn;
    
    sort(elecRecoPairIndices.begin(), elecRecoPairIndices.end(), [branchElectron](const pair<int,int> lhs, const pair<int,int> rhs) {
    	return fabs(((((Electron*)branchElectron->At(lhs.first))->P4() + ((Electron*)branchElectron->At(lhs.second))->P4())).M() -91 ) <
	  fabs( ((((Electron*)branchElectron->At(rhs.first))->P4() + ((Electron*)branchElectron->At(rhs.second))->P4()).M()) -91 ) ; 
      });

    remove_overlaps(elecRecoPairIndices);

    // muons 
    vector< pair<int,int>> muRecoPairIndicesIn;
    vector< pair<int,int>> muRecoPairIndices;
    vector <vector<int>> muRecoPairIndices_;
    if(goodMu_reco_indices.size() > 1) 
      muRecoPairIndices_=combinationsNoRepetitionAndOrderDoesNotMatter(2,goodMu_reco_indices);

    // remove all combinations not satifing criteria 
    for(int i=0; i<(int)muRecoPairIndices_.size(); i++){
      int elecRecoIndex=muRecoPairIndices_[i].at(0);
      int elecRecoIndex2=muRecoPairIndices_[i].at(1);

      Muon *el1_reco= (Muon*) branchMuon->At(elecRecoIndex);
      Muon *el2_reco=(Muon*) branchMuon->At(elecRecoIndex2);
      if( el1_reco->Charge ==  el2_reco->Charge ) continue;
      muRecoPairIndicesIn.push_back(make_pair(elecRecoIndex,elecRecoIndex2));
    }

    muRecoPairIndices=muRecoPairIndicesIn;

    sort(muRecoPairIndices.begin(), muRecoPairIndices.end(), [branchMuon](const pair<int,int> lhs, const pair<int,int> rhs) {
    	return fabs(((((Muon*)branchMuon->At(lhs.first))->P4() + ((Muon*)branchMuon->At(lhs.second))->P4())).M() -91 ) <
	  fabs( ((((Muon*)branchMuon->At(rhs.first))->P4() + ((Muon*)branchMuon->At(rhs.second))->P4()).M()) -91 ) ; 
      });
    
    //sort(muRecoPairIndices.begin(),muRecoPairIndices.end(), [branchMuon]( pair<int,int>   & lhs,  pair<int,int>   & rhs) {
    //	    int index11_reco=(lhs).first;
    //	    int index12_reco=(lhs).second;
    //	    int index21_reco=(rhs).first;
    //	int index22_reco=(rhs).second;
    //	    return fabs(((((Muon*)branchMuon->At(index11_reco))->P4() + ((Muon*)branchMuon->At(index12_reco))->P4())).M() - 91) <
    //	  fabs( ((((Muon*)branchMuon->At(index21_reco))->P4() + ((Muon*)branchMuon->At(index22_reco))->P4()).M()) -91);
    //});
    
    remove_overlaps(muRecoPairIndices);


    //increaseCount(cutFlowMap_reco,"at least two lep pairs",weight);
    
    // GB: Isn't this a cut ??
    // order quads by mZ1 and mZ2
    //if( muRecoPairIndices.size() + elecRecoPairIndices.size() < 2) continue; // no candidate found;
    //if( switchVal_reco==0 && ( muRecoPairIndices.size() + elecRecoPairIndices.size() > 1)) {
    //increaseCount(cutFlowMap_reco,"at least two lep pairs",weight);
    ////cutFlowMap_reco["at least two lep pairs"] = {cutVal_reco,cutValW_reco};
    //}
    //else switchVal_reco=1; 
    
    int thisRecoEventType=-1; 
    
    vector<pair<int,pair<int,int>>> RecoPairIndices; // 0 for electron 1 for muon
    for(int i=0; i<(int) elecRecoPairIndices.size(); i++)
      RecoPairIndices.push_back(make_pair(0,elecRecoPairIndices.at(i)));
    for(int i=0; i<(int) muRecoPairIndices.size(); i++)
      RecoPairIndices.push_back(make_pair(1,muRecoPairIndices.at(i)));

    // sort all of the indices by closeness to mZ
    sort(RecoPairIndices.begin(), RecoPairIndices.end(), [branchMuon,branchElectron]  ( const pair<int,pair<int,int>>  lhs , const pair<int,pair<int,int>>   rhs ){

	double mass1_reco=0 ;
    	double mass2_reco=0;
	
	if( lhs.first==0 ) mass1_reco = (((Electron*)branchElectron->At(lhs.second.first))->P4() + ((Electron*)branchElectron->At(lhs.second.second))->P4()).M();
	else mass1_reco = (((Muon*)branchMuon->At(lhs.second.first))->P4() + ((Muon*)branchMuon->At(lhs.second.second))->P4()).M();
	
	if( rhs.first==0 ) mass2_reco = (((Electron*)branchElectron->At(rhs.second.first))->P4() + ((Electron*)branchElectron->At(rhs.second.second))->P4()).M();
	else mass2_reco = (((Muon*)branchMuon->At(rhs.second.first))->P4() + ((Muon*)branchMuon->At(rhs.second.second))->P4()).M();
	
	return fabs(mass1_reco -91.0 )  <  fabs(mass2_reco -91.0); 
	
      });

    //cout<<"Reco pair indices "<<RecoPairIndices.size()<<endl;
    //for( int i=0; i<(int)RecoPairIndices.size(); i++){
    //if(RecoPairIndices[i].first==0) 
    //	cout<<" "<<"electrons "<< (((Electron*)branchElectron->At(RecoPairIndices[i].second.first))->P4() + ((Electron*)branchElectron->At(RecoPairIndices[i].second.second))->P4()).M()<<" "<<i<<endl;
    //else
    //	cout<<" "<<"muons "<< (((Muon*)branchMuon->At(RecoPairIndices[i].second.first))->P4() + ((Muon*)branchMuon->At(RecoPairIndices[i].second.second))->P4()).M()<<"  "<<i<<endl;
    //}
    
    if(RecoPairIndices.size() < 2 ) switchVal_reco=1;
    else if (switchVal_reco==0) increaseCount(cutFlowMap_reco,"OSFL",weight);
    ///cutFlowMap_reco["OSFL"] = {cutVal_reco,cutValW_reco};
    
    
    if( switchVal_reco==0){
      if( RecoPairIndices[0].first == 1 && RecoPairIndices[1].first == 1) thisRecoEventType=0;
      else if( RecoPairIndices[0].first == 0 && RecoPairIndices[1].first == 0) thisRecoEventType=1;
      else if( RecoPairIndices[0].first == 1 && RecoPairIndices[1].first == 0) thisRecoEventType=2;
      else if( RecoPairIndices[0].first == 0 && RecoPairIndices[1].first == 1) thisRecoEventType=3;
    }

    // case 4mu 
    if( thisRecoEventType==0 ) {
      
      // take first two muons
      Muon *muon1_reco = (Muon *) branchMuon->At( RecoPairIndices[0].second.first);
      Muon *muon2_reco = (Muon *) branchMuon->At( RecoPairIndices[0].second.second);
      // take first two muons
      Muon *muon3_reco = (Muon *) branchMuon->At( RecoPairIndices[1].second.first);
      Muon *muon4_reco = (Muon *) branchMuon->At( RecoPairIndices[1].second.second);
      
      q1_reco = muon1_reco->Charge;
      q2_reco = muon2_reco->Charge;
      q3_reco = muon3_reco->Charge;
      q4_reco = muon4_reco->Charge;
      
      l1_reco=muon1_reco->P4();
      l2_reco=muon2_reco->P4();
      l3_reco=muon3_reco->P4();
      l4_reco=muon4_reco->P4(); 
      
    }
    
  // case 4e
  else if( thisRecoEventType==1 ) {
    
    // take first two electrons
    Electron *muon1_reco = (Electron *) branchElectron->At( RecoPairIndices[0].second.first);
    Electron *muon2_reco = (Electron *) branchElectron->At( RecoPairIndices[0].second.second);
    // take first two electrons
    Electron *muon3_reco = (Electron *) branchElectron->At( RecoPairIndices[1].second.first);
    Electron *muon4_reco = (Electron *) branchElectron->At( RecoPairIndices[1].second.second);

    q1_reco = muon1_reco->Charge;
    q2_reco = muon2_reco->Charge;
    q3_reco = muon3_reco->Charge;
    q4_reco = muon4_reco->Charge;

    l1_reco=muon1_reco->P4();
    l2_reco=muon2_reco->P4();
    l3_reco=muon3_reco->P4();
    l4_reco=muon4_reco->P4(); 
  }

  // case 2mu2e 
  if( thisRecoEventType==2 ) {
    
    // take first two muons
    Muon *muon1_reco = (Muon *) branchMuon->At( RecoPairIndices[0].second.first);
    Muon *muon2_reco = (Muon *) branchMuon->At( RecoPairIndices[0].second.second);
    // take first two electrons
    Electron *muon3_reco = (Electron *) branchElectron->At( RecoPairIndices[1].second.first);
    Electron *muon4_reco = (Electron *) branchElectron->At( RecoPairIndices[1].second.second);

    q1_reco = muon1_reco->Charge;
    q2_reco = muon2_reco->Charge;
    q3_reco = muon3_reco->Charge;
    q4_reco = muon4_reco->Charge;

    l1_reco=muon1_reco->P4();
    l2_reco=muon2_reco->P4();
    l3_reco=muon3_reco->P4();
    l4_reco=muon4_reco->P4(); 
  }

  // case 2e2mu 
   else if( thisRecoEventType==3 ) {
    
    // take first two electrons
    Electron *muon1_reco = (Electron *) branchElectron->At( RecoPairIndices[0].second.first);
    Electron *muon2_reco = (Electron *) branchElectron->At( RecoPairIndices[0].second.second);
    // take first two muons
    Muon *muon3_reco = (Muon *) branchMuon->At( RecoPairIndices[1].second.first);
    Muon *muon4_reco = (Muon *) branchMuon->At( RecoPairIndices[1].second.second);

    q1_reco = muon1_reco->Charge;
    q2_reco = muon2_reco->Charge;
    q3_reco = muon3_reco->Charge;
    q4_reco = muon4_reco->Charge;

    l1_reco=muon1_reco->P4();
    l2_reco=muon2_reco->P4();
    l3_reco=muon3_reco->P4();
    l4_reco=muon4_reco->P4(); 

   }

  if (thisRecoEventType==-1)  switchVal_reco=1;
 
  
  
  if(  switchVal_reco==0){
    z1_reco=l1_reco + l2_reco;
    z2_reco=l3_reco + l4_reco;
    fourl_reco=l1_reco + l2_reco + l3_reco + l4_reco;
    
    l1l2deltaPhireco=(l1_reco.Phi() > l2_reco.Phi() ? -1:+1)*TMath::Abs(l2_reco.Phi() - l1_reco.Phi());
    l3l4deltaPhireco=(l3_reco.Phi() > l4_reco.Phi() ? -1:+1)*TMath::Abs(l4_reco.Phi() - l3_reco.Phi());
    
    l1l2deltaEtareco=(l1_reco.Eta() > l2_reco.Eta() ? -1:+1)*TMath::Abs(l2_reco.Eta() - l1_reco.Eta());
    l3l4deltaEtareco=(l3_reco.Eta() > l4_reco.Eta() ? -1:+1)*TMath::Abs(l4_reco.Eta() - l3_reco.Eta());
  
    l1l2deltaRreco=sqrt((l1l2deltaPhireco*l1l2deltaPhireco)+(l1l2deltaEtareco*l1l2deltaEtareco));
    l3l4deltaRreco=sqrt((l3l4deltaPhireco*l3l4deltaPhireco)+(l3l4deltaEtareco*l3l4deltaEtareco));
  
    l1cosThetareco=l1_reco.CosTheta();
    l2cosThetareco=l2_reco.CosTheta();
    l3cosThetareco=l3_reco.CosTheta();
    l4cosThetareco=l4_reco.CosTheta();
    fourlcosThetareco=fourl_reco.CosTheta();
  
    l1_reco.Boost(-z1_reco.BoostVector());
    l2_reco.Boost(-z1_reco.BoostVector());
    l3_reco.Boost(-z2_reco.BoostVector());
    l4_reco.Boost(-z2_reco.BoostVector());
    l1cosThetaBoostreco=l1_reco.CosTheta();
    l2cosThetaBoostreco=l2_reco.CosTheta();
    l3cosThetaBoostreco=l3_reco.CosTheta();
    l4cosThetaBoostreco=l4_reco.CosTheta();
    fourlcosThetaBoostreco=fourl_reco.CosTheta();
    l1l2deltaPhiBoostreco=(l1_reco.Phi() > l2_reco.Phi() ? -1:+1)*TMath::Abs(l2_reco.Phi() - l1_reco.Phi());
    l3l4deltaPhiBoostreco=(l3_reco.Phi() > l4_reco.Phi() ? -1:+1)*TMath::Abs(l4_reco.Phi() - l3_reco.Phi());
    l1l2deltaEtaBoostreco=(l1_reco.Eta() > l2_reco.Eta() ? -1:+1)*TMath::Abs(l2_reco.Eta() - l1_reco.Eta());
    l3l4deltaEtaBoostreco=(l3_reco.Eta() > l4_reco.Eta() ? -1:+1)*TMath::Abs(l4_reco.Eta() - l3_reco.Eta());
    l1_reco.Boost(z1_reco.BoostVector());
    l2_reco.Boost(z1_reco.BoostVector());
    l3_reco.Boost(z2_reco.BoostVector());
    l4_reco.Boost(z2_reco.BoostVector());
  
    // collins soper frame
    l1l2CScosThetareco=(q1_reco > q2_reco ? -1:+1)*TMath::Abs(2*(l2_reco.Pz()*l1_reco.E()-l1_reco.Pz()*l2_reco.E())/(z1_reco.M()*sqrt(z1_reco.M()*z1_reco.M()+z1_reco.Pt()*z1_reco.Pt())));
    l3l4CScosThetareco=(q3_reco > q4_reco ? -1:+1)*TMath::Abs(2*(l4_reco.Pz()*l3_reco.E()-l3_reco.Pz()*l4_reco.E())/(z2_reco.M()*sqrt(z2_reco.M()*z2_reco.M()+z2_reco.Pt()*z2_reco.Pt())));
  
    zzdeltaPhireco=(z1_reco.Phi() > z2_reco.Phi() ? -1:+1)*TMath::Abs(z2_reco.Phi() - z1_reco.Phi());
    zzdeltaEtareco=(z1_reco.Eta() > z2_reco.Eta() ? -1:+1)*TMath::Abs(z2_reco.Eta() - z1_reco.Eta());
    zzdeltaRreco=sqrt((zzdeltaPhireco*zzdeltaPhireco)+(zzdeltaEtareco*zzdeltaEtareco));
  }
  
  //cout<<"Done reco analysis  Event passed up to cut "<<switchVal_reco<<endl;
  
  
  //------------------------------------------------------------------------------------------------------------------------------------------------------------
  // particle
  //------------------------------------------------------------------------------------------------------------------------------------------------------------
  
  // higgs + jets

/*
    int bjets = 0;
    for(int i=0; i<(int)branchGenJet->GetEntries(); i++){
      Jet *genjet=(Jet*) branchGenJet->At(i);
      if (ghost_btag(branchGenParticle, genjet)){
        bjets += 1; 
      if (bjets == 1) {
        b1_particle = genjet->P4();
      } else if (bjets == 2){
        b2_particle = genjet->P4();
      } else break;
    }
  }


    h_particle = b1_particle + b2_particle;

    double bbdeltaPhiparticle=(b1_particle.Phi() > b2_particle.Phi() ? -1:+1)*TMath::Abs(b2_particle.Phi() - b1_particle.Phi());
    double bbdeltaEtaparticle=(b1_particle.Eta() > b2_particle.Eta() ? -1:+1)*TMath::Abs(b2_particle.Eta() - b1_particle.Eta());
    double bbdeltaRparticle=sqrt((bbdeltaPhiparticle*bbdeltaPhiparticle)+(bbdeltaEtaparticle*bbdeltaEtaparticle));



    int jetsParticle = 0;
    for(int i=0; i<(int)branchGenJet->GetEntries(); i++){
      Jet *genjet=(Jet*) branchGenJet->At(i);
      if (!ghost_btag(branchGenParticle, genjet)){
       jetsParticle += 1;
       if (jetsParticle == 1) {
       j1_particle = genjet->P4();
       } else if (jetsParticle == 2){
       j2_particle = genjet->P4();
       } else break;
       }
       }

       double jjdeltaPhiparticle=(j1_particle.Phi() > j2_particle.Phi() ? -1:+1)*TMath::Abs(j2_particle.Phi() - j1_particle.Phi());
       double jjdeltaEtaparticle=(j1_particle.Eta() > j2_particle.Eta() ? -1:+1)*TMath::Abs(j2_particle.Eta() - j1_particle.Eta());
       double jjdeltaRparticle=sqrt((jjdeltaPhiparticle*jjdeltaPhiparticle)+(jjdeltaEtaparticle*jjdeltaEtaparticle));
*/

//---------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------

// higgs
  double bbdeltaPhiparticle = -9999;
  double bbdeltaEtaparticle = -9999;
  double bbdeltaRparticle= -9999;
  Jet *jet1_particle = nullptr;
  Jet *jet2_particle = nullptr;
  double jjdeltaPhiparticle =  -9999;
  double jjdeltaEtaparticle =   -9999;
  double jjdeltaRparticle   =  -9999;
  
  double l1l2deltaPhiparticle=-9999;
  double l3l4deltaPhiparticle=-9999;
  
  double l1l2deltaEtaparticle=-9999;
  double l3l4deltaEtaparticle=-9999;

  double l1l2deltaRparticle=-9999;
  double l3l4deltaRparticle=-9999;

  double l1cosThetaparticle=-9999;
  double l2cosThetaparticle=-9999;
  double l3cosThetaparticle=-9999;
  double l4cosThetaparticle=-9999;
  double fourlcosThetaparticle=-9999;
  double l1cosThetaBoostparticle=-9999;
  double l2cosThetaBoostparticle=-9999;
  double l3cosThetaBoostparticle=-9999;
  double l4cosThetaBoostparticle=-9999;
  double fourlcosThetaBoostparticle=-9999;
  double l1l2CScosThetaparticle=-9999;
  double l3l4CScosThetaparticle=-9999;

  double zzdeltaPhiparticle=-9999;
  double zzdeltaEtaparticle=-9999;
  double zzdeltaRparticle=-9999;
  
  int switchVal_particle = 0;
  

    vector <int> btagIndexParticle;
    vector <int> noBtagParticle;
    vector <int> goodJetIndexParticle;

    for(int i=0; i<(int)branchGenJet->GetEntries(); i++){
      Jet *jet=(Jet*) branchGenJet->At(i);
    //if( jet->PT < 20) continue;
    //if (fabs(jet->Eta) > 4.4) continue; 
      // this is not correct as Btag>0 is not particle level.. 
      //if( jet->BTag>0) btagIndexParticle.push_back(i);
      if( ghost_btag(branchGenParticle, jet,1)) btagIndexParticle.push_back(i);
      else noBtagParticle.push_back(i);
           goodJetIndexParticle.push_back(i);
    }

    increaseCount(cutFlowMap_particle,"initial particle",weight);
    //cutFlowMap_particle["initial particle"] = {cutVal_particle,cutValW_particle};

    // at least one b tag 
    if(switchVal_particle == 0 && btagIndexParticle.size() > 0) increaseCount(cutFlowMap_particle,"1 btag particle",weight);
    else  switchVal_particle = 1;

    // at least two jets
    if(switchVal_particle == 0 && goodJetIndexParticle.size() > 1) increaseCount(cutFlowMap_particle,"2 good j particle",weight);
    else switchVal_particle = 1;
     
    
    
    sort(btagIndexParticle.begin(), btagIndexParticle.end(), [branchGenJet](const int& lhs, const int& rhs) {
	return ((Jet*)branchGenJet->At(lhs))->PT > ((Jet*)branchGenJet->At(rhs))->PT;
      });
    sort(noBtagParticle.begin(), noBtagParticle.end(), [branchGenJet](const int& lhs, const int& rhs) {
	return ((Jet*)branchGenJet->At(lhs))->PT > ((Jet*)branchGenJet->At(rhs))->PT;
      });
    sort(goodJetIndexParticle.begin(), goodJetIndexParticle.end(), [branchGenJet](const int& lhs, const int& rhs) {
	return ((Jet*)branchGenJet->At(lhs))->PT > ((Jet*)branchGenJet->At(rhs))->PT;
      });
    
    Jet *b1Particle=nullptr;
    Jet *b2Particle=nullptr;
    
    vector<pair<int,int>> bJetPairsParticle;
    vector<vector <int>> bJetPairsCombParticle;
    if(goodJetIndexParticle.size() > 1)
      bJetPairsCombParticle=combinationsNoRepetitionAndOrderDoesNotMatter(2,goodJetIndexParticle);
//  vector<vector <int>> bJetPairsCombParticle=combinationsNoRepetitionAndOrderDoesNotMatter(2,btagIndexParticle);
    
    //if( bJetPairsCombParticle.size() < 1) continue; // need at least two good jets;  
    if(switchVal_particle==0 && bJetPairsCombParticle.size() >0 ) 
      increaseCount(cutFlowMap_particle,"2 b-like jet pairs part",weight);
    else switchVal_particle=1; 

    for(int i=0; i<(int)bJetPairsCombParticle.size(); i++)
      bJetPairsParticle.push_back(make_pair(bJetPairsCombParticle[i][0],bJetPairsCombParticle[i][1]));
    
    if( bJetPairsParticle.size() > 1) 
      sort(bJetPairsParticle.begin(), bJetPairsParticle.end(), [branchGenJet](const pair<int,int> lhs, const pair<int,int> rhs) {
	  return fabs(((((Jet*)branchGenJet->At(lhs.first))->P4() + ((Jet*)branchGenJet->At(lhs.second))->P4())).M() - 125 ) <
	    fabs( ((((Jet*)branchGenJet->At(rhs.first))->P4() + ((Jet*)branchGenJet->At(rhs.second))->P4()).M()) - 125 ) ; 
	});
    
    pair <int,int> higgsbbcandidateParticle;
    bool foundBjetParticle=false; 
    for(int i=0; i<(int) bJetPairsParticle.size(); i++){
      b1Particle=(Jet*)branchGenJet->At(bJetPairsParticle[i].first);
      b2Particle=(Jet*)branchGenJet->At(bJetPairsParticle[i].second);
      if( ghost_btag(branchGenParticle, b1Particle) || ghost_btag(branchGenParticle, b2Particle)) {
	higgsbbcandidateParticle=bJetPairsParticle[i];
	foundBjetParticle=true;
	break;
      }
    }
    
    
    if(switchVal_particle == 0 && foundBjetParticle) { // b pair
      increaseCount(cutFlowMap_particle,"found bb particle",weight);
   } else  switchVal_particle = 1;
    
   if(switchVal_particle == 0 && b1Particle !=nullptr && b2Particle !=nullptr && foundBjetParticle){
    
    b1_particle = b1Particle->P4();
    b2_particle = b2Particle->P4();
    h_particle = b1_particle + b2_particle;
    
    
    if (b1_particle.Eta() > b2_particle.Eta()) {
      bbdeltaPhiparticle = remainder( b1_particle.Phi() - b2_particle.Phi(), 2*M_PI );
      bbdeltaEtaparticle= b1_particle.Eta() - b2_particle.Eta();
    }
    else{
      bbdeltaPhiparticle = remainder( b2_particle.Phi() - b1_particle.Phi(), 2*M_PI );
      bbdeltaEtaparticle= b2_particle.Eta() - b1_particle.Eta();
    }

    // double bbdeltaPhiparticle=(b1_particle.Phi() > b2_particle.Phi() ? -1:+1)*TMath::Abs(b2_particle.Phi() - b1_particle.Phi());
    // double bbdeltaEtaparticle=(b1_particle.Eta() > b2_particle.Eta() ? -1:+1)*TMath::Abs(b2_particle.Eta() - b1_particle.Eta());
    // double bbdeltaPhiparticle= (b1_particle.Phi() - b2_particle.Phi());
    // double bbdeltaEtaparticle= (b1_particle.Eta() - b2_particle.Eta());
    bbdeltaRparticle= sqrt((bbdeltaPhiparticle*bbdeltaPhiparticle)+(bbdeltaEtaparticle*bbdeltaEtaparticle));
   }
   
    // jets
   vector <int> nonHiggsJetParticle;
    
    for(int i=0; i<(int)goodJetIndexParticle.size(); i++){
      if( goodJetIndexParticle[i] == higgsbbcandidateParticle.first  || goodJetIndexParticle[i] == higgsbbcandidateParticle.second) continue;
      nonHiggsJetParticle.push_back(i);
    }

    if(nonHiggsJetParticle.size() > 1) {
      sort(nonHiggsJetParticle.begin(), nonHiggsJetParticle.end(), [branchGenJet](const int& lhs, const int& rhs) {
	return ((Jet*)branchGenJet->At(lhs))->PT > ((Jet*)branchGenJet->At(rhs))->PT;
	});}
    
    // at least 2 vbf jets 
    if(switchVal_particle == 0 && nonHiggsJetParticle.size() > 1) 
      increaseCount(cutFlowMap_particle,"2 vbfj particle",weight);
    else switchVal_particle = 1;
    

    vector<pair<int,int>> vbfJetIndexParticle;
    vector<vector <int>> vbfJetIndexCombParticle;
    if(switchVal_particle==0 && nonHiggsJet.size()>1){
      vbfJetIndexCombParticle=combinationsNoRepetitionAndOrderDoesNotMatter(2,nonHiggsJetParticle);
      increaseCount(cutFlowMap_particle,"comb vbf part",weight);
    }
    else switchVal_particle=1; 
     
    for(int i=0; i<(int)vbfJetIndexCombParticle.size(); i++){
      //cout<<" pair "<<i<<" Jet 1 pT "<<((Jet*)branchGenJet->At(vbfJetIndexCombParticle[i][0]))->PT<<" 2 "<<((Jet*)branchGenJet->At(vbfJetIndexCombParticle[i][1]))->PT<<endl;
      vbfJetIndexParticle.push_back(make_pair(vbfJetIndexCombParticle[i][0],vbfJetIndexCombParticle[i][1]));
    }
    //cout<<endl;
    //  if( branchMuon->GetEntries() + branchElectron->GetEntries() < 4) continue;

    //if( switchVal_particle==0 && vbfJetIndexParticle.size() > 1) 
    //sort(vbfJetIndexParticle.begin(), vbfJetIndexParticle.end(), [branchGenJet](const pair<int,int> lhs, const pair<int,int> rhs) {
    //	  return fabs((((Jet*)branchGenJet->At(lhs.first))->Eta - ((Jet*)branchGenJet->At(lhs.second))->Eta) ) >
    //	    fabs((((Jet*)branchGenJet->At(rhs.first))->Eta - ((Jet*)branchGenJet->At(rhs.second))->Eta) ) ; 
    //	});

    //cout<<"Here "<<endl;
    //cout<<" number of jet "<<vbfJetIndexParticle.size()<<endl;

    
    //if( switchVal_particle==0 && vbfJetIndexParticle.size() > 1) 
    //sort(vbfJetIndexParticle.begin(), vbfJetIndexParticle.end(), [branchGenJet](const pair<int,int> lhs, const pair<int,int> rhs) {
    //	  return ((((Jet*)branchGenJet->At(lhs.first))->PT)) < ((Jet*)branchGenJet->At(lhs.second))->PT ; 
    //	});
    
    
    int vbfJetsIndexParticleCandidate = -1;
    vector<pair<int,int>> vbfJetIndexParticle2;
    // loop and take first w eta > 2.5
    for (int i=0; i<(int)vbfJetIndexParticle.size(); i++) {
      if( fabs((((Jet*)branchGenJet->At(vbfJetIndexParticle[i].first))->Eta - ((Jet*)branchGenJet->At(vbfJetIndexParticle[i].second))->Eta)) <= 2.5 ) {
	continue; 
      }
      else {
	vbfJetIndexParticle2.push_back(vbfJetIndexParticle.at(i));
	vbfJetsIndexParticleCandidate = i;
      }
    }
   
    vbfJetIndexParticle.clear();
    vbfJetIndexParticle=vbfJetIndexParticle2;
    if(switchVal_particle==0 && vbfJetIndexParticle.size()>0)
      increaseCount(cutFlowMap_particle,"2.5 deltaEta vbf particle",weight);
    else
      switchVal_particle=1; 
  
    //cutFlowMap_particle["2.5 deltaEta vbf particle"] = {cutVal_particle,cutValW_particle}; //last particle cut
    
    //Jet *jet1 = (Jet*) branchGenJet->At(vbfJetIndex[vbfJetsIndexCandidate].first);
    // Jet *jet2 = (Jet*) branchGenJet->At(vbfJetIndex[vbfJetsIndexCandidate].second);
    
    
    // Re sort by eta 
    if( switchVal_particle==0 && vbfJetIndexParticle.size() > 1) 
      sort(vbfJetIndexParticle.begin(), vbfJetIndexParticle.end(), [branchGenJet](const pair<int,int> lhs, const pair<int,int> rhs) {
    	  return fabs((((Jet*)branchGenJet->At(lhs.first))->Eta - ((Jet*)branchGenJet->At(lhs.second))->Eta) ) >
    	    fabs((((Jet*)branchGenJet->At(rhs.first))->Eta - ((Jet*)branchGenJet->At(rhs.second))->Eta) ) ; 
	});
    
    if( switchVal_particle==0){
      jet1_particle = (Jet*) branchGenJet->At(vbfJetIndexParticle[0].first);
      jet2_particle = (Jet*) branchGenJet->At(vbfJetIndexParticle[0].second);
      
      j1_particle=jet1_particle->P4();
      j2_particle=jet2_particle->P4();
      
      if (j1_particle.Eta() > j2_particle.Eta()) {
	jjdeltaPhiparticle = remainder( j1_particle.Phi() - j2_particle.Phi(), 2*M_PI );
      }
      else{
	jjdeltaPhiparticle = remainder( j2_particle.Phi() - j1_particle.Phi(), 2*M_PI );
      }
      
      // double jjdeltaPhiparticle=(j1_particle.Phi() > j2_particle.Phi() ? -1:+1)*TMath::Abs(j2_particle.Phi() - j1_particle.Phi());
      // double jjdeltaEtaparticle=(j1_particle.Eta() > j2_particle.Eta() ? -1:+1)*TMath::Abs(j2_particle.Eta() - j1_particle.Eta());
      // double jjdeltaPhiparticle= (j1_particle.Phi() - j2_particle.Phi());
      jjdeltaEtaparticle= (j1_particle.Eta() - j2_particle.Eta());
      jjdeltaRparticle=sqrt((jjdeltaPhiparticle*jjdeltaPhiparticle)+(jjdeltaEtaparticle*jjdeltaEtaparticle));
    }
//---------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------

// z + leptons
vector <int> goodE_particle_indices; 
for(int i=0; i<(int)branchGenParticle->GetEntries(); i++){
    GenParticle *particle=(GenParticle*) branchGenParticle->At(i);  
    if( particle->Status !=1 ) continue;
    if( fabs(particle->PID) != 11) continue;
    if( fabs(particle->Eta)>= 2.5 || particle->PT <= 5) continue; 
    GenParticle *parent=find_parent(branchGenParticle,particle,particle->PID);
    //if(parent == nullptr) continue; 
    //cout<<" --> Parent "<<parent->PID<<endl;
    if( abs(parent->PID) != 22 &&  abs(parent->PID)!=23  &&  abs(parent->PID)!= 25 ) continue; 
      
    goodE_particle_indices.push_back(i); 
}
 
vector <int> goodMu_particle_indices; 
for(int i=0; i<(int)branchGenParticle->GetEntries(); i++){
    GenParticle *particle=(GenParticle*) branchGenParticle->At(i);  
    if( particle->Status !=1 ) continue; 
    // muons 
    if( fabs(particle->PID) != 13) continue;
    if( fabs(particle->Eta)>=2.5 || particle->PT <= 5.0) continue; 
      
    GenParticle *parent=find_parent(branchGenParticle,particle,particle->PID);
    //if(parent == nullptr) continue; 
    //cout<<" --> Parent "<<parent->PID<<endl;
    if( abs(parent->PID) != 22 &&  abs(parent->PID)!=23  &&  abs(parent->PID)!= 25 ) continue; 

    goodMu_particle_indices.push_back(i); 
 }
 

    // sort the indices by pT ;
    sort(goodE_particle_indices.begin(), goodE_particle_indices.end(), [branchGenParticle](const int& lhs, const int& rhs) {
	return ((GenParticle*) branchGenParticle->At(lhs))->PT > ((GenParticle*) branchGenParticle->At(rhs))->PT;
      });

    // sort the indices by pT ;
    sort(goodMu_particle_indices.begin(), goodMu_particle_indices.end(), [branchGenParticle](const int& lhs, const int& rhs) {
		return ((GenParticle*) branchGenParticle->At(lhs))->PT > ((GenParticle*) branchGenParticle->At(rhs))->PT;
      });


    //cout<<"Particles Electrons "<<endl;
    //for(int i=0;i<(int)goodE_particle_indices.size(); i++)
    //cout<<" e Pt: "<<((GenParticle*)branchGenParticle->At(goodE_particle_indices[i]))->PT<<" eta "<<((GenParticle*)branchGenParticle->At(goodE_particle_indices[i]))->Eta<<" i "<<i<<endl;
    //cout<<"Particles Muons "<<endl;
    //for(int i=0;i<(int)goodMu_particle_indices.size(); i++)
    //cout<<" mu Pt: "<<((GenParticle*)branchGenParticle->At(goodMu_particle_indices[i]))->PT<<" eta "<<((GenParticle*)branchGenParticle->At(goodMu_particle_indices[i]))->Eta<<" i "<<i<<endl;
    
    
    // form pairs 
    // electrons 
    vector< pair<int,int>> elecParticlePairIndices;
    vector< pair<int,int>> elecParticlePairIndicesIn;
    vector <vector<int>> elecParticlePairIndices_;
    if( goodE_particle_indices.size() > 1 )
    elecParticlePairIndices_=combinationsNoRepetitionAndOrderDoesNotMatter(2,goodE_particle_indices);

    // remove all combinations not satisfying criteria 
    for(int i=0;i<(int)elecParticlePairIndices_.size(); i++){
      int elecParticleIndex=elecParticlePairIndices_[i].at(0);
      int elecParticleIndex2=elecParticlePairIndices_[i].at(1);
      
      GenParticle *el1_particle=(GenParticle*) branchGenParticle->At(elecParticleIndex);
      GenParticle *el2_particle=(GenParticle*) branchGenParticle->At(elecParticleIndex2);
      if( el1_particle->Charge == el2_particle->Charge ) continue;
      
      elecParticlePairIndicesIn.push_back(make_pair(elecParticleIndex,elecParticleIndex2));
    }

    elecParticlePairIndices=elecParticlePairIndicesIn;
    
    sort(elecParticlePairIndices.begin(), elecParticlePairIndices.end(), [branchGenParticle](const pair<int,int> lhs, const pair<int,int> rhs) {
    	return fabs(((((GenParticle*) branchGenParticle->At(lhs.first))->P4() + ((GenParticle*) branchGenParticle->At(lhs.second))->P4())).M() -91 ) <
	  fabs( ((((GenParticle*) branchGenParticle->At(rhs.first))->P4() + ((GenParticle*) branchGenParticle->At(rhs.second))->P4()).M()) -91 ) ; 
      });

    remove_overlaps(elecParticlePairIndices);

    // muons 
    vector< pair<int,int>> muParticlePairIndicesIn;
    vector< pair<int,int>> muParticlePairIndices;
    vector <vector<int>> muParticlePairIndices_;
    if(goodMu_particle_indices.size() > 1) 
      muParticlePairIndices_=combinationsNoRepetitionAndOrderDoesNotMatter(2,goodMu_particle_indices);

    // remove all combinations not satifing criteria 
    for(int i=0; i<(int)muParticlePairIndices_.size(); i++){
      int elecParticleIndex=muParticlePairIndices_[i].at(0);
      int elecParticleIndex2=muParticlePairIndices_[i].at(1);

      GenParticle *el1_particle=(GenParticle*) branchGenParticle->At(elecParticleIndex);
      GenParticle *el2_particle=(GenParticle*) branchGenParticle->At(elecParticleIndex2);
      if( el1_particle->Charge ==  el2_particle->Charge ) continue;
      muParticlePairIndicesIn.push_back(make_pair(elecParticleIndex,elecParticleIndex2));
    }

    muParticlePairIndices=muParticlePairIndicesIn;

    sort(muParticlePairIndices.begin(), muParticlePairIndices.end(), [branchGenParticle](const pair<int,int> lhs, const pair<int,int> rhs) {
    	return fabs(((((GenParticle*) branchGenParticle->At(lhs.first))->P4() + ((GenParticle*) branchGenParticle->At(lhs.second))->P4())).M() -91 ) <
	  fabs( ((((GenParticle*) branchGenParticle->At(rhs.first))->P4() + ((GenParticle*) branchGenParticle->At(rhs.second))->P4()).M()) -91 ) ; 
      });
    
    //sort(muParticlePairIndices.begin(),muParticlePairIndices.end(), [branchGenParticle]( pair<int,int>   & lhs,  pair<int,int>   & rhs) {
    //	int index11_particle=(lhs).first;
    //	    int index12_particle=(lhs).second;
    //	    int index21_particle=(rhs).first;
    //	int index22_particle=(rhs).second;
    //	    return fabs(((((GenParticle*) branchGenParticle->At(index11_particle))->P4() + ((GenParticle*) branchGenParticle->At(index12_particle))->P4())).M() - 91) <
    //	  fabs( ((((GenParticle*) branchGenParticle->At(index21_particle))->P4() + ((GenParticle*) branchGenParticle->At(index22_particle))->P4()).M()) -91);
    //});

    remove_overlaps(muParticlePairIndices);

    //increaseCount(cutFlowMap_particle,"at least two lep pairs",weight);
    
    // GB: This is a cut !!! 
    // order quads by mZ1 and mZ2 
    //if( muParticlePairIndices.size() + elecParticlePairIndices.size() < 2) continue; // no candidate found;
    //if( muParticlePairIndices.size() + elecParticlePairIndices.size() < 2) {
    //switchVal_particle=1; 
    //}
    //else {
    //increaseCount(cutFlowMap_particle,"at least two lep pairs",weight);
    ////cutFlowMap_particle["at least two lep pairs"] = {cutVal_particle,cutValW_particle};
    //}

    int thisParticleEventType=-1; 

    vector<pair<int,pair<int,int>>> ParticlePairIndices; // 0 for electron 1 for muon
    for(int i=0; i<(int) elecParticlePairIndices.size(); i++)
      ParticlePairIndices.push_back(make_pair(0,elecParticlePairIndices.at(i)));
    for(int i=0; i<(int) muParticlePairIndices.size(); i++)
      ParticlePairIndices.push_back(make_pair(1,muParticlePairIndices.at(i)));
    
    // sort all of the indices by closeness to mZ
    sort(ParticlePairIndices.begin(), ParticlePairIndices.end(), [branchGenParticle]  ( const pair<int,pair<int,int>>  lhs , const pair<int,pair<int,int>>   rhs ){
	
	double mass1_particle=0;
    	double mass2_particle=0;
	
      if( lhs.first==0 ) mass1_particle = (((GenParticle*) branchGenParticle->At(lhs.second.first))->P4() + ((GenParticle*) branchGenParticle->At(lhs.second.second))->P4()).M();
      else mass1_particle= ((((GenParticle*) branchGenParticle->At(lhs.second.first))->P4() + ((GenParticle*) branchGenParticle->At(lhs.second.second))->P4())).M();
          
      if( rhs.first==0 ) mass2_particle = (((GenParticle*) branchGenParticle->At(rhs.second.first))->P4() + ((GenParticle*) branchGenParticle->At(rhs.second.second))->P4()).M();
      else mass2_particle = (((GenParticle*) branchGenParticle->At(rhs.second.first))->P4() + ((GenParticle*) branchGenParticle->At(rhs.second.second))->P4()).M();
          
      return fabs(mass1_particle -91.0 )  <  fabs(mass2_particle -91.0); 
    });
    
    //cout<<"Particle pair indices "<<ParticlePairIndices.size()<<endl;
    //for( int i=0; i<(int)ParticlePairIndices.size(); i++){
    //if(ParticlePairIndices[i].first==0) 
    //	cout<<" "<<"electrons "<< (((GenParticle*) branchGenParticle->At(ParticlePairIndices[i].second.first))->P4() + ((GenParticle*) branchGenParticle->At(ParticlePairIndices[i].second.second))->P4()).M()<<" "<<i<<endl;
    //else
    //	cout<<" "<<"muons "<< (((GenParticle*) branchGenParticle->At(ParticlePairIndices[i].second.first))->P4() + ((GenParticle*) branchGenParticle->At(ParticlePairIndices[i].second.second))->P4()).M()<<" "<<i<<endl;
    //}
    //cout<<endl;

    if( ParticlePairIndices.size() <2)  switchVal_particle=1;
    else if (switchVal_particle ==0)  increaseCount(cutFlowMap_particle,"OSFL",weight);
        
    if(switchVal_particle==0){
      if( ParticlePairIndices[0].first == 1 && ParticlePairIndices[1].first == 1) thisParticleEventType=0;
      else if( ParticlePairIndices[0].first == 0 && ParticlePairIndices[1].first == 0) thisParticleEventType=1;
      else if( ParticlePairIndices[0].first == 1 && ParticlePairIndices[1].first == 0) thisParticleEventType=2;
      else if( ParticlePairIndices[0].first == 0 && ParticlePairIndices[1].first == 1) thisParticleEventType=3;
    }
  // case 4mu 
  if( thisParticleEventType==0 ) {
    
    // take first two muons
    GenParticle *muon1_particle = (GenParticle*) branchGenParticle->At( ParticlePairIndices[0].second.first);
    GenParticle *muon2_particle = (GenParticle*) branchGenParticle->At( ParticlePairIndices[0].second.second);
    // take first two muons
    GenParticle *muon3_particle = (GenParticle*) branchGenParticle->At( ParticlePairIndices[1].second.first);
    GenParticle *muon4_particle = (GenParticle*) branchGenParticle->At( ParticlePairIndices[1].second.second);
    
    q1_particle = muon1_particle->Charge;
    q2_particle = muon2_particle->Charge;
    q3_particle = muon3_particle->Charge;
    q4_particle = muon4_particle->Charge;

    l1_particle=muon1_particle->P4();
    l2_particle=muon2_particle->P4();
    l3_particle=muon3_particle->P4();
    l4_particle=muon4_particle->P4(); 

  }

  // case 4e
  else if( thisParticleEventType==1 ) {
    
    // take first two electrons
    GenParticle *muon1_particle = (GenParticle*) branchGenParticle->At( ParticlePairIndices[0].second.first);
    GenParticle *muon2_particle = (GenParticle*) branchGenParticle->At( ParticlePairIndices[0].second.second);
    // take first two electrons
    GenParticle *muon3_particle = (GenParticle*) branchGenParticle->At( ParticlePairIndices[1].second.first);
    GenParticle *muon4_particle = (GenParticle*) branchGenParticle->At( ParticlePairIndices[1].second.second);

    q1_particle = muon1_particle->Charge;
    q2_particle = muon2_particle->Charge;
    q3_particle = muon3_particle->Charge;
    q4_particle = muon4_particle->Charge;

    l1_particle=muon1_particle->P4();
    l2_particle=muon2_particle->P4();
    l3_particle=muon3_particle->P4();
    l4_particle=muon4_particle->P4(); 
  }

  // case 2mu2e 
  if( thisParticleEventType==2 ) {
    
    // take first two muons
    GenParticle *muon1_particle = (GenParticle*) branchGenParticle->At( ParticlePairIndices[0].second.first);
    GenParticle *muon2_particle = (GenParticle*) branchGenParticle->At( ParticlePairIndices[0].second.second);
    // take first two electrons
    GenParticle *muon3_particle = (GenParticle*) branchGenParticle->At( ParticlePairIndices[1].second.first);
    GenParticle *muon4_particle = (GenParticle*) branchGenParticle->At( ParticlePairIndices[1].second.second);
  
    q1_particle = muon1_particle->Charge;
    q2_particle = muon2_particle->Charge;
    q3_particle = muon3_particle->Charge;
    q4_particle = muon4_particle->Charge;

    l1_particle=muon1_particle->P4();
    l2_particle=muon2_particle->P4();
    l3_particle=muon3_particle->P4();
    l4_particle=muon4_particle->P4(); 
  }

  // case 2e2mu 
   else if( thisParticleEventType==3 ) {
    
    // take first two electrons
    GenParticle *muon1_particle = (GenParticle*) branchGenParticle->At( ParticlePairIndices[0].second.first);
    GenParticle *muon2_particle = (GenParticle*) branchGenParticle->At( ParticlePairIndices[0].second.second);
    // take first two muons
    GenParticle *muon3_particle = (GenParticle*) branchGenParticle->At( ParticlePairIndices[1].second.first);
    GenParticle *muon4_particle = (GenParticle*) branchGenParticle->At( ParticlePairIndices[1].second.second);

    q1_particle = muon1_particle->Charge;
    q2_particle = muon2_particle->Charge;
    q3_particle = muon3_particle->Charge;
    q4_particle = muon4_particle->Charge;

    l1_particle=muon1_particle->P4();
    l2_particle=muon2_particle->P4();
    l3_particle=muon3_particle->P4();
    l4_particle=muon4_particle->P4(); 

   }

    z1_particle=l1_particle + l2_particle;
    z2_particle=l3_particle + l4_particle;
    fourl_particle=l1_particle + l2_particle + l3_particle + l4_particle;


    l1l2deltaPhiparticle=(l1_particle.Phi() > l2_particle.Phi() ? -1:+1)*TMath::Abs(l2_particle.Phi() - l1_particle.Phi());
    l3l4deltaPhiparticle=(l3_particle.Phi() > l4_particle.Phi() ? -1:+1)*TMath::Abs(l4_particle.Phi() - l3_particle.Phi());

    l1l2deltaEtaparticle=(l1_particle.Eta() > l2_particle.Eta() ? -1:+1)*TMath::Abs(l2_particle.Eta() - l1_particle.Eta());
    l3l4deltaEtaparticle=(l3_particle.Eta() > l4_particle.Eta() ? -1:+1)*TMath::Abs(l4_particle.Eta() - l3_particle.Eta());

    l1l2deltaRparticle=sqrt((l1l2deltaPhiparticle*l1l2deltaPhiparticle)+(l1l2deltaEtaparticle*l1l2deltaEtaparticle));
    l3l4deltaRparticle=sqrt((l3l4deltaPhiparticle*l3l4deltaPhiparticle)+(l3l4deltaEtaparticle*l3l4deltaEtaparticle));

    l1cosThetaparticle=l1_particle.CosTheta();
    l2cosThetaparticle=l2_particle.CosTheta();
    l3cosThetaparticle=l3_particle.CosTheta();
    l4cosThetaparticle=l4_particle.CosTheta();
    fourlcosThetaparticle=fourl_particle.CosTheta();

    l1_particle.Boost(-z1_particle.BoostVector());
    l2_particle.Boost(-z1_particle.BoostVector());
    l3_particle.Boost(-z2_particle.BoostVector());
    l4_particle.Boost(-z2_particle.BoostVector());
    l1cosThetaBoostparticle=l1_particle.CosTheta();
    l2cosThetaBoostparticle=l2_particle.CosTheta();
    l3cosThetaBoostparticle=l3_particle.CosTheta();
    l4cosThetaBoostparticle=l4_particle.CosTheta();
    fourlcosThetaBoostparticle=fourl_particle.CosTheta();
    l1_particle.Boost(z1_particle.BoostVector());
    l2_particle.Boost(z1_particle.BoostVector());
    l3_particle.Boost(z2_particle.BoostVector());
    l4_particle.Boost(z2_particle.BoostVector());

  // collins soper frame
    l1l2CScosThetaparticle=(q1_particle > q2_particle ? -1:+1)*TMath::Abs(2*(l2_particle.Pz()*l1_particle.E()-l1_particle.Pz()*l2_particle.E())/(z1_particle.M()*sqrt(z1_particle.M()*z1_particle.M()+z1_particle.Pt()*z1_particle.Pt())));
    l3l4CScosThetaparticle=(q3_particle > q4_particle ? -1:+1)*TMath::Abs(2*(l4_particle.Pz()*l3_particle.E()-l3_particle.Pz()*l4_particle.E())/(z2_particle.M()*sqrt(z2_particle.M()*z2_particle.M()+z2_particle.Pt()*z2_particle.Pt())));

    zzdeltaPhiparticle=(z1_particle.Phi() > z2_particle.Phi() ? -1:+1)*TMath::Abs(z2_particle.Phi() - z1_particle.Phi());
    zzdeltaEtaparticle=(z1_particle.Eta() > z2_particle.Eta() ? -1:+1)*TMath::Abs(z2_particle.Eta() - z1_particle.Eta());
    zzdeltaRparticle=sqrt((zzdeltaPhiparticle*zzdeltaPhiparticle)+(zzdeltaEtaparticle*zzdeltaEtaparticle));

// neutrinos 
// jets we use the GenJet container

    //cout<<"Done particle level "<<endl;
    
//------------------------------------------------------------------------------------------------------------------------------------------------------------
// parton
//------------------------------------------------------------------------------------------------------------------------------------------------------------

// higgs

    increaseCount(cutFlowMap_parton,"initial parton",weight);

    int switchVal_parton = 0;

    bool HiggsRecord=FillHiggsTruthRecord(branchGenParticle,h_parton,b1_parton,b2_parton,j1_parton,j2_parton);
    
    if(switchVal_parton == 0 && HiggsRecord) { 
      increaseCount(cutFlowMap_parton,"Higgs Candidate",weight);
    } else  switchVal_parton = 1;
    
    if(HiggsRecord){
      if (b1_parton.Eta() > b2_parton.Eta()) {
	bbdeltaPhiparton = remainder( b1_parton.Phi() - b2_parton.Phi(), 2*M_PI );
	bbdeltaEtaparton= b1_parton.Eta() - b2_parton.Eta();
      }
      else{
	bbdeltaPhiparton = remainder( b2_parton.Phi() - b1_parton.Phi(), 2*M_PI );
	bbdeltaEtaparton= b2_parton.Eta() - b1_parton.Eta();
    }
      
      // double bbdeltaPhiparton=(b1_parton.Phi() > b2_parton.Phi() ? -1:+1)*TMath::Abs(b2_parton.Phi() - b1_parton.Phi());
      // double bbdeltaEtaparton=(b1_parton.Eta() > b2_parton.Eta() ? -1:+1)*TMath::Abs(b2_parton.Eta() - b1_parton.Eta());
      // double bbdeltaPhiparton= (b1_parton.Phi() - b2_parton.Phi());
      // double bbdeltaEtaparton= (b1_parton.Eta() - b2_parton.Eta());
      bbdeltaRparton=sqrt((bbdeltaPhiparton*bbdeltaPhiparton)+(bbdeltaEtaparton*bbdeltaEtaparton));

      if (j1_parton.Eta() > j2_parton.Eta()) {
	jjdeltaPhiparton = remainder( j1_parton.Phi() - j2_parton.Phi(), 2*M_PI );
	jjdeltaEtaparton= j1_parton.Eta() - j2_parton.Eta();
      }
      else{
	jjdeltaPhiparton = remainder( j2_parton.Phi() - j1_parton.Phi(), 2*M_PI );
	jjdeltaEtaparton= j2_parton.Eta() - j1_parton.Eta();
      }
    }
    
// z + leptons
   
    vector <int> genZBosons; 
    for(int i=0; i<(int)branchGenParticle->GetEntries(); i++){
      GenParticle *particle=(GenParticle*)branchGenParticle->At(i);
      int d1_pid = 9999;
      if (particle->D1 != -1) {
        GenParticle *daughter1 = (GenParticle*) branchGenParticle->At(particle->D1);
        d1_pid = daughter1 -> PID;
      }
      // find a Z.
      if( particle->PID == 23 && d1_pid != 23)
        genZBosons.push_back(i);
    }

    if(genZBosons.size() >1){
      sort(genZBosons.begin(), genZBosons.end(), [branchGenParticle](const int& lhs, const int& rhs) {
	  return ((GenParticle*)branchGenParticle->At(lhs))->P4().M() > (((GenParticle*)branchGenParticle->At(rhs))->P4()).M(); 
	});
    }

    bool foundZZ=genZBosons.size() > 1;
    if(switchVal_parton == 0 && foundZZ)
      increaseCount(cutFlowMap_parton,"ZZ parton",weight);
    else  switchVal_parton = 1;
  
     if(foundZZ){
    z1_parton = ((GenParticle*)branchGenParticle->At(genZBosons[0]))->P4(); 
    z2_parton = ((GenParticle*)branchGenParticle->At(genZBosons[1]))->P4(); 
    
    vector <int> Z1children;
    Z1children.push_back( ((GenParticle*)branchGenParticle->At(genZBosons[0]))->D1 ); 
    Z1children.push_back( ((GenParticle*)branchGenParticle->At(genZBosons[0]))->D2 ); 
    
    vector <int> Z2children; 
    Z2children.push_back( ((GenParticle*)branchGenParticle->At(genZBosons[1]))->D1 );
    Z2children.push_back( ((GenParticle*)branchGenParticle->At(genZBosons[1]))->D2 ); 
    
    l1_parton = ((GenParticle*) branchGenParticle->At(Z1children[0]))->P4();
    l2_parton = ((GenParticle*) branchGenParticle->At(Z1children[1]))->P4();
    l3_parton = ((GenParticle*) branchGenParticle->At(Z2children[0]))->P4();
    l4_parton = ((GenParticle*) branchGenParticle->At(Z2children[1]))->P4();
    
    q1_parton = ((GenParticle*) branchGenParticle->At(Z1children[0]))->Charge;
    q2_parton = ((GenParticle*) branchGenParticle->At(Z1children[1]))->Charge;
    q3_parton = ((GenParticle*) branchGenParticle->At(Z2children[0]))->Charge;
    q4_parton = ((GenParticle*) branchGenParticle->At(Z2children[1]))->Charge;
    
    //z1_parton=l1_parton + l2_parton;
    //z2_parton=l3_parton + l4_parton;
    fourl_parton=l1_parton + l2_parton + l3_parton + l4_parton;
    
    l1l2deltaPhiparton=(l1_parton.Phi() > l2_parton.Phi() ? -1:+1)*TMath::Abs(l2_parton.Phi() - l1_parton.Phi());
    l3l4deltaPhiparton=(l3_parton.Phi() > l4_parton.Phi() ? -1:+1)*TMath::Abs(l4_parton.Phi() - l3_parton.Phi());
    
    l1l2deltaEtaparton=(l1_parton.Eta() > l2_parton.Eta() ? -1:+1)*TMath::Abs(l2_parton.Eta() - l1_parton.Eta());
    l3l4deltaEtaparton=(l3_parton.Eta() > l4_parton.Eta() ? -1:+1)*TMath::Abs(l4_parton.Eta() - l3_parton.Eta());
    
    l1l2deltaRparton=sqrt((l1l2deltaPhiparton*l1l2deltaPhiparton)+(l1l2deltaEtaparton*l1l2deltaEtaparton));
    l3l4deltaRparton=sqrt((l3l4deltaPhiparton*l3l4deltaPhiparton)+(l3l4deltaEtaparton*l3l4deltaEtaparton));
    
    l1cosThetaparton=l1_parton.CosTheta();
    l2cosThetaparton=l2_parton.CosTheta();
    l3cosThetaparton=l3_parton.CosTheta();
    l4cosThetaparton=l4_parton.CosTheta();
    fourlcosThetaparton=fourl_parton.CosTheta();

    l1_parton.Boost(-z1_parton.BoostVector());
    l2_parton.Boost(-z1_parton.BoostVector());
    l3_parton.Boost(-z2_parton.BoostVector());
    l4_parton.Boost(-z2_parton.BoostVector());
    l1cosThetaBoostparton=l1_parton.CosTheta();
    l2cosThetaBoostparton=l2_parton.CosTheta();
    l3cosThetaBoostparton=l3_parton.CosTheta();
    l4cosThetaBoostparton=l4_parton.CosTheta();
    fourlcosThetaBoostparton=fourl_parton.CosTheta();
    l1_parton.Boost(z1_parton.BoostVector());
    l2_parton.Boost(z1_parton.BoostVector());
    l3_parton.Boost(z2_parton.BoostVector());
    l4_parton.Boost(z2_parton.BoostVector());

  // collins soper frame
    double l1l2CScosThetaparton=(q1_parton > q2_parton ? -1:+1)*TMath::Abs(2*(l2_parton.Pz()*l1_parton.E()-l1_parton.Pz()*l2_parton.E())/(z1_parton.M()*sqrt(z1_parton.M()*z1_parton.M()+z1_parton.Pt()*z1_parton.Pt())));
    double l3l4CScosThetaparton=(q3_parton > q4_parton ? -1:+1)*TMath::Abs(2*(l4_parton.Pz()*l3_parton.E()-l3_parton.Pz()*l4_parton.E())/(z2_parton.M()*sqrt(z2_parton.M()*z2_parton.M()+z2_parton.Pt()*z2_parton.Pt())));

    double zzdeltaPhiparton=(z1_parton.Phi() > z2_parton.Phi() ? -1:+1)*TMath::Abs(z2_parton.Phi() - z1_parton.Phi());
    double zzdeltaEtaparton=(z1_parton.Eta() > z2_parton.Eta() ? -1:+1)*TMath::Abs(z2_parton.Eta() - z1_parton.Eta());
    double zzdeltaRparton=sqrt((zzdeltaPhiparton*zzdeltaPhiparton)+(zzdeltaEtaparton*zzdeltaEtaparton));
  }
    // Print cutflow intermideiate

    if( entry % 1000 == 0 ){
      cout<<"Processed "<<entry<< " / " <<numberOfEntries <<" "<< entry/ numberOfEntries *100 <<" %"<<endl;
      cout<<" Reco CutF Flow "<<endl;
      PrintCutFlow(cutFlowMap_reco,cutList_reco,"Reco");
      cout<<" Particle  CutF Flow "<<endl;
      PrintCutFlow(cutFlowMap_particle,cutList_particle, "Particle");
      cout<<" Parton  CutF Flow "<<endl;
      PrintCutFlow(cutFlowMap_parton,cutList_parton, "Parton");
    }
//------------------------------------------------------------------------------------------------------------------------------------------------------------
// fill histos
//------------------------------------------------------------------------------------------------------------------------------------------------------------

    

// higgs + b
    // reco
    if(switchVal_reco==0){
    // pT + m
      hHpTreco -> Fill(h_reco.Pt(),weight);
      hHmreco -> Fill(h_reco.M(),weight); 
      hb1pTreco -> Fill(b1_reco.Pt(),weight);
      hb1mreco -> Fill(b1_reco.M(),weight); 
      hb2pTreco -> Fill(b2_reco.Pt(),weight);
      hb2mreco -> Fill(b2_reco.M(),weight); 
      // phi
      hHphireco -> Fill(h_reco.Phi(),weight);
      hb1phireco -> Fill(b1_reco.Phi(),weight);
      hb2phireco -> Fill(b2_reco.Phi(),weight);
      hbbdeltaPhireco -> Fill(bbdeltaPhireco,weight);
      // eta
      hHetareco -> Fill(h_reco.Eta(),weight);
      hb1etareco -> Fill(b1_reco.Eta(),weight);
      hb2etareco -> Fill(b2_reco.Eta(),weight);
      hbbdeltaEtareco -> Fill(bbdeltaEtareco,weight);
      // R = sqrt[phi^2 + eta^2]^1/2
      hHRreco -> Fill(sqrt((h_reco.Phi()*h_reco.Phi())+(h_reco.Eta()*h_reco.Eta())),weight);
      hb1Rreco -> Fill(sqrt((b1_reco.Phi()*b1_reco.Phi())+(b1_reco.Eta()*b1_reco.Eta())),weight);
      hb2Rreco -> Fill(sqrt((b2_reco.Phi()*b2_reco.Phi())+(b2_reco.Eta()*b2_reco.Eta())),weight);
      hbbdeltaRreco -> Fill(bbdeltaRreco,weight);
    }
    // particle
    if(switchVal_particle==0){
      // pT + m
      hHpTparticle -> Fill(h_particle.Pt(), weight);
      hHmparticle -> Fill(h_particle.M(), weight);
      hb1pTparticle -> Fill(b1_particle.Pt(),weight);
      hb1mparticle -> Fill(b1_particle.M(),weight); 
      hb2pTparticle -> Fill(b2_particle.Pt(),weight);
      hb2mparticle -> Fill(b2_particle.M(),weight); 
      // phi
      hHphiparticle -> Fill(h_particle.Phi(),weight);
      hb1phiparticle -> Fill(b1_particle.Phi(),weight);
      hb2phiparticle -> Fill(b2_particle.Phi(),weight);
      hbbdeltaPhiparticle->Fill(bbdeltaPhiparticle,weight);
      // eta
      hHetaparticle -> Fill(h_particle.Eta(),weight);
      hb1etaparticle -> Fill(b1_particle.Eta(),weight);
      hb2etaparticle -> Fill(b2_particle.Eta(),weight);
      hbbdeltaEtaparticle -> Fill(bbdeltaEtaparticle,weight);
      // R = sqrt[phi^2 + eta^2]^1/2
      hHRparticle -> Fill(sqrt((h_particle.Phi()*h_particle.Phi())+(h_particle.Eta()*h_particle.Eta())),weight);
      hb1Rparticle -> Fill(sqrt((b1_particle.Phi()*b1_particle.Phi())+(b1_particle.Eta()*b1_particle.Eta())),weight);
      hb2Rparticle -> Fill(sqrt((b2_particle.Phi()*b2_particle.Phi())+(b2_particle.Eta()*b2_particle.Eta())),weight);
      hbbdeltaRparticle -> Fill(bbdeltaRparticle,weight);
    }
    // parton
    if(switchVal_parton==0){
    // pT + m
      hHpTparton -> Fill(h_parton.Pt(), weight);
      hHmparton -> Fill(h_parton.M(), weight);
      hb1pTparton-> Fill(b1_parton.Pt(),weight);
      hb1mparton -> Fill(b1_parton.M(),weight); 
      hb2pTparton -> Fill(b2_parton.Pt(),weight);
      hb2mparton -> Fill(b2_parton.M(),weight); 
      // phi
      hHphiparton -> Fill(h_parton.Phi(),weight);
      hb1phiparton -> Fill(b1_parton.Phi(),weight);
      hb2phiparton -> Fill(b2_parton.Phi(),weight);
      hbbdeltaPhiparton->Fill(bbdeltaPhiparton,weight);
      // eta
      hHetaparton -> Fill(h_parton.Eta(),weight);
      hb1etaparton -> Fill(b1_parton.Eta(),weight);
      hb2etaparton -> Fill(b2_parton.Eta(),weight);
      hbbdeltaEtaparton -> Fill(bbdeltaEtaparton,weight);
      // R = sqrt[phi^2 + eta^2]^1/2
      hHRparton -> Fill(sqrt((h_parton.Phi()*h_parton.Phi())+(h_parton.Eta()*h_parton.Eta())),weight);
      hb1Rparton -> Fill(sqrt((b1_parton.Phi()*b1_parton.Phi())+(b1_parton.Eta()*b1_parton.Eta())),weight);
      hb2Rparton -> Fill(sqrt((b2_parton.Phi()*b2_parton.Phi())+(b2_parton.Eta()*b2_parton.Eta())),weight);
      hbbdeltaRparton -> Fill(bbdeltaRparton,weight);
    }

// jets
    // reco
      if(switchVal_reco==0){
      // pT + m
      hjjpTreco->Fill(j1_reco.Pt()+j2_reco.Pt(),weight);
        hj1pTreco->Fill(j1_reco.Pt(),weight);
        hj2pTreco->Fill(j2_reco.Pt(),weight);
        // phi
        hj1phireco->Fill(j1_reco.Phi(),weight); 
        hj2phireco->Fill(j2_reco.Phi(),weight); 
        hjjdeltaPhireco->Fill(jjdeltaPhireco,weight); 
        // eta
        hj1etareco->Fill(j1_reco.Eta(),weight); 
        hj2etareco->Fill(j2_reco.Eta(),weight);
        hjjdeltaEtareco->Fill(jjdeltaEtareco, weight);
        // R = sqrt[phi^2 + eta^2]^1/2
        hj1Rreco -> Fill(sqrt((j1_reco.Phi()*j1_reco.Phi())+(j1_reco.Eta()*j1_reco.Eta())),weight);
        hj2Rreco -> Fill(sqrt((j2_reco.Phi()*j2_reco.Phi())+(j2_reco.Eta()*j2_reco.Eta())),weight);
        hjjdeltaRreco -> Fill(jjdeltaRreco,weight);
      }

    // particle
      // pT + m
      if(switchVal_particle==0){
        hjjpTparticle->Fill(j1_particle.Pt()+j2_particle.Pt(),weight);
        hj1pTparticle->Fill(j1_particle.Pt(),weight);
        hj2pTparticle->Fill(j2_particle.Pt(),weight);
        // phi
        hj1phiparticle->Fill(j1_particle.Phi(),weight); 
        hj2phiparticle->Fill(j2_particle.Phi(),weight); 
        hjjdeltaPhiparticle->Fill(jjdeltaPhiparticle,weight); 
        // eta
        hj1etaparticle->Fill(j1_particle.Eta(),weight); 
        hj2etaparticle->Fill(j2_particle.Eta(),weight);
        hjjdeltaEtaparticle->Fill(jjdeltaEtaparticle, weight);
        // R = sqrt[phi^2 + eta^2]^1/2
        hj1Rparticle -> Fill(sqrt((j1_particle.Phi()*j1_particle.Phi())+(j1_particle.Eta()*j1_particle.Eta())),weight);
        hj2Rparticle -> Fill(sqrt((j2_particle.Phi()*j2_particle.Phi())+(j2_particle.Eta()*j2_particle.Eta())),weight);
        hjjdeltaRparticle -> Fill(jjdeltaRparticle,weight);
      }

// z
      //reco
      if(switchVal_reco==0){
        // pT + m
        hz1pTreco->Fill(z1_reco.Pt(),weight);
        hz2pTreco->Fill(z2_reco.Pt(),weight);
        hz1mreco->Fill(z1_reco.M(),weight);
        hz2mreco->Fill(z2_reco.M(),weight);
        // cos
        hz1cosThetareco->Fill(z1_reco.CosTheta(),weight);
        hz2cosThetareco->Fill(z2_reco.CosTheta(),weight);
        // phi
        hz1phireco->Fill(z1_reco.Phi(),weight); 
        hz2phireco->Fill(z2_reco.Phi(),weight); 
        hzzdeltaPhireco->Fill(zzdeltaPhireco,weight); 
        // eta
        hz1etareco->Fill(z1_reco.Eta(),weight); 
        hz2etareco->Fill(z2_reco.Eta(),weight);
        hzzdeltaEtareco->Fill(zzdeltaEtareco, weight);
        // R = sqrt[phi^2 + eta^2]^1/2
        hz1Rreco -> Fill(sqrt((z1_reco.Phi()*z1_reco.Phi())+(z1_reco.Eta()*z1_reco.Eta())),weight);
        hz2Rreco -> Fill(sqrt((z2_reco.Phi()*z2_reco.Phi())+(z2_reco.Eta()*z2_reco.Eta())),weight);
        hzzdeltaRreco -> Fill(zzdeltaRreco,weight);
      }

    // particle
      if(switchVal_particle==0){
        // pT + m
        hz1pTparticle->Fill(z1_particle.Pt(), weight);
        hz2pTparticle->Fill(z2_particle.Pt(), weight);
        hz1mparticle -> Fill(z1_particle.M(), weight);
        hz2mparticle -> Fill(z2_particle.M(), weight);
        // phi
        hz1phiparticle->Fill(z1_particle.Phi(),weight); 
        hz2phiparticle->Fill(z2_particle.Phi(),weight); 
        hzzdeltaPhiparticle->Fill(zzdeltaPhiparticle,weight); 
        // eta
        hz1etaparticle->Fill(z1_particle.Eta(),weight); 
        hz2etaparticle->Fill(z2_particle.Eta(),weight);
        hzzdeltaEtaparticle->Fill(zzdeltaEtaparticle, weight);
        // R = sqrt[phi^2 + eta^2]^1/2
        hz1Rparticle -> Fill(sqrt((z1_particle.Phi()*z1_particle.Phi())+(z1_particle.Eta()*z1_particle.Eta())),weight);
        hz2Rparticle -> Fill(sqrt((z2_particle.Phi()*z2_particle.Phi())+(z2_particle.Eta()*z2_particle.Eta())),weight);
        hzzdeltaRparticle -> Fill(zzdeltaRparticle,weight);
      }

      // parton
      // pT + m
      if(switchVal_parton==0 ){
      hz1pTparton->Fill(z1_parton.Pt(),weight);
      hz2pTparton->Fill(z2_parton.Pt(),weight);
      hz1mparton->Fill(z1_parton.M(),weight);
      hz2mparton->Fill(z2_parton.M(),weight);
      // phi
      hz1phiparton->Fill(z1_parton.Phi(),weight); 
      hz2phiparton->Fill(z2_parton.Phi(),weight); 
      hzzdeltaPhiparton->Fill(zzdeltaPhiparton,weight); 
      // eta
      hz1etaparton->Fill(z1_parton.Eta(),weight); 
      hz2etaparton->Fill(z2_parton.Eta(),weight);
      hzzdeltaEtaparton->Fill(zzdeltaEtaparton, weight);
      // R = sqrt[phi^2 + eta^2]^1/2
      hz1Rparton -> Fill(sqrt((z1_parton.Phi()*z1_parton.Phi())+(z1_parton.Eta()*z1_parton.Eta())),weight);
      hz2Rparton -> Fill(sqrt((z2_parton.Phi()*z2_parton.Phi())+(z2_parton.Eta()*z2_parton.Eta())),weight);
      hzzdeltaRparton -> Fill(zzdeltaRparton,weight);
      }

//leptons
    // reco
	if(switchVal_reco==0){
        // pT + m
        hl1pTreco->Fill(l1_reco.Pt(),weight);
        hl2pTreco->Fill(l2_reco.Pt(),weight);
        hl3pTreco->Fill(l3_reco.Pt(),weight);
        hl4pTreco->Fill(l4_reco.Pt(),weight);
        // phi
        hl1phireco->Fill(l1_reco.Phi() ,weight);
        hl2phireco->Fill(l2_reco.Phi() ,weight);
        hl3phireco->Fill(l3_reco.Phi() ,weight);
        hl4phireco->Fill(l4_reco.Phi() ,weight);
        hl1l2deltaPhireco->Fill(l1l2deltaPhireco ,weight);
        hl3l4deltaPhireco->Fill(l3l4deltaPhireco ,weight);
        hl1l2deltaPhiBoostreco->Fill(l1l2deltaPhiBoostreco ,weight);
        hl3l4deltaPhiBoostreco->Fill(l3l4deltaPhiBoostreco ,weight);
        // eta
        hl1etareco->Fill(l1_reco.Eta() ,weight);
        hl2etareco->Fill(l2_reco.Eta() ,weight);
        hl3etareco->Fill(l3_reco.Eta() ,weight);
        hl4etareco->Fill(l4_reco.Eta() ,weight);
        hl1l2deltaEtareco->Fill(l1l2deltaEtareco ,weight);
        hl3l4deltaEtareco->Fill(l3l4deltaEtareco ,weight);
        hl1l2deltaEtaBoostreco->Fill(l1l2deltaEtaBoostreco ,weight);
        hl3l4deltaEtaBoostreco->Fill(l3l4deltaEtaBoostreco ,weight);
        // R
        hl1Rreco->Fill(sqrt((l1_reco.Phi()*l1_reco.Phi())+(l1_reco.Eta()*l1_reco.Eta())) ,weight);
        hl2Rreco->Fill(sqrt((l2_reco.Phi()*l2_reco.Phi())+(l2_reco.Eta()*l2_reco.Eta())) ,weight);
        hl3Rreco->Fill(sqrt((l3_reco.Phi()*l3_reco.Phi())+(l3_reco.Eta()*l3_reco.Eta())) ,weight);
        hl4Rreco->Fill(sqrt((l4_reco.Phi()*l4_reco.Phi())+(l4_reco.Eta()*l4_reco.Eta())) ,weight);
        hl1l2deltaRreco->Fill(l1l2deltaRreco ,weight);
        hl3l4deltaRreco->Fill(l3l4deltaRreco ,weight);
        // cos
        hl1cosThetareco->Fill(l1cosThetareco,weight);
        hl2cosThetareco->Fill(l2cosThetareco,weight);
        hl3cosThetareco->Fill(l3cosThetareco,weight);
        hl4cosThetareco->Fill(l4cosThetareco,weight);
        hfourlcosThetareco->Fill(fourlcosThetareco,weight);
        hl1cosThetaBoostreco->Fill(l1cosThetaBoostreco,weight);
        hl2cosThetaBoostreco->Fill(l2cosThetaBoostreco,weight);
        hl3cosThetaBoostreco->Fill(l3cosThetaBoostreco,weight);
        hl4cosThetaBoostreco->Fill(l4cosThetaBoostreco,weight);
        hfourlcosThetaBoostreco->Fill(fourlcosThetaBoostreco,weight);
        hl1l2CScosThetareco->Fill(l1l2CScosThetareco,weight);
        hl3l4CScosThetareco->Fill(l3l4CScosThetareco,weight);
	}
      // particle
	if(switchVal_particle==0){
	
        // pT + m
        hl1pTparticle->Fill(l1_particle.Pt(),weight);
        hl2pTparticle->Fill(l2_particle.Pt(),weight);
        hl3pTparticle->Fill(l3_particle.Pt(),weight);
        hl4pTparticle->Fill(l4_particle.Pt(),weight);
        // phi
        hl1phiparticle->Fill(l1_particle.Phi() ,weight);
        hl2phiparticle->Fill(l2_particle.Phi() ,weight);
        hl3phiparticle->Fill(l3_particle.Phi() ,weight);
        hl4phiparticle->Fill(l4_particle.Phi() ,weight);
        hl1l2deltaPhiparticle->Fill(l1l2deltaPhiparticle ,weight);
        hl3l4deltaPhiparticle->Fill(l3l4deltaPhiparticle ,weight);
        // eta
        hl1etaparticle->Fill(l1_particle.Eta() ,weight);
        hl2etaparticle->Fill(l2_particle.Eta() ,weight);
        hl3etaparticle->Fill(l3_particle.Eta() ,weight);
        hl4etaparticle->Fill(l4_particle.Eta() ,weight);
        hl1l2deltaEtaparticle->Fill(l1l2deltaEtaparticle ,weight);
        hl3l4deltaEtaparticle->Fill(l3l4deltaEtaparticle ,weight);
        // R
        hl1Rparticle->Fill(sqrt((l1_particle.Phi()*l1_particle.Phi())+(l1_particle.Eta()*l1_particle.Eta())) ,weight);
        hl2Rparticle->Fill(sqrt((l2_particle.Phi()*l2_particle.Phi())+(l2_particle.Eta()*l2_particle.Eta())) ,weight);
        hl3Rparticle->Fill(sqrt((l3_particle.Phi()*l3_particle.Phi())+(l3_particle.Eta()*l3_particle.Eta())) ,weight);
        hl4Rparticle->Fill(sqrt((l4_particle.Phi()*l4_particle.Phi())+(l4_particle.Eta()*l4_particle.Eta())) ,weight);
        hl1l2deltaRparticle->Fill(l1l2deltaRparticle ,weight);
        hl3l4deltaRparticle->Fill(l3l4deltaRparticle ,weight);
        // cos
        hl1cosThetaparticle->Fill(l1cosThetaparticle,weight);
        hl2cosThetaparticle->Fill(l2cosThetaparticle,weight);
        hl3cosThetaparticle->Fill(l3cosThetaparticle,weight);
        hl4cosThetaparticle->Fill(l4cosThetaparticle,weight);
        hfourlcosThetaparticle->Fill(fourlcosThetaparticle,weight);
        hl1cosThetaBoostparticle->Fill(l1cosThetaBoostparticle,weight);
        hl2cosThetaBoostparticle->Fill(l2cosThetaBoostparticle,weight);
        hl3cosThetaBoostparticle->Fill(l3cosThetaBoostparticle,weight);
        hl4cosThetaBoostparticle->Fill(l4cosThetaBoostparticle,weight);
        hfourlcosThetaBoostparticle->Fill(fourlcosThetaBoostparticle,weight);
        hl1l2CScosThetaparticle->Fill(l1l2CScosThetaparticle,weight);
        hl3l4CScosThetaparticle->Fill(l3l4CScosThetaparticle,weight);
	}
      // parton
	if(switchVal_parton==0 ){
        // pT + m
        hl1pTparton->Fill(l1_parton.Pt(),weight);
        hl2pTparton->Fill(l2_parton.Pt(),weight);
        hl3pTparton->Fill(l3_parton.Pt(),weight);
        hl4pTparton->Fill(l4_parton.Pt(),weight);
        // phi
        hl1phiparton->Fill(l1_parton.Phi() ,weight);
        hl2phiparton->Fill(l2_parton.Phi() ,weight);
        hl3phiparton->Fill(l3_parton.Phi() ,weight);
        hl4phiparton->Fill(l4_parton.Phi() ,weight);
        hl1l2deltaPhiparton->Fill(l1l2deltaPhiparton ,weight);
        hl3l4deltaPhiparton->Fill(l3l4deltaPhiparton ,weight);
        // eta
        hl1etaparton->Fill(l1_parton.Eta() ,weight);
        hl2etaparton->Fill(l2_parton.Eta() ,weight);
        hl3etaparton->Fill(l3_parton.Eta() ,weight);
        hl4etaparton->Fill(l4_parton.Eta() ,weight);
        hl1l2deltaEtaparton->Fill(l1l2deltaEtaparton ,weight);
        hl3l4deltaEtaparton->Fill(l3l4deltaEtaparton ,weight);
        // R
        hl1Rparton->Fill(sqrt((l1_parton.Phi()*l1_parton.Phi())+(l1_parton.Eta()*l1_parton.Eta())) ,weight);
        hl2Rparton->Fill(sqrt((l2_parton.Phi()*l2_parton.Phi())+(l2_parton.Eta()*l2_parton.Eta())) ,weight);
        hl3Rparton->Fill(sqrt((l3_parton.Phi()*l3_parton.Phi())+(l3_parton.Eta()*l3_parton.Eta())) ,weight);
        hl4Rparton->Fill(sqrt((l4_parton.Phi()*l4_parton.Phi())+(l4_parton.Eta()*l4_parton.Eta())) ,weight);
        hl1l2deltaRparton->Fill(l1l2deltaRparton ,weight);
        hl3l4deltaRparton->Fill(l3l4deltaRparton ,weight);
        // cos
        hl1cosThetaparton->Fill(l1cosThetaparton,weight);
        hl2cosThetaparton->Fill(l2cosThetaparton,weight);
        hl3cosThetaparton->Fill(l3cosThetaparton,weight);
        hl4cosThetaparton->Fill(l4cosThetaparton,weight);
        hfourlcosThetaparton->Fill(fourlcosThetaparton,weight);
        hl1cosThetaBoostparton->Fill(l1cosThetaBoostparton,weight);
        hl2cosThetaBoostparton->Fill(l2cosThetaBoostparton,weight);
        hl3cosThetaBoostparton->Fill(l3cosThetaBoostparton,weight);
        hl4cosThetaBoostparton->Fill(l4cosThetaBoostparton,weight);
        hfourlcosThetaBoostparton->Fill(fourlcosThetaBoostparton,weight);
        hl1l2CScosThetaparton->Fill(l1l2CScosThetaparton,weight);
        hl3l4CScosThetaparton->Fill(l3l4CScosThetaparton,weight);
	}

// comp
// 2D - parton(1) particle(2) reco(3)
	if(switchVal_parton==0 && switchVal_particle==0 ){
	  hHpT12Comp -> Fill(h_parton.Pt(), h_particle.Pt(), weight);
	  hHm12Comp -> Fill(h_parton.M(), h_particle.M(), weight);
	  hbbdeltaPhi12Comp -> Fill(bbdeltaPhiparton, bbdeltaPhiparticle, weight);
	  hbbdeltaEta12Comp -> Fill(bbdeltaEtaparton, bbdeltaEtaparticle, weight);
	  hjjpT12Comp -> Fill(j1_parton.Pt()+j2_parton.Pt(),j1_particle.Pt()+j2_particle.Pt(), weight);
	  hjjdeltaPhi12Comp -> Fill(jjdeltaPhiparton, jjdeltaPhiparticle, weight);
	  hz1pT23Comp->Fill(z1_parton.Pt(), z1_particle.Pt(), weight);
	  hz1m23Comp->Fill(z1_parton.M(), z1_particle.M(), weight);
	  hz2pT23Comp->Fill(z2_parton.Pt(), z2_particle.Pt(), weight);
	  hz2m23Comp->Fill(z2_parton.M(), z2_particle.M(), weight);
	}
	if(switchVal_reco ==0  && switchVal_particle==0 ){
	  hHpT23Comp -> Fill(h_particle.Pt(), h_reco.Pt(), weight);
	  hHm23Comp -> Fill(h_particle.M(), h_reco.M(), weight);
	  hbbdeltaPhi23Comp -> Fill(bbdeltaPhiparticle, bbdeltaPhireco, weight);
	  hbbdeltaEta23Comp -> Fill(bbdeltaEtaparticle, bbdeltaEtareco, weight);
	  hjjpT23Comp -> Fill(j1_particle.Pt()+j2_particle.Pt(),j1_reco.Pt()+j2_reco.Pt(), weight);
	  hjjdeltaPhi23Comp -> Fill(jjdeltaPhiparticle, jjdeltaPhireco, weight);
	  hj1Phi23Comp -> Fill(j1_particle.Phi(), j1_reco.Phi(), weight);
	  hj2Phi23Comp -> Fill(j2_particle.Phi(), j2_reco.Phi(), weight);
	  hz1pT12Comp->Fill(z1_particle.Pt(), z1_reco.Pt(), weight);
	  hz1m12Comp->Fill(z1_particle.M(), z1_reco.M(), weight);
	  hz2pT12Comp->Fill(z2_particle.Pt(), z2_reco.Pt(), weight);
	  hz2m12Comp->Fill(z2_particle.M(), z2_reco.M(), weight);
	}
	if(switchVal_reco==0  && switchVal_parton==0 ){
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
	  // l vs h
	  hl1l2deltaPhiHpTcompreco->Fill(l1l2deltaPhiparton, h_reco.Pt(), weight);
	  hl3l4deltaPhiHpTcompreco->Fill(l3l4deltaPhiparton, h_reco.Pt(), weight);
	}
	
  
	
	if(switchVal_reco==0){
	  hHz1pTcompreco->Fill(h_reco.Pt(), z1_reco.Pt(), weight);
	  hHz2pTcompreco->Fill(h_reco.Pt(), z2_reco.Pt(), weight);
	  hHpTz1etacompreco->Fill(h_reco.Pt(), z1_reco.Eta(), weight);
	  hHpTz2etacompreco->Fill(h_reco.Pt(), z2_reco.Eta(), weight);
	  hHzzpTcompreco->Fill(h_reco.Pt(), z1_reco.Pt() + z2_reco.Pt(), weight);
	}
	if(switchVal_particle==0){
	  hHz1pTcompparticle->Fill(h_particle.Pt(), z1_particle.Pt(), weight);
	  hHz2pTcompparticle->Fill(h_particle.Pt(), z2_particle.Pt(), weight);
	  hHpTz1etacompparticle->Fill(h_particle.Pt(), z1_particle.Eta(), weight);
	  hHpTz2etacompparticle->Fill(h_particle.Pt(), z2_particle.Eta(), weight);
	  hHzzpTcompparticle->Fill(h_particle.Pt(), z1_particle.Pt() + z2_particle.Pt(), weight);
	}
	if(switchVal_parton){
	  hHz1pTcompparton->Fill(h_parton.Pt(), z1_parton.Pt(), weight);
	  hHz2pTcompparton->Fill(h_parton.Pt(), z2_parton.Pt(), weight);
	  hHpTz1etacompparton->Fill(h_parton.Pt(), z1_parton.Eta(), weight);
	  hHpTz2etacompparton->Fill(h_parton.Pt(), z2_parton.Eta(), weight);
	  hHzzpTcompparton->Fill(h_parton.Pt(), z1_parton.Pt() + z2_parton.Pt(), weight);
	}

    
    
  nPassed+=weight;
  nPassedRaw++;

  cutVal_reco++; cutValW_reco+=weight;
  cutVal_particle++; cutValW_particle+=weight;
  //cutVal_parton++; cutValW_parton+=weight;
  }
  

//------------------------------------------------------------------------------------------------------------------------------------------------------------
// write histos
//------------------------------------------------------------------------------------------------------------------------------------------------------------


   cout<<" Reco CutF Flow "<<endl;
   PrintCutFlow(cutFlowMap_reco,cutList_reco,  "Reco");
   cout<<" Particle  CutF Flow "<<endl;
   PrintCutFlow(cutFlowMap_particle,cutList_particle, "Particle");
   cout<<" Parton  CutF Flow "<<endl;
   PrintCutFlow(cutFlowMap_parton,cutList_parton, "Parton");
   
   for(std::vector<string>::iterator it=selType.begin(); it!=selType.end(); it++){
     FillCutFlow(cutFlowHists[(*it)],cutFlowEffs[(*it)],*cutFlowMapAll[(*it)],cutFlowMByType[(*it)], (*it));
   }
  
  
   TFile *hists= new TFile(outputFile,"recreate");
   hists->cd();
   for(std::vector<TH1F*>::iterator h=listOfTH1.begin(); h!=listOfTH1.end(); h++){
     ( *h)->Write();
     delete (*h); 
   }
   for(std::vector<TH2F*>::iterator h=listOfTH2.begin(); h!=listOfTH2.end(); h++){
     ( *h)->Write();
     delete (*h); 
   }

   for(std::vector<TProfile*>::iterator h=listOfTPorifles.begin(); h!=listOfTPorifles.end(); h++){
     ( *h)->Write();
     delete (*h); 
   }
   
   kappaLambda -> Write();
   delete kappaLambda;
   
   hists->Close();
   
  std::cout<<"Total number of entries "<<numberOfEntries<<" Passed "<<nPassed<<" raw "<<nPassedRaw<<std::endl;


  
}

int main(int argc, char* argv[]) {
  const char *inputFileName = argv[1];
  const char *outputFileName = argv[2];
  // O: for ZZ H JJ, 1: (->H) ZZ  jj: 2: (->H) ZZ, 3: Hjj, 4: WW H JJ, 5: WW (->H) jj: 6: (->H) WW,
  int analysisType=0;
  if( argc > 2 )  analysisType=atoi(argv[3]);
  //  const char *process_name = argv[3];                                                                                                                                                                                                                                                                                                                                   
  zAnalyzer(inputFileName, outputFileName,analysisType);
  return 1;
}
