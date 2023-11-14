#ifdef __CLING__
R__LOAD_LIBRARY("libDelphes")
#endif

#include "classes/DelphesClasses.h"
#include "classes/DelphesLHEFReader.h"
#include "ExRootAnalysis/ExRootTreeReader.h"
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

using namespace std;

//------------------------------------------------------------------------------------------------------------------------------------------------------------
// initialize combinations
//------------------------------------------------------------------------------------------------------------------------------------------------------------

vector <string> types={"4mu","4e", "2mu2e", "2e2mu"};
std::vector <string> cutList_reco={"initial reco", "1 btag reco", "2 good j reco", "2 b-like jet pairs reco", "found bb reco", "2 vbfj reco", "vbfj pairs","OSFL"};
std::vector <string> cutList_particle={"initial particle", "1 btag particle", "2 good j particle", "2 b-like jet pairs part", "found bb particle", "2 vbfj particle", "comb vbf part","OSFL"};

void PrintCutFlow(std::map<string, std::pair<int,double>> cutFlowMap, std::vector <string> cutList, string label){

std::cout<<std::left<<std::setw(25)<<label<<" Cut"<<std::setw(10)<<label<<" Passed"<<std::setw(15)<<" Rel Eff "<< std::setw(15)<<label <<" Efficiency"<< std::endl;
    for(int i=0; i<(int) cutList.size(); i++) {
        const std::string cutName = cutList[i];
        double passed_reco =  cutFlowMap[cutName].first;
        double efficiency_reco = 100.00 * cutFlowMap[cutName].second / cutFlowMap[cutList[0]].second;

	double relEff= 100;;
	if(i>0)
	  relEff = 100.00 * cutFlowMap[cutName].second / cutFlowMap[cutList.at(i-1)].second; 
	std::cout<<std::left<<std::setw(25)<<cutName<<std::setw(10)<<passed_reco<< std::setw(15)<<relEff <<std::setw(15)<<efficiency_reco<< std::endl;
    }

    std::cout<<std::left<<std::setw(25)<<label<<" Cut"<<std::setw(10)<<label<<" Passed"<<std::setw(15)<< "Rel Eff "<<std::setw(15)<<label<<" Efficiency"<< std::endl;
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

void draw_hist2(TH2 *histo, const char *name, const char *title, const char *xaxistitle, const char *yaxistitle) {
  TCanvas *c = new TCanvas(name, title, 1500, 1200);
  histo->GetXaxis()->SetTitle(xaxistitle);
  histo->GetYaxis()->SetTitle(yaxistitle);
  //histo->SetStatX(0.875);
  //histo->SetStatY(0.875);
  histo->Draw("COLZ");
  PrintCanvas(c, name);
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
    return NULL;
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
void zAnalyzer(const char *inputFile,const char *outputFile) {

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
        TH2 *hHpT12Comp = new TH2F("H_pT_comp_12", "p_{T}^{hbb}", pTBins, hpTmin, hpTmax, pTBins, hpTmin, hpTmax);
        TH2 *hHm12Comp = new TH2F("H_m_comp_12", "m_{hbb}", mBins, hmmin, hmmax, mBins, hmmin, hmmax);
        TH2 *hHpT23Comp = new TH2F("H_pT_comp_23", "p_{T}^{hbb}", pTBins, hpTmin, hpTmax, pTBins, hpTmin, hpTmax);
        TH2 *hHm23Comp = new TH2F("H_m_comp_23", "m_{hbb}", mBins, hmmin, hmmax, mBins, hmmin, hmmax);
        TH2 *hHpT13Comp = new TH2F("H_pT_comp_13", "p_{T}^{hbb}", pTBins, hpTmin, hpTmax, pTBins, hpTmin, hpTmax);
        TH2 *hHm13Comp = new TH2F("H_m_comp_13", "m_{hbb}", mBins, hmmin, hmmax, mBins, hmmin, hmmax);
        TH2 *hbbdeltaPhi12Comp = new TH2F("bb_#Delta#phi_comp_12", "#Delta#phi_{bb}", phiBins, -TMath::Pi(),+TMath::Pi(), phiBins, -TMath::Pi(),+TMath::Pi());
        TH2 *hbbdeltaPhi23Comp = new TH2F("bb_#Delta#phi_comp_23", "#Delta#phi_{bb}", phiBins, -TMath::Pi(),+TMath::Pi(), phiBins, -TMath::Pi(),+TMath::Pi());
        TH2 *hbbdeltaPhi13Comp = new TH2F("bb_#Delta#phi_comp_13", "#Delta#phi_{bb}", phiBins, -TMath::Pi(),+TMath::Pi(), phiBins, -TMath::Pi(),+TMath::Pi());
        TH2 *hbbdeltaEta12Comp = new TH2F("bb_#Delta#eta_comp_12", "#Delta#eta_{bb}", etaBins, hetamin, hetamax, etaBins, hetamin, hetamax);
        TH2 *hbbdeltaEta23Comp = new TH2F("bb_#Delta#eta_comp_23", "#Delta#eta_{bb}", etaBins, hetamin, hetamax, etaBins, hetamin, hetamax);
        TH2 *hbbdeltaEta13Comp = new TH2F("bb_#Delta#eta_comp_13", "#Delta#eta_{bb}", etaBins, hetamin, hetamax, etaBins, hetamin, hetamax);
    // jets
        TH2 *hjjpT12Comp = new TH2F("jj_pT_comp_12", "p_{T}^{jj}", pTBins, jpTmin, jpTmax, pTBins, jpTmin, jpTmax);
        TH2 *hjjpT23Comp = new TH2F("jj_pT_comp_23", "p_{T}^{jj}", pTBins, jpTmin, jpTmax, pTBins, jpTmin, jpTmax);
        TH2 *hjjpT13Comp = new TH2F("jj_pT_comp_13", "p_{T}^{jj}", pTBins, jpTmin, jpTmax, pTBins, jpTmin, jpTmax);
        TH2 *hj1Phi23Comp = new TH2F("j1_#phi_comp_23", "#phi_{j1}", phiBins, -TMath::Pi(),+TMath::Pi(), phiBins, -TMath::Pi(),+TMath::Pi());
        TH2 *hj2Phi23Comp = new TH2F("j2_#phi_comp_23", "#phi_{j2}", phiBins, -TMath::Pi(),+TMath::Pi(), phiBins, -TMath::Pi(),+TMath::Pi());
        TH2 *hjjdeltaPhi12Comp = new TH2F("jj_#Delta#phi_comp_12", "#Delta#phi_{jj}", phiBins, -TMath::Pi(),+TMath::Pi(), phiBins, -TMath::Pi(),+TMath::Pi());
        TH2 *hjjdeltaPhi23Comp = new TH2F("jj_#Delta#phi_comp_23", "#Delta#phi_{jj}", phiBins, -TMath::Pi(),+TMath::Pi(), phiBins, -TMath::Pi(),+TMath::Pi());        
        TH2 *hjjdeltaPhi13Comp = new TH2F("jj_#Delta#phi_comp_13", "#Delta#phi_{jj}", phiBins, -TMath::Pi(),+TMath::Pi(), phiBins, -TMath::Pi(),+TMath::Pi());
    //z
      // pT + m
        TH2 *hz1pT12Comp = new TH2F("z1_pT_comp_12", "p^{T}_{z1}", pTBins, zpTmin, zpTmax, pTBins, zpTmin, zpTmax);
        TH2 *hz1m12Comp = new TH2F("z1_m_comp_12", "m_{z1}", mBins, zmmin, zmmax, mBins, zmmin, zmmax);
        TH2 *hz2pT12Comp = new TH2F("z2_pT_comp_12", "p^{T}_{z2}", pTBins, zpTmin, zpTmax, pTBins, zpTmin, zpTmax);
        TH2 *hz2m12Comp = new TH2F("z2_m_comp_12", "m_{z2}", mBins, zmmin, zmmax, mBins, zmmin, zmmax);
        TH2 *hz1pT23Comp = new TH2F("z1_pT_comp_23", "p^{T}_{z1}", pTBins, zpTmin, zpTmax, pTBins, zpTmin, zpTmax);
        TH2 *hz1m23Comp = new TH2F("z1_m_comp_23", "m_{z1}", mBins, zmmin, zmmax, mBins, zmmin, zmmax);
        TH2 *hz2pT23Comp = new TH2F("z2_pT_comp_23", "p^{T}_{z2}", pTBins, zpTmin, zpTmax, pTBins, zpTmin, zpTmax);
        TH2 *hz2m23Comp = new TH2F("z2_m_comp_23", "m_{z2}", mBins, zmmin, zmmax, mBins, zmmin, zmmax);
        TH2 *hz1pT13Comp = new TH2F("z1_pT_comp_13", "p^{T}_{z1}", pTBins, zpTmin, zpTmax, pTBins, zpTmin, zpTmax);
        TH2 *hz1m13Comp = new TH2F("z1_m_comp_13", "m_{z1}", mBins, zmmin, zmmax, mBins, zmmin, zmmax);
        TH2 *hz2pT13Comp = new TH2F("z2_pT_comp_13", "p^{T}_{z2}", pTBins, zpTmin, zpTmax, pTBins, zpTmin, zpTmax);
        TH2 *hz2m13Comp = new TH2F("z2_m_comp_13", "m_{z2}", mBins, zmmin, zmmax, mBins, zmmin, zmmax);

        TH2 *hl1l2deltaPhiHpTcompreco = new TH2F("l1l2_delta#phi_H_pT_comp_reco", "l1l2_delta#phi_H_pT_comp_reco", phiBins, -TMath::Pi(),+TMath::Pi(), pTBins, hpTmin, hpTmax);
        TH2 *hl3l4deltaPhiHpTcompreco = new TH2F("l3l4_delta#phi_H_pT_comp_reco", "l1l2_delta#phi_H_pT_comp_reco", phiBins, -TMath::Pi(),+TMath::Pi(), pTBins, hpTmin, hpTmax);

        TH2 *hHz1pTcompreco = new TH2F("h_z1_pT_comp_reco", "h_z1_pT_comp_reco", pTBins, hpTmin, hpTmax, pTBins, zpTmin, zpTmax);
        TH2 *hHz2pTcompreco = new TH2F("h_z2_pT_comp_reco", "h_z2_pT_comp_reco", pTBins, hpTmin, hpTmax, pTBins, zpTmin, zpTmax);
        TH2 *hHz1pTcompparticle = new TH2F("h_z1_pT_comp_particle", "h_z1_pT_comp_particle", pTBins, hpTmin, hpTmax, pTBins, zpTmin, zpTmax);
        TH2 *hHz2pTcompparticle = new TH2F("h_z2_pT_comp_particle", "h_z2_pT_comp_particle", pTBins, hpTmin, hpTmax, pTBins, zpTmin, zpTmax);
        TH2 *hHz1pTcompparton = new TH2F("h_z1_pT_comp_parton", "h_z1_pT_comp_parton", pTBins, hpTmin, hpTmax, pTBins, zpTmin, zpTmax);
        TH2 *hHz2pTcompparton = new TH2F("h_z2_pT_comp_parton", "h_z2_pT_comp_parton", pTBins, hpTmin, hpTmax, pTBins, zpTmin, zpTmax);

        TH2 *hHpTz1etacompreco = new TH2F("h_pT_z1_eta_comp_reco", "h_pT_z1_eta_comp_reco", pTBins, hpTmin, hpTmax, etaBins, zetamin, zetamax);
        TH2 *hHpTz2etacompreco = new TH2F("h_pT_z2_eta_comp_reco", "h_pT_z2_eta_comp_reco", pTBins, hpTmin, hpTmax, etaBins, zetamin, zetamax);
        TH2 *hHpTz1etacompparticle = new TH2F("h_pT_z1_eta_comp_particle", "h_pT_z1_eta_comp_particle", pTBins, hpTmin, hpTmax, etaBins, zetamin, zetamax);
        TH2 *hHpTz2etacompparticle = new TH2F("h_pT_z2_eta_comp_particle", "h_pT_z2_eta_comp_particle", pTBins, hpTmin, hpTmax, etaBins, zetamin, zetamax);
        TH2 *hHpTz1etacompparton = new TH2F("h_pT_z1_eta_comp_parton", "h_pT_z1_eta_comp_parton", pTBins, hpTmin, hpTmax, etaBins, zetamin, zetamax);
        TH2 *hHpTz2etacompparton= new TH2F("h_pT_z2_eta_comp_parton", "h_pT_z2_eta_comp_parton", pTBins, hpTmin, hpTmax, etaBins, zetamin, zetamax);

        TH2 *hHzzpTcompreco = new TH2F("h_zz_pT_comp_reco", "h_zz_pT_comp_reco", pTBins, hpTmin, 2*hpTmax, pTBins, zpTmin, 2*zpTmax);
        TH2 *hHzzpTcompparticle = new TH2F("h_zz_pT_comp_particle", "h_zz_pT_comp_particle", pTBins, hpTmin, 2*hpTmax, pTBins, zpTmin, 2*zpTmax);
        TH2 *hHzzpTcompparton = new TH2F("h_zz_pT_comp_parton", "h_zz_pT_comp_parton", pTBins, hpTmin, 2*hpTmax, pTBins, zpTmin, 2*zpTmax);

  TProfile *kappaLambda = new TProfile("kappaLambda", "kappaLambda", 40, -20, 20);
  kappaLambda -> GetXaxis() -> SetTitle("#kappa_{#lambda}");

  double  nPassed=0;
  double Lumi=3e3;
  double totWeightedEntries=0;
  int  nPassedRaw=0;
  int nQuads=0;

  GenParticle *daughter1;
  GenParticle *daughter2;

  GenParticle *particle1;
  GenParticle *particle2;

  GenParticle *p1daughter1;
  GenParticle *p1daughter2;
  GenParticle *p2daughter1;
  GenParticle *p2daughter2;

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
std::vector <string> cutList_parton={};
std::map<string, std::pair<int,double>> cutFlowMap_parton;
for(int i=0; i<(int) cutList_parton.size(); i++) { 
  cutFlowMap_parton[cutList_parton.at(i)] = make_pair(0,0.0); 
}

//------------------------------------------------------------------------------------------------------------------------------------------------------------
// loop
//------------------------------------------------------------------------------------------------------------------------------------------------------------
  for(Int_t entry = 0; entry < numberOfEntries; ++entry) {

    treeReader->ReadEntry(entry);
    std::map<int,double> kappaLambdaWeights;
    HepMCEvent *event = (HepMCEvent*) branchEvent -> At(0);

    Float_t weight = event->Weight/numberOfEntries*Lumi;
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
      if( jet->BTag>0) btagIndex.push_back(i);
      else noBtag.push_back(i);
           goodJetIndex.push_back(i);
    }

    increaseCount(cutFlowMap_reco,"initial reco",weight);
    //cutFlowMap_reco["initial reco"] = {cutVal_reco,cutValW_reco};

   sort(btagIndex.begin(), btagIndex.end(), [branchJet](const int& lhs, const int& rhs) {
       return ((Jet*)branchJet->At(lhs))->PT < ((Jet*)branchJet->At(rhs))->PT;
     });
   sort(noBtag.begin(), noBtag.end(), [branchJet](const int& lhs, const int& rhs) {
       return ((Jet*)branchJet->At(lhs))->PT < ((Jet*)branchJet->At(rhs))->PT;
     });
   sort(goodJetIndex.begin(), goodJetIndex.end(), [branchJet](const int& lhs, const int& rhs) {
       return ((Jet*)branchJet->At(lhs))->PT < ((Jet*)branchJet->At(rhs))->PT;
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

   if(switchVal_reco == 0 )
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
     if( b1->BTag>0 && b2->BTag>0) {
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
    
  //        sort(nonHiggsJet.begin(), nonHiggsJet.end(), [branchJet](const int& lhs, const int& rhs) {
  //	    return ((Jet*)branchJet->At(lhs))->PT < ((Jet*)branchJet->At(rhs))->PT;
  //    });

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

  if(switchVal_reco==0)    
    for(int i=0; i<(int)vbfJetIndexComb.size(); i++)
    vbfJetIndex.push_back(make_pair(vbfJetIndexComb[i][0],vbfJetIndexComb[i][1]));
  
  // GB: What is this ? 
  //  if( branchMuon->GetEntries() + branchElectron->GetEntries() < 4) continue;
  
  if( switchVal_reco==0 && vbfJetIndex.size() > 1) 
    sort(vbfJetIndex.begin(), vbfJetIndex.end(), [branchJet](const pair<int,int> lhs, const pair<int,int> rhs) {
	return fabs((((Jet*)branchJet->At(lhs.first))->Eta - ((Jet*)branchJet->At(lhs.second))->Eta) ) >
	  fabs((((Jet*)branchJet->At(rhs.first))->Eta - ((Jet*)branchJet->At(rhs.second))->Eta) ) ; 
      });
  
  /*
    int vbfJetsIndexCandidate = -1;
    
    // loop and take first w eta > 2.5
    for (int i=0; i<(int)vbfJetIndex.size(); i++) {
    if( fabs((((Jet*)branchGenJet->At(vbfJetIndex[i].first))->Eta - ((Jet*)branchGenJet->At(vbfJetIndex[i].second))->Eta)) > 2.5 ) {
    vbfJetsIndexCandidate = i;
    break;
    }
    }
    
    cutFlowMap_reco["2.5 deltaEta vbf reco"] = {cutVal_reco,cutValW_reco}; //last reco cut
    
    Jet *jet1 = (Jet*) branchJet->At(vbfJetIndex[vbfJetsIndexCandidate].first);
    Jet *jet2 = (Jet*) branchJet->At(vbfJetIndex[vbfJetsIndexCandidate].second);
  */


  // GB, well what happens here if we don't have two jets in the event?
  //Jet *jet1 = (Jet*) branchJet->At(vbfJetIndex[0].first);
  //Jet *jet2 = (Jet*) branchJet->At(vbfJetIndex[0].second);

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
      if( el_reco->PT > 1 && fabs(el_reco->Eta) < 2.5) goodE_reco_indices.push_back(i);
      
    }

    // sort the indices by pT ;
    sort(goodE_reco_indices.begin(), goodE_reco_indices.end(), [branchElectron](const int& lhs, const int& rhs) {
	return ((Electron*)branchElectron->At(lhs))->PT < ((Electron*)branchElectron->At(rhs))->PT;
      });

    vector <int> goodMu_reco_indices; 
    for(int i=0; i<(int)branchMuon->GetEntries(); i++){
      Muon *mu_reco = (Muon *) branchMuon->At(i);
    
    // pT and eta cuts 
      if( mu_reco->PT > 1 && fabs(mu_reco->Eta) < 2.5) goodMu_reco_indices.push_back(i);

    }

    // sort the indices by pT ;
    sort(goodMu_reco_indices.begin(), goodMu_reco_indices.end(), [branchMuon](const int& lhs, const int& rhs) {
	return ((Muon*)branchMuon->At(lhs))->PT < ((Muon*)branchMuon->At(rhs))->PT;
      });
    
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

    sort(muRecoPairIndices.begin(),muRecoPairIndices.end(), [branchMuon]( pair<int,int>   & lhs,  pair<int,int>   & rhs) {
	    int index11_reco=(lhs).first;
	    int index12_reco=(lhs).second;
	    int index21_reco=(rhs).first;
    	int index22_reco=(rhs).second;
	    return fabs(((((Muon*)branchMuon->At(index11_reco))->P4() + ((Muon*)branchMuon->At(index12_reco))->P4())).M() - 91) <
	  fabs( ((((Muon*)branchMuon->At(index21_reco))->P4() + ((Muon*)branchMuon->At(index22_reco))->P4()).M()) -91);
      });

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
	
	return fabs(mass1_reco ) <  fabs(mass2_reco) ;
	
      });


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
      if( jet->BTag>0) btagIndexParticle.push_back(i);
      else noBtagParticle.push_back(i);
           goodJetIndexParticle.push_back(i);
    }

    increaseCount(cutFlowMap_particle,"initial particle",weight);
    //cutFlowMap_particle["initial particle"] = {cutVal_particle,cutValW_particle};

    // at least one b tag 
    if(switchVal_particle == 0 && btagIndex.size() < 1) increaseCount(cutFlowMap_particle,"1 btag particle",weight);
    else  switchVal_particle = 1;

    // at least two jets
    if(switchVal_particle == 0 && goodJetIndexParticle.size() > 1) increaseCount(cutFlowMap_particle,"2 good j particle",weight);
    else switchVal_particle = 1;
     
    
    
    sort(btagIndexParticle.begin(), btagIndexParticle.end(), [branchGenJet](const int& lhs, const int& rhs) {
	return ((Jet*)branchGenJet->At(lhs))->PT < ((Jet*)branchGenJet->At(rhs))->PT;
      });
    sort(noBtagParticle.begin(), noBtagParticle.end(), [branchGenJet](const int& lhs, const int& rhs) {
	return ((Jet*)branchGenJet->At(lhs))->PT < ((Jet*)branchGenJet->At(rhs))->PT;
      });
    sort(goodJetIndexParticle.begin(), goodJetIndexParticle.end(), [branchGenJet](const int& lhs, const int& rhs) {
	return ((Jet*)branchGenJet->At(lhs))->PT < ((Jet*)branchGenJet->At(rhs))->PT;
      });
    
    Jet *b1Particle=nullptr;
    Jet *b2Particle=nullptr;
    
    vector<pair<int,int>> bJetPairsParticle;
    vector<vector <int>> bJetPairsCombParticle=combinationsNoRepetitionAndOrderDoesNotMatter(2,goodJetIndexParticle);
//  vector<vector <int>> bJetPairsCombParticle=combinationsNoRepetitionAndOrderDoesNotMatter(2,btagIndexParticle);
    
    //if( bJetPairsCombParticle.size() < 1) continue; // need at least two good jets;  
    if(switchVal_particle==0 && bJetPairsCombParticle.size() >0 ) 
      increaseCount(cutFlowMap_particle,"2 b-like jet pairs part",weight);
    else
      switchVal_particle=1; 

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

    // b pair 
    if(switchVal_particle == 0 && !foundBjetParticle)  
      increaseCount(cutFlowMap_particle,"found bb particle",weight);
    else switchVal_particle = 1;
    
    // jets

    vector <int> nonHiggsJetParticle;

    for(int i=0; i<(int)goodJetIndexParticle.size(); i++){
      if( goodJetIndexParticle[i] == higgsbbcandidateParticle.first  || goodJetIndexParticle[i] == higgsbbcandidateParticle.second) continue;
      nonHiggsJetParticle.push_back(i);
    }

    //        sort(nonHiggsJet.begin(), nonHiggsJet.end(), [branchGenJet](const int& lhs, const int& rhs) {
    //	    return ((Jet*)branchGenJet->At(lhs))->PT < ((Jet*)branchGenJet->At(rhs))->PT;
    //    });

    // at least 2 vbf jets 
    if(switchVal_particle == 0 && nonHiggsJetParticle.size() < 2) 
      increaseCount(cutFlowMap_particle,"2 vbfj particle",weight);
    else switchVal_particle = 1;
    

    vector<pair<int,int>> vbfJetIndexParticle;
    vector<vector <int>> vbfJetIndexCombParticle;
    if(switchVal_particle==0)
      combinationsNoRepetitionAndOrderDoesNotMatter(2,nonHiggsJetParticle);
    
    if( vbfJetIndexCombParticle.size() < 1 ) switchVal_particle=1; 
    else increaseCount(cutFlowMap_particle,"comb vbf part",weight);
     
 
    for(int i=0; i<(int)vbfJetIndexCombParticle.size(); i++)
      vbfJetIndexParticle.push_back(make_pair(vbfJetIndexCombParticle[i][0],vbfJetIndexCombParticle[i][1]));
    //  if( branchMuon->GetEntries() + branchElectron->GetEntries() < 4) continue;

    if( vbfJetIndexParticle.size() > 1) 
      sort(vbfJetIndexParticle.begin(), vbfJetIndexParticle.end(), [branchGenJet](const pair<int,int> lhs, const pair<int,int> rhs) {
	  return fabs((((Jet*)branchGenJet->At(lhs.first))->Eta - ((Jet*)branchGenJet->At(lhs.second))->Eta) ) >
	    fabs((((Jet*)branchGenJet->At(rhs.first))->Eta - ((Jet*)branchGenJet->At(rhs.second))->Eta) ) ; 
	});

    /*
      int vbfJetsIndexCandidate = -1;
  
      // loop and take first w eta > 2.5
      for (int i=0; i<(int)vbfJetIndex.size(); i++) {
      if( fabs((((Jet*)branchGenJet->At(vbfJetIndex[i].first))->Eta - ((Jet*)branchGenJet->At(vbfJetIndex[i].second))->Eta)) > 2.5 ) {
      vbfJetsIndexCandidate = i;
      break;
      }
      }

  
 cutFlowMap_particle["2.5 deltaEta vbf particle"] = {cutVal_particle,cutValW_particle}; //last particle cut

    Jet *jet1 = (Jet*) branchGenJet->At(vbfJetIndex[vbfJetsIndexCandidate].first);
    Jet *jet2 = (Jet*) branchGenJet->At(vbfJetIndex[vbfJetsIndexCandidate].second);
*/

    
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
    if( fabs(particle->PID) == 11 && fabs(particle->Eta)< 2.5 && particle->PT > 2.5) 
            goodE_particle_indices.push_back(i); 
}
 
vector <int> goodMu_particle_indices; 
for(int i=0; i<(int)branchGenParticle->GetEntries(); i++){
    GenParticle *particle=(GenParticle*) branchGenParticle->At(i);  
    if( particle->Status !=1 ) continue; 
    // electrons 
    if( fabs(particle->PID) == 13 && fabs(particle->Eta)< 2.5 && particle->PT > 2.5) 
            goodMu_particle_indices.push_back(i); 
}

    // sort the indices by pT ;
    sort(goodE_particle_indices.begin(), goodE_particle_indices.end(), [branchGenParticle](const int& lhs, const int& rhs) {
	return ((GenParticle*) branchGenParticle->At(lhs))->PT < ((GenParticle*) branchGenParticle->At(rhs))->PT;
      });

    // sort the indices by pT ;
    sort(goodMu_particle_indices.begin(), goodMu_particle_indices.end(), [branchGenParticle](const int& lhs, const int& rhs) {
		return ((GenParticle*) branchGenParticle->At(lhs))->PT < ((GenParticle*) branchGenParticle->At(rhs))->PT;
      });
    
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

    sort(muParticlePairIndices.begin(),muParticlePairIndices.end(), [branchGenParticle]( pair<int,int>   & lhs,  pair<int,int>   & rhs) {
	    int index11_particle=(lhs).first;
	    int index12_particle=(lhs).second;
	    int index21_particle=(rhs).first;
    	int index22_particle=(rhs).second;
	    return fabs(((((GenParticle*) branchGenParticle->At(index11_particle))->P4() + ((GenParticle*) branchGenParticle->At(index12_particle))->P4())).M() - 91) <
	  fabs( ((((GenParticle*) branchGenParticle->At(index21_particle))->P4() + ((GenParticle*) branchGenParticle->At(index22_particle))->P4()).M()) -91);
      });

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
          
	    return fabs(mass1_particle ) <  fabs(mass2_particle) ;
    });


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


//------------------------------------------------------------------------------------------------------------------------------------------------------------
// parton
//------------------------------------------------------------------------------------------------------------------------------------------------------------

// higgs


//cutFlowMap_parton["initial parton"] = {cutVal_parton,cutValW_parton};

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
        // check for b parton children
        if (abs(d1_pid) == 5 && abs(d2_pid) == 5) {
          if (daughter1 -> PT > daughter2 -> PT) {
            b1_parton = daughter1 -> P4();
            b2_parton = daughter2 -> P4();
          } else {
            b1_parton = daughter2 -> P4();
            b2_parton = daughter1 -> P4();
          }
        }
      } else if ((particle->PID == 1 || particle->PID == 2 || particle->PID == 3 || particle->PID == 4 || particle->PID == 6) && d1_pid != particle->PID && d2_pid != particle->PID){
          if (daughter1 -> PT > daughter2 -> PT) {
            j1_parton = daughter1 -> P4();
            j2_parton = daughter2 -> P4();
          } else {
            j1_parton = daughter2 -> P4();
            j2_parton = daughter1 -> P4();
          }
        }
    }

    double bbdeltaPhiparton =  9999;
    double bbdeltaEtaparton = 9999;

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
    double bbdeltaRparton=sqrt((bbdeltaPhiparton*bbdeltaPhiparton)+(bbdeltaEtaparton*bbdeltaEtaparton));
  
    double jjdeltaPhiparton = 9999;
    double jjdeltaEtaparton = 9999;

	  if (j1_parton.Eta() > j2_parton.Eta()) {
	    jjdeltaPhiparton = remainder( j1_parton.Phi() - j2_parton.Phi(), 2*M_PI );
      jjdeltaEtaparton= j1_parton.Eta() - j2_parton.Eta();
	  }
	  else{
	    jjdeltaPhiparton = remainder( j2_parton.Phi() - j1_parton.Phi(), 2*M_PI );
      jjdeltaEtaparton= j2_parton.Eta() - j1_parton.Eta();
	  }

// z + leptons

vector <int> genZBosons; 
for(int i=0; i<(int)branchGenParticle->GetEntries(); i++){
    GenParticle *particle=(GenParticle*)branchGenParticle->At(i);
    int d1_pid = 9999;
    if (particle->D1 != -1) {
        daughter1 = (GenParticle*) branchGenParticle->At(particle->D1);
        d1_pid = daughter1 -> PID;
      }
    // find a Z.
        if( particle->PID == 23 && d1_pid != 23) 
        genZBosons.push_back(i);
}

    sort(genZBosons.begin(), genZBosons.end(), [branchGenParticle](const int& lhs, const int& rhs) {
	return ((GenParticle*)branchGenParticle->At(lhs))->P4().M() < (((GenParticle*)branchGenParticle->At(rhs))->P4()).M(); 
    });

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

int q1_parton = ((GenParticle*) branchGenParticle->At(Z1children[0]))->Charge;
int q2_parton = ((GenParticle*) branchGenParticle->At(Z1children[1]))->Charge;
int q3_parton = ((GenParticle*) branchGenParticle->At(Z2children[0]))->Charge;
int q4_parton = ((GenParticle*) branchGenParticle->At(Z2children[1]))->Charge;

    //z1_parton=l1_parton + l2_parton;
    //z2_parton=l3_parton + l4_parton;
    fourl_parton=l1_parton + l2_parton + l3_parton + l4_parton;

    double l1l2deltaPhiparton=(l1_parton.Phi() > l2_parton.Phi() ? -1:+1)*TMath::Abs(l2_parton.Phi() - l1_parton.Phi());
    double l3l4deltaPhiparton=(l3_parton.Phi() > l4_parton.Phi() ? -1:+1)*TMath::Abs(l4_parton.Phi() - l3_parton.Phi());

    double l1l2deltaEtaparton=(l1_parton.Eta() > l2_parton.Eta() ? -1:+1)*TMath::Abs(l2_parton.Eta() - l1_parton.Eta());
    double l3l4deltaEtaparton=(l3_parton.Eta() > l4_parton.Eta() ? -1:+1)*TMath::Abs(l4_parton.Eta() - l3_parton.Eta());

    double l1l2deltaRparton=sqrt((l1l2deltaPhiparton*l1l2deltaPhiparton)+(l1l2deltaEtaparton*l1l2deltaEtaparton));
    double l3l4deltaRparton=sqrt((l3l4deltaPhiparton*l3l4deltaPhiparton)+(l3l4deltaEtaparton*l3l4deltaEtaparton));

    double l1cosThetaparton=l1_parton.CosTheta();
    double l2cosThetaparton=l2_parton.CosTheta();
    double l3cosThetaparton=l3_parton.CosTheta();
    double l4cosThetaparton=l4_parton.CosTheta();
    double fourlcosThetaparton=fourl_parton.CosTheta();

    l1_parton.Boost(-z1_parton.BoostVector());
    l2_parton.Boost(-z1_parton.BoostVector());
    l3_parton.Boost(-z2_parton.BoostVector());
    l4_parton.Boost(-z2_parton.BoostVector());
    double l1cosThetaBoostparton=l1_parton.CosTheta();
    double l2cosThetaBoostparton=l2_parton.CosTheta();
    double l3cosThetaBoostparton=l3_parton.CosTheta();
    double l4cosThetaBoostparton=l4_parton.CosTheta();
    double fourlcosThetaBoostparton=fourl_parton.CosTheta();
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

    // Print cutflow intermideiate

    if( entry % 1000 == 0 ){
      cout<<"Processed "<<entry<< " / " <<numberOfEntries <<" "<< entry/ numberOfEntries *100 <<" %"<<endl;
      cout<<" Reco CutF Flow "<<endl;
      PrintCutFlow(cutFlowMap_reco,cutList_reco,"Reco");
      cout<<" Particle  CutF Flow "<<endl;
      PrintCutFlow(cutFlowMap_particle,cutList_particle, "Particle");
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

// comp
// 2D - parton(1) particle(2) reco(3)

  // higgs + b       
    // pT + m
	
    hHpT12Comp -> Fill(h_parton.Pt(), h_particle.Pt(), weight);
    hHm12Comp -> Fill(h_parton.M(), h_particle.M(), weight);
    hHpT23Comp -> Fill(h_particle.Pt(), h_reco.Pt(), weight);
    hHm23Comp -> Fill(h_particle.M(), h_reco.M(), weight);
    hHpT13Comp -> Fill(h_parton.Pt(), h_reco.Pt(), weight);
    hHm13Comp -> Fill(h_parton.M(), h_reco.M(), weight);
    hbbdeltaPhi12Comp -> Fill(bbdeltaPhiparton, bbdeltaPhiparticle, weight);
    hbbdeltaPhi23Comp -> Fill(bbdeltaPhiparticle, bbdeltaPhireco, weight);
    hbbdeltaPhi13Comp -> Fill(bbdeltaPhiparton, bbdeltaPhireco, weight);
    hbbdeltaEta12Comp -> Fill(bbdeltaEtaparton, bbdeltaEtaparticle, weight);
    hbbdeltaEta23Comp -> Fill(bbdeltaEtaparticle, bbdeltaEtareco, weight);
    hbbdeltaEta13Comp -> Fill(bbdeltaEtaparton, bbdeltaEtareco, weight);

  // jets
    hjjpT12Comp -> Fill(j1_parton.Pt()+j2_parton.Pt(),j1_particle.Pt()+j2_particle.Pt(), weight);
    hjjpT23Comp -> Fill(j1_particle.Pt()+j2_particle.Pt(),j1_reco.Pt()+j2_reco.Pt(), weight);
    hjjpT13Comp -> Fill(j1_parton.Pt()+j2_parton.Pt(),j1_reco.Pt()+j2_reco.Pt(), weight);

    hjjdeltaPhi12Comp -> Fill(jjdeltaPhiparton, jjdeltaPhiparticle, weight);
    hjjdeltaPhi23Comp -> Fill(jjdeltaPhiparticle, jjdeltaPhireco, weight);
    hjjdeltaPhi13Comp -> Fill(jjdeltaPhiparton, jjdeltaPhireco, weight);

    hj1Phi23Comp -> Fill(j1_particle.Phi(), j1_reco.Phi(), weight);
    hj2Phi23Comp -> Fill(j2_particle.Phi(), j2_reco.Phi(), weight);

  // z 
    hz1pT13Comp->Fill(z1_parton.Pt(), z1_reco.Pt(), weight);
    hz1m13Comp->Fill(z1_parton.M(), z1_reco.M(), weight);
    hz2pT13Comp->Fill(z2_parton.Pt(), z2_reco.Pt(), weight);
    hz2m13Comp->Fill(z2_parton.M(), z2_reco.M(), weight);

    hz1pT12Comp->Fill(z1_particle.Pt(), z1_reco.Pt(), weight);
    hz1m12Comp->Fill(z1_particle.M(), z1_reco.M(), weight);
    hz2pT12Comp->Fill(z2_particle.Pt(), z2_reco.Pt(), weight);
    hz2m12Comp->Fill(z2_particle.M(), z2_reco.M(), weight);

    hz1pT23Comp->Fill(z1_parton.Pt(), z1_particle.Pt(), weight);
    hz1m23Comp->Fill(z1_parton.M(), z1_particle.M(), weight);
    hz2pT23Comp->Fill(z2_parton.Pt(), z2_particle.Pt(), weight);
    hz2m23Comp->Fill(z2_parton.M(), z2_particle.M(), weight);

  // l vs h
    hl1l2deltaPhiHpTcompreco->Fill(l1l2deltaPhiparton, h_reco.Pt(), weight);
    hl3l4deltaPhiHpTcompreco->Fill(l3l4deltaPhiparton, h_reco.Pt(), weight);

    hHz1pTcompreco->Fill(h_reco.Pt(), z1_reco.Pt(), weight);
    hHz2pTcompreco->Fill(h_reco.Pt(), z2_reco.Pt(), weight);
    hHz1pTcompparticle->Fill(h_particle.Pt(), z1_particle.Pt(), weight);
    hHz2pTcompparticle->Fill(h_particle.Pt(), z2_particle.Pt(), weight);
    hHz1pTcompparton->Fill(h_parton.Pt(), z1_parton.Pt(), weight);
    hHz2pTcompparton->Fill(h_parton.Pt(), z2_parton.Pt(), weight);

    hHpTz1etacompreco->Fill(h_reco.Pt(), z1_reco.Eta(), weight);
    hHpTz2etacompreco->Fill(h_reco.Pt(), z2_reco.Eta(), weight);
    hHpTz1etacompparticle->Fill(h_particle.Pt(), z1_particle.Eta(), weight);
    hHpTz2etacompparticle->Fill(h_particle.Pt(), z2_particle.Eta(), weight);
    hHpTz1etacompparton->Fill(h_parton.Pt(), z1_parton.Eta(), weight);
    hHpTz2etacompparton->Fill(h_parton.Pt(), z2_parton.Eta(), weight);

    hHzzpTcompreco->Fill(h_reco.Pt(), z1_reco.Pt() + z2_reco.Pt(), weight);
    hHzzpTcompparticle->Fill(h_particle.Pt(), z1_particle.Pt() + z2_particle.Pt(), weight);
    hHzzpTcompparton->Fill(h_parton.Pt(), z1_parton.Pt() + z2_parton.Pt(), weight);
    
  nPassed+=weight;
  nPassedRaw++;

  cutVal_reco++; cutValW_reco+=weight;
  cutVal_particle++; cutValW_particle+=weight;
  //cutVal_parton++; cutValW_parton+=weight;
  }
  

//------------------------------------------------------------------------------------------------------------------------------------------------------------
// write histos
//------------------------------------------------------------------------------------------------------------------------------------------------------------
   
   TFile *hists= new TFile(outputFile,"recreate");
   hists->cd();
   for(std::vector<TH1F*>::iterator h=listOfTH1.begin(); h!=listOfTH1.end(); h++){
     ( *h)->Write();
     delete (*h); 
   }
   /*
// higgs

    hHpTreco -> Write();
    hHmreco -> Write();
    hb1pTreco -> Write();
    hb1mreco -> Write();
    hb2pTreco -> Write();
    hb2mreco -> Write();
    hHphireco -> Write();
    hb1phireco -> Write();
    hb2phireco -> Write();
    hbbdeltaPhireco -> Write();
    hHetareco -> Write();
    hb1etareco -> Write();
    hb2etareco -> Write();
    hbbdeltaEtareco -> Write();
    hHRreco -> Write();
    hb1Rreco -> Write();
    hb2Rreco -> Write();
    hbbdeltaRreco -> Write();

    hHpTparticle -> Write();
    hHmparticle -> Write();
    hb1pTparticle -> Write();
    hb1mparticle -> Write();
    hb2pTparticle -> Write();
    hb2mparticle -> Write();
    hHphiparticle -> Write();
    hb1phiparticle -> Write();
    hb2phiparticle -> Write();
    hbbdeltaPhiparticle -> Write();
    hHetaparticle -> Write();
    hb1etaparticle -> Write();
    hb2etaparticle -> Write();
    hbbdeltaEtaparticle -> Write();
    hHRparticle -> Write();
    hb1Rparticle -> Write();
    hb2Rparticle -> Write();
    hbbdeltaRparticle -> Write();

    hHpTparton -> Write();
    hHmparton -> Write();
    hb1pTparton -> Write();
    hb1mparton -> Write();
    hb2pTparton -> Write();
    hb2mparton -> Write();
    hHphiparton -> Write();
    hb1phiparton -> Write();
    hb2phiparton -> Write();
    hbbdeltaPhiparton -> Write();
    hHetaparton -> Write();
    hb1etaparton -> Write();
    hb2etaparton -> Write();
    hbbdeltaEtaparton -> Write();
    hHRparton -> Write();
    hb1Rparton -> Write();
    hb2Rparton -> Write();
    hbbdeltaRparton -> Write();

// jets

    hjjpTreco -> Write();
    hj1pTreco -> Write();
    hj2pTreco -> Write();
    hj1phireco -> Write();
    hj2phireco -> Write();
    hjjdeltaPhireco -> Write();
    hj1etareco -> Write();
    hj2etareco -> Write();
    hjjdeltaEtareco -> Write();
    hj1Rreco -> Write();
    hj2Rreco -> Write();
    hjjdeltaRreco -> Write();

    hjjpTparticle -> Write();
    hj1pTparticle -> Write();
    hj2pTparticle -> Write();
    hj1phiparticle -> Write();
    hj2phiparticle -> Write();
    hjjdeltaPhiparticle -> Write();
    hj1etaparticle -> Write();
    hj2etaparticle -> Write();
    hjjdeltaEtaparticle -> Write();
    hj1Rparticle -> Write();
    hj2Rparticle -> Write();
    hjjdeltaRparticle -> Write();


// z 

    hz1pTreco -> Write();
    hz2pTreco -> Write();
    hz1mreco -> Write();
    hz2mreco -> Write();
    hz1cosThetareco -> Write();
    hz2cosThetareco -> Write();
    hz1phireco -> Write();
    hz2phireco -> Write();
    hzzdeltaPhireco -> Write();
    hz1etareco -> Write();
    hz2etareco -> Write();
    hzzdeltaEtareco -> Write();
    hz1Rreco -> Write();
    hz2Rreco -> Write();
    hzzdeltaRreco -> Write();

    hz1pTparticle -> Write();
    hz2pTparticle -> Write();
    hz1mparticle -> Write();
    hz2mparticle -> Write();
    hz1phiparticle -> Write();
    hz2phiparticle -> Write();
    hzzdeltaPhiparticle -> Write();
    hz1etaparticle -> Write();
    hz2etaparticle -> Write();
    hzzdeltaEtaparticle -> Write();
    hz1Rparticle -> Write();
    hz2Rparticle -> Write();
    hzzdeltaRparticle -> Write();

    hz1pTparton -> Write();
    hz2pTparton -> Write();
    hz1mparton -> Write();
    hz2mparton -> Write();
    hz1phiparton -> Write();
    hz2phiparton -> Write();
    hzzdeltaPhiparton -> Write();
    hz1etaparton -> Write();
    hz2etaparton -> Write();
    hzzdeltaEtaparton -> Write();
    hz1Rparton -> Write();
    hz2Rparton -> Write();
    hzzdeltaRparton -> Write();

// leptons 

    hl1pTreco -> Write();
    hl2pTreco -> Write();
    hl3pTreco -> Write();
    hl4pTreco -> Write();
    hl1phireco -> Write();
    hl2phireco -> Write();
    hl3phireco -> Write();
    hl4phireco -> Write();
    hl1l2deltaPhireco -> Write();
    hl3l4deltaPhireco -> Write();
    hl1l2deltaPhiBoostreco -> Write();
    hl3l4deltaPhiBoostreco -> Write();
    hl1etareco -> Write();
    hl2etareco -> Write();
    hl3etareco -> Write();
    hl4etareco -> Write();
    hl1l2deltaEtareco -> Write();
    hl3l4deltaEtareco -> Write();
    hl1l2deltaEtaBoostreco -> Write();
    hl3l4deltaEtaBoostreco -> Write();
    hl1Rreco -> Write();
    hl2Rreco -> Write();
    hl3Rreco -> Write();
    hl4Rreco -> Write();
    hl1l2deltaRreco -> Write();
    hl3l4deltaRreco -> Write();
    hl1cosThetareco -> Write();
    hl2cosThetareco -> Write();
    hl3cosThetareco -> Write();
    hl4cosThetareco -> Write();
    hfourlcosThetareco -> Write();
    hl1cosThetaBoostreco -> Write();
    hl2cosThetaBoostreco -> Write();
    hl3cosThetaBoostreco -> Write();
    hl4cosThetaBoostreco -> Write();
    hfourlcosThetaBoostreco -> Write();
    hl1l2CScosThetareco -> Write();
    hl3l4CScosThetareco -> Write();

    hl1pTparticle -> Write();
    hl2pTparticle -> Write();
    hl3pTparticle -> Write();
    hl4pTparticle -> Write();
    hl1phiparticle -> Write();
    hl2phiparticle -> Write();
    hl3phiparticle -> Write();
    hl4phiparticle -> Write();
    hl1l2deltaPhiparticle -> Write();
    hl3l4deltaPhiparticle -> Write();
    hl1l2deltaPhiBoostparticle -> Write();
    hl3l4deltaPhiBoostparticle -> Write();
    hl1etaparticle -> Write();
    hl2etaparticle -> Write();
    hl3etaparticle -> Write();
    hl4etaparticle -> Write();
    hl1l2deltaEtaparticle -> Write();
    hl3l4deltaEtaparticle -> Write();
    hl1l2deltaEtaBoostparticle -> Write();
    hl3l4deltaEtaBoostparticle -> Write();
    hl1Rparticle -> Write();
    hl2Rparticle -> Write();
    hl3Rparticle -> Write();
    hl4Rparticle -> Write();
    hl1l2deltaRparticle -> Write();
    hl3l4deltaRparticle -> Write();
    hl1cosThetaparticle -> Write();
    hl2cosThetaparticle -> Write();
    hl3cosThetaparticle -> Write();
    hl4cosThetaparticle -> Write();
    hfourlcosThetaparticle -> Write();
    hl1cosThetaBoostparticle -> Write();
    hl2cosThetaBoostparticle -> Write();
    hl3cosThetaBoostparticle -> Write();
    hl4cosThetaBoostparticle -> Write();
    hfourlcosThetaBoostparticle -> Write();
    hl1l2CScosThetaparticle -> Write();
    hl3l4CScosThetaparticle -> Write();

    hl1pTparton -> Write();
    hl2pTparton -> Write();
    hl3pTparton -> Write();
    hl4pTparton -> Write();
    hl1phiparton -> Write();
    hl2phiparton -> Write();
    hl3phiparton -> Write();
    hl4phiparton -> Write();
    hl1l2deltaPhiparton -> Write();
    hl3l4deltaPhiparton -> Write();
    hl1l2deltaPhiBoostparton -> Write();
    hl3l4deltaPhiBoostparton -> Write();
    hl1etaparton -> Write();
    hl2etaparton -> Write();
    hl3etaparton -> Write();
    hl4etaparton -> Write();
    hl1l2deltaEtaparton -> Write();
    hl3l4deltaEtaparton -> Write();
    hl1l2deltaEtaBoostparton -> Write();
    hl3l4deltaEtaBoostparton -> Write();
    hl1Rparton -> Write();
    hl2Rparton -> Write();
    hl3Rparton -> Write();
    hl4Rparton -> Write();
    hl1l2deltaRparton -> Write();
    hl3l4deltaRparton -> Write();
    hl1cosThetaparton -> Write();
    hl2cosThetaparton -> Write();
    hl3cosThetaparton -> Write();
    hl4cosThetaparton -> Write();
    hfourlcosThetaparton -> Write();
    hl1cosThetaBoostparton -> Write();
    hl2cosThetaBoostparton -> Write();
    hl3cosThetaBoostparton -> Write();
    hl4cosThetaBoostparton -> Write();
    hfourlcosThetaBoostparton -> Write();
    hl1l2CScosThetaparton -> Write();
    hl3l4CScosThetaparton -> Write();
   */

// comp

    hHpT12Comp -> Write();
    hHm12Comp  -> Write();
    hHpT23Comp -> Write();
    hHm23Comp -> Write();
    hHpT13Comp -> Write();
    hHm13Comp -> Write();
    hbbdeltaPhi12Comp -> Write();
    hbbdeltaPhi23Comp -> Write();
    hbbdeltaPhi13Comp -> Write();
    hbbdeltaEta12Comp -> Write();
    hbbdeltaEta23Comp -> Write();
    hbbdeltaEta13Comp -> Write();

    hjjpT12Comp -> Write();
    hjjpT23Comp -> Write();
    hjjpT13Comp -> Write();

    hjjdeltaPhi12Comp -> Write();
    hjjdeltaPhi23Comp -> Write();
    hjjdeltaPhi13Comp -> Write();

    hj1Phi23Comp -> Write();
    hj2Phi23Comp -> Write();
    

    hz1pT13Comp -> Write();
    hz1m13Comp -> Write();
    hz2pT13Comp -> Write();
    hz2m13Comp -> Write();
    hz1pT23Comp -> Write();
    hz1m23Comp -> Write();
    hz2pT23Comp -> Write();
    hz2m23Comp -> Write();
    hz1pT12Comp -> Write();
    hz1m12Comp -> Write();
    hz2pT12Comp -> Write();
    hz2m12Comp -> Write();

    hl1l2deltaPhiHpTcompreco -> Write();
    hl3l4deltaPhiHpTcompreco -> Write();

    hHz1pTcompreco -> Write();
    hHz2pTcompreco -> Write();
    hHz1pTcompparticle -> Write();
    hHz2pTcompparticle -> Write();
    hHz1pTcompparton -> Write();
    hHz2pTcompparton -> Write();

    hHpTz1etacompreco -> Write();
    hHpTz2etacompreco -> Write();
    hHpTz1etacompparticle -> Write();
    hHpTz2etacompparticle -> Write();
    hHpTz1etacompparton -> Write();
    hHpTz2etacompparton -> Write();

    hHzzpTcompreco -> Write();
    hHzzpTcompparticle -> Write();
    hHzzpTcompparton -> Write();

    kappaLambda -> Write();

    hists->Close();
    // gROOT->SetBatch(kTRUE);

//------------------------------------------------------------------------------------------------------------------------------------------------------------
// draw histos
//------------------------------------------------------------------------------------------------------------------------------------------------------------

// higgs
/*
    draw_hist(hHpTreco,"hbb_pT_reco", "p^{T}_{hbb}_reco", "pT (GeV)");
    draw_hist(hHmreco,"hbb_m_reco", "m_{hbb}_reco", "mass (GeV)");
    draw_hist(hb1pTreco,"b1_pT_reco", "p^{T}_{b1}_reco", "pT (GeV)");
    draw_hist(hb1mreco,"b1_m_reco", "m_{b1}_reco", "mass (GeV)");
    draw_hist(hb2pTreco,"b2_pT_reco", "p^{T}_{b2}_reco", "pT (GeV)");
    draw_hist(hb2mreco,"b2_m_reco", "m_{b2}_reco", "mass (GeV)");
    draw_hist(hHphireco,"hbb_phi_reco", "#phi_{hbb}_reco", "#phi");
    draw_hist(hb1phireco,"b1_phi_reco", "#phi_{b1}_reco", "#phi");
    draw_hist(hb2phireco,"b2_phi_reco", "#phi_{b2}_reco", "#phi");
    draw_hist(hbbdeltaPhireco,"bb_Deltaphi_reco", "#Delta#phi_{bb}_reco", "#Delta#phi");
    draw_hist(hHetareco,"hbb_eta_reco", "#eta_{hbb}_reco", "#eta");
    draw_hist(hb1etareco,"b1_eta_reco", "#eta_{b1}_reco", "#eta");
    draw_hist(hb2etareco,"b2_eta_reco", "#eta_{b2}_reco", "#eta");
    draw_hist(hbbdeltaEtareco,"bb_Deltaeta_reco", "#Delta#eta_{bb}_reco", "#Delta#eta");
    draw_hist(hHRreco,"hbb_R_reco", "R_{hbb}_reco", "R");
    draw_hist(hb1Rreco,"b1_R_reco", "R_{b1}_reco", "R");
    draw_hist(hb2Rreco,"b2_R_reco", "R_{b2}_reco", "R");
    draw_hist(hbbdeltaRreco,"bb_DeltaR_reco", "#DeltaR_{bb}_reco", "#DeltaR");

    draw_hist(hHpTparticle,"hbb_pT_particle", "p^{T}_{hbb}_particle", "pT (GeV)");
    draw_hist(hHmparticle,"hbb_m_particle", "m_{hbb}_particle", "mass (GeV)");
    draw_hist(hb1pTparticle,"b1_pT_particle", "p^{T}_{b1}_particle", "pT (GeV)");
    draw_hist(hb1mparticle,"b1_m_particle", "m_{b1}_particle", "mass (GeV)");
    draw_hist(hb2pTparticle,"b2_pT_particle", "p^{T}_{b2}_particle", "pT (GeV)");
    draw_hist(hb2mparticle,"b2_m_particle", "m_{b2}_particle", "mass (GeV)");
    draw_hist(hHphiparticle,"hbb_#phi_particle", "#phi_{hbb}_particle", "#phi");
    draw_hist(hb1phiparticle,"b1_#phi_particle", "#phi_{b1}_particle", "#phi");
    draw_hist(hb2phiparticle,"b2_#phi_particle", "#phi_{b2}_particle", "#phi");
    draw_hist(hbbdeltaPhiparticle,"bb_#Delta#phi_particle", "#Delta#phi_{bb}_particle", "#Delta#phi");
    draw_hist(hHetaparticle,"hbb_#eta_particle", "#eta_{hbb}_particle", "#eta");
    draw_hist(hb1etaparticle,"b1_#eta_particle", "#eta_{b1}_particle", "#eta");
    draw_hist(hb2etaparticle,"b2_#eta_particle", "#eta_{b2}_particle", "#eta");
    draw_hist(hbbdeltaEtaparticle,"bb_#Delta#eta_particle", "#Delta#eta_{bb}_particle", "#Delta#eta");
    draw_hist(hHRparticle,"hbb_R_particle", "R_{hbb}_particle", "R");
    draw_hist(hb1Rparticle,"b1_R_particle", "R_{b1}_particle", "R");
    draw_hist(hb2Rparticle,"b2_R_particle", "R_{b2}_particle", "R");
    draw_hist(hbbdeltaRparticle,"bb_#DeltaR_particle", "#DeltaR_{bb}_particle", "#DeltaR");

    draw_hist(hHpTparton,"hbb_pT_parton", "p^{T}_{hbb}_parton", "pT (GeV)");
    draw_hist(hHmparton,"hbb_m_parton", "m_{hbb}_parton", "mass (GeV)");
    draw_hist(hb1pTparton,"b1_pT_parton", "p^{T}_{b1}_parton", "pT (GeV)");
    draw_hist(hb1mparton,"b1_m_parton", "m_{b1}_parton", "mass (GeV)");
    draw_hist(hb2pTparton,"b2_pT_parton", "p^{T}_{b2}_parton", "pT (GeV)");
    draw_hist(hb2mparton,"b2_m_parton", "m_{b2}_parton", "mass (GeV)");
    draw_hist(hHphiparton,"hbb_#phi_parton", "#phi_{hbb}_parton", "#phi");
    draw_hist(hb1phiparton,"b1_#phi_parton", "#phi_{b1}_parton", "#phi");
    draw_hist(hb2phiparton,"b2_#phi_parton", "#phi_{b2}_parton", "#phi");
    draw_hist(hbbdeltaPhiparton,"bb_#Delta#phi_parton", "#Delta#phi_{bb}_parton", "#Delta#phi");
    draw_hist(hHetaparton,"hbb_#eta_parton", "#eta_{hbb}_parton", "#eta");
    draw_hist(hb1etaparton,"b1_#eta_parton", "#eta_{b1}_parton", "#eta");
    draw_hist(hb2etaparton,"b2_#eta_parton", "#eta_{b2}_parton", "#eta");
    draw_hist(hbbdeltaEtaparton,"bb_#Delta#eta_parton", "#Delta#eta_{bb}_parton", "#Delta#eta");
    draw_hist(hHRparton,"hbb_R_parton", "R_{hbb}_parton", "R");
    draw_hist(hb1Rparton,"b1_R_parton", "R_{b1}_parton", "R");
    draw_hist(hb2Rparton,"b2_R_parton", "R_{b2}_parton", "R");
    draw_hist(hbbdeltaRparton,"bb_#DeltaR_parton", "#DeltaR_{bb}_parton", "#DeltaR");

// jets

    draw_hist(hjjpTreco,"jj_pT_reco", "p^{T}_{jj}_reco", "pT (GeV)");
    draw_hist(hj1pTreco,"j1_pT_reco", "p^{T}_{j1}_reco", "pT (GeV)");
    draw_hist(hj2pTreco,"j2_pT_reco", "p^{T}_{j2}_reco", "pT (GeV)");

    draw_hist(hj1phireco,"j1_#phi_reco", "#phi_{j1}_reco", "#phi");
    draw_hist(hj2phireco,"j2_#phi_reco", "#phi_{j2}_reco", "#phi");

    draw_hist(hjjdeltaPhireco,"jj_#Delta#phi_reco", "#Delta#phi_{jj}_reco", "#Delta#phi");
    draw_hist(hj1etareco,"j1_#eta_reco", "#eta_{j1}_reco", "#eta");
    draw_hist(hj2etareco,"j2_#eta_reco", "#eta_{j2}_reco", "#eta");
    draw_hist(hjjdeltaEtareco,"jj_#Delta#eta_reco", "#Delta#eta_{jj}_reco", "#Delta#eta");
    draw_hist(hj1Rreco,"j1_R_reco", "R_{j1}_reco", "R");
    draw_hist(hj2Rreco,"j2_R_reco", "R_{j2}_reco", "R");
    draw_hist(hjjdeltaRreco,"jj_#DeltaR_reco", "#DeltaR_{jj}_reco", "#DeltaR");

    draw_hist(hjjpTparticle,"jj_pT_particle", "p^{T}_{jj}_particle", "pT (GeV)");
    draw_hist(hj1pTparticle,"j1_pT_particle", "p^{T}_{j1}_particle", "pT (GeV)");
    draw_hist(hj2pTparticle,"j2_pT_particle", "p^{T}_{j2}_particle", "pT (GeV)");

    draw_hist(hj1phiparticle,"j1_#phi_particle", "#phi_{j1}_particle", "#phi");
    draw_hist(hj2phiparticle,"j2_#phi_particle", "#phi_{j2}_particle", "#phi");

    draw_hist(hjjdeltaPhiparticle,"jj_#Delta#phi_particle", "#Delta#phi_{jj}_particle", "#Delta#phi");
    draw_hist(hj1etaparticle,"j1_#eta_particle", "#eta_{j1}_particle", "#eta");
    draw_hist(hj2etaparticle,"j2_#eta_particle", "#eta_{j2}_particle", "#eta");
    draw_hist(hjjdeltaEtaparticle,"jj_#Delta#eta_particle", "#Delta#eta_{jj}_particle", "#Delta#eta");
    draw_hist(hj1Rparticle,"j1_R_particle", "R_{j1}_particle", "R");
    draw_hist(hj2Rparticle,"j2_R_particle", "R_{j2}_particle", "R");
    draw_hist(hjjdeltaRparticle,"jj_#DeltaR_particle", "#DeltaR_{jj}_particle", "#DeltaR");

// z 

    draw_hist(hz1pTreco,"z1_pT_reco", "p^{T}_{z1}_reco", "pT (GeV)");
    draw_hist(hz2pTreco,"z2_pT_reco", "p^{T}_{z2}_reco", "pT (GeV)");
    draw_hist(hz1mreco,"z1_m_reco", "m_{z1}_reco", "pT (GeV)");
    draw_hist(hz2mreco,"z2_m_reco", "m_{z2}_reco", "pT (GeV)");
    draw_hist(hz1phireco,"z1_#phi_reco", "#phi_{z1}_reco", "#phi");
    draw_hist(hz2phireco,"z2_#phi_reco", "#phi_{z2}_reco", "#phi");
    draw_hist(hzzdeltaPhireco,"zz_#Delta#phi_reco", "#Delta#phi_{zz}_reco", "#Delta#phi");
    draw_hist(hz1etareco,"z1_#eta_reco", "#eta_{z1}_reco", "#eta");
    draw_hist(hz2etareco,"z2_#eta_reco", "#eta_{z2}_reco", "#eta");
    draw_hist(hzzdeltaEtareco,"zz_#Delta#eta_reco", "#Delta#eta_{zz}_reco", "#Delta#eta");
    draw_hist(hz1Rreco,"z1_R_reco", "R_{z1}_reco", "R");
    draw_hist(hz2Rreco,"z2_R_reco", "R_{z2}_reco", "R");
    draw_hist(hzzdeltaRreco,"zz_#DeltaR_reco", "#DeltaR_{zz}_reco", "#DeltaR");
    draw_hist(hz1cosThetareco,"z1_cos#theta_reco", "cos#theta_{z1}_reco", "cos#theta");
    draw_hist(hz2cosThetareco,"z2_cos#theta_reco", "cos#theta_{z1}_reco", "cos#theta");

    draw_hist(hz1pTparticle,"z1_pT_particle", "p^{T}_{z1}_particle", "pT (GeV)");
    draw_hist(hz2pTparticle,"z2_pT_particle", "p^{T}_{z2}_particle", "pT (GeV)");
    draw_hist(hz1mparticle,"z1_m_particle", "m_{z1}_particle", "mass (GeV)");
    draw_hist(hz2mparticle,"z2_m_particle", "m_{z2}_particle", "mass (GeV)");
    draw_hist(hz1phiparticle,"z1_#phi_particle", "#phi_{z1}_particle", "#phi");
    draw_hist(hz2phiparticle,"z2_#phi_particle", "#phi_{z2}_particle", "#phi");
    draw_hist(hzzdeltaPhiparticle,"zz_#Delta#phi_particle", "#Delta#phi_{zz}_particle", "#Delta#phi");
    draw_hist(hz1etaparticle,"z1_#eta_particle", "#eta_{z1}_particle", "#eta");
    draw_hist(hz2etaparticle,"z2_#eta_particle", "#eta_{z2}_particle", "#eta");
    draw_hist(hzzdeltaEtaparticle,"zz_#Delta#eta_particle", "#Delta#eta_{zz}_particle", "#Delta#eta");
    draw_hist(hz1Rparticle,"z1_R_particle", "R_{z1}_particle", "R");
    draw_hist(hz2Rparticle,"z2_R_particle", "R_{z2}_particle", "R");
    draw_hist(hzzdeltaRparticle,"zz_#DeltaR_particle", "#DeltaR_{zz}_particle", "#DeltaR");

    draw_hist(hz1pTparton,"z1_pT_parton", "p^{T}_{z1}_parton", "pT (GeV)");
    draw_hist(hz2pTparton,"z2_pT_parton", "p^{T}_{z2}_parton", "pT (GeV)");
    draw_hist(hz1mparton,"z1_m_parton", "m_{z1}_parton", "mass (GeV)");
    draw_hist(hz2mparton,"z2_m_parton", "m_{z2}_parton", "mass (GeV)");
    draw_hist(hz1phiparton,"z1_#phi_parton", "#phi_{z1}_parton", "#phi");
    draw_hist(hz2phiparton,"z2_#phi_parton", "#phi_{z2}_parton", "#phi");
    draw_hist(hzzdeltaPhiparton,"zz_#Delta#phi_parton", "#Delta#phi_{zz}_parton", "#Delta#phi");
    draw_hist(hz1etaparton,"z1_#eta_parton", "#eta_{z1}_parton", "#eta");
    draw_hist(hz2etaparton,"z2_#eta_parton", "#eta_{z2}_parton", "#eta");
    draw_hist(hzzdeltaEtaparton,"zz_#Delta#eta_parton", "#Delta#eta_{zz}_parton", "#Delta#eta");
    draw_hist(hz1Rparton,"z1_R_parton", "R_{z1}_parton", "R");
    draw_hist(hz2Rparton,"z2_R_parton", "R_{z2}_parton", "R");
    draw_hist(hzzdeltaRparton,"zz_#DeltaR_parton", "#DeltaR_{zz}_parton", "#DeltaR");

// leptons

    draw_hist(hl1pTreco,"l1_pT_reco", "p^{T}_{l1}_reco", "pT (GeV)");
    draw_hist(hl2pTreco,"l2_pT_reco", "p^{T}_{l2}_reco", "pT (GeV)");
    draw_hist(hl3pTreco,"l3_pT_reco", "p^{T}_{l3}_reco", "pT (GeV)");
    draw_hist(hl4pTreco,"l4_pT_reco", "p^{T}_{l4}_reco", "pT (GeV)");
    draw_hist(hl1phireco,"l1_#phi_reco", "#phi_{l1}_reco", "#phi");
    draw_hist(hl2phireco,"l2_#phi_reco", "#phi_{l2}_reco", "#phi");
    draw_hist(hl3phireco,"l3_#phi_reco", "#phi_{l3}_reco", "#phi");
    draw_hist(hl4phireco,"l4_#phi_reco", "#phi_{l4}_reco", "#phi");
    draw_hist(hl1l2deltaPhireco,"l1l2_#Delta#phi_reco", "#Delta#phi_{l1l2}_reco", "#Delta#phi");
    draw_hist(hl3l4deltaPhireco,"l3l4_#Delta#phi_reco", "#Delta#phi_{l3l4}_reco", "#Delta#phi");
    draw_hist(hl1l2deltaPhiBoostreco,"l1l2_#Delta#phi_Boost_reco", "#Delta#phi_{l1l2}_Boost_reco", "#Delta#phi");
    draw_hist(hl3l4deltaPhiBoostreco,"l3l4_#Delta#phi_Boost_reco", "#Delta#phi_{l3l4}_Boost_reco", "#Delta#phi");
    draw_hist(hl1etareco,"l1_#eta_reco", "#eta_{l1}_reco", "#eta");
    draw_hist(hl2etareco,"l2_#eta_reco", "#eta_{l2}_reco", "#eta");
    draw_hist(hl3etareco,"l3_#eta_reco", "#eta_{l3}_reco", "#eta");
    draw_hist(hl4etareco,"l4_#eta_reco", "#eta_{l4}_reco", "#eta");
    draw_hist(hl1l2deltaEtareco,"l1l2_#Delta#eta_reco", "#Delta#eta_{l1l2}_reco", "#Delta#eta");
    draw_hist(hl3l4deltaEtareco,"l3l4_#Delta#eta_reco", "#Delta#eta_{l3l4}_reco", "#Delta#eta");
    draw_hist(hl1l2deltaEtaBoostreco,"l1l2_#Delta#eta_Boost_reco", "#Delta#eta_{l1l2}_Boost_reco", "#Delta#eta");
    draw_hist(hl3l4deltaEtaBoostreco,"l3l4_#Delta#eta_Boost_reco", "#Delta#eta_{l3l4}_Boost_reco", "#Delta#eta");
    draw_hist(hl1Rreco,"l1_R_reco", "R_{l1}_reco", "R");
    draw_hist(hl2Rreco,"l2_R_reco", "R_{l2}_reco", "R");
    draw_hist(hl3Rreco,"l3_R_reco", "R_{l3}_reco", "R");
    draw_hist(hl4Rreco,"l4_R_reco", "R_{l4}_reco", "R");
    draw_hist(hl1l2deltaRreco,"l1l2_#Delta R_reco", "#DeltaR_{l1l2}_reco", "#DeltaR");
    draw_hist(hl3l4deltaRreco,"l3l4_#Delta R_reco", "#DeltaR_{l3l4}_reco", "#DeltaR");
    draw_hist(hl1l2CScosThetareco,"l1l2_cos#theta_{CS}_reco", "cos#theta_{CS l1l2}_reco", "cos#theta_{CS}");
    draw_hist(hl3l4CScosThetareco,"l3l4_cos#theta_{CS}_reco", "cos#theta_{CS l3l4}_reco", "cos#theta_{CS}");
    draw_hist(hl1cosThetareco,"l1_cos#theta_reco", "cos#theta_{l1}_reco", "cos#theta");
    draw_hist(hl2cosThetareco,"l2_cos#theta_reco", "cos#theta_{l2}_reco", "cos#theta");
    draw_hist(hl3cosThetareco,"l3_cos#theta_reco", "cos#theta_{l3}_reco", "cos#theta");
    draw_hist(hl4cosThetareco,"l4_cos#theta_reco", "cos#theta_{l4}_reco", "cos#theta");
    draw_hist(hfourlcosThetareco,"four_cos#theta_reco", "cos#theta_{fourl}_reco", "cos#theta");
    draw_hist(hl1cosThetaBoostreco,"l1_cos#theta_Boost_reco", "cos#theta_{l1}_Boost_reco", "cos#theta");
    draw_hist(hl2cosThetaBoostreco,"l2_cos#theta_Boost_reco", "cos#theta_{l2}_Boost_reco", "cos#theta");
    draw_hist(hl3cosThetaBoostreco,"l3_cos#theta_Boost_reco", "cos#theta_{l3}_Boost_reco", "cos#theta");
    draw_hist(hl4cosThetaBoostreco,"l4_cos#theta_Boost_reco", "cos#theta_{l4}_Boost_reco", "cos#theta");
    draw_hist(hfourlcosThetaBoostreco,"four_cos#theta_Boost_reco", "cos#theta_{fourl}_Boost_reco", "cos#theta");

    draw_hist(hl1pTparticle,"l1_pT_particle", "p^{T}_{l1}_particle", "pT (GeV)");
    draw_hist(hl2pTparticle,"l2_pT_particle", "p^{T}_{l2}_particle", "pT (GeV)");
    draw_hist(hl3pTparticle,"l3_pT_particle", "p^{T}_{l3}_particle", "pT (GeV)");
    draw_hist(hl4pTparticle,"l4_pT_particle", "p^{T}_{l4}_particle", "pT (GeV)");
    draw_hist(hl1phiparticle,"l1_#phi_particle", "#phi_{l1}_particle", "#phi");
    draw_hist(hl2phiparticle,"l2_#phi_particle", "#phi_{l2}_particle", "#phi");
    draw_hist(hl3phiparticle,"l3_#phi_particle", "#phi_{l3}_particle", "#phi");
    draw_hist(hl4phiparticle,"l4_#phi_particle", "#phi_{l4}_particle", "#phi");
    draw_hist(hl1l2deltaPhiparticle,"l1l2_#Delta#phi_particle", "#Delta#phi_{l1l2}_particle", "#Delta#phi");
    draw_hist(hl3l4deltaPhiparticle,"l3l4_#Delta#phi_particle", "#Delta#phi_{l3l4}_particle", "#Delta#phi");
    draw_hist(hl1l2deltaPhiBoostparticle,"l1l2_#Delta#phi_Boost_particle", "#Delta#phi_{l1l2}_Boost_particle", "#Delta#phi");
    draw_hist(hl3l4deltaPhiBoostparticle,"l3l4_#Delta#phi_Boost_particle", "#Delta#phi_{l3l4}_Boost_particle", "#Delta#phi");
    draw_hist(hl1etaparticle,"l1_#eta_particle", "#eta_{l1}_particle", "#eta");
    draw_hist(hl2etaparticle,"l2_#eta_particle", "#eta_{l2}_particle", "#eta");
    draw_hist(hl3etaparticle,"l3_#eta_particle", "#eta_{l3}_particle", "#eta");
    draw_hist(hl4etaparticle,"l4_#eta_particle", "#eta_{l4}_particle", "#eta");
    draw_hist(hl1l2deltaEtaparticle,"l1l2_#Delta#eta_particle", "#Delta#eta_{l1l2}_particle", "#Delta#eta");
    draw_hist(hl3l4deltaEtaparticle,"l3l4_#Delta#eta_particle", "#Delta#eta_{l3l4}_particle", "#Delta#eta");
    draw_hist(hl1l2deltaEtaBoostparticle,"l1l2_#Delta#eta_Boost_particle", "#Delta#eta_{l1l2}_Boost_particle", "#Delta#eta");
    draw_hist(hl3l4deltaEtaBoostparticle,"l3l4_#Delta#eta_Boost_particle", "#Delta#eta_{l3l4}_Boost_particle", "#Delta#eta");
    draw_hist(hl1Rparticle,"l1_R_particle", "R_{l1}_particle", "R");
    draw_hist(hl2Rparticle,"l2_R_particle", "R_{l2}_particle", "R");
    draw_hist(hl3Rparticle,"l3_R_particle", "R_{l3}_particle", "R");
    draw_hist(hl4Rparticle,"l4_R_particle", "R_{l4}_particle", "R");
    draw_hist(hl1l2deltaRparticle,"l1l2_#DeltaR_particle", "#DeltaR_{l1l2}_particle", "#Delta_R");
    draw_hist(hl3l4deltaRparticle,"l3l4_#DeltaR_particle", "#DeltaR_{l3l4}_particle", "#Delta_R");
    draw_hist(hl1l2CScosThetaparticle,"l1l2_cos#theta_{CS}_particle", "cos#theta_{CSl1l2}_particle", "cos#theta_{CS}");
    draw_hist(hl3l4CScosThetaparticle,"l3l4_cos#theta_{CS}_particle", "cos#theta_{CSl3l4}_particle", "cos#theta_{CS}");
    draw_hist(hl1cosThetaparticle,"l1_cos#theta_particle", "cos#theta_{l1}_particle", "cos#theta");
    draw_hist(hl2cosThetaparticle,"l2_cos#theta_particle", "cos#theta_{l2}_particle", "cos#theta");
    draw_hist(hl3cosThetaparticle,"l3_cos#theta_particle", "cos#theta_{l3}_particle", "cos#theta");
    draw_hist(hl4cosThetaparticle,"l4_cos#theta_particle", "cos#theta_{l4}_particle", "cos#theta");
    draw_hist(hfourlcosThetaparticle,"four_cos#theta_particle", "cos#theta_{fourl}_particle", "cos#theta");
    draw_hist(hl1cosThetaBoostparticle,"l1_cos#theta_Boost_particle", "cos#theta_{l1}_Boost_particle", "cos#theta");
    draw_hist(hl2cosThetaBoostparticle,"l2_cos#theta_Boost_particle", "cos#theta_{l2}_Boost_particle", "cos#theta");
    draw_hist(hl3cosThetaBoostparticle,"l3_cos#theta_Boost_particle", "cos#theta_{l3}_Boost_particle", "cos#theta");
    draw_hist(hl4cosThetaBoostparticle,"l4_cos#theta_Boost_particle", "cos#theta_{l4}_Boost_particle", "cos#theta");
    draw_hist(hfourlcosThetaBoostparticle,"four_cos#theta_Boost_particle", "cos#theta_{fourl}_Boost_particle", "cos#theta");

    draw_hist(hl1pTparton,"l1_pT_parton", "p^{T}_{l1}_parton", "pT (GeV)");
    draw_hist(hl2pTparton,"l2_pT_parton", "p^{T}_{l2}_parton", "pT (GeV)");
    draw_hist(hl3pTparton,"l3_pT_parton", "p^{T}_{l3}_parton", "pT (GeV)");
    draw_hist(hl4pTparton,"l4_pT_parton", "p^{T}_{l4}_parton", "pT (GeV)");
    draw_hist(hl1phiparton,"l1_#phi_parton", "#phi_{l1}_parton", "#phi");
    draw_hist(hl2phiparton,"l2_#phi_parton", "#phi_{l2}_parton", "#phi");
    draw_hist(hl3phiparton,"l3_#phi_parton", "#phi_{l3}_parton", "#phi");
    draw_hist(hl4phiparton,"l4_#phi_parton", "#phi_{l4}_parton", "#phi");
    draw_hist(hl1l2deltaPhiparton,"l1l2_#Delta#phi_parton", "#Delta#phi_{l1l2}_parton", "#Delta#phi");
    draw_hist(hl3l4deltaPhiparton,"l3l4_#Delta#phi_parton", "#Delta#phi_{l3l4}_parton", "#Delta#phi");
    draw_hist(hl1l2deltaPhiBoostparton,"l1l2_#Delta#phi_Boost_parton", "#Delta#phi_{l1l2}_Boost_parton", "#Delta#phi");
    draw_hist(hl3l4deltaPhiBoostparton,"l3l4_#Delta#phi_Boost_parton", "#Delta#phi_{l3l4}_Boost_parton", "#Delta#phi");
    draw_hist(hl1etaparton,"l1_#eta_parton", "#eta_{l1}_parton", "#eta");
    draw_hist(hl2etaparton,"l2_#eta_parton", "#eta_{l2}_parton", "#eta");
    draw_hist(hl3etaparton,"l3_#eta_parton", "#eta_{l3}_parton", "#eta");
    draw_hist(hl4etaparton,"l4_#eta_parton", "#eta_{l4}_parton", "#eta");
    draw_hist(hl1l2deltaEtaparton,"l1l2_#Delta#eta_parton", "#Delta#eta_{l1l2}_parton", "#Delta#eta");
    draw_hist(hl3l4deltaEtaparton,"l3l4_#Delta#eta_parton", "#Delta#eta_{l3l4}_parton", "#Delta#eta");
    draw_hist(hl1l2deltaEtaBoostparton,"l1l2_#Delta#eta_Boost_parton", "#Delta#eta_{l1l2}_Boost_parton", "#Delta#eta");
    draw_hist(hl3l4deltaEtaBoostparton,"l3l4_#Delta#eta_Boost_parton", "#Delta#eta_{l3l4}_Boost_parton", "#Delta#eta");
    draw_hist(hl1Rparton,"l1_R_parton", "R_{l1}_parton", "R");
    draw_hist(hl2Rparton,"l2_R_parton", "R_{l2}_parton", "R");
    draw_hist(hl3Rparton,"l3_R_parton", "R_{l3}_parton", "R");
    draw_hist(hl4Rparton,"l4_R_parton", "R_{l4}_parton", "R");
    draw_hist(hl1l2deltaRparton,"l1l2_#Delta R_parton", "#Delta_R_{l1l2}_parton", "#Delta_R");
    draw_hist(hl3l4deltaRparton,"l3l4_#Delta R_parton", "#Delta_R_{l3l4}_parton", "#Delta_R");
    draw_hist(hl1l2CScosThetaparton,"l1l2_cos#theta_{CS}_parton", "cos#theta_{CSl1l2}_parton", "cos#theta_{CS}");
    draw_hist(hl3l4CScosThetaparton,"l3l4_cos#theta_{CS}_parton", "cos#theta_{CSl3l4}_parton", "cos#theta_{CS}");
    draw_hist(hl1cosThetaparton,"l1_cos#theta_parton", "cos#theta_{l1}_parton", "cos#theta");
    draw_hist(hl2cosThetaparton,"l2_cos#theta_parton", "cos#theta_{l2}_parton", "cos#theta");
    draw_hist(hl3cosThetaparton,"l3_cos#theta_parton", "cos#theta_{l3}_parton", "cos#theta");
    draw_hist(hl4cosThetaparton,"l4_cos#theta_parton", "cos#theta_{l4}_parton", "cos#theta");
    draw_hist(hfourlcosThetaparton,"four_cos#theta_parton", "cos#theta_{fourl}_parton", "cos#theta");
    draw_hist(hl1cosThetaBoostparton,"l1_cos#theta_Boost_parton", "cos#theta_{l1}_Boost_parton", "cos#theta");
    draw_hist(hl2cosThetaBoostparton,"l2_cos#theta_Boost_parton", "cos#theta_{l2}_Boost_parton", "cos#theta");
    draw_hist(hl3cosThetaBoostparton,"l3_cos#theta_Boost_parton", "cos#theta_{l3}_Boost_parton", "cos#theta");
    draw_hist(hl4cosThetaBoostparton,"l4_cos#theta_Boost_parton", "cos#theta_{l4}_Boost_parton", "cos#theta");
    draw_hist(hfourlcosThetaBoostparton,"four_cos#theta_Boost_parton", "cos#theta_{fourl}_Boost_parton", "cos#theta");


// comp 
// 2D - parton(1) particle(2) reco(3)

    draw_hist2(hHpT12Comp, "H_pT_comp_12", "p^{T}_{h}", "parton level pT (GeV)", "particle level pT (GeV)");
    draw_hist2(hHm12Comp, "H_m_comp_12", "m_{h}", "parton level mass (GeV)", "particle level mass (GeV)");
    draw_hist2(hHpT23Comp, "H_pT_comp_23", "p^{T}_{h}", "particle level pT (GeV)", "reco level pT (GeV)");
    draw_hist2(hHm23Comp, "H_m_comp_23", "m_{h}", "particle level mass (GeV)", "reco level mass (GeV)");
    draw_hist2(hHpT13Comp, "H_pT_comp_13", "p^{T}_{h}", "parton level pT (GeV)", "reco level pT (GeV)");
    draw_hist2(hHm13Comp, "H_m_comp_13", "m_{h}", "parton level mass (GeV)", "reco level mass (GeV)");

    draw_hist2(hbbdeltaPhi12Comp, "bb_delta#phi_comp_12", "#Delta#phi_{bb}", "parton level #Delta#phi", "particle level #Delta#phi");
    draw_hist2(hbbdeltaPhi23Comp, "bb_delta#phi_comp_23", "#Delta#phi_{bb}", "particle level #Delta#phi", "reco level #Delta#phi");
    draw_hist2(hbbdeltaPhi13Comp, "bb_delta#phi_comp_13", "#Delta#phi_{bb}", "parton level #Delta#phi", "reco level #Delta#phi");
    draw_hist2(hbbdeltaEta12Comp, "bb_delta#eta_comp_12", "#Delta#eta_{bb}", "parton level #Delta#eta", "particle level #Delta#eta");
    draw_hist2(hbbdeltaEta23Comp, "bb_delta#eta_comp_23", "#Delta#eta_{bb}", "particle level #Delta#eta", "reco level #Delta#eta");
    draw_hist2(hbbdeltaEta13Comp, "bb_delta#eta_comp_13", "#Delta#eta_{bb}", "parton level #Delta#eta", "reco level #Delta#eta");

    draw_hist2(hjjpT12Comp, "jj_pT_comp_12", "p^{T}_{jj}", "parton level pT (GeV)", "particle level pT (GeV)");
    draw_hist2(hjjpT23Comp, "jj_pT_comp_23", "p^{T}_{jj}", "particle level pT (GeV)", "reco level pT (GeV)");
    draw_hist2(hjjpT13Comp, "jj_pT_comp_13", "p^{T}_{jj}", "partonlevel pT (GeV)", "reco level pT (GeV)");
  
    draw_hist2(hjjdeltaPhi12Comp, "jj_#Delta#phi_comp_12", "#Delta#phi_{jj}", "parton level #Delta#phi", "particle level #Delta#phi");
    draw_hist2(hjjdeltaPhi23Comp, "jj_#Delta#phi_comp_23", "#Delta#phi_{jj}", "particle level #Delta#phi", "reco level #Delta#phi");
    draw_hist2(hjjdeltaPhi13Comp, "jj_#Delta#phi_comp_13", "#Delta#phi_{jj}", "parton level #Delta#phi", "reco level #Delta#phi");
    
    draw_hist2(hj1Phi23Comp, "j1_#phi_comp_23", "#phi_{j1}", "particle level #phi", "reco level #phi");
    draw_hist2(hj2Phi23Comp, "j2_#phi_comp_23", "#phi_{j2}", "particle level #phi", "reco level #phi");

    draw_hist2(hl1l2deltaPhiHpTcompreco, "l1l2_delta#phi_H_pT_comp_reco", " ", "#Delta#phi_l1l2", "p^{T}_{h}_reco");
    draw_hist2(hl3l4deltaPhiHpTcompreco, "l3l4_delta#phi_H_pT_comp_reco", " ", "#Delta#phi_l1l2", "p^{T}_{h}_reco");
    draw_hist2(hz1pT12Comp, "z1_pT_comp_12", "p^{T}_{z1}", "parton level pT (GeV)", "particle level pT (GeV)");
    draw_hist2(hz1m12Comp, "z1_m_comp_12", "m_{z1}", "parton level mass (GeV)", "particle level mass (GeV)");
    draw_hist2(hz2pT12Comp, "z2_pT_comp_12", "p^{T}_{z2}", "parton level pT (GeV)", "particle level pT (GeV)");
    draw_hist2(hz2m12Comp, "z2_m_comp_12", "m_{z2}", "parton level mass (GeV)", "particle level mass (GeV)");
    draw_hist2(hz1pT23Comp, "z1_pT_comp_23", "p^{T}_{z1}", "particle level pT (GeV)", "reco level pT (GeV)");
    draw_hist2(hz1m23Comp, "z1_m_comp_23", "m_{z1}", "particle level mass (GeV)", "reco level mass (GeV)");
    draw_hist2(hz2pT23Comp, "z2_pT_comp_23", "p^{T}_{z2}", "particle level pT (GeV)", "reco level pT (GeV)");
    draw_hist2(hz2m23Comp, "z2_m_comp_23", "m_{z2}", "particle level mass (GeV)", "reco level mass (GeV)");
    draw_hist2(hz1pT13Comp, "z1_pT_comp_13", "p^{T}_{z1}", "parton level pT (GeV)", "reco level pT (GeV)");
    draw_hist2(hz1m13Comp, "z1_m_comp_13", "m_{z1}", "parton level mass (GeV)", "reco  level mass (GeV)");
    draw_hist2(hz2pT13Comp, "z2_pT_comp_13", "p^{T}_{z2}", "parton level pT (GeV)", "reco level pT (GeV)");
    draw_hist2(hz2m13Comp, "z2_m_comp_13", "m_{z2}", "parton level mass (GeV)", "reco  level mass (GeV)");

    draw_hist2(hHz1pTcompreco, "h_z1_pT_comp_reco", "h_z1_pT_comp_reco", "reco level h pT (GeV)", "reco level z1 pT (GeV)");
    draw_hist2(hHz2pTcompreco, "h_z2_pT_comp_reco", "h_z2_pT_comp_reco", "reco level h pT (GeV)", "reco level z2 pT (GeV)");
    draw_hist2(hHz1pTcompparticle, "h_z1_pT_comp_particle", "h_z1_pT_comp_particle", "particle level h pT (GeV)", "particle level z1 pT (GeV)");
    draw_hist2(hHz2pTcompparticle, "h_z2_pT_comp_particle", "h_z2_pT_comp_particle", "particle level h pT (GeV)", "particle level z2 pT (GeV)");
    draw_hist2(hHz1pTcompparton, "h_z1_pT_comp_parton", "h_z1_pT_comp_parton", "parton level h pT (GeV)", "parton level z1 pT (GeV)");
    draw_hist2(hHz2pTcompparton, "h_z2_pT_comp_parton", "h_z2_pT_comp_parton", "parton level h pT (GeV)", "parton level z2 pT (GeV)");

    draw_hist2(hHpTz1etacompreco, "h_pT_z1_#eta_comp_reco", "h_pT_z1_#eta_comp_reco", "reco level h pT (GeV)", "reco level z1 #eta");
    draw_hist2(hHpTz2etacompreco, "h_pT_z2_#eta_comp_reco", "h_pT_z2_#eta_comp_reco", "reco level h pT (GeV)", "reco level z2 #eta");
    draw_hist2(hHpTz1etacompparticle, "h_pT_z1_#eta_comp_particle", "h_pT_z1_#eta_comp_particle", "particle level h pT (GeV)", "particle level z1 #eta");
    draw_hist2(hHpTz2etacompparticle, "h_pT_z2_#eta_comp_particle", "h_pT_z2_#eta_comp_particle", "particle level h pT (GeV)", "particle level z2 #eta");
    draw_hist2(hHpTz1etacompparton, "h_pT_z1_#eta_comp_parton", "h_pT_z1_#eta_comp_parton", "parton level h pT (GeV)", "parton level z1 #eta");
    draw_hist2(hHpTz2etacompparton, "h_pT_z2_#eta_comp_parton", "h_pT_z2_#eta_comp_parton", "parton level h pT (GeV)", "parton level z2 #eta");

    draw_hist2(hHzzpTcompreco, "h_zz_pT_comp_reco", "h_zz_pT_comp_reco", "reco level h pT (GeV)", "reco level zz pT (GeV)");
    draw_hist2(hHzzpTcompparticle, "h_zz_pT_comp_particle", "h_zz_pT_comp_particle", "particle level h pT (GeV)", "particle level zz pT (GeV)");
    draw_hist2(hHzzpTcompparton , "h_zz_pT_comp_parton", "h_zz_pT_comp_parton", "parton level h pT (GeV)", "parton level zz pT (GeV)");


  draw_hist(kappaLambda," ", "crossX", "kappa_lambda");
*/
//------------------------------------------------------------------------------------------------------------------------------------------------------------
// delete ?
//------------------------------------------------------------------------------------------------------------------------------------------------------------

// higgs
    /*
  hHpTreco -> Clear();
  hHmreco -> Clear();
  hb1pTreco -> Clear();
  hb1mreco -> Clear();
  hb2pTreco -> Clear();
  hb2mreco -> Clear();
  hHphireco -> Clear();
  hb1phireco -> Clear();
  hb2phireco -> Clear();
  hbbdeltaPhireco -> Clear();
  hHetareco -> Clear();
  hb1etareco -> Clear();
  hb2etareco -> Clear();
  hbbdeltaEtareco -> Clear();
  hHRreco -> Clear();
  hb1Rreco -> Clear();
  hb2Rreco -> Clear();
  hbbdeltaRreco -> Clear();

  hHpTparticle -> Clear();
  hHmparticle -> Clear();
  hb1pTparticle -> Clear();
  hb1mparticle -> Clear();
  hb2pTparticle -> Clear();
  hb2mparticle -> Clear();
  hHphiparticle -> Clear();
  hb1phiparticle -> Clear();
  hb2phiparticle -> Clear();
  hbbdeltaPhiparticle -> Clear();
  hHetaparticle -> Clear();
  hb1etaparticle -> Clear();
  hb2etaparticle -> Clear();
  hbbdeltaEtaparticle -> Clear();
  hHRparticle -> Clear();
  hb1Rparticle -> Clear();
  hb2Rparticle -> Clear();
  hbbdeltaRparticle -> Clear();

  hHpTparton -> Clear();
  hHmparton -> Clear();
  hb1pTparton -> Clear();
  hb1mparton -> Clear();
  hb2pTparton -> Clear();
  hb2mparton -> Clear();
  hHphiparton -> Clear();
  hb1phiparton -> Clear();
  hb2phiparton -> Clear();
  hbbdeltaPhiparton -> Clear();
  hHetaparton -> Clear();
  hb1etaparton -> Clear();
  hb2etaparton -> Clear();
  hbbdeltaEtaparton -> Clear();
  hHRparton -> Clear();
  hb1Rparton -> Clear();
  hb2Rparton -> Clear();
  hbbdeltaRparton -> Clear();

// jets
  hjjpTreco -> Clear();
  hj1pTreco -> Clear();
  hj2pTreco -> Clear();
  hj1phireco -> Clear();
  hj2phireco -> Clear();
  hjjdeltaPhireco -> Clear();
  hj1etareco -> Clear();
  hj2etareco -> Clear();
  hjjdeltaEtareco -> Clear();
  hj1Rreco -> Clear();
  hj2Rreco -> Clear();
  hjjdeltaRreco -> Clear();

  hjjpTparticle -> Clear();
  hj1pTparticle -> Clear();
  hj2pTparticle -> Clear();
  hj1phiparticle -> Clear();
  hj2phiparticle -> Clear();
  hjjdeltaPhiparticle -> Clear();
  hj1etaparticle -> Clear();
  hj2etaparticle -> Clear();
  hjjdeltaEtaparticle -> Clear();
  hj1Rparticle -> Clear();
  hj2Rparticle -> Clear();
  hjjdeltaRparticle -> Clear();

// z 
  hz1pTreco -> Clear();
  hz2pTreco -> Clear();
  hz1mreco -> Clear();
  hz2mreco -> Clear();

  hz1phireco -> Clear();
  hz2phireco -> Clear();
  hzzdeltaPhireco -> Clear();

  hz1etareco -> Clear();
  hz2etareco -> Clear();
  hzzdeltaEtareco -> Clear();

  hz1Rreco -> Clear();
  hz2Rreco -> Clear();
  hzzdeltaRreco -> Clear();

  hz1cosThetareco -> Clear();
  hz2cosThetareco -> Clear();

  hz1pTparticle -> Clear();
  hz2pTparticle -> Clear();
  hz1mparticle -> Clear();
  hz2mparticle -> Clear();

  hz1phiparticle -> Clear();
  hz2phiparticle -> Clear();
  hzzdeltaPhiparticle -> Clear();

  hz1etaparticle -> Clear();
  hz2etaparticle -> Clear();
  hzzdeltaEtaparticle -> Clear();

  hz1Rparticle -> Clear();
  hz2Rparticle -> Clear();
  hzzdeltaRparticle -> Clear();

  hz1pTparton -> Clear();
  hz2pTparton -> Clear();
  hz1mparton -> Clear();
  hz2mparton -> Clear();

  hz1phiparton -> Clear();
  hz2phiparton -> Clear();
  hzzdeltaPhiparton -> Clear();

  hz1etaparton -> Clear();
  hz2etaparton -> Clear();
  hzzdeltaEtaparton -> Clear();

  hz1Rparton -> Clear();
  hz2Rparton -> Clear();
  hzzdeltaRparton -> Clear();

// leptons
  hl1pTreco -> Clear();
  hl2pTreco -> Clear();
  hl3pTreco -> Clear();
  hl4pTreco -> Clear();
  hl1phireco -> Clear();
  hl2phireco -> Clear();
  hl3phireco -> Clear();
  hl4phireco -> Clear();
  hl1l2deltaPhireco -> Clear();
  hl3l4deltaPhireco -> Clear();
  hl1l2deltaPhiBoostreco -> Clear();
  hl3l4deltaPhiBoostreco -> Clear();

  hl1etareco -> Clear();
  hl2etareco -> Clear();
  hl3etareco -> Clear();
  hl4etareco -> Clear();
  hl1l2deltaEtareco -> Clear();
  hl3l4deltaEtareco -> Clear();
  hl1l2deltaEtaBoostreco -> Clear();
  hl3l4deltaEtaBoostreco -> Clear();

  hl1Rreco -> Clear();
  hl2Rreco -> Clear();
  hl3Rreco -> Clear();
  hl4Rreco -> Clear();
  hl1l2deltaRreco -> Clear();
  hl3l4deltaRreco -> Clear();

  hl1cosThetareco -> Clear();
  hl2cosThetareco -> Clear();
  hl3cosThetareco -> Clear();
  hl4cosThetareco -> Clear();
  hfourlcosThetareco -> Clear();
  hl1cosThetaBoostreco -> Clear();
  hl2cosThetaBoostreco -> Clear();
  hl3cosThetaBoostreco -> Clear();
  hl4cosThetaBoostreco -> Clear();
  hfourlcosThetaBoostreco -> Clear();
  hl1l2CScosThetareco -> Clear();
  hl3l4CScosThetareco -> Clear();

  // particle
  hl1pTparticle -> Clear();
  hl2pTparticle -> Clear();
  hl3pTparticle -> Clear();
  hl4pTparticle -> Clear();

  hl1phiparticle -> Clear();
  hl2phiparticle -> Clear();
  hl3phiparticle -> Clear();
  hl4phiparticle -> Clear();
  hl1l2deltaPhiparticle -> Clear();
  hl3l4deltaPhiparticle -> Clear();
  hl1l2deltaPhiBoostparticle -> Clear();
  hl3l4deltaPhiBoostparticle -> Clear();

  hl1etaparticle -> Clear();
  hl2etaparticle -> Clear();
  hl3etaparticle -> Clear();
  hl4etaparticle -> Clear();
  hl1l2deltaEtaparticle -> Clear();
  hl3l4deltaEtaparticle -> Clear();
  hl1l2deltaEtaBoostparticle -> Clear();
  hl3l4deltaEtaBoostparticle -> Clear();

  hl1Rparticle -> Clear();
  hl2Rparticle -> Clear();
  hl3Rparticle -> Clear();
  hl4Rparticle -> Clear();
  hl1l2deltaRparticle -> Clear();
  hl3l4deltaRparticle -> Clear();

  hl1cosThetaparticle -> Clear();
  hl2cosThetaparticle -> Clear();
  hl3cosThetaparticle -> Clear();
  hl4cosThetaparticle -> Clear();
  hfourlcosThetaparticle -> Clear();
  hl1cosThetaBoostparticle -> Clear();
  hl2cosThetaBoostparticle -> Clear();
  hl3cosThetaBoostparticle -> Clear();
  hl4cosThetaBoostparticle -> Clear();
  hfourlcosThetaBoostparticle -> Clear();
  hl1l2CScosThetaparticle -> Clear();
  hl3l4CScosThetaparticle -> Clear();

// parton
  hl1pTparton -> Clear();
  hl2pTparton -> Clear();
  hl3pTparton -> Clear();
  hl4pTparton -> Clear();

  hl1phiparton -> Clear();
  hl2phiparton -> Clear();
  hl3phiparton -> Clear();
  hl4phiparton -> Clear();
  hl1l2deltaPhiparton -> Clear();
  hl3l4deltaPhiparton -> Clear();
  hl1l2deltaPhiBoostparton -> Clear();
  hl3l4deltaPhiBoostparton -> Clear();

  hl1etaparton -> Clear();
  hl2etaparton -> Clear();
  hl3etaparton -> Clear();
  hl4etaparton -> Clear();
  hl1l2deltaEtaparton -> Clear();
  hl3l4deltaEtaparton -> Clear();
  hl1l2deltaEtaBoostparton -> Clear();
  hl3l4deltaEtaBoostparton -> Clear();

  hl1Rparton -> Clear();
  hl2Rparton -> Clear();
  hl3Rparton -> Clear();
  hl4Rparton -> Clear();
  hl1l2deltaRparton -> Clear();
  hl3l4deltaRparton -> Clear();

  hl1cosThetaparton -> Clear();
  hl2cosThetaparton -> Clear();
  hl3cosThetaparton -> Clear();
  hl4cosThetaparton -> Clear();
  hfourlcosThetaparton -> Clear();
  hl1cosThetaBoostparton -> Clear();
  hl2cosThetaBoostparton -> Clear();
  hl3cosThetaBoostparton -> Clear();
  hl4cosThetaBoostparton -> Clear();
  hfourlcosThetaBoostparton -> Clear();
  hl1l2CScosThetaparton -> Clear();
  hl3l4CScosThetaparton -> Clear();

  hHpT12Comp -> Clear();
  hHm12Comp -> Clear();
  hHpT23Comp -> Clear();
  hHm23Comp -> Clear();
  hHpT13Comp -> Clear();
  hHm13Comp -> Clear();
  hbbdeltaPhi12Comp -> Clear();
  hbbdeltaPhi23Comp -> Clear();
  hbbdeltaPhi13Comp -> Clear();
  hbbdeltaEta12Comp -> Clear();
  hbbdeltaEta23Comp -> Clear();
  hbbdeltaEta13Comp -> Clear();

  hjjpT12Comp -> Clear();
  hjjpT23Comp -> Clear();
  hjjpT13Comp -> Clear();
  hj1Phi23Comp -> Clear();
  hj2Phi23Comp -> Clear();
  hjjdeltaPhi12Comp -> Clear();
  hjjdeltaPhi23Comp -> Clear();
  hjjdeltaPhi13Comp -> Clear();

  hz1pT12Comp -> Clear();
  hz1m12Comp -> Clear();
  hz2pT12Comp -> Clear();
  hz2m12Comp -> Clear();
  hz1pT23Comp -> Clear();
  hz1m23Comp -> Clear();
  hz2pT23Comp -> Clear();
  hz2m23Comp -> Clear();
  hz1pT13Comp -> Clear();
  hz1m13Comp -> Clear();
  hz2pT13Comp -> Clear();
  hz2m13Comp -> Clear();

  hl1l2deltaPhiHpTcompreco -> Clear();
  hl3l4deltaPhiHpTcompreco -> Clear();

  hHz1pTcompreco -> Clear();
  hHz2pTcompreco -> Clear();
  hHz1pTcompparticle -> Clear();
  hHz2pTcompparticle -> Clear();
  hHz1pTcompparton -> Clear();
  hHz2pTcompparton -> Clear();

  hHpTz1etacompreco -> Clear();
  hHpTz2etacompreco -> Clear();
  hHpTz1etacompparticle -> Clear();
  hHpTz2etacompparticle -> Clear();
  hHpTz1etacompparton -> Clear();
  hHpTz2etacompparton -> Clear();

  hHzzpTcompreco -> Clear();
  hHzzpTcompparticle -> Clear();
  hHzzpTcompparton -> Clear();

  kappaLambda -> Clear();
    */



/*
  std::vector<std::pair<std::string, std::pair<int, double>>> sortedCutFlow(cutFlowMap.begin(), cutFlowMap.end());
  std::sort(sortedCutFlow.begin(), sortedCutFlow.end(), [](const auto& a, const auto& b) {
    return a.second.first > b.second.first;
  });
*/

   cout<<" Reco CutF Flow "<<endl;
   PrintCutFlow(cutFlowMap_reco,cutList_reco,  "Reco");
   cout<<" Particle  CutF Flow "<<endl;
   PrintCutFlow(cutFlowMap_particle,cutList_particle, "Particle");
   
	/*
  std::cout<<std::left<<std::setw(25)<<"Reco Cut"<<std::setw(10)<<"Reco Passed"<<std::setw(15)<<" Rel Eff "<< std::setw(15)<<"Reco Efficiency"<< std::endl;
    for(int i=0; i<(int) cutList_reco.size(); i++) {
        const std::string cutName_reco = cutList_reco[i];
        double passed_reco =  cutFlowMap_reco[cutName_reco].first;
        double efficiency_reco = 100.00 * cutFlowMap_reco[cutName_reco].second / cutFlowMap_reco["initial reco"].second;

	double relEff= 100;;
	if(i>0)
	  relEff = 100.00 * cutFlowMap_reco[cutName_reco].second / cutFlowMap_reco[cutList_reco.at(i-1)].second; 
	std::cout<<std::left<<std::setw(25)<<cutName_reco<<std::setw(10)<<passed_reco<< std::setw(15)<<relEff <<std::setw(15)<<efficiency_reco<< std::endl;
    }

    std::cout<<std::left<<std::setw(25)<<"Particle Cut"<<std::setw(10)<<"Particle Passed"<<std::setw(15)<< "Rel Eff "<<std::setw(15)<<"Particle Efficiency"<< std::endl;
    cout<<endl;
    
    for(int i=0; i<(int) cutList_particle.size(); i++) {
        const std::string cutName_particle = cutList_particle[i];
        double passed_particle =  cutFlowMap_particle[cutName_particle].first ;
        double efficiency_particle = 100.00 * cutFlowMap_particle[cutName_particle].second / cutFlowMap_particle["initial particle"].second;

	double relEff= 100;;
	if(i>0)
	  relEff = 100.00 * cutFlowMap_particle[cutName_particle].second / cutFlowMap_particle[cutList_particle.at(i-1)].second; 
	
        std::cout<<std::left<<std::setw(25)<<cutName_particle<<std::setw(10)<<passed_particle<<std::setw(15)<<relEff<<std::setw(15)<<efficiency_particle<< std::endl;
  }
	*/

  std::cout<<"Total number of entries "<<numberOfEntries<<" Passed "<<nPassed<<" raw "<<nPassedRaw<<std::endl;


  
}

int main(int argc, char* argv[]) {
  const char *inputFileName = argv[1];
  const char *outputFileName = argv[2];
  //  const char *process_name = argv[3];                                                                                                                                                                                                                                                                                                                                   
  zAnalyzer(inputFileName, outputFileName);
  return 1;
}
