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

//#define MDEBUG
#define MSEED 1234 

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

#include "includes/cutflow_include.h"
#include "includes/crossx_include.h"
#include "includes/hist_include.h"
#include "includes/weights_include.h"
#include "includes/selections_include.h"
#include "includes/kinematics_include.h"

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


// void zAnalyzer(const char *inputFile,const char *outputFile, int kappaVal = 8) {
void zAnalyzer(const char *inputFile, const char *outputFile, const char *process_name, string analysis="HZZJJ"){
    
  #ifdef __CLING__
    gSystem->Load("libDelphes");
  #endif

// cutflow table

  vector <TH1F*> listOfTH1;
  vector <TH2F*> listOfTH2;
  vector <TProfile*> listOfTProfiles;

  DefineSelections();
  
  cutList_reco=cutSelectionProcessReco[analysis];
  cutList_particle=cutSelectionProcessParticle[analysis];
  cutList_parton=cutSelectionProcessParton[analysis];

  std::map<string, bool> enableCutReco;
  std::map<string, bool> enableCutParticle;
  std::map<string, bool> enableCutParton;
 
  for(vector<string>::iterator it=cutSelectionProcessReco["all"].begin(); it!=cutSelectionProcessReco["all"].end(); it++){
    enableCutReco[ (*it)] = hasCut(cutList_reco, (*it));
  }
  for(vector<string>::iterator it=cutSelectionProcessParticle["all"].begin(); it!=cutSelectionProcessParticle["all"].end(); it++){
    enableCutParticle[ (*it)] = hasCut(cutList_particle, (*it));
  }
  for(vector<string>::iterator it=cutSelectionProcessParton["all"].begin(); it!=cutSelectionProcessParton["all"].end(); it++){
    enableCutParton[ (*it)] = hasCut(cutList_parton, (*it));
  }

  int cutVal_reco = 0;
  double cutValW_reco = 0;
 
  std::map<string, std::pair<int,double>> cutFlowMap_reco;
  for(int i=0; i<(int) cutList_reco.size(); i++) { 
    cutFlowMap_reco[cutList_reco.at(i)] = make_pair(0,0.0); 
  }
 
  int cutVal_particle = 0;
  double cutValW_particle = 0;
 
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

  typedef std::map<std::string, std::pair<int,double>> cutFlowMapDef;
  std::map<string, cutFlowMapDef* > cutFlowMapAll;
    cutFlowMapAll["reco"] =  & cutFlowMap_reco;
    cutFlowMapAll["particle"] = & cutFlowMap_particle;
    cutFlowMapAll["parton"] = & cutFlowMap_parton;
  
  for(std::vector<string>::iterator it=selType.begin(); it!=selType.end(); it++){
    cutFlowHists[(*it)]=new TH1F(Form("hSel_%s",(*it).c_str()),"",cutFlowMByType[(*it)].size(),0,cutFlowMByType[(*it)].size()+1);
    cutFlowEffs[(*it)]=new TProfile(Form("hEff_%s",(*it).c_str()),"",cutFlowMByType[(*it)].size(),0,cutFlowMByType[(*it)].size()+1);
    
    listOfTH1.push_back(cutFlowHists[(*it)]);
    listOfTProfiles.push_back((cutFlowEffs[(*it)]));
  }
  
// delphes

  TChain chain("Delphes");
  chain.Add(inputFile);

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
  TClonesArray *branchWeight   = treeReader->UseBranch("Weight");

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

  TFile *hists= new TFile(outputFile,"recreate");

//------------------------------------------------------------------------------------------------------------------------------------------------------------
// BOOK HISTOGRAMS
//------------------------------------------------------------------------------------------------------------------------------------------------------------


// def bins + range

  bool fill_1D = true;
  bool fill_2D = true;

  const double mBins = 10;
  const double pTBins = 10;
  const double phiBins = 10;
  const double etaBins = 10;
  const double RBins = 10;
  const double cosBins = 10;

  const double hpTmin = 0;
  const double jpTmin = 0;
  const double zpTmin = 0;
  const double wpTmin = 0;
  const double lpTmin = 0;
  const double METpTmin = 0;
  const double hpTmax = 500;
  const double jpTmax = 500;
  const double zpTmax = 500;
  const double wpTmax = 500;
  const double lpTmax = 100;
  const double METpTmax = 500;

  const double hmmin = 0;
  const double jmmin = 0;
  const double zmmin = 0;
  const double wmmin = 0;
  const double hmmax = 200;
  const double jmmax = 200;
  const double zmmax = 200;
  const double wmmax = 200;

  const double hetamin = -2.5;
  const double jetamin = -2.5;
  const double zetamin = -2.5;
  const double wetamin = -2.5;
  const double letamin = -2.5;
  const double METetamin = -2.5;
  const double hetamax = 2.5;
  const double jetamax = 2.5;
  const double zetamax = 2.5;
  const double wetamax = 2.5;
  const double letamax = 2.5;
  const double METetamax = 2.5;

  const double hRmin = 0;
  const double jRmin = 0;
  const double zRmin = 0;
  const double wRmin = 0;
  const double lRmin = 0;
  const double METRmin = 0;
  const double hRmax = 5;
  const double jRmax = 5;
  const double zRmax = 5;
  const double wRmax = 5;
  const double lRmax = 5;
  const double METRmax = 5;

// 1D

  // higgs
    TH1F *hHpTreco = new TH1F("Hbb_pT_reco", "p^{T}_{Hbb} Reco", pTBins, hpTmin, hpTmax); 	listOfTH1.push_back(hHpTreco);
    TH1F *hHmreco = new TH1F("Hbb_m_reco", "m_{Hbb} Reco", mBins, hmmin, hmmax); listOfTH1.push_back(hHmreco);
    TH1F *hbbdeltaPhireco = new TH1F("bb_deltaPhi_reco", "#Delta#phi_{bb} Reco", phiBins, -TMath::Pi(), +TMath::Pi()); listOfTH1.push_back(hbbdeltaPhireco);
    TH1F *hbbdeltaEtareco = new TH1F("bb_deltaEta_reco", "#Delta#eta_{bb} Reco", etaBins, hetamin, hetamax); listOfTH1.push_back(hbbdeltaEtareco);
    TH1F *hbbdeltaRreco = new TH1F("bb_deltaR_reco", "#DeltaR_{bb} Reco", RBins, hRmin, hRmax); listOfTH1.push_back(hbbdeltaRreco);

    TH1F *hHpTparticle = new TH1F("Hbb_pT_particle", "p^{T}_{Hbb} Particle", pTBins, hpTmin, hpTmax); listOfTH1.push_back(hHpTparticle);
    TH1F *hHmparticle = new TH1F("Hbb_m_particle", "m_{Hbb} Particle", mBins, hmmin, hmmax); listOfTH1.push_back(hHmparticle);
    TH1F *hbbdeltaPhiparticle = new TH1F("bb_deltaPhi_particle", "#Delta#phi_{bb} Particle", phiBins, -TMath::Pi(), +TMath::Pi()); listOfTH1.push_back(hbbdeltaPhiparticle);
    TH1F *hbbdeltaEtaparticle = new TH1F("bb_deltaEta_particle", "#Delta#eta_{bb} Particle", etaBins, hetamin, hetamax); listOfTH1.push_back(hbbdeltaEtaparticle);
    TH1F *hbbdeltaRparticle = new TH1F("bb_deltaR_particle", "#DeltaR_{bb} Particle", RBins, hRmin, hRmax); listOfTH1.push_back(hbbdeltaRparticle);

    TH1F *hHpTparton = new TH1F("Hbb_pT_parton", "p^{T}_{hbb} Parton", pTBins, hpTmin, hpTmax); listOfTH1.push_back(hHpTparton);
    TH1F *hHmparton = new TH1F("Hbb_m_parton", "m_{hbb} Parton", mBins, hmmin,  hmmax); listOfTH1.push_back(hHmparton);
    TH1F *hbbdeltaPhiparton = new TH1F("bb_deltaPhi_parton", "#Delta#phi_{bb} Parton", phiBins, -TMath::Pi(), +TMath::Pi()); listOfTH1.push_back(hbbdeltaPhiparton);
    TH1F *hbbdeltaEtaparton = new TH1F("bb_deltaEta_parton", "#Delta#eta_{bb} Parton", etaBins, hetamin, hetamax); listOfTH1.push_back(hbbdeltaEtaparton);
    TH1F *hbbdeltaRparton = new TH1F("bb_deltaR_parton", "#DeltaR_{bb} Parton", RBins, hRmin, hRmax); listOfTH1.push_back(hbbdeltaRparton);

  // vbf jets
    TH1F *hjjpTreco = new TH1F("jj_pT_reco", "p^{T}_{jj} Reco", pTBins, jpTmin, jpTmax); listOfTH1.push_back(hjjpTreco);
    TH1F *hjjdeltaPhireco = new TH1F("jj_deltaPhi_reco", "#Delta#phi_{jj} Reco", phiBins, -TMath::Pi(), +TMath::Pi()); listOfTH1.push_back(hjjdeltaPhireco);
    TH1F *hjjdeltaEtareco = new TH1F("jj_deltaEta_reco", "#Delta#eta_{jj} Reco", etaBins, jetamin, jetamax); listOfTH1.push_back(hjjdeltaEtareco);
    TH1F *hjjdeltaRreco = new TH1F("jj_deltaR_reco", "#DeltaR_{jj} Reco", RBins, jRmin, jRmax); listOfTH1.push_back(hjjdeltaRreco);

    TH1F *hjjpTparticle = new TH1F("jj_pT_particle", "p^{T}_{jj} Particle", pTBins, jpTmin, jpTmax); listOfTH1.push_back(hjjpTparticle);
    TH1F *hjjdeltaPhiparticle = new TH1F("jj_deltaPhi_particle", "#Delta#phi_{jj} Particle", phiBins, -TMath::Pi(), +TMath::Pi()); listOfTH1.push_back(hjjdeltaPhiparticle);
    TH1F *hjjdeltaEtaparticle = new TH1F("jj_deltaEta_particle", "#Delta#eta_{jj} Particle", etaBins, jetamin, jetamax); listOfTH1.push_back(hjjdeltaEtaparticle);
    TH1F *hjjdeltaRparticle = new TH1F("jj_deltaR_particle", "#DeltaR_{jj} Particle", RBins, jRmin, jRmax); listOfTH1.push_back(hjjdeltaRparticle);

  // z1
    TH1F *hZ1pTreco = new TH1F("Z1_pT_reco", "p^{T}_{Z1} Reco", pTBins, zpTmin, zpTmax); listOfTH1.push_back(hZ1pTreco);
    TH1F *hZ1mreco = new TH1F("Z1_m_reco", "m_{Z1} Reco", mBins, zmmin, zmmax); listOfTH1.push_back(hZ1mreco);

    TH1F *hZ1pTparticle = new TH1F("Z1_pT_particle", "p^{T}_{Z1} Particle", pTBins, zpTmin, zpTmax); listOfTH1.push_back(hZ1pTparticle);
    TH1F *hZ1mparticle = new TH1F("Z1_m_particle", "m_{Z1} Particle", mBins, zmmin, zmmax); listOfTH1.push_back(hZ1mparticle);

    TH1F *hZ1pTparton = new TH1F("Z1_pT_parton", "p^{T}_{Z1} Parton", pTBins, zpTmin, zpTmax); listOfTH1.push_back(hZ1pTparton);
    TH1F *hZ1mparton = new TH1F("Z1_m_parton", "m_{Z1} Parton", mBins, zmmin, zmmax); listOfTH1.push_back(hZ1mparton);

  // z2
    TH1F *hZ2pTreco = new TH1F("Z2_pT_reco", "p^{T}_{Z2} Reco", pTBins, zpTmin, zpTmax); listOfTH1.push_back(hZ2pTreco);
    TH1F *hZ2mreco = new TH1F("Z2_m_reco", "m_{Z2} Reco", mBins, zmmin, zmmax); listOfTH1.push_back(hZ2mreco);

    TH1F *hZ2pTparticle = new TH1F("Z2_pT_particle", "p^{T}_{Z2} Particle", pTBins, zpTmin, zpTmax); listOfTH1.push_back(hZ2pTparticle);
    TH1F *hZ2mparticle = new TH1F("Z2_m_particle", "m_{Z2} Particle", mBins, zmmin, zmmax); listOfTH1.push_back(hZ2mparticle);

    TH1F *hZ2pTparton = new TH1F("Z2_pT_parton", "p^{T}_{Z2} Parton", pTBins, zpTmin, zpTmax); listOfTH1.push_back(hZ2pTparton);
    TH1F *hZ2mparton = new TH1F("Z2_m_parton", "m_{Z2} Parton", mBins, zmmin, zmmax); listOfTH1.push_back(hZ2mparton);

  // zz
    TH1F *hZZdeltaPhireco = new TH1F("ZZ_#Delta#phi_reco", "#Delta#phi_{ZZ} Reco", phiBins, -TMath::Pi(), +TMath::Pi()); listOfTH1.push_back(hZZdeltaPhireco);
    TH1F *hZZdeltaEtareco = new TH1F("ZZ_#Delta#eta_reco", "#Delta#eta_{ZZ} Reco", etaBins, zetamin, zetamax);listOfTH1.push_back(hZZdeltaEtareco);
    TH1F *hZZdeltaRreco = new TH1F("ZZ_#DeltaR_reco", "#DeltaR_{ZZ} Reco", RBins, zRmin, zRmax); listOfTH1.push_back(hZZdeltaRreco);

    TH1F *hZZdeltaPhiparticle = new TH1F("ZZ_#Delta#phi_particle", "#Delta#phi_{ZZ} Particle", phiBins, -TMath::Pi(), +TMath::Pi()); listOfTH1.push_back(hZZdeltaPhiparticle);
    TH1F *hZZdeltaEtaparticle = new TH1F("ZZ_#Delta#eta_particle", "#Delta#eta_{ZZ} Particle", etaBins, zetamin, zetamax);listOfTH1.push_back(hZZdeltaEtaparticle);
    TH1F *hZZdeltaRparticle = new TH1F("ZZ_#DeltaR_particle", "#DeltaR_{ZZ} Particle", RBins, zRmin, zRmax); listOfTH1.push_back(hZZdeltaRparticle);

    TH1F *hZZdeltaPhiparton = new TH1F("ZZ_#Delta#phi_parton", "#Delta#phi_{ZZ} Parton", phiBins, -TMath::Pi(), +TMath::Pi()); listOfTH1.push_back(hZZdeltaPhiparton);
    TH1F *hZZdeltaEtaparton = new TH1F("ZZ_#Delta#eta_parton", "#Delta#eta_{ZZ} Parton", etaBins, zetamin, zetamax);listOfTH1.push_back(hZZdeltaEtaparton);
    TH1F *hZZdeltaRparton = new TH1F("ZZ_#DeltaR_parton", "#DeltaR_{ZZ} Parton", RBins, zRmin, zRmax); listOfTH1.push_back(hZZdeltaRparton);

  // w1
    TH1F *hW1pTreco = new TH1F("W1_pT_reco", "p^{T}_{W1} Reco", pTBins, wpTmin, wpTmax); listOfTH1.push_back(hW1pTreco);
    TH1F *hW1mreco = new TH1F("W1_m_reco", "m_{W1} Reco", mBins, wmmin, wmmax); listOfTH1.push_back(hW1mreco);

    TH1F *hW1pTparticle = new TH1F("W1_pT_particle", "p^{T}_{W1} Particle", pTBins, wpTmin, wpTmax); listOfTH1.push_back(hW1pTparticle);
    TH1F *hW1mparticle = new TH1F("W1_m_particle", "m_{W1} Particle", mBins, wmmin, wmmax); listOfTH1.push_back(hW1mparticle);

    TH1F *hW1pTparton = new TH1F("W1_pT_parton", "p^{T}_{W1} Parton", pTBins, wpTmin, wpTmax); listOfTH1.push_back(hW1pTparton);
    TH1F *hW1mparton = new TH1F("W1_m_parton", "m_{W1} Parton", mBins, wmmin, wmmax); listOfTH1.push_back(hW1mparton);

  // w2
    TH1F *hW2pTreco = new TH1F("W2_pT_reco", "p^{T}_{W2} Reco", pTBins, wpTmin, wpTmax); listOfTH1.push_back(hW2pTreco);
    TH1F *hW2mreco = new TH1F("W2_m_reco", "m_{W2} Reco", mBins, wmmin, wmmax); listOfTH1.push_back(hW2mreco);

    TH1F *hW2pTparticle = new TH1F("W2_pT_particle", "p^{T}_{W2} Particle", pTBins, wpTmin, wpTmax); listOfTH1.push_back(hW2pTparticle);
    TH1F *hW2mparticle = new TH1F("W2_m_particle", "m_{W2} Particle", mBins, wmmin, wmmax); listOfTH1.push_back(hW2mparticle);

    TH1F *hW2pTparton = new TH1F("W2_pT_parton", "p^{T}_{W2} Parton", pTBins, wpTmin, wpTmax); listOfTH1.push_back(hW2pTparton);
    TH1F *hW2mparton = new TH1F("W2_m_parton", "m_{W2} Parton", mBins, wmmin, wmmax); listOfTH1.push_back(hW2mparton);

  // ww
    TH1F *hWWdeltaPhireco = new TH1F("WW_#Delta#phi_reco", "#Delta#phi_{WW} Reco", phiBins, -TMath::Pi(), +TMath::Pi()); listOfTH1.push_back(hWWdeltaPhireco);
    TH1F *hWWdeltaEtareco = new TH1F("WW_#Delta#eta_reco", "#Delta#eta_{WW} Reco", etaBins, wetamin, wetamax);listOfTH1.push_back(hWWdeltaEtareco);
    TH1F *hWWdeltaRreco = new TH1F("WW_#DeltaR_reco", "#DeltaR_{WW} Reco", RBins, wRmin, wRmax); listOfTH1.push_back(hWWdeltaRreco);

    TH1F *hWWdeltaPhiparticle = new TH1F("WW_#Delta#phi_particle", "#Delta#phi_{WW} Particle", phiBins, -TMath::Pi(), +TMath::Pi()); listOfTH1.push_back(hWWdeltaPhiparticle);
    TH1F *hWWdeltaEtaparticle = new TH1F("WW_#Delta#eta_particle", "#Delta#eta_{WW} Particle", etaBins, wetamin, wetamax);listOfTH1.push_back(hWWdeltaEtaparticle);
    TH1F *hWWdeltaRparticle = new TH1F("WW_#DeltaR_particle", "#DeltaR_{WW} Particle", RBins, wRmin, wRmax); listOfTH1.push_back(hWWdeltaRparticle);

    TH1F *hWWdeltaPhiparton = new TH1F("WW_#Delta#phi_parton", "#Delta#phi_{WW} Parton", phiBins, -TMath::Pi(), +TMath::Pi()); listOfTH1.push_back(hWWdeltaPhiparton);
    TH1F *hWWdeltaEtaparton = new TH1F("WW_#Delta#eta_parton", "#Delta#eta_{WW} Parton", etaBins, wetamin, wetamax);listOfTH1.push_back(hWWdeltaEtaparton);
    TH1F *hWWdeltaRparton = new TH1F("WW_#DeltaR_parton", "#DeltaR_{WW} Parton", RBins, wRmin, wRmax); listOfTH1.push_back(hWWdeltaRparton);

  // leptons

    TH1F *lead_e_pt_reco = new TH1F("lead_e_pT_reco", "p^{T}_{e1}_reco", pTBins, lpTmin, lpTmax); listOfTH1.push_back(lead_e_pt_reco);
    TH1F *lead_e_eta_reco = new TH1F("lead_e_#eta_reco", "#eta_{e1}_reco", etaBins, letamin, letamax); listOfTH1.push_back(lead_e_eta_reco);
    TH1F *lead_e_phi_reco = new TH1F("lead_e_#phi_reco", "#phi_{e1}_reco", phiBins, -TMath::Pi(), TMath::Pi()); listOfTH1.push_back(lead_e_phi_reco);
    TH1F *sublead_e_pt_reco = new TH1F("sublead_e_pT_reco", "p^{T}_{e2}_reco", pTBins, lpTmin, lpTmax); listOfTH1.push_back(sublead_e_pt_reco);
    TH1F *sublead_e_eta_reco = new TH1F("sublead_e_#eta_reco", "#eta_{e2}_reco", etaBins, letamin, letamax); listOfTH1.push_back(sublead_e_eta_reco);
    TH1F *sublead_e_phi_reco = new TH1F("sublead_e_#phi_reco", "#phi_{e2}_reco", phiBins, -TMath::Pi(), TMath::Pi()); listOfTH1.push_back(sublead_e_phi_reco);
    TH1F *goodE_size_reco = new TH1F("goodE_size_reco", "size", 5, 0, 5); listOfTH1.push_back(goodE_size_reco);

    TH1F *lead_mu_pt_reco = new TH1F("lead_mu_pT_reco", "p^{T}_{mu1}_reco", pTBins, lpTmin, lpTmax); listOfTH1.push_back(lead_mu_pt_reco);
    TH1F *lead_mu_eta_reco = new TH1F("lead_mu_#eta_reco", "#eta_{mu1}_reco", etaBins, letamin, letamax); listOfTH1.push_back(lead_mu_eta_reco);
    TH1F *lead_mu_phi_reco = new TH1F("lead_mu_#phi_reco", "#phi_{mu1}_reco", phiBins, -TMath::Pi(), TMath::Pi()); listOfTH1.push_back(lead_mu_phi_reco);
    TH1F *sublead_mu_pt_reco = new TH1F("sublead_mu_pT_reco", "p^{T}_{mu2}_reco", pTBins, lpTmin, lpTmax); listOfTH1.push_back(sublead_mu_pt_reco);
    TH1F *sublead_mu_eta_reco = new TH1F("sublead_mu_#eta_reco", "#eta_{mu2}_reco", etaBins, letamin, letamax); listOfTH1.push_back(sublead_mu_eta_reco);
    TH1F *sublead_mu_phi_reco = new TH1F("sublead_mu_#phi_reco", "#phi_{mu2}_reco", phiBins, -TMath::Pi(), TMath::Pi()); listOfTH1.push_back(sublead_mu_phi_reco);
    TH1F *goodMu_size_reco = new TH1F("goodMu_size_reco", "size", 5, 0, 5); listOfTH1.push_back(goodMu_size_reco);

    TH1F *lead_e_pt_particle = new TH1F("lead_e_pT_particle", "p^{T}_{e1}_particle", pTBins, lpTmin, lpTmax); listOfTH1.push_back(lead_e_pt_particle);
    TH1F *lead_e_eta_particle = new TH1F("lead_e_#eta_particle", "#eta_{e1}_particle", etaBins, letamin, letamax); listOfTH1.push_back(lead_e_eta_particle);
    TH1F *lead_e_phi_particle = new TH1F("lead_e_#phi_particle", "#phi_{e1}_particle", phiBins, -TMath::Pi(), TMath::Pi()); listOfTH1.push_back(lead_e_phi_particle);
    TH1F *sublead_e_pt_particle = new TH1F("sublead_e_pT_particle", "p^{T}_{e2}_particle", pTBins, lpTmin, lpTmax); listOfTH1.push_back(sublead_e_pt_particle);
    TH1F *sublead_e_eta_particle = new TH1F("sublead_e_#eta_particle", "#eta_{e2}_particle", etaBins, letamin, letamax); listOfTH1.push_back(sublead_e_eta_particle);
    TH1F *sublead_e_phi_particle = new TH1F("sublead_e_#phi_particle", "#phi_{e2}_particle", phiBins, -TMath::Pi(), TMath::Pi()); listOfTH1.push_back(sublead_e_phi_particle);
    TH1F *goodE_size_particle = new TH1F("goodE_size_particle", "size", 5, 0, 5); listOfTH1.push_back(goodE_size_particle);

    TH1F *lead_mu_pt_particle = new TH1F("lead_mu_pT_particle", "p^{T}_{mu1}_particle", pTBins, lpTmin, lpTmax); listOfTH1.push_back(lead_mu_pt_particle);
    TH1F *lead_mu_eta_particle = new TH1F("lead_mu_#eta_particle", "#eta_{mu1}_particle", etaBins, letamin, letamax); listOfTH1.push_back(lead_mu_eta_particle);
    TH1F *lead_mu_phi_particle = new TH1F("lead_mu_#phi_particle", "#phi_{mu1}_particle", phiBins, -TMath::Pi(), TMath::Pi()); listOfTH1.push_back(lead_mu_phi_particle);
    TH1F *sublead_mu_pt_particle = new TH1F("sublead_mu_pT_particle", "p^{T}_{mu2}_particle", pTBins, lpTmin, lpTmax); listOfTH1.push_back(sublead_mu_pt_particle);
    TH1F *sublead_mu_eta_particle = new TH1F("sublead_mu_#eta_particle", "#eta_{mu2}_particle", etaBins, letamin, letamax); listOfTH1.push_back(sublead_mu_eta_particle);
    TH1F *sublead_mu_phi_particle = new TH1F("sublead_mu_#phi_particle", "#phi_{mu2}_particle", phiBins, -TMath::Pi(), TMath::Pi()); listOfTH1.push_back(sublead_mu_phi_particle);
    TH1F *goodMu_size_particle = new TH1F("goodMu_size_particle", "size", 5, 0, 5); listOfTH1.push_back(goodMu_size_particle);

    TH1F *goodE_size_parton = new TH1F("goodE_size_parton", "size", 5, 0, 5); listOfTH1.push_back(goodE_size_parton);

    TH1F *goodMu_size_parton = new TH1F("goodMu_size_parton", "size", 5, 0, 5); listOfTH1.push_back(goodMu_size_parton);

// 2D - parton(1) particle(2) reco(3)

  // paired
    TH2F *hPJsize23Comp = new TH2F("PAIReD_jet_size_comp_23", "size", 5, 0, 5, 5, 0, 5); listOfTH2.push_back(hPJsize23Comp);
    TH2F *hPJBsize23Comp = new TH2F("PAIReD_b_jet_size_comp_23", "size", 5, 0, 5, 5, 0, 5); listOfTH2.push_back(hPJBsize23Comp);

  // higgs
    TH2F *hHpT12Comp = new TH2F("H_pT_comp_12", "p_{T}^{hbb} Parton vs. Particle", pTBins, hpTmin, hpTmax, pTBins, hpTmin, hpTmax); listOfTH2.push_back(hHpT12Comp);
    TH2F *hHpT23Comp = new TH2F("H_pT_comp_23", "p_{T}^{hbb} Particle vs. Reco", pTBins, hpTmin, hpTmax, pTBins, hpTmin, hpTmax); listOfTH2.push_back(hHpT23Comp);
    TH2F *hHpT13Comp = new TH2F("H_pT_comp_13", "p_{T}^{hbb} Parton vs. Reco", pTBins, hpTmin, hpTmax, pTBins, hpTmin, hpTmax); listOfTH2.push_back(hHpT13Comp);
    TH2F *hHm12Comp = new TH2F("H_m_comp_12", "m_{hbb} Parton vs. Particle", mBins, hmmin, hmmax, mBins, hmmin, hmmax); listOfTH2.push_back(hHm12Comp);
    TH2F *hHm23Comp = new TH2F("H_m_comp_23", "m_{hbb} Particle vs. Reco", mBins, hmmin, hmmax, mBins, hmmin, hmmax); listOfTH2.push_back(hHm23Comp);      
    TH2F *hHm13Comp = new TH2F("H_m_comp_13", "m_{hbb} Parton vs. Reco", mBins, hmmin, hmmax, mBins, hmmin, hmmax); listOfTH2.push_back(hHm13Comp);
    TH2F*hbbdeltaPhi12Comp = new TH2F("bb_#Delta#phi_comp_12", "#Delta#phi_{bb} Parton vs. Particle", phiBins, -TMath::Pi(),+TMath::Pi(), phiBins, -TMath::Pi(),+TMath::Pi()); listOfTH2.push_back(hbbdeltaPhi12Comp);
    TH2F*hbbdeltaPhi23Comp = new TH2F("bb_#Delta#phi_comp_23", "#Delta#phi_{bb} Particle vs. Reco", phiBins, -TMath::Pi(),+TMath::Pi(), phiBins, -TMath::Pi(),+TMath::Pi()); listOfTH2.push_back(hbbdeltaPhi23Comp);
    TH2F*hbbdeltaPhi13Comp = new TH2F("bb_#Delta#phi_comp_13", "#Delta#phi_{bb} Parton vs. Reco", phiBins, -TMath::Pi(),+TMath::Pi(), phiBins, -TMath::Pi(),+TMath::Pi()); listOfTH2.push_back(hbbdeltaPhi13Comp); 
    TH2F*hbbdeltaEta12Comp = new TH2F("bb_#Delta#eta_comp_12", "#Delta#eta_{bb} Parton vs. Particle", etaBins, hetamin, hetamax, etaBins, hetamin, hetamax); listOfTH2.push_back(hbbdeltaEta12Comp);
    TH2F*hbbdeltaEta23Comp = new TH2F("bb_#Delta#eta_comp_23", "#Delta#eta_{bb} Particle vs. Reco", etaBins, hetamin, hetamax, etaBins, hetamin, hetamax); listOfTH2.push_back(hbbdeltaEta23Comp);
    TH2F*hbbdeltaEta13Comp = new TH2F("bb_#Delta#eta_comp_13", "#Delta#eta_{bb} Parton vs. Reco", etaBins, hetamin, hetamax, etaBins, hetamin, hetamax); listOfTH2.push_back(hbbdeltaEta13Comp);

  // vbfj
    TH2F*hjjpT12Comp = new TH2F("jj_pT_comp_12", "p_{T}^{jj} Parton vs. Particle", pTBins, jpTmin, jpTmax, pTBins, jpTmin, jpTmax); listOfTH2.push_back(hjjpT12Comp);
    TH2F*hjjpT23Comp = new TH2F("jj_pT_comp_23", "p_{T}^{jj} Particle vs. Reco", pTBins, jpTmin, jpTmax, pTBins, jpTmin, jpTmax); listOfTH2.push_back(hjjpT23Comp);
    TH2F*hjjpT13Comp = new TH2F("jj_pT_comp_13", "p_{T}^{jj} Parton vs. Reco", pTBins, jpTmin, jpTmax, pTBins, jpTmin, jpTmax); listOfTH2.push_back(hjjpT13Comp);
    TH2F*hjjdeltaPhi12Comp = new TH2F("jj_#Delta#phi_comp_12", "#Delta#phi_{jj} Parton vs. Particle", phiBins, -TMath::Pi(),+TMath::Pi(), phiBins, -TMath::Pi(),+TMath::Pi()); listOfTH2.push_back(hjjdeltaPhi12Comp);
    TH2F*hjjdeltaPhi23Comp = new TH2F("jj_#Delta#phi_comp_23", "#Delta#phi_{jj} Particle vs. Reco", phiBins, -TMath::Pi(),+TMath::Pi(), phiBins, -TMath::Pi(),+TMath::Pi()); listOfTH2.push_back(hjjdeltaPhi23Comp);
    TH2F*hjjdeltaPhi13Comp = new TH2F("jj_#Delta#phi_comp_13", "#Delta#phi_{jj} Parton vs. Reco", phiBins, -TMath::Pi(),+TMath::Pi(), phiBins, -TMath::Pi(),+TMath::Pi()); listOfTH2.push_back(hjjdeltaPhi13Comp);

  // leptons
    TH2F*hl1pT12Comp = new TH2F("l1_pT_comp_12", "p^{T}_{l1}  Parton vs. Particle", pTBins, lpTmin, lpTmax, pTBins, lpTmin, lpTmax); listOfTH2.push_back(hl1pT12Comp);
    TH2F*hl1pT23Comp = new TH2F("l1_pT_comp_23", "p^{T}_{l1} Particle vs. Reco", pTBins, lpTmin, lpTmax, pTBins, lpTmin, lpTmax); listOfTH2.push_back(hl1pT23Comp);
    TH2F*hl1pT13Comp = new TH2F("l1_pT_comp_13", "p^{T}_{l1} Parton vs. Reco", pTBins, lpTmin, lpTmax, pTBins, lpTmin, lpTmax); listOfTH2.push_back(hl1pT13Comp);
    TH2F*hl2pT12Comp = new TH2F("l2_pT_comp_12", "p^{T}_{l2}  Parton vs. Particle", pTBins, lpTmin, lpTmax, pTBins, lpTmin, lpTmax); listOfTH2.push_back(hl2pT12Comp);
    TH2F*hl2pT23Comp = new TH2F("l2_pT_comp_23", "p^{T}_{l2} Particle vs. Reco", pTBins, lpTmin, lpTmax, pTBins, lpTmin, lpTmax); listOfTH2.push_back(hl2pT23Comp);
    TH2F*hl2pT13Comp = new TH2F("l2_pT_comp_13", "p^{T}_{l2} Parton vs. Reco", pTBins, lpTmin, lpTmax, pTBins, lpTmin, lpTmax); listOfTH2.push_back(hl2pT13Comp);

  // z comps
    TH2F*hz1pT12Comp = new TH2F("z1_pT_comp_12", "p^{T}_{z1}  Parton vs. Particle", pTBins, zpTmin, zpTmax, pTBins, zpTmin, zpTmax); listOfTH2.push_back(hz1pT12Comp);
    TH2F*hz1pT23Comp = new TH2F("z1_pT_comp_23", "p^{T}_{z1} Particle vs. Reco", pTBins, zpTmin, zpTmax, pTBins, zpTmin, zpTmax); listOfTH2.push_back(hz1pT23Comp);
    TH2F*hz1pT13Comp = new TH2F("z1_pT_comp_13", "p^{T}_{z1} Parton vs. Reco", pTBins, zpTmin, zpTmax, pTBins, zpTmin, zpTmax); listOfTH2.push_back(hz1pT13Comp);
    TH2F*hz2pT12Comp = new TH2F("z2_pT_comp_12", "p^{T}_{z2}  Parton vs. Particle", pTBins, zpTmin, zpTmax, pTBins, zpTmin, zpTmax); listOfTH2.push_back(hz2pT12Comp);
    TH2F*hz2pT23Comp = new TH2F("z2_pT_comp_23", "p^{T}_{z2} Particle vs. Reco", pTBins, zpTmin, zpTmax, pTBins, zpTmin, zpTmax); listOfTH2.push_back(hz2pT23Comp);
    TH2F*hz2pT13Comp = new TH2F("z2_pT_comp_13", "p^{T}_{z2} Parton vs. Reco", pTBins, zpTmin, zpTmax, pTBins, zpTmin, zpTmax); listOfTH2.push_back(hz2pT13Comp);
    TH2F*hz1m12Comp = new TH2F("z1_m_comp_12", "m_{z1}  Parton vs. Particle", mBins, zmmin, zmmax, mBins, zmmin, zmmax); listOfTH2.push_back(hz1m12Comp);
    TH2F*hz1m23Comp = new TH2F("z1_m_comp_23", "m_{z1} Particle vs. Reco", mBins, zmmin, zmmax, mBins, zmmin, zmmax); listOfTH2.push_back(hz1m23Comp);      
    TH2F*hz1m13Comp = new TH2F("z1_m_comp_13", "m_{z1} Parton vs. Reco", mBins, zmmin, zmmax, mBins, zmmin, zmmax); listOfTH2.push_back(hz1m13Comp);
    TH2F*hz2m12Comp = new TH2F("z2_m_comp_12", "m_{z2}  Parton vs. Particle", mBins, zmmin, zmmax, mBins, zmmin, zmmax); listOfTH2.push_back(hz2m12Comp);
    TH2F*hz2m23Comp = new TH2F("z2_m_comp_23", "m_{z2} Particle vs. Reco", mBins, zmmin, zmmax, mBins, zmmin, zmmax); listOfTH2.push_back(hz2m23Comp);
    TH2F*hz2m13Comp = new TH2F("z2_m_comp_13", "m_{z2} Parton vs. Reco", mBins, zmmin, zmmax, mBins, zmmin, zmmax); listOfTH2.push_back(hz2m13Comp);

  // w comps
    TH2F*hw1pT12Comp = new TH2F("w1_pT_comp_12", "p^{T}_{w1}  Parton vs. Particle", pTBins, wpTmin, wpTmax, pTBins, wpTmin, wpTmax); listOfTH2.push_back(hw1pT12Comp);
    TH2F*hw1pT23Comp = new TH2F("w1_pT_comp_23", "p^{T}_{w1} Particle vs. Reco", pTBins, wpTmin, wpTmax, pTBins, wpTmin, wpTmax); listOfTH2.push_back(hw1pT23Comp);
    TH2F*hw1pT13Comp = new TH2F("w1_pT_comp_13", "p^{T}_{w1} Parton vs. Reco", pTBins, wpTmin, wpTmax, pTBins, wpTmin, wpTmax); listOfTH2.push_back(hw1pT13Comp);
    TH2F*hw2pT12Comp = new TH2F("w2_pT_comp_12", "p^{T}_{w2}  Parton vs. Particle", pTBins, wpTmin, wpTmax, pTBins, wpTmin, wpTmax); listOfTH2.push_back(hw2pT12Comp);
    TH2F*hw2pT23Comp = new TH2F("w2_pT_comp_23", "p^{T}_{w2} Particle vs. Reco", pTBins, wpTmin, wpTmax, pTBins, wpTmin, wpTmax); listOfTH2.push_back(hw2pT23Comp);
    TH2F*hw2pT13Comp = new TH2F("w2_pT_comp_13", "p^{T}_{w2} Parton vs. Reco", pTBins, wpTmin, wpTmax, pTBins, wpTmin, wpTmax); listOfTH2.push_back(hw2pT13Comp);
    TH2F*hw1m12Comp = new TH2F("w1_m_comp_12", "m_{w1}  Parton vs. Particle", mBins, wmmin, wmmax, mBins, wmmin, wmmax); listOfTH2.push_back(hw1m12Comp);
    TH2F*hw1m23Comp = new TH2F("w1_m_comp_23", "m_{w1} Particle vs. Reco", mBins, wmmin, wmmax, mBins, wmmin, wmmax); listOfTH2.push_back(hw1m23Comp);      
    TH2F*hw1m13Comp = new TH2F("w1_m_comp_13", "m_{w1} Parton vs. Reco", mBins, wmmin, wmmax, mBins, wmmin, wmmax); listOfTH2.push_back(hw1m13Comp);
    TH2F*hw2m12Comp = new TH2F("w2_m_comp_12", "m_{w2}  Parton vs. Particle", mBins, wmmin, wmmax, mBins, wmmin, wmmax); listOfTH2.push_back(hw2m12Comp);
    TH2F*hw2m23Comp = new TH2F("w2_m_comp_23", "m_{w2} Particle vs. Reco", mBins, wmmin, wmmax, mBins, wmmin, wmmax); listOfTH2.push_back(hw2m23Comp);
    TH2F*hw2m13Comp = new TH2F("w2_m_comp_13", "m_{w2} Parton vs. Reco", mBins, wmmin, wmmax, mBins, wmmin, wmmax); listOfTH2.push_back(hw2m13Comp);


// 2d mixed
      
  // higgs dphi vs. vbfj dphi
    TH2F*hbbjjdeltaPhicompreco = new TH2F("bb_jj_#Delta#phi_comp_reco", "bb_jj_#Delta#phi_comp_reco", phiBins, -TMath::Pi(), +TMath::Pi(), phiBins, -TMath::Pi(), +TMath::Pi()); listOfTH2.push_back(hbbjjdeltaPhicompreco);
    TH2F*hbbjjdeltaPhicompparticle = new TH2F("bb_jj_#Delta#phi_comp_particle", "bb_jj_#Delta#phi_comp_particle", phiBins, -TMath::Pi(), +TMath::Pi(), phiBins, -TMath::Pi(), +TMath::Pi()); listOfTH2.push_back(hbbjjdeltaPhicompparticle);
    TH2F*hbbjjdeltaPhicompparton = new TH2F("bb_jj_#Delta#phi_comp_parton", "bb_jj_#Delta#phi_comp_parton", phiBins, -TMath::Pi(), +TMath::Pi(), phiBins, -TMath::Pi(), +TMath::Pi()); listOfTH2.push_back(hbbjjdeltaPhicompparton);

  // higgs deta vs. vbfj deta
    TH2F*hbbdeltaEtajjdeltaPhicompreco = new TH2F("bb_#Delta#eta_jj_#Delta#phi_comp_reco", "bb_#Delta#eta_jj_#Delta#phi_comp_reco", etaBins, hetamin, hetamax, phiBins, -TMath::Pi(), +TMath::Pi()); listOfTH2.push_back(hbbdeltaEtajjdeltaPhicompreco);
    TH2F*hbbdeltaEtajjdeltaPhicompparticle = new TH2F("bb_#Delta#eta_jj_#Delta#phi_comp_particle", "bb_#Delta#eta_jj_#Delta#phi_comp_particle", etaBins, hetamin, hetamax, phiBins, -TMath::Pi(), +TMath::Pi()); listOfTH2.push_back(hbbdeltaEtajjdeltaPhicompparticle);
    TH2F*hbbdeltaEtajjdeltaPhicompparton = new TH2F("bb_#Delta#eta_jj_#Delta#phi_comp_parton", "bb_#Delta#eta_jj_#Delta#phi_comp_parton", etaBins, hetamin, hetamax, phiBins, -TMath::Pi(), +TMath::Pi()); listOfTH2.push_back(hbbdeltaEtajjdeltaPhicompparton);

  // higgs pT vs. lepton dphi 
    TH2F*hHpTl1l2deltaPhicompreco = new TH2F("h_pT_l1l2_delta#phi_comp_reco", "h_pT_l1l2_delta#phi_comp_reco", pTBins, hpTmin, hpTmax, phiBins, -TMath::Pi(), +TMath::Pi()); listOfTH2.push_back(hHpTl1l2deltaPhicompreco);
    TH2F*hHpTl3l4deltaPhicompreco = new TH2F("h_pT_l3l4_delta#phi_comp_reco", "h_pT_l1l2_delta#phi_comp_reco", pTBins, hpTmin, hpTmax, phiBins, -TMath::Pi(), +TMath::Pi()); listOfTH2.push_back(hHpTl3l4deltaPhicompreco);
    TH2F*hHpTl1l2deltaPhicompparticle = new TH2F("h_pT_l1l2_delta#phi_comp_particle", "h_pT_l1l2_delta#phi_comp_particle", pTBins, hpTmin, hpTmax, phiBins, -TMath::Pi(), +TMath::Pi()); listOfTH2.push_back(hHpTl1l2deltaPhicompparticle);
    TH2F*hHpTl3l4deltaPhicompparticle = new TH2F("h_pT_l3l4_delta#phi_comp_particle", "h_pT_l1l2_delta#phi_comp_particle", pTBins, hpTmin, hpTmax, phiBins, -TMath::Pi(), +TMath::Pi()); listOfTH2.push_back(hHpTl3l4deltaPhicompparticle);
    TH2F*hHpTl1l2deltaPhicompparton = new TH2F("h_pT_l1l2_delta#phi_comp_parton", "h_pT_l1l2_delta#phi_comp_parton", pTBins, hpTmin, hpTmax, phiBins, -TMath::Pi(), +TMath::Pi()); listOfTH2.push_back(hHpTl1l2deltaPhicompparton);
    TH2F*hHpTl3l4deltaPhicompparton = new TH2F("h_pT_l3l4_delta#phi_comp_parton", "h_pT_l1l2_delta#phi_comp_parton", pTBins, hpTmin, hpTmax, phiBins, -TMath::Pi(), +TMath::Pi()); listOfTH2.push_back(hHpTl3l4deltaPhicompparton);

  // higgs pT vs. z pT
    TH2F*hHz1pTcompreco = new TH2F("h_z1_pT_comp_reco", "h_z1_pT_comp_reco", pTBins, hpTmin, hpTmax, pTBins, zpTmin, zpTmax); listOfTH2.push_back(hHz1pTcompreco);
    TH2F*hHz2pTcompreco = new TH2F("h_z2_pT_comp_reco", "h_z2_pT_comp_reco", pTBins, hpTmin, hpTmax, pTBins, zpTmin, zpTmax); listOfTH2.push_back(hHz2pTcompreco);
    TH2F*hHz1pTcompparticle = new TH2F("h_z1_pT_comp_particle", "h_z1_pT_comp_particle", pTBins, hpTmin, hpTmax, pTBins, zpTmin, zpTmax); listOfTH2.push_back(hHz1pTcompparticle);
    TH2F*hHz2pTcompparticle = new TH2F("h_z2_pT_comp_particle", "h_z2_pT_comp_particle", pTBins, hpTmin, hpTmax, pTBins, zpTmin, zpTmax); listOfTH2.push_back(hHz2pTcompparticle);
    TH2F*hHz1pTcompparton = new TH2F("h_z1_pT_comp_parton", "h_z1_pT_comp_parton", pTBins, hpTmin, hpTmax, pTBins, zpTmin, zpTmax); listOfTH2.push_back(hHz1pTcompparton);
    TH2F*hHz2pTcompparton = new TH2F("h_z2_pT_comp_parton", "h_z2_pT_comp_parton", pTBins, hpTmin, hpTmax, pTBins, zpTmin, zpTmax); listOfTH2.push_back(hHz2pTcompparton);

  // higgs pT vs. zz pT      
    TH2F*hHzzpTcompreco = new TH2F("h_zz_pT_comp_reco", "h_zz_pT_comp_reco", pTBins, hpTmin, 2*hpTmax, pTBins, zpTmin, 2*zpTmax); listOfTH2.push_back(hHzzpTcompreco);
    TH2F*hHzzpTcompparticle = new TH2F("h_zz_pT_comp_particle", "h_zz_pT_comp_particle", pTBins, hpTmin, 2*hpTmax, pTBins, zpTmin, 2*zpTmax); listOfTH2.push_back(hHzzpTcompparticle);
    TH2F*hHzzpTcompparton = new TH2F("h_zz_pT_comp_parton", "h_zz_pT_comp_parton", pTBins, hpTmin, 2*hpTmax, pTBins, zpTmin, 2*zpTmax); listOfTH2.push_back(hHzzpTcompparton);

  // higgs pT vs. zz dphi
    TH2F*hHpTzzdeltaPhicompreco = new TH2F("h_pT_zz_delta#phi_comp_reco", "h_pT_zz_delta#phi_comp_reco", pTBins, hpTmin, hpTmax, phiBins, -TMath::Pi(), +TMath::Pi()); listOfTH2.push_back(hHpTzzdeltaPhicompreco);
    TH2F*hHpTzzdeltaPhicompparticle = new TH2F("h_pT_zz_delta#phi_comp_particle", "h_pT_zz_delta#phi_comp_particle", pTBins, hpTmin, hpTmax, phiBins, -TMath::Pi(), +TMath::Pi()); listOfTH2.push_back(hHpTzzdeltaPhicompparticle);
    TH2F*hHpTzzdeltaPhicompparton = new TH2F("h_pT_zz_delta#phi_comp_parton", "h_pT_zz_delta#phi_comp_parton", pTBins, hpTmin, hpTmax, phiBins, -TMath::Pi(), +TMath::Pi()); listOfTH2.push_back(hHpTzzdeltaPhicompparton);

  // higgs pT vs. zz deta    
    TH2F*hHpTzzdeltaEtacompreco = new TH2F("h_pT_zz_delta#eta_comp_reco", "h_pT_zz_delta#eta_comp_reco", pTBins, hpTmin, hpTmax, etaBins, zetamin, zetamax); listOfTH2.push_back(hHpTzzdeltaEtacompreco);
    TH2F*hHpTzzdeltaEtacompparticle = new TH2F("h_pT_zz_delta#eta_comp_particle", "h_pT_zz_delta#eta_comp_particle", pTBins, hpTmin, hpTmax, etaBins, zetamin, zetamax); listOfTH2.push_back(hHpTzzdeltaEtacompparticle);
    TH2F*hHpTzzdeltaEtacompparton = new TH2F("h_pT_zz_delta#eta_comp_parton", "h_pT_zz_delta#eta_comp_parton", pTBins, hpTmin, hpTmax, etaBins, zetamin, zetamax); listOfTH2.push_back(hHpTzzdeltaEtacompparton);

  // higgs dphi vs. zz dphi
    TH2F*hbbzzdeltaPhicompreco = new TH2F("bb_zz_#Delta#phi_comp_reco", "bb_zz_#Delta#phi_comp_reco", phiBins, -TMath::Pi(), +TMath::Pi(), phiBins, -TMath::Pi(), +TMath::Pi()); listOfTH2.push_back(hbbzzdeltaPhicompreco);
    TH2F*hbbzzdeltaPhicompparticle = new TH2F("bb_zz_#Delta#phi_comp_particle", "bb_zz_#Delta#phi_comp_particle", phiBins, -TMath::Pi(), +TMath::Pi(), phiBins, -TMath::Pi(), +TMath::Pi()); listOfTH2.push_back(hbbzzdeltaPhicompparticle);
    TH2F*hbbzzdeltaPhicompparton = new TH2F("bb_zz_#Delta#phi_comp_parton", "bb_zz_#Delta#phi_comp_parton", phiBins, -TMath::Pi(), +TMath::Pi(), phiBins, -TMath::Pi(), +TMath::Pi()); listOfTH2.push_back(hbbzzdeltaPhicompparton);

  // higgs deta vs. zz deta  
    TH2F*hbbzzdeltaEtacompreco = new TH2F("bb_zz_delta#eta_comp_reco", "bb_zz_delta#eta_comp_reco", etaBins, hetamin, hetamax, etaBins, zetamin, zetamax); listOfTH2.push_back(hbbzzdeltaEtacompreco);
    TH2F*hbbzzdeltaEtacompparticle = new TH2F("bb_zz_delta#eta_comp_particle", "bb_zz_delta#eta_comp_particle", etaBins, hetamin, hetamax, etaBins, zetamin, zetamax); listOfTH2.push_back(hbbzzdeltaEtacompparticle);
    TH2F*hbbzzdeltaEtacompparton = new TH2F("bb_zz_delta#eta_comp_parton", "bb_zz_delta#eta_comp_parton", etaBins, hetamin, hetamax, etaBins, zetamin, zetamax); listOfTH2.push_back(hbbzzdeltaEtacompparton);


// misc

  // leptonic event types
    TH1F *recoET = new TH1F("reco_event_type", "ET", 5, -1, 4); listOfTH1.push_back(recoET);
    TH1F *particleET = new TH1F("particle_event_type", "ET", 5, -1, 4); listOfTH1.push_back(particleET);
    TH1F *partonET = new TH1F("parton_event_type", "ET", 5, -1, 4); listOfTH1.push_back(partonET);
  
  // w - reco by event type
    TH1F *hllpTET0reco = new TH1F("ll_ET0_pT_reco", "p^{T}_ET0_{ll}_reco", pTBins, wpTmin, wpTmax); listOfTH1.push_back(hllpTET0reco);
    TH1F *hllmET0reco = new TH1F("ll_ET0_m_reco", "m_{ll}_ET0_reco", mBins, wmmin, wmmax); listOfTH1.push_back(hllmET0reco);
    TH1F *hllpTET1reco = new TH1F("ll_ET1_pT_reco", "p^{T}_ET1_{ll}_reco", pTBins, wpTmin, wpTmax); listOfTH1.push_back(hllpTET1reco);
    TH1F *hllmET1reco = new TH1F("ll_ET1_m_reco", "m_{ll}_ET1_reco", mBins, wmmin, wmmax); listOfTH1.push_back(hllmET1reco);
    TH1F *hllpTET2reco = new TH1F("ll_ET2_pT_reco", "p^{T}_ET2_{ll}_reco", pTBins, wpTmin, wpTmax); listOfTH1.push_back(hllpTET2reco);
    TH1F *hllmET2reco = new TH1F("ll_ET2_m_reco", "m_{ll}_ET2_reco", mBins, wmmin, wmmax); listOfTH1.push_back(hllmET2reco);
    TH1F *hllpTET3reco = new TH1F("ll_ET3_pT_reco", "p^{T}_ET3_{ll}_reco", pTBins, wpTmin, wpTmax); listOfTH1.push_back(hllpTET3reco);
    TH1F *hllmET3reco = new TH1F("ll_ET3_m_reco", "m_{ll}_ET3_reco", mBins, wmmin, wmmax); listOfTH1.push_back(hllmET3reco);

  // w - particle by event type
    TH1F *hllpTET0particle = new TH1F("ll_ET0_pT_particle", "p^{T}_ET0_{ll}_particle", pTBins, wpTmin, wpTmax); listOfTH1.push_back(hllpTET0particle);
    TH1F *hllmET0particle = new TH1F("ll_ET0_m_particle", "m_{ll}_ET0_particle", mBins, wmmin, wmmax); listOfTH1.push_back(hllmET0particle);
    TH1F *hllpTET1particle = new TH1F("ll_ET1_pT_particle", "p^{T}_ET1_{ll}_particle", pTBins, wpTmin, wpTmax); listOfTH1.push_back(hllpTET1particle);
    TH1F *hllmET1particle = new TH1F("ll_ET1_m_particle", "m_{ll}_ET1_particle", mBins, wmmin, wmmax); listOfTH1.push_back(hllmET1particle);
    TH1F *hllpTET2particle = new TH1F("ll_ET2_pT_particle", "p^{T}_ET2_{ll}_particle", pTBins, wpTmin, wpTmax); listOfTH1.push_back(hllpTET2particle);
    TH1F *hllmET2particle = new TH1F("ll_ET2_m_particle", "m_{ll}_ET2_particle", mBins, wmmin, wmmax); listOfTH1.push_back(hllmET2particle);
    TH1F *hllpTET3particle = new TH1F("ll_ET3_pT_particle", "p^{T}_ET3_{ll}_particle", pTBins, wpTmin, wpTmax); listOfTH1.push_back(hllpTET3particle);
    TH1F *hllmET3particle = new TH1F("ll_ET3_m_particle", "m_{ll}_ET3_particle", mBins, wmmin, wmmax); listOfTH1.push_back(hllmET3particle);

  // mll [0,15]
    TH1F *hmll_0_15_reco = new TH1F("ll_m_reco_0_15", "m_{ll}_reco_0_15", 75, 0, 15); listOfTH1.push_back(hmll_0_15_reco);
    TH1F *hllm_0_15_particle = new TH1F("ll_m_particle_0_15", "m_{ll}_particle_0_15", 75, 0, 15); listOfTH1.push_back(hllm_0_15_particle);

  // ?
    vector <TH1F*> lepPT_partonV;
    for(int i=0; i<4; i++){
      lepPT_partonV.push_back(new TH1F(Form("lepPT_partonV_%d",i),"",50,0,2e2));
        listOfTH1.push_back(lepPT_partonV.at(i));
    }



//------------------------------------------------------------------------------------------------------------------------------------------------------------
// DEF QUANTITIES
//------------------------------------------------------------------------------------------------------------------------------------------------------------


  double  nPassed=0;
  double totWeightedEntries=0;
  int  nPassedRaw=0;

  double Lumi=1000*300; // 1000 to convert from pb to fb, 300 for end of run 3
  cout << "LUMI: "<< Lumi << endl;

// TLorentzVectors

  // higgs
  TLorentzVector b1_reco, b1_particle,  b1_parton;
  TLorentzVector b2_reco, b2_particle,  b2_parton;
  TLorentzVector h_reco, h_parton, h_particle;

  // vbfj
  TLorentzVector j1_reco, j1_particle,  j1_parton;
  TLorentzVector j2_reco, j2_particle,  j2_parton;

  // z
  TLorentzVector z1_reco, z1_particle,  z1_parton;
  TLorentzVector z2_reco, z2_particle,  z2_parton;

  // w
  TLorentzVector w1_reco, w1_particle,  w1_parton;
  TLorentzVector w2_reco, w2_particle,  w2_parton;
  TLorentzVector met;

  // leptons
  TLorentzVector l1_reco, l1_particle,  l1_parton;
  TLorentzVector l2_reco, l2_particle,  l2_parton;
  TLorentzVector l3_reco, l3_particle,  l3_parton;
  TLorentzVector l4_reco, l4_particle,  l4_parton;

  TLorentzVector e1_reco, e1_particle,  e1_parton;
  TLorentzVector e2_reco, e2_particle,  e2_parton;
  TLorentzVector e3_reco, e3_particle,  e3_parton;
  TLorentzVector e4_reco, e4_particle,  e4_parton;

  TLorentzVector m1_reco, m1_particle,  m1_parton;
  TLorentzVector m2_reco, m2_particle,  m2_parton;
  TLorentzVector m3_reco, m3_particle,  m3_parton;
  TLorentzVector m4_reco, m4_particle,  m4_parton;

// kinematic quantities
/*
  // higgs
  double bbdeltaPhireco = 9999;
  double bbdeltaEtareco = 9999;
  double bbdeltaRreco = -9999;

  double bbdeltaPhiparticle = -9999;
  double bbdeltaEtaparticle = -9999;
  double bbdeltaRparticle= -9999;

  double bbdeltaPhiparton = -9999;
  double bbdeltaEtaparton = -9999;
  double bbdeltaRparton= -9999;

  // vbfj
  double jjdeltaPhireco =  -9999;
  double jjdeltaEtareco= -9999;
  double jjdeltaRreco = -9999;

  double jjdeltaPhiparticle =  -9999;
  double jjdeltaEtaparticle =   -9999;
  double jjdeltaRparticle   =  -9999;

  double jjdeltaPhiparton = -9999;
  double jjdeltaEtaparton = -9999;

  double jet1_pt = -9999;
  double jet2_pt = -9999;
  double jet1_energy = -9999;
  double jet2_energy = -9999;
  double pairedbbmass = -9999;

  // leps
  double l1l2deltaPhireco= -99999;
  double l3l4deltaPhireco=-99999;
  double l1l2deltaEtareco= -99999;
  double l3l4deltaEtareco=-99999; 
  double l1l2deltaRreco=-99999;
  double l3l4deltaRreco=-99999;

  double l1l2deltaPhiBoostreco=-99999;
  double l3l4deltaPhiBoostreco=-99999;
  double l1l2deltaEtaBoostreco=-99999;
  double l3l4deltaEtaBoostreco=-99999;

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

  double l1l2CScosThetareco=-99999;
  double l3l4CScosThetareco=-99999;

  double l1l2deltaPhiparticle=-9999;
  double l3l4deltaPhiparticle=-9999;
  double l1l2deltaEtaparticle=-9999;
  double l3l4deltaEtaparticle=-9999;
  double l1l2deltaRparticle=-9999;
  double l3l4deltaRparticle=-9999;

  double l1l2deltaPhiBoostparticle=-99999;
  double l3l4deltaPhiBoostparticle=-99999;
  double l1l2deltaEtaBoostparticle=-99999;
  double l3l4deltaEtaBoostparticle=-99999;

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

  double l1l2deltaPhiparton=-9999;
  double l3l4deltaPhiparton=-9999;
  double l1l2deltaEtaparton=-9999;
  double l3l4deltaEtaparton=-9999;
  double l1l2deltaRparton=-9999;
  double l3l4deltaRparton=-9999;

  // add dphi deta boost

  double l1cosThetaparton=-9999;
  double l2cosThetaparton=-9999;
  double l3cosThetaparton=-9999;
  double l4cosThetaparton=-9999;
  double fourlcosThetaparton=-9999;

  double l1cosThetaBoostparton=-9999;
  double l2cosThetaBoostparton=-9999;
  double l3cosThetaBoostparton=-9999;
  double l4cosThetaBoostparton=-9999;
  double fourlcosThetaBoostparton=-9999;

  double l1l2CScosThetaparton=-9999;
  double l3l4CScosThetaparton=-9999;

  // z
  double zzdeltaPhireco= -99999;
  double zzdeltaEtareco= -99999;
  double zzdeltaRreco=-99999;

  double zzdeltaPhiparticle=-9999;
  double zzdeltaEtaparticle=-9999;
  double zzdeltaRparticle=-9999;
      
  double zzdeltaPhiparton=-9999;
  double zzdeltaEtaparton=-9999;
  double zzdeltaRparton=-9999;

  // w
  double wwdeltaPhireco= -99999;
  double wwdeltaEtareco= -99999;
  double wwdeltaRreco=-99999;

  double wwdeltaPhiparticle=-9999;
  double wwdeltaEtaparticle=-9999;
  double wwdeltaRparticle=-9999;
      
  double wwdeltaPhiparton=-9999;
  double wwdeltaEtaparton=-9999;
  double wwdeltaRparton=-9999;

  // lep charge
  int q1_reco=0;
  int q2_reco=0;
  int q3_reco=0;
  int q4_reco=0;

  int q1_particle=0;
  int q2_particle=0;
  int q3_particle=0;
  int q4_particle=0;

  int q1_parton = -9999;
  int q2_parton = -9999;
  int q3_parton = -9999;
  int q4_parton = -9999;
*/


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


#ifdef MDEBUG
    numberOfEntries=1000;
#endif 

  for(Int_t entry = 0; entry < numberOfEntries; ++entry) {

    treeReader->ReadEntry(entry);
    HepMCEvent *event = (HepMCEvent*) branchEvent -> At(0);
    Float_t weight = event->Weight*Lumi*cross_section*numberOfEntries/(numEntries*totalWeight);
    Float_t test_weight = event->Weight*cross_section*numberOfEntries/(numEntries*totalWeight);
    hWeight -> Fill(event->Weight, test_weight);


  //------------------------------------------------------------------------------------------------------------------------------------------------------------
  // RECO - HIGGS 
  //------------------------------------------------------------------------------------------------------------------------------------------------------------

    int switchVal_reco = 0;
    if(enableCutReco["initial - reco"]) increaseCount(cutFlowMap_reco,"initial - reco",weight);
    

    bool foundHiggs_reco = false;

    vector <int> goodJetIndex=GoodJetIndices(branchJet);

    if(enableCutReco["jet pT > 20 - reco"]) {
      if(switchVal_reco == 0 && goodJetIndex.size() > 0) increaseCount(cutFlowMap_reco,"jet pT > 20 - reco", weight);
      else switchVal_reco = 1;
    }
    
    std::vector<std::pair< std::map<TString, float>, std::map<TString, std::vector<float>>>>  pairedJet=paired::PAIReDjointEvent(branchGenParticle,branchPFCand,branchJet,0.4,false,false,true,1.0,false);
    //cout<<"PAIRED lables bb "<<pairedJet.first["label_bb"]<<" cc "<<pairedJet.first["label_cc"]<<" ll "<<pairedJet.first["label_ll"]<<" indices 1: "<<pairedJet.first["jet1_index"]<<" 2: "<<pairedJet.first["jet1_index"]<<endl;
    std::vector<std::pair< std::map<TString, float>, std::map<TString, std::vector<float>>>>  pairedJetB;

    if(enableCutReco["1 PAIReD jet - reco"]) {
      if(switchVal_reco == 0 && pairedJet.size() > 0) increaseCount(cutFlowMap_reco,"1 PAIReD jet - reco",weight);
      else switchVal_reco = 1;
    }

    for(int i=0; i<(int)pairedJet.size(); i++){
      std::pair< std::map<TString, float>, std::map<TString, std::vector<float>>> thisPaired=pairedJet.at(i);
      if( thisPaired.first["isbtagged"] > 0) pairedJetB.push_back(thisPaired);
    }

    vector <int> btagIndex;
    int pairedJetSize_reco = pairedJet.size();
    int pairedBJetSize_reco = pairedJetB.size();

    if(enableCutReco["1 bb PAIReD jet - reco"]) {
        if(switchVal_reco == 0 && pairedJetB.size()>0){
          increaseCount(cutFlowMap_reco,"1 bb PAIReD jet - reco",weight);
          foundHiggs_reco = true;
        } else switchVal_reco = 1;
    }

    std::map<TString, float> paired_jet;

    if (switchVal_reco == 0 && foundHiggs_reco){

      paired_jet = pairedJetB.at(0).first;

      btagIndex.push_back(paired_jet["jet1_index"]);
      btagIndex.push_back(paired_jet["jet2_index"]);

      b1_reco.SetPtEtaPhiM(paired_jet["jet1_pt"],paired_jet["jet1_eta"],paired_jet["jet1_phi"],paired_jet["jet1_mass"]);
      b2_reco.SetPtEtaPhiM(paired_jet["jet2_pt"],paired_jet["jet2_eta"],paired_jet["jet2_phi"],paired_jet["jet2_mass"]);

      h_reco = b1_reco + b2_reco; // dijet

      double bbdeltaPhireco = deltaPhi(b1_reco, b2_reco)
      double bbdeltaEtareco = deltaEta(b1_reco, b2_reco)
      double bbdeltaRreco = deltaR(b1_reco, b2_reco)

    }


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
      }*/
    
    
    
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


  //------------------------------------------------------------------------------------------------------------------------------------------------------------
  // RECO - VBF JETS
  //------------------------------------------------------------------------------------------------------------------------------------------------------------

    bool foundVBF_reco = False;

    vector <int> nonHiggsJet;
    vector<pair<int,int>> vbfJetIndex;
    vector<vector <int>> vbfJetIndexComb;
    vector<pair<int,int>> vbfJetIndex_dEta;
    int vbfJetIndexCandidate = -1;

    // check that they do not belong to higgs + sort by pT
    for(int i=0; i<(int)branchJet->GetEntries(); i++) {
      if( goodJetIndex[i] == paired_jet["jet1_index"] || goodJetIndex[i] == paired_jet["jet2_index"] ) continue;
    nonHiggsJet.push_back(goodJetIndex[i]);
    SortByPtIndices(nonHiggsJet,branchJet);
    }

    // check that there are at least two + make combinations
    if(enableCutReco["2 VBF jet - reco"]) {
      if(switchVal_reco==0 && nonHiggsJet.size() > 1 ) {
	      increaseCount(cutFlowMap_reco,"2 VBF jet - reco",weight);
	      vbfJetIndexComb=combinationsNoRepetitionAndOrderDoesNotMatter(2,nonHiggsJet);
      }
      else switchVal_reco=1;
    }

    // convert vector of vector of ints to vector of pairs of ints
    vbfJetIndex=GetvbfJetIndex(vbfJetIndexComb);

    // loop and take those w dEta > 2.5
    for (int i=0; i<(int)vbfJetIndex.size(); i++) {
      if( fabs((((Jet*)branchJet->At(vbfJetIndex[i].first))->Eta - ((Jet*)branchJet->At(vbfJetIndex[i].second))->Eta)) <= 2.5 ) continue; 
      else { 
	      vbfJetIndex_dEta.push_back(vbfJetIndex.at(i));
	      vbfJetIndexCandidate = i;
      }
    }

    if(enableCutReco["2.5 deltaEta VBF jet - reco"]) {
      if(switchVal_reco==0 && vbfJetIndex_dEta.size()>0) {
        increaseCount(cutFlowMap_reco,"2.5 deltaEta VBF jet - reco",weight);
        foundVBF_reco = True;
      } else switchVal_reco=1;
    }

    // sort them again by eta for leading/subleading
    SortByEtaIndices(vbfJetIndex_dEta,branchJet); 

    Jet *jet1 =nullptr;
    Jet *jet2 =nullptr;
  
    if(switchVal_reco==0 && foundVBF_reco) {

      jet1 = (Jet*) branchJet->At(vbfJetIndex_dEta[0].first);
      jet2 = (Jet*) branchJet->At(vbfJetIndex_dEta[0].second);
      j1_reco=jet1->P4();
      j2_reco=jet2->P4();

      double jjdeltaPhireco = deltaPhi(j1_reco, j2_reco)
      double jjdeltaEtareco = deltaEta(j1_reco, j2_reco)
      double jjdeltaRreco = deltaR(j1_reco, j2_reco)

    }

  //------------------------------------------------------------------------------------------------------------------------------------------------------------
  // RECO - LEPTONS
  //------------------------------------------------------------------------------------------------------------------------------------------------------------

  int thisRecoEventType=-1;

  double e_pT_min = 15.0;
  double e_eta_max = 2.5;
  double mu_pT_min = 15.0;
  double mu_eta_max = 2.5;
    
  // get e+e- mu+mu-
  vector <int> goodE_min_reco_indices = get_good_reco_lepton_indices(branchElectron, e_pT_min, e_eta_max, analysis, "electron", -1);
  vector <int> goodE_plus_reco_indices = get_good_reco_lepton_indices(branchElectron, e_pT_min, e_eta_max, analysis, "electron", 1);
  vector <int> goodMu_min_reco_indices = get_good_reco_lepton_indices(branchMuon, mu_pT_min, mu_eta_max, analysis, "muon", -1);
  vector <int> goodMu_plus_reco_indices = get_good_reco_lepton_indices(branchMuon, mu_pT_min, mu_eta_max, analysis, "muon", 1);

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
  ConcatenateIndices(goodE_min_reco_indices, goodE_reco_indices);
  ConcatenateIndices(goodE_plus_reco_indices, goodE_reco_indices);
  ConcatenateIndices(goodMu_min_reco_indices, goodMu_reco_indices);
  ConcatenateIndices(goodMu_plus_reco_indices, goodMu_reco_indices);

  if(enableCutReco["lep pT & eta cut - reco"]){
      if(switchVal_reco == 0 && ((goodE_reco_indices.size() + goodMu_reco_indices.size()) > 0)) increaseCount(cutFlowMap_reco,"lep pT & eta cut - reco",weight);
      else switchVal_reco = 1;
  }

  // re-sort
  sort_by_pT(goodE_reco_indices, branchElectron, "Electron"); 
  sort_by_pT(goodMu_reco_indices, branchMuon, "Muon");

  goodE_size_reco->Fill(goodE_reco_indices.size(),weight);
  goodMu_size_reco->Fill(goodMu_reco_indices.size(),weight);

  // lep details
  if (goodE_reco_indices.size()>0) e1_reco = ((Electron *) branchElectron->At(goodE_reco_indices[0]))->P4();
  if (goodE_reco_indices.size()>1) e2_reco = ((Electron *) branchElectron->At(goodE_reco_indices[1]))->P4();
  if (goodMu_reco_indices.size()>0) m1_reco = ((Muon *) branchMuon->At(goodMu_reco_indices[0]))->P4();
  if (goodMu_reco_indices.size()>1) m2_reco = ((Muon *) branchMuon->At(goodMu_reco_indices[1]))->P4();

  lead_e_pt_reco->Fill(e1_reco.Pt(),weight);
  lead_e_eta_reco->Fill(e1_reco.Eta(),weight);
  lead_e_phi_reco->Fill(e1_reco.Phi(),weight);
  sublead_e_pt_reco->Fill(e2_reco.Pt(),weight);
  sublead_e_eta_reco->Fill(e2_reco.Eta(),weight);
  sublead_e_phi_reco->Fill(e2_reco.Phi(),weight);

  lead_mu_pt_reco->Fill(m1_reco.Pt(),weight);
  lead_mu_eta_reco->Fill(m1_reco.Eta(),weight);
  lead_mu_phi_reco->Fill(m1_reco.Phi(),weight);
  sublead_mu_pt_reco->Fill(m2_reco.Pt(),weight);
  sublead_mu_eta_reco->Fill(m2_reco.Eta(),weight);
  sublead_mu_phi_reco->Fill(m2_reco.Phi(),weight);
  
// V details

  vector<pair<int,pair<int,int>>> ZRecoPairIndices;
  vector <int> wleps;

  if(analysis == "HZZJJ"){
  
    // form pairs for each flavour
    vector< pair<int,int>> elecZRecoPairIndices=GetelecRecoPairIndices(branchElectron,goodE_reco_indices); 
    vector< pair<int,int>> muZRecoPairIndices=GetmuRecoPairIndices(branchMuon,goodMu_reco_indices);
    
    //increaseCount(cutFlowMap_reco,"at least two lep pairs",weight); 
    
    ZRecoPairIndices=GetRecoPairIndices(elecZRecoPairIndices,muZRecoPairIndices,branchElectron,branchMuon); // 0 for electron 1 for muon

    if(enableCutReco["OSSF - reco"]){
      if (switchVal_reco==0 && ZRecoPairIndices.size()>=2) increaseCount(cutFlowMap_reco,"OSSF - reco",weight);
      else switchVal_reco=1;
    }

    if( switchVal_reco==0 && ZRecoPairIndices.size()>=2){
      if( ZRecoPairIndices[0].first == 1 && ZRecoPairIndices[1].first == 1) thisRecoEventType=0;
      else if( ZRecoPairIndices[0].first == 0 && ZRecoPairIndices[1].first == 0) thisRecoEventType=1;
      else if( ZRecoPairIndices[0].first == 1 && ZRecoPairIndices[1].first == 0) thisRecoEventType=2;
      else if( ZRecoPairIndices[0].first == 0 && ZRecoPairIndices[1].first == 1) thisRecoEventType=3;
    }

    recoET->Fill(thisRecoEventType, weight);

    getRecoZLeps(thisRecoEventType, ZRecoPairIndices, branchElectron, branchMuon, l1_reco, l2_reco, l3_reco, l4_reco, q1_reco, q2_reco, q3_reco, q4_reco);

    if( switchVal_reco == 0 && thisRecoEventType != -1 && ZRecoPairIndices.size() >= 2 ){
      z1_reco=l1_reco + l2_reco;
      z2_reco=l3_reco + l4_reco;
      zzdeltaPhireco = deltaPhi(z1_reco, z2_reco)
      zzdeltaEtareco = deltaEta(z1_reco, z2_reco)
      zzdeltaRreco = deltaR(z1_reco, z2_reco)
    }

    // WW 

   } else if(analysis == "HWWJJ") {
      
      // FOR OSOF SWITCH mu mu / e e EVENT TYPE TO -1

      if((goodE_reco_indices.size() + goodMu_reco_indices.size()) >= 2 ){
        if (goodMu_min_reco_indices.size() > 0 && goodMu_plus_reco_indices.size() > 0) { // case mu- mu+
            thisRecoEventType = 0;
            wleps.push_back(goodMu_min_reco_indices[0]);
            wleps.push_back(goodMu_plus_reco_indices[0]);
            sort2_by_pT(wleps, branchElectron, branchMuon, "Muon", "Muon");
            l1_reco = ((Muon *) branchMuon->At(wleps[0]))->P4();
            l2_reco = ((Muon *) branchMuon->At(wleps[1]))->P4();
        } else if (goodE_min_reco_indices.size() > 0 && goodE_plus_reco_indices.size() > 0) {  // case e- e+ 
            thisRecoEventType = 1; 
            wleps.push_back(goodE_min_reco_indices[0]);
            wleps.push_back(goodE_plus_reco_indices[0]);
            sort2_by_pT(wleps, branchElectron, branchMuon, "Electron", "Electron");
            l1_reco = ((Electron *) branchElectron->At(wleps[0]))->P4();
            l2_reco = ((Electron *) branchElectron->At(wleps[1]))->P4();
        } else if (goodMu_min_reco_indices.size() > 0 && goodE_plus_reco_indices.size() > 0) { // case mu- e+
            thisRecoEventType = 2;
            wleps.push_back(goodMu_min_reco_indices[0]);
            wleps.push_back(goodE_plus_reco_indices[0]);
            sort2_by_pT(wleps, branchElectron, branchMuon, "Muon", "Electron");
            l1_reco = ((Muon *) branchMuon->At(wleps[0]))->P4();
            l2_reco = ((Electron *) branchElectron->At(wleps[1]))->P4();
        } else if (goodE_min_reco_indices.size() > 0 && goodMu_plus_reco_indices.size() > 0) {  // case e- mu+
            thisRecoEventType = 3; 
            wleps.push_back(goodE_min_reco_indices[0]);
            wleps.push_back(goodMu_plus_reco_indices[0]);
            sort2_by_pT(wleps, branchElectron, branchMuon, "Electron", "Muon");
            l1_reco = ((Electron *) branchElectron->At(wleps[0]))->P4();
            l2_reco = ((Muon *) branchMuon->At(wleps[1]))->P4();
        }
    }

      recoET->Fill(thisRecoEventType);

      if(enableCutReco["OSOF - reco"]){
        if (switchVal_reco==0 && thisRecoEventType != -1) increaseCount(cutFlowMap_reco, "OSOF - reco", weight);
        else switchVal_reco=1;
      }

      hmll_0_15_reco->Fill((l1_reco+l2_reco).M(),weight);

     if(enableCutReco["mll > 10 - reco"]) {
          if(switchVal_reco == 0 && (l1_reco+l2_reco).M() >= 10) increaseCount(cutFlowMap_reco,"mll > 10 - reco",weight);
          else switchVal_reco = 1;
      }

      met = ((MissingET*)branchMissingET->At(0))->P4();

      double w1mass = calculate_mT(l1_reco.Pt(), met.Pt(), l1_reco.Phi() - met.Phi());
      double w2mass = calculate_mT(l2_reco.Pt(), met.Pt(), l2_reco.Phi() - met.Phi());

      if( switchVal_reco == 0){
        w1_reco=l1_reco + met;
        w2_reco=l2_reco + met;
        double wwdeltaPhireco = deltaPhi(w1_reco, w2_reco)
        double wwdeltaEtareco = deltaEta(w1_reco, w2_reco)
        double wwdeltaRreco = deltaR(w1_reco, w2_reco)
      }

    }

  //------------------------------------------------------------------------------------------------------------------------------------------------------------
  // PARTICLE - HIGGS
  //------------------------------------------------------------------------------------------------------------------------------------------------------------

    int switchVal_particle = 0;
    bool foundBjetParticle = false;

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
        foundBjetParticle = true;
      } else switchVal_particle = 1;
    }

    std::map<TString, float> paired_jet_particle;

    if (switchVal_particle == 0 && pairedJetBParticle.size()>0){

      foundBjetParticle = true;
      paired_jet_particle = pairedJetBParticle.at(0).first;

      btagIndexParticle.push_back(paired_jet_particle["jet1_index"]);
      btagIndexParticle.push_back(paired_jet_particle["jet2_index"]);

      b1_particle.SetPtEtaPhiM(paired_jet_particle["jet1_pt"],paired_jet_particle["jet1_eta"],paired_jet_particle["jet1_phi"],paired_jet_particle["jet1_mass"]);
      b2_particle.SetPtEtaPhiM(paired_jet_particle["jet2_pt"],paired_jet_particle["jet2_eta"],paired_jet_particle["jet2_phi"],paired_jet_particle["jet2_mass"]);

      h_particle = b1_particle + b2_particle; // dijet

      double bbdeltaPhiparticle = deltaPhi(b1_particle, b2_particle)
      double bbdeltaEtarparticle = deltaEta(b1_particle, b2_particle)
      double bbdeltaRrparticle = deltaR(b1_particle, b2_particle)
 
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

      double jjdeltaPhiparticle = deltaPhi(j1_particle, j2_particle)
      double jjdeltaEtarparticle = deltaEta(j1_particle, j2_particle)
      double jjdeltaRrparticle = deltaR(j1_particle, j2_particle)
  
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
  ConcatenateIndices( goodE_min_particle_indices,goodE_particle_indices);
  ConcatenateIndices( goodE_plus_particle_indices,goodE_particle_indices);
  ConcatenateIndices( goodMu_min_particle_indices,goodMu_particle_indices);
  ConcatenateIndices( goodMu_plus_particle_indices,goodMu_particle_indices);

  vector<int> goodLep_particle_indices;
  ConcatenateIndices( goodE_particle_indices,goodLep_particle_indices);
  ConcatenateIndices( goodMu_particle_indices,goodLep_particle_indices);

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
      double zzdeltaPhiparticle = deltaPhi(z1_particle, z2_particle)
      double zzdeltaEtarparticle = deltaEta(z1_particle, z2_particle)
      double zzdeltaRrparticle = deltaR(z1_particle, z2_particle)
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
        double wwdeltaPhiparticle = deltaPhi(w1_particle, w2_particle)
        double wwdeltaEtarparticle = deltaEta(w1_particle, w2_particle)
        double wwdeltaRrparticle = deltaR(w1_particle, w2_particle)

      }

    }

    //------------------------------------------------------------------------------------------------------------------------------------------------------------
    // PARTON - HIGGS + JETS
    //------------------------------------------------------------------------------------------------------------------------------------------------------------
 
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
      double bbdeltaPhiparton = deltaPhi(b1_parton, b2_parton)
      double bbdeltaEtarparton = deltaEta(b1_parton, b2_parton)
      double bbdeltaRrparton = deltaR(b1_parton, b2_parton)
    }

    //------------------------------------------------------------------------------------------------------------------------------------------------------------
    // PARTON - LEPTONS
    //------------------------------------------------------------------------------------------------------------------------------------------------------------

    int thisPartonEventType=-1;

    vector <int> goodE_parton_indices  = GoodElectronPartonIndices(branchGenParticle, analysis);
    vector <int> goodMu_parton_indices = GoodMuonPartonIndices(branchGenParticle, analysis);
    goodE_size_parton->Fill(goodE_parton_indices.size(),weight);
    goodMu_size_parton->Fill(goodMu_parton_indices.size(),weight);

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
        double zzdeltaPhiparton = deltaPhi(z1_parton, z2_parton)
        double zzdeltaEtarparton = deltaEta(z1_parton, z2_parton)
        double zzdeltaRrparton = deltaR(z1_parton, z2_parton)
      }


    } if(analysis == "HWWJJ") { 

      WPartonIndices = GetWPartonIndices(branchGenParticle, analysis);

      if(WPartonIndices.size() > 1) foundWW = true;

      if(enableCutParton["WW parton"]){
        if(switchVal_parton == 0 && foundWW) increaseCount(cutFlowMap_parton,"WW parton",weight);
        else switchVal_parton = 1;
      }

      if(switchVal_parton == 0) getPartonWLeps(thisPartonEventType, WPartonIndices, branchGenParticle, w1_parton, w2_parton, l1_parton, l2_parton, q1_parton, q2_parton);

      partonET->Fill(thisPartonEventType,weight);

      if(foundWW) {
	
        double wwdeltaPhiparton = deltaPhi(w1_parton, w2_parton)
        double wwdeltaEtarparton = deltaEta(w1_parton, w2_parton)
        double wwdeltaRrparton = deltaR(w1_parton, w2_parton)

        lepPT_partonV.at(0)->Fill(l1_parton.Pt());
        lepPT_partonV.at(1)->Fill(l2_parton.Pt());
	
      }
    }

  //------------------------------------------------------------------------------------------------------------------------------------------------------------
  // PRINT CFT 
  //------------------------------------------------------------------------------------------------------------------------------------------------------------


  if( entry % 1000 == 0 ){
      cout<<"Processed "<<entry<< " / " <<numberOfEntries <<" "<< entry/ numberOfEntries *100 <<" %"<<endl;
      PrintCutFlow(cutFlowMap_reco,cutList_reco,"Reco");
      PrintCutFlow(cutFlowMap_particle,cutList_particle, "Particle");
      PrintCutFlow(cutFlowMap_parton,cutList_parton, "Parton");
    }


  //------------------------------------------------------------------------------------------------------------------------------------------------------------
  // FILL HISTOGRAMS - RECO
  //------------------------------------------------------------------------------------------------------------------------------------------------------------

  // 1D

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
      if(foundHiggs_parton){
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
        //hllpTreco->Fill((l1_reco+l2_reco).Pt(),weight);
        //hllmreco->Fill((l1_reco+l2_reco).M(),weight);
        hW1pTreco->Fill(w1_reco.Pt(),weight);
        hW1mreco->Fill(massTransverse(l1_reco, met),weight);
        hW2pTreco->Fill(w2_reco.Pt(),weight);
        hW2mreco->Fill(massTransverse(l2_reco, met),weight);
        hWWpTreco->Fill((w1_reco + w2_reco).Pt(),weight);
        hWWmreco->Fill((w1_reco + w2_reco).M(),weight);

        hWWdeltaPhireco->Fill(wwdeltaPhireco,weight); 
        hWWdeltaEtareco->Fill(wwdeltaEtareco, weight);
        hWWdeltaRreco -> Fill(wwdeltaRreco,weight);
      }
    }

    if(switchVal_particle==0){
      if(thisParticleEventType!=-1 && WParticlePairIndices.size()>=2){
        //hllpTparticle->Fill((l1_particle+l2_particle).Pt(),weight);
        //hllmparticle->Fill((l1_particle+l2_particle).M(),weight);
        hW1pTparticle->Fill(w1_particle.Pt(),weight);
        hW1mparticle->Fill(massTransverse(l1_particle, met),weight);
        hW2pTparticle->Fill(w2_particle.Pt(),weight);
        hW2mparticle->Fill(massTransverse(l2_particle, met),weight);
        hWWpTparticle->Fill((w1_particle + w2_sparticle).Pt(),weight);
        hWWmparticle->Fill((w1_particle + w2_particle).M(),weight);

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

  hPJsize23Comp -> Fill(pairedJetSize_particle, pairedJetSize_reco, weight);
  hPJBsize23Comp -> Fill(pairedBJetSize_particle, pairedBJetSize_reco, weight);

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


  //------------------------------------------------------------------------------------------------------------------------------------------------------------
  // END OF EVENT LOOP
  //------------------------------------------------------------------------------------------------------------------------------------------------------------


    nPassed+=weight;
    nPassedRaw++;

    cutVal_reco++; cutValW_reco+=weight;
    cutVal_particle++; cutValW_particle+=weight;
    cutVal_parton++; cutValW_parton+=weight;

  }
    

//------------------------------------------------------------------------------------------------------------------------------------------------------------
// FILL+PRINT CUTFLOW
//------------------------------------------------------------------------------------------------------------------------------------------------------------


  cout << " " << endl;
  cout << "RECO CUT FLOW" << endl;
    PrintCutFlow(cutFlowMap_reco,cutList_reco,  "Reco");
  cout << "PARTICLE CUT FLOW" << endl;
    PrintCutFlow(cutFlowMap_particle,cutList_particle, "Particle");
  cout << "PARTON CUT FLOW" << endl;
    PrintCutFlow(cutFlowMap_parton,cutList_parton, "Parton");
   
  for(std::vector<string>::iterator it=selType.begin(); it!=selType.end(); it++){
    FillCutFlow(cutFlowHists[(*it)],cutFlowEffs[(*it)],*cutFlowMapAll[(*it)],cutFlowMByType[(*it)], (*it));
  }
  
//------------------------------------------------------------------------------------------------------------------------------------------------------------
// WRITE HISTOGRAMS
//------------------------------------------------------------------------------------------------------------------------------------------------------------

  hists->cd();

  for(std::vector<TH1F*>::iterator h=listOfTH1.begin(); h!=listOfTH1.end(); h++){
    ( *h)->Write();
    delete (*h); 
  }

  for(std::vector<TH2F*>::iterator h=listOfTH2.begin(); h!=listOfTH2.end(); h++){
    ( *h)->Write();
    delete (*h); 
  }

  for(std::vector<TProfile*>::iterator h=listOfTProfiles.begin(); h!=listOfTProfiles.end(); h++){
    ( *h)->Write();
    delete (*h); 
  }
  
  //kappaLambda -> Write();
  //delete kappaLambda;
   
  hists->Close();

}

//------------------------------------------------------------------------------------------------------------------------------------------------------------
// MAIN
//------------------------------------------------------------------------------------------------------------------------------------------------------------


int main(int argc, char* argv[]){

  const char *inputFileName = argv[1];
  const char *outputFileName = argv[2];
  const char *process_name = argv[3];
  
  string analysisType="HZZJJ";

  // O: for ZZ H JJ, 1: (->H) ZZ  jj: 2: (->H) ZZ, 3: Hjj, 4: WW H JJ, 5: WW (->H) jj: 6: (->H) WW

  if( argc > 2 )  analysisType=string(argv[4]);

  cout << "RUNNING ANALYSIS: " << analysisType << endl;

  zAnalyzer(inputFileName, outputFileName, process_name, analysisType);

  return 0;
}
