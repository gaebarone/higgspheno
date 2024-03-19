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
  const double lpTmax = 500;
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

  // higgs - reco
  TH1F *hHpTreco = new TH1F("hbb_pT_reco", "p^{T}_{hbb}_reco", pTBins, hpTmin, hpTmax); 	listOfTH1.push_back(hHpTreco);
  TH1F *hHmreco = new TH1F("hbb_m_reco", "m_{hbb}_reco", mBins, hmmin, hmmax); listOfTH1.push_back(hHmreco);
  TH1F *hbbdeltaPhireco = new TH1F("bb_#Delta#phi_reco", "#Delta#phi_{bb}_reco", phiBins, -TMath::Pi(), +TMath::Pi()); listOfTH1.push_back(hbbdeltaPhireco);
  TH1F *hbbdeltaEtareco = new TH1F("bb_#Delta#eta_reco", "#Delta#eta_{bb}_reco", etaBins, hetamin, hetamax); listOfTH1.push_back(hbbdeltaEtareco);
  TH1F *hbbdeltaRreco = new TH1F("bb_#DeltaR_reco", "#DeltaR_{bb}_reco", RBins, hRmin, hRmax); listOfTH1.push_back(hbbdeltaRreco);

  // higgs - particle
  TH1F *hHpTparticle = new TH1F("hbb_pT_particle", "p^{T}_{hbb}_particle", pTBins, hpTmin, hpTmax); listOfTH1.push_back(hHpTparticle);
  TH1F *hHmparticle = new TH1F("hbb_m_particle", "m_{hbb}_particle", mBins, hmmin, hmmax); listOfTH1.push_back(hHmparticle);
  TH1F *hbbdeltaPhiparticle = new TH1F("bb_#Delta#phi_particle", "#Delta#phi_{bb}_particle", phiBins, -TMath::Pi(), +TMath::Pi()); listOfTH1.push_back(hbbdeltaPhiparticle);
  TH1F *hbbdeltaEtaparticle = new TH1F("bb_#Delta#eta_particle", "#Delta#eta_{bb}_particle", etaBins, hetamin, hetamax); listOfTH1.push_back(hbbdeltaEtaparticle);
  TH1F *hbbdeltaRparticle = new TH1F("bb_#DeltaR_particle", "#DeltaR_{bb}_particle", RBins, hRmin, hRmax); listOfTH1.push_back(hbbdeltaRparticle);

  // higgs - parton
  TH1F *hHpTparton = new TH1F("hbb_pT_parton", "p^{T}_{hbb}_parton", pTBins, hpTmin, hpTmax); listOfTH1.push_back(hHpTparton);
  TH1F *hHmparton = new TH1F("hbb_m_parton", "m_{hbb}_parton", mBins, hmmin,  hmmax); listOfTH1.push_back(hHmparton);
  TH1F *hbbdeltaPhiparton = new TH1F("bb_#Delta#phi_parton", "#Delta#phi_{bb}_parton", phiBins, -TMath::Pi(), +TMath::Pi()); listOfTH1.push_back(hbbdeltaPhiparton);
  TH1F *hbbdeltaEtaparton = new TH1F("bb_#Delta#eta_parton", "#Delta#eta_{bb}_parton", etaBins, hetamin, hetamax); listOfTH1.push_back(hbbdeltaEtaparton);
  TH1F *hbbdeltaRparton = new TH1F("bb_#DeltaR_parton", "#DeltaR_{bb}_parton", RBins, hRmin, hRmax); listOfTH1.push_back(hbbdeltaRparton);

  // vbfj - reco
  TH1F *hjjpTreco = new TH1F("jj_pT_reco", "p^{T}_{jj}_reco", pTBins, jpTmin, jpTmax); listOfTH1.push_back(hjjpTreco);
  TH1F *hjjdeltaPhireco = new TH1F("jj_#Delta#phi_reco", "#Delta#phi_{jj}_reco", phiBins, -TMath::Pi(), +TMath::Pi()); listOfTH1.push_back(hjjdeltaPhireco);
  TH1F *hjjdeltaEtareco = new TH1F("jj_#Delta#eta_reco", "#Delta#eta_{jj}_reco", etaBins, jetamin, jetamax); listOfTH1.push_back(hjjdeltaEtareco);
  TH1F *hjjdeltaRreco = new TH1F("jj_#DeltaR_reco", "#DeltaR_{jj}_reco", RBins, jRmin, jRmax); listOfTH1.push_back(hjjdeltaRreco);

  // vbfj - particle
  TH1F *hjjpTparticle = new TH1F("jj_pT_particle", "p^{T}_{jj}_particle", pTBins, jpTmin, jpTmax); listOfTH1.push_back(hjjpTparticle);
  TH1F *hjjdeltaPhiparticle = new TH1F("jj_#Delta#phi_particle", "#Delta#phi_{jj}_particle", phiBins, -TMath::Pi(), +TMath::Pi()); listOfTH1.push_back(hjjdeltaPhiparticle);
  TH1F *hjjdeltaEtaparticle = new TH1F("jj_#Delta#eta_particle", "#Delta#eta_{jj}_particle", etaBins, jetamin, jetamax); listOfTH1.push_back(hjjdeltaEtaparticle);
  TH1F *hjjdeltaRparticle = new TH1F("jj_#DeltaR_particle", "#DeltaR_{jj}_particle", RBins, jRmin, jRmax); listOfTH1.push_back(hjjdeltaRparticle);

  // z - reco
  TH1F *hz1pTreco = new TH1F("z1_pT_reco", "p^{T}_{z1}_reco", pTBins, zpTmin, zpTmax); listOfTH1.push_back(hz1pTreco);
  TH1F *hz2pTreco = new TH1F("z2_pT_reco", "p^{T}_{z2}_reco", pTBins, zpTmin, zpTmax); listOfTH1.push_back(hz2pTreco);
  TH1F *hz1mreco = new TH1F("z1_m_reco", "m_{z1}_reco", mBins, zmmin, zmmax); listOfTH1.push_back(hz1mreco);
  TH1F *hz2mreco = new TH1F("z2_m_reco", "m_{z2}_reco", mBins, zmmin, zmmax); listOfTH1.push_back(hz2mreco);
  TH1F *hzzdeltaPhireco = new TH1F("zz_#Delta#phi_reco", "#Delta#phi_{zz}_reco", phiBins, -TMath::Pi(), +TMath::Pi()); listOfTH1.push_back(hzzdeltaPhireco);
  TH1F *hzzdeltaEtareco = new TH1F("zz_#Delta#eta_reco", "#Delta#eta_{zz}_reco", etaBins, zetamin, zetamax);listOfTH1.push_back(hzzdeltaEtareco);
  TH1F *hzzdeltaRreco = new TH1F("zz_#DeltaR_reco", "#DeltaR_{zz}_reco", RBins, zRmin, zRmax); listOfTH1.push_back(hzzdeltaRreco);
  // cosTheta
  TH1F *hz1cosThetareco = new TH1F("z1_cos#theta_reco", "cos#theta_{z1}_reco", cosBins, -1, 1); listOfTH1.push_back(hz1cosThetareco);
  TH1F *hz2cosThetareco = new TH1F("z2_cos#theta_reco", "cos#theta_{z2}_reco", cosBins, -1, 1); listOfTH1.push_back(hz2cosThetareco);

  // z - particle
  TH1F *hz1pTparticle = new TH1F("z1_pT_particle", "p^{T}_{z1}_particle", pTBins, zpTmin, zpTmax);listOfTH1.push_back(hz1pTparticle);
  TH1F *hz2pTparticle = new TH1F("z2_pT_particle", "p^{T}_{z2}_particle", pTBins, zpTmin, zpTmax);listOfTH1.push_back(hz2pTparticle);
  TH1F *hz1mparticle = new TH1F("z1_m_particle", "m_{z1}_particle", mBins, zmmin, zmmax);listOfTH1.push_back(hz1mparticle);
  TH1F *hz2mparticle = new TH1F("z2_m_particle", "m_{z2}_particle", mBins, zmmin, zmmax);listOfTH1.push_back(hz2mparticle);
  TH1F *hzzdeltaPhiparticle = new TH1F("zz_#Delta#phi_particle", "#Delta#phi_{zz}_particle", phiBins, -TMath::Pi(), +TMath::Pi()); listOfTH1.push_back(hzzdeltaPhiparticle);
  TH1F *hzzdeltaEtaparticle = new TH1F("zz_#Delta#eta_particle", "#Delta#eta_{zz}_particle", etaBins, zetamin, zetamax); listOfTH1.push_back(hzzdeltaEtaparticle);
  TH1F *hzzdeltaRparticle = new TH1F("zz_#DeltaR_particle", "#DeltaR_{zz}_particle", RBins, zRmin, zRmax); listOfTH1.push_back(hzzdeltaRparticle);
  // cosTheta
  TH1F *hz1cosThetaparticle = new TH1F("z1_cos#theta_particle", "cos#theta_{z1}_particle", cosBins, -1, 1); listOfTH1.push_back(hz1cosThetaparticle);
  TH1F *hz2cosThetaparticle = new TH1F("z2_cos#theta_particle", "cos#theta_{z2}_particle", cosBins, -1, 1); listOfTH1.push_back(hz2cosThetaparticle);

  // z - parton
  TH1F *hz1pTparton = new TH1F("z1_pT_parton", "p^{T}_{z1}_parton", pTBins, zpTmin, zpTmax); listOfTH1.push_back(hz1pTparton);
  TH1F *hz2pTparton = new TH1F("z2_pT_parton", "p^{T}_{z2}_parton", pTBins, zpTmin, zpTmax); listOfTH1.push_back(hz2pTparton);
  TH1F *hz1mparton = new TH1F("z1_m_parton", "m_{z1}_parton", mBins, zmmin, zmmax); listOfTH1.push_back(hz1mparton);
  TH1F *hz2mparton = new TH1F("z2_m_parton", "m_{z2}_parton", mBins, zmmin, zmmax); listOfTH1.push_back(hz2mparton);
  TH1F *hzzdeltaPhiparton = new TH1F("zz_#Delta#phi_parton", "#Delta#phi_{zz}_parton", phiBins, -TMath::Pi(), +TMath::Pi()); listOfTH1.push_back(hzzdeltaPhiparton);
  TH1F *hzzdeltaEtaparton = new TH1F("zz_#Delta#eta_parton", "#Delta#eta_{zz}_parton", etaBins, zetamin, zetamax); listOfTH1.push_back(hzzdeltaEtaparton);
  TH1F *hzzdeltaRparton = new TH1F("zz_#DeltaR_parton", "#DeltaR_{zz}_parton", RBins, zRmin, zRmax); listOfTH1.push_back(hzzdeltaRparton);

  // w - reco
  TH1F *hllpTreco = new TH1F("ll_pT_reco", "p^{T}_{ll}_reco", pTBins, wpTmin, wpTmax); listOfTH1.push_back(hllpTreco);
  TH1F *hllmreco = new TH1F("ll_m_reco", "m_{ll}_reco", mBins, wmmin, wmmax); listOfTH1.push_back(hllmreco);

  TH1F *hllm_0_15_reco = new TH1F("ll_m_reco_0_15", "m_{ll}_reco_0_15", 75, 0, 15); listOfTH1.push_back(hllm_0_15_reco);

  TH1F *hw1pTreco = new TH1F("w1_pT_reco", "p^{T}_{w1}_reco", pTBins, wpTmin, wpTmax); listOfTH1.push_back(hw1pTreco);
  TH1F *hw1mreco = new TH1F("w1_m_reco", "m_{w1}_reco", mBins, wmmin, wmmax); listOfTH1.push_back(hw1mreco);
  TH1F *hw2pTreco = new TH1F("w2_pT_reco", "p^{T}_{w2}_reco", pTBins, wpTmin, wpTmax); listOfTH1.push_back(hw2pTreco);
  TH1F *hw2mreco = new TH1F("w2_m_reco", "m_{w2}_reco", mBins, wmmin, wmmax); listOfTH1.push_back(hw2mreco);
  TH1F *hwwpTreco = new TH1F("ww_pT_reco", "p^{T}_{ww}_reco", pTBins, wpTmin, wpTmax); listOfTH1.push_back(hwwpTreco);
  TH1F *hwwmreco = new TH1F("ww_m_reco", "m_{ww}_reco", mBins, wmmin, wmmax); listOfTH1.push_back(hwwmreco);
  TH1F *hwwdeltaPhireco = new TH1F("ww_#Delta#phi_reco", "#Delta#phi_{ww}_reco", phiBins, -TMath::Pi(), +TMath::Pi()); listOfTH1.push_back(hwwdeltaPhireco);
  TH1F *hwwdeltaEtareco = new TH1F("ww_#Delta#eta_reco", "#Delta#eta_{ww}_reco", etaBins, wetamin, wetamax);listOfTH1.push_back(hwwdeltaEtareco);
  TH1F *hwwdeltaRreco = new TH1F("ww_#DeltaR_reco", "#DeltaR_{ww}_reco", RBins, wRmin, wRmax); listOfTH1.push_back(hwwdeltaRreco);
  
  // w - reco by event type
  TH1F *hllpTET0reco = new TH1F("ll_ET0_pT_reco", "p^{T}_ET0_{ll}_reco", pTBins, wpTmin, wpTmax); listOfTH1.push_back(hllpTET0reco);
  TH1F *hllmET0reco = new TH1F("ll_ET0_m_reco", "m_{ll}_ET0_reco", mBins, wmmin, wmmax); listOfTH1.push_back(hllmET0reco);
  TH1F *hllpTET1reco = new TH1F("ll_ET1_pT_reco", "p^{T}_ET1_{ll}_reco", pTBins, wpTmin, wpTmax); listOfTH1.push_back(hllpTET1reco);
  TH1F *hllmET1reco = new TH1F("ll_ET1_m_reco", "m_{ll}_ET1_reco", mBins, wmmin, wmmax); listOfTH1.push_back(hllmET1reco);
  TH1F *hllpTET2reco = new TH1F("ll_ET2_pT_reco", "p^{T}_ET2_{ll}_reco", pTBins, wpTmin, wpTmax); listOfTH1.push_back(hllpTET2reco);
  TH1F *hllmET2reco = new TH1F("ll_ET2_m_reco", "m_{ll}_ET2_reco", mBins, wmmin, wmmax); listOfTH1.push_back(hllmET2reco);
  TH1F *hllpTET3reco = new TH1F("ll_ET3_pT_reco", "p^{T}_ET3_{ll}_reco", pTBins, wpTmin, wpTmax); listOfTH1.push_back(hllpTET3reco);
  TH1F *hllmET3reco = new TH1F("ll_ET3_m_reco", "m_{ll}_ET3_reco", mBins, wmmin, wmmax); listOfTH1.push_back(hllmET3reco);

  // w - particle
  TH1F *hllpTparticle = new TH1F("ll_pT_particle", "p^{T}_{ll}_particle", pTBins, wpTmin, wpTmax); listOfTH1.push_back(hllpTparticle);
  TH1F *hllmparticle = new TH1F("ll_m_particle", "m_{ll}_particle", mBins, wmmin, wmmax); listOfTH1.push_back(hllmparticle);

  TH1F *hllm_0_15_particle = new TH1F("ll_m_particle_0_15", "m_{ll}_particle_0_15", 75, 0, 15); listOfTH1.push_back(hllm_0_15_particle);

  TH1F *hw1pTparticle = new TH1F("w1_pT_particle", "p^{T}_{w1}_particle", pTBins, wpTmin, wpTmax); listOfTH1.push_back(hw1pTparticle);
  TH1F *hw1mparticle = new TH1F("w1_m_particle", "m_{w1}_particle", mBins, wmmin, wmmax); listOfTH1.push_back(hw1mparticle);
  TH1F *hw2pTparticle = new TH1F("w2_pT_particle", "p^{T}_{w2}_particle", pTBins, wpTmin, wpTmax); listOfTH1.push_back(hw2pTparticle);
  TH1F *hw2mparticle = new TH1F("w2_m_particle", "m_{w2}_particle", mBins, wmmin, wmmax); listOfTH1.push_back(hw2mparticle);
  TH1F *hwwpTparticle = new TH1F("ww_pT_particle", "p^{T}_{ww}_particle", pTBins, wpTmin, wpTmax); listOfTH1.push_back(hwwpTparticle);
  TH1F *hwwmparticle = new TH1F("ww_m_particle", "m_{ww}_particle", mBins, wmmin, wmmax); listOfTH1.push_back(hwwmparticle);
  TH1F *hwwdeltaPhiparticle = new TH1F("ww_#Delta#phi_particle", "#Delta#phi_{ww}_particle", phiBins, -TMath::Pi(), +TMath::Pi()); listOfTH1.push_back(hwwdeltaPhiparticle);
  TH1F *hwwdeltaEtaparticle = new TH1F("ww_#Delta#eta_particle", "#Delta#eta_{ww}_particle", etaBins, wetamin, wetamax); listOfTH1.push_back(hwwdeltaEtaparticle);
  TH1F *hwwdeltaRparticle = new TH1F("ww_#DeltaR_particle", "#DeltaR_{ww}_particle", RBins, wRmin, wRmax); listOfTH1.push_back(hwwdeltaRparticle);

  // w - particle by event type
  TH1F *hllpTET0particle = new TH1F("ll_ET0_pT_particle", "p^{T}_ET0_{ll}_particle", pTBins, wpTmin, wpTmax); listOfTH1.push_back(hllpTET0particle);
  TH1F *hllmET0particle = new TH1F("ll_ET0_m_particle", "m_{ll}_ET0_particle", mBins, wmmin, wmmax); listOfTH1.push_back(hllmET0particle);
  TH1F *hllpTET1particle = new TH1F("ll_ET1_pT_particle", "p^{T}_ET1_{ll}_particle", pTBins, wpTmin, wpTmax); listOfTH1.push_back(hllpTET1particle);
  TH1F *hllmET1particle = new TH1F("ll_ET1_m_particle", "m_{ll}_ET1_particle", mBins, wmmin, wmmax); listOfTH1.push_back(hllmET1particle);
  TH1F *hllpTET2particle = new TH1F("ll_ET2_pT_particle", "p^{T}_ET2_{ll}_particle", pTBins, wpTmin, wpTmax); listOfTH1.push_back(hllpTET2particle);
  TH1F *hllmET2particle = new TH1F("ll_ET2_m_particle", "m_{ll}_ET2_particle", mBins, wmmin, wmmax); listOfTH1.push_back(hllmET2particle);
  TH1F *hllpTET3particle = new TH1F("ll_ET3_pT_particle", "p^{T}_ET3_{ll}_particle", pTBins, wpTmin, wpTmax); listOfTH1.push_back(hllpTET3particle);
  TH1F *hllmET3particle = new TH1F("ll_ET3_m_particle", "m_{ll}_ET3_particle", mBins, wmmin, wmmax); listOfTH1.push_back(hllmET3particle);

  // w - parton
  TH1F *hw1pTparton = new TH1F("w1_pT_parton", "p^{T}_{w1}_parton", pTBins, wpTmin, wpTmax); listOfTH1.push_back(hw1pTparton);
  TH1F *hw2pTparton = new TH1F("w2_pT_parton", "p^{T}_{w2}_parton", pTBins, wpTmin, wpTmax); listOfTH1.push_back(hw2pTparton);
  TH1F *hw1mparton = new TH1F("w1_m_parton", "m_{w1}_parton", mBins, wmmin, wmmax); listOfTH1.push_back(hw1mparton);
  TH1F *hw2mparton = new TH1F("w2_m_parton", "m_{w2}_parton", mBins, wmmin, wmmax); listOfTH1.push_back(hw2mparton);
  TH1F *hwwdeltaPhiparton = new TH1F("ww_#Delta#phi_parton", "#Delta#phi_{ww}_parton", phiBins, -TMath::Pi(), +TMath::Pi()); listOfTH1.push_back(hwwdeltaPhiparton);
  TH1F *hwwdeltaEtaparton = new TH1F("ww_#Delta#eta_parton", "#Delta#eta_{ww}_parton", etaBins, wetamin, wetamax); listOfTH1.push_back(hwwdeltaEtaparton);
  TH1F *hwwdeltaRparton = new TH1F("ww_#DeltaR_parton", "#DeltaR_{ww}_parton", RBins, wRmin, wRmax); listOfTH1.push_back(hwwdeltaRparton);

  // lepton - reco
  TH1F *hl1l2deltaPhireco = new TH1F("l1l2_#Delta#phi_reco", "#Delta#phi_{l1l2}_reco", phiBins, -TMath::Pi(), +TMath::Pi()); listOfTH1.push_back(hl1l2deltaPhireco);
  TH1F *hl3l4deltaPhireco = new TH1F("l3l4_#Delta#phi_reco", "#Delta#phi_{l3l4}_reco", phiBins, -TMath::Pi(), +TMath::Pi()); listOfTH1.push_back(hl3l4deltaPhireco);
  TH1F *hl1l2deltaPhiBoostreco = new TH1F("l1l2_#Delta#phi_Boost_reco", "#Delta#phi_{l1l2}_reco", phiBins, -TMath::Pi(), +TMath::Pi()); listOfTH1.push_back(hl1l2deltaPhiBoostreco);
  TH1F *hl3l4deltaPhiBoostreco = new TH1F("l3l4_#Delta#phi_Boost_reco", "#Delta#phi_{l3l4}_reco", phiBins, -TMath::Pi(), +TMath::Pi()); listOfTH1.push_back(hl3l4deltaPhiBoostreco);
  TH1F *hl1l2deltaEtareco = new TH1F("l1l2_#Delta#eta_reco", "#Delta#eta_{l1l2}_reco", etaBins, letamin, letamax); listOfTH1.push_back(hl1l2deltaEtareco);
  TH1F *hl3l4deltaEtareco = new TH1F("l3l4_#Delta#eta_reco", "#Delta#eta_{l3l4}_reco", etaBins, letamin, letamax); listOfTH1.push_back(hl3l4deltaEtareco);
  TH1F *hl1l2deltaEtaBoostreco = new TH1F("l1l2_#Delta#eta_Boost_reco", "#Delta#eta_{l1l2}_Boost_reco", etaBins, letamin, letamax); listOfTH1.push_back(hl1l2deltaEtaBoostreco);
  TH1F *hl3l4deltaEtaBoostreco = new TH1F("l3l4_#Delta#eta_Boost_reco", "#Delta#eta_{l3l4}_Boost_reco", etaBins, letamin, letamax); listOfTH1.push_back(hl3l4deltaEtaBoostreco);
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

  // lepton - particle
  TH1F *hl1l2deltaPhiparticle = new TH1F("l1l2_#Delta#phi_particle", "#Delta#phi_{l1l2}_particle", phiBins, -TMath::Pi(), +TMath::Pi()); listOfTH1.push_back(hl1l2deltaPhiparticle);
  TH1F *hl3l4deltaPhiparticle = new TH1F("l3l4_#Delta#phi_particle", "#Delta#phi_{l3l4}_particle", phiBins, -TMath::Pi(), +TMath::Pi()); listOfTH1.push_back(hl3l4deltaPhiparticle);
  TH1F *hl1l2deltaPhiBoostparticle = new TH1F("l1l2_#Delta#phi_Boost_particle", "#Delta#phi_{l1l2}_particle", phiBins, -TMath::Pi(), +TMath::Pi()); listOfTH1.push_back(hl1l2deltaPhiBoostparticle);
  TH1F *hl3l4deltaPhiBoostparticle = new TH1F("l3l4_#Delta#phi_Boost_particle", "#Delta#phi_{l3l4}_particle", phiBins, -TMath::Pi(), +TMath::Pi()); listOfTH1.push_back(hl3l4deltaPhiBoostparticle);
  TH1F *hl1l2deltaEtaparticle = new TH1F("l1l2_#Delta#eta_particle", "#Delta#eta_{l1l2}_particle", etaBins, letamin, letamax); listOfTH1.push_back(hl1l2deltaEtaparticle);
  TH1F *hl3l4deltaEtaparticle = new TH1F("l3l4_#Delta#eta_particle", "#Delta#eta_{l3l4}_particle", etaBins, letamin, letamax); listOfTH1.push_back(hl3l4deltaEtaparticle);
  TH1F *hl1l2deltaEtaBoostparticle = new TH1F("l1l2_#Delta#eta_Boost_particle", "#Delta#eta_{l1l2}_Boost_particle", etaBins, letamin, letamax); listOfTH1.push_back(hl1l2deltaEtaBoostparticle);
  TH1F *hl3l4deltaEtaBoostparticle = new TH1F("l3l4_#Delta#eta_Boost_particle", "#Delta#eta_{l3l4}_Boost_particle", etaBins, letamin, letamax);listOfTH1.push_back(hl3l4deltaEtaBoostparticle);
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

  // lepton - parton
  TH1F *hl1l2deltaPhiparton = new TH1F("l1l2_#Delta#phi_parton", "#Delta#phi_{l1l2}_parton", phiBins, -TMath::Pi(), +TMath::Pi()); listOfTH1.push_back(hl1l2deltaPhiparton);
  TH1F *hl3l4deltaPhiparton = new TH1F("l3l4_#Delta#phi_parton", "#Delta#phi_{l3l4}_parton", phiBins, -TMath::Pi(), +TMath::Pi()); listOfTH1.push_back(hl3l4deltaPhiparton);
  TH1F *hl1l2deltaPhiBoostparton = new TH1F("l1l2_#Delta#phi_Boost_parton", "#Delta#phi_{l1l2}_parton", phiBins, -TMath::Pi(), +TMath::Pi()); listOfTH1.push_back(hl1l2deltaPhiBoostparton);
  TH1F *hl3l4deltaPhiBoostparton = new TH1F("l3l4_#Delta#phi_Boost_parton", "#Delta#phi_{l3l4}_parton", phiBins, -TMath::Pi(), +TMath::Pi()); listOfTH1.push_back(hl3l4deltaPhiBoostparton);
  TH1F *hl1l2deltaEtaparton = new TH1F("l1l2_#Delta#eta_parton", "#Delta#eta_{l1l2}_parton", etaBins, letamin, letamax); listOfTH1.push_back(hl1l2deltaEtaparton);
  TH1F *hl3l4deltaEtaparton = new TH1F("l3l4_#Delta#eta_parton", "#Delta#eta_{l3l4}_parton", etaBins, letamin, letamax); listOfTH1.push_back(hl3l4deltaEtaparton);
  TH1F *hl1l2deltaEtaBoostparton = new TH1F("l1l2_#Delta#eta_Boost_parton", "#Delta#eta_{l1l2}_Boost_parton", etaBins, letamin, letamax); listOfTH1.push_back(hl1l2deltaEtaBoostparton);
  TH1F *hl3l4deltaEtaBoostparton = new TH1F("l3l4_#Delta#eta_Boost_parton", "#Delta#eta_{l3l4}_Boost_parton", etaBins, letamin, letamax); listOfTH1.push_back(hl3l4deltaEtaBoostparton);
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


// 2D - parton(1) particle(2) reco(3)

  // paired
  TH2F *hPJsize23Comp = new TH2F("PAIReD_jet_size_comp_23", "size", 5, 0, 5, 5, 0, 5); listOfTH2.push_back(hPJsize23Comp);
  TH2F *hPJBsize23Comp = new TH2F("PAIReD_b_jet_size_comp_23", "size", 5, 0, 5, 5, 0, 5); listOfTH2.push_back(hPJBsize23Comp);

  // higgs
  TH2F *hHpT12Comp = new TH2F("H_pT_comp_12", "p_{T}^{hbb}", pTBins, hpTmin, hpTmax, pTBins, hpTmin, hpTmax); listOfTH2.push_back(hHpT12Comp);
  TH2F *hHpT23Comp = new TH2F("H_pT_comp_23", "p_{T}^{hbb}", pTBins, hpTmin, hpTmax, pTBins, hpTmin, hpTmax); listOfTH2.push_back(hHpT23Comp);
  TH2F *hHpT13Comp = new TH2F("H_pT_comp_13", "p_{T}^{hbb}", pTBins, hpTmin, hpTmax, pTBins, hpTmin, hpTmax); listOfTH2.push_back(hHpT13Comp);
  TH2F *hHm12Comp = new TH2F("H_m_comp_12", "m_{hbb}", mBins, hmmin, hmmax, mBins, hmmin, hmmax); listOfTH2.push_back(hHm12Comp);
  TH2F *hHm23Comp = new TH2F("H_m_comp_23", "m_{hbb}", mBins, hmmin, hmmax, mBins, hmmin, hmmax); listOfTH2.push_back(hHm23Comp);      
  TH2F *hHm13Comp = new TH2F("H_m_comp_13", "m_{hbb}", mBins, hmmin, hmmax, mBins, hmmin, hmmax); listOfTH2.push_back(hHm13Comp);
  TH2F*hbbdeltaPhi12Comp = new TH2F("bb_#Delta#phi_comp_12", "#Delta#phi_{bb}", phiBins, -TMath::Pi(),+TMath::Pi(), phiBins, -TMath::Pi(),+TMath::Pi()); listOfTH2.push_back(hbbdeltaPhi12Comp);
  TH2F*hbbdeltaPhi23Comp = new TH2F("bb_#Delta#phi_comp_23", "#Delta#phi_{bb}", phiBins, -TMath::Pi(),+TMath::Pi(), phiBins, -TMath::Pi(),+TMath::Pi()); listOfTH2.push_back(hbbdeltaPhi23Comp);
  TH2F*hbbdeltaPhi13Comp = new TH2F("bb_#Delta#phi_comp_13", "#Delta#phi_{bb}", phiBins, -TMath::Pi(),+TMath::Pi(), phiBins, -TMath::Pi(),+TMath::Pi()); listOfTH2.push_back(hbbdeltaPhi13Comp); 
  TH2F*hbbdeltaEta12Comp = new TH2F("bb_#Delta#eta_comp_12", "#Delta#eta_{bb}", etaBins, hetamin, hetamax, etaBins, hetamin, hetamax); listOfTH2.push_back(hbbdeltaEta12Comp);
  TH2F*hbbdeltaEta23Comp = new TH2F("bb_#Delta#eta_comp_23", "#Delta#eta_{bb}", etaBins, hetamin, hetamax, etaBins, hetamin, hetamax); listOfTH2.push_back(hbbdeltaEta23Comp);
  TH2F*hbbdeltaEta13Comp = new TH2F("bb_#Delta#eta_comp_13", "#Delta#eta_{bb}", etaBins, hetamin, hetamax, etaBins, hetamin, hetamax); listOfTH2.push_back(hbbdeltaEta13Comp);

  // vbfj
  TH2F*hjjpT12Comp = new TH2F("jj_pT_comp_12", "p_{T}^{jj}", pTBins, jpTmin, jpTmax, pTBins, jpTmin, jpTmax); listOfTH2.push_back(hjjpT12Comp);
  TH2F*hjjpT23Comp = new TH2F("jj_pT_comp_23", "p_{T}^{jj}", pTBins, jpTmin, jpTmax, pTBins, jpTmin, jpTmax); listOfTH2.push_back(hjjpT23Comp);
  TH2F*hjjpT13Comp = new TH2F("jj_pT_comp_13", "p_{T}^{jj}", pTBins, jpTmin, jpTmax, pTBins, jpTmin, jpTmax); listOfTH2.push_back(hjjpT13Comp);
  TH2F*hjjdeltaPhi12Comp = new TH2F("jj_#Delta#phi_comp_12", "#Delta#phi_{jj}", phiBins, -TMath::Pi(),+TMath::Pi(), phiBins, -TMath::Pi(),+TMath::Pi()); listOfTH2.push_back(hjjdeltaPhi12Comp);
  TH2F*hjjdeltaPhi23Comp = new TH2F("jj_#Delta#phi_comp_23", "#Delta#phi_{jj}", phiBins, -TMath::Pi(),+TMath::Pi(), phiBins, -TMath::Pi(),+TMath::Pi()); listOfTH2.push_back(hjjdeltaPhi23Comp);
  TH2F*hjjdeltaPhi13Comp = new TH2F("jj_#Delta#phi_comp_13", "#Delta#phi_{jj}", phiBins, -TMath::Pi(),+TMath::Pi(), phiBins, -TMath::Pi(),+TMath::Pi()); listOfTH2.push_back(hjjdeltaPhi13Comp);

  // leps
  TH2F*hl1pT12Comp = new TH2F("l1_pT_comp_12", "p^{T}_{l1}", pTBins, lpTmin, lpTmax, pTBins, lpTmin, lpTmax); listOfTH2.push_back(hl1pT12Comp);
  TH2F*hl1pT23Comp = new TH2F("l1_pT_comp_23", "p^{T}_{l1}", pTBins, lpTmin, lpTmax, pTBins, lpTmin, lpTmax); listOfTH2.push_back(hl1pT23Comp);
  TH2F*hl1pT13Comp = new TH2F("l1_pT_comp_13", "p^{T}_{l1}", pTBins, lpTmin, lpTmax, pTBins, lpTmin, lpTmax); listOfTH2.push_back(hl1pT13Comp);
  TH2F*hl2pT12Comp = new TH2F("l2_pT_comp_12", "p^{T}_{l2}", pTBins, lpTmin, lpTmax, pTBins, lpTmin, lpTmax); listOfTH2.push_back(hl2pT12Comp);
  TH2F*hl2pT23Comp = new TH2F("l2_pT_comp_23", "p^{T}_{l2}", pTBins, lpTmin, lpTmax, pTBins, lpTmin, lpTmax); listOfTH2.push_back(hl2pT23Comp);
  TH2F*hl2pT13Comp = new TH2F("l2_pT_comp_13", "p^{T}_{l2}", pTBins, lpTmin, lpTmax, pTBins, lpTmin, lpTmax); listOfTH2.push_back(hl2pT13Comp);

  // z comps
  TH2F*hz1pT12Comp = new TH2F("z1_pT_comp_12", "p^{T}_{z1}", pTBins, zpTmin, zpTmax, pTBins, zpTmin, zpTmax); listOfTH2.push_back(hz1pT12Comp);
  TH2F*hz1pT23Comp = new TH2F("z1_pT_comp_23", "p^{T}_{z1}", pTBins, zpTmin, zpTmax, pTBins, zpTmin, zpTmax); listOfTH2.push_back(hz1pT23Comp);
  TH2F*hz1pT13Comp = new TH2F("z1_pT_comp_13", "p^{T}_{z1}", pTBins, zpTmin, zpTmax, pTBins, zpTmin, zpTmax); listOfTH2.push_back(hz1pT13Comp);
  TH2F*hz2pT12Comp = new TH2F("z2_pT_comp_12", "p^{T}_{z2}", pTBins, zpTmin, zpTmax, pTBins, zpTmin, zpTmax); listOfTH2.push_back(hz2pT12Comp);
  TH2F*hz2pT23Comp = new TH2F("z2_pT_comp_23", "p^{T}_{z2}", pTBins, zpTmin, zpTmax, pTBins, zpTmin, zpTmax); listOfTH2.push_back(hz2pT23Comp);
  TH2F*hz2pT13Comp = new TH2F("z2_pT_comp_13", "p^{T}_{z2}", pTBins, zpTmin, zpTmax, pTBins, zpTmin, zpTmax); listOfTH2.push_back(hz2pT13Comp);
  TH2F*hz1m12Comp = new TH2F("z1_m_comp_12", "m_{z1}", mBins, zmmin, zmmax, mBins, zmmin, zmmax); listOfTH2.push_back(hz1m12Comp);
  TH2F*hz1m23Comp = new TH2F("z1_m_comp_23", "m_{z1}", mBins, zmmin, zmmax, mBins, zmmin, zmmax); listOfTH2.push_back(hz1m23Comp);      
  TH2F*hz1m13Comp = new TH2F("z1_m_comp_13", "m_{z1}", mBins, zmmin, zmmax, mBins, zmmin, zmmax); listOfTH2.push_back(hz1m13Comp);
  TH2F*hz2m12Comp = new TH2F("z2_m_comp_12", "m_{z2}", mBins, zmmin, zmmax, mBins, zmmin, zmmax); listOfTH2.push_back(hz2m12Comp);
  TH2F*hz2m23Comp = new TH2F("z2_m_comp_23", "m_{z2}", mBins, zmmin, zmmax, mBins, zmmin, zmmax); listOfTH2.push_back(hz2m23Comp);
  TH2F*hz2m13Comp = new TH2F("z2_m_comp_13", "m_{z2}", mBins, zmmin, zmmax, mBins, zmmin, zmmax); listOfTH2.push_back(hz2m13Comp);

  // w comps
  TH2F*hw1pT12Comp = new TH2F("w1_pT_comp_12", "p^{T}_{w1}", pTBins, wpTmin, wpTmax, pTBins, wpTmin, wpTmax); listOfTH2.push_back(hw1pT12Comp);
  TH2F*hw1pT23Comp = new TH2F("w1_pT_comp_23", "p^{T}_{w1}", pTBins, wpTmin, wpTmax, pTBins, wpTmin, wpTmax); listOfTH2.push_back(hw1pT23Comp);
  TH2F*hw1pT13Comp = new TH2F("w1_pT_comp_13", "p^{T}_{w1}", pTBins, wpTmin, wpTmax, pTBins, wpTmin, wpTmax); listOfTH2.push_back(hw1pT13Comp);
  TH2F*hw2pT12Comp = new TH2F("w2_pT_comp_12", "p^{T}_{w2}", pTBins, wpTmin, wpTmax, pTBins, wpTmin, wpTmax); listOfTH2.push_back(hw2pT12Comp);
  TH2F*hw2pT23Comp = new TH2F("w2_pT_comp_23", "p^{T}_{w2}", pTBins, wpTmin, wpTmax, pTBins, wpTmin, wpTmax); listOfTH2.push_back(hw2pT23Comp);
  TH2F*hw2pT13Comp = new TH2F("w2_pT_comp_13", "p^{T}_{w2}", pTBins, wpTmin, wpTmax, pTBins, wpTmin, wpTmax); listOfTH2.push_back(hw2pT13Comp);
  TH2F*hw1m12Comp = new TH2F("w1_m_comp_12", "m_{w1}", mBins, wmmin, wmmax, mBins, wmmin, wmmax); listOfTH2.push_back(hw1m12Comp);
  TH2F*hw1m23Comp = new TH2F("w1_m_comp_23", "m_{w1}", mBins, wmmin, wmmax, mBins, wmmin, wmmax); listOfTH2.push_back(hw1m23Comp);      
  TH2F*hw1m13Comp = new TH2F("w1_m_comp_13", "m_{w1}", mBins, wmmin, wmmax, mBins, wmmin, wmmax); listOfTH2.push_back(hw1m13Comp);
  TH2F*hw2m12Comp = new TH2F("w2_m_comp_12", "m_{w2}", mBins, wmmin, wmmax, mBins, wmmin, wmmax); listOfTH2.push_back(hw2m12Comp);
  TH2F*hw2m23Comp = new TH2F("w2_m_comp_23", "m_{w2}", mBins, wmmin, wmmax, mBins, wmmin, wmmax); listOfTH2.push_back(hw2m23Comp);
  TH2F*hw2m13Comp = new TH2F("w2_m_comp_13", "m_{w2}", mBins, wmmin, wmmax, mBins, wmmin, wmmax); listOfTH2.push_back(hw2m13Comp);


  // mixed
      
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

  //TProfile *kappaLambda = new TProfile("kappaLambda", "kappaLambda", 40, -20, 20);
  // kappaLambda -> GetXaxis() -> SetTitle("#kappa_{#lambda}");

  TH1F *recoET = new TH1F("reco_event_type", "ET", 5, -1, 3); listOfTH1.push_back(recoET);
  TH1F *particleET = new TH1F("particle_event_type", "ET", 5, -1, 3); listOfTH1.push_back(particleET);
  TH1F *partonET = new TH1F("parton_event_type", "ET", 5, -1, 3); listOfTH1.push_back(partonET);




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

  // leps
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

  TLorentzVector fourl_reco, fourl_particle, fourl_parton;

  // z
  TLorentzVector z1_reco, z1_particle,  z1_parton;
  TLorentzVector z2_reco, z2_particle,  z2_parton;

  TLorentzVector w1_reco, w1_particle,  w1_parton;
  TLorentzVector w2_reco, w2_particle,  w2_parton;
  TLorentzVector v1_reco, v1_particle,  v1_parton;
  TLorentzVector v2_reco, v2_particle,  v2_parton;
  TLorentzVector MET, met1, met2;

// kinematic quantities

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
  
  for(Int_t entry = 0; entry < numberOfEntries; ++entry) {

    //if( entry > 10) break;

    treeReader->ReadEntry(entry);
    // std::map<int,double> kappaLambdaWeights;
    HepMCEvent *event = (HepMCEvent*) branchEvent -> At(0);
    Float_t weight = event->Weight*Lumi*cross_section*numberOfEntries/(numEntries*totalWeight);
    Float_t test_weight = event->Weight*cross_section*numberOfEntries/(numEntries*totalWeight);
    hWeight -> Fill(event->Weight, test_weight);

  //------------------------------------------------------------------------------------------------------------------------------------------------------------
  // RECO - HIGGS 
  //------------------------------------------------------------------------------------------------------------------------------------------------------------


    int switchVal_reco = 0;
    bool foundBjet = false;

    if(enableCutReco["initial - reco"]){
      increaseCount(cutFlowMap_reco,"initial - reco",weight);
    }
    
    // vector <int> btagIndex;
    // vector <int> noBtag;
    // vector <int> goodJetIndex=GoodJetIndices(btagIndex,noBtag,branchJet,branchGenParticle);

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
      else  switchVal_reco = 1;
    }

    for(int i=0; i<(int)pairedJet.size(); i++){
    std::pair< std::map<TString, float>, std::map<TString, std::vector<float>>> thisPaired=pairedJet.at(i);
      if( thisPaired.first["isbtagged"] > 0) {
        pairedJetB.push_back(thisPaired);
      }
    }

    vector <int> btagIndex;
    int pairedJetSize_reco = pairedJet.size();
    int pairedBJetSize_reco = pairedJetB.size();

    if(enableCutReco["1 bb PAIReD jet - reco"]) {
        if(switchVal_reco == 0 && pairedJetB.size()>0){
          increaseCount(cutFlowMap_reco,"1 bb PAIReD jet - reco",weight);
          foundBjet = true;
        } else switchVal_reco = 1;
    }

    std::map<TString, float> paired_jet;

    if (switchVal_reco == 0 && pairedJetB.size()>0){

      foundBjet = true;
      paired_jet = pairedJetB.at(0).first;

      btagIndex.push_back(paired_jet["jet1_index"]);
      btagIndex.push_back(paired_jet["jet2_index"]);

      b1_reco.SetPtEtaPhiM(paired_jet["jet1_pt"],paired_jet["jet1_eta"],paired_jet["jet1_phi"],paired_jet["jet1_mass"]);
      b2_reco.SetPtEtaPhiM(paired_jet["jet2_pt"],paired_jet["jet2_eta"],paired_jet["jet2_phi"],paired_jet["jet2_mass"]);

      h_reco = b1_reco + b2_reco; // dijet

      if (b1_reco.Eta() > b2_reco.Eta()) {

        bbdeltaPhireco = remainder( b1_reco.Phi() - b2_reco.Phi(), 2*M_PI );
        bbdeltaEtareco = b1_reco.Eta() - b2_reco.Eta();
        } else {
        bbdeltaPhireco = remainder( b2_reco.Phi() - b1_reco.Phi(), 2*M_PI );
        bbdeltaEtareco = b2_reco.Eta() - b1_reco.Eta();

      }

      bbdeltaRreco=sqrt((bbdeltaPhireco*bbdeltaPhireco)+(bbdeltaEtareco*bbdeltaEtareco));

      //Double_t bbdeltaPhireco = b1_reco.DeltaPhi(b2_reco);
      //Double_t bbdeltaEtareco = b1_reco.DeltaEta(b2_reco);
      //Double_t bbdeltaRreco = b1_reco.DeltaR(b2_reco, kFALSE);

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
      if( fabs((((Jet*)branchJet->At(vbfJetIndex[i].first))->Eta - ((Jet*)branchJet->At(vbfJetIndex[i].second))->Eta)) <= 2.5 ) {
	      continue; 
      } else { 
	      vbfJetIndex_dEta.push_back(vbfJetIndex.at(i));
	      vbfJetIndexCandidate = i;
      }
    }

    if(enableCutReco["2.5 deltaEta VBF jet - reco"]) {
      if(switchVal_reco==0 && vbfJetIndex_dEta.size()>0) increaseCount(cutFlowMap_reco,"2.5 deltaEta VBF jet - reco",weight);
      else switchVal_reco=1;
    }

    // sort them again by eta for leading/subleading
    SortByEtaIndices(vbfJetIndex_dEta,branchJet); 

    Jet *jet1 =nullptr;
    Jet *jet2 =nullptr;
  
    if(switchVal_reco==0 && vbfJetIndex_dEta.size()>0) {
      jet1 = (Jet*) branchJet->At(vbfJetIndex_dEta[0].first);
      jet2 = (Jet*) branchJet->At(vbfJetIndex_dEta[0].second);
      j1_reco=jet1->P4();
      j2_reco=jet2->P4();
    
      if (j1_reco.Eta() > j2_reco.Eta()) {
	      jjdeltaPhireco = remainder( j1_reco.Phi() - j2_reco.Phi(), 2*M_PI );
      } else {
	      jjdeltaPhireco = remainder( j2_reco.Phi() - j1_reco.Phi(), 2*M_PI );
      }

      // double jjdeltaPhireco=(j1_reco.Phi() > j2_reco.Phi() ? -1:+1)*TMath::Abs(j2_reco.Phi() - j1_reco.Phi());
      // double jjdeltaEtareco=(j1_reco.Eta() > j2_reco.Eta() ? -1:+1)*TMath::Abs(j2_reco.Eta() - j1_reco.Eta());
      // double jjdeltaPhireco= (j1_reco.Phi() - j2_reco.Phi());
      jjdeltaEtareco= (j1_reco.Eta() - j2_reco.Eta());
      jjdeltaRreco=sqrt((jjdeltaPhireco*jjdeltaPhireco)+(jjdeltaEtareco*jjdeltaEtareco));
    }

  //------------------------------------------------------------------------------------------------------------------------------------------------------------
  // RECO - LEPTONS
  //------------------------------------------------------------------------------------------------------------------------------------------------------------

  int thisRecoEventType=-1;

  vector<pair<int,pair<int,int>>> ZRecoPairIndices;
  vector<pair<pair<int,int>,int>> WRecoPairIndices;

  vector <int> goodE_reco_indices  = GoodElectronRecoIndices(branchElectron, analysis);
  vector <int> goodMu_reco_indices = GoodMuonRecoIndices(branchMuon, analysis);

  if(analysis == "HZZJJ"){

    if(enableCutReco["lep pT > 5 & eta < 2.5 - reco"]){
      if (switchVal_reco==0) {
       if (goodE_reco_indices.size() > 0 || goodMu_reco_indices.size() > 0) increaseCount(cutFlowMap_reco,"lep pT > 5 & eta < 2.5 - reco",weight);
      } 
      else  switchVal_reco=1;
    }

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

    recoET->Fill(thisRecoEventType,weight);

    getRecoZLeps(thisRecoEventType, ZRecoPairIndices, branchElectron, branchMuon, l1_reco, l2_reco, l3_reco, l4_reco, q1_reco, q2_reco, q3_reco, q4_reco);

    if( switchVal_reco == 0 && thisRecoEventType != -1 && ZRecoPairIndices.size() >= 2 ){
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

    // WW 

   } else if(analysis == "HWWJJ") {

      if(enableCutReco["lep pT > 15 & eta < 2.5 - reco"]){
        if (switchVal_reco==0) {
          if (goodE_reco_indices.size() > 0 || goodMu_reco_indices.size() > 0) increaseCount(cutFlowMap_reco,"lep pT > 15 & eta < 2.5 - reco",weight); 
          else switchVal_reco=1;
        }
      }
      
      WRecoPairIndices = GetWRecoPairIndices(goodE_reco_indices, goodMu_reco_indices, branchElectron, branchMuon, branchMissingET);

      // FOR OFOS SWITCH mu mu / e e EVENT TYPE TO -1
      if( switchVal_reco==0 && WRecoPairIndices.size()>=2){
          if( WRecoPairIndices[0].first.first == 1 && WRecoPairIndices[1].first.first == 1) thisRecoEventType=0; // mu mu
          else if( WRecoPairIndices[0].first.first == 0 && WRecoPairIndices[1].first.first == 0) thisRecoEventType=1; // e e
          else if( WRecoPairIndices[0].first.first == 1 && WRecoPairIndices[1].first.first == 0) thisRecoEventType=2; // mu e
          else if( WRecoPairIndices[0].first.first == 0 && WRecoPairIndices[1].first.first == 1) thisRecoEventType=3; // e mu
      }

    recoET->Fill(thisRecoEventType,weight);

    // cout << "reco event type: " << thisRecoEventType << endl; 

      if(enableCutReco["OSOF - reco"]){
        if (switchVal_reco==0 && thisRecoEventType != -1) increaseCount(cutFlowMap_reco, "OSOF - reco", weight);
        else switchVal_reco=1;
      }

      getWReco(thisRecoEventType, WRecoPairIndices, branchElectron, branchMuon, branchMissingET, l1_reco, l2_reco, q1_reco, q2_reco, met1, met2);

      //cout << "mll_reco: " << (l1_reco+l2_reco).M() << endl;

      hllm_0_15_reco->Fill((l1_reco+l2_reco).M(),weight);

      if(enableCutReco["mll > 12 - reco"]) {
          if(switchVal_reco == 0 && (l1_reco+l2_reco).M() >= 12) increaseCount(cutFlowMap_reco,"mll > 12 - reco",weight);
          else switchVal_reco = 1;
      }

    if( switchVal_reco == 0){

      w1_reco=l1_reco + met1;
      w2_reco=l2_reco + met2;
    
      l1l2deltaPhireco=(l1_reco.Phi() > l2_reco.Phi() ? -1:+1)*TMath::Abs(l2_reco.Phi() - l1_reco.Phi());
      l1l2deltaEtareco=(l1_reco.Eta() > l2_reco.Eta() ? -1:+1)*TMath::Abs(l2_reco.Eta() - l1_reco.Eta());
      l1l2deltaRreco=sqrt((l1l2deltaPhireco*l1l2deltaPhireco)+(l1l2deltaEtareco*l1l2deltaEtareco));
  
      l1cosThetareco=l1_reco.CosTheta();
      l2cosThetareco=l2_reco.CosTheta();
  
      l1_reco.Boost(-w1_reco.BoostVector());
      l2_reco.Boost(-w1_reco.BoostVector());
      l1cosThetaBoostreco=l1_reco.CosTheta();
      l2cosThetaBoostreco=l2_reco.CosTheta();
      l1l2deltaPhiBoostreco=(l1_reco.Phi() > l2_reco.Phi() ? -1:+1)*TMath::Abs(l2_reco.Phi() - l1_reco.Phi());
      l1l2deltaEtaBoostreco=(l1_reco.Eta() > l2_reco.Eta() ? -1:+1)*TMath::Abs(l2_reco.Eta() - l1_reco.Eta());
      l1_reco.Boost(w1_reco.BoostVector());
      l2_reco.Boost(w1_reco.BoostVector());

      // collins soper frame 

      l1l2CScosThetareco=(q1_reco > q2_reco ? -1:+1)*TMath::Abs(2*(l2_reco.Pz()*l1_reco.E()-l1_reco.Pz()*l2_reco.E())/(z1_reco.M()*sqrt(z1_reco.M()*z1_reco.M()+z1_reco.Pt()*z1_reco.Pt())));  
    
      wwdeltaPhireco=(w1_reco.Phi() > w2_reco.Phi() ? -1:+1)*TMath::Abs(w2_reco.Phi() - w1_reco.Phi());
      wwdeltaEtareco=(w1_reco.Eta() > w2_reco.Eta() ? -1:+1)*TMath::Abs(w2_reco.Eta() - w1_reco.Eta());
      wwdeltaRreco=sqrt((wwdeltaPhireco*wwdeltaPhireco)+(wwdeltaEtareco*wwdeltaEtareco));
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
    
    // vector <int> btagIndexParticle;
    // vector <int> noBtagParticle;
    // vector <int> goodJetIndexParticle=GoodJetIndices(btagIndexParticle,noBtagParticle,branchGenJet,branchGenParticle);

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

      if (b1_particle.Eta() > b2_particle.Eta()) {

        bbdeltaPhiparticle = remainder( b1_particle.Phi() - b2_particle.Phi(), 2*M_PI );
        bbdeltaEtaparticle = b1_particle.Eta() - b2_particle.Eta();
        } else {
        bbdeltaPhiparticle = remainder( b2_particle.Phi() - b1_particle.Phi(), 2*M_PI );
        bbdeltaEtaparticle = b2_particle.Eta() - b1_particle.Eta();

      }

      bbdeltaRparticle=sqrt((bbdeltaPhiparticle*bbdeltaPhiparticle)+(bbdeltaEtaparticle*bbdeltaEtaparticle));

      //Double_t bbdeltaPhireco = b1_reco.DeltaPhi(b2_reco);
      //Double_t bbdeltaEtareco = b1_reco.DeltaEta(b2_reco);
      //Double_t bbdeltaRreco = b1_reco.DeltaR(b2_reco, kFALSE);

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
    
      if (j1_particle.Eta() > j2_particle.Eta()) {
	      jjdeltaPhiparticle = remainder( j1_particle.Phi() - j2_particle.Phi(), 2*M_PI );
      } else {
	      jjdeltaPhiparticle = remainder( j2_particle.Phi() - j1_particle.Phi(), 2*M_PI );
      }

      // double jjdeltaPhireco=(j1_reco.Phi() > j2_reco.Phi() ? -1:+1)*TMath::Abs(j2_reco.Phi() - j1_reco.Phi());
      // double jjdeltaEtareco=(j1_reco.Eta() > j2_reco.Eta() ? -1:+1)*TMath::Abs(j2_reco.Eta() - j1_reco.Eta());
      // double jjdeltaPhireco= (j1_reco.Phi() - j2_reco.Phi());
      jjdeltaEtaparticle= (j1_particle.Eta() - j2_particle.Eta());
      jjdeltaRparticle=sqrt((jjdeltaPhiparticle*jjdeltaPhiparticle)+(jjdeltaEtaparticle*jjdeltaEtaparticle));
    }


  //------------------------------------------------------------------------------------------------------------------------------------------------------------
  // PARTICLE - LEPTONS
  //------------------------------------------------------------------------------------------------------------------------------------------------------------

  int thisParticleEventType=-1;

  vector<pair<int,pair<int,int>>> ZParticlePairIndices;
  vector<pair<pair<int,int>,int>> WParticlePairIndices;

  vector <int> goodE_particle_indices  = GoodElectronParticleIndices(branchGenParticle, analysis);
  vector <int> goodMu_particle_indices = GoodMuonParticleIndices(branchGenParticle, analysis);

  if(analysis == "HZZJJ"){

    if(enableCutParticle["lep pT > 5 & eta < 2.5 - particle"]){
      if (switchVal_particle==0) {
       if (goodE_particle_indices.size() > 0 || goodMu_particle_indices.size() > 0) increaseCount(cutFlowMap_particle,"lep pT > 5 & eta < 2.5 - particle",weight);
      } 
      else  switchVal_particle=1;
    }

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

    particleET->Fill(thisParticleEventType,weight);

    getParticleZLeps(thisParticleEventType, ZParticlePairIndices, branchGenParticle, l1_particle, l2_particle, l3_particle, l4_particle, q1_particle, q2_particle, q3_particle, q4_particle);

    if( switchVal_particle == 0 && thisParticleEventType != -1 && ZParticlePairIndices.size() >= 2 ){
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
      l1l2deltaPhiBoostparticle=(l1_particle.Phi() > l2_particle.Phi() ? -1:+1)*TMath::Abs(l2_particle.Phi() - l1_particle.Phi());
      l3l4deltaPhiBoostparticle=(l3_particle.Phi() > l4_particle.Phi() ? -1:+1)*TMath::Abs(l4_particle.Phi() - l3_particle.Phi());
      l1l2deltaEtaBoostparticle=(l1_particle.Eta() > l2_particle.Eta() ? -1:+1)*TMath::Abs(l2_particle.Eta() - l1_particle.Eta());
      l3l4deltaEtaBoostparticle=(l3_particle.Eta() > l4_particle.Eta() ? -1:+1)*TMath::Abs(l4_particle.Eta() - l3_particle.Eta());
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
    }

    // WW 

   } else if(analysis == "HWWJJ") {

      if(enableCutParticle["lep pT > 15 & eta < 2.5 - particle"]){
        if (switchVal_particle==0) {
          if (goodE_particle_indices.size() > 0 || goodMu_particle_indices.size() > 0) increaseCount(cutFlowMap_particle,"lep pT > 15 & eta < 2.5 - particle",weight);
          else  switchVal_particle=1; 
        } 
      }

      WParticlePairIndices = GetWParticlePairIndices(goodE_particle_indices, goodMu_particle_indices, branchGenParticle, branchMissingET);

      // FOR OFOS SWITCH mu mu / e e EVENT TYPE TO -1
      if( switchVal_particle==0 && WParticlePairIndices.size()>=2){
          if( WParticlePairIndices[0].first.first == 1 && WParticlePairIndices[1].first.first == 1) thisParticleEventType = 0; // mu mu
          else if( WParticlePairIndices[0].first.first == 0 && WParticlePairIndices[1].first.first == 0) thisParticleEventType = 1; // e e
          else if( WParticlePairIndices[0].first.first == 1 && WParticlePairIndices[1].first.first == 0) thisParticleEventType = 2; // mu e
          else if( WParticlePairIndices[0].first.first == 0 && WParticlePairIndices[1].first.first == 1) thisParticleEventType = 3; // e mu
      }

      particleET->Fill(thisParticleEventType,weight);

      // cout << "particle event type: " << thisParticleEventType << endl; 

      if(enableCutParticle["OSOF - particle"]){
        if (switchVal_particle==0 && thisParticleEventType != -1) increaseCount(cutFlowMap_particle,"OSOF - particle",weight);
        else switchVal_particle=1;
      }

      getWParticle(thisParticleEventType, WParticlePairIndices, branchGenParticle,branchMissingET, l1_particle, l2_particle, q1_particle, q2_particle, met1, met2);

      // cout << "mll_particle: " << (l1_particle+l2_particle).M() << endl;

      hllm_0_15_particle->Fill((l1_particle+l2_particle).M(),weight);

      if(enableCutParticle["mll > 12 - particle"]) {
        if(switchVal_particle == 0 &&  (l1_particle+l2_particle).M() >= 12) increaseCount(cutFlowMap_particle,"mll > 12 - particle",weight);
        else switchVal_particle = 1;
      }

      if( switchVal_particle == 0 ) {

        w1_particle=l1_particle + met1;
        w2_particle=l2_particle + met2;

        l1l2deltaPhiparticle=(l1_particle.Phi() > l2_particle.Phi() ? -1:+1)*TMath::Abs(l2_particle.Phi() - l1_particle.Phi());
        l1l2deltaEtaparticle=(l1_particle.Eta() > l2_particle.Eta() ? -1:+1)*TMath::Abs(l2_particle.Eta() - l1_particle.Eta());
        l1l2deltaRparticle=sqrt((l1l2deltaPhiparticle*l1l2deltaPhiparticle)+(l1l2deltaEtaparticle*l1l2deltaEtaparticle));
    
        l1cosThetaparticle=l1_particle.CosTheta();
        l2cosThetaparticle=l2_particle.CosTheta();
    
        l1_particle.Boost(-w1_particle.BoostVector());
        l2_particle.Boost(-w1_particle.BoostVector());
        l1cosThetaBoostparticle=l1_particle.CosTheta();
        l2cosThetaBoostparticle=l2_particle.CosTheta();
        l1l2deltaPhiBoostparticle=(l1_particle.Phi() > l2_particle.Phi() ? -1:+1)*TMath::Abs(l2_particle.Phi() - l1_particle.Phi());
        l1l2deltaEtaBoostparticle=(l1_particle.Eta() > l2_particle.Eta() ? -1:+1)*TMath::Abs(l2_particle.Eta() - l1_particle.Eta());
        l1_particle.Boost(w1_particle.BoostVector());
        l2_particle.Boost(w1_particle.BoostVector());

        // collins soper frame 

        l1l2CScosThetaparticle=(q1_particle > q2_particle ? -1:+1)*TMath::Abs(2*(l2_particle.Pz()*l1_particle.E()-l1_particle.Pz()*l2_particle.E())/(z1_particle.M()*sqrt(z1_particle.M()*z1_particle.M()+z1_particle.Pt()*z1_particle.Pt())));  

        wwdeltaPhiparticle=(w1_particle.Phi() > w2_particle.Phi() ? -1:+1)*TMath::Abs(w2_particle.Phi() - w1_particle.Phi());
        wwdeltaEtaparticle=(w1_particle.Eta() > w2_particle.Eta() ? -1:+1)*TMath::Abs(w2_particle.Eta() - w1_particle.Eta());
        wwdeltaRparticle=sqrt((wwdeltaPhiparticle*wwdeltaPhiparticle)+(wwdeltaEtaparticle*wwdeltaEtaparticle));

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
      if (b1_parton.Eta() > b2_parton.Eta()) {
	      bbdeltaPhiparton = remainder( b1_parton.Phi() - b2_parton.Phi(), 2*M_PI );
	      bbdeltaEtaparton= b1_parton.Eta() - b2_parton.Eta();
      } else {
	      bbdeltaPhiparton = remainder( b2_parton.Phi() - b1_parton.Phi(), 2*M_PI );
	      bbdeltaEtaparton= b2_parton.Eta() - b1_parton.Eta();
      }

      bbdeltaRparton=sqrt((bbdeltaPhiparton*bbdeltaPhiparton)+(bbdeltaEtaparton*bbdeltaEtaparton));

      if (j1_parton.Eta() > j2_parton.Eta()) {
	      jjdeltaPhiparton = remainder( j1_parton.Phi() - j2_parton.Phi(), 2*M_PI );
	      jjdeltaEtaparton= j1_parton.Eta() - j2_parton.Eta();
      } else {
	      jjdeltaPhiparton = remainder( j2_parton.Phi() - j1_parton.Phi(), 2*M_PI );
	      jjdeltaEtaparton= j2_parton.Eta() - j1_parton.Eta();
      }
    }
    
    //------------------------------------------------------------------------------------------------------------------------------------------------------------
    // PARTON - LEPTONS
    //------------------------------------------------------------------------------------------------------------------------------------------------------------

    int thisPartonEventType=-1;

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
      
      getPartonZLeps(thisPartonEventType, ZPartonIndices, branchGenParticle, z1_parton, z2_parton, l1_parton, l2_parton, l3_parton, l4_parton, q1_parton, q2_parton, q3_parton, q4_parton);

      if (foundZZ){
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

    } if(analysis == "HWWJJ") { 

      WPartonIndices = GetWPartonIndices(branchGenParticle, analysis);

      if(WPartonIndices.size() > 1) foundWW = true;

      if(enableCutParton["WW parton"]){
        if(switchVal_parton == 0 && foundWW) increaseCount(cutFlowMap_parton,"WW parton",weight);
        else switchVal_parton = 1;
      }

      getPartonWLeps(thisPartonEventType, WPartonIndices, branchGenParticle, w1_parton, w2_parton, l1_parton, l2_parton, q1_parton, q2_parton);

      partonET->Fill(thisPartonEventType,weight);

      if(foundWW) {
        l1l2deltaPhiparton=(l1_parton.Phi() > l2_parton.Phi() ? -1:+1)*TMath::Abs(l2_parton.Phi() - l1_parton.Phi());
        l1l2deltaEtaparton=(l1_parton.Eta() > l2_parton.Eta() ? -1:+1)*TMath::Abs(l2_parton.Eta() - l1_parton.Eta());
        l1l2deltaRparton=sqrt((l1l2deltaPhiparton*l1l2deltaPhiparton)+(l1l2deltaEtaparton*l1l2deltaEtaparton));
      
        l1cosThetaparton=l1_parton.CosTheta();
        l2cosThetaparton=l2_parton.CosTheta();

        l1_parton.Boost(-w1_parton.BoostVector());
        l2_parton.Boost(-w1_parton.BoostVector());
        l1cosThetaBoostparton=l1_parton.CosTheta();
        l2cosThetaBoostparton=l2_parton.CosTheta();
        l1_parton.Boost(w1_parton.BoostVector());
        l2_parton.Boost(w1_parton.BoostVector());

        // collins soper frame
        double l1l2CScosThetaparton=(q1_parton > q2_parton ? -1:+1)*TMath::Abs(2*(l2_parton.Pz()*l1_parton.E()-l1_parton.Pz()*l2_parton.E())/(w1_parton.M()*sqrt(w1_parton.M()*w1_parton.M()+w1_parton.Pt()*w1_parton.Pt())));

        double wwdeltaPhiparton=(w1_parton.Phi() > w2_parton.Phi() ? -1:+1)*TMath::Abs(w2_parton.Phi() - w1_parton.Phi());
        double wwdeltaEtaparton=(w1_parton.Eta() > w2_parton.Eta() ? -1:+1)*TMath::Abs(w2_parton.Eta() - w1_parton.Eta());
        double wwdeltaRparton=sqrt((wwdeltaPhiparton*zzdeltaPhiparton)+(wwdeltaEtaparton*zzdeltaEtaparton));
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

    // higgs - reco
    if(switchVal_reco==0){
      if(foundBjet){
        hHpTreco -> Fill(h_reco.Pt(),weight);
        hHmreco -> Fill(h_reco.M(),weight); 
	      hbbdeltaPhireco -> Fill(bbdeltaPhireco,weight);
	      hbbdeltaEtareco -> Fill(bbdeltaEtareco,weight);
	      hbbdeltaRreco -> Fill(bbdeltaRreco,weight);
      }
    }

    // vbfj - reco
    if(switchVal_reco==0){
      if(vbfJetIndex.size()>0){
        hjjpTreco->Fill(j1_reco.Pt()+j2_reco.Pt(),weight);
        hjjdeltaPhireco->Fill(jjdeltaPhireco,weight); 
        hjjdeltaEtareco->Fill(jjdeltaEtareco, weight);
        hjjdeltaRreco -> Fill(jjdeltaRreco,weight);
      }
    }

    // z - reco
    if(switchVal_reco==0){
      if(thisRecoEventType!=-1 && ZRecoPairIndices.size()>=2){
        hz1pTreco->Fill(z1_reco.Pt(),weight);
        hz2pTreco->Fill(z2_reco.Pt(),weight);
        hz1mreco->Fill(z1_reco.M(),weight);
        hz2mreco->Fill(z2_reco.M(),weight);
        hzzdeltaPhireco->Fill(zzdeltaPhireco,weight); 
        hzzdeltaEtareco->Fill(zzdeltaEtareco, weight);
        hzzdeltaRreco -> Fill(zzdeltaRreco,weight);

        hz1cosThetareco->Fill(z1_reco.CosTheta(),weight);
        hz2cosThetareco->Fill(z2_reco.CosTheta(),weight);
      }
    }

    // w - reco
    if(switchVal_reco==0){
      if(thisRecoEventType!=-1 && WRecoPairIndices.size()>=2){
        hllpTreco->Fill((l1_reco+l2_reco).Pt(),weight);
        hllmreco->Fill((l1_reco+l2_reco).M(),weight);
        hw1pTreco->Fill(w1_reco.Pt(),weight);
        hw1mreco->Fill(w1_reco.M(),weight);
        hw2pTreco->Fill(w2_reco.Pt(),weight);
        hw2mreco->Fill(w2_reco.M(),weight);
        hwwpTreco->Fill((w1_reco + w2_reco).Pt(),weight);
        hwwmreco->Fill((w1_reco + w2_reco).M(),weight);

        hwwdeltaPhireco->Fill(wwdeltaPhireco,weight); 
        hwwdeltaEtareco->Fill(wwdeltaEtareco, weight);
        hwwdeltaRreco -> Fill(wwdeltaRreco,weight);
      }
    }

    // w ET- reco 
    if(switchVal_reco==0){
      if(thisRecoEventType==0 && WRecoPairIndices.size()>=2){
        hllpTET0reco->Fill((l1_reco+l2_reco).Pt(),weight);
        hllmET0reco->Fill((l1_reco+l2_reco).M(),weight);
      } else if(thisRecoEventType==1 && WRecoPairIndices.size()>=2){
        hllpTET1reco->Fill((l1_reco+l2_reco).Pt(),weight);
        hllmET1reco->Fill((l1_reco+l2_reco).M(),weight);
      } else if(thisRecoEventType==2 && WRecoPairIndices.size()>=2){
        hllpTET2reco->Fill((l1_reco+l2_reco).Pt(),weight);
        hllmET2reco->Fill((l1_reco+l2_reco).M(),weight);
      } else if(thisRecoEventType==3 && WRecoPairIndices.size()>=2){
        hllpTET3reco->Fill((l1_reco+l2_reco).Pt(),weight);
        hllmET3reco->Fill((l1_reco+l2_reco).M(),weight);
      }
    }

    // leptons - reco
    if(switchVal_reco==0){
      if(thisRecoEventType!=-1 && ZRecoPairIndices.size()>=2){
        hl1l2deltaPhireco->Fill(l1l2deltaPhireco ,weight);
        hl3l4deltaPhireco->Fill(l3l4deltaPhireco ,weight);
        hl1l2deltaPhiBoostreco->Fill(l1l2deltaPhiBoostreco ,weight);
        hl3l4deltaPhiBoostreco->Fill(l3l4deltaPhiBoostreco ,weight);
        hl1l2deltaEtareco->Fill(l1l2deltaEtareco ,weight);
        hl3l4deltaEtareco->Fill(l3l4deltaEtareco ,weight);
        hl1l2deltaEtaBoostreco->Fill(l1l2deltaEtaBoostreco ,weight);
        hl3l4deltaEtaBoostreco->Fill(l3l4deltaEtaBoostreco ,weight);
        hl1l2deltaRreco->Fill(l1l2deltaRreco ,weight);
        hl3l4deltaRreco->Fill(l3l4deltaRreco ,weight);

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
      }
  

  //------------------------------------------------------------------------------------------------------------------------------------------------------------
  // FILL HISTOGRAMS - PARTICLE
  //------------------------------------------------------------------------------------------------------------------------------------------------------------


      // higgs - particle
    if(switchVal_particle==0){
      if(foundBjetParticle){
        hHpTparticle -> Fill(h_particle.Pt(), weight);
        hHmparticle -> Fill(h_particle.M(), weight);
        hbbdeltaPhiparticle->Fill(bbdeltaPhiparticle,weight);
        hbbdeltaEtaparticle -> Fill(bbdeltaEtaparticle,weight);
        hbbdeltaRparticle -> Fill(bbdeltaRparticle,weight);
      }
    }

    // vbfj - particle
    if(switchVal_particle==0){
      if(vbfJetIndexParticle.size() > 0){
        hjjpTparticle->Fill(j1_particle.Pt()+j2_particle.Pt(),weight);
        hjjdeltaPhiparticle->Fill(jjdeltaPhiparticle,weight); 
        hjjdeltaEtaparticle->Fill(jjdeltaEtaparticle, weight);
        hjjdeltaRparticle -> Fill(jjdeltaRparticle,weight);
      }
    }

    // z - particle
    if(switchVal_particle==0){
      if(thisParticleEventType!=-1 && ZParticlePairIndices.size()>=2){
        hz1pTparticle->Fill(z1_particle.Pt(),weight);
        hz2pTparticle->Fill(z2_particle.Pt(),weight);
        hz1mparticle->Fill(z1_particle.M(),weight);
        hz2mparticle->Fill(z2_particle.M(),weight);
        hzzdeltaPhiparticle->Fill(zzdeltaPhiparticle,weight); 
        hzzdeltaEtaparticle->Fill(zzdeltaEtaparticle, weight);
        hzzdeltaRparticle -> Fill(zzdeltaRparticle,weight);

        hz1cosThetaparticle->Fill(z1_particle.CosTheta(),weight);
        hz2cosThetaparticle->Fill(z2_particle.CosTheta(),weight);
      }
    }

    // w - particle
    if(switchVal_particle==0){
      if(thisParticleEventType!=-1 && WParticlePairIndices.size()>=2){
        hllpTparticle->Fill((l1_particle+l2_particle).Pt(),weight);
        hllmparticle->Fill((l1_particle+l2_particle).M(),weight);
        hw1pTparticle->Fill(w1_particle.Pt(),weight);
        hw1mparticle->Fill(w1_particle.M(),weight);
        hw2pTparticle->Fill(w2_particle.Pt(),weight);
        hw2mparticle->Fill(w2_particle.M(),weight);
        hwwpTparticle->Fill((w1_particle + w2_particle).Pt(),weight);
        hwwmparticle->Fill((w1_particle + w2_particle).M(),weight);

        hwwdeltaPhiparticle->Fill(wwdeltaPhiparticle,weight); 
        hwwdeltaEtaparticle->Fill(wwdeltaEtaparticle, weight);
        hwwdeltaRparticle -> Fill(wwdeltaRparticle,weight);
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

    // leptons - particle
      if(thisParticleEventType!=-1){
        hl1l2deltaPhiparticle->Fill(l1l2deltaPhiparticle ,weight);
        hl3l4deltaPhiparticle->Fill(l3l4deltaPhiparticle ,weight);
        hl1l2deltaPhiBoostparticle->Fill(l1l2deltaPhiBoostparticle ,weight);
        hl3l4deltaPhiBoostparticle->Fill(l3l4deltaPhiBoostparticle ,weight);
        hl1l2deltaEtaparticle->Fill(l1l2deltaEtaparticle ,weight);
        hl3l4deltaEtaparticle->Fill(l3l4deltaEtaparticle ,weight);
        hl1l2deltaEtaBoostparticle->Fill(l1l2deltaEtaBoostparticle ,weight);
        hl3l4deltaEtaBoostparticle->Fill(l3l4deltaEtaBoostparticle ,weight);
        hl1l2deltaRparticle->Fill(l1l2deltaRparticle ,weight);
        hl3l4deltaRparticle->Fill(l3l4deltaRparticle ,weight);

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
        hl3l4CScosThetaparticle->Fill(l3l4CScosThetareco,weight);
      }
    

  //------------------------------------------------------------------------------------------------------------------------------------------------------------
  // FILL HISTOGRAMS - PARTON
  //------------------------------------------------------------------------------------------------------------------------------------------------------------

    // higgs - parton
    if(switchVal_parton==0){
      if(HiggsRecord){
        hHpTparton -> Fill(h_parton.Pt(), weight);
        hHmparton -> Fill(h_parton.M(), weight);
        hbbdeltaPhiparton->Fill(bbdeltaPhiparton,weight);
        hbbdeltaEtaparton -> Fill(bbdeltaEtaparton,weight);
        hbbdeltaRparton -> Fill(bbdeltaRparton,weight);
      }
    }

    // z - parton
    if(switchVal_parton==0 ){
      if(foundZZ){
        hz1pTparton->Fill(z1_parton.Pt(),weight);
        hz2pTparton->Fill(z2_parton.Pt(),weight);
        hz1mparton->Fill(z1_parton.M(),weight);
        hz2mparton->Fill(z2_parton.M(),weight);
        hzzdeltaPhiparton->Fill(zzdeltaPhiparton,weight); 
        hzzdeltaEtaparton->Fill(zzdeltaEtaparton, weight);
        hzzdeltaRparton -> Fill(zzdeltaRparton,weight);
      }
    }

    // w - parton
    if(switchVal_parton==0 ){
      if(foundWW){
        hw1pTparton->Fill(w1_parton.Pt(),weight);
        hw2pTparton->Fill(w2_parton.Pt(),weight);
        hw1mparton->Fill(w1_parton.M(),weight);
        hw2mparton->Fill(w2_parton.M(),weight);
        hwwdeltaPhiparton->Fill(wwdeltaPhiparton,weight); 
        hwwdeltaEtaparton->Fill(wwdeltaEtaparton, weight);
        hwwdeltaRparton -> Fill(wwdeltaRparton,weight);
      }
    }

    // leptons - parton
    if(switchVal_parton==0 ){
      if(foundZZ){
        hl1l2deltaPhiparton->Fill(l1l2deltaPhiparton ,weight);
        hl3l4deltaPhiparton->Fill(l3l4deltaPhiparton ,weight);
        hl1l2deltaEtaparton->Fill(l1l2deltaEtaparton ,weight);
        hl3l4deltaEtaparton->Fill(l3l4deltaEtaparton ,weight);
        hl1l2deltaRparton->Fill(l1l2deltaRparton ,weight);
        hl3l4deltaRparton->Fill(l3l4deltaRparton ,weight);

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
