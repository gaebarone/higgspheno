#include "classes/DelphesClasses.h"
#include "ExRootAnalysis/ExRootTreeReader.h"
#include "ghost_tagging.h"
#include <iostream>
#include <fstream>
#include <string>
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
#include <vector>
#include "TClonesArray.h"
#include "get_cross_section.h"
//------------------------------------------------------------------------------

Long64_t get_total_events(const char *process_name) {
  std::string inputFileName = std::string(process_name) + "_inputs.txt";
  std::ifstream inputFile(inputFileName.c_str());
  std::string line;
  TChain chain("Delphes");
  Long64_t total = 0;
  while (std::getline(inputFile, line)) {
    chain.Add(line.c_str());
  }
  ExRootTreeReader *treeReader = new ExRootTreeReader(&chain);
  Long64_t numEntries = treeReader->GetEntries();
  delete treeReader;
  return numEntries;
}

// be able to calculate total weight of events in a process
/*Float_t get_file_weight(const char *inputName) {
  TChain chain("Delphes");
  chain.Add(inputName);
  ExRootTreeReader *treeReader = new ExRootTreeReader(&chain);
  Float_t file_weight = 0;
  Long64_t numberOfEntries = treeReader->GetEntries();
  TClonesArray *branchEvent = (TClonesArray*)treeReader->UseBranch("Event");
  for(Int_t entry = 0; entry < numberOfEntries; ++entry){
    // Load selected branches with data from specified event
    treeReader->ReadEntry(entry);
    HepMCEvent *event = (HepMCEvent*) branchEvent -> At(0);
    file_weight += event->Weight;
  }
  return file_weight;
}

Float_t get_total_weight(const char *process_name) {
  std::string inputFileName = std::string(process_name) + "_inputs.txt";
  std::ifstream inputFile(inputFileName.c_str());
  std::string line;
  TChain chain("Delphes");
  Float_t total = 0;
  while (std::getline(inputFile, line)) {
    total += get_file_weight(line.c_str());
  }
  return total;
}*/
// make a ton of plots for zhbb events (z -> l l) (h -> b b)
void zhbb_analyze(const char *inputFile, const char *outputFile, const char *process_name) {
  // SET CUTS
  const double e_pt_cut_lead = 27;
  const double mu_pt_cut_lead = 20;
  const double e_pt_cut_sub = 15;
  const double mu_pt_cut_sub = 12;
  const double eta_cut = 2.5;
  const double jet_pt_cut = 20;
  const double higgs_mass_cut_low = 0;
  const double higgs_mass_cut_hi = 250;
  const double higgs_pt_cut = 1e-6;
  const double z_mass_cut_low = 80;
  const double z_mass_cut_hi = 100;
  const double z_pt_cut = 20;
  const double dphi_zh_cut = 2.5;
  // get process cross section
  double Lumi=200e3;
  // Create chain of and append the file
  TChain chain("Delphes");
  chain.Add(inputFile);
  // Create object of class ExRootTreeReader
  ExRootTreeReader *treeReader = new ExRootTreeReader(&chain);
  Long64_t numberOfEntries = treeReader->GetEntries();
  Long64_t numEntries = get_total_events(process_name);
  double cross_section = get_cross_section(process_name);
  Float_t totalWeight = 0.0;
  // Get pointers to branches used in this analysis
  TClonesArray *branchJet = (TClonesArray*)treeReader->UseBranch("Jet");
  TClonesArray *branchElectron = (TClonesArray*)treeReader->UseBranch("Electron");
  TClonesArray *branchMuon = (TClonesArray*)treeReader->UseBranch("Muon");
  TClonesArray *branchEvent = (TClonesArray*)treeReader->UseBranch("Event");
  TClonesArray *branchGenParticle = (TClonesArray*)treeReader->UseBranch("Particle");
  TClonesArray *branchGenJet = (TClonesArray*) treeReader->UseBranch("GenJet");
  // Book histograms
  TH1 *hWeight = new TH1F("weights", "weight", 50, 0.0, 1.0);
  // mass
  TH1 *hMZH = new TH1F("mass_ZH", "m_{ZH}", 50, 120.0, 670.0);
  TH1 *hMbbR = new TH1F("mass_bb_reco", "Reco m_{bb}", 50, higgs_mass_cut_low, higgs_mass_cut_hi);
  TH1 *hMllR = new TH1F("mass_ll_reco", "Reco m_{ll}", 50, z_mass_cut_low, z_mass_cut_hi);
  TH1 *hMbbllR = new TH1F("mass_bbll_reco", "Reco m_{bbll}", 50, 120.0, 670.0);
  TH1 *hMbbP = new TH1F("mass_bb_particle", "m_{bb}", 50, higgs_mass_cut_low, higgs_mass_cut_hi);
  TH1 *hMllP = new TH1F("mass_ll_particle", "m_{ll}", 50, z_mass_cut_low, z_mass_cut_hi);
  TH1 *hMbbllP = new TH1F("mass_bbll_particle", "m_{bbll}", 50, 120.0, 670.0);
  TH2 *hMbbComp = new TH2F("mass_bb_Comp", "m_{bb}", 50, higgs_mass_cut_low, higgs_mass_cut_hi, 50, higgs_mass_cut_low, higgs_mass_cut_hi);
  TH2 *hMllComp = new TH2F("mass_ll_Comp", "m_{ll}", 50, z_mass_cut_low, z_mass_cut_hi, 50, z_mass_cut_low, z_mass_cut_hi);
  TH2 *hMZHComp = new TH2F("mass_ZH_Comp", "m_{ZH}", 50, 120.0, 670.0, 50, 120.0, 670.0);
  TH2 *hMZHCompP = new TH2F("mass_ZH_Comp_particle", "m_{ZH}", 50, 120.0, 670.0, 50, 120.0, 670.0);
  // pt parton level
  TH1 *hPtH = new TH1F("pt_H", "p_{T}^{H}", 50, 0.0, 400.0);
  TH1 *hPtZ = new TH1F("pt_Z", "p_{T}^{Z}", 50, 0.0, 400.0);
  TH1 *hPtLE = new TH1F("pt_LE",  "Lead p_{T}^{e}", 50, 20.0, 200.0);
  TH1 *hPtLM = new TH1F("pt_LM",  "Lead p_{T}^{#mu}", 50, 20.0, 200.0);
  TH1 *hPtSLE = new TH1F("pt_SLE",  "Sublead p_{T}^{e}", 50, 20.0, 200.0);
  TH1 *hPtSLM = new TH1F("pt_SLM",  "Sublead p_{T}^{#mu}", 50, 20.0, 200.0);
  TH1 *hPtLB = new TH1F("pt_LB", "Lead b Parton p_{T}", 50, 20.0, 200.0);
  TH1 *hPtSLB = new TH1F("pt_SLB", "Sublead b Parton p_{T}", 50, 20.0, 200.0);
  TH1 *hPtZH = new TH1F("pt_ZH", "p_{T}^{ZH}", 50, 0.0, 400.0);
  // pt reco level
  TH1 *hPtBBR = new TH1F("pt_bb_reco", "Reco p_{T}^{bb}", 50, 0.0, 400.0);
  TH1 *hPtLLR = new TH1F("pt_ll_reco", "Reco p_{T}^{ll}", 50, 0.0, 400.0);
  TH1 *hPtLER = new TH1F("pt_LE_reco",  "Lead Reco p_{T}^{e}", 50, 20.0, 200.0);
  TH1 *hPtLMR = new TH1F("pt_LM_reco",  "Lead Reco p_{T}^{#mu}", 50, 20.0, 200.0);
  TH1 *hPtSLER = new TH1F("pt_SLE_reco",  "Sublead Reco p_{T}^{e}", 50, 20.0, 200.0);
  TH1 *hPtSLMR = new TH1F("pt_SLM_reco",  "Sublead Reco p_{T}^{#mu}", 50, 20.0, 200.0);
  TH1 *hPtLBR = new TH1F("pt_LB_reco", "Lead Reco b Jet p_{T}", 50, 20.0, 200.0);
  TH1 *hPtSLBR = new TH1F("pt_SLB_reco", "Sublead Reco b Jet p_{T}", 50, 20.0, 200.0);
  TH1 *hPtBBLLR = new TH1F("pt_bbll_reco", "Reco p_{T}^{bbll}", 50, 0.0, 400.0);
  // pt particle level
  TH1 *hPtBBP = new TH1F("pt_bb_particle", "p_{T}^{bb}", 50, 0.0, 400.0);
  TH1 *hPtLLP = new TH1F("pt_ll_particle", "p_{T}^{ll}", 50, 0.0, 400.0);
  TH1 *hPtLEP = new TH1F("pt_LE_particle",  "Lead p_{T}^{e}", 50, 20.0, 200.0);
  TH1 *hPtLMP = new TH1F("pt_LM_particle",  "Lead p_{T}^{#mu}", 50, 20.0, 200.0);
  TH1 *hPtSLEP = new TH1F("pt_SLE_particle",  "Sublead p_{T}^{e}", 50, 20.0, 200.0);
  TH1 *hPtSLMP = new TH1F("pt_SLM_particle",  "Sublead p_{T}^{#mu}", 50, 20.0, 200.0);
  TH1 *hPtLBP = new TH1F("pt_LB_particle", "Lead b Jet p_{T}", 50, 20.0, 200.0);
  TH1 *hPtSLBP = new TH1F("pt_SLB_particle", "Sublead b Jet p_{T}", 50, 20.0, 200.0);
  TH1 *hPtBBLLP = new TH1F("pt_bbll_particle", "p_{T}^{bbll}", 50, 0.0, 400.0);
  // pt parton-reco comparison
  TH2 *hPtHComp = new TH2F("pt_H_Comp", "p_{T}^{H}", 50, 0.0, 400.0, 50, 0.0, 400.0);
  TH2 *hPtZComp = new TH2F("pt_Z_Comp", "p_{T}^{Z}", 50, 0.0, 400.0, 50, 0.0, 400.0);
  TH2 *hPtLEComp = new TH2F("pt_LE_Comp", "Lead p_{T}^{e}", 50, 20.0, 200.0, 50, 20.0, 200.0);
  TH2 *hPtLMComp = new TH2F("pt_LM_Comp", "Lead p_{T}^{#mu}", 50, 20.0, 200.0, 50, 20.0, 200.0);
  TH2 *hPtSLEComp = new TH2F("pt_SLE_Comp", "Sublead p_{T}^{e}", 50, 20.0, 200.0, 50, 20.0, 200.0);
  TH2 *hPtSLMComp = new TH2F("pt_SLM_Comp", "Sublead p_{T}^{#mu}", 50, 20.0, 200.0, 50, 20.0, 200.0);
  TH2 *hPtLBComp = new TH2F("pt_LB_Comp", "Lead p_{T}^{b}", 50, 20.0, 200.0, 50, 20.0, 200.0);
  TH2 *hPtSLBComp = new TH2F("pt_SLB_Comp", "Sublead p_{T}^{b}", 50, 20.0, 200.0, 50, 20.0, 200.0);
  TH2 *hPtZHComp = new TH2F("pt_ZH_Comp", "Sublead p_{T}^{ZH}", 50, 0.0, 400.0, 50, 0.0, 400.0);
  // pt particle-reco comparison
  TH2 *hPtHCompP = new TH2F("pt_H_Comp_particle", "p_{T}^{H}", 50, 0.0, 400.0, 50, 0.0, 400.0);
  TH2 *hPtZCompP = new TH2F("pt_Z_Comp_particle", "p_{T}^{Z}", 50, 0.0, 400.0, 50, 0.0, 400.0);
  TH2 *hPtLECompP = new TH2F("pt_LE_Comp_particle", "Lead p_{T}^{e}", 50, 20.0, 200.0, 50, 20.0, 200.0);
  TH2 *hPtLMCompP = new TH2F("pt_LM_Comp_particle", "Lead p_{T}^{#mu}", 50, 20.0, 200.0, 50, 20.0, 200.0);
  TH2 *hPtSLECompP = new TH2F("pt_SLE_Comp_particle", "Sublead p_{T}^{e}", 50, 20.0, 200.0, 50, 20.0, 200.0);
  TH2 *hPtSLMCompP = new TH2F("pt_SLM_Comp_particle", "Sublead p_{T}^{#mu}", 50, 20.0, 200.0, 50, 20.0, 200.0);
  TH2 *hPtLBCompP = new TH2F("pt_LB_Comp_particle", "Lead p_{T}^{b}", 50, 20.0, 200.0, 50, 20.0, 200.0);
  TH2 *hPtSLBCompP = new TH2F("pt_SLB_Comp_particle", "Sublead p_{T}^{b}", 50, 20.0, 200.0, 50, 20.0, 200.0);
  TH2 *hPtZHCompP = new TH2F("pt_ZH_Comp_particle", "Sublead p_{T}^{ZH}", 50, 0.0, 400.0, 50, 0.0, 400.0);
  // eta parton level
  TH1 *hEtaH = new TH1F("eta_H", "#eta_{H}", 50, -3.0, 3.0);
  TH1 *hEtaZ = new TH1F("eta_Z", "#eta_{Z}", 50, -3.0, 3.0);
  TH1 *hEtaLE = new TH1F("eta_LE",  "Lead #eta_{e}", 50, -3.0, 3.0);
  TH1 *hEtaLM = new TH1F("eta_LM",  "Lead #eta_{#mu}", 50, -3.0, 3.0);
  TH1 *hEtaSLE = new TH1F("eta_SLE",  "Sublead #eta_{e}", 50, -3.0, 3.0);
  TH1 *hEtaSLM = new TH1F("eta_SLM",  "Sublead #eta_{#mu}", 50, -3.0, 3.0);
  TH1 *hEtaLB = new TH1F("eta_LB", "Lead #eta_{b}", 50, -3.0, 3.0);
  TH1 *hEtaSLB = new TH1F("eta_SLB", "Sublead #eta_{b}", 50, -3.0, 3.0);
  TH1 *hEtaZH = new TH1F("eta_ZH", "#eta_{ZH}", 50, -3.0, 3.0);
  TH1 *hDEtaEE = new TH1F("deta_EE", "#Delta#eta_{EE}", 50, -5.0, 5.0);
  TH1 *hDEtaMM = new TH1F("deta_MM", "#Delta#eta_{MM}", 50, -5.0, 5.0);
  // eta reco level
  TH1 *hEtaBBR = new TH1F("eta_bb_reco", "Reco #eta_{bb}", 50, -3.0, 3.0);
  TH1 *hEtaLLR = new TH1F("eta_ll_reco", "Reco #eta_{ll}", 50, -3.0, 3.0);
  TH1 *hEtaLER = new TH1F("eta_LE_reco",  "Lead Reco #eta_{e}", 50, -3.0, 3.0);
  TH1 *hEtaLMR = new TH1F("eta_LM_reco",  "Lead Reco #eta_{#mu}", 50, -3.0, 3.0);
  TH1 *hEtaSLER = new TH1F("eta_SLE_reco",  "Sublead Reco #eta_{e}", 50, -3.0, 3.0);
  TH1 *hEtaSLMR = new TH1F("eta_SLM_reco",  "Sublead Reco #eta_{#mu}", 50, -3.0, 3.0);
  TH1 *hEtaLBR = new TH1F("eta_LB_reco", "Lead Reco b Jet #eta", 50, -3.0, 3.0);
  TH1 *hEtaSLBR = new TH1F("eta_SLB_reco", "Sublead Reco b Jet #eta", 50, -3.0, 3.0);
  TH1 *hEtaBBLLR = new TH1F("eta_bbll_reco", "Reco #eta_{bbll}", 50, -3.0, 3.0);
  TH1 *hDEtaEER = new TH1F("deta_EE_reco", "Reco #Delta#eta_{EE}", 50, -5.0, 5.0);
  TH1 *hDEtaMMR = new TH1F("deta_MM_reco", "Reco #Delta#eta_{MM}", 50, -5.0, 5.0);
  // eta particle level
  TH1 *hEtaBBP = new TH1F("eta_bb_particle", "#eta_{bb}", 50, -3.0, 3.0);
  TH1 *hEtaLLP = new TH1F("eta_ll_particle", "#eta_{ll}", 50, -3.0, 3.0);
  TH1 *hEtaLEP = new TH1F("eta_LE_particle",  "Lead #eta_{e}", 50, -3.0, 3.0);
  TH1 *hEtaLMP = new TH1F("eta_LM_particle",  "Lead #eta_{#mu}", 50, -3.0, 3.0);
  TH1 *hEtaSLEP = new TH1F("eta_SLE_particle",  "Sublead #eta_{e}", 50, -3.0, 3.0);
  TH1 *hEtaSLMP = new TH1F("eta_SLM_particle",  "Sublead #eta_{#mu}", 50, -3.0, 3.0);
  TH1 *hEtaLBP = new TH1F("eta_LB_particle", "Lead b Jet #eta", 50, -3.0, 3.0);
  TH1 *hEtaSLBP = new TH1F("eta_SLB_particle", "Sublead b Jet #eta", 50, -3.0, 3.0);
  TH1 *hEtaBBLLP = new TH1F("eta_bbll_particle", "#eta_{bbll}", 50, -3.0, 3.0);
  TH1 *hDEtaEEP = new TH1F("deta_EE_particle", "#Delta#eta_{EE}", 50, -5.0, 5.0);
  TH1 *hDEtaMMP = new TH1F("deta_MM_particle", "#Delta#eta_{MM}", 50, -5.0, 5.0);
  // eta parton-reco comparison
  TH2 *hEtaHComp = new TH2F("eta_H_Comp", "#eta_{H}", 50, -3.0, 3.0, 50, -3.0, 3.0);
  TH2 *hEtaZComp = new TH2F("eta_Z_Comp", "#eta_{Z}", 50, -3.0, 3.0, 50, -3.0, 3.0);
  TH2 *hEtaLEComp = new TH2F("eta_LE_Comp", "Lead #eta_{e}", 50,-3.0, 3.0, 50, -3.0, 3.0);
  TH2 *hEtaLMComp = new TH2F("eta_LM_Comp", "Lead #eta_{#mu}", 50,-3.0, 3.0, 50, -3.0, 3.0);
  TH2 *hEtaSLEComp = new TH2F("eta_SLE_Comp", "Sublead #eta_{e}", 50,-3.0, 3.0, 50, -3.0, 3.0);
  TH2 *hEtaSLMComp = new TH2F("eta_SLM_Comp", "Sublead #eta_{#mu}", 50,-3.0, 3.0, 50, -3.0, 3.0);
  TH2 *hEtaLBComp = new TH2F("eta_LB_Comp", "Lead #eta_{b}", 50,-3.0, 3.0, 50, -3.0, 3.0);
  TH2 *hEtaSLBComp = new TH2F("eta_SLB_Comp", "Sublead #eta_{b}", 50,-3.0, 3.0, 50, -3.0, 3.0);
  TH2 *hEtaZHComp = new TH2F("eta_ZH_Comp", "Sublead #eta_{ZH}", 50, -3.0, 3.0, 50, -3.0, 3.0);
  TH2 *hDEtaEEComp = new TH2F("deta_EE_Comp", "#Delta#eta_{e}", 50, -5.0, 5.0, 50, -5.0, 5.0);
  TH2 *hDEtaMMComp = new TH2F("deta_MM_Comp", "#Delta#eta_{#mu}", 50, -5.0, 5.0, 50, -5.0, 5.0);
  // eta particle-reco comparison
  TH2 *hEtaHCompP = new TH2F("eta_H_Comp_particle", "#eta_{H}", 50, -3.0, 3.0, 50, -3.0, 3.0);
  TH2 *hEtaZCompP = new TH2F("eta_Z_Comp_particle", "#eta_{Z}", 50, -3.0, 3.0, 50, -3.0, 3.0);
  TH2 *hEtaLECompP = new TH2F("eta_LE_Comp_particle", "Lead #eta_{e}", 50,-3.0, 3.0, 50, -3.0, 3.0);
  TH2 *hEtaLMCompP = new TH2F("eta_LM_Comp_particle", "Lead #eta_{#mu}", 50,-3.0, 3.0, 50, -3.0, 3.0);
  TH2 *hEtaSLECompP = new TH2F("eta_SLE_Comp_particle", "Sublead #eta_{e}", 50,-3.0, 3.0, 50, -3.0, 3.0);
  TH2 *hEtaSLMCompP = new TH2F("eta_SLM_Comp_particle", "Sublead #eta_{#mu}", 50,-3.0, 3.0, 50, -3.0, 3.0);
  TH2 *hEtaLBCompP = new TH2F("eta_LB_Comp_particle", "Lead #eta_{b}", 50,-3.0, 3.0, 50, -3.0, 3.0);
  TH2 *hEtaSLBCompP = new TH2F("eta_SLB_Comp_particle", "Sublead #eta_{b}", 50,-3.0, 3.0, 50, -3.0, 3.0);
  TH2 *hEtaZHCompP = new TH2F("eta_ZH_Comp_particle", "Sublead #eta_{ZH}", 50, -3.0, 3.0, 50, -3.0, 3.0);
  TH2 *hDEtaEECompP = new TH2F("deta_EE_Comp_particle", "#Delta#eta_{e}", 50, -5.0, 5.0, 50, -5.0, 5.0);
  TH2 *hDEtaMMCompP = new TH2F("deta_MM_Comp_particle", "#Delta#eta_{#mu}", 50, -5.0, 5.0, 50, -5.0, 5.0);
  // cos theta parton 1evel
  TH1 *hCosThetaLE = new TH1F("costheta_LE",  "Lead cos#theta_{e}", 50, -1.0, 1.0);
  TH1 *hCosThetaLM = new TH1F("costheta_LM",  "Lead cos#theta_{#mu}", 50, -1.0, 1.0);
  TH1 *hCosThetaSLE = new TH1F("costheta_SLE",  "Sublead cos#theta_{e}", 50, -1.0, 1.0);
  TH1 *hCosThetaSLM = new TH1F("costheta_SLM",  "Sublead cos#theta_{#mu}", 50, -1.0, 1.0);
  TH1 *hCosThetaLB = new TH1F("costheta_LB", "Lead cos#theta_{b}", 50, -1.0, 1.0);
  TH1 *hCosThetaSLB = new TH1F("costheta_SLB", "Sublead cos#theta_{b}", 50, -1.0, 1.0);
  // cos theta reco 1evel
  TH1 *hCosThetaLER = new TH1F("costheta_LE_reco",  "Lead Reco cos#theta_{e}", 50, -1.0, 1.0);
  TH1 *hCosThetaLMR = new TH1F("costheta_LM_reco",  "Lead Reco cos#theta_{#mu}", 50, -1.0, 1.0);
  TH1 *hCosThetaSLER = new TH1F("costheta_SLE_reco",  "Sublead Reco cos#theta_{e}", 50, -1.0, 1.0);
  TH1 *hCosThetaSLMR = new TH1F("costheta_SLM_reco",  "Sublead Reco cos#theta_{#mu}", 50, -1.0, 1.0);
  TH1 *hCosThetaLBR = new TH1F("costheta_LB_reco", "Lead Reco b Jet cos#theta", 50, -1.0, 1.0);
  TH1 *hCosThetaSLBR = new TH1F("costheta_SLB_reco", "Sublead Reco b Jet cos#theta", 50, -1.0, 1.0);
  // cos theta particle 1evel
  TH1 *hCosThetaLEP = new TH1F("costheta_LE_particle",  "Lead cos#theta_{e}", 50, -1.0, 1.0);
  TH1 *hCosThetaLMP = new TH1F("costheta_LM_particle",  "Lead cos#theta_{#mu}", 50, -1.0, 1.0);
  TH1 *hCosThetaSLEP = new TH1F("costheta_SLE_particle",  "Sublead cos#theta_{e}", 50, -1.0, 1.0);
  TH1 *hCosThetaSLMP = new TH1F("costheta_SLM_particle",  "Sublead cos#theta_{#mu}", 50, -1.0, 1.0);
  TH1 *hCosThetaLBP = new TH1F("costheta_LB_particle", "Lead b Jet cos#theta", 50, -1.0, 1.0);
  TH1 *hCosThetaSLBP = new TH1F("costheta_SLB_particle", "Sublead b Jet cos#theta", 50, -1.0, 1.0);
  // cos theta parton-reco comparison
  TH2 *hCosThetaLEComp = new TH2F("costheta_LE_Comp", "Lead cos#theta_{e}", 50,-1.0, 1.0, 50, -1.0, 1.0);
  TH2 *hCosThetaLMComp = new TH2F("costheta_LM_Comp", "Lead cos#theta_{#mu}", 50,-1.0, 1.0, 50, -1.0, 1.0);
  TH2 *hCosThetaSLEComp = new TH2F("costheta_SLE_Comp", "Sublead cos#theta_{e}", 50,-1.0, 1.0, 50, -1.0, 1.0);
  TH2 *hCosThetaSLMComp = new TH2F("costheta_SLM_Comp", "Sublead cos#theta_{#mu}", 50,-1.0, 1.0, 50, -1.0, 1.0);
  TH2 *hCosThetaLBComp = new TH2F("costheta_LB_Comp", "Lead cos#theta_{b}", 50,-1.0, 1.0, 50, -1.0, 1.0);
  TH2 *hCosThetaSLBComp = new TH2F("costheta_SLB_Comp", "Sublead cos#theta_{b}", 50,-1.0, 1.0, 50, -1.0, 1.0);
  // cos theta particle-reco comparison
  TH2 *hCosThetaLECompP = new TH2F("costheta_LE_Comp_particle", "Lead cos#theta_{e}", 50,-1.0, 1.0, 50, -1.0, 1.0);
  TH2 *hCosThetaLMCompP = new TH2F("costheta_LM_Comp_particle", "Lead cos#theta_{#mu}", 50,-1.0, 1.0, 50, -1.0, 1.0);
  TH2 *hCosThetaSLECompP = new TH2F("costheta_SLE_Comp_particle", "Sublead cos#theta_{e}", 50,-1.0, 1.0, 50, -1.0, 1.0);
  TH2 *hCosThetaSLMCompP = new TH2F("costheta_SLM_Comp_particle", "Sublead cos#theta_{#mu}", 50,-1.0, 1.0, 50, -1.0, 1.0);
  TH2 *hCosThetaLBCompP = new TH2F("costheta_LB_Comp_particle", "Lead cos#theta_{b}", 50,-1.0, 1.0, 50, -1.0, 1.0);
  TH2 *hCosThetaSLBCompP = new TH2F("costheta_SLB_Comp_particle", "Sublead cos#theta_{b}", 50,-1.0, 1.0, 50, -1.0, 1.0);
  // phi parton level
  TH1 *hPhiH = new TH1F("phi_H", "#phi_{H}", 50, -3.15, 3.15);
  TH1 *hPhiZ = new TH1F("phi_Z", "#phi_{Z}", 50, -3.15, 3.15);
  TH1 *hPhiLE = new TH1F("phi_LE",  "Lead #phi_{e}", 50, -3.15, 3.15);
  TH1 *hPhiLM = new TH1F("phi_LM",  "Lead #phi_{#mu}", 50, -3.15, 3.15);
  TH1 *hPhiSLE = new TH1F("phi_SLE",  "Sublead #phi_{e}", 50, -3.15, 3.15);
  TH1 *hPhiSLM = new TH1F("phi_SLM",  "Sublead #phi_{#mu}", 50, -3.15, 3.15);
  TH1 *hPhiLB = new TH1F("phi_LB", "Lead #phi_{b}", 50, -3.15, 3.15);
  TH1 *hPhiSLB = new TH1F("phi_SLB", "Sublead #phi_{b}", 50, -3.15, 3.15);
  TH1 *hPhiZH = new TH1F("phi_ZH", "#phi_{ZH}", 50, -3.15, 3.15);
  TH1 *hDPhiEE = new TH1F("dphi_EE", "#Delta#phi_{EE}", 50, -3.15, 3.15);
  TH1 *hDPhiMM = new TH1F("dphi_MM", "#Delta#phi_{MM}", 50, -3.15, 3.15);
  TH1 *hDPhiZH = new TH1F("dphi_ZH", "#Delta#phi_{ZH}", 50, -3.15, 3.15);
  // phi reco level
  TH1 *hPhiBBR = new TH1F("phi_bb_reco", "Reco #phi_{bb}", 50, -3.15, 3.15);
  TH1 *hPhiLLR = new TH1F("phi_ll_reco", "Reco #phi_{ll}", 50, -3.15, 3.15);
  TH1 *hPhiLER = new TH1F("phi_LE_reco",  "Lead Reco #phi_{e}", 50, -3.15, 3.15);
  TH1 *hPhiLMR = new TH1F("phi_LM_reco",  "Lead Reco #phi_{#mu}", 50, -3.15, 3.15);
  TH1 *hPhiSLER = new TH1F("phi_SLE_reco",  "Sublead Reco #phi_{e}", 50, -3.15, 3.15);
  TH1 *hPhiSLMR = new TH1F("phi_SLM_reco",  "Sublead Reco #phi_{#mu}", 50, -3.15, 3.15);
  TH1 *hPhiLBR = new TH1F("phi_LB_reco", "Lead Reco b Jet #phi", 50, -3.15, 3.15);
  TH1 *hPhiSLBR = new TH1F("phi_SLB_reco", "Sublead Reco b Jet #phi", 50, -3.15, 3.15);
  TH1 *hPhiBBLLR = new TH1F("phi_bbll_reco", "Reco #phi_{bbll}", 50, -3.15, 3.15);
  TH1 *hDPhiEER = new TH1F("dphi_EE_reco", "Reco #Delta#phi_{EE}", 50, -3.15, 3.15);
  TH1 *hDPhiMMR = new TH1F("dphi_MM_reco", "Reco #Delta#phi_{MM}", 50, -3.15, 3.15);
  TH1 *hDPhiZHR = new TH1F("dphi_ZH_reco", "Reco #Delta#phi_{bb,ll}", 50, -3.15, 3.15);
  // phi particle level
  TH1 *hPhiBBP = new TH1F("phi_bb_particle", "#phi_{bb}", 50, -3.15, 3.15);
  TH1 *hPhiLLP = new TH1F("phi_ll_particle", "#phi_{ll}", 50, -3.15, 3.15);
  TH1 *hPhiLEP = new TH1F("phi_LE_particle",  "Lead #phi_{e}", 50, -3.15, 3.15);
  TH1 *hPhiLMP = new TH1F("phi_LM_particle",  "Lead #phi_{#mu}", 50, -3.15, 3.15);
  TH1 *hPhiSLEP = new TH1F("phi_SLE_particle",  "Sublead #phi_{e}", 50, -3.15, 3.15);
  TH1 *hPhiSLMP = new TH1F("phi_SLM_particle",  "Sublead #phi_{#mu}", 50, -3.15, 3.15);
  TH1 *hPhiLBP = new TH1F("phi_LB_particle", "Lead b Jet #phi", 50, -3.15, 3.15);
  TH1 *hPhiSLBP = new TH1F("phi_SLB_particle", "Sublead b Jet #phi", 50, -3.15, 3.15);
  TH1 *hPhiBBLLP = new TH1F("phi_bbll_particle", "#phi_{bbll}", 50, -3.15, 3.15);
  TH1 *hDPhiEEP = new TH1F("dphi_EE_particle", "#Delta#phi_{EE}", 50, -3.15, 3.15);
  TH1 *hDPhiMMP = new TH1F("dphi_MM_particle", "#Delta#phi_{MM}", 50, -3.15, 3.15);
  TH1 *hDPhiZHP = new TH1F("dphi_ZH_particle", "#Delta#phi_{bb,ll}", 50, -3.15, 3.15);
  // phi comparison
  TH2 *hPhiHComp = new TH2F("phi_H_Comp", "#phi_{H}", 50, -3.15, 3.15, 50, -3.15, 3.15);
  TH2 *hPhiZComp = new TH2F("phi_Z_Comp", "#phi_{Z}", 50, -3.15, 3.15, 50, -3.15, 3.15);
  TH2 *hPhiLEComp = new TH2F("phi_LE_Comp", "Lead #phi_{e}", 50,-3.15, 3.15, 50, -3.15, 3.15);
  TH2 *hPhiLMComp = new TH2F("phi_LM_Comp", "Lead #phi_{#mu}", 50,-3.15, 3.15, 50, -3.15, 3.15);
  TH2 *hPhiSLEComp = new TH2F("phi_SLE_Comp", "Sublead #phi_{e}", 50,-3.15, 3.15, 50, -3.15, 3.15);
  TH2 *hPhiSLMComp = new TH2F("phi_SLM_Comp", "Sublead #phi_{#mu}", 50,-3.15, 3.15, 50, -3.15, 3.15);
  TH2 *hPhiLBComp = new TH2F("phi_LB_Comp", "Lead #phi_{b}", 50,-3.15, 3.15, 50, -3.15, 3.15);
  TH2 *hPhiSLBComp = new TH2F("phi_SLB_Comp", "Sublead #phi_{b}", 50,-3.15, 3.15, 50, -3.15, 3.15);
  TH2 *hPhiZHComp = new TH2F("phi_ZH_Comp", "Sublead #phi_{ZH}", 50, -3.15, 3.15, 50, -3.15, 3.15);
  TH2 *hDPhiEEComp = new TH2F("dphi_EE_Comp", "#Delta#phi_{e}", 50, 0, 3.15, 50, 0, 3.15);
  TH2 *hDPhiMMComp = new TH2F("dphi_MM_Comp", "#Delta#phi_{#mu}", 50, 0, 3.15, 50, 0, 3.15);
  TH2 *hDPhiEECompP = new TH2F("dphi_EE_Comp_particle", "#Delta#phi_{e}", 50, -3.15, 3.15, 50, -3.15, 3.15);
  TH2 *hDPhiMMCompP = new TH2F("dphi_MM_Comp_particle", "#Delta#phi_{#mu}", 50, -3.15, 3.15, 50, -3.15, 3.15);
  // initialize variables needed for filling histograms
  // info for totals
  double  nPassed=0;
  int  nPassedRaw=0;
  // variables to keep track of event content
  int bjets = 0;
  int num_elec_reco;
  int num_mu_reco;
  int num_elec_particle;
  int num_mu_particle;
  bool elec_ev_parton;
  // vectors to store event constituents
  TLorentzVector b1_parton;
  TLorentzVector b2_parton;
  TLorentzVector b1_reco;
  TLorentzVector b2_reco;
  TLorentzVector b1_particle;
  TLorentzVector b2_particle;
  TLorentzVector e1_parton;
  TLorentzVector e2_parton;
  TLorentzVector e1_particle;
  TLorentzVector e2_particle;
  TLorentzVector m1_parton;
  TLorentzVector m2_parton;
  TLorentzVector m1_particle;
  TLorentzVector m2_particle;
  TLorentzVector higgsvec;
  TLorentzVector elecvec1;
  TLorentzVector elecvec2;
  TLorentzVector muvec1;
  TLorentzVector muvec2;
  TLorentzVector zvec;
  TLorentzVector sysvec;
  TLorentzVector zpartonvec;
  TLorentzVector zparticlevec;
  TLorentzVector hpartonvec;
  TLorentzVector hparticlevec;
  // particle pointers
  Electron *elec;
  Electron *elec1;
  Electron *elec2;
  Muon *muon;
  Muon *muon1;
  Muon *muon2;
  GenParticle *daughter1;
  GenParticle *daughter2;
  // calculate total weight
  for(Int_t entry = 0; entry < numberOfEntries; ++entry){
    // Load selected branches with data from specified event
    treeReader->ReadEntry(entry);
    HepMCEvent *event = (HepMCEvent*) branchEvent -> At(0);
    totalWeight += event->Weight;
  }
  // Loop over events
  for(Int_t entry = 0; entry < numberOfEntries; ++entry){
    // Load selected branches with data from specified event
    treeReader->ReadEntry(entry);
    HepMCEvent *event = (HepMCEvent*) branchEvent -> At(0);
    Float_t weight = event->Weight*Lumi*cross_section*numberOfEntries/(numEntries*totalWeight);
    Float_t test_weight = event->Weight*cross_section*numberOfEntries/(numEntries*totalWeight);
    hWeight -> Fill(event->Weight, test_weight);
    bool fill_reco = true;
    bool fill_parton = true;
    bool fill_particle = true;
    // Loop over jets in event and save the b jets in 2b events as lorentz vectors to reconstruct the higgs
    bjets = 0;
    for(int i=0; i<(int)branchJet->GetEntries(); i++){
      Jet *jet=(Jet*) branchJet->At(i);
      if (isMyBTag(jet, branchGenParticle) && abs(jet->Eta) < eta_cut && jet->PT > jet_pt_cut) {
        bjets += 1;
        if (bjets == 1) {
          b1_reco = jet->P4();
        } else if (bjets == 2){
          b2_reco = jet->P4();
          break;
        }
      }
    }
    if (bjets != 2) fill_reco = false;
    higgsvec = b1_reco + b2_reco;
    if (higgsvec.M() > higgs_mass_cut_hi || higgsvec.M() < higgs_mass_cut_low || higgsvec.Pt() < higgs_pt_cut) {
        fill_reco = false;
    }
    // Check for two electrons meeting requirements 
    num_elec_reco = 0;
    num_mu_reco = 0;
    if (fill_reco) {
      if (branchElectron->GetEntries() > 1) {
        for(int i=0; i<(int)branchElectron->GetEntries(); i++) {
          elec = (Electron *) branchElectron->At(i);
          if (abs(elec->Eta) < eta_cut && num_elec_reco == 0 && elec->PT > e_pt_cut_lead) {
            elec1 = elec;
            num_elec_reco += 1;
          } else if (abs(elec->Eta) < eta_cut && num_elec_reco == 1 && elec->PT > e_pt_cut_sub) {
            elec2 = elec;
            num_elec_reco += 1;
            if ((elec1->Charge + elec2->Charge) != 0) {
              fill_reco = false;
              break;
            }
          } else if (abs(elec->Eta) < eta_cut && num_elec_reco == 2 && elec->PT > e_pt_cut_sub) {
            fill_reco = false;
            break;
          }
        }
      }
      // Check for two muons meeting requirements 
      if (branchMuon->GetEntries() > 1) {
        for(int i=0; i<(int)branchMuon->GetEntries(); i++) {
          muon = (Muon *) branchMuon->At(i);
          if (abs(muon->Eta) < eta_cut && num_mu_reco == 0 && muon->PT > mu_pt_cut_lead) {
            muon1 = muon;
            num_mu_reco += 1;
          } else if (abs(muon->Eta) < eta_cut && num_mu_reco == 1 && muon->PT > mu_pt_cut_sub) {
            muon2 = muon;
            num_mu_reco += 1;
            if ((muon1->Charge + muon2->Charge) != 0) {
              fill_reco = false;
              break;
            }
          } else if (abs(muon->Eta) < eta_cut && num_mu_reco == 2 && muon->PT > mu_pt_cut_sub) {
            fill_reco = false;
            break;
          }
        }
      }
    }
    // fill electron data
    if (fill_reco && num_elec_reco == 2 && num_mu_reco == 0) {
      elecvec1 = elec1->P4();
      elecvec2 = elec2->P4();
      zvec = elecvec1 + elecvec2;

      if (zvec.M() > z_mass_cut_hi || zvec.M() < z_mass_cut_low || zvec.Pt() < z_pt_cut || abs(deltaPhi(zvec.Phi(), higgsvec.Phi())) < dphi_zh_cut) {
        fill_reco = false;
      } else {
        hPtLER -> Fill(elec1->PT, weight);
        hEtaLER -> Fill(elec1->Eta, weight);
        hPhiLER -> Fill(elec1->Phi, weight);
        hPtSLER -> Fill(elec2->PT, weight);
        hEtaSLER -> Fill(elec2->Eta, weight);
        hPhiSLER -> Fill(elec2->Phi, weight);
        hDPhiEER -> Fill(deltaPhi(elec1->Phi, elec2->Phi), weight); 
        hDEtaEER -> Fill(elec1->Eta - elec2->Eta, weight); 
        elecvec1.Boost(-zvec.BoostVector());
        elecvec2.Boost(-zvec.BoostVector());
        hCosThetaLER -> Fill(elecvec1.CosTheta(), weight);
        hCosThetaSLER -> Fill(elecvec2.CosTheta(), weight);
        elecvec1.Boost(zvec.BoostVector());
        elecvec2.Boost(zvec.BoostVector());
      }
    // fill muon data
    } else if (fill_reco && num_elec_reco == 0 && num_mu_reco == 2) {
      muvec1 = muon1->P4();
      muvec2 = muon2->P4();
      zvec = muvec1 + muvec2;
      if (zvec.M() > z_mass_cut_hi || zvec.M() < z_mass_cut_low || zvec.Pt() < z_pt_cut || abs(deltaPhi(zvec.Phi(), higgsvec.Phi())) < dphi_zh_cut) {
        fill_reco = false;
      } else {
        hPtLMR -> Fill(muon1->PT, weight);
        hEtaLMR -> Fill(muon1->Eta, weight);
        hPhiLMR -> Fill(muon1->Phi, weight);
        hPtSLMR -> Fill(muon2->PT, weight);
        hEtaSLMR -> Fill(muon2->Eta, weight);
        hPhiSLMR -> Fill(muon2->Phi, weight);
        hDPhiMMR -> Fill(deltaPhi(muon1->Phi, muon2->Phi), weight);
        hDEtaMMR -> Fill(muon1->Eta - muon2->Eta, weight);
        muvec1.Boost(-zvec.BoostVector());
        muvec2.Boost(-zvec.BoostVector());
        hCosThetaLMR -> Fill(muvec1.CosTheta(), weight);
        hCosThetaSLMR -> Fill(muvec2.CosTheta(), weight);
        muvec1.Boost(zvec.BoostVector());
        muvec2.Boost(zvec.BoostVector());
      } 
    } else {
      fill_reco = false;
    }
    // fill mass reco histograms
    sysvec = higgsvec + zvec;
    if (fill_reco) {
      hMllR -> Fill(zvec.M(), weight);
      hMbbllR -> Fill(sysvec.M(), weight);
      hMbbR -> Fill(higgsvec.M(), weight);
      // fill pt reco histograms 
      hPtBBR -> Fill(higgsvec.Pt(), weight);
      hPtLLR -> Fill(zvec.Pt(), weight);
      hPtLBR -> Fill(b1_reco.Pt(), weight);
      hPtSLBR -> Fill(b2_reco.Pt(), weight);
      hPtBBLLR -> Fill(sysvec.Pt(), weight);
      // fill eta reco histograms
      hEtaBBR -> Fill(higgsvec.Eta(), weight);
      hEtaLLR -> Fill(zvec.Eta(), weight);
      hEtaLBR -> Fill(b1_reco.Eta(), weight);
      hEtaSLBR -> Fill(b2_reco.Eta(), weight);
      hEtaBBLLR -> Fill(sysvec.Eta(), weight);
      // fill phi reco histograms
      hPhiBBR -> Fill(higgsvec.Phi(), weight);
      hPhiLLR -> Fill(zvec.Phi(), weight);
      hPhiLBR -> Fill(b1_reco.Phi(), weight);
      hPhiSLBR -> Fill(b2_reco.Phi(), weight);
      hPhiBBLLR -> Fill(sysvec.Phi(), weight);
      hDPhiZHR -> Fill(deltaPhi(zvec.Phi(), higgsvec.Phi()), weight);
      // fill cos theta reco histograms
      b1_reco.Boost(-higgsvec.BoostVector());
      b2_reco.Boost(-higgsvec.BoostVector());
      hCosThetaLBR -> Fill(b1_reco.CosTheta(), weight);
      hCosThetaSLBR -> Fill(b2_reco.CosTheta(), weight);
      b1_reco.Boost(higgsvec.BoostVector());
      b2_reco.Boost(higgsvec.BoostVector());
    }
    // loop over true particles and fill those histograms
    b1_parton.SetPtEtaPhiM(0,0,0,0);
    b2_parton.SetPtEtaPhiM(0,0,0,0);
    e1_parton.SetPtEtaPhiM(0,0,0,0);
    e2_parton.SetPtEtaPhiM(0,0,0,0);
    m1_parton.SetPtEtaPhiM(0,0,0,0);
    m2_parton.SetPtEtaPhiM(0,0,0,0);
    e1_particle.SetPtEtaPhiM(0,0,0,0);
    e2_particle.SetPtEtaPhiM(0,0,0,0);
    m1_particle.SetPtEtaPhiM(0,0,0,0);
    m2_particle.SetPtEtaPhiM(0,0,0,0);
    hpartonvec.SetPtEtaPhiM(0,0,0,0);
    zpartonvec.SetPtEtaPhiM(0,0,0,0);
    num_elec_particle = 0;
    num_mu_particle = 0;
    for(int i=0; i<(int)branchGenParticle->GetEntries(); i++){
      // safely access mother and daughter particles
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
      // Higgs parton
      if (particle->PID == 25 && d1_pid != 25 && d2_pid != 25) {
        hpartonvec = particle->P4();
        // check for b parton children
        if (abs(d1_pid) == 5 && abs(d2_pid) == 5) {
          if (daughter1 -> PT > daughter2 -> PT) {
            b1_parton = daughter1 -> P4();
            b2_parton = daughter2 -> P4();
          } else {
            b1_parton = daughter2 -> P4();
            b2_parton = daughter1 -> P4();
          }
          if (abs(b1_parton.Eta()) > eta_cut || abs(b2_parton.Eta()) > eta_cut) fill_parton = false;
          if (b1_parton.Pt() < jet_pt_cut || b2_parton.Pt() < jet_pt_cut) fill_parton = false;
        } else fill_parton = false;
      // Z parton
      } else if ((particle->PID == 23) && (d1_pid != 23) && (d2_pid != 23)) {
        zpartonvec = particle->P4();
        // check for electron daughters
        if (abs(d1_pid) == 11 && abs(d2_pid) == 11) {
          if (daughter1 -> PT > daughter2 -> PT) {
            e1_parton = daughter1 -> P4();
            e2_parton = daughter2 -> P4();
          } else {
            e1_parton = daughter2 -> P4();
            e2_parton = daughter1 -> P4();
          }
          elec_ev_parton = true;
          if (abs(e1_parton.Eta()) > eta_cut || abs(e2_parton.Eta()) > eta_cut) {
            fill_parton = false;
          } else if (e1_parton.Pt() < e_pt_cut_lead || e1_parton.Pt() < e_pt_cut_sub) {
            fill_parton = false;
          }
        // check for muon daughters
        } else if (abs(d1_pid) == 13 && abs(d2_pid) == 13) {
          if (daughter1 -> PT > daughter2 -> PT) {
            m1_parton = daughter1 -> P4();
            m2_parton = daughter2 -> P4();
          } else {
            m1_parton = daughter2 -> P4();
            m2_parton = daughter1 -> P4();
          }
          elec_ev_parton = false;
          if (abs(m1_parton.Eta()) > eta_cut || abs(m2_parton.Eta()) > eta_cut) {
            fill_parton = false;
          } else if (m1_parton.Pt() < mu_pt_cut_lead || m1_parton.Pt() < mu_pt_cut_sub) {
            fill_parton = false;
          }
        } else fill_parton = false;
      // now look for status 1 (particle level) electrons
      } else if (abs(particle -> PID) == 11 && particle -> Status == 1 && abs(particle -> Eta) < eta_cut) {
        if (num_elec_particle > 1 && particle->PT > e_pt_cut_sub) {
          fill_particle = false;
          continue;
        }
        if (particle->PT > e_pt_cut_lead && particle->PT > e1_particle.Pt()) {
          num_elec_particle += 1;
          e2_particle.SetPtEtaPhiM(e1_particle.Pt(), e1_particle.Eta(), e1_particle.Phi(), e1_particle.M());
          e1_particle = particle -> P4();
        } else if (particle->PT > e_pt_cut_sub) {
          num_elec_particle += 1;
          e2_particle = particle -> P4();
        }
      // now look for status 1 (particle level) muons
      } else if (abs(particle -> PID) == 13 && particle -> Status == 1 && abs(particle -> Eta) < eta_cut) {
        if (num_mu_particle > 1 && particle->PT > mu_pt_cut_sub) {
          fill_particle = false;
          continue;
        }
        if (particle->PT > mu_pt_cut_lead && particle -> PT > m1_particle.Pt()) {
          num_mu_particle += 1;
          m2_particle.SetPtEtaPhiM(m1_particle.Pt(), m1_particle.Eta(), m1_particle.Phi(), m1_particle.M());
          m1_particle = particle -> P4();
        } else if (particle->PT > mu_pt_cut_sub) {
          num_mu_particle += 1;
          m2_particle = particle -> P4();
        }
      }
    }
    if (hpartonvec.Pt() < higgs_pt_cut || zpartonvec.Pt() < z_pt_cut || hpartonvec.Pt() < higgs_pt_cut  || abs(deltaPhi(zpartonvec.Phi(), hpartonvec.Phi())) < dphi_zh_cut) fill_parton = false;
    if (fill_parton) {
      // fill electron parton histograms
      if (elec_ev_parton) {
        hPtLE -> Fill(e1_parton.Pt(), weight);
        hEtaLE -> Fill(e1_parton.Eta(), weight);
        hPhiLE -> Fill(e1_parton.Phi(), weight);
        hPtSLE -> Fill(e2_parton.Pt(), weight);
        hEtaSLE -> Fill(e2_parton.Eta(), weight);
        hPhiSLE -> Fill(e2_parton.Phi(), weight);
        hDPhiEE -> Fill(deltaPhi(e1_parton.Phi(), e2_parton.Phi()), weight);
        hDEtaEE -> Fill(e1_parton.Eta()-e2_parton.Eta(), weight);
        // boost to z frame for cos theta
        e1_parton.Boost(-zpartonvec.BoostVector());
        e2_parton.Boost(-zpartonvec.BoostVector());
        hCosThetaLE -> Fill(e1_parton.CosTheta(), weight);
        hCosThetaSLE -> Fill(e2_parton.CosTheta(), weight);
        // fill cos theta comparison histograms (boosting back to z frame for reco as well)
        if (fill_reco && num_elec_reco == 2) {
          elecvec1.Boost(-zvec.BoostVector());
          elecvec2.Boost(-zvec.BoostVector());
          hCosThetaLEComp -> Fill(e1_parton.CosTheta(), elecvec1.CosTheta(), weight);
          hCosThetaSLEComp -> Fill(e2_parton.CosTheta(), elecvec2.CosTheta(), weight);
          // boost back to lab frame
          elecvec1.Boost(zvec.BoostVector());
          elecvec2.Boost(zvec.BoostVector());
        }
        e1_parton.Boost(zpartonvec.BoostVector());
        e2_parton.Boost(zpartonvec.BoostVector());
        // fill other electron comparison histograms
        if (fill_reco && num_elec_reco == 2) {
          hPtLEComp -> Fill(e1_parton.Pt(), elecvec1.Pt(), weight);
          hEtaLEComp -> Fill(e1_parton.Eta(), elecvec1.Eta(), weight);
          hPhiLEComp -> Fill(e1_parton.Phi(), elecvec1.Phi(), weight);
          hPtSLEComp -> Fill(e2_parton.Pt(), elecvec2.Pt(), weight);
          hEtaSLEComp -> Fill(e2_parton.Eta(), elecvec2.Eta(), weight);
          hPhiSLEComp -> Fill(e2_parton.Phi(), elecvec2.Phi(), weight);
          hDPhiEEComp -> Fill(deltaPhi(e1_parton.Phi(), e2_parton.Phi()), deltaPhi(elecvec1.Phi(), elecvec2.Phi()), weight);
          hDEtaEEComp -> Fill(e1_parton.Eta()-e2_parton.Eta(), elecvec1.Eta()-elecvec2.Eta(), weight);
        }
      } else {
        // fill muon parton histograms
        hPtLM -> Fill(m1_parton.Pt(), weight);
        hEtaLM -> Fill(m1_parton.Eta(), weight);
        hPhiLM -> Fill(m1_parton.Phi(), weight);
        hPtSLM -> Fill(m2_parton.Pt(), weight);
        hEtaSLM -> Fill(m2_parton.Eta(), weight);
        hPhiSLM -> Fill(m2_parton.Phi(), weight);
        hDPhiMM -> Fill(deltaPhi(m1_parton.Phi(), m2_parton.Phi()), weight);
        hDEtaMM -> Fill(m1_parton.Eta()-m2_parton.Eta(), weight);
        // boost to z frame for cos theta
        m1_parton.Boost(-zpartonvec.BoostVector());
        m2_parton.Boost(-zpartonvec.BoostVector());
        hCosThetaLM -> Fill(m1_parton.CosTheta(), weight);
        hCosThetaSLM -> Fill(m2_parton.CosTheta(), weight);
        // fill cos theta comparison histograms (boosting to z frame for reco as well)
        if (fill_reco && num_mu_reco == 2) {
          muvec1.Boost(-zvec.BoostVector());
          muvec2.Boost(-zvec.BoostVector());
          hCosThetaLMComp -> Fill(m1_parton.CosTheta(), muvec1.CosTheta(), weight);
          hCosThetaSLMComp -> Fill(m2_parton.CosTheta(), muvec2.CosTheta(), weight);
          //boost back to lab frame
          muvec1.Boost(zvec.BoostVector());
          muvec2.Boost(zvec.BoostVector());
        }
        m1_parton.Boost(zpartonvec.BoostVector());
        m2_parton.Boost(zpartonvec.BoostVector());
        // fill other muon comparison histograms
        if (fill_reco && num_mu_reco == 2) {
          hPtSLMComp -> Fill(m2_parton.Pt(), muvec2.Pt(), weight);
          hEtaSLMComp -> Fill(m2_parton.Eta(), muvec2.Eta(), weight);
          hPhiSLMComp -> Fill(m2_parton.Phi(), muvec2.Phi(), weight);
          hPtLMComp -> Fill(m1_parton.Pt(), muvec1.Pt(), weight);
          hEtaLMComp -> Fill(m1_parton.Eta(), muvec1.Eta(), weight);
          hPhiLMComp -> Fill(m1_parton.Phi(), muvec1.Phi(), weight);
          hDPhiMMComp -> Fill(deltaPhi(m1_parton.Phi(), m2_parton.Phi()), deltaPhi(muvec1.Phi(), muvec2.Phi()), weight);
          hDEtaMMComp -> Fill(m1_parton.Eta()-m2_parton.Eta(), muvec1.Eta()-muvec2.Eta(), weight);
        }
      }
      // fill higgs and z parton histograms
      hPtH -> Fill(hpartonvec.Pt(), weight);
      hEtaH -> Fill(hpartonvec.Eta(), weight);
      hPhiH -> Fill(hpartonvec.Phi(), weight);
      hPtZ -> Fill(zpartonvec.Pt(), weight);
      hEtaZ -> Fill(zpartonvec.Eta(), weight);
      hPhiZ -> Fill(zpartonvec.Phi(), weight);
      if (fill_reco) {
        // fill higgs and z comparison histograms
        hPtHComp -> Fill(hpartonvec.Pt(), higgsvec.Pt(), weight);
        hEtaHComp -> Fill(hpartonvec.Eta(), higgsvec.Eta(), weight);
        hPhiHComp -> Fill(hpartonvec.Phi(), higgsvec.Phi(), weight);
        hPtZComp -> Fill(zpartonvec.Pt(), zvec.Pt(), weight);
        hEtaZComp -> Fill(zpartonvec.Eta(), zvec.Eta(), weight);
        hPhiZComp -> Fill(zpartonvec.Phi(), zvec.Phi(), weight);
      }
      // fill b parton histograms
      hPtLB -> Fill(b1_parton.Pt(), weight);
      hPtSLB -> Fill(b2_parton.Pt(), weight);
      hEtaLB -> Fill(b1_parton.Eta(), weight);
      hEtaSLB -> Fill(b2_parton.Eta(), weight);
      hPhiLB -> Fill(b1_parton.Phi(), weight);
      hPhiSLB -> Fill(b2_parton.Phi(), weight);
      // boost to higgs frame for cos theta
      b1_parton.Boost(-hpartonvec.BoostVector());
      b2_parton.Boost(-hpartonvec.BoostVector());
      hCosThetaLB -> Fill(b1_parton.CosTheta(), weight);
      hCosThetaSLB -> Fill(b2_parton.CosTheta(), weight);
      // fill cos theta comparison histograms (boosting to higgs frame for reco as well)
      if (fill_reco) {
        b1_reco.Boost(-higgsvec.BoostVector());
        b2_reco.Boost(-higgsvec.BoostVector());
        hCosThetaLMComp -> Fill(b1_parton.CosTheta(), b1_reco.CosTheta(), weight);
        hCosThetaSLMComp -> Fill(b2_parton.CosTheta(), b2_reco.CosTheta(), weight);
        // boost back to lab frame
        b1_reco.Boost(higgsvec.BoostVector());
        b2_reco.Boost(higgsvec.BoostVector());
      }
      b1_parton.Boost(hpartonvec.BoostVector());
      b2_parton.Boost(hpartonvec.BoostVector());
      // fill other b comparison histograms
      if (fill_reco) {
        hPtLBComp -> Fill(b1_parton.Pt(), b1_reco.Pt(), weight);
        hPtSLBComp -> Fill(b2_parton.Pt(), b2_reco.Pt(), weight);
        hEtaLBComp -> Fill(b1_parton.Eta(), b1_reco.Eta(), weight);
        hEtaSLBComp -> Fill(b2_parton.Eta(), b2_reco.Eta(), weight);
        hPhiLBComp -> Fill(b1_parton.Phi(), b1_reco.Phi(), weight);
        hPhiSLBComp -> Fill(b2_parton.Phi(), b2_reco.Phi(), weight);
      }
      // fill parton level composite histograms
      hMZH -> Fill((zpartonvec+hpartonvec).M(), weight);
      hPtZH -> Fill((zpartonvec+hpartonvec).Pt(), weight);
      hEtaZH -> Fill((zpartonvec+hpartonvec).Eta(), weight);
      hPhiZH -> Fill((zpartonvec+hpartonvec).Phi(), weight);
      hDPhiZH -> Fill(deltaPhi(zpartonvec.Phi(), hpartonvec.Phi()), weight);
      if (fill_reco) {
        hMZHComp -> Fill((zpartonvec+hpartonvec).M(), sysvec.M(), weight);
        hPtZHComp -> Fill((zpartonvec+hpartonvec).Pt(), sysvec.Pt(), weight);
        hEtaZHComp -> Fill((zpartonvec+hpartonvec).Eta(), sysvec.Eta(), weight);
        hPhiZHComp -> Fill((zpartonvec+hpartonvec).Phi(), sysvec.Phi(), weight);
      }
    }
    // gather generated b jets if the particle level leptons look good
    if (!((num_elec_particle == 2 && num_mu_particle == 0) || (num_elec_particle == 0 && num_mu_particle == 2))) {
      fill_particle = false;
    }
    if (fill_particle) {
      bjets = 0;
      for(int i=0; i<(int)branchGenJet->GetEntries(); i++){
        Jet *genjet=(Jet*) branchGenJet->At(i);
        if (ghost_btag(branchGenParticle, genjet) && abs(genjet->Eta) < eta_cut && genjet->PT > jet_pt_cut) {
          bjets += 1;
          if (bjets == 1) {
            b1_particle = genjet->P4();
          } else if (bjets == 2) {
            b2_particle = genjet->P4();
          } else break;
        }
      }
      if (bjets != 2) fill_particle = false;
      hparticlevec = b1_particle + b2_particle;
      if (hparticlevec.M() > higgs_mass_cut_hi || hparticlevec.M() < higgs_mass_cut_low || hparticlevec.Pt() < higgs_pt_cut) fill_particle = false;
    }
    if (fill_particle) {
      if (num_elec_particle == 2) {
      // fill electron particle histograms
        zparticlevec = e1_particle + e2_particle;
        if (zparticlevec.M() > z_mass_cut_hi || zparticlevec.M() < z_mass_cut_low || zparticlevec.Pt() < z_pt_cut  || abs(deltaPhi(zparticlevec.Phi(), hparticlevec.Phi())) < dphi_zh_cut) {
          fill_particle = false;
        } else {
          hPtLEP -> Fill(e1_particle.Pt(), weight);
          hEtaLEP -> Fill(e1_particle.Eta(), weight);
          hPhiLEP -> Fill(e1_particle.Phi(), weight);
          hPtSLEP -> Fill(e2_particle.Pt(), weight);
          hEtaSLEP -> Fill(e2_particle.Eta(), weight);
          hPhiSLEP -> Fill(e2_particle.Phi(), weight);
          hDPhiEEP -> Fill(deltaPhi(e1_particle.Phi(), e2_particle.Phi()), weight);
          hDEtaEEP -> Fill(e1_particle.Eta()-e2_particle.Eta(), weight);
          // boost to z frame and for cos theta
          e1_particle.Boost(-zparticlevec.BoostVector());
          e2_particle.Boost(-zparticlevec.BoostVector());
          hCosThetaLEP -> Fill(e1_particle.CosTheta(), weight);
          hCosThetaSLEP -> Fill(e2_particle.CosTheta(), weight);
        // fill cos theta comparison histograms (boosting back to z frame for reco as well)
          if (fill_reco && num_elec_reco == 2) {
            elecvec1.Boost(-zvec.BoostVector());
            elecvec2.Boost(-zvec.BoostVector());
            hCosThetaLECompP -> Fill(e1_particle.CosTheta(), elecvec1.CosTheta(), weight);
            hCosThetaSLECompP -> Fill(e2_particle.CosTheta(), elecvec2.CosTheta(), weight);
            // boost back to lab frame
            elecvec1.Boost(zvec.BoostVector());
            elecvec2.Boost(zvec.BoostVector());
          }
          e1_particle.Boost(zparticlevec.BoostVector());
          e2_particle.Boost(zparticlevec.BoostVector());
          // fill other electron particle-reco comparison histograms
          if (fill_reco && num_elec_reco == 2) {
            hPtLECompP -> Fill(e1_particle.Pt(), elecvec1.Pt(), weight);
            hEtaLECompP -> Fill(e1_particle.Eta(), elecvec1.Eta(), weight);
            hPtSLECompP -> Fill(e2_particle.Pt(), elecvec2.Pt(), weight);
            hEtaSLECompP -> Fill(e2_particle.Eta(), elecvec2.Eta(), weight);
            hDPhiEECompP -> Fill(deltaPhi(e1_particle.Phi(), e2_particle.Phi()), deltaPhi(elecvec1.Phi(), elecvec2.Phi()), weight);
            hDEtaEECompP -> Fill(e1_particle.Eta()-e2_particle.Eta(), elecvec1.Eta()-elecvec2.Eta(), weight);
          }
        }
      } else if (num_mu_particle == 2){
        // LEFT OFF HERE
        // fill muon particle histograms
        zparticlevec = m1_particle + m2_particle;
        if (zparticlevec.M() > z_mass_cut_hi || zparticlevec.M() < z_mass_cut_low || zparticlevec.Pt() < z_pt_cut  || abs(deltaPhi(zparticlevec.Phi(), hparticlevec.Phi())) < dphi_zh_cut) {
          fill_particle = false;
        } else {
          hPtLMP -> Fill(m1_particle.Pt(), weight);
          hEtaLMP -> Fill(m1_particle.Eta(), weight);
          hPhiLMP -> Fill(m1_particle.Phi(), weight);
          hPtSLMP -> Fill(m2_particle.Pt(), weight);
          hEtaSLMP -> Fill(m2_particle.Eta(), weight);
          hPhiSLMP -> Fill(m2_particle.Phi(), weight);
          hDPhiMMP -> Fill(deltaPhi(m1_particle.Phi(), m2_particle.Phi()), weight);
          hDEtaMMP -> Fill(m1_particle.Eta()-m2_particle.Eta(), weight);
          // boost to z frame and for cos theta
          m1_particle.Boost(-zparticlevec.BoostVector());
          m2_particle.Boost(-zparticlevec.BoostVector());
          hCosThetaLMP -> Fill(m1_particle.CosTheta(), weight);
          hCosThetaSLMP -> Fill(m2_particle.CosTheta(), weight);
          // fill cos theta comparison histograms (boosting back to z frame for reco as well)
          if (fill_reco && num_mu_reco == 2) {
            muvec1.Boost(-zvec.BoostVector());
            muvec2.Boost(-zvec.BoostVector());
            hCosThetaLMCompP -> Fill(m1_particle.CosTheta(), muvec1.CosTheta(), weight);
            hCosThetaSLMCompP -> Fill(m2_particle.CosTheta(), muvec2.CosTheta(), weight);
            // boost back to lab frame
            muvec1.Boost(zvec.BoostVector());
            muvec2.Boost(zvec.BoostVector());
          }
          m1_particle.Boost(zparticlevec.BoostVector());
          m2_particle.Boost(zparticlevec.BoostVector());
          // fill other muon particle-reco comparison histograms
          if (fill_reco && num_mu_reco == 2) {
            hPtSLMCompP -> Fill(m2_particle.Pt(), muvec2.Pt(), weight);
            hEtaSLMCompP -> Fill(m2_particle.Eta(), muvec2.Eta(), weight);
            hPtLMCompP -> Fill(m1_particle.Pt(), muvec1.Pt(), weight);
            hEtaLMCompP -> Fill(m1_particle.Eta(), muvec1.Eta(), weight);
            hDPhiMMCompP -> Fill(deltaPhi(m1_particle.Phi(), m2_particle.Phi()), deltaPhi(muvec1.Phi(), muvec2.Phi()), weight);
            hDEtaMMCompP -> Fill(m1_particle.Eta()-m2_particle.Eta(), muvec1.Eta()-muvec2.Eta(), weight);
          }
        }
      }
    }
    if (fill_particle) {
      // fill particle level b jet histograms
      hPtLBP -> Fill(b1_particle.Pt(), weight);
      hPtSLBP -> Fill(b2_particle.Pt(), weight);
      hEtaLBP -> Fill(b1_particle.Eta(), weight);
      hEtaSLBP -> Fill(b2_particle.Eta(), weight);
      hPhiLBP -> Fill(b1_particle.Phi(), weight);
      hPhiSLBP -> Fill(b2_particle.Phi(), weight);
      // boost to higgs frame and for cos theta
      b1_particle.Boost(-hparticlevec.BoostVector());
      b2_particle.Boost(-hparticlevec.BoostVector());
      hCosThetaLBP -> Fill(b1_particle.CosTheta(), weight);
      hCosThetaSLBP -> Fill(b2_particle.CosTheta(), weight);
      // fill cos theta comparison histograms (boosting back to higgs frame for reco as well)
      if (fill_reco) {
        b1_reco.Boost(-higgsvec.BoostVector());
        b2_reco.Boost(-higgsvec.BoostVector());
        hCosThetaLBCompP -> Fill(b1_particle.CosTheta(), b1_reco.CosTheta(), weight);
        hCosThetaSLBCompP -> Fill(b2_particle.CosTheta(), b2_reco.CosTheta(), weight);
        // boost back to lab frame
        b1_reco.Boost(higgsvec.BoostVector());
        b2_reco.Boost(higgsvec.BoostVector());
      }
      b1_particle.Boost(hparticlevec.BoostVector());
      b2_particle.Boost(hparticlevec.BoostVector());
      // fill other b jet comparison histograms
      if (fill_reco) {
        hPtLBCompP -> Fill(b1_particle.Pt(), b1_reco.Pt(), weight);
        hPtSLBCompP -> Fill(b2_particle.Pt(), b2_reco.Pt(), weight);
        hEtaLBCompP -> Fill(b1_particle.Eta(), b1_reco.Eta(), weight);
        hEtaSLBCompP -> Fill(b2_particle.Eta(), b2_reco.Eta(), weight);
      }
      // fill mass histograms
      hMbbP -> Fill(hparticlevec.M(), weight);
      hMllP -> Fill(zparticlevec.M(), weight);
      hMbbllP -> Fill((hparticlevec+zparticlevec).M(), weight);
      if (fill_reco) {
        hMZHCompP -> Fill((zparticlevec+hparticlevec).M(), sysvec.M(), weight);
        hMbbComp -> Fill(hparticlevec.M(), higgsvec.M(), weight);
        hMllComp -> Fill(zparticlevec.M(), zvec.M(), weight);
      }
      // fill higgs and z histograms
      hPtBBP -> Fill(hparticlevec.Pt(), weight);
      hPtLLP -> Fill(zparticlevec.Pt(), weight);
      hEtaBBP -> Fill(hparticlevec.Eta(), weight);
      hEtaLLP -> Fill(zparticlevec.Eta(), weight);
      hPhiBBP -> Fill(hparticlevec.Phi(), weight);
      hPhiLLP -> Fill(zparticlevec.Phi(), weight);
      // fill higgs and z particle-reco comparison histograms
      if (fill_reco) {
        hPtHCompP -> Fill(hparticlevec.Pt(), higgsvec.Pt(), weight);
        hEtaHCompP -> Fill(hparticlevec.Eta(), higgsvec.Eta(), weight);
        hPtZCompP -> Fill(zparticlevec.Pt(), zvec.Pt(), weight);
        hEtaZCompP -> Fill(zparticlevec.Eta(), zvec.Eta(), weight);
      }
      // fill total event histograms
      hPtBBLLP -> Fill((hparticlevec+zparticlevec).Pt(), weight);
      hEtaBBLLP -> Fill((hparticlevec+zparticlevec).Eta(), weight);
      hPhiBBLLP -> Fill((hparticlevec+zparticlevec).Phi(), weight);
      hDPhiZHP -> Fill(deltaPhi(zparticlevec.Phi(), hparticlevec.Phi()), weight);
      if (fill_reco) {
        hPtZHCompP -> Fill((zparticlevec+hparticlevec).Pt(), sysvec.Pt(), weight);
        hEtaZHCompP -> Fill((zparticlevec+hparticlevec).Eta(), sysvec.Eta(), weight);
      }
    }
    nPassed+=weight;
    nPassedRaw++;
  }
  // Fit higgs mass distribution (can also do Z but W is less insightful due to MET)
  /*TF1 *fitMHR = new TF1("higgs_dscb_reco", DSCBf, 20, 250, 7);
  fitMHR -> SetParameters(1, 2, 1, 2, 125, 10, 1);
  fitMHR -> SetParNames("alpha_low", "n_low", "alpha_high", "n_high", "mean", "sigma", "N");
  TFitResultPtr higgs_results_reco = hMbbR -> Fit("higgs_dscb_reco", "S");
  std::cout<<"chi2 for higgs reco fit "<< higgs_results_reco->Chi2()<<std::endl;
  TF1 *fitMHP = new TF1("higgs_dscb_particle", DSCBf, 20, 250, 7);
  fitMHP -> SetParameters(1, 2, 1, 2, 125, 10, 1);
  fitMHP -> SetParNames("alpha_low", "n_low", "alpha_high", "n_high", "mean", "sigma", "N");
  TFitResultPtr higgs_results_particle = hMbbP -> Fit("higgs_dscb_particle", "S");
  std::cout<<"chi2 for higgs particle fit "<< higgs_results_particle->Chi2()<<std::endl;

  TF1 *fitMZR = new TF1("z_dscb_reco", DSCBf, 20, 250, 7);
  fitMZR -> SetParameters(1, 2, 1, 2, 90, 5, 1);
  fitMZR -> SetParNames("alpha_low", "n_low", "alpha_high", "n_high", "mean", "sigma", "N");
  TFitResultPtr z_results_reco = hMllR -> Fit("z_dscb_reco", "S");
  std::cout<<"chi2 for Z reco fit "<< z_results_reco->Chi2()<<std::endl;
  TF1 *fitMZP = new TF1("z_dscb_particle", DSCBf, 20, 250, 7);
  fitMZP -> SetParameters(1, 2, 1, 2, 90, 5, 1);
  fitMZP -> SetParNames("alpha_low", "n_low", "alpha_high", "n_high", "mean", "sigma", "N");
  TFitResultPtr z_results_particle = hMllP -> Fit("z_dscb_particle", "S");
  std::cout<<"chi2 for Z particle fit "<< z_results_particle->Chi2()<<std::endl;*/

  // Write histograms
  TFile *hists= new TFile(outputFile,"recreate");
  hists->cd();
  //mass
  hWeight->Write();
  hMZH->Write();
  hMbbR->Write();
  hMllR->Write();
  hMbbllR->Write();
  hMbbP->Write();
  hMllP->Write();
  hMbbllP->Write();
  hMbbComp->Write();
  hMllComp->Write();
  hMZHComp->Write();
  hMZHCompP->Write();
  // pt
  hPtH->Write();
  hPtZ->Write();
  hPtLE->Write();
  hPtLM->Write();
  hPtSLE->Write();
  hPtSLM->Write();
  hPtLB->Write();
  hPtSLB->Write();
  hPtZH->Write();
  hPtBBR->Write();
  hPtLLR->Write();
  hPtLER->Write();
  hPtLMR->Write();
  hPtSLER->Write();
  hPtSLMR->Write();
  hPtLBR->Write();
  hPtSLBR->Write();
  hPtBBLLR->Write();
  hPtBBP->Write();
  hPtLLP->Write();
  hPtLEP->Write();
  hPtLMP->Write();
  hPtSLEP->Write();
  hPtSLMP->Write();
  hPtLBP->Write();
  hPtSLBP->Write();
  hPtBBLLP->Write();
  hPtHComp->Write();
  hPtZComp->Write();
  hPtLEComp->Write();
  hPtLMComp->Write();
  hPtSLEComp->Write();
  hPtSLMComp->Write();
  hPtLBComp->Write();
  hPtSLBComp->Write();
  hPtZHComp->Write();
  hPtHCompP->Write();
  hPtZCompP->Write();
  hPtLECompP->Write();
  hPtLMCompP->Write();
  hPtSLECompP->Write();
  hPtSLMCompP->Write();
  hPtLBCompP->Write();
  hPtSLBCompP->Write();
  hPtZHCompP->Write();
  //eta
  hEtaH->Write();
  hEtaZ->Write();
  hEtaLE->Write();
  hEtaLM->Write();
  hEtaSLE->Write();
  hEtaSLM->Write();
  hEtaLB->Write();
  hEtaSLB->Write();
  hEtaZH->Write();
  hEtaBBR->Write();
  hEtaLLR->Write();
  hEtaLER->Write();
  hEtaLMR->Write();
  hEtaSLER->Write();
  hEtaSLMR->Write();
  hEtaLBR->Write();
  hEtaSLBR->Write();
  hEtaBBLLR->Write();
  hEtaBBP->Write();
  hEtaLLP->Write();
  hEtaLEP->Write();
  hEtaLMP->Write();
  hEtaSLEP->Write();
  hEtaSLMP->Write();
  hEtaLBP->Write();
  hEtaSLBP->Write();
  hEtaBBLLP->Write();
  hEtaHComp->Write();
  hEtaZComp->Write();
  hEtaLEComp->Write();
  hEtaLMComp->Write();
  hEtaSLEComp->Write();
  hEtaSLMComp->Write();
  hEtaLBComp->Write();
  hEtaSLBComp->Write();
  hEtaZHComp->Write();
  hEtaHCompP->Write();
  hEtaZCompP->Write();
  hEtaLECompP->Write();
  hEtaLMCompP->Write();
  hEtaSLECompP->Write();
  hEtaSLMCompP->Write();
  hEtaLBCompP->Write();
  hEtaSLBCompP->Write();
  hEtaZHCompP->Write();
  hDEtaEE->Write();
  hDEtaMM->Write();
  hDEtaEER->Write();
  hDEtaMMR->Write();
  hDEtaEEP->Write();
  hDEtaMMP->Write();
  hDEtaEEComp->Write();
  hDEtaEECompP->Write();
  hDEtaMMComp->Write();
  hDEtaMMCompP->Write();
   //cos theta
  hCosThetaLE->Write();
  hCosThetaLM->Write();
  hCosThetaSLE->Write();
  hCosThetaSLM->Write();
  hCosThetaLB->Write();
  hCosThetaSLB->Write();
  hCosThetaLER->Write();
  hCosThetaLMR->Write();
  hCosThetaSLER->Write();
  hCosThetaSLMR->Write();
  hCosThetaLBR->Write();
  hCosThetaSLBR->Write();
  hCosThetaLEP->Write();
  hCosThetaLMP->Write();
  hCosThetaSLEP->Write();
  hCosThetaSLMP->Write();
  hCosThetaLBP->Write();
  hCosThetaSLBP->Write();
  hCosThetaLEComp->Write();
  hCosThetaLMComp->Write();
  hCosThetaSLEComp->Write();
  hCosThetaSLMComp->Write();
  hCosThetaLBComp->Write();
  hCosThetaSLBComp->Write();
  hCosThetaLECompP->Write();
  hCosThetaLMCompP->Write();
  hCosThetaSLECompP->Write();
  hCosThetaSLMCompP->Write();
  hCosThetaLBCompP->Write();
  hCosThetaSLBCompP->Write();
  // phi
  hPhiH->Write();
  hPhiZ->Write();
  hPhiLE->Write();
  hPhiLM->Write();
  hPhiSLE->Write();
  hPhiSLM->Write();
  hDPhiEE->Write();
  hDPhiMM->Write();
  hDPhiZH->Write();
  hPhiLB->Write();
  hPhiSLB->Write();
  hPhiZH->Write();
  hPhiBBR->Write();
  hPhiLLR->Write();
  hPhiLER->Write();
  hPhiLMR->Write();
  hPhiSLER->Write();
  hPhiSLMR->Write();
  hDPhiEER->Write();
  hDPhiMMR->Write();
  hDPhiZHR->Write();
  hPhiLBR->Write();
  hPhiSLBR->Write();
  hPhiBBLLR->Write();
  hPhiBBP->Write();
  hPhiLLP->Write();
  hPhiLEP->Write();
  hPhiLMP->Write();
  hPhiSLEP->Write();
  hPhiSLMP->Write();
  hDPhiEEP->Write();
  hDPhiMMP->Write();
  hDPhiZHP->Write();
  hPhiLBP->Write();
  hPhiSLBP->Write();
  hPhiBBLLP->Write();
  hPhiHComp->Write();
  hPhiZComp->Write();
  hPhiLEComp->Write();
  hPhiLMComp->Write();
  hPhiSLEComp->Write();
  hPhiSLMComp->Write();
  hDPhiEEComp->Write();
  hDPhiEECompP->Write();
  hDPhiMMComp->Write();
  hDPhiMMCompP->Write();
  hPhiLBComp->Write();
  hPhiSLBComp->Write();
  hPhiZHComp->Write();
  hists->Close();
  // clear all the histograms
  //mass
  hWeight->Clear();
  hMZH->Clear();
  hMbbR->Clear();
  hMllR->Clear();
  hMbbllR->Clear();
  hMbbP->Clear();
  hMllP->Clear();
  hMbbllP->Clear();
  hMbbComp->Clear();
  hMllComp->Clear();
  hMZHComp->Clear();
  hMZHCompP->Clear();
  // pt
  hPtH->Clear();
  hPtZ->Clear();
  hPtLE->Clear();
  hPtLM->Clear();
  hPtSLE->Clear();
  hPtSLM->Clear();
  hPtLB->Clear();
  hPtSLB->Clear();
  hPtZH->Clear();
  hPtBBR->Clear();
  hPtLLR->Clear();
  hPtLER->Clear();
  hPtLMR->Clear();
  hPtSLER->Clear();
  hPtSLMR->Clear();
  hPtLBR->Clear();
  hPtSLBR->Clear();
  hPtBBLLR->Clear();
  hPtBBP->Clear();
  hPtLLP->Clear();
  hPtLEP->Clear();
  hPtLMP->Clear();
  hPtSLEP->Clear();
  hPtSLMP->Clear();
  hPtLBP->Clear();
  hPtSLBP->Clear();
  hPtBBLLP->Clear();
  hPtHComp->Clear();
  hPtZComp->Clear();
  hPtLEComp->Clear();
  hPtLMComp->Clear();
  hPtSLEComp->Clear();
  hPtSLMComp->Clear();
  hPtLBComp->Clear();
  hPtSLBComp->Clear();
  hPtZHComp->Clear();
  hPtHCompP->Clear();
  hPtZCompP->Clear();
  hPtLECompP->Clear();
  hPtLMCompP->Clear();
  hPtSLECompP->Clear();
  hPtSLMCompP->Clear();
  hPtLBCompP->Clear();
  hPtSLBCompP->Clear();
  hPtZHCompP->Clear();
  //eta
  hEtaH->Clear();
  hEtaZ->Clear();
  hEtaLE->Clear();
  hEtaLM->Clear();
  hEtaSLE->Clear();
  hEtaSLM->Clear();
  hEtaLB->Clear();
  hEtaSLB->Clear();
  hEtaZH->Clear();
  hEtaBBR->Clear();
  hEtaLLR->Clear();
  hEtaLER->Clear();
  hEtaLMR->Clear();
  hEtaSLER->Clear();
  hEtaSLMR->Clear();
  hEtaLBR->Clear();
  hEtaSLBR->Clear();
  hEtaBBLLR->Clear();
  hEtaBBP->Clear();
  hEtaLLP->Clear();
  hEtaLEP->Clear();
  hEtaLMP->Clear();
  hEtaSLEP->Clear();
  hEtaSLMP->Clear();
  hEtaLBP->Clear();
  hEtaSLBP->Clear();
  hEtaBBLLP->Clear();
  hEtaHComp->Clear();
  hEtaZComp->Clear();
  hEtaLEComp->Clear();
  hEtaLMComp->Clear();
  hEtaSLEComp->Clear();
  hEtaSLMComp->Clear();
  hEtaLBComp->Clear();
  hEtaSLBComp->Clear();
  hEtaZHComp->Clear();
  hEtaHCompP->Clear();
  hEtaZCompP->Clear();
  hEtaLECompP->Clear();
  hEtaLMCompP->Clear();
  hEtaSLECompP->Clear();
  hEtaSLMCompP->Clear();
  hEtaLBCompP->Clear();
  hEtaSLBCompP->Clear();
  hEtaZHCompP->Clear();
  hDEtaEE->Clear();
  hDEtaMM->Clear();
  hDEtaEER->Clear();
  hDEtaMMR->Clear();
  hDEtaEEP->Clear();
  hDEtaMMP->Clear();
  hDEtaEEComp->Clear();
  hDEtaEECompP->Clear();
  hDEtaMMComp->Clear();
  hDEtaMMCompP->Clear();
   //cos theta
  hCosThetaLE->Clear();
  hCosThetaLM->Clear();
  hCosThetaSLE->Clear();
  hCosThetaSLM->Clear();
  hCosThetaLB->Clear();
  hCosThetaSLB->Clear();
  hCosThetaLER->Clear();
  hCosThetaLMR->Clear();
  hCosThetaSLER->Clear();
  hCosThetaSLMR->Clear();
  hCosThetaLBR->Clear();
  hCosThetaSLBR->Clear();
  hCosThetaLEP->Clear();
  hCosThetaLMP->Clear();
  hCosThetaSLEP->Clear();
  hCosThetaSLMP->Clear();
  hCosThetaLBP->Clear();
  hCosThetaSLBP->Clear();
  hCosThetaLEComp->Clear();
  hCosThetaLMComp->Clear();
  hCosThetaSLEComp->Clear();
  hCosThetaSLMComp->Clear();
  hCosThetaLBComp->Clear();
  hCosThetaSLBComp->Clear();
  hCosThetaLECompP->Clear();
  hCosThetaLMCompP->Clear();
  hCosThetaSLECompP->Clear();
  hCosThetaSLMCompP->Clear();
  hCosThetaLBCompP->Clear();
  hCosThetaSLBCompP->Clear();
  // phi
  hPhiH->Clear();
  hPhiZ->Clear();
  hPhiLE->Clear();
  hPhiLM->Clear();
  hPhiSLE->Clear();
  hPhiSLM->Clear();
  hDPhiEE->Clear();
  hDPhiMM->Clear();
  hDPhiZH->Clear();
  hPhiLB->Clear();
  hPhiSLB->Clear();
  hPhiZH->Clear();
  hPhiBBR->Clear();
  hPhiLLR->Clear();
  hPhiLER->Clear();
  hPhiLMR->Clear();
  hPhiSLER->Clear();
  hPhiSLMR->Clear();
  hDPhiEER->Clear();
  hDPhiMMR->Clear();
  hDPhiZHR->Clear();
  hPhiLBR->Clear();
  hPhiSLBR->Clear();
  hPhiBBLLR->Clear();
  hPhiBBP->Clear();
  hPhiLLP->Clear();
  hPhiLEP->Clear();
  hPhiLMP->Clear();
  hPhiSLEP->Clear();
  hPhiSLMP->Clear();
  hDPhiEEP->Clear();
  hDPhiMMP->Clear();
  hDPhiZHP->Clear();
  hPhiLBP->Clear();
  hPhiSLBP->Clear();
  hPhiBBLLP->Clear();
  hPhiHComp->Clear();
  hPhiZComp->Clear();
  hPhiLEComp->Clear();
  hPhiLMComp->Clear();
  hPhiSLEComp->Clear();
  hPhiSLMComp->Clear();
  hDPhiEEComp->Clear();
  hDPhiEECompP->Clear();
  hDPhiMMComp->Clear();
  hDPhiMMCompP->Clear();
  hPhiLBComp->Clear();
  hPhiSLBComp->Clear();
  hPhiZHComp->Clear();
}

int main(int argc, char* argv[]) {
  const char *inputFileName = argv[1];
  const char *outputFileName = argv[2];
  const char *process_name = argv[3];
  zhbb_analyze(inputFileName, outputFileName, process_name);
  return 1;
}
