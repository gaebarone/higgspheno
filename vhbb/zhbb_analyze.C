// root -l -b
// .x examples/whbb_analyze.C("../wplus_hbb/Events/run_01/tag_1_delphes_events.root", "../wminus_hbb/Events/run_01/tag_1_delphes_events.root")
#ifdef __CLING__
R__LOAD_LIBRARY(libDelphes)
#include "classes/DelphesClasses.h"
#include "external/ExRootAnalysis/ExRootTreeReader.h"
#include "DSCBf.h"
#include "ghost_tagging.h"
#endif

//------------------------------------------------------------------------------
// save plots as different file types
void PrintCanvas(TCanvas *c=nullptr,string name="default"){
  std::vector <string> types={"jpg"}; 
  for(std::vector<string>::iterator it=types.begin(); it!=types.end(); it++) {
    c->Print(Form("zhbb_plots/%s.%s",name.c_str(),(*it).c_str()),(*it).c_str());
  }
}

// draw 1d histograms
void draw_hist(TH1 *histo, const char *name, const char *title, const char *axistitle) {
  TCanvas *c = new TCanvas(name, title, 1500,1200);
  c->cd();
  histo->GetXaxis()->SetTitle(axistitle);
  histo->Draw();
  PrintCanvas(c, name);
}

// draw 2d histograms
void draw_hist2(TH2 *histo, const char *name, const char *title, const char *xaxistitle, const char *yaxistitle) {
  TCanvas *c = new TCanvas(name, title, 1500, 1200);
  histo->GetXaxis()->SetTitle(xaxistitle);
  histo->GetYaxis()->SetTitle(yaxistitle);
  histo->Draw("COLZ");
  PrintCanvas(c, name); 
}

// make a ton of plots for zhbb events (z -> l l)
void zhbb_analyze(const char *inputFile) {
  gSystem->Load("libDelphes");
  // Create chain of root trees
  TChain chain("Delphes");
  chain.Add(inputFile);
  // Create object of class ExRootTreeReader
  ExRootTreeReader *treeReader = new ExRootTreeReader(&chain);
  Long64_t numberOfEntries = treeReader->GetEntries();
  // Get pointers to branches used in this analysis
  TClonesArray *branchJet = treeReader->UseBranch("Jet");
  TClonesArray *branchElectron = treeReader->UseBranch("Electron");
  TClonesArray *branchMuon = treeReader->UseBranch("Muon");
  TClonesArray *branchEvent = treeReader->UseBranch("Event");
  TClonesArray *branchGenParticle = treeReader->UseBranch("Particle");
  TClonesArray *branchGenJet = treeReader->UseBranch("GenJet");
  // SET CUTS
  const double e_pt_cut_lead = 27;
  const double mu_pt_cut_lead = 20;
  const double e_pt_cut_sub = 15;
  const double mu_pt_cut_sub = 12;
  const double eta_cut = 2.5;
  // Book histograms
  // mass
  TH1 *hMZH = new TH1F("mass_ZH", "m_{ZH}", 100, 120.0, 670.0);
  TH1 *hMbbR = new TH1F("mass_bb_reco", "Reco m_{bb}", 100, 0.0, 250.0);
  TH1 *hMllR = new TH1F("mass_ll_reco", "Reco m_{ll}", 100, 40.0, 140.0);
  TH1 *hMbbllR = new TH1F("mass_bbll_reco", "Reco m_{bbll}", 100, 120.0, 670.0);
  TH1 *hMbbP = new TH1F("mass_bb_particle", "m_{bb}", 100, 0.0, 250.0);
  TH1 *hMllP = new TH1F("mass_ll_particle", "m_{ll}", 100, 40.0, 140.0);
  TH1 *hMbbllP = new TH1F("mass_bbll_particle", "m_{bbll}", 100, 120.0, 670.0);
  TH2 *hMbbComp = new TH2F("mass_bb_Comp", "m_{bb}", 50, 0.0, 250.0, 50, 0.0, 250.0);
  TH2 *hMllComp = new TH2F("mass_ll_Comp", "m_{ll}", 50, 40.0, 140.0, 50, 40.0, 140.0);
  TH2 *hMZHComp = new TH2F("mass_ZH_Comp", "m_{ZH}", 50, 120.0, 670.0, 50, 120.0, 670.0);
  TH2 *hMZHCompP = new TH2F("mass_ZH_Comp_particle", "m_{ZH}", 50, 120.0, 670.0, 50, 120.0, 670.0);
  // pt parton level
  TH1 *hPtH = new TH1F("pt_H", "p_{T}^{H}", 100, 0.0, 400.0);
  TH1 *hPtZ = new TH1F("pt_Z", "p_{T}^{Z}", 100, 0.0, 400.0);
  TH1 *hPtLE = new TH1F("pt_LE",  "Lead p_{T}^{e}", 100, 20.0, 200.0);
  TH1 *hPtLM = new TH1F("pt_LM",  "Lead p_{T}^{#mu}", 100, 20.0, 200.0);
  TH1 *hPtSLE = new TH1F("pt_SLE",  "Sublead p_{T}^{e}", 100, 20.0, 200.0);
  TH1 *hPtSLM = new TH1F("pt_SLM",  "Sublead p_{T}^{#mu}", 100, 20.0, 200.0);
  TH1 *hPtLB = new TH1F("pt_LB", "Lead b Parton p_{T}", 100, 20.0, 200.0);
  TH1 *hPtSLB = new TH1F("pt_SLB", "Sublead b Parton p_{T}", 100, 20.0, 200.0);
  TH1 *hPtZH = new TH1F("pt_ZH", "p_{T}^{ZH}", 100, 0.0, 400.0);
  // pt reco level
  TH1 *hPtBBR = new TH1F("pt_bb_reco", "Reco p_{T}^{bb}", 100, 0.0, 400.0);
  TH1 *hPtLLR = new TH1F("pt_ll_reco", "Reco p_{T}^{ll}", 100, 0.0, 400.0);
  TH1 *hPtLER = new TH1F("pt_LE_reco",  "Lead Reco p_{T}^{e}", 100, 20.0, 200.0);
  TH1 *hPtLMR = new TH1F("pt_LM_reco",  "Lead Reco p_{T}^{#mu}", 100, 20.0, 200.0);
  TH1 *hPtSLER = new TH1F("pt_SLE_reco",  "Sublead Reco p_{T}^{e}", 100, 20.0, 200.0);
  TH1 *hPtSLMR = new TH1F("pt_SLM_reco",  "Sublead Reco p_{T}^{#mu}", 100, 20.0, 200.0);
  TH1 *hPtLBR = new TH1F("pt_LB_reco", "Lead Reco b Jet p_{T}", 100, 20.0, 200.0);
  TH1 *hPtSLBR = new TH1F("pt_SLB_reco", "Sublead Reco b Jet p_{T}", 100, 20.0, 200.0);
  TH1 *hPtBBLLR = new TH1F("pt_bbll_reco", "Reco p_{T}^{bbll}", 100, 0.0, 400.0);
  // pt particle level
  TH1 *hPtBBP = new TH1F("pt_bb_particle", "p_{T}^{bb}", 100, 0.0, 400.0);
  TH1 *hPtLLP = new TH1F("pt_ll_particle", "p_{T}^{ll}", 100, 0.0, 400.0);
  TH1 *hPtLEP = new TH1F("pt_LE_particle",  "Lead p_{T}^{e}", 100, 20.0, 200.0);
  TH1 *hPtLMP = new TH1F("pt_LM_particle",  "Lead p_{T}^{#mu}", 100, 20.0, 200.0);
  TH1 *hPtSLEP = new TH1F("pt_SLE_particle",  "Sublead p_{T}^{e}", 100, 20.0, 200.0);
  TH1 *hPtSLMP = new TH1F("pt_SLM_particle",  "Sublead p_{T}^{#mu}", 100, 20.0, 200.0);
  TH1 *hPtLBP = new TH1F("pt_LB_particle", "Lead b Jet p_{T}", 100, 20.0, 200.0);
  TH1 *hPtSLBP = new TH1F("pt_SLB_particle", "Sublead b Jet p_{T}", 100, 20.0, 200.0);
  TH1 *hPtBBLLP = new TH1F("pt_bbll_particle", "p_{T}^{bbll}", 100, 0.0, 400.0);
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
  TH1 *hEtaH = new TH1F("eta_H", "#eta_{H}", 100, -3.0, 3.0);
  TH1 *hEtaZ = new TH1F("eta_Z", "#eta_{Z}", 100, -3.0, 3.0);
  TH1 *hEtaLE = new TH1F("eta_LE",  "Lead #eta_{e}", 100, -3.0, 3.0);
  TH1 *hEtaLM = new TH1F("eta_LM",  "Lead #eta_{#mu}", 100, -3.0, 3.0);
  TH1 *hEtaSLE = new TH1F("eta_SLE",  "Sublead #eta_{e}", 100, -3.0, 3.0);
  TH1 *hEtaSLM = new TH1F("eta_SLM",  "Sublead #eta_{#mu}", 100, -3.0, 3.0);
  TH1 *hEtaLB = new TH1F("eta_LB", "Lead #eta_{b}", 100, -3.0, 3.0);
  TH1 *hEtaSLB = new TH1F("eta_SLB", "Sublead #eta_{b}", 100, -3.0, 3.0);
  TH1 *hEtaZH = new TH1F("eta_ZH", "#eta_{ZH}", 100, -3.0, 3.0);
  // eta reco level
  TH1 *hEtaBBR = new TH1F("eta_bb_reco", "Reco #eta_{bb}", 100, -3.0, 3.0);
  TH1 *hEtaLLR = new TH1F("eta_ll_reco", "Reco #eta_{ll}", 100, -3.0, 3.0);
  TH1 *hEtaLER = new TH1F("eta_LE_reco",  "Lead Reco #eta_{e}", 100, -3.0, 3.0);
  TH1 *hEtaLMR = new TH1F("eta_LM_reco",  "Lead Reco #eta_{#mu}", 100, -3.0, 3.0);
  TH1 *hEtaSLER = new TH1F("eta_SLE_reco",  "Sublead Reco #eta_{e}", 100, -3.0, 3.0);
  TH1 *hEtaSLMR = new TH1F("eta_SLM_reco",  "Sublead Reco #eta_{#mu}", 100, -3.0, 3.0);
  TH1 *hEtaLBR = new TH1F("eta_LB_reco", "Lead Reco b Jet #eta", 100, -3.0, 3.0);
  TH1 *hEtaSLBR = new TH1F("eta_SLB_reco", "Sublead Reco b Jet #eta", 100, -3.0, 3.0);
  TH1 *hEtaBBLLR = new TH1F("eta_bbll_reco", "Reco #eta_{bbll}", 100, -3.0, 3.0);
  // eta particle level
  TH1 *hEtaBBP = new TH1F("eta_bb_particle", "#eta_{bb}", 100, -3.0, 3.0);
  TH1 *hEtaLLP = new TH1F("eta_ll_particle", "#eta_{ll}", 100, -3.0, 3.0);
  TH1 *hEtaLEP = new TH1F("eta_LE_particle",  "Lead #eta_{e}", 100, -3.0, 3.0);
  TH1 *hEtaLMP = new TH1F("eta_LM_particle",  "Lead #eta_{#mu}", 100, -3.0, 3.0);
  TH1 *hEtaSLEP = new TH1F("eta_SLE_particle",  "Sublead #eta_{e}", 100, -3.0, 3.0);
  TH1 *hEtaSLMP = new TH1F("eta_SLM_particle",  "Sublead #eta_{#mu}", 100, -3.0, 3.0);
  TH1 *hEtaLBP = new TH1F("eta_LB_particle", "Lead b Jet #eta", 100, -3.0, 3.0);
  TH1 *hEtaSLBP = new TH1F("eta_SLB_particle", "Sublead b Jet #eta", 100, -3.0, 3.0);
  TH1 *hEtaBBLLP = new TH1F("eta_bbll_particle", "#eta_{bbll}", 100, -3.0, 3.0);
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
  // phi parton level
  TH1 *hPhiH = new TH1F("phi_H", "#phi_{H}", 100, -3.15, 3.15);
  TH1 *hPhiZ = new TH1F("phi_Z", "#phi_{Z}", 100, -3.15, 3.15);
  TH1 *hPhiLE = new TH1F("phi_LE",  "Lead #phi_{e}", 100, -3.15, 3.15);
  TH1 *hPhiLM = new TH1F("phi_LM",  "Lead #phi_{#mu}", 100, -3.15, 3.15);
  TH1 *hPhiSLE = new TH1F("phi_SLE",  "Sublead #phi_{e}", 100, -3.15, 3.15);
  TH1 *hPhiSLM = new TH1F("phi_SLM",  "Sublead #phi_{#mu}", 100, -3.15, 3.15);
  TH1 *hPhiLB = new TH1F("phi_LB", "Lead #phi_{b}", 100, -3.15, 3.15);
  TH1 *hPhiSLB = new TH1F("phi_SLB", "Sublead #phi_{b}", 100, -3.15, 3.15);
  TH1 *hPhiZH = new TH1F("phi_ZH", "#phi_{ZH}", 100, -3.15, 3.15);
  // phi reco level
  TH1 *hPhiBBR = new TH1F("phi_bb_reco", "Reco #phi_{bb}", 100, -3.15, 3.15);
  TH1 *hPhiLLR = new TH1F("phi_ll_reco", "Reco #phi_{ll}", 100, -3.15, 3.15);
  TH1 *hPhiLER = new TH1F("phi_LE_reco",  "Lead Reco #phi_{e}", 100, -3.15, 3.15);
  TH1 *hPhiLMR = new TH1F("phi_LM_reco",  "Lead Reco #phi_{#mu}", 100, -3.15, 3.15);
  TH1 *hPhiSLER = new TH1F("phi_SLE_reco",  "Sublead Reco #phi_{e}", 100, -3.15, 3.15);
  TH1 *hPhiSLMR = new TH1F("phi_SLM_reco",  "Sublead Reco #phi_{#mu}", 100, -3.15, 3.15);
  TH1 *hPhiLBR = new TH1F("phi_LB_reco", "Lead Reco b Jet #phi", 100, -3.15, 3.15);
  TH1 *hPhiSLBR = new TH1F("phi_SLB_reco", "Sublead Reco b Jet #phi", 100, -3.15, 3.15);
  TH1 *hPhiBBLLR = new TH1F("phi_bbll_reco", "Reco #phi_{bbll}", 100, -3.15, 3.15);
  // phi particle level
  TH1 *hPhiBBP = new TH1F("phi_bb_particle", "#phi_{bb}", 100, -3.15, 3.15);
  TH1 *hPhiLLP = new TH1F("phi_ll_particle", "#phi_{ll}", 100, -3.15, 3.15);
  TH1 *hPhiLEP = new TH1F("phi_LE_particle",  "Lead #phi_{e}", 100, -3.15, 3.15);
  TH1 *hPhiLMP = new TH1F("phi_LM_particle",  "Lead #phi_{#mu}", 100, -3.15, 3.15);
  TH1 *hPhiSLEP = new TH1F("phi_SLE_particle",  "Sublead #phi_{e}", 100, -3.15, 3.15);
  TH1 *hPhiSLMP = new TH1F("phi_SLM_particle",  "Sublead #phi_{#mu}", 100, -3.15, 3.15);
  TH1 *hPhiLBP = new TH1F("phi_LB_particle", "Lead b Jet #phi", 100, -3.15, 3.15);
  TH1 *hPhiSLBP = new TH1F("phi_SLB_particle", "Sublead b Jet #phi", 100, -3.15, 3.15);
  TH1 *hPhiBBLLP = new TH1F("phi_bbll_particle", "#phi_{bbll}", 100, -3.15, 3.15);
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

  // initialize variables needed for filling histograms
  double  nPassed=0;
  int  nPassedRaw=0;
  double Lumi=3e3;
  int bjets = 0;
  int num_elec_reco;
  int num_mu_reco;
  int num_elec_particle;
  int num_mu_particle;
  bool elec_ev_parton;
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
  Electron *elec;
  Electron *elec1;
  Electron *elec2;
  Muon *muon;
  Muon *muon1;
  Muon *muon2;
  GenParticle *daughter1;
  GenParticle *daughter2;
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

  // Loop over events
  for(Int_t entry = 0; entry < numberOfEntries; ++entry){
    // Load selected branches with data from specified event
    treeReader->ReadEntry(entry);
    HepMCEvent *event = (HepMCEvent*) branchEvent -> At(0);
    Float_t weight = event->Weight/numberOfEntries*Lumi;
    bool fill_reco = true;
    bool fill_parton = true;
    bool fill_particle = true;
    // Loop over jets in event and save the b jets in 2b events as lorentz vectors to reconstruct the higgs
    bjets = 0;
    for(int i=0; i<(int)branchJet->GetEntries(); i++){
      Jet *jet=(Jet*) branchJet->At(i);
      if (jet -> BTag == 1 && abs(jet -> Eta) < eta_cut) {
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
      hPtLER -> Fill(elec1->PT, weight);
      hEtaLER -> Fill(elec1->Eta, weight);
      hPhiLER -> Fill(elec1->Phi, weight);
      hPtSLER -> Fill(elec2->PT, weight);
      hEtaSLER -> Fill(elec2->Eta, weight);
      hPhiSLER -> Fill(elec2->Phi, weight);
    // fill muon data
    } else if (fill_reco && num_elec_reco == 0 && num_mu_reco == 2) {
      muvec1 = muon1->P4();
      muvec2 = muon2->P4();
      zvec = muvec1 + muvec2;
      hPtLMR -> Fill(muon1->PT, weight);
      hEtaLMR -> Fill(muon1->Eta, weight);
      hPhiLMR -> Fill(muon1->Phi, weight);
      hPtSLMR -> Fill(muon2->PT, weight);
      hEtaSLMR -> Fill(muon2->Eta, weight);
      hPhiSLMR -> Fill(muon2->Phi, weight);
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
        } else fill_parton = false;
      // Z parton
      } else if ((particle->PID == 23) && (d1_pid != 23) && (d2_pid != 23)) {
        zpartonvec = particle->P4();
        // check for electron daughters
        if (abs(d1_pid) == 11) {
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
        } else if (abs(d1_pid) == 13) {
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
        }
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
    if (fill_parton) {
      // fill electron parton histograms
      if (elec_ev_parton) {
        hPtLE -> Fill(e1_parton.Pt(), weight);
        hEtaLE -> Fill(e1_parton.Eta(), weight);
        hPhiLE -> Fill(e1_parton.Phi(), weight);
        hPtSLE -> Fill(e2_parton.Pt(), weight);
        hEtaSLE -> Fill(e2_parton.Eta(), weight);
        hPhiSLE -> Fill(e2_parton.Phi(), weight);
        if (fill_reco && num_elec_reco == 2) {
          // fill electron comparison histograms
          hPtLEComp -> Fill(e1_parton.Pt(), elecvec1.Pt(), weight);
          hEtaLEComp -> Fill(e1_parton.Eta(), elecvec1.Eta(), weight);
          hPhiLEComp -> Fill(e1_parton.Phi(), elecvec1.Phi(), weight);
          hPtSLEComp -> Fill(e2_parton.Pt(), elecvec2.Pt(), weight);
          hEtaSLEComp -> Fill(e2_parton.Eta(), elecvec2.Eta(), weight);
          hPhiSLEComp -> Fill(e2_parton.Phi(), elecvec2.Phi(), weight);
        }
      } else {
        // fill muon parton histograms
        hPtLM -> Fill(m1_parton.Pt(), weight);
        hEtaLM -> Fill(m1_parton.Eta(), weight);
        hPhiLM -> Fill(m1_parton.Phi(), weight);
        hPtSLM -> Fill(m2_parton.Pt(), weight);
        hEtaSLM -> Fill(m2_parton.Eta(), weight);
        hPhiSLM -> Fill(m2_parton.Phi(), weight);
        if (fill_reco && num_mu_reco == 2) {
          // fill muon comparison histograms
          hPtSLMComp -> Fill(m2_parton.Pt(), muvec2.Pt(), weight);
          hEtaSLMComp -> Fill(m2_parton.Eta(), muvec2.Eta(), weight);
          hPhiSLMComp -> Fill(m2_parton.Phi(), muvec2.Phi(), weight);
          hPtLMComp -> Fill(m1_parton.Pt(), muvec1.Pt(), weight);
          hEtaLMComp -> Fill(m1_parton.Eta(), muvec1.Eta(), weight);
          hPhiLMComp -> Fill(m1_parton.Phi(), muvec1.Phi(), weight);
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
      if (fill_reco) {
        // fill b comparison histograms
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
        if (ghost_btag(branchGenParticle, genjet) && abs(genjet->Eta) < eta_cut) {
          bjets += 1;
          if (bjets == 1) {
            b1_particle = genjet->P4();
          } else if (bjets == 2) {
            b2_particle = genjet->P4();
          } else break;
        }
      }
      if (bjets != 2) fill_particle = false;
    }
    if (fill_particle) {
      if (num_elec_particle == 2) {
      // fill electron particle histograms
        zparticlevec = e1_particle + e2_particle;
        hPtLEP -> Fill(e1_particle.Pt(), weight);
        hEtaLEP -> Fill(e1_particle.Eta(), weight);
        hPhiLEP -> Fill(e1_particle.Phi(), weight);
        hPtSLEP -> Fill(e2_particle.Pt(), weight);
        hEtaSLEP -> Fill(e2_particle.Eta(), weight);
        hPhiSLEP -> Fill(e2_particle.Phi(), weight);
        // fill electron particle-reco comparison histograms
        if (fill_reco && num_elec_particle == 2) {
          hPtLECompP -> Fill(e1_particle.Pt(), elecvec1.Pt(), weight);
          hEtaLECompP -> Fill(e1_particle.Eta(), elecvec1.Eta(), weight);
          hPtSLECompP -> Fill(e2_particle.Pt(), elecvec2.Pt(), weight);
          hEtaSLECompP -> Fill(e2_particle.Eta(), elecvec2.Eta(), weight);
        }
      } else if (num_mu_particle == 2){
        // fill muon particle histograms
        zparticlevec = m1_particle + m2_particle;
        hPtLMP -> Fill(m1_particle.Pt(), weight);
        hEtaLMP -> Fill(m1_particle.Eta(), weight);
        hPhiLMP -> Fill(m1_particle.Phi(), weight);
        hPtSLMP -> Fill(m2_particle.Pt(), weight);
        hEtaSLMP -> Fill(m2_particle.Eta(), weight);
        hPhiSLMP -> Fill(m2_particle.Phi(), weight);
        // fill muon particle-reco comparison histograms
        if (fill_reco && num_mu_particle == 2) {
          hPtSLMCompP -> Fill(m2_particle.Pt(), muvec2.Pt(), weight);
          hEtaSLMCompP -> Fill(m2_particle.Eta(), muvec2.Eta(), weight);
          hPtLMCompP -> Fill(m1_particle.Pt(), muvec1.Pt(), weight);
          hEtaLMCompP -> Fill(m1_particle.Eta(), muvec1.Eta(), weight);
        }
      }
      // fill particle level b jet histograms
      hPtLBP -> Fill(b1_particle.Pt(), weight);
      hPtSLBP -> Fill(b2_particle.Pt(), weight);
      hEtaLBP -> Fill(b1_particle.Eta(), weight);
      hEtaSLBP -> Fill(b2_particle.Eta(), weight);
      hPhiLBP -> Fill(b1_particle.Phi(), weight);
      hPhiSLBP -> Fill(b2_particle.Phi(), weight);
      hparticlevec = b1_particle + b2_particle;
      // fill b jet comparison histograms
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
      if (fill_reco) {
        hPtZHCompP -> Fill((zparticlevec+hparticlevec).Pt(), sysvec.Pt(), weight);
        hEtaZHCompP -> Fill((zparticlevec+hparticlevec).Eta(), sysvec.Eta(), weight);
      }
    }
    nPassed+=weight;
    nPassedRaw++;
  }

  // Fit higgs mass distribution (can also do Z but W is less insightful due to MET)
  TF1 *fitMHR = new TF1("higgs_dscb_reco", DSCBf, 20, 250, 7);
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
  std::cout<<"chi2 for Z particle fit "<< z_results_particle->Chi2()<<std::endl;

  // Write histograms
  TFile *hists= new TFile("zhbb_hists.root","recreate");
  hists->cd();
  //mass
  hMZH -> Write();
  hMbbR->Write();
  hMllR->Write();
  hMbbllR -> Write();
  hMbbP->Write();
  hMllP->Write();
  hMbbllP -> Write();
  hMbbComp -> Write();
  hMllComp -> Write();
  hMZHComp -> Write();
  hMZHCompP -> Write();
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
  // phi
  hPhiH->Write();
  hPhiZ->Write();
  hPhiLE->Write();
  hPhiLM->Write();
  hPhiSLE->Write();
  hPhiSLM->Write();
  hPhiLB->Write();
  hPhiSLB->Write();
  hPhiZH->Write();
  hPhiBBR->Write();
  hPhiLLR->Write();
  hPhiLER->Write();
  hPhiLMR->Write();
  hPhiSLER->Write();
  hPhiSLMR->Write();
  hPhiLBR->Write();
  hPhiSLBR->Write();
  hPhiBBLLR->Write();
  hPhiBBP->Write();
  hPhiLLP->Write();
  hPhiLEP->Write();
  hPhiLMP->Write();
  hPhiSLEP->Write();
  hPhiSLMP->Write();
  hPhiLBP->Write();
  hPhiSLBP->Write();
  hPhiBBLLP->Write();
  hPhiHComp->Write();
  hPhiZComp->Write();
  hPhiLEComp->Write();
  hPhiLMComp->Write();
  hPhiSLEComp->Write();
  hPhiSLMComp->Write();
  hPhiLBComp->Write();
  hPhiSLBComp->Write();
  hPhiZHComp->Write();
  hists->Close();
  // draw histograms
  gROOT->SetBatch(kTRUE);
  // mass
  draw_hist(hMZH, "mass_ZH", "m_{T}^{ZH}", "mass (GeV)");
  draw_hist(hMbbR, "mass_bb_reco", "Reco m_{bb}", "mass (GeV)");
  draw_hist(hMllR, "mass_ll_reco", "Reco m_{ll}", "mass (GeV)");
  draw_hist(hMbbllR, "mass_bbll_reco", "Reco m_{bbll}", "mass (GeV)");
  draw_hist(hMbbP, "mass_bb_particle", "m_{bb}", "mass (GeV)");
  draw_hist(hMllP, "mass_ll_particle", "m_{ll}", "mass (GeV)");
  draw_hist(hMbbllP, "mass_bbll_particle", "m_{bbll}", "mass (GeV)");
  draw_hist2(hMbbComp, "mass_bb_Comp", "m_{bb}", "Particle Level mass (GeV)", "Reco Level mass (GeV)");
  draw_hist2(hMllComp, "mass_ll_Comp", "m_{ll}", "Particle Level mass (GeV)", "Reco Level mass (GeV)");
  draw_hist2(hMZHComp, "mass_ZH_Comp", "m_{T}^{ZH}", "Parton Level mass (GeV)", "Reco Level mass (GeV)");
  draw_hist2(hMZHCompP, "mass_ZH_Comp_particle", "m_{T}^{ZH}", "Particle Level mass (GeV)", "Reco Level mass (GeV)");
  // pt
  draw_hist(hPtH, "pt_H", "p_{T}^{H}", "p_{T}");
  draw_hist(hPtZ, "pt_Z", "p_{T}^{Z}", "p_{T}");
  draw_hist(hPtLE, "pt_LE",  "Lead p_{T}^{e}", "p_{T}");
  draw_hist(hPtLM, "pt_LM",  "Lead p_{T}^{#mu}", "p_{T}");
  draw_hist(hPtSLE, "pt_SLE",  "Sublead p_{T}^{e}", "p_{T}");
  draw_hist(hPtSLM, "pt_SLM",  "Sublead p_{T}^{#mu}", "p_{T}");
  draw_hist(hPtLB, "pt_LB", "Lead b Parton p_{T}", "p_{T}");
  draw_hist(hPtSLB, "pt_SLB", "Sublead b Parton p_{T}", "p_{T}");
  draw_hist(hPtZH, "pt_ZH", "p_{T}^{ZH}", "p_{T}");
  draw_hist(hPtBBR, "pt_bb_reco", "Reco p_{T}^{bb}", "p_{T}");
  draw_hist(hPtLLR, "pt_ll_reco", "Reco p_{T}^{ll}", "p_{T}");
  draw_hist(hPtLER, "pt_LE_reco",  "Lead Reco p_{T}^{e}", "p_{T}");
  draw_hist(hPtLMR, "pt_LM_reco",  "Lead Reco p_{T}^{#mu}", "p_{T}");
  draw_hist(hPtSLER, "pt_SLE_reco",  "Sublead Reco p_{T}^{e}", "p_{T}");
  draw_hist(hPtSLMR, "pt_SLM_reco",  "Sublead Reco p_{T}^{#mu}", "p_{T}");
  draw_hist(hPtLBR, "pt_LB_reco", "Lead Reco b Jet p_{T}", "p_{T}");
  draw_hist(hPtSLBR, "pt_SLB_reco", "Sublead Reco b Jet p_{T}", "p_{T}");
  draw_hist(hPtBBLLR, "pt_bbll_reco", "Reco p_{T}^{bbll}", "p_{T}");
  draw_hist(hPtBBP, "pt_bb_particle", "p_{T}^{bb}", "p_{T}");
  draw_hist(hPtLLP, "pt_ll_particle", "p_{T}^{ll}", "p_{T}");
  draw_hist(hPtLEP, "pt_LE_particle",  "Lead p_{T}^{e}", "p_{T}");
  draw_hist(hPtLMP, "pt_LM_particle",  "Lead p_{T}^{#mu}", "p_{T}");
  draw_hist(hPtSLEP, "pt_SLE_particle",  "Sublead p_{T}^{e}", "p_{T}");
  draw_hist(hPtSLMP, "pt_SLM_particle",  "Sublead p_{T}^{#mu}", "p_{T}");
  draw_hist(hPtLBP, "pt_LB_particle", "Lead b Jet p_{T}", "p_{T}");
  draw_hist(hPtSLBP, "pt_SLB_particle", "Sublead b Jet p_{T}", "p_{T}");
  draw_hist(hPtBBLLP, "pt_bbll_particle", "p_{T}^{bbll}", "p_{T}");
  // pt comparison
  draw_hist2(hPtHComp, "pt_H_Comp", "p_{T}^{H}", "Parton Level P_{T} (GeV)", "Reco Level P_{T} (GeV)");
  draw_hist2(hPtZComp, "pt_Z_Comp", "p_{T}^{Z}", "Parton Level P_{T} (GeV)", "Reco Level P_{T} (GeV)");
  draw_hist2(hPtLEComp, "pt_LE_Comp", "Lead p_{T}^{e}", "Parton Level P_{T} (GeV)", "Reco Level P_{T} (GeV)");
  draw_hist2(hPtLMComp, "pt_LM_Comp", "Lead p_{T}^{#mu}", "Parton Level P_{T} (GeV)", "Reco Level P_{T} (GeV)");
  draw_hist2(hPtSLEComp, "pt_SLE_Comp", "Sublead p_{T}^{e}", "Parton Level P_{T} (GeV)", "Reco Level P_{T} (GeV)");
  draw_hist2(hPtSLMComp, "pt_SLM_Comp", "Sublead p_{T}^{#mu}", "Parton Level P_{T} (GeV)", "Reco Level P_{T} (GeV)");
  draw_hist2(hPtLBComp, "pt_LB_Comp", "Lead p_{T}^{b}", "Parton Level P_{T} (GeV)", "Reco Level P_{T} (GeV)");
  draw_hist2(hPtSLBComp, "pt_SLB_Comp", "Sublead p_{T}^{b}", "Parton Level P_{T} (GeV)", "Reco Level P_{T} (GeV)");
  draw_hist2(hPtZHComp, "pt_ZH_Comp", "Sublead p_{T}^{ZH}", "Parton Level P_{T} (GeV)", "Reco Level P_{T} (GeV)");
  draw_hist2(hPtHCompP, "pt_H_Comp_particle", "p_{T}^{H}", "Particle Level P_{T} (GeV)", "Reco Level P_{T} (GeV)");
  draw_hist2(hPtZCompP, "pt_Z_Comp_particle", "p_{T}^{Z}", "Particle Level P_{T} (GeV)", "Reco Level P_{T} (GeV)");
  draw_hist2(hPtLECompP, "pt_LE_Comp_particle", "Lead p_{T}^{e}", "Particle Level P_{T} (GeV)", "Reco Level P_{T} (GeV)");
  draw_hist2(hPtLMCompP, "pt_LM_Comp_particle", "Lead p_{T}^{#mu}", "Particle Level P_{T} (GeV)", "Reco Level P_{T} (GeV)");
  draw_hist2(hPtSLECompP, "pt_SLE_Comp_particle", "Sublead p_{T}^{e}", "Particle Level P_{T} (GeV)", "Reco Level P_{T} (GeV)");
  draw_hist2(hPtSLMCompP, "pt_SLM_Comp_particle", "Sublead p_{T}^{#mu}", "Particle Level P_{T} (GeV)", "Reco Level P_{T} (GeV)");
  draw_hist2(hPtLBCompP, "pt_LB_Comp_particle", "Lead p_{T}^{b}", "Particle Level P_{T} (GeV)", "Reco Level P_{T} (GeV)");
  draw_hist2(hPtSLBCompP, "pt_SLB_Comp_particle", "Sublead p_{T}^{b}", "Particle Level P_{T} (GeV)", "Reco Level P_{T} (GeV)");
  draw_hist2(hPtZHCompP, "pt_ZH_Comp_particle", "Sublead p_{T}^{ZH}", "Particle Level P_{T} (GeV)", "Reco Level P_{T} (GeV)");
  // eta
  draw_hist(hEtaH, "eta_H", "#eta_{H}", "#eta");
  draw_hist(hEtaZ, "eta_Z", "#eta_{Z}", "#eta");
  draw_hist(hEtaLE, "eta_LE",  "Lead #eta_{e}", "#eta");
  draw_hist(hEtaLM, "eta_LM",  "Lead #eta_{#mu}", "#eta");
  draw_hist(hEtaSLE, "eta_SLE",  "Sublead #eta_{e}", "#eta");
  draw_hist(hEtaSLM, "eta_SLM",  "Sublead #eta_{#mu}", "#eta");
  draw_hist(hEtaLB, "eta_LB", "Lead #eta_{b}", "#eta");
  draw_hist(hEtaSLB, "eta_SLB", "Sublead #eta_{b}", "#eta");
  draw_hist(hEtaZH, "eta_ZH", "#eta_{ZH}", "#eta");
  draw_hist(hEtaBBR, "eta_bb_reco", "Reco #eta_{bb}", "#eta");
  draw_hist(hEtaLLR, "eta_ll_reco", "Reco #eta_{ll}", "#eta");
  draw_hist(hEtaLER, "eta_LE_reco",  "Lead Reco #eta_{e}", "#eta");
  draw_hist(hEtaLMR, "eta_LM_reco",  "Lead Reco #eta_{#mu}", "#eta");
  draw_hist(hEtaSLER, "eta_SLE_reco",  "Sublead Reco #eta_{e}", "#eta");
  draw_hist(hEtaSLMR, "eta_SLM_reco",  "Sublead Reco #eta_{#mu}", "#eta");
  draw_hist(hEtaLBR, "eta_LB_reco", "Lead Reco b Jet #eta", "#eta");
  draw_hist(hEtaSLBR, "eta_SLB_reco", "Sublead Reco b Jet #eta", "#eta");
  draw_hist(hEtaBBLLR, "eta_bbll_reco", "Reco #eta_{bbll}", "#eta");
  draw_hist(hEtaBBP, "eta_bb_particle", "#eta_{bb}", "#eta");
  draw_hist(hEtaLLP, "eta_ll_particle", "#eta_{ll}", "#eta");
  draw_hist(hEtaLEP, "eta_LE_particle",  "Lead #eta_{e}", "#eta");
  draw_hist(hEtaLMP, "eta_LM_particle",  "Lead #eta_{#mu}", "#eta");
  draw_hist(hEtaSLEP, "eta_SLE_particle",  "Sublead #eta_{e}", "#eta");
  draw_hist(hEtaSLMP, "eta_SLM_particle",  "Sublead #eta_{#mu}", "#eta");
  draw_hist(hEtaLBP, "eta_LB_particle", "Lead b Jet #eta", "#eta");
  draw_hist(hEtaSLBP, "eta_SLB_particle", "Sublead b Jet #eta", "#eta");
  draw_hist(hEtaBBLLP, "eta_bbll_particle", "#eta_{bbll}", "#eta");
  // eta comparison
  draw_hist2(hEtaHComp, "eta_H_Comp", "#eta_{H}", "Parton Level #eta", "Reco Level #eta");
  draw_hist2(hEtaZComp, "eta_Z_Comp", "#eta_{Z}", "Parton Level #eta", "Reco Level #eta");
  draw_hist2(hEtaLEComp, "eta_LE_Comp", "Lead #eta_{e}","Parton Level #eta", "Reco Level #eta");
  draw_hist2(hEtaLMComp, "eta_LM_Comp", "Lead #eta_{#mu}","Parton Level #eta", "Reco Level #eta");
  draw_hist2(hEtaSLEComp, "eta_SLE_Comp", "Sublead #eta_{e}","Parton Level #eta", "Reco Level #eta");
  draw_hist2(hEtaSLMComp, "eta_SLM_Comp", "Sublead #eta_{#mu}","Parton Level #eta", "Reco Level #eta");
  draw_hist2(hEtaLBComp, "eta_LB_Comp", "Lead #eta_{b}","Parton Level #eta", "Reco Level #eta");
  draw_hist2(hEtaSLBComp, "eta_SLB_Comp", "Sublead #eta_{b}","Parton Level #eta", "Reco Level #eta");
  draw_hist2(hEtaZHComp, "eta_ZH_Comp", "Sublead #eta_{ZH}", "Parton Level #eta", "Reco Level #eta");
  draw_hist2(hEtaHCompP, "eta_H_Comp_particle", "#eta_{H}", "Particle Level #eta", "Reco Level #eta");
  draw_hist2(hEtaZCompP, "eta_Z_Comp_particle", "#eta_{Z}", "Particle Level #eta", "Reco Level #eta");
  draw_hist2(hEtaLECompP, "eta_LE_Comp_particle", "Lead #eta_{e}","Particle Level #eta", "Reco Level #eta");
  draw_hist2(hEtaLMCompP, "eta_LM_Comp_particle", "Lead #eta_{#mu}","Particle Level #eta", "Reco Level #eta");
  draw_hist2(hEtaSLECompP, "eta_SLE_Comp_particle", "Sublead #eta_{e}","Particle Level #eta", "Reco Level #eta");
  draw_hist2(hEtaSLMCompP, "eta_SLM_Comp_particle", "Sublead #eta_{#mu}","Particle Level #eta", "Reco Level #eta");
  draw_hist2(hEtaLBCompP, "eta_LB_Comp_particle", "Lead #eta_{b}","Particle Level #eta", "Reco Level #eta");
  draw_hist2(hEtaSLBCompP, "eta_SLB_Comp_particle", "Sublead #eta_{b}","Particle Level #eta", "Reco Level #eta");
  draw_hist2(hEtaZHCompP, "eta_ZH_Comp_particle", "Sublead #eta_{ZH}", "Particle Level #eta", "Reco Level #eta");
  // phi parton level
  draw_hist(hPhiH, "phi_H", "#phi_{H}", "#phi (Rad)");
  draw_hist(hPhiZ, "phi_Z", "#phi_{Z}", "#phi (Rad)");
  draw_hist(hPhiLE, "phi_LE",  "Lead #phi_{e}", "#phi (Rad)");
  draw_hist(hPhiLM, "phi_LM",  "Lead #phi_{#mu}", "#phi (Rad)");
  draw_hist(hPhiSLE, "phi_SLE",  "Sublead #phi_{e}", "#phi (Rad)");
  draw_hist(hPhiSLM, "phi_SLM",  "Sublead #phi_{#mu}", "#phi (Rad)");
  draw_hist(hPhiLB, "phi_LB", "Lead #phi_{b}", "#phi (Rad)");
  draw_hist(hPhiSLB, "phi_SLB", "Sublead #phi_{b}", "#phi (Rad)");
  draw_hist(hPhiZH, "phi_ZH", "#phi_{ZH}", "#phi (Rad)");
  draw_hist(hPhiBBR, "phi_bb_reco", "Reco #phi_{bb}", "#phi (Rad)");
  draw_hist(hPhiLLR, "phi_ll_reco", "Reco #phi_{ll}", "#phi (Rad)");
  draw_hist(hPhiLER, "phi_LE_reco",  "Lead Reco #phi_{e}", "#phi (Rad)");
  draw_hist(hPhiLMR, "phi_LM_reco",  "Lead Reco #phi_{#mu}", "#phi (Rad)");
  draw_hist(hPhiSLER, "phi_SLE_reco",  "Sublead Reco #phi_{e}", "#phi (Rad)");
  draw_hist(hPhiSLMR, "phi_SLM_reco",  "Sublead Reco #phi_{#mu}", "#phi (Rad)");
  draw_hist(hPhiLBR, "phi_LB_reco", "Lead Reco b Jet #phi", "#phi (Rad)");
  draw_hist(hPhiSLBR, "phi_SLB_reco", "Sublead Reco b Jet #phi", "#phi (Rad)");
  draw_hist(hPhiBBLLR, "phi_bbll_reco", "Reco #phi_{bbll}", "#phi (Rad)");
  draw_hist(hPhiBBP, "phi_bb_particle", "#phi_{bb}", "#phi (Rad)");
  draw_hist(hPhiLLP, "phi_ll_particle", "#phi_{ll}", "#phi (Rad)");
  draw_hist(hPhiLEP, "phi_LE_particle",  "Lead #phi_{e}", "#phi (Rad)");
  draw_hist(hPhiLMP, "phi_LM_particle",  "Lead #phi_{#mu}", "#phi (Rad)");
  draw_hist(hPhiSLEP, "phi_SLE_particle",  "Sublead #phi_{e}", "#phi (Rad)");
  draw_hist(hPhiSLMP, "phi_SLM_particle",  "Sublead #phi_{#mu}", "#phi (Rad)");
  draw_hist(hPhiLBP, "phi_LB_particle", "Lead b Jet #phi", "#phi (Rad)");
  draw_hist(hPhiSLBP, "phi_SLB_particle", "Sublead b Jet #phi", "#phi (Rad)");
  draw_hist(hPhiBBLLP, "phi_bbll_particle", "#phi_{bbll}", "#phi (Rad)");
  // phi comparison
  draw_hist2(hPhiHComp, "phi_H_Comp", "#phi_{H}", "Parton Level #phi (Rad)", "Reco Level #phi (Rad)");
  draw_hist2(hPhiZComp, "phi_Z_Comp", "#phi_{Z}", "Parton Level #phi (Rad)", "Reco Level #phi (Rad)");
  draw_hist2(hPhiLEComp, "phi_LE_Comp", "Lead #phi_{e}","Parton Level #phi (Rad)", "Reco Level #phi (Rad)");
  draw_hist2(hPhiLMComp, "phi_LM_Comp", "Lead #phi_{#mu}","Parton Level #phi (Rad)", "Reco Level #phi (Rad)");
  draw_hist2(hPhiSLEComp, "phi_SLE_Comp", "Sublead #phi_{e}","Parton Level #phi (Rad)", "Reco Level #phi (Rad)");
  draw_hist2(hPhiSLMComp, "phi_SLM_Comp", "Sublead #phi_{#mu}","Parton Level #phi (Rad)", "Reco Level #phi (Rad)");
  draw_hist2(hPhiLBComp, "phi_LB_Comp", "Lead #phi_{b}","Parton Level #phi (Rad)", "Reco Level #phi (Rad)");
  draw_hist2(hPhiSLBComp, "phi_SLB_Comp", "Sublead #phi_{b}","Parton Level #phi (Rad)", "Reco Level #phi (Rad)");
  draw_hist2(hPhiZHComp, "phi_ZH_Comp", "Sublead #phi_{ZH}", "Parton Level #phi (Rad)", "Reco Level #phi (Rad)");
  // print some info
  std::cout<<"Integral of m_H hist "<< hMbbR->Integral()<<" integral of m_Z; "<<hMllR->Integral()<<std::endl;
  std::cout<<"Total number of entries "<<numberOfEntries<<" Passed "<<nPassed<<" raw "<<nPassedRaw<<std::endl;
}
