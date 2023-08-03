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
    c->Print(Form("whbb_plots/%s.%s",name.c_str(),(*it).c_str()),(*it).c_str());
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

// search for the final daughter given a particle and the id of the daughter to search for
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

// make a ton of plots from whbb events (w -> l v)
void whbb_analyze(const char *inputFilePlus, const char *inputFileMinus) {
  gSystem->Load("libDelphes");
  // Create chain of root trees
  TChain chain("Delphes");
  chain.Add(inputFilePlus);
  // Create object of class ExRootTreeReader
  ExRootTreeReader *treeReader = new ExRootTreeReader(&chain);
  Long64_t numberOfPlusEntries = treeReader->GetEntries();
  chain.Add(inputFileMinus);
  Long64_t numberOfEntries = treeReader->GetEntries();
  // Get pointers to branches used in this analysis
  TClonesArray *branchJet = treeReader->UseBranch("Jet");
  TClonesArray *branchElectron = treeReader->UseBranch("Electron");
  TClonesArray *branchMuon = treeReader->UseBranch("Muon");
  TClonesArray *branchMissingET = treeReader->UseBranch("MissingET");
  TClonesArray *branchEvent = treeReader->UseBranch("Event");
  TClonesArray *branchGenParticle = treeReader->UseBranch("Particle");
  TClonesArray *branchGenJet = treeReader->UseBranch("GenJet");
  TClonesArray *branchGenMissingET = treeReader->UseBranch("GenMissingET");
  // SET CUTS
  const double e_pt_cut = 34;
  const double mu_pt_cut = 30;
  const double eta_cut = 2.5;
  // Book histograms
  // mass
  TH1 *hMTW = new TH1F("massT_W", "m_{T}^{W}", 100, 40.0, 140.0);
  TH1 *hMTWH = new TH1F("massT_WH", "m_{T}^{WH}", 100, 120.0, 900.0);
  TH1 *hMbbR = new TH1F("mass_bb_reco", "Reco m_{bb}", 100, 0.0, 250.0);
  TH1 *hMTlvR = new TH1F("massT_lv_reco", "Reco m_{T}^{l#nu}", 100, 40.0, 300.0);
  TH1 *hMTbblvR = new TH1F("massT_bblv_reco", "Reco m_{T}^{bbl#nu}", 100, 120.0, 900.0);
  TH1 *hMbbP = new TH1F("mass_bb_particle", "m_{bb}", 100, 0.0, 250.0);
  TH1 *hMTlvP = new TH1F("massT_lv_particle", "m_{T}^{l#nu}", 100, 40.0, 300.0);
  TH1 *hMTbblvP = new TH1F("massT_bblv_particle", "m_{T}^{bbl#nu}", 100, 120.0, 900.0);
  TH2 *hMbbComp = new TH2F("mass_bb_Comp", "m_{bb}", 50, 0.0, 250.0, 50, 0.0, 250.0);
  TH2 *hMTWComp = new TH2F("massT_W_Comp", "m_{T}^{W}", 50, 40.0, 300.0, 50, 40.0, 300.0);
  TH2 *hMTWHComp = new TH2F("massT_WH_Comp", "m_{T}^{WH}", 50, 120.0, 900.0, 50, 120.0, 900.0);
  TH2 *hMTWCompP = new TH2F("massT_W_Comp_particle", "m_{T}^{W}", 50, 50.0, 600.0, 50, 50.0, 600.0);
  TH2 *hMTWHCompP = new TH2F("massT_WH_Comp_particle", "m_{T}^{WH}", 50, 120.0, 900.0, 50, 120.0, 900.0);
  // pt parton level
  TH1 *hPtH = new TH1F("pt_H", "p_{T}^{H}", 100, 0.0, 400.0);
  TH1 *hPtW = new TH1F("pt_W", "p_{T}^{W}", 100, 0.0, 400.0);
  TH1 *hPtE = new TH1F("pt_e",  "p_{T}^{e}", 100, 20.0, 200.0);
  TH1 *hPtNu = new TH1F("pt_nu",  "p_{T}^{MET}", 100, 20.0, 200.0);
  TH1 *hPtMu = new TH1F("pt_mu",  "p_{T}^{#mu}", 100, 20.0, 200.0);
  TH1 *hPtLB = new TH1F("pt_LB", "Lead b Parton p_{T}", 100, 20.0, 200.0);
  TH1 *hPtSLB = new TH1F("pt_SLB", "Sublead b Parton p_{T}", 100, 20.0, 200.0);
  TH1 *hPtWH = new TH1F("pt_WH", "p_{T}^{WH}", 100, 0.0, 400.0);
  // pt reco level
  TH1 *hPtBBR = new TH1F("pt_bb_reco", "Reco p_{T}^{bb}", 100, 0.0, 400.0);
  TH1 *hPtLVR = new TH1F("pt_lv_reco", "Reco p_{T}^{l#nu}", 100, 0.0, 400.0);
  TH1 *hPtER = new TH1F("pt_e_reco",  "Reco p_{T}^{e}", 100, 20.0, 200.0);
  TH1 *hPtMetR = new TH1F("pt_met_reco",  "Reco p_{T}^{MET}", 100, 20.0, 200.0);
  TH1 *hPtMuR = new TH1F("pt_mu_reco",  "Reco p_{T}^{#mu}", 100, 20.0, 200.0);
  TH1 *hPtLBR = new TH1F("pt_LB_reco", "Lead Reco b Jet p_{T}", 100, 20.0, 200.0);
  TH1 *hPtSLBR = new TH1F("pt_SLB_reco", "Sublead Reco b Jet p_{T}", 100, 20.0, 200.0);
  TH1 *hPtBBLVR = new TH1F("pt_bblv_reco", "Reco p_{T}^{bbl#nu}", 100, 0.0, 400.0);
  // pt particle level
  TH1 *hPtBBP = new TH1F("pt_bb_particle", "p_{T}^{bb}", 100, 0.0, 400.0);
  TH1 *hPtLVP = new TH1F("pt_lv_particle", "p_{T}^{l#nu}", 100, 0.0, 400.0);
  TH1 *hPtEP = new TH1F("pt_e_particle",  "p_{T}^{e}", 100, 20.0, 200.0);
  TH1 *hPtMetP = new TH1F("pt_met_particle", "p_{T}^{MET}", 100, 20.0, 200.0);
  TH1 *hPtMuP = new TH1F("pt_mu_particle",  "p_{T}^{#mu}", 100, 20.0, 200.0);
  TH1 *hPtLBP = new TH1F("pt_LB_particle", "Lead b Jet p_{T}", 100, 20.0, 200.0);
  TH1 *hPtSLBP = new TH1F("pt_SLB_particle", "Sublead b Jet p_{T}", 100, 20.0, 200.0);
  TH1 *hPtBBLVP = new TH1F("pt_bblv_particle", "p_{T}^{bbl#nu}", 100, 0.0, 400.0);
  // pt parton-reco comparison
  TH2 *hPtHComp = new TH2F("pt_H_Comp", "p_{T}^{H}", 50, 0.0, 400.0, 50, 0.0, 400.0);
  TH2 *hPtWComp = new TH2F("pt_W_Comp", "p_{T}^{W}", 50, 0.0, 400.0, 50, 0.0, 400.0);
  TH2 *hPtEComp = new TH2F("pt_e_Comp", "p_{T}^{e}", 50, 20.0, 200.0, 50, 20.0, 200.0);
  TH2 *hPtMetComp = new TH2F("pt_nu_Comp", "p_{T}^{#nu}", 50, 20.0, 200.0, 50, 20.0, 200.0);
  TH2 *hPtMuComp = new TH2F("pt_mu_Comp", "p_{T}^{#mu}", 50, 20.0, 200.0, 50, 20.0, 200.0);
  TH2 *hPtLBComp = new TH2F("pt_LB_Comp", "Lead p_{T}^{b}", 50, 20.0, 200.0, 50, 20.0, 200.0);
  TH2 *hPtSLBComp = new TH2F("pt_SLB_Comp", "Sublead p_{T}^{b}", 50, 20.0, 200.0, 50, 20.0, 200.0);
  TH2 *hPtWHComp = new TH2F("pt_WH_Comp", "Sublead p_{T}^{WH}", 50, 0.0, 400.0, 50, 0.0, 400.0);
  // pt particle-reco comparison
  TH2 *hPtHCompP = new TH2F("pt_H_Comp_particle", "p_{T}^{H}", 50, 0.0, 400.0, 50, 0.0, 400.0);
  TH2 *hPtWCompP = new TH2F("pt_W_Comp_particle", "p_{T}^{W}", 50, 0.0, 400.0, 50, 0.0, 400.0);
  TH2 *hPtECompP = new TH2F("pt_e_Comp_particle", "p_{T}^{e}", 50, 20.0, 200.0, 50, 20.0, 200.0);
  TH2 *hPtMetCompP = new TH2F("pt_met_Comp_particle", "MET", 50, 20.0, 200.0, 50, 20.0, 200.0);
  TH2 *hPtMuCompP = new TH2F("pt_mu_Comp_particle", "p_{T}^{#mu}", 50, 20.0, 200.0, 50, 20.0, 200.0);
  TH2 *hPtLBCompP = new TH2F("pt_LB_Comp_particle", "Lead p_{T}^{b}", 50, 20.0, 200.0, 50, 20.0, 200.0);
  TH2 *hPtSLBCompP = new TH2F("pt_SLB_Comp_particle", "Sublead p_{T}^{b}", 50, 20.0, 200.0, 50, 20.0, 200.0);
  TH2 *hPtWHCompP = new TH2F("pt_WH_Comp_particle", "Sublead p_{T}^{WH}", 50, 0.0, 400.0, 50, 0.0, 400.0);
  // eta parton level
  TH1 *hEtaH = new TH1F("eta_H", "#eta_{H}", 100, -3.0, 3.0);
  TH1 *hEtaW = new TH1F("eta_W", "#eta_{W}", 100, -3.0, 3.0);
  TH1 *hEtaE = new TH1F("eta_e",  "#eta_{e}", 100, -3.0, 3.0);
  TH1 *hEtaNu = new TH1F("eta_nu",  "#eta_{MET}", 100, -3.0, 3.0);
  TH1 *hEtaMu = new TH1F("eta_mu",  "#eta_{#mu}", 100, -3.0, 3.0);
  TH1 *hEtaLB = new TH1F("eta_LB", "Lead #eta_{b}", 100, -3.0, 3.0);
  TH1 *hEtaSLB = new TH1F("eta_SLB", "Sublead #eta_{b}", 100, -3.0, 3.0);
  TH1 *hEtaWH = new TH1F("eta_WH", "#eta_{WH}", 100, -3.0, 3.0);
  // eta reco level
  TH1 *hEtaBBR = new TH1F("eta_bb_reco", "Reco #eta_{bb}", 100, -3.0, 3.0);
  TH1 *hEtaLVR = new TH1F("eta_lv_reco", "Reco #eta_{l#nu}", 100, -3.0, 3.0);
  TH1 *hEtaER = new TH1F("eta_e_reco",  "Reco #eta_{e}", 100, -3.0, 3.0);
  TH1 *hEtaMetR = new TH1F("eta_met_reco",  "Reco #eta_{MET}", 100, -3.0, 3.0);
  TH1 *hEtaMuR = new TH1F("eta_mu_reco",  "Reco #eta_{#mu}", 100, -3.0, 3.0);
  TH1 *hEtaLBR = new TH1F("eta_LB_reco", "Lead Reco b Jet #eta", 100, -3.0, 3.0);
  TH1 *hEtaSLBR = new TH1F("eta_SLB_reco", "Sublead Reco b Jet #eta", 100, -3.0, 3.0);
  TH1 *hEtaBBLVR = new TH1F("eta_bblv_reco", "Reco #eta_{bbl#nu}", 100, -3.0, 3.0);
  // eta particle level
  TH1 *hEtaBBP = new TH1F("eta_bb_particle", "#eta_{bb}", 100, -3.0, 3.0);
  TH1 *hEtaLVP = new TH1F("eta_lv_particle", "#eta_{l#nu}", 100, -3.0, 3.0);
  TH1 *hEtaEP = new TH1F("eta_e_particle",  "#eta_{e}", 100, -3.0, 3.0);
  TH1 *hEtaMetP = new TH1F("eta_met_particle",  "#eta_{MET}", 100, -3.0, 3.0);
  TH1 *hEtaMuP = new TH1F("eta_mu_particle",  "#eta_{#mu}", 100, -3.0, 3.0);
  TH1 *hEtaLBP = new TH1F("eta_LB_particle", "Lead b Jet #eta", 100, -3.0, 3.0);
  TH1 *hEtaSLBP = new TH1F("eta_SLB_particle", "Sublead b Jet #eta", 100, -3.0, 3.0);
  TH1 *hEtaBBLVP = new TH1F("eta_bblv_particle", "#eta_{bbl#nu}", 100, -3.0, 3.0);
  // eta parton-reco comparison
  TH2 *hEtaHComp = new TH2F("eta_H_Comp", "#eta_{H}", 50, -3.0, 3.0, 50, -3.0, 3.0);
  TH2 *hEtaWComp = new TH2F("eta_W_Comp", "#eta_{W}", 50, -3.0, 3.0, 50, -3.0, 3.0);
  TH2 *hEtaEComp = new TH2F("eta_e_Comp", "#eta_{e}", 50,-3.0, 3.0, 50, -3.0, 3.0);
  TH2 *hEtaMetComp = new TH2F("eta_nu_Comp", "#eta_{#nu}", 50,-3.0, 3.0, 50, -3.0, 3.0);
  TH2 *hEtaMuComp = new TH2F("eta_mu_Comp", "#eta_{#mu}", 50,-3.0, 3.0, 50, -3.0, 3.0);
  TH2 *hEtaLBComp = new TH2F("eta_LB_Comp", "Lead #eta_{b}", 50,-3.0, 3.0, 50, -3.0, 3.0);
  TH2 *hEtaSLBComp = new TH2F("eta_SLB_Comp", "Sublead #eta_{b}", 50,-3.0, 3.0, 50, -3.0, 3.0);
  TH2 *hEtaWHComp = new TH2F("eta_WH_Comp", "Sublead #eta_{WH}", 50, -3.0, 3.0, 50, -3.0, 3.0);
  // eta particle-reco comparison
  TH2 *hEtaHCompP = new TH2F("eta_H_Comp_particle", "#eta_{H}", 50, -3.0, 3.0, 50, -3.0, 3.0);
  TH2 *hEtaWCompP = new TH2F("eta_W_Comp_particle", "#eta_{W}", 50, -3.0, 3.0, 50, -3.0, 3.0);
  TH2 *hEtaECompP = new TH2F("eta_e_Comp_particle", "#eta_{e}", 50,-3.0, 3.0, 50, -3.0, 3.0);
  TH2 *hEtaMetCompP = new TH2F("eta_nu_Comp_particle", "#eta_{#nu}", 50,-3.0, 3.0, 50, -3.0, 3.0);
  TH2 *hEtaMuCompP = new TH2F("eta_mu_Comp_particle", "#eta_{#mu}", 50,-3.0, 3.0, 50, -3.0, 3.0);
  TH2 *hEtaLBCompP = new TH2F("eta_LB_Comp_particle", "Lead #eta_{b}", 50,-3.0, 3.0, 50, -3.0, 3.0);
  TH2 *hEtaSLBCompP = new TH2F("eta_SLB_Comp_particle", "Sublead #eta_{b}", 50,-3.0, 3.0, 50, -3.0, 3.0);
  TH2 *hEtaWHCompP = new TH2F("eta_WH_Comp_particle", "Sublead #eta_{WH}", 50, -3.0, 3.0, 50, -3.0, 3.0);
  // phi parton level
  TH1 *hPhiH = new TH1F("phi_H", "#phi_{H}", 100, -3.15, 3.15);
  TH1 *hPhiW = new TH1F("phi_W", "#phi_{W}", 100, -3.15, 3.15);
  TH1 *hPhiE = new TH1F("phi_e",  "#phi_{e}", 100, -3.15, 3.15);
  TH1 *hPhiNu = new TH1F("phi_nu",  "#phi_{MET}", 100, -3.15, 3.15);
  TH1 *hPhiMu = new TH1F("phi_mu",  "#phi_{#mu}", 100, -3.15, 3.15);
  TH1 *hPhiLB = new TH1F("phi_LB", "Lead #phi_{b}", 100, -3.15, 3.15);
  TH1 *hPhiSLB = new TH1F("phi_SLB", "Sublead #phi_{b}", 100, -3.15, 3.15);
  TH1 *hPhiWH = new TH1F("phi_WH", "#phi_{WH}", 100, -3.15, 3.15);
  // phi reco level
  TH1 *hPhiBBR = new TH1F("phi_bb_reco", "Reco #phi_{bb}", 100, -3.15, 3.15);
  TH1 *hPhiLVR = new TH1F("phi_lv_reco", "Reco #phi_{l#nu}", 100, -3.15, 3.15);
  TH1 *hPhiER = new TH1F("phi_e_reco",  "Reco #phi_{e}", 100, -3.15, 3.15);
  TH1 *hPhiMetR = new TH1F("phi_met_reco",  "Reco #phi_{MET}", 100, -3.15, 3.15);
  TH1 *hPhiMuR = new TH1F("phi_mu_reco",  "Reco #phi_{#mu}", 100, -3.15, 3.15);
  TH1 *hPhiLBR = new TH1F("phi_LB_reco", "Lead Reco b Jet #phi", 100, -3.15, 3.15);
  TH1 *hPhiSLBR = new TH1F("phi_SLB_reco", "Sublead Reco b Jet #phi", 100, -3.15, 3.15);
  TH1 *hPhiBBLVR = new TH1F("phi_bblv_reco", "Reco #phi_{bbl#nu}", 100, -3.15, 3.15);
  // phi particle level
  TH1 *hPhiBBP = new TH1F("phi_bb_particle", "#phi_{bb}", 100, -3.15, 3.15);
  TH1 *hPhiLVP = new TH1F("phi_lv_particle", "#phi_{l#nu}", 100, -3.15, 3.15);
  TH1 *hPhiEP = new TH1F("phi_e_particle",  "#phi_{e}", 100, -3.15, 3.15);
  TH1 *hPhiMetP = new TH1F("phi_met_particle",  "#phi_{MET}", 100, -3.15, 3.15);
  TH1 *hPhiMuP = new TH1F("phi_mu_particle",  "#phi_{#mu}", 100, -3.15, 3.15);
  TH1 *hPhiLBP = new TH1F("phi_LB_particle", "Lead b Jet #phi", 100, -3.15, 3.15);
  TH1 *hPhiSLBP = new TH1F("phi_SLB_particle", "Sublead b Jet #phi", 100, -3.15, 3.15);
  TH1 *hPhiBBLVP = new TH1F("phi_bblv_particle", "#phi_{bbl#nu}", 100, -3.15, 3.15);
  // phi comparison
  TH2 *hPhiHComp = new TH2F("phi_H_Comp", "#phi_{H}", 50, -3.15, 3.15, 50, -3.15, 3.15);
  TH2 *hPhiWComp = new TH2F("phi_W_Comp", "#phi_{W}", 50, -3.15, 3.15, 50, -3.15, 3.15);
  TH2 *hPhiEComp = new TH2F("phi_e_Comp", "#phi_{e}", 50,-3.15, 3.15, 50, -3.15, 3.15);
  TH2 *hPhiMetComp = new TH2F("phi_nu_Comp", "#phi_{#nu}", 50,-3.15, 3.15, 50, -3.15, 3.15);
  TH2 *hPhiMuComp = new TH2F("phi_mu_Comp", "#phi_{#mu}", 50,-3.15, 3.15, 50, -3.15, 3.15);
  TH2 *hPhiLBComp = new TH2F("phi_LB_Comp", "Lead #phi_{b}", 50,-3.15, 3.15, 50, -3.15, 3.15);
  TH2 *hPhiSLBComp = new TH2F("phi_SLB_Comp", "Sublead #phi_{b}", 50,-3.15, 3.15, 50, -3.15, 3.15);
  TH2 *hPhiWHComp = new TH2F("phi_WH_Comp", "Sublead #phi_{WH}", 50, -3.15, 3.15, 50, -3.15, 3.15);
  // initialize variables needed for filling histograms
  double  nPassed=0;
  int  nPassedRaw=0;
  double Lumi=3e3;
  int bjets = 0;
  bool elec_ev_parton;
  TLorentzVector b1_parton;
  TLorentzVector b2_parton;
  TLorentzVector b1_reco;
  TLorentzVector b2_reco;
  TLorentzVector b1_particle;
  TLorentzVector b2_particle;
  TLorentzVector higgsvec;
  TLorentzVector elecvec;
  TLorentzVector muvec;
  Electron *elec;
  Electron *temp_elec;
  Muon *muon;
  Muon *temp_muon;
  MissingET *met;
  MissingET *genmet;
  GenParticle *daughter1;
  GenParticle *daughter2;
  GenParticle *daughter;
  GenParticle *daughternu;
  TLorentzVector metvec;
  TLorentzVector genmetvec;
  TLorentzVector wvec;
  TLorentzVector sysvec;
  TLorentzVector wpartonvec;
  TLorentzVector wparticlevec;
  TLorentzVector hpartonvec;
  TLorentzVector hparticlevec;
  TLorentzVector eparticlevec;
  TLorentzVector mparticlevec;
  bool fill_reco = true;
  bool fill_parton = true;
  bool fill_particle = true;

  // Loop over events
  for(Int_t entry = 0; entry < numberOfEntries; ++entry){
    fill_reco = true;
    fill_parton = true;
    fill_particle = true;
    // Load selected branches with data from specified event
    treeReader->ReadEntry(entry);
    HepMCEvent *event = (HepMCEvent*) branchEvent -> At(0);
    Float_t weight = event->Weight/numberOfEntries*Lumi;
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
    // Now check for events with 1 charged lepton and save its information
    // Also combine it with MET data to reconstruct the W
    met = (MissingET *) branchMissingET->At(0);
    metvec = met ->P4();
    elec = NULL;
    muon = NULL;
    if (fill_reco) {
      // check for one electron meeting requirements
      if (branchElectron->GetEntries() > 0) {
        for(int i=0; i<(int)branchElectron->GetEntries(); i++) {
          temp_elec = (Electron *) branchElectron->At(i);
          if (abs(temp_elec->Eta) < eta_cut && temp_elec->PT > e_pt_cut) {
            if (elec == NULL) {
              elec = temp_elec;
            } else {
              fill_reco = false;
              break;
            }
          }
        }
      }
      // Check for one muon meeting requirements 
      if (branchMuon->GetEntries() > 0) {
        for(int i=0; i<(int)branchMuon->GetEntries(); i++) {
          temp_muon = (Muon *) branchMuon->At(i);
          if (abs(temp_muon->Eta) < eta_cut && temp_muon->PT > e_pt_cut) {
            if (muon == NULL) {
              muon = temp_muon;
            } else {
              fill_reco = false;
              break;
            }
          }
        }
      }
    }
    // fill electron histograms
    if (fill_reco && elec != NULL && muon == NULL) {
      elecvec = elec->P4();
      wvec = elecvec + metvec;
      hPtER -> Fill(elec->PT, weight);
      hEtaER -> Fill(elec->Eta, weight);
      hPhiER -> Fill(elec->Phi, weight);
    // fill muon histograms
    } else if (fill_reco && elec == NULL && muon != NULL) {
      muvec = muon->P4();
      wvec = muvec + metvec;
      hPtMuR -> Fill(muon->PT, weight);
      hEtaMuR -> Fill(muon->Eta, weight);
      hPhiMuR -> Fill(muon->Phi, weight);
    } else {
      fill_reco = false;
    }
    // fill other reco histograms
    sysvec = higgsvec + wvec;
    if (fill_reco) {
      // fill met reco histograms
      hPtMetR -> Fill(metvec.Pt(), weight);
      hEtaMetR -> Fill(metvec.Eta(), weight);
      hPhiMetR -> Fill(metvec.Phi(), weight);
      // fill mass reco histograms
      hMTlvR -> Fill(wvec.Mt(), weight);
      hMTbblvR -> Fill(sysvec.Mt(), weight);
      hMbbR -> Fill(higgsvec.M(), weight);
      // fill pt reco histograms 
      hPtBBR -> Fill(higgsvec.Pt(), weight);
      hPtLVR -> Fill(wvec.Pt(), weight);
      hPtLBR -> Fill(b1_reco.Pt(), weight);
      hPtSLBR -> Fill(b2_reco.Pt(), weight);
      hPtBBLVR -> Fill(sysvec.Pt(), weight);
      // fill eta reco histograms
      hEtaBBR -> Fill(higgsvec.Eta(), weight);
      hEtaLVR -> Fill(wvec.Eta(), weight);
      hEtaLBR -> Fill(b1_reco.Eta(), weight);
      hEtaSLBR -> Fill(b2_reco.Eta(), weight);
      hEtaBBLVR -> Fill(sysvec.Eta(), weight);
      // fill phi reco histograms
      hPhiBBR -> Fill(higgsvec.Phi(), weight);
      hPhiLVR -> Fill(wvec.Phi(), weight);
      hPhiLBR -> Fill(b1_reco.Phi(), weight);
      hPhiSLBR -> Fill(b2_reco.Phi(), weight);
      hPhiBBLVR -> Fill(sysvec.Phi(), weight);
    }
    // loop over true particles and fill those histograms
    b1_parton.SetPtEtaPhiM(0,0,0,0);
    b2_parton.SetPtEtaPhiM(0,0,0,0);
    eparticlevec.SetPtEtaPhiM(0,0,0,0);
    mparticlevec.SetPtEtaPhiM(0,0,0,0);
    genmet = (MissingET *) branchGenMissingET->At(0);
    genmetvec = genmet ->P4();
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
        // b partons
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
      // W parton
      } else if ((particle -> PID == 24 && d1_pid != 24 && d2_pid != 24) || (particle -> PID == -24 && d1_pid != -24 && d1_pid != -24)) {
        wpartonvec = particle->P4();
        // check for electron daughter 1
        if (abs(d1_pid) == 11 || abs(d2_pid) == 11) {
          elec_ev_parton = true;
          // parton level
          if (abs(d1_pid) == 11) {
            daughter = daughter1;
            daughternu = daughter2;
          } else {
            daughter = daughter2;
            daughternu = daughter1;
          }
          if (daughter->PT < e_pt_cut || abs(daughter->Eta) > eta_cut) fill_parton = false;
        // check for muon daughter 1
        } else if (abs(d1_pid) == 13 || abs(d2_pid) == 13) {
          if (abs(d1_pid) == 13) {
            daughter = daughter1;
            daughternu = daughter2;
          } else {
            daughter = daughter2;
            daughternu = daughter1;
          }
          if (daughter->PT < mu_pt_cut || abs(daughter->Eta) > eta_cut) fill_parton = false;
        }
      // check for status 1 (particle level) electron
      } else if (abs(particle -> PID) == 11 && particle -> Status == 1 && particle->PT > e_pt_cut && abs(particle->Eta) < eta_cut) {
        if (eparticlevec.Pt() == 0) {
          eparticlevec = particle -> P4();
          wparticlevec = genmetvec + eparticlevec;
        } else {
          fill_particle = false;
        }
      } else if (abs(particle -> PID) == 13 && particle -> Status == 1 && particle->PT > mu_pt_cut && abs(particle->Eta) < eta_cut) {
        if (mparticlevec.Pt() == 0) {
          mparticlevec = particle -> P4();
          wparticlevec = genmetvec + mparticlevec;
        } else {
          fill_particle = false;
        }
      }
    }
    // check only one charged lepton
    if ((eparticlevec.Pt() == 0 && mparticlevec.Pt() == 0) || (eparticlevec.Pt() != 0 && mparticlevec.Pt() != 0)) fill_particle = false;
    // fill parton level histograms
    if (fill_parton) {
      // fill higgs and w boson parton level histograms
      hPtH -> Fill(hpartonvec.Pt(), weight);
      hEtaH -> Fill(hpartonvec.Eta(), weight);
      hPhiH -> Fill(hpartonvec.Phi(), weight);
      hPtW -> Fill(wpartonvec.Pt(), weight);
      hEtaW -> Fill(wpartonvec.Eta(), weight);
      hPhiW -> Fill(wpartonvec.Phi(), weight);
      hMTW -> Fill(wpartonvec.Mt(), weight);
      if (fill_reco) {
        // fill higgs and w boson parton-reco comparison histograms
        hPtHComp -> Fill(hpartonvec.Pt(), higgsvec.Pt(), weight);
        hEtaHComp -> Fill(hpartonvec.Eta(), higgsvec.Eta(), weight);
        hPhiHComp -> Fill(hpartonvec.Phi(), higgsvec.Phi(), weight);
        hPtWComp -> Fill(wpartonvec.Pt(), wvec.Pt(), weight);
        hEtaWComp -> Fill(wpartonvec.Eta(), wvec.Eta(), weight);
        hPhiWComp -> Fill(wpartonvec.Phi(), wvec.Phi(), weight);
        hMTWComp -> Fill(wpartonvec.Mt(), sysvec.Mt(), weight);
      }
      // fill electron parton histograms
      if (elec_ev_parton) {
        hPtE -> Fill(daughter -> PT, weight);
        hEtaE -> Fill(daughter -> Eta, weight);
        hPhiE -> Fill(daughter -> Phi, weight);
        if (fill_reco && elec != NULL) {
          hPtEComp -> Fill(daughter -> PT, elecvec.Pt(), weight);
          hEtaEComp -> Fill(daughter -> Eta, elecvec.Eta(), weight);
          hPhiEComp -> Fill(daughter -> Phi, elecvec.Phi(), weight);
        }
      // fill muon parton histograms
      } else {
        hPtMu -> Fill(daughter -> PT, weight);
        hEtaMu -> Fill(daughter -> Eta, weight);
        hPhiMu -> Fill(daughter -> Phi, weight);
        if (fill_reco && elec == NULL) {
          hPtMuComp -> Fill(daughter -> PT, muvec.Pt(), weight);
          hEtaMuComp -> Fill(daughter -> Eta, muvec.Eta(), weight);
          hPhiMuComp -> Fill(daughter -> Phi, muvec.Phi(), weight);
        }
      }
      // fill neutrino parton histograms
      hPtNu -> Fill(daughternu -> PT, weight);
      hEtaNu -> Fill(daughternu -> Eta, weight);
      hPhiNu -> Fill(daughternu -> Phi, weight);
      if (fill_reco) {
        hPtMetComp -> Fill(daughternu -> PT, metvec.Pt(), weight);
        hEtaMetComp -> Fill(daughternu -> Eta, metvec.Eta(), weight);
        hPhiMetComp -> Fill(daughternu -> Phi, metvec.Phi(), weight);
      }
      // fill b parton histograms
      hPtLB -> Fill(b1_parton.Pt(), weight);
      hPtSLB -> Fill(b2_parton.Pt(), weight);
      hEtaLB -> Fill(b1_parton.Eta(), weight);
      hEtaSLB -> Fill(b2_parton.Eta(), weight);
      hPhiLB -> Fill(b1_parton.Phi(), weight);
      hPhiSLB -> Fill(b2_parton.Phi(), weight);
      if (fill_reco) {
        hPtLBComp -> Fill(b1_parton.Pt(), b1_reco.Pt(), weight);
        hPtSLBComp -> Fill(b2_parton.Pt(), b2_reco.Pt(), weight);
        hEtaLBComp -> Fill(b1_parton.Eta(), b1_reco.Eta(), weight);
        hEtaSLBComp -> Fill(b2_parton.Eta(), b2_reco.Eta(), weight);
        hPhiLBComp -> Fill(b1_parton.Phi(), b1_reco.Phi(), weight);
        hPhiSLBComp -> Fill(b2_parton.Phi(), b2_reco.Phi(), weight);
      }
      // fill parton level composite histograms
      hMTWH -> Fill((wpartonvec+hpartonvec).Mt(), weight);
      hPtWH -> Fill((wpartonvec+hpartonvec).Pt(), weight);
      hEtaWH -> Fill((wpartonvec+hpartonvec).Eta(), weight);
      hPhiWH -> Fill((wpartonvec+hpartonvec).Phi(), weight);
      if (fill_reco) {
        hPtWHComp -> Fill((wpartonvec+hpartonvec).Pt(), sysvec.Pt(), weight);
        hEtaWHComp -> Fill((wpartonvec+hpartonvec).Eta(), sysvec.Eta(), weight);
        hPhiWHComp -> Fill((wpartonvec+hpartonvec).Phi(), sysvec.Phi(), weight);
        hMTWHComp -> Fill((wpartonvec+hpartonvec).Mt(), sysvec.Mt(), weight);
      }
    }
    // gather particle level jets
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
    // fill particle level histograms
    if (fill_particle) {
      // fill particle level electron histograms
      if (eparticlevec.Pt() != 0) {
        hPtEP -> Fill(eparticlevec.Pt(), weight);
        hEtaEP -> Fill(eparticlevec.Eta(), weight);
        hPhiEP -> Fill(eparticlevec.Phi(), weight);
        if (fill_reco && elec != NULL) {
          hPtECompP -> Fill(eparticlevec.Pt(), elecvec.Pt(), weight);
          hEtaECompP -> Fill(eparticlevec.Eta(), elecvec.Eta(), weight);
        }
      // fill particle level muon histograms
      } else {
        hPtMuP -> Fill(mparticlevec.Pt(), weight);
        hEtaMuP -> Fill(mparticlevec.Eta(), weight);
        hPhiMuP -> Fill(mparticlevec.Phi(), weight);
        if (fill_reco && muon != NULL) {
          hPtMuCompP -> Fill(mparticlevec.Pt(), muvec.Pt(), weight);
          hEtaMuCompP -> Fill(mparticlevec.Eta(), muvec.Eta(), weight);
        }
      }
      // fill particle level met histograms
      hPtMetP -> Fill(genmetvec.Pt(), weight);
      hEtaMetP -> Fill(genmetvec.Eta(), weight);
      hPhiMetP -> Fill(genmetvec.Phi(), weight);
      if (fill_reco) {
        hPtMetCompP -> Fill(genmetvec.Pt(), metvec.Pt(), weight);
        hEtaMetCompP -> Fill(genmetvec.Eta(), metvec.Eta(), weight);
      }
      // fill particle level b jet histograms
      hPtLBP -> Fill(b1_particle.Pt(), weight);
      hPtSLBP -> Fill(b2_particle.Pt(), weight);
      hEtaLBP -> Fill(b1_particle.Eta(), weight);
      hEtaSLBP -> Fill(b2_particle.Eta(), weight);
      hPhiLBP -> Fill(b1_particle.Phi(), weight);
      hPhiSLBP -> Fill(b2_particle.Phi(), weight);
      if (fill_reco) {
        hPtLBCompP -> Fill(b1_particle.Pt(), b1_reco.Pt(), weight);
        hPtSLBCompP -> Fill(b2_particle.Pt(), b2_reco.Pt(), weight);
        hEtaLBCompP -> Fill(b1_particle.Eta(), b1_reco.Eta(), weight);
        hEtaSLBCompP -> Fill(b2_particle.Eta(), b2_reco.Eta(), weight);
      }
      // fill particle level higgs and w histograms
      hparticlevec = b1_particle + b2_particle;
      hMbbP -> Fill(hparticlevec.M(), weight);
      hMTlvP -> Fill(wparticlevec.Mt(), weight);
      hPtBBP -> Fill(hparticlevec.Pt(), weight);
      hPtLVP -> Fill(wparticlevec.Pt(), weight);
      hEtaBBP -> Fill(hparticlevec.Eta(), weight);
      hEtaLVP -> Fill(wparticlevec.Eta(), weight);
      hPhiBBP -> Fill(hparticlevec.Phi(), weight);
      hPhiLVP -> Fill(wparticlevec.Phi(), weight);
      // fill higgs and w particle-reco comparison histograms
      if (fill_reco) {
        hPtHCompP -> Fill(hparticlevec.Pt(), higgsvec.Pt(), weight);
        hEtaHCompP -> Fill(hparticlevec.Eta(), higgsvec.Eta(), weight);
        hPtWCompP -> Fill(wparticlevec.Pt(), wvec.Pt(), weight);
        hEtaWCompP -> Fill(wparticlevec.Eta(), wvec.Eta(), weight);
        hMTWCompP -> Fill(wparticlevec.Mt(), sysvec.Mt(), weight);
        hMbbComp -> Fill(hparticlevec.M(), higgsvec.M(), weight);
      }
      // fill particle level total event histograms
      hMTbblvP -> Fill((hparticlevec+wparticlevec).Mt(), weight);
      hPtBBLVP -> Fill((hparticlevec+wparticlevec).Pt(), weight);
      hEtaBBLVP -> Fill((hparticlevec+wparticlevec).Eta(), weight);
      hPhiBBLVP -> Fill((hparticlevec+wparticlevec).Phi(), weight);
      if (fill_reco) {
        hPtWHCompP -> Fill((wparticlevec+hparticlevec).Pt(), sysvec.Pt(), weight);
        hEtaWHCompP -> Fill((wparticlevec+hparticlevec).Eta(), sysvec.Eta(), weight);
        hMTWHCompP -> Fill((wparticlevec+hparticlevec).Mt(), sysvec.Mt(), weight);
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
  // Write histograms
  TFile *hists= new TFile("whbb_hists.root","recreate");
  hists->cd();
  //mass
  hMTW -> Write();
  hMTWH -> Write();
  hMbbR->Write();
  hMTlvR->Write();
  hMTbblvR -> Write();
  hMbbP->Write();
  hMTlvP->Write();
  hMTbblvP -> Write();
  hMbbComp -> Write();
  hMTWComp -> Write();
  hMTWHComp -> Write();
  hMTWCompP -> Write();
  hMTWHCompP -> Write();
  // pt
  hPtH->Write();
  hPtW->Write();
  hPtE->Write();
  hPtNu->Write();
  hPtMu->Write();
  hPtLB->Write();
  hPtSLB->Write();
  hPtWH->Write();
  hPtBBR->Write();
  hPtLVR->Write();
  hPtER->Write();
  hPtMetR->Write();
  hPtMuR->Write();
  hPtLBR->Write();
  hPtSLBR->Write();
  hPtBBLVR -> Write();
  hPtBBP->Write();
  hPtLVP->Write();
  hPtEP->Write();
  hPtMetP->Write();
  hPtMuP->Write();
  hPtLBP->Write();
  hPtSLBP->Write();
  hPtBBLVP->Write();
  hPtHComp->Write();
  hPtWComp->Write();
  hPtEComp->Write();
  hPtMetComp->Write();
  hPtMuComp->Write();
  hPtLBComp->Write();
  hPtSLBComp->Write();
  hPtWHComp->Write();
  hPtHCompP->Write();
  hPtWCompP->Write();
  hPtECompP->Write();
  hPtMetCompP->Write();
  hPtMuCompP->Write();
  hPtLBCompP->Write();
  hPtSLBCompP->Write();
  hPtWHCompP->Write();
  //eta
  hEtaH->Write();
  hEtaW->Write();
  hEtaE->Write();
  hEtaNu->Write();
  hEtaMu->Write();
  hEtaLB->Write();
  hEtaSLB->Write();
  hEtaWH->Write();
  hEtaBBR->Write();
  hEtaLVR->Write();
  hEtaER->Write();
  hEtaMetR->Write();
  hEtaMuR->Write();
  hEtaLBR->Write();
  hEtaSLBR->Write();
  hEtaBBLVR->Write();
  hEtaBBP->Write();
  hEtaLVP->Write();
  hEtaEP->Write();
  hEtaMetP->Write();
  hEtaMuP->Write();
  hEtaLBP->Write();
  hEtaSLBP->Write();
  hEtaBBLVP->Write();
  hEtaHComp->Write();
  hEtaWComp->Write();
  hEtaEComp->Write();
  hEtaMetComp->Write();
  hEtaMuComp->Write();
  hEtaLBComp->Write();
  hEtaSLBComp->Write();
  hEtaWHComp->Write();
  hEtaHCompP->Write();
  hEtaWCompP->Write();
  hEtaECompP->Write();
  hEtaMetCompP->Write();
  hEtaMuCompP->Write();
  hEtaLBCompP->Write();
  hEtaSLBCompP->Write();
  hEtaWHCompP->Write();
  // phi
  hPhiH->Write();
  hPhiW->Write();
  hPhiE->Write();
  hPhiNu->Write();
  hPhiMu->Write();
  hPhiLB->Write();
  hPhiSLB->Write();
  hPhiWH->Write();
  hPhiBBR->Write();
  hPhiLVR->Write();
  hPhiER->Write();
  hPhiMetR->Write();
  hPhiMuR->Write();
  hPhiLBR->Write();
  hPhiSLBR->Write();
  hPhiBBLVR->Write();
  hPhiBBP->Write();
  hPhiLVP->Write();
  hPhiEP->Write();
  hPhiMetP->Write();
  hPhiMuP->Write();
  hPhiLBP->Write();
  hPhiSLBP->Write();
  hPhiBBLVP->Write();
  hPhiHComp->Write();
  hPhiWComp->Write();
  hPhiEComp->Write();
  hPhiMetComp->Write();
  hPhiMuComp->Write();
  hPhiLBComp->Write();
  hPhiSLBComp->Write();
  hPhiWHComp->Write();
  hists->Close();
  // draw histograms
  gROOT->SetBatch(kTRUE);
  // mass
  draw_hist(hMTW, "massT_W", "m_{T}^{W}", "mass (GeV)");
  draw_hist(hMTWH, "massT_WH", "m_{T}^{WH}", "mass (GeV)");
  draw_hist(hMbbR, "mass_bb_reco", "Reco m_{bb}", "mass (GeV)");
  draw_hist(hMTlvR, "massT_lv_reco", "Reco m_{T}^{l#nu}", "mass (GeV)");
  draw_hist(hMTbblvR, "massT_bblv_reco", "Reco m_{T}^{bbl#nu}", "mass (GeV)");
  draw_hist(hMbbP, "mass_bb_particle", "m_{bb}", "mass (GeV)");
  draw_hist(hMTlvP, "massT_lv_particle", "m_{T}^{l#nu}", "mass (GeV)");
  draw_hist(hMTbblvP, "massT_bblv_particle", "m_{T}^{bbl#nu}", "mass (GeV)");
  draw_hist2(hMbbComp, "mass_bb_Comp", "m_{bb}", "Particle Level mass (GeV)", "Reco Level mass (GeV)");
  draw_hist2(hMTWComp, "massT_W_Comp", "m_{T}^{W}", "Parton Level mass (GeV)", "Reco Level mass (GeV)");
  draw_hist2(hMTWHComp, "massT_WH_Comp", "m_{T}^{WH}", "Parton Level mass (GeV)", "Reco Level mass (GeV)");
  draw_hist2(hMTWCompP, "massT_W_Comp_particle", "m_{T}^{W}", "Particle Level mass (GeV)", "Reco Level mass (GeV)");
  draw_hist2(hMTWHCompP, "massT_WH_Comp_particle", "m_{T}^{WH}", "Particle Level mass (GeV)", "Reco Level mass (GeV)");
  // pt
  draw_hist(hPtH, "pt_H", "p_{T}^{H}", "p_{T}");
  draw_hist(hPtW, "pt_W", "p_{T}^{W}", "p_{T}");
  draw_hist(hPtE, "pt_e",  "p_{T}^{e}", "p_{T}");
  draw_hist(hPtNu, "pt_nu",  "p_{T}^{MET}", "p_{T}");
  draw_hist(hPtMu, "pt_mu",  "p_{T}^{#mu}", "p_{T}");
  draw_hist(hPtLB, "pt_LB", "Lead b Parton p_{T}", "p_{T}");
  draw_hist(hPtSLB, "pt_SLB", "Sublead b Parton p_{T}", "p_{T}");
  draw_hist(hPtWH, "pt_WH", "p_{T}^{WH}", "p_{T}");
  draw_hist(hPtBBR, "pt_bb_reco", "Reco p_{T}^{bb}", "p_{T}");
  draw_hist(hPtLVR, "pt_lv_reco", "Reco p_{T}^{l#nu}", "p_{T}");
  draw_hist(hPtER, "pt_e_reco",  "Reco p_{T}^{e}", "p_{T}");
  draw_hist(hPtMetR, "pt_met_reco",  "Reco p_{T}^{MET}", "p_{T}");
  draw_hist(hPtMuR, "pt_mu_reco",  "Reco p_{T}^{#mu}", "p_{T}");
  draw_hist(hPtLBR, "pt_LB_reco", "Lead Reco b Jet p_{T}", "p_{T}");
  draw_hist(hPtSLBR, "pt_SLB_reco", "Sublead Reco b Jet p_{T}", "p_{T}");
  draw_hist(hPtBBLVR, "pt_bblv_reco", "Reco p_{T}^{bbl#nu}", "p_{T}");
  draw_hist(hPtBBP, "pt_bb_particle", "p_{T}^{bb}", "p_{T}");
  draw_hist(hPtLVP, "pt_lv_particle", "p_{T}^{l#nu}", "p_{T}");
  draw_hist(hPtEP, "pt_e_particle",  "p_{T}^{e}", "p_{T}");
  draw_hist(hPtMetP, "pt_met_particle", "p_{T}^{MET}", "p_{T}");
  draw_hist(hPtMuP, "pt_mu_particle",  "p_{T}^{#mu}", "p_{T}");
  draw_hist(hPtLBP, "pt_LB_particle", "Lead b Jet p_{T}", "p_{T}");
  draw_hist(hPtSLBP, "pt_SLB_particle", "Sublead b Jet p_{T}", "p_{T}");
  draw_hist(hPtBBLVP, "pt_bblv_particle", "p_{T}^{bbl#nu}", "p_{T}");
  // pt comparison
  draw_hist2(hPtHComp, "pt_H_Comp", "p_{T}^{H}", "Parton Level P_{T} (GeV)", "Reco Level P_{T} (GeV)");
  draw_hist2(hPtWComp, "pt_W_Comp", "p_{T}^{W}", "Parton Level P_{T} (GeV)", "Reco Level P_{T} (GeV)");
  draw_hist2(hPtEComp, "pt_e_Comp", "p_{T}^{e}", "Parton Level P_{T} (GeV)", "Reco Level P_{T} (GeV)");
  draw_hist2(hPtMetComp, "pt_nu_Comp", "p_{T}^{#nu}", "Parton Level P_{T} (GeV)", "Reco Level P_{T} (GeV)");
  draw_hist2(hPtMuComp, "pt_mu_Comp", "p_{T}^{#mu}", "Parton Level P_{T} (GeV)", "Reco Level P_{T} (GeV)");
  draw_hist2(hPtLBComp, "pt_LB_Comp", "Lead p_{T}^{b}", "Parton Level P_{T} (GeV)", "Reco Level P_{T} (GeV)");
  draw_hist2(hPtSLBComp, "pt_SLB_Comp", "Sublead p_{T}^{b}", "Parton Level P_{T} (GeV)", "Reco Level P_{T} (GeV)");
  draw_hist2(hPtWHComp, "pt_WH_Comp", "p_{T}^{WH}", "Parton Level P_{T} (GeV)", "Reco Level P_{T} (GeV)");

  draw_hist2(hPtHCompP, "pt_H_Comp_particle", "p_{T}^{H}", "Particle Level P_{T} (GeV)", "Reco Level P_{T} (GeV)");
  draw_hist2(hPtWCompP, "pt_W_Comp_particle", "p_{T}^{W}", "Particle Level P_{T} (GeV)", "Reco Level P_{T} (GeV)");
  draw_hist2(hPtECompP, "pt_e_Comp_particle", "p_{T}^{e}", "Particle Level P_{T} (GeV)", "Reco Level P_{T} (GeV)");
  draw_hist2(hPtMetCompP, "pt_nu_Comp_particle", "p_{T}^{#nu}", "Particle Level P_{T} (GeV)", "Reco Level P_{T} (GeV)");
  draw_hist2(hPtMuCompP, "pt_mu_Comp_particle", "p_{T}^{#mu}", "Particle Level P_{T} (GeV)", "Reco Level P_{T} (GeV)");
  draw_hist2(hPtLBCompP, "pt_LB_Comp_particle", "Lead p_{T}^{b}", "Particle Level P_{T} (GeV)", "Reco Level P_{T} (GeV)");
  draw_hist2(hPtSLBCompP, "pt_SLB_Comp_particle", "Sublead p_{T}^{b}", "Particle Level P_{T} (GeV)", "Reco Level P_{T} (GeV)");
  draw_hist2(hPtWHCompP, "pt_WH_Comp_particle", "p_{T}^{WH}", "Particle Level P_{T} (GeV)", "Reco Level P_{T} (GeV)");
  // eta
  draw_hist(hEtaH, "eta_H", "#eta_{H}", "#eta");
  draw_hist(hEtaW, "eta_W", "#eta_{W}", "#eta");
  draw_hist(hEtaE, "eta_e",  "#eta_{e}", "#eta");
  draw_hist(hEtaNu, "eta_nu",  "#eta_{MET}", "#eta");
  draw_hist(hEtaMu, "eta_mu",  "#eta_{#mu}", "#eta");
  draw_hist(hEtaLB, "eta_LB", "Lead #eta_{b}", "#eta");
  draw_hist(hEtaSLB, "eta_SLB", "Sublead #eta_{b}", "#eta");
  draw_hist(hEtaWH, "eta_WH", "#eta_{WH}", "#eta");
  draw_hist(hEtaBBR, "eta_bb_reco", "Reco #eta_{bb}", "#eta");
  draw_hist(hEtaLVR, "eta_lv_reco", "Reco #eta_{l#nu}", "#eta");
  draw_hist(hEtaER, "eta_e_reco",  "Reco #eta_{e}", "#eta");
  draw_hist(hEtaMetR, "eta_met_reco",  "Reco #eta_{MET}", "#eta");
  draw_hist(hEtaMuR, "eta_mu_reco",  "Reco #eta_{#mu}", "#eta");
  draw_hist(hEtaLBR, "eta_LB_reco", "Lead Reco b Jet #eta", "#eta");
  draw_hist(hEtaSLBR, "eta_SLB_reco", "Sublead Reco b Jet #eta", "#eta");
  draw_hist(hEtaBBLVR, "eta_bblv_reco", "Reco #eta_{bbl#nu}", "#eta");
  draw_hist(hEtaBBP, "eta_bb_particle", "#eta_{bb}", "#eta");
  draw_hist(hEtaLVP, "eta_lv_particle", "#eta_{l#nu}", "#eta");
  draw_hist(hEtaEP, "eta_e_particle",  "#eta_{e}", "#eta");
  draw_hist(hEtaMetP, "eta_met_particle",  "#eta_{MET}", "#eta");
  draw_hist(hEtaMuP, "eta_mu_particle",  "#eta_{#mu}", "#eta");
  draw_hist(hEtaLBP, "eta_LB_particle", "Lead b Jet #eta", "#eta");
  draw_hist(hEtaSLBP, "eta_SLB_particle", "Sublead b Jet #eta", "#eta");
  draw_hist(hEtaBBLVP, "eta_bblv_particle", "#eta_{bbl#nu}", "#eta");
  // eta comparison
  draw_hist2(hEtaHComp, "eta_H_Comp", "#eta_{H}", "Parton Level #eta", "Reco Level #eta");
  draw_hist2(hEtaWComp, "eta_W_Comp", "#eta_{W}", "Parton Level #eta", "Reco Level #eta");
  draw_hist2(hEtaEComp, "eta_e_Comp", "#eta_{e}","Parton Level #eta", "Reco Level #eta");
  draw_hist2(hEtaMetComp, "eta_nu_Comp", "#eta_{#nu}","Parton Level #eta", "Reco Level #eta");
  draw_hist2(hEtaMuComp, "eta_mu_Comp", "#eta_{#mu}","Parton Level #eta", "Reco Level #eta");
  draw_hist2(hEtaLBComp, "eta_LB_Comp", "Lead #eta_{b}","Parton Level #eta", "Reco Level #eta");
  draw_hist2(hEtaSLBComp, "eta_SLB_Comp", "Sublead #eta_{b}","Parton Level #eta", "Reco Level #eta");
  draw_hist2(hEtaWHComp, "eta_WH_Comp", "Sublead #eta_{WH}", "Parton Level #eta", "Reco Level #eta");

  draw_hist2(hEtaHCompP, "eta_H_Comp_particle", "#eta_{H}", "Particle Level #eta", "Reco Level #eta");
  draw_hist2(hEtaWCompP, "eta_W_Comp_particle", "#eta_{W}", "Particle Level #eta", "Reco Level #eta");
  draw_hist2(hEtaECompP, "eta_e_Comp_particle", "#eta_{e}","Particle Level #eta", "Reco Level #eta");
  draw_hist2(hEtaMetCompP, "eta_nu_Comp_particle", "#eta_{#nu}","Particle Level #eta", "Reco Level #eta");
  draw_hist2(hEtaMuCompP, "eta_mu_Comp_particle", "#eta_{#mu}","Particle Level #eta", "Reco Level #eta");
  draw_hist2(hEtaLBCompP, "eta_LB_Comp_particle", "Lead #eta_{b}","Particle Level #eta", "Reco Level #eta");
  draw_hist2(hEtaSLBCompP, "eta_SLB_Comp_particle", "Sublead #eta_{b}","Particle Level #eta", "Reco Level #eta");
  draw_hist2(hEtaWHCompP, "eta_WH_Comp_particle", "Sublead #eta_{WH}", "Particle Level #eta", "Reco Level #eta");
  // phi parton level
  draw_hist(hPhiH, "phi_H", "#phi_{H}", "#phi (Rad)");
  draw_hist(hPhiW, "phi_W", "#phi_{W}", "#phi (Rad)");
  draw_hist(hPhiE, "phi_e",  "#phi_{e}", "#phi (Rad)");
  draw_hist(hPhiNu, "phi_nu",  "#phi_{MET}", "#phi (Rad)");
  draw_hist(hPhiMu, "phi_mu",  "#phi_{#mu}", "#phi (Rad)");
  draw_hist(hPhiLB, "phi_LB", "Lead #phi_{b}", "#phi (Rad)");
  draw_hist(hPhiSLB, "phi_SLB", "Sublead #phi_{b}", "#phi (Rad)");
  draw_hist(hPhiWH, "phi_WH", "#phi_{WH}", "#phi (Rad)");
  draw_hist(hPhiBBR, "phi_bb_reco", "Reco #phi_{bb}", "#phi (Rad)");
  draw_hist(hPhiLVR, "phi_lv_reco", "Reco #phi_{l#nu}", "#phi (Rad)");
  draw_hist(hPhiER, "phi_e_reco",  "Reco #phi_{e}", "#phi (Rad)");
  draw_hist(hPhiMetR, "phi_met_reco",  "Reco #phi_{MET}", "#phi (Rad)");
  draw_hist(hPhiMuR, "phi_mu_reco",  "Reco #phi_{#mu}", "#phi (Rad)");
  draw_hist(hPhiLBR, "phi_LB_reco", "Lead Reco b Jet #phi", "#phi (Rad)");
  draw_hist(hPhiSLBR, "phi_SLB_reco", "Sublead Reco b Jet #phi", "#phi (Rad)");
  draw_hist(hPhiBBLVR, "phi_bblv_reco", "Reco #phi_{bbl#nu}", "#phi (Rad)");
  draw_hist(hPhiBBP, "phi_bb_particle", "#phi_{bb}", "#phi (Rad)");
  draw_hist(hPhiLVP, "phi_lv_particle", "#phi_{l#nu}", "#phi (Rad)");
  draw_hist(hPhiEP, "phi_e_particle",  "#phi_{e}", "#phi (Rad)");
  draw_hist(hPhiMetP, "phi_met_particle",  "#phi_{MET}", "#phi (Rad)");
  draw_hist(hPhiMuP, "phi_mu_particle",  "#phi_{#mu}", "#phi (Rad)");
  draw_hist(hPhiLBP, "phi_LB_particle", "Lead b Jet #phi", "#phi (Rad)");
  draw_hist(hPhiSLBP, "phi_SLB_particle", "Sublead b Jet #phi", "#phi (Rad)");
  draw_hist(hPhiBBLVP, "phi_bblv_particle", "#phi_{bbl#nu}", "#phi (Rad)");
  // phi comparison
  draw_hist2(hPhiHComp, "phi_H_Comp", "#phi_{H}", "Parton Level #phi (Rad)", "Reco Level #phi (Rad)");
  draw_hist2(hPhiWComp, "phi_W_Comp", "#phi_{W}", "Parton Level #phi (Rad)", "Reco Level #phi (Rad)");
  draw_hist2(hPhiEComp, "phi_e_Comp", "#phi_{e}","Parton Level #phi (Rad)", "Reco Level #phi (Rad)");
  draw_hist2(hPhiMetComp, "phi_nu_Comp", "#phi_{#nu}","Parton Level #phi (Rad)", "Reco Level #phi (Rad)");
  draw_hist2(hPhiMuComp, "phi_mu_Comp", "#phi_{#mu}","Parton Level #phi (Rad)", "Reco Level #phi (Rad)");
  draw_hist2(hPhiLBComp, "phi_LB_Comp", "Lead #phi_{b}","Parton Level #phi (Rad)", "Reco Level #phi (Rad)");
  draw_hist2(hPhiSLBComp, "phi_SLB_Comp", "Sublead #phi_{b}","Parton Level #phi (Rad)", "Reco Level #phi (Rad)");
  draw_hist2(hPhiWHComp, "phi_WH_Comp", "Sublead #phi_{WH}", "Parton Level #phi (Rad)", "Reco Level #phi (Rad)");
  // print some info
  std::cout<<"Integral of m_bb hist "<< hMbbR->Integral()<<" integral of mT_lv "<<hMTlvR->Integral()<<std::endl;
  std::cout<<"Total number of entries "<<numberOfEntries<<" Passed "<<nPassed<<" raw "<<nPassedRaw<<std::endl;
}
