// root -l -b
// .x examples/whbb_analyze.C("../wplus_hbb/Events/run_01/tag_1_delphes_events.root", "../wminus_hbb/Events/run_01/tag_1_delphes_events.root")
#ifdef __CLING__
R__LOAD_LIBRARY(libDelphes)
#include "classes/DelphesClasses.h"
#include "external/ExRootAnalysis/ExRootTreeReader.h"
#include "DSCBf.h"
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

GenParticle* find_status1_child(TClonesArray *branchGenParticle, GenParticle *particle, int target) {
  if (abs(particle -> PID) == target && particle -> Status == 1) {
    return particle;
  }
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
  if (abs(d1_pid) == target) {
    return find_status1_child(branchGenParticle, daughter1, target);
  } else if (abs(d2_pid) == target) {
    return find_status1_child(branchGenParticle, daughter2, target);
  } else {
    return NULL;
  }
}


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
  // Book histograms
  // mass
  TH1 *hMTW = new TH1F("massT_W", "m_{T}^{W}", 100, 40.0, 140.0);
  TH1 *hMTWH = new TH1F("massT_WH", "m_{T}^{WH}", 100, 120.0, 670.0);
  TH1 *hMbbR = new TH1F("mass_bb_reco", "Reco m_{bb}", 100, 0.0, 250.0);
  TH1 *hMTlvR = new TH1F("massT_lv_reco", "Reco m_{T}^{l#nu}", 100, 40.0, 140.0);
  TH1 *hMTbblvR = new TH1F("massT_bblv_reco", "Reco m_{T}^{bbl#nu}", 100, 120.0, 670.0);
  TH1 *hMbbP = new TH1F("mass_bb_particle", "m_{bb}", 100, 0.0, 250.0);
  TH1 *hMTlvP = new TH1F("massT_lv_particle", "m_{T}^{l#nu}", 100, 40.0, 140.0);
  TH1 *hMTbblvP = new TH1F("massT_bblv_particle", "m_{T}^{bbl#nu}", 100, 120.0, 670.0);
  TH2 *hMbbComp = new TH2F("mass_bb_Comp", "m_{bb}", 50, 0.0, 250.0, 50, 0.0, 250.0);
  TH2 *hMTWComp = new TH2F("massT_W_Comp", "m_{T}^{W}", 50, 40.0, 140.0, 50, 40.0, 140.0);
  TH2 *hMTWHComp = new TH2F("massT_WH_Comp", "m_{T}^{WH}", 50, 120.0, 670.0, 50, 120.0, 670.0);
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
  // pt comparison
  TH2 *hPtHComp = new TH2F("pt_H_Comp", "p_{T}^{H}", 50, 0.0, 400.0, 50, 0.0, 400.0);
  TH2 *hPtWComp = new TH2F("pt_W_Comp", "p_{T}^{W}", 50, 0.0, 400.0, 50, 0.0, 400.0);
  TH2 *hPtEComp = new TH2F("pt_e_Comp", "p_{T}^{e}", 50, 20.0, 200.0, 50, 20.0, 200.0);
  TH2 *hPtMetComp = new TH2F("pt_nu_Comp", "p_{T}^{#nu}", 50, 20.0, 200.0, 50, 20.0, 200.0);
  TH2 *hPtMuComp = new TH2F("pt_mu_Comp", "p_{T}^{#mu}", 50, 20.0, 200.0, 50, 20.0, 200.0);
  TH2 *hPtLBComp = new TH2F("pt_LB_Comp", "Lead p_{T}^{b}", 50, 20.0, 200.0, 50, 20.0, 200.0);
  TH2 *hPtSLBComp = new TH2F("pt_SLB_Comp", "Sublead p_{T}^{b}", 50, 20.0, 200.0, 50, 20.0, 200.0);
  TH2 *hPtWHComp = new TH2F("pt_WH_Comp", "Sublead p_{T}^{WH}", 50, 0.0, 400.0, 50, 0.0, 400.0);
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
  // eta comparison
  TH2 *hEtaHComp = new TH2F("eta_H_Comp", "#eta_{H}", 50, -3.0, 3.0, 50, -3.0, 3.0);
  TH2 *hEtaWComp = new TH2F("eta_W_Comp", "#eta_{W}", 50, -3.0, 3.0, 50, -3.0, 3.0);
  TH2 *hEtaEComp = new TH2F("eta_e_Comp", "#eta_{e}", 50,-3.0, 3.0, 50, -3.0, 3.0);
  TH2 *hEtaMetComp = new TH2F("eta_nu_Comp", "#eta_{#nu}", 50,-3.0, 3.0, 50, -3.0, 3.0);
  TH2 *hEtaMuComp = new TH2F("eta_mu_Comp", "#eta_{#mu}", 50,-3.0, 3.0, 50, -3.0, 3.0);
  TH2 *hEtaLBComp = new TH2F("eta_LB_Comp", "Lead #eta_{b}", 50,-3.0, 3.0, 50, -3.0, 3.0);
  TH2 *hEtaSLBComp = new TH2F("eta_SLB_Comp", "Sublead #eta_{b}", 50,-3.0, 3.0, 50, -3.0, 3.0);
  TH2 *hEtaWHComp = new TH2F("eta_WH_Comp", "Sublead #eta_{WH}", 50, -3.0, 3.0, 50, -3.0, 3.0);
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
  Muon *muon;
  MissingET *met;
  MissingET *genmet;
  GenParticle *daughter1;
  GenParticle *daughter2;
  GenParticle *daughter;
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

  // Loop over events
  for(Int_t entry = 0; entry < numberOfEntries; ++entry){
    // Load selected branches with data from specified event
    treeReader->ReadEntry(entry);
    HepMCEvent *event = (HepMCEvent*) branchEvent -> At(0);
    Float_t weight = event->Weight/numberOfEntries*Lumi;
    // Loop over jets in event and save the b jets in 2b events as lorentz vectors to reconstruct the higgs
    bjets = 0;
    for(int i=0; i<(int)branchJet->GetEntries(); i++){
      Jet *jet=(Jet*) branchJet->At(i);
      if (jet -> Flavor == 5) bjets += 1; // add additional cuts?
      else continue;
      if (bjets == 1) {
        b1_reco = jet->P4();
      } else if (bjets == 2){
        b2_reco = jet->P4();
        break;
      }
    }
    if (bjets != 2) continue;
    higgsvec = b1_reco + b2_reco;
    // Now check for events with 1 charged lepton and save its information
    // Also combine it with MET data to reconstruct the W
    met = (MissingET *) branchMissingET->At(0);
    metvec = met ->P4();
    if (branchElectron->GetEntries() == 1 && branchMuon->GetEntries() == 0){
      elec = (Electron *) branchElectron->At(0);
      if ((elec->Charge == 1 && entry < numberOfPlusEntries) || (elec->Charge == -1 && entry >= numberOfPlusEntries)) {
        elecvec = elec->P4();
        wvec = elecvec + metvec;
        hPtER -> Fill(elec->PT, weight);
        hEtaER -> Fill(elec->Eta, weight);
        hPhiER -> Fill(elec->Phi, weight);
      }
    } else if (branchMuon->GetEntries() == 1 && branchElectron->GetEntries() == 0) {
      muon = (Muon *) branchMuon->At(0);
      if ((muon->Charge == 1 && entry < numberOfPlusEntries) || (muon->Charge == -1 && entry >= numberOfPlusEntries)) {
        muvec = muon->P4();
        wvec = muvec + metvec;
        hPtMuR -> Fill(muon->PT, weight);
        hEtaMuR -> Fill(muon->Eta, weight);
        hPhiMuR -> Fill(muon->Phi, weight);
      }
    } else {
      continue;
    }
    sysvec = higgsvec + wvec;
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
    // loop over true particles and fill those histograms
    b1_parton.SetPtEtaPhiM(0,0,0,0);
    b2_parton.SetPtEtaPhiM(0,0,0,0);
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
      if (particle->PID == 25 && d1_pid != 25) {
        hPtH -> Fill(particle -> PT, weight);
        hEtaH -> Fill(particle -> Eta, weight);
        hPhiH -> Fill(particle -> Phi, weight);
        hPtHComp -> Fill(particle -> PT, higgsvec.Pt(), weight);
        hEtaHComp -> Fill(particle -> Eta, higgsvec.Eta(), weight);
        hPhiHComp -> Fill(particle -> Phi, higgsvec.Phi(), weight);
        hpartonvec = particle->P4();
        // b partons
        if (abs(d1_pid) == 5) {
          if (daughter1 -> PT > daughter2 -> PT) {
            b1_parton = daughter1 -> P4();
            b2_parton = daughter2 -> P4();
          } else {
            b1_parton = daughter2 -> P4();
            b2_parton = daughter1 -> P4();
          }
        }
      // W parton
      } else if ((particle -> PID == 24 && d1_pid != 24) || (particle -> PID == -24 && d1_pid != -24)) {
        hPtW -> Fill(particle -> PT, weight);
        hEtaW -> Fill(particle -> Eta, weight);
        hPhiW -> Fill(particle -> Phi, weight);
        hPtWComp -> Fill(particle -> PT, wvec.Pt(), weight);
        hEtaWComp -> Fill(particle -> Eta, wvec.Eta(), weight);
        hPhiWComp -> Fill(particle -> Phi, wvec.Phi(), weight);
        wpartonvec = particle->P4();
        hMTW -> Fill(wpartonvec.Mt(), weight);
        hMTWComp -> Fill(wpartonvec.Mt(), sysvec.Mt(), weight);
        // check for electron daughter 1
        if (abs(d1_pid) == 11 || abs(d2_pid) == 11) {
          // parton level
          if (abs(d1_pid) == 11) {
            daughter = daughter1;
          } else {
            daughter = daughter2;
          }
          hPtE -> Fill(daughter -> PT, weight);
          hEtaE -> Fill(daughter -> Eta, weight);
          hPhiE -> Fill(daughter -> Phi, weight);
          hPtEComp -> Fill(daughter -> PT, elecvec.Pt(), weight);
          hEtaEComp -> Fill(daughter -> Eta, elecvec.Eta(), weight);
          hPhiEComp -> Fill(daughter -> Phi, elecvec.Phi(), weight);
          // particle level
          eparticlevec = find_status1_child(branchGenParticle, daughter, 11) -> P4();
          hPtEP -> Fill(eparticlevec.Pt(), weight);
          hEtaEP -> Fill(eparticlevec.Eta(), weight);
          hPhiEP -> Fill(eparticlevec.Phi(), weight);
          wparticlevec = genmetvec + eparticlevec;
        // check for muon daughter 1
        } else if (abs(d1_pid) == 13 || abs(d2_pid) == 13) {
          if (abs(d1_pid) == 13) {
            daughter = daughter1;
          } else {
            daughter = daughter2;
          }
          // parton level
          hPtMu -> Fill(daughter -> PT, weight);
          hEtaMu -> Fill(daughter -> Eta, weight);
          hPhiMu -> Fill(daughter -> Phi, weight);
          hPtMuComp -> Fill(daughter -> PT, muvec.Pt(), weight);
          hEtaMuComp -> Fill(daughter -> Eta, muvec.Eta(), weight);
          hPhiMuComp -> Fill(daughter -> Phi, muvec.Phi(), weight);
          // particle level
          mparticlevec = find_status1_child(branchGenParticle, daughter, 13) -> P4();
          hPtMuP -> Fill(mparticlevec.Pt(), weight);
          hEtaMuP -> Fill(mparticlevec.Eta(), weight);
          hPhiMuP -> Fill(mparticlevec.Phi(), weight);
          wparticlevec = genmetvec + mparticlevec;
        // check for neutrino daughter 1
        } else if (abs(d1_pid) == 12 || abs(d1_pid) == 14 || abs(d2_pid) == 12 || abs(d2_pid) == 14) {
          if (abs(d1_pid) == 12 || abs(d1_pid) == 14) {
            daughter = daughter1;
          } else {
            daughter = daughter2;
          }
          hPtNu -> Fill(daughter -> PT, weight);
          hEtaNu -> Fill(daughter -> Eta, weight);
          hPhiNu -> Fill(daughter -> Phi, weight);
          hPtMetComp -> Fill(daughter -> PT, metvec.Pt(), weight);
          hEtaMetComp -> Fill(daughter -> Eta, metvec.Eta(), weight);
          hPhiMetComp -> Fill(daughter -> Phi, metvec.Phi(), weight);
        }
      }
    }
    // fill b parton histograms
    hPtLB -> Fill(b1_parton.Pt(), weight);
    hPtSLB -> Fill(b2_parton.Pt(), weight);
    hEtaLB -> Fill(b1_parton.Eta(), weight);
    hEtaSLB -> Fill(b2_parton.Eta(), weight);
    hPhiLB -> Fill(b1_parton.Phi(), weight);
    hPhiSLB -> Fill(b2_parton.Phi(), weight);
    hPtLBComp -> Fill(b1_parton.Pt(), b1_reco.Pt(), weight);
    hPtSLBComp -> Fill(b2_parton.Pt(), b2_reco.Pt(), weight);
    hEtaLBComp -> Fill(b1_parton.Eta(), b1_reco.Eta(), weight);
    hEtaSLBComp -> Fill(b2_parton.Eta(), b2_reco.Eta(), weight);
    hPhiLBComp -> Fill(b1_parton.Phi(), b1_reco.Phi(), weight);
    hPhiSLBComp -> Fill(b2_parton.Phi(), b2_reco.Phi(), weight);
    // fill parton level composite histograms
    hMTWH -> Fill((wpartonvec+hpartonvec).Mt(), weight);
    hMTWHComp -> Fill((wpartonvec+hpartonvec).Mt(), sysvec.Mt(), weight);
    hPtWH -> Fill((wpartonvec+hpartonvec).Pt(), weight);
    hEtaWH -> Fill((wpartonvec+hpartonvec).Eta(), weight);
    hPhiWH -> Fill((wpartonvec+hpartonvec).Phi(), weight);
    hPtWHComp -> Fill((wpartonvec+hpartonvec).Pt(), sysvec.Pt(), weight);
    hEtaWHComp -> Fill((wpartonvec+hpartonvec).Eta(), sysvec.Eta(), weight);
    hPhiWHComp -> Fill((wpartonvec+hpartonvec).Phi(), sysvec.Phi(), weight);
    // fill particle level met histograms
    hPtMetP -> Fill(genmetvec.Pt(), weight);
    hEtaMetP -> Fill(genmetvec.Eta(), weight);
    hPhiMetP -> Fill(genmetvec.Phi(), weight);
    // gather particle level jets
    bjets = 0;
    for(int i=0; i<(int)branchGenJet->GetEntries(); i++){
      Jet *genjet=(Jet*) branchGenJet->At(i);
      if (genjet -> Flavor == 5) bjets += 1; // add additional cuts?
      else continue;
      if (bjets == 1) {
        b1_particle = genjet->P4();
      } else if (bjets == 2){
        b2_particle = genjet->P4();
        break;
      }
    }
    // fill particle level b jet histograms
    hPtLBP -> Fill(b1_particle.Pt(), weight);
    hPtSLBP -> Fill(b2_particle.Pt(), weight);
    hEtaLBP -> Fill(b1_particle.Eta(), weight);
    hEtaSLBP -> Fill(b2_particle.Eta(), weight);
    hPhiLBP -> Fill(b1_particle.Phi(), weight);
    hPhiSLBP -> Fill(b2_particle.Phi(), weight);
    // fill particle level composite histograms
    hparticlevec = b1_particle + b2_particle;
    hMbbP -> Fill(hparticlevec.M(), weight);
    hMTlvP -> Fill(wparticlevec.Mt(), weight);
    hMTbblvP -> Fill((hparticlevec+wparticlevec).Mt(), weight);
    hMbbComp -> Fill(hparticlevec.M(), higgsvec.M(), weight);
    hPtBBP -> Fill(hparticlevec.Pt(), weight);
    hPtLVP -> Fill(wparticlevec.Pt(), weight);
    hPtBBLVP -> Fill((hparticlevec+wparticlevec).Pt(), weight);
    hEtaBBP -> Fill(hparticlevec.Eta(), weight);
    hEtaLVP -> Fill(wparticlevec.Eta(), weight);
    hEtaBBLVP -> Fill((hparticlevec+wparticlevec).Eta(), weight);
    hPhiBBP -> Fill(hparticlevec.Phi(), weight);
    hPhiLVP -> Fill(wparticlevec.Phi(), weight);
    hPhiBBLVP -> Fill((hparticlevec+wparticlevec).Phi(), weight);
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
