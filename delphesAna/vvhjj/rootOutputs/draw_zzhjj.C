#ifdef __CLING__
R__LOAD_LIBRARY(libDelphes)
#include <iostream>
#include <fstream>
#include <cmath>
#endif

/*
// be able to calculate total number of events in a process
Long64_t get_file_num_events(const char *inputName) {
  TChain chain("Delphes");
  chain.Add(inputName);
  ExRootTreeReader *treeReader = new ExRootTreeReader(&chain);
  return treeReader->GetEntries();
}

Float_t get_total_events(const char *process_name) {
  std::string inputFileName = std::string(process_name) + "_inputs.txt";
  std::ifstream inputFile(inputFileName.c_str());
  std::string line;
  TChain chain("Delphes");
  Long64_t total = 0;
  while (std::getline(inputFile, line)) {
    total += get_file_num_events(line.c_str());
  }
  return static_cast<Float_t>(total);
}

// be able to calculate total weight of events in a process
Float_t get_file_weight(const char *inputName) {
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
}

// store the scales for the different processes as global variables
Float_t sig_scale = 0;
Float_t ttbar_scale = 0;
Float_t ttHbb_scale = 0;
Float_t drellyan_scale = 0;
Float_t diboson_scale = 0;
*/

// save plots as different file types
void PrintCanvas(TCanvas *c=nullptr, string name="default", string outputFolder="default", string subFolder="default"){
  std::vector <string> types={"png"}; 
  for(std::vector<string>::iterator it=types.begin(); it!=types.end(); it++) {
    c->Print(Form("%s/%s/%s.%s", outputFolder.c_str(), subFolder.c_str(), name.c_str(), (*it).c_str()), (*it).c_str());
    c->SetLogy(); 
    c->Print(Form("%s/%s/%s_log.%s", outputFolder.c_str(), subFolder.c_str(), name.c_str(), (*it).c_str()), (*it).c_str());
  }
}


// draw 1d histograms
void draw_stack(TFile *sig_file, TFile *ttbar_file, TFile *ttHbb_file, TFile *diboson_file, TFile *drellyan_file, string name, string title, string axistitle, string outputFolder) {

  // get the histograms from the files
  TH1F *sig_hist = (TH1F*)sig_file->Get(name.c_str());
  TH1F *ttbar_hist = (TH1F*)ttbar_file->Get(name.c_str());
  TH1F *ttHbb_hist = (TH1F*)ttHbb_file->Get(name.c_str());
  TH1F *diboson_hist = (TH1F*)diboson_file->Get(name.c_str());
  TH1F *drellyan_hist = (TH1F*)drellyan_file->Get(name.c_str());

  TH1F *sigClone=(TH1F*)sig_hist->Clone("sigClone");

  double lumiScaling = 10;
  
  // change lumi if desired
  sig_hist->Scale(lumiScaling*0.000013); //*0.000013 for W
  ttbar_hist->Scale(lumiScaling);
  diboson_hist->Scale(lumiScaling);
  ttHbb_hist->Scale(lumiScaling);
  drellyan_hist->Scale(lumiScaling);
  sigClone->Scale(lumiScaling*0.000013); //*0.000013 for W    

  // rebin if desired
  sig_hist->Rebin(1);
  ttbar_hist->Rebin(1);
  ttHbb_hist->Rebin(1);
  diboson_hist->Rebin(1);
  drellyan_hist->Rebin(1);

  // total hist 
  TH1F *total=(TH1F*)ttbar_hist->Clone("total");
  total->Add(sig_hist);
  total->Add(ttHbb_hist);
  total->Add(drellyan_hist); 
  total->Add(diboson_hist);

  // fill color
  sig_hist->SetFillColor(kCyan);
  ttbar_hist->SetFillColor(kOrange);
  ttHbb_hist->SetFillColor(kMagenta);
  diboson_hist->SetFillColor(kRed);
  drellyan_hist->SetFillColor(kYellow);

  /*
  sig_hist->Scale(sig_scale);
  ttbar_hist->Scale(ttbar_scale);
  ttHbb_hist->Scale(ttHbb_scale);
  drellyan_hist->Scale(drellyan_scale);
  diboson_hist->Scale(diboson_scale);
  */

  // form the histogram stack
  THStack *stack = new THStack(name.c_str(), title.c_str());
  stack->Add(ttbar_hist);
  stack->Add(drellyan_hist);
  stack->Add(diboson_hist);
  stack->Add(ttHbb_hist);
  stack->Add(sig_hist);
  
  // make a legend
  TLegend *legend = new TLegend(0.775, 0.775, 0.875, 0.875);
  legend->AddEntry(sig_hist, "Signal", "f");
  legend->AddEntry(drellyan_hist, "drell-yan", "f");
  legend->AddEntry(ttbar_hist, "ttbar", "f");
  legend->AddEntry(ttHbb_hist, "ttHbb", "f");
  legend->AddEntry(diboson_hist, "diboson", "f");
  
  // make a canvas and draw on it
  TCanvas *c = new TCanvas(name.c_str(), title.c_str(), 1500, 1200);
  c->cd();
  stack->Draw("hist");
  total->Draw("hist same pe");
  stack->GetXaxis()->SetTitle(axistitle.c_str());
  legend->Draw();

  // sig on top
  sigClone->Scale(100000);
  sigClone->SetLineColor(kBlack);
  //sigClone->SetFillColor(kWhite);
  //legend->AddEntry(sigClone, "Signal x 1000000", "l");
  sigClone->SetLineWidth(2);
  //  sigClone->Draw("hist same");
  cout <<" Output folder "<<outputFolder<<endl;
  
  TString nameOut(name);
  nameOut.ReplaceAll(" ","_");
  nameOut.ReplaceAll("#","");
  nameOut.ReplaceAll("^","");
  nameOut.ReplaceAll("/","");
  nameOut.ReplaceAll("{","");
  nameOut.ReplaceAll("}","");
  cout<<" Name "<<nameOut.Data()<<endl;
  PrintCanvas(c,nameOut.Data(), outputFolder, "stacks");
  c->Close();
  delete c; 
}

// draw 2d histograms
void draw_hist2(TFile *file, string name, string title, string xaxistitle,  string yaxistitle, string outputFolder, string subFolder) {
  TH1 *histo = dynamic_cast<TH1*>(file->Get(name.c_str()));
  TCanvas *c = new TCanvas(name.c_str(), title.c_str(), 1500, 1200);
  histo->GetXaxis()->SetTitle(xaxistitle.c_str());
  histo->GetYaxis()->SetTitle(yaxistitle.c_str());
  histo->Draw("COLZ");
  PrintCanvas(c, name, outputFolder, subFolder);
  c->Close();
}



void draw_zzhjj_(string sig_filename = "signal.root", string bkg_filename = "all_bkg.root",string ttbar_filename = "ttbar.root", string ttHbb_filename = "ttHbb.root", string diboson_filename = "diboson.root", string drellyan_filename = "drellyan.root", string outputFolder = "../histograms/") {

  // open files
  TFile *sig_file = TFile::Open(sig_filename.c_str(), "READ");
  TFile *bkg_file = TFile::Open(bkg_filename.c_str(), "READ");
  TFile *ttbar_file = TFile::Open(ttbar_filename.c_str(), "READ");
  TFile *ttHbb_file = TFile::Open(ttHbb_filename.c_str(), "READ");
  TFile *diboson_file = TFile::Open(diboson_filename.c_str(), "READ");
  TFile *drellyan_file = TFile::Open(drellyan_filename.c_str(), "READ");

  // draw histograms
  //gROOT->SetBatch(kTRUE);

  // cft
  draw_stack(sig_file, ttbar_file, ttHbb_file, diboson_file, drellyan_file, "hSel_reco", "", "", outputFolder);
  draw_stack(sig_file, ttbar_file, ttHbb_file, diboson_file, drellyan_file, "hSel_particle", "", "", outputFolder);
  draw_stack(sig_file, ttbar_file, ttHbb_file, diboson_file, drellyan_file, "hSel_parton", "", "", outputFolder);
  draw_stack(sig_file, ttbar_file, ttHbb_file, diboson_file, drellyan_file, "hEff_reco", "", "", outputFolder);
  draw_stack(sig_file, ttbar_file, ttHbb_file, diboson_file, drellyan_file, "hEff_particle", "", "", outputFolder);
  draw_stack(sig_file, ttbar_file, ttHbb_file, diboson_file, drellyan_file, "hEff_parton", "", "", outputFolder);

  // higgs - reco
  draw_stack(sig_file, ttbar_file, ttHbb_file, diboson_file, drellyan_file, "hbb_pT_reco", "p^{T}_{hbb}_reco", "p_{T}", outputFolder);
  draw_stack(sig_file, ttbar_file, ttHbb_file, diboson_file, drellyan_file, "hbb_m_reco", "m_{hbb}_reco", "mass (GeV)", outputFolder);
  draw_stack(sig_file, ttbar_file, ttHbb_file, diboson_file, drellyan_file, "bb_#Delta#phi_reco", "#Delta#phi_{bb}_reco", "#phi", outputFolder);
  draw_stack(sig_file, ttbar_file, ttHbb_file, diboson_file, drellyan_file, "bb_#Delta#eta_reco", "#Delta#eta_{bb}_reco", "#eta", outputFolder);
  draw_stack(sig_file, ttbar_file, ttHbb_file, diboson_file, drellyan_file, "bb_#DeltaR_reco", "#DeltaR_{bb}_reco", "R", outputFolder);  
  // higgs - particle
  draw_stack(sig_file, ttbar_file, ttHbb_file, diboson_file, drellyan_file, "hbb_pT_particle", "p^{T}_{hbb}_particle", "p_{T}", outputFolder);
  draw_stack(sig_file, ttbar_file, ttHbb_file, diboson_file, drellyan_file, "hbb_m_particle", "m_{hbb}_particle", "mass (GeV)", outputFolder);
  draw_stack(sig_file, ttbar_file, ttHbb_file, diboson_file, drellyan_file, "bb_#Delta#phi_particle", "#Delta#phi_{bb}_particle", "#phi", outputFolder);
  draw_stack(sig_file, ttbar_file, ttHbb_file, diboson_file, drellyan_file, "bb_#Delta#eta_particle", "#Delta#eta_{bb}_particle", "#eta", outputFolder);
  draw_stack(sig_file, ttbar_file, ttHbb_file, diboson_file, drellyan_file, "bb_#DeltaR_particle", "#DeltaR_{bb}_particle", "R", outputFolder);
  // higgs - parton
  draw_stack(sig_file, ttbar_file, ttHbb_file, diboson_file, drellyan_file, "hbb_pT_parton", "p^{T}_{hbb}_parton", "p_{T}", outputFolder);
  draw_stack(sig_file, ttbar_file, ttHbb_file, diboson_file, drellyan_file, "hbb_m_parton", "m_{hbb}_parton", "mass (GeV)", outputFolder);
  draw_stack(sig_file, ttbar_file, ttHbb_file, diboson_file, drellyan_file, "bb_#Delta#phi_parton", "#Delta#phi_{bb}_parton", "#phi", outputFolder);
  draw_stack(sig_file, ttbar_file, ttHbb_file, diboson_file, drellyan_file, "bb_#Delta#eta_parton", "#Delta#eta_{bb}_parton", "#eta", outputFolder);
  draw_stack(sig_file, ttbar_file, ttHbb_file, diboson_file, drellyan_file, "bb_#DeltaR_parton", "#DeltaR_{bb}_parton", "R", outputFolder);
  // vbfj - reco
  draw_stack(sig_file, ttbar_file, ttHbb_file, diboson_file, drellyan_file, "jj_pT_reco", "p^{T}_{jj}_reco", "p_{T}", outputFolder);
  draw_stack(sig_file, ttbar_file, ttHbb_file, diboson_file, drellyan_file, "jj_#Delta#phi_reco", "#Delta#phi_{jj}_reco", "#phi", outputFolder);
  draw_stack(sig_file, ttbar_file, ttHbb_file, diboson_file, drellyan_file, "jj_#Delta#eta_reco", "#Delta#eta_{jj}_reco", "#eta", outputFolder);
  draw_stack(sig_file, ttbar_file, ttHbb_file, diboson_file, drellyan_file, "jj_#DeltaR_reco", "#DeltaR_{jj}_reco", "R", outputFolder);
  // vbfj - particle
  draw_stack(sig_file, ttbar_file, ttHbb_file, diboson_file, drellyan_file, "jj_pT_particle", "p^{T}_{jj}_particle", "p_{T}", outputFolder);
  draw_stack(sig_file, ttbar_file, ttHbb_file, diboson_file, drellyan_file, "jj_#DeltaR_particle", "#DeltaR_{jj}_particle", "R", outputFolder);
  draw_stack(sig_file, ttbar_file, ttHbb_file, diboson_file, drellyan_file, "jj_#Delta#eta_particle", "#Delta#eta_{jj}_particle", "#eta", outputFolder);
  draw_stack(sig_file, ttbar_file, ttHbb_file, diboson_file, drellyan_file, "jj_#DeltaR_particle", "#DeltaR_{jj}_particle", "R", outputFolder);
  // z - reco
  draw_stack(sig_file, ttbar_file, ttHbb_file, diboson_file, drellyan_file, "z1_pT_reco", "p^{T}_{z1}_reco", "p_{T}", outputFolder);
  draw_stack(sig_file, ttbar_file, ttHbb_file, diboson_file, drellyan_file, "z2_pT_reco", "p^{T}_{z2}_reco", "p_{T}", outputFolder);
  draw_stack(sig_file, ttbar_file, ttHbb_file, diboson_file, drellyan_file, "z1_m_reco", "m_{z1}_reco", "mass (GeV)", outputFolder);
  draw_stack(sig_file, ttbar_file, ttHbb_file, diboson_file, drellyan_file, "z2_m_reco", "m_{z2}_reco", "mass (GeV)", outputFolder);
  draw_stack(sig_file, ttbar_file, ttHbb_file, diboson_file, drellyan_file, "zz_#Delta#phi_reco", "#Delta#phi_{zz}_reco", "#phi", outputFolder);
  draw_stack(sig_file, ttbar_file, ttHbb_file, diboson_file, drellyan_file, "zz_#Delta#eta_reco", "#Delta#eta_{zz}_reco", "#eta", outputFolder);
  draw_stack(sig_file, ttbar_file, ttHbb_file, diboson_file, drellyan_file, "zz_#DeltaR_reco", "#DeltaR_{zz}_reco", "R", outputFolder);
  // z - particle
  draw_stack(sig_file, ttbar_file, ttHbb_file, diboson_file, drellyan_file, "z1_pT_particle", "p^{T}_{z1}_particle", "p_{T}", outputFolder);
  draw_stack(sig_file, ttbar_file, ttHbb_file, diboson_file, drellyan_file, "z2_pT_particle", "p^{T}_{z2}_particle", "p_{T}", outputFolder);
  draw_stack(sig_file, ttbar_file, ttHbb_file, diboson_file, drellyan_file, "z1_m_particle", "m_{z1}_particle", "mass (GeV)", outputFolder);
  draw_stack(sig_file, ttbar_file, ttHbb_file, diboson_file, drellyan_file, "z2_m_particle", "m_{z2}_particle", "mass (GeV)", outputFolder);
  draw_stack(sig_file, ttbar_file, ttHbb_file, diboson_file, drellyan_file, "zz_#Delta#phi_particle", "#Delta#phi_{zz}_particle", "#phi", outputFolder);
  draw_stack(sig_file, ttbar_file, ttHbb_file, diboson_file, drellyan_file, "zz_#Delta#eta_particle", "#Delta#eta_{zz}_particle", "#eta", outputFolder);
  draw_stack(sig_file, ttbar_file, ttHbb_file, diboson_file, drellyan_file, "zz_#DeltaR_particle", "#DeltaR_{zz}_particle", "R", outputFolder);
  // z - parton
  draw_stack(sig_file, ttbar_file, ttHbb_file, diboson_file, drellyan_file, "z1_pT_parton", "p^{T}_{z1}_parton", "p_{T}", outputFolder);
  draw_stack(sig_file, ttbar_file, ttHbb_file, diboson_file, drellyan_file, "z2_pT_parton", "p^{T}_{z2}_parton", "p_{T}", outputFolder);
  draw_stack(sig_file, ttbar_file, ttHbb_file, diboson_file, drellyan_file, "z1_m_parton", "m_{z1}_parton", "mass (GeV)", outputFolder);
  draw_stack(sig_file, ttbar_file, ttHbb_file, diboson_file, drellyan_file, "z2_m_parton", "m_{z2}_parton", "mass (GeV)", outputFolder);
  draw_stack(sig_file, ttbar_file, ttHbb_file, diboson_file, drellyan_file, "zz_#Delta#phi_parton", "#Delta#phi_{zz}_parton", "#phi", outputFolder);
  draw_stack(sig_file, ttbar_file, ttHbb_file, diboson_file, drellyan_file, "zz_#Delta#eta_parton", "#Delta#eta_{zz}_parton", "#eta", outputFolder);
  draw_stack(sig_file, ttbar_file, ttHbb_file, diboson_file, drellyan_file, "zz_#DeltaR_parton", "#DeltaR_{zz}_parton", "R", outputFolder);
  // w - reco
  draw_stack(sig_file, ttbar_file, ttHbb_file, diboson_file, drellyan_file, "ll_pT_reco", "p^{T}_{ll}_reco", "p_{T}", outputFolder);
  draw_stack(sig_file, ttbar_file, ttHbb_file, diboson_file, drellyan_file, "ll_m_reco", "m_{ll}_reco", "mass (GeV)", outputFolder);
  draw_stack(sig_file, ttbar_file, ttHbb_file, diboson_file, drellyan_file, "w1_pT_reco", "p^{T}_{w1}_reco", "p_{T}", outputFolder);
  draw_stack(sig_file, ttbar_file, ttHbb_file, diboson_file, drellyan_file, "w1_m_reco", "m_{w1}_reco", "mass (GeV)", outputFolder);
  draw_stack(sig_file, ttbar_file, ttHbb_file, diboson_file, drellyan_file, "w2_pT_reco", "p^{T}_{w2}_reco", "p_{T}", outputFolder);
  draw_stack(sig_file, ttbar_file, ttHbb_file, diboson_file, drellyan_file, "w2_m_reco", "m_{w2}_reco", "mass (GeV)", outputFolder);
  draw_stack(sig_file, ttbar_file, ttHbb_file, diboson_file, drellyan_file, "ww_pT_reco", "p^{T}_{ww}_reco", "p_{T}", outputFolder);
  draw_stack(sig_file, ttbar_file, ttHbb_file, diboson_file, drellyan_file, "ww_m_reco", "m_{ww}_reco", "mass (GeV)", outputFolder);
  // w - particle 
  draw_stack(sig_file, ttbar_file, ttHbb_file, diboson_file, drellyan_file, "ll_pT_particle", "p^{T}_{ll}_particle", "p_{T}", outputFolder);
  draw_stack(sig_file, ttbar_file, ttHbb_file, diboson_file, drellyan_file, "ll_m_particle", "m_{ll}_particle", "mass (GeV)", outputFolder);
  draw_stack(sig_file, ttbar_file, ttHbb_file, diboson_file, drellyan_file, "w1_pT_particle", "p^{T}_{w1}_particle", "p_{T}", outputFolder);
  draw_stack(sig_file, ttbar_file, ttHbb_file, diboson_file, drellyan_file, "w1_m_particle", "m_{w1}_particle", "mass (GeV)", outputFolder);
  draw_stack(sig_file, ttbar_file, ttHbb_file, diboson_file, drellyan_file, "w2_pT_particle", "p^{T}_{w2}_particle", "p_{T}", outputFolder);
  draw_stack(sig_file, ttbar_file, ttHbb_file, diboson_file, drellyan_file, "w2_m_particle", "m_{w2}_particle", "mass (GeV)", outputFolder);
  draw_stack(sig_file, ttbar_file, ttHbb_file, diboson_file, drellyan_file, "ww_pT_particle", "p^{T}_{ww}_particle", "p_{T}", outputFolder);
  draw_stack(sig_file, ttbar_file, ttHbb_file, diboson_file, drellyan_file, "ww_m_particle", "m_{ww}_particle", "mass (GeV)", outputFolder);
  // leps - reco
  draw_stack(sig_file, ttbar_file, ttHbb_file, diboson_file, drellyan_file, "l1l2_#Delta#phi_reco", "#Delta#phi_{l1l2}_reco", "#phi", outputFolder);
  draw_stack(sig_file, ttbar_file, ttHbb_file, diboson_file, drellyan_file, "l3l4_#Delta#phi_reco", "#Delta#phi_{l3l4}_reco", "#phi", outputFolder);
  draw_stack(sig_file, ttbar_file, ttHbb_file, diboson_file, drellyan_file, "l1l2_#Delta#eta_reco", "#Delta#eta_{l1l2}_reco", "#eta", outputFolder);
  draw_stack(sig_file, ttbar_file, ttHbb_file, diboson_file, drellyan_file, "l3l4_#Delta#eta_reco", "#Delta#eta_{l3l4}_reco", "#eta", outputFolder);
  draw_stack(sig_file, ttbar_file, ttHbb_file, diboson_file, drellyan_file, "l1l2_#Delta R_reco", "#Delta R_{l1l2}_reco", "R", outputFolder);
  draw_stack(sig_file, ttbar_file, ttHbb_file, diboson_file, drellyan_file, "l3l4_#Delta R_reco", "#Delta R_{l3l4}_reco", "R", outputFolder);
  // leps - particle
  draw_stack(sig_file, ttbar_file, ttHbb_file, diboson_file, drellyan_file, "l1l2_#Delta#phi_particle", "#Delta#phi_{l1l2}_particle", "#phi", outputFolder);
  draw_stack(sig_file, ttbar_file, ttHbb_file, diboson_file, drellyan_file, "l3l4_#Delta#phi_particle", "#Delta#phi_{l3l4}_particle", "#phi", outputFolder);
  draw_stack(sig_file, ttbar_file, ttHbb_file, diboson_file, drellyan_file, "l1l2_#Delta#eta_particle", "#Delta#eta_{l1l2}_particle", "#eta", outputFolder);
  draw_stack(sig_file, ttbar_file, ttHbb_file, diboson_file, drellyan_file, "l3l4_#Delta#eta_particle", "#Delta#eta_{l3l4}_particle", "#eta", outputFolder);
  draw_stack(sig_file, ttbar_file, ttHbb_file, diboson_file, drellyan_file, "l1l2_#Delta R_particle", "#Delta R_{l1l2}_particle", "R", outputFolder);
  draw_stack(sig_file, ttbar_file, ttHbb_file, diboson_file, drellyan_file, "l3l4_#Delta R_particle", "#Delta R_{l3l4}_particle", "R", outputFolder);
  // leps - parton
  draw_stack(sig_file, ttbar_file, ttHbb_file, diboson_file, drellyan_file, "l1l2_#Delta#phi_parton", "#Delta#phi_{l1l2}_parton", "#phi", outputFolder);
  draw_stack(sig_file, ttbar_file, ttHbb_file, diboson_file, drellyan_file, "l3l4_#Delta#phi_parton", "#Delta#phi_{l3l4}_parton", "#phi", outputFolder);
  draw_stack(sig_file, ttbar_file, ttHbb_file, diboson_file, drellyan_file, "l1l2_#Delta#eta_parton", "#Delta#eta_{l1l2}_parton", "#eta", outputFolder);
  draw_stack(sig_file, ttbar_file, ttHbb_file, diboson_file, drellyan_file, "l3l4_#Delta#eta_parton", "#Delta#eta_{l3l4}_parton", "#eta", outputFolder);
  draw_stack(sig_file, ttbar_file, ttHbb_file, diboson_file, drellyan_file, "l1l2_#Delta R_parton", "#Delta R_{l1l2}_parton", "R", outputFolder);
  draw_stack(sig_file, ttbar_file, ttHbb_file, diboson_file, drellyan_file, "l3l4_#Delta R_parton", "#Delta R_{l3l4}_parton", "R", outputFolder);
  /* 
  // higgs comps
  draw_hist2(sig_file, "H_pT_comp_12", "p_{T}^{hbb}", "Parton Level p_{T}", "Particle Level p_{T}", outputFolder, "signal2D");
  draw_hist2(sig_file, "H_pT_comp_23", "p_{T}^{hbb}", "Particle Level p_{T}", "Reco Level p_{T}", outputFolder, "signal2D");
  draw_hist2(sig_file, "H_pT_comp_13", "p_{T}^{hbb}", "Parton Level p_{T}", "Reco Level p_{T}", outputFolder, "signal2D");
  draw_hist2(sig_file, "H_m_comp_12", "m_{hbb}", "Parton Level Mass (GeV)", "Particle Level Mass (GeV)", outputFolder, "signal2D");
  draw_hist2(sig_file, "H_m_comp_23", "m_{hbb}", "Particle Level Mass (GeV)", "Reco Level Mass (GeV)", outputFolder, "signal2D");
  draw_hist2(sig_file, "H_m_comp_13", "m_{hbb}", "Parton Level Mass (GeV)", "Reco Level Mass (GeV)", outputFolder, "signal2D");
  draw_hist2(sig_file, "bb_#Delta#phi_comp_12", "#Delta#phi_{bb}", "Parton Level #Delta#phi", "Particle Level #Delta#phi", outputFolder, "signal2D");
  draw_hist2(sig_file, "bb_#Delta#phi_comp_23", "#Delta#phi_{bb}", "Particle Level #Delta#phi", "Reco Level #Delta#phi", outputFolder, "signal2D");
  draw_hist2(sig_file, "bb_#Delta#phi_comp_13", "#Delta#phi_{bb}", "Parton Level #Delta#phi", "Reco Level #Delta#phi", outputFolder, "signal2D");
  draw_hist2(sig_file, "bb_#Delta#eta_comp_12", "#Delta#eta_{bb}", "Parton Level #Delta#eta", "Particle Level #Delta#eta", outputFolder, "signal2D");
  draw_hist2(sig_file, "bb_#Delta#eta_comp_23", "#Delta#eta_{bb}", "Particle Level #Delta#eta", "Reco Level #Delta#eta", outputFolder, "signal2D");
  draw_hist2(sig_file, "bb_#Delta#eta_comp_13", "#Delta#eta_{bb}", "Parton Level #Delta#eta", "Reco Level #Delta#eta", outputFolder, "signal2D");
  // vbfj comps
  draw_hist2(sig_file, "jj_pT_comp_12", "p_{T}^{jj}", "Parton Level p_{T}", "Particle Level p_{T}", outputFolder, "signal2D");
  draw_hist2(sig_file, "jj_pT_comp_23", "p_{T}^{jj}", "Particle Level p_{T}", "Reco Level p_{T}", outputFolder, "signal2D");
  draw_hist2(sig_file, "jj_pT_comp_13", "p_{T}^{jj}", "Parton Level p_{T}", "Reco Level p_{T}", outputFolder, "signal2D");
  draw_hist2(sig_file, "jj_#Delta#phi_comp_12", "#Delta#phi_{jj}", "Parton Level #Delta#phi", "Particle Level #Delta#phi", outputFolder, "signal2D");
  draw_hist2(sig_file, "jj_#Delta#phi_comp_23", "#Delta#phi_{jj}", "Particle Level #Delta#phi", "Reco Level #Delta#phi", outputFolder, "signal2D");
  draw_hist2(sig_file, "jj_#Delta#phi_comp_13", "#Delta#phi_{jj}", "Parton Level #Delta#phi", "Reco Level #Delta#phi", outputFolder, "signal2D");
  */

  draw_stack(sig_file, ttbar_file, ttHbb_file, diboson_file, drellyan_file, "ll_m_reco_0_15", "m_{ll}_reco_0_15", "mass (GeV)", outputFolder);
  draw_stack(sig_file, ttbar_file, ttHbb_file, diboson_file, drellyan_file, "ll_m_particle_0_15", "m_{ll}_particle_0_15", "mass (GeV)", outputFolder);
  draw_stack(sig_file, ttbar_file, ttHbb_file, diboson_file, drellyan_file, "reco_event_type", "ET", "ET", outputFolder);
  draw_stack(sig_file, ttbar_file, ttHbb_file, diboson_file, drellyan_file, "particle_event_type", "ET", "ET", outputFolder);
  // print signal significance
  
  std::cout<<"Reco Level Contributions:"<<std::endl;
  
  TH1 *sig_mbb = dynamic_cast<TH1*>(sig_file->Get("hbb_m_reco"));
  TH1 *bkg_mbb = dynamic_cast<TH1*>(bkg_file->Get("hbb_m_reco"));
  TH1 *ttbar_mbb = dynamic_cast<TH1*>(ttbar_file->Get("hbb_m_reco"));
  TH1 *ttHbb_mbb = dynamic_cast<TH1*>(ttHbb_file->Get("hbb_m_reco"));
  TH1 *diboson_mbb = dynamic_cast<TH1*>(diboson_file->Get("hbb_m_reco"));
  TH1 *drellyan_mbb = dynamic_cast<TH1*>(drellyan_file->Get("hbb_m_reco"));

  std::cout<<"ttbar: "<<ttbar_mbb->Integral()<<std::endl;
  std::cout<<"ttHbb: "<<ttHbb_mbb->Integral()<<std::endl;
  std::cout<<"diboson: "<<diboson_mbb->Integral()<<std::endl;
  std::cout<<"drell-yan: "<<drellyan_mbb->Integral()<<std::endl;
  std::cout<<"signal: "<<sig_mbb->Integral()<<std::endl;
  
  double sig_sig = sig_mbb->Integral()/sqrt(bkg_mbb->Integral());
  std::cout<<"Signal Significance: "<<sig_sig<<std::endl;

  std::cout<<"Particle Level Contributions:"<<std::endl;
  
  TH1 *sig_mbb_p = dynamic_cast<TH1*>(sig_file->Get("hbb_m_particle"));
  TH1 *bkg_mbb_p = dynamic_cast<TH1*>(bkg_file->Get("hbb_m_particle"));
  TH1 *ttbar_mbb_p = dynamic_cast<TH1*>(ttbar_file->Get("hbb_m_particle"));
  TH1 *ttHbb_mbb_p = dynamic_cast<TH1*>(ttHbb_file->Get("hbb_m_particle"));
  TH1 *diboson_mbb_p = dynamic_cast<TH1*>(diboson_file->Get("hbb_m_particle"));
  TH1 *drellyan_mbb_p = dynamic_cast<TH1*>(drellyan_file->Get("hbb_m_particle"));

  std::cout<<"ttbar: "<<ttbar_mbb_p->Integral()<<std::endl;
  std::cout<<"ttHbb: "<<ttHbb_mbb_p->Integral()<<std::endl;
  std::cout<<"diboson: "<<diboson_mbb_p->Integral()<<std::endl;
  std::cout<<"drell-yan: "<<drellyan_mbb_p->Integral()<<std::endl;
  std::cout<<"signal: "<<sig_mbb_p->Integral()<<std::endl;
  
  double sig_sig_p = sig_mbb_p->Integral()/sqrt(bkg_mbb_p->Integral());
  std::cout<<"Signal Significance: "<<sig_sig_p<<std::endl;

  std::cout<<"Integral before cuts with Lumi=300:"<<std::endl;
  
  TH1 *sig_w = dynamic_cast<TH1*>(sig_file->Get("weights"));
  TH1 *bkg_w = dynamic_cast<TH1*>(bkg_file->Get("weights"));
  TH1 *ttbar_w = dynamic_cast<TH1*>(ttbar_file->Get("weights"));
  TH1 *ttHbb_w = dynamic_cast<TH1*>(ttHbb_file->Get("weights"));
  TH1 *diboson_w = dynamic_cast<TH1*>(diboson_file->Get("weights"));
  TH1 *drellyan_w = dynamic_cast<TH1*>(drellyan_file->Get("weights"));

  std::cout<<"ttbar: "<<ttbar_w->Integral()<<std::endl;
  std::cout<<"ttHbb: "<<ttHbb_w->Integral()<<std::endl;
  std::cout<<"diboson: "<<diboson_w->Integral()<<std::endl;
  std::cout<<"drell-yan: "<<drellyan_w->Integral()<<std::endl;
  std::cout<<"signal: "<<sig_w->Integral()<<std::endl;

  /*TFile *test_file = TFile::Open("../outputs/histograms/ttbar012j/delphes_704728_5.root", "READ");
  TH1 *test_w = dynamic_cast<TH1*>(test_file->Get("weights"));
  std::cout<<"single file ttbar integral: "<<test_w->Integral()<<std::endl;
  test_file->Close();*/
  
  sig_file->Close();
  bkg_file->Close();
  ttbar_file->Close();
  ttHbb_file->Close();
  diboson_file->Close();
  drellyan_file->Close();

}


void draw_zzhjj(string sel=""){
  std::map<string,string> selMapToHistsSignal;
  selMapToHistsSignal["HWWJJ"] = "wpwmhqq_HWWJJ.root";
  selMapToHistsSignal["HZZJJ"] = "zzhqq_HZZJJ.root";

  system(Form("mkdir -p ../histograms/%s/stacks",sel.c_str()));

  draw_zzhjj_(selMapToHistsSignal[sel],
	     Form("all_bkg_%s.root",sel.c_str()),
	     Form("ttbar_%s.root",sel.c_str()),
	     Form("ttHbb_%s.root",sel.c_str()),
	     Form("diboson_%s.root",sel.c_str()),
	     Form("drellyan_%s.root",sel.c_str()),
	     Form("../histograms/%s",sel.c_str()));
	     
	     
}

