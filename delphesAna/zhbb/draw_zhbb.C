#ifdef __CLING__
R__LOAD_LIBRARY(libDelphes)
#include <iostream>
#include <fstream>
#include <cmath>
#endif

// save plots as different file types
void PrintCanvas(TCanvas *c=nullptr, string name="default", string outputFolder="default", string subFolder="default"){
  std::vector <string> types={"jpg"}; 
  for(std::vector<string>::iterator it=types.begin(); it!=types.end(); it++) {
    c->Print(Form("%s/%s/%s.%s", outputFolder.c_str(), subFolder.c_str(), name.c_str(), (*it).c_str()), (*it).c_str());
  }
}

// draw 1d histograms
void draw_stack(TFile *sig_file, TFile *ttbar_file, TFile *ttHbb_file, TFile *diboson_file, TFile *drellyan_file, const char *name, const char *title, const char *axistitle, const char *outputFolder) {
  // get the histograms from the files
  TH1F *sig_hist = (TH1F*)sig_file->Get(name);
  TH1F *ttbar_hist = (TH1F*)ttbar_file->Get(name);
  TH1F *ttHbb_hist = (TH1F*)ttHbb_file->Get(name);
  TH1F *diboson_hist = (TH1F*)diboson_file->Get(name);
  TH1F *drellyan_hist = (TH1F*)drellyan_file->Get(name);
  //std::cout<<"Name "<<name<<std::endl;
  sig_hist->SetFillColor(kCyan);
  ttbar_hist->SetFillColor(kOrange);
  ttHbb_hist->SetFillColor(kMagenta);
  diboson_hist->SetFillColor(kRed);
  drellyan_hist->SetFillColor(kYellow);
  // rebin
  sig_hist->Rebin(2);
  ttbar_hist->Rebin(2);
  ttHbb_hist->Rebin(2);
  diboson_hist->Rebin(2);
  drellyan_hist->Rebin(2);
  // form the histogram stack
  THStack *stack = new THStack(name, title);
  stack->Add(diboson_hist);
  stack->Add(drellyan_hist);
  stack->Add(ttbar_hist);
  stack->Add(ttHbb_hist);
  stack->Add(sig_hist);
  // make a legend
  TLegend *legend = new TLegend(0.8, 0.8, 0.9, 0.9);
  legend->AddEntry(sig_hist, "Signal", "f");
  legend->AddEntry(drellyan_hist, "drell-yan", "f");
  legend->AddEntry(ttbar_hist, "ttbar", "f");
  legend->AddEntry(ttHbb_hist, "ttH", "f");
  legend->AddEntry(diboson_hist, "diboson", "f");
  // make a canvas and draw on it
  TCanvas *c = new TCanvas(name, title, 1500, 1200);
  c->cd();
  stack->Draw("hist e");
  stack->GetXaxis()->SetTitle(axistitle);
  legend->Draw();
  PrintCanvas(c, name, outputFolder, "stacks");
  c->Close();
}

// draw 2d histograms
void draw_hist2(TFile *file, const char *name, const char *title, const char *xaxistitle, const char *yaxistitle, const char *outputFolder, const char *subFolder) {
  TH1 *histo = dynamic_cast<TH1*>(file->Get(name));
  TCanvas *c = new TCanvas(name, title, 1500, 1200);
  histo->GetXaxis()->SetTitle(xaxistitle);
  histo->GetYaxis()->SetTitle(yaxistitle);
  histo->Draw("COLZ");
  PrintCanvas(c, name, outputFolder, subFolder);
  c->Close();
}

void draw_zhbb(const char *sig_filename, const char *bkg_filename, const char *ttbar_filename, const char *ttHbb_filename, const char *diboson_filename, const char *drellyan_filename, const char *outputFolder) {\
  TFile *sig_file = TFile::Open(sig_filename, "READ");
  TFile *bkg_file = TFile::Open(bkg_filename, "READ");
  TFile *ttbar_file = TFile::Open(ttbar_filename, "READ");
  TFile *ttHbb_file = TFile::Open(ttHbb_filename, "READ");
  TFile *diboson_file = TFile::Open(diboson_filename, "READ");
  TFile *drellyan_file = TFile::Open(drellyan_filename, "READ");
  // draw histograms
  gROOT->SetBatch(kTRUE);
  // mass
  draw_stack(sig_file, ttbar_file, ttHbb_file, diboson_file, drellyan_file, "mass_ZH", "m_{T}^{ZH}", "mass (GeV)", outputFolder);
  draw_stack(sig_file, ttbar_file, ttHbb_file, diboson_file, drellyan_file, "mass_bb_reco", "Reco m_{bb}", "mass (GeV)", outputFolder);
  draw_stack(sig_file, ttbar_file, ttHbb_file, diboson_file, drellyan_file, "mass_ll_reco", "Reco m_{ll}", "mass (GeV)", outputFolder);
  draw_stack(sig_file, ttbar_file, ttHbb_file, diboson_file, drellyan_file, "mass_bbll_reco", "Reco m_{bbll}", "mass (GeV)", outputFolder);
  draw_stack(sig_file, ttbar_file, ttHbb_file, diboson_file, drellyan_file, "mass_bb_particle", "m_{bb}", "mass (GeV)", outputFolder);
  draw_stack(sig_file, ttbar_file, ttHbb_file, diboson_file, drellyan_file, "mass_ll_particle", "m_{ll}", "mass (GeV)", outputFolder);
  draw_stack(sig_file, ttbar_file, ttHbb_file, diboson_file, drellyan_file, "mass_bbll_particle", "m_{bbll}", "mass (GeV)", outputFolder);
  // pt
  draw_stack(sig_file, ttbar_file, ttHbb_file, diboson_file, drellyan_file, "pt_H", "p_{T}^{H}", "p_{T}", outputFolder);
  draw_stack(sig_file, ttbar_file, ttHbb_file, diboson_file, drellyan_file, "pt_Z", "p_{T}^{Z}", "p_{T}", outputFolder);
  draw_stack(sig_file, ttbar_file, ttHbb_file, diboson_file, drellyan_file, "pt_LE",  "Lead p_{T}^{e}", "p_{T}", outputFolder);
  draw_stack(sig_file, ttbar_file, ttHbb_file, diboson_file, drellyan_file, "pt_LM",  "Lead p_{T}^{#mu}", "p_{T}", outputFolder);
  draw_stack(sig_file, ttbar_file, ttHbb_file, diboson_file, drellyan_file, "pt_SLE",  "Sublead p_{T}^{e}", "p_{T}", outputFolder);
  draw_stack(sig_file, ttbar_file, ttHbb_file, diboson_file, drellyan_file, "pt_SLM",  "Sublead p_{T}^{#mu}", "p_{T}", outputFolder);
  draw_stack(sig_file, ttbar_file, ttHbb_file, diboson_file, drellyan_file, "pt_LB", "Lead b Parton p_{T}", "p_{T}", outputFolder);
  draw_stack(sig_file, ttbar_file, ttHbb_file, diboson_file, drellyan_file, "pt_SLB", "Sublead b Parton p_{T}", "p_{T}", outputFolder);
  draw_stack(sig_file, ttbar_file, ttHbb_file, diboson_file, drellyan_file, "pt_ZH", "p_{T}^{ZH}", "p_{T}", outputFolder);
  draw_stack(sig_file, ttbar_file, ttHbb_file, diboson_file, drellyan_file, "pt_bb_reco", "Reco p_{T}^{bb}", "p_{T}", outputFolder);
  draw_stack(sig_file, ttbar_file, ttHbb_file, diboson_file, drellyan_file, "pt_ll_reco", "Reco p_{T}^{ll}", "p_{T}", outputFolder);
  draw_stack(sig_file, ttbar_file, ttHbb_file, diboson_file, drellyan_file, "pt_LE_reco",  "Lead Reco p_{T}^{e}", "p_{T}", outputFolder);
  draw_stack(sig_file, ttbar_file, ttHbb_file, diboson_file, drellyan_file, "pt_LM_reco",  "Lead Reco p_{T}^{#mu}", "p_{T}", outputFolder);
  draw_stack(sig_file, ttbar_file, ttHbb_file, diboson_file, drellyan_file, "pt_SLE_reco",  "Sublead Reco p_{T}^{e}", "p_{T}", outputFolder);
  draw_stack(sig_file, ttbar_file, ttHbb_file, diboson_file, drellyan_file, "pt_SLM_reco",  "Sublead Reco p_{T}^{#mu}", "p_{T}", outputFolder);
  draw_stack(sig_file, ttbar_file, ttHbb_file, diboson_file, drellyan_file, "pt_LB_reco", "Lead Reco b Jet p_{T}", "p_{T}", outputFolder);
  draw_stack(sig_file, ttbar_file, ttHbb_file, diboson_file, drellyan_file, "pt_SLB_reco", "Sublead Reco b Jet p_{T}", "p_{T}", outputFolder);
  draw_stack(sig_file, ttbar_file, ttHbb_file, diboson_file, drellyan_file, "pt_bbll_reco", "Reco p_{T}^{bbll}", "p_{T}", outputFolder);
  draw_stack(sig_file, ttbar_file, ttHbb_file, diboson_file, drellyan_file, "pt_bb_particle", "p_{T}^{bb}", "p_{T}", outputFolder);
  draw_stack(sig_file, ttbar_file, ttHbb_file, diboson_file, drellyan_file, "pt_ll_particle", "p_{T}^{ll}", "p_{T}", outputFolder);
  draw_stack(sig_file, ttbar_file, ttHbb_file, diboson_file, drellyan_file, "pt_LE_particle",  "Lead p_{T}^{e}", "p_{T}", outputFolder);
  draw_stack(sig_file, ttbar_file, ttHbb_file, diboson_file, drellyan_file, "pt_LM_particle",  "Lead p_{T}^{#mu}", "p_{T}", outputFolder);
  draw_stack(sig_file, ttbar_file, ttHbb_file, diboson_file, drellyan_file, "pt_SLE_particle",  "Sublead p_{T}^{e}", "p_{T}", outputFolder);
  draw_stack(sig_file, ttbar_file, ttHbb_file, diboson_file, drellyan_file, "pt_SLM_particle",  "Sublead p_{T}^{#mu}", "p_{T}", outputFolder);
  draw_stack(sig_file, ttbar_file, ttHbb_file, diboson_file, drellyan_file, "pt_LB_particle", "Lead b Jet p_{T}", "p_{T}", outputFolder);
  draw_stack(sig_file, ttbar_file, ttHbb_file, diboson_file, drellyan_file, "pt_SLB_particle", "Sublead b Jet p_{T}", "p_{T}", outputFolder);
  draw_stack(sig_file, ttbar_file, ttHbb_file, diboson_file, drellyan_file, "pt_bbll_particle", "p_{T}^{bbll}", "p_{T}", outputFolder);
  // eta
  draw_stack(sig_file, ttbar_file, ttHbb_file, diboson_file, drellyan_file, "eta_H", "#eta_{H}", "#eta", outputFolder);
  draw_stack(sig_file, ttbar_file, ttHbb_file, diboson_file, drellyan_file, "eta_Z", "#eta_{Z}", "#eta", outputFolder);
  draw_stack(sig_file, ttbar_file, ttHbb_file, diboson_file, drellyan_file, "eta_LE",  "Lead #eta_{e}", "#eta", outputFolder);
  draw_stack(sig_file, ttbar_file, ttHbb_file, diboson_file, drellyan_file, "eta_LM",  "Lead #eta_{#mu}", "#eta", outputFolder);
  draw_stack(sig_file, ttbar_file, ttHbb_file, diboson_file, drellyan_file, "eta_SLE",  "Sublead #eta_{e}", "#eta", outputFolder);
  draw_stack(sig_file, ttbar_file, ttHbb_file, diboson_file, drellyan_file, "eta_SLM",  "Sublead #eta_{#mu}", "#eta", outputFolder);
  draw_stack(sig_file, ttbar_file, ttHbb_file, diboson_file, drellyan_file, "eta_LB", "Lead #eta_{b}", "#eta", outputFolder);
  draw_stack(sig_file, ttbar_file, ttHbb_file, diboson_file, drellyan_file, "eta_SLB", "Sublead #eta_{b}", "#eta", outputFolder);
  draw_stack(sig_file, ttbar_file, ttHbb_file, diboson_file, drellyan_file, "eta_ZH", "#eta_{ZH}", "#eta", outputFolder);
  draw_stack(sig_file, ttbar_file, ttHbb_file, diboson_file, drellyan_file, "deta_EE", "#Delta#eta_{ee}", "#Delta#eta", outputFolder);
  draw_stack(sig_file, ttbar_file, ttHbb_file, diboson_file, drellyan_file, "deta_MM", "#Delta#eta_{#mu#mu}", "#Delta#eta", outputFolder);
  draw_stack(sig_file, ttbar_file, ttHbb_file, diboson_file, drellyan_file, "eta_bb_reco", "Reco #eta_{bb}", "#eta", outputFolder);
  draw_stack(sig_file, ttbar_file, ttHbb_file, diboson_file, drellyan_file, "eta_ll_reco", "Reco #eta_{ll}", "#eta", outputFolder);
  draw_stack(sig_file, ttbar_file, ttHbb_file, diboson_file, drellyan_file, "eta_LE_reco",  "Lead Reco #eta_{e}", "#eta", outputFolder);
  draw_stack(sig_file, ttbar_file, ttHbb_file, diboson_file, drellyan_file, "eta_LM_reco",  "Lead Reco #eta_{#mu}", "#eta", outputFolder);
  draw_stack(sig_file, ttbar_file, ttHbb_file, diboson_file, drellyan_file, "eta_SLE_reco",  "Sublead Reco #eta_{e}", "#eta", outputFolder);
  draw_stack(sig_file, ttbar_file, ttHbb_file, diboson_file, drellyan_file, "eta_SLM_reco",  "Sublead Reco #eta_{#mu}", "#eta", outputFolder);
  draw_stack(sig_file, ttbar_file, ttHbb_file, diboson_file, drellyan_file, "eta_LB_reco", "Lead Reco b Jet #eta", "#eta", outputFolder);
  draw_stack(sig_file, ttbar_file, ttHbb_file, diboson_file, drellyan_file, "eta_SLB_reco", "Sublead Reco b Jet #eta", "#eta", outputFolder);
  draw_stack(sig_file, ttbar_file, ttHbb_file, diboson_file, drellyan_file, "eta_bbll_reco", "Reco #eta_{bbll}", "#eta", outputFolder);
  draw_stack(sig_file, ttbar_file, ttHbb_file, diboson_file, drellyan_file, "deta_EE_reco", "Reco #Delta#eta_{ee}", "#Delta#eta", outputFolder);
  draw_stack(sig_file, ttbar_file, ttHbb_file, diboson_file, drellyan_file, "deta_MM_reco", "Reco #Delta#eta_{#mu#mu}", "#Delta#eta", outputFolder);
  draw_stack(sig_file, ttbar_file, ttHbb_file, diboson_file, drellyan_file, "eta_bb_particle", "#eta_{bb}", "#eta", outputFolder);
  draw_stack(sig_file, ttbar_file, ttHbb_file, diboson_file, drellyan_file, "eta_ll_particle", "#eta_{ll}", "#eta", outputFolder);
  draw_stack(sig_file, ttbar_file, ttHbb_file, diboson_file, drellyan_file, "eta_LE_particle",  "Lead #eta_{e}", "#eta", outputFolder);
  draw_stack(sig_file, ttbar_file, ttHbb_file, diboson_file, drellyan_file, "eta_LM_particle",  "Lead #eta_{#mu}", "#eta", outputFolder);
  draw_stack(sig_file, ttbar_file, ttHbb_file, diboson_file, drellyan_file, "eta_SLE_particle",  "Sublead #eta_{e}", "#eta", outputFolder);
  draw_stack(sig_file, ttbar_file, ttHbb_file, diboson_file, drellyan_file, "eta_SLM_particle",  "Sublead #eta_{#mu}", "#eta", outputFolder);
  draw_stack(sig_file, ttbar_file, ttHbb_file, diboson_file, drellyan_file, "eta_LB_particle", "Lead b Jet #eta", "#eta", outputFolder);
  draw_stack(sig_file, ttbar_file, ttHbb_file, diboson_file, drellyan_file, "eta_SLB_particle", "Sublead b Jet #eta", "#eta", outputFolder);
  draw_stack(sig_file, ttbar_file, ttHbb_file, diboson_file, drellyan_file, "eta_bbll_particle", "#eta_{bbll}", "#eta", outputFolder);
  draw_stack(sig_file, ttbar_file, ttHbb_file, diboson_file, drellyan_file, "deta_EE_particle", "#Delta#eta_{ee}", "#Delta#eta", outputFolder);
  draw_stack(sig_file, ttbar_file, ttHbb_file, diboson_file, drellyan_file, "deta_MM_particle", "#Delta#eta_{#mu#mu}", "#Delta#eta", outputFolder);
  // cos theta
  draw_stack(sig_file, ttbar_file, ttHbb_file, diboson_file, drellyan_file, "costheta_LE",  "Lead cos#theta_{e}", "cos#theta", outputFolder);
  draw_stack(sig_file, ttbar_file, ttHbb_file, diboson_file, drellyan_file, "costheta_LM",  "Lead cos#theta_{#mu}", "cos#theta", outputFolder);
  draw_stack(sig_file, ttbar_file, ttHbb_file, diboson_file, drellyan_file, "costheta_SLE",  "Sublead cos#theta_{e}", "cos#theta", outputFolder);
  draw_stack(sig_file, ttbar_file, ttHbb_file, diboson_file, drellyan_file, "costheta_SLM",  "Sublead cos#theta_{#mu}", "cos#theta", outputFolder);
  draw_stack(sig_file, ttbar_file, ttHbb_file, diboson_file, drellyan_file, "costheta_LB", "Lead cos#theta_{b}", "cos#theta", outputFolder);
  draw_stack(sig_file, ttbar_file, ttHbb_file, diboson_file, drellyan_file, "costheta_SLB", "Sublead cos#theta_{b}", "cos#theta", outputFolder);
  draw_stack(sig_file, ttbar_file, ttHbb_file, diboson_file, drellyan_file, "costheta_LE_reco",  "Lead Reco cos#theta_{e}", "cos#theta", outputFolder);
  draw_stack(sig_file, ttbar_file, ttHbb_file, diboson_file, drellyan_file, "costheta_LM_reco",  "Lead Reco cos#theta_{#mu}", "cos#theta", outputFolder);
  draw_stack(sig_file, ttbar_file, ttHbb_file, diboson_file, drellyan_file, "costheta_SLE_reco",  "Sublead Reco cos#theta_{e}", "cos#theta", outputFolder);
  draw_stack(sig_file, ttbar_file, ttHbb_file, diboson_file, drellyan_file, "costheta_SLM_reco",  "Sublead Reco cos#theta_{#mu}", "cos#theta", outputFolder);
  draw_stack(sig_file, ttbar_file, ttHbb_file, diboson_file, drellyan_file, "costheta_LB_reco", "Lead Reco b Jet cos#theta", "cos#theta", outputFolder);
  draw_stack(sig_file, ttbar_file, ttHbb_file, diboson_file, drellyan_file, "costheta_SLB_reco", "Sublead Reco b Jet cos#theta", "cos#theta", outputFolder);
  draw_stack(sig_file, ttbar_file, ttHbb_file, diboson_file, drellyan_file, "costheta_LE_particle",  "Lead cos#theta_{e}", "cos#theta", outputFolder);
  draw_stack(sig_file, ttbar_file, ttHbb_file, diboson_file, drellyan_file, "costheta_LM_particle",  "Lead cos#theta_{#mu}", "cos#theta", outputFolder);
  draw_stack(sig_file, ttbar_file, ttHbb_file, diboson_file, drellyan_file, "costheta_SLE_particle",  "Sublead cos#theta_{e}", "cos#theta", outputFolder);
  draw_stack(sig_file, ttbar_file, ttHbb_file, diboson_file, drellyan_file, "costheta_SLM_particle",  "Sublead cos#theta_{#mu}", "cos#theta", outputFolder);
  draw_stack(sig_file, ttbar_file, ttHbb_file, diboson_file, drellyan_file, "costheta_LB_particle", "Lead b Jet cos#theta", "cos#theta", outputFolder);
  draw_stack(sig_file, ttbar_file, ttHbb_file, diboson_file, drellyan_file, "costheta_SLB_particle", "Sublead b Jet cos#theta", "cos#theta", outputFolder);
  // phi
  draw_stack(sig_file, ttbar_file, ttHbb_file, diboson_file, drellyan_file, "phi_H", "#phi_{H}", "#phi (Rad)", outputFolder);
  draw_stack(sig_file, ttbar_file, ttHbb_file, diboson_file, drellyan_file, "phi_Z", "#phi_{Z}", "#phi (Rad)", outputFolder);
  draw_stack(sig_file, ttbar_file, ttHbb_file, diboson_file, drellyan_file, "phi_LE",  "Lead #phi_{e}", "#phi (Rad)", outputFolder);
  draw_stack(sig_file, ttbar_file, ttHbb_file, diboson_file, drellyan_file, "phi_LM",  "Lead #phi_{#mu}", "#phi (Rad)", outputFolder);
  draw_stack(sig_file, ttbar_file, ttHbb_file, diboson_file, drellyan_file, "phi_SLE",  "Sublead #phi_{e}", "#phi (Rad)", outputFolder);
  draw_stack(sig_file, ttbar_file, ttHbb_file, diboson_file, drellyan_file, "phi_SLM",  "Sublead #phi_{#mu}", "#phi (Rad)", outputFolder);
  draw_stack(sig_file, ttbar_file, ttHbb_file, diboson_file, drellyan_file, "dphi_EE", "#Delta#phi_{ee}", "#Delta#phi (Rad)", outputFolder);
  draw_stack(sig_file, ttbar_file, ttHbb_file, diboson_file, drellyan_file, "dphi_MM", "#Delta#phi_{#mu#mu}", "#Delta#phi (Rad)", outputFolder);
  draw_stack(sig_file, ttbar_file, ttHbb_file, diboson_file, drellyan_file, "phi_LB", "Lead #phi_{b}", "#phi (Rad)", outputFolder);
  draw_stack(sig_file, ttbar_file, ttHbb_file, diboson_file, drellyan_file, "phi_SLB", "Sublead #phi_{b}", "#phi (Rad)", outputFolder);
  draw_stack(sig_file, ttbar_file, ttHbb_file, diboson_file, drellyan_file, "phi_ZH", "#phi_{ZH}", "#phi (Rad)", outputFolder);
  draw_stack(sig_file, ttbar_file, ttHbb_file, diboson_file, drellyan_file, "phi_bb_reco", "Reco #phi_{bb}", "#phi (Rad)", outputFolder);
  draw_stack(sig_file, ttbar_file, ttHbb_file, diboson_file, drellyan_file, "phi_ll_reco", "Reco #phi_{ll}", "#phi (Rad)", outputFolder);
  draw_stack(sig_file, ttbar_file, ttHbb_file, diboson_file, drellyan_file, "phi_LE_reco",  "Lead Reco #phi_{e}", "#phi (Rad)", outputFolder);
  draw_stack(sig_file, ttbar_file, ttHbb_file, diboson_file, drellyan_file, "phi_LM_reco",  "Lead Reco #phi_{#mu}", "#phi (Rad)", outputFolder);
  draw_stack(sig_file, ttbar_file, ttHbb_file, diboson_file, drellyan_file, "phi_SLE_reco",  "Sublead Reco #phi_{e}", "#phi (Rad)", outputFolder);
  draw_stack(sig_file, ttbar_file, ttHbb_file, diboson_file, drellyan_file, "phi_SLM_reco",  "Sublead Reco #phi_{#mu}", "#phi (Rad)", outputFolder);
  draw_stack(sig_file, ttbar_file, ttHbb_file, diboson_file, drellyan_file, "dphi_EE_reco", "Reco #Delta#phi_{ee}", "#Delta#phi (Rad)", outputFolder);
  draw_stack(sig_file, ttbar_file, ttHbb_file, diboson_file, drellyan_file, "dphi_MM_reco", "Reco #Delta#phi_{#mu#mu}", "#Delta#phi (Rad)", outputFolder);
  draw_stack(sig_file, ttbar_file, ttHbb_file, diboson_file, drellyan_file, "phi_LB_reco", "Lead Reco b Jet #phi", "#phi (Rad)", outputFolder);
  draw_stack(sig_file, ttbar_file, ttHbb_file, diboson_file, drellyan_file, "phi_SLB_reco", "Sublead Reco b Jet #phi", "#phi (Rad)", outputFolder);
  draw_stack(sig_file, ttbar_file, ttHbb_file, diboson_file, drellyan_file, "phi_bbll_reco", "Reco #phi_{bbll}", "#phi (Rad)", outputFolder);
  draw_stack(sig_file, ttbar_file, ttHbb_file, diboson_file, drellyan_file, "phi_bb_particle", "#phi_{bb}", "#phi (Rad)", outputFolder);
  draw_stack(sig_file, ttbar_file, ttHbb_file, diboson_file, drellyan_file, "phi_ll_particle", "#phi_{ll}", "#phi (Rad)", outputFolder);
  draw_stack(sig_file, ttbar_file, ttHbb_file, diboson_file, drellyan_file, "phi_LE_particle",  "Lead #phi_{e}", "#phi (Rad)", outputFolder);
  draw_stack(sig_file, ttbar_file, ttHbb_file, diboson_file, drellyan_file, "phi_LM_particle",  "Lead #phi_{#mu}", "#phi (Rad)", outputFolder);
  draw_stack(sig_file, ttbar_file, ttHbb_file, diboson_file, drellyan_file, "phi_SLE_particle",  "Sublead #phi_{e}", "#phi (Rad)", outputFolder);
  draw_stack(sig_file, ttbar_file, ttHbb_file, diboson_file, drellyan_file, "phi_SLM_particle",  "Sublead #phi_{#mu}", "#phi (Rad)", outputFolder);
  draw_stack(sig_file, ttbar_file, ttHbb_file, diboson_file, drellyan_file, "dphi_EE_particle", "#Delta#phi_{ee}", "#Delta#phi (Rad)", outputFolder);
  draw_stack(sig_file, ttbar_file, ttHbb_file, diboson_file, drellyan_file, "dphi_MM_particle", "#Delta#phi_{#mu#mu}", "#Delta#phi (Rad)", outputFolder);
  draw_stack(sig_file, ttbar_file, ttHbb_file, diboson_file, drellyan_file, "phi_LB_particle", "Lead b Jet #phi", "#phi (Rad)", outputFolder);
  draw_stack(sig_file, ttbar_file, ttHbb_file, diboson_file, drellyan_file, "phi_SLB_particle", "Sublead b Jet #phi", "#phi (Rad)", outputFolder);
  draw_stack(sig_file, ttbar_file, ttHbb_file, diboson_file, drellyan_file, "phi_bbll_particle", "#phi_{bbll}", "#phi (Rad)", outputFolder);
  // signal comparison histograms
  // mass comparison
  draw_hist2(sig_file, "mass_bb_Comp", "m_{bb}", "Particle Level mass (GeV)", "Reco Level mass (GeV)", outputFolder, "signal2D");
  draw_hist2(sig_file, "mass_ll_Comp", "m_{ll}", "Particle Level mass (GeV)", "Reco Level mass (GeV)", outputFolder, "signal2D");
  draw_hist2(sig_file, "mass_ZH_Comp", "m_{T}^{ZH}", "Parton Level mass (GeV)", "Reco Level mass (GeV)", outputFolder, "signal2D");
  draw_hist2(sig_file, "mass_ZH_Comp_particle", "m_{T}^{ZH}", "Particle Level mass (GeV)", "Reco Level mass (GeV)", outputFolder, "signal2D");
  // pt comparison
  draw_hist2(sig_file, "pt_H_Comp", "p_{T}^{H}", "Parton Level P_{T} (GeV)", "Reco Level P_{T} (GeV)", outputFolder, "signal2D");
  draw_hist2(sig_file, "pt_Z_Comp", "p_{T}^{Z}", "Parton Level P_{T} (GeV)", "Reco Level P_{T} (GeV)", outputFolder, "signal2D");
  draw_hist2(sig_file, "pt_LE_Comp", "Lead p_{T}^{e}", "Parton Level P_{T} (GeV)", "Reco Level P_{T} (GeV)", outputFolder, "signal2D");
  draw_hist2(sig_file, "pt_LM_Comp", "Lead p_{T}^{#mu}", "Parton Level P_{T} (GeV)", "Reco Level P_{T} (GeV)", outputFolder, "signal2D");
  draw_hist2(sig_file, "pt_SLE_Comp", "Sublead p_{T}^{e}", "Parton Level P_{T} (GeV)", "Reco Level P_{T} (GeV)", outputFolder, "signal2D");
  draw_hist2(sig_file, "pt_SLM_Comp", "Sublead p_{T}^{#mu}", "Parton Level P_{T} (GeV)", "Reco Level P_{T} (GeV)", outputFolder, "signal2D");
  draw_hist2(sig_file, "pt_LB_Comp", "Lead p_{T}^{b}", "Parton Level P_{T} (GeV)", "Reco Level P_{T} (GeV)", outputFolder, "signal2D");
  draw_hist2(sig_file, "pt_SLB_Comp", "Sublead p_{T}^{b}", "Parton Level P_{T} (GeV)", "Reco Level P_{T} (GeV)", outputFolder, "signal2D");
  draw_hist2(sig_file, "pt_ZH_Comp", "Sublead p_{T}^{ZH}", "Parton Level P_{T} (GeV)", "Reco Level P_{T} (GeV)", outputFolder, "signal2D");
  draw_hist2(sig_file, "pt_H_Comp_particle", "p_{T}^{H}", "Particle Level P_{T} (GeV)", "Reco Level P_{T} (GeV)", outputFolder, "signal2D");
  draw_hist2(sig_file, "pt_Z_Comp_particle", "p_{T}^{Z}", "Particle Level P_{T} (GeV)", "Reco Level P_{T} (GeV)", outputFolder, "signal2D");
  draw_hist2(sig_file, "pt_LE_Comp_particle", "Lead p_{T}^{e}", "Particle Level P_{T} (GeV)", "Reco Level P_{T} (GeV)", outputFolder, "signal2D");
  draw_hist2(sig_file, "pt_LM_Comp_particle", "Lead p_{T}^{#mu}", "Particle Level P_{T} (GeV)", "Reco Level P_{T} (GeV)", outputFolder, "signal2D");
  draw_hist2(sig_file, "pt_SLE_Comp_particle", "Sublead p_{T}^{e}", "Particle Level P_{T} (GeV)", "Reco Level P_{T} (GeV)", outputFolder, "signal2D");
  draw_hist2(sig_file, "pt_SLM_Comp_particle", "Sublead p_{T}^{#mu}", "Particle Level P_{T} (GeV)", "Reco Level P_{T} (GeV)", outputFolder, "signal2D");
  draw_hist2(sig_file, "pt_LB_Comp_particle", "Lead p_{T}^{b}", "Particle Level P_{T} (GeV)", "Reco Level P_{T} (GeV)", outputFolder, "signal2D");
  draw_hist2(sig_file, "pt_SLB_Comp_particle", "Sublead p_{T}^{b}", "Particle Level P_{T} (GeV)", "Reco Level P_{T} (GeV)", outputFolder, "signal2D");
  draw_hist2(sig_file, "pt_ZH_Comp_particle", "Sublead p_{T}^{ZH}", "Particle Level P_{T} (GeV)", "Reco Level P_{T} (GeV)", outputFolder, "signal2D");
  // eta comparison
  draw_hist2(sig_file, "eta_H_Comp", "#eta_{H}", "Parton Level #eta", "Reco Level #eta", outputFolder, "signal2D");
  draw_hist2(sig_file, "eta_Z_Comp", "#eta_{Z}", "Parton Level #eta", "Reco Level #eta", outputFolder, "signal2D");
  draw_hist2(sig_file, "eta_LE_Comp", "Lead #eta_{e}","Parton Level #eta", "Reco Level #eta", outputFolder, "signal2D");
  draw_hist2(sig_file, "eta_LM_Comp", "Lead #eta_{#mu}","Parton Level #eta", "Reco Level #eta", outputFolder, "signal2D");
  draw_hist2(sig_file, "eta_SLE_Comp", "Sublead #eta_{e}","Parton Level #eta", "Reco Level #eta", outputFolder, "signal2D");
  draw_hist2(sig_file, "eta_SLM_Comp", "Sublead #eta_{#mu}","Parton Level #eta", "Reco Level #eta", outputFolder, "signal2D");
  draw_hist2(sig_file, "eta_LB_Comp", "Lead #eta_{b}","Parton Level #eta", "Reco Level #eta", outputFolder, "signal2D");
  draw_hist2(sig_file, "eta_SLB_Comp", "Sublead #eta_{b}","Parton Level #eta", "Reco Level #eta", outputFolder, "signal2D");
  draw_hist2(sig_file, "eta_ZH_Comp", "Sublead #eta_{ZH}", "Parton Level #eta", "Reco Level #eta", outputFolder, "signal2D");
  draw_hist2(sig_file, "eta_H_Comp_particle", "#eta_{H}", "Particle Level #eta", "Reco Level #eta", outputFolder, "signal2D");
  draw_hist2(sig_file, "eta_Z_Comp_particle", "#eta_{Z}", "Particle Level #eta", "Reco Level #eta", outputFolder, "signal2D");
  draw_hist2(sig_file, "eta_LE_Comp_particle", "Lead #eta_{e}","Particle Level #eta", "Reco Level #eta", outputFolder, "signal2D");
  draw_hist2(sig_file, "eta_LM_Comp_particle", "Lead #eta_{#mu}","Particle Level #eta", "Reco Level #eta", outputFolder, "signal2D");
  draw_hist2(sig_file, "eta_SLE_Comp_particle", "Sublead #eta_{e}","Particle Level #eta", "Reco Level #eta", outputFolder, "signal2D");
  draw_hist2(sig_file, "eta_SLM_Comp_particle", "Sublead #eta_{#mu}","Particle Level #eta", "Reco Level #eta", outputFolder, "signal2D");
  draw_hist2(sig_file, "eta_LB_Comp_particle", "Lead #eta_{b}","Particle Level #eta", "Reco Level #eta", outputFolder, "signal2D");
  draw_hist2(sig_file, "eta_SLB_Comp_particle", "Sublead #eta_{b}","Particle Level #eta", "Reco Level #eta", outputFolder, "signal2D");
  draw_hist2(sig_file, "eta_ZH_Comp_particle", "Sublead #eta_{ZH}", "Particle Level #eta", "Reco Level #eta", outputFolder, "signal2D");
  draw_hist2(sig_file, "deta_EE_Comp", "#Delta#eta_{ee}","Parton Level #Delta#eta", "Reco Level #Delta#eta", outputFolder, "signal2D");
  draw_hist2(sig_file, "deta_MM_Comp", "#Delta#eta_{#mu#mu}","Parton Level #Delta#eta", "Reco Level #Delta#eta", outputFolder, "signal2D");
  draw_hist2(sig_file, "deta_EE_Comp_particle", "#Delta#eta_{ee}","Particle Level #Delta#eta", "Reco Level #Delta#eta", outputFolder, "signal2D");
  draw_hist2(sig_file, "deta_MM_Comp_particle", "#Delta#eta_{#mu#mu}","Particle Level #Delta#eta", "Reco Level #Delta#eta", outputFolder, "signal2D");
  // cos theta comparison
  draw_hist2(sig_file, "costheta_LE_Comp", "Lead cos#theta_{e}","Parton Level cos#theta", "Reco Level cos#theta", outputFolder, "signal2D");
  draw_hist2(sig_file, "costheta_LM_Comp", "Lead cos#theta_{#mu}","Parton Level cos#theta", "Reco Level cos#theta", outputFolder, "signal2D");
  draw_hist2(sig_file, "costheta_SLE_Comp", "Sublead cos#theta_{e}","Parton Level cos#theta", "Reco Level cos#theta", outputFolder, "signal2D");
  draw_hist2(sig_file, "costheta_SLM_Comp", "Sublead cos#theta_{#mu}","Parton Level cos#theta", "Reco Level cos#theta", outputFolder, "signal2D");
  draw_hist2(sig_file, "costheta_LB_Comp", "Lead cos#theta_{b}","Parton Level cos#theta", "Reco Level cos#theta", outputFolder, "signal2D");
  draw_hist2(sig_file, "costheta_SLB_Comp", "Sublead cos#theta_{b}","Parton Level cos#theta", "Reco Level cos#theta", outputFolder, "signal2D");
  draw_hist2(sig_file, "costheta_LE_Comp_particle", "Lead cos#theta_{e}","Particle Level cos#theta", "Reco Level cos#theta", outputFolder, "signal2D");
  draw_hist2(sig_file, "costheta_LM_Comp_particle", "Lead cos#theta_{#mu}","Particle Level cos#theta", "Reco Level cos#theta", outputFolder, "signal2D");
  draw_hist2(sig_file, "costheta_SLE_Comp_particle", "Sublead cos#theta_{e}","Particle Level cos#theta", "Reco Level cos#theta", outputFolder, "signal2D");
  draw_hist2(sig_file, "costheta_SLM_Comp_particle", "Sublead cos#theta_{#mu}","Particle Level cos#theta", "Reco Level cos#theta", outputFolder, "signal2D");
  draw_hist2(sig_file, "costheta_LB_Comp_particle", "Lead cos#theta_{b}","Particle Level cos#theta", "Reco Level cos#theta", outputFolder, "signal2D");
  draw_hist2(sig_file, "costheta_SLB_Comp_particle", "Sublead cos#theta_{b}","Particle Level cos#theta", "Reco Level cos#theta", outputFolder, "signal2D");
  // phi comparison
  draw_hist2(sig_file, "phi_H_Comp", "#phi_{H}", "Parton Level #phi (Rad)", "Reco Level #phi (Rad)", outputFolder, "signal2D");
  draw_hist2(sig_file, "phi_Z_Comp", "#phi_{Z}", "Parton Level #phi (Rad)", "Reco Level #phi (Rad)", outputFolder, "signal2D");
  draw_hist2(sig_file, "phi_LE_Comp", "Lead #phi_{e}","Parton Level #phi (Rad)", "Reco Level #phi (Rad)", outputFolder, "signal2D");
  draw_hist2(sig_file, "phi_LM_Comp", "Lead #phi_{#mu}","Parton Level #phi (Rad)", "Reco Level #phi (Rad)", outputFolder, "signal2D");
  draw_hist2(sig_file, "phi_SLE_Comp", "Sublead #phi_{e}","Parton Level #phi (Rad)", "Reco Level #phi (Rad)", outputFolder, "signal2D");
  draw_hist2(sig_file, "phi_SLM_Comp", "Sublead #phi_{#mu}","Parton Level #phi (Rad)", "Reco Level #phi (Rad)", outputFolder, "signal2D");
  draw_hist2(sig_file, "dphi_EE_Comp", "#Delta#phi_{ee}","Parton Level #Delta#phi (Rad)", "Reco Level #Delta#phi (Rad)", outputFolder, "signal2D");
  draw_hist2(sig_file, "dphi_MM_Comp", "#Delta#phi_{#mu#mu}","Parton Level #Delta#phi (Rad)", "Reco Level #Delta#phi (Rad)", outputFolder, "signal2D");
  draw_hist2(sig_file, "dphi_EE_Comp_particle", "#Delta#phi_{ee}","Particle Level #Delta#phi (Rad)", "Reco Level #Delta#phi (Rad)", outputFolder, "signal2D");
  draw_hist2(sig_file, "dphi_MM_Comp_particle", "#Delta#phi_{#mu#mu}","Particle Level #Delta#phi (Rad)", "Reco Level #Delta#phi (Rad)", outputFolder, "signal2D");
  draw_hist2(sig_file, "phi_LB_Comp", "Lead #phi_{b}","Parton Level #phi (Rad)", "Reco Level #phi (Rad)", outputFolder, "signal2D");
  draw_hist2(sig_file, "phi_SLB_Comp", "Sublead #phi_{b}","Parton Level #phi (Rad)", "Reco Level #phi (Rad)", outputFolder, "signal2D");
  draw_hist2(sig_file, "phi_ZH_Comp", "Sublead #phi_{ZH}", "Parton Level #phi (Rad)", "Reco Level #phi (Rad)", outputFolder, "signal2D");
  // bkg comparison histograms
  // mass comparison
  draw_hist2(bkg_file, "mass_bb_Comp", "m_{bb}", "Particle Level mass (GeV)", "Reco Level mass (GeV)", outputFolder, "bkg2D");
  draw_hist2(bkg_file, "mass_ll_Comp", "m_{ll}", "Particle Level mass (GeV)", "Reco Level mass (GeV)", outputFolder, "bkg2D");
  draw_hist2(bkg_file, "mass_ZH_Comp", "m_{T}^{ZH}", "Parton Level mass (GeV)", "Reco Level mass (GeV)", outputFolder, "bkg2D");
  draw_hist2(bkg_file, "mass_ZH_Comp_particle", "m_{T}^{ZH}", "Particle Level mass (GeV)", "Reco Level mass (GeV)", outputFolder, "bkg2D");
  // pt comparison
  draw_hist2(bkg_file, "pt_H_Comp", "p_{T}^{H}", "Parton Level P_{T} (GeV)", "Reco Level P_{T} (GeV)", outputFolder, "bkg2D");
  draw_hist2(bkg_file, "pt_Z_Comp", "p_{T}^{Z}", "Parton Level P_{T} (GeV)", "Reco Level P_{T} (GeV)", outputFolder, "bkg2D");
  draw_hist2(bkg_file, "pt_LE_Comp", "Lead p_{T}^{e}", "Parton Level P_{T} (GeV)", "Reco Level P_{T} (GeV)", outputFolder, "bkg2D");
  draw_hist2(bkg_file, "pt_LM_Comp", "Lead p_{T}^{#mu}", "Parton Level P_{T} (GeV)", "Reco Level P_{T} (GeV)", outputFolder, "bkg2D");
  draw_hist2(bkg_file, "pt_SLE_Comp", "Sublead p_{T}^{e}", "Parton Level P_{T} (GeV)", "Reco Level P_{T} (GeV)", outputFolder, "bkg2D");
  draw_hist2(bkg_file, "pt_SLM_Comp", "Sublead p_{T}^{#mu}", "Parton Level P_{T} (GeV)", "Reco Level P_{T} (GeV)", outputFolder, "bkg2D");
  draw_hist2(bkg_file, "pt_LB_Comp", "Lead p_{T}^{b}", "Parton Level P_{T} (GeV)", "Reco Level P_{T} (GeV)", outputFolder, "bkg2D");
  draw_hist2(bkg_file, "pt_SLB_Comp", "Sublead p_{T}^{b}", "Parton Level P_{T} (GeV)", "Reco Level P_{T} (GeV)", outputFolder, "bkg2D");
  draw_hist2(bkg_file, "pt_ZH_Comp", "Sublead p_{T}^{ZH}", "Parton Level P_{T} (GeV)", "Reco Level P_{T} (GeV)", outputFolder, "bkg2D");
  draw_hist2(bkg_file, "pt_H_Comp_particle", "p_{T}^{H}", "Particle Level P_{T} (GeV)", "Reco Level P_{T} (GeV)", outputFolder, "bkg2D");
  draw_hist2(bkg_file, "pt_Z_Comp_particle", "p_{T}^{Z}", "Particle Level P_{T} (GeV)", "Reco Level P_{T} (GeV)", outputFolder, "bkg2D");
  draw_hist2(bkg_file, "pt_LE_Comp_particle", "Lead p_{T}^{e}", "Particle Level P_{T} (GeV)", "Reco Level P_{T} (GeV)", outputFolder, "bkg2D");
  draw_hist2(bkg_file, "pt_LM_Comp_particle", "Lead p_{T}^{#mu}", "Particle Level P_{T} (GeV)", "Reco Level P_{T} (GeV)", outputFolder, "bkg2D");
  draw_hist2(bkg_file, "pt_SLE_Comp_particle", "Sublead p_{T}^{e}", "Particle Level P_{T} (GeV)", "Reco Level P_{T} (GeV)", outputFolder, "bkg2D");
  draw_hist2(bkg_file, "pt_SLM_Comp_particle", "Sublead p_{T}^{#mu}", "Particle Level P_{T} (GeV)", "Reco Level P_{T} (GeV)", outputFolder, "bkg2D");
  draw_hist2(bkg_file, "pt_LB_Comp_particle", "Lead p_{T}^{b}", "Particle Level P_{T} (GeV)", "Reco Level P_{T} (GeV)", outputFolder, "bkg2D");
  draw_hist2(bkg_file, "pt_SLB_Comp_particle", "Sublead p_{T}^{b}", "Particle Level P_{T} (GeV)", "Reco Level P_{T} (GeV)", outputFolder, "bkg2D");
  draw_hist2(bkg_file, "pt_ZH_Comp_particle", "Sublead p_{T}^{ZH}", "Particle Level P_{T} (GeV)", "Reco Level P_{T} (GeV)", outputFolder, "bkg2D");
  // eta comparison
  draw_hist2(bkg_file, "eta_H_Comp", "#eta_{H}", "Parton Level #eta", "Reco Level #eta", outputFolder, "bkg2D");
  draw_hist2(bkg_file, "eta_Z_Comp", "#eta_{Z}", "Parton Level #eta", "Reco Level #eta", outputFolder, "bkg2D");
  draw_hist2(bkg_file, "eta_LE_Comp", "Lead #eta_{e}","Parton Level #eta", "Reco Level #eta", outputFolder, "bkg2D");
  draw_hist2(bkg_file, "eta_LM_Comp", "Lead #eta_{#mu}","Parton Level #eta", "Reco Level #eta", outputFolder, "bkg2D");
  draw_hist2(bkg_file, "eta_SLE_Comp", "Sublead #eta_{e}","Parton Level #eta", "Reco Level #eta", outputFolder, "bkg2D");
  draw_hist2(bkg_file, "eta_SLM_Comp", "Sublead #eta_{#mu}","Parton Level #eta", "Reco Level #eta", outputFolder, "bkg2D");
  draw_hist2(bkg_file, "eta_LB_Comp", "Lead #eta_{b}","Parton Level #eta", "Reco Level #eta", outputFolder, "bkg2D");
  draw_hist2(bkg_file, "eta_SLB_Comp", "Sublead #eta_{b}","Parton Level #eta", "Reco Level #eta", outputFolder, "bkg2D");
  draw_hist2(bkg_file, "eta_ZH_Comp", "Sublead #eta_{ZH}", "Parton Level #eta", "Reco Level #eta", outputFolder, "bkg2D");
  draw_hist2(bkg_file, "eta_H_Comp_particle", "#eta_{H}", "Particle Level #eta", "Reco Level #eta", outputFolder, "bkg2D");
  draw_hist2(bkg_file, "eta_Z_Comp_particle", "#eta_{Z}", "Particle Level #eta", "Reco Level #eta", outputFolder, "bkg2D");
  draw_hist2(bkg_file, "eta_LE_Comp_particle", "Lead #eta_{e}","Particle Level #eta", "Reco Level #eta", outputFolder, "bkg2D");
  draw_hist2(bkg_file, "eta_LM_Comp_particle", "Lead #eta_{#mu}","Particle Level #eta", "Reco Level #eta", outputFolder, "bkg2D");
  draw_hist2(bkg_file, "eta_SLE_Comp_particle", "Sublead #eta_{e}","Particle Level #eta", "Reco Level #eta", outputFolder, "bkg2D");
  draw_hist2(bkg_file, "eta_SLM_Comp_particle", "Sublead #eta_{#mu}","Particle Level #eta", "Reco Level #eta", outputFolder, "bkg2D");
  draw_hist2(bkg_file, "eta_LB_Comp_particle", "Lead #eta_{b}","Particle Level #eta", "Reco Level #eta", outputFolder, "bkg2D");
  draw_hist2(bkg_file, "eta_SLB_Comp_particle", "Sublead #eta_{b}","Particle Level #eta", "Reco Level #eta", outputFolder, "bkg2D");
  draw_hist2(bkg_file, "eta_ZH_Comp_particle", "Sublead #eta_{ZH}", "Particle Level #eta", "Reco Level #eta", outputFolder, "bkg2D");
  draw_hist2(bkg_file, "deta_EE_Comp", "#Delta#eta_{ee}","Parton Level #Delta#eta", "Reco Level #Delta#eta", outputFolder, "bkg2D");
  draw_hist2(bkg_file, "deta_MM_Comp", "#Delta#eta_{#mu#mu}","Parton Level #Delta#eta", "Reco Level #Delta#eta", outputFolder, "bkg2D");
  draw_hist2(bkg_file, "deta_EE_Comp_particle", "#Delta#eta_{ee}","Particle Level #Delta#eta", "Reco Level #Delta#eta", outputFolder, "bkg2D");
  draw_hist2(bkg_file, "deta_MM_Comp_particle", "#Delta#eta_{#mu#mu}","Particle Level #Delta#eta", "Reco Level #Delta#eta", outputFolder, "bkg2D");
  // cos theta comparison
  draw_hist2(bkg_file, "costheta_LE_Comp", "Lead cos#theta_{e}","Parton Level cos#theta", "Reco Level cos#theta", outputFolder, "bkg2D");
  draw_hist2(bkg_file, "costheta_LM_Comp", "Lead cos#theta_{#mu}","Parton Level cos#theta", "Reco Level cos#theta", outputFolder, "bkg2D");
  draw_hist2(bkg_file, "costheta_SLE_Comp", "Sublead cos#theta_{e}","Parton Level cos#theta", "Reco Level cos#theta", outputFolder, "bkg2D");
  draw_hist2(bkg_file, "costheta_SLM_Comp", "Sublead cos#theta_{#mu}","Parton Level cos#theta", "Reco Level cos#theta", outputFolder, "bkg2D");
  draw_hist2(bkg_file, "costheta_LB_Comp", "Lead cos#theta_{b}","Parton Level cos#theta", "Reco Level cos#theta", outputFolder, "bkg2D");
  draw_hist2(bkg_file, "costheta_SLB_Comp", "Sublead cos#theta_{b}","Parton Level cos#theta", "Reco Level cos#theta", outputFolder, "bkg2D");
  draw_hist2(bkg_file, "costheta_LE_Comp_particle", "Lead cos#theta_{e}","Particle Level cos#theta", "Reco Level cos#theta", outputFolder, "bkg2D");
  draw_hist2(bkg_file, "costheta_LM_Comp_particle", "Lead cos#theta_{#mu}","Particle Level cos#theta", "Reco Level cos#theta", outputFolder, "bkg2D");
  draw_hist2(bkg_file, "costheta_SLE_Comp_particle", "Sublead cos#theta_{e}","Particle Level cos#theta", "Reco Level cos#theta", outputFolder, "bkg2D");
  draw_hist2(bkg_file, "costheta_SLM_Comp_particle", "Sublead cos#theta_{#mu}","Particle Level cos#theta", "Reco Level cos#theta", outputFolder, "bkg2D");
  draw_hist2(bkg_file, "costheta_LB_Comp_particle", "Lead cos#theta_{b}","Particle Level cos#theta", "Reco Level cos#theta", outputFolder, "bkg2D");
  draw_hist2(bkg_file, "costheta_SLB_Comp_particle", "Sublead cos#theta_{b}","Particle Level cos#theta", "Reco Level cos#theta", outputFolder, "bkg2D");
  // phi comparison
  draw_hist2(bkg_file, "phi_H_Comp", "#phi_{H}", "Parton Level #phi (Rad)", "Reco Level #phi (Rad)", outputFolder, "bkg2D");
  draw_hist2(bkg_file, "phi_Z_Comp", "#phi_{Z}", "Parton Level #phi (Rad)", "Reco Level #phi (Rad)", outputFolder, "bkg2D");
  draw_hist2(bkg_file, "phi_LE_Comp", "Lead #phi_{e}","Parton Level #phi (Rad)", "Reco Level #phi (Rad)", outputFolder, "bkg2D");
  draw_hist2(bkg_file, "phi_LM_Comp", "Lead #phi_{#mu}","Parton Level #phi (Rad)", "Reco Level #phi (Rad)", outputFolder, "bkg2D");
  draw_hist2(bkg_file, "phi_SLE_Comp", "Sublead #phi_{e}","Parton Level #phi (Rad)", "Reco Level #phi (Rad)", outputFolder, "bkg2D");
  draw_hist2(bkg_file, "phi_SLM_Comp", "Sublead #phi_{#mu}","Parton Level #phi (Rad)", "Reco Level #phi (Rad)", outputFolder, "bkg2D");
  draw_hist2(bkg_file, "dphi_EE_Comp", "#Delta#phi_{ee}","Parton Level #Delta#phi (Rad)", "Reco Level #Delta#phi (Rad)", outputFolder, "bkg2D");
  draw_hist2(bkg_file, "dphi_MM_Comp", "#Delta#phi_{#mu#mu}","Parton Level #Delta#phi (Rad)", "Reco Level #Delta#phi (Rad)", outputFolder, "bkg2D");
  draw_hist2(bkg_file, "dphi_EE_Comp_particle", "#Delta#phi_{ee}","Particle Level #Delta#phi (Rad)", "Reco Level #Delta#phi (Rad)", outputFolder, "bkg2D");
  draw_hist2(bkg_file, "dphi_MM_Comp_particle", "#Delta#phi_{#mu#mu}","Particle Level #Delta#phi (Rad)", "Reco Level #Delta#phi (Rad)", outputFolder, "bkg2D");
  draw_hist2(bkg_file, "phi_LB_Comp", "Lead #phi_{b}","Parton Level #phi (Rad)", "Reco Level #phi (Rad)", outputFolder, "bkg2D");
  draw_hist2(bkg_file, "phi_SLB_Comp", "Sublead #phi_{b}","Parton Level #phi (Rad)", "Reco Level #phi (Rad)", outputFolder, "bkg2D");
  draw_hist2(bkg_file, "phi_ZH_Comp", "Sublead #phi_{ZH}", "Parton Level #phi (Rad)", "Reco Level #phi (Rad)", outputFolder, "bkg2D");
  // print signal significance
  TH1 *sig_mbb = dynamic_cast<TH1*>(sig_file->Get("mass_bb_reco"));
  TH1 *bkg_mbb = dynamic_cast<TH1*>(sig_file->Get("mass_bb_reco"));
  TH1 *ttbar_mbb = dynamic_cast<TH1*>(ttbar_file->Get("mass_bb_reco"));
  TH1 *ttHbb_mbb = dynamic_cast<TH1*>(ttHbb_file->Get("mass_bb_reco"));
  TH1 *diboson_mbb = dynamic_cast<TH1*>(diboson_file->Get("mass_bb_reco"));
  TH1 *drellyan_mbb = dynamic_cast<TH1*>(drellyan_file->Get("mass_bb_reco"));

  std::cout<<"ttbar contribution: "<<ttbar_mbb->Integral()<<std::endl;
  std::cout<<"ttHbb contribution: "<<ttHbb_mbb->Integral()<<std::endl;
  std::cout<<"diboson contribution: "<<diboson_mbb->Integral()<<std::endl;
  std::cout<<"drell-yan contribution: "<<drellyan_mbb->Integral()<<std::endl;
  
  double sig_sig = sig_mbb->Integral()/sqrt(bkg_mbb->Integral());
  std::cout<<"Signal Significance: "<<sig_sig<<std::endl;
  sig_file->Close();
  bkg_file->Close();
}