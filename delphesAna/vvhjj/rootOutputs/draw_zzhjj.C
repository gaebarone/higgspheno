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

  // sig_hist->SetFillColor(kCyan);
  ttbar_hist->SetFillColor(kOrange);
  ttHbb_hist->SetFillColor(kViolet);
  diboson_hist->SetFillColor(kRed);
  drellyan_hist->SetFillColor(kYellow);

  // rebin
  
  sig_hist->Rebin(2);
  ttbar_hist->Rebin(2);
  ttHbb_hist->Rebin(2);
  diboson_hist->Rebin(2);
  drellyan_hist->Rebin(2);
  
  /*  
  ttbar_hist->Scale(ttbar_scale);
  ttHbb_hist->Scale(ttHbb_scale);
  drellyan_hist->Scale(drellyan_scale);
  diboson_hist->Scale(diboson_scale);
  */

  // form the histogram stack
  THStack *stack = new THStack(name, title);
  stack->Add(ttbar_hist);
  stack->Add(drellyan_hist);
  stack->Add(diboson_hist);
  stack->Add(ttHbb_hist);
  //stack->Add(sig_hist);

  // make a legend
  TLegend *legend = new TLegend(0.725, 0.725, 0.875, 0.875);
  //legend->AddEntry(sig_hist, "  signal", "f");
  legend->AddEntry(drellyan_hist, "  drellyan", "f");
  legend->AddEntry(diboson_hist, "  diboson", "f");
  legend->AddEntry(ttbar_hist, "  ttbar", "f");
  legend->AddEntry(ttHbb_hist, "  ttHbb", "f");
  
  // make a canvas and draw on it
  TCanvas *c = new TCanvas(name, title, 1500, 1200);
  c->cd();
  stack->Draw("hist e");
  stack->GetXaxis()->SetTitle(axistitle);
  legend->Draw();

  TH1F *summed_hist = (TH1F*)stack->GetHistogram()->Clone();
  summed_hist->SetStats(kTRUE); // Enable stat box
  c->Update(); // Update canvas to make stat box appear

  // sig on top
  sig_hist->Scale(10000);
  sig_hist->SetLineColor(kBlack);
  legend->AddEntry(sig_hist, "sig x 10000", "l");
  sig_hist->SetLineWidth(2);
  sig_hist->Draw("hist same");
  
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

//void draw_zzhjj(const char *sig_filename = "../../../outputs/histograms/signal.root", const char *bkg_filename = "../../../outputs/histograms/all_bkg.root", const char *ttbar_filename = "../../../outputs/histograms/ttbar.root", const char *ttHbb_filename = "../../../outputs/histograms/ttHbb.root", const char *diboson_filename = "../../../outputs/histograms/diboson.root", const char *drellyan_filename = "../../../outputs/histograms/drellyan.root", const char *outputFolder = "../../../outputs/plots") {

void draw_zzhjj(const char *sig_filename, const char *bkg_filename, const char *ttbar_filename, const char *ttHbb_filename, const char *diboson_filename, const char *drellyan_filename, const char *outputFolder){

  // set scales according to total number of events
  /*sig_scale = 1.0/get_total_events("zh_zll_hbb_012j");
  ttbar_scale = 1.0/get_total_events("ttbar012j");
  ttHbb_scale = 1.0/get_total_events("ttHbb");
  drellyan_scale = 1.0/get_total_events("zll_123j");
  diboson_scale = 1.0/(get_total_events("wwjj_j")+get_total_events("wzjj_j")+get_total_events("wz_wjj_123j")+get_total_events("zzjj_j")+get_total_events("zz_zjj_123j"));
  */

  // open files
  TFile *sig_file = TFile::Open(sig_filename, "READ");
  TFile *bkg_file = TFile::Open(bkg_filename, "READ");
  TFile *ttbar_file = TFile::Open(ttbar_filename, "READ");
  TFile *ttHbb_file = TFile::Open(ttHbb_filename, "READ");
  TFile *diboson_file = TFile::Open(diboson_filename, "READ");
  TFile *drellyan_file = TFile::Open(drellyan_filename, "READ");
  // draw histograms
  gROOT->SetBatch(kTRUE);

// higgs - reco
  draw_stack(sig_file, ttbar_file, ttHbb_file, diboson_file, drellyan_file, "hbb_pT_reco", "Reco Level p^{T}_{hbb}", "p_{T}", outputFolder);
  draw_stack(sig_file, ttbar_file, ttHbb_file, diboson_file, drellyan_file, "hbb_m_reco", "Reco Level m_{hbb}", "mass (GeV)", outputFolder);
  draw_stack(sig_file, ttbar_file, ttHbb_file, diboson_file, drellyan_file, "bb_#Delta#phi_reco", "Reco Level #Delta#phi_{bb}", "#phi", outputFolder);
  draw_stack(sig_file, ttbar_file, ttHbb_file, diboson_file, drellyan_file, "bb_#Delta#eta_reco", "Reco Level #Delta#eta_{bb}", "#eta", outputFolder);
  draw_stack(sig_file, ttbar_file, ttHbb_file, diboson_file, drellyan_file, "bb_#DeltaR_reco", "Reco Level #DeltaR_{bb}", "R", outputFolder);  
// higgs - particle
  draw_stack(sig_file, ttbar_file, ttHbb_file, diboson_file, drellyan_file, "hbb_pT_particle", "Particle Level p^{T}_{hbb}", "p_{T}", outputFolder);
  draw_stack(sig_file, ttbar_file, ttHbb_file, diboson_file, drellyan_file, "hbb_m_particle", "Particle Level m_{hbb}", "mass (GeV)", outputFolder);
  draw_stack(sig_file, ttbar_file, ttHbb_file, diboson_file, drellyan_file, "bb_#Delta#phi_particle", "Particle Level #Delta#phi_{bb}", "#phi", outputFolder);
  draw_stack(sig_file, ttbar_file, ttHbb_file, diboson_file, drellyan_file, "bb_#Delta#eta_particle", "Particle Level #Delta#eta_{bb}", "#eta", outputFolder);
  draw_stack(sig_file, ttbar_file, ttHbb_file, diboson_file, drellyan_file, "bb_#DeltaR_particle", "Particle Level #DeltaR_{bb}", "R", outputFolder);
// higgs - parton
  draw_stack(sig_file, ttbar_file, ttHbb_file, diboson_file, drellyan_file, "hbb_pT_parton", "Parton Level p^{T}_{hbb}", "p_{T}", outputFolder);
  draw_stack(sig_file, ttbar_file, ttHbb_file, diboson_file, drellyan_file, "hbb_m_parton", "Parton Level m_{hbb}", "mass (GeV)", outputFolder);
  draw_stack(sig_file, ttbar_file, ttHbb_file, diboson_file, drellyan_file, "bb_#Delta#phi_parton", "Parton Level #Delta#phi_{bb}", "#phi", outputFolder);
  draw_stack(sig_file, ttbar_file, ttHbb_file, diboson_file, drellyan_file, "bb_#Delta#eta_parton", "Parton Level #Delta#eta_{bb}", "#eta", outputFolder);
  draw_stack(sig_file, ttbar_file, ttHbb_file, diboson_file, drellyan_file, "bb_#DeltaR_parton", "Parton Level #DeltaR_{bb}", "R", outputFolder);
// vbfj - reco
  draw_stack(sig_file, ttbar_file, ttHbb_file, diboson_file, drellyan_file, "jj_pT_reco", "Reco Level p^{T}_{jj}", "p_{T}", outputFolder);
  draw_stack(sig_file, ttbar_file, ttHbb_file, diboson_file, drellyan_file, "jj_#Delta#phi_reco", "Reco Level #Delta#phi_{jj}", "#phi", outputFolder);
  draw_stack(sig_file, ttbar_file, ttHbb_file, diboson_file, drellyan_file, "jj_#Delta#eta_reco", "Reco Level #Delta#eta_{jj}", "#eta", outputFolder);
  draw_stack(sig_file, ttbar_file, ttHbb_file, diboson_file, drellyan_file, "jj_#DeltaR_reco", "Reco Level #DeltaR_{jj}", "R", outputFolder);
// vbfj - particle
  draw_stack(sig_file, ttbar_file, ttHbb_file, diboson_file, drellyan_file, "jj_pT_particle", "Particle Level p^{T}_{jj}", "p_{T}", outputFolder);
  draw_stack(sig_file, ttbar_file, ttHbb_file, diboson_file, drellyan_file, "jj_#DeltaR_particle", "Particle Level #DeltaR_{jj}", "R", outputFolder);
  draw_stack(sig_file, ttbar_file, ttHbb_file, diboson_file, drellyan_file, "jj_#Delta#eta_particle", "Particle Level #Delta#eta_{jj}", "#eta", outputFolder);
  draw_stack(sig_file, ttbar_file, ttHbb_file, diboson_file, drellyan_file, "jj_#DeltaR_particle", "Particle Level #DeltaR_{jj}", "R", outputFolder);

// cutFlow - reco
  draw_stack(sig_file, ttbar_file, ttHbb_file, diboson_file, drellyan_file,"hSel_reco", "Reco Level Selection", " ", outputFolder);
  draw_stack(sig_file, ttbar_file, ttHbb_file, diboson_file, drellyan_file,"hEff_reco", "Reco Level Selection Efficiency", " ", outputFolder);
// cutFlow - particle
  draw_stack(sig_file, ttbar_file, ttHbb_file, diboson_file, drellyan_file,"hSel_particle", "Particle Level Selection", " ", outputFolder);
  draw_stack(sig_file, ttbar_file, ttHbb_file, diboson_file, drellyan_file,"hEff_particle", "Particle Level Selection Efficiency", " ", outputFolder);

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


  // print signal significance
  std::cout<<"Reco Level Contrubutions:"<<std::endl;
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

  std::cout<<"Integral before cuts with Lumi=1:"<<std::endl;
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
