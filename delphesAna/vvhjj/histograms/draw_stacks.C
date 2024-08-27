#ifdef __CLING__
R__LOAD_LIBRARY(libDelphes)
#include <iostream>
#include <fstream>
#include <cmath>
#include <cstring>
#include <string>
#endif


std::string to_string_without_decimals(double value) {
  std::ostringstream oss;
  oss << std::fixed << std::setprecision(0) << value;
  return oss.str();
}

// Helper function to check if a name starts with a given prefix
auto starts_with = [](const std::string& str, const std::string& prefix) {
  return str.rfind(prefix, 0) == 0;
};


//-------------------------------------------------------------------------------


// save plots as different file types
void PrintCanvas(TCanvas *c=nullptr, string name="default", string outputFolder="default"){
  std::vector <string> types={"png"}; 
  for(std::vector<string>::iterator it=types.begin(); it!=types.end(); it++) {
    c->Print(Form("%s/%s.%s", outputFolder.c_str(), name.c_str(), (*it).c_str()), (*it).c_str());
    //c->SetLogy(); 
    //c->Print(Form("%s/%s_log.%s", outputFolder.c_str(), name.c_str(), (*it).c_str()), (*it).c_str());
  }
}

//-------------------------------------------------------------------------------     


// draw 1d histograms
void draw_stack(TFile *sig_file, TFile *ttbar_file, TFile *ttHbb_file, TFile *diboson_file, TFile *drellyan_file, string name, string title, string axistitle, double sig_scale, string sigName, string outputFolder) {

  // get the histograms from the files
  TH1F *sig_hist = (TH1F*)sig_file->Get(name.c_str());
  TH1F *ttbar_hist = (TH1F*)ttbar_file->Get(name.c_str());
  TH1F *ttHbb_hist = (TH1F*)ttHbb_file->Get(name.c_str());
  TH1F *diboson_hist = (TH1F*)diboson_file->Get(name.c_str());
  TH1F *drellyan_hist = (TH1F*)drellyan_file->Get(name.c_str());

  if ( starts_with(name, "hEff") || starts_with(name, "hClosure") || starts_with(name, "weight") ) {
    std::cout << "Skipping histogram: " << name << std::endl;
    return; // Exit the function early
  }

  TH1F *sigClone=(TH1F*)sig_hist->Clone("sigClone");

  // Add a small constant to each bin to avoid log(0)
  /* double epsilon = 1e-12;
  for (int i = 1; i <= sig_hist->GetNbinsX(); ++i) {
    sig_hist->SetBinContent(i, sig_hist->GetBinContent(i) + epsilon);
    ttbar_hist->SetBinContent(i, ttbar_hist->GetBinContent(i) + epsilon);
    ttHbb_hist->SetBinContent(i, ttHbb_hist->GetBinContent(i) + epsilon);
    diboson_hist->SetBinContent(i, diboson_hist->GetBinContent(i) + epsilon);
    drellyan_hist->SetBinContent(i, drellyan_hist->GetBinContent(i) + epsilon);
    sigClone->SetBinContent(i, sigClone->GetBinContent(i) + epsilon);
  }
  */

  double lumiScaling = 1;
  
  // change lumi if desired
  sig_hist->Scale(lumiScaling);
  ttbar_hist->Scale(lumiScaling);
  diboson_hist->Scale(lumiScaling);
  ttHbb_hist->Scale(lumiScaling);
  drellyan_hist->Scale(lumiScaling);
  sigClone->Scale(lumiScaling);

  sig_hist->Scale(sig_scale);
  ttbar_hist->Scale(1);
  ttHbb_hist->Scale(1);
  drellyan_hist->Scale(1);
  diboson_hist->Scale(1);

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
  ttbar_hist->SetFillColor(kRed);
  ttHbb_hist->SetFillColor(kMagenta);
  diboson_hist->SetFillColor(kOrange);
  drellyan_hist->SetFillColor(kYellow);

  
  // form the histogram stack
  THStack *stack = new THStack(name.c_str(), title.c_str());
  stack->Add(ttbar_hist);
  stack->Add(drellyan_hist);
  stack->Add(diboson_hist);
  stack->Add(ttHbb_hist);
  stack->Add(sig_hist);

  stack->SetMinimum(0);

  // make a legend
  TLegend *legend = new TLegend(0.775, 0.775, 0.875, 0.875);
  legend->AddEntry(sig_hist, (sigName + " x " + to_string_without_decimals(sig_scale)).c_str(), "f"); 
  legend->AddEntry(drellyan_hist, "Drell-Yan", "f");
  legend->AddEntry(ttbar_hist, "ttbar", "f");
  legend->AddEntry(ttHbb_hist, "ttHbb", "f");
  legend->AddEntry(diboson_hist, "Diboson", "f");
  
  // make a canvas and draw on it
  TCanvas *cstack = new TCanvas(name.c_str(), title.c_str(), 1500, 1200);
  cstack->cd();
  stack->Draw("hist");
  total->Draw("hist same pe");
  legend->Draw();

  // sig on top
  
  //int scale = 100;  
  //sigClone->Scale(scale);
  //sigClone->SetLineColor(kBlack);
  //legend->AddEntry(sigClone, ("Signal x " + to_string(scale)).c_str(), "l");
  //sigClone->Draw("hist same");
    
  TString nameOut(name);
  nameOut.ReplaceAll(" ","_");
  nameOut.ReplaceAll("#","");
  nameOut.ReplaceAll("^","");
  nameOut.ReplaceAll("/","");
  nameOut.ReplaceAll("{","");
  nameOut.ReplaceAll("}","");
  cout<<" Name "<<nameOut.Data()<<endl;
  PrintCanvas(cstack,nameOut.Data(), outputFolder);
  cstack->Close();
  delete cstack; 
}

//-------------------------------------------------------------------------------     

// draw 1d histograms
void draw_overlay_normalized(TFile *sig_file, TFile *ttbar_file, TFile *ttHbb_file, TFile *diboson_file, TFile *drellyan_file, string name, string title, string axistitle, string sigName, string outputFolder) {

  // get the histograms from the files
  TH1F *sig_hist = (TH1F*)sig_file->Get(name.c_str());
  TH1F *ttbar_hist = (TH1F*)ttbar_file->Get(name.c_str());
  TH1F *ttHbb_hist = (TH1F*)ttHbb_file->Get(name.c_str());
  TH1F *diboson_hist = (TH1F*)diboson_file->Get(name.c_str());
  TH1F *drellyan_hist = (TH1F*)drellyan_file->Get(name.c_str());

  if (starts_with(name, "hEff") || starts_with(name, "hClosure") || starts_with(name, "weight")) {
    std::cout << "Skipping histogram: " << name << std::endl;
    return; // Exit the function early
  }

  // Normalize the histograms by their integrals
  if (sig_hist->Integral() != 0) sig_hist->Scale(1.0 / sig_hist->Integral());
  if (ttbar_hist->Integral() != 0) ttbar_hist->Scale(1.0 / ttbar_hist->Integral());
  if (ttHbb_hist->Integral() != 0) ttHbb_hist->Scale(1.0 / ttHbb_hist->Integral());
  if (diboson_hist->Integral() != 0) diboson_hist->Scale(1.0 / diboson_hist->Integral());
  if (drellyan_hist->Integral() != 0) drellyan_hist->Scale(1.0 / drellyan_hist->Integral());

  // Set line colors for each histogram
  sig_hist->SetLineColor(kBlack);
  ttbar_hist->SetLineColor(kRed);
  ttHbb_hist->SetLineColor(kGreen);
  diboson_hist->SetLineColor(kOrange);
  drellyan_hist->SetLineColor(kBlue);

  // Set line width for better visibility
  sig_hist->SetLineWidth(1);
  ttbar_hist->SetLineWidth(1);
  ttHbb_hist->SetLineWidth(1);
  diboson_hist->SetLineWidth(1);
  drellyan_hist->SetLineWidth(1);

  sig_hist->SetMinimum(0);
  sig_hist->SetMaximum(1);
  sig_hist->SetStats(false);

  // Make a legend
  TLegend *legend = new TLegend(0.725, 0.725, 0.875, 0.875);
  legend->AddEntry(sig_hist, sigName.c_str(), "l");
  legend->AddEntry(ttbar_hist, "ttbar", "l");
  legend->AddEntry(ttHbb_hist, "ttHbb", "l");
  legend->AddEntry(diboson_hist, "Diboson", "l");
  legend->AddEntry(drellyan_hist, "Drell-Yan", "l");

  // Make a canvas and draw on it
  TCanvas *c = new TCanvas(name.c_str(), title.c_str(), 1500, 1200);
  c->cd();

  // Draw histograms
  sig_hist->Draw("hist");
  ttbar_hist->Draw("hist same");
  ttHbb_hist->Draw("hist same");
  diboson_hist->Draw("hist same");
  drellyan_hist->Draw("hist same");

  // Draw the legend
  legend->Draw();

  // Add axis titles
  sig_hist->GetXaxis()->SetTitle(axistitle.c_str());
  sig_hist->GetYaxis()->SetTitle("Normalized Events");

  // Save the canvas
  TString nameOut(name);
  nameOut.ReplaceAll(" ","_");
  nameOut.ReplaceAll("#","");
  nameOut.ReplaceAll("^","");
  nameOut.ReplaceAll("/","");
  nameOut.ReplaceAll("{","");
  nameOut.ReplaceAll("}","");
  cout<<" Name "<<nameOut.Data()<<endl;
  PrintCanvas(c, nameOut.Data(), outputFolder);
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
  PrintCanvas(c, name, outputFolder);
  c->Close();
}


//-------------------------------------------------------------------------------     

void draw_stacks_(string sig_filename = "../outputs/wpwmhqq.root", string bkg_filename = "../outputs/all_bkg.root",string ttbar_filename = "../outputs/ttbar.root", string ttHbb_filename = "../outputs/ttHbb.root", string diboson_filename = "../outputs/diboson.root", string drellyan_filename = "../outputs/drellyan.root", double sig_scale = 1, string sigName = "zzhjj", string outputFolder = ".") {

  // open files
  TFile *sig_file = TFile::Open(sig_filename.c_str(), "READ");
  TFile *bkg_file = TFile::Open(bkg_filename.c_str(), "READ");
  TFile *ttbar_file = TFile::Open(ttbar_filename.c_str(), "READ");
  TFile *ttHbb_file = TFile::Open(ttHbb_filename.c_str(), "READ");
  TFile *diboson_file = TFile::Open(diboson_filename.c_str(), "READ");
  TFile *drellyan_file = TFile::Open(drellyan_filename.c_str(), "READ");

  // draw histograms
  TIter next(sig_file->GetListOfKeys());
  TKey *key;
  while ((key = (TKey*)next())) {
       
      const char* keyName = key->GetName();
      if (strncmp(keyName, sigName.c_str(), sigName.length()) != 0) continue;

      TObject* obj = key->ReadObj();

      if (obj->IsA()->InheritsFrom("TH2F")) continue;
      if (obj->IsA()->InheritsFrom("TH1F")) {
          TH1* hist = (TH1*)obj;
          const char* name = hist->GetName();
          const char* title = hist->GetTitle();
          const char* xlabel = hist->GetXaxis()->GetTitle(); // Set your x-axis label here
          draw_stack(sig_file, ttbar_file, ttHbb_file, diboson_file, drellyan_file, name, title, xlabel, sig_scale, sigName, outputFolder);
     }
  }
  
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
  
  sig_file->Close();
  bkg_file->Close();
  ttbar_file->Close();
  ttHbb_file->Close();
  diboson_file->Close();
  drellyan_file->Close();

}


//-------------------------------------------------------------------------------  


void draw_stacks(string sig="", int sig_scale = 1){

  system(Form("mkdir -p stack_plots/%s",sig.c_str()));

  draw_stacks_( "../outputs/"+sig+".root","../outputs/all_bkg.root","../outputs/ttbar.root","../outputs/ttHbb.root","../outputs/diboson.root","../outputs/drellyan.root", sig_scale, sig, "stack_plots/"+sig );

}

