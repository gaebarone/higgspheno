#ifdef __CLING__
R__LOAD_LIBRARY(libDelphes)
#include <iostream>
#include <fstream>
#include <cmath>
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
void draw_overlay_normalized(string name, string title, string axistitle, string outputFolder, string wwhjj, string zzhjj, string wwjj, string zzjj, string hwwjj, string hzzjj) {

  // open backgrounds
  TFile *ttbar_file = TFile::Open("../outputs/ttbar.root", "READ");
  TFile *ttHbb_file = TFile::Open("../outputs/ttHbb.root", "READ");
  TFile *diboson_file = TFile::Open("../outputs/diboson.root", "READ");
  TFile *drellyan_file = TFile::Open("../outputs/drellyan.root", "READ");
  
  // open signals
  TFile *wwhjj_file = TFile::Open(("../outputs/"+wwhjj+".root").c_str(), "READ");
  TFile *zzhjj_file = TFile::Open(("../outputs/"+zzhjj+".root").c_str(), "READ");
  TFile *wwjj_file = TFile::Open(("../outputs/"+wwjj+".root").c_str(), "READ");
  TFile *zzjj_file = TFile::Open(("../outputs/"+zzjj+".root").c_str(), "READ");
  TFile *hwwjj_file = TFile::Open(("../outputs/"+hwwjj+".root").c_str(), "READ");
  TFile *hzzjj_file = TFile::Open(("../outputs/"+hzzjj+".root").c_str(), "READ");

  // get the histograms from the files
  TH1F *ttbar_hist = (TH1F*)ttbar_file->Get(name.c_str());
  TH1F *ttHbb_hist = (TH1F*)ttHbb_file->Get(name.c_str());
  TH1F *diboson_hist = (TH1F*)diboson_file->Get(name.c_str());
  TH1F *drellyan_hist = (TH1F*)drellyan_file->Get(name.c_str());

  TH1F *wwhjj_hist = (TH1F*)wwhjj_file->Get(name.c_str());
  TH1F *zzhjj_hist = (TH1F*)zzhjj_file->Get(name.c_str());
  TH1F *wwjj_hist = (TH1F*)wwjj_file->Get(name.c_str());
  TH1F *zzjj_hist = (TH1F*)zzjj_file->Get(name.c_str());
  TH1F *hwwjj_hist = (TH1F*)hwwjj_file->Get(name.c_str());
  TH1F *hzzjj_hist = (TH1F*)hzzjj_file->Get(name.c_str());

  if (starts_with(name, "hEff") || starts_with(name, "hClosure") || starts_with(name, "weight")) {
    std::cout << "Skipping histogram: " << name << std::endl;
    return; // Exit the function early
  }

  // Add a small constant to each bin to avoid log(0)                                                                                                                                                                                                                          
  /*
  double epsilon = 1e-12;
  for (int i = 1; i <= ttbar_hist->GetNbinsX(); ++i) {
    ttbar_hist->SetBinContent(i, ttbar_hist->GetBinContent(i) + epsilon);
    ttHbb_hist->SetBinContent(i, ttHbb_hist->GetBinContent(i) + epsilon);
    diboson_hist->SetBinContent(i, diboson_hist->GetBinContent(i) + epsilon);
    drellyan_hist->SetBinContent(i, drellyan_hist->GetBinContent(i) + epsilon);
    wwhjj_hist->SetBinContent(i, wwhjj_hist->GetBinContent(i) + epsilon);
    zzhjj_hist->SetBinContent(i, zzhjj_hist->GetBinContent(i) + epsilon);
    //wwjj_hist->SetBinContent(i, wwjj_hist->GetBinContent(i) + epsilon);
    //zzjj_hist->SetBinContent(i, zzjj_hist->GetBinContent(i) + epsilon);
    //hwwjj_hist->SetBinContent(i, hwwjj_hist->GetBinContent(i) + epsilon);
    //hzzjj_hist->SetBinContent(i, hzzjj_hist->GetBinContent(i) + epsilon);
  }
  */

  // Normalize the histograms by their integrals
  if (ttbar_hist->Integral() != 0) ttbar_hist->Scale(1.0 / ttbar_hist->Integral());
  if (ttHbb_hist->Integral() != 0) ttHbb_hist->Scale(1.0 / ttHbb_hist->Integral());
  if (diboson_hist->Integral() != 0) diboson_hist->Scale(1.0 / diboson_hist->Integral());
  if (drellyan_hist->Integral() != 0) drellyan_hist->Scale(1.0 / drellyan_hist->Integral());

  if (wwhjj_hist->Integral() != 0) wwhjj_hist->Scale(1.0 / wwhjj_hist->Integral());
  if (zzhjj_hist->Integral() != 0) zzhjj_hist->Scale(1.0 / zzhjj_hist->Integral());
  if (wwjj_hist->Integral() != 0) wwjj_hist->Scale(1.0 / wwjj_hist->Integral());
  if (zzjj_hist->Integral() != 0) zzjj_hist->Scale(1.0 / zzjj_hist->Integral());
  if (hwwjj_hist->Integral() != 0) hwwjj_hist->Scale(1.0 / hwwjj_hist->Integral());
  if (hzzjj_hist->Integral() != 0) hzzjj_hist->Scale(1.0 / hzzjj_hist->Integral());

  // Set line colors for each histogram
  ttbar_hist->SetLineColor(kRed+1);
  ttHbb_hist->SetLineColor(kRed-4);
  diboson_hist->SetLineColor(kOrange+8);
  drellyan_hist->SetLineColor(kOrange-2);

  wwhjj_hist->SetLineColor(kBlue);
  zzhjj_hist->SetLineColor(kGreen);
  wwjj_hist->SetLineColor(kBlue+2);
  zzjj_hist->SetLineColor(kGreen+2);
  hwwjj_hist->SetLineColor(kMagenta);
  hzzjj_hist->SetLineColor(kViolet);


  ttbar_hist->SetMinimum(0);
  ttbar_hist->SetMaximum(1.1);
  ttbar_hist->SetStats(false);

  wwhjj_hist->SetMinimum(0);
  wwhjj_hist->SetMaximum(1.1);
  wwhjj_hist ->SetStats(false);

  // Make a legend
  TLegend *legend = new TLegend(0.725, 0.725, 0.875, 0.875);
  
  legend->AddEntry(ttbar_hist, "ttbar", "l");
  legend->AddEntry(ttHbb_hist, "ttHbb", "l");
  legend->AddEntry(diboson_hist, "Diboson", "l");
  legend->AddEntry(drellyan_hist, "Drell-Yan", "l");
  
  legend->AddEntry(wwhjj_hist, wwhjj.c_str(), "l");
  legend->AddEntry(zzhjj_hist, zzhjj.c_str(), "l");
  legend->AddEntry(wwjj_hist, wwjj.c_str(), "l");
  legend->AddEntry(zzjj_hist, zzjj.c_str(), "l");
  legend->AddEntry(hwwjj_hist, hwwjj.c_str(), "l");
  legend->AddEntry(hzzjj_hist, hzzjj.c_str(), "l");

  // Make a canvas and draw on it
  TCanvas *c = new TCanvas(name.c_str(), title.c_str(), 1500, 1200);
  c->cd();

  // Draw histograms
  ttbar_hist->Draw("hist");
  ttHbb_hist->Draw("hist same");
  diboson_hist->Draw("hist same");
  drellyan_hist->Draw("hist same");

  wwhjj_hist->Draw("hist same");
  zzhjj_hist->Draw("hist same");
  wwjj_hist->Draw("hist same");
  zzjj_hist->Draw("hist same");
  hwwjj_hist->Draw("hist same");
  hzzjj_hist->Draw("hist same");

  // Draw the legend
  legend->Draw();

  // Add axis titles
  ttbar_hist->GetXaxis()->SetTitle(axistitle.c_str());
  ttbar_hist->GetYaxis()->SetTitle("Normalized Events");

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

  wwhjj_file->Close();
  zzhjj_file->Close();
  wwjj_file->Close();
  zzjj_file->Close();
  hwwjj_file->Close();
  hzzjj_file->Close();

  ttbar_file->Close();
  ttHbb_file->Close();
  diboson_file->Close();
  drellyan_file->Close();

  delete wwhjj_file;
  delete zzhjj_file;
  delete wwjj_file;
  delete zzjj_file;
  delete hwwjj_file;
  delete hzzjj_file;

  delete ttbar_file;
  delete ttHbb_file;
  delete diboson_file;
  delete drellyan_file;

}


//-------------------------------------------------------------------------------     


void draw_overlay_normalized_(string outputFolder, string wwhjj, string zzhjj, string wwjj, string zzjj, string hwwjj, string hzzjj) {

  TFile *ttbar_file = TFile::Open("../outputs/ttbar.root", "READ");

  // draw hists
  TIter next(ttbar_file->GetListOfKeys());
  TKey *key;
  while ((key = (TKey*)next())) {

      TObject* obj = key->ReadObj();

      const char* keyName = key->GetName();
      if (strncmp(keyName, "nocuts_", 7) != 0) continue;

      if (obj->IsA()->InheritsFrom("TH2F")) continue;

      if (obj->IsA()->InheritsFrom("TH1F")) {
          TH1* hist = (TH1*)obj;
          const char* name = hist->GetName();
          const char* title = hist->GetTitle();
          const char* xlabel = hist->GetXaxis()->GetTitle(); // Set your x-axis label here
	  draw_overlay_normalized(name, title, xlabel, outputFolder, wwhjj, zzhjj, wwjj, zzjj, hwwjj, hzzjj);
	    }
  }
  
}


//-------------------------------------------------------------------------------  


void draw_normalized_dists( string wwhjj, string zzhjj, string wwjj, string zzjj, string hwwjj, string hzzjj) {

  system("mkdir -p normalized_dists");

  draw_overlay_normalized_( "normalized_dists", wwhjj, zzhjj, wwjj, zzjj, hwwjj, hzzjj );

}

