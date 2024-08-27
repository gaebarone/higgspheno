#ifdef __CLING__
R__LOAD_LIBRARY(libDelphes)
#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#endif


void plot_eff_table( string sig1_filename, string sig2_filename, string sig3_filename, string sig4_filename, string sig5_filename, string sig6_filename, string ttbar_filename, string ttHbb_filename, string diboson_filename, string drellyan_filename, string all_bkg_filename, string hEff, string outputFolder) {

  TFile *sig1_file = TFile::Open(sig1_filename.c_str(), "READ");
  TFile *sig2_file = TFile::Open(sig2_filename.c_str(), "READ");
  TFile *sig3_file = TFile::Open(sig3_filename.c_str(), "READ");
  TFile *sig4_file = TFile::Open(sig4_filename.c_str(), "READ");
  TFile *sig5_file = TFile::Open(sig5_filename.c_str(), "READ");
  TFile *sig6_file = TFile::Open(sig6_filename.c_str(), "READ");

  TFile *all_bkg_file = TFile::Open(all_bkg_filename.c_str(), "READ");
  TFile *ttbar_file = TFile::Open(ttbar_filename.c_str(), "READ");
  TFile *ttHbb_file = TFile::Open(ttHbb_filename.c_str(), "READ");
  TFile *diboson_file = TFile::Open(diboson_filename.c_str(), "READ");
  TFile *drellyan_file = TFile::Open(drellyan_filename.c_str(), "READ");

  // Create a 2D histogram with bins for the x-axis and the y-axis representing files
  int num_files = 10;
  const int num_bins = ((TH1F*)sig1_file->Get(hEff.c_str()))->GetNbinsX();
  TH2F *h2d = new TH2F( (hEff + "_2D").c_str(), (hEff + "_2D").c_str(), num_bins, 0, num_bins, num_files, 0, num_files);
  
  // Set the y-axis labels for the different files
  h2d->GetYaxis()->SetBinLabel(1, "(H>ZZ)jj");
  h2d->GetYaxis()->SetBinLabel(2, "(H>W+W-)jj");
  h2d->GetYaxis()->SetBinLabel(3, "ZZjj");
  h2d->GetYaxis()->SetBinLabel(4, "W+W-jj");
  h2d->GetYaxis()->SetBinLabel(5, "ZZHjj");
  h2d->GetYaxis()->SetBinLabel(6, "W+W-Hjj");
  //h2d->GetYaxis()->SetBinLabel(3, "All Bkg");
  h2d->GetYaxis()->SetBinLabel(7, "Drell-Yan");
  h2d->GetYaxis()->SetBinLabel(8, "Diboson");
  h2d->GetYaxis()->SetBinLabel(9, "ttbar");
  h2d->GetYaxis()->SetBinLabel(10, "ttHbb");

  TH1F *ref_hist = (TH1F*)sig1_file->Get(hEff.c_str());
  if (!ref_hist) {
    std::cerr << "Histogram " << hEff << " not found in the signal file." << std::endl;
    return;
  }
  for (int bin = 1; bin <= num_bins; ++bin) {
    TString bin_label = ref_hist->GetXaxis()->GetBinLabel(bin);
    h2d->GetXaxis()->SetBinLabel(bin, bin_label);
  }
  h2d->GetXaxis()->SetLabelSize(0.02);

  // File names and their corresponding histograms
  std::vector<std::pair<std::string, TH1F*>> file_histograms = {
    {"(H>ZZ)jj", (TH1F*)sig6_file->Get(hEff.c_str())},
    {"(H>W+W-)jj", (TH1F*)sig5_file->Get(hEff.c_str())},
    {"ZZjj", (TH1F*)sig4_file->Get(hEff.c_str())},
    {"W+W-jj", (TH1F*)sig3_file->Get(hEff.c_str())},
    {"ZZHjj", (TH1F*)sig2_file->Get(hEff.c_str())},
    {"W+W-Hjj", (TH1F*)sig1_file->Get(hEff.c_str())},
    //{"All Bkgs", (TH1F*)all_bkg_file->Get(hEff.c_str())},
    {"Drell-Yan", (TH1F*)drellyan_file->Get(hEff.c_str())},
    {"Diboson", (TH1F*)diboson_file->Get(hEff.c_str())},
    {"ttbar", (TH1F*)ttbar_file->Get(hEff.c_str())},
    {"ttHbb", (TH1F*)ttHbb_file->Get(hEff.c_str())}
  };

  for (int i = 0; i < num_files; ++i) {
    const auto& [label, hist] = file_histograms[i];
    if (hist && hist->GetNbinsX() == num_bins) { // Ensure the histograms have the expected number of bins
      for (int bin = 1; bin <= num_bins; ++bin) {
        double bin_content = hist->GetBinContent(bin);
        // Explicitly set the bin content to zero if it's zero
	if (bin_content==0) {
	  h2d->SetBinContent(bin, i + 1, 0.01);
	} else {
	  h2d->SetBinContent(bin, i + 1, bin_content);
	}
      }
    } else {
      std::cerr << "Histogram for " << label << " is not found or has unexpected number of bins." << std::endl;
    }
  }

  // Draw the 2D histogram
  TCanvas *c = new TCanvas(("c_"+hEff+"_2D").c_str(), (hEff+" 2D Histogram").c_str(), 2000, 1600);
  c->cd();
  h2d->SetStats(false);
  h2d->SetTitle(hEff.c_str());
  h2d->Draw("COLZ");
  gStyle->SetPaintTextFormat("1.0f"); // Set the format of the numbers
  h2d->SetMarkerSize(1);
  h2d->Draw("COLZ TEXT");

  double x_min = h2d->GetXaxis()->GetXmin();
  double x_max = h2d->GetXaxis()->GetXmax();
  TLine *line = new TLine(x_min, num_files-4, x_max, num_files-4);
  TLine *line2 = new TLine(x_max-num_files+4, 0, x_max-num_files+4, num_files);
  line->SetLineColor(kRed);
  line2->SetLineColor(kRed);
  line->SetLineWidth(3);
  line2->SetLineWidth(3);
  line->Draw();
  line2->Draw();

  // Save the canvas
  TString nameOut( hEff + "_2D" );
  nameOut.ReplaceAll(" ", "_");
  nameOut.ReplaceAll("#", "");
  nameOut.ReplaceAll("^", "");
  nameOut.ReplaceAll("/", "");
  nameOut.ReplaceAll("{", "");
  nameOut.ReplaceAll("}", "");

  std::string outputFileName = outputFolder + "/" + nameOut.Data() + ".png";
  std::cout << "Saving histogram to " << outputFileName << std::endl;
  c->SaveAs(outputFileName.c_str());

  // Clean up

  delete line;
  delete line2;

  delete c;
  delete h2d;

  sig1_file->Close();
  sig2_file->Close();
  sig3_file->Close();
  sig4_file->Close();
  sig5_file->Close();
  sig6_file->Close();
  //all_bkg_file->Close();
  ttbar_file->Close();
  ttHbb_file->Close();
  diboson_file->Close();
  drellyan_file->Close();

  delete sig1_file;
  delete sig2_file;
  delete sig3_file;
  delete sig4_file;
  delete sig5_file;
  delete sig6_file;
  //delete all_bkg_file;
  delete ttbar_file;
  delete ttHbb_file;
  delete diboson_file;
  delete drellyan_file;

}

//-------------------------------------------------------------------------------  


void draw_eff_table(string sig1="", string sig2="", string sig3="", string sig4="", string sig5="", string sig6=""){

  system("mkdir -p eff_tables");

  plot_eff_table( "../outputs/"+sig1+".root","../outputs/"+sig2+".root","../outputs/"+sig3+".root","../outputs/"+sig4+".root", "../outputs/"+sig5+".root","../outputs/"+sig6+".root", "../outputs/ttbar.root","../outputs/ttHbb.root","../outputs/diboson.root","../outputs/drellyan.root","../outputs/all_bkg.root", "hEff_reco",  "eff_tables" );
  plot_eff_table( "../outputs/"+sig1+".root","../outputs/"+sig2+".root","../outputs/"+sig3+".root","../outputs/"+sig4+".root", "../outputs/"+sig5+".root","../outputs/"+sig6+".root", "../outputs/ttbar.root","../outputs/ttHbb.root","../outputs/diboson.root","../outputs/drellyan.root","../outputs/all_bkg.root", "hEff_particle",  "eff_tables" );
  plot_eff_table( "../outputs/"+sig1+".root","../outputs/"+sig2+".root","../outputs/"+sig3+".root","../outputs/"+sig4+".root", "../outputs/"+sig5+".root","../outputs/"+sig6+".root", "../outputs/ttbar.root","../outputs/ttHbb.root","../outputs/diboson.root","../outputs/drellyan.root","../outputs/all_bkg.root", "hEff_parton", "eff_tables" );

}

