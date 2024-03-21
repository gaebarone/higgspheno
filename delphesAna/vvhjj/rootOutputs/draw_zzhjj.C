#ifdef __CLING__
R__LOAD_LIBRARY(libDelphes)
#include <iostream>
#include <fstream>
#include <cmath>
#endif

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
  sig_hist->Scale(lumiScaling); //*0.000013 for W
  ttbar_hist->Scale(lumiScaling);
  diboson_hist->Scale(lumiScaling);
  ttHbb_hist->Scale(lumiScaling);
  drellyan_hist->Scale(lumiScaling);
  sigClone->Scale(lumiScaling); //*0.000013 for W    

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
  TIter next(sig_file->GetListOfKeys());
  TKey *key;
  while ((key = (TKey*)next())) {
      TObject* obj = key->ReadObj();
      if (obj->IsA()->InheritsFrom("TH1")) {
          TH1* hist = (TH1*)obj;
          const char* name = hist->GetName();
          const char* title = hist->GetTitle();
          const char* xlabel = ""; // Set your x-axis label here
          draw_stack(sig_file, ttbar_file, ttHbb_file, diboson_file, drellyan_file, name, title, xlabel, outputFolder);
      }
  }
  
/* 
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
*/

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

