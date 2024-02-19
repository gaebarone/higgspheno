#ifndef HIST_INCLUDE_H
#define HIST_INCLUDE_H

//------------------------------------------------------------------------------------------------------------------------------------------------------------
// HISTOGRAM DEFINITIONS
//------------------------------------------------------------------------------------------------------------------------------------------------------------

void PrintCanvas(TCanvas *c=nullptr,string name="default"){
  std::vector <string> types={"jpg"}; 
  for(std::vector<string>::iterator it=types.begin(); it!=types.end(); it++) {
    c->Print(Form("jpg/%s.%s",name.c_str(),(*it).c_str()),(*it).c_str());
  }
}

void draw_hist(TH1 *histo, const char *name, const char *title, const char *axistitle) {
  TCanvas *c = new TCanvas(name, title, 1500,1200);
  c->cd();
  histo->GetXaxis()->SetTitle(axistitle);
  histo->SetMinimum(0.0);
  histo->Draw("hist e");
  PrintCanvas(c, name);
}

void draw_hist2(TH2F*histo, const char *name, const char *title, const char *xaxistitle, const char *yaxistitle) {
  TCanvas *c = new TCanvas(name, title, 1500, 1200);
  histo->GetXaxis()->SetTitle(xaxistitle);
  histo->GetYaxis()->SetTitle(yaxistitle);
  //histo->SetStatX(0.875);
  //histo->SetStatY(0.875);
  histo->Draw("COLZ");
  PrintCanvas(c, name);
}


//------------------------------------------------------------------------------------------------------------------------------------------------------------
// BOOK HISTOGRAMS
//------------------------------------------------------------------------------------------------------------------------------------------------------------

/*
const double mBins = 20;
const double pTBins = 20;
const double phiBins = 10;
const double etaBins = 10;
const double RBins = 10;
const double cosBins = 10;

const double hpTmin = 0;
const double jpTmin = 0;
const double zpTmin = 0;
const double lpTmin = 0;
const double hpTmax = 500;
const double jpTmax = 750;
const double zpTmax = 500;
const double lpTmax = 750;

const double hmmin = 0;
const double jmmin = 0;
const double zmmin = 0;
const double hmmax = 200;
const double jmmax = 200;
const double zmmax = 200;

const double hetamin = -2.5;
const double jetamin = 0;
const double zetamin = 0;
const double letamin = 0;
const double hetamax = 2.5;
const double jetamax = 5;
const double zetamax = 5;
const double letamax = 5;

const double hRmin = 0;
const double jRmin = 0;
const double zRmin = 0;
const double lRmin = 0;
const double hRmax = 5;
const double jRmax = 5;
const double zRmax = 5;
const double lRmax = 5;

  vector <TH1F*> listOfTH1;
  vector <TH2F*> listOfTH2;

void book_histograms() {
// 1D

  // higgs - reco
    // pT + m
      TH1F *hHpTreco = new TH1F("hbb_pT_reco", "p^{T}_{hbb}_reco", pTBins, hpTmin, hpTmax); 	listOfTH1.push_back(hHpTreco);
      TH1F *hHmreco = new TH1F("hbb_m_reco", "m_{hbb}_reco", mBins, hmmin, hmmax); listOfTH1.push_back(hHmreco);
      //TH1F *hb1pTreco = new TH1F("b1_pT_reco", "p^{T}_{b1}_reco", pTBins, hpTmin, hpTmax); listOfTH1.push_back(hb1pTreco);
      //TH1F *hb1mreco = new TH1F("b1_m_reco", "m_{b1}_reco", mBins, hmmin, hmmax);listOfTH1.push_back(hb1mreco);
      //TH1F *hb2pTreco = new TH1F("b2_pT_reco", "p^{T}_{b2}_reco", pTBins, hpTmin, hpTmax); listOfTH1.push_back(hb2pTreco);
      //TH1F *hb2mreco = new TH1F("b2_m_reco", "m_{b2}_reco", mBins, hmmin, hmmax); listOfTH1.push_back(hb2mreco);
    // phi
      //TH1F *hHphireco = new TH1F("hbb_#phi_reco", "#phi_{hbb}_reco", phiBins, -TMath::Pi(), +TMath::Pi()); listOfTH1.push_back(hHphireco);
      //TH1F *hb1phireco = new TH1F("b1_#phi_reco", "#phi_{b1}_reco", phiBins, -TMath::Pi(), +TMath::Pi()); listOfTH1.push_back(hb1phireco);
      //TH1F *hb2phireco = new TH1F("b2_#phi_reco", "#phi_{b2}_reco", phiBins, -TMath::Pi(), +TMath::Pi()); listOfTH1.push_back(hb2phireco);
      TH1F *hbbdeltaPhireco = new TH1F("bb_#Delta#phi_reco", "#Delta#phi_{bb}_reco", phiBins, -TMath::Pi(), +TMath::Pi()); listOfTH1.push_back(hbbdeltaPhireco);
    // eta
      //TH1F *hHetareco = new TH1F("hbb_#eta_reco", "#eta_{hbb}_reco", etaBins, hetamin, hetamax); listOfTH1.push_back(hHetareco);
      //TH1F *hb1etareco = new TH1F("b1_#eta_reco", "#eta_{b1}_reco", etaBins, hetamin, hetamax); listOfTH1.push_back(hb1etareco);
      //TH1F *hb2etareco = new TH1F("b2_#eta_reco", "#eta_{b2}_reco", etaBins, hetamin, hetamax); listOfTH1.push_back(hb2etareco);
      TH1F *hbbdeltaEtareco = new TH1F("bb_#Delta#eta_reco", "#Delta#eta_{bb}_reco", etaBins, hetamin, hetamax); listOfTH1.push_back(hbbdeltaEtareco);
     // R
      //TH1F *hHRreco = new TH1F("hbb_R_reco", "R_{hbb}_reco", RBins, hRmin, hRmax); listOfTH1.push_back(hHRreco);
      //TH1F *hb1Rreco = new TH1F("b1_R_reco", "R_{b1}_reco", RBins, hRmin, hRmax); listOfTH1.push_back(hb1Rreco);
      //TH1F *hb2Rreco = new TH1F("b2_R_reco", "R_{b2}_reco", RBins, hRmin, hRmax); listOfTH1.push_back(hb2Rreco);
      TH1F *hbbdeltaRreco = new TH1F("bb_#DeltaR_reco", "#DeltaR_{bb}_reco", RBins, hRmin, hRmax); listOfTH1.push_back(hbbdeltaRreco);
  
  // higgs - particle
    // pT + m
      TH1F *hHpTparticle = new TH1F("hbb_pT_particle", "p^{T}_{hbb}_particle", pTBins, hpTmin, hpTmax); listOfTH1.push_back(hHpTparticle);
      TH1F *hHmparticle = new TH1F("hbb_m_particle", "m_{hbb}_particle", mBins, hmmin, hmmax); listOfTH1.push_back(hHmparticle);
      //TH1F *hb1pTparticle = new TH1F("b1_pT_particle", "p^{T}_{b1}_particle", pTBins, hpTmin, hpTmax); listOfTH1.push_back(hb1pTparticle);
      //TH1F *hb1mparticle = new TH1F("b1_m_particle", "m_{b1}_particle", mBins, hmmin, hmmax); listOfTH1.push_back(hb1mparticle);
      //TH1F *hb2pTparticle = new TH1F("b2_pT_particle", "p^{T}_{b2}_particle", pTBins, hpTmin, hpTmax); listOfTH1.push_back(hb2pTparticle);
      //TH1F *hb2mparticle = new TH1F("b2_m_particle", "m_{b2}_particle", mBins, hmmin, hmmax); listOfTH1.push_back(hb2mparticle);
    // phi
      //TH1F *hHphiparticle = new TH1F("hbb_#phi_particle", "#phi_{hbb}_particle", phiBins, -TMath::Pi(), +TMath::Pi()); listOfTH1.push_back(hHphiparticle);
      //TH1F *hb1phiparticle = new TH1F("b1_#phi_particle", "#phi_{b1}_particle", phiBins, -TMath::Pi(), +TMath::Pi()); listOfTH1.push_back(hb1phiparticle);
      //TH1F *hb2phiparticle = new TH1F("b2_#phi_particle", "#phi_{b2}_particle", phiBins, -TMath::Pi(), +TMath::Pi()); listOfTH1.push_back(hb2phiparticle);
      TH1F *hbbdeltaPhiparticle = new TH1F("bb_#Delta#phi_particle", "#Delta#phi_{bb}_particle", phiBins, -TMath::Pi(), +TMath::Pi()); listOfTH1.push_back(hbbdeltaPhiparticle);
    // eta
      //TH1F *hHetaparticle = new TH1F("hbb_#eta_particle", "#eta_{hbb}_particle", etaBins, hetamin, hetamax); listOfTH1.push_back(hHetaparticle);
      //TH1F *hb1etaparticle = new TH1F("b1_#eta_particle", "#eta_{b1}_particle", etaBins, hetamin, hetamax); listOfTH1.push_back(hb1etaparticle);
      //TH1F *hb2etaparticle = new TH1F("b2_#eta_particle", "#eta_{b2}_particle", etaBins, hetamin, hetamax); listOfTH1.push_back(hb2etaparticle);
      TH1F *hbbdeltaEtaparticle = new TH1F("bb_#Delta#eta_particle", "#Delta#eta_{bb}_particle", etaBins, hetamin, hetamax); listOfTH1.push_back(hbbdeltaEtaparticle);
   // R
      //TH1F *hHRparticle = new TH1F("hbb_R_particle", "R_{hbb}_particle", RBins, hRmin, hRmax); listOfTH1.push_back(hHRparticle);
      //TH1F *hb1Rparticle = new TH1F("b1_R_particle", "R_{b1}_particle", RBins, hRmin, hRmax); listOfTH1.push_back(hb1Rparticle);
      //TH1F *hb2Rparticle = new TH1F("b2_R_particle", "R_{b2}_particle", RBins, hRmin, hRmax); listOfTH1.push_back(hb2Rparticle);
      TH1F *hbbdeltaRparticle = new TH1F("bb_#DeltaR_particle", "#DeltaR_{bb}_particle", RBins, hRmin, hRmax); listOfTH1.push_back(hbbdeltaRparticle);
  
  // higgs - parton
    // pT + m
      TH1F *hHpTparton = new TH1F("hbb_pT_parton", "p^{T}_{hbb}_parton", pTBins, hpTmin, hpTmax); listOfTH1.push_back(hHpTparton);
      TH1F *hHmparton = new TH1F("hbb_m_parton", "m_{hbb}_parton", mBins, hmmin,  hmmax); listOfTH1.push_back(hHmparton);
      //TH1F *hb1pTparton = new TH1F("b1_pT_parton", "p^{T}_{b1}_parton", pTBins, hpTmin, hpTmax); listOfTH1.push_back(hb1pTparton);
      //TH1F *hb1mparton = new TH1F("b1_m_parton", "m_{b1}_parton", mBins, hmmin, hmmax); listOfTH1.push_back(hb1mparton);
      //TH1F *hb2pTparton = new TH1F("b2_pT_parton", "p^{T}_{b2}_parton", pTBins, hpTmin, hpTmax); listOfTH1.push_back(hb2pTparton);
      //TH1F *hb2mparton = new TH1F("b2_m_parton", "m_{b2}_parton", mBins, hmmin, hmmax); listOfTH1.push_back(hb2mparton);
    // phi
      //TH1F *hHphiparton = new TH1F("hbb_#phi_parton", "#phi_{hbb}_parton", phiBins, -TMath::Pi(), +TMath::Pi()); listOfTH1.push_back(hHphiparton);
      //TH1F *hb1phiparton = new TH1F("b1_#phi_parton", "#phi_{b1}_parton", phiBins, -TMath::Pi(), +TMath::Pi()); listOfTH1.push_back(hb1phiparton);
      //TH1F *hb2phiparton = new TH1F("b2_#phi_parton", "#phi_{b2}_parton", phiBins, -TMath::Pi(), +TMath::Pi()); listOfTH1.push_back(hb2phiparton);
      TH1F *hbbdeltaPhiparton = new TH1F("bb_#Delta#phi_parton", "#Delta#phi_{bb}_parton", phiBins, -TMath::Pi(), +TMath::Pi()); listOfTH1.push_back(hbbdeltaPhiparton);
    // eta
      //TH1F *hHetaparton = new TH1F("hbb_#eta_parton", "#eta_{hbb}_parton", etaBins, hetamin, hetamax); listOfTH1.push_back(hHetaparton);
      //TH1F *hb1etaparton = new TH1F("b1_#eta_parton", "#eta_{b1}_parton", etaBins, hetamin, hetamax); listOfTH1.push_back(hb1etaparton);
      //TH1F *hb2etaparton = new TH1F("b2_#eta_parton", "#eta_{b2}_parton", etaBins, hetamin, hetamax); listOfTH1.push_back(hb2etaparton);
      TH1F *hbbdeltaEtaparton = new TH1F("bb_#Delta#eta_parton", "#Delta#eta_{bb}_parton", etaBins, hetamin, hetamax); listOfTH1.push_back(hbbdeltaEtaparton);
    // R
      //TH1F *hHRparton = new TH1F("hbb_R_parton", "R_{hbb}_parton", RBins, hRmin, hRmax); listOfTH1.push_back(hHRparton);
      //TH1F *hb1Rparton = new TH1F("b1_R_parton", "R_{b1}_parton", RBins, hRmin, hRmax); listOfTH1.push_back(hb1Rparton);
      //TH1F *hb2Rparton = new TH1F("b2_R_parton", "R_{b2}_parton", RBins, hRmin, hRmax); listOfTH1.push_back(hb2Rparton);
      TH1F *hbbdeltaRparton = new TH1F("bb_#DeltaR_parton", "#DeltaR_{bb}_parton", RBins, hRmin, hRmax); listOfTH1.push_back(hbbdeltaRparton);
  
  // vbfj - reco
    // pT + m
      TH1F *hjjpTreco = new TH1F("jj_pT_reco", "p^{T}_{jj}_reco", pTBins, jpTmin, jpTmax); listOfTH1.push_back(hjjpTreco);
      //TH1F *hj1pTreco = new TH1F("j1_pT_reco", "p^{T}_{j1}_reco", pTBins, jpTmin, jpTmax); listOfTH1.push_back(hj1pTreco);
      //TH1F *hj2pTreco = new TH1F("j2_pT_reco", "p^{T}_{j2}_reco", pTBins, jpTmin, jpTmax); listOfTH1.push_back(hj2pTreco);
    // phi
      //TH1F *hj1phireco = new TH1F("j1_#phi_reco", "#phi_{j1}_reco", phiBins, -TMath::Pi(), +TMath::Pi()); listOfTH1.push_back(hj1phireco);
      //TH1F *hj2phireco = new TH1F("j2_#phi_reco", "#phi_{j2}_reco", phiBins, -TMath::Pi(), +TMath::Pi()); listOfTH1.push_back(hj2phireco);
      TH1F *hjjdeltaPhireco = new TH1F("jj_#Delta#phi_reco", "#Delta#phi_{jj}_reco", phiBins, -TMath::Pi(), +TMath::Pi()); listOfTH1.push_back(hjjdeltaPhireco);
    // eta
      //TH1F *hj1etareco = new TH1F("j1_#eta_reco", "#eta_{j1}_reco", etaBins, jetamin, jetamax); listOfTH1.push_back(hj1etareco);
      //TH1F *hj2etareco = new TH1F("j2_#eta_reco", "#eta_{j2}_reco", etaBins, jetamin, jetamax); listOfTH1.push_back(hj2etareco);
      TH1F *hjjdeltaEtareco = new TH1F("jj_#Delta#eta_reco", "#Delta#eta_{jj}_reco", etaBins, jetamin, jetamax); listOfTH1.push_back(hjjdeltaEtareco);
    // R
      //TH1F *hj1Rreco = new TH1F("j1_R_reco", "R_{j1}_reco", RBins, jRmin, jRmax); listOfTH1.push_back(hj1Rreco);
      //TH1F *hj2Rreco = new TH1F("j2_R_reco", "R_{j2}_reco", RBins, jRmin, jRmax); listOfTH1.push_back(hj2Rreco);
      TH1F *hjjdeltaRreco = new TH1F("jj_#DeltaR_reco", "#DeltaR_{jj}_reco", RBins, jRmin, jRmax); listOfTH1.push_back(hjjdeltaRreco);
  
  // vbfj - particle
    // pT + m
      TH1F *hjjpTparticle = new TH1F("jj_pT_particle", "p^{T}_{jj}_particle", pTBins, jpTmin, jpTmax); listOfTH1.push_back(hjjpTparticle);
      //TH1F *hj1pTparticle = new TH1F("j1_pT_particle", "p^{T}_{j1}_particle", pTBins, jpTmin, jpTmax); listOfTH1.push_back(hj1pTparticle);
      //TH1F *hj2pTparticle = new TH1F("j2_pT_particle", "p^{T}_{j2}_particle", pTBins, jpTmin, jpTmax); listOfTH1.push_back(hj2pTparticle);
    // phi
      //TH1F *hj1phiparticle = new TH1F("j1_#phi_particle", "#phi_{j1}_particle", phiBins, -TMath::Pi(), +TMath::Pi()); listOfTH1.push_back(hj1phiparticle);
      //TH1F *hj2phiparticle = new TH1F("j2_#phi_particle", "#phi_{j2}_particle", phiBins, -TMath::Pi(), +TMath::Pi()); listOfTH1.push_back(hj2phiparticle);
      TH1F *hjjdeltaPhiparticle = new TH1F("jj_#Delta#phi_particle", "#Delta#phi_{jj}_particle", phiBins, -TMath::Pi(), +TMath::Pi()); listOfTH1.push_back(hjjdeltaPhiparticle);
    // eta
      //TH1F *hj1etaparticle = new TH1F("j1_#eta_particle", "#eta_{j1}_particle", etaBins, jetamin, jetamax); listOfTH1.push_back(hj1etaparticle);
      //TH1F *hj2etaparticle = new TH1F("j2_#eta_particle", "#eta_{j2}_particle", etaBins, jetamin, jetamax); listOfTH1.push_back(hj2etaparticle);
      TH1F *hjjdeltaEtaparticle = new TH1F("jj_#Delta#eta_particle", "#Delta#eta_{jj}_particle", etaBins, jetamin, jetamax); listOfTH1.push_back(hjjdeltaEtaparticle);
    // R
      //TH1F *hj1Rparticle = new TH1F("j1_R_particle", "R_{j1}_particle", RBins, jRmin, jRmax); listOfTH1.push_back(hj1Rparticle);
      //TH1F *hj2Rparticle = new TH1F("j2_R_particle", "R_{j2}_particle", RBins, jRmin, jRmax); listOfTH1.push_back(hj2Rparticle);
      TH1F *hjjdeltaRparticle = new TH1F("jj_#DeltaR_particle", "#DeltaR_{jj}_particle", RBins, jRmin, jRmax); listOfTH1.push_back(hjjdeltaRparticle);
      
  // z - reco
    // pT + m
      TH1F *hz1pTreco = new TH1F("z1_pT_reco", "p^{T}_{z1}_reco", pTBins, zpTmin, zpTmax); listOfTH1.push_back(hz1pTreco);
      TH1F *hz2pTreco = new TH1F("z2_pT_reco", "p^{T}_{z2}_reco", pTBins, zpTmin, zpTmax); listOfTH1.push_back(hz2pTreco);
      TH1F *hz1mreco = new TH1F("z1_m_reco", "m_{z1}_reco", mBins, zmmin, zmmax); listOfTH1.push_back(hz1mreco);
      TH1F *hz2mreco = new TH1F("z2_m_reco", "m_{z2}_reco", mBins, zmmin, zmmax); listOfTH1.push_back(hz2mreco);
    // phi
      //TH1F *hz1phireco = new TH1F("z1_#phi_reco", "#phi_{z1}_reco", phiBins, -TMath::Pi(), +TMath::Pi()); listOfTH1.push_back(hz1phireco);
      //TH1F *hz2phireco = new TH1F("z2_#phi_reco", "#phi_{z2}_reco", phiBins, -TMath::Pi(), +TMath::Pi()); listOfTH1.push_back(hz2phireco);
      TH1F *hzzdeltaPhireco = new TH1F("zz_#Delta#phi_reco", "#Delta#phi_{zz}_reco", phiBins, -TMath::Pi(), +TMath::Pi()); listOfTH1.push_back(hzzdeltaPhireco);
    // eta
      //TH1F *hz1etareco = new TH1F("z1_#eta_reco", "#eta_{z1}_reco", etaBins, zetamin, zetamax);listOfTH1.push_back(hz1etareco);
      //TH1F *hz2etareco = new TH1F("z2_#eta_reco", "#eta_{z2}_reco", etaBins, zetamin, zetamax);listOfTH1.push_back(hz2etareco);
      TH1F *hzzdeltaEtareco = new TH1F("zz_#Delta#eta_reco", "#Delta#eta_{zz}_reco", etaBins, zetamin, zetamax);listOfTH1.push_back(hzzdeltaEtareco);
    // R
      //TH1F *hz1Rreco = new TH1F("z1_R_reco", "R_{z1}_reco", RBins, zRmin, zRmax); listOfTH1.push_back(hz1Rreco);
      //TH1F *hz2Rreco = new TH1F("z2_R_reco", "R_{z2}_reco", RBins, zRmin, zRmax); listOfTH1.push_back(hz2Rreco);
      TH1F *hzzdeltaRreco = new TH1F("zz_#DeltaR_reco", "#DeltaR_{zz}_reco", RBins, zRmin, zRmax); listOfTH1.push_back(hzzdeltaRreco);
    // cos 
      TH1F *hz1cosThetareco = new TH1F("z1_cos#theta_reco", "cos#theta_{z1}_reco", cosBins, -1, 1); listOfTH1.push_back(hz1cosThetareco);
      TH1F *hz2cosThetareco = new TH1F("z2_cos#theta_reco", "cos#theta_{z2}_reco", cosBins, -1, 1); listOfTH1.push_back(hz2cosThetareco);
  
  // z - particle
    // pT + m
      TH1F *hz1pTparticle = new TH1F("z1_pT_particle", "p^{T}_{z1}_particle", pTBins, zpTmin, zpTmax);listOfTH1.push_back(hz1pTparticle);
      TH1F *hz2pTparticle = new TH1F("z2_pT_particle", "p^{T}_{z2}_particle", pTBins, zpTmin, zpTmax);listOfTH1.push_back(hz2pTparticle);
      TH1F *hz1mparticle = new TH1F("z1_m_particle", "m_{z1}_particle", mBins, zmmin, zmmax);listOfTH1.push_back(hz1mparticle);
      TH1F *hz2mparticle = new TH1F("z2_m_particle", "m_{z2}_particle", mBins, zmmin, zmmax);listOfTH1.push_back(hz2mparticle);
    // phi
      //TH1F *hz1phiparticle = new TH1F("z1_#phi_particle", "#phi_{z1}_particle", phiBins, -TMath::Pi(), +TMath::Pi());listOfTH1.push_back(hz1phiparticle);
      //TH1F *hz2phiparticle = new TH1F("z2_#phi_particle", "#phi_{z2}_particle", phiBins, -TMath::Pi(), +TMath::Pi());listOfTH1.push_back(hz2phiparticle);
      TH1F *hzzdeltaPhiparticle = new TH1F("zz_#Delta#phi_particle", "#Delta#phi_{zz}_particle", phiBins, -TMath::Pi(), +TMath::Pi()); listOfTH1.push_back(hzzdeltaPhiparticle);
    // eta
      //TH1F *hz1etaparticle = new TH1F("z1_#eta_particle", "#eta_{z1}_particle", etaBins, zetamin, zetamax); listOfTH1.push_back(hz1etaparticle);
      //TH1F *hz2etaparticle = new TH1F("z2_#eta_particle", "#eta_{z2}_particle", etaBins, zetamin, zetamax); listOfTH1.push_back(hz2etaparticle);
      TH1F *hzzdeltaEtaparticle = new TH1F("zz_#Delta#eta_particle", "#Delta#eta_{zz}_particle", etaBins, zetamin, zetamax); listOfTH1.push_back(hzzdeltaEtaparticle);
    // R
      //TH1F *hz1Rparticle = new TH1F("z1_R_particle", "R_{z1}_particle", RBins, zRmin, zRmax); listOfTH1.push_back(hz1Rparticle);
      //TH1F *hz2Rparticle = new TH1F("z2_R_particle", "R_{z2}_particle", RBins, zRmin, zRmax); listOfTH1.push_back(hz2Rparticle);
      TH1F *hzzdeltaRparticle = new TH1F("zz_#DeltaR_particle", "#DeltaR_{zz}_particle", RBins, zRmin, zRmax); listOfTH1.push_back(hzzdeltaRparticle);
  
  // z - parton
    // pT + m
      TH1F *hz1pTparton = new TH1F("z1_pT_parton", "p^{T}_{z1}_parton", pTBins, zpTmin, zpTmax); listOfTH1.push_back(hz1pTparton);
      TH1F *hz2pTparton = new TH1F("z2_pT_parton", "p^{T}_{z2}_parton", pTBins, zpTmin, zpTmax); listOfTH1.push_back(hz2pTparton);
      TH1F *hz1mparton = new TH1F("z1_m_parton", "m_{z1}_parton", mBins, zmmin, zmmax); listOfTH1.push_back(hz1mparton);
      TH1F *hz2mparton = new TH1F("z2_m_parton", "m_{z2}_parton", mBins, zmmin, zmmax); listOfTH1.push_back(hz2mparton);
    // phi
      //TH1F *hz1phiparton = new TH1F("z1_#phi_parton", "#phi_{z1}_parton", phiBins, -TMath::Pi(), +TMath::Pi()); listOfTH1.push_back(hz1phiparton);
      //TH1F *hz2phiparton = new TH1F("z2_#phi_parton", "#phi_{z2}_parton", phiBins, -TMath::Pi(), +TMath::Pi()); listOfTH1.push_back(hz2phiparton);
      TH1F *hzzdeltaPhiparton = new TH1F("zz_#Delta#phi_parton", "#Delta#phi_{zz}_parton", phiBins, -TMath::Pi(), +TMath::Pi()); listOfTH1.push_back(hzzdeltaPhiparton);
    // eta
      //TH1F *hz1etaparton = new TH1F("z1_#eta_parton", "#eta_{z1}_parton", etaBins, zetamin, zetamax); listOfTH1.push_back(hz1etaparton);
      //TH1F *hz2etaparton = new TH1F("z2_#eta_parton", "#eta_{z2}_parton", etaBins, zetamin, zetamax); listOfTH1.push_back(hz2etaparton);
      TH1F *hzzdeltaEtaparton = new TH1F("zz_#Delta#eta_parton", "#Delta#eta_{zz}_parton", etaBins, zetamin, zetamax); listOfTH1.push_back(hzzdeltaEtaparton);
    // R
      //TH1F *hz1Rparton = new TH1F("z1_R_parton", "R_{z1}_parton", RBins, zRmin, zRmax); listOfTH1.push_back(hz1Rparton);
      //TH1F *hz2Rparton = new TH1F("z2_R_parton", "R_{z2}_parton", RBins, zRmin, zRmax); listOfTH1.push_back(hz2Rparton);
      TH1F *hzzdeltaRparton = new TH1F("zz_#DeltaR_parton", "#DeltaR_{zz}_parton", RBins, zRmin, zRmax); listOfTH1.push_back(hzzdeltaRparton);

  // lepton - reco
    // pT
      //TH1F *hl1pTreco = new TH1F("l1_pT_reco", "p^{T}_{l1}_reco", pTBins, lpTmin, lpTmax); listOfTH1.push_back(hl1pTreco);
      //TH1F *hl2pTreco = new TH1F("l2_pT_reco", "p^{T}_{l2}_reco", pTBins, lpTmin, lpTmax); listOfTH1.push_back(hl2pTreco);
      //TH1F *hl3pTreco = new TH1F("l3_pT_reco", "p^{T}_{l3}_reco", pTBins, lpTmin, lpTmax); listOfTH1.push_back(hl3pTreco);
      //TH1F *hl4pTreco = new TH1F("l4_pT_reco", "p^{T}_{l4}_reco", pTBins, lpTmin, lpTmax);  listOfTH1.push_back(hl4pTreco);
    // phi
      //TH1F *hl1phireco = new TH1F("l1_#phi_reco", "#phi_{l1}_reco", phiBins, -TMath::Pi(), +TMath::Pi()); listOfTH1.push_back(hl1phireco);
      //TH1F *hl2phireco = new TH1F("l2_#phi_reco", "#phi_{l2}_reco", phiBins, -TMath::Pi(), +TMath::Pi()); listOfTH1.push_back(hl2phireco);
      //TH1F *hl3phireco = new TH1F("l3_#phi_reco", "#phi_{l3}_reco", phiBins, -TMath::Pi(), +TMath::Pi()); listOfTH1.push_back(hl3phireco);
      //TH1F *hl4phireco = new TH1F("l4_#phi_reco", "#phi_{l4}_reco", phiBins, -TMath::Pi(), +TMath::Pi()); listOfTH1.push_back(hl4phireco);
      TH1F *hl1l2deltaPhireco = new TH1F("l1l2_#Delta#phi_reco", "#Delta#phi_{l1l2}_reco", phiBins, -TMath::Pi(), +TMath::Pi()); listOfTH1.push_back(hl1l2deltaPhireco);
      TH1F *hl3l4deltaPhireco = new TH1F("l3l4_#Delta#phi_reco", "#Delta#phi_{l3l4}_reco", phiBins, -TMath::Pi(), +TMath::Pi()); listOfTH1.push_back(hl3l4deltaPhireco);
      TH1F *hl1l2deltaPhiBoostreco = new TH1F("l1l2_#Delta#phi_Boost_reco", "#Delta#phi_{l1l2}_reco", phiBins, -TMath::Pi(), +TMath::Pi()); listOfTH1.push_back(hl1l2deltaPhiBoostreco);
      TH1F *hl3l4deltaPhiBoostreco = new TH1F("l3l4_#Delta#phi_Boost_reco", "#Delta#phi_{l3l4}_reco", phiBins, -TMath::Pi(), +TMath::Pi()); listOfTH1.push_back(hl3l4deltaPhiBoostreco);
    // eta
      //TH1F *hl1etareco = new TH1F("l1_#eta_reco", "#eta_{l1}_reco", etaBins, letamin, letamax);listOfTH1.push_back(hl1etareco);
      //TH1F *hl2etareco = new TH1F("l2_#eta_reco", "#eta_{l2}_reco", etaBins, letamin, letamax); listOfTH1.push_back(hl2etareco);
      //TH1F *hl3etareco = new TH1F("l3_#eta_reco", "#eta_{l3}_reco", etaBins, letamin, letamax); listOfTH1.push_back(hl3etareco);
      //TH1F *hl4etareco = new TH1F("l4_#eta_reco", "#eta_{l4}_reco", etaBins, letamin, letamax); listOfTH1.push_back(hl4etareco);
      TH1F *hl1l2deltaEtareco = new TH1F("l1l2_#Delta#eta_reco", "#Delta#eta_{l1l2}_reco", etaBins, letamin, letamax); listOfTH1.push_back(hl1l2deltaEtareco);
      TH1F *hl3l4deltaEtareco = new TH1F("l3l4_#Delta#eta_reco", "#Delta#eta_{l3l4}_reco", etaBins, letamin, letamax); listOfTH1.push_back(hl3l4deltaEtareco);
      TH1F *hl1l2deltaEtaBoostreco = new TH1F("l1l2_#Delta#eta_Boost_reco", "#Delta#eta_{l1l2}_Boost_reco", etaBins, letamin, letamax); listOfTH1.push_back(hl1l2deltaEtaBoostreco);
      TH1F *hl3l4deltaEtaBoostreco = new TH1F("l3l4_#Delta#eta_Boost_reco", "#Delta#eta_{l3l4}_Boost_reco", etaBins, letamin, letamax); listOfTH1.push_back(hl3l4deltaEtaBoostreco);
    // R
      //TH1F *hl1Rreco = new TH1F("l1_R_reco", "R_{l1}_reco", RBins, lRmin, lRmax); listOfTH1.push_back(hl1Rreco);
      //TH1F *hl2Rreco = new TH1F("l2_R_reco", "R_{l2}_reco", RBins, lRmin, lRmax); listOfTH1.push_back(hl2Rreco);
      //TH1F *hl3Rreco = new TH1F("l3_R_reco", "R_{l3}_reco", RBins, lRmin, lRmax); listOfTH1.push_back(hl3Rreco);
      //TH1F *hl4Rreco = new TH1F("l4_R_reco", "R_{l4}_reco", RBins, lRmin, lRmax); listOfTH1.push_back(hl4Rreco);
      TH1F *hl1l2deltaRreco = new TH1F("l1l2_#Delta R_reco", "#Delta R_{l1l2}_reco", RBins, lRmin, lRmax); listOfTH1.push_back(hl1l2deltaRreco);
      TH1F *hl3l4deltaRreco = new TH1F("l3l4_#Delta R_reco", "#Delta R_{l3l4}_reco", RBins, lRmin, lRmax); listOfTH1.push_back(hl3l4deltaRreco);
    // cos
      TH1F *hl1cosThetareco = new TH1F("l1_cos#theta_reco", "cos#theta_{l1}_reco", cosBins, -1, 1); listOfTH1.push_back(hl1cosThetareco);
      TH1F *hl2cosThetareco = new TH1F("l2_cos#theta_reco", "cos#theta_{l2}_reco", cosBins, -1, 1); listOfTH1.push_back(hl2cosThetareco);
      TH1F *hl3cosThetareco = new TH1F("l3_cos#theta_reco", "cos#theta_{l3}_reco", cosBins, -1, 1); listOfTH1.push_back(hl3cosThetareco);
      TH1F *hl4cosThetareco = new TH1F("l4_cos#theta_reco", "cos#theta_{l4}_reco", cosBins, -1, 1); listOfTH1.push_back(hl4cosThetareco);
      TH1F *hfourlcosThetareco = new TH1F("fourl_cos#theta_reco", "cos#theta_{fourl}_reco", cosBins, -1, 1); listOfTH1.push_back(hfourlcosThetareco);
      TH1F *hl1cosThetaBoostreco = new TH1F("l1_cos#theta_Boost_reco", "cos#theta_{l1}_Boost_reco", cosBins, -1, 1); listOfTH1.push_back(hl1cosThetaBoostreco);
      TH1F *hl2cosThetaBoostreco = new TH1F("l2_cos#theta_Boost_reco", "cos#theta_{l2}_Boost_reco", cosBins, -1, 1); listOfTH1.push_back(hl2cosThetaBoostreco);
      TH1F *hl3cosThetaBoostreco = new TH1F("l3_cos#theta_Boost_reco", "cos#theta_{l3}_Boost_reco", cosBins, -1, 1); listOfTH1.push_back(hl3cosThetaBoostreco);
      TH1F *hl4cosThetaBoostreco = new TH1F("l4_cos#theta_Boost_reco", "cos#theta_{l4}_Boost_reco", cosBins, -1, 1); listOfTH1.push_back(hl4cosThetaBoostreco);
      TH1F *hfourlcosThetaBoostreco = new TH1F("fourl_cos#theta_Boost_reco", "cos#theta_{fourl}_Boost_reco", cosBins, -1, 1); listOfTH1.push_back(hfourlcosThetaBoostreco);
      TH1F *hl1l2CScosThetareco = new TH1F("l1l2_cos#theta_{CS}_reco", "cos#theta_{CSl1l2}_reco", cosBins, -1, 1); listOfTH1.push_back(hl1l2CScosThetareco);
      TH1F *hl3l4CScosThetareco = new TH1F("l3l4_cos#theta_{CS}_reco", "cos#theta_{CSl3l4}_reco", cosBins, -1, 1); listOfTH1.push_back(hl3l4CScosThetareco);

  // lepton - particle
    // pT
      //TH1F *hl1pTparticle = new TH1F("l1_pT_particle", "p^{T}_{l1}_particle", pTBins, lpTmin, lpTmax); listOfTH1.push_back(hl1pTparticle);
      //TH1F *hl2pTparticle = new TH1F("l2_pT_particle", "p^{T}_{l2}_particle", pTBins, lpTmin, lpTmax); listOfTH1.push_back(hl2pTparticle);
      //TH1F *hl3pTparticle = new TH1F("l3_pT_particle", "p^{T}_{l3}_particle", pTBins, lpTmin, lpTmax); listOfTH1.push_back(hl3pTparticle);
      //TH1F *hl4pTparticle = new TH1F("l4_pT_particle", "p^{T}_{l4}_particle", pTBins, lpTmin, lpTmax); listOfTH1.push_back(hl4pTparticle);
    // phi
      //TH1F *hl1phiparticle = new TH1F("l1_#phi_particle", "#phi_{l1}_particle", phiBins, -TMath::Pi(), +TMath::Pi()); listOfTH1.push_back(hl1phiparticle);
      //TH1F *hl2phiparticle = new TH1F("l2_#phi_particle", "#phi_{l2}_particle", phiBins, -TMath::Pi(), +TMath::Pi()); listOfTH1.push_back(hl2phiparticle);
      //TH1F *hl3phiparticle = new TH1F("l3_#phi_particle", "#phi_{l3}_particle", phiBins, -TMath::Pi(), +TMath::Pi()); listOfTH1.push_back(hl3phiparticle);
      //TH1F *hl4phiparticle = new TH1F("l4_#phi_particle", "#phi_{l4}_particle", phiBins, -TMath::Pi(), +TMath::Pi()); listOfTH1.push_back(hl4phiparticle);
      TH1F *hl1l2deltaPhiparticle = new TH1F("l1l2_#Delta#phi_particle", "#Delta#phi_{l1l2}_particle", phiBins, -TMath::Pi(), +TMath::Pi()); listOfTH1.push_back(hl1l2deltaPhiparticle);
      TH1F *hl3l4deltaPhiparticle = new TH1F("l3l4_#Delta#phi_particle", "#Delta#phi_{l3l4}_particle", phiBins, -TMath::Pi(), +TMath::Pi()); listOfTH1.push_back(hl3l4deltaPhiparticle);
      TH1F *hl1l2deltaPhiBoostparticle = new TH1F("l1l2_#Delta#phi_Boost_particle", "#Delta#phi_{l1l2}_particle", phiBins, -TMath::Pi(), +TMath::Pi()); listOfTH1.push_back(hl1l2deltaPhiBoostparticle);
      TH1F *hl3l4deltaPhiBoostparticle = new TH1F("l3l4_#Delta#phi_Boost_particle", "#Delta#phi_{l3l4}_particle", phiBins, -TMath::Pi(), +TMath::Pi()); listOfTH1.push_back(hl3l4deltaPhiBoostparticle);
    // eta
      //TH1F *hl1etaparticle = new TH1F("l1_#eta_particle", "#eta_{l1}_particle", etaBins, letamin, letamax); listOfTH1.push_back(hl1etaparticle);
      //TH1F *hl2etaparticle = new TH1F("l2_#eta_particle", "#eta_{l2}_particle", etaBins, letamin, letamax); listOfTH1.push_back(hl2etaparticle);
      //TH1F *hl3etaparticle = new TH1F("l3_#eta_particle", "#eta_{l3}_particle", etaBins, letamin, letamax); listOfTH1.push_back(hl3etaparticle);
      //TH1F *hl4etaparticle = new TH1F("l4_#eta_particle", "#eta_{l4}_particle", etaBins, letamin, letamax); listOfTH1.push_back(hl4etaparticle);
      TH1F *hl1l2deltaEtaparticle = new TH1F("l1l2_#Delta#eta_particle", "#Delta#eta_{l1l2}_particle", etaBins, letamin, letamax); listOfTH1.push_back(hl1l2deltaEtaparticle);
      TH1F *hl3l4deltaEtaparticle = new TH1F("l3l4_#Delta#eta_particle", "#Delta#eta_{l3l4}_particle", etaBins, letamin, letamax); listOfTH1.push_back(hl3l4deltaEtaparticle);
      TH1F *hl1l2deltaEtaBoostparticle = new TH1F("l1l2_#Delta#eta_Boost_particle", "#Delta#eta_{l1l2}_Boost_particle", etaBins, letamin, letamax); listOfTH1.push_back(hl1l2deltaEtaBoostparticle);
      TH1F *hl3l4deltaEtaBoostparticle = new TH1F("l3l4_#Delta#eta_Boost_particle", "#Delta#eta_{l3l4}_Boost_particle", etaBins, letamin, letamax);listOfTH1.push_back(hl3l4deltaEtaBoostparticle);
    // R
      //TH1F *hl1Rparticle = new TH1F("l1_R_particle", "R_{l1}_particle", RBins, lRmin, lRmax); listOfTH1.push_back(hl1Rparticle);
      //TH1F *hl2Rparticle = new TH1F("l2_R_particle", "R_{l2}_particle", RBins, lRmin, lRmax); listOfTH1.push_back(hl2Rparticle);
      //TH1F *hl3Rparticle = new TH1F("l3_R_particle", "R_{l3}_particle", RBins, lRmin, lRmax); listOfTH1.push_back(hl3Rparticle);
      //TH1F *hl4Rparticle = new TH1F("l4_R_particle", "R_{l4}_particle", RBins, lRmin, lRmax); listOfTH1.push_back(hl4Rparticle);
      TH1F *hl1l2deltaRparticle = new TH1F("l1l2_#Delta R_particle", "#Delta R_{l1l2}_particle", RBins, lRmin, lRmax); listOfTH1.push_back(hl1l2deltaRparticle);
      TH1F *hl3l4deltaRparticle = new TH1F("l3l4_#Delta R_particle", "#Delta R_{l3l4}_particle", RBins, lRmin, lRmax); listOfTH1.push_back(hl3l4deltaRparticle);
    // cos 
      TH1F *hl1cosThetaparticle = new TH1F("l1_cos#theta_particle", "cos#theta_{l1}_particle", cosBins, -1, 1); listOfTH1.push_back(hl1cosThetaparticle);
      TH1F *hl2cosThetaparticle = new TH1F("l2_cos#theta_particle", "cos#theta_{l2}_particle", cosBins, -1, 1); listOfTH1.push_back(hl2cosThetaparticle);
      TH1F *hl3cosThetaparticle = new TH1F("l3_cos#theta_particle", "cos#theta_{l3}_particle", cosBins, -1, 1); listOfTH1.push_back(hl3cosThetaparticle);
      TH1F *hl4cosThetaparticle = new TH1F("l4_cos#theta_particle", "cos#theta_{l4}_particle", cosBins, -1, 1); listOfTH1.push_back(hl4cosThetaparticle);
      TH1F *hfourlcosThetaparticle = new TH1F("fourl_cos#theta_particle", "cos#theta_{fourl}_particle", cosBins, -1, 1); listOfTH1.push_back(hfourlcosThetaparticle);
      TH1F *hl1cosThetaBoostparticle = new TH1F("l1_cos#theta_Boost_particle", "cos#theta_{l1}_Boost_particle", cosBins, -1, 1); listOfTH1.push_back(hl1cosThetaBoostparticle);
      TH1F *hl2cosThetaBoostparticle = new TH1F("l2_cos#theta_Boost_particle", "cos#theta_{l2}_Boost_particle", cosBins, -1, 1); listOfTH1.push_back(hl2cosThetaBoostparticle);
      TH1F *hl3cosThetaBoostparticle = new TH1F("l3_cos#theta_Boost_particle", "cos#theta_{l3}_Boost_particle", cosBins, -1, 1); listOfTH1.push_back(hl3cosThetaBoostparticle);
      TH1F *hl4cosThetaBoostparticle = new TH1F("l4_cos#theta_Boost_particle", "cos#theta_{l4}_Boost_particle", cosBins, -1, 1); listOfTH1.push_back(hl4cosThetaBoostparticle);
      TH1F *hfourlcosThetaBoostparticle = new TH1F("fourl_cos#theta_Boost_particle", "cos#theta_{fourl}_Boost_particle", cosBins, -1, 1); listOfTH1.push_back(hfourlcosThetaBoostparticle);
      TH1F *hl1l2CScosThetaparticle = new TH1F("l1l2_cos#theta_{CS}_particle", "cos#theta_{CSl1l2}_particle", cosBins, -1, 1); listOfTH1.push_back(hl1l2CScosThetaparticle);
      TH1F *hl3l4CScosThetaparticle = new TH1F("l3l4_cos#theta_{CS}_particle", "cos#theta_{CSl3l4}_particle", cosBins, -1, 1); listOfTH1.push_back(hl3l4CScosThetaparticle);

  // lepton - parton
    // pT
      //TH1F *hl1pTparton = new TH1F("l1_pT_parton", "p^{T}_{l1}_parton", pTBins, lpTmin, lpTmax); listOfTH1.push_back(hl1pTparton);
      //TH1F *hl2pTparton = new TH1F("l2_pT_parton", "p^{T}_{l2}_parton", pTBins, lpTmin, lpTmax); listOfTH1.push_back(hl2pTparton);
      //TH1F *hl3pTparton = new TH1F("l3_pT_parton", "p^{T}_{l3}_parton", pTBins, lpTmin, lpTmax);  listOfTH1.push_back(hl3pTparton);
      //TH1F *hl4pTparton = new TH1F("l4_pT_parton", "p^{T}_{l4}_parton", pTBins, lpTmin, lpTmax); listOfTH1.push_back(hl4pTparton);
    // phi
      //TH1F *hl1phiparton = new TH1F("l1_#phi_parton", "#phi_{l1}_parton", phiBins, -TMath::Pi(), +TMath::Pi()); listOfTH1.push_back(hl1phiparton);
      //TH1F *hl2phiparton = new TH1F("l2_#phi_parton", "#phi_{l2}_parton", phiBins, -TMath::Pi(), +TMath::Pi()); listOfTH1.push_back(hl2phiparton);
      //TH1F *hl3phiparton = new TH1F("l3_#phi_parton", "#phi_{l3}_parton", phiBins, -TMath::Pi(), +TMath::Pi()); listOfTH1.push_back(hl3phiparton);
      //TH1F *hl4phiparton = new TH1F("l4_#phi_parton", "#phi_{l4}_parton", phiBins, -TMath::Pi(), +TMath::Pi()); listOfTH1.push_back(hl4phiparton);
      TH1F *hl1l2deltaPhiparton = new TH1F("l1l2_#Delta#phi_parton", "#Delta#phi_{l1l2}_parton", phiBins, -TMath::Pi(), +TMath::Pi()); listOfTH1.push_back(hl1l2deltaPhiparton);
      TH1F *hl3l4deltaPhiparton = new TH1F("l3l4_#Delta#phi_parton", "#Delta#phi_{l3l4}_parton", phiBins, -TMath::Pi(), +TMath::Pi()); listOfTH1.push_back(hl3l4deltaPhiparton);
      TH1F *hl1l2deltaPhiBoostparton = new TH1F("l1l2_#Delta#phi_Boost_parton", "#Delta#phi_{l1l2}_parton", phiBins, -TMath::Pi(), +TMath::Pi()); listOfTH1.push_back(hl1l2deltaPhiBoostparton);
      TH1F *hl3l4deltaPhiBoostparton = new TH1F("l3l4_#Delta#phi_Boost_parton", "#Delta#phi_{l3l4}_parton", phiBins, -TMath::Pi(), +TMath::Pi()); listOfTH1.push_back(hl3l4deltaPhiBoostparton);
    // eta
      //TH1F *hl1etaparton = new TH1F("l1_#eta_parton", "#eta_{l1}_parton", etaBins, letamin, letamax); listOfTH1.push_back(hl1etaparton);
      //TH1F *hl2etaparton = new TH1F("l2_#eta_parton", "#eta_{l2}_parton", etaBins, letamin, letamax); listOfTH1.push_back(hl2etaparton);
      //TH1F *hl3etaparton = new TH1F("l3_#eta_parton", "#eta_{l3}_parton", etaBins, letamin, letamax); listOfTH1.push_back(hl3etaparton);
      //TH1F *hl4etaparton = new TH1F("l4_#eta_parton", "#eta_{l4}_parton", etaBins, letamin, letamax); listOfTH1.push_back(hl4etaparton);
      TH1F *hl1l2deltaEtaparton = new TH1F("l1l2_#Delta#eta_parton", "#Delta#eta_{l1l2}_parton", etaBins, letamin, letamax); listOfTH1.push_back(hl1l2deltaEtaparton);
      TH1F *hl3l4deltaEtaparton = new TH1F("l3l4_#Delta#eta_parton", "#Delta#eta_{l3l4}_parton", etaBins, letamin, letamax); listOfTH1.push_back(hl3l4deltaEtaparton);
      TH1F *hl1l2deltaEtaBoostparton = new TH1F("l1l2_#Delta#eta_Boost_parton", "#Delta#eta_{l1l2}_Boost_parton", etaBins, letamin, letamax); listOfTH1.push_back(hl1l2deltaEtaBoostparton);
      TH1F *hl3l4deltaEtaBoostparton = new TH1F("l3l4_#Delta#eta_Boost_parton", "#Delta#eta_{l3l4}_Boost_parton", etaBins, letamin, letamax); listOfTH1.push_back(hl3l4deltaEtaBoostparton);
    // R
      //TH1F *hl1Rparton = new TH1F("l1_R_parton", "R_{l1}_parton", RBins, lRmin, lRmax); listOfTH1.push_back(hl1Rparton);
      //TH1F *hl2Rparton = new TH1F("l2_R_parton", "R_{l2}_parton", RBins, lRmin, lRmax); listOfTH1.push_back(hl2Rparton);
      //TH1F *hl3Rparton = new TH1F("l3_R_parton", "R_{l3}_parton", RBins, lRmin, lRmax); listOfTH1.push_back(hl3Rparton);
      //TH1F *hl4Rparton = new TH1F("l4_R_parton", "R_{l4}_parton", RBins, lRmin, lRmax); listOfTH1.push_back(hl4Rparton);
      TH1F *hl1l2deltaRparton = new TH1F("l1l2_#Delta R_parton", "#Delta R_{l1l2}_parton", RBins, lRmin, lRmax); listOfTH1.push_back(hl1l2deltaRparton);
      TH1F *hl3l4deltaRparton = new TH1F("l3l4_#Delta R_parton", "#Delta R_{l3l4}_parton", RBins, lRmin, lRmax); listOfTH1.push_back(hl3l4deltaRparton);
    // cos
      TH1F *hl1cosThetaparton = new TH1F("l1_cos#theta_parton", "cos#theta_{l1}_parton", cosBins, -1, 1); listOfTH1.push_back(hl1cosThetaparton);
      TH1F *hl2cosThetaparton = new TH1F("l2_cos#theta_parton", "cos#theta_{l2}_parton", cosBins, -1, 1); listOfTH1.push_back(hl2cosThetaparton);
      TH1F *hl3cosThetaparton = new TH1F("l3_cos#theta_parton", "cos#theta_{l3}_parton", cosBins, -1, 1); listOfTH1.push_back(hl3cosThetaparton);
      TH1F *hl4cosThetaparton = new TH1F("l4_cos#theta_parton", "cos#theta_{l4}_parton", cosBins, -1, 1); listOfTH1.push_back(hl4cosThetaparton);
      TH1F *hfourlcosThetaparton = new TH1F("fourl_cos#theta_parton", "cos#theta_{fourl}_parton", cosBins, -1, 1); listOfTH1.push_back(hfourlcosThetaparton);
      TH1F *hl1cosThetaBoostparton = new TH1F("l1_cos#theta_Boost_parton", "cos#theta_{l1}_Boost_parton", cosBins, -1, 1); listOfTH1.push_back(hl1cosThetaBoostparton);
      TH1F *hl2cosThetaBoostparton = new TH1F("l2_cos#theta_Boost_parton", "cos#theta_{l2}_Boost_parton", cosBins, -1, 1); listOfTH1.push_back(hl2cosThetaBoostparton);
      TH1F *hl3cosThetaBoostparton = new TH1F("l3_cos#theta_Boost_parton", "cos#theta_{l3}_Boost_parton", cosBins, -1, 1); listOfTH1.push_back(hl3cosThetaBoostparton);
      TH1F *hl4cosThetaBoostparton = new TH1F("l4_cos#theta_Boost_parton", "cos#theta_{l4}_Boost_parton", cosBins, -1, 1); listOfTH1.push_back(hl4cosThetaBoostparton);
      TH1F *hfourlcosThetaBoostparton = new TH1F("fourl_cos#theta_Boost_parton", "cos#theta_{fourl}_Boost_parton", cosBins, -1, 1); listOfTH1.push_back(hfourlcosThetaBoostparton);
      TH1F *hl1l2CScosThetaparton = new TH1F("l1l2_cos#theta_{CS}_parton", "cos#theta_{CSl1l2}_parton", cosBins, -1, 1); listOfTH1.push_back(hl1l2CScosThetaparton);
      TH1F *hl3l4CScosThetaparton = new TH1F("l3l4_cos#theta_{CS}_parton", "cos#theta_{CSl3l4}_parton", cosBins, -1, 1); listOfTH1.push_back(hl3l4CScosThetaparton); 

// 2D comp - parton(1) particle(2) reco(3)

  // higgs comps
    // pT + m
      TH2F *hHpT12Comp = new TH2F("H_pT_comp_12", "p_{T}^{hbb}", pTBins, hpTmin, hpTmax, pTBins, hpTmin, hpTmax); listOfTH2.push_back(hHpT12Comp);
      TH2F *hHpT23Comp = new TH2F("H_pT_comp_23", "p_{T}^{hbb}", pTBins, hpTmin, hpTmax, pTBins, hpTmin, hpTmax); listOfTH2.push_back(hHpT23Comp);
      TH2F *hHpT13Comp = new TH2F("H_pT_comp_13", "p_{T}^{hbb}", pTBins, hpTmin, hpTmax, pTBins, hpTmin, hpTmax); listOfTH2.push_back(hHpT13Comp);
      TH2F *hHm12Comp = new TH2F("H_m_comp_12", "m_{hbb}", mBins, hmmin, hmmax, mBins, hmmin, hmmax); listOfTH2.push_back(hHm12Comp);
      TH2F *hHm23Comp = new TH2F("H_m_comp_23", "m_{hbb}", mBins, hmmin, hmmax, mBins, hmmin, hmmax); listOfTH2.push_back(hHm23Comp);      
      TH2F *hHm13Comp = new TH2F("H_m_comp_13", "m_{hbb}", mBins, hmmin, hmmax, mBins, hmmin, hmmax); listOfTH2.push_back(hHm13Comp);
    // phi
      TH2F*hbbdeltaPhi12Comp = new TH2F("bb_#Delta#phi_comp_12", "#Delta#phi_{bb}", phiBins, -TMath::Pi(),+TMath::Pi(), phiBins, -TMath::Pi(),+TMath::Pi()); listOfTH2.push_back(hbbdeltaPhi12Comp);
      TH2F*hbbdeltaPhi23Comp = new TH2F("bb_#Delta#phi_comp_23", "#Delta#phi_{bb}", phiBins, -TMath::Pi(),+TMath::Pi(), phiBins, -TMath::Pi(),+TMath::Pi()); listOfTH2.push_back(hbbdeltaPhi23Comp);
      TH2F*hbbdeltaPhi13Comp = new TH2F("bb_#Delta#phi_comp_13", "#Delta#phi_{bb}", phiBins, -TMath::Pi(),+TMath::Pi(), phiBins, -TMath::Pi(),+TMath::Pi()); listOfTH2.push_back(hbbdeltaPhi13Comp);
    // eta  
      TH2F*hbbdeltaEta12Comp = new TH2F("bb_#Delta#eta_comp_12", "#Delta#eta_{bb}", etaBins, hetamin, hetamax, etaBins, hetamin, hetamax); listOfTH2.push_back(hbbdeltaEta12Comp);
      TH2F*hbbdeltaEta23Comp = new TH2F("bb_#Delta#eta_comp_23", "#Delta#eta_{bb}", etaBins, hetamin, hetamax, etaBins, hetamin, hetamax); listOfTH2.push_back(hbbdeltaEta23Comp);
      TH2F*hbbdeltaEta13Comp = new TH2F("bb_#Delta#eta_comp_13", "#Delta#eta_{bb}", etaBins, hetamin, hetamax, etaBins, hetamin, hetamax); listOfTH2.push_back(hbbdeltaEta13Comp);
    
  // vbfj comps
    // pT  
      TH2F*hjjpT12Comp = new TH2F("jj_pT_comp_12", "p_{T}^{jj}", pTBins, jpTmin, jpTmax, pTBins, jpTmin, jpTmax); listOfTH2.push_back(hjjpT12Comp);
      TH2F*hjjpT23Comp = new TH2F("jj_pT_comp_23", "p_{T}^{jj}", pTBins, jpTmin, jpTmax, pTBins, jpTmin, jpTmax); listOfTH2.push_back(hjjpT23Comp);
      TH2F*hjjpT13Comp = new TH2F("jj_pT_comp_13", "p_{T}^{jj}", pTBins, jpTmin, jpTmax, pTBins, jpTmin, jpTmax); listOfTH2.push_back(hjjpT13Comp);
    // phi 
      //TH2F*hj1Phi23Comp = new TH2F("j1_#phi_comp_23", "#phi_{j1}", phiBins, -TMath::Pi(),+TMath::Pi(), phiBins, -TMath::Pi(),+TMath::Pi()); listOfTH2.push_back(hj1Phi23Comp);
      //TH2F*hj2Phi23Comp = new TH2F("j2_#phi_comp_23", "#phi_{j2}", phiBins, -TMath::Pi(),+TMath::Pi(), phiBins, -TMath::Pi(),+TMath::Pi()); listOfTH2.push_back(hj2Phi23Comp);
      TH2F*hjjdeltaPhi12Comp = new TH2F("jj_#Delta#phi_comp_12", "#Delta#phi_{jj}", phiBins, -TMath::Pi(),+TMath::Pi(), phiBins, -TMath::Pi(),+TMath::Pi()); listOfTH2.push_back(hjjdeltaPhi12Comp);
      TH2F*hjjdeltaPhi23Comp = new TH2F("jj_#Delta#phi_comp_23", "#Delta#phi_{jj}", phiBins, -TMath::Pi(),+TMath::Pi(), phiBins, -TMath::Pi(),+TMath::Pi()); listOfTH2.push_back(hjjdeltaPhi23Comp);
      TH2F*hjjdeltaPhi13Comp = new TH2F("jj_#Delta#phi_comp_13", "#Delta#phi_{jj}", phiBins, -TMath::Pi(),+TMath::Pi(), phiBins, -TMath::Pi(),+TMath::Pi()); listOfTH2.push_back(hjjdeltaPhi13Comp);
    
  // z comps
    // pT + m
      TH2F*hz1pT12Comp = new TH2F("z1_pT_comp_12", "p^{T}_{z1}", pTBins, zpTmin, zpTmax, pTBins, zpTmin, zpTmax); listOfTH2.push_back(hz1pT12Comp);
      TH2F*hz1pT23Comp = new TH2F("z1_pT_comp_23", "p^{T}_{z1}", pTBins, zpTmin, zpTmax, pTBins, zpTmin, zpTmax); listOfTH2.push_back(hz1pT23Comp);
      TH2F*hz1pT13Comp = new TH2F("z1_pT_comp_13", "p^{T}_{z1}", pTBins, zpTmin, zpTmax, pTBins, zpTmin, zpTmax); listOfTH2.push_back(hz1pT13Comp);
      TH2F*hz2pT12Comp = new TH2F("z2_pT_comp_12", "p^{T}_{z2}", pTBins, zpTmin, zpTmax, pTBins, zpTmin, zpTmax); listOfTH2.push_back(hz2pT12Comp);
      TH2F*hz2pT23Comp = new TH2F("z2_pT_comp_23", "p^{T}_{z2}", pTBins, zpTmin, zpTmax, pTBins, zpTmin, zpTmax); listOfTH2.push_back(hz2pT23Comp);
      TH2F*hz2pT13Comp = new TH2F("z2_pT_comp_13", "p^{T}_{z2}", pTBins, zpTmin, zpTmax, pTBins, zpTmin, zpTmax); listOfTH2.push_back(hz2pT13Comp);
      TH2F*hz1m12Comp = new TH2F("z1_m_comp_12", "m_{z1}", mBins, zmmin, zmmax, mBins, zmmin, zmmax); listOfTH2.push_back(hz1m12Comp);
      TH2F*hz1m23Comp = new TH2F("z1_m_comp_23", "m_{z1}", mBins, zmmin, zmmax, mBins, zmmin, zmmax); listOfTH2.push_back(hz1m23Comp);      
      TH2F*hz1m13Comp = new TH2F("z1_m_comp_13", "m_{z1}", mBins, zmmin, zmmax, mBins, zmmin, zmmax); listOfTH2.push_back(hz1m13Comp);
      TH2F*hz2m12Comp = new TH2F("z2_m_comp_12", "m_{z2}", mBins, zmmin, zmmax, mBins, zmmin, zmmax); listOfTH2.push_back(hz2m12Comp);
      TH2F*hz2m23Comp = new TH2F("z2_m_comp_23", "m_{z2}", mBins, zmmin, zmmax, mBins, zmmin, zmmax); listOfTH2.push_back(hz2m23Comp);
      TH2F*hz2m13Comp = new TH2F("z2_m_comp_13", "m_{z2}", mBins, zmmin, zmmax, mBins, zmmin, zmmax); listOfTH2.push_back(hz2m13Comp);
    

  // mixed comps
    
    // higgs dphi vs. vbfj dphi
      TH2F*hbbjjdeltaPhicompreco = new TH2F("bb_jj_#Delta#phi_comp_reco", "bb_jj_#Delta#phi_comp_reco", phiBins, -TMath::Pi(), +TMath::Pi(), phiBins, -TMath::Pi(), +TMath::Pi()); listOfTH2.push_back(hbbjjdeltaPhicompreco);
      TH2F*hbbjjdeltaPhicompparticle = new TH2F("bb_jj_#Delta#phi_comp_particle", "bb_jj_#Delta#phi_comp_particle", phiBins, -TMath::Pi(), +TMath::Pi(), phiBins, -TMath::Pi(), +TMath::Pi()); listOfTH2.push_back(hbbjjdeltaPhicompparticle);
      TH2F*hbbjjdeltaPhicompparton = new TH2F("bb_jj_#Delta#phi_comp_parton", "bb_jj_#Delta#phi_comp_parton", phiBins, -TMath::Pi(), +TMath::Pi(), phiBins, -TMath::Pi(), +TMath::Pi()); listOfTH2.push_back(hbbjjdeltaPhicompparton);

    // higgs deta vs. vbfj deta
      TH2F*hbbdeltaEtajjdeltaPhicompreco = new TH2F("bb_#Delta#eta_jj_#Delta#phi_comp_reco", "bb_#Delta#eta_jj_#Delta#phi_comp_reco", etaBins, hetamin, hetamax, phiBins, -TMath::Pi(), +TMath::Pi()); listOfTH2.push_back(hbbdeltaEtajjdeltaPhicompreco);
      TH2F*hbbdeltaEtajjdeltaPhicompparticle = new TH2F("bb_#Delta#eta_jj_#Delta#phi_comp_particle", "bb_#Delta#eta_jj_#Delta#phi_comp_particle", etaBins, hetamin, hetamax, phiBins, -TMath::Pi(), +TMath::Pi()); listOfTH2.push_back(hbbdeltaEtajjdeltaPhicompparticle);
      TH2F*hbbdeltaEtajjdeltaPhicompparton = new TH2F("bb_#Delta#eta_jj_#Delta#phi_comp_parton", "bb_#Delta#eta_jj_#Delta#phi_comp_parton", etaBins, hetamin, hetamax, phiBins, -TMath::Pi(), +TMath::Pi()); listOfTH2.push_back(hbbdeltaEtajjdeltaPhicompparton);

    // higgs pT vs. lepton dphi 
      TH2F*hHpTl1l2deltaPhicompreco = new TH2F("h_pT_l1l2_delta#phi_comp_reco", "h_pT_l1l2_delta#phi_comp_reco", pTBins, hpTmin, hpTmax, phiBins, -TMath::Pi(), +TMath::Pi()); listOfTH2.push_back(hHpTl1l2deltaPhicompreco);
      TH2F*hHpTl3l4deltaPhicompreco = new TH2F("h_pT_l3l4_delta#phi_comp_reco", "h_pT_l1l2_delta#phi_comp_reco", pTBins, hpTmin, hpTmax, phiBins, -TMath::Pi(), +TMath::Pi()); listOfTH2.push_back(hHpTl3l4deltaPhicompreco);
      TH2F*hHpTl1l2deltaPhicompparticle = new TH2F("h_pT_l1l2_delta#phi_comp_particle", "h_pT_l1l2_delta#phi_comp_particle", pTBins, hpTmin, hpTmax, phiBins, -TMath::Pi(), +TMath::Pi()); listOfTH2.push_back(hHpTl1l2deltaPhicompparticle);
      TH2F*hHpTl3l4deltaPhicompparticle = new TH2F("h_pT_l3l4_delta#phi_comp_particle", "h_pT_l1l2_delta#phi_comp_particle", pTBins, hpTmin, hpTmax, phiBins, -TMath::Pi(), +TMath::Pi()); listOfTH2.push_back(hHpTl3l4deltaPhicompparticle);
      TH2F*hHpTl1l2deltaPhicompparton = new TH2F("h_pT_l1l2_delta#phi_comp_parton", "h_pT_l1l2_delta#phi_comp_parton", pTBins, hpTmin, hpTmax, phiBins, -TMath::Pi(), +TMath::Pi()); listOfTH2.push_back(hHpTl1l2deltaPhicompparton);
      TH2F*hHpTl3l4deltaPhicompparton = new TH2F("h_pT_l3l4_delta#phi_comp_parton", "h_pT_l1l2_delta#phi_comp_parton", pTBins, hpTmin, hpTmax, phiBins, -TMath::Pi(), +TMath::Pi()); listOfTH2.push_back(hHpTl3l4deltaPhicompparton);

    // higgs pT vs. z pT
      TH2F*hHz1pTcompreco = new TH2F("h_z1_pT_comp_reco", "h_z1_pT_comp_reco", pTBins, hpTmin, hpTmax, pTBins, zpTmin, zpTmax); listOfTH2.push_back(hHz1pTcompreco);
      TH2F*hHz2pTcompreco = new TH2F("h_z2_pT_comp_reco", "h_z2_pT_comp_reco", pTBins, hpTmin, hpTmax, pTBins, zpTmin, zpTmax); listOfTH2.push_back(hHz2pTcompreco);
      TH2F*hHz1pTcompparticle = new TH2F("h_z1_pT_comp_particle", "h_z1_pT_comp_particle", pTBins, hpTmin, hpTmax, pTBins, zpTmin, zpTmax); listOfTH2.push_back(hHz1pTcompparticle);
      TH2F*hHz2pTcompparticle = new TH2F("h_z2_pT_comp_particle", "h_z2_pT_comp_particle", pTBins, hpTmin, hpTmax, pTBins, zpTmin, zpTmax); listOfTH2.push_back(hHz2pTcompparticle);
      TH2F*hHz1pTcompparton = new TH2F("h_z1_pT_comp_parton", "h_z1_pT_comp_parton", pTBins, hpTmin, hpTmax, pTBins, zpTmin, zpTmax); listOfTH2.push_back(hHz1pTcompparton);
      TH2F*hHz2pTcompparton = new TH2F("h_z2_pT_comp_parton", "h_z2_pT_comp_parton", pTBins, hpTmin, hpTmax, pTBins, zpTmin, zpTmax); listOfTH2.push_back(hHz2pTcompparton);

    // higgs pT vs. zz pT      
      TH2F*hHzzpTcompreco = new TH2F("h_zz_pT_comp_reco", "h_zz_pT_comp_reco", pTBins, hpTmin, 2*hpTmax, pTBins, zpTmin, 2*zpTmax); listOfTH2.push_back(hHzzpTcompreco);
      TH2F*hHzzpTcompparticle = new TH2F("h_zz_pT_comp_particle", "h_zz_pT_comp_particle", pTBins, hpTmin, 2*hpTmax, pTBins, zpTmin, 2*zpTmax); listOfTH2.push_back(hHzzpTcompparticle);
      TH2F*hHzzpTcompparton = new TH2F("h_zz_pT_comp_parton", "h_zz_pT_comp_parton", pTBins, hpTmin, 2*hpTmax, pTBins, zpTmin, 2*zpTmax); listOfTH2.push_back(hHzzpTcompparton);

    // higgs pT vs. zz dphi
      TH2F*hHpTzzdeltaPhicompreco = new TH2F("h_pT_zz_delta#phi_comp_reco", "h_pT_zz_delta#phi_comp_reco", pTBins, hpTmin, hpTmax, phiBins, -TMath::Pi(), +TMath::Pi()); listOfTH2.push_back(hHpTzzdeltaPhicompreco);
      TH2F*hHpTzzdeltaPhicompparticle = new TH2F("h_pT_zz_delta#phi_comp_particle", "h_pT_zz_delta#phi_comp_particle", pTBins, hpTmin, hpTmax, phiBins, -TMath::Pi(), +TMath::Pi()); listOfTH2.push_back(hHpTzzdeltaPhicompparticle);
      TH2F*hHpTzzdeltaPhicompparton = new TH2F("h_pT_zz_delta#phi_comp_parton", "h_pT_zz_delta#phi_comp_parton", pTBins, hpTmin, hpTmax, phiBins, -TMath::Pi(), +TMath::Pi()); listOfTH2.push_back(hHpTzzdeltaPhicompparton);
    
    // higgs pT vs. zz deta    
      TH2F*hHpTzzdeltaEtacompreco = new TH2F("h_pT_zz_delta#eta_comp_reco", "h_pT_zz_delta#eta_comp_reco", pTBins, hpTmin, hpTmax, etaBins, zetamin, zetamax); listOfTH2.push_back(hHpTzzdeltaEtacompreco);
      TH2F*hHpTzzdeltaEtacompparticle = new TH2F("h_pT_zz_delta#eta_comp_particle", "h_pT_zz_delta#eta_comp_particle", pTBins, hpTmin, hpTmax, etaBins, zetamin, zetamax); listOfTH2.push_back(hHpTzzdeltaEtacompparticle);
      TH2F*hHpTzzdeltaEtacompparton = new TH2F("h_pT_zz_delta#eta_comp_parton", "h_pT_zz_delta#eta_comp_parton", pTBins, hpTmin, hpTmax, etaBins, zetamin, zetamax); listOfTH2.push_back(hHpTzzdeltaEtacompparton);

    // higgs pT vs. z eta
      //TH2F*hHpTz1etacompreco = new TH2F("h_pT_z1_eta_comp_reco", "h_pT_z1_eta_comp_reco", pTBins, hpTmin, hpTmax, etaBins, zetamin, zetamax); listOfTH2.push_back(hHpTz1etacompreco);
      //TH2F*hHpTz2etacompreco = new TH2F("h_pT_z2_eta_comp_reco", "h_pT_z2_eta_comp_reco", pTBins, hpTmin, hpTmax, etaBins, zetamin, zetamax); listOfTH2.push_back(hHpTz2etacompreco);
      //TH2F*hHpTz1etacompparticle = new TH2F("h_pT_z1_eta_comp_particle", "h_pT_z1_eta_comp_particle", pTBins, hpTmin, hpTmax, etaBins, zetamin, zetamax); listOfTH2.push_back(hHpTz1etacompparticle);
      //TH2F*hHpTz2etacompparticle = new TH2F("h_pT_z2_eta_comp_particle", "h_pT_z2_eta_comp_particle", pTBins, hpTmin, hpTmax, etaBins, zetamin, zetamax); listOfTH2.push_back(hHpTz2etacompparticle);
      //TH2F*hHpTz1etacompparton = new TH2F("h_pT_z1_eta_comp_parton", "h_pT_z1_eta_comp_parton", pTBins, hpTmin, hpTmax, etaBins, zetamin, zetamax); listOfTH2.push_back(hHpTz1etacompparton);
      //TH2F*hHpTz2etacompparton= new TH2F("h_pT_z2_eta_comp_parton", "h_pT_z2_eta_comp_parton", pTBins, hpTmin, hpTmax, etaBins, zetamin, zetamax); listOfTH2.push_back(hHpTz2etacompparton);
    
    // higgs dphi vs. zz dphi
      TH2F*hbbzzdeltaPhicompreco = new TH2F("bb_zz_#Delta#phi_comp_reco", "bb_zz_#Delta#phi_comp_reco", phiBins, -TMath::Pi(), +TMath::Pi(), phiBins, -TMath::Pi(), +TMath::Pi()); listOfTH2.push_back(hbbzzdeltaPhicompreco);
      TH2F*hbbzzdeltaPhicompparticle = new TH2F("bb_zz_#Delta#phi_comp_particle", "bb_zz_#Delta#phi_comp_particle", phiBins, -TMath::Pi(), +TMath::Pi(), phiBins, -TMath::Pi(), +TMath::Pi()); listOfTH2.push_back(hbbzzdeltaPhicompparticle);
      TH2F*hbbzzdeltaPhicompparton = new TH2F("bb_zz_#Delta#phi_comp_parton", "bb_zz_#Delta#phi_comp_parton", phiBins, -TMath::Pi(), +TMath::Pi(), phiBins, -TMath::Pi(), +TMath::Pi()); listOfTH2.push_back(hbbzzdeltaPhicompparton);

    // higgs deta vs. zz deta  
      TH2F*hbbzzdeltaEtacompreco = new TH2F("bb_zz_delta#eta_comp_reco", "bb_zz_delta#eta_comp_reco", etaBins, hetamin, hetamax, etaBins, zetamin, zetamax); listOfTH2.push_back(hbbzzdeltaEtacompreco);
      TH2F*hbbzzdeltaEtacompparticle = new TH2F("bb_zz_delta#eta_comp_particle", "bb_zz_delta#eta_comp_particle", etaBins, hetamin, hetamax, etaBins, zetamin, zetamax); listOfTH2.push_back(hbbzzdeltaEtacompparticle);
      TH2F*hbbzzdeltaEtacompparton = new TH2F("bb_zz_delta#eta_comp_parton", "bb_zz_delta#eta_comp_parton", etaBins, hetamin, hetamax, etaBins, zetamin, zetamax); listOfTH2.push_back(hbbzzdeltaEtacompparton);

  TProfile *kappaLambda = new TProfile("kappaLambda", "kappaLambda", 40, -20, 20);
  kappaLambda -> GetXaxis() -> SetTitle("#kappa_{#lambda}");

    // b tag
      TH1F *leadbscorereco = new TH1F("lead_bscore_reco", "leading b score reco", 100, 0, 1); listOfTH1.push_back(leadbscorereco);
      TH1F *leadbscoreparticle = new TH1F("lead_bscore_particle", "leading b score particle", 100, 0, 1); listOfTH1.push_back(leadbscoreparticle);

      TH1F *subleadbscorereco = new TH1F("sublead_bscore_reco", "sub-leading b score reco", 100, 0, 1); listOfTH1.push_back(subleadbscorereco);
      TH1F *subleadbscoreparticle = new TH1F("sublead_bscore_particle", "sub-leading b score particle", 100, 0, 1); listOfTH1.push_back(subleadbscoreparticle); 

      //TH2F *leadbscore23 = new TH2F("lead_bscore_comp_23", "leading b score", 5, 0, 1, 5, 0, 1); listOfTH2.push_back(leadbscore23);
      //TH2F *subleadbscore23 = new TH2F("sublead_bscore_comp_23", "sub-leading b score", 5, 0, 1, 5, 0, 1); listOfTH2.push_back(subleadbscore23);

      TH1F *pairedJetsize = new TH1F("pairedJet_size", "pairedJet size", 10, 0, 5); listOfTH1.push_back(pairedJetsize);
      TH1F *pairedJetBsize = new TH1F("pairedJetB_size", "pairedJetB size", 5, 0, 5); listOfTH1.push_back(pairedJetBsize);
      TH2F *pairedJetbbsizemass = new TH2F("pairedJetbb_size_mass", "pairedJetbb_size_mass", 5, 0, 5, 5, hmmin, hmmax); listOfTH2.push_back(pairedJetbbsizemass);
}
*/



#endif