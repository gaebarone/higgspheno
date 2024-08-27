#ifndef HIST_INCLUDE_H
#define HIST_INCLUDE_H

#include <map>
#include <string>
#include <TH1F.h>
#include <TH2F.h>
#include <TMath.h>

#include "helperfunctions_include.h"

using namespace std;

class histograms {

public: 

  histograms(string name = " ") : sel_name( name ) {}
  string sel_name;

  std::map<std::string, TH1F*> bbhists;
  std::map<std::string, TH1F*> jjhists;
  std::map<std::string, TH1F*> wwhists;
  std::map<std::string, TH1F*> zzhists;
  std::map<std::string, TH1F*> lephists;
  std::map<std::string, TH1F*> methists;

  std::map<std::string, TH2F*> bbhists2D;
  std::map<std::string, TH2F*> jjhists2D;
  std::map<std::string, TH2F*> wwhists2D;
  std::map<std::string, TH2F*> zzhists2D;
  std::map<std::string, TH2F*> lephists2D;

  void initialize_bb() {

  // 1D

    bbhists[sel_name + "_bb_pT_reco"] = new TH1F( (sel_name + "_bb_pT_reco").c_str() , "p^{T}_{bb} Reco ; p^{T} (GeV) ; Events", 20, 0, 500);
    bbhists[sel_name + "_bb_m_reco"] = new TH1F( (sel_name + "_bb_m_reco").c_str() , "m_{bb} Reco ; Mass (GeV) ; Events", 20, 0, 500);
    bbhists[sel_name + "_bb_dphi_reco"] = new TH1F( (sel_name + "_bb_dphi_reco").c_str() , "#Delta#phi_{bb} Reco ; #Delta#phi ; Events", 10, -TMath::Pi(), +TMath::Pi());
    bbhists[sel_name + "_bb_deta_reco"] = new TH1F( (sel_name + "_bb_deta_reco").c_str() , "#Delta#eta_{bb} Reco ; #Delta#eta ; Events", 10, -10, 10);
    bbhists[sel_name + "_bb_dr_reco"] = new TH1F( (sel_name + "_bb_dr_reco").c_str() , "#DeltaR_{bb} Reco ; #DeltaR ; Events", 10, 0, 5);

    bbhists[sel_name + "_bb_pT_particle"] = new TH1F( (sel_name + "_bb_pT_particle").c_str() , "p^{T}_{bb} Particle ; p^{T} (GeV) ; Events", 20, 0, 500);
    bbhists[sel_name + "_bb_m_particle"] = new TH1F( (sel_name + "_bb_m_particle").c_str() , "m_{bb} Particle ; Mass (GeV) ; Events", 20, 0, 500);
    bbhists[sel_name + "_bb_dphi_particle"] = new TH1F( (sel_name + "_bb_dphi_particle").c_str() , "#Delta#phi_{bb} Particle ; #Delta#phi ; Events", 10, -TMath::Pi(), +TMath::Pi());
    bbhists[sel_name + "_bb_deta_particle"] = new TH1F( (sel_name + "_bb_deta_particle").c_str() , "#Delta#eta_{bb} Particle ; #Delta#eta ; Events", 10, -10, 10);
    bbhists[sel_name + "_bb_dr_particle"] = new TH1F( (sel_name + "_bb_dr_particle").c_str() , "#DeltaR_{bb} Particle ; #DeltaR ; Events", 10, 0, 5);

    bbhists[sel_name + "_bb_pT_parton"] = new TH1F( (sel_name + "_bb_pT_parton").c_str() , "p^{T}_{bb} Parton ; p^{T} (GeV) ; Events", 20, 0, 500);
    bbhists[sel_name + "_bb_m_parton"] = new TH1F( (sel_name + "_bb_m_parton").c_str() , "m_{bb} Parton ; Mass (GeV) ; Events", 20, 0, 500);
    bbhists[sel_name + "_bb_dphi_parton"] = new TH1F( (sel_name + "_bb_dphi_parton").c_str() , "#Delta#phi_{bb} Parton ; #Delta#phi ; Events", 10, -TMath::Pi(), +TMath::Pi());
    bbhists[sel_name + "_bb_deta_parton"] = new TH1F( (sel_name + "_bb_deta_parton").c_str() , "#Delta#eta_{bb} Parton ; #Delta#eta ; Events", 10, -10, 10);
    bbhists[sel_name + "_bb_dr_parton"] = new TH1F( (sel_name + "_bb_dr_parton").c_str() , "#DeltaR_{bb} Parton ; #DeltaR ; Events", 10, 0, 5);
  
  // 2D

    bbhists2D[sel_name + "_bb_pT_reco_particle"] = new TH2F( (sel_name + "_bb_pT_reco_particle").c_str() , "p^{T}_{bb} Reco Particle ; Reco p^{T} (GeV) ; Particle p^{T} (GeV)", 20, 0, 500, 20, 0, 500);
    bbhists2D[sel_name + "_bb_m_reco_particle"] = new TH2F( (sel_name + "_bb_m_reco_particle").c_str() , "m_{bb} Reco Particle ; Reco Mass (GeV) ; Particle Mass (GeV)", 20, 0, 500, 20, 0, 500);
    bbhists2D[sel_name + "_bb_dphi_reco_particle"] = new TH2F( (sel_name + "_bb_dphi_reco_particle").c_str() , "#Delta#phi_{bb} Reco Particle ; Reco #Delta#phi ; Particle #Delta#phi", 10, -TMath::Pi(), +TMath::Pi(), 10, -TMath::Pi(), +TMath::Pi());
    bbhists2D[sel_name + "_bb_deta_reco_particle"] = new TH2F( (sel_name + "_bb_deta_reco_particle").c_str() , "#Delta#eta_{bb} Reco Particle ; Reco #Delta#eta ; Particle #Delta#eta", 10, -10, 10, 10, -10, 10);
    bbhists2D[sel_name + "_bb_dr_reco_particle"] = new TH2F( (sel_name + "_bb_dr_reco_particle").c_str() , "#DeltaR_{bb} Reco Particle ; Reco #DeltaR ; Particle #DeltaR", 10, 0, 5, 10, 0, 5);

    bbhists2D[sel_name + "_bb_pT_particle_parton"] = new TH2F( (sel_name + "_bb_pT_particle_parton").c_str() , "p^{T}_{bb} Particle Parton ; Particle p^{T} (GeV) ; Parton p^{T} (GeV", 20, 0, 500, 20, 0, 500);
    bbhists2D[sel_name + "_bb_m_particle_parton"] = new TH2F( (sel_name + "_bb_m_particle_parton").c_str() , "m_{bb} Particle Parton ; Particle Mass (GeV) ; Parton Mass (GeV)", 20, 0, 500, 20, 0, 500);
    bbhists2D[sel_name + "_bb_dphi_particle_parton"] = new TH2F( (sel_name + "_bb_dphi_particle_parton").c_str() , "#Delta#phi_{bb} Particle Parton ; Particle #Delta#phi ; Parton #Delta#phi", 10, -TMath::Pi(), +TMath::Pi(), 10, -TMath::Pi(), +TMath::Pi());
    bbhists2D[sel_name + "_bb_deta_particle_parton"] = new TH2F( (sel_name + "_bb_deta_particle_parton").c_str() , "#Delta#eta_{bb} Particle Parton ; Particle #Delta#eta ; Parton #Delta#eta", 10, -10, 10, 10, -10, 10);
    bbhists2D[sel_name + "_bb_dr_particle_parton"] = new TH2F( (sel_name + "_bb_dr_particle_parton").c_str() , "#DeltaR_{bb} Particle Parton ; Particle #DeltaR ; Parton #DeltaR", 10, 0, 5, 10, 0, 5);

  }

  void initialize_jj() {

    jjhists[sel_name + "_jj_pT_reco"] = new TH1F( (sel_name + "_jj_pT_reco").c_str() , "p^{T}_{jj} Reco ; p^{T} (GeV) ; Events", 20, 0, 500);
    jjhists[sel_name + "_jj_m_reco"] = new TH1F( (sel_name + "_jj_m_reco").c_str() , "m_{jj} Reco ; Mass (GeV) ; Events", 40, 0, 2000);
    jjhists[sel_name + "_jj_dphi_reco"] = new TH1F( (sel_name + "_jj_dphi_reco").c_str() , "#Delta#phi_{jj} Reco ; #Delta#phi ; Events", 10, -TMath::Pi(), +TMath::Pi());
    jjhists[sel_name + "_jj_deta_reco"] = new TH1F( (sel_name + "_jj_deta_reco").c_str() , "|#Delta#eta_{jj}| Reco ; |#Delta#eta| ; Events", 20, 0, 10);
    jjhists[sel_name + "_jj_dr_reco"] = new TH1F( (sel_name + "_jj_dr_reco").c_str() , "#DeltaR_{jj} Reco ; #DeltaR ; Events", 10, 0, 5);

    jjhists[sel_name + "_jj_pT_particle"] = new TH1F( (sel_name + "_jj_pT_particle").c_str() , "p^{T}_{jj} Particle ; p^{T} (GeV) ; Events", 20, 0, 500);
    jjhists[sel_name + "_jj_m_particle"] = new TH1F( (sel_name + "_jj_m_particle").c_str() , "m_{jj} Particle ; Mass (GeV) ; Events", 40, 0, 2000);
    jjhists[sel_name + "_jj_dphi_particle"] = new TH1F( (sel_name + "_jj_dphi_particle").c_str() , "#Delta#phi_{jj} Particle ; #Delta#phi ; Events", 10, -TMath::Pi(), +TMath::Pi());
    jjhists[sel_name + "_jj_deta_particle"] = new TH1F( (sel_name + "_jj_deta_particle").c_str() , "|#Delta#eta_{jj}| Particle ; |#Delta#eta| ; Events", 20, 0, 10);
    jjhists[sel_name + "_jj_dr_particle"] = new TH1F( (sel_name + "_jj_dr_particle").c_str() , "#DeltaR_{jj} Particle ; #DeltaR ; Events", 10, 0, 5);

    jjhists[sel_name + "_jj_pT_parton"] = new TH1F( (sel_name + "_jj_pT_parton").c_str() , "p^{T}_{jj} Parton ; p^{T} (GeV) ; Events", 20, 0, 500);
    jjhists[sel_name + "_jj_m_parton"] = new TH1F( (sel_name + "_jj_m_parton").c_str() , "m_{jj} Parton ; Mass (GeV) ; Events", 40, 0, 2000);
    jjhists[sel_name + "_jj_dphi_parton"] = new TH1F( (sel_name + "_jj_dphi_parton").c_str() , "#Delta#phi_{jj} Parton ; #Delta#phi ; Events", 10, -TMath::Pi(), +TMath::Pi());
    jjhists[sel_name + "_jj_deta_parton"] = new TH1F( (sel_name + "_jj_deta_parton").c_str() , "|#Delta#eta_{jj}| Parton ; |#Delta#eta| ; Events", 20, 0, 10);
    jjhists[sel_name + "_jj_dr_parton"] = new TH1F( (sel_name + "_jj_dr_parton").c_str() , "#DeltaR_{jj} Parton ; #DeltaR ; Events", 10, 0, 5);

  // 2D

    jjhists2D[sel_name + "_jj_pT_reco_particle"] = new TH2F( (sel_name + "_jj_pT_reco_particle").c_str() , "p^{T}_{jj} Reco Particle ; Reco p^{T} (GeV) ; Particle p^{T} (GeV)", 20, 0, 500, 20, 0, 500);
    jjhists2D[sel_name + "_jj_m_reco_particle"] = new TH2F( (sel_name + "_jj_m_reco_particle").c_str() , "m_{jj} Reco Particle ; Reco Mass (GeV) ; Particle Mass (GeV)", 20, 0, 2000, 20, 0, 2000);
    jjhists2D[sel_name + "_jj_dphi_reco_particle"] = new TH2F( (sel_name + "_jj_dphi_reco_particle").c_str() , "#Delta#phi_{jj} Reco Particle ; Reco #Delta#phi ; Particle #Delta#phi", 10, -TMath::Pi(), +TMath::Pi(), 10, -TMath::Pi(), +TMath::Pi());
    jjhists2D[sel_name + "_jj_deta_reco_particle"] = new TH2F( (sel_name + "_jj_deta_reco_particle").c_str() , "|#Delta#eta_{jj}| Reco Particle ; Reco |#Delta#eta| ; Particle #Delta#eta", 20, 0, 10, 20, 0, 10);
    jjhists2D[sel_name + "_jj_dr_reco_particle"] = new TH2F( (sel_name + "_jj_dr_reco_particle").c_str() , "#DeltaR_{jj} Reco Particle ; Reco #DeltaR ; Particle #DeltaR", 10, 0, 5, 10, 0, 5);

    jjhists2D[sel_name + "_jj_pT_particle_parton"] = new TH2F( (sel_name + "_jj_pT_particle_parton").c_str() , "p^{T}_{jj} Particle Parton ; Particle p^{T} (GeV) ; Parton p^{T} (GeV)", 20, 0, 500, 20, 0, 500);
    jjhists2D[sel_name + "_jj_m_particle_parton"] = new TH2F( (sel_name + "_jj_m_particle_parton").c_str() , "m_{jj} Particle Parton ; Particle Mass (GeV) ; Parton Mass (GeV)", 20, 0, 2000, 20, 0, 2000);
    jjhists2D[sel_name + "_jj_dphi_particle_parton"] = new TH2F( (sel_name + "_jj_dphi_particle_parton").c_str() , "#Delta#phi_{jj} Particle Parton ; Particle #Delta#phi ; Parton #Delta#phi", 10, -TMath::Pi(), +TMath::Pi(), 10, -TMath::Pi(), +TMath::Pi());
    jjhists2D[sel_name + "_jj_deta_particle_parton"] = new TH2F( (sel_name + "_jj_deta_particle_parton").c_str() , "|#Delta#eta_{jj}| Particle Parton ; Particle |#Delta#eta| ; Parton #Delta#eta", 20, 0, 10, 20, 0, 10);
    jjhists2D[sel_name + "_jj_dr_particle_parton"] = new TH2F( (sel_name + "_jj_dr_particle_parton").c_str() , "#DeltaR_{jj} Particle Parton ; Particle #DeltaR ; Parton #DeltaR", 10, 0, 5, 10, 0, 5);

  }

  void initialize_ww() {

    wwhists[sel_name + "_w1_pT_reco"] = new TH1F( (sel_name + "_w1_pT_reco").c_str() , "p^{T}_{w1} Reco ; p^{T} (GeV) ; Events", 20, 0, 500);
    wwhists[sel_name + "_w1_m_reco"] = new TH1F( (sel_name + "_w1_m_reco").c_str() , "m_{w1} Reco ; Mass (GeV) ; Events", 20, 0, 500); wwhists[sel_name + "_w1_mT_reco"] = new TH1F( (sel_name + "_w1_mT_reco").c_str() , "m^{T}_{w1} Reco ; M^{T} (GeV) ; Events", 20, 0, 500);
    wwhists[sel_name + "_w2_pT_reco"] = new TH1F( (sel_name + "_w2_pT_reco").c_str() , "p^{T}_{w2} Reco ; p^{T} (GeV) ; Events", 20, 0, 500);
    wwhists[sel_name + "_w2_m_reco"] = new TH1F( (sel_name + "_w2_m_reco").c_str() , "m_{w2} Reco ; Mass (GeV) ; Events", 20, 0, 500); wwhists[sel_name + "_w2_mT_reco"] = new TH1F( (sel_name + "_w2_mT_reco").c_str() , "m^{T}_{w2} Reco ; M^{T} (GeV) ; Events", 20, 0, 500);
    wwhists[sel_name + "_ww_pT_reco"] = new TH1F( (sel_name + "_ww_pT_reco").c_str() , "p^{T}_{ww} Reco ; p^{T} (GeV) ; Events", 20, 0, 500);
    wwhists[sel_name + "_ww_m_reco"] = new TH1F( (sel_name + "_ww_m_reco").c_str() , "m_{ww} Reco ; Mass (GeV) ; Events", 20, 0, 500); wwhists[sel_name + "_ww_mT_reco"] = new TH1F( (sel_name + "_ww_mT_reco").c_str() , "m^{T}_{ww} Reco ; M^{T} (GeV) ; Events", 20, 0, 500);
    wwhists[sel_name + "_ww_dphi_reco"] = new TH1F( (sel_name + "_ww_dphi_reco").c_str() , "#Delta#phi_{ww} Reco ; #Delta#phi ; Events", 10, -TMath::Pi(), +TMath::Pi());
    wwhists[sel_name + "_ww_deta_reco"] = new TH1F( (sel_name + "_ww_deta_reco").c_str() , "#Delta#eta_{ww} Reco ; #Delta#eta ; Events", 10, -10, 10);
    wwhists[sel_name + "_ww_dr_reco"] = new TH1F( (sel_name + "_ww_dr_reco").c_str() , "#DeltaR_{ww} Reco ; #DeltaR ; Events", 10, 0, 5);
    wwhists[sel_name + "_ww_et_reco"] = new TH1F( (sel_name + "_ww_et_reco").c_str() , "Event Type Reco", 5, -1, 4);

    wwhists[sel_name + "_w1_pT_particle"] = new TH1F( (sel_name + "_w1_pT_particle").c_str() , "p^{T}_{w1} Particle ; p^{T} (GeV) ; Events", 20, 0, 500);
    wwhists[sel_name + "_w1_m_particle"] = new TH1F( (sel_name + "_w1_m_particle").c_str() , "m_{w1} Particle ; Mass (GeV) ; Events", 20, 0, 500); wwhists[sel_name + "_w1_mT_particle"] = new TH1F( (sel_name + "_w1_mT_particle").c_str() , "m^{T}_{w1} Particle ; M^{T} (GeV) ; Events", 20, 0, 500);
    wwhists[sel_name + "_w2_pT_particle"] = new TH1F( (sel_name + "_w2_pT_particle").c_str() , "p^{T}_{w2} Particle ; p^{T} (GeV) ; Events", 20, 0, 500);
    wwhists[sel_name + "_w2_m_particle"] = new TH1F( (sel_name + "_w2_m_particle").c_str() , "m_{w2} Particle ; Mass (GeV) ; Events", 20, 0, 500); wwhists[sel_name + "_w2_mT_particle"] = new TH1F( (sel_name + "_w2_mT_particle").c_str() , "m^{T}_{w2} Particle ; M^{T} (GeV) ; Events", 20, 0, 500);
    wwhists[sel_name + "_ww_pT_particle"] = new TH1F( (sel_name + "_ww_pT_particle").c_str() , "p^{T}_{ww} Particle ; p^{T} (GeV) ; Events", 20, 0, 500);
    wwhists[sel_name + "_ww_m_particle"] = new TH1F( (sel_name + "_ww_m_particle").c_str(), "m_{ww} Particle ; Mass (GeV) ; Events", 20, 0, 500); wwhists[sel_name + "_ww_mT_particle"] = new TH1F( (sel_name + "_ww_mT_particle").c_str(), "m^{T}_{ww} Particle ; M^{T} (GeV) ; Events", 20, 0, 500);
    wwhists[sel_name + "_ww_dphi_particle"] = new TH1F( (sel_name + "_ww_dphi_particle").c_str() , "#Delta#phi_{ww} Particle ; #Delta#phi ; Events", 10, -TMath::Pi(), +TMath::Pi());
    wwhists[sel_name + "_ww_deta_particle"] = new TH1F( (sel_name + "_ww_deta_particle").c_str() , "#Delta#eta_{ww} Particle ; #Delta#eta ; Events", 10, -10, 10);
    wwhists[sel_name + "_ww_dr_particle"] = new TH1F( (sel_name + "_ww_dr_particle").c_str() , "#DeltaR_{ww} Particle ; #DeltaR ; Events", 10, 0, 5);
    wwhists[sel_name + "_ww_et_particle"] = new TH1F( (sel_name + "_ww_et_particle").c_str() , "Event Type Particle", 5, -1, 4);

    wwhists[sel_name + "_w1_pT_parton"] = new TH1F( (sel_name + "_w1_pT_parton").c_str() , "p^{T}_{w1} Parton ; p^{T} (GeV) ; Events", 20, 0, 500);
    wwhists[sel_name + "_w1_m_parton"] = new TH1F( (sel_name + "_w1_m_parton").c_str() , "m_{w1} Parton ; Mass (GeV) ; Events", 20, 0, 500); wwhists[sel_name + "_w1_mT_parton"] = new TH1F( (sel_name + "_w1_mT_parton").c_str() , "m^{T}_{w1} Parton ; M^{T} (GeV) ; Events", 20, 0, 500);
    wwhists[sel_name + "_w2_pT_parton"] = new TH1F( (sel_name + "_w2_pT_parton").c_str() , "p^{T}_{w2} Parton ; p^{T} (GeV) ; Events", 20, 0, 500);
    wwhists[sel_name + "_w2_m_parton"] = new TH1F( (sel_name + "_w2_m_parton").c_str() , "m_{w2} Parton ; Mass (GeV) ; Events", 20, 0, 500); wwhists[sel_name + "_w2_mT_parton"] = new TH1F( (sel_name + "_w2_mT_parton").c_str() , "m^{T}_{w2} Parton ; M^{T} (GeV) ; Events", 20, 0, 500);
    wwhists[sel_name + "_ww_pT_parton"] = new TH1F( (sel_name + "_ww_pT_parton").c_str() , "p^{T}_{ww} Parton ; p^{T} (GeV) ; Events", 20, 0, 500);
    wwhists[sel_name + "_ww_m_parton"] = new TH1F( (sel_name + "_ww_m_parton").c_str(), "m_{ww} Parton ; Mass (GeV) ; Events", 20, 0, 500); wwhists[sel_name + "_ww_mT_parton"] = new TH1F( (sel_name + "_ww_mT_parton").c_str(), "m^{T}_{ww} Parton ; M^{T} (GeV) ; Events", 20, 0, 500);
    wwhists[sel_name + "_ww_dphi_parton"] = new TH1F( (sel_name + "_ww_dphi_parton").c_str() , "#Delta#phi_{ww} Parton ; #Delta#phi ; Events", 10, -TMath::Pi(), +TMath::Pi());
    wwhists[sel_name + "_ww_deta_parton"] = new TH1F( (sel_name + "_ww_deta_parton").c_str() , "#Delta#eta_{ww} Parton ; #Delta#eta ; Events", 10, -10, 10);
    wwhists[sel_name + "_ww_dr_parton"] = new TH1F( (sel_name + "_ww_dr_parton").c_str() , "#DeltaR_{ww} Parton ; #DeltaR ; Events", 10, 0, 5);
    wwhists[sel_name + "_ww_et_parton"] = new TH1F( (sel_name + "_ww_et_parton").c_str() , "Event Type Parton", 5, -1, 4);

  // 2D

    wwhists2D[sel_name + "_w1_pT_reco_particle"] = new TH2F( (sel_name + "_w1_pT_reco_particle").c_str() , "p^{T}_{w1} Reco Particle ; Reco p^{T} (GeV) ; Particle p^{T} (GeV)", 20, 0, 500, 20, 0, 500);
    wwhists2D[sel_name + "_w1_m_reco_particle"] = new TH2F( (sel_name + "_w1_m_reco_particle").c_str() , "m_{w1} Reco Particle ; Reco Mass (GeV) ; Particle Mass (GeV)", 20, 0, 500, 20, 0, 500); wwhists2D[sel_name + "_w1_mT_reco_particle"] = new TH2F( (sel_name + "_w1_mT_reco_particle").c_str() , "m^{T}_{w1} Reco Particle ; Reco M^{T} (GeV) ; Particle M^{T} (GeV)", 20, 0, 500, 20, 0, 500);
    wwhists2D[sel_name + "_w2_pT_reco_particle"] = new TH2F( (sel_name + "_w2_pT_reco_particle").c_str() , "p^{T}_{w2} Reco Particle ; Reco p^{T} (GeV) ; Particle p^{T} (GeV)", 20, 0, 500, 20, 0, 500);
    wwhists2D[sel_name + "_w2_m_reco_particle"] = new TH2F( (sel_name + "_w2_m_reco_particle").c_str() , "m_{w2} Reco Particle ; Reco Mass (GeV) ; Particle Mass (GeV)", 20, 0, 500, 20, 0, 500); wwhists2D[sel_name + "_w2_mT_reco_particle"] = new TH2F( (sel_name + "_w2_mT_reco_particle").c_str() , "m^{T}_{w2} Reco Particle ; Reco M^{T} (GeV) ; Particle M^{T} (GeV)", 20, 0, 500, 20, 0, 500);
    wwhists2D[sel_name + "_ww_pT_reco_particle"] = new TH2F( (sel_name + "_ww_pT_reco_particle").c_str() , "p^{T}_{ww} Reco Particle ; Reco p^{T} (GeV) ; Particle p^{T} (GeV)", 20, 0, 500, 20, 0, 500);
    wwhists2D[sel_name + "_ww_m_reco_particle"] = new TH2F( (sel_name + "_ww_m_reco_particle").c_str() , "m_{ww} Reco Particle ; Reco Mass (GeV) ; Particle Mass (GeV)", 20, 0, 500, 20, 0, 500); wwhists2D[sel_name + "_ww_mT_reco_particle"] = new TH2F( (sel_name + "_ww_mT_reco_particle").c_str() , "m^{T}_{ww} Reco Particle ; Reco M^{T} (GeV) ; Particle M^{T} (GeV)", 20, 0, 500, 20, 0, 500);
    wwhists2D[sel_name + "_ww_dphi_reco_particle"] = new TH2F( (sel_name + "_ww_dphi_reco_particle").c_str() , "#Delta#phi_{ww} Reco Particle ; Reco #Delta#phi ; Particle #Delta#phi", 10, -TMath::Pi(), +TMath::Pi(), 10, -TMath::Pi(), +TMath::Pi());
    wwhists2D[sel_name + "_ww_deta_reco_particle"] = new TH2F( (sel_name + "_ww_deta_reco_particle").c_str() , "#Delta#eta_{ww} Reco Particle ; Reco #Delta#eta ; Particle #Delta#eta", 10, -10, 10, 10, -10, 10);
    wwhists2D[sel_name + "_ww_dr_reco_particle"] = new TH2F( (sel_name + "_ww_dr_reco_particle").c_str() , "#DeltaR_{ww} Reco Particle ; Reco #DeltaR ; Particle #DeltaR", 10, 0, 5, 10, 0, 5);

    wwhists2D[sel_name + "_w1_pT_particle_parton"] = new TH2F( (sel_name + "_w1_pT_particle_parton").c_str() , "p^{T}_{w1} Particle Parton ; Particle p^{T} (GeV) ; Parton p^{T} (GeV)", 20, 0, 500, 20, 0, 500);
    wwhists2D[sel_name + "_w1_m_particle_parton"] = new TH2F( (sel_name + "_w1_m_particle_parton").c_str() , "m_{w1} Particle Parton ; Particle Mass (GeV) ; Parton Mass (GeV)", 20, 0, 500, 20, 0, 500); wwhists2D[sel_name + "_w1_mT_particle_parton"] = new TH2F( (sel_name + "_w1_mT_particle_parton").c_str() , "m^{T}_{w1} Particle Parton ; Particle M^{T} (GeV) ; Parton M^{T} (GeV)", 20, 0, 500, 20, 0, 500);
    wwhists2D[sel_name + "_w2_pT_particle_parton"] = new TH2F( (sel_name + "_w2_pT_particle_parton").c_str() , "p^{T}_{w2} Particle Parton ; Particle p^{T} (GeV) ; Parton p^{T} (GeV)", 20, 0, 500, 20, 0, 500);
    wwhists2D[sel_name + "_w2_m_particle_parton"] = new TH2F( (sel_name + "_w2_m_particle_parton").c_str() , "m_{w2} Particle Parton ; Particle Mass (GeV) ; Parton Mass (GeV)", 20, 0, 500, 20, 0, 500); wwhists2D[sel_name + "_w2_mT_particle_parton"] = new TH2F( (sel_name + "_w2_mT_particle_parton").c_str() , "m^{T}_{w2} Particle Parton ; Particle M^{T} (GeV) ; Parton M^{T} (GeV)", 20, 0, 500, 20, 0, 500);
    wwhists2D[sel_name + "_ww_pT_particle_parton"] = new TH2F( (sel_name + "_ww_pT_particle_parton").c_str() , "p^{T}_{ww} Particle Parton ; Particle p^{T} (GeV) ; Parton p^{T} (GeV)", 20, 0, 500, 20, 0, 500);
    wwhists2D[sel_name + "_ww_m_particle_parton"] = new TH2F( (sel_name + "_ww_m_particle_parton").c_str() , "m_{ww} Particle Parton ; Particle Mass (GeV) ; Parton Mass (GeV)", 20, 0, 500, 20, 0, 500); wwhists2D[sel_name + "_ww_mT_particle_parton"] = new TH2F( (sel_name + "_ww_mT_particle_parton").c_str() , "m^{T}_{ww} Particle Parton ; Particle M^{T} (GeV) ; Parton M^{T} (GeV)", 20, 0, 500, 20, 0, 500);
    wwhists2D[sel_name + "_ww_dphi_particle_parton"] = new TH2F( (sel_name + "_ww_dphi_particle_parton").c_str() , "#Delta#phi_{ww} Particle Parton ; Particle #Delta#phi ; Parton #Delta#phi", 10, -TMath::Pi(), +TMath::Pi(), 10, -TMath::Pi(), +TMath::Pi());
    wwhists2D[sel_name + "_ww_deta_particle_parton"] = new TH2F( (sel_name + "_ww_deta_particle_parton").c_str() , "#Delta#eta_{ww} Particle Parton ; Particle #Delta#eta ; Parton #Delta#eta", 10, -10, 10, 10, -10, 10);
    wwhists2D[sel_name + "_ww_dr_particle_parton"] = new TH2F( (sel_name + "_ww_dr_particle_parton").c_str() , "#DeltaR_{ww} Particle Parton ; Particle #DeltaR ; Parton #DeltaR", 10, 0, 5, 10, 0, 5);

    }

  void initialize_zz() {

    zzhists[sel_name + "_z1_pT_reco"] = new TH1F( (sel_name + "_z1_pT_reco").c_str() , "p^{T}_{z1} Reco ; p^{T} (GeV) ; Events", 20, 0, 500);
    zzhists[sel_name + "_z1_m_reco"] = new TH1F( (sel_name + "_z1_m_reco").c_str() , "m_{z1} Reco ; Mass (GeV) ; Events", 20, 0, 250);
    zzhists[sel_name + "_z2_pT_reco"] = new TH1F( (sel_name + "_z2_pT_reco").c_str() , "p^{T}_{z2} Reco ; p^{T} (GeV) ; Events", 20, 0, 500);
    zzhists[sel_name + "_z2_m_reco"] = new TH1F( (sel_name + "_z2_m_reco").c_str() , "m_{z2} Reco ; Mass (GeV) ; Events", 20, 0, 250);
    zzhists[sel_name + "_zz_pT_reco"] = new TH1F( (sel_name + "_zz_pT_reco").c_str() , "p^{T}_{zz} Reco ; p^{T} (GeV) ; Events", 20, 0, 500);
    zzhists[sel_name + "_zz_m_reco"] = new TH1F( (sel_name + "_zz_m_reco").c_str() , "m_{zz} Reco ; Mass (GeV) ; Events", 20, 0, 500);
    zzhists[sel_name + "_zz_dphi_reco"] = new TH1F( (sel_name + "_zz_dphi_reco").c_str() , "#Delta#phi_{zz} Reco ; #Delta#phi ; Events", 10, -TMath::Pi(), +TMath::Pi());
    zzhists[sel_name + "_zz_deta_reco"] = new TH1F( (sel_name + "_zz_deta_reco").c_str() , "#Delta#eta_{zz} Reco ; #Delta#eta ; Events", 10, -10, 10);
    zzhists[sel_name + "_zz_dr_reco"] = new TH1F( (sel_name + "_zz_dr_reco").c_str() , "#DeltaR_{zz} Reco ; #DeltaR ; Events", 10, 0, 5);
    zzhists[sel_name + "_zz_et_reco"] = new TH1F( (sel_name + "_zz_et_reco").c_str() , "Event Type Reco", 5, -1, 4);

    zzhists[sel_name + "_z1_pT_particle"] = new TH1F( (sel_name + "_z1_pT_particle").c_str() , "p^{T}_{z1} Particle ; p^{T} (GeV) ; Events", 20, 0, 500);
    zzhists[sel_name + "_z1_m_particle"] = new TH1F( (sel_name + "_z1_m_particle").c_str() , "m_{z1} Particle ; Mass (GeV) ; Events", 20, 0, 250);
    zzhists[sel_name + "_z2_pT_particle"] = new TH1F( (sel_name + "_z2_pT_particle").c_str() , "p^{T}_{z2} Particle ; p^{T} (GeV) ; Events", 20, 0, 500);
    zzhists[sel_name + "_z2_m_particle"] = new TH1F( (sel_name + "_z2_m_particle").c_str() , "m_{z2} Particle ; Mass (GeV) ; Events", 20, 0, 250);
    zzhists[sel_name + "_zz_pT_particle"] = new TH1F( (sel_name + "_zz_pT_particle").c_str() , "p^{T}_{zz} Particle ; p^{T} (GeV) ; Events", 20, 0, 250);
    zzhists[sel_name + "_zz_m_particle"] = new TH1F( (sel_name + "_zz_m_particle").c_str() , "m_{zz} Particle ; Mass (GeV) ; Events", 20, 0, 500);
    zzhists[sel_name + "_zz_dphi_particle"] = new TH1F( (sel_name + "_zz_dphi_particle").c_str() , "#Delta#phi_{zz} Particle ; #Delta#phi ; Events", 10, -TMath::Pi(), +TMath::Pi());
    zzhists[sel_name + "_zz_deta_particle"] = new TH1F( (sel_name + "_zz_deta_particle").c_str() , "#Delta#eta_{zz} Particle ; #Delta#eta ; Events", 10, -10, 10);
    zzhists[sel_name + "_zz_dr_particle"] = new TH1F( (sel_name + "_zz_dr_particle").c_str() , "#DeltaR_{zz} Particle ; #DeltaR ; Events", 10, 0, 5);
    zzhists[sel_name + "_zz_et_particle"] = new TH1F( (sel_name + "_zz_et_particle").c_str() , "Event Type Particle", 5, -1, 4);

    zzhists[sel_name + "_z1_pT_parton"] = new TH1F( (sel_name + "_z1_pT_parton").c_str() , "p^{T}_{z1} Parton ; p^{T} (GeV) ; Events", 20, 0, 500);
    zzhists[sel_name + "_z1_m_parton"] = new TH1F( (sel_name + "_z1_m_parton").c_str() , "m_{z1} Parton ; Mass (GeV) ; Events", 20, 0, 250);
    zzhists[sel_name + "_z2_pT_parton"] = new TH1F( (sel_name + "_z2_pT_parton").c_str() , "p^{T}_{z2} Parton ; p^{T} (GeV) ; Events", 20, 0, 500);
    zzhists[sel_name + "_z2_m_parton"] = new TH1F( (sel_name + "_z2_m_parton").c_str() , "m_{z2} Parton ; Mass (GeV) ; Events", 20, 0, 250);
    zzhists[sel_name + "_zz_pT_parton"] = new TH1F( (sel_name + "_zz_pT_parton").c_str() , "p^{T}_{zz} Parton ; p^{T} (GeV) ; Events", 20, 0, 500);
    zzhists[sel_name + "_zz_m_parton"] = new TH1F( (sel_name + "_zz_m_parton").c_str() , "m_{zz} Parton ; Mass (GeV) ; Events", 20, 0, 500);
    zzhists[sel_name + "_zz_dphi_parton"] = new TH1F( (sel_name + "_zz_dphi_parton").c_str() , "#Delta#phi_{zz} Parton ; #Delta#phi ; Events", 10, -TMath::Pi(), +TMath::Pi());
    zzhists[sel_name + "_zz_deta_parton"] = new TH1F( (sel_name + "_zz_deta_parton").c_str() , "#Delta#eta_{zz} Parton ; #Delta#eta ; Events", 10, -10, 10);
    zzhists[sel_name + "_zz_dr_parton"] = new TH1F( (sel_name + "_zz_dr_parton").c_str() , "#DeltaR_{zz} Parton ; #DeltaR ; Events", 10, 0, 5);
    zzhists[sel_name + "_zz_et_parton"] = new TH1F( (sel_name + "_zz_et_parton").c_str() , "Event Type Parton", 5, -1, 4);

  // 2D

    zzhists2D[sel_name + "_z1_pT_reco_particle"] = new TH2F( (sel_name + "_z1_pT_reco_particle").c_str() , "p^{T}_{z1} Reco Particle ; Reco p^{T} (GeV) ; Particle p^{T} (GeV)", 20, 0, 500, 20, 0, 500);
    zzhists2D[sel_name + "_z1_m_reco_particle"] = new TH2F( (sel_name + "_z1_m_reco_particle").c_str() , "m_{z1} Reco Particle ; Reco Mass (GeV) ; Particle Mass (GeV)", 20, 0, 250, 20, 0, 250);
    zzhists2D[sel_name + "_z2_pT_reco_particle"] = new TH2F( (sel_name + "_z2_pT_reco_particle").c_str() , "p^{T}_{z2} Reco Particle ; Reco p^{T} (GeV) ; Particle p^{T} (GeV)", 20, 0, 500, 20, 0, 500);
    zzhists2D[sel_name + "_z2_m_reco_particle"] = new TH2F( (sel_name + "_z2_m_reco_particle").c_str() , "m_{z2} Reco Particle ; Reco Mass (GeV) ; Particle Mass (GeV)", 20, 0, 250, 20, 0, 250);
    zzhists2D[sel_name + "_zz_pT_reco_particle"] = new TH2F( (sel_name + "_zz_pT_reco_particle").c_str() , "p^{T}_{zz} Reco Particle ; Reco p^{T} (GeV) ; Particle p^{T} (GeV)", 20, 0, 500, 20, 0, 500);
    zzhists2D[sel_name + "_zz_m_reco_particle"] = new TH2F( (sel_name + "_zz_m_reco_particle").c_str() , "m_{zz} Reco Particle ; Reco Mass (GeV) ; Particle Mass (GeV)", 20, 0, 500, 20, 0, 500);
    zzhists2D[sel_name + "_zz_dphi_reco_particle"] = new TH2F( (sel_name + "_zz_dphi_reco_particle").c_str() , "#Delta#phi_{zz} Reco Particle ; Reco #Delta#phi ; Particle #Delta#phi", 10, -TMath::Pi(), +TMath::Pi(), 10, -TMath::Pi(), +TMath::Pi());
    zzhists2D[sel_name + "_zz_deta_reco_particle"] = new TH2F( (sel_name + "_zz_deta_reco_particle").c_str() , "#Delta#eta_{zz} Reco Particle ; Reco #Delta#eta ; Particle #Delta#eta", 10, -10, 10, 10, -10, 10);
    zzhists2D[sel_name + "_zz_dr_reco_particle"] = new TH2F( (sel_name + "_zz_dr_reco_particle").c_str() , "#DeltaR_{zz} Reco Particle ; Reco #DeltaR ; Particle #DeltaR", 10, 0, 5, 10, 0, 5);

    zzhists2D[sel_name + "_z1_pT_particle_parton"] = new TH2F( (sel_name + "_z1_pT_particle_parton").c_str() , "p^{T}_{z1} Particle Parton ; Particle p^{T} (GeV) ; Parton p^{T} (GeV)", 20, 0, 500, 20, 0, 500);
    zzhists2D[sel_name + "_z1_m_particle_parton"] = new TH2F( (sel_name + "_z1_m_particle_parton").c_str() , "m_{z1} Particle Parton ; Particle Mass (GeV) ; Parton Mass (GeV)", 20, 0, 250, 20, 0, 250);
    zzhists2D[sel_name + "_z2_pT_particle_parton"] = new TH2F( (sel_name + "_z2_pT_particle_parton").c_str() , "p^{T}_{z2} Particle Parton ; Particle p^{T} (GeV) ; Parton p^{T} (GeV)", 20, 0, 500, 20, 0, 500);
    zzhists2D[sel_name + "_z2_m_particle_parton"] = new TH2F( (sel_name + "_z2_m_particle_parton").c_str() , "m_{z2} Particle Parton ; Particle Mass (GeV) ; Parton Mass (GeV)", 20, 0, 250, 20, 0, 250);
    zzhists2D[sel_name + "_zz_pT_particle_parton"] = new TH2F( (sel_name + "_zz_pT_particle_parton").c_str() , "p^{T}_{zz} Particle Parton ; Particle p^{T} (GeV) ; Parton p^{T} (GeV)", 20, 0, 500, 20, 0, 500);
    zzhists2D[sel_name + "_zz_m_particle_parton"] = new TH2F( (sel_name + "_zz_m_particle_parton").c_str() , "m_{zz} Particle Parton ; Particle Mass (GeV) ; Parton Mass (GeV)", 20, 0, 500, 20, 0, 500);
    zzhists2D[sel_name + "_zz_dphi_particle_parton"] = new TH2F( (sel_name + "_zz_dphi_particle_parton").c_str() , "#Delta#phi_{zz} Particle Parton ; Particle #Delta#phi ; Parton #Delta#phi", 10, -TMath::Pi(), +TMath::Pi(), 10, -TMath::Pi(), +TMath::Pi());
    zzhists2D[sel_name + "_zz_deta_particle_parton"] = new TH2F( (sel_name + "_zz_deta_particle_parton").c_str() , "#Delta#eta_{zz} Particle Parton ; Particle #Delta#eta ; Parton #Delta#eta", 10, -10, 10, 10, -10, 10);
    zzhists2D[sel_name + "_zz_dr_particle_parton"] = new TH2F( (sel_name + "_zz_dr_particle_parton").c_str() , "#DeltaR_{zz} Particle Parton ; Particle #DeltaR ; Parton #DeltaR", 10, 0, 5, 10, 0, 5);

  }

  void initialize_2l() {

    lephists[sel_name + "_l1_pT_reco"] = new TH1F( (sel_name + "_l1_pT_reco").c_str() , "p^{T}_{l1} Reco ; p^{T} (GeV) ; Events", 50, 0, 250);
    lephists[sel_name + "_l2_pT_reco"] = new TH1F( (sel_name + "_l2_pT_reco").c_str() , "p^{T}_{l2} Reco ; p^{T} (GeV) ; Events", 50, 0, 250);
    lephists[sel_name + "_l1_eta_reco"] = new TH1F( (sel_name + "_l1_eta_reco").c_str() , "#eta_{l1} Reco ; #eta ; Events", 40, -10, 10);
    lephists[sel_name + "_l2_eta_reco"] = new TH1F( (sel_name + "_l2_eta_reco").c_str() , "#eta_{l2} Reco ; #eta ; Events", 40, -10, 10);
    lephists[sel_name + "_l1l2_pT_reco"] = new TH1F( (sel_name + "_l1l2_pT_reco").c_str() , "p^{T}_{l1l2} Reco ; p^{T} (GeV) ; Events", 10, 0, 250);
    lephists[sel_name + "_l1l2_m_reco"] = new TH1F( (sel_name + "_l1l2_m_reco").c_str() , "m_{l1l2} Reco ; Mass (GeV) ; Events", 20, 0, 250);
    lephists[sel_name + "_l1l2_dphi_reco"] = new TH1F( (sel_name + "_l1l2_dphi_reco").c_str() , "#Delta#phi_{l1l2} Reco ; #Delta#phi ; Events", 10, -TMath::Pi(), +TMath::Pi());
    lephists[sel_name + "_l1l2_deta_reco"] = new TH1F( (sel_name + "_l1l2_deta_reco").c_str() , "#Delta#eta_{l1l2} Reco ; #Delta#eta ; Events", 10, -10, 10);
    lephists[sel_name + "_l1l2_dr_reco"] = new TH1F( (sel_name + "_l1l2_dr_reco").c_str() , "#DeltaR_{l1l2} Reco ; #DeltaR ; Events", 10, 0, 5);

    lephists[sel_name + "_l1_pT_particle"] = new TH1F( (sel_name + "_l1_pT_particle").c_str() , "p^{T}_{l1} Particle ; p^{T} (GeV) ; Events", 50, 0, 250);
    lephists[sel_name + "_l2_pT_particle"] = new TH1F( (sel_name + "_l2_pT_particle").c_str() , "p^{T}_{l2} Particle ; p^{T} (GeV) ; Events", 50, 0, 250);
    lephists[sel_name + "_l1_eta_particle"] = new TH1F( (sel_name + "_l1_eta_particle").c_str() , "#eta_{l1} Particle ; #eta ; Events", 40, -10, 10);
    lephists[sel_name + "_l2_eta_particle"] = new TH1F( (sel_name + "_l2_eta_particle").c_str() , "#eta_{l2} Particle ; #eta ; Events", 40, -10, 10);
    lephists[sel_name + "_l1l2_pT_particle"] = new TH1F( (sel_name + "_l1l2_pT_particle").c_str() , "p^{T}_{l1l2} Particle ; p^{T} (GeV) ; Events", 10, 0, 250);
    lephists[sel_name + "_l1l2_m_particle"] = new TH1F( (sel_name + "_l1l2_m_particle").c_str() , "m_{l1l2} Particle ; Mass (GeV) ; Events", 20, 0, 250);
    lephists[sel_name + "_l1l2_dphi_particle"] = new TH1F( (sel_name + "_l1l2_dphi_particle").c_str() , "#Delta#phi_{l1l2} Particle ; #Delta#phi ; Events", 10, -TMath::Pi(), +TMath::Pi());
    lephists[sel_name + "_l1l2_deta_particle"] = new TH1F( (sel_name + "_l1l2_deta_particle").c_str() , "#Delta#eta_{l1l2} Particle ; #Delta#eta ; Events", 10, -10, 10);
    lephists[sel_name + "_l1l2_dr_particle"] = new TH1F( (sel_name + "_l1l2_dr_particle").c_str() , "#DeltaR_{l1l2} Particle ; #DeltaR ; Events", 10, 0, 5);

    lephists[sel_name + "_l1_pT_parton"] = new TH1F( (sel_name + "_l1_pT_parton").c_str() , "p^{T}_{l1} Parton ; p^{T} (GeV) ; Events", 50, 0, 250);
    lephists[sel_name + "_l2_pT_parton"] = new TH1F( (sel_name + "_l2_pT_parton").c_str() , "p^{T}_{l2} Parton ; p^{T} (GeV) ; Events", 50, 0, 250);
    lephists[sel_name + "_l1_eta_parton"] = new TH1F( (sel_name + "_l1_eta_parton").c_str() , "#eta_{l1} Parton ; #eta ; Events", 40, -10, 10);
    lephists[sel_name + "_l2_eta_parton"] = new TH1F( (sel_name + "_l2_eta_parton").c_str() , "#eta_{l2} Parton ; #eta ; Events", 40, -10, 10);
    lephists[sel_name + "_l1l2_pT_parton"] = new TH1F( (sel_name + "_l1l2_pT_parton").c_str() , "p^{T}_{l1l2} Parton ; p^{T} (GeV) ; Events", 10, 0, 250);
    lephists[sel_name + "_l1l2_m_parton"] = new TH1F( (sel_name + "_l1l2_m_parton").c_str() , "m_{l1l2} Parton ; Mass (GeV) ; Events", 20, 0, 250);
    lephists[sel_name + "_l1l2_dphi_parton"] = new TH1F( (sel_name + "_l1l2_dphi_parton").c_str() , "#Delta#phi_{l1l2} Parton ; #Delta#phi ; Events", 10, -TMath::Pi(), +TMath::Pi());
    lephists[sel_name + "_l1l2_deta_parton"] = new TH1F( (sel_name + "_l1l2_deta_parton").c_str() , "#Delta#eta_{l1l2} Parton ; #Delta#eta ; Events", 10, -10, 10);
    lephists[sel_name + "_l1l2_dr_parton"] = new TH1F( (sel_name + "_l1l2_dr_parton").c_str() , "#DeltaR_{l1l2} Parton ; #DeltaR ; Events", 10, 0, 5);

  // 2D

    lephists2D[sel_name + "_l1_pT_reco_particle"] = new TH2F( (sel_name + "_l1_pT_reco_particle").c_str() , "p^{T}_{l1} Reco Particle ; Reco p^{T} (GeV) ; Particle p^{T} (GeV)", 50, 0, 250, 50, 0, 250);
    lephists2D[sel_name + "_l2_pT_reco_particle"] = new TH2F( (sel_name + "_l2_pT_reco_particle").c_str() , "p^{T}_{l2} Reco Particle ; Reco p^{T} (GeV) ; Particle p^{T} (GeV)", 50, 0, 250, 50, 0, 250);
    lephists2D[sel_name + "_l1_eta_reco_particle"] = new TH2F( (sel_name + "_l1_eta_reco_particle").c_str() , "#eta_{l1} Reco Particle ; Reco #eta ; Particle #eta", 40, -10, 10, 40, -10, 10);
    lephists2D[sel_name + "_l2_eta_reco_particle"] = new TH2F( (sel_name + "_l2_eta_reco_particle").c_str() , "#eta_{l2} Reco Particle ; Reco #eta ; Particle #eta", 40, -10, 10, 40, -10, 10);
    lephists2D[sel_name + "_l1l2_pT_reco_particle"] = new TH2F( (sel_name + "_l1l2_pT_reco_particle").c_str() , "p^{T}_{l1l2} Reco Particle ; Reco p^{T} (GeV) ; Particle p^{T} (GeV)", 10, 0, 250, 10, 0, 250);
    lephists2D[sel_name + "_l1l2_m_reco_particle"] = new TH2F( (sel_name + "_l1l2_m_reco_particle").c_str() , "m_{l1l2} Reco Particle ; Reco Mass (GeV) ; Particle Mass (GeV)", 20, 0, 250, 20, 0, 250);
    lephists2D[sel_name + "_l1l2_dphi_reco_particle"] = new TH2F( (sel_name + "_l1l2_dphi_reco_particle").c_str() , "#Delta#phi_{l1l2} Reco Particle ; Reco #Delta#phi ; Particle #Delta#phi", 10, -TMath::Pi(), +TMath::Pi(), 10, -TMath::Pi(), +TMath::Pi());
    lephists2D[sel_name + "_l1l2_deta_reco_particle"] = new TH2F( (sel_name + "_l1l2_deta_reco_particle").c_str() , "#Delta#eta_{l1l2} Reco Particle ; Reco #Delta#eta ; Particle #Delta#eta", 10, -10, 10, 10, -10, 10);
    lephists2D[sel_name + "_l1l2_dr_reco_particle"] = new TH2F( (sel_name + "_l1l2_dr_reco_particle").c_str() , "#DeltaR_{l1l2} Reco Particle ; Reco #DeltaR ; Particle #DeltaR", 10, 0, 5, 10, 0, 5);

    lephists2D[sel_name + "_l1_pT_particle_parton"] = new TH2F( (sel_name + "_l1_pT_particle_parton").c_str() , "p^{T}_{l1} Particle Parton ; Particle p^{T} (GeV) ; Parton p^{T} (GeV)", 50, 0, 250, 50, 0, 250);
    lephists2D[sel_name + "_l2_pT_particle_parton"] = new TH2F( (sel_name + "_l2_pT_particle_parton").c_str() , "p^{T}_{l2} Particle Parton ; Particle p^{T} (GeV) ; Parton p^{T} (GeV)", 50, 0, 250, 50, 0, 250);
    lephists2D[sel_name + "_l1_eta_particle_parton"] = new TH2F( (sel_name + "_l1_eta_particle_parton").c_str() , "#eta_{l1} Particle Parton ; Particle #eta ; Parton #eta", 40, -10, 10, 40, -10, 10);
    lephists2D[sel_name + "_l2_eta_particle_parton"] = new TH2F( (sel_name + "_l2_eta_particle_parton").c_str() , "#eta_{l2} Particle Parton ; Particle #eta ; Parton #eta", 40, -10, 10, 40, -10, 10);
    lephists2D[sel_name + "_l1l2_pT_particle_parton"] = new TH2F( (sel_name + "_l1l2_pT_particle_parton").c_str() , "p^{T}_{l1l2} Particle Parton ; Particle p^{T} (GeV) ; Parton p^{T} (GeV)", 10, 0, 250, 10, 0, 250);
    lephists2D[sel_name + "_l1l2_m_particle_parton"] = new TH2F( (sel_name + "_l1l2_m_particle_parton").c_str() , "m_{l1l2} Particle Parton ; Particle Mass (GeV) ; Parton Mass (GeV)", 20, 0, 250, 20, 0, 250);
    lephists2D[sel_name + "_l1l2_dphi_particle_parton"] = new TH2F( (sel_name + "_l1l2_dphi_particle_parton").c_str() , "#Delta#phi_{l1l2} Particle Parton ; Particle #Delta#phi ; Parton #Delta#phi", 10, -TMath::Pi(), +TMath::Pi(), 10, -TMath::Pi(), +TMath::Pi());
    lephists2D[sel_name + "_l1l2_deta_particle_parton"] = new TH2F( (sel_name + "_l1l2_deta_particle_parton").c_str() , "#Delta#eta_{l1l2} Particle Parton ; Particle #Delta#eta ; Parton #Delta#eta", 10, -10, 10, 10, -10, 10);
    lephists2D[sel_name + "_l1l2_dr_particle_parton"] = new TH2F( (sel_name + "_l1l2_dr_particle_parton").c_str() , "#DeltaR_{l1l2} Particle Parton ; Particle #DeltaR ; Parton #DeltaR", 10, 0, 5, 10, 0, 5);

  } 

  void initialize_4l() {

    lephists[sel_name + "_l3_pT_reco"] = new TH1F( (sel_name + "_l3_pT_reco").c_str() , "p^{T}_{l3} Reco ; p^{T} (GeV) ; Events", 50, 0, 250);
    lephists[sel_name + "_l4_pT_reco"] = new TH1F( (sel_name + "_l4_pT_reco").c_str() , "p^{T}_{l4} Reco ; p^{T} (GeV) ; Events", 50, 0, 250);
    lephists[sel_name + "_l3_eta_reco"] = new TH1F( (sel_name + "_l3_eta_reco").c_str() , "#eta_{l3} Reco ; #eta ; Events", 40, -10, 10);
    lephists[sel_name + "_l4_eta_reco"] = new TH1F( (sel_name + "_l4_eta_reco").c_str() , "#eta_{l4} Reco ; #eta ; Events", 40, -10, 10);
    lephists[sel_name + "_l3l4_pT_reco"] = new TH1F( (sel_name + "_l3l4_pT_reco").c_str() , "p^{T}_{l3l4} Reco ; p^{T} (GeV) ; Events", 10, 0, 250);
    lephists[sel_name + "_l3l4_m_reco"] = new TH1F( (sel_name + "_l3l4_m_reco").c_str() , "m_{l3l4} Reco ; Mass (GeV) ; Events", 20, 0, 250);
    lephists[sel_name + "_l3l4_dphi_reco"] = new TH1F( (sel_name + "_l3l4_dphi_reco").c_str() , "#Delta#phi_{l3l4} Reco ; #Delta#phi ; Events", 10, -TMath::Pi(), +TMath::Pi());
    lephists[sel_name + "_l3l4_deta_reco"] = new TH1F( (sel_name + "_l3l4_deta_reco").c_str() , "#Delta#eta_{l3l4} Reco ; #Delta#eta ; Events", 10, -10, 10);
    lephists[sel_name + "_l3l4_dr_reco"] = new TH1F( (sel_name + "_l3l4_dr_reco").c_str() , "#DeltaR_{l3l4} Reco ; #DeltaR ; Events", 10, 0, 5);

    lephists[sel_name + "_l3_pT_particle"] = new TH1F( (sel_name + "_l3_pT_particle").c_str() , "p^{T}_{l3} Particle ; p^{T} (GeV) ; Events", 50, 0, 250);
    lephists[sel_name + "_l4_pT_particle"] = new TH1F( (sel_name + "_l4_pT_particle").c_str() , "p^{T}_{l4} Particle ; p^{T} (GeV) ; Events", 50, 0, 250);
    lephists[sel_name + "_l3_eta_particle"] = new TH1F( (sel_name + "_l3_eta_particle").c_str() , "#eta_{l3} Particle ; #eta ; Events", 40, -10, 10);
    lephists[sel_name + "_l4_eta_particle"] = new TH1F( (sel_name + "_l4_eta_particle").c_str() , "#eta_{l4} Particle ; #eta ; Events", 40, -10, 10);
    lephists[sel_name + "_l3l4_pT_particle"] = new TH1F( (sel_name + "_l3l4_pT_particle").c_str() , "p^{T}_{l3l4} Particle ; p^{T} (GeV) ; Events", 10, 0, 250);
    lephists[sel_name + "_l3l4_m_particle"] = new TH1F( (sel_name + "_l3l4_m_particle").c_str() , "m_{l3l4} Particle ; Mass (GeV) ; Events", 20, 0, 250);
    lephists[sel_name + "_l3l4_dphi_particle"] = new TH1F( (sel_name + "_l3l4_dphi_particle").c_str() , "#Delta#phi_{l3l4} Particle ; #Delta#phi ; Events", 10, -TMath::Pi(), +TMath::Pi());
    lephists[sel_name + "_l3l4_deta_particle"] = new TH1F( (sel_name + "_l3l4_deta_particle").c_str() , "#Delta#eta_{l3l4} Particle ; #Delta#eta ; Events", 10, -10, 10);
    lephists[sel_name + "_l3l4_dr_particle"] = new TH1F( (sel_name + "_l3l4_dr_particle").c_str() , "#DeltaR_{l3l4} Particle ; #DeltaR ; Events", 10, 0, 5);

    lephists[sel_name + "_l3_pT_parton"] = new TH1F( (sel_name + "_l3_pT_parton").c_str() , "p^{T}_{l3} Parton ; p^{T} (GeV) ; Events", 50, 0, 250);
    lephists[sel_name + "_l4_pT_parton"] = new TH1F( (sel_name + "_l4_pT_parton").c_str() , "p^{T}_{l4} Parton ; p^{T} (GeV) ; Events", 50, 0, 250);
    lephists[sel_name + "_l3_eta_parton"] = new TH1F( (sel_name + "_l3_eta_parton").c_str() , "#eta_{l3} Parton ; #eta ; Events", 40, -10, 10);
    lephists[sel_name + "_l4_eta_parton"] = new TH1F( (sel_name + "_l4_eta_parton").c_str() , "#eta_{l4} Parton ; #eta ; Events", 40, -10, 10);
    lephists[sel_name + "_l3l4_pT_parton"] = new TH1F( (sel_name + "_l3l4_pT_parton").c_str() , "p^{T}_{l3l4} Parton ; p^{T} (GeV) ; Events", 10, 0, 250);
    lephists[sel_name + "_l3l4_m_parton"] = new TH1F( (sel_name + "_l3l4_m_parton").c_str() , "m_{l3l4} Parton ; Mass (GeV) ; Events", 20, 0, 250);
    lephists[sel_name + "_l3l4_dphi_parton"] = new TH1F( (sel_name + "_l3l4_dphi_parton").c_str() , "#Delta#phi_{l3l4} Parton ; #Delta#phi ; Events", 10, -TMath::Pi(), +TMath::Pi());
    lephists[sel_name + "_l3l4_deta_parton"] = new TH1F( (sel_name + "_l3l4_deta_parton").c_str() , "#Delta#eta_{l3l4} Parton ; #Delta#eta ; Events", 10, -10, 10);
    lephists[sel_name + "_l3l4_dr_parton"] = new TH1F( (sel_name + "_l3l4_dr_parton").c_str() , "#DeltaR_{l3l4} Parton ; #DeltaR ; Events", 10, 0, 5);

  // 2D

    lephists2D[sel_name + "_l3_pT_reco_particle"] = new TH2F( (sel_name + "_l3_pT_reco_particle").c_str() , "p^{T}_{l3} Reco Particle ; Reco p^{T} (GeV) ; Particle p^{T} (GeV)", 50, 0, 250, 50, 0, 250);
    lephists2D[sel_name + "_l4_pT_reco_particle"] = new TH2F( (sel_name + "_l4_pT_reco_particle").c_str() , "p^{T}_{l4} Reco Particle ; Reco p^{T} (GeV) ; Particle p^{T} (GeV)", 50, 0, 250, 50, 0, 250);
    lephists2D[sel_name + "_l3_eta_reco_particle"] = new TH2F( (sel_name + "_l3_eta_reco_particle").c_str() , "#eta_{l3} Reco Particle ; Reco #eta ; Particle #eta", 40, -10, 10, 40, -10, 10);
    lephists2D[sel_name + "_l4_eta_reco_particle"] = new TH2F( (sel_name + "_l4_eta_reco_particle").c_str() , "#eta_{l4} Reco Particle ; Reco #eta ; Particle #eta", 40, -10, 10, 40, -10, 10);
    lephists2D[sel_name + "_l3l4_pT_reco_particle"] = new TH2F( (sel_name + "_l3l4_pT_reco_particle").c_str() , "p^{T}_{l3l4} Reco Particle ; Reco p^{T} (GeV) ; Particle p^{T} (GeV)", 10, 0, 250, 10, 0, 250);
    lephists2D[sel_name + "_l3l4_m_reco_particle"] = new TH2F( (sel_name + "_l3l4_m_reco_particle").c_str() , "m_{l3l4} Reco Particle ; Reco Mass (GeV) ; Particle Mass (GeV)", 20, 0, 250, 20, 0, 250);
    lephists2D[sel_name + "_l3l4_dphi_reco_particle"] = new TH2F( (sel_name + "_l3l4_dphi_reco_particle").c_str() , "#Delta#phi_{l3l4} Reco Particle ; Reco #Delta#phi ; Particle #Delta#phi", 10, -TMath::Pi(), +TMath::Pi(), 10, -TMath::Pi(), +TMath::Pi());
    lephists2D[sel_name + "_l3l4_deta_reco_particle"] = new TH2F( (sel_name + "_l3l4_deta_reco_particle").c_str() , "#Delta#eta_{l3l4} Reco Particle ; Reco #Delta#eta ; Particle #Delta#eta", 10, -10, 10, 10, -10, 10);
    lephists2D[sel_name + "_l3l4_dr_reco_particle"] = new TH2F( (sel_name + "_l3l4_dr_reco_particle").c_str() , "#DeltaR_{l3l4} Reco Particle ; Reco #DeltaR ; Particle #DeltaR", 10, 0, 5, 10, 0, 5);

    lephists2D[sel_name + "_l3_pT_particle_parton"] = new TH2F( (sel_name + "_l3_pT_particle_parton").c_str() , "p^{T}_{l3} Particle Parton ; Particle p^{T} (GeV) ; Parton p^{T} (GeV)", 50, 0, 250, 50, 0, 250);
    lephists2D[sel_name + "_l4_pT_particle_parton"] = new TH2F( (sel_name + "_l4_pT_particle_parton").c_str() , "p^{T}_{l4} Particle Parton ; Particle p^{T} (GeV) ; Parton p^{T} (GeV)", 50, 0, 250, 50, 0, 250);
    lephists2D[sel_name + "_l3_eta_particle_parton"] = new TH2F( (sel_name + "_l3_eta_particle_parton").c_str() , "#eta_{l3} Particle Parton ; Particle #eta ; Parton #eta", 40, -10, 10, 40, -10, 10);
    lephists2D[sel_name + "_l4_eta_particle_parton"] = new TH2F( (sel_name + "_l4_eta_particle_parton").c_str() , "#eta_{l4} Particle Parton ; Particle #eta ; Parton #eta", 40, -10, 10, 40, -10, 10);
    lephists2D[sel_name + "_l3l4_pT_particle_parton"] = new TH2F( (sel_name + "_l3l4_pT_particle_parton").c_str() , "p^{T}_{l3l4} Particle Parton ; Particle p^{T} (GeV) ; Parton p^{T} (GeV)", 10, 0, 250, 10, 0, 250);
    lephists2D[sel_name + "_l3l4_m_particle_parton"] = new TH2F( (sel_name + "_l3l4_m_particle_parton").c_str() , "m_{l3l4} Particle Parton ; Particle Mass (GeV) ; Parton Mass (GeV)", 20, 0, 250, 20, 0, 250);
    lephists2D[sel_name + "_l3l4_dphi_particle_parton"] = new TH2F( (sel_name + "_l3l4_dphi_particle_parton").c_str() , "#Delta#phi_{l3l4} Particle Parton ; Particle #Delta#phi ; Parton #Delta#phi", 10, -TMath::Pi(), +TMath::Pi(), 10, -TMath::Pi(), +TMath::Pi());
    lephists2D[sel_name + "_l3l4_deta_particle_parton"] = new TH2F( (sel_name + "_l3l4_deta_particle_parton").c_str() , "#Delta#eta_{l3l4} Particle Parton ; Particle #Delta#eta ; Parton #Delta#eta", 10, -10, 10, 10, -10, 10);
    lephists2D[sel_name + "_l3l4_dr_particle_parton"] = new TH2F( (sel_name + "_l3l4_dr_particle_parton").c_str() , "#DeltaR_{l3l4} Particle Parton ; Particle #DeltaR ; Parton #DeltaR", 10, 0, 5, 10, 0, 5);

  } 

  void initialize_met() {

    methists[sel_name + "_MET_reco"] = new TH1F( (sel_name + "_MET_reco").c_str() , "E^{T}_{Missing} Reco ; E^{T}_{Missing} (GeV) ; Events", 50, 0, 500);
    methists[sel_name + "_MET_particle"] = new TH1F( (sel_name + "_MET_particle").c_str() , "E^{T}_{Missing} Particle ; E^{T}_{Missing} (GeV) ; Events", 50, 0, 500);

  }

  void fill_bb( string analysis_type, double weight, TLorentzVector &b1, TLorentzVector &b2, string process_name ) {

    double pT = -99; double m = -99; double dphi = -99; double deta = -99; double dr = -99; 

    if ( !(b1.Pt() == 0 && b1.Eta() == 0 && b1.Phi() == 0) && !(b2.Pt() == 0 && b2.Eta() == 0 && b2.Phi() == 0) ) {

      pT = ( b1 + b2 ).Pt();
      m = ( b1 + b2 ).M();

      dphi = delta_phi( b1 , b2 );
      deta = delta_eta( b1 , b2 );
      dr = delta_r( b1 , b2 );
    
    }

    if ( analysis_type == "reco" ){

      bbhists[ sel_name + "_bb_pT_reco" ]->Fill(pT, weight);
      bbhists[ sel_name + "_bb_m_reco" ]->Fill(m, weight);
      bbhists[ sel_name + "_bb_dphi_reco" ]->Fill(dphi, weight);
      bbhists[ sel_name + "_bb_deta_reco" ]->Fill(deta, weight);
      bbhists[ sel_name + "_bb_dr_reco" ]->Fill(dr, weight);

    } else if ( analysis_type == "particle" ){

      bbhists[ sel_name + "_bb_pT_particle" ]->Fill(pT, weight);
      bbhists[ sel_name + "_bb_m_particle" ]->Fill(m, weight);
      bbhists[ sel_name + "_bb_dphi_particle" ]->Fill(dphi, weight);
      bbhists[ sel_name + "_bb_deta_particle" ]->Fill(deta, weight);
      bbhists[ sel_name + "_bb_dr_particle" ]->Fill(dr, weight);

    } else if ( analysis_type == "parton" ){

      bbhists[ sel_name + "_bb_pT_parton" ]->Fill(pT, weight);
      bbhists[ sel_name + "_bb_m_parton" ]->Fill(m, weight);
      bbhists[ sel_name + "_bb_dphi_parton" ]->Fill(dphi, weight);
      bbhists[ sel_name + "_bb_deta_parton" ]->Fill(deta, weight);
      bbhists[ sel_name + "_bb_dr_parton" ]->Fill(dr, weight);

    }

  }

  void fill_bb_2D( string analysis_type1, string analysis_type2, double weight, TLorentzVector &at1_b1, TLorentzVector &at1_b2, TLorentzVector &at2_b1, TLorentzVector &at2_b2 ) {

    double pT1 = -99; double m1 = -99; double dphi1 = -99; double deta1 = -99; double dr1 = -99;
    double pT2 = -99; double m2 = -99; double dphi2 = -99; double deta2 = -99; double dr2 = -99;

    if ( !(at1_b1.Pt() == 0 && at1_b1.Eta() == 0 && at1_b1.Phi() == 0) && !(at1_b2.Pt() == 0 && at1_b2.Eta() == 0 && at1_b2.Phi() == 0) && !(at2_b1.Pt() == 0 && at2_b1.Eta() == 0 && at2_b1.Phi() == 0) && !(at2_b2.Pt() == 0 && at2_b2.Eta() == 0 && at2_b2.Phi() == 0) ) {

    pT1 = ( at1_b1 + at1_b2 ).Pt(); pT2 = ( at2_b1 + at2_b2 ).Pt();
    m1 = ( at1_b1 + at1_b2 ).M(); m2 = ( at2_b1 + at2_b2 ).M();

    dphi1 = delta_phi( at1_b1 , at1_b2 ); dphi2 = delta_phi( at2_b1 , at2_b2 );
    deta1 = delta_eta( at1_b1 , at1_b2 ); deta2 = delta_eta( at2_b1 , at2_b2 );
    dr1 = delta_r( at1_b1 , at1_b2 ); dr2 = delta_r( at2_b1 , at2_b2 );

    }

    if ( analysis_type1 == "reco" && analysis_type2 == "particle" ){

      bbhists2D[ sel_name + "_bb_pT_reco_particle" ]->Fill(pT1, pT2, weight);
      bbhists2D[ sel_name + "_bb_m_reco_particle" ]->Fill(m1, m2, weight);
      bbhists2D[ sel_name + "_bb_dphi_reco_particle" ]->Fill(dphi1, dphi2, weight);
      bbhists2D[ sel_name + "_bb_deta_reco_particle" ]->Fill(deta1, deta2, weight);
      bbhists2D[ sel_name + "_bb_dr_reco_particle" ]->Fill(dr1, dr2, weight);

    } else if ( analysis_type1 == "particle" && analysis_type2 == "parton" ){

      bbhists2D[ sel_name + "_bb_pT_particle_parton" ]->Fill(pT1, pT2, weight);
      bbhists2D[ sel_name + "_bb_m_particle_parton" ]->Fill(m1, m2, weight);
      bbhists2D[ sel_name + "_bb_dphi_particle_parton" ]->Fill(dphi1, dphi2, weight);
      bbhists2D[ sel_name + "_bb_deta_particle_parton" ]->Fill(deta1, deta2, weight);
      bbhists2D[ sel_name + "_bb_dr_particle_parton" ]->Fill(dr1, dr2, weight);

    }

  }


  void fill_jj( string analysis_type, double weight, TLorentzVector &j1, TLorentzVector &j2, string process_name ) {


    double pT = -99; double m = -99; double dphi = -99; double deta = -99; double dr = -99; 

    if ( !(j1.Pt() == 0 && j1.Eta() == 0 && j1.Phi() == 0) && !(j2.Pt() == 0 && j2.Eta() == 0 && j2.Phi() == 0) ) {

      pT = ( j1 + j2 ).Pt();
      m = ( j1 + j2 ).M();

      dphi = delta_phi( j1 , j2 );
      deta = abs(delta_eta( j1 , j2 ));
      dr = delta_r( j1 , j2 );

    }

    if ( analysis_type == "reco" ){

      jjhists[ sel_name + "_jj_pT_reco" ]->Fill(pT, weight);
      jjhists[ sel_name + "_jj_m_reco" ]->Fill(m, weight);
      jjhists[ sel_name + "_jj_dphi_reco" ]->Fill(dphi, weight);
      jjhists[ sel_name + "_jj_deta_reco" ]->Fill(deta, weight);
      jjhists[ sel_name + "_jj_dr_reco" ]->Fill(dr, weight);

    } else if ( analysis_type == "particle" ){

      jjhists[ sel_name + "_jj_pT_particle" ]->Fill(pT, weight);
      jjhists[ sel_name + "_jj_m_particle" ]->Fill(m, weight);
      jjhists[ sel_name + "_jj_dphi_particle" ]->Fill(dphi, weight);
      jjhists[ sel_name + "_jj_deta_particle" ]->Fill(deta, weight);
      jjhists[ sel_name + "_jj_dr_particle" ]->Fill(dr, weight);

    } else if ( analysis_type == "parton" ){

      jjhists[ sel_name + "_jj_pT_parton" ]->Fill(pT, weight);
      jjhists[ sel_name + "_jj_m_parton" ]->Fill(m, weight);
      jjhists[ sel_name + "_jj_dphi_parton" ]->Fill(dphi, weight);
      jjhists[ sel_name + "_jj_deta_parton" ]->Fill(deta, weight);
      jjhists[ sel_name + "_jj_dr_parton" ]->Fill(dr, weight);

    }

  }

  void fill_jj_2D( string analysis_type1, string analysis_type2, double weight, TLorentzVector &at1_j1, TLorentzVector &at1_j2, TLorentzVector &at2_j1, TLorentzVector &at2_j2 ) {

    double pT1 = -99; double m1 = -99; double dphi1 = -99; double deta1 = -99; double dr1 = -99;
    double pT2 = -99; double m2 = -99; double dphi2 = -99; double deta2 = -99; double dr2 = -99;

    if ( !(at1_j1.Pt() == 0 && at1_j1.Eta() == 0 && at1_j1.Phi() == 0) && !(at1_j2.Pt() == 0 && at1_j2.Eta() == 0 && at1_j2.Phi() == 0) && !(at2_j1.Pt() == 0 && at2_j1.Eta() == 0 && at2_j1.Phi() == 0) && !(at2_j2.Pt() == 0 && at2_j2.Eta() == 0 && at2_j2.Phi() == 0) ) {

    pT1 = ( at1_j1 + at1_j2 ).Pt(); pT2 = ( at2_j1 + at2_j2 ).Pt();
    m1 = ( at1_j1 + at1_j2 ).M(); m2 = ( at2_j1 + at2_j2 ).M();

    dphi1 = delta_phi( at1_j1 , at1_j2 ); dphi2 = delta_phi( at2_j1 , at2_j2 );
    deta1 = abs(delta_eta( at1_j1 , at1_j2 )); deta2 = abs(delta_eta( at2_j1 , at2_j2 ));
    dr1 = delta_r( at1_j1 , at1_j2 ); dr2 = delta_r( at2_j1 , at2_j2 );

    }

    if ( analysis_type1 == "reco" && analysis_type2 == "particle" ){

      jjhists2D[ sel_name + "_jj_pT_reco_particle" ]->Fill(pT1, pT2, weight);
      jjhists2D[ sel_name + "_jj_m_reco_particle" ]->Fill(m1, m2, weight);
      jjhists2D[ sel_name + "_jj_dphi_reco_particle" ]->Fill(dphi1, dphi2, weight);
      jjhists2D[ sel_name + "_jj_deta_reco_particle" ]->Fill(deta1, deta2, weight);
      jjhists2D[ sel_name + "_jj_dr_reco_particle" ]->Fill(dr1, dr2, weight);

    } else if ( analysis_type1 == "particle" && analysis_type2 == "parton" ){

      jjhists2D[ sel_name + "_jj_pT_particle_parton" ]->Fill(pT1, pT2, weight);
      jjhists2D[ sel_name + "_jj_m_particle_parton" ]->Fill(m1, m2, weight);
      jjhists2D[ sel_name + "_jj_dphi_particle_parton" ]->Fill(dphi1, dphi2, weight);
      jjhists2D[ sel_name + "_jj_deta_particle_parton" ]->Fill(deta1, deta2, weight);
      jjhists2D[ sel_name + "_jj_dr_particle_parton" ]->Fill(dr1, dr2, weight);

    }

  }


  void fill_ww( string analysis_type, double weight, int et, TLorentzVector &w1, TLorentzVector &w2, string process_name ) {

    double pT1 = -99; double m1 = -99; double mT1 = -99; 
    double pT2 = -99; double m2 = -99; double mT2 = -99;
    double pT = -99; double m = -99; double mT = -99;

    double dphi = -99; double deta = -99; double dr = -99; 

    if ( !(w1.Pt() == 0 && w1.Eta() == 0 && w1.Phi() == 0) && !(w2.Pt() == 0 && w2.Eta() == 0 && w2.Phi() == 0) ) {

      pT1 = w1.Pt(); pT2 = w2.Pt(); pT = ( w1 + w2 ).Pt();

      m1 = w1.M(); m2 = w2.M(); m = ( w1 + w2 ).M();
      mT1 = w1.Mt(); mT2 = w2.Mt(); mT = ( w1 + w2 ).Mt();

      dphi = delta_phi( w1 , w2 );
      deta = delta_eta( w1 , w2 );
      dr = delta_r( w1 , w2 );

    }

    if ( analysis_type == "reco" ){

      wwhists[ sel_name + "_w1_pT_reco" ]->Fill(pT1, weight);
      wwhists[ sel_name + "_w1_m_reco" ]->Fill(m1, weight); wwhists[ sel_name + "_w1_mT_reco" ]->Fill(mT1, weight);
      wwhists[ sel_name + "_w2_pT_reco" ]->Fill(pT2, weight);
      wwhists[ sel_name + "_w2_m_reco" ]->Fill(m2, weight); wwhists[ sel_name + "_w2_mT_reco" ]->Fill(mT2, weight);
      wwhists[ sel_name + "_ww_pT_reco" ]->Fill(pT, weight);
      wwhists[ sel_name + "_ww_m_reco" ]->Fill(m, weight); wwhists[ sel_name + "_ww_mT_reco" ]->Fill(mT, weight);

      wwhists[ sel_name + "_ww_dphi_reco" ]->Fill(dphi, weight);
      wwhists[ sel_name + "_ww_deta_reco" ]->Fill(deta, weight);
      wwhists[ sel_name + "_ww_dr_reco" ]->Fill(dr, weight);

      wwhists[ sel_name + "_ww_et_reco" ]->Fill(et, weight);

    } else if ( analysis_type == "particle" ){

      wwhists[ sel_name + "_w1_pT_particle" ]->Fill(pT1, weight);
      wwhists[ sel_name + "_w1_m_particle" ]->Fill(m1, weight); wwhists[ sel_name + "_w1_mT_particle" ]->Fill(mT1, weight);
      wwhists[ sel_name + "_w2_pT_particle" ]->Fill(pT2, weight);
      wwhists[ sel_name + "_w2_m_particle" ]->Fill(m2, weight); wwhists[ sel_name + "_w2_mT_particle" ]->Fill(mT2, weight);
      wwhists[ sel_name + "_ww_pT_particle" ]->Fill(pT, weight);
      wwhists[ sel_name + "_ww_m_particle" ]->Fill(m, weight); wwhists[ sel_name + "_ww_mT_particle" ]->Fill(mT, weight);

      wwhists[ sel_name + "_ww_dphi_particle" ]->Fill(dphi, weight);
      wwhists[ sel_name + "_ww_deta_particle" ]->Fill(deta, weight);
      wwhists[ sel_name + "_ww_dr_particle" ]->Fill(dr, weight);

      wwhists[ sel_name + "_ww_et_particle" ]->Fill(et, weight);

    } else if ( analysis_type == "parton" ){

      wwhists[ sel_name + "_w1_pT_parton" ]->Fill(pT1, weight);
      wwhists[ sel_name + "_w1_m_parton" ]->Fill(m1, weight); wwhists[ sel_name + "_w1_mT_parton" ]->Fill(mT1, weight);
      wwhists[ sel_name + "_w2_pT_parton" ]->Fill(pT2, weight);
      wwhists[ sel_name + "_w2_m_parton" ]->Fill(m2, weight); wwhists[ sel_name + "_w2_mT_parton" ]->Fill(mT2, weight);
      wwhists[ sel_name + "_ww_pT_parton" ]->Fill(pT, weight);
      wwhists[ sel_name + "_ww_m_parton" ]->Fill(m, weight); wwhists[ sel_name + "_ww_mT_parton" ]->Fill(mT, weight);

      wwhists[ sel_name + "_ww_dphi_parton" ]->Fill(dphi, weight);
      wwhists[ sel_name + "_ww_deta_parton" ]->Fill(deta, weight);
      wwhists[ sel_name + "_ww_dr_parton" ]->Fill(dr, weight);

      wwhists[ sel_name + "_ww_et_parton" ]->Fill(et, weight);

    }

  }


  void fill_ww_2D( string analysis_type1, string analysis_type2, double weight, TLorentzVector &at1_w1, TLorentzVector &at1_w2, TLorentzVector &at2_w1, TLorentzVector &at2_w2 ) {

    double at1_pT1 = -99; double at2_pT1 = -99;
    double at1_pT2 = -99; double at2_pT2 = -99;
    double at1_pT = -99; double at2_pT = -99;

    double at1_m1 = -99; double at2_m1 = -99;
    double at1_m2 = -99; double at2_m2 = -99;
    double at1_m = -99; double at2_m = -99;

    double at1_mT1 = -99; double at2_mT1 = -99;
    double at1_mT2 = -99; double at2_mT2 = -99;
    double at1_mT = -99; double at2_mT = -99;

    double dphi1 = -99; double dphi2 = -99;
    double deta1 = -99; double deta2 = -99;
    double dr1 = -99; double dr2 = -99;

    if (!(at1_w1.Pt() == 0 && at1_w1.Eta() == 0 && at1_w1.Phi() == 0) && !(at1_w2.Pt() == 0 && at1_w2.Eta() == 0 && at1_w2.Phi() == 0) && 
        !(at2_w1.Pt() == 0 && at2_w1.Eta() == 0 && at2_w1.Phi() == 0) && !(at2_w2.Pt() == 0 && at2_w2.Eta() == 0 && at2_w2.Phi() == 0)) {

      at1_pT1 = at1_w1.Pt(); at2_pT1 = at2_w1.Pt();
      at1_pT2 = at1_w2.Pt(); at2_pT2 = at2_w2.Pt();
      at1_pT = (at1_w1 + at1_w2).Pt(); at2_pT = (at2_w1 + at2_w2).Pt();

      at1_m1 = at1_w1.M(); at2_m1 = at2_w1.M();
      at1_m2 = at1_w2.M(); at2_m2 = at2_w2.M();
      at1_m = (at1_w1 + at1_w2).M(); at2_m = (at2_w1 + at2_w2).M();

      at1_mT1 = at1_w1.Mt(); at2_mT1 = at2_w1.Mt();
      at1_mT2 = at1_w2.Mt(); at2_mT2 = at2_w2.Mt();
      at1_mT = (at1_w1 + at1_w2).Mt(); at2_mT = (at2_w1 + at2_w2).Mt();

      dphi1 = delta_phi(at1_w1, at1_w2); dphi2 = delta_phi(at2_w1, at2_w2);
      deta1 = delta_eta(at1_w1, at1_w2); deta2 = delta_eta(at2_w1, at2_w2);
      dr1 = delta_r(at1_w1, at1_w2); dr2 = delta_r(at2_w1, at2_w2);

    }


    if ( analysis_type1 == "reco" && analysis_type2 == "particle" ){

      wwhists2D[ sel_name + "_w1_pT_reco_particle" ]->Fill(at1_pT1, at2_pT2, weight);
      wwhists2D[ sel_name + "_w1_m_reco_particle" ]->Fill(at1_m1, at2_m2, weight); wwhists2D[ sel_name + "_w1_mT_reco_particle" ]->Fill(at1_mT1, at2_mT2, weight); 
      wwhists2D[ sel_name + "_w2_pT_reco_particle" ]->Fill(at1_pT1, at2_pT2, weight);
      wwhists2D[ sel_name + "_w2_m_reco_particle" ]->Fill(at1_m1, at2_m2, weight); wwhists2D[ sel_name + "_w2_mT_reco_particle" ]->Fill(at1_mT1, at2_mT2, weight);
      wwhists2D[ sel_name + "_ww_pT_reco_particle" ]->Fill(at1_pT, at2_pT, weight);
      wwhists2D[ sel_name + "_ww_m_reco_particle" ]->Fill(at1_m1, at2_m2, weight); wwhists2D[ sel_name + "_ww_mT_reco_particle" ]->Fill(at1_mT1, at2_mT2, weight);
      wwhists2D[ sel_name + "_ww_dphi_reco_particle" ]->Fill(dphi1, dphi2, weight);
      wwhists2D[ sel_name + "_ww_deta_reco_particle" ]->Fill(deta1, deta2, weight);
      wwhists2D[ sel_name + "_ww_dr_reco_particle" ]->Fill(dr1, dr2, weight);

    } else if ( analysis_type1 == "particle" && analysis_type2 == "parton" ){

      wwhists2D[ sel_name + "_w1_pT_particle_parton" ]->Fill(at1_pT1, at2_pT2, weight);
      wwhists2D[ sel_name + "_w1_m_particle_parton" ]->Fill(at1_m1, at2_m2, weight); wwhists2D[ sel_name + "_w1_mT_particle_parton" ]->Fill(at1_mT1, at2_mT2, weight);
      wwhists2D[ sel_name + "_w2_pT_particle_parton" ]->Fill(at1_pT1, at2_pT2, weight);
      wwhists2D[ sel_name + "_w2_m_particle_parton" ]->Fill(at1_m1, at2_m2, weight); wwhists2D[ sel_name + "_w2_mT_particle_parton" ]->Fill(at1_mT1, at2_mT2, weight);
      wwhists2D[ sel_name + "_ww_pT_particle_parton" ]->Fill(at1_pT1, at2_pT2, weight);
      wwhists2D[ sel_name + "_ww_m_particle_parton" ]->Fill(at1_m1, at2_m2, weight); wwhists2D[ sel_name + "_ww_mT_particle_parton" ]->Fill(at1_mT1, at2_mT2, weight);
      wwhists2D[ sel_name + "_ww_dphi_particle_parton" ]->Fill(dphi1, dphi2, weight);
      wwhists2D[ sel_name + "_ww_deta_particle_parton" ]->Fill(deta1, deta2, weight);
      wwhists2D[ sel_name + "_ww_dr_particle_parton" ]->Fill(dr1, dr2, weight);

    }

  }


  void fill_zz( string analysis_type, double weight, int et, TLorentzVector &z1, TLorentzVector &z2, string process_name ) {

    double pT1 = -99; double m1 = -99; 
    double pT2 = -99; double m2 = -99;
    double pT = -99; double m = -99;

    double dphi = -99; double deta = -99; double dr = -99; 

    if ( !(z1.Pt() == 0 && z1.Eta() == 0 && z1.Phi() == 0) && !(z2.Pt() == 0 && z2.Eta() == 0 && z2.Phi() == 0) ) {

      pT1 = z1.Pt(); pT2 = z2.Pt(); pT = ( z1 + z2 ).Pt();
      m1 = z1.M(); m2 = z2.M(); m = ( z1 + z2 ).M();

      dphi = delta_phi( z1 , z2 );
      deta = delta_eta( z1 , z2 );
      dr = delta_r( z1 , z2 );

    }

    if ( analysis_type == "reco" ){

      zzhists[ sel_name + "_z1_pT_reco"]->Fill(pT1, weight);
      zzhists[ sel_name + "_z1_m_reco"]->Fill(m1, weight);
      zzhists[ sel_name + "_z2_pT_reco"]->Fill(pT2, weight);
      zzhists[ sel_name + "_z2_m_reco"]->Fill(m2, weight);
      zzhists[ sel_name + "_zz_pT_reco"]->Fill(pT, weight);
      zzhists[ sel_name + "_zz_m_reco"]->Fill(m, weight);

      zzhists[ sel_name + "_zz_dphi_reco"]->Fill(dphi, weight);
      zzhists[ sel_name + "_zz_deta_reco"]->Fill(deta, weight);
      zzhists[ sel_name + "_zz_dr_reco"]->Fill(dr, weight);

      zzhists[ sel_name + "_zz_et_reco" ]->Fill(et, weight);

    } else if ( analysis_type == "particle" ){

      zzhists[ sel_name + "_z1_pT_particle" ]->Fill(pT1, weight);
      zzhists[ sel_name + "_z1_m_particle" ]->Fill(m1, weight);
      zzhists[ sel_name + "_z2_pT_particle" ]->Fill(pT2, weight);
      zzhists[ sel_name + "_z2_m_particle" ]->Fill(m2, weight);
      zzhists[ sel_name + "_zz_pT_particle" ]->Fill(pT, weight);
      zzhists[ sel_name + "_zz_m_particle" ]->Fill(m, weight);

      zzhists[ sel_name + "_zz_dphi_particle" ]->Fill(dphi, weight);
      zzhists[ sel_name + "_zz_deta_particle" ]->Fill(deta, weight);
      zzhists[ sel_name + "_zz_dr_particle" ]->Fill(dr, weight);

      zzhists[ sel_name + "_zz_et_particle" ]->Fill(et, weight);

    } else if ( analysis_type == "parton" ){

      zzhists[ sel_name + "_z1_pT_parton" ]->Fill(pT1, weight);
      zzhists[ sel_name + "_z1_m_parton" ]->Fill(m1, weight);
      zzhists[ sel_name + "_z2_pT_parton" ]->Fill(pT2, weight);
      zzhists[ sel_name + "_z2_m_parton" ]->Fill(m2, weight);
      zzhists[ sel_name + "_zz_pT_parton" ]->Fill(pT, weight);
      zzhists[ sel_name + "_zz_m_parton" ]->Fill(m, weight);

      zzhists[ sel_name + "_zz_dphi_parton" ]->Fill(dphi, weight);
      zzhists[ sel_name + "_zz_deta_parton" ]->Fill(deta, weight);
      zzhists[ sel_name + "_zz_dr_parton" ]->Fill(dr, weight);

      zzhists[ sel_name + "_zz_et_parton" ]->Fill(et, weight);

    }

  }


  void fill_zz_2D( string analysis_type1, string analysis_type2, double weight, TLorentzVector &at1_z1, TLorentzVector &at1_z2, TLorentzVector &at2_z1, TLorentzVector &at2_z2 ) {

    double at1_pT1 = -99; double at2_pT1 = -99;
    double at1_pT2 = -99; double at2_pT2 = -99;
    double at1_pT = -99; double at2_pT = -99;

    double at1_m1 = -99; double at2_m1 = -99;
    double at1_m2 = -99; double at2_m2 = -99;
    double at1_m = -99; double at2_m = -99;

    double dphi1 = -99; double dphi2 = -99;
    double deta1 = -99; double deta2 = -99;
    double dr1 = -99; double dr2 = -99;

    if (!(at1_z1.Pt() == 0 && at1_z1.Eta() == 0 && at1_z1.Phi() == 0) && !(at1_z2.Pt() == 0 && at1_z2.Eta() == 0 && at1_z2.Phi() == 0) && 
        !(at2_z1.Pt() == 0 && at2_z1.Eta() == 0 && at2_z1.Phi() == 0) && !(at2_z2.Pt() == 0 && at2_z2.Eta() == 0 && at2_z2.Phi() == 0)) {

      at1_pT1 = at1_z1.Pt(); at2_pT1 = at2_z1.Pt();
      at1_pT2 = at1_z2.Pt(); at2_pT2 = at2_z2.Pt();
      at1_pT = (at1_z1 + at1_z2).Pt(); at2_pT = (at2_z1 + at2_z2).Pt();

      at1_m1 = at1_z1.M(); at2_m1 = at2_z1.M();
      at1_m2 = at1_z2.M(); at2_m2 = at2_z2.M();
      at1_m = (at1_z1 + at1_z2).M(); at2_m = (at2_z1 + at2_z2).M();

      dphi1 = delta_phi(at1_z1, at1_z2); dphi2 = delta_phi(at2_z1, at2_z2);
      deta1 = delta_eta(at1_z1, at1_z2); deta2 = delta_eta(at2_z1, at2_z2);
      dr1 = delta_r(at1_z1, at1_z2); dr2 = delta_r(at2_z1, at2_z2);
      
    }

    if ( analysis_type1 == "reco" && analysis_type2 == "particle" ){

      zzhists2D[ sel_name + "_z1_pT_reco_particle" ]->Fill(at1_pT1, at2_pT2, weight);
      zzhists2D[ sel_name + "_z1_m_reco_particle" ]->Fill(at1_m1, at2_m2, weight); 
      zzhists2D[ sel_name + "_z2_pT_reco_particle" ]->Fill(at1_pT1, at2_pT2, weight);
      zzhists2D[ sel_name + "_z2_m_reco_particle" ]->Fill(at1_m1, at2_m2, weight);
      zzhists2D[ sel_name + "_zz_pT_reco_particle" ]->Fill(at1_pT, at2_pT, weight);
      zzhists2D[ sel_name + "_zz_m_reco_particle" ]->Fill(at1_m1, at2_m2, weight);
      zzhists2D[ sel_name + "_zz_dphi_reco_particle" ]->Fill(dphi1, dphi2, weight);
      zzhists2D[ sel_name + "_zz_deta_reco_particle" ]->Fill(deta1, deta2, weight);
      zzhists2D[ sel_name + "_zz_dr_reco_particle" ]->Fill(dr1, dr2, weight);

    } else if ( analysis_type1 == "particle" && analysis_type2 == "parton" ){

      zzhists2D[ sel_name + "_z1_pT_particle_parton" ]->Fill(at1_pT1, at2_pT2, weight);
      zzhists2D[ sel_name + "_z1_m_particle_parton" ]->Fill(at1_m1, at2_m2, weight);
      zzhists2D[ sel_name + "_z2_pT_particle_parton" ]->Fill(at1_pT1, at2_pT2, weight);
      zzhists2D[ sel_name + "_z2_m_particle_parton" ]->Fill(at1_m1, at2_m2, weight);
      zzhists2D[ sel_name + "_zz_pT_particle_parton" ]->Fill(at1_pT1, at2_pT2, weight);
      zzhists2D[ sel_name + "_zz_m_particle_parton" ]->Fill(at1_m1, at2_m2, weight);
      zzhists2D[ sel_name + "_zz_dphi_particle_parton" ]->Fill(dphi1, dphi2, weight);
      zzhists2D[ sel_name + "_zz_deta_particle_parton" ]->Fill(deta1, deta2, weight);
      zzhists2D[ sel_name + "_zz_dr_particle_parton" ]->Fill(dr1, dr2, weight);

    }

  }


  void fill_2l( string analysis_type, double weight, TLorentzVector &l1, TLorentzVector &l2, string process_name ) {

    double pT1 = -99; double eta1 = -99; double pT2 = -99; double eta2 = -99; 
    double pT = -99; double m = -99;
    double dphi = -99; double deta = -99; double dr = -99; 

    if ( !(l1.Pt() == 0 && l1.Eta() == 0 && l1.Phi() == 0) && !(l2.Pt() == 0 && l2.Eta() == 0 && l2.Phi() == 0) ) {

      pT1 = l1.Pt(); pT2 = l2.Pt();
      eta1 = l1.Eta(); eta2 = l2.Eta();

      pT = ( l1 + l2 ).Pt(); m = ( l1 + l2 ).M();

      dphi = delta_phi( l1 , l2 );
      deta = delta_eta( l1 , l2 );
      dr = delta_r( l1 , l2 );

    }

    if ( analysis_type == "reco" ){

      lephists[ sel_name + "_l1_pT_reco"]->Fill(pT1, weight);
      lephists[ sel_name + "_l2_pT_reco"]->Fill(pT2, weight);
      lephists[ sel_name + "_l1_eta_reco"]->Fill(eta1, weight);
      lephists[ sel_name + "_l2_eta_reco"]->Fill(eta2, weight);
      lephists[ sel_name + "_l1l2_pT_reco"]->Fill(pT, weight);
      lephists[ sel_name + "_l1l2_m_reco"]->Fill(m, weight);

      lephists[ sel_name + "_l1l2_dphi_reco"]->Fill(dphi, weight);
      lephists[ sel_name + "_l1l2_deta_reco"]->Fill(deta, weight);
      lephists[ sel_name + "_l1l2_dr_reco"]->Fill(dr, weight);

    } else if ( analysis_type == "particle" ){

      lephists[ sel_name + "_l1_pT_particle" ]->Fill(pT1, weight);
      lephists[ sel_name + "_l2_pT_particle" ]->Fill(pT2, weight);
      lephists[ sel_name + "_l1_eta_particle"]->Fill(eta1, weight);
      lephists[ sel_name + "_l2_eta_particle"]->Fill(eta2, weight);
      lephists[ sel_name + "_l1l2_pT_particle" ]->Fill(pT, weight);
      lephists[ sel_name + "_l1l2_m_particle" ]->Fill(m, weight);

      lephists[ sel_name + "_l1l2_dphi_particle" ]->Fill(dphi, weight);
      lephists[ sel_name + "_l1l2_deta_particle" ]->Fill(deta, weight);
      lephists[ sel_name + "_l1l2_dr_particle" ]->Fill(dr, weight);

    } else if ( analysis_type == "parton" ){

      lephists[ sel_name + "_l1_pT_parton" ]->Fill(pT1, weight);
      lephists[ sel_name + "_l2_pT_parton" ]->Fill(pT2, weight);
      lephists[ sel_name + "_l1_eta_parton"]->Fill(eta1, weight);
      lephists[ sel_name + "_l2_eta_parton"]->Fill(eta2, weight);
      lephists[ sel_name + "_l1l2_pT_parton" ]->Fill(pT, weight);
      lephists[ sel_name + "_l1l2_m_parton" ]->Fill(m, weight);

      lephists[ sel_name + "_l1l2_dphi_parton" ]->Fill(dphi, weight);
      lephists[ sel_name + "_l1l2_deta_parton" ]->Fill(deta, weight);
      lephists[ sel_name + "_l1l2_dr_parton" ]->Fill(dr, weight);

    }

  }


  void fill_2l_2D( string analysis_type1, string analysis_type2, double weight, TLorentzVector &at1_l1, TLorentzVector &at1_l2, TLorentzVector &at2_l1, TLorentzVector &at2_l2 ) {

    double at1_pT1 = -99; double at2_pT1 = -99;
    double at1_pT2 = -99; double at2_pT2 = -99;
    double at1_eta1 = at1_l1.Eta(); double at2_eta1 = at2_l1.Eta();
    double at1_eta2 = at1_l2.Eta(); double at2_eta2 = at2_l2.Eta();

    double at1_pT = -99; double at2_pT = -99;
    double at1_m = -99; double at2_m = -99;

    double dphi1 = -99; double dphi2 = -99;
    double deta1 = -99; double deta2 = -99;
    double dr1 = -99; double dr2 = -99;

    if (!(at1_l1.Pt() == 0 && at1_l1.Eta() == 0 && at1_l1.Phi() == 0) && !(at1_l2.Pt() == 0 && at1_l2.Eta() == 0 && at1_l2.Phi() == 0) && 
        !(at2_l1.Pt() == 0 && at2_l1.Eta() == 0 && at2_l1.Phi() == 0) && !(at2_l2.Pt() == 0 && at2_l2.Eta() == 0 && at2_l2.Phi() == 0)) {

      at1_pT1 = at1_l1.Pt(); at2_pT1 = at2_l1.Pt();
      at1_pT2 = at1_l2.Pt(); at2_pT2 = at2_l2.Pt();
      at1_eta1 = at1_l1.Eta(); at2_eta1 = at2_l1.Eta();
      at1_eta2 = at1_l2.Eta(); at2_eta2 = at2_l2.Eta();

      at1_pT = (at1_l1 + at1_l2).Pt(); at2_pT = (at2_l1 + at2_l2).Pt();
      at1_m = (at1_l1 + at1_l2).M(); at2_m = (at2_l1 + at2_l2).M();

      dphi1 = delta_phi(at1_l1, at1_l2); dphi2 = delta_phi(at2_l1, at2_l2);
      deta1 = delta_eta(at1_l1, at1_l2); deta2 = delta_eta(at2_l1, at2_l2);
      dr1 = delta_r(at1_l1, at1_l2); dr2 = delta_r(at2_l1, at2_l2);
      
    }

    if ( analysis_type1 == "reco" && analysis_type2 == "particle" ){

      lephists2D[ sel_name + "_l1_pT_reco_particle" ]->Fill(at1_pT1, at2_pT2, weight);
      lephists2D[ sel_name + "_l2_pT_reco_particle" ]->Fill(at1_pT1, at2_pT2, weight);
      lephists2D[ sel_name + "_l1_eta_reco_particle" ]->Fill(at1_eta1, at2_eta1, weight);
      lephists2D[ sel_name + "_l2_eta_reco_particle" ]->Fill(at1_eta2, at2_eta2, weight);
      lephists2D[ sel_name + "_l1l2_pT_reco_particle" ]->Fill(at1_pT, at2_pT, weight);
      lephists2D[ sel_name + "_l1l2_m_reco_particle" ]->Fill(at1_m, at2_m, weight);
      lephists2D[ sel_name + "_l1l2_dphi_reco_particle" ]->Fill(dphi1, dphi2, weight);
      lephists2D[ sel_name + "_l1l2_deta_reco_particle" ]->Fill(deta1, deta2, weight);
      lephists2D[ sel_name + "_l1l2_dr_reco_particle" ]->Fill(dr1, dr2, weight);

    } else if ( analysis_type1 == "particle" && analysis_type2 == "parton" ){

      lephists2D[ sel_name + "_l1_pT_particle_parton" ]->Fill(at1_pT1, at2_pT2, weight);
      lephists2D[ sel_name + "_l2_pT_particle_parton" ]->Fill(at1_pT1, at2_pT2, weight);
      lephists2D[ sel_name + "_l1_eta_particle_parton" ]->Fill(at1_eta1, at2_eta1, weight);
      lephists2D[ sel_name + "_l2_eta_particle_parton" ]->Fill(at1_eta2, at2_eta2, weight);
      lephists2D[ sel_name + "_l1l2_pT_particle_parton" ]->Fill(at1_pT, at2_pT, weight);
      lephists2D[ sel_name + "_l1l2_m_particle_parton" ]->Fill(at1_m, at2_m, weight);
      lephists2D[ sel_name + "_l1l2_dphi_particle_parton" ]->Fill(dphi1, dphi2, weight);
      lephists2D[ sel_name + "_l1l2_deta_particle_parton" ]->Fill(deta1, deta2, weight);
      lephists2D[ sel_name + "_l1l2_dr_particle_parton" ]->Fill(dr1, dr2, weight);

    }

  }

  void fill_4l( string analysis_type, double weight, TLorentzVector &l3, TLorentzVector &l4, string process_name ) {

    double pT3 = -99; double eta3 = -99; double pT4 = -99; double eta4 = -99; 
    double pT = -99; double m = -99;
    double dphi = -99; double deta = -99; double dr = -99; 

    if ( !(l3.Pt() == 0 && l3.Eta() == 0 && l3.Phi() == 0) && !(l4.Pt() == 0 && l4.Eta() == 0 && l4.Phi() == 0) ) {

      pT3 = l3.Pt(); pT4 = l4.Pt();
      eta3 = l3.Eta(); eta4 = l4.Eta();

      pT = ( l3 + l4 ).Pt(); m = ( l3 + l4 ).M();

      dphi = delta_phi( l3 , l4 );
      deta = delta_eta( l3 , l4 );
      dr = delta_r( l3 , l4 );

    }

    if ( analysis_type == "reco" ){

      lephists[ sel_name + "_l3_pT_reco"]->Fill(pT3, weight);
      lephists[ sel_name + "_l4_pT_reco"]->Fill(pT4, weight);
      lephists[ sel_name + "_l3_eta_reco"]->Fill(eta3, weight);
      lephists[ sel_name + "_l4_eta_reco"]->Fill(eta4, weight);
      lephists[ sel_name + "_l3l4_pT_reco"]->Fill(pT, weight);
      lephists[ sel_name + "_l3l4_m_reco"]->Fill(m, weight);

      lephists[ sel_name + "_l3l4_dphi_reco"]->Fill(dphi, weight);
      lephists[ sel_name + "_l3l4_deta_reco"]->Fill(deta, weight);
      lephists[ sel_name + "_l3l4_dr_reco"]->Fill(dr, weight);

    } else if ( analysis_type == "particle" ){

      lephists[ sel_name + "_l3_pT_particle" ]->Fill(pT3, weight);
      lephists[ sel_name + "_l4_pT_particle" ]->Fill(pT4, weight);
      lephists[ sel_name + "_l3_eta_particle"]->Fill(eta3, weight);
      lephists[ sel_name + "_l4_eta_particle"]->Fill(eta4, weight);
      lephists[ sel_name + "_l3l4_pT_particle" ]->Fill(pT, weight);
      lephists[ sel_name + "_l3l4_m_particle" ]->Fill(m, weight);

      lephists[ sel_name + "_l3l4_dphi_particle" ]->Fill(dphi, weight);
      lephists[ sel_name + "_l3l4_deta_particle" ]->Fill(deta, weight);
      lephists[ sel_name + "_l3l4_dr_particle" ]->Fill(dr, weight);

    } else if ( analysis_type == "parton" ){

      lephists[ sel_name + "_l3_pT_parton" ]->Fill(pT3, weight);
      lephists[ sel_name + "_l4_pT_parton" ]->Fill(pT4, weight);
      lephists[ sel_name + "_l3_eta_parton"]->Fill(eta3, weight);
      lephists[ sel_name + "_l4_eta_parton"]->Fill(eta4, weight);
      lephists[ sel_name + "_l3l4_pT_parton" ]->Fill(pT, weight);
      lephists[ sel_name + "_l3l4_m_parton" ]->Fill(m, weight);

      lephists[ sel_name + "_l3l4_dphi_parton" ]->Fill(dphi, weight);
      lephists[ sel_name + "_l3l4_deta_parton" ]->Fill(deta, weight);
      lephists[ sel_name + "_l3l4_dr_parton" ]->Fill(dr, weight);

    }

  }


  void fill_4l_2D( string analysis_type1, string analysis_type2, double weight, TLorentzVector &at1_l3, TLorentzVector &at1_l4, TLorentzVector &at2_l3, TLorentzVector &at2_l4 ) {

    double at1_pT3 = -99; double at2_pT3 = -99;
    double at1_pT4 = -99; double at2_pT4 = -99;
    double at1_eta3 = at1_l3.Eta(); double at2_eta3 = at2_l3.Eta();
    double at1_eta4 = at1_l4.Eta(); double at2_eta4 = at2_l4.Eta();

    double at1_pT = -99; double at2_pT = -99;
    double at1_m = -99; double at2_m = -99;

    double dphi1 = -99; double dphi2 = -99;
    double deta1 = -99; double deta2 = -99;
    double dr1 = -99; double dr2 = -99;

    if (!(at1_l3.Pt() == 0 && at1_l3.Eta() == 0 && at1_l3.Phi() == 0) && !(at1_l4.Pt() == 0 && at1_l4.Eta() == 0 && at1_l4.Phi() == 0) && 
        !(at2_l3.Pt() == 0 && at2_l3.Eta() == 0 && at2_l3.Phi() == 0) && !(at2_l4.Pt() == 0 && at2_l4.Eta() == 0 && at2_l4.Phi() == 0)) {

      at1_pT3 = at1_l3.Pt(); at2_pT3 = at2_l3.Pt();
      at1_pT4 = at1_l4.Pt(); at2_pT4 = at2_l4.Pt();
      at1_eta3 = at1_l3.Eta(); at2_eta3 = at2_l3.Eta();
      at1_eta4 = at1_l4.Eta(); at2_eta4 = at2_l4.Eta();
      
      at1_pT = (at1_l3 + at1_l4).Pt(); at2_pT = (at2_l3 + at2_l4).Pt();
      at1_m = (at1_l3 + at1_l4).M(); at2_m = (at2_l3 + at2_l4).M();

      dphi1 = delta_phi(at1_l3, at1_l4); dphi2 = delta_phi(at2_l3, at2_l4);
      deta1 = delta_eta(at1_l3, at1_l4); deta2 = delta_eta(at2_l3, at2_l4);
      dr1 = delta_r(at1_l3, at1_l4); dr2 = delta_r(at2_l3, at2_l4);
      
    }

    if ( analysis_type1 == "reco" && analysis_type2 == "particle" ){

      lephists2D[ sel_name + "_l3_pT_reco_particle" ]->Fill(at1_pT3, at2_pT3, weight);
      lephists2D[ sel_name + "_l4_pT_reco_particle" ]->Fill(at1_pT4, at2_pT4, weight);
      lephists2D[ sel_name + "_l3_eta_reco_particle" ]->Fill(at1_eta3, at2_eta3, weight);
      lephists2D[ sel_name + "_l4_eta_reco_particle" ]->Fill(at1_eta4, at2_eta4, weight);
      lephists2D[ sel_name + "_l3l4_pT_reco_particle" ]->Fill(at1_pT, at2_pT, weight);
      lephists2D[ sel_name + "_l3l4_m_reco_particle" ]->Fill(at1_m, at2_m, weight);
      lephists2D[ sel_name + "_l3l4_dphi_reco_particle" ]->Fill(dphi1, dphi2, weight);
      lephists2D[ sel_name + "_l3l4_deta_reco_particle" ]->Fill(deta1, deta2, weight);
      lephists2D[ sel_name + "_l3l4_dr_reco_particle" ]->Fill(dr1, dr2, weight);

    } else if ( analysis_type1 == "particle" && analysis_type2 == "parton" ){

      lephists2D[ sel_name + "_l3_pT_particle_parton" ]->Fill(at1_pT3, at2_pT3, weight);
      lephists2D[ sel_name + "_l4_pT_particle_parton" ]->Fill(at1_pT4, at2_pT4, weight);
      lephists2D[ sel_name + "_l3_eta_particle_parton" ]->Fill(at1_eta3, at2_eta3, weight);
      lephists2D[ sel_name + "_l4_eta_particle_parton" ]->Fill(at1_eta4, at2_eta4, weight);
      lephists2D[ sel_name + "_l3l4_pT_particle_parton" ]->Fill(at1_pT, at2_pT, weight);
      lephists2D[ sel_name + "_l3l4_m_particle_parton" ]->Fill(at1_m, at2_m, weight);
      lephists2D[ sel_name + "_l3l4_dphi_particle_parton" ]->Fill(dphi1, dphi2, weight);
      lephists2D[ sel_name + "_l3l4_deta_particle_parton" ]->Fill(deta1, deta2, weight);
      lephists2D[ sel_name + "_l3l4_dr_particle_parton" ]->Fill(dr1, dr2, weight);

    }

  }

  void fill_met( string analysis_type, double weight, TLorentzVector &met, string process_name ) {

    double MET = -99;

    if (!(met.Et() == 0 && met.Eta() == 0 && met.Phi() == 0)) MET = met.Et();

    if ( analysis_type == "reco" ) methists[sel_name + "_MET_reco"]->Fill(MET, weight);
    if ( analysis_type == "particle" ) methists[sel_name + "_MET_particle"]->Fill(MET, weight);

  }

void write_all_hist() {

    for (auto& pair : bbhists) {
      // histogram->SetDirectory(0)
      // cd 
      pair.second->Write();
    }
    for (auto& pair : bbhists2D) {
      pair.second->Write();
    }

    for (auto& pair : jjhists) {
      pair.second->Write();
    }
    for (auto& pair : jjhists2D) {
      pair.second->Write();
    }

    for (auto& pair : wwhists) {
      pair.second->Write();
    }
    for (auto& pair : wwhists2D) {
      pair.second->Write();
    }

    for (auto& pair : zzhists) {
      pair.second->Write();
    }
    for (auto& pair : zzhists2D) {
      pair.second->Write();
    }

    for (auto& pair : lephists) {
      pair.second->Write();
    }
    for (auto& pair : lephists2D) {
      pair.second->Write();
    }
    for (auto& pair : methists) {
      pair.second->Write();
    }
}

~histograms() {

  for (auto& pair : bbhists) {
    // delete pair.second;
  }
  for (auto& pair : bbhists2D) {
    // delete pair.second;
  }
  for (auto& pair : jjhists) {
  }
  for (auto& pair : jjhists2D) {
  }
  for (auto& pair : wwhists) {
  }
  for (auto& pair : wwhists2D) {
  }  
  for (auto& pair : zzhists) {
  }
  for (auto& pair : zzhists2D) {
  }
  for (auto& pair : lephists) {
  }
  for (auto& pair : lephists2D) {
  }
  for (auto& pair : methists) {
  }
  
}

};


// OLD



/*
void initialize_histograms() {

  vector <TH1F*> listOfTH1;
  vector <TH2F*> listOfTH2;
  vector <TProfile*> listOfTProfiles;

  int mBins = 10, pTBins = 10, phiBins = 10,  etaBins = 10, rBins = 10, cosBins = 10;

  double hpTmin = 0, hpTmax = 500;
  double jpTmin = 0, jpTmax = 500;
  double zpTmin = 0, zpTmax = 500;
  double wpTmin = 0, wpTmax = 500;
  double lpTmin = 0, lpTmax = 500;
  double metpTmin = 0, metpTmax = 500;

  const double hmmin = 0;
  const double jmmin = 0;
  const double zmmin = 0;
  const double wmmin = 0;
  const double hmmax = 200;
  const double jmmax = 200;
  const double zmmax = 200;
  const double wmmax = 200;

  const double hetamin = -2.5;
  const double jetamin = -2.5;
  const double zetamin = -2.5;
  const double wetamin = -2.5;
  const double letamin = -2.5;
  const double METetamin = -2.5;
  const double hetamax = 2.5;
  const double jetamax = 2.5;
  const double zetamax = 2.5;
  const double wetamax = 2.5;
  const double letamax = 2.5;
  const double METetamax = 2.5;

  const double hRmin = 0;
  const double jRmin = 0;
  const double zRmin = 0;
  const double wRmin = 0;
  const double lRmin = 0;
  const double METRmin = 0;
  const double hRmax = 5;
  const double jRmax = 5;
  const double zRmax = 5;
  const double wRmax = 5;
  const double lRmax = 5;
  const double METRmax = 5;

// 1D

  TH1F *hWeight = new TH1F("weights", "weight", 50, 0.0, 1.0);  listOfTH1.push_back(hWeight);

  // higgs
    TH1F *hHpTreco = new TH1F("Hbb_pT_reco", "p^{T}_{Hbb} Reco", pTBins, hpTmin, hpTmax); 	listOfTH1.push_back(hHpTreco);
    TH1F *hHmreco = new TH1F("Hbb_m_reco", "m_{Hbb} Reco", mBins, hmmin, hmmax); listOfTH1.push_back(hHmreco);
    TH1F *hbbdeltaPhireco = new TH1F("bb_deltaPhi_reco", "#Delta#phi_{bb} Reco", phiBins, -TMath::Pi(), +TMath::Pi()); listOfTH1.push_back(hbbdeltaPhireco);
    TH1F *hbbdeltaEtareco = new TH1F("bb_deltaEta_reco", "#Delta#eta_{bb} Reco", etaBins, hetamin, hetamax); listOfTH1.push_back(hbbdeltaEtareco);
    TH1F *hbbdeltaRreco = new TH1F("bb_deltaR_reco", "#DeltaR_{bb} Reco", rBins, hRmin, hRmax); listOfTH1.push_back(hbbdeltaRreco);

    TH1F *hHpTparticle = new TH1F("Hbb_pT_particle", "p^{T}_{Hbb} Particle", pTBins, hpTmin, hpTmax); listOfTH1.push_back(hHpTparticle);
    TH1F *hHmparticle = new TH1F("Hbb_m_particle", "m_{Hbb} Particle", mBins, hmmin, hmmax); listOfTH1.push_back(hHmparticle);
    TH1F *hbbdeltaPhiparticle = new TH1F("bb_deltaPhi_particle", "#Delta#phi_{bb} Particle", phiBins, -TMath::Pi(), +TMath::Pi()); listOfTH1.push_back(hbbdeltaPhiparticle);
    TH1F *hbbdeltaEtaparticle = new TH1F("bb_deltaEta_particle", "#Delta#eta_{bb} Particle", etaBins, hetamin, hetamax); listOfTH1.push_back(hbbdeltaEtaparticle);
    TH1F *hbbdeltaRparticle = new TH1F("bb_deltaR_particle", "#DeltaR_{bb} Particle", rBins, hRmin, hRmax); listOfTH1.push_back(hbbdeltaRparticle);

    // TH1F *hHpTparton = new TH1F("Hbb_pT_parton", "p^{T}_{hbb} Parton", pTBins, hpTmin, hpTmax); listOfTH1.push_back(hHpTparton);
    // TH1F *hHmparton = new TH1F("Hbb_m_parton", "m_{hbb} Parton", mBins, hmmin,  hmmax); listOfTH1.push_back(hHmparton);
    // TH1F *hbbdeltaPhiparton = new TH1F("bb_deltaPhi_parton", "#Delta#phi_{bb} Parton", phiBins, -TMath::Pi(), +TMath::Pi()); listOfTH1.push_back(hbbdeltaPhiparton);
    // TH1F *hbbdeltaEtaparton = new TH1F("bb_deltaEta_parton", "#Delta#eta_{bb} Parton", etaBins, hetamin, hetamax); listOfTH1.push_back(hbbdeltaEtaparton);
    // TH1F *hbbdeltaRparton = new TH1F("bb_deltaR_parton", "#DeltaR_{bb} Parton", rBins, hRmin, hRmax); listOfTH1.push_back(hbbdeltaRparton);

  // vbf jets
    TH1F *hjjpTreco = new TH1F("jj_pT_reco", "p^{T}_{jj} Reco", pTBins, jpTmin, jpTmax); listOfTH1.push_back(hjjpTreco);
    TH1F *hjjdeltaPhireco = new TH1F("jj_deltaPhi_reco", "#Delta#phi_{jj} Reco", phiBins, -TMath::Pi(), +TMath::Pi()); listOfTH1.push_back(hjjdeltaPhireco);
    TH1F *hjjdeltaEtareco = new TH1F("jj_deltaEta_reco", "#Delta#eta_{jj} Reco", etaBins, jetamin, jetamax); listOfTH1.push_back(hjjdeltaEtareco);
    TH1F *hjjdeltaRreco = new TH1F("jj_deltaR_reco", "#DeltaR_{jj} Reco", rBins, jRmin, jRmax); listOfTH1.push_back(hjjdeltaRreco);

    TH1F *hjjpTparticle = new TH1F("jj_pT_particle", "p^{T}_{jj} Particle", pTBins, jpTmin, jpTmax); listOfTH1.push_back(hjjpTparticle);
    TH1F *hjjdeltaPhiparticle = new TH1F("jj_deltaPhi_particle", "#Delta#phi_{jj} Particle", phiBins, -TMath::Pi(), +TMath::Pi()); listOfTH1.push_back(hjjdeltaPhiparticle);
    TH1F *hjjdeltaEtaparticle = new TH1F("jj_deltaEta_particle", "#Delta#eta_{jj} Particle", etaBins, jetamin, jetamax); listOfTH1.push_back(hjjdeltaEtaparticle);
    TH1F *hjjdeltaRparticle = new TH1F("jj_deltaR_particle", "#DeltaR_{jj} Particle", rBins, jRmin, jRmax); listOfTH1.push_back(hjjdeltaRparticle);

  // z1
    TH1F *hZ1pTreco = new TH1F("Z1_pT_reco", "p^{T}_{Z1} Reco", pTBins, zpTmin, zpTmax); listOfTH1.push_back(hZ1pTreco);
    TH1F *hZ1mreco = new TH1F("Z1_m_reco", "m_{Z1} Reco", mBins, zmmin, zmmax); listOfTH1.push_back(hZ1mreco);

    TH1F *hZ1pTparticle = new TH1F("Z1_pT_particle", "p^{T}_{Z1} Particle", pTBins, zpTmin, zpTmax); listOfTH1.push_back(hZ1pTparticle);
    TH1F *hZ1mparticle = new TH1F("Z1_m_particle", "m_{Z1} Particle", mBins, zmmin, zmmax); listOfTH1.push_back(hZ1mparticle);

    // TH1F *hZ1pTparton = new TH1F("Z1_pT_parton", "p^{T}_{Z1} Parton", pTBins, zpTmin, zpTmax); listOfTH1.push_back(hZ1pTparton);
    // TH1F *hZ1mparton = new TH1F("Z1_m_parton", "m_{Z1} Parton", mBins, zmmin, zmmax); listOfTH1.push_back(hZ1mparton);

  // z2
    TH1F *hZ2pTreco = new TH1F("Z2_pT_reco", "p^{T}_{Z2} Reco", pTBins, zpTmin, zpTmax); listOfTH1.push_back(hZ2pTreco);
    TH1F *hZ2mreco = new TH1F("Z2_m_reco", "m_{Z2} Reco", mBins, zmmin, zmmax); listOfTH1.push_back(hZ2mreco);

    TH1F *hZ2pTparticle = new TH1F("Z2_pT_particle", "p^{T}_{Z2} Particle", pTBins, zpTmin, zpTmax); listOfTH1.push_back(hZ2pTparticle);
    TH1F *hZ2mparticle = new TH1F("Z2_m_particle", "m_{Z2} Particle", mBins, zmmin, zmmax); listOfTH1.push_back(hZ2mparticle);

    // TH1F *hZ2pTparton = new TH1F("Z2_pT_parton", "p^{T}_{Z2} Parton", pTBins, zpTmin, zpTmax); listOfTH1.push_back(hZ2pTparton);
    // TH1F *hZ2mparton = new TH1F("Z2_m_parton", "m_{Z2} Parton", mBins, zmmin, zmmax); listOfTH1.push_back(hZ2mparton);

  // zz
    TH1F *hZZdeltaPhireco = new TH1F("ZZ_#Delta#phi_reco", "#Delta#phi_{ZZ} Reco", phiBins, -TMath::Pi(), +TMath::Pi()); listOfTH1.push_back(hZZdeltaPhireco);
    TH1F *hZZdeltaEtareco = new TH1F("ZZ_#Delta#eta_reco", "#Delta#eta_{ZZ} Reco", etaBins, zetamin, zetamax);listOfTH1.push_back(hZZdeltaEtareco);
    TH1F *hZZdeltaRreco = new TH1F("ZZ_#DeltaR_reco", "#DeltaR_{ZZ} Reco", rBins, zRmin, zRmax); listOfTH1.push_back(hZZdeltaRreco);

    TH1F *hZZdeltaPhiparticle = new TH1F("ZZ_#Delta#phi_particle", "#Delta#phi_{ZZ} Particle", phiBins, -TMath::Pi(), +TMath::Pi()); listOfTH1.push_back(hZZdeltaPhiparticle);
    TH1F *hZZdeltaEtaparticle = new TH1F("ZZ_#Delta#eta_particle", "#Delta#eta_{ZZ} Particle", etaBins, zetamin, zetamax);listOfTH1.push_back(hZZdeltaEtaparticle);
    TH1F *hZZdeltaRparticle = new TH1F("ZZ_#DeltaR_particle", "#DeltaR_{ZZ} Particle", rBins, zRmin, zRmax); listOfTH1.push_back(hZZdeltaRparticle);

    // TH1F *hZZdeltaPhiparton = new TH1F("ZZ_#Delta#phi_parton", "#Delta#phi_{ZZ} Parton", phiBins, -TMath::Pi(), +TMath::Pi()); listOfTH1.push_back(hZZdeltaPhiparton);
    // TH1F *hZZdeltaEtaparton = new TH1F("ZZ_#Delta#eta_parton", "#Delta#eta_{ZZ} Parton", etaBins, zetamin, zetamax);listOfTH1.push_back(hZZdeltaEtaparton);
    // TH1F *hZZdeltaRparton = new TH1F("ZZ_#DeltaR_parton", "#DeltaR_{ZZ} Parton", rBins, zRmin, zRmax); listOfTH1.push_back(hZZdeltaRparton);

  // w1
    TH1F *hW1pTreco = new TH1F("W1_pT_reco", "p^{T}_{W1} Reco", pTBins, wpTmin, wpTmax); listOfTH1.push_back(hW1pTreco);
    TH1F *hW1mreco = new TH1F("W1_m_reco", "m_{W1} Reco", mBins, wmmin, wmmax); listOfTH1.push_back(hW1mreco);

    TH1F *hW1pTparticle = new TH1F("W1_pT_particle", "p^{T}_{W1} Particle", pTBins, wpTmin, wpTmax); listOfTH1.push_back(hW1pTparticle);
    TH1F *hW1mparticle = new TH1F("W1_m_particle", "m_{W1} Particle", mBins, wmmin, wmmax); listOfTH1.push_back(hW1mparticle);

    // TH1F *hW1pTparton = new TH1F("W1_pT_parton", "p^{T}_{W1} Parton", pTBins, wpTmin, wpTmax); listOfTH1.push_back(hW1pTparton);
    // TH1F *hW1mparton = new TH1F("W1_m_parton", "m_{W1} Parton", mBins, wmmin, wmmax); listOfTH1.push_back(hW1mparton);

  // w2
    TH1F *hW2pTreco = new TH1F("W2_pT_reco", "p^{T}_{W2} Reco", pTBins, wpTmin, wpTmax); listOfTH1.push_back(hW2pTreco);
    TH1F *hW2mreco = new TH1F("W2_m_reco", "m_{W2} Reco", mBins, wmmin, wmmax); listOfTH1.push_back(hW2mreco);

    TH1F *hW2pTparticle = new TH1F("W2_pT_particle", "p^{T}_{W2} Particle", pTBins, wpTmin, wpTmax); listOfTH1.push_back(hW2pTparticle);
    TH1F *hW2mparticle = new TH1F("W2_m_particle", "m_{W2} Particle", mBins, wmmin, wmmax); listOfTH1.push_back(hW2mparticle);

    // TH1F *hW2pTparton = new TH1F("W2_pT_parton", "p^{T}_{W2} Parton", pTBins, wpTmin, wpTmax); listOfTH1.push_back(hW2pTparton);
    // TH1F *hW2mparton = new TH1F("W2_m_parton", "m_{W2} Parton", mBins, wmmin, wmmax); listOfTH1.push_back(hW2mparton);

  // ww
    TH1F *hWWdeltaPhireco = new TH1F("WW_#Delta#phi_reco", "#Delta#phi_{WW} Reco", phiBins, -TMath::Pi(), +TMath::Pi()); listOfTH1.push_back(hWWdeltaPhireco);
    TH1F *hWWdeltaEtareco = new TH1F("WW_#Delta#eta_reco", "#Delta#eta_{WW} Reco", etaBins, wetamin, wetamax);listOfTH1.push_back(hWWdeltaEtareco);
    TH1F *hWWdeltaRreco = new TH1F("WW_#DeltaR_reco", "#DeltaR_{WW} Reco", rBins, wRmin, wRmax); listOfTH1.push_back(hWWdeltaRreco);

    TH1F *hWWdeltaPhiparticle = new TH1F("WW_#Delta#phi_particle", "#Delta#phi_{WW} Particle", phiBins, -TMath::Pi(), +TMath::Pi()); listOfTH1.push_back(hWWdeltaPhiparticle);
    TH1F *hWWdeltaEtaparticle = new TH1F("WW_#Delta#eta_particle", "#Delta#eta_{WW} Particle", etaBins, wetamin, wetamax);listOfTH1.push_back(hWWdeltaEtaparticle);
    TH1F *hWWdeltaRparticle = new TH1F("WW_#DeltaR_particle", "#DeltaR_{WW} Particle", rBins, wRmin, wRmax); listOfTH1.push_back(hWWdeltaRparticle);

    // TH1F *hWWdeltaPhiparton = new TH1F("WW_#Delta#phi_parton", "#Delta#phi_{WW} Parton", phiBins, -TMath::Pi(), +TMath::Pi()); listOfTH1.push_back(hWWdeltaPhiparton);
    // TH1F *hWWdeltaEtaparton = new TH1F("WW_#Delta#eta_parton", "#Delta#eta_{WW} Parton", etaBins, wetamin, wetamax);listOfTH1.push_back(hWWdeltaEtaparton);
    // TH1F *hWWdeltaRparton = new TH1F("WW_#DeltaR_parton", "#DeltaR_{WW} Parton", rBins, wRmin, wRmax); listOfTH1.push_back(hWWdeltaRparton);

  // leptons

    TH1F *lead_e_pt_reco = new TH1F("lead_e_pT_reco", "p^{T}_{e1}_reco", pTBins, lpTmin, lpTmax); listOfTH1.push_back(lead_e_pt_reco);
    TH1F *lead_e_eta_reco = new TH1F("lead_e_#eta_reco", "#eta_{e1}_reco", etaBins, letamin, letamax); listOfTH1.push_back(lead_e_eta_reco);
    TH1F *lead_e_phi_reco = new TH1F("lead_e_#phi_reco", "#phi_{e1}_reco", phiBins, -TMath::Pi(), TMath::Pi()); listOfTH1.push_back(lead_e_phi_reco);
    TH1F *sublead_e_pt_reco = new TH1F("sublead_e_pT_reco", "p^{T}_{e2}_reco", pTBins, lpTmin, lpTmax); listOfTH1.push_back(sublead_e_pt_reco);
    TH1F *sublead_e_eta_reco = new TH1F("sublead_e_#eta_reco", "#eta_{e2}_reco", etaBins, letamin, letamax); listOfTH1.push_back(sublead_e_eta_reco);
    TH1F *sublead_e_phi_reco = new TH1F("sublead_e_#phi_reco", "#phi_{e2}_reco", phiBins, -TMath::Pi(), TMath::Pi()); listOfTH1.push_back(sublead_e_phi_reco);
    TH1F *goodE_size_reco = new TH1F("goodE_size_reco", "size", 5, 0, 5); listOfTH1.push_back(goodE_size_reco);

    TH1F *lead_mu_pt_reco = new TH1F("lead_mu_pT_reco", "p^{T}_{mu1}_reco", pTBins, lpTmin, lpTmax); listOfTH1.push_back(lead_mu_pt_reco);
    TH1F *lead_mu_eta_reco = new TH1F("lead_mu_#eta_reco", "#eta_{mu1}_reco", etaBins, letamin, letamax); listOfTH1.push_back(lead_mu_eta_reco);
    TH1F *lead_mu_phi_reco = new TH1F("lead_mu_#phi_reco", "#phi_{mu1}_reco", phiBins, -TMath::Pi(), TMath::Pi()); listOfTH1.push_back(lead_mu_phi_reco);
    TH1F *sublead_mu_pt_reco = new TH1F("sublead_mu_pT_reco", "p^{T}_{mu2}_reco", pTBins, lpTmin, lpTmax); listOfTH1.push_back(sublead_mu_pt_reco);
    TH1F *sublead_mu_eta_reco = new TH1F("sublead_mu_#eta_reco", "#eta_{mu2}_reco", etaBins, letamin, letamax); listOfTH1.push_back(sublead_mu_eta_reco);
    TH1F *sublead_mu_phi_reco = new TH1F("sublead_mu_#phi_reco", "#phi_{mu2}_reco", phiBins, -TMath::Pi(), TMath::Pi()); listOfTH1.push_back(sublead_mu_phi_reco);
    TH1F *goodMu_size_reco = new TH1F("goodMu_size_reco", "size", 5, 0, 5); listOfTH1.push_back(goodMu_size_reco);

    TH1F *lead_e_pt_particle = new TH1F("lead_e_pT_particle", "p^{T}_{e1}_particle", pTBins, lpTmin, lpTmax); listOfTH1.push_back(lead_e_pt_particle);
    TH1F *lead_e_eta_particle = new TH1F("lead_e_#eta_particle", "#eta_{e1}_particle", etaBins, letamin, letamax); listOfTH1.push_back(lead_e_eta_particle);
    TH1F *lead_e_phi_particle = new TH1F("lead_e_#phi_particle", "#phi_{e1}_particle", phiBins, -TMath::Pi(), TMath::Pi()); listOfTH1.push_back(lead_e_phi_particle);
    TH1F *sublead_e_pt_particle = new TH1F("sublead_e_pT_particle", "p^{T}_{e2}_particle", pTBins, lpTmin, lpTmax); listOfTH1.push_back(sublead_e_pt_particle);
    TH1F *sublead_e_eta_particle = new TH1F("sublead_e_#eta_particle", "#eta_{e2}_particle", etaBins, letamin, letamax); listOfTH1.push_back(sublead_e_eta_particle);
    TH1F *sublead_e_phi_particle = new TH1F("sublead_e_#phi_particle", "#phi_{e2}_particle", phiBins, -TMath::Pi(), TMath::Pi()); listOfTH1.push_back(sublead_e_phi_particle);
    TH1F *goodE_size_particle = new TH1F("goodE_size_particle", "size", 5, 0, 5); listOfTH1.push_back(goodE_size_particle);

    TH1F *lead_mu_pt_particle = new TH1F("lead_mu_pT_particle", "p^{T}_{mu1}_particle", pTBins, lpTmin, lpTmax); listOfTH1.push_back(lead_mu_pt_particle);
    TH1F *lead_mu_eta_particle = new TH1F("lead_mu_#eta_particle", "#eta_{mu1}_particle", etaBins, letamin, letamax); listOfTH1.push_back(lead_mu_eta_particle);
    TH1F *lead_mu_phi_particle = new TH1F("lead_mu_#phi_particle", "#phi_{mu1}_particle", phiBins, -TMath::Pi(), TMath::Pi()); listOfTH1.push_back(lead_mu_phi_particle);
    TH1F *sublead_mu_pt_particle = new TH1F("sublead_mu_pT_particle", "p^{T}_{mu2}_particle", pTBins, lpTmin, lpTmax); listOfTH1.push_back(sublead_mu_pt_particle);
    TH1F *sublead_mu_eta_particle = new TH1F("sublead_mu_#eta_particle", "#eta_{mu2}_particle", etaBins, letamin, letamax); listOfTH1.push_back(sublead_mu_eta_particle);
    TH1F *sublead_mu_phi_particle = new TH1F("sublead_mu_#phi_particle", "#phi_{mu2}_particle", phiBins, -TMath::Pi(), TMath::Pi()); listOfTH1.push_back(sublead_mu_phi_particle);
    TH1F *goodMu_size_particle = new TH1F("goodMu_size_particle", "size", 5, 0, 5); listOfTH1.push_back(goodMu_size_particle);

    // TH1F *goodE_size_parton = new TH1F("goodE_size_parton", "size", 5, 0, 5); listOfTH1.push_back(goodE_size_parton);
    // TH1F *goodMu_size_parton = new TH1F("goodMu_size_parton", "size", 5, 0, 5); listOfTH1.push_back(goodMu_size_parton);


// 2D - parton(1) particle(2) reco(3)

  // higgs
    TH2F *hHpT12Comp = new TH2F("H_pT_comp_12", "p_{T}^{hbb} Parton vs. Particle", pTBins, hpTmin, hpTmax, pTBins, hpTmin, hpTmax); listOfTH2.push_back(hHpT12Comp);
    TH2F *hHpT23Comp = new TH2F("H_pT_comp_23", "p_{T}^{hbb} Particle vs. Reco", pTBins, hpTmin, hpTmax, pTBins, hpTmin, hpTmax); listOfTH2.push_back(hHpT23Comp);
    TH2F *hHpT13Comp = new TH2F("H_pT_comp_13", "p_{T}^{hbb} Parton vs. Reco", pTBins, hpTmin, hpTmax, pTBins, hpTmin, hpTmax); listOfTH2.push_back(hHpT13Comp);
    TH2F *hHm12Comp = new TH2F("H_m_comp_12", "m_{hbb} Parton vs. Particle", mBins, hmmin, hmmax, mBins, hmmin, hmmax); listOfTH2.push_back(hHm12Comp);
    TH2F *hHm23Comp = new TH2F("H_m_comp_23", "m_{hbb} Particle vs. Reco", mBins, hmmin, hmmax, mBins, hmmin, hmmax); listOfTH2.push_back(hHm23Comp);      
    TH2F *hHm13Comp = new TH2F("H_m_comp_13", "m_{hbb} Parton vs. Reco", mBins, hmmin, hmmax, mBins, hmmin, hmmax); listOfTH2.push_back(hHm13Comp);
    TH2F*hbbdeltaPhi12Comp = new TH2F("bb_#Delta#phi_comp_12", "#Delta#phi_{bb} Parton vs. Particle", phiBins, -TMath::Pi(),+TMath::Pi(), phiBins, -TMath::Pi(),+TMath::Pi()); listOfTH2.push_back(hbbdeltaPhi12Comp);
    TH2F*hbbdeltaPhi23Comp = new TH2F("bb_#Delta#phi_comp_23", "#Delta#phi_{bb} Particle vs. Reco", phiBins, -TMath::Pi(),+TMath::Pi(), phiBins, -TMath::Pi(),+TMath::Pi()); listOfTH2.push_back(hbbdeltaPhi23Comp);
    TH2F*hbbdeltaPhi13Comp = new TH2F("bb_#Delta#phi_comp_13", "#Delta#phi_{bb} Parton vs. Reco", phiBins, -TMath::Pi(),+TMath::Pi(), phiBins, -TMath::Pi(),+TMath::Pi()); listOfTH2.push_back(hbbdeltaPhi13Comp); 
    TH2F*hbbdeltaEta12Comp = new TH2F("bb_#Delta#eta_comp_12", "#Delta#eta_{bb} Parton vs. Particle", etaBins, hetamin, hetamax, etaBins, hetamin, hetamax); listOfTH2.push_back(hbbdeltaEta12Comp);
    TH2F*hbbdeltaEta23Comp = new TH2F("bb_#Delta#eta_comp_23", "#Delta#eta_{bb} Particle vs. Reco", etaBins, hetamin, hetamax, etaBins, hetamin, hetamax); listOfTH2.push_back(hbbdeltaEta23Comp);
    TH2F*hbbdeltaEta13Comp = new TH2F("bb_#Delta#eta_comp_13", "#Delta#eta_{bb} Parton vs. Reco", etaBins, hetamin, hetamax, etaBins, hetamin, hetamax); listOfTH2.push_back(hbbdeltaEta13Comp);

  // vbfj
    TH2F*hjjpT12Comp = new TH2F("jj_pT_comp_12", "p_{T}^{jj} Parton vs. Particle", pTBins, jpTmin, jpTmax, pTBins, jpTmin, jpTmax); listOfTH2.push_back(hjjpT12Comp);
    TH2F*hjjpT23Comp = new TH2F("jj_pT_comp_23", "p_{T}^{jj} Particle vs. Reco", pTBins, jpTmin, jpTmax, pTBins, jpTmin, jpTmax); listOfTH2.push_back(hjjpT23Comp);
    TH2F*hjjpT13Comp = new TH2F("jj_pT_comp_13", "p_{T}^{jj} Parton vs. Reco", pTBins, jpTmin, jpTmax, pTBins, jpTmin, jpTmax); listOfTH2.push_back(hjjpT13Comp);
    TH2F*hjjdeltaPhi12Comp = new TH2F("jj_#Delta#phi_comp_12", "#Delta#phi_{jj} Parton vs. Particle", phiBins, -TMath::Pi(),+TMath::Pi(), phiBins, -TMath::Pi(),+TMath::Pi()); listOfTH2.push_back(hjjdeltaPhi12Comp);
    TH2F*hjjdeltaPhi23Comp = new TH2F("jj_#Delta#phi_comp_23", "#Delta#phi_{jj} Particle vs. Reco", phiBins, -TMath::Pi(),+TMath::Pi(), phiBins, -TMath::Pi(),+TMath::Pi()); listOfTH2.push_back(hjjdeltaPhi23Comp);
    TH2F*hjjdeltaPhi13Comp = new TH2F("jj_#Delta#phi_comp_13", "#Delta#phi_{jj} Parton vs. Reco", phiBins, -TMath::Pi(),+TMath::Pi(), phiBins, -TMath::Pi(),+TMath::Pi()); listOfTH2.push_back(hjjdeltaPhi13Comp);

  // leptons
    TH2F*hl1pT12Comp = new TH2F("l1_pT_comp_12", "p^{T}_{l1}  Parton vs. Particle", pTBins, lpTmin, lpTmax, pTBins, lpTmin, lpTmax); listOfTH2.push_back(hl1pT12Comp);
    TH2F*hl1pT23Comp = new TH2F("l1_pT_comp_23", "p^{T}_{l1} Particle vs. Reco", pTBins, lpTmin, lpTmax, pTBins, lpTmin, lpTmax); listOfTH2.push_back(hl1pT23Comp);
    TH2F*hl1pT13Comp = new TH2F("l1_pT_comp_13", "p^{T}_{l1} Parton vs. Reco", pTBins, lpTmin, lpTmax, pTBins, lpTmin, lpTmax); listOfTH2.push_back(hl1pT13Comp);
    TH2F*hl2pT12Comp = new TH2F("l2_pT_comp_12", "p^{T}_{l2}  Parton vs. Particle", pTBins, lpTmin, lpTmax, pTBins, lpTmin, lpTmax); listOfTH2.push_back(hl2pT12Comp);
    TH2F*hl2pT23Comp = new TH2F("l2_pT_comp_23", "p^{T}_{l2} Particle vs. Reco", pTBins, lpTmin, lpTmax, pTBins, lpTmin, lpTmax); listOfTH2.push_back(hl2pT23Comp);
    TH2F*hl2pT13Comp = new TH2F("l2_pT_comp_13", "p^{T}_{l2} Parton vs. Reco", pTBins, lpTmin, lpTmax, pTBins, lpTmin, lpTmax); listOfTH2.push_back(hl2pT13Comp);

  // z comps
    TH2F*hz1pT12Comp = new TH2F("z1_pT_comp_12", "p^{T}_{z1}  Parton vs. Particle", pTBins, zpTmin, zpTmax, pTBins, zpTmin, zpTmax); listOfTH2.push_back(hz1pT12Comp);
    TH2F*hz1pT23Comp = new TH2F("z1_pT_comp_23", "p^{T}_{z1} Particle vs. Reco", pTBins, zpTmin, zpTmax, pTBins, zpTmin, zpTmax); listOfTH2.push_back(hz1pT23Comp);
    TH2F*hz1pT13Comp = new TH2F("z1_pT_comp_13", "p^{T}_{z1} Parton vs. Reco", pTBins, zpTmin, zpTmax, pTBins, zpTmin, zpTmax); listOfTH2.push_back(hz1pT13Comp);
    TH2F*hz2pT12Comp = new TH2F("z2_pT_comp_12", "p^{T}_{z2}  Parton vs. Particle", pTBins, zpTmin, zpTmax, pTBins, zpTmin, zpTmax); listOfTH2.push_back(hz2pT12Comp);
    TH2F*hz2pT23Comp = new TH2F("z2_pT_comp_23", "p^{T}_{z2} Particle vs. Reco", pTBins, zpTmin, zpTmax, pTBins, zpTmin, zpTmax); listOfTH2.push_back(hz2pT23Comp);
    TH2F*hz2pT13Comp = new TH2F("z2_pT_comp_13", "p^{T}_{z2} Parton vs. Reco", pTBins, zpTmin, zpTmax, pTBins, zpTmin, zpTmax); listOfTH2.push_back(hz2pT13Comp);
    TH2F*hz1m12Comp = new TH2F("z1_m_comp_12", "m_{z1}  Parton vs. Particle", mBins, zmmin, zmmax, mBins, zmmin, zmmax); listOfTH2.push_back(hz1m12Comp);
    TH2F*hz1m23Comp = new TH2F("z1_m_comp_23", "m_{z1} Particle vs. Reco", mBins, zmmin, zmmax, mBins, zmmin, zmmax); listOfTH2.push_back(hz1m23Comp);      
    TH2F*hz1m13Comp = new TH2F("z1_m_comp_13", "m_{z1} Parton vs. Reco", mBins, zmmin, zmmax, mBins, zmmin, zmmax); listOfTH2.push_back(hz1m13Comp);
    TH2F*hz2m12Comp = new TH2F("z2_m_comp_12", "m_{z2}  Parton vs. Particle", mBins, zmmin, zmmax, mBins, zmmin, zmmax); listOfTH2.push_back(hz2m12Comp);
    TH2F*hz2m23Comp = new TH2F("z2_m_comp_23", "m_{z2} Particle vs. Reco", mBins, zmmin, zmmax, mBins, zmmin, zmmax); listOfTH2.push_back(hz2m23Comp);
    TH2F*hz2m13Comp = new TH2F("z2_m_comp_13", "m_{z2} Parton vs. Reco", mBins, zmmin, zmmax, mBins, zmmin, zmmax); listOfTH2.push_back(hz2m13Comp);

  // w comps
    TH2F*hw1pT12Comp = new TH2F("w1_pT_comp_12", "p^{T}_{w1}  Parton vs. Particle", pTBins, wpTmin, wpTmax, pTBins, wpTmin, wpTmax); listOfTH2.push_back(hw1pT12Comp);
    TH2F*hw1pT23Comp = new TH2F("w1_pT_comp_23", "p^{T}_{w1} Particle vs. Reco", pTBins, wpTmin, wpTmax, pTBins, wpTmin, wpTmax); listOfTH2.push_back(hw1pT23Comp);
    TH2F*hw1pT13Comp = new TH2F("w1_pT_comp_13", "p^{T}_{w1} Parton vs. Reco", pTBins, wpTmin, wpTmax, pTBins, wpTmin, wpTmax); listOfTH2.push_back(hw1pT13Comp);
    TH2F*hw2pT12Comp = new TH2F("w2_pT_comp_12", "p^{T}_{w2}  Parton vs. Particle", pTBins, wpTmin, wpTmax, pTBins, wpTmin, wpTmax); listOfTH2.push_back(hw2pT12Comp);
    TH2F*hw2pT23Comp = new TH2F("w2_pT_comp_23", "p^{T}_{w2} Particle vs. Reco", pTBins, wpTmin, wpTmax, pTBins, wpTmin, wpTmax); listOfTH2.push_back(hw2pT23Comp);
    TH2F*hw2pT13Comp = new TH2F("w2_pT_comp_13", "p^{T}_{w2} Parton vs. Reco", pTBins, wpTmin, wpTmax, pTBins, wpTmin, wpTmax); listOfTH2.push_back(hw2pT13Comp);
    TH2F*hw1m12Comp = new TH2F("w1_m_comp_12", "m_{w1}  Parton vs. Particle", mBins, wmmin, wmmax, mBins, wmmin, wmmax); listOfTH2.push_back(hw1m12Comp);
    TH2F*hw1m23Comp = new TH2F("w1_m_comp_23", "m_{w1} Particle vs. Reco", mBins, wmmin, wmmax, mBins, wmmin, wmmax); listOfTH2.push_back(hw1m23Comp);      
    TH2F*hw1m13Comp = new TH2F("w1_m_comp_13", "m_{w1} Parton vs. Reco", mBins, wmmin, wmmax, mBins, wmmin, wmmax); listOfTH2.push_back(hw1m13Comp);
    TH2F*hw2m12Comp = new TH2F("w2_m_comp_12", "m_{w2}  Parton vs. Particle", mBins, wmmin, wmmax, mBins, wmmin, wmmax); listOfTH2.push_back(hw2m12Comp);
    TH2F*hw2m23Comp = new TH2F("w2_m_comp_23", "m_{w2} Particle vs. Reco", mBins, wmmin, wmmax, mBins, wmmin, wmmax); listOfTH2.push_back(hw2m23Comp);
    TH2F*hw2m13Comp = new TH2F("w2_m_comp_13", "m_{w2} Parton vs. Reco", mBins, wmmin, wmmax, mBins, wmmin, wmmax); listOfTH2.push_back(hw2m13Comp);


// 2d mixed
      
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

  // higgs dphi vs. zz dphi
    TH2F*hbbzzdeltaPhicompreco = new TH2F("bb_zz_#Delta#phi_comp_reco", "bb_zz_#Delta#phi_comp_reco", phiBins, -TMath::Pi(), +TMath::Pi(), phiBins, -TMath::Pi(), +TMath::Pi()); listOfTH2.push_back(hbbzzdeltaPhicompreco);
    TH2F*hbbzzdeltaPhicompparticle = new TH2F("bb_zz_#Delta#phi_comp_particle", "bb_zz_#Delta#phi_comp_particle", phiBins, -TMath::Pi(), +TMath::Pi(), phiBins, -TMath::Pi(), +TMath::Pi()); listOfTH2.push_back(hbbzzdeltaPhicompparticle);
    TH2F*hbbzzdeltaPhicompparton = new TH2F("bb_zz_#Delta#phi_comp_parton", "bb_zz_#Delta#phi_comp_parton", phiBins, -TMath::Pi(), +TMath::Pi(), phiBins, -TMath::Pi(), +TMath::Pi()); listOfTH2.push_back(hbbzzdeltaPhicompparton);

  // higgs deta vs. zz deta  
    TH2F*hbbzzdeltaEtacompreco = new TH2F("bb_zz_delta#eta_comp_reco", "bb_zz_delta#eta_comp_reco", etaBins, hetamin, hetamax, etaBins, zetamin, zetamax); listOfTH2.push_back(hbbzzdeltaEtacompreco);
    TH2F*hbbzzdeltaEtacompparticle = new TH2F("bb_zz_delta#eta_comp_particle", "bb_zz_delta#eta_comp_particle", etaBins, hetamin, hetamax, etaBins, zetamin, zetamax); listOfTH2.push_back(hbbzzdeltaEtacompparticle);
    TH2F*hbbzzdeltaEtacompparton = new TH2F("bb_zz_delta#eta_comp_parton", "bb_zz_delta#eta_comp_parton", etaBins, hetamin, hetamax, etaBins, zetamin, zetamax); listOfTH2.push_back(hbbzzdeltaEtacompparton);


// misc

  // leptonic event types
    TH1F *recoET = new TH1F("reco_event_type", "ET", 5, -1, 4); listOfTH1.push_back(recoET);
    TH1F *particleET = new TH1F("particle_event_type", "ET", 5, -1, 4); listOfTH1.push_back(particleET);
    TH1F *partonET = new TH1F("parton_event_type", "ET", 5, -1, 4); listOfTH1.push_back(partonET);
  
  // w - reco by event type
    TH1F *hllpTET0reco = new TH1F("ll_ET0_pT_reco", "p^{T}_ET0_{ll}_reco", pTBins, wpTmin, wpTmax); listOfTH1.push_back(hllpTET0reco);
    TH1F *hllmET0reco = new TH1F("ll_ET0_m_reco", "m_{ll}_ET0_reco", mBins, wmmin, wmmax); listOfTH1.push_back(hllmET0reco);
    TH1F *hllpTET1reco = new TH1F("ll_ET1_pT_reco", "p^{T}_ET1_{ll}_reco", pTBins, wpTmin, wpTmax); listOfTH1.push_back(hllpTET1reco);
    TH1F *hllmET1reco = new TH1F("ll_ET1_m_reco", "m_{ll}_ET1_reco", mBins, wmmin, wmmax); listOfTH1.push_back(hllmET1reco);
    TH1F *hllpTET2reco = new TH1F("ll_ET2_pT_reco", "p^{T}_ET2_{ll}_reco", pTBins, wpTmin, wpTmax); listOfTH1.push_back(hllpTET2reco);
    TH1F *hllmET2reco = new TH1F("ll_ET2_m_reco", "m_{ll}_ET2_reco", mBins, wmmin, wmmax); listOfTH1.push_back(hllmET2reco);
    TH1F *hllpTET3reco = new TH1F("ll_ET3_pT_reco", "p^{T}_ET3_{ll}_reco", pTBins, wpTmin, wpTmax); listOfTH1.push_back(hllpTET3reco);
    TH1F *hllmET3reco = new TH1F("ll_ET3_m_reco", "m_{ll}_ET3_reco", mBins, wmmin, wmmax); listOfTH1.push_back(hllmET3reco);

  // w - particle by event type
    TH1F *hllpTET0particle = new TH1F("ll_ET0_pT_particle", "p^{T}_ET0_{ll}_particle", pTBins, wpTmin, wpTmax); listOfTH1.push_back(hllpTET0particle);
    TH1F *hllmET0particle = new TH1F("ll_ET0_m_particle", "m_{ll}_ET0_particle", mBins, wmmin, wmmax); listOfTH1.push_back(hllmET0particle);
    TH1F *hllpTET1particle = new TH1F("ll_ET1_pT_particle", "p^{T}_ET1_{ll}_particle", pTBins, wpTmin, wpTmax); listOfTH1.push_back(hllpTET1particle);
    TH1F *hllmET1particle = new TH1F("ll_ET1_m_particle", "m_{ll}_ET1_particle", mBins, wmmin, wmmax); listOfTH1.push_back(hllmET1particle);
    TH1F *hllpTET2particle = new TH1F("ll_ET2_pT_particle", "p^{T}_ET2_{ll}_particle", pTBins, wpTmin, wpTmax); listOfTH1.push_back(hllpTET2particle);
    TH1F *hllmET2particle = new TH1F("ll_ET2_m_particle", "m_{ll}_ET2_particle", mBins, wmmin, wmmax); listOfTH1.push_back(hllmET2particle);
    TH1F *hllpTET3particle = new TH1F("ll_ET3_pT_particle", "p^{T}_ET3_{ll}_particle", pTBins, wpTmin, wpTmax); listOfTH1.push_back(hllpTET3particle);
    TH1F *hllmET3particle = new TH1F("ll_ET3_m_particle", "m_{ll}_ET3_particle", mBins, wmmin, wmmax); listOfTH1.push_back(hllmET3particle);

  // mll [0,15]
    TH1F *hmll_0_15_reco = new TH1F("ll_m_reco_0_15", "m_{ll}_reco_0_15", 75, 0, 15); listOfTH1.push_back(hmll_0_15_reco);
    TH1F *hllm_0_15_particle = new TH1F("ll_m_particle_0_15", "m_{ll}_particle_0_15", 75, 0, 15); listOfTH1.push_back(hllm_0_15_particle);

}
*/

//------------------------------------------------------------------------------------------------------------------------------------------------------------
// printing defs
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



#endif