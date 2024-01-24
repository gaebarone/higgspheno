//#ifdef __CLING__
//R__LOAD_LIBRARY(libDelphes)
#ifndef GHOST_TAGGING_H
#define GHOST_TAGGING_H
#include "TClonesArray.h"
#include "classes/DelphesClasses.h"
//#endif
#include "inEllipse.h"
#include "TRandom.h"
#include "TRandom3.h"
#include <iostream>
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"

// see https://pdg.lbl.gov/2007/reviews/montecarlorpp.pdf
int get_bottom_number(int pid) {
  pid = abs(pid) ;
  //if( pid % 10 == 5 ) return pid; 
  if (pid == 5) {
    return 1;
  } else if (pid / 100 == 55) {
    return 2;
  } else if ((pid / 10) % 100 == 55) {
    return 2;
  } else if (pid / 1000 == 5) {
    return 1;
  } else if ((pid / 100) % 10 == 5) {
    return 1;
  } else {
    return 0;
  }
}





bool ghost_btag(TClonesArray *branchGenParticle, Jet *jet,double jet_radius = 0.4) {

  // loop through particles and check if in jet
  int nEntries =((TClonesArray*)branchGenParticle)->GetEntries();
  for(int i=0; i<nEntries; i++){
    GenParticle *particle=(GenParticle*) branchGenParticle->At(i);
    //std::cout<<"Checking particle "<<i<<" delta R "<<jet->P4().DeltaR(particle->P4())<< " PID "<<particle->PID<<" function "<<get_bottom_number(particle->PID)<<std::endl;
    if (jet->P4().DeltaR(particle->P4()) > jet_radius) continue;
    // check for b hadrons
    // rivet equivalent 
    if (Rivet::PID::hasBottom(particle->PID)) { /*std::cout<< "return b"<<std::endl;*/ return true;}
    //if (get_bottom_number(particle->PID) > 0) { /*std::cout<< "return b"<<std::endl;*/ return true;}
    //if(   isBottomMeson )
  }
  // no b, return false
  return false;
}

bool ghost_bbtag(TClonesArray *branchGenParticle, Jet *jet1, Jet *jet2) {
  int num_bs = 0;
  //int particle_bs = 0;
  int  particle_bs = 0;
  // loop through particles and check if in jet
  for(int i=0; i<(int) ((TClonesArray*)branchGenParticle)->GetEntries(); i++){
    GenParticle *particle=(GenParticle*) (branchGenParticle->At(i));
    if (!inEllipse(jet1, jet2, particle)) continue;
    // check for b hadrons
    //particle_bs = get_bottom_number(particle->PID);
    if (Rivet::PID::hasBottom(particle->PID)) {
      // safely check that daughters are not more bs
      int d1_pid = 9999;
      if (particle->D1 != -1) {
       GenParticle* daughter1 = (GenParticle*) branchGenParticle->At(particle->D1);
        d1_pid = daughter1 -> PID;
      }
      int d2_pid = 9999;
      if (particle->D2 != -1) {
        GenParticle* daughter2 = (GenParticle*) branchGenParticle->At(particle->D2);
        d2_pid = daughter2 -> PID;
      }
      //if (get_bottom_number(d1_pid) == 0 && get_bottom_number(d2_pid) == 0) {
      if( Rivet::PID::hasBottom(d1_pid) && Rivet::PID::hasBottom(d2_pid)){
	// add to total number of bs and return true if hit 2
        //num_bs += particle_bs;
	num_bs++; 
	if (num_bs > 1) return true;
      }
    }
  }
  // < 2 bs, return false
  return false;
}


bool isMyBTag (Jet *jet, TClonesArray *branchGenParticle=nullptr,int seed=0,double jet_radius = 0.4, double effWrk=0.9,double fake_eff=0.01){
  // tracker acceptance
  
  bool isB=ghost_btag(branchGenParticle,jet,jet_radius);
  TRandom3 randEff,randIneff; 
  randEff.SetSeed(seed); 
  randIneff.SetSeed(seed);

  //bool passEff=randEff.Uniform(0,1) < effWrk; 
  //bool passFake=randIneff.Uniform(0,1) < fake_eff;

  bool passEff=randEff.Binomial(1,effWrk) > 0; 
  bool passFake=randIneff.Binomial(1,fake_eff) > 0;
  
  if( isB) return passEff;
  else return passFake; 
}

double ghost_btagPseudoRecoScore(TClonesArray *branchGenParticle, Jet *jet,double jet_radius = 0.4) {
  
  // loop through particles and check if in jet
  int nEntries =((TClonesArray*)branchGenParticle)->GetEntries();
  for(int i=0; i<nEntries; i++){
    GenParticle *particle=(GenParticle*) branchGenParticle->At(i);
    if (jet->P4().DeltaR(particle->P4()) > jet_radius) continue;
    
    if (Rivet::PID::hasBottom(particle->PID)) {
      double score = (1 - jet->P4().DeltaR(particle->P4())/jet_radius);
    }
  }
  // no b, return false
  return 0.0; 
}

double ghost_btagPseudoRecoScoreSmeared(TClonesArray *branchGenParticle, Jet *jet,int seed=0,double jet_radius = 0.4, double effWrk=0.9,double fake_eff=0.01, double smear=0.01){
  TRandom3 randEff,randIneff,smearRand;
  randEff.SetSeed(seed); 
  randIneff.SetSeed(seed);
  smearRand.SetSeed(seed);
  double true_score=ghost_btagPseudoRecoScore(branchGenParticle,jet,jet_radius);
  double smeared_score=smearRand.Gaus(true_score,true_score*smear);
  return smeared_score;
}



#endif
