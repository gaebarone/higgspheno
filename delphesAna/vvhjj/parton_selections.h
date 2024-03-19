#ifndef PARTON_SELECTIONS_H
#define PARTON_SELECTIONS_H

#include "HepMC/GenParticle.h"
#include "classes/DelphesClasses.h"
#include "classes/DelphesLHEFReader.h"
#include "external/ExRootAnalysis/ExRootTreeReader.h"
#include "../common_includes/ghost_tagging.h"
#include "../common_includes/combinations.h"
#include "../common_includes/get_cross_section.h"
#include <iostream>
#include <fstream>
#include <string>
#include <map>
#include "TH1F.h"
#include "TH2F.h"
#include "TClonesArray.h"
#include "TTree.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <TROOT.h>
#include <TLorentzVector.h>
#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>
#include <TMath.h>
#include <Rtypes.h>
#include <TString.h>
#include <TRandom.h>
#include <TRandom3.h>
#include "TParticle.h"
#include "TCanvas.h"
#include "TProfile.h"
#include <vector>
#include <iomanip>

GenParticle* find_parent(TClonesArray *branchGenParticle, GenParticle *particle, int target) {
  
  // check particle itself
  if (particle -> PID != target){
    //cout<<endl;
    //if( particle != nullptr) cout<<"Done first Parent is "<<particle->PID<<endl;
    return particle;
  }
  
  
  // safely access daughters
  int d1_pid = 9999;
  int d2_pid = 9999;
  GenParticle *parent1;
  GenParticle *parent2;
  if (particle->M1 != -1) {
    parent1 = (GenParticle*) branchGenParticle->At(particle->M1);
    d1_pid = parent1 -> PID;
  }
  if (particle->M2 != -1) {
    parent2 = (GenParticle*) branchGenParticle->At(particle->M2);
    d2_pid = parent2 -> PID;
  }
  //cout<<"Particle "<<particle->PID<<" parent 1 "<<parent1->PID<<" 2 "<<parent2->PID<<endl;

  if (d1_pid == target ) {
    //cout<<"-> same particle  for parent 1 "<<endl;
    return find_parent(branchGenParticle, parent1,target);
  }
  else
    return parent1;
  /*
  // recursive call on parents or return null
  if ( d1_pid != target && d2_pid==target){
  //cout<<"Done Parent 1 is "<<parent1->PID<<endl;
  return parent1; 
  }
  else if( d2_pid !=target && d1_pid==target){
  //cout<<"Done Parent 2 is "<<parent2->PID<<endl;
  return parent2;
  }
  else if (d1_pid == target ) {
  //cout<<"-> same particle  for parent 1 "<<endl;
  return find_parent(branchGenParticle, parent1,target);
  }
  else if (d2_pid == target) {
  //cout<<"-> same particle  for parent 2 "<<endl;
  return find_parent(branchGenParticle, parent2, target);
  }
  else{
  //cout<<"Done Parent is "<<particle->PID<<endl;
  //cout<<endl;
  return particle;
  }
  */
}



//------------------------------------------------------------------------------------------------------------------------------------------------------------
// define find_status1_child
//------------------------------------------------------------------------------------------------------------------------------------------------------------

GenParticle* find_status1_child(TClonesArray *branchGenParticle, GenParticle *particle, int target) {
  // check particle itself
  if (abs(particle -> PID) == target && particle -> Status == 1) {
    return particle;
  }
  // safely access daughters
  int d1_pid = 9999;
  int d2_pid = 9999;
  GenParticle *daughter1;
  GenParticle *daughter2;
  if (particle->D1 != -1) {
    daughter1 = (GenParticle*) branchGenParticle->At(particle->D1);
    d1_pid = daughter1 -> PID;
  }
  if (particle->D2 != -1) {
    daughter2 = (GenParticle*) branchGenParticle->At(particle->D2);
    d2_pid = daughter2 -> PID;
  }
  // recursive call on daughters or return null
  if (abs(d1_pid) == target) {
    return find_status1_child(branchGenParticle, daughter1, target);
  } else if (abs(d2_pid) == target) {
    return find_status1_child(branchGenParticle, daughter2, target);
  } else {
    return nullptr;
  }
}


bool FillHiggsTruthRecord(TClonesArray *branchGenParticle, TLorentzVector & h_parton, TLorentzVector & b1_parton, TLorentzVector & b2_parton, TLorentzVector & j1_parton, TLorentzVector &j2_parton){
    
  bool foundHiggsOrinTree=false;
  bool foundHiggsToBbb=false;
  bool foundHiggsToBbbJets=false;

  GenParticle *daughter1=nullptr;
  GenParticle *daughter2=nullptr;

  for(int i=0; i<(int)branchGenParticle->GetEntries(); i++) {

    GenParticle *particle=(GenParticle*) branchGenParticle->At(i);

    int d1_pid = 9999;
    int d2_pid = 9999;

    if (particle->D1 != -1) {
      daughter1 = (GenParticle*) branchGenParticle->At(particle->D1);
      d1_pid = daughter1 -> PID;
    }

    if (particle->D2 != -1) {
      daughter2 = (GenParticle*) branchGenParticle->At(particle->D2);
      d2_pid = daughter2 -> PID;
    } 

    if (particle->PID == 25 && d1_pid != 25 && d2_pid != 25) {

      h_parton = particle->P4();
      foundHiggsOrinTree=true;
    
      if (abs(d1_pid) == 5 && abs(d2_pid) == 5) {

        foundHiggsToBbb=true; 

        if (daughter1 -> PT > daughter2 -> PT) {
          b1_parton = daughter1 -> P4();
          b2_parton = daughter2 -> P4();
        } else {
          b1_parton = daughter2 -> P4();
          b2_parton = daughter1 -> P4();
        }
      }

    } else if ((particle->PID == 1 || particle->PID == 2 || particle->PID == 3 || particle->PID == 4 || particle->PID == 6) && d1_pid != particle->PID && d2_pid != particle->PID) {

      foundHiggsToBbbJets=true;

      if (daughter1 -> PT > daughter2 -> PT) {
	      j1_parton = daughter1 -> P4();
	      j2_parton = daughter2 -> P4();
      } else {
	      j1_parton = daughter2 -> P4();
	      j2_parton = daughter1 -> P4();
      }
    }
  }

  return (foundHiggsOrinTree && foundHiggsToBbb && foundHiggsToBbbJets); 

}





#endif
