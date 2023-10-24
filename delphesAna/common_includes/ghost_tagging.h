//#ifdef __CLING__
//R__LOAD_LIBRARY(libDelphes)
#include "TClonesArray.h"
#include "classes/DelphesClasses.h"
//#endif
#include "inEllipse.h"
// see https://pdg.lbl.gov/2007/reviews/montecarlorpp.pdf
int get_bottom_number(int pid) {
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

bool ghost_btag(TClonesArray *branchGenParticle, Jet *jet) {
  double jet_radius = 0.4;
  // loop through particles and check if in jet
  int nEntries =((TClonesArray*)branchGenParticle)->GetEntries();
  for(int i=0; i<nEntries; i++){
    GenParticle *particle=(GenParticle*) branchGenParticle->At(i);
    if (jet->P4().DeltaR(particle->P4()) > jet_radius) continue;
    // check for b hadrons
    if (get_bottom_number(particle->PID) > 0) return true;
  }
  // no b, return false
  return false;
}

bool ghost_bbtag(TClonesArray *branchGenParticle, Jet *jet1, Jet *jet2) {
  int num_bs = 0;
  int particle_bs = 0;
  // loop through particles and check if in jet
  for(int i=0; i<(int) ((TClonesArray*)branchGenParticle)->GetEntries(); i++){
    GenParticle *particle=(GenParticle*) (branchGenParticle->At(i));
    if (!inEllipse(jet1, jet2, particle)) continue;
    // check for b hadrons
    particle_bs = get_bottom_number(particle->PID);
    if (particle_bs > 0) {
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
      if (get_bottom_number(d1_pid) == 0 && get_bottom_number(d2_pid) == 0) {
        // add to total number of bs and return true if hit 2
        num_bs += particle_bs;
        if (num_bs > 1) return true;
      }
    }
  }
  // < 2 bs, return false
  return false;
}
