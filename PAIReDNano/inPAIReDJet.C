// #include <iostream>
// #include <unordered_set>
// #include <utility>
// #include "TClonesArray.h"
// #include "Math/LorentzVector.h"
// #include "classes/DelphesClasses.h"
// #include "external/ExRootAnalysis/ExRootTreeReader.h"

#ifdef __CLING__
R__LOAD_LIBRARY(libDelphes)
#include "classes/DelphesClasses.h"
#include "external/ExRootAnalysis/ExRootTreeReader.h"
#include "external/ExRootAnalysis/ExRootResult.h"
#include <cstdlib>
#include <iostream>
#include <functional>
#include <time.h>
#else
class ExRootTreeReader;
#endif

const double pi = TMath::Pi();

double deltaPhi(double phi1, double phi2) { return TVector2::Phi_mpi_pi(phi1 - phi2); }

double deltaR(double eta1, double phi1, double eta2, double phi2) {
  double deta = eta1 - eta2;
  double dphi = deltaPhi(phi1, phi2);
  return std::hypot(deta, dphi);
}

double shiftpi(double phi, double shift, double lim) {
  if (shift == 0) return phi;
  if (shift > 0) {
    if (phi < lim) return phi + shift;
  } else if (phi > lim) {
    return phi + shift;
  }
  return phi;
}

template <class T1, class T2>
double deltaR(const T1 &a, const T2 &b) {
  return deltaR(a->Eta, a->Phi, b->Eta, b->Phi);
}



struct ParticleInfo {
  // initialize generated particle
  ParticleInfo(const GenParticle *particle, int src = -1) {
    pt = particle->PT;
    eta = particle->Eta;
    phi = particle->Phi;
    mass = particle->Mass;
    p4 = ROOT::Math::PtEtaPhiMVector(pt, eta, phi, mass);
    px = p4.px();
    py = p4.py();
    pz = p4.pz();
    energy = p4.energy();
    charge = particle->Charge;
    pid = particle->PID;
    fUniqueID = particle->GetUniqueID();
    source = src;
  }

  // initialize reconstructed particle
  ParticleInfo(const ParticleFlowCandidate *particle, int src = -1) {
    pt = particle->PT;
    eta = particle->Eta;
    phi = particle->Phi;
    mass = particle->Mass;
    p4 = ROOT::Math::PtEtaPhiMVector(pt, eta, phi, mass);
    px = p4.px();
    py = p4.py();
    pz = p4.pz();
    energy = p4.energy();
    charge = particle->Charge;
    pid = particle->PID;
    d0 = particle->D0;
    d0err = particle->ErrorD0;
    dz = particle->DZ;
    dzerr = particle->ErrorDZ;
    source = src;
  }

  double pt;
  double eta;
  double phi;
  double mass;
  double px;
  double py;
  double pz;
  double energy;
  ROOT::Math::PtEtaPhiMVector p4;

  int charge;
  int pid;
  int fUniqueID;
  int source;

  float d0 = 0;
  float d0err = 0;
  float dz = 0;
  float dzerr = 0;
};

template <typename T>
std::pair<bool,int> isInEllipse(const Jet *jet1, const Jet *jet2, T const *p, float semimajoradd = 1.0) {
  float eta1 = jet1->Eta;
  float eta2 = jet2->Eta;
  float phi1 = jet1->Phi;
  float phi2 = jet2->Phi;
  float djet, semimajor, focus, eta_center, phi_center, eta_f1,phi_f1,eta_f2,phi_f2, eta0,phi0, f1dist, f2dist, distsum; //, phi_m1, phi_m2;
  double semiminor;

  djet = deltaR(jet1,jet2);
  semiminor = 1.5;
  semimajor = std::max({semiminor,double(djet/2+semimajoradd)});
  focus = pow(pow(semimajor,2)-pow(semiminor,2),0.5);  // Distance of focus to center

  eta_center = (eta1+eta2)/2;
  phi_center = TVector2::Phi_mpi_pi(phi1 + TVector2::Phi_mpi_pi(phi2 - phi1)/2);
 
  //focus 1
  eta_f1 = eta_center + focus/(djet/2)  * (eta1-eta_center);
  phi_f1 = TVector2::Phi_mpi_pi(phi_center + focus/(djet/2)  *TVector2::Phi_mpi_pi(phi1-phi_center));
  
  //focus 2
  eta_f2 = eta_center + focus/(djet/2)  * (eta2-eta_center);
  phi_f2 = TVector2::Phi_mpi_pi(phi_center + focus/(djet/2)  *TVector2::Phi_mpi_pi(phi2-phi_center));

  // Two ends of major axis. This is necesssary to make sure that the point p is not in between the foci on the wrong side of the phi axis
  // phi_m1 = TVector2::Phi_mpi_pi(phi_center + semimajor/(djet/2)  *TVector2::Phi_mpi_pi(phi1-phi_center));
  // phi_m2 = TVector2::Phi_mpi_pi(phi_center + semimajor/(djet/2)  *TVector2::Phi_mpi_pi(phi2-phi_center));

  double shift = 0, lim = 0;
  // if (phi_center > phi_m1 && phi_center > phi_m2) shift = 2*pi;
  // if (phi_center < phi_m1 && phi_center < phi_m2) shift = -2*pi;
  
  if (phi_center >= 0) {
    shift = 2*pi;
    lim = phi_center - pi;
  }
  else {
    shift = -2*pi;
    lim = phi_center + pi;
  }

  // if (abs(phi1-phi2) > 3.4) cout  << "(" << eta1 << "," << phi1 << "), " << "(" << eta2 << "," << phi2 << "), " << "(" << eta_f1 << "," << phi_f1 << "), "<< "(" << eta_f2 << "," << phi_f2 << "), " << "(" << eta_center << "," << phi_center << ")" << endl;

  eta0 = p->Eta;
  phi0 = shiftpi(p->Phi,shift,lim);

  f1dist = std::hypot(eta0-eta_f1,phi0-shiftpi(phi_f1,shift,lim));
  f2dist = std::hypot(eta0-eta_f2,phi0-shiftpi(phi_f2,shift,lim));
  distsum = f1dist+f2dist;

  //if in ellipse, the sum of the distances will be less than 2*semimajor
  if (distsum < 2*semimajor) return make_pair(true,1);
  else return make_pair(false,-1);
}