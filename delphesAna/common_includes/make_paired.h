#ifndef MAKE_PAIRED_H
#define MAKE_PAIRED_H
// #include <iostream>
// #include <unordered_set>
// #include <utility>
// #include "TClonesArray.h"
// #include "Math/LorentzVector.h"
// #include "classes/DelphesClasses.h"
// #include "external/ExRootAnalysis/ExRootTreeReader.h"


#include "classes/DelphesClasses.h"
#include "external/ExRootAnalysis/ExRootTreeReader.h"
#include "external/ExRootAnalysis/ExRootResult.h"
#include <cstdlib>
#include <iostream>
#include <functional>
#include <time.h>
#include <Math/GenVector/PtEtaPhiM4D.h>
#include <Math/Vector4D.h>
#include <Math/GenVector/LorentzVector.h>

//class ExRootTreeReader;
namespace paired

{

const double pi = TMath::Pi();

 double deltaPhi(double phi1, double phi2) { return TVector2::Phi_mpi_pi(phi1 - phi2); }
 
 double deltaR_M(double eta1, double phi1, double eta2, double phi2) {
   double deta = eta1 - eta2;
   double dphi = paired::deltaPhi(phi1, phi2);
   return std::hypot(deta, dphi);
 }
 
double shiftpi(double phi, double shift, double lim) {
  if (shift == 0) return phi;
  if (shift > 0) {
    if (phi < lim) return phi + shift;
  }
  else {
    if (phi > lim) return phi + shift;
  }
  return phi;
}

template <class T1, class T2>
double deltaR(const T1 &a, const T2 &b) {
  return paired::deltaR_M(a->Eta, a->Phi, b->Eta, b->Phi);
}



struct ParticleInfo {
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
std::pair<bool,int> isInBridgedJet(const Jet *jet1, const Jet *jet2, T const *p, float jetR, bool bridge=false) {
  if (std::abs(p->Eta) > 5 || p->PT <= 0) return make_pair(false,-1);
  if (!p) return make_pair(false,-1);

  float eta1 = jet1->Eta;
  float eta2 = jet2->Eta;
  float phi1 = jet1->Phi;
  float phi2 = jet2->Phi;
  float eta0,phi0;
  float num, den, dist, jet1d, jet2d;
  den = paired::deltaR(jet1,jet2);

  jet1d = paired::deltaR(jet1,p);
  jet2d = paired::deltaR(jet2,p);
  if (jet1d < jetR) return make_pair(true,1);
  if (jet2d < jetR) return make_pair(true,2);

  if (!bridge) {
    return make_pair(false,-1);
  }

  //bridging

  if (pow(jet1d,2)+pow(jet2d,2) > pow(den,2)) return make_pair(false,-1); //obtuse angle triangle, point lies outside
  eta0 = p->Eta;
  phi0 = p->Phi;
  num = (eta2-eta1)*paired::deltaPhi(phi1,phi0) - (eta1-eta0)*paired::deltaPhi(phi2,phi1);
  dist = num/den;
  if (dist > jetR) return make_pair(false,-1);
  else return make_pair(true,3);

}

template <typename T>
std::pair<bool,int> isInEllipse(const Jet *jet1, const Jet *jet2, T const *p, float semimajoradd = 1.0) {
  float eta1 = jet1->Eta;
  float eta2 = jet2->Eta;
  float phi1 = jet1->Phi;
  float phi2 = jet2->Phi;
  float djet, semimajor, focus, eta_center, phi_center, eta_f1,phi_f1,eta_f2,phi_f2, eta0,phi0, f1dist, f2dist, distsum; //, phi_m1, phi_m2;
  double semiminor;

  djet = paired::deltaR(jet1,jet2);
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

template <typename T>
std::pair<bool,int> isInJet(const Jet *jet1, const Jet *jet2, T const *p, float jetR, bool bridge=false, bool ellipse = false, float semimajoradd = 1.0) {
  if (ellipse) return paired::isInEllipse(jet1,jet2,p,semimajoradd);
  else return paired::isInBridgedJet(jet1,jet2,p,jetR,bridge);
}

//------------------------------------------------------------------------------
std::pair< std::map<TString, float>, std::map<TString, std::vector<float>> > processBridge(const Jet *jet1, const Jet *jet2, const TClonesArray *branchPF, const TClonesArray *branchParticle, float jetR,  bool bridge=true, bool ellipse = true, float semimajoradd = 1.0) {
  
  std::map<TString, float> _floatVars;
  std::map<TString, std::vector<float>> _arrayVars;
  
  // Reco PF particles 
  std::vector<paired::ParticleInfo> particles;
  ParticleFlowCandidate *p;
  int ncands = branchPF->GetEntriesFast();
  bool isparticle1 = true;
  TLorentzVector jetp4, genjetp4;

  for (Int_t i = 0; i < ncands; ++i) {
    p = (ParticleFlowCandidate *)branchPF->At(i);
    auto retpair = paired::isInJet(jet1,jet2,p,jetR,bridge,ellipse,semimajoradd);
    if (retpair.first) {
      particles.push_back(paired::ParticleInfo(p,retpair.second));
      if (isparticle1) {
        isparticle1 = false;
        jetp4 = p->P4();
      }
      else {
        jetp4 += p->P4();
      }
    //  cout << "(" << jet1->Eta << ","<< jet1->Phi << ") ("<< jet2->Eta << ","<< jet2->Phi << ") ("<< p->Eta << ","<< p->Phi << ") "<< retpair.second << endl;
    }
  }
  _floatVars["jet1_pt"] = jet1->PT;
  _floatVars["jet1_eta"] = jet1->Eta;
  _floatVars["jet1_phi"] = jet1->Phi;
  _floatVars["jet1_energy"] = jet1->P4().Energy();
  _floatVars["jet1_nparticles"] =jet1->Constituents.GetEntriesFast();
  _floatVars["jet2_pt"] = jet2->PT;
  _floatVars["jet2_eta"] = jet2->Eta;
  _floatVars["jet2_phi"] = jet2->Phi;
  _floatVars["jet2_energy"] = jet2->P4().Energy();
  _floatVars["jet2_nparticles"] =jet2->Constituents.GetEntriesFast();

  _floatVars["label_cc"] = 0;
  _floatVars["label_bb"] = 0;
  _floatVars["label_ll"] = 0;
  _floatVars["parton_location"] = 0;
  _floatVars["higgs_pt"] = 0;
  _floatVars["drqq"] = 0;
  _floatVars["z_pt"] = 0;
  _floatVars["gen_flav"] = 0;

  _floatVars["dijet_mass"] = jetp4.M();
  _floatVars["dijet_pt"] = jetp4.Pt();
  _floatVars["dijet_eta"] = jetp4.Eta();
  _floatVars["dijet_phi"] = jetp4.Phi();
  
  std::sort(particles.begin(), particles.end(), [](const auto &a, const auto &b) { return a.pt > b.pt; });

  _floatVars["jet_nparticles"] = particles.size();

  for (const auto &p : particles) {
    if (std::abs(p.pz) > 10000 || std::abs(p.eta) > 5 || p.pt <= 0) continue;
    _arrayVars["part_px"].push_back(p.px);
    _arrayVars["part_py"].push_back(p.py);
    _arrayVars["part_pz"].push_back(p.pz);
    _arrayVars["part_phi"].push_back(p.phi);
    _arrayVars["part_eta"].push_back(p.eta);
    _arrayVars["part_energy"].push_back(p.energy);
    _arrayVars["part_pt"].push_back(p.pt);
    _arrayVars["part_deta1"].push_back((jet1->Eta > 0 ? 1 : -1) * (p.eta - jet1->Eta));
    _arrayVars["part_dphi1"].push_back(paired::deltaPhi(p.phi, jet1->Phi));
    _arrayVars["part_deta"].push_back((jet1->Eta > 0 ? 1 : -1) * (p.eta - jet1->Eta));
    _arrayVars["part_dphi"].push_back(paired::deltaPhi(p.phi, jet1->Phi));
    _arrayVars["part_deta2"].push_back((jet2->Eta > 0 ? 1 : -1) * (p.eta - jet2->Eta));
    _arrayVars["part_dphi2"].push_back(paired::deltaPhi(p.phi, jet2->Phi));
    _arrayVars["part_charge"].push_back(p.charge);
    _arrayVars["part_pid"].push_back(p.pid);
    _arrayVars["part_d0val"].push_back(p.d0);
    _arrayVars["part_d0err"].push_back(p.d0err);
    _arrayVars["part_dzval"].push_back(p.dz);
    _arrayVars["part_dzerr"].push_back(p.dzerr);  
    _arrayVars["part_source"].push_back(p.source);  
  }


  // Gen Particles
  std::vector<const GenParticle *> genParticles_;
  for (Int_t i = 0; i < branchParticle->GetEntriesFast(); ++i) {
    genParticles_.push_back((GenParticle *)branchParticle->At(i));
  }

  // Higgs
  bool isHiggs = false;
  float hpt=0, zpt=0, drqq=0, hmass=0;
  const GenParticle *d1, *d2;
  int nc=0,nb=0,nbrid=0;
  std::pair<bool,int> retpair1,retpair2,thisret;
  for (const auto *thisgp : genParticles_) {
    if( thisgp==nullptr || thisgp->D1<0) continue;
    d1 = genParticles_[thisgp->D1];
    if ((thisgp->PID == 25 && d1->PID != 25) || (thisgp->PID == 23 && (abs(d1->PID) == 4 || abs(d1->PID) == 5))) {
      isHiggs = true;
      hpt = thisgp->PT;
      hmass = thisgp->Mass;
      if( thisgp==nullptr || thisgp->D2<0) continue;
      d2 = genParticles_[thisgp->D2];
      drqq = paired::deltaR(d1,d2);
      _floatVars["gen_flav"] = thisgp->PID*10 + abs(d1->PID);
      break;
    }    
  }

  if (isHiggs) {
    retpair1 = paired::isInJet(jet1,jet2,d1,jetR,bridge,ellipse,semimajoradd);
    retpair2 = paired::isInJet(jet1,jet2,d2,jetR,bridge,ellipse,semimajoradd);
    if (retpair1.first && retpair2.first) {
      if (abs(d1->PID)==5 && abs(d2->PID)==5) _floatVars["label_bb"] = 1;
      if (abs(d1->PID)==4 && abs(d2->PID)==4) _floatVars["label_cc"] = 1;

      if (retpair1.second+retpair2.second==3) _floatVars["parton_location"] = 1; //one parton in each jet
      if ((retpair1.second==3 && retpair2.second<3) || (retpair1.second<3 && retpair2.second==3)) _floatVars["parton_location"] = 2; //one parton in jet and one in bridge
      if (retpair1.second==3 && retpair2.second==3) _floatVars["parton_location"] = 3; //both partons in bridge
    }
  }
  // if (_floatVars["label_bb"] || _floatVars["label_cc"]) cout << _floatVars["label_bb"] << " " << _floatVars["label_cc"] << " " << retpair1.second << " " << retpair2.second << " " << _floatVars["parton_location"] << endl;

  // Not Higgs
  else {
    std::vector<const GenParticle *> partons;
    for (const auto *thisgp : genParticles_) {
        if (thisgp->M1 == 0 && thisgp->PID!=23) {
          partons.push_back(thisgp);
          thisret = paired::isInJet(jet1,jet2,thisgp,jetR,bridge,ellipse,semimajoradd);
          if (thisret.first) {
            if (abs(thisgp->PID) == 5) {
              nb++;
              if (thisret.second==3) nbrid++;
            }
            else if (abs(thisgp->PID) == 4) {
              nc++;
              if (thisret.second==3) nbrid++;
            }
          }
        }
    }
    if (partons.size() >=2) {
        sort(partons.begin(),partons.end(),
              [](const GenParticle * lhs, const GenParticle * rhs)  {
                  return lhs->PT > rhs->PT;
              });
        auto dipart = partons.at(0)->P4() + partons.at(1)->P4();
        hpt = dipart.Pt();
        hmass = dipart.M();
        drqq = paired::deltaR(partons.at(0),partons.at(1));

        if (nb>=2) _floatVars["label_bb"] = 1;
        else if (nc>=2) _floatVars["label_cc"] = 1;

        if (nbrid>=2) _floatVars["parton_location"] = 3;
        else if (nbrid > 0) _floatVars["parton_location"] = 2;
        else _floatVars["parton_location"] = 1;
    }
  }

  // Keep Z information
  for (const auto *thisgp : genParticles_) {
    if( thisgp==nullptr || thisgp->D1<0) continue;
    d1 = genParticles_[thisgp->D1];
    if (thisgp->PID == 23 && d1->PID != 23) {
      zpt = thisgp->PT;
      break;
    }
  }

  // make genJet
  isparticle1 = true;
  for (const auto *thisgp : genParticles_) {
    if (thisgp->Status != 1) continue;
    auto retpair = paired::isInJet(jet1,jet2,thisgp,jetR,bridge,ellipse,semimajoradd);
    if (retpair.first) {
      if (isparticle1) {
        isparticle1 = false;
        genjetp4 = thisgp->P4();
      }
      else {
        genjetp4 += thisgp->P4();
      }
    //  cout << "(" << jet1->Eta << ","<< jet1->Phi << ") ("<< jet2->Eta << ","<< jet2->Phi << ") ("<< p->Eta << ","<< p->Phi << ") "<< retpair.second << endl;
    }
  }
  _floatVars["gendijet_mass"] = genjetp4.M();
  _floatVars["gendijet_pt"] = genjetp4.Pt();
  _floatVars["gendijet_eta"] = genjetp4.Eta();
  _floatVars["gendijet_phi"] = genjetp4.Phi();


  // cout << _floatVars["gen_flav"] << " " << zpt << " " << hpt << endl;
  _floatVars["higgs_pt"] = hpt;
  _floatVars["hmass"] = hmass;
  _floatVars["drqq"] = drqq;
  _floatVars["z_pt"] = zpt;
  if (_floatVars["label_bb"] == 0 && _floatVars["label_cc"] == 0)
    _floatVars["label_ll"] = 1;

  return std::make_pair(_floatVars,_arrayVars);  
}

std::vector<float> getEventInfo(const TClonesArray *branchParticle) {
  bool isHiggs = false;
  float hpt=0, zpt=0, gen_flav=0,dRqq=0, hmass=0;
  const GenParticle *d1, *d2;

  std::vector<const GenParticle *> genParticles_;
  for (Int_t i = 0; i < branchParticle->GetEntriesFast(); ++i) {
    genParticles_.push_back((GenParticle *)branchParticle->At(i));
  }

  for (const auto *thisgp : genParticles_) {
    if(thisgp==nullptr) continue;
    if( thisgp->D1 < 0 ) continue; 
    d1 = genParticles_[thisgp->D1];
    if ((thisgp->PID == 25 && d1->PID != 25) || (thisgp->PID == 23 && (abs(d1->PID) == 4 || abs(d1->PID) == 5))) {
      isHiggs = true;
      hpt = thisgp->PT;
      hmass = thisgp->Mass;
      gen_flav = thisgp->PID*10 + abs(d1->PID);
      if( thisgp==nullptr || thisgp->D2<0) continue;
      d2 = genParticles_[thisgp->D2];
      dRqq = paired::deltaR(d1,d2);
    }    
    if (thisgp->PID == 23 && (abs(d1->PID) == 11 || abs(d1->PID) == 13) ) {
      zpt = thisgp->PT;
    }
    if (hpt>0 && zpt>0) break;
  }


  if (!isHiggs) {
    std::vector<const GenParticle *> partons;
    for (const auto *thisgp : genParticles_) {
      if( thisgp==nullptr) continue;
        if (thisgp->M1 == 0 && thisgp->PID!=23) partons.push_back(thisgp);
    }
    if (partons.size() >=2) {
        sort(partons.begin(),partons.end(),
              [](const GenParticle * lhs, const GenParticle * rhs)  {
                  return lhs->PT > rhs->PT;
              });
        auto dipart = partons.at(0)->P4() + partons.at(1)->P4();
        dRqq = paired::deltaR(partons.at(0),partons.at(1));
        hpt = dipart.Pt();
        hmass = dipart.M();
    }
  }

  std::vector<float> ret{hpt,zpt,gen_flav,dRqq,hmass};
  return ret;

}


 std::pair< std::map<TString, float>, std::map<TString, std::vector<float>>>

   PAIReDjointEvent( TClonesArray *branchParticle = nullptr,
					  TClonesArray *branchPFCand = nullptr,
					  TClonesArray *branchJet = nullptr,
					  float jetR = 0.4,
					  bool forwardjet = false, bool bridge=false,
					  bool ellipse = false, float semimajoradd = 1.0, 
					  bool sigonly=false){
  // define branches
  std::map<TString, float> floatVars;
  floatVars["label_bb"] = 0;
  floatVars["label_cc"] = 0;
  floatVars["label_ll"] = 0;

  floatVars["event"] = 0;
  floatVars["jetno"] = 0;

  floatVars["jet1_index"] = 0;
  floatVars["jet1_pt"] = 0;
  floatVars["jet1_eta"] = 0;
  floatVars["jet1_phi"] = 0;
  floatVars["jet1_energy"] = 0;
  floatVars["jet1_nparticles"] = 0;

  floatVars["jet2_index"] = 0;
  floatVars["jet2_pt"] = 0;
  floatVars["jet2_eta"] = 0;
  floatVars["jet2_phi"] = 0;
  floatVars["jet2_energy"] = 0;
  floatVars["jet2_nparticles"] = 0;
  floatVars["jet_nparticles"] = 0;

  floatVars["higgs_pt"] = 0;
  floatVars["hmass"] = 0;
  floatVars["drqq"] = 0;
  floatVars["z_pt"] = 0;
  floatVars["gen_flav"] = 0;
  floatVars["parton_location"] = 0;  

  floatVars["dijet_mass"] = 0;
  floatVars["dijet_pt"] = 0;
  floatVars["dijet_eta"] = 0;
  floatVars["dijet_phi"] = 0;
  floatVars["gendijet_mass"] = 0;
  floatVars["gendijet_pt"] = 0;
  floatVars["gendijet_eta"] = 0;
  floatVars["gendijet_phi"] = 0;

  floatVars["dijet_mass_nobridge"] = 0;
  floatVars["dijet_pt_nobridge"] = 0;
  floatVars["dijet_eta_nobridge"] = 0;
  floatVars["dijet_phi_nobridge"] = 0;

  floatVars["n_c"] = 0;

  std::map<TString, std::vector<float>> arrayVars;
  arrayVars["part_px"];
  arrayVars["part_py"];
  arrayVars["part_pz"];
  arrayVars["part_phi"];
  arrayVars["part_eta"];
  arrayVars["part_energy"];
  arrayVars["part_pt"];
  arrayVars["part_deta1"];
  arrayVars["part_dphi1"];
  arrayVars["part_deta"];
  arrayVars["part_dphi"];
  arrayVars["part_deta2"];
  arrayVars["part_dphi2"];
  arrayVars["part_charge"];
  arrayVars["part_pid"];
  arrayVars["part_d0val"];
  arrayVars["part_d0err"];
  arrayVars["part_dzval"];
  arrayVars["part_dzerr"];
  arrayVars["part_source"];

  float higgspt,genflav,zpt,dRqq,hmass;
  int num_processed = 0;
  int jetno = 0;

  // analyze
  //TClonesArray *branchParticle = treeReader->UseBranch("Particle");
  //TClonesArray *branchPFCand = treeReader->UseBranch("ParticleFlowCandidate");
  //TClonesArray *branchJet = treeReader->UseBranch(jetBranch);
  

  std::pair<std::map<TString, float>, std::map<TString, std::vector<float>>> jetbrid;

  auto eventVect = paired::getEventInfo(branchParticle);
  higgspt = eventVect.at(0);
  zpt = eventVect.at(1);
  genflav = eventVect.at(2);
  dRqq = eventVect.at(3);
  hmass = eventVect.at(4);


  int n_c = 0;
    for (Int_t i = 0; i < branchJet->GetEntriesFast(); ++i) {
      const Jet *jet = (Jet *)branchJet->At(i);
      if (!forwardjet) {
        if (jet->PT < 20 || std::abs(jet->Eta) > 2.5)  continue;
      }
      if (jet->Flavor==4) n_c++;
    }           
    // Loop over all jets in event
    for (Int_t i = 0; i < branchJet->GetEntriesFast(); ++i) {
      const Jet *jet = (Jet *)branchJet->At(i);
      if (!forwardjet) {
        if (jet->PT < 20 || std::abs(jet->Eta) > 2.5)  continue;
      }
      else {
        if (jet->PT < 20)  continue;
      }

      std::map<TString, float> floats;
      std::map<TString, std::vector<float>> arrays;
      Int_t j;

      for (j = i+1; j < branchJet->GetEntriesFast(); ++j) {
        const Jet *jet2 = (Jet *)branchJet->At(j);
        if (!forwardjet) {
          if (jet2->PT < 20 || std::abs(jet2->Eta) > 2.5)  continue;
        }
        else {
          if (jet2->PT < 20) continue;
          if (std::abs(jet->Eta) > 2.5 && std::abs(jet2->Eta) > 2.5)  continue;
          // if (std::abs(jet->Eta) > 2.5 || std::abs(jet2->Eta) > 2.5)  
          //     cout << "Recovered a forward jet." << endl;
        }

        for (auto &v : floatVars) {
          v.second = 0;
        }
        for (auto &v : arrayVars) {
          v.second.clear();
        }

        jetbrid = paired::processBridge(jet, jet2, branchPFCand, branchParticle, jetR, bridge, ellipse, semimajoradd);
        floats = jetbrid.first;
        arrays = jetbrid.second;

        for (auto &p : arrays) {
          arrayVars[p.first] = p.second;
        }
        for (auto &p : floats) {
          floatVars[p.first] = p.second;
        }
        ///floatVars["event"] = entry;
        floatVars["n_c"] = n_c;

        if (sigonly) {
          if (floatVars["label_ll"]==1) continue;
        }

        auto dijet = jet->P4() + jet2->P4();
        floatVars["dijet_mass_nobridge"] = dijet.M();
        floatVars["dijet_pt_nobridge"] = dijet.Pt();
        floatVars["dijet_eta_nobridge"] = dijet.Eta();
        floatVars["dijet_phi_nobridge"] = dijet.Phi();
        floatVars["jetno"] = jetno;

        ++jetno;
        ++num_processed;
      }
    }

    std::pair< std::map<TString, float>, std::map<TString, std::vector<float>>>  output=make_pair(floatVars,arrayVars);
    
    return output;
}




}
//------------------------------------------------------------------------------
#endif
