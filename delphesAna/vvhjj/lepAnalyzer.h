#ifndef LEPANALYZER_H
#define LEPANALYZER_H

#include "HepMC/GenParticle.h"
#include "classes/DelphesClasses.h"
#include "classes/DelphesLHEFReader.h"
#include "external/ExRootAnalysis/ExRootTreeReader.h"
// #include "../../common_includes/ghost_tagging.h"
#include "../common_includes/combinations.h"
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

using namespace std;


//------------------------------------------------------------------------------------------------------------------------------------------------------------
// Z + W
//------------------------------------------------------------------------------------------------------------------------------------------------------------


void remove_overlaps(vector< pair<int,int>> muPairIndices){
  for( vector< pair<int,int>>::iterator it=muPairIndices.begin(); it!=muPairIndices.end(); it++){
    pair<int,int> one=(*it);
    for(vector< pair<int,int>>::iterator it2=it+1; it2!=muPairIndices.end(); it2++){
      pair<int,int> two=(*it2);
      if( one.first == two.first || one.first == two.second || one.second == two.second || one.second == two.second) muPairIndices.erase(it2--);
    }
  }
}


vector <int> GoodElectronRecoIndices(TClonesArray *branchElectron=nullptr, string analysis="HZZJJ"){

    vector <int> goodE_reco_indices;

    if( branchElectron == nullptr) return goodE_reco_indices;

    for(int i=0; i<(int)branchElectron->GetEntries(); i++){

      Electron *el_reco = (Electron *) branchElectron->At(i);

      // pT and eta cuts 
      if( analysis == "HZZJJ" || analysis == "ZZJJ" ){

        if( el_reco->PT > 5.0 && fabs(el_reco->Eta) < 2.5) goodE_reco_indices.push_back(i);

      } else if( analysis == "HWWJJ" || analysis == "WWJJ" ){

        if( el_reco->PT > 15.0 && fabs(el_reco->Eta) < 2.5) goodE_reco_indices.push_back(i);

      }
    }

    // sort the indices by pT ;
    sort(goodE_reco_indices.begin(), goodE_reco_indices.end(), [branchElectron](const int& lhs, const int& rhs) {
      return ((Electron*)branchElectron->At(lhs))->PT > ((Electron*)branchElectron->At(rhs))->PT;
    });
  
    return goodE_reco_indices; 
  
}

vector <int> GoodElectronParticleIndices(TClonesArray *branchGenParticle=nullptr, string analysis="HZZJJ") {

    vector <int> goodE_particle_indices; 

    for(int i=0; i<(int)branchGenParticle->GetEntries(); i++){

      GenParticle *particle=(GenParticle*) branchGenParticle->At(i); 

      if( particle->Status !=1 ) continue;
      if( fabs(particle->PID) != 11) continue;

      if( analysis == "HZZJJ" || analysis == "ZZJJ "){

        if( particle->PT > 5 && fabs(particle->Eta) < 2.5 ) goodE_particle_indices.push_back(i);

      } else if( analysis == "HWWJJ" || analysis == "WWJJ" ) {
        
        if( particle->PT > 15 && fabs(particle->Eta) < 2.5 ) goodE_particle_indices.push_back(i);

      }  

      //GenParticle *parent=find_parent(branchGenParticle,particle,particle->PID);
      //if(parent == nullptr) continue; 
      //cout<<" --> Parent "<<parent->PID<<endl;
      //if( abs(parent->PID) != 22 &&  abs(parent->PID)!=23  &&  abs(parent->PID)!= 25 ) continue; 
       
      // sort the indices by pT ;
      sort(goodE_particle_indices.begin(), goodE_particle_indices.end(), [branchGenParticle](const int& lhs, const int& rhs) {
	      return ((GenParticle*) branchGenParticle->At(lhs))->PT > ((GenParticle*) branchGenParticle->At(rhs))->PT;
      });
    }

    return goodE_particle_indices;

}


vector <int> GoodMuonRecoIndices(TClonesArray *branchMuon=nullptr, string analysis="HZZJJ"){

    vector <int> goodMu_reco_indices;

    if( branchMuon==nullptr) return goodMu_reco_indices;

    for(int i=0; i<(int)branchMuon->GetEntries(); i++){

      Muon *mu_reco = (Muon *) branchMuon->At(i);

      // pT and eta cuts 
      if( analysis == "HZZJJ" || analysis == "ZZJJ" ){

        if( mu_reco->PT > 5.0 && fabs(mu_reco->Eta) < 2.5) goodMu_reco_indices.push_back(i);

      } else if( analysis == "HWWJJ" || analysis == "WWJJ" ){

        if( mu_reco->PT > 15.0 && fabs(mu_reco->Eta) < 2.5) goodMu_reco_indices.push_back(i);

      }
    }

    // sort the indices by pT ;
    sort(goodMu_reco_indices.begin(), goodMu_reco_indices.end(), [branchMuon](const int& lhs, const int& rhs) {
      return ((Muon*)branchMuon->At(lhs))->PT > ((Muon*)branchMuon->At(rhs))->PT;
    });
  
    return goodMu_reco_indices; 
  
}

vector <int> GoodMuonParticleIndices(TClonesArray *branchGenParticle=nullptr, string analysis="HZZJJ"){

    vector <int> goodMu_particle_indices; 

    for(int i=0; i<(int)branchGenParticle->GetEntries(); i++){

      GenParticle *particle=(GenParticle*) branchGenParticle->At(i); 

      if( particle->Status !=1 ) continue;
      if( fabs(particle->PID) != 13) continue;

      if( analysis == "HZZJJ" || analysis == "ZZJJ "){

        if( particle->PT > 5 && fabs(particle->Eta) < 2.5 ) goodMu_particle_indices.push_back(i);

      } else if( analysis == "HWWJJ" || analysis == "WWJJ" ) {
        
        if( particle->PT > 15 && fabs(particle->Eta) < 2.5 ) goodMu_particle_indices.push_back(i);

      }  

      //GenParticle *parent=find_parent(branchGenParticle,particle,particle->PID);
      //if(parent == nullptr) continue; 
      //cout<<" --> Parent "<<parent->PID<<endl;
      //if( abs(parent->PID) != 22 &&  abs(parent->PID)!=23  &&  abs(parent->PID)!= 25 ) continue; 
       
      // sort the indices by pT ;
      sort(goodMu_particle_indices.begin(), goodMu_particle_indices.end(), [branchGenParticle](const int& lhs, const int& rhs) {
	      return ((GenParticle*) branchGenParticle->At(lhs))->PT > ((GenParticle*) branchGenParticle->At(rhs))->PT;
      });
    }

    return goodMu_particle_indices;
    
}



//------------------------------------------------------------------------------------------------------------------------------------------------------------
// Z 
//------------------------------------------------------------------------------------------------------------------------------------------------------------


vector< pair<int,int>> GetelecRecoPairIndices(TClonesArray *branchElectron=nullptr, vector<int> goodE_reco_indices=vector<int>(0)){
  vector< pair<int,int>> elecRecoPairIndices;
  vector< pair<int,int>> elecRecoPairIndicesIn;
  vector <vector<int>> elecRecoPairIndices_;
  if( goodE_reco_indices.size() > 1 )
    elecRecoPairIndices_=combinationsNoRepetitionAndOrderDoesNotMatter(2,goodE_reco_indices);

    // remove all combinations not satisfying criteria 
    for(int i=0;i<(int)elecRecoPairIndices_.size(); i++){
      int elecRecoIndex=elecRecoPairIndices_[i].at(0);
      int elecRecoIndex2=elecRecoPairIndices_[i].at(1);
      
      Electron *el1_reco=(Electron*) branchElectron->At(elecRecoIndex);
      Electron *el2_reco=(Electron*) branchElectron->At(elecRecoIndex2);

      if( el1_reco->Charge == el2_reco->Charge ) continue;
      
      elecRecoPairIndicesIn.push_back(make_pair(elecRecoIndex,elecRecoIndex2));
    }

    elecRecoPairIndices=elecRecoPairIndicesIn;
    
    sort(elecRecoPairIndices.begin(), elecRecoPairIndices.end(), [branchElectron](const pair<int,int> lhs, const pair<int,int> rhs) {
    	return fabs(((((Electron*)branchElectron->At(lhs.first))->P4() + ((Electron*)branchElectron->At(lhs.second))->P4())).M() -91 ) <
	  fabs( ((((Electron*)branchElectron->At(rhs.first))->P4() + ((Electron*)branchElectron->At(rhs.second))->P4()).M()) -91 ) ; 
      });

    remove_overlaps(elecRecoPairIndices);
 
    return elecRecoPairIndices;
  }

vector< pair<int,int>> GetelecParticlePairIndices(TClonesArray *branchGenParticle=nullptr, vector<int> goodE_particle_indices=vector<int>(0)){
  vector< pair<int,int>> elecParticlePairIndices;
  vector< pair<int,int>> elecParticlePairIndicesIn;
  vector <vector<int>> elecParticlePairIndices_;
  if( goodE_particle_indices.size() > 1 )
    elecParticlePairIndices_=combinationsNoRepetitionAndOrderDoesNotMatter(2,goodE_particle_indices);

    // remove all combinations not satisfying criteria 
    for(int i=0;i<(int)elecParticlePairIndices_.size(); i++){
      int elecParticleIndex=elecParticlePairIndices_[i].at(0);
      int elecParticleIndex2=elecParticlePairIndices_[i].at(1);
      
      GenParticle *el1_particle=(GenParticle*) branchGenParticle->At(elecParticleIndex);
      GenParticle *el2_particle=(GenParticle*) branchGenParticle->At(elecParticleIndex2);

      if( el1_particle->Charge == el2_particle->Charge ) continue;
      
      elecParticlePairIndicesIn.push_back(make_pair(elecParticleIndex,elecParticleIndex2));
    }

    elecParticlePairIndices=elecParticlePairIndicesIn;
    
    sort(elecParticlePairIndices.begin(), elecParticlePairIndices.end(), [branchGenParticle](const pair<int,int> lhs, const pair<int,int> rhs) {
    	return fabs(((((GenParticle*)branchGenParticle->At(lhs.first))->P4() + ((GenParticle*)branchGenParticle->At(lhs.second))->P4())).M() -91 ) <
	  fabs( ((((GenParticle*)branchGenParticle->At(rhs.first))->P4() + ((GenParticle*)branchGenParticle->At(rhs.second))->P4()).M()) -91 ) ; 
      });

    remove_overlaps(elecParticlePairIndices);
 
    return elecParticlePairIndices;
}


vector< pair<int,int>> GetmuRecoPairIndices(TClonesArray *branchMuon=nullptr, vector<int> goodMu_reco_indices=vector<int>(0)){
    vector< pair<int,int>> muRecoPairIndices;
    vector< pair<int,int>> muRecoPairIndicesIn;
    vector <vector<int>> muRecoPairIndices_;
    if(goodMu_reco_indices.size() > 1) 
      muRecoPairIndices_=combinationsNoRepetitionAndOrderDoesNotMatter(2,goodMu_reco_indices);
    
      // remove all combinations not satifing criteria 
    for(int i=0; i<(int)muRecoPairIndices_.size(); i++){
      int muRecoIndex=muRecoPairIndices_[i].at(0);
      int muRecoIndex2=muRecoPairIndices_[i].at(1);
      
      Muon *el1_reco= (Muon*) branchMuon->At(muRecoIndex);
      Muon *el2_reco=(Muon*) branchMuon->At(muRecoIndex2);
      if( el1_reco->Charge ==  el2_reco->Charge ) continue;
      muRecoPairIndicesIn.push_back(make_pair(muRecoIndex,muRecoIndex2));
    }

    muRecoPairIndices=muRecoPairIndicesIn;

    sort(muRecoPairIndices.begin(), muRecoPairIndices.end(), [branchMuon](const pair<int,int> lhs, const pair<int,int> rhs) {
    	return fabs(((((Muon*)branchMuon->At(lhs.first))->P4() + ((Muon*)branchMuon->At(lhs.second))->P4())).M() -91 ) <
	  fabs( ((((Muon*)branchMuon->At(rhs.first))->P4() + ((Muon*)branchMuon->At(rhs.second))->P4()).M()) -91 ) ; 
      });
    
    //sort(muRecoPairIndices.begin(),muRecoPairIndices.end(), [branchMuon]( pair<int,int>   & lhs,  pair<int,int>   & rhs) {
    //	    int index11_reco=(lhs).first;
    //	    int index12_reco=(lhs).second;
    //	    int index21_reco=(rhs).first;
    //	int index22_reco=(rhs).second;
    //	    return fabs(((((Muon*)branchMuon->At(index11_reco))->P4() + ((Muon*)branchMuon->At(index12_reco))->P4())).M() - 91) <
    //	  fabs( ((((Muon*)branchMuon->At(index21_reco))->P4() + ((Muon*)branchMuon->At(index22_reco))->P4()).M()) -91);
    //});
    
    remove_overlaps(muRecoPairIndices);

    return muRecoPairIndices;
}

vector< pair<int,int>> GetmuParticlePairIndices(TClonesArray *branchGenParticle=nullptr, vector<int> goodMu_particle_indices=vector<int>(0)){
  vector< pair<int,int>> muParticlePairIndices;
  vector< pair<int,int>> muParticlePairIndicesIn;
  vector <vector<int>> muParticlePairIndices_;
  if( goodMu_particle_indices.size() > 1 )
    muParticlePairIndices_=combinationsNoRepetitionAndOrderDoesNotMatter(2,goodMu_particle_indices);

    // remove all combinations not satisfying criteria 
    for(int i=0;i<(int)muParticlePairIndices_.size(); i++){
      int muParticleIndex=muParticlePairIndices_[i].at(0);
      int muParticleIndex2=muParticlePairIndices_[i].at(1);
      
      GenParticle *el1_particle=(GenParticle*) branchGenParticle->At(muParticleIndex);
      GenParticle *el2_particle=(GenParticle*) branchGenParticle->At(muParticleIndex2);

      if( el1_particle->Charge == el2_particle->Charge ) continue;
      
      muParticlePairIndicesIn.push_back(make_pair(muParticleIndex,muParticleIndex2));
    }

    muParticlePairIndices=muParticlePairIndicesIn;
    
    sort(muParticlePairIndices.begin(), muParticlePairIndices.end(), [branchGenParticle](const pair<int,int> lhs, const pair<int,int> rhs) {
    	return fabs(((((GenParticle*)branchGenParticle->At(lhs.first))->P4() + ((GenParticle*)branchGenParticle->At(lhs.second))->P4())).M() -91 ) <
	  fabs( ((((GenParticle*)branchGenParticle->At(rhs.first))->P4() + ((GenParticle*)branchGenParticle->At(rhs.second))->P4()).M()) -91 ) ; 
      });

    remove_overlaps(muParticlePairIndices);
 
    return muParticlePairIndices;
}

/*
vector< pair<int,int>> GetelecmuRecoPairIndices(TClonesArray *branchMuon=nullptr, vector<int> goodMu_reco_indices=vector<int>(0)){
    vector< pair<int,int>> elecmuRecoPairIndices;
    vector< pair<int,int>> elecmuRecoPairIndicesIn;
    vector <vector<int>> elecmuRecoPairIndices_;
    if(goodE_reco_indices.size() > 0 && goodMu_reco_indices.size() > 0) 
      elecmuRecoPairIndices_=combinationsNoRepetitionAndOrderDoesNotMatter(2,goodMu_reco_indices);
    
      // remove all combinations not satifing criteria 
    for(int i=0; i<(int)muRecoPairIndices_.size(); i++){
      int elecRecoIndex=muRecoPairIndices_[i].at(0);
      int elecRecoIndex2=muRecoPairIndices_[i].at(1);
      
      Muon *el1_reco= (Muon*) branchMuon->At(elecRecoIndex);
      Muon *el2_reco=(Muon*) branchMuon->At(elecRecoIndex2);
      if( el1_reco->Charge ==  el2_reco->Charge ) continue;
      muRecoPairIndicesIn.push_back(make_pair(elecRecoIndex,elecRecoIndex2));
    }

    muRecoPairIndices=muRecoPairIndicesIn;

    sort(muRecoPairIndices.begin(), muRecoPairIndices.end(), [branchMuon](const pair<int,int> lhs, const pair<int,int> rhs) {
    	return fabs(((((Muon*)branchMuon->At(lhs.first))->P4() + ((Muon*)branchMuon->At(lhs.second))->P4())).M() -91 ) <
	  fabs( ((((Muon*)branchMuon->At(rhs.first))->P4() + ((Muon*)branchMuon->At(rhs.second))->P4()).M()) -91 ) ; 
      });
    
    //sort(muRecoPairIndices.begin(),muRecoPairIndices.end(), [branchMuon]( pair<int,int>   & lhs,  pair<int,int>   & rhs) {
    //	    int index11_reco=(lhs).first;
    //	    int index12_reco=(lhs).second;
    //	    int index21_reco=(rhs).first;
    //	int index22_reco=(rhs).second;
    //	    return fabs(((((Muon*)branchMuon->At(index11_reco))->P4() + ((Muon*)branchMuon->At(index12_reco))->P4())).M() - 91) <
    //	  fabs( ((((Muon*)branchMuon->At(index21_reco))->P4() + ((Muon*)branchMuon->At(index22_reco))->P4()).M()) -91);
    //});
    
    remove_overlaps(muRecoPairIndices);

    return muRecoPairIndices;
}
*/

vector<pair<int,pair<int,int>>> GetRecoPairIndices(vector< pair<int,int>> elecRecoPairIndices, vector< pair<int,int>> muRecoPairIndices, TClonesArray *branchElectron=nullptr,TClonesArray *branchMuon=nullptr){

  vector<pair<int,pair<int,int>>> RecoPairIndices;


for(int i=0; i<(int) elecRecoPairIndices.size(); i++)
      RecoPairIndices.push_back(make_pair(0,elecRecoPairIndices.at(i)));
    for(int i=0; i<(int) muRecoPairIndices.size(); i++)
      RecoPairIndices.push_back(make_pair(1,muRecoPairIndices.at(i)));

    // sort all of the indices by closeness to mZ
    sort(RecoPairIndices.begin(), RecoPairIndices.end(), [branchMuon,branchElectron]  ( const pair<int,pair<int,int>>  lhs , const pair<int,pair<int,int>>   rhs ){

	double mass1_reco=0;
  double mass2_reco=0;
	
	if( lhs.first==0 ) mass1_reco = (((Electron*)branchElectron->At(lhs.second.first))->P4() + ((Electron*)branchElectron->At(lhs.second.second))->P4()).M();
	else mass1_reco = (((Muon*)branchMuon->At(lhs.second.first))->P4() + ((Muon*)branchMuon->At(lhs.second.second))->P4()).M();
	
	if( rhs.first==0 ) mass2_reco = (((Electron*)branchElectron->At(rhs.second.first))->P4() + ((Electron*)branchElectron->At(rhs.second.second))->P4()).M();
	else mass2_reco = (((Muon*)branchMuon->At(rhs.second.first))->P4() + ((Muon*)branchMuon->At(rhs.second.second))->P4()).M();
	
	return fabs(mass1_reco -91.0 )  <  fabs(mass2_reco -91.0); 
	
      });
    
    return RecoPairIndices;
}

vector<pair<int,pair<int,int>>> GetParticlePairIndices(vector< pair<int,int>> elecParticlePairIndices, vector< pair<int,int>> muParticlePairIndices, TClonesArray *branchGenParticle=nullptr){

  vector<pair<int,pair<int,int>>> ParticlePairIndices;


  for(int i=0; i<(int) elecParticlePairIndices.size(); i++)
      ParticlePairIndices.push_back(make_pair(0,elecParticlePairIndices.at(i)));
  for(int i=0; i<(int) muParticlePairIndices.size(); i++)
      ParticlePairIndices.push_back(make_pair(1,muParticlePairIndices.at(i)));

  // sort all of the indices by closeness to mZ
  sort(ParticlePairIndices.begin(), ParticlePairIndices.end(), [branchGenParticle]  ( const pair<int,pair<int,int>>  lhs , const pair<int,pair<int,int>>   rhs ){

    //double mass1_particle=0 ;
    //double mass2_particle=0;
    
    //if( lhs.first==0 ) mass1_particle = (((GenParticle*)branchGenParticle->At(lhs.second.first))->P4() + ((GenParticle*)branchGenParticle->At(lhs.second.second))->P4()).M();
    //else mass1_particle = (((GenParticle*)branchGenParticle->At(lhs.second.first))->P4() + ((GenParticle*)branchGenParticle->At(lhs.second.second))->P4()).M();
    
    //if( rhs.first==0 ) mass2_particle = (((GenParticle*)branchGenParticle->At(rhs.second.first))->P4() + ((GenParticle*)branchGenParticle->At(rhs.second.second))->P4()).M();
    //else mass2_particle = (((GenParticle*)branchGenParticle->At(rhs.second.first))->P4() + ((GenParticle*)branchGenParticle->At(rhs.second.second))->P4()).M();
    
    double mass1_particle = (((GenParticle*)branchGenParticle->At(lhs.second.first))->P4() + ((GenParticle*)branchGenParticle->At(lhs.second.second))->P4()).M();
    double mass2_particle = (((GenParticle*)branchGenParticle->At(rhs.second.first))->P4() + ((GenParticle*)branchGenParticle->At(rhs.second.second))->P4()).M();

    return fabs(mass1_particle -91.0 )  <  fabs(mass2_particle -91.0); 
	
  });

    return ParticlePairIndices;
}

vector <int> GetZPartonIndices(TClonesArray *branchGenParticle=nullptr, string analysis="HZZJJ") {

    vector <int> ZPartonIndices; 

    for(int i=0; i<(int)branchGenParticle->GetEntries(); i++){

      GenParticle *particle=(GenParticle*)branchGenParticle->At(i);

      int d1_pid = 9999;
      int d2_pid = 9999;

      if (particle->D1 != -1) {
        GenParticle *daughter1 = (GenParticle*) branchGenParticle->At(particle->D1);
        d1_pid = daughter1 -> PID;
      }

      if (particle->D2 != -1) {
        GenParticle *daughter2 = (GenParticle*) branchGenParticle->At(particle->D2);
        d2_pid = daughter2 -> PID;
      }

      if( particle->PID == 23 && d1_pid != 23) ZPartonIndices.push_back(i);
    
    }

    // sort by mass
    if(ZPartonIndices.size() > 1){
      sort(ZPartonIndices.begin(), ZPartonIndices.end(), [branchGenParticle](const int& lhs, const int& rhs) {
        return ((GenParticle*)branchGenParticle->At(lhs))->P4().M() > (((GenParticle*)branchGenParticle->At(rhs))->P4()).M(); 
      });
    }


    return ZPartonIndices;

}

void getRecoZLeps(int& thisRecoEventType,  const vector<pair<int,pair<int,int>>> RecoPairIndices, TClonesArray *branchElectron, TClonesArray *branchMuon,  TLorentzVector &l1_reco,  TLorentzVector &l2_reco,  TLorentzVector& l3_reco,  TLorentzVector& l4_reco, int& q1_reco, int& q2_reco, int& q3_reco, int& q4_reco){
    
  if( thisRecoEventType==-1 ) return;

  if( thisRecoEventType==0 ) { // case 4mu 
    
    // take first two muons
    Muon *muon1_reco = (Muon *) branchMuon->At( RecoPairIndices[0].second.first);
    Muon *muon2_reco = (Muon *) branchMuon->At( RecoPairIndices[0].second.second);
    // take first two muons
    Muon *muon3_reco = (Muon *) branchMuon->At( RecoPairIndices[1].second.first);
    Muon *muon4_reco = (Muon *) branchMuon->At( RecoPairIndices[1].second.second);
    
    q1_reco = muon1_reco->Charge;
    q2_reco = muon2_reco->Charge;
    q3_reco = muon3_reco->Charge;
    q4_reco = muon4_reco->Charge;
    
    l1_reco=muon1_reco->P4();
    l2_reco=muon2_reco->P4();
    l3_reco=muon3_reco->P4();
    l4_reco=muon4_reco->P4(); 
    
  } else if( thisRecoEventType==1 ) { // case 4e
  
    // take first two electrons
    Electron *muon1_reco = (Electron *) branchElectron->At( RecoPairIndices[0].second.first);
    Electron *muon2_reco = (Electron *) branchElectron->At( RecoPairIndices[0].second.second);
    // take first two electrons
    Electron *muon3_reco = (Electron *) branchElectron->At( RecoPairIndices[1].second.first);
    Electron *muon4_reco = (Electron *) branchElectron->At( RecoPairIndices[1].second.second);

    q1_reco = muon1_reco->Charge;
    q2_reco = muon2_reco->Charge;
    q3_reco = muon3_reco->Charge;
    q4_reco = muon4_reco->Charge;

    l1_reco=muon1_reco->P4();
    l2_reco=muon2_reco->P4();
    l3_reco=muon3_reco->P4();
    l4_reco=muon4_reco->P4(); 

  } else if( thisRecoEventType==2 ) { // case 2mu2e 
  
    // take first two muons
    Muon *muon1_reco = (Muon *) branchMuon->At( RecoPairIndices[0].second.first);
    Muon *muon2_reco = (Muon *) branchMuon->At( RecoPairIndices[0].second.second);
    // take first two electrons
    Electron *muon3_reco = (Electron *) branchElectron->At( RecoPairIndices[1].second.first);
    Electron *muon4_reco = (Electron *) branchElectron->At( RecoPairIndices[1].second.second);

    q1_reco = muon1_reco->Charge;
    q2_reco = muon2_reco->Charge;
    q3_reco = muon3_reco->Charge;
    q4_reco = muon4_reco->Charge;

    l1_reco=muon1_reco->P4();
    l2_reco=muon2_reco->P4();
    l3_reco=muon3_reco->P4();
    l4_reco=muon4_reco->P4(); 

  } else if( thisRecoEventType==3 ) { // case 2e2mu 
  
    // take first two electrons
    Electron *muon1_reco = (Electron *) branchElectron->At( RecoPairIndices[0].second.first);
    Electron *muon2_reco = (Electron *) branchElectron->At( RecoPairIndices[0].second.second);
    // take first two muons
    Muon *muon3_reco = (Muon *) branchMuon->At( RecoPairIndices[1].second.first);
    Muon *muon4_reco = (Muon *) branchMuon->At( RecoPairIndices[1].second.second);

    q1_reco = muon1_reco->Charge;
    q2_reco = muon2_reco->Charge;
    q3_reco = muon3_reco->Charge;
    q4_reco = muon4_reco->Charge;

    l1_reco=muon1_reco->P4();
    l2_reco=muon2_reco->P4();
    l3_reco=muon3_reco->P4();
    l4_reco=muon4_reco->P4(); 

  }

}


void getParticleZLeps(int& thisParticleEventType,  const vector<pair<int,pair<int,int>>> ParticlePairIndices, TClonesArray *branchGenParticle,  TLorentzVector &l1_particle,  TLorentzVector &l2_particle,  TLorentzVector& l3_particle,  TLorentzVector& l4_particle, int& q1_particle, int& q2_particle, int& q3_particle, int& q4_particle){
  
  if( thisParticleEventType==-1 ) return;

  //if( thisParticleEventType==0 ) { // case 4mu 
  
    // take first two muons
    GenParticle *muon1_particle = (GenParticle*) branchGenParticle->At( ParticlePairIndices[0].second.first);
    GenParticle *muon2_particle = (GenParticle*) branchGenParticle->At( ParticlePairIndices[0].second.second);
    // take first two muons
    GenParticle *muon3_particle = (GenParticle*) branchGenParticle->At( ParticlePairIndices[1].second.first);
    GenParticle *muon4_particle = (GenParticle*) branchGenParticle->At( ParticlePairIndices[1].second.second);
  
    q1_particle = muon1_particle->Charge;
    q2_particle = muon2_particle->Charge;
    q3_particle = muon3_particle->Charge;
    q4_particle = muon4_particle->Charge;

    l1_particle=muon1_particle->P4();
    l2_particle=muon2_particle->P4();
    l3_particle=muon3_particle->P4();
    l4_particle=muon4_particle->P4(); 
  /*
  } else if( thisParticleEventType==1 ) { // case 4e
  
    // take first two electrons
    GenParticle *muon1_particle = (GenParticle*) branchGenParticle->At( ParticlePairIndices[0].second.first);
    GenParticle *muon2_particle = (GenParticle*) branchGenParticle->At( ParticlePairIndices[0].second.second);
    // take first two electrons
    GenParticle *muon3_particle = (GenParticle*) branchGenParticle->At( ParticlePairIndices[1].second.first);
    GenParticle *muon4_particle = (GenParticle*) branchGenParticle->At( ParticlePairIndices[1].second.second);

    q1_particle = muon1_particle->Charge;
    q2_particle = muon2_particle->Charge;
    q3_particle = muon3_particle->Charge;
    q4_particle = muon4_particle->Charge;

    l1_particle=muon1_particle->P4();
    l2_particle=muon2_particle->P4();
    l3_particle=muon3_particle->P4();
    l4_particle=muon4_particle->P4(); 

  } else if( thisParticleEventType==2 ) { // case 2mu2e 
  
    // take first two muons
    GenParticle *muon1_particle = (GenParticle*) branchGenParticle->At( ParticlePairIndices[0].second.first);
    GenParticle *muon2_particle = (GenParticle*) branchGenParticle->At( ParticlePairIndices[0].second.second);
    // take first two electrons
    GenParticle *muon3_particle = (GenParticle*) branchGenParticle->At( ParticlePairIndices[1].second.first);
    GenParticle *muon4_particle = (GenParticle*) branchGenParticle->At( ParticlePairIndices[1].second.second);

    q1_particle = muon1_particle->Charge;
    q2_particle = muon2_particle->Charge;
    q3_particle = muon3_particle->Charge;
    q4_particle = muon4_particle->Charge;

    l1_particle=muon1_particle->P4();
    l2_particle=muon2_particle->P4();
    l3_particle=muon3_particle->P4();
    l4_particle=muon4_particle->P4(); 

  } else if( thisParticleEventType==3 ) { // case 2e2mu
  
    // take first two electrons
    GenParticle *muon1_particle = (GenParticle*) branchGenParticle->At( ParticlePairIndices[0].second.first);
    GenParticle *muon2_particle = (GenParticle*) branchGenParticle->At( ParticlePairIndices[0].second.second);
    // take first two muons
    GenParticle *muon3_particle = (GenParticle*) branchGenParticle->At( ParticlePairIndices[1].second.first);
    GenParticle *muon4_particle = (GenParticle*) branchGenParticle->At( ParticlePairIndices[1].second.second);

    q1_particle = muon1_particle->Charge;
    q2_particle = muon2_particle->Charge;
    q3_particle = muon3_particle->Charge;
    q4_particle = muon4_particle->Charge;

    l1_particle=muon1_particle->P4();
    l2_particle=muon2_particle->P4();
    l3_particle=muon3_particle->P4();
    l4_particle=muon4_particle->P4(); 
  }
  */
}

void getPartonZLeps(int thisPartonEventType, vector <int> ZPartonIndices, TClonesArray *branchGenParticle,   TLorentzVector &z1_parton,  TLorentzVector &z2_parton, TLorentzVector &l1_parton,  TLorentzVector &l2_parton,  TLorentzVector& l3_parton,  TLorentzVector& l4_parton, int& q1_parton, int& q2_parton, int& q3_parton, int& q4_parton){
  
  vector <int> Z1children;
  Z1children.push_back( ((GenParticle*)branchGenParticle->At(ZPartonIndices[0]))->D1 ); 
  Z1children.push_back( ((GenParticle*)branchGenParticle->At(ZPartonIndices[0]))->D2 ); 

  vector <int> Z2children; 
  Z2children.push_back( ((GenParticle*)branchGenParticle->At(ZPartonIndices[1]))->D1 );
  Z2children.push_back( ((GenParticle*)branchGenParticle->At(ZPartonIndices[1]))->D2 ); 

  if (abs(((GenParticle*) branchGenParticle->At(Z1children[0]))->PID) == 13 && abs(((GenParticle*) branchGenParticle->At(Z1children[1]))->PID) == 13 && abs(((GenParticle*) branchGenParticle->At(Z2children[0]))->PID) == 13 && abs(((GenParticle*) branchGenParticle->At(Z2children[1]))->PID) == 13) thisPartonEventType = 0;
  else if (abs(((GenParticle*) branchGenParticle->At(Z1children[0]))->PID) == 11 && abs(((GenParticle*) branchGenParticle->At(Z1children[1]))->PID) == 11 && abs(((GenParticle*) branchGenParticle->At(Z2children[0]))->PID) == 11 && abs(((GenParticle*) branchGenParticle->At(Z2children[1]))->PID) == 11) thisPartonEventType = 1;
  else if (abs(((GenParticle*) branchGenParticle->At(Z1children[0]))->PID) == 13 && abs(((GenParticle*) branchGenParticle->At(Z1children[1]))->PID) == 13 && abs(((GenParticle*) branchGenParticle->At(Z2children[0]))->PID) == 11 && abs(((GenParticle*) branchGenParticle->At(Z2children[1]))->PID) == 11) thisPartonEventType = 2;
  else if (abs(((GenParticle*) branchGenParticle->At(Z1children[0]))->PID) == 11 && abs(((GenParticle*) branchGenParticle->At(Z1children[1]))->PID) == 11 && abs(((GenParticle*) branchGenParticle->At(Z2children[0]))->PID) == 13 && abs(((GenParticle*) branchGenParticle->At(Z2children[2]))->PID) == 13) thisPartonEventType = 3;
  else return;

  z1_parton = ((GenParticle*)branchGenParticle->At(ZPartonIndices[0]))->P4(); 
  z2_parton = ((GenParticle*)branchGenParticle->At(ZPartonIndices[1]))->P4(); 

  l1_parton = ((GenParticle*) branchGenParticle->At(Z1children[0]))->P4();
  l2_parton = ((GenParticle*) branchGenParticle->At(Z1children[1]))->P4();
  l3_parton = ((GenParticle*) branchGenParticle->At(Z2children[0]))->P4();
  l4_parton = ((GenParticle*) branchGenParticle->At(Z2children[1]))->P4();

  q1_parton = ((GenParticle*) branchGenParticle->At(Z1children[0]))->Charge;
  q2_parton = ((GenParticle*) branchGenParticle->At(Z1children[1]))->Charge;
  q3_parton = ((GenParticle*) branchGenParticle->At(Z2children[0]))->Charge;
  q4_parton = ((GenParticle*) branchGenParticle->At(Z2children[1]))->Charge;
    
}

//------------------------------------------------------------------------------------------------------------------------------------------------------------
// W
//------------------------------------------------------------------------------------------------------------------------------------------------------------

/*
   vector< pair<int,int>> GetWRecoIndices(TClonesArray *branchElectron=nullptr, TClonesArray *branchMuon=nullptr, vector<int> goodE_reco_indices=vector<int>(0), vector<int> goodMu_reco_indices=vector<int>(0)){
 
    vector< pair<int,int>> lepRecoIndices;

    for(int i=0; i<(int) goodE_reco_indices.size(); i++) lepRecoIndices.push_back(make_pair(0,goodE_reco_indices.at(i)));
    for(int i=0; i<(int) goodMu_reco_indices.size(); i++) lepRecoIndices.push_back(make_pair(1,goodMu_reco_indices.at(i)));

    sort(lepRecoIndices.begin(), lepRecoIndices.end(), [branchElectron, branchMuon](const pair<int, int>& lhs, const pair<int, int>& rhs) {
      if (lhs.first == 1 && rhs.first == 1) return ((Muon*)branchMuon->At(lhs.second))->PT > ((Muon*)branchMuon->At(rhs.second))->PT; // mu mu
      else if (lhs.first == 0 && rhs.first == 0) return ((Electron*)branchElectron->At(lhs.second))->PT > ((Electron*)branchElectron->At(rhs.second))->PT; // e e
      else if (lhs.first == 1 && rhs.first == 0) return ((Muon*)branchMuon->At(lhs.second))->PT > ((Electron*)branchElectron->At(rhs.second))->PT; // mu e
      else if (lhs.first == 0 && rhs.first == 1) return ((Electron*)branchElectron->At(lhs.second))->PT > ((Muon*)branchMuon->At(rhs.second))->PT; // e mu
      return false; 
    });

    return  lepRecoIndices;

  }
*/

/*
  vector< pair<int,int>> GetWParticleIndices(TClonesArray *branchGenParticle=nullptr, vector<int> goodE_particle_indices=vector<int>(0), vector<int> goodMu_particle_indices=vector<int>(0)){

    vector< pair<int,int>> lepParticleIndices;

    for(int i=0; i<(int) goodE_particle_indices.size(); i++) lepParticleIndices.push_back(make_pair(0,goodE_particle_indices.at(i)));
    for(int i=0; i<(int) goodMu_particle_indices.size(); i++) lepParticleIndices.push_back(make_pair(1,goodMu_particle_indices.at(i)));

    sort(lepParticleIndices.begin(), lepParticleIndices.end(), [branchGenParticle](const pair<int, int>& lhs, const pair<int, int>& rhs) {
      if (lhs.first == 1 && rhs.first == 1) return ((GenParticle*)branchGenParticle->At(lhs.second))->PT > ((GenParticle*)branchGenParticle->At(rhs.second))->PT; // mu mu
      else if (lhs.first == 0 && rhs.first == 0) return ((GenParticle*)branchGenParticle->At(lhs.second))->PT > ((GenParticle*)branchGenParticle->At(rhs.second))->PT; // e e
      else if (lhs.first == 1 && rhs.first == 0) return ((GenParticle*)branchGenParticle->At(lhs.second))->PT > ((GenParticle*)branchGenParticle->At(rhs.second))->PT; // mu e
      else if (lhs.first == 0 && rhs.first == 1) return ((GenParticle*)branchGenParticle->At(lhs.second))->PT > ((GenParticle*)branchGenParticle->At(rhs.second))->PT; // e mu
      return false;
    });

    return  lepParticleIndices;

  }
*/


vector<pair<pair<int,int>, int>> GetWRecoPairIndices(vector<int> goodE_reco_indices=vector<int>(0), vector<int> goodMu_reco_indices=vector<int>(0), TClonesArray *branchElectron=nullptr, TClonesArray *branchMuon=nullptr, TClonesArray *branchMissingET=nullptr){

  vector< pair<int,int>> lepRecoIndices;
  vector<pair<pair<int,int>, int>> WRecoPairIndices;

  for(int i=0; i<(int) goodE_reco_indices.size(); i++) lepRecoIndices.push_back(make_pair(0,goodE_reco_indices.at(i)));
  for(int i=0; i<(int) goodMu_reco_indices.size(); i++) lepRecoIndices.push_back(make_pair(1,goodMu_reco_indices.at(i)));

  //for(int i=0; i<(int) min(lepRecoIndices.size(),goodMET_reco_indices.size()); i++) WRecoPairIndices.push_back(make_pair(lepRecoIndices.at(i), goodMET_reco_indices.at(i)));

  // sort all of the indices by closeness to mW

  int METIndex=0;

  sort(lepRecoIndices.begin(), lepRecoIndices.end(), [branchMuon,branchElectron,branchMissingET,&METIndex]  ( const pair<int,int>   & lhs , const pair<int,int>  &rhs ){

    double massW1_reco=0;
    double massW2_reco=0;

    if( lhs.first == 0 ) massW1_reco = (((Electron*)branchElectron->At(lhs.second))->P4() + ((MissingET*)branchMissingET->At(METIndex))->P4()).Mt();
    else massW1_reco = (((Muon*)branchMuon->At(lhs.second))->P4() + ((MissingET*)branchMissingET->At(METIndex))->P4()).Mt();
    
    if( rhs.first == 0 ) massW2_reco = (((Electron*)branchElectron->At(rhs.second))->P4() + ((MissingET*)branchMissingET->At(METIndex))->P4()).Mt();
    else massW2_reco = (((Muon*)branchMuon->At(rhs.second))->P4() + ((MissingET*)branchMissingET->At(METIndex))->P4()).Mt();
    
    return fabs(massW1_reco - 80.0 )  <  fabs(massW2_reco - 80.0); 
    
  });

  for( int i=0; i<(int) lepRecoIndices.size(); i++) WRecoPairIndices.push_back(make_pair(make_pair(lepRecoIndices[i].first,lepRecoIndices[i].second),METIndex));
  return WRecoPairIndices;

}


vector<pair<pair<int,int>, int>> GetWParticlePairIndices(vector<int> goodE_particle_indices=vector<int>(0), vector<int> goodMu_particle_indices=vector<int>(0), TClonesArray *branchGenParticle=nullptr, TClonesArray *branchMissingET=nullptr){

  vector< pair<int,int>> lepParticleIndices;
  vector<pair<pair<int,int>, int>> WParticlePairIndices;

  for(int i=0; i<(int) goodE_particle_indices.size(); i++) lepParticleIndices.push_back(make_pair(0,goodE_particle_indices.at(i)));
  for(int i=0; i<(int) goodMu_particle_indices.size(); i++) lepParticleIndices.push_back(make_pair(1,goodMu_particle_indices.at(i)));

  //for(int i=0; i<(int) min(lepParticleIndices.size(),goodMET_reco_indices.size()); i++) WParticlePairIndices.push_back(make_pair(lepParticleIndices.at(i), goodMET_reco_indices.at(i)));

  // sort all of the indices by closeness to mW

  int METIndex=0;

  sort(lepParticleIndices.begin(), lepParticleIndices.end(), [branchGenParticle,branchMissingET,&METIndex]  ( const pair<int,int>   & lhs , const pair<int,int>  &rhs ){

    double massW1_particle=0;
    double massW2_particle=0;

    if( lhs.first == 0 ) massW1_particle = (((GenParticle*)branchGenParticle->At(lhs.second))->P4() + ((MissingET*)branchMissingET->At(METIndex))->P4()).Mt();
    else massW1_particle = (((GenParticle*)branchGenParticle->At(lhs.second))->P4() + ((MissingET*)branchMissingET->At(METIndex))->P4()).Mt();
    
    if( rhs.first == 0 ) massW2_particle = (((GenParticle*)branchGenParticle->At(rhs.second))->P4() + ((MissingET*)branchMissingET->At(METIndex))->P4()).Mt();
    else massW2_particle = (((GenParticle*)branchGenParticle->At(rhs.second))->P4() + ((MissingET*)branchMissingET->At(METIndex))->P4()).Mt();
    
    return fabs(massW1_particle - 80.0 )  <  fabs(massW2_particle - 80.0); 
    
  });

  for( int i=0; i<(int) lepParticleIndices.size(); i++) WParticlePairIndices.push_back(make_pair(make_pair(lepParticleIndices[i].first,lepParticleIndices[i].second),METIndex));
  return WParticlePairIndices;

}

vector <int> GetWPartonIndices(TClonesArray *branchGenParticle=nullptr, string analysis="HWWJJ") {

    vector <int> WPartonIndices; 

    for(int i=0; i<(int)branchGenParticle->GetEntries(); i++){

      GenParticle *particle=(GenParticle*)branchGenParticle->At(i);
      int d1_pid = 9999;

      if (particle->D1 != -1) {
        GenParticle *daughter1 = (GenParticle*) branchGenParticle->At(particle->D1);
        d1_pid = daughter1 -> PID;
      }

      if( abs(particle->PID) == 24 && abs(d1_pid) != 24) WPartonIndices.push_back(i);
      
    
    }

    // sort by mass
    if(WPartonIndices.size() > 1){
      sort(WPartonIndices.begin(), WPartonIndices.end(), [branchGenParticle](const int& lhs, const int& rhs) {
        return ((GenParticle*)branchGenParticle->At(lhs))->P4().M() > (((GenParticle*)branchGenParticle->At(rhs))->P4()).M(); 
      });
    }

    return WPartonIndices;

}

void getWReco(int& thisRecoEventType,  const vector<pair<pair<int,int>, int>> WRecoPairIndices, TClonesArray *branchElectron, TClonesArray *branchMuon, TClonesArray *branchMissingET, TLorentzVector &l1_reco, TLorentzVector &l2_reco, int& q1_reco, int& q2_reco,  TLorentzVector& met1, TLorentzVector& met2){
    
    if( thisRecoEventType==-1 ) return;

    MissingET *Met1 = (MissingET *) branchMissingET->At( WRecoPairIndices[0].second);
    MissingET *Met2 = (MissingET *) branchMissingET->At( WRecoPairIndices[1].second);

    met1.SetPtEtaPhiM((Met1->MET)/2,0,Met1->Phi,0);
    met2.SetPtEtaPhiM((Met2->MET)/2,0,Met2->Phi,0);

    if( thisRecoEventType==0 ) { // case 2mu 
      
      Muon *muon1_reco = (Muon *) branchMuon->At( WRecoPairIndices[0].first.second);
      Muon *muon2_reco = (Muon *) branchMuon->At( WRecoPairIndices[1].first.second);
      
      q1_reco = muon1_reco->Charge;
      q2_reco = muon2_reco->Charge;
      
      l1_reco=muon1_reco->P4();
      l2_reco=muon2_reco->P4();
      
    } else if( thisRecoEventType==1 ) { // case 2e
    
      Electron *muon1_reco = (Electron *) branchElectron->At( WRecoPairIndices[0].first.second);
      Electron *muon2_reco = (Electron *) branchElectron->At( WRecoPairIndices[1].first.second);

      q1_reco = muon1_reco->Charge;
      q2_reco = muon2_reco->Charge;

      l1_reco=muon1_reco->P4();
      l2_reco=muon2_reco->P4();

    } else if( thisRecoEventType==2 ) { // case 1mu1e 
    
      Muon *muon1_reco = (Muon *) branchMuon->At( WRecoPairIndices[0].first.second);
      Electron *muon2_reco = (Electron *) branchElectron->At( WRecoPairIndices[1].first.second);

      q1_reco = muon1_reco->Charge;
      q2_reco = muon2_reco->Charge;

      l1_reco=muon1_reco->P4();
      l2_reco=muon2_reco->P4();

    } else if( thisRecoEventType==3 ) { // case 1e1mu 
    
      Electron *muon1_reco = (Electron *) branchElectron->At( WRecoPairIndices[0].first.second);
      Muon *muon2_reco = (Muon *) branchMuon->At( WRecoPairIndices[1].first.second);

      q1_reco = muon1_reco->Charge;
      q2_reco = muon2_reco->Charge;

      l1_reco=muon1_reco->P4();
      l2_reco=muon2_reco->P4();

    }

}


void getWParticle(int& thisParticleEventType,  const vector<pair<pair<int,int>,int>> WParticlePairIndices, TClonesArray *branchGenParticle, TClonesArray *branchMissingET, TLorentzVector &l1_particle, TLorentzVector &l2_particle, int& q1_particle, int& q2_particle,  TLorentzVector& met1, TLorentzVector& met2){    
    
    if( thisParticleEventType==-1 ) return;

    MissingET *Met1 = (MissingET *) branchMissingET->At( WParticlePairIndices[0].second);
    MissingET *Met2 = (MissingET *) branchMissingET->At( WParticlePairIndices[1].second);

    met1.SetPtEtaPhiM((Met1->MET)/2,0,Met1->Phi,0);
    met2.SetPtEtaPhiM((Met2->MET)/2,0,Met2->Phi,0);

  //  if( thisParticleEventType==0 ) { // case 2mu 
      
      GenParticle *muon1_particle = (GenParticle *) branchGenParticle->At( WParticlePairIndices[0].first.second);
      GenParticle *muon2_particle = (GenParticle *) branchGenParticle->At( WParticlePairIndices[1].first.second);
      
      q1_particle = muon1_particle->Charge;
      q2_particle = muon2_particle->Charge;
      
      l1_particle=muon1_particle->P4();
      l2_particle=muon2_particle->P4();
  /*   
   } else if( thisParticleEventType==1 ) { // case 2e
    
      GenParticle *muon1_particle = (GenParticle *) branchGenParticle->At( WParticlePairIndices[0].first.second);
      GenParticle *muon2_particle = (GenParticle *) branchGenParticle->At( WParticlePairIndices[0].first.second);

      q1_particle = muon1_particle->Charge;
      q2_particle = muon2_particle->Charge;

      l1_particle=muon1_particle->P4();
      l2_particle=muon2_particle->P4();

    } else if( thisParticleEventType==2 ) { // case 1mu1e 
    
      GenParticle *muon1_particle = (GenParticle *) branchGenParticle->At( WParticlePairIndices[0].first.second);
      GenParticle *muon2_particle = (GenParticle *) branchGenParticle->At( WParticlePairIndices[0].first.second);

      q1_particle = muon1_particle->Charge;
      q2_particle = muon2_particle->Charge;

      l1_particle=muon1_particle->P4();
      l2_particle=muon2_particle->P4();

    } else if( thisParticleEventType==3 ) { // case 1e1mu 
    
      GenParticle *muon1_particle = (GenParticle *) branchGenParticle->At( WParticlePairIndices[0].first.second);
      GenParticle *muon2_particle = (GenParticle *) branchGenParticle->At( WParticlePairIndices[0].first.second);

      q1_particle = muon1_particle->Charge;
      q2_particle = muon2_particle->Charge;

      l1_particle=muon1_particle->P4();
      l2_particle=muon2_particle->P4();
    
    }
  */
}


void getPartonWLeps(int& thisPartonEventType, vector <int> WPartonIndices, TClonesArray *branchGenParticle,   TLorentzVector &w1_parton,  TLorentzVector &w2_parton, TLorentzVector &l1_parton,  TLorentzVector &l2_parton, int& q1_parton, int& q2_parton){
  
  vector <int> W1childrenIndices;
  vector <int> W2childrenIndices;

  // get w
  GenParticle *w1=(GenParticle*) branchGenParticle->At(WPartonIndices[0]);
  GenParticle *w2=(GenParticle*) branchGenParticle->At(WPartonIndices[1]);

  // get w daughters + pid
  GenParticle *w1daughter1 = (GenParticle*) branchGenParticle->At(w1->D1);
  GenParticle *w1daughter2 = (GenParticle*) branchGenParticle->At(w1->D2);
  GenParticle *w2daughter1 = (GenParticle*) branchGenParticle->At(w2->D1);
  GenParticle *w2daughter2 = (GenParticle*) branchGenParticle->At(w2->D2);
  int w1d1_pid = w1daughter1 -> PID;
  int w1d2_pid = w1daughter2 -> PID;
  int w2d1_pid = w2daughter1 -> PID;
  int w2d2_pid = w2daughter2 -> PID;

  // get the daughter that is not a neutrino
  if(abs(w1d1_pid) == 11 || abs(w1d1_pid) == 13) W1childrenIndices.push_back(w1->D1);
  else if(abs(w1d2_pid) == 11 || abs(w1d2_pid) == 13) W1childrenIndices.push_back(w1->D2);
  else return;

  if(abs(w2d1_pid) == 11 || abs(w2d1_pid) == 13) W2childrenIndices.push_back(w2->D1);
  else if(abs(w2d2_pid) == 11 || abs(w2d2_pid) == 13) W2childrenIndices.push_back(w2->D2);
  else return;

  // sort here ?

  // get event type
  GenParticle *w1lep = (GenParticle*) branchGenParticle->At(W1childrenIndices[0]);
  GenParticle *w2lep = (GenParticle*) branchGenParticle->At(W2childrenIndices[0]);
  int w1lep_pid = w1lep -> PID;
  int w2lep_pid = w2lep -> PID;

  if (abs(w1lep_pid) == 13 && abs(w2lep_pid) == 13) thisPartonEventType = 0;
  else if (abs(w1lep_pid) == 11 && abs(w2lep_pid) == 11) thisPartonEventType = 1;
  else if (abs(w1lep_pid) == 13 && abs(w2lep_pid) == 11) thisPartonEventType = 2;
  else if (abs(w1lep_pid) == 11 && abs(w2lep_pid) == 13) thisPartonEventType = 3;
  else return;

  w1_parton = ((GenParticle*)branchGenParticle->At(WPartonIndices[0]))->P4(); 
  w2_parton = ((GenParticle*)branchGenParticle->At(WPartonIndices[1]))->P4(); 

  l1_parton = ((GenParticle*) branchGenParticle->At(W1childrenIndices[0]))->P4();
  l2_parton = ((GenParticle*) branchGenParticle->At(W2childrenIndices[0]))->P4();

  q1_parton = ((GenParticle*) branchGenParticle->At(W1childrenIndices[0]))->Charge;
  q2_parton = ((GenParticle*) branchGenParticle->At(W2childrenIndices[0]))->Charge;

}

#endif