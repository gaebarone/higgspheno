//#ifndef LEPANALYZER_H
//#define LEPANALYZER_H

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

vector <int> GoodElectronIndices(TClonesArray *branchElectron=nullptr){
  vector <int> goodE_reco_indices;
  if( branchElectron==nullptr) return goodE_reco_indices;
  for(int i=0; i<(int)branchElectron->GetEntries(); i++){
      Electron *el_reco = (Electron *) branchElectron->At(i);
      // pT and eta cuts 
      if( el_reco->PT > 5.0 && fabs(el_reco->Eta) < 2.5) goodE_reco_indices.push_back(i);
  }
  // sort the indices by pT ;
  sort(goodE_reco_indices.begin(), goodE_reco_indices.end(), [branchElectron](const int& lhs, const int& rhs) {
      return ((Electron*)branchElectron->At(lhs))->PT > ((Electron*)branchElectron->At(rhs))->PT;
    });
  
  return goodE_reco_indices; 
  

}

vector <int> GoodMuonIndices(TClonesArray *branchMuon=nullptr){
  vector <int> goodMu_reco_indices; 
  if( branchMuon==nullptr) return goodMu_reco_indices; 

  
  for(int i=0; i<(int)branchMuon->GetEntries(); i++){
      Muon *mu_reco = (Muon *) branchMuon->At(i);
    
      // pT and eta cuts 
      if( mu_reco->PT > 5.0 && fabs(mu_reco->Eta) < 2.5) goodMu_reco_indices.push_back(i);

    }

    // sort the indices by pT ;
    sort(goodMu_reco_indices.begin(), goodMu_reco_indices.end(), [branchMuon](const int& lhs, const int& rhs) {
	return ((Muon*)branchMuon->At(lhs))->PT > ((Muon*)branchMuon->At(rhs))->PT;
      });
    return goodMu_reco_indices; 
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


vector< pair<int,int>> GetmuRecoPairIndices(TClonesArray *branchMuon=nullptr, vector<int> goodMu_reco_indices=vector<int>(0)){
    vector< pair<int,int>> muRecoPairIndices;
    vector< pair<int,int>> muRecoPairIndicesIn;
    vector <vector<int>> muRecoPairIndices_;
    if(goodMu_reco_indices.size() > 1) 
      muRecoPairIndices_=combinationsNoRepetitionAndOrderDoesNotMatter(2,goodMu_reco_indices);
    
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

vector<pair<int,pair<int,int>>> GetRecoPairIndices(vector< pair<int,int>> elecRecoPairIndices, vector< pair<int,int>> muRecoPairIndices,int &thisRecoEventType,TClonesArray *branchElectron=nullptr,TClonesArray *branchMuon=nullptr){

  vector<pair<int,pair<int,int>>> RecoPairIndices;


for(int i=0; i<(int) elecRecoPairIndices.size(); i++)
      RecoPairIndices.push_back(make_pair(0,elecRecoPairIndices.at(i)));
    for(int i=0; i<(int) muRecoPairIndices.size(); i++)
      RecoPairIndices.push_back(make_pair(1,muRecoPairIndices.at(i)));

    // sort all of the indices by closeness to mZ
    sort(RecoPairIndices.begin(), RecoPairIndices.end(), [branchMuon,branchElectron]  ( const pair<int,pair<int,int>>  lhs , const pair<int,pair<int,int>>   rhs ){

	double mass1_reco=0 ;
  double mass2_reco=0;
	
	if( lhs.first==0 ) mass1_reco = (((Electron*)branchElectron->At(lhs.second.first))->P4() + ((Electron*)branchElectron->At(lhs.second.second))->P4()).M();
	else mass1_reco = (((Muon*)branchMuon->At(lhs.second.first))->P4() + ((Muon*)branchMuon->At(lhs.second.second))->P4()).M();
	
	if( rhs.first==0 ) mass2_reco = (((Electron*)branchElectron->At(rhs.second.first))->P4() + ((Electron*)branchElectron->At(rhs.second.second))->P4()).M();
	else mass2_reco = (((Muon*)branchMuon->At(rhs.second.first))->P4() + ((Muon*)branchMuon->At(rhs.second.second))->P4()).M();
	
	return fabs(mass1_reco -91.0 )  <  fabs(mass2_reco -91.0); 
	
      });

     //cout<<"Reco pair indices "<<RecoPairIndices.size()<<endl;
    //for( int i=0; i<(int)RecoPairIndices.size(); i++){
    //if(RecoPairIndices[i].first==0) 
    //	cout<<" "<<"electrons "<< (((Electron*)branchElectron->At(RecoPairIndices[i].second.first))->P4() + ((Electron*)branchElectron->At(RecoPairIndices[i].second.second))->P4()).M()<<" "<<i<<endl;
    //else
    //	cout<<" "<<"muons "<< (((Muon*)branchMuon->At(RecoPairIndices[i].second.first))->P4() + ((Muon*)branchMuon->At(RecoPairIndices[i].second.second))->P4()).M()<<"  "<<i<<endl;
    //}
    
    return RecoPairIndices;
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

//------------------------------------------------------------------------------------------------------------------------------------------------------------
// W
//------------------------------------------------------------------------------------------------------------------------------------------------------------

   vector< pair<int,int>> GetLepRecoIndices(TClonesArray *branchElectron=nullptr, TClonesArray *branchMuon=nullptr, vector<int> goodE_reco_indices=vector<int>(0), vector<int> goodMu_reco_indices=vector<int>(0)){
 
    vector< pair<int,int>> lepRecoIndices;

    for(int i=0; i<(int) goodE_reco_indices.size(); i++) lepRecoIndices.push_back(make_pair(0,goodE_reco_indices.at(i)));
    for(int i=0; i<(int) goodMu_reco_indices.size(); i++) lepRecoIndices.push_back(make_pair(1,goodMu_reco_indices.at(i)));

    sort(lepRecoIndices.begin(), lepRecoIndices.end(), [branchElectron, branchMuon](const pair<int, int>& lhs, const pair<int, int>& rhs) {
      if (lhs.first == 1 && rhs.first == 1) return ((Muon*)branchMuon->At(lhs.second))->PT > ((Muon*)branchMuon->At(rhs.second))->PT; // mu mu
      else if (lhs.first == 0 && rhs.first == 0) return ((Electron*)branchElectron->At(lhs.second))->PT > ((Electron*)branchElectron->At(rhs.second))->PT; // e e
      else if (lhs.first == 1 && rhs.first == 0) return ((Muon*)branchMuon->At(lhs.second))->PT > ((Electron*)branchElectron->At(rhs.second))->PT; // mu e
      else if (lhs.first == 0 && rhs.first == 1) return ((Electron*)branchElectron->At(lhs.second))->PT > ((Muon*)branchMuon->At(rhs.second))->PT; // e mu
    });

    return  lepRecoIndices;

  }

void getRecoWLeps(int& thisRecoEventType,  const vector< pair<int,int>> WRecoIndices, TClonesArray *branchElectron, TClonesArray *branchMuon,  TLorentzVector &l1_reco,  TLorentzVector &l2_reco, int& q1_reco, int& q2_reco){
    
    if( thisRecoEventType==-1 ) return;

    if( thisRecoEventType==0 ) { // case 2mu 
      
      Muon *muon1_reco = (Muon *) branchMuon->At( WRecoIndices[0].second);
      Muon *muon2_reco = (Muon *) branchMuon->At( WRecoIndices[1].second);
      
      q1_reco = muon1_reco->Charge;
      q2_reco = muon2_reco->Charge;
      
      l1_reco=muon1_reco->P4();
      l2_reco=muon2_reco->P4();
      
    } else if( thisRecoEventType==1 ) { // case 2e
    
      Electron *muon1_reco = (Electron *) branchElectron->At( WRecoIndices[0].second);
      Electron *muon2_reco = (Electron *) branchElectron->At( WRecoIndices[0].second);

      q1_reco = muon1_reco->Charge;
      q2_reco = muon2_reco->Charge;

      l1_reco=muon1_reco->P4();
      l2_reco=muon2_reco->P4();

    } else if( thisRecoEventType==2 ) { // case 1mu1e 
    
      Muon *muon1_reco = (Muon *) branchMuon->At( WRecoIndices[0].second);
      Electron *muon2_reco = (Electron *) branchElectron->At( WRecoIndices[0].second);

      q1_reco = muon1_reco->Charge;
      q2_reco = muon2_reco->Charge;

      l1_reco=muon1_reco->P4();
      l2_reco=muon2_reco->P4();

    } else if( thisRecoEventType==3 ) { // case 1e1mu 
    
      Electron *muon1_reco = (Electron *) branchElectron->At( WRecoIndices[0].second);
      Muon *muon2_reco = (Muon *) branchMuon->At( WRecoIndices[0].second);

      q1_reco = muon1_reco->Charge;
      q2_reco = muon2_reco->Charge;

      l1_reco=muon1_reco->P4();
      l2_reco=muon2_reco->P4();

    }

}