#ifndef SELECTIONS_H
#define SELECTIONS_H
#include "HepMC/GenParticle.h"
#include "classes/DelphesClasses.h"
#include "classes/DelphesLHEFReader.h"
#include "external/ExRootAnalysis/ExRootTreeReader.h"
#include "../common_includes/ghost_tagging.h"
#include "../common_includes/get_cross_section.h"
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







#endif
