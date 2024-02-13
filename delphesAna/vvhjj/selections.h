#ifndef SELECTIONS_H
#define SELECTIONS_H
#include "HepMC/GenParticle.h"
#include "classes/DelphesClasses.h"
#include "classes/DelphesLHEFReader.h"
#include "external/ExRootAnalysis/ExRootTreeReader.h"
#include "../common_includes/ghost_tagging.h"
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


double btagEff=0.85;
double fakeEff=0.01;


void remove_overlaps(vector< pair<int,int>> muPairIndices){
  for( vector< pair<int,int>>::iterator it=muPairIndices.begin(); it!=muPairIndices.end(); it++){
    pair<int,int> one=(*it);
    for(vector< pair<int,int>>::iterator it2=it+1; it2!=muPairIndices.end(); it2++){
      pair<int,int> two=(*it2);
      if( one.first == two.first || one.first == two.second || one.second == two.second || one.second == two.second) muPairIndices.erase(it2--);
    }
  }
}


pair<int,int> GethiggsbbcandidateNoMass(vector<vector <int>>  &bJetPairsComb,
				   vector<pair<int,int>> &bJetPairs,
				   pair<int,int> &b12pos,
				   bool & foundBjet, 
				   TClonesArray *branchJet=nullptr,
				   TClonesArray *branchGenParticle=nullptr){
  
  pair <int,int> higgsbbcandidate;
 
  for(int i=0; i<(int)bJetPairsComb.size(); i++)
    bJetPairs.push_back(make_pair(bJetPairsComb[i][0],bJetPairsComb[i][1]));
  
  // sort the pairs of b-tags by pT;
  //if( bJetPairs.size() > 1){
  //sort(bJetPairs.begin(), bJetPairs.end(), [branchJet](const pair<int,int> lhs, const pair<int,int> rhs) {
  //  return ((Jet*)branchJet->At(lhs.first))->PT > ((Jet*)branchJet->At(rhs.first))->PT; 
  //});
  //}
  
  for(int i=0; i<(int) bJetPairs.size(); i++){
    higgsbbcandidate=bJetPairs[i];
    // Attention
    Jet *b1=(Jet*)branchJet->At(b12pos.first);
    Jet *b2=(Jet*)branchJet->At(b12pos.second);
    foundBjet=true;
    break; 
  }
  
    return higgsbbcandidate;
  

}

pair <int,int> Gethiggsbbcandidate(vector<vector <int>>  &bJetPairsComb,
				   vector<pair<int,int>> &bJetPairs,
				   pair<int,int> &b12pos,
				   bool & foundBjet, 
				   TClonesArray *branchJet=nullptr,
				   TClonesArray *branchGenParticle=nullptr){
  
  pair <int,int> higgsbbcandidate;
 
    for(int i=0; i<(int)bJetPairsComb.size(); i++)
      bJetPairs.push_back(make_pair(bJetPairsComb[i][0],bJetPairsComb[i][1]));
   
    if( bJetPairs.size() > 1) 
      sort(bJetPairs.begin(), bJetPairs.end(), [branchJet](const pair<int,int> lhs, const pair<int,int> rhs) {
	  return fabs(((((Jet*)branchJet->At(lhs.first))->P4() + ((Jet*)branchJet->At(lhs.second))->P4())).M() - 125 ) <
	    fabs( ((((Jet*)branchJet->At(rhs.first))->P4() + ((Jet*)branchJet->At(rhs.second))->P4()).M()) - 125 ) ; 
	});
    
    for(int i=0; i<(int) bJetPairs.size(); i++){
      //afs/  b12pos=make_pair(bJetPairs[i].first,(bJetPairs[i].second));
      //*b1=(Jet*)branchJet->At(bJetPairs[i].first);
      //*b2=(Jet*)branchJet->At(bJetPairs[i].second);
      //if( b1->BTag>0 && b2->BTag>0) {
      // Attention
      Jet *b1=(Jet*)branchJet->At(b12pos.first);
      Jet *b2=(Jet*)branchJet->At(b12pos.second);
      
      if( isMyBTag(b1, branchGenParticle,0,0.4,btagEff,fakeEff) && abs(b1->Eta) < 2.5 || (isMyBTag(b2, branchGenParticle,0,0.4,btagEff,fakeEff) && abs(b2->Eta)<2.5) ) {
	higgsbbcandidate=bJetPairs[i];
	foundBjet=true;
	break;
      }
    }
  
    return higgsbbcandidate;
}


vector<pair <int,double>> JetBtagScoreIndexParicle(vector <int> goodJetIndex,
						   TClonesArray *branchGenJet=nullptr,
						   TClonesArray *branchGenParticle=nullptr){
  
  
  vector<pair <int,double>> scores;
  for( unsigned int i=0; i<goodJetIndex.size(); i++){
    
    scores.push_back( make_pair(goodJetIndex[i],
				ghost_btagPseudoRecoScore(branchGenParticle,
							  dynamic_cast <Jet*> (branchGenJet->At(i)),
							  0.4)));
    
   
  }
  

  sort( scores.begin(), scores.end(), [branchGenJet](const pair<int,double> lhs, const pair <int,double> rhs){
    return lhs.second > rhs.second;
  });
  return scores; 
}

vector<pair <int,double>>  JetBtagScoreIndex(std::vector <int> goodJetIndex,
					    TClonesArray *branchJet=nullptr,
					    TClonesArray *branchGenParticle=nullptr){
  
  
  vector<pair <int,double>> scores;;
  for( unsigned int i=0; i<goodJetIndex.size(); i++){
    
    scores.push_back( make_pair(goodJetIndex[i],
					 ghost_btagPseudoRecoScoreSmeared(branchGenParticle,
									  (Jet*)branchJet->At(i),
									  0,
									  0.4,
									  0.9,
									  0.01,
									  0.01)));
    
  }
  

  sort( scores.begin(), scores.end(), [branchJet](const pair<int,double> lhs, const pair <int,double> rhs){
      return lhs.second > rhs.second;
  });
  return scores; 
}


vector <int> GoodJetIndices( vector <int> & btagIndex, 
			     vector <int> & noBtag,
			     TClonesArray *branchJet=nullptr,
			     TClonesArray *branchGenParticle=nullptr){
  vector <int> goodJetIndex;

  
    for(int i=0; i<(int)branchJet->GetEntries(); i++){
      Jet *jet=(Jet*) branchJet->At(i);
      //if( jet->PT < 20) continue;
      //if (fabs(jet->Eta) > 4.4) continue; 
      //  if( jet->BTag>0) btagIndex.push_back(i);
      if( isMyBTag(jet, branchGenParticle,0,0.4,btagEff,fakeEff) && abs(jet->Eta) < 2.5 ) btagIndex.push_back(i); 

      else noBtag.push_back(i);
      goodJetIndex.push_back(i);
    }
    sort(btagIndex.begin(), btagIndex.end(), [branchJet](const int& lhs, const int& rhs) {
	return ((Jet*)branchJet->At(lhs))->PT > ((Jet*)branchJet->At(rhs))->PT;
      });
    sort(noBtag.begin(), noBtag.end(), [branchJet](const int& lhs, const int& rhs) {
	return ((Jet*)branchJet->At(lhs))->PT > ((Jet*)branchJet->At(rhs))->PT;
      });
    sort(goodJetIndex.begin(), goodJetIndex.end(), [branchJet](const int& lhs, const int& rhs) {
	return ((Jet*)branchJet->At(lhs))->PT > ((Jet*)branchJet->At(rhs))->PT;
      });
    
    
  return goodJetIndex; 
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
void SortByPtIndices(vector <int> nonHiggsJet, TClonesArray *branchJet){
  if( nonHiggsJet.size() > 1 ){
    sort(nonHiggsJet.begin(), nonHiggsJet.end(), [branchJet](const int& lhs, const int& rhs) {
	return ((Jet*)branchJet->At(lhs))->PT > ((Jet*)branchJet->At(rhs))->PT;
      });}  
}

vector<pair<int,int>> GetvbfJetIndex( vector<vector <int>> vbfJetIndexComb){
  vector<pair<int,int>> vbfJetIndex;
  for(int i=0; i<(int)vbfJetIndexComb.size(); i++)
    vbfJetIndex.push_back(make_pair(vbfJetIndexComb[i][0],vbfJetIndexComb[i][1]));
  return vbfJetIndex;
}



#endif
