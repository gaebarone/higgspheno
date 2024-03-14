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
