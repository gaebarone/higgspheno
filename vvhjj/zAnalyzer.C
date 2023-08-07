#ifdef __CLING__
R__LOAD_LIBRARY(libDelphes)
#include "classes/DelphesClasses.h"
#include "external/ExRootAnalysis/ExRootTreeReader.h"
#include "combinations.h"     
#include "DSCBf.h"
#include "TH1.h"
#include "TCanvas.h"
#endif

//------------------------------------------------------------------------------------------------------------------------------------------------------------
// initialize combinations
//------------------------------------------------------------------------------------------------------------------------------------------------------------

vector <string> types={"4mu","4e", "2mu2e", "2e2mu"};

void remove_overlaps(vector< pair<int,int>> muPairIndices){
  for( vector< pair<int,int>>::iterator it=muPairIndices.begin(); it!=muPairIndices.end(); it++){
    pair<int,int> one=(*it);
    for(vector< pair<int,int>>::iterator it2=it+1; it2!=muPairIndices.end(); it2++){
      pair<int,int> two=(*it2);
      if( one.first == two.first || one.first == two.second || one.second == two.second || one.second == two.second) muPairIndices.erase(it2--);
    }
  }
}

std::vector<std::vector<int>> combinationsNoRepetitionAndOrderDoesNotMatter (int subsetSize, std::vector<int> setOfNumbers){
  std::vector<std::vector<int> > subsets{};
  subsets.reserve (count_each_combination (setOfNumbers.begin (), setOfNumbers.begin () + subsetSize, setOfNumbers.end ()));
  for_each_combination (setOfNumbers.begin (), setOfNumbers.begin () + subsetSize, setOfNumbers.end (), [&subsets] (auto first, auto last) {
    subsets.push_back (std::vector<int>{ first, last });
    return false;
  });
  return subsets;
}

std::vector <std::vector <int>> getAllCombinations(vector<int> inputVector, int k){
  std::vector<vector<int>> combinations; 
  std::vector<int> selector(inputVector.size());
    std::fill(selector.begin(), selector.begin() + k, 1);

    do {
        std::vector<int> selectedIds;
        std::vector<int> selectedVectorElements;
        for (int i = 0; i < inputVector.size(); i++) {
            if (selector[i]) {
                selectedIds.push_back(i);
            }
        }
        for (auto& id : selectedIds) {
            selectedVectorElements.push_back(inputVector[id]);
        }
        combinations.push_back(selectedVectorElements);
    } while (std::prev_permutation(selector.begin(), selector.end()));

    return combinations;
}

//------------------------------------------------------------------------------------------------------------------------------------------------------------
// histogram settings
//------------------------------------------------------------------------------------------------------------------------------------------------------------

void PrintCanvas(TCanvas *c=nullptr,string name="default"){
  std::vector <string> types={"jpg"}; 
  for(std::vector<string>::iterator it=types.begin(); it!=types.end(); it++) {
    c->Print(Form("zzhjj_plots/%s.%s",name.c_str(),(*it).c_str()),(*it).c_str());
  }
}

void draw_hist(TH1 *histo, const char *name, const char *title, const char *axistitle) {
  TCanvas *c = new TCanvas(name, title, 1500,1200);
  c->cd();
  histo->GetXaxis()->SetTitle(axistitle);
  histo->Draw("axis");
//  PrintCanvas(c, name);
}

void draw_hist2(TH2 *histo, const char *name, const char *title, const char *xaxistitle, const char *yaxistitle) {
  TCanvas *c = new TCanvas(name, title, 1500, 1200);
  histo->GetXaxis()->SetTitle(xaxistitle);
  histo->GetYaxis()->SetTitle(yaxistitle);
  histo->Draw("COLZ");
//  PrintCanvas(c, name);
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
    return NULL;
  }
}

//------------------------------------------------------------------------------------------------------------------------------------------------------------
// z analyzer
//------------------------------------------------------------------------------------------------------------------------------------------------------------


void zAnalyzer(const char *inputFile,const char *outputFile) {
  gSystem->Load("libDelphes");

  TChain chain("Delphes");
  chain.Add(inputFile);

  ExRootTreeReader *treeReader = new ExRootTreeReader(&chain);
  Long64_t numberOfEntries = treeReader->GetEntries();

  TClonesArray *branchJet = treeReader->UseBranch("Jet");
  TClonesArray *branchElectron = treeReader->UseBranch("Electron");
  TClonesArray *branchMuon = treeReader->UseBranch("Muon");
  TClonesArray *branchEvent = treeReader->UseBranch("Event");
  TClonesArray *branchGenParticle = treeReader->UseBranch("Particle");
  TClonesArray *branchGenJet = treeReader->UseBranch("GenJet");

// higgs
    // 1D
    // reco 
        // pT + m
        TH1 *hHpTreco = new TH1F("H_pT_reco", "p^{T}_{h}", 100, 0, 500);
        TH1 *hHmreco = new TH1F("H_m_reco", "m_{h}", 50, 0, 250);
    // particle
        // pT + m
        TH1 *hHpTparticle = new TH1F("H_pT_particle", "p^{T}_{h}", 100, 0, 500);
        TH1 *hHmparticle = new TH1F("H_m_particle", "m_{h}", 50, 0, 250);
    // parton
        // pT + m
        TH1 *hHpTparton = new TH1F("H_pT_parton", "p^{T}_{h}", 100, 0, 500);
        TH1 *hHmparton = new TH1F("H_m_parton", "m_{h}", 50, 0,  250);
    //2D
        TH2 *hHpTComp = new TH2F("H_pT_comp", "p_{T}^{h}", 20, 0.0, 500, 20, 0.0, 500);
        TH2 *hHmComp = new TH2F("H_m_comp", "m_{h}", 20, 0.0, 500, 20, 0.0, 500);


// jets
    // reco 
        // pT + m
        TH1 *hjjv1pTreco = new TH1F("jj_pT_reco_v1", "p^{T}_{jj}", 100, 0, 1000);
        TH1 *hjjv1mreco = new TH1F("jj_m_reco_v1", "m_{jj}", 100, 0, 500);
        TH1 *hjjv2pTreco = new TH1F("jj_pT_reco_v2", "p^{T}_{jj}", 100, 0, 1000);
        TH1 *hjjv2mreco = new TH1F("jj_m_reco_v2", "m_{jj}", 100, 0, 500);
        TH1 *hj1pTreco = new TH1F("j1_pT_reco", "p^{T}_{j1}", 100, 0, 1000);
        TH1 *hj1mreco = new TH1F("j1_m_reco", "m_{j1}", 100, 0, 100);
        TH1 *hj2pTreco = new TH1F("j2_pT_reco", "p^{T}_{j2}", 100, 0, 1000);
        TH1 *hj2mreco = new TH1F("j2_m_reco", "m_{j2}", 100, 0, 100);
        // eta
        TH1 *hj1etareco = new TH1F("j1_eta_reco", "eta_{j1}", 100, -5, 5);
        TH1 *hj2etareco = new TH1F("j2_eta_reco", "eta_{j2}", 100, -5, 5);
        // phi
        TH1F *hjjdeltaPhireco = new TH1F("jj_#phi_reco", "#phi_{jj}",10,-TMath::Pi(),+TMath::Pi());

// z 
    // reco
        // pT + m
        TH1 *hz1pTreco = new TH1F("z1_pT_reco", "p^{T}_{z1}", 100, 0, 1000);
        TH1 *hz2pTreco = new TH1F("z2_pT_reco", "p^{T}_{z2}", 100, 0, 1000);
        TH1 *hz1mreco = new TH1F("z1_m_reco", "m_{z1}", 100, 0, 250);
        TH1 *hz2mreco = new TH1F("z2_m_reco", "m_{z2}", 100, 0, 250);
        // cos
        TH1 *hz1cosThetareco = new TH1F("z1_cos#theta_reco", "cos#theta_{z1}",10,-1,1);
        TH1 *hz2cosThetareco = new TH1F("z2_cos#theta_reco", "cos#theta_{z1}",10,-1,1);
        // phi
        TH1 *hzzdeltaPhireco = new TH1F("zz_#phi_reco", "#phi_{zz}",10,-TMath::Pi(),+TMath::Pi());
    // particle
        // pT + m
        TH1 *hz1pTparticle = new TH1F("z1_pT_particle", "p^{T}_{z1}", 100, 0, 1000);
        TH1 *hz2pTparticle = new TH1F("z2_pT_particle", "p^{T}_{z2}", 100, 0, 1000);
        TH1 *hz1mparticle = new TH1F("z1_m_particle", "m_{z1}", 100, 0, 250);
        TH1 *hz2mparticle = new TH1F("z2_m_particle", "m_{z2}", 100, 0, 250);
    // parton
        // pT + m
        TH1 *hz1pTparton = new TH1F("z1_pT_parton", "p^{T}_{z1}", 100, 0, 1000);
        TH1 *hz2pTparton = new TH1F("z2_pT_parton", "p^{T}_{z2}", 100, 0, 1000);
        TH1 *hz1mparton = new TH1F("z1_m_parton", "m_{z1}", 100, 0, 250);
        TH1 *hz2mparton = new TH1F("z2_m_parton", "m_{z2}", 100, 0, 250);
    //2D
        TH2 *hz1pTComp = new TH2F("z1_pT_comp", "p^{T}_{z1}", 20, 0.0, 500, 20, 0.0, 500);
        TH2 *hz1mComp = new TH2F("z1_m_comp", "m_{z1}", 20, 0.0, 500, 20, 0.0, 500);
        TH2 *hz2pTComp = new TH2F("z2_pT_comp", "p^{T}_{z2}", 20, 0.0, 500, 20, 0.0, 500);
        TH2 *hz2mComp = new TH2F("z2_m_comp", "m_{z2}", 20, 0.0, 500, 20, 0.0, 500);

// leptons
    // reco
        // pT + m
        TH1 *hl1pTreco = new TH1F("l1_pT_reco", "p^{T}_{l1}", 100, 0, 1000);
        TH1 *hl2pTreco = new TH1F("l2_pT_reco", "p^{T}_{l2}", 100, 0, 1000);
        TH1 *hl3pTreco = new TH1F("l3_pT_reco", "p^{T}_{l3}", 100, 0, 1000);
        TH1 *hl4pTreco = new TH1F("l4_pT_reco", "p^{T}_{l4}", 100, 0, 1000);
        // phi
        TH1 *hl1l2deltaPhireco = new TH1F("l1l2_#phi_reco", "#phi_{l1l2}",10,-TMath::Pi(),+TMath::Pi());
        TH1 *hl3l4deltaPhireco = new TH1F("l3l4_#phi_reco", "#phi_{l3l4}",10,-TMath::Pi(),+TMath::Pi());
  
  double  nPassed=0;
  double Lumi=3e3;
  double totWeightedEntries=0;
  int  nPassedRaw=0;
  int bjets = 0;
  int nQuads=0;
  bool electronEvent;

  GenParticle *daughter1;
  GenParticle *daughter2;

  GenParticle *particle1;
  GenParticle *particle2;

  GenParticle *p1daughter1;
  GenParticle *p1daughter2;
  GenParticle *p2daughter1;
  GenParticle *p2daughter2;

  TLorentzVector j1_reco, j1_particle,  j1_parton;
  TLorentzVector j2_reco, j2_particle,  j2_parton;

  TLorentzVector b1_reco, b1_particle,  b1_parton;
  TLorentzVector b2_reco, b2_particle,  b2_parton;

  TLorentzVector h_reco, h_parton, h_particle;

  TLorentzVector e1_reco, e1_particle,  e1_parton;
  TLorentzVector e2_reco, e2_particle,  e2_parton;
  TLorentzVector e3_reco, e3_particle,  e3_parton;
  TLorentzVector e4_reco, e4_particle,  e4_parton;

  TLorentzVector m1_reco, m1_particle,  m1_parton;
  TLorentzVector m2_reco, m2_particle,  m2_parton;
  TLorentzVector m3_reco, m3_particle,  m3_parton;
  TLorentzVector m4_reco, m4_particle,  m4_parton;

  TLorentzVector l1_reco, l1_particle,  l1_parton;
  TLorentzVector l2_reco, l2_particle,  l2_parton;
  TLorentzVector l3_reco, l3_particle,  l3_parton;
  TLorentzVector l4_reco, l4_particle,  l4_parton;

  TLorentzVector z1_reco, z1_particle,  z1_parton;
  TLorentzVector z2_reco, z2_particle,  z2_parton;

  TLorentzVector fourl_reco, fourl_particle, fourl_parton; 

//------------------------------------------------------------------------------------------------------------------------------------------------------------
// reco
//------------------------------------------------------------------------------------------------------------------------------------------------------------

  for(Int_t entry = 0; entry < numberOfEntries; ++entry) {

    treeReader->ReadEntry(entry);
    HepMCEvent *event = (HepMCEvent*) branchEvent -> At(0);
    
    Float_t weight = event->Weight/numberOfEntries*Lumi;
    totWeightedEntries+=weight;

// higgs

    vector <int> btagIndex;
    vector <int> noBtag;
    vector <int> goodJetIndex;

    for(int i=0; i<(int)branchJet->GetEntries(); i++){
      Jet *jet=(Jet*) branchJet->At(i);
    //if( jet->PT < 20) continue;
    //if (fabs(jet->Eta) > 4.4) continue; 
      if( jet->BTag>0) btagIndex.push_back(i);
      else noBtag.push_back(i);
           goodJetIndex.push_back(i);
    }
    
    sort(btagIndex.begin(), btagIndex.end(), [branchJet](const int& lhs, const int& rhs) {
	    return ((Jet*)branchJet->At(lhs))->PT < ((Jet*)branchJet->At(rhs))->PT;
    });
    sort(noBtag.begin(), noBtag.end(), [branchJet](const int& lhs, const int& rhs) {
	    return ((Jet*)branchJet->At(lhs))->PT < ((Jet*)branchJet->At(rhs))->PT;
    });
    sort(goodJetIndex.begin(), goodJetIndex.end(), [branchJet](const int& lhs, const int& rhs) {
	    return ((Jet*)branchJet->At(lhs))->PT < ((Jet*)branchJet->At(rhs))->PT;
    });

    if(btagIndex.size() <1) continue ; // at least one b tag 
    if(goodJetIndex.size() < 2) continue ; // at least two jets 

    Jet *b1=nullptr;
    Jet *b2=nullptr;

    vector<pair<int,int>> bJetPairs;
    vector<vector <int>> bJetPairsComb=combinationsNoRepetitionAndOrderDoesNotMatter(2,goodJetIndex);
//  vector<vector <int>> bJetPairsComb=combinationsNoRepetitionAndOrderDoesNotMatter(2,btagIndex);

    if( bJetPairsComb.size() < 1) continue; // need at least two good jets; 
    for(int i=0; i<(int)bJetPairsComb.size(); i++)
      bJetPairs.push_back(make_pair(bJetPairsComb[i][0],bJetPairsComb[i][1]));


    if( bJetPairs.size() > 1) 
      sort(bJetPairs.begin(), bJetPairs.end(), [branchJet](const pair<int,int> lhs, const pair<int,int> rhs) {
	  return fabs(((((Jet*)branchJet->At(lhs.first))->P4() + ((Jet*)branchJet->At(lhs.second))->P4())).M() - 125 ) <
	    fabs( ((((Jet*)branchJet->At(rhs.first))->P4() + ((Jet*)branchJet->At(rhs.second))->P4()).M()) - 125 ) ; 
	});

    
    pair <int,int> higgsbbcandidate;
    bool foundBjet=false; 
    for(int i=0; i<(int) bJetPairs.size(); i++){
      b1=(Jet*)branchJet->At(bJetPairs[i].first);
      b2=(Jet*)branchJet->At(bJetPairs[i].second);
      if( b1->BTag>0 || b2->BTag>0) {
    higgsbbcandidate=bJetPairs[i];
    foundBjet=true;
	  break;
      }
    }

    h_reco = b1->P4() + b2->P4();

    if(!foundBjet) continue; 
    
// jets

    vector <int> nonHiggsJet;

    for(int i=0; i<(int)goodJetIndex.size(); i++){
      if( goodJetIndex[i] == higgsbbcandidate.first  || goodJetIndex[i] == higgsbbcandidate.second) continue;
      nonHiggsJet.push_back(i);
    }

    if(nonHiggsJet.size() < 2) continue ; 
    vector<pair<int,int>> vbfJetIndex;
    vector<vector <int>> vbfJetIndexComb=combinationsNoRepetitionAndOrderDoesNotMatter(2,nonHiggsJet);
    if( vbfJetIndexComb.size() < 1 ) continue; 
    for(int i=0; i<(int)vbfJetIndexComb.size(); i++)
      vbfJetIndex.push_back(make_pair(vbfJetIndexComb[i][0],vbfJetIndexComb[i][1]));
    if( branchMuon->GetEntries() + branchElectron->GetEntries() < 4) continue;

    Jet *jet1 = (Jet*) branchJet->At(vbfJetIndex[0].first);
    Jet *jet2 = (Jet*) branchJet->At(vbfJetIndex[0].second);

    j1_reco=jet1->P4();
    j2_reco=jet2->P4();

    double jjdeltaPhireco=(j1_reco.Phi() > j2_reco.Phi() ? -1:+1)*TMath::Abs(j1_reco.Phi() - j2_reco.Phi());

// leptons + z

    vector <int> goodE_reco_indices; 
    for(int i=0; i<(int)branchElectron->GetEntries(); i++){
      Electron *el_reco = (Electron *) branchElectron->At(i);

      // pT and eta cuts 
      if( el_reco->PT > 1 && fabs(el_reco->Eta) < 2.5) goodE_reco_indices.push_back(i);
    }

    // sort the indices by pT ;
    sort(goodE_reco_indices.begin(), goodE_reco_indices.end(), [branchElectron](const int& lhs, const int& rhs) {
	return ((Electron*)branchElectron->At(lhs))->PT < ((Electron*)branchElectron->At(rhs))->PT;
      });
    
    vector <int> goodMu_reco_indices; 
    for(int i=0; i<(int)branchMuon->GetEntries(); i++){
    // pT and eta cuts 
      Muon *mu_reco = (Muon *) branchMuon->At(i);
      if( mu_reco->PT > 1 && fabs(mu_reco->Eta) < 2.5) goodMu_reco_indices.push_back(i);
    }

    // sort the indices by pT ;
    sort(goodMu_reco_indices.begin(), goodMu_reco_indices.end(), [branchMuon](const int& lhs, const int& rhs) {
	return ((Muon*)branchMuon->At(lhs))->PT < ((Muon*)branchMuon->At(rhs))->PT;
      });
    
    // form pairs for each flavour

    // electrons 
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

    // muons 
    vector< pair<int,int>> muRecoPairIndicesIn;
    vector< pair<int,int>> muRecoPairIndices;
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

    sort(muRecoPairIndices.begin(),muRecoPairIndices.end(), [branchMuon]( pair<int,int>   & lhs,  pair<int,int>   & rhs) {
	    int index11_reco=(lhs).first;
	    int index12_reco=(lhs).second;
	    int index21_reco=(rhs).first;
    	int index22_reco=(rhs).second;
	    return fabs(((((Muon*)branchMuon->At(index11_reco))->P4() + ((Muon*)branchMuon->At(index12_reco))->P4())).M() - 91) <
	  fabs( ((((Muon*)branchMuon->At(index21_reco))->P4() + ((Muon*)branchMuon->At(index22_reco))->P4()).M()) -91);
      });

    remove_overlaps(muRecoPairIndices);
    
    // order quads by mZ1 and mZ2
    if( muRecoPairIndices.size() + elecRecoPairIndices.size() < 2) continue; // no candidate found;
   
    int thisRecoEventType=-1; 

    vector<pair<int,pair<int,int>>> RecoPairIndices; // 0 for electron 1 for muon
    for(int i=0; i<(int) elecRecoPairIndices.size(); i++)
      RecoPairIndices.push_back(make_pair(0,elecRecoPairIndices.at(i)));
    for(int i=0; i<(int) muRecoPairIndices.size(); i++)
      RecoPairIndices.push_back(make_pair(1,muRecoPairIndices.at(i)));

    // sort all of the indices by closeness to mZ
    sort(RecoPairIndices.begin(), RecoPairIndices.end(), [branchMuon,branchElectron]  ( const pair<int,pair<int,int>>  lhs , const pair<int,pair<int,int>>   rhs ){

	    double mass1_reco=0 ;
    	double mass2_reco=0;
	
      if( lhs.first==0 ) mass1_reco = (((Electron*)branchElectron->At(lhs.second.first))->P4() + ((Electron*)branchElectron->At(lhs.second.second))->P4()).M();
      else mass1_reco = ((((Muon*)branchMuon->At(lhs.second.first))->P4() + ((Muon*)branchMuon->At(lhs.second.second))->P4())).M();
          
      if( rhs.first==0 ) mass2_reco = (((Electron*)branchElectron->At(rhs.second.first))->P4() + ((Electron*)branchElectron->At(rhs.second.second))->P4()).M();
      else mass2_reco = (((Muon*)branchMuon->At(rhs.second.first))->P4() + ((Muon*)branchMuon->At(rhs.second.second))->P4()).M();
          
	    return fabs(mass1_reco ) <  fabs(mass2_reco) ;
    });
  
  if( RecoPairIndices[0].first == 1 && RecoPairIndices[1].first == 1) thisRecoEventType=0;
  else if( RecoPairIndices[0].first == 0 && RecoPairIndices[1].first == 0) thisRecoEventType=1;
  else if( RecoPairIndices[0].first == 1 && RecoPairIndices[1].first == 0) thisRecoEventType=2;
  else if( RecoPairIndices[0].first == 0 && RecoPairIndices[1].first == 1) thisRecoEventType=3;
  
  // case 4mu 
  if( thisRecoEventType==0 ) {
    
    // take first two muons
    Muon *muon1_reco = (Muon *) branchMuon->At( RecoPairIndices[0].second.first);
    Muon *muon2_reco = (Muon *) branchMuon->At( RecoPairIndices[0].second.second);
    // take first two muons
    Muon *muon3_reco = (Muon *) branchMuon->At( RecoPairIndices[1].second.first);
    Muon *muon4_reco = (Muon *) branchMuon->At( RecoPairIndices[1].second.second);

    fourl_reco=muon1_reco->P4() + muon2_reco->P4() + muon3_reco->P4() + muon4_reco->P4(); 
    
    z1_reco=muon1_reco->P4() + muon2_reco->P4() ;
    z2_reco=muon3_reco->P4() + muon4_reco->P4() ;

    l1_reco=muon1_reco->P4();
    l2_reco=muon2_reco->P4();
    l3_reco=muon3_reco->P4();
    l4_reco=muon4_reco->P4(); 

  }

  // case 4e
  else if( thisRecoEventType==1 ) {
    
    // take first two electrons
    Electron *muon1_reco = (Electron *) branchElectron->At( RecoPairIndices[0].second.first);
    Electron *muon2_reco = (Electron *) branchElectron->At( RecoPairIndices[0].second.second);
    // take first two electrons
    Electron *muon3_reco = (Electron *) branchElectron->At( RecoPairIndices[1].second.first);
    Electron *muon4_reco = (Electron *) branchElectron->At( RecoPairIndices[1].second.second);

    fourl_reco=muon1_reco->P4() + muon2_reco->P4() + muon3_reco->P4() + muon4_reco->P4(); 

    z1_reco=muon1_reco->P4() + muon2_reco->P4() ;
    z2_reco=muon3_reco->P4() + muon4_reco->P4() ;

    l1_reco=muon1_reco->P4();
    l2_reco=muon2_reco->P4();
    l3_reco=muon3_reco->P4();
    l4_reco=muon4_reco->P4(); 
  }

  // case 2mu2e 
  if( thisRecoEventType==2 ) {
    
    // take first two muons
    Muon *muon1_reco = (Muon *) branchMuon->At( RecoPairIndices[0].second.first);
    Muon *muon2_reco = (Muon *) branchMuon->At( RecoPairIndices[0].second.second);
    // take first two electrons
    Electron *muon3_reco = (Electron *) branchElectron->At( RecoPairIndices[1].second.first);
    Electron *muon4_reco = (Electron *) branchElectron->At( RecoPairIndices[1].second.second);

    fourl_reco=muon1_reco->P4() + muon2_reco->P4() + muon3_reco->P4() + muon4_reco->P4(); 

    z1_reco=muon1_reco->P4() + muon2_reco->P4() ;
    z2_reco=muon3_reco->P4() + muon4_reco->P4() ;

    l1_reco=muon1_reco->P4();
    l2_reco=muon2_reco->P4();
    l3_reco=muon3_reco->P4();
    l4_reco=muon4_reco->P4(); 
  }

  // case 2e2mu 
   else if( thisRecoEventType==3 ) {
    
    // take first two electrons
    Electron *muon1_reco = (Electron *) branchElectron->At( RecoPairIndices[0].second.first);
    Electron *muon2_reco = (Electron *) branchElectron->At( RecoPairIndices[0].second.second);
    // take first two muons
    Muon *muon3_reco = (Muon *) branchMuon->At( RecoPairIndices[1].second.first);
    Muon *muon4_reco = (Muon *) branchMuon->At( RecoPairIndices[1].second.second);

    fourl_reco=muon1_reco->P4() + muon2_reco->P4() + muon3_reco->P4() + muon4_reco->P4(); 

    z1_reco=muon1_reco->P4() + muon2_reco->P4() ;
    z2_reco=muon3_reco->P4() + muon4_reco->P4() ;

    l1_reco=muon1_reco->P4();
    l2_reco=muon2_reco->P4();
    l3_reco=muon3_reco->P4();
    l4_reco=muon4_reco->P4(); 

   }

    double zzdeltaPhireco=(z2_reco.Phi() > z1_reco.Phi() ? -1:+1)*TMath::Abs(z1_reco.Phi() - z2_reco.Phi());
    double l1l2deltaPhireco=(l2_reco.Phi() > l1_reco.Phi() ? -1:+1)*TMath::Abs(l1_reco.Phi() - l2_reco.Phi());
    double l3l4deltaPhireco=(l4_reco.Phi() > l3_reco.Phi() ? -1:+1)*TMath::Abs(l3_reco.Phi() - l4_reco.Phi());

//------------------------------------------------------------------------------------------------------------------------------------------------------------
// particle
//------------------------------------------------------------------------------------------------------------------------------------------------------------

// higgs

    for(int i=0; i<(int)branchGenJet->GetEntries(); i++){
      Jet *genjet=(Jet*) branchGenJet->At(i);
      if (genjet -> Flavor == 5) bjets += 1; 
      else continue;
      if (bjets == 1) {
        b1_particle = genjet->P4();
      } else if (bjets == 2){
        b2_particle = genjet->P4();
        break;
      }
    }

    h_particle = b1_particle + b2_particle;

// z + leptons

// electrons 
vector <int> goodE_particle_indices; 
for(int i=0; i<(int)branchGenParticle->GetEntries(); i++){
    GenParticle *particle=(GenParticle*) branchGenParticle->At(i);  
    if( particle->Status !=1 ) continue;
    if( fabs(particle->PID) == 11 && fabs(particle->Eta)< 2.5 && particle->PT > 2.5) 
            goodE_particle_indices.push_back(i); 
}

// muons  
vector <int> goodMu_particle_indices; 
for(int i=0; i<(int)branchGenParticle->GetEntries(); i++){
    GenParticle *particle=(GenParticle*) branchGenParticle->At(i);  
    if( particle->Status !=1 ) continue; 
    // electrons 
    if( fabs(particle->PID) == 13 && fabs(particle->Eta)< 2.5 && particle->PT > 2.5) 
            goodMu_particle_indices.push_back(i); 
// neutrinos 
// jets we use the GenJet container.. 
}

    // sort the indices by pT ;
    sort(goodE_particle_indices.begin(), goodE_particle_indices.end(), [branchGenParticle](const int& lhs, const int& rhs) {
	return ((GenParticle*) branchGenParticle->At(lhs))->PT < ((GenParticle*) branchGenParticle->At(rhs))->PT;
      });

    // sort the indices by pT ;
    sort(goodMu_particle_indices.begin(), goodMu_particle_indices.end(), [branchGenParticle](const int& lhs, const int& rhs) {
		return ((GenParticle*) branchGenParticle->At(lhs))->PT < ((GenParticle*) branchGenParticle->At(rhs))->PT;
      });
    
    // form pairs 
    // electrons 
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
    	return fabs(((((GenParticle*) branchGenParticle->At(lhs.first))->P4() + ((GenParticle*) branchGenParticle->At(lhs.second))->P4())).M() -91 ) <
	  fabs( ((((GenParticle*) branchGenParticle->At(rhs.first))->P4() + ((GenParticle*) branchGenParticle->At(rhs.second))->P4()).M()) -91 ) ; 
      });

    remove_overlaps(elecParticlePairIndices);

    // muons 
    vector< pair<int,int>> muParticlePairIndicesIn;
    vector< pair<int,int>> muParticlePairIndices;
    vector <vector<int>> muParticlePairIndices_;
    if(goodMu_particle_indices.size() > 1) 
      muParticlePairIndices_=combinationsNoRepetitionAndOrderDoesNotMatter(2,goodMu_particle_indices);

    // remove all combinations not satifing criteria 
    for(int i=0; i<(int)muParticlePairIndices_.size(); i++){
      int elecParticleIndex=muParticlePairIndices_[i].at(0);
      int elecParticleIndex2=muParticlePairIndices_[i].at(1);

      GenParticle *el1_particle=(GenParticle*) branchGenParticle->At(elecParticleIndex);
      GenParticle *el2_particle=(GenParticle*) branchGenParticle->At(elecParticleIndex2);
      if( el1_particle->Charge ==  el2_particle->Charge ) continue;
      muParticlePairIndicesIn.push_back(make_pair(elecParticleIndex,elecParticleIndex2));
    }

    muParticlePairIndices=muParticlePairIndicesIn;

    sort(muParticlePairIndices.begin(),muParticlePairIndices.end(), [branchGenParticle]( pair<int,int>   & lhs,  pair<int,int>   & rhs) {
	    int index11_particle=(lhs).first;
	    int index12_particle=(lhs).second;
	    int index21_particle=(rhs).first;
    	int index22_particle=(rhs).second;
	    return fabs(((((GenParticle*) branchGenParticle->At(index11_particle))->P4() + ((GenParticle*) branchGenParticle->At(index12_particle))->P4())).M() - 91) <
	  fabs( ((((GenParticle*) branchGenParticle->At(index21_particle))->P4() + ((GenParticle*) branchGenParticle->At(index22_particle))->P4()).M()) -91);
      });

    remove_overlaps(muParticlePairIndices);
    
    // order quads by mZ1 and mZ2
    if( muParticlePairIndices.size() + elecParticlePairIndices.size() < 2) continue; // no candidate found;
   
    int thisParticleEventType=-1; 

    vector<pair<int,pair<int,int>>> ParticlePairIndices; // 0 for electron 1 for muon
    for(int i=0; i<(int) elecParticlePairIndices.size(); i++)
      ParticlePairIndices.push_back(make_pair(0,elecParticlePairIndices.at(i)));
    for(int i=0; i<(int) muParticlePairIndices.size(); i++)
      ParticlePairIndices.push_back(make_pair(1,muParticlePairIndices.at(i)));

    // sort all of the indices by closeness to mZ
    sort(ParticlePairIndices.begin(), ParticlePairIndices.end(), [branchGenParticle]  ( const pair<int,pair<int,int>>  lhs , const pair<int,pair<int,int>>   rhs ){

	    double mass1_particle=0 ;
    	double mass2_particle=0;
	
      if( lhs.first==0 ) mass1_particle = (((GenParticle*) branchGenParticle->At(lhs.second.first))->P4() + ((GenParticle*) branchGenParticle->At(lhs.second.second))->P4()).M();
      else mass1_particle= ((((GenParticle*) branchGenParticle->At(lhs.second.first))->P4() + ((GenParticle*) branchGenParticle->At(lhs.second.second))->P4())).M();
          
      if( rhs.first==0 ) mass2_particle = (((GenParticle*) branchGenParticle->At(rhs.second.first))->P4() + ((GenParticle*) branchGenParticle->At(rhs.second.second))->P4()).M();
      else mass2_particle = (((GenParticle*) branchGenParticle->At(rhs.second.first))->P4() + ((GenParticle*) branchGenParticle->At(rhs.second.second))->P4()).M();
          
	    return fabs(mass1_particle ) <  fabs(mass2_particle) ;
    });
  
  if( ParticlePairIndices[0].first == 1 && ParticlePairIndices[1].first == 1) thisParticleEventType=0;
  else if( ParticlePairIndices[0].first == 0 && ParticlePairIndices[1].first == 0) thisParticleEventType=1;
  else if( ParticlePairIndices[0].first == 1 && ParticlePairIndices[1].first == 0) thisParticleEventType=2;
  else if( ParticlePairIndices[0].first == 0 && ParticlePairIndices[1].first == 1) thisParticleEventType=3;
  
  // case 4mu 
  if( thisParticleEventType==0 ) {
    
    // take first two muons
    GenParticle *muon1_particle = (GenParticle*) branchGenParticle->At( ParticlePairIndices[0].second.first);
    GenParticle *muon2_particle = (GenParticle*) branchGenParticle->At( ParticlePairIndices[0].second.second);
    // take first two muons
    GenParticle *muon3_particle = (GenParticle*) branchGenParticle->At( ParticlePairIndices[1].second.first);
    GenParticle *muon4_particle = (GenParticle*) branchGenParticle->At( ParticlePairIndices[1].second.second);

    fourl_particle=muon1_particle->P4() + muon2_particle->P4() + muon3_particle->P4() + muon4_particle->P4(); 
    
    z1_particle=muon1_particle->P4() + muon2_particle->P4() ;
    z2_particle=muon3_particle->P4() + muon4_particle->P4() ;

    l1_particle=muon1_particle->P4();
    l2_particle=muon2_particle->P4();
    l3_particle=muon3_particle->P4();
    l4_particle=muon4_particle->P4(); 

  }

  // case 4e
  else if( thisParticleEventType==1 ) {
    
    // take first two electrons
    GenParticle *muon1_particle = (GenParticle*) branchGenParticle->At( ParticlePairIndices[0].second.first);
    GenParticle *muon2_particle = (GenParticle*) branchGenParticle->At( ParticlePairIndices[0].second.second);
    // take first two electrons
    GenParticle *muon3_particle = (GenParticle*) branchGenParticle->At( ParticlePairIndices[1].second.first);
    GenParticle *muon4_particle = (GenParticle*) branchGenParticle->At( ParticlePairIndices[1].second.second);

    fourl_particle=muon1_particle->P4() + muon2_particle->P4() + muon3_particle->P4() + muon4_particle->P4(); 

    z1_particle=muon1_particle->P4() + muon2_particle->P4() ;
    z2_particle=muon3_particle->P4() + muon4_particle->P4() ;

    l1_particle=muon1_particle->P4();
    l2_particle=muon2_particle->P4();
    l3_particle=muon3_particle->P4();
    l4_particle=muon4_particle->P4(); 
  }

  // case 2mu2e 
  if( thisParticleEventType==2 ) {
    
    // take first two muons
    GenParticle *muon1_particle = (GenParticle*) branchGenParticle->At( ParticlePairIndices[0].second.first);
    GenParticle *muon2_particle = (GenParticle*) branchGenParticle->At( ParticlePairIndices[0].second.second);
    // take first two electrons
    GenParticle *muon3_particle = (GenParticle*) branchGenParticle->At( ParticlePairIndices[1].second.first);
    GenParticle *muon4_particle = (GenParticle*) branchGenParticle->At( ParticlePairIndices[1].second.second);

    fourl_particle=muon1_particle->P4() + muon2_particle->P4() + muon3_particle->P4() + muon4_particle->P4(); 

    z1_particle=muon1_particle->P4() + muon2_particle->P4() ;
    z2_particle=muon3_particle->P4() + muon4_particle->P4() ;

    l1_particle=muon1_particle->P4();
    l2_particle=muon2_particle->P4();
    l3_particle=muon3_particle->P4();
    l4_particle=muon4_particle->P4(); 
  }

  // case 2e2mu 
   else if( thisParticleEventType==3 ) {
    
    // take first two electrons
    GenParticle *muon1_particle = (GenParticle*) branchGenParticle->At( ParticlePairIndices[0].second.first);
    GenParticle *muon2_particle = (GenParticle*) branchGenParticle->At( ParticlePairIndices[0].second.second);
    // take first two muons
    GenParticle *muon3_particle = (GenParticle*) branchGenParticle->At( ParticlePairIndices[1].second.first);
    GenParticle *muon4_particle = (GenParticle*) branchGenParticle->At( ParticlePairIndices[1].second.second);

    fourl_particle=muon1_particle->P4() + muon2_particle->P4() + muon3_particle->P4() + muon4_particle->P4(); 

    z1_particle=muon1_particle->P4() + muon2_particle->P4() ;
    z2_particle=muon3_particle->P4() + muon4_particle->P4() ;

    l1_particle=muon1_particle->P4();
    l2_particle=muon2_particle->P4();
    l3_particle=muon3_particle->P4();
    l4_particle=muon4_particle->P4(); 

   }

//------------------------------------------------------------------------------------------------------------------------------------------------------------
// parton
//------------------------------------------------------------------------------------------------------------------------------------------------------------

// higgs

    for(int i=0; i<(int)branchGenParticle->GetEntries(); i++){
      GenParticle *particle=(GenParticle*) branchGenParticle->At(i);
      int d1_pid = 9999;
      if (particle->D1 != -1) {
        daughter1 = (GenParticle*) branchGenParticle->At(particle->D1);
        d1_pid = daughter1 -> PID;
      }
      int d2_pid = 9999;
      if (particle->D2 != -1) {
        daughter2 = (GenParticle*) branchGenParticle->At(particle->D2);
        d2_pid = daughter2 -> PID;
      } 
      // higgs parton
      if (particle->PID == 25 && d1_pid != 25) {
        h_parton = particle->P4();
        // check for b parton children
        if (abs(d1_pid) == 5) {
          if (daughter1 -> PT > daughter2 -> PT) {
            b1_parton = daughter1 -> P4();
            b2_parton = daughter2 -> P4();
          } else {
            b1_parton = daughter2 -> P4();
            b2_parton = daughter1 -> P4();
          }
        }
      }
    }

// z + leptons

vector <int> genZBosons; 
for(int i=0; i<(int)branchGenParticle->GetEntries(); i++){
    GenParticle *particle=(GenParticle*)branchGenParticle->At(i);
    int d1_pid = 9999;
    if (particle->D1 != -1) {
        daughter1 = (GenParticle*) branchGenParticle->At(particle->D1);
        d1_pid = daughter1 -> PID;
      }
    // find a Z.
        if( particle->PID == 23 && d1_pid != 23) 
        genZBosons.push_back(i);
}

    sort(genZBosons.begin(), genZBosons.end(), [branchGenParticle](const int& lhs, const int& rhs) {
	return ((GenParticle*)branchGenParticle->At(lhs))->P4().M() < (((GenParticle*)branchGenParticle->At(rhs))->P4()).M(); 
    });
/**
    vector< pair<int,int>> ZPairIndices;
    vector< pair<int,int>> ZPairIndicesIn;
    vector <vector<int>> ZPairIndices_;
    if( genZBosons.size() > 1 )
    ZPairIndices_=combinationsNoRepetitionAndOrderDoesNotMatter(2,genZBosons);

    for(int i=0;i<(int)ZPairIndices_.size(); i++){
      int Z1Index=ZPairIndices_[i].at(0);
      int Z2Index=ZPairIndices_[i].at(1);
      
      GenParticle *particle1=(GenParticle*)branchGenParticle->At(Z1Index);
      GenParticle *particle2=(GenParticle*)branchGenParticle->At(Z2Index);
      
      ZPairIndicesIn.push_back(make_pair(Z1Index,Z2Index));
    }

    ZPairIndices=ZPairIndicesIn;
*/
/*
    sort(ZPairIndices.begin(), ZPairIndices.end(), [branchGenParticle](const pair<int,int> lhs, const pair<int,int> rhs) {
    	return fabs(((((GenParticle*)branchGenParticle->At(lhs.first))->P4() + ((GenParticle*)branchGenParticle->At(lhs.second))->P4())).M() -182 ) <
	  fabs( ((((GenParticle*)branchGenParticle->At(rhs.first))->P4() + ((GenParticle*)branchGenParticle->At(rhs.second))->P4()).M()) -182 ) ; 
      });
*/
/*
     // sort all of the indices by closeness to mZ
    sort(ZPairIndices.begin(), ZPairIndices.end(), [branchGenParticle]  ( const pair<int,pair<int,int>>  lhs , const pair<int,pair<int,int>>   rhs ){

	    double z1mass=0 ;
    	double z2mass=0;
	
      if( lhs.first==0 ) z1mass = (((GenParticle*)branchGenParticle->At(lhs.first))->P4()).M();
      else z1mass = (((GenParticle*)branchGenParticle->At(lhs.second))->P4()).M();;
          
      if( rhs.first==0 ) z2mass = (((GenParticle*)branchGenParticle->At(rhs.first))->P4()).M();
      else z2mass = (((GenParticle*)branchGenParticle->At(rhs.first))->P4()).M();
          
	    return fabs(z1mass) < fabs(z2mass) ;
    });

    remove_overlaps(ZPairIndices);
*/

//for(int i=0; i<(int)branchGenParticle->GetEntries(); i++){
  //  GenParticle *particle=(GenParticle*)branchGenParticle->At(i);

z1_parton = ((GenParticle*)branchGenParticle->At(genZBosons[0]))->P4(); 
z2_parton = ((GenParticle*)branchGenParticle->At(genZBosons[1]))->P4(); 

vector <int> Z1children;
Z1children.push_back( ((GenParticle*)branchGenParticle->At(genZBosons[0]))->D1 ); 
Z1children.push_back( ((GenParticle*)branchGenParticle->At(genZBosons[0]))->D2 ); 

vector <int> Z2children; 
Z2children.push_back( ((GenParticle*)branchGenParticle->At(genZBosons[1]))->D1 );
Z2children.push_back( ((GenParticle*)branchGenParticle->At(genZBosons[1]))->D2 ); 

p1daughter1 = (GenParticle*) branchGenParticle->At(Z1children[0]);
p1daughter2 = (GenParticle*) branchGenParticle->At(Z1children[1]);
/*
p1_d1_pid = p1daughter1 -> PID;
p1_d2_pid = p1daughter2 -> PID;
*/
p2daughter1 = (GenParticle*) branchGenParticle->At(Z2children[0]);
p2daughter2 = (GenParticle*) branchGenParticle->At(Z2children[1]);
/*
p2_d1_pid = p2daughter1 -> PID;
p2_d2_pid = p2daughter2 -> PID;
*/

/*
      int p1_d1_pid = 9999;
      if (particle1->D1 != -1) {
        p1daughter1 = (GenParticle*) branchGenParticle->At(particle1->D1);
        p1_d1_pid = p1daughter1 -> PID;
      }
      int p1_d2_pid = 9999;
      if (particle1->D2 != -1) {
        p1daughter2 = (GenParticle*) branchGenParticle->At(particle1->D2);
        p1_d2_pid = p2daughter2 -> PID;
      }
      int p2_d1_pid = 9999;
      if (particle2->D1 != -1) {
        p2daughter1 = (GenParticle*) branchGenParticle->At(particle2->D1);
        p2_d1_pid = p2daughter1 -> PID;
      }
      int p2_d2_pid = 9999;
      if (particle2->D2 != -1) {
        p2daughter2 = (GenParticle*) branchGenParticle->At(particle2->D2);
        p2_d2_pid = p2daughter2 -> PID;
      } 
*/

/*
        if ((abs(p1_d1_pid) == 11) && (abs(p1_d2_pid) == 11)) {
          if (p1daughter1 -> PT > p1daughter2 -> PT) {
            e1_parton = p1daughter1 -> P4();
            e2_parton = p1daughter2 -> P4();
            e1_particle = find_status1_child(branchGenParticle, p1daughter1, 11) -> P4();
            e2_particle = find_status1_child(branchGenParticle, p1daughter2, 11) -> P4();
          } else {
            e1_parton = p1daughter2 -> P4();
            e2_parton = p1daughter1 -> P4();
            e1_particle = find_status1_child(branchGenParticle, p1daughter2, 11) -> P4();
            e2_particle = find_status1_child(branchGenParticle, p1daughter1, 11) -> P4();
          }
        } else if (abs(p2_d1_pid) == 11) {
          if (p2daughter1 -> PT > p2daughter2 -> PT) {
            e1_parton = p2daughter1 -> P4();
            e2_parton = p2daughter2 -> P4();
            e1_particle = find_status1_child(branchGenParticle, p2daughter1, 11) -> P4();
            e2_particle = find_status1_child(branchGenParticle, p2daughter2, 11) -> P4();
          } else {
            e1_parton = p2daughter2 -> P4();
            e2_parton = p2daughter1 -> P4();
            e1_particle = find_status1_child(branchGenParticle, p2daughter2, 11) -> P4();
            e2_particle = find_status1_child(branchGenParticle, p2daughter1, 11) -> P4();
          }
        // check for muon daughters
        } else if (abs(p1_d1_pid) == 13) {
          if (p1daughter1 -> PT > p1daughter2 -> PT) {
            m1_parton = p1daughter1 -> P4();
            m2_parton = p1daughter2 -> P4();
            m1_particle = find_status1_child(branchGenParticle, p1daughter1, 13) -> P4();
            m2_particle = find_status1_child(branchGenParticle, p1daughter2, 13) -> P4();
          } else {
            m1_parton = p1daughter2 -> P4();
            m2_parton = p1daughter1 -> P4();
            m1_particle = find_status1_child(branchGenParticle, p1daughter2, 13) -> P4();
            m2_particle = find_status1_child(branchGenParticle, p1daughter1, 13) -> P4();
          }
        } else if (abs(p2_d1_pid) == 13) {
          if (p2daughter1 -> PT > p2daughter2 -> PT) {
            m1_parton = p2daughter1 -> P4();
            m2_parton = p2daughter2 -> P4();
            m1_particle = find_status1_child(branchGenParticle, p2daughter1, 13) -> P4();
            m2_particle = find_status1_child(branchGenParticle, p2daughter2, 13) -> P4();
          } else {
            m1_parton = p2daughter2 -> P4();
            m2_parton = p2daughter1 -> P4();
            m1_particle = find_status1_child(branchGenParticle, p2daughter2, 13) -> P4();
            m2_particle = find_status1_child(branchGenParticle, p2daughter1, 13) -> P4();
          }
        }
*/

    nPassed+=weight;
    nPassedRaw++;

//------------------------------------------------------------------------------------------------------------------------------------------------------------
// fill histos
//------------------------------------------------------------------------------------------------------------------------------------------------------------

// higgs
    // reco
    hHpTreco -> Fill(h_reco.Pt(),weight);
    hHmreco -> Fill(h_reco.M(),weight); 
    // particle
    hHpTparticle -> Fill(h_particle.Pt(), weight);
    hHmparticle -> Fill(h_particle.M(), weight);
    // parton
    hHpTparton -> Fill(h_parton.Pt(), weight);
    hHmparton -> Fill(h_parton.M(), weight);
    // comp
    hHpTComp -> Fill(h_parton.Pt(), h_reco.Pt(), weight);
    hHmComp -> Fill(h_parton.M(), h_reco.M(), weight);

// jets
    // reco
        // pT + m
        hjjv1pTreco->Fill((j1_reco+j2_reco).Pt(),weight);
        hjjv1mreco->Fill((j1_reco+j2_reco).M(),weight);
        hjjv2pTreco->Fill(j1_reco.Pt()+j2_reco.Pt(),weight);
        hjjv2mreco->Fill(j1_reco.M()+j2_reco.M(),weight);
        hj1pTreco->Fill(j1_reco.Pt(),weight);
        hj1mreco->Fill(j1_reco.M(),weight);
        hj2pTreco->Fill(j2_reco.Pt(),weight);
        hj2mreco->Fill(j2_reco.M(),weight);
        // eta
        hj1etareco->Fill(j1_reco.Eta(), weight); 
        hj2etareco->Fill(j2_reco.Eta(), weight); 
        // phi
        hjjdeltaPhireco->Fill(jjdeltaPhireco,weight); 
        // cos

// z
    //reco
        // pT + m
        hz1pTreco->Fill(z1_reco.Pt(),weight);
        hz2pTreco->Fill(z2_reco.Pt(),weight);
        hz1mreco->Fill(z1_reco.M(),weight);
        hz2mreco->Fill(z2_reco.M(),weight);
        // cos
        hz1cosThetareco->Fill(z1_reco.CosTheta(),weight);
        hz2cosThetareco->Fill(z2_reco.CosTheta(),weight);
        // phi
        hzzdeltaPhireco->Fill(zzdeltaPhireco,weight);
    // particle
        // pT + m
        hz1pTparticle->Fill(z1_particle.Pt(), weight);
        hz2pTparticle->Fill(z2_particle.Pt(), weight);
        hz1mparticle -> Fill(z1_particle.M(), weight);
        hz2mparticle -> Fill(z2_particle.M(), weight);
    // parton
        hz1pTparton->Fill(z1_parton.Pt(),weight);
        hz2pTparton->Fill(z2_parton.Pt(),weight);
        hz1mparton->Fill(z1_parton.M(),weight);
        hz2mparton->Fill(z2_parton.M(),weight);
    // comp
        hz1pTComp->Fill(z1_parton.Pt(), z1_reco.Pt(), weight);
        hz1mComp->Fill(z1_parton.M(), z1_reco.M(), weight);
        hz2pTComp->Fill(z2_parton.Pt(), z2_reco.Pt(), weight);
        hz2mComp->Fill(z2_parton.M(), z2_reco.M(), weight);


//leptons
    // reco
        // pT + m
        hl1pTreco->Fill(l1_reco.Pt(),weight);
        hl2pTreco->Fill(l2_reco.Pt(),weight);
        hl3pTreco->Fill(l3_reco.Pt(),weight);
        hl4pTreco->Fill(l4_reco.Pt(),weight);
        // eta
        // phi
        hl1l2deltaPhireco->Fill(l1l2deltaPhireco,weight);
        hl3l4deltaPhireco->Fill(l3l4deltaPhireco,weight);
        // cos

}

//------------------------------------------------------------------------------------------------------------------------------------------------------------
// write+draw histos
//------------------------------------------------------------------------------------------------------------------------------------------------------------
   
   TFile *hists= new TFile(outputFile,"recreate");
   hists->cd();

// higgs
    hHpTreco -> Write();
    hHmreco -> Write();
    hHpTparticle -> Write();
    hHmparticle -> Write();
    hHpTparton -> Write();
    hHmparton -> Write();
    hHpTComp -> Write();
    hHmComp -> Write();

// jets
    hjjv1pTreco -> Write();
    hjjv1mreco -> Write();
    hjjv2pTreco -> Write();
    hjjv2mreco -> Write();

    hj1pTreco -> Write();
    hj1mreco -> Write();
    hj2pTreco -> Write();
    hj2mreco -> Write();
    hj1etareco -> Write();
    hj2etareco -> Write();

    hjjdeltaPhireco -> Write();

// z 

    hz1pTreco -> Write();
    hz2pTreco -> Write();
    hz1mreco -> Write();
    hz2mreco -> Write();

    hz1cosThetareco -> Write();
    hz2cosThetareco -> Write();

    hzzdeltaPhireco -> Write();

    hz1pTparticle -> Write();
    hz2pTparticle -> Write();
    hz1mparticle -> Write();
    hz2mparticle -> Write();

    hz1pTparton -> Write();
    hz2pTparton -> Write();
    hz1mparton -> Write();
    hz2mparton -> Write();
    hz1pTComp -> Write();
    hz1mComp -> Write();
    hz2pTComp -> Write();
    hz2mComp -> Write();

// leptons 

    hl1pTreco -> Write();
    hl2pTreco -> Write();
    hl3pTreco -> Write();
    hl4pTreco -> Write();
    hl1l2deltaPhireco -> Write();
    hl3l4deltaPhireco -> Write();

    gROOT->SetBatch(kTRUE);

    draw_hist(hHpTreco,"H_pT_reco", "p^{T}_{h}", "pT (GeV)");
    draw_hist(hHpTparticle,"H_pT_particle", "p^{T}_{h}", "pT (GeV)");
    draw_hist(hHpTparton,"H_pT_parton", "p^{T}_{h}", "pT (GeV)");
    draw_hist(hHmreco,"H_m_reco", "m_{h}", "mass (GeV)");
    draw_hist(hHmparticle,"H_m_particle", "m_{h}", "mass (GeV)");
    draw_hist(hHmparton,"H_m_parton", "m_{h}", "mass (GeV)");
    draw_hist2(hHpTComp, "H_pT_comp", "p^{T}_{h}", "parton level pT (GeV)", "reco level pT (GeV)");
    draw_hist2(hHmComp, "H_m_comp", "m_{h}", "parton level mass (GeV)", "reco level mass (GeV)");

    draw_hist(hjjv1pTreco,"jj_pT_reco_v1", "p^{T}_{jj}", "pT (GeV)");
    draw_hist(hjjv1mreco,"jj_m_reco_v1", "m_{jj}", "mass (GeV)");
    draw_hist(hjjv2pTreco,"jj_pT_reco_v2", "p^{T}_{jj}", "pT (GeV)");
    draw_hist(hjjv2mreco,"jj_m_reco_v2", "m_{jj}", "mass (GeV)");
    draw_hist(hj1pTreco,"j1_pT_reco", "p^{T}_{j1}", "pT (GeV)");
    draw_hist(hj1mreco,"j1_m_reco", "m_{j1}", "mass (GeV)");
    draw_hist(hj2pTreco,"j2_pT_reco", "p^{T}_{j2}", "pT (GeV)");
    draw_hist(hj2mreco,"j2_m_reco", "m_{j2}", "mass (GeV)");
    draw_hist(hj1etareco,"j1_eta_reco", "eta_{j1}", "eta");
    draw_hist(hj2etareco,"j2_eta_reco", "eta_{j2}", "eta");
    draw_hist(hjjdeltaPhireco,"jj#phi_reco", "#phi_{jj}", "#phi");
    
    draw_hist(hz1pTreco,"z1_pT_reco", "p^{T}_{z1}", "pT (GeV)");
    draw_hist(hz2pTreco,"z2_pT_reco", "p^{T}_{z2}", "pT (GeV)");
    draw_hist(hz1mreco,"z1_m_reco", "m_{z1}", "pT (GeV)");
    draw_hist(hz2mreco,"z2_m_reco", "m_{z2}", "pT (GeV)");

    draw_hist(hz1pTparticle,"z1_pT_particle", "p^{T}_{z1}", "pT (GeV)");
    draw_hist(hz2pTparticle,"z2_pT_particle", "p^{T}_{z2}", "pT (GeV)");
    draw_hist(hz1mparticle,"z1_m_particle", "m_{z1}", "mass (GeV)");
    draw_hist(hz2mparticle,"z2_m_particle", "m_{z2}", "mass (GeV)");

    draw_hist(hz1cosThetareco,"z1_cos#theta_reco", "cos#theta_{z1}", "cos#theta");
    draw_hist(hz2cosThetareco,"z2_cos#theta_reco", "cos#theta_{z1}", "cos#theta");

    draw_hist(hzzdeltaPhireco,"zz_#phi_reco", "#phi_{zz}", "#phi");

    draw_hist(hz1pTparton,"z1_pT_parton", "p^{T}_{z1}", "pT (GeV)");
    draw_hist(hz2pTparton,"z2_pT_parton", "p^{T}_{z2}", "pT (GeV)");
    draw_hist(hz1mparton,"z1_m_parton", "m_{z1}", "mass (GeV)");
    draw_hist(hz2mparton,"z2_m_parton", "m_{z2}", "mass (GeV)");
    draw_hist2(hz1pTComp, "z1_pT_comp", "p^{T}_{z1}", "parton level pT (GeV)", "reco level pT (GeV)");
    draw_hist2(hz1mComp, "z1_m_comp", "m_{z1}", "parton level mass (GeV)", "reco level mass (GeV)");
    draw_hist2(hz2pTComp, "z2_pT_comp", "p^{T}_{z2}", "parton level pT (GeV)", "reco level pT (GeV)");
    draw_hist2(hz2mComp, "z2_m_comp", "m_{z2}", "parton level mass (GeV)", "reco level mass (GeV)");

    draw_hist(hl1pTreco,"l1_pT_reco", "p^{T}_{l1}", "pT (GeV)");
    draw_hist(hl2pTreco,"l2_pT_reco", "p^{T}_{l2}", "pT (GeV)");
    draw_hist(hl3pTreco,"l3_pT_reco", "p^{T}_{l3}", "pT (GeV)");
    draw_hist(hl4pTreco,"l4_pT_reco", "p^{T}_{l4}", "pT (GeV)");
    draw_hist(hl1l2deltaPhireco,"l1l2_#phi_reco", "#phi_{l1l2}", "#phi");
    draw_hist(hl3l4deltaPhireco,"l3l4_#phi_reco", "#phi_{l3l4}", "#phi");
  
    hists->Close();
}
