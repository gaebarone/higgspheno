#ifndef CUTFLOW_INCLUDE_H
#define CUTFLOW_INCLUDE_H

// typedef std::map<std::string, std::pair<int,double>> cutFlowMapDef;

std::vector <string> cutList_reco;  
std::vector <string> cutList_particle;  
std::vector <string> cutList_parton;  

std::map<std::string, std::vector<string> > cutSelectionProcessReco;
std::map<std::string, std::vector<string> > cutSelectionProcessParticle;
std::map<std::string, std::vector<string> > cutSelectionProcessParton;

void DefineSelections(){

    cutSelectionProcessReco["all"]={"initial - reco", "found higgs - reco", "found vbfjets - reco", "found ww - reco", "found zz - reco"};
    cutSelectionProcessParticle["all"]={"initial - particle", "found higgs - particle", "found vbfjets - particle", "found ww - particle", "found zz - particle"};

    cutSelectionProcessParton["all"]={"initial parton", "Higgs Candidate", "ZZ parton", "WW parton"};

}


bool hasCut(vector<string> cutList, string cut) {
    for(int i=0; i<(int)cutList.size(); i++) {
        if( cut.compare(cutList.at(i))==0 ) return true; 
    }
        return false; 
}

void FillCutFlow(TH1F* hSel, TProfile *hEff,std::map<string, std::pair<int,double>> cutFlowMap, std::vector <string> cutList, string label) {
    for(int i=0; i<(int) cutList.size(); i++) {
        const std::string cutName = cutList[i];
        double passed_reco =  cutFlowMap[cutName].second;
        double efficiency_reco = 100.00 * cutFlowMap[cutName].second / cutFlowMap[cutList[0]].second;

        hSel->GetXaxis()->SetBinLabel(i+1,cutName.c_str());
        hEff->GetXaxis()->SetBinLabel(i+1,cutName.c_str());
        hSel->SetBinContent(i+1,passed_reco);
        hEff->Fill(i+1.0,efficiency_reco);
    } 
}

void PrintCutFlow(std::map<std::string, std::pair<int, double>> cutFlowMap, std::vector<std::string> cutList, std::string label) {
    int width = 15;

    auto printLine = [&]() {    
        std::cout << std::setw(140) << std::setfill('-') << "" << std::setfill(' ') << std::endl;
    };

    // auto printRow = [&](const std::string& name, int passed, double relEff, double efficiency, double normpassed) {
    auto printRow = [&](const std::string& name, int passed, double efficiency, double normpassed) {
        std::cout << "| " << std::setw(width*3) << std::left << name << "|";
        std::cout << std::setw(width*2) << std::left << passed << "|";
        // std::cout << std::setw(width*2) << std::left << relEff << "|";
        std::cout << std::setw(width*2) << std::left << efficiency << "|";
        std::cout << std::setw(width*2) << std::left << normpassed << "|" << std::endl;
    };

    printLine();

    std::cout << "| " << std::setw(width*3) << std::left << label + " Cut" << "|";
    std::cout << std::setw(width*2) << std::left << label + " Passed" << "|";
    // std::cout << std::setw(width*2) << std::left << " Rel Eff " << "|";
    std::cout << std::setw(width*2) << std::left << label + " Efficiency" << "|" ;
    std::cout << std::setw(width*2) << std::left << label + " Norm Count" << "|" << std::endl;

    printLine();

    for (const std::string& cutName : cutList) {
        double passed_reco = cutFlowMap[cutName].first;
        double efficiency_reco = 100.00 * cutFlowMap[cutName].second / cutFlowMap[cutList[0]].second;
        // double relEff = (cutList.size() > 1 && &cutName != &cutList[0]) ? 100.00 * cutFlowMap[cutName].second / cutFlowMap[cutList.at(&cutName - &cutList[1])].second : 100;
        double passedNorm=cutFlowMap[cutName].second;

        // printRow(cutName, passed_reco, relEff, efficiency_reco,passedNorm);
        printRow(cutName, passed_reco, efficiency_reco, passedNorm);
    }

    printLine();
}

void update_cft(std::map<string, std::pair<int,double>> & cft_event, std::map<string, std::pair<int,double>> & cft_total, string cut_name, double weight) {
    cft_event[cut_name] = make_pair(1,weight);
    cft_total[cut_name] = make_pair(cft_total[cut_name].first+1,cft_total[cut_name].second+weight);
}

void increaseCount(std::map<string, std::pair<int,double>> & cutFlowMap, string cutName, double weight) {
    cutFlowMap[cutName]=make_pair(cutFlowMap[cutName].first+1,cutFlowMap[cutName].second+weight);
}

void increaseAllCounts(std::vector<string,std::map<string, std::pair<int,double>>> & allCutFlows,std::map<string,vector<string>> &cutSelectionProcessReco_){
  for(std::map<string,vector<string>>::iterator it=cutSelectionProcessReco_.begin(); it!=cutSelectionProcessReco_.end(); it++){
  }
}

/*
class Histograms {
    initializeHists; 
    fillHists; 
    writeHists 

} 



class Selection : public bla {

    public: 
bool  passSelection(Event); 
}; 

// initialize stuff; 
for (selectionConfiguration in selectionsConfigurartionsWeant )
    selectionMap[selectionConfiguration]=new Slection; 
    selectionMap[selectionConfiguration]->configureAsNeeded; 
    Histograms[selectionConfiguration]->initializeHists; 


// Event Loop 
for ( n in nEvents )

cutflowMap; 
for (selectionConfiguration in selectionsConfigurartionsWeant )
        mapPass[selectionConfiguration]=selectionMap[selectionConfiguration]-passSelection(n)
        FillSelectionAmp(scutflowMap[selectionConfiguration],mapPass[selectionConfiguration]); 
        Histograms[selectionConfiguration]->Fill(n);

    ); 



*/


#endif
