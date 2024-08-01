#ifndef CUTFLOW_INCLUDE_H
#define CUTFLOW_INCLUDE_H

// typedef std::map<std::string, std::pair<int,double>> cutFlowMapDef;

std::vector <string> cut_list_reco;  
std::vector <string> cut_list_particle;  
std::vector <string> cut_list_parton;  

std::map<std::string, std::vector<string> > cut_sel_process_reco;
std::map<std::string, std::vector<string> > cut_sel_process_particle;
std::map<std::string, std::vector<string> > cut_sel_process_parton;

void DefineSelections(){

    cut_sel_process_reco["all"]={"initial - reco", "pT_j > 20 && eta_j < 5 - reco", "absdetajj > 2.5 - reco", "no bb vbfj - reco", "found jj - reco", "found bb - reco", "no bb veto - reco", "pT_l > 15 && eta_l < 2.5 - reco", "2 passing leps - reco", "4 passing leps - reco", "mll > 10 - reco", "2l OS - reco", "4l OS - reco", "met > 20 - reco", "found ww - reco", "found zz - reco",  "final wwhjj - reco",  "final zzhjj - reco"};
    cut_sel_process_particle["all"]={"initial - particle", "pT_j > 20 && eta_j < 5 - particle", "absdetajj > 2.5 - particle", "no bb vbfj - particle", "found jj - particle", "found bb - particle", "no bb veto - particle", "pT_l > 15 && eta_l < 2.5 - particle", "2 passing leps - particle", "4 passing leps - particle", "mll > 10 - particle", "2l OS - particle", "4l OS - particle", "met > 20 - particle", "found ww - particle", "found zz - particle",  "final wwhjj - particle",  "final zzhjj - particle"};
    cut_sel_process_parton["all"]={"initial - parton", "found bb - parton", "found jj - parton", "found ww - parton", "found zz - parton", "final wwhjj - parton",  "final zzhjj - parton"};

    cut_sel_process_reco["wwhjj"]={"initial - reco", "pT_j > 20 && eta_j < 5 - reco", "absdetajj > 2.5 - reco", "no bb vbfj - reco", "found jj - reco", "found bb - reco", "no bb veto - reco", "pT_l > 15 && eta_l < 2.5 - reco", "2 passing leps - reco", "mll > 10 - reco", "2l OS - reco", "met > 20 - reco", "found ww - reco",  "final wwhjj - reco"};
    cut_sel_process_particle["wwhjj"]={"initial - particle", "pT_j > 20 && eta_j < 5 - particle", "absdetajj > 2.5 - particle", "no bb vbfj - particle", "found jj - particle", "found bb - particle", "no bb veto - particle", "pT_l > 15 && eta_l < 2.5 - particle", "2 passing leps - particle", "mll > 10 - particle", "2l OS - particle", "met > 20 - particle", "found ww - particle",  "final wwhjj - particle"};
    cut_sel_process_parton["wwhjj"]={"initial - parton", "found bb - parton", "found jj - parton", "found ww - parton", "final wwhjj - parton"};

    cut_sel_process_reco["zzhjj"]={"initial - reco", "pT_j > 20 && eta_j < 5 - reco", "absdetajj > 2.5 - reco", "no bb vbfj - reco", "found jj - reco", "found bb - reco", "no bb veto - reco", "pT_l > 15 && eta_l < 2.5 - reco", "4 passing leps - reco", "mll > 10 - reco", "4l OS - reco", "found zz - reco",  "final zzhjj - reco"};
    cut_sel_process_particle["zzhjj"]={"initial - particle", "pT_j > 20 && eta_j < 5 - particle", "absdetajj > 2.5 - particle", "no bb vbfj - particle", "found jj - particle", "found bb - particle", "no bb veto - particle", "pT_l > 15 && eta_l < 2.5 - particle", "4 passing leps - particle", "mll > 10 - particle", "4l OS - particle", "found bb - particle", "found jj - particle", "found zz - particle",  "final zzhjj - particle"};
    cut_sel_process_parton["zzhjj"]={"initial - parton", "found bb - parton", "found jj - parton", "found zz - parton", "final zzhjj - parton"};

    // cut_sel_process_reco["wwjj"]={"initial - reco", "lep pT > 15 && lep eta < 2.5 - reco", "detajj > 2.5 - reco", "no bb vbfj - reco", "2 passing leps - reco", "jet pT > 20 - reco", "mll > 10 -reco", "met > 15 - reco"}; //, "found jj - reco", "found ww - reco"};
    // cut_sel_process_particle["wwjj"]={"initial - particle", "lep pT > 15 && lep eta < 2.5 - particle", "detajj > 2.5 - particle", "no bb vbfj - particle", "2 passing leps - particle", "jet pT > 20 - particle", "mll > 10 - particle", "met > 15 - particle"}; //, "found jj - particle", "found ww - particle"};
    // cut_sel_process_parton["wwjj"]={"initial - parton", "found jj - parton", "found ww - parton", "final - parton"};

    // cut_sel_process_reco["zzjj"]={"initial - reco", "lep pT > 15 && lep eta < 2.5 - reco", "detajj > 2.5 - reco", "no bb vbfj - reco", "4 passing leps - reco", "jet pT > 20 - reco"}; //, "found jj - reco", "found zz - reco"};
    // cut_sel_process_particle["zzjj"]={"initial - particle", "lep pT > 15 && lep eta < 2.5 - particle", "detajj > 2.5 - particle", "no bb vbfj - particle", "4 passing leps - particle", "jet pT > 20 - particle"}; //, "found jj - particle", "found zz - particle"};
    // cut_sel_process_parton["zzjj"]={"initial - parton", "found jj - parton", "found zz - parton", "final - parton"};
    
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
