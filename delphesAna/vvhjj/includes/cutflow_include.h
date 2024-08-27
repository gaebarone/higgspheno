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

    cut_sel_process_reco["all"]={"initial - reco", "pT_j > 20 && eta_j < 5 - reco", "2 jets - reco", "4 jets - reco", "no bb vbfj - reco", "found jj - reco", "absdetajj > 2.5 - reco", "mjj > 250 - reco",  "no bb found - reco", "found bb - reco", "bb != jj - reco", "pT_l > 15 && eta_l < 2.5 - reco", "2 leps - reco", "4 leps - reco", "2l OS - reco", "4l OS - reco", "mvv < 150 - reco", "mvv > 150 - reco", "met < 30 - reco",  "met > 30 - reco", "found ww - reco", "found zz - reco",  "final wpwmhjj - reco",  "final zzhjj - reco",  "final wpwmjj - reco",  "final zzjj - reco", "final hwpwmjj - reco",  "final hzzjj - reco"};
    cut_sel_process_particle["all"]={"initial - particle", "pT_j > 20 && eta_j < 5 - particle", "2 jets - particle", "4 jets - particle", "no bb vbfj - particle", "found jj - particle", "absdetajj > 2.5 - particle", "mjj > 250 - particle", "no bb found - particle", "found bb - particle", "bb != jj - particle", "pT_l > 15 && eta_l < 2.5 - particle", "2 leps - particle", "4 leps - particle", "2l OS - particle", "4l OS - particle", "mvv < 150 - particle", "mvv > 150 - particle", "met < 30 - particle", "met > 30 - particle", "found ww - particle", "found zz - particle",  "final wpwmhjj - particle",  "final zzhjj - particle",  "final wpwmjj - particle",  "final zzjj - particle",  "final hwpwmjj - particle",  "final hzzjj - particle"};
    cut_sel_process_parton["all"]={"initial - parton", "found bb - parton", "no bb found - parton", "found jj - parton", "found ww - parton", "found zz - parton", "found hww - parton", "found hzz - parton", "final wpwmhjj - parton",  "final zzhjj - parton", "final wpwmjj - parton",  "final zzjj - parton", "final hwpwmjj - parton",  "final hzzjj - parton"};

    // vvhjj
    cut_sel_process_reco["wpwmhjj"]={"initial - reco", "pT_j > 20 && eta_j < 5 - reco", "4 jets - reco", "absdetajj > 2.5 - reco", "mjj > 250 - reco", "no bb vbfj - reco", "found jj - reco", "found bb - reco", "bb != jj - reco", "pT_l > 15 && eta_l < 2.5 - reco", "2 leps - reco", "2l OS - reco", "met > 30 - reco", "mvv > 150 - reco", "found ww - reco",  "final wpwmhjj - reco"};
    cut_sel_process_particle["wpwmhjj"]={"initial - particle", "pT_j > 20 && eta_j < 5 - particle", "4 jets - particle", "absdetajj > 2.5 - particle", "mjj > 250 - particle", "no bb vbfj - particle", "found jj - particle", "found bb - particle", "bb != jj - particle", "pT_l > 15 && eta_l < 2.5 - particle", "2 leps - particle", "2l OS - particle", "met > 30 - particle", "mvv > 150 - particle", "found ww - particle",  "final wpwmhjj - particle"};
    cut_sel_process_parton["wpwmhjj"]={"initial - parton", "found bb - parton", "found jj - parton", "found ww - parton", "final wpwmhjj - parton"};

    cut_sel_process_reco["zzhjj"]={"initial - reco", "pT_j > 20 && eta_j < 5 - reco", "4 jets - reco", "absdetajj > 2.5 - reco", "mjj > 250 - reco", "no bb vbfj - reco", "found jj - reco", "found bb - reco", "bb != jj - reco", "pT_l > 15 && eta_l < 2.5 - reco", "4 leps - reco", "4l OS - reco", "met < 30 - reco", "mvv > 150 - reco", "found zz - reco",  "final zzhjj - reco"};
    cut_sel_process_particle["zzhjj"]={"initial - particle", "pT_j > 20 && eta_j < 5 - particle", "4 jets - particle", "absdetajj > 2.5 - particle", "mjj > 250 - particle", "no bb vbfj - particle", "found jj - particle", "found bb - particle", "bb != jj - particle", "pT_l > 15 && eta_l < 2.5 - particle", "4 leps - particle", "4l OS - particle", "met < 30 - particle", "mvv > 150 - particle", "found zz - particle",  "final zzhjj - particle"};
    cut_sel_process_parton["zzhjj"]={"initial - parton", "found bb - parton", "found jj - parton", "found zz - parton", "final zzhjj - parton"};

    // vvjj
    cut_sel_process_reco["wpwmjj"]={"initial - reco", "pT_j > 20 && eta_j < 5 - reco", "2 jets - reco", "absdetajj > 2.5 - reco", "mjj > 250 - reco", "no bb vbfj - reco", "found jj - reco", "no bb found - reco", "pT_l > 15 && eta_l < 2.5 - reco", "2 leps - reco", "2l OS - reco", "met > 30 - reco", "mvv > 150 - reco", "found ww - reco",  "final wpwmjj - reco"};
    cut_sel_process_particle["wpwmjj"]={"initial - particle", "pT_j > 20 && eta_j < 5 - particle", "2 jets - particle", "absdetajj > 2.5 - particle", "mjj > 250 - particle", "no bb vbfj - particle", "found jj - particle", "no bb found - particle", "pT_l > 15 && eta_l < 2.5 - particle", "2 leps - particle", "2l OS - particle", "mvv > 150 - particle", "met > 30 - particle", "found ww - particle",  "final wpwmjj - particle"};
    cut_sel_process_parton["wpwmjj"]={"initial - parton", "found jj - parton", "no bb found - parton", "found ww - parton", "final wpwmjj - parton"};

    cut_sel_process_reco["zzjj"]={"initial - reco", "pT_j > 20 && eta_j < 5 - reco", "2 jets - reco", "absdetajj > 2.5 - reco", "mjj > 250 - reco", "no bb vbfj - reco", "found jj - reco", "no bb found - reco", "pT_l > 15 && eta_l < 2.5 - reco", "4 leps - reco", "4l OS - reco", "met < 30 - reco", "mvv > 150 - reco", "found zz - reco",  "final zzjj - reco"};
    cut_sel_process_particle["zzjj"]={"initial - particle", "pT_j > 20 && eta_j < 5 - particle", "2 jets - particle", "absdetajj > 2.5 - particle", "mjj > 250 - particle", "no bb vbfj - particle", "found jj - particle", "no bb found - particle", "pT_l > 15 && eta_l < 2.5 - particle", "4 leps - particle", "4l OS - particle", "met < 30 - particle", "mvv > 150 - particle", "found zz - particle",  "final zzjj - particle"};
    cut_sel_process_parton["zzjj"]={"initial - parton", "found jj - parton", "no bb found - parton", "found zz - parton", "final zzjj - parton"};

    // hvvjj
    cut_sel_process_reco["hwpwmjj"]={"initial - reco", "pT_j > 20 && eta_j < 5 - reco", "2 jets - reco", "absdetajj > 2.5 - reco", "mjj > 250 - reco", "no bb vbfj - reco", "found jj - reco", "no bb found - reco", "pT_l > 15 && eta_l < 2.5 - reco", "2 leps - reco", "2l OS - reco", "met > 30 - reco", "mvv < 150 - reco","found ww - reco", "final hwpwmjj - reco"};
    cut_sel_process_particle["hwpwmjj"]={"initial - particle", "pT_j > 20 && eta_j < 5 - particle", "2 jets - particle", "absdetajj > 2.5 - particle", "mjj > 250 - particle", "no bb vbfj - particle", "found jj - particle", "no bb found - particle", "pT_l > 15 && eta_l < 2.5 - particle", "2 leps - particle", "2l OS - particle", "met > 30 - particle", "mvv < 150 - particle", "found ww - particle", "final hwpwmjj - particle"};
    cut_sel_process_parton["hwpwmjj"]={"initial - parton", "found jj - parton", "no bb found - parton", "found hww - parton", "final hwpwmjj - parton"};

    cut_sel_process_reco["hzzjj"]={"initial - reco", "pT_j > 20 && eta_j < 5 - reco", "2 jets - reco", "absdetajj > 2.5 - reco", "mjj > 250 - reco", "no bb vbfj - reco", "found jj - reco", "no bb found - reco", "pT_l > 15 && eta_l < 2.5 - reco", "4 leps - reco", "4l OS - reco", "met < 30 - reco", "mvv < 150 - reco", "found zz - reco", "final hzzjj - reco"};
    cut_sel_process_particle["hzzjj"]={"initial - particle", "pT_j > 20 && eta_j < 5 - particle", "2 jets - particle", "absdetajj > 2.5 - particle", "mjj > 250 - particle", "no bb vbfj - particle", "found jj - particle", "no bb found - particle", "pT_l > 15 && eta_l < 2.5 - particle", "4 leps - particle", "4l OS - particle", "met < 30 - particle", "mvv < 150 - particle", "found zz - particle", "final hzzjj - particle"};
    cut_sel_process_parton["hzzjj"]={"initial - parton", "found jj - parton", "no bb found - parton", "found hzz - parton", "final hzzjj - parton"};

    // no cuts
    cut_sel_process_reco["nocuts"]={"initial - reco"};
    cut_sel_process_particle["nocuts"]={"initial - particle"};
    cut_sel_process_parton["nocuts"]={"initial - parton"};

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
