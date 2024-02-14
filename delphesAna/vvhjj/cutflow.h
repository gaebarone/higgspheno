// typedef std::map<std::string, std::pair<int,double>> cutFlowMapDef;

std::vector <string> cutList_reco;  //{"initial reco", "1 PAIReD jet", "1 bb PAIReD jet", "vbfj pairs","2.5 deltaEta vbf reco","OSFL"};
std::vector <string> cutList_particle;  //{"initial particle", "1 btag particle", "2 good j particle", "2 b-like jet pairs part", "found bb particle", "2 vbfj particle", "comb vbf part","2.5 deltaEta vbf particle","OSFL"};
std::vector <string> cutList_parton;  //{"initial parton", "Higgs Candidate", "ZZ parton"};

std::map<std::string, std::vector<string> > cutSelectionProcessReco;
std::map<std::string, std::vector<string> > cutSelectionProcessParticle;
std::map<std::string, std::vector<string> > cutSelectionProcessParton;

void DefineSelections(){
  // Cut selections to consider
  // ghost seleciton
  cutSelectionProcessReco["all"]={"initial reco", "1 btag reco", "2 good j reco", "2 b-like jet pairs reco", "found bb reco", "2 vbfj reco", "vbfj pairs","2.5 deltaEta vbf reco","OSFL"}; 
  cutSelectionProcessParticle["all"]={"initial particle", "1 btag particle", "2 good j particle", "2 b-like jet pairs part", "found bb particle", "2 vbfj particle", "comb vbf part","2.5 deltaEta vbf particle","OSFL"};
  cutSelectionProcessParton["all"]={"initial parton", "Higgs Candidate", "ZZ parton"};
  
  cutSelectionProcessReco["HZZJJ"]={"initial reco", "1 btag reco", "2 good j reco", "2 b-like jet pairs reco", "found bb reco", "2 vbfj reco", "vbfj pairs","2.5 deltaEta vbf reco","OSFL"}; 
  cutSelectionProcessParticle["HZZJJ"]={"initial particle", "1 btag particle", "2 good j particle", "2 b-like jet pairs part", "found bb particle", "2 vbfj particle", "comb vbf part","2.5 deltaEta vbf particle","OSFL"};
  cutSelectionProcessParton["HZZJJ"]={"initial parton", "Higgs Candidate", "ZZ parton"};
  
  cutSelectionProcessReco["ZZJJ"]={"initial reco", "2 vbfj reco", "vbfj pairs","2.5 deltaEta vbf reco","OSFL"}; 
  cutSelectionProcessParticle["ZZJJ"]={"initial particle", "2 vbfj particle", "comb vbf part","2.5 deltaEta vbf particle","OSFL"};
  cutSelectionProcessParton["ZZJJ"]={"initial parton","ZZ parton"};
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
    int nameWidth = 30;
    int valueWidth = 10;

    auto printLine = [&]() {    
        std::cout << std::setw(nameWidth + valueWidth * 3 + 7) << std::setfill('-') << "" << std::setfill(' ') << std::endl;
    };

    auto printRow = [&](const std::string& name, int passed, double relEff, double efficiency, double normpassed) {
        std::cout << "| " << std::setw(nameWidth) << std::left << name << "|";
        std::cout << std::setw(valueWidth) << std::left << passed << "|";
        std::cout << std::setw(valueWidth) << std::left << relEff << "|";
        std::cout << std::setw(valueWidth) << std::left << efficiency << "|";
        std::cout << std::setw(20) << std::left << normpassed << "|" << std::endl;
    };

    printLine();

    std::cout << "| " << std::setw(nameWidth) << std::left << label + " Cut" << "|";
    std::cout << std::setw(valueWidth) << std::left << label + " Passed" << "|";
    std::cout << std::setw(valueWidth) << std::left << " Rel Eff " << "|";
    std::cout << std::setw(valueWidth) << std::left << label + " Efficiency" << "|" ;
    std::cout << std::setw(valueWidth) << std::left << label + " Norm Count" << "|" << std::endl;

    printLine();

    for (const std::string& cutName : cutList) {
        double passed_reco = cutFlowMap[cutName].first;
        double efficiency_reco = 100.00 * cutFlowMap[cutName].second / cutFlowMap[cutList[0]].second;
        double relEff = (cutList.size() > 1 && &cutName != &cutList[0]) ? 100.00 * cutFlowMap[cutName].second / cutFlowMap[cutList.at(&cutName - &cutList[1])].second : 100;
        double passedNorm=cutFlowMap[cutName].second;

        printRow(cutName, passed_reco, relEff, efficiency_reco,passedNorm);
    }

    printLine();
}

void increaseCount(std::map<string, std::pair<int,double>> & cutFlowMap, string cutName, double weight) {
    cutFlowMap[cutName]=make_pair(cutFlowMap[cutName].first+1,cutFlowMap[cutName].second+weight);
}

void increaseAllCounts(std::vector<string,std::map<string, std::pair<int,double>>> & allCutFlows,std::map<string,vector<string>> &cutSelectionProcessReco_){
  for(std::map<string,vector<string>>::iterator it=cutSelectionProcessReco_.begin(); it!=cutSelectionProcessReco_.end(); it++){
  }
}