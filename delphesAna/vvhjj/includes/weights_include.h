#ifndef WEIGHTS_INCLUDE_H
#define WEIGHTS_INCLUDE_H


Long64_t get_num_entries(const char *inputName) {
  TChain chain("Delphes");
  chain.Add(inputName);
  ExRootTreeReader *treeReader = new ExRootTreeReader(&chain);
  return treeReader->GetEntries();
}

Long64_t get_total_num_entries(const char *process_name) {
  std::string inputFileName = std::string(process_name) + "_inputs.txt";
  std::ifstream inputFile(inputFileName.c_str());
  if(inputFile.fail()){
    cout<<" File "<<inputFileName<<" does not exist returing -1 total num entries"<<endl;
    return -1;
  }
  std::string line;
  TChain chain("Delphes");
  Long64_t total=0;
  while (std::getline(inputFile, line)) {
    total += get_num_entries(line.c_str());
  }
  return total;
}

Long64_t get_total_events(const char *process_name) {
  std::string inputFileName = std::string(process_name) + "_inputs.txt";
  std::ifstream inputFile(inputFileName.c_str());
  if(inputFile.fail()){
    cout<<" File "<<inputFileName<<" does not exist returing -1 total events"<<endl;
    return -1;
  }
  std::string line;
  TChain chain("Delphes");
  Long64_t total = 0;
  while (std::getline(inputFile, line)) {
    chain.Add(line.c_str());
  }
  ExRootTreeReader *treeReader = new ExRootTreeReader(&chain);
  Long64_t numEntries = treeReader->GetEntries();
  delete treeReader;
  return numEntries;
}


#endif