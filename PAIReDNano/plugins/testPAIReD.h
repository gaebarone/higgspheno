#include <cstdlib>
#include <iostream>
#include <functional>
#include <time.h>
#include <iostream>
#include <fstream>
#include <string>
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
#include <vector>

void saveParams(std::vector<std::vector<float>> data, std::vector<std::string> data_names) {
  TFile *hists= new TFile("PAIReDtests.root","recreate");
  hists -> cd();
  for (size_t i = 0; i < data.size(); ++i) {
    float min_val = data[i][0];
    float max_val = data[i][0];
    for (size_t j = 1; j < data[i].size(); ++j) {
      if (data[i][j] < min_val) min_val = data[i][j];
      if (data[i][j] > max_val) max_val = data[i][j];
    }
    TH1 *hData = new TH1F(data_names[i].c_str(), data_names[i].c_str(), 20, min_val, max_val);
    for (size_t j = 0; j < data[i].size(); ++j) {
      hData -> Fill(data[i][j], 1);
    }
    hData -> Write();
    hData -> Clear();
  }
  hists -> Close();
}