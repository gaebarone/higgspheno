#include "../../common_includes/tools.h"

void drawXsec(){
  
  // kappa lambda for VVHjj [pb] 
  vector <double> xsecZZhjj= { 0.00000257238, 0.000257238 , 0.0000643096, 0.0000102895, 0.00000257238, 0.00000257238, 0.000257238, 0.0000643096, 0.0000102895 } ; 
  vector <double> kappa_lambdaZZhjj={ 1.0,-10.0,-5,-2,-1,1,10,5,2};
      
      
  vector <double> xsecZZhjjRatiotoSM;
  for(int i=0; i<(int)xsecZZhjj.size(); i++){
    xsecZZhjjRatiotoSM.push_back(xsecZZhjj[i]/xsecZZhjj[0]);
  }


 
  TLegend *leg=DefLeg();
  
  

  TGraph *graphxsecZZhjjRatiotoSM=new TGraph(xsecZZhjjRatiotoSM.size(),&kappa_lambdaZZhjj[0],&xsecZZhjjRatiotoSM[0]);
  graphxsecZZhjjRatiotoSM->GetXaxis()->SetTitle("#it{#lambda}^{BSM}/#it{#lambda}^{SM}");
  graphxsecZZhjjRatiotoSM->GetYaxis()->SetTitle("#it{#sigma}^{BSM}/#it{#sigma}^{SM}");
  
  graphxsecZZhjjRatiotoSM->SetMarkerStyle(8);
  graphxsecZZhjjRatiotoSM->SetMarkerColor(kRed);
  //graphxsecZZhjjRatiotoSM->SetMarkerSize(2);
  graphxsecZZhjjRatiotoSM->SetTitle("");

  TF1 *fgraphxsecZZhjjRatiotoSM=new TF1("fgraphxsecZZhjjRatiotoSM","[0]*x*x",-10,+10);
  fgraphxsecZZhjjRatiotoSM->SetParameter(1,1);
  fgraphxsecZZhjjRatiotoSM->SetLineColor(kRed);
  graphxsecZZhjjRatiotoSM->Fit(fgraphxsecZZhjjRatiotoSM);

  fgraphxsecZZhjjRatiotoSM->SetLineWidth(4);


  // kappaV ZZhjj
  vector <double> xsecZZhjjKappaV={0.025729,0.00160807,0.0000411665,0.0000025729,0.00000000025729,0.00000000025729,0.0000025729,0.0000411665,0.00160807,0.025729};
  vector <double> kappa_VZhjj={-10,-5,-2,-1,0.1,0.1,1,2,5,10};
  vector <double> xsecKappaVZZhjjRatiotoSM;
  for(int i=0; i<(int)xsecZZhjjKappaV.size(); i++){
    xsecKappaVZZhjjRatiotoSM.push_back(xsecZZhjjKappaV[i]/0.0000025729);
  }


  TGraph *graphxsecKappaVZZhjjRatiotoSM=new TGraph(xsecKappaVZZhjjRatiotoSM.size(),&kappa_VZhjj[0],&xsecKappaVZZhjjRatiotoSM[0]);
  graphxsecKappaVZZhjjRatiotoSM->GetXaxis()->SetTitle("#it{#kappa}");
  graphxsecKappaVZZhjjRatiotoSM->GetYaxis()->SetTitle("#it{#sigma}^{BSM}/#it{#sigma}^{SM}");
  
  graphxsecKappaVZZhjjRatiotoSM->SetMarkerStyle(23);
  graphxsecKappaVZZhjjRatiotoSM->SetMarkerColor(kBlue);
  graphxsecKappaVZZhjjRatiotoSM->SetMarkerSize(2);
  graphxsecKappaVZZhjjRatiotoSM->SetTitle("");

  TF1 *fgraphxsecKappaVZZhjjRatiotoSM=new TF1("fgraphxsecKappaVZZhjjRatiotoSM","pol4",-10,+10);
  fgraphxsecKappaVZZhjjRatiotoSM->SetParameter(1,1);
  fgraphxsecKappaVZZhjjRatiotoSM->SetLineColor(kBlue);
  graphxsecKappaVZZhjjRatiotoSM->Fit(fgraphxsecKappaVZZhjjRatiotoSM);

  fgraphxsecKappaVZZhjjRatiotoSM->SetLineWidth(4);

  leg->SetHeader("Modified couplings");
  leg->AddEntry(fgraphxsecZZhjjRatiotoSM,"Higgs self coupling");
  leg->AddEntry(fgraphxsecKappaVZZhjjRatiotoSM,"Higgs to Vector bosons");
  

  TH1F *hempty=new TH1F("hempty","",20,-15,+15);
  SetDef(hempty);
  //hempty->GetXaxis()->SetTitle("#it{#lambda}^{BSM}/#it{#lambda}^{SM}");
  hempty->GetXaxis()->SetTitle("coupling ratio to SM");
  hempty->GetYaxis()->SetTitle("#it{#sigma}^{BSM}/#it{#sigma}^{SM}");
  hempty->GetYaxis()->SetTitleOffset(1.2);
  hempty->SetLineColor(kWhite);
  hempty->SetFillColor(kWhite);
  hempty->GetYaxis()->SetRangeUser(0,110);
  
  TCanvas *cgraphxsecZZhjjRatiotoSM=new TCanvas("graphxsecZZhjjRatiotoSM","");{
    SetCanvasDefaults(cgraphxsecZZhjjRatiotoSM);
    cgraphxsecZZhjjRatiotoSM->cd();
    hempty->Draw("hsit");
    graphxsecZZhjjRatiotoSM->Draw("P same");
    graphxsecKappaVZZhjjRatiotoSM->Draw("P same");
    leg->Draw();
  }
  
  
}
