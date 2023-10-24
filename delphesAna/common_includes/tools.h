/*---------------------------------------------------------------------------------------

tools.h

The purpose of this header file is to collate common code snippets for other macros. 

The snippets include: 

  FUNCTION         | PURPOSE
  ----------------------------------------------------------------------------------------
  DivideBinW       | Divide the contents of a histogram by their bin widths
  GiveSigHist      | Given two input histograms, calculate the significance of 
                   | histogram 1 relative to histogram 2
  DefLeg           | Initialize a legend with some default properties
  RadioUncert      | Determine the uncertainty in a ratio of x and y
  Round            | ?
  channeLegend     | Create a legend labelled according to the different channels
  Lumi             | Calculate a label reporting the luminosity
  ATLASLeg         | Create a legend bearing the ATLAS header, with or without the 
                   | 'simulation' label
  ATLASLabel       | Create the ATLAS plot label with sqrt(s) and the luminosity
  ReadableName     | Convert strings containing underscores to strings with spaces,
                   | ie, make "human readable"
  RemoveBlanks     | Remove blanks/spaces from a string and replace them with underscores
                   | ie, make "computer readable"
  SetCanvasDefaults| Assign default values for canvas margins
  PtConeVar        | Create a string for the ptvarcone20 variable
  GetVar           | ?
  RootToLatex      | Convert ROOT's version of a Latex string into .tex syntax
  MoveH            | Shift histogram contents to a new central value
  MoveHNO          | Shift histogram contents to a new central value with a new error and
                   | and a new name
  RemoveEtaRad     | ?
---------------------------------------------------------------------------------------*/

#include "TApplication.h"
#include "TSystem.h"
#include "TROOT.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <string>
#include "TFile.h"
#include "TH1.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TLine.h"
#include "TF1.h"
#include "TFile.h"
#include "TROOT.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TSystem.h"
#include "TF1.h"

#include "RooWorkspace.h"
#include "RooAbsData.h"
#include "RooStats/NumberCountingUtils.h"
#include "RooStats/ModelConfig.h"
#include "RooStats/FeldmanCousins.h"
#include "RooStats/ToyMCSampler.h"
#include "RooStats/PointSetInterval.h"
#include "RooStats/ConfidenceBelt.h"
#include "RooDataHist.h"
#include "RooGenericPdf.h"
#include "RooRealVar.h"
#include "RooAbsPdf.h"
#include "RooHistPdf.h"
#include "RooStats/ProfileLikelihoodCalculator.h"
#include "RooStats/LikelihoodInterval.h"
#include "RooStats/LikelihoodIntervalPlot.h"
#include "RooStats/HypoTestResult.h"
#include "RooStats/ConfInterval.h"
#include "RooStats/ProfileLikelihoodTestStat.h"
#include "RooStats/SamplingDistribution.h"
#include "RooStats/SamplingDistPlot.h"
#include "RooStats/RooStatsUtils.h"
#include "RooStats/HypoTestResult.h"
#include "RooStats/SamplingDistPlot.h"
#include "RooStats/ModelConfig.h"
#include "TMath.h"
#include <vector>
//#include <stdio>
//#include <time>

//#include " /Users/gaetano/Atlas/programs/releases/AnalysisBase/2.3.28/RootCore/scripts/load_packages.C"

using namespace RooFit;
using namespace RooStats;
using namespace std;

void DivideBinW(TH1F *h,bool derr=true){
//##:::PURPOSE -- Divide the contents of a histogram by their bin widths.
 
 if(h==NULL) return;
  h->Sumw2(true);

  double entr=h->GetEntries();
  double inte=h->Integral();
  for( int b=1; b<(int)h->GetNbinsX()+1; b++)
    {
      double cont=h->GetBinContent(b);
      double nEvt=h->GetBinContent(b)/inte;
      if( h->GetBinWidth(b)!=0){
	h->SetBinContent(b,h->GetBinContent(b)/h->GetBinWidth(b));
	//double err=sqrt(nEvt*entr)*h->GetBinContent(b)/cont/h->GetBinWidth(b);
	//h->SetBinError(b,err);
	h->SetBinError(b,h->GetBinError(b)/h->GetBinWidth(b));
      }
      else {
	h->SetBinContent(b,0);
	h->SetBinError(b,0);
      }
    }
}


TH1F *GiveSigHist(TH1F *h1,TH1F *h2,double stat=1){
//##::: PURPOSE -- Given two input histograms, calculate significance of histogram 1 relative to histogram 2
  
  TH1F *h=(TH1F*)h1->Clone("h");
  h->GetYaxis()->SetTitle("signifiance [#sigma]");
  h->GetXaxis()->SetTitle(h1->GetXaxis()->GetTitle());
  
  for(int i=1;i<(int)h1->GetNbinsX()+1;i++){
    // double sig=NumberCountingUtils::BinomialObsZ(h1->GetBinContent(i),h2->GetBinContent(i),
    //h2->GetBinError(i)/h2->GetBinContent(i) );
    //double sig=(h1->GetBinContent(i)/stat - h2->GetBinContent(i)/stat)/sqrt(h2->GetBinContent(i)*stat);
    double sig=(h1->GetBinContent(i)*stat - h2->GetBinContent(i)*stat)/sqrt(h2->GetBinContent(i)*stat+h1->GetBinError(i)*stat);
    h->SetBinContent(i,sig);
  }
  return h;  
}

TLegend *DefLeg(double x1=0.6,double x2=0.7,double x3=.75,double x4=0.85,double size=22){
//##::: PURPOSE -- initialize a legend with some default properties
  
  TLegend* legend=new TLegend(x1,x2,x3,x4);
  legend->SetFillColor(kWhite);
  legend->SetBorderSize(0);
  legend->SetTextSize(22);
  legend->SetTextFont(43);

  return legend;
}

double RatioUncert(double x=0,double xe=0,double y=1,double ye=0){
//##::: PURPOSE -- Determine the uncertainty on a ratio
  
  double ret=0;
  if(x==0) x=1;
  ret+=pow(xe/y,2);
  ret+=pow(x/(y*y)*ye,2);
  ret=sqrt(ret);
  return ret;
}

int Round( double x , double base = 5 ){
  x = x * pow(10, base);
  x = x + 0.5;
  x = floor(x) / pow(10, base);
  return x;

  /*

//##::: PURPOSE --
  int k=0;
  int max = 10 ; 
  
  if( x < pow(10,-base) ){
    for( k = 0 ; k < max ; k++)    
      if (  x*pow(10,k) > pow(10,base) && x * pow(10,k)  <= pow(10,base+1)     ) 
	return -k;
  }
      else {
        for( k = max ; k > base ; k--)  
          if (  x / pow(10,k)  >  1   && x/pow(10,k) <  10  ) 
            return fabs(k);
      }
  return 0;
  */
}

TLegend *channeLegend(bool el=true){
//##::: PURPOSE -- Create legend labelled according to the different channels
  
  TLegend* chLegendM=new TLegend(0.24,0.64,0.25,0.7);
  chLegendM->SetFillColor(kWhite);
  chLegendM->SetBorderSize(0);
  chLegendM->SetTextSize(0.06);
  chLegendM->SetTextFont(42);
  /*
  if(el)
    chLegendM->AddEntry( (TObject*)NULL,"#it{e}+jets","");
  else 
    chLegendM->AddEntry( (TObject*)NULL,"#it{#mu}+jets","");
  */

 if(el)
    chLegendM->AddEntry( (TObject*)NULL,"#bf{Electron channel}","");
  else 
    chLegendM->AddEntry( (TObject*)NULL,"#bf{Muon channel}","");
  return chLegendM;
}

void Lumi(TCanvas &c1,double lumi=4.59,double s=7,double x=0.1,double y=0.8){
//##::: PURPOSE -- calculate a label reporting the luminosity
  
  TLatex t;
  string dataSet="#scale[0.5]{ #sqrt{s} = 13 TeV #int#it{L}dt = ";
  char slumi[100];
  sprintf(slumi,"%.2f",lumi);
  dataSet+=slumi;
  dataSet+=" fb^{-1} }";
  c1.cd();
  t.DrawLatex(x,y,dataSet.c_str());
}

TLegend *ATLASLeg(bool simulation=false,int Wk=1,double lumi=4.59,
		  double x1=0.17, double y1=0.74,double x2=0.2,double y2=0.89){
//##::: PURPOSE -- create a legend bearing the ATLAS header, with or without 'simulation' label

  TLegend *ATLAS=new TLegend(x1,y1,x2,y2);
  ATLAS->SetFillColor(kWhite);
  ATLAS->SetBorderSize(0);
  //ATLAS->SetTextSize(0.04);
  ATLAS->SetTextFont(43);
  
  TLegend *ATLASIM=(TLegend*)ATLAS->Clone("ATLASIM");
  
  string s;
  if(Wk==0)
    s="#it{#bf{ATLAS}}";
  if(Wk==1)
    s="#it{#bf{ATLAS}} Internal";
  if(Wk==2)
    s="#it{#bf{ATLAS}} Internal";

  if(Wk!=-1)
    ATLAS->AddEntry((TObject*)NULL,s.c_str(),"");

  if(Wk==1 && simulation){
    ATLASIM->AddEntry((TObject*)NULL,"#font[43]{ #it{#bf{ATLAS}} Simulation,}","");
    ATLASIM->AddEntry((TObject*)NULL,"#font[43]{ #scale[0.8]{Internal}}","");
  }
  else if(Wk!=-1)
    ATLASIM->AddEntry((TObject*)NULL,"#font[43]{ #it{#bf{ATLAS}} Simulation Internal}","");
  
  if(lumi!=0){
    /*
    string dataSet="#scale[0.5]{ #sqrt{s} = 7 TeV #int#it{L}dt = ";
    char slumi[100];
    //sprintf(slumi,"%.2f",4.66);
    sprintf(slumi,"%.2f",lumi);
    dataSet+=slumi;
    dataSet+=" fb^{-1} }";
    */
    char dataSet[1000];
    sprintf(dataSet,"#sqrt{s} = 13 TeV, %.1f fb^{-1}",lumi);
    //sprintf(dataSet,"#scale[0.7]{ #it{#sqrt{s}} = 13 TeV, %.1f fb^{-1} }",lumi);
    //sprintf(dataSet,"#scale[0.7]{#it{#sqrt{s}} = 8 TeV #int#it{L}dt = %.2f fb^{-1} }",lumi);

    ATLAS->AddEntry((TObject*)NULL,dataSet,"");
  }
  if(!simulation) return ATLAS;
  else return ATLASIM;
}


void ATLASLabel(TLatex *ATLAS=NULL,double x1=0.17,double y1=0.74,bool simulation=false,int Wk=1,double lumi=4.59){
//##::: PURPOSE -- Create the ATLAS plot label with sqrt(s) and the luminosity 
  
  if(ATLAS == NULL ) return; 
  ATLAS->SetTextSize(22);
  ATLAS->SetTextFont(43);
  
  char *s;
  if( Wk == 0 )
    s = ( char* ) "#font[42]{ #it{#bf{ATLAS}}}";
  if( Wk==1)
    s = ( char* ) "#font[42]{ #it{#bf{ATLAS}} Internal}";
  if( Wk==2 )
    s = ( char* ) "#font[42]{ #it{#bf{ATLAS}} Internal}";
  
  if( simulation )
    s = ( char* ) "#font[42]{ #it{#bf{ATLAS}} Simulation Internal}";
  ATLAS->DrawLatex( x1,y1,s );
  char dataSet[ 100 ];
  //sprintf(dataSet,"#scale[1.0]{ #sqrt{s} = 13 TeV, %.1f fb^{-1} }",lumi);
  sprintf(dataSet,"#scale[1.0]{ #sqrt{#it{s}} = 14 TeV}");
  //sprintf(dataSet,"#scale[0.5]{ #sqrt{s} = 8 TeV #int#it{L}dt = %.2f fb^{-1} }",lumi);
  //if(!simulation)
  ATLAS->DrawLatex(x1,y1*0.9,dataSet);
  
}


string ReadableName(string r){
//##::: PURPOSE -- Convert strings containing underscores to strings with spaces, ie "human readable"
 
 string s=r;
  while (s.find("_")!=string::npos){
    size_t f=s.find("_");
    s.replace(f,std::string("_").length()," ");
  }
  return s;
}

string RemoveBlanks(string r){
//##::: PURPOSE -- Remove blanks/spaces from a string and replace them with underscores, ie make "computer readable"
  string s=r;
  while (s.find(" ")!=string::npos){
    size_t f=s.find(" ");
    s.replace(f,std::string(" ").length(),"_");
  }
  return s;
}

void SetCanvasDefaults(TCanvas *c1){
//##::: PURPOSE -- assign default values for canvas margins
  c1->SetTicks();
  c1->SetTopMargin(0.015);
  c1->SetRightMargin(0.020);
  c1->SetBottomMargin(0.15);
  c1->SetLeftMargin(0.12);
}

char *PtConeVar(){
//##::: PURPOSE -- create a string for the ptvarcone20 variable 
  return ( char* ) "#it{p}_{T}^{20} [GeV]";
}

string GetVar(TH1F *h=NULL){
//##::: PURPOSE --
  string b=h->GetXaxis()->GetTitle();
  b.erase(b.find("[")-1,b.size());
  return b;
}

string  RootToLatex(string s){
//##::: PURPOSE -- Convert ROOT's version of a Latex string into .tex syntax
  while (s.find("#")!=string::npos){
    size_t f=s.find("#");
    s.replace(f, std::string("#").length(),"\\");
  }

  
  
  return s;
}

TH1F *MoveH(TH1F *h, double center=-9999){
//##::: PURPOSE -- Shift histogram contents to a new central value
  if(center ==-9999) 
    center=h->GetMean();

  char newname[100];
  sprintf(newname,"%s%s",h->GetName(),"_mv");
  TH1F *hn=new TH1F(newname,h->GetTitle(),h->GetNbinsX(),h->GetBinLowEdge(1)-center,
		    h->GetBinLowEdge(h->GetNbinsX()+1)-center);
  hn->GetXaxis()->SetTitle(h->GetXaxis()->GetTitle());
  hn->GetYaxis()->SetTitle(h->GetYaxis()->GetTitle());
  hn->SetMarkerStyle(h->GetMarkerStyle());
  hn->SetMarkerSize(h->GetMarkerSize());
  hn->SetMarkerColor(h->GetMarkerColor());
  hn->SetFillColor(h->GetFillColor());
  hn->SetLineStyle(h->GetLineStyle());
  hn->SetLineWidth(h->GetLineWidth());
  hn->SetLineColor(h->GetLineColor());

  for(int i=1;i<(int)h->GetNbinsX()+1;i++){
    double cont=h->GetBinContent(i);
    double oldCent=h->GetBinCenter(i);
    int newBinCent=hn->FindBin(oldCent-center);
    hn->SetBinContent(newBinCent,hn->GetBinContent(newBinCent)+cont);
  }
  return hn;
}



void MoveHNO(TH1F *h, double center=-9999){
//##::: PURPOSE -- Shift histogram contents to a new central value
  if(center ==-9999) 
    center=h->GetMean();

  char newname[100];
  sprintf(newname,"%s%s",h->GetName(),"_mv");
  TH1F *hn=new TH1F(newname,h->GetTitle(),h->GetNbinsX(),h->GetBinLowEdge(1)-center,
		    h->GetBinLowEdge(h->GetNbinsX()+1)-center);
  hn->GetXaxis()->SetTitle(h->GetXaxis()->GetTitle());
  hn->GetYaxis()->SetTitle(h->GetYaxis()->GetTitle());
  hn->SetMarkerStyle(h->GetMarkerStyle());
  hn->SetMarkerSize(h->GetMarkerSize());
  hn->SetMarkerColor(h->GetMarkerColor());
  hn->SetFillColor(h->GetFillColor());
  hn->SetLineStyle(h->GetLineStyle());
  hn->SetLineWidth(h->GetLineWidth());
  hn->SetLineColor(h->GetLineColor());

  for(int i=1;i<(int)h->GetNbinsX()+1;i++){
    double cont=h->GetBinContent(i);
    double oldCent=h->GetBinCenter(i);
    int newBinCent=hn->FindBin(oldCent-center);
    hn->SetBinContent(newBinCent,hn->GetBinContent(newBinCent)+cont);
    hn->SetBinError(newBinCent,h->GetBinError(i));
  }
  hn->SetName(h->GetName());
  h->SetName(newname);
  h=hn;
}


void RemoveEtaRad(TH1F *h){
//##::: PURPOSE -- 
  string hName=h->GetName();
  
  if(hName.find("Phi")!=string::npos || 
     hName.find("phi")!=string::npos || 
     hName.find("PHI")!=string::npos) 
    return ;

  //if( !(hName.find("eta")!=string::npos || hName.find("Eta")!=string::npos || 
  //hName.find("ETA")!=string::npos)) 
  // return ;

  
  vector <string> ifind;
  ifind.push_back("[rad]"); ifind.push_back("rad");

  string iaxis=h->GetXaxis()->GetTitle();
  for( vector <string>::iterator it=ifind.begin(); it!=ifind.end(); it++){
    string is=*it;
    size_t f = iaxis.find(is);
    if(f!=string::npos)
      iaxis.replace(f, is.length(),""); 
  }
  h->GetXaxis()->SetTitle(iaxis.c_str());

  string yaxis=h->GetYaxis()->GetTitle();
  for( vector <string>::iterator it=ifind.begin(); it!=ifind.end(); it++){
    string is=*it;
    size_t f = yaxis.find(is);
    if(f!=string::npos)
      yaxis.replace(f, is.length(),""); 
  }
  h->GetYaxis()->SetTitle(yaxis.c_str());
  
  if(string(h->GetName()).compare("hPhPtCone20")==0) 
    h->GetXaxis()->SetTitle("#it{p}^{iso}_{T} [GeV]");

  if(string(h->GetName()).compare("hPhotonPt")==0) 
    h->GetXaxis()->SetTitle("#it{E}_{T}(#gamma) [GeV]");

}


void SetDef(TH1 *h){
  h->SetStats(0);
  h->SetLineWidth(2);
  h->Sumw2(true);
  
  h->GetYaxis()->SetTitleFont(43);
  h->GetYaxis()->SetTitleSize(25);
  h->GetYaxis()->SetLabelFont(43);
  h->GetYaxis()->SetLabelSize(20);
  h->GetYaxis()->SetTitleOffset(1.2);

  h->GetXaxis()->SetTitleFont(43);
  h->GetXaxis()->SetTitleSize(25);
  h->GetXaxis()->SetLabelFont(43);
  h->GetXaxis()->SetLabelSize(20);
  h->GetYaxis()->SetTitleOffset(0.9);
  //h->GetXaxis()->SetTitleOffset(2.4);
  h->SetNdivisions(505);

  h->SetLineColor(kBlack);
}


void GlobalShutUp(){   
  RooMsgService::instance().setSilentMode(kTRUE);
  RooMsgService::instance().setGlobalKillBelow(RooFit::FATAL);      
}
