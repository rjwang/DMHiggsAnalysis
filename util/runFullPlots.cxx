//
//  runplotter.cc
//
//
//  Created by RENJIE WANG on 2/17/15.
//
//

#include <iostream>
#include <sstream>
#include <fstream>
#include <map>
#include <algorithm>
#include <vector>
#include <set>

#include "TSystem.h"
#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TMath.h"
#include "TStyle.h"
#include "TH2F.h"
#include "TProfile.h"
#include "TROOT.h"
#include "TString.h"
#include "TList.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TLegend.h"
#include "TPaveText.h"
#include "TObjArray.h"
#include "THStack.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TLatex.h"
#include "TAxis.h"

#include<iostream>


using namespace std;

std::string m_outputFolder;
void showOptions(){



  std::cout << "====================================================================== " << std::endl;
  std::cout << "================ Options for runFullPlots.cxx : ===================== " << std::endl;
  std::cout << "====================================================================== " << std::endl;
  std::cout << " " << std::endl;
  std::cout << " " << std::endl;
  std::cout << " " << std::endl;
  std::cout << " -o  output folder name " << std::endl;
  std::cout << " -h  histograms name " << std::endl;

}


bool compareMax(TH1F* hist1, TH1F* hist2){

  return (hist1->Integral() < hist2->Integral());
}

void addText(double x1, double x2, double y1, double y2, TString TEXT, Color_t color, Float_t angle = 0)
{
    TPaveText* T = new TPaveText(x1,y1,x2,y2, "NDC");
    T->SetFillColor(0);
    T->SetFillStyle(0);
    T->SetLineColor(0);
    T->SetTextAlign(22);
    T->SetTextFont(42);
    T->SetTextColor(color);
    TText *text = T->AddText(TEXT);
    text->SetTextAngle(angle);
    text->SetTextAlign(22);
    T->SetTextFont(42);
    T->Draw("same");
    T->SetBorderSize(0);
}


void mergeOverFlow(TH1F *h_)
{
  int binsize = h_->GetXaxis()->GetNbins();
  //underflow bin
  h_->SetBinContent(1,h_->GetBinContent(1)+h_->GetBinContent(0));
  h_->SetBinError(1,sqrt( h_->GetBinError(1)*h_->GetBinError(1)+h_->GetBinError(0)*h_->GetBinError(0) ) );
  
    //overflow bin
  h_->SetBinContent(binsize,h_->GetBinContent(binsize)+h_->GetBinContent(binsize+1));
  h_->SetBinError(binsize,sqrt( h_->GetBinError(binsize)*h_->GetBinError(binsize)+h_->GetBinError(binsize+1)*h_->GetBinError(binsize+1) ) );
  
}

void addHist(TString InputDir = ".root", TString hist="metsig_sel", bool showRatio = false, TString SaveName = "SaveName", bool plotSignal = false, bool isDataBlind = false)
{

  if(isDataBlind) showRatio = false;

  gStyle->SetPadRightMargin (0.03);

  gStyle->SetPadTopMargin(0.07);

  TCanvas *canv = new TCanvas("canv", "limits canvas", 800., 800.);
  canv->cd();
  canv->SetGridx(0);
  canv->SetGridy(0);
  TPad* t1 = new TPad("t1","t1", 0.0, 0., 1.0, 1.0);
  t1->SetRightMargin(0.05);
  TPad* t2;
  if(showRatio) {
    t1 = new TPad("t1","t1", 0.0, 0.25, 1.0, 1.0);
    t1->SetBottomMargin(0.);
    t2 = new TPad("t1","t1", 0.0, 0.0, 1.0, 0.25);
    t2->SetBottomMargin(0.4);
    t2->SetGridy(true);
    t2->Draw();
  }
  t1->Draw();
  t1->cd();
  t1->SetLogy(true);
  //if(hist.Contains("mt_")) t1->SetLogy(false);
  t1->SetLogx(false);

  std::vector<TString> Dirs;
  //  if(InputDir.Contains("plotter.root")){
  Dirs.push_back("V#gamma#gamma");
  Dirs.push_back("V#gamma");
  Dirs.push_back("#gamma#gamma");
  Dirs.push_back("#gamma+jets");
  Dirs.push_back("SMHiggs");
  Dirs.push_back("Mx(1)MZ(200)");
  Dirs.push_back("M_{Z}(1000)M_{A}(200)");

    //  }
  TFile *infile = TFile::Open(InputDir, "READ");

  if( !infile ) {
    std::cout << "Cannot open the input file : " << InputDir << std::endl;
    exit(-1);
  }

  // TList* fileKeys = (TList*) infile->GetListOfKeys();
  // for( TObject* obj : *fileKeys ){
  //   if( )

  // }



  // data
  TH1F* data = (TH1F*) infile->Get("Data/"+hist);
  mergeOverFlow(data);

  if( !data ) {
    std::cout << "Input file : " << InputDir << " doesn't contain the requested histogram Data/" << hist << std::endl;
    exit(-1);
  }
    
  //
  // all background
  //
  TH1F* hist_smh = (TH1F*) data->Clone();
  hist_smh->Reset();


  TH1F* hist_vgammagamma = (TH1F*) data->Clone();
  hist_vgammagamma->Reset();

  TH1F* hist_vgamma           = (TH1F*) data->Clone();
  hist_vgamma->Reset();

  TH1F* hist_gj        = (TH1F*) data->Clone();
  hist_gj->Reset();

  TH1F* hist_gg           = (TH1F*) data->Clone();
  hist_gg->Reset();

  TH1F* hist_totalbkg     = (TH1F*) data->Clone();
  hist_totalbkg->Reset();

  TH1F* h_zprime=NULL;
  TH1F* h_2hdm=NULL;

  if(plotSignal) {
    //signal
    h_zprime  = (TH1F*) infile->Get("Mx(1)MZ(200)/"+hist);
    h_2hdm  = (TH1F*) infile->Get("M_{Z}(1000)M_{A}(200)/"+hist);

    std::cout << " Problem not with 2hdm " << std::endl;
    mergeOverFlow(h_zprime);
    mergeOverFlow(h_2hdm);
    h_zprime->SetLineColor(kRed);
    h_zprime->SetLineWidth(2);


    h_2hdm->SetLineColor(kBlue);
    h_2hdm->SetLineStyle(2);
    h_2hdm->SetLineWidth(2);



  }


  for(size_t j=0; j<Dirs.size(); j++) {
    TH1F* h_ = NULL;
    h_ = (TH1F*) infile->Get(Dirs[j]+"/"+hist);
    if(h_==NULL) continue;
    mergeOverFlow(h_);

    if(!Dirs[j].Contains("Data") 
       && !Dirs[j].Contains("M_{A}_(200)") && !Dirs[j].Contains("MZ(200)")
       ) {
      hist_totalbkg->Add(h_);
    }

    if(Dirs[j].Contains("SMHiggs")) hist_smh->Add(h_);
    if(Dirs[j].Contains("V#gamma#gamma")) hist_vgammagamma->Add(h_);
    if(Dirs[j].Contains("#gamma#gamma")) 	hist_gg->Add(h_);
    if(Dirs[j].Contains("V#gamma")) 	hist_vgamma->Add(h_);
    if(Dirs[j].Contains("#gamma+jets"))           hist_gj->Add(h_);

  }

  for(int ibin=1; ibin<=data->GetXaxis()->GetNbins(); ibin++) {
    if(data->GetBinContent(ibin)==0) {
      int binsize = data->GetBinWidth(ibin);
      data->SetBinError(ibin,1.8/binsize);
    }
  }
	   
	   
  //
  // Error band
  //
	   
	   
  TFile* er_file_sm = TFile::Open("Systematics_SMHiggs.root","READ");
  TGraphAsymmErrors* error_sm = (TGraphAsymmErrors*) er_file_sm->Get("SysBand");
  TFile* er_file_gg = TFile::Open("Systematics_gg.root","READ");
  TGraphAsymmErrors* error_gg = (TGraphAsymmErrors*) er_file_gg->Get("SysBand");
  TGraphAsymmErrors* errors = new TGraphAsymmErrors();

  for(int ibin=1; ibin<=hist_totalbkg->GetXaxis()->GetNbins(); ++ibin) {

    double systematicErrorUp  = 0 ;
    double systematicErrorLow  = 0 ;
    if( hist.Contains("metsig_sel") ){
      int point_sm = TMath::Nint( hist_totalbkg->GetXaxis()->GetBinCenter(ibin)/((error_sm->GetXaxis()->GetXmax() - error_sm->GetXaxis()->GetXmin() )/error_sm->GetN()) );
      int point_gg = TMath::Nint( hist_totalbkg->GetXaxis()->GetBinCenter(ibin)/((error_gg->GetXaxis()->GetXmax() - error_gg->GetXaxis()->GetXmin() )/error_gg->GetN()) );
      
      systematicErrorUp += TMath::Power(error_gg->GetErrorYhigh(point_gg)*hist_gg->GetBinContent(ibin),2);
      systematicErrorUp += TMath::Power(error_sm->GetErrorYhigh(point_sm)*hist_smh->GetBinContent(ibin),2);
      systematicErrorLow += TMath::Power(error_gg->GetErrorYlow(point_gg)*hist_gg->GetBinContent(ibin),2);
      systematicErrorLow += TMath::Power(error_sm->GetErrorYlow(point_sm)*hist_smh->GetBinContent(ibin),2);
    }
    double statisticalError = hist_totalbkg->GetBinError(ibin);
    double errorUp = TMath::Sqrt(TMath::Power( statisticalError , 2 ) +  systematicErrorUp  );  
    double errorDown = TMath::Sqrt(TMath::Power( statisticalError , 2 ) + systematicErrorLow );  
    errors->SetPoint(ibin-1,hist_totalbkg->GetXaxis()->GetBinCenter(ibin), hist_totalbkg->GetBinContent(ibin) );
    errors->SetPointError(ibin-1,hist_totalbkg->GetXaxis()->GetBinWidth(ibin)/2.0,hist_totalbkg->GetXaxis()->GetBinWidth(ibin)/2.0, errorDown , errorUp );
  }

  errors->SetFillStyle(3254); //RJ
  errors->SetFillColor(kGray+3);
  errors->SetLineStyle(1);

  THStack *stack = new THStack("stack","");
  std::vector<TH1F*> orderTH1;
  orderTH1.push_back(hist_vgamma);
  orderTH1.push_back(hist_gg);
  orderTH1.push_back(hist_vgammagamma);
  orderTH1.push_back(hist_gj);
  orderTH1.push_back(hist_smh);

  std::sort(orderTH1.begin() , orderTH1.end() , compareMax );

  for( TH1F* hist : orderTH1)
    stack->Add(hist,"HIST");

  hist_smh->SetFillColor(kCyan-7);
  hist_vgammagamma->SetFillColor(kOrange-2);
  hist_vgamma->SetFillColor(kAzure-4);
  hist_gg->SetFillColor(kTeal+2);
  hist_gj->SetFillColor(kPink+1);

  stack->SetMaximum(stack->GetMaximum()*1000);
  stack->SetMinimum(stack->GetMinimum()*0.01);
  stack->Draw("");


  data->Draw("E1 sames");
  gStyle->SetEndErrorSize(0);
  errors->Draw("2 sames");

  if(plotSignal) {
    h_zprime->Draw("hist same");
    h_2hdm->Draw("hist sames");
  }

  t1->RedrawAxis();

  //Set XTitle
  if(showRatio) {
    stack->GetXaxis()->SetTitle();
    stack->GetXaxis()->SetLabelSize(0);
  } else {
    stack->GetXaxis()->SetTitle("E^{miss}_{T} Significance [#sqrt{GeV}]");
  }
  //Set YTitle
  stack->GetYaxis()->SetTitle("Events / 1 GeV");
  stack->GetYaxis()->SetTitleSize(0.06);
  stack->GetYaxis()->SetTitleOffset(0.8);



  TLegend* leg  = new TLegend(0.5,0.7,0.94,0.95, "NDC");

  if(!showRatio) leg  = new TLegend(0.5,0.7,0.92,0.95, "NDC");

  leg->AddEntry(data,"Data","lep");
  leg->SetNColumns(2);
  leg->SetTextSize(0.032);
  leg->AddEntry(hist_gj,"#gamma+jets","F");
  leg->AddEntry(hist_vgammagamma,"V#gamma#gamma","F");
  leg->AddEntry(hist_gg,"#gamma#gamma","F");
  leg->AddEntry(hist_vgamma,"V#gamma","F");
  leg->AddEntry(hist_smh,"SMHiggs","F");
  leg->AddEntry(errors,"Syst. #oplus Stat. Unc.","F");

  if(plotSignal) {
    leg->AddEntry(h_zprime,"Mx(1)MZ(200)","l");
    leg->AddEntry(h_2hdm,"M_{Z}(1000)M_{A}(200)","l");


  }
  leg->SetFillColor(0);
  leg->SetFillStyle(0);
  leg->SetLineColor(0);
  leg->SetHeader("");
  leg->Draw("same");

  TLatex latex;
  latex.SetTextSize(0.04);
  latex.SetNDC();
  latex.SetTextAlign(22);
  latex.SetTextFont(42);
  latex.SetTextSizePixels(22);


  latex.DrawLatex(0.3,0.85,"#bf{#it{ATLAS}} Internal");
  latex.DrawLatex(0.3,0.75,"#sqrt{s} = 13TeV, #int #it{L}dt = 15.4 fb^{-1}");

  if(showRatio) {
    t2->cd();
    //t2->SetLogy(1);
    TH1F* h_ratio1 = (TH1F*) data->Clone();
    //cout << h_ratio1->GetBinError(19) << endl;
    h_ratio1->Divide(hist_totalbkg);

    TGraphAsymmErrors* h_ratio = new TGraphAsymmErrors(errors->GetN());

    for(int ibin=0; ibin<errors->GetN(); ++ibin) {
      double xvalue = 0 ;
      double yvalue = 0 ;
      errors->GetPoint(ibin,xvalue,yvalue);
      double errx=h_ratio1->GetXaxis()->GetBinWidth(h_ratio1->GetXaxis()->FindBin(xvalue))/2.;
      double errydown = yvalue !=0  ? errors->GetErrorYlow(ibin)/yvalue : 0 ;
      double erryup = yvalue !=0  ? errors->GetErrorYhigh(ibin)/yvalue : 0 ;
      h_ratio->SetPoint(ibin,xvalue,1);
      h_ratio->SetPointError(ibin,errx,errx,errydown,erryup);
    } 

    h_ratio1->Draw("");
    h_ratio->Draw("E1");
    gStyle->SetOptStat(0);


    h_ratio1->SetMinimum(-0.0001);
    h_ratio1->SetMaximum(2.1);


    h_ratio1->GetXaxis()->SetTitle("E^{miss}Significance_{T} [#sqrt{GeV}]");
    h_ratio1->GetYaxis()->SetTitleOffset(0.44);
    h_ratio1->GetYaxis()->SetLabelSize(0.15);
    h_ratio1->GetYaxis()->SetNdivisions(/*5*/3);
    h_ratio1->GetYaxis()->SetTickLength(0.02);
    h_ratio1->GetYaxis()->SetTitleSize(/*0.17*/0.15);
    h_ratio1->GetYaxis()->SetTitleOffset(/*0.40*/0.3);

    h_ratio1->GetXaxis()->SetTitleOffset(0.8);
    h_ratio1->GetXaxis()->SetLabelSize(0.16);
    h_ratio1->GetXaxis()->SetTitleSize(/*0.20*/0.21);
    h_ratio1->GetXaxis()->SetTickLength(0.12);


    h_ratio1->GetXaxis()->SetNdivisions(5 + 100*2 + 10000*5);
    h_ratio1->GetYaxis()->SetNdivisions(3 + 100*5 + 10000*0);

    h_ratio1->GetYaxis()->SetTitle("Data / MC");


    // ============== Plot stat + sys uncertainties ===================


    TGraphAsymmErrors *denRelUnc=new TGraphAsymmErrors(errors->GetN());
    TH1D *denRelUncH=(TH1D *) hist_totalbkg->Clone("mcrelunc");
    double llmin=-999;
    double llmax=999;

    for(int ibin=0; ibin<errors->GetN(); ++ibin) {
      if(ibin == 1) llmin=denRelUncH->GetXaxis()->GetBinLowEdge(ibin);
      llmax=denRelUncH->GetXaxis()->GetBinUpEdge(ibin);
      double xvalue = 0 ;
      double yvalue = 0 ;
      errors->GetPoint(ibin,xvalue,yvalue);
      double errx=denRelUncH->GetXaxis()->GetBinWidth(denRelUncH->GetXaxis()->FindBin(xvalue))/2.;
      double errydown = yvalue !=0  ? errors->GetErrorYlow(ibin)/yvalue : 0 ;
      double erryup = yvalue !=0  ? errors->GetErrorYhigh(ibin)/yvalue : 0 ;
      denRelUnc->SetPoint(ibin,xvalue,1);
      denRelUnc->SetPointError(ibin,errx,errx,errydown,erryup);
    } 

    denRelUnc->SetLineColor(1);
    denRelUnc->SetFillStyle(3254);
    denRelUnc->SetFillColor(kRed);
    denRelUnc->SetMarkerColor(1);
    denRelUnc->SetMarkerStyle(1);
    denRelUnc->Draw("2");

    TLegend* leg2  = new TLegend(0.17,0.70,0.5,1.05, "NDC");
    leg2->AddEntry(denRelUnc,"Stat. #oplus syst. unc.","F");
    
    leg2->SetFillColor(0);
    leg2->SetFillStyle(0);
    leg2->SetLineColor(0);
    leg2->SetHeader("");
    leg2->Draw("same");
    leg2->SetTextFont(42);

    TLine *llbase;
    llbase = new TLine(llmin,1.,llmax,1.);
    llbase->SetLineWidth(1);
    llbase->SetLineStyle(1);
    llbase->SetLineColor(kBlue);
    llbase->Draw("same");

  }

  infile->Close();

  canv->SaveAs(TString(m_outputFolder)+SaveName+".pdf");
  canv->SaveAs(TString(m_outputFolder)+SaveName+".png");
  delete canv;

}

int main(int argc, char** argv)
{

  if( argc < 2 ){
    std::cout << "Input arguments are zero. This algorithm works with two inputs: one the output folder and the second, the filename containing histograms for systematics. Exiting." << std::endl;
    return 0;
  }



  std::vector<std::string> list;
  std::map<std::string,std::vector<std::string> > options;

  for( int i = 0 ;  i < argc ; ++i ){


    std::string argStr=argv[i];
    if( TString(argStr).BeginsWith("-") )
      argStr=argv[i+1];
    else
      continue;
    for (size_t argI=0,n; argI <= argStr.length(); argI=n+1){
      n = argStr.find_first_of(',',argI);
      if (n == std::string::npos)
	n = argStr.length();
      string tmp = argStr.substr(argI,n-argI);
      list.push_back(tmp);
    }
    std::string option = argv[i]; 
    options[option.substr(1,1)]=list;
    list.clear();
  }


  std::string outputFolder;
  std::string inputFiles;
  std::vector<std::string> histogramStem;

  for( auto itMap : options){

    if( itMap.first == "o" )
      outputFolder = itMap.second[0];
    else if( itMap.first == "h" )
      histogramStem = itMap.second;
    else if( itMap.first == "f" )
      inputFiles = itMap.second[0];
    else 
      showOptions();
  }


  m_outputFolder = TString(outputFolder).EndsWith("/") ? outputFolder : outputFolder+"/" ;

  if( inputFiles == "" )
    inputFiles = "/afs/cern.ch/work/a/alopezso/private/DMAnalysis20.7/Plots/monohiggs_paper_comb1615_data16_h013_v1_Lite/plotter.root";
  bool isDataBlind = true;

  for( std::string histoName : histogramStem )
    addHist(TString(inputFiles),TString(histoName),true,TString(histoName),true);



  return 0;

}
