#include <iostream>

#include "HGamAnalysisFramework/RunUtils.h"
#include "HGamAnalysisFramework/Config.h"
#include "DMHiggsAnalysis/DataEvtSummaryHandler.h"
#include "DMHiggsAnalysis/PhysicsEvent.h"

#include "TSystem.h"
#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TGraphAsymmErrors.h"
#include "TGraph.h"
#include "TMultiGraph.h"
#include "TProfile.h"
#include "TEventList.h"
#include "TROOT.h"
#include "TMath.h"

#include "DMHiggsAnalysis/SmartSelectionMonitor.h"
#include "DMHiggsAnalysis/deltaR.h"

std::string m_outputFolder;

using namespace std;

bool symemtricCuts(double met, double ptyy,int cut){

  return   met >= cut*10 && ptyy >= cut*10;

}
bool Categorisation(std::vector<double> variables,int category){


  //  double met = variables[0];
  double metSig = variables[0];
  double ptyy = variables[1];
  double ptHard = variables[2];


  switch(category){


  case 0 :

    return true;
    break;
  case 4 :

    if( ptyy >= 15 && ( ptyy <= 25  ||  metSig <= 4 ) )
      return true;


    break;
  case 3 :
    if(  ptyy > 25  && ( metSig > 4 && metSig <= 7 ) )
      return true;

    
    break;
  case 2 :
    
    if(  ptyy <= 90  && ( metSig > 7 ) )
      return true;

    break;
  case 1 :

    if(  ptyy > 90  && ( metSig > 7 ) )
      return true;    
    break;

  default :

    std::cout << " Category number " << category << " is not defined. Currently there are 4 different categories." << std::endl;
    
    break;

  }



  return false;
}

int getCategory(std::vector<double> variables){


  //  double met = variables[0];
  double metSig = variables[0];
  double ptyy = variables[1];
  //  double ptHard = variables[2];


  if(  metSig > 7  ){
    if( ptyy > 90 )
      return 1; 
    else
      return 2;
  } else{
    if( metSig > 4 && ptyy > 25 )
      return 3;
    else
      if( ptyy > 15 )
	return 4;
  }

  return 0;
}

void getGraphFromHist(TH2F* hist , TGraph* graph, TString projected_axis){

  TH1D* projection ;
  if( projected_axis == "x")
    projection = hist->ProjectionX(TString(hist->GetName())+"_projx");

  if( projected_axis == "y")
    projection = hist->ProjectionX(TString(hist->GetName())+"_projy");

  for( int j = 1 ; j <= projection->GetNbinsX();++j)
    graph->SetPoint(j-1,projection->GetXaxis()->GetBinCenter(j),projection->GetBinContent(j));



}
void getGraphFromHist(TH2F* hist_up , TH2F* hist_down, TGraphAsymmErrors* graph,TGraph* nominal, TString projected_axis){

  TH1D* projection_up ;
  TH1D* projection_down ;
  if( projected_axis == "x"){
    projection_up = hist_up->ProjectionX(TString(hist_up->GetName())+"_projx");
    projection_down = hist_down->ProjectionX(TString(hist_down->GetName())+"_projx");
  }
  if( projected_axis == "y"){
    projection_up = hist_up->ProjectionX(TString(hist_up->GetName())+"_projy");
    projection_down = hist_down->ProjectionX(TString(hist_down->GetName())+"_projy");

  }
  for( int j = 1 ; j <= projection_up->GetNbinsX();++j){
    double x , y;
    nominal->GetPoint(j-1,x,y);
    graph->SetPoint(j-1,x,y);
    double err_high = projection_down->GetBinContent(j) >= projection_up->GetBinContent(j) ?  projection_down->GetBinContent(j)-y : projection_up->GetBinContent(j)-y;
    double err_low = projection_down->GetBinContent(j) <= projection_up->GetBinContent(j) ?  projection_down->GetBinContent(j)-y : projection_up->GetBinContent(j)-y;
    //    graph->SetPointError(j-1,0,0,projection_down->GetBinContent(j),projection_up->GetBinContent(j));
    graph->SetPointError(j-1,0,0,err_low,err_high);
  }


}
void plotBand(TFile* fileHist,TString canvName,TString sys_name,TString x_axis_name,TString y_axis_name){



  TGraph* nominal = new TGraph();

  getGraphFromHist((TH2F*)fileHist->Get("METSIG"),nominal,"y");
  TGraphAsymmErrors* graph = new TGraphAsymmErrors();

  TString HistName = "METSIG_"+sys_name;
  getGraphFromHist( (TH2F*)fileHist->Get(HistName+"__1up") , (TH2F*)fileHist->Get(HistName+"__1down") ,graph,nominal,"y");


  gStyle->SetPadTopMargin(0.06);
  gStyle->SetPadRightMargin (0.03);

  TCanvas *canv = new TCanvas("canv", "limits canvas", 700., 550.);
  canv->cd();
  canv->SetGridx(0);
  canv->SetGridy(0);
  TPad* t1 = new TPad("t1","t1", 0.0, 0., 1.0, 1.0);
  t1->Draw();
  t1->cd();
  //  t1->SetLogy(true);
  //t1->SetLogx(true);

  TMultiGraph *mg = new TMultiGraph();


  TGraph *grExp;
  grExp = new TGraph(graph->GetN());

  for( int i = 0 ; i < graph->GetN();++i){
    double x;
    double y;
    graph->GetPoint(i,x,y);
    grExp->SetPoint(i,x,y);

  }

  grExp->SetLineWidth(2);
  grExp->SetLineStyle(2);
  grExp->SetMarkerSize(1); //RJ
  grExp->SetMarkerStyle(24);


  graph->SetFillColor(kGreen);



  mg->Add(graph);
  mg->Add(grExp,"PL");
  mg->Draw("a3");

  mg->GetYaxis()->SetTitle("Systematics (%)");
  mg->GetYaxis()->SetTitleSize(0.055);
  mg->GetYaxis()->SetTitleOffset(1.35);


  mg->GetXaxis()->SetTitle(x_axis_name);
  mg->GetXaxis()->SetTitleSize(0.055);
  mg->GetXaxis()->SetTitleOffset(0.8);
  mg->GetYaxis()->SetTitle(y_axis_name);
  mg->GetYaxis()->SetTitleSize(0.055);
  mg->GetYaxis()->SetTitleOffset(0.7);


  float posx1 = 0.18;
  float posx2 = 0.58;
  float posy1 = 0.65-0.04;
  float posy2 = 0.75;
  TLegend *leg = new TLegend(posx1, posy1, posx2, posy2);
  leg->SetFillColor(0);
  leg->SetFillStyle(0);
  leg->SetLineColor(0);
  leg->SetTextFont(42);
  leg->SetBorderSize(0);

  leg->AddEntry(grExp, "Nominal value", "PL");
  leg->AddEntry(graph, "Systematic variation", "F");


  TLatex latex2;
  latex2.SetTextSize(0.04);
  latex2.SetNDC();
  latex2.DrawLatex(0.20+0.4, 0.74, "#sqrt{s} = 13 TeV, 10.1 fb^{-1}");



  latex2.DrawLatex(0.20+0.4, 0.8, "#bf{#it{ATLAS}} Internal");
  latex2.DrawLatex(0.20+0.4, 0.77,sys_name);

  leg->Draw("SAME");
  canv->SaveAs(m_outputFolder+canvName+".pdf");
  canv->SaveAs(m_outputFolder+canvName+".png");
  delete canv;
}

void plotEachSyst(TFile* fileHist,TString canvName,std::vector<TString> syslist){

  int sizeHist = 0;
  for( unsigned int i = 0 ; i < syslist.size() ;++i){

      if( !syslist[i].Contains("_1down") )
	continue;

      ++sizeHist;
  }


  TH1D* graph_up = new TH1D("graph_up","Systematic Variation",sizeHist,0,sizeHist);
  TH1D* graph = new TH1D("graph","Systematic Variation",sizeHist,0,sizeHist);
  TH1D* graph_down = new TH1D("graph_down","Systematic Variation",sizeHist,0,sizeHist);
  int point = 0 ;
  double nominal = ((TH2F*)fileHist->Get("METSIG"))->Integral();

  for( unsigned int i = 0 ; i < syslist.size() ;++i){

    TString name_stem;
      if( syslist[i].Contains("_1down"))
	name_stem = syslist[i].Remove(syslist[i].Index("__1down"));
      else
	continue;

      ++point;
      double up = ((TH2F*)fileHist->Get("METSIG_"+name_stem+"__1up"))->Integral();
      double down = ((TH2F*)fileHist->Get("METSIG_"+name_stem+"__1down"))->Integral();
      double err_high = (up - nominal)/nominal;
      double err_low = (down - nominal)/nominal;
      graph->SetBinContent(point,0);
      graph_up->SetBinContent(point,err_high*100);
      graph_down->SetBinContent(point,err_low*100);

      graph->GetXaxis()->SetBinLabel(point,name_stem);
      graph_up->GetXaxis()->SetBinLabel(point,name_stem);
      graph_down->GetXaxis()->SetBinLabel(point,name_stem);

  }

  gStyle->SetPadTopMargin(0.06);
  gStyle->SetPadRightMargin (0.03);

  TCanvas *canv = new TCanvas("canv", "limits canvas", 700., 550.);
  canv->cd();
  canv->SetGridx(0);
  canv->SetGridy(0);
  TPad* t1 = new TPad("t1","t1", 0.0, 0., 1.0, 1.0);
  t1->Draw();
  t1->cd();
  t1->SetLogy(true);
  //t1->SetLogx(true);


  graph->SetStats(0);

  graph->SetLineColor(kBlack);
  graph->SetLineWidth(3);

  graph_up->SetLineColor(kRed);
  graph_up->SetLineWidth(3);

  graph_down->SetLineColor(kGreen);
  graph_down->SetLineWidth(3);


  graph->GetYaxis()->SetTitle("");
  graph->GetYaxis()->SetTitle("Systematics (%)");
  graph->GetYaxis()->SetTitleSize(0.055);
  graph->GetYaxis()->SetTitleOffset(1.35);

  // graph->GetXaxis()->SetTitleSize(0.055);
  // graph->GetXaxis()->SetTitleOffset(0.8);
  // graph->GetYaxis()->SetTitleSize(0.055);
  // graph->GetYaxis()->SetTitleOffset(0.7);

  graph->SetMaximum(10);
  graph->SetMinimum(-10);


  canv->cd();
  graph->Draw();
  graph_up->Draw("SAME");
  graph_down->Draw("SAME");




  float posx1 = 0.18;
  float posx2 = 0.58;
  float posy1 = 0.65-0.04;
  float posy2 = 0.75;

  TLegend *leg = new TLegend(posx1, posy1, posx2, posy2);
  leg->SetFillColor(0);
  leg->SetFillStyle(0);
  leg->SetLineColor(0);
  leg->SetTextFont(42);
  leg->SetBorderSize(0);

  leg->AddEntry(graph, "Nominal");
  leg->AddEntry(graph_up, "Systematic variation up (%)");
  leg->AddEntry(graph_down, "Systematic variation down (%)");


  TLatex latex2;
  latex2.SetTextSize(0.04);
  latex2.SetNDC();
  latex2.DrawLatex(0.20+0.4, 0.74, "#sqrt{s} = 13 TeV, 10.1 fb^{-1}");
  latex2.DrawLatex(0.20+0.4, 0.8, "#bf{#it{ATLAS}} Internal");
  //  latex2.DrawLatex(0.20+0.4, 0.77,sys_name);

  leg->Draw("SAME");
  canv->SaveAs(m_outputFolder+canvName+".pdf");
  canv->SaveAs(m_outputFolder+canvName+".png");
  delete canv;
}

void calculateTotalBand(TFile* fileHist,TString canvName,std::vector<TString> syslist){

  int sizeHist = 0;
  for( unsigned int i = 0 ; i < syslist.size() ;++i){

      if( !syslist[i].Contains("_1down") )
	continue;

      ++sizeHist;
  }

  std::map<std::string,TH1D*> vector_syst_up;
  std::map<std::string,TH1D*> vector_syst_down;

  TH1D* nominal =  (TH1D*)((TH2F*)fileHist->Get("METSIG"))->ProjectionX("METSIG_projx");

  for( unsigned int i = 0 ; i < syslist.size() ;++i){

    std::string name_stem;
    if( syslist[i].Contains("_1down")){
      name_stem = syslist[i].Data();
      name_stem = (TString(name_stem).Remove(TString(name_stem).Index("__1down"))).Data();
    }else
      continue;
    // TH1D* up = ((TH2F*)fileHist->Get("METSIG_"+name_stem+"__1up"))->ProjectionX("METSIG_"+name_stem+"__1up_projx");
    // TH1D* down = ((TH2F*)fileHist->Get("METSIG_"+name_stem+"__1down"))->ProjectionX("METSIG_"+name_stem+"__1down_projx");
    vector_syst_up[name_stem] =  (TH1D*) ((TH2F*)fileHist->Get(TString("METSIG_"+name_stem+"__1up")))->ProjectionX(TString("METSIG_"+name_stem+"__1up_projx"));
    vector_syst_down[name_stem] = (TH1D*) ((TH2F*)fileHist->Get(TString("METSIG_"+name_stem+"__1down")))->ProjectionX(TString("METSIG_"+name_stem+"__1down_projx"));
  }

  TGraphAsymmErrors* graphAsymm = new TGraphAsymmErrors();
  TGraph* nominal_graph = new TGraph();


  for( int i = 1 ; i <= nominal->GetNbinsX() ;++i){
    double err2_up = 0 ;
    double err2_down = 0;

    // for( auto itMap : vector_syst_up)
    // 	err2_up += nominal->GetBinContent(i) !=0  ? TMath::Power((itMap.second->GetBinContent(i)-nominal->GetBinContent(i))/nominal->GetBinContent(i),2) : 0 ;

    // for( auto itMap : vector_syst_down)
    //   err2_down += nominal->GetBinContent(i) !=0  ? TMath::Power((itMap.second->GetBinContent(i)-nominal->GetBinContent(i))/nominal->GetBinContent(i),2) : 0 ;

    for( auto itMap : vector_syst_up){
      double upVar = itMap.second->GetBinContent(i);
      double downVar = vector_syst_down[itMap.first]->GetBinContent(i);
      double nom = nominal->GetBinContent(i);

      if( nom != 0 ){
	if( upVar > nom && downVar > nom )
	  err2_up += upVar >= downVar ? pow((upVar-nom)/nom,2) : pow((downVar-nom)/nom ,2);
	else if( upVar < nom && downVar < nom )
	  err2_down += upVar >= downVar ? pow((downVar-nom)/nom,2) : pow((upVar-nom)/nom,2) ;
	else if( upVar < nom && downVar > nom ){
	  err2_down += pow((upVar-nom)/nom,2) ;
	  err2_up += pow((downVar-nom)/nom,2) ;
	}else if( upVar > nom && downVar < nom ){
	  err2_up += pow((upVar-nom)/nom,2) ;
	  err2_down += pow((downVar-nom)/nom,2) ;
	}else{

	}
      }

    }


    nominal_graph->SetPoint(i-1,nominal->GetXaxis()->GetBinCenter(i),0);
    graphAsymm->SetPoint(i-1,nominal->GetXaxis()->GetBinCenter(i),0);
    graphAsymm->SetPointError(i-1,0,0,TMath::Sqrt(err2_down)*100,TMath::Sqrt(err2_up)*100);
    std::cout << " X is : " << nominal->GetXaxis()->GetBinCenter(i) << " Error down " << (-1)*TMath::Sqrt(err2_down)*100 << std::endl;
  }

  // for( int i = 0 ; i < graphAsymm->GetN() ; ++i ){
  //   double errlx = graphAsymm->GetErrorXlow(i);
  //   double errly = graphAsymm->GetErrorYlow(i);
  //   double errhx = graphAsymm->GetErrorXhigh(i);
  //   double errhy = graphAsymm->GetErrorYhigh(i);
  //   std::cout << "X " << graphAsymm->GetXaxis()->GetBinCenter(i) << " Error Xlow : " << errlx  << " Error Xhigh : " << errhx  << " Error Ylow : " << errly  << " Error Yhigh : " << errhy <<std::endl;


  // }



  gStyle->SetPadTopMargin(0.06);
  gStyle->SetPadRightMargin (0.03);

  TCanvas *canv = new TCanvas("canv", "limits canvas", 700., 550.);
  canv->cd();
  canv->SetGridx(0);
  canv->SetGridy(0);
  TPad* t1 = new TPad("t1","t1", 0.0, 0., 1.0, 1.0);
  t1->Draw();
  t1->cd();
  t1->SetLogy(true);
  //t1->SetLogx(true);

  //  nominal_graph->SetStats(0);

  nominal_graph->SetLineColor(kBlack);
  //  nominal_graph->SetLineWidth(3);

  graphAsymm->SetFillColor(kOrange+1);

  nominal_graph->GetXaxis()->SetTitle("E^{miss}_{T} Significance [#sqrt{GeV}]");
  nominal_graph->GetYaxis()->SetTitle("Systematics (%)");
  nominal_graph->GetYaxis()->SetTitleSize(0.055);
  nominal_graph->GetYaxis()->SetTitleOffset(1.35);

  // nominal_graph->GetXaxis()->SetTitleSize(0.055);
  // nominal_graph->GetXaxis()->SetTitleOffset(0.8);
  // nominal_graph->GetYaxis()->SetTitleSize(0.055);
  // nominal_graph->GetYaxis()->SetTitleOffset(0.7);

  nominal_graph->SetMaximum(20);
  nominal_graph->SetMinimum(-20);
  graphAsymm->SetMaximum(20);
  graphAsymm->SetMinimum(-20);


  canv->cd();

  TMultiGraph *mg = new TMultiGraph();
  mg->Add(graphAsymm);
  mg->Add(nominal_graph,"PL");

  //  nominal_graph->Draw("AC*");
  //  graphAsymm->Draw("ACF");
  mg->Draw("a3");

  float posx1 = 0.18;
  float posx2 = 0.58;
  float posy1 = 0.75;
  float posy2 = 0.9;

  TLegend *leg = new TLegend(posx1, posy1, posx2, posy2);
  leg->SetFillColor(0);
  leg->SetFillStyle(0);
  leg->SetLineColor(0);
  leg->SetTextFont(42);
  leg->SetBorderSize(0);

  leg->AddEntry(nominal_graph, "Nominal","PL");
  leg->AddEntry(graphAsymm, "Systematic variation (%)","F");

  TLatex latex2;
  latex2.SetTextSize(0.04);
  latex2.SetNDC();
  latex2.DrawLatex(0.20+0.4, 0.74, "#sqrt{s} = 13 TeV");
  latex2.DrawLatex(0.20+0.4, 0.8, "#bf{#it{ATLAS}} Internal");

  leg->Draw("SAME");

  canv->SaveAs(TString(m_outputFolder+canvName+"_TotalBand.pdf"),"RECREATE");
  canv->SaveAs(TString(m_outputFolder+canvName+"_TotalBand.png"),"RECREATE");

  TFile * file = TFile::Open(TString(m_outputFolder+canvName+"_SystBand.root"),"RECREATE");
  file->cd();
  graphAsymm->SetName("SysBand");
  graphAsymm->Write();
  file->Close();
  fileHist->cd();

  delete canv;
}


int main(int argc, char **argv)
{

  if( argc < 2 ){
    std::cout << "Input arguments are zero. This algorithm works with two inputs: one the output folder and the second, the filename containing histograms for systematics. Exiting." << std::endl;
    return 0;
  }
  std::vector<TString> syslist;
  TString sysListDir = "/afs/cern.ch/work/a/alopezso/private/DMAnalysis20.7/DMHiggsAnalysis/data/SysList_data16_h013.txt";
  std::ifstream sysListFile (sysListDir.Data());
  syslist.clear();
  syslist.push_back(""); //norminal

  if (sysListFile.is_open()) {
    std::string line;
    while ( getline (sysListFile,line) ) {
      if(line == "")continue;
      TString tempSys = line;
      syslist.push_back(tempSys);
    }
  } else {
    std::cout << "Unable to open systematic list: " << sysListDir <<std::endl;
    abort();
  }

  sysListFile.close();

  std::string argStr=argv[1];
  std::vector<std::string> fileList;
  for (size_t argI=0,n; argI <= argStr.length(); argI=n+1){
    n = argStr.find_first_of(',',argI);
    if (n == std::string::npos)
      n = argStr.length();
    string tmp = argStr.substr(argI,n-argI);
    fileList.push_back(tmp);
  }
  std::cout << TString(fileList[0]) << std::endl;
  m_outputFolder = TString(fileList[0]).EndsWith("/") ? fileList[0] : fileList[0]+"/" ;

  TString hadd = "hadd -f tmp_MxAOD_Sys.root ";
  for(unsigned int i = 1 ; i < fileList.size() ; ++i )    
    hadd += TString("  "+fileList[i]);
  
  std::cout << hadd << std::endl;
  system(hadd.Data());
  TFile* fileHist;
  if( fileList.size() == 2)
    fileHist = TFile::Open(TString(fileList[1]),"READ");
  else
    fileHist = TFile::Open(TString("tmp_MxAOD_Sys.root"),"READ");

  std::string filename_cut = fileList[1].substr(fileList[1].find("mc15c.")+6);

  if(filename_cut == "") {
    std::cout << " Filename : " << fileList[1] << " has not a correct name. Exiting" << std::endl;
    return 0;
  }

  filename_cut.replace(filename_cut.find(".SYSHIST"),-1,"");
  // std::string tag = fileList[1].substr(fileList[1].find("SYSHIST")+8);
  // tag.replace(tag.find(".root"),-1,"");
  // tag = tag != "" ? "."+tag : "";

  for( unsigned int i = 0 ; i < syslist.size() ;++i){

    TString name_stem;
    if( syslist[i].Contains("_1down")){
      name_stem = syslist[i];
      name_stem = name_stem.Remove(name_stem.Index("__1down"));
    }else
	continue;
    //      plotBand(fileHist,TString(filename_cut)+tag+"_"+name_stem,name_stem,"E^{miss}_{T} Significance [#sqrt{GeV}]","Weights");
      plotBand(fileHist,TString(filename_cut)+"_"+name_stem,name_stem,"E^{miss}_{T} Significance [#sqrt{GeV}]","Weights");
    
  }
  //  plotEachSyst(fileHist,TString(filename_cut)+tag+"_AllSyst",syslist);
  plotEachSyst(fileHist,TString(filename_cut)+"_AllSyst",syslist);
  //  calculateTotalBand(fileHist,TString(filename_cut)+tag,syslist);
  calculateTotalBand(fileHist,TString(filename_cut),syslist);
    return 0;
}
