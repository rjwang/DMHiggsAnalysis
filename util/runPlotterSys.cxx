#include <iostream>
#include <dirent.h>
#include <stdio.h>

#include "HGamAnalysisFramework/RunUtils.h"
#include "HGamAnalysisFramework/Config.h"
#include "DMHiggsAnalysis/DataEvtSummaryHandler.h"
#include "DMHiggsAnalysis/PhysicsEvent.h"

#include "TSystem.h"
#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TH1.h"
#include "TObject.h"
#include "TObjArray.h"
#include "TList.h"
#include "TH2F.h"
#include "TGraphAsymmErrors.h"
#include "TGraph.h"
#include "TMultiGraph.h"
#include "TProfile.h"
#include "TEventList.h"
#include "TROOT.h"
#include "TMath.h"
#include "THStack.h"

#include "DMHiggsAnalysis/SmartSelectionMonitor.h"
#include "DMHiggsAnalysis/deltaR.h"

#include <vector>
#include <map>

#include "./JSONWrapper.cc"

std::string m_outputFolder;

using namespace std;


bool compareMax(TH1* hist1, TH1* hist2){

  return (hist1->Integral() < hist2->Integral());
}

void showOptions(){



  std::cout << "====================================================================== " << std::endl;
  std::cout << "================ Options for runPlotterSys.cxx : ===================== " << std::endl;
  std::cout << "====================================================================== " << std::endl;
  std::cout << " " << std::endl;
  std::cout << " " << std::endl;
  std::cout << " " << std::endl;
  std::cout << " -o  output folder name " << std::endl;
  std::cout << " -j  jsonFile containing the layout for each file histogramoutput folder name " << std::endl;
  std::cout << " -f  input files that will be stacked " << std::endl; 
  std::cout << " -h  histogram name " << std::endl;

}


std::map<std::string,std::vector<int> > configureHistogramsLayout(JSONWrapper::Object jsonFile ){

  std::map<std::string,std::vector<int> > tmp_map;

  std::vector<JSONWrapper::Object> objects = jsonFile["proc"].daughters();

  for( unsigned int i = 0 ; i < objects.size() ; ++i ){

    std::vector<int> layout;
    if( objects[i].isTag("color") )
      layout.push_back(objects[i]["color"].toInt());
    else
      layout.push_back(2);
    if( objects[i].isTag("marker") )
      layout.push_back(objects[i]["marker"].toInt());
    else
      layout.push_back(2);
    if( objects[i].isTag("lcolor") )
      layout.push_back(objects[i]["lcolor"].toInt());
    else
      layout.push_back(2);

    tmp_map[objects[i]["tag"].toString()] = layout;
  }

  return tmp_map;

}
void addHistogram(std::map<std::string,TH1*>& mapHist,TH1* hist,TString name,JSONWrapper::Object jsonFile,std::string tagName=""){

  std::vector<JSONWrapper::Object> objects = jsonFile["proc"].daughters();
  for( unsigned int i = 0 ; i < objects.size() ; ++i ){
    if( tagName != "" && tagName !=  objects[i]["tag"].toString() )
      continue;
    std::vector<JSONWrapper::Object> samples = objects[i]["data"].daughters();
    for( unsigned int j = 0 ; j < samples.size() ; ++j){

      if( !TString(samples[j]["dtag"].toString()).Contains(name) )
	continue;
      if( name.Contains("_filt2") && objects[i]["tag"].toString() == "#gamma#gamma" )
	continue;


      if( mapHist.find(objects[i]["tag"].toString()) == mapHist.end() )
	mapHist[objects[i]["tag"].toString()] = (TH1*) hist->Clone(TString(objects[i]["tag"].toString()));
      else
	mapHist[objects[i]["tag"].toString()]->Add( hist );
    }
  }

}
void addHistStacker(THStack* stack, std::map<std::string,TH1*> mapHists,std::map<std::string,std::vector<int> > mapLayout){


  std::vector<TH1*> orderTH1;

  for( auto hist : mapHists ){ 
    
    hist.second->SetLineColor((int)mapLayout[hist.first][2]);
    hist.second->SetFillColor((int)mapLayout[hist.first][0]);
    hist.second->SetFillStyle(1001);
    hist.second->SetMarkerColor((int)mapLayout[hist.first][0]);
    hist.second->SetMarkerStyle((int)mapLayout[hist.first][1]);
    hist.second->SetStats(0);
    orderTH1.push_back(hist.second);
  }


  std::sort(orderTH1.begin() , orderTH1.end() , compareMax );

  for( TH1* hist : orderTH1 ) 
    stack->Add(hist);






}
void addSystematicsGraph(TGraphAsymmErrors* graph,TH1* hist,TString name){


  const char* folder = "/afs/cern.ch/work/a/alopezso/private/DMAnalysis20.7/Plots/";
  std::vector<TString> directories;


  DIR *dir = opendir(folder);

  struct dirent *entry = readdir(dir);

  while (entry != NULL)
    {

      directories.push_back(TString(entry->d_name));
      entry = readdir(dir);

    }

  closedir(dir);

  if( !hist )
    return;
  TString corresp_dir="";
  for( TString dirName : directories ){
    std::string dir_string = dirName.Data();
      if( dir_string == "." || dir_string == ".."  )
	continue;
    if( TString(hist->GetTitle()).Contains(dirName) ){
      corresp_dir = dirName.EndsWith("/") ? dirName : dirName + "/";
      break;
    }
  }


  std::string check = corresp_dir.Data();
  if( check == "" || check == "." ) 
    return;

  TString filename=TString(corresp_dir);
  filename.Remove(corresp_dir.Sizeof()-2);
  filename += "_SystBand.root";

  TFile* fileSyst= TFile::Open(folder + corresp_dir + filename,"READ");
  TGraphAsymmErrors* tmp_graph = (TGraphAsymmErrors*)fileSyst->Get("SysBand");


  double xvalue = 0 ;
  double yvalue = 0 ;
  double errxlow = 0 ;
  double errxhigh = 0 ;
  double errylow = 0 ;
  double erryhigh = 0 ;

  for( int i = 0 ; i < tmp_graph->GetN(); ++i ){

    tmp_graph->GetPoint(i,xvalue,yvalue);
    graph->SetPoint(i,xvalue,1);

  }
  for( int i = 0 ; i < tmp_graph->GetN(); ++i ){

    graph->GetPoint(i,xvalue,yvalue);
    errxlow=graph->GetErrorXlow(i);
    errxhigh=graph->GetErrorXhigh(i);
    errylow=graph->GetErrorYlow(i);
    erryhigh=graph->GetErrorYhigh(i);

    double bin_content = hist->GetBinContent(hist->GetXaxis()->FindBin(xvalue));
    errylow += xvalue <= hist->GetXaxis()->GetXmax() ?TMath::Power( bin_content*tmp_graph->GetErrorYlow(i),2 ) : 0;
    erryhigh += xvalue <= hist->GetXaxis()->GetXmax() ?TMath::Power( bin_content*tmp_graph->GetErrorYhigh(i),2 ) : 0 ;

    graph->SetPoint(i,xvalue,1);
    graph->SetPointError(i,errxlow,errxhigh,errylow,erryhigh);

  }
  fileSyst->Close();
  
  if( graph->GetN() == 0  )
    return;


  for( int i = 0 ; i < graph->GetN(); ++i ){

    double xvalue = 0 ;
    double yvalue = 0 ;
    graph->GetPoint(i,xvalue,yvalue);
    double errxlow=graph->GetErrorXlow(i);
    double errxhigh=graph->GetErrorXhigh(i);
    double errylow=graph->GetErrorYlow(i);
    double erryhigh=graph->GetErrorYhigh(i);

    double bin_content = hist->GetBinContent(hist->GetXaxis()->FindBin(xvalue));
    errylow = xvalue <= hist->GetXaxis()->GetXmax() ? TMath::Sqrt( errylow + pow(hist->GetBinError( hist->GetXaxis()->FindBin(xvalue) ),2) ) : 0 ;
    erryhigh = xvalue <= hist->GetXaxis()->GetXmax() ? TMath::Sqrt( erryhigh + pow(hist->GetBinError( hist->GetXaxis()->FindBin(xvalue) ),2) ) : 0 ;

    graph->SetPoint(i,xvalue,1);
    graph->SetPointError(i,errxlow,errxhigh,errylow,erryhigh);

  }

}
void storeSystematicsGraph(std::vector<std::string> inputFiles,std::map<std::string,TH1*> histMap,TString histogramName,JSONWrapper::Object jsonFile){


  const char* folder = "/afs/cern.ch/work/a/alopezso/private/DMAnalysis20.7/Plots/";
  std::vector<TString> directories;


  DIR *dir = opendir(folder);

  struct dirent *entry = readdir(dir);

  while (entry != NULL)
    {

      directories.push_back(TString(entry->d_name));
      entry = readdir(dir);

    }

  closedir(dir);

  std::map<std::string,TGraphAsymmErrors*> errorsMap;

  std::vector<JSONWrapper::Object> objects = jsonFile["proc"].daughters();
  for( unsigned int i = 0 ; i < objects.size() ; ++i )
    errorsMap[objects[i]["tag"].toString()] = new TGraphAsymmErrors();


  for( unsigned int j = 0 ; j < inputFiles.size(); ++j ){
    std::string inputFile = inputFiles[j];
    if( inputFile == "" || inputFile == ".")
      continue;

    TFile* fileHist = TFile::Open(TString(inputFile),"READ");
    std::string filename_cut = TString(inputFile).Contains("mc") ? inputFile.substr(inputFile.find("mc15c.")+6) : inputFile.substr(inputFile.find("13TeV.")+6);

    if( filename_cut == "" || filename_cut == "." )
      continue;
    if( TString(inputFile).Contains("data16_13TeV") || TString(inputFile).Contains("data15_13TeV") ) 
      continue;

    filename_cut.replace(filename_cut.find(".root"),-1,"");

    TH1* histogram = (TH1*) fileHist->Get(histogramName);
    histogram->SetName(histogramName+"_"+TString(filename_cut) );
    histogram->SetTitle(histogramName+"_"+TString(filename_cut) );


    for( unsigned int i = 0 ; i < objects.size() ; ++i ){

      std::vector<JSONWrapper::Object> samples = objects[i]["data"].daughters();
      for( unsigned int j = 0 ; j < samples.size() ; ++j){

	if( !TString(samples[j]["dtag"].toString()).Contains(filename_cut) )
	  continue;
	if( TString(filename_cut).Contains("_filt2") && objects[i]["tag"].toString() == "#gamma#gamma" )
	  continue;

	addSystematicsGraph(errorsMap[objects[i]["tag"].toString()],histogram,filename_cut);
      }
    }

  }


  for( auto itMap : errorsMap){
    std::cout << itMap.first << std::endl;
    // if( itMap.first == "" || itMap.first == "Data" )
    //   continue;
    TH1* sumMC = (TH1*) histMap[itMap.first];
    if( !sumMC )
      continue;
    for( int i = 0 ; i < itMap.second->GetN(); ++i ){
      double xvalue = 0 ;
      double yvalue = 0 ;
      itMap.second->GetPoint(i,xvalue,yvalue);
      itMap.second->SetPoint(i,xvalue,sumMC->GetBinContent( sumMC->GetXaxis()->FindBin(xvalue) ));
      double errxlow=itMap.second->GetErrorXlow(i);
      double errxhigh=itMap.second->GetErrorXhigh(i);
      double errylow = sumMC->GetBinContent( sumMC->GetXaxis()->FindBin(xvalue) ) !=0 ?  itMap.second->GetErrorYlow(i)/sumMC->GetBinContent( sumMC->GetXaxis()->FindBin(xvalue) ) : 0 ;
      double erryhigh = sumMC->GetBinContent( sumMC->GetXaxis()->FindBin(xvalue) ) !=0 ?  itMap.second->GetErrorYhigh(i)/sumMC->GetBinContent( sumMC->GetXaxis()->FindBin(xvalue) ) : 0 ;
      itMap.second->SetPoint(i,xvalue,0);
      itMap.second->SetPointError(i,errxlow,errxhigh,errylow,erryhigh);
    }

    std::cout << " Problem in file declaration " << std::endl;

    TString filename = itMap.first; 

    if( itMap.first == "#gamma#gamma" )
      filename = "gg";
    if( itMap.first == "V#gamma#gamma" )
      filename = "Vgg";
    if( itMap.first == "V#gamma" )
      filename = "Vg";
    if( itMap.first == "#gamma+jets" )
      filename = "gj";

    TFile* fileError = TFile::Open(TString("Systematics_"+filename+".root"),"RECREATE");
    fileError->cd();
    itMap.second->SetName("SysBand");
    itMap.second->Write();
    fileError->Close();
    std::cout << " Problem outside file declaration " << std::endl;

  }




  

}
void plotStacker(TH1* data,THStack* stack,std::map<std::string,TH1*> dmHists,TGraphAsymmErrors* errorGraph,TCanvas* canvas, std::map<std::string,std::vector<int> > mapLayout){

  TH1* sumMC = (TH1*)stack->GetStack()->Last();

  TGraphAsymmErrors* errorGraphRatio = new TGraphAsymmErrors();

  for( int i = 0 ; i < errorGraph->GetN(); ++i ){
    double xvalue = 0 ;
    double yvalue = 0 ;
    errorGraph->GetPoint(i,xvalue,yvalue);
    errorGraph->SetPoint(i,xvalue,sumMC->GetBinContent( sumMC->GetXaxis()->FindBin(xvalue) ));
    double errxlow=errorGraph->GetErrorXlow(i);
    double errxhigh=errorGraph->GetErrorXhigh(i);
    double errylow = sumMC->GetBinContent( sumMC->GetXaxis()->FindBin(xvalue) ) !=0 ?  errorGraph->GetErrorYlow(i)/sumMC->GetBinContent( sumMC->GetXaxis()->FindBin(xvalue) ) : 0 ;
    double erryhigh = sumMC->GetBinContent( sumMC->GetXaxis()->FindBin(xvalue) ) !=0 ?  errorGraph->GetErrorYhigh(i)/sumMC->GetBinContent( sumMC->GetXaxis()->FindBin(xvalue) ) : 0 ;
    errorGraphRatio->SetPoint(i,xvalue,1);
    errorGraphRatio->SetPointError(i,errxlow,errxhigh,errylow,erryhigh);
  }

  canvas->Clear();

  TPad* pad_Up = new TPad("pad_Up","pad_Up",0.1,0.35,0.9,1,0);
  TPad* pad_Down = new TPad("pad_Down","pad_Down",0.1,0.1,0.9,0.35,0);

  pad_Up->SetBottomMargin(0);
  pad_Up->SetBorderMode(0);
  pad_Down->SetTopMargin(0);
  pad_Down->SetBottomMargin( 0.3 );
  pad_Up->Clear();
  pad_Down->Clear();
  
  canvas->cd();
  pad_Up->Draw();
  pad_Up->cd();
  gPad->SetTicks(1,1);

  data->SetStats(0);

  gPad->SetLogy(1);
  data->GetYaxis()->SetRangeUser(0.001,data->GetBinContent(data->GetMaximumBin())*1000);
  stack->SetMaximum(data->GetBinContent(data->GetMaximumBin())*1000);
  stack->SetMinimum(0.001);

  if( !stack ){std::cout << " mc Histogram is not filled  " << std::endl;exit(-1);}


  if( data )data->SetStats(0);
   if( stack )
     stack->Draw("hist");
  else
    std::cout << "No stack histogram filled. Please check the code." << std::endl;
  if( data )
    data->Draw("SAME");
  else
    std::cout << "No data histogram filled. Please check the code." << std::endl;
  
  for( auto hist : dmHists){
    hist.second->SetLineColor(mapLayout[hist.first][0]);
    hist.second->SetLineWidth(3);
  }
  

  for(auto hist : dmHists){
    TH1* hist2 = hist.second;
    hist2->Draw("hist""SAME");
  }


  errorGraph->SetFillStyle(3002);
  errorGraph->Draw("3""SAME");

  TLegend* leg = new TLegend(0.45,0.6,0.9,0.85);
  leg->SetFillStyle(0);
  leg->SetNColumns(3);
  leg->SetTextSize(0.032);
  leg->SetBorderSize(0);
  leg->AddEntry(data,"Data","P");
  leg->AddEntry(errorGraph,"Stat #oplus syst");

  for( TObject* obj : *(stack->GetStack()) ){
    TH1* hist = (TH1*) obj;
    if( !hist )
      continue;

    hist->SetLineColor((int)mapLayout[hist->GetName()][2]);
    hist->SetFillColor((int)mapLayout[hist->GetName()][0]);
    hist->SetFillStyle(1001);
    hist->SetMarkerStyle((int)mapLayout[hist->GetName()][1]);
    hist->SetMarkerColor((int)mapLayout[hist->GetName()][0]);
    
    leg->AddEntry(hist,hist->GetName(),"f");
  }

  for( auto hist : dmHists)
    leg->AddEntry(hist.second,TString(hist.first),"L");

  if( !leg ){std::cout << " Legend is not filled  " << std::endl;exit(-1);}

  leg->Draw("SAME");

  TLatex latex;
  latex.SetTextSize(0.04);
  latex.SetNDC();
  latex.SetTextAlign(22);
  latex.SetTextFont(42);
  latex.SetTextSizePixels(22);


  latex.DrawLatex(0.2,0.85,"#bf{#it{ATLAS}} Internal");
  latex.DrawLatex(0.2,0.75,"#sqrt{s} = 13TeV, #int #it{L}dt = 15.4 fb^{-1}");

  stack->GetXaxis()->SetTitle(data->GetXaxis()->GetTitle());
  stack->GetYaxis()->SetTitle(data->GetYaxis()->GetTitle());
  stack->SetTitle(data->GetTitle());
  stack->GetXaxis()->SetTitleOffset(0.8);
  stack->GetYaxis()->SetTitleOffset(0.8);
  stack->GetXaxis()->SetTitleSize(0.05);
  stack->GetYaxis()->SetTitleSize(0.05);

  pad_Up->Modified();
  canvas->cd();
  pad_Down->Draw();
  pad_Down->cd();
  gPad->SetTicks(1,1);


  TLine* line = new TLine( data->GetXaxis()->GetXmin(), 1 ,data->GetXaxis()->GetXmax(), 1 );
  line->SetLineColor(kBlue);
  line->SetLineWidth(2);
  line->SetLineStyle(2);

  TH1* ratio= (TH1*)data->Clone("ratio");
  ratio->SetStats(0);
  //  if(isNorm) (*mcAddition)=(*mcAddition)*(data->Integral()/(mcAddition->Integral()));
  ratio->Divide(sumMC);
  ratio->GetXaxis()->SetTitleSize(0.15);
  ratio->GetXaxis()->SetTitleOffset(0.8);
  ratio->GetYaxis()->SetLabelSize(0.09);
  ratio->GetYaxis()->SetTitleSize(0.15);
  ratio->GetYaxis()->SetTitleOffset(0.3);
  ratio->GetYaxis()->SetTitle("Data/MC");
  ratio->GetXaxis()->SetLabelSize(0.15);
  ratio->GetYaxis()->SetRangeUser(0.,2.);
  errorGraphRatio->GetXaxis()->SetTitleSize(0.15);
  errorGraphRatio->GetXaxis()->SetTitleOffset(0.8);
  errorGraphRatio->GetYaxis()->SetLabelSize(0.09);
  errorGraphRatio->GetYaxis()->SetTitleSize(0.15);
  errorGraphRatio->GetYaxis()->SetTitleOffset(0.3);
  errorGraphRatio->GetYaxis()->SetTitle("Data/MC");
  errorGraphRatio->GetXaxis()->SetLabelSize(0.15);
  errorGraphRatio->GetYaxis()->SetRangeUser(0.,2.);
  errorGraphRatio->GetYaxis()->SetRangeUser(0.,2.);
  errorGraphRatio->GetXaxis()->SetRangeUser(ratio->GetXaxis()->GetXmin(),ratio->GetXaxis()->GetXmax());
  errorGraphRatio->SetFillColor(kOrange+2); 
  errorGraphRatio->Draw("a3");
  ratio->Draw("ep""SAME");



  ratio->SetTitle("");
  gPad->SetGridy();
  gStyle->SetAxisColor(kBlack);
  line->Draw();
  gPad->Modified();

}



int main(int argc, char **argv)
{

  if( argc < 2 ){
    std::cout << "Input arguments are zero. This algorithm works with two inputs: one the output folder and the second, the filename containing histograms for systematics. Exiting." << std::endl;
    return 0;
  }



  std::vector<std::string> fileList;
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
      fileList.push_back(tmp);
    }
    std::string option = argv[i]; 
    options[option.substr(1,1)]=fileList;
    fileList.clear();
  }

  std::string outputFolder;
  std::vector<std::string> inputFiles;
  std::vector<TString> histograms;
  std::vector<std::string> histogramStem;
  std::string jsonFile;

  for( auto itMap : options){

    if( itMap.first == "o" )
      outputFolder = itMap.second[0];
    else if( itMap.first == "j" )
      jsonFile = itMap.second[0];
    else if( itMap.first == "f" )
      inputFiles = itMap.second;
    else if( itMap.first == "h" )
      histogramStem = itMap.second;
    else 
      showOptions();
  }

  m_outputFolder = TString(outputFolder).EndsWith("/") ? outputFolder : outputFolder+"/" ;

  JSONWrapper::Object json(jsonFile,true);
  std::map<std::string,std::vector<int> > histogram_colors = configureHistogramsLayout( json );
  TFile* fileTest = TFile::Open(TString(inputFiles[0]),"READ");
  TList* list = fileTest->GetListOfKeys();
  
  for( TObject* obj : *list ) {
    TH1* tmp_hist = (TH1*) obj;
    if( !tmp_hist )
      continue;

    bool match=false;
    for( unsigned int z = 0 ; z < histogramStem.size() ; ++z ){
      std::string objName = obj->GetName(); 
      if( objName == histogramStem[z] )
	match = true;  
      continue;
    }

    if( match  )
	histograms.push_back(TString(obj->GetName()));

  }

  fileTest->Close();

  for( unsigned int j = 0 ; j < inputFiles.size(); ++j ){
    std::string inputFile = inputFiles[j];
    if( inputFile == "" || inputFile == ".")
      continue;
    std::cout << " Running on file : " << inputFile << std::endl;
  }


  for( unsigned int i = 0 ; i < histograms.size(); ++i ){
    TString histogramName = histograms[i];

    THStack* stackPlot = nullptr;
    TH1* data = nullptr;
    std::map<std::string,TH1*> darkMatter;
    std::map<std::string,TH1*> allHists;

    TGraphAsymmErrors* graph = new TGraphAsymmErrors();

    for( unsigned int j = 0 ; j < inputFiles.size(); ++j ){
      std::string inputFile = inputFiles[j];
      if( inputFile == "" || inputFile == ".")
	continue;

      TFile* fileHist = TFile::Open(TString(inputFile),"READ");
      std::string filename_cut = TString(inputFile).Contains("mc") ? inputFile.substr(inputFile.find("mc15c.")+6) : inputFile.substr(inputFile.find("13TeV.")+6);

      if( filename_cut == "" || filename_cut == "." )
	continue;
      if( TString(inputFile).Contains("data16_13TeV") || TString(inputFile).Contains("data15_13TeV") ) {
	

	if( data != NULL  ){	
	  data->Add( (TH1*)fileHist->Get(histogramName) );
	  data->SetDirectory(0);
	}else{
	  data = (TH1*)fileHist->Get(histogramName);
	  data->SetDirectory(0);

	}
	data->SetDirectory(0);
	fileHist->Close();
	continue;

      }

      filename_cut.replace(filename_cut.find(".root"),-1,"");

      TH1* histogram = (TH1*) fileHist->Get(histogramName);
      histogram->SetName(histogramName+"_"+TString(filename_cut) );
      histogram->SetTitle(histogramName+"_"+TString(filename_cut) );

      if( TString(inputFile).Contains("hxx") || TString(inputFile).Contains("2hdmxx") ) {
	addHistogram(darkMatter,histogram,filename_cut,json);
	continue;
      }


      addHistogram(allHists,histogram,filename_cut,json);
      addSystematicsGraph(graph,histogram,filename_cut);
    }

    if( !stackPlot )
      stackPlot = new THStack(histogramName,histogramName);
    addHistStacker(stackPlot,allHists,histogram_colors);
    TObjArray* list = stackPlot->GetStack();
    for( TObject* obj : *list){
      TH1* hist = (TH1*) obj;
      if( !hist )
    	continue;
    }
    TCanvas* canvas = new TCanvas("canvas","canvas",1000,700);
    plotStacker(data,stackPlot,darkMatter,graph,canvas,histogram_colors);

    canvas->SaveAs(m_outputFolder+histogramName+".pdf","RECREATE");
    canvas->SaveAs(m_outputFolder+histogramName+".png","RECREATE");
    storeSystematicsGraph(inputFiles,allHists,histogramName,json);
  }









  return 0;
}
