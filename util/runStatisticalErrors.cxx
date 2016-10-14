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


bool compareMax(TH1F* hist1, TH1F* hist2){

  return (hist1->Integral() > hist2->Integral());
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

    if( match )
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


    for( unsigned int j = 0 ; j < inputFiles.size(); ++j ){
      std::string inputFile = inputFiles[j];
      if( inputFile == "" || inputFile == ".")
	continue;

      TFile* fileHist = TFile::Open(TString(inputFile),"READ");
      std::string filename_cut = TString(inputFile).Contains("mc") ? inputFile.substr(inputFile.find("mc15c.")+6) : inputFile.substr(inputFile.find("13TeV.")+6);
      filename_cut = TString(inputFile).Contains("/2hdm") ?  inputFile.substr(inputFile.find("/2hdm")+1): filename_cut;
      filename_cut = TString(inputFile).Contains("/hxx") ?  inputFile.substr(inputFile.find("/hxx")+1): filename_cut;
      filename_cut.replace(filename_cut.find(".NTuple"),-1,"");

      TH1F* myhist = (TH1F*)fileHist->Get(histogramName);
      myhist->SetDirectory(0);

      TH1F* statHist = new TH1F("statHist","Statistical Errors",myhist->GetNbinsX(),myhist->GetXaxis()->GetXmin(),myhist->GetXaxis()->GetXmax());

      for( int j = 1 ; j <= myhist->GetNbinsX() ; ++j)
	statHist->Fill(myhist->GetBinCenter(j),myhist->GetBinContent(j) != 0 ? myhist->GetBinError(j)/myhist->GetBinContent(j) : 0 );
    
      TCanvas* canvas = new TCanvas("canvas","canvas");
      canvas->cd();
      statHist->SetStats(0);
      statHist->GetXaxis()->SetTitle(myhist->GetXaxis()->GetTitle());
      statHist->GetYaxis()->SetTitle("StatError/BinContent");

      statHist->SetFillColor(kBlue);
      statHist->Draw("hist""BAR");

      canvas->SaveAs(m_outputFolder+histogramName+"_"+filename_cut+"_StatisticalErrors.pdf","RECREATE");
      canvas->SaveAs(m_outputFolder+histogramName+"_"+filename_cut+"_StatisticalErrors.png","RECREATE");

      fileHist->Close();
    }

  }

  return 0;
}
