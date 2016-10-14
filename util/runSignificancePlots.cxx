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
  std::cout << " -s  signal input files" << std::endl; 
  std::cout << " -b  background input files " << std::endl; 
  std::cout << " -h  histogram name " << std::endl;

}

void getSignificancePlot2D(TH2F* Significance, TH2F* Signal, TH2F* Background ){

 
  for( int i = 1 ; i <= Signal->GetNbinsX() ; ++i ){
    for( int j = 1 ; j <= Signal->GetNbinsY() ; ++j ){
      double signalIntegral = Signal->Integral(i,Signal->GetNbinsX(),j,Signal->GetNbinsY());
      double backgroundIntegral = Background->Integral(i,Background->GetNbinsX(),j,Background->GetNbinsY());
      double significance =  backgroundIntegral > 0 ? signalIntegral/TMath::Sqrt(backgroundIntegral) : 0 ;
      //      Significance->Fill( Signal->GetXaxis()->GetBinCenter(i),Signal->GetYaxis()->GetBinCenter(j),significance);
      Significance->SetBinContent(i,j,significance);
    
      if( backgroundIntegral < 0){
	std::cout << " Background integral is lower than zero. Filling significance with 0." << std::endl;
	std::cout << " However, check your event yields. A significance study should use histgrams on which event yields are divided, ergo at most background integral can be zero" << std::endl;
      }
    }
  }
}

void getSignificancePlot1D(TH1F* Significance, TH1F* Signal, TH1F* Background ){

  for( int i = 1 ; i <= Signal->GetNbinsX() ; ++i ){
    double signalIntegral = Signal->Integral(i,Signal->GetNbinsX());
    double backgroundIntegral = Background->Integral(i,Background->GetNbinsX());
    double significance =  backgroundIntegral > 0 ? signalIntegral/TMath::Sqrt(backgroundIntegral) : 0 ;
    Significance->Fill( Signal->GetXaxis()->GetBinCenter(i), significance  );

    if( backgroundIntegral < 0){
      std::cout << " Background integral is lower than zero. Filling significance with 0." << std::endl;
      std::cout << " However, check your event yields. A significance study should use histgrams on which event yields are divided, ergo at most background integral can be zero" << std::endl;
    }
  }
  
}

void plotSignificance(TObject* signal , TObject* background, TString outCanvName ){

  TH1F* test1D_Signal =(TH1F*) signal;
  TH2F* test2D_Signal =(TH2F*) signal;
  TH1F* test1D_Back =(TH1F*) background;
  TH2F* test2D_Back =(TH2F*) background;


  TCanvas* canvas = new TCanvas("canvas","canvas",1000,700);

  if(test2D_Signal && test2D_Back)
    {

      TH2F* Signal =(TH2F*) signal->Clone("Signal");
      TH2F* Background = (TH2F*)background->Clone("Background");
      TH2F* Significance = new TH2F(TString(Signal->GetName())+"_Significance","Significance",Signal->GetNbinsX(),Signal->GetXaxis()->GetXmin(),Signal->GetXaxis()->GetXmax(),Signal->GetNbinsY(),Signal->GetYaxis()->GetXmin(),Signal->GetYaxis()->GetXmax());
      getSignificancePlot2D(Significance,Signal,Background);
      
      if( !Significance )
	std::cout << "Significance was not correctly filled. "<< std::endl;

      canvas->cd();
      Significance->SetStats(0);
      Significance->GetXaxis()->SetTitle(Signal->GetXaxis()->GetTitle());
      Significance->GetYaxis()->SetTitle(Signal->GetYaxis()->GetTitle());
      Significance->GetZaxis()->SetTitle(Signal->GetZaxis()->GetTitle());
      Significance->Draw("COLZ");
      canvas->SaveAs(outCanvName+"_Significance.pdf","RECREATE");
      canvas->SaveAs(outCanvName+"_Significance.png","RECREATE");

  }
  else if( test1D_Signal && test1D_Back)
    {
      TH1F* Signal =(TH1F*) signal->Clone("Signal");
      TH1F* Background = (TH1F*)background->Clone("Background");
      TH1F* Significance = new TH1F(TString(Signal->GetName())+"_Significance","Significance",Signal->GetNbinsX(),Signal->GetXaxis()->GetXmin(),Signal->GetXaxis()->GetXmax());
      getSignificancePlot1D(Significance,Signal,Background);
      
      if( !Significance )
	std::cout << "Significance was not correctly filled. "<< std::endl;

      canvas->cd();
      Significance->SetStats(0);
      Significance->GetXaxis()->SetTitle(Signal->GetXaxis()->GetTitle());
      Significance->GetYaxis()->SetTitle(Signal->GetYaxis()->GetTitle());
      Significance->Draw();
      canvas->SaveAs(outCanvName+"_Significance.pdf","RECREATE");
      canvas->SaveAs(outCanvName+"_Significance.png","RECREATE");
    }
  else{

    if( !test1D_Signal && test1D_Back )
      std::cout << " 1D histogram : signal is not defined " << std::endl;

    else if( test1D_Signal && !test1D_Back )
      std::cout << " 1D histogram : background is not defined " << std::endl;

    else if( !test2D_Signal && test2D_Back )
      std::cout << " 2D histogram : signal is not defined " << std::endl;

    else if( test2D_Signal && !test2D_Back )
      std::cout << " 2D histogram : background is not defined " << std::endl;
    else
      std::cout << "Unused type of histograms : " << signal->GetName() << std::endl;

  }

  delete canvas;
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
  std::vector<std::string> SignalinputFiles;
  std::vector<std::string> BackgroundinputFiles;
  std::vector<TString> histograms2D;
  std::vector<TString> histograms1D;
  std::vector<std::string> histogramStem;

  for( auto itMap : options){

    if( itMap.first == "o" )
      outputFolder = itMap.second[0];
    else if( itMap.first == "s" )
      SignalinputFiles = itMap.second;
    else if( itMap.first == "b" )
      BackgroundinputFiles = itMap.second;
    else if( itMap.first == "h" )
      histogramStem = itMap.second;
    else 
      showOptions();
  }

  m_outputFolder = TString(outputFolder).EndsWith("/") ? outputFolder : outputFolder+"/" ;

  TFile* fileTest = TFile::Open(TString(SignalinputFiles[0]),"READ");
  TList* list = fileTest->GetListOfKeys();
  

  unsigned int idenHistograms = 0 ; 

  for( unsigned int z = 0 ; z < histogramStem.size() ; ++z ){

    TH1F* hist1D = nullptr ;
    TH2F* hist2D = nullptr; 

    if( hist1D )
      hist1D->Clear();
    if( hist2D )
      hist2D->Clear(); 
    
    bool match=false;

    for( TObject* obj : *list ) {
  
      hist2D = (TH2F*) obj;
      hist1D=(TH1F*) obj;

      if( !hist1D && !hist2D )
	continue;

      std::string objName = obj->GetName();
      if( objName == histogramStem[z] )
	match = true;
      
    }

    if( match ){
      if( hist2D )
	histograms2D.push_back(TString(histogramStem[z]));
      else{
	if( hist1D )
	  histograms1D.push_back(TString(histogramStem[z]));
	else
	  continue;
      }
    }else{
      std::cout << " Histogram "+ histogramStem[z] + " is not in currentFile. Please check the correct name. Exiting " << std::endl;
      exit(-1);
    }

  }
  fileTest->Close();

  std::cout << " Histograms 2D : " << histograms2D.size() << std::endl;
  std::cout << " Histograms 1D : " << histograms1D.size() << std::endl;

  for( unsigned int i = 0 ; i < histograms2D.size(); ++i ){
    TString histogramName = histograms2D[i];

    TH2F* Background = nullptr;
    TH2F* Signal = nullptr;

    for( unsigned int j = 0 ; j < BackgroundinputFiles.size(); ++j ){
      std::string inputFile = BackgroundinputFiles[j];
      if( inputFile == "" || (inputFile == "." || inputFile == ".." ) )
	continue;

      TFile* fileHist = TFile::Open(TString(inputFile),"READ");
      if( !Background )
	Background = (TH2F*) fileHist->Get(histogramName);
      else
	Background->Add( (TH2F*) fileHist->Get(histogramName) );
      Background->SetDirectory(0);
      fileHist->Close();
    }

    for( unsigned int j = 0 ; j < SignalinputFiles.size(); ++j ){
      std::string inputFile = SignalinputFiles[j];
      if( inputFile == "" || (inputFile == "." || inputFile == ".." ))
	continue;

      TFile* fileHist = TFile::Open(TString(inputFile),"READ");

      TH2F* hist = (TH2F*) fileHist->Get(histogramName);
      std::string filename_cut = TString(inputFile).Contains("mc") ? inputFile.substr(inputFile.find("mc15c.")+6) : inputFile.substr(inputFile.find("13TeV.")+6);

      if( !hist ){
	std::cout << " Histogram " << histogramName << " does not exist in input file " << inputFile << std::endl;
	continue;
      }
      gSystem->Exec("mkdir -p "+m_outputFolder+TString(filename_cut));
      plotSignificance(hist , Background, TString(m_outputFolder+TString(filename_cut+"/")+histogramName) );

      // TH2F* Significance = new TH2F(TString(Signal->GetName())+"_Significance","Significance",Signal->GetNbinsX(),Signal->GetXaxis()->GetXmin(),Signal->GetXaxis()->GetXmax(),Signal->GetNbinsY(),Signal->GetYaxis()->GetXmin(),Signal->GetYaxis()->GetXmax());
      // getSignificancePlot2D(Significance,Signal,Background);

      // Significance->GetXaxis()->SetTitle(Signal->GetXaxis()->GetTitle());
      // Significance->GetYaxis()->SetTitle(Signal->GetYaxis()->GetTitle());
      // Significance->GetZaxis()->SetTitle(Signal->GetZaxis()->GetTitle());

      // TCanvas* canvas = new TCanvas("canvas_"+histogramName,"canvas");
      // canvas->cd();
      // Significance->SetStats(0);

      // if( !Significance )
      // 	std::cout << "Significance was not correctly filled. "<< std::endl;

      // Significance->Draw("COLZ");
      // canvas->SaveAs(m_outputFolder+TString(filename_cut+"/")+histogramName+"_Significance.pdf","RECREATE");
      // canvas->SaveAs(m_outputFolder+TString(filename_cut+"/")+histogramName+"_Significance.png","RECREATE");
      if( !Signal )
	Signal = (TH2F*) hist->Clone("AllSignals");
      else
	Signal->Add( hist );
      Signal->SetDirectory(0);
      fileHist->Close();
    }

    if( !Signal || ! Background){
      std::cout << "Signal or background histograms are empty. Please, check input files or this code for errors." << std::endl;
      exit(-1);
    }
 
    plotSignificance(Signal , Background, TString(m_outputFolder+histogramName+"AllSignals") );
    // TH2F* Significance = new TH2F(TString(Signal->GetName())+"_Significance","Significance",Signal->GetNbinsX(),Signal->GetXaxis()->GetXmin(),Signal->GetXaxis()->GetXmax(),Signal->GetNbinsY(),Signal->GetYaxis()->GetXmin(),Signal->GetYaxis()->GetXmax());
    // getSignificancePlot2D(Significance,Signal,Background);

    // Significance->GetXaxis()->SetTitle(Signal->GetXaxis()->GetTitle());
    // Significance->GetYaxis()->SetTitle(Signal->GetYaxis()->GetTitle());
    // Significance->GetZaxis()->SetTitle(Signal->GetZaxis()->GetTitle());

    // TCanvas* canvas = new TCanvas("canvas_"+histogramName,"canvas");
    // canvas->cd();
    // Significance->SetStats(0);

    // if( !Significance )
    //   std::cout << "Significance was not correctly filled. "<< std::endl;

    // Significance->Draw("COLZ");
    // canvas->SaveAs(m_outputFolder+histogramName+"AllSignals_Significance.pdf","RECREATE");
    // canvas->SaveAs(m_outputFolder+histogramName+"AllSignals_Significance.png","RECREATE");
  }

  for( unsigned int i = 0 ; i < histograms1D.size(); ++i ){
    TString histogramName = histograms1D[i];

    TH1F* Background = nullptr;
    TH1F* Signal = nullptr;


    for( unsigned int j = 0 ; j < BackgroundinputFiles.size(); ++j ){
      std::string inputFile = BackgroundinputFiles[j];
      if( inputFile == "" || inputFile == ".")
	continue;

      TFile* fileHist = TFile::Open(TString(inputFile),"READ");
      if( !Background )
	Background = (TH1F*) fileHist->Get(histogramName);
      else
	Background->Add( (TH1F*) fileHist->Get(histogramName) );

      Background->SetDirectory(0);
      fileHist->Close();
    }

    for( unsigned int j = 0 ; j < SignalinputFiles.size(); ++j ){
      std::string inputFile = SignalinputFiles[j];
      if( inputFile == "" || (inputFile == "." || inputFile == ".." ))
	continue;

      TFile* fileHist = TFile::Open(TString(inputFile),"READ");

      TH1F* hist = (TH1F*) fileHist->Get(histogramName);
      if( !hist ){
	std::cout << " Histogram " << histogramName << " does not exist in input file " << inputFile << std::endl;
	continue;
      }

      std::string filename_cut = TString(inputFile).Contains("mc") ? inputFile.substr(inputFile.find("mc15c.")+6) : inputFile.substr(inputFile.find("13TeV.")+6);

      gSystem->Exec("mkdir -p "+m_outputFolder+TString(filename_cut));
      plotSignificance(hist , Background, TString(m_outputFolder+TString(filename_cut+"/")+histogramName) );



      // TH1F* Significance = new TH1F(TString(Signal->GetName())+"_Significance","Significance",Signal->GetNbinsX(),Signal->GetXaxis()->GetXmin(),Signal->GetXaxis()->GetXmax());
      // getSignificancePlot1D(Significance,Signal,Background);

      // if( !Significance )
      // 	std::cout << "Significance was not correctly filled. "<< std::endl;


      // TCanvas* canvas = new TCanvas("canvas","canvas",1000,700);
      // canvas->cd();
      // Significance->SetStats(0);
      // Significance->GetXaxis()->SetTitle(Signal->GetXaxis()->GetTitle());
      // Significance->GetYaxis()->SetTitle(Signal->GetYaxis()->GetTitle());
      // Significance->Draw();
      // canvas->SaveAs(m_outputFolder+TString(filename_cut+"/")+histogramName+"_Significance.pdf","RECREATE");
      // canvas->SaveAs(m_outputFolder+TString(filename_cut+"/")+histogramName+"_Significance.png","RECREATE");

      if( !Signal )
	Signal = (TH1F*) hist->Clone("AllSignals");
      else
	Signal->Add( hist );

      Signal->SetDirectory(0);
      fileHist->Close();




    }

    if( !Signal || ! Background){
      std::cout << "Signal or background histograms are empty. Please, check input files or this code for errors." << std::endl;
      exit(-1);
    }
    Signal->SetName(TString(histogramName+"AllSignals"));
    plotSignificance(Signal , Background, TString(m_outputFolder+histogramName+"AllSignals") );
    // TH1F* Significance = new TH1F(TString(Signal->GetName())+"_Significance","Significance",Signal->GetNbinsX(),Signal->GetXaxis()->GetXmin(),Signal->GetXaxis()->GetXmax());
    // getSignificancePlot1D(Significance,Signal,Background);

    // if( !Significance )
    //   std::cout << "Significance was not correctly filled. "<< std::endl;


    // TCanvas* canvas = new TCanvas("canvas","canvas",1000,700);
    // canvas->cd();
    // Significance->SetStats(0);
    // Significance->GetXaxis()->SetTitle(Signal->GetXaxis()->GetTitle());
    // Significance->GetYaxis()->SetTitle(Signal->GetYaxis()->GetTitle());
    // Significance->Draw();
    // canvas->SaveAs(m_outputFolder+histogramName+"AllSignals_Significance.pdf","RECREATE");
    // canvas->SaveAs(m_outputFolder+histogramName+"AllSignals_Significance.png","RECREATE");
  }
  return 0;
}
