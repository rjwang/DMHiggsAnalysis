#include <string>
#include <vector>
#include <iostream>
#include <list>
#include <iterator>
#include <algorithm>
#include <unordered_map>
#include <sstream>

#include "TROOT.h"
#include "TFile.h"
#include "TDirectory.h"
#include "TChain.h"
#include "TObject.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TLegend.h"
#include "TGraph.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TTree.h"
#include "TF1.h"
#include "TCutG.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TMultiGraph.h"
#include "TPaveText.h"
#include "THStack.h"
#include "Rtypes.h"
#include "TList.h"
#include "TGaxis.h"

#include "TLatex.h"
#include "TPad.h"


#include "./tdrstyle.C"
#include "./JSONWrapper.cc"

int cutIndex=-1;
string cutIndexStr="";
double iLumi = 2007;
double iEcm=8;
bool showChi2 = false;
bool showUnc=false;
double baseRelUnc=0.044;
bool noLog=false;
bool isSim=false;
bool isinProgress=false;
bool isDataBlind=false;
bool do2D  = true;
bool do1D  = true;
bool doTex = true;
bool StoreInFile = true;
bool doPlot = true;
bool splitCanvas = false;
bool onlyCutIndex = false;
bool nosig = false;
bool noratio = false;
double scaleSignal=1.0;
double scaleYMax = 1.0;
double scaleYMin = 1.0;
bool useDataMinusMC = false;
string inDir   = "OUTNew/";
string jsonFile = "../../data/beauty-samples.json";
string outDir  = "Img/";
string plotExt = ".png";
string plotExt2 = ".pdf";
string outFile = "plotter.root";
string cutflowhisto = "all_cutflow";
string ctrlSample = "";
int rebin = 1; //RJ
float cutLine1, cutLine2;
float leftcutLine, rightcutLine;


struct stSampleInfo{ double PURescale_up; double PURescale_down; double initialNumberOfEvents;};
std::unordered_map<string, stSampleInfo> sampleInfoMap;
std::unordered_map<string, bool> FileExist;
TString getChannelName(std::string SaveName);
void savePath(TCanvas* c,std::string outDir,std::string SaveName,std::string plotExt);

struct NameAndType{
   std::string name;
   bool type;
   bool isIndexPlot;
   NameAndType(std::string name_,  bool type_, bool isIndexPlot_){name = name_; type = type_; isIndexPlot = isIndexPlot_;}

   bool operator==(const NameAndType& a){ return a.name == name;}
   bool operator< (const NameAndType& a){ return a.name < name;}
 };


TObject* GetObjectFromPath(TDirectory* File, std::string Path, bool GetACopy=false)
{
   size_t pos = Path.find("/");
   if(pos < 256){
      std::string firstPart = Path.substr(0,pos);
      std::string endPart   = Path.substr(pos+1,Path.length());
      TDirectory* TMP = (TDirectory*)File->Get(firstPart.c_str());
      if(TMP!=NULL){
         TObject* TMP2 =  GetObjectFromPath(TMP,endPart,GetACopy);
         return TMP2;
      }

      return NULL;
   }else{
      TObject* TMP = File->Get(Path.c_str());

      if(GetACopy){
         return TMP->Clone();
      }else{
         return TMP;
      }
   }

}

void GetListOfObject(JSONWrapper::Object& Root, std::string RootDir, std::list<NameAndType>& histlist, TDirectory* dir=NULL, std::string parentPath="")
{
  if(parentPath=="" && !dir)
    {
      int dataProcessed = 0;
      int signProcessed = 0;
      int bckgProcessed = 0;

      std::vector<JSONWrapper::Object> Process = Root["proc"].daughters();
      //loop on all Proc
      for(size_t ip=0; ip<Process.size(); ip++){
          if(Process[ip]["isinvisible"].toBool())continue;
	  bool isData (  Process[ip]["isdata"].toBool()  );
          bool isSign ( !isData &&  Process[ip].isTag("spimpose") && Process[ip]["spimpose"].toBool());
  	  bool isMC   = !isData && !isSign;
	  string filtExt("");
	  if(Process[ip].isTag("mctruthmode") ) { char buf[255]; sprintf(buf,"_filt%d",(int)Process[ip]["mctruthmode"].toInt()); filtExt += buf; }

          //just to make it faster, only consider the first 3 sample of a same kind
          if(isData){if(dataProcessed>=1){continue;}else{dataProcessed++;}}
          if(isSign){if(signProcessed>=2){continue;}else{signProcessed++;}}
          if(isMC  ){if(bckgProcessed>8) {continue;}else{bckgProcessed++;}}

	  std::vector<JSONWrapper::Object> Samples = (Process[ip])["data"].daughters();
	  for(size_t id=0; id<Samples.size(); id++){
	      int split = 1;
	      if(Samples[id].isTag("split"))split = Samples[id]["split"].toInt();
	      string segmentExt;
	      if(split>1) { char buf[255]; sprintf(buf,"_%i",0); segmentExt += buf; }
              string FileName = RootDir + (Samples[id])["dtag"].toString() +  (Samples[id].isTag("suffix")?(Samples[id])["suffix"].toString():string("")) + segmentExt + filtExt + ".root";

              printf("\033[33mAdding all objects from %25s to the list of considered objects\033[0m\n",  FileName.c_str());
	      TFile* file = new TFile(FileName.c_str());
	      if(file->IsZombie())
		{
		  if(isData) dataProcessed--;
		  if(isSign) signProcessed--;
		  if(isMC) bckgProcessed--;
		}
	      GetListOfObject(Root,RootDir,histlist,(TDirectory*)file,"" );
	      file->Close();
	    }
	}

      return;
   }

   if(dir==NULL)return;

   TList* list = dir->GetListOfKeys();
   for(int i=0;i<list->GetSize();i++){
      TObject* tmp = GetObjectFromPath(dir,list->At(i)->GetName(),false);
      if(tmp->InheritsFrom("TTree")) continue;

      if(tmp->InheritsFrom("TDirectory")){
         GetListOfObject(Root,RootDir,histlist,(TDirectory*)tmp,parentPath+ list->At(i)->GetName()+"/" );
      }else if(tmp->InheritsFrom("TH1")){
        bool isTH1 = !(tmp->InheritsFrom("TH2") || tmp->InheritsFrom("TH3"));
        bool hasIndex = string(((TH1*)tmp)->GetXaxis()->GetTitle()).find("cut index")<string::npos;
        if(hasIndex){isTH1=true;}
	histlist.push_back(NameAndType(parentPath+list->At(i)->GetName(), isTH1, hasIndex ) );
      }

      delete tmp;
   }


}

void checkFileExist(JSONWrapper::Object& Root, std::string RootDir, NameAndType HistoProperties){

   std::vector<JSONWrapper::Object> Process = Root["proc"].daughters();
   for(unsigned int i=0;i<Process.size();i++){
      if(!StoreInFile && Process[i]["isinvisible"].toBool())continue;
      string filtExt("");
      if(Process[i].isTag("mctruthmode") ) { char buf[255]; sprintf(buf,"_filt%d",(int)Process[i]["mctruthmode"].toInt()); filtExt += buf; }
      std::vector<JSONWrapper::Object> Samples = (Process[i])["data"].daughters();
      for(unsigned int j=0;j<Samples.size();j++){
	 int split = 1;
	 if(Samples[j].isTag("split"))split = Samples[j]["split"].toInt();

	 TH1* tmphist = NULL;
         for(int s=0;s<split;s++){
	   string segmentExt; if(split>1) { char buf[255]; sprintf(buf,"_%i",s); segmentExt += buf; }
            string FileName = RootDir + (Samples[j])["dtag"].toString() + ((Samples[j].isTag("suffix"))?(Samples[j])["suffix"].toString():string("")) + segmentExt + filtExt + ".root";
            TFile* File = new TFile(FileName.c_str());
            bool& fileExist = FileExist[FileName];
            if(!File || File->IsZombie() || !File->IsOpen() || File->TestBit(TFile::kRecovered) ){fileExist=false;  continue; }else{fileExist=true;}

            TH1* tmptmphist = (TH1*) GetObjectFromPath(File,HistoProperties.name);
            if(tmptmphist)
              {
                if(!tmphist){gROOT->cd(); tmphist = (TH1*)tmptmphist->Clone(tmptmphist->GetName());}else{tmphist->Add(tmptmphist);}
                delete tmptmphist;
              }
            delete File;
         }
	 if(!tmphist)continue;
	 stSampleInfo& sampleInfo = sampleInfoMap[(Samples[j])["dtag"].toString()];
	 sampleInfo.initialNumberOfEvents = tmphist->GetBinContent(1);

	 delete tmphist;

      }
   }
}



void SavingToFile(JSONWrapper::Object& Root, std::string RootDir, NameAndType HistoProperties, TFile* OutputFile){
   std::vector<TObject*> ObjectToDelete;

   std::vector<JSONWrapper::Object> Process = Root["proc"].daughters();
   for(unsigned int i=0;i<Process.size();i++){
      TH1* hist = NULL;
      std::vector<JSONWrapper::Object> Samples = (Process[i])["data"].daughters();
      for(unsigned int j=0;j<Samples.size();j++){
         double Weight = 1.0;
/*
         if(!Process[i]["isdata"].toBool() && !Process[i]["isdatadriven"].toBool())Weight*= iLumi;
         if(Samples[j].isTag("xsec")     ) Weight*= Samples[j]["xsec"].toDouble();
	 if(Samples[j].isTag("BR"))             Weight*= Samples[j]["BR"].toDouble();
	 if(Samples[j].isTag("filtEff"))        Weight*= Samples[j]["filtEff"].toDouble();
	 if(Samples[j].isTag("kfactor"))        Weight*= Samples[j]["kfactor"].toDouble();
	 //if(Samples[j].isTag("nevts")     ) Weight/= Samples[j]["nevts"].toDouble();
	 stSampleInfo& sampleInfo = sampleInfoMap[(Samples[j])["dtag"].toString()];
         Weight /= sampleInfo.initialNumberOfEvents;
*/
	 string filtExt("");
	 if(Process[i].isTag("mctruthmode") ) { char buf[255]; sprintf(buf,"_filt%d",(int)Process[i]["mctruthmode"].toInt()); filtExt += buf; }

	 //std::vector<JSONWrapper::Object> BR = Samples[j]["br"].daughters();for(unsigned int b=0;b<BR.size();b++){Weight*=BR[b].toDouble();}

         if(HistoProperties.name.find("optim_cut")!=string::npos){Weight=1.0;}

         int split = 1;
         if(Samples[j].isTag("split"))split = Samples[j]["split"].toInt();
         TH1* tmphist = NULL;
         for(int s=0;s<split;s++){
	   string segmentExt; if(split>1){ char buf[255]; sprintf(buf,"_%i",s); segmentExt += buf;}

            string FileName = RootDir + (Samples[j])["dtag"].toString() + ((Samples[j].isTag("suffix"))?(Samples[j])["suffix"].toString():string("")) +  segmentExt + filtExt + ".root";
            if(!FileExist[FileName])continue;
            TFile* File = new TFile(FileName.c_str());
            if(!File || File->IsZombie() || !File->IsOpen() || File->TestBit(TFile::kRecovered) )continue;
            TH1* tmptmphist = (TH1*) GetObjectFromPath(File,HistoProperties.name);
            if(!tmptmphist){delete File;continue;}
            if(!tmphist){gROOT->cd(); tmphist = (TH1*)tmptmphist->Clone(tmptmphist->GetName());}else{tmphist->Add(tmptmphist);}

            delete tmptmphist;
            delete File;
         }
         if(!tmphist)continue;
         if(!hist){gROOT->cd(); hist = (TH1*)tmphist->Clone(tmphist->GetName());hist->Scale(Weight);}else{hist->Add(tmphist,Weight);}
         delete tmphist;
      }
      if(!hist)continue;

      string dirName = Process[i]["tag"].c_str();while(dirName.find("/")!=std::string::npos)dirName.replace(dirName.find("/"),1,"-");
      OutputFile->cd();
      TDirectory* subdir = OutputFile->GetDirectory(dirName.c_str());
      if(!subdir || subdir==OutputFile) subdir = OutputFile->mkdir(dirName.c_str());
      subdir->cd();
      hist->Write();
      delete hist;
   }
}



void fixExtremities(TH1* h,bool addOverflow, bool addUnderflow)
{
  if(h==0) return;

  if(addUnderflow){
      double fbin  = h->GetBinContent(0) + h->GetBinContent(1);
      double fbine = sqrt(h->GetBinError(0)*h->GetBinError(0)
                          + h->GetBinError(1)*h->GetBinError(1));
      h->SetBinContent(1,fbin);
      h->SetBinError(1,fbine);
      h->SetBinContent(0,0);
      h->SetBinError(0,0);
    }

  if(addOverflow){
      int nbins = h->GetNbinsX();
      double fbin  = h->GetBinContent(nbins) + h->GetBinContent(nbins+1);
      double fbine = sqrt(h->GetBinError(nbins)*h->GetBinError(nbins)
                          + h->GetBinError(nbins+1)*h->GetBinError(nbins+1));
      h->SetBinContent(nbins,fbin);
      h->SetBinError(nbins,fbine);
      h->SetBinContent(nbins+1,0);
      h->SetBinError(nbins+1,0);
    }
}


void Draw2DHistogramSplitCanvas(JSONWrapper::Object& Root, std::string RootDir, NameAndType HistoProperties){
   if(HistoProperties.isIndexPlot && cutIndex<0)return;

   std::string SaveName = "";

   TPaveText* T = new TPaveText(0.2,0.995,0.84,0.95, "NDC");
   T->SetFillColor(0);
   T->SetFillStyle(0);  T->SetLineColor(0);
   T->SetTextAlign(32);
   char Buffer[1024];
   if(isSim)	sprintf(Buffer, "#it{#bf{ATLAS}} Simulation, #sqrt{s}=%.1f TeV, L=%.1f fb^{-1}", iEcm, iLumi/1000);
   else sprintf(Buffer, "#it{#bf{ATLAS}} Internal, #sqrt{s}=%.1f TeV, L=%.1f fb^{-1}", iEcm, iLumi/1000);
   T->AddText(Buffer);

   std::vector<JSONWrapper::Object> Process = Root["proc"].daughters();
   std::vector<TObject*> ObjectToDelete;
   for(unsigned int i=0;i<Process.size();i++){
      if(Process[i]["isinvisible"].toBool())continue;
      string filtExt("");
      if(Process[i].isTag("mctruthmode") ) { char buf[255]; sprintf(buf,"_filt%d",(int)Process[i]["mctruthmode"].toInt()); filtExt += buf; }

      TCanvas* c1 = new TCanvas("c1","c1",500,500);
      c1->SetLogz(true);

      TH1* hist = NULL;
      std::vector<JSONWrapper::Object> Samples = (Process[i])["data"].daughters();
      for(unsigned int j=0;j<Samples.size();j++){
         double Weight = 1.0;
/*
         if(!Process[i]["isdata"].toBool()  && !Process[i]["isdatadriven"].toBool())Weight*= iLumi;
         if(Samples[j].isTag("xsec"))  		Weight*= Samples[j]["xsec"].toDouble();
	 if(Samples[j].isTag("BR"))             Weight*= Samples[j]["BR"].toDouble();
	 if(Samples[j].isTag("filtEff"))  	Weight*= Samples[j]["filtEff"].toDouble();
	 if(Samples[j].isTag("kfactor"))        Weight*= Samples[j]["kfactor"].toDouble();
	 //if(Samples[j].isTag("nevts")) 		Weight/= Samples[j]["nevts"].toDouble();
         stSampleInfo& sampleInfo = sampleInfoMap[(Samples[j])["dtag"].toString()];
         Weight /= sampleInfo.initialNumberOfEvents;
*/
         //std::vector<JSONWrapper::Object> BR = Samples[j]["br"].daughters();for(unsigned int b=0;b<BR.size();b++){Weight*=BR[b].toDouble();}


         int split = 1;
         if(Samples[j].isTag("split"))split = Samples[j]["split"].toInt();
         TH1* tmphist = NULL;
         for(int s=0;s<split;s++){
	   string segmentExt; if(split>1) { char buf[255]; sprintf(buf,"_%i",s); segmentExt += buf; }

            string FileName = RootDir + (Samples[j])["dtag"].toString() + ((Samples[j].isTag("suffix"))?(Samples[j])["suffix"].toString():string("")) + segmentExt + filtExt + ".root";
            if(!FileExist[FileName])continue;
            TFile* File = new TFile(FileName.c_str());
            if(!File || File->IsZombie() || !File->IsOpen() || File->TestBit(TFile::kRecovered) )continue;
            TH1* tmptmphist = (TH1*) GetObjectFromPath(File,HistoProperties.name);
            if(!tmptmphist){delete File;continue;}
            if(!tmphist){gROOT->cd(); tmphist = (TH1*)tmptmphist->Clone(tmptmphist->GetName());}else{tmphist->Add(tmptmphist);}
            delete tmptmphist;
            delete File;
         }
         if(!tmphist)continue;
         if(!hist){gROOT->cd(); hist = (TH1*)tmphist->Clone(tmphist->GetName());hist->Scale(Weight);}else{hist->Add(tmphist,Weight);}
         delete tmphist;
      }
      if(!hist)continue;

      SaveName = hist->GetName();
      ObjectToDelete.push_back(hist);
      hist->SetTitle("");
      hist->SetStats(kFALSE);

      gStyle->SetPalette(1);
      hist->Draw("COLZ");

      TPaveText* leg = new TPaveText(0.75,0.35,0.90,0.20, "NDC");
      leg->SetFillColor(0);
      leg->SetFillStyle(0);  leg->SetLineColor(0);
      leg->SetTextAlign(12);
      leg->AddText(Process[i]["tag"].c_str());
      ObjectToDelete.push_back(leg);
      T->Draw("same");


          //RJ adding ee or mumu channel number!
          TPaveText* T2 = new TPaveText(0.2+0.55,0.93,0.36+0.55,0.83, "NDC");
          T2->SetFillColor(0);
          T2->SetFillStyle(0);  T2->SetLineColor(0);
          T2->SetTextAlign(22);
          char Buffer2[1024];
          TString ssName = SaveName;
          if(SaveName.find("eeeq0jets") != string::npos)        sprintf(Buffer2,"#it{ee channel, 0 jets}");
          if(SaveName.find("eeeq1jets") != string::npos)        sprintf(Buffer2,"#it{ee channel, 1 jets}");
          if(SaveName.find("mumueq0jets") != string::npos)      sprintf(Buffer2,"#it{#mu#mu channel, 0 jets}");
          if(SaveName.find("mumueq1jets") != string::npos)      sprintf(Buffer2,"#it{#mu#mu channel, 1 jets}");
          if(SaveName.find("emueq0jets") != string::npos)       sprintf(Buffer2,"#it{e#mu channel, 0 jets}");
	  if(SaveName.find("emueq1jets") != string::npos)       sprintf(Buffer2,"#it{e#mu channel, 1 jets}");
          T2->AddText(Buffer2);

          if(SaveName.find("eeeq0jets") != string::npos || SaveName.find("eeeq1jets") != string::npos ||
                SaveName.find("mumueq0jets") != string::npos || SaveName.find("mumueq1jets") != string::npos ||
                SaveName.find("emueq0jets") != string::npos || SaveName.find("emueq1jets") != string::npos )
          {
                T2->Draw("same");
          }



      savePath(c1,outDir,SaveName+"_"+(Process[i])["tag"].toString(),plotExt);
      savePath(c1,outDir,SaveName+"_"+(Process[i])["tag"].toString(),".pdf");
      delete c1;
   }

   for(unsigned int d=0;d<ObjectToDelete.size();d++){delete ObjectToDelete[d];}ObjectToDelete.clear();
   delete T;
}


void Draw2DHistogram(JSONWrapper::Object& Root, std::string RootDir, NameAndType HistoProperties){
   if(HistoProperties.isIndexPlot && cutIndex<0)return;

   std::string SaveName = "";
   std::vector<JSONWrapper::Object> Process = Root["proc"].daughters();
   int NSampleToDraw = 0;
   for(unsigned int i=0;i<Process.size();i++){
      if(Process[i]["isinvisible"].toBool())continue;
      std::string procName = Process[i]["tag"].c_str();
      if(procName.find("Multijets")!=std::string::npos || procName.find("Z#rightarrow #tau#tau")!=std::string::npos) continue;
      NSampleToDraw++;
   }
   int CanvasX = 4;
   int CanvasY = ceil(NSampleToDraw/CanvasX);
   TCanvas* c1 = new TCanvas("c1","c1",CanvasX*350,CanvasY*350);
   c1->Divide(CanvasX,CanvasY,0,0);


   std::vector<TObject*> ObjectToDelete;
   int ipad = 0;
   for(unsigned int i=0;i<Process.size();i++){
      if(Process[i]["isinvisible"].toBool())continue;
      std::string procName = Process[i]["tag"].c_str();
      if(procName.find("Multijets")!=std::string::npos || procName.find("Z#rightarrow #tau#tau")!=std::string::npos) {
		continue;
      }

      TVirtualPad* pad = c1->cd(ipad+1);
      pad->SetLogz(true);
      pad->SetTopMargin(0.13); pad->SetBottomMargin(0.15);  pad->SetRightMargin(0.15); pad->SetLeftMargin(0.15);
      if(i%CanvasX == 0) pad->SetLeftMargin(0.2);

      TH1* hist = NULL;
      std::vector<JSONWrapper::Object> Samples = (Process[i])["data"].daughters();
      for(unsigned int j=0;j<Samples.size();j++){
         double Weight = 1.0;
//         if(!Process[i]["isdata"].toBool()  && !Process[i]["isdatadriven"].toBool())Weight*= iLumi;
	  string filtExt("");
	  if(Process[i].isTag("mctruthmode") ) { char buf[255]; sprintf(buf,"_filt%d",(int)Process[i]["mctruthmode"].toInt()); filtExt += buf; }

/*
         if(Samples[j].isTag("xsec")     )Weight*= Samples[j]["xsec"].toDouble();
	 if(Samples[j].isTag("BR")) 	        Weight*= Samples[j]["BR"].toDouble();
	 if(Samples[j].isTag("filtEff"))        Weight*= Samples[j]["filtEff"].toDouble();
	 if(Samples[j].isTag("kfactor"))        Weight*= Samples[j]["kfactor"].toDouble();
	 //if(Samples[j].isTag("nevts")     ) Weight/= Samples[j]["nevts"].toDouble();
         stSampleInfo& sampleInfo = sampleInfoMap[(Samples[j])["dtag"].toString()];
         Weight /= sampleInfo.initialNumberOfEvents;
*/

         //std::vector<JSONWrapper::Object> BR = Samples[j]["br"].daughters();for(unsigned int b=0;b<BR.size();b++){Weight*=BR[b].toDouble();}


         int split = 1;
         if(Samples[j].isTag("split"))split = Samples[j]["split"].toInt();
         TH1* tmphist = NULL;
         for(int s=0;s<split;s++){
	   string segmentExt; if(split>1) { char buf[255]; sprintf(buf,"_%i",s); segmentExt += buf; }

            string FileName = RootDir + (Samples[j])["dtag"].toString() + ((Samples[j].isTag("suffix"))?(Samples[j])["suffix"].toString():string("")) + segmentExt + filtExt + ".root";
            if(!FileExist[FileName])continue;
            TFile* File = new TFile(FileName.c_str());
            if(!File || File->IsZombie() || !File->IsOpen() || File->TestBit(TFile::kRecovered) )continue;
            TH1* tmptmphist = (TH1*) GetObjectFromPath(File,HistoProperties.name);
            if(!tmptmphist){delete File;continue;}
            if(!tmphist){gROOT->cd(); tmphist = (TH1*)tmptmphist->Clone(tmptmphist->GetName());}else{tmphist->Add(tmptmphist);}
            //if(Process[i]["isdata"].toBool())printf("%s --> %f\n",(Samples[j])["dtag"].toString().c_str(), tmptmphist->Integral());
            delete tmptmphist;
            delete File;
         }
         if(!tmphist)continue;
         if(!hist){gROOT->cd(); hist = (TH1*)tmphist->Clone(tmphist->GetName());hist->Scale(Weight);}else{hist->Add(tmphist,Weight);}
         delete tmphist;
      }
      if(!hist)continue;

      SaveName = hist->GetName();
      //std::string procName = Process[i]["tag"].c_str();
      //if(procName.find("Multijets")!=std::string::npos || procName.find("Z#rightarrow #tau#tau")!=std::string::npos) continue;
      ObjectToDelete.push_back(hist);
      hist->GetZaxis()->SetTitle("");
      hist->SetStats(kFALSE);
      hist->GetXaxis()->SetTitleOffset(1.3);
      hist->GetYaxis()->SetTitleOffset(2.0);
      hist->Draw("COL Z");

      //TString savename =getChannelName(SaveName);
      TPaveText* leg = new TPaveText(0.20,0.993,0.95,0.88, "NDC");
      leg->SetFillColor(0);
      leg->SetFillStyle(0);  leg->SetLineColor(0);
      leg->SetTextAlign(31);
      char Buffer[1024]; sprintf(Buffer, "#it{#bf{ATLAS}} Internal, #sqrt{s}=%.1f TeV, L=%.1f fb^{-1}", iEcm,iLumi/1000);
      leg->AddText(Buffer);
      TString processName = Process[i]["tag"].c_str();
      //processName += ", "+savename;
      //leg->AddText(Process[i]["tag"].c_str());
      leg->AddText(processName);
      leg->Draw("same");
      ObjectToDelete.push_back(leg);
      ipad++;
   }
   c1->cd(0);

   savePath(c1,outDir,SaveName,plotExt);
   savePath(c1,outDir,SaveName,".pdf");

   for(unsigned int d=0;d<ObjectToDelete.size();d++){delete ObjectToDelete[d];}ObjectToDelete.clear();
   delete c1;
}





void savePath(TCanvas* c,std::string outDir,std::string SaveName,std::string plotExt)
{
   string SavePath = SaveName + plotExt;
   while(SavePath.find("*")!=std::string::npos)SavePath.replace(SavePath.find("*"),1,"");
   while(SavePath.find("#")!=std::string::npos)SavePath.replace(SavePath.find("#"),1,"");
   while(SavePath.find("{")!=std::string::npos)SavePath.replace(SavePath.find("{"),1,"");
   while(SavePath.find("}")!=std::string::npos)SavePath.replace(SavePath.find("}"),1,"");
   while(SavePath.find("(")!=std::string::npos)SavePath.replace(SavePath.find("("),1,"");
   while(SavePath.find(")")!=std::string::npos)SavePath.replace(SavePath.find(")"),1,"");
   while(SavePath.find("^")!=std::string::npos)SavePath.replace(SavePath.find("^"),1,"");
   while(SavePath.find("/")!=std::string::npos)SavePath.replace(SavePath.find("/"),1,"-");
   if(outDir.size()) SavePath = outDir +"/"+ SavePath;
   system(string(("rm -f ") + SavePath).c_str());
   c->SaveAs(SavePath.c_str());
}



//
void Draw1DHistogram(JSONWrapper::Object& Root, std::string RootDir, NameAndType HistoProperties){
   if(HistoProperties.isIndexPlot && cutIndex<0)return;

   //gStyle->SetLineWidth(2);
   TCanvas* c1 = new TCanvas("c1","c1",700,700);
   TPad* t1 = new TPad();
   if(isDataBlind || noratio) t1 = new TPad("t1","t1", 0.0, 0.0, 1.0, 1.0);
   else t1 = new TPad("t1","t1", 0.0, 0.3, 1.0, 1.0);


   t1->Draw();
   t1->cd();
   if(!isDataBlind || noratio) t1->SetBottomMargin(0);
   TString name_denRelUncH="";
   if(!noLog) t1->SetLogy(true);
   float maximumFound(noLog);


   TLegend* legA  = new TLegend();
   if(!(isDataBlind||noratio)) legA  = new TLegend(0.50,0.80,0.95,0.98, "NDC");
   else legA = new TLegend(0.45,0.78-0.03,0.9,0.98-0.03, "NDC");

   legA->SetNColumns(3);
   THStack* stack = new THStack("MC","MC");
   TH1 *     mc   = NULL;
   TH1 *     mcPlusRelUnc = NULL;
   TH1 *     mctotalUnc = NULL;
   std::vector<TH1 *> spimpose;
   std::vector<TString> spimposeOpts;
   TH1*     data = NULL;
   std::vector<TObject*> ObjectToDelete;
   std::string SaveName = "";
   std::vector<JSONWrapper::Object> Process = Root["proc"].daughters();
   std::vector<TH1*> fakeStack;

   // // ===================  Alvaro has added  ==================
   // std::map<std::string,TGraphAsymmErrors*> vector_SystErrors;
   // TGraphAsymmErrors* sys_error = new TGraphAsymmErrors();
   // // =========================================================

   for(unsigned int i=0;i<Process.size();i++){
      if(Process[i]["isinvisible"].toBool())continue;
      TH1* hist = NULL;
      
      //      string FileSys;
      std::vector<JSONWrapper::Object> Samples = (Process[i])["data"].daughters();
      for(unsigned int j=0;j<Samples.size();j++){
         double Weight = 1.0;
//         if(!Process[i]["isdata"].toBool() && !Process[i]["isdatadriven"].toBool() )Weight*= iLumi;
	  string filtExt("");
	 if(Process[i].isTag("mctruthmode") ) { char buf[255]; sprintf(buf,"_filt%d",(int)Process[i]["mctruthmode"].toInt()); filtExt += buf; }

/*
         if(Samples[j].isTag("xsec")     )    Weight*= Samples[j]["xsec"].toDouble();
	 if(Samples[j].isTag("BR"))             Weight*= Samples[j]["BR"].toDouble();
	 if(Samples[j].isTag("filtEff"))        Weight*= Samples[j]["filtEff"].toDouble();
	 if(Samples[j].isTag("kfactor"))        Weight*= Samples[j]["kfactor"].toDouble();
	 //if(Samples[j].isTag("nevts")     )   Weight/= Samples[j]["nevts"].toDouble();
         stSampleInfo& sampleInfo = sampleInfoMap[(Samples[j])["dtag"].toString()];
         Weight /= sampleInfo.initialNumberOfEvents;
*/
         //std::vector<JSONWrapper::Object> BR = Samples[j]["br"].daughters();for(unsigned int b=0;b<BR.size();b++){Weight*=BR[b].toDouble();}

         if(HistoProperties.name.find("optim_cut")!=string::npos){Weight=1.0;}

         int split = 1;
         if(Samples[j].isTag("split"))split = Samples[j]["split"].toInt();
         TH1* tmphist = NULL;
         for(int s=0;s<split;s++){
	   string segmentExt; if(split>1) { char buf[255]; sprintf(buf,"_%i",s); segmentExt += buf; }

	    string FileName = RootDir + (Samples[j])["dtag"].toString() + ((Samples[j].isTag("suffix"))?(Samples[j])["suffix"].toString():string("")) + segmentExt + filtExt + ".root";
	    if(!FileExist[FileName]){continue;}
            TFile* File = new TFile(FileName.c_str());
            if(!File || File->IsZombie() || !File->IsOpen() || File->TestBit(TFile::kRecovered) )continue;
            TH1* tmptmphist = NULL;
            if(HistoProperties.isIndexPlot && cutIndex>=0){
               TH2* tmp2D = (TH2*) GetObjectFromPath(File,HistoProperties.name);
               if(tmp2D){tmptmphist = tmp2D->ProjectionY((string(tmp2D->GetName())+cutIndexStr).c_str(),cutIndex,cutIndex); delete tmp2D;}
//            }else if(HistoProperties.name.find("mt_shape")){
//               TH2* tmp2D = (TH2*) GetObjectFromPath(File,HistoProperties.name);
//               if(tmp2D){tmp2D->GetXaxis()->SetTitle("#events"); tmptmphist = tmp2D->ProjectionX((string(tmp2D->GetName())+cutIndexStr).c_str()); delete tmp2D;}
            }else if(!HistoProperties.isIndexPlot){
               tmptmphist = (TH1*) GetObjectFromPath(File,HistoProperties.name);
	    }
	    if(!tmptmphist){delete File;continue;}
	    if(tmptmphist->GetEntries() == 0) {delete File;continue;}
            if(!tmphist){gROOT->cd(); tmphist = (TH1*)tmptmphist->Clone(tmptmphist->GetName());}else{tmphist->Add(tmptmphist);}
            //if(Process[i]["isdata"].toBool())printf("%s --> %f\n",(Samples[j])["dtag"].toString().c_str(), tmptmphist->Integral());
	    //cout << ">>>>>>>>>>>>> debug tmphist: " << tmptmphist->GetNbinsX() << " weight: " << Weight << " File: " << FileName.c_str() << endl;

	    // ======================= Alvaro Added =============================
	    // if( !TString(FileName).Contains("data") )
	    //   FileSys = "/afs/cern.ch/work/a/alopezso/private/DMAnalysis20.7/Plots/"+(FileName.substr(FileName.find("mc15c.")+6)).replace(FileName.find(".NTuple.root"),-1,"")+"/"+(FileName.substr(FileName.find("mc15c.")+6)).replace(FileName.find(".NTuple.root"),-1,"_SysBand.root");
	    // else
	    //   FileSys="";
	    // ==================================================================
            delete tmptmphist;
            delete File;
         }
         if(!tmphist)continue;
         if(!hist){gROOT->cd(); hist = (TH1*)tmphist->Clone(tmphist->GetName());hist->Scale(Weight);}else{hist->Add(tmphist,Weight);}
         delete tmphist;
      }
      if(!hist)continue;
      // // ===================  Alvaro has added  ==================
      // if( FileSys != ""  || !TFile::Open( FileSys.c_str(),"READ")->IsOpen() ){
      // 	TFile* fileSyst = TFile::Open( FileSys.c_str(),"READ"); 
      // 	vector_SystErrors[Process[i]["tag"].c_str()] = (TGraphAsymmErrors*)fileSyst->Get("SysBand");

      // 	for( int j = 0 ; j < hist->GetNbinsX() ; ++j){
      // 	  double err_up = vector_SystErrors[Process[i]["tag"].c_str()]->GetErrorYhigh(j);
      // 	  double err_down = vector_SystErrors[Process[i]["tag"].c_str()]->GetErrorYlow(j);
      // 	  double err2_up =  sys_error ? sys_error->GetErrorYhigh(j) : 0;
      // 	  double err2_down = sys_error ? sys_error->GetErrorYlow(j) : 0;

      // 	  err2_up += pow(err_up*hist->GetBinContent(j+1),2);
      // 	  err2_down += pow(err_down*hist->GetBinContent(j+1),2);

      // 	  vector_SystErrors[Process[i]["tag"].c_str()]->SetPointError(j,0,0,err_down*hist->GetBinContent(j+1),err_up*hist->GetBinContent(j+1));
      // 	  sys_error->SetPoint(j,hist->GetXaxis()->GetBinCenter(j+1),1);	
      // 	  sys_error->SetPointError(j,0,0,err2_down,err2_up);	
      // 	}
      // }
      // // =========================================================

      SaveName = hist->GetName();
      TString postfix(""); postfix+=i;
      hist->SetName(SaveName+postfix);
      hist->Rebin(rebin);
      if(Process[i].isTag("color" ) )hist->SetLineColor  ((int)Process[i][ "color"].toDouble()); else hist->SetLineColor  (1);
      if(Process[i].isTag("color" ) )hist->SetMarkerColor((int)Process[i][ "color"].toDouble()); else hist->SetMarkerColor(1);
      if(Process[i].isTag("color" ) )hist->SetFillColor  ((int)Process[i][ "color"].toDouble()); else hist->SetFillColor  (0);
      if(Process[i].isTag("lcolor") )hist->SetLineColor  ((int)Process[i]["lcolor"].toDouble());
      if(Process[i].isTag("mcolor") )hist->SetMarkerColor((int)Process[i]["mcolor"].toDouble());
      if(Process[i].isTag("fcolor") )hist->SetFillColor  ((int)Process[i]["fcolor"].toDouble());
      if(Process[i].isTag("lwidth") )hist->SetLineWidth  ((int)Process[i]["lwidth"].toDouble());// else hist->SetLineWidth  (1);
      if(Process[i].isTag("lstyle") )hist->SetLineStyle  ((int)Process[i]["lstyle"].toDouble());// else hist->SetLinStyle  (1);
      if(Process[i].isTag("fill"  ) )hist->SetFillColor  ((int)Process[i]["fill"  ].toDouble());
      if(Process[i].isTag("marker") )hist->SetMarkerStyle((int)Process[i]["marker"].toDouble());// else hist->SetMarkerStyle(1);

      //fixExtremities(hist,true,true);
      hist->SetTitle("");
      hist->SetStats(kFALSE);
      hist->SetMinimum(2e-2*scaleYMin);//5e-1);//2e-2);
      //hist->SetMaximum(1E6);
      //hist->SetMaximum(hist->GetBinContent(hist->GetMaximumBin())*1.10);
      //TString tSaveName = hist->GetName();

      //cout <<hist->GetNbinsX() << endl;

      ObjectToDelete.push_back(hist);
      if(Process[i].isTag("normto")) hist->Scale( Process[i]["normto"].toDouble()/hist->Integral() );
      if((!Process[i].isTag("spimpose") || !Process[i]["spimpose"].toBool()) && !Process[i]["isdata"].toBool()){
         //Add to Stack
	 //cout << "------------" << endl;
	 //cout << "stack: " << Process[i]["tag"].c_str() << "\t GetNbinsX: " << hist->GetNbinsX()<< endl;
	 if(Process[i]["tag"].c_str() == ctrlSample)
		stack->Add(hist, "HIST");
	 else
		fakeStack.push_back(hist);

         legA->AddEntry(hist, Process[i]["tag"].c_str(), "F");
         if(!mc){mc = (TH1D*)hist->Clone("mc");}else{mc->Add(hist);}
      }
      else if(Process[i].isTag("spimpose") && Process[i]["spimpose"].toBool())
	{
   	  //legB->AddEntry(hist, Process[i]["tag"].c_str(), "L");
	  if(!nosig) {
		if(scaleSignal==1) legA->AddEntry(hist, Process[i]["tag"].c_str(), Process[i]["isdata"].toBool() ? "P" : "L" );
		else{
                	std::ostringstream strs;
                	strs << scaleSignal;
                	std::string _str_ = strs.str();
                	TString str_scalesignal="#times"+_str_;
                	TString newleg = Process[i]["tag"].c_str()+str_scalesignal;
                	legA->AddEntry(hist, newleg, Process[i]["isdata"].toBool() ? "P" : "L" );
		}
	  }

	  spimposeOpts.push_back( Process[i]["isdata"].toBool() ? "e1" : "hist" );
	  hist->Scale(scaleSignal); //scale the signal, will not affect the root file
	  spimpose.push_back(hist);
	  if(maximumFound<hist->GetMaximum()) maximumFound=hist->GetMaximum()*5.0;//1.1;
	}
      else{
	if(Process[i]["isdata"].toBool()){
	  if(!data){
	    data = hist;
	    legA->AddEntry(hist, Process[i]["tag"].c_str(), "P");
	  }
	  else data->Add(hist);
	}
      }
   }
//   if(mc &&  maximumFound<mc->GetMaximum()) maximumFound=mc->GetMaximum()*5.0;//1.1;
//   if(data && maximumFound<data->GetMaximum()) maximumFound=data->GetMaximum()*5.0;//1.1;

   //cout << mc->GetNbinsX() << endl;
   if(mc &&  maximumFound<mc->GetBinContent(mc->GetMaximumBin()) ) maximumFound=mc->GetBinContent(mc->GetMaximumBin())*5.0;//1.1;
   if(data &&  maximumFound<data->GetBinContent(data->GetMaximumBin()) ) maximumFound=data->GetBinContent(data->GetMaximumBin())*5.0;//1.1;

   bool canvasIsFilled(false);

   for(unsigned int j=0; j<fakeStack.size(); j++){
	stack->Add(fakeStack[j], "HIST");
   }

   if(stack && stack->GetStack() && stack->GetStack()->GetEntriesFast()>0){
     stack->Draw("");

     stack->GetYaxis()->SetTitleOffset(0.9);
     stack->GetYaxis()->SetLabelSize(0.06);
     stack->GetYaxis()->SetTitleSize(0.06);

     if(noratio){
	stack->GetXaxis()->SetTitleOffset(0.9);
	stack->GetXaxis()->SetLabelSize(0.04);
	stack->GetXaxis()->SetTitleSize(0.05);

	stack->GetYaxis()->SetTitleOffset(1.1);
	stack->GetYaxis()->SetLabelSize(0.04);
	stack->GetYaxis()->SetTitleSize(0.05);
     }

     TH1 *hist=(TH1*)stack->GetStack()->At(0);
     if(stack->GetXaxis())
       {
	 if(isDataBlind || noratio) stack->GetXaxis()->SetTitle(hist->GetXaxis()->GetTitle());
	 name_denRelUncH = hist->GetXaxis()->GetTitle();

	 TString tSaveName = hist->GetName();
	 float binsize = hist->GetBinWidth(1);
	 float binsizeN = hist->GetBinWidth(hist->GetXaxis()->GetNbins());
	 if(binsize!=binsizeN) binsize=1;
	 std::ostringstream strs;
	 strs.precision(2);
	 strs << binsize;
	 TString binSize = strs.str();
	 //if(tSaveName.Contains("pfmet_dataDY") || tSaveName.Contains("qt_dataDY") || tSaveName.Contains("pfmetType2_dataDY")) binSize = "1";
	 if(name_denRelUncH.Contains("multiplicity")) stack->GetYaxis()->SetTitle("Events");
	 else if(!name_denRelUncH.Contains("1 GeV")){
	 	if(name_denRelUncH.Contains("GeV") &&
			!name_denRelUncH.Contains(">") && !name_denRelUncH.Contains("<"))
				stack->GetYaxis()->SetTitle("Events / "+binSize+" GeV");
	 	else 			stack->GetYaxis()->SetTitle("Events / "+binSize);
	 }


	 double minimumFound = hist->GetMinimum();
	 stack->SetMinimum(minimumFound);

	 if(tSaveName.Contains("njets_raw")) maximumFound*= 50;
	 if(tSaveName.Contains("dphiZMET_presel")) maximumFound*= 50;
	 if(tSaveName.Contains("balancedif_presel")) maximumFound*= 50;
         if(tSaveName.Contains("axialpfmet_presel")) maximumFound*= 100;

      	 if(tSaveName.Contains("nvtx")) maximumFound*= 50;
      	 if(tSaveName.Contains("nleptons")) maximumFound*= 50;
         if(tSaveName.Contains("npfjets")) maximumFound*= 50;
         if(tSaveName.Contains("npfbjets")) maximumFound*= 50;
	 if(tSaveName.Contains("_WWCtrl")) maximumFound*= 50;
	 if(tSaveName.Contains("zmassType2_WWCtrl")) maximumFound = 1e+5;
	 if(tSaveName.Contains("Eta_raw")) maximumFound*= 40;
	 if(tSaveName.Contains("wmt_raw")) maximumFound*= 20;
	 if(tSaveName.Contains("pfmet_raw")) maximumFound*= 100;
	 if(tSaveName.Contains("DPhiZMET_dataDY")) maximumFound*= 20;
	 if(tSaveName.Contains("eleLooseFakePt")) maximumFound*= 20;
	 if(tSaveName.Contains("eleTightFakePt")) maximumFound*= 20;
	 if(tSaveName.Contains("FakeEta")) maximumFound*= 40;
	 if(tSaveName.Contains("eventflow")) maximumFound*= 100;

	 if(tSaveName.Contains("pho_phi_sel")) maximumFound*= 50;
         if(tSaveName.Contains("pho_eta_sel")) maximumFound*= 50;
         if(tSaveName.Contains("pho_pt_sel")) maximumFound*= 50;
         if(tSaveName.Contains("met_sel")) maximumFound*= 50;
         if(tSaveName.Contains("diphoton_pt_sel") || tSaveName.Contains("pthard")) maximumFound*= 50;
         if(tSaveName.Contains("diphoton_mass_sel")) maximumFound*= 50;
         if(tSaveName.Contains("balancedif_sel")) maximumFound*= 50;
         if(tSaveName.Contains("dphi")) maximumFound*= 50;
	 if(tSaveName.Contains("nphotons") || tSaveName.Contains("nelectrons") || tSaveName.Contains("nmuons") || tSaveName.Contains("njets") ) maximumFound*= 100;
	 if(tSaveName.Contains("metphi_sel")) maximumFound*= 50;
	 if(tSaveName.Contains("pu_sel") || tSaveName.Contains("puwgt_sel")) maximumFound*= 50;
         if(tSaveName.Contains("yields")) maximumFound*= 100;
	 if(tSaveName.Contains("metsig")) maximumFound*= 100;


	 maximumFound *= scaleYMax;

	 stack->SetMaximum(maximumFound);

/*
         if(tSaveName.Contains("pfmet")){
                stack->SetMinimum(5e-5);
                //if(tSaveName.Contains("eq0jets")) stack->SetMaximum(1e3);
                //if(tSaveName.Contains("eq1jets")) stack->SetMaximum(2e4);
         }

*/
         if(tSaveName.Contains("mt_final")){
                stack->SetMaximum(100);
         }


	 if(noLog)
	   {
	     stack->SetMaximum(maximumFound);
	   }
       }
     ObjectToDelete.push_back(stack);
     canvasIsFilled=true;


     if(showUnc && mc)
       {
	 mcPlusRelUnc = (TH1 *) mc->Clone("totalmcwithunc"); mcPlusRelUnc->SetDirectory(0);
	 for(int ibin=1; ibin<=mcPlusRelUnc->GetXaxis()->GetNbins(); ibin++)
	   {
	     Double_t error=sqrt(pow(mcPlusRelUnc->GetBinError(ibin),2)+pow(mcPlusRelUnc->GetBinContent(ibin)*baseRelUnc,2));
	     mcPlusRelUnc->SetBinError(ibin,error);
	   }
	 mcPlusRelUnc->SetFillStyle(3427);
	 mcPlusRelUnc->SetFillColor(kGray+1);
	 mcPlusRelUnc->SetMarkerStyle(1);
       }
   }
   if(data){
       if(!isDataBlind) data->Draw(canvasIsFilled ? "E1 same" : "E1");
       canvasIsFilled=true;
   }
   for(size_t ip=0; ip<spimpose.size(); ip++){
     TString opt=spimposeOpts[ip];
     if(!nosig) spimpose[ip]->Draw(opt + (canvasIsFilled ? "same": "") );
     canvasIsFilled=true;
   }

   if(mc)
   {
	mctotalUnc=(TH1D *) mc->Clone("mcrelunc");
       	for(int xbin=1; xbin<=mctotalUnc->GetXaxis()->GetNbins(); xbin++)
          {
          	if(mctotalUnc->GetBinContent(xbin)==0) continue;
           	mctotalUnc->SetBinContent(xbin,mctotalUnc->GetBinContent(xbin));
           	mctotalUnc->SetBinError(xbin,mctotalUnc->GetBinError(xbin));
          }
        mctotalUnc->SetDirectory(0);

	TGraphErrors *mcgr=new TGraphErrors;
        //TGraphAsymmErrors *mcgr=new TGraphAsymmErrors;
        for(int ibin=1; ibin<=mctotalUnc->GetXaxis()->GetNbins(); ibin++)
          {
	    // // Alvaro includes 

	    // double err_up = sys_error->GetErrorYhigh(ibin-1);
	    // double err_down = sys_error->GetErrorYlow(ibin-1) ;	     
	    // double err_stat = mctotalUnc->GetBinError(ibin);
	    // //

        	mcgr->SetPoint(ibin-1,mctotalUnc->GetXaxis()->GetBinCenter(ibin),mctotalUnc->GetBinContent(ibin));
		mcgr->SetPointError(ibin-1,mctotalUnc->GetXaxis()->GetBinWidth(ibin)/2,mctotalUnc->GetBinError(ibin));
                //mcgr->SetPointError(ibin-1 , 0 , 0 , mctotalUnc->GetBinContent(ibin)-sqrt(pow(err_stat,2) + err_down ) , mctotalUnc->GetBinContent(ibin)+sqrt(pow(err_stat,2) + err_up ) );
          }
       	mcgr->SetFillStyle(3254);
       	mcgr->SetFillColor(1);
        mcgr->SetMarkerStyle(1);
       	mcgr->Draw("2 same");
	legA->AddEntry(mcgr,"Stat. Unc.", "F");
   }

   //compare data and MC
   if(showChi2 && data && mc && data->Integral()>0 && mc->Integral()>0)
     {
       TPaveText *pave = new TPaveText(0.6,0.85,0.8,0.9,"NDC");
       pave->SetBorderSize(0);
       pave->SetFillStyle(0);
       pave->SetTextAlign(32);
       pave->SetTextFont(42);
       char buf[100];
       sprintf(buf,"#chi^{2}/ndof : %3.2f", data->Chi2Test(mc,"WWCHI2/NDF") );
       pave->AddText(buf);
       pave->Draw();
     }

//   TPaveText* T = new TPaveText(0.7,0.994,0.95,0.94, "NDC");
   TPaveText* T = new TPaveText(0.18,0.80,0.43,0.85, "NDC");
   if(isDataBlind || noratio) T = new TPaveText(0.7,0.994,0.95,0.95, "NDC");
   T->SetFillColor(0);
   T->SetFillStyle(0);  T->SetLineColor(0);
   T->SetTextAlign(22);
   char Buffer[1024];
   //if(iLumi>1000) sprintf(Buffer, "%.1f fb^{-1} (%.0f TeV)", iLumi/1000,iEcm);
   //else		  sprintf(Buffer, "%.1f pb^{-1} (%.0f TeV)", iLumi, iEcm);
   if(iLumi>1000) sprintf(Buffer, "#sqrt{s} = %.0f TeV, %.1f fb^{-1}", iEcm, iLumi/1000);
   else 	  sprintf(Buffer, "#sqrt{s} = %.0f TeV, %.1f pb^{-1}", iEcm, iLumi);


   T->AddText(Buffer);
   T->Draw("same");
   T->SetBorderSize(0);

   //T = new TPaveText(0.06,0.994,0.58, 0.935, "NDC");
   T = new TPaveText(0.18, 0.83, 0.43, 0.96, "NDC");
   if(isDataBlind || noratio) T = new TPaveText(0.05,0.994,0.45, 0.95, "NDC");
   if(isSim) T->AddText("#it{#bf{ATLAS}} Simulation");
   else if(isinProgress) T->AddText("#it{#bf{ATLAS}} Work in Progress");
   else      T->AddText("#it{#bf{ATLAS}} Internal");

   T->Draw("same");
   T->SetBorderSize(0);

/*
   TPaveText* T2 = new TPaveText(0.2,0.93,0.39,0.8, "NDC");
   if(noratio) T2 = new TPaveText(0.2,0.93,0.39,0.88, "NDC");
   T2->SetFillColor(0);
   T2->SetFillStyle(0);  T2->SetLineColor(0);
   T2->SetTextAlign(22);
   TString savename = getChannelName(SaveName);
   T2->AddText(savename);
   if(savename!="") T2->Draw("same");
*/
   //RJ
   TGraph *leftcut;
   if(leftcutLine != 0){
   	const int n = 35;
	double x1[n],y1[n];
	for (int i=0;i<n;i++) {
		x1[i] = leftcutLine;
		double ratio = (double) i/n;
     		y1[i] = maximumFound*ratio;
   	}
   	leftcut = new TGraph(n,x1,y1);
   	leftcut->SetLineColor(6);
	leftcut->SetFillColor(6);
   	leftcut->SetLineWidth(603);
   	leftcut->SetFillStyle(3005);
	leftcut->Draw("same");
   }


   TGraph *rightcut;
   if(rightcutLine != 0){
   	const int n = 35;
	double x1[n],y1[n];
	for (int i=0;i<n;i++) {
		x1[i] = rightcutLine;
		double ratio = (double) i/n;
     		y1[i] = maximumFound*ratio;
   	}
   	rightcut = new TGraph(n,x1,y1);
   	rightcut->SetLineColor(6);
	rightcut->SetFillColor(6);
   	rightcut->SetLineWidth(-603);
   	rightcut->SetFillStyle(3005);
	rightcut->Draw("same");
   }

   //RJ
   TLine *llcut1;
   if (cutLine1 != 0){
	llcut1 = new TLine(cutLine1, 0, cutLine1, maximumFound);
	llcut1->SetLineWidth(2);
	llcut1->SetLineStyle(2);
	llcut1->SetLineColor(kBlack);
	llcut1->Draw("same");
  }

   TLine *llcut2;
   if (cutLine2 != 0){
        llcut2 = new TLine(cutLine2, 0, cutLine2, maximumFound);
        llcut2->SetLineWidth(2);
        llcut2->SetLineStyle(2);
        llcut2->SetLineColor(kBlack);
        llcut2->Draw("same");
  }


   legA->SetFillColor(0); legA->SetFillStyle(0); legA->SetLineColor(0);
   legA->SetHeader("");
   legA->Draw("same");
   legA->SetTextFont(42);
   //    legB->SetFillColor(0); legB->SetFillStyle(0); legB->SetLineColor(0);
   //    legB->SetHeader("");
   //    legB->Draw("same");
   //   legB->SetTextFont(42);

   //
   std::vector<TH1 *> compDists;
   if(data)                   compDists.push_back(data);
   else if(spimpose.size()>0) compDists=spimpose;
   if(mc && compDists.size() && !isDataBlind && !noratio)
     {
       t1->SetTopMargin(0.03);
       t1->SetRightMargin(0.03);
       c1->cd();
       TPad* t2 = new TPad("t2","t2", 0.0, 0.0, 1.0, 0.3);
       t2->Draw();
       t2->cd();
       t2->SetGridy(true);
       //t2->SetTopMargin(0);
       t2->SetTopMargin(0.05);
       t2->SetRightMargin(0.03);
       t2->SetBottomMargin(0.5);
       //mc stats
       TH1D *denRelUncH=0;
       if(mcPlusRelUnc) denRelUncH=(TH1D *) mcPlusRelUnc->Clone("mcrelunc");
       else             denRelUncH=(TH1D *) mc->Clone("mcrelunc");
       for(int xbin=1; xbin<=denRelUncH->GetXaxis()->GetNbins(); xbin++)
	 {
	   if(denRelUncH->GetBinContent(xbin)==0) continue;
	   Double_t err=denRelUncH->GetBinError(xbin)/denRelUncH->GetBinContent(xbin);
	   denRelUncH->SetBinContent(xbin,1);
	   denRelUncH->SetBinError(xbin,err);
	 }
       //TGraphErrors *denRelUnc=new TGraphErrors(denRelUncH);
       TGraphErrors *denRelUnc=new TGraphErrors(denRelUncH->GetXaxis()->GetNbins());
       int icutgRJ=0;
       double llmin=-999;
       double llmax=999;
       for(int ibin=1; ibin<=denRelUncH->GetXaxis()->GetNbins(); ibin++){
           	//double err=denRelUncH->GetBinError(ibin)/denRelUncH->GetBinContent(ibin);
		if(ibin == 1) llmin=denRelUncH->GetXaxis()->GetBinLowEdge(ibin);
		llmax=denRelUncH->GetXaxis()->GetBinUpEdge(ibin);

		denRelUnc->SetPoint(icutgRJ,denRelUncH->GetXaxis()->GetBinCenter(ibin),1);
		if(useDataMinusMC) denRelUnc->SetPoint(icutgRJ,denRelUncH->GetXaxis()->GetBinCenter(ibin),0);
		if(denRelUncH->GetBinContent(ibin) < 2e-2) continue;
		denRelUnc->SetPointError(icutgRJ,denRelUncH->GetXaxis()->GetBinWidth(ibin)/2.0, denRelUncH->GetBinError(ibin)/denRelUncH->GetBinContent(ibin));
          	icutgRJ++;
       } //RJ


       denRelUnc->SetLineColor(1);
       denRelUnc->SetFillStyle(3254);
       denRelUnc->SetFillColor(kRed);
       denRelUnc->SetMarkerColor(1);
       denRelUnc->SetMarkerStyle(1);
       denRelUncH->Reset("ICE");
       denRelUncH->Draw();
       denRelUnc->Draw("2");

       TLine *llbase;
       llbase = new TLine(llmin,1.,llmax,1.);
       if(useDataMinusMC) llbase = new TLine(llmin,0.,llmax,0.);
       llbase->SetLineWidth(1);
       llbase->SetLineStyle(1);
       llbase->SetLineColor(kBlue);
       llbase->Draw("same");


       float yscale = (1.0-0.2)/(0.18-0);
       denRelUncH->GetYaxis()->SetTitle("Data/MC");
       if(useDataMinusMC) denRelUncH->GetYaxis()->SetTitle("#frac{Data-MC}{MC}");
       denRelUncH->SetMinimum(0.09);
       denRelUncH->SetMaximum(1.99);
       if(useDataMinusMC) denRelUncH->SetMinimum(-1.99);
       if(useDataMinusMC) denRelUncH->SetMaximum(1.99);
       //denRelUncH->GetXaxis()->SetTitle("");
       denRelUncH->GetXaxis()->SetTitle(name_denRelUncH); //RJ
       //denRelUncH->SetMinimum(0);
       //denRelUncH->SetMaximum(data->GetBinContent(data->GetMaximumBin())*1.10);
       denRelUncH->GetXaxis()->SetTitleOffset(1.2/*1.3*/);
       denRelUncH->GetXaxis()->SetLabelSize(0.03*yscale);
       denRelUncH->GetXaxis()->SetTitleSize(/*0.03*/0.04*yscale);
       denRelUncH->GetXaxis()->SetTickLength(0.03*yscale);
       denRelUncH->GetYaxis()->SetTitleOffset(0.39/*0.33*/);
       denRelUncH->GetYaxis()->SetNdivisions(5);
       denRelUncH->GetYaxis()->SetLabelSize(0.03*yscale);
       denRelUncH->GetYaxis()->SetTitleSize(0.03*yscale);


       //add comparisons
       for(size_t icd=0; icd<compDists.size(); icd++)
	 {
	   TString name("CompHistogram"); name+=icd;
	   TH1D *dataToObsH = (TH1D*)compDists[icd]->Clone(name);
	   if(useDataMinusMC) dataToObsH->Add(mc,-1);
	   dataToObsH->Divide(mc);
	   dataToObsH->Draw("E1 same");
	 }
     }
   else
     {
       //if not comparison resize the canvas
       c1->SetWindowSize(800,700);
       c1->SetCanvasSize(800,700);
       t1->SetPad(0,0,1,1);
       t1->SetBottomMargin(0.12);
       t1->SetTopMargin(0.05);
       t1->SetLeftMargin(0.15);
       if(!stack) stack->GetXaxis()->SetTitle(name_denRelUncH);
     }

   t1->RedrawAxis();

   c1->Modified();
   c1->Update();
   c1->cd();
   savePath(c1,outDir,SaveName,plotExt);
   savePath(c1,outDir,SaveName,".pdf");

   delete c1;
   for(unsigned int d=0;d<ObjectToDelete.size();d++){delete ObjectToDelete[d];}ObjectToDelete.clear();
   delete legA;
   delete T;
}



TString getChannelName(std::string SaveName){
   TString Buffer2;
   if(SaveName.find("eeeq0jets") != string::npos)       {Buffer2="#it{ee, 0 jets channel}";}
   if(SaveName.find("eeeq1jets") != string::npos)       {Buffer2="#it{ee, 1 jets channel}";}
   if(SaveName.find("eelesq1jets") != string::npos)     {Buffer2="#it{ee, #leq 1 jets channel}";}
   if(SaveName.find("eegeq2jets") != string::npos)      {Buffer2="#it{ee, #geq 2 jets channel}";}
   if(SaveName.find("mumueq0jets") != string::npos)     {Buffer2="#it{#mu#mu, 0 jets channel}";}
   if(SaveName.find("mumueq1jets") != string::npos)     {Buffer2="#it{#mu#mu, 1 jets channel}";}
   if(SaveName.find("mumulesq1jets") != string::npos)   {Buffer2="#it{#mu#mu, #leq 1 jets channel}";}
   if(SaveName.find("mumugeq2jets") != string::npos)    {Buffer2="#it{#mu#mu, #geq 2 jets channel}";}
   if(SaveName.find("emueq0jets") != string::npos)      {Buffer2="#it{e#mu, 0 jets channel}";}
   if(SaveName.find("emueq1jets") != string::npos)      {Buffer2="#it{e#mu, 1 jets channel}";}
   if(SaveName.find("emulesq1jets") != string::npos)    {Buffer2="#it{e#mu, #leq 1 jets channel}";}
   if(SaveName.find("emugeq2jets") != string::npos)     {Buffer2="#it{e#mu, #geq 2 jets channel}";}
   if(SaveName.find("lleq0jets") != string::npos)       {Buffer2="#it{ee,#mu#mu, 0 jets channel}";}
   if(SaveName.find("lleq1jets") != string::npos)       {Buffer2="#it{ee,#mu#mu, 1 jets channel}";}
   if(SaveName.find("lllesq1jets") != string::npos)     {Buffer2="#it{ee,#mu#mu, #leq 1 jets channel}";}
   if(SaveName.find("llgeq2jets") != string::npos)      {Buffer2="#it{ee,#mu#mu, #geq 2 jets channel}";}
   if(SaveName.find("ee") != string::npos && SaveName.find("eq0jets") == string::npos
        && SaveName.find("eq1jets") == string::npos && SaveName.find("lesq1jets") == string::npos
        && SaveName.find("geq2jets") == string::npos )  {Buffer2="#it{ee, #geq 0 jets channel}";}
   if(SaveName.find("mumu") != string::npos && SaveName.find("eq0jets") == string::npos
        && SaveName.find("eq1jets") == string::npos && SaveName.find("lesq1jets") == string::npos
        && SaveName.find("geq2jets") == string::npos )  {Buffer2="#it{#mu#mu, #geq 0 jets channel}";}
   if(SaveName.find("emu") != string::npos && SaveName.find("eq0jets") == string::npos
        && SaveName.find("eq1jets") == string::npos && SaveName.find("lesq1jets") == string::npos
        && SaveName.find("geq2jets") == string::npos )  {Buffer2="#it{e#mu, #geq 0 jets channel}";}
   if(SaveName.find("ll") != string::npos && SaveName.find("eq0jets") == string::npos
        && SaveName.find("eq1jets") == string::npos && SaveName.find("lesq1jets") == string::npos
        && SaveName.find("geq2jets") == string::npos )  {Buffer2="#it{ee,#mu#mu, #geq 0 jets channel}";}

   return Buffer2;
}



int main(int argc, char* argv[]){

   setTDRStyle();
   gStyle->SetPadTopMargin   (0.06);
   gStyle->SetPadBottomMargin(0.12);
   gStyle->SetPadRightMargin (0.06);
   gStyle->SetPadLeftMargin  (0.14);
   gStyle->SetTitleSize(0.04, "XYZ");
   gStyle->SetTitleXOffset(1.1);
   gStyle->SetTitleYOffset(1.45);
   gStyle->SetPalette(1);
   gStyle->SetNdivisions(505);

   std::vector<string> histoNameMask;
   std::vector<string> histoNameMaskStart;
   std::vector<string> nohistoNameMask;
   std::vector<string> histoChannelName;

   for(int i=1;i<argc;i++){
     string arg(argv[i]);

     if(arg.find("--help")!=string::npos){
	printf("\033[35m");
	printf("%s\n",std::string(100,'-').c_str());
        printf("--%18s --> print this helping text\n","help");

	printf("--%18s --> integrated luminosity to be used for the MC rescale\n","iLumi");
	printf("--%18s --> center of mass energy (TeV) = 8 TeV by default\n","iEcm");
        printf("--%18s --> print ATLAS Simulation instead of the standard title\n","isSim");
        printf("--%18s --> path to the directory containing the .root files to process\n","inDir");
        printf("--%18s --> path of the directory that will contains the output plots and tables\n","outDir");
        printf("--%18s --> path of the output summary .root file\n","outFile");
        printf("--%18s --> containing list of process (and associated style) to process to process\n","json");
        printf("--%18s --> processing only the objects containing the following argument in their name\n","only");
        printf("--%18s --> processing only the objects starting with the following argument in their name\n","OnlyStartWith");
	printf("--%18s --> processing only the objects starting without the following argument in their name\n","NoPlotwith");
        printf("--%18s --> processing only the objects with channel name\n","channel");
        printf("--%18s --> will do the projection on that index for histos of type cutIndex\n","index");
        printf("--%18s --> show the data/MC chi^2\n","chi2");
        printf("--%18s --> show stat uncertainty (if number is given use it as relative bin by bin uncertainty (e.g. lumi)\n","showUnc");
	printf("--%18s --> use linear scale\n","noLog");
	printf("--%18s --> use Ratio: (Data-MC)/MC\n","useDataMinusMC");
        printf("--%18s --> Skip processing of 1D objects\n","no1D");
        printf("--%18s --> Skip processing of 2D objects\n","no2D");
        printf("--%18s --> Do not create latex table (when possible)\n","noTex");
        printf("--%18s --> Do not make a summary .root file\n","noRoot");
        printf("--%18s --> Do not creates plot files (useful to speedup processing)\n","noPlot");
	printf("--%18s --> extension to save\n","plotExt");
	printf("--%18s --> name of the histogram with the original number of events (cutflow by default)\n","cutflow");
        printf("--%18s --> (only for 2D plots) save all the samples in separated pltos\n","splitCanvas");
	printf("--%18s --> number of bin (default = 1) for rebinning histograms, combine --only option to rebin some specific histograms\n","rebin");
	printf("--%18s --> Blind the Data point from 1D Histograms\n","isDataBlind");
	printf("--%18s --> No ratio plots\n","noratio");
	printf("--%18s --> Add cut Line 1\n","addcutLine1");
	printf("--%18s --> Add cut Line 2\n","addcutLine2");
	printf("--%18s --> Add left cut Line \n","addleftcutLine");
	printf("--%18s --> Add right cut Line \n","addrightcutLine");
	printf("--%18s --> Scale the signal with a number \n","scaleSignal");
	printf("--%18s --> Scale Ymax for 1D Histogram \n","scaleYMax");
	printf("--%18s --> Scale Ymin for 1D Histogram \n","scaleYMix");
	printf("%s\n",std::string(100,'-').c_str());
	printf("\033[0m");
	printf("\033[33m");
        printf("command line example: runPlotter --json ../data/beauty-samples.json --iLumi 2007 --inDir OUT/ --outDir OUT/plots/ --outFile plotter.root --noRoot --noPlot\n");
	printf("\033[0m");
	return 0;
     }

     if(arg.find("--iLumi"  )!=string::npos && i+1<argc){ sscanf(argv[i+1],"%lf",&iLumi); i++; printf("Lumi = %f\n", iLumi); }
     if(arg.find("--iEcm"  )!=string::npos && i+1<argc){ sscanf(argv[i+1],"%lf",&iEcm); i++; printf("Ecm = %f TeV\n", iEcm); }

     if(arg.find("--inDir"  )!=string::npos && i+1<argc){ inDir    = argv[i+1];  i++;  printf("inDir = %s\n", inDir.c_str());  }
     if(arg.find("--outDir" )!=string::npos && i+1<argc){ outDir   = argv[i+1];  i++;  printf("outDir = %s\n", outDir.c_str());  }
     if(arg.find("--outFile")!=string::npos && i+1<argc){ outFile  = argv[i+1];  i++; printf("output file = %s\n", outFile.c_str()); }
     if(arg.find("--json"   )!=string::npos && i+1<argc){ jsonFile = argv[i+1];  i++;  }
     if(arg.find("--OnlyStartWith"   )!=string::npos && i+1<argc){ histoNameMaskStart.push_back(argv[i+1]); printf("Only process Histo starting with '%s'\n", argv[i+1]); i++;  }
     if(arg.find("--only"   )!=string::npos && i+1<argc)         { histoNameMask.push_back(argv[i+1]); printf("Only process Histo containing '%s'\n", argv[i+1]); i++;  }
     if(arg.find("--channel")!=string::npos && i+1<argc)        { histoChannelName.push_back(argv[i+1]);printf("Only process Histo with Channel '%s'\n", argv[i+1]); i++;}
     if(arg.find("--NoPlotwith")!=string::npos && i+1<argc)	{ nohistoNameMask.push_back(argv[i+1]); printf("Only process Histo not containing '%s'\n", argv[i+1]); i++;  }
     if(arg.find("--index"  )!=string::npos && i+1<argc)         { sscanf(argv[i+1],"%d",&cutIndex); i++; onlyCutIndex=(cutIndex>=0); printf("index = %i\n", cutIndex);  }
     if(arg.find("--chi2"  )!=string::npos)                      { showChi2 = true;  }
     if(arg.find("--showUnc") != string::npos) {
       showUnc=true;
       if(i+1<argc) {
	 string nextArg(argv[i+1]);
	 if(nextArg.find("--")==string::npos)
	   {
	     sscanf(argv[i+1],"%lf",&baseRelUnc);
	     i++;
	   }
       }
       printf("Uncertainty band will be included for MC with base relative uncertainty of: %3.2f",baseRelUnc);
     }
     if(arg.find("--isSim")!=string::npos){ isSim = true;    }
     if(arg.find("--isinProgress")!=string::npos){ isinProgress = true;    }
     if(arg.find("--isDataBlind")!=string::npos){ isDataBlind = true; printf("isDataBlind\n");} //RJ
     if(arg.find("--noratio")!=string::npos){ noratio = true; printf("noratio\n");}
     if(arg.find("--nosig")!=string::npos){nosig = true; printf("noSig plot=1\n");}//RJ
     if(arg.find("--noLog")!=string::npos){ noLog = true;    }
     if(arg.find("--useDataMinusMC")!=string::npos){ useDataMinusMC = true; }
     if(arg.find("--no2D"  )!=string::npos){ do2D = false;    }
     if(arg.find("--no1D"  )!=string::npos){ do1D = false;    }
     if(arg.find("--noTex" )!=string::npos){ doTex= false;    }
     if(arg.find("--noRoot")!=string::npos){ StoreInFile = false;    }
     if(arg.find("--noPlot")!=string::npos){ doPlot = false;    }
     if(arg.find("--plotExt" )!=string::npos && i+1<argc){ plotExt   = argv[i+1];  i++;  printf("saving plots as = %s\n", plotExt.c_str());  }
     if(arg.find("--cutflow" )!=string::npos && i+1<argc){ cutflowhisto   = argv[i+1];  i++;  printf("Normalizing from 1st bin in = %s\n", cutflowhisto.c_str());  }
     if(arg.find("--splitCanvas")!=string::npos){ splitCanvas = true;    }
     if(arg.find("--addcutLine1"  )!=string::npos && i+1<argc) { sscanf(argv[i+1],"%f",&cutLine1); i++; printf("cutLine1 = %f\n", cutLine1); } //RJ
     if(arg.find("--addcutLine2"  )!=string::npos && i+1<argc) { sscanf(argv[i+1],"%f",&cutLine2); i++; printf("cutLine2 = %f\n", cutLine2); } //RJ
     if(arg.find("--addleftcutLine"  )!=string::npos && i+1<argc) { sscanf(argv[i+1],"%f",&leftcutLine); i++; printf("leftcutLine = %f\n", leftcutLine); } //RJ
     if(arg.find("--addrightcutLine"  )!=string::npos && i+1<argc) { sscanf(argv[i+1],"%f",&rightcutLine); i++; printf("rightcutLine = %f\n", rightcutLine); } //RJ
     if(arg.find("--ctrlSample" )!=string::npos && i+1<argc){ ctrlSample   = argv[i+1];  i++;  printf("Control Sample = %s\n", ctrlSample.c_str());  }
     if(arg.find("--rebin"  )!=string::npos && i+1<argc) { sscanf(argv[i+1],"%d",&rebin); i++; printf("rebin = %d\n", rebin); } //RJ
     if(arg.find("--scaleSignal")!=string::npos && i+1<argc){ sscanf(argv[i+1],"%lf",&scaleSignal); i++; printf("scaleSignal = %f\n", scaleSignal); }
     if(arg.find("--scaleYMax")  !=string::npos && i+1<argc){ sscanf(argv[i+1],"%lf",&scaleYMax); i++; printf("scaleYMax = %f\n", scaleYMax); }
     if(arg.find("--scaleYMin")  !=string::npos && i+1<argc){ sscanf(argv[i+1],"%lf",&scaleYMin); i++; printf("scaleYMin = %f\n", scaleYMin); }

   }
   system( (string("mkdir -p ") + outDir).c_str());
   if(ctrlSample.find("_")!=std::string::npos)ctrlSample.replace(ctrlSample.find("_"),1," ");

   char buf[255];
   sprintf(buf, "_Index%d", cutIndex);
   cutIndexStr = buf;

   TString nameJsonFile = jsonFile;


   JSONWrapper::Object Root(jsonFile, true);
   checkFileExist(Root,inDir,NameAndType(cutflowhisto,true, false));  //Used to get the rescale factor based on the total number of events geenrated
   std::list<NameAndType> histlist;
   GetListOfObject(Root,inDir,histlist);
   histlist.sort();
   histlist.unique();


   TFile* OutputFile = NULL;
   if(StoreInFile) OutputFile = new TFile(outFile.c_str(),"RECREATE");
   if(!doPlot){
   	printf("\033[32m Progressing Bar              :0%%       20%%       40%%       60%%       80%%       100%%\n \033[0m");
  	printf("\033[32m                             : \033[0m");
   }
   int TreeStep = histlist.size()/50;if(TreeStep==0)TreeStep=1;
   string csvFile(outDir +"/histlist.csv");
   system(("echo \"\" > " + csvFile).c_str());


   int ictr(0);
   for(std::list<NameAndType>::iterator it= histlist.begin(); it!= histlist.end(); it++,ictr++)
     {
       if(!doPlot){ if(ictr%TreeStep==0){printf(".");fflush(stdout);} }
       bool passMasking(false);
       for(unsigned int i=0;i<histoNameMask.size();i++){

	if(histoNameMaskStart.size()==0){
		if(it->name.find(histoNameMask[i])!=std::string::npos)
		passMasking = true;
	} else {
		for(unsigned int j=0;j<histoNameMaskStart.size();j++){

			if(it->name.find(histoNameMask[i])!=std::string::npos &&
				it->name.find(histoNameMaskStart[j])!=std::string::npos) passMasking = true;
		}
	}
       }
       if(histoNameMask.size()==0 && histoNameMaskStart.size()==0)passMasking = true;

	//excluding some hists
       if(passMasking && nohistoNameMask.size()!=0) {
		for(unsigned int j=0;j<nohistoNameMask.size();j++)
		{
			if(it->name.find(nohistoNameMask[j])!=std::string::npos){
				passMasking = false;
				break;}
		}
       }


       if(!passMasking)continue;


       bool skipChannel(true);
       for(size_t j=0; j<histoChannelName.size(); j++)
       {
	if(it->name.find(histoChannelName[j])!=std::string::npos)
		skipChannel=false;
       }
       if(skipChannel) continue;

       system(("echo \"" + it->name + "\" >> " + csvFile).c_str());
       TString itname = it->name;

       if(doPlot && do2D  && !it->type){                      if(!splitCanvas){ Draw2DHistogram(Root,inDir,*it);}else{Draw2DHistogramSplitCanvas(Root,inDir,*it);}}
       if(doPlot && do1D  &&  it->type){                                        Draw1DHistogram(Root,inDir,*it); }

       if(StoreInFile && do2D  && !it->type){                                  	SavingToFile(Root,inDir,*it, OutputFile); }
       if(StoreInFile && do1D  &&  it->type){					SavingToFile(Root,inDir,*it, OutputFile); }

     }


   printf("\n");
   if(StoreInFile) OutputFile->Close();

   system(("generateJSONplotterFromList.py -i " + csvFile + " -o "+outDir+"/plotter.json").c_str());
   //system(("rm " + csvFile).c_str());
   system(("generatehtml.py -o "+outDir).c_str());
   printf("\033[31mYou can browse the results using: firefox %sindex.html\033[0m\n",outDir.c_str());
}

