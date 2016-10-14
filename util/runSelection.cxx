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
#include "TProfile.h"
#include "TEventList.h"
#include "TROOT.h"
#include "TMath.h"
#include "TList.h"

#include "DMHiggsAnalysis/SmartSelectionMonitor.h"
#include "DMHiggsAnalysis/deltaR.h"



using namespace std;

bool symemtricCuts(double met, double ptyy,int cut){

  return   met >= cut*10 && ptyy >= cut*10;

}
bool Categorisation(std::vector<double> variables,int category){


  //  double met = variables[0];
  double metSig = variables[0];
  double ptyy = variables[1];
  double ptHard = variables[2];
  double Nlep = variables[3];


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

  case 5 :

    if(  ptyy > 90  && ( metSig > 7 && Nlep == 0 ) )
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



int main(int argc, char *argv[])
{

    // check arguments
    if(argc<2) {
        std::cout << "Usage : " << argv[0] << " parameters_cfg.py" << std::endl;
        exit(0);
    }

    DataEvtSummaryHandler summaryHandler_;

    HG::Config runProcess;
    HG::StrV Files = HG::parseArguments(&runProcess, argc, argv);

    //tree info
    int evStart     = runProcess.getInt("evStart");
    int evEnd       = runProcess.getInt("evEnd");
    TString dirname = runProcess.getStr("dirName");

    bool isMC       = runProcess.getBool("isMC");
    TString inputUrl = runProcess.getStr("InputFile","InputFile");
    TString SampleName = runProcess.getStr("SampleName","SampleName");

    int mctruthmode = runProcess.getInt("mctruthmode");
    bool isSignal   = runProcess.getBool("issignal");

    bool isMC_H125  = isMC && (SampleName=="SMHiggs");
    bool isMC_gamjet = isMC && (SampleName=="#gamma+jets");
    bool isMC_gamgam =  isMC && (SampleName=="#gamma#gamma");
    bool isMC_vgam = isMC && (SampleName=="V#gamma");



    TFile *inputFile = TFile::Open(inputUrl);
    printf("Looping on %s\n",inputUrl.Data());
    if(inputFile==0) return -1;
    if(inputFile->IsZombie()) return -1;
    if( !summaryHandler_.attachToTree( (TTree *)inputFile->Get(dirname) ) ) {
        inputFile->Close();
        return -1;
    }

    const int totalEntries= summaryHandler_.getEntries();
    cout << "totalEntries: " << totalEntries << endl;
    if(evEnd<0 || evEnd>summaryHandler_.getEntries() ) evEnd=totalEntries;
    if(evStart > evEnd ) {
        inputFile->Close();
        return -1;
    }


    TString outdir=runProcess.getStr("OutputDir");
    TString outUrl( outdir );
    gSystem->Exec("mkdir -p " + outUrl);

    TString outFileUrl(gSystem->BaseName(inputUrl));
    outFileUrl.ReplaceAll(".root","");
    if(mctruthmode!=0) {
        outFileUrl += "_filt";
        outFileUrl += mctruthmode;
    }

    TString outTxtUrl_final= outUrl + "/" + outFileUrl + "_FinalList.txt";
    FILE* outTxtFile_final = NULL;
    outTxtFile_final = fopen(outTxtUrl_final.Data(), "w");
    printf("TextFile URL = %s\n",outTxtUrl_final.Data());
    fprintf(outTxtFile_final,"run event diphotonmass diphotonpt met\n");


    TString outTreeUrl = outUrl + "/";
    outTreeUrl += outFileUrl + "_tree.root";
    bool ifsaveEvents = (!isMC || isSignal || isMC_H125 || isMC_gamgam || isMC_gamjet);
    //always save tree
    ifsaveEvents = true;

    TFile *otreefile=TFile::Open(outTreeUrl, "recreate");

    /*
        TTree *myEvents = new TTree("tree","tree");
        double mgg;
        //double weight_t, weight_all;
        myEvents->Branch("mgg",&mgg,"mgg/D");
        //myEvents->Branch("weight_t", &weight_t,"weight_t/D");
        //myEvents->Branch("weight_all", &weight_all,"weight_all/D");
    */

    TTree *myEvents_bin1 = new TTree("tree_bin1","tree_bin1");
    double mgg_bin1(0.), weight_bin1_t(1.);
    myEvents_bin1->Branch("mgg",&mgg_bin1,"mgg/D");
    myEvents_bin1->Branch("weight_t", &weight_bin1_t,"weight_t/D");


    TTree *myEvents_bin2 = new TTree("tree_bin2","tree_bin2");
    double mgg_bin2(0.), weight_bin2_t(1.);
    myEvents_bin2->Branch("mgg",&mgg_bin2,"mgg/D");
    myEvents_bin2->Branch("weight_t", &weight_bin2_t,"weight_t/D");


    TTree *myEvents_bin3 = new TTree("tree_bin3","tree_bin3");
    double mgg_bin3(0.), weight_bin3_t(1.);
    myEvents_bin3->Branch("mgg",&mgg_bin3,"mgg/D");
    myEvents_bin3->Branch("weight_t", &weight_bin3_t,"weight_t/D");


    TTree *myEvents_bin4 = new TTree("tree_bin4","tree_bin4");
    double mgg_bin4(0.), weight_bin4_t(1.);
    myEvents_bin4->Branch("mgg",&mgg_bin4,"mgg/D");
    myEvents_bin4->Branch("weight_t", &weight_bin4_t,"weight_t/D");


    TTree *myEvents_bin5 = new TTree("tree_bin5","tree_bin5");
    double mgg_bin5(0.), weight_bin5_t(1.);
    myEvents_bin5->Branch("mgg",&mgg_bin5,"mgg/D");
    myEvents_bin5->Branch("weight_t", &weight_bin5_t,"weight_t/D");


    TTree *myEvents_bin6 = new TTree("tree_bin6","tree_bin6");
    double mgg_bin6(0.), weight_bin6_t(1.);
    myEvents_bin6->Branch("mgg",&mgg_bin6,"mgg/D");
    myEvents_bin6->Branch("weight_t", &weight_bin6_t,"weight_t/D");


    TTree *myEvents_bin7 = new TTree("tree_bin7","tree_bin7");
    double mgg_bin7(0.), weight_bin7_t(1.);
    myEvents_bin7->Branch("mgg",&mgg_bin7,"mgg/D");
    myEvents_bin7->Branch("weight_t", &weight_bin7_t,"weight_t/D");


    TTree *myEvents_bin8 = new TTree("tree_bin8","tree_bin8");
    double mgg_bin8(0.), weight_bin8_t(1.);
    myEvents_bin8->Branch("mgg",&mgg_bin8,"mgg/D");
    myEvents_bin8->Branch("weight_t", &weight_bin8_t,"weight_t/D");


    TTree *myEvents_bin9 = new TTree("tree_bin9","tree_bin9");
    double mgg_bin9(0.), weight_bin9_t(1.);
    myEvents_bin9->Branch("mgg",&mgg_bin9,"mgg/D");
    myEvents_bin9->Branch("weight_t", &weight_bin9_t,"weight_t/D");


    TTree *myEvents_bin10 = new TTree("tree_bin10","tree_bin10");
    double mgg_bin10(0.), weight_bin10_t(1.);
    myEvents_bin10->Branch("mgg",&mgg_bin10,"mgg/D");
    myEvents_bin10->Branch("weight_t", &weight_bin10_t,"weight_t/D");


    TTree *myEvents_bin11 = new TTree("tree_bin11","tree_bin11");
    double mgg_bin11(0.), weight_bin11_t(1.);
    myEvents_bin11->Branch("mgg",&mgg_bin11,"mgg/D");
    myEvents_bin11->Branch("weight_t", &weight_bin11_t,"weight_t/D");


    TTree *myEvents_bin12 = new TTree("tree_bin12","tree_bin12");
    double mgg_bin12(0.), weight_bin12_t(1.);
    myEvents_bin12->Branch("mgg",&mgg_bin12,"mgg/D");
    myEvents_bin12->Branch("weight_t", &weight_bin12_t,"weight_t/D");


    TTree *myEvents_bin13 = new TTree("tree_bin13","tree_bin13");
    double mgg_bin13(0.), weight_bin13_t(1.);
    myEvents_bin13->Branch("mgg",&mgg_bin13,"mgg/D");
    myEvents_bin13->Branch("weight_t", &weight_bin13_t,"weight_t/D");


    TTree *myEvents_bin14 = new TTree("tree_bin14","tree_bin14");
    double mgg_bin14(0.), weight_bin14_t(1.);
    myEvents_bin14->Branch("mgg",&mgg_bin14,"mgg/D");
    myEvents_bin14->Branch("weight_t", &weight_bin14_t,"weight_t/D");


    TTree *myEvents_bin15 = new TTree("tree_bin15","tree_bin15");
    double mgg_bin15(0.), weight_bin15_t(1.);
    myEvents_bin15->Branch("mgg",&mgg_bin15,"mgg/D");
    myEvents_bin15->Branch("weight_t", &weight_bin15_t,"weight_t/D");


    TTree *myEvents_bin16 = new TTree("tree_bin16","tree_bin16");
    double mgg_bin16(0.), weight_bin16_t(1.);
    myEvents_bin16->Branch("mgg",&mgg_bin16,"mgg/D");
    myEvents_bin16->Branch("weight_t", &weight_bin16_t,"weight_t/D");



    TTree *myEvents_allbins = new TTree("tree_allbins","tree_allbins");
    double mgg_allbins(0.), weight_allbins_t(1.);
    int mycat_idx(-1);
    myEvents_allbins->Branch("mgg",&mgg_allbins,"mgg/D");
    myEvents_allbins->Branch("weight_t", &weight_allbins_t,"weight_t/D");
    myEvents_allbins->Branch("mycat", &mycat_idx,"mycat_idx/I");


    //##################################################################################
    //##########################    INITIATING HISTOGRAMS     ##########################
    //##################################################################################

    SmartSelectionMonitor mon;


    mon.addHistogram( new TH1F( "pu_sel", ";Pileup;Events", 50,0,50) );
    mon.addHistogram( new TH1F( "puwgt_sel", ";Pileup;Events", 50,0,50) );

    mon.addHistogram( new TH1F( "nvtx_sel",     ";Vertices;Events",50,0,50) );
    mon.addHistogram( new TH1F( "nvtxwgt_sel",  ";Vertices;Events",50,0,50) );

    //for MC normalization (to 1/pb)
    TH1F* Hcutflow  = (TH1F*) mon.addHistogram(  new TH1F ("cutflow"    , "cutflow"    ,6,0,6) ) ;
    TH1F* Heventflow  = (TH1F*) mon.addHistogram(  new TH1F ("eventflow"    , "eventflow"    ,15+6+2,0,15+6+2) ) ;
    TH1F* Heventflow_raw  = (TH1F*) mon.addHistogram(  new TH1F ("eventflow_raw"    , "eventflow_raw"    ,15+6+2,0,15+6+2) ) ;



    double PTBins[]= {0,10,20,30,40,50,60,70,80,90,100,120,140,160,180,200,250,300,350,400,600};
    const int nPTBins = sizeof(PTBins)/sizeof(double) - 1;

    double METSigBins[]= {0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.2,2.4,2.6,2.8,3.0,3.2,3.4,3.6,3.8,4.4,4.8,5.2,5.6,6,7,8,10};
    const int nMETSigBins = sizeof(METSigBins)/sizeof(double) - 1;

    mon.addHistogram( new TH1F( "diphoton_mass_sel",    ";#it{m}_{#gamma#gamma} [GeV];Events", 62,105,160) );
    mon.addHistogram( new TH1F( "diphoton_pt_sel",    ";#it{p}_{T}^{#gamma#gamma} [GeV];Events", 100,0,600) );
    mon.addHistogram( new TH1F( "diphoton_pt_rebin_sel", ";#it{p}_{T}^{#gamma#gamma} [GeV];Events", nPTBins, PTBins) );
    mon.addHistogram( new TH1F( "met_sel",    	      ";#it{E}_{T}^{miss} [GeV];Events", 100,0,600) );
    mon.addHistogram( new TH1F( "sumet_sel",            ";sumE_{T} [GeV];Events", 100,50,1050) );
    mon.addHistogram( new TH1F( "metsig_sel",         ";E_{T}^{miss} Significance [GeV];Events", 100,0,10) );
    mon.addHistogram( new TH1F( "metsig2_sel",         ";E_{T}^{miss} Significance [GeV];Events", 40,0,20) );
    mon.addHistogram( new TH1F( "metsig_selbeg6",         ";E_{T}^{miss} Significance [GeV];Events", 1,0,1) );
    mon.addHistogram( new TH1F( "metsig_selbeg7",         ";E_{T}^{miss} Significance [GeV];Events", 1,0,1) );

    mon.addHistogram( new TH1F( "metsig_rebin_sel",         ";E_{T}^{miss} Significance [GeV];Events", nMETSigBins,METSigBins) );
    mon.addHistogram( new TH1F( "balancedif_sel",     ";|E_{T}^{miss}-#it{p}_{T}^{#gamma#gamma}|/#it{p}_{T}^{#gamma#gamma};Events", 5,0,1.0) );
    mon.addHistogram( new TH1F( "metphi_sel",         ";#phi(#it{E}_{T}^{miss}) [rad];Events", 50,-TMath::Pi(),TMath::Pi()) );
    mon.addHistogram( new TH1F( "pthard_sel",         ";#it{p}_{T}^{hard} [GeV];Events", 100,0,600) );
    mon.addHistogram( new TH1F( "dphiGamGamMET_sel",   ";#Delta#it{#phi}(#it{p}_{T}^{#gamma#gamma},E_{T}^{miss}) [rad];Events", 50,0,TMath::Pi()) );

    mon.addHistogram( new TH1F( "leadingjetpt_sel",         ";#it{p}_{T}^{j} [GeV];Events", nPTBins, PTBins) );
    mon.addHistogram( new TH1F( "subleadingjetpt_sel",         ";#it{p}_{T}^{j} [GeV];Events", nPTBins, PTBins) );



    mon.addHistogram( new TH1F( "leadingpho_pt_sel",    ";Leading #it{p}_{T}^{#gamma} [GeV];Events", 100,0,600) );
    mon.addHistogram( new TH1F( "trailingpho_pt_sel",    ";Trailing #it{p}_{T}^{#gamma} [GeV];Events", 100,0,600) );
    mon.addHistogram( new TH1F( "leadingpho_eta_sel",    ";Leading #it{#eta}^{#gamma};Events", 50,-2.4,2.4) );
    mon.addHistogram( new TH1F( "trailingpho_eta_sel",    ";Trailing #it{#eta}^{#gamma};Events", 50,-2.4,2.4) );
    mon.addHistogram( new TH1F( "leadingpho_phi_sel",    ";Leading #it{#phi}^{#gamma} [rad];Events", 50,-1.*TMath::Pi(),TMath::Pi()) );
    mon.addHistogram( new TH1F( "trailingpho_phi_sel",    ";Trailing #it{#phi}^{#gamma} [rad];Events", 50,-1.*TMath::Pi(),TMath::Pi()) );

    mon.addHistogram( new TH1F( "diphoton_mass_3jetctrl",    ";#it{m}_{#gamma#gamma} [GeV];Events", 62,105,160) );

    // mon.addHistogram( new TH1F( "dphiGamGamMET_bin1",   ";#Delta#it{#phi}(#it{p}_{T}^{#gamma#gamma},E_{T}^{miss}) [rad];Events", 50,0,TMath::Pi()) );
    // mon.addHistogram( new TH1F( "metsig_bin1",         ";E_{T}^{miss} Significance [GeV];Events", 100,0,10) );
    // mon.addHistogram( new TH1F( "balancedif_bin1",     ";|E_{T}^{miss}-#it{p}_{T}^{#gamma#gamma}|/#it{p}_{T}^{#gamma#gamma};Events", 5,0,1.0) );

    // mon.addHistogram( new TH1F( "diphoton_mass_bin1",    ";#it{m}_{#gamma#gamma} [GeV];Events", 62,105,160) );
    // mon.addHistogram( new TH1F( "diphoton_mass_bin2",    ";#it{m}_{#gamma#gamma} [GeV];Events", 62,105,160) );
    // mon.addHistogram( new TH1F( "diphoton_mass_bin3",    ";#it{m}_{#gamma#gamma} [GeV];Events", 62,105,160) );
    // mon.addHistogram( new TH1F( "diphoton_mass_bin4",    ";#it{m}_{#gamma#gamma} [GeV];Events", 62,105,160) );

    // mon.addHistogram( new TH1F( "yields_diphoton_mass_bin0",    ";#it{m}_{#gamma#gamma} [GeV];Events", 1,0,1) );
    // mon.addHistogram( new TH1F( "yields_diphoton_mass_bin1",    ";#it{m}_{#gamma#gamma} [GeV];Events", 1,0,1) );
    // mon.addHistogram( new TH1F( "yields_diphoton_mass_bin2",    ";#it{m}_{#gamma#gamma} [GeV];Events", 1,0,1) );
    // mon.addHistogram( new TH1F( "yields_diphoton_mass_bin3",    ";#it{m}_{#gamma#gamma} [GeV];Events", 1,0,1) );
    // mon.addHistogram( new TH1F( "yields_diphoton_mass_bin4",    ";#it{m}_{#gamma#gamma} [GeV];Events", 1,0,1) );


    mon.addHistogram( new TH1F( "met_mu0_10_sel",            ";#it{E}_{T}^{miss} [GeV];Events", 100,0,600) );
    mon.addHistogram( new TH1F( "sumet_mu0_10_sel",            ";sumE_{T} [GeV];Events", 100,50,1050) );
    mon.addHistogram( new TH1F( "metsig_mu0_10_sel",         ";E_{T}^{miss} Significance [GeV];Events", 100,0,10) );
    mon.addHistogram( new TH1F( "met_mu10_20_sel",            ";#it{E}_{T}^{miss} [GeV];Events", 100,0,600) );
    mon.addHistogram( new TH1F( "sumet_mu10_20_sel",            ";sumE_{T} [GeV];Events", 100,50,1050) );
    mon.addHistogram( new TH1F( "metsig_mu10_20_sel",         ";E_{T}^{miss} Significance [GeV];Events", 100,0,10) );
    mon.addHistogram( new TH1F( "met_mu20_30_sel",            ";#it{E}_{T}^{miss} [GeV];Events", 100,0,600) );
    mon.addHistogram( new TH1F( "sumet_mu20_30_sel",            ";sumE_{T} [GeV];Events", 100,50,1050) );
    mon.addHistogram( new TH1F( "metsig_mu20_30_sel",         ";E_{T}^{miss} Significance [GeV];Events", 100,0,10) );




    //mass vs METSig cut
    mon.addHistogram( new TH1F( "diphoton_mass_metsigbeg1",    ";#it{m}_{#gamma#gamma} [GeV];Events", 62,105,160) );
    mon.addHistogram( new TH1F( "diphoton_mass_metsigbeg2",    ";#it{m}_{#gamma#gamma} [GeV];Events", 62,105,160) );
    mon.addHistogram( new TH1F( "diphoton_mass_metsigbeg3",    ";#it{m}_{#gamma#gamma} [GeV];Events", 62,105,160) );
    mon.addHistogram( new TH1F( "diphoton_mass_metsigbeg4",    ";#it{m}_{#gamma#gamma} [GeV];Events", 62,105,160) );
    mon.addHistogram( new TH1F( "diphoton_mass_metsigbeg5",    ";#it{m}_{#gamma#gamma} [GeV];Events", 62,105,160) );
    mon.addHistogram( new TH1F( "diphoton_mass_metsigbeg6",    ";#it{m}_{#gamma#gamma} [GeV];Events", 62,105,160) );

    //pt vs METSig cut
    mon.addHistogram( new TH1F( "diphoton_pt_metsigbeg1",    ";#it{p}_{T}^{#gamma#gamma} [GeV];Events", 100,0,600) );
    mon.addHistogram( new TH1F( "diphoton_pt_metsigbeg2",    ";#it{p}_{T}^{#gamma#gamma} [GeV];Events", 100,0,600) );
    mon.addHistogram( new TH1F( "diphoton_pt_metsigbeg3",    ";#it{p}_{T}^{#gamma#gamma} [GeV];Events", 100,0,600) );
    mon.addHistogram( new TH1F( "diphoton_pt_metsigbeg4",    ";#it{p}_{T}^{#gamma#gamma} [GeV];Events", 100,0,600) );
    mon.addHistogram( new TH1F( "diphoton_pt_metsigbeg5",    ";#it{p}_{T}^{#gamma#gamma} [GeV];Events", 100,0,600) );
    mon.addHistogram( new TH1F( "diphoton_pt_metsigbeg6",    ";#it{p}_{T}^{#gamma#gamma} [GeV];Events", 100,0,600) );


    mon.addHistogram( new TH2F( "pv_dist_vs_met_sig",    ";z^{hardest}_{PV} - z^{#gamma#gamma}_{PV} (mm);E^{miss}_{T} significance;Events", 100,-1.5,1.5,20,0,20) );
    mon.addHistogram( new TH1F( "pv_dist",    ";z^{hardest}_{PV} - z^{#gamma#gamma}_{PV} (mm);E^{miss}_{T} significance;Events", 100,-1.5,1.5) );
    mon.addHistogram( new TH1F( "pv_dist_metsigbeg1",    ";z^{hardest}_{PV} - z^{#gamma#gamma}_{PV} (mm);Events", 100,-1.5,1.5) );
    mon.addHistogram( new TH1F( "pv_dist_metsigbeg2",    ";z^{hardest}_{PV} - z^{#gamma#gamma}_{PV} (mm);Events", 100,-1.5,1.5) );
    mon.addHistogram( new TH1F( "pv_dist_metsigbeg3",    ";z^{hardest}_{PV} - z^{#gamma#gamma}_{PV} (mm);Events", 100,-1.5,1.5) );
    mon.addHistogram( new TH1F( "pv_dist_metsigbeg4",    ";z^{hardest}_{PV} - z^{#gamma#gamma}_{PV} (mm);Events", 100,-1.5,1.5) );
    mon.addHistogram( new TH1F( "pv_dist_metsigbeg5",    ";z^{hardest}_{PV} - z^{#gamma#gamma}_{PV} (mm);Events", 100,-1.5,1.5) );
    mon.addHistogram( new TH1F( "pv_dist_metsigbeg6",    ";z^{hardest}_{PV} - z^{#gamma#gamma}_{PV} (mm);Events", 100,-1.5,1.5) );

 
    mon.addHistogram( new TH2F( "allpart_phi_vs_met_sig",    ";|#Delta#phi(hardObj,E^{miss}_{T})|;E^{miss}_{T} significance;Events", 50,0,TMath::Pi(),20,0,20) );
    mon.addHistogram( new TH1F( "allpart_phi",    ";|#Delta#phi(hardObj,E^{miss}_{T})|;Events", 50,0,TMath::Pi()) );
    mon.addHistogram( new TH1F( "allpart_phi_metsigbeg1",    ";|#Delta#phi(hardObj,E^{miss}_{T})|;Events", 50,0,TMath::Pi()) );
    mon.addHistogram( new TH1F( "allpart_phi_metsigbeg2",    ";|#Delta#phi(hardObj,E^{miss}_{T})|;Events", 50,0,TMath::Pi()) );
    mon.addHistogram( new TH1F( "allpart_phi_metsigbeg3",    ";|#Delta#phi(hardObj,E^{miss}_{T})|;Events", 50,0,TMath::Pi()) );
    mon.addHistogram( new TH1F( "allpart_phi_metsigbeg4",    ";|#Delta#phi(hardObj,E^{miss}_{T})|;Events", 50,0,TMath::Pi()) );
    mon.addHistogram( new TH1F( "allpart_phi_metsigbeg5",    ";|#Delta#phi(hardObj,E^{miss}_{T})|;Events", 50,0,TMath::Pi()) );
    mon.addHistogram( new TH1F( "allpart_phi_metsigbeg6",    ";|#Delta#phi(hardObj,E^{miss}_{T})|;Events", 50,0,TMath::Pi()) );



    for( unsigned int i = 0 ; i < 6 ; ++i){

      mon.addHistogram( new TH1F( "allpart_phi_cat"+TString::Format("%d",i),    ";|#Delta#phi(hardObj,E^{miss}_{T})|;Events", 50,0,TMath::Pi()) );
      mon.addHistogram( new TH1F( "pv_dist_cat"+TString::Format("%d",i),    ";z^{hardest}_{PV} - z^{#gamma#gamma}_{PV} (mm);Events", 100,-1.5,1.5) );
      mon.addHistogram( new TH1F( "diphoton_mass_bin"+TString::Format("%d",i),    ";#it{m}_{#gamma#gamma} [GeV];Events", 62,105,160) );
      mon.addHistogram( new TH1F( "dphiGamGamMET_bin"+TString::Format("%d",i),   ";#Delta#it{#phi}(#it{p}_{T}^{#gamma#gamma},E_{T}^{miss}) [rad];Events", 50,0,TMath::Pi()) );
      mon.addHistogram( new TH1F( "metsig_bin"+TString::Format("%d",i),         ";E_{T}^{miss} Significance [GeV];Events", 100,0,10) );
      mon.addHistogram( new TH1F( "balancedif_bin"+TString::Format("%d",i),     ";|E_{T}^{miss}-#it{p}_{T}^{#gamma#gamma}|/#it{p}_{T}^{#gamma#gamma};Events", 5,0,1.0) );
      mon.addHistogram( new TH1F( "yields_diphoton_mass_bin"+TString::Format("%d",i),    ";#it{m}_{#gamma#gamma} [GeV];Events", 1,0,1) );


    }






    //optimization
    std::vector<double> optim_Cuts_METSig;
    std::vector<double> optim_Cuts_MET;
    std::vector<double> optim_Cuts_DphiyyMET;
    std::vector<double> optim_Cuts_yypt;


    for(double metsig=0; metsig<=10; metsig+=0.5) {
        optim_Cuts_METSig     .push_back(metsig);
    }
    for(double met=0; met<=200; met+=10) {
        optim_Cuts_MET     .push_back(met);
    }
    for(double dphi=0; dphi<=3.1; dphi+=0.31) {
        optim_Cuts_DphiyyMET  .push_back(dphi);
    }
    for(double yypt=0; yypt<=150; yypt+=10) {
        optim_Cuts_yypt  .push_back(yypt);
    }

    size_t nOptims_METSig = optim_Cuts_METSig.size();
    size_t nOptims_MET = optim_Cuts_MET.size();
    size_t nOptims_DphiyyMET = optim_Cuts_DphiyyMET.size();
    size_t nOptims_yypt = optim_Cuts_yypt.size();

    TH2F *h2 = (TH2F*) mon.addHistogram( new TH2F( "yields_dphi_vs_metsig",";E_{T}^{miss} Significance; #Delta#it{#phi}(#it{p}_{T}^{#gamma#gamma},E_{T}^{miss}) [rad];Events", nOptims_METSig,0,nOptims_METSig, nOptims_DphiyyMET,0,nOptims_DphiyyMET) );
    TH2F *h3 = (TH2F*) mon.addHistogram( new TH2F( "yields_dphi_vs_met",";E_{T}^{miss} [GeV]; #Delta#it{#phi}(#it{p}_{T}^{#gamma#gamma},E_{T}^{miss}) [rad];Events", nOptims_MET,0,nOptims_MET, nOptims_DphiyyMET,0,nOptims_DphiyyMET) );
    TH2F *h4 = (TH2F*) mon.addHistogram( new TH2F( "yields_metsig_vs_yypt", ";#it{p}_{T}^{#gamma#gamma}; E_{T}^{miss} Significance; Events", nOptims_yypt,0,nOptims_yypt, nOptims_METSig,0,nOptims_METSig) );

    for(int ibin=1; ibin<=h2->GetXaxis()->GetNbins(); ibin++) {
        if(ibin%4 != 0 ) continue;
        TString label("");
        label +="> ";
        std::ostringstream buff;
        buff << fixed  << setprecision(0) << optim_Cuts_METSig[ibin-1];
        std::string str = buff.str();
        TString title = str;
        label += title;
        h2->GetXaxis()->SetBinLabel(ibin,label);
    }
    for(int ibin=1; ibin<=h2->GetYaxis()->GetNbins(); ibin++) {
        if(ibin%2 !=0 ) continue;
        TString label("");
        label +="> ";
        std::ostringstream buff;
        buff << fixed  << setprecision(1) << optim_Cuts_DphiyyMET[ibin-1];
        std::string str = buff.str();
        TString title = str;
        label += title;
        h2->GetYaxis()->SetBinLabel(ibin,label);
    }


    for(int ibin=1; ibin<=h3->GetXaxis()->GetNbins(); ibin++) {
        if(ibin%4 != 0 ) continue;
        TString label("");
        label +="> ";
        std::ostringstream buff;
        buff << fixed  << setprecision(0) << optim_Cuts_MET[ibin-1];
        std::string str = buff.str();
        TString title = str;
        label += title;
        h3->GetXaxis()->SetBinLabel(ibin,label);
    }
    for(int ibin=1; ibin<=h3->GetYaxis()->GetNbins(); ibin++) {
        if(ibin%2 !=0 ) continue;
        TString label("");
        label +="> ";
        std::ostringstream buff;
        buff << fixed  << setprecision(1) << optim_Cuts_DphiyyMET[ibin-1];
        std::string str = buff.str();
        TString title = str;
        label += title;
        h3->GetYaxis()->SetBinLabel(ibin,label);
    }


    for(int ibin=1; ibin<=h4->GetXaxis()->GetNbins(); ibin++) {
        if(ibin%4 != 0 ) continue;
        TString label("");
        label +="> ";
        std::ostringstream buff;
        buff << fixed  << setprecision(0) << optim_Cuts_yypt[ibin-1];
        std::string str = buff.str();
        TString title = str;
        label += title;
        h4->GetXaxis()->SetBinLabel(ibin,label);
    }
    for(int ibin=1; ibin<=h4->GetYaxis()->GetNbins(); ibin++) {
        if(ibin%2 !=0 ) continue;
        TString label("");
        label +="> ";
        std::ostringstream buff;
        buff << fixed  << setprecision(1) <<optim_Cuts_METSig[ibin-1];
        std::string str = buff.str();
        TString title = str;
        label += title;
        h4->GetYaxis()->SetBinLabel(ibin,label);
    }







    TH1F *h1 = (TH1F*) mon.addHistogram( new TH1F( "nphotons_sel", ";Photon multiplicity;Events", 3,2,5) );
    for(int ibin=1; ibin<=h1->GetXaxis()->GetNbins(); ibin++) {
        TString label("");
        label +="= ";
        label += (ibin+1);
        h1->GetXaxis()->SetBinLabel(ibin,label);
    }


    h1 = (TH1F*) mon.addHistogram( new TH1F( "nelectrons_sel", ";Electron multiplicity;Events", 3,0,3) );
    for(int ibin=1; ibin<=h1->GetXaxis()->GetNbins(); ibin++) {
        TString label("");
        label +="= ";
        label += (ibin-1);
        h1->GetXaxis()->SetBinLabel(ibin,label);
    }


    h1 = (TH1F*) mon.addHistogram( new TH1F( "nmuons_sel", ";Muon multiplicity;Events", 3,0,3) );
    for(int ibin=1; ibin<=h1->GetXaxis()->GetNbins(); ibin++) {
        TString label("");
        label +="= ";
        label += (ibin-1);
        h1->GetXaxis()->SetBinLabel(ibin,label);
    }



    h1 = (TH1F*) mon.addHistogram( new TH1F( "njets_sel",";Jet multiplicity;Events", 5,0,5) );
    for(int ibin=1; ibin<=h1->GetXaxis()->GetNbins(); ibin++) {
        TString label("");
        label +="= ";
        label += (ibin-1);
        h1->GetXaxis()->SetBinLabel(ibin,label);
    }

    h1 = (TH1F*) mon.addHistogram( new TH1F( "yields_finalsel",";;Events", 17,0,17) );
    h1->GetXaxis()->SetBinLabel(1,"High #it{E}_{T}^{miss}, high #it{p}_{T}^{#gamma#gamma}");
    h1->GetXaxis()->SetBinLabel(2,"High #it{E}_{T}^{miss}, low #it{p}_{T}^{#gamma#gamma}");
    h1->GetXaxis()->SetBinLabel(3,"Intermediate #it{E}_{T}^{miss}");
    h1->GetXaxis()->SetBinLabel(4,"Rest category");
    //model independent limits
    h1->GetXaxis()->SetBinLabel(5,"#it{E}_{T}^{miss} > 10, #it{p}_{T}^{#gamma#gamma} > 10");
    h1->GetXaxis()->SetBinLabel(6,"#it{E}_{T}^{miss} > 20, #it{p}_{T}^{#gamma#gamma} > 20");
    h1->GetXaxis()->SetBinLabel(7,"#it{E}_{T}^{miss} > 30, #it{p}_{T}^{#gamma#gamma} > 30");
    h1->GetXaxis()->SetBinLabel(8,"#it{E}_{T}^{miss} > 40, #it{p}_{T}^{#gamma#gamma} > 40");
    h1->GetXaxis()->SetBinLabel(9,"#it{E}_{T}^{miss} > 50, #it{p}_{T}^{#gamma#gamma} > 50");
    h1->GetXaxis()->SetBinLabel(10,"#it{E}_{T}^{miss} > 60, #it{p}_{T}^{#gamma#gamma} > 60");
    h1->GetXaxis()->SetBinLabel(11,"#it{E}_{T}^{miss} > 70, #it{p}_{T}^{#gamma#gamma} > 70");
    h1->GetXaxis()->SetBinLabel(12,"#it{E}_{T}^{miss} > 80, #it{p}_{T}^{#gamma#gamma} > 80");
    h1->GetXaxis()->SetBinLabel(13,"#it{E}_{T}^{miss} > 90, #it{p}_{T}^{#gamma#gamma} > 90");
    h1->GetXaxis()->SetBinLabel(14,"#it{E}_{T}^{miss} > 100, #it{p}_{T}^{#gamma#gamma} > 100");
    h1->GetXaxis()->SetBinLabel(15,"#it{E}_{T}^{miss} > 110, #it{p}_{T}^{#gamma#gamma} > 110");
    h1->GetXaxis()->SetBinLabel(16,"#it{E}_{T}^{miss} > 120, #it{p}_{T}^{#gamma#gamma} > 120");
    h1->GetXaxis()->SetBinLabel(17,"High #it{E}_{T}^{miss}, high #it{p}_{T}^{#gamma#gamma}, N_{e^{-}} = 0");

    h1 = (TH1F*) mon.addHistogram( new TH1F( "yields_finalsel5cat",";;Events", 5,0,5) );
    h1->GetXaxis()->SetBinLabel(1,"High #it{E}_{T}^{miss}, high #it{p}_{T}^{#gamma#gamma}");
    h1->GetXaxis()->SetBinLabel(2,"High #it{E}_{T}^{miss}, low #it{p}_{T}^{#gamma#gamma}");
    h1->GetXaxis()->SetBinLabel(3,"Intermediate #it{E}_{T}^{miss}");
    h1->GetXaxis()->SetBinLabel(4,"Rest category");
    h1->GetXaxis()->SetBinLabel(5,"High #it{E}_{T}^{miss}, high #it{p}_{T}^{#gamma#gamma}, N_{e^{-}} = 0");



    //Gen plot
    TH1F* Nevent_gen = (TH1F*) mon.addHistogram( new TH1F( "Nevent_gen",            ";Gen;Events", 12,0,12) );
    Nevent_gen->GetXaxis()->SetBinLabel(1,"Total event");
    Nevent_gen->GetXaxis()->SetBinLabel(2,"Acceptance Bin 1");
    Nevent_gen->GetXaxis()->SetBinLabel(3,"Efficiency Bin 1");
    Nevent_gen->GetXaxis()->SetBinLabel(4,"Acceptance Bin 2");
    Nevent_gen->GetXaxis()->SetBinLabel(5,"Efficiency Bin 2");
    Nevent_gen->GetXaxis()->SetBinLabel(6,"Acceptance Bin 3");
    Nevent_gen->GetXaxis()->SetBinLabel(7,"Efficiency Bin 3");
    Nevent_gen->GetXaxis()->SetBinLabel(8,"Acceptance Bin 4");
    Nevent_gen->GetXaxis()->SetBinLabel(9,"Efficiency Bin 4");
    Nevent_gen->GetXaxis()->SetBinLabel(10,"Acceptance Bin 5");
    Nevent_gen->GetXaxis()->SetBinLabel(11,"Efficiency Bin 5");






    mon.addHistogram( new TH1F( "diphoton_mass_gen",    ";#it{m}_{#gamma#gamma} [GeV];Events", 62,105,160) );
    mon.addHistogram( new TH1F( "diphoton_pt_gen",    ";#it{p}_{T}^{#gamma#gamma} [GeV];Events", 100,0,600) );
    mon.addHistogram( new TH1F( "met_gen",            ";#it{E}_{T}^{miss} [GeV];Events", 100,0,600) );
    mon.addHistogram( new TH1F( "met_gen_dmPart",            ";#it{E}_{T}^{miss} [GeV];Events", 100,0,600) );
    mon.addHistogram( new TH1F( "sumet_gen",            ";sumE_{T} [GeV];Events", 100,50,1050) );
    mon.addHistogram( new TH1F( "metsig_gen",         ";E_{T}^{miss} Significance [GeV];Events", 100,0,10) );


    mon.addHistogram( new TH1F( "Nevent",            ";;Events", 1,0,1) );





    //MC normalization (to 1/pb)
    double sumInitialEvents=1.0;
    double skim_eff=1.0;
    if(isMC && !inputUrl.Contains("_reweighting") ) {

      TList* listObjInFile = inputFile->GetListOfKeys();
      bool match = false;
      int countHistos = 0 ;  
      for( TObject* obj : *listObjInFile){
	TString objName = obj->GetName();
	if( objName.Contains("CutFlow_noDalitz_weighted")){
	  std::cout << "Input file contains CutFlow_noDalitz_weighted" << std::endl;
	  ++countHistos;
	}else if( objName.Contains("CutFlow_weighted")){
	  std::cout << "Input file contains CutFlow_weighted" << std::endl;
	  ++countHistos;
	}

      }
	if( countHistos == 2 )
	  match = true;
	else
	  std::cout << "File doens't contain needed histograms. Please check your input File production." << std::endl;

	TH1F* H_noDalitz_weighted = match ? (TH1F *) inputFile->Get("CutFlow_noDalitz_weighted") : nullptr;
	TH1F* H_weighted 	  = match ? (TH1F *) inputFile->Get("CutFlow_weighted") : nullptr ;
	sumInitialEvents = match ? H_noDalitz_weighted->GetBinContent(3) : 1.0 ;

        // Hard-coding to bin number 1,2
	double NxAOD      = match ? H_weighted->GetBinContent(1) : 1;
	double NDxAOD     = match ? H_weighted->GetBinContent(2) : 1;
	skim_eff = NDxAOD / NxAOD;

        //printf("sumInitialEvents = %f, DxAOD skimming efficiency = %f\n",sumInitialEvents, skim_eff);


        for(int bin=1; bin < Heventflow->GetNbinsX()+1; bin++) {
	  if(bin<16) {
	    Heventflow -> SetBinContent(bin, match ? H_noDalitz_weighted->GetBinContent(bin) : 1);
	    //  Heventflow -> SetBinContent(bin, H_weighted->GetBinContent(bin));
	    Heventflow -> GetXaxis()->SetBinLabel(bin,match ? H_noDalitz_weighted->GetXaxis()->GetBinLabel(bin) : "");
	    //  Heventflow -> GetXaxis()->SetBinLabel(bin,H_weighted->GetXaxis()->GetBinLabel(bin));

	    Heventflow_raw -> SetBinContent(bin, match ? H_noDalitz_weighted->GetBinContent(bin) : 1);
	    Heventflow_raw -> GetXaxis()->SetBinLabel(bin,match ? H_noDalitz_weighted->GetXaxis()->GetBinLabel(bin) : "" );
	    //  Heventflow_raw -> SetBinContent(bin, H_weighted->GetBinContent(bin));
	    //  Heventflow_raw -> GetXaxis()->SetBinLabel(bin,H_weighted->GetXaxis()->GetBinLabel(bin));
	  }
        }

        Heventflow->GetXaxis()->SetBinLabel(16+1,"m_{#gamma#gamma}");
        Heventflow->GetXaxis()->SetBinLabel(16+2,"High #it{E}_{T}^{miss}, high #it{p}_{T}^{#gamma#gamma}");
        Heventflow->GetXaxis()->SetBinLabel(17+2,"High #it{E}_{T}^{miss}, low #it{p}_{T}^{#gamma#gamma}");
        Heventflow->GetXaxis()->SetBinLabel(18+2,"Intermediate #it{E}_{T}^{miss}");
        Heventflow->GetXaxis()->SetBinLabel(19+2,"Rest category");

        Heventflow_raw->GetXaxis()->SetBinLabel(16+1,"m_{#gamma#gamma}");
        Heventflow_raw->GetXaxis()->SetBinLabel(16+2,"High #it{E}_{T}^{miss}, high #it{p}_{T}^{#gamma#gamma}");
        Heventflow_raw->GetXaxis()->SetBinLabel(17+2,"High #it{E}_{T}^{miss}, low #it{p}_{T}^{#gamma#gamma}");
        Heventflow_raw->GetXaxis()->SetBinLabel(18+2,"Intermediate #it{E}_{T}^{miss}");
        Heventflow_raw->GetXaxis()->SetBinLabel(19+2,"Rest category");
        //
        if(isSignal) Nevent_gen -> SetBinContent(1,match ?  H_noDalitz_weighted->GetBinContent(1) : 1 );
        //if(isSignal) Nevent_gen -> SetBinContent(1, H_weighted->GetBinContent(1));
      }

      Hcutflow->SetBinContent(1,sumInitialEvents/skim_eff);

      //
      TFile *GamjetsWeights_File = TFile::Open("/afs/cern.ch/work/a/alopezso/private/DMAnalysis20.7/DMHiggsAnalysis/data/etmissSigWeights.root");
      TH1F* h_gamjetweights = NULL;
      if(isMC_gamjet)  h_gamjetweights = (TH1F *) GamjetsWeights_File->Get("etmissSigWeights");

      //
    TFile *GamGamWeights_File = TFile::Open("/afs/cern.ch/work/a/alopezso/private/DMAnalysis20.7/DMHiggsAnalysis/data/wgt_yymass_shape.root");
    TH1F* h_gamgamweights = NULL;
    if(isMC_gamgam)  h_gamgamweights = (TH1F *) GamGamWeights_File->Get("wgt_yymass_shape");



    //####################################################################################################################
    //###########################################           EVENT LOOP         ###########################################
    //####################################################################################################################


    float lumiXsecWeight(1.0);
    // loop on all the events
    int treeStep = (evEnd-evStart)/50;
    if(treeStep==0)treeStep=1;
    printf("Progressing Bar     :0%%       20%%       40%%       60%%       80%%       100%%\n");
    printf("Scanning the ntuple :");

    for( int iev=evStart; iev<evEnd; iev++) {
        if((iev-evStart)%treeStep==0) {
            printf("#");
            fflush(stdout);
        }



        //prepare the tag's vectors for histo filling
        std::vector<TString> tags(1,"all");

        //##############################################   EVENT LOOP STARTS   ##############################################
        //load the event content from tree
        summaryHandler_.getEntry(iev);
        DataEvtSummary_t &ev=summaryHandler_.getEvent();

        //if(!isMC && ev.RunNumber>301973) continue;
        //all the weights, include XS*BR*EFF, genWeights * PU * PVz
        float weight = 1.0;

        if(isMC) weight *= ev.totWeight;
        if(isMC) lumiXsecWeight = ev.lumiXsecWeight;
        if(isMC) weight *= 15.365; //13.783; //lumi weight, from h013, ntuple is done with 1fb-1

        //remove the pileup weights
        //float pileupWeight = ev.pileupWeight;
        //if(isMC && pileupWeight!=0) weight /= pileupWeight;

        //cout << "ev.totWeight: " << ev.totWeight << " ev.lumiXsecWeight: " << ev.lumiXsecWeight << endl;
        //cout << "lumiXsecWeight: " << lumiXsecWeight << endl;


        //special for gamgam3jet sample
        if(isMC_gamgam) {
            weight = 1.0;
            weight *= ev.evtWeight;

	    //            weight *= 122956*0.8;
            weight *= 136806*0.8;
            weight /= 3.1000714e+07;
        }

        if(isMC_gamjet) {
            weight = 1.0;
            weight *= ev.evtWeight;

            weight *= 136806*0.2;
            weight /= 3.0740822e+07;
        }

        // remove the overlap of Vgam and Vgamgam
        if(isMC_vgam) {
            vector<int> truth_type, truth_origin;
            truth_type.clear();
            truth_origin.clear();
            for(int igen=0; igen<ev.nmcparticles; igen++) {
                if(ev.mc_id[igen] != 22) continue;
                truth_type.push_back(ev.mc_type[igen]);
                truth_origin.push_back(ev.mc_origin[igen]);
            }

            if(truth_type.size()>1) {

                bool remove_Vyy = (truth_type[0]>12&&truth_type[0]<17);   //y1 match to true photon
                remove_Vyy &= (truth_type[1]>12&&truth_type[1]<17);        //y1 match to true photon
                remove_Vyy &= (  truth_origin[0]==3 &&  truth_origin[1]==3);   //y1 and y2 have origin single photon
                if(remove_Vyy) continue;
            }
        }


        float gen_sumet = ev.mc_sumet;
        float gen_met(0), gen_yypt(0);
        float gen_metsig(0);

        bool passGencuts_bin1(true);
        bool passGencuts_bin2(true);
        bool passGencuts_bin3(true);
        bool passGencuts_bin4(true);
        bool passGencuts_bin5(true);

        if( isSignal || isMC_H125 ) {

            LorentzVector P4yy(0,0,0,0);
            LorentzVector P4met(0,0,0,0);

            std::vector<LorentzVector> GenPhotons;
            std::vector<LorentzVector> GenElectrons;

            for(int igen=0; igen<ev.nmcparticles; igen++) {
	      if(ev.mc_id[igen] != 22 || fabs(ev.mc_id[igen]) != 11 ) continue;
	      if( ev.mc_id[igen] == 22 ){
                LorentzVector P4(ev.mc_px[igen],ev.mc_py[igen],ev.mc_pz[igen],ev.mc_en[igen]);
                GenPhotons.push_back(P4);
                P4yy += P4;
	      } 
	      if(fabs( ev.mc_id[igen] ) == 11 ){
                LorentzVector P4(ev.mc_px[igen],ev.mc_py[igen],ev.mc_pz[igen],ev.mc_en[igen]);
                GenElectrons.push_back(P4);
	      }
            }

	    if( ev.ntruthDarkMatters  < 2 ){
	      for(int igen=0; igen<ev.nmcparticles; igen++) {
                if(ev.mc_id[igen] != 0) continue; // met
                LorentzVector P4(ev.mc_px[igen],ev.mc_py[igen],ev.mc_pz[igen],ev.mc_en[igen]);
                P4met += P4;
	      }
	    }else{
	      for(int igen=0; igen<ev.ntruthDarkMatters; igen++) {
                LorentzVector P4(ev.DarkMatterTruth_Px[igen]/1000,ev.DarkMatterTruth_Py[igen]/1000,ev.DarkMatterTruth_Pz[igen]/1000,ev.DarkMatterTruth_E[igen]/1000);
                P4met += P4;
	      }
	    }

            gen_met = P4met.pt();
            gen_yypt = P4yy.pt();
            if(gen_sumet!=0) gen_metsig =  gen_met/sqrt(gen_sumet) ;

            mon.fillHisto("diphoton_mass_gen", tags, P4yy.mass(), 1.);
            mon.fillHisto("diphoton_pt_gen", tags, P4yy.pt(), 1.);
            mon.fillHisto("met_gen", tags, P4met.pt(), 1.);
            mon.fillHisto("sumet_gen", tags, gen_sumet, 1.);
            mon.fillHisto("metsig_gen", tags, gen_metsig, 1.);

            bool passGencuts(true);
            passGencuts &= (GenPhotons.size()==2);
            passGencuts &= (P4yy.mass()>124 && P4yy.mass() <126);

	    if( GenPhotons.size() >= 2 ){
	      double leading_ypt = GenPhotons[0].pt();
	      double trailing_ypt = GenPhotons[1].pt();
	      if(trailing_ypt>leading_ypt) {
                double tmp = leading_ypt;
                leading_ypt = trailing_ypt;
                trailing_ypt = tmp;
	      }

	      passGencuts &= P4yy.mass() != 0 ? (leading_ypt/P4yy.mass() > 0.35) : false ;
	      passGencuts &= P4yy.mass() != 0 ? (trailing_ypt/P4yy.mass() > 0.25) : false ;

	      passGencuts_bin1 &= (passGencuts && gen_metsig>7 && P4yy.pt()>90);
	      passGencuts_bin2 &= (passGencuts && gen_metsig>7 && P4yy.pt()<90);
	      passGencuts_bin3 &= (passGencuts && gen_metsig<7 && (P4yy.pt()>25 && gen_metsig>4) );
	      passGencuts_bin4 &= (passGencuts && gen_metsig<7 && P4yy.pt()>15 && (P4yy.pt()<25 || gen_metsig<4) );
	      passGencuts_bin5 &= (passGencuts_bin1 && GenElectrons.size() == 0 );

	      if(passGencuts) mon.fillHisto("Nevent_gen",tags, 1, 1);
	      if(passGencuts_bin1) mon.fillHisto("Nevent_gen",tags, 3, 1);
	      if(passGencuts_bin2) mon.fillHisto("Nevent_gen",tags, 5, 1);
	      if(passGencuts_bin3) mon.fillHisto("Nevent_gen",tags, 7, 1);
	      if(passGencuts_bin4) mon.fillHisto("Nevent_gen",tags, 9, 1);
	      if(passGencuts_bin5) mon.fillHisto("Nevent_gen",tags, 11, 1);
	      // if(passGencuts_bin1) mon.fillHisto("Nevent_gen",tags, 1, 1);
	      // if(passGencuts_bin2) mon.fillHisto("Nevent_gen",tags, 3, 1);
	      // if(passGencuts_bin3) mon.fillHisto("Nevent_gen",tags, 5, 1);
	      // if(passGencuts_bin4) mon.fillHisto("Nevent_gen",tags, 7, 1);
	      // if(passGencuts_bin5) mon.fillHisto("Nevent_gen",tags, 11, 1);
	    }

        }


        // add PhysicsEvent_t class, get all tree to physics objects
        PhysicsEvent_t phys=getPhysicsEventFrom(ev);




        // looping photons
        int nGoodPhotons(0);
        std::vector<std::pair<int,LorentzVector> > GoodPhotons;
        for(size_t ipho=0; ipho<phys.photons.size(); ipho++) {
            LorentzVector pho=phys.photons[ipho];
            int phoid = phys.photons[ipho].id;
            std::pair <int,LorentzVector> goodpho;
            goodpho = std::make_pair(phoid,pho);
            GoodPhotons.push_back(goodpho);
            nGoodPhotons++;
        }

        if(nGoodPhotons<2) continue; // 2 tight photons



        LorentzVector pho1 = GoodPhotons[0].second;
        LorentzVector pho2 = GoodPhotons[1].second;
        LorentzVector diphoton(pho1+pho2);

        bool passLeadingPhoton (pho1.pt()/diphoton.mass() > 0.35);
        bool passTrailingPhoton (pho2.pt()/diphoton.mass() > 0.25);
        bool passmassWindow(diphoton.mass()>105 && diphoton.mass()<160);

        LorentzVector met = phys.met; // diphoton vertex based
        //LorentzVector met = phys.met_hv; // hardest vertex based
        double dphiGamGamMET=fabs(deltaPhi(diphoton.phi(),met.phi()));
        double balanceDif = fabs(1-met.pt()/diphoton.pt());


        //
        // Electrons
        //
        std::vector<std::pair<int,LorentzVector> > GoodElectrons;
        for(size_t iele=0; iele<phys.electrons.size(); iele++) {
            LorentzVector ele=phys.electrons[iele];
            int eleid = phys.electrons[iele].id;
            std::pair <int,LorentzVector> goodlep;
            goodlep = std::make_pair(eleid,ele);
            GoodElectrons.push_back(goodlep);
        }


        //
        // Muons
        //

        std::vector<std::pair<int,LorentzVector> > GoodMuons;
        for(size_t imu=0; imu<phys.muons.size(); imu++) {
            LorentzVector mu=phys.muons[imu];
            int muid = phys.muons[imu].id;
            std::pair <int,LorentzVector> goodlep;
            goodlep = std::make_pair(muid,mu);
            GoodMuons.push_back(goodlep);
        }



        //
        // Jets
        //
        LorentzVector leading_jet(0,0,0,0);
        LorentzVector subleading_jet(0,0,0,0);

        PhysicsObjectJetCollection GoodJets;
        for(size_t ijet=0; ijet<phys.jets.size(); ijet++) {
            GoodJets.push_back(phys.jets[ijet]);
            if(phys.jets[ijet].pt()>leading_jet.pt()) {
                LorentzVector tmpjet = leading_jet;
                leading_jet = phys.jets[ijet];
                subleading_jet = tmpjet;
            } else if(phys.jets[ijet].pt()>subleading_jet.pt()) {
                subleading_jet = phys.jets[ijet];
            }
        }

	if(GoodJets.size()==3) mon.fillHisto("diphoton_mass_3jetctrl", tags, diphoton.mass(), weight);



        LorentzVector hardsum(0,0,0,0);
        for(size_t ijet=0; ijet<GoodJets.size(); ijet++) {
            hardsum += GoodJets[ijet];
        }
        for(size_t ipho=0; ipho<GoodPhotons.size(); ipho++) {
            hardsum += GoodPhotons[ipho].second;
        }

        /*
            	fprintf(outTxtFile_final,"   %.5f %.5f %.5f %.5f\n",pho1.pt(),pho1.eta(),pho1.phi(),pho1.E());
            	fprintf(outTxtFile_final,"   %.5f %.5f %.5f %.5f\n",pho2.pt(),pho2.eta(),pho2.phi(),pho2.E());
            	for(size_t ipho=2; ipho<GoodPhotons.size(); ipho++) {
                        LorentzVector pho = GoodPhotons[ipho].second;
            	    fprintf(outTxtFile_final,"   %.5f %.5f %.5f %.5f\n",pho.pt(),pho.eta(),pho.phi(),pho.E());
            	}
        */


        if(!passLeadingPhoton || !passTrailingPhoton || !passmassWindow) continue;

        if(!isMC) fprintf(outTxtFile_final,"%d %d %.5f %.5f %.5f\n",ev.RunNumber, ev.EventNumber, diphoton.mass(), diphoton.pt(), met.pt());



        float sum_et = ev.sumet/1000.;
        float met_sig = met.pt()/sqrt(sum_et);


        if(isMC_gamjet) weight *= getSFfrom1DHist(met_sig, h_gamjetweights);
        //if(isMC_gamgam) weight *= getSFfrom1DHist(diphoton.mass(), h_gamgamweights);


        mon.fillHisto("pu_sel",   tags, ev.mu,      1.0);
        mon.fillHisto("puwgt_sel",tags, ev.mu,      weight);
        //# of PV
        mon.fillHisto("nvtx_sel",   tags, phys.nvtx,      1.0);
        mon.fillHisto("nvtxwgt_sel",tags, phys.nvtx,      weight);


        mon.fillHisto("nphotons_sel",tags, nGoodPhotons, weight);
        mon.fillHisto("diphoton_pt_sel", tags, diphoton.pt(), weight);
        mon.fillHisto("diphoton_pt_rebin_sel", tags, diphoton.pt(), weight, true);
        mon.fillHisto("met_sel", tags, met.pt(), weight);
        mon.fillHisto("sumet_sel", tags, sum_et, weight);
        mon.fillHisto("metsig_sel", tags, met_sig, weight);
        mon.fillHisto("metsig2_sel", tags, met_sig, weight);
        if(met_sig>6) mon.fillHisto("metsig_selbeg6", tags, 0, weight);
        if(met_sig>7) mon.fillHisto("metsig_selbeg7", tags, 0, weight);

        if(phys.nvtx>0 && phys.nvtx<10) {
            mon.fillHisto("met_mu0_10_sel", tags, met.pt(), weight);
            mon.fillHisto("sumet_mu0_10_sel", tags, sum_et, weight);
            mon.fillHisto("metsig_mu0_10_sel", tags, met_sig, weight);
        }
        if(phys.nvtx>10 && phys.nvtx<20) {
            mon.fillHisto("met_mu10_20_sel", tags, met.pt(), weight);
            mon.fillHisto("sumet_mu10_20_sel", tags, sum_et, weight);
            mon.fillHisto("metsig_mu10_20_sel", tags, met_sig, weight);
        }
        if(phys.nvtx>20 && phys.nvtx<30) {
            mon.fillHisto("met_mu20_30_sel", tags, met.pt(), weight);
            mon.fillHisto("sumet_mu20_30_sel", tags, sum_et, weight);
            mon.fillHisto("metsig_mu20_30_sel", tags, met_sig, weight);
        }





        mon.fillHisto("metsig_rebin_sel", tags, met_sig, weight,true);
        mon.fillHisto("balancedif_sel"           ,tags, balanceDif, weight);
        mon.fillHisto("metphi_sel", tags, met.phi(), weight);
        mon.fillHisto("dphiGamGamMET_sel",tags, dphiGamGamMET, weight);

        mon.fillHisto("leadingpho_pt_sel", tags, pho1.pt(), weight);
        mon.fillHisto("trailingpho_pt_sel", tags, pho2.pt(), weight);
        mon.fillHisto("leadingpho_eta_sel", tags, pho1.eta(), weight);
        mon.fillHisto("trailingpho_eta_sel", tags, pho2.eta(), weight);
        mon.fillHisto("leadingpho_phi_sel", tags, pho1.phi(), weight);
        mon.fillHisto("trailingpho_phi_sel", tags, pho2.phi(), weight);

        mon.fillHisto("nelectrons_sel",tags, GoodElectrons.size(), weight);
        mon.fillHisto("nmuons_sel",tags, GoodMuons.size(), weight);
        mon.fillHisto("njets_sel",tags, GoodJets.size(), weight);
        mon.fillHisto("pthard_sel",tags, hardsum.pt(), weight);

        if(leading_jet.pt()>0) mon.fillHisto("leadingjetpt_sel",tags, leading_jet.pt(), weight, true);
        if(subleading_jet.pt()>0) mon.fillHisto("subleadingjetpt_sel",tags, subleading_jet.pt(), weight, true);




        if(met_sig>1) mon.fillHisto("diphoton_mass_metsigbeg1", tags, diphoton.mass(), weight);
        if(met_sig>2) mon.fillHisto("diphoton_mass_metsigbeg2", tags, diphoton.mass(), weight);
        if(met_sig>3) mon.fillHisto("diphoton_mass_metsigbeg3", tags, diphoton.mass(), weight);
        if(met_sig>4) mon.fillHisto("diphoton_mass_metsigbeg4", tags, diphoton.mass(), weight);
        if(met_sig>5) mon.fillHisto("diphoton_mass_metsigbeg5", tags, diphoton.mass(), weight);
        if(met_sig>6) mon.fillHisto("diphoton_mass_metsigbeg6", tags, diphoton.mass(), weight);

        if(met_sig>1) mon.fillHisto("diphoton_pt_metsigbeg1", tags, diphoton.pt(), weight);
        if(met_sig>2) mon.fillHisto("diphoton_pt_metsigbeg2", tags, diphoton.pt(), weight);
        if(met_sig>3) mon.fillHisto("diphoton_pt_metsigbeg3", tags, diphoton.pt(), weight);
        if(met_sig>4) mon.fillHisto("diphoton_pt_metsigbeg4", tags, diphoton.pt(), weight);
        if(met_sig>5) mon.fillHisto("diphoton_pt_metsigbeg5", tags, diphoton.pt(), weight);
        if(met_sig>6) mon.fillHisto("diphoton_pt_metsigbeg6", tags, diphoton.pt(), weight);

	mon.fillHisto("pv_dist", tags, phys.pv_hard_z-phys.pv_diphot_z, weight);
	mon.fillHisto("pv_dist_vs_met_sig", tags, phys.pv_hard_z-phys.pv_diphot_z, met_sig, weight);

        if(met_sig>1) mon.fillHisto("pv_dist_metsigbeg1", tags, phys.pv_hard_z-phys.pv_diphot_z, weight);
        if(met_sig>2) mon.fillHisto("pv_dist_metsigbeg2", tags, phys.pv_hard_z-phys.pv_diphot_z, weight);
        if(met_sig>3) mon.fillHisto("pv_dist_metsigbeg3", tags, phys.pv_hard_z-phys.pv_diphot_z, weight);
        if(met_sig>4) mon.fillHisto("pv_dist_metsigbeg4", tags, phys.pv_hard_z-phys.pv_diphot_z, weight);
        if(met_sig>5) mon.fillHisto("pv_dist_metsigbeg5", tags, phys.pv_hard_z-phys.pv_diphot_z, weight);
        if(met_sig>6) mon.fillHisto("pv_dist_metsigbeg6", tags, phys.pv_hard_z-phys.pv_diphot_z, weight);






        mon.fillHisto("Nevent",   tags, 0, weight);
        mon.fillHisto("diphoton_mass_sel", tags, diphoton.mass(), weight);

//        if(diphoton.mass()>122 && diphoton.mass()<128) {

        for(size_t metsig=0; metsig<optim_Cuts_METSig.size(); metsig++) {
            for(size_t dphi=0; dphi<optim_Cuts_DphiyyMET.size(); dphi++) {

                if(met.pt()>50 && diphoton.pt()>50 && dphiGamGamMET>optim_Cuts_DphiyyMET[dphi] && met_sig > optim_Cuts_METSig[metsig]) {
                    mon.fillHisto("yields_dphi_vs_metsig",tags, metsig, dphi, weight);
                }
            }
        }

        for(size_t imet=0; imet<optim_Cuts_MET.size(); imet++) {
            for(size_t dphi=0; dphi<optim_Cuts_DphiyyMET.size(); dphi++) {

                if(met.pt()>optim_Cuts_MET[imet] && diphoton.pt()>optim_Cuts_MET[imet] && dphiGamGamMET>optim_Cuts_DphiyyMET[dphi]) {
                    mon.fillHisto("yields_dphi_vs_met",tags, imet, dphi, weight);
                }
            }
        }



        for(size_t ipt=0; ipt<optim_Cuts_yypt.size(); ipt++) {
            for(size_t metsig=0; metsig<optim_Cuts_METSig.size(); metsig++) {

                if( diphoton.pt()>optim_Cuts_yypt[ipt] && met_sig > optim_Cuts_METSig[metsig]) {
                    mon.fillHisto("yields_metsig_vs_yypt",tags, ipt, metsig, weight);
                }
            }
        }



//        }



        mon.fillHisto("eventflow",tags, 15+1, weight);
        mon.fillHisto("eventflow_raw",tags, 15+1, 1);
        mon.fillHisto("yields_diphoton_mass_bin0", tags, 0, weight);
        if(met_sig>7) {
            if(diphoton.pt()>90) {
	      if( GoodElectrons.size() == 0 ){
                if(isSignal && passGencuts_bin5) mon.fillHisto("Nevent_gen",tags, 2, 1);
                mon.fillHisto("diphoton_mass_bin5", tags, diphoton.mass(), weight);
                mon.fillHisto("yields_diphoton_mass_bin5", tags, 0, weight);
                mon.fillHisto("yields_finalsel",tags, 16, weight);
		mon.fillHisto("yields_finalsel5cat",tags, 4, weight);
                mon.fillHisto("eventflow",tags, 19+2, weight);
                mon.fillHisto("eventflow_raw",tags, 19+2, 1);



                mon.fillHisto("metsig_bin5", tags, met_sig, weight);
                mon.fillHisto("balancedif_bin5", tags, balanceDif, weight);
                mon.fillHisto("dphiGamGamMET_bin5",tags, dphiGamGamMET, weight);
	      }
	      if(isSignal && passGencuts_bin1) mon.fillHisto("Nevent_gen",tags, 2, 1);
	      mon.fillHisto("diphoton_mass_bin1", tags, diphoton.mass(), weight);
	      mon.fillHisto("yields_diphoton_mass_bin1", tags, 0, weight);
	      mon.fillHisto("yields_finalsel",tags, 0, weight);
	      mon.fillHisto("yields_finalsel5cat",tags, 0, weight);
	      mon.fillHisto("eventflow",tags, 15+2, weight);
	      mon.fillHisto("eventflow_raw",tags, 15+2, 1);



	      mon.fillHisto("metsig_bin1", tags, met_sig, weight);
	      mon.fillHisto("balancedif_bin1", tags, balanceDif, weight);
	      mon.fillHisto("dphiGamGamMET_bin1",tags, dphiGamGamMET, weight);
            } else {
                if(isSignal && passGencuts_bin2) mon.fillHisto("Nevent_gen",tags, 4, 1);
                mon.fillHisto("diphoton_mass_bin2", tags, diphoton.mass(), weight);
                mon.fillHisto("yields_diphoton_mass_bin2", tags, 0, weight);
                mon.fillHisto("yields_finalsel",tags, 1, weight);
                mon.fillHisto("yields_finalsel5cat",tags, 1, weight);
                mon.fillHisto("eventflow",tags, 16+2, weight);
                mon.fillHisto("eventflow_raw",tags, 16+2, 1);



            }
        } else if(met_sig>4 && diphoton.pt()>25 ) {
            if(isSignal && passGencuts_bin3) mon.fillHisto("Nevent_gen",tags, 6, 1);
            mon.fillHisto("diphoton_mass_bin3", tags, diphoton.mass(), weight);
            mon.fillHisto("yields_diphoton_mass_bin3", tags, 0, weight);
            mon.fillHisto("yields_finalsel",tags, 2, weight);
            mon.fillHisto("yields_finalsel5cat",tags, 2, weight);
            mon.fillHisto("eventflow",tags, 17+2, weight);
            mon.fillHisto("eventflow_raw",tags, 17+2, 1);



        } else if( diphoton.pt()>15) {
            if(isSignal && passGencuts_bin4) mon.fillHisto("Nevent_gen",tags, 8, 1);
            mon.fillHisto("diphoton_mass_bin4", tags, diphoton.mass(), weight);
            mon.fillHisto("yields_diphoton_mass_bin4", tags, 0, weight);
            mon.fillHisto("yields_finalsel",tags, 3, weight);
            mon.fillHisto("yields_finalsel5cat",tags, 3, weight);
            mon.fillHisto("eventflow",tags, 18+2, weight);
            mon.fillHisto("eventflow_raw",tags, 18+2, 1);


        }

        if(met.pt()>10 && diphoton.pt()>10) mon.fillHisto("yields_finalsel",tags, 4, weight);
        if(met.pt()>20 && diphoton.pt()>20) mon.fillHisto("yields_finalsel",tags, 5, weight);
        if(met.pt()>30 && diphoton.pt()>30) mon.fillHisto("yields_finalsel",tags, 6, weight);
        if(met.pt()>40 && diphoton.pt()>40) mon.fillHisto("yields_finalsel",tags, 7, weight);
        if(met.pt()>50 && diphoton.pt()>50) mon.fillHisto("yields_finalsel",tags, 8, weight);
        if(met.pt()>60 && diphoton.pt()>60) mon.fillHisto("yields_finalsel",tags, 9, weight);
        if(met.pt()>70 && diphoton.pt()>70) mon.fillHisto("yields_finalsel",tags, 10, weight);
        if(met.pt()>80 && diphoton.pt()>80) mon.fillHisto("yields_finalsel",tags, 11, weight);
        if(met.pt()>90 && diphoton.pt()>90) mon.fillHisto("yields_finalsel",tags, 12, weight);
        if(met.pt()>100 && diphoton.pt()>100) mon.fillHisto("yields_finalsel",tags, 13, weight);
        if(met.pt()>110 && diphoton.pt()>110) mon.fillHisto("yields_finalsel",tags, 14, weight);
        if(met.pt()>120 && diphoton.pt()>120) mon.fillHisto("yields_finalsel",tags, 15, weight);




	LorentzVector allpart;

	for( auto photon : GoodPhotons)
	  allpart += photon.second;
	for( auto jet : GoodJets)
	  allpart += jet;
	for( auto electron : GoodElectrons)
	  allpart += electron.second;
	for( auto muon : GoodMuons)
	  allpart += muon.second;

	mon.fillHisto("allpart_phi", tags, fabs(deltaPhi(allpart.phi(),met.phi())), weight);
	mon.fillHisto("allpart_phi_vs_met_sig", tags, fabs(deltaPhi(allpart.phi(),met.phi())),met_sig, weight);



	if(met_sig>1) mon.fillHisto("allpart_phi_metsigbeg1", tags, fabs(deltaPhi(allpart.phi(),met.phi())), weight);
        if(met_sig>2) mon.fillHisto("allpart_phi_metsigbeg2", tags, fabs(deltaPhi(allpart.phi(),met.phi())), weight);
        if(met_sig>3) mon.fillHisto("allpart_phi_metsigbeg3", tags, fabs(deltaPhi(allpart.phi(),met.phi())), weight);
        if(met_sig>4) mon.fillHisto("allpart_phi_metsigbeg4", tags, fabs(deltaPhi(allpart.phi(),met.phi())), weight);
        if(met_sig>5) mon.fillHisto("allpart_phi_metsigbeg5", tags, fabs(deltaPhi(allpart.phi(),met.phi())), weight);
        if(met_sig>6) mon.fillHisto("allpart_phi_metsigbeg6", tags, fabs(deltaPhi(allpart.phi(),met.phi())), weight);
	




	for( unsigned int i = 0 ; i < 5 ; ++i ){

	  std::vector<double> var;
	  var.push_back(met_sig);
	  var.push_back(diphoton.pt());
	  var.push_back(0);
	  var.push_back(GoodElectrons.size());

	  if(Categorisation(var,i)){
	    mon.fillHisto("allpart_phi_cat"+TString::Format("%d",i), tags, fabs(deltaPhi(allpart.phi(),met.phi())), weight);
	    mon.fillHisto("pv_dist_cat"+TString::Format("%d",i), tags, phys.pv_hard_z-phys.pv_diphot_z, weight);
	  }

	}


        if( ifsaveEvents /*&& (diphoton.mass()<120 || diphoton.mass()>130)*/ ) {

            //if(!isMC && diphoton.mass()>120 && diphoton.mass()<130) continue;

            if(met_sig>7) {
                if(diphoton.pt()>90) {

                    mgg_bin1 = diphoton.mass();
                    weight_bin1_t = weight;
                    myEvents_bin1->Fill();

		    mycat_idx=0;

                } else {

                    mgg_bin2 = diphoton.mass();
                    weight_bin2_t = weight;
                    myEvents_bin2->Fill();

		    mycat_idx=1;

                }
            } else if(met_sig>4 && diphoton.pt()>25) {

                mgg_bin3 = diphoton.mass();
                weight_bin3_t = weight;
                myEvents_bin3->Fill();

		mycat_idx=2;

            } else if(diphoton.pt()>15) {

                mgg_bin4 = diphoton.mass();
                weight_bin4_t = weight;
                myEvents_bin4->Fill();

		mycat_idx=3;

            }else{
		mycat_idx=-1;
	    }



            //if(met.pt()>10 && diphoton.pt()>10 ) {
	    //test
	    if(diphoton.pt()>10 ) {
                mgg_bin5 = diphoton.mass();
                weight_bin5_t = weight;
                myEvents_bin5->Fill();
            }
            if(diphoton.pt()>20 ) {
                mgg_bin6 = diphoton.mass();
                weight_bin6_t = weight;
                myEvents_bin6->Fill();
            }
            if(diphoton.pt()>30 ) {
                mgg_bin7 = diphoton.mass();
                weight_bin7_t = weight;
                myEvents_bin7->Fill();
            }
            if(diphoton.pt()>40 ) {
                mgg_bin8 = diphoton.mass();
                weight_bin8_t = weight;
                myEvents_bin8->Fill();
            }
            if(diphoton.pt()>50 ) {
                mgg_bin9 = diphoton.mass();
                weight_bin9_t = weight;
                myEvents_bin9->Fill();
            }
            if(diphoton.pt()>60 ) {
                mgg_bin10 = diphoton.mass();
                weight_bin10_t = weight;
                myEvents_bin10->Fill();
            }

            if(diphoton.pt()>70 ) {
                mgg_bin11 = diphoton.mass();
                weight_bin11_t = weight;
                myEvents_bin11->Fill();
            }
            if(diphoton.pt()>80 ) {
                mgg_bin12 = diphoton.mass();
                weight_bin12_t = weight;
                myEvents_bin12->Fill();
            }
            if(diphoton.pt()>90 ) {
                mgg_bin13 = diphoton.mass();
                weight_bin13_t = weight;
                myEvents_bin13->Fill();
            }
            if(diphoton.pt()>100 ) {
                mgg_bin14 = diphoton.mass();
                weight_bin14_t = weight;
                myEvents_bin14->Fill();
            }
            if(diphoton.pt()>110 ) {
                mgg_bin15 = diphoton.mass();
                weight_bin15_t = weight;
                myEvents_bin15->Fill();
            }
            if(diphoton.pt()>120 ) {
                mgg_bin16 = diphoton.mass();
                weight_bin16_t = weight;
                myEvents_bin16->Fill();
            }


            mgg_allbins = diphoton.mass();
            weight_allbins_t = weight;
            myEvents_allbins->Fill();


        }





    } // loop on all events END


    if(isMC) {

        for(int bin=1; bin < Heventflow->GetNbinsX()+1; bin++) {
            if(bin<16) {
                float bincontent = Heventflow->GetBinContent(bin);
                Heventflow -> SetBinContent(bin, bincontent * lumiXsecWeight);
            }
        }
    }

    //##############################################
    //########     SAVING HISTO TO FILE     ########
    //##############################################

    printf("\n");
    //save control plots to file
    outUrl += "/";
    outUrl += outFileUrl + ".root";
    printf("Results saved in %s\n", outUrl.Data());

    //save all to the file
    TFile *ofile=TFile::Open(outUrl, "recreate");
    mon.Write();
    ofile->Close();

    // saving hgg tree to output file
    if(ifsaveEvents) {
        otreefile->cd();
        myEvents_bin1->Write();
        myEvents_bin2->Write();
        myEvents_bin3->Write();
        myEvents_bin4->Write();
        //model indepedent
        myEvents_bin5->Write();
        myEvents_bin6->Write();
        myEvents_bin7->Write();
        myEvents_bin8->Write();
        myEvents_bin9->Write();
        myEvents_bin10->Write();
        myEvents_bin11->Write();
        myEvents_bin12->Write();
        myEvents_bin13->Write();
        myEvents_bin14->Write();
        myEvents_bin15->Write();
        myEvents_bin16->Write();

        myEvents_allbins->Write();
        printf("Tree saved in %s\n", outTreeUrl.Data());
    }
    otreefile->Close();
    if(!ifsaveEvents) gSystem->Exec("rm " + outTreeUrl);

    printf("\n");
    inputFile->Close();

    GamjetsWeights_File->Close();
    GamGamWeights_File->Close();

    if(outTxtFile_final)fclose(outTxtFile_final);

    return 0;
}
