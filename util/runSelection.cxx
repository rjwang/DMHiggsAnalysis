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

#include "DMHiggsAnalysis/SmartSelectionMonitor.h"
#include "DMHiggsAnalysis/deltaR.h"



using namespace std;

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
    bool isSignal   = runProcess.getBool("issignal");
    bool isMC_H125  = isMC && (inputUrl.Contains("H125"));
    bool isMC_gamgam = isMC && ( inputUrl.Contains("_gamgam_") || inputUrl.Contains("Pythia8_1Hard1BremDP40") );
    bool isMC_gamjet = isMC && (inputUrl.Contains("_gamjet_"));


    TFile *inputFile = TFile::Open(inputUrl);
    printf("Looping on %s\n",inputUrl.Data());
    if(inputFile==0) return -1;
    if(inputFile->IsZombie()) return -1;
    if( !summaryHandler_.attachToTree( (TTree *)inputFile->Get(dirname) ) ) {
        inputFile->Close();
        return -1;
    }

    const int totalEntries= summaryHandler_.getEntries();
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

    TString outTxtUrl_final= outUrl + "/" + outFileUrl + "_FinalList.txt";
    FILE* outTxtFile_final = NULL;
    outTxtFile_final = fopen(outTxtUrl_final.Data(), "w");
    printf("TextFile URL = %s\n",outTxtUrl_final.Data());
    fprintf(outTxtFile_final,"run event diphotonmass categorybin evtWeight lumiWeight pthard diphotonpt met\n");


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
    myEvents_allbins->Branch("mgg",&mgg_allbins,"mgg/D");
    myEvents_allbins->Branch("weight_t", &weight_allbins_t,"weight_t/D");










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

    mon.addHistogram( new TH1F( "diphoton_mass_sel",    ";#it{m}_{#gamma#gamma} [GeV];Events", 62,105,160) );
    mon.addHistogram( new TH1F( "diphoton_pt_sel",    ";#it{p}_{T}^{#gamma#gamma} [GeV];Events", 100,0,600) );
    mon.addHistogram( new TH1F( "met_sel",    	      ";#it{E}_{T}^{miss} [GeV];Events", 100,0,600) );
    mon.addHistogram( new TH1F( "sumet_sel",            ";sumE_{T} [GeV];Events", 100,50,1050) );
    mon.addHistogram( new TH1F( "metsig_sel",         ";E_{T}^{miss} Significance [GeV];Events", 100,0,10) );
    mon.addHistogram( new TH1F( "balancedif_sel",     ";|E_{T}^{miss}-#it{p}_{T}^{#gamma#gamma}|/#it{p}_{T}^{#gamma#gamma};Events", 5,0,1.0) );
    mon.addHistogram( new TH1F( "metphi_sel",         ";#phi(#it{E}_{T}^{miss}) [rad];Events", 50,-TMath::Pi(),TMath::Pi()) );
    mon.addHistogram( new TH1F( "pthard_sel",         ";#it{p}_{T}^{hard} [GeV];Events", 100,0,600) );
    mon.addHistogram( new TH1F( "dphiGamGamMET_sel",   ";#Delta#it{#phi}(#it{p}_{T}^{#gamma#gamma},E_{T}^{miss}) [rad];Events", 50,0,TMath::Pi()) );

    mon.addHistogram( new TH1F( "leadingpho_pt_sel",    ";Leading #it{p}_{T}^{#gamma} [GeV];Events", 100,0,600) );
    mon.addHistogram( new TH1F( "trailingpho_pt_sel",    ";Trailing #it{p}_{T}^{#gamma} [GeV];Events", 100,0,600) );
    mon.addHistogram( new TH1F( "leadingpho_eta_sel",    ";Leading #it{#eta}^{#gamma};Events", 50,-2.4,2.4) );
    mon.addHistogram( new TH1F( "trailingpho_eta_sel",    ";Trailing #it{#eta}^{#gamma};Events", 50,-2.4,2.4) );
    mon.addHistogram( new TH1F( "leadingpho_phi_sel",    ";Leading #it{#phi}^{#gamma} [rad];Events", 50,-1.*TMath::Pi(),TMath::Pi()) );
    mon.addHistogram( new TH1F( "trailingpho_phi_sel",    ";Trailing #it{#phi}^{#gamma} [rad];Events", 50,-1.*TMath::Pi(),TMath::Pi()) );

    mon.addHistogram( new TH1F( "dphiGamGamMET_bin1",   ";#Delta#it{#phi}(#it{p}_{T}^{#gamma#gamma},E_{T}^{miss}) [rad];Events", 50,0,TMath::Pi()) );
    mon.addHistogram( new TH1F( "metsig_bin1",         ";E_{T}^{miss} Significance [GeV];Events", 100,0,10) );
    mon.addHistogram( new TH1F( "balancedif_bin1",     ";|E_{T}^{miss}-#it{p}_{T}^{#gamma#gamma}|/#it{p}_{T}^{#gamma#gamma};Events", 5,0,1.0) );


    //optimization
    std::vector<double> optim_Cuts_METSig;
    std::vector<double> optim_Cuts_MET;
    std::vector<double> optim_Cuts_DphiyyMET;


    for(double metsig=0; metsig<=10; metsig+=0.5) {
        optim_Cuts_METSig     .push_back(metsig);
    }
    for(double met=0; met<=200; met+=10) {
        optim_Cuts_MET     .push_back(met);
    }
    for(double dphi=0; dphi<=3.1; dphi+=0.31) {
        optim_Cuts_DphiyyMET  .push_back(dphi);
    }

    size_t nOptims_METSig = optim_Cuts_METSig.size();
    size_t nOptims_MET = optim_Cuts_MET.size();
    size_t nOptims_DphiyyMET = optim_Cuts_DphiyyMET.size();

    TH2F *h2 = (TH2F*) mon.addHistogram( new TH2F( "yields_dphi_vs_metsig",";E_{T}^{miss} Significance; #Delta#it{#phi}(#it{p}_{T}^{#gamma#gamma},E_{T}^{miss}) [rad];Events", nOptims_METSig,0,nOptims_METSig, nOptims_DphiyyMET,0,nOptims_DphiyyMET) );

    TH2F *h3 = (TH2F*) mon.addHistogram( new TH2F( "yields_dphi_vs_met",";E_{T}^{miss} [GeV]; #Delta#it{#phi}(#it{p}_{T}^{#gamma#gamma},E_{T}^{miss}) [rad];Events", nOptims_MET,0,nOptims_MET, nOptims_DphiyyMET,0,nOptims_DphiyyMET) );


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

    h1 = (TH1F*) mon.addHistogram( new TH1F( "yields_finalsel",";;Events", 16,0,16) );
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














    //MC normalization (to 1/pb)
    double sumInitialEvents=1.0;
    double skim_eff=1.0;
    if(isMC) {
        TH1F* H_noDalitz_weighted = (TH1F *) inputFile->Get("CutFlow_noDalitz_weighted");
        TH1F* H_weighted 	  = (TH1F *) inputFile->Get("CutFlow_weighted");
        sumInitialEvents = H_noDalitz_weighted->GetBinContent(3);

        // Hard-coding to bin number 1,2
        double NxAOD      = H_weighted->GetBinContent(1);
        double NDxAOD     = H_weighted->GetBinContent(2);
        skim_eff = NDxAOD / NxAOD;

        printf("sumInitialEvents = %f, DxAOD skimming efficiency = %f\n",sumInitialEvents, skim_eff);


        for(int bin=1; bin < Heventflow->GetNbinsX()+1; bin++) {
            if(bin<16) {
                Heventflow -> SetBinContent(bin, H_noDalitz_weighted->GetBinContent(bin));
                Heventflow -> GetXaxis()->SetBinLabel(bin,H_noDalitz_weighted->GetXaxis()->GetBinLabel(bin));
            }
        }

        Heventflow->GetXaxis()->SetBinLabel(16+1,"m_{#gamma#gamma}");
        Heventflow->GetXaxis()->SetBinLabel(16+2,"High #it{E}_{T}^{miss}, high #it{p}_{T}^{#gamma#gamma}");
        Heventflow->GetXaxis()->SetBinLabel(17+2,"High #it{E}_{T}^{miss}, low #it{p}_{T}^{#gamma#gamma}");
        Heventflow->GetXaxis()->SetBinLabel(18+2,"Intermediate #it{E}_{T}^{miss}");
        Heventflow->GetXaxis()->SetBinLabel(19+2,"Rest category");





    }

    Hcutflow->SetBinContent(1,sumInitialEvents/skim_eff);



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

        //all the weights, include XS*BR*EFF, genWeights * PU * PVz
        float weight = 1.0;
        if(isMC) weight *= ev.totWeight;
        if(isMC) lumiXsecWeight = ev.lumiXsecWeight;
        // add PhysicsEvent_t class, get all tree to physics objects
        PhysicsEvent_t phys=getPhysicsEventFrom(ev);





        // looping photons
        int nGoodPhotons(0);
        std::vector<std::pair<int,LorentzVector> > GoodPhotons;
        for(size_t ipho=0; ipho<phys.photons.size(); ipho++) {
            LorentzVector pho=phys.photons[ipho];
            int phoid = phys.photons[ipho].id;
            /*
                        if(pho.pt()<25) continue;
                        if(fabs(pho.eta())>2.37) continue;
                        if(fabs(pho.eta())>1.37 && fabs(pho.eta())<1.52) continue;

                        //Iso, ID
                        //bool hasTightIdandIso(true);
                        //hasTightIdandIso &= phys.photons[ipho].isTightID;
                        //hasTightIdandIso &= phys.photons[ipho].isTightIso;
                        //if(!hasTightIdandIso) continue;

            */
            std::pair <int,LorentzVector> goodpho;
            goodpho = std::make_pair(phoid,pho);
            GoodPhotons.push_back(goodpho);
            nGoodPhotons++;
        }

        if(nGoodPhotons<2) continue; // 2 tight photons



        LorentzVector pho1 = GoodPhotons[0].second;
        LorentzVector pho2 = GoodPhotons[1].second;
        /*
                //find the leading and trailing photon
                if(pho1.pt()<pho2.pt()) {
                    LorentzVector pho_ = pho1;
                    pho1 = pho2;
                    pho2 = pho_;
                }

                for(size_t ipho=2; ipho<GoodPhotons.size(); ipho++) {
                    LorentzVector pho_ = GoodPhotons[ipho].second;
                    if( pho_.pt() > pho1.pt() ) {
                        pho2 = pho1;
                        pho1 = pho_;
                    }
                    if( pho_.pt() > pho2.pt() && pho_.pt() < pho1.pt() ) {
                        pho2 = pho_;
                    }
                }
        */

        LorentzVector diphoton(pho1+pho2);

        bool passLeadingPhoton (pho1.pt()/diphoton.mass() > 0.35);
        bool passTrailingPhoton (pho2.pt()/diphoton.mass() > 0.25);
        bool passmassWindow(diphoton.mass()>105 && diphoton.mass()<160);


        LorentzVector met = phys.met;
        double dphiGamGamMET=fabs(deltaPhi(diphoton.phi(),met.phi()));
        double balanceDif = fabs(1-met.pt()/diphoton.pt());


        //
        // Electrons
        //
        std::vector<std::pair<int,LorentzVector> > GoodElectrons;
        for(size_t iele=0; iele<phys.electrons.size(); iele++) {
            LorentzVector ele=phys.electrons[iele];
            int eleid = phys.electrons[iele].id;
            /*
                        if(ele.pt()<10) continue;
                        if(fabs(ele.eta())>2.47) continue;
                        if(fabs(ele.eta())>1.37 && fabs(ele.eta())<1.52) continue;

                        if(!phys.electrons[iele].isMediumID) continue;
                        if(!phys.electrons[iele].isLooseIso) continue;

                        bool overlapremoval(false);
                        for(size_t ipho=0; ipho<GoodPhotons.size(); ipho++) {
                            if( reco::deltaR(ele,GoodPhotons[ipho].second) < 0.4) {
                                overlapremoval |=true;
                                break;
                            }
                        }
                        if(overlapremoval) continue;
            */
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
            /*
                        if(mu.pt()<10) continue;
                        if(fabs(mu.eta())>2.7) continue;

                        if(!phys.muons[imu].isIsoGradientLoose) continue;

                        bool overlapremoval(false);
                        for(size_t ipho=0; ipho<GoodPhotons.size(); ipho++) {
                            if( reco::deltaR(mu,GoodPhotons[ipho].second) < 0.4) {
                                overlapremoval |=true;
                                break;
                            }
                        }
                        if(overlapremoval) continue;

                        for(size_t ijet=0; ijet<GoodJets.size(); ijet++) {
                            if( reco::deltaR(mu,GoodJets[ijet]) < 0.4) {
                                overlapremoval |=true;
                                break;
                            }
                        }
                        if(overlapremoval) continue;

            */
            std::pair <int,LorentzVector> goodlep;
            goodlep = std::make_pair(muid,mu);
            GoodMuons.push_back(goodlep);
        }



        //
        // Jets
        //
        PhysicsObjectJetCollection GoodJets;
        for(size_t ijet=0; ijet<phys.jets.size(); ijet++) {

            /*
                        if(phys.jets[ijet].pt()<25) continue;
                        if(fabs(phys.jets[ijet].eta())>4.4) continue;

                        if(phys.jets[ijet].pt()<50 && phys.jets[ijet].pt()>20 && fabs(phys.jets[ijet].eta())<2.4) {
                            if(phys.jets[ijet].JVT > 0.64) continue;
                        }

                        bool overlapremoval(false);
                        for(size_t ipho=0; ipho<GoodPhotons.size(); ipho++) {
                            if( reco::deltaR(phys.jets[ijet],GoodPhotons[ipho].second) < 0.4) {
                                overlapremoval |=true;
                                break;
                            }
                        }
                        if(overlapremoval) continue;

                        for(size_t iele=0; iele<goodElectrons.size(); iele++) {
                            if( reco::deltaR(phys.jets[ijet],goodElectrons[iele].second) < 0.2) {
                                overlapremoval |=true;
                                break;
                            }
                        }
                        if(overlapremoval) continue;
            */

            GoodJets.push_back(phys.jets[ijet]);
        }




        /*
                std::vector<std::pair<int,LorentzVector> > GoodElectrons;
                for(size_t iele=0; iele<goodElectrons.size(); iele++) {
                    bool overlapremoval(false);
                    for(size_t ijet=0; ijet<GoodJets.size(); ijet++) {
                        if( reco::deltaR(goodElectrons[iele].second, GoodJets[ijet]) < 0.4) {
                            overlapremoval |=true;
                            break;
                        }
                    }
                    if(overlapremoval) continue;

                    GoodElectrons.push_back(goodElectrons[iele]);
                }
        */



        LorentzVector hardsum(0,0,0,0);
        for(size_t ijet=0; ijet<GoodJets.size(); ijet++) {
            hardsum += GoodJets[ijet];
        }
        for(size_t ipho=0; ipho<GoodPhotons.size(); ipho++) {
            hardsum += GoodPhotons[ipho].second;
        }



        /*
        	fprintf(outTxtFile_final,"%d %d %.5f %d %.5f %.5f %.5f %.5f %.5f\n",ev.RunNumber, ev.EventNumber, diphoton.mass(), 0, ev.evtWeight, ev.lumiXsecWeight, hardsum.pt(), diphoton.pt(), met.pt());
        	fprintf(outTxtFile_final,"   %.5f %.5f %.5f %.5f\n",pho1.pt(),pho1.eta(),pho1.phi(),pho1.E());
        	fprintf(outTxtFile_final,"   %.5f %.5f %.5f %.5f\n",pho2.pt(),pho2.eta(),pho2.phi(),pho2.E());
        	for(size_t ipho=2; ipho<GoodPhotons.size(); ipho++) {
                    LorentzVector pho = GoodPhotons[ipho].second;
        	    fprintf(outTxtFile_final,"   %.5f %.5f %.5f %.5f\n",pho.pt(),pho.eta(),pho.phi(),pho.E());

        	}

        */


        if(!passLeadingPhoton || !passTrailingPhoton || !passmassWindow) continue;



        float sum_et = ev.sumet/1000.;
        float met_sig = met.pt()/sqrt(sum_et);

        mon.fillHisto("pu_sel",   tags, ev.mu,      1.0);
        mon.fillHisto("puwgt_sel",tags, ev.mu,      weight);
        //# of PV
        mon.fillHisto("nvtx_sel",   tags, phys.nvtx,      1.0);
        mon.fillHisto("nvtxwgt_sel",tags, phys.nvtx,      weight);


        mon.fillHisto("nphotons_sel",tags, nGoodPhotons, weight);
        mon.fillHisto("diphoton_pt_sel", tags, diphoton.pt(), weight);
        mon.fillHisto("met_sel", tags, met.pt(), weight);
        mon.fillHisto("sumet_sel", tags, sum_et, weight);
        mon.fillHisto("metsig_sel", tags, met_sig, weight);
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

        if(isMC) {
            mon.fillHisto("diphoton_mass_sel", tags, diphoton.mass(), weight);
        }
        if(!isMC && (diphoton.mass()<120 || diphoton.mass()>130)) {
            mon.fillHisto("diphoton_mass_sel", tags, diphoton.mass(), weight);
        }

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




//        }



        mon.fillHisto("eventflow",tags, 15+1, weight);
        if(met.pt()>100) {
            if(isMC) {
                if(diphoton.pt()>100) {
                    //if(diphoton.mass()>122 && diphoton.mass()<128) {
                    mon.fillHisto("yields_finalsel",tags, 0, weight);
                    mon.fillHisto("eventflow",tags, 15+2, weight);

                    mon.fillHisto("metsig_bin1", tags, met_sig, weight);
                    mon.fillHisto("balancedif_bin1", tags, balanceDif, weight);
                    mon.fillHisto("dphiGamGamMET_bin1",tags, dphiGamGamMET, weight);

                    fprintf(outTxtFile_final,"%d %d %.5f %d %.5f %.5f %.5f %.5f %.5f\n",ev.RunNumber, ev.EventNumber, diphoton.mass(), 1, ev.evtWeight, ev.lumiXsecWeight, hardsum.pt(), diphoton.pt(), met.pt());
                    //}
                } else {
                    //if(diphoton.mass()>122 && diphoton.mass()<128) {
                    mon.fillHisto("yields_finalsel",tags, 1, weight);
                    mon.fillHisto("eventflow",tags, 16+2, weight);
                    fprintf(outTxtFile_final,"%d %d %.5f %d %.5f %.5f %.5f %.5f %.5f\n",ev.RunNumber, ev.EventNumber, diphoton.mass(), 2, ev.evtWeight, ev.lumiXsecWeight, hardsum.pt(), diphoton.pt(), met.pt());
                    //}
                }
            }
        } else if(met.pt()>50 && hardsum.pt()>40) {
            //if(diphoton.mass()>122 && diphoton.mass()<128) {
            mon.fillHisto("yields_finalsel",tags, 2, weight);
            mon.fillHisto("eventflow",tags, 17+2, weight);
            fprintf(outTxtFile_final,"%d %d %.5f %d %.5f %.5f %.5f %.5f %.5f\n",ev.RunNumber, ev.EventNumber, diphoton.mass(), 3, ev.evtWeight, ev.lumiXsecWeight, hardsum.pt(), diphoton.pt(), met.pt());
            //}
        } else if( diphoton.pt()>15) {
            //if(diphoton.mass()>122 && diphoton.mass()<128) {
            mon.fillHisto("yields_finalsel",tags, 3, weight);
            mon.fillHisto("eventflow",tags, 18+2, weight);
            fprintf(outTxtFile_final,"%d %d %.5f %d %.5f %.5f %.5f %.5f %.5f\n",ev.RunNumber, ev.EventNumber, diphoton.mass(), 4, ev.evtWeight, ev.lumiXsecWeight, hardsum.pt(), diphoton.pt(), met.pt());
            //}
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





















        if( ifsaveEvents /*&& (diphoton.mass()<120 || diphoton.mass()>130)*/ ) {

            if(!isMC && diphoton.mass()>120 && diphoton.mass()<130) continue;

            if(met.pt()>100) {
                if(diphoton.pt()>100) {

                    mgg_bin1 = diphoton.mass();
                    weight_bin1_t = weight;
                    myEvents_bin1->Fill();


                } else {

                    mgg_bin2 = diphoton.mass();
                    weight_bin2_t = weight;
                    myEvents_bin2->Fill();

                }
            } else if(met.pt()>50 && hardsum.pt()>40) {


                mgg_bin3 = diphoton.mass();
                weight_bin3_t = weight;
                myEvents_bin3->Fill();


            } else if(diphoton.pt()>15) {

                mgg_bin4 = diphoton.mass();
                weight_bin4_t = weight;
                myEvents_bin4->Fill();

            }



            if(met.pt()>10 && diphoton.pt()>10 ) {
                mgg_bin5 = diphoton.mass();
                weight_bin5_t = weight;
                myEvents_bin5->Fill();
            }
            if(met.pt()>20 && diphoton.pt()>20 ) {
                mgg_bin6 = diphoton.mass();
                weight_bin6_t = weight;
                myEvents_bin6->Fill();
            }
            if(met.pt()>30 && diphoton.pt()>30 ) {
                mgg_bin7 = diphoton.mass();
                weight_bin7_t = weight;
                myEvents_bin7->Fill();
            }
            if(met.pt()>40 && diphoton.pt()>40 ) {
                mgg_bin8 = diphoton.mass();
                weight_bin8_t = weight;
                myEvents_bin8->Fill();
            }
            if(met.pt()>50 && diphoton.pt()>50 ) {
                mgg_bin9 = diphoton.mass();
                weight_bin9_t = weight;
                myEvents_bin9->Fill();
            }
            if(met.pt()>60 && diphoton.pt()>60 ) {
                mgg_bin10 = diphoton.mass();
                weight_bin10_t = weight;
                myEvents_bin10->Fill();
            }

            if(met.pt()>70 && diphoton.pt()>70 ) {
                mgg_bin11 = diphoton.mass();
                weight_bin11_t = weight;
                myEvents_bin11->Fill();
            }
            if(met.pt()>80 && diphoton.pt()>80 ) {
                mgg_bin12 = diphoton.mass();
                weight_bin12_t = weight;
                myEvents_bin12->Fill();
            }
            if(met.pt()>90 && diphoton.pt()>90 ) {
                mgg_bin13 = diphoton.mass();
                weight_bin13_t = weight;
                myEvents_bin13->Fill();
            }
            if(met.pt()>100 && diphoton.pt()>100 ) {
                mgg_bin14 = diphoton.mass();
                weight_bin14_t = weight;
                myEvents_bin14->Fill();
            }
            if(met.pt()>110 && diphoton.pt()>110 ) {
                mgg_bin15 = diphoton.mass();
                weight_bin15_t = weight;
                myEvents_bin15->Fill();
            }
            if(met.pt()>120 && diphoton.pt()>120 ) {
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


    if(outTxtFile_final)fclose(outTxtFile_final);

    return 0;
}
