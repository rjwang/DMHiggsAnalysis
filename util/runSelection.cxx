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


    TString inputUrl = runProcess.getStr("InputFile","InputFile");

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



    //##################################################################################
    //##########################    INITIATING HISTOGRAMS     ##########################
    //##################################################################################

    SmartSelectionMonitor mon;


    mon.addHistogram( new TH1F( "pu_raw", ";pileup;Events", 50,0,50) );
    mon.addHistogram( new TH1F( "puwgt_raw", ";pileup;Events", 50,0,50) );

    mon.addHistogram( new TH1F( "nvtx_raw",     ";Vertices;Events",50,0,50) );
    mon.addHistogram( new TH1F( "nvtxwgt_raw",  ";Vertices;Events",50,0,50) );


    mon.addHistogram( new TH2F( "pt2mass_pho2vspho1_raw",";#it{p}_{T}^{#gamma1}/#it{m}_{#gamma#gamma};#it{p}_{T}^{#gamma2}/#it{m}_{#gamma#gamma};Events",100,0,1., 100,0,1.) );

    mon.addHistogram( new TH1F( "diphoton_mass_raw",    ";#it{m}_{#gamma#gamma} [GeV];Events", 62,105,160) );
    mon.addHistogram( new TH1F( "diphoton_pt_raw",    ";#it{p}_{T}^{#gamma#gamma} [GeV];Events", 100,0,600) );

    mon.addHistogram( new TH1F( "leadingpho_pt_raw",    ";Leading #it{p}_{T}^{#gamma} [GeV];Events", 100,0,600) );
    mon.addHistogram( new TH1F( "trailingpho_pt_raw",    ";Trailing #it{p}_{T}^{#gamma} [GeV];Events", 100,0,600) );
    mon.addHistogram( new TH1F( "leadingpho_eta_raw",    ";Leading #it{#eta}^{#gamma};Events", 100,-2.4,2.4) );
    mon.addHistogram( new TH1F( "trailingpho_eta_raw",    ";Trailing #it{#eta}^{#gamma};Events", 100,-2.4,2.4) );
    mon.addHistogram( new TH1F( "leadingpho_phi_raw",    ";Leading #it{#phi}^{#gamma};Events", 50,-1.*TMath::Pi(),TMath::Pi()) );
    mon.addHistogram( new TH1F( "trailingpho_phi_raw",    ";Trailing #it{#phi}^{#gamma};Events", 50,-1.*TMath::Pi(),TMath::Pi()) );




    TH1F *h1 = (TH1F*) mon.addHistogram( new TH1F( "nphotons_raw", ";Photon multiplicity;Events", 3,2,5) );
    for(int ibin=1; ibin<=h1->GetXaxis()->GetNbins(); ibin++) {
        TString label("");
        label +="=";
        label += (ibin+1);
        h1->GetXaxis()->SetBinLabel(ibin,label);
    }

    //####################################################################################################################
    //###########################################           EVENT LOOP         ###########################################
    //####################################################################################################################


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

        //bad event filter
        if(!ev.passQualityCuts) continue;


        //weights:  genWeights * PU * PVz
        float weight = 1.0;
        bool isMC = ev.isMC;
        if(isMC) weight *= ev.initWeight;

        // add PhysicsEvent_t class, get all tree to physics objects
        PhysicsEvent_t phys=getPhysicsEventFrom(ev);


        mon.fillHisto("pu_raw",   tags, ev.mu,      1.0);
        mon.fillHisto("puwgt_raw",tags, ev.mu,      weight);

        //# of PV
        mon.fillHisto("nvtx_raw",   tags, phys.nvtx,      1.0);
        mon.fillHisto("nvtxwgt_raw",tags, phys.nvtx,      weight);




        // looping photons
        int nGoodPhotons(0);
        std::vector<std::pair<int,LorentzVector> > goodPhotons;
        for(size_t ipho=0; ipho<phys.photons.size(); ipho++) {
            LorentzVector pho=phys.photons[ipho];
            int phoid = phys.photons[ipho].id;
            if(pho.pt()<25) continue;
            if(fabs(pho.eta())>2.37) continue;
            if(fabs(pho.eta())>1.37 && fabs(pho.eta())<1.52) continue;

            //Iso, ID
            bool hasTightIdandIso(true);
            hasTightIdandIso &= phys.photons[ipho].isTightID;
            hasTightIdandIso &= phys.photons[ipho].isTightIso;
            if(!hasTightIdandIso) continue;

            nGoodPhotons++;
            std::pair <int,LorentzVector> goodpho;
            goodpho = std::make_pair(phoid,pho);
            goodPhotons.push_back(goodpho);
        }

        if(nGoodPhotons<2) continue; // 2 tight photons



        LorentzVector pho1 = goodPhotons[0].second;
        LorentzVector pho2 = goodPhotons[1].second;

	//find the leading and trailing photon
        if(pho1.pt()<pho2.pt()) {
            LorentzVector pho_ = pho1;
            pho1 = pho2;
            pho2 = pho_;
        }

        for(size_t ipho=2; ipho<goodPhotons.size(); ipho++) {
            LorentzVector pho_ = goodPhotons[ipho].second;
            if( pho_.pt() > pho1.pt() ) {
                pho2 = pho1;
                pho1 = pho_;
            }
            if( pho_.pt() > pho2.pt() && pho_.pt() < pho1.pt() ) {
                pho2 = pho_;
            }
        }

        LorentzVector diphoton(pho1+pho2);

        bool passLeadingPhoton (pho1.pt()/diphoton.mass() > 0.35);
        bool passTrailingPhoton (pho2.pt()/diphoton.mass() > 0.25);
        bool passmassWindow(diphoton.mass()>105 && diphoton.mass()<160);


	mon.fillHisto("pt2mass_pho2vspho1_raw",tags, pho1.pt()/diphoton.mass(), pho2.pt()/diphoton.mass(), weight);

        if(passLeadingPhoton && passTrailingPhoton && passmassWindow) {

	    mon.fillHisto("nphotons_raw",tags, nGoodPhotons, weight);
            mon.fillHisto("diphoton_mass_raw", tags, diphoton.mass(), weight);
	    mon.fillHisto("diphoton_pt_raw", tags, diphoton.pt(), weight);

	    mon.fillHisto("leadingpho_pt_raw", tags, pho1.pt(), weight);
	    mon.fillHisto("trailingpho_pt_raw", tags, pho2.pt(), weight);
	    mon.fillHisto("leadingpho_eta_raw", tags, pho1.eta(), weight);
	    mon.fillHisto("trailingpho_eta_raw", tags, pho2.eta(), weight);
	    mon.fillHisto("leadingpho_phi_raw", tags, pho1.phi(), weight);
	    mon.fillHisto("trailingpho_phi_raw", tags, pho2.phi(), weight);



        }














    } // loop on all events END




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


    printf("\n");
    inputFile->Close();



    return 0;
}
