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

    //for MC normalization (to 1/pb)
    TH1F* Hcutflow  = (TH1F*) mon.addHistogram(  new TH1F ("cutflow"    , "cutflow"    ,6,0,6) ) ;

//    mon.addHistogram( new TH2F( "pt2mass_pho2vspho1_raw",";#it{p}_{T}^{#gamma1}/#it{m}_{#gamma#gamma};#it{p}_{T}^{#gamma2}/#it{m}_{#gamma#gamma};Events",100,0,1., 100,0,1.) );

    mon.addHistogram( new TH1F( "diphoton_mass_raw",    ";#it{m}_{#gamma#gamma} [GeV];Events", 62,105,160) );
    mon.addHistogram( new TH1F( "diphoton_pt_raw",    ";#it{p}_{T}^{#gamma#gamma} [GeV];Events", 100,0,600) );
    mon.addHistogram( new TH1F( "met_raw",    	      ";#it{E}_{T}^{miss} [GeV];Events", 100,0,600) );
    mon.addHistogram( new TH1F( "metphi_raw",         ";#phi(#it{E}_{T}^{miss}) [rad];Events", 50,-TMath::Pi(),TMath::Pi()) );

    mon.addHistogram( new TH1F( "dphiGamGamMET_raw",   ";#Delta#it{#phi}(#it{p}_{T}^{#gamma#gamma},E_{T}^{miss});Events", 50,0,TMath::Pi()) );


    mon.addHistogram( new TH1F( "leadingpho_pt_raw",    ";Leading #it{p}_{T}^{#gamma} [GeV];Events", 100,0,600) );
    mon.addHistogram( new TH1F( "trailingpho_pt_raw",    ";Trailing #it{p}_{T}^{#gamma} [GeV];Events", 100,0,600) );
    mon.addHistogram( new TH1F( "leadingpho_eta_raw",    ";Leading #it{#eta}^{#gamma};Events", 50,-2.4,2.4) );
    mon.addHistogram( new TH1F( "trailingpho_eta_raw",    ";Trailing #it{#eta}^{#gamma};Events", 50,-2.4,2.4) );
    mon.addHistogram( new TH1F( "leadingpho_phi_raw",    ";Leading #it{#phi}^{#gamma};Events", 50,-1.*TMath::Pi(),TMath::Pi()) );
    mon.addHistogram( new TH1F( "trailingpho_phi_raw",    ";Trailing #it{#phi}^{#gamma};Events", 50,-1.*TMath::Pi(),TMath::Pi()) );



    TH1F *h1 = (TH1F*) mon.addHistogram( new TH1F( "nphotons_raw", ";Photon multiplicity;Events", 3,2,5) );
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



    Hcutflow->SetBinContent(1,sumInitialEvents/skim_eff);



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
        if(isMC) weight *= ev.initWeight;

        // add PhysicsEvent_t class, get all tree to physics objects
        PhysicsEvent_t phys=getPhysicsEventFrom(ev);





        // looping photons
        int nGoodPhotons(0);
        std::vector<std::pair<int,LorentzVector> > GoodPhotons;
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
            GoodPhotons.push_back(goodpho);
        }

        if(nGoodPhotons<2) continue; // 2 tight photons



        LorentzVector pho1 = GoodPhotons[0].second;
        LorentzVector pho2 = GoodPhotons[1].second;

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

        LorentzVector diphoton(pho1+pho2);

        bool passLeadingPhoton (pho1.pt()/diphoton.mass() > 0.35);
        bool passTrailingPhoton (pho2.pt()/diphoton.mass() > 0.25);
        bool passmassWindow(diphoton.mass()>105 && diphoton.mass()<160);
        LorentzVector met = phys.met;
        double dphiGamGamMET=fabs(deltaPhi(diphoton.phi(),met.phi()));


        mon.fillHisto("trailingpho_phi_raw", tags, pho2.phi(), weight);


        if(!passLeadingPhoton || !passTrailingPhoton || !passmassWindow) continue;

        mon.fillHisto("pu_raw",   tags, ev.mu,      1.0);
        mon.fillHisto("puwgt_raw",tags, ev.mu,      weight);
        //# of PV
        mon.fillHisto("nvtx_raw",   tags, phys.nvtx,      1.0);
        mon.fillHisto("nvtxwgt_raw",tags, phys.nvtx,      weight);


        mon.fillHisto("nphotons_raw",tags, nGoodPhotons, weight);
        mon.fillHisto("diphoton_mass_raw", tags, diphoton.mass(), weight);
        mon.fillHisto("diphoton_pt_raw", tags, diphoton.pt(), weight);
        mon.fillHisto("met_raw", tags, met.pt(), weight);
        mon.fillHisto("metphi_raw", tags, met.phi(), weight);
        mon.fillHisto("dphiGamGamMET_raw",tags, dphiGamGamMET, weight);

        mon.fillHisto("leadingpho_pt_raw", tags, pho1.pt(), weight);
        mon.fillHisto("trailingpho_pt_raw", tags, pho2.pt(), weight);
        mon.fillHisto("leadingpho_eta_raw", tags, pho1.eta(), weight);
        mon.fillHisto("trailingpho_eta_raw", tags, pho2.eta(), weight);
        mon.fillHisto("leadingpho_phi_raw", tags, pho1.phi(), weight);


        //
        // Electrons
        //
        std::vector<std::pair<int,LorentzVector> > goodElectrons;
        for(size_t iele=0; iele<phys.electrons.size(); iele++) {
            LorentzVector ele=phys.electrons[iele];
            int eleid = phys.electrons[iele].id;
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

            std::pair <int,LorentzVector> goodlep;
            goodlep = std::make_pair(eleid,ele);
            goodElectrons.push_back(goodlep);
        }


        //
        // Jets
        //
        PhysicsObjectJetCollection GoodJets;
        for(size_t ijet=0; ijet<phys.jets.size(); ijet++) {

            if(phys.jets[ijet].pt()<25) continue;
            if(fabs(phys.jets[ijet].eta())>4.4) continue;

            //jet ID
            if(fabs(phys.jets[ijet].eta()) < 2.4 && phys.jets[ijet].JVT < 0.64) continue;

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

            GoodJets.push_back(phys.jets[ijet]);
        }




        //
        // Muons
        //

        std::vector<std::pair<int,LorentzVector> > GoodMuons;
        for(size_t imu=0; imu<phys.muons.size(); imu++) {
            LorentzVector mu=phys.muons[imu];
            int muid = phys.muons[imu].id;
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


            std::pair <int,LorentzVector> goodlep;
            goodlep = std::make_pair(muid,mu);
            GoodMuons.push_back(goodlep);
        }




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








        mon.fillHisto("nelectrons_sel",tags, GoodElectrons.size(), weight);
        mon.fillHisto("nmuons_sel",tags, GoodMuons.size(), weight);
        mon.fillHisto("njets_sel",tags, GoodJets.size(), weight);
















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
