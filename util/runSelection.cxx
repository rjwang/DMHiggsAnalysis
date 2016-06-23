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



    double xsec     = runProcess.getNum("xsec");
    double BR     = runProcess.getNum("BR");
    double filtEff     = runProcess.getNum("filtEff");
    double kfactor     = runProcess.getNum("kfactor");


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


    TString outTreeUrl = outUrl + "/";
    outTreeUrl += outFileUrl + "_tree.root";
    //bool ifsaveEvents = (!isMC || isSignal || isMC_H125 || isMC_gamgam || isMC_gamjet);
    bool ifsaveEvents = true;

    TFile *otreefile=TFile::Open(outTreeUrl, "recreate");

    /*
        TTree *myEvents = new TTree("tree","tree");
        double mgg;
        //double weight_t, weight_all;
        myEvents->Branch("mgg",&mgg,"mgg/D");
        //myEvents->Branch("weight_t", &weight_t,"weight_t/D");
        //myEvents->Branch("weight_all", &weight_all,"weight_all/D");
    */

    TTree *myEvents_bin0 = new TTree("tree_bin0","tree_bin0");
    double mgg_bin0(0.), weight_bin0_t(1.), weight_bin0_all(1.);
    myEvents_bin0->Branch("mgg",&mgg_bin0,"mgg/D");
    myEvents_bin0->Branch("weight_t", &weight_bin0_t,"weight_t/D");
    myEvents_bin0->Branch("weight_all", &weight_bin0_all,"weight_all/D");


    TTree *myEvents_bin1 = new TTree("tree_bin1","tree_bin1");
    double mgg_bin1(0.), weight_bin1_t(1.), weight_bin1_all(1.);
    myEvents_bin1->Branch("mgg",&mgg_bin1,"mgg/D");
    myEvents_bin1->Branch("weight_t", &weight_bin1_t,"weight_t/D");
    myEvents_bin1->Branch("weight_all", &weight_bin1_all,"weight_all/D");


    TTree *myEvents_bin2 = new TTree("tree_bin2","tree_bin2");
    double mgg_bin2(0.), weight_bin2_t(1.), weight_bin2_all(1.);
    myEvents_bin2->Branch("mgg",&mgg_bin2,"mgg/D");
    myEvents_bin2->Branch("weight_t", &weight_bin2_t,"weight_t/D");
    myEvents_bin2->Branch("weight_all", &weight_bin2_all,"weight_all/D");


    TTree *myEvents_bin3 = new TTree("tree_bin3","tree_bin3");
    double mgg_bin3(0.), weight_bin3_t(1.), weight_bin3_all(1.);
    myEvents_bin3->Branch("mgg",&mgg_bin3,"mgg/D");
    myEvents_bin3->Branch("weight_t", &weight_bin3_t,"weight_t/D");
    myEvents_bin3->Branch("weight_all", &weight_bin3_all,"weight_all/D");


    TTree *myEvents_allbins = new TTree("tree_allbins","tree_allbins");
    double mgg_allbins(0.), weight_allbins_t(1.), weight_allbins_all(1.);
    myEvents_allbins->Branch("mgg",&mgg_allbins,"mgg/D");
    myEvents_allbins->Branch("weight_t", &weight_allbins_t,"weight_t/D");
    myEvents_allbins->Branch("weight_all", &weight_allbins_all,"weight_all/D");










    //##################################################################################
    //##########################    INITIATING HISTOGRAMS     ##########################
    //##################################################################################

    SmartSelectionMonitor mon;


    mon.addHistogram( new TH1F( "pu_sel", ";pileup;Events", 50,0,50) );
    mon.addHistogram( new TH1F( "puwgt_sel", ";pileup;Events", 50,0,50) );

    mon.addHistogram( new TH1F( "nvtx_sel",     ";Vertices;Events",50,0,50) );
    mon.addHistogram( new TH1F( "nvtxwgt_sel",  ";Vertices;Events",50,0,50) );

    //for MC normalization (to 1/pb)
    TH1F* Hcutflow  = (TH1F*) mon.addHistogram(  new TH1F ("cutflow"    , "cutflow"    ,6,0,6) ) ;

    mon.addHistogram( new TH1F( "diphoton_mass_sel",    ";#it{m}_{#gamma#gamma} [GeV];Events", 62,105,160) );
    mon.addHistogram( new TH1F( "diphoton_pt_sel",    ";#it{p}_{T}^{#gamma#gamma} [GeV];Events", 100,0,600) );
    mon.addHistogram( new TH1F( "met_sel",    	      ";#it{E}_{T}^{miss} [GeV];Events", 100,0,600) );
    mon.addHistogram( new TH1F( "balancedif_sel",     ";|E_{T}^{miss}-#it{p}_{T}^{#gamma#gamma}|/#it{p}_{T}^{#gamma#gamma};Events", 5,0,1.0) );
    mon.addHistogram( new TH1F( "metphi_sel",         ";#phi(#it{E}_{T}^{miss}) [rad];Events", 50,-TMath::Pi(),TMath::Pi()) );
    mon.addHistogram( new TH1F( "pthard_sel",         ";#it{p}_{T}^{hard} [GeV];Events", 100,0,600) );

    mon.addHistogram( new TH1F( "dphiGamGamMET_sel",   ";#Delta#it{#phi}(#it{p}_{T}^{#gamma#gamma},E_{T}^{miss});Events", 50,0,TMath::Pi()) );

    mon.addHistogram( new TH1F( "leadingpho_pt_sel",    ";Leading #it{p}_{T}^{#gamma} [GeV];Events", 100,0,600) );
    mon.addHistogram( new TH1F( "trailingpho_pt_sel",    ";Trailing #it{p}_{T}^{#gamma} [GeV];Events", 100,0,600) );
    mon.addHistogram( new TH1F( "leadingpho_eta_sel",    ";Leading #it{#eta}^{#gamma};Events", 50,-2.4,2.4) );
    mon.addHistogram( new TH1F( "trailingpho_eta_sel",    ";Trailing #it{#eta}^{#gamma};Events", 50,-2.4,2.4) );
    mon.addHistogram( new TH1F( "leadingpho_phi_sel",    ";Leading #it{#phi}^{#gamma};Events", 50,-1.*TMath::Pi(),TMath::Pi()) );
    mon.addHistogram( new TH1F( "trailingpho_phi_sel",    ";Trailing #it{#phi}^{#gamma};Events", 50,-1.*TMath::Pi(),TMath::Pi()) );


    mon.addHistogram( new TH1F( "diphoton_mass_bin0_sel",    ";#it{m}_{#gamma#gamma} [GeV];Events", 55,105,160) );
    mon.addHistogram( new TH1F( "diphoton_mass_bin1_sel",    ";#it{m}_{#gamma#gamma} [GeV];Events", 55,105,160) );
    mon.addHistogram( new TH1F( "diphoton_mass_bin2_sel",    ";#it{m}_{#gamma#gamma} [GeV];Events", 55,105,160) );
    mon.addHistogram( new TH1F( "diphoton_mass_bin3_sel",    ";#it{m}_{#gamma#gamma} [GeV];Events", 55,105,160) );



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

    h1 = (TH1F*) mon.addHistogram( new TH1F( "yields_finalsel",";;Events", 4,0,4) );
    h1->GetXaxis()->SetBinLabel(1,"High #it{E}_{T}^{miss}, high #it{p}_{T}^{#gamma#gamma}");
    h1->GetXaxis()->SetBinLabel(2,"High #it{E}_{T}^{miss}, low #it{p}_{T}^{#gamma#gamma}");
    h1->GetXaxis()->SetBinLabel(3,"Intermediate #it{E}_{T}^{miss}");
    h1->GetXaxis()->SetBinLabel(4,"Rest");


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

        //weights:  genWeights * PU * PVz
        float weight = 1.0;
        if(isMC) weight *= ev.initWeight;
        /*
        	//normalization weights
        	if(isMC){
        		weight *= (xsec*BR*filtEff*kfactor);
        		weight /= (sumInitialEvents/skim_eff);
        	}
        */
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
//        bool passmassWindow(diphoton.mass()>105 && diphoton.mass()<160);
        bool passmassWindow(diphoton.mass()>105 && diphoton.mass()<180);


        LorentzVector met = phys.met;
        double dphiGamGamMET=fabs(deltaPhi(diphoton.phi(),met.phi()));
        double balanceDif = fabs(1-met.pt()/diphoton.pt());




        if(!passLeadingPhoton || !passTrailingPhoton || !passmassWindow) continue;


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


        mon.fillHisto("pu_sel",   tags, ev.mu,      1.0);
        mon.fillHisto("puwgt_sel",tags, ev.mu,      weight);
        //# of PV
        mon.fillHisto("nvtx_sel",   tags, phys.nvtx,      1.0);
        mon.fillHisto("nvtxwgt_sel",tags, phys.nvtx,      weight);


        mon.fillHisto("nphotons_sel",tags, nGoodPhotons, weight);
        mon.fillHisto("diphoton_mass_sel", tags, diphoton.mass(), weight);
        mon.fillHisto("diphoton_pt_sel", tags, diphoton.pt(), weight);
        mon.fillHisto("met_sel", tags, met.pt(), weight);
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


        if(met.pt()>100) {
            if(diphoton.pt()>100) {
                if(diphoton.mass()>122 && diphoton.mass()<128) mon.fillHisto("yields_finalsel",tags, 0, weight);
                mon.fillHisto("diphoton_mass_bin0_sel", tags, diphoton.mass(), weight);
            } else {
                if(diphoton.mass()>122 && diphoton.mass()<128) mon.fillHisto("yields_finalsel",tags, 1, weight);
                mon.fillHisto("diphoton_mass_bin1_sel", tags, diphoton.mass(), weight);
            }
        } else if(met.pt()>50 && hardsum.pt()>40) {
            if(diphoton.mass()>122 && diphoton.mass()<128) mon.fillHisto("yields_finalsel",tags, 2, weight);
            mon.fillHisto("diphoton_mass_bin2_sel", tags, diphoton.mass(), weight);
        } else if( diphoton.pt()>15) {
            if(diphoton.mass()>122 && diphoton.mass()<128) mon.fillHisto("yields_finalsel",tags, 3, weight);
            mon.fillHisto("diphoton_mass_bin3_sel", tags, diphoton.mass(), weight);
        }












        if( ifsaveEvents && (diphoton.mass()<120 || diphoton.mass()>130) ) {

            if(met.pt()>100) {
                if(diphoton.pt()>100) {

                    mgg_bin0 = diphoton.mass();
                    weight_bin0_t = weight;
                    weight_bin0_all = weight;
                    if(isMC) {
                        weight_bin0_all *= (xsec*BR*filtEff*kfactor);
                        weight_bin0_all /= (sumInitialEvents/skim_eff);
                    }
                    myEvents_bin0->Fill();

                } else {

                    mgg_bin1 = diphoton.mass();
                    weight_bin1_t = weight;
                    weight_bin1_all = weight;
                    if(isMC) {
                        weight_bin1_all *= (xsec*BR*filtEff*kfactor);
                        weight_bin1_all /= (sumInitialEvents/skim_eff);
                    }
                    myEvents_bin1->Fill();

                }
            } else if(met.pt()>50 && hardsum.pt()>40) {


                mgg_bin2 = diphoton.mass();
                weight_bin2_t = weight;
                weight_bin2_all = weight;
                if(isMC) {
                    weight_bin2_all *= (xsec*BR*filtEff*kfactor);
                    weight_bin2_all /= (sumInitialEvents/skim_eff);
                }
                myEvents_bin2->Fill();


            } else if(diphoton.pt()>15) {

                mgg_bin3 = diphoton.mass();
                weight_bin3_t = weight;
                weight_bin3_all = weight;
                if(isMC) {
                    weight_bin3_all *= (xsec*BR*filtEff*kfactor);
                    weight_bin3_all /= (sumInitialEvents/skim_eff);
                }
                myEvents_bin3->Fill();

            }



                mgg_allbins = diphoton.mass();
                weight_allbins_t = weight;
                weight_allbins_all = weight;
                if(isMC) {
                    weight_allbins_all *= (xsec*BR*filtEff*kfactor);
                    weight_allbins_all /= (sumInitialEvents/skim_eff);
                }
                myEvents_allbins->Fill();


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

    // saving hgg tree to output file
    if(ifsaveEvents) {
        otreefile->cd();
        myEvents_bin0->Write();
        myEvents_bin1->Write();
        myEvents_bin2->Write();
        myEvents_bin3->Write();
        myEvents_allbins->Write();
        printf("Tree saved in %s\n", outTreeUrl.Data());
    }
    otreefile->Close();
    if(!ifsaveEvents) gSystem->Exec("rm " + outTreeUrl);

    printf("\n");
    inputFile->Close();



    return 0;
}
