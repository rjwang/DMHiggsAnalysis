#include "DMHiggsAnalysis/DMHiggsSystematics.h"
#include "HGamAnalysisFramework/HgammaIncludes.h"




//Egamma includes
#include "xAODEgamma/PhotonxAODHelpers.h"
#include "xAODEgamma/EgammaDefs.h"
#include "xAODEgamma/EgammaEnums.h"
#include "xAODEgamma/EgammaxAODHelpers.h"

// Trigger includes
#include "TrigConfxAOD/xAODConfigTool.h"
#include "TrigDecisionTool/TrigDecisionTool.h"

// ETmissHandler
#include "xAODBase/IParticleHelpers.h"

// Truth Particle
#include "xAODTruth/TruthParticleContainer.h"
#include "xAODTruth/TruthParticleAuxContainer.h"
#include <assert.h>


#if !defined(__CINT__) || defined(__MAKECINT__)

#include <string.h>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <set>
#include <cmath>

#include "Math/LorentzVector.h"
#include "TMath.h"
#include "TVector2.h"
#include "TTree.h"
#include "TLorentzVector.h"

#endif

typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > LorentzVector;
typedef std::vector<LorentzVector> LorentzVectorCollection;



// this is needed to distribute the algorithm to the workers
ClassImp(DMHiggsSystematics)



DMHiggsSystematics::DMHiggsSystematics(const char *name)
    : HgammaAnalysis(name)
{
    // Here you put any code for the base initialization of variables,
    // e.g. initialize all pointers to 0.  Note that you should only put
    // the most basic initialization here, since this method will be
    // called on both the submission and the worker node.  Most of your
    // initialization code will go into histInitialize() and
    // initialize().

}



DMHiggsSystematics::~DMHiggsSystematics()
{
    // Here you delete any memory you allocated during your analysis.
}



EL::StatusCode DMHiggsSystematics::createOutput()
{
    // Here you setup the histograms needed for you analysis. This method
    // gets called after the Handlers are initialized, so that the systematic
    // registry is already filled.



    //Create a TTree
    //  TFile *outfile = wk()->getOutputFile("MxAOD");



    return EL::StatusCode::SUCCESS;
}



EL::StatusCode DMHiggsSystematics::initialize()
{

    HgammaAnalysis::initialize();
    std::string inputfileName = wk()->inputFileName();
    //currentfilename = inputfileName;
    inputfileName.replace(inputfileName.find(".MxAOD") , -1, "") ;
    TString TagName = inputfileName;
    inputfileName.append(".NTuple.root");

    m_outputFile = TFile::Open(inputfileName.c_str(),"RECREATE");





    //#############################
    // doing systematics
    //#############################

    // Set the function type:
//    TString sysListDir = "/afs/cern.ch/work/r/rewang/monoHiggs/hgam_012_000247/DMHiggsAnalysis/data/SysList_data16_h012_rel20p7.txt";
    TString sysListDir = "/afs/cern.ch/work/r/rewang/monoHiggs/hgam_012_000247/DMHiggsAnalysis/data/SysList_data16_h012_rel20p7_noMuonElev2.txt";
    std::ifstream sysListFile (sysListDir.Data());
    syslist.clear();
    syslist.push_back(""); //norminal

    if (sysListFile.is_open()) {
        std::string line;
        while ( getline (sysListFile,line) ) {
            if(line == "")continue;
            TString tempSys = line;
            std::cout << "using sys: " << tempSys << std::endl;
            syslist.push_back(tempSys);
            //break; //for testing
        }
    } else {
        std::cout << "Unable to open systematic list: " << sysListDir <<std::endl;
        abort();
    }
    sysListFile.close();

    SYSHIST = syslist.size();
    for(int hist=0; hist<SYSHIST; hist++) {

        TString Hname = "mgg";
        if(syslist[hist]!="") Hname += ("_"+syslist[hist]);

        myhisto[hist] = new TH1F( Hname.Data() , ";;Events" , 14, 0, 14);

        myhisto[hist]->GetXaxis()->SetBinLabel(1,"High #it{E}_{T}^{miss}, high #it{p}_{T}^{#gamma#gamma}");
        myhisto[hist]->GetXaxis()->SetBinLabel(2,"High #it{E}_{T}^{miss}, low #it{p}_{T}^{#gamma#gamma}");
        myhisto[hist]->GetXaxis()->SetBinLabel(3,"Intermediate #it{E}_{T}^{miss}");
        myhisto[hist]->GetXaxis()->SetBinLabel(4,"Rest Category");
        //model independent limits
        myhisto[hist]->GetXaxis()->SetBinLabel(5,"#it{E}_{T}^{miss} > 10, #it{p}_{T}^{#gamma#gamma} > 10");
        myhisto[hist]->GetXaxis()->SetBinLabel(6,"#it{E}_{T}^{miss} > 20, #it{p}_{T}^{#gamma#gamma} > 20");
        myhisto[hist]->GetXaxis()->SetBinLabel(7,"#it{E}_{T}^{miss} > 30, #it{p}_{T}^{#gamma#gamma} > 30");
        myhisto[hist]->GetXaxis()->SetBinLabel(8,"#it{E}_{T}^{miss} > 40, #it{p}_{T}^{#gamma#gamma} > 40");
        myhisto[hist]->GetXaxis()->SetBinLabel(9,"#it{E}_{T}^{miss} > 50, #it{p}_{T}^{#gamma#gamma} > 50");
        myhisto[hist]->GetXaxis()->SetBinLabel(10,"#it{E}_{T}^{miss} > 60, #it{p}_{T}^{#gamma#gamma} > 60");
        myhisto[hist]->GetXaxis()->SetBinLabel(11,"#it{E}_{T}^{miss} > 70, #it{p}_{T}^{#gamma#gamma} > 70");
        myhisto[hist]->GetXaxis()->SetBinLabel(12,"#it{E}_{T}^{miss} > 80, #it{p}_{T}^{#gamma#gamma} > 80");
        myhisto[hist]->GetXaxis()->SetBinLabel(13,"#it{E}_{T}^{miss} > 90, #it{p}_{T}^{#gamma#gamma} > 90");
        myhisto[hist]->GetXaxis()->SetBinLabel(14,"#it{E}_{T}^{miss} > 100, #it{p}_{T}^{#gamma#gamma} > 100");
        myhisto[hist]->GetXaxis()->SetBinLabel(15,"#it{E}_{T}^{miss} > 110, #it{p}_{T}^{#gamma#gamma} > 110");
        myhisto[hist]->GetXaxis()->SetBinLabel(16,"#it{E}_{T}^{miss} > 120, #it{p}_{T}^{#gamma#gamma} > 120");


    }



    return EL::StatusCode::SUCCESS;
}

EL::StatusCode DMHiggsSystematics::execute()
{
    // Here you do everything that needs to be done on every single
    // events, e.g. read input variables, apply cuts, and fill
    // histograms and trees.  This is where most of your actual analysis
    // code will go.

    // Important to keep this, so that internal tools / event variables
    // are filled properly.
    using namespace std;

    HgammaAnalysis::execute();



    SG::AuxElement::Accessor<unsigned int> runNumber("runNumber");
    SG::AuxElement::Accessor<unsigned int> lumiBlock("lumiBlock");
    SG::AuxElement::Accessor<unsigned long long> eventNumber("eventNumber");
    SG::AuxElement::Accessor< std::vector<float> > mcEventWeights("mcEventWeights");

    SG::AuxElement::Accessor<int> NPV("numberOfPrimaryVertices");
    SG::AuxElement::Accessor<float> mu("mu");
    SG::AuxElement::Accessor<float> initWeight("weightInitial");
    SG::AuxElement::Accessor<char> isPassed("isPassed");
    SG::AuxElement::Accessor<char> isPassedJetEventClean("isPassedJetEventClean");
    SG::AuxElement::Accessor<char> isDalitz("isDalitz");
    SG::AuxElement::Accessor<float> met_TST("met_TST");
    SG::AuxElement::Accessor<float> met_sumet("sumet_TST");
    SG::AuxElement::Accessor<float> met_phi("phi_TST");




    // Systematics loop

    for(int hist=0; hist<SYSHIST; hist++) {

        TString HGamEventInfoTag = "HGamEventInfo";
        if(syslist[hist]!="") HGamEventInfoTag += ("_"+syslist[hist]);



        // photons
        const xAOD::EventInfo* HGameventInfo = 0 ;
        if(event()->retrieve(HGameventInfo,HGamEventInfoTag.Data()).isFailure() )
            HG::fatal("Cannot retrieve event Info .");

        const xAOD::EventInfo* eventInfo = 0 ;
        if(event()->retrieve(eventInfo,"EventInfo").isFailure() )
            HG::fatal("Cannot retrieve event Info .");

        isPassed_var = isPassed(*HGameventInfo) == 1 ? 1 : 0 ;
        isPassedJetEventClean_var = isPassedJetEventClean(*HGameventInfo) == 1 ? 1 : 0 ;
        if(!isPassed_var || !isPassedJetEventClean_var) return EL::StatusCode::SUCCESS;

        if(isMC()) {
            isDalitz_var = isDalitz(*HGameventInfo) == 1 ? 1 : 0 ;
            if(isDalitz_var) return EL::StatusCode::SUCCESS;
        }



        /*
                RunNumber = runNumber( *eventInfo );
                LumiBlock = lumiBlock( *eventInfo );
                EventNumber = eventNumber( *eventInfo );
                initWeight_var = initWeight(*HGameventInfo);
                NPV_var = NPV(*HGameventInfo);
                mu_var = mu(*HGameventInfo);
        */

        xAOD::PhotonContainer photons_H = photonHandler()->getCorrectedContainer() ;
        xAOD::ElectronContainer electrons_H = electronHandler()->getCorrectedContainer() ;
        xAOD::MuonContainer muons_H = muonHandler()->getCorrectedContainer() ;
        xAOD::JetContainer jets_H = jetHandler()->getCorrectedContainer() ;


        xAOD::PhotonContainer photons = photonHandler()->applySelection(photons_H);
        xAOD::ElectronContainer electrons = electronHandler()->applySelection(electrons_H);
        xAOD::MuonContainer muons = muonHandler()->applySelection(muons_H);
        xAOD::JetContainer jets = jetHandler()->applySelection(jets_H);


        //move this overlap checking to selection level
        overlapHandler()->removeOverlap(photons,jets,electrons, muons);


        //if(!HgammaAnalysis::pass( &photons, &electrons, &muons, &jets)) return EL::StatusCode::SUCCESS;
        //if(!HgammaAnalysis::passJetEventCleaning() ) return EL::StatusCode::SUCCESS;



        HgammaAnalysis::setSelectedObjects(&photons, &electrons, &muons, &jets);

        //weights:  genWeights * PU * PVz
        float weight = 1.0;
        //if(isMC()) weight *= HgammaAnalysis::weightFinal();
        if(isMC()) weight *= HgammaAnalysis::weight();
        // mcWeight_var * pileupWeight_var * vertexWeight_var * weigth from Photon SF






        //
        // Photon
        //

        if(photons.size()<2) continue; // 2 tight photons
        std::vector<LorentzVector> GoodPhotons;
        for(size_t gn=0; gn<photons.size(); gn++) {

            double pt = photons[gn]->p4().Pt();
            double eta = photons[gn]->p4().Eta();
            double phi = photons[gn]->p4().Phi();
            double e = photons[gn]->p4().E();

            TLorentzVector pho;
            pho.Clear();
            pho.SetPtEtaPhiE(pt,eta,phi,e);

            LorentzVector P4( pho.Px()/1000., pho.Py()/1000., pho.Pz()/1000., pho.E()/1000.);
            GoodPhotons.push_back(P4);
        }


        LorentzVector pho1 = GoodPhotons[0];
        LorentzVector pho2 = GoodPhotons[1];

        /*
            //find the leading and trailing photon
            if(pho1.pt()<pho2.pt()) {
                LorentzVector pho_ = pho1;
                pho1 = pho2;
                pho2 = pho_;
            }
        */

        LorentzVector diphoton(pho1+pho2);
        bool passLeadingPhoton (pho1.pt()/diphoton.mass() > 0.35);
        bool passTrailingPhoton (pho2.pt()/diphoton.mass() > 0.25);
        bool passmassWindow(diphoton.mass()>105 && diphoton.mass()<160);

        if(!passLeadingPhoton || !passTrailingPhoton || !passmassWindow) continue;



        //
        // Electron
        //
//                nElectrons=0;
//                for(size_t gn=0; gn<electrons.size(); gn++) {
//                    electron_Px[nElectrons] = electrons[gn]->p4().Px();
//                    electron_Py[nElectrons] = electrons[gn]->p4().Py();
//                    electron_Pz[nElectrons] = electrons[gn]->p4().Pz();
//                    electron_E[nElectrons] = electrons[gn]->e();
//                    nElectrons++;
//                }
//
//
//
//                //
//                // Muon
//                //
//                nMuons=0;
//                for(size_t gn=0; gn<muons.size(); gn++) {
//                    muon_Px[nMuons] = muons[gn]->p4().Px();
//                    muon_Py[nMuons] = muons[gn]->p4().Py();
//                    muon_Pz[nMuons] = muons[gn]->p4().Pz();
//                    muon_E[nMuons] = muons[gn]->e();
//
//                    nMuons++;
//                }
//
//
        //
        // Jets
        //
        std::vector<LorentzVector> GoodJets;
        for(size_t ijet=0; ijet<jets.size(); ijet++) {

            double pt = jets[ijet]->p4().Pt();
            double eta = jets[ijet]->p4().Eta();
            double phi = jets[ijet]->p4().Phi();
            double e = jets[ijet]->p4().E();

            TLorentzVector pho;
            pho.Clear();
            pho.SetPtEtaPhiE(pt,eta,phi,e);

            LorentzVector P4( pho.Px()/1000., pho.Py()/1000., pho.Pz()/1000., pho.E()/1000.);
            GoodJets.push_back(P4);
        }




        LorentzVector hardsum(0,0,0,0);
        for(size_t ijet=0; ijet<GoodJets.size(); ijet++) {
            hardsum += GoodJets[ijet];
        }
        for(size_t ipho=0; ipho<GoodPhotons.size(); ipho++) {
            hardsum += GoodPhotons[ipho];
        }


        //
        // MET
        //
        met = met_TST(*HGameventInfo);
        sumet = met_sumet(*HGameventInfo);
        phi_met = met_phi(*HGameventInfo);


        LorentzVector goodMET =  LorentzVector( met*cos(phi_met)/1000., met*sin(phi_met)/1000., 0, met/1000. );

        if(goodMET.pt()>100) {
            if(diphoton.pt()>100) {
                myhisto[hist] -> Fill(0.,weight);
            } else {
                myhisto[hist] -> Fill(1.,weight);
            }
        } else if(goodMET.pt()>50 && hardsum.pt()>40) {
            myhisto[hist] -> Fill(2.,weight);
        } else if( diphoton.pt()>15) {
            myhisto[hist] -> Fill(3.,weight);
        }


        if(goodMET.pt()>10 && diphoton.pt()>10) myhisto[hist] -> Fill(4, weight);
        if(goodMET.pt()>20 && diphoton.pt()>20) myhisto[hist] -> Fill(5, weight);
        if(goodMET.pt()>30 && diphoton.pt()>30) myhisto[hist] -> Fill(6, weight);
        if(goodMET.pt()>40 && diphoton.pt()>40) myhisto[hist] -> Fill(7, weight);
        if(goodMET.pt()>50 && diphoton.pt()>50) myhisto[hist] -> Fill(8, weight);
        if(goodMET.pt()>60 && diphoton.pt()>60) myhisto[hist] -> Fill(9, weight);
        if(goodMET.pt()>70 && diphoton.pt()>70) myhisto[hist] -> Fill(10, weight);
        if(goodMET.pt()>80 && diphoton.pt()>80) myhisto[hist] -> Fill(11, weight);
        if(goodMET.pt()>90 && diphoton.pt()>90) myhisto[hist] -> Fill(12, weight);
        if(goodMET.pt()>100 && diphoton.pt()>100) myhisto[hist] -> Fill(13, weight);
        if(goodMET.pt()>110 && diphoton.pt()>110) myhisto[hist] -> Fill(14, weight);
        if(goodMET.pt()>120 && diphoton.pt()>120) myhisto[hist] -> Fill(15, weight);







    }//systematics loop


    return EL::StatusCode::SUCCESS;
}



EL::StatusCode DMHiggsSystematics::finalize()
{


    m_outputFile->cd();

    for(int hist=0; hist<SYSHIST; hist++) {
        myhisto[hist] ->Write();
    }
    m_outputFile->Close();


    HgammaAnalysis::finalize();

    return EL::StatusCode::SUCCESS;

}
