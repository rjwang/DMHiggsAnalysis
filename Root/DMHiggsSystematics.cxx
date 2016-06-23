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


void DMHiggsSystematics::declareVariables()
{
    myEvents->Branch("RunNumber",&RunNumber,"RunNumber/I");
    myEvents->Branch("LumiBlock",&LumiBlock,"LumiBlock/I");
    myEvents->Branch("EventNumber",&EventNumber,"EventNumber/I");
    myEvents->Branch("NPV",&NPV_var,"NPV_var/I");
    myEvents->Branch("mu",&mu_var,"mu_var/F");
    myEvents->Branch("initWeight",&initWeight_var,"initWeight_var/F");

    myEvents->Branch("nPhotons", &nPhotons,"nPhotons/I");
    myEvents->Branch("photon_Px", photon_Px,"photon_Px[nPhotons]/F");
    myEvents->Branch("photon_Py", photon_Py,"photon_Py[nPhotons]/F");
    myEvents->Branch("photon_Pz", photon_Pz,"photon_Pz[nPhotons]/F");
    myEvents->Branch("photon_E",  photon_E, "photon_E[nPhotons]/F");

    myEvents->Branch("nElectrons", &nElectrons,"nElectrons/I");
    myEvents->Branch("electron_Px", electron_Px,"electron_Px[nElectrons]/F");
    myEvents->Branch("electron_Py", electron_Py,"electron_Py[nElectrons]/F");
    myEvents->Branch("electron_Pz", electron_Pz,"electron_Pz[nElectrons]/F");
    myEvents->Branch("electron_E",  electron_E, "electron_E[nElectrons]/F");

    myEvents->Branch("nMuons", &nMuons,"nMuons/I");
    myEvents->Branch("muon_Px", muon_Px,"muon_Px[nMuons]/F");
    myEvents->Branch("muon_Py", muon_Py,"muon_Py[nMuons]/F");
    myEvents->Branch("muon_Pz", muon_Pz,"muon_Pz[nMuons]/F");
    myEvents->Branch("muon_E",  muon_E, "muon_E[nMuons]/F");

    myEvents->Branch("nJets", &nJets,"nJets/I");
    myEvents->Branch("jet_Px", jet_Px,"jet_Px[nJets]/F");
    myEvents->Branch("jet_Py", jet_Py,"jet_Py[nJets]/F");
    myEvents->Branch("jet_Pz", jet_Pz,"jet_Pz[nJets]/F");
    myEvents->Branch("jet_E", jet_E,"jet_E[nJets]/F");

    myEvents->Branch("met", &met,"met/F");
    myEvents->Branch("phi_met", &phi_met,"phi_met/F");
    myEvents->Branch("sumet", &sumet,"sumet/F");

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

    syslist.clear();
    syslist.push_back(""); //norminal
    syslist.push_back("EG_RESOLUTION_ALL__1down");
/*
    syslist.push_back("EG_RESOLUTION_ALL__1up");
    syslist.push_back("EG_SCALE_ALL__1down");
    syslist.push_back("EG_SCALE_ALL__1up");
    syslist.push_back("EL_EFF_ID_TotalCorrUncertainty__1down");
    syslist.push_back("EL_EFF_ID_TotalCorrUncertainty__1up");
    syslist.push_back("EL_EFF_Iso_TotalCorrUncertainty__1down");
    syslist.push_back("EL_EFF_Iso_TotalCorrUncertainty__1up");
    syslist.push_back("EL_EFF_Reco_TotalCorrUncertainty__1down");
    syslist.push_back("EL_EFF_Reco_TotalCorrUncertainty__1up");
    syslist.push_back("EL_EFF_TriggerEff_TotalCorrUncertainty__1down");
    syslist.push_back("EL_EFF_TriggerEff_TotalCorrUncertainty__1up");
    syslist.push_back("EL_EFF_Trigger_TotalCorrUncertainty__1down");
    syslist.push_back("EL_EFF_Trigger_TotalCorrUncertainty__1up");
    syslist.push_back("FT_EFF_Eigen_B_0__1down");
    syslist.push_back("FT_EFF_Eigen_B_0__1up");
    syslist.push_back("FT_EFF_Eigen_B_1__1down");
    syslist.push_back("FT_EFF_Eigen_B_1__1up");
    syslist.push_back("FT_EFF_Eigen_B_2__1down");
    syslist.push_back("FT_EFF_Eigen_B_2__1up");
    syslist.push_back("FT_EFF_Eigen_B_3__1down");
    syslist.push_back("FT_EFF_Eigen_B_3__1up");
    syslist.push_back("FT_EFF_Eigen_B_4__1down");
    syslist.push_back("FT_EFF_Eigen_B_4__1up");
    syslist.push_back("FT_EFF_Eigen_C_0__1down");
    syslist.push_back("FT_EFF_Eigen_C_0__1up");
    syslist.push_back("FT_EFF_Eigen_C_1__1down");
    syslist.push_back("FT_EFF_Eigen_C_1__1up");
    syslist.push_back("FT_EFF_Eigen_C_2__1down");
    syslist.push_back("FT_EFF_Eigen_C_2__1up");
    syslist.push_back("FT_EFF_Eigen_C_3__1down");
    syslist.push_back("FT_EFF_Eigen_C_3__1up");
    syslist.push_back("FT_EFF_Eigen_Light_0__1down");
    syslist.push_back("FT_EFF_Eigen_Light_0__1up");
    syslist.push_back("FT_EFF_Eigen_Light_10__1down");
    syslist.push_back("FT_EFF_Eigen_Light_10__1up");
    syslist.push_back("FT_EFF_Eigen_Light_11__1down");
    syslist.push_back("FT_EFF_Eigen_Light_11__1up");
    syslist.push_back("FT_EFF_Eigen_Light_1__1down");
    syslist.push_back("FT_EFF_Eigen_Light_1__1up");
    syslist.push_back("FT_EFF_Eigen_Light_2__1down");
    syslist.push_back("FT_EFF_Eigen_Light_2__1up");
    syslist.push_back("FT_EFF_Eigen_Light_3__1down");
    syslist.push_back("FT_EFF_Eigen_Light_3__1up");
    syslist.push_back("FT_EFF_Eigen_Light_4__1down");
    syslist.push_back("FT_EFF_Eigen_Light_4__1up");
    syslist.push_back("FT_EFF_Eigen_Light_5__1down");
    syslist.push_back("FT_EFF_Eigen_Light_5__1up");
    syslist.push_back("FT_EFF_Eigen_Light_6__1down");
    syslist.push_back("FT_EFF_Eigen_Light_6__1up");
    syslist.push_back("FT_EFF_Eigen_Light_7__1down");
    syslist.push_back("FT_EFF_Eigen_Light_7__1up");
    syslist.push_back("FT_EFF_Eigen_Light_8__1down");
    syslist.push_back("FT_EFF_Eigen_Light_8__1up");
    syslist.push_back("FT_EFF_Eigen_Light_9__1down");
    syslist.push_back("FT_EFF_Eigen_Light_9__1up");
    syslist.push_back("FT_EFF_extrapolation_from_charm__1down");
    syslist.push_back("FT_EFF_extrapolation_from_charm__1up");
    syslist.push_back("FT_EFF_extrapolation__1down");
    syslist.push_back("FT_EFF_extrapolation__1up");
    syslist.push_back("JET_JER_SINGLE_NP__1up");
    syslist.push_back("MET_JetTrk_ScaleDown");
    syslist.push_back("MET_JetTrk_ScaleUp");
    syslist.push_back("MET_SoftTrk_ResoPara");
    syslist.push_back("MET_SoftTrk_ResoPerp");
    syslist.push_back("MET_SoftTrk_ScaleDown");
    syslist.push_back("MET_SoftTrk_ScaleUp");
    syslist.push_back("MUONS_ID__1down");
    syslist.push_back("MUONS_ID__1up");
    syslist.push_back("MUONS_MS__1down");
    syslist.push_back("MUONS_MS__1up");
    syslist.push_back("MUONS_SCALE__1down");
    syslist.push_back("MUONS_SCALE__1up");
    syslist.push_back("MUON_EFF_STAT_LOWPT__1down");
    syslist.push_back("MUON_EFF_STAT_LOWPT__1up");
    syslist.push_back("MUON_EFF_STAT__1down");
    syslist.push_back("MUON_EFF_STAT__1up");
    syslist.push_back("MUON_EFF_SYS_LOWPT__1down");
    syslist.push_back("MUON_EFF_SYS_LOWPT__1up");
    syslist.push_back("MUON_EFF_SYS__1down");
    syslist.push_back("MUON_EFF_SYS__1up");
    syslist.push_back("MUON_EFF_TrigStatUncertainty__1down");
    syslist.push_back("MUON_EFF_TrigStatUncertainty__1up");
    syslist.push_back("MUON_EFF_TrigSystUncertainty__1down");
    syslist.push_back("MUON_EFF_TrigSystUncertainty__1up");
    syslist.push_back("MUON_ISO_STAT__1down");
    syslist.push_back("MUON_ISO_STAT__1up");
    syslist.push_back("MUON_ISO_SYS__1down");
    syslist.push_back("MUON_ISO_SYS__1up");
    syslist.push_back("PH_EFF_ID_Uncertainty__1down");
    syslist.push_back("PH_EFF_ID_Uncertainty__1up");
    syslist.push_back("PH_EFF_TRKISO_Uncertainty__1down");
    syslist.push_back("PH_EFF_TRKISO_Uncertainty__1up");
    syslist.push_back("PH_Iso_DDonoff");
    syslist.push_back("PRW_DATASF__1down");
    syslist.push_back("PRW_DATASF__1up");
*/

    SYSHIST = syslist.size();


    for(int hist=0; hist<SYSHIST; hist++) {

        TString Hname = "mgg";
        if(syslist[hist]!="") Hname += ("_"+syslist[hist]);

        myhisto[hist] = new TH1F( Hname.Data() , ";;Events" , 4, 0, 4);

        myhisto[hist]->GetXaxis()->SetBinLabel(1,"High #it{E}_{T}^{miss}, high #it{p}_{T}^{#gamma#gamma}");
        myhisto[hist]->GetXaxis()->SetBinLabel(2,"High #it{E}_{T}^{miss}, low #it{p}_{T}^{#gamma#gamma}");
        myhisto[hist]->GetXaxis()->SetBinLabel(3,"Intermediate #it{E}_{T}^{miss}");
        myhisto[hist]->GetXaxis()->SetBinLabel(4,"Rest");

    }


    // tree
    //myEvents = new TTree("DMHiggsSystematics","DMHiggsSystematics");
    //declareVariables();

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

        RunNumber = runNumber( *eventInfo );
        LumiBlock = lumiBlock( *eventInfo );
        EventNumber = eventNumber( *eventInfo );
        initWeight_var = initWeight(*HGameventInfo);
        NPV_var = NPV(*HGameventInfo);
        mu_var = mu(*HGameventInfo);


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


        //weights:  genWeights * PU * PVz
        float weight = 1.0;
        if(isMC()) weight *= HgammaAnalysis::weightFinal();//initWeight_var;


        //
        // Photon
        //

        if(photons.size()<2) continue; // 2 tight photons
        std::vector<LorentzVector> GoodPhotons;
        for(size_t gn=0; gn<photons.size(); gn++) {
            LorentzVector P4( photons[gn]->p4().Px()/1000.,  photons[gn]->p4().Py()/1000., photons[gn]->p4().Pz()/1000., photons[gn]->p4().E()/1000.);
            GoodPhotons.push_back(P4);
        }


        LorentzVector pho1 = GoodPhotons[0];
        LorentzVector pho2 = GoodPhotons[1];

        //find the leading and trailing photon
        if(pho1.pt()<pho2.pt()) {
            LorentzVector pho_ = pho1;
            pho1 = pho2;
            pho2 = pho_;
        }

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
            LorentzVector P4( jets[ijet]->p4().Px()/1000.,  jets[ijet]->p4().Py()/1000., jets[ijet]->p4().Pz()/1000., jets[ijet]->p4().E()/1000.);
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
                if(diphoton.mass()>122 && diphoton.mass()<128) myhisto[hist] -> Fill(0.,weight);
            } else {
                if(diphoton.mass()>122 && diphoton.mass()<128) myhisto[hist] -> Fill(1.,weight);
            }
        } else if(goodMET.pt()>50 && hardsum.pt()>40) {
            if(diphoton.mass()>122 && diphoton.mass()<128) myhisto[hist] -> Fill(2.,weight);
        } else if( diphoton.pt()>15) {
            if(diphoton.mass()>122 && diphoton.mass()<128) myhisto[hist] -> Fill(3.,weight);
        }






    }//systematics loop


    //myEvents->Fill();
    return EL::StatusCode::SUCCESS;
}



EL::StatusCode DMHiggsSystematics::finalize()
{


    m_outputFile->cd();
    //myEvents->Write();

    for(int hist=0; hist<SYSHIST; hist++) {
        myhisto[hist] ->Write();
    }

    m_outputFile->Close();



    HgammaAnalysis::finalize();

    return EL::StatusCode::SUCCESS;

}
