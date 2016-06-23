#include "DMHiggsAnalysis/DMHiggsAnalysis.h"
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




// this is needed to distribute the algorithm to the workers
ClassImp(DMHiggsAnalysis)



DMHiggsAnalysis::DMHiggsAnalysis(const char *name)
    : HgammaAnalysis(name)
{
    // Here you put any code for the base initialization of variables,
    // e.g. initialize all pointers to 0.  Note that you should only put
    // the most basic initialization here, since this method will be
    // called on both the submission and the worker node.  Most of your
    // initialization code will go into histInitialize() and
    // initialize().

}



DMHiggsAnalysis::~DMHiggsAnalysis()
{
    // Here you delete any memory you allocated during your analysis.
}



EL::StatusCode DMHiggsAnalysis::createOutput()
{
    // Here you setup the histograms needed for you analysis. This method
    // gets called after the Handlers are initialized, so that the systematic
    // registry is already filled.



    //Create a TTree
    //  TFile *outfile = wk()->getOutputFile("MxAOD");



    return EL::StatusCode::SUCCESS;
}


void DMHiggsAnalysis::declareVariables()
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


EL::StatusCode DMHiggsAnalysis::initialize()
{

    using namespace std;
    HgammaAnalysis::initialize();
    std::string inputfileName = wk()->inputFileName();
    //currentfilename = inputfileName;
    inputfileName.replace(inputfileName.find(".MxAOD") , -1, "") ;
    TString TagName = inputfileName;
    TagName.ReplaceAll("mc15c.","");
    inputfileName.append(".NTuple.root");
    TString InputfileName = inputfileName;

    m_outputFile = TFile::Open(InputfileName.Data(),"RECREATE");

    if(TagName.Contains("Pythia8_WH125") || TagName.Contains("Pythia8_ZH125") ) TagName += "_big";
    if(TagName.Contains("Pythia8_1Hard1BremDP40")) TagName = "Pythia81Hard1BremDP40";

    // copy cutflow hists
    if( isMC() ) {
        CutFlow_ = (TH1F*) getCutFlowHistogram("CutFlow_"+TagName,"")->Clone("CutFlow");
        CutFlow_noDalitz_ = (TH1F*) getCutFlowHistogram("CutFlow_"+TagName,"_noDalitz")->Clone("CutFlow_noDalitz");
        CutFlow_weighted_ = (TH1F*) getCutFlowHistogram("CutFlow_"+TagName,"_weighted")->Clone("CutFlow_weighted");
        CutFlow_noDalitz_weighted_ = (TH1F*) getCutFlowHistogram("CutFlow_"+TagName,"_noDalitz_weighted")->Clone("CutFlow_noDalitz_weighted");
    }


    // tree
    myEvents = new TTree("DMHiggsAnalysis","DMHiggsAnalysis");
    declareVariables();

    return EL::StatusCode::SUCCESS;
}

EL::StatusCode DMHiggsAnalysis::execute()
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





    // photons
    const xAOD::EventInfo* HGameventInfo = 0 ;
    if(event()->retrieve(HGameventInfo,"HGamEventInfo").isFailure() )
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
    if(isMC()) initWeight_var = HgammaAnalysis::weightFinal(); //initWeight(*HGameventInfo);
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


/*
    if(!HgammaAnalysis::pass( &photons, &electrons, &muons, &jets)) return EL::StatusCode::SUCCESS;
    if(!HgammaAnalysis::passJetEventCleaning() ) return EL::StatusCode::SUCCESS;
*/


    //move this overlap checking to selection level
    overlapHandler()->removeOverlap(photons,jets,electrons, muons);


    //
    // Photon
    //
    nPhotons=0;
    for(size_t gn=0; gn<photons.size(); gn++) {
        photon_Px[nPhotons] = photons[gn]->p4().Px();
        photon_Py[nPhotons] = photons[gn]->p4().Py();
        photon_Pz[nPhotons] = photons[gn]->p4().Pz();
        photon_E[nPhotons] = photons[gn]->p4().E();
        nPhotons++;

    }


    //
    // Electron
    //
    nElectrons=0;
    for(size_t gn=0; gn<electrons.size(); gn++) {
        electron_Px[nElectrons] = electrons[gn]->p4().Px();
        electron_Py[nElectrons] = electrons[gn]->p4().Py();
        electron_Pz[nElectrons] = electrons[gn]->p4().Pz();
        electron_E[nElectrons] = electrons[gn]->e();
        nElectrons++;
    }



    //
    // Muon
    //
    nMuons=0;
    for(size_t gn=0; gn<muons.size(); gn++) {
        muon_Px[nMuons] = muons[gn]->p4().Px();
        muon_Py[nMuons] = muons[gn]->p4().Py();
        muon_Pz[nMuons] = muons[gn]->p4().Pz();
        muon_E[nMuons] = muons[gn]->e();

        nMuons++;
    }


    //
    // Jets
    //
    nJets=0;
    for(size_t gn=0; gn<jets.size(); gn++) {
        jet_Px[nJets] = jets[gn]->px();
        jet_Py[nJets] = jets[gn]->py();
        jet_Pz[nJets] = jets[gn]->pz();
        jet_E[nJets]  = jets[gn]->e();

        nJets++;
    }



    //
    // MET
    //
    met = met_TST(*HGameventInfo);
    sumet = met_sumet(*HGameventInfo);
    phi_met = met_phi(*HGameventInfo);



    myEvents->Fill();
    return EL::StatusCode::SUCCESS;
}



EL::StatusCode DMHiggsAnalysis::finalize()
{


    m_outputFile->cd();
    if( isMC() ) {
        CutFlow_->Write();
        CutFlow_noDalitz_->Write();
        CutFlow_weighted_->Write();
        CutFlow_noDalitz_weighted_->Write();
    }
    myEvents->Write();
    m_outputFile->Close();



    HgammaAnalysis::finalize();

    return EL::StatusCode::SUCCESS;

}
