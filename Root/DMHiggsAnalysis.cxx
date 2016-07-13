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


    //histoStore()->createTH1F("mu", 50, 0, 50);

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
    myEvents->Branch("totWeight",&totWeight_var,"totWeight_var/F");
    myEvents->Branch("lumiXsecWeight",&lumiXsecWeight_var,"lumiXsecWeight_var/F");
    myEvents->Branch("evtWeight",&evtWeight_var,"evtWeight_var/F");
    myEvents->Branch("mcWeight",&mcWeight_var,"mcWeight_var/F");
    myEvents->Branch("pileupWeight",&pileupWeight_var,"pileupWeight_var/F");
    myEvents->Branch("vertexWeight",&vertexWeight_var,"vertexWeight_var/F");

    myEvents->Branch("qscale",&qscale,"qscale/F");
    myEvents->Branch("x1",&x1,"x1/F");
    myEvents->Branch("x2",&x2,"x2/F");
    myEvents->Branch("id1",&id1,"id1/I");
    myEvents->Branch("id2",&id2,"id2/I");





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
    myEvents->Branch("electron_charge",electron_charge,"electron_charge[nElectrons]/I");

    myEvents->Branch("nMuons", &nMuons,"nMuons/I");
    myEvents->Branch("muon_Px", muon_Px,"muon_Px[nMuons]/F");
    myEvents->Branch("muon_Py", muon_Py,"muon_Py[nMuons]/F");
    myEvents->Branch("muon_Pz", muon_Pz,"muon_Pz[nMuons]/F");
    myEvents->Branch("muon_E",  muon_E, "muon_E[nMuons]/F");
    myEvents->Branch("muon_charge",muon_charge,"muon_charge[nMuons]/I");


    myEvents->Branch("nJets", &nJets,"nJets/I");
    myEvents->Branch("jet_Px", jet_Px,"jet_Px[nJets]/F");
    myEvents->Branch("jet_Py", jet_Py,"jet_Py[nJets]/F");
    myEvents->Branch("jet_Pz", jet_Pz,"jet_Pz[nJets]/F");
    myEvents->Branch("jet_E", jet_E,"jet_E[nJets]/F");

    myEvents->Branch("met", &met,"met/F");
    myEvents->Branch("phi_met", &phi_met,"phi_met/F");
    myEvents->Branch("sumet", &sumet,"sumet/F");

    myEvents->Branch("met_hv", &met_hv,"met_hv/F");
    myEvents->Branch("phi_met_hv", &phi_met_hv,"phi_met_hv/F");
    myEvents->Branch("sumet_hv", &sumet_hv,"sumet_hv/F");

    //mc truth
    myEvents->Branch("nmcparticles",  &nmcparticles,   "nmcparticles/I");
    myEvents->Branch("mc_px",         mc_px,           "mc_px[nmcparticles]/F");
    myEvents->Branch("mc_py",         mc_py,           "mc_py[nmcparticles]/F");
    myEvents->Branch("mc_pz",         mc_pz,           "mc_pz[nmcparticles]/F");
    myEvents->Branch("mc_en",         mc_en,           "mc_en[nmcparticles]/F");
    myEvents->Branch("mc_id",         mc_id,           "mc_id[nmcparticles]/I");
    myEvents->Branch("mc_type",       mc_type,         "mc_type[nmcparticles]/I");
    myEvents->Branch("mc_origin",     mc_origin,       "mc_origin[nmcparticles]/I");


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
    SG::AuxElement::Accessor<char> isPassed("isPassed");
    SG::AuxElement::Accessor<char> isPassedJetEventClean("isPassedJetEventClean");
    SG::AuxElement::Accessor<char> isDalitz("isDalitz");

    // from h013, TST is diphoton vertex based
    SG::AuxElement::Accessor<float> met_TST("met_TST");
    SG::AuxElement::Accessor<float> met_sumet("sumet_TST");
    SG::AuxElement::Accessor<float> met_phi("phi_TST");

    SG::AuxElement::Accessor<float> met_hardVertexTST("met_hardVertexTST");
    SG::AuxElement::Accessor<float> sumet_hardVertexTST("sumet_hardVertexTST");
    SG::AuxElement::Accessor<float> phi_hardVertexTST("phi_hardVertexTST");

    SG::AuxElement::Accessor<int> truthType("truthType");
    SG::AuxElement::Accessor<int> truthOrigin("truthOrigin");



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
    if(isMC()) {
        isDalitz_var = isDalitz(*HGameventInfo) == 1 ? 1 : 0 ;
        if(isDalitz_var) return EL::StatusCode::SUCCESS;
    }


    RunNumber = runNumber( *eventInfo );
    LumiBlock = lumiBlock( *eventInfo );
    EventNumber = eventNumber( *eventInfo );
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


    //histoStore()->fillTH1F("mu", eventHandler()->mu(), weight());

    if(isMC()) {
        HgammaAnalysis::setSelectedObjects(&photons, &electrons, &muons, &jets);
        //lumiXsecWeight_var * evtWeight_var
        totWeight_var = HgammaAnalysis::weightFinal();

        lumiXsecWeight_var = HgammaAnalysis::lumiXsecWeight();
        // mcWeight_var * pileupWeight_var * vertexWeight_var * weigth from Photon SF
        evtWeight_var = HgammaAnalysis::weight();

        mcWeight_var = eventHandler()->mcWeight();
        pileupWeight_var = eventHandler()->pileupWeight();
        vertexWeight_var = eventHandler()->vertexWeight();

        const xAOD::TruthEventContainer *truthEvents = nullptr;
        if (event()->retrieve(truthEvents, "TruthEvents").isFailure())      HG::fatal("Can't access TruthEvents");
        xAOD::TruthEventContainer seltruth =*truthEvents;
        qscale  = seltruth[0]->auxdataConst<float>("Q")*seltruth[0]->auxdataConst<float>("Q");
        x1  = seltruth[0]->auxdataConst<float>("X1");
        x2  = seltruth[0]->auxdataConst<float>("X2");
        id1 = seltruth[0]->auxdataConst<int>("PDGID1");
        id2 = seltruth[0]->auxdataConst<int>("PDGID2");
    }





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

        double pt = photons[gn]->p4().Pt();
        double eta = photons[gn]->p4().Eta();
        double phi = photons[gn]->p4().Phi();
        double e = photons[gn]->p4().E();

        TLorentzVector pho;
        pho.Clear();
        pho.SetPtEtaPhiE(pt,eta,phi,e);

        photon_Px[nPhotons] = pho.Px();
        photon_Py[nPhotons] = pho.Py();
        photon_Pz[nPhotons] = pho.Pz();
        photon_E[nPhotons] =  pho.E();

        nPhotons++;
    }


    //
    // Electron
    //
    nElectrons=0;
    for(size_t gn=0; gn<electrons.size(); gn++) {

        double pt = electrons[gn]->p4().Pt();
        double eta = electrons[gn]->p4().Eta();
        double phi = electrons[gn]->p4().Phi();
        double e = electrons[gn]->p4().E();

        TLorentzVector pho;
        pho.Clear();
        pho.SetPtEtaPhiE(pt,eta,phi,e);


        electron_Px[nElectrons] = pho.Px();
        electron_Py[nElectrons] = pho.Py();
        electron_Pz[nElectrons] = pho.Pz();
        electron_E[nElectrons]  = pho.E();

        electron_charge[nElectrons] = electrons[gn]->charge();
        nElectrons++;
    }



    //
    // Muon
    //
    nMuons=0;
    for(size_t gn=0; gn<muons.size(); gn++) {

        double pt = muons[gn]->p4().Pt();
        double eta = muons[gn]->p4().Eta();
        double phi = muons[gn]->p4().Phi();
        double e = muons[gn]->p4().E();

        TLorentzVector pho;
        pho.Clear();
        pho.SetPtEtaPhiE(pt,eta,phi,e);

        muon_Px[nMuons] = pho.Px();
        muon_Py[nMuons] = pho.Py();
        muon_Pz[nMuons] = pho.Pz();
        muon_E[nMuons]  = pho.E();

        muon_charge[nMuons] = muons[gn]->charge();
        nMuons++;
    }


    //
    // Jets
    //
    nJets=0;
    for(size_t gn=0; gn<jets.size(); gn++) {

        double pt = jets[gn]->p4().Pt();
        double eta = jets[gn]->p4().Eta();
        double phi = jets[gn]->p4().Phi();
        double e = jets[gn]->p4().E();

        TLorentzVector pho;
        pho.Clear();
        pho.SetPtEtaPhiE(pt,eta,phi,e);

        jet_Px[nJets] = pho.Px();
        jet_Py[nJets] = pho.Py();
        jet_Pz[nJets] = pho.Pz();
        jet_E[nJets]  = pho.E();

        nJets++;
    }



    //
    // MET
    //
    met = met_TST(*HGameventInfo);
    sumet = met_sumet(*HGameventInfo);
    phi_met = met_phi(*HGameventInfo);

    std::string inputfileName = wk()->inputFileName();
    TString TagName = inputfileName;
    // diphoton vertex start to default metTST from h013
    if(!TagName.Contains("h012")) {
        met_hv = met_hardVertexTST(*HGameventInfo);
        sumet_hv = sumet_hardVertexTST(*HGameventInfo);
        phi_met_hv = phi_hardVertexTST(*HGameventInfo);
    }



    //
    // Generator-level information
    //
    nmcparticles = 0;
    if( isMC() ) {

        xAOD::TruthParticleContainer    truthPhotons = truthHandler()->getPhotons() ;
        xAOD::TruthParticleContainer    truthElectrons = truthHandler()->getElectrons();
        xAOD::TruthParticleContainer    truthMuons = truthHandler()->getMuons();
        xAOD::JetContainer              truthJets = truthHandler()->getJets();
        xAOD::MissingETContainer        truthMET = truthHandler()->getMissingET();


        for( xAOD::TruthParticle* truthpart : truthPhotons) {

            mc_px[nmcparticles] = truthpart->p4().Px();
            mc_py[nmcparticles] = truthpart->p4().Py();
            mc_pz[nmcparticles] = truthpart->p4().Pz();
            mc_en[nmcparticles] = truthpart->p4().E();
            mc_origin[nmcparticles] = truthOrigin( *truthpart );
            mc_type[nmcparticles] = truthType( *truthpart );
            mc_id[nmcparticles] = 22;
            nmcparticles++;

        }





        mc_px[nmcparticles] = truthMET["NonInt"]->mpx();
        mc_py[nmcparticles] = truthMET["NonInt"]->mpy();
        mc_pz[nmcparticles] = 0;
        mc_en[nmcparticles] = truthMET["NonInt"]->met();
        mc_origin[nmcparticles] = -1;
        mc_type[nmcparticles] = -1;
        mc_id[nmcparticles] = 0; //met pdgid assgined as 0
        nmcparticles++;


    }





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
