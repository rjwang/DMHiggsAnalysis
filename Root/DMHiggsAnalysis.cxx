#include "DMHiggsAnalysis/DMHiggsAnalysis.h"
#include "HGamAnalysisFramework/HgammaIncludes.h"

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

    //histoStore()->createTH1F("m_yy", 60, 110, 140);


    //Create a TTree
    //  TFile *outfile = wk()->getOutputFile("MxAOD");


    return EL::StatusCode::SUCCESS;
}


void DMHiggsAnalysis::declareVariables()
{

    myEvents->Branch("RunNumber",&RunNumber,"RunNumber/I");
    myEvents->Branch("LumiBlock",&LumiBlock,"LumiBlock/I");
    myEvents->Branch("EventNumber",&EventNumber,"EventNumber/I");
    myEvents->Branch("mcEventWeights",&mcEventWeights_var,"mcEventWeights_var/F");

    myEvents->Branch("mcID",&mcID_var,"mcID_var/I");
    myEvents->Branch("NPV",&NPV_var,"NPV_var/I");
    myEvents->Branch("mu",&mu_var,"mu_var/F");
    myEvents->Branch("isMC",&isMC_var,"isMC_var/I");
    myEvents->Branch("initWeight",&initWeight_var,"initWeight_var/F");
    myEvents->Branch("xsectionBRfilterEff",&xsecBrFilterEff_var,"xsecBrFilterEff_var/F");
    // myEvents->Branch("myy",&myy_var,"/I");
    // myEvents->Branch("phi_yy_met",&phi_yy_met_var,"/F");
    // myEvents->Branch("phi_yyj_met",&phi_yyj_met_var,"/F");
    // myEvents->Branch("phi_yyjj_met",&phi_yyjj_met_var,"/F");
    // myEvents->Branch("phi_jj_met",&phi_jj_met_var,"/F");
    // myEvents->Branch("phi_yy_jj",&phi_yy_jj_var,"/F");
    // myEvents->Branch("phi_y1_j",&phi_y1_j_var,"/F");
    // myEvents->Branch("phi_y2_j",&phi_y2_j_var,"/F");
    // myEvents->Branch("phi_y1_met",&phi_y1_met_var,"/F");
    // myEvents->Branch("phi_y2_met",&phi_y2_met_var,"/F");
    // myEvents->Branch("phi_y1_y2",&phi_y1_y2_var,"/F");
    // myEvents->Branch("phi_j1_j2",&phi_j1_j2_var,"/F");
    // myEvents->Branch("phi_allparticles_met",&phi_allparticles_met_var,"/F");
    myEvents->Branch("metSig",&metSig_var,"metSig_var/F");
    myEvents->Branch("passVertex",&passVertex_var,"passVertex_var/I");
    myEvents->Branch("passHiggsSelection",&passHiggsSelection_var,"passHiggsSelection_var/I");
    myEvents->Branch("passQualityCuts",&passQualityCuts_var,"passQualityCuts_var/I");


    myEvents->Branch("nPhotons", &nPhotons,"nPhotons/I");
    //  myEvents->Branch("photonAuthor", photonAuthor,"photonAuthor[nPhotons]/I");
    myEvents->Branch("photon_Px", photon_Px,"photon_Px[nPhotons]/F");
    myEvents->Branch("photon_Py", photon_Py,"photon_Py[nPhotons]/F");
    myEvents->Branch("photon_Pz", photon_Pz,"photon_Pz[nPhotons]/F");
    myEvents->Branch("photon_E",  photon_E, "photon_E[nPhotons]/F");
    myEvents->Branch("photon_Eps",photon_Eps,"photon_Eps[nPhotons]/F");
    myEvents->Branch("photon_E1",photon_E1,"photon_E1[nPhotons]/F");
    myEvents->Branch("photon_E2",photon_E2,"photon_E2[nPhotons]/F");
    myEvents->Branch("photon_E3",photon_E3,"photon_E3[nPhotons]/F");
    myEvents->Branch("photon_conversion",photon_conversionType,"photon_conversionType[nPhotons]/I");
    myEvents->Branch("photon_conversionRadius",photon_conversionRadius,"photon_conversionRadius[nPhotons]/F");
    myEvents->Branch("photon_isTight", photon_isTight,"photon_isTight[nPhotons]/I");
    myEvents->Branch("photon_isLoose", photons_isLoose,"photons_isLoose[nPhotons]/I");
    myEvents->Branch("photon_isLoosePrime2", photons_isLoosePrime2,"photons_isLoosePrime2[nPhotons]/I");
    myEvents->Branch("photon_isLoosePrime3", photons_isLoosePrime3,"photons_isLoosePrime3[nPhotons]/I");
    myEvents->Branch("photon_isLoosePrime4", photons_isLoosePrime4,"photons_isLoosePrime4[nPhotons]/I");
    myEvents->Branch("photon_isLoosePrime5", photons_isLoosePrime5,"photons_isLoosePrime5[nPhotons]/I");
    myEvents->Branch("photon_isIsoFixedCutTight",  photon_isIsoFixedCutTight, "photon_isIsoFixedCutTight[nPhotons]/I");
    myEvents->Branch("photon_isIsoFixedCutLoose",photon_isIsoFixedCutLoose,"photon_isIsoFixedCutLoose[nPhotons]/I");
    myEvents->Branch("photon_isIsoFixedCutLooseCaloOnly",photon_isIsoFixedCutLooseCaloOnly,"photon_isIsoFixedCutLooseCaloOnly[nPhotons]/I");
    myEvents->Branch("photon_isIsoFixedCutTightCaloOnly",photon_isIsoFixedCutTightCaloOnly,"photon_isIsoFixedCutTightCaloOnly[nPhotons]/I");
    myEvents->Branch("photon_topoCone20",photon_topoCone20,"photon_topoCone20[nPhotons]/F");
    myEvents->Branch("photon_topoCone40",photon_topoCone40,"photon_topoCone40[nPhotons]/F");
    myEvents->Branch("photon_Cone20",photon_Cone20,"photon_Cone20[nPhotons]/F");
    myEvents->Branch("photon_Cone40",photon_Cone40,"photon_Cone40[nPhotons]/F");

    myEvents->Branch("nElectrons", &nElectrons,"nElectrons/I");
    //  myEvents->Branch("electronAuthor", electronAuthor,"electronAuthor[nElectrons]/I");
    myEvents->Branch("electron_Px", electron_Px,"electron_Px[nElectrons]/F");
    myEvents->Branch("electron_Py", electron_Py,"electron_Py[nElectrons]/F");
    myEvents->Branch("electron_Pz", electron_Pz,"electron_Pz[nElectrons]/F");
    myEvents->Branch("electron_E",  electron_E, "electron_E[nElectrons]/F");
    myEvents->Branch("electron_Eps",electron_Eps,"electron_Eps[nElectrons]/F");
    myEvents->Branch("electron_E1",electron_E1,"electron_E1[nElectrons]/F");
    myEvents->Branch("electron_E2",electron_E2,"electron_E2[nElectrons]/F");
    myEvents->Branch("electron_E3",electron_E3,"electron_E3[nElectrons]/F");
    myEvents->Branch("electron_charge",electron_charge,"electron_charge[nElectrons]/F");
    myEvents->Branch("electron_isTight", electron_isTight,"electron_isTight[nElectrons]/I");
    myEvents->Branch("electron_isMedium", electron_isMedium,"electron_isMedium[nElectrons]/I");
    myEvents->Branch("electron_isIsoLoose", electron_isIsoLoose,"electron_isIsoLoose[nElectrons]/I");
    myEvents->Branch("electron_ptvarCone20",electron_ptvarCone20,"electron_ptvarCone20[nElectrons]/F");
    myEvents->Branch("electron_topoCone20",electron_topoCone20,"electron_topoCone20[nElectrons]/F");


    myEvents->Branch("nMuons", &nMuons,"nMuons/I");
    //  myEvents->Branch("muonAuthor", muonAuthor,"muonAuthor[nMuons]/I");
    myEvents->Branch("muon_Px", muon_Px,"muon_Px[nMuons]/F");
    myEvents->Branch("muon_Py", muon_Py,"muon_Py[nMuons]/F");
    myEvents->Branch("muon_Pz", muon_Pz,"muon_Pz[nMuons]/F");
    myEvents->Branch("muon_E",  muon_E, "muon_E[nMuons]/F");
    myEvents->Branch("muon_charge",muon_charge,"muon_charge[nMuons]/F");
    myEvents->Branch("muon_passIPcut", muon_passIPcut,"muon_passIPcut[nMuons]/I");
    // myEvents->Branch("muon_isTight", muon_isTight,"muon_isTight[nMuons]/I");
    // myEvents->Branch("muon_isMedium", muon_isMedium,"muon_isMedium[nMuons]/I");
    // myEvents->Branch("muon_isLoose", muon_isLoose,"muon_isLoose[nMuons]/I");
    myEvents->Branch("muon_isIsoGradientLoose", muon_isIsoGradientLoose,"muon_isIsoGradientLoose[nMuons]/I");
    myEvents->Branch("muon_isIsoGradient", muon_isIsoGradient,"muon_isIsoGradient[nMuons]/I");
    myEvents->Branch("muon_isIsoLoose", muon_isIsoLoose,"muon_isIsoLoose[nMuons]/I");
    myEvents->Branch("muon_ptvarCone20",muon_ptvarCone20,"muon_ptvarCone20[nMuons]/F");
    myEvents->Branch("muon_topoCone20",muon_topoCone20,"muon_topoCone20[nMuons]/F");


    myEvents->Branch("nJets", &nJets,"nJets/I");
    myEvents->Branch("jet_Px", jet_Px,"jet_Px[nJets]/F");
    myEvents->Branch("jet_Py", jet_Py,"jet_Py[nJets]/F");
    myEvents->Branch("jet_Pz", jet_Pz,"jet_Pz[nJets]/F");
    myEvents->Branch("jet_E", jet_E,"jet_E[nJets]/F");
    myEvents->Branch("jet_Jvt",  jet_Jvt, "jet_Jvt[nJets]/F");
    myEvents->Branch("jet_PassSelection",  jet_PassSelection, "jet_PassSelection[nJets]/I");

    myEvents->Branch("met", &met,"met/F");
    myEvents->Branch("sumet", &sumet,"sumet/F");
    myEvents->Branch("phi_met", &phi_met,"phi_met/F");


    myEvents->Branch("ntruthPhotons", &ntruthPhotons,"ntruthPhotons/I");
    myEvents->Branch("photonTruthPx", photonTruthPx,"photonTruthPx[ntruthPhotons]/F");
    myEvents->Branch("photonTruthPy", photonTruthPy,"photonTruthPy[ntruthPhotons]/F");
    myEvents->Branch("photonTruthPz", photonTruthPz,"photonTruthPz[ntruthPhotons]/F");
    myEvents->Branch("photonTruthE",  photonTruthE, "photonTruthE[ntruthPhotons]/F");
    myEvents->Branch("photonTruth_ptcone20",photonTruth_ptcone20,"photonTruth_ptcone20[ntruthPhotons]/F");
    myEvents->Branch("photonTruth_ptcone40",photonTruth_ptcone40,"photonTruth_ptcone40[ntruthPhotons]/F");
    myEvents->Branch("photonTruth_etcone20",photonTruth_etcone20,"photonTruth_etcone20[ntruthPhotons]/F");
    myEvents->Branch("photonTruth_etcone40",photonTruth_etcone40,"photonTruth_etcone40[ntruthPhotons]/F");
    myEvents->Branch("photonTruth_truthOrigin",photonTruth_truthOrigin,"photonTruth_truthOrigin[ntruthPhotons]/I");
    myEvents->Branch("photonTruth_truthType",photonTruth_truthType,"photonTruth_truthType[ntruthPhotons]/I");

    myEvents->Branch("ntruthElectrons", &ntruthElectrons,"ntruthElectrons/I");
    myEvents->Branch("electronTruthPx", electronTruthPx,"electronTruthPx[ntruthElectrons]/F");
    myEvents->Branch("electronTruthPy", electronTruthPy,"electronTruthPy[ntruthElectrons]/F");
    myEvents->Branch("electronTruthPz", electronTruthPz,"electronTruthPz[ntruthElectrons]/F");
    myEvents->Branch("electronTruthE",  electronTruthE, "electronTruthE[ntruthElectrons]/F");

    myEvents->Branch("ntruthMuons", &ntruthMuons,"ntruthMuons/I");
    myEvents->Branch("muonTruthPx", muonTruthPx,"muonTruthPx[ntruthMuons]/F");
    myEvents->Branch("muonTruthPy", muonTruthPy,"muonTruthPy[ntruthMuons]/F");
    myEvents->Branch("muonTruthPz", muonTruthPz,"muonTruthPz[ntruthMuons]/F");
    myEvents->Branch("muonTruthE",  muonTruthE, "muonTruthE[ntruthMuons]/F");

    myEvents->Branch("ntruthJets", &ntruthJets,"ntruthJets/I");
    myEvents->Branch("jetTruthPx", jetTruthPx,"jetTruthPx[ntruthJets]/F");
    myEvents->Branch("jetTruthPy", jetTruthPy,"jetTruthPy[ntruthJets]/F");
    myEvents->Branch("jetTruthPz", jetTruthPz,"jetTruthPz[ntruthJets]/F");
    myEvents->Branch("jetTruthE",  jetTruthE, "jetTruthE[ntruthJets]/F");

    myEvents->Branch("mpxTruthInt", &mpxTruthInt,"mpxTruthInt/F");
    myEvents->Branch("mpyTruthInt", &mpyTruthInt,"mpyTruthInt/F");
    myEvents->Branch("metTruthInt", &metTruthInt,"metTruthInt/F");
    myEvents->Branch("sumetTruthInt",  &sumetTruthInt, "sumetTruthInt/F");
    myEvents->Branch("mpxTruthNonInt", &mpxTruthNonInt,"mpxTruthNonInt/F");
    myEvents->Branch("mpyTruthNonInt", &mpyTruthNonInt,"mpyTruthNonInt/F");
    myEvents->Branch("metTruthNonInt", &metTruthNonInt,"metTruthNonInt/F");
    myEvents->Branch("sumetTruthNonInt",  &sumetTruthNonInt, "sumetTruthNonInt/F");
}


void DMHiggsAnalysis::clearVectors()
{


    for( int iparticle = 0 ; iparticle < MAXPARTICLES ; ++iparticle) {

        //    photonAuthor[iparticle] =  - 9999;
        photon_Px[iparticle] =  - 9999;
        photon_E[iparticle] =  - 9999;
        photon_Py[iparticle] =  - 9999;
        photon_Pz[iparticle] =  - 9999;
        photon_Eps[iparticle] =  - 9999;
        photon_E1[iparticle] =  - 9999;
        photon_E2[iparticle] =  - 9999;
        photon_E3[iparticle] =  - 9999;

        photon_conversionType[iparticle] =  - 9999;
        photon_isTight[iparticle] =  - 9999;
        photons_isLoose[iparticle] =  - 9999;
        photons_isLoosePrime2[iparticle] =  - 9999;
        photons_isLoosePrime3[iparticle] =  - 9999;
        photons_isLoosePrime4[iparticle] =  - 9999;
        photons_isLoosePrime5[iparticle] =  - 9999;

        photon_isIsoFixedCutTight[iparticle] =  - 9999;
        photon_isIsoFixedCutTightCaloOnly[iparticle] =  - 9999;
        photon_isIsoFixedCutLooseCaloOnly[iparticle] =  - 9999;
        photon_isIsoFixedCutLoose[iparticle] =  - 9999;
        photon_Cone20[iparticle] =  - 9999;
        photon_Cone40[iparticle] =  - 9999;
        photon_topoCone20[iparticle] =  - 9999;
        photon_topoCone40[iparticle] =  - 9999;

        //electronAuthor[iparticle] =  - 9999;
        electron_Px[iparticle] =  - 9999;
        electron_E[iparticle] =  - 9999;
        electron_Py[iparticle] =  - 9999;
        electron_Pz[iparticle] =  - 9999;
        electron_Eps[iparticle] =  - 9999;
        electron_E1[iparticle] =  - 9999;
        electron_E2[iparticle] =  - 9999;
        electron_E3[iparticle] =  - 9999;

        electron_isTight[iparticle] =  - 9999;
        electron_isMedium[iparticle] =  - 9999;
        electron_isIsoLoose[iparticle] =  - 9999;
        electron_topoCone20[iparticle] =  - 9999;
        electron_ptvarCone20[iparticle] =  - 9999;

        //muonAuthor[iparticle] =  - 9999;
        muon_Px[iparticle] =  - 9999;
        muon_E[iparticle] =  - 9999;
        muon_Py[iparticle] =  - 9999;
        muon_Pz[iparticle] =  - 9999;

        muon_passIPcut[iparticle] =  - 9999;
        muon_topoCone20[iparticle] =  - 9999;
        muon_ptvarCone20[iparticle] =  - 9999;
        // muon_isLoose[iparticle] =  - 9999;
        // muon_isTight[iparticle] =  - 9999;
        // muon_isMedium[iparticle] =  - 9999;
        muon_isIsoGradientLoose[iparticle] =  - 9999;
        muon_isIsoGradient[iparticle] =  - 9999;
        muon_isIsoLoose[iparticle] =  - 9999;

        //jetAuthor[iparticle] =  - 9999;
        jet_Px[iparticle] =  - 9999;
        jet_Py[iparticle] =  - 9999;
        jet_Pz[iparticle] =  - 9999;
        jet_E[iparticle] =  - 9999;
        jet_Jvt[iparticle] =  - 9999;
        jet_PassSelection[iparticle] =  - 9999;



        photonTruthPx[iparticle] =  - 9999;
        photonTruthE[iparticle] =  - 9999;
        photonTruthPy[iparticle] =  - 9999;
        photonTruthPz[iparticle] =  - 9999;
        photonTruth_etcone20[iparticle] =  - 9999;
        photonTruth_etcone40[iparticle] =  - 9999;
        photonTruth_ptcone20[iparticle] =  - 9999;
        photonTruth_ptcone40[iparticle] =  - 9999;
        photonTruth_truthOrigin[iparticle] =  - 9999;
        photonTruth_truthType[iparticle] =  - 9999;


        electronTruthPx[iparticle] =  - 9999;
        electronTruthE[iparticle] =  - 9999;
        electronTruthPy[iparticle] =  - 9999;
        electronTruthPz[iparticle] =  - 9999;

        muonTruthPx[iparticle] =  - 9999;
        muonTruthE[iparticle] =  - 9999;
        muonTruthPy[iparticle] =  - 9999;
        muonTruthPz[iparticle] =  - 9999;


        jetTruthPx[iparticle] =  - 9999;
        jetTruthE[iparticle] =  - 9999;
        jetTruthPy[iparticle] =  - 9999;
        jetTruthPz[iparticle] =  - 9999;

        mpxTruthInt =  - 9999;
        mpyTruthInt =  - 9999;
        metTruthInt =  - 9999;
        sumetTruthInt =  - 9999;
        mpxTruthNonInt =  - 9999;
        mpyTruthNonInt =  - 9999;
        metTruthNonInt =  - 9999;
        sumetTruthNonInt =  - 9999;
    }

}

EL::StatusCode DMHiggsAnalysis::initialize()
{


    HgammaAnalysis::initialize();
    std::string inputfileName = wk()->inputFileName();
    currentfilename = inputfileName;
    inputfileName.replace(inputfileName.find(".MxAOD") , -1, "") ;
    inputfileName.append(".NTuple.root");


    m_outputFile = TFile::Open(inputfileName.c_str(),"RECREATE");
    myEvents = new TTree("DMHiggsAnalysis","DMHiggsAnalysis");
    declareVariables();

    m_muonTightSelectionTool = new CP::MuonSelectionTool("MuonTightSelectionTool");
    m_muonMediumSelectionTool = new CP::MuonSelectionTool("MuonMediumSelectionTool");
    m_muonLooseSelectionTool = new CP::MuonSelectionTool("MuonLooseSelectionTool");

    CP_CHECK( "initialize()" , m_muonTightSelectionTool->setProperty("MuQuality", int(xAOD::Muon::Tight)) );
    CP_CHECK( "initialize()" , m_muonMediumSelectionTool->setProperty("MuQuality", int(xAOD::Muon::Medium)) );
    CP_CHECK( "initialize()" , m_muonLooseSelectionTool->setProperty("MuQuality", int(xAOD::Muon::Loose)) );

    if( m_muonTightSelectionTool->initialize().isFailure() )
        HG::fatal("Couldn't initalize MuonSelectionTool for tight WP. Exiting");
    if( m_muonMediumSelectionTool->initialize().isFailure() )
        HG::fatal("Couldn't initalize MuonSelectionTool for medium WP. Exiting");
    if( m_muonLooseSelectionTool->initialize().isFailure() )
        HG::fatal("Couldn't initalize MuonSelectionTool for loose WP. Exiting");




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

    HgammaAnalysis::execute();


    std::string inputfileName = wk()->inputFileName();

    if( currentfilename != inputfileName ) {

        currentfilename = inputfileName;

        inputfileName.replace(inputfileName.find(".MxAOD") , -1, "") ;
        inputfileName.append(".NTuple.root");

        m_outputFile->cd();

        myEvents->Write();
        m_outputFile->Close();

        m_outputFile = TFile::Open(inputfileName.c_str(),"RECREATE");
        myEvents = new TTree("DMHiggsAnalysis","DMHiggsAnalysis");

        declareVariables();

        inputfileName = wk()->inputFileName();


    }


    SG::AuxElement::Accessor<unsigned int> runNumber("runNumber");
    SG::AuxElement::Accessor<unsigned int> lumiBlock("lumiBlock");
    SG::AuxElement::Accessor<unsigned long long> eventNumber("eventNumber");
    SG::AuxElement::Accessor< std::vector<float> > mcEventWeights("mcEventWeights");

    SG::AuxElement::Accessor<int> NPV("numberOfPrimaryVertices");
    SG::AuxElement::Accessor<float> mu("mu");
    SG::AuxElement::Accessor<float> initWeight("weightInitial");
    SG::AuxElement::Accessor<float> crossSectionBRflterEff("crossSectionBRfilterEff");
    //  SG::AuxElement::Accessor<float> TotalWeight("TotalWeight");
    SG::AuxElement::Accessor<float> myy("m_yy");
    SG::AuxElement::Accessor<float> met_sumet("sumet_TST");
    SG::AuxElement::Accessor<float> met_phi("phi_TST");
    SG::AuxElement::Accessor<float> met_TST("met_TST");
    SG::AuxElement::Accessor<float> phi_allparticles_met("DeltaPhi_allPart_MET");
    SG::AuxElement::Accessor<float> photon_ePS("E0_raw");
    SG::AuxElement::Accessor<float> photon_e1("E1_raw");
    SG::AuxElement::Accessor<float> photon_e2("E2_raw");
    SG::AuxElement::Accessor<float> photon_e3("E3_raw");
    SG::AuxElement::Accessor<int> photon_converted("conversionType");
    SG::AuxElement::Accessor<float> photon_convertedRadius("conversionRadius");
    SG::AuxElement::Accessor<float> hardestVertex("hardestVertexZ");
    SG::AuxElement::Accessor<float> selectedVertex("selectedVertexZ");
    SG::AuxElement::Accessor<char> passHiggsSelection("isPassed");
    SG::AuxElement::Accessor<char> passQualityCuts("isPassedBasic");
    SG::AuxElement::Accessor<char> passDoubleHiggsSelection("isPassedLowHighMyy");

    SG::AuxElement::Accessor<float> ptvarCone20("ptvarcone20");
    SG::AuxElement::Accessor<float> topoetCone20("topoetcone20");
    SG::AuxElement::Accessor<float> topoetCone40("topoetcone40");
    SG::AuxElement::Accessor<float> etCone20("topoetcone20");
    SG::AuxElement::Accessor<float> etCone40("topoetcone40");
    SG::AuxElement::Accessor<float> ptCone20("ptcone20");
    SG::AuxElement::Accessor<float> ptCone40("ptcone40");
    SG::AuxElement::Accessor<int> truthType("truthType");
    SG::AuxElement::Accessor<int> truthOrigin("truthOrigin");
    SG::AuxElement::Accessor<char> isisoFixedCutTight("isIsoFixedCutTight");
    SG::AuxElement::Accessor<char> isisoFixedCutLoose("isIsoFixedCutLoose");
    SG::AuxElement::Accessor<char> isisoFixedCutTightCaloOnly("isIsoFixedCutTightCaloOnly");
    SG::AuxElement::Accessor<char> isisoFixedCutLooseCaloOnly("isIsoFixedCutLooseCaloOnly");
    SG::AuxElement::Accessor<char> isTight("isTight");
    SG::AuxElement::Accessor<char> isMedium("isMedium");
    SG::AuxElement::Accessor<char> isLoose("isLoose");
    SG::AuxElement::Accessor<char> passIPcut("passIPCut");

    SG::AuxElement::Accessor<float> Jvt("Jvt");



    // photons


    const xAOD::EventInfo* HGameventInfo = 0 ;

    if(event()->retrieve(HGameventInfo,"HGamEventInfo").isFailure() )
        HG::fatal("Cannot retrieve event Info .");



    const xAOD::EventInfo* eventInfo = 0 ;

    if(event()->retrieve(eventInfo,"EventInfo").isFailure() )
        HG::fatal("Cannot retrieve event Info .");


    RunNumber = runNumber( *eventInfo );
    LumiBlock = lumiBlock( *eventInfo );
    EventNumber = eventNumber( *eventInfo );

    TLorentzVector etmissVector;

    mcEventWeights_var = isMC() ?  mcEventWeights( *eventInfo )[0] : 1.;
    xsecBrFilterEff_var = crossSectionBRflterEff.isAvailable( *HGameventInfo ) ?  crossSectionBRflterEff( *HGameventInfo )  : 1.;
    mcID_var = isMC() ? eventInfo->mcChannelNumber() : -999;
    NPV_var = NPV(*HGameventInfo);
    mu_var = mu(*HGameventInfo);
    isMC_var = isMC();
    initWeight_var = initWeight(*HGameventInfo);
    //myy_var = myy(*HGameventInfo);
    metSig_var = met_TST(*HGameventInfo)/met_sumet(*HGameventInfo);
    passVertex_var = fabs(hardestVertex(*HGameventInfo) - selectedVertex(*HGameventInfo)) < 0.3 ? 1 : 0;
    passHiggsSelection_var = passHiggsSelection(*HGameventInfo) == 1 ?  1 : 0;
    passQualityCuts_var = passQualityCuts(*HGameventInfo) == 1 ? 1 : 0 ;



    std::string cutFlowName ;

    cutFlowName = "CutFlow_" + inputfileName;
    cutFlowName.replace(cutFlowName.find(".MxAOD") , -1, "") ;
    cutFlowName.append("_weighted");

    if( m_eventCounter == 1 ) m_histCutFlow[inputfileName] = (TH1F*) HG::getHistogramFromFile(cutFlowName,inputfileName);

    xAOD::PhotonContainer photons = photonHandler()->getCorrectedContainer() ;
    xAOD::ElectronContainer electrons = electronHandler()->getCorrectedContainer() ;
    xAOD::MuonContainer muons = muonHandler()->getCorrectedContainer() ;
    xAOD::JetContainer jets = jetHandler()->getCorrectedContainer() ;

    //move this overlap checking to selection level
    //overlapHandler()->removeOverlap(photons,jets,electrons, muons);


    //
    // Photon
    //
    nPhotons=0;
    for(size_t gn=0; gn<photons.size(); gn++) {
        //    photonAuthor[nPhotons] = photons[gn]->author();
        photon_Px[nPhotons] = photons[gn]->p4().Px();
        photon_Py[nPhotons] = photons[gn]->p4().Py();
        photon_Pz[nPhotons] = photons[gn]->p4().Pz();
        photon_E[nPhotons] = photons[gn]->p4().E();
        photon_Eps[nPhotons] = photon_ePS( *photons[gn] );
        photon_E1[nPhotons] = photon_e1( *photons[gn] );
        photon_E2[nPhotons] = photon_e2( *photons[gn] );
        photon_E3[nPhotons] = photon_e3( *photons[gn] );
        photon_conversionType[nPhotons] = photon_converted( *photons[gn] );
        photon_conversionRadius[nPhotons] = photon_convertedRadius( *photons[gn] );
        photon_topoCone20[nPhotons] = topoetCone20( *photons[gn] );
        photon_topoCone40[nPhotons] = topoetCone40( *photons[gn] );
        photon_Cone20[nPhotons] = ptCone20( *photons[gn] );
        photon_Cone40[nPhotons] = ptCone40( *photons[gn] );
        photon_isIsoFixedCutTight[nPhotons] = isisoFixedCutTight( *photons[gn] ) == 1 ? 1 : 0;
        photon_isIsoFixedCutTightCaloOnly[nPhotons] = isisoFixedCutTightCaloOnly( *photons[gn] ) == 1 ? 1 : 0;
        photon_isIsoFixedCutLooseCaloOnly[nPhotons] = isisoFixedCutLooseCaloOnly( *photons[gn] ) == 1 ? 1 : 0;
        photon_isIsoFixedCutLoose[nPhotons] = isisoFixedCutLoose( *photons[gn] ) == 1 ? 1 : 0;
        photon_isTight[nPhotons] = isTight( *photons[gn] ) == 1 ? 1 : 0;
        photons_isLoose[nPhotons] = photonHandler()->passPIDCut( photons[gn] , egammaPID::IsEMLoose  ) ? 1 : 0;
        photons_isLoosePrime2[nPhotons] = photonHandler()->passLoosePrime( photons[gn] , 2 ) == 1 ? 1 : 0;
        photons_isLoosePrime3[nPhotons] = photonHandler()->passLoosePrime( photons[gn] , 3 ) == 1 ? 1 : 0;
        photons_isLoosePrime4[nPhotons] = photonHandler()->passLoosePrime( photons[gn] , 4 ) == 1 ? 1 : 0;
        photons_isLoosePrime5[nPhotons] = photonHandler()->passLoosePrime( photons[gn] , 5 ) == 1 ? 1 : 0;
        nPhotons++;
    }


    //
    // Electron
    //
    nElectrons=0;
    for(size_t gn=0; gn<electrons.size(); gn++) {
        //    electronAuthor[nElectrons] = electrons[gn]->author();
        electron_Px[nElectrons] = electrons[gn]->p4().Px();
        electron_Py[nElectrons] = electrons[gn]->p4().Py();
        electron_Pz[nElectrons] = electrons[gn]->p4().Pz();
        electron_E[nElectrons] = electrons[gn]->e();
        electron_Eps[nElectrons] = electrons[gn]->caloCluster() != nullptr ? electrons[gn]->caloCluster()->energyBE(0) : 0. ;
        electron_E1[nElectrons] = electrons[gn]->caloCluster() != nullptr ? electrons[gn]->caloCluster()->energyBE(1) : 0. ;
        electron_E2[nElectrons] = electrons[gn]->caloCluster() != nullptr ? electrons[gn]->caloCluster()->energyBE(2) : 0. ;
        electron_E3[nElectrons] = electrons[gn]->caloCluster() != nullptr ? electrons[gn]->caloCluster()->energyBE(3) : 0. ;

        electron_topoCone20[nElectrons] = topoetCone20( *electrons[gn] );
        electron_ptvarCone20[nElectrons] = ptvarCone20( *electrons[gn] );
        electron_isTight[nElectrons] = isTight( *electrons[gn] ) == 1 ? 1 : 0;
        electron_isMedium[nElectrons] = electronHandler()->passPIDCut( electrons[gn] , "Medium" )  ?  1 : 0;
        electron_isIsoLoose[nElectrons] = electronHandler()->passIsoCut( electrons[gn] , HG::Iso::Loose )  ? 1 : 0;
        electron_isMedium[nElectrons] = electronHandler()->passPIDCut( electrons[gn] , "Medium" )  ? 1 : 0;
        electron_charge[nElectrons] = electrons[gn]->charge() ;
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

        muon_topoCone20[nMuons] = topoetCone20( *muons[gn] );
        muon_ptvarCone20[nMuons] = ptvarCone20( *muons[gn] );
        muon_charge[nMuons] = muons[gn]->charge() ;
        muon_passIPcut[nMuons] = passIPcut( *muons[gn] );
        muon_isIsoGradientLoose[nMuons] = muonHandler()->passIsoCut( muons[gn], HG::Iso::GradientLoose ) ? 1 : 0;
        muon_isIsoGradient[nMuons] = muonHandler()->passIsoCut( muons[gn], HG::Iso::Gradient ) ? 1 : 0;
        muon_isIsoLoose[nMuons] = muonHandler()->passIsoCut( muons[gn], HG::Iso::Loose ) ? 1 : 0;

        // muon_isLoose[nMuons] = m_muonLooseSelectionTool->accept( muons[gn] ) ? 1 : 0;
        // muon_isMedium[nMuons] = m_muonMediumSelectionTool->accept( muons[gn] ) ? 1 : 0;
        // muon_isTight[nMuons] = m_muonTightSelectionTool->accept( muons[gn] ) ? 1 : 0;


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
        jet_Jvt[nJets] = Jvt( *jets[gn] );
        jet_PassSelection[nJets] = ( jetHandler()->passPtEtaCuts( jets[gn]) && jetHandler()->passJVFCut( jets[gn]) ) && jetHandler()->passJVTCut( jets[gn] ) ? 1 : 0;

        nJets++;
    }




    //
    // MET
    //
    met = met_TST(*HGameventInfo);
    sumet = met_sumet(*HGameventInfo);
    phi_met = met_phi(*HGameventInfo);


    ntruthPhotons = 0 ;
    ntruthElectrons = 0 ;
    ntruthMuons = 0 ;
    ntruthJets = 0 ;

    //
    // Generator-level information
    //
    //HG::TruthParticles* 		truthParticles = truthHandler()->getTruthParticles();

    if( isMC() ) {
        xAOD::TruthParticleContainer 	truthPhotons = truthHandler()->getPhotons() ;
        xAOD::TruthParticleContainer 	truthElectrons = truthHandler()->getElectrons();
        xAOD::TruthParticleContainer 	truthMuons = truthHandler()->getMuons();
        xAOD::JetContainer 		truthJets = truthHandler()->getJets();
        xAOD::MissingETContainer 	truthMET = truthHandler()->getMissingET();


        ntruthPhotons = 0 ;

        for( xAOD::TruthParticle* truthpart : truthPhotons) {

            photonTruthPx[ntruthPhotons] = truthpart->p4().Px();
            photonTruthPy[ntruthPhotons] = truthpart->p4().Py();
            photonTruthPz[ntruthPhotons] = truthpart->p4().Pz();
            photonTruthE[ntruthPhotons] = truthpart->p4().E();
            photonTruth_ptcone20[ntruthPhotons] = ptCone20( *truthpart );
            photonTruth_ptcone40[ntruthPhotons] = ptCone40( *truthpart );
            photonTruth_etcone20[ntruthPhotons] = etCone20( *truthpart );
            photonTruth_etcone40[ntruthPhotons] = etCone40( *truthpart );
            photonTruth_truthOrigin[ntruthPhotons] = truthOrigin( *truthpart );
            photonTruth_truthType[ntruthPhotons] = truthType( *truthpart );


            ++ntruthPhotons;
        }


        ntruthElectrons = 0 ;
        for( xAOD::TruthParticle* truthpart : truthElectrons) {
            electronTruthPx[ntruthElectrons] = truthpart->p4().Px();
            electronTruthPy[ntruthElectrons] = truthpart->p4().Py();
            electronTruthPz[ntruthElectrons] = truthpart->p4().Pz();
            electronTruthE[ntruthElectrons] = truthpart->p4().E();
            ++ntruthElectrons;
        }




        ntruthMuons = 0 ;
        for( xAOD::TruthParticle* truthpart : truthMuons) {
            muonTruthPx[ntruthMuons] = truthpart->p4().Px();
            muonTruthPy[ntruthMuons] = truthpart->p4().Py();
            muonTruthPz[ntruthMuons] = truthpart->p4().Pz();
            muonTruthE[ntruthMuons] = truthpart->p4().E();
            ++ntruthMuons;
        }


        ntruthJets = 0 ;
        for( xAOD::Jet* truthpart : truthJets) {

            jetTruthPx[ntruthJets] = truthpart->p4().Px();
            jetTruthPy[ntruthJets] = truthpart->p4().Py();
            jetTruthPz[ntruthJets] = truthpart->p4().Pz();
            jetTruthE[ntruthJets] = truthpart->p4().E();
            ++ntruthJets;
        }



        // mpxTruthInt = ((*truthMET["Int"])+(*truthMET["IntMuons"])).mpx();
        // mpyTruthInt = ((*truthMET["Int"])+(*truthMET["IntMuons"])).mpy();
        // metTruthInt = ((*truthMET["Int"])+(*truthMET["IntMuons"])).met();
        // sumetTruthInt = ((*truthMET["Int"])+(*truthMET["IntMuons"])).sumet();

        mpxTruthInt = truthMET["Int"]->mpx();
        mpyTruthInt = truthMET["Int"]->mpy();
        metTruthInt = truthMET["Int"]->met();
        sumetTruthInt = truthMET["Int"]->sumet();
        mpxTruthNonInt = truthMET["NonInt"]->mpx();
        mpyTruthNonInt = truthMET["NonInt"]->mpy();
        metTruthNonInt = truthMET["NonInt"]->met();
        sumetTruthNonInt = truthMET["NonInt"]->sumet();

    }
    //  for( xAOD::TruthParticle* )




    myEvents->Fill();
    return EL::StatusCode::SUCCESS;
}



EL::StatusCode DMHiggsAnalysis::finalize()
{


    m_outputFile->cd();
    myEvents->Write();
    //  if( m_histCutFlow ) m_histCutFlow[it.first]->Write();
    m_outputFile->Close();



    HgammaAnalysis::finalize();

    return EL::StatusCode::SUCCESS;

}
