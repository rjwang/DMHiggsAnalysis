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

  histoStore()->createTH1F("m_yy", 60, 110, 140);
  
  
  //Create a TTree
  //  TFile *outfile = wk()->getOutputFile("MxAOD");
  
  
    return EL::StatusCode::SUCCESS;
}
void DMHiggsAnalysis::declareVariables(){

  myEvents->Branch("RunNumber",&RunNumber,"RunNumber/I"); 
  myEvents->Branch("LumiBlock",&LumiBlock,"LumiBlock/I"); 
  myEvents->Branch("EventNumber",&EventNumber,"EventNumber/I"); 
  myEvents->Branch("mcEventWeights",&mcEventWeights_var,"mcEventWeights_var/F");

  myEvents->Branch("mcID",&mcID_var,"mcID_var/I"); 
  myEvents->Branch("NPV",&NPV_var,"NPV_var/I"); 
  myEvents->Branch("mu",&mu_var,"mu_var/I"); 
  myEvents->Branch("isMC",&isMC_var,"isMC_var/I"); 
  myEvents->Branch("initWeight",&initWeight_var,"initWeight_var/I"); 
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
  myEvents->Branch("photonPt", photonPt,"photonPt[nPhotons]/F");
  myEvents->Branch("photonEta", photonEta,"photonEta[nPhotons]/F");
  myEvents->Branch("photonPhi", photonPhi,"photonPhi[nPhotons]/F");
  myEvents->Branch("photonE",  photonE, "photonE[nPhotons]/F");
  myEvents->Branch("photons_Eps",photons_Eps,"photon_Eps[nPhotons]/F");
  myEvents->Branch("photons_E1",photons_E1,"photons_E1[nPhotons]/F");
  myEvents->Branch("photons_E2",photons_E2,"photons_E2[nPhotons]/F");
  myEvents->Branch("photons_E3",photons_E3,"photons_E3[nPhotons]/F");
  myEvents->Branch("photons_conversion",photons_conversionType,"photons_conversionType[nPhotons]/I");
  myEvents->Branch("photons_conversionRadius",photons_conversionRadius,"photons_conversionRadius[nPhotons]/F");
  myEvents->Branch("photon_isTight", photons_isTight,"photons_isTight[nPhotons]/F");
  myEvents->Branch("photon_isFixedCutTight",  photons_isFixedCutTight, "photons_isFixedCutTight[nPhotons]/F");
  myEvents->Branch("photons_isFixedCutCaloOnly",photons_isFixedCutLoose,"photons_isFixedCutCaloOnly[nPhotons]/F");
  myEvents->Branch("photons_isFixedCutLoose",photons_isFixedCutLoose,"photons_isFixedCutLoose[nPhotons]/F");
  myEvents->Branch("photons_topoCone20",photons_topoCone20,"photons_topoCone20[nPhotons]/F");
  myEvents->Branch("photons_topoCone40",photons_topoCone40,"photons_topoCone40[nPhotons]/F");
  myEvents->Branch("photons_Cone20",photons_Cone20,"photons_Cone20[nPhotons]/F");
  myEvents->Branch("photons_Cone40",photons_Cone40,"photons_Cone40[nPhotons]/F");

  myEvents->Branch("nElectrons", &nElectrons,"nElectrons/I");
  //  myEvents->Branch("electronAuthor", electronAuthor,"electronAuthor[nElectrons]/I");
  myEvents->Branch("electronPt", electronPt,"electronPt[nElectrons]/F");
  myEvents->Branch("electronEta", electronEta,"electronEta[nElectrons]/F");
  myEvents->Branch("electronPhi", electronPhi,"electronPhi[nElectrons]/F");
  myEvents->Branch("electronE",  electronE, "electronE[nElectrons]/F");
  myEvents->Branch("electrons_Eps",electrons_Eps,"electron_Eps[nElectrons]/F");
  myEvents->Branch("electrons_E1",electrons_E1,"electrons_E1[nElectrons]/F");
  myEvents->Branch("electrons_E2",electrons_E2,"electrons_E2[nElectrons]/F");
  myEvents->Branch("electrons_E3",electrons_E3,"electrons_E3[nElectrons]/F");
  myEvents->Branch("electrons_charge",electrons_charge,"electrons_charge[nElectrons]/F");
  myEvents->Branch("electron_isTight", electrons_isTight,"electrons_isTight[nElectrons]/I");
  myEvents->Branch("electrons_ptvarCone20",electrons_ptvarCone20,"electrons_ptvarCone20[nElectrons]/F");
  myEvents->Branch("electrons_topoCone20",electrons_topoCone20,"electrons_topoCone20[nElectrons]/F");


  myEvents->Branch("nMuons", &nMuons,"nMuons/I");
  //  myEvents->Branch("muonAuthor", muonAuthor,"muonAuthor[nMuons]/I");
  myEvents->Branch("muonPt", muonPt,"muonPt[nMuons]/F");
  myEvents->Branch("muonEta", muonEta,"muonEta[nMuons]/F");
  myEvents->Branch("muonPhi", muonPhi,"muonPhi[nMuons]/F");
  myEvents->Branch("muonE",  muonE, "muonE[nMuons]/F");
  myEvents->Branch("muons_charge",muons_charge,"muons_charge[nMuons]/F");
  myEvents->Branch("muon_passIPcut", muons_passIPcut,"muons_passIPcut[nMuons]/I");
  myEvents->Branch("muons_ptvarCone20",muons_ptvarCone20,"muons_ptvarCone20[nMuons]/F");
  myEvents->Branch("muons_topoCone20",muons_topoCone20,"muons_topoCone20[nMuons]/F");


  myEvents->Branch("nJets", &nJets,"nJets/I");
  myEvents->Branch("jetPt", jetPt,"jetPt[nJets]/F");
  myEvents->Branch("jetEta", jetEta,"jetEta[nJets]/F");
  myEvents->Branch("jetPhi", jetPhi,"jetPhi[nJets]/F");
  myEvents->Branch("jetJvt",  jetJvt, "jetJvt[nJets]/F");
  myEvents->Branch("jetPassSelection",  jetPassSelection, "jetPassSelection[nJets]/I");

  myEvents->Branch("met", &met,"met/F");
  myEvents->Branch("sumet", &sumet,"sumet/F");
  myEvents->Branch("phi_met", &phi_met,"phi_met/F");
}

void DMHiggsAnalysis::clearVectors(){


  for( int iparticle = 0 ; iparticle < MAXPARTICLES ; ++iparticle){

    //    photonAuthor[iparticle] =  - 9999;
    photonPt[iparticle] =  - 9999;
    photonE[iparticle] =  - 9999;
    photonEta[iparticle] =  - 9999;
    photonPhi[iparticle] =  - 9999;
    photons_Eps[iparticle] =  - 9999;
    photons_E1[iparticle] =  - 9999;
    photons_E2[iparticle] =  - 9999;
    photons_E3[iparticle] =  - 9999;

    photons_conversionType[iparticle] =  - 9999;
    photons_isTight[iparticle] =  - 9999;
    //    photons_isLoose[iparticle] =  - 9999;
    //    photons_isLoosePrime[iparticle] =  - 9999;

    photons_isFixedCutTight[iparticle] =  - 9999;
    photons_isFixedCutTightCaloOnly[iparticle] =  - 9999;
    photons_isFixedCutLoose[iparticle] =  - 9999;
    photons_Cone20[iparticle] =  - 9999;
    photons_Cone40[iparticle] =  - 9999;
    photons_topoCone20[iparticle] =  - 9999;
    photons_topoCone40[iparticle] =  - 9999;

    //electronAuthor[iparticle] =  - 9999;
    electronPt[iparticle] =  - 9999;
    electronE[iparticle] =  - 9999;
    electronEta[iparticle] =  - 9999;
    electronPhi[iparticle] =  - 9999;
    electrons_Eps[iparticle] =  - 9999;
    electrons_E1[iparticle] =  - 9999;
    electrons_E2[iparticle] =  - 9999;
    electrons_E3[iparticle] =  - 9999;

    electrons_isTight[iparticle] =  - 9999;
    electrons_topoCone20[iparticle] =  - 9999;
    electrons_ptvarCone20[iparticle] =  - 9999;

    //muonAuthor[iparticle] =  - 9999;
    muonPt[iparticle] =  - 9999;
    muonE[iparticle] =  - 9999;
    muonEta[iparticle] =  - 9999;
    muonPhi[iparticle] =  - 9999;

    muons_passIPcut[iparticle] =  - 9999;
    muons_topoCone20[iparticle] =  - 9999;
    muons_ptvarCone20[iparticle] =  - 9999;

    //jetAuthor[iparticle] =  - 9999;
    jetPt[iparticle] =  - 9999;
    jetJvt[iparticle] =  - 9999;
    jetEta[iparticle] =  - 9999;
    jetPhi[iparticle] =  - 9999;
    jetPassSelection[iparticle] =  - 9999;
  }

}

EL::StatusCode DMHiggsAnalysis::initialize()
{


  HgammaAnalysis::initialize();
  outFile = TFile::Open("NTuple.root","RECREATE");



  myEvents = new TTree("DMHiggsAnalysis","DMHiggsAnalysis");

  wk()->addOutput(outFile);


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

  HgammaAnalysis::execute();


  SG::AuxElement::Accessor<unsigned int> runNumber("runNumber");
  SG::AuxElement::Accessor<unsigned int> lumiBlock("lumiBlock");
  SG::AuxElement::Accessor<unsigned long long> eventNumber("eventNumber");
  SG::AuxElement::Accessor< std::vector<float> > mcEventWeights("mcEventWeights");

  SG::AuxElement::Accessor<int> NPV("numberOfPrimaryVertices");
  SG::AuxElement::Accessor<float> mu("mu");
  SG::AuxElement::Accessor<float> initWeight("weightInitial");
  SG::AuxElement::Accessor<float> crossSectionBRflterEff("crossSectionBRflterEff");

  //  SG::AuxElement::Accessor<float> TotalWeight("TotalWeight");
  SG::AuxElement::Accessor<float> myy("m_yy");
  SG::AuxElement::Accessor<float> met_sumet("sumet_TST");
  SG::AuxElement::Accessor<float> met_phi("phi_TST");
  SG::AuxElement::Accessor<float> met_TST("met_TST");
  SG::AuxElement::Accessor<float> phi_allparticles_met("DeltaPhi_allPart_MET");
  SG::AuxElement::Accessor<float> photons_ePS("E0_raw");
  SG::AuxElement::Accessor<float> photons_e1("E1_raw");
  SG::AuxElement::Accessor<float> photons_e2("E2_raw");
  SG::AuxElement::Accessor<float> photons_e3("E3_raw");
  SG::AuxElement::Accessor<int> photons_converted("conversionType");
  SG::AuxElement::Accessor<float> photons_convertedRadius("conversionRadius");
  SG::AuxElement::Accessor<float> hardestVertex("hardestVertexZ");
  SG::AuxElement::Accessor<float> selectedVertex("selectedVertexZ");
  SG::AuxElement::Accessor<char> passHiggsSelection("isPassed");
  SG::AuxElement::Accessor<char> passQualityCuts("isPassedBasic");
  SG::AuxElement::Accessor<char> passDoubleHiggsSelection("isPassedLowHighMyy");

  SG::AuxElement::Accessor<float> ptvarCone20("ptvarcone20");
  SG::AuxElement::Accessor<float> topoetCone20("topoetcone20");
  SG::AuxElement::Accessor<float> topoetCone40("topoetcone40");
  SG::AuxElement::Accessor<float> ptCone20("ptcone20");
  SG::AuxElement::Accessor<float> ptCone40("ptcone40");
  SG::AuxElement::Accessor<char> isisoFixedCutTight("isisoFixedCutTight");
  SG::AuxElement::Accessor<char> isisoFixedCutLoose("isisoFixedCutLoose");
  SG::AuxElement::Accessor<char> isisoFixedCutTightCaloOnly("isisoFixedCutTightCaloOnly");
  SG::AuxElement::Accessor<char> isisoFixedCutLooseCaloOnly("isisoFixedCutLooseCaloOnly");
  SG::AuxElement::Accessor<char> isTight("isTight");
  SG::AuxElement::Accessor<char> passIPcut("passIPCut");

  SG::AuxElement::Accessor<float> Jvt("Jvt");



  // photons


  //if (photons.size() < 2) return EL::StatusCode::SUCCESS;
  //TLorentzVector h = photons[0]->p4() + photons[1]->p4();
  //histoStore()->fillTH1F("m_yy", h.M()/HG::GeV);


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
  myy_var = myy(*HGameventInfo);
  metSig_var = met_TST(*HGameventInfo)/met_sumet(*HGameventInfo);
  passVertex_var = fabs(hardestVertex(*HGameventInfo) - selectedVertex(*HGameventInfo)) < 0.3 ? 1 : 0;
  passHiggsSelection_var = passHiggsSelection(*HGameventInfo) == 1 ?  1 : 0;
  passQualityCuts_var = passQualityCuts(*HGameventInfo) == 1 ? 1 : 0 ;


  std::string inputfileName = wk()->inputFileName();
  std::string cutFlowName ;

  cutFlowName = "CutFlow_" + inputfileName;
  cutFlowName.replace(cutFlowName.find(".MxAOD") , -1, "") ;
  cutFlowName.append("_weighted");

  if( m_eventCounter == 1 ) m_histCutFlow[inputfileName] = (TH1F*) HG::getHistogramFromFile(cutFlowName,inputfileName);

  xAOD::PhotonContainer photons = photonHandler()->getCorrectedContainer() ;
  xAOD::ElectronContainer electrons = electronHandler()->getCorrectedContainer() ;
  xAOD::MuonContainer muons = muonHandler()->getCorrectedContainer() ;
  xAOD::JetContainer jets = jetHandler()->getCorrectedContainer() ;

  overlapHandler()->removeOverlap(photons,jets,electrons, muons);


  TLorentzVector photon_lead;
  TLorentzVector photon_sublead;
  TLorentzVector jet_lead;
  TLorentzVector jet_sublead;

  if( photons.size() > 0 ) photon_lead = (photons)[0]->p4();
  if( photons.size() > 1 ) photon_sublead = (photons)[1]->p4();


  if( jets.size() > 0 ) jet_lead = jets[0]->p4();
  if( jets.size() > 1 ) jet_sublead = jets[1]->p4();

  TLorentzVector allparticles;

  for(xAOD::Jet* jet : jets){
    TLorentzVector test=jet->p4();
    allparticles = allparticles+test;
  }
  for(xAOD::Photon* photon: photons){
    TLorentzVector test = photon->p4();
    allparticles = allparticles +test;
  }
  for(xAOD::Electron* electron: electrons){
    TLorentzVector test = electron->p4();
    allparticles = allparticles +test;
  }
  for(xAOD::Muon* muon: muons){
    TLorentzVector test = muon->p4();
    allparticles = allparticles +test;
  }

  // switch( photons.size() ) {


  // case 0 :

  //   break;
  // case 1 :


  //   phi_y1_met_var = (photon_lead).DeltaPhi(etmissVector);
  //   if( jets.size() > 0 )phi_y1_j_var = (photon_lead).DeltaPhi(jet_lead);

  //   break;
  // case 2 : 

  //   phi_y2_met_var = (photon_sublead).DeltaPhi(etmissVector);
  //   phi_yy_met_var = (photon_lead + photon_sublead).DeltaPhi(etmissVector);
  //   phi_y1_y2_var = (photon_lead).DeltaPhi(photon_sublead);


  //   if( jets.size() > 0 ){
  //     phi_yyj_met_var = (photon_lead + photon_sublead + jet_lead).DeltaPhi(etmissVector);
  //     phi_y2_j_var = (photon_sublead).DeltaPhi(jet_lead);
  //   }
  //   if( jets.size() > 1 ){
  //     phi_yyjj_met_var = (photon_lead + photon_sublead + jet_lead + jet_sublead).DeltaPhi(etmissVector);
  //     phi_yy_jj_var = (photon_lead + photon_sublead).DeltaPhi(jet_lead + jet_sublead);
  //   }



  // default : 

  //   phi_y2_met_var = (photon_sublead).DeltaPhi(etmissVector);
  //   phi_yy_met_var = (photon_lead + photon_sublead).DeltaPhi(etmissVector);
  //   phi_y1_y2_var = (photon_lead).DeltaPhi(photon_sublead);


  //   if( jets.size() > 0 ){
  //     phi_yyj_met_var = (photon_lead + photon_sublead + jet_lead).DeltaPhi(etmissVector);
  //     phi_y2_j_var = (photon_sublead).DeltaPhi(jet_lead);
  //   }
  //   if( jets.size() > 1 ){
  //     phi_yyjj_met_var = (photon_lead + photon_sublead + jet_lead + jet_sublead).DeltaPhi(etmissVector);
  //     phi_yy_jj_var = (photon_lead + photon_sublead).DeltaPhi(jet_lead + jet_sublead);
  //   }


  //   break;

  // }



  // phi_allparticles_met_var = allparticles.DeltaPhi(etmissVector);


  nPhotons=0;


  for(size_t gn=0; gn<photons.size(); gn++) {
    //    photonAuthor[nPhotons] = photons[gn]->author();
    photonPt[nPhotons] = photons[gn]->p4().Pt();
    photonEta[nPhotons] = photons[gn]->p4().Eta();
    photonPhi[nPhotons] = photons[gn]->p4().Phi();
    photonE[nPhotons] = photons[gn]->p4().E();
    photons_Eps[nPhotons] = photons_ePS( *photons[gn] ); 
    photons_E1[nPhotons] = photons_e1( *photons[gn] ); 
    photons_E2[nPhotons] = photons_e2( *photons[gn] ); 
    photons_E3[nPhotons] = photons_e3( *photons[gn] ); 
    photons_conversionType[nPhotons] = photons_converted( *photons[gn] ); 
    photons_conversionRadius[nPhotons] = photons_convertedRadius( *photons[gn] ); 
    photons_topoCone20[nPhotons] = topoetCone20( *photons[gn] ); 
    photons_topoCone40[nPhotons] = topoetCone40( *photons[gn] ); 

    photons_Cone20[nPhotons] = ptCone20( *photons[gn] ); 
    photons_Cone40[nPhotons] = ptCone40( *photons[gn] ); 
    photons_isFixedCutTight[nPhotons] = isisoFixedCutTight( *photons[gn] ) == 1 ? 1 : 0; 
    photons_isFixedCutTightCaloOnly[nPhotons] = isisoFixedCutTightCaloOnly( *photons[gn] ) == 1 ? 1 : 0; 
    photons_isFixedCutLooseCaloOnly[nPhotons] = isisoFixedCutLooseCaloOnly( *photons[gn] ) == 1 ? 1 : 0; 
    photons_isFixedCutLoose[nPhotons] = isisoFixedCutLoose( *photons[gn] ) == 1 ? 1 : 0;
    photons_isTight[nPhotons] = isTight( *photons[gn] ) == 1 ? 1 : 0;    

    nPhotons++;
  }

  nElectrons=0;


  for(size_t gn=0; gn<electrons.size(); gn++) {
    //    electronAuthor[nElectrons] = electrons[gn]->author();
    electronPt[nElectrons] = electrons[gn]->pt();
    electronEta[nElectrons] = electrons[gn]->eta();
    electronPhi[nElectrons] = electrons[gn]->phi();
    electronE[nElectrons] = electrons[gn]->e();
    electrons_Eps[nElectrons] = electrons[gn]->caloCluster() != nullptr ? electrons[gn]->caloCluster()->energyBE(0) : 0. ; 
    electrons_E1[nElectrons] = electrons[gn]->caloCluster() != nullptr ? electrons[gn]->caloCluster()->energyBE(1) : 0. ; 
    electrons_E2[nElectrons] = electrons[gn]->caloCluster() != nullptr ? electrons[gn]->caloCluster()->energyBE(2) : 0. ; 
    electrons_E3[nElectrons] = electrons[gn]->caloCluster() != nullptr ? electrons[gn]->caloCluster()->energyBE(3) : 0. ;
    
    electrons_topoCone20[nElectrons] = topoetCone20( *electrons[gn] ); 
    electrons_ptvarCone20[nElectrons] = ptvarCone20( *electrons[gn] );
    electrons_isTight[nElectrons] = isTight( *electrons[gn] ) == 1 ? 1 : 0; 
    electrons_charge[nElectrons] = electrons[gn]->charge() ;   

    nElectrons++;
  }


  nMuons=0;


  for(size_t gn=0; gn<muons.size(); gn++) {

    muonPt[nMuons] = muons[gn]->pt();
    muonEta[nMuons] = muons[gn]->eta();
    muonPhi[nMuons] = muons[gn]->phi();
    muonE[nMuons] = muons[gn]->e();
    
    muons_topoCone20[nMuons] = topoetCone20( *muons[gn] ); 
    muons_ptvarCone20[nMuons] = ptvarCone20( *muons[gn] );
    muons_charge[nMuons] = muons[gn]->charge() ;   
    muons_passIPcut[nMuons] = passIPcut( *muons[gn] );   

    nMuons++;
  }


  nJets=0;


  for(size_t gn=0; gn<jets.size(); gn++) {
    //    jetAuthor[nJets] = jets[gn]->author();
    jetPt[nJets] = jets[gn]->pt();
    jetEta[nJets] = jets[gn]->eta();
    jetPhi[nJets] = jets[gn]->phi();
    jetJvt[nJets] = Jvt( *jets[gn] );
    jetPassSelection[nJets] = ( jetHandler()->passPtEtaCuts( jets[gn]) && jetHandler()->passJVFCut( jets[gn]) ) && jetHandler()->passJVTCut( jets[gn] ) ? 1 : 0;
    nJets++;
  }




  met = met_TST(*HGameventInfo);
  sumet = met_sumet(*HGameventInfo);
  phi_met = met_phi(*HGameventInfo);

  myEvents->Fill();
  return EL::StatusCode::SUCCESS;
}





EL::StatusCode DMHiggsAnalysis::finalize() {


  outFile->cd();
  myEvents->Write();

  for( auto it : m_histCutFlow ){
    outFile->cd();
    std::cout << " Writing the cutflow histogram of sample  : " + TString(it.first) << std::endl;
    if( !it.second ) continue;
      //HG::fatal("Cut flow histogram is empty.");
    it.second->Write();
  }

  outFile->Close();
  



  HgammaAnalysis::finalize();

  return EL::StatusCode::SUCCESS;

}
