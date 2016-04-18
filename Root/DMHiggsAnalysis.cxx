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


  myEvents->Branch("mcID",&mcID_var,"mcID/I"); 
  myEvents->Branch("NPV",&NPV_var,"/I"); 
  myEvents->Branch("mu",&mu_var,"/I"); 
  myEvents->Branch("isMC",&isMC_var,"/I"); 
  myEvents->Branch("initWeight",&initWeight_var,"/I"); 
  myEvents->Branch("myy",&myy_var,"/I"); 
  myEvents->Branch("phi_yy_met",&phi_yy_met_var,"/F"); 
  myEvents->Branch("phi_yyj_met",&phi_yyj_met_var,"/F"); 
  myEvents->Branch("phi_yyjj_met",&phi_yyjj_met_var,"/F"); 
  myEvents->Branch("phi_jj_met",&phi_jj_met_var,"/F"); 
  myEvents->Branch("phi_yy_jj",&phi_yy_jj_var,"/F"); 
  myEvents->Branch("phi_y1_j",&phi_y1_j_var,"/F"); 
  myEvents->Branch("phi_y2_j",&phi_y2_j_var,"/F"); 
  myEvents->Branch("phi_y1_met",&phi_y1_met_var,"/F"); 
  myEvents->Branch("phi_y2_met",&phi_y2_met_var,"/F");  
  myEvents->Branch("phi_y1_y2",&phi_y1_y2_var,"/F");  
  myEvents->Branch("phi_j1_j2",&phi_j1_j2_var,"/F");  
  myEvents->Branch("phi_allparticles_met",&phi_allparticles_met_var,"/F"); 
  myEvents->Branch("metSig",&metSig_var,"/F"); 
  myEvents->Branch("passVertex",&passVertex_var,"/I"); 
  myEvents->Branch("passHiggsSelection",&passHiggsSelection_var,"/I"); 
  myEvents->Branch("passQualityCuts",&passQualityCuts_var,"/I"); 
  myEvents->Branch("passDoubleHiggsSelection",&passDoubleHiggsSelection_var,"/I"); 


  myEvents->Branch("nPhotons", &nPhotons,"nPhotons/I");
  myEvents->Branch("photonPt", photonPt,"photonPt[nPhotons]/F");
  myEvents->Branch("photonEta", photonEta,"photonEta[nPhotons]/F");
  myEvents->Branch("photonPhi", photonPhi,"photonPhi[nPhotons]/F");
  myEvents->Branch("photonE",  photonE, "photonE[nPhotons]/F");
  myEvents->Branch("photons_Eps",&photons_Eps,"/F");
  myEvents->Branch("photons_E1",&photons_E1,"/F");
  myEvents->Branch("photons_E2",&photons_E2,"/F");
  myEvents->Branch("photons_E3",&photons_E3,"/F");
  myEvents->Branch("photons_conversion",&photons_conversionType,"/I");


}


void DMHiggsAnalysis::clearVectors(){


  for( int iparticle = 0 ; iparticle < MAXPARTICLES ; ++iparticle){

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
    photons_isLoose[iparticle] =  - 9999;
    photons_isLoosePrime[iparticle] =  - 9999;

    photons_isFixedCutTight[iparticle] =  - 9999;
    photons_isFixedCutTightCaloOnly[iparticle] =  - 9999;
    photons_isFixedCutLoose[iparticle] =  - 9999;
    photons_isCone20[iparticle] =  - 9999;
    photons_isCone40[iparticle] =  - 9999;
    photons_isTopocone20[iparticle] =  - 9999;
    photons_isTopocone40[iparticle] =  - 9999;

  }

}

EL::StatusCode DMHiggsAnalysis::initialize()
{


  HgammaAnalysis::initialize();
  outFile = TFile::Open("NTuple.root","RECREATE");



  myEvents = new TTree("DMHiggsAnalysis","DMHiggsAnalysis");


  wk()->addOutput(outFile);

  
  //myEvents->SetDirectory(outFile);

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


  SG::AuxElement::Accessor<int> NPV("numberOfPrimaryVertices");
  SG::AuxElement::Accessor<float> mu("mu");
  //  SG::AuxElement::Accessor<int> isMC("isMC");
  SG::AuxElement::Accessor<float> initWeight("weightInitial");
  //  SG::AuxElement::Accessor<float> XsecLumiEffKWeight("XsecLumiEffKWeight");
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
  SG::AuxElement::Accessor<float> hardestVertex("hardestVertexZ");
  SG::AuxElement::Accessor<float> selectedVertex("selectedVertexZ");
  SG::AuxElement::Accessor<char> passHiggsSelection("isPassed");
  SG::AuxElement::Accessor<char> passQualityCuts("isPassedBasic");
  SG::AuxElement::Accessor<char> passDoubleHiggsSelection("isPassedLowHighMyy");




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



  TLorentzVector etmissVector;


  mcID_var = isMC() ? eventInfo->mcChannelNumber() : -999;
  NPV_var = NPV(*HGameventInfo);
  mu_var = mu(*HGameventInfo);
  isMC_var = isMC();
  initWeight_var = initWeight(*HGameventInfo);
  //  XsecLumiEffKWeight_var = XsecLumiEffKWeight(*HGameventInfo);
  //  TotalWeight_var =TotalWeight( ;
  myy_var = myy(*HGameventInfo);
  metSig_var = met_TST(*HGameventInfo)/met_sumet(*HGameventInfo);
  passVertex_var = fabs(hardestVertex(*HGameventInfo) - selectedVertex(*HGameventInfo)) < 0.3 ? 1 : 0;
  passHiggsSelection_var = passHiggsSelection(*HGameventInfo) == 1 ?  1 : 0;
  passQualityCuts_var = passQualityCuts(*HGameventInfo) == 1 ? 1 : 0 ;


  xAOD::PhotonContainer photons = photonHandler()->getCorrectedContainer() ;
  xAOD::ElectronContainer electrons = electronHandler()->getCorrectedContainer() ;
  xAOD::MuonContainer muons = muonHandler()->getCorrectedContainer() ;
  xAOD::JetContainer jets = jetHandler()->getCorrectedContainer() ;

  xAOD::MissingETContainer etmiss = etmissHandler()->getCorrectedContainer(&photons,&jets,&electrons,&muons);


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



  switch( photons.size() ) {


  case 0 :

    break;
  case 1 :


    phi_y1_met_var = (photon_lead).DeltaPhi(etmissVector);
    if( jets.size() > 0 )phi_y1_j_var = (photon_lead).DeltaPhi(jet_lead);

    break;
  case 2 : 

    phi_y2_met_var = (photon_sublead).DeltaPhi(etmissVector);
    phi_yy_met_var = (photon_lead + photon_sublead).DeltaPhi(etmissVector);
    phi_y1_y2_var = (photon_lead).DeltaPhi(photon_sublead);


    if( jets.size() > 0 ){
      phi_yyj_met_var = (photon_lead + photon_sublead + jet_lead).DeltaPhi(etmissVector);
      phi_y2_j_var = (photon_sublead).DeltaPhi(jet_lead);
    }
    if( jets.size() > 1 ){
      phi_yyjj_met_var = (photon_lead + photon_sublead + jet_lead + jet_sublead).DeltaPhi(etmissVector);
      phi_yy_jj_var = (photon_lead + photon_sublead).DeltaPhi(jet_lead + jet_sublead);
    }



  default : 

    phi_y2_met_var = (photon_sublead).DeltaPhi(etmissVector);
    phi_yy_met_var = (photon_lead + photon_sublead).DeltaPhi(etmissVector);
    phi_y1_y2_var = (photon_lead).DeltaPhi(photon_sublead);


    if( jets.size() > 0 ){
      phi_yyj_met_var = (photon_lead + photon_sublead + jet_lead).DeltaPhi(etmissVector);
      phi_y2_j_var = (photon_sublead).DeltaPhi(jet_lead);
    }
    if( jets.size() > 1 ){
      phi_yyjj_met_var = (photon_lead + photon_sublead + jet_lead + jet_sublead).DeltaPhi(etmissVector);
      phi_yy_jj_var = (photon_lead + photon_sublead).DeltaPhi(jet_lead + jet_sublead);
    }


    break;

  }




  phi_allparticles_met_var = allparticles.DeltaPhi(etmissVector);


  nPhotons=0;
  for(size_t gn=0; gn<photons.size(); gn++) {
    photonPt[nPhotons] = photons[gn]->p4().Pt();
    photonEta[nPhotons] = photons[gn]->p4().Eta();
    photonPhi[nPhotons] = photons[gn]->p4().Phi();
    photonE[nPhotons] = photons[gn]->p4().E();
      
    photons_Eps[nPhotons] = photons_ePS( *photons[gn] ); 
    photons_E1[nPhotons] = photons_e1( *photons[gn] ); 
    photons_E2[nPhotons] = photons_e2( *photons[gn] ); 
    photons_E3[nPhotons] = photons_e3( *photons[gn] ); 
    photons_conversionType[nPhotons] = photons_converted( *photons[gn] ); 
    nPhotons++;
  }


  myEvents->Fill();
  return EL::StatusCode::SUCCESS;
}





EL::StatusCode DMHiggsAnalysis::finalize() {


  outFile->cd();
  myEvents->Write();
  outFile->Close();
  



  HgammaAnalysis::finalize();

  return EL::StatusCode::SUCCESS;

}
