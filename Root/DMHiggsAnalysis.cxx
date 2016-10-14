#include "DMHiggsAnalysis/DMHiggsAnalysis.h"
#include "HGamAnalysisFramework/HgammaIncludes.h"
#include "HGamAnalysisFramework/HGamVariables.h"



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

bool DMHiggsAnalysis::isTruthPhoton(const xAOD::TruthParticle* truth){

  if( !truth ) return false;
  if( fabs(truth->pdgId()) != 22 ) return false;

  return true;    
    
}


bool DMHiggsAnalysis::isDarkMatter(const xAOD::TruthParticle* truth){
    
    
  if( !truth ) return false;
  if( truth->status() != 1 ) return false;
  if( fabs(truth->pdgId()) != 1000022 ) return false;

  return true;  
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

  myEvents->Branch("mc_sumet", &mc_sumet,"mc_sumet/F");


  myEvents->Branch("vertexZ", &vertexZ,"vertexZ/F");
  myEvents->Branch("vertexZ_hv", &vertexZ_hv,"vertexZ_hv/F");

  // myEvents->Branch("metTruthInt", &metTruthInt,"metTruthInt/F");
  // myEvents->Branch("metTruthNonInt", &metTruthNonInt,"metTruthNonInt/F");
  // myEvents->Branch("metPhiTruthInt", &metPhiTruthInt,"metPhiTruthInt/F");
  // myEvents->Branch("metPhiTruthNonInt", &metPhiTruthNonInt,"metPhiTruthNonInt/F");



  myEvents->Branch("ntruthPhotons", &ntruthPhotons,"ntruthPhotons/I");
  myEvents->Branch("photonTruthPx", photonTruthPx,"photonTruthPx[ntruthPhotons]/F");
  myEvents->Branch("photonTruthPy", photonTruthPy,"photonTruthPy[ntruthPhotons]/F");
  myEvents->Branch("photonTruthPz", photonTruthPz,"photonTruthPz[ntruthPhotons]/F");
  myEvents->Branch("photonTruthE",  photonTruthE, "photonTruthE[ntruthPhotons]/F");

  myEvents->Branch("ntruthDarkMatters", &ntruthDarkMatters,"ntruthDarkMatters/I");
  myEvents->Branch("DarkMatterTruthPx", DarkMatterTruthPx,"DarkMatterTruthPx[ntruthDarkMatters]/F");
  myEvents->Branch("DarkMatterTruthPy", DarkMatterTruthPy,"DarkMatterTruthPy[ntruthDarkMatters]/F");
  myEvents->Branch("DarkMatterTruthPz", DarkMatterTruthPz,"DarkMatterTruthPz[ntruthDarkMatters]/F");
  myEvents->Branch("DarkMatterTruthE",  DarkMatterTruthE, "DarkMatterTruthE[ntruthDarkMatters]/F");


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
  std::string numberFile = inputfileName;
  //currentfilename = inputfileName;
  if( config()->getBool("IsMxAOD") )
    m_isMxAOD=true;
  else
    m_isMxAOD=false;

  std::cout << " Filename : " << inputfileName << std::endl;
  if( m_isMxAOD )
    inputfileName.replace(inputfileName.find(".MxAOD") , -1, "") ;
  else
    inputfileName=config()->getStr("NTupleName");

  TString TagName = inputfileName;
  TagName.ReplaceAll("mc15c.","");

  if( m_isMxAOD )
    inputfileName.append(".NTuple.root");
  else
    inputfileName.append(".root");

  TString InputfileName = inputfileName;

  if(TString(inputfileName).Contains("165_3jets") && m_isMxAOD )
    {
      numberFile.replace(0,numberFile.find("h013."),"");
      numberFile.replace(numberFile.find(".root"),-1,"");
      inputfileName.replace(inputfileName.find(".NTuple.root") , -1, "");
      InputfileName = inputfileName;
      InputfileName+="."+numberFile+".NTuple.root";
    }

  m_outputFile = TFile::Open(InputfileName.Data(),"RECREATE");

  //if(TagName.Contains("Pythia8_WH125") || TagName.Contains("Pythia8_ZH125") ) TagName += "_big";
  //    if(TagName.Contains("Pythia8_1Hard1BremDP40")) TagName = "Pythia81Hard1BremDP40";

  // copy cutflow hists


  if( isMC() && m_isMxAOD) {
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
  SG::AuxElement::Accessor<float> myy("m_yy");
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

  if( m_isMxAOD){
    RunNumber = runNumber( *eventInfo );
    LumiBlock = lumiBlock( *eventInfo );
    EventNumber = eventNumber( *eventInfo );
    NPV_var = NPV(*HGameventInfo);
    mu_var = mu(*HGameventInfo);
    myy_var = myy(*HGameventInfo);

  

    isPassed_var =  isPassed(*HGameventInfo) == 1 ? 1 : 0 ;
    //isPassed_var=0;
    //if(var::cutFlow()>13 && 105 <= myy_var/1e3 && myy_var/1e3 < 160) isPassed_var=1;


    isPassedJetEventClean_var = isPassedJetEventClean(*HGameventInfo) == 1 ? 1 : 0 ;
    if(!isPassed_var || !isPassedJetEventClean_var) return EL::StatusCode::SUCCESS;

  }else{

    NPV_var = eventHandler()->numberOfPrimaryVertices();
    mu_var = eventHandler()->mu();
  


    //isPassed_var=0;
    //if(var::cutFlow()>13 && 105 <= myy_var/1e3 && myy_var/1e3 < 160) isPassed_var=1;




  }


  if(isMC() && m_isMxAOD) {
    isDalitz_var = isDalitz(*HGameventInfo) == 1 ? 1 : 0 ;
    if(isDalitz_var) return EL::StatusCode::SUCCESS;
  }



  xAOD::PhotonContainer photons_H = photonHandler()->getCorrectedContainer() ;
  xAOD::ElectronContainer electrons_H = electronHandler()->getCorrectedContainer() ;
  xAOD::MuonContainer muons_H = muonHandler()->getCorrectedContainer() ;
  xAOD::JetContainer jets_H = jetHandler()->getCorrectedContainer() ;


  xAOD::PhotonContainer photons = photonHandler()->applySelection(photons_H);
  xAOD::ElectronContainer electrons = electronHandler()->applySelection(electrons_H);
  xAOD::MuonContainer muons = muonHandler()->applySelection(muons_H);
  xAOD::JetContainer jets = jetHandler()->applySelection(jets_H);

  xAOD::MissingETContainer etmiss = etmissHandler()->getCorrectedContainer(&photons,&jets_H,&electrons,&muons);

  if( !m_isMxAOD ){
    isPassed_var =  HgammaAnalysis::pass(&photons,&electrons,&muons,&jets) ? 1 : 0 ;
    isPassedJetEventClean_var = HgammaAnalysis::passJetEventCleaning() == 1 ? 1 : 0 ;

  }

  if(!isPassed_var || !isPassedJetEventClean_var) return EL::StatusCode::SUCCESS;

  //histoStore()->fillTH1F("mu", eventHandler()->mu(), weight());

  if(isMC() && m_isMxAOD) {
    HgammaAnalysis::setSelectedObjects(&photons, &electrons, &muons, &jets);
    //lumiXsecWeight_var * evtWeight_var
    totWeight_var = HgammaAnalysis::weightFinal() * eventHandler()->getVar<float>("weightJvt");
    //eventHandler()->getVar<float>("weight")  * eventHandler()->getVar<float>("weightJvt");

    lumiXsecWeight_var = HgammaAnalysis::lumiXsecWeight();
    // mcWeight_var * pileupWeight_var * vertexWeight_var * weigth from Photon SF
    //        evtWeight_var = eventHandler()->getVar<float>("weight") * eventHandler()->getVar<float>("weightJvt");
    evtWeight_var = HgammaAnalysis::weight() * eventHandler()->getVar<float>("weightJvt");

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
  overlapHandler()->removeOverlap(photons,jets,electrons,muons);


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


  met = m_isMxAOD ? met_TST(*HGameventInfo) : etmiss["TST"]->met()/HG::GeV;
  sumet = m_isMxAOD ? met_sumet(*HGameventInfo) : etmiss["TST"]->sumet()/HG::GeV;
  phi_met = m_isMxAOD ? met_phi(*HGameventInfo) : etmiss["TST"]->phi()/HG::GeV;

  std::string inputfileName = wk()->inputFileName();
  TString TagName = inputfileName;
  // diphoton vertex start to default metTST from h013
  if(!TagName.Contains("h012")) {
    met_hv = m_isMxAOD ? met_hardVertexTST(*HGameventInfo) : etmiss["hardVertexTST"]->met()/HG::GeV;
    sumet_hv = m_isMxAOD ? sumet_hardVertexTST(*HGameventInfo) : etmiss["hardVertexTST"]->sumet()/HG::GeV;
    phi_met_hv = m_isMxAOD ? phi_hardVertexTST(*HGameventInfo) : etmiss["hardVertexTST"]->phi()/HG::GeV;
  }


  //Vertex
  vertexZ     = eventHandler()->selectedVertexZ();
  vertexZ_hv  = eventHandler()->hardestVertexZ();

  //
  // Generator-level information
  //
  nmcparticles = 0;
  ntruthDarkMatters = 0 ;
  ntruthPhotons = 0 ;




  if( isMC() ) {

    xAOD::TruthParticleContainer    truthPhotons = truthHandler()->getPhotons() ;
    xAOD::TruthParticleContainer    truthElectrons = truthHandler()->getElectrons();
    xAOD::TruthParticleContainer    truthMuons = truthHandler()->getMuons();
    xAOD::JetContainer              truthJets = truthHandler()->getJets();
    xAOD::MissingETContainer        truthMET = truthHandler()->getMissingET();

    const xAOD::TruthParticleContainer*   truthParticles;

    if( !m_isMxAOD ) 
      truthParticles = truthHandler()->getTruthParticles();


    for( xAOD::TruthParticle* truthpart : truthPhotons) {

      mc_px[nmcparticles] = truthpart->p4().Px()/1000.;
      mc_py[nmcparticles] = truthpart->p4().Py()/1000.;
      mc_pz[nmcparticles] = truthpart->p4().Pz()/1000.;
      mc_en[nmcparticles] = truthpart->p4().E()/1000.;
      mc_origin[nmcparticles] = m_isMxAOD ? truthOrigin( *truthpart ) : truthpart->parent()->pdgId() ;
      mc_type[nmcparticles] = m_isMxAOD ? truthType(*truthpart) : truthpart->type();
      mc_id[nmcparticles] = 22;
      nmcparticles++;


      if( !truthpart) continue;
      if( fabs(truthType(*truthpart)) != 22 ) continue;
      if( !truthpart->parent() ) continue;
      if( truthOrigin(*(truthpart->parent())) != 25 ) continue;
  
      photonTruthPx[ntruthPhotons] = truthpart->p4().Px();
      photonTruthPy[ntruthPhotons] = truthpart->p4().Py();
      photonTruthPz[ntruthPhotons] = truthpart->p4().Pz();
      photonTruthE[ntruthPhotons] = truthpart->p4().E();
      ++ntruthPhotons;
	    
    }

    ntruthPhotons = 0 ;

    if( !m_isMxAOD ){

      for( const xAOD::TruthParticle* truthpart : *truthParticles){
 
	if( !truthpart ) continue;
	//      if( truthpart->status() != 1 ) continue;
	if( fabs(truthpart->pdgId()) != 1000022 ) continue;
	DarkMatterTruthPx[ntruthDarkMatters] = truthpart->p4().Px();
	DarkMatterTruthPy[ntruthDarkMatters] = truthpart->p4().Py();
	DarkMatterTruthPz[ntruthDarkMatters] = truthpart->p4().Pz();
	DarkMatterTruthE[ntruthDarkMatters] = truthpart->p4().E();
	++ntruthDarkMatters;

      }

    }
    mc_sumet = truthMET["NonInt"]->sumet()/1000.;
    mc_px[nmcparticles] = truthMET["NonInt"]->mpx()/1000.;
    mc_py[nmcparticles] = truthMET["NonInt"]->mpy()/1000.;
    mc_pz[nmcparticles] = 0;
    mc_en[nmcparticles] = truthMET["NonInt"]->met()/1000.;
    mc_origin[nmcparticles] = -1;
    mc_type[nmcparticles] = -1;
    mc_id[nmcparticles] = 0; //met pdgid assgined as 0
    nmcparticles++;

    //	metTruthInt =  (*truthMET["Int"]+*truthMET["IntMuons"]).met() /1000.; 
    metTruthInt =  truthMET["Int"]->met() /1000.; 
    metTruthNonInt = truthMET["NonInt"]->met()/1000.;
    //	metPhiTruthInt =  (*truthMET["Int"]+*truthMET["IntMuons"]).phi()/1000.;
    metPhiTruthInt =  truthMET["Int"]->phi()/1000.;
    metPhiTruthNonInt  = truthMET["NonInt"]->phi()/1000.;

  }

  myEvents->Fill();
  return EL::StatusCode::SUCCESS;
}



EL::StatusCode DMHiggsAnalysis::finalize()
{


  m_outputFile->cd();
  if( isMC() && m_isMxAOD ) {
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
