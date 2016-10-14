#include "DMHiggsAnalysis/DataEvtSummaryHandler.h"

using namespace std;


//
DataEvtSummaryHandler::DataEvtSummaryHandler()
{
}

//
DataEvtSummaryHandler::~DataEvtSummaryHandler()
{
}

/*
//
bool DataEvtSummaryHandler::initTree(TTree *t)
{
    if(t==0) return false;
    t_ = t;


    t_->Branch("nPhotons", &evSummary_.nPhotons,"nPhotons/I");
    t_->Branch("photonPx", evSummary_.photonPx,"photonPx[nPhotons]/F");
    t_->Branch("photonPy", evSummary_.photonPy,"photonPy[nPhotons]/F");
    t_->Branch("photonPz", evSummary_.photonPz,"photonPz[nPhotons]/F");
    t_->Branch("photonE",  evSummary_.photonE, "photonE[nPhotons]/F");

    return true;
}
*/

//
bool DataEvtSummaryHandler::attachToTree(TTree *t)
{
    if(t==0) return false;
    t_ = t;


    //event info
    t_->SetBranchAddress("RunNumber",                      	&evSummary_.RunNumber);
    t_->SetBranchAddress("EventNumber",                      	&evSummary_.EventNumber);
    t_->SetBranchAddress("LumiBlock",                      	&evSummary_.LumiBlock);
    t_->SetBranchAddress("NPV",                      		&evSummary_.NPV);
    t_->SetBranchAddress("mu",                      		&evSummary_.mu);
    t_->SetBranchAddress("totWeight",             		&evSummary_.totWeight);
    t_->SetBranchAddress("lumiXsecWeight",             		&evSummary_.lumiXsecWeight);
    t_->SetBranchAddress("evtWeight",             		&evSummary_.evtWeight);
    t_->SetBranchAddress("pileupWeight",             		&evSummary_.pileupWeight);
//    t_->SetBranchAddress("passVertex",          		&evSummary_.passVertex);



    t_->SetBranchAddress("nPhotons",                      	&evSummary_.nPhotons);
    t_->SetBranchAddress("photon_Px",                   	evSummary_.photon_Px);
    t_->SetBranchAddress("photon_Py",                   	evSummary_.photon_Py);
    t_->SetBranchAddress("photon_Pz",                   	evSummary_.photon_Pz);
    t_->SetBranchAddress("photon_E",                   		evSummary_.photon_E);
/*
    t_->SetBranchAddress("photon_isTight",                   	evSummary_.photon_isTight);
    t_->SetBranchAddress("photon_isIsoFixedCutTight",   	evSummary_.photon_isIsoFixedCutTight);
    t_->SetBranchAddress("photon_isIsoFixedCutLoose",   	evSummary_.photon_isIsoFixedCutLoose);
    t_->SetBranchAddress("photon_isIsoFixedCutTightCaloOnly",  evSummary_.photon_isIsoFixedCutTightCaloOnly);
    t_->SetBranchAddress("photon_isIsoFixedCutLooseCaloOnly",  evSummary_.photon_isIsoFixedCutLooseCaloOnly);
*/


    t_->SetBranchAddress("nElectrons",                      	&evSummary_.nElectrons);
    t_->SetBranchAddress("electron_Px",                   	evSummary_.electron_Px);
    t_->SetBranchAddress("electron_Py",                   	evSummary_.electron_Py);
    t_->SetBranchAddress("electron_Pz",                   	evSummary_.electron_Pz);
    t_->SetBranchAddress("electron_E",                   	evSummary_.electron_E);
    t_->SetBranchAddress("electron_charge",                     evSummary_.electron_charge);
/*
    t_->SetBranchAddress("electron_isTight",                   	evSummary_.electron_isTight);
    t_->SetBranchAddress("electron_isMedium",                   evSummary_.electron_isMedium);
    t_->SetBranchAddress("electron_isIsoLoose",                 evSummary_.electron_isIsoLoose);
*/


    t_->SetBranchAddress("nMuons",                     		&evSummary_.nMuons);
    t_->SetBranchAddress("muon_Px",                   		evSummary_.muon_Px);
    t_->SetBranchAddress("muon_Py",                   		evSummary_.muon_Py);
    t_->SetBranchAddress("muon_Pz",                   		evSummary_.muon_Pz);
    t_->SetBranchAddress("muon_E",                   		evSummary_.muon_E);
    t_->SetBranchAddress("muon_charge",                         evSummary_.muon_charge);
/*
    t_->SetBranchAddress("muon_isIsoGradientLoose",             evSummary_.muon_isIsoGradientLoose);
    t_->SetBranchAddress("muon_isIsoGradient",                  evSummary_.muon_isIsoGradient);
    t_->SetBranchAddress("muon_isIsoLoose",                 	evSummary_.muon_isIsoLoose);
*/


    t_->SetBranchAddress("nJets",                      		&evSummary_.nJets);
    t_->SetBranchAddress("jet_Px",                   		evSummary_.jet_Px);
    t_->SetBranchAddress("jet_Py",                   		evSummary_.jet_Py);
    t_->SetBranchAddress("jet_Pz",                   		evSummary_.jet_Pz);
    t_->SetBranchAddress("jet_E",                   		evSummary_.jet_E);

    t_->SetBranchAddress("ntruthDarkMatters",                      		&evSummary_.ntruthDarkMatters);
    t_->SetBranchAddress("DarkMatterTruthPx",                   		evSummary_.DarkMatterTruth_Px);
    t_->SetBranchAddress("DarkMatterTruthPy",                   		evSummary_.DarkMatterTruth_Py);
    t_->SetBranchAddress("DarkMatterTruthPz",                   		evSummary_.DarkMatterTruth_Pz);
    t_->SetBranchAddress("DarkMatterTruthE",                   		evSummary_.DarkMatterTruth_E);
//    t_->SetBranchAddress("jet_Jvt",                   		evSummary_.jet_Jvt);

    t_->SetBranchAddress("vertexZ",                                 &evSummary_.vertexZ);
    t_->SetBranchAddress("vertexZ_hv",                               &evSummary_.vertexZ_hv);




    t_->SetBranchAddress("met",                      		&evSummary_.met);
    t_->SetBranchAddress("sumet",                    		&evSummary_.sumet);
    t_->SetBranchAddress("phi_met",                  		&evSummary_.phi_met);

    t_->SetBranchAddress("met_hv",                                      &evSummary_.met_hv);
    t_->SetBranchAddress("sumet_hv",                                    &evSummary_.sumet_hv);
    t_->SetBranchAddress("phi_met_hv",                                  &evSummary_.phi_met_hv);

    t_->SetBranchAddress("mc_sumet",                               &evSummary_.mc_sumet);


    //mc truth
    t_->SetBranchAddress("nmcparticles",  &evSummary_.nmcparticles);
    t_->SetBranchAddress("mc_px",         evSummary_.mc_px);
    t_->SetBranchAddress("mc_py",         evSummary_.mc_py);
    t_->SetBranchAddress("mc_pz",         evSummary_.mc_pz);
    t_->SetBranchAddress("mc_en",         evSummary_.mc_en);
    t_->SetBranchAddress("mc_id",         evSummary_.mc_id);
    t_->SetBranchAddress("mc_type",       evSummary_.mc_type);
    t_->SetBranchAddress("mc_origin",     evSummary_.mc_origin);



    return true;
}


//
void DataEvtSummaryHandler::resetStruct()
{
    evSummary_.nPhotons=0;
    evSummary_.nElectrons=0;
    evSummary_.nMuons=0;
    evSummary_.nJets=0;
    evSummary_.nmcparticles=0;
}

/*
//
void DataEvtSummaryHandler::fillTree()
{
    if(t_) t_->Fill();
}
*/
