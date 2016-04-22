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

    t_->SetBranchAddress("isMC",                      	&evSummary_.isMC);
    t_->SetBranchAddress("NPV",                      	&evSummary_.NPV);
    t_->SetBranchAddress("mu",                      	&evSummary_.mu);
    t_->SetBranchAddress("initWeight",             		&evSummary_.initWeight);
    t_->SetBranchAddress("passQualityCuts",          	&evSummary_.passQualityCuts);






    t_->SetBranchAddress("nPhotons",                      	&evSummary_.nPhotons);
    t_->SetBranchAddress("photonPx",                   		evSummary_.photonPx);
    t_->SetBranchAddress("photonPy",                   		evSummary_.photonPy);
    t_->SetBranchAddress("photonPz",                   		evSummary_.photonPz);
    t_->SetBranchAddress("photonE",                   		evSummary_.photonE);
    t_->SetBranchAddress("photons_isTight",                   	evSummary_.photons_isTight);
    t_->SetBranchAddress("photons_isFixedCutTight",   		evSummary_.photons_isFixedCutTight);




    t_->SetBranchAddress("nElectrons",                      	&evSummary_.nElectrons);
    t_->SetBranchAddress("electronPx",                   	evSummary_.electronPx);
    t_->SetBranchAddress("electronPy",                   	evSummary_.electronPy);
    t_->SetBranchAddress("electronPz",                   	evSummary_.electronPz);
    t_->SetBranchAddress("electronE",                   	evSummary_.electronE);



    t_->SetBranchAddress("nMuons",                     		&evSummary_.nMuons);
    t_->SetBranchAddress("muonPx",                   		evSummary_.muonPx);
    t_->SetBranchAddress("muonPy",                   		evSummary_.muonPy);
    t_->SetBranchAddress("muonPz",                   		evSummary_.muonPz);
    t_->SetBranchAddress("muonE",                   		evSummary_.muonE);



    t_->SetBranchAddress("nJets",                      		&evSummary_.nJets);
    t_->SetBranchAddress("jetPx",                   		evSummary_.jetPx);
    t_->SetBranchAddress("jetPy",                   		evSummary_.jetPy);
    t_->SetBranchAddress("jetPz",                   		evSummary_.jetPz);
    t_->SetBranchAddress("jetE",                   		evSummary_.jetE);




    t_->SetBranchAddress("met",                      		&evSummary_.met);
    t_->SetBranchAddress("sumet",                    		&evSummary_.sumet);
    t_->SetBranchAddress("phi_met",                  		&evSummary_.phi_met);




    return true;
}


//
void DataEvtSummaryHandler::resetStruct()
{
    evSummary_.nPhotons=0;
    evSummary_.nElectrons=0;
    evSummary_.nMuons=0;
    evSummary_.nJets=0;
}

/*
//
void DataEvtSummaryHandler::fillTree()
{
    if(t_) t_->Fill();
}
*/
