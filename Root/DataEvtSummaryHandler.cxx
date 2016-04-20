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
    t_->SetBranchAddress("RunNumber",                      &evSummary_.RunNumber);
    t_->SetBranchAddress("EventNumber",                      &evSummary_.EventNumber);
    t_->SetBranchAddress("LumiBlock",                      &evSummary_.LumiBlock);

    t_->SetBranchAddress("isMC_var",                      &evSummary_.isMC_var);
    t_->SetBranchAddress("NPV_var",                      &evSummary_.NPV_var);
    t_->SetBranchAddress("mu_var",                      &evSummary_.mu_var);
    t_->SetBranchAddress("initWeight_var",                      &evSummary_.initWeight_var);




    t_->SetBranchAddress("nPhotons",                      &evSummary_.nPhotons);
    t_->SetBranchAddress("photonPx",                   evSummary_.photonPx);
    t_->SetBranchAddress("photonPy",                   evSummary_.photonPy);
    t_->SetBranchAddress("photonPz",                   evSummary_.photonPz);
    t_->SetBranchAddress("photonE",                   evSummary_.photonE);

    t_->SetBranchAddress("nElectrons",                      &evSummary_.nElectrons);
    t_->SetBranchAddress("electronPx",                   evSummary_.electronPx);
    t_->SetBranchAddress("electronPy",                   evSummary_.electronPy);
    t_->SetBranchAddress("electronPz",                   evSummary_.electronPz);
    t_->SetBranchAddress("electronE",                   evSummary_.electronE);

    t_->SetBranchAddress("nMuons",                      &evSummary_.nMuons);
    t_->SetBranchAddress("muonPx",                   evSummary_.muonPx);
    t_->SetBranchAddress("muonPy",                   evSummary_.muonPy);
    t_->SetBranchAddress("muonPz",                   evSummary_.muonPz);
    t_->SetBranchAddress("muonE",                   evSummary_.muonE);

    t_->SetBranchAddress("nJets",                      &evSummary_.nJets);
    t_->SetBranchAddress("jetPx",                   evSummary_.jetPx);
    t_->SetBranchAddress("jetPy",                   evSummary_.jetPy);
    t_->SetBranchAddress("jetPz",                   evSummary_.jetPz);
    t_->SetBranchAddress("jetE",                   evSummary_.jetE);

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
