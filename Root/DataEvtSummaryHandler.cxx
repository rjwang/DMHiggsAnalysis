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

//
bool DataEvtSummaryHandler::attachToTree(TTree *t)
{
    if(t==0) return false;
    t_ = t;


    //event info

    t_->SetBranchAddress("nPhotons",                      &evSummary_.nPhotons);
    t_->SetBranchAddress("photonPx",                   evSummary_.photonPx);
    t_->SetBranchAddress("photonPy",                   evSummary_.photonPy);
    t_->SetBranchAddress("photonPz",                   evSummary_.photonPz);
    t_->SetBranchAddress("photonE",                   evSummary_.photonE);

    return true;
}


//
void DataEvtSummaryHandler::resetStruct()
{
    evSummary_.nPhotons=0;

}

//
void DataEvtSummaryHandler::fillTree()
{
    if(t_) t_->Fill();
}

