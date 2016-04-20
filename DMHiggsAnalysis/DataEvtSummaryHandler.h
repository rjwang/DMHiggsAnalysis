#ifndef dataevtsummaryhandler_h
#define dataevtsummaryhandler_h

#if !defined(__CINT__) || defined(__MAKECINT__)

#include <string.h>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <set>
#include <cmath>

#include "Math/LorentzVector.h"
#include "TMath.h"
#include "TVector2.h"
#include "TTree.h"
#include "TLorentzVector.h"

#endif

typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > LorentzVector;
typedef std::vector<LorentzVector> LorentzVectorCollection;


#define MAXPARTICLES 99



struct DataEvtSummary_t {

  int RunNumber;
  int EventNumber;
  int LumiBlock;

  int isMC_var;
  int NPV_var;
  float mu_var;
  float initWeight_var;

  //photons
  int nPhotons;
  float photonPx[MAXPARTICLES];
  float photonPy[MAXPARTICLES];
  float photonPz[MAXPARTICLES];
  float photonE[MAXPARTICLES];

  //Electrons
  int nElectrons;
  float electronPx[MAXPARTICLES];
  float electronPy[MAXPARTICLES];
  float electronPz[MAXPARTICLES];
  float electronE[MAXPARTICLES];

  //Muons
  int nMuons;
  float muonPx[MAXPARTICLES];
  float muonPy[MAXPARTICLES];
  float muonPz[MAXPARTICLES];
  float muonE[MAXPARTICLES];

  //Jets
  int nJets;
  float jetPx[MAXPARTICLES];
  float jetPy[MAXPARTICLES];
  float jetPz[MAXPARTICLES];
  float jetE[MAXPARTICLES];


};

class DataEvtSummaryHandler {
public:
    //
    DataEvtSummaryHandler();
    virtual ~DataEvtSummaryHandler();


    //current event
    DataEvtSummary_t evSummary_;
    DataEvtSummary_t &getEvent() {
        return evSummary_;
    }

    //write mode
    //bool initTree(TTree *t);
    //void fillTree();

    //read mode
    bool attachToTree(TTree *t);
    int getEntries() { return (t_ ? t_->GetEntriesFast() : 0); }
    void getEntry(int ientry) {
    	resetStruct();
    	if(t_) t_->GetEntry(ientry);
    }

    void resetStruct();

private:
    //the tree
    TTree *t_;
};




#endif
