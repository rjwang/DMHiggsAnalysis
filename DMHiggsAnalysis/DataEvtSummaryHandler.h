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

  int NPV;
  float mu;
  float totWeight;
  float lumiXsecWeight;
  float evtWeight;
  float pileupWeight;


//  int passVertex;

  //photons
  int nPhotons;
  float photon_Px[MAXPARTICLES];
  float photon_Py[MAXPARTICLES];
  float photon_Pz[MAXPARTICLES];
  float photon_E[MAXPARTICLES];
/*
  int photon_isTight[MAXPARTICLES];
  int photon_isIsoFixedCutTight[MAXPARTICLES];
  int photon_isIsoFixedCutLoose[MAXPARTICLES];
  int photon_isIsoFixedCutTightCaloOnly[MAXPARTICLES];
  int photon_isIsoFixedCutLooseCaloOnly[MAXPARTICLES];
*/

  //Electrons
  int nElectrons;
  float electron_Px[MAXPARTICLES];
  float electron_Py[MAXPARTICLES];
  float electron_Pz[MAXPARTICLES];
  float electron_E[MAXPARTICLES];
  int electron_charge[MAXPARTICLES];
/*
  int electron_isTight[MAXPARTICLES];
  int electron_isMedium[MAXPARTICLES];
  int electron_isIsoLoose[MAXPARTICLES];
*/

  //Muons
  int nMuons;
  float muon_Px[MAXPARTICLES];
  float muon_Py[MAXPARTICLES];
  float muon_Pz[MAXPARTICLES];
  float muon_E[MAXPARTICLES];
  int muon_charge[MAXPARTICLES];
/*
  int muon_isIsoGradientLoose[MAXPARTICLES];
  int muon_isIsoGradient[MAXPARTICLES];
  int muon_isIsoLoose[MAXPARTICLES];
*/

  //Jets
  int nJets;
  float jet_Px[MAXPARTICLES];
  float jet_Py[MAXPARTICLES];
  float jet_Pz[MAXPARTICLES];
  float jet_E[MAXPARTICLES];
//  float jet_Jvt[MAXPARTICLES];


  //MET
  float met;
  float sumet;
  float phi_met;

  float met_hv;
  float sumet_hv;
  float phi_met_hv;

  float mc_sumet;

    //gen level event
    int nmcparticles;
    float mc_px[MAXPARTICLES],mc_py[MAXPARTICLES],mc_pz[MAXPARTICLES],mc_en[MAXPARTICLES];
    int mc_id[MAXPARTICLES], mc_type[MAXPARTICLES], mc_origin[MAXPARTICLES];




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
