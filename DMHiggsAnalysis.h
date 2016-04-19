#ifndef DMHiggsAnalysis_DMHiggsAnalysis_H
#define DMHiggsAnalysis_DMHiggsAnalysis_H

#include "HGamAnalysisFramework/HgammaAnalysis.h"


#define MAXPARTICLES 99

class DMHiggsAnalysis : public HgammaAnalysis
{
  // put your configuration variables here as public variables.
  // that way they can be set directly from CINT and python.
public:
  // float cutValue;



  // variables that don't get filled at submission time should be
  // protected from being send from the submission node to the worker
  // node (done by the //!)
private:

  TFile* outFile;//!
  TTree *myEvents; //!
  // TH1 *myHist; //!

  float mcID_var;
  int NPV_var;
  float mu_var;
  int isMC_var;
  float initWeight_var;
  //  float XsecLumiEffKWeight_var = XsecLumiEffKWeight(*eventInfo);
  //  float TotalWeight_var =TotalWeight(;
  float myy_var;
  float phi_yy_met_var;
  float phi_yyj_met_var;
  float phi_jj_met_var;
  float phi_yyjj_met_var;
  float phi_yy_jj_var;
  float phi_y1_j_var;
  float phi_y2_j_var; 
  float phi_y1_met_var;
  float phi_y2_met_var;
  float phi_y1_y2_var;
  float phi_j1_j2_var;
  float phi_allparticles_met_var;
  float metSig_var;
  int passVertex_var;
  int passHiggsSelection_var;
  int passQualityCuts_var;
  int passDoubleHiggsSelection_var;

  //photons
  int nPhotons;
  float photonPt[MAXPARTICLES];
  float photonEta[MAXPARTICLES];
  float photonPhi[MAXPARTICLES];
  float photonE[MAXPARTICLES];

  float photons_Eps[MAXPARTICLES];
  float photons_E1[MAXPARTICLES];
  float photons_E2[MAXPARTICLES];
  float photons_E3[MAXPARTICLES];


  int photons_conversionType[MAXPARTICLES];
  int photons_isTight[MAXPARTICLES];
  int photons_isLoose[MAXPARTICLES];
  int photons_isLoosePrime[MAXPARTICLES];

  int photons_isFixedCutTight[MAXPARTICLES];
  int photons_isFixedCutTightCaloOnly[MAXPARTICLES];
  int photons_isFixedCutLoose[MAXPARTICLES];
  int photons_isCone20[MAXPARTICLES];
  int photons_isCone40[MAXPARTICLES];
  int photons_isTopocone20[MAXPARTICLES];
  int photons_isTopocone40[MAXPARTICLES];

public:
  // this is a standard constructor
  DMHiggsAnalysis() { }
  DMHiggsAnalysis(const char *name);
  virtual ~DMHiggsAnalysis();

  // these are the functions inherited from HgammaAnalysis

  void declareVariables();
  void clearVectors();

  virtual EL::StatusCode createOutput();
  virtual EL::StatusCode initialize();
  virtual EL::StatusCode execute();
  virtual EL::StatusCode finalize();


  // this is needed to distribute the algorithm to the workers
  ClassDef(DMHiggsAnalysis, 1);
};

#endif // DMHiggsAnalysis_DMHiggsAnalysis_H