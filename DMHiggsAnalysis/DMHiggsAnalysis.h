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



  // General sample information 
  int RunNumber;
  int EventNumber;
  int LumiBlock;
  float mcEventWeights_var;

  float mcID_var;
  int NPV_var;
  float mu_var;
  int isMC_var;
  float initWeight_var;
  float xsecBrFilterEff_var;
  //  float TotalWeight_var =TotalWeight(;
  float myy_var;
  /* float phi_yy_met_var; */
  /* float phi_yyj_met_var; */
  /* float phi_jj_met_var; */
  /* float phi_yyjj_met_var; */
  /* float phi_yy_jj_var; */
  /* float phi_y1_j_var; */
  /* float phi_y2_j_var;  */
  /* float phi_y1_met_var; */
  /* float phi_y2_met_var; */
  /* float phi_y1_y2_var; */
  /* float phi_j1_j2_var; */
  /* float phi_allparticles_met_var; */
  float metSig_var;
  int passVertex_var;
  int passHiggsSelection_var;
  int passQualityCuts_var;

  //photons
  int nPhotons;
  int photonAuthor[MAXPARTICLES];
  float photonPt[MAXPARTICLES];
  float photonEta[MAXPARTICLES];
  float photonPhi[MAXPARTICLES];
  float photonE[MAXPARTICLES];
  float photons_Eps[MAXPARTICLES];
  float photons_E1[MAXPARTICLES];
  float photons_E2[MAXPARTICLES];
  float photons_E3[MAXPARTICLES];
  int photons_conversionType[MAXPARTICLES];
  float photons_conversionRadius[MAXPARTICLES];
  int photons_isTight[MAXPARTICLES];
  //   int photons_isLoose[MAXPARTICLES];
  int photons_isFixedCutTight[MAXPARTICLES];
  int photons_isFixedCutTightCaloOnly[MAXPARTICLES];
  int photons_isFixedCutLooseCaloOnly[MAXPARTICLES];
  int photons_isFixedCutLoose[MAXPARTICLES];
  float photons_Cone20[MAXPARTICLES];
  float photons_Cone40[MAXPARTICLES];
  float photons_topoCone20[MAXPARTICLES];
  float photons_topoCone40[MAXPARTICLES];

  //electrons
  int nElectrons;
  int electronAuthor[MAXPARTICLES];
  float electronPt[MAXPARTICLES];
  float electronEta[MAXPARTICLES];
  float electronPhi[MAXPARTICLES];
  float electronE[MAXPARTICLES];
  float electrons_Eps[MAXPARTICLES];
  float electrons_E1[MAXPARTICLES];
  float electrons_E2[MAXPARTICLES];
  float electrons_E3[MAXPARTICLES];
  float electrons_charge[MAXPARTICLES];
  int electrons_isTight[MAXPARTICLES];
  //   int electrons_isLoose[MAXPARTICLES];

  float electrons_topoCone20[MAXPARTICLES];
  float electrons_ptvarCone20[MAXPARTICLES];


  //muons
  int nMuons;
  float muonPt[MAXPARTICLES];
  float muonEta[MAXPARTICLES];
  float muonPhi[MAXPARTICLES];
  float muonE[MAXPARTICLES];
  float muons_Eps[MAXPARTICLES];
  float muons_E1[MAXPARTICLES];
  float muons_E2[MAXPARTICLES];
  float muons_E3[MAXPARTICLES];
  float muons_charge[MAXPARTICLES];
  int muons_passIPcut[MAXPARTICLES];


  float muons_topoCone20[MAXPARTICLES];
  float muons_ptvarCone20[MAXPARTICLES];


  //jets
  int nJets;
  float jetPt[MAXPARTICLES];
  float jetEta[MAXPARTICLES];
  float jetPhi[MAXPARTICLES];
  float jetJvt[MAXPARTICLES];
  int jetPassSelection[MAXPARTICLES];


  // Met
  float met;
  float sumet;
  float phi_met;



  std::map<std::string,TH1F*> m_histCutFlow; //!


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
