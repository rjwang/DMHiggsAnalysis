#ifndef DMHiggsAnalysis_DMHiggsAnalysis_H
#define DMHiggsAnalysis_DMHiggsAnalysis_H

#include "HGamAnalysisFramework/HgammaAnalysis.h"

#include "PATInterfaces/SystematicRegistry.h"

#define MAXPARTICLES 99

class DMHiggsAnalysis : public HgammaAnalysis {
    // put your configuration variables here as public variables.
    // that way they can be set directly from CINT and python.
public:
    // float cutValue;



    // variables that don't get filled at submission time should be
    // protected from being send from the submission node to the worker
    // node (done by the //!)
private:

    TFile* m_outputFile;//!
    TH1F *CutFlow_; //!
    TH1F *CutFlow_noDalitz_; //!
    TH1F *CutFlow_weighted_; //!
    TH1F *CutFlow_noDalitz_weighted_; //!
    TTree *myEvents; //!

    //  ===========================  General sample information  =========================== //
    int RunNumber;
    int EventNumber;
    int LumiBlock;
//  float mcEventWeights_var;
//  float mcID_var;
    int NPV_var;
    float mu_var;
//  int isMC_var;
    float totWeight_var;
    float lumiXsecWeight_var;
    float evtWeight_var;
    float mcWeight_var;
    float pileupWeight_var;
    float vertexWeight_var;


    float qscale, x1,x2;
    int id1,id2;

//  float xsecBrFilterEff_var;
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
    /*
      int passVertex_var;
      int passHiggsSelection_var;
      int passQualityCuts_var;
    */

    int isPassed_var;
    int isPassedJetEventClean_var;
    int isDalitz_var;
    //  ===========================  Reco objects  information  =========================== //

    //
    //Photons
    //
    int nPhotons;
    float photon_Px[MAXPARTICLES];
    float photon_Py[MAXPARTICLES];
    float photon_Pz[MAXPARTICLES];
    float photon_E[MAXPARTICLES];

    /*
      float photon_Eps[MAXPARTICLES];
      float photon_E1[MAXPARTICLES];
      float photon_E2[MAXPARTICLES];
      float photon_E3[MAXPARTICLES];
      int photon_conversionType[MAXPARTICLES];
      float photon_conversionRadius[MAXPARTICLES];
      int photon_isTight[MAXPARTICLES];
      int photons_isLoose[MAXPARTICLES];
      int photons_isLoosePrime2[MAXPARTICLES];
      int photons_isLoosePrime3[MAXPARTICLES];
      int photons_isLoosePrime4[MAXPARTICLES];
      int photons_isLoosePrime5[MAXPARTICLES];
      int photon_isIsoFixedCutTight[MAXPARTICLES];
      int photon_isIsoFixedCutTightCaloOnly[MAXPARTICLES];
      int photon_isIsoFixedCutLooseCaloOnly[MAXPARTICLES];
      int photon_isIsoFixedCutLoose[MAXPARTICLES];
      float photon_Cone20[MAXPARTICLES];
      float photon_Cone40[MAXPARTICLES];
      float photon_topoCone20[MAXPARTICLES];
      float photon_topoCone40[MAXPARTICLES];
    */
    //Electrons
    int nElectrons;
    float electron_Px[MAXPARTICLES];
    float electron_Py[MAXPARTICLES];
    float electron_Pz[MAXPARTICLES];
    float electron_E[MAXPARTICLES];
    int electron_charge[MAXPARTICLES];
    /*
      int electron_Author[MAXPARTICLES];
      float electron_Eps[MAXPARTICLES];
      float electron_E1[MAXPARTICLES];
      float electron_E2[MAXPARTICLES];
      float electron_E3[MAXPARTICLES];
      int electron_isTight[MAXPARTICLES];
      int electron_isMedium[MAXPARTICLES];
      int electron_isIsoLoose[MAXPARTICLES];
      float electron_topoCone20[MAXPARTICLES];
      float electron_ptvarCone20[MAXPARTICLES];
    */

    //Muons
    int nMuons;
    float muon_Px[MAXPARTICLES];
    float muon_Py[MAXPARTICLES];
    float muon_Pz[MAXPARTICLES];
    float muon_E[MAXPARTICLES];
    int muon_charge[MAXPARTICLES];

    /*
      float muon_Eps[MAXPARTICLES];
      float muon_E1[MAXPARTICLES];
      float muon_E2[MAXPARTICLES];
      float muon_E3[MAXPARTICLES];
      int   muon_passIPcut[MAXPARTICLES];
      float muon_topoCone20[MAXPARTICLES];
      float muon_ptvarCone20[MAXPARTICLES];
    */
    /* int muon_isTight[MAXPARTICLES]; */
    /* int muon_isMedium[MAXPARTICLES]; */
    /* int muon_isLoose[MAXPARTICLES]; */
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
    /*
      float jet_Jvt[MAXPARTICLES];
      int jet_PassSelection[MAXPARTICLES];
    */

    // Met
    float met;
    float sumet;
    float phi_met;

    float met_hv;
    float sumet_hv;
    float phi_met_hv;

    float mc_sumet;
//float metSig_var;


    float vertexZ;
    float vertexZ_hv;


    // ===========================  Truth information  =========================== //
    //gen level event
    int nmcparticles;
    float mc_px[MAXPARTICLES],mc_py[MAXPARTICLES],mc_pz[MAXPARTICLES],mc_en[MAXPARTICLES];
    int mc_id[MAXPARTICLES], mc_type[MAXPARTICLES], mc_origin[MAXPARTICLES];


    // Photons

    int ntruthPhotons;
    float photonTruthPx[MAXPARTICLES];
    float photonTruthPy[MAXPARTICLES];
    float photonTruthPz[MAXPARTICLES];
    float photonTruthE[MAXPARTICLES];
    float photonTruth_ptcone20[MAXPARTICLES];
    float photonTruth_ptcone40[MAXPARTICLES];
    float photonTruth_etcone20[MAXPARTICLES];
    float photonTruth_etcone40[MAXPARTICLES];

    int photonTruth_truthOrigin[MAXPARTICLES];
    int photonTruth_truthType[MAXPARTICLES];


    // Electrons

    int ntruthElectrons;
    float electronTruthPx[MAXPARTICLES];
    float electronTruthPy[MAXPARTICLES];
    float electronTruthPz[MAXPARTICLES];
    float electronTruthE[MAXPARTICLES];

    // Muons

    int ntruthMuons;
    float muonTruthPx[MAXPARTICLES];
    float muonTruthPy[MAXPARTICLES];
    float muonTruthPz[MAXPARTICLES];
    float muonTruthE[MAXPARTICLES];


    // Jets

    int ntruthJets;
    float jetTruthPx[MAXPARTICLES];
    float jetTruthPy[MAXPARTICLES];
    float jetTruthPz[MAXPARTICLES];
    float jetTruthE[MAXPARTICLES];

    // MissingET

    float mpxTruthInt;
    float mpyTruthInt;
    float metTruthInt;
    float sumetTruthInt;
    float mpxTruthNonInt;
    float mpyTruthNonInt;
    float metTruthNonInt;
    float sumetTruthNonInt;

    float metPhiTruthInt ;
    float metPhiTruthNonInt ;

    int ntruthDarkMatters;
    float DarkMatterTruthPx[MAXPARTICLES];
    float DarkMatterTruthPy[MAXPARTICLES];
    float DarkMatterTruthPz[MAXPARTICLES];
    float DarkMatterTruthE[MAXPARTICLES];



    //std::map<std::string,TH1F*> m_histCutFlow; //!
    /* std::map<std::string,TFile*> m_outputFiles; //! */
    /* std::map<std::string,TTree*> m_outputTTree; //! */

    //std::string currentfilename;

    /*
      CP::MuonSelectionTool              *m_muonLooseSelectionTool; //!
      CP::MuonSelectionTool              *m_muonMediumSelectionTool; //!
      CP::MuonSelectionTool              *m_muonTightSelectionTool; //!
    */

public:
    // this is a standard constructor
    DMHiggsAnalysis() { }
    DMHiggsAnalysis(const char *name);
    virtual ~DMHiggsAnalysis();

    // these are the functions inherited from HgammaAnalysis

    void declareVariables();

    bool isTruthPhoton(const xAOD::TruthParticle* truth);
    bool isDarkMatter(const xAOD::TruthParticle* truth);

    bool m_isMxAOD;

    virtual EL::StatusCode createOutput();
    virtual EL::StatusCode initialize();
    virtual EL::StatusCode execute();
    virtual EL::StatusCode finalize();


    // this is needed to distribute the algorithm to the workers
    ClassDef(DMHiggsAnalysis, 1);
};

#endif // DMHiggsAnalysis_DMHiggsAnalysis_H
