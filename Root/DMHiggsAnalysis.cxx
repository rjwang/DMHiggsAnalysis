#include "DMHiggsAnalysis/DMHiggsAnalysis.h"
#include "HGamAnalysisFramework/HgammaIncludes.h"

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



EL::StatusCode DMHiggsAnalysis::createOutput()
{
    // Here you setup the histograms needed for you analysis. This method
    // gets called after the Handlers are initialized, so that the systematic
    // registry is already filled.

    histoStore()->createTH1F("m_yy", 60, 110, 140);


  //Create a TTree
  TFile *outfile = wk()->getOutputFile("MxAOD");
  myEvents = new TTree("DMHiggsAnalysis","DMHiggsAnalysis");
  myEvents->SetDirectory(outfile);
  myEvents->Branch("nPhotons", &nPhotons,"nPhotons/I");
  myEvents->Branch("photonPx", photonPx,"photonPx[nPhotons]/F");
  myEvents->Branch("photonPy", photonPy,"photonPy[nPhotons]/F");
  myEvents->Branch("photonPz", photonPz,"photonPz[nPhotons]/F");
  myEvents->Branch("photonE",  photonE, "photonE[nPhotons]/F");



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
    HgammaAnalysis::execute();

    // photons
    xAOD::PhotonContainer photons = photonHandler()->getCorrectedContainer();

    nPhotons=0;
    for(size_t gn=0; gn<photons.size(); gn++) {
        photonPx[nPhotons] = photons[0]->p4().Px();
        photonPy[nPhotons] = photons[0]->p4().Py();
        photonPz[nPhotons] = photons[0]->p4().Pz();
        photonE[nPhotons] = photons[0]->p4().E();

        nPhotons++;
    }

    //if (photons.size() < 2) return EL::StatusCode::SUCCESS;
    //TLorentzVector h = photons[0]->p4() + photons[1]->p4();
    //histoStore()->fillTH1F("m_yy", h.M()/HG::GeV);




    myEvents->Fill();
    return EL::StatusCode::SUCCESS;
}
