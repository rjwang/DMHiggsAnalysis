#ifndef DMHiggsAnalysis_DMHiggsAnalysis_H
#define DMHiggsAnalysis_DMHiggsAnalysis_H

#include "HGamAnalysisFramework/HgammaAnalysis.h"

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
  // Tree *myTree; //!
  // TH1 *myHist; //!



public:
  // this is a standard constructor
  DMHiggsAnalysis() { }
  DMHiggsAnalysis(const char *name);
  virtual ~DMHiggsAnalysis();

  // these are the functions inherited from HgammaAnalysis
  virtual EL::StatusCode createOutput();
  virtual EL::StatusCode execute();


  // this is needed to distribute the algorithm to the workers
  ClassDef(DMHiggsAnalysis, 1);
};

#endif // DMHiggsAnalysis_DMHiggsAnalysis_H
