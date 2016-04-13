#include "DMHiggsAnalysis/DMHiggsAnalysis.h"
#include "HGamAnalysisFramework/RunUtils.h"

int main(int argc, char *argv[])
{
  // Set up the job for xAOD access
  xAOD::Init().ignore();

  // Create our algorithm
  DMHiggsAnalysis *alg = new DMHiggsAnalysis("DMHiggsAnalysis");

  // Use helper to start the job
  HG::runJob(alg, argc, argv);

  return 0;
}
