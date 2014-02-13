#include "prism/RNASeqParameters.h"


int main(int argc, char* argv[]) {

  CommandLineParser clp;
  prism::RNASeqParameters params;

  prism::SetupCommandLine(clp, params);

  clp.ParseCommanLine(argc, argv, params.readsFileNames);

  






}
