#ifndef MULTI_THREADED_PARAMETERS_H_
#define MULTI_THREADED_PARAMETERS_H_


#include "BaseParameters.h"


class MultiThreadParameters : public BaseParameters { 
 public:
  int nProc;

  void Init() {
    nProc = 1;
  }

  void RegisterCommandLineOptions(CommandLineParser &clp) {
  	clp.RegisterIntOption("nproc", &params.nProc, "Number of cores to use.", CommandLineParser::PositiveInteger);
  }


};


#endif
