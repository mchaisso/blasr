#ifndef BASE_PARAMETERS_H_
#define BASE_PARAMETERS_H_

#include "CommandLineParser.h"

class BaseParameters {
 public:
  virtual void Init() {}
  virtual void RegisterCommandLineOptions(CommandLineParser &clp) {}
  virtual void MakeSane();
};


#endif
