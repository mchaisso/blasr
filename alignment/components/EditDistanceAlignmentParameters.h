#ifndef EDIT_DISTANCE_ALIGNMENT_PARAMETERS_H_
#define EDIT_DISTANCE_ALIGNMENT_PARAMETERS_H_

#include "CommandLineParser.h"
#include "BaseParameters.h"


class EditDistanceAlignmentParameters  : public BaseParameters {
 public:
  int insertion;
  int deletion;
  int mismatch;
	int match;
  void Init() {

  }

  void RegisterCommandLineOptions(CommandLineParser &clp) {
    clp.RegisterIntOption("match", &match, "", CommandLineParser::Integer);
    clp.RegisterIntOption("mismatch", &mismatch, "", CommandLineParser::Integer);
    clp.RegisterIntOption("indel", &indel, "", CommandLineParser::Integer);
    clp.RegisterIntOption("insertion", &insertion, "", CommandLineParser::Integer);
    clp.RegisterIntOption("deletion", &deletion, "", CommandLineParser::Integer);
  }

  void MakeSane() {

  }
};
  




};


#endif
