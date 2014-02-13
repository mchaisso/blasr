#ifndef ALIGNMETN_FILTER_PARAMETERS_H_
#define ALIGNMETN_FILTER_PARAMETERS_H_

#include "CommandLineParser.h"
#include "BaseParameters.h"


class AlignmentFilterParameters : public BaseParamters {
 public:

	int  maxScore;
  bool placeRandomly;
	float minPctIdentity;
  float maxPctIdentity;
	VectorIndex nBest;

  void Init() {
    maxScore = -200;
    placeRandomly = false;
		minPctIdentity = 0;
    maxPctIdentity = 100.1;
		nBest = 10;
  }

  void RegisterCommandLineOptions(CommandLineParser &clp) {
    clp.RegisterIntOption("maxScore", &maxScore, "", CommandLineParser::Integer);
    clp.RegisterFlagOption("placeRandomly", &placeRandomly, "");
    clp.RegisterFloatOption("minPctIdentity", &minPctIdentity, "", CommandLineParser::NonNegativeFloat);
    clp.RegisterFloatOption("maxPctIdentity", &maxPctIdentity, "", CommandLineParser::NonNegativeFloat);
    clp.RegisterIntOption("bestn", (int*) &nBest, "", CommandLineParser::PositiveInteger);
  }

  void MakeSane() {
  }
};


#endif
