#ifndef MATCH_ANCHOR_PARAMETERS_H_
#define MATCH_ANCHOR_PARAMETERS_H_

#include "CommandLineParser.h"
#include "BaseParameters.h"

class MatchAnchorParamters  : public BaseParamters {
 public:

	int minMatchLength, maxMatchLength;
	int maxExpand, minExpand;
  bool stopMappingOnceUnique;
  int maxAnchorsPerPosition;
	int advanceExactMatches;
  void Init() {
    //		anchorParameters.minMatchLength = minMatchLength = 14;
		minMatchLength = 14;
    maxMatchLength = 0;
		maxExpand = 0;
		minExpand = 0;
    //		anchorParameters.stopMappingOnceUnique = true;
    stopMappingOnceUnique = true;
    maxAnchorsPerPosition = 1000;
    //		anchorParameters.advanceExactMatches = advanceExactMatches = 0;
		advanceExactMatches = 0;
  }

  void RegisterCommandLineOptions(CommandLineParser &clp) {
    clp.RegisterIntOption("minMatch", &minMatchLength, "", CommandLineParser::PositiveInteger);
    clp.RegisterIntOption("maxMatch", &maxMatchLength, "", CommandLineParser::PositiveInteger);
    clp.RegisterIntOption("maxExpand", &maxExpand, "", CommandLineParser::PositiveInteger);
    clp.RegisterIntOption("minExpand", &minExpand, "", CommandLineParser::NonNegativeInteger);
    clp.RegisterIntOption("stopMappingOnceUnique", (int*) &anchorParameters.stopMappingOnceUnique, "", CommandLineParser::NonNegativeInteger);
    clp.RegisterIntOption("maxAnchorsPerPosition", &anchorParameters.maxAnchorsPerPosition, "", CommandLineParser::PositiveInteger);
    clp.RegisterIntOption("advanceExactMatches", &anchorParameters.advanceExactMatches, "", CommandLineParser::NonNegativeInteger);
  }

  void MakeSane() {
    //    if (anchorParameters.maxLCPLength != 0 and anchorParameters.maxLCPLength < anchorParameters.minMatchLength) {
    if (maxLCPLength != 0 and maxLCPLength < minMatchLength) {
      cout << "ERROR: maxMatch is less than minMatch, which will result in no hits." << endl;
    }
  }
};





#endif
