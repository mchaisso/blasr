#ifndef COMPONENTS_SDP_PARAMETERS_H_
#define COMPONENTS_SDP_PARAMETERS_H_

#include "CommandLineParser.h"


class SDPParameters {
 public:
	int sdpIndel;
  int sdpIns, sdpDel;
	int sdpTupleSize;

  void Init() {
		sdpIndel = 5;
    sdpIns   = 5;
    sdpDel   = 10;
		sdpTupleSize = 11;
  }


  void RegisterCommandLineOptions(CommandLineParser &clp) {

    clp.RegisterIntOption("sdpTupleSize", &params.sdpTupleSize, "", CommandLineParser::PositiveInteger);
    clp.RegisterIntOption("sdpindel", &params.sdpIndel, "", CommandLineParser::Integer);
    clp.RegisterIntOption("sdpIns", &params.sdpIns, "", CommandLineParser::Integer);
    clp.RegisterIntOption("sdpDel", &params.sdpDel, "", CommandLineParser::Integer);

  }  



};


#endif
