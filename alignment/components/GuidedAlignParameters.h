#ifndef GUIDED_ALIGN_PARAMETERS_H_
#define GUIDED_ALIGN_PARAMETERS_H_

#include "CommandLineParser.h"
#include "BaseParameters.h"


class GuidedAlignParameters  : public BaseParameters {
 public:
	bool useGuidedAlign;
  int  guidedAlignBandSize;

  void Init() {
		useGuidedAlign = true;
    guidedAlignBandSize = 16;
  }

  void RegisterCommandLineOptions(CommandLineParser &clp) {
    clp.RegisterFlagOption("guidedAlign", (bool*)&useGuidedAlign, "");
    clp.RegisterFlagOption("useGuidedAlign", (bool*)&trashbinBool, "");
    clp.RegisterFlagOption("noUseGuidedAlign", (bool*)&useGuidedAlign, "");
    clp.RegisterIntOption("guidedAlignBandSize", &guidedAlignBandSize, "", CommandLineParser::PositiveInteger);	
  }

  void MakeSane() {
  }
};



#endif
