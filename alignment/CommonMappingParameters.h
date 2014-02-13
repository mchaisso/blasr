#ifndef COMMON_MAPPING_PARAMETERS_H_
#define COMMON_MAPPING_PARAMETERS_H_

#include "algorithms/alignment/AlignmentFormats.h"

class CommonMappingParameters {
 public:
	vector<string> readsFileNames;
	bool useTitleTable;
	string titleTableName;
	string unalignedFileName;
	bool printUnaligned;
	bool refineAlignments;
  bool  printVersion;
	string genomeFileName;
	bool refineBetweenAnchorsOnly;
	bool printDiscussion;
  bool printVerboseHelp;
	int verbosity;
	bool detailedSDPAlignment, nouseDetailedSDPAlignment;
  
  void Init() {
		titleTableName  = "";
		useTitleTable = false;
		printUnaligned = false;
		unalignedFileName = "";
		refineAlignments = true;
    printVersion = false;
		genomeFileName = "";
		refineBetweenAnchorsOnly = false;
		printDiscussion = false;
    maxReadIndex = -1;
    verbosity = 0;
    printVerboseHelp = false;
		detailedSDPAlignment = true;
		nouseDetailedSDPAlignment = false;
  }


  void MakeSane() {

		if (titleTableName != "") {
			useTitleTable = true;
		}
		if (unalignedFileName != "") {
			printUnaligned = true;
		}
		if (sdpFilterType == 0) {
			detailedSDPAlignment = true;
      nouseDetailedSDPAlignment = false;
		}
		if (detailedSDPAlignment == false) {
			sdpFilterType = 1;
		}
		if (nouseDetailedSDPAlignment == true) {
			detailedSDPAlignment = false;
		}
		if (nouseDetailedSDPAlignment == false) {
			detailedSDPAlignment = true;
		}
  }

  void RegisterCommandLineOptions(CommandLineParser &clp) {
    clp.RegisterStringOption("unaligned", &unalignedFileName, "");
    clp.RegisterFlagOption("noRefineAlign", (bool*) &refineAlign, "");
    clp.RegisterFlagOption("version", (bool*)&printVersion, "");
    clp.RegisterFlagOption("rbao", &refineBetweenAnchorsOnly, "");
    clp.RegisterFlagOption("help", &printDiscussion, "");
    clp.RegisterFlagOption("h", &printVerboseHelp, "");
    clp.RegisterFlagOption("v", (bool*) &verbosity, "");
    clp.RegisterIntOption("V", &verbosity, "Specify a level of verbosity.", CommandLineParser::NonNegativeInteger);
    clp.RegisterFlagOption("useDetailedSDP", &detailedSDPAlignment, "");
    clp.RegisterFlagOption("nouseDetailedSDP", &trashbinBool, "");
  }


};


#endif
