#ifndef INDEX_PARAMETERS_H_
#define INDEX_PARAMETERS_H_

#include "CommandLineParser.h"
#include "components/BaseParameters.h"

class IndexParameters {
 public:
  int useSuffixArray;
	int useBwt;
	string suffixArrayFileName;
	string bwtFileName;
	int useSeqDB;
	string seqDBName;
	int lookupTableLength;
	bool usePrefixLookupTable;

  void Init() {
		useSuffixArray = false;
		useBwt = false;
		suffixArrayFileName= "";
		bwtFileName = "";
		useSeqDB = 0;
		lookupTableLength = 8;
		seqDBName = "";
  }

  void RegisterCommandLineOptions(CommandLineParser &clp) {
    clp.RegisterStringOption("sa", &params.suffixArrayFileName, "Suffix array of genome.");
    clp.RegisterStringOption("bwt", &params.bwtFileName, "BWT-FM index of a genome.");
    clp.RegisterStringOption("seqdb",  &params.seqDBName, "");
    clp.RegisterIntOption("saLookupTableLength", &lookupTableLength, "", CommandLineParser::PositiveInteger);
  }

  void MakeSane() {
		if (suffixArrayFileName != "") {
			useSuffixArray = true;
		}
		if (bwtFileName != "") {
			useBwt = true;
		}
		if (useBwt and useSuffixArray) {
			cout << "ERROR, sa and bwt must be used independently." << endl;
			exit(1);
		}
  }
};

#endif
