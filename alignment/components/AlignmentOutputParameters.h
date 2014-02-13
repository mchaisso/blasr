#ifndef ALIGMENT_OUTPUT_PARAMETERS_H_
#define ALIGMENT_OUTPUT_PARAMETERS_H_

#include "algorithms/alignment/printers/SAMPrinter.h"
#include "CommandLineParser.h"
#include "BaseAlignment.h"

class AlignmentOutputParameters : public BaseAlignment {
 public:
  bool printHeader;
  string clippingString;
  SAMOutput::Clipping clipping;
  bool  forPicard;
  bool  separateGaps;
	int printFormat;
  bool printSAM;
	string outFileName;

  void Init() {
    printHeader = false;
    clipping = SAMOutput::none;
    clippingString = "";
    forPicard = false;
    separateGaps = false;
		printFormat = SummaryPrint;
		outFileName = "";
    printSAM = false;
  }

  void RegisterCommandLineOptions(CommandLineParser &clp) {
    clp.RegisterFlagOption("header", (bool*)&printHeader, "For column formatted output, print a line "
                           "containing the names of each column.");
    clp.RegisterFlagOption("forPicard", &forPicard, "");
        clp.RegisterStringOption("clipping", &clippingString, "");
    clp.RegisterFlagOption("sam", &printSAM, "");
    clp.RegisterIntOption("m", &printFormat, "", CommandLineParser::NonNegativeInteger);
    clp.RegisterStringOption("out", &outFileName, "");
  }
  
  void MakeSane() {

    if (printSAM) {
      printFormat = SAM;
      forPicard = true;
    }

    //
    // Parse the clipping.
    //
    if (clippingString == "soft") {
      clipping = SAMOutput::soft;
    }
    else if (clippingString == "hard") {
      clipping = SAMOutput::hard;
    }
    else if (clippingString == "none") {
      clipping = SAMOutput::none;
    }
  }
};


#endif
