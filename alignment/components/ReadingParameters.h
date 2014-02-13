#ifndef COMPONENTS_READING_PARAMETERS_H_
#define COMPONENTS_READING_PARAMETERS_H_

//
// This class is used to define the parameters for how to interpret IO
// of SMRT,FASTQ, and FASTA sequences.

#include "CommandLineParser.h"

class ReadingParameters {
 public:
	bool useRegionTable;
	bool useHQRegionTable;
	bool unrollCcs;
	bool useCcs;
	bool useAllSubreadsInCcs;
	bool useCcsOnly;
	string regionTableFileName;
  int   maxReadIndex;
	int startRead;
	int stride;
	vector<string> regionTableFileNames;
	bool readSeparateRegionTable;
	float subsample;
  int readIndex;
  bool useRandomSeed;
  int  randomSeed;
	bool byAdapter;
	int readsFileIndex;
	int  minReadLength;
  int  maxReadLength;
  
  void Init() {
		useRegionTable = true;
		useHQRegionTable=true;
		unrollCcs  = false;
		useCcs     = false;
		useCcsOnly = false;
		useAllSubreadsInCcs = false;
		readSeparateRegionTable = false;
		regionTableFileName = "";
		startRead = 0;
		stride    = 1;
		subsample = 1.1;
    readIndex = -1;
    useRandomSeed = false;
    randomSeed = 0;
		byAdapter = false;
		readsFileIndex = 0;
		minReadLength = 50;
    maxReadLength = 0; // means no max read length
  }
  
  void RegisterCommandLineOptions(CommandLineParser &clp) {
    clp.RegisterStringOption("regionTable", &params.regionTableFileName, "");
    clp.RegisterFlagOption("unrollCcs", &unrollCcs, "");
    clp.RegisterFlagOption("useccs", &useCcs, "");
    clp.RegisterFlagOption("useccsdenovo", &useCcsOnly, "");
    clp.RegisterFlagOption("useccsall", &useAllSubreadsInCcs, "");
    clp.RegisterFlagOption("ignoreRegions", &useRegionTable, "");
    clp.RegisterFlagOption("ignoreHQRegions", &useHQRegionTable, "");
    clp.RegisterIntOption("start", &startRead, "", CommandLineParser::NonNegativeInteger);
    clp.RegisterIntOption("stride", &stride, "", CommandLineParser::NonNegativeInteger);
    clp.RegisterFloatOption("subsample", &subsample, "", CommandLineParser::PositiveFloat);
    clp.RegisterIntOption("readIndex", &readIndex, "", CommandLineParser::NonNegativeInteger);
    clp.RegisterIntOption("maxReadIndex", &maxReadIndex, "", CommandLineParser::NonNegativeInteger);
    clp.RegisterIntOption("randomSeed", &randomSeed, "", CommandLineParser::Integer);
    clp.RegisterFlagOption("divideByAdapter", &byAdapter, "");
    clp.RegisterIntOption("minReadLength", &minReadLength, "", CommandLineParser::NonNegativeInteger);
    clp.RegisterIntOption("maxReadLength", &maxReadLength, "", CommandLineParser::NonNegativeInteger);
  }

  void MakeSane() {
		if (regionTableFileName != "") {
			useRegionTable = true;
			readSeparateRegionTable = true;
		}

		if (subsample < 1 and stride > 1) {
			cout << "ERROR, subsample and stride must be used independently." << endl;
			exit(1);
		}

    if (randomSeed != 0) {
      useRandomSeed = true;
    }
  }

};




#endif
