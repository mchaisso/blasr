#ifndef BLASR_PARAMETERS_H_
#define BLASR_PARAMETERS_H_
#include <vector>

#include "tuples/TupleMetrics.h"
#include "datastructures/anchoring/AnchorParameters.h"
#include "qvs/QualityValue.h"
#include "CommonMappingParameters.h"
#include "CommandLineParser.h"
#include "components/BaseParameters.h"

class BlasrParameters : public BaseParameters, CommonMappingParameters {
 public:
  //
  // Parameters for global substitution, insertion, and deletion priors.
  //
  float minFractionToBeConsideredOverlapping;
	float indelRate;
	int indel;
	int refineAlign;
	int globalChainType;
	string tupleListName;
  string anchorFileName;
  string clusterFileName;
	int useTupleList;
	int useCountTable;
	string countTableName;
	int listTupleSize;
	int pValueType;
	int sortRefinedAlignments;
  bool storeMapQV;
	bool doSensitiveSearch;
	bool emulateNucmer;
  bool extendDenovoCCSSubreads;
	TupleMetrics saTupleMetrics;
	TupleMetrics sdpTupleMetrics;
	AnchorParameters anchorParameters;
	bool storeMetrics;
	bool ignoreQualities;
	bool extendFrontAlignment;
	bool extendAlignments;
  int  maxExtendDropoff;

  //
  // These may need to go to some prefiltering step.
  //
	int  minSubreadLength;
	int  minAvgQual;

	bool overlap;
	bool advanceHalf; 
	float approximateMaxInsertionRate;
	int nCandidates;
	bool doGlobalAlignment;
	float averageMismatchScore;
	bool mapSubreadsSeparately;
	string metricsFileName;
  string lcpBoundsFileName;
  string fullMetricsFileName;
	bool printSubreadTitle;


	int  chunkSize;
	int  subreadMapType;
	int  sdpFilterType;
	int  bandSize;
  int  extendBandSize;
	bool useQVScore;
	int  scoreType;
	float sdpBypassThreshold;
  bool computeAlignProbability;
  float readAccuracyPrior;
  int   substitutionPrior;
  int   globalDeletionPrior;
  bool  outputByThread;
  int   recurseOver;
  int   branchExpand;

	void Init() {
    CommonMappingParameters::Init();
    minFractionToBeConsideredOverlapping = 0.10;
		indelRate = 0.3;
		indel    = 5;
		refineAlign = 1;
		readsFileNames;
		tupleListName;
    anchorFileName = "";
		nCandidates = 10;
		pValueType = 0;
		useTupleList = 0;
		useCountTable = 0;
		countTableName = "";
		listTupleSize = 6;
		sortRefinedAlignments = 1;
		saTupleMetrics.Initialize(listTupleSize);
		sdpTupleMetrics.Initialize(sdpTupleSize);
    storeMapQV = true;
    extendDenovoCCSSubreads = false;
		storeMetrics = false;
		ignoreQualities = false;
		extendFrontAlignment = false;
		extendAlignments = false;
    maxExtendDropoff = 10;
		minSubreadLength = 0;
		minAvgQual = 0;
		overlap = false;
		advanceHalf = false;
		approximateMaxInsertionRate = 1.30;
		doGlobalAlignment = false;
		mapSubreadsSeparately=true;
		globalChainType = 0;
		metricsFileName = "";
    fullMetricsFileName = "";
		doSensitiveSearch = false;
		emulateNucmer = false;
		printSubreadTitle = true;
		subreadMapType = 0;
		chunkSize = 10000000;
		sdpFilterType = 0;
		bandSize = 0;
    extendBandSize = 10;
		useQVScore = false;
		sdpBypassThreshold = 1000000.0;
		scoreType = 0;
    computeAlignProbability = false;    
    readAccuracyPrior = 0.85;
    substitutionPrior = 20;
    globalDeletionPrior = 13;
    outputByThread = false;
    recurseOver = 10000;
    branchExpand = 1;
	}

	MappingParameters() {
		Init();
	}
	
	void MakeSane(){ 
    CommonMappingParameters::MakeSane();

		//
		// Fix all logical incompatibilities with parameters.
		//
		
		if (nCandidates < nBest) {
      nCandidates = nBest;
		}


    if (placeRandomly and nBest == 1) {
      cout << "ERROR. When attempting to select equivalently scoring reads at random " << endl
           << "the bestn parameter must be greater than one." << endl;
      exit(1);
    }
		if (sdpFilterType > 1) {
			cout << "Warning: using new filter method for SDP alignments.  The parameter is " << endl
					 << "either 0 or 1, but " << sdpFilterType << " was specified." << endl;
			sdpFilterType = 1;
		}
    //		anchorParameters.minMatchLength = minMatchLength;
		if (countTableName != "") {
			useCountTable = true;
		}
		if (metricsFileName != "" or fullMetricsFileName != "") {
			storeMetrics = true;
		}
		if (useCcsOnly) {
			useCcs = true;
		}
		if (useAllSubreadsInCcs == true) {
			useCcs = true;
		}
    
    if (subreadMapType < 0 or subreadMapType > 1) {
      cout << "Error, subreadImplType must be 0 or 1" << endl;
      exit(1);
    }

		if (emulateNucmer) {
      SetEmulateNucmer();
		}


    else if (clippingString != "") {
      cout << "ERROR, clipping should either be soft, hard, or none." << endl;
      exit(1);
    }
  }

  void SetEmulateNucmer() {
    anchorParameters.stopMappingOnceUnique = true;
    anchorParameters.advanceExactMatches   = 30;
    anchorParameters.maxAnchorsPerPosition = 1;
    sdpBypassThreshold                     = 0.75;
    sdpTupleSize                           = 15;
    anchorParameters.minMatchLength        = 30;
    useGuidedAlign                         = true;
    refineAlignments                       = false;
  }

  void SetForSensitivity() {
    advanceExactMatches = 0;
    anchorParameters.numBranches = 1;
    anchorParameters.maxAnchorsPerPosition = 10000;
  }

  void RegisterCommandLineOptions(CommandLineParser &clp) {
    clp.RegisterFlagOption("onegap", &separateGaps, "");
    clp.RegisterStringOption("ctab", &countTableName, "" );
    clp.RegisterFlagOption("extend", &extendAlignments, "");
    clp.RegisterIntOption("branchExpand", &branchExpand, "", CommandLineParser::NonNegativeInteger);
    clp.RegisterIntOption("maxExtendDropoff", &maxExtendDropoff, "", CommandLineParser::NonNegativeInteger);
    clp.RegisterFlagOption("nucmer", &emulateNucmer, "");
    clp.RegisterStringOption("anchors",  &anchorFileName, "");
    clp.RegisterStringOption("clusters", &clusterFileName, "");
    clp.RegisterFlagOption("noStoreMapQV", &storeMapQV, "");
    clp.RegisterIntOption("subreadImplType", &subreadMapType, "", CommandLineParser::PositiveInteger);
    clp.RegisterIntOption("bandSize", &bandSize, "", CommandLineParser::PositiveInteger);	
    clp.RegisterIntOption("extendBandSize", &extendBandSize, "", CommandLineParser::PositiveInteger);	
    clp.RegisterFloatOption("indelRate", &indelRate, "", CommandLineParser::NonNegativeFloat);
    clp.RegisterFloatOption("sdpbypass", &sdpBypassThreshold, "", CommandLineParser::NonNegativeFloat);
    clp.RegisterFloatOption("minFrac", &trashbinFloat, "", CommandLineParser::NonNegativeFloat);
    clp.RegisterIntOption("pvaltype", &pValueType, "", CommandLineParser::NonNegativeInteger);
    clp.RegisterFlagOption("sortRefinedAlignments",(bool*) &sortRefinedAlignments, "");
    clp.RegisterIntOption("contextAlignLength", &anchorParameters.contextAlignLength, "", CommandLineParser::PositiveInteger);
    clp.RegisterFlagOption("skipLookupTable", &anchorParameters.useLookupTable, "");
    clp.RegisterStringOption("metrics", &metricsFileName, "");
    clp.RegisterStringOption("fullMetrics", &fullMetricsFileName, "");
    clp.RegisterFlagOption("ignoreQuality", &ignoreQualities, "");
    clp.RegisterFlagOption("noFrontAlign", &extendFrontAlignment, "");
    clp.RegisterIntOption("minSubreadLength", &minSubreadLength, "", CommandLineParser::NonNegativeInteger);
    clp.RegisterIntOption("minAvgQual", &minAvgQual, "", CommandLineParser::Integer);
    clp.RegisterFlagOption("advanceHalf", &advanceHalf, "");
    clp.RegisterIntOption("maxLCPLength", &anchorParameters.maxLCPLength, "", CommandLineParser::NonNegativeInteger);
    clp.RegisterFlagOption("extendDenovoCCSSubreads", &extendDenovoCCSSubreads, "");
    clp.RegisterFlagOption("noRefineAlignments", &refineAlignments, "");
    clp.RegisterIntOption("nCandidates", &nCandidates, "", CommandLineParser::NonNegativeInteger);
    clp.RegisterFlagOption("noSplitSubreads", &mapSubreadsSeparately, "");
    clp.RegisterIntOption("subreadMapType", &subreadMapType, "", CommandLineParser::NonNegativeInteger);
    clp.RegisterStringOption("titleTable", &titleTableName, "");
    clp.RegisterFlagOption("useSensitiveSearch", &doSensitiveSearch, "");
    clp.RegisterFlagOption("computeAlignProbability", &computeAlignProbability, "");
    clp.RegisterFlagOption("global", &doGlobalAlignment, "");
    clp.RegisterIntOption("globalChainType", &globalChainType, "", CommandLineParser::NonNegativeInteger);
    clp.RegisterFlagOption("noPrintSubreadTitle", (bool*) &printSubreadTitle, "");
    clp.RegisterIntOption("sdpFilterType", &sdpFilterType, "", CommandLineParser::NonNegativeInteger);
    clp.RegisterIntOption("scoreType", &scoreType, "", CommandLineParser::NonNegativeInteger);
    clp.RegisterFloatOption("accuracyPrior",    &readAccuracyPrior, "", CommandLineParser::NonNegativeFloat);
    clp.RegisterIntOption("substitutionPrior",  &substitutionPrior, "", CommandLineParser::NonNegativeInteger);
    clp.RegisterIntOption("deletionPrior",  &globalDeletionPrior, "", CommandLineParser::NonNegativeInteger);
    clp.RegisterIntOption("recurseOver", &recurseOver, "", CommandLineParser::NonNegativeInteger);

  }
};


#endif
