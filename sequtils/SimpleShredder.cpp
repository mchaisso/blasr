#include "FASTAReader.h"
#include "FASTQSequence.h"
#include "FASTASequence.h"
#include "CommandLineParser.h"
#include "algorithms/metagenomics/FindRandomSequence.h"
#include "utils.h"
#include "statistics/statutils.h"
#include <string>
#include <sstream>
using namespace std;

int main(int argc, char* argv[]) {
	string inFileName, readsFileName;
	DNALength readLength;
	float coverage = 0;
  bool noRandInit = false;
	int numReads = -1;
	CommandLineParser clp;
  int qualityValue = 20;
  bool printFastq = false;
  int stratify = 0;
  string titleType = "pacbio";
  string fastqType = "illumina"; // or "sanger"
	clp.RegisterStringOption("inFile", &inFileName, "Reference sequence", 0);
	clp.RegisterPreviousFlagsAsHidden();
	clp.RegisterIntOption("readLength", (int*) &readLength, "The length of reads to simulate.  The length is fixed.",
												CommandLineParser::PositiveInteger, "Length of every read.", 0);
	clp.RegisterFloatOption("coverage", &coverage, "Total coverage (from which the number of reads is calculated",
													CommandLineParser::PositiveFloat, 0);
  clp.RegisterFlagOption("nonRandInit", &noRandInit, "Skip initializing the random number generator with time.");
	clp.RegisterIntOption("nReads", &numReads, "Total number of reads (from which coverage is calculated)", CommandLineParser::PositiveInteger, 0);
	clp.RegisterStringOption("readsFile", &readsFileName, "Reads output file", 0);
  clp.RegisterFlagOption("fastq", &printFastq, "Fake fastq output with constant quality value (20)");
  clp.RegisterIntOption("quality", &qualityValue, "Value to use for fastq quality", CommandLineParser::PositiveInteger);
  clp.RegisterIntOption("stratify", &stratify, "Sample a read every 'stratify' bases, rather than randomly.", CommandLineParser::PositiveInteger);
  clp.RegisterStringOption("titleType", &titleType, "Set the name of the title: 'pacbio'|'illumina'");
  clp.RegisterStringOption("fastqType", &fastqType, "Set the type of fastq: 'illumina'|'sanger'");
	vector<string> leftovers;
	clp.ParseCommandLine(argc, argv, leftovers);

  if (!noRandInit) {
    InitializeRandomGeneratorWithTime();
  }

	FASTAReader inReader;
	inReader.Init(inFileName);
	vector<FASTASequence> reference;
	
	inReader.ReadAllSequences(reference);
	ofstream readsFile;
  if (readsFileName == "") {
    cout << "ERROR.  You must specify a reads file." << endl;
    exit(1);
  }
	CrucialOpen(readsFileName, readsFile, std::ios::out);
  
  ofstream sangerFastqFile;
  if (fastqType == "sanger") {
    string sangerFastqFileName = readsFileName + ".fastq";
    CrucialOpen(sangerFastqFileName, sangerFastqFile, std::ios::out);
  }

	DNALength refLength = 0;
	int i;
	for (i = 0; i < reference.size(); i++) {
		refLength += reference[i].length;
	}
	if (numReads == -1 and coverage == 0 and stratify == 0) {
		cout << "ERROR, you must specify either coverage, nReads, or stratify." << endl;
		exit(1);
	}
	else if (numReads == -1) {
		numReads = (refLength / readLength) * coverage;
	}

  if (stratify) {
    if (!readLength) {
      cout << "ERROR. If you are using stratification, a read length must be specified." << endl;
      exit(1);
    }
  }

	DNASequence sampleSeq;
	sampleSeq.length = readLength;
	int maxRetry = 10000000;
	int retryNumber = 0;
  DNALength seqIndex, seqPos;
  if (stratify) {
    seqIndex = 0;
    seqPos   = 0;
  }
  DNALength origReadLength = readLength;
	for (i = 0; stratify or i < numReads; i++) {
    if (stratify == 0) {
			int positionFound;
      positionFound = FindRandomPos(reference, seqIndex, seqPos, readLength );
			if (positionFound == false) {
				cout << "Could not simulate a position for a read of length " << readLength << endl;
				exit(1);
			}

    }
    else {
      //
      // find the next start pos, or bail if done
      //
      if (seqPos >= reference[seqIndex].length) {
        if (seqIndex == reference.size() - 1) {
          break;
        }
        else {
          seqIndex = seqIndex + 1;
          seqPos   = 0;
          continue;
        }
      }
      readLength = min(reference[seqIndex].length - seqPos, origReadLength);
    }
		sampleSeq.seq = &reference[seqIndex].seq[seqPos];
		int j;
		int gappedRead = 0;
    string title;
    stringstream titleStrm;
    if (titleType == "pacbio") {
      titleStrm << i << "|"<< reference[seqIndex].GetName() << "|" << seqPos << "|" << seqPos + readLength;
    }
    else if (titleType == "illumina") {
      titleStrm << "SE_" << i << "_0@" << seqPos << "-"<<seqPos+readLength <<"/1";
    }
    else {
      cout << "ERROR. Bad title type " << titleType << endl;
      exit(1);
    }
    title = titleStrm.str();
    sampleSeq.length = readLength;
    if (!printFastq) {
      readsFile << ">" << title << endl;
      sampleSeq.PrintSeq(readsFile);
    }
    else {
      FASTQSequence fastqSampleSeq;
      fastqSampleSeq.CopyTitle(title);
      fastqSampleSeq.seq = sampleSeq.seq;
      fastqSampleSeq.length = sampleSeq.length;
      fastqSampleSeq.qual.data = new unsigned char[sampleSeq.length];
      fill(fastqSampleSeq.qual.data, fastqSampleSeq.qual.data + sampleSeq.length, qualityValue);
      if (fastqType == "illumina") {
        fastqSampleSeq.PrintFastq(readsFile, fastqSampleSeq.length+1);
      }
      else {
        fastqSampleSeq.PrintSeq(readsFile);
        fastqSampleSeq.PrintQual(sangerFastqFile);
      }
      delete[] fastqSampleSeq.qual.data;
      delete[] fastqSampleSeq.title;
    }
    
    if (stratify) {
      seqPos += readLength;
    }

	}
	return 0;
}
