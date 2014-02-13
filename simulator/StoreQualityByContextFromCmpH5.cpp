#include "files/ReaderAgglomerate.h"
#include "SMRTSequence.h"
#include "utils/FileOfFileNames.h"
#include "simulator/ContextSet.h"
#include "simulator/OutputSampleListSet.h"
#include "datastructures/alignment/CmpFile.h"
#include "data/hdf/HDFCmpFile.h"
#include "algorithms/alignment/printers/StickAlignmentPrinter.h"

class ScoredLength {
public:
  int score, length;
  int operator<(const ScoredLength &rhs) const { 
    return score < rhs.score;
  }
  ScoredLength(int s, int l) : score(s), length(l) {}
  ScoredLength() {}
};

void PrintUsage() {
	cout << "cmpH5StoreQualityByContext - grab quality values from cmp.h5 files until minimum requirements for the number of times a context has been sampled are met." << endl;
	cout << "usage: cmpH5StoreQualityByContext aligned_reads.cmp.h5  output.qbc  [options] " << endl;
	cout << "options: " << endl
			 << " -contextlength L  The length of the context to sample" << endl << endl
			 << " -minSamples S(500)Report pass if all contexts are sampled" <<endl
			 << "                   at least S times." << endl << endl
			 << " -maxSamples S(1000)Stop sampling a context once it has reached"
			 << "                   S samples." << endl << endl
       << " -onlyMaxLength" <<endl
       << "                   Store only the length of the longest subread as part of the length model." << endl;
}

int main(int argc, char* argv[]) {
	string outFileName;
	int contextLength = 5;
	int minSamples = 500;
	int maxSamples = 1000;
	if (argc < 3) {
		PrintUsage();
		exit(1);
	}
	
	int argi = 1;
  string cmpH5FileName;
	cmpH5FileName = argv[argi++];
	outFileName   = argv[argi++];
	int minAverageQual = 0;
  bool onlyMaxLength = false;

	while (argi < argc) {
		if (strcmp(argv[argi], "-contextLength") == 0) {
			contextLength = atoi(argv[++argi]);
		}
		else if (strcmp(argv[argi], "-minSamples") == 0) {
			minSamples = atoi(argv[++argi]);
		}
		else if (strcmp(argv[argi], "-maxSamples") == 0) {
			maxSamples = atoi(argv[++argi]);
		}
    else if (strcmp(argv[argi], "-onlyMaxLength") == 0) {
      onlyMaxLength = true;
    }
		else {
			PrintUsage();
			cout << "ERROR, bad option: " << argv[argi] << endl;
			exit(1);
		}
		++argi;
	}
  map<string, ScoredLength> maxLengthMap;
  OutputSampleListSet samples(contextLength);
	SMRTSequence read;

	ofstream sampleOut;
	CrucialOpen(outFileName, sampleOut, std::ios::out|std::ios::binary);
	int fileNameIndex;
	
	int numContextsReached = 0;
	int numContexts = 1 << (contextLength*2);
	ReaderAgglomerate reader;
	samples.keyLength = contextLength;
	HDFCmpFile<CmpAlignment> cmpReader;
  cmpReader.IncludeField("QualityValue");
  cmpReader.IncludeField("DeletionQV");
  cmpReader.IncludeField("InsertionQV");
  cmpReader.IncludeField("SubstitutionQV");
  cmpReader.IncludeField("SubstitutionTag");
  cmpReader.IncludeField("DeletionTag");
  cmpReader.IncludeField("PulseIndex");
  cmpReader.IncludeField("WidthInFrames");
  cmpReader.IncludeField("PreBaseFrames");

	if (cmpReader.Initialize(cmpH5FileName, H5F_ACC_RDWR) == 0) {
		cout << "ERROR, could not open the cmp file." << endl;
		exit(1);
	}
	cout << "Reading cmp file." << endl;

	CmpFile cmpFile;

  cmpReader.ReadAlignmentDescriptions(cmpFile);
  cmpReader.ReadStructure(cmpFile);
  cout << "done reading structure."<<endl;
  int alignmentIndex;
  int nAlignments = cmpReader.alnInfoGroup.GetNAlignments();
  vector<int> alignmentToBaseMap;

	for (alignmentIndex = 0; 
       alignmentIndex < nAlignments and
         !samples.Sufficient();
       alignmentIndex++) {
    //
    // For ease of use, store the length of the alignment to make another model.
    //

    ByteAlignment alignmentArray;
    cmpReader.ReadAlignmentArray(alignmentIndex, alignmentArray);
    Alignment alignment;
    ByteAlignmentToAlignment(alignmentArray, alignment);
    string readSequence, refSequence;
    readSequence.resize(alignmentArray.size());
    refSequence.resize(alignmentArray.size());
    DNASequence readDNA, refDNA;

    ByteAlignmentToQueryString(&alignmentArray[0], alignmentArray.size(), &readSequence[0]);
    ByteAlignmentToRefString(&alignmentArray[0], alignmentArray.size(), &refSequence[0]);				
    RemoveGaps(readSequence, readSequence);
    RemoveGaps(refSequence, refSequence);

    readDNA.seq = (Nucleotide*) readSequence.c_str();
    readDNA.length = readSequence.size();
    refDNA.seq = (Nucleotide*) refSequence.c_str();
    refDNA.length = refSequence.size();
    CmpAlignment cmpAlignment;

    cmpReader.ImportReadFromCmpH5(alignmentIndex, cmpAlignment, read);

    CreateAlignmentToSequenceMap(alignmentArray, alignmentToBaseMap);

    if (read.length < contextLength) {
      continue;
    }
    int subreadLength = (cmpFile.alnInfo.alignments[alignmentIndex].GetQueryEnd() - 
                         cmpFile.alnInfo.alignments[alignmentIndex].GetQueryStart());
    if (onlyMaxLength == false) {
      samples.lengths.push_back(subreadLength);
    }
    else {
      int score = (cmpAlignment.GetNMatch() - 
                   cmpAlignment.GetNMismatch() - 
                   cmpAlignment.GetNInsertions() - 
                   cmpAlignment.GetNDeletions());
      stringstream nameStrm;
      nameStrm << cmpAlignment.GetMovieId() << "_" << cmpAlignment.GetHoleNumber();
      string nameStr = nameStrm.str();
      if (maxLengthMap.find(nameStr) == maxLengthMap.end()) {
        maxLengthMap[nameStr] = ScoredLength(score, subreadLength);
      }
    }

    int sampleEnd = alignmentArray.size() - contextLength/2;
    int a;
    for (a = contextLength/2; a < sampleEnd; a++) {

      // Make sure the context begins on a real nucleotide.
      while (a < sampleEnd and 
             ((RefChar[alignmentArray[a]] == ' '))) {a++;}

			//
			// Move ab back to an index where there are contextLength/2 non-gap characters, counted by nb
			//
      int ab; // num bases
			int ae; // alignment end
      ab = a-1;
      int nb = 0, ne = 0;
      while (true) {
        if (RefChar[alignmentArray[ab]] != ' ') {
          nb++;
        }
        if (ab == 0 or nb == contextLength/2) break;
        ab--;
      }

			//
			// Advance ae to an index where there are contextLength/2 non-gap characters, counted by ne.
			//
      ae = a + 1;
      while (ae < alignmentArray.size() and ne < contextLength/ 2) {
        if (RefChar[alignmentArray[ae]] != ' ') {
          ne++;
        }
        ae++;
      }

			//
			// Make sure there are no edge effects that prevent a context of the correct length from being assigned.
			//
      if (nb + ne + 1 != contextLength) {
        continue;
      }
      int ai;
      string context;
      for (ai = ab; ai < ae; ai++) {
        if (RefChar[alignmentArray[ai]] != ' ') {
          context.push_back(RefChar[alignmentArray[ai]]);
        }
      }
      assert(context.size() == contextLength);
      //      cout << "got context: " << context << endl;
      //
      // Now create the context.
      //
      OutputSample sample;

      //
      // This context is a deletion, create that.
      //
      sample.type = OutputSample::Deletion;

      //
      // This context is either an insertion or substitution
      //
      // Look to see if the previous aligned position was an
      // insertion, and move back as far as the insertion extends.
      int aq = a-1;
      int sampleLength;
			
      if (QueryChar[alignmentArray[a]] == ' ') {
        sample.type = OutputSample::Deletion;
				sampleLength = 0;
      }
			else if (RefChar[alignmentArray[aq]] == ' ') {
        while (aq > 0 
               and RefChar[alignmentArray[aq]] == ' ' 
               and QueryChar[alignmentArray[aq]] != ' ') {
          aq--;
        }
        sample.type = OutputSample::Insertion;
				sampleLength = a - aq;
      }
      else if (QueryChar[alignmentArray[a]] == RefChar[alignmentArray[aq]]) {
        sample.type = OutputSample::Match;
				sampleLength = 1;
      }
      else {
        sample.type = OutputSample::Substitution;
				sampleLength = 1;
      }
        

      sample.Resize(sampleLength);
			if (sampleLength > 0) {
				int seqPos = alignmentToBaseMap[aq];
				if (seqPos < read.length) {
					sample.CopyFromSeq(read, seqPos, sampleLength);
					string nucs;
					int n;
					for (n = 0; n < sample.nucleotides.size(); n++) { 
						char c = sample.nucleotides[n];
						assert(c == 'A' or c == 'T' or c == 'G' or c == 'C');
						nucs.push_back(sample.nucleotides[n]); 
					}
				}
			}
			samples.AppendOutputSample(context, sample);
    }
    read.Free();
  }

  if (onlyMaxLength) {
    map<string, ScoredLength>::iterator maxScoreIt;
    for (maxScoreIt = maxLengthMap.begin(); maxScoreIt != maxLengthMap.end(); ++maxScoreIt) {
      cout << maxScoreIt->second.length << endl;
      samples.lengths.push_back(maxScoreIt->second.length); 
    }
  }


	samples.Write(sampleOut);

	return 0;
}

	
