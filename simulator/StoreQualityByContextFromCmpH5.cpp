#include "files/ReaderAgglomerate.h"
#include "SMRTSequence.h"
#include "utils/FileOfFileNames.h"
#include "simulator/ContextSet.h"
#include "simulator/OutputSampleListSet.h"
#include "algorithms/alignment/readers/sam/SAMReader.h"
#include "algorithms/alignment/AlignmentUtils.h"
#include "datastructures/alignmentset/SAMToAlignmentCandidateAdapter.h"

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
	cout << "storeQualityByContext - grab quality values from sam files until minimum requirements for the number of times a context has been sampled are met." << endl;
	cout << "usage: storeQualityByContext input.sam  output.qbc  [options] " << endl;
	cout << "options: " << endl
			 << " -contextlength L  The length of the context to sample" << endl << endl
		   << " -reference  ref   The reference file. Necessary when reading from SAM" << endl << endl
			 << " -minSamples S(500)Report pass if all contexts are sampled" <<endl
			 << "                   at least S times." << endl << endl
			 << " -maxSamples S(1000)Stop sampling a context once it has reached"
			 << "                   S samples." << endl << endl
       << " -onlyMaxLength" <<endl
       << "                   Store only the length of the longest subread as part of the length model." << endl;
}

void RemoveGaps(string &src, string &dest) {
	int i;
	for (i = 0; i< src.size(); i++) { if (src[i] != '-') dest.push_back(src[i]);}
}
int CopyFieldValues(unsigned char* dest, string src) {
	if (src.size() == 0) {
		return 0;
	}
	else {
		memcpy(dest, src.c_str(), src.size());
		return src.size();
	}
}

int AssignData( unsigned char* &field, string src) {
	if (field == NULL) {
		if (src.size() > 0) {
			cout << "ERROR. Copying into unallocated field." << endl;
			assert(0);
		}
		else {
			return 0;
		}
	}
	int res = CopyFieldValues(field, src);
	return res;
}

int AssignQualityData( unsigned char* &field, string src) {
	int res = AssignData(field, src);

	if (res == 0) {
		field = NULL;
	}		
	else {
		QualityStringToStored(field, src.size());
	}
	return res;
}

int SamReadAlignment(  SAMReader<SAMFullReferenceSequence, SAMFullReadGroup, SAMPosAlignment> &samReader,
											 vector<FASTASequence> &references,
											 map<string, int> &refNameToIndex,
											 string &refAln, string &queryAln,
											 SMRTSequence &read) {

  SAMAlignment samAlignment;

	if (samReader.GetNextAlignment(samAlignment) == false) {
		return 0;
	}

	read.Allocate(samAlignment.seq.size());
	read.length = samAlignment.seq.size();
	AssignData(read.seq, samAlignment.seq);
	AssignQualityData(read.qual.data, samAlignment.qual);
	AssignQualityData(read.deletionQV.data, samAlignment.qd);
	AssignQualityData(read.insertionQV.data, samAlignment.qi);
	AssignQualityData(read.substitutionQV.data, samAlignment.qs);
	AssignQualityData(read.mergeQV.data, samAlignment.qm);
	AssignData((unsigned char*&) read.substitutionTag, samAlignment.ts);
	AssignData((unsigned char*&) read.deletionTag, samAlignment.td);

	read.CopyTitle(samAlignment.qName);
	string baseName;
	baseName = GetBaseTitle(samAlignment.qName);
	read.CopyTitle(baseName);
	//
	// To fit with other methods, this has to read a vector of alignments, even
	// though the first is the only real alignment.
	//
	vector<AlignmentCandidate<> > convertedAlignments;
	SAMAlignmentsToCandidates(samAlignment, 
														references, refNameToIndex,
														convertedAlignments);

	if (convertedAlignments.size() == 0) {
		return 0;
	}

	string alignStr;

	CreateAlignmentStrings(convertedAlignments[0],
												 convertedAlignments[0].qAlignedSeq, convertedAlignments[0].tAlignedSeq,
												 refAln, alignStr, queryAln);

	return refAln.size();

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
  string alignmentFileName;
	alignmentFileName = argv[argi++];
	outFileName   = argv[argi++];
	int minAverageQual = 0;
  bool onlyMaxLength = false;

	int SAM_FILE   = 1;
		
	string referenceName = "";
	
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
		else if (strcmp(argv[argi], "-reference") == 0) {
			referenceName = argv[++argi];
		}
		else {
			PrintUsage();
			cout << "ERROR, bad option: " << argv[argi] << endl;
			exit(1);
		}
		++argi;
	}
  map<string, int> maxLengthMap;
	map<string, int> nameToIndex;
  OutputSampleListSet samples(contextLength);
	SMRTSequence read;
	vector<FASTASequence> reference;
	if (referenceName != "") {
		FASTAReader reader;
		reader.Initialize(referenceName);
		reader.storeName = true;
		reader.ReadAllSequences(reference);
		int i;
		for (i = 0; i < reference.size(); i++) {
			nameToIndex[reference[i].title] = i;
		}
	}
	ofstream sampleOut;
	CrucialOpen(outFileName, sampleOut, std::ios::out|std::ios::binary);
	int fileNameIndex;
	
	int numContextsReached = 0;
	int numContexts = 1 << (contextLength*2);
	ReaderAgglomerate reader;
	samples.keyLength = contextLength;

  SAMReader<SAMFullReferenceSequence, SAMFullReadGroup, SAMPosAlignment> samReader;
  AlignmentSet<SAMFullReferenceSequence, SAMFullReadGroup, SAMPosAlignment> alignmentSet;

	samReader.Initialize(alignmentFileName);
	samReader.ReadHeader(alignmentSet);

	
  int alignmentIndex;

  vector<int> alignmentToBaseMap;

	string refSeq, querySeq;
	string refAln, queryAln;
	for (alignmentIndex = 0; 
			 samples.Sufficient() == false;
       alignmentIndex++) {
    //
    // For ease of use, store the length of the alignment to make another model.
    //
		int readSize = 0;		
		readSize = SamReadAlignment(samReader,
																reference,
																nameToIndex,
																refAln, queryAln, read);

		RemoveGaps(refAln, refSeq);
		RemoveGaps(queryAln, querySeq);
		alignmentToBaseMap.resize(queryAln.size());
		
		int i=0, p=0;
		for (i =0;i < queryAln.size(); i++) {
			alignmentToBaseMap[i] = p;			
			if (queryAln[i] != ' ') { p++; }
		}
				
		if (readSize == 0) {
			break;
		}
    int sampleEnd = refAln.size() - contextLength/2;
    int a;
    for (a = contextLength/2; a < sampleEnd; a++) {

      // Make sure the context begins on a real nucleotide.
      while (a < sampleEnd and
						 refAln[a] == ' ') {a++;}


			//
			// Move ab back to an index where there are contextLength/2 non-gap characters, counted by nb
			//
      int ab; // num bases
			int ae; // alignment end
      ab = a-1;
      int nb = 0, ne = 0;
      while (true) {
        if (refAln[ab] != ' ') {
          nb++;
        }
        if (ab == 0 or nb == contextLength/2) break;
        ab--;
      }

			//
			// Advance ae to an index where there are contextLength/2 non-gap characters, counted by ne.
			//
      ae = a + 1;
      while (ae < refAln.size() and ne < contextLength/ 2) {
        if (refAln[ae] != ' ') {
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
        if (refAln[ai] != ' ') {
          context.push_back(refAln[ai]);
        }
      }
      assert(context.size() == contextLength);

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
			
      if (queryAln[a] == ' ') {
        sample.type = OutputSample::Deletion;
				sampleLength = 0;
      }
			else if (refAln[aq] == ' ') {
        while (aq > 0 
               and refAln[aq] == ' ' 
               and queryAln[aq] != ' ') {
          aq--;
        }
        sample.type = OutputSample::Insertion;
				sampleLength = a - aq;
      }
      else if (queryAln[a] == refAln[aq]) {
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
				}
			}
			samples.AppendOutputSample(context, sample);

			
    }

		if (onlyMaxLength) {
			if (maxLengthMap.find(read.title) == maxLengthMap.end()){
				maxLengthMap[read.title] = refSeq.size();
			}
		}
		else {
			samples.lengths.push_back(refSeq.size());
		}
    read.Free();		
  }

  if (onlyMaxLength) {
    map<string, int>::iterator maxScoreIt;
    for (maxScoreIt = maxLengthMap.begin(); maxScoreIt != maxLengthMap.end(); ++maxScoreIt) {
      samples.lengths.push_back(maxScoreIt->second); 
    }
  }


	samples.Write(sampleOut);

	return 0;
}

	
