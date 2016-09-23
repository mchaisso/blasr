#include "files/ReaderAgglomerate.h"
#include "SMRTSequence.h"
#include "utils/FileOfFileNames.h"
#include "simulator/ContextSet.h"
#include "simulator/OutputSampleListSet.h"
#include "algorithms/alignment/readers/sam/SAMReader.h"
#include "algorithms/alignment/AlignmentUtils.h"
#include "datastructures/alignmentset/SAMToAlignmentCandidateAdapter.h"
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
	cout << "storeQualityByContext - grab quality values from sam files until minimum requirements for the number of times a context has been sampled are met." << endl;
	cout << "usage: storeQualityByContext input.sam reference.fasta output.qbc  [options] " << endl;
	cout << "options: " << endl
			 << " -contextLength L  The length of the context to sample" << endl << endl
			 << " -minSamples S(500)Report pass if all contexts are sampled" <<endl
			 << "                   at least S times." << endl << endl
			 << " -maxSamples S(1000)Stop sampling a context once it has reached"
			 << "                   S samples." << endl << endl
       << " -onlyMaxLength" <<endl
       << "                   Store only the length of the longest subread as part of the length model." << endl
			 << " -minLengthSamples N  Continue reading alignments until N read lengths are observed, even if the model is complete." << endl;
}

void RemoveGaps(string &src, string &dest) {
	int i;
	dest = "";
	for (i = 0; i< src.size(); i++) { if (src[i] != '-') dest.push_back(src[i]);}
}
int CopyFieldValues(unsigned char* dest, string src, int strand) {
	if (src.size() == 0) {
		return 0;
	}
	else {
		memcpy(dest, src.c_str(), src.size());
		if (strand == 1) {
			std::reverse(&dest[0], &dest[src.size()]);
		}
		return src.size();
	}
}

int AssignData( unsigned char* &field, string src, int strand) {
	if (field == NULL) {
		if (src.size() > 0) {
			cout << "ERROR. Copying into unallocated field." << endl;
			assert(0);
		}
		else {
			return 0;
		}
	}
	int res = CopyFieldValues(field, src, strand);
	return res;
}

int AssignQualityData( unsigned char* &field, string src, int strand) {
	int res = AssignData(field, src, strand);

	if (res == 0) {
		field = NULL;
	}		
	else {
		//
		// Transform what was printed to a numerical value.
		//
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
	int strand = samAlignment.flag & 0x10 >> 2;
	AssignData(read.seq, samAlignment.seq, strand);
	AssignQualityData(read.qual.data, samAlignment.qual, strand);
	AssignQualityData(read.deletionQV.data, samAlignment.qd, strand);
	AssignQualityData(read.insertionQV.data, samAlignment.qi, strand);
	AssignQualityData(read.substitutionQV.data, samAlignment.qs, strand);
	AssignQualityData(read.mergeQV.data, samAlignment.qm, strand);
	AssignData((unsigned char*&) read.substitutionTag, samAlignment.ts, strand);
	AssignData((unsigned char*&) read.deletionTag, samAlignment.td, strand);

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

	refAln = "";
	queryAln = "";
	CreateAlignmentStrings(convertedAlignments[0],
												 convertedAlignments[0].qAlignedSeq, convertedAlignments[0].tAlignedSeq,
												 refAln, alignStr, queryAln);

	/*	StickPrintAlignment(convertedAlignments[0],
											convertedAlignments[0].qAlignedSeq, convertedAlignments[0].tAlignedSeq, cout);
	*/
	string refGap = "";
	//
	// qAlignedSeqPos is from the beginning of the full read (full offset), but the sam string is just the subread.
	// To account for the position, subtract where the subread starts from the full offset.
	//

	// samAlignment.xs is 1-based, but qAlignedSeqPos is 0-based. Add+1 to account for this
	int readAlignStart = convertedAlignments[0].qAlignedSeqPos - samAlignment.xs + 1;
	refGap.resize(readAlignStart, '-');
	refAln = refGap + refAln;
	string newQueryAln = samAlignment.seq.substr(0, readAlignStart) + queryAln;
	queryAln = newQueryAln;
	return refAln.size();

}

int main(int argc, char* argv[]) {
	string outFileName;
	int contextLength = 5;
	int minSamples = 500;
	int maxSamples = 1000;
	int minLengthSamples = 1000;
	if (argc < 3) {
		PrintUsage();
		exit(1);
	}
	
	int argi = 1;
  string alignmentFileName;
	string referenceName = "";
	
	alignmentFileName = argv[argi++];
	referenceName     = argv[argi++];
	outFileName       = argv[argi++];
	int minAverageQual = 0;
  bool onlyMaxLength = false;

	int SAM_FILE   = 1;
		
		
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
		else if (strcmp(argv[argi], "-minLengthSamples") == 0) {
			minLengthSamples = atoi(argv[++argi]);
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

	FASTAReader reader;
	reader.Initialize(referenceName);
	reader.storeName = true;
	reader.ReadAllSequences(reference);
	int i;
	for (i = 0; i < reference.size(); i++) {
		nameToIndex[reference[i].title] = i;
	}

	ofstream sampleOut;
	CrucialOpen(outFileName, sampleOut, std::ios::out|std::ios::binary);
	int fileNameIndex;
	
	int numContextsReached = 0;
	int numContexts = 1 << (contextLength*2);

	samples.keyLength = contextLength;
	if (minSamples > maxSamples) {
		maxSamples = minSamples;
	}
	samples.minSamples = minSamples;
	samples.maxSamples = maxSamples;
  SAMReader<SAMFullReferenceSequence, SAMFullReadGroup, SAMPosAlignment> samReader;
  AlignmentSet<SAMFullReferenceSequence, SAMFullReadGroup, SAMPosAlignment> alignmentSet;

	samReader.Initialize(alignmentFileName);
	samReader.ReadHeader(alignmentSet);

	
  int alignmentIndex;

  vector<int> alignmentToBaseMap;
	long nMatch = 0, nMisMatch = 0, nInsertion=0, nDeletion=0;
	for (alignmentIndex = 0; 
			 samples.Sufficient() == false;
       alignmentIndex++) {
    //
    // For ease of use, store the length of the alignment to make another model.
    //
		string refSeq, querySeq;
		string refAln, queryAln;
		
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
			if (queryAln[i] != '-') { p++; }
		}
				
		if (readSize == 0) {
			break;
		}
    int sampleEnd = refAln.size() - contextLength/2;
    int a;

		a = contextLength/2;
		// Make sure the context begins on a real nucleotide.
		while (a < sampleEnd and
					 refAln[a] == '-') {a++;}

    for (; a < sampleEnd; a++) {


			//
			// Move ab back to an index where there are contextLength/2 non-gap characters, counted by nb
			//
      int ab; // num bases
			int ae; // alignment end
      ab = a;
      int nb = 0, ne = 0;
      while (ab < refAln.size() ) {
        if (refAln[ab] != '-') {
          nb++;
        }
        if (ab == 0 or nb == contextLength/2) break;
        ab--;
      }

			//
			// Advance ae to an index where there are contextLength/2 non-gap characters, counted by ne.
			//
      ae = a + 1;
      while (ae < refAln.size() and ne < contextLength/ 2 + 1) {
        if (refAln[ae] != '-') {
          ne++;
        }
        ae++;
      }

			//
			// Make sure there are no edge effects that prevent a context of the correct length from being assigned.
			//
      if (nb + ne != contextLength) {
        continue;
      }
      int ai;
      string context;
      for (ai = ab; ai < ae; ai++) {
        if (refAln[ai] != '-') {
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
      int sampleLength;
			int queryBasePos= -1;

      if (queryAln[a] == '-') {

        sample.type = OutputSample::Deletion;

				//
				// Move pointers to the boundaries of the deletion
				//
				int ds=a,de=a;
				while (ds > 0 and queryAln[ds] == '-') {
					ds-=1;
				}
				while (de < queryAln.size() and queryAln[de] == '-') {
					de+=1;
				}
				//
				// The sample length is the length of the deletion.
				//
				sampleLength = de-ds;
      }
			else if (refAln[a] == '-') {
				int aIns = a;
        while (aIns < refAln.size()
               and refAln[aIns] == '-' 
               and queryAln[aIns] != '-') {
          aIns++;
        }
        sample.type  = OutputSample::Insertion;
				sampleLength = aIns - a;
				assert(sampleLength > 0);
				sample.nNuc  = sampleLength;
      }
      else if (queryAln[a] == refAln[a]) {
        sample.type  = OutputSample::Match;
				sampleLength = 1;				
      }
      else {
        sample.type  = OutputSample::Substitution;
				sampleLength = 1;
      }

			if (sample.type != OutputSample::Deletion) {
				sample.Resize(sampleLength);
				if (sampleLength > 0) {
					int seqPos = alignmentToBaseMap[a];
					assert(seqPos + sampleLength < read.length);
					sample.CopyFromSeq(read, seqPos, sampleLength);
				}
				sample.nNuc = sampleLength;
			}
			samples.AppendOutputSample(context, sample);
			if (sample.type == OutputSample::Substitution) {
				assert(sample.nNuc == 1 and sample.nucleotides[0] != '\0');
			}

			
			a+=sampleLength;

			if (sample.type == OutputSample::Deletion){
				nDeletion += sampleLength;
			}
			else if (sample.type == OutputSample::Insertion) {
				nInsertion+= sampleLength;
			}
			else if (sample.type == OutputSample::Match) {
				nMatch += sampleLength;
			}
			else if (sample.type == OutputSample::Substitution) {
				nMisMatch += sampleLength;
			}
    }

		if (onlyMaxLength) {
			if (maxLengthMap.find(read.title) == maxLengthMap.end()){
				maxLengthMap[read.title] = refSeq.size();
			}
			else {
				if (maxLengthMap[read.title] < refSeq.size()) {
					maxLengthMap[read.title] = refSeq.size();
				}
			}
		}
		else {
			samples.lengths.push_back(refSeq.size());
		}
    read.Free();		
  }

	while (alignmentIndex < minLengthSamples) {
		int readSize;
		string refSeq, querySeq;
		string refAln, queryAln;
		
		readSize = SamReadAlignment(samReader,
																reference,
																nameToIndex,
																refAln, queryAln, read);
		if (readSize == 0) {
			break;
		}

		RemoveGaps(refAln, refSeq);
		RemoveGaps(queryAln, querySeq);
		
		if (maxLengthMap.find(read.title) == maxLengthMap.end()){
				maxLengthMap[read.title] = refSeq.size();
				alignmentIndex ++;				
		}
		else {
			if (maxLengthMap[read.title] < refSeq.size()) {
				maxLengthMap[read.title] = refSeq.size();
			}
		}

    read.Free();		
		refAln = queryAln = "";
	}


	float total = nInsertion + nDeletion + nMisMatch + nMatch + 1;

	cerr << "Match     " << nMatch << endl
			 << "MisMatch  " << nMisMatch << endl
			 << "Insertion " << nInsertion << endl
			 << "Deletion  " << nDeletion << endl
			 << "Accurcy " << nMatch / total << endl;

  if (onlyMaxLength) {
    map<string, int>::iterator maxLengthIt;
    for (maxLengthIt = maxLengthMap.begin(); maxLengthIt != maxLengthMap.end(); ++maxLengthIt) {
			samples.lengths.push_back(maxLengthIt->second); 
    }
  }


	samples.Write(sampleOut);

	for (int i = 0; i < reference.size(); i++) {
		reference[i].Free();
	}
	return 0;
}

	
