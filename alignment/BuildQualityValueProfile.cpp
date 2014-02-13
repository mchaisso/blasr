#include "../common/data/hdf/HDFBasReader.h"
#include "../common/datastructures/reads/ReadList.h"
#include "../common/datastructures/alignment/Alignment.h"
#include "../common/qvs/QualityValueProfile.h"
#include "../common/algorithms/alignment/readers/CompareSequencesAlignmentReader.h"
#include "../common/utils/FileOfFileNames.h"
#include "../common/files/ReaderAgglomerate.h"
#include "../common/CommandLineParser.h"
#include "../common/FASTQSequence.h"
#include "../common/SMRTSequence.h"

using namespace std;
int main(int argc, char* argv[]) {

	CommandLineParser clp;
	string alignmentsFileName;
	vector<string> readsFileNames;
	string readsFileName;
	int wordSize;
	string matchProfileName, mismatchProfileName;
	clp.RegisterStringOption("alignment_file", &alignmentsFileName, "File containing alignments of reads in HDF format.");
	clp.RegisterStringOption("reads_file", &readsFileName, "Reads in .bas.h5 format, or a file containing a list of file names.");
	clp.RegisterIntOption("word_size", &wordSize, "Length of keyword", CommandLineParser::PositiveInteger);
	clp.RegisterStringOption("match_profile", &matchProfileName, "A file containing the profile of quality values of matching nucleotides.");
	clp.RegisterStringOption("mismatch_profile", &mismatchProfileName, "A file containing the profile of quality values of mismatching nucleotides.");
	clp.RegisterPreviousFlagsAsHidden();
	vector<string> opts;
	clp.ParseCommandLine(argc, argv, opts);

	
  ReadList<SMRTSequence> reads;
	if (FileOfFileNames::IsFOFN(readsFileName)) {
		FileOfFileNames::FOFNToList(readsFileName, readsFileNames);
	}
	else {
		readsFileNames.push_back(readsFileName);
	}
	//
	// Open the alignments file now to make sure it exists.  Try to open
	// as many files as possible, so open this and the output file
	// before doing anything to them. 
	//
	ifstream alignmentsIn;
	CrucialOpen(alignmentsFileName, alignmentsIn);

	ofstream matchProfileOut;
	ofstream mismatchProfileOut;
	CrucialOpen(matchProfileName, matchProfileOut, std::ios::out);
	CrucialOpen(mismatchProfileName, mismatchProfileOut, std::ios::out);

	int readsFileIndex;
	for (readsFileIndex = 0; readsFileIndex < readsFileNames.size()-2; readsFileIndex++) {
		ReaderAgglomerate reader;
		reader.Initialize(readsFileNames[readsFileIndex]);
		SMRTSequence read, *readPtr;
		cout << "reading from " << readsFileNames[readsFileIndex] << endl;
		while(reader.GetNext(read)){ 
			readPtr = new SMRTSequence;
			readPtr->Assign(read);
			reads.Add(readPtr);
		}
	}
	// make the reads list searchable.
	reads.Order();

	// Read in alignments.
	std::vector<CompSeqAlignment*> alignmentPtrs;
	CompSeqAlignment *alignmentPtr;
	while(1) {
		alignmentPtr = new CompSeqAlignment;
		if (CompareSequencesAlignmentReader<CompSeqAlignment>::Read(alignmentsIn, *alignmentPtr) == 0) {
			delete alignmentPtr;
			break;
		}
		alignmentPtrs.push_back(alignmentPtr);
	}

	//
	// Now process all of the alignments.
	//
  
  int alignmentIndex;
  VectorIndex readIndex;
	QualityValueProfile matchQualityProfile(wordSize, 256), mismatchQualityProfile(wordSize, 256);
	Nucleotide *word = new Nucleotide[wordSize+1];
	word[wordSize] = '\0';

	for (alignmentIndex = 0; alignmentIndex < alignmentPtrs.size(); alignmentIndex++) {
		if (reads.LookupReadByName(alignmentPtrs[alignmentIndex]->qName, readIndex)) {
			//
			// Found the read corresponding to this alignment, add its quality values to the profile.
			//
			int alignPos;
			CompSeqAlignment* alignment = alignmentPtrs[alignmentIndex];
			int wordPos = 0;
			int queryPos = alignment->qPos;
			for (alignPos = 0; alignPos < alignment->qString.size() - wordSize + 1; alignPos++ ){
				//
				// Update the word that is being used for the quality value profile.
				//
				if (alignment->tString[alignPos] != '-') {
					if (wordPos == wordSize) {
						int w;
						for (w = 0; w < wordSize - 1; w++) {
							word[w] = word[w+1];
						}
						word[wordPos-1] = alignment->tString[alignPos];
					}
					else {
						word[wordPos] = alignment->tString[alignPos];
						wordPos++;
					}
				}
				
				if (alignment->qString[alignPos] != '-' and alignment->qString[alignPos] == alignment->tString[alignPos]) {
					if (wordPos == wordSize) {
						matchQualityProfile.Update(word, reads.reads[readIndex]->qual[queryPos]);
					}
				}
				else {
					if (wordPos == wordSize) {
						mismatchQualityProfile.Update(word, reads.reads[readIndex]->qual[queryPos]);
					}
				}
				// 
				// Advance count in the query if 
				if (alignment->qString[alignPos] != '-') {
					queryPos++;
				}
			}
		}
	}
	delete[] word;
	// Done computing quality profiles, exit.
	matchQualityProfile.ProfileToCDF();
	mismatchQualityProfile.ProfileToCDF();
	matchQualityProfile.Print(matchProfileOut);
	mismatchQualityProfile.Print(mismatchProfileOut);

	return 0;
}
