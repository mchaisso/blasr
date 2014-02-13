#include "../common/algorithms/alignment/AffineKBandAlign.h"
#include "../common/algorithms/alignment/SWAlign.h"
#include "../common/algorithms/alignment/AlignmentUtils.h"
#include "../common/algorithms/alignment/QualityValueScoreFunction.h"
#include "../common/algorithms/alignment/DistanceMatrixScoreFunction.h"
#include "../common/algorithms/alignment/AlignmentPrinter.h"
#include "../common/datastructures/alignment/Path.h"
#include "../common/datastructures/alignment/AlignmentCandidate.h"
#include "../common/algorithms/alignment/ScoreMatrices.h"
#include "../common/utils/FileOfFileNames.h"
#include "../common/FASTQSequence.h"
#include "../common/FASTAReader.h"
#include "../common/FASTASequence.h"
#include "../common/tuples/DNATuple.h"
#include "../common/tuples/TupleMetrics.h"
#include "../common/CommandLineParser.h"
#include "../common/files/ReaderAgglomerate.h"
#include "../common/data/hdf/HDFRegionTableReader.h"

#include <string>
#include <semaphore.h>

using namespace std;


class ReadKeyword {
public:
	DNATuple tuple;
	int readIndex;
	int readPos;
	int operator<(const ReadKeyword &rhs) const {
		return tuple < rhs.tuple;
	}
	ReadKeyword & operator=(const ReadKeyword &rhs ) {
		tuple     = rhs.tuple;
		readIndex = rhs.readIndex;
		readPos   = rhs.readPos;
		return *this;
	}
	
};


typedef AlignmentCandidate<FASTQSequence, FASTQSequence> FastqAlignment;


class Data {
public:
	//	vector<ReadKeyword> *readKeywords;
	std::vector<int> *prevAlignedGenomePos;
	std::vector<int> *readOptScore;
	std::vector<FastqAlignment > *optAlignment;
	std::vector<int> *optGenomeAlignPos;
	std::vector<int> *optGenomeAlignLength;
	std::vector<ReadKeyword> *keywords;
	std::vector<FASTQSequence> *reads;
	float insRate;
	FASTQSequence *genome;
	TupleMetrics *tm;
};

	
void KeywordSeededAlignment(Data *data) {
	FASTQSequence genomeSubstring;
	DNATuple genomeTuple;
  DNALength genomePos;
	ReadKeyword genomeKeyword;
	std::vector<ReadKeyword>::iterator keyIt, upKeyIt;

	//
	// Scan the genome.
	//

	vector<int> scoreMat;
	vector<Arrow> pathMat;
	vector<Arrow> hpInsPathMat, insPathMat;
	vector<int> hpInsScoreMat, insScoreMat;
	DistanceMatrixScoreFunction<DNASequence, FASTQSequence> distanceMatrixScoreFn;
	distanceMatrixScoreFn.InitializeScoreMatrix(SMRTDistanceMatrix);
	distanceMatrixScoreFn.del = 6;
	distanceMatrixScoreFn.ins = 6;
	genomePos = 0;


	for (genomePos = 0; genomePos < data->genome->length - data->tm->tupleSize + 1; genomePos++) {
		genomeKeyword.tuple.FromStringLR(&data->genome->seq[genomePos], *data->tm);
		keyIt = lower_bound(data->keywords->begin(), data->keywords->end(), genomeKeyword);
		upKeyIt  = upper_bound(data->keywords->begin(), data->keywords->end(), genomeKeyword);
		//
		// Find all the reads and all the positions in reads that
		// have this keyword.
		for (; keyIt != upKeyIt; keyIt++ ){
			DNALength prefixLength = (*keyIt).readPos * data->insRate;
			DNALength substringLength = (*data->reads)[(*keyIt).readIndex].length * data->insRate;
			DNALength substringPos;
			if (genomePos < substringLength) {
				substringPos = 0;
			}
			else {
				substringPos = genomePos - prefixLength;
			}
			//
			// Do not bother aligning the read again if it aligns to the same position.
			//
  		if ((*data->prevAlignedGenomePos)[(*keyIt).readIndex] == substringPos)
				continue;

			if (substringPos + substringLength > data->genome->length) {
				substringLength = data->genome->length - substringPos;
			}

			genomeSubstring.seq = &data->genome->seq[substringPos];
			genomeSubstring.length = substringLength;
			FastqAlignment alignment;
			int readIndex = (*keyIt).readIndex;
			int alignScore;
			alignScore = KBandAlign((*data->reads)[readIndex], genomeSubstring, SMRTDistanceMatrix, 
															6, // ins
															6, // del
															0.30*(*data->reads)[readIndex].length,
															insScoreMat, insPathMat,
															alignment, distanceMatrixScoreFn, QueryFit);
			if (alignScore < (*data->readOptScore)[readIndex]) {
				(*data->readOptScore)[readIndex] = alignScore;
				(*data->optAlignment)[readIndex] = alignment;
				(*data->optAlignment)[readIndex].tAlignedSeqPos = substringPos;
				(*data->optGenomeAlignPos)[readIndex] = substringPos;
				(*data->optGenomeAlignLength)[readIndex] = substringLength;
			}
			(*data->prevAlignedGenomePos)[readIndex] = substringPos;
			/*
				cout << genomePos << " read: " << readIndex 
				<< " readpos: " << (*keyIt).readPos << " score " << alignScore << endl;
			*/
		}
		if (genomePos % 1000 == 0) {
			cerr << genomePos << endl;
		}
  }
}

int main(int argc, char* argv[]) {
	
	string genomeFileName, readsFileName;
	TupleMetrics tm;
	float insRate = 0.10;
	tm.tupleSize = 8;
	CommandLineParser clp;
	int nProcessors = 1;
	clp.SetProgramName("exhalign");
	clp.SetProgramSummary("Count the number of occurrences of every k-mer in a file.");
	clp.RegisterStringOption("genome", &genomeFileName, "The file of the genome to align to.");
	clp.RegisterStringOption("reads",  &readsFileName,  "The reads to align.");
	clp.RegisterPreviousFlagsAsHidden();
	clp.RegisterIntOption("wordsize", &tm.tupleSize, "Size of words to count", 
												CommandLineParser::NonNegativeInteger);
	clp.RegisterFloatOption("insrate", &insRate, "Roughly the insertion rate (10%)", 
													CommandLineParser::NonNegativeFloat);
	clp.RegisterIntOption("nProc", &nProcessors, "Number of processors to use", CommandLineParser::NonNegativeInteger);
	clp.ParseCommandLine(argc, argv);

	insRate+=1.0;
	//
	// Process the reads into a vector of read keywords
	//
	

	vector<string> readsFileNames;
	vector<FASTQSequence> reads;
	vector<vector<ReadKeyword> > keywords;
	SMRTSequence seq, seqRC;
	ReadKeyword keyword;
	int readIndex = 0;

	if (FileOfFileNames::IsFOFN(readsFileName)) {
		FileOfFileNames::FOFNToList(readsFileName, readsFileNames);
	}
	else {
		readsFileNames.push_back(readsFileName);
	}

	ReaderAgglomerate genomeReader;	
	HDFRegionTableReader regionTableReader;
	genomeReader.Initialize(genomeFileName);
	FASTQSequence genome;
	genomeReader.GetNext(genome);
	SubreadIterator subreadIterator;

	keywords.resize(nProcessors);
	RegionTable  regionTable, *regionTablePtr;

	int readsFileIndex;
	for (readsFileIndex = 0; readsFileIndex < readsFileNames.size(); readsFileIndex++ ) {
		
		ReaderAgglomerate reader;
		reader.Initialize(readsFileNames[readsFileIndex]);
		regionTalePtr = NULL;
		
		if (reader.fileType == HDFPulse or
				reader.fileType == HDFBase) {
			regionTableReader.Initialize(readsFileNames[readsFileIndex]);
			regionTableReader.Read(regionTable);
			regionTablePtr = &regionTable;
		}
		else {
			regionTablePtr = NULL;
		}
		SMRTSequence fullSequence;
		while(reader.GetNext(fullSequence)) {

			subreadIterator.Initialize(&fullSequence, regionTablePtr);
			
			SMRTSequence seq;
			while (subreadIterator.GetNext(seq)) {
				DNALength pos;
				if (seq.length < tm.tupleSize) 
					continue;
				reads.push_back(seq);
				for (pos = 0; pos < seq.length - tm.tupleSize + 1; pos++) {
					keyword.tuple.FromStringLR(&seq.seq[pos], tm);
					keyword.readPos = pos;
					keyword.readIndex = readIndex;
					keywords[(readIndex/2)%nProcessors].push_back(keyword);
				}
				readIndex++;
				seq.MakeRC(seqRC);
				reads.push_back(seqRC);
				for (pos = 0; pos < seqRC.length - tm.tupleSize + 1; pos++) {
					keyword.tuple.FromStringLR(&seqRC.seq[pos], tm);
					keyword.readPos = pos;
					keyword.readIndex = readIndex;
					keywords[(readIndex/2)%nProcessors].push_back(keyword);
				}
				readIndex++;
				//				seq.Free();
				seqRC.Free();
			}
			fullSequence.Free();
		}
	}
	int procIndex;
	for (procIndex = 0; procIndex < nProcessors; procIndex++) {
		std::sort(keywords[procIndex].begin(), keywords[procIndex].end());
	}


  std::vector<int> prevAlignedGenomePos;
  std::vector<int> readOptScore;
  std::vector<FastqAlignment > optAlignment;
	std::vector<int> optGenomeAlignPos;
	std::vector<int> optGenomeAlignLength;

  prevAlignedGenomePos.resize(reads.size());
  readOptScore.resize(reads.size());
  optAlignment.resize(reads.size());
	optGenomeAlignPos.resize(reads.size());
	optGenomeAlignLength.resize(reads.size());
	vector<Data> tdata;
	tdata.resize(nProcessors);
  std::fill(prevAlignedGenomePos.begin(), prevAlignedGenomePos.end(), -1);
	for (procIndex = 0; procIndex < nProcessors; procIndex++) {
		tdata[procIndex].prevAlignedGenomePos = &prevAlignedGenomePos;
		tdata[procIndex].readOptScore         = &readOptScore;
		tdata[procIndex].optAlignment         = &optAlignment;
		tdata[procIndex].optGenomeAlignPos    = &optGenomeAlignPos;
		tdata[procIndex].optGenomeAlignLength = &optGenomeAlignLength;
		tdata[procIndex].keywords             = &keywords[procIndex];
		tdata[procIndex].genome               = &genome;
		tdata[procIndex].insRate              = insRate;
		tdata[procIndex].reads                = &reads;
		tdata[procIndex].tm                   = &tm;
	}
	if (nProcessors == 1) {
		KeywordSeededAlignment(&tdata[0]);
	}
	else {
		pthread_t *threads = new pthread_t[nProcessors];
		pthread_attr_t *threadAttr = new pthread_attr_t[nProcessors];
		for (procIndex = 0; procIndex < nProcessors; procIndex++) {
			pthread_attr_init(&threadAttr[procIndex]);			
			pthread_create(&threads[procIndex], &threadAttr[procIndex], (void*(*)(void*))KeywordSeededAlignment, &tdata[procIndex]);
		}
		for (procIndex = 0; procIndex < nProcessors; procIndex++) {
			pthread_join(threads[procIndex], NULL);
		}

	}
	VectorIndex i;
	//	cout << "printing alignments for " << reads.size() << " reads." << endl;
	for (readIndex = 0; readIndex < readOptScore.size(); readIndex +=2 ){
		int optIndex = readIndex;
		if (readOptScore[readIndex] > readOptScore[readIndex+1]) {
			optIndex= readIndex + 1;
		}
		FASTQSequence genomeSubstring;
		genomeSubstring.seq = &genome.seq[optGenomeAlignPos[optIndex]];
		genomeSubstring.length =  optGenomeAlignLength[optIndex];
		if (prevAlignedGenomePos[optIndex] >= 0) {
			optAlignment[optIndex].qName.assign(reads[optIndex].title, reads[optIndex].titleLength);
			optAlignment[optIndex].tName.assign(genome.GetName());
			ComputeAlignmentStats(optAlignment[optIndex], reads[optIndex].seq, genomeSubstring.seq, SMRTDistanceMatrix, 6, 6);
			if (optAlignment[optIndex].blocks.size() > 0) {
				PrintCompareSequencesAlignment(optAlignment[optIndex], reads[optIndex], genomeSubstring,cout);
			}
			/*			StickPrintAlignment(optAlignment[optIndex],
													reads[optIndex],
													genomeSubstring, cout, 0, optGenomeAlignPos[optIndex]);
			*/
						
		}
	}
	for (readIndex = 0; readIndex < readOptScore.size(); readIndex++ ) {
		reads[readIndex].Free();
	}

	return 0;
}
