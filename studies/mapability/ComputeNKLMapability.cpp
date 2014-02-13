#include "FASTASequence.h"
#include "FASTAReader.h"
#include "tuples/DNATupleList.h"
#include "algorithms/anchoring/GlobalChain.h"
#include "algorithms/metagenomics/FindRandomSequence.h"
#include "datastructures/suffixarray/SuffixArray.h"
#include "datastructures/suffixarray/SuffixArrayTypes.h"
#include "utils/FileOfFileNames.h"
#include <string>
#include <set>
#include <stdlib.h>
#include <semaphore.h>
#include "statistics/statutils.h"
#include "tuples/TupleMetrics.h"
#include <limits.h>
#include "Types.h"

sem_t largeMemorySemaphore;

bool ReadContainsN(DNASequence &read) {
	DNALength i;
	for (i = 0; i < read.length; i++ ){
		if (read.seq[i] == 'N') return true;
	}
	return false;
}

class CompareTuplePosWithTuple {
	int operator()(const PositionDNATuple &posTuple, const DNATuple &nonPosTuple) const {
		return posTuple.tuple < nonPosTuple.tuple;
	}
	int operator()(const DNATuple &nonPosTuple, const PositionDNATuple &posTuple) const {
		return nonPosTuple.tuple < posTuple.tuple;
	}
};

typedef vector<PositionDNATuple> PosTupleList;
typedef pair<PosTupleList::iterator, PosTupleList::iterator> PosTupleListBounds;
typedef vector<PosTupleList> PosTupleMatrix;



class Anchor : public PositionDNATuple{
public:
	static TupleMetrics *tm;
	UInt GetT() {
		return pos;
	}
	UInt GetQ() {
		return (UInt) readPos;
	}
	UInt GetW() {
		return tm->tupleSize;
	}
	int readPos;
	Anchor& operator=(const Anchor& rhs) {
		PositionDNATuple::operator=(rhs);
		readPos = rhs.readPos;
		return *this;
	}
};

TupleMetrics* Anchor::tm = NULL;


class Fragment {
public:
	PositionDNATuple tuple;
	vector<int> previous;
	int prevMax;
	int count;
	int readPos;
	Fragment() {
		count = 0;
	}
};




void CopyMaxChainAnchors(vector<Anchor> &anchors, vector<int> &maxChainIndices, int indexOffset,
												 int maxScore, int maxIndex, vector<Anchor> &maxChain) {
	// maxScore is the same as the number of anchors in the chain.
	maxChain.resize(maxScore);
	int i;
	int chainIndex = maxIndex;
	int maxChainIndex = maxScore-1;
	for (i = 0; i < maxScore; i++) {
		maxChain[maxChainIndex] = anchors[chainIndex + indexOffset];
		--maxChainIndex;
		//
		// max chain indices contains the index of the previous anchor on 
		// this chain.  When at the first in the chain, chainIndex ==
		// maxChainInddices[chainIndex].
		//
		chainIndex = maxChainIndices[chainIndex];
	}
}
		
									
int FindMaxPath(vector<Anchor> &anchors,
								int start, int end, 
								int tupleSize, float fracThreshold, vector<int> &maxPath, int &globalMaxScore, int &globalMaxIndex) {
	int listLength = end - start;
	maxPath.resize(listLength);
	vector<int> maxScore;
	maxScore.resize(listLength);
	std::fill(maxScore.begin(), maxScore.end(), 0);
	int i, j;
	globalMaxScore = 0;
	globalMaxIndex = 0;
	//
	// Anchors are sorted in order of the genome.  For each anchor in
	// the genome, look at all other anchors that land before it to see
	// if they are in the valid alignment area.
	// 
	for (i = 0; i < listLength; i++ ) { 
		int maxPrevScore = 0;
		int maxPrevIndex = i;
		for (j = 0; j < i; j++ ) {
			//
			// Are anchors i and j in order?
			//
			float readDist, genomeDist;
			if ((readDist =  anchors[i+start].readPos - (anchors[j+start].readPos + tupleSize)) >= 0  and
					(genomeDist = anchors[i+start].pos - (anchors[j+start].pos + tupleSize)) >= 0) {
				if (min(readDist,genomeDist) / max(readDist,genomeDist) > fracThreshold) {
					if (maxScore[j] > maxPrevScore) {
						maxPrevScore = maxScore[j];
						maxPrevIndex = j;
					}
				}
			}
		}
		maxScore[i] = maxPrevScore + 1;
		maxPath[i]  = maxPrevIndex;
		if (maxScore[i] > globalMaxScore) {
			globalMaxScore = maxScore[i];
			globalMaxIndex = i;
		}
	}
	return globalMaxScore;
}



void BuildTuplePosLists(DNASequence &seq, DNALength genomePosition, DNASequence &genome,
			TupleMetrics &tm, 
			DNASuffixArray &sa,
			vector<Anchor> &tupleMatches, bool &lockedMemorySemaphore, int maxMatches, int semValue, int procIndex, int seqIndex, int simSeqPos) {
	int numPositions = 0;
	DNALength seqPos;
	if (seq.length < tm.tupleSize) {
		tupleMatches.clear();
		return;
	}
	PositionDNATuple queryTuple;
	int tooLarge = 0;
	for (seqPos = 0; !tooLarge and seqPos < seq.length - tm.tupleSize + 1; seqPos++) {
		SAIndex low, high;
		SAIndex matchLength;
		matchLength = sa.StoreLCPBounds(genome.seq, genome.length, &seq.seq[seqPos], tm.tupleSize, low, high);

		//		assert(matchLength == tm.tupleSize);
		if (matchLength == tm.tupleSize) {
			SAIndex sai;
			Anchor tuple;
			tuple.FromStringLR(&seq.seq[seqPos], tm);
			tuple.readPos = seqPos;
			if (maxMatches == 0 or high - low < maxMatches) {
				for (sai = low; !tooLarge and sai < high; sai++ ){
					//
					// Don't add tuples corresponding to the read.
					//
					if (sa.index[sai] != genomePosition + seqPos) {
						tuple.pos = sa.index[sai];
						tupleMatches.push_back(tuple);
						if (genome.seq[sa.index[sai]] == 'N') {
						  cout << "ERROR, for some reason the match includes a 'N'" <<endl;
						  cout << "chr: " << seqIndex << " pos " << simSeqPos << endl;
						  assert(genome.seq[sa.index[sai]] != 'N');
						}
						if (tupleMatches.size() == 1000000 and semValue > 0) {
							//
							// large memory, wait on a semaphore to proceed.
							//
							if (!lockedMemorySemaphore) {
								sem_wait(&largeMemorySemaphore);
								lockedMemorySemaphore = true;
							}
						}
						if (tupleMatches.size() * sizeof(Anchor) > 100000000) {
						  tooLarge = true;
						  break;
						}
						assert(tupleMatches[tupleMatches.size()-1].readPos < seq.length);
					}
				}
			}
		}
	}
	sort(tupleMatches.begin(), tupleMatches.end(), OrderPositionDNATuplesByPosition());
}

class NKLData {
public:
	vector<DNASuffixArray> sa;
	vector<FASTASequence>  genome;
	int n, L, k;
	int L1, L2, k1, k2;
	float floatFraction;
	ostream *out;
	int max;
	int nproc;
	int semLimit;
  int procIndex;
};

sem_t outputSemaphore;

void	FindNKLMatchesFull(vector<FASTASequence> &genome, vector<DNASuffixArray> &sa, int n, int L1, int L2, int k1, int k2, float floatFraction, ostream *out, int semValue, int maxNAnchors=0, int procIndex=0) {
	TupleList<PositionDNATuple> tupleList;
	DNALength pos;
	
	DNASequence read;
	int i = 0;
	DNALength maxPos = 0;
	PosTupleList mappingList;
	OrderPositionDNATuplesByTuple orderByTuple;
	PosTupleList readTupleList;
	vector<int> maxChainIndices;		
	vector<Anchor> tupleMatches;
	set<vector<UInt> > optChainSet;
	vector<UInt> optChainIndices;
	vector<UInt> scores;
	vector<UInt> prevOpt;
	bool lockedMemorySemaphore = false;
	int numSamples = 0;
	while (numSamples < n) {
		int nIts = 0;
		UInt seqIndex, pos;
		FindRandomPos(genome, seqIndex, pos, L2);
		read.seq = &genome[seqIndex].seq[pos];
		int L, k;
		for (L = max(L1,k1); L <= L2; L+=100 ){ 
			read.length = L2;
			for (k = k2; k >= k1; k-=5) {
				map<int,int> kMapFreq;
				TupleMetrics tm;
				tm.tupleSize = k;
				Anchor::tm = &tm;
				readTupleList.clear();
				int p;
				Anchor tuple;
				for (p = 0; p < read.length - tm.tupleSize + 1; p++) {
					tuple.FromStringLR(&read.seq[p], tm);
					tuple.pos = p;
					readTupleList.push_back(tuple);
				}

				int refIndex;
				for (refIndex = 0; refIndex < genome.size(); refIndex++) {
					tupleMatches.clear();
					optChainSet.clear();
					BuildTuplePosLists(read, pos, genome[refIndex], tm, sa[refIndex], tupleMatches, lockedMemorySemaphore, maxNAnchors,semValue, procIndex, seqIndex,pos);

					DNALength genomeMatchIndex = 0;

					for (genomeMatchIndex = 0; genomeMatchIndex < tupleMatches.size(); genomeMatchIndex++) {
						DNALength matchEnd = genomeMatchIndex;
						while(matchEnd < tupleMatches.size() and
									tupleMatches[matchEnd].pos - tupleMatches[genomeMatchIndex].pos < read.length) {
							matchEnd++;
						}
						maxChainIndices.clear();
						int maxChainScore, maxChainIndex;
						optChainIndices.clear();
						maxChainScore = RestrictedGlobalChain(&tupleMatches[genomeMatchIndex], 
																									matchEnd - genomeMatchIndex, floatFraction, optChainIndices,
																									scores, prevOpt);
			
						//
						// Transform the indices into the matches into genome coordinate indices.
						//

						int i;
						DNALength nextMatchIndex = optChainIndices[optChainIndices.size()-1] + genomeMatchIndex;
						
						for (i = 0; i < optChainIndices.size(); i++ ){
							optChainIndices[i] = tupleMatches[optChainIndices[i] + genomeMatchIndex].GetT();
						}

						//
						// Check this chain for uniqueness from other chains that are found.
						//
						
						if (optChainIndices.size() > 0 and 
								optChainSet.find(optChainIndices) == optChainSet.end()) {
							optChainSet.insert(optChainIndices);
							kMapFreq[optChainIndices.size()]++;
						}

						set<vector<UInt> >::iterator setIt, nextIt;
						setIt = optChainSet.begin();
						while (setIt != optChainSet.end()) {
							// 
							// Assume we are not adding blank sets here.
							//
							assert((*setIt).size() != 0);
							int lastChainIndex = (*setIt).size() - 1;
							if ((*setIt)[lastChainIndex] < tupleMatches[genomeMatchIndex].GetT()) {
								nextIt = setIt;
								++nextIt;
								optChainSet.erase(setIt);
								setIt = nextIt;
							}
							else {
								++setIt;
							}
						}
						//
						// Advance past this match to save some time and not double-count hits.
						//
						if (nextMatchIndex > genomeMatchIndex) {
							genomeMatchIndex = nextMatchIndex;
						}
						//
						// Free up some of the opt chain set
						//
						
					}
				}

				if (lockedMemorySemaphore) {
					//
					// Do some memory clean up.
					//
					vector<Anchor>().swap(tupleMatches);
					sem_post(&largeMemorySemaphore);
					lockedMemorySemaphore= false;
				}
				map<int,int>::iterator freqIt, freqEnd;
				sem_wait(&outputSemaphore);
				*out << seqIndex << " " << pos << " " << L << " " << k << " ";
				for (freqIt = kMapFreq.begin(); freqIt != kMapFreq.end(); ++freqIt) {
					*out << freqIt->first << " " << freqIt->second << ", ";
				}
				*out << endl;
				sem_post(&outputSemaphore);
			}
		}
		++numSamples;
	}
}

void FindNKLMatches(NKLData *data) {
  FindNKLMatchesFull(data->genome, data->sa, data->n, data->L1, data->L2, data->k1, data->k2, data->floatFraction, data->out, data->semLimit, data->max, data->procIndex);
}	

void PrintUsage() {
		cout << "usage: computeNLKMapability genomeFileName suffixarrayFileNmame n l_start l_end k_start k_end f outFile"<<endl;
		cout << " n  Number of reads." << endl
				 << " l_start l_end  Length of read. Iterate by 100 from l_start to l_end" << endl
				 << " k_start k_end  K-mer sample size.  Iterate by 10." << endl
				 << " f  Fractional indel rate (0.20 for example)" << endl;
		cout << " -max   m   Do not add anchors if there are more than m at a position " << endl
				 << " -nproc n   Fork off n processes" << endl;
		cout << " -norandinit Do not do random number initialization" << endl
				 << " -sem n(5)  Place value of n on semaphore." << endl;
}

int main(int argc, char* argv[]) {

	if (argc < 8) {
		PrintUsage();
		exit(1);
	}
	int nProc = 1;
	int argi = 1;
	string genomeFileName = argv[argi++];
	string saFileName     = argv[argi++];
	NKLData data;
	int semValue = 0;
	data.max = 0;
	data.n = atoi(argv[argi++]);
	data.L1 = atoi(argv[argi++]);
	data.L2 = atoi(argv[argi++]);
	data.k1 = atoi(argv[argi++]);
	data.k2 = atoi(argv[argi++]);
	data.floatFraction = atof(argv[argi++]);
	string outFileName = argv[argi++];
	
	bool randInit = true;
	while (argi < argc) {
		if (strcmp(argv[argi], "-nproc") == 0){ 
			nProc = atoi(argv[++argi]);
		}
		else if (strcmp(argv[argi], "-norandinit") == 0){ 
			randInit = false;
		}
		else if (strcmp(argv[argi], "-sem") == 0) {
			semValue = atoi(argv[++argi]);
		}
		else if (strcmp(argv[argi], "-max") == 0) {
			data.max = atoi(argv[++argi]);
		}
		else {
			PrintUsage();
			exit(1);
		}
		++argi;
	}
	data.semLimit = semValue;
	ofstream outFile;
	CrucialOpen(outFileName, outFile, std::ios::out);
	data.out = &outFile;
	//	data.sa.Read(saFileName);
	vector<string> genomeFileNames;
	vector<string> saFileNames;
	if (randInit) {
		void InitializeRandomGeneratorWithTime();
	}

	if (FileOfFileNames::IsFOFN(genomeFileName)) {
		FileOfFileNames::FOFNToList(genomeFileName, genomeFileNames);
	}	
	else {
		genomeFileNames.push_back(genomeFileName);
	}

	if (FileOfFileNames::IsFOFN(saFileName)) {
		FileOfFileNames::FOFNToList(saFileName, saFileNames);
	}	
	else {
		saFileNames.push_back(saFileName);
	}
	
	data.genome.resize(genomeFileNames.size());
	assert(genomeFileNames.size() == saFileNames.size());
	data.sa.resize(saFileNames.size());
	int refIndex;
	for(refIndex = 0 ; refIndex < genomeFileNames.size(); refIndex++) {
		FASTAReader reader;
		reader.Initialize(genomeFileNames[refIndex]);
		reader.ReadAllSequencesIntoOne(data.genome[refIndex]);
		data.sa[refIndex].Read(saFileNames[refIndex]);
	}
	
	//	FASTAReader reader;
	//	reader.Init(genomeFileName);
	//	reader.ReadAllSequencesIntoOne(data.genome);
	if (sem_init(&outputSemaphore, 0,1) == -1) {
		cout << "ERROR, cout not initialize output semaphore." << endl;
		exit(1);
	}
	if (sem_init(&largeMemorySemaphore, 0, semValue) == -1) {
		cout << "ERROR, cout not initialize output semaphore." << endl;
		exit(1);
	}

	int procIndex;
	if (nProc ==1) {
		FindNKLMatches(&data);
	}
	else {
		pthread_t *threads = new pthread_t[nProc];
		pthread_attr_t *threadAttr = new pthread_attr_t[nProc];
		for (procIndex = 0; procIndex < nProc; procIndex++) {
			pthread_attr_init(&threadAttr[procIndex]);	
			data.procIndex = procIndex;
			pthread_create(&threads[procIndex], &threadAttr[procIndex], (void*(*)(void*))FindNKLMatches, &data);
			
		}
		for (procIndex = 0; procIndex < nProc; procIndex++) {
			pthread_join(threads[procIndex], NULL);
		}
	}
	
}
