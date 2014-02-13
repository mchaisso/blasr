#include "datastructures/suffixarray/SuffixArray.h"
#include "FASTAReader.h"
#include "FASTASequence.h"
#include "datastructures/suffixarray/SuffixArrayTypes.h"
#include "datastructures/matrix/Matrix.h"
#include <iostream>
#include <string>
#include <vector>
#include <semaphore.h>


sem_t sem_thread_list;
sem_t sem_thread_available;
sem_t sem_result_matrix;
class IPCData {
public:
	DNASuffixArray *sarray;
	FASTASequence  *seq;
	int n;
	Matrix<float>  *matrix;
	vector<int> ns;
	int ni;
	int lMin, lMax;
	int threadIndex;
	vector<bool>   unique;
	vector<char>   *threadIsRunning;
};

using namespace std;

void SearchWordUniqueness(IPCData *data) {
	DNALength i;
	DNALength numUnique = 0;
	DNALength s;
	DNALength l;
	int n = data->n;
	DNASuffixArray &sarray = *data->sarray;
	FASTASequence  &seq    = *data->seq;
  Matrix<float>  &matrix = *data->matrix;
	vector<char>   &threadIsRunning = *data->threadIsRunning;

	int lMin, lMax, nMin;
	nMin = data->ns[data->ni];
	lMin = data->lMin;
	lMax = data->lMax;

	std::fill(data->unique.begin(), data->unique.end(), false);
	int prevNumUnique = 0;
	for (l = lMin; l <= lMax; l++) {
		cerr << "l: " << l << endl;
		for (i = 0; i < seq.length - l; ) {
			DNALength j = i + 1;
			//
			// This word is already unique at a shorter length
			//
			prevNumUnique = numUnique;
			if (data->unique[i] == true) {
				i++;
				continue;
			}
			
			if (seq.length - sarray.index[i] < l) {
				i++;
				continue;
			}
 			//
			// Look to see how many words past position i are the same as i for up to 'l' characters.
			//
				
			// 
			//  At this point, the suffix at position i is not l-unique,
			//  therefore the first l-1 characters are equal. 
			//
			DNALength i2 = i+1;
			while(i2 < seq.length -l) {
				DNALength l1, l2;
				l1 = l2 = l;
				if (seq.length - sarray.index[i] < l1) {
					l1 = seq.length - sarray.index[i];
				}
				if (seq.length - sarray.index[i2] < l2) {
					l2 = seq.length - sarray.index[i2];
  			}
				DNALength lMin = std::min(l1,l2);
				if (lMin < l or (seq.seq[sarray.index[i]+l-1] != seq.seq[sarray.index[i2]+l-1])) {
					break;
				}
				/*				if (strncmp((const char*) &seq.seq[sarray.index[i]], 
										(const char*) &seq.seq[sarray.index[i2]], lMin) != 0) {
					break;
					}*/
				i2++;
			}
			if (i2 - i <= n) {
				//
				// There are less than 'n' instances of words of length l from i to i2-1.  Count them all unique.
				//
				numUnique += i2-i;
				DNALength iu;
				// Mark these positions as unique since if they are unique for 'l' bases they will be unique for 'l+1'.
				for (iu = i; iu < i2; iu++) {
					data->unique[iu] = true;
				}
			}
			i = i2;
		}
		DNALength il;
		sem_wait(&sem_result_matrix);
		matrix[data->ni][l-lMin] = numUnique / (1.0*seq.length);
		sem_post(&sem_result_matrix);
	}
	sem_wait(&sem_thread_list);
	threadIsRunning[data->threadIndex] = 0;
	sem_post(&sem_thread_list);
	sem_post(&sem_thread_available)		;
	pthread_exit(NULL); 
}


int main(int argc, char* argv[]) {
	string genomeFileName;
	string suffixArrayFileName;
	if (argc < 6) {
		cout << "Usage: countWordUniqueness genome suffixArray lenMin lenMax nproc n1 n2 ... nM"  << endl;
		exit(1);
	}
	int nMin, nMax, lMin, lMax;
	genomeFileName = argv[1];
	suffixArrayFileName = argv[2];
	lMin = atoi(argv[3]);
	lMax = atoi(argv[4]);
	int nProc= atoi(argv[5]);	
	int argi = 6;
	vector<int> ns;
	while(argi < argc) {
		ns.push_back(atoi(argv[argi]));
		argi++; 
	}
  cout << "checking mult <= " << ns[0] << endl;
	Matrix<float> matrix;

	// Get the ref sequence.
	FASTAReader reader;
	reader.Init(genomeFileName);
	FASTASequence seq;
	reader.ReadAllSequencesIntoOne(seq);
	seq.ToUpper();
	nMin = 0; nMax = ns.size() - 1;

	matrix.Resize(nMax-nMin+1, lMax - lMin + 1);
	DNASuffixArray sarray;
	sarray.Read(suffixArrayFileName);
	vector<IPCData> ipc;
	vector<int>  procSlots;
	procSlots.resize(nMax-nMin+1);
	std::fill(procSlots.begin(), procSlots.end(), -1);
	vector<pthread_t> threads;
	vector<pthread_attr_t> threadAttr;
	vector<char>      threadIsRunning;
	threadAttr.resize(nProc);
	threadIsRunning.resize(nProc);
	std::fill(threadIsRunning.begin(), threadIsRunning.end(), 0);

	ipc.resize(nProc);
	int i;
	seq.ToUpper();
	for (i = 0; i < nProc;i++) {
		ipc[i].sarray = &sarray;
		ipc[i].seq    = &seq;
		ipc[i].matrix  = &matrix;
		ipc[i].lMin   = lMin;
		ipc[i].lMax   = lMax;
		ipc[i].ns     = ns;

		ipc[i].threadIsRunning = &threadIsRunning;
		ipc[i].threadIndex = i;
		ipc[i].unique.resize(seq.length);
	}
	int ki;
	int n;
	int nComplete = 0;
	int nCounts   = nMax - nMin + 1;

	sem_init(&sem_thread_list, 0, 1);
	sem_init(&sem_thread_available, 0, 0);
	sem_init(&sem_result_matrix, 0, 1);
	threads.resize(nProc);
	int threadIndex = 0;
    cout << "nmin: " << nMin << " nmax " << nMax << endl;
	for (n = nMin; n <= nMax; n++ ){
		if (threadIndex < nProc) {
			ipc[threadIndex].n = ns[n];
			ipc[threadIndex].ni = n;
			pthread_attr_init(&threadAttr[threadIndex]);
			sem_wait(&sem_thread_list);
			threadIsRunning[threadIndex] = 1;
			sem_post(&sem_thread_list);
			pthread_create(&threads[threadIndex], &threadAttr[threadIndex], 
										 (void* (*)(void*))SearchWordUniqueness, &ipc[threadIndex]);			
			++threadIndex;
		}
		else {
			// Need to process more rows. Launch threads to do this.
			sem_wait(&sem_thread_available);
			
			// Wait for other threads to not be writing to the availability list
			sem_wait(&sem_thread_list);
			int ti;
			for (ti = 0; ti < nProc; ti++) {
				if (threadIsRunning[ti] == 0) {
					// Found a slot where a thread may be started. Do that.
					threadIsRunning[threadIndex] = 1;
					ipc[ti].n = n;
					sem_post(&sem_thread_list);
					pthread_create(&threads[ti], &threadAttr[ti], 
												 (void* (*)(void*))SearchWordUniqueness, &ipc[ti]);			
					break;
				}					
			}
			if (ti >= nProc) {
				cout << "ERROR! A thread should be open but no slots were found!" << endl;
			}
		}
	}
	//
	// Wait or all tasks to finish.
	DNALength l;
	for (threadIndex = 0; threadIndex < min(nMax-nMin+1, nProc); threadIndex++) {
		pthread_join(threads[threadIndex], NULL);
	}
	
	for (n = nMin; n <= nMax; n++) {
		for (l = lMin; l <= lMax; l++) {
			cout << matrix[n-nMin][l-lMin] << ", ";
		}
		cout << endl;
	}
}
