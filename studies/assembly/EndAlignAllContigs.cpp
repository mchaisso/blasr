#include <string>
#include <vector>
#include <iostream>
#include "FASTAReader.h"
#include "FASTASequence.h"
#include "algorithms/alignment/SWAlign.h"
#include "algorithms/alignment/ScoreMatrices.h"
#include "algorithms/alignment/printers/StickAlignmentPrinter.h"
#include "datastructures/matrix/Matrix.h"
#include "utils.h"

using namespace std;


int main(int argc, char* argv[]) {


	FASTAReader reader;
	FASTASequence read;
	int minLength = 100;
	int maxAlignLength = 500;
	if (argc < 3) {
		cout << "usage: endAlignAllContigs inFile outMatrix [-minLength L( " << minLength << "]" 
				 << " [-maxAlignLength L ] " << endl;
		exit(1);
	}
	string readsFileName, matrixFileName;
	readsFileName = argv[1];
	matrixFileName = argv[2];
	int argi = 3;
	while (argi < argc) {
		if (strcmp(argv[argi], "-minLength") == 0) {
			minLength = atoi(argv[++argi]);
		}
		if (strcmp(argv[argi], "-maxAlignLength") == 0) {
			maxAlignLength = atoi(argv[++argi]);
		}
		++argi;
	}
	vector<FASTASequence *> reads;
	vector<FASTASequence *> rcReads;
	reader.Init(readsFileName);
	ofstream matrixOut;
	CrucialOpen(matrixFileName, matrixOut);
	while(reader.GetNext(read)) {
		if (read.length > minLength) {
			FASTASequence *readCopy = new FASTASequence;
			(*readCopy) = read;
			reads.push_back(readCopy);
			readCopy    = new FASTASequence;
			read.MakeRC(*readCopy);
			rcReads.push_back(readCopy);
		}
	}

	Matrix<int> alignScores;
	Matrix<float> alignIdentities;
	alignScores.Resize(reads.size(), reads.size());
	alignIdentities.Resize(reads.size(), reads.size());
	vector<int> scoreMat;
	vector<Arrow> pathMat;
	int i, j;
	int alignScore;
	FASTASequence readi, readj;
	FASTASequence rcReadi, rcReadj;
	
	for (i = 0; i < reads.size(); i++) {
		float maxFrontIdent, maxEndIdent;
		int   maxFrontIdentIndex, maxEndIdentIndex;
		maxFrontIdent = 0; maxEndIdent = 0;
		maxFrontIdentIndex = 0;
		maxEndIdentIndex   = 0;
		int maxFrontIdentLength = 0;
		int maxEndIdentLength  = 0;
		int maxFrontLength     = 0;
		int	maxEndLength       = 0;
		int nmaxFrontLengthIndex = 0;
		int maxEndLengthIndex  = 0;
		float maxFrontLengthIdent = 0;
		float maxEndLengthIdent = 0;
		int maxFrontLengthIndex = 0;
		for (j = 0; j < reads.size(); j++ ){
			// 
			// Store the two ends of the alignment.
			//
			alignScore = 0;
			int rcAlignScore;
			Alignment alignment;
			Alignment rcAlignment;
			Alignment *optAlignment;
			if (j < i) {
				if (reads[i]->length < maxAlignLength) {
					readi.length = reads[i]->length;
					readi.seq    = &reads[i]->seq[0];
					rcReadi.length = rcReads[i]->length;
					rcReadi.seq    = &rcReads[i]->seq[0];
				}
				else {
					readi.length = maxAlignLength;
					readi.seq    = &reads[i]->seq[reads[i]->length - maxAlignLength];
					rcReadi.length = maxAlignLength;
					rcReadi.seq    = &rcReads[i]->seq[rcReads[i]->length - maxAlignLength];
				}				
				if (reads[j]->length < maxAlignLength) {
					readj.length = reads[j]->length;
					readj.seq    = &reads[j]->seq[0];
					rcReadj.length = rcReads[j]->length;
					rcReadj.seq    = &rcReads[j]->seq[0];
				}
				else {
					readj.length = maxAlignLength;
					readj.seq    = &reads[j]->seq[reads[j]->length - maxAlignLength];
					rcReadj.length = maxAlignLength;
					rcReadj.seq    = &rcReads[j]->seq[rcReads[j]->length - maxAlignLength];
				}				
				alignScore = SWAlign(readi, readj, SMRTDistanceMatrix, 3, scoreMat, pathMat, alignment, EndAnchored);
				rcAlignScore = SWAlign(rcReadi, rcReadj, SMRTDistanceMatrix, 3, scoreMat, pathMat, rcAlignment, EndAnchored);
			}
			
			if (j > i) {
				if (reads[i]->length < maxAlignLength) {
					readi.length = reads[i]->length;
					readi.seq    = &reads[i]->seq[0];
					rcReadi.length = rcReads[i]->length;
					rcReadi.seq    = &rcReads[i]->seq[0];

				}
				else {
					readi.length = maxAlignLength;
					readi.seq    = &reads[i]->seq[0];
					rcReadi.length = maxAlignLength;
					rcReadi.seq    = &rcReads[i]->seq[0];
				}				
				if (reads[j]->length < maxAlignLength) {
					readj.length = reads[j]->length;
					readj.seq    = &reads[j]->seq[0];
					rcReadj.length = rcReads[j]->length;
					rcReadj.seq    = &rcReads[j]->seq[0];
				}
				else {
					readj.length = maxAlignLength;
					readj.seq    = &reads[j]->seq[0];
					rcReadj.length = maxAlignLength;
					rcReadj.seq    = &rcReads[j]->seq[0];
				}				
				alignScore = SWAlign(readi, readj, SMRTDistanceMatrix, 3, scoreMat, pathMat, alignment, FrontAnchored);				
				rcAlignScore = SWAlign(rcReadi, rcReadj, SMRTDistanceMatrix, 3, scoreMat, pathMat, rcAlignment, FrontAnchored);
			}
			ComputeAlignmentStats(alignment, readi.seq, readj.seq, SMRTDistanceMatrix, 3,3 );
			ComputeAlignmentStats(rcAlignment, rcReadi.seq, rcReadj.seq, SMRTDistanceMatrix, 3,3 );

			if (alignScore < rcAlignScore){ 
				optAlignment = &alignment;
			}
			else {
				optAlignment = &rcAlignment;
			}
			int nBlocks = optAlignment->blocks.size();
			
			if (nBlocks > 0) {
				int qAlignedLength = optAlignment->blocks[nBlocks-1].qPos - optAlignment->blocks[0].qPos;
				/*				cout << alignScore << " " << optAlignment.pctSimilarity << " " << qAlignedLength << " "
						 << " " << i << " " << j << " " << readi.length << " " << readj.length << endl;
				*/
				alignIdentities[i][j] = optAlignment->pctSimilarity;
				//				StickPrintOptAlignment(optAlignment, readi, readj, cout, 0,0);
				if (j < i ) {
					if (optAlignment->pctSimilarity > maxEndIdent) {
						maxEndIdent = optAlignment->pctSimilarity;
						maxEndIdentIndex = j;
						maxEndIdentLength = qAlignedLength;
					}
					if (qAlignedLength > maxEndLength) {
						maxEndLength = qAlignedLength;
						maxEndLengthIdent = optAlignment->pctSimilarity;
						maxEndLengthIndex = j;
					}
				}
				else if (j > i) {
					if (optAlignment->pctSimilarity > maxFrontIdent) {
						maxFrontIdent       = optAlignment->pctSimilarity;
						maxFrontIdentIndex  = j;
						maxFrontIdentLength = qAlignedLength;
					}
					if (qAlignedLength > maxFrontLength) {
						maxFrontLength = qAlignedLength;
						maxFrontLengthIdent = optAlignment->pctSimilarity;
						maxFrontLengthIndex = j;
					}
				}
			}
			else {
				alignScore = 0;
				alignIdentities[i][j] = 0;				
			}
			alignScores[i][j] = alignScore;

		}
		cout << i << " " << reads[i]->title << " "
				 << " " << maxFrontIdent     << " " << maxFrontIdentLength << " " << maxFrontIdentIndex << " " << reads[maxFrontIdentIndex]->title << " "
				 << " " << maxEndIdent << " " << maxEndIdentLength << " "  << " " << maxEndIdentIndex << " " << reads[maxEndIdentIndex]->title << " "
				 << " " << maxFrontLength    << " " << maxFrontLengthIdent << " " << maxEndLengthIndex << " " << reads[maxEndLengthIndex]->title << " " 
				 << maxEndLength << " " << maxEndLengthIdent << " " << maxEndLengthIndex << " " << reads[maxEndLengthIndex]->title << endl;
	}

	alignIdentities.Print(matrixOut);
	return 0;
}
