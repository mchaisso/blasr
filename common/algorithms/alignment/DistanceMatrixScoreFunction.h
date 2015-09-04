#ifndef ALGORITHMS_ALIGNMENT_DISTANCE_MATRIX_SCORE_FUNCTION_H_
#define ALGORITHMS_ALIGNMENT_DISTANCE_MATRIX_SCORE_FUNCTION_H_

#include "ScoreMatrices.h"
#include "BaseScoreFunction.h"

#include "FASTASequence.h"
#include "FASTQSequence.h"


template<typename T_RefSequence, typename T_QuerySequence>
	class DistanceMatrixScoreFunction : public BaseScoreFunction {
 public:
 DistanceMatrixScoreFunction() {
   ins = 0;
   del = 0;
 }

 DistanceMatrixScoreFunction(int scoreMatrixP[5][5], int insertionP, int deletionP) {
   InitializeScoreMatrix(scoreMatrixP);
   ins = insertionP;
   del = deletionP;
 }

 int scoreMatrix[5][5];
 void InitializeScoreMatrix(int scoreMatrixP[5][5]) {
   int i, j;
   for (i = 0; i < 5; i++ ){ 
     for (j = 0; j < 5; j++ ){
       scoreMatrix[i][j] = scoreMatrixP[i][j];
     }
   }
 }
 int Deletion(T_RefSequence &seq, DNALength pos, T_QuerySequence &querySeq, DNALength queryPos) {
   return del;
 }
 int Insertion(T_RefSequence &seq, DNALength pos, T_QuerySequence &querySeq, DNALength queryPos) {
   return ins;
 }
 int Deletion(T_RefSequence &seq, DNALength pos) {
   return del;
 }
 int Match(T_RefSequence &ref, DNALength refPos, T_QuerySequence &query, DNALength queryPos) {
   return scoreMatrix[ThreeBit[ref[refPos]]][ThreeBit[query[queryPos]]];		
 }

 //
 // Define the score function on dereferenced pointers for speed.
 //
 int Match(Nucleotide ref, Nucleotide query) {
   return scoreMatrix[ThreeBit[ref]][ThreeBit[query]];
 }
											 
 int Insertion(T_QuerySequence &seq, DNALength pos) {
   return ins;
 }

 float NormalizedMatch(T_RefSequence &ref, DNALength refPos, T_QuerySequence &query, DNALength queryPos) {return 0;}
 float NormalizedInsertion(T_RefSequence &ref, DNALength refPos, T_QuerySequence &query, DNALength queryPos) {return 0;}
 float NormalizedDeletion(T_RefSequence &ref, DNALength refPos, T_QuerySequence &query, DNALength queryPos) {return 0;}
  
};



/*
 * Define all specializations for a FASTA reference and FASTQSequence for the query, or FASTA sequence for query.
 */

template<>
int DistanceMatrixScoreFunction<FASTASequence, FASTASequence>::Deletion(FASTASequence &ref, DNALength pos) {
	return del;
}

template<>
int DistanceMatrixScoreFunction<DNASequence, FASTQSequence>::Deletion(DNASequence &ref, DNALength pos) {
	return del;
}

template<>
int DistanceMatrixScoreFunction<FASTASequence, FASTASequence>::Insertion(FASTASequence &query, DNALength pos) {
	// positive value for quality value penalizes the alignment.
	return ins;
}

template<>
int DistanceMatrixScoreFunction<DNASequence, FASTQSequence>::Insertion(FASTQSequence &query, DNALength pos) {
	// positive value for quality value penalizes the alignment.
	return ins;
}

template<>
int DistanceMatrixScoreFunction<FASTASequence, FASTASequence>::Match(FASTASequence &ref, DNALength refPos, 
                                                                     FASTASequence &query, DNALength queryPos) {
	// positive value for quality value penalizes the alignment.
	return scoreMatrix[ThreeBit[query.seq[queryPos]]][ThreeBit[ref.seq[refPos]]];
}


template<>
int DistanceMatrixScoreFunction<DNASequence, FASTQSequence>::Match(DNASequence &ref, DNALength refPos, 
																																	 FASTQSequence &query, DNALength queryPos) {
	// positive value for quality value penalizes the alignment.
	return scoreMatrix[ThreeBit[query.seq[queryPos]]][ThreeBit[ref.seq[refPos]]];
}



#endif
