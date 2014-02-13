#ifndef ALGORITHMS_ALIGNMENT_QUALITY_VALUE_SCORE_FUNCTION_H_
#define ALGORITHMS_ALIGNMENT_QUALITY_VALUE_SCORE_FUNCTION_H_
#include "../../FASTASequence.h"
#include "../../FASTQSequence.h"
#include "../../NucConversion.h"
#include "ScoreMatrices.h"
#include "BaseScoreFunction.h"

template<typename T_RefSequence, typename T_QuerySequence>
	class QualityValueScoreFunction : public BaseScoreFunction {
 public:
	
		int Deletion(T_RefSequence &seq, DNALength refPos, T_QuerySequence &querySeq, DNALength queryPos) {
		cout << "For now, deletion must be specialized with FASTQ or FASTA Sequences. " << endl;
		exit(1);
		return 0;
		}
		
	int Deletion(T_RefSequence &seq, DNALength pos) {
		cout << "(2) For now, deletion must be specialized with FASTQ or FASTA Sequences. " << endl;
		exit(1);
		return 0;
	}
	int Match(T_RefSequence &ref, DNALength refPos, T_QuerySequence &query, DNALength queryPos) {
		cout << "Match For now, this function must be specialized with either FASTQ or FASTA sequences" << endl;
		exit(1);
		return 0;
	}

	int Insertion(T_RefSequence &ref, DNALength refPos, T_QuerySequence &seq, DNALength pos) {
		cout << "For now, this function must be specialized with either FASTQ or FASTA sequences" << endl;
		exit(1);
		return 0;
	}

	int Insertion(T_QuerySequence &seq, DNALength pos) {
		cout << "For now, this function must be specialized with either FASTQ or FASTA sequences" << endl;
		exit(1);
		return 0;
	}
};


/*
 * Define all specializations for a FASTA reference and FASTQSequence for the query, or FASTA sequence for query.
 */
template<>
int QualityValueScoreFunction<FASTASequence, FASTQSequence>::Deletion(FASTASequence &ref, DNALength pos) {
	return del; // For now there is no global deletion quality value.
}

template<>
int QualityValueScoreFunction<DNASequence, FASTQSequence>::Deletion(DNASequence &ref, DNALength pos) {
	return del; // For now there is no global deletion quality value.
}

template<>
int QualityValueScoreFunction<DNASequence, FASTQSequence>::Deletion(DNASequence &seq, DNALength refPos, FASTQSequence &querySeq, DNALength queryPos) {
	return Deletion(seq, refPos);
}



template<>
int QualityValueScoreFunction<DNASequence, FASTQSequence>::Insertion(FASTQSequence &query, DNALength pos) {
	// positive value for quality value penalizes the alignment.
	return query.qual[pos];
}

template<>
int QualityValueScoreFunction<DNASequence, FASTQSequence>::Insertion(DNASequence &ref, DNALength refPos,
																																		 FASTQSequence &query, DNALength pos) {
	// positive value for quality value penalizes the alignment.
	//	return query.qual[pos];
	//	return Insertion(query, pos);
	return ins;
	
}

template<>
int QualityValueScoreFunction<DNASequence, FASTQSequence>::Match(DNASequence &ref, DNALength refPos, 
																																 FASTQSequence &query, DNALength queryPos) {
	// positive value for quality value penalizes the alignment.
	return QVDistanceMatrix[ThreeBit[query.seq[queryPos]]][ThreeBit[ref.seq[refPos]]] * query.qual[queryPos];
}








#endif
