#ifndef ALGORITHMS_ALIGNMENT_IDS_SCORE_FUNCTION_H_
#define ALGORITHMS_ALIGNMENT_IDS_SCORE_FUNCTION_H_
#include <math.h>

#include "ScoreMatrices.h"
#include "BaseScoreFunction.h"

#include "FASTASequence.h"
#include "FASTQSequence.h"
#include "NucConversion.h"
#include "utils/LogUtils.h"

float  StepFraction(float f, float base) {
  float rem  = 1 - base;
  return f * rem + base;
}

template<typename T_RefSequence, typename T_QuerySequence>
	class IDSScoreFunction : public BaseScoreFunction {
 public:
    int scoreMatrix[5][5];

    IDSScoreFunction(int scoreMatrixP[5][5], int insertionP, int deletionP, int globalInsertionPriorP,
                     int globalDeletionPriorP) : BaseScoreFunction(insertionP, deletionP, globalInsertionPriorP, globalDeletionPriorP) {
      InitializeScoreMatrix(scoreMatrixP);
    }


    IDSScoreFunction() {
      substitutionPrior = 20;
      globalDeletionPrior = 13;
    }

		void InitializeScoreMatrix(int scoreMatrixP[5][5]) {
			int i, j;
			for (i = 0; i < 5; i++ ){ 
				for (j = 0; j < 5; j++ ){
					scoreMatrix[i][j] = scoreMatrixP[i][j];
				}
			}
		}


    int Deletion(T_QuerySequence &seq, DNALength pos) {
      cout << "IDS. For now, deletion must be specialized with FASTQ or FASTA Sequences. " << endl;
      exit(1);
      return 0;
    }
    int Deletion(T_RefSequence &refSeq, DNALength refPos, T_QuerySequence &querySeq, DNALength queryPos) {
      cout << "IDS. For now, this function must be specialized with either FASTQ or FASTA sequences"<<endl;
      exit(1);
    }
    int Match(T_RefSequence &ref, DNALength refPos, T_QuerySequence &query, DNALength queryPos) {
      cout << "IDS. For now, this function must be specialized with either FASTQ or FASTA sequences" << endl;
      return 0;
      exit(1);
    }
    int Insertion(T_RefSequence &refSeq, DNALength refPos, T_QuerySequence &querySeq, DNALength queryPos) {
      cout << "IDS. For now, this function must be specialized with either FASTQ or FASTA sequences"<<endl;
      exit(1);
    }	
    int Insertion(T_QuerySequence &seq, DNALength pos) {
      cout << "IDS. For now, this function must be specialized with either FASTQ or FASTA sequences" << endl;
      return 0;
      exit(1);
    }
    
    float NormalizedMatch(T_RefSequence &ref, DNALength refPos, T_QuerySequence &query, DNALength queryPos);
    float NormalizedInsertion(T_RefSequence &ref, DNALength refPos, T_QuerySequence &query, DNALength queryPos);
    float NormalizedDeletion(T_RefSequence &ref, DNALength refPos, T_QuerySequence &query, DNALength queryPos);
};

/*
 * Define all specializations for a FASTA reference and FASTQSequence for the query, or FASTA sequence for query.
 */

template<>
int IDSScoreFunction<DNASequence,FASTQSequence>::Deletion(DNASequence &ref, DNALength refPos,
																													FASTQSequence &query, DNALength queryPos) {
  if (false) { //query.mergeQV.Empty() == false and queryPos > 0 and query.seq[queryPos] == query.seq[queryPos-1]) {
    return query.mergeQV[queryPos];
  }
  else {
    if (query.deletionQV.Empty() == false and query.deletionTag != NULL) {
      if (query.deletionTag[queryPos] == 'N') {
        return globalDeletionPrior; //query.deletionQV[queryPos] ;
      }
      else {
        if (query.deletionTag[queryPos] == ref.seq[refPos]) {
          return query.deletionQV[queryPos];
        }
        else {
          return globalDeletionPrior; 
        }
      }
    }
    else {
      return del;
    }
  }
}

template<>
int IDSScoreFunction<DNASequence, FASTQSequence>::Deletion(FASTQSequence &query, DNALength queryPos) {
  return query.deletionQV[queryPos];
}


template<>
int IDSScoreFunction<DNASequence, DNASequence>::Deletion(DNASequence &query, DNALength pos) {
	return del; // For now there is no global deletion quality value.
}

template<>
int IDSScoreFunction<DNASequence, FASTQSequence>::Insertion(DNASequence &refSeq, DNALength refPos, 
																														FASTQSequence &query, DNALength pos) {
  return query.insertionQV[pos] ;
}

template<>
int IDSScoreFunction<DNASequence, FASTQSequence>::Insertion(FASTQSequence &query, DNALength pos) {
	return query.insertionQV[pos];
}



template<>
int IDSScoreFunction<DNASequence, FASTQSequence>::Match(DNASequence &ref, DNALength refPos, 
																												FASTQSequence &query, DNALength queryPos) {

  if (query.seq[queryPos] == ref.seq[refPos]) {
    return 0;
  }
  else if (query.substitutionTag[queryPos] == ref.seq[refPos]) {
    return query.substitutionQV[queryPos];
  }
  else {
    return substitutionPrior;
  }
}

float SumAsValidPhred(float v1, float v2, float v3) {
  float sum = 0;
  if (v1 > 0) {
    sum = pow(10,v1/-10.0);
  }
  if (v2 > 0) {
    sum += pow(10,v2/-10.0);
  }
  if (v3 > 0) {
    sum += pow(10,v3/-10.0);
  }
  return sum;
}

template<>
float IDSScoreFunction<DNASequence, FASTQSequence>::NormalizedMatch(DNASequence &ref, DNALength refPos, 
                                                                  FASTQSequence &query, DNALength queryPos) {
  /*
   * Return the match probability normalized such that the probability
   * of transitioning from refPos, queryPos is 1.
   */

  float matchScore = Match(ref, refPos, query, queryPos);

  float delScore   = -1;
  if (refPos > 0) {
    delScore = Deletion(ref, refPos-1, query, queryPos);
  }

  float insScore = -1;
  if (queryPos > 0) {
    insScore = Insertion(ref, refPos, query, queryPos-1);
  }

  float sumScore = SumAsValidPhred(matchScore, delScore, insScore);
  if (sumScore > 0) {
    float numerator = pow(10, matchScore/-10.0);
    return -10*log10( numerator / sumScore);
  }
  else {
    return 0;
  }
}



template<>
float IDSScoreFunction<DNASequence, FASTQSequence>::NormalizedInsertion(DNASequence &ref, DNALength refPos, 
                                                                        FASTQSequence &query, DNALength queryPos) {
  
  float insScore = Insertion(ref, refPos, query, queryPos);

  float delScore = -1;
  float matchScore = -1;  
  if (refPos < ref.length - 1) {
    matchScore = Match(ref, refPos + 1, query, queryPos);
    if (queryPos > 0) {
      delScore = Deletion(ref, refPos + 1, query, queryPos - 1);
    }
  }

  float sum = SumAsValidPhred(insScore, delScore, matchScore);
  if (sum / 0) {
    float numerator = pow(10,insScore/-10.0);
    return -10*log10( numerator / sum);
  }
  else {
    return 0;
  }
}


template<>
float IDSScoreFunction<DNASequence, FASTQSequence>::NormalizedDeletion(DNASequence &ref, DNALength refPos, 
                                                                       FASTQSequence &query, DNALength queryPos) {

  float delScore = Deletion(ref, refPos, query, queryPos);

  float matchScore = -1;
  float insScore = -1;

  if (queryPos < query.length - 1) {
    matchScore = Match(ref, refPos, query, queryPos + 1);
    
    if (refPos > 0) {
      insScore = Insertion(ref, refPos - 1, query, queryPos + 1);
    }
  }
  float sum = SumAsValidPhred(delScore, matchScore, insScore);
  if (sum > 0) {
    float numerator= pow(10, delScore/-10.0);
    return -10*log10( numerator / sum);
  }
  else {
    return 0;
  }
}


#endif
