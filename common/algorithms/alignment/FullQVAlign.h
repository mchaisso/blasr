#ifndef ALGORITHMS_ALIGNMENT_FULL_QV_ALIGN_H_
#define ALGORITHMS_ALIGNMENT_FULL_QV_ALIGN_H_

#include "../../datastructures/matrix/Matrix.h"
#include "../../FASTQSequence.h"
#include "../../FASTASequence.h"

template<typename T_Query, typename T_Reference>
double FullQVAlign(T_Query       &query,
									T_Reference   &target,
									Matrix<double> &alignProb) {
	
	alignProb.Resize(query.length + 1, target.length + 1);
	alignProb.Initialize(0);

	DNALength q, t;

	if (query.length == 0 or target.length == 0) { return 0; }

	// Initialize boundaries of ins/del/match probability matrices.
	q = 0;
	VectorIndex numCols = target.length + 1;
	VectorIndex numRows = query.length + 1;
	alignProb[0][0] = 1;
	for (t = 1; t < numCols; t++ ) {
		// cannot match to a gap
		alignProb[0][t] = log(target.GetInsertionQV(t-1)) + alignProb[0][t-1];
	}	

	for (q = 1; q < numRows; q++) {
		alignProb[q][0] = log(query.GetInsertionQV(q-1))  + alignProb[q-1][0];
	}
	
	// Now compute probability of alignment with the Forward algorithm.
	for (q = 1; q < numRows; q++) {
		for (t = 1; t < numCols; t++) {
			// First compute p_ins[q,t] as transitions from match matrix

			double logMatchedPulseProb  = 0;
			double logInsertedPulseProb = 0;
			double logDeletedPulseProb  = 0;


			//
			// Use inefficient coding for now.
			//
			
			// Compute match, the bases are either the same, in which case
			// this is simply the product of the probabilities of the two
			// matches.  Otherwise, either one of the pulses may be correct,
			// and the probability is the union of the two cases.
			// 
			
			double matchedPulseProb;
			
			if (query.seq[q-1] == target.seq[t-1]) {
				matchedPulseProb = (1-query.GetSubstitutionQV(q-1)) * (1-target.GetSubstitutionQV(t-1));
			}
			else {
				matchedPulseProb = (query.GetSubstitutionQV(q-1)/3.0)*(1-target.GetSubstitutionQV(t-1)) +
					((1-query.GetSubstitutionQV(q-1)))*(target.GetSubstitutionQV(t-1)/3.0);
			}

			matchedPulseProb = exp(alignProb[q-1][t-1])*matchedPulseProb;
			// 
			// An insertion in the query can be either a normal extra base
			// in the query, or a deletion in the reference.
			//
			//			logInsertedPulseProb = uery.GetInsertionQV(q-1)) + alignProb[q-1][t];
			
			double insertedPulseProb = 0;
			if (target.GetDeletionTag(t-1) != 'N') {
				//
				// The target has a pulse that was not strong enough to call a
				// real incorporation.  For now assume that the weak pulse is
				// the previous nucleotide in the query.  So the likelihood of
				// the weak pulse is influenced by the likelihood of the
				// previous nucleotide in the query. 
				//
				// Also, we only consider the previous base to be a missed
				// weak pulse if the current base is a match. 
				//
	
				if (q > 1) {
					insertedPulseProb = 
						(target.GetPreBaseDeletionQV(t-1, query.seq[q-2]) *target.GetDeletionQV(t-1)
						 + query.GetInsertionQV(q-1)) 
						* exp(alignProb[q-1][t]);
				}
				else {
					//
					// There can be no pre-base deletion tag here (could probably be an assert statement).
					//
					insertedPulseProb = query.GetInsertionQV(q-1)  * exp(alignProb[q-1][t]);
				}
			}
			else {
				insertedPulseProb = (query.GetInsertionQV(q-1) + target.GetDeletionQV(t-1))*exp(alignProb[q-1][t]);
			}
			
			//
			// An insertion in the target may be either a normal extra base
			// in the target, or a deletion in the query.
			//			logDeletedPulseProb = target.GetInsertionQV(t-1)) + alignProb[q][t-1];
			double deletedPulseProb = 0;
			if (query.GetDeletionTag(q-1) != 'N') {
				if (t > 1) {
					deletedPulseProb = (query.GetPreBaseDeletionQV(q-1, target.seq[t-2]) * query.GetDeletionQV(q-1) 
															+ target.GetInsertionQV(t-1))*exp(alignProb[q][t-1]);
				}
				else {
					// There was a dropped pulse before this position, but nothing to align it to.  
					deletedPulseProb = target.GetInsertionQV(t-1) * exp(alignProb[q][t-1]);
				}
			}
			else {
				deletedPulseProb =  (target.GetInsertionQV(t-1) + query.GetDeletionQV(q-1)) *exp(alignProb[q][t-1]);
			}

			// Determine the total probability of reaching this position.
			/*			cout << "align prob " << q << " " << t << " " <<  logMatchedPulseProb << " " 
							<<  logInsertedPulseProb << " " <<  logDeletedPulseProb << endl;*/
			alignProb[q][t] = log(matchedPulseProb + insertedPulseProb + deletedPulseProb);
		}
	}
	double fullAlignProb = alignProb[numRows-1][numCols-1];
	alignProb.Free();
	return fullAlignProb;
};


#endif
