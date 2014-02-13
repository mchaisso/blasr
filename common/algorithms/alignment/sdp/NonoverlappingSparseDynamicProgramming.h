#ifndef NONOVERLAPPING_SPARSE_DYNAMIC_PROGRAMMING_H_
#define NONOVERLAPPING_SPARSE_DYNAMIC_PROGRAMMING_H_

/*
 *
 * Compute the subset of fragmentSet of largest weight that is not overlapping.
 *
 */
#include "FragmentSort.h"
#include "SDPSet.h"
#include "SDPFragment.h"
#include "SDPColumn.h"
#include "../AlignmentUtils.h"
#include "../../../datastructures/alignment/Alignment.h"


template<typename T_Fragment>
int SDPHeaviestSubsequence(int queryLength,
													 vector<T_Fragment> &fragmentSet,
													 int indel, int match,
													 vector<int> &maxFragmentChain,
													 AlignmentType alignType=Global) {

	maxFragmentChain.clear();

	if (fragmentSet.size() < 1)
		return 0;
		
	std::sort(fragmentSet.begin(), fragmentSet.end(), LexicographicFragmentSort<T_Fragment>());

	SDPSet<SDPColumn> activeSet;

	int sweepX;
	int trailRow;
	int fSweep;
	int fi;
	for (fi = 0; fi < fragmentSet.size(); fi++) {
		fragmentSet[fi].index = fi;
	}

	sweepX = fragmentSet[0].x;

	fSweep = 0;
	int maxChainLength = 0;
	int maxChainFragment = -1;
	int minFragmentCost, minFragmentIndex;
	minFragmentCost = INF_INT;
	minFragmentIndex = -1;
	for (sweepX = 0; sweepX < queryLength; sweepX++) {
		//
		// For each fragment ending in sweepX, attempt to compute its cost
		// by looking for the best previous match.
		//
		int fSweepStart = fSweep;
		while (fSweep < fragmentSet.size() and 
					 fragmentSet[fSweep].x == sweepX) {

			//
			// Compute the cost of every fragment in the sweep.
			//
			int cp = INF_INT, cl = INF_INT, ca = INF_INT;
			SDPColumn curCol, predCol;
			curCol.col = fragmentSet[fSweep].y;
			//
			// Search preceeding fragments.
			//

			//
			// Compute the cost of fragment_f
			int foundPrev = 0;
			T_Fragment fragStartPoint;
			fragStartPoint.x = fragmentSet[fSweep].x - fragmentSet[fSweep].GetLength();
			fragStartPoint.y = fragmentSet[fSweep].y - fragmentSet[fSweep].GetLength();
			T_Fragment predFrag;
			
			int predFragScore;
			
			if (activeSet.Predecessor(fragStartPoint, predFrag)) {
				//
				// predCol points to the fragment with greatest value less than curCol.
				// 
				predFragScore = fragmentSet[predCol.optFragment].cost + 
					abs((fragmentSet[fSweep].x - fragmentSet[fSweep].y)
							- (fragmentSet[predCol.optFragment].x - fragmentSet[predCol.optFragment].y)) * indel;
				foundPrev = 1;
			}

			//
			// Now assign the score for the fragment.  When doing a local alignment
			// start a new chain if the score is not good.
			
			if (foundPrev and 
					(alignType == Global or
					 (alignType == Local and predFragScore < 0))) {				
				fragmentSet[fSweep].chainPrev = predCol.optFragment;
				fragmentSet[fSweep].cost = predFragScore - fragmentSet[fSweep].weight;
				fragmentSet[fSweep].chainLength = fragmentSet[fragmentSet[fSweep].chainPrev].chainLength + 1;
			}
			else {
				fragmentSet[fSweep].chainPrev = -1;
				fragmentSet[fSweep].cost = fragmentSet[fSweep].GetLength() * match - fragmentSet[fSweep].weight;
				fragmentSet[fSweep].chainLength = 1;
			}
			
			if (minFragmentCost > fragmentSet[fSweep].cost) {
				minFragmentCost  = fragmentSet[fSweep].cost;
				minFragmentIndex = fSweep;
				//	maxChainLength = fragmentSet[fSweep].chainLength;
			}
			
			if (fragmentSet[fSweep].chainLength > maxChainLength) {
				maxChainLength = fragmentSet[fSweep].chainLength;
				maxChainFragment = fSweep;
			}

			fSweep++;
		}
		
		//
		// Each column must contain the highest scoring hit in that column. 
		//
		fSweep = fSweepStart;
		while (fSweep < fragmentSet.size() and
					 fragmentSet[fSweep].x == sweepX) {
			//
			// These elements are removed from the sweep set since they are done being processed.
			// If they are the lowest cost in the value, update activeSet
			//
			SDPColumn col;
			int storeCol = 0;
			col.col = fragmentSet[fSweep].y;

			if (activeSet.Member(col)) {
				if (fragmentSet[col.optFragment].cost < fragmentSet[fSweep].cost) {
					storeCol = 1;
				}
			}
			else {
				storeCol = 1;
			}
			if (storeCol) {
				col.col = fragmentSet[fSweep].y;
				col.optFragment = fSweep;
				// 
				// Insert new column or replace col with a more optimal one.
				//
				activeSet.Insert(col);

				// 
				// The invariant structure of the activeSet is that
				// after inserting a fragment of score S at column col, 
				// the score of all columns greater than 'col' in col set
				// must be less than col. 
				//
				// To preserve this invariant, when an element is inserted
				// at 'col', look to columns greater.  As long as any columns
				// have scores that are greater than col, remove them.
				// Once a column col_next has been found that has a score less than S
				// by the structure of the loop invariant, all columns greater than col_next
				// are guaranteed to have lower score than S, so we can continue searching
				// through this loop.
				//
				// Since fragments are processed at most once, this remains O(M).
				
				SDPColumn successorCol = col;

				while (activeSet.Successor(col, successorCol) and
							 fragmentSet[successorCol.optFragment].cost > fragmentSet[fSweep].cost) {
					activeSet.Delete(successorCol);
				}
			}

			//
			// Now remove this fragment, it is at the end of the sweep line.
			//
			int deleted;
			deleted = activeSet.Delete(fragmentSet[fSweep]);
			assert(deleted);
			++fSweep;
		}
	}
	if (alignType == Local) {
		maxChainFragment = minFragmentIndex;
	}
	while (maxChainFragment != -1) {
		maxFragmentChain.push_back(maxChainFragment);
		maxChainFragment = fragmentSet[maxChainFragment].chainPrev;
	}
	std::reverse(maxFragmentChain.begin(), maxFragmentChain.end());
	return maxFragmentChain.size();
}

#endif
