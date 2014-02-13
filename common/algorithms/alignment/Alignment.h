#ifndef ALIGNMENT_TYPE_H_
#define ALIGNMENT_TYPE_H_

#include <vector>
using namespace std;
#include "../../datastructures/alignment/AlignmentMap.h"
#include "../../datastructures/alignment/AlignmentBlock.h"
#include "../../datastructures/alignment/Path.h"
#include "../../datastructures/alignment/AlignmentStats.h"

class Alignment : public AlignmentStats {
 public:
	int qPos, tPos;
	// The positions in every block are relative to qPos and tPos;
	vector<Block> blocks;
	int size() {
		return blocks.size();
	}

	void AppendAlignment(Alignment &next) {
		int qOffset = next.qPos - qPos;
		int tOffset = next.tPos - tPos;
		Block tempBlock;
		int n;
		for (n = 0; n < next.blocks.size(); n++ ) {
			tempBlock = next.blocks[n];
			tempBlock.qPos += qOffset;
			tempBlock.tPos += tOffset;
			blocks.push_back(tempBlock);
		}
	}

	void ArrowPathToAlignment(vector<Arrow> &optPath) {
		int q, t;
		int a = 1;
		q = 0; t = 0;
		Block b;
		a = 0;
		while (a < optPath.size()) {
			if (optPath[a] == Diagonal) {
				// Start of a block;
				b.qPos = q;
				b.tPos = t;
				b.length = 0;
				while(a < optPath.size() and optPath[a] == Diagonal) {
					b.length++;
					a++;
					t++;
					q++;
				}
				blocks.push_back(b);
			}
			while(a < optPath.size() and optPath[a] == Left) {
				t++;
				a++;
			}
			while(a < optPath.size() and optPath[a] == Up) {
				q++;
				a++;
			}
		}
	}

	void AlignmentToMap(AlignmentMap &map) {
		map.tPos = tPos;
		map.qPos = qPos;
		if (blocks.size() == 0) {
			map.alignPos.resize(0);
			return;
		}

		int lastQueryPos = blocks[blocks.size()- 1].qPos + blocks[blocks.size() - 1].length;
		map.alignPos.resize(lastQueryPos - qPos + 1);
		int b;
		std::fill(map.alignPos.begin(),
							map.alignPos.end(), 
							-1);
	}

	void ComputeStats(char *tSeq, char* qSeq) {
		int b, bp;
		int qp, tp;
		qp = qPos, tp=tPos;
		nMatch = nMismatch = nIns = nDel = 0;
		for (b = 0; b < blocks.size(); b++) {
			//
			// Count matches in a fixed-length block.
			//
			for (bp = 0; bp < blocks[b].length; bp++ ){
				if (TwoBit[tSeq[tp]] == TwoBit[qSeq[qp]]) {
					nMatch++;
				}
				else {
					nMismatch++;
				}
				tp++;
				qp++;
			}

			//
			// Count the matches in the random drift portion of a common gap.
			//
			if (b < blocks.size() - 1) {
				int qGap = blocks[b+1].qPos - (blocks[b].qPos + blocks[b].length);
				int tGap = blocks[b+1].tPos - (blocks[b].tPos + blocks[b].length);
				int commonGap = qGap;
				if (tGap < qGap) {
					commonGap = tGap;
				}
				for (bp = 0; bp < commonGap; bp++) {
					if (TwoBit[tSeq[tp]] == TwoBit[qSeq[qp]]) {
						nMatch++;
					}
					else {
						nMismatch++;
					}
					tp++;
					qp++;
				}
				nDel += tGap;
				nIns += qGap;
				tp += tGap;
				qp += qGap;
			}
		}
		if (tp + qp > 0) {
			pctSimilarity = (nMatch*2) / (tp + qp);
		}
		else {
			pctSimilarity = 0;
		}
	}
};


#endif
