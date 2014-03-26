#ifndef SAMSUPPLEMENTALQVLIST
#define SAMSUPPLEMENTALQVLIST
#include "SMRTSequence.h"
#include "datastructures/alignmentset/SAMQVConversion.h"

class SupplementalQVList {
 public:
	enum QVList {Insertion=0x1, Deletion=0x2, Substitution=0x4, Merge=0x8, SubstitutionTag=0x10, DeletionTag=0x20};
	enum QVIndex {I_Insertion=1,I_Deletion=2,I_Substitution=3,I_Merge=4,I_SubstitutionTag=5,I_DeletionTag=6};
	unsigned int useqv;
	void SetDefaultQV() {
		useqv = Insertion | Deletion | Substitution | Merge | SubstitutionTag | DeletionTag;
	}
	static const char* qvTags[];
	static const char* qvNames[];
	static int nqvTags;
	static int nTags;

	int UseQV(vector<string> &qvList) {
		int i;
		useqv = 0;
		for (i = 0; i < qvList.size(); i++) {
			int j;
			for (j = 0; j < nTags; j++) {
				if (qvList[i] == qvNames[j]) {
					useqv |= 1 << j;
					break;
				}
			}
			if (j == nTags) {
				return 1;
			}
		}
		return 0;
	}
	void FormatQVOptionalFields(SMRTSequence &alignedSubsequence) {
		int i;
		for (i = 0; i < nqvTags; i++) {
			if (alignedSubsequence.GetQVPointerByIndex(i+1)->data == NULL) {
				// mask off this quality value since it does not exist
				useqv = useqv & ~(1 << i);
			}
		}

		for (i = 0; i < nqvTags; i++) {
			if (useqv & (1 << i)) {
				QualityVectorToPrintable(alignedSubsequence.GetQVPointerByIndex(i+1)->data, alignedSubsequence.length);
			}
		}
	}
	
	void PrintQVOptionalFields(SMRTSequence &alignedSubsequence, ostream &out) {
		int i = 0;
		for (i = 0; i < nqvTags; i++) {
			if (alignedSubsequence.GetQVPointerByIndex(i+1)->data == NULL) {
				// mask off this quality value since it does not exist
				useqv = useqv & ~(1 << i);
			}
		}


		for (i = 0; i < nTags; i++) {
			if (alignedSubsequence.GetQVPointerByIndex(i) != NULL and (useqv & (1 << i)) ) {
				out << "\t" << qvTags[i] << ":Z:";
				alignedSubsequence.PrintAsciiRichQuality(out, i + 1, 0);
			}
		}

		if (alignedSubsequence.substitutionTag != NULL and (useqv & SubstitutionTag)) {
			out << "\t" << qvTags[I_SubstitutionTag-1] << ":Z:";
			alignedSubsequence.PrintAsciiRichQuality(out, I_SubstitutionTag, 0);
		}

		if (alignedSubsequence.deletionTag != NULL and (useqv & DeletionTag)) {
			out << "\t" << qvTags[I_DeletionTag-1] << ":Z:";
			alignedSubsequence.PrintAsciiRichQuality(out, I_DeletionTag, 0);
		}
		
	}
};

const char* SupplementalQVList::qvNames[] = {"Insertion", "Deletion", "Substitution", "Merge", "SubstitutionTag", "DeletionTag"};
const char* SupplementalQVList::qvTags[] = {"qi", "qd", "qs", "qm", "ts", "td"};

// Only the first 4 tags are quality values.
int SupplementalQVList::nqvTags = 4;
int SupplementalQVList::nTags = 6;


#endif
