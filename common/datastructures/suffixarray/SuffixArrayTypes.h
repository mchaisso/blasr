#ifndef SUFFIX_ARRAY_TYPES_H_
#define SUFFIX_ARRAY_TYPES_H_

#include "SuffixArray.h"
#include "SharedSuffixArray.h"

#include "../../cmpseq/CompressedSequence.h"
#include "../../algorithms/compare/Compare4BitCompressed.h"
#include "../../FASTASequence.h"

typedef SuffixArray<Nucleotide, vector<int> > DNASuffixArray;
typedef SuffixArray<Nucleotide, vector<int>, 
	                  Compare4BitCompressed<Nucleotide>,
	                  CompressedDNATuple<FASTASequence> >       CompressedDNASuffixArray;

#endif
