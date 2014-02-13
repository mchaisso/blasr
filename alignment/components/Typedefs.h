#ifndef ALIGNMENT_COMPONENTS_TYPEDEFS_H_
#define ALIGNMENT_COMPONENTS_TYPEDEFS_H_

typedef SMRTSequence T_Sequence;
typedef FASTASequence T_GenomeSequence;
typedef DNASuffixArray T_SuffixArray;
typedef DNATuple T_Tuple;
typedef LISPValueWeightor<T_GenomeSequence, DNATuple, vector<ChainedMatchPos> >  PValueWeightor;
typedef LISSMatchFrequencyPValueWeightor<T_GenomeSequence, DNATuple, vector<ChainedMatchPos> >  MultiplicityPValueWeightor;
typedef MappingData<T_SuffixArray, T_GenomeSequence, T_Tuple> MappingIPC;


#endif
