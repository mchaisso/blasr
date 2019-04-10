#include <string>
#include <iostream>
#include <vector>
#include <set>
#include <sstream>
#include <pthread.h>
#include <stdlib.h>
#include <time.h>
#include <signal.h>
#include <execinfo.h>
#include <algorithm>

#include "MappingIPC.h"
#include "MappingSemaphores.h"

#include "CCSSequence.h"
#include "SMRTSequence.h"
#include "FASTASequence.h"
#include "FASTAReader.h"
#include "SeqUtils.h"
#include "defs.h"
#include "utils.h"

#include "tuples/DNATuple.h"
#include "tuples/HashedTupleList.h"
#include "algorithms/compare/CompareStrings.h"
#include "algorithms/alignment.h"
#include "algorithms/alignment/AffineKBandAlign.h"
#include "algorithms/alignment/GuidedAlign.h"
#include "algorithms/alignment/AffineGuidedAlign.h"
#include "algorithms/alignment/FullQVAlign.h"
#include "algorithms/alignment/ExtendAlign.h"
#include "algorithms/alignment/OneGapAlignment.h"
#include "algorithms/alignment/AlignmentUtils.h"
#include "algorithms/alignment/DistanceMatrixScoreFunction.h"
#include "algorithms/alignment/StringToScoreMatrix.h"
#include "algorithms/alignment/AlignmentFormats.h"
#include "algorithms/anchoring/LISPValue.h"
#include "algorithms/anchoring/LISPValueWeightor.h"
#include "algorithms/anchoring/LISSizeWeightor.h"
#include "algorithms/anchoring/LISQValueWeightor.h"
#include "algorithms/anchoring/FindMaxInterval.h" 
#include "algorithms/anchoring/MapBySuffixArray.h"
#include "algorithms/anchoring/ClusterProbability.h"
#include "algorithms/anchoring/BWTSearch.h"
#include "datastructures/metagenome/SequenceIndexDatabase.h"
#include "datastructures/metagenome/TitleTable.h"
#include "datastructures/suffixarray/SharedSuffixArray.h"
#include "datastructures/suffixarray/SuffixArrayTypes.h"
#include "datastructures/tuplelists/TupleCountTable.h" 
#include "datastructures/anchoring/WeightedInterval.h"
#include "datastructures/anchoring/AnchorParameters.h"
#include "datastructures/alignment/AlignmentCandidate.h"
#include "datastructures/alignment/AlignmentBlock.h"
#include "datastructures/alignment/AlignmentContext.h"

#include "datastructures/mapping/MappingMetrics.h"
#include "datastructures/reads/ReadInterval.h"
#include "utils/FileOfFileNames.h"
#include "utils/RegionUtils.h"
#include "qvs/QualityTransform.h"
#include "files/ReaderAgglomerate.h"
#include "files/CCSIterator.h"
#include "files/FragmentCCSIterator.h"
#include "data/hdf/HDFRegionTableReader.h"
#include "datastructures/bwt/BWT.h"
#include "datastructures/sequence/PackedDNASequence.h"
#include "CommandLineParser.h" 
#include "qvs/QualityValue.h"
#include "statistics/pdfs.h"
#include "statistics/cdfs.h"
#include "statistics/statutils.h"
#include "statistics/LookupAnchorDistribution.h"

#include "algorithms/anchoring/PiecewiseMatch.h"

#define MAX_PHRED_SCORE 254
#define MAPQV_END_ALIGN_WIGGLE 5

#ifdef USE_GOOGLE_PROFILER
#include "google/profiler.h"
#endif

using namespace std;

MappingSemaphores semaphores;

/*
 * Declare global structures that are shared between threads.
 */

ostream *outFilePtr;

static bool seen;
HDFRegionTableReader *regionTableReader;

typedef SMRTSequence T_Sequence;
typedef FASTASequence T_GenomeSequence;
typedef DNASuffixArray T_SuffixArray;
typedef DNATuple T_Tuple;

typedef LISPValueWeightor<T_GenomeSequence, DNATuple, vector<ChainedMatchPos> >  PValueWeightor;
//typedef LISSMatchFrequencyPValueWeightor<T_GenomeSequence, DNATuple, vector<ChainedMatchPos> >  MultiplicityPValueWeightor;

ReaderAgglomerate *reader;

typedef MappingData<T_SuffixArray, T_GenomeSequence, T_Tuple> MappingIPC;

long totalBases = 0;
int  totalReads = 0;

class ClusterInformation {
public:
  int maxClusterSize;
  float meanAnchorBasesPerRead;
  float sdAnchorBasesPerRead;
  int score;
  float pctSimilarity;
  int readLength;
  float nStdDev ;
  int numSignificant;
  int numForward, numReverse;
};


class ReadAlignments {
public:
  /*
    This class stores the alignments from a read.  A read may be
    aligned in several different modes:
    1. noSplitSureads - Treat the read as a unit from start to end
    2. subreads       - Align each subread independently
    3. denovo         - Only align the CCS sequence from a read
    4. allpass        - Align the de novo ccs sequences and then the
    subreads to where the denovo ccs aligned.
    5. fullpass       - Same as allpass, except using only complete
    subreads.
   
    The alignments are a raggad array of n sequences; n is 1 for cases 
    1 and 3, the number of subreads for cases 2 and 4, and the number
    of full length passes for case 5.

    A ReadAligments class must only have alignments for a single type
    of read in it.

  */

  vector<vector<T_AlignmentCandidate*> > subreadAlignments;
  vector<SMRTSequence> subreads;
  AlignMode alignMode;
  SMRTSequence read;
  int GetNAlignedSeq() {
    return subreadAlignments.size();
  }

  bool AllSubreadsHaveAlignments() {
    int i, nAlignedSeq;
    nAlignedSeq = subreadAlignments.size();
    for (i = 0; i < nAlignedSeq; i++) {
      if (subreadAlignments[i].size() == 0) {
        return false;
      }
    }
    return true;
  }

  void Clear() {
    int i;
    int nAlignedSeq;
    for (i = 0, nAlignedSeq = subreadAlignments.size(); i < nAlignedSeq; i++) {
      int nAlignments;
      int a;
      for (a = 0, nAlignments = subreadAlignments[i].size(); a < nAlignments; a++) {
        delete subreadAlignments[i][a];
      }
      subreadAlignments[i].clear();
    }

    for (i = 0, nAlignedSeq = subreads.size(); i< nAlignedSeq; i++) {
      subreads[i].FreeIfControlled();
      if (subreads[i].title != NULL) {
        delete[] subreads[i].title;
        subreads[i].title = NULL;
      }
    }

    subreadAlignments.clear();
    read.Free();
  }

  void Resize(int nSeq) {
    subreadAlignments.resize(nSeq);
    subreads.resize(nSeq);
  }
  
  void CheckSeqIndex(int seqIndex) {
    if ( seqIndex < 0 or seqIndex >= subreads.size() ) {
        cout << "ERROR, adding a sequence to an unallocated position." << endl;
        assert(0);
    }
  }

  void SetSequence(int seqIndex, SMRTSequence &seq) {
    CheckSeqIndex(seqIndex);
    subreads[seqIndex] = seq;
  }

  void AddAlignmentForSeq(int seqIndex, T_AlignmentCandidate *alignmentPtr) {
    CheckSeqIndex(seqIndex);
    subreadAlignments[seqIndex].push_back(alignmentPtr);
  }

  void AddAlignmentsForSeq(int seqIndex, vector<T_AlignmentCandidate*> &seqAlignmentPtrs) {
    CheckSeqIndex(seqIndex);
    subreadAlignments[seqIndex].insert(subreadAlignments[seqIndex].end(), seqAlignmentPtrs.begin(), seqAlignmentPtrs.end());
  }

  ~ReadAlignments() {
    read.FreeIfControlled();
  }
};


class SortAlignmentPointersByScore {
public:
  int operator()(T_AlignmentCandidate *lhs, T_AlignmentCandidate* rhs) {
    if (lhs->score == rhs->score) {
      return lhs->tPos + lhs->tAlignedSeqPos < rhs->tPos + rhs->tAlignedSeqPos;
    }
    else {
      return lhs->score < rhs->score;
    }
  }
};

class SortAlignmentPointersByMapQV {
public:
  int operator()(T_AlignmentCandidate *lhs, T_AlignmentCandidate* rhs) {
    if (lhs->mapQV == rhs->mapQV) {
      if (lhs->score == rhs->score) {
        return lhs->tPos + lhs->tAlignedSeqPos < rhs->tPos + rhs->tAlignedSeqPos;
      }
      else {
        return lhs->score < rhs->score;
      }
    }
    else {
      return lhs->mapQV > rhs->mapQV;
    }
  }
};

string GetMajorVersion() {
  return "1.MC.rc46";
}

void GetVersion(string &version) {
  string perforceVersionString("$Change$");
  version = GetMajorVersion();
  if (perforceVersionString.size() > 12) {
    version.insert(version.size(), ".");
    version.insert(version.size(), perforceVersionString, 9, perforceVersionString.size() - 11);
  }
}

void GetSoftwareVersion(string &changelistId, string &softwareVersion) {
	
	int nDots = 0;
	int i;
	softwareVersion = "";
	for (i = 0; i < changelistId.size(); i++) {
		if ( changelistId[i]== '.') {
			nDots +=1;
		}
		if (nDots == 2) {
			break;
		}
	}
	softwareVersion = changelistId.substr(0,i);
}

void MakeSAMHDString(string &hdString) {
  stringstream out;
  out << "@HD\t" << "VN:" << GetMajorVersion() << "\t" << "pb:3.0.1";
  hdString = out.str();
}

void MakeSAMPGString(string &commandLineString, string &pgString) {
  stringstream out;
  string version;
  GetVersion(version);
  out << "@PG\t" << "ID:BLASR\tVN:"<< version<<"\tCL:"<<commandLineString;
  pgString = out.str();
}


void ParseChipIdFromMovieName(string &movieName, string &chipId) {
  vector<string> movieNameFields;
  Tokenize(movieName, "_", movieNameFields);
  if (movieNameFields.size() == 1) {
    chipId = movieNameFields[0];
  }
  else if (movieNameFields.size() > 4) {
    chipId = movieNameFields[3];
  }
  else {
    chipId = "NO_CHIP_ID";
  }
}

//
// Define a list of buffers that are meant to grow to high-water
// marks, and not shrink down past that.   The memory is reused rather
// than having multiple calls to new.
//
class MappingBuffers {
public:
  vector<int> hpInsScoreMat, insScoreMat;
  vector<int> kbandScoreMat;
  vector<Arrow> hpInsPathMat, insPathMat;
  vector<Arrow> kbandPathMat;
  vector<int>   scoreMat;
  vector<Arrow> pathMat;
  vector<int>  affineScoreMat;
  vector<Arrow> affinePathMat;
  vector<ChainedMatchPos> matchPosList;
  vector<ChainedMatchPos> rcMatchPosList;
  vector<BasicEndpoint<ChainedMatchPos> > globalChainEndpointBuffer;
  vector<Fragment> sdpFragmentSet, sdpPrefixFragmentSet, sdpSuffixFragmentSet;
  TupleList<PositionDNATuple> sdpCachedTargetTupleList;
  TupleList<PositionDNATuple> sdpCachedTargetPrefixTupleList;
  TupleList<PositionDNATuple> sdpCachedTargetSuffixTupleList;
  std::vector<int> sdpCachedMaxFragmentChain;
  vector<double> probMat;
  vector<double> optPathProbMat;
  vector<float>  lnSubPValueMat;
  vector<float>  lnInsPValueMat;
  vector<float>  lnDelPValueMat;
  vector<float>  lnMatchPValueMat;

  void Reset() {
    vector<int>().swap(hpInsScoreMat);
    vector<int>().swap(insScoreMat);
    vector<int>().swap(kbandScoreMat);
    vector<Arrow>().swap(hpInsPathMat);
    vector<Arrow>().swap(insPathMat);
    vector<Arrow>().swap(kbandPathMat);
    vector<int>().swap(scoreMat);
    vector<Arrow>().swap(pathMat);
    vector<ChainedMatchPos>().swap(matchPosList);
    vector<ChainedMatchPos>().swap(rcMatchPosList);
    vector<BasicEndpoint<ChainedMatchPos> >().swap(globalChainEndpointBuffer);
    vector<Fragment>().swap(sdpFragmentSet);
    vector<Fragment>().swap(sdpPrefixFragmentSet);
    vector<Fragment>().swap(sdpSuffixFragmentSet);
    sdpCachedTargetTupleList.Reset();
    sdpCachedTargetPrefixTupleList.Reset();
    sdpCachedTargetSuffixTupleList.Reset();
    vector<int>().swap(sdpCachedMaxFragmentChain);
    vector<double>().swap(probMat);
    vector<double>().swap(optPathProbMat);
    vector<float>().swap(lnSubPValueMat);
    vector<float>().swap(lnInsPValueMat);
    vector<float>().swap(lnDelPValueMat);
    vector<float>().swap(lnMatchPValueMat);
  }
};

void SetHelp(string &str) {
  stringstream helpStream;
  helpStream << "   Options for blasr " << endl
             << "   Basic usage: 'blasr reads.{fasta,bas.h5} genome.fasta [-options] " << endl
             << "   option\tDescription (default_value)." << endl << endl
             << " Input Files." << endl
             << "   reads.fasta is a multi-fasta file of reads.  While any fasta file is valid input, " <<endl
             << "               it is preferable to use pls.h5 or bas.h5 files because they contain" << endl
             << "               more rich quality value information." << endl << endl
             << "   reads.bas.h5|reads.pls.h5 Is the native output format in Hierarchical Data Format of " <<endl
             << "               SMRT reads. This is the preferred input to blasr because rich quality" << endl
             << "               value (insertion,deletion, and substitution quality values) information is " << endl
             << "               maintained.  The extra quality information improves variant detection and mapping"<<endl
             << "               speed." << endl <<endl
             << "   -sa suffixArrayFile"<< endl
             << "               Use the suffix array 'sa' for detecting matches" << endl 
             << "               between the reads and the reference.  The suffix" << endl 
             << "               array has been prepared by the sawriter program." << endl << endl 
             << "   -ctab tab "<<endl
             << "               A table of tuple counts used to estimate match significance.  This is " << endl
             << "               by the program 'printTupleCountTable'.  While it is quick to generate on " << endl
             << "               the fly, if there are many invocations of blasr, it is useful to"<<endl
             << "               precompute the ctab." <<endl << endl
             << "   -regionTable table" << endl
             << "               Read in a read-region table in HDF format for masking portions of reads." << endl
             << "               This may be a single table if there is just one input file. " << endl
             << "               or a fofn.  When a region table is specified, any region table inside " << endl
             << "               the reads.bas.h5 or reads.bas.h5 files are ignored."<< endl
             << endl 
             << " Options for modifying reads. There is ancilliary information about substrings of reads " << endl
             << "               that is stored in a 'region table' for each read file.  Because " << endl
             << "               HDF is used, the region table may be part of the .bas.h5 or .pls.h5 file," << endl
             << "               or a separate file.  A contiguously read substring from the template " << endl
             << "               is a subread, and any read may contain multiple subreads. The boundaries " << endl
             << "               of the subreads may be inferred from the region table either directly or " <<endl
             << "               by definition of adapter boundaries.  Typically region tables also" << endl
             << "               contain information for the location of the high and low quality regions of"<<endl
             << "               reads.  Reads produced by spurrious reads from empty ZMWs have a high"<<endl
             << "               quality start coordinate equal to high quality end, making no usable read." <<endl
             << "   -useccs   " << endl
             << "               Align the circular consensus sequence (ccs), then report alignments" << endl
             << "               of the ccs subreads to the window that the ccs was mapped to.  Only " << endl
             << "               alignments of the subreads are reported." << endl  
             << "   -useccsall"<<endl 
             << "               Similar to -useccs, except all subreads are aligned, rather than just" << endl
             << "               the subreads used to call the ccs.  This will include reads that only"<<endl
             << "               cover part of the template." << endl
             << "   -useccsdenovo" << endl
             << "               Align the circular consensus, and report only the alignment of the ccs"<<endl
             << "               sequence." << endl
             << "   -noSplitSubreads (false)" <<endl
             << "               Do not split subreads at adapters.  This is typically only " << endl
             << "               useful when the genome in an unrolled version of a known template, and " << endl
             << "               contains template-adapter-reverse_template sequence." << endl
             << "   -ignoreRegions(false)" << endl
             << "               Ignore any information in the region table." << endl
             << "   -ignoreHQRegions (false)Ignore any hq regions in the region table." << endl
             << endl
             << " Alignment Output." << endl
             << "   -bestn n (10)" <<endl
			 << "               Report the top 'n' alignments." << endl
             << "   -sam        Write output in SAM format." << endl
             << "   -clipping [none|hard|soft|subread] (none)" << endl
             << "               Use no/hard/soft clipping for SAM output."<< endl
             << "   -out out (terminal)  " << endl
             << "               Write output to 'out'." << endl
             << "   -unaligned file" << endl
             << "               Output reads that are not aligned to 'file'" << endl
             << "   -m t           " << endl
             << "               If not printing SAM, modify the output of the alignment." << endl
             << "                t=" << StickPrint <<   " Print blast like output with |'s connecting matched nucleotides." << endl 
             << "                  " << SummaryPrint << " Print only a summary: score and pos." << endl 
             << "                  " << CompareXML <<   " Print in Compare.xml format." << endl 
             << "                  " << Vulgar <<       " Print in vulgar format (deprecated)." << endl
             << "                  " << Interval <<     " Print a longer tabular version of the alignment." << endl 
             << "                  " << CompareSequencesParsable  << " Print in a machine-parsable format that is read by compareSequences.py." << endl
             << "   -noSortRefinedAlignments (false) " << endl
             << "               Once candidate alignments are generated and scored via sparse dynamic "<< endl
             << "               programming, they are rescored using local alignment that accounts " << endl
             << "               for different error profiles." <<endl
             << "               Resorting based on the local alignment may change the order the hits are returned." << endl
             << "   -allowAdjacentIndels " << endl
             << "               When specified, adjacent insertion or deletions are allowed. Otherwise, adjacent " << endl
             << "               insertion and deletions are merged into one operation.  Using quality values " << endl
             << "               to guide pairwise alignments may dictate that the higher probability alignment "<<endl
             << "               contains adjacent insertions or deletions.  Current tools such as GATK do not permit" << endl
             << "               this and so they are not reported by default." << endl
             << "   -header" <<endl
             << "               Print a header as the first line of the output file describing the contents of each column."<<endl
             << "   -titleTable tab (NULL) " << endl
             << "               Construct a table of reference sequence titles.  The reference sequences are " << endl
             << "               enumerated by row, 0,1,...  The reference index is printed in alignment results" << endl
             << "               rather than the full reference name.  This makes output concise, particularly when" <<endl
             << "               very verbose titles exist in reference names."<<endl
             << "   -minPctIdentity p (0)"<< endl
             << "               Only report alignments if they are greater than p percent identity."<<endl
             << "   -minAlignLength l (0)"<< endl
             << "               Only report alignments if they are greater than l bases on the genome."<<endl
             << "   -unaligned file" << endl
             << "               Output reads that are not aligned to 'file'." << endl
             << endl 
             << " Options for anchoring alignment regions. This will have the greatest effect on speed and sensitivity." << endl
             << "   -minMatch m (10) " << endl
             << "               Minimum seed length.  Higher minMatch will speed up alignment, " << endl
             << "               but decrease sensitivity." << endl
             << "   -maxMatch l (inf)" << endl
             << "               Stop mapping a read to the genome when the lcp length reaches l.  " << endl
             << "               This is useful when the query is part of the reference, for example when " <<endl
             << "               constructing pairwise alignments for de novo assembly."<<endl
             << "   -maxLCPLength l (inf)" << endl
             << "               The same as -maxMatch." << endl
             << "   -maxAnchorsPerPosition m (inf) " << endl
             << "               Do not add anchors from a position if it matches to more than 'm' locations in the target." << endl
             << "   -advanceExactMatches E (0)" << endl
             << "               Another trick for speeding up alignments with match - E fewer anchors.  Rather than" << endl 
             << "               finding anchors between the read and the genome at every position in the read, " <<endl
             << "               when an anchor is found at position i in a read of length L, the next position " << endl
             << "               in a read to find an anchor is at i+L-E." << endl
             << "               Use this when alignining already assembled contigs." << endl
             << "   -nCandidates n (10)" << endl 
             << "               Keep up to 'n' candidates for the best alignment.  A large value of n will slow mapping" << endl
             << "               because the slower dynamic programming steps are applied to more clusters of anchors" <<endl
             << "               which can be a rate limiting step when reads are very long."<<endl
             << endl
             << "  Options for Refining Hits." << endl
             << "   -sdpTupleSize K (11)" << endl
             << "               Use matches of length K to speed dynamic programming alignments.  This controls" <<endl
             << "               accuracy of assigning gaps in pairwise alignments once a mapping has been found,"<<endl
             << "               rather than mapping sensitivity itself."<<endl
             << "   -scoreMatrix \"score matrix string\" " << endl
             << "               Specify an alternative score matrix for scoring fasta reads.  The matrix is " << endl
             << "               in the format " << endl
             << "                  ACGTN" << endl
             << "                A abcde" << endl
             << "                C fghij" << endl
             << "                G klmno" << endl
             << "                T pqrst" << endl
             << "                N uvwxy" << " . The values a...y should be input as a quoted space separated " << endl
             << "               string: \"a b c ... y\". Lower scores are better, so matches should be less " << endl
             << "               than mismatches e.g. a,g,m,s = -5 (match), mismatch = 6. " << endl
             << "   -affineOpen value (10) " << endl
             << "               Set the penalty for opening an affine alignment." << endl
             << "   -affineExtend value (0)" << endl
             << "               Change affine (extension) gap penalty. Lower value allows more gaps." << endl
		         << "   -alignContigs" << endl
						 << "               Configure mapping parameters to map very long (10s of Mbp) contigs against a reference." << endl << endl
             << " Options for overlap/dynamic programming alignments and pairwise overlap for de novo assembly. " << endl
             << "   -useQuality (false)" << endl
             << "               Use substitution/insertion/deletion/merge quality values to score gap and " << endl
             << "               mismatch penalties in pairwise alignments.  Because the insertion and deletion" << endl
             << "               rates are much higher than substitution, this will make many alignments " <<endl
             << "               favor an insertion/deletion over a substitution.  Naive consensus calling methods "<<endl
             << "               will then often miss substitution polymorphisms. This option should be " << endl
             << "               used when calling consensus using the Quiver method.  Furthermore, when " << endl
             << "               not using quality values to score alignments, there will be a lower consensus " << endl
             << "               accuracy in homolymer regions." << endl
             << "   -affineAlign (true)" << endl
             << "               Refine alignment using affine guided align." << endl << endl
             << " Options for filtering reads." << endl
             << "   -minReadLength l(50)" << endl
             << "               Skip reads that have a full length less than l. Subreads may be shorter." << endl 
             << "   -minSubreadLength l(0)" << endl
             << "               Do not align subreads of length less than l." << endl
             << "   -maxScore m (0)" << endl
             << "               Maximum score to output (high is bad, negative good)." << endl << endl
             << " Options for parallel alignment." << endl
             << "   -nproc N (1)" << endl
             << "               Align using N processes.  All large data structures such as the suffix array and " << endl
             << "               tuple count table are shared."<<endl
             << "   -start S (0)" << endl
             << "               Index of the first read to begin aligning. This is useful when multiple instances " << endl
             << "               are running on the same data, for example when on a multi-rack cluster."<<endl
             << "   -stride S (1)" << endl
             << "               Align one read every 'S' reads." << endl << endl
             << " Options for subsampling reads." << endl
             << "   -subsample (0)" << endl
             << "               Proportion of reads to randomly subsample (expressed as a decimal) and align." << endl
             << endl
//             << " Options for dynamic programming alignments. " << endl << endl
//             << "   -ignoreQuality" << endl
//             << "                 Ignore quality values when computing alignments (they still may be used." << endl 
//             << "                 when mapping)." << endl << endl
             << " -v            Print some verbose information." << endl 
             << " -V 2          Make verbosity more verbose.  Probably only useful for development." << endl
             << " -h            Print this help file." << endl << endl
             << "To cite BLASR, please use: Chaisson M.J., and Tesler G., Mapping " << endl
             << "single molecule sequencing reads using Basic Local Alignment with " << endl
             << "Successive Refinement (BLASR): Theory and Application, BMC " << endl
             << "Bioinformatics 2012, 13:238 ." << endl << endl;
  str = helpStream.str();
}

void SetConciseHelp(string &conciseHelp) {
  stringstream strm;
  strm << "blasr - a program to map reads to a genome" << endl
       << " usage: blasr reads genome " << endl
       << " Run with -h for a list of commands " << endl
       << "          -help for verbose discussion of how to run blasr." << endl;
  conciseHelp = strm.str();
}

void PrintDiscussion() {
  cout << "NAME"<<endl;
  cout << "         blasr - Map SMRT Sequences to a reference genome."<< endl << endl;
  cout << "SYNOPSIS" << endl
       << "         blasr reads.fasta genome.fasta " << endl << endl
       << "         blasr reads.fasta genome.fasta -sa genome.fasta.sa" << endl << endl
       << "         blasr reads.bas.h5 genome.fasta [-sa genome.fasta.sa] " << endl << endl
       << "         blasr reads.bas.h5 genome.fasta -sa genome.fasta.sa -maxScore -100 -minMatch 15 ... " << endl << endl
       << "         blasr reads.bas.h5 genome.fasta -sa genome.fasta.sa -nproc 24 -out alignment.out ... " << endl << endl
       << "DESCRIPTION " << endl
       << "  blasr is a read mapping program that maps reads to positions " << endl
       << "  in a genome by clustering short exact matches between the read and" << endl
       << "  the genome, and scoring clusters using alignment. The matches are" << endl
       << "  generated by searching all suffixes of a read against the genome" << endl
       << "  using a suffix array. Global chaining methods are used to score " << endl
       << "  clusters of matches." << endl << endl
       << "  The only required inputs to blasr are a file of reads and a" << endl
       << "  reference genome.  It is exremely useful to have read filtering" << endl
       << "  information, and mapping runtime may decrease substantially when a" << endl
       << "  precomputed suffix array index on the reference sequence is" << endl
       << "  specified." << endl
       << "  " << endl
       << "  Although reads may be input in FASTA format, the recommended input is HDF" << endl
       << "  bas.h5 and pls.h5 files because these contain qualtiy value" << endl
       << "  information that is used in the alignment and produces higher quality" << endl
       << "  variant detection.  " << endl
       << "  " << endl
       << "  Read filtering information is contained in the .bas.h5 input files as" << endl
       << "  well as generated by other post-processing programs with analysis of" << endl
       << "  pulse files and read in from a separate .region.h5 file.  The current" << endl
       << "  set of filters that are applied to reads are high quality region" << endl
       << "  filtering, and adapter filtering.  Regions outside high-quality" << endl
       << "  regions are ignored in mapping.  Reads that contain regions annotated" << endl
       << "  as adapter are split into non-adapter (template) regions, and mapped" << endl
       << "  separately." << endl
       << "  " << endl
       << "  When suffix array index of a genome is not specified, the suffix array is" << endl
       << "  built before producing alignment.   This may be prohibitively slow" << endl
       << "  when the genome is large (e.g. Human).  It is best to precompute the" << endl
       << "  suffix array of a genome using the program sawriter, and then specify" << endl
       << "  the suffix array on the command line using -sa genome.fa.sa." << endl
       << "  " << endl
       << "  The optional parameters are roughly divided into three categories:" << endl
       << "  control over anchoring, alignment scoring, and output. " << endl
       << "  " << endl
       << "  The default anchoring parameters are optimal for small genomes and" << endl
       << "  samples with up to 5% divergence from the reference genome.  The main" << endl
       << "  parameter governing speed and sensitivity is the -minMatch parameter." << endl
       << "  For human genome alignments, a value of 11 or higher is recommended.  " << endl
       << "  Several methods may be used to speed up alignments, at the expense of" << endl
       << "  possibly decreasing sensitivity.  " << endl
       << "  " << endl
       << "  Regions that are too repetitive may be ignored during mapping by" << endl
       << "  limiting the number of positions a read maps to with the" << endl
       << "  -maxAnchorsPerPosition option.  Values between 500 and 1000 are effective" << endl
       << "  in the human genome." << endl
	   << "  " << endl
       << "  For small genomes such as bacterial genomes or BACs, the default parameters " << endl
       << "  are sufficient for maximal sensitivity and good speed." << endl
       << endl << endl;
}

int CountZero(unsigned char *ptr, int length) {
  int i;
  int nZero = 0;
  for (i = 0; i < length; i++) {
    if (ptr[i] == 0) { ++nZero; }
  }
  return nZero;
}


bool ReadHasMeaningfulQualityValues(FASTQSequence &sequence) {
  if (sequence.qual.Empty() == true) {
    return 0;
  }
  else {
    int q;
    int numZero=0, numNonZero=0;
    if (sequence.qual.data == NULL) {
      return false;
    }
    numZero = CountZero(sequence.qual.data, sequence.length);
    numNonZero = sequence.length - numZero;
    int subNumZero = 0, subNonZero = 0;

    if (sequence.substitutionQV.data == NULL) {
      return false;
    }
    subNumZero = CountZero(sequence.substitutionQV.data, sequence.length);
    subNonZero = sequence.length - subNumZero;

    if (numZero < 0.5*numNonZero and subNumZero < 0.5 * subNonZero) {
       return true;
     }
    else {
      return false;
    }
  }
}

template<typename T_RefSequence, typename T_Sequence>
void PairwiseLocalAlign(T_Sequence &qSeq, T_RefSequence &tSeq, 
                        int k, 
                        MappingParameters &params, T_AlignmentCandidate &alignment, 
                        MappingBuffers &mappingBuffers,
                        AlignmentType alignType=Global) {
  //
  // Perform a pairwise alignment between qSeq and tSeq, but choose
  // the pairwise alignment method based on the parameters.  The
  // options for pairwise alignment are:
  //  - Affine KBanded alignment: usually used for sequences with no
  //                              quality information.
  //  - KBanded alignment: For sequences with quality information.
  //                       Gaps are scored with quality values.
  //  
  DistanceMatrixScoreFunction<DNASequence, FASTQSequence> distScoreFn;
	params.InitializeScoreFunction(distScoreFn);
  int kbandScore;
  int qvAwareScore;
  if (params.ignoreQualities || qSeq.qual.Empty() || !ReadHasMeaningfulQualityValues(qSeq) ) {

    kbandScore = AffineKBandAlign(qSeq, tSeq, SMRTDistanceMatrix, 
                                  params.indel+2, params.indel - 3, // homopolymer insertion open and extend
                                  params.indel+2, params.indel - 1, // any insertion open and extend
                                  params.indel, // deletion
                                  k*1.2,
                                  mappingBuffers.scoreMat, mappingBuffers.pathMat, 
                                  mappingBuffers.hpInsScoreMat, mappingBuffers.hpInsPathMat,
                                  mappingBuffers.insScoreMat, mappingBuffers.insPathMat,
                                  alignment, Global);

    alignment.score = kbandScore;
    if (params.verbosity >= 2) {
      cout << "align score: " << kbandScore << endl;
    }
  }
  else {

        
    if (qSeq.insertionQV.Empty() == false) {
      qvAwareScore = KBandAlign(qSeq, tSeq, SMRTDistanceMatrix, 
                                params.indel+2, // ins
                                params.indel+2, // del
                                k,
                                mappingBuffers.scoreMat, mappingBuffers.pathMat,
                                alignment, distScoreFn, alignType);
      if (params.verbosity >= 2) {
        cout << "ids score fn score: " << qvAwareScore << endl;
      }
    }
    else {
      qvAwareScore = KBandAlign(qSeq, tSeq, SMRTDistanceMatrix, 
                                params.indel+2, // ins
                                params.indel+2, // del
                                k,
                                mappingBuffers.scoreMat, mappingBuffers.pathMat,
                                alignment, distScoreFn, alignType);
      if (params.verbosity >= 2) {
        cout << "qv score fn score: " << qvAwareScore << endl;
      }
    }
    alignment.sumQVScore = qvAwareScore;
    alignment.score = qvAwareScore;
    alignment.probScore = 0;
  }
  // Compute stats and assign a default alignment score using an edit distance.
  ComputeAlignmentStats(alignment, qSeq.seq, tSeq.seq, distScoreFn);

  if (params.scoreType == 1) {
    alignment.score = alignment.sumQVScore;
  }
  
}


int SDPAlignLite(DNASequence &query, DNASequence &target, int wordSize, int maxMatchesPerPosition, 
									vector<Fragment> &fragmentSet ) {


	TupleList<PositionDNATuple> targetTupleList;

	/*
		Collect a set of matching fragments between query and target.
		Since this function is an inner-loop for alignment, anything to
		speed it up will help.  One way to speed it up is to re-use the
		vectors that contain the sdp matches. 
	*/

	TupleMetrics tm, tmSmall;
	tm.Initialize(wordSize);

	SequenceToTupleList(target, tm, targetTupleList);

  targetTupleList.Sort();
			 
  //
  // Store in fragmentSet the tuples that match between the target
  // and query.
  //

	StoreMatchingPositions(query, tm, targetTupleList, fragmentSet, maxMatchesPerPosition);
	VectorIndex f;  
	for (f = 0; f < fragmentSet.size(); f++) {
		fragmentSet[f].weight = tm.tupleSize;
		fragmentSet[f].length = tm.tupleSize;

  }
	
	std::sort(fragmentSet.begin(), fragmentSet.end(), LexicographicFragmentSort<Fragment>());
	//
	// Assume inversion will be in rc max frament chain set.
	//
	std::vector<int> maxFragmentChain;
	SDPLongestCommonSubsequence(query.length, fragmentSet, tm.tupleSize, 2, 2, -5, maxFragmentChain, Local);
  //GlobalChainFragmnetList(fragmentSet, maxFragmentChain);

	for (f = 0; f < maxFragmentChain.size(); f++) {
		fragmentSet[f] = fragmentSet[maxFragmentChain[f]];
	}
	fragmentSet.resize(f);
	return maxFragmentChain.size();

 }

template<typename T_RefSequence, typename T_Sequence>
void RefineAlignment(T_Sequence &query,
                     T_RefSequence &genome,
                     T_AlignmentCandidate  &alignmentCandidate, MappingParameters &params,
                     MappingBuffers &mappingBuffers) {


  FASTQSequence qSeq;
  DNASequence   tSeq;
  DistanceMatrixScoreFunction<DNASequence, FASTQSequence> distScoreFn;
	params.InitializeScoreFunction(distScoreFn);
  distScoreFn.InitializeScoreMatrix(SMRTDistanceMatrix);

  if (params.doGlobalAlignment) {
    SMRTSequence subread;
    ((FASTQSequence*)&subread)->ReferenceSubstring(query,
																									 query.subreadStart,
																									 query.subreadEnd - query.subreadStart);

    int drift = ComputeDrift(alignmentCandidate);
    T_AlignmentCandidate refinedAlignment;

    KBandAlign(subread, alignmentCandidate.tAlignedSeq, SMRTDistanceMatrix, 
               params.insertion, params.deletion,
                              drift,
                              mappingBuffers.scoreMat, mappingBuffers.pathMat,
                              refinedAlignment, distScoreFn, Global);
    refinedAlignment.RemoveEndGaps();
    
    alignmentCandidate.blocks = refinedAlignment.blocks;
    alignmentCandidate.gaps   = refinedAlignment.gaps;
    alignmentCandidate.tPos   = refinedAlignment.tPos;
    alignmentCandidate.qPos   = refinedAlignment.qPos + query.subreadStart;
    alignmentCandidate.score  = refinedAlignment.score;
  }
  else if (params.useGuidedAlign) {
    T_AlignmentCandidate refinedAlignment;
    int lastBlock = alignmentCandidate.blocks.size() - 1;
    

    if (alignmentCandidate.blocks.size() > 0) {

      /*
       * Refine the alignment without expanding past the current
       * boundaries of the sequences that are already aligned.
       */

      //
      // NOTE** this only makes sense when 
      // alignmentCandidate.blocks[0].tPos == 0. Otherwise the length
      // of the sequence is not correct.
      //
      tSeq.Copy(alignmentCandidate.tAlignedSeq, 
                alignmentCandidate.tPos,
                (alignmentCandidate.blocks[lastBlock].tPos + 
                 alignmentCandidate.blocks[lastBlock].length));
    
      //      qSeq.ReferenceSubstring(alignmentCandidate.qAlignedSeq,
      qSeq.ReferenceSubstring(query,
                              alignmentCandidate.qAlignedSeqPos + alignmentCandidate.qPos, 
                              (alignmentCandidate.blocks[lastBlock].qPos +
                               alignmentCandidate.blocks[lastBlock].length));


			if (params.affineAlign) {
				AffineGuidedAlign(qSeq, tSeq, alignmentCandidate, 
													distScoreFn, params.bandSize,
													mappingBuffers, 
													refinedAlignment, Global, false);
			}
			else {
        GuidedAlign(qSeq, tSeq, alignmentCandidate, 
                    distScoreFn, params.guidedAlignBandSize,
                    mappingBuffers,
                    refinedAlignment, Global, false);
			}

      ComputeAlignmentStats(refinedAlignment, 
                            qSeq.seq,
                            tSeq.seq, 
                            distScoreFn, params.affineAlign);

      //
      // Copy the refine alignment, which may be a subsequence of the
      // alignmentCandidate into the alignment candidate.  
      //

      // First copy the alignment block and gap (the description of
      // the base by base alignment).

      alignmentCandidate.blocks.clear();
      alignmentCandidate.blocks = refinedAlignment.blocks;

      alignmentCandidate.CopyStats(refinedAlignment);

      alignmentCandidate.gaps   = refinedAlignment.gaps;
      alignmentCandidate.score  = refinedAlignment.score;
      alignmentCandidate.nCells = refinedAlignment.nCells;

      // Next copy the information that describes what interval was
      // aligned.  Since the reference sequences of the alignment
      // candidate have been modified, they are reassigned.
      alignmentCandidate.tAlignedSeq.TakeOwnership(tSeq);
      alignmentCandidate.ReassignQSequence(qSeq);
      alignmentCandidate.tAlignedSeqPos    += alignmentCandidate.tPos; 
      alignmentCandidate.qAlignedSeqPos    += alignmentCandidate.qPos;

      //
      // tPos and qPos are the positions within the interval where the
      // alignment begins. The refined alignment has adifferent tPos
      // and qPos from the alignment candidate.
      alignmentCandidate.tPos = refinedAlignment.tPos;
      alignmentCandidate.qPos = refinedAlignment.qPos;

      // The lengths of the newly aligned sequences may differ, update those.
      alignmentCandidate.tAlignedSeqLength = tSeq.length;
      alignmentCandidate.qAlignedSeqLength = qSeq.length;
    }
  }
  else {


    //
    // This assumes an SDP alignment has been performed to create 'alignmentCandidate'. 
  
    //
    // Recompute the alignment using a banded smith waterman to
    // get rid of any spurious effects of usign the seeded gaps.
    //

    //
    // The k-banded alignment is over a subsequence of the first
    // (sparse dynamic programming, SDP) alignment.  The SDP
    // alignment is over a large window that may contain the
    // candidate sequence.  The k-band alignment is over a tighter
    // region.  

    int drift = ComputeDrift(alignmentCandidate);
          
    //
    // Rescore the alignment with a banded alignment that has a
    // better model of sequencing error.
    //

    if (alignmentCandidate.blocks.size() == 0 ){ 
      alignmentCandidate.score = 0;
      return;
    }
    int lastBlock = alignmentCandidate.blocks.size() - 1;

    //
    // Assign the sequences that are going to be realigned using
    // banded alignment.  The SDP alignment does not give that great
    // of a score, but it does do a good job at finding a backbone
    // alignment that closely defines the sequence that is aligned.
    // Reassign the subsequences for alignment with a tight bound
    // around the beginning and ending of each sequence, so that
    // global banded alignment may be performed.
    //
  
    //
    // This section needs to be cleaned up substantially.  Right now it
    // copies a substring from the ref to a temp, then from the temp
    // back to the ref.  It may be possible to just keep one pointer per
    // read to the memory that was allocated, then allow the seq
    // parameter to float around.  The reason for all the copying is
    // that in case there is a compressed version of the genome the
    // seqences must be transformed before alignment.
    //

    if (alignmentCandidate.qIsSubstring) {
      qSeq.ReferenceSubstring(query,  // the original sequence
                              alignmentCandidate.qPos + alignmentCandidate.qAlignedSeqPos, 
                              alignmentCandidate.blocks[lastBlock].qPos + alignmentCandidate.blocks[lastBlock].length);
    }
    else {
      qSeq.ReferenceSubstring(alignmentCandidate.qAlignedSeq, // the subsequence that the alignment points to
                              alignmentCandidate.qPos  + alignmentCandidate.qAlignedSeqPos, 
                              alignmentCandidate.blocks[lastBlock].qPos + alignmentCandidate.blocks[lastBlock].length - alignmentCandidate.blocks[0].qPos);
    }
      
    tSeq.Copy(alignmentCandidate.tAlignedSeq, // the subsequence the alignment points to
              alignmentCandidate.tPos, // ofset into the subsequence
              alignmentCandidate.blocks[lastBlock].tPos + alignmentCandidate.blocks[lastBlock].length - alignmentCandidate.blocks[0].tPos);

    T_AlignmentCandidate refinedAlignment;

    //
    // When the parameter bandSize is 0, set the alignment band size
    // to the drift off the diagonal, plus a little more for wiggle
    // room.  When the parameteris nonzero, use that as a fixed band.
    //
    int k;
    if (params.bandSize == 0) {
      k = abs(drift) * 1.5;
    }
    else {
      k = params.bandSize;
    }
    if (params.verbosity > 0) {
      cout << "drift: " << drift << " qlen: " << alignmentCandidate.qAlignedSeq.length << " tlen: " << alignmentCandidate.tAlignedSeq.length << " k: " << k << endl;
      cout << "aligning in " << k << " * " << alignmentCandidate.tAlignedSeq.length << " " << k * alignmentCandidate.tAlignedSeq.length << endl;
    }
    if (k < 10) {
      k = 10;
    }

    alignmentCandidate.tAlignedSeqPos    += alignmentCandidate.tPos; 
    
    VectorIndex lastSDPBlock = alignmentCandidate.blocks.size() - 1;

    if (alignmentCandidate.blocks.size() > 0) {
      DNALength prevLength =  alignmentCandidate.tAlignedSeqLength -= alignmentCandidate.tPos;
      alignmentCandidate.tAlignedSeqLength = (alignmentCandidate.blocks[lastSDPBlock].tPos 
                                              + alignmentCandidate.blocks[lastSDPBlock].length 
                                              - alignmentCandidate.blocks[0].tPos);
    }
    else {
      alignmentCandidate.tAlignedSeqLength = 0;
    }

    alignmentCandidate.tPos              = 0;
    alignmentCandidate.qAlignedSeqPos    += alignmentCandidate.qPos;

    if (alignmentCandidate.blocks.size() > 0) {
      DNALength prevLength =  alignmentCandidate.qAlignedSeqLength -= alignmentCandidate.qPos; 
      alignmentCandidate.qAlignedSeqLength = (alignmentCandidate.blocks[lastSDPBlock].qPos 
                                              + alignmentCandidate.blocks[lastSDPBlock].length
                                              - alignmentCandidate.blocks[0].qPos);
    }
    else {
      alignmentCandidate.qAlignedSeqLength = 0;
    }
    alignmentCandidate.qPos                = 0;

    alignmentCandidate.blocks.clear();
    alignmentCandidate.tAlignedSeq.Free();
    alignmentCandidate.tAlignedSeq.TakeOwnership(tSeq);
    alignmentCandidate.ReassignQSequence(qSeq);

    if (params.verbosity >= 2) {
      cout << "refining target: " << endl;
      alignmentCandidate.tAlignedSeq.PrintSeq(cout);
      cout << "refining query: " << endl;
      ((DNASequence*)&alignmentCandidate.qAlignedSeq)->PrintSeq(cout);
      cout << endl;
    }
    PairwiseLocalAlign(qSeq, tSeq, k, params, alignmentCandidate, mappingBuffers, Fit);
  }
}


int AlignSubstring(DNASequence &tSeq, int tPos, int tLength, 
									 FASTQSequence &qSeq, int qPos, int qLength, 
									 MappingParameters &params,
									 int sdpTupleSize,
									 DistanceMatrixScoreFunction<DNASequence, FASTQSequence> &distScoreFn,
									 MappingBuffers &refinementBuffers,
									 T_AlignmentCandidate &alignment) {

	DNASequence tSubSeq;
	FASTQSequence qSubSeq;
	tSubSeq.ReferenceSubstring(tSeq, tPos, tLength);
	qSubSeq.ReferenceSubstring(qSeq, qPos, qLength);

	int alignScore;
	if ((tLength < params.bandSize or qLength < params.bandSize) and 
			((tLength >= qLength and tLength < 1.5 * qLength) or
			 (tLength <= qLength and tLength > 0.5 * qLength))) {
		alignScore = AffineKBandAlign(qSubSeq, tSubSeq, distScoreFn.scoreMatrix, 
																	params.indel+2, params.indel - 3, // homopolymer insertion open and extend
																	params.indel+2, params.indel - 1, // any insertion open and extend
																	params.indel, // deletion
																	params.bandSize,
																	refinementBuffers.scoreMat, refinementBuffers.pathMat, 
																	refinementBuffers.hpInsScoreMat, refinementBuffers.hpInsPathMat,
																	refinementBuffers.insScoreMat, refinementBuffers.insPathMat,
																			alignment, Global);
	}
	else {
		AlignmentType alignType = Global;
		
		alignScore = SDPAlign(qSubSeq, tSubSeq, distScoreFn, sdpTupleSize, 
													params.sdpIns, params.sdpDel, 0.25, 
													alignment, refinementBuffers, alignType,
													params.detailedSDPAlignment, 
													true,
													params.sdpPrefix,
													params.recurse,
													params.recurseOver, 
													params.sdpMaxAnchorsPerPosition, Global);
		
		alignment.qAlignedSeq.ReferenceSubstring(qSubSeq);
		alignment.tAlignedSeq.ReferenceSubstring(tSubSeq);

		if (alignment.blocks.size() > 0) {
			int last = alignment.blocks.size()-1;
			if (alignment.blocks[last].tPos + alignment.blocks[last].length < tSubSeq.length and 
					alignment.blocks[last].qPos + alignment.blocks[last].length < qSubSeq.length) {
				Block end;
				end.tPos = tSubSeq.length - 1;
				end.qPos = qSubSeq.length - 1;
				end.length = 1;
				alignment.blocks.push_back(end);
			}
			
			if (alignment.blocks[0].tPos > 0 and alignment.blocks[0].qPos > 0) {
				//
				// Pin the alignment to the beginning of the sequences.
				//
				Block start;
				start.tPos = 0;
				start.qPos = 0;
				start.length = 1;
				alignment.blocks.insert(alignment.blocks.begin(), start);
			}
		}

		RefineAlignment(qSubSeq, tSubSeq, 
										alignment, 
										params, 
										refinementBuffers);
	}
	return alignScore;
}

void AppendSubseqAlignment(T_AlignmentCandidate &alignment, 
													 T_AlignmentCandidate &alignmentInGap, 
													 DNALength qPos, 
													 DNALength tPos) {
	if (alignmentInGap.blocks.size() > 0) {
		int b;
		//
		// Configure this block to be relative to the beginning
		// of the aligned substring.  
		//
		for (b = 0; b < alignmentInGap.size(); b++) {
			alignmentInGap.blocks[b].tPos += tPos + alignmentInGap.tPos + alignmentInGap.tAlignedSeqPos;
			alignmentInGap.blocks[b].qPos += qPos + alignmentInGap.qPos + alignmentInGap.qAlignedSeqPos;

			assert(alignmentInGap.blocks[b].tPos < alignment.tAlignedSeq.length);
			assert(alignmentInGap.blocks[b].qPos < alignment.qAlignedSeq.length);
		}
	}
	// Add the blocks for the refined aignment.
	alignment.blocks.insert(alignment.blocks.end(),
													alignmentInGap.blocks.begin(),
													alignmentInGap.blocks.end());
}

bool LengthInBounds(DNALength len, DNALength min, DNALength max) {
	return (len > min and len < max);
}

template<typename T_TargetSequence, typename T_QuerySequence, typename TDBSequence>
void AlignIntervals(T_TargetSequence &genome, T_QuerySequence &read, T_QuerySequence &rcRead,
                    WeightedIntervalSet<ChainedMatchPos> &weightedIntervals,
                    int mutationCostMatrix[][5], 
                    int ins, int del, int sdpTupleSize,
                    int useSeqDB, SequenceIndexDatabase<TDBSequence> &seqDB,
                    vector<T_AlignmentCandidate*> &alignments,
                    MappingParameters &params,
                    int useScoreCutoff, int maxScore,
                    MappingBuffers &mappingBuffers,
                    int procId=0) {
                        
  vector<T_QuerySequence*> forrev;
  forrev.resize(2);
  forrev[Forward] = &read;
  forrev[Reverse] = &rcRead;

  //
  // Use an edit distance scoring function instead of IDS.  Although
  // the IDS should be more accurate, it is more slow, and it is more
  // important at this stage to have faster alignments than accurate,
  // since all alignments are rerun using GuidedAlignment later on.
  //
  DistanceMatrixScoreFunction<DNASequence, FASTQSequence> distScoreFn(SMRTDistanceMatrix, params.insertion, params.deletion);
	distScoreFn.affineOpen = params.affineOpen;
	distScoreFn.affineExtend = params.affineExtend;
  
  //
  // Assume there is at least one interval.
  //
  if (weightedIntervals.size() == 0) 
    return;

  WeightedIntervalSet<ChainedMatchPos>::iterator intvIt = weightedIntervals.begin();
  int alignmentIndex = 0;
  
	if (params.progress) {
		cout << "Mapping intervals " << endl;
	}
  do {

    T_AlignmentCandidate *alignment = alignments[alignmentIndex];
    alignment->clusterWeight= (*intvIt).totalAnchorSize;
    alignment->clusterScore = (*intvIt).pValue;


    //
    // Advance references.  Intervals are stored in reverse order, so
    // go backwards in the list, and alignments are in forward order.
    // That should probably be changed.
    //
    ++alignmentIndex;

    // 
    // Try aligning the read to the genome.
    //
    DNALength matchIntervalStart, matchIntervalEnd;
    matchIntervalStart = (*intvIt).start;
    matchIntervalEnd   = (*intvIt).end;

    bool readOverlapsContigStart    = false;
    bool readOverlapsContigEnd      = false;
    int  startOverlappedContigIndex = 0;
    int  endOverlappedContigIndex   = 0;
    if (params.verbosity > 0) {
      cout << "aligning interval : " << read.length << " " << (*intvIt).start << " " 
           << (*intvIt).end  << " " << (*intvIt).qStart << " " << (*intvIt).qEnd
           << " " << matchIntervalStart << " to " << matchIntervalEnd << " " 
           << params.approximateMaxInsertionRate << " "  << endl;
    }
    assert(matchIntervalEnd >= matchIntervalStart);

    //
    // If using a sequence database, check to make sure that the
    // boundaries of the sequence windows do not overlap with 
    // the boundaries of the reads.  If the beginning is before
    // the boundary, move the beginning up to the start of the read.
    // If the end is past the end boundary of the read, similarly move 
    // the window boundary to the end of the read boundary.

    DNALength tAlignedContigStart = 0;
    int seqDBIndex = 0;


    //
    // Stretch the alignment interval so that it is close to where
    // the read actually starts.
    //
    DNALength subreadStart = read.subreadStart;
    DNALength subreadEnd   = read.subreadEnd;
    if ((*intvIt).GetStrandIndex() == Reverse) {
      subreadEnd   = read.MakeRCCoordinate(read.subreadStart) + 1;
      subreadStart = read.MakeRCCoordinate(read.subreadEnd-1);
    }

		if (params.extendEnds) {
			DNALength lengthBeforeFirstMatch = ((*intvIt).qStart - subreadStart) * params.approximateMaxInsertionRate ;
			DNALength lengthAfterLastMatch   = (subreadEnd - (*intvIt).qEnd) * params.approximateMaxInsertionRate;
			if (matchIntervalStart < lengthBeforeFirstMatch  or params.doGlobalAlignment) {
				matchIntervalStart = 0;
			}
			else {
				matchIntervalStart -= lengthBeforeFirstMatch;
			}
			
			if (genome.length < matchIntervalEnd + lengthAfterLastMatch or params.doGlobalAlignment) {
				matchIntervalEnd = genome.length;
			}
			else {
				matchIntervalEnd += lengthAfterLastMatch;
			}
		}

			

		bool trimMatches = false;
    DNALength intervalContigStartPos, intervalContigEndPos;
    if (useSeqDB) {
      //
      // The sequence db index is the one where the actual match is
      // contained. The matchIntervalStart might be before the sequence
      // index boundary due to the extrapolation of alignment start by
      // insertion rate.  If this is the case, bump up the
      // matchIntervalStart to be at the beginning of the boundary. 
      // Modify bounds similarly for the matchIntervalEnd and the end
      // of a boundary.
      //
      seqDBIndex = seqDB.SearchForIndex((*intvIt).start);
      intervalContigStartPos = seqDB.seqStartPos[seqDBIndex];
      if (intervalContigStartPos > matchIntervalStart) {
        matchIntervalStart = intervalContigStartPos;
      }
      intervalContigEndPos = seqDB.seqStartPos[seqDBIndex+1] - 1;
      if (intervalContigEndPos < matchIntervalEnd) {
        matchIntervalEnd = intervalContigEndPos;
				trimMatches = true;
      }
      alignment->tName    = seqDB.GetSpaceDelimitedName(seqDBIndex);
      alignment->tLength  = intervalContigEndPos - intervalContigStartPos;
      //
      // When there are multiple sequences in the database, store the
      // index of this sequence.  This lets one compare the contigs
      // that reads are mapped to, for instance.
      //
      alignment->tIndex   = seqDBIndex;
    }
    else {
      alignment->tLength     = genome.length;
      alignment->tName       = genome.GetName();
      intervalContigStartPos = 0;
      intervalContigEndPos   = genome.length;
    }

    alignment->qName = read.title;
    int alignScore;
    alignScore = 0;
    alignment->tAlignedSeqPos     = matchIntervalStart;
    alignment->tAlignedSeqLength  = matchIntervalEnd - matchIntervalStart;

    if ((*intvIt).GetStrandIndex() == Forward) {
      alignment->tAlignedSeq.Copy(genome, alignment->tAlignedSeqPos, alignment->tAlignedSeqLength);
      alignment->tStrand = Forward;
    }
    else {
      DNALength rcAlignedSeqPos = genome.MakeRCCoordinate(alignment->tAlignedSeqPos + alignment->tAlignedSeqLength - 1);
      genome.CopyAsRC(alignment->tAlignedSeq, rcAlignedSeqPos, alignment->tAlignedSeqLength);
      // Map forward coordinates into reverse complement.

      intervalContigStartPos    = genome.MakeRCCoordinate(intervalContigStartPos) + 1;
      intervalContigEndPos      = genome.MakeRCCoordinate(intervalContigEndPos - 1);
      swap(intervalContigStartPos, intervalContigEndPos);
      alignment->tAlignedSeqPos = rcAlignedSeqPos;
      alignment->tStrand        = Reverse;
    }
		
		//
    // Reference the subread of the query, always in the forward
    // strand for refining alignments.
		//

		alignment->qAlignedSeqPos = read.subreadStart;
		alignment->qAlignedSeq.ReferenceSubstring(read, read.subreadStart, subreadEnd - subreadStart);
    alignment->qAlignedSeqLength = alignment->qAlignedSeq.length;
    alignment->qLength           = read.length;
    alignment->qStrand           = 0;
		
		//
    // First count how much of the read matches the genome exactly.  
		// It may be too low to bother refining.
		// 
    int intervalSize = 0;
    for (int m = 0; m < intvIt->matches.size(); m++) {
			intervalSize += intvIt->matches[m].l;
		} 

		//		cout << "Interval size: " << intervalSize << endl;
    int subreadLength = forrev[(*intvIt).GetStrandIndex()]->subreadEnd - forrev[(*intvIt).GetStrandIndex()]->subreadStart;
    if ((1.0*intervalSize) / subreadLength < params.sdpBypassThreshold and !params.emulateNucmer) {
      //
      // Not enough of the read maps to the genome, need to use
      // sdp alignment to define the regions of the read that map.
      //
      if (params.refineBetweenAnchorsOnly) {

        //
        // Run SDP alignment only between the genomic anchors,
        // including the genomic anchors as part of the alignment.
        //
        int m;

				vector<ChainedMatchPos> *matches;
        Alignment anchorsOnly;
				//
				// The strand bookkeeping is a bit confusing, so hopefully
				// this will set things straight.
				//
				// If the alignment is forward strand, the coordinates of the
				// blocks are relative to the forward read, starting at 0, not
				// the subread start.
				// If the alignment is reverse strand, the coordinates of the
				// blocks are relative to the reverse strand, starting at the
				// position of the subread on the reverse strand.
				// 
				// The coordinates of the blocks in the genome are always
				// relative to the forward strand on the genome, starting at
				// 0.  
				//

				//
				// The first step to refining between anchors only is to make
				// the anchors relative to the tAlignedSeq.
				//

				matches=(vector<ChainedMatchPos>*)&intvIt->matches;

				//
				// Flip the target.
				//

				DNALength targetOffset = alignment->tAlignedSeqPos;
				if (alignment->tStrand == 1) {
					targetOffset = genome.MakeRCCoordinate(alignment->tAlignedSeqPos+ alignment->tAlignedSeq.length- 1);
				}
				for (m = 0; m < matches->size(); m++) {
					(*matches)[m].t -= targetOffset;
				}

				//
				// If the alignment was on the reverse strand, swap the coordinates so they are on the forward strand of the read.
				//
				if (alignment->tStrand == 1) {
					for (m = 0; m < (*matches).size(); m++) {
						(*matches)[m].q = alignment->qAlignedSeq.length - ((*matches)[m].q + (*matches)[m].l);
						(*matches)[m].t = alignment->tAlignedSeq.length - ((*matches)[m].t + (*matches)[m].l);						
					}
					std::reverse(matches->begin(), matches->end());
				}

				//
				// Remove overlapping matches.
				//
					
				vector<bool> toRemove(matches->size(), false);
				if (matches->size() > 0) {
					m = 0;
					while (m < matches->size()-1) {
						int n = m+1;
						while (n < matches->size() and 
									 (((*matches)[m].q+(*matches)[m].l > (*matches)[n].q) or
										((*matches)[m].t+(*matches)[m].l > (*matches)[n].t))) {
							toRemove[n] = true;
							n++;
						}
						m=n;
					}

					
					m=0;
					int n=0;
					for(n=0;n<matches->size();n++) {
						if (toRemove[n] == false) {
							(*matches)[m] = (*matches)[n];
							m++;
						}
					}

					matches->resize(m);
				}

				  

				if (matches->size() > 0) {
					toRemove.resize(matches->size());
					std::fill(toRemove.begin(), toRemove.end(), false);					
					m = 0;
					int nRemoved =0;
					int blockSize = 0;
					int MAX_GAP=25;
					int MIN_BLOCK=50;
					vector<int> bStart, bEnd, bSize, bScore, bPrev;
					bStart.push_back(0);
					bEnd.push_back(0);
					bSize.push_back(0);
					bScore.push_back(0);
					bPrev.push_back(0);
		      while ( m < matches->size() - 1 ) { 			
						blockSize = 0;
						int blockStart = m;
						int qNextGap = 0;
						int tNextGap = 0;
						int qPrevGap = 0;
						int tPrevGap = 0;
						while (m < matches->size() - 1 and abs(tNextGap - qNextGap) < MAX_GAP) {
							qNextGap = (*matches)[m+1].q - (*matches)[m].q+(*matches)[m].l;
							tNextGap = (*matches)[m+1].t - (*matches)[m].t+(*matches)[m].l;
							blockSize+= (*matches)[m].l;
							m++;
						}
						bStart.push_back(blockStart);
						bEnd.push_back(m-1);
						bSize.push_back(blockSize);
						bScore.push_back(0);
						bPrev.push_back(0);
					}
					bStart.push_back(matches->size()-1);
					bEnd.push_back(matches->size()-1);
					bSize.push_back((*matches)[matches->size()-1].l);
					bScore.push_back(0);
					bPrev.push_back(0);

					int mi, mj;
					float GAP=-0.1;
					int MATCH=1;
					for (mi = 0; mi < bStart.size(); mi++) {
						for (mj = 0; mj < mi; mj++) {
							int p = bEnd[mj];
							int c = bStart[mi];

							int tGap = (*matches)[c].t - (*matches)[p].t;
							int qGap = (*matches)[c].q - (*matches)[p].q;
							int gap = abs(tGap - qGap);
							/*
							if (mj == mi - 1) {
								cerr << mi <<"\tprev\t" << gap << endl;
								}*/
							int score = bScore[mj] + gap*GAP + bSize[mj];
							if (bScore[mi] < score) {
								bScore[mi] = score;
								bPrev[mi] = mj;
							}
						}
					}
					int maxScore = 0;
					int maxIndex = 0;
					
					for (mi = 0; mi < bScore.size(); mi++) {
						if (bScore[mi] > maxScore) {
							maxScore = bScore[mi];
							maxIndex = mi;
						}
					}
					/*
					for (mi = 0; mi < bScore.size(); mi++) {
						cerr << mi << "\t" << bStart[mi] << "\t" << bEnd[mi] << "\t" << bSize[mi] << "\t" << bScore[mi] << "\t" << bPrev[mi] << endl;
					}
					cerr << endl;
					*/
					mi = maxIndex;
					while (mi > 0) {
						bScore[mi] = -1;
						mi= bPrev[mi];
					}
					
					for (mi = 1; mi < bScore.size(); mi++) {
						if (bScore[mi] != -1) {
							for (mj=bStart[mi]; mj <= bEnd[mi]; mj++) {
								toRemove[mj] = true;
							}
						}
					}
					
					/*
				if (seen == false){ 


					ofstream mergedOut("merged.txt");

					for (mi = 0; mi < (*matches).size(); mi++) {
						mergedOut << (*matches)[mi].q << "\t" << (*matches)[mi].t << "\t";
						if (toRemove[mi]) {
							mergedOut << "1";
						}
						else {
							mergedOut << "0";
						}
						mergedOut << endl;
					}
					mergedOut.close();
					seen = true;
				}
					*/

					m=0;
					int n=0;
					for(n=0;n<matches->size();n++) {
						if (toRemove[n] == false) {
							(*matches)[m] = (*matches)[n];
							m++;
						}
					}

					matches->resize(m);
				}

				DNASequence tSubSeq;
				FASTQSequence qSubSeq, mSubSeq;
				
				//
				// Refine alignments between matches.
				//

				MappingBuffers refinementBuffers;
				int f = 0, r = matches->size();
				int tGap, qGap;				

			
				for (m = 0; matches->size() > 0 and m < matches->size()-1; m++) {
          Block block;
          block.qPos = (*matches)[m].q;
          block.tPos = (*matches)[m].t;
          block.length = (*matches)[m].l;
          
          //
          // Find the lengths of the gaps between anchors.
          //

          tGap = (*matches)[m+1].t - ((*matches)[m].t + (*matches)[m].l);
          qGap = (*matches)[m+1].q - ((*matches)[m].q + (*matches)[m].l);

          float gapRatio = (1.0*tGap)/qGap;

					//
					// Add the original block
					//
					alignment->blocks.push_back(block);
					anchorsOnly.blocks.push_back(block);
						
          if (tGap > 0 and qGap > 0) {
						//						cerr << "refining between " << (*matches)[m].q << "-" << (*matches)[m+1].q << "(" << qGap << ") , " << (*matches)[m].t << "-" << (*matches)[m+1].t << "(" << tGap << ")" << endl;
						T_AlignmentCandidate alignmentInGap;

            DNALength tPos, qPos;
            tPos = block.tPos + block.length;
            qPos = block.qPos + block.length;

						//
						// The following is to prevent the attempt to align massive gaps.
						//
						float tRatio = ((float)tGap)/qGap;
						float qRatio = 1/tRatio;
						int maxGap = 100000;
						int adaptSDPTupleSize = sdpTupleSize;
						if (alignment->tAlignedSeq.length  < 1000 && alignment->qAlignedSeq.length < 1000) {
							adaptSDPTupleSize = 7;
						}
						if (tRatio > 0.001 and qRatio > 0.001 and tGap < maxGap and qGap< maxGap) {
							AlignSubstring(alignment->tAlignedSeq, tPos, tGap,
														 alignment->qAlignedSeq, qPos, qGap,
														 params,
														 adaptSDPTupleSize,
														 distScoreFn,
														 refinementBuffers,
														 alignmentInGap);


							//
							// Now, splice the fragment alignment into the current
							// alignment. 
							//
							AppendSubseqAlignment(*alignment, alignmentInGap, qPos, tPos);
						}
					}
        }

				//
				// Add the last block.
				//

				if (matches->size() > 0) {
					Block block;
					block.qPos = (*matches)[matches->size()-1].q;
					block.tPos = (*matches)[matches->size()-1].t;
					if (block.tPos > alignment->tAlignedSeq.length) {
						cout << "ERROR mapping " << read.title << endl;
						read.PrintSeq(cout);
					}
					assert(block.tPos <= alignment->tAlignedSeq.length);
					assert(block.qPos <= alignment->qAlignedSeq.length);
					
					block.length = (*matches)[m].l;
					alignment->blocks.push_back(block);        
					anchorsOnly.blocks.push_back(block);
					}

				//
				// By convention, blocks start at 0, and the
				// alignment->tPos,qPos give the start of the alignment.
				// Modify the block positions so that they are offset by 0.
				//
				if (alignment->blocks.size() > 0) {
					alignment->tPos = alignment->blocks[0].tPos;
					alignment->qPos = alignment->blocks[0].qPos;
					int b;
					int blocksSize = alignment->blocks.size();
					for (b = 0; b < blocksSize ; b++) {
						assert(alignment->tPos <= alignment->blocks[b].tPos);
						assert(alignment->qPos <= alignment->blocks[b].qPos);
						alignment->blocks[b].tPos -= alignment->tPos;
						alignment->blocks[b].qPos -= alignment->qPos;
					}
					for (b = 0; b < anchorsOnly.blocks.size(); b++) {
						anchorsOnly.blocks[b].tPos -= alignment->tPos;
						anchorsOnly.blocks[b].qPos -= alignment->qPos;
					}
					anchorsOnly.tPos = alignment->tPos;
					anchorsOnly.qPos = alignment->qPos;
					
					//
					// Adjacent blocks without gaps (e.g. 50M50M50M) are not yet merged.  Do the merging now.
					//
					int cur = 0, next = 1;
					while (next < alignment->blocks.size()) {
						while (next < alignment->blocks.size() and 
									 alignment->blocks[cur].tPos + alignment->blocks[cur].length == alignment->blocks[next].tPos and 
									 alignment->blocks[cur].qPos + alignment->blocks[cur].length == alignment->blocks[next].qPos) {
						alignment->blocks[cur].length += alignment->blocks[next].length;
						alignment->blocks[next].length = 0;
						next +=1;
						}
						cur = next;
						next += 1;
					}
					cur = 0; next = 0;
					while (next < alignment->blocks.size()) {
						while (next < alignment->blocks.size() and alignment->blocks[next].length == 0) {
							next +=1;
						}
						if (next < alignment->blocks.size()) {
						alignment->blocks[cur] = alignment->blocks[next];
						cur+=1;
						next+=1;
						}
					}
					alignment->blocks.resize(cur);
					
					alignment->gaps.clear();


				
					ComputeAlignmentStats(*alignment, alignment->qAlignedSeq.seq, alignment->tAlignedSeq.seq,
																distScoreFn, params.affineAlign);
				}
			}
      else {
        alignScore = SDPAlign(alignment->qAlignedSeq, alignment->tAlignedSeq, distScoreFn, 
                              sdpTupleSize, params.sdpIns, params.sdpDel, params.indelRate*3, 
                              *alignment, mappingBuffers, 
                              Local, 
                              params.detailedSDPAlignment, 
															params.extendFrontAlignment,
															params.sdpPrefix, params.recurse, params.recurseOver, params.sdpMaxAnchorsPerPosition);
        ComputeAlignmentStats(*alignment, alignment->qAlignedSeq.seq, alignment->tAlignedSeq.seq,
                              distScoreFn, params.affineAlign);
      }
    }

    else {
      //
      // The anchors used to anchor the sequence are sufficient to extend the alignment.
      //
      int m;
      for (m = 0; m < (*intvIt).matches.size(); m++ ){
        Block block;
        block.qPos = (*intvIt).matches[m].q - alignment->qAlignedSeqPos;
        block.tPos = (*intvIt).matches[m].t - alignment->tAlignedSeqPos;
        block.length = (*intvIt).matches[m].l;
        alignment->blocks.push_back(block);
      }
    }



    //
    //  The anchors/sdp alignments may leave portions of the read
    //  unaligned at the beginning and end.  If the parameters
    //  specify extending alignments, try and align extra bases at
    //  the beginning and end of alignments.
    if (params.extendAlignments) {

      //
      // Modify the alignment so that the start and end of the
      // alignment strings are at the alignment boundaries.
      //
      // Since the query sequence is pointing at a subsequence of the
      // read (and is always in the forward direction), just reference
      // a new portion of the read.
      alignment->qAlignedSeqPos = alignment->qAlignedSeqPos + alignment->qPos;
      alignment->qAlignedSeqLength = alignment->QEnd();
      alignment->qAlignedSeq.ReferenceSubstring(read, alignment->qAlignedSeqPos, alignment->qAlignedSeqLength );
      alignment->qPos = 0; 

      //
      // Since the target sequence may be on the forward or reverse
      // strand, a copy of the subsequence is made, and the original
      // sequence free'd.
      //
      DNASequence tSubseq;
      alignment->tAlignedSeqPos = alignment->tAlignedSeqPos + alignment->tPos;
      alignment->tAlignedSeqLength = alignment->TEnd();
      tSubseq.Copy(alignment->tAlignedSeq, alignment->tPos, alignment->tAlignedSeqLength);      
      alignment->tPos = 0;

      alignment->tAlignedSeq.Free();
      alignment->tAlignedSeq.TakeOwnership(tSubseq);

          
      DNALength maximumExtendLength = 500;

      if (alignment->blocks.size() > 0 ){
        int lastAlignedBlock = alignment->blocks.size() - 1;
        DNALength lastAlignedQPos  = alignment->blocks[lastAlignedBlock].QEnd() + alignment->qPos + alignment->qAlignedSeqPos;
        DNALength lastAlignedTPos  = alignment->blocks[lastAlignedBlock].TEnd() + alignment->tPos + alignment->tAlignedSeqPos;
        T_AlignmentCandidate extendedAlignmentForward, extendedAlignmentReverse;
        int forwardScore, reverseScore;

        SMRTSequence  readSuffix;
        DNALength     readSuffixLength;
        DNASequence   genomeSuffix;
        DNALength     genomeSuffixLength;
        
        SMRTSequence   readPrefix;
        DNALength     readPrefixLength;
        DNASequence   genomePrefix;
        DNALength     genomePrefixLength;

        //
        // Align the entire end of the read if it is short enough.
        //
        readSuffixLength = min(read.length - lastAlignedQPos, maximumExtendLength);
        if (readSuffixLength > 0) {
          readSuffix.ReferenceSubstring(read, lastAlignedQPos, readSuffixLength);
        }
        else {
          readSuffix.length = 0;
        }
        
        //
        // Align The entire end of the genome up to the maximum extend length;
        //
        genomeSuffixLength = min(intervalContigEndPos - lastAlignedTPos, maximumExtendLength);
        if (genomeSuffixLength > 0) {
          if (alignment->tStrand == Forward) {
            genomeSuffix.Copy(genome, lastAlignedTPos, genomeSuffixLength);
          }
          else {
            ((DNASequence)genome).CopyAsRC(genomeSuffix, lastAlignedTPos, genomeSuffixLength);
          }
        }
        else {
          genomeSuffix.length = 0;
        }
        forwardScore = 0;
        if (readSuffix.length > 0 and genomeSuffix.length > 0) {
          forwardScore = ExtendAlignmentForward(readSuffix, 0,
                                                genomeSuffix, 0,
                                                params.extendBandSize, 
                                                // Reuse buffers to speed up alignment
                                                mappingBuffers.scoreMat,
                                                mappingBuffers.pathMat,
                                                // Do the alignment in the forward direction.
                                                extendedAlignmentForward,
                                                distScoreFn,
                                                1, // don't bother attempting
                                                // to extend the alignment
                                                // if one of the sequences
                                                // is less than 1 base long
                                                params.maxExtendDropoff);
        }
        
        if ( forwardScore < 0 ) {
          //
          // The extended alignment considers the whole genome, but
          // should be modified to be starting at the end of where 
          // the original alignment left off.
          //
          if (params.verbosity > 0) {
            cout << "forward extended an alignment of score " << alignment->score << " with score " << forwardScore << " by " << extendedAlignmentForward.blocks.size() << " blocks and length " << extendedAlignmentForward.blocks[extendedAlignmentForward.blocks.size()-1].qPos << endl;
          }
          extendedAlignmentForward.tAlignedSeqPos = lastAlignedTPos;

          extendedAlignmentForward.qAlignedSeqPos = lastAlignedQPos;

          genomeSuffix.length = extendedAlignmentForward.tPos + extendedAlignmentForward.TEnd();
          alignment->tAlignedSeq.Append(genomeSuffix);
          alignment->qAlignedSeq.length += extendedAlignmentForward.qPos + extendedAlignmentForward.QEnd();
          assert(alignment->qAlignedSeq.length <= read.length);
          alignment->AppendAlignment(extendedAlignmentForward);
        }
        genomeSuffix.Free();

        DNALength firstAlignedQPos = alignment->qPos + alignment->qAlignedSeqPos;
        DNALength firstAlignedTPos = alignment->tPos + alignment->tAlignedSeqPos;
        
        
        readPrefixLength = min(firstAlignedQPos, maximumExtendLength);
        if (readPrefixLength > 0) {
          readPrefix.ReferenceSubstring(read, firstAlignedQPos-readPrefixLength, readPrefixLength);
        }
        else {
          readPrefix.length = 0;
        }
        
        genomePrefixLength = min(firstAlignedTPos - intervalContigStartPos, maximumExtendLength);
        if (genomePrefixLength > 0) {
          if (alignment->tStrand == 0) {
            genomePrefix.Copy(genome, firstAlignedTPos - genomePrefixLength, genomePrefixLength);
          }
          else {
            ((DNASequence)genome).MakeRC(genomePrefix, firstAlignedTPos - genomePrefixLength, genomePrefixLength);
          }
        }
        reverseScore = 0;
        if (readPrefix.length > 0 and genomePrefix.length > 0) {
          reverseScore = ExtendAlignmentReverse(readPrefix, readPrefix.length-1,
                                                genomePrefix, genomePrefixLength - 1,
                                                params.extendBandSize, //k
                                                mappingBuffers.scoreMat,
                                                mappingBuffers.pathMat,
                                                extendedAlignmentReverse,
                                                distScoreFn,
                                                1, // don't bother attempting
                                                // to extend the alignment
                                                // if one of the sequences
                                                // is less than 1 base long
                                                params.maxExtendDropoff);
        }
      
        if (reverseScore < 0 ) {
          //
          // Make alignment->tPos relative to the beginning of the
          // extended alignment so that when it is appended, the
          // coordinates match correctly.
          if (params.verbosity > 0) {
            cout << "reverse extended an alignment of score " << alignment->score << " with score " << reverseScore << " by " << extendedAlignmentReverse.blocks.size() << " blocks and length " << extendedAlignmentReverse.blocks[extendedAlignmentReverse.blocks.size()-1].qPos << endl;
          }
          extendedAlignmentReverse.tAlignedSeqPos = firstAlignedTPos - genomePrefixLength;
          extendedAlignmentReverse.qAlignedSeqPos = firstAlignedQPos - readPrefixLength;
          extendedAlignmentReverse.AppendAlignment(*alignment);

          genomePrefix.Append(alignment->tAlignedSeq, genomePrefix.length - alignment->tPos);
          alignment->tAlignedSeq.Free();
          alignment->tAlignedSeq.TakeOwnership(genomePrefix);
          
          alignment->blocks = extendedAlignmentReverse.blocks;
          
          alignment->tAlignedSeqPos = extendedAlignmentReverse.tAlignedSeqPos;
          alignment->tPos = extendedAlignmentReverse.tPos;


          alignment->qAlignedSeqPos     = extendedAlignmentReverse.qAlignedSeqPos; 
          alignment->qAlignedSeq.length = readPrefix.length + alignment->qAlignedSeq.length;
          alignment->qPos               = extendedAlignmentReverse.qPos;
          alignment->qAlignedSeq.seq    = readPrefix.seq;
          //
          // Make sure the two ways of accounting for aligned sequence
          // length are in sync.  This needs to go.
          //
          if (alignment->blocks.size() > 0) {
            int lastBlock = alignment->blocks.size() - 1;
            alignment->qAlignedSeqLength = alignment->qAlignedSeq.length;
            alignment->tAlignedSeqLength = alignment->tAlignedSeq.length;
          }                
          else {
            alignment->qAlignedSeqLength = alignment->qAlignedSeq.length = 0;
            alignment->tAlignedSeqLength = alignment->tAlignedSeq.length = 0;
          }
        }
        else {
          genomePrefix.Free();
        }
      }
    }

    ComputeAlignmentStats(*alignment, 
                          alignment->qAlignedSeq.seq,
                          alignment->tAlignedSeq.seq, distScoreFn, params.affineAlign);


    intvIt++;
		if (params.progress) {
			cout << "done with interval " << alignmentIndex << endl;
		}
		if (params.verbosity > 0) {
      cout << "interval align score: " << alignScore << endl;
      StickPrintAlignment(*alignment,
                          (DNASequence&) alignment->qAlignedSeq,
                          (DNASequence&) alignment->tAlignedSeq,
                          cout,
                          0, alignment->tAlignedSeqPos);

    }

  } while (intvIt != weightedIntervals.end());
}


bool FirstContainsSecond(DNALength aStart, DNALength aEnd, DNALength bStart, DNALength bEnd) {
  return ((bStart > aStart and bEnd <= aEnd) or
          (bStart >= aStart and bEnd < aEnd));
}

template<typename T_Sequence>
bool CheckForSufficientMatch(T_Sequence &read, vector<T_AlignmentCandidate*> &alignmentPtrs, MappingParameters &params) {
  if (alignmentPtrs.size() > 0 and alignmentPtrs[0]->score < params.maxScore) {
    return true;
  }
  else {
    return false;
  }
}

void DeleteAlignments(vector<T_AlignmentCandidate*> &alignmentPtrs, int startIndex=0) {
  int i;
  for (i = startIndex; i < alignmentPtrs.size(); i++ ) {
    delete alignmentPtrs[i];
  }
  alignmentPtrs.resize(0);
}




int RemoveLowQualitySDPAlignments(int readLength, vector<T_AlignmentCandidate*> &alignmentPtrs, MappingParameters &params) {
  // Just a hack.  For now, assume there is at least 1 match per 50 bases.
  int totalBasesMatched = 0;
  int a;
  for (a = 0; a < alignmentPtrs.size(); a++) {
    int b;
    for (b = 0; b < alignmentPtrs[a]->blocks.size(); b++) {
      totalBasesMatched += alignmentPtrs[a]->blocks[b].length;
    }
    int expectedMatches = params.sdpTupleSize/50.0 * readLength;
    if (totalBasesMatched < expectedMatches) {
      delete alignmentPtrs[a];
      alignmentPtrs[a] = NULL;
    }
  }
  int packedAlignmentIndex = 0;
  for (a = 0; a < alignmentPtrs.size(); a++) {
    if (alignmentPtrs[a] != NULL) {
      alignmentPtrs[packedAlignmentIndex] = alignmentPtrs[a];
      packedAlignmentIndex++;
    }
  }
  alignmentPtrs.resize(packedAlignmentIndex);
  return packedAlignmentIndex;
}


template<typename T_Sequence>
int RemoveLowQualityAlignments(T_Sequence &read, vector<T_AlignmentCandidate*> &alignmentPtrs, MappingParameters &params) {
  if (params.verbosity > 0) {
    cout << "checking at least " << alignmentPtrs.size() << " alignments to see if they are accurate." << endl;
  }
  UInt i;
  for (i = 0; i < MIN(params.nCandidates, alignmentPtrs.size()); i++) { 
    if (params.verbosity > 0) {
      cout << "Quality check  " << i << " " << alignmentPtrs[i]->score << endl;
    }
	if (params.maxGap > 0) {
	  UInt j;
	  for (j = 1; j < alignmentPtrs[i]->blocks.size(); j++) {
		int gq, gt;
		gq = alignmentPtrs[i]->blocks[j].qPos - alignmentPtrs[i]->blocks[j-1].qPos - alignmentPtrs[i]->blocks[j].length;
		gt = alignmentPtrs[i]->blocks[j].tPos - alignmentPtrs[i]->blocks[j-1].tPos - alignmentPtrs[i]->blocks[j].length;
		if ((gq >= params.maxGap and gt < params.maxGap) or 
			(gt > params.maxGap and gq < params.maxGap) ) {
		  // 
		  // Flag alignment to be deleted.
		  //
		  alignmentPtrs[i]->score = params.maxScore;
		  break;
		}
	  }
	}
    if (alignmentPtrs[i]->blocks.size() == 0 or
        alignmentPtrs[i]->score > params.maxScore) {
      //
      // Since the alignments are sorted according to alignment
      // score, once one of the alignments is too low of a score,
      // all remaining alignments are also too low, and should be
      // removed as well.  Do that all at once.
      //
      if (alignmentPtrs[i]->blocks.size() == 0 and params.verbosity > 0) {
        cout << "Removing empty alignment " << alignmentPtrs[i]->qName << endl;
      }
      if (params.verbosity  > 0) { 
        cout << alignmentPtrs[i]->qName << " alignment " << i << " is too low of a score." << alignmentPtrs[i]->score << endl;
      }
	  
      int deletedIndex = i;
      for (; deletedIndex < alignmentPtrs.size(); deletedIndex++) {
        delete alignmentPtrs[deletedIndex];
        alignmentPtrs[deletedIndex] = NULL;
      }
      alignmentPtrs.erase(i + alignmentPtrs.begin(), alignmentPtrs.end());
      break;
    }
    else {
      if (params.verbosity > 0) {
        cout << "Keeping alignment " << i << " " << alignmentPtrs[i]->qPos << " " << alignmentPtrs[i]->qPos + alignmentPtrs[i]->qLength
             << " " << alignmentPtrs[i]->tName << " " << alignmentPtrs[i]->tPos << " " << alignmentPtrs[i]->tPos + alignmentPtrs[i]->tLength << " " << alignmentPtrs[i]->tStrand
             << " from score: " << alignmentPtrs[i]->score << endl;
      }
    }
  }
  return alignmentPtrs.size();
}

int RemoveOverlappingAlignments(vector<T_AlignmentCandidate*> &alignmentPtrs, MappingParameters &params) {
  vector<unsigned char> alignmentIsContained;
  alignmentIsContained.resize(alignmentPtrs.size());
  std::fill(alignmentIsContained.begin(), alignmentIsContained.end(), false);

  int j;
  int numContained = 0;
  int curNotContained = 0;
    
  if (alignmentPtrs.size() > 0) {
    UInt i;
    for (i = 0; i < alignmentPtrs.size()-1; i++ ){
      T_AlignmentCandidate *aref = alignmentPtrs[i];
      if (aref->pctSimilarity < params.minPctIdentity) {
        continue;
      }
      for (j = i + 1; j < alignmentPtrs.size(); j++ ){
        //
        // Make sure this alignment isn't already removed.
        //
        if (alignmentIsContained[j]) {
          continue;
        }
            
        //
        // Only check for containment if the two sequences are from the same contig.
        //
        if (alignmentPtrs[i]->tIndex != alignmentPtrs[j]->tIndex) {
          continue;
        }

        // 
        // Check for an alignment that is fully overlapping another 
        // alignment.
        if (aref->GenomicTBegin() <= alignmentPtrs[j]->GenomicTBegin() and
            aref->GenomicTEnd() >= alignmentPtrs[j]->GenomicTEnd() and 
            alignmentPtrs[i]->tStrand == alignmentPtrs[j]->tStrand) {
          //
          // Alignment i is contained in j is only true if it has a worse score.
          //
          if (aref->score <= alignmentPtrs[j]->score) {
            alignmentIsContained[j] = true;
          }
          if (params.verbosity >= 2) {
            cout << "alignment " << i << " is contained in " << j << endl;
            cout << aref->tAlignedSeqPos << " " <<  alignmentPtrs[j]->tAlignedSeqPos << " "
                 << aref->tAlignedSeqPos + aref->tAlignedSeqLength << " " 
                 << alignmentPtrs[j]->tAlignedSeqPos + alignmentPtrs[j]->tAlignedSeqLength << endl;
          }
        }
        else if (alignmentPtrs[j]->GenomicTBegin() <= aref->GenomicTBegin() and
                 alignmentPtrs[j]->GenomicTEnd()   >= aref->GenomicTEnd() and 
                 alignmentPtrs[i]->tStrand == alignmentPtrs[j]->tStrand) {
          if (params.verbosity >= 2) {
            cout << "ALIGNMENT " << j << " is contained in " << i << endl;
            cout << alignmentPtrs[j]->tAlignedSeqPos << " " <<  aref->tAlignedSeqPos << " "
                 << alignmentPtrs[j]->tAlignedSeqPos + alignmentPtrs[j]->tAlignedSeqLength << " " 
                 <<  aref->tAlignedSeqPos + aref->tAlignedSeqLength << endl;
          }
          if (alignmentPtrs[j]->score <= aref->score) {
            alignmentIsContained[i] = true;
          }
        }
      }
    }
    for (i = 0; i < alignmentPtrs.size(); i++) {
      T_AlignmentCandidate *aref = alignmentPtrs[i];
      if (alignmentIsContained[i]) {
        delete alignmentPtrs[i];
        alignmentPtrs[i] = NULL;
        numContained++;
      }
      else {
        alignmentPtrs[curNotContained] = aref;
        ++curNotContained;
      }
    }
    alignmentPtrs.resize(alignmentPtrs.size() - numContained);
  } 
  return alignmentPtrs.size();
}

template<typename T_RefSequence, typename T_Sequence>
void RefineAlignments(vector<T_Sequence*> &bothQueryStrands,
                      T_RefSequence &genome,
                      vector<T_AlignmentCandidate*> &alignmentPtrs, MappingParameters &params, MappingBuffers &mappingBuffers) {

  
  UInt i;
  for (i = 0; i < alignmentPtrs.size(); i++ ) {
    RefineAlignment(*bothQueryStrands[0], genome, *alignmentPtrs[i], params, mappingBuffers);
  }
  //
  // It's possible the alignment references change their order after running
  // the local alignments.  This is made into a parameter rather than resorting
  // every time so that the performance gain by resorting may be measured.
  //
  if (params.sortRefinedAlignments) {
    std::sort(alignmentPtrs.begin(), alignmentPtrs.end(), SortAlignmentPointersByScore());
  }
}
    



void AssignRefContigLocation(T_AlignmentCandidate &alignment, SequenceIndexDatabase<FASTQSequence> &seqdb, DNASequence &genome) {
    //
    // If the sequence database is used, the start position of
    // the alignment is relative to the start of the chromosome,
    // not the entire index.  Subtract off the start position of
    // the chromosome to get the true position.
    //
  DNALength forwardTPos;
  int seqDBIndex;
  if (alignment.tStrand == 0) {
    forwardTPos = alignment.tAlignedSeqPos;
    seqDBIndex = seqdb.SearchForIndex(forwardTPos);
    alignment.tAlignedSeqPos -= seqdb.seqStartPos[seqDBIndex];
  }
  else {
    //
    // Flip coordinates into forward strand in order to find the boundaries 
    // of the contig, then reverse them in order to find offset.
    //

    // Find the reverse complement coordinate of the index of the last aligned base.
    assert(alignment.tAlignedSeqLength > 0);
    forwardTPos = genome.MakeRCCoordinate(alignment.tAlignedSeqPos + alignment.tAlignedSeqLength - 1);
    seqDBIndex  = seqdb.SearchForIndex(forwardTPos);

    
    //
    // Find the reverse comlement coordinate of the last base of this
    // sequence.  This would normally be the start of the next contig
    // -1 to get the length, but since an 'N' is added between every
    // pair of sequences, this is -2.
    //
    DNALength reverseTOffset;
    reverseTOffset = genome.MakeRCCoordinate(seqdb.seqStartPos[seqDBIndex+1]-2);
    alignment.tAlignedSeqPos -= reverseTOffset;
  }
}

void AssignRefContigLocations(vector<T_AlignmentCandidate*> &alignmentPtrs, SequenceIndexDatabase<FASTQSequence> &seqdb, DNASequence &genome) {
  
  UInt i;
  for (i = 0; i < alignmentPtrs.size(); i++) {
    T_AlignmentCandidate *aref = alignmentPtrs[i];
    AssignRefContigLocation(*aref, seqdb, genome);
  }
}


template<typename T_RefSequence>
void AssignGenericRefContigName(vector<T_AlignmentCandidate*> &alignmentPtrs, T_RefSequence &genome) {
  UInt i;
  for (i = 0; i < alignmentPtrs.size(); i++) {
    T_AlignmentCandidate *aref = alignmentPtrs[i];
    aref->tName = genome.title;
  }
}

int basesAligned = 0;
int nReads = 0;

template<typename T_Sequence, typename T_RefSequence, typename T_SuffixArray, typename T_TupleCountTable>
void MapRead(T_Sequence &read, T_Sequence &readRC, T_RefSequence &genome, 
             T_SuffixArray &sarray, 
             BWT &bwt,
             SeqBoundaryFtr<FASTQSequence> &seqBoundary, 
             T_TupleCountTable &ct,
             SequenceIndexDatabase<FASTQSequence> &seqdb,
             MappingParameters &params,
             MappingMetrics    &metrics,
             vector<T_AlignmentCandidate*> &alignmentPtrs, 
             MappingBuffers &mappingBuffers,
             MappingIPC *mapData) {


  bool matchFound;
  WeightedIntervalSet<ChainedMatchPos> topIntervals(params.nCandidates);
  int numKeysMatched=0, rcNumKeysMatched=0;
  int expand = params.minExpand;
  metrics.clocks.total.Tick();
  int nTotalCells = 0;
  int forwardNumBasesMatched = 0, reverseNumBasesMatched = 0;
  do {
    matchFound = false;
    mappingBuffers.matchPosList.clear();
    mappingBuffers.rcMatchPosList.clear();
    alignmentPtrs.clear();
    topIntervals.clear();
    params.anchorParameters.expand = expand;

    metrics.clocks.mapToGenome.Tick();
		if (params.verbosity > 1) {
			cerr << "mapping to genome." << endl;
		}
    if (params.useSuffixArray) {
      params.anchorParameters.lcpBoundsOutPtr = mapData->lcpBoundsOutPtr;
      numKeysMatched   = 
        MapReadToGenome(genome, sarray, read, 
						params.lookupTableLength, 
						mappingBuffers.matchPosList,
                        params.anchorParameters);
      
      //
      // Only print values for the read in forward direction (and only
      // the first read). 
      //
      mapData->lcpBoundsOutPtr = NULL;
      if (!params.forwardOnly) {
        rcNumKeysMatched = 
          MapReadToGenome(genome, sarray, readRC, params.lookupTableLength, mappingBuffers.rcMatchPosList, 
                          params.anchorParameters);
      }
    }
    else if (params.useBwt){ 
      numKeysMatched   = MapReadToGenome(bwt, read, read.subreadStart, read.subreadEnd, 
                                         mappingBuffers.matchPosList, params.anchorParameters, forwardNumBasesMatched);
      if (!params.forwardOnly) {
        rcNumKeysMatched = MapReadToGenome(bwt, readRC, readRC.subreadStart, readRC.subreadEnd, 
                                           mappingBuffers.rcMatchPosList, params.anchorParameters, reverseNumBasesMatched); 
      }
    }

    //
    // Look to see if only the anchors are printed.
    if (params.anchorFileName != "") {
      int i;
      if (params.nProc > 1) {
#ifdef __APPLE__
        sem_wait(semaphores.writer);
#else
        sem_wait(&semaphores.writer);
#endif
      }
			cout << "writing anchors to " << params.anchorFileName << endl;
      for (i = 0; i < mappingBuffers.matchPosList.size(); i++) {
        *mapData->anchorFilePtr << mappingBuffers.matchPosList[i].q << " " << mappingBuffers.matchPosList[i].t << " " << mappingBuffers.matchPosList[i].l << " " << 0 << endl;
      }
      for (i = 0; i < mappingBuffers.rcMatchPosList.size(); i++) {
        *mapData->anchorFilePtr << read.length - mappingBuffers.rcMatchPosList[i].q << " " << mappingBuffers.rcMatchPosList[i].t << " " << mappingBuffers.rcMatchPosList[i].l << " " << 1 << endl;
      }
      
      if (params.nProc > 1) {
#ifdef __APPLE__
        sem_post(semaphores.writer);
#else
        sem_post(&semaphores.writer);
#endif
      }
    }

    metrics.totalAnchors += mappingBuffers.matchPosList.size() + mappingBuffers.rcMatchPosList.size();
    metrics.clocks.mapToGenome.Tock();

    metrics.clocks.sortMatchPosList.Tick();
    SortMatchPosList(mappingBuffers.matchPosList);
    SortMatchPosList(mappingBuffers.rcMatchPosList);
    metrics.clocks.sortMatchPosList.Tock();
    PValueWeightor lisPValue(read, genome, ct.tm, &ct);
		//    MultiplicityPValueWeightor lisPValueByWeight(genome);

    LISSumOfLogPWeightor<T_GenomeSequence,vector<ChainedMatchPos> > lisPValueByLogSum(genome);

    LISSizeWeightor<vector<ChainedMatchPos> > lisWeightFn;

    IntervalSearchParameters intervalSearchParameters;
		params.InitializeIntervalSearchParameters(intervalSearchParameters);
    //
    // If specified, only align a band from the anchors.
    //
    DNALength squareRefLength = read.length * 1.25 + params.limsAlign;    
    if (params.limsAlign != 0) {
      int fi;
      for (fi = 0; fi < mappingBuffers.matchPosList.size(); fi++) {
        if (mappingBuffers.matchPosList[fi].t >= squareRefLength) { break; }
      }
      if (fi < mappingBuffers.matchPosList.size()) {
        mappingBuffers.matchPosList.resize(fi);
      }
    }

	
    metrics.clocks.findMaxIncreasingInterval.Tick();
    
    //
    // For now say that something that has a 50% chance of happening
    // by chance is too high of a p value. This is probably many times
    // the size.
    //
    intervalSearchParameters.maxPValue = log(0.5); 


    //
    // Remove anchors that are fully encompassed by longer ones.  This
    // speeds up limstemplate a lot.
    //

    RemoveOverlappingAnchors(mappingBuffers.matchPosList);
    RemoveOverlappingAnchors(mappingBuffers.rcMatchPosList);

	if (params.progress != 0) {
	  cout << "anchors " << mappingBuffers.matchPosList.size() << " " << mappingBuffers.rcMatchPosList.size() << endl;
	}

    if (params.pValueType == 0) {
      int original = mappingBuffers.matchPosList.size();

      int numMerged = 0;
      if (params.printDotPlots) {
        ofstream dotPlotOut;
        string dotPlotName = string(read.title) + ".anchors";
        CrucialOpen(dotPlotName, dotPlotOut, std::ios::out);
        int mp;
        for (mp = 0; mp < mappingBuffers.matchPosList.size(); mp++ ){
          dotPlotOut << mappingBuffers.matchPosList[mp].q << " " << mappingBuffers.matchPosList[mp].t << " " << mappingBuffers.matchPosList[mp].l << " " << endl;
        }
        dotPlotOut.close();
      }

	  DNALength maxGapLength = 50000;
	  DNALength intervalLength = (read.subreadEnd - read.subreadStart) * (1 + params.indelRate);
	  if (intervalLength  - (read.subreadEnd - read.subreadStart) > maxGapLength) {
		intervalLength = (read.subreadEnd - read.subreadStart) + maxGapLength;
	  }
	  if ( params.piecewiseMatch) {
			LISPValueWeightor<T_GenomeSequence, DNATuple, vector<DirMatch > > pwPValue(read, genome, ct.tm, &ct);
			LISSizeWeightor<vector<DirMatch> > pwWeightFn;			
			PiecewiseMatch(mappingBuffers.matchPosList,
										 mappingBuffers.rcMatchPosList,
										 params.nCandidates,
										 seqBoundary,
										 pwPValue,
										 pwWeightFn,
										 topIntervals, genome, read, intervalSearchParameters);


	  }
	  else {
			FindMaxIncreasingInterval(Forward,
															mappingBuffers.matchPosList,
															// allow for indels to stretch out the mapping of the read.
															intervalLength, params.nCandidates,
															seqBoundary,
															lisPValue,//lisPValue2,
															lisWeightFn,
															topIntervals, genome, read, intervalSearchParameters,
															&mappingBuffers.globalChainEndpointBuffer, 
															&mappingBuffers.sdpFragmentSet, read.title);
		// Uncomment when the version of the weight functor needs the sequence.
		
		FindMaxIncreasingInterval(Reverse, mappingBuffers.rcMatchPosList,
															intervalLength, params.nCandidates, 
															seqBoundary,
															lisPValue,//lisPValue2
															lisWeightFn,
															topIntervals, genome, readRC, intervalSearchParameters,
															&mappingBuffers.globalChainEndpointBuffer,
															&mappingBuffers.sdpFragmentSet, read.title);
	  }
	}
			/*    else if (params.pValueType == 1) {
      FindMaxIncreasingInterval(Forward,
                                mappingBuffers.matchPosList,
                                // allow for indels to stretch out the mapping of the read.
                                (DNALength) ((read.subreadEnd - read.subreadStart) * (1 + params.indelRate)), params.nCandidates,
                                seqBoundary,
                                lisPValueByWeight,
                                lisWeightFn,
                                topIntervals, genome, read, intervalSearchParameters,
                                &mappingBuffers.globalChainEndpointBuffer,
                                read.title);


      FindMaxIncreasingInterval(Reverse, mappingBuffers.rcMatchPosList,
                                (DNALength) ((read.subreadEnd - read.subreadStart) * (1 + params.indelRate)), params.nCandidates, 
                                seqBoundary,
                                lisPValueByWeight,
                                lisWeightFn,
                                topIntervals, genome, readRC, intervalSearchParameters,
                                &mappingBuffers.globalChainEndpointBuffer,
                                read.title);
																}*/
    else if (params.pValueType == 2) {
      FindMaxIncreasingInterval(Forward,
                                mappingBuffers.matchPosList,
                                // allow for indels to stretch out the mapping of the read.
                                (DNALength) ((read.subreadEnd - read.subreadStart) * (1 + params.indelRate)), params.nCandidates,
                                seqBoundary,
                                lisPValueByLogSum,
                                lisWeightFn,
                                topIntervals, genome, read, intervalSearchParameters,
                                &mappingBuffers.globalChainEndpointBuffer,
                                &mappingBuffers.sdpFragmentSet, 
                                read.title);

      FindMaxIncreasingInterval(Reverse, mappingBuffers.rcMatchPosList,
                                (DNALength) ((read.subreadEnd - read.subreadStart) * (1 + params.indelRate)), params.nCandidates, 
                                seqBoundary,
                                lisPValueByLogSum,
                                lisWeightFn,
                                topIntervals, genome, readRC, intervalSearchParameters,
                                &mappingBuffers.globalChainEndpointBuffer,
                                &mappingBuffers.sdpFragmentSet, 
                                read.title);
    }

    metrics.clocks.findMaxIncreasingInterval.Tock();    

    //
    // Print verbose output.
    //

    WeightedIntervalSet<ChainedMatchPos>::iterator topIntIt, topIntEnd;
    topIntEnd = topIntervals.end();
		
	if (params.removeContainedIntervals) {
	  topIntervals.RemoveContained();
	}
	
    //
    // Allocate candidate alignments on the stack.  Each interval is aligned.
    //
    alignmentPtrs.resize(topIntervals.size());

    if (params.verbosity > 0) {
      int topintind = 0;
      cout << " intv: index start end qstart qend seq_boundary_start seq_boundary_end strand pvalue " << endl;
      for (topIntIt = topIntervals.begin();topIntIt != topIntEnd ; ++topIntIt) {
        cout << " intv: " << topintind << " " << (*topIntIt).start << " " 
             << (*topIntIt).end << " " 
             << (*topIntIt).qStart << " " << (*topIntIt).qEnd << " "
             << seqBoundary((*topIntIt).start) << " " << seqBoundary((*topIntIt).end) << " "
						 << (*topIntIt).GetStrandIndex() << " " 
             << (*topIntIt).pValue << endl;
        if (params.verbosity > 2) {
          int m;
          for (m = 0; m < (*topIntIt).matches.size(); m++) {
            cout << " (" << (*topIntIt).matches[m].q << ", " << (*topIntIt).matches[m].t << ", " << (*topIntIt).matches[m].l << ") ";
          }
          cout << endl;
        }
        ++topintind;
      }
    }

    UInt i;
    for (i = 0; i < alignmentPtrs.size(); i++ ) {
      alignmentPtrs[i] = new T_AlignmentCandidate;
    }
    metrics.clocks.alignIntervals.Tick();
    AlignIntervals( genome, read, readRC,
                    topIntervals,
                    SMRTDistanceMatrix,
                    params.indel, params.indel, 
                    params.sdpTupleSize, 
                    params.useSeqDB, seqdb,
                    alignmentPtrs,
                    params,
                    params.useScoreCutoff, params.maxScore,
                    mappingBuffers,
                    params.startRead );


    std::sort(alignmentPtrs.begin(), alignmentPtrs.end(), SortAlignmentPointersByScore());
    metrics.clocks.alignIntervals.Tock();

    //
    // Evalutate the matches that are found for 'good enough'.
    //
      
    matchFound = CheckForSufficientMatch(read, alignmentPtrs, params);

    //
    // When no proper alignments are found, the loop will resume.
    // Delete all alignments because they are bad.
    // 
    if (expand < params.maxExpand and matchFound == false) {
      DeleteAlignments(alignmentPtrs, 0);
    }


    //
    // Record some metrics that show how long this took to run per base.
    //

    if (alignmentPtrs.size() > 0) {
      metrics.RecordNumAlignedBases(read.length);
      metrics.RecordNumCells(alignmentPtrs[0]->nCells);
    }

    if (matchFound == true) {
      metrics.totalAnchorsForMappedReads += mappingBuffers.matchPosList.size() + mappingBuffers.rcMatchPosList.size();
    }
    ++expand;
  } while ( expand <= params.maxExpand and matchFound == false);
  metrics.clocks.total.Tock();
  UInt i;
  int totalCells = 0;
  for (i = 0; i< alignmentPtrs.size(); i++) {
    totalCells += alignmentPtrs[i]->nCells;
  }
  //
  //  Some of the alignments are to spurious regions. Delete the
  //  references that have too small of a score.
  //

  int effectiveReadLength = 0;
  for (i = 0; i< read.length; i++) {
    if (read.seq[i] != 'N') effectiveReadLength++;
  }
  if (params.sdpFilterType == 0) {
    RemoveLowQualityAlignments(read, alignmentPtrs, params);
  }
  else if (params.sdpFilterType == 1) {
    RemoveLowQualitySDPAlignments(effectiveReadLength, alignmentPtrs, params);
  }

  //
  // Now remove overlapping alignments.
  //

  vector<T_Sequence*> bothQueryStrands;
  bothQueryStrands.resize(2);
  bothQueryStrands[Forward] = &read;
  bothQueryStrands[Reverse] = &readRC;


  //
  // Possibly use banded dynamic programming to refine the columns
  // of an alignment and the alignment score.
  //
  if (params.refineAlignments) {
    RefineAlignments(bothQueryStrands, genome, alignmentPtrs, params, mappingBuffers);
	}
	RemoveLowQualityAlignments(read,alignmentPtrs, params);
	RemoveOverlappingAlignments(alignmentPtrs, params);

  if (params.forPicard) {
    int a;
    for (a = 0; a < alignmentPtrs.size(); a++ ) {
      alignmentPtrs[a]->OrderGapsByType();
    }
  }

  //
  // Assign the query name and strand for each alignment.
  //

  for (i = 0; i < alignmentPtrs.size(); i++) { 
    T_AlignmentCandidate *aref = alignmentPtrs[i];
    //    aref->qStrand = aref->readIndex;
    if (aref->tStrand == 0) {
      aref->qName = read.GetName();
    }
    else {
      aref->qName = readRC.GetName();
    }
  }

  AssignRefContigLocations(alignmentPtrs, seqdb, genome);
}


void SumMismatches(SMRTSequence &read,
                   T_AlignmentCandidate &alignment,
									 int mismatchScore,
                   int fullIntvStart, int fullIntvEnd,
                   int &sum) {
  int alnStart, alnEnd;
  alignment.GetQIntervalOnForwardStrand(alnStart, alnEnd);
  int p;
  sum = 0;

  if (read.substitutionQV.Empty() == false) {
    for (p = fullIntvStart; p < alnStart; p++) {
      sum += read.substitutionQV[p];
    }
    for (p = alnEnd; p < fullIntvEnd; p++) {
      sum += read.substitutionQV[p];
    }
  }
	else {
		sum += mismatchScore * ((alnStart - fullIntvStart) + (fullIntvEnd - alnEnd));
	}
}

int FindMaxLengthAlignment(vector<T_AlignmentCandidate*> alignmentPtrs,
                           int &maxLengthIndex) {
  int i;
  int maxLength = 0;
  maxLengthIndex = -1;

  for (i = 0; i < alignmentPtrs.size(); i++) {
    int qStart, qEnd;
    alignmentPtrs[i]->GetQInterval(qStart, qEnd);
    if (qEnd - qStart > maxLength) {
      maxLengthIndex = i;
      maxLength = qEnd - qStart;
    }
  }
  return (maxLength != -1);
}

bool AlignmentsOverlap(T_AlignmentCandidate &alnA, T_AlignmentCandidate &alnB, float minPercentOverlap) {
  int alnAStart, alnAEnd, alnBStart, alnBEnd;
  bool useForwardStrand=true;
  alnA.GetQInterval(alnAStart, alnAEnd, useForwardStrand);
  alnB.GetQInterval(alnBStart, alnBEnd, useForwardStrand);
  // Look if one alignment encompasses the other
  int ovp = 0;
  if (alnAStart <= alnBStart and alnAEnd >= alnBEnd) {
    return true;
  }
  else if (alnBStart <= alnAStart and alnBEnd >= alnAEnd) {
    return true;
  }
  else {
    //
    // Look to see if the alignments overlap
    //

    if (alnAEnd >= alnBStart and alnAEnd <= alnBEnd) {
      ovp = alnAEnd - alnBStart;
    }
    else if (alnAStart >= alnBStart and alnAStart <= alnBEnd) {
      ovp = alnBEnd - alnAStart;
    }
  }
	float ovpPercent;
	if (alnAEnd - alnAStart > 0 and alnBEnd - alnBStart > 0) {
		ovpPercent = max(float(ovp) / float((alnAEnd - alnAStart)), 
										 float(ovp) /  float((alnBEnd - alnBStart)));
	}
	else {
		ovpPercent = 0;
	}
  // returns true when an overlap is found.
  return ( ovpPercent > minPercentOverlap);

}


void PartitionOverlappingAlignments(vector<T_AlignmentCandidate*> &alignmentPtrs,
                                    vector<set<int> > &partitions,
                                    float minOverlap) {
  if (alignmentPtrs.size() == 0) {
    partitions.clear();
    return;
  }

  set<int>::iterator setIt, setEnd;
  int i, p;
  bool overlapFound = false;
  for (i = 0; i < alignmentPtrs.size(); i++) {
    overlapFound = false;
    for (p = 0; p < partitions.size() and overlapFound == false; p++) {
      setEnd = partitions[p].end();
      for (setIt = partitions[p].begin(); setIt != partitions[p].end() and overlapFound == false; ++setIt) {
        if (AlignmentsOverlap(*alignmentPtrs[i], *alignmentPtrs[*setIt], minOverlap) or
            (alignmentPtrs[i]->QAlignStart() <= alignmentPtrs[*setIt]->QAlignStart()) and
            (alignmentPtrs[i]->QAlignEnd()  > alignmentPtrs[*setIt]->QAlignEnd())) {
          partitions[p].insert(i);
          overlapFound = true;
        }
      }
    }
    // 
    // If this alignment does not overlap any other, create a
    // partition with it as the first element.
    //
    if (overlapFound == false) {
      partitions.push_back(set<int>());
      partitions[partitions.size()-1].insert(i);
    }
  }
}

void StoreMapQVs(SMRTSequence &read,
                 vector<T_AlignmentCandidate*> &alignmentPtrs, 
                 MappingParameters &params, 
                 MappingBuffers &mappingBuffers, 
                 string title="") {
  
  //
  // Only weight alignments for mapqv against eachother if they are overlapping.
  //
  int a;
  vector<set<int> > partitions; // Each set contains alignments that overlap on the read.
  DistanceMatrixScoreFunction<DNASequence, FASTQSequence> distScoreFn;
	
	params.InitializeScoreFunction(distScoreFn);
  distScoreFn.InitializeScoreMatrix(SMRTLogProbMatrix);
	
  //
  // Rescore the alignment so that it uses probabilities. 
  //
  for (a = 0; a < alignmentPtrs.size(); a++) {
		alignmentPtrs[a]->probScore = -ComputeAlignmentScore(*alignmentPtrs[a],
																												 alignmentPtrs[a]->qAlignedSeq,
																												 alignmentPtrs[a]->tAlignedSeq,
																												 distScoreFn, params.affineAlign) / 10.0;
  }
  PartitionOverlappingAlignments(alignmentPtrs, partitions, params.minFractionToBeConsideredOverlapping);
  
  int p;
  set<int>::iterator partIt, partEnd;
  
  //
  // For each partition, store where on the read it begins, and where
  // it ends. 
  //
  vector<int> partitionBeginPos, partitionEndPos, partitionScore;
	
  partitionBeginPos.resize(partitions.size(), -1);
  partitionEndPos.resize(partitions.size(), -1);
	partitionScore.resize(partitions.size(), 0);
  vector<char> assigned;
  assigned.resize( alignmentPtrs.size());
  fill(assigned.begin(), assigned.end(), false);
	
  for (p = 0; p < partitions.size(); p++) {
    partEnd = partitions[p].end();
    int alnStart, alnEnd;

    if (partitions[p].size() > 0) {
      partIt = partitions[p].begin();
      alignmentPtrs[*partIt]->GetQInterval(alnStart, alnEnd);
      partitionBeginPos[p] = alnStart; 
      partitionEndPos[p]   = alnEnd;
			partitionScore[p] = alignmentPtrs[*partIt]->nMatch;
			
      ++partIt;
      partEnd = partitions[p].end();
      for (; partIt != partEnd; ++partIt) {
        alignmentPtrs[*partIt]->GetQInterval(alnStart, alnEnd);
				// 
        if (alnEnd - alnStart > partitionEndPos[p] - partitionBeginPos[p] and alignmentPtrs[*partIt]->nMatch *1.20 >= partitionScore[p]) {
          partitionBeginPos[p] = alnStart;
          partitionEndPos[p]   = alnEnd;
					partitionScore[p]    = alignmentPtrs[*partIt]->nMatch;
        }
      }
    }
  }
  
  //
  // For each partition, determine the widest parts of the read that
  // are aligned in the partition.  All alignments will be extended to
  // the end of the widest parts of the partition.
  //
  const static bool convertToForwardStrand = true;

  UInt i; 

  //
  // For now, just use the alignment score as the probability score.
  // Although it is possible to use the full forward probability, for
  // the most part it is pretty much the same as the Vitterbi
  // probability, but it takes a lot longer to compute.
  //

  //
  // Now estimate what the alignment scores would be if they were
  // extended past the ends of their current alignment.
  //

  for (p = 0; p < partitions.size(); p++) {
    partEnd = partitions[p].end();
    int alnStart, alnEnd;
    for (partIt = partitions[p].begin(); partitions[p].size() > 0 and partIt != partEnd; ++partIt) {
      int frontSum = 0, endSum = 0;
      int alnStart, alnEnd;
      int mismatchSum = 0;
      alignmentPtrs[*partIt]->GetQInterval(alnStart, alnEnd, convertToForwardStrand);

      if (alnStart - partitionBeginPos[p] > MAPQV_END_ALIGN_WIGGLE or
          partitionEndPos[p] - alnEnd > MAPQV_END_ALIGN_WIGGLE) {
        SumMismatches(read, *alignmentPtrs[*partIt], 15,
											partitionBeginPos[p], partitionEndPos[p], mismatchSum);
      }
      //
      // Random sequence can be aligned with about 50% similarity due
      // to optimization, so weight the qv sum 
      //
      alignmentPtrs[*partIt]->probScore += -(mismatchSum) * 0.5;
    }
  }

  //                                           
  // Determine mapqv by summing qvscores in partitions

  float mapQVDenominator = 0;
  for (p = 0; p < partitions.size(); p++) {
    set<int>::iterator nextIt;
    if (partitions[p].size() == 0) {
      continue;
    }
    int index = *partitions[p].begin();

    mapQVDenominator = alignmentPtrs[index]->probScore;
		int maxNMatch = alignmentPtrs[index]->nMatch;
		
    if (partitions[p].size() > 1) {
      partIt  = partitions[p].begin();
      partEnd = partitions[p].end();
      ++partIt;

      for (; partIt != partEnd; ++partIt) {
        index = *partIt;
				if (alignmentPtrs[index]->nMatch > maxNMatch) {
					maxNMatch = alignmentPtrs[index]->nMatch;
				}
				if (alignmentPtrs[index]->nMatch * 1.2 >= maxNMatch) {
					mapQVDenominator = LogSumOfTwo(mapQVDenominator, alignmentPtrs[index]->probScore);
				}
      }
    }
    
    
    for (partIt = partitions[p].begin(); 
         partIt != partitions[p].end(); ++partIt) {
      //
      // If only one alignment is found, assume maximum mapqv.
      //
      assigned[*partIt] = true;
      if (partitions[p].size() == 1) {
        alignmentPtrs[*partIt]->mapQV = MAX_PHRED_SCORE;
      }
      
      //
      // Look for overflow.
      //
      else if (alignmentPtrs[*partIt]->probScore - mapQVDenominator < -20) {
        alignmentPtrs[*partIt]->mapQV = 0;
      }
      else {
        double log10 = log(10);
        double sub   = alignmentPtrs[*partIt]->probScore - mapQVDenominator;
        double expo = exp(log10*sub);
        double diff = 1.0 - expo;
        int phredValue;
        
        if (expo == 0) {
          phredValue = 0;
        }
        else if (diff == 0) {
          phredValue = MAX_PHRED_SCORE;
        }
        else {
          phredValue = Phred(diff);
        }
        if (phredValue > MAX_PHRED_SCORE) {
          phredValue = MAX_PHRED_SCORE;
        }
        if (phredValue < 0) {
						phredValue = 0;
				}
        alignmentPtrs[*partIt]->mapQV = phredValue;
        assigned[*partIt]=true;
      }
    }
  }

  for (i = 0; i < assigned.size(); i++) {
    assert(assigned[i]);
  }
}

                            

//
// The full read is not the subread, and does not have masked off characters.
//
void PrintAlignment(T_AlignmentCandidate &alignment, SMRTSequence &fullRead, MappingParameters &params, AlignmentContext &alignmentContext, ostream &outFile) {
  if (alignment.score > params.maxScore) {
		if (params.verbosity > 0) {
			cout << "Not using " << alignment.qAlignedSeqPos << " " << alignment.tAlignedSeqPos << " because score: " << alignment.score << " is too low (" << params.maxScore  << ")" << endl;
		}
    return;
  }
  if (alignment.pctSimilarity < params.minPctIdentity) {
		if (params.verbosity > 0) {
			cout << "Not using " << alignment.qAlignedSeqPos << " " << alignment.tAlignedSeqPos << " because identity: " << alignment.pctSimilarity << " is too low (" << params.minPctIdentity  << ")" << endl;
		}
    return;
  }
	if (alignment.tAlignedSeq.length < params.minAlignLength) {
		if (params.verbosity > 0) {
			cout << "Not using " << alignment.qAlignedSeqPos << " " << alignment.tAlignedSeqPos << " because length: " << alignment.tAlignedSeq.length << " is too short (" << params.minAlignLength  << ")" << endl;
		}
		return;
	}
	if (alignment.mapQV < params.minMapQV) {
		return;
	}
  try {
    int lastBlock = alignment.blocks.size() - 1;
    if (params.printFormat == StickPrint) {
      PrintAlignmentStats(alignment, outFile);
      StickPrintAlignment(alignment,
                          (DNASequence&) alignment.qAlignedSeq,
                          (DNASequence&) alignment.tAlignedSeq,
                          outFile,
                          alignment.qAlignedSeqPos, alignment.tAlignedSeqPos);
    }
    else if (params.printFormat == SAM) {
      SAMOutput::PrintAlignment(alignment, fullRead, outFile, alignmentContext, params.samQVList, params.clipping, params.passthrough);
    }
    else if (params.printFormat == CompareXML) {
      CompareXMLPrintAlignment(alignment,
                               (DNASequence&) alignment.qAlignedSeq, (DNASequence&) alignment.tAlignedSeq,
                               outFile,
                               alignment.qAlignedSeqPos, alignment.tAlignedSeqPos);
    }
    else if (params.printFormat == Vulgar) {
      PrintAlignmentStats(alignment, outFile);
      VulgarPrintAlignment(alignment, outFile);
    }
    else if (params.printFormat == CompareSequencesParsable) {
      PrintCompareSequencesAlignment(alignment, alignment.qAlignedSeq, alignment.tAlignedSeq, outFile);
    }
    else if (params.printFormat == Interval) {
      if (alignment.blocks.size() > 0) {
        IntervalAlignmentPrinter::Print(alignment, outFile);
      }
    }
    else if (params.printFormat == SummaryPrint) {
      if (alignment.blocks.size() > 0) {
        SummaryAlignmentPrinter::Print(alignment, outFile);
      }
    }
  }
  catch (ostream::failure f) {
    cout << "ERROR writing to output file. The output drive may be full, or you  " << endl;
    cout << "may not have proper write permissions." << endl;
    exit(1);
  }
  
}


void PrintAlignments(vector<T_AlignmentCandidate*> alignmentPtrs,
                     SMRTSequence &read,
                     MappingParameters &params, ostream &outFile, 
                     AlignmentContext alignmentContext) {
  //
  // Print all alignments, unless the parameter placeRandomly is set.
  // In this case only one read is printed and it is selected from all
  // equally top scoring hits.
  //
  
  UInt i;
  int optScore;
  int nOpt = 0;
	if (params.verbosity > 0) {
		cout <<  "Printing " << alignmentPtrs.size() << " alignments." << endl;
	}
  if (params.placeRandomly) {
    if (alignmentPtrs.size() > 0) {
      optScore = alignmentPtrs[0]->score;
      nOpt = 1;
      // First find the minimum score, and count how many times it
      // exists
      for (i = 1; i < alignmentPtrs.size(); i++) { 
        if (alignmentPtrs[i]->score < optScore) {
          optScore = alignmentPtrs[i]->score;
          nOpt = 1;
        }
        else if (alignmentPtrs[i]->score == optScore) {
          nOpt++;
        }
      }
    }
  }

	int prev=totalBases / 10000000;
	totalBases += read.length;
	totalReads += 1;
	if (totalBases / 10000000 > prev) {
		cerr << "Processed " << totalReads << " (" << totalBases << " bp)" << endl;
	}
  
  if (params.nProc > 1) {
#ifdef __APPLE__
    sem_wait(semaphores.writer);
#else
    sem_wait(&semaphores.writer);


#endif
  }

  int optIndex = 0;
  int startIndex = 0;
  int endIndex = 0;
  if (params.placeRandomly) {
    startIndex = RandomInt(nOpt);
    assert(startIndex < nOpt);
    endIndex = startIndex + 1;
  }
  else {
    startIndex = 0;
    endIndex   = MIN(params.nBest, alignmentPtrs.size());
  }
  for (i = startIndex; i < endIndex; i++) { 
    T_AlignmentCandidate *aref = alignmentPtrs[i];      
      
    if (aref->blocks.size() == 0) {

      //
      // If the SDP alignment finds nothing, there will be no
      // blocks.  This may happen if the sdp block size is larger
      // than the anchor size found with the suffix array.  When no
      // blocks are found there is no alignment, so zero-out the
      // score and continue.
      //
      aref->score = 0;
			if (params.verbosity > 0) {
				cout << "Zero blocks found for " << aref->qName << " " << aref->qAlignedSeqPos << " " << aref->tAlignedSeqPos << endl;
			}
      continue;
    }
    
    //
    // Configure some of the alignment context before printing.
    //
    if (i > 0 and params.placeRandomly == false) {
      alignmentContext.isPrimary = false;
    }
    else {
      alignmentContext.isPrimary = true;
    }
    
    PrintAlignment(*alignmentPtrs[i], read, params, alignmentContext, outFile);
  }

  if (params.nProc > 1) {
#ifdef __APPLE__
    sem_post(semaphores.writer);
#else
    sem_post(&semaphores.writer);
#endif
  }

}

template<typename T_Sequence>
bool GetNextReadThroughSemaphore(ReaderAgglomerate &reader, MappingParameters &params, T_Sequence &read, AlignmentContext &context) {

  //
  // Grab the value of the semaphore for debugging purposes.
  //
  // uncomment when needed
  // int semvalue;
  // if (params.nProc > 1) {
  //  sem_getvalue(&semaphores.reader, &semvalue);
  // }

  //
  // Wait on a semaphore
  if (params.nProc > 1) {
#ifdef __APPLE__
    sem_wait(semaphores.reader);
#else
    sem_wait(&semaphores.reader);
#endif
  }

  bool returnValue = true;
  //
  // CCS Reads are read differently from other reads.  Do static casting here
  // of this.
  //
  if (reader.GetNext(read) == 0) {
    returnValue = false;
  }

  //
  // Set the read group id before releasing the semaphore, since other
  // threads may change the reader object to a new read group before
  // sending this alignment out to printing. 
  context.readGroupId = reader.readGroupId;
  
  if (params.nProc > 1) {
#ifdef __APPLE__
    sem_post(semaphores.reader);
#else
    sem_post(&semaphores.reader);
#endif
  }


  return returnValue;
}

void AssignMapQV(vector<T_AlignmentCandidate*> &alignmentPtrs) {
  int i;
  int mapQV = 1;
  if (alignmentPtrs.size() > 1 and alignmentPtrs[0]->score == alignmentPtrs[1]->score) {
    // the top two alignments have the same score, don't consider them as mapped.
    mapQV = 0;
  }
  
  for (i = 0; i < alignmentPtrs.size(); i++) {
    alignmentPtrs[i]->mapQV = mapQV;
  }
}

//template<typename T_SuffixArray, typename T_GenomeSequence, typename T_Tuple>
void MapReads(MappingData<T_SuffixArray, T_GenomeSequence, T_Tuple> *mapData) { 
  
  //
  // Step 1, initialize local pointers to map data 
  // for programming shorthand.
  //

  MappingParameters params = mapData->params;


  DNASuffixArray sarray;
  TupleCountTable<T_GenomeSequence, DNATuple> ct;
  SequenceIndexDatabase<FASTQSequence> seqdb;
  FASTQSequence fastaGenome;
  T_GenomeSequence    genome;
  BWT *bwtPtr;
  
  mapData->ShallowCopySuffixArray(sarray);
  mapData->ShallowCopyReferenceSequence(genome);
  mapData->ShallowCopySequenceIndexDatabase(seqdb);
  mapData->ShallowCopyTupleCountTable(ct);

  bwtPtr = mapData->bwtPtr;
  SeqBoundaryFtr<FASTQSequence> seqBoundary(&seqdb);

  VectorIndex i, j;
  if (params.match != 0) {
  }    

  int numAligned = 0;
  
  SMRTSequence smrtRead, smrtReadRC;
  SMRTSequence unrolledReadRC;
  CCSSequence  ccsRead;
  RegionAnnotation annotation;
  T_Sequence read;
  int readIndex = -1;
  int readsFileIndex;
  bool readIsCCS = false;

  //
  // Reuse the following buffers during alignment.  Since these keep
  // storage contiguous, hopefully this will decrease memory
  // fragmentation.
  //
  MappingBuffers mappingBuffers;
  while (true) {

    //
    // Scan the next read from input.  This may either be a CCS read,
    // or regular read (though this may be aligned in whole, or by
    // subread).
    //

    AlignmentContext alignmentContext;
    if (mapData->reader->GetFileType() == HDFCCS) {
      if (GetNextReadThroughSemaphore(*mapData->reader, params, ccsRead, alignmentContext) == false) {
        break;
      }
      else {
        if (params.unrollCcs == false) {
          readIsCCS = true;
          smrtRead.Copy(ccsRead);
          ccsRead.zmwData = smrtRead.zmwData = ccsRead.unrolledRead.zmwData;
          ccsRead.SetQVScale(params.qvScaleType);
        }
        else {
          smrtRead.Copy(ccsRead.unrolledRead);
        }
        ++readIndex;
      }
    }
    else {
      if (GetNextReadThroughSemaphore(*mapData->reader, params, smrtRead, alignmentContext) == false) {
        break;
      }
      else {
        ++readIndex;
        smrtRead.SetQVScale(params.qvScaleType);
      }
    }

    //
    // Only normal (non-CCS) reads should be masked.  Since CCS reads store the raw read, that is masked.
    //
    bool readHasGoodRegion = true;
    if (params.useRegionTable and params.useHQRegionTable) {
      if (readIsCCS) {
        readHasGoodRegion = MaskRead(ccsRead.unrolledRead, ccsRead.unrolledRead.zmwData, *mapData->regionTablePtr);       
      }
      else {
        readHasGoodRegion = MaskRead(smrtRead, smrtRead.zmwData, *mapData->regionTablePtr);
      }
      //
      // Store the high quality start and end of this read for masking purposes when printing.
      //
      int hqStart, hqEnd;
      int score;
      LookupHQRegion(smrtRead.zmwData.holeNumber, *mapData->regionTablePtr, hqStart, hqEnd, score);
      smrtRead.lowQualityPrefix = hqStart;
      smrtRead.lowQualitySuffix = smrtRead.length - hqEnd;
    }
    else {
      smrtRead.lowQualityPrefix = 0;
      smrtRead.lowQualitySuffix = 0;
    }

    //
    // Give the opportunity to align a subset of reads.
    //
    if (params.maxReadIndex >= 0 and smrtRead.zmwData.holeNumber > params.maxReadIndex) {
      smrtRead.Free();
      break;
    }
		
		if (params.SkipRead(smrtRead.zmwData.holeNumber)) {
      if (readIsCCS) {
        ccsRead.Free();
      }
      smrtRead.Free();
			continue;
		}

    if (readHasGoodRegion == false or (params.readIndex >= 0 and smrtRead.zmwData.holeNumber != params.readIndex)) {

      //
      // Nothing to do with this read. Skip aligning it entirely.
      //
      if (readIsCCS) {
        ccsRead.Free();
      }
      smrtRead.Free();
      // Stop processing once the specified read index is reached.
      // Eventually this will change to just seek to readIndex, and
      // just align one read anyway.
      if (smrtRead.zmwData.holeNumber > params.readIndex) {
        break;
      }
      continue;
    }
    
    

    // 
    // Discard reads that are too small, or not labeled as having any
    // useable/good sequence. 
    //
    if (smrtRead.length < params.minReadLength or readHasGoodRegion == false or 
                          (params.maxReadLength != 0 and 
                           smrtRead.length > params.maxReadLength)) {
      if (readIsCCS) {
        ccsRead.Free();
      }
      smrtRead.Free();
      continue;
    }

    if (smrtRead.qual.Empty() == false and smrtRead.GetAverageQuality() < params.minAvgQual) {
      if (readIsCCS) {
        ccsRead.Free();
      }
      smrtRead.Free();
      continue;
    }

    smrtRead.MakeRC(smrtReadRC);
    
    if (readIsCCS) {
      ccsRead.unrolledRead.MakeRC(unrolledReadRC);
    }

    //
    // When aligning subreads separately, iterate over each subread, and print the alignments for these.
    //

    ReadAlignments allReadAlignments;
    allReadAlignments.read = smrtRead;

    

    if (readIsCCS == false and params.mapSubreadsSeparately) {
      vector<ReadInterval> subreadIntervals;
      //
      // Determine endpoints of this subread in the main read.
      //
      if (params.useRegionTable == false) {
        //
        // When there is no region table, the subread is the entire
        // read.
        //
        ReadInterval wholeRead(0, smrtRead.length);
        // The set of subread intervals is just the entire read.
        subreadIntervals.push_back(wholeRead);
      }
      else {
        // 
        // Grab the subread intervals from the entire region table to
        // iterate over.
        //
        CollectSubreadIntervals(smrtRead, mapData->regionTablePtr, subreadIntervals, params.byAdapter);
      }
      
      //
      // Make room for alignments.
      //
      allReadAlignments.Resize(subreadIntervals.size());
      allReadAlignments.alignMode = Subread;

      int endIndex = subreadIntervals.size();
      DNALength highQualityStartPos = smrtRead.lowQualityPrefix;
      DNALength highQualityEndPos   = smrtRead.length - smrtRead.lowQualitySuffix;
      DNALength intvIndex;
      for (intvIndex = 0; intvIndex < endIndex; intvIndex++) {
        SMRTSequence subreadSequence, subreadSequenceRC;
        //
        // There is a user specified minimum read length.  Don't
        // bother with reads shorter than this.
        //

        if ( (subreadIntervals[intvIndex].end < highQualityStartPos) or
              (subreadIntervals[intvIndex].start > highQualityEndPos) ) {
          continue;
        }
        
        if (subreadIntervals[intvIndex].end - subreadIntervals[intvIndex].start < params.minSubreadLength) {
          continue;
        }
        
        //
        // subreadMapType is a way of limiting the portion of the read
        // that is aligned.  The output is similar, but the
        // computation is slightly different.  The subreadMapType 0
        // was written first, and just creates a hard mask over the
        // regions that are not to be aligned.  The subreadMapType is
        // slightly more formal mode where a new read is pointed to
        // the subread then aligned.
        //
        // subreadMapType of 0 is always used, however eventually it
        // may be faster to go to 1, just 1 isn't tested thoroughly
        // yet. 
        // 
        // Note, for proper SAM printing, subreadMaptype of 0 is needed.
        //
        if (params.subreadMapType == 0) {
          smrtRead.MakeSubreadAsMasked(subreadSequence, subreadIntervals[intvIndex].start, subreadIntervals[intvIndex].end);
        }
        else if (params.subreadMapType == 1) {
          smrtRead.MakeSubreadAsReference(subreadSequence, subreadIntervals[intvIndex].start, subreadIntervals[intvIndex].end);
        }
        
        //
        // Trim the boundaries of the read so that only the hq regions are aligned, and no N's.
        //
        if ((subreadSequence.subreadStart < highQualityStartPos) and
             (subreadSequence.subreadEnd   > highQualityStartPos)) {
          subreadSequence.subreadStart = highQualityStartPos;
        }
        if (subreadSequence.subreadStart  < highQualityEndPos and
            subreadSequence.subreadEnd  > highQualityEndPos) {
          subreadSequence.subreadEnd = highQualityEndPos;
        }
        if (!params.preserveReadTitle) {
          smrtRead.SetSubreadTitle(subreadSequence, subreadSequence.subreadStart, subreadSequence.subreadEnd);
        }
        else {
          subreadSequence.CopyTitle(smrtRead.title);
        }

        subreadSequence.MakeRC(subreadSequenceRC);
        subreadSequenceRC.subreadStart = smrtRead.length - subreadSequence.subreadEnd;
        subreadSequenceRC.subreadEnd   = smrtRead.length - subreadSequence.subreadStart;

        //
        // Store the sequence that is being mapped in case no hits are
        // found, and missing sequences are printed.
        //
        allReadAlignments.SetSequence(intvIndex, subreadSequence);

        vector<T_AlignmentCandidate*> alignmentPtrs;
        mapData->metrics.numReads++;
        //
        // Try default and fast parameters to map the read.
        //
        subreadSequence.zmwData.holeNumber = smrtRead.zmwData.holeNumber;

        MapRead(subreadSequence, subreadSequenceRC, 
                genome,           // possibly multi fasta file read into one sequence
                sarray, *bwtPtr,  // The suffix array, and the bwt-fm index structures
                seqBoundary,      // Boundaries of contigs in the
                                  // genome, alignments do not span
                                  // the ends of boundaries.
                ct,               // Count table to use word frequencies in the genome to weight matches.
                seqdb,            // Information about the names of
                                  // chromosomes in the genome, and
                                  // where their sequences are in the genome.
                params,           // A huge list of parameters for
                                  // mapping, only compile/command
                                  // line values set.
                mapData->metrics, // Keep track of time/ hit counts,
                                  // etc.. Not fully developed, but
                                  // should be.
                alignmentPtrs,    // Where the results are stored.
                mappingBuffers,   // A class of buffers for structurs
                                  // like dyanmic programming
                                  // matrices, match lists, etc., that are not
                                  // reallocated between calls to
                                  // MapRead.  They are cleared though.
                mapData           // Some values that are shared
                                  // across threads.
                );

        //
        // No alignments were found, sometimes parameters are
        // specified to try really hard again to find an alignment.
        // This sets some parameters that use a more sensitive search
        // at the cost of time.
        //

        if ((alignmentPtrs.size() == 0 or alignmentPtrs[0]->pctSimilarity < 80) and params.doSensitiveSearch) {
          MappingParameters sensitiveParams = params;
          sensitiveParams.SetForSensitivity();
          MapRead(subreadSequence, subreadSequenceRC, genome, 
                  sarray, *bwtPtr, 
                  seqBoundary, ct, seqdb,
                  sensitiveParams, mapData->metrics, 
                  alignmentPtrs, mappingBuffers, 
                  mapData);
        }

        //
        // Store the mapping quality values.
        //
        if (alignmentPtrs.size() > 0 and 
            alignmentPtrs[0]->score < params.maxScore and 
            params.storeMapQV) {
          StoreMapQVs(subreadSequence, alignmentPtrs, params, mappingBuffers, subreadSequence.title);
        }
        
        allReadAlignments.AddAlignmentsForSeq(intvIndex, alignmentPtrs);


        //
        // Move reference from subreadSequence, which will be freed at
        // the end of this loop to the smrtRead, which exists for the
        // duration of aligning all subread of the smrtRead.
        //
        if (params.subreadMapType == 0) {
          int a;
          for (a = 0; a < alignmentPtrs.size(); a++) {
            if (alignmentPtrs[a]->qStrand == 0) {
              alignmentPtrs[a]->qAlignedSeq.ReferenceSubstring(smrtRead,
                                                               alignmentPtrs[a]->qAlignedSeq.seq - subreadSequence.seq,
                                                               alignmentPtrs[a]->qAlignedSeqLength);
            }
            else {
              alignmentPtrs[a]->qAlignedSeq.ReferenceSubstring(smrtReadRC,
                                                               alignmentPtrs[a]->qAlignedSeq.seq - subreadSequenceRC.seq, 
                                                               alignmentPtrs[a]->qAlignedSeqLength);
            }
          }
          subreadSequence.Free();
          subreadSequenceRC.Free();
        }

      }
    }
    else {
      //
      // The read The read must be mapped as a whole, even it it contains subreads.
      //
      vector<T_AlignmentCandidate*> alignmentPtrs;
      mapData->metrics.numReads++;
      smrtRead.subreadStart   = 0; smrtRead.subreadEnd   = smrtRead.length;
      smrtReadRC.subreadStart = 0; smrtReadRC.subreadEnd = smrtRead.length;

      MapRead(smrtRead, smrtReadRC, 
              genome, sarray, *bwtPtr, seqBoundary, ct, seqdb, params, mapData->metrics,
              alignmentPtrs, mappingBuffers, mapData);

      //
      // Just one sequence is aligned.  There is one primary hit, and
      // all other are secondary.
      //

      if (readIsCCS == false or params.useCcsOnly) {
        //
        // Record some information for proper SAM Annotation.
        //
        allReadAlignments.Resize(1);
        allReadAlignments.AddAlignmentsForSeq(0, alignmentPtrs);
        if (params.useCcsOnly) {
          allReadAlignments.alignMode = CCSDeNovo;
        }
        else {
          allReadAlignments.alignMode = Fullread;
        }
        allReadAlignments.SetSequence(0, smrtRead);
      }
      else if (readIsCCS) {
        //
        // Align ccs reads.
        //

        //
        // Align the ccs subread to where the denovo sequence mapped (explode).
        //
        SMRTSequence readRC;
        
        CCSIterator ccsIterator;
        FragmentCCSIterator fragmentCCSIterator;
        CCSIterator *subreadIterator;
        
        //
        // Choose a different iterator over subreads depending on the
        // alignment mode.  When the mode is allpass, include the
        // framgents that are not necessarily full pass.
        //
        if (params.useAllSubreadsInCcs) {
          // 
          // Use all subreads even if they are not full pass
          fragmentCCSIterator.Initialize(&ccsRead, mapData->regionTablePtr);
          subreadIterator = &fragmentCCSIterator;
          allReadAlignments.alignMode = CCSAllPass;
        }
        else {
          // Use only full pass reads.
          ccsIterator.Initialize(&ccsRead);
          subreadIterator = &ccsIterator;
          allReadAlignments.alignMode = CCSFullPass;
        }

        allReadAlignments.Resize(subreadIterator->GetNumPasses());
          
        int passDirection, passStartBase, passNumBases;
        SMRTSequence subread;
        
        //
        // The read was previously set to the smrtRead, which was the
        // de novo ccs sequence.  Since the alignments of exploded
        // reads are reported, the unrolled read should be used as the
        // reference when printing.
        //
        allReadAlignments.read = ccsRead.unrolledRead;
        subreadIterator->Reset();
        int subreadIndex;

        //
        // Realign all subreads.
        //
        for (subreadIndex = 0; subreadIndex < subreadIterator->GetNumPasses(); subreadIndex++) {
          subreadIterator->GetNext(passDirection, passStartBase, passNumBases);
          if (passNumBases <= 0) { continue; }

          subread.ReferenceSubstring(ccsRead.unrolledRead, passStartBase, passNumBases-1);
          subread.CopyTitle(ccsRead.title);
          //allReadAlignments.SetSequence(subreadIndex, subread);
          // The unrolled alignment should be relative to the entire read.
          allReadAlignments.SetSequence(subreadIndex, ccsRead.unrolledRead);
		  
          int alignmentIndex;
          
          //
          // Align all subreads to the positions that the de novo
          // sequence has aligned to.
          //
          for (alignmentIndex = 0; alignmentIndex < alignmentPtrs.size(); alignmentIndex++) {

            T_AlignmentCandidate *alignment = alignmentPtrs[alignmentIndex];          
            if (alignment->score > params.maxScore) break;
            //
            // Determine where in the genome the subread has mapped.
            //
            DNASequence ccsAlignedForwardRefSequence, ccsAlignedReverseRefSequence;
            
            if (alignment->tStrand == 0) {
              // This needs to be changed -- copy copies RHS into LHS,
              // CopyAsRC copies LHS into RHS
              ccsAlignedForwardRefSequence.Copy(alignment->tAlignedSeq);
              alignment->tAlignedSeq.CopyAsRC(ccsAlignedReverseRefSequence);
            }
            else {
              alignment->tAlignedSeq.CopyAsRC(ccsAlignedForwardRefSequence);
              ccsAlignedReverseRefSequence.Copy(alignment->tAlignedSeq);
            }

            DistanceMatrixScoreFunction<DNASequence, FASTQSequence> distScoreFn;
            distScoreFn.del = params.deletion;
            distScoreFn.ins = params.insertion;
            distScoreFn.InitializeScoreMatrix(SMRTDistanceMatrix);
          
            /*
             * Determine the strand to align the subread strand to.
             */

            if (alignment->tStrand == passDirection) {
              //
              // The alignment is in the same direction as the subread.  Easy.
              //

              T_AlignmentCandidate exploded;
              int explodedScore;
              bool computeProbIsFalse = false; 
              // need to make this more compact...
              explodedScore = GuidedAlign(subread, ccsAlignedForwardRefSequence,
                                          distScoreFn, 12,
                                          params.sdpIns, params.sdpDel, params.indelRate,
                                          mappingBuffers,
                                          exploded,  
                                          // Use a small tuple size
                                          Local, computeProbIsFalse, 6);

              if (params.verbosity > 0) {
                StickPrintAlignment(exploded,
                                    (DNASequence&) subread,
                                    (DNASequence&) ccsAlignedForwardRefSequence,
                                    cout,
                                    exploded.qAlignedSeqPos, exploded.tAlignedSeqPos);
              }
                
              if (exploded.blocks.size() > 0) {
                ComputeAlignmentStats(exploded, subread.seq, ccsAlignedForwardRefSequence.seq, distScoreFn); //SMRTDistanceMatrix, params.indel, params.indel );
                if (exploded.score <= params.maxScore) {

                  //
                  // Reset the coordinates of the alignment so that
                  // they are relative to the genome, not the aligned
                  // substring.
                  //
                  exploded.qStrand = 0;
                  exploded.tStrand = alignment->tStrand;

                  exploded.qLength = ccsRead.unrolledRead.length;
                  exploded.tLength = alignment->tLength;
                  exploded.tAlignedSeq.Copy(ccsAlignedForwardRefSequence);
                  exploded.tAlignedSeqPos = alignment->tAlignedSeqPos;
                  exploded.tAlignedSeqLength = alignment->tAlignedSeqLength;
                  exploded.qAlignedSeq.ReferenceSubstring(subread);
                  exploded.qAlignedSeqPos = passStartBase;
                  exploded.qAlignedSeqLength = passNumBases;
                  exploded.mapQV = alignment->mapQV;
                
                  //
                  // Assign the name of this subread.
                  //
                  exploded.tName = alignment->tName;
                  stringstream namestrm;
                  namestrm << "/" << passStartBase << "_" << passStartBase + passNumBases;
                  exploded.qName = string(ccsRead.unrolledRead.title) + namestrm.str();
                
                  //
                  // Assign the proper chromosome coordiantes.
                  //
                  AssignRefContigLocation(exploded, seqdb, genome);

                  //
                  // Save this alignment for printing later.
                  //
                  T_AlignmentCandidate *alignmentPtr = new T_AlignmentCandidate;
                  *alignmentPtr = exploded;
                  allReadAlignments.AddAlignmentForSeq(subreadIndex, alignmentPtr);
                }
              }
            }
            else {
              int pos = 0;


              T_AlignmentCandidate explodedrc;
              int explodedRCScore;
              bool computeProbIsFalse = false;
              explodedRCScore = GuidedAlign(subread, ccsAlignedReverseRefSequence,
                                            distScoreFn, 10,
                                            params.sdpIns, params.sdpDel, params.indelRate,
                                            mappingBuffers,
                                            explodedrc,  
                                            Global, computeProbIsFalse, 4);


              int explodedrcscore = explodedrc.score;
              if (params.verbosity > 0) {
                cout << "subread: " << endl;
                ((DNASequence) subread).PrintSeq(cout);
                cout << endl;
                cout << "ccsAlignedSequence " << endl;
                ((DNASequence) ccsAlignedReverseRefSequence).PrintSeq(cout);
                StickPrintAlignment(explodedrc,
                                    (DNASequence&) subread,
                                    (DNASequence&) ccsAlignedReverseRefSequence,
                                    cout,
                                    explodedrc.qAlignedSeqPos, explodedrc.tAlignedSeqPos);
              }
              if (explodedrc.blocks.size() > 0) {
                ComputeAlignmentStats(explodedrc, subread.seq, ccsAlignedReverseRefSequence.seq, distScoreFn);//SMRTDistanceMatrix, params.indel, params.indel );

                if (explodedrc.score <= params.maxScore) {
                  explodedrc.qStrand = 0;
                  explodedrc.tStrand = 1;
                  explodedrc.qLength = ccsRead.unrolledRead.length;
                  explodedrc.tLength = alignment->tLength;
                  explodedrc.tAlignedSeq.Copy(ccsAlignedReverseRefSequence);
                  explodedrc.qAlignedSeq.ReferenceSubstring(subread);
                  explodedrc.tAlignedSeqPos = genome.MakeRCCoordinate(alignment->tAlignedSeqPos + alignment->tAlignedSeqLength);
                  explodedrc.tAlignedSeqLength = alignment->tAlignedSeqLength;
                  explodedrc.qAlignedSeqPos = passStartBase;
                  explodedrc.mapQV = alignment->mapQV;

                  explodedrc.tName = alignment->tName;
                  stringstream namestrm;
                  namestrm << "/" << passStartBase << "_" << passStartBase + passNumBases;
                  explodedrc.qName = string(ccsRead.unrolledRead.title) + namestrm.str();

                
                  AssignRefContigLocation(explodedrc, seqdb, genome);
                
                  T_AlignmentCandidate *alignmentPtr = new T_AlignmentCandidate;
                  *alignmentPtr = explodedrc;
                  allReadAlignments.AddAlignmentForSeq(subreadIndex, alignmentPtr);
                }
              }
            }
          }
        }
			}
		}

    int subreadIndex;
    int nAlignedSubreads = allReadAlignments.GetNAlignedSeq();

    //
    // Initialize the alignemnt context with information applicable to SAM output.
    //
    alignmentContext.alignMode = allReadAlignments.alignMode;
    for (subreadIndex = 0; subreadIndex < nAlignedSubreads; subreadIndex++) {
      if (allReadAlignments.subreadAlignments[subreadIndex].size() > 0) {
        alignmentContext.numProperlyAlignedSubreads++;
      }
    }
    if (alignmentContext.numProperlyAlignedSubreads == allReadAlignments.subreadAlignments.size()) {
      alignmentContext.allSubreadsProperlyAligned = true;
    }
    alignmentContext.nSubreads = nAlignedSubreads;

    for (subreadIndex = 0; subreadIndex < nAlignedSubreads; subreadIndex++) {
      alignmentContext.subreadIndex = subreadIndex;
      if (subreadIndex < nAlignedSubreads-1 and allReadAlignments.subreadAlignments[subreadIndex+1].size() > 0) {
        alignmentContext.nextSubreadPos = allReadAlignments.subreadAlignments[subreadIndex+1][0]->QAlignStart();
        alignmentContext.nextSubreadDir = allReadAlignments.subreadAlignments[subreadIndex+1][0]->qStrand;
        alignmentContext.rNext = allReadAlignments.subreadAlignments[subreadIndex+1][0]->tName;
        alignmentContext.hasNextSubreadPos = true;
      }
      else {
        alignmentContext.nextSubreadPos = 0;
        alignmentContext.nextSubreadDir = 0;
        alignmentContext.rNext = "";
        alignmentContext.hasNextSubreadPos = false;
      }
        
      if (allReadAlignments.subreadAlignments[subreadIndex].size() > 0) {
        int alnIndex;

        PrintAlignments(allReadAlignments.subreadAlignments[subreadIndex], 
//                        smrtRead, // the source read
                        allReadAlignments.subreads[subreadIndex], // the source read
                        // for these alignments
                        params, *mapData->outFilePtr,
                        alignmentContext);   
      }
      else {
        //
        // Print the unaligned sequences.
        //
        if (params.printFormat == SAM) {
          if (params.nProc == 1) {
            SAMOutput::PrintUnalignedRead(allReadAlignments.subreads[subreadIndex], *mapData->outFilePtr, alignmentContext, params.samQVList, params.clipping);
          }
          else {
#ifdef __APPLE__
            sem_wait(semaphores.writer);
#else
            sem_wait(&semaphores.writer);
#endif
            SAMOutput::PrintUnalignedRead(allReadAlignments.subreads[subreadIndex], *mapData->outFilePtr, alignmentContext, params.samQVList, params.clipping);
#ifdef __APPLE__
            sem_post(semaphores.writer);
#else
            sem_post(&semaphores.writer);
#endif
          }
        }

        if (params.printUnaligned == true) {
          if (params.nProc == 1) {
            allReadAlignments.subreads[subreadIndex].PrintSeq(*mapData->unalignedFilePtr);
          }
          else {
#ifdef __APPLE__
            sem_wait(semaphores.unaligned);
#else
            sem_wait(&semaphores.unaligned);
#endif
            allReadAlignments.subreads[subreadIndex].PrintSeq(*mapData->unalignedFilePtr);
#ifdef __APPLE__
            sem_post(semaphores.unaligned);
#else
            sem_post(&semaphores.unaligned);
#endif
          }
        }
      }

    }
    
    allReadAlignments.Clear();
		smrtReadRC.Free();
		smrtRead.Free();

		if (readIsCCS) {
			ccsRead.Free();
      unrolledReadRC.Free();
		}
    numAligned++;
    if(numAligned % 200 == 0) {
      mappingBuffers.Reset();
    }
	}
	if (params.nProc > 1) {
#ifdef __APPLE__
		sem_wait(semaphores.reader);
		sem_post(semaphores.reader);
#else
		sem_wait(&semaphores.reader);
		sem_post(&semaphores.reader);
#endif
	}
	if (params.nProc > 1) {
		pthread_exit(NULL); 
	}
}

float ComputePMatch(float accuracy, int anchorLength) {
  assert(anchorLength >= 0);
  if (anchorLength == 0) { 
    return 0;
  }
  else {
    return pow(accuracy,anchorLength);
  }
}

//
// Assume the number of mismatches in a row follow a geometric distribution.
//
void GeometricDistributionSummaryStats(float pSuccess,
                                       float &mean, float &variance) {
  mean = 1/pSuccess;
  variance = (1-pSuccess)/ (pow(pSuccess,2));
}

int ComputeExpectedWaitingBases(float mean, float variance, float certainty) {
  float nStdDev;
  assert(FindQNorm(certainty, nStdDev) != false);
  return mean + sqrt(variance) * nStdDev;
}




int main(int argc, char* argv[]) {
  seen=false;
	//
	// Configure parameters for refining alignments.
	//
	MappingParameters params;
	ReverseCompressIndex index;
	pid_t parentPID;
	pid_t *pids;
	
	CommandLineParser clp;
  string commandLine;
	string helpString;
	SetHelp(helpString);

	string conciseHelpString;
	SetConciseHelp(conciseHelpString);
	
	stringstream usageSStrm;
	usageSStrm << "   Basic usage: 'blasr reads.{fasta,bas.h5} genome.fasta [-options] " << endl
						 << " [option]\tDescription (default_value)." << endl << endl
						 << " Input Files." << endl
						 << "   reads.fasta is a multi-fasta file of reads.  While any fasta file is valid input, " 
						 "it is preferable to use pls.h5 or bas.h5 files because they contain "
						 "more rich quality value information." << endl
						 << "   reads.bas.h5|reads.pls.h5 Is the native output format in Hierarchical Data Format of "
		"SMRT reads. This is the preferred input to blasr because rich quality"
		"value (insertion,deletion, and substitution quality values) information is "
		"maintained.  The extra quality information improves variant detection and mapping"<<
		"speed." << endl << endl;


	clp.SetHelp(helpString);
	clp.SetConciseHelp(conciseHelpString);
	clp.SetProgramSummary(usageSStrm.str());
	clp.SetProgramName("blasr");

	//
	// Make the default arguments.
	//
	bool required=true;
	bool optional=false;
	int  trashbinInt;
	float trashbinFloat;
  bool trashbinBool;
	//	clp.RegisterStringListOption("input_files", &params.readsFileNames, "Read files followed by genome.");
	//	clp.RegisterPreviousFlagsAsHidden();
	bool printVerboseHelp = false;
	bool printLongHelp    = false;
	clp.RegisterStringOption("sa", &params.suffixArrayFileName, "");
	clp.RegisterStringOption("ctab", &params.countTableName, "" );
	clp.RegisterStringOption("regionTable", &params.regionTableFileName, "");
	clp.RegisterIntOption("bestn", (int*) &params.nBest, "", CommandLineParser::PositiveInteger);
  clp.RegisterIntOption("limsAlign", &params.limsAlign, "", CommandLineParser::PositiveInteger);
  clp.RegisterFlagOption("printOnlyBest", &params.printOnlyBest, "");
	clp.RegisterFlagOption("outputByThread", &params.outputByThread, "");
	clp.RegisterFlagOption("rbao", &params.refineBetweenAnchorsOnly, "");
  clp.RegisterFlagOption("allowAdjacentIndels", &params.forPicard, "");
  clp.RegisterFlagOption("onegap", &params.separateGaps, "");
	clp.RegisterFlagOption("overlap", &params.overlap, "");
  clp.RegisterFlagOption("allowAdjacentIndels", &params.forPicard, "");
  clp.RegisterFlagOption("placeRepeatsRandomly", &params.placeRandomly, "");
  clp.RegisterIntOption("randomSeed", &params.randomSeed, "", CommandLineParser::Integer);
	clp.RegisterFlagOption("extend", &params.extendAlignments, "");
  clp.RegisterIntOption("branchExpand", &params.anchorParameters.branchExpand, "", CommandLineParser::NonNegativeInteger);
  clp.RegisterIntOption("maxExtendDropoff", &params.maxExtendDropoff, "", CommandLineParser::NonNegativeInteger);
	clp.RegisterFlagOption("nucmer", &params.emulateNucmer, "");
	clp.RegisterIntOption("maxExpand", &params.maxExpand, "", CommandLineParser::PositiveInteger);
	clp.RegisterIntOption("minExpand", &params.minExpand, "", CommandLineParser::NonNegativeInteger);
	clp.RegisterStringOption("seqdb",  &params.seqDBName, "");
	clp.RegisterStringOption("anchors",  &params.anchorFileName, "");
  clp.RegisterStringOption("clusters", &params.clusterFileName, "");
  clp.RegisterFlagOption("noStoreMapQV", &params.storeMapQV, "");
	clp.RegisterFlagOption("noRefineAlign", (bool*) &params.refineAlign, "");
	clp.RegisterFlagOption("guidedAlign", (bool*)&params.useGuidedAlign, "");
  clp.RegisterFlagOption("useGuidedAlign", (bool*)&trashbinBool, "");
  clp.RegisterFlagOption("noUseGuidedAlign", (bool*)&params.useGuidedAlign, "");
  clp.RegisterFlagOption("header", (bool*)&params.printHeader, "");
  clp.RegisterIntOption("subreadImplType", &params.subreadMapType, "", CommandLineParser::PositiveInteger);
	clp.RegisterIntOption("bandSize", &params.bandSize, "", CommandLineParser::PositiveInteger);	
	clp.RegisterIntOption("extendBandSize", &params.extendBandSize, "", CommandLineParser::PositiveInteger);	
	clp.RegisterIntOption("guidedAlignBandSize", &params.guidedAlignBandSize, "", CommandLineParser::PositiveInteger);	
	clp.RegisterIntOption("maxAnchorsPerPosition", &params.anchorParameters.maxAnchorsPerPosition, "", CommandLineParser::PositiveInteger);
	clp.RegisterIntOption("sdpMaxAnchorsPerPosition", &params.sdpMaxAnchorsPerPosition, "", CommandLineParser::PositiveInteger);
	clp.RegisterIntOption("stopMappingOnceUnique", (int*) &params.anchorParameters.stopMappingOnceUnique, "", CommandLineParser::NonNegativeInteger);
	clp.RegisterStringOption("out", &params.outFileName, "");
	clp.RegisterIntOption("match", &params.match, "", CommandLineParser::Integer);
	clp.RegisterIntOption("mismatch", &params.mismatch, "", CommandLineParser::Integer);
	clp.RegisterIntOption("minMatch", &params.minMatchLength, "", CommandLineParser::PositiveInteger);
	clp.RegisterIntOption("maxMatch", &params.anchorParameters.maxLCPLength, "", CommandLineParser::NonNegativeInteger);
	clp.RegisterIntOption("maxLCPLength", &params.anchorParameters.maxLCPLength, "", CommandLineParser::NonNegativeInteger);
	clp.RegisterIntOption("indel", &params.indel, "", CommandLineParser::Integer);
	clp.RegisterIntOption("insertion", &params.insertion, "", CommandLineParser::Integer);
	clp.RegisterIntOption("deletion", &params.deletion, "", CommandLineParser::Integer);
	clp.RegisterIntOption("idsIndel", &params.idsIndel, "", CommandLineParser::Integer);
	clp.RegisterIntOption("sdpindel", &params.sdpIndel, "", CommandLineParser::Integer);
	clp.RegisterIntOption("sdpIns", &params.sdpIns, "", CommandLineParser::Integer);
	clp.RegisterIntOption("sdpDel", &params.sdpDel, "", CommandLineParser::Integer);
	clp.RegisterFloatOption("indelRate", &params.indelRate, "", CommandLineParser::NonNegativeFloat);
	clp.RegisterFloatOption("minRatio", &params.minRatio, "", CommandLineParser::NonNegativeFloat); 
	clp.RegisterFloatOption("sdpbypass", &params.sdpBypassThreshold, "", CommandLineParser::NonNegativeFloat);
	clp.RegisterFloatOption("minFrac", &trashbinFloat, "", CommandLineParser::NonNegativeFloat);
	clp.RegisterIntOption("maxScore", &params.maxScore, "", CommandLineParser::Integer);
	clp.RegisterStringOption("bwt", &params.bwtFileName, "");
	clp.RegisterIntOption("m", &params.printFormat, "", CommandLineParser::NonNegativeInteger);
  clp.RegisterFlagOption("sam", &params.printSAM, "");
  clp.RegisterStringOption("clipping", &params.clippingString, "");
	clp.RegisterIntOption("sdpTupleSize", &params.sdpTupleSize, "", CommandLineParser::PositiveInteger);
	clp.RegisterIntOption("sdpPrefix", &params.sdpPrefix, "", CommandLineParser::NonNegativeInteger);
	clp.RegisterIntOption("pvaltype", &params.pValueType, "", CommandLineParser::NonNegativeInteger);
	clp.RegisterIntOption("start", &params.startRead, "", CommandLineParser::NonNegativeInteger);
	clp.RegisterIntOption("stride", &params.stride, "", CommandLineParser::NonNegativeInteger);
	clp.RegisterFloatOption("subsample", &params.subsample, "", CommandLineParser::PositiveFloat);
	clp.RegisterIntOption("nproc", &params.nProc, "", CommandLineParser::PositiveInteger);
	clp.RegisterFlagOption("sortRefinedAlignments",(bool*) &params.sortRefinedAlignments, "");
	clp.RegisterIntOption("quallc", &params.qualityLowerCaseThreshold, "", CommandLineParser::Integer);
	clp.RegisterFlagOption("p", (bool*) &params.progress, "");
	clp.RegisterFlagOption("v", (bool*) &params.verbosity, "");
	clp.RegisterIntOption("V", &params.verbosity, "Specify a level of verbosity.", CommandLineParser::NonNegativeInteger);
	clp.RegisterIntOption("contextAlignLength", &params.anchorParameters.contextAlignLength, "", CommandLineParser::PositiveInteger);
	clp.RegisterFlagOption("skipLookupTable", &params.anchorParameters.useLookupTable, "");
	clp.RegisterStringOption("metrics", &params.metricsFileName, "");
	clp.RegisterStringOption("lcpBounds", &params.lcpBoundsFileName, "");
  clp.RegisterStringOption("fullMetrics", &params.fullMetricsFileName, "");
	clp.RegisterIntOption("nbranch", &params.anchorParameters.numBranches, "", CommandLineParser::NonNegativeInteger);
	clp.RegisterFlagOption("divideByAdapter", &params.byAdapter, "");
	clp.RegisterFlagOption("useQuality", &params.ignoreQualities, "");
	clp.RegisterFlagOption("noFrontAlign", &params.extendFrontAlignment, "");
	clp.RegisterIntOption("minReadLength", &params.minReadLength, "", CommandLineParser::NonNegativeInteger);
	clp.RegisterIntOption("minAlignLength", &params.minAlignLength, "", CommandLineParser::NonNegativeInteger);
	clp.RegisterIntOption("maxReadLength", &params.maxReadLength, "", CommandLineParser::NonNegativeInteger);
	clp.RegisterIntOption("minSubreadLength", &params.minSubreadLength, "", CommandLineParser::NonNegativeInteger);
	clp.RegisterIntOption("minAvgQual", &params.minAvgQual, "", CommandLineParser::Integer);
	clp.RegisterFlagOption("advanceHalf", &params.advanceHalf, "");
	clp.RegisterIntOption("advanceExactMatches", &params.anchorParameters.advanceExactMatches, "", CommandLineParser::NonNegativeInteger);
	clp.RegisterFlagOption("unrollCcs", &params.unrollCcs, "");
	clp.RegisterFlagOption("useccs", &params.useCcs, "");
	clp.RegisterFlagOption("useccsdenovo", &params.useCcsOnly, "");
	clp.RegisterFlagOption("useccsall", &params.useAllSubreadsInCcs, "");
  clp.RegisterFlagOption("extendDenovoCCSSubreads", &params.extendDenovoCCSSubreads, "");
	clp.RegisterFlagOption("noRefineAlignments", &params.refineAlignments, "");
	clp.RegisterFloatOption("minPctIdentity", &params.minPctIdentity, "", CommandLineParser::NonNegativeFloat);
	clp.RegisterFloatOption("maxPctIdentity", &params.maxPctIdentity, "", CommandLineParser::NonNegativeFloat);
	clp.RegisterIntOption("nCandidates", &params.nCandidates, "", CommandLineParser::NonNegativeInteger);
	clp.RegisterFlagOption("useTemp", (bool*) &params.tempDirectory, "");
	clp.RegisterFlagOption("noSplitSubreads", &params.mapSubreadsSeparately, "");
  clp.RegisterIntOption("subreadMapType", &params.subreadMapType, "", CommandLineParser::NonNegativeInteger);
	clp.RegisterStringOption("titleTable", &params.titleTableName, "");
	clp.RegisterFlagOption("useSensitiveSearch", &params.doSensitiveSearch, "");
	clp.RegisterFlagOption("ignoreRegions", &params.useRegionTable, "");
	clp.RegisterFlagOption("ignoreHQRegions", &params.useHQRegionTable, "");
  clp.RegisterFlagOption("computeAlignProbability", &params.computeAlignProbability, "");
	clp.RegisterStringOption("unaligned", &params.unalignedFileName, "");
  clp.RegisterFlagOption("global", &params.doGlobalAlignment, "");
	clp.RegisterIntOption("globalChainType", &params.globalChainType, "", CommandLineParser::NonNegativeInteger);
	clp.RegisterFlagOption("noPrintSubreadTitle", (bool*) &params.printSubreadTitle, "");
	clp.RegisterIntOption("saLookupTableLength", &params.lookupTableLength, "", CommandLineParser::PositiveInteger);
	clp.RegisterFlagOption("useDetailedSDP", &params.detailedSDPAlignment, "");
	clp.RegisterFlagOption("nouseDetailedSDP", &trashbinBool, "");
	clp.RegisterIntOption("sdpFilterType", &params.sdpFilterType, "", CommandLineParser::NonNegativeInteger);
	clp.RegisterIntOption("scoreType", &params.scoreType, "", CommandLineParser::NonNegativeInteger);
	clp.RegisterFlagOption("help", &params.printDiscussion, "");
	clp.RegisterFlagOption("h", &printVerboseHelp, "");
  clp.RegisterFloatOption("accuracyPrior",    &params.readAccuracyPrior, "", CommandLineParser::NonNegativeFloat);
  clp.RegisterIntOption("readIndex", &params.readIndex, "", CommandLineParser::NonNegativeInteger);
	clp.RegisterIntListOption("readIndices", &params.readIndices, "", false);
  clp.RegisterIntOption("maxReadIndex", &params.maxReadIndex, "", CommandLineParser::NonNegativeInteger);
	clp.RegisterFlagOption("version", (bool*)&params.printVersion, "");
  clp.RegisterIntOption("substitutionPrior",  &params.substitutionPrior, "", CommandLineParser::NonNegativeInteger);
  clp.RegisterIntOption("deletionPrior",  &params.globalDeletionPrior, "", CommandLineParser::NonNegativeInteger);
  clp.RegisterIntOption("recurseOver", &params.recurseOver, "", CommandLineParser::NonNegativeInteger);
  clp.RegisterStringOption("scoreMatrix", &params.scoreMatrixString, "");
  clp.RegisterFlagOption("printDotPlots", &params.printDotPlots, "");
  clp.RegisterFlagOption("preserveReadTitle", &params.preserveReadTitle,"");
  clp.RegisterFlagOption("forwardOnly", &params.forwardOnly,"");
  clp.RegisterFlagOption("affineAlign", &params.affineAlign, "");
  clp.RegisterIntOption("affineOpen", &params.affineOpen, "", CommandLineParser::NonNegativeInteger);
  clp.RegisterIntOption("affineExtend", &params.affineExtend, "", CommandLineParser::NonNegativeInteger);
	clp.RegisterFlagOption("alignContigs", &params.alignContigs, "", false);
	clp.RegisterStringOption("findex", &params.findex, "", false);
	clp.RegisterIntOption("minInterval", &params.minInterval, "", CommandLineParser::NonNegativeInteger);
	clp.RegisterStringListOption("samqv", &params.samqv, "", false);
	clp.RegisterIntOption("minMapQV", &params.minMapQV, "", CommandLineParser::NonNegativeInteger);
	clp.RegisterIntOption("maxRefine", &params.maxRefine, "", CommandLineParser::NonNegativeInteger);
	clp.RegisterFlagOption("removeContained", &params.removeContainedIntervals, "", false);
	clp.RegisterFlagOption("piecewise", &params.piecewiseMatch, "", false);
	clp.RegisterFlagOption("noSelf", &params.noSelf, "", false);
	clp.RegisterIntOption("maxAnchorGap", &params.maxAnchorGap, "", CommandLineParser::NonNegativeInteger);
	clp.RegisterIntOption("maxGap", &params.maxGap, "", CommandLineParser::NonNegativeInteger);
	clp.RegisterStringOption("fileType", &params.fileType, "");
	clp.RegisterFlagOption("streaming", &params.streaming, "", false);
	clp.RegisterFlagOption("passthrough", &params.passthrough, "", false);	
	clp.ParseCommandLine(argc, argv, params.readsFileNames);
  clp.CommandLineToString(argc, argv, commandLine);
	
  if (params.printVersion) {
    string version;
    GetVersion(version);
    cout << version << endl;
    exit(0);
  }
    

	if (printVerboseHelp) {
		cout << helpString << endl;
		exit(0);
	}

	if (params.printDiscussion) {
		PrintDiscussion();
		exit(0);
	}
	if (argc < 3) {
		cout << conciseHelpString;
		exit(1);
	}
  
  int a, b;
  for (a = 0; a < 5; a++ ){
    for (b = 0; b < 5; b++ ){
      if (a != b) {
        SMRTDistanceMatrix[a][b] += params.mismatch;
      }
      else {
        SMRTDistanceMatrix[a][b] += params.match;
      }
    }
  }
  
  if (params.scoreMatrixString != "") {
    if (StringToScoreMatrix(params.scoreMatrixString, SMRTDistanceMatrix) == false) {
      cout << "ERROR. The string " << endl
           << params.scoreMatrixString << endl
           << "is not a valid format.  It should be a quoted, space separated string of " << endl
           << "integer values.  The matrix: " << endl
           << "    A  C  G  T  N" << endl
           << " A  1  2  3  4  5" << endl
           << " C  6  7  8  9 10" << endl
           << " G 11 12 13 14 15" << endl
           << " T 16 17 18 19 20" << endl
           << " N 21 22 23 24 25" << endl
           << " should be specified as \"1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25\"" << endl;
      exit(1);
    }
  }
  
  FileOfFileNames::ExpandFileNameList(params.readsFileNames);

	// The last reads file is the genome
	if (params.readsFileNames.size() > 0) {
		params.genomeFileName = params.readsFileNames[params.readsFileNames.size()-1];
	}

	if (params.printDiscussion == true) {
		PrintDiscussion();
		exit(0);
	}
	params.MakeSane();
	params.anchorParameters.verbosity = params.verbosity; 

  //
  // The random number generator is used for subsampling for debugging
  // and testing consensus.
  //
  if (params.useRandomSeed == true) {
    InitializeRandomGenerator(params.randomSeed);
  }
  else {
    InitializeRandomGeneratorWithTime();
  }
  
  //
  // Various aspects of timing are stored here.  However this isn't
  // quite finished.
  //
  MappingMetrics metrics;

  ofstream fullMetricsFile;
  if (params.fullMetricsFileName != "") {
    CrucialOpen(params.fullMetricsFileName, fullMetricsFile, std::ios::out);
    metrics.SetStoreList();
  }

	/*
	 * If reading a separate region table, there is a 1-1 correspondence
	 * between region table and bas file.
	 */
	if (params.readSeparateRegionTable) {
		if (FileOfFileNames::IsFOFN(params.regionTableFileName)) {
			FileOfFileNames::FOFNToList(params.regionTableFileName, params.regionTableFileNames);
		}
		else {
			params.regionTableFileNames.push_back(params.regionTableFileName);
		}
	}

  if (params.regionTableFileNames.size() != 0 and 
      params.regionTableFileNames.size() != params.readsFileNames.size() - 1) {
    cout << "Error, there are not the same number of region table files as input files." << endl;
    exit(1);
  }

	//	params.readsFileNames.pop_back();
	if (params.readsFileNames.size() < 2) {
		cout << "Error, you must provide at least one reads file and a genome file." <<endl;
		exit(1);
	}

	//  The input reads files must have file extensions.
	for (int i = 0; i < params.readsFileNames.size()-1; i++) {
		if (params.readsFileNames[i] != "/dev/stdin") {
			size_t dotPos = params.readsFileNames[i].find_last_of('.');
			if (dotPos == string::npos) {
				cout<<"ERROR, the input reads files must include file extensions."<<endl;
				exit(1);
			}
		}
	}

	// -useQuality can not be used in combination with a fasta input
	if (!params.ignoreQualities) {
		for (int i = 0; i < params.readsFileNames.size()-1; i++) {
			size_t dotPos = params.readsFileNames[i].find_last_of('.');
			assert (dotPos != string::npos); //dotPos should have been checked above
			string suffix = params.readsFileNames[i].substr(dotPos+1);
			if (suffix == "fasta") {
				cout<<"ERROR, you can not use -useQuality option when any of the input reads files are in multi-fasta format."<<endl; 
				exit(1);
			}
		}
	}

	if (params.nProc > 1 and params.outFileName == "") {
		cout << "ERROR: You must explicitly specify an output file name with the option -out "<<endl
				 << "when aligning in parallel mode. "<<endl;
		exit(1);
	}

	parentPID = getpid();


	SequenceIndexDatabase<FASTASequence> seqdb;
	SeqBoundaryFtr<FASTASequence> seqBoundary(&seqdb);

	//
	// Initialize the sequence index database if it used. If it is not
	// specified, it is initialized by default when reading a multiFASTA
	// file.
	//
	if (params.useSeqDB) {
		ifstream seqdbin;
		CrucialOpen(params.seqDBName, seqdbin);
		seqdb.ReadDatabase(seqdbin);
	}

	//
	// Make sure the reads file exists and can be opened before
	// trying to read any of the larger data structures.
	//
	

	FASTASequence   fastaGenome;
	T_Sequence      genome;
	FASTAReader     genomeReader;

	// 
	// The genome is in normal FASTA, or condensed (lossy homopolymer->unipolymer) 
	// format.  Both may be read in using a FASTA reader.
	//
	if (!genomeReader.Init(params.genomeFileName)) {
		cout << "Could not open genome file " << params.genomeFileName << endl;
		exit(1);
	}

  if (params.printSAM) {
    genomeReader.computeMD5 = true;
  }
	//
	// If no sequence title database is supplied, initialize one when
	// reading in the reference, and consider a seqdb to be present.
	//
	if (!params.useSeqDB) {
		genomeReader.ReadAllSequencesIntoOne(fastaGenome, &seqdb);
		params.useSeqDB = true;
		if (params.noSelf) {
		  seqdb.BuildNameMaps();
		}
	}
	else {
		genomeReader.ReadAllSequencesIntoOne(fastaGenome);
	}
	genomeReader.Close();
	//
	// The genome may have extra spaces in the fasta name. Get rid of those.
	//
	VectorIndex t;
	for (t = 0; t < fastaGenome.titleLength; t++ ){
		if (fastaGenome.title[t] == ' ') {
			fastaGenome.titleLength = t;
			fastaGenome.title[t] = '\0';
			break;
		}
	}
	genome.seq = fastaGenome.seq;
	genome.length = fastaGenome.length;
	genome.title = fastaGenome.title;
	genome.titleLength = fastaGenome.titleLength;
	genome.ToUpper();


	DNASuffixArray sarray;
	TupleCountTable<T_GenomeSequence, DNATuple> ct;

	int listTupleSize;
	
	ofstream outFile;
  outFile.exceptions(ostream::failbit);
	ofstream unalignedOutFile;
	BWT bwt;

	//
	// Look to see if:
	//  1. no indices have been specified on the command line
	//  2. If so, if a .bwt exitsts
	//     2. If not, if a .sa exists
	
	
	
	if (params.useBwt == false and params.useSuffixArray == false) {
		params.bwtFileName = params.genomeFileName + ".bwt";
		ifstream in(params.bwtFileName.c_str());

		if (in.good()) {
			params.useBwt = true;
		}
		else {
			params.bwtFileName = "";
			params.suffixArrayFileName = params.genomeFileName + ".sa";
			ifstream saIn(params.suffixArrayFileName.c_str());
			if (saIn) {
				params.useSuffixArray = true;
			}
			else {
				params.suffixArrayFileName = "";
			}
		}
		in.close();
	}
	
	
	if (params.useBwt) {

		if (bwt.Read(params.bwtFileName) == 0) {
			cout << "ERROR! Could not read the BWT file. " << params.bwtFileName << endl;
			exit(1);
		}
	}
	else {
		if (!params.useSuffixArray) {
			//
			// There was no explicit specification of a suffix
			// array on the command line, so build it on the fly here.
			//
			genome.ToThreeBit();		
			vector<int> alphabet;
			sarray.InitThreeBitDNAAlphabet(alphabet);
			sarray.LarssonBuildSuffixArray(genome.seq, genome.length, alphabet);
			if (params.minMatchLength > 0) {
				if (params.anchorParameters.useLookupTable == true) {
          if (params.lookupTableLength > params.minMatchLength) {
            params.lookupTableLength = params.minMatchLength;
          }
					sarray.BuildLookupTable(genome.seq, genome.length, params.lookupTableLength);
				}
			}
			genome.ConvertThreeBitToAscii();
			params.useSuffixArray = 1;
    }
		else if (params.useSuffixArray) {
			if (sarray.Read(params.suffixArrayFileName)) {
        if (params.minMatchLength != 0) {
          params.listTupleSize = min(8, params.minMatchLength);
        }
        else {
          params.listTupleSize = sarray.lookupPrefixLength;
        }
        if (params.minMatchLength < sarray.lookupPrefixLength) {
          cerr << "WARNING. The value of -minMatch " << params.minMatchLength << " is less than the smallest searched length of " << sarray.lookupPrefixLength << ".  Setting -minMatch to " << sarray.lookupPrefixLength << "." << endl;
          params.minMatchLength = sarray.lookupPrefixLength;
        }
      }
      else {
        cout << "ERROR. " << params.suffixArrayFileName << " is not a valid suffix array. " << endl
             << " Make sure it is generated with the latest version of sawriter." << endl;
        exit(1);
      }
		}
	}
	
  if (params.minMatchLength < sarray.lookupPrefixLength) {
    cerr << "WARNING. The value of -minMatch " << params.minMatchLength << " is less than the smallest searched length of " << sarray.lookupPrefixLength << ".  Setting -minMatch to " << sarray.lookupPrefixLength << "." << endl;
    params.minMatchLength = sarray.lookupPrefixLength;
  }

	//
	// It is required to have a tuple count table
	// for estimating the background frequencies
	// for word matching. 
	// If one is specified on the command line, simply read
	// it in.  If not, this is operating under the mode 
	// that everything is computed from scratch.
	//
  long l;
  TupleMetrics saLookupTupleMetrics;

	if (params.useCountTable == false) {
		params.countTableName = params.genomeFileName + ".ctab";
		ifstream ctab(params.countTableName.c_str());
		if (ctab.good() == true) {
			params.useCountTable = true;
		}
		else {
			params.countTableName = "";
		}
	}

	if (params.useCountTable) {
		ifstream ctIn;
		CrucialOpen(params.countTableName, ctIn, std::ios::in | std::ios::binary);
		ct.Read(ctIn);
		saLookupTupleMetrics = ct.tm;

	} else {
		saLookupTupleMetrics.Initialize(params.lookupTableLength);
		ct.InitCountTable(saLookupTupleMetrics);
		ct.AddSequenceTupleCountsLR(genome);
	}
	TitleTable titleTable;

	if (params.useTitleTable) {
		ofstream titleTableOut;
		CrucialOpen(params.titleTableName, titleTableOut);
		//
		// When using a sequence index database, the title table is simply copied 
		// from the sequencedb. 
		//
		if (params.useSeqDB) {
			titleTable.Copy(seqdb.names, seqdb.nSeqPos-1);
			titleTable.ResetTableToIntegers(seqdb.names, seqdb.nameLengths, seqdb.nSeqPos-1);
		}
		else {
			//
			// No seqdb, so there is just one sequence. Still the user specified a title
			// table, so just the first sequence in the fasta file should be used. 
			//
			titleTable.Copy(&fastaGenome.title, 1);
			titleTable.ResetTableToIntegers(&genome.title, &genome.titleLength, 1);
			fastaGenome.titleLength = strlen(genome.title);
		}
		titleTable.Write(titleTableOut);
	}
	else {
		if (params.useSeqDB) {
			//
			// When using a sequence index database, but not the titleTable,
			// it is necessary to truncate the titles at the first space to
			// be compatible with the way other alignment programs interpret
			// fasta titles.  When printing the title table, there is all
			// sorts of extra storage space, so the full line is stored.
			//
			seqdb.SequenceTitleLinesToNames();
		}
	}

	ostream  *outFilePtr = &cout;
	ofstream outFileStrm;
	ofstream unalignedFile;
	ostream *unalignedFilePtr = NULL;
	ofstream metricsOut, lcpBoundsOut;
  ofstream anchorFileStrm;
  ofstream clusterOut, *clusterOutPtr;
 
  if (params.anchorFileName != "") {
    CrucialOpen(params.anchorFileName, anchorFileStrm, std::ios::out);
  }

  if (params.clusterFileName != "") {
    CrucialOpen(params.clusterFileName, clusterOut, std::ios::out);
    clusterOutPtr = &clusterOut;
    clusterOut << "total_size p_value n_anchors read_length align_score read_accuracy anchor_probability min_exp_anchors seq_length" << endl;
  }
  else {
    clusterOutPtr = NULL;
  }

	if (params.outFileName != "") {
		CrucialOpen(params.outFileName, outFileStrm, std::ios::out);
		outFilePtr = &outFileStrm;
	}

  if (params.printHeader) {
    switch(params.printFormat) {
    case(SummaryPrint):
      SummaryAlignmentPrinter::PrintHeader(*outFilePtr);
      break;
    case(Interval):
      IntervalAlignmentPrinter::PrintHeader(*outFilePtr);
      break;
    }
  }

	if (params.printUnaligned == true) {
		CrucialOpen(params.unalignedFileName, unalignedFile, std::ios::out);
		unalignedFilePtr = &unalignedFile;
	}
	
	if (params.metricsFileName != "") {
		CrucialOpen(params.metricsFileName, metricsOut);
	}

  if (params.lcpBoundsFileName != "") {
    CrucialOpen(params.lcpBoundsFileName, lcpBoundsOut);
    //    lcpBoundsOut << "pos depth width lnwidth" << endl;
  }
	
	//
	// Configure the mapping database.
	//

	MappingData<T_SuffixArray, T_GenomeSequence, T_Tuple> *mapdb = new MappingData<T_SuffixArray, T_GenomeSequence, T_Tuple>[params.nProc];

	int procIndex;
	pthread_attr_t *threadAttr = new pthread_attr_t[params.nProc];
	//  MappingSemaphores semaphores;
	//
	// When there are multiple processes running along, sometimes there
	// are semaphores to worry about.
	//

	if (params.nProc > 1) {
		semaphores.InitializeAll();
	}
	for (procIndex = 0; procIndex < params.nProc; procIndex++ ){
		pthread_attr_init(&threadAttr[procIndex]);
	}

	//
	// Start the mapping jobs.
	//
	int readsFileIndex = 0;
	if (params.subsample < 1) {
		InitializeRandomGeneratorWithTime();
		reader = new ReaderAgglomerate(params.subsample);
	}
	else {
		reader = new ReaderAgglomerate(params.startRead, params.stride);
	}
  //  In case the input is fasta, make all bases in upper case.
  reader->SetToUpper();


	regionTableReader = new HDFRegionTableReader;
	RegionTable regionTable;
  //
  // Store lists of how long it took to map each read.
  //
  metrics.clocks.SetStoreList(true);
	if (params.useCcs) {
		reader->UseCCS();
	}
	
	vector<int > holeNumbers;
	if (params.readIndex != -1 or params.readIndices.size() > 0) {
		if (params.readIndices.size() > 0) {
			holeNumbers = params.readIndices;
		}
		else {
			holeNumbers.push_back(params.readIndex);
		}
	}
	
  if (params.printSAM) {
    string hdString, sqString, rgString, pgString;
		if (params.passthrough == false) {
			MakeSAMHDString(hdString);
			*outFilePtr << hdString << endl;
		}
		set<string> readGroups;
    for (readsFileIndex = 0; readsFileIndex < params.readsFileNames.size() - 1; readsFileIndex++ ) {    
      reader->SetReadFileName(params.readsFileNames[readsFileIndex], params.fileType);
      reader->Initialize();
			if (params.passthrough == false) {
				string movieNameMD5;
				ScanData scanData;
				reader->GetScanData(scanData);
				MakeMD5(scanData.movieName, movieNameMD5, 10);
				string chipId;
				ParseChipIdFromMovieName(scanData.movieName, chipId);
				//
				// Each movie may only be represented once in the header.
				//
				string changelistId;
				reader->GetChangelistId(changelistId);
				string softwareVersion;
				GetSoftwareVersion(changelistId, softwareVersion);
				
				if (readGroups.find(scanData.movieName) == readGroups.end()) {
					*outFilePtr << "@RG\t" 
											<< "ID:" << movieNameMD5 << "\t" 
											<< "PU:"<< scanData.movieName << "\t" 
											<< "SM:"<< chipId << "\t" 
											<< "PL:PACBIO" << "\t"
											<< "DS:READTYPE=SUBREAD;" 
											<< "CHANGELISTID="<<changelistId <<";"
											<< "BINDINGKIT=" << scanData.bindingKit << ";" 
											<< "SEQUENCINGKIT=" << scanData.sequencingKit << ";"
											<< "FRAMERATEHZ=100;" 
											<< "BASECALLERVERSION=" << softwareVersion;

					int q;
					for (q = 0; q < params.samQVList.nTags; q++){ 
						if (params.samQVList.useqv & params.samQVList.qvFlagIndex[q]) {
							*outFilePtr << ";" << params.samQVList.qvNames[q] << "=" << params.samQVList.qvTags[q];
						}
					}
					*outFilePtr << endl;
					readGroups.insert(scanData.movieName);					
				}
			}
			else {
				//
					// When passing through sam values, do not try and guess the
					// read groups etd., just emit whatever header lines were
					// present.
					//
				int i;
				for (i = 0; i < reader->samReader.allHeaders.size(); i++ ) {
					*outFilePtr << reader->samReader.allHeaders[i] << endl;
				}
			}
			

			if (params.streaming == false) {
				reader->Close();
			}
    }
    seqdb.MakeSAMSQString(sqString);
    *outFilePtr << sqString; // this already outputs endl
		
    string commandLineString;
    clp.CommandLineToString(argc, argv, commandLineString);
    MakeSAMPGString(commandLineString, pgString);
    *outFilePtr << pgString << endl;
  }

	
	for (readsFileIndex = 0; readsFileIndex < params.readsFileNames.size() - 1; readsFileIndex++ ){ 
		params.readsFileIndex = readsFileIndex;

		//
		// Configure the reader to use the correct read and region
		// file names.
		//
		reader->SetReadFileName(params.readsFileNames[params.readsFileIndex], params.fileType);
		if (holeNumbers.size() > 0) {
			reader->InitializeHoleNumbers(holeNumbers);
		}
		//
		// Initialize using already set file names.
		//
		
		int initReturnValue = 1;
		//
		// Not the best approach.  When writing SAM output, all the input
		// files must be scanned for read groups so that the header may be
		// written.  When streaming, we don't want to close this since the
		// streaming pipe will be killed, so during header initialization
		// of streamed files, the reader is not closed.  This is ok since
		// only one file is ever read from (stdin).
		//
		if (params.streaming == false or reader->IsInitialized() == false) {
			 initReturnValue = reader->Initialize();
		} else {
			initReturnValue = reader->IsInitialized();
		}
    string changeListIdString;
    reader->hdfBasReader.GetChangeListID(changeListIdString);
    ChangeListID changeListId(changeListIdString);
    params.qvScaleType = DetermineQVScaleFromChangeListID(changeListId);
		if (reader->FileHasZMWInformation() and params.useRegionTable) {
			if (params.readSeparateRegionTable) {
				if (regionTableReader->Initialize(params.regionTableFileNames[params.readsFileIndex]) == 0) {
					cout << "ERROR! Could not read the region table " << params.regionTableFileNames[params.readsFileIndex] <<endl;
					exit(1);
				}
				params.useRegionTable = true;
			}
			else {
				if (reader->HasRegionTable()) {
					if (regionTableReader->Initialize(params.readsFileNames[params.readsFileIndex]) == 0) {
						cout << "ERROR! Could not read the region table " << params.regionTableFileNames[params.readsFileIndex] <<endl;
						exit(1);
					}
					params.useRegionTable = true;
				}
				else {
					params.useRegionTable = false;
				}
			}
		}
		else {
			params.useRegionTable = false;
		}

		/*
		 * Check to see if there is a region table. If there is a separate
		 * region table, use that (over the region table in the bas
		 * file).  If there is a region table in the bas file, use that,
		 * without having to specify a region table on the command line. 
		 */
		if (params.useRegionTable) {
			regionTable.Reset();
			regionTableReader->ReadTable(regionTable);
			regionTableReader->Close();
			regionTable.SortTableByHoleNumber();
		}

#ifdef USE_GOOGLE_PROFILER
    char *profileFileName = getenv("CPUPROFILE");
    if (profileFileName != NULL) {
      ProfilerStart(profileFileName);
    }
    else {
      ProfilerStart("google_profile.txt");
    }
#endif

		if (initReturnValue > 0) {
			if (params.nProc == 1) {
				mapdb[0].Initialize(&sarray, &genome, &seqdb, &ct, &index, params, reader, &regionTable, 
                            outFilePtr, unalignedFilePtr, &anchorFileStrm, clusterOutPtr);
				mapdb[0].bwtPtr = &bwt;
        if (params.fullMetricsFileName != "") {
          mapdb[0].metrics.SetStoreList(true);
        }
        if (params.lcpBoundsFileName != "") {
          mapdb[0].lcpBoundsOutPtr = &lcpBoundsOut;
        }
        else {
          mapdb[0].lcpBoundsOutPtr = NULL;
        }
				MapReads(&mapdb[0]);
				metrics.Collect(mapdb[0].metrics);
			}
			else {
				pthread_t *threads = new pthread_t[params.nProc];
				for (procIndex = 0; procIndex < params.nProc; procIndex++ ){ 
					//
					// Initialize thread-specific parameters.
					//
            
					mapdb[procIndex].Initialize(&sarray, &genome, &seqdb, &ct, &index, params, reader, &regionTable, 
                                      outFilePtr, unalignedFilePtr, &anchorFileStrm, clusterOutPtr);
					mapdb[procIndex].bwtPtr      = &bwt;
          if (params.fullMetricsFileName != "") {
            mapdb[procIndex].metrics.SetStoreList(true);
          }
          if (params.lcpBoundsFileName != "") {
            mapdb[procIndex].lcpBoundsOutPtr = &lcpBoundsOut;
          }
          else {
            mapdb[procIndex].lcpBoundsOutPtr = NULL;
          }

          if (params.outputByThread) {
            ofstream *outPtr =new ofstream;
            mapdb[procIndex].outFilePtr = outPtr;
            stringstream outNameStream;
            outNameStream << params.outFileName << "." << procIndex;
            mapdb[procIndex].params.outFileName = outNameStream.str();
            CrucialOpen(mapdb[procIndex].params.outFileName, *outPtr, std::ios::out);
          }
					pthread_create(&threads[procIndex], &threadAttr[procIndex], (void* (*)(void*))MapReads, &mapdb[procIndex]);
				}
				for (procIndex = 0; procIndex < params.nProc; procIndex++) {
					pthread_join(threads[procIndex], NULL);
				}
				for (procIndex = 0; procIndex < params.nProc; procIndex++) {
					metrics.Collect(mapdb[procIndex].metrics);
          if (params.outputByThread) {
            delete mapdb[procIndex].outFilePtr;
          }
				}
			}
		}
		reader->Close();
	}
  
  delete reader;

	fastaGenome.Free();
#ifdef USE_GOOGLE_PROFILER
  ProfilerStop();
#endif

	if (mapdb != NULL) {
		delete[] mapdb;
	}
	if (threadAttr != NULL) {
		delete[] threadAttr;
	}
	seqdb.FreeDatabase();
	delete regionTableReader;
	if (params.metricsFileName != "") {
		metrics.PrintSummary(metricsOut);
	}
  if (params.fullMetricsFileName != "") {
    metrics.PrintFullList(fullMetricsFile);
  }
	if (params.outFileName != "") {
		outFileStrm.close();
	}

	cerr << "Finished mapping " << totalReads << " (" << totalBases << " bp)" << endl;


	return 0;
}

