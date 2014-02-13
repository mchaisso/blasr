#include <cmath>
#include <math.h>
#include <cstdlib>

/******************************************************************************/
/* randn()
code yanked from http://www.dreamincode.net/code/snippet1446.htm
 * 
 * Normally (Gaussian) distributed random numbers, using the Box-Muller 
 * transformation.  This transformation takes two uniformly distributed deviates
 * within the unit circle, and transforms them into two independently 
 * distributed normal deviates.  Utilizes the internal rand() function; this can
 * easily be changed to use a better and faster RNG.
 * 
 * The parameters passed to the function are the mean and standard deviation of 
 * the desired distribution.  The default values used, when no arguments are 
 * passed, are 0 and 1 - the standard normal distribution.
 * 
 * 
 * Two functions are provided:
 * 
 * The first uses the so-called polar version of the B-M transformation, using
 * multiple calls to a uniform RNG to ensure the initial deviates are within the
 * unit circle.  This avoids making any costly trigonometric function calls.
 * 
 * The second makes only a single set of calls to the RNG, and calculates a 
 * position within the unit circle with two trigonometric function calls.
 * 
 * The polar version is generally superior in terms of speed; however, on some
 * systems, the optimization of the math libraries may result in better 
 * performance of the second.  Try it out on the target system to see which
 * works best for you.  On my test machine (Athlon 3800+), the non-trig version
 * runs at about 3x10^6 calls/s; while the trig version runs at about
 * 1.8x10^6 calls/s (-O2 optimization).
 * 
 * 
 * Example calls:
 * randn_notrig();	//returns normal deviate with mean=0.0, std. deviation=1.0
 * randn_notrig(5.2,3.0);	//returns deviate with mean=5.2, std. deviation=3.0
 * 
 * 
 * Dependencies - requires <cmath> for the sqrt(), sin(), and cos() calls, and a
 * #defined value for PI.
 */

/******************************************************************************/
//	"Polar" version without trigonometric calls
double randn_notrig(double mu=0.0, double sigma=1.0) {
	static bool deviateAvailable=false;	//	flag
	static float storedDeviate;			//	deviate from previous calculation
	double polar, rsquared, var1, var2;
	
	//	If no deviate has been stored, the polar Box-Muller transformation is 
	//	performed, producing two independent normally-distributed random
	//	deviates.  One is stored for the next round, and one is returned.
	if (!deviateAvailable) {
		
		//	choose pairs of uniformly distributed deviates, discarding those 
		//	that don't fall within the unit circle
		do {
			var1=2.0*( double(rand())/double(RAND_MAX) ) - 1.0;
			var2=2.0*( double(rand())/double(RAND_MAX) ) - 1.0;
			rsquared=var1*var1+var2*var2;
		} while ( rsquared>=1.0 || rsquared == 0.0);
		
		//	calculate polar tranformation for each deviate
		polar=sqrt(-2.0*log(rsquared)/rsquared);
		
		//	store first deviate and set flag
		storedDeviate=var1*polar;
		deviateAvailable=true;
		
		//	return second deviate
		return var2*polar*sigma + mu;
	}
	
	//	If a deviate is available from a previous call to this function, it is
	//	returned, and the flag is set to false.
	else {
		deviateAvailable=false;
		return storedDeviate*sigma + mu;
	}
}


/******************************************************************************/
//	Standard version with trigonometric calls
#define PI 3.14159265358979323846

double randn_trig(double mu=0.0, double sigma=1.0) {
	static bool deviateAvailable=false;	//	flag
	static float storedDeviate;			//	deviate from previous calculation
	double dist, angle;
	
	//	If no deviate has been stored, the standard Box-Muller transformation is 
	//	performed, producing two independent normally-distributed random
	//	deviates.  One is stored for the next round, and one is returned.
	if (!deviateAvailable) {
		
		//	choose a pair of uniformly distributed deviates, one for the
		//	distance and one for the angle, and perform transformations
		dist=sqrt( -2.0 * log(double(rand()) / double(RAND_MAX)) );
		angle=2.0 * PI * (double(rand()) / double(RAND_MAX));
		
		//	calculate and store first deviate and set flag
		storedDeviate=dist*cos(angle);
		deviateAvailable=true;
		
		//	calcaulate return second deviate
		return dist * sin(angle) * sigma + mu;
	}
	
	//	If a deviate is available from a previous call to this function, it is
	//	returned, and the flag is set to false.
	else {
		deviateAvailable=false;
		return storedDeviate*sigma + mu;
	}
}

#include "FASTAReader.h"
#include "FASTASequence.h"

#include "statistics/statutils.h"
#include "utils.h"
#include "algorithms/metagenomics/FindRandomSequence.h"
#include "CommandLineParser.h"

#include <string>

int main(int argc, char* argv[]) {
  string genomeFileName;
  int nBasesToSimulate = 0;

  CommandLineParser clp;

  int fragmentLengthMean = 0;
  int fragmentLengthSD = 0;
  int numFragmentsMean = 0;
  int numFragmentsSD = 0;

  string adapterFileName = "";
  string outputFileName = "";
  clp.SetProgramSummary("sageify - sample a bunch of sequences from a genome, and concatenate them.");

  //
  // Required files.
  //
  clp.RegisterStringOption("genome", &genomeFileName, "The reference genome (multi FASTA)");
  clp.RegisterIntOption("nBases",    &nBasesToSimulate, "The number of bases to simulate.", 
                        CommandLineParser::PositiveInteger);
  clp.RegisterIntOption("fragLen",   &fragmentLengthMean, "The average length of a fragment.", 
                        CommandLineParser::PositiveInteger, false);
  clp.RegisterIntOption("fragNum", &numFragmentsMean, "The number of fragments that are ligated on average.",
                        CommandLineParser::PositiveInteger, false);
  clp.RegisterStringOption("out", &outputFileName, "Output file name. Two files will be created, output, and output.pos"
                           "where output.pos contains the original locations of all of the sampled reads.", true);
  clp.RegisterIntOption("fragSD",    &fragmentLengthSD, "The standard deviation of fragment length.", 
                        CommandLineParser::PositiveInteger, false);
  clp.RegisterIntOption("fragNumSD", &numFragmentsSD, "The standard deviation of number of fragments (sqrt(fragNum))", 
                        CommandLineParser::NonNegativeInteger);
  clp.RegisterStringOption("adapter", &adapterFileName, "Possible adapter sequence between fragments.", false);

  clp.ParseCommandLine(argc, argv);

  ofstream outputFile;

  //
  // Check some of the input.
  //
  if (fragmentLengthMean == 0) {
    cout << "The fragment length must be greater than 0." << endl;
    exit(1);
  }
  if (numFragmentsMean == 0) {
    cout << "The mean number of fragments must be greater than 0." << endl;
    exit(1);
  }

  // Make sure writing to the output file works.
  CrucialOpen(outputFileName, outputFile, std::ios::out);

  //
  // Read in the input genome to sample from.
  //
  vector<FASTASequence> genome;
  FASTAReader reader;
  reader.Initialize(genomeFileName);
  reader.ReadAllSequences(genome);
  
  int nSampledBases = 0;

  FASTASequence adapter;
  if (adapterFileName != "") {
    reader.Close();
    reader.Initialize(adapterFileName);
    reader.GetNext(adapter);
  }
  
  if (numFragmentsSD == 0) {
    numFragmentsSD = sqrt(numFragmentsMean);
  }

  if (fragmentLengthSD == 0) {
    fragmentLengthSD = sqrt(fragmentLengthMean);
  }
  int seqIndex = 0;
  
  ofstream fragmentPosFile;
  string fragmentPosFileName = outputFileName + ".pos";
  CrucialOpen(fragmentPosFileName, fragmentPosFile, std::ios::out);
  fragmentPosFile << "ligated fragment chromIdx chromTitle start length" << endl;

  while (nSampledBases < nBasesToSimulate) {
    unsigned int seqIndex;
    unsigned int seqPos;
    
    int nFragments = randn_trig(numFragmentsMean, numFragmentsSD);

    int f;
    FASTASequence ligatedSequence;
    stringstream titleStrm;
    titleStrm << seqIndex;
    if (adapter.length > 0) {
      ligatedSequence.Append(adapter);
    }

    for (f = 0; f < nFragments; f++) {
      int fragLen = randn_trig(fragmentLengthMean, fragmentLengthSD);
      FindRandomPos(genome, seqIndex, seqPos, fragLen);
      DNASequence fragment;
      fragment.seq = &genome[seqIndex].seq[seqPos];
      fragment.length = fragLen;
      ligatedSequence.Append(fragment);
      titleStrm << " " << genome[seqIndex].title << "_" << seqPos << "_" << fragLen;
      fragmentPosFile << seqIndex << " " << f << " " << seqIndex << " " << genome[seqIndex].title << " " << seqPos << " " << fragLen << endl;
      if (adapter.length > 0) {
        ligatedSequence.Append(adapter);
      }
    }
    ligatedSequence.CopyTitle(titleStrm.str());
    ligatedSequence.PrintSeq(outputFile);
    nSampledBases += ligatedSequence.length;
  }
  outputFile.close();
}
