#include "simulator/ContextOutputList.h"
#include <string>
#include <sstream>
#include <iostream>
#include "utils.h"
#include "DNASequence.h"
#include "FASTAReader.h"
#include "FASTASequence.h"
#include "CommandLineParser.h"
#include "algorithms/metagenomics/FindRandomSequence.h"
#include "utils.h"
#include "statistics/statutils.h"
#include "simulator/LengthHistogram.h"
#include "simulator/OutputSampleListSet.h"
#include "data/hdf/HDFBasWriter.h"
#include "data/hdf/HDFRegionTableWriter.h"
#include "Enumerations.h"
#include "datastructures/metagenome/TitleTable.h"
using namespace std;

long ParseNumber(string number) {
	if (number.size() == 0) {
		return 0;
	}
	
	long n = atoi(number.c_str());

	int l = number.size()-1;
	char s = number[l];
	
	if (s == 'G' or s=='g') {
		n*=1000000000;
	}
	if (s == 'M' or s == 'm') {
		n*=1000000;
	}
	else if (s =='k' or s=='K') {
		n*=1000;
	}
	return n;
}
		

void SetHelp(string & str) {
    stringstream helpStream;
    helpStream 
        << "usage: alchemy outputModel [ options ]" << endl
        << " options: " << endl
        << "  -genome genome.fasta" << endl
        << "            Simulate reads from the reference genome 'genome.fasta'." << endl << endl
        << "  -basesPerFile N" << endl
        << "            Limit the number of bases per output file to this. Allow G/M/K format" << endl << endl
			  << "  -hdf    Output bas.h5 and not fasta." << endl << endl
        << "  -sourceReads filename " << endl
        << "            When set, simulate reads by reading from 'filename', " << endl
        << "            rather than simulating from a genome." << endl 
        << "            The format of the fasta titles should be >read_index|chr|start_pos|end_pos" << endl << endl
        << "  -lengthModel" << endl
        << "            Use lengths from the alchemy model, rather than the read length.  This " << endl
        << "            is used in conjunction with the sourceReadsFile, to modulate the lenghts" << endl
        << "            of the reads." << endl << endl
        << "  -fixedLength length " << endl
        << "            Set simulated read length to a fixed value of 'length', rather than " << endl
        << "            sampling from a length mode." << endl 
        << "  -movieName name (\"simulated_movie\")" << endl 
        << "            Use 'name' for movies rather than m000_000..." << endl << endl
        << "  -titleTable name" <<endl 
        << "            Read in the titleTable to assign chromosome indices from " << endl
        << "            simulated reads." << endl << endl
        << "  -baseFileName name (\"simulated\")" << endl 
        << "            Use an alternative name for the output file, rather than 'simulated'" << endl << endl 
        << "  -nFiles N (1)" << endl
        << "            The number of files to simulate. " << endl << endl 
        << "  -meanLength L(0)" << endl
        << "            When set, scales the length of the average read to L." << endl  << endl 
        << "  -posMap   filename  " << endl
        << "            Use this when running alignment through compareSequences.py " << endl
        << "            and writing to cmp.h5. Specify a map between movie names and " << endl
        << "            chromosome/positions. " << endl << endl 
        << "            When set, the simulated positions are not sored in " << endl
        << "            the bas.h5 files and instead printed to 'filename'" << endl << endl 
        << "  -printPercentRepeat" << endl
        << "            Add to the title table a field that has the percent " << endl
        << "            repeat content of the read shown by lower case in " << endl
        << "            the reference." << endl << endl
				<< "  -ignoreUnplaced When it is impossible to simulate a read in the genome, " << endl
				<< "            skip and move to the next read.  When this is not set, the program will exit." << endl << endl;
		
    str = helpStream.str();
}

int main(int argc, char* argv[]) {
    string refGenomeFileName = "";
    string lengthModelFileName = "";
    string outputModelFileName = "";
    DNALength basesPerFile = 0;
    string sourceReadsFileName = "";
    string titleTableFileName = "";
    int numBasH5Files = 1;
    string basH5BaseFileName = "simulated";
    string movieName = "m000000_000000_00000_cSIMULATED_s";
    bool   doRandGenInit = false;
    bool   usePosMap     = false;
    bool   printPercentRepeat = false;
    string posMapFileName = "";
    vector<string> movieNames;
    bool useLengthModel = false;
    bool useFixedLength = false;
    ofstream posMapFile;
    int scaledLength = 0;
    int fixedLength = 0;
    int nBasFiles = 1;
    bool useLengthsModel = true;
    bool printHelp = false;

		bool printFasta = true;
    CommandLineParser clp;
    string commandLine;
    string helpString;
    SetHelp(helpString);
    vector<string> fns;
		string basesPerFileString;
		bool ignoreUnplaced = false;
    clp.RegisterStringOption("genome", &refGenomeFileName, "");
    clp.RegisterStringOption("basesPerFile", &basesPerFileString, "",
            CommandLineParser::PositiveInteger);
		clp.RegisterFlagOption("hdf", &printFasta, "");
    clp.RegisterStringOption("sourceReads", &sourceReadsFileName, "");
    clp.RegisterStringOption("lengthModel", &lengthModelFileName, "");
    clp.RegisterIntOption("fixedLength", &fixedLength, "",
            CommandLineParser::PositiveInteger);
    clp.RegisterFlagOption("lengthModel", &useLengthModel, "");
    clp.RegisterStringOption("movieName", &movieName, "");
    clp.RegisterStringOption("titleTable", &titleTableFileName, "");
    clp.RegisterStringOption("baseFileName", &basH5BaseFileName, "");
    clp.RegisterIntOption("nFiles", &nBasFiles, "",
            CommandLineParser::PositiveInteger);
    clp.RegisterIntOption("meanLength", &scaledLength, "",
            CommandLineParser::PositiveInteger);
    clp.RegisterStringOption("posMap", &posMapFileName, "");
    clp.RegisterFlagOption("printPercentRepeat", &printPercentRepeat, "");
		clp.RegisterFlagOption("ignoreUnplaced", &ignoreUnplaced, "");
		clp.RegisterFlagOption("initRandomGenerator", &doRandGenInit, "");
    clp.RegisterFlagOption("h", &printHelp, "");

    clp.SetHelp(helpString);
    clp.ParseCommandLine(argc, argv, fns);
    clp.CommandLineToString(argc, argv, commandLine);

    clp.SetProgramName("alchemy");

		basesPerFile = ParseNumber(basesPerFileString);
    outputModelFileName = fns[0];
    if (argc <= 1 or printHelp or outputModelFileName == "") {
        cout << helpString << endl;
        exit(0);
    }

    if (usePosMap) {
			CrucialOpen(posMapFileName, posMapFile, std::ios::out);
    }

    if (sourceReadsFileName == "" and fixedLength == 0) {
        useLengthModel = true;
    }

    if (useLengthModel and fixedLength != 0) {
        cout << "ERROR! You must either use a length model or a fixed length." << endl;
        exit(1);
    }

    if (sourceReadsFileName == "" and basesPerFile == 0) {
        cout << "ERROR! You must specify either a set of read to use as " << endl
             << "original reads for simulation or the total number of bases " << endl
             << "to simulate in each bas.h5 file." << endl;
        exit(1);
    }
 
    if (sourceReadsFileName == "" and refGenomeFileName == "") {
        cout << "ERROR! You must specify a genome to sample reads from or a set of read "<<endl
            << "to use as original reads for simulation." << endl;
        exit(1);
    }

    if (fixedLength != 0 and refGenomeFileName == "") {
        cout << "ERROR! You must specify a genome file if using a fixed length." << endl;
        exit(1);
    }

    if ((fixedLength != 0 or scaledLength != 0) and sourceReadsFileName != "") {
        cout << "ERROR! You cannot specify a fixed length nor mean length with a source " << endl
            << "reads file.  The read lengths are taken from the source reads or the length model." << endl;
        exit(1);
    }

    LengthHistogram   lengthHistogram;
		lengthHistogram.interpolate = true;
    OutputSampleListSet   outputModel(0);
    TitleTable titleTable;

    if (doRandGenInit) {
        InitializeRandomGeneratorWithTime();
    }

    //
    // Read models.
    //
    if (titleTableFileName != "") {
        titleTable.Read(titleTableFileName);
    }


    outputModel.Read(outputModelFileName);

    if (useLengthModel) {
        lengthHistogram.BuildFromAlignmentLengths(outputModel.lengths);
    }


    vector<int> alignmentLengths;
    int meanAlignmentLength;


    if (scaledLength != 0 and useLengthModel) {
        //
        // Scale the histogram so that the average length is 'scaledLength'.
        //

        // 1. Integrate histogram
        long totalLength = 0;
        long totalSamples = 0;
        int hi;
        for (hi = 0; hi < lengthHistogram.lengthHistogram.cdf.size()-1; hi++) {
            int ni;
            ni = lengthHistogram.lengthHistogram.cdf[hi+1] - lengthHistogram.lengthHistogram.cdf[hi];
            totalLength += ni * lengthHistogram.lengthHistogram.data[hi];
        }
        totalSamples = lengthHistogram.lengthHistogram.cdf[lengthHistogram.lengthHistogram.cdf.size()-1];

        float meanSampleLength = totalLength / (1.0*totalSamples);
        float fractionIncrease = scaledLength / meanSampleLength;

        for (hi = 0; hi < lengthHistogram.lengthHistogram.cdf.size(); hi++) {
            lengthHistogram.lengthHistogram.data[hi] *= fractionIncrease;
        }
    }

    FASTAReader inReader, seqReader;
    vector<FASTASequence> reference;
    DNALength refLength = 0;
    int i;
    if (refGenomeFileName != "") {
        inReader.Init(refGenomeFileName);
        inReader.ReadAllSequences(reference);

        for (i = 0; i < reference.size(); i++) {
            refLength += reference[i].length;
        }
    }

    if (sourceReadsFileName !=  "") {
        seqReader.Init(sourceReadsFileName);
    }

    ofstream readsFile;

    //
    // Create and simulate bas.h5 files.
    //
    int baseFileIndex;
    bool readsRemain = true;

    for (baseFileIndex = 0; ((sourceReadsFileName == "" and baseFileIndex < nBasFiles)  // case 1 is reads are generated by file
                or (sourceReadsFileName != "" and readsRemain)); // case 2 is reads are generated by an input file.
            baseFileIndex++) {
        //
        // Prep the base file for writing.
        //
        stringstream fileNameStrm, movieNameStrm;
        movieNameStrm << movieName << "." << baseFileIndex;
        string fullMovieName = movieNameStrm.str();
				ofstream fastaOutFile;
        HDFBasWriter basWriter;
        HDFRegionTableWriter regionWriter;
				RegionTable regionTable;				
				
				if (printFasta == false) {
					fileNameStrm  << fullMovieName <<  ".bas.h5";

					regionTable.CreateDefaultAttributes();


					basWriter.SetPlatform(Springfield);
					//
					// Use a fixed set of fields for now.
					//

					// These are all pulled from the outputModel.
					basWriter.IncludeField("Basecall");
					basWriter.IncludeField("QualityValue");
					basWriter.IncludeField("SubstitutionQV");
					basWriter.IncludeField("SubstitutionTag");
					basWriter.IncludeField("InsertionQV");
					basWriter.IncludeField("DeletionQV");
					basWriter.IncludeField("MergeQV");
					basWriter.IncludeField("DeletionTag");
					basWriter.IncludeField("WidthInFrames");
					basWriter.IncludeField("PreBaseFrames");
					basWriter.IncludeField("PulseIndex");
					

					// Just go from 0 .. hole Number
					basWriter.IncludeField("HoleNumber");
					// Fixed to 0.
					basWriter.IncludeField("HoleXY");
					if (usePosMap == false) {
            basWriter.IncludeField("SimulatedSequenceIndex");
            basWriter.IncludeField("SimulatedCoordinate");
					}
					string changelistId = "2.3.0.0.140640";					
					basWriter.Initialize(fileNameStrm.str(), movieNameStrm.str(), changelistId);
					regionWriter.Initialize(basWriter.pulseDataGroup);					
					//					basWriter.SetChangeListID(changelistId);
					string bindingKit    = "100356300";
					string sequencingKit = "100356200";

					basWriter.AddBindingKit(bindingKit);
					basWriter.AddSequencingKit(sequencingKit);
					
				}
				else {
					fileNameStrm << fullMovieName << ".fasta";
					string fastaFileName = fileNameStrm.str();
					CrucialOpen(fastaFileName, fastaOutFile, std::ios::out);
				}


        //
        // This is mainly used to create the atributes.
        //
        vector<unsigned char> qualityValue, substitutionQV, substitutionTag, insertionQV, deletionQV, mergeQV, deletionTag;
        vector<HalfWord> widthInFrames, preBaseFrames, pulseIndex;

        DNALength numSimulatedBases  = 0;
        FASTASequence sampleSeq;
        //sampleSeq.length = readLength;
        int maxRetry = 10000000;
        int retryNumber = 0;
        int numReads = 0;
        int readLength = 0;

        while (basesPerFile == 0 or numSimulatedBases < basesPerFile) {
            DNALength seqIndex, seqPos;
            if (useLengthModel or fixedLength) {
                if (useLengthModel) {
                    lengthHistogram.GetRandomLength(readLength);
                }
                else {
                    readLength = fixedLength;
                }
            }
            if (refGenomeFileName != "") {
							int readPositionFound;
							readPositionFound = FindRandomPos(reference, seqIndex, seqPos, readLength + (outputModel.keyLength - 1));
							if (readPositionFound == true) {
                sampleSeq.seq    = &reference[seqIndex].seq[seqPos];
                sampleSeq.length = readLength + (outputModel.keyLength - 1);
                assert(reference[seqIndex].length >= sampleSeq.length);
							}
							else {
								if (ignoreUnplaced) {
									continue;
								}
							}
            }
            else if (sourceReadsFileName != "") {
                if (seqReader.GetNext(sampleSeq) == false) {
                    readsRemain = false;
                    break;
                }
                if (sampleSeq.length < outputModel.keyLength) {
                    continue;
                }
                //
                // Now attempt to parse the position from the fasta title.
                //

                if (useLengthModel) {
                    int tryNumber = 0;
                    readLength = 0;
                    int maxNTries = 1000;
                    int tryBuffer[5] = {-1,-1,-1,-1,-1};
                    while (tryNumber < maxNTries and readLength < outputModel.keyLength) {
                        lengthHistogram.GetRandomLength(readLength);
                        readLength = sampleSeq.length = min(sampleSeq.length, (unsigned int) readLength);
                        tryBuffer[tryNumber%5] = readLength;
                        tryNumber++;
                    }
                    if (tryNumber >= maxNTries) {
                        cout << "ERROR. Could not generate a read length greater than the " << outputModel.keyLength << " requried " <<endl
                            << "minimum number of bases using the length model specified in the alchemy." <<endl
                            << "model.  Something is either wrong with the model or the context length is too large." <<endl;
                        cout << "The last few tries were: " << tryBuffer[0] << " " << tryBuffer[1] << " " << tryBuffer[2] << " " << tryBuffer[3] << " " << tryBuffer[4] << endl;
                        exit(1);
                    }
                }

                readLength = sampleSeq.length;
                vector<string> tokens;
                Tokenize(sampleSeq.title, "|", tokens);
                if (tokens.size() == 4) {
                    seqPos = atoi(tokens[2].c_str());
                    if (titleTableFileName == "") {
                        seqIndex = 0;
                    }
                    else {
                        int index;
                        titleTable.Lookup(tokens[1], index);
                        seqIndex = index;
                    }
                }
                else {
                    seqPos   = 0;
                }
            }

            numSimulatedBases += readLength;

            int p;
            // create the sample sequence
            int contextLength = outputModel.keyLength;
            int contextMiddle = contextLength / 2;
            string outputString;

            int nDel = 0;
            int nIns = 0;

            //
            // Simulate to beyond the sample length.
            //
            qualityValue.clear(); 
            substitutionQV.clear(); 
            substitutionTag.clear(); 
            insertionQV.clear(); 
            deletionQV.clear(); 
            deletionTag.clear();
            pulseIndex.clear();
            widthInFrames.clear();
            preBaseFrames.clear();
            assert(sampleSeq.length > contextMiddle + 1);

            for (p = contextMiddle;
                    p < sampleSeq.length - contextMiddle - 1; p++) {
                string refContext;
                refContext.assign((const char*) &sampleSeq.seq[p-contextMiddle], contextLength);

                string outputContext;
                int    contextWasFound;
                OutputSample sample;
                int i;
                for (i = 0; i < refContext.size(); i++) { refContext[i] = toupper(refContext[i]);}
                outputModel.SampleRandomSample(refContext, sample);
                if (sample.type == OutputSample::Deletion ) {
                    //
                    // There was a deletion.  Advance in reference, then output
                    // the base after the deletion.
                    //
                    p += sample.nNuc;
                    ++nDel;
                }

								else {
									int cp;
									//
									// Add the sampled context, possibly multiple characters because of an insertion.
									//
									for (i = 0; i < sample.nucleotides.size(); i++) {
                    outputString.push_back(sample.nucleotides[i]);
                    qualityValue.push_back(sample.qualities[i].qv[0]);
                    deletionQV.push_back(sample.qualities[i].qv[1]);
                    insertionQV.push_back(sample.qualities[i].qv[2]);
                    substitutionQV.push_back(sample.qualities[i].qv[3]);
                    mergeQV.push_back(sample.qualities[i].qv[4]);
                    deletionTag.push_back(sample.qualities[i].tags[0]);
                    substitutionTag.push_back(sample.qualities[i].tags[1]);
                    pulseIndex.push_back(sample.qualities[i].frameValues[0]);
                    preBaseFrames.push_back(sample.qualities[i].frameValues[1]);
                    widthInFrames.push_back(sample.qualities[i].frameValues[2]);
									}
									nIns += sample.qualities.size() - 1;
								}
            }
            if (outputString.find('N') != outputString.npos or
                    outputString.find('n') != outputString.npos) {
                cerr << "WARNING!  The sampled string " << endl << outputString << endl
                    << "should not contain N's, but it seems to.  This is being ignored "<<endl
                    << "for now so that simulation may continue, but this shouldn't happen"<<endl
                    << "and is really a bug." << endl;
                numSimulatedBases -= readLength;
                continue;
            }
            //
            // Ok, done creating the read, now time to create some quality values!!!!!
            //
            SMRTSequence read;
            read.length = outputString.size();
            read.Allocate(read.length);
            memcpy(read.seq, outputString.c_str(), read.length * sizeof(unsigned char));
            assert(qualityValue.size() == read.length * sizeof(unsigned char));
            memcpy(read.qual.data, &qualityValue[0], read.length * sizeof(unsigned char));
            memcpy(read.deletionQV.data, &deletionQV[0], read.length * sizeof(unsigned char));
            memcpy(read.insertionQV.data, &insertionQV[0], read.length * sizeof(unsigned char));
            memcpy(read.substitutionQV.data, &substitutionQV[0], read.length * sizeof(unsigned char));
            memcpy(read.mergeQV.data, &mergeQV[0], read.length * sizeof(unsigned char));						
            memcpy(read.deletionTag, &deletionTag[0], read.length * sizeof(unsigned char));
            memcpy(read.substitutionTag, &substitutionTag[0], read.length * sizeof(unsigned char));
            memcpy(read.pulseIndex, &pulseIndex[0], read.length * sizeof(int));
            memcpy(read.preBaseFrames, &preBaseFrames[0], read.length * sizeof(HalfWord));
            memcpy(read.widthInFrames, &widthInFrames[0], read.length * sizeof(HalfWord));

            //
            // The pulse index for now is just fake data.
            //
            int i;
            for (i = 0; i < read.length; i++) {
                read.pulseIndex[i] = 1;
            }
            read.xy[0] = seqIndex;
            read.xy[1] = seqPos;
            read.zmwData.holeNumber = numReads;

						if (printFasta == false) {
							basWriter.Write(read);
							
							// Record where this was simulated from.
							if (usePosMap == false) {
                basWriter.WriteSimulatedCoordinate(seqPos);
                basWriter.WriteSimulatedSequenceIndex(seqIndex);
							}
							else {
                posMapFile << fullMovieName << "/" << numReads << "/0_" << read.length << " " << seqIndex << " "<< seqPos;
                if (printPercentRepeat) {
									DNALength nRepeat = sampleSeq.GetRepeatContent();
									posMapFile << " " << nRepeat*1.0/sampleSeq.length;
                }
                posMapFile << endl;
							}
							
						}
						else {
							string title;
							stringstream titleStream;
							titleStream << fullMovieName << "/" << numReads;
							title = titleStream.str();
							read.CopyTitle(title);
							read.PrintSeq(fastaOutFile);
						}
				
            RegionAnnotation region;
            region.row[0] = read.zmwData.holeNumber;
            region.row[1] = 1;
            region.row[2] = 0;
            region.row[3] = read.length;
            region.row[4] = 1000; // Should be enough.
            regionWriter.Write(region);
            region.row[1] = 2; // Rewrite for hq region encompassing everything.
            regionWriter.Write(region);      
            if (sourceReadsFileName != "") {
                sampleSeq.Free();
            }
            read.Free();
            ++numReads;
        }

				if (printFasta == false) {
					regionWriter.Finalize(regionTable.columnNames, 
																regionTable.regionTypes,
																regionTable.regionDescriptions,
																regionTable.regionSources);
					basWriter.Close();
					numReads = 0;
					//
					// The bas writer should automatically flush on closing.
					//
				}
    }
    if (usePosMap) {
        posMapFile.close();
    }

    for (i = 0; i < reference.size(); i++) {
        reference[i].Free();
    }
}

