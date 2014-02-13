#ifndef FILES_ALIGNMENT_WRITER_AGGLOMERATE
#define FILES_ALIGNMENT_WRITER_AGGLOMERATE

#include <iostream>
#include <fstream>
#include "../utils.h"

#include "BaseSequenceIO.h"

#include "../FASTAReader.h"
#include "../FASTQReader.h"
#include "../data/hdf/HDFBasWriter.h"
#include "../Enumerations.h"

template <typename T_Sequence>
class WriterAgglomerate :public BaseSequenceIO {
	ofstream     faOut;
 public:
	HDFBasWriter hdfWriter;
	int Initialize(string &pFileName) {
		if (DetermineFileTypeByExtension(pFileName, fileType)) {
			fileName = pFileName;
			return Initialize("","");
		}
	}

	int Initialize(FileType &pFileType, string &pFileName) {
		SetFiles(pFileType, pFileName);
		return Initialize();
	}

	int Initialize(string movieName, string runCode) {
		switch(fileType) {
		case Fasta:
		case Fastq:
			CrucialOpen(fileName, faOut, std::ios::out);
			break;
		case HDFBase:
		case HDFPulse:
			hdfWriter.Initialize(fileName, movieName, runCode);
			break;
		}
		return 1;
	}

	
	int Write(T_Sequence &seq) {
		switch(fileType) {
		case Fastq:
		case Fasta:
				seq.PrintSeq(faOut);
				return 1;
			break;
		case HDFBase:
		case HDFPulse:
				return hdfWriter.Write(seq);
			break;
		}
	}
};

#endif
