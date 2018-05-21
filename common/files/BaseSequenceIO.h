#ifndef FILES_DATA_FILE_TYPE
#define FILES_DATA_FILE_TYPE

#include <string>
#include <iostream>
#include <assert.h>
#include "Enumerations.h"

using namespace std;
class BaseSequenceIO {
 protected:
	FileType fileType;
	string fileName;
 public:

	void SetFiles(FileType &pFileType, string &pFileName) {
		fileType = pFileType;
		fileName = pFileName;
	}
	
	FileType GetFileType() {
		return fileType;
	}

	static int DetermineFileTypeByExtension(string fileName, FileType &type, bool exitOnFailure=true) {
		string::size_type dotPos = fileName.rfind(".");

		if (dotPos != string::npos) {
			string extension;
			extension.assign(fileName, dotPos+1, fileName.size() - (dotPos+1));
			if (extension == "fasta" or
					extension == "fa" or
					extension == "fas" or
					extension == "fsta" or
					extension == "screen"
					) {
				type = Fasta;
				return 1;
			}
			else if (extension == "h5") {
				dotPos = fileName.rfind(".", dotPos-1); 
				extension.assign(fileName, dotPos+1, fileName.size() - (dotPos + 1));
				if (extension == "pls.h5" or extension == "plx.h5" ) {
					type = HDFPulse;
					return 1;
				}
				else if (extension == "bas.h5" or extension == "bax.h5") {
					type = HDFBase;
					return 1;
				}
                else if (extension == "ccs.h5") {
                    type = HDFCCSONLY;
                    return 1;
                }
				else {
					type = None;
					return 0;
				}
			}
			else if (extension == "fastq" or
							 extension == "fq") {
				type = Fastq;
				return 1;
			}
			else if (extension == "4bit" or
							 extension == "fourbit") {
				type = Fourbit;
				assert("Four bit reading is not yet implemented for the reader agglomerate!" == 0);
				return 1;
			}
			else if (extension == "sam") {
				type = SAM_READ;
				return 1;
			}
			
			else {
				type = None;
				if (exitOnFailure) {
					cout << "ERROR, file type '." << extension << "' is not understood to be one of pls.h5, fasta, nor fastq. " << endl;
					exit(1);
				}
				return 0;
			}
			return 0;
		}
    return 0;
	}
};	
#endif
