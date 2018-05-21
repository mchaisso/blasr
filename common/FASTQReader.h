#ifndef FASTQ_READER_H_
#define FASTQ_READER_H_
#include "FASTAReader.h"
#include "FASTQSequence.h"
#include "qvs/QualityValue.h"
#include <stdio.h>
#include <limits.h>

class FASTQReader : public FASTAReader {
 public:

    FASTQReader() : FASTAReader() {
        endOfReadDelim = '+';
    }

    long GetNext(FASTASequence &seq) {
			return ((FASTAReader*)this)->GetNext(seq, '@');
    }

    unsigned char phredQVtoPacbioQV(unsigned char phredQV){
        int qual = floor(100.0 * log10(pow(10.0, phredQV/10.0) - 1.0) + 0.5); 
        qual = qual > 250 ? 250 : qual;
        qual = qual < 1   ? 1   : qual;
        return (unsigned char) qual;
    }

    int GetNext(FASTQSequence &seq) {
        char c;
				if (file.eof() == true) {
					return false;
				}
        AdvanceToTitleStart('@');
				if (file.eof() == true) {
					return false;
				}
        CheckValidTitleStart('@');
        ReadTitle(seq.title, seq.titleLength);
        // Title ends on '\n', consume that;
				int res;
				res=FASTAReader::GetNext(seq, '@');
				if (res == 0) {
					return false;
				}
				
				string line;
				// move past '+'
				getline(file, line);
				if (line[0] != '+') {
					cout << "ERROR in fastq format. Expcting '+', got " << line[0] << endl;
					exit(0);
				}
        if (seq.length > 0) {
          seq.qual.Allocate(seq.length);
					char buf[seq.length+1];
					// using only seq.length stops 1 base early for \0.
					file.get((char*) buf, seq.length+1);
					memcpy(seq.qual.data, buf, seq.length);
					file.get();// get the \n
        }
        else {
          seq.qual.data = NULL;
        }
				seq.deleteOnExit = true;
				return true;
    }



    int Advance(int nSteps) {
			while (nSteps > 0 and file.eof() == false) {
				int i;
				for (i = 0; i < 4; i++) {
					string line;
					getline(file,line);
				}
			}
    }
};


#endif
