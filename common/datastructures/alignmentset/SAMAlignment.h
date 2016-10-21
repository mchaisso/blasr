#ifndef ALIGNMENT_SAM_ALIGNMENT_H_
#define ALIGNMENT_SAM_ALIGNMENT_H_

#include <string>
#include <sstream>

#include "SAMKeywordValuePair.h"
#include "datastructures/alignment/CigarString.h"
#include "datastructures/alignmentset/SAMSupplementalQVList.h"
#include "utils/StringUtils.h"

class SAMAlignment {
 public:
  enum SAMAlignmentRequiredFields {
    S_QNAME,
    S_FLAG,
    S_RNAME,
    S_POS,
    S_MAPQV,
    S_CIGAR,
    S_RNEXT,
    S_PNEXT,
    S_TLEN,
    S_SEQ,
    S_QUAL};

  static const char* SAMAlignmentRequiredFieldNames[];

  string qName;
  unsigned int flag;
  string rName;
  unsigned int pos;
  short mapQV;
  CigarString cigar;
  string rNext;
  int pNext;
  int tLen;
	int zmw;
	int numPasses;
	float readQuality;
	int subreadContext;
  string seq;
  string qual;

    // Optional tags defined in blasr:
    // "RG" read group Id
    // "AS" alignment score
    // "XS" read alignment start position without counting previous soft clips (1 based) 
    // "XE" read alignment end position without counting previous soft clips (1 based) 
    // "XL" aligned read length 
    // "XQ" query read length
    // "XT" # of continues reads, always 1 for blasr 
    // "NM" # of subreads 
    // "FI" read alignment start position (1 based) 
    //
  float score;
  int   xs, xe;
	//
	// qs and qe are the positions in the original sequence, though the
	// read record should not contain the entire orginal
	// "polymerase-read".
	int qStart, qEnd;
  int xl;
  int xq;

  string rg;
  int   as;
  int   xt;
  int   nm;
  int   fi;
  string optTagStr; 

	//
	// Quality values.
	//
  string qi, qd, qs, qm, ts, td;	

	float snr[4];
  SAMAlignment() {
    //
    // Initialize all optional fields.  Required fields will be
    // assigned a value later.
    //
    score = xs = xe = as = xt = xq = nm = fi = xl = 0;
    rg = optTagStr = "";
		readQuality = 0;
		numPasses = 0;
		snr[0] = snr[1] = snr[2] = snr[3] = 0;
  }

  void PrintSAMAlignment(ostream & out = cout) {
      out << qName << "\t" << flag  << "\t" << rName << "\t"
          << pos   << "\t" << mapQV << "\t" << cigar << "\t"
          << rNext << "\t" << pNext << "\t" << tLen  << "\t"
          << seq   << "\t" << qual  << "\t" << optTagStr 
          << "\n";
  }

  // Find position of the nth character in a string.
  int FindPosOfNthChar(string str, int n, char c) {
    if (n < 1) {
      cout << "Nth should be a positive integer." << endl;
      exit(0);
    }
    int count = 1;
    int pos = str.find(c, 0);
    // pos is the position of the first character c;
    while(count < n and pos != string::npos) {
        pos = str.find(c, pos+1);
        count ++;
    }
    return pos;
  }

  // Trim the end '\n\r' characters from a string.
  string TrimStringEnd(string str) {
    string newStr = str;
    while(newStr[newStr.size()-1] == '\r' or
          newStr[newStr.size()-1] == '\n') {
        newStr.erase(newStr.size()-1, 1);
    }
    return newStr;
  }

	bool StoreSNR(istream &strm) {
		strm.get();
		strm.get();
		if (!(strm >> snr[0])) return false;
		strm.get();
		if (!(strm >> snr[1])) return false;
		strm.get();
		if (!(strm >> snr[2])) return false;
		strm.get();
		if (!(strm >> snr[3])) return false;

		return true;
	}
	
  bool StoreValues(string &line,  int lineNumber=0, bool allowUnaligned=false) {
    stringstream strm(line);
    vector<bool> usedFields;
    usedFields.resize(S_QUAL);
    fill(usedFields.begin(), usedFields.end(), false);
    string kvPair;
    int i;
    bool parseError = false;
    SAMAlignmentRequiredFields field;
    //
    // Define a temporary mapqv value that gets over a GMAP bug that prints a mapqv < 0.
    //
    int tmpMapQV;
    if (!(strm >> qName)) {
      parseError = true;
      field = S_QNAME;
    }
    else if (! (strm >> flag) ){ 
      parseError = true;
      field = S_FLAG;
    }
    else if (! (strm >> rName) ) {
      parseError = true;
      field = S_RNAME;
    }
    else if (! (strm >> pos) ) {
      parseError = true;
      field = S_POS;
    }
    else if (! (strm >> tmpMapQV)) {
      parseError = true; field = S_MAPQV;
    }
    else if (! (strm >> cigar)) {
      parseError = true; field = S_CIGAR;
    }

    else if (! (strm >> rNext)) {
      parseError = true; field = S_RNEXT;
    }

    else if (! (strm >> pNext)) {
      parseError = true; field = S_PNEXT;
    }
    else if (! (strm >> tLen)) {
      parseError = true; field = S_TLEN;
    }
    else if (! (strm >> seq)) {
      parseError = true; field = S_SEQ;
    }
    else if (! (strm >> qual)) {
      parseError = true; field = S_QUAL;
    }

    mapQV = (unsigned char) tmpMapQV;

    // Find posisition of the 11th tab.
    int optTagsStartPos = FindPosOfNthChar(strm.str(), 11, '\t');
    // Save all optional tags.
    if (optTagsStartPos != string::npos) {
        optTagStr = strm.str().substr(optTagsStartPos+1);
        optTagStr = TrimStringEnd(optTagStr);
    } 

    //
    // If not aligned, stop trying to read in elements from the sam string.
    //
    if (allowUnaligned == false and rName == "*") {
      return true;
    }

    if (parseError) {
      cout << "Error parsing alignment line " << lineNumber << ". Missing or error in field " << SAMAlignmentRequiredFieldNames[field] << endl;
			cout << line << endl;
      exit(1);
    }
    
    //
    // Now parse optional data.
    //
    while (strm) {
      string kvName, kvType, kvValue;
      string typedKVPair;
      if ((strm >> typedKVPair) == 0) {
        break;
      }
      if (TypedKeywordValuePair::Separate(typedKVPair, kvName, kvType, kvValue)) {
        stringstream strm(kvValue);
        if (kvName == "RG") {
          rg = kvValue;
        }
        else if (kvName == "AS") {
          strm >> as;
        }
        else if (kvName == "XS") {
          strm >> xs;
        }
        else if (kvName == "XE") {
          strm >> xe;
        }
        else if (kvName == "XL") {
          strm >> xl;
        }
        else if (kvName == "XT") {
          strm >> xt;
        }
        else if (kvName == "NM") {
          strm >> nm;
        }
        else if (kvName == "FI") {
          strm >> fi;
        }
        else if (kvName == "XQ") {
          strm >> xq;
        }
				else if (kvName == "zm") {
					strm >> zmw;
				}
				else if (kvName == "qs") {
					strm >> qStart;
				}
				else if (kvName == "qe") {
					strm >> qEnd;
				}
				else if (kvName == "np") {
					strm >> numPasses;
				}
				else if (kvName == "rq") {
					strm >> readQuality;
				}
				else if (kvName == "cx") {
					strm >> subreadContext;
				}
				else if (kvName == "sn") {
					if (StoreSNR(strm) == false) {
						cout << "ERROR. Could not store SNR " << typedKVPair << endl;;
					}
				}
				
				
				else if (kvName == SupplementalQVList::qvTags[SupplementalQVList::I_Insertion-1]) {
					strm >> qi;
				}
				else if (kvName == SupplementalQVList::qvTags[SupplementalQVList::I_Deletion-1]) {
					strm >> qd;
				}
				else if (kvName == SupplementalQVList::qvTags[SupplementalQVList::I_Substitution-1]) {
			    strm >> qs;
				}
				else if (kvName == SupplementalQVList::qvTags[SupplementalQVList::I_Merge-1]) {
					strm >> qm;
				}
				else if (kvName == SupplementalQVList::qvTags[SupplementalQVList::I_SubstitutionTag-1]) {
					strm >> ts;
				}
				else if (kvName == SupplementalQVList::qvTags[SupplementalQVList::I_DeletionTag-1]) {
					strm >> td;
				}
      }
      else {
        cout << "ERROR.  Could not parse typed keyword value " << typedKVPair << endl;
      }
    }
		if (zmw == -1) {
			//
			// Parse the zmw from the qname,since it was not specified on the command line
			//
			vector<string> qNameTokens;
			Tokenize(qName, "/", qNameTokens);
			if (qNameTokens.size() == 3) {
				zmw = atoi(qNameTokens[1].c_str());
			}
		}
	}
};

const char* SAMAlignment::SAMAlignmentRequiredFieldNames[] = { "QNAME", "FLAG", 
                                                               "RNAME", "POS", 
                                                               "MAPQV", "CIGAR", 
                                                               "RNEXT", "PNEXT", 
                                                               "TLEN", "SEQ", "QUAL"} ;


using namespace std;

class SAMPosAlignment : public SAMAlignment {
 public:
  unsigned int qStart, qEnd;
  unsigned int tStart, tEnd;
  int qStrand, tStrand;
};

/* 
 * Write the full one later
 */




  

#endif
