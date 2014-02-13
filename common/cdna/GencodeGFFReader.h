#ifndef CDNA_GENCODE_GFF_READER_H_
#define CDNA_GENCODE_GFF_READER_H_

#include "GencodeGFFEntry.h"
#include "utils/StringUtils.h"

template<typename T_Value>
void StringToValue(string valueStr, T_Value &value) {
  stringstream strm(valueStr);
  strm >> value;
}

void StringToValue(string valueStr, string &value) {
  //
  // The string value just has extra quotes around it.
  //
  
  // Make sure the quotes are there.
  if (valueStr.size() < 2) {
    value = "";
    return;
  }
  if (valueStr[0] != '"' or valueStr[valueStr.size()-1] != '"') {
    value = "";
    return;
  }
  
  value = valueStr.substr(1,valueStr.size()-2);
}


void KWStringToPair(string &kwPair, string &key, string &value) {
  stringstream strm(kwPair);
  strm >> key >> value;
}

bool ReadGencodeGFFLine(istream &in, GencodeGFFEntry &entry) {
  // Here is a sample line:
  //  chr1	HAVANA	exon	35721	36073	.	-	.	gene_id "ENSG00000237613.2"; transcript_id "ENST00000461467.1"; gene_type "protein_coding"; gene_status "KNOWN"; gene_name "FAM138A"; transcript_type "processed_transcript"; transcript_status "KNOWN"; transcript_name "FAM138A-002"; level 2; havana_gene "OTTHUMG00000000960.1"; havana_transcript "OTTHUMT00000002843.1";
  string line;
  getline(in, line);
  stringstream strm(line);
  if (! (strm >> entry.chr >> entry.annotationSource >> entry.genomicLocusType >> entry.start >> entry.end
         >> entry.scoreNOTUSED >> entry.strand >> entry.phase) ) {
    return false;
  }

  string kwPairsLine;
  getline(strm, kwPairsLine);
  
  vector<string> kwPairs;
  ParseSeparatedList(kwPairsLine, kwPairs, ';');
  int i;
  for (i = 0; i < kwPairs.size(); i++) {
    string key, value;
    KWStringToPair(kwPairs[i], key, value);
    if (key == "gene_id") {
      StringToValue(value, entry.geneId);
    }
    else if (key == "transcript_id") {
      StringToValue(value, entry.transcriptId);
    }
    else if (key == "gene_type") {
      StringToValue(value, entry.geneType);
    }
    else if (key == "gene_status") {
      StringToValue(value, entry.geneStatus);
    }
    else if (key == "gene_name") {
      StringToValue(value, entry.geneName);
    }
    else if (key == "transcript_type") {
      StringToValue(value, entry.transcriptType);
    }
    else if (key == "transcript_status") {
      StringToValue(value, entry.transcriptStatus);
    }
    else if (key == "transcript_name") {
      StringToValue(value, entry.transcriptName);
    }
    else if (key == "level") {
      StringToValue(value, entry.level);
    }
    else if (key == "havana_gene" or 
             key == "havana_transcript" or 
             key == "tag" or 
             key == "ccdsid" or 
             key == "ont") {
    }
    else {
      cout <<" ERROR.  Keyword-value " << kwPairs[i] << " is not supported." << endl;
      cout << "The offending line is " << line << endl;
      return false;
    }
  }
  return true;
}


#endif
