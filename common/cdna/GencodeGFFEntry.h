#ifndef CDNA_GENCODE_GFF_ENTRY
#define CDNA_GENCODE_GFF_ENTRY

#include <string>
using namespace std;

class GencodeGFFEntry {
 public:
  string chr;
  string annotationSource;
  string genomicLocusType;
  unsigned int start, end;
  char scoreNOTUSED;
  char strand;
  char phase;
  string geneId;
  string transcriptId;
  string geneType;
  string geneStatus;
  string geneName;
  string transcriptType;
  string transcriptStatus;
  string transcriptName;
  int level;
  /*
    From: http://www.gencodegenes.org/gencodeformat.html
    gene_id 	ENSGXXXXXXXXXXX *
    transcript_id 	ENSTXXXXXXXXXXX *
    gene_type 	list of biotypes
    gene_status 	{KNOWN, NOVEL, PUTATIVE}
    gene_name 	string
    transcript_type 	list of biotypes
    transcript_status 	{KNOWN, NOVEL, PUTATIVE}
    transcript_name 	string
    level 	1 (verified loci),
    2 (manually annotated loci),
  */

};


/*
 * A GencodeGFFGene is a collection of GencodeGFFEntries, all linked by the same transcript id.
 */



#endif
