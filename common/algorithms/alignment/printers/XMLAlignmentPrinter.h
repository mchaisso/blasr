#ifndef XML_ALIGNMENT_PRINTER_H_
#define XML_ALIGNMENT_PRINTER_H_

#include "../../../formatting/SimpleXMLUtils.h"

template<typename T_Alignment, typename T_Sequence> 
void	CompareXMLPrintAlignment(T_Alignment &alignment, 
													 T_Sequence &query, T_Sequence &text, ostream &out,
													 int qPrintStart = 0,
													 int tPrintStart = 0,
													 int maxPrintLength = 50) {
		/*
		 * Sample alignment:
		 *
		 <hit name="x15_y33_1220208-0008_m081205_152444_Uni_p2_b15" unalignedLength="1301" start="1051" end="1016" strand="-" targetStart="1" targetEnd="44" targetStrand="+">
		 <zScore value="-6.091"/>
		 <nInsert value="1" percent="2.86" />
		 <nDelete value="9" percent="25.71" />
		 <nMismatch value="1" percent="2.86" />
		 <nCorrect value="24" percent="68.57" />
		 <alignment><query>
		 AG--CGTTCC-TATGG-TG-GGGTCGTTA-ACT---GTCGCCAG
		 </query><target>
		 AGCCCG-TCCTTATGGTTGAGGGTTGTTACACTTCGGTCGCCAG
		 </target></alignment>
		 </hit>
		*/
		char strand[2] = {'+', '-'};
		string tAlignStr, alignStr, qAlignStr;
		CreateAlignmentStrings(alignment, query.seq, text.seq, tAlignStr, alignStr, qAlignStr);
    int alignLength = tAlignStr.size();
    if (alignLength == 0) {
      alignLength = 1; // Make sure there are no divide by zero.
      alignment.nIns = 0;
      alignment.nDel = 0;
      alignment.nMismatch = 0;
      alignment.nMatch = 0;
    }
    int lastBlock = alignment.blocks.size()-1;
		out << BeginDataEntry(string("hit"),
													CreateKeywordValuePair(string("name"), alignment.qName) +
													CreateKeywordValuePair(string("unalignedLength"), alignment.qLength) +
													CreateKeywordValuePair(string("start"), alignment.qAlignedSeqPos + alignment.blocks[0].qPos) + 
													CreateKeywordValuePair(string("end"),  alignment.qAlignedSeqPos + alignment.blocks[lastBlock].qPos + alignment.blocks[lastBlock].length) +
													CreateKeywordValuePair(string("strand"), strand[alignment.qStrand]) + 
													CreateKeywordValuePair(string("targetStart"), alignment.tAlignedSeqPos + alignment.blocks[0].tPos ) +													 
													CreateKeywordValuePair(string("targetEnd"), alignment.tAlignedSeqPos + alignment.blocks[lastBlock].tPos  + alignment.blocks[lastBlock].length) + 
													CreateKeywordValuePair(string("targetStrand"), strand[alignment.tStrand])) << endl;
		out << CreateDataEntry(string("zScore"),
													 CreateKeywordValuePair(string("value"), alignment.zScore)) << endl;
		out << CreateDataEntry(string("nInsert"),
													 CreateKeywordValuePair(string("value"), alignment.nIns) + " " +
													 CreateKeywordValuePair(string("percent"), alignment.nIns*0.5/alignLength)) 
				<< endl;
		out << CreateDataEntry(string("nDelete"),
													 CreateKeywordValuePair(string("value"), alignment.nDel) + " " + 
													 CreateKeywordValuePair(string("percent"), alignment.nDel*0.5/alignLength)) 
				<< endl;
		out << CreateDataEntry(string("nMismatch"),
													 CreateKeywordValuePair(string("value"), alignment.nMismatch) + " " + 
													 CreateKeywordValuePair(string("percent"), alignment.nMismatch*0.5/alignLength)) 
				<< endl;
		out << CreateDataEntry(string("nCorrect"),
													 CreateKeywordValuePair(string("value"), alignment.nMatch) + " " +
													 CreateKeywordValuePair(string("percent"), alignment.nMatch*0.5/alignLength)) 
				<< endl;
		
		out << CreateStartEntry(string("alignment"), string("")) << CreateStartEntry(string("query"), string("")) << endl;
		out << qAlignStr << endl;
		out << CreateEndEntry(string("query")) << endl;
		out << CreateStartEntry(string("target"), string("")) << endl;
		out << tAlignStr << endl;
		out << CreateEndEntry(string("target")) << endl;
		out << CreateEndEntry(string("alignment")) << endl;
		out << CreateEndEntry(string("hit")) << endl;
	};




#endif
