#ifndef STICK_ALIGNMENT_PRINTER
#define STICK_ALIGNMENT_PRINTER

#include <string>
#include <iostream>
#include <fstream>
#include "algorithms/alignment/Alignment.h"
#include "algorithms/alignment/AlignmentUtils.h"
template<typename T_Alignment, typename T_QuerySequence, typename T_TargetSequence>
void StickPrintAlignment(T_Alignment &alignment, 
												 T_QuerySequence &query, T_TargetSequence &text, ostream &out,
												 unsigned int qPrintStart = 0,
												 unsigned int tPrintStart = 0,
												 int maxPrintLength = 50) {
	//
	// Eventually hack to add options for alignment format.
	//
	string queryStr, alignStr, textStr;
	DNALength q, t;
	q = alignment.qPos;
	t = alignment.tPos;
	char gapChar = '-';
	//
	// [q/t]PrintStart are the beginning/ending of the subsequence that
	// is being printed however the actual alignment begins at
	// [q/t]PrintStart + alignment->[q/t]Pos, since there may be a small
	// gap between the beginning of the alignment and the beginning of
	// the substring that is aligned.
	//
	out << "         Query: " << alignment.qName << endl;
	out << "        Target: " << alignment.tName << endl;
	out << "         Model: a hybrid of global/local non-affine alignment" <<endl;
	out << "     Raw score: " << alignment.score << endl;
	out << "        Map QV: " << alignment.mapQV << endl;
	out << "  Query strand: " << alignment.qStrand << endl;
	out << " Target strand: " << alignment.tStrand << endl;
	if (alignment.blocks.size() > 0) {
		out << "   QueryRange: " << qPrintStart + alignment.qPos //alignment.blocks[0].qPos 
				<< " -> " 
				<< (qPrintStart + alignment.qPos + 
						+ alignment.blocks[alignment.blocks.size() - 1 ].qPos
						+ alignment.blocks[alignment.blocks.size() - 1 ].length) << " of " << alignment.qLength << endl;
		out << "  TargetRange: " << tPrintStart + alignment.tPos //alignment.blocks[0].tPos 
				<< " -> " 
			//				<< (tPrintStart + alignment.blocks[alignment.blocks.size() - 1].tPos +
				<< (tPrintStart + alignment.tPos +
						alignment.blocks[alignment.blocks.size() - 1 ].tPos + 
						alignment.blocks[alignment.blocks.size() - 1].length) << " of " << alignment.tLength << endl;
	}
	else {
		out << "   QueryRange: NONE " << endl;
		out << "  TargetRange: NONE " << endl;
	}

	if (alignment.blocks.size() == 0) 
		return;
	CreateAlignmentStrings(alignment, query.seq, text.seq, textStr, alignStr, queryStr, query.length, text.length);

	//
	// Process any gaps at the beginning of 
	// the alignment.
	//
	
	DNALength qPos, tPos;

	int coordsPrintWidth = MAX(GetNumberWidth(alignment.qPos + qPrintStart + query.length), 
														 GetNumberWidth(alignment.tPos + tPrintStart + query.length));

	// 
	// Now print the alignment
	VectorIndex lineLength = maxPrintLength;
	if (lineLength > textStr.size())
		lineLength = textStr.size();

	string textLine, queryLine, alignLine;
	DNALength pos = 0;
	tPos = alignment.tPos;
	qPos = alignment.qPos;
	int spacePad = 2;
	string coordPadding(spacePad, ' ');
	string alignStrPadding(coordsPrintWidth + spacePad, ' ');

	while (lineLength > 0) {
		textLine.assign(textStr, pos, lineLength);
		alignLine.assign(alignStr, pos, lineLength);
		queryLine.assign(queryStr, pos, lineLength);
		out.width(coordsPrintWidth);
		out << qPrintStart + qPos;
		out.width(0);
		out << coordPadding << queryLine << endl;
		out << alignStrPadding <<  alignLine << endl;
		out.width(coordsPrintWidth);
		out << tPrintStart + tPos ;
		out.width(0);
		out << coordPadding << textLine  << endl;
		out << endl;
		
		pos += lineLength;
		tPos += lineLength;
		qPos += lineLength;
		//
		// Adjust tPos and qPos for indels.
		//
		VectorIndex p;
		for (p = 0; p < textLine.size(); p++) {
			if (textLine[p] == gapChar)
				tPos--;
		}
		for (p = 0; p < queryLine.size(); p++) {
			if (queryLine[p] == gapChar)
				qPos--;
		}
		
		lineLength = maxPrintLength;
		if (textStr.size() < (pos +lineLength)) {
			// sets lineLength to 0 on the last iteration
			lineLength = textStr.size() - pos;
		}
		else {
			lineLength = maxPrintLength;
		}
	}
}


#endif
