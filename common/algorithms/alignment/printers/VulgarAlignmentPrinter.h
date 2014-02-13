#ifndef VULGAR_ALIGNMENT_PRINTER_H_
#define VULGAR_ALIGNMENT_PRINTER_H_

#include <string>
#include <sstream>

template <typename T_Alignment>
int CreateVulgarString(T_Alignment &alignment, string &vstring) {
	stringstream vstream;
	VectorIndex b;
	int tGap, qGap, cGap;
	if (alignment.blocks.size() == 0) {
		vstring = "";
		return 1;
	}
	for (b = 0; b < alignment.blocks.size() - 1; b++) {
		tGap = (alignment.blocks[b+1].tPos 
						- (alignment.blocks[b].length + 
							 alignment.blocks[b].tPos));
		qGap = (alignment.blocks[b+1].qPos 
						- (alignment.blocks[b].length + 
							 alignment.blocks[b].qPos));
		if (tGap > 0 and qGap > 0)
			cGap = abs(tGap - qGap);
		else 
			cGap = 0;
		tGap -= cGap;
		qGap -= cGap;
		vstream << " M " << alignment.blocks[b].length + cGap;
		if (tGap > 0) {
			vstream << " D " << tGap;
		}
		else {
			vstream << " I " << qGap;
		}
	}
	if (alignment.blocks.size() > 0) {
		vstream << " M " << alignment.blocks[alignment.blocks.size() - 1].length;
	}
	vstring = vstream.str();
	return 1;
}

template <typename T_Alignment>
void VulgarPrintAlignment(T_Alignment &alignment, ostream &out) {
	string vstring;
	CreateVulgarString(alignment, vstring);
	out << vstring;
}

#endif
