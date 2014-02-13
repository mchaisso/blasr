#ifndef MENAGE_ALIGNER_H_
#define MENAGE_ALIGNER_H_

#include "../../datastructures/alignment/Alignment.h"
#include "../../datastructures/matrix/FlatMatrix.h"



template<typename T_Sequence>
int BellEndAlign(T_Sequence &seqX, T_Sequence &seqY, T_Sequence &seqZ,
								int matchMat[5][5], int gap,
								FlatMatrix3D<int> &scoreMat,
								FlatMatrix3D<Arrow> &pathMat,
								vector<Arrow> &optAlignment) {

	//
	// Three dimensional semi-global alignment.
	// seqY is the only sequence assumed to be an entire sequencing of a template
	// seqX is a suffix of the template, and seqZ is a prefix of the template.
	// 
	
	int x, y, z;
	scoreMat.Grow(seqX.length + 1, seqY.length + 1,seqZ.length + 1);
	pathMat.Grow(seqX.length + 1, seqY.length + 1,seqZ.length + 1);


	/*
	 * When inserting sequences along the x axis, only penalize for 
	 * gapping x against y (x*gap), instead of x against y and z (x*2*gap).
	 */
	for (x = 0; x < seqX.length + 1; x++ ) {
		scoreMat.Set(x,0,0,x*gap);
		pathMat.Set(x,0,0,InsertX);
	}

	/*
	 * Ditto for gapping y.
	 */

	for (y = 0; y < seqY.length + 1; y++ ) {
		scoreMat.Set(0,y,0,y*gap);
		pathMat.Set(0,y,0, InsertY);
	}

	/*
	 * If Z is aligned at the beginning against x and y, then the suffix z overlaps the prefix x,
	 * therefore x and z both fully overlap y, therfore inserting z must gap both x and y.
	 */
	for (z = 0; z < seqZ.length + 1; z++ ){
		scoreMat.Set(0,0,z,z*2*gap);
		pathMat.Set(0,0,z,InsertZ);
	}

	
	int minScore= 0;
	Arrow arrow = NoArrow;
	
	// align y to z, gapping x.
	int yIns, zIns, xIns;
	int yzMatch;
	for (y = 1; y < seqY.length + 1; y++ ){
		for (z = 1; z < seqZ.length + 1; z++ ) {
			yzMatch = scoreMat.Get(0,y-1,z-1) + matchMat[seqY.seq[y-1]][seqZ.seq[z-1]] + gap; // gap in x
			yIns    = scoreMat.Get(0,y-1,z) + gap*2; // penalize for gap in X as well.
			zIns    = scoreMat.Get(0,y,z-1) + gap*2; // ditto
			minScore = MIN(yzMatch, MIN(yIns, zIns));

      arrow = NoArrow;
			if (minScore == yzMatch) {
				arrow  = DiagonalYZ;
			}
			else if (minScore == yIns) {
				arrow  = InsertY;
			}
			else if (minScore == zIns) {
				arrow = InsertZ;
			}
			scoreMat.Set(0,y,z, minScore);
			pathMat.Set(0,y,z, arrow);
		}
	}

	// Align x to z, gapping y.
	int xzMatch;
	for (x = 1; x < seqX.length + 1; x++) {
		for (z = 1; z <seqZ.length + 1; z++) {
			xzMatch = scoreMat.Get(x-1,0,z-1) + matchMat[seqX.seq[x-1]][seqZ.seq[z-1]] + gap; // gap penalty for y
			xIns    = scoreMat.Get(x-1,0,z) + gap*2;
			zIns    = scoreMat.Get(x,0,z-1) + gap*2;
			minScore = MIN(xzMatch, MIN(xIns, zIns));
			if (minScore == xzMatch) {
				arrow = DiagonalXZ;
			}
			else if (minScore == xIns) {
				arrow = InsertX;
			}
			else {
				arrow = InsertZ;
			}
			scoreMat.Set(x,0,z,minScore);
			pathMat.Set(x,0,z,arrow);
		}
	}

	/*
	 * align x an y, gapping z.
	 */
	int xyMatch;
	for (x = 1; x < seqX.length + 1; x++ ) {
		for (y = 1; y < seqY.length + 1; y++) {
			xyMatch = scoreMat.Get(x-1,y-1,0) + matchMat[seqX.seq[x-1]][seqY.seq[y-1]]; // no gap penalty for z;
			xIns = scoreMat.Get(x-1, y, 0) + gap; // gap y, since all of y must be aligned.

			// Only penalize for gapping x to prefixes of y. 
			// 
			yIns = scoreMat.Get(x,y-1,0); 
			if (x < seqX.length) {
				yIns += gap; // The end suffix of y not 
              			 // matching x is not gapped.
			}

			minScore = MIN(xyMatch, MIN(xIns, yIns));
			if (minScore == xyMatch) {
				arrow = DiagonalXY;
			}
			else if (minScore == xIns) {
				arrow = InsertX;
			}
			else {
				arrow = InsertY;
			}
			scoreMat.Set(x,y,0,minScore);
			pathMat.Set(x,y,0,arrow);
		}
	}


	for (x = 1; x < seqX.length + 1; x++) {
		for (y = 1; y < seqY.length + 1; y++ ) {
			for (z = 1; z < seqZ.length + 1; z++ ) {
				//
				// Use sum-of-pairs scoring
				//
				//						 DiagonalXYZ,
				//						 InsertX,InsertY,InsertZ, // imply diagonal yz/xz/xy
				//						 DiagonalXY, DiagonalYZ, DiagonalXZ  // imply insertion of Z/X/Y				
				int xyMatch, xzMatch, yzMatch;
				int xIns, yIns, zIns;

				xyMatch = matchMat[seqX.seq[x-1]][seqY.seq[y-1]];
				xzMatch = matchMat[seqX.seq[x-1]][seqZ.seq[z-1]];
				yzMatch = matchMat[seqY.seq[y-1]][seqZ.seq[z-1]];

				int diagonalXYZ, insertX, insertY, insertZ, diagonalXY, diagonalYZ, diagonalXZ;

				diagonalXYZ = xyMatch + xzMatch + yzMatch + scoreMat.Get(x-1,y-1,z-1);
				diagonalXY = scoreMat.Get(x-1,y-1,z) + xyMatch + gap; // gap on z
				diagonalXZ = scoreMat.Get(x-1,y,z-1) + xzMatch + gap; // gap on y
				
				insertX = scoreMat.Get(x-1,y,z)   + 2*gap; // gap z,y


				/*
				 * The following three conditions cover alignments on the YZ face when x = xmax.
				 * In this plane, the entire x sequence has been aligned, so there should be no
				 * penalty for inserting X here.
				 */

				diagonalYZ = scoreMat.Get(x,y-1,z-1) + yzMatch; // gap on x				
				if (x < seqX.length)
					diagonalYZ += gap;

				insertY = scoreMat.Get(x,y-1,z) + gap;
				if (x < seqX.length) 
					insertY += gap;

				insertZ = scoreMat.Get(x,y,z-1)   + gap; // gap x,y				
				if (x < seqX.length)
					insertZ += gap;

				minScore = MIN(diagonalXYZ,
											 MIN(diagonalXY,
													 MIN(diagonalXZ,
															 MIN(diagonalYZ,
																	 MIN(insertX,
																			 MIN(insertY, insertZ))))));

				scoreMat.Set(x,y,z, minScore);
				if (minScore == diagonalXYZ) {
					pathMat.Set(x,y,z, DiagonalXYZ);
				}
				else if (minScore == diagonalXY){ 
					pathMat.Set(x,y,z,DiagonalXY );
				}
				else if (minScore == diagonalXZ){ 
					pathMat.Set(x,y,z,DiagonalXZ );
				}
				else if (minScore == diagonalYZ ){ 
					pathMat.Set(x,y,z, DiagonalYZ);
				}
				else if (minScore == insertX){
					pathMat.Set(x,y,z, InsertX );
				}
				else if (minScore == insertY){
					pathMat.Set(x,y,z, InsertY);
				}
				else {
					// by default: minScore == insertZ
					pathMat.Set(x,y,z, InsertZ);
				}

			}
		}
	}

	for (x = 1; x < seqX.length + 1; x++) {
		for (y = 1; y < seqY.length + 1; y++ ) {
			for (z = 1; z < seqZ.length + 1; z++ ) {
				assert(pathMat.Get(x,y,z) != Diagonal);
			}
		}
	}
	//
	// Now, find the optimal alignment score allowing for Z sequence to align to an arbitrary prefix of Y.
	//

	x = seqX.length;
	y = seqY.length;
	z = seqZ.length;
	
	// Add gaps after z to make the alignment complete.

	while  (x > 0 or y > 0 or z > 0) {
		arrow = pathMat.Get(x,y,z);
		optAlignment.push_back(arrow);
		if (arrow == DiagonalXYZ) {
			x--; y--; z--;
		}
		else if (arrow == DiagonalXY) {
			x--; y--;
		}
		else if (arrow == DiagonalXZ) {
			x--; z--;
		}
		else if (arrow == DiagonalYZ) {
			y--; z--;
		}
		else if (arrow == InsertX) {
			x--;
		}
		else if (arrow == InsertY) {
			y--;
		}
		else if (arrow == InsertZ) {
			z--;
		}
	}
	if (optAlignment.size() > 1) 
		std::reverse(optAlignment.begin(), optAlignment.end());
	
	return scoreMat.Get(seqX.length, seqY.length, seqZ.length);
}

template<typename T_Sequence>
void CreateBellEndAlignmentStrings(T_Sequence &seqX, T_Sequence &seqY, T_Sequence &seqZ,
																	vector<Arrow> &alignment,
																	string &stringX, string &stringY, string &stringZ) {

	seqX.ToAscii();
	seqY.ToAscii();
	seqZ.ToAscii();
	int a;
	int x, y, z;
	x = 0, y = 0, z = 0;
	for (a = 0; a < alignment.size(); a++) {
		if (alignment[a] == DiagonalXYZ) {
			stringX.push_back(seqX.seq[x]);
			stringY.push_back(seqY.seq[y]);
			stringZ.push_back(seqZ.seq[z]);
			x++, y++, z++;
		}
		else if (alignment[a] == DiagonalXZ) {
			stringX.push_back(seqX.seq[x]);
			stringY.push_back('-');
			stringZ.push_back(seqZ.seq[z]);
			x++, z++;
		}
		else if (alignment[a] == DiagonalXY) {
			stringX.push_back(seqX.seq[x]);
			stringY.push_back(seqY.seq[y]);
			stringZ.push_back('-');
			x++, y++;
		}
		else if (alignment[a] == DiagonalYZ) {
			stringX.push_back('-');
			stringY.push_back(seqY.seq[y]);
			stringZ.push_back(seqZ.seq[z]);
			y++, z++;
		}
		else if (alignment[a] == InsertX) {
			stringX.push_back(seqX.seq[x]);
			stringY.push_back('-');
			stringZ.push_back('-');
			x++;
		}
		else if (alignment[a] == InsertY) {
			stringX.push_back('-');
			stringY.push_back(seqY.seq[y]);
			stringZ.push_back('-');
			y++;
		}
		else if (alignment[a] == InsertZ) {
			stringX.push_back('-');
			stringY.push_back('-');
			stringZ.push_back(seqZ.seq[z]);
			z++;
		}
		else {
			cout << "ERROR! Alignment sequence has an arrow" << alignment[a] << " that is not valid!" << endl;
			assert(0);
		}
	}
}


void PrintMAlignStrings(string &x, string &y, string &z, ostream &out) {

	int lineLength = 70;
	int spos = 0;
	assert(x.length() == y.length());
	assert(y.length() == z.length());
	while (spos < x.length()) {
		int substrLength = lineLength;
		if (spos +lineLength > x.length()) {
			lineLength = x.length() - spos;
		}
		out << x.substr(spos, lineLength) << endl
				<< y.substr(spos, lineLength) << endl
				<< z.substr(spos, lineLength) << endl << endl;
		spos += lineLength;
	}
}

#endif
