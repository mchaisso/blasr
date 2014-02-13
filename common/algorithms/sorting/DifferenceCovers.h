#ifndef ALGORITHMS_SORTING_DIFFERENCE_COVERS_H_
#define ALGORITHMS_SORTING_DIFFERENCE_COVERS_H_


#include <vector>
#include "../../Types.h"

#define N_COVERS 5
const UInt diffCovers[N_COVERS][60] = {
	{7,  3,  1, 2, 4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
	{32, 7,  1, 2, 3, 4, 8,12,20, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
	{64, 9,  1, 2, 3, 6,15,17,35,43,60, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
	{111,12, 1, 2, 3, 6,13,28,37,39,45,53,66,94, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
	{2281, 58, 0, 1,2,3,4,5,6,7,8,9,19,38,57,76,95,114,133,152,171,190,229,268,307,346,385,424,463,502,541,580,619,658,697,736,775,814,853,892,931,951,971,991,1011,1031,1051,1071,1091,1111,1131,1132,1133,1134,1135,1136,1137,1138,1139, 1140}};

int InitializeDifferenceCover(int diffCoverSize, UInt &diffCoverLength, UInt *&diffCover) {
	UInt index;
	for (index = 0; index < N_COVERS; index++) {
		if (diffCovers[index][0] == diffCoverSize) {
			diffCoverLength = diffCovers[index][1];
			diffCover = new UInt[diffCoverLength];
			memcpy(diffCover, &diffCovers[index][2], sizeof(UInt)*diffCoverLength);
			return 1;
		}
	}
	return 0;
}

#endif
