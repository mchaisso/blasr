#ifndef DATASTRUCTURES_READS_HOLE_XY_H_
#define DATASTRUCTURES_READS_HOLE_XY_H_

class HoleXY {
 public:
	int16_t xy[2];
	bool operator<(const HoleXY &rhs) const {
		return *this < rhs.xy;
	}
	bool operator<(const int16_t xyP[2]) const {
		if (xy[0] == xyP[0]) {
			return xy[1] < xyP[1];
		}
		else {
			return xy[0] < xyP[0];
		}
	}		
};

#endif
