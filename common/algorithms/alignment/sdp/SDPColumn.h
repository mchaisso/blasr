#ifndef SDP_COLUMN_H_
#define SDP_COLUMN_H_

class SDPColumn {
 public:
	int col;
	int optFragment;
	int operator<(const SDPColumn & rhs) const {
		return col < rhs.col;
	}
	int operator<(const int y) const {
		return col < y;
	}
	int operator==(const SDPColumn &rhs) const {
		return col == rhs.col;
	}
	int operator>(const SDPColumn &rhs) const {
		return (!(*this < rhs) && !(*this == rhs));
	}
};

#endif
