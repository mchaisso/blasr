#ifndef SDP_FRAGMENT_H_
#define SDP_FRAGMENT_H_


class Fragment {
 public:
	unsigned int x;
	unsigned int y;
	unsigned int weight;
	unsigned int length;
	int index;
	int chainPrev;
	int cost;
	int above;
	unsigned int chainLength;

	unsigned int GetX() const {
		return x;
	}
	unsigned int GetY() const {
		return y;
	}

	Fragment(unsigned int px, unsigned int py, int pweight=0) {
		x = px;
		y = py;
		weight = pweight;
		chainPrev = 0;
    index = length = 0;
    cost = 0;
		above = -1;
	}
	
	int SetAbove(int a) {
		above = a;
		return true;
	}

	int GetAbove(int &a) {
		if (above >= 0) {
			a = above;
			return true;
		}
		else {
			a = -1;
			return false;
		}
	}
	//
	// Provide default constructor that will
	// give bad results if members are not properly initialized
  // later on.
	//
	Fragment() {
		x = -1;
		y = -1;
		chainPrev = 0;
		chainLength = length;
	}
	int LessThanXY(const Fragment &f) const {
		if (x < f.x)
			return 1;
		else if (x == f.x) 
			return y < f.y;
		else 
			return 0;
	}

	int LessThanYX(const Fragment &f) const {
		if (y < f.y)
			return 1;
		else if (y == f.y) 
			return x < f.x;
		else 
			return 0;
	}

	int operator<(const Fragment &f) const {
		// 
		// Sort fragments by diagonal:
		//
		int diag, fDiag;
		diag = (int)(y - x);
		fDiag = (int)f.y - (int)f.x;
		if (diag < fDiag)
			return 1;
		else if (diag == fDiag)
			return (x < f.x);
		else
			return 0;
	}
	Fragment& operator=(const Fragment &rhs) {
		x           = rhs.x;
		y           = rhs.y;
		index       = rhs.index;
		cost        = rhs.cost;
		weight      = rhs.weight;
    length      = rhs.length;
		chainLength = rhs.chainLength;
		chainPrev   = rhs.chainPrev;
		above       = rhs.above;
		return *this;
	}
		
	int operator==(const Fragment &f) const {
		return (x == f.x and y == f.y);
	}
	int operator>(const Fragment &f) const {
		return (!(*this < f) &&  !(*this == f));
	}
	int GetLength() {
		return length;
	}
	int GetXLength() {
		return length;
	}
	int GetYLength() {
		return length;
	}

	void SetLength(int _length) {
		length = _length;
	}
};


#endif
