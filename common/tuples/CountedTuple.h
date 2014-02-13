#ifndef COUNTED_TUPLE_H_
#define COUNTED_TUPLE_H_

#include "tuples/DNATuple.h"

class CountedDNATuple : public DNATuple {
 public:
	int count;
	int IncrementCount() {
		count = count + 1;
	}
	
};


#endif
