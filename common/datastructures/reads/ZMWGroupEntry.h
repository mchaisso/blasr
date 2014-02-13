#ifndef DATASTRUCTURES_READS_ZMW_GROUP_ENTRY_H_
#define DATASTRUCTURES_READS_ZMW_GROUP_ENTRY_H_

class ZMWGroupEntry {
 public:
	UInt holeNumber;
	UInt x;
	UInt y;
	int  numEvents;
    unsigned char holeStatus;
    ZMWGroupEntry() {
        holeNumber = x = y = 0;
        numEvents = 0;
        holeStatus = '0';
    }
};


#endif
