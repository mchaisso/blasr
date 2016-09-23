#ifndef INTERVAL_SEARCH_PARAMETERS_H_
#define INTERVAL_SEARCH_PARAMETERS_H_

class IntervalSearchParameters {
 public:
	int   globalChainType;
	float maxPValue;
	bool  overlap;
	int   minMatch;
	int   minInterval;
	int   maxAnchorGap;
	bool  noSelf;
	IntervalSearchParameters() {
		globalChainType = 0;
		maxPValue       = log(0.1);
		overlap         = false;
		minMatch        = 0;
		minInterval     = 0;
		maxAnchorGap    = 0;
		noSelf          = false;
	}
};


#endif
