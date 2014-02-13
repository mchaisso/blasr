#ifndef DEFS_H_
#define DEFS_H_

#ifndef MIN
#define MIN(a, b) ((a) < (b) ? (a) : (b))
#endif


#ifndef MAX
#define MAX(a, b) ((a) > (b) ? (a) : (b))
#endif

#define SWAP(a, b) (((a) == (b)) ? 0 : ((a) ^= (b), (b) ^= (a), (a) ^= (b)))

#ifndef INF_INT
#define INF_INT INT_MAX
#endif

#endif
