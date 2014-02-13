#ifndef TUPLES_TUPLE_MASK
#define TUPLES_TUPLE_MASK
//
// Each f is conveniently 2 bases.
//
//                     
static unsigned long TupleMask[] = {0, 0xf, 0xff, 0xfff, 0xffff, 0xfffff, 0xffffff,  // 0 - 6
                                  0xfffffff, 0xfffffff, 0xffffffff, 0xfffffffff,   // 7 - 9
                                  0xffffffffff, 0xfffffffffff, 0xffffffffffff,     // 10-12
                                  0xfffffffffffff, 0xffffffffffffff, 0xfffffffffffffff, // 13-15
                                  0xffffffffffffffff}; // 16

#endif
