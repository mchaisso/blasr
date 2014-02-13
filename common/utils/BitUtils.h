#ifndef UTILS_BIT_UTILS_H_
#define UTILS_BIT_UTILS_H_
#include <stdint.h>
#include "../Types.h"
#include <bitset>

int32_t CountBits(uint32_t v) {
	/* 
	 * Attribute Sean Anderson for this method.
	 * The page: http://graphics.stanford.edu/~seander/bithacks.html
	 * contains this code.
	 */
  return bitset<sizeof(uint32_t)*8>(v).count();
/*
	unsigned long long c;
	c =  ((v & 0xfff) * 0x1001001001001ULL & 0x84210842108421ULL) % 0x1f;
	c += (((v & 0xfff000) >> 12) * 0x1001001001001ULL & 0x84210842108421ULL) % 
		0x1f;
	c += ((v >> 24) * 0x1001001001001ULL & 0x84210842108421ULL) % 0x1f;
	return c;*/
}


int GetSetBitPosition64(uint64_t v) {
	unsigned int s;      // Output: Resulting position of bit with rank r [1-64]
	uint64_t a, b, c, d; // Intermediate temporaries for bit count.
	unsigned int t;      // Bit count temporary.
		
	unsigned int r=1;
	// Do a normal parallel bit count for a 64-bit integer,                     
	// but store all intermediate steps.                                        
	// a = (v & 0x5555...) + ((v >> 1) & 0x5555...);
	a =  v - ((v >> 1) & ~0UL/3);
	// b = (a & 0x3333...) + ((a >> 2) & 0x3333...);
	b = (a & ~0UL/5) + ((a >> 2) & ~0UL/5);
	// c = (b & 0x0f0f...) + ((b >> 4) & 0x0f0f...);
	c = (b + (b >> 4)) & ~0UL/0x11;
	// d = (c & 0x00ff...) + ((c >> 8) & 0x00ff...);
	d = (c + (c >> 8)) & ~0UL/0x101;
	t = (d >> 32) + (d >> 48);
	// Now do branchless select!                                                
	s  = 64;
	// if (r > t) {s -= 32; r -= t;}
	s -= ((t - r) & 256) >> 3; r -= (t & ((t - r) >> 8));
	t  = (d >> (s - 16)) & 0xff;
	// if (r > t) {s -= 16; r -= t;}
	s -= ((t - r) & 256) >> 4; r -= (t & ((t - r) >> 8));
	t  = (c >> (s - 8)) & 0xf;
	// if (r > t) {s -= 8; r -= t;}
	s -= ((t - r) & 256) >> 5; r -= (t & ((t - r) >> 8));
	t  = (b >> (s - 4)) & 0xf;
	// if (r > t) {s -= 4; r -= t;}
	s -= ((t - r) & 256) >> 6; r -= (t & ((t - r) >> 8));
	t  = (a >> (s - 2)) & 0x3;
	// if (r > t) {s -= 2; r -= t;}
	s -= ((t - r) & 256) >> 7; r -= (t & ((t - r) >> 8));
	t  = (v >> (s - 1)) & 0x1;
	// if (r > t) s--;
	s -= ((t - r) & 256) >> 8;
	//	s = 64 - s;
	s = s - 1;
	return s;
}

unsigned int GetSetBitPosition32(UInt v) {
	unsigned int r;      // Input: bit's desired rank [1-64].
	unsigned int s;      // Output: Resulting position of bit with rank r [1-64]
	uint32_t a, b, c, d; // Intermediate temporaries for bit count.
	unsigned int t;      // Bit count temporary.
	r = 1;
	// Do a normal parallel bit count for a 64-bit integer,                     
	// but store all intermediate steps.                                        
	// a = (v & 0x5555...) + ((v >> 1) & 0x5555...);
	a =  v - ((v >> 1) & ~0U/3);
	// b = (a & 0x3333...) + ((a >> 2) & 0x3333...);
	b = (a & ~0U/5) + ((a >> 2) & ~0U/5);
	// c = (b & 0x0f0f...) + ((b >> 4) & 0x0f0f...);
	c = (b + (b >> 4)) & ~0U/0x11;
	// d = (c & 0x00ff...) + ((c >> 8) & 0x00ff...);
	d = (c + (c >> 8)) & ~0U/0x101;
	t = (d >> 16) + (d >> 24);


	// Now do branchless select!                                                
	s  = 32;
	// if (r > t) {s -= 16; r -= t;}
	s -= ((t - r) & 256) >> 4; r -= (t & ((t - r) >> 8));
	t  = (c >> (s - 8)) & 0xf;
	// if (r > t) {s -= 8; r -= t;}
	s -= ((t - r) & 256) >> 5; r -= (t & ((t - r) >> 8));
	t  = (b >> (s - 4)) & 0xf;
	// if (r > t) {s -= 4; r -= t;}
	s -= ((t - r) & 256) >> 6; r -= (t & ((t - r) >> 8));
	t  = (a >> (s - 2)) & 0x3;
	// if (r > t) {s -= 2; r -= t;}
	s -= ((t - r) & 256) >> 7; r -= (t & ((t - r) >> 8));
	t  = (v >> (s - 1)) & 0x1;
	// if (r > t) s--;
	s -= ((t - r) & 256) >> 8;
	//	s = 32 - s;
	s = s - 1;
	return s;
}

#endif
