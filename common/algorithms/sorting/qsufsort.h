#ifndef QSUFSORT_H_
#define QSUFSORT_H_
#include <assert.h>

void suffixsort(int *x, int *p, int n, int k, int l);


/* qsufsort.c
   Copyright 1999, N. Jesper Larsson, all rights reserved.

   This file contains an implementation of the algorithm presented in "Faster
   Suffix Sorting" by N. Jesper Larsson (jesper@cs.lth.se) and Kunihiko
   Sadakane (sada@is.s.u-tokyo.ac.jp).

   This software may be used freely for any purpose. However, when distributed,
   the original source must be clearly stated, and, when the source code is
   distributed, the copyright notice must be retained and any alterations in
   the code must be clearly marked. No warranty is given regarding the quality
   of this software.*/

#include <limits.h>

#define KEY(p)          (V[*(p)+(h)])


#define MED3(a, b, c)   (KEY(a)<KEY(b) ?                                \
                         (KEY(b)<KEY(c) ? (b) : KEY(a)<KEY(c) ? (c) : (a)) \
                         : (KEY(b)>KEY(c) ? (b) : KEY(a)>KEY(c) ? (c) : (a)))

/* Subroutine for select_sort_split and sort_split. Sets group numbers for a
   group whose lowest position in I is pl and highest position is pm.*/

#include <vector>
#include <algorithm>

using namespace std;

template<typename T_Index, long T_Index_MAX=0>
  class LarssonSuffixSort {
 private:
  T_Index *I;
  T_Index *V;
  T_Index r;
  T_Index h;
  T_Index SWAPPtr(T_Index *p, T_Index *q) { T_Index tmp=*(p); *(p)=*(q); *(q)=tmp; return *(q);}
	vector<char> boundaries;

 public:
	T_Index INDEX_MAX;
	LarssonSuffixSort() {
		INDEX_MAX = T_Index_MAX;
        r = h = 0;
        I = V = NULL;
	}
  void operator()(T_Index *x, T_Index *p, T_Index n, T_Index k, T_Index l)  {
    suffixsort(x,p,n,k,l);
  }
  
  void update_group(T_Index *pl, T_Index *pm)
  {
    int g;
    
    g=pm-I;                      /* group number.*/
    V[*pl]=g;                    /* update group number of first position.*/
    if (pl==pm) {
			/*MC*/
			assert(pl -I >= 0);
			boundaries[pl - I] = 1;
			//      *pl=-1;                   /* one element, sorted group.*/
		}
    else
      do                        /* more than one element, unsorted group.*/
        V[*++pl]=g;            /* update group numbers.*/
      while (pl<pm);
  }
  
  /* Quadratic sorting method to use for small subarrays. To be able to update
     group numbers consistently, a variant of selection sorting is used.*/

  void select_sort_split(T_Index *p, T_Index n) {
    T_Index *pa, *pb, *pi, *pn;
    T_Index f, v;

    pa=p;                        /* pa is start of group being picked out.*/
    pn=p+n-1;                    /* pn is last position of subarray.*/
    while (pa<pn) {
      for (pi=pb=pa+1, f=KEY(pa); pi<=pn; ++pi)
        if ((v=KEY(pi))<f) {
          f=v;                /* f is smallest key found.*/
          SWAPPtr(pi, pa);       /* place smallest element at beginning.*/
          pb=pa+1;            /* pb is position for elements equal to f.*/
        } else if (v==f) {     /* if equal to smallest key.*/
          SWAPPtr(pi, pb);       /* place next to other smallest elements.*/
          ++pb;
        }
      update_group(pa, pb-1);   /* update group values for new group.*/
      pa=pb;                    /* continue sorting rest of the subarray.*/
    }
    if (pa==pn) {                /* check if last part is single element.*/
			assert(pa - I >= 0);
			//			assert(boundaries[pa-I] == 0 );
      V[*pa]=pa-I;
			/*MC*/
			//      *pa=-1;                   /* sorted group.*/
			boundaries[pa - I] = 1; 
    }
  }

  /* Subroutine for sort_split, algorithm by Bentley & McIlroy.*/

  T_Index choose_pivot(T_Index *p, T_Index n) {
    T_Index *pl, *pm, *pn;
    T_Index s;
   
    pm=p+(n>>1);                 /* small arrays, middle element.*/
    if (n>7) {
      pl=p;
      pn=p+n-1;
      if (n>40) {               /* big arrays, pseudomedian of 9.*/
        s=n>>3;
        pl=MED3(pl, pl+s, pl+s+s);
        pm=MED3(pm-s, pm, pm+s);
        pn=MED3(pn-s-s, pn-s, pn);
      }
      pm=MED3(pl, pm, pn);      /* midsize arrays, median of 3.*/
    }
    return KEY(pm);
  }

  /* Sorting routine called for each unsorted group. Sorts the array of integers
     (suffix numbers) of length n starting at p. The algorithm is a ternary-split
     quicksort taken from Bentley & McIlroy, "Engineering a Sort Function",
     Software -- Practice and Experience 23(11), 1249-1265 (November 1993). This
     function is based on Program 7.*/

  void sort_split(T_Index *p, T_Index n)
  {
    T_Index *pa, *pb, *pc, *pd, *pl, *pm, *pn;
    T_Index f, v, s, t;

    if (n<7) {                   /* multi-selection sort smallest arrays.*/
      select_sort_split(p, n);
      return;
    }

    v=choose_pivot(p, n);
    pa=pb=p;
    pc=pd=p+n-1;
    while (1) {                  /* split-end partition.*/
      while (pb<=pc && (f=KEY(pb))<=v) {
        if (f==v) {
          SWAPPtr(pa, pb);
          ++pa;
        }
        ++pb;
      }
      while (pc>=pb && (f=KEY(pc))>=v) {
        if (f==v) {
          SWAPPtr(pc, pd);
          --pd;
        }
        --pc;
      }
      if (pb>pc)
        break;
      SWAPPtr(pb, pc);
      ++pb;
      --pc;
    }
    pn=p+n;
    if ((s=pa-p)>(t=pb-pa))
      s=t;
    for (pl=p, pm=pb-s; s; --s, ++pl, ++pm)
      SWAPPtr(pl, pm);
    if ((s=pd-pc)>(t=pn-pd-1))
      s=t;
    for (pl=pb, pm=pn-s; s; --s, ++pl, ++pm)
      SWAPPtr(pl, pm);

    s=pb-pa;
    t=pd-pc;
    if (s>0)
      sort_split(p, s);
    update_group(p+s, p+n-t-1);
    if (t>0)
      sort_split(p+n-t, t);
  }

  /* Bucketsort for first iteration.

     Input: x[0...n-1] holds integers in the range 1...k-1, all of which appear
     at least once. x[n] is 0. (This is the corresponding output of transform.) k
     must be at most n+1. p is array of size n+1 whose contents are disregarded.

     Output: x is V and p is I after the initial sorting stage of the refined
     suffix sorting algorithm.
  */
      
  void bucketsort(T_Index *x, T_Index *p, T_Index n, T_Index k)
  {
    T_Index *pi, i, c, d, g;
  
    for (pi=p; pi<p+k; ++pi) {
			//      *pi=-1;                   /* mark linked lists empty.*/
			assert(pi - p >= 0);
			assert(pi - p == pi - I);
			//			boundaries[pi-p] = 0;
		}
		int *buckets = new int[k];
		T_Index *starts  = new T_Index[k];
		/*MC+1*/
		for (i = 0; i < k; i++ ){
			buckets[i] = -1;
		}
		/*MC-1*/
		for (i=0; i<=n; ++i) {
			/*MC+2*/
			if (buckets[x[i]] == -1) {
				starts[x[i]] = i;
			}
      x[i]=buckets[c=x[i]];           /* insert in linked list.*/
      buckets[c]=i;
			//*x[i]=pi[c=x[i]]
			p[c] = i;
			/*MC-2*/
    }
    for (pi=p+k-1, i=n; pi>=p; --pi) {
      d=x[c=*pi];               /* c is position, d is next in list.*/
      x[c]=g=i;                 /* last position equals group number.*/
			//			if (d>=0) {               /* if more than one element in group.*/
			if (c != starts[pi - p]) {
        p[i--]=c;              /* p is permutation for the sorted x.*/
        do {
          d=x[c=d];           /* next in linked list.*/
          x[c]=g;             /* group number in x.*/
          p[i--]=c;           /* permutation in p.*/
        } 
				while (c != starts[pi - p]);
				//while (d>=0);
      } else {
				/*MC*/
				boundaries[i--]=true;
				//        p[i--]=-1;             /* one element, sorted group.*/
			}
    }
		delete[] starts;
		delete[] buckets;
  }

  /* Transforms the alphabet of x by attempting to aggregate several symbols into
     one, while preserving the suffix order of x. The alphabet may also be
     compacted, so that x on output comprises all integers of the new alphabet
     with no skipped numbers.

     Input: x is an array of size n+1 whose first n elements are positive
     integers in the range l...k-1. p is array of size n+1, used for temporary
     storage. q controls aggregation and compaction by defining the maximum value
     for any symbol during transformation: q must be at least k-l; if q<=n,
     compaction is guaranteed; if k-l>n, compaction is never done; if q is
     INT_MAX, the maximum number of symbols are aggregated into one.
   
     Output: Returns an integer j in the range 1...q representing the size of the
     new alphabet. If j<=n+1, the alphabet is compacted. The global variable r is
     set to the number of old symbols grouped into one. Only x[n] is 0.*/

  T_Index transform(T_Index *x, T_Index *p, T_Index n, T_Index k, T_Index l, T_Index q)
  {
    T_Index b, c, d, e, i, j, m, s;
    T_Index *pi, *pj;
   
    for (s=0, i=k-l; i; i>>=1)
      ++s;                      /* s is number of bits in old symbol.*/
    e=INDEX_MAX>>s;                /* e is for overflow checking.*/
    for (b=d=r=0; r<n && d<=e && (c=d<<s|(k-l))<=q; ++r) {
      b=b<<s|(x[r]-l+1);        /* b is start of x in chunk alphabet.*/
      d=c;                      /* d is max symbol in chunk alphabet.*/
    }
    m=(1<<(r-1)*s)-1;            /* m masks off top old symbol from chunk.*/
    x[n]=l-1;                    /* emulate zero terminator.*/
    if (d<=n) {                  /* if bucketing possible, compact alphabet.*/
      for (pi=p; pi<=p+d; ++pi)
        *pi=0;                 /* zero transformation table.*/
      for (pi=x+r, c=b; pi<=x+n; ++pi) {
        p[c]=1;                /* mark used chunk symbol.*/
        c=(c&m)<<s|(*pi-l+1);  /* shift in next old symbol in chunk.*/
      }
      for (i=1; i<r; ++i) {     /* handle last r-1 positions.*/
        p[c]=1;                /* mark used chunk symbol.*/
        c=(c&m)<<s;            /* shift in next old symbol in chunk.*/
      }
      for (pi=p, j=1; pi<=p+d; ++pi)
        if (*pi)
          *pi=j++;            /* j is new alphabet size.*/
      for (pi=x, pj=x+r, c=b; pj<=x+n; ++pi, ++pj) {
        *pi=p[c];              /* transform to new alphabet.*/
        c=(c&m)<<s|(*pj-l+1);  /* shift in next old symbol in chunk.*/
      }
      while (pi<x+n) {          /* handle last r-1 positions.*/
        *pi++=p[c];            /* transform to new alphabet.*/
        c=(c&m)<<s;            /* shift right-end zero in chunk.*/
      }
    } else {                     /* bucketing not possible, don't compact.*/
      for (pi=x, pj=x+r, c=b; pj<=x+n; ++pi, ++pj) {
        *pi=c;                 /* transform to new alphabet.*/
        c=(c&m)<<s|(*pj-l+1);  /* shift in next old symbol in chunk.*/
      }
      while (pi<x+n) {          /* handle last r-1 positions.*/
        *pi++=c;               /* transform to new alphabet.*/
        c=(c&m)<<s;            /* shift right-end zero in chunk.*/
      }
      j=d+1;                    /* new alphabet size.*/
    }
    x[n]=0;                      /* end-of-string symbol is zero.*/
    return j;                    /* return new alphabet size.*/
  }

  /* Makes suffix array p of x. x becomes inverse of p. p and x are both of size
     n+1. Contents of x[0...n-1] are integers in the range l...k-1. Original
     contents of x[n] is disregarded, the n-th symbol being regarded as
     end-of-string smaller than all other symbols.*/

  void suffixsort(T_Index *x, T_Index *p, T_Index n, T_Index k, T_Index l)
  {
    T_Index *pi, *pk;
    T_Index i, j, s, sl;
		boundaries.resize(n+1);
		fill(boundaries.begin(), boundaries.end(), 0);

    V=x;                         /* set global values.*/
    I=p;
   
    if (n>=k-l) {                /* if bucketing possible,*/
      j=transform(V, I, n, k, l, n);
      bucketsort(V, I, n, j);   /* bucketsort on first r positions.*/
    } else {
      transform(V, I, n, k, l, INDEX_MAX);
      for (i=0; i<=n; ++i)
        I[i]=i;                /* initialize I with suffix numbers.*/
      h=0;
      sort_split(I, n+1);       /* quicksort on first r positions.*/
    }
    h=r;                         /* number of symbols aggregated by transform.*/
   
    while (*I<=n) {
      pi=I;                     /* pi is first position of group.*/
      sl=0;                     /* sl is negated length of sorted groups.*/
      do {
				/*MC-3*/
				//        if ((s=*pi)<0) {
				s = *pi;
				if (boundaries[pi - I] == 1) {
					//          pi-=s;              /* skip over sorted group.*/
					s = 1;
					pi+=s;
          sl+=s;              /* add negated length to sl.*/
        } else {
          if (sl) {
						//            *(pi+sl)=sl;     /* combine sorted groups before pi.*/
						assert(boundaries[(pi-sl) - I] == 1);
						assert(pi - sl >= I);
            *(pi-sl)=sl;     /* combine sorted groups before pi.*/
            sl=0;
          }
          pk=I+V[s]+1;        /* pk-1 is last position of unsorted group.*/
          sort_split(pi, pk-pi);
          pi=pk;              /* next group.*/
        }
      } while (pi<=I+n);
      if (sl)                   /* if the array ends with a sorted group.*/
        *(pi-sl)=sl;           /* combine sorted groups at end of I.*/
      h=2*h;                    /* double sorted-depth.*/
    }

    for (i=0; i<=n; ++i)         /* reconstruct suffix array from inverse.*/
      I[V[i]]=i;
  }
};
#endif
