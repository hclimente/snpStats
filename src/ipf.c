/* IPF for 2^K contingency table (K<=14) */

#include <stdlib.h>
#include <math.h>
#include <R.h>
#include <Rinternals.h>
#include "ipf.h"

/* 
   Count bits set in an integer (up to 14 bit integer)
   See http://graphics.stanford.edu/~seander/bithacks.html
*/

inline int bitcount(unsigned int x) {
  return((x * 0x200040008001ULL & 0x111111111111111ULL) % 0xf);
}

/* Extract specified bits (mask) and shift to end of word */

unsigned int bitxtr(unsigned int x, unsigned int mask) {
  unsigned int res = 0;
  unsigned int add = 01;
  while (x) {
    if (mask & 0x01) {
      if (x & 0x01) 
	res = res|add;
      add = add<<1;
    }
    x = x>>1;
    mask = mask>>1;
  }
  return(res);
}
  
	
/* 
   K                number of dimensions 
   observed         input table, first subscript varying fastest
   nterms           number of terms (fitted margins)
   terms            array of terms; each term is represented by a bit pattern 
                    in which each bit indicates whether the corresponding 
                    factor appears in the term or not
   expected         output table of fitted freqencies
                    if expected[0]>=0, this should also hold an initial fit
   maxit            maximum number of steps
   eps              convergence criterion
*/

int ipf(int K, const double *observed, 
	const int nterms, const unsigned int *terms,  double *expected,
	const int maxit, const double eps) {

  /* Size of observed and expected arrays */
  int size = (1<<K);
  /* Initialize */
  
  if (expected[0]<0.0) 
    for (int i=0; i<size; i++)
      expected[i] = 1.0;

  /* Calculate work space needed and allocate */

  int maxsize = 0;
  for (int j=0; j<nterms; j++) {
    unsigned int mask = terms[j];
    int mdim = bitcount(mask);
    int msize = 1<<mdim;
    if (msize>maxsize)
      maxsize = msize;
  }
  double *mexp = (double *) Calloc(maxsize, double);
  double *mobs = (double *) Calloc(maxsize, double);

  /* IPF */
  
  int it = 0; /* step counter */
  double test = 0.0; /* convergence test */
  while (it<maxit) {
    /* Loop over terms */
    for (int j=0; j<nterms; j++) {
      unsigned int mask = terms[j];
      int mdim = bitcount(mask);
      int msize = 1<<mdim;
      int mbytes = msize*sizeof(double);
      memset(mexp, 0, mbytes);
      memset(mobs, 0, mbytes);
     /* Calculate observed and expected margins for current fit */
      for (unsigned int i=0; i<size; i++) {
	int mi = bitxtr(i, mask);
	mobs[mi] += observed[i];
	mexp[mi] += expected[i];
      }
      /* Scaling factors and convergence test */
      for (int i=0; i<msize; i++) {
	double mexpi = mexp[i];
	if (mexpi) {
	  double scale = mobs[i]/mexpi;
	  double ti = fabs(scale - 1.0);
	  if (ti>test)
	  test = ti;
	  mexp[i] = scale;
	}
      }
      /* Scale fitted distribution */
      for (unsigned int i=0; i<size; i++) {
	int mi = bitxtr(i, mask);
	expected[i] *= mexp[mi];
      }
    }
    if (test<eps) {
      Free(mobs);
      Free(mexp);
      return(0); /* Convergence */
    }
    it++;
  }
  Free(mobs);
  Free(mexp);
  return(1); /* Maximum steps reached */
}
