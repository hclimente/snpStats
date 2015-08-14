/* Modified for R_xlen_t 26/6/2015 */

#include <R.h>
#include <Rmath.h>
#include <Rdefines.h>
#include <Rinternals.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "uncertain.h"
#include "Rmissing.h"

SEXP test_switch(const SEXP Snps, const SEXP Snps2, const SEXP Split, 
		 const SEXP Prior) {
  
  int *female = NULL;
  SEXP cl = GET_CLASS(Snps);
  if (TYPEOF(cl) != STRSXP) {
    cl = R_data_class(Snps, FALSE); /* S4 way of getting class attribute */
  }
  SEXP diploid = NULL;
  if (!strcmp(CHAR(STRING_ELT(cl, 0)), "XSnpMatrix")) {
    diploid = R_do_slot(Snps, mkString("diploid"));
    female = LOGICAL(diploid);
  }
  unsigned char *snps = RAW(Snps); 
  int N1 = nrows(Snps);
  int M = ncols(Snps);
  int *split = NULL;
  int N2 = 0;
  unsigned char *snps2 = NULL;
  int *female2 = NULL;
  if (TYPEOF(Snps2)==NILSXP) 
    split = INTEGER(Split);
  else {
    N2 = nrows(Snps2);
    snps2 = RAW(Snps2);
    if (female) {
      diploid = R_do_slot(Snps2, mkString("diploid"));
      female2 = LOGICAL(diploid);
    }
  }
    
  double eta = *REAL(Prior);
  
  SEXP Result;
  PROTECT(Result = allocVector(REALSXP, M));
  double *result = REAL(Result);

  double logten = log(10.0);
  for (int j=0; j<M; j++) {
    int n0=0, n1=0, np0=0, np1=0, g=1, N=N1;
    unsigned char *snpsg = snps;
    int *fg = female;
    while (1) {
      R_xlen_t ij = (R_xlen_t)N*(R_xlen_t)j;
      for (int i=0; i<N; i++, ij++) {
	const int s = (int) snpsg[ij];
	if (s) {
	  if (split) 
	    g = split[i];
	  if (g!=NA_INTEGER) {
	    if (female && !female[i]) {
	      if (g==2) {
		n1++;
		np1 += (s-1)/2;
	      }
	      else {
		n0++;
		np0 += (s-1)/2;
	      }
	    }
	    else {
	      if (g==2) {
		n1 += 2;
		np1 += (s-1);
	      }
	      else {
		n0 += 2;
		np0 += (s-1);
	      }
	    }
	  }
	}
      }
      if (split || g==2)
	break;
      else {
	snpsg = snps2;
	if (female) 
	  fg = female2;
	N = N2;
	g = 2;
      }
    }
    /* log base 10 Bayes factor for switch */
    result[j] = logten*(lbeta(eta + (double) (np0+n1-np1),
			      eta + (double) (n0-np0+np1)) -
			lbeta(eta + (double) (np0+np1), 
			      eta + (double) (n0+n1-np0-np1)));
  }
  UNPROTECT(1);
  return Result;
}
	
	
SEXP smat_switch(SEXP X, SEXP Switch) {
  SEXP Result = duplicate(X);
  unsigned char *res = RAW(Result); 
  int N = nrows(Result);
  int nsw = length(Switch);
  int *sw = INTEGER(Switch);
  for (int i=0; i<nsw; i++) {
    unsigned char *resij = res + N*(sw[i] - 1);
    for (int j=0; j<N; j++, resij++) {
      unsigned char g = *resij;
      if (g) {
	if (g>3) {
	  double p0, p1, p2;
	  g2post(g, &p0, &p1, &p2);
	  *resij = post2g(p1, p0);
	}
	else
	  *resij = 4 - g;
      }
    }
  }
  return Result;
}

