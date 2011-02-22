#include <R.h>
#include <Rdefines.h>
#include <Rinternals.h>
#include <stdlib.h>
#include <stdio.h>
#include "uncertain.h"

/* As.numeric function */

SEXP asnum(SEXP Snps) {
  
  int N, M;
  SEXP Result;
  const unsigned char *snps = RAW(Snps);
  SEXP names = getAttrib(Snps, R_DimNamesSymbol);
  if (TYPEOF(names)==NILSXP) {
    N = length(Snps);
    M = 1;
    PROTECT(Result = allocVector(REALSXP, N));
    names = getAttrib(Snps, R_NamesSymbol);
    setAttrib(Result, R_NamesSymbol, names);
  }
  else {
    N = nrows(Snps);
    M = ncols(Snps);
    PROTECT(Result = allocMatrix(REALSXP, N, M));
    SEXP Dimnames;
    PROTECT(Dimnames = allocVector(VECSXP, 2));
    SET_VECTOR_ELT(Dimnames, 0, VECTOR_ELT(names, 0));
    SET_VECTOR_ELT(Dimnames, 1, VECTOR_ELT(names, 1));
    setAttrib(Result, R_DimNamesSymbol, Dimnames);
    UNPROTECT(1);
  }
  double *result = REAL(Result);
  for (int j=0, ij=0; j<M; j++) {
    for (int i=0; i<N; i++, ij++) {
      unsigned char g = snps[ij];
      if (g) {
	double p0, p1, p2;
	g2post(g, &p0, &p1, &p2);
	result[ij] = p1 + 2*p2;
      }
      else
	result[ij] = NA_REAL;
    }
  }
  UNPROTECT(1);
  return Result;
}
    
