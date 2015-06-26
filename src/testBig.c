#include <R.h>
#include <Rinternals.h>
#include <stdio.h>
#include <stdlib.h>

/* Generate a (possibly large) object of class SnpMatrix */

SEXP snp_big(SEXP nrow, SEXP ncol) {
  int N = *INTEGER(nrow);
  int M = *INTEGER(ncol);

  R_xlen_t size = (R_xlen_t)N * (R_xlen_t)M;
  SEXP Result, Class, Rnames, Cnames, Dnames, Package;
  PROTECT(Result = allocMatrix(RAWSXP, N, M));
  PROTECT(Class = allocVector(STRSXP, 1));
  SET_STRING_ELT(Class, 0, mkChar("SnpMatrix"));
  PROTECT(Package = allocVector(STRSXP, 1));
  SET_STRING_ELT(Package, 0, mkChar("snpStats"));
  setAttrib(Class, install("package"), Package);
  classgets(Result, Class);
  SET_S4_OBJECT(Result);

  PROTECT(Dnames = allocVector(VECSXP, 2));
  setAttrib(Result, R_DimNamesSymbol, Dnames);
  PROTECT(Rnames = allocVector(STRSXP, N));
  PROTECT(Cnames = allocVector(STRSXP, M));
  SET_VECTOR_ELT(Dnames, 0, Rnames);
  SET_VECTOR_ELT(Dnames, 1, Cnames);
  for (int i=0; i<N; i++) {
    char name[16];
    sprintf(name, "Subject%d", i+1);
    SET_STRING_ELT(Rnames, i, mkChar(name));
  }
  for (int i=0; i<M; i++) {
    char name[16];
    sprintf(name, "Snp%d", i+1);
    SET_STRING_ELT(Cnames, i, mkChar(name));
  }
  char *Ri = RAW(Result);
  for (R_xlen_t i; i<size; i++) {
    int x = 1+rand()%2+rand()%2;
    Ri[i] = x;
  }
  UNPROTECT(6);
  return Result;
}

