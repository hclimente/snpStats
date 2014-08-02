#include <R.h>
#include <Rinternals.h>
#include "uncertain.h"

SEXP Rmean2g(SEXP Mean, SEXP MaxE) {
  if (TYPEOF(MaxE)!=LGLSXP)
    error("maxE argument not of type logical");
  int maxE = *LOGICAL(MaxE);
  if (TYPEOF(Mean)!=REALSXP)
    error("argument is not of type numeric");
  int N = length(Mean);
  double *mean = REAL(Mean);

  SEXP Res;
  PROTECT(Res = allocVector(RAWSXP, N));
  unsigned char *res = RAW(Res);

  for (int i=0; i<N; i++)
    res[i] = mean2g(mean[i], maxE);

  UNPROTECT(1);
  return(Res);
}

SEXP Rg2post(SEXP G, SEXP Trans) {
  if (TYPEOF(G)!=RAWSXP)
    error("argument is not of type raw");
  int N = length(G);
  unsigned char *g = RAW(G);
  if (TYPEOF(Trans)!=LGLSXP)
    error("transpose argument not of type logical");
  int *trans = LOGICAL(Trans);

  SEXP Res;
  if (*trans) {
    PROTECT(Res = allocMatrix(REALSXP, 3, N));
    double *res = REAL(Res);
    double *rij = res; 
    for (int i=0; i<N; i++) {
      double *raa = rij++;
      double *rab = rij++;
      double *rbb = rij++;
      if (!g2post(g[i], raa, rab, rbb))
	*raa = *rab = *rbb = NA_REAL;
    }
  }
  else {
    PROTECT(Res = allocMatrix(REALSXP, N, 3));
    double *res = REAL(Res);
    for (int i=0; i<N; i++) {
      double *raa = res+i;
      double *rab = raa+N;
      double *rbb = rab+N;
      if (!g2post(g[i], raa, rab, rbb))
	*raa = *rab = *rbb = NA_REAL;
    }
  }
  UNPROTECT(1);
  return(Res);
}

SEXP Rpost2g(SEXP Pos, SEXP Trans) {
  if (TYPEOF(Pos)!=REALSXP || !isMatrix(Pos))
    error("argument is not a numeric matrix");
  double *pos = REAL(Pos);
  if (TYPEOF(Trans)!=LGLSXP)
    error("transpose argument is not of type logical");
  int *trans = LOGICAL(Trans);

  SEXP Res;
  if (*trans) {
    if (nrows(Pos)!=3)
      error("matrix does not have 3 rows");
    int N = ncols(Pos);
    PROTECT(Res = allocVector(RAWSXP, N));
    unsigned char *res = RAW(Res);
    double *pij = pos;
    for (int i=0; i<N; i++) {
      pij++;
      double pAB = *(pij++);
      double pBB = *(pij++);
      res[i] = post2g(pAB, pBB);
    }
  }
  else {
    if (ncols(Pos)!=3)
      error("matrix does not have 3 columns");
    int N = nrows(Pos);
    PROTECT(Res = allocVector(RAWSXP, N));
    unsigned char *res = RAW(Res);
    for (int i=0; i<N; i++) {
      double *paa = pos+i;
      double *pab = paa+N;
      double *pbb = pab+N;
      res[i] = post2g(*pab, *pbb);
    }
  }
    
  UNPROTECT(1);
  return(Res);
  }










