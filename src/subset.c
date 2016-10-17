#include <R.h>
#include <Rinternals.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* Subset of a SnpMatrix or XSnpMatrix object */

SEXP subset(SEXP X, SEXP Rows, SEXP Cols) {
  /* Attributes of input object */
  const char *classX = CHAR(STRING_ELT(getAttrib(X, R_ClassSymbol), 0));
  int *diploid = NULL;
  SEXP Diploid = R_NilValue;
  if (!strcmp(classX, "XSnpMatrix")) {
    Diploid = R_do_slot(X, mkString("diploid"));
    diploid = LOGICAL(Diploid);
  }
  int *dim = INTEGER(getAttrib(X, R_DimSymbol));
  int N = dim[0];
  int M = dim[1];
  SEXP Dimnames = getAttrib(X, R_DimNamesSymbol);
  SEXP Rownames = VECTOR_ELT(Dimnames, 0);
  SEXP Colnames = VECTOR_ELT(Dimnames, 1);
  int nrows = LENGTH(Rows);
  int *rows = NULL;
  if (nrows)
    rows = INTEGER(Rows);
  else
    nrows = N;
  int ncols = LENGTH(Cols);
  int *cols = NULL;
  if (ncols)
    cols = INTEGER(Cols);
  else 
    ncols = M;
  if (!cols && !rows) {
    warning("No selection made");
    return X;
  }

  const Rbyte *x = RAW(X);

  /* Result */

  int nprotect = 5;
  SEXP Result, Rdim, Rdimnames, Rrownames, Rcolnames, Rdiploid, Rclass, Package;
  PROTECT(Result = allocMatrix(RAWSXP, nrows, ncols));
  Rbyte *r = RAW(Result);
  PROTECT(Rclass = allocVector(STRSXP, 1));
  if (diploid)
    SET_STRING_ELT(Rclass, 0, mkChar("XSnpMatrix"));
  else
    SET_STRING_ELT(Rclass, 0, mkChar("SnpMatrix"));
  PROTECT(Package = allocVector(STRSXP, 1));
  SET_STRING_ELT(Package, 0, mkChar("snpStats"));
  setAttrib(Rclass, install("package"), Package);
  setAttrib(Result, R_ClassSymbol, Rclass);
  SET_S4_OBJECT(Result);
  PROTECT(Rdim = allocVector(INTSXP, 2));
  int *rdim = INTEGER(Rdim);
  rdim[0] = nrows;
  rdim[1] = ncols;
  setAttrib(Result, R_DimSymbol, Rdim);
  PROTECT(Rdimnames = allocVector(VECSXP, 2));
  if (rows) {
    PROTECT(Rrownames = allocVector(STRSXP, nrows));
    nprotect++;
    SET_VECTOR_ELT(Rdimnames, 0, Rrownames);
  }
  else {
    SET_VECTOR_ELT(Rdimnames, 0, duplicate(VECTOR_ELT(Dimnames, 0)));
  }
  if (cols) {
    PROTECT(Rcolnames = allocVector(STRSXP, ncols));
    nprotect++;
    SET_VECTOR_ELT(Rdimnames, 1, Rcolnames);
  }
  else {
    SET_VECTOR_ELT(Rdimnames, 1, duplicate(VECTOR_ELT(Dimnames, 1)));
  }
  setAttrib(Result, R_DimNamesSymbol, Rdimnames);
  int* rdiploid = NULL;
  if (diploid) {
    if (rows) {
      PROTECT(Rdiploid = allocVector(LGLSXP, nrows));
      nprotect++;
      rdiploid = LOGICAL(Rdiploid);
      R_do_slot_assign(Result, mkString("diploid"), Rdiploid);
    }
    else {
      R_do_slot_assign(Result, mkString("diploid"), duplicate(Diploid));
    }
  }

  /* Populate the new object */
  
  R_xlen_t uv1 = -N, ij = 0;
  for (int j=0; j<ncols; j++) {
    if (cols) {
      int v = cols[j]-1;
      uv1 = (R_xlen_t)N * (R_xlen_t)v;
      SET_STRING_ELT(Rcolnames, j, duplicate(STRING_ELT(Colnames, v))); 
    } 
    else {
      uv1 += N;
    }
    if (rows) {
      for (int i=0; i<nrows; i++) {
        int u = rows[i]-1;
        R_xlen_t uv = uv1 + u;
        r[ij++] = x[uv];
      }
    }
    else {
      memcpy(r+ij, x+uv1, N);
      ij += N;
    }
  }
  if (rows) {
    for (int i=0; i<nrows; i++) {
      int u = rows[i]-1;
      SET_STRING_ELT(Rrownames, i, duplicate(STRING_ELT(Rownames, u)));
      if (diploid)
        rdiploid[i] = diploid[u];
    }
  }
 
  if (ij > R_LEN_T_MAX)
    warning("Output  SnpMatrix has more than 2^31-1 elements. Many functions do not support such objects");

  UNPROTECT(nprotect);
  return Result;
  
}
