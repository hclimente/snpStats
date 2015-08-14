#include <string.h>
#include <stdio.h>
#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include "uncertain.h"

/* Fst calculation -- as used in HapMap, but with option to be more like AOV */

SEXP Fst(SEXP Snps, SEXP Group, SEXP HapMap) {
  
  /* Process Snps argument */

  const char *classS = NULL;

  classS = CHAR(STRING_ELT(getAttrib(Snps, R_ClassSymbol), 0));
  int ifX = 0; /* default Not X */
  if (!strcmp(classS, "SnpMatrix"))
    ifX = 0;
  else if (!strcmp(classS, "XSnpMatrix"))
    ifX = 1;
  else {
    ifX = 0; /* to avoid warning message */
    error("Argument error - class(Snps)");
  }
  if(!IS_S4_OBJECT(Snps)) {
    error("Argument error - Snps is not S4 object");
  }

  if (TYPEOF(Snps)!=RAWSXP)
    error("Argument error - Snps");
  const unsigned char *snps = RAW(Snps);
  int N, nsnp;
  N = nrows(Snps);
  nsnp = ncols(Snps);
  
     
  const int *diploid = NULL;
  if (ifX) {
    SEXP Diploid = R_do_slot(Snps, mkString("diploid")); 
    diploid = LOGICAL(Diploid);
  }

  /* Process Group argument */

  classS = CHAR(STRING_ELT(getAttrib(Group, R_ClassSymbol), 0));
  if (strcmp(classS, "factor"))
    error("Argument error - class(Group)");
  if (LENGTH(Group)!=N)
    error("Non-conformant arguments");
  int ngrp = nlevels(Group);
  int *group = INTEGER(Group);

  /* Process HapMap argument */

  if (TYPEOF(HapMap)!=LGLSXP)
    error("Argument error - typeof(HapMap)");
  int hapmap = *LOGICAL(HapMap);

  /* Setup result object */

  SEXP Fst, Weight;
  PROTECT(Fst = allocVector(REALSXP, nsnp));
  PROTECT(Weight = allocVector(REALSXP, nsnp));
  double *fst = REAL(Fst);
  double *weight = REAL(Weight);

  /* Work arrays */

  int *na2 = (int *) Calloc(ngrp, int);
  int *na = (int *) Calloc(ngrp, int);
  double *gwts = (double *) Calloc(ngrp, double);
  
  /* Calculate group weights */

  memset(na, 0x00, ngrp*sizeof(int));
  for (int i=0; i<N; i++) {
    int gi = group[i];
    if (gi != NA_INTEGER) {
      gi--;
      if (ifX) {
	int fi = diploid[i];
	  na[gi] += fi? 2: 1;
      }
      else {
	na[gi] += 2;
      }	
    }
  }

  double sgw = 0.0;
  for (int g=0; g<ngrp; g++) {
    double w = (double)na[g];
    if (hapmap)
      w = w*(w-1.0);
    sgw += w;
    gwts[g] = w;
  }
  for (int g=0; g<ngrp; g++) 
    gwts[g] /= sgw;


  /* Calculate Fst for each SNP, plus its weight in the overall estimate  */


  SEXP Result, Names;
  PROTECT(Result = allocVector(VECSXP, 2));
  PROTECT(Names = allocVector(STRSXP, 2));
  SET_STRING_ELT(Names, 0, mkChar("Fst"));
  SET_STRING_ELT(Names, 1, mkChar("weight"));
  setAttrib(Result, R_NamesSymbol, Names);
  SET_VECTOR_ELT(Result, 0, Fst);
  SET_VECTOR_ELT(Result, 1, Weight);
  
  UNPROTECT(4);
  return Result;
}

      
