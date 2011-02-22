#include <R.h>
#include <Rinternals.h>
#include "Rmissing.h"
#include <math.h>
#include <stdio.h>
#include <string.h>
#include "hash_index.h"
#include "imputation.h"
#include "uncertain.h"

SEXP X_snp_summary(const SEXP Snps, const SEXP Rules, const SEXP Uncertain) {

  /* SNPs ---- an XSnpMatrix */

  int *ifFemale;
  SEXP Female = R_do_slot(Snps, mkString("Female"));
  if (TYPEOF(Female)!=LGLSXP)
    error("Argument error -  Female slot wrong type");
  ifFemale = LOGICAL(Female);
  if (TYPEOF(Snps)!=RAWSXP)
    error("Argument error - Snps wrong type");
  if (Snps == R_NilValue) {
    error("Argument error - Snps = NULL");
  }
  if(!IS_S4_OBJECT(Snps)) {
    error("Argument error - Snps is not S4 object");
  }
  const unsigned char *snps = RAW(Snps);
  int N, M;
  int *dim = INTEGER(getAttrib(Snps, R_DimSymbol));
  N = dim[0];
  M = dim[1];
  SEXP names = getAttrib(Snps, R_DimNamesSymbol);
  if (names == R_NilValue) {
    error("Argument error - Snps object has no names");
  }
  SEXP snpNames = VECTOR_ELT(names, 1);
  if (snpNames == R_NilValue) {
    error("Argument error - Snps object has no snp names");
  }
   
  /* Rules object */

  int nrules = 0;
  if (!isNull(Rules)) {
    const char *classR = NULL;
    if (TYPEOF(R_data_class(Rules, FALSE)) == STRSXP) {
      classR = CHAR(STRING_ELT(R_data_class(Rules, FALSE), 0));
    } else {
      classR = CHAR(STRING_ELT(getAttrib(Rules, R_ClassSymbol), 0));
    }
    if (strcmp(classR, "ImputationRules")!=0) 
      error("Argument error - Rules");
    nrules = LENGTH(Rules);
  }

  /* Handling of uncertain genotypes */

  if (TYPEOF(Uncertain) != LGLSXP)
    error("Argument error: Uncertain is wrong type");
  int uncert = *LOGICAL(Uncertain);

  /* Output object */

  int ncol = nrules==0? M: nrules;

  SEXP Result, Calls, Call_rate, Certain_call_rate, MAF, 
    P_AA, P_AB, P_BB, P_AY, P_BY, Z_HWE, Calls_female;
  PROTECT(Result = allocVector(VECSXP, 11));
  PROTECT(Calls = allocVector(INTSXP, ncol));
  SET_VECTOR_ELT(Result, 0, Calls);
  PROTECT(Call_rate = allocVector(REALSXP, ncol));
  SET_VECTOR_ELT(Result, 1, Call_rate);
  PROTECT(Certain_call_rate = allocVector(REALSXP, ncol));
  SET_VECTOR_ELT(Result, 2, Certain_call_rate);
  PROTECT(MAF = allocVector(REALSXP, ncol));
  SET_VECTOR_ELT(Result, 3, MAF);
  PROTECT(P_AA = allocVector(REALSXP, ncol));
  SET_VECTOR_ELT(Result, 4, P_AA);
  PROTECT(P_AB = allocVector(REALSXP, ncol));
  SET_VECTOR_ELT(Result, 5, P_AB);
  PROTECT(P_BB = allocVector(REALSXP, ncol));
  SET_VECTOR_ELT(Result, 6, P_BB);
  PROTECT(P_AY = allocVector(REALSXP, ncol));
  SET_VECTOR_ELT(Result, 7, P_AY);
  PROTECT(P_BY = allocVector(REALSXP, ncol));
  SET_VECTOR_ELT(Result, 8, P_BY);
  PROTECT(Z_HWE = allocVector(REALSXP, ncol));
  SET_VECTOR_ELT(Result, 9, Z_HWE);
  PROTECT(Calls_female = allocVector(INTSXP, ncol));
  SET_VECTOR_ELT(Result, 10, Calls_female);
  SEXP Names;
  PROTECT(Names = allocVector(STRSXP, 11));
  SET_STRING_ELT(Names, 0, mkChar("Calls"));
  SET_STRING_ELT(Names, 1, mkChar("Call.rate"));
  SET_STRING_ELT(Names, 2, mkChar("Certain.calls"));
  SET_STRING_ELT(Names, 3, mkChar("MAF"));
  SET_STRING_ELT(Names, 4, mkChar("P.AA"));
  SET_STRING_ELT(Names, 5, mkChar("P.AB"));
  SET_STRING_ELT(Names, 6, mkChar("P.BB"));
  SET_STRING_ELT(Names, 7, mkChar("P.AY"));
  SET_STRING_ELT(Names, 8, mkChar("P.BY"));
  SET_STRING_ELT(Names, 9, mkChar("z.HWE"));
  SET_STRING_ELT(Names, 10, mkChar("Calls.female"));
  int *calls = INTEGER(Calls);
  double *call_rate = REAL(Call_rate);
  double *certain = REAL(Certain_call_rate);
  double *maf = REAL(MAF);
  double *p_aa = REAL(P_AA);
  double *p_ab = REAL(P_AB);
  double *p_bb = REAL(P_BB);
  double *p_ay = REAL(P_AY);
  double *p_by = REAL(P_BY);
  double *z_hwe = REAL(Z_HWE);
  int *calls_female = INTEGER(Calls_female);
  setAttrib(Result, R_NamesSymbol, Names);
  setAttrib(Result, R_RowNamesSymbol, snpNames);
  SEXP dfClass;
  PROTECT(dfClass = allocVector(STRSXP, 1));
  SET_STRING_ELT(dfClass, 0, mkChar("data.frame"));
  setAttrib(Result, R_ClassSymbol, dfClass);

  /* If imputation rules, set work arrays and hash index */

  index_db name_index = NULL;
  GTYPE **gt2ht = NULL;
  double *add, *dom;
  int pmax;
  if (nrules!=0) {
    name_index = create_name_index(snpNames);
    add = (double *)Calloc(N, double);
    dom = (double *)Calloc(N, double);
    pmax = *INTEGER(getAttrib(Rules, install("Max.predictors")));
    gt2ht = (GTYPE **)Calloc(pmax, GTYPE *);
    for (int i=0; i<pmax; i++)
      gt2ht[i] = create_gtype_table(i+1);
  }

  /* Calculations */

  int *obs = (int *) Calloc(N, int);
  int i, j, ij;
  for (i=0; i<N; i++)
    obs[i] = 0;

  for (j=0, ij=0; j<ncol; j++) { 
    int aa = 0, ab = 0, bb = 0, ay=0, by=0;
    double aap = 0.0, abp = 0.0, bbp = 0.0, ayp=0.0, byp=0.0;
    int ncall=0, ncertain=0;
    if (nrules!=0) {
      SEXP Rule =  VECTOR_ELT(Rules, j);
      if (!isNull(Rule)){ /* Monomorphic */
	do_impute(Snps, N, NULL, NULL, N, name_index, Rule, gt2ht, 
		  add, dom);
	for (int i=0; i<N; i++) {
	  double addi = add[i];
	  double domi = dom[i];
	  if (!ISNA(addi) && !ISNA(domi)) {
	    ncall++;
	    obs[i] = 1;
	    double p2 = domi;
	    double p1 = addi-2.0*p2;
	    double p0 = 1.0-p1-p2;
	    if (ifFemale[i]) {
	      aap += p0;
	      abp += p1;
	      bbp += p2;
	    }
	    else {
	      aap += p0;
	      bbp += p2;
	    }
	  }
	}
      }
    } else {
      for (i=0; i<N; i++) {
	int g = (int) snps[ij++];
	if (g) {
	  ncall++;
	  obs[i] = 1;
	  if (g<4) {
	    ncertain++;
	    if (ifFemale[i]) 
	      switch (g) {
	      case 1: aa++; break;
	      case 2: ab++; break;
	      case 3: bb++; 
	      }	  
	    else 
	      switch (g) {
	      case 1: ay++; break;
	      case 3: by++; 
	      }
	  }
	  else if (uncert) {
	    double p0, p1, p2;
	    g2post(g, &p0, &p1, &p2);
	    if (ifFemale[i]) {
	      aap += p0;
	      abp += p1;
	      bbp += p2;
	    }
	    else {
	      aap += p0;
	      bbp += p2;
	    }
	  }
	}
      }
    }
    
    /* HWE test only involves certain assignments in females */

    double nv = aa + ab + bb;
    double na = 2*aa + ab;
    double p = na/(2.0*nv);
    double q = 1.0 - p;
    double den = 2*p*q*sqrt(nv);
    double z = den>0.0? (ab - 2*p*q*nv)/den: NA_REAL;
    z_hwe[j] = z;

    /* Call rate stuff */

    calls[j] = ncall;
    calls_female[j] = nv;
    call_rate[j] = (double) ncall/(double) N;
    certain[j] = ncall>0? (double) ncertain /(double) ncall: NA_REAL;

    /* Frequencies */

    nv += aap + abp + bbp;
    na += 2*aap + abp;
    double ny = ay + by + ayp + byp;
    double nc = 2*nv + ny;
    p = (na + ay + ayp)/nc;
    if (p>0.5)
      p = 1.0 - p;
    maf[j] = nc>0? p: NA_REAL;
    p_aa[j] = nv>0? ((double) aa + aap)/nv: NA_REAL;
    p_ab[j] = nv>0? ((double) ab + abp)/nv: NA_REAL;
    p_bb[j] = nv>0? ((double) bb + bbp)/nv: NA_REAL;
    p_ay[j] = ny>0? ((double) ay + ayp)/ny: NA_REAL;
    p_by[j] = ny>0? ((double) by + ayp)/ny: NA_REAL;
  }

  /* Call rate */

  int Nobs = 0;
  for (i=0; i<N; i++) 
    Nobs += obs[i];
  if (Nobs < N) {
    warning("%d rows were empty - ignored when calculating call rates", 
	    N - Nobs);
    double infl = (double) N / (double) Nobs;
    if (Nobs) {
      for (j=0; j<ncol; j++)
	call_rate[j] *= infl;
    }
    else {
      error("Empty matrix");
    } 
  }
      
  /* Tidy up */

  UNPROTECT(14);

  Free(obs);
  if (nrules!=0) {
    index_destroy(name_index);
    Free(add);
    Free(dom);
    for (int i=0; i<pmax; i++) 
      destroy_gtype_table(gt2ht[i], i+1);
    Free(gt2ht);
  }

  return Result;
}

SEXP snp_summary(const SEXP Snps, const SEXP Rules, const SEXP Uncertain) {

  /* SNPs ---- a SnpMatrix */

  if (TYPEOF(Snps)!=RAWSXP)
    error("Argument error - Snps wrong type");
  if (Snps == R_NilValue) {
    error("Argument error - Snps = NULL");
  }
  if(!IS_S4_OBJECT(Snps)) {
    error("Argument error - Snps is not S4 object");
  }
  const unsigned char *snps = RAW(Snps);
  int N, M;
  int *dim = INTEGER(getAttrib(Snps, R_DimSymbol));
  N = dim[0];
  M = dim[1];
  SEXP names = getAttrib(Snps, R_DimNamesSymbol);
  if (names == R_NilValue) {
    error("Argument error - Snps object has no names");
  }
  SEXP snpNames = VECTOR_ELT(names, 1);
  if (snpNames == R_NilValue) {
    error("Argument error - Snps object has no snp names");
  }
  
  /* Rules object */

  int nrules = 0;
  if (!isNull(Rules)) {
    const char *classR = NULL;
    if (TYPEOF(R_data_class(Rules, FALSE)) == STRSXP) {
      classR = CHAR(STRING_ELT(R_data_class(Rules, FALSE), 0));
    } else {
      classR = CHAR(STRING_ELT(getAttrib(Rules, R_ClassSymbol), 0));
    }
    if (strcmp(classR, "ImputationRules")!=0) 
      error("Argument error - Rules");
    nrules = LENGTH(Rules);
  }

  /* Handling of uncertain genotypes */

  if (TYPEOF(Uncertain) != LGLSXP)
    error("Argument error: Uncertain is wrong type");
  int uncert = *LOGICAL(Uncertain);

  /* Output object */

  int ncol = nrules==0? M: nrules;

  SEXP Result, Calls, Call_rate, Certain_call_rate, MAF, 
    P_AA, P_AB, P_BB, Z_HWE;
  PROTECT(Result = allocVector(VECSXP, 8));
  PROTECT(Calls = allocVector(INTSXP, ncol));
  SET_VECTOR_ELT(Result, 0, Calls);
  PROTECT(Call_rate = allocVector(REALSXP, ncol));
  SET_VECTOR_ELT(Result, 1, Call_rate);
  PROTECT(Certain_call_rate = allocVector(REALSXP, ncol));
  SET_VECTOR_ELT(Result, 2, Certain_call_rate);
  PROTECT(MAF = allocVector(REALSXP, ncol));
  SET_VECTOR_ELT(Result, 3, MAF);
  PROTECT(P_AA = allocVector(REALSXP, ncol));
  SET_VECTOR_ELT(Result, 4, P_AA);
  PROTECT(P_AB = allocVector(REALSXP, ncol));
  SET_VECTOR_ELT(Result, 5, P_AB);
  PROTECT(P_BB = allocVector(REALSXP, ncol));
  SET_VECTOR_ELT(Result, 6, P_BB);
  PROTECT(Z_HWE = allocVector(REALSXP, ncol));
  SET_VECTOR_ELT(Result, 7, Z_HWE);

  SEXP Names;
  PROTECT(Names = allocVector(STRSXP, 8));
  SET_STRING_ELT(Names, 0, mkChar("Calls"));
  SET_STRING_ELT(Names, 1, mkChar("Call.rate"));
  SET_STRING_ELT(Names, 2, mkChar("Certain.calls"));
  SET_STRING_ELT(Names, 3, mkChar("MAF"));
  SET_STRING_ELT(Names, 4, mkChar("P.AA"));
  SET_STRING_ELT(Names, 5, mkChar("P.AB"));
  SET_STRING_ELT(Names, 6, mkChar("P.BB"));
  SET_STRING_ELT(Names, 7, mkChar("z.HWE"));

  int *calls = INTEGER(Calls);
  double *call_rate = REAL(Call_rate);
  double *certain = REAL(Certain_call_rate);
  double *maf = REAL(MAF);
  double *p_aa = REAL(P_AA);
  double *p_ab = REAL(P_AB);
  double *p_bb = REAL(P_BB);
  double *z_hwe = REAL(Z_HWE);

  setAttrib(Result, R_NamesSymbol, Names);
  setAttrib(Result, R_RowNamesSymbol, duplicate(snpNames));
  SEXP dfClass;
  PROTECT(dfClass = allocVector(STRSXP, 1));
  SET_STRING_ELT(dfClass, 0, mkChar("data.frame"));
  setAttrib(Result, R_ClassSymbol, dfClass);

  /* If imputation rules, set work arrays and hash index */

  index_db name_index = NULL;
  GTYPE **gt2ht = NULL;
  double *add, *dom;
  int pmax;
  if (nrules!=0) {
    name_index = create_name_index(snpNames);
    add = (double *)Calloc(N, double);
    dom = (double *)Calloc(N, double);
    pmax = *INTEGER(getAttrib(Rules, install("Max.predictors")));
    gt2ht = (GTYPE **)Calloc(pmax, GTYPE *);
    for (int i=0; i<pmax; i++)
      gt2ht[i] = create_gtype_table(i+1);
   }

  /* Calculations */

  int *obs = (int *) Calloc(N, int);
  int i, j, ij;
  for (i=0; i<N; i++)
    obs[i] = 0;
  for (j=0, ij=0; j<ncol; j++) {
    int ncall=0, ncertain=0; 
    int aa = 0, ab = 0, bb = 0;
    double aap=0.0, abp=0.0, bbp=0.0;
    if (nrules!=0) {
      SEXP Rule =  VECTOR_ELT(Rules, j);
      if (!isNull(Rule)){ /* Monomorphic */
	do_impute(Snps, N, NULL, NULL, N, name_index, Rule, gt2ht, 
		  add, dom);
	for (int i=0; i<N; i++) {
	  double addi = add[i];
	  double domi = dom[i];
	  if(!ISNA(addi) && !ISNA(domi)) {
	    ncall++;
	    obs[i] = 1;
	    double p2 = domi;
	    double p1 = addi-2.0*p2;
	    double p0 = 1.0-p1-p2;
	    aap += p0;
	    abp += p1;
	    bbp += p2;
	  }
	}
      }
    } else {
      for (i=0; i<N; i++) {
	int g = (int) snps[ij++];
	if (g) {
	  ncall++;
	  obs[i] = 1;
	  if (g<4) {
	    ncertain++;
	    switch (g) {
	    case 1: aa++; break;
	    case 2: ab++; break;
	    case 3: bb++; 	
	    } 
	  }
	  else if (uncert) {
	    double p0, p1, p2;
	    g2post(g, &p0, &p1, &p2);
	    aap += p0;
	    abp += p1;
	    bbp += p2;
	  }
	}
      }
    }
    
    /* HWE test in certain assignments only */

    double nv = aa + ab + bb;
    double na = 2*aa + ab;
    double p =  na/(2*nv);
    double q = 1.0 - p;
    double den = 2*p*q*sqrt(nv);
    double z = den>0.0? (ab - 2*p*q*nv)/den: NA_REAL;
    z_hwe[j] = z;

    /* Call rates */

    calls[j] = ncall;
    call_rate[j] = (double) ncall/(double) N;
    certain[j] = ncall>0? (double) ncertain /(double) ncall: NA_REAL;

    /* Frequencies */

    nv += aap + abp + bbp;
    na += 2*aap + abp;
    p =  na/(2*nv);
    if (p>0.5)
      p = 1.0 - p;
    maf[j] = nv>0? p: NA_REAL;
    p_aa[j] = nv>0? ((double) aa + aap)/nv: NA_REAL;
    p_ab[j] = nv>0? ((double) ab + abp)/nv: NA_REAL;
    p_bb[j] = nv>0? ((double) bb + bbp)/nv: NA_REAL;
  }

  /* Call rate */

  int Nobs = 0;
  for (i=0; i<N; i++) 
    Nobs += obs[i];
  if (Nobs < N) {
    warning("%d rows were empty - ignored when calculating call rates", 
	    N - Nobs);
    double infl = (double) N / (double) Nobs;
    if (Nobs) {
      for (j=0; j<ncol; j++)
	call_rate[j] *= infl;
    }
    else {
      error("Empty matrix");
    }
  }

  /* Tidy up */

  UNPROTECT(11);

  Free(obs);
  if (nrules!=0) {
    index_destroy(name_index);
    Free(add);
    Free(dom);
    for (int i=0; i<pmax; i++) 
      destroy_gtype_table(gt2ht[i], i+1);
    Free(gt2ht);
  }

  return Result;
}


SEXP row_summary(const SEXP Snps) {

  if (TYPEOF(Snps)!=RAWSXP)
    error("Argument error - Snps wrong type");
  if (Snps == R_NilValue) {
    error("Argument error - Snps = NULL");
  }
  if(!IS_S4_OBJECT(Snps)) {
    error("Argument error - Snps is not S4 object");
  }
  const unsigned char *snps = RAW(Snps);
  int N, M;
  int *dim = INTEGER(getAttrib(Snps, R_DimSymbol));
  N = dim[0];
  M = dim[1];
  SEXP names = getAttrib(Snps, R_DimNamesSymbol);
  if (names == R_NilValue) {
    error("Argument error - Snps object has no names");
  }
  SEXP rowNames = VECTOR_ELT(names, 0);
  if (rowNames == R_NilValue) {
    error("Argument error - Snps object has no row names");
  }
  
  /* Output object */

  SEXP Result, Call_rate, Certain_call_rate, Het;
  PROTECT(Result = allocVector(VECSXP, 3));
  PROTECT(Call_rate = allocVector(REALSXP, N));
  SET_VECTOR_ELT(Result, 0, Call_rate);
  PROTECT(Certain_call_rate = allocVector(REALSXP, N));
  SET_VECTOR_ELT(Result, 1, Certain_call_rate);
  PROTECT(Het = allocVector(REALSXP, N));
  SET_VECTOR_ELT(Result, 2, Het);

  SEXP Names;
  PROTECT(Names = allocVector(STRSXP, 3));
  SET_STRING_ELT(Names, 0, mkChar("Call.rate"));
  SET_STRING_ELT(Names, 1, mkChar("Certain.calls"));
  SET_STRING_ELT(Names, 2, mkChar("Heterozygosity"));

  double *call_rate = REAL(Call_rate);
  double *certain = REAL(Certain_call_rate);
  double *het = REAL(Het);

  setAttrib(Result, R_NamesSymbol, Names);
  setAttrib(Result, R_RowNamesSymbol, duplicate(rowNames));
  SEXP dfClass;
  PROTECT(dfClass = allocVector(STRSXP, 1));
  SET_STRING_ELT(dfClass, 0, mkChar("data.frame"));
  setAttrib(Result, R_ClassSymbol, dfClass);

  /* Calculations */
  int i, j, ij;
  for (i=0; i<N; i++) {
    int ncall = 0, nhet=0, ncertain=0;
    for (j=0, ij=i; j<M; j++, ij+=N) {
      unsigned char g = snps[ij];
      if (g) {
	ncall++;
	if (g<4) {
	  ncertain++;
	  if (g == 0x02) /* Het */
	    nhet++;
	}
      }
    }
    call_rate[i] = (double) ncall/ (double) M;
    certain[i] = ncall>0? (double) ncertain/(double) ncall: NA_REAL;
    het[i] = ncall>0? (double) nhet/ (double) ncall: NA_REAL;
  }
  UNPROTECT(6);
  return Result;
}
