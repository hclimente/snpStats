/* 
  If A is a SnpMatrix this routine calculates the matrix B.B-transpose, 
  where B is derived from A by normalising columns to have zero mean and 
  unit standard deviation under HWE. That is, if p is the allele frequency 
  for one column of A, the corresponding column of B (if elements are
  coded 0, 1, or 2) is obtained by subtracting the mean, 2*p and dividing by 
  the SD, sqrt(2*p*(1-p)). For male samples and the X chromosome, codes are 
  0 and 2 so that the mean is again 2*p, but the SD is now 2*sqrt(p*(1-p)).
  Missing genotypes score zero (equivalent to replacing missing values by 
  the mean in the original matrix.

  Missing data treatment:

  If Correct_for_missing is FALSE, missing genotypes are replaced by 
  their (marginal) expectations - i.e. twice the allele frequency

  If TRUE, contributions to the output matrix are weighted using inverse
  probability weights. The (small) probability that locus k is missing in 
  subject i is assumed to be mu*alpha_i*beta_k where alpha and beta are 
  are vectors with mean 1. Then the probability that subject i is observed 
  at locus k is 1-mu*alpha_i*beta_k, and the probability that locus k 
  is observed in both subject i and subject j is 
  (1-mu*alpha_i*beta_k)*(1-mu*alpha_j*beta_k) 

  alpha_i and beta_k are estimated from the observed numbers of missing calls 
  as follows:

  T_ik = 1 if call (i,k) is missing, 0 otherwise

  alpha_i = N*T_i./T_..
  beta_k = M*T_.k/T..
  mu = T_../(N*M)

  so that

  mu*alpha_i*beta_k = T_i.*T_.k/T_..

  If Uncertain==TRUE, uncertain genotypes are replace by posterior 
  expectations. Otherwise they are treated as missing.

*/

#include <R.h>
#include <Rmath.h>
#include <Rdefines.h>
#include <Rinternals.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "uncertain.h"
#include "Rmissing.h"

SEXP xxt(const SEXP Snps, const SEXP Strata, const SEXP Correct_for_missing, 
	 const SEXP Lower_only, const SEXP Uncertain) {
  
  if (TYPEOF(Correct_for_missing)!=LGLSXP)
    error("Argument error - Correct_for_missing wrong type");
  const int correct = *LOGICAL(Correct_for_missing);

  if (TYPEOF(Lower_only)!=LGLSXP)
    error("Argument error - Lower_only wrong type");
  const int lower = *LOGICAL(Lower_only);

  int *ifFemale = NULL;
  SEXP cl = GET_CLASS(Snps);
  if (TYPEOF(cl) != STRSXP) {
    cl = R_data_class(Snps, FALSE); /* S4 way of getting class attribute */
  }
  if (!strcmp(CHAR(STRING_ELT(cl, 0)), "XSnpMatrix")) {
    SEXP Female = R_do_slot(Snps, mkString("Female"));
    if (TYPEOF(Female)!=LGLSXP)
      error("Argument error -  Female slot wrong type");
    ifFemale = LOGICAL(Female);
  }
  else if (strcmp(CHAR(STRING_ELT(cl, 0)), "SnpMatrix")) {
    error("Argument error - Snps wrong type");
  }    

  const unsigned char *snps = RAW(Snps);
  int N, M;
  int *dim = INTEGER(getAttrib(Snps, R_DimSymbol));
  N = dim[0];
  M = dim[1];

  int *strata = NULL;
  double *mu = NULL, *sd = NULL;
  int *count = NULL; 
  double *acount = NULL;
  int nstrata = 0;
  if (!isNull(Strata)) {
    if (LENGTH(Strata)!=N)
      error("Argument error - Strata argument wrong length");
    strata = INTEGER(Strata);
    nstrata = nlevels(Strata);
    count = Calloc(nstrata, int);
    acount = (double *)Calloc(nstrata, double);
    mu = Calloc(nstrata, double);
    sd = Calloc(nstrata, double);
  }
  double mean=0.0, sd2=0.0; 

  /* Weights for missing data correction */

  int *Ti=NULL, *Tk=NULL, T=0;
  if (correct) {
    warning("With correct.for.missing option set, result may not be a positive semi-definite matrix");
    T = 0;
    Ti = Calloc(N, int);
    memset(Ti, 0x00, N*sizeof(int));
    Tk = Calloc(M, int);
    memset(Tk, 0x00, M*sizeof(int));
    for (int k=0, ik=0; k<M; k++) {
      for (int i=0; i<N; i++) {
	int sik = (int) snps[ik++];
	if (!sik){
	  T++;
	  Ti[i]++;
	  Tk[k]++;
	}
      }
    } 
  }


  /* Handling of uncertain genotypes */

  if (TYPEOF(Uncertain) != LGLSXP)
    error("Argument error: Uncertainty is wrong type");
  int uncert = *LOGICAL(Uncertain);

  /* Result matrix */

  SEXP Result;
  PROTECT(Result = allocMatrix(REALSXP, N, N));
  double *result = REAL(Result);
  memset(result, 0x00, N*N*sizeof(double));

  /* Update result matrix for each locus in turn */

  const double rt2 = sqrt(2.0);

  for (int ik=0, k=0; k<M; k++) {

    /* Calculate allele frequency */

    int s1=0, s2=0, polymorphic=0;
    if (strata) {
      memset(count, 0x00, nstrata*sizeof(int));
      memset(acount, 0x00, nstrata*sizeof(double));
      for (int ki=ik, i=0; i<N; i++) {
	unsigned char w = snps[ki++];
	if (w&&((w<4)|uncert)) {
	  int si = strata[i]-1;
	  double gm = g2mean(w);
	  if (ifFemale && !ifFemale[i]) {
	    count[si]++;
	    acount[si] += gm/2.0;
	  }
	  else {
	    count[si] += 2;
	    acount[si] += gm;
	  }
	}
      }
      for (int i=0; i<nstrata; i++) {
	int s1 = count[i];
	double s2 = acount[i];
	if (s1 && (s2>0.0) && (s2<s1)){
	  polymorphic = 1;
	  double afk = s2/(double) s1;
	  mu[i] = 2.0*afk + 1.0;
	  sd[i] = sqrt(2.0*afk*(1.0-afk));
	}
	else 
	  mu[i] = sd[i] = 0.0;
      }
    }
    else {
      for (int ki=ik, i=0; i<N; i++) {
	unsigned char w =  snps[ki++];
	if (w&&((w<4)|uncert)) {
	  double gm = g2mean(w);
	  if (ifFemale && !ifFemale[i]) {
	    s1++;
	    s2 += gm/2.0;
	  }
	  else {
	    s1 += 2;
	    s2 += gm;
	  }
	}
      }
      polymorphic = (s1 && (s2>0.0) && (s2<s1));
      if (polymorphic){
	double afk = ((double) s2) / ((double) s1);
	mean = 2.0*afk;
	sd2 = sqrt(2.0*afk*(1.0-afk));
      }
      else
	mean = sd2 = 0.0;
    }
    
    /* If polymorphic, add contribution */

    if (polymorphic) {
      double tk=0.0, ipw=0.0;
      if (correct)
	tk = T? (double) Tk[k] / (double) T: 0.0;
    
      /* Update X.X-transpose matrix */

      for (int i=0, ij=0; i<N; i++, ik++) {
	if (strata) {
	  int si = strata[i]-1;
	  mean = mu[si];
	  sd2 = sd[si];
	}
	unsigned char sik = snps[ik];
	if (sik && sd2!=0.0) {
	  double xik = g2mean(sik) - mean;
	  xik /= (ifFemale && !ifFemale[i])? (rt2*sd2): sd2;
	  if (correct) {
	    ipw = 1.0/(1.0 - tk*(double)Ti[i]); /* IPW for diagonal */
	  }
	  ij += i;
	  for (int jk=ik, j=i; j<N; j++, ij++) {
	    if (strata) {
	      int sj = strata[j]-1;
	      mean = mu[sj];
	      sd2 = sd[sj];
	    }
	    unsigned char sjk = snps[jk++]; 
	    if (sjk && sd2!=0.0) {
	      double xjk = g2mean(sjk) - mean;
	      xjk/=  (ifFemale && !ifFemale[j])? (rt2*sd2): sd2;
	      if (correct) 
		result[ij] += xik*xjk*((i==j)?ipw: ipw/(1.0-tk*(double)Ti[j]));
	      else 
		result[ij] += xik*xjk;
	    }
	  }
	}
	else {
	  ij += N;
	}
      }
    }
    else {
      ik += N;
    }
  }

  /* Copy lower triangle to upper triangle */

  if (!lower) {
    for (int i=0, ij=0; i<N; i++) {
      ij += (i+1);
      for (int j=i+1, ji=ij-1+N; j<N; j++, ij++, ji+=N) {
	result[ji] = result[ij];
      }
    }
  }

  /* Return work space */

  if (correct) {
    Free(Tk);
    Free(Ti);
  }
  if (strata) {
    Free(acount);
    Free(count);
    Free(mu);
    Free(sd);
  }
      
  UNPROTECT(1);
  return(Result);
}
  

/* 
   Correlations between columns of a snpmatrix and columns of a normal matrix
*/

 
SEXP corsm(const SEXP Snps, const SEXP X, const SEXP Uncertain) { 

  if (!inherits(Snps, "SnpMatrix"))
    error("Argument error - Snps wrong type");
  const unsigned char *snps = RAW(Snps);
  int *dim = INTEGER(getAttrib(Snps, R_DimSymbol));
  int N = dim[0];
  int M = dim[1];

  if (TYPEOF(X)!=REALSXP)
    error("Argument error - X wrong type");
  if (X == R_NilValue) {
    error("Argument error - X = NULL");
  }
  const double *x = REAL(X);
  dim = INTEGER(getAttrib(X, R_DimSymbol));
  if (dim[0] != N) {
    error("Unequal numbers of rows");
  }
  int P = dim[1];

  /* Handling of uncertain genotypes */

  if (TYPEOF(Uncertain) != LGLSXP)
    error("Argument error: Uncertain is wrong type");
  int uncert = *LOGICAL(Uncertain);

  SEXP Result;
  PROTECT(Result = allocMatrix(REALSXP, M, P));
  double *result = REAL(Result);

  for (int j=0, ij=0, jks=0; j<P; j++, jks+=N) {
    for (int i=0, ik=0; i<M; i++, ij++) {
      double sg=0.0, sgg=0, sx=0.0, sxx=0.0, sgx=0.0;
      int s=0;
      for (int k=0, jk=jks; k<N; k++) {
	unsigned char g = snps[ik++];
	double xk = x[jk++];
	if (g && ((g<4)|uncert) && !ISNA(xk)) {
	  double gm = g2mean(g);
	  s++;
	  sg += gm;
	  sgg += gm*gm;
	  sx += xk;
	  sxx += xk*xk;
	  sgx += gm*xk;
	}
      }
      if (s) {
	sgg -= sg*sg/(double)s;
	sxx -= sx*sx/(double)s;
	sgx -= sg*sx/(double)s;
	if (sgg>0.0 && sxx>0.0)
	  result[ij] = sgx/sqrt(sgg*sxx);
	else
	  result[ij] = NA_REAL;
      }
      else {
	result[ij] = NA_REAL;
      }
    }
  }
  
  UNPROTECT(1);
  return(Result);
}

/*
  IBS matrix.

  Calculates NxN matrix of type integer, whose upper triangle counts number 
  non-missing pairs of chromosomes and whose lower triangle counts number 
  IBS. The diagonal counts non-missing calls for each subject

  For autosomes, each pair of non-missing SNPs add 4 above the diagonal, and 
  0, 2, or 4 below the diagonal, according to IBS state:
                      AA AB BB
		   AA  4  2  0
                   AB  2  2  2
                   BB  0  2  4
  For an X SNP, pairs of females count similarly but pairs of males count 1 
  above the diagonal and a male/femal pair counts 2. Below the diagonal we 
  have:
       AY BY               AA AB BB
    AY  1  0      or    AY  2  1  0
    BY  0  1            BY  0  1  2

  If Uncertain==TRUE, contributions are weighted by posterior probability.
  Otherwise they are treated as missing.

*/


SEXP ibs_count(const SEXP Snps, const SEXP Uncertain) { 

  double lutab[3][3] = {{4., 2., 0.}, {2., 2., 2.}, {0., 2., 4.}};

  int *ifFemale = NULL;
  SEXP cl = GET_CLASS(Snps);
  if (TYPEOF(cl) != STRSXP) {
    cl = R_data_class(Snps, FALSE); /* S4 way of getting class attribute */
  }
  if (!strcmp(CHAR(STRING_ELT(cl, 0)), "XSnpMatrix")) {
    SEXP Female = R_do_slot(Snps, mkString("Female"));
    if (TYPEOF(Female)!=LGLSXP)
      error("Argument error -  Female slot wrong type");
    ifFemale = LOGICAL(Female);
  }
  else if (strcmp(CHAR(STRING_ELT(cl, 0)), "SnpMatrix")) {
    error("Argument error - Snps wrong type");
  }    
  SEXP names = getAttrib(Snps, R_DimNamesSymbol);
  if (names == R_NilValue) {
    error("Argument error - Snps object has no names");
  }
  SEXP rowNames = VECTOR_ELT(names, 0);
  if (rowNames == R_NilValue) {
    error("Argument error - Snps object has no row names");
  }

  const unsigned char *snps = RAW(Snps);
  int N, M;
  int *dim = INTEGER(getAttrib(Snps, R_DimSymbol));
  N = dim[0];
  M = dim[1];

  /* Handling of uncertain genotypes */

  if (TYPEOF(Uncertain) != LGLSXP)
    error("Argument error: Uncertain is wrong type");
  int uncert = *LOGICAL(Uncertain);

  /* Result matrix */

  SEXP Result, dimnames;
  PROTECT(Result = allocMatrix(REALSXP, N, N));
  PROTECT(dimnames = allocVector(VECSXP, 2));
  SET_VECTOR_ELT(dimnames, 0, duplicate(rowNames));    
  SET_VECTOR_ELT(dimnames, 1, duplicate(rowNames));    
  setAttrib(Result, R_DimNamesSymbol, dimnames);    
  double *result = REAL(Result);
  memset(result, 0x00, N*N*sizeof(double));

  /* Update result matrix for each locus in turn */

  for (int ik=0, k=0; k<M; k++) {

    /* Update IBS matrix */

    int N1 = N+1;
    for (int i=0, ii=0; i<N; i++, ii+=N1) {
      int base_div;
      if (ifFemale && !ifFemale[i])
	base_div = 2;
      else
	base_div = 1;
      unsigned char sik = snps[ik++];
      if (sik&&((sik<4)|uncert)) {
	result[ii]++;
	double pi[3];
	g2post(sik, pi, pi+1, pi+2);
	for (int j=i+1, jk=ik, ji=ii+1, ij=ii+N; 
	     j<N; j++, ji++, ij+=N) {
	  int div = base_div ;
	  if (ifFemale && !ifFemale[j])
	    div *= 2;
	  unsigned char sjk = snps[jk++]; 
	  if (sjk&&((sjk<4)|uncert)) {
	    double pj[3];
	    g2post(sjk, pj, pj+1, pj+2);
	    double add = 0;
	    for (int ii=0; ii<2; ii++) {
	      double pii = pi[ii];
	      if (pii) {
		for (int jj=0; jj<3; jj++) {
		  double pjj = pj[jj];
		  if (pjj)
		    add += pii*pjj*lutab[ii][jj];
		}
	      }
	    }
	    result[ij] += add/div;
	    result[ji] += 4/div;
	  }
	}
      }
    }
  }
  
  UNPROTECT(2);
  return(Result);
}

/*
  Distance matrix based on IBS counts
*/

SEXP ibs_dist(const SEXP Ibsc) {
  
  if (!isReal(Ibsc))
    error("Input object is not a real array");

  const double *ibsc = REAL(Ibsc);
  int N, M;
  int *dim = INTEGER(getAttrib(Ibsc, R_DimSymbol));
  N = dim[0];
  M = dim[1];
  if (!N || N!=M) 
    error("Input object is not a square matrix");
  SEXP names = getAttrib(Ibsc, R_DimNamesSymbol);
  if (names == R_NilValue) {
    error("Argument error - no names");
  }
  SEXP rowNames = VECTOR_ELT(names, 0);
  if (rowNames == R_NilValue) {
    error("Argument error - no sample identifiers");
  }

  /* Result matrix */

  R_len_t Nout = (N*(N-1))/2;
  SEXP Result, Size, Class;
  PROTECT(Result = allocVector(REALSXP, Nout));
  PROTECT(Size = allocVector(INTSXP, 1));
  INTEGER(Size)[0] = N;
  PROTECT(Class = allocVector(STRSXP, 1));
  SET_STRING_ELT(Class, 0, mkChar("dist"));
  setAttrib(Result, install("Size"), Size);
  setAttrib(Result, install("Labels"), duplicate(rowNames));
  classgets(Result, Class);
  double *result = REAL(Result);
  memset(result, 0x00, Nout*sizeof(double));
  for (int i=1, ii=0, k=0; i<N; i++, ii+=(N+1)) {
    for (int j=i, ji=ii+1, ij=ii+N; j<N; j++, ji++, ij+=N){
      result[k++] = (ibsc[ji]-ibsc[ij])/ibsc[ji];
    }
  }

  UNPROTECT(3);
  return(Result);
}


