#include <R.h>
#include <Rinternals.h>
#include <stdio.h>
#include <string.h>

#include "Rmissing.h"

/* #include "glm_test.h" */
#include "hash_index.h"
#include "imputation.h"
#include "invert.h"
#include "uncertain.h"

#define MAX_NAME 128
              
/* Routines for SSq and products. In X.Y case, Y variable varies fastest 
   in result array. That is, ssp is an My*Mx matrix */

/* Complete lines only */

void ssqprod_c(const int N, const int Mx, const double *X,
	       const int My, const double *Y, 
	       const int *stratum, const int *order,
	       double *ssp, int *df) {
  double *sumx = (double *)Calloc(Mx, double);
  memset(sumx, 0x00, Mx*sizeof(double));
  double *sumy = NULL;
  int nssp;
  if (My) {
    sumy = (double *)Calloc(My, double);
    memset(sumy, 0x00, My*sizeof(double));
    nssp = Mx*My;
  }
  else {
    nssp = (Mx*(Mx+1))/2;
  }
  memset(ssp, 0x00, nssp*sizeof(double));
  *df = 0;
  int Ns = 0;
  int sprev=NA_INTEGER;
  for (int s=0; s<N; s++) {
    int is = order[s]-1;
    if (is<0){  /* order[i]==0 means skip this row */
      continue;
    }
    if (stratum && (stratum[is]!=sprev) && Ns) {
      sprev = stratum[is];
      for (int i=0, ij=0; i<Mx; i++) {
	double sxi = sumx[i];
	if (My) {
	  for (int j=0; j<My; j++, ij++)
	    ssp[ij] -= (sxi*sumy[j]/(double)Ns);
	}
	else {
	  for (int j=0; j<=i; j++, ij++)
	    ssp[ij] -= (sxi*sumx[j]/(double)Ns);
	}
      }
      *df += (Ns - 1);
      Ns = 0;
      memset(sumx, 0x00, Mx*sizeof(double));
      if (My)
	memset(sumy, 0x00, My*sizeof(double));
    }
    Ns++;
    const double *Xi = X+is;
    for (int i=0, ij=0; i<Mx; i++, Xi+=N) {
      sumx[i] += (*Xi);
      if (My) {
	const double *Yj = Y+is;
	for (int j=0; j<My; j++, ij++, Yj+=N) {
	  if (!i)
	    sumy[j] += (*Yj);
	  ssp[ij] += (*Xi)*(*Yj);
	}
      }
      else {
	const double *Xj = X+is;
	for (int j=0; j<=i; j++, ij++, Xj+=N)
	  ssp[ij] += (*Xi)*(*Xj);
      }
    }
  }
  if (Ns) {
    for (int i=0, ij=0; i<Mx; i++) {
      double sxi = sumx[i];
      if (My) {
	for (int j=0; j<My; j++, ij++)
	  ssp[ij] -= (sxi*sumy[j]/(double)Ns); 
      }
      else {
	for (int j=0; j<=i; j++, ij++)
	  ssp[ij] -= (sxi*sumx[j]/(double)Ns);
      }
    }
    *df += (double)(Ns - 1); 
  }
  Free(sumx);
  if (My)
    Free(sumy);
}	  

/* Incomplete lines */

void ssqprod_i(const int N, const int Mx, const double *X,
	       const int My, const double *Y, 
	       const int *stratum, const int *order,
	       double *ssp, int *df) {
  const double *Xi = X;
  for (int i=0, ij=0; i<Mx; i++, Xi+=N) {
    const double *Yj;
    int jto;
    if (My) {
      Yj = Y;
      jto = My;
    }
    else {
      Yj = X;
      jto = i+1;
    }
    int sprev = NA_INTEGER;
    for (int j=0; j<jto; j++, ij++, Yj+=N) {
      double wxy = 0.0, wx = 0.0, wy = 0.0;
      int wdf = 0;
      int Ns = 0;
      for (int s=0; s<N; s++) {
	int is = order[s]-1;
	if (is<0) /* order[i]==0 means skip this row */
	  continue;
	if (stratum && (stratum[is]!=sprev)) {
	  sprev = stratum[is];
	  wxy -= (wx*wy/(double)Ns);
	  wx = wy = 0.0;
	  wdf += (Ns - 1);
	  Ns = 0;
	}
	double xis = Xi[is];
	double yjs = Yj[is];
	if (!(ISNA(xis)) || (ISNA(yjs))){
	  Ns++;
	  wx += xis;
	  wy += yjs;
	  wxy += xis*yjs;
	}
      }
      wxy -= (wx*wy/(double)Ns);      
      wdf += (Ns - 1);
      ssp[ij] = wxy;
      df[ij] = wdf;
    }
  }
}

/* Multiplier for var-cov matrix of covariances between X and Y */

void mcov_i(const int N, const int Mx, const double *X,
	       const int My, const double *Y, 
	       const int *stratum, const int *order,
	       double *multiplier) {
}

	    

SEXP mvphen(const SEXP Pheno, const SEXP Snps, const SEXP Rules, 
	    const SEXP Stratum, const SEXP Order, const SEXP Tests, 
	    const SEXP Complete, const SEXP Uncertain, const SEXP Score){
  
  Rprintf("**** Warning: this code is under development\n");
  /* Pheno should be a matrix or a vector */

  if (TYPEOF(Pheno) != REALSXP)
    error("illegal storage mode for phenotype");
  double *Y = NULL;
  int N = 0, M = 0;
  SEXP Pnames = R_NilValue;
  if (isMatrix(Pheno)) {
    N = nrows(Pheno);
    M = ncols(Pheno);
    Y = REAL(Pheno);
    Pnames = VECTOR_ELT(getAttrib(Pheno, R_DimNamesSymbol), 1);
  }
  else {
    M = 1;
    N = length(Pheno);
  }

  /* Snp should be a SnpMatrix or an XSnpMatrix */

  const char *classSnps = CHAR(STRING_ELT(getAttrib(Snps, R_ClassSymbol), 0));
  int ifX = 0;
  int *female=NULL;
  if (!strcmp(classSnps, "SnpMatrix"))
    ifX = 0;
  else if (!strcmp(classSnps, "XSnpMatrix")) {
    ifX = 1;
    error("tests for XSnpMatrix not yet implemented");
  } 
  else 
    error("Argument error - class(Snps): %s", classSnps);
  if (nrows(Snps) != N)
    error("unequal numbers of rows in phenotype and snp matrices");
  int nsnp = ncols(Snps);
  unsigned char *snps = RAW(Snps);
  
  /* If imputation involved, calculate snp name index and lookup tables */
  
  index_db name_index = NULL;
  SEXP Snp_names =  VECTOR_ELT(getAttrib(Snps, R_DimNamesSymbol), 1);
  SEXP Rule_names = R_NilValue;	
  int pmax = 0;
  GTYPE **gt2ht = NULL;
  if (TYPEOF(Rules)!=NILSXP) {
    name_index = create_name_index(Snp_names);
    Rule_names = getAttrib(Rules, R_NamesSymbol);
    pmax = *INTEGER(getAttrib(Rules, install("Max.predictors")));
    gt2ht = (GTYPE **)Calloc(pmax, GTYPE *);
    for (int i=0; i<pmax; i++)
      gt2ht[i] = create_gtype_table(i+1);
  }
  
  /* Tests. Snps referred to by col number not name */

  int test_size = 1, test_list = 0;
  int *test_int = NULL;
  int ntest = nsnp;
  SEXPTYPE tests_type = TYPEOF(Tests);
  if (tests_type==VECSXP) {
    ntest = LENGTH(Tests);
    test_list = 1;
    for (int i=0; i<ntest; i++) {
      SEXP testi =  VECTOR_ELT(Tests, i);
      if (TYPEOF(testi)!=INTSXP)
	error("Non-integer list - Tests[i]");
      int leni = LENGTH(testi);
      if (leni>test_size)
	test_size = leni;
    }
    if (test_size==1)
      error("List of multi-SNP tests doesn't contain a multi-SNP spec");
  }
  else if (tests_type==INTSXP) {
    test_int = INTEGER(Tests);
    ntest = LENGTH(Tests);
  }
  else if (tests_type!=NILSXP)
    error("Argument error - sets");
  
  /* Stratum */
  
  int *stratum = NULL;
  if (TYPEOF(Stratum)==INTSXP) {
    if (LENGTH(Stratum)!=N) 
      error("Dimension error - Stratum"); 
    stratum = INTEGER(Stratum);
  }
  else if (TYPEOF(Stratum)!=NILSXP)
    error("Argument error - Stratum not integer");
  
  /* Order */

  int *order = NULL;
  int *lorder = Calloc(N, int);
  if (TYPEOF(Order)==INTSXP) {
    order = INTEGER(Order);
    for (int i=0; i<N; i++)
      lorder[i] = order[i];
  }
  else if (TYPEOF(Order)==NILSXP) {
    for (int i=0; i<N; i++)
      lorder[i] = i+1;
  }
  else 
    error("Argument error - Order argument not integer");

    
  /* Complete lines only? */

  if (TYPEOF(Complete)!=LGLSXP) 
    error("illegal complete argument");
  int complete = *LOGICAL(Complete);

  /* Use uncertain genotypes? */
  
  if (TYPEOF(Uncertain)!=LGLSXP) 
    error("illegal uncertain argument");
  int uncertain = *LOGICAL(Uncertain);
  
  /* Save U and V? */
  
  if (TYPEOF(Score)!=LGLSXP) 
    error("illegal score argument");
  int score = *LOGICAL(Score);

  /* Output list */
  
  SEXP Result, TestNames = R_NilValue, Chisq, Df, Nused, Namelist; 
  PROTECT(Result = allocS4Object());
  PROTECT(Chisq = allocVector(REALSXP, ntest));
  double *chisq = REAL(Chisq);
  PROTECT(Df = allocVector(INTSXP, ntest));
  int *df = INTEGER(Df);
  PROTECT(Nused = allocVector(INTSXP, ntest));
  int *nused = INTEGER(Nused);
  if (test_list)
    PROTECT(Namelist = allocVector(VECSXP, ntest));
  else
    PROTECT(Namelist = allocVector(STRSXP, ntest));
  R_do_slot_assign(Result, mkString("chisq"), Chisq);
  R_do_slot_assign(Result, mkString("df"), Df);
  R_do_slot_assign(Result, mkString("N"), Nused);
  R_do_slot_assign(Result, mkString("snp.names"), Namelist);
  R_do_slot_assign(Result, mkString("var.names"), Pnames);
		   
  /* If score and score variance to be saved */

  SEXP Save_score = R_NilValue, UVnames = R_NilValue;
  if (score) {
    PROTECT(Save_score = allocVector(VECSXP, ntest));
    R_do_slot_assign(Result, mkString("score"), Save_score);
    PROTECT(UVnames = allocVector(STRSXP, 2));
    SET_STRING_ELT(UVnames, 0, mkChar("U"));
    SET_STRING_ELT(UVnames, 1, mkChar("V"));
  }
		   
  /* Test names */

  if (test_list) {
    SEXP NameTests = getAttrib(Tests, R_NamesSymbol);
    setAttrib(Namelist, R_NamesSymbol, NameTests);
  }

  
  /* Do tests */

  int st;
  int *sts = &st;

   /* Require Y rows complete -- regardless of Complete argument */

  int *xmiss = (int *)Calloc(N, int);
  memset(xmiss, 0x00, N*sizeof(int));
  for (int i=0; i<N; i++) {
    for (int j=0, ij=i; j<M; j++, ij+=N) {
      if (ISNA(Y[ij])) {
	xmiss[i] = 1;
	break;
      }
    }
  }
  int Ncomplete = 0;
  for (int i=0; i<N; i++) {
    int is = lorder[i] -1;
    if (xmiss[is])
      lorder[i] = 0;
    else
      Ncomplete++;
  }
		  	   
  for (int test=0; test<ntest; test++) {
    SEXP Snames = R_NilValue;
    int nsts = 1;
    if (test_list) { /* List of multi-SNP tests */
      SEXP St = VECTOR_ELT(Tests, test);
      nsts = LENGTH(St);
      sts = INTEGER(St);
      Snames = allocVector(STRSXP, nsts);
      SET_VECTOR_ELT(Namelist, test, Snames); /* Should protect Snames */
    }
    else { /* Single snps */
      if (test_int)  /* Selected list */
	sts = test_int+test;
      else
	st = test + 1;
    }

    /* Arrays to hold score and score variance */
  
    int nu = nsts*M;
    int nv = (nu*(nu+1)/2);
    SEXP U = R_NilValue, V = R_NilValue;
    double *u, *v;
    if (score) {
      if (M == 1) 
	PROTECT(U = allocVector(REALSXP, nsts));
      else
	PROTECT(U = allocMatrix(REALSXP, M, nsts));
      PROTECT(V = allocVector(REALSXP, nv));
      u = REAL(U);
      v = REAL(V);
    }
    else {
      u = (double *) Calloc(nu, double);
      v = (double *) Calloc(nv, double);
    }
 
    int *lorder_t = lorder;
    if (complete) {
      lorder_t = (int *) Calloc(N, int);
      memset(lorder_t, 0x00, N*sizeof(int));
    }

    /* Extract X matrix */
    
    int space = MAX_NAME-1;
    char testname[MAX_NAME];
    double *X = (double *)Calloc(N*nsts, double);
    memset(xmiss, 0x00, N*sizeof(int));
    for (int j=0, ij=0; j<nsts; j++) {
      int stsj = sts[j];
      if (stsj>0) {
	stsj--;
	SEXP Snp_namej =  STRING_ELT(Snp_names, stsj);
	if (test_list)
	  SET_STRING_ELT(Snames, j, Snp_namej);	  
	else
	  SET_STRING_ELT(Namelist, test, Snp_namej);
	const unsigned char *gj = snps + N*stsj;
	for (int i=0; i<N; i++, ij++) {
	  unsigned char gij = gj[i];
	  if (gij && (uncertain||(gij<4))) {
	    X[ij] = g2mean(gij);
	  }
	  else {
	    X[ij] = NA_REAL;
	    if (complete)
	      xmiss[i] = 1;
	  }
	}
      }
      else {
	stsj = -(1+stsj);
	SEXP Rule_namej = STRING_ELT(Rule_names, stsj);
	if (test_list)
	  SET_STRING_ELT(Snames, j, Rule_namej);	  
	else
	  SET_STRING_ELT(Namelist, test, Rule_namej);
	SEXP Rule =  VECTOR_ELT(Rules, stsj);
	if (!isNull(Rule)){ /* Not monomorphic */
	  do_impute(Snps, N, female, NULL, N, name_index, Rule, gt2ht, X+ij, 
		    NULL);
	  ij += N;
	}
	else {
	  X[ij++] = NA_REAL;
	  if (complete) {
	    for (int i=0; i<N; i++)
	      xmiss[i] = 1;
	  }
	}   
      }
    }
    if (complete) {
      Ncomplete = 0;
      for (int i=0; i<N; i++) {
	int is = lorder[i] - 1;
	if (is<0 || xmiss[is])
	  lorder_t[i] = 0;
	else {
	  lorder_t[i] = order[i];
	  Ncomplete++;
	}
      }
    }
    if (!Ncomplete) {
      warning("No data -- test skipped");
      if (!score) {
	Free(u);
	Free(v);
      }
      if (complete)
	Free(lorder_t);
      Free(X);
    }
    else {

      /* Do calculations */

      int nsts2 = (nsts*(nsts+1))/2;
      double *XX = (double *)Calloc(nsts2, double);
      int M2 = (M*(M+1))/2;
      double *YY = (double *)Calloc(M2, double);
      int df_r;
      int *df_rXX = NULL, *df_rYY = NULL, *df_rXY = NULL;
      double *df_rV; 
      if (complete) {
	ssqprod_c(N, nsts, X, 0, NULL, stratum, lorder_t, XX, &df_r); 
	ssqprod_c(N, M, Y, 0, NULL, stratum, lorder_t, YY, &df_r);
	ssqprod_c(N, nsts, X, M, Y, stratum, lorder_t, u, &df_r);
      }
      else {
	error("Incomplete data option not yet implemented");
	df_rXX = (int *)Calloc(nsts2, int);
	df_rYY = (int *)Calloc(M2, int);
	df_rXY = (int *)Calloc(nu, int);
	df_rV = NULL; /* To be corrected */
	ssqprod_i(N, nsts, X, 0, NULL, stratum, lorder, XX, df_rXX); 
	ssqprod_i(N, M, Y, 0, NULL, stratum, lorder, YY, df_rYY);
	ssqprod_i(N, nsts, X, M, Y, stratum, lorder, u, df_rXY);
      }
      Free(X); 
    
      /* Make v matrix  */
  
      for (int i=0, ips=0, iv=0; i<nsts; i++, ips+=i) {
	for (int j=0, jqs=0; j<M; j++, jqs+=j) {
	  for (int p=0, ip=ips; p<=i; p++, ip++) {
	    double XXip=XX[ip];
	    int qto = (p==i)? j+1: M;
	    for (int q=0, jq=jqs; q<qto; q++, jq+=(q<j)?1:q, iv++) {	      
	      v[iv] = YY[jq]*XXip/(double)df_r;
	    }
	  }
	}
      }

      /* Calculate chi-squared and degrees of freedom */

      double chi2 = NA_REAL, ifault;
      int dft = NA_INTEGER;
      
      if (score) {
	double *w = (double *)Calloc(nv, double);
	ifault = qform(nu, u, v, w, &chi2, &dft);
	Free(w);
      }
      else {
	ifault = qform(nu, u, v, v, &chi2, &dft);
      }
     
      if (ifault) {
	warning("Variance matrix not positive semi-definite");
	chi2 = NA_REAL;
	dft = NA_INTEGER;
      }

      chisq[test] = chi2;
      df[test] = dft;
      nused[test] = Ncomplete;

      /* Tidy up at end of this test */


      if (score) {
	SEXP Scoret, Udnames;
	PROTECT(Scoret = allocVector(VECSXP, 2));
	setAttrib(Scoret, R_NamesSymbol, UVnames);
	SET_VECTOR_ELT(Scoret, 0, U);
	SET_VECTOR_ELT(Scoret, 1, V);
	SET_VECTOR_ELT(Save_score, test, Scoret);
	UNPROTECT(3);
      }
      else {
	Free(u);
	Free(v);
      }
      Free(XX);
      Free(YY);
      if(complete) {
	Free(lorder_t);
      }
      else {
	Free(df_rXX);
	Free(df_rYY);
	Free(df_rXY);
	Free(df_rV);
      }
    }
  }

  /* Return hash table memory and gt->ht tables */

  Free(lorder);
  if (name_index) {
    index_destroy(name_index);
    for (int i=0; i<pmax; i++)
      destroy_gtype_table(gt2ht[i], i+1);
    Free(gt2ht);
  }
  Free(lorder);
  Free(xmiss);

  SEXP Class, Package;
  PROTECT(Class = allocVector(STRSXP, 1));
  SET_STRING_ELT(Class, 0, mkChar(score?  "GlmtestsScore": "GlmTests"));
  PROTECT(Package = allocVector(STRSXP, 1));
  SET_STRING_ELT(Package, 0, mkChar("snpStats"));
  setAttrib(Class, install("package"), Package);
  classgets(Result, Class);

  if (score)
    UNPROTECT(9); /* Save_score, UVnames PLUS ... */
  else
    UNPROTECT(7); /* Package, Class, Nused, Df, Chisq, Namelist, Result */
  return(Result);
}
