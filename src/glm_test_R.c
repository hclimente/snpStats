#include <R.h>
#include <Rinternals.h>
#include <stdio.h>
#include <string.h>

#include "glm_test.h"
#include "hash_index.h"
#include "imputation.h"
#include "invert.h"
#include "uncertain.h"
#include "Rmissing.h"

#define MAX_NAME_LENGTH 81

SEXP snp_lhs_score(const SEXP Y, const SEXP X, const SEXP Stratum, 
		   const SEXP Z, const SEXP Snp_subset, 
		   const SEXP Robust, const SEXP Cluster, const SEXP Uncertain,
		   const SEXP Control, const SEXP If_score) {

  /* Y should be a SnpMatrix or an XSnpMatrix */
  const char *classY = NULL;
  if (TYPEOF(R_data_class(Y, FALSE)) == STRSXP) {
    classY = CHAR(STRING_ELT(R_data_class(Y, FALSE), 0));
  } else {
    classY = CHAR(STRING_ELT(getAttrib(Y, R_ClassSymbol), 0));
  }
  if(!IS_S4_OBJECT(Y)) {
    error("Y in snp_lhs_score is missing S4 Object bit");
  }
  int ifX = 0;
  if (!strcmp(classY, "SnpMatrix"))
    ifX = 0;
  else if (!strcmp(classY, "XSnpMatrix"))
    ifX = 1;
  else 
    error("Argument error - class(Y)");

  /* Y and its dimensions */

  if (TYPEOF(Y)!=RAWSXP)
    error("Argument error - Y");
  const unsigned char *y = RAW(Y);
  int *dim;
  int N, nsnp;
  if (strlen(classY)>5) {
    dim = INTEGER(getAttrib(Y, R_DimSymbol));
    N = dim[0];
    nsnp = dim[1];
  }
  else {
    N = LENGTH(Y);
    nsnp = 1;
  }

  /* SNP subset */

  int *snp_subset = NULL;
  int ntest = nsnp;
  SEXPTYPE sstype = TYPEOF(Snp_subset);
  if (sstype==INTSXP) {
    if (LENGTH(Snp_subset)>nsnp)
      error("Dimension error - Snp_subset");
    snp_subset = INTEGER(Snp_subset);
    ntest = LENGTH(Snp_subset);
  }
  else if (sstype!=NILSXP)
    error("Argument error - Snp.subset");
  

  /* If X-chromosome, indicators of diploid genotype */

  const int *diploid;
  if (ifX) {
    SEXP Diploid = R_do_slot(Y, mkString("diploid")); 
    diploid = LOGICAL(Diploid);
  }
  else
    diploid = NULL;
  
  /* X and its dimensions */

  double *x = NULL;
  int M = 0;
  if (TYPEOF(X)==REALSXP) {
    x = REAL(X);
    const int *dim = INTEGER(getAttrib(X, R_DimSymbol));
    M = dim[1];
    if (dim[0]!=N)
      error("Dimension error - X"); 
  }
  else if (TYPEOF(X)!=NILSXP)
    error("Argument error - X");


  /* Stratum */

  int S = 1;
  const int *stratum = NULL;
  if (TYPEOF(Stratum)==INTSXP) {
    if (LENGTH(Stratum)!=N) 
      error("Dimension error - Stratum"); 
    stratum = INTEGER(Stratum);
    for (int i=0; i<N; i++) 
      if (stratum[i]>S) S = stratum[i];
  }
  else if (TYPEOF(Stratum)!=NILSXP)
    error("Argument error - Stratum");

  /* Z and its dimensions */

  if (TYPEOF(Z)!=REALSXP)
    error("Argument error - Z");
  double *z = REAL(Z);
  dim = INTEGER(getAttrib(Z, R_DimSymbol));
  int P = dim[1];
  if (dim[0]!=N)
    error("Dimension error - Z"); 


  /* Robust */

  if (TYPEOF(Robust)!=LGLSXP)
    error("Argument error - Robust");
  int robust = LOGICAL(Robust)[0];

  /* Cluster */

  int C = robust? 1: 0;
  int *cluster = NULL;
  if (TYPEOF(Cluster)==INTSXP) {
    if (LENGTH(Cluster)!=N)
      error("Dimension error - Cluster"); 
    cluster = INTEGER(Cluster);
    for (int i=0; i<N; i++)
      if (cluster[i]>C) C = cluster[i];
  }
  else if (TYPEOF(Cluster)!=NILSXP)
    error("Argument error - Cluster");


  /* Handling of uncertain genotypes */

  if (TYPEOF(Uncertain) != LGLSXP)
    error("Argument error: Uncertain is wrong type");
  int uncert = *LOGICAL(Uncertain);
  /* Force robust variance if uncertain option used */
  robust = robust || uncert; 

  
  /* Control object */

  if (TYPEOF(Control)!=VECSXP || LENGTH(Control)!=3) 
    error("Argument error - Control");
  SEXP Maxit = VECTOR_ELT(Control, 0);
  if (TYPEOF(Maxit)!=INTSXP || length(Maxit)!=1)
    error("Argument error - Maxit");
  int maxit = INTEGER(Maxit)[0];
  SEXP Epsilon = VECTOR_ELT(Control, 1);
  if (TYPEOF(Epsilon)!=REALSXP || LENGTH(Epsilon)!=1)
    error("Argument error - Epsilon");
  double epsilon = REAL(Epsilon)[0];
  SEXP R2Max = VECTOR_ELT(Control, 2);
  if (TYPEOF(R2Max)!=REALSXP || LENGTH(R2Max)!=1)
    error("Argument error - R2Max");
  double r2Max = REAL(R2Max)[0];
  
  /* Are scores and score variances required? */

  int if_score = *LOGICAL(If_score);

  /* Work arrays */

  double *yd = (double *) R_alloc(N, sizeof(double));
  double *fitted = (double *) R_alloc(N, sizeof(double));
  double *prior = (double *) R_alloc(N, sizeof(double));
  double *resid = (double *) R_alloc(N, sizeof(double));
  double *weights = (double *) R_alloc(N, sizeof(double));
  double *xb;
  if (x)
    xb = (double *) R_alloc(N*M, sizeof(double));
  else 
    xb = NULL;

  /* Output arrays */

  SEXP Result, Rnames, Vnames, Chisq, Df, Nused, 
    Score = R_NilValue, UVnames = R_NilValue, Varnames = R_NilValue;
  PROTECT(Result = allocS4Object());

  SEXP snpNames = VECTOR_ELT(getAttrib(Y, R_DimNamesSymbol), 1);
  SEXPTYPE sntype = TYPEOF(snpNames);
  if (sntype!=STRSXP)
    error("SNP names not string variables");
  PROTECT(Rnames = allocVector(STRSXP, ntest));
  PROTECT(Vnames = VECTOR_ELT(getAttrib(Z, R_DimNamesSymbol), 1));
  PROTECT(Chisq = allocVector(REALSXP, ntest));
  double *chisq = REAL(Chisq);
  PROTECT(Df = allocVector(INTSXP, ntest));
  int *df = INTEGER(Df);
  PROTECT(Nused = allocVector(INTSXP, ntest));
  int *nused = INTEGER(Nused);
  PROTECT(Varnames = allocVector(STRSXP, 1));
  /* Should put something in as a place holder */
  SET_STRING_ELT(Varnames, 0, mkChar("X-variable"));
  R_do_slot_assign(Result, mkString("var.names"), Varnames);
  R_do_slot_assign(Result, mkString("snp.names"), Rnames);
  R_do_slot_assign(Result, mkString("chisq"), Chisq);
  R_do_slot_assign(Result, mkString("df"), Df);
  R_do_slot_assign(Result, mkString("N"), Nused);

  if (if_score) {
    PROTECT(Score = allocVector(VECSXP, ntest));
    R_do_slot_assign(Result, mkString("score"), Score);
    PROTECT(UVnames = allocVector(STRSXP, 2));
    SET_STRING_ELT(UVnames, 0, mkChar("U"));
    SET_STRING_ELT(UVnames, 1, mkChar("V"));
  }

  /* Do calculations */

  for (int t=0; t<ntest; t++) {

    SEXP U = R_NilValue, V = R_NilValue;
    double *u, *v;
    if (if_score) {
      PROTECT(U = allocVector(REALSXP, P));
      PROTECT(V = allocVector(REALSXP, (P*(P+1))/2));
      u = REAL(U);
      v = REAL(V);
    }
    else {
      u = (double *) Calloc(P, double);
      v = (double *) Calloc((P*(P+1))/2, double);
    }

    int j = snp_subset? snp_subset[t] - 1: t; 
    const unsigned char *yj = y + N*j;
    int mono = 1, yv = 0, nu = 0;
 

    /* Load SNP as Binomial y-variate, with prior weights */

    for (int i=0; i<N; i++) {
      int yij = (int) yj[i];
      if (yij && (uncert||yij<4)) {
	if (mono) {
	  if (yij<4) {
	    mono = !yv || (yij==yv);
	    yv = yij;
	  }
	  else 
	    mono = 0;
	}
	prior[i] = (!diploid || diploid[i])? 2.0: 1.0;
	yd[i] = g2mean(yij)/2.0;
	nu++;
      }
      else {
	prior[i] = yd[i] = 0.0;
      }
    }
    nused[t] = nu;

    if (mono) { /* Monomorphic SNP */
      warning("Monomorphic SNP: %d", t+1);
      memset(u, 0x00, P*sizeof(double));
      memset(v, 0x00, sizeof(double)*(P*(P+1))/2); 
    }
    else {
      
      /* Fit base model */
    
      int rank, dfr;
      double scale;
      double rc = glm_fit(BINOMIAL, LOGIT, N, M, 0, S, yd, prior, x, stratum, 
			  maxit, epsilon, r2Max, 0, 
			  &rank, xb, fitted, resid, weights, &scale, &dfr,
			  NULL, NULL, NULL, NULL);
      if (rc) 
	warning("Failure to converge while fitting base model for SNP %d",j+1);
    
      /* Score test */
      
      if (dfr) {
	glm_score_test(N, rank, S, stratum, P, z, C, cluster,
		       resid, weights, xb, scale,
		       r2Max, u, v);
      }
      else {
	memset(u, 0x00, P*sizeof(double));
	memset(v, 0x00, sizeof(double)*(P*(P+1))/2);
      }
    }
    int nunit = 0;
    for (int i=0; i<N; i++) 
      if (weights[i]) nunit++;
    nused[t] = nunit;
        
    if (P>1) {
      if (qform(P, u, v, NULL, chisq+t, df+t)) {
	warning("Matrix not positive semi-definite in test ", t+1);
	chisq[t] = NA_REAL;
	df[t] = NA_INTEGER;
      }
    }
    else {
      if ((*v)>0) {
	chisq[t] = (*u)*(*u)/(*v);
	df[t] = 1;
      }
      else {
	chisq[t] = NA_REAL;
	df[t] = NA_INTEGER;
      }
    }
    if (if_score) {
      SEXP Scoret;
      PROTECT(Scoret = allocVector(VECSXP, 2));
      setAttrib(Scoret, R_NamesSymbol, UVnames);
      SET_VECTOR_ELT(Scoret, 0, U);
      SET_VECTOR_ELT(Scoret, 1, V);
      SET_VECTOR_ELT(Score, t, Scoret);
      UNPROTECT(3);
    }
    
    /* Name of test */

    SET_STRING_ELT(Rnames, t, mkChar(CHAR(STRING_ELT(snpNames,j))));

    /* Return work space */

    if (!if_score) {
      Free(u);
      Free(v);
    }
  }
  
  SEXP Class, Package;
  PROTECT(Class = allocVector(STRSXP, 1));
  if (if_score)
    SET_STRING_ELT(Class, 0, mkChar("GlmTestsScore"));
  else 
    SET_STRING_ELT(Class, 0, mkChar("GlmTests"));
  PROTECT(Package = allocVector(STRSXP, 1));
  SET_STRING_ELT(Package, 0, mkChar("snpStats"));
  setAttrib(Class, install("package"), Package);
  classgets(Result, Class);

  if (if_score)
    UNPROTECT(11);
  else
    UNPROTECT(9);

  return(Result);
}


SEXP snp_rhs_score(SEXP Y, SEXP family, SEXP link, 
		   SEXP X, SEXP Stratum, SEXP Z, SEXP Rules, 
		   SEXP Prior, SEXP Tests, SEXP Robust, SEXP Cluster, 
		   SEXP Uncertain, 
		   SEXP Control, SEXP MissAllow, SEXP If_score) {

  char testname[MAX_NAME_LENGTH];
  int max_name_length =  MAX_NAME_LENGTH -1;
  testname[max_name_length] = (char) 0; /* Make sure its null terminated */

  /* Y and its dimensions */

  if (TYPEOF(Y)!=REALSXP)
    error("Argument error - Y");
  int N = LENGTH(Y);
  double *y = REAL(Y);

  /* Family */

  if (TYPEOF(family)!=INTSXP || LENGTH(family)!=1)
    error("Argument error - family");
  int fam = INTEGER(family)[0];
  if (fam<1 || fam>4)
    error("Illegal family argument");
  
  /* Link */

  if (TYPEOF(link)!=INTSXP || LENGTH(link)!=1)
    error("Argument error - link");
  int lnk = INTEGER(link)[0];
  if (lnk<1 || lnk>4)
    error("Illegal link argument");

  /* X and its dimensions */

  double *x = NULL;
  int M = 0;
  int *dim = NULL;
  if (TYPEOF(X)==REALSXP) {
    x = REAL(X);
    dim = INTEGER(getAttrib(X, R_DimSymbol));
    M = dim[1];
    if (dim[0]!=N)
      error("Dimension error - X"); 
  }
  else if (TYPEOF(X)==NILSXP) {
    M = 0; 
    x = NULL;
  }
  else
    error("Argument error - X");


  /* Stratum */

  int S = 1;
  int *stratum = NULL;
  if (TYPEOF(Stratum)==INTSXP) {
    if (LENGTH(Stratum)!=N) 
      error("Dimension error - Stratum"); 
    stratum = INTEGER(Stratum);
    for (int i=0; i<N; i++) 
      if (stratum[i]>S) S = stratum[i];
  }
  else if (TYPEOF(Stratum)!=NILSXP)
    error("Argument error - Stratum");

  /* Z should be a SnpMatrix or an XSnpMatrix */

  const char *classZ = NULL;
  if (TYPEOF(R_data_class(Z, FALSE)) == STRSXP) {
    classZ = CHAR(STRING_ELT(R_data_class(Z, FALSE), 0));
  } else {
    classZ = CHAR(STRING_ELT(getAttrib(Z, R_ClassSymbol), 0));
  }
  int ifX = 0;
  if (!strcmp(classZ, "SnpMatrix"))
    ifX = 0;
  else if (!strcmp(classZ, "XSnpMatrix"))
    ifX = 1;
  else 
    error("Argument error - class(Z)");
  if (!IS_S4_OBJECT(Z)) {
    error("Z in snp_rhs_score is missing S4 Object bit");
  }
 
  /* Z and its dimensions */

  int Nz;
  if (TYPEOF(Z)!=RAWSXP)
    error("Argument error - Z");
  const unsigned char *z = RAW(Z);
  int nsnp;
  if (strlen(classZ)>5) {
    dim = INTEGER(getAttrib(Z, R_DimSymbol));
    Nz = dim[0];
    nsnp = dim[1];
  }
  else {
    Nz = LENGTH(Z);
    nsnp = 1;
  }
  if (Nz!=N)
    error("Dimension error - Z");


  /* If X-chromosome, indicators of diploid genotype */

  int *diploid;
  if (ifX) {
    SEXP Diploid = R_do_slot(Z, mkString("diploid")); 
    diploid = LOGICAL(Diploid);
  }
  else
    diploid = NULL;

  /* If imputation involved, calculate snp name index and lookup tables */

  index_db name_index = NULL;
  SEXP Snp_names =  VECTOR_ELT(getAttrib(Z, R_DimNamesSymbol), 1);
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

  
  /* Prior weights */

  double *prior_all = NULL;
  if (TYPEOF(Prior)==REALSXP) {
    if (LENGTH(Prior)!=N)
      error("Dimension error - Prior");
    prior_all = REAL(Prior);
  }
  else if (TYPEOF(Prior)!=NILSXP)
    error("Argument error - Prior");

  /* Tests. For now snps referred to by col number not name */

  int test_size = 1;
  int *test_int = NULL;
  int ntest = nsnp;
  SEXPTYPE tests_type = TYPEOF(Tests);
  if (tests_type==VECSXP) {
    ntest = LENGTH(Tests);
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
    error("Argument error - tests");

  /* Robust */

  if (TYPEOF(Robust)!=LGLSXP)
    error("Argument error - Robust");
  int robust = LOGICAL(Robust)[0];

  /* Cluster */

  int C = robust? 1: 0;
  int *cluster = NULL;
  if (TYPEOF(Cluster)==INTSXP) {
    if (LENGTH(Cluster)!=N)
      error("Dimension error - Cluster"); 
    cluster = INTEGER(Cluster);
    for (int i=0; i<N; i++)
      if (cluster[i]>C) C = cluster[i];
  }
  else if (TYPEOF(Cluster)!=NILSXP)
    error("Argument error - Cluster");
  
  /* Handling of uncertain genotypes */

  if (TYPEOF(Uncertain) != LGLSXP)
    error("Argument error: Uncertain is wrong type");
  int uncert = *LOGICAL(Uncertain);

  
/* Control object */

  if (TYPEOF(Control)!=VECSXP || LENGTH(Control)!=3) 
    error("Argument error - Control");
  SEXP Maxit = VECTOR_ELT(Control, 0);
  if (TYPEOF(Maxit)!=INTSXP || length(Maxit)!=1)
    error("Argument error - Maxit");
  int maxit = INTEGER(Maxit)[0];
  SEXP Epsilon = VECTOR_ELT(Control, 1);
  if (TYPEOF(Epsilon)!=REALSXP || LENGTH(Epsilon)!=1)
    error("Argument error - Epsilon");
  double epsilon = REAL(Epsilon)[0];
  SEXP R2Max = VECTOR_ELT(Control, 2);
  if (TYPEOF(R2Max)!=REALSXP || LENGTH(R2Max)!=1)
    error("Argument error - R2Max");
  double r2Max = REAL(R2Max)[0];
  
  /* Allowable missingness */

  if (TYPEOF(MissAllow)!=REALSXP || LENGTH(MissAllow)!=1)
    error("Argument error - MissAllow");
  int mallowed = N * REAL(MissAllow)[0];

  /* Score vectors to be saved? */

  int if_score = *LOGICAL(If_score);


  /* Work arrays */

  double *zw = (double *) R_alloc(N*test_size, sizeof(double));
  double *prior = (double *) R_alloc(N, sizeof(double));
  double *resid = (double *) R_alloc(N, sizeof(double));
  double *resid_all = (double *) R_alloc(N, sizeof(double));
  double *weights = (double *) R_alloc(N, sizeof(double));
  double *weights_all = (double *) R_alloc(N, sizeof(double));
  double *fitted = (double *) R_alloc(N, sizeof(double));
  double *fitted_all = (double *) R_alloc(N, sizeof(double));
  double *xb, *xb_all;
  if (x){
    xb = (double *) R_alloc(N*M, sizeof(double));
    xb_all = (double *) R_alloc(N*M, sizeof(double));
  }
  else 
    xb = xb_all =  NULL;
  
  /* Fit base model to full data */

  int rank, rank_all, df_r, df_r_all;
  double scale, scale_all;
  int err = glm_fit(fam, lnk, N, M, 0, S, y, prior_all, x, stratum, 
		    maxit, epsilon, r2Max, 0, 
		    &rank_all, xb_all, 
		    fitted_all, resid_all, weights_all, 
		    &scale_all, &df_r_all, NULL, NULL, NULL, NULL);
  if (err) 
    error("failure to converge at initial fitting of base model");


  /* Output list */

  SEXP Result, TestNames = R_NilValue, Chisq, Df, Nused, Score = R_NilValue, 
    Varnames = R_NilValue;
  PROTECT(Result = allocS4Object());
  PROTECT(Chisq = allocVector(REALSXP, ntest));
  double *chisq = REAL(Chisq);
  PROTECT(Df = allocVector(INTSXP, ntest));
  int *df = INTEGER(Df);
  PROTECT(Nused = allocVector(INTSXP, ntest));
  int *nused = INTEGER(Nused);
  PROTECT(Varnames = allocVector(STRSXP, 1));
  /* Should put something in as a place holder */
  SET_STRING_ELT(Varnames, 0, mkChar("Y-variable"));
  if (test_size>1) 
    PROTECT(TestNames=allocVector(VECSXP, ntest));
  else {
    PROTECT(TestNames=allocVector(STRSXP, ntest));
    setAttrib(TestNames, R_NamesSymbol, getAttrib(Tests, R_NamesSymbol));
  }
  R_do_slot_assign(Result, mkString("snp.names"), TestNames);
  R_do_slot_assign(Result, mkString("var.names"), Varnames);
  R_do_slot_assign(Result, mkString("chisq"), Chisq);
  R_do_slot_assign(Result, mkString("df"), Df);
  R_do_slot_assign(Result, mkString("N"), Nused);
  int nprot = 6;  
  SEXP UVnames=R_NilValue;
  if (if_score) {
    PROTECT(Score = allocVector(VECSXP, ntest));
    R_do_slot_assign(Result, mkString("score"), Score);
    PROTECT(UVnames = allocVector(STRSXP, 2));
    SET_STRING_ELT(UVnames, 0, mkChar("U"));
    SET_STRING_ELT(UVnames, 1, mkChar("V"));
    nprot = 8;
  }
     
  /* Do tests */

  int snp;
  int *snps = &snp;
  for (int test=0; test<ntest; test++) {
    SEXP Snames=R_NilValue;
    int nsnpt = 1;
    if (test_size==1) { /* Single SNP tests */
      if (test_int)  /* Selected list */
	snps = test_int+test;
      else
	snp = test + 1;
      
    }
    else  { /* List of multi-SNP tests */
      SEXP Snps = VECTOR_ELT(Tests, test);
      nsnpt = LENGTH(Snps);
      snps = INTEGER(Snps);
      PROTECT(Snames = allocVector(STRSXP, nsnpt));
      SET_VECTOR_ELT(TestNames, test, Snames);
    }

    /* Arrays to hold score and score variance */

    SEXP U = R_NilValue, V = R_NilValue;
    double *u, *v;
    if (if_score) {
      PROTECT(U = allocVector(REALSXP, nsnpt));
      PROTECT(V = allocVector(REALSXP, (nsnpt*(nsnpt+1))/2));
      u = REAL(U);
      v = REAL(V);
    }
    else {
      u = (double *) Calloc(nsnpt, double);
      v = (double *) Calloc((nsnpt*(nsnpt+1))/2, double);
    }

    for (int i=0; i<N; i++)
      prior[i] = prior_all? prior_all[i]: 1.0;

    /* Set up Z matrix, tracking incomplete cases  */

    int missing = 0;
    for (int j=0, ij=0; j<nsnpt; j++) {
      /*      int len = 0; unused */
      int snpsj = snps[j];
      if (snpsj>0) {
	snpsj--;
	SEXP Snp_namej =  STRING_ELT(Snp_names, snpsj);
	if (test_size==1)
	  SET_STRING_ELT(TestNames, test, Snp_namej);
	else
	  SET_STRING_ELT(Snames, j, Snp_namej);
	const unsigned char *zj = z + N*snpsj;
	for (int i=0; i<N; i++, ij++) {
	  unsigned char zij = zj[i];
	  if (zij && (uncert||(zij<4))) {
	    zw[ij] = g2mean(zij);
	  }
	  else {
	    zw[ij] = 0.0;
	    if (prior[i]) {
	      missing ++;
	      prior[i] = 0.0;
	    }
	  }
	}
      }
      else {
	snpsj = -(1+snpsj);
	SEXP Rule_namej = STRING_ELT(Rule_names, snpsj);
	if (test_size==1)
	  SET_STRING_ELT(TestNames, test, Rule_namej);
	else
	  SET_STRING_ELT(Snames, j, Rule_namej);
	SEXP Rule =  VECTOR_ELT(Rules, snpsj);
	if (!isNull(Rule)){ /* Not monomorphic */
	  do_impute(Z, N, diploid, NULL, N, name_index, Rule, gt2ht, zw+ij, 
		    NULL);
	  for (int i=0; i<N; i++, ij++) {
	    if (ISNA(zw[ij])) {
	      zw[ij] = 0.0;
	      if (prior[i]) {
		missing ++;
		prior[i] = 0.0;
	      }
	    }
	  }   
	}
	else {
	  for (int i=0; i<N; i++, ij++)
	    zw[ij] = 0.0;
	}
      }
    }

    if (missing == N) { 

      /* No data */
      
      warning("No valid data for test %d", test+1);
      memset(u, 0x00, nsnpt*sizeof(double));
      memset(v, 0x00, sizeof(double)*(nsnpt*(nsnpt+1))/2);

    }
    else if (missing > mallowed) {

      /* Refit the model using current fitted values as start */

      for (int i=0; i<N; i++) 
	fitted[i] = fitted_all[i];
      int err = glm_fit(fam, lnk, N, M, 0, S, y, prior, x, stratum, 
			maxit, epsilon, r2Max, 1,
			&rank, xb, fitted, resid, weights, 
			&scale, &df_r, NULL, NULL, NULL, NULL);
      if (err) 
	warning("No convergence while fitting base model for test %d", test+1);
      glm_score_test(N, rank, S, stratum, nsnpt, zw, C, cluster,
		     resid, weights, xb, scale,
		     r2Max, u, v);
 
    }
    else {

      /* Use  GLM fit to total data */

      df_r = df_r_all - missing;
      rank = rank_all;
      for (int i=0; i<N; i++) 
	weights[i] = prior[i]? weights_all[i]: 0.0;
      glm_score_test(N, rank, S, stratum, nsnpt, zw, C, cluster,
		     resid_all, weights, xb_all, scale_all,
		     r2Max, u, v);
    }
 
    int nu = 0;
    for (int i=0; i<N; i++)
      if(weights[i]) nu++;
    nused[test] = nu;
    if (nsnpt>1) {
      if (qform(nsnpt, u, v, NULL, chisq+test, df+test)) {
	warning("Matrix not positive semi-definite in test ", test+1);
	chisq[test] = NA_REAL;
	df[test] = NA_INTEGER;
      }
      else if (!df[test])
	chisq[test] = NA_REAL;
    }
    else {
      if ((*v)>0.0) {
	chisq[test] = (*u)*(*u)/(*v);
	df[test] = 1;
      }
      else {
	chisq[test] = NA_REAL;
	df[test] = 0;
      }
    }
    if (if_score) {
      SEXP Scoret;
      PROTECT(Scoret = allocVector(VECSXP, 2));
      setAttrib(Scoret, R_NamesSymbol, UVnames);
      SET_VECTOR_ELT(Scoret, 0, U);
      SET_VECTOR_ELT(Scoret, 1, V);
      SET_VECTOR_ELT(Score, test, Scoret);
      UNPROTECT(3); /* Scoret, U and V */
    }
    else {
      Free(u);
      Free(v);
    }
    if (test_size>1)
      UNPROTECT(1); /* Snames */
  }

  /* Return hash table memory and gt->ht tables */

  if (name_index) {
    index_destroy(name_index);
    for (int i=0; i<pmax; i++)
      destroy_gtype_table(gt2ht[i], i+1);
    Free(gt2ht);
  }

  SEXP Class, Package;
  PROTECT(Class = allocVector(STRSXP, 1));
  if (if_score) 
    SET_STRING_ELT(Class, 0, mkChar("GlmTestsScore"));
  else
    SET_STRING_ELT(Class, 0, mkChar("GlmTests"));
  PROTECT(Package = allocVector(STRSXP, 1));
  SET_STRING_ELT(Package, 0, mkChar("snpStats"));
  setAttrib(Class, install("package"), Package);
  classgets(Result, Class);
  nprot += 2;

  UNPROTECT(nprot);

  return(Result);
  
}

/* Pooling two objects */

SEXP pool2_glm(SEXP X, SEXP Y, SEXP If_score) {
  
  SEXP Xs = R_do_slot(X, mkString("score"));
  SEXP Ys = R_do_slot(Y, mkString("score"));
  SEXP Xn = R_do_slot(X, mkString("N"));
  int *xn = INTEGER(Xn);
  SEXP Yn = R_do_slot(Y, mkString("N"));
  int *yn = INTEGER(Yn);
  SEXP Test_names = R_do_slot(X, mkString("snp.names"));
  SEXP Var_names = R_do_slot(X, mkString("var.names"));
  int N = LENGTH(Xs);
  if (N!=LENGTH(Ys))
    error("pool2_glm: unequal length arguments");
  int if_score = *LOGICAL(If_score);

  SEXP Result, Chisq, Df, Nused, Score = R_NilValue, UVnames=R_NilValue;
  PROTECT(Result = allocS4Object());
  PROTECT(Chisq = allocVector(REALSXP, N));
  double *chisq = REAL(Chisq);
  PROTECT(Df = allocVector(INTSXP, N));
  int *df = INTEGER(Df);
  PROTECT(Nused = allocVector(INTSXP, N));
  int *nused = INTEGER(Nused);
  int nprot = 4;  
  
  if (if_score) {
    PROTECT(Score = allocVector(VECSXP, N));
    setAttrib(Score, R_NamesSymbol, Test_names);
    PROTECT(UVnames = allocVector(STRSXP, 2));
    SET_STRING_ELT(UVnames, 0, mkChar("U"));
    SET_STRING_ELT(UVnames, 1, mkChar("V"));
    nprot += 2;
  }

  for (int i=0; i<N; i++) {
    SEXP Xsi = VECTOR_ELT(Xs, i);
    SEXP Ysi = VECTOR_ELT(Ys, i);
    SEXP XiU = VECTOR_ELT(Xsi, 0);
    double *xiu = REAL(XiU);
    SEXP XiV = VECTOR_ELT(Xsi, 1);
    double *xiv = REAL(XiV);
    SEXP YiU = VECTOR_ELT(Ysi, 0);
    double *yiu = REAL(YiU);
    SEXP YiV = VECTOR_ELT(Ysi, 1);
    double *yiv = REAL(YiV);
    int nu = LENGTH(XiU);
    int nv = LENGTH(XiV);
    if (LENGTH(YiU)!=nu)
      error("attempt to pool tests on unequal numbers of parameters");
    SEXP RU = R_NilValue, RV = R_NilValue;
    double *ru, *rv;
    if (if_score) {
      PROTECT(RU = allocVector(REALSXP, nu));
      ru = REAL(RU);
      PROTECT(RV = allocVector(REALSXP, nv));
      rv = REAL(RV);
    }
    else {
      ru = (double *) Calloc(nu, double);
      rv = (double *) Calloc(nv, double);
    }
    memset(ru, 0x00, nu*sizeof(double));
    memset(rv, 0x00, nv*sizeof(double));
    for (int j=0; j<nu; j++) 
      ru[j] = xiu[j] + yiu[j];
    for (int j=0; j<nv; j++) 
      rv[j] = xiv[j] + yiv[j];
    if (nu>1) {
      if (qform(nu, ru, rv, NULL, chisq+i, df+i)) {
	warning("Matrix not positive semi-definite in pooled test ", i+1);
	chisq[i] = NA_REAL;
	df[i] = NA_INTEGER;
      }
      else if (df[i]==0)
	chisq[i] = NA_REAL;
    }
    else {
      if ((*rv)==0) {
	df[i] = 0;
	chisq[i] = NA_REAL;
      }
      else {
	df[i] = 1;
	chisq[i] = (*ru)*(*ru)/(*rv);
      }
    }
    nused[i] = xn[i] + yn[i];
    if (if_score) {
      SEXP Scorei;
      PROTECT(Scorei = allocVector(VECSXP, 2));
      setAttrib(Scorei, R_NamesSymbol, UVnames);
      /* SEXP Unames = getAttrib(XiU, R_NamesSymbol); unused */
      SET_VECTOR_ELT(Scorei, 0, RU);
      SET_VECTOR_ELT(Scorei, 1, RV);
      SET_VECTOR_ELT(Score, i, Scorei);
      UNPROTECT(3); /* RU, RV, Scorei */
    }
    else {
      Free(ru);
      Free(rv);
    }
  }
  R_do_slot_assign(Result, mkString("snp.names"), Test_names);
  R_do_slot_assign(Result, mkString("var.names"), Var_names);
  R_do_slot_assign(Result, mkString("chisq"), Chisq);
  R_do_slot_assign(Result, mkString("df"), Df);
  R_do_slot_assign(Result, mkString("N"), Nused);
  
  SEXP Class, Package;
  PROTECT(Class = allocVector(STRSXP, 1));
  if (if_score) {
    R_do_slot_assign(Result, mkString("score"), Score);
    SET_STRING_ELT(Class, 0, mkChar("GlmTestsScore"));
  }
  else
    SET_STRING_ELT(Class, 0, mkChar("GlmTests"));
  PROTECT(Package = allocVector(STRSXP, 1));
  SET_STRING_ELT(Package, 0, mkChar("snpStats"));
  setAttrib(Class, install("package"), Package);
  classgets(Result, Class);
  nprot += 2;

  UNPROTECT(nprot);


  return(Result);
}

    
/* Estimation */

SEXP snp_lhs_estimate(const SEXP Y, const SEXP X, const SEXP Stratum, 
		      const SEXP Snp_subset, const SEXP P_est, 
		      const SEXP Robust, const SEXP Cluster, 
		      const SEXP Uncertain, const SEXP Control){

  /* Y should be a SnpMatrix or an XSnpMatrix */
  const char *classY = NULL;
  if (TYPEOF(R_data_class(Y, FALSE)) == STRSXP) {
    classY = CHAR(STRING_ELT(R_data_class(Y, FALSE), 0));
  } else {
    classY = CHAR(STRING_ELT(getAttrib(Y, R_ClassSymbol), 0));
  }
  if(!IS_S4_OBJECT(Y)) {
    error("Y in snp_lhs_score is missing S4 Object bit");
  }
  int ifX = 0;
  if (!strcmp(classY, "SnpMatrix"))
    ifX = 0;
  else if (!strcmp(classY, "XSnpMatrix"))
    ifX = 1;
  else 
    error("Argument error - class(Y)");

  /* Y and its dimensions */

  if (TYPEOF(Y)!=RAWSXP)
    error("Argument error - Y");
  const unsigned char *y = RAW(Y);
  const int *dim;
  int N, nsnp;
  if (strlen(classY)>5) {
    dim = INTEGER(getAttrib(Y, R_DimSymbol));
    N = dim[0];
    nsnp = dim[1];
  }
  else {
    N = LENGTH(Y);
    nsnp = 1;
  }

  /* SNP subset */

  int *snp_subset = NULL;
  int nest = nsnp;
  SEXPTYPE sstype = TYPEOF(Snp_subset);
  if (sstype==INTSXP) {
    if (LENGTH(Snp_subset)>nsnp)
      error("Dimension error - Snp_subset");
    snp_subset = INTEGER(Snp_subset);
    nest = LENGTH(Snp_subset);
  }
  else if (sstype!=NILSXP)
    error("Argument error - Snp.subset");
  

  /* If X-chromosome, indicators of diploid genotype */

  const int *diploid;
  if (ifX) {
    SEXP Diploid = R_do_slot(Y, mkString("diploid")); 
    diploid = LOGICAL(Diploid);
  }
  else
    diploid = NULL;
  
  /* X and its dimensions */

  if (TYPEOF(X)!=REALSXP)
    error("Argument error - X");
  double *x = REAL(X);
  dim = INTEGER(getAttrib(X, R_DimSymbol));
  int M = dim[1];
  if (dim[0]!=N)
    error("Dimension error - X"); 
  SEXP Xnames = VECTOR_ELT(getAttrib(X, R_DimNamesSymbol), 1);

  /* Stratum */

  int S = 1;
  const int *stratum = NULL;
  if (TYPEOF(Stratum)==INTSXP) {
    if (LENGTH(Stratum)!=N) 
      error("Dimension error - Stratum"); 
    stratum = INTEGER(Stratum);
    for (int i=0; i<N; i++) 
      if (stratum[i]>S) S = stratum[i];
  }
  else if (TYPEOF(Stratum)!=NILSXP)
    error("Argument error - Stratum");

  /* P */

  if (TYPEOF(P_est)!=INTSXP)
    error("Argument error - P");
  int P = *INTEGER(P_est);
  if (P<=0 || P>M)
    error("Invalid number of parameters to be estimated");
  
  /* Robust */

  if (TYPEOF(Robust)!=LGLSXP)
    error("Argument error - Robust");
  int robust = LOGICAL(Robust)[0];

  /* Cluster */

  int C = robust? 1: 0;
  int *cluster = NULL;
  if (TYPEOF(Cluster)==INTSXP) {
    if (LENGTH(Cluster)!=N)
      error("Dimension error - Cluster"); 
    cluster = INTEGER(Cluster);
    for (int i=0; i<N; i++)
      if (cluster[i]>C) C = cluster[i];
  }
  else if (TYPEOF(Cluster)!=NILSXP)
    error("Argument error - Cluster");

  /* Handling of uncertain genotypes */

  if (TYPEOF(Uncertain) != LGLSXP)
    error("Argument error: Uncertain is wrong type");
  int uncert = *LOGICAL(Uncertain);
  /* Force robust variance if uncertain option used */
  robust = robust || uncert; 
  
  /* Control object */

  if (TYPEOF(Control)!=VECSXP || LENGTH(Control)!=3) 
    error("Argument error - Control");
  SEXP Maxit = VECTOR_ELT(Control, 0);
  if (TYPEOF(Maxit)!=INTSXP || length(Maxit)!=1)
    error("Argument error - Maxit");
  int maxit = INTEGER(Maxit)[0];
  SEXP Epsilon = VECTOR_ELT(Control, 1);
  if (TYPEOF(Epsilon)!=REALSXP || LENGTH(Epsilon)!=1)
    error("Argument error - Epsilon");
  double epsilon = REAL(Epsilon)[0];
  SEXP R2Max = VECTOR_ELT(Control, 2);
  if (TYPEOF(R2Max)!=REALSXP || LENGTH(R2Max)!=1)
    error("Argument error - R2Max");
  double r2Max = REAL(R2Max)[0];
  
  /* Work arrays */

  double *yd = (double *) R_alloc(N, sizeof(double));
  double *fitted = (double *) R_alloc(N, sizeof(double));
  double *prior = (double *) R_alloc(N, sizeof(double));
  double *resid = (double *) R_alloc(N, sizeof(double));
  double *weights = (double *) R_alloc(N, sizeof(double));
  double *xb = (double *) R_alloc(N*M, sizeof(double));
  
  int *which = (int *) R_alloc(P, sizeof(int));
  double *beta = (double *) R_alloc(P, sizeof(double));
  double *tri = (double *) R_alloc((P*(P+1))/2, sizeof(double));

  /* Output arrays */

  SEXP Estimates, Rnames, BVnames;
  PROTECT(Estimates = allocVector(VECSXP, nest));
  
  SEXP snpNames = VECTOR_ELT(getAttrib(Y, R_DimNamesSymbol), 1);
  SEXPTYPE sntype = TYPEOF(snpNames);
  PROTECT(Rnames = allocVector(sntype, nest));
  
  PROTECT(BVnames = allocVector(STRSXP, 4));
  SET_STRING_ELT(BVnames, 0, mkChar("beta"));
  SET_STRING_ELT(BVnames, 1, mkChar("Var.beta"));
  SET_STRING_ELT(BVnames, 2, mkChar("N"));
  SET_STRING_ELT(BVnames, 3, mkChar("Y.var"));

  /* Do calculations */


  for (int t=0; t<nest; t++) {
    int j = snp_subset? snp_subset[t] - 1: t;

     /* Name of model and yvar */

    SEXP Yname;
    PROTECT(Yname = allocVector(STRSXP, 1));
    SET_STRING_ELT(Yname, 0, mkChar(CHAR(STRING_ELT(snpNames,j))));

    if (sntype==INTSXP) 
      INTEGER(Rnames)[t] = INTEGER(snpNames)[j];
    else 
      SET_STRING_ELT(Rnames, t, STRING_ELT(Yname, 0));

    const unsigned char *yj = y + N*j;
    int mono = 1, yv = 0;
 

    /* Load SNP as Binomial y-variate, with prior weights */

    for (int i=0; i<N; i++) {
      int yij = (int) yj[i];
      if (yij && (uncert||yij<4)) {
	if (mono) {
	  if (yij<4) {
	    mono = !yv || (yij==yv);
	    yv = yij;
	  }
	  else 
	    mono = 0;
	}
	prior[i] = (!diploid || diploid[i])? 2.0: 1.0;
	yd[i] = g2mean(yij)/2.0;
      }
      else {
	prior[i] = yd[i] = 0.0;
      }
    }

    if (!mono) { 
      
    /* Fit model */

      int rank, dfr, pest;
      double scale;
      double rc = glm_fit(BINOMIAL, LOGIT, N, M, P, S, yd, prior, x, stratum, 
			  maxit, epsilon, r2Max, 0, 
			  &rank, xb, fitted, resid, weights, &scale, &dfr,
			  &pest, which, beta, tri);
      if (rc) 
	warning("Failure to converge while fitting model for SNP %d",j+1);
      if (!pest) {
	warning("No estimable parameters for test %d", t+1);
	SET_VECTOR_ELT(Estimates, t, R_NilValue);
      }
      else {
	SEXP B, Bnames, V, Nu;
	double *b, *v;
	PROTECT(B = allocVector(REALSXP, pest));
	PROTECT(Bnames = allocVector(STRSXP, pest));
	PROTECT(V = allocVector(REALSXP, (pest*(pest+1))/2));
	PROTECT(Nu = allocVector(INTSXP, 1));
	b = REAL(B);
	v = REAL(V);
	if (C) {
	  meat_matrix(N, pest, C, cluster, xb+N*(rank-pest), resid, weights, v);
	  glm_est(pest, beta, tri, scale, v, b, v);
	}
	else   
	  glm_est(pest, beta, tri, scale, NULL, b, v);
	int nunit = 0;
	for (int i=0; i<N; i++) 
	  if (weights[i]) nunit++;
	*INTEGER(Nu) = nunit;
      
	/* Names for B */
	
	for (int i=0; i<pest; i++) 
	  SET_STRING_ELT(Bnames, i, STRING_ELT(Xnames, which[i]));
	setAttrib(B, R_NamesSymbol, Bnames);
	
	SEXP Estt;
	PROTECT(Estt = allocVector(VECSXP, 4));
	setAttrib(Estt, R_NamesSymbol, BVnames);
	SET_VECTOR_ELT(Estt, 0, B);
	SET_VECTOR_ELT(Estt, 1, V);	
	SET_VECTOR_ELT(Estt, 2, Nu);
	SET_VECTOR_ELT(Estt, 3, Yname);
	SET_VECTOR_ELT(Estimates, t, Estt);
	UNPROTECT(5);
      }
    }
    else {
      warning("Monomorphic SNP: %d", t+1);
      SET_VECTOR_ELT(Estimates, t, R_NilValue);
    }
    UNPROTECT(1);
  }
	  
  setAttrib(Estimates, R_NamesSymbol, Rnames);  
  SEXP Class, Package, Lhs;
  PROTECT(Class = allocVector(STRSXP, 1));
  SET_STRING_ELT(Class, 0, mkChar("GlmEstimates"));
  PROTECT(Package = allocVector(STRSXP, 1));
  SET_STRING_ELT(Package, 0, mkChar("snpStats"));
  setAttrib(Class, install("package"), Package);
  classgets(Estimates, Class);
  PROTECT(Lhs = allocVector(LGLSXP, 1));
  LOGICAL(Lhs)[0] = 1;
  setAttrib(Estimates, install("snpLHS"), Lhs);
  UNPROTECT(6);
  SET_S4_OBJECT(Estimates);
  return(Estimates);
}

SEXP snp_rhs_estimate(SEXP Y, SEXP family, SEXP link, 
		      SEXP X, SEXP Stratum, SEXP Z, SEXP Rules, 
		      SEXP Prior, SEXP Sets, SEXP Robust, SEXP Cluster, 
		      SEXP Uncertain, SEXP Control) {

  char setname[MAX_NAME_LENGTH];
  int max_name_length =  MAX_NAME_LENGTH -1;
  setname[max_name_length] = (char) 0; /* Make sure its null terminated */

  /* Y and its dimensions and label */

  if (TYPEOF(Y)!=REALSXP)
    error("Argument error - Y");
  double *y = REAL(Y);
  int N = LENGTH(Y);
  SEXP Yname;
  PROTECT(Yname = getAttrib(Y, install("label")));

  /* Family */

  if (TYPEOF(family)!=INTSXP || LENGTH(family)!=1)
    error("Argument error - family");
  int fam = INTEGER(family)[0];
  if (fam<1 || fam>4)
    error("Illegal family argument");
  
  /* Link */

  if (TYPEOF(link)!=INTSXP || LENGTH(link)!=1)
    error("Argument error - link");
  int lnk = INTEGER(link)[0];
  if (lnk<1 || lnk>4)
    error("Illegal link argument");

  /* X and its dimensions */

  double *x = NULL;
  int M = 0;
  int *dim = NULL;
  if (TYPEOF(X)==REALSXP) {
    x = REAL(X);
    dim = INTEGER(getAttrib(X, R_DimSymbol));
    M = dim[1];
    if (dim[0]!=N)
      error("Dimension error - X"); 
  }
  else if (TYPEOF(X)==NILSXP) {
    M = 0; 
    x = NULL;
  }
  else
    error("Argument error - X");

  /* Stratum */

  int S = 1;
  int *stratum = NULL;
  if (TYPEOF(Stratum)==INTSXP) {
    if (LENGTH(Stratum)!=N) 
      error("Dimension error - Stratum"); 
    stratum = INTEGER(Stratum);
    for (int i=0; i<N; i++) 
      if (stratum[i]>S) S = stratum[i];
  }
  else if (TYPEOF(Stratum)!=NILSXP)
    error("Argument error - Stratum");

/* Z should be a SnpMatrix or an XSnpMatrix */

  const char *classZ = NULL;
  if (TYPEOF(R_data_class(Z, FALSE)) == STRSXP) {
    classZ = CHAR(STRING_ELT(R_data_class(Z, FALSE), 0));
  } else {
    classZ = CHAR(STRING_ELT(getAttrib(Z, R_ClassSymbol), 0));
  }
  if(!IS_S4_OBJECT(Z)) {
    error("Z in snp_rhs_score is missing S4 Object bit");
  }
  int ifX = 0;
  if (!strcmp(classZ, "SnpMatrix"))
    ifX = 0;
  else if (!strcmp(classZ, "XSnpMatrix"))
    ifX = 1;
  else 
    error("Argument error - class(Z)");

  /* Z and its dimensions */

  int Nz;
  if (TYPEOF(Z)!=RAWSXP)
    error("Argument error - Z");
  const unsigned char *z = RAW(Z);
  int nsnp;
  if (strlen(classZ)>5) {
    dim = INTEGER(getAttrib(Z, R_DimSymbol));
    Nz = dim[0];
    nsnp = dim[1];
  }
  else {
    Nz = LENGTH(Z);
    nsnp = 1;
  }
  if (Nz!=N)
    error("Dimension error - Z");

  /* If X-chromosome, indicators of diploid genotype */

  int *diploid;
  if (ifX) {
    SEXP Diploid = R_do_slot(Z, mkString("diploid")); 
    diploid = LOGICAL(Diploid);
  }
  else
    diploid = NULL;
  
  /* If imputation involved, calculate snp name index and lookup tables */

  index_db name_index = NULL;
  SEXP Snp_names =  VECTOR_ELT(getAttrib(Z, R_DimNamesSymbol), 1);
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

  /* Prior weights */

  double *prior = NULL;
  if (TYPEOF(Prior)==REALSXP) {
    if (LENGTH(Prior)!=N)
      error("Dimension error - Prior");
    prior = REAL(Prior);
  }
  else if (TYPEOF(Prior)!=NILSXP)
    error("Argument error - Prior");

  /* Estimates. For now snps referred to by col number not name */

  int set_size = 1, gen_names = 1;
  int *set_int = NULL;
  int nset = nsnp;
  SEXP NameSets = R_NilValue;
  SEXPTYPE sets_type = TYPEOF(Sets);
  if (sets_type==VECSXP) {
    nset = LENGTH(Sets);
    NameSets = getAttrib(Sets, R_NamesSymbol);
    gen_names = (NameSets==R_NilValue);
    for (int i=0; i<nset; i++) {
      SEXP seti =  VECTOR_ELT(Sets, i);
      if (TYPEOF(seti)!=INTSXP)
	error("Non-integer list - Sets[i]");
      int leni = LENGTH(seti);
      if (leni>set_size)
	set_size = leni;
    }
    if (set_size==1)
      error("List of multi-SNP sets doesn't contain a multi-SNP spec");
  }
  else if (sets_type==INTSXP) {
    set_int = INTEGER(Sets);
    nset = LENGTH(Sets);
  }
  else if (sets_type!=NILSXP)
    error("Argument error - Sets");

  /* Robust */

  if (TYPEOF(Robust)!=LGLSXP)
    error("Argument error - Robust");
  int robust = LOGICAL(Robust)[0];

  /* Cluster */

  int C = robust? 1: 0;
  int *cluster = NULL;
  if (TYPEOF(Cluster)==INTSXP) {
    if (LENGTH(Cluster)!=N)
      error("Dimension error - Cluster"); 
    cluster = INTEGER(Cluster);
    for (int i=0; i<N; i++)
      if (cluster[i]>C) C = cluster[i];
  }
  else if (TYPEOF(Cluster)!=NILSXP)
    error("Argument error - Cluster");
  
  /* Handling of uncertain genotypes */

  if (TYPEOF(Uncertain) != LGLSXP)
    error("Argument error: Uncertain is wrong type");
  int uncert = *LOGICAL(Uncertain);

  /* Control object */

  if (TYPEOF(Control)!=VECSXP || LENGTH(Control)!=3) 
    error("Argument error - Control");
  SEXP Maxit = VECTOR_ELT(Control, 0);
  if (TYPEOF(Maxit)!=INTSXP || length(Maxit)!=1)
    error("Argument error - Maxit");
  int maxit = INTEGER(Maxit)[0];
  SEXP Epsilon = VECTOR_ELT(Control, 1);
  if (TYPEOF(Epsilon)!=REALSXP || LENGTH(Epsilon)!=1)
    error("Argument error - Epsilon");
  double epsilon = REAL(Epsilon)[0];
  SEXP R2Max = VECTOR_ELT(Control, 2);
  if (TYPEOF(R2Max)!=REALSXP || LENGTH(R2Max)!=1)
    error("Argument error - R2Max");
  double r2Max = REAL(R2Max)[0];

  /* Work arrays */

  int M_max = M+set_size;
  double *new_prior = (double *) R_alloc(N, sizeof(double));
  double *resid = (double *) R_alloc(N, sizeof(double));
  double *weights = (double *) R_alloc(N, sizeof(double));
  double *fitted = (double *) R_alloc(N, sizeof(double));
  double *xb = (double *) R_alloc(N*M_max, sizeof(double));
  double *xz = (double *) R_alloc(N*M_max, sizeof(double));

  double *beta =  (double *) R_alloc(M_max, sizeof(double));
  int *which = (int *) R_alloc(M_max, sizeof(int));
  double *tri = (double *) R_alloc((M_max*(M_max+1))/2, sizeof(double));
  double *Uc;

  if (C>1) 
    Uc = (double *) R_alloc(C*M_max, sizeof(double));

  /* Copy base model */

  int NM = N*M;
  for (int ij=0; ij<NM; ij++)
    xz[ij] = x[ij];
    
  /* Output list */

  SEXP Estimates, BVnames;
  PROTECT(Estimates = allocVector(VECSXP, nset));
  
  PROTECT(BVnames = allocVector(STRSXP, 4));
  SET_STRING_ELT(BVnames, 0, mkChar("beta"));
  SET_STRING_ELT(BVnames, 1, mkChar("Var.beta"));
  SET_STRING_ELT(BVnames, 2, mkChar("N"));
  SET_STRING_ELT(BVnames, 3, mkChar("Y.var"));

  SEXP SetNames = R_NilValue;
  if (gen_names) {
    PROTECT(SetNames = allocVector(STRSXP, nset));
    setAttrib(Estimates, R_NamesSymbol, SetNames);
  }
  else 
    setAttrib(Estimates, R_NamesSymbol, NameSets);
  
  /* Do calculations */

  int snp;
  int *snps = &snp;
  for (int set=0; set<nset; set++) {

    int nsnp_set = 1;
    if (set_size==1) { /* Single SNP sets */
      if (set_int)  /* Selected list */
	snps = set_int+set;
      else
	snp = set + 1;
    }
    else  { /* List of multi-SNP sets */
      SEXP Snps = VECTOR_ELT(Sets, set);
      nsnp_set = LENGTH(Snps);
      snps = INTEGER(Snps);
    }
    for (int i=0; i<N; i++) new_prior[i] = prior? prior[i]: 1.0;
 
    /* Add to XZ matrix, tracking incomplete cases  */

    int missing = 0;
    for (int j=0, ij=NM; j<nsnp_set; j++) {
      int snpsj = snps[j];
      const char *srname;
      if (snpsj>0) {
	snpsj--; 
	srname = CHAR(STRING_ELT(Snp_names, snpsj));
	const unsigned char *zj = z + N*snpsj;
	for (int i=0; i<N; i++, ij++) {
	  unsigned char zij = zj[i];
	  if (zij && (uncert||(zij<4))) {
	    xz[ij] = g2mean(zij);
	  }
	  else {
	    xz[ij] = 0.0;
	    if (new_prior[i]){
	      missing ++;
	      new_prior[i] = 0.0;
	    }
	  }
	}
      }	      		      
      else {
	snpsj = -(1+snpsj);
	SEXP Rule_namej = STRING_ELT(Rule_names, snpsj);
	srname = CHAR(Rule_namej);
	SEXP Rule =  VECTOR_ELT(Rules, snpsj);
	if (!isNull(Rule)){ /* Not monomorphic */
	  do_impute(Z, N, diploid, NULL, N, name_index, Rule, gt2ht, xz+ij, 
		    NULL);
	  for (int i=0; i<N; i++, ij++) {
	    if (ISNA(xz[ij])) {
	      xz[ij] = 0.0;
	      if (new_prior[i]) {
		missing ++;
		new_prior[i] = 0.0;
	      }
	    }
	  } 
	}
	else {
	  for (int i=0; i<N; i++, ij++)
	    xz[ij] = 0.0;
	}
      }
      if (gen_names) {
	if (!j){
	  strncpy(setname, srname, max_name_length);
	}
	else if ((nsnp_set>1) && (j==(nsnp_set-1))) {
	  int len = strlen(setname);
	  int space = max_name_length - len;
	  if (space>2) {
	    setname[len++] = '.';
	    setname[len++] = '.';
	  }
	  strncpy(setname+len, srname, space-2);
	}
      }
    }
    if (gen_names)
      SET_STRING_ELT(SetNames, set, mkChar(setname));
 
    if (missing == N) { /* No data */
      warning("No valid data for set %d", set+1);
      SET_VECTOR_ELT(Estimates, set, R_NilValue);
    }
    else {
      int rank, pest, df_r;
      double scale;
      int err = glm_fit(fam, lnk, N, M+nsnp_set, nsnp_set, S, 
			y, new_prior, xz, stratum, 
			maxit, epsilon, r2Max, 0,
			&rank, xb, fitted, resid, weights, 
			&scale, &df_r, &pest, which, beta, tri);
      if (err) 
	warning("No convergence while fitting model for set %d", set+1);
      if (!pest) {
	warning("No estimable parameters for set %d", set+1);
	SET_VECTOR_ELT(Estimates, set, R_NilValue);
      }
      else {
	SEXP B = R_NilValue, V = R_NilValue, Nu= R_NilValue, 
	  Snames = R_NilValue;
	PROTECT(B = allocVector(REALSXP, pest));
	PROTECT(V = allocVector(REALSXP, (pest*(pest+1))/2));
	PROTECT(Nu = allocVector(INTSXP, pest));
	double *b = REAL(B);
	double *v = REAL(V);
	PROTECT(Snames = allocVector(STRSXP, pest));
	if (C) {
	  meat_matrix(N, pest, C, cluster, xb+N*(rank-pest), resid, weights, v);
	  glm_est(pest, beta, tri, scale, v, b, v);
	}
	else   
	  glm_est(pest, beta, tri, scale, NULL, b, v);
	int nunit = 0;
	for (int i=0; i<N; i++) 
	  if (weights[i]) nunit++;
	*INTEGER(Nu) = nunit;
      

	/* Names for B */
	for (int i=0; i<pest; i++) {
	  int snpswi = snps[which[i]-M];
	  if (snpswi>0) 
	    SET_STRING_ELT(Snames, i, 
			 STRING_ELT(Snp_names, snpswi-1));
	  else
	    SET_STRING_ELT(Snames, i, 
			   STRING_ELT(Rule_names, -(snpswi+1)));
	}
	setAttrib(B, R_NamesSymbol, Snames);
	SEXP Estt;
	PROTECT(Estt = allocVector(VECSXP, 4));
	setAttrib(Estt, R_NamesSymbol, BVnames);
	SET_VECTOR_ELT(Estt, 0, B);
	SET_VECTOR_ELT(Estt, 1, V);	
	SET_VECTOR_ELT(Estt, 2, Nu);
	SET_VECTOR_ELT(Estt, 3, Yname);
	SET_VECTOR_ELT(Estimates, set, Estt);
	UNPROTECT(5);
      }
    }
  }

  /* Return hash table memory and gt->ht tables */

  if (name_index) {
    index_destroy(name_index);
    for (int i=0; i<pmax; i++)
      destroy_gtype_table(gt2ht[i], i+1);
    Free(gt2ht);
  }

  SEXP Class, Package, Lhs;
  PROTECT(Class = allocVector(STRSXP, 1));
  SET_STRING_ELT(Class, 0, mkChar("GlmEstimates"));
  PROTECT(Package = allocVector(STRSXP, 1));
  SET_STRING_ELT(Package, 0, mkChar("snpStats"));
  setAttrib(Class, install("package"), Package);
  classgets(Estimates, Class);
  PROTECT(Lhs = allocVector(LGLSXP, 1));
  LOGICAL(Lhs)[0] = 0;
  setAttrib(Estimates, install("snpLHS"), Lhs);

  UNPROTECT(gen_names? 7: 6);

  SET_S4_OBJECT(Estimates);
  return(Estimates);
  
}

/* Wald tests from estimates list element */


SEXP wald(const SEXP From) {
  SEXP LHS = getAttrib(From, install("snpLHS"));
  int lhs = LOGICAL(LHS)[0];
  int ntest = LENGTH(From);
  int maxsnp, maxwork;
  SEXP To, Snp_names, Var_names, Chisq, Df, N;

  /* Find first valid entry */
  SEXP F1 = R_NilValue;
  for (int i=0; i<ntest; i++) {
    F1 = VECTOR_ELT(From, i);
    if (!isNull(F1))
      break;
  }

  SEXP From_names = getAttrib(From, R_NamesSymbol);
  if (lhs) {
    maxsnp = 1;
    Snp_names = From_names;
    Var_names = getAttrib(VECTOR_ELT(F1, 0), R_NamesSymbol);
    maxwork = LENGTH(Var_names);
  }
  else {
    /* Maximum size */
    maxsnp = 0;
    for (int i=0; i<ntest; i++) {
      SEXP Fi = VECTOR_ELT(From, i);
      if (!isNull(Fi)) {
	int nb = LENGTH(VECTOR_ELT(Fi, 0));
	if (nb>maxsnp)
	  maxsnp = nb;
      }
    }
    if (maxsnp>1) {
      PROTECT(Snp_names = allocVector(VECSXP, ntest));
      setAttrib(Snp_names, R_NamesSymbol,  getAttrib(From, R_NamesSymbol));
    }
    else
      Snp_names = From_names;
    Var_names = VECTOR_ELT(F1, 3);
    maxwork = maxsnp;
  }
  PROTECT(Chisq = allocVector(REALSXP, ntest));
  double *chisq = REAL(Chisq);
  PROTECT(Df = allocVector(INTSXP, ntest));
  int *df = INTEGER(Df);
  PROTECT(N = allocVector(INTSXP, ntest));
  int *n = INTEGER(N);
  PROTECT(To = R_do_new_object(R_getClassDef("GlmTests")));
  /* setAttrib(To, R_NamesSymbol, getAttrib(From, R_NamesSymbol)); */
  R_do_slot_assign(To, mkString("snp.names"), Snp_names);
  R_do_slot_assign(To, mkString("var.names"), Var_names);
  R_do_slot_assign(To, mkString("chisq"), Chisq);
  R_do_slot_assign(To, mkString("df"), Df);
  R_do_slot_assign(To, mkString("N"), N);
  
  /* Loop over fits */

  double *work = NULL;
  if (maxwork>1) 
    work = (double *) Calloc((maxwork*(maxwork+1))/2, double);
  for (int i=0; i<ntest; i++) {
    SEXP Fi = VECTOR_ELT(From, i);
    double wald = NA_REAL;
    int rank = NA_INTEGER, ni = NA_INTEGER; 
    if (!isNull(Fi)) {
      SEXP B = VECTOR_ELT(Fi, 0);
      double *b = REAL(B);
      SEXP B_names = getAttrib(B, R_NamesSymbol);
      if (!lhs) {
	if (maxsnp>1) 
	  SET_VECTOR_ELT(Snp_names, i, B_names);
      }
      int nb = LENGTH(B);
      double *v = REAL(VECTOR_ELT(Fi, 1));
      if (maxwork>1) {
	int res = qform(nb, b, v, work, &wald, &rank);
	if (!res) {
	  wald = NA_REAL;
	  rank = NA_INTEGER;
	}
      }
      else {
	if (v[0]>0.0) {
	  wald = b[0]*b[0]/v[0];
	  rank = 1;
	}
      }
      ni = INTEGER(VECTOR_ELT(Fi, 2))[0];
    }
    chisq[i] = wald;
    df[i] = rank;
    n[i] = ni;
  }

  if (maxwork>1)
    Free(work);
  if (lhs || maxsnp==1)
    UNPROTECT(4);
  else 
    UNPROTECT(5);
 
  return(To);
}

  
