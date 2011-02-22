/* 
   Single SNP tests - 1df and 2df, controlling for a stratification 
   
   Calculates score test and permutation variance (variance under 
   random permutations of phenotype vector within strata) and then the 
   conventional U^2/V chi-squared test
 
*/

#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "hash_index.h"
#include "imputation.h"
#include "uncertain.h"
#include "Rmissing.h"

SEXP score_single(const SEXP Phenotype, const SEXP Stratum, const SEXP Snps, 
		  const SEXP Rules, const SEXP Subset, const SEXP Snp_subset, 
		  const SEXP Uncertain){
  double vec[4] = {1.0, 0.0, 0.0, 0.0};
  /* int i, j, k, km, m, t, su; */

  /* Phenotype */

  if (TYPEOF(Phenotype)!=REALSXP)
    error("Argument error - Phenotype");
  const double *phenotype = REAL(Phenotype);
  int n = LENGTH(Phenotype);
  
  /* Stratum */

  const int *stratum = NULL;
  int nstrata = 1;
  SEXPTYPE stype = TYPEOF(Stratum);
  if (stype==INTSXP) {
    if (LENGTH(Stratum)!=n)
      error("Dimension error - Stratum");
    stratum = INTEGER(Stratum);
    for (int i=0; i<n; i++) {
      int si = stratum[i];
      if (si>nstrata) nstrata = si;
    }
  }
  else if (stype!=NILSXP)
    error("Argument error - Stratum");

  /* SNPs ---- should be a SnpMatrix or an XSnpMatrix */

  const char *classS = NULL;
  if (TYPEOF(R_data_class(Snps, FALSE)) == STRSXP) {
    classS = CHAR(STRING_ELT(R_data_class(Snps, FALSE), 0));
  } else {
    classS = CHAR(STRING_ELT(getAttrib(Snps, R_ClassSymbol), 0));
  }
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
  if (strlen(classS)>5) {
    int *dim = INTEGER(getAttrib(Snps, R_DimSymbol));
    N = dim[0];
    nsnp = dim[1];
  }
  else {
    N = LENGTH(Snps);
    nsnp = 1;
  }
  if (N!=n)
    error("Dimension error - Snps");
  SEXP snp_names = VECTOR_ELT(getAttrib(Snps, R_DimNamesSymbol), 1);
  index_db name_index = create_name_index(snp_names);

  /* Subset */

  int *subset = NULL;
  int nsubj = n;
  SEXPTYPE sutype = TYPEOF(Subset);
  if (sutype==INTSXP) {
    if (LENGTH(Subset)>n)
      error("Dimenion error - Subset");
    subset = INTEGER(Subset);
    nsubj = LENGTH(Subset);
  }
  else if (sutype!=NILSXP)
    error("Argument error - Subset");

  /* Rules */

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

  /* SNP subset */

  int *snp_subset = NULL;
  int ntest = nrules? nrules: nsnp;
  SEXPTYPE sstype = TYPEOF(Snp_subset);
  if (sstype==INTSXP) {
    snp_subset = INTEGER(Snp_subset);
    ntest = LENGTH(Snp_subset);
  }
  else if (sstype!=NILSXP)
    error("Argument error - Snp.subset");

  /* Sex */

  const int *female = NULL;
  if (ifX) {
    SEXP Female = R_do_slot(Snps, mkString("Female")); 
    female = LOGICAL(Female);
  }

  /* Handling of uncertain genotypes */

  if (TYPEOF(Uncertain) != LGLSXP)
    error("Argument error: Uncertain is wrong type");
  int uncert = *LOGICAL(Uncertain);

 
  /* Output objects */
 
  SEXP Result, Used, U, V;
  PROTECT(Result = allocVector(VECSXP, 4));
  
  if (female) {
    PROTECT(U =  allocMatrix(REALSXP, ntest, 3) ); /* Scores */
    PROTECT(V =  allocMatrix(REALSXP, ntest, 4) ); /* Score variances */
  }
  else {
    PROTECT(U =  allocMatrix(REALSXP, ntest, 2) ); /* Scores */
    PROTECT(V =  allocMatrix(REALSXP, ntest, 3) ); /* Score variances */
  }
  SET_VECTOR_ELT(Result, 0, U);
  double *umat = REAL(U); 
  SET_VECTOR_ELT(Result, 1, V);
  double *vmat = REAL(V);

  PROTECT(Used =  allocVector(INTSXP, ntest) ); /* N used */
  int *Nused = INTEGER(Used);
  SET_VECTOR_ELT(Result, 2, Used);

  SEXP Nr2;
  double *nr2;
  if (isNull(Rules)) {
    PROTECT(Nr2 = allocVector(REALSXP, 0));
    nr2 = NULL;
  }
  else {
    PROTECT(Nr2 = allocVector(REALSXP, ntest)); /* N used x R-squared */
    nr2 = REAL(Nr2);
  }
  SET_VECTOR_ELT(Result, 3, Nr2);
  
  /* Space to hold imputed values */

  double *xadd = NULL;
  double *xdom = NULL;
  GTYPE **gt2ht = NULL;
  int pmax = 0;
  if (nrules) {
    xadd = (double *) Calloc(nsubj, double);
    xdom = (double *) Calloc(nsubj, double);
    pmax = *INTEGER(getAttrib(Rules, install("Max.predictors")));
    gt2ht = (GTYPE **)Calloc(pmax, GTYPE *);
    for (int i=0; i<pmax; i++)
      gt2ht[i] = create_gtype_table(i+1);
  }



  /* Do calculations */
  
  if (!ifX) {

    /* Work arrays */
    
    double **UV = (double **) Calloc(nstrata, double *);
    for (int i=0; i<nstrata; i++) 
      UV[i] = (double *) Calloc(10, double);

    for (int t=0; t<ntest; t++) {
      int i = snp_subset? snp_subset[t] - 1: t; 
      const unsigned char *snpsi = NULL;
      double r2 = 0.0;
      if (nrules) {
	if (i >= nrules)
	  error("snp_subset out of range");
	SEXP Rule =  VECTOR_ELT(Rules, i);
	if (isNull(Rule)){ /* Monomorphic */
	  for (int j=0; j<nsubj; j++)
	    xadd[j] = xdom[j] = 0.0;
	}
	else {
	  do_impute(Snps, n, NULL, subset, nsubj, name_index, Rule, gt2ht, 
		    xadd, xdom);
	  r2 = *REAL(VECTOR_ELT(Rule, 1));
	}
      }
      else {
	if (i >= nsnp)
	  error("snp_subset out of range");
	snpsi = snps + n*i;
      }
      /* Initialise score and score variance array */
      for (int j=0; j<nstrata; j++) {
	double *uvj = UV[j];
	for (int k=0; k<10; k++)
	  uvj[k] = 0.0;
      }

      int nu = 0;
      for (int su=0; su<nsubj; su++) {
	int j = subset? subset[su] - 1: su;
	if (j>=n) 
	  error("subset out of range");
	int strat = (nstrata==1)? 0: (stratum[j]-1);
	vec[1] = phenotype[j];
	double *uv = UV[strat];
	if (nrules) {
	  double ax = xadd[su];
	  vec[2] = ax;
	  if (!ISNA(ax)) {
	    nu++;
	    vec[3] = xdom[su];
	    for (int k=0, km=0; k<4; k++) {
	      double vk = vec[k];
	      for (int m=0; m<=k; m++)
		uv[km++] += vk*vec[m];
	    }
	  }
	}
	else {
	  unsigned char zij = snpsi[j];
	  if (zij&&((zij<4)|uncert)) {
	    nu++;
	    g2ad(zij, vec+2, vec+3);
	    for (int k=0, km=0; k<4; k++) {
	      double vk = vec[k];
	      for (int m=0; m<=k; m++)
		uv[km++] += vk*vec[m];
	    }
	  }
	}
      }
      double u1=0.0, u2=0.0, v11=0.0, v12=0.0, v22=0.0; 
      for (int j=0; j<nstrata; j++) {
	double *uvj = UV[j];
	double N = uvj[0];
	if (N>1.0) {
	  u1 += (uvj[4] - uvj[3]*uvj[1]/N);
	  u2 += (uvj[7] - uvj[6]*uvj[1]/N);
	  double vy = (uvj[2] - uvj[1]*uvj[1]/N)/(N-1.0);
	  v11 += vy*(uvj[5] - uvj[3]*uvj[3]/N);
	  v12 += vy*(uvj[8] - uvj[3]*uvj[6]/N);
	  v22 += vy*(uvj[9] - uvj[6]*uvj[6]/N);
	}
      }
      Nused[t] = nu;
      if (nrules)
	nr2[t] = nu*r2;
      umat[t] = u1;
      umat[ntest+t] = u2;
      vmat[t] = v11;
      vmat[ntest+t] = v12;
      vmat[2*ntest+t] = v22;
    }

    /* Return work arrays */

    for (int i=0; i<nstrata; i++) 
      Free(UV[i]);
    Free(UV);
  }
  else {

    /* Work arrays */
    
    double **UVM = (double **) Calloc(nstrata, double *);
    double **UVF = (double **) Calloc(nstrata, double *);
    for (int i=0; i<nstrata; i++) {
      UVM[i] = (double *) Calloc(10, double);
      UVF[i] = (double *) Calloc(10, double);
    }
    
    for (int t=0; t<ntest; t++) {
      int i = snp_subset? snp_subset[t] - 1: t; 
      const unsigned char *snpsi = NULL;
      double r2 = 0.0;
      if (nrules) {
	if (i >= nrules)
	  error("snp_subset out of range");
	SEXP Rule = VECTOR_ELT(Rules, i);
	if (isNull(Rule)) {
	  for (int j=0; j<nsubj; j++)
	    xadd[j] = xdom[j] = 0.0;
	}
	else {
	  do_impute(Snps, n, female, subset, nsubj, name_index, Rule, gt2ht, 
		    xadd, xdom);
	  r2 = *REAL(VECTOR_ELT(Rule, 1));
	}
      }
      else {
	if (i >= nsnp)
	  error("snp_subset out of range");
	snpsi = snps + n*i;
      }
      /* Initialise score and score variance arrays */
      for (int j=0; j<nstrata; j++) {
	double *uvmj = UVM[j];
	double *uvfj = UVF[j];
	for (int k=0; k<10; k++)
	  uvmj[k] = uvfj[k] = 0.0;
      }
      int nu = 0;
      for (int su=0; su<nsubj; su++) {
	int j = subset? subset[su] - 1: su;
	if (j>=n) 
	  error("subset out of range");
	int strat = (nstrata==1)? 0: (stratum[j]-1);
	vec[1] = phenotype[j];
	double *uv = female[j] ? UVF[strat]: UVM[strat];
	if (nrules) {
	  vec[2] = xadd[su];
	  if (!ISNA(vec[2])) {
	    nu++;
	    vec[3] = xdom[su];
	    for (int k=0, km=0; k<4; k++) {
	      double vk = vec[k];
	      for (int m=0; m<=k; m++)
		uv[km++] += vk*vec[m];
	    }
	  }
	}
	else {
	  unsigned char zij = snpsi[j];
	  if (zij&&((zij<4)|uncert)) {
	    nu++;
	    g2ad(zij, vec+2, vec+3);
	    for (int k=0, km=0; k<4; k++) {
	      double vk = vec[k];
	      for (int m=0; m<=k; m++)
		uv[km++] += vk*vec[m];
	    }
	  }
	}
      }
      /* rewrite this loop */
      double u=0, v=0, u1=0.0, u2=0.0, v11=0.0, v12=0.0, v22=0.0;
      for (int j=0; j<nstrata; j++) {
	double *uvmj = UVM[j];
	double Nm = uvmj[0];
	double *uvfj = UVF[j];
	double Nf = uvfj[0];
	double Nt = Nf + Nm;
	if (Nt>0) {
	  double ybar = (uvmj[1] + uvfj[1])/Nt;
	  double yb2 = ybar*ybar;
	  u += uvmj[4] + uvfj[4] - ybar*(uvmj[3]+uvfj[3]);
	  /* Variance of male contribution (Bernoulli) */
	  if (Nm>0) {
	    double ssy = uvmj[2] - 2*ybar*uvmj[1] + Nm*yb2;
	    if (nrules) {
	      /* I think this should be right but needs checking */
	      v += ssy*(uvmj[5] - uvmj[3]*uvmj[3]/Nm)/(Nm - 1.0);
	    }
	    else {
	      double af = (uvfj[3] + uvmj[3]/2)/(2*Nf+Nm);
	      v += 4*af*(1.0-af)*ssy;
	    }
	  }
	  /* Female contribution */
	  if (Nf>1) {
	    u1 += uvfj[4] - uvfj[3]*uvfj[1]/Nf;
	    u2 += uvfj[7] - uvfj[6]*uvfj[1]/Nf;
	    double ssy = (uvfj[2] - 2*ybar*uvfj[1] + Nf*yb2)/(Nf-1.0);
	    double w = (uvfj[5] - uvfj[3]*uvfj[3]/Nf);
	    v += ssy*w;
	    ssy = (uvfj[2] - uvfj[1]*uvfj[1]/Nf)/(Nf-1.0);
	    v11 += ssy*w;
	    v12 += ssy*(uvfj[8] - uvfj[3]*uvfj[6]/Nf);
	    v22 += ssy*(uvfj[9] - uvfj[6]*uvfj[6]/Nf);
	  }
	}
      }
      Nused[t] = nu;
      if (nrules)
	nr2[t] = nu*r2;
      umat[t] = u;
      umat[ntest+t] = u1;
      umat[2*ntest+t] = u2;
      vmat[t] = v;
      vmat[ntest+t] = v11;
      vmat[2*ntest+t] = v12;
      vmat[3*ntest+t] = v22;
    }

    /* Return work arrays */

    for (int i=0; i<nstrata; i++) {
      Free(UVM[i]);
      Free(UVF[i]);
    }
    Free(UVM);
    Free(UVF);
  }

  /* Tidy up */

  index_destroy(name_index);
  if (nrules) {
    Free(xadd);
    Free(xdom);
    for (int i=0; i<pmax; i++) 
      destroy_gtype_table(gt2ht[i], i+1);
    Free(gt2ht);
  }

  /* Attributes of output object */

  SEXP Names;
  PROTECT(Names = allocVector(STRSXP, 4));
  SET_STRING_ELT(Names, 0, mkChar("U"));
  SET_STRING_ELT(Names, 1, mkChar("V"));
  SET_STRING_ELT(Names, 2, mkChar("N"));
  SET_STRING_ELT(Names, 3, mkChar("N.r2"));

  setAttrib(Result, R_NamesSymbol, Names);
  UNPROTECT(6);

  return(Result);
}

SEXP chisq_single(const SEXP Scores) {
  const double tol=1.e-2;
  SEXP U = VECTOR_ELT(Scores, 0);
  SEXP V = VECTOR_ELT(Scores, 1);
  int ntest = nrows(U);
  double *umat = REAL(U);
  double *vmat = REAL(V);
  SEXP Result;
  PROTECT(Result = allocMatrix(REALSXP, ntest, 2));
  double *chisq = REAL(Result);
  if (ncols(U)==3) { /* X chromosome */
    for (int i=0; i<ntest; i++) {
     /* 1 df test */
      double u = umat[i], u1 = umat[ntest+i], u2 = umat[2*ntest+i];
      double v = vmat[i], v11 = vmat[ntest+i], v12 = vmat[2*ntest+i],
	v22 = vmat[3*ntest+i];
      chisq[i] =  v>0.0? u*u/v: NA_REAL;
      /* 2 df test */
      double r2 = v12*v12/(v11*v22);
      if (v11 <= 0.0 || v22 <= 0.0 || (1.0-r2) < tol) 
	chisq[ntest+i] = NA_REAL; /* Not positive definite --- enough! */
      else 
	chisq[ntest+i] = chisq[i] + 
	  (r2*u1*u1/v11 + u2*u2/v22 - 2.0*r2*u1*u2/v12)/(1.0-r2);
    }
  }
  else { /* Autosome */
    for (int i=0; i<ntest; i++) {
      double u1 = umat[i], u2 = umat[ntest+i];
      double v11 = vmat[i], v12 = vmat[ntest+i], v22 = vmat[2*ntest+i];    
      /* 1 df test */
      chisq[i] =  v11>0.0? u1*u1/v11: NA_REAL;
      /* 2 df test */
      double r2 = v12*v12/(v11*v22);
      if (v11 <= 0.0 || v22 <= 0.0 || (1.0-r2) < tol) 
	chisq[ntest+i] = NA_REAL; /* Not positive definite --- enough! */
      else 
	chisq[ntest+i] = (u1*u1/v11 + u2*u2/v22 - 2.0*r2*u1*u2/v12)/(1.0-r2);
    }
  }
  SEXP Dimnames, Colnames;
  PROTECT(Dimnames = allocVector(VECSXP, 2));
  PROTECT(Colnames = allocVector(STRSXP, 2));
  SET_STRING_ELT(Colnames, 0, mkChar("1 df"));
  SET_STRING_ELT(Colnames, 1, mkChar("2 df"));
  SET_VECTOR_ELT(Dimnames, 0, R_NilValue);
  SET_VECTOR_ELT(Dimnames, 1, Colnames);
  setAttrib(Result, R_DimNamesSymbol, Dimnames);
  UNPROTECT(3);
  return Result;
}

 
    
