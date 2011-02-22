/* 
   Single trio-based SNP tests - 1df and 2df
   
   Calculates score test and variance 
 
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

SEXP score_tdt(const SEXP Proband, const SEXP Father, const SEXP Mother, 
	       const SEXP Cluster, const SEXP Snps,  const SEXP Rules, 
               const SEXP Snp_subset, const SEXP Check, const SEXP Robust,
	       const SEXP Uncertain){
 
  int mendelian[36] = {1, 0, 0, 1, 1, 0, 0, 1, 0,
		       1, 1, 0, 1, 1, 1, 0, 1, 1,
		       0, 1, 0, 0, 1, 1, 0, 0, 1,
                       1, 0, 0, 1, 0, 1, 0, 0, 1};
   
  /* Trios */

  if (TYPEOF(Proband)!=INTSXP || TYPEOF(Father)!=INTSXP || 
      TYPEOF(Mother)!=INTSXP || TYPEOF(Cluster)!=INTSXP )
    error("Argument error - Proband|Father|Mother|Cluster");
  const int *proband = INTEGER(Proband);
  const int *father = INTEGER(Father);
  const int *mother = INTEGER(Mother);
  const int *cluster = INTEGER(Cluster);

  int ntrio = LENGTH(Proband);
  
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
  int nsubj, nsnp;
  if (strlen(classS)>5) {
    int *dim = INTEGER(getAttrib(Snps, R_DimSymbol));
    nsubj = dim[0];
    nsnp = dim[1];
  }
  else {
    nsubj = LENGTH(Snps);
    nsnp = 1;
  }
  SEXP snp_names = VECTOR_ELT(getAttrib(Snps, R_DimNamesSymbol), 1);
  index_db name_index = create_name_index(snp_names);

  /* Rules */

  int nrules = 0;
  int pmax = 0;
  GTYPE **gt2ht = NULL;
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
    pmax = *INTEGER(getAttrib(Rules, install("Max.predictors")));
    gt2ht = (GTYPE **)Calloc(pmax, GTYPE *);
    for (int i=0; i<pmax; i++) 
      gt2ht[i] = create_gtype_table(i+1);
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

  /* Check inheritance */

  if (TYPEOF(Check)!=LGLSXP)
    error("Argument error `Check'");
  int check = *LOGICAL(Check);
 
  /* Robust option */

  if (TYPEOF(Robust)!=LGLSXP)    
    error("Argument error `Robust'");
  int robust = *LOGICAL(Robust);

  /* Handling of uncertain genotypes */

  if (TYPEOF(Uncertain) != LGLSXP)
    error("Argument error: Uncertainty is wrong type");
  int uncert = *LOGICAL(Uncertain);
  /* Force robust variance if uncertain option used */
  robust = robust || uncert; 


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
  if (nrules) {
    xadd = (double *) Calloc(nsubj, double);
    xdom = (double *) Calloc(nsubj, double);
  }

  /* Do calculations */
  
  int nonmend=0, Xerrors=0, sexerrors=0;
  
  const double thrsix = 3.0/16.0;
  const double quart = 1.0/4.0;
  const double half = 1.0/2.0;

  for (int t=0; t<ntest; t++) {
    int i = snp_subset? snp_subset[t] - 1: t; 
    const unsigned char *snpsi = NULL;
    double r2 = 0.0;
    if (nrules) {
      if (i >= nrules)
	error("snp_subset out of range");
      SEXP Rule =  VECTOR_ELT(Rules, i);
      if (isNull(Rule)){ /* Monomorphic */
	for (int j=0; j<nsubj; j++){
	  xadd[j] = xdom[j] = 0.0;
	}
      }
      else {
	do_impute(Snps, nsubj, female, NULL, nsubj, name_index, Rule, gt2ht, 
		  xadd, xdom);
	r2 = *REAL(VECTOR_ELT(Rule, 0));
      }
    }
    else {
      if (i >= nsnp)
	error("snp_subset out of range");
      snpsi = snps + nsubj*i;
    }

    /* Initialise score and score variance array */
    
    double u1=0.0, u2=0.0, v11=0.0, v12=0.0, v22=0.0; 
    double up1=0.0, up2=0.0; /* Cluster score contributions */

    int nu = 0;
    int last_clust = 0;
    for (int j=0; j<ntrio; j++) {
      int pj = proband[j] - 1;
      int fj = father[j] - 1;
      int mj = mother[j] - 1;
      int Xmalep = 0;  /* X and male proband */
      int se=0, xe=0;
      if (ifX) {
	Xmalep = !female[pj];
	int se = female[fj] || !female[mj];
	sexerrors += se;
      }   
      if (nrules) {
	double xap = xadd[pj];
	double xaf = xadd[fj];
	double xam = xadd[mj];
	double xdp = xdom[pj];
	double xdf = xdom[fj];
	double xdm = xdom[mj];
	/* Skip if any missing data */
	if (!(ISNA(xap) || ISNA(xaf) || ISNA(xam))) {
	  nu++;
	  if (Xmalep) {
	    up1 += xap - xam;
	    up2 += xdp - xdm;
	  }
	  else {
	    up1 += xap - (xaf+xam)/2.0;
	    up2 += xdp - (xdf+xdm)/2.0;
	  }
	}
      }
      else {
	unsigned char sp = snpsi[pj];
	unsigned char sf = snpsi[fj];
	unsigned char sm = snpsi[mj];
	if ((sp && sf && sm)) { /* Skip if any missing data */
	  int cert = (sp<4 && sf<4 && sm<4);
	  /* Test for mendelian inheritance */
	  int jm = cert && mendelian[sp + 3*sm + 9*sf - 13];
	  nonmend += !jm;
	  if (ifX){
	    int xe = (sp==2 && Xmalep) || (sf==2);
	    Xerrors += xe;
	  }
	  if ((cert|uncert) && (!check||jm) && (!(xe||se))) {
	    double xap, xaf, xam, xdp, xdm, xdf;
	    g2ad(sp, &xap, &xdp);
	    g2ad(sf, &xaf, &xdf);
	    g2ad(sm, &xam, &xdm);
	    nu++;
	    if (Xmalep) {
	      up1 += xap - xam;
	      up2 += xdp - xdm;
	    }
	    else {
	      up1 += xap - (xaf+xam)/2.0;
	      up2 += xdp - (xdf+xdm)/2.0;
	    }
	    if (!robust) { 
	      
	      /* Theoretical variance calculation */
	      
	      if (!Xmalep) {
		if (sm==2) { /* Mother heterozygous */
		  if (sf==2) { /* Father heterozygous */
		      v22 += thrsix;
		      v12 += quart;
		      v11 += half;
		  }
		  else {
		    v11 += quart;
		    if (sf==3) {
		      v22 += quart;
		      v12 += quart;
		    }
		  }
		}
		else if (sf==2) { /* Father heterozygous*/
		  v11 += quart;
		  if (sm==3) {
		    v22 += quart;
		    v12 += quart;
		  }
		}
	      }
	      else if (sm==2) { /* X, heterozygous mother, male offspring */
		v11++;
		v12 += half;
		v22 += quart;
	      }
	    }
	  }
	}
      }

      /* When cluster is complete ... */

      int cj = cluster[j];
      if (cj!=last_clust || j==ntrio) {
	if (robust) {
	  
	  /* Robust variance calculation */
	  
	  v11 += up1*up1;
	  v12 += up1*up2;
	  v22 += up2*up2;
	}

	/* Update scores and zero cluster constributions */

	u1 += up1;
	u2 += up2;
	up1 = up2 = 0.0;
	last_clust = cj;
      }
    }
 
    /* Store results in output object */
    
    Nused[t] = nu;
    if (nrules)
      nr2[t] = nu*r2;
    umat[t] = u1;
    umat[ntest+t] = u2;
    vmat[t] = v11;
    vmat[ntest+t] = v12;
    vmat[2*ntest+t] = v22;
  }
 
   /* Warnings */

  if (Xerrors) 
    warning("%d instances of a male coded as heterozygous at an X locus", 
	    Xerrors);
  if (sexerrors)
    warning("%d instances of sex inconsistent with parental status", 
	    sexerrors);
  if (nonmend)
    warning("%d misinheritances were detected", nonmend);

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
