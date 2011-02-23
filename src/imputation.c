#include <R.h>
#include <Rmath.h>
#include <Rdefines.h>
#include <Rinternals.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "covwin.h"
#include "hash_index.h"
#include "imputation.h"
#include "uncertain.h"
#include "ipf.h"
#include "Rmissing.h"

#define min2(x, y) (x<y? x: y)


int bin_search(const double *sorted, const int len, const double value);
int nearest_N(const double *sorted, const int len, const double value, 
	      const int N);
double covariances(int i, int j, va_list ap);
double snpcov(const unsigned char *x, const unsigned char *y, 
	      const int *female, const int N, const int phase, 
	      const double minA);
double snpmean(const unsigned char *x, 
	       const int *female, const int N);
void utinv(double *, const int);


SEXP snp_impute(const SEXP X, const SEXP Y, const SEXP Xord, const SEXP Yord,
		const SEXP Xpos, const SEXP Ypos, const SEXP Phase, 
		const SEXP Try, const SEXP Stop,
		const SEXP Hapcontr, const SEXP EMcontr,
		const SEXP MinA){

  int try =  *INTEGER(Try);   /* Number to search */
  if (LENGTH(Stop)!=3)
    error("Stop argument not of length 3");
  double r2stop = REAL(Stop)[0]; /* R^2 value to stop inclusion */
  if (r2stop>1.0)
    r2stop = 1.0;
  int pmax = (int) REAL(Stop)[1];   /* Maximum number of predictor variables */
  if (pmax>try)
    pmax = try;
  double dR2 = REAL(Stop)[2]; /* Minimum increase in R^2 to include */
  if (dR2<0)
    dR2 = NA_REAL;
  int phase = *LOGICAL(Phase);   /* haploid or diploid computation */
  double minA = *REAL(MinA);     /* min paired data test */

  /* Fully phased haplotype-based  imputation control args */

  double hapr2 = REAL(Hapcontr)[0]; /* R^2 value to force try */
  /* Required gain in R^2 to persist */
  double hapimp = 0.0;
  if (LENGTH(Hapcontr)>1)
    hapimp = REAL(Hapcontr)[1]; 
  if (LENGTH(EMcontr)!=4)
    error("EMcontr argument not of length 4");
  int maxit_em = REAL(EMcontr)[0];  /* Max EM iterations */
  double tol_em = REAL(EMcontr)[1]; /* EM convergence tolerance */
  int maxit_ipf = REAL(EMcontr)[2]; /* Max iterations for initial IPF */
  int tol_ipf = REAL(EMcontr)[2]; /* Max iterations for initial IPF */


  const double *xpos = REAL(Xpos);  /* Sorted list of X positions */
  int nx = LENGTH(Xpos);      /* Number of X s */
  int *xord = INTEGER(Xord);  /* Corresponding columns in X */
  const double *ypos = REAL(Ypos);  /* Sorted list of Y positions  */
  int ny = LENGTH(Ypos);      /* Number of Y s */
  int *yord = INTEGER(Yord);  /* Corresponding columns in Y */
  int nsubject = nrows(X);
  unsigned char *x = RAW(X);
  unsigned char *y = RAW(Y);
  SEXP Xsnpnames = VECTOR_ELT(getAttrib(X, R_DimNamesSymbol), 1); 
  int *female = NULL;
  SEXP cl = GET_CLASS(X);
  if (TYPEOF(cl) != STRSXP) {
    cl = R_data_class(X, FALSE); /* S4 way of getting class attribute */
  }
  if (!strcmp(CHAR(STRING_ELT(cl, 0)), "XSnpMatrix")) {
    SEXP Female = R_do_slot(X, mkString("Female"));
    female = LOGICAL(Female);
  }

  /* Work arrays */

  double *xy = (double *)Calloc(try, double);       /* XY covariances */
  double *xxd = (double *)Calloc(try, double);      /* Diagonals of XX */
  double *xxi = (double *)Calloc(try*pmax, double); /* Row of XX covars */
  int *sel = (int *)Calloc(pmax, int);              /* Selected SNPs */
  double *coef =(double *)Calloc((pmax*(pmax+1))/2, double);/* Coefficients*/ 
  double *ycoef = (double *)Calloc(pmax, double);   /* Y coefficients */
  COV_WIN_PTR cache = new_window(try, 0);            /* Covariance cache */

  /* Work arrays for haplotype phasing etc. */

  int tmax = (1 << 2*(pmax+1));  /* Space for 4x4x..x4 table */
  int hmax = (1 << (pmax+1));    /* Space for 2x2x..x2 table */
  int *tcell = (int *)Calloc(nsubject, int); /* addresses in table */
  int *contin = (int *)Calloc(tmax, int); 
  int *hcontin=NULL;
  if (female) 
    hcontin = (int *)Calloc(tmax, int);
  double *phap = (double *)Calloc(hmax, double);
  double *phap2 = (double *)Calloc(hmax, double);
  /* gtype->htype lookup tables */
  GTYPE **tables = (GTYPE **)Calloc(pmax+1, GTYPE *);
  for (int i=0; i<=pmax; i++)
    tables[i] = create_gtype_table(i+1);
  /* This is only big enough for first order model */
  unsigned int *llmod = (unsigned int *) Calloc(pmax, unsigned int);

  /* Result */
 
  SEXP Result;
  PROTECT(Result = allocVector(VECSXP, ny));
  setAttrib(Result, R_NamesSymbol, 
	    VECTOR_ELT(getAttrib(Y, R_DimNamesSymbol),1));

  /* Main loop */

  int n_em_fail=0, n_one=0, n_sat=0, n_mod=0, maxpred=0;
  for (int i=0; i<ny; i++) {
    unsigned char *yi = y + nsubject*(yord[i]-1);
    /* Minor allele frequency */
    int ng=0, na=0;
    for (int j=0; j<nsubject; j++) {
      int yij = (int) yi[j];
      if (yij && (yij<4)) {
	ng++;
	na += yij;
      }
    }
    if (ng>0) {
      double maf = (double) (na - ng)/ (double) (2*ng);
      if (maf>0.5)
	maf = 1.0 - maf;
      double yy = snpcov(yi, yi, female, nsubject, phase, minA);
      if (!ISNA(yy)) {
	int start = nearest_N(xpos, nx, ypos[i], try);
	for (int j=0; j<try; j++) { 
	  int jx = nsubject*(xord[start+j]-1);
	  xy[j] = snpcov(x+jx, yi, female, nsubject, phase, minA);
	}
	move_window(cache, start);
	get_diag(cache, xxd, covariances, x, nsubject, xord, female, phase, 
		 minA);
	double resid = yy;
	double rsq = 0.0;
	int nregr=0, ic = 0;
	while (1) {
	  /* Find next snp to be included */
	  double max_due = 0.0;
	  int best = -1;
	  for (int j=0; j<try; j++) {
	    double xxj = xxd[j];
	    double xyj = xy[j];
	    if (xxj==0.0 || ISNA(xxj) || ISNA(xyj))
	      continue;
	    double xyj2 = xyj*xyj;
	    if (xyj2>(xxj*resid)) { /* r^2 > 1 */
	      xy[j] = xyj>0.0? sqrt(xxj*resid): -sqrt(xxj*resid);
	      best = j;
	      max_due = resid;
	    }
	    else {
	      double due = xyj2/xxj;
	      if (due>max_due) {
		max_due = due;
		best = j;
	      }
	    }
	  }
	  if (best<0)
	    break;
	  double *xxin = xxi + try*nregr;
	  get_row(cache, start+best, xxin, 
		  covariances, x, nsubject, xord, female, phase, minA);
	  double *xxik = xxi;
	  int reject = 0;
	  for (int k=0; k<nregr; k++, xxik+=try) {
	    int selk = sel[k];
	    double xxink = xxin[selk], xxikk = xxik[selk];
	    if (ISNA(xxink) || ISNA(xxikk)) {
	      reject = 1;
	      ic -= k;
	      break;
	    }
	    double ck = xxink/xxikk;
	    coef[ic++] = ck;
	    for (int j=0; j<try; j++) {
	      double w1 = xxin[j], w2 = xxik[j];
	      xxin[j] = ISNA(w1) || ISNA(w2)? NA_REAL: w1 - ck*w2;
	    }
	  }
	  if (reject) {
	    xxd[best] = 0.0;
	    continue;
	  }
	  sel[nregr] = best; /* Save index */
	  double bestc = xy[best]/xxd[best]; 
	  ycoef[nregr] = bestc; /* Save regression coefficient */
	  if (ISNA(dR2))
	    dR2 = 2/ (double) (nsubject-nregr-2);
	  double deltaR2 =  max_due/resid;
	  resid -= max_due;
	  rsq = 1.0 - resid/yy; 
	  nregr++;
	  int stop = (rsq>=r2stop)||(nregr==pmax)||(deltaR2<dR2); 
	  if (stop) 
	    break;

	  double vn = xxd[best];
	  xxd[best] = 0.0;
	  for (int j=0; j<try; j++) {
	    double w = xxin[j];
	    if (!ISNA(w)) {
	      xy[j] -=  bestc*w;
	      w = xxd[j]-w*w/vn;
	      xxd[j] = w>0.0? w: 0.0;
	    }
	  }
	}

	if (nregr>0) {	

	  /* Unphased multilocus genotype as 4x4x... table */

	  if (nregr>maxpred)
	    maxpred = nregr;
	  for (int j=0; j<nsubject; j++) {
	    unsigned char yij = yi[j];
	    if (yij>4)
	      yij = 0;
	    tcell[j] = (int) yij;
	  }
	  for (int k=0, sh=2; k<nregr; k++, sh+=2) {
	    unsigned char *xk = x + nsubject*(xord[start+sel[k]]-1);
	    for (int j=0; j<nsubject; j++) {
	      unsigned int xkj = (unsigned int) xk[j];
	      if (xkj>4)
		xkj = 0;
	      tcell[j] = tcell[j] | (xkj << sh);
	    }
	  }
	  int dim = nregr+1;
	  int tsize = 1<<(2*dim);
	  memset(contin, 0x00, tsize*sizeof(int));
	  if (hcontin)
	    memset(hcontin, 0x00, tsize*sizeof(int));
	  for (int j=0; j<nsubject; j++) {
	    if (female && !female[j])
	      hcontin[tcell[j]]++;
	    else
	      contin[tcell[j]]++;
	  }

	  /* Phased haplotype frequencies by EM and IPF algorithms */

	  phap[0] = -1.0; 
	  double r2_mod=-1.0, r2_sat=-1.0, r2=-1.0;
	  double *whichp = NULL;

	  /* Phase using saturated model */
	  
	  int em_fail_sat = emhap(dim, contin, hcontin, tables[nregr], 
			  maxit_em, tol_em, phap, 0, NULL); 
	  if (em_fail_sat>=0) { 
	    r2_sat = hap_r2(nregr, phap);
	    if (nregr == 1) {
	      n_one++;
	      r2 = r2_sat;
	      whichp = phap;
	    }
	    else if (hapr2 < 1.0) {
	    
	      /* First order log-linear model */

	      unsigned int bit = 0x02;
	      unsigned int all = 0x00;
	      for (int i=0; i<nregr; i++, bit=bit<<1) {
		all = all | bit;
		llmod[i] = bit & 0x01;
	      }
	      /* Get starting value by IPF on haplotype frequencies */
	      phap2[0] = -1;
	      ipf(nregr+1, phap, nregr, llmod, phap2, maxit_ipf, tol_ipf);
	      int em_fail = emhap(dim, contin, hcontin, tables[nregr], 
			      maxit_em, tol_em, phap2, 1+nregr, llmod);
	      if (em_fail>=0) 
		r2_mod = hap_r2(nregr, phap2);	  
	    }
	    
	    /* Only use saturated model if enough  improvedment */

	    if (r2_mod<0.0 || (r2_sat-r2_mod)>hapimp) {
	      whichp = phap;
	      r2 = r2_sat;
	      n_sat++;
	    }
	    else { 
	      whichp = phap2;
	      r2 = r2_mod;
	      n_mod++;
	    }

	    /* Save imputation rule */

	    if (r2>0.0) {
	      SEXP Rule, Rlnames, Maf, R2, Pnames, Coefs;
	      PROTECT(Rule = allocVector(VECSXP, 4));
	    
	      PROTECT(Rlnames = allocVector(STRSXP, 4));
	      SET_STRING_ELT(Rlnames, 0, mkChar("maf"));
	      SET_STRING_ELT(Rlnames, 1, mkChar("r.squared"));
	      SET_STRING_ELT(Rlnames, 2, mkChar("snps"));
	      SET_STRING_ELT(Rlnames, 3, mkChar("hap.probs"));
	      
	      PROTECT(Maf = allocVector(REALSXP, 1));
	      *REAL(Maf) = maf;
	      PROTECT(R2 = allocVector(REALSXP, 1));
	      PROTECT(Pnames = allocVector(STRSXP, nregr));
	      for (int j=0; j<nregr; j++) {
	      int xsnp =  xord[start+sel[j]]-1;
	      SET_STRING_ELT(Pnames, j, STRING_ELT(Xsnpnames, xsnp));
	      }
	      *REAL(R2) = r2;
	      int lenp = 1 << (nregr+1);
	      PROTECT(Coefs = allocVector(REALSXP, lenp));
	      double *coefs = REAL(Coefs);
	      for (int j=0; j<lenp; j++) 
		coefs[j] = whichp[j];

	      SET_VECTOR_ELT(Rule, 0, Maf);
	      SET_VECTOR_ELT(Rule, 1, R2);
	      SET_VECTOR_ELT(Rule, 2, Pnames);
	      SET_VECTOR_ELT(Rule, 3, Coefs);
	      setAttrib(Rule, R_NamesSymbol, Rlnames);
	      SET_VECTOR_ELT(Result, yord[i]-1, Rule);
	      UNPROTECT(6);
	    }
	    else {
	      /* Failed prediction algorithm */
	      SET_VECTOR_ELT(Result, yord[i]-1, R_NilValue);
	    }
	  }
	  else {
	    warning("phasing algorithm failed");
	    n_em_fail ++; 
	  }
	}
	else {
	  /* No valid predictor */
	  SET_VECTOR_ELT(Result, yord[i]-1, R_NilValue);
	}
      }
      else {
	/*MAF too low  */
	SET_VECTOR_ELT(Result, yord[i]-1, R_NilValue);
      }
    }
    else {
      /* No data */
      SET_VECTOR_ELT(Result, yord[i]-1, R_NilValue);
    }
  }
  if (n_em_fail)
    Rprintf("Failures of phasing algorithm: %d\n", n_em_fail);
  if (n_one)
    Rprintf("SNPs tagged by a single SNP: %d\n", n_one);
  n_sat -= n_one;
  if (n_mod)
    Rprintf("SNPs tagged by multiple tag haplotypes (log-linear model): %d\n", 
	    n_mod);
 if (n_sat)
    Rprintf("SNPs tagged by multiple tag haplotypes (saturated model): %d\n", 
	    n_sat);
 
  SEXP IrClass, Package;
  PROTECT(IrClass = allocVector(STRSXP, 1));
  SET_STRING_ELT(IrClass, 0, mkChar("ImputationRules"));
  PROTECT(Package = allocVector(STRSXP, 1));
  SET_STRING_ELT(Package, 0, mkChar("snpStats"));
  setAttrib(IrClass, install("package"), Package);
  classgets(Result, IrClass);
  SEXP Maxpred;
  PROTECT(Maxpred = allocVector(INTSXP, 1));
  INTEGER(Maxpred)[0] = maxpred;
  setAttrib(Result, install("Max.predictors"), Maxpred);
  SET_S4_OBJECT(Result);

  /* Tidy up */

  Free(xy);
  Free(xxi);
  Free(xxd);
  Free(sel);
  Free(coef);
  Free(ycoef);
  free_window(cache);
  Free(contin);
  if (hcontin)
    Free(hcontin);
  Free(phap);
  Free(phap2);
  Free(tcell);
  for (int i=0; i<=pmax; i++)
    destroy_gtype_table(tables[i], i+1);
  Free(tables);
  Free(llmod);

  UNPROTECT(4);
  return Result;
}

double covariances(int i, int j, va_list ap) {
  unsigned char *snps = va_arg(ap, unsigned char *);
  int N = va_arg(ap, int);
  int *cols = va_arg(ap, int *);
  int *female = va_arg(ap, int *);
  int phase = va_arg(ap, int); 
  double minA = va_arg(ap, double);
  int ik = N*(cols[i]-1), jk = N*(cols[j]-1);
  return snpcov(snps+ik, snps+jk, female, N, phase, minA);
}


double snpcov(const unsigned char *x, const unsigned char *y, 
	  const int *female, const int N, const int phase, const double minA) {
  int n1=0, n2=0, nt=0, sx=0, sy=0, sxy=0;
  double cov, n11;
  if (phase) {
    if (female)
      error("phase=TRUE not yet implemented for the X chromosome");
    error("phase=TRUE not yet implemented");
    return NA_REAL;
  }
  else {
    if (female) {
      for (int k=0; k<N; k++) {
	int xk = (int) *(x++);
	int yk = (int) *(y++);
	if (xk && (xk<4) && yk && (yk<4)) { 
	  xk--;
	  yk--;
	  if (female[k]) {
	    n2++;
	  }
	  else {
	    n1++;
	    xk/=2;
	    yk/=2;
	  }
	  sx += xk;
	  sy += yk;
	  sxy += xk*yk;
	}
      }
      nt = 2*n2 + n1;
      if (nt<2)
	return NA_REAL;
      int nt_1 = nt-1;
      double p2 = (double)(2*n2)/(double)nt;
      double ps = (double)sx * (double)sy;
      cov = ((double)sxy - (1.0+p2)*ps/(double)nt)/
	((double)nt_1 - p2);
      n11 = (double)nt_1*(sxy - p2*ps/(double)nt_1)/((double)nt_1-p2);
    }
    else {
      for (int k=0; k<N; k++) {
	int xk = (int) *(x++);
	int yk = (int) *(y++);
	if (xk && (xk<4) && yk && (yk<4)) {
	  xk--;
	  yk--;
	  n2++;
	  sx += xk;
	  sy += yk;
	  sxy += xk*yk;
	}
      }
      if (n2 < 2)
	return NA_REAL;
      nt = 2*n2;
      double ps = (double)sx * (double)sy;
      double n_1 = n2 - 1;
      cov = 0.5*((double)sxy - ps/(double)n2)/(double)n_1;
      double twon_1 = nt - 1;
      n11 = (double)twon_1*((double)sxy - ps/(double)twon_1)/
	(2.0*(double)n_1);
    }
    double test = (cov > 0.0)?
      min2(n11, nt - sx - sy + n11):
      min2(sx - n11, sy - n11);
    /*    printf("n11 = %lf, test = %lf, ", n11, test); */
    if (test < minA)
      return NA_REAL;
    return cov;
  }
}

/*  Routine for testing snpcov */

SEXP snpcov_test(const SEXP X, const SEXP i, const SEXP j, const SEXP minA) {
  int ii = *INTEGER(i) - 1;
  int jj = *INTEGER(j) - 1;
  int N = nrows(X);
  double ma = *REAL(minA);
  unsigned char *x = RAW(X);
  double mycov = snpcov(x+N*ii, x+N*jj, 0, N, 0, ma);
  Rprintf("N = %d, cov = ", N);
  if (ISNA(mycov))
    Rprintf("NA_REAL\n");
  else
    Rprintf("%lf\n", mycov);
  SEXP Result = allocVector(REALSXP, 1);
  *REAL(Result) = mycov;
  return Result;
}
 
double snpmean(const unsigned char *x, const int *female, const int N) {
  int sum=0, sumx=0;
  if (female) {
    for (int i=0; i<N; i++) {
      int wt = female[i]? 2: 1;
      int w = (int) *(x++);
      if (w && (w<4)) {
	sum += wt;
	sumx += wt*w;
      }
    }
  }
  else {
    for (int i=0; i<N; i++) {
      int w = (int) *(x++);
      if (w && (w<4)) {
	sum++;
	sumx += w;
      }
    }
  }
  if (sum)
    return (double) sumx/(double) sum - 1.0;
  else
    return NA_REAL;
}

/* Inverse of unit triangular matrix -- diagonal not stored */
void utinv(double *mat, const int N){
  if (N<2)
    return;
  for (int j=1, ij=0; j<N; j++) {
    for (int i=0, is=0; i<j; i++, ij++) {
      double w = mat[ij];
      if (ISNA(w))
	warning("Bug: NAs in triangular coefficients matrix");
      for (int k=(i+1), k1=ij+1, k2=is; k<j; k++){ 
	w += mat[k1]*mat[k2];
	k1++;
	k2+= (k+1);
      }
      mat[ij] = (-w);
      is += (i+2);
    }
  }
}
  

int bin_search(const double *sorted, const int len, const double value) {
  int low = 0, high = len - 1;
  int mid = (low + high)/2;
  while (low<mid) {
    if (sorted[mid] > value)
      high = mid;
    else if (sorted[mid] < value)
      low = mid;
    else
      return mid;
    mid = (low + high)/2;
  }
  return low;
}

int nearest_N(const double *sorted, const int len, const double value, 
	      const int N) {
  int last = len - N;
  int res = bin_search(sorted, len, value) - N/2;
  if (res<0)
    res = 0;
  if (res>last)
    res = last;
  if ((value - sorted[res])>(sorted[res+N-1] - value)) {
    while (res<last) {
      res++;
      if ((value - sorted[res])<=(sorted[res+N-1] - value))
	return res;
    }
  }
  else {
    while (res>0) {
      res--;
      if ((value - sorted[res])>=(sorted[res+N-1] - value))
	return res;
    }
  }
  return res;
}

/* Create a hash index from a set of names */

index_db create_name_index(const SEXP names) {
  if (TYPEOF(names)!=STRSXP) 
    error("Names not character variable");
  int N = LENGTH(names);
  index_db res = index_create(N);
  for (int i=0; i<N; i++) {
    if (index_insert(res, CHAR(STRING_ELT(names, i)), i)!=0)
      error("Duplicate names");
  }
  return res;
}

/* Do an imputation (on selected rows) */

void do_impute(const SEXP Obs_snps, const int nrow, 
	       const int *female,
	       const int *rows, int nuse, 
	       index_db snp_names,
	       SEXP Rule, GTYPE **gt2ht, 
	       double *value_a, double *value_d) {
  unsigned char *snps = RAW(Obs_snps);
  SEXP Snps = VECTOR_ELT(Rule, 2);
  int nsnp = LENGTH(Snps);
  SEXP Coefs = VECTOR_ELT(Rule, 3);
  int ncoefs = LENGTH(Coefs);
  double *coefs = REAL(Coefs);
  if (!rows)
    nuse = nrow;

  if (ncoefs==(nsnp+1))  /* Regression imputation */
    error("Old imputation rule; not supported by this version");
  else { /* Imputation from phased haplotypes */
    int *gt = (int *)Calloc(nuse, int);
    int *fem = NULL;
    if (female)
      fem = (int *)Calloc(nuse, int);
    memset(gt, 0x00, nuse*sizeof(int));
    /* Calculate predictor genotypes */
    for (int j=0, sh=0; j<nsnp; j++, sh+=2) {
      int jj = index_lookup(snp_names, CHAR(STRING_ELT(Snps, j)));
      if (jj<0)
	error("Couldn't match snp name: %s", CHAR(STRING_ELT(Snps, j)));
      for (int r=0, ist=nrow*jj; r<nuse; r++) {
	int i = rows? rows[r]-1: r;
	int sij = (int)snps[ist+i];
	gt[r] = gt[r] | (sij << sh);
	if (fem)
	  fem[r] = female[i];
      }
    }
    /* 
       Score genotypes
       Perhaps more efficient to compute a lookup table at outset
    */

    const GTYPE *gtab = gt2ht[nsnp-1];
    for (int i=0; i<nuse; i++) {
      double posterior[3];
      int gti = gt[i];
      if (gti) {
	int mX = fem? (!fem[i]): 0; /* X and male? */
	predict_gt(nsnp, gti, mX, coefs, gtab, posterior);
	int ispna = ISNA(posterior[0]);
	value_a[i] = ispna? NA_REAL: posterior[1]+2.0*posterior[2];
	if (value_d)
	  value_d[i] = ispna? NA_REAL: posterior[2];
      }
      else {
	value_a[i] = NA_REAL;
	if (value_d)
	  value_d[i] = NA_REAL;
      }
    }
    Free(gt);
    if (fem)
      Free(fem);
  }
}
  
SEXP impute_snps(const SEXP Rules, const SEXP Snps, const SEXP Subset, 
		 const SEXP As_numeric) { 
  int *female_in=NULL, *female=NULL;
  SEXP cl = GET_CLASS(Snps);
  if (TYPEOF(cl) != STRSXP) {
    cl = R_data_class(Snps, FALSE); /* S4 way of getting class attribute */
  }
  if (!strcmp(CHAR(STRING_ELT(cl, 0)), "XSnpMatrix")) {
    SEXP Female = R_do_slot(Snps, mkString("Female"));
    female_in = LOGICAL(Female);
  }

  SEXP names = getAttrib(Snps, R_DimNamesSymbol);
  index_db name_index = create_name_index(VECTOR_ELT(names, 1));
  int N = nrows(Snps);
    int M = LENGTH(Rules);
  int *subset = NULL;
  int nsubj = N;
  SEXPTYPE sutype = TYPEOF(Subset);
  if (sutype==INTSXP) { 
    if (LENGTH(Subset)>N)
      error("Dimension error - Subset");
    subset = INTEGER(Subset);
    nsubj = LENGTH(Subset);
  }
  else if (sutype!=NILSXP)
    error("Argument error - Subset");


  double *w1 = (double *)Calloc(nsubj, double);
  double *w2 = (double *)Calloc(nsubj, double);
  SEXP Result, Female, Dimnames, Class, Package;
  double *dresult = NULL;
  unsigned char *rresult = NULL;

  int as_numeric = *LOGICAL(As_numeric);
  if (as_numeric) {
    PROTECT(Result = allocMatrix(REALSXP, nsubj, M));
    dresult = REAL(Result);
  }
  else {
    PROTECT(Result = allocMatrix(RAWSXP, nsubj, M));
    rresult = RAW(Result);
    PROTECT(Class = allocVector(STRSXP, 1));
    if (female_in) {
      PROTECT(Female = allocVector(REALSXP, nsubj));
      R_do_slot_assign(Result, mkString("Female"), Female);
      SET_STRING_ELT(Class, 0, mkChar("XSnpMatrix"));
      female = LOGICAL(Female);
    }
    else {
      SET_STRING_ELT(Class, 0, mkChar("SnpMatrix"));
    }
    PROTECT(Package = allocVector(STRSXP, 1));
    SET_STRING_ELT(Package, 0, mkChar("snpStats"));
    setAttrib(Class, install("package"), Package);
    classgets(Result, Class);
    SET_S4_OBJECT(Result);
  }
  PROTECT(Dimnames = allocVector(VECSXP, 2));
  SET_VECTOR_ELT(Dimnames, 0, VECTOR_ELT(names, 0));
  SET_VECTOR_ELT(Dimnames, 1, getAttrib(Rules, R_NamesSymbol));
  setAttrib(Result, R_DimNamesSymbol, Dimnames);
  int pmax = *INTEGER(getAttrib(Rules, install("Max.predictors")));
  GTYPE **gt2ht = (GTYPE **)Calloc(pmax, GTYPE *); 
  for (int i=0; i<pmax; i++)
    gt2ht[i] = create_gtype_table(i+1);
  for (int j=0, ji=0; j<M; j++) {
    SEXP Rule = VECTOR_ELT(Rules, j);
    if (isNull(Rule)) {
      if (as_numeric) {
	for (int i=0; i<nsubj; i++, ji++)
	  dresult[ji] = NA_REAL;
      }
      else {
	for (int i=0; i<nsubj; i++, ji++)
	  rresult[ji] = 0;
      }
    }
    else {
      do_impute(Snps, N, female_in, subset, nsubj, name_index, Rule, gt2ht, 
		w1, w2); 
      if (as_numeric) {
	for (int i=0; i<nsubj; i++, ji++){
	  dresult[ji] = w1[i];
	}
      }
      else {
	for (int i=0; i<nsubj; i++, ji++) {
	  double w1i = w1[i];
	  double w2i = w2[i];
	  rresult[ji] = ISNA(w1i)? 0: post2g(w1i-2.0*w2i, w2i);
	}
	if (female) {
	  for (int i=0; i<nsubj; i++) {
	    int ii = subset? subset[i]-1: i;
	    female[i] = female_in[ii];
	  }
	}
      }
    }
  }
  
  index_destroy(name_index);
  for (int i=0; i<pmax; i++) 
    destroy_gtype_table(gt2ht[i], i+1);
  Free(gt2ht);
  if (as_numeric)
    UNPROTECT(2);
  else {
    if (female_in)
      UNPROTECT(5);
    else 
      UNPROTECT(4);
  } 
  Free(w1);
  Free(w2);
  return Result;
}

/* Extract r-square and  number of tag snps */

SEXP r2_impute(const SEXP Rules) {
  int M = LENGTH(Rules);
  SEXP Result;
  PROTECT(Result = allocMatrix(REALSXP, M, 2));
  double *result = REAL(Result);
  for (int i=0; i<M; i++) {
    SEXP Rule = VECTOR_ELT(Rules, i);
    if (TYPEOF(Rule)==NILSXP){
      result[i] = NA_REAL;
      result[i+M] = NA_REAL;
    }
    else {
      result[i] = *REAL(VECTOR_ELT(Rule, 1));
      result[i+M] = LENGTH(VECTOR_ELT(Rule, 2));
    }
  }
  UNPROTECT(1);
  return Result;
}

