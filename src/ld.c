#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>
#include <Rdefines.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#define EPS 1.0e-10

/* Declare GSL function */

int gsl_poly_solve_cubic (double a, double b, double c,
                      double *x0, double *x1, double *x2);

/* Phase function */

int phase(const int, const unsigned char *, const unsigned char *, 
	  const int *, double *, double *, double*);

/* Calculate LD stats from phased haplotype frequencies */

void set_arrays(const double *, const double *, double,  double **, int);
  
/* R interface */
 
SEXP ld(SEXP X, SEXP Y, SEXP Depth, SEXP Stats, SEXP Symmetric) {
  char *statnames[7] = {"LLR", "OR", "Q", "Covar", "D.prime", "R.squared", "R"};
  double *arrays[7];

  int depth = *INTEGER(Depth);
  
  /* Stats to calculate */

  int *stats = LOGICAL(Stats);
  if (length(Stats)!=7)
    error("Stats argument");
  int nstats = 0;
  for (int i=0; i<7; i++)
    if (stats[i]) nstats++;

  /* X ---- should be a SnpMatrix or an XSnpMatrix */

  const char *classX = NULL;
  classX = CHAR(STRING_ELT(GET_CLASS(X), 0));
  int *diploid = NULL; /* default Not X */
  if (!strcmp(classX, "XSnpMatrix")) {
    SEXP Diploid = GET_SLOT(X, mkString("diploid"));
    diploid = LOGICAL(Diploid);
  }
  else if (strcmp(classX, "SnpMatrix")) {
    error("Argument error - class(Snps)");
  }
  unsigned char *x = RAW(X);
  int N = nrows(X);
  int MX = ncols(X);

  int XY = 0, MY = 0;
  unsigned char *y;
  if (!isNull(Y)) {
    y= RAW(Y);
    XY = 1;
    const char *classY = NULL;
    classY = CHAR(STRING_ELT(getAttrib(Y, R_ClassSymbol), 0));
    if (!strcmp(classY, "SnpMatrix")) {
      if (diploid)
	error("X and Y are incompatible types");
    }
    else if (!strcmp(classY, "XSnpMatrix")) {
      if (!diploid)
	error("X and Y are incompatible types");
    }
    else {
      error("Argument error - class(Y)");
    }
    if (nrows(Y)!=N)
      error("unequal number of rows in X and Y matrices");
    MY = ncols(Y);
  }
   
  int symmetric = *LOGICAL(Symmetric);
  SEXP Result, Snames;
  int NR = XY? MX*MY: (MX*(MX-1) - (MX-depth)*(MX-depth-1))/2;
  if (nstats>1) {
    PROTECT(Result = allocVector(VECSXP, nstats));
    PROTECT(Snames = allocVector(STRSXP, nstats));
  }
  SEXP Xnames = VECTOR_ELT(GET_DIMNAMES(X), 1);
  if (XY) {
    SEXP LDarrayNames;
    PROTECT(LDarrayNames = allocVector(VECSXP, 2));
    SEXP Ynames = VECTOR_ELT(GET_DIMNAMES(Y), 1);
    SET_VECTOR_ELT(LDarrayNames, 0, Xnames);
    SET_VECTOR_ELT(LDarrayNames, 1, Ynames);
    for (int i=0, is=0; i<7; i++) {
      if (stats[i]) {
	SEXP LDstat;
	PROTECT(LDstat = allocMatrix(REALSXP, MX, MY));
	SET_DIMNAMES(LDstat, LDarrayNames);
	arrays[i] = REAL(LDstat);
	if (nstats>1) {
	  SET_VECTOR_ELT(Result, is, LDstat);
	  SET_STRING_ELT(Snames, is, mkChar(statnames[i]));
	}
	else 
	  Result = LDstat;
	is++;
	UNPROTECT(1); /* LDstat */
      }
      else
	arrays[i] = NULL;
    }
    UNPROTECT(1); /* LDarrayNames */
  }
  else {
    /* Make other arrays for band matrix slots */
    SEXP i_slot, p_slot, Dim_slot, Dimnames_slot, factors_slot, uplo_slot,
      mout_class;
    PROTECT(i_slot = allocVector(INTSXP, NR));
    int *iv = INTEGER(i_slot);
    PROTECT(p_slot = allocVector(INTSXP, MX+1));
    int *pv = INTEGER(p_slot);
    pv[0] = pv[1] = 0;
    for (int j=1, ij=0; j<MX; j++) {
      int ifr = j<depth? 0: j - depth;
      for (int i=ifr; i<j; i++, ij++) 
	iv[ij] = i;
      pv[j+1] = ij;
    }
    PROTECT(Dim_slot = allocVector(INTSXP, 2));
    int *dim = INTEGER(Dim_slot);
    dim[0] = dim[1] = MX;
    PROTECT(Dimnames_slot = allocVector(VECSXP, 2));
    SET_VECTOR_ELT(Dimnames_slot, 0, Xnames);
    SET_VECTOR_ELT(Dimnames_slot, 1, Xnames);
    PROTECT(factors_slot = allocVector(VECSXP, 0));
    PROTECT(uplo_slot = mkString("U"));
    if (symmetric)
      PROTECT(mout_class = MAKE_CLASS("dsCMatrix"));
    else
      PROTECT(mout_class = MAKE_CLASS("dgCMatrix"));
    for (int i=0, is=0; i<7; i++) {
      if (stats[i]){
	SEXP LDstat, x_slot;
	PROTECT(x_slot = allocVector(REALSXP, NR));
	arrays[i] = REAL(x_slot);
	PROTECT(LDstat = NEW_OBJECT(mout_class)); 
	SET_SLOT(LDstat, mkString("x"), x_slot);
	SET_SLOT(LDstat, mkString("i"), i_slot);
	SET_SLOT(LDstat, mkString("p"), p_slot);
	SET_SLOT(LDstat, mkString("Dim"), Dim_slot);
	SET_SLOT(LDstat, mkString("Dimnames"), Dimnames_slot);
	SET_SLOT(LDstat, mkString("factors"), factors_slot);
	if (symmetric)
	  SET_SLOT(LDstat, mkString("uplo"), uplo_slot);
	if (nstats>1) {
	  SET_VECTOR_ELT(Result, is, LDstat);
	  SET_STRING_ELT(Snames, is, mkChar(statnames[i]));
	}
	else 
	  Result = LDstat;
	is++;
	UNPROTECT(2); /* LDstat and x_slot */
      }
    }
    UNPROTECT(7); /* i, p, Dim, Dimnames, factors, uplo, dsCMatrix */
  }

  /* Calculate LD statistics */

  double hapfreqs[4], margins[4], LLR;
  if (XY) {
    double ld[7];
    unsigned char*yj = y;
    for (int j=0, ij=0; j<MY; j++, yj+=N) {
      unsigned char *xi = x;
      for (int i=0; i<MX; i++, xi+=N, ij++) {
	int pr = phase(N, xi, yj, diploid, hapfreqs, margins, &LLR);
	if (pr) {
	  for (int i=0; i<7; i++)
	    if (stats[i]) (arrays[i])[ij] = NA_REAL;
	}
	else 
	  set_arrays(hapfreqs, margins, LLR, arrays, ij);
      }
    }
  }
  else {
    double ld[7];
    unsigned char *xj = x;
    for (int j=1, ij=0; j<MX; j++, xj+=N) {
      int ifr = j<depth? 0: j - depth;
      unsigned char *xi = x+N*ifr; 
      for (int i=ifr; i<j; i++, xi+=N, ij++) {
	int pr = phase(N, xi, xj, diploid, hapfreqs, margins, &LLR);
	if (pr) {
	  for (int i=0; i<7; i++)
	    if (stats[i]) (arrays[i])[ij] = NA_REAL;
	}
	else 
	  set_arrays(hapfreqs, margins, LLR, arrays, ij);
      }
    }
  }
  if (nstats>1) {
    setAttrib(Result, R_NamesSymbol, Snames);
    UNPROTECT(2); /* Result, Snames */
  }
  return Result;
}


/* Function to calculate phased haplotype frequencies of two SNPs */
 
int phase(const int N, const unsigned char *x, const unsigned char *y, 
	  const int *diploid, double *hapfreq, double *margins, double *LLR) {
  int T[4]={0, 0, 0, 0}, G[9]={0, 0, 0, 0, 0, 0, 0, 0, 0};
  /* Defaults (monomorphic SNP) */
  for (int i=0; i<N; i++) {
    unsigned char xi = x[i], yi = y[i];
    if (!xi || !yi || xi>3 || yi>3)
      continue;
    int xyi = (xi-1)*3 + yi-1;
    if (!diploid || diploid[i]) 
      G[xyi]++;
    else 
      switch (xyi) {
      case 0: T[0]++; break;
      case 2: T[1]++; break;
      case 6: T[2]++; break;
      case 8: T[3]++; break;
      default: return 2; /* Heterozygous haploid genotype on X */
      }
  }
  T[0]+= 2*G[0]+G[1]+G[3];
  T[1]+= 2*G[2]+G[1]+G[5];
  T[2]+= 2*G[6]+G[3]+G[7];
  T[3]+= 2*G[8]+G[7]+G[5];
  int Dh = G[4];
  double Nhap = T[0]+T[1]+T[2]+T[3]+2*Dh;
  double E1[4];
  /* Allele frequencies */
  margins[0] = (double)(T[0]+T[1]+Dh)/Nhap;
  margins[1] = 1.0-margins[0];
  margins[2] = (double)(T[0]+T[2]+Dh)/Nhap;
  margins[3] = 1.0-margins[2];
  if (!(margins[0]&&margins[1]&&margins[2]&&margins[3])) {
    return 1; /* At least one SNP is monomorphic */
  }
  /* Solve cubic equation */
  double  roots[3];
  int nroot=1;
  double w2 = (double)(T[0]*T[3]);
  double w3 = (double)(T[1]*T[2]);
  if (!Dh) {
    roots[0] = w2/(w2+w3);
  }
  else {
    double w1 = (double)(T[0]+T[3]-T[1]-T[2])/(double)Dh;
    double Dh2 = Dh*Dh;
    w2 /= Dh2;
    w3 /= Dh2;
    double a = (w1-3.0)/2.0;
    double b = (w2 + w3 - w1 + 1.0)/2.0;
    double c = -w2/2.0;
    nroot = gsl_poly_solve_cubic(a, b, c, roots, roots+1, roots+2);
  }
  double llh, p=-1.0;
  if (LLR || (nroot>1)) {
    for (int i=0; i<nroot; i++) {
      double pi = roots[i];
      if (pi<(-EPS) || pi>(1.0+EPS))
	continue;
      if (pi<0.0)
	pi = 0.0;
      if (pi>1.0)
	pi = 1.0;
      double Dhpi = Dh*pi;
      double Dhqi = Dh-Dhpi;
      E1[0] = ((double)T[0] + Dhpi)/Nhap;
      E1[1] = ((double)T[1] + Dhqi)/Nhap;
      E1[2] = ((double)T[2] + Dhqi)/Nhap;
      E1[3] = ((double)T[3] + Dhpi)/Nhap;
      /* Probability of double het (divided by 2)*/
      double EDh1 = (E1[0]*E1[3] + E1[1]*E1[2])/2.0;
      double llhi = 0.0;
      if (Dh) 
	llhi += Dh*log(EDh1);
      for (int j=0; j<4; j++) 
	if (T[j]) llhi += T[j]*log(E1[j]);
      if (!i || (llhi>llh)) {
	llh = llhi;
	p = pi;
      }
    }
  }
  else 
    p = roots[0];
  if (p==-1.0)
    return 3;
  
  /* Normal return */

  double Dhp = Dh*p;
  double Dhq = Dh-Dhp;
  hapfreq[0] = ((double)T[0] + Dhp)/Nhap;
  hapfreq[1] = ((double)T[1] + Dhq)/Nhap;
  hapfreq[2] = ((double)T[2] + Dhq)/Nhap;
  hapfreq[3] = ((double)T[3] + Dhp)/Nhap;
  if (LLR) {
    double llh0=0.0;
    for (int i=0; i<4; i++)
      llh0 += margins[i]*log(margins[i]);
    *LLR = llh - Nhap*llh0;
  }
  return 0;
}

void set_arrays(const double *hapfreqs, const double *margins, double LLR, 
		double **arrays, int ij) {
  /* LLR */
  if (arrays[0]) (arrays[0])[ij] = LLR;
  /* OR */
  double OR =  hapfreqs[0]*hapfreqs[3]/(hapfreqs[1]*hapfreqs[2]);
  if (arrays[1]) (arrays[1])[ij] = OR;
  /* Yules Q */
  if (arrays[2]) (arrays[2])[ij] = (OR-1.0)/(OR+1.0);
  /* Covariance */
  double covar = hapfreqs[0] - margins[0]*margins[2];
  if (arrays[3]) (arrays[3])[ij] = covar;
  /* D-prime */
  if (arrays[4]) {
    if (covar>0) {
      double P1Q2 = margins[0]*margins[3];
      double P2Q1 = margins[1]*margins[2];
      (arrays[4])[ij] = covar/(P1Q2<P2Q1? P1Q2: P2Q1);
    }
    else {
      double P1Q1 = margins[0]*margins[2];
      double P2Q2 = margins[1]*margins[3];
      (arrays[4])[ij] = -covar/(P1Q1<P2Q2? P1Q1: P2Q2);
    }
  }
  /* R-squared */
  double mprod = margins[0]*margins[1]*margins[2]*margins[3];
  if (arrays[5]) (arrays[5])[ij] = covar*covar/mprod;
  /* R */
  if (arrays[6]) (arrays[6])[ij] = covar/sqrt(mprod);
}
