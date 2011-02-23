#include <stdio.h>
#include <ctype.h>
#include <R.h>
#include <Rinternals.h>
#include "uncertain.h"
#include "zlib.h"

#define BUFFERSIZE 126
/* wc function */

void gzwc(const gzFile infile, const int nline, 
	  int *chars, int *words, int *lines){
  int ch;
  int sp=1;

  *chars = 0;
  *words = 0;
  *lines = 0;
  while((ch=gzgetc(infile))!=EOF && (!nline || (*lines)<nline)) {
    ++(*chars);
    if(isspace(ch)) 
      sp=1;
    else if(sp) {
      ++(*words);
      sp=0;
    }
    if (ch=='\n') 
      ++(*lines);
  }
  gzrewind(infile);
  return;
}

/* Next white-space delimited field */

void gznext(const gzFile infile, char *buffer, const int len) {
  char ch;
  while (isspace(ch=gzgetc(infile))){} /* Skip leading white space */
  int i = 0, maxi = len-2;
  do {
    if (i>maxi)
      error("input field exceeds buffer length");
    buffer[i] = ch;
    i++;
  } while (!(isspace(ch=gzgetc(infile))));
  buffer[i] = 00;
}

/* MACH MLPROB file */

SEXP read_mach(const SEXP Filename, const SEXP Colnames, const SEXP Nsubject) {
  int nsubj=0;
  SEXPTYPE tns = TYPEOF(Nsubject);
  if (tns!=NILSXP) {
    if (tns==INTSXP)
      nsubj = *INTEGER(Nsubject);
    else if (tns==REALSXP)
      nsubj = (int) *REAL(Nsubject);
    else
      error("illegal type for nrow argument");
  }
  if (TYPEOF(Filename)!=STRSXP || length(Filename)>1)
    error("Argument type error: Filename");
  const char *filename = CHAR(STRING_ELT(Filename, 0));
  Rprintf("Reading MACH data from file %s\n", filename);
  gzFile infile =  gzopen(filename, "rb");
  if (!infile)
    error("Could not open input file");
  int chars, words, lines, ncol;
  if (!nsubj) {
    gzwc(infile, 0, &chars, &words, &lines);
    if (words%lines)
      error("Number of fields is not a multiple of number of lines");
    ncol = words/lines - 2;
  }
  else {
    gzwc(infile, 1,  &chars, &words, &lines);
    lines = nsubj;
    ncol = words -2;
  }
  if (ncol<=0)
    error("No loci to read");
  if (ncol%2)
    error("Odd number of fields");
  ncol/=2;
  if (TYPEOF(Colnames)!=NILSXP) {
    if (TYPEOF(Colnames)!=STRSXP)
      error("column names are not of character type");
    if (length(Colnames)!=ncol) 
      error("Number of entries on file does not correspond with column names");
  }

  Rprintf("Reading SnpMatrix with %d rows and %d columns\n", lines, ncol);

  /* Build output object */ 

  SEXP Result, Dimnames, Rnames=R_NilValue, Package, Class;
  PROTECT(Result = allocMatrix(RAWSXP, lines, ncol));
  unsigned char *result = RAW(Result);
  memset(result, 0x00, lines*ncol);
  PROTECT(Dimnames = allocVector(VECSXP, 2));
  PROTECT(Rnames = allocVector(STRSXP, lines));
  SET_VECTOR_ELT(Dimnames, 0, Rnames);
  if (TYPEOF(Colnames)!=NILSXP) {
    SET_VECTOR_ELT(Dimnames, 1, Colnames);
  }
  else {
    SEXP Cnames;
    PROTECT(Cnames = allocVector(STRSXP, ncol));
    char id[BUFFERSIZE];
    for (int i=0; i<ncol; i++) {
      sprintf(id, "SNP%d", i+1);
      SET_STRING_ELT(Cnames, i, mkChar(id));
    }
    SET_VECTOR_ELT(Dimnames, 1, Cnames);
    UNPROTECT(1); /* Cnames should be protected via Dimnames */
  }
  setAttrib(Result, R_DimNamesSymbol, Dimnames);
  UNPROTECT(2); /* Dimnames and Rnames should be protected via Result */

  /* Class */

  PROTECT(Class = allocVector(STRSXP, 1));
  SET_STRING_ELT(Class, 0, mkChar("SnpMatrix"));
  PROTECT(Package = allocVector(STRSXP, 1));
  SET_STRING_ELT(Package, 0, mkChar("snpStats"));
  setAttrib(Class, install("package"), Package);
  classgets(Result, Class);
  SET_S4_OBJECT(Result);
  UNPROTECT(2);
 
  /* Read in data */

  char buffer[BUFFERSIZE];
  for (int i=0; i<lines; i++) {
    gznext(infile, buffer, BUFFERSIZE);
    SET_STRING_ELT(Rnames, i, mkChar(buffer));
    gznext(infile, buffer, BUFFERSIZE);
    if (strcmp(buffer, "ML_PROB")!=0)
      error("file does not appear to be an MLPROB file (field 2=%s)", buffer);
    for (int j=0, ij=i; j<ncol; j++, ij+=lines) {
      double pAA, pAB;
      gznext(infile, buffer, BUFFERSIZE);
      if (sscanf(buffer, "%lf", &pAA)!=1)
	error("read error at line %d, SNP %d: %s", i, j, buffer);
      gznext(infile, buffer, BUFFERSIZE);
      if (sscanf(buffer, "%lf", &pAB)!=1)
	error("read error at line %d, SNP %d: %s", i, j, buffer);
      double pBB = 1.0 - pAA - pAB;
      /* Deal with rounding error */
      if (pBB<0.0){
	pBB = 0.0;
	double tot = pAA+pAB;
	pAA/=tot;
	pAB/=tot;
      }
      result[ij] = post2g(pAB, pBB);
    }
  }
  UNPROTECT(1);
  return(Result);
}

/* IMPUTE or BEAGLE output */

SEXP read_impute(const SEXP Filename, const SEXP Rownames, const SEXP Nsnp, 
		 const SEXP Snpcol, const SEXP Header) {
  int nsnp=0;
  SEXPTYPE tns = TYPEOF(Nsnp);
  if (tns!=NILSXP) {
    if (tns==INTSXP)
      nsnp = *INTEGER(Nsnp);
    else if (tns==REALSXP)
      nsnp = (int) *REAL(Nsnp);
    else
      error("illegal type for nsnp argument");
  }
  /* For impute which snp id used. zero or NULL means BEAGLE output */ 
  int snpcol=0; 
  tns = TYPEOF(Snpcol);
  if (tns!=NILSXP) {
    if (tns==INTSXP)
      snpcol = *INTEGER(Snpcol);
    else if (tns==REALSXP)
      snpcol = (int) *REAL(Snpcol);
    else
      error("illegal type for snpcol argument");
  }
  if (snpcol<0 || snpcol>2)
    error("illegal snpcol argument");
  int ncol_skip = snpcol? 5: 3;
      
  if (TYPEOF(Header) != LGLSXP) 
    error("illegal header argument");
  int header = *LOGICAL(Header);

  if (TYPEOF(Filename)!=STRSXP || length(Filename)>1)
    error("Argument type error: Filename");
  const char *filename = CHAR(STRING_ELT(Filename, 0));
  Rprintf("Reading IMPUTE data from file %s\n", filename);
  gzFile infile =  gzopen(filename, "rb");
  if (!infile)
    error("Could not open input file");
  int chars, words, lines, N;
  if (!nsnp) {
    gzwc(infile, 0, &chars, &words, &lines);
    if (words%lines)
      error("Number of fields is not a multiple of number of lines");
    N = words/lines - ncol_skip;
    nsnp = lines;
  }
  else {
    gzwc(infile, 1,  &chars, &words, &lines);
    N = words - ncol_skip;
  }
  if (N<=0)
    error("No loci to read");
  if (N%3)
    error("Number of probabilities is not a multiple of 3");
  N/=3;
  int norownames = 1;
  if (TYPEOF(Rownames)!=NILSXP) {
    norownames = 0;
    if (TYPEOF(Rownames)!=STRSXP)
      error("row names are not of character type");
    if (length(Rownames)!=N) 
      error("Number of entries on file does not correspond with row names");
  }

  Rprintf("Reading SnpMatrix with %d rows and %d columns\n", N, nsnp);

  /* Build output object */ 

  SEXP Result, Dimnames, Colnames, Package, Class, Rnames=R_NilValue;
  PROTECT(Result = allocMatrix(RAWSXP, N, nsnp));
  unsigned char *result = RAW(Result);
  memset(result, 0x00, N*nsnp);
  PROTECT(Dimnames = allocVector(VECSXP, 2));
  PROTECT(Colnames = allocVector(STRSXP, nsnp));
  SET_VECTOR_ELT(Dimnames, 1, Colnames);
  if (TYPEOF(Rownames)!=NILSXP) {
    SET_VECTOR_ELT(Dimnames, 0, Rownames);
  }
  else {
    PROTECT(Rnames = allocVector(STRSXP, N));
    char id[BUFFERSIZE];
    if (!header) {
      for (int i=0; i<N; i++) {
	sprintf(id, "Sample%d", i+1);
	SET_STRING_ELT(Rnames, i, mkChar(id));
      }
    }
    SET_VECTOR_ELT(Dimnames, 0, Rnames);
    UNPROTECT(1); /* Rnames should be protected via Dimnames */
  }
  setAttrib(Result, R_DimNamesSymbol, Dimnames);
  UNPROTECT(2); /* Dimnames and Colnames should be protected via Result */

  /* Class */

  PROTECT(Class = allocVector(STRSXP, 1));
  SET_STRING_ELT(Class, 0, mkChar("SnpMatrix"));
  PROTECT(Package = allocVector(STRSXP, 1));
  SET_STRING_ELT(Package, 0, mkChar("snpStats"));
  setAttrib(Class, install("package"), Package);
  classgets(Result, Class);
  SET_S4_OBJECT(Result);
  UNPROTECT(2);

  char buffer[BUFFERSIZE];

  /* If BEAGLE, read in header line */

  if (!snpcol && header) {
    gznext(infile, buffer, BUFFERSIZE);
    if (strcmp(buffer, "marker"))
      error("Header line not compatible with BEAGLE output format");
    gznext(infile, buffer, BUFFERSIZE);
    gznext(infile, buffer, BUFFERSIZE);
    for (int i=0; i<N; i++) {
      gznext(infile, buffer, BUFFERSIZE);
      if (norownames)
	SET_STRING_ELT(Rnames, i, mkChar(buffer));
      gznext(infile, buffer, BUFFERSIZE);
      gznext(infile, buffer, BUFFERSIZE);
    }
  }
  if (snpcol)
    snpcol--;

  /* Read in data */

  for (int i=0, ij=0; i<nsnp; i++) {
    for (int j=0; j<ncol_skip; j++) {
      gznext(infile, buffer, BUFFERSIZE);
      if (j==snpcol) {
	SET_STRING_ELT(Colnames, i, mkChar(buffer));
      }
    }
    for (int j=0; j<N; j++, ij++) {
      double pAA, pAB, pBB;
      gznext(infile, buffer, BUFFERSIZE);
      if (sscanf(buffer, "%lf", &pAA)!=1)
	error("read error at line %d, sample %d: %s", i, j, buffer);
      gznext(infile, buffer, BUFFERSIZE);
      if (sscanf(buffer, "%lf", &pAB)!=1)
	error("read error at line %d, sample %d: %s", i, j, buffer);
      gznext(infile, buffer, BUFFERSIZE);
      if (sscanf(buffer, "%lf", &pBB)!=1)
	error("read error at line %d, sample %d: %s", i, j, buffer);
      double psum = pAA+pAB+pBB;
      if (psum>0.0) {
	pAB/=psum;
	pBB/=psum;
	result[ij] = post2g(pAB, pBB);
	
      }
      else { 
	result[ij] = 0x00;
      }
    }
  }
  UNPROTECT(1);
  return(Result);
}


