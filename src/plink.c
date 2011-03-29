/* 
   Read a plink .bed file as a SnpMatrix
 
*/


#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
/* #include "Rmissing.h" */


/* this may not be fast but seems reliable */

void skip(FILE *in, int lines, int size) {
  if (lines) {
    for (int i=0; i<lines; i++)
      for (int j=0; j<size; j++)
	if (fgetc(in)==EOF) error("unexpected end of file");
  }    
}


SEXP readbed(SEXP Bed, SEXP Id, SEXP Snps, SEXP Rsel, SEXP Csel) {
  const unsigned char recode[4] = {'\x01', '\x00', '\x02', '\x03'};
  const unsigned char mask = '\x03';
  int nrow = LENGTH(Id);
  int ncol = LENGTH(Snps);
  const char *file = CHAR(STRING_ELT(Bed, 0));
  FILE *in = fopen(file, "rb");
  if (!in)
    error("Couln't open input file: %s", file);
  unsigned char start[3];
  if (fread(start, 1, 3, in)!=3)
    error("Failed to read first 3 bytes");
  if (start[0]!='\x6C' || start[1]!='\x1B')
    error("Input file does not appear to be a .bed file (%X, %X)", 
	  start[0], start[1]);

  /* Create output object */

  SEXP Result, Dimnames, Package, Class;
  PROTECT(Result = allocMatrix(RAWSXP, nrow, ncol));
  PROTECT(Dimnames = allocVector(VECSXP, 2));
  SET_VECTOR_ELT(Dimnames, 0, Id);
  SET_VECTOR_ELT(Dimnames, 1, Snps);
  setAttrib(Result, R_DimNamesSymbol, Dimnames);
  PROTECT(Class = allocVector(STRSXP, 1));
  SET_STRING_ELT(Class, 0, mkChar("SnpMatrix"));
  PROTECT(Package = allocVector(STRSXP, 1));
  SET_STRING_ELT(Package, 0, mkChar("snpStats"));
  setAttrib(Class, install("package"), Package);
  classgets(Result, Class);
  SET_S4_OBJECT(Result);
  
  unsigned char *result = RAW(Result); 
  int ncell = nrow*ncol;
  memset(result, 0x00, ncell);

  /* Read in data */

  int *seek = NULL;
  int snp_major = start[2];
  int nbyte = 0;
  if (snp_major) {
    if (!isNull(Rsel))
      error("can't select by rows when .bed file is in SNP-major order");
    if (!isNull(Csel)) {
      seek = INTEGER(Csel);
      nbyte = 1 + (nrow-1)/4; 
    }
  }
  else {
    if (!isNull(Csel))
      error("can't select by columns when .bed file is in individual-major order");
    if (!isNull(Rsel)) {
      seek = INTEGER(Rsel);
      nbyte = 1 + (ncol-1)/4;
    }
  }

  if (seek) 
    skip(in, seek[0]-1, nbyte);
  int part=0, ij=0, i=0, j=0;
  unsigned char byte = 0x00;
  while (1) {
    if (!part) {
      byte = (unsigned char) fgetc(in);
      if (feof(in))
	error("Unexpected end of file reached");
      part = 4;
    }
    unsigned char code = byte & mask;
    byte = byte >> 2;
    part--;
    result[ij] = recode[code];
    if (snp_major) {
      ij++;
      i++;
      if (i==nrow) {
	i = part = 0;
	j++;
	if (seek)
	  skip(in, seek[j]-seek[j-1]-1, nbyte);
	if (j==ncol)
	  break;
      }
    }	
    else {
      ij += nrow;
      j++;
      if (j==ncol){
	j = part = 0;
	i++;
	if (seek)
	  skip(in, seek[i]-seek[i-1]-1, nbyte);
	if (i==nrow)
	  break;
	ij = i;
      }
    }
  }

  UNPROTECT(4);
  return Result;
}

SEXP writebed(const SEXP Snps, const SEXP File, const SEXP SnpMajor) {
  const char *file = CHAR(STRING_ELT(File, 0));
  FILE *out = fopen(file, "wb");
  if (!out)
    error("Couldn't open output file: %s", file);
  /* Magic number */
  fputc(0x6C, out); fputc(0x1B, out);
  /* Order */
  int snpmaj = *LOGICAL(SnpMajor);
  /* Snps */
  int N = nrows(Snps);
  int M = ncols(Snps);
  unsigned char *snps = RAW(Snps);
  unsigned char recode[4] = {0x01, 0x00, 0x02, 0x03};
  unsigned char byte = 0x00;
  if (snpmaj) {
    fputc(0x01, out);
    for (int j=0, ij=0; j<M; j++) {
      for (int i=0; i<N; i++, ij++) {
	unsigned char s = snps[ij];
	if (s>3)
	  error("Uncertain genotype [%d,%d]: cannot be dealt with by this version", i, j);
	int part = i%4;
	if (!part && i) {
	  fputc(byte, out);
	  byte = 0x00;
	}
	byte = byte | (recode[s] << 2*part);
      }
      fputc(byte, out);
      byte = 0x00;
    }
  }
  else {
    fputc(0x00, out);
    for (int i=0; i<N; i++) {
      int part;
      for (int j=0, ij=i; j<M; j++, ij+=N) {
	unsigned char s = snps[ij];
	if (s>3)
	  error("Uncertain genotype [%d,%d]: cannot be dealt with by this version", i, j);
	part = j%4;
	if (!part && j) {
	  fputc(byte, out);
	  byte = 0x00;
	}
	byte = byte | (recode[s]<<2*part);
      }
      fputc(byte, out);
      byte = 0x00;
    }
  }
  fclose(out);
  return R_NilValue;
}
  
 
