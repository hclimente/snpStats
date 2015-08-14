#include <R.h>
#include <Rinternals.h>
#include <stdio.h>
#include <string.h>
#include "hash_index.h"
#include "Rmissing.h"

/* rbind and cbind for SnpMatrix and XSnpMatrix classes */

SEXP snp_rbind(SEXP args) {
  int nb = length(args) - 1;
  SEXP argsin = args;
  const char *class = NULL;
  SEXP Cnames = R_NilValue;
  SEXP Class = R_NilValue; 
  int nr = 0, nc=0;
  for (int i=0; i<nb; i++) {
    args = CDR(args);
    const SEXP This = CAR(args);
    Class = getAttrib(This, R_ClassSymbol);
    if (TYPEOF(Class) != STRSXP) {
      Class = R_data_class(This, FALSE);
    }
    const char *cli = CHAR(STRING_ELT(Class, 0));
    if(!IS_S4_OBJECT(This)) {
      warning("rbinding SnpMatrix object without S4 object bit");
    }
    int nci = ncols(This);
    nr += nrows(This);
    SEXP dn = getAttrib(This, R_DimNamesSymbol);
    if (dn==R_NilValue)
      error("Missing dimnames attribute in SnpMatrix object");
    SEXP cni = VECTOR_ELT(dn, 1);
    if (cni==R_NilValue)
      error("Missing column names in SnpMatrix object");
    SEXP rni = VECTOR_ELT(dn, 0);
    if (rni==R_NilValue)
      error("Missing row names in SnpMatrix object");
    if (!i) { /* First arg */
      class = cli;
      if ((strcmp(class,"SnpMatrix")!=0) && (strcmp(class,"XSnpMatrix")!=0))
	error("argument not a SnpMatrix");
      nc = nci;
      Cnames = cni;
    }
    else {
      if (strcmp(class, cli) != 0)
	error("arguments have incompatible classes");
      if (nci != nc)
	error("matrices have unequal number of columns");
      if (cni != R_NilValue) {
	if (Cnames == R_NilValue) 
	  Cnames = cni;
	else for (int j=0; j<nc; j++) {
	  const char *one = CHAR(STRING_ELT(Cnames, j));
	  const char *other = CHAR(STRING_ELT(cni, j));
	  if (strcmp(one, other) != 0)
	    error("column names do not match");
	}
      }
    }
  }

  /* Result matrix */

  SEXP Result, Rnames, Dnames, Diploid = R_NilValue;
  PROTECT(Result = allocMatrix(RAWSXP, nr, nc));
  classgets(Result, duplicate(Class));
  SET_S4_OBJECT(Result);
  PROTECT(Rnames = allocVector(STRSXP, nr));
  PROTECT(Dnames = allocVector(VECSXP, 2));
  SET_VECTOR_ELT(Dnames, 0, Rnames);
  SET_VECTOR_ELT(Dnames, 1, duplicate(Cnames));
  setAttrib(Result, R_DimNamesSymbol, Dnames);
  int *diploid = NULL;
  int X = (strcmp(class, "XSnpMatrix") == 0);
  if (X) {
    PROTECT(Diploid = allocVector(LGLSXP, nr));
    R_do_slot_assign(Result, mkString("diploid"), Diploid);
    diploid = LOGICAL(Diploid);
  }
  unsigned char *result = RAW(Result);

  /* Fill Result matrix */

  args = argsin;
  int nri;
  index_db row_index = index_create(nr);
  int rows_done=0;
  for (int i=0; i<nb; i++, rows_done += nri) {
    unsigned char *r = result + rows_done;
    args = CDR(args);
    const SEXP This = CAR(args);
    nri = nrows(This);
    /* Copy matrix body */
    unsigned char *this = RAW(This);
    unsigned char *rj = r;
    for (int j=0; j<nc; j++, rj+=nr) {
      unsigned char *rjk = rj;
      for (int k=0; k<nri; k++) 
	*(rjk++) = *(this++);
    }
    /* Copy row names */
    SEXP dn = getAttrib(This, R_DimNamesSymbol);
    if( dn != R_NilValue) {
      SEXP rni = VECTOR_ELT(dn, 0);
      if( rni != R_NilValue) {
	for (int k=0, nk=rows_done; k<nri; k++, nk++) {
	  SEXP rnik = STRING_ELT(rni, k);
	  if (rnik != R_NilValue) {
	    SET_STRING_ELT(Rnames, nk, rnik);
	    if (index_insert(row_index, CHAR(rnik), k)!=0)
	      warning("Duplicated row name at row %d overall from row %d of object %d", nk+1, k+1, i+1);
	  }
	}
      }
    }
    /* Copy diploid indicators */
    if (X) {
      SEXP Di = R_do_slot(This, mkString("diploid"));
      int *di = LOGICAL(Di);
      for (int k=0, nk=rows_done; k<nri; k++, nk++)
	diploid[nk] = di[k];
    }
  }
  if (X) {
    /* copying row names to diploid slot names;
       the input row names should agree with the input slot names,
       so the worst case is that output slots have names which 
       input slots don't have */
    setAttrib(Diploid, R_NamesSymbol, duplicate(Rnames));
  }
  index_destroy(row_index);
  UNPROTECT(X? 4:3);
  return(Result);
}

SEXP snp_cbind(SEXP args) {
  int X = FALSE, nb = length(args) - 1;
  SEXP argsin = args;
  const char *class = NULL;
  SEXP Diploid = R_NilValue;
  int *diploid = NULL; 
  SEXP Rnames = R_NilValue, Class = R_NilValue;
  int nr = 0, nc=0;
  for (int i=0; i<nb; i++) {
    args = CDR(args);
    SEXP This = CAR(args);
    Class = getAttrib(This, R_ClassSymbol);
    if (TYPEOF(Class) != STRSXP) {
      Class = R_data_class(This, FALSE);
    }
    const char *cli = CHAR(STRING_ELT(Class, 0));
    if(!IS_S4_OBJECT(This)) {
      warning("cbinding SnpMatrix object without S4 object bit");
    }
    SEXP Di = R_NilValue;
    int *di = NULL;
    X = (strcmp(cli, "XSnpMatrix")==0);
    if (X) {
      Di = R_do_slot(This, mkString("diploid"));
      di = LOGICAL(Di);
    }
    int nri = nrows(This);
    nc += ncols(This);
    /* Copy column names, check row names */
    SEXP dn = getAttrib(This, R_DimNamesSymbol);
    if (dn==R_NilValue)
      error("Missing dimnames attribute in SnpMatrix object");
    SEXP cni = VECTOR_ELT(dn, 1);
    if (cni==R_NilValue)
      error("Missing column names in SnpMatrix object");
    SEXP rni = VECTOR_ELT(dn, 0);
    if (rni==R_NilValue)
      error("Missing row names in SnpMatrix object");
    if (!i) { /* First arg */
      class = cli;
      if ((strcmp(class,"SnpMatrix")!=0) && (strcmp(class,"XSnpMatrix")!=0))
	error("argument not a SnpMatrix");
      Rnames = rni;
      nr = nri;
      if (X) {
	Diploid = Di;
	diploid = di;
      }
    }
    else {
      if (strcmp(class, cli) != 0)
	error("incompatible argument classes");
      if (nri != nr)
	error("unequal number of rows");
      for (int j=0; j<nr; j++) {
	const char *one = CHAR(STRING_ELT(Rnames, j));
	const char *other = CHAR(STRING_ELT(rni, j));
	if (strcmp(one, other) != 0)
	  error("row names do not match");
	if (X && (diploid[j]!=di[j]))
	  error("inconsistent ploidy in row %d", j+1);
      }
    }
  }

  /* Result matrix */

  SEXP Result, Cnames, Dnames;
  PROTECT(Result = allocMatrix(RAWSXP, nr, nc));
  classgets(Result, duplicate(Class));
  SET_S4_OBJECT(Result);
  PROTECT(Dnames = allocVector(VECSXP, 2));
  setAttrib(Result, R_DimNamesSymbol, Dnames);
  PROTECT(Cnames = allocVector(STRSXP, nc));
  SET_VECTOR_ELT(Dnames, 0, duplicate(Rnames));
  SET_VECTOR_ELT(Dnames, 1, Cnames);
  if (X) 
    R_do_slot_assign(Result, mkString("diploid"), duplicate(Diploid));
  unsigned char *result;
  result = RAW(Result);

  /* Fill Result matrix */

  args = argsin;
  index_db col_index = index_create(nc);
  R_xlen_t ij = 0;
  for (int i=0; i<nb; i++) {
    args = CDR(args);
    SEXP This = CAR(args);
    /* Copy matrix body */
    unsigned char *this = RAW(This);
    int nci = ncols(This);
    R_xlen_t len = length(This);
    for (R_xlen_t j=0; j<len; j++) 
      *(result++) = *(this++);
    /* Copy column names */
    SEXP dn = getAttrib(This, R_DimNamesSymbol);
    if(dn != R_NilValue) {
      SEXP cni = VECTOR_ELT(dn, 1);
      if(cni != R_NilValue) {
	for (int j=0; j<nci; j++, ij++) {
	  SEXP cnij = STRING_ELT(cni, j);
	  if (cnij != R_NilValue) {
	    SET_STRING_ELT(Cnames, ij, cnij);
	    if (index_insert(col_index, CHAR(cnij), ij)!=0)
	      error("Duplicated column name at column %d overall from column %d of object %d", ij+1, j+1, i+1);
	  }
	}
      }
    } else {
      Rprintf("names empty\n");
    }
  }
  index_destroy(col_index);
  /* in cbind we never create a new diploid slot, but only copy 
     from the first one so unlike rbind() we don't need to decide 
     unprotect level */ 
  UNPROTECT(3); 
  return(Result);
}


