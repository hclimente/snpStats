#include <R.h>
#include <Rinternals.h>
#include <stdio.h>
#include <ctype.h>

#define MAX_ID 128    /* Max space to store an id */

#define SET_TO_NA_IF_ZERO(x) if (!(x)) { (x) = NA_INTEGER ; }

#define GTYPE_N 0x00
#define GTYPE_A 0x01
#define GTYPE_H 0x02
#define GTYPE_B 0x03
#define GTYPE_INVALID       0x80
#define GTYPE_PART_MISSING  0x10
#define GTYPE_NON_DIALLELIC 0x20

int count_fields(FILE *f) {
  int nbl=0, tab=0, cont=1, nf=0;
  while (cont) {
    int ci = fgetc(f);
    switch (ci) {
    case EOF:
      return 0;
    case ' ':
      if(nbl) {
	nf++;
	nbl = 0;
      }
      break;
    case '\t':
      if (nbl || tab) 
	nf++;
      tab = 1;
      break;
    case '\n':
      cont = 0;
      if (nbl || tab) 
	nf++;
      break;
    default:
      nbl = 1;
      tab = 0;
    }
  }
  return nf;
}

/* Count non-empty lines */

int count_lines(FILE *f) {
  int nl = 0, nw = 0, cont = 1;
  while (cont) {
    int c = fgetc(f);
    switch (c) {
    case EOF:
      cont = 0;
      if (nw)
	nl++;
      break;
    case '\n':
      nl++;
      nw = 0;
      break;
    default:
      nw = 1;
    }
  }
  return nl;
}

/* Compute genotype code and, if necessary, update allele codes */

unsigned char gcode(unsigned char acodes[2], char a1, char a2, char miss, int ifX, int male) {
  /* Check for valid missing data pattern */
  if (a1==miss) {
    if (a2==miss)
      return GTYPE_N;
    else {
      if (ifX && male)
	a1 = a2;
      else
	return (GTYPE_INVALID | GTYPE_PART_MISSING);
    }
  }
  else {
    if (a2==miss) {
      if (ifX && male)
	a2 = a1;
      else
	return (GTYPE_INVALID | GTYPE_PART_MISSING);
    }
  }

  if (acodes[0]==miss) {  /* No codes yet set */
    acodes[0] = a1;
    if (a2==a1) 
      return GTYPE_A;
    else {
      acodes[1] = a2;
      return GTYPE_H;
    }
  }
  else if (acodes[1]==miss) {  /* One code only set */
    if (a1==acodes[0]) {
      if (a1==a2)
	return GTYPE_A;
      else {
	acodes[1] = a2;
	return GTYPE_H;
      }
    }
    else if (a2==acodes[0]) {
      acodes[1] = a1;
      return GTYPE_H;
    }
    else if (a1==a2) {
      acodes[1] = a1;
      return GTYPE_B;
    }
    else
      return (GTYPE_INVALID | GTYPE_NON_DIALLELIC) ;
  }
  else { /* Both codes set */
    unsigned char g = GTYPE_A;
    if (a1==acodes[1])
      g++;
    else if (a1!=acodes[0])
      return (GTYPE_INVALID | GTYPE_NON_DIALLELIC);
    if (a2==acodes[1])
      g++;
    else if (a2!=acodes[0])
      return (GTYPE_INVALID | GTYPE_NON_DIALLELIC);
    return g;
  }
}

/* The R-callable function */

SEXP readped(SEXP filename, SEXP snp_names, SEXP missing, SEXP X, SEXP sep) {
  const char *fname = CHAR(STRING_ELT(filename, 0));
  FILE *infile = fopen(fname, "r");
  if (!infile)
    error("Failed to open file %s for read", fname);
  int nf = count_fields(infile);
  if (nf % 2)
    error("Odd number of fields per line");
  int nc = (nf - 6)/2;
  int nr = count_lines(infile)+1;
  rewind(infile);
  int ifX = LOGICAL(X)[0];
  char mval = *CHAR(STRING_ELT(missing, 0));
  char sepchar = *CHAR(STRING_ELT(sep, 0));
  int i=0, j=0;

  if (length(snp_names) && nc != length(snp_names))
    error("Length of snp.names array conflicts with the pedfile");

  /* Set up R objects to hold ped data */

  SEXP Family, Member, Father, Mother, Sex, Affected;
  PROTECT(Family   = allocVector(STRSXP, nr));

  PROTECT(Member   = allocVector(INTSXP, nr));
  PROTECT(Father   = allocVector(INTSXP, nr));
  PROTECT(Mother   = allocVector(INTSXP, nr));
  PROTECT(Sex      = allocVector(INTSXP, nr));
  PROTECT(Affected = allocVector(INTSXP, nr));
  int protected = 6;

  int *member   = INTEGER(Member);
  int *father   = INTEGER(Father);
  int *mother   = INTEGER(Mother);
  int *sex      = INTEGER(Sex);
  int *affected = INTEGER(Affected);

  /* snp matrix Data part */

  SEXP Smat;
  PROTECT(Smat = allocMatrix(RAWSXP, nr, nc));
  protected++;
  unsigned char *smat = RAW(Smat);
  
  /* R objects to hold SnpMatrix attributes etc. */
  
  int *female = NULL; 
  SEXP Rnames, Acodes, diploid = R_NilValue, FSlot = R_NilValue;
  PROTECT(Rnames = allocVector(STRSXP, nr));
  PROTECT(Acodes = allocMatrix(RAWSXP, nc, 2));
  unsigned char *acodes = RAW(Acodes);
  protected += 2;
  for (i = 0 ; i < 2 * nc ; i++)
    acodes[i] = mval;
  if(ifX) {
    PROTECT(diploid = allocVector(LGLSXP, nr));
    PROTECT(FSlot  = allocVector(STRSXP, 1));
    protected += 2;
    SET_STRING_ELT(FSlot, 0, mkChar("diploid"));
    female = LOGICAL(diploid);
  }

  /* Array to indicate any errors */

  int *skip = Calloc(nc, int);
  for ( j = 0 ; j < nc ; j++)
    skip[j] = 0;
  
  /* Read data */

  int warn_col = 0;     /* whether there is at least one warning column */
  int warn_missing = 0; /* whether there is at leat one missing genotype */
  char fid[MAX_ID], fmid[MAX_ID];

  for (i = 0; i < nr ; i++) {
    /* Pedfile mandatory 6 fields */
    if (fscanf(infile, " %70s", fid) != 1)
      error("Error while reading family id on row %d", i+1 );
    SET_STRING_ELT(Family, i, mkChar(fid));
    if(fscanf(infile, " %d %d %d %d %d", 
	   member+i, father+i, mother+i, sex+i, affected+i) != 5 )
      error("Error reading pedigree information on row %d", i+1 );
    /* male sex */
    int mi = (sex[i]==1);
    if (ifX)
      female[i] = !mi;
    SET_TO_NA_IF_ZERO(sex[i]);
    SET_TO_NA_IF_ZERO(affected[i]);
    SET_TO_NA_IF_ZERO(father[i]);
    SET_TO_NA_IF_ZERO(mother[i]);

    /* Generate row name for SnpMatrix */
    int memi = member[i];
    snprintf(fmid, MAX_ID, "%s%c%d", fid, sepchar, memi);
    if (!memi)
      member[i] = NA_INTEGER;
    SET_STRING_ELT(Rnames, i, mkChar(fmid));
    /* Read genotypes */
    unsigned char *acodesj = acodes;
    unsigned char *smatj = smat + i;
    for (j = 0; j < nc ; j++) {
      char a1, a2;
      if (fscanf(infile, " %c %c", &a1, &a2) != 2)
	error("Error reading genotype %d on row %d", j+1, i+1);
      unsigned char g = gcode( acodesj , a1, a2, mval, ifX, mi); 
      if ( g & GTYPE_INVALID ) {
	warn_missing += 1;
	g = GTYPE_N;
      }
      else if ( g & GTYPE_NON_DIALLELIC ) {
	if(!skip[j]) {
	  warn_col += 1;
	  skip[j] = 1;
	}
      }
      *smatj = g;
      smatj += nr;
      acodesj += 2;
    }
  }
  
  if(warn_missing)
    warning("%i of partially missing genotype treated as missing", warn_missing);
  if (warn_col)
    warning("More than 2 alleles for %i columns", warn_col);

  /* Set any skipped columns to NA */

  unsigned char* smatij = smat;
  for (j = 0 ; j < nc ; j++) {
    if (skip[j]) {
      for(i = 0 ; i < nr ; i++, smatij++)
	*smatij = GTYPE_N;
    } else {
      smatij += nr;
    }
  }
  Free(skip);
  /* Subject support frame */

  SEXP Support, DFClass, DFNames;
  PROTECT(DFNames = allocVector(STRSXP, 6));
  PROTECT(Support = allocVector(VECSXP, 6));
  SET_VECTOR_ELT(Support, 0, Family);
  SET_VECTOR_ELT(Support, 1, Member);
  SET_VECTOR_ELT(Support, 2, Father);
  SET_VECTOR_ELT(Support, 3, Mother);
  SET_STRING_ELT(DFNames, 0, mkChar("Family"));
  SET_STRING_ELT(DFNames, 1, mkChar("Member"));
  SET_STRING_ELT(DFNames, 2, mkChar("Father"));
  SET_STRING_ELT(DFNames, 3, mkChar("Mother"));
  SET_STRING_ELT(DFNames, 4, mkChar("Sex"));
  SET_STRING_ELT(DFNames, 5, mkChar("Affected"));

  SEXP Sexlevs;
  PROTECT(Sexlevs = allocVector(STRSXP, 2));
  SET_STRING_ELT(Sexlevs, 0, mkChar("Male"));
  SET_STRING_ELT(Sexlevs, 1, mkChar("Female"));
  setAttrib(Sex, R_LevelsSymbol, Sexlevs);
  SEXP Fclass;
  PROTECT(Fclass = allocVector(STRSXP, 1));
  SET_STRING_ELT(Fclass, 0, mkChar("factor"));
  classgets(Sex, Fclass);
  SET_VECTOR_ELT(Support, 4, Sex);
  SEXP Afflevs;
  PROTECT(Afflevs = allocVector(STRSXP, 2));
  SET_STRING_ELT(Afflevs, 0, mkChar("Unaffected"));
  SET_STRING_ELT(Afflevs, 1, mkChar("Affected"));
  setAttrib(Affected, R_LevelsSymbol, Afflevs);
  classgets(Affected, duplicate(Fclass));
  SET_VECTOR_ELT(Support, 5, Affected);
  setAttrib(Support, R_RowNamesSymbol, Rnames);
  setAttrib(Support, R_NamesSymbol, DFNames);
  PROTECT(DFClass = allocVector(STRSXP, 1));
  SET_STRING_ELT(DFClass, 0, mkChar("data.frame"));
  classgets(Support, DFClass);
  protected += 6;

  /* SNP matrix */

  SEXP SMClass, DimNames, Cnames;
  PROTECT(SMClass = allocVector(STRSXP, 1));
  if (ifX) {
    SET_STRING_ELT(SMClass, 0, mkChar("XSnpMatrix"));
    R_do_slot_assign(Smat, FSlot, diploid); 
  }
  else
    SET_STRING_ELT(SMClass, 0, mkChar("SnpMatrix"));
  classgets(Smat, SMClass);
  SET_S4_OBJECT(Smat);
  PROTECT(DimNames = allocVector(VECSXP, 2));
  SET_VECTOR_ELT(DimNames, 0, duplicate(Rnames));
  protected += 2;
  if (length(snp_names)) 
    SET_VECTOR_ELT(DimNames, 1, duplicate(snp_names));
  else {
    PROTECT(Cnames = allocVector(STRSXP, nc));
    protected++;
    for (j=0; j<nc; j++) {
      sprintf(fmid,"%d", j+1);
      SET_STRING_ELT(Cnames, j, mkChar(fmid));
    }
    SET_VECTOR_ELT(DimNames, 1, Cnames);
  }
  setAttrib(Smat, R_DimNamesSymbol, DimNames);
  
  /* Result list */

  SEXP Result, ResNames;
  PROTECT(Result = allocVector(VECSXP, 2));
  SET_VECTOR_ELT(Result, 0, Smat);
  SET_VECTOR_ELT(Result, 1, Support);
  PROTECT(ResNames = allocVector(STRSXP, 2));
  SET_STRING_ELT(ResNames, 0, mkChar("snps"));
  SET_STRING_ELT(ResNames, 1, mkChar("subject.support"));
  setAttrib(Result, R_NamesSymbol, ResNames);
  protected += 2;

  UNPROTECT(protected);
  return (Result);
}

