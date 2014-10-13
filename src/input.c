/* Read genotype data in long format */

#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <stdio.h>
#include <R.h>
#include <Rmath.h>
#include <Rdefines.h>
#include <Rinternals.h>
#include "hash_index.h"
#include "zlib.h"

#define MAX_FLD 128

/* Declare routines called but defined later */

int str_match(const char *a, const char *b, const int forward);
int str_inlist(const SEXP strlist, const char *target);
int nucleotide(char c);
SEXP simplify_names(const SEXP x);
int next_field(const gzFile infile, const char sep, const char comment, 
	       const char replace, char *field, const int len_field);
int skip_to_eol(const gzFile infile);
int recode_snp(unsigned char *matrix, const int N, const int M);

/* Function defined in imputation.c */

index_db create_name_index(const SEXP names);


/* Main C function */

SEXP insnp_new(const SEXP Filenames, const SEXP Sample_id, const SEXP Snp_id,
	       const SEXP Diploid, const SEXP Fields,  
	       const SEXP Codes, const SEXP Threshold, const SEXP Lower, 
	       const SEXP Sep, const SEXP Comment, const SEXP Skip, 
	       const SEXP Simplify, const SEXP Verbose, 
	       const SEXP In_order, const SEXP Every){

  /* Process arguments */ 
  
  if (TYPEOF(Verbose)!=LGLSXP)
    error("Argument type error: Verbose");
  if (length(Verbose)>1)
    warning("Only first element of argument used: Verbose");
  int verbose = *INTEGER(Verbose);

  if (TYPEOF(In_order)!=LGLSXP)
    error("Argument type error: In_order");
  if (length(In_order)>1)
    warning("Only first element of argument used: In_order");
  int in_order = *INTEGER(In_order);

  if (TYPEOF(Filenames)!=STRSXP)
    error("Argument type error: Filenames");
  int Nfile = length(Filenames);

  int Nsample=0;
  if (TYPEOF(Sample_id)==STRSXP)
    Nsample = length(Sample_id);
  else if (TYPEOF(Sample_id)!=NILSXP)
    error("Argument type error: Sample_id");
  
  int Nsnp=0;
  if (TYPEOF(Snp_id)==STRSXP)
    Nsnp = length(Snp_id);
  else if (TYPEOF(Snp_id)!=NILSXP)
    error("Argument type error: Snp_id");
  
  /* file is 1 = sample, 2 = snp, 0 = irrelevant */

  int file_is = 0;
  if (!Nsample) {
    if (Nsnp) {
      Nsample = Nfile;
      Rprintf("Each file is assumed to concern a single sample\n"); 
      Rprintf("(Sample IDs are assumed to be included in filenames)\n");
      file_is = 1;
    }
    else 
      error("No sample or SNP IDs specified");
  }
  else if (!Nsnp) {
    Nsnp = Nfile;
    Rprintf("Each file is assumed to concern a single SNP\n");
    Rprintf("(SNP IDs are assumed to be included in filenames)\n");
    file_is = 2;
  }

  /* If not in order, set up hash tables */

  index_db sample_index = NULL;
  index_db snp_index = NULL;
  if (!in_order) {
    if (file_is != 1) 
      sample_index = create_name_index(Sample_id);
    if (file_is != 2) 
      snp_index = create_name_index(Snp_id);
  }

  int *diploid=NULL;
  if (TYPEOF(Diploid)==LGLSXP) {
    if (length(Diploid)!=Nsample)
      error("Argument length error: diploid argument");
    diploid = LOGICAL(Diploid);
  }
  else if (TYPEOF(Diploid)!=NILSXP)
    error("Argument type error: diploid argument");
  
  if (TYPEOF(Fields)!=INTSXP) 
    error("Argument type error: Fields");
  int *fields = INTEGER(Fields);
  int fsamp=0, fsnp=0, fgt=0, fa1=0, fa2=0, fconf=0;
  SEXP Fnames = getAttrib(Fields, R_NamesSymbol);
  if (TYPEOF(Fnames)==NILSXP) 
    error("Argument error: Fields argument has no names");
  int fmax = 0;
  int Nfield = length(Fields);
  for (int i=0; i<Nfield; i++) {
    const char *fname = CHAR(STRING_ELT(Fnames, i));
    int fi = fields[i];
    if (!strcmp(fname, "sample"))
      fsamp = fi;
    else if (!strcmp(fname, "snp"))
      fsnp = fi;
    else if (!strcmp(fname, "genotype"))
      fgt = fi;
    else if (!strcmp(fname, "allele1"))
      fa1 = fi;
    else if (!strcmp(fname, "allele2"))
      fa2 = fi;
    else if (!strcmp(fname, "confidence"))
      fconf = fi;
    else
      error("Unrecognized input field name: %s", fname);
    if (fi>fmax) 
      fmax = fi;
  }
  if (verbose) {
    Rprintf("Reading one call per input line\n");
    if (fsamp) {
      Rprintf("   Sample id is in field %d", fsamp);
      if (file_is==1) 
	Rprintf(" (ignored)\n");
      else
	Rprintf("\n");
    }
    if (fsnp) {
      Rprintf("   SNP id is in field %d", fsnp);
      if (file_is==2) 
	Rprintf(" (ignored)\n");
      else
	Rprintf("\n");
    }
    if (fgt)
      Rprintf("   Genotype is in field %d\n", fgt);
    if (fa1)
      Rprintf("   Allele 1 is in field %d\n", fa1);
    if (fa2)
      Rprintf("   Allele 2 is in field %d\n", fa2);
    if (fconf)
      Rprintf("   Confidence score is in field %d\n", fconf);
  }
  
  if (file_is==1)
    fsamp = 0;
  else if (file_is==2)
    fsnp = 0;


  /* Allele or genotype coding? */

  int gcoding;
  if (fgt) {
    if (fa1 || fa2) 
      error("Coding must be by genotype OR allele, not both");
    gcoding = 1;
  }
  else {
    if (!(fa1 && fa2)) 
      error("No genotype or allele field(s) specified");
    if (!(fa1 && fa2)) 
      error("Field positions for both alleles must be specified"); 
    gcoding = 0;
  }

  int nuc = 0;
  if (TYPEOF(Codes)!=STRSXP)
    error("Argument type error: Codes");
  if (length(Codes)==1) {
    SEXP Code = STRING_ELT(Codes, 0);
    const char *code = CHAR(Code);
    if (!strcmp(code, "nucleotide"))
      nuc = 1;
    else
      error("Unrecognized coding: %s", code);
  }
  else {
    int ncode = length(Codes);
    if (gcoding) {
      if (diploid) {
	if (ncode!=5) {
	  if (ncode==3)
	    warning("Genotype coding for X: haploid genotypes are assumed to be coded as homozygous");
	  else
	    error("Genotype coding for X.snp: three or five genotype codes must be specified");
	}
      }
      else {
	if (ncode!=3)
	  error("Genotype coding: three genotype codes must be specified");
      }
    }
    else {
      if (ncode!=2) 
	error("Allele coding: two allele codes must be specified");
    }
  }

  if (TYPEOF(Threshold)==NILSXP && !fconf) 
    error("Argument type error: no threshold argument");
  if (TYPEOF(Threshold)!=REALSXP)
    error("Argument type error: Threshold");
  double threshold = *REAL(Threshold);
  if (fconf && threshold==NA_REAL)
    error("Confidence score is read but no threshold is set");

  if (TYPEOF(Lower)!=LGLSXP)
    error("Argument type error: Lower");
  if (length(Lower)>1)
    warning("Only first element of argument used: Lower");
  int lower = *INTEGER(Lower);

  char sep = ' ';
  if (TYPEOF(Sep)==STRSXP) {
    if (length(Sep)>1)
      warning("Only first element of argument used: Sep");
    const char *c = CHAR(STRING_ELT(Sep, 0));
    if (strlen(c)>1) 
      warning("Only first character used: Sep");
    sep = c[0];
  }
  else if (TYPEOF(Sep)!=NILSXP) 
    error("Argument type error: Sep");

  char comment = (char) 0;
  if (TYPEOF(Comment)==STRSXP) {
    if (length(Sep)>1)
      warning("Only first element of argument used: Comment");
    const char *c = CHAR(STRING_ELT(Comment, 0));
    if (strlen(c)>1) 
      warning("Only first character used: Comment");
    comment = c[0];
  }
  else if (TYPEOF(Comment)!=NILSXP) 
    error("Argument type error: Comment");

  int skip = 0;
  if (TYPEOF(Skip)==INTSXP) {
    if (length(Skip)>1)
      warning("Only first element used: Skip");
    skip = INTEGER(Skip)[0];
  }
  else if (TYPEOF(Skip)!=NILSXP) 
    error("Argument type error: Skip");
 
  if (TYPEOF(Simplify)!=LGLSXP)
    error("Argument type error: Simplify");
  if (length(Simplify)>2)
    error("Argument length error: Simplify");
  int *simplify = INTEGER(Simplify);
  if (length(Simplify)==1)
    simplify[1] = simplify[0];

  int every=0;
  if (TYPEOF(Every)==INTSXP) {
    if (length(Every)>1) 
      warning("Only first element used: Every");
    every = INTEGER(Every)[0];
  }
  else if (TYPEOF(Every)!=NILSXP)
    error("Argument type error: Every");
    

  /* Create output object and initialise to zero */

  if (verbose) {
    if (diploid)
      Rprintf("Reading XSnpMatrix with %d rows and %d columns\n", 
	      Nsample, Nsnp);
    else
      Rprintf("Reading SnpMatrix with %d rows and %d columns\n", 
	      Nsample, Nsnp);
  }
  SEXP Result, Dimnames, Package, Class;
  PROTECT(Result = allocMatrix(RAWSXP, Nsample, Nsnp));
  PROTECT(Dimnames = allocVector(VECSXP, 2));
  if (simplify[0]) {
    SET_VECTOR_ELT(Dimnames, 0, 
		   simplify_names(file_is==1? Filenames: Sample_id));
  }
  else {
    SET_VECTOR_ELT(Dimnames, 0, 
		   duplicate(file_is==1? Filenames: Sample_id));
  }
  if (simplify[1]) {
    SET_VECTOR_ELT(Dimnames, 1, 
		   simplify_names(file_is==2? Filenames: Snp_id));
  }
  else {
    SET_VECTOR_ELT(Dimnames, 1, 
		   duplicate(file_is==2? Filenames: Snp_id));
  }
  setAttrib(Result, R_DimNamesSymbol, Dimnames);

  /* Class */

  PROTECT(Class = allocVector(STRSXP, 1));
  if (diploid) {
    R_do_slot_assign(Result, mkString("diploid"), Diploid);
    SET_STRING_ELT(Class, 0, mkChar("XSnpMatrix"));
  }
  else {
    SET_STRING_ELT(Class, 0, mkChar("SnpMatrix"));
  }
  PROTECT(Package = allocVector(STRSXP, 1));
  SET_STRING_ELT(Package, 0, mkChar("snpStats"));
  setAttrib(Class, install("package"), Package);
  classgets(Result, Class);
  SET_S4_OBJECT(Result);
  unsigned char *result = RAW(Result);
  memset(result, 0x00, Nsample*Nsnp);

  /* Read in data */

  char field[MAX_FLD];
  int Naccept = 0, Nreject = 0, Nocall = 0, Nskipped = 0, Nxerror = 0;
  int i_this = 0, j_this = 0;
  const char *this_sample=NULL, *this_snp=NULL;
  if (fsamp) {
    this_sample = CHAR(STRING_ELT(Sample_id, 0)); 
  }
  if (fsnp) {
    this_snp = CHAR(STRING_ELT(Snp_id, 0));
  }
  if (verbose) {
    Rprintf("                             Cumulative totals\n");
    Rprintf("                    -----------------------------------\n");
    Rprintf("    File     Line   Accepted Rejected  No call  Skipped    File name\n");
  }
  
  /* slowest varying 0 = don't know, 1 = sample, 2 = snp */

  int slowest = file_is, last=Nsample*Nsnp-1;
  int advance = 0, finished=0;

  for (int f=0; f<Nfile; f++) {
    /* Open input file */
    const char *filename = CHAR(STRING_ELT(Filenames, f));
    if (verbose) {
      int lfn = strlen(filename); 
      if (lfn > 20) {
	Rprintf("%59s...%-17s\r", "", filename+lfn-17);
      }
      else
	Rprintf("%59s%-20s\r", "", filename);
    }
    gzFile infile = gzopen(filename, "rb");
    if (!infile) {
      warning("Failure to open input file: %s", filename);
      continue;
    }
    int fterm = 2, line = 0, found_in_file = 0;
    /* Skip any header lines */
    for (int i=0; i<skip; i++) {
      line++;
      if (skip_to_eol(infile)==3)
	error("End-of-file reached on line %d", line);
    }
    Nskipped += skip;
    /* Read data lines */
 
    while (fterm!=3) {

      /* Read a line */

      line++;
      if (verbose && every && !(line % every)) 
	Rprintf("%8d %8d %10d %8d %8d %8d\r", 
		f+1, line, Naccept, Nreject, Nocall, Nskipped);
      int genotype=0, allele1=0, allele2=0;
      /* wanted is coded as:
	 1  if this call is to be accepted
	 0  if it is not wanted (or a comment line)
	 -1 if rejected due to insufficient confidence
	 -2 coded as no-call
      */
      int wanted = 1; 
      char sampid[MAX_FLD], snpid[MAX_FLD];
      char gtype1[MAX_FLD], gtype2[MAX_FLD];
      char cscore[MAX_FLD];
      sampid[0] = snpid[0] = cscore[0] = (char) 0;
      for (int r=1; (r<=fmax); r++) {
	fterm = next_field(infile, sep, comment, '_', field, MAX_FLD);
	if (!fterm) 
	  error("Field overflow: line %d, field %d", line, r);
	if ((fterm>1) && (r<fmax)) {
	  if (r==1) {
	    if(!field[0]) {
	    /* Empty line or comment line */
	      wanted = 0;
	      if (fterm==2) {
		Nskipped++;
		continue;
	      }
	      else
		break;
	    }
	  }
	  error("Incomplete line: %d (last field read: %d = %s)", 
		line, r, field);
	}
	/* Save fields */
	if (r==fsamp) {
	  strncpy(sampid, field, MAX_FLD-1);
	}
	else if (r==fsnp) {
	  strncpy(snpid, field, MAX_FLD-1);
	}
	else if (r==fconf) {
	  strncpy(cscore, field, MAX_FLD-1);
	}
	else if (r==fgt) 
	  strncpy(gtype1, field, MAX_FLD-1);
        else if (r==fa1) 
	  strncpy(gtype1, field, MAX_FLD-1);
	else if (r==fa2) 
	  strncpy(gtype2, field, MAX_FLD-1);
	else {
	  /* skip field */
	}
      } /* Matches: for (int r=1; (r<=fmax); r++) { */

      if (!wanted) 
	continue;

      /* Discard any further fields */

      if (fterm<2) {
	fterm = skip_to_eol(infile);
      }

      /* Find next target and check matches */

      int match_sample, match_snp;
      if (in_order) {

	/* Advance to next target read */

	if (advance) {
	  
	  /* 
	     If unknown, determine sort order by seeing which indicator 
	     has changed
	  */

	  if (!slowest) {
	    if (strcmp(this_sample, sampid)) 
	      slowest = 2; /* sample fastest, SNP slowest */
	    else if (strcmp(this_snp, snpid))
	      slowest = 1; /* SNP fastest, sample slowest */
	    else
	      error("Error in input file sort order");
	  }
	  
	  /* Now advance fastest varying indicator */
	  
	  if (slowest==1) {
	    j_this++;
	    if (j_this==Nsnp) {
	      j_this = 0;
	      if (fsamp) { 
		i_this++;
		if (i_this==Nsample) {
		  finished = 1;
		  break;
		}
		else
		  this_sample =  CHAR(STRING_ELT(Sample_id, i_this)); 
	      }
	    }
	    this_snp =  CHAR(STRING_ELT(Snp_id, j_this));
	  }
	  else {
	    i_this++;
	    if (i_this==Nsample) {
	      i_this = 0;
	      if (fsnp) {
		j_this++;
		if (j_this==Nsnp) {
		  finished = 1;
		  break;
		}
		else 
		  this_snp =  CHAR(STRING_ELT(Snp_id, j_this));
	      }
	    }
	    this_sample =  CHAR(STRING_ELT(Sample_id, i_this));
	  }
	}
	
	/* Does current line match current target? */
	
	match_sample = file_is==1 || !strcmp(this_sample, sampid);
	match_snp = file_is==2 || !strcmp(this_snp, snpid);
      }

      else { /* Not in order */

	if (file_is!=1) {
	  if (this_sample && strcmp(this_sample, sampid))
	    i_this = index_lookup(sample_index, sampid);
	  if (i_this >= 0) {
	    this_sample =  CHAR(STRING_ELT(Sample_id, i_this)); 
	    match_sample = 1;
	  }
	  else {
	    this_sample = NULL;
	    match_sample = 0;
	  }
	}
	else {
	  match_sample = 1;
	  i_this = f;
	}
	if (file_is!=2) {
	  if (this_snp && strcmp(this_snp, snpid))
	    j_this = index_lookup(snp_index, snpid);
	  if (j_this >= 0) {
	    this_snp = CHAR(STRING_ELT(Snp_id, j_this));
	    match_snp = 1;
	  }
	  else {
	    this_snp = NULL;
	    match_snp = 0;
	  }
	}
	else {
	  match_snp = 1;
	  j_this = f;
	}
      }

      if (match_sample && match_snp) {

	/* Next target read found in file(s) */

	found_in_file++;
	int ij_this = j_this*Nsample + i_this;
	finished = (ij_this==last);
	
	/* Check confidence score */
	
	if (fconf) {
	  double conf;
	  if (sscanf(cscore, "%lf", &conf)!=1) 
	    error("Failure to read confidence score: line %d", line);
	  if ((lower && conf<threshold) || (!lower && conf>threshold)) {
	    wanted = -1;
	    Nreject++;
	    if (finished)
	      break;
	    else
	      advance = 1;
	    continue;
	  }
	}
	
	/* Decode genotype */
	
	int which;
	if (gcoding) {
	  if (nuc) {
	    switch (strlen(gtype1)) {
	    case 0:
	      allele1 = allele2 = 0;
	      break;
	    case 1:
	      allele1 = allele2 = nucleotide(gtype1[0]);
	      break;
	    case 2:
	      allele1 = nucleotide(gtype1[0]);
	      allele2 = nucleotide(gtype1[1]);
	      break;
	    default:
	      error("Nucleotide coded genotype should be 2 character string: line %d", line);
	    }
	  }
	  else {
	    which = str_inlist(Codes, gtype1);
	    if (!which)
	      genotype = 0;
	    else if (which>3) 
	      genotype = 2*which - 7;
	    else
	      genotype = which;
	  }
	}
	else {
	  if (nuc) {
	    allele1 = nucleotide(gtype1[0]);
	    allele2 = nucleotide(gtype2[0]);
	  }
	  else {
	    allele1 = str_inlist(Codes, gtype1);
	    allele2 = str_inlist(Codes, gtype2);
	  }
	}

	/* Successful read, store genotype in result[i_this, j_this] */
	
	if (nuc || !gcoding) {
	  if (allele2 < allele1) {
	    genotype = allele2;
	    allele2 = allele1;
	    allele1 = genotype;
	  }
	  if (allele1 && allele2) 
	    genotype = allele1 + (allele2*(allele2-1))/2;	
	  else
	    genotype = 0;
	}
	if (genotype) {
	  if (diploid && !diploid[i_this] && (genotype==2))
	    Nxerror++;
	  else {
	    Naccept++;
	    result[ij_this] = (unsigned char) genotype;
	  }
	}
	else {
	  wanted = -2;
	  Nocall++;
	}
	
	/* Flag need to advance to next target */

	advance = 1;
      } /* matches if (match_sample && match_snp) { */
      else {
	Nskipped++;

	/* Flag no advance to next target */

	advance = 0;	
      }
      if (finished)
	break;
    } /* matches: while (fterm!=3) { */
    if (file_is==1)  
      i_this++;
    if (file_is==2) 
      j_this++;
    if (verbose) {
      Rprintf("%8d %8d %10d %8d %8d %8d\r", f+1, line, Naccept, Nreject, Nocall, Nskipped);
    }
    if(!found_in_file)
      warning("No calls found in file %s", filename);
    gzclose(infile);
    if (finished)
      break;
  }

  /* Warnings */

  if (in_order && !finished) 
    warning("End of data reached before search completed");
  if (Nxerror) 
    warning("%d haploid genotypes were coded as heterozygous; set to NA", Nxerror);
  if (Nskipped)
    warning("%d lines of input file(s) were skipped", Nskipped);

  /* Report */

  if (verbose)
    Rprintf("\n");
  Rprintf("%d genotypes successfully read\n", Naccept);
  if (Nreject)
    Rprintf("%d genotypes were rejected due to low confidence\n", Nreject);
  if (Nocall)
    Rprintf("%d genotypes were not called\n", Nocall);
  if (Nsample*Nsnp > Naccept+Nxerror+Nreject+Nocall)
    Rprintf("%d genotypes could not be found on input file(s)\n",
	    Nsample*Nsnp - Naccept - Nreject - Nxerror-Nocall);
  if (nuc) {
    if (verbose)
      Rprintf("Recasting and checking nucleotide coding\n");
    int none_snps = recode_snp(result, Nsample, Nsnp);
    if (none_snps) {
      Rprintf("%d polymorphisms were not SNPs and have been set to NA ", 
	      none_snps);
      Rprintf("(see warnings for details)\n");
    }
  }
  UNPROTECT(4);
  
  /* Destroy hash indexes */

  if (sample_index)
    index_destroy(sample_index);
  if (snp_index)
    index_destroy(snp_index);

  return Result;
}

/* 
   String match. Returns number of characters matching in forward or reverse
   directions
*/

int str_match(const char *a, const char *b, const int forward) {
  int len = 0;
  if (forward) {
    while (*a && *b && (*a==*b)) {
      len++;
      a++;
      b++;
    }
  }
  else {
    int la = strlen(a), lb = strlen(b);
    while(la-- && lb--) {
      if (a[la] == b[lb])
	len++;
      else
	break;
    }
  }
  return len;
}

/* 
   Find string in a list of strings. In a list of length L returns 
   1:L if found and 0 if not
*/

int str_inlist(const SEXP strlist, const char *target){
  int L = length(strlist);
  for (int i=0; i<L; i++) {
    if (!strcmp(target, CHAR(STRING_ELT(strlist, i)))) 
      return i+1;
  }
  return 0;
}

/* Nucleotide A, C, G, T -> 1, 2, 3, 4. 0 for no match */

int nucleotide(char c) {
  char cu = toupper(c);
  switch (cu) {
  case 'A': return 1;
  case 'C': return 2;
  case 'G': return 3;
  case 'T': return 4;
  default: return 0;
  }
}
  

/* 
   Copy and simplify identifying character vector by removal of any 
   repeated leading or trailing character sequence
*/

SEXP simplify_names(const SEXP x) {
  if (TYPEOF(x)!=STRSXP)
    error("simplify: argument type error");
  int len = length(x);
  int lenf=0, lenb=0;
  if (len>1) {
    char front[MAX_FLD], back[MAX_FLD];
    strncpy(front, CHAR(STRING_ELT(x, 0)), MAX_FLD-1);
    strncpy(back, front, MAX_FLD-1);
    lenf = lenb = strlen(front);
    char *end_back = back + lenb; /* pointer to terminating 0 */
    char *start_back = back;
    for (int i=1; i<len; i++) {
      const char *xi = CHAR(STRING_ELT(x, i));
      if (lenf) {
	lenf = str_match(front, xi, 1);
	front[lenf] = (char) 0;
      }
      if (lenb) {
	lenb = str_match(start_back, xi, 0);
	start_back = end_back - lenb ;
      }
    }
  }
  
  SEXP Result;
  PROTECT(Result = allocVector(STRSXP, len));
  char id[MAX_FLD];
  for (int i=0; i<len; i++) {
    const char *xi = CHAR(STRING_ELT(x, i));
    int lenx = strlen(xi);
    int ncp = lenx - lenf - lenb;
    if (ncp>=MAX_FLD) 
      error("simplify: id length overflow");
    strncpy(id, xi+lenf, ncp);
    *(id+ncp) = (char) 0;
    SET_STRING_ELT(Result, i, mkChar(id));
  }
  UNPROTECT(1);
  return Result;
}
    

/*
  Read to next separator character OR end-of-line OR end-of-file
  Returns 1, 2 or 3 according to termination. 
  Leading blanks are skipped
  After this, characters are copied to <field>
  Any characters after <comment> are skipped and end-of-line (1) returned
  On return <field> has trailing blanks stripped
  Blanks within the returned string are replaced by <replace>
  Returns 0 if <field> overflows
*/

  
int next_field(const gzFile infile, const char sep, const char comment, 
	       const char replace, char *field, const int len_field) {
  const int blank = (int)' ';
  const int delim = (int) sep;
  const int eol = (int) '\n';
  const int eof = -1;
  const int com = (int) comment;

  /* skip to first none-blank character */

  field[0] = (char) 0;
  int c;
  /* Skip leading blanks */
  while ((c=gzgetc(infile))==blank) {};
  int non_blank = 1;
  char *term = field + len_field;
  int len = 0, result = 0;
  do {
    if (c==eol) {
      result = 2;
    }
    else if (c==eof) {
      result = 3;
    }
    else if (com && c==com) {
      /* Skip to end of line */
      while (1) {
	c=gzgetc(infile);
        if (c==eol) {
	  result = 2;
	  break;
	}
	if (c==eof) {
	  result = 3;
	  break;
	}
      }
    }
    else if (c==delim) {
      result = 1;
    }
    else if (c==blank) {
      result = 0;
      if (non_blank) {
	non_blank = 0;
	term = field;
      }
      len++;
      if (len<len_field)
	*(field++) = replace;
    }
    else {
      result = 0;
      non_blank = 1;
      len++;
      if (len>=len_field)
	return 0;
      *(field++) = (char) c;
    }
    if (!result)
      c = gzgetc(infile);
  }
  while (!result);
  if (!non_blank)
    *term = (char) 0;
  else
    *field = (char) 0;
return result;
}

/* Skip to end of line (or end of file) */

int skip_to_eol(const gzFile infile){
  const int eol = (int) '\n';
  const int eof = -1;
  while (1) {
    int c = gzgetc(infile);
    if (c==eol)
      return 2;
    if (c==eof)
      return 3;
  }
}

/* Recast nucleotide coding (1:10) to snp coding (1:3) */

int recode_snp(unsigned char *matrix, const int N, const int M) {
  int none_snps = 0;
  for (int s=0; s<M; s++, matrix+=N) {
    int count[11], lookup[11];
    memset(count, 0x00, 11*sizeof(int));
    memset(lookup, 0x00, 11*sizeof(int));
    for (int i=0; i<N; i++) {
      int m = (int) matrix[i];
      count[m]++;
    }
    int first=0, second=0, error=0;
    for (int j=1, ij=1; j<5; j++, ij++) {
      for (int i=1; i<j; i++, ij++) {
	if (count[ij]) { /* Hets observed */
	  if (second || (first && (first!=i))) {
	    error = 1;
	    break;
	  }
	  first = i;
	  second = j;
	  lookup[ij] = 2;
	}
      }
      if (error)
	break;
      if (count[ij]) {  /* Homs observed */
	if (!first) {
	  first = j;
	  lookup[ij] = 1;
	}
	else {
	  if (second && (second!=j)) {
	    error = 1;
	    break;
	  }
	  second = j;
	  lookup[ij] = 3;
	}
      }
    }
    if (error) {
      memset(matrix, 0x00, N*sizeof(unsigned char));
      none_snps++;
      warning("None-SNP in column %d", s+1);
    }
    else {
      for (int i=0; i<N; i++) {
	int m = (int) matrix[i];
	matrix[i] = (unsigned char) lookup[m];
      }
    }
  }
  return none_snps;
}
