#include <stdio.h>

void put_name(FILE *file, char *name, int quote) {
  if (quote) {
    fputc('\"', file);
    fputs(name, file);
    fputc('\"', file);
  }
  else
    fputs(name, file);
}

void write_as_matrix(char **file, char *x, int *N, int *M, 
		     char **rnames, char **cnames, int *asalleles, int *append, 
		     int *quote, char **sep, char **eol, char **na, 
		     int *row_names, int *col_names, int *iferror) {
  int zerom = (int) '0' -1 ; /* Character code for zero ... minus 1 */
  int nrow = *N;
  int ncol = *M;
  FILE *  outfile;
  int i=0, j=0, ij=0;
  if (*append)
    outfile = fopen(*file, "a");
  else
    outfile = fopen(*file, "w");
  if (!outfile) {
    *iferror = 1;
    return;
  }
  if (*col_names) {
    for (i=0; i<ncol; i++) {
      if (i)
	fputs(*sep, outfile);
      put_name(outfile, cnames[i], *quote);
    }
    fputs(*eol, outfile);
  }
  for (i=0; i<nrow; i++) { 
    if (*row_names) {
      put_name(outfile, rnames[i], *quote);
      fputs(*sep, outfile);
    }
    for (j=0, ij=i; j<ncol; j++, ij+=nrow) {
      if (j)
	fputs(*sep, outfile);
      int  g = (int) x[ij];
      if (*asalleles) {
	if (!g) {
	  fputs(*na, outfile);
	  fputs(*sep, outfile);
	  fputs(*na, outfile);
	}
	else {
	  fputc(g<3? '1': '2', outfile);
	  fputs(*sep, outfile);
	  fputc(g<2? '1': '2', outfile);
	}
      }
      else {
	if (!g)
	  fputs(*na, outfile);
	else {
	  g += zerom;
	  fputc((char) g, outfile);
	}
      }
    }
    fputs(*eol, outfile);
  }
  fclose(outfile);
  *iferror = 0;
  return;
}

    

