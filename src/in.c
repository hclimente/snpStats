/* Read genotype data in long format

filename    Name of input text file (4 fields per row) 
tmpdir      Temporary directory for sorting (no trailing /)
threshold   Quality threshold for inclusion of data
nchip       Number of chips
chip_id     Array of chip identifier strings
nsnp        Number of snps
snp_id      Array of snp identifier strings
gtypes      Array of char's, length nchip*nsnp, to hold genotype data

*/

#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#define MAX_ID 128
#define MAX_GT 16

void insnp(char *filename, char *tmpdir,  
	   int *nchip, char **chip_id, int *nsnps, char **snp_id,
	   char *codes[3], double *threshold, char *gtypes, 
	   int counts[2], int *iferror) {
  /* Sort with chips varying fastest */
  char sort_command[160];
  int i = 0, j = 0;
  sprintf(sort_command, 
	  "sort  -k 2,2 -k 1,1 -T \"%s\" -o \"%s\" \"%s\"",
          tmpdir, filename, filename);
  printf("%s\n", sort_command);
  int error = system(sort_command);
  if (error) goto sort_error;
  FILE *infile = fopen(filename, "r");
  if (!infile) 
    goto open_error;

  int not_called = 0;
  int called = 0;
  char *code_aa = codes[0];
  char *code_ab = codes[1];
  char *code_bb = codes[2];
  char chip_in[MAX_ID], snp_in[MAX_ID], gt_in[MAX_GT];
  double thr_in;
  if (fscanf(infile, " %s %s %s %lf", chip_in, snp_in, gt_in, &thr_in)!=4)
    goto read_error;
  int ij=0;
  for (j=0; j<(*nsnps); j++) {
    char *snp_target = snp_id[j];
    int jcmp;
    while ((jcmp = strcmp(snp_in, snp_target))<0) {
      int scanned = fscanf(infile, " %s %s %s %lf", 
		       chip_in, snp_in, gt_in, &thr_in);
      if (scanned==EOF) goto normal;
      else if (scanned!=4)
	goto read_error;  
    }   
    for (i=0; i<(*nchip); i++, ij++) {
      char * chip_target = chip_id[i];
      int icmp; 
      if (!jcmp) {
	while ((icmp = strcmp(chip_in, chip_target))<0) { 
	  int scanned = fscanf(infile, " %s %s %s %lf", 
			       chip_in, snp_in, gt_in, &thr_in);
	  if (scanned == EOF) 
	    goto normal;
	  else if (scanned!=4)
	    goto read_error;
	}
	if (!icmp) {
	  /* Assign genotype */
	  if (!strcmp(code_aa, gt_in)) {
	    gtypes[ij] = (char) 1;
	    called ++;
	  }
	  else if  (!strcmp(code_ab, gt_in)) {
	    gtypes[ij] = (char) 2;
	    called ++;
	  }
	  else if  (!strcmp(code_bb, gt_in)) {
	    gtypes[ij] = (char) 3;
	    called ++;
	  }
	  else {
	    gtypes[ij] = (char) 0;
	    not_called ++;
	  }
	}
      }
      else {
	gtypes[ij] = (char) 0;
      }
    }
  }

  /* Normal return */

 normal: {
  int full = (*nsnps)*(*nchip);
  while (ij<full)
    gtypes[ij++] = (char)0;
  counts[0] = called;
  counts[1] = not_called;
  *iferror = 0;
  return;
  }
  

    /* Error conditions */ 
  
 sort_error: *iferror = 1; return;
 open_error: *iferror = 2; return;
 read_error: *iferror = 3; return;
 
}

int main() {
  char* chips[2] = {"1", "2"};
  char* snps[3] = {"a", "b", "c"};
  char* codes[3] = {"aa", "ab", "bb"};
  int counts[2] = {0,0};
  char res[6];
  int nchips = 2;
  int nsnps = 3;
  int iferror;
  double thresh=0.8;
  int i = 0;
  insnp("test.txt", "~/temp", &nchips, chips, &nsnps, snps, codes, &thresh, 
	res, counts, &iferror); 
  printf("iferror = %d, counts = %d, %d\n", iferror, counts[0], counts[1]);
  for (i=0; i<6; i++)
    printf("%-2o\n", res[i]);
  
  exit(0);
}
