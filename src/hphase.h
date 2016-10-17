/* Element of genotype->haplotype lookup table */

typedef struct {
  int nphase;  /* The number of possible phases */
  int *haps;   /* 2*nphase array of haplotype-pair assignments */
} GTYPE;

/* 
   Genotypes are numbered with first allele varying fastest starting 
   from 1,0,0... e.g.
   1,0,0,... = 0
   2,0,0,... = 1
   3,0,0,... = 2
   0,1,0,... = 3
   For m SNPs, genotypes are coded 0:(4^m-2) with 0 meaning missing; i.e.
   g = -1 + alleles[0] + 4*alleles[1] + 16*alleles[2] + ...
   and -1 is ignored (representing completely missing data).

   Haplotypes are also numbered with first SNP varying fastest. For 
   m SNPs, haplotypes are coded 0:(2^m-1)
*/

/* Create lookup table */
GTYPE *create_gtype_table(const int nsnp); 
/* Destroy lookup table */
void destroy_gtype_table(GTYPE *gtt, const int nsnp); 
/* EM algorithm */
int emhap(const int nsnp, const int *gtable, const int *htable, 
	  GTYPE *gtypes, const int maxit, const double tol,
	  double *hprob, const int nllm, const unsigned int *llm);
/* Predicted genotype scores */
void predict_gt(const int npr, const int g, const int mX, const double *hprob, 
		const GTYPE *gtypes, double *pred);
/* Predict allele on a haplotype */
double predict_ht(const int h, const double *hprob);
/* Expected r^2 for prediction of SNP 1 from SNPs 2,...,nsnp */
double hap_r2(const int npr, const double *hprob);
double gen_r2(const int npr, double *hprob, const GTYPE *gtypes);
