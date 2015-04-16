#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <R.h>

#include "hphase.h"

int ipf(int K, const double *observed, 
	const int nterms, const unsigned int *terms,  double *expected,
	const int maxit, const double eps);

/* Generate genotype->haplotype lookup table */

GTYPE *create_gtype_table(const int nsnp) {
  int ngt = (1 << (2*nsnp)) - 1;  /* 4^nsnp - 1 */
  GTYPE *result = (GTYPE *)Calloc(ngt, GTYPE);
  if (!result)
    return NULL;
  int *alleles = (int *)Calloc(nsnp, int);
  memset(alleles, 0x00, nsnp*sizeof(int));
  int igt=0;
  while(1) {
    /* Advance genotype by one */
    int adv = 0;
    while (adv<nsnp) {
      if (++(alleles[adv])==4) 
	alleles[adv++] = 0;
      else 
	break;
    } 
    if (adv==nsnp)
      break;

    /* Calculate number of phased assignments for this genotype */
  
    int hom=1, het=0;
    for (int i=0; i<nsnp; i++) {
      int ai = alleles[i];
      if (!ai) {
	het = hom + 4*het; /* Effect of adding 1/2 genotype at this SNP */
	hom = 2*hom; /* Effect of adding 1/1 and 2/2 genotypes */
      }
      else if (ai==2) {
	het = hom + 2*het; /* Effect of adding 1/2 genotype at this SNP */
	hom = 0; /* No homozygous genotypes left */
      }
    }
    int nph = hom + het;

    /* Create element in output array */

    result[igt].nphase = nph;
    int *haps = (int *)Calloc(2*nph, int);
    if (!haps)
      return NULL;
    result[igt].haps = haps;
    igt++;

    /* Compute haplotype assignments */

    haps[0] = haps[1] = 0;
    for (int i=0, iph=1, one=1; i<nsnp; i++) {
      int ai = alleles[i];
      /* iph is the number of assignments for the partial genotype */
      int add = 0;
      for (int j=0, j1=0, k=2*iph; j<iph; j++) {
	int h1 = haps[j1];
	int h2 = haps[j1+1];
	if (!ai) {
	  haps[j1++] = h1;
	  haps[j1++] = h2;
	  haps[k++] = h1 | one;
	  haps[k++] = h2 | one;
	  haps[k++] = h1;
	  haps[k++] = h2 | one;
	  add+=2;
	  if (h1!=h2) {
	    add++;
	    haps[k++] = h1 | one;
	    haps[k++] = h2;
	  }	    
	}
	else if (ai==1) {
	  haps[j1++] = h1;
	  haps[j1++] = h2;
	}
	else if (ai==2) {
	  haps[j1++] = h1;
	  haps[j1++] = h2 | one;
	  if (h1!=h2) {
	    add++;
	    haps[k++] = h1 | one;
	    haps[k++] = h2;
	  }
	}
	else if (ai==3) {
	  haps[j1++] = h1 | one;
	  haps[j1++] = h2 | one;
	}    
      }
      iph += add;
      one = (one<<1);
    }
  }
  Free(alleles); 
  return result;
}  

/* Destroy genotype->haplotype lookup table */
      
void destroy_gtype_table(GTYPE *gtt, const int nsnp){
  int ngt = (1 << (2*nsnp)) - 1;  /* 4^nsnp - 1 */
  for (int i=0; i<ngt; i++) {
    Free(gtt[i].haps);
  }
  Free(gtt);
}


/* EM algorithm to calculate haplotype frequencies from tables of unphased 
   genotype frequencies  plus an optional table of phased haplotype frequencies

   nsnp    Number of SNPs
   gtable  Table of dimension 4^nsnp containing unphased genotype frequencies
   htable  Table of same dimension containing male X-chromosome data
   gtypes  (Optional) genotype->haplotype lookup table 
   maxit   Maximum number of EM iterations
   tol     Tolerance for proportional change in fitted haplotype probabilities
   hprob   Output table of dimension 2^nsnp containing haplotype probabilities
   nllm    Number of terms in log-linear smoothing model
   llm     Model terms (as bit patterns)
   
*/


int emhap(const int nsnp, const int *gtable, const int *htable, 
	  GTYPE *gtypes, const int maxit, const double tol,
	  double *hprob, const int nllm, const unsigned int *llm){

  GTYPE *lookup;

  /* Minimum number of EM steps before warning messages */

  int minit = 3;
  
  /* IPF control */

  const int ipfsteps = 10; /* Max number of IPF steps in each EM step */
  const double ipfeps = 0.001; /* IPF criterion (relative change in expected) */

  /* If lookup table not supplied, create one */
  
  if(gtypes) 
    lookup = gtypes;
  else
    lookup = create_gtype_table(nsnp);


  /* Table sizes */

  int ngt = (1 << (2*nsnp));
  int nht = (1 << nsnp);


  /* Count of haplotypes  */

  int npu = 0, npk = 0;
  for (int i=1; i<ngt; i++) {
    npu += gtable[i];
    if (htable)
      npk += htable[i];
  }
  npu *= 2;

  double total = npu + npk;
  if (!total) 
    return -1;

  /* Maximum number of possible haplotype assignments */

  int maxhaps = (1 << 2*(nsnp - 1));

  /* Work arrays */
  
  double *sum = (double *)Calloc(nht, double);
  double *prg = (double *)Calloc(maxhaps, double);
  double *prh = NULL;
  if (htable)
    prh = (double *)Calloc(maxhaps, double);

  /* If no starting values, initialize haplotype frequency vector */

  if (hprob[0]<0.0) {
    double maxp = 1.0/ (double)nht; 
    for (int i=0; i<nht; i++)
      hprob[i] = maxp;
  }

  /* EM algorithm */

  int it = 0, result = 0;
  double logL_prev = 0.0;
  while(1) {
    memset(sum, 0x00, nht*sizeof(double));
    double logL = 0.0;
    for (int i=1; i<ngt; i++) {
      int gti = gtable[i];
      int hti = htable? htable[i]: 0;
      if (gti || hti) {
	double psumg = 0.0, psumh = 0.0;
	GTYPE *lupi = lookup + i - 1;
	int nph = lupi->nphase;

	/* Posterior probabilities */

	for (int j=0, jj=0; j<nph; j++) {
	  int h1 = lupi->haps[jj++];
	  int h2 = lupi->haps[jj++];
	  if (gti) {
	    double p = hprob[h1]*hprob[h2];
	    if (h1!=h2)
	      p *= 2;
	    prg[j] = p;
	    psumg += p;
	  }
	  if (hti && (h1==h2)) { /* Males are coded as homozygous on X */
	    double p = hprob[h1];
	    prh[j] = p;
	    psumh += p;
	  }
	}
	if (gti)
	  logL += gti*log(psumg);
	if (hti)
	  logL += hti*log(psumh);

	/* Increment haplotype table with expected frequencies */
	
	double fgi = psumg? gtable[i]/psumg: 0.0;
	double fhi = psumh? htable[i]/psumh: 0.0;
	if (fgi || fhi) {
	  for (int j=0, jj=0; j<nph; j++){
	    int h1 = lupi->haps[jj++];
	    int h2 = lupi->haps[jj++];
	    if (fgi) {
	      double fij = fgi*prg[j];
	      sum[h1] += fij;
	      sum[h2] += fij;
	    }
	    if (fhi && (h1==h2)) {
	      double fij = fhi*prh[j];
	      sum[h1] += fij;
	    }
	  }
	}
      }
    }

    /* New estimates of haplotype frequencies */

    int ipfault = 0;
    if (nllm) { /* Log-linear smoothing model */
      for (int i=0; i<nht; i++) {
	sum[i] /= total;
	/* do ipfsteps of IPF algorithm */
	ipfault = ipf(nsnp, sum, nllm, llm, hprob, ipfsteps, ipfeps);
      }
    }
    else { /* Saturated model */
      for (int i=0; i<nht; i++) 
      hprob[i] = sum[i]/total;
    }

    /* Convergence test */

    double ctest = logL - logL_prev;
    logL_prev = logL;
    if (it++) {
      if (it>minit && ctest<0.0) {
	warning("Log likelihood decreased in EM algorithm at iteration %d", it);
	result = -2;
	break;
      }
      else if (it>maxit) {
	result = 1;
	break;
      }
      else if (ctest < tol) {
	break;
      }
    }
  }

  /* Return work arrays */

  if(!gtypes)
    destroy_gtype_table(lookup, nsnp);
  Free(sum);
  Free(prg);
  if (prh)
    Free(prh);

  return result;
}

/* Predict SNP from  SNP genotype

   npr    number of predictor SNPs
   g      predictor genotype (single code, npr SNPs)
   mX     1 if haploid, otherwise 0
   hprob  haplotype probs (table of dimension npr+1 with SNP to be 
          predicted varying fastest)
   gtypes genotype->haplotype lookup table for genotype table of dimension npr

   pred   (length 3) posterior probabilities for genotype == 0, 1, 2
*/

void predict_gt(const int npr, const int g, const int mX, const double *hprob, 
		const GTYPE *gtypes, double *pred) {
  if (!g) {
    pred[0] = pred[1] = pred[2] = NA_REAL;
    return;
  }
  int nhaps = gtypes[g-1].nphase;
  int *haps = gtypes[g-1].haps;
  /* Loop over possible phased haplotype assignments */
  double psum = 0.0, qsum = 0.0, qprod = 0.0;
  for (int i=0, ii=0; i<nhaps; i++) {
    int h1 = 2*haps[ii++];
    int h2 = 2*haps[ii++];
    if (mX) {
      if (h1==h2) { /* Only consider homozygous codings */
	double p0 = hprob[h1];
	double p1 = hprob[h1+1];
	double p = p0+p1;
	psum += p;
	qsum += p1;
      }
    }
    else {
      double p10 = hprob[h1];
      double p11 = hprob[h1+1];
      double p1 = p10+p11;
      double q1 = p11/p1; /* Conditional probability h1 carries 1 */
      double p20 = hprob[h2];
      double p21 = hprob[h2+1];
      double p2 = p20+p21;
      double q2 = p21/p2;/* Conditional probability h2 carries 1 */
      /* Probability of this haplotype assignment */
      double p = p1*p2;
      if (h1!=h2)
	p *= 2;
      psum += p;
      if (p) {
	qsum += p*(q1+q2);
	qprod += p*q1*q2;
      }
    }
  }
  /* Posterior probabilities of genotypes */
  if (psum > 0.0) {
    if (mX) {
      double pg2 = qsum/psum;
      pred[2] = pg2; 
      pred[1] = 0.0;
      pred[0] = 1.0 - pg2;
    }
    double pg1 = (qsum - 2.0*qprod)/psum;
    double pg2 =  qprod/psum;
    pred[0] = 1.0 - pg1 - pg2;
    pred[1] = pg1;
    pred[2] = pg2;
  }
  else {
    pred[0] = pred[1] = pred[2] = NA_REAL;
  }
}


/* Predict allele on a haplotype 

   npr    number of predictor SNPs
   h      predictor haplotype (npr SNPs)
   hprob  haplotype probs (table of dimension npr+1 with SNP to be 
          predicted varying fastest)

   Returns probability that predicted SNP has the second allele

*/


double predict_ht(const int h, const double *hprob) {
  int h0 = 2*h;
  double p1 = hprob[h0+1];
  double den = hprob[h0] + p1;
  return den>0? p1/den: NA_REAL;
}

/* r^2 for prediction of first SNP. hprob is assumed to sum to 1.0 */

double hap_r2(const int npr, const double *hprob){
  if (npr<1)
    return -1;
  int nhp = (1 << npr);
  double s0 = 0.0, s1 = 0.0;
  for (int i=0, ii=0; i<nhp; i++) {
    double f0 = hprob[ii++];
    double f1 = hprob[ii++];
    double f = f0 + f1;
    if (f) {
      s0 += f1;
      s1 += f1*f1/f;
    }
  }
  double var1 = s0*(1.0-s0);
  double var2 = (s1 - s0*s0);
  return var2/var1;
}

/* Genotype version. gtypes is gtype->htype table for nsnp-1 SNPs */

double gen_r2(const int npr, double *hprob, const GTYPE *gtypes) {

  int *alleles = (int *)Calloc(npr, int);
  int g = 0;
  for (int i=0, inc=1; i<npr; i++) {
    alleles[i] = 0;
    g += inc;
    inc = inc << 2;
  }
  /* Cycle through complete genotypes */
  double sp = 0.0, mu =0.0, vp=0.0;
  double pr[3];
  while(1) {    
    predict_gt(npr, g, 0, hprob, gtypes, pr);
    double yp = pr[1] + 2.0*pr[2];
    double p = pr[0];
    sp += p;
    if (p) {
      double pyp = p*yp;
      mu += pyp;
      vp += pyp*yp;
    }
    /* Advance the genotype */
    int adv = 0, inc=1;
    g++;
    while (adv<npr) {
      if (++(alleles[adv])==3) {
	g += inc;
	alleles[adv] = 0;
	adv++;
	inc = inc << 2;
      }
      else 
	break;
    } 
    if (adv==npr)
      break;
  }
  Free(alleles);
  mu /= sp;
  vp = vp/sp - mu*mu;
  double vo = mu * (1.0 - mu/2.0);
  return vp/vo;
}

