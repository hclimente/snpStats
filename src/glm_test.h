/* Family */

#define BINOMIAL  1
#define POISSON   2
#define GAUSSIAN  3
#define GAMMA     4

/* Link */

#define LOGIT     1
#define LOG       2
#define IDENTITY  3
#define INVERSE   4

/* GLM definition functions */

double varfun(int, double);
int muvalid(int, double);
double linkfun(int, double);
double invlink(int, double);
double dlink(int, double);



/* Fit a base model */

int glm_fit(int family, int link, int N, int M, int P, int S,
	    const double *y, const double *prior, const double *X, 
	    const int *stratum, int maxit, double conv, int init, 
	    int *rank, double *Xb, 
	    double *fitted, double *resid, double *weights, 
	    double *scale, int *df_resid,
	    int *P_est, int *which, double *beta, double *tri);

/* Score test for additional terms */

void glm_score_test(int N, int M, int S, const int *stratum, 
		    int P, const double *Z, int C, const int *cluster,
		    const double *resid, const double *weights, 
		    const double *Xb, double scale,
		    double max_r2, double *U, double *V);

/* Parameter estimation */

void glm_est(int P_est, const double *betaQ, double *tri, 
	     double scale, const double *meatrix, 
	     double *beta, double *var_beta);
 
void meat_matrix(int N, int P, int C, const int *cluster,
		 const double *Xb, const double *resid, const double *weights,
		 double *meatrix);

