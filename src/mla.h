int wcenter(const double *y, int n, const double *weight, const int *stratum, 
	    int nstrata, int resid, double *ynew);

double wresid(const double *y, int n, const double *weight, const double *x, 
	   double *ynew);

double wssq(const double *y, int n, const double *weight);

double wsum(const double *y, int n, const double *weight);

double wspr(const double *y, const double *x, int n, const double *weight);
