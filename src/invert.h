/* Cholesky factorization */

int chol(double *, int, double *, int *, double *);

/* Inversion of positive semi-definite symmetric matrix */

int syminv(double *, int, double *, double *, int *, double *);

/* Inversion of triangular matrix */

int trinv(int, double *, double *);

/* Quadratic form U-transpose.V-inverse.U */

int qform(int, double *, double *, double *, double *, int *);
