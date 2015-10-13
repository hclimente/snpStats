/* Note that these functions deal with triangular matrices stored in the
   following order:
                           0  1  3  6  ... etc.
                              2  4  7
                                 5  8
                                    9                               */
 
#include <math.h>
#include <stdlib.h>


int chol(double a[], int n, double u[], int *nullty, double *log_det)
/* Based on Algorithm AS6 J.R.S.S.(C), 1968, Vol 17, No 2. */
{
     double  ldet, w, etaak, eta=0.000001;
     int nty, i, j, k, l, m, icol, irow;

     if (n<=0)
	 return(1);
     for (ldet=0.0, nty=j=k=icol=0; icol<n; icol++, j=k)
         {
         for (l=irow=0; irow<=icol; irow++, k++)
             {
             for (w=a[k], m=j, i=0; i<irow; i++, l++, m++)
                 w=w-u[l]*u[m];
             if (irow==icol)
                 {
                 etaak=eta*a[k];
		 if (w>etaak)
		     {
		     ldet += log(w);
		     u[k]=sqrt(w);
		     }
                 else if (w<-etaak)
                     return(2);
                 else
                     {
                     u[k]=0.0;
                     nty++;
                     }
                 }
             else
                 {
                 if (u[l]==0.0)
                     u[k]=0.0;
                 else
                     u[k]=w/u[l];
                 }
             l++;
             }
         }
     *nullty = nty;
     *log_det = ldet;
     return(0);
}


int syminv(double a[], int nrow, double c[], double w[], int *nullty,
       double *log_det)
/* Based on Algorithm AS7 J.R.S.S.(C), 1968, Vol 17, No 2. */
{
    double x;
    int i, j, k, l, irow, icol, jcol, mdiag, ndiag, last_cell, last_row;

    if (nrow<=0)
        return(1);
    i=chol(a,nrow,c,nullty,log_det);
    if (i!=0)
        return(i);
    last_cell=((nrow*(nrow+1))/2) - 1;
    last_row=nrow-1;
    for (ndiag=last_cell, irow=last_row; irow>=0; ndiag-=irow, ndiag--, irow--)
        {
        if (c[ndiag]!=0.0)
            {
            for (l=ndiag, i=irow; i<nrow; i++, l+=i)
                w[i]=c[l];
            for (icol=last_row, mdiag=jcol=last_cell; icol>=irow;
                      mdiag-=icol ,mdiag-- , icol--, jcol--)
                {
                if (icol==irow)
                    x=1.0/w[irow];
                else
                    x=0.0;
                for (k=last_row, l=jcol; k>irow; k--)
                    {
                    x-= (w[k]*c[l]);
                    if (l>mdiag)
                        l-=k;
                    else
                        l--;
                    }
                c[l]=x/w[irow];
                }
            }
        else
            for (j=irow, l=ndiag; j<nrow; j++, l+=j)
                c[l]=0.0;
        }
    return(0);
}

int trinv(int n, double *u, double *c)
/* Invert n*n upper triangular matrix, u[], store in c[] (c may overwrite u).
   Returns Nullity - the number of zero diagonals. */
{
   int nullty, i, j, k, ij, ik, jj, jk;
   double uii, cij;

   for (nullty=i=ij=0; i<n; i++){
      uii = u[ij+i];
      if (uii != 0.0){
         for(j=jj=0; j<i; jj += ++j+1){
            for(k=j, ik=ij, jk=jj, cij=0.0; k<i; ik++, jk += ++k)
               cij += c[jk]*u[ik];
            c[ij++] = -cij/uii;
         }
         c[ij++] = 1/uii;
      }else{
         for(j=0; j<=i;j++, ij++) c[ij] = 0.0;
         nullty ++;
      }
   }
   return nullty;
}
   
/* Quadratic form U-t.V-inv.U. W is work space same size as V but may 
   coincide with V (which is then destroyed) */

int qform(int N, double *U, double *V, double *W, double *quad, int *rank)
/* Quadratic form U-transpose.V-inverse.U */
{
  int dynW = 0;
  if (!W) {
    dynW = 1;
    W = (double *) calloc((N*(N+1))/2, sizeof(double));
  }
  int nullty;
  double log_det;
  
  int ifault = chol(V, N, W, &nullty, &log_det);
  if (ifault)
    return(ifault);
  nullty = trinv(N, W, W);
  double res=0.0;
  for (int i=0, ij=0; i<N; i++) {
    double w=0.0;
    for (int j=0; j<=i; j++, ij++) {
      w += U[j]*W[ij];
    }
    res += w*w;
  }
  *rank = N - nullty;
  *quad = res;
  if (dynW)
    free(W);
  return(0);
}
	
