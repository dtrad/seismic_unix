#include "su.h"
#include "segy.h"
#include "Complex.h"
// Add to  a square  matrix A(nxn) + Diag(n,n)
float maxmax(complex **A,int nr, int nc)
{
     int i, j;
     float cmax;
     cmax=abs(A[0][0]);
     for (i=0;i<nr;i++)
       for (j=0;j<nc;j++)
	         cmax=MAX(abs(A[i][j]),abs(cmax));
     return(cmax);
}
float maxmax(complex *A,int n)
{
     int i;
     float cmax;
     cmax=abs(A[0]);
     for (i=0;i<n;i++)
	         cmax=MAX(abs(A[i]),abs(cmax));
     return(cmax);
}
