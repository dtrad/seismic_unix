#include "/usr/local/cwp/include/su.h"
#include "/usr/local/cwp/include/segy.h"
#include "/usr/local/cwp/src/Complex/include/Complex.h"
// Add to  a square  matrix A(nxn) + Diag(n,n)
void AplusDiag(complex **C,complex **A,complex *Diag,int n)
{
     int i;
     
     for (i=0;i<n;i++)
	         C[i][i]=A[i][i]+Diag[i];
     return;
}
