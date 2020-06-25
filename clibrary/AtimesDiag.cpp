#include "su.h"
#include "segy.h"
#include "Complex.h"
// Multiply a rectangular matrix A(nrxnc) * Diag(nc,nc)

void AtimesDiag(complex **C,complex **A,float *Diag,int nr,int nc)
{
     int i,j,k;
     
     for (j=0;j<nc;j++)
         for (i=0;i<nr;i++)
	         C[i][j]=A[i][j]*Diag[j];
     return;
}

void AtimesDiag(complex **C,complex **A,complex *Diag,int nr,int nc)
{
     int i,j,k;
     
     for (j=0;j<nc;j++)
         for (i=0;i<nr;i++)
	         C[i][j]=A[i][j]*Diag[j];
     return;
}

