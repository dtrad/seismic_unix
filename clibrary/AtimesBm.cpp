#include "su.h"
#include "segy.h"
#include "Complex.h"
void AtimesBm(complex **C,complex **A,complex **B,int nr,int ni,int nc)
{
     int i,j,k;
     for (i=0;i<nr;i++){
         for (j=0;j<nc;j++){
             C[i][j].r=0;
             C[i][j].i=0;
             for (k=0;k<ni;k++)
	         C[i][j]=C[i][j]+A[i][k]*B[k][j];
          }
     } 
     return;
}

void diagxtimesA(complex **B, float *x, complex **A, int nr, int nc)
{
     int i,j;
     for (i=0;i<nr;i++)
       for (j=0;j<nc;j++)
	 B[i][j]=x[i]*A[i][j];
     return;
}



void Atimes_conjB_elem(complex **C,complex **A,complex **B,int n1, int n2)
{
     int i,j;
     for (i=0;i<n2;i++)
         for (j=0;j<n1;j++)
             C[i][j]=A[i][j]*conjg(B[i][j]);

     return;
}

void Atimes_B_elem(complex **C,complex **A,complex **B,int n1, int n2)
{
     int i,j;
     for (i=0;i<n2;i++)
         for (j=0;j<n1;j++)
             C[i][j]=A[i][j]*B[i][j];

     return;
}


void AtimesBm(complex **C,complex **A,complex **B,int nr,int ni,int nc, 
	      char ch)
{
     int i,j,k;
     // This version computes only lower
     //if (ch=='t')
     for (i=0;i<nr;i++){
         for (j=i;j<nc;j++){
             C[i][j].r=0;
             C[i][j].i=0;
             for (k=0;k<ni;k++)
	         C[i][j]=C[i][j]+A[i][k]*B[k][j];
          }
     } 
     return;
}
void AtimesBm(complex **C,complex **A,complex **B,int nr,int ni,int nc, 
	      complex *D)
{
     int i,j,k;
     // This version computes only lower
     //if (ch=='t')
     for (i=0;i<nr;i++){
         for (j=0;j<nc;j++){
             C[i][j].r=0;
             C[i][j].i=0;
             for (k=0;k<ni;k++)
	         C[i][j]=C[i][j]+A[i][k]*D[k]*B[k][j];
          }
     } 
     return;
}
void AtimesBm(complex **C,complex **A,complex **B,int nr,int ni,int nc, 
	      complex *D, char ch)
{
     int i,j,k;
     // This version computes only lower
     //if (ch=='t')
     for (i=0;i<nr;i++){
         for (j=i;j<nc;j++){
             C[i][j].r=0;
             C[i][j].i=0;
             for (k=0;k<ni;k++)
	         C[i][j]=C[i][j]+A[i][k]*D[k]*B[k][j];
          }
     } 
     return;
}

void AtimesBH(complex **C,complex **A,complex **B,int nr1, int nc1, int nr2, int nc2)
{
  // C= A x B^H (H means Hermitian) 
  int i,j,k, ni=0;
  int nr=nr1;
  int nc=nr2;
  if (nc1==nc2) ni=nc1;
  else err("nr1 must be equal to nr2");
  

  for (i=0;i<nr;i++){
    for (j=0;j<nc;j++){
      C[i][j].r=0;
      C[i][j].i=0;
      for (k=0;k<ni;k++)
	C[i][j]=C[i][j]+A[i][k]*conjg(B[j][k]);
    }
  } 
  return;
}


void AHtimesB(complex **C, complex **A, complex **B, int nr1, int nc1, int nr2, int nc2)
{
  // C= A^H x B (H means Hermitian) 
  int i,j,k, ni=0;
  int nr=nc1;
  int nc=nc2;
  if (nr1==nr2) ni=nr1;
  else err("nr1 must be equal to nr2");
  


  for (i=0;i<nr;i++){
    for (j=0;j<nc;j++){
      C[i][j].r=0;
      C[i][j].i=0;
      for (k=0;k<ni;k++)
	C[i][j]=C[i][j]+conjg(A[k][i])*B[k][j];
    }
  } 
  return;
}

void transpose(complex **F, complex **FT, int nr, int nc)
{
  int i,j;

  for (i=0;i<nr;i++)
    for (j=0;j<nc;j++)
      FT[j][i]=conjg(F[i][j]);
  return;
}


