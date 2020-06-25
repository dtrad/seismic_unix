#include "/usr/local/cwp/include/su.h"
#include "/usr/local/cwp/include/segy.h"
#include "/usr/local/cwp/src/Complex/include/Complex.h"
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
