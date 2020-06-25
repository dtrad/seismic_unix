#include "/usr/local/cwp/include/su.h"
#include "/usr/local/cwp/include/segy.h"
#include "/usr/local/cwp/src/Complex/include/Complex.h"
void displayA(complex *A,int nr,int nc)
{
     int i,j,k;
     for (i=0;i<nr;i++){
         for (j=0;j<nc;j++){
           fprintf(stderr,"A.r[%d][%d]=%e,\t A.i[%d][%d]=%e\n",i,j,A[i][j].r
              ,i,j,A[i][j].i);
           fprintf(stderr,"A.r[%d][%d]=%e,\t A.i[%d][%d]=%e\n",j,i,A[j][i].r
              ,j,i,A[j][i].i);
          }
     } 
     return;
}






