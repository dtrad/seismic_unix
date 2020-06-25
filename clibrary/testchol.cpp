#include "su.h"
#include "segy.h"
#include "Complex.h"
void displayA(complex **A,int nr,int nc);
float testchol(complex **A,double *p,int nr,int nc)
{
     int i,j,k;
     complex **L, **B, temp, czero;
     czero.r=czero.i=0;

     if ((L=alloc2complex(nc,nr))==NULL)
          err("cannot allocate memory for L\n");
     if ((B=alloc2complex(nc,nr))==NULL)
          err("cannot allocate memory for L\n");
     int nq=nr;
     for (i=0;i<nq;i++){
         for (j=0;j<i;j++){
             L[i][j]=A[i][j];
             L[j][i]=czero;
             B[j][i]=A[j][i];
             B[i][j]=conjg(A[j][i]);            
          }
     }
     for (i=0;i<nq;i++)
	     L[i][i]=p[i];
     
     //displayA(L,nq,nq);
     
  
     float res=0;
     for (i=0;i<nq;i++){
         for (j=0;j<nq;j++){
             temp=czero;
             for (k=0;k<nq;k++)
                   temp=temp+L[i][k]*conjg(L[j][k]);
             res=res+abs(B[i][j]-temp);//abs(B[i][j]);
             /*fprintf(stderr,"B[%d][%d].r=%e\n",i,j,B[i][j].r);
             fprintf(stderr,"temp.r=%e\n",temp.r);
             fprintf(stderr,"B[%d][%d].i=%e\n",i,j,B[i][j].i);
             fprintf(stderr,"temp.i=%e\n",temp.i);*/          
          }
     }

     res/=(nq*nq/2);

     free2complex(B);
     free2complex(L);

     return(res);
}





