#include "su.h"
#include "Complex.h"
#include "clibrary.h"
void chol_all(complex *R,complex **ll,complex *Qp,int n,complex *b,
              complex *x,double  *diagll)    
{
  // construct ll from R
       int i, j;
       for (i=0;i<n;i++)
	 for (j=i;j<n;j++)
            ll[i][j]=R[j-i];
      
	 
       for (i=0;i<n;i++) {ll[i][i]+=Qp[i];ll[i][i].i=0;}
       choldc(ll,n,diagll);

       cholsl(ll,n,diagll,b,x);
       return;
}





