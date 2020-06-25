#include "su.h"
void interval(float *x,int lx,float *pmx, float *pax)
{
   int i;
   float dx, mx, ax;		
   mx=0;
   ax=0;	
     for (i=1;i<lx;i++){
	 dx=fabs(x[i]-x[i-1]);
	 if (dx>mx) 
	        mx=dx;
	 ax=ax+dx;
     }
     ax=ax/(lx-1);
     *pmx=mx;
     *pax=ax;
     return;
}
