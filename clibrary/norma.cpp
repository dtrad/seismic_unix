#include "math.h"
//float abs(float);
void norma(float **x, int nt, int nh)
{
       register int it, ih;	 
       float  xmax;
       xmax=x[1][1];
       for (it=0;it<nt;it++)
	  for (ih=0;ih<nh;ih++)     
              if(fabs(x[it][ih])>xmax) xmax=fabs(x[it][ih]);
       for (it=0;it<nt;it++)
	  for (ih=0;ih<nh;ih++)     
              x[it][ih]=x[it][ih]/xmax;                  

       return;
}

