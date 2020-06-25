#include "su.h"
float dot(int n, float *a, float *b);

void adaptive_subtract(float **data, float **datapred, int nh, int nt)
{
  /* Calculate a scale factor that minimize the difference between
     data and datapred in a least squares sense */
  
  float dpd=0;
  float dpdp=0;
  float scale=0;
  int ih,it;

  fprintf(stderr,"substracting multiples \n");
  dpd=dot(nh*nt,data[0],datapred[0]);
  dpdp=dot(nh*nt,datapred[0],datapred[0]);
  scale=dpd/dpdp;
  fprintf(stderr,"scale ===>%f\n",scale);
  for (ih=0;ih<nh;ih++) 
    for (it=0;it<nt;it++)
      data[ih][it]=data[ih][it]-scale*datapred[ih][it];

  return;
}



 
 


