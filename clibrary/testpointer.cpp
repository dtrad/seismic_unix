#include "su.h"
#include "clibrarytd.h"
void testpointer(float **d, int nt, int nh)
{
  float **dp;
  int ih, it;

  if ((dp=(float **) alloc1(nt,sizeof(float*)))==NULL) 
    err("Cannot assign mem to dp"); 
  dp[0]=(float *)d[0];

  //for (ih=0;ih<nh;ih++) for (it=1;it<nt;++it) dp[it][ih]=&d[ih][it];

  // for (ih=0;ih<nh;ih++)
  //for (it=0;it<nt;it++)
      //fprintf(stderr,"dp[it][ih]=%f,d[ih][it]=%f\n",dp[it][0],d[0][it]);
      

  free((void **) dp);  
  return;

}

