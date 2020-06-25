#include <time.h>
#include "dan.h"
#include "radonfkcube.h"

float adjteststoltzcube(int nt, int nh, int nv, float *t, float *h, float *vel)
{
  
  float **dr1;
  float ***mr1;
  float **dr2;
  float ***mr2;
  float dp1;
  float dp2;
  int it,ih, iv;
  float test;
  int ny=nt*nh;
  int nx=nt*nh*nv;
  unsigned int seed;	

  dr1=ealloc2float(nt,nh);
  dr2=ealloc2float(nt,nh);
  mr1=ealloc3float(nt,nh,nv);
  mr2=ealloc3float(nt,nh,nv);
  

  // Set seed
  if (-1 == (seed = (unsigned int) time((time_t *) NULL))) {
    err("time() failed to set seed");
  }
  srannor(seed);

  for (ih=0;ih<nh;ih++) for (it=0;it<nt;it++) dr1[ih][it]=frannor();

  // Change seed
  if (-1 == (seed = (unsigned int) time((time_t *) NULL))) {
    err("time() failed to set seed");
  }
  srannor(seed);

  for (iv=0;iv<nv;iv++)
    for (ih=0;ih<nh;ih++)
      for (it=0;it<nt;it++) 
	mr1[iv][ih][it]=frannor();


  stoltzop2(dr1,mr2,nt,nh,nv,t,h,vel,1);
  stoltzop2(dr2,mr1,nt,nh,nv,t,h,vel,0);

  dp1=dot(ny,dr1[0],dr2[0]);
  dp2=dot(nx,*mr1[0],*mr2[0]);
  fprintf(stderr,"dp1 = %f, dp2=%f \n",dp1,dp2);

  if (dp2!=0) test=dp1/dp2;
  else test=0;

  fprintf(stderr,"++++Test adjoint = %f \n",test);

  free3float(mr1);
  free3float(mr2);
  free2float(dr1);
  free2float(dr2);
  return(test);  

}






