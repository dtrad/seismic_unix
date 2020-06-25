#include <time.h>
#include "dan.h"
#include "interpfk.h"

float adjteststoltz(int nt, int nh, float *t, float *h, float vel,float **data1, 
		    float **model1, complex **F, complex **FLS)
{
  
  float **dr1;
  float **mr1;
  float **dr2;
  float **mr2;
  float dp1;
  float dp2;
  int it,ih,iq;
  float test;
  int nq=nh;
  int ny=nt*nh;
  int nx=nt*nq;
  unsigned int seed;	
  dr1=ealloc2float(nt,nh);
  dr2=ealloc2float(nt,nh);
  mr1=ealloc2float(nt,nq);
  mr2=ealloc2float(nt,nq);
  

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


  for (iq=0;iq<nq;iq++) for (it=0;it<nt;it++) mr1[iq][it]=frannor();

  //for (ih=0;ih<nh;ih++) for (it=0;it<nt;it++) dr1[ih][it]=data1[ih][it];
  //for (iq=0;iq<nq;iq++) for (it=0;it<nt;it++) mr1[iq][it]=model1[iq][it];

  stoltzop(dr1,mr2,nt,nh,t,h,vel,F,FLS,1);
  stoltzop(mr1,dr2,nt,nh,t,h,vel,F,FLS,0);

  dp1=dot(ny,dr1[0],dr2[0]);
  dp2=dot(nx,mr1[0],mr2[0]);
  fprintf(stderr,"dp1 = %f, dp2=%f \n",dp1,dp2);

  if (dp2!=0) test=dp1/dp2;
  else test=0;

  fprintf(stderr,"++++Test adjoint = %f \n",test);

  free2float(mr1);
  free2float(mr2);
  free2float(dr1);
  free2float(dr2);
  return(test);  

}






