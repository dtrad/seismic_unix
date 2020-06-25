//#include "su.h"
//#include "cwp.h"
#include <time.h>
#include "dan.h"
#include "interpfk.h"

void contruc_2(int conj,int add, int nx, float *xx, int nb, float *bb,float *yy);

void convgather(float *d, float *m, int nh, int nt, float *wavelet, int nw, int adj)
{
  int ih;//, it;
  float *dtemp;
  float *dtemp2;

    
  fprintf(stderr,"nw=%d........\n",nw);
  dtemp=ealloc1float(nt+nw);
  dtemp2=ealloc1float(nt+nw);
  if (!adj){
    if (nw<=1) {
      for (ih=0;ih<nh;ih++)
	memcpy((void *) &d[ih*nt],(const void *) &m[ih*nt],nt*sizeof(float));
    }
    else{
      for (ih=0;ih<nh;ih++){
	memset((void *) dtemp,(int)'\0',(nt+nw)*FSIZE); 
	contruc_2(0,0,nw,wavelet,nt,&m[ih*nt],dtemp);
	memcpy((void *) &d[ih*nt],(const void *) dtemp,nt*sizeof(float));
      }
    }
  }
  else{ // if adj
    if (nw<=1) {
      for (ih=0;ih<nh;ih++)
	memcpy((void *) &m[ih*nt],(const void *) &d[ih*nt],nt*sizeof(float));
    }
    else{
      for (ih=0;ih<nh;ih++){
	memset((void *) dtemp,(int)'\0',(nt+nw)*FSIZE); 
	memset((void *) dtemp2,(int)'\0',(nt+nw)*FSIZE); 
	memcpy((void *) dtemp2,(const void *) &d[ih*nt],nt*sizeof(float));
	//for (it=0;it<nw;it++) dtemp2[it]=0;
	contruc_2(1,0,nw,wavelet,nt,dtemp,dtemp2);
	memcpy((void *) &m[ih*nt],(const void *) dtemp,nt*sizeof(float));
      }
    }
  }
  free1float(dtemp2);
  free1float(dtemp);
  return;
}


float adjtestconvgather(int nh, int nt, float *wavelet, int nw)
{
  
  float **dr1;
  float **mr1;
  float **dr2;
  float **mr2;
  float dp1;
  float dp2;
  int it,ih,iq;
  float test;
  int nq=nh; // for this operator nq=nh always
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


  convgather(dr1[0],mr2[0],nh,nt,wavelet,nw,1);
  convgather(dr2[0],mr1[0],nh,nt,wavelet,nw,0);

  dp1=dot(ny,dr1[0],dr2[0]);
  dp2=dot(nx,mr1[0],mr2[0]);
  fprintf(stderr,"dp1 = %f, dp2=%f \n",dp1,dp2);

  if (dp2!=0) test=dp1/dp2;
  else test=0;

  fprintf(stderr,"++++Test adjoint convgather = %f \n",test);

  free2float(mr1);
  free2float(mr2);
  free2float(dr1);
  free2float(dr2);
  return(test);  

}













