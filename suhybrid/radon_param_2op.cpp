#include "su.h"
#include "stddef.h"



void radon_param(float fmax,float *,int, float dx_av, float qmin, 
   float *qmaxt, float *qmax, float *dq, int nq, int rtmethod, float factor);
void interval(float *pos, int nh, float *dx_max, float *dx_av);


void radon_param_2op(float fmax1, float fmax2, float *h, int nh, float *q, float nq, 
		      float qmin1, float qmin2, float *q1, 
		      float *q2, int nq1, int nq2, float depth1, 
		      float depth2, int  rtmethod1, int rtmethod2, float factor1,
		      float factor2,float *pdq1, float *pdq2, int symmetricq1)
/*

   Given the field geometry (dx, xmin, xmax and fmax) and the chosen qmin,    
   it computes dq, maximum allowable qmax, nq.
    Hence, the NMO must be adjust such that q < qmax.	
    The dx can be the average or maximum, it is under research.
    rtmethod=1 LRT
    rtmethod=2 PRT 
    This version allows to define a different factor to undersample the q space
    For example when we want to go beyond the qmax defined by Nyquist
    Daniel Trad- UBC- 16-2-99

*/
{	
  float qmaxt;
  float qmax;
  float dx_max;
  float dx_av;
  float dq1;
  float dq2;
  int iq;

  fprintf(stderr,"qmin1=%e,qmin2=%e\n",qmin1,qmin2);

  interval(h,nh,&dx_max,&dx_av);
  fprintf(stderr,"dx_max=%f, dx_av=%f\n", dx_max, dx_av);
    
  radon_param(fmax1,h,nh,dx_av,qmin1,&qmaxt,&qmax,&dq1,nq1,rtmethod1,factor1);
    
  fprintf(stderr,"q max=%e,qmax used=%e\n", qmaxt,qmax);
  fprintf(stderr,"freq max=%f,dq=%e\n", fmax1,dq1);

  if (symmetricq1){
    for (iq=0;iq<nq1/2;iq++){
      q1[iq]=qmin1+iq*dq1;
      fprintf(stderr,"q[%d]=%f\n",iq,q1[iq]);
    }
    //fprintf(stderr,"q[%d]=%f\n",nq1/2-1,q1[nq1/2-1]);
    q1[nq1/2]=-q1[nq1/2-1];
    fprintf(stderr,"q[%d]=%f\n",nq1/2,q1[nq1/2]);
    for (iq=nq1/2+1;iq<nq1;iq++){
      q1[iq]=q1[iq-1]+dq1;
      fprintf(stderr,"q[%d]=%f\n",iq,q1[iq]);
    }
  }
  else{
    //qmin2=q1[nq1-1]+dq1;
    for (iq=0;iq<nq1;iq++){
      q1[iq]=qmin1+iq*dq1;
      fprintf(stderr,"q1[%d]=%f\n",iq,q1[iq]);
    }
  }
  
  radon_param(fmax2,h,nh,dx_av,qmin2,&qmaxt,&qmax,&dq2,nq1,rtmethod2,factor2);

  for (iq=0;iq<nq2;iq++) q2[iq]=qmin2+iq*dq2;

  for (iq=0;iq<nq1;iq++) q[iq]=q1[iq];
  for (iq=0;iq<nq2;iq++){
    q[iq+nq1]=q2[iq];
    fprintf(stderr,"q2[%d]=%e\n",iq,q2[iq]);
  }
  *pdq1=dq1;
  *pdq2=dq2;

  return;
}










