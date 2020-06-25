#include <math.h>
#include <stdio.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
// Matlab
#include "C:\MATLAB5\EXTERN\INCLUDE\mex.h"

#define printf mexPrintf

void rstack( double *x, double *t, double *h, double *p, double *y, 
			int nt, int nh, int np, double dt) ;
double sinc(double arg);

void mexFunction( int nlhs, mxArray *plhs[],
				 int nrhs, const mxArray *prhs[])

{
	if (nrhs==0) {
		printf("[v]=radonopi(u,t,h,p)\n") ;
		printf("L hyperbolic radon operator\n") ;	
		return;
	}
	
	double *x, *z;
	int nt,nh,np;
	double *t,*p,*h, dt;

///////////////////////////////////
// Get input from MATLAB
	nt=mxGetM(prhs[0]);printf("nt=%d\t",nt);
	nh=mxGetN(prhs[0]);printf("nh=%d\t",nh);
	z=mxGetPr(prhs[0]);
	t=mxGetPr(prhs[1]);
	h=mxGetPr(prhs[2]);
	p=mxGetPr(prhs[3]);
	np=mxGetN(prhs[3]);printf("np=%d\n",np);
	dt=t[1]-t[0];printf("dt=%f\n",dt);
    //output			
	plhs[0]=mxCreateDoubleMatrix(nt,np, mxREAL);
	x=mxGetPr(plhs[0]);
////////////////////////////////	

	rstack(x,t,h,p,z,nt,nh,np,dt);
    	
	return;
}


void rstack( double *m, double *t, double *h, double *q, double *d, 
			int nt, int nh, int nq, double dt) 
{
  int i,k,ih,iq,itau;
  unsigned int j;
  double time,it;
  int nx, ny;
  nx=nt*nq;
  ny=nt*nh;
    
  for (i=0;i<(nq*nt);i++) m[i]=0;

  for (ih=0;ih<nh;ih++){
    for (iq=0;iq<nq;iq++){
	  for (itau=0;itau<nt;itau++){
	    k=iq*nt+itau;
	    time=sqrt(pow(t[itau],2)+pow(h[ih],2)*q[iq]);
		it=time/dt;
	    j=ih*nt+(int) floor(it);
	    if ((it!=nt)&&(j<ny)&&(k<nx)) 
			m[k]=m[k]+d[j];//*sincin(it-floor(it));
		
	  }
	  
    }
  }
  return;
}

