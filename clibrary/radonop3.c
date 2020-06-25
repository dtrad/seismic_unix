#include <math.h>
#include <stdio.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
/* Matlab */
#include "mex.h"

#define printf mexPrintf

void rstack( double *x, double *t, double *h, double *p, double *y, 
			int nt, int nh, int np, double dt, double *ww) ;
double sinc(double arg);

void mexFunction( int nlhs, mxArray *plhs[],
				 int nrhs, const mxArray *prhs[])

{
	double *x, *z, *ww;
	int nt,nh,np,flag;
	double *t,*p,*h, dt;
	if (nrhs==0) {
		printf("[v]=radonopi(u,t,h,p)\n") ;
		printf("LH Inverse hyperbolic radon operator\n") ;
		printf("u and v can be input and output as 1D vectors\n") ;
		printf("or as 2D vectors. u(nt,nh) gives v(nt,np) \n") ;
		printf("u(ntxnh,1) gives v(ntxnp,1) \n");
		printf("(time is fast dimension)\n") ;
		printf("to use in CG algorithms\n") ;
		printf("Daniel Trad - UBC\n") ;
		return;
	}
	


	/*/////////////////////////////////
	  // Get input from MATLAB*/
	/*nt=mxGetM(prhs[0]);printf("nt=%d\t",nt);*/
	nh=mxGetN(prhs[0]);printf("nh=%d\t",nh);
        if (nh==1) flag=1; else flag=0;
	z=mxGetPr(prhs[0]);
	t=mxGetPr(prhs[1]);
        nt=mxGetN(prhs[1]);
	h=mxGetPr(prhs[2]);
        nh=mxGetN(prhs[2]);
	p=mxGetPr(prhs[3]);
	np=mxGetN(prhs[3]);printf("np=%d\n",np);
	dt=t[1]-t[0];printf("dt=%f\n",dt);
	ww=mxGetPr(prhs[4]);
	/*output*/
        if (flag==1) plhs[0]=mxCreateDoubleMatrix(nt*np,1, mxREAL);
        else if (flag==0) plhs[0]=mxCreateDoubleMatrix(nt,np, mxREAL);

	x=mxGetPr(plhs[0]);
	/*////////////////////////////*/	

	rstack(x,t,h,p,z,nt,nh,np,dt,ww);
    	
	return;
}


void rstack( double *m, double *t, double *h, double *q, double *d, 
			int nt, int nh, int nq, double dt, double *ww) 
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
		if (fabs(ww[iq])>1e-2){
			for (itau=0;itau<nt;itau++){
				k=iq*nt+itau;
			      /*time=sqrt(pow(t[itau],2)+pow(h[ih],2)*q[iq]);*/
                                time=t[itau]+pow(h[ih],2)*q[iq];
				it=time/dt;
				j=ih*nt+(int) floor(it);
				if ((it!=nt)&&(j<ny)&&(k<nx)) 
					m[k]=m[k]+d[j];/**sincin(it-floor(it));*/
		
			}
		}
	  
	}
  }
  return;
}






