#include <math.h>
#include <stdio.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
/* Matlab */
#include "mex.h"

#define printf mexPrintf

void rstack( double *x, double *t, double *h, double *p, double *y, 
			int nt, int nh, int np, double dt, double theta) ;
double sinc(double arg);

void mexFunction( int nlhs, mxArray *plhs[],
				 int nrhs, const mxArray *prhs[])

{
	double *x, *z;
	int nt,nh,np,flag,nrow,ncol;
	double *t,*p,*h, dt, theta;
	if (nrhs==0) {
		printf("[u]=radonop(v,t,h,p,theta)\n") ;
		printf("L hyperbolic radon operator with dip\n") ;	
		printf("v and u can be input and output as 1D vectors\n") ;
		printf("or as 2D vectors. v(nt,np) gives u(nt,nh) \n") ;
		printf("v(ntxnp,1) gives u(ntxnh,1) \n");
		printf("(time is fast dimension)\n") ;
		printf("to use in CG algorithms\n") ;
		printf("Theta is the dip angle \n") ;
		printf("Daniel Trad - UBC\n") ;
		return;
	}
	


	/*/////////////////////////////////
	  // Get input from MATLAB */

	/*nt=mxGetM(prhs[0]);printf("nt=%d\t",nt);*/

	np=mxGetN(prhs[0]);/*printf("np=%d\t",np);*/
	if (np==1) flag=1; else flag=0;
	x=mxGetPr(prhs[0]);

	t=mxGetPr(prhs[1]);
	nrow=mxGetM(prhs[1]);
	ncol=mxGetN(prhs[1]);
	nt=(nrow>ncol)? nrow : ncol; 
	printf("nt=%d\n",nt);
        
	h=mxGetPr(prhs[2]);
	nrow=mxGetM(prhs[2]);
	ncol=mxGetN(prhs[2]);
	nh=(nrow>ncol)? nrow : ncol;
	printf("nh=%d\n",nh);

	p=mxGetPr(prhs[3]);
	nrow=mxGetM(prhs[3]);
	ncol=mxGetN(prhs[3]);
	np=(nrow>ncol)? nrow : ncol;
	printf("np=%d\n",np);

        
        theta=mxGetScalar(prhs[4]);
	dt=t[1]-t[0];/*printf("dt=%f\n",dt);*/
	/*output*/
	if (flag==1) plhs[0]=mxCreateDoubleMatrix(nt*nh,1,mxREAL);
	else if (flag==0) plhs[0]=mxCreateDoubleMatrix(nt,nh, mxREAL);
	z=mxGetPr(plhs[0]);
	/*////////////////////////////*/	

	rstack(x,t,h,p,z,nt,nh,np,dt,theta);
    	
	return;
}


void rstack( double *m, double *t, double *h, double *q, double *d, 
			int nt, int nh, int nq, double dt, double theta) 
{
  int i,k,ih,iq,itau;
  unsigned int j;
  double time,it,hxh,qxhxh,qxhxhxs,sintheta,a1,a2;
  int nx, ny,itint;
  nx=nt*nq;
  ny=nt*nh;
  theta=theta/180.*acos(-1.);
  printf("theta in radians %e\n",theta);
  sintheta=sin(theta)*sin(theta);
  printf("sin(theta)^2  %e\n",sintheta);
  for (i=0;i<(nh*nt);i++) d[i]=0;

  for (ih=0;ih<nh;ih++){
    hxh=h[ih]*h[ih];
    for (iq=0;iq<nq;iq++){
      qxhxh=q[iq]*hxh;
      qxhxhxs=qxhxh*sintheta;
      for (itau=0;itau<nt;itau++){
	k=iq*nt+itau;
	if (fabs(m[k])>1e-2){
	  time=sqrt(pow(t[itau],2)+qxhxh-qxhxhxs);
	  it=time/dt;
	  itint=(int) floor(it);
	  j=ih*nt+itint;
	  a1=1-(it-itint);
	  a2=it-itint;
	  if ((it!=nt)&&(j<ny-1)&&(k<nx)) {
	    d[j]=d[j]+a1*m[k];
	    d[j+1]=d[j+1]+a2*m[k];
	  }
	}
      }
      
    }
  }
  return;
}



double sinc(double arg)
{
	double y;
	double pi = 3.1415926535;
	if (arg<1e-1) return(y=1.);
	y=sin(pi*arg)/(pi*arg);
	return(y);
}


