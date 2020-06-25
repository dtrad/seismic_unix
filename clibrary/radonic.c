#include <math.h>
#include <stdio.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
/* Matlab */
#include "mex.h"

#define printf mexPrintf

void rstack( double *x, double *t, double *h, double *p, double *y, 
			int nt, int nh, int np, double dt) ;

void mexFunction( int nlhs, mxArray *plhs[],
				 int nrhs, const mxArray *prhs[])

{
	double *x, *z, *ww;
	int nt,nh,np,flag,nrow,ncol,flag2;
	double *t,*p,*h, dt;
	if (nrhs==0) {
		printf("[v]=radonic(u,t,h,p)\n") ;
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
	nh=mxGetN(prhs[0]);/*printf("nh=%d\t",nh);*/
        if (nh==1) flag=1; else flag=0;
	z=mxGetPr(prhs[0]);

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


	
	dt=t[1]-t[0];/*printf("dt=%f\n",dt);*/
        flag2=0;
        if (mxIsEmpty(prhs[4])) flag2=1;
	else  ww=mxGetPr(prhs[4]);
        printf("No Weight matrix flag2=%d;\n",flag2);
        
	/*output*/
        if (flag==1) plhs[0]=mxCreateDoubleMatrix(nt*np+1,1, mxREAL);
        else if (flag==0) plhs[0]=mxCreateDoubleMatrix(nt+1,np+1, mxREAL);

	x=mxGetPr(plhs[0]);
	/*////////////////////////////*/	

        rstack(x,t,h,p,z,nt,nh,np,dt);
    	
	return;
}



void rstack( double *m, double *t, double *h, double *q, double *d, 
			int nt, int nh, int nq, double dt) 
{
  int i,k,ih,iq,itau;
  unsigned int j;
  double time,it,hxh,pxhxh,a1,a2,itint;
  int nx, ny,iqxnt,ihxnt;
  
  nx=nt*nq;
  ny=nt*nh;
    
  for (i=0;i<(nq*nt);i++) m[i]=0;

  for (ih=0;ih<nh;ih++){
    hxh=h[ih]*h[ih];
    ihxnt=ih*nt;
    for (iq=0;iq<nq;iq++){
      pxhxh=hxh*q[iq];
      iqxnt=iq*nt;
      for (itau=0;itau<nt;itau++){
	k=iqxnt+itau;
	time=sqrt(pow(t[itau],2)+pxhxh);       
    	it=time/dt;
	itint=(int) floor(it);
	j=ihxnt+itint;

	if ((it!=nt)&&(j<ny-1)&&(k<nx)) 
	  m[k]+=d[j]; 
	
      }
    }  
  }
  return;
}













