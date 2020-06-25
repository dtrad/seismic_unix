#include <math.h>
#include <stdio.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
/* Matlab */
#include "mex.h"

#define printf mexPrintf

void semblance( double *x, double *t, double *h, double *p, double *y, 
			int nt, int nh, int np, double dt) ;


void mexFunction( int nlhs, mxArray *plhs[],
				 int nrhs, const mxArray *prhs[])

{
	double *x, *z;
	int nt,nh,np,flag,nrow,ncol;
	double *t,*p,*h, dt;
	if (nrhs==0) {
		printf("[v]=semblance(u,t,h,p)\n") ;
		printf("Semblance analysis as a function of q=1/V^2\n") ;
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

	dt=t[1]-t[0];printf("dt=%f\n",dt);
	/*output*/
        if (flag==1) plhs[0]=mxCreateDoubleMatrix(nt*np,1, mxREAL);
        else if (flag==0) plhs[0]=mxCreateDoubleMatrix(nt,np, mxREAL);

	x=mxGetPr(plhs[0]);
	/*////////////////////////////*/	

	semblance(x,t,h,p,z,nt,nh,np,dt);
    	
	return;
}


void semblance( double *m, double *t, double *h, double *q, double *d, 
		int nt, int nh, int nq, double dt) 
{
  int i,k,ih,iq,itau;
  unsigned int j;
  double time,it,*den,*suma;
  int nx, ny;
  int nwind=3;
  nx=nt*nq;
  ny=nt*nh;
  printf("nx=%d,ny=%d,\n",nx,ny);
  suma=(double*) malloc((size_t) (nx+1)*sizeof(double));
  den=(double*) malloc((size_t) (nx+1)*sizeof(double));
  if (!den) printf("allocation failure in vector, semblance ");
  for (i=0;i<(nq*nt);i++) m[i]=0;
  
  for (ih=0;ih<nh;ih++){
    for (iq=0;iq<nq;iq++){
      for (itau=0;itau<nt;itau++){
	k=iq*nt+itau;
	time=sqrt(pow(t[itau],2)+pow(h[ih],2)*q[iq]);
	it=time/dt;
	j=ih*nt+(int) floor(it);
	if ((it!=nt)&&(j<ny)&&(k<nx)){ 
	  m[k]=m[k]+d[j];
	  den[k]=den[k]+d[j]*d[j]; /**sincin(it-floor(it));*/
	}
      }
      
    }
  }

  for (i=0;i<nx;i++) suma[i]=0;
  for (k=nwind;k<nx-nwind;k++)
    for (i=k-nwind;i<k+nwind;i++) suma[i]+=(m[i]*m[i]);
  for (k=0;k<nx;k++)  m[k]=suma[k];
  
  for (i=0;i<nx;i++) suma[i]=0;
  for (k=nwind;k<nx-nwind;k++) 
    for (i=k-nwind;i<k+nwind;i++) suma[i]+=(den[i]);
  for (k=0;k<nx;k++)  den[k]=suma[k];
  
  for (k=0;k<nx;k++) if (den[k]>1e-33) m[k]=(m[k])/(nh*den[k]); else m[k]==0;
  free(den);
  free(suma);
  return;
}









