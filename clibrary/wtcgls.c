#include <math.h>
#include <stdio.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
/* Matlab */
#include "mex.h"
#include "nrutil.c"
#include "nrutil.h"
#include "nr.h"
#define printf mexPrintf
void materr(char error_text[]);
void radonpi( double *m, double *t, double *h, double *q, double *d, 
			int nt, int nh, int nq, double dt, double *ww);
void radonp( double *m, double *t, double *h, double *q, double *d, 
			int nt, int nh, int nq, double dt);

void wtcgls(double *t, double *qaxis, double *h,double *x,double *b,double *Qp,double tol,
	    int itercg,int nx, int ny, int nt, int nh, int nq, double step);
double dot(int n, double *a, double *b);

double sinc(double arg);

void mexFunction( int nlhs, mxArray *plhs[],
				 int nrhs, const mxArray *prhs[])

{
	double *x, *z, *ww;
	int nt,nh,np,flag;
	double *t,*p,*h, dt, step, tol;
        int itercg;
        int nx,ny;

  
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
        itercg=mxGetScalar(prhs[5]);
        tol=mxGetScalar(prhs[6]);
        step=mxGetScalar(prhs[7]);

	/*output*/
        if (flag==1) plhs[0]=mxCreateDoubleMatrix(nt*np,1, mxREAL);
        else if (flag==0) plhs[0]=mxCreateDoubleMatrix(nt,np, mxREAL);
        nx=np*nt;
        ny=nh*nt;
        printf("nx=%d,ny=%d\n",nx,ny);
	x=mxGetPr(plhs[0]);
	/*////////////////////////////*/	
        wtcgls(t,p,h,x,z,ww,tol,itercg,nx,ny,nt,nh,np,step);   	
	return;
}

void wtcgls(double *t, double *qaxis, double *h,double *x,double *b,double *Qp,double tol,
int itercg,int nx, int ny, int nt, int nh, int nq, double step)
  
{
  double normb,dq,dq2,nit,beta,betanum, betaden, alpha, alphanum, alphaden;
  int k,i,j,in,num;
  double *q,*q1,*s,*x1,*z,*z1,*r,*Az,*eta,*rho,*gcv;
  double dt=t[2]-t[1];
  printf("nx=%d,ny=%d\n",nx,ny);
  printf("dt=%f,nt=%d,nh=%d,nq=%d,step=%f,tol=%f,itercg=%d\n",dt,nt,nh,nq,step,tol,itercg);
  

  
  if ((q=dvector(0,nx))==NULL) 
    materr("cannot allocate memory for qcg\n");
  if ((q1=dvector(0,nx))==NULL)
    materr("cannot allocate memory for q1cg\n");
  if ((s=dvector(0,nx))==NULL)
    materr("cannot allocate memory for scg\n");
  if ((x1=dvector(0,nx))==NULL)
    materr("cannot allocate memory for x1cg\n");
  if ((z=dvector(0,nx))==NULL)
    materr("cannot allocate memory for zcg\n");
  if ((z1=dvector(0,nx))==NULL)
    materr("cannot allocate memory for z1cg\n");
  if ((r=dvector(0,ny))==NULL)
    materr("cannot allocate memory for rcg\n");
  if ((Az=dvector(0,ny))==NULL)
    materr("cannot allocate memory for Azcg\n");
  if ((eta=dvector(0,nx))==NULL)
    materr("cannot allocate memory for eta\n");
  if ((rho=dvector(0,nx))==NULL)
    materr("cannot allocate memory for rho\n");  
  if ((gcv=dvector(0,nx))==NULL)
     materr("cannot allocate memory for gcv\n");   
  /*for (i=0;i<nx;i++) Qp[i]=1;*/
    
    /*if (sqrt(fabs(x[i]))>1e-3) Qp[i]=1./sqrt(fabs(x[i]));
    else Qp[i]=1./1e-3;
  */
  for (i=0;i<nx;i++) x[i]=0.;
  normb=sqrt(dot(ny,b,b));
  /*printf("normb=%e\n",normb);*/
  for (i=0;i<ny;i++) r[i]=b[i];
  radonpi(s,t,h,qaxis,r,nt,nh,nq,dt,Qp);
  /*for (i=0;i<10;i++) printf("s[i]=%e\n",s[i]);*/

  nit=itercg;
  for(i=0;i<nx;i++){
    q1[i]=s[i]/Qp[i];
    q[i]=q1[i]/Qp[i];
  }
  for (i=0;i<nx;i++) z[i]=q[i];
  dq=dot(nx,s,q);
  for (i=0;i<nx;i++) z1[i]=q1[i];
  for(i=0;i<nx;i++) x1[i]=0.;       
  for (j=0;j<nit;j++){
    /*rstack(t,qaxis,h,z,Az,0);*/
    radonp(z,t,h,qaxis,Az,nt,nh,nq,dt);
    alphanum=dq;
    alphaden=dot(ny,Az,Az);
    printf("alphaden=%e\n",alphaden);
    if (alphaden < 0.) printf("alphaden=%e\n",alphaden);
    if (alphaden < tol ){ 
      printf("alphanum=%e,alphaden=%e,j=%d\n",
	      alphanum,alphaden,j);
      /*break;*/
    }
    alpha=alphanum/alphaden;
    alpha*=step;
    /*printf("j=%d,alpha=%e\n",j,alpha); */        
    /* Update model u and residuals*/

    for(i=0;i<nx;i++) x[i]=x[i]+alpha*z[i];
    for(i=0;i<ny;i++) r[i]=r[i]-alpha*Az[i];  
    
    /*resold=resid;*/
    /*rstack(t,qaxis,h,s,r,1);*/
    radonpi(s,t,h,qaxis,r,nt,nh,nq,dt,Qp);
    for(i=0;i<nx;i++){
      q1[i]=s[i]/Qp[i];
      q[i]=q1[i]/Qp[i];
    }
    dq2=dot(nx,s,q);
    beta=dq2/dq;
    /*printf("j=%d,beta=%e\n",j,beta);*/
    dq=dq2;
    for (i=0;i<nx;i++) z[i]=q[i]+beta*z[i];
    rho[j] = sqrt(dot(ny,r,r))/normb;
    printf("rho[%d]=%e\n",j,rho[j]);      
    for (i=0;i<nx;i++) {
      x1[i]=x1[i]+alpha*z1[i]; 
      z1[i]=q1[i]+beta*z1[i];
    }
    eta[j]=sqrt(dot(nx,x1,x1));
    if ((tol==0) && (j>2)){ /* GCV criteria */
       in = j-1;
       for (i=1;i<=in;i++){
       num=(ny-(i-1))*(ny-(i-1)); 
       gcv[i]=(rho[i]*rho[i])/num;
       }      
       if (gcv[j-2]<gcv[j-1]){ 
         printf("GCV Criteria was reached in iteration %d\n",j-1);
         nit = j-1;
         return;
       } 
       
       else if ((tol!=0) && (rho[j] < tol)){ 
        printf("Convergence have been acheived at iteration # %d\n",j);
        return;
       }
    }          
  }
  printf("j=%d\n",j);
  free_dvector(gcv,0,nx);
  free_dvector(rho,0,nx);
  free_dvector(eta,0,nx);
  free_dvector(Az,0,ny);
  free_dvector(r,0,ny);
  free_dvector(z1,0,nx);
  free_dvector(z,0,nx);
  free_dvector(x1,0,nx);
  free_dvector(s,0,nx);
  free_dvector(q,0,nx);
  free_dvector(q1,0,nx); 

  return;
}


double dot(int n, double *a, double *b)
/********************************************************************  
return the  dot product
*********************************************************************/
{
	int j;
	double sum=0.;
	for(j=0;j<n;j++) sum += a[j]*b[j];
	return(sum);
}


void radonpi( double *m, double *t, double *h, double *q, double *d, 
			int nt, int nh, int nq, double dt, double *ww) 
{
  int i,k,ih,iq,itau;
  unsigned int j;
  double time,it,pxhxh,hxh;
  int nx, ny;
  nx=nt*nq;
  ny=nt*nh;
    
  for (i=0;i<(nq*nt);i++) m[i]=0;

  for (ih=0;ih<nh;ih++){
	hxh=h[ih]*h[ih];
    for (iq=0;iq<nq;iq++){
		pxhxh=q[iq]*hxh;
		if (fabs(ww[iq])>1e-2){
			for (itau=0;itau<nt;itau++){
				k=iq*nt+itau;
				time=sqrt(pow(t[itau],2)+pxhxh);
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


void radonp( double *m, double *t, double *h, double *q, double *d, 
			int nt, int nh, int nq, double dt) 
{
  int i,k,ih,iq,itau;
  unsigned int j;
  double time,it,pxhxh,hxh;
  int nx, ny;
  nx=nt*nq;
  ny=nt*nh;
  
    
  for (i=0;i<(nh*nt);i++) d[i]=0;

  for (ih=0;ih<nh;ih++){
    hxh=h[ih]*h[ih];
    for (iq=0;iq<nq;iq++){
      pxhxh=q[iq]*hxh;
	  for (itau=0;itau<nt;itau++){
	    k=iq*nt+itau;
		if (fabs(m[k])>1e-10){
			time=sqrt(pow(t[itau],2)+pxhxh);
			it=time/dt;
			j=ih*nt+(int) floor(it);
			if ((it!=nt)&&(j<ny)&&(k<nx)) {
			d[j]=d[j]+m[k];/**sinc(it-floor(it));*/			
			}
		}
	  }
	  
    }
  }
  return;
}

void materr(char error_text[])
/* Numerical Recipes standard error handler */
{
	printf("Numerical Recipes run-time error...\n");
	printf("%s\n",error_text);
	printf("...now exiting to system...\n");
	exit(1);
}



