void nlcg1_interface(complex *d,complex **L,complex *m, float nh, float nq, float ftol, int itercg, )
{
  float *dr;
  float *mr;
  float **LR;
  float (*func)(float []);
  void (*dfunc)(float [], float []);

  func=costfunc;
  dfunc=gradient;

  // Transform the complex system of equation in a real one 
  dr=ealloc1float(2*nh);
  mr=ealloc1float(2*nq);
  LR=ealloc2float(2*nq,2*nh);

  complex2real(d[freq],L,dr,LR,nh,nq);
  
  // test
  if (1) matmul0 (dr,LR,mr,1,2*nh,2*nq);
     
  nlcg(d2r,LR,m2r,2*nq,2*nh,iter_end,itercg);  
  frprmn(mr,2*nq,ftol, itercg,&fret,func,dfunc);
  matmul0 (d2r,LR,m2r,0,2*nh,2*nq);
  real2complex(d2[freq],d2r,m2[freq],m2r,nh,nq);  

  free2float(LR);
  free1float(mr);
  free1float(dr);

}



void dfunc1(float *x, float *xc);
float func1(float *x);
void linmin(float p[], float xi[], int n, float *fret,float (*func)(float []));

#include "nrutil.h"
#include "Complex.h"
#include "su.h"
#include "dan.h"




void frprmn(float p[], int n, float ftol, int *iter, float *fret,
	float (*func)(float []), void (*dfunc)(float [], float []))
{
	void linmin(float p[], float xi[], int n, float *fret,
		float (*func)(float []));
	void dlinmin(float p[], float xi[], int n, float *fret, 
		    float (*func)(float []),void (*dfunc)(float [], float []));
	int j,its;
	float gg,gam,fp,dgg;
	float *g,*h,*xi;

	g=vector(1,n);
	h=vector(1,n);
	xi=vector(1,n);
 	fp=(*func)(p);
	(*dfunc)(p,xi);
	for (j=1;j<=n;j++) {
		g[j] = -xi[j];
		xi[j]=h[j]=g[j];
	}
 
	for (its=1;its<=ITMAX;its++) {
		*iter=its;
		dlinmin(p,xi,n,fret,func,dfunc);
		if (2.0*fabs(*fret-fp) <= ftol*(fabs(*fret)+fabs(fp)+EPS)) {
			FREEALL
			return;
		}
		fp=(*func)(p);
		(*dfunc)(p,xi);
		dgg=gg=0.0;
		for (j=1;j<=n;j++) {
			gg += g[j]*g[j];
			dgg += (xi[j]+g[j])*xi[j];
		}
		if (gg == 0.0) {
			FREEALL
			return;
		}
		gam=dgg/gg;
		for (j=1;j<=n;j++) {
			g[j] = -xi[j];
			xi[j]=h[j]=g[j]+gam*h[j];
		}
	}
	nrerror("Too many iterations in frprmn");
}
#undef ITMAX
#undef EPS
#undef FREEALL
#undef NRANSI
/* (C) Copr. 1986-92 Numerical Recipes Software . */
#define NRANSI
#include "nrutil.h"
#define TOL 2.0e-4

int ncom;
float *pcom,*xicom,(*nrfunc)(float []);
void (*nrdfun)(float [], float []);

void dlinmin(float p[], float xi[], int n, float *fret, float (*func)(float []),void (*dfunc)(float [], float []))
{
	float dbrent(float ax, float bx, float cx,
		float (*f)(float), float (*df)(float), float tol, float *xmin);
	float f1dim(float x);
	float df1dim(float x);
	void mnbrak(float *ax, float *bx, float *cx, float *fa, float *fb,
		float *fc, float (*func)(float));
	int j;
	float xx,xmin,fx,fb,fa,bx,ax;

	ncom=n;
	pcom=vector(1,n);
	xicom=vector(1,n);
	nrfunc=func;
	nrdfun=dfunc;
	for (j=1;j<=n;j++) {
		pcom[j]=p[j];
		xicom[j]=xi[j];
	}
	ax=0.0;
	xx=1.0;
	mnbrak(&ax,&xx,&bx,&fa,&fx,&fb,f1dim);
	*fret=dbrent(ax,xx,bx,f1dim,df1dim,TOL,&xmin);
	for (j=1;j<=n;j++) {
		xi[j] *= xmin;
		p[j] += xi[j];
	}
	free_vector(xicom,1,n);
	free_vector(pcom,1,n);
}
#undef TOL
#undef NRANSI
/* (C) Copr. 1986-92 Numerical Recipes Software . */
#define NRANSI
#include "nrutil.h"

extern int ncom;
extern float *pcom,*xicom,(*nrfunc)(float []);
extern void (*nrdfun)(float [], float []);

float df1dim(float x)
{
	int j;
	float df1=0.0;
	float *xt,*df;

	xt=vector(1,ncom);
	df=vector(1,ncom);
	for (j=1;j<=ncom;j++) xt[j]=pcom[j]+x*xicom[j];
	(*nrdfun)(xt,df);
	for (j=1;j<=ncom;j++) df1 += df[j]*xicom[j];
	free_vector(df,1,ncom);
	free_vector(xt,1,ncom);
	return df1;
}
#undef NRANSI
/* (C) Copr. 1986-92 Numerical Recipes Software . */
#include <math.h>
#define NRANSI
#include "nrutil.h"
#define ITMAX 100
#define ZEPS 1.0e-10
#define MOV3(a,b,c, d,e,f) (a)=(d);(b)=(e);(c)=(f);

float dbrent(float ax, float bx, float cx, float (*f)(float),
	float (*df)(float), float tol, float *xmin)
{
	int iter,ok1,ok2;
	float a,b,d,d1,d2,du,dv,dw,dx,e=0.0;
	float fu,fv,fw,fx,olde,tol1,tol2,u,u1,u2,v,w,x,xm;

	a=(ax < cx ? ax : cx);
	b=(ax > cx ? ax : cx);
	x=w=v=bx;
	fw=fv=fx=(*f)(x);
	dw=dv=dx=(*df)(x);
	for (iter=1;iter<=ITMAX;iter++) {
		xm=0.5*(a+b);
		tol1=tol*fabs(x)+ZEPS;
		tol2=2.0*tol1;
		if (fabs(x-xm) <= (tol2-0.5*(b-a))) {
			*xmin=x;
			return fx;
		}
		if (fabs(e) > tol1) {
			d1=2.0*(b-a);
			d2=d1;
			if (dw != dx) d1=(w-x)*dx/(dx-dw);
			if (dv != dx) d2=(v-x)*dx/(dx-dv);
			u1=x+d1;
			u2=x+d2;
			ok1 = (a-u1)*(u1-b) > 0.0 && dx*d1 <= 0.0;
			ok2 = (a-u2)*(u2-b) > 0.0 && dx*d2 <= 0.0;
			olde=e;
			e=d;
			if (ok1 || ok2) {
				if (ok1 && ok2)
					d=(fabs(d1) < fabs(d2) ? d1 : d2);
				else if (ok1)
					d=d1;
				else
					d=d2;
				if (fabs(d) <= fabs(0.5*olde)) {
					u=x+d;
					if (u-a < tol2 || b-u < tol2)
						d=SIGN(tol1,xm-x);
				} else {
					d=0.5*(e=(dx >= 0.0 ? a-x : b-x));
				}
			} else {
				d=0.5*(e=(dx >= 0.0 ? a-x : b-x));
			}
		} else {
			d=0.5*(e=(dx >= 0.0 ? a-x : b-x));
		}
		if (fabs(d) >= tol1) {
			u=x+d;
			fu=(*f)(u);
		} else {
			u=x+SIGN(tol1,d);
			fu=(*f)(u);
			if (fu > fx) {
				*xmin=x;
				return fx;
			}
		}
		du=(*df)(u);
		if (fu <= fx) {
			if (u >= x) a=x; else b=x;
			MOV3(v,fv,dv, w,fw,dw)
			MOV3(w,fw,dw, x,fx,dx)
			MOV3(x,fx,dx, u,fu,du)
		} else {
			if (u < x) a=u; else b=u;
			if (fu <= fw || w == x) {
				MOV3(v,fv,dv, w,fw,dw)
				MOV3(w,fw,dw, u,fu,du)
			} else if (fu < fv || v == x || v == w) {
				MOV3(v,fv,dv, u,fu,du)
			}
		}
	}
	nrerror("Too many iterations in routine dbrent");
	return 0.0;
}
#undef ITMAX
#undef ZEPS
#undef MOV3
#undef NRANSI
/* (C) Copr. 1986-92 Numerical Recipes Software . */
#include <math.h>
#define NRANSI
#include "nrutil.h"
#define GOLD 1.618034
#define GLIMIT 100.0
#define TINY 1.0e-20
#define SHFT(a,b,c,d) (a)=(b);(b)=(c);(c)=(d);

void mnbrak(float *ax, float *bx, float *cx, float *fa, float *fb, float *fc,
	float (*func)(float))
{
	float ulim,u,r,q,fu,dum;

	*fa=(*func)(*ax);
	*fb=(*func)(*bx);
	if (*fb > *fa) {
		SHFT(dum,*ax,*bx,dum)
		SHFT(dum,*fb,*fa,dum)
	}
	*cx=(*bx)+GOLD*(*bx-*ax);
	*fc=(*func)(*cx);
	while (*fb > *fc) {
		r=(*bx-*ax)*(*fb-*fc);
		q=(*bx-*cx)*(*fb-*fa);
		u=(*bx)-((*bx-*cx)*q-(*bx-*ax)*r)/
			(2.0*SIGN(FMAX(fabs(q-r),TINY),q-r));
		ulim=(*bx)+GLIMIT*(*cx-*bx);
		if ((*bx-u)*(u-*cx) > 0.0) {
			fu=(*func)(u);
			if (fu < *fc) {
				*ax=(*bx);
				*bx=u;
				*fa=(*fb);
				*fb=fu;
				return;
			} else if (fu > *fb) {
				*cx=u;
				*fc=fu;
				return;
			}
			u=(*cx)+GOLD*(*cx-*bx);
			fu=(*func)(u);
		} else if ((*cx-u)*(u-ulim) > 0.0) {
			fu=(*func)(u);
			if (fu < *fc) {
				SHFT(*bx,*cx,u,*cx+GOLD*(*cx-*bx))
				SHFT(*fb,*fc,fu,(*func)(u))
			}
		} else if ((u-ulim)*(ulim-*cx) >= 0.0) {
			u=ulim;
			fu=(*func)(u);
		} else {
			u=(*cx)+GOLD*(*cx-*bx);
			fu=(*func)(u);
		}
		SHFT(*ax,*bx,*cx,u)
		SHFT(*fa,*fb,*fc,fu)
	}
}
#undef GOLD
#undef GLIMIT
#undef TINY
#undef SHFT
#undef NRANSI
/* (C) Copr. 1986-92 Numerical Recipes Software . */
#define NRANSI
#include "nrutil.h"

extern int ncom;
extern float *pcom,*xicom,(*nrfunc)(float []);

float f1dim(float x)
{
  int j;
  float f,*xt;
  xt=vector(1,ncom);
  for (j=1;j<=ncom;j++) xt[j]=pcom[j]+x*xicom[j];
  f=(*nrfunc)(xt);
  free_vector(xt,1,ncom);
  return f;
}
#undef NRANSI
/* (C) Copr. 1986-92 Numerical Recipes Software . */



//#include "clibrary.h"
#include <math.h>
#define NRANSI
#include "nrutil.h"
float func1(float *x)
{
  extern int nh, nq;
  float Jt,Jd,Jm;
  int i,j,nnh,nnq;
  extern float **LR, *CD, *CM, *RHS, *DC,*RES, lambda, sigma;  
  
  nnh=2*nh;
  nnq=2*nq;

  Atimesx(nnh,nnq,DC,LR,x);  
  for (i=1;i<=nnh;i++) RES[i]=-RHS[i]+DC[i];
  for (Jd=0.,i=1;i<=nnh;i++) Jd+=RES[i]*RES[i];
  //for (Jd=0.,i=1;i<=nnh;i++) Jd+=abs(RES[i]);
  for (Jm=0.,i=1;i<=nnq;i++) Jm+=ABS(x[i]);
  //for (Jm=0.,i=1;i<=nnq;i++) Jm+=log(x[i]*x[i]+sigma);
  Jt=Jd+lambda*Jm;
  //fprintf(stderr,"Jt=%e,Jd=%e,Jm=%e\n",Jt,Jd,lambda*Jm);
  //for (i=1;i<nnq;i++) fprintf(stderr,"RES[%d]=%e\n",i,RES[i]);

  return(Jt);
}


#include <math.h>
#define NRANSI
#include "nrutil.h"
//#include "clibrary.h"


void dfunc1(float *x, float *xc)
{
  extern int  nq, nh;
  extern float **LR, **LRH, *CM, *RHS, *DC, *XXC, *RES, lambda, sigma;  
  int i, nnq, nnh;
  
  
  nnh=2*nh;
  nnq=2*nq;
  //fprintf(stderr,"###############Inside dfunc1\n");
  //Atimesx(nnh,nnq,DC,LR,x);
//for (i=1;i<=nnh;i++) dc[i]/=Cd[i];
  //for (i=1;i<=nnh;i++) RES[i]=RHS[i]-DC[i];
  Atimesx(nnq,nnh,xc,LRH,RES);
  //for (i=1;i<=nnq;i++) {fprintf(stderr,"xc[%d]=%e\n",i,xc[i]);
  for (i=1;i<=nnq;i++) XXC[i]=lambda*SGN(x[i]);
  //for (i=1;i<=nnq;i++) XXC[i]=lambda*2*x[i]/(x[i]*x[i]+sigma);
  for (i=1;i<=nnq;i++) xc[i]+=XXC[i];
  //fprintf(stderr,"*********sigma=%f***\n",sigma);
  return;
}







































































