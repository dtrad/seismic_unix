#ifndef RADONTOEP_H
#define RADONTOEP_H

#define floatprint(expr) fprintf(stderr,#expr " = %f\n", expr)
#define expfprint(expr) fprintf(stderr,#expr " = %g\n", expr)
#define intprint(expr) fprintf(stderr,#expr " = %d\n", expr)

#define LOOKFAC	2	/* Look ahead factor for npfaro	  */
#define PFA_MAX	720720	

#ifndef TRACE
#define TRACE fprintf(stderr,"%s @ %u\n",__FILE__,__LINE__)
#endif

#include"su.h"
#include"dan.h"
#include "inversion_par.h"


void radontoepf(float **data, float *h, int nh,float *t, int nt, float dt, float **model, float *q, int nq, float *vel, float eps1, float factor, float smute, float nmofactor, int rtmethod, float depth, float fmax);


void toepradon(float *pos, float **d, float dt,float **m, float *q, float dq,
	       float eps1, float eps2, float eps, float fmax, int nt, int nh,
	       int nq, int rtmethod, float depth);

int ctoeplitz( int n, complex *r, complex *a, complex *b,
		 complex *f, complex *g );

void hrrti(float **d,float *pos,float dt,float **m,float *q, float fmax, 
	   int nt, int  nh, int nq, int rtmethod);

void radon_matrix(complex *R, complex **l,complex **lh,float *h,float *q,int nh,int nq,float w,
		  float *dh);

void radon_matrix(complex **l, float *h,float *q,int nh,int nq,float w);

void interpovv (int nt, int ncdp, float *cdp, float **ovv, float cdpt, float *ovvt);

void plotgather(float **d, int nh, int nt, float dt);

void plotgather(float **d, int nh, int nt, float dt, const char *s);

void save_gather(float **d, int nh, int nt, float dt, char* name);

int taper(float **data, int nt, int nh, int ntaper,int flag);

void nmo(float *d,float *m,float *t,float h,float *vel,int adj,int nt,float dt,
	 float smute);

void getvelocities(float dt, int nt, int ncdp, float *cdp, float **ovv);

void radon_param(float fmax,float *,int, float dx_av, float qmin, 
   float *qmaxt, float *qmax, float *dq, int nq, int rtmethod, float factor);

void interval(float *pos, int nh, float *dx_max, float *dx_av);

void fftgo0(int sign,float **,complex **, int nh, int nt, float dt, int *nf0);

int fftback0(int sign,float **,complex  **, int nh, int nt, float dt, int nf0);

#endif


