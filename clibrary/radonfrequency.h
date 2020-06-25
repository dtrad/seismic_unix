#ifndef RADONPARABOLIC_H
#define RADONPARABOLIC_H

#define LOOKFAC	2	/* Look ahead factor for npfaro	  */
#define PFA_MAX	720720	

#ifndef TRACE
#define TRACE fprintf(stderr,"%s @ %u\n",__FILE__,__LINE__)
#endif

#define floatprint(expr) fprintf(stderr,#expr " = %f\n", expr)
#define expfprint(expr) fprintf(stderr,#expr " = %g\n", expr)
#define intprint(expr) fprintf(stderr,#expr " = %d\n", expr)

#include"su.h"
#include"dan.h"
#include "inversion_par.h"


void radonsolver0_subtract(float **data, float *h, int nh,float *t, int nt, float dt, float **model, float *q, int nq, float *vel, int itercg, int iter_end, float step, float eps2, float eps1, float quantil, int norm, float factor, float smute, float nmofactor, int rtmethod, float depth, float fmax, char *solver, float **M, int mute);


void radonsolver(float *pos, int nh, float **data, float *t, int nt, float dt,float **model, float *q, int nq, float dq, float eps1, float eps2, float eps, float fmax, float *Wd, int itercg, int iter_end,int norm,float step, int testadj, int rtmethod, float depth, char *solver);

void radon_matrix(complex **l, float *h,float *q,int nh,int nq,float w);

void radon_moveout(float *h, float *g, int nh, int rtmethod, float depth);

void hrrti(float **d,float *pos,float dt,float **m,float *q, float fmax, 
	   int nt, int  nh, int nq, int rtmethod, float depth);

int taper(float **data, int nt, int nh, int ntaper,int flag);

void adaptive_subtract(float **data1, float **data2, int nh, int nt);

void muting(float **model, int nq, int nt, float *q, float *t, float parmute, float t0, int side);


int taper(float *data, int nh, int ntaper);

void radon_param(float fmax,float *,int, float dx_av, float qmin, 
   float *qmaxt, float *qmax, float *dq, int nq, int rtmethod, float factor);

void interval(float *pos, int nh, float *dx_max, float *dx_av);

void interpovv (int nt, int ncdp, float *cdp, float **ovv, float cdpt, float *ovvt);

void getvelocities(float dt, int nt, int ncdp, float *cdp, float **ovv);

void nmo(float *d,float *m,float *t,float h,float *vel,int adj,int nt,float dt,
	 float smute);

void weights_inv(complex *m, int nx, int norm, float sigmam, float *Wm, int iter);

void deviations(complex *m, int nx, complex *d, int ny, int norm, float quantil1, 
		float quantil2, float *sigmam, float *sigmad);

void dataweigths(float *pos, int nh, float *Wd, int add);

void radon_matrix_cgfft(complex **L,complex *RC, int nh, int nq, int nf2, float *Cd);

void  circ_mult(int n,complex *RC,complex *f,complex *g,int nf, complex *fc, complex *gc);

float wtcgls(complex *b,complex **L, complex *x,float *Wm,
	     float *Wd,int nh,int nq, float tol, float step, int itercg);

float radon_cgfft(complex *d2, complex **L, complex *RC, complex *m2, int nh, int nq, int nf2, float *Wm, float *Wd, float eps, int itercg, float step);

float radon_toeplitz(complex *d, complex **L, complex *m, float sigmad, int nh, int nq);

float radon_cgtoep(complex *d, complex **L, complex *m, float sigmad, int nh, int nq);

void weights_window_inv(complex **m, int buffer, int nq, int freq, int norm, float sigmam, float *Wm, int iter);

int  stopcriteria(float *gcv, float *rho, float normdata, int nh, int nq, int iter, int gcvbool);

/************************************************************************/

void AtimesB(float **A, float **B, int nh, int nt);

typedef struct {	/* Parameters for inversion  */
  float tmin; 
  float tmax;
  int ihmin;
  int ihmax;
  float slope;
  float threshold;
} mutemask_par;

void mutemask2(float **M, int nh, int nt, float dt, mutemask_par mutepar);

void write_curve_mute(mutemask_par par);

/************************************************************************/

#endif









