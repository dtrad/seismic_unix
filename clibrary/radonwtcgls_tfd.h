#ifndef RADONWTCGLS_TFD_H
#define RADONWTCGLS_TFD_H

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

void radonwtcgls0_tfd(float **data, float *h, int nh,float *t, int nt, float dt, float **model, float *q, int nq, float *vel, float fmax, int itercg, int iter_end, float step, float eps2, float eps1, float quantil, int norm, float factor, float smute, float nmofactor, int rtmethod, float depth);

void radonwtcgls_tfd(float *pos, int nh, float **data, float *t, int nt, float dt,
		     float **model, float *q, int nq, float dq, float eps1, float eps2,
		     float eps, float fmax, float **Wd, int itercg, int iter_end,
		     int norm,float step, int testadj, float quantil, int rtmethod, 
		     float depth);

void hrrti(float **d,float *pos,float dt,float **m,float *q, float fmax, 
	   int nt, int  nh, int nq, int rtmethod);

void radon_matrix(complex **l, float *h,float *q,int nh,int nq,float w);

void radon_moveout(float *h, float *g, int nh, int rtmethod, float depth);

void radon_param(float fmax,float *,int, float dx_av, float qmin, 
   float *qmaxt, float *qmax, float *dq, int nq, int rtmethod, float factor);

void interpovv (int nt, int ncdp, float *cdp, float **ovv, float cdpt, float *ovvt);

void plotgather(float **d, int nh, int nt, float dt);

void plotgather(float **d, int nh, int nt, float dt, const char *s);

void save_gather(float **d, int nh, float *h, int nt, float dt, char* name);

void save_gather(float *d, int nh, int nt, float dt, char* name);

void save_gather(float **d, int nh, int nt, float dt, char* name);

void save_gather(complex **d, int nh, float *h, int nt, float dt, int nfft, char* name);

int taper(float **data, int nt, int nh, int ntaper,int flag);

void nmo(float *d,float *m,float *t,float h,float *vel,int adj,int nt,float dt,
	 float smute);

void getvelocities(float dt, int nt, int ncdp, float *cdp, float **ovv);

void interval(float *pos, int nh, float *dx_max, float *dx_av);

void fftgo0(int sign,float **,complex **, int nh, int nt, float dt, int *nf0);

void fftgo1(int sign,float *,complex **, int nh, int nt, float dt, int nf0);

int fftback0(int sign,float **,complex  **, int nh, int nt, float dt, int nf0);

int fftback1(int sign,float *,complex  **, int nh, int nt, float dt, int nf0);

void deviations(float *m, int nx, float *d, int ny, int norm, float quantil1, float quantil2, 
		float *sigmam, float *sigmad);

void weights(float *m, int nx, int norm, float eps1, float *Wm, int iter);

#endif


