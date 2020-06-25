#ifndef HRFFT2_H
#define HRFFT2_H

#define LOOKFAC	2	/* Look ahead factor for npfaro	  */
#define PFA_MAX	720720	

#ifndef TRACE
#define TRACE fprintf(stderr,"%s @ %u\n",__FILE__,__LINE__)
#endif

#include"su.h"
#include"dan.h"
#include "inversion_par.h"

#define floatprint(expr) fprintf(stderr,#expr " = %f\n", expr)
#define expfprint(expr) fprintf(stderr,#expr " = %g\n", expr)
#define intprint(expr) fprintf(stderr,#expr " = %d\n", expr)

void dft2_interface(float **data, complex **model, float *h, int nh, float *t, int nt, float *k, 
		    int nk, float *vel, float smute, float nmofactor, float fmax, 
		    float **data2, float *h2, int nh2, float **Wd, inv_par inv, int method, float& df);

void dft2_toep(float **d, complex **m, float *t, float *pos, float *k,
	       int nt, int nh, int nk, float eps1, float fmax, int *pnf0);

void dft2_wtcgls(float **d, complex **m2, float *t, float *pos, float *k,
	    int nt, int nh, int nk, float **Wd, inv_par inv, float fmax, int *pnf0);

void dft2_inv(float **d, complex **m2, float *t, float *pos, float *k,
	     int nt, int nh, int nk, float fmax, int nf0);

void dft_matrix(complex *R, complex **F,complex **FH,float *h,float *k,int nh,int nk, float *dh);

void dft_matrix(complex **F,float *h,float *k,int nh,int nk);

void modelweight_inv(complex *m, int nx, int norm, float eps1, float *Wm);

void wtcgls(complex *b,complex **L, complex *x,float *Wm,
	  float *Wd,int nh,int nq,complex *q,complex *q1,complex *s,
	  complex *x1,complex *z,complex *z1,complex *r,complex *Az,
	  float *eta,float *rho, float *gcv, float tol, float step, int itercg);

int ctoeplitz( int n, complex *r, complex *a, complex *b, complex *f, complex *g );

void interpovv (int nt, int ncdp, float *cdp, float **ovv, float cdpt, float *ovvt);

void plotgather(float **d, int nh, int nt, float dt);

void plotgather(float **d, int nh, int nt, float dt, const char *s);

void save_gather(float **d, int nh, int nt, float dt, char* name);

int taper(float **data, int nt, int nh, int ntaper,int flag);

void nmo(float *d,float *m,float *t,float h,float *vel,int adj,int nt,float dt, float smute);

void getvelocities(float dt, int nt, int ncdp, float *cdp, float **ovv);

void fftgo0(int sign, float **d, complex  **m, int nh, int nt, float dt, int *nf0);

int fftback0(int sign,float **,complex  **, int nh, int nt, float dt, int nf0);

int read_ascii_file(const char *name,float *x);

#endif


