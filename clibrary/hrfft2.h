#ifndef HRFFT2_H
#define HRFFT2_H

#define LOOKFAC	2	/* Look ahead factor for npfaro	  */
#define PFA_MAX	720720	

#ifndef TRACE
#define TRACE fprintf(stderr,"%s @ %u\n",__FILE__,__LINE__)
#endif

#include"su.h"
#include"dan.h"

#define floatprint(expr) fprintf(stderr,#expr " = %f\n", expr)
#define expfprint(expr) fprintf(stderr,#expr " = %g\n", expr)
#define intprint(expr) fprintf(stderr,#expr " = %d\n", expr)


void dft2_interface(float **data, complex **model, float *h, int nh, float *t, int nt, float *k, 
		    int nk, float *vel, float eps1, float smute, float nmofactor, float fmax, 
		    float **data, float *h2, int nh2);

void dft2_toep(float **d, complex **m, float *t, float *pos, float *k,
	       int nt, int nh, int nk, float eps1, float fmax, int *pnf0);

void dft2_inv(float **d, complex **m2, float *t, float *pos, float *k, int nt, int nh,int nk, 
	     float fmax, int nf0);

void fft_matrix(complex *R, complex **F,complex **FH,float *h,float *k,int nh,int nk, float *dh);

void interpovv (int nt, int ncdp, float *cdp, float **ovv, float cdpt, float *ovvt);

void plotgather(float **d, int nh, int nt, float dt);

void plotgather(float **d, int nh, int nt, float dt, const char *s);

void save_gather(float **d, int nh, int nt, float dt, char* name);

int taper(float **data, int nt, int nh, int ntaper,int flag);

void nmo(float *d,float *m,float *t,float h,float *vel,int adj,int nt,float dt, float smute);

void getvelocities(float dt, int nt, int ncdp, float *cdp, float **ovv);

void fftgo0(int sign, float **d, complex  **m, int nh, int nt, float dt, int *nf0);

int read_ascii_file(const char *name,float *x);

#endif


