#ifndef FFT2_H
#define FFT2_H

#ifndef TRACE
#define TRACE fprintf(stderr,"%s @ %u\n",__FILE__,__LINE__)
#endif

#include"su.h"
#include"dan.h"

#define floatprint(expr) fprintf(stderr,#expr " = %f\n", expr)
#define expfprint(expr) fprintf(stderr,#expr " = %g\n", expr)
#define intprint(expr) fprintf(stderr,#expr " = %d\n", expr)

void fft2_0(float **data, complex **model, float *h, int nh, float *t, int nt, 
	    float *k, int nk, float *w, int nw, float *vel, float eps1, 
	    float smute, float nmofactor, float fmax, float **data2, float *h2, int nh2);

void fft2(float **dataxt, complex **datafk, float *t, float *x, int nt, int nx, 
	  float *w, float *k, int *pnw, int *pnk, float **data2xt);

void interpovv (int nt, int ncdp, float *cdp, float **ovv, float cdpt, float *ovvt);

void plotgather(float **d, int nh, int nt, float dt);

void plotgather(float **d, int nh, int nt, float dt, const char *s);

void save_gather(float **d, int nh, int nt, float dt, char* name);

int taper(float **data, int nt, int nh, int ntaper,int flag);

void nmo(float *d,float *m,float *t,float h,float *vel,int adj,int nt,float dt, float smute);

void getvelocities(float dt, int nt, int ncdp, float *cdp, float **ovv);

int read_ascii_file(const char *name,float *x);

#endif


