#ifndef RADONLINETFD_H
#define RADONLINETFD_H

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

/* suradoncgfft_tfd uses this */
void radoncgfft0_tfd(float **data, float *h, int nh,float *t, int nt, float dt, 
		     float **model, float *q, int nq, float *vel, float fmax,  
                     float quantil, float factor, float smute, float nmofactor, 
                     int rtmethod, float depth, inv_par inv);
//////////////////////////////////////////////////////////////////////////////////////
/* suradoncgfft_tfd2 uses this */
void radoncgfft0_tfd_mute_conv(float **data, float *h, int nh,float *t, int nt, float dt, 
			       float **model, float *q, int nq, float *vel, float fmax,  
			       float quantil, float factor, float smute, float nmofactor, 
			       int rtmethod, float depth, inv_par inv, float parmute1, 
			       float parmute2, int mute1, int mute2, float t0, 
			       float *wavelet, int nw, int plot);
//////////////////////////////////////////////////////////////////////////////////////

void radoncgfft_tfd_conv(float *pos, int nh, float **data, float *t, int nt, float dt,
			 float **model, float *q, int nq, float dq,  float fmax, float **Wd, 
			 int testadj, float quantil, int rtmethod,  float depth, inv_par inv, 
			 float *wavelet, int nw);

float cgfft(float *b,  complex ***L, float *x, complex **d2, complex **RC, complex **m2, 
	    int nt, int nh, int nq, int nf2, float dt,  int nfft, float fmax, 
	    float *Wm, float *Wd, inv_par inv, float *wavelet, int nw);

void radonf(float *model, float *data, complex ***L,int nt, int nh, int nq, 
	    float dt,complex **m2, complex **d2, int nfft, float fmax, int adj);

void Atimesx_by_fft_loop(float *y, complex **A, float *x, complex **y2, complex **x2,
			 int nt, int n, float dt, int nfft, float fmax, int nf2);

void  Atimesx_by_fft(complex *y,complex *A,complex *x, int n, int nf, complex *yc, complex *xc);

void hrrti(float **d,float *pos,float dt,float **m,float *q, float fmax, 
	   int nt, int nh, int nq, int rtmethod, float depth);

/* radoncgfft_tfd_mute uses doublemute  */
void doublemute(float **data, float **model, int nq, int nh, int nt, float *q, float *h, float *t, float parmute1, float parmute2, int mute1, int mute2, float fmax, int rtmethod, float depth, float t0);
/* radoncgfft_tfd_mute_conv uses doublemute with comvolution */
void doublemute(float **data, float **model, int nq, int nh, int nt, float *q, float *h, float *t, float parmute1, float parmute2, int mute1, int mute2, float fmax, int rtmethod, float depth, float t0, float *wavelet, int nw, float *vel, float nmofactor, float smute, int plot);


//////////////////////////////////////////////////////////////////////////////////////

void radon_matrix(complex **l, float *h,float *q,int nh,int nq,float w);

void radon_matrix_cgfft(complex **L, complex *RC, int nh, int nq, int nf2, float *Cd);

void radon_moveout(float *h, float *g, int nh, int rtmethod, float depth);

void radon_param(float fmax,float *,int, float dx_av, float qmin, 
   float *qmaxt, float *qmax, float *dq, int nq, int rtmethod, float factor);

void interpovv (int nt, int ncdp, float *cdp, float **ovv, float cdpt, float *ovvt);

int taper(float **data, int nt, int nh, int ntaper,int flag);

void nmo(float *d,float *m,float *t,float h,float *vel,int adj,int nt,float dt,
	 float smute);

void getvelocities(float dt, int nt, int ncdp, float *cdp, float **ovv);

void interval(float *pos, int nh, float *dx_max, float *dx_av);

void deviations(float *m, int nx, float *d, int ny, int norm, float quantil1, float quantil2, 
		float *sigmam, float *sigmad);

void weights(float *m, int nx, int norm, float eps1, float *Wm, int iter);

void weights_cgfft(float *m, int nx, int norm, float sigmam, float *Wm, int iter);

float testadjop(void (*oper) (float *model, float *data, complex ***L,int nt, int nh, int nq, 
			      float dt,complex **m2, complex **d2, int nfft, float fmax, int adj)
		,complex ***L,int nt, int nh, int nq, float dt, complex **m2, complex **d2, 
		int nfft, float fmax, int adj);

int taper(float **data, int nt, int nh, int ntaper,int flag);

int  stopcriteria(float *gcv, float *rho, float normdata, int nh, int nq, int iter, int gcvbool);

//////////////////////////////////////////////////////////////////////////////////////
/***  Routines for convolution and correlation ***/
void contruc_2(int conj,int add, int nx, float *xx, int nb, float *bb,float *yy);

void conv_corr(float *din, float *dout, int nt, int nh, float *wavelet, int nw, int adj);

int get_wavelet(float *wavelet,const char *name,int nw, int type, float dt, float fpeak);

#endif




