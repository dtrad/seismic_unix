#ifndef RADONHYBRIDINV_H
#define RADONHYBRIDINV_H

#define floatprint(expr) fprintf(stderr,#expr " = %f\n", expr)
#define expfprint(expr) fprintf(stderr,#expr " = %g\n", expr)
#define intprint(expr) fprintf(stderr,#expr " = %d\n", expr)

#include"su.h"
#include"dan.h"
#include "inversion_par.h"

void radonhybridinv(float **data, float *h, int nh,float *t, int nt, float dt, float **model, float *q, int nq, float *vel, int nq1, int nq2, float smute, float nmofactor, int rtmethod1, int rtmethod2, float depth1, float depth2, float fmax1, float fmax2, int filtering, int npoly, float *ffilter, float *amps);

void polygonalFilter(float *f, float *amps, int npoly,
				int nfft, float dt, float *filter);


void radoninv_2op(float **data, float *ph1, float *ph2, int nh, float *t, int nt, 
		  float **model, float *q, int nq, float *q1, int nq1, float *q2, int nq2, 
		  float fmax1, float fmax2, int filtering, int npoly, float *f, float *amps);

void radon_param_2op(float fmax1, float fmax2, float *h, int nh, float *q, float nq, 
		      float qmin1, float qmin2, float *q1, 
		      float *q2, int nq1, int nq2, float depth1, 
		      float depth2, int  rtmethod1, int rtmethod2, float factor1,
		      float factor2,float *pdq1, float *pdq2);

void radon_matrix(complex **l, float *g,float *q,int nh,int nq,float w);

void radon_matrix_irrq(complex **L, float *h, float *q,int nh,int nq,float w);

void radon_matrix_2op(complex **L, complex **L1, complex **L2, int nh, int nq1, int nq2);

void radon_moveout(float *h, float *g, int nh, int rtmethod, float depth);

void interpovv (int nt, int ncdp, float *cdp, float **ovv, float cdpt, float *ovvt);

int taper(float **data, int nt, int nh, int ntaper,int flag);

void nmo(float *d,float *m,float *t,float h,float *vel,int adj,int nt,float dt,
	 float smute);

void getvelocities(float dt, int nt, int ncdp, float *cdp, float **ovv);


#endif













