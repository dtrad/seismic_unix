#ifndef RADONSOLVER_H
#define RADONSOLVER_H

#define LOOKFAC	2	/* Look ahead factor for npfaro	  */
#define PFA_MAX	720720	

#ifndef TRACE
#define TRACE fprintf(stderr,"%s @ %u\n",__FILE__,__LINE__)
#endif

#define floatprint(expr) fprintf(stderr,#expr " = %f\n", expr)
#define expfprint(expr) fprintf(stderr,#expr " = %g\n", expr)
#define intprint(expr) fprintf(stderr,#expr " = %d\n", expr)

#include"su.h"
#include"segy.h"
#include"dan.h"
#include "inversion_par.h"



void radonline_stack0(float **data, complex **oper,float *trace, float *h, int nh,float *t, 
		      int nt, float dt, float *vel, float smute, float nmofactor,  float fmax, 
		      int nfft, float df);

void radonline_LS(float **data, complex ***RT, complex ***LT, float **model, float *trace, 
		  float *h, int nh,float *t, 
		  int nt, float *q, int nq, float dt, float *vel, float smute, float nmofactor,
		  float fmax, int nfft, float df, float *mutevector, float *stackvector);

void test2(segy tr, int nh, int nt);

void line_stack0(float **data, complex **oper,float *trace, float *h, int nh,float *t, int nt, float dt, float *vel, float smute, float nmofactor,  float fmax, int nfft, float df);


/*** stack_rt_operator ***/
void stack_rt_operator(complex **F, int maxfreq, float df, float *h, int nh, 
		       float *q, int nq, float *t, 
		       int nt, float eps, float parmute);

void stack_rt_LS_operators(complex ***LS, complex ***L, int maxfreq, float df, 
			   float *h, int nh, float *q, int nq, float *t, int nt, float eps, 
			   float parmute, float *mutevector, float *stackop);

void LSmatrix(complex **F, complex **FLST, int nh, int nk, float eps);

void mute_vector(float *mute,float *q, int nq, float parmute);

void stack_vector(float *stackop, int nh);

void inverse_matrix_multiply (int nrows1, complex **matrix1, int ncols2,
			      int nrows2, complex **matrix2, complex **out_matrix);

void testop(complex **F, int maxfreq, int nh);

/**************************/
void radon_matrix(complex **l, float *h,float *q,int nh,int nq,float w);

void radon_moveout(float *h, float *g, int nh, int rtmethod, float depth);

//void hrrti(float **d,float *pos,float dt,float **m,float *q, float fmax, 
//	   int nt, int  nh, int nq, int rtmethod, float depth);

int taper(float **data, int nt, int nh, int ntaper,int flag);

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

void doublemute(float **data, float **model, int nq, int nh, int nt, float *q, float *h, float *t, float parmute1, float parmute2, int mute1, int mute2, float fmax, int rtmethod, float depth, float t0mute);

void weights_window_inv(complex **m, int buffer, int nq, int freq, int norm, float sigmam, float *Wm, int iter);


#endif









