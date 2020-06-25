#ifndef RADONLOGBAR_H
#define RADONLOGBAR_H

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

/* Encapsulae all parameters for the log barrier method */
typedef struct{

  int itercg;	/* Number maximum of iterations in CG*/
  int iter_end;	/* Number maximum of External Iterations */
  float PI_t;    /* Primal infeasibility tolerance */
  float DI_t;     /* Dual infeasibility tolerance */
  float GI_t;   /* GAP infeasibility */
  float delta; /* Perturbation hyperparameter */ 
  float gamma; /* Residual hyperparameter */

} logb_param;



void radonlogbar0(float **data, float *h, int nh,float *t, int nt, float dt, float **model, float *q, int nq, float *vel, float factor, float smute, float nmofactor, int rtmethod, float depth, float fmax, char *solver, logb_param par);

void radonlogbar(float *pos, int nh, float **data, float *t, int nt, float dt,float **model, float *q, int nq, float dq, float fmax, float *Wd,  int rtmethod, float depth, char *solver, logb_param par);

int logbarrier_interface(complex *d,complex **L,complex *m, int nh, int nq, logb_param par);

int logbarrier_method(float *d, float **L, float *m, int nq, int nh, logb_param par);

int PD_lp_bm(float *b, float **A, float *x, int nx, int ny, logb_param par);  

int  CG_solver(float **A,float *D,int np,int nd,
	       float *x, float *b, float delta, int iter_end);

void matmul0(float *y, float **A, float *x, int adj, int ny, int nx);

float min3(float a, float b, float c);

void complex2real(complex *b, complex **A, float *br, float **AR, int nb, int nx);

void complex2real(complex *b, complex **A, complex *m, float *br, float **AR, float *mr, 
		  int nb, int nx);

void real2complex(complex *b, float *br, complex *x, float *xr, int nb, int nx);

void radon_matrix(complex **l, float *h,float *q,int nh,int nq,float w);

void radon_moveout(float *h, float *g, int nh, int rtmethod, float depth);

void hrrti(float **d,float *pos,float dt,float **m,float *q, float fmax, 
	   int nt, int  nh, int nq, int rtmethod, float depth);

int taper(float **data, int nt, int nh, int ntaper,int flag);

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



#endif









