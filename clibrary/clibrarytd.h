#ifndef CLIBRARYTD_H
#define CLIBRARYTD_H
//#include "inversion_par.h"

#define floatprint(expr) fprintf(stderr,#expr " = %f\n", expr)
#define expfprint(expr) fprintf(stderr,#expr " = %g\n", expr)
#define intprint(expr) fprintf(stderr,#expr " = %d\n", expr)
// Hyperbolic Radon Transform in time domain
void radtd(float *t, float *q, float *h, float *x,float *y,float eps, 
	   float qmin, float qmax, float fmax);
void radtd2(float *t, float *q, float *h, float *x,float *y,float eps, 
	   float qmin, float qmax, float fmax, float *vel);
void radtd2(float *t, float *q, float *h, float *x,float *y,float eps, 
	    float qmin, float qmax, float fmax, float *vel, float dperv, 
	    float pervmin);
void radint(float *t, float *q, float *h, float *h2, float *m,float *d,float *dint,float eps,float qmin,float qmax, float fmax, float *vel, float dperv, float pervmin);
void radint(float *t, float *q, float *h, float *h2, float *m,float *d,
	    float *dint,float eps,float qmin,float qmax, float fmax, 
	    float *vel, float dperv, float pervmin,
	    float alum, float blum, float clum, int norm);

float dot(int n, float *a, float *b);
void lsqr(float *t,float *qaxis,float *h,float *x,float *d,float  tol,
	  int reorth);
void semblance( float *m, float *t, float *h, float *q, float *d, 
			int nt, int nh, int nq, float dt);
void radonopi( float *m, float *t, float *h, float *q, float *d); 
void radonopi( float *m, float *t, float *h, float *q, float *d, float *ww);
void radonopi( float *m, float *t, float *h, float *q, float *d, float *ww, float theta); 
void radonop(float *m, float *t, float *h, float *q, float *d);
void radonop(float *m, float *t, float *h, float *q, float *d, float theta);
void wtcglstd(float *t, float *qaxis, float *h,float *x,float *b,float *Qp, float tol, int rtmethod);
void wtcglstd(float *t, float *qaxis, float *h,float *x,float *b,float *Qp, float tol, float theta);
void wtcgls(void (*oper) (float *,float *,float *,float *,float *,float **,int ,int ,int,int), int nt, int nh, int np, float *t, float *h, float *p, float *x,float *b,float *Wm, float *Wd, float **vel,float tol, float step, float eps1, float eps2, int itercg);

// Velocity Analysis
void velanrt(float *t, float *q, float *h, float *x,float *y,float eps, 
	   float Vmin, float Vmax);
void wtcglsvel(float *t, float *qaxis, float *h,float *x,float *b,float *Qp, float tol);
void velopi( float *m, float *t, float *h, float *q, float *d); 
void velopi( float *m, float *t, float *h, float *q, float *d, float *ww);
void velopi( float *m, float *t, float *h, float *q, float *d, float *ww, float theta); 
void velop(float *m, float *t, float *h, float *q, float *d);
void velop(float *m, float *t, float *h, float *q, float *d, float theta);
void radonopi_id(float *m, float *t, float *h, float *q, float *d); 
void rho_filter (int npoints, int nt, float dt, float *rho);

void radonopi_lin(float *m, float *t, float *h, float *q, float *d);
void radonopi_par(float *m, float *t, float *h, float *q, float *d);

void radonop_lin(float *m, float *t, float *h, float *q, float *d);
void radonop_par(float *m, float *t, float *h, float *q, float *d);
void radon_param(float fmax,float *, int, float dx_av, float qmin, 
   float *qmaxt, float *qmax, float *dq, int nq, int rtmethod);
void radon_param(float fmax,float *, int,float dx_av, float qmin, 
   float *qmaxt, float *qmax, float *dq, int nq, int rtmethod, float factor);
void interval(float *pos, int nh, float *dx_max, float *dx_av);
void smoothing(float *d,int nt,int nx,int nl,int nr,int flag);

float testadjrad(float *t, float *h, float *q);

float testadjop(void (*oper) (float *,float *,float *,float *,float *,int ,int ,int,int), float *t, float *h, float *q, int nt, int nh, int nq);

float testadjop(void (*oper) (float *,float *,float *,float *,float *,float *,int ,int ,int,int), float *t, float *h, float *q, float *vel,int nt, int nh, int nq);

float testadjop(void (*oper) (float *,float *,float *,float *,float *,float **,int ,int ,int,int), float *t, float *h, float *q, float **vel,int nt, int nh, int nq);

float testadjop(void (*oper) (float *,float *,unsigned short **,int, int ,int, int, int),unsigned short  **index,int nt, int nh, int nq, int nsparse);

float testadjop(void (*oper) (float *,float *,unsigned int **,int, int ,int, int, int),unsigned int **index,int nt, int nh, int nq, int nsparse);

float testadjop(void (*oper) (float *,float *,float **,int ,int ,int, int),float **index,int nt, int nh, int nq);

float testadjop(void (*oper) (float *,float *,float **,int ,int ,int, int, int),float **index,int nt, int nh, int nq, int ns);


float testadjop(void (*oper) (float *,float *,float **,int ,int ,int, int, int, float *wavelet),float **index,int nt, int nh, int nq, int ns, float *wavelet, int nw);

float testadjop(void (*oper) (float *,float *,unsigned int **,int ,int ,int, int, int, float *wavelet, int nw),unsigned int **index,int nt, int nh, int nq, int nsparse, float *wavelet,int nw);

void wtcgls(void (*oper) (float *,float *,float *,float *,float *,int ,int ,int,int), int nt, int nh, int np, float *t, float *h, float *p, float *x,float *b,float *Wm, float *Wd, float tol, float step, float eps1, float eps2, int itercg);

void wtcgls(void (*oper) (float *,float *,float *,float *,float *,float *,int ,int ,int,int), int nt, int nh, int np, float *t, float *h, float *p, float *x,float *b,float *Wm, float *Wd, float *vel,float tol, float step, float eps1, float eps2, int itercg);

void wtcgls(void (*oper) (float *,float *,float *,float *,float *,float **,int ,int ,int,int), int nt, int nh, int np, float *t, float *h, float *p, float *x,float *b,float *Wm, float *Wd, float **vel,float tol, float step, float eps1, float eps2, int itercg);

void wtcgls(void (*oper) (float *,float *,unsigned short **,int ,int ,int, int, int), int nt, int nh, int np, int nsparse, float *x,float *b,float *Wm, float *Wd, unsigned short **index, float tol, float step, float eps1, float eps2, int itercg);

void radonhyp(float *m, float *t, float *h, float *q, float *d, int adj, int nt, int nh, int nq); 

void radonhyp(float *m, float *d, unsigned int **index, int adj, int nt, int nh, int nq, int nsparse, float *wavelet, int nw);

void radonhyp(float *m, float *t, float *h, float *q, float *d, float *vel, int adj, int nt, int nh, int nq);

void radonhyp(float *m, float *t, float *h, float *q, float *d, float **vel, int adj, int nt, int nh, int nq);
 
void radonhyp(float *m, float *d, unsigned int **index, int adj, int nt, int nh, int nq, int nsparse);

void radonhyp_crude(float *m, float *t, float *h, float *q, float *d, int adj, int nt, int nh, int nq);

void radonhyp_crude(float *m,float *t, float *h, float *q, float *d, float **vel,int adj,int nt, int nh, int nq);

void radonhyp_inv(float *m,float *t, float *h, float *q, float *d, float **vel,int adj,int nt, int nh, int nq);

void radonhyp_crude_inv(float *m,float *t, float *h, float *q, float *d, float **vel,int adj,int nt, int nh, int nq);

void radonpar(float *m, float *t, float *h, float *q, float *d, int adj, int nt, int nh, int nq); 

void radonlin(float *m, float *t, float *h, float *q, float *d, int adj, int nt, int nh, int nq); 

void enmo1(float *m,float *t,float *p,float *v,float *d, float *vel,int adj, int nt, int np, int nv);

int radhypsp(float *t,float *h, float *q, float **vel, int nt, int nh, int nq);

void radhypsp(float *t, float *h, float *q, float **vel, int nt, int nh, int nq,float **index);
 
void radhypsp(float *t, float *h, float *q, float **vel, int nt, int nh, int nq, unsigned int **index); 

void radonhyp(float *m, float *d, float **index, int adj, int nt, int nh, int nq);

void radonhyp(float *m, float *d, float **index, int adj, int nt, int nh, int nq);

void radint(float *t, float *q, float *h, float *h2, float *m,float *d,float *dint,float eps,float qmin,float qmax, float fmax, float *vel, float dperv, float pervmin);

void radint(float *t, float *q, float *h, float *h2, float *m,float *d,float *dint,float eps,float qmin,float qmax, float fmax, float *vel, float dperv, float pervmin, float alum, float blum, float clum, int norm, float t0);

void radint(float *t, float *q, float *h, float *h2, float *m,float *d,float *dint,float eps,float qmin,float qmax, float fmax, float *vel, float dperv, float pervmin, float alum, float blum, float clum, int norm, float t0, int mute, float parmute);

void lpcgnr(void (*oper) (float *,float *,unsigned short **, int, int ,int ,int,int), int nt, int nh, int np, int nsparse, float *x,float *b,float *Wm, float *Wd, float *M, unsigned short **index, float tol, float step, float eps1, float eps2, int itercg,int restart);

void mpcgnr(void (*oper) (float *,float *,unsigned short **, int, int ,int ,int,int), int nt, int nh, int np, int nsparse, float *x,float *b,float *Wm, float *Wd, float *M, unsigned short **index, float tol, float step, float eps1, float eps2, int itercg,int restart);

void Lumley_precond(float *M,float *t, float *h, int nt, int nh, float alum, 
		    float blum, float clum);
void rhofilt(float *d, float *t, float *h, int nt, int nh, float dt);

void save_vector(float *d, int n, char* name);
void save_gather(float **d, int nh, int nt, float dt, char* name);

float wpcgnr(void (*oper)(float *,float *,float *,float *,float *,float **,int ,int ,int,int),int nt, int nh, int np, float *t, float *h, float *p, float *x,float *b,float *Wd, float *M, float **vel, float tol, float step, int itercg,int restart);

float wpcgnr(void (*oper) (float *,float *,unsigned short **, int, int ,int ,int,int), int nt, int nh, int np, int nsparse, float *x,float *b,float *Wd, float *M, unsigned short **index, float tol, float step,int itercg,int restart);

float wpcgnr(void (*oper) (float *,float *,unsigned int **, int, int ,int ,int,int), int nt, int nh, int np, int nsparse, float *x,float *b,float *Wd, float *M, unsigned int **index, float tol, float step,int itercg,int restart);

float wpcgnr(void (*oper) (float *,float *,float **, int, int ,int ,int), int nt, int nh, int np, float *x,float *b,float *Wd, float *M, float **index, float tol, float step,int itercg,int restart);

float wpcgnr(void (*oper) (float *,float *,float **, int, int, int ,int ,int), int nt, int nh, int np, int ns, float *x,float *b,float *Wd, float *M, float **index, float tol, float step,int itercg,int restart);

void modelweight(float *m, int nx, int norm, float eps1, float *Wm);
void modelweight_inv(float *m, int nx, int norm, float eps1, float *Wm);

float mpcgne(void (*oper) (float *,float *,unsigned short **, int, int ,int ,int,int), int nt, int nh, int np, int nsparse, float *x,float *b,float *Wd, float *M, unsigned short **index, float tol, float step,int itercg,int restart);

float pcgnr(void (*oper) (float *,float *,unsigned short **, int, int ,int ,int,int), int nt, int nh, int np, int nsparse, float *x,float *b,float *Wd, float *Wmi, unsigned short **index, float tol, float step,int itercg,int restart);

void irreg_velaxis(int nq, int nt, float pervmin, float dperv, float *t, float *vel,float *q, float **vgrid);

void irreg_velaxis(int nq, int nt, float pervmin, float dperv, float *t, float *vel,float *q, float **vgrid, int centralq);

void irreg_slowness_axis(int nq, int nt, float pervmin, float dperv, float *t, float *vel,float *q, float **vgrid, int centralq);

void reg_slowness_axis(int nq, int nt, float dperv, float *t, float *vel,float *q, float **vgrid);

void build_index_rad(float *t, float *h, float *q, float **vel, int nt, int nh, int nq,unsigned int **index, int it0);

void build_index_rad_inv(float *t, float *h, float *q, float **vel, int nt, int nh, int nq,unsigned int **index, int it0);

void build_index_inv(float *t, float *h, float **vgrid, int nt, int nh, 
		       int nq,unsigned int **index, int nonzero, float t0);

int findgaps(float *h,float *h2,int nh,float tol, float hmin, float hmax);

int taper(float **data, int nt, int nh, int ntaper,int flag);
int taper(float *data, int nt, int nh, int ntaper,int flag);



void sethflag(int *hflag, int nhflag, float *h, int nh, float *h2, int nh2);

void sethflag2(int *hflag, int nhflag, float *h, int nh, float *h2, int nh2, float under, float over);

void plotgather(float **d, int nh, int nt, float dt, const char *s);
void plotvector(float *d, int nt, const char *s);
void contran(int conj,int sum, int nx, float *xx, int nb, float *bb,float *yy);
void contruc(int conj,int add, int lag, int nx, float *xx, int nb, float *bb,int ny,float *yy);
void convin(int conj,int sum, int nx, float *xx, int nb, float *bb,float *yy);
float testadjop_conv(void (*oper) (int, int, int, float *,int, float *,float *), 
		     int nb, float *bb, int nx, int ny);
segy cleansegy(segy tr);
#endif























