#ifndef CLIBRARY_H
#define CLIBRARY_H
#ifndef TRACE
#define TRACE fprintf(stderr,"%s @ %u\n",__FILE__,__LINE__)
#endif


#define floatprint(expr) fprintf(stderr,#expr " = %f\n", expr)
#define expfprint(expr) fprintf(stderr,#expr " = %g\n", expr)
#define intprint(expr) fprintf(stderr,#expr " = %d\n", expr)
void hrrtf(float *pos, float **d ,float dt,float eps1,float eps2, float qmin, 
	   float **m, float *q, float dq, float freq, float eps, float *bad);
void hrrtf(float *pos, float **d ,float dt,float eps1,float eps2, float qmin, 
	   float **m, float *q, float dq, float freq, float eps, float *bad,
	   int taper);
void AtimesDiag(complex **C,complex **A,complex *Diag,int nr,int nc);
void AtimesDiag(complex **C,complex **A,float *Diag,int nr,int nc);
void Atimesx(complex *b,complex **A,complex *x,int nr,int nc);
void Atimesx(float *b,float **A,float *x,int nr,int nc);
void Atimesx(int nr,int nc,float *b,float **A,float *x);
void Atimesx(complex *b,complex **A,complex *x,int nr,int nc, int adj);
void htmul(int n, complex *a, complex *x, complex *y);
void xminusy(complex *z,complex *x,complex *y,int n);
void xminusy(float *z,float *x,float *y,int n);
complex mycdot(complex *x,complex *y,int n);
float rcdot(int n, complex *a, complex *b);
void xequaly(complex *x,complex *y,int n);
void xtimesy(complex *z,complex *x,complex *y,int n);
void xtimesy(float *z,float *x,float *y,int n);
float modgrad(complex *u, int nq, int norm, float powd, complex *Qp);
float modgrad(complex *u, int nq, int norm,float powd, complex *Qp, int method);
void xplusy(complex *z,complex *x,complex *y,int n);
complex cdot(complex *x,complex *y,int n);
float costfunc(complex *d, complex *dc, complex *u,int nh,int np,
float sigma,int norma);
void displayA(complex **A,int nr,int nc);
void displayA(complex *A,int n);
void displayA(float *A,int n);
float testchol(complex **A,double *p,int nr, int nc);
float maxmax(complex **A,int nr, int nc);
float maxmax(complex *A,int n);

void fftgo(int sign,float **,complex **, int nh, int nt, float dt, int *nf0);
void fftgo0(int sign,float **,complex **, int nh, int nt, float dt, int *nf0);
void fftgo1(int sign,float *,complex **, int nh, int nt, float dt, int nf0);

int fftback(int sign,float **,complex  **, int nh, int nt, float dt, int nf0);
int fftback0(int sign,float **,complex  **, int nh, int nt, float dt, int nf0);
int fftback1(int sign,float *,complex  **, int nh, int nt, float dt, int nf0);

float real(complex);
void cholsl(complex **a, int n,double p[],complex b[],complex x[]);
void cholsl(float **a, int n, float p[], float b[], float x[]);
void choldc(float **a, int n, float p[]);
void choldc(complex **a, int n, double p[]);

void AtimesBm(complex **C,complex **A,complex **B,int nr,int ni,int nc);
void AtimesBm(complex **C,complex **A,complex **B,int nr,int ni,int nc,
 char ch);
void AtimesBm(complex **C,complex **A,complex **B,int nr,int ni,int nc,
 complex *D);
void AtimesBm(complex **C,complex **A,complex **B,int nr,int ni,int nc, 
	      complex *D, char ch);
//void p_stack(float *pos, float **d,float dt,float **m,float *q, float dq,
//	     float eps1, float eps2, float eps, float fmax);
void toepradon(float *pos, float **d,float dt,float **m,float *q, float dq,
	     float eps1, float eps2, float eps, float fmax);
void toepradoncg(float *pos, float **d,float dt,float **m,float *q, float dq,
	     float eps1, float eps2, float eps, float fmax);
void cholund(float *pos, float **d,float dt,float **m,float *q, float dq,
	     float eps1, float eps2, float eps, float fmax);
void cholover(float *pos, float **d,float dt,float **m,float *q, float dq,
	     float eps1, float eps2, float eps, float fmax);
void norma(float **m, int nt, int nq);
void radon_param(float fmax,float *,int, float dx_av, float qmin, 
   float *qmaxt, float *qmax, float *dq, int nq, int rtmethod);
void radon_param(float fmax,float *,int, float dx_av, float qmin, 
   float *qmaxt, float *qmax, float *dq, int nq, int rtmethod, float factor);

void interval(float *pos, int nh, float *dx_max, float *dx_av);
int ctoeplitz( int n, complex *r, complex *a, complex *b,
		 complex *f, complex *g );
int ctoephcg( int niter, int n, complex *a, complex *x, complex *y, 
	complex *s, complex *ss, complex *g, complex *rr);
int cghrrt( int niter, int n, complex **a, complex *x, complex *y, 
	complex *s, complex *ss, complex *g, complex *rr);
void hrrti(float **d,float *pos,float dt,float **m,float *q, float fmax);
float misfit(complex *a,complex *b, float eps1,float *dh, int n);
float misfit(complex *a,complex *b, float eps1,float *dh, int n, int norm);
float modnorm(complex *u,float eps2,float dq, int n);
float freqweight(int j, float df, float f1, float f2);

//void conjgradrt(float *pos,float **data,float dt,float **model,float *q,
//		float dq,float eps1, float fmax);

//void conjgradrt2(float *pos,float **data,float dt,float **model,float *q,
//		 float dq,float eps1, float fmax);

void chol_all(complex *R,complex **LL,complex *Qp,int n,complex *b,
              complex *x,double *diagll);

void choleskyrt(float *pos, float **d, float dt,float **m, float *q, float dq,
	     float eps1, float eps2, float fmax);
void matrix_4(complex *R, complex **l,complex **lh,float *pos,float *q,int nh,int nq,float w,float *dh, float dq, int rtmethod, float *g);

void matrix_5(complex *R, complex **l,complex **lh,float *pos,float *q,int nh,int nq,float w,float *dh, float dq, int rtmethod, float *g);

void matrix_3(complex *R, complex **l,complex **lh,float *pos,
	     float *q,int nh,int nq,float w,float *dh, float dq, int rtmethod);
void matrix_3(complex *R, complex **l,complex **lh,float *pos,
	      float *q,int nh,int nq,float w,float *dh, float dq, int rtmethod, float *Cd);
void matrix_3(complex **l,complex **lh,float *pos,float *q,int nh,int nq,float w,float *dh, float dq, int rtmethod);
void matrix_3(complex **l,float *pos,float *q,int nh,int nq,float w, int rtmethod);
void matrix_3(complex **l,float *pos,float *q,int nh,int nq,float w, int rtmethod, float *g);
void matrix_3(complex **l,complex **lh,float *pos,float *q,int nh,int nq,float w,float *dh, float dq, int rtmethod, int flag);
void matrix_3(complex **l,complex **lh,float *pos,float *q,int nh,int nq,float w,float *dh, float *Wd, float dq, int rtmethod, float *g);
void conjgradrt3(float *pos,float **data,float dt,float **model,float *q,
		 float dq,float eps1,float eps2,float eps,float fmax);
void conjgradrt4(float *pos,float **data,float dt,float **model,float *q,
		 float dq,float eps1,float eps2,float eps,float fmax, 
		 float *Wd);
void  circ_mult(int n,complex *RC,complex *f,complex *g,int nf,
		complex *fc, complex *gc);
void svdcmp(float **a, int m, int n, float w[], float **v);
void svbksb(float **u, float w[], float **v, int m, int n, float b[], float x[]);
void svdover(float *pos, float **d, float dt,float **m, float *q, float dq,
	     float eps1, float eps2, float eps, float fmax);
     //#include "/home/dtrad/radon/clibrary/nrutil.h"
void conjgrad3r(float *pos, float **data, float dt,float **model, float *q, float dq,float eps1, float eps2, float eps, float fmax);
void linbcg1(double **A, double x[], double b[],int n,int itol,double tol,
	int itmax, int *iter, double *err);

void cgls(complex *d,complex **L, complex **LH,complex *u,complex *Qp,
	  int nh,int nq,complex *q,complex *r,complex *raux,complex *s,
	  complex *p, float eps1,float eps2, float eps);
void cgls0(float *pos, float **data, float dt,float **model, float *q, float dq,float eps1, float eps2, float eps, float fmax);
void wtcgls0(float *pos, float **data, float dt,float **model, float *q, float dq,float eps1, float eps2, float eps, float fmax, float *bad);
void wtcgls0(float *pos, int nh, float **data, float *t, int nt, float dt,float **model, float *q, int nq, float dq, float eps1, float eps2, float eps, float fmax, float *Wd, int itercg, int iter_end,int norm,float step, int testadj);

void wtcgls0_td(float *pos, int nh, float **data, float *t, int nt, float dt,float **model, float *q, int nq, float dq, float eps1, float eps2, float eps, float fmax, float *Wd, int itercg, int iter_end,int norm,float step, int testadj);

void radoncg_levinson(float *pos, int nh, float **data, float *t,int nt,float dt,float **model, float *q, int nq, float dq,float eps1, float eps2, float eps, float fmax, float *Wd, int itercg, int iter_end, int norm, float step, int testadj,float quantil, int rtmethod, float depth);

void radon_cholund(float *pos, int nh, float **d,float *t,int nt, float dt,float **m, float *q, int nq,float dq,float eps1, float eps2, float eps, float fmax,float *Wd, int itercg, int iter_end, int norm, float step, int testadj);

void wtcgls(complex *b,complex **L, complex **LH,complex *x,complex *Wm,
          float *Wd, int nh,int nq,complex *q,complex *q1,complex *s,
	  complex *x1,complex *z,complex *z1,complex *r,complex *Az,
	  float *eta,float *rho, float *gcv, float tol);
void wtcgls(complex *b,complex **L, complex **LH,complex *x,complex *Wm,
	  float *Wd,int nh,int nq,complex *q,complex *q1,complex *s,
	  complex *x1,complex *z,complex *z1,complex *r,complex *Az,
	  float *eta,float *rho, float *gcv, float tol, float step, int itercg);


void frprmn(float p[], int n, float ftol, int *iter, float *fret,
	float (*func)(float []), void (*dfunc)(float [], float []));

//void dfunc1(float *x, float *xc);
//float func1(float *x);
void cgnl0(float *pos, float **d, float dt,float **m, float *q, float dq,
	     float eps1, float eps2, float eps, float fmax);

void lsqr0(float *pos, float **data, float dt,float **model, float *q, float dq,float eps1, float eps2, float eps, float fmax);
void lsqrc(complex **L, complex **LH,complex *x,complex *b,float tol,
float eps2,int reorth,complex **U,complex **B,complex **V,complex *Qp,complex *u,complex *v,complex *w,complex *q,complex *p, complex *r,complex *temp2,float *eta,
float *rho,float *gcv);

void matrix_3(complex **l,complex **lh,float *pos,float *q,int nh,int nq,float w,float *dh, float *dh2, float dq, int rtmethod, int flag);

float testadj_rad_f(complex **L, complex **LH);

int wpcgnr(void (*oper) (complex *,complex *, complex  **, int, int, int), int nh, int np, complex *x,complex *b,float *Wd, float *M, complex  **L, float tol, float step,int itercg,int restart);

void radon_freq (complex *,complex *, complex  **, int, int, int);

void rad_wpcgnr(float *pos, float **data, float dt,float **model, float *q, float dq,float eps1, float eps2, float eps, float fmax, float *Wd);

void modelweight(complex *m, int nx, int norm, float eps1, float *Wm);
void savesudata(float **d, int nt, int nh, float dt, char *name, int shape);
void savesudata(float *d, int nt, int nh, float dt, char *name);

void sethflag2(int *hflag, int nhflag, float *h, int nh, float *h2, int nh2, float under, float over);

void save_gather(float **d, int nh, float *h, int nt, float dt, char* name);

void save_gather(float *d, int nh, int nt, float dt, char* name);

void save_gather(float **d, int nh, int nt, float dt, char* name);

void save_gather(complex **d, int nh, float *h, int nt, float dt, int nfft, char* name);

void radon_matrix(complex *R, complex **l,complex **lh,float *h,float *q,int nh,int nq,float w,
		  float *dh);

void radon_matrix(complex **l, float *h,float *q,int nh,int nq,float w);

#endif
