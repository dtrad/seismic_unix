#ifndef DAN_H
#define DAN_H

#define TRACE fprintf(stderr,"%s @ %u\n",__FILE__,__LINE__)

#define TRACE2 { \
  fprintf(stderr,"%s @ %u\n",__FILE__,__LINE__); \
  fflush(stderr); \
}
#define TRACEF(a) { \
  fprintf(stderr,"%s @ %u\n",__FILE__,__LINE__); \
  fprintf a; \
  fflush(stderr); \
}

#include"su.h"
#include"cwp.h"
#include"segy.h"
segy cleansegy(segy tr);
void xtimesA(complex *b,complex *x,complex **A,int nr,int nc);
void xtimesA(complex *b,float *x,complex **A,int nr,int nc);
void Atimesx(complex *b,complex **A,complex *x,int nr,int nc);
void Atimesx(float *b,float **A,float *x,int nr,int nc);
void Atimesx(int nr,int nc,float *b,float **A,float *x);
void Atimesx(complex *b,complex **A,complex *x,int nr,int nc, int adj);
void Atimesx(complex **A, float x, int n1, int n2);
void Atimesx(float **A, float x, int n1, int n2);
void Atimesx(float *y,float **A,float *x, int adj, int ny,int nx); /* Old not use */
void Atimesx2(float *y,float **A,float *x, int ny,int nx, int adj); /* For v(0,n-1) */
void Atimesxnr(float *y,float **A,float *x, int ny,int nx, int adj); /* For v(1,n) */
void xplusy(complex *z,complex *x,complex *y,int n);
void AequalB(float **A,float **B,int n1, int n2);
void AequalB(complex **A,complex **B,int n1, int n2);
void Aplus_equalB(float **A,float **B,int n1, int n2);
void Aplus_equalB(complex **A,complex **B,int n1, int n2);
void xminusy(complex *z,complex *x,complex *y,int n);
void xminusy(float *z,float *x,float *y,int n);
complex cdot(complex *x,complex *y,int n);
float rcdot(int n, complex *a, complex *b);
void xequaly(complex *x,complex *y,int n);
void xequaly(float *x,float *y,int n);
void xtimesy(complex *z,complex *x,complex *y,int n);
void xtimesy(float *z,float *x,float *y,int n);
void xtimesy(complex *z,float *x,complex *y,int n);
void xtimesy(complex *y,complex *x,float a,int n);
void xtimesy(complex *z,complex *x,float *y,int n);
/* AtimesBm.cpp */
void AtimesBm(complex **C,complex **A,complex **B,int nr,int ni,int nc);
void AtimesBm(complex **C,complex **A,complex **B,int nr,int ni,int nc,
 char ch);
void AtimesBm(complex **C,complex **A,complex **B,int nr,int ni,int nc,
 complex *D);
void AtimesBm(complex **C,complex **A,complex **B,int nr,int ni,int nc, 
	      complex *D, char ch);
void diagxtimesA(complex **B, float *x, complex **A, int nr, int nc);
void Atimes_conjB_elem(complex **C,complex **A,complex **B,int n1, int n2);
void Atimes_B_elem(complex **C,complex **A,complex **B,int n1, int n2);
void AtimesBH(complex **C,complex **A,complex **B,int nr1, int nc1, int nr2, int nc2);
void AHtimesB(complex **C, complex **A, complex **B, int nr1, int nc1, int nr2, 
	      int nc2);

void transpose(complex **F, complex **FT, int nr, int nc);
/********************/
float dot(int n, float *a, float *b);
float freqweight(int j, float df, float f1, float f2);
void freqweighting(complex **d, int nfreq, int nh, float df, float fmin, float fmax);
float normalize(int n, float *a);
float normalize(int n, complex *a);
void plotcurves(float *d, int n1, int n2, const char *s);
void plotcurves(float **d, int n1, int n2, int f2, int l2, const char *s);
void plotgather(float *d, int n1, int n2, const char *s);
void plotgather(float **d, int n1, int n2, const char *s, float perc);
void plotgather(float **d, int nh, int nt, float dt);
void plotgather(float **d, int nh, int nt, float dt, const char *s);
void gather_pipe(float **d, int n2, int n1, const char *s);
void plotgather_pipe(float **d, int n2, int n1, const char *s);
void plotgather_pipe(float *d, int n2, int n1, const char *s);
void save_gather(float ***d, int nv, float *vel, int nh, float *h, int nt, float *t,
		 char* name);
void save_gather(float ***d, int nv, float *vel, int nh, float *h, int nt, float dt,
		 char* name);
void save_gather(float ***d, int nx, int nh, int nt, float dt, char* name);
void save_gather(float **d, int nh, int nt, float dt, char* name);
void save_gather(float **d, int nh, float *h, int nt, float dt, char* name);
void save_gather(float *d, int nh, float *h, int nt, float dt, char* name);

void save_gather_tx(float **d, int nh, float *h, int nt, float dt, char* name);
void save_vector(float *d, int n, char* name);
void save1dfile(float *d, int n1, const char *s);
void save2dfile(float **d, int n1, int n2, const char *s, int transpose);
void save2dfile(float **d, int nh, int nt, float dt, const char *s);
void save2dfile(complex **d, int n1, int n2, const char *s);
void save3dfile(complex ***d, int n1, int n2, int n3, const char *s);
void read1dfile(float *d, int n1, const char *s);
void read2dfile(complex **d, int n1, int n2, const char *s);
void read3dfile(complex ***d, int n1, int n2, int n3, const char *s);

void zero_array(complex **d,int n2, int n1);
void zero_array(float **d,int n2, int n1);
void zero_array(complex ***d,int n3, int n2, int n1);
void zero_array(float ***d, int n3, int n2, int n1);

void zero_vector(complex *d,int n);
void zero_vector(float *d,int n);

int read_ascii_file(const char *name,float *x);
/***********************************************************************************/
/* File: fft.cpp */
/* 1D fft */
void fft_parameters(int nt, float dt, int *pnfft, int *pnf, float *pdf);
/*** Gather ***/
void fftgo_xt2fx(int sign,float **d,complex **m,int nh,int nt, float dt, int nfft, int nf);
int fftback_fx2xt(int sign,float **d,complex **m,int nh,int nt,float dt,int nfft, int nf0);
void fftgo_xt2fx(int sign,float *d,complex **m,int nh,int nt, float dt, int nfft, int nf);
int fftback_fx2xt(int sign,float *d,complex **m,int nh,int nt,float dt,int nfft, int nf0);
/** Single trace **/
void fftgo_xt2fx(int sign,float *d,complex  *m, int nt, float dt, int nfft, int nf);
int fftback_fx2xt(int sign,float *d,complex  *m, int nt, float dt, int nfft, int nf);

/* 2D FFT */
void fft2_parameters(int nt, float dt, int nx, float dx, int *pnw, float *pdw,
		     int *pnk, float *pdk, float *w, float *k);

void fft2_xt2kf(float **dataxt, complex **datakf, int nt, int nx, int nw, int nk);
void fft2_kf2xt(float **dataxt, complex **datakf, int nt, int nx, int nw, int nk);
void fksize(int nt, int nx, int *pnw, int *pnk);
/***********************************************************************************/
/* Max min : normalize.cpp*/
float max(int n, float *a);
float max(int n, complex *a);
float min(int n, complex *a);
float min(int n, float *a);
float min_with_sign(int n, float *a);
void maxmin(float *x,int lx, float *pmin, float *pmax);

/* Windowing */
float **window(float **a, int a2f, int a2l, int a1f, int a1l);
complex **window(complex **a, int a2f, int a2l, int a1f, int a1l);
void window_cpy(float **a, int a2f, int a2l, int a1f, int a1l, float **m, int go);
float **eallocwindow(int a2f, int a2l, int a1f, int a1l);

/* Filtering */
void filtoutliers(float **data,int nh, int nt, float quantil);
void filtering(complex **d, int n2, int n1, int nfft, float dt, float *f, float *amps, int npoly);
void polygonalFilter(float *f, float *amps, int npoly, int nfft, float dt, float *filter);

int findpower2(int nt);

void padding_zeros(float *d, int n2, int n1, int n1z);

//float quest(float p, int n, float x[]);
#endif











