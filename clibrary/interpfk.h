#ifndef INTERPFK_H
#define INTERPFK_H

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

void save_gather(float **d, int nh, int nt, float dt, char* name);

void save2d(float **d, int nh, int nt, float dt, char* name);

void stolt1k (float k, float v, float s, float fmax, int nt, float dt,
		     complex *p, complex *q);
void stolt1kadj (float k, float v, float s, float fmax, int nt, float dt,
			complex *p, complex *q);
int taper(float **data, int nt, int nh, int ntaper,int flag);
void taper (int lxtaper, int lbtaper, 
		   int nx, int ix, int nt, float *trace);
void plotAmpSpec(complex **pk, int nk, int nt, float dt);
void dft_matrix(complex **F,float *h,float *k,int nh,int nk);
void save_gather_freq(complex **d, int nh, float *h, int nt, float dt, int nfft, 
		      char* name);

void stoltzop2(float **data, float **model, int nt, int nh, int nx, float *t, 
	       float *h, float *x, float vel, complex **F, complex **F2, float *wavelet, 
	       int nw, int adj);

void stoltzopinv2(float **data, float **model, int nt, int nx, float *t, 
		  float *x, float vel, complex **F, float *wavelet, int nw);

float adjteststoltz(int nt, int nh, int nq, float *t, float *h, float *q, float vel,
		    complex **F, complex **F2, float *wavelet, int nw);
void Atimesx_DFT(complex *y,complex **A,complex *x,int ny,int nx, int adj);

void stoltz_wtcgls(float **data, float **model, float *h, int nh,  float *t, int nt,
		   float *x, int nx, float vel, inv_par inv, complex **F, complex **F2, 
		   float *wavelet, int nw);

void stoltz_wtcgls(float **data, float **model, float **Wd, float *h, int nh,  float *t, int nt,
		   float *x, int nx, float vel, inv_par inv, complex **F, complex **F2, 
		   float *wavelet, int nw);

float wpcgnr_mig(void (*oper)  (float **data, float **model, int nt, int nh, int nq, 
				float *t,float *h, float *q, float vel, complex **F, 
				complex **F2, float *wavelet,int nw,int adj),
		 float **x, float **b,int nt, int nh, int nq, float *t, float *h, float *q, 
		 float vel, float **Wm, float **Wd, inv_par inv, complex **F, complex **F2,
		 float *wavelet, int nw);

void deviations(float *m, int nx, float *d, int ny, int norm, float quantil1, float quantil2, 
		float *sigmam, float *sigmad);

void weights(float *m, int nx, int norm, float sigmam, float *Wm, int iter);

void CequalAxB(float **C,float **A,float **B,int n1,int n2);


void CequalApluskxB(float **C, float **A, float k, float **B, int n1,int n2);

void LU_decomposition (int nrows, complex **matrix, int *idx, float *d);

void backward_substitution (int nrows, complex **matrix, int *idx, complex *b);

void inverse_matrix_multiply (int nrows1, complex **matrix1, int ncols2,
			      int nrows2, complex **matrix2, complex **out_matrix);
void dftls(complex **F, complex **FLST, int nh, int nk, float eps);

void AequalB(float **A,float **B,int n1, int n2);

void AequalB(complex **A, complex **B,int n1, int n2);

void FTmatrix(complex **F, complex **FLS, float *h, float *k, int nh, int nk, float eps);

void kaxis(int nx, float vmax, float dt, int nt, float dx, int *pnk, float *pdk);

void kaxis(int nx, float vmax, float dt, int nt, float dx, int *pnk, float *pdk, int dft);

void dft_matrix(complex **F,float *h,float *k,int nh,int nk);

void Atimesx_DFT(complex *y,complex **A,complex *x,int ny,int nx, int adj);

void save_gather_test(float **d, int nh, int nt, float dt, char* name);

/* stretching */
void makev (int nmig, float *tmig, float *vmig, float vscale,
	    int nt, float dt, float ft, float **v, float *vmin, float *vmax);
void makeut (float vstolt, float fmax, float *vrms,
	     int nt, float dt, float **ut, int *nu, float *du, float **tu);
void makeu (float vstolt, float *v, int nt, float dt, float *u);

void stretch(float **data, float **datas, int nt, int nu, int nh, float *t, float *tu, float *ut,
	     float dt, float du, int str);

/****************************************************************************************/
 

void contruc_2(int conj,int add, int nx, float *xx, int nb, float *bb,float *yy);

int get_wavelet(float *wavelet,const char *name,int nw, int type, float dt, float fpeak);

void convgather(float *d, float *m, int nh, int nt, float *wavelet, int nw, int adj);

float adjtestconvgather(int nh, int nt, float *wavelet, int nw);
void plotvector(float *d, int nt, const char *s);

int fft2_zeropad(float **datain, float **dataout, float **Wd, int nt, int nx, 
		 float interpfact);

int add_zerotraces(float **datain, float **dataout, float **Wd, int nt, int nx, 
		   float interpfact);

int add_zerotraces(float **datain, float **dataout, float **Wd, int nt, int nx, int nx2, 
		   float *h, float *h2);

void AtimesB(float **A, float **B, int nh, int nt);

typedef struct {	/* Parameters for inversion  */
  float tmin; 
  float tmax;
  int ihmin;
  int ihmax;
  float slope;
  float threshold;
} mutemask_par;

void mutemask(float **M, float **d, int nh, int nt, float dt, mutemask_par mutepar);


#endif



