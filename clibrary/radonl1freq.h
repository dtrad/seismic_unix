#ifndef RADONL1FREQ_H
#define RADONL1FREQ_H

#include"su.h"
#include"dan.h"
#include "inversion_par.h"

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

typedef struct {	/* parameters for the double Radon operator */
  int rtmethod1;
  int rtmethod2;
  int nq1;
  int nq2;  
  float factor1;
  float factor2;
  float depth1;
  float depth2;   
  float dq1; 
  float dq2; 
  float qmin1;
  float qmin2;
  float fmax1;
  float fmax2;
  int mute1;
  int mute2;
} op_param;


void radonl1freq(float **data, float *h, int nh,float *t, int nt, float dt, float **model, 
		 float *q, float *vel, int itercg, int iter_end, op_param op2, float smute,
		 float nmofactor, float *ffilter, float *amps);


void radonl1freq_loop(float **data, float **model, float *t, float *q, float *h, 
		      int nt, int nq, int nh, float dt, int itercg, int iter_end, 
		      op_param op2);

void matrix(complex **L, float *h, float *q,int nh,int nq,float w, float dq, int rtmethod);

void radonl1freq_karmarkar(float *d, float **L, float *m, int nq, int nh, 
			   int iter_end, int itercg);

int PD_lp_bm(float *b, float **A, float *x, int nx, int ny,   
	     int iter_ext, float delta, float gamma, 
	     float PrI, float DI, float DG,int iter_end, int itercg);

int  CG_solver(float **A,float *D,int np,int nd,
	       float *x, float *b, float delta, int iter_end);

void matmul0(float *y, float **A, float *x, int adj, int ny, int nx);

float min3(float a, float b, float c);

void radon_param_beam(float fmax1, float fmax2, float *h, int nh,  float *q, float nq, 
		      float qmin1, float qmin2, float *q1, 
		      float *q2, float *ph1, float *ph2, int nq1, int nq2, float depth1, 
		      float depth2, int  rtmethod1, int rtmethod2, float factor1, 
		      float factor2, float *pdq1, float *pdq2);

void complex2real(complex *b, complex **A, float *br, float **AR, int nb, int nx);

void real2complex(complex *b, float *br, complex *x, float *xr, int nb, int nx);

void plotgather(float **d, int nh, int nt, float dt, const char *s);
void plotgatherb(float **d, int nh, int nt, float dt, const char *s);

#endif






