#ifndef RADONL1_H
#define RADONL1_H


void radon_PD_lp_bm(float **d, float **m, float *t, float *q, float *h, int nt, 
		    int nq, int nh, int iter_end, int itercg, int method);

int PD_lp_bm(int nx, int ny, int na, float *c, float *A, float *x, float *b, 
	     int *irow,int *icol,int iter_ext, float delta, float gamma, float PrI, 
	     float DI, float DG,int iter_end, int itercg);

void matmul0(int conj,int nb,float *bb, int nx, float *x, int ny, float *y, 
	     int *irow, int *icol);

void matmul1(int conj,int nb,float *bb, int nx, float *x, int ny, float *y, 
	     int *irow, int *icol);

float min3(float a, float b, float c);

int  CG_solver(float *A,float *D,int na,int np,int nd,int *irow,int *icol,
	       float *x, float *b, float delta, int iter_end);


int radonop(int nh, float *h, int nt, float *t, int nv, float *v, float *A, 
	    int *icol,int *irow, int vel);


#endif
