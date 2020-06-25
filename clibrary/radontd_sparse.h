#include "inversion_par.h"

#ifndef RADONTD_SPARSE_H
#define RADONTD_SPARSE_H

#ifndef TRACE
#define TRACE fprintf(stderr,"%s @ %u\n",__FILE__,__LINE__)
#endif


/* FUNCTION PROTOTYPES */

void radontd_sparse(float *t, float *q, float *h, float **m,float **d, int nt, int nh, int nq, float dt, float *vel, float dperv, float pervmin, float t0, inv_par inv, int centralq, int filtout, int nw, float fpeak, int typewav, int LI, int nreg, float parmute, int mute);

void radontd_sparse_win(float *t, float *q, float *h, float **m,float **d, int nt, int nh, int nq, float dt, float **vgrid, inv_par inv, int dataprec, int nw, float fpeak, int typewav, int LI, float parmute, int mute);


void build_index_slowness(float *t, float *h, float *q, float **vel, int nt, int nh, int nq,unsigned int **index, int it0);

void irreg_slowness_axis(int nq, int nt, float pervmin, float dperv, float *t, float *vel,float *q, float **vgrid, int centralq);

void irreg_slowness_axis(int nq, int nt, float pervmin, float dperv, float *t, float *vel,float *q, float **vgrid, int centralq, int nreg);

float wpcgnr(void (*oper) (float *,float *,unsigned int **, int, int ,int ,int,int), int nt, int nh, int np, int nsparse, float *x,float *b,float *Wd, float *M, unsigned int **index, inv_par inv);

float wpcgnr(void (*oper) (float *,float *,unsigned int **, int, int ,int ,int,int, float *wavelet, int nw), int nt, int nh, int np, int nsparse, float *x,float *b,float *Wd, float *M, unsigned int **index, inv_par inv, float *wavelet, int nw);

float wpcgnr(void (*oper)(float *,float *,float *,float *,float *,float **,int ,int ,int,int),int nt, int nh, int np, float *t, float *h, float *p, float *x,float *b,float *Wd, float *M, float **vel, inv_par inv);

void plotgather(float **d, int nh, int nt, float dt);

int read_ascii_file(const char *name,float *x);

void getvelocities(float dt, int nt, int ncdp, float *cdp, float **ovv);

void save_gather(float **d, int nh, int nt, float dt, char* name);

void save2dfile(float **d, int nh, int nt, float dt, const char *s);

void interpovv (int nt, int ncdp, float *cdp, float **ovv, float cdpt, float *ovvt);

void radonhyp(float *m,float *t, float *h, float *q, float *d, float **vel,int adj,int nt, int nh, int nq);

void radonhyp_sinc(float *m,float *t, float *h, float *q, float *d, float **vel,int adj,int nt, int nh, int nq);

/********************* Protoypes for sinc **************************/

void radontd_sparse(float *t, float *q, float *h, float **m,float **d, int nt, int nh, int nq, float dt, float *vel, float dperv, float pervmin, float t0, inv_par inv, int centralq, int filtout, int nw, float fpeak, int typewav, int LI);

void radonhyp(float *m, float *d, float **index, int adj, int nt, int nh, int nq, int ns);

void build_index_slowness(float *t, float *h, float *q, float **slow2, int nt, int nh, int nq,
			  float **index, int it0);

float wpcgnr(void (*oper) (float *, float *,float **, int, int ,int ,int,int), int nt, int nh, int np, int nsparse, float *x,float *b,float *Wd, float *M, float **index, inv_par inv);

float testadjop(void (*oper) (float *,float *,float **,int ,int ,int, int, int),float **index,int nt, int nh, int nq, int nsparse);

float wpcgnr_m0(void (*oper) (float *,float *,unsigned int **, int, int ,int ,int,int), int nt, int nh, int np, int nsparse, float *x,float *b,float *Wd, float *M, unsigned int **index, inv_par inv);

/****************************************************************************************/

/********************* Protoypes for LI **************************/


void build_index_slowness_li(float *t, float *h, float *q, float  **vel, int nt, int nh, int nq, unsigned int **index, int it0);


void radonhyp_li(float *m, float *d, unsigned int **index, int adj, int nt, int nh, int nq, int ns);

void radonhyp_li(float *m, float *d, unsigned int **index, int adj, int nt, int nh, int nq, int ns, float *wavelet, int nw);

/****************************************************************************************/


void contruc_2(int conj,int add, int nx, float *xx, int nb, float *bb,float *yy);

int get_wavelet(float *wavelet,const char *name,int nw, int type, float dt, float fpeak);

/***********************Protoptypes fo radontd_sparse_interp ****************************/

void radontd_sparse(float *t, float *q, float *h, float **m,float **d, int nt, int nh, int nq, float dt, float *vel, float dperv, float pervmin, float t0,  inv_par inv, int centralq, int filtout, int nw, float fpeak, int typewav, int LI);

float wpcgnr(void (*oper) (float *,float *,float **, int, int ,int ,int,int, float *wavelet, int nw), int nt, int nh, int np, int nsparse, float *x,float *b,float *Wd, float *M, float **index, inv_par inv, float *wavelet, int nw);

float testadjop(void (*oper) (float *,float *,float **,int ,int ,int, int, int), float **index,int nt, int nh, int nq, int nsparse);

float testadjop(void (*oper) (float *,float *,float **,int ,int ,int, int, int, 
			      float *wavelet, int nw),float **index,int nt, int nh, 
		int nq, int nsparse, float *wavelet, int nw);

void radonhyp(float *m, float *d, float **index, int adj, int nt, int nh, int nq, int nsparse);

void radonhyp(float *m, float *d, float **index, int adj, int nt, int nh, int nq, int nsparse, float *wavelet, int nw);

/****************************************************************************************/

#endif





