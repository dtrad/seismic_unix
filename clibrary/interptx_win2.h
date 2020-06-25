#include "inversion_par.h"
#include "dan.h"

#ifndef RADONTD_WIN_H
#define RADONTD_WIN_H

#ifndef TRACE
#define TRACE fprintf(stderr,"%s @ %u\n",__FILE__,__LINE__)
#endif


/* FUNCTION PROTOTYPES */

void interptx_sparse(float *t, float *q, float *h, float **m,float **d, float **d,int nt, int nh, int nq, float dt, float *vel, float dperv, float pervmin, float t0, inv_par inv, int centralq, int filtout, int nw, float fpeak, int typewav, int LI, int nreg, float parmute, int mute, int itm, int plot);

void interptx_sparse_win(float *t, float *q, float *h, float **m,float **d, float **dp, int nt, int nh, int nq, float dt, float *vel, inv_par inv, int dataprec, int nw, float fpeak, int typewav, int LI, float parmute, int mute, int itm, int plot);

void build_index_slowness(float *t, float *h, float *q, float *vel, int nt, int nh, int nq,unsigned int **index);

void irreg_slowness_axis(int nq, int nt, float pervmin, float dperv, float *t, float *vel,float *q, float **vgrid, int centralq);

void irreg_slowness_axis(int nq, int nt, float pervmin, float dperv, float *t, float *vel,float *q, float **vgrid, int centralq, int nreg);

float wpcgnr(void (*oper) (float *,float *,unsigned int **, int, int ,int ,int,int), int nt, int nh, int np, int nsparse, float *x,float *b,float *Wd, float *M, unsigned int **index, inv_par inv);

float wpcgnr(void (*oper) (float *,float *,unsigned int **, int, int ,int ,int,int, float *wavelet, int nw), int nt, int nh, int np, int nsparse, float *x,float *b,float *Wd, float *M, unsigned int **index, inv_par inv, float *wavelet, int nw);


void plotgather(float **d, int nh, int nt, float dt);

int read_ascii_file(const char *name,float *x);

void getvelocities(float dt, int nt, int ncdp, float *cdp, float **ovv);

void save_gather(float **d, int nh, int nt, float dt, char* name);

void save2dfile(float **d, int nh, int nt, float dt, const char *s);

void interpovv (int nt, int ncdp, float *cdp, float **ovv, float cdpt, float *ovvt);

void radonhyp(float *m,float *t, float *h, float *q, float *d, float *vel, int adj,int nt, int nh, int nq);

void radonhyp_sinc(float *m,float *t, float *h, float *q, float *d, float *vel,int adj,int nt, int nh, int nq);

void rest_multiples(void  (*oper) (float *, float *, unsigned int **, int , int , int, 
				   int, int, float *, int), 
		    float **data, float **model, unsigned int **index, int adj, int nt, 
		    int nh, int nq, int nsparse, float *wavelet, int nw, float parmute, 
		    float *q, float *h, int itm, int plot);


void rest_multiples(void  (*oper) (float *, float *, unsigned int **, int , int , int, 
				   int, int), 
		    float **data, float **model, unsigned int **index, int adj, int nt, 
		    int nh, int nq, int nsparse, float parmute, float *q, float *h, int itm,
		    int plot);

/********************* Protoypes for sinc **************************/

void radontd_sparse(float *t, float *q, float *h, float **m,float **d, int nt, int nh, int nq, float dt, float *vel, float dperv, float pervmin, float t0, inv_par inv, int centralq, int filtout, int nw, float fpeak, int typewav, int LI);

void radonhyp(float *m, float *d, float **index, int adj, int nt, int nh, int nq, int ns);

void build_index_slowness(float *t, float *h, float *q, float **slow2, int nt, int nh, int nq,
			  float **index);

float wpcgnr(void (*oper) (float *, float *,float **, int, int ,int ,int,int), int nt, int nh, int np, int nsparse, float *x,float *b,float *Wd, float *M, float **index, inv_par inv);

float testadjop(void (*oper) (float *,float *,float **,int ,int ,int, int, int),float **index,int nt, int nh, int nq, int nsparse);

float wpcgnr_m0(void (*oper) (float *,float *,unsigned int **, int, int ,int ,int,int), int nt, int nh, int np, int nsparse, float *x,float *b,float *Wd, float *M, unsigned int **index, inv_par inv);

/****************************************************************************************/

/********************* Protoypes for LI **************************/


void build_index_slowness_li(float *t, float *h, float *q, float  *vel, int nt, int nh, int nq, unsigned int **index);


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
void smute_gather(float **d, int nt, int nh, float *t, float *h, float *vel, float smute);

void weights_td(float *m, int nx, int norm, float sigmam, float *Wm, int iter);

void deviations(float *m, int nx, float *d, int ny, int norm, float quantil1, float quantil2, float *sigmam, float *sigmad);

typedef struct {	/* Parameters for inversion  */
  float tmin; 
  float tmax;
  int ihmin;
  int ihmax;
  float slope;
  float threshold;
} mutemask_par;

void mutemask(float **M, float **d, int nh, int nt, float dt, mutemask_par mutepar);

void AtimesB(float **A, float **B, int nh, int nt);


void taper (int lxtaper, int lbtaper, int nx, int ix, int nt, float *trace, int inv);

#endif







