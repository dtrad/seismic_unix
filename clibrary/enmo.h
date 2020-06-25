#include "inversion_par.h"

#ifndef ENMO_H
#define ENMO_H

#ifndef TRACE
#define TRACE fprintf(stderr,"%s @ %u\n",__FILE__,__LINE__)
#endif


/* FUNCTION PROTOTYPES */
void irreg_velocity_grid(int nv, int nt, float pervmin, float dperv, float *t, float *vel,
			 float **vgrid, int centralq, int nreg);

int rad_ellip_sp(float *t,float *p,float *v, float **vgrid,int nt, int np, int nv,  unsigned short **nullspace, float *critangle);

void rad_ellip_sp(float *t,float *p,float *v, float **vgrid,int nt, int np, int nv, unsigned int **index, int nsparse, float *critangle);

void rad_ellip_sp_LI(float *t,float *p,float *v, float **vgrid,int nt, int np, int nv, unsigned int **index, int nsparse, float *critangle);

void rad_ellip_sp_old(float *t,float *p,float *v, float **vgrid,int nt, int np, int nv, unsigned int **index, int nsparse, unsigned short **nullspace, float *critangle);

void radonellip(float *m, float *d, unsigned int **index, int adj, int nt, int nh, int nq, int nsparse);

void radonellip_LI(float *m, float *d, unsigned int **index, int adj, int nt, int nh, int nq, int nsparse);

float testadjop(void (*oper) (float *,float *,float *,float *,float *,float **,int ,int ,int,int), float *t, float *h, float *q, float **vel,int nt, int nh, int nq);

float dot(int n, float *a, float *b);

void enmo1(float *m,float *t,float *p,float *v,float *d, float *vel,int adj, int nt, int np, int nv);

void enmo1(float *m,float *t,float *p,float *v,float *d, float **vgrid,int adj, int nt, int np, int nv);

void rms2intvel(int n,float *t0, float *vs, float *v, float *h);

void myintlin (int nin, float xin[], float yin[], float yinl, float yinr, 
	int nout, float xout[], float yout[]);

void plotgather(float **d, int nh, int nt, float dt, const char *s);

int tapercentre(float **data, int nh, int nt, int ntaper, int ntapertime, int centre);

void create_elliptical_taper(float *t,float *p,float *v, float **vgrid,int nt, int np, int nv, 
			     float **Wd, float *critangle);

void enmo_make_win(float **d, float **m, float *t, float *p, float *v, int nt, int np, int nv, 
		   float **vgrid, float *critangle, float **Wd, int LI, float t0, float tmax, 
		   inv_par inv, int itercgfin);

void enmo_win(float **d, float **m, float *t, float *p, float *v, int nt, int np, int nv, 
	      float **vgrid, float *critangle, float **Wd, int LI, inv_par inv, int itercgfin);

void create_taper(float **Wd, float *t, float *p, int nt, int np, int taper, int ntaper, int ntapertime);

void critical_angle(float *t, int nt, float *tvel, int nvel, float *critang, float *critangle);


/****************************************************************************************/

#endif





