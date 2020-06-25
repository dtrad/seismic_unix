#ifndef RADON_H
#define RADON_H

#include "inversion_par.h"

#ifndef TRACE
#define TRACE fprintf(stderr,"%s @ %u\n",__FILE__,__LINE__)
#endif


#define floatprint(expr) fprintf(stderr,#expr " = %f\n", expr)
#define expfprint(expr) fprintf(stderr,#expr " = %g\n", expr)
#define intprint(expr) fprintf(stderr,#expr " = %d\n", expr)

#define LOOKFAC	2	/* Look ahead factor for npfaro	  */
#define PFA_MAX	720720	/* Largest allowed nfft	          */

void hrrti(float **d,float *pos,float dt,float **m,float *q, float fmax, 
	   int nt, int  nh, int nq, int rtmethod);

void radonfreqband(float *pos, int nh, float **data, float *t, int nt, float dt,
		   float **model, float *q, int nq, float dq, float fmax, float **Wd,
		   int norm, int testadj, float quantil, int rtmethod);

float wpcgnr(void (*oper)  (float *model, float *data, complex ***L,int nt, int nh, 
			    int nq, float dt,complex **m2, complex **d2, int nf0, 
			    float fmax, int adj),
	     float *x, float *b,  complex ***L, int nt, int nh, 
	     int nq, float dt, complex **m2, complex **d2, int nf0, float fmax, 
	     float *Wd, float *Wm, inv_par inv);

void radonf(float *model, float *data, complex ***L,int nt,int nh,int nq,
	    float dt,complex **m2, complex **d2, int nf0,float fmax, int adj);

void interpovv (int nt, int ncdp, float *cdp, float **ovv, float cdpt, float *ovvt);

void plotgather(float **d, int nh, int nt, float dt);

void plotgather(float **d, int nh, int nt, float dt, const char *s);

void save_gather(float **d, int nh, int nt, float dt, char* name);

int taper(float **data, int nt, int nh, int ntaper,int flag);

void nmo(float *d,float *m,float *t,float h,float *vel,int adj,int nt,float dt,
	 float smute);

void getvelocities(float dt, int nt, int ncdp, float *cdp, float **ovv);

void modelweight(float *m, int nx, int norm, float eps1, float *Wm);

#endif


