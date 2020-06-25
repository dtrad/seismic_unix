#ifndef KMIG_H
#define KMIG_H

#ifndef TRACE
#define TRACE fprintf(stderr,"%s @ %u\n",__FILE__,__LINE__)
#endif

#define floatprint(expr) fprintf(stderr,#expr " = %f\n", expr)
#define expfprint(expr) fprintf(stderr,#expr " = %g\n", expr)
#define intprint(expr) fprintf(stderr,#expr " = %d\n", expr)

#include "su.h"
#include "segy.h"
#include "dan.h"
#include "header.h"
#include <signal.h>

void rjwfilter(float **d,int nt,int nh, float dt);
void filt(float *trace,int nt,float dt,float fmax,int ls,int m,float *trf);
void interpovv(int nt, int ncdp, float *cdp, float **ovv, float cdpt, 
	       float *ovvt);
void kmig3(float *d, float cdp, float h,float **m, float *t, float *x, float **vel, int nt, int nx, float dt, float cdpspace, float aper);
void kmig3(float *d, float cdp, float h,float **m, float *t, float *x, float **vel, int nt, int nx, float dt); 

void getvelocities(float dt, int nt, int ncdp, float *cdp, float **ovv);


#endif
