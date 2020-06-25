#ifndef RADONWAVELET_H
#define RADONWAVELET_H

#define floatprint(expr) fprintf(stderr,#expr " = %f\n", expr)
#define expfprint(expr) fprintf(stderr,#expr " = %g\n", expr)
#define intprint(expr) fprintf(stderr,#expr " = %d\n", expr)

#ifndef TRACE
#define TRACE fprintf(stderr,"%s @ %u\n",__FILE__,__LINE__)
#endif

#include"su.h"
#include"dan.h"
#include "inversion_par.h"
#include "nrutil.h"
#include "nr.h"

void dwt(float **data, int nh, int nt, int sign);
void radonwavelet0(float **data, float *t, int nt, float *h, int nh, float **model, 
		   float *q, int nq, int wavn, float quantil, int scalefilt);

void mrgathers(float **data, int nt, int nh, int scale,int *itd, float **datatemp);

void mrcoefficients(float **data, int nt, int nh, int scale, float quantil, int scalefilt);

void radonwindow(float **data, float *t, int nt, float *h, int nh, float **model, 
		 float *q, int nq, int scale, float quantil, int scalefilt);
float dyad_timeaxis(int nt, int ntd, float *td, float dt, float t0);

int dyad(int nt, int *id);

#endif













