#ifndef PWDF_H
#define PWDF_H

#define floatprint(expr) fprintf(stderr,#expr " = %f\n", expr)
#define expfprint(expr) fprintf(stderr,#expr " = %g\n", expr)
#define intprint(expr) fprintf(stderr,#expr " = %d\n", expr)

#ifndef TRACE
#define TRACE fprintf(stderr,"%s @ %u\n",__FILE__,__LINE__)
#endif

#include"su.h"
#include"dan.h"
#include "inversion_par.h"

void pwdf0(float **data, float *t, int nt, float *h, int nh, float **model, float *q, int nq);
void hpwavelet(float **wav,float sx, float st,float tau,int nt,int nx,float dt,float dx);

#endif













