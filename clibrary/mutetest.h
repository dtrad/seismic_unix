#ifndef MUTETEST_H
#define MUTETEST_H

#define LOOKFAC 2       /* Look ahead factor for npfaro   */
#define PFA_MAX 720720  

#ifndef TRACE
#define TRACE fprintf(stderr,"%s @ %u\n",__FILE__,__LINE__)
#endif

#define floatprint(expr) fprintf(stderr,#expr " = %f\n", expr)
#define expfprint(expr) fprintf(stderr,#expr " = %g\n", expr)
#define intprint(expr) fprintf(stderr,#expr " = %d\n", expr)

#include "su.h"
#include "dan.h"


void radonhyper(TGather &data, TGather &model);
void AtimesB(float **A, float **B, int nh, int nt);


typedef struct {        /* Parameters for inversion  */
  float tmin; 
  float tmax;
  int ihmin;
  int ihmax;
  float slope;
  float threshold;
} mutemask_par;


void mutemask(TGather &data, mutemask_par par);

#endif
