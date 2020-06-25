#ifndef EOMIG_H
#define EOMIG_H


#ifndef TRACE
#define TRACE fprintf(stderr,"%s @ %u\n",__FILE__,__LINE__)
#endif

#include"su.h"
#include"cwp.h"
#include"segy.h"
#include"dan.h"


void eomig1(float **csp, float *he,float *m, float *t, float *vel,float t0, int nt, int nhcsp);
void rjwfilter(float **d,int nt,int nh, float dt);

void filt(float *trace,int nt,float dt,float fmax,int ls,int m,float *trf);

void equiv_offset(segy tr,float **csp,unsigned short **cspfold,float *t,
		  float xcsp,float *velint,
		  float dt,float *he, float dhe, int nt, int nhcsp, float *tc,
		  float *hec,int precise,float beta,float cdpspace,float scalefold);

/* segy arguments give problems, use members of the struture instead */ 
void equiv_offset_1side(float *trace, int offset, int cdp, float **csp, 
			unsigned short **cspfold, float *t,float xcsp,float *velint,
			float dt,float *he, float dhe, int nt, int nhmax, float *tc,
			float *hec, int precise, float beta, float cdpspace, 
			float scalefold);

void equiv_offset_test(float *trace, int offset, int cdp, float **csp,unsigned short **cspfold,
		       float *t,float xcsp,float *velint,float dt,float *he, 
		       float dhe, int nt, int nhmax,float beta,float cdpspace,
		       float scalefold);

// Function for reading and interpolate velocities.

void interpovv (int nt, int ncdp, float *cdp, float **ovv, float cdpt, float *ovvt);

void getvelocities(float dt, int nt, int ncdp, float *cdp, float **ovv);

#endif













