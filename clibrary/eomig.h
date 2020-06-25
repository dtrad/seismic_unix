#ifndef EOMIG_H
#define EOMIG_H

/* Maximum value an `unsigned short int' can hold.  (Minimum is 0).  */
#ifndef USHRT_MAX
#define USHRT_MAX 65535
#endif

#ifndef TRACE
#define TRACE fprintf(stderr,"%s @ %u\n",__FILE__,__LINE__)
#endif

#include"su.h"
#include"cwp.h"
#include"segy.h"
#include"dan.h"
#include "inversion_par.h"

//////////////////////////////////////////////////////////////////
// Function Prototypes
// ls migration for small data
void eomigls_sparse(float *m,float *mm,float *t, float *h, float *q, 
		    float *vel,float **data, int nt, int nh, int nq, 
		    float fmax,float step, float eps, float eps1, float eps2, 
		    int itercg, int iter_end, int smooth, int testadj, 
		    int buffer, float pervmin,float dperv, int norm, 
		    float parmute, unsigned short mute, float t0, 
		    float *Wm,float *Wd,int centralq);

void eomgetvel(float *m,float *mm,float *t, float *h, float *q, float *vel,
	       float **data, int nt, int nh, int nq, float fmax, 
	       float step, float eps, float eps1, float eps2, int itercg, 
	       int iter_end, int smooth, int testadj, int buffer, 
	       float pervmin, float dperv, int norm, int nitervel, float t0,
	       float *Wm, float *Wd);


// ls migration for large data
void eomigls_sparse_large(float *m,float *mm,float *t, float *h, float *q, 
		    float *vel,float **data, int nt, int nh, int nq, 
		    float fmax, float step, float eps, float eps1, 
		    float eps2, int itercg, int iter_end, int smooth, 
		    int testadj, int buffer, float pervmin, float dperv,
		    int norm, float parmute, unsigned short mute, 
		    float t0, float *Wm, float *Wd, int centralq);

void eomigls_sparse_sinc(float *m,float *mm,float *t, float *h, float *q, 
		    float *vel,float **data, int nt, int nh, int nq, 
		    float fmax,float step, float eps, float eps1, float eps2, 
		    int itercg, int iter_end, int smooth, int testadj, 
		    int buffer, float pervmin,float dperv, int norm, 
		    float parmute, unsigned short mute, float t0, 
		    float *Wm,float *Wd,int centralq);



void eomgetvel_large(float *m,float *mm,float *t, float *h, float *q, 
		     float *vel,float **data, int nt, int nh, int nq, 
		     float fmax, float step, float eps, float eps1, float eps2,
		     int itercg, 
		     int iter_end, int smooth, int testadj, int buffer, 
		     float pervmin, float dperv, int norm, int nitervel, 
		     float t0, float *Wm, float *Wd);

void eomgetvel(float *m,float *mm,float *t, float *h, float *q, float *vel,
	       float **data, int nt, int nh, int nq, float fmax, 
	       float step, float eps, float eps1, float eps2, int itercg, 
	       int iter_end, int smooth, int testadj, int buffer, 
	       float pervmin, float dperv, int norm, int nitervel);

void eomig1(float **csp, float *he, float *m, float *t,float *velint,float t0);
void eomig1(float **csp, float *he,float *m, float *t, float *vel,float t0, int nt, int nhcsp);
void rjwfilter(float **d,int nt,int nh, float dt);

void filt(float *trace,int nt,float dt,float fmax,int ls,int m,float *trf);

void equiv_offset(segy tr,float **csp,unsigned short **cspfold,float *t,
		  float xcsp,float *velint,
		  float dt,float *he, float dhe, int nt, int nhcsp, float *tc,
		  float *hec,int precise,float beta,float cdpspace,float scalefold);

void equiv_offset_1side(segy tr,float **csp,unsigned short **cspfold,
			float *t,float xcsp,float *velint,
			float dt,float *he, float dhe, int nt, int nhmax, 
			float *tc, float *hec, int precise, float beta, 
			float cdpspace,float scalefold);
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

void equiv_offset_test(segy tr,float **csp,unsigned short **cspfold,
		       float *t,float xcsp,float *velint,float dt,float *he, 
		       float dhe, int nt, int nhmax,float beta,float cdpspace,
		       float scalefold);

// Function for reading and inetrpolate velocities.

void interpovv (int nt, int ncdp, float *cdp, float **ovv, float cdpt, float *ovvt);

void interpovv(int nt, int ncdp, float *cdp, float **ovv, float **oa1, 
	       float **oa2, float cdpt, float *ovvt, float *oa1t, float *oa2t);

void getvelocities(float dt, int nt, int ncdp, float *cdp, float **ovv);

void getvelocities(float dt, int nt, int ncdp, float *cdp, float **ovv, float **oa1, float **oa2);

//////////////////////////////////////////////////////////////////////

void eomigls(float *t, float *vel, float *h, float *m,float **data,
	   float fmax, int nvel, float step, float eps1, float eps2, 
	   int itercg, int iter_end);

//////////////////////////////////////////////////////////////////////
// Prototypes
void irreg_velaxis(int nq, int nt, float pervmin, float dperv, float *t, float *vel,float *q, float **vgrid);
 
int count_nsparse(float *t, float *h, float *q, float **vel, 
		  int nt, int nh, int nq, float t0);

void build_index(float *t, float *h, float **vel, int nt, int nh, 
		 int nq,unsigned short **index, int nonzero, float t0);

void build_index_inv(float *t, float *h, float **vel, int nt, int nh, 
		 int nq,unsigned short **index, int nonzero, float t0);

void build_index(float *t, float *h, float **vel, int nt, int nh, 
		 int nq,unsigned int **index, int nonzero, float t0);

void build_index_inv(float *t, float *h, float **vel, int nt, int nh, 
		 int nq,unsigned int **index, int nonzero, float t0);

void build_index(float *t, float *h, float **vel, int nt, int nh, int nq,float **index, int ns, float t0 );


void build_index_inv(float *t, float *h, float **vel, int nt, int nh, int nq,float **index, int ns, float t0);

void get_vel(float *m, float *t, float *q, float *vel, float *vel2, 
	     int nt, int nq, int buffer); 

void eomig(float *m, float *d, unsigned short **index, int adj, int nt, int nh, int nq, int nsparse);

void eomig(float *m, float *d, unsigned int **index, int adj, int nt, int nh, int nq, int nsparse);

void eomig(float *m, float *d, float **index, int adj, int nt, int nh, int nq, int nsparse);

void eomig (float *,float *,float *,float *,float *,float *,int ,int ,
	    int,int);

float testadjop(void (*oper) (float *,float *,float *,float *,float *,
			      float *,int ,int ,int,int), 
		float *t, float *h, float *q, float *vel,int nt,int nh,int nq);

void wtcgls(void (*oper) (float *,float *,float *,float *,float *,float *,int ,int ,int,int), int nt, int nh, int np, float *t, float *h, float *p, float *x,float *b,float *Wm, float *Wd, float *vel,float tol, float step, float eps1, float eps2, int itercg);

float dot(int n, float *a, float *b);

void filt(float *trace,int nt,float dt,float fmax,int ls,int m,float *trf);

void eomig1(float **csp, float *he,float *m, float *t, float *vel);

void modelweight(float *m, int nx, int norm, float eps1, float *Wm);

void plotgather(float **d, int nh, int nt, float dt);
void plotgather(float *d, int n1, int n2, const char *s);
void plotgather_pipe(float *d, int n2, int n1, const char *s);
void plotgather_pipe(float **d, int n2, int n1, const char *s);
void save_vector(float *d, int n, char* name);

void radonfreq(float *h, int nh, float **data ,float *t, int nt, float dt, float *vel, float eps1,float eps2, float eps, float **m, float *q, int nq, int itercg, int iter_end, float step, int lsmethod, float factor,float ssmute,float nmofactor, unsigned short mute, float parmute, float quantil);

void radonfreqint(float *h, int nh, float **data ,float *t, int nt, float dt, float *vel, float eps1,float eps2, float eps, float **model, float *q, int nq, int itercg, int iter_end, float step, int lsmethod, float factor, float smute,float nmofactor,unsigned short  mute, float parmute, float *h2, int nh, float quantil);

void nmo(float *d,float *m,float *t,float h,float *vel,int adj,int nt,float dt,float smute);

int taper(float **data, int nt, int nh, int ntaper,int flag);

void doublemute(float **data, float **model, int nq, int nh, int nt, float *q, float *h, float *t, float parmute1, float parmute2, int mute1, int mute2, float fmax, int rtmethod, float depth, char *solver);

#endif













