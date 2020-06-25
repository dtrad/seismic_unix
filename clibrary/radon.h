#ifndef RADON_H
#define RADON_H

#define floatprint(expr) fprintf(stderr,#expr " = %f\n", expr)
#define expfprint(expr) fprintf(stderr,#expr " = %g\n", expr)
#define intprint(expr) fprintf(stderr,#expr " = %d\n", expr)

void radontoepf(float **data, float *h, int nh,float *t, int nt, float dt, float **model, float *q, int nq, float *vel, float eps1, float factor, float smute, float nmofactor, int rtmethod, float depth, float fmax);


void toepradon(float *pos, float **d, float dt,float **m, float *q, float dq,
	       float eps1, float eps2, float eps, float fmax, int nt, int nh,
	       int nq, int rtmethod);

void radoncgfft(float **data, float *h, int nh,float *t, int nt, float dt, float **model, float *q, int nq, float *vel, int itercg, int iter_end, float step, float eps2, float eps1, float quantil, int norm, float factor, float smute, float nmofactor, int rtmethod, float depth);

void radoncg_levinson(float *pos, int nh, float **data, float *t,int nt,float dt,float **model, float *q, int nq, float dq,float eps1, float eps2, float eps, float fmax, float *Wd, int itercg, int iter_end, int norm, float step, int testadj, float quantil, int rtmethod, float depth);

void radonwtcgls0(float **data, float *h, int nh,float *t, int nt, float dt, float **model, float *q, int nq, float *vel, int itercg, int iter_end, float step, float eps2, float eps1, float quantil, int norm, float factor, float smute, float nmofactor, int rtmethod, float depth, float fmax);

void wtcgls0(float *pos, int nh, float **data, float *t, int nt, float dt,float **model, float *q, int nq, float dq, float eps1, float eps2, float eps, float fmax, float *Wd, int itercg, int iter_end,int norm,float step, int testadj, int rtmethod);

void hrrti(float **d,float *pos,float dt,float **m,float *q, float fmax, 
	   int nt, int  nh, int nq, int rtmethod);

void radonhybrid(float **data, float *h, int nh,float *t, int nt, float dt, float **model, float *q, int nq, float *vel, int itercg, int iter_end, float step, float eps2, float eps1, float quantil, int norm, float qmin1, float qmin2,  int nq1, int nq2, float factor1, float factor2, float smute, float nmofactor, int rtmethod1, int rtmethod2, float depth1, float depth2, float fmax1, float fmax2,float *ffilter, float *amps);

void radonwtcgls_beam(float *h, int nh, float **data, float *t, int nt, float dt,float **model, float *q, int nq, float dq, float eps1, float eps2, float eps, float fmax, float *Wd, int itercg, int iter_end,int norm,float step, int testadj, float quantil, int nq1, int nq2, int rtmethod1, int rtmethod2, float depth1, float depth2, float qmin1, float qmin2, float factor1, float factor2);

void radonwtcgls_beam2(float *h, int nh, float **data, float *t, int nt, float dt,float **model, float *q, int nq, float dq, float eps1, float eps2, float eps, float fmax, float *Wd, int itercg, int iter_end,int norm,float step, int testadj, float quantil, int nq1, int nq2, int rtmethod1, int rtmethod2, float depth1, float depth2, float qmin1, float qmin2, float factor1, float factor2, float fmax1, float fmax2);

void radoninv_beam(float *h, int nh, float **data, float *t, int nt, float dt,
		   float **model, float *q,  int nq, int nq1, int nq2, float qmin1, 
		   float qmin2,  int rtmethod1, int rtmethod2, float depth1, 
		   float depth2, float factor1, float factor2, float fmax1, float fmax2,
  		   float *f, float *amps);

void radonwtcgls0_tfd(float **data, float *h, int nh,float *t, int nt, float dt, float **model, float *q, int nq, float *vel, int itercg, int iter_end, float step, float eps2, float eps1, float quantil, int norm, float factor, float smute, float nmofactor, int rtmethod, float depth);

void radonwtcgls_tfd(float *pos, int nh, float **data, float *t, int nt, float dt,
		     float **model, float *q, int nq, float dq, float eps1, float eps2,
		     float eps, float fmax, float **Wd, int itercg, int iter_end,
		     int norm,float step, int testadj, float quantil, int rtmethod);


void interpovv (int nt, int ncdp, float *cdp, float **ovv, float cdpt, float *ovvt);

void interpovv (int nt, int ncdp, float *cdp, float **ovv, float **oa1, 
	float **oa2, float cdpt, float *ovvt, float *oa1t, float *oa2t);


void plotgather(float **d, int nh, int nt, float dt);

void plotgather(float **d, int nh, int nt, float dt, const char *s);

void save_gather(float **d, int nh, int nt, float dt, char* name);

int taper(float **data, int nt, int nh, int ntaper,int flag);

void nmo(float *d,float *m,float *t,float h,float *vel,int adj,int nt,float dt,
	 float smute);

void getvelocities(float dt, int nt, int ncdp, float *cdp, float **ovv);

void getvelocities(float dt, int nt, int ncdp, float *cdp, float **ovv, float **oa1, float **oa2);



#endif


