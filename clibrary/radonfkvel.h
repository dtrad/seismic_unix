#ifndef RADONFK_H
#define RADONFK_H

#define LOOKFAC	2	/* Look ahead factor for npfaro	  */
#define PFA_MAX	720720	

#ifndef TRACE
#define TRACE fprintf(stderr,"%s @ %u\n",__FILE__,__LINE__)
#endif

#define floatprint(expr) fprintf(stderr,#expr " = %f\n", expr)
#define expfprint(expr) fprintf(stderr,#expr " = %g\n", expr)
#define intprint(expr) fprintf(stderr,#expr " = %d\n", expr)

#include"su.h"
#include"dan.h"
#include "inversion_par.h"

void save_gather(float **d, int nh, int nt, float dt, char* name);

void save2d(float **d, int nh, int nt, float dt, char* name);

void stolt1k (float k, float v, float s, float fmax, int nt, float dt,
		     complex *p, complex *q);
void stolt1kadj (float k, float v, float s, float fmax, int nt, float dt,
			complex *p, complex *q);
int taper(float **data, int nt, int nh, int ntaper,int flag);

void taper (int lxtaper, int lbtaper, 
		   int nx, int ix, int nt, float *trace);

void plotAmpSpec(complex **pk, int nk, int nt, float dt);

void save_gather_freq(complex **d, int nh, float *h, int nt, float dt, int nfft, 
		      char* name);

void stoltzop2(float **data, float **model, int nt, int nh, float *t, 
	       float *h, float vel, int adj);

void stoltzopinv2(float **data, float **model, int nt, int nx, float *t, 
		  float *x, float vel);

float adjteststoltz(int nt, int nh, float *t, float *h, float vel);

void stoltz_wtcgls(float **data, float **model, float **Wd, float *h, int nh,  float *t, int nt,
		   float vel, inv_par inv);

float wpcgnr_mig(void (*oper)  (float **data, float **model, int nt, int nh, float *t,
				      float *h, float vel, int adj),
		       float **x, float **b,int nt, int nh, float *t, float *h,  
		       float vel, float **Wm, float **Wd, inv_par inv);

/*************             Weights               *****************/
void deviations(float *m, int nx, float *d, int ny, int norm, float quantil1, float quantil2, 
		float *sigmam, float *sigmad);

void weights(float *m, int nx, int norm, float sigmam, float *Wm, int iter);

void wdmask(float **Wd, int nt, int nh);

void wmmask(float **Wm, int nt, int nq);
/******************************************************************/
void CequalAxB(float **C,float **A,float **B,int n1,int n2);


void CequalApluskxB(float **C, float **A, float k, float **B, int n1,int n2);

void AequalB(float **A,float **B,int n1, int n2);

void AequalB(complex **A, complex **B,int n1, int n2);

void kaxis(int nx, float vmax, float dt, int nt, float dx, int *pnk, float *pdk);

void save_gather_test(float **d, int nh, int nt, float dt, char* name);

/* PLots */
void xplotgather(float **d, int nh, int nt, float dt, char *s, int num, char *s2);
void ximageplotgather(float **d, int nh, int nt, float dt, char *s, int num, char *s2);

void plot_after_stretch(float **d, int nh, int nt, int dt, char *s1, char *s2);
/* stretching */
void makev (int nmig, float *tmig, float *vmig, float vscale,
	    int nt, float dt, float ft, float **v, float *vmin, float *vmax);
void makeut (float vstolt, float fmax, float *vrms,
	     int nt, float dt, float **ut, int *nu, float *du, float **tu);
void makeu (float vstolt, float *v, int nt, float dt, float *u);

void stretch(float **data, float **datas, int nt, int nu, int nh, float *t, float *tu, float *ut,
	     float dt, float du, int str);

/****************************************************************************************/
 

void plotvector(float *d, int nt, const char *s);

void AtimesB(float **A, float **B, int nh, int nt);

typedef struct {	/* Parameters for inversion  */
  float tmin; 
  float tmax;
  int ihmin;
  int ihmax;
  float slope;
  float threshold;
} mutemask_par;

void mutemask(float **M, float **d, int nh, int nt, float dt, mutemask_par mutepar);

void write_curve_mute(mutemask_par par);

#endif















