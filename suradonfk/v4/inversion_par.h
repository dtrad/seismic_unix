#ifndef INVERSION_PAR_H
#define INVERSION_PAR_H

#define MAXITER 100
/* TYPEDEFS */
typedef struct {	/* Parameters for inversion  */

  int itercg;	/* Number maximum of iterations in CG*/

  int iter_end;	/* Number maximum of External Iterations */

  float eps1;	/* noise covariance    */

  float eps2 ;	/* model covariance */

  float eps;	/* small number for tolerance */

  float step;	/* step length factor  */

  int norm;     /* norm to use in model weights; */

  int restart;  /* always set to 1 for now */

  int taperflag; /* If set applies taper to 5 outer traces */

  int mute;      /* if set applies mute */

  float parmute; /* if mute is set it controls the length of muting */

  float J[MAXITER];
  
} inv_par;


typedef struct {	/* Parameters for inversion  */

  int taperflag; /* If set applies taper to 5 outer traces */

  int mute;      /* if set applies mute */

  float parmute; /* if mute is set it controls the length of muting */

  int testadj;   /* if set test adjoint */

} options_par;

/* FUNCTION PROTOTYPES */
void radint_sparse4(float *t, float *q, float *h, float *h2, float **m,float **d,float **dint,int nt, int nh, int nq, int nh2, float dt, float *vel, float dperv, float pervmin, float t0, options_par opt, inv_par inv,int flagvel, int centralq, int filtout);

float wpcgnr(void (*oper) (float *,float *,unsigned int **, int, int ,int ,int,int), int nt, int nh, int np, int nsparse, float *x,float *b,float *Wd, float *M, unsigned int **index, inv_par inv);

float wpcgnr(void (*oper) (float *,float *,unsigned int **, int, int ,int ,int,int, float *wavelet, int nw), int nt, int nh, int np, int nsparse, float *x,float *b,float *Wd, float *M, unsigned int **index, inv_par inv, float *wavelet, int nw);

float wpcgnr(void (*oper)(float *,float *,float *,float *,float *,float **,int ,int ,int,int),int nt, int nh, int np, float *t, float *h, float *p, float *x,float *b,float *Wd, float *M, float **vel, inv_par inv);


#endif
























