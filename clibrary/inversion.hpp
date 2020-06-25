#ifndef INVERSION_PAR_H
#define INVERSION_PAR_H


/* TYPEDEFS */
class inv_par {	/* Parameters for inversion  */
public:
  int itercg=20 ;	/* Number maximum of iterations in CG*/
  int iter_end=1;	/* Number maximum of External Iterations */
  float eps1=1e-3;	/* noise covariance    */
  float eps2=1e-3;	/* model covariance */
  float eps=1e-7;	/* small number for tolerance */
  float step=0.95;	/* step length factor  */
  int norm=0;     /* norm to use in model weights; */
  int restart=1;  /* always set to 1 for now */
  int taperflag=0; /* If set applies taper to 5 outer traces */
  int mute=0;      /* if set applies mute */
  float parmute=0; /* if mute is set it controls the length of muting */
  void speak(){ fprintf(stderr," here my first class \n") ; return; }
} ;


/* FUNCTION PROTOTYPES */

#endif
























