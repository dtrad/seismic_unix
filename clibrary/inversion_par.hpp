#include <iostream.h>
#ifndef INVERSION_PAR_H
#define INVERSION_PAR_H


/* TYPEDEFS */
class inversion {	/* Parameters for inversion  */
public:
  inversion();
  ~inversion();
  void SetItercg(int itercg);
  int GetItercg() const;
  int iter_end;	/* Number maximum of External Iterations */
  float eps1;	/* noise covariance    */
  float eps2;	/* model covariance */
  float eps;	/* small number for tolerance */
  float step;	/* step length factor  */
  int norm;     /* norm to use in model weights; */
  int restart;  /* always set to 1 for now */
  int taperflag; /* If set applies taper to 5 outer traces */
  int mute;      /* if set applies mute */
  float parmute; /* if mute is set it controls the length of muting */
private:
  int itercg;	/* Number maximum of iterations in CG*/  

} ;

inversion::inversion()
{ 
  fprintf(stderr," here my first class constructor\n") ; 
}

inversion::~inversion()
{
 fprintf(stderr," here my first class destructor \n") ; 
}

void inversion::SetItercg(int i)
{ 
  itercg=i;
  fprintf(stderr," here my first class constructor\n") ; 
}

int inversion::GetItercg() const
{ 
  fprintf(stderr," here GetItercg\n") ;
  return(itercg);
}


/* FUNCTION PROTOTYPES */

#endif
























