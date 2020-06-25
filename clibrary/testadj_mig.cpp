/*
This function tests that two given operators are adjoint each other.
For an opeartor L defined as  d=Lm , the adjoint LH, such that  ma=LHd 
fullfill the property
(Lm,d)=(m,LHd)
In this test m, and d are ramdom vectors and L and LH are applied to verify
this property. 
Ref: Claerbout, Processing vs Inversion, 1987, Blackwell Pub. chapter 5
Daniel Trad- December 1999

 */
#include <math.h>
#include "su.h"

float dot(int n, float *a, float *b);
void migration(float *m, float *t, float *h, float *vel, float *d); 
void modelling(float *m, float *t, float *h, float *vel, float *d);


float testadj_mig(float *t, float *h, float *vel)
{
  float *dr1;
  float *mr1;
  float *dr2;
  float *mr2;
  float dp1;
  float dp2;

  extern int nt;
  extern int nhcsp;
  int ny=nt*nhcsp;
  int it;
  int ih;

  if ((dr1=alloc1float(ny))==NULL)
    fprintf(stderr,"***Sorry, space for vrand1 could not be allocated\n");
  if ((mr1=alloc1float(nt))==NULL)
    fprintf(stderr,"***Sorry, space for vrand2 could not be allocated\n");
  if ((dr2=alloc1float(ny))==NULL)
    fprintf(stderr,"***Sorry, space for vrand1 could not be allocated\n");
  if ((mr2=alloc1float(nt))==NULL)
    fprintf(stderr,"***Sorry, space for vrand2 could not be allocated\n");
 
  for (it=0;it<ny;it++) dr1[it]=frannor();
  for (ih=0;ih<nhcsp;ih++) for (it=0;it<100;it++) dr1[ih*nt+it]=0;
  for (it=0;it<nt;it++) mr1[it]=frannor();
  for (it=nt-100;it<nt;it++) mr1[it]=0;
  migration(mr2,t,h,vel,dr1);

  modelling(mr1,t,h,vel,dr2); 

  dp1=dot(ny,dr1,dr2);

  dp2=dot(nt,mr1,mr2);
  if (dp2!=0){
    fprintf(stderr,"testadj=%e\n",dp1/dp2);
    return(dp1/dp2);
  }
  else 
    return(0.);
  

  free1float(mr2);
  free1float(dr2);
  free1float(mr1);
  free1float(dr1);

}










