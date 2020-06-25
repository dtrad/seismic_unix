#include <math.h>
#include "su.h"
float dot(int n, float *a, float *b);
void radonop(float *m, float *t, float *h, float *q, float *d);
void radonopi(float *m, float *t, float *h, float *q, float *d);


float testadjrad(float *t, float *h, float *q)
{
  float *dr1;
  float *mr1;
  float *dr2;
  float *mr2;
  float dp1;
  float dp2;
  int it;
  extern int ny;
  extern int nx;
  


  if ((dr1=alloc1float(ny))==NULL)
    fprintf(stderr,"***Sorry, space for vrand1 could not be allocated\n");
  if ((mr1=alloc1float(nx))==NULL)
    fprintf(stderr,"***Sorry, space for vrand2 could not be allocated\n");
  if ((dr2=alloc1float(ny))==NULL)
    fprintf(stderr,"***Sorry, space for vrand1 could not be allocated\n");
  if ((mr2=alloc1float(nx))==NULL)
    fprintf(stderr,"***Sorry, space for vrand2 could not be allocated\n");
 
  for (it=0;it<ny;it++) dr1[it]=frannor();
  for (it=0;it<nx;it++) mr1[it]=frannor();
  
  radonopi(mr2,t,h,q,dr1);

  radonop(mr1,t,h,q,dr2); 

  dp1=dot(ny,dr1,dr2);
  dp2=dot(nx,mr1,mr2);

  if (dp2!=0)

  return(dp1/dp2);
  
  free1float(mr2);
  free1float(dr2);
  free1float(mr1);
  free1float(dr1);

}



