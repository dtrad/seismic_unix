/* Copyright (c) University of British Columbia, 1999.*/
/* All rights reserved.                       */
/* EQUIV_OFFSET:  $Date: October 1999  */
/* Compute equivalent offset time by time to test other methods
   This method must be precise but slow, so must be used only for test 
   hei  is the equivalent offset at bin i
   ihei   is the corresponding index
   t0   is the time at zero offset
   it0  is the corresponding offset
   he2  is he[i]*he[i]
   xm2  is xm*xm
   h2   is h*h
   crossterm is 4*xm2*h2/(vel*vel);
*/

#include "su.h"
#include "segy.h"

void equiv_offset_test(segy tr,float **csp,float *t,float xcsp,float *velint,float dt,float *he, float dhe, int nt, int nhmax, float beta, float cdpspace)
{
  int it00;
  float t0;     /* time at zero offset */
  int it0;      /* Index */
  float hei;    /* equivalent offset for every trace */
  int   ihe ;   /* equivalent offset index */   
  float he2;    /* he*he  */
  int   nhec;   /* Number of offset bins at every trace */
  float cdp;    /* midpoint coordinate */
  float xm;     /* midpoint distance from trace to scp = cdp_tr - x_scp */
  float xm2;             // xm*xm
  float h;               /* halfoffset */
  float h2;              // h*h
  float crossterm;       // 4*xm^2*h*2/vel
  float hecmax;          /* Maximum possible eq.  offset for this trace */
  float vel1;    // Temporal variable for velocity
  int i,j,k;
  register int it;
  float t00;
  float sinbeta=sin(beta);
  int sign;
  int verbose2=0; 
  float gx;

  if (tr.offset > 0) sign=1;
  else sign=-1;
  if (verbose2) fprintf(stderr,"dhe=%f\n",dhe);
  
  h=fabs((float) tr.offset);
  h/=2;  // halfoffset
  //cdp=tr.sx+tr.offset/2.;

  cdp=tr.cdp;
  cdp*=cdpspace;
  xm=fabs(cdp-xcsp);

  xm2=xm*xm;
  h2=h*h;

  // First time to use for a ray coming horizontally
  t00=2*xm/(sinbeta*velint[0]);  
  it00=(int) (t00/dt);
  if (it00>=nt) it00=nt-1;
  vel1=velint[0];
  hei=xm2+h2-4*xm2*h2/((t[it00]*vel1)*(t[it00]*vel1));
  if (hei > 0) hei=sign*sqrt(hei);
  else hei=0;
  ihe=(int) ((hei-he[0])/dhe);
  if (sign==-1) ihe+=1;
  if (ihe>=0 &&  ihe < nhmax){
    if (it00 < 10) it00=10;
    for (it=it00;it< nt;it++){
      // For every t[it] estimate vel0 and equiv offset
      // The velocity must be estimated at the apex, i.e., t0
      t0=t[it]*t[it]/4-(hei*hei)/(vel1*vel1);
      if (t0 > 0) t0=2*sqrt(t0);
      else t0=0; 
      it0=int (t0/dt);
      if (it0<nt) vel1=velint[it0];
      else vel1=velint[nt-1];
      // xm2+h2 would be enough for horizontal layers
      // the crossterm corrects for xm != 0
      hei=xm2+h2-4*xm2*h2/((t[it]*vel1)*(t[it]*vel1));
      if (hei > 0) hei=sign*sqrt(hei);
      else hei=0;
      ihe=(int) ((hei-he[0])/dhe);
      if (sign==-1) ihe+=1;
      if (ihe >=0 && ihe < nhmax) csp[it][ihe]+=tr.data[it];      
    }
  }
  
  return;
}  

  


