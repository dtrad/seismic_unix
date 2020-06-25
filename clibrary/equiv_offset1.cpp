/* Copyright (c) University of British Columbia, 1999.*/
/* All rights reserved.                       */
/* EQUIV_OFFSET:  $Date: October 1999  */
/* Precompute times for changing offset bin 
   hec[i] is the equivalent offset at bin i
   ihec   is the corresponding index
   tc[i]  is the time at which euivalent offset bin changes
   itc is the corresponding offset
   hec2 is hec[i]*hec[i]
   xm2  is xm*xm
   h2   is h*h
   crossterm is 4*xm2*h2/(vel*vel);
*/

#include "su.h"
#include "segy.h"
#include <math.h>

void equiv_offset(segy tr,float **csp,float *t,float xcsp,float *velint,float dt,float *he, float dhe, int nt, int nhmax, float *tc, float *hec, int precise, float beta, float cdpspace)
{
  //float  *hec;           /* equivalent offset for every trace */
  //float  *tc;            /* time when equivalent offset bin increases */

  int   ihec;   /* equivalent offset index */
  int   itc;    /* corresponding index for tc[i] */   
  float hec2;   /* he*he  */
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
  float tc0;
  float vel2;
  int itc2;
  float sinbeta=sin(beta);
  int sign;
  int ihec0; // Index that corresponds to he=0
  int verbose2=0; 
  float gx;
  

  if (tr.offset > 0) sign=1;
  else sign=-1;
  dhe=sign*dhe;
  if (verbose2) fprintf(stderr,"dhe=%f,cdpspace=%f\n",dhe,cdpspace);

  for (i=0;i<nhmax;i++){
    tc[i]=0;
    hec[i]=0;
  }

  h=fabs((float) tr.offset);
  h/=2;  // halfoffset
  
  //gx=tr.sx + tr.offset;
  //cdp=tr.sx+tr.offset/2.;
  cdp=tr.cdp;
  xm=fabs(cdp-xcsp);
  xm*=cdpspace;
  xm2=xm*xm;
  h2=h*h;
  hecmax=sign*sqrt(xm2+h2);

  /* To eliminate particular offsets or midpoints use */
  //if (h==0) return; 
  
  if (xm != 0){
    i=0;
    tc[i]=2*xm/(sinbeta*velint[0]);  // First time to use for a ray coming horizontally
    itc=int (tc[i]/dt);if (itc<nt) vel1=velint[itc];else vel1=velint[nt-1];
    hec[i]=xm2+h2-4*xm2*h2/((tc[i]*vel1)*(tc[i]*vel1));

    if (hec[i] > 0) hec[i]=sign*sqrt(hec[i]);
    else hec[i]=0;

    

 
    //fprintf(stderr,"cdp=%f,xcsp=%f\n",cdp,xcsp);    
    do{
      i++;
      hec[i]=hec[i-1]+dhe;
      //fprintf(stderr,"hec[%d]=%f\n",i,hec[i]);
      hec2=hec[i]*hec[i];
      // xm2+h2 would be enough for horizontal layers
      crossterm=xm2+h2-hec2;
      if (crossterm>0) {
	crossterm=sqrt(crossterm);
	if (precise){
	  /* First estimate of time uses previous velocity */
	  for (k=0;k<3;k++){
	    tc[i]=2*xm*h/(vel1*crossterm);
	    tc0=tc[i]*tc[i]/4-hec2/(vel1*vel1);
	    if (tc0 > 0) tc0=2*sqrt(tc0);
	    else tc0=0; 
	    itc=int (tc0/dt);
	    if (itc<nt) vel1=velint[itc];
	    else vel1=velint[nt-1];
	  }
	}
	else {
	  tc[i]=2*xm*h/(vel1*crossterm);
	  itc=(int) (tc[i]/dt);
          if (itc<nt) vel1=velint[itc];
	  else vel1=velint[nt-1];
	}
      }
      else tc[i]=t[nt-2];   
    }while((tc[i] < t[nt-2]) && (i<nhmax) && (fabs(hec[i])<fabs(hecmax)));
    nhec=i++;  
  }
  else if (xm==0) {  // xm=0 is a CMP, hence he=h; h=0 is a zero offset trace
    tc[0]=2*(xm+h)/(sinbeta*velint[0]);
    hec[0]=sign*sqrt(xm2+h2);
    nhec=1;
  }
  tc[nhec]=t[nt-2];
  hec[nhec]=hecmax;


  it=int (tc[0]/dt);
  if (sign==1) ihec=(int) ((hec[0]-he[0])/fabs(dhe));
  else  ihec=(int) ((hec[0]-he[0])/fabs(dhe))+1;

  if (verbose2)
  for (i=0;i<=nhec;i++){
    j=ihec+sign*i; 
    if (j >=0 && j < nhmax) fprintf(stderr,"tc[%d]=%f,hec[%d]=%f,he[%d]=%f,dhe=%f\n",i,tc[i],i,hec[i],j,he[j],dhe);
  }


  if (verbose2) if (ihec >=0 && ihec < nhmax) 
    fprintf(stderr,"hec[0]=%f,he[%d]=%f\n",hec[0],ihec,he[ihec]);
  for (i=0;i<nhec;i++){
    itc=(int) (tc[i+1]/dt);
    if (ihec >= 0 && ihec < nhmax ){    
      do{
	it++;
	if (verbose2) fprintf(stderr,"he[%d]=%f,it=%d\n",ihec,he[ihec],it);
       	csp[it][ihec]+=tr.data[it];	    
      }while(it<=itc);
      ihec=ihec+sign*1;
    }
    else break;
  }
  return;
}  

  












