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

void equiv_offset(segy tr,float **csp,float *t,float xcsp,float *velint,float dt,float *he, float dhe, int nt, int nhmax, float *tc, float *hec)
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
  int sign;

  if (tr.offset > 0) sign=1;
  else sign=-1;
  dhe=sign*dhe;
  
  h=fabs((float) tr.offset);

  h/=2;  // halfoffset
  cdp=(float) tr.cdp;
  xm=fabs(cdp-xcsp);

  xm2=xm*xm;
  h2=h*h;
  /* To eliminate particular offsets or midpoints use 
  if (h==0) return; 
  */
  if ((xm!=0)&&(h!=0)){
    i=0;
    tc[i]=2*xm/velint[0];  // First time to use
    itc=int (tc[i]/dt);if (itc<nt) vel1=velint[itc];else vel1=velint[nt-1];
    hec[i]=sign*xm;
    hecmax=sign*sqrt(xm2+h2);
        
    do{
      i++;
      hec[i]=hec[i-1]+dhe;
      hec2=hec[i]*hec[i];
      crossterm=xm2+h2-hec2;
      if (crossterm>0) {
	crossterm=sqrt(crossterm);
	/* First estimate of time uses previous velocity */
	for (k=0;k<3;k++){
	  tc[i]=2*xm*h/(vel1*crossterm);
	  itc=int (tc[i]/dt);
	  if (itc<nt) vel1=velint[itc];
	  else vel1=velint[nt-1];
	  //fprintf(stderr,"itc=%d,k=%d,i=%d,tc[i]=%f\n",itc,k,i,tc[i]);
	}
      }
      else tc[i]=t[nt-2];
      //fprintf(stderr,"Inside do loop, i=%d,nhmax=%d\n",i,nhmax);
      //fprintf(stderr,"Inside do loop, tc[%d]=%f,t[nt-2]=%f\n",i,tc[i],t[nt-2]);
      //fprintf(stderr,"Inside do loop, hec[%d]=%f,hecmax=%f\n",i,hec[i],hecmax);   
    }while((tc[i] < t[nt-2]) && (i<nhmax) && (hec[i]<hecmax));
    nhec=i++;
    //fprintf(stderr,"Inside equiv_offset, nt=%d,nhmax=%d,hecmax=%f\n",nt,nhmax,hecmax);    
  }
  else{
    tc[0]=2*(xm+h)/velint[0];
    hec[0]=sign*sqrt(xm2+h2);
    nhec=1;
  }
  tc[nhec]=t[nt-2];
  hec[nhec]=hecmax;
          
  it=int (tc[0]/dt);
  ihec=(int) fabs((hec[0]/dhe));
  for (i=0;i<nhec;i++){
    itc=(int) (tc[i+1]/dt);    
    if (ihec < nhmax){    
      do{
	it++;
	csp[it][ihec]+=tr.data[it];	    
      }while(it<=itc);
      ihec=ihec+1;
    }
    else break;
  }
return;
}  

  


