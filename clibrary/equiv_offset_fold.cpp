/* Copyright (c) University of British Columbia, 1999.*/
/* All rights reserved.                       */
/* EQUIV_OFFSET:  $Date: October 1999  */
/* Map this trace to the current csp 
   Three options:
   equiv_offset precomputes times, and uses + and - offsets
   equiv_offset_test does all computations, (precise and slow)
   equiv_offset_1side precomputes times, and uses 1 side of offset

   Precompute times for changing offset bin 
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
#define TRACE fprintf(stderr,"%s @ %u\n",__FILE__,__LINE__)

int t2t0(float t,float h,float v, float dt);

void equiv_offset(segy tr,float **csp,unsigned short **cspfold,float *t,float xcsp,float *velint,float dt,float *he, float dhe, int nt, int nhmax, float *tc, float *hec, int precise, float beta, float cdpspace,float scalefold)
{
  //float  *hec;           /* equivalent offset for every trace */
  //float  *tc;            /* time when equivalent offset bin increases */

  int   ihec;   /* equivalent offset index */
  int   itc;    /* corresponding index for tc[i] */   
  int   it0;     /* corresponding index for t0[i] */   
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
  int i,j;
  register int it;
  float sinbeta=sin(beta);
  int sign;
  int verbose2=0; 
  
  // If the trace has negative offset all its samples will map to the
  // negative side of the CSP
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

  // Options to calculate the cdp coordinate.
  // If cdp is set right and is integer just use tr.cdp *cdpspace
  // Otherwise compute from sx and gx the actual coordinate  
  //gx=tr.sx + tr.offset;
  //cdp=tr.sx+tr.offset/2.;
  //cdp=(tr.gx-tr.sx)/2.+tr.sx;

  cdp=tr.cdp;
  xm=fabs(cdp-xcsp);
  //fprintf(stderr,"cdp=%f,xcsp=%f,xm=%f\n",cdp,xcsp,xm);
  xm*=cdpspace; // cdpspace is the physical distance between cdps
  xm2=xm*xm;
  h2=h*h;
  hecmax=sign*sqrt(xm2+h2);

 // Beta ==0 is a code for processing only cdps
  if (beta==0) sinbeta=1;

  /* To eliminate particular offsets or midpoints use */
  //if (h==0) return; 
  //if (xm==0) return;
  // Velocity must be estimated at t
  if ((xm != 0)&&(beta!=0)){
    i=0;
    // First time to use for a ray coming horizontally
    tc[i]=2*MAX(xm,h)/(sinbeta*velint[0]);  
    itc=int (tc[i]/dt+0.5);
    it0=t2t0(t[MIN(itc,nt-1)],MAX(xm,h),velint[MIN(itc,nt-1)],dt);
    vel1=velint[MIN(it0,nt-1)];
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
	tc[i]=2*xm*h/(vel1*crossterm);
	itc=(int) (tc[i]/dt+0.5);
	it0=t2t0(t[MIN(itc,nt-1)],hec[i-1],velint[MIN(itc,nt-1)],dt);
	vel1=velint[MIN(it0,nt-1)];
      }
      else tc[i]=t[nt-2];   
    }while((tc[i] < t[nt-2]) && (i<nhmax) && (fabs(hec[i])<fabs(hecmax)));
    nhec=i++;  
  }
  else{  // xm=0 is a CMP, hence he=h; h=0 is a zero offset trace
    tc[0]=2*h/(sinbeta*velint[0]);
    hec[0]=sign*h;
    nhec=1;
  }
  tc[nhec]=t[nt-2];
  hec[nhec]=hecmax;

  it=int (tc[0]/dt+0.5);
  ihec=(int) ((hec[0]-he[0])/fabs(dhe)+sign*0.5);

  if (verbose2)
  for (i=0;i<=nhec;i++){
    j=ihec+sign*i; 
    if (j >=0 && j < nhmax) fprintf(stderr,"tc[%d]=%f,hec[%d]=%f,he[%d]=%f,dhe=%f\n",i,tc[i],i,hec[i],j,he[j],dhe);
  }


  if (verbose2) if (ihec >=0 && ihec < nhmax) 
    fprintf(stderr,"hec[0]=%f,he[%d]=%f\n",hec[0],ihec,he[ihec]);
  for (i=0;i<nhec;i++){
    itc=(int) (tc[i+1]/dt+0.5);
    if (ihec >= 0 && ihec < nhmax ){    
      do{
	it++;
	//if (verbose2) fprintf(stderr,"he[%d]=%f,it=%d\n",ihec,he[ihec],it);
       	csp[ihec][it]+=tr.data[it];
	cspfold[ihec][it]=cspfold[ihec][it]+(int) fabs(scalefold*tr.data[it]);
	//cspfold[ihec][it]++;
      }while(it<=itc);
      ihec=ihec+sign;
    }
    else break;
  }
  return;
}  




void equiv_offset_1side(segy tr,float **csp, unsigned short **cspfold, float *t,float xcsp,float *velint,float dt,float *he, float dhe, int nt, int nhmax, float *tc, float *hec, int precise, float beta, float cdpspace, float scalefold)
{
  //float  *hec;           /* equivalent offset for every trace */
  //float  *tc;            /* time when equivalent offset bin increases */

  int   ihec;   /* equivalent offset index */
  int   itc;    /* corresponding index for tc[i] */   
  int   it0;     /* corresponding index for t0[i] */   
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
  int i;
  register int it;

  //fprintf(stderr,"beta=%f\n",beta);
  float sinbeta=sin(beta);
  //fprintf(stderr,"sin(beta)=%f\n",sinbeta);

  for (i=0;i<nhmax;i++){
    tc[i]=0;
    hec[i]=0;
  }

  h=fabs((float) tr.offset);
  h/=2;  // halfoffset
  //cdp=tr.sx+tr.offset/2.;
  
  cdp=tr.cdp;

  xm=fabs(cdp-xcsp);
  // test for reconstructing the cdp from csp.
  // if (!xm) return;
  xm*=cdpspace;

  xm2=xm*xm;
  h2=h*h;
  hecmax=sqrt(xm2+h2);
  /* To eliminate particular offsets or midpoints use */
  //if (h==0) return; 
  //if (xm==0) return;

 // Beta ==0 is a code for processing only cdps
  if (beta==0) sinbeta=1;

  if ( xm && beta ){
    i=0;
    tc[i]=2*MAX(xm,h)/(sinbeta*velint[0]);  // First time to use for a ray coming horizontally
    //if (i==0 && tc[0]<0.1) fprintf(stderr,"xm=%f,h=%f\n",xm,h);
    itc=int (tc[i]/dt);
    it0=t2t0(t[MIN(itc,nt-1)],MAX(xm,h),velint[MIN(itc,nt-1)],dt);
    vel1=velint[MIN(it0,nt-1)];
    hec[i]=xm2+h2-4*xm2*h2/((tc[i]*vel1)*(tc[i]*vel1));

    if (hec[i] > 0) hec[i]=sqrt(hec[i]);
    else hec[i]=0;
    
    hecmax=sqrt(xm2+h2);
 
    //fprintf(stderr,"cdp=%f,xcsp=%f\n",cdp,xcsp);    
    do{
      i++;
      hec[i]=hec[i-1]+dhe;
      hec2=hec[i]*hec[i];
      // xm2+h2 would be enough for horizontal layers
      crossterm=xm2+h2-hec2;
      if (crossterm>0) {
	crossterm=sqrt(crossterm);
	tc[i]=2*xm*h/(vel1*crossterm);
	itc=(int) (tc[i]/dt+0.5);
	it0=t2t0(t[MIN(itc,nt-1)],hec[i-1],velint[MIN(itc,nt-1)],dt);
	vel1=velint[MIN(it0,nt-1)];
      }
      else tc[i]=t[nt-2];   
    }while((tc[i] < t[nt-2]) && (i<nhmax) && (hec[i]<hecmax));
    nhec=i++;  
  }
  else if (!xm){  // xm=0 is a CMP, hence he=h; h=0 is a zero offset trace

    tc[0]=0;//2*h/(sinbeta*velint[0]);
    hec[0]=h;
    nhec=1;
  }
  else return;

  tc[nhec]=t[nt-2];
  hec[nhec]=hecmax;
       
  it=int (tc[0]/dt+0.5);
  ihec=(int) ((hec[0]-he[0])/dhe);
  //int inithec=ihec; /* For printing */
  //fprintf(stderr,"ihec=%d\n",ihec);
  
  if (ihec>0)
    for (i=0;i<nhec;i++){
      itc=(int) (tc[i+1]/dt+0.5);
      if ((ihec < nhmax)&&(itc<nt)){    
	do{
	  it++;
	  //fprintf(stderr,"ihec=%d,it=%d,nhmax=%d,nt=%d\n",ihec,it,nhmax,nt);
	  csp[ihec][it]+=tr.data[it];
	  cspfold[ihec][it]=cspfold[ihec][it]+( int) fabs(scalefold*tr.data[it]);
	  //if (fabs(tr.data[it])>scalefold) cspfold[ihec][it]++;	    
	}while(it<=itc);
	ihec++;
      }
      else break;
    }
  //fprintf(stderr,"xm=%f, h=%f inithec=%d, nhmax=%d\n", xm, h, inithec, nhmax);
  return;
}  

void equiv_offset_1side(float *trace, int offset, int cdp, float **csp, unsigned short **cspfold, float *t,float xcsp,float *velint,float dt,float *he, float dhe, int nt, int nhmax, float *tc, float *hec, int precise, float beta, float cdpspace, float scalefold)
{
  //float  *hec;           /* equivalent offset for every trace */
  //float  *tc;            /* time when equivalent offset bin increases */

  int   ihec;   /* equivalent offset index */
  int   itc;    /* corresponding index for tc[i] */   
  int   it0;     /* corresponding index for t0[i] */   
  float hec2;   /* he*he  */
  int   nhec;   /* Number of offset bins at every trace */
  //float cdp;    /* midpoint coordinate */
  float xm;     /* midpoint distance from trace to scp = cdp_tr - x_scp */
  float xm2;             // xm*xm
  float h;               /* halfoffset */
  float h2;              // h*h
  float crossterm;       // 4*xm^2*h*2/vel
  float hecmax;          /* Maximum possible eq.  offset for this trace */
  float vel1;    // Temporal variable for velocity
  int i;
  register int it;

  //fprintf(stderr,"beta=%f\n",beta);
  float sinbeta=sin(beta);
  //fprintf(stderr,"sin(beta)=%f\n",sinbeta);

  for (i=0;i<nhmax;i++){
    tc[i]=0;
    hec[i]=0;
  }

  h=(float) abs(offset);
  h/=2;  // halfoffset
  //cdp=tr.sx+tr.offset/2.;
  
  //cdp=cdpnumber;

  xm=fabs(cdp-xcsp);
  // test for reconstructing the cdp from csp.
  // if (!xm) return;
  xm*=cdpspace;

  xm2=xm*xm;
  h2=h*h;
  hecmax=sqrt(xm2+h2);
  /* To eliminate particular offsets or midpoints use */
  //if (h==0) return; 
  //if (xm==0) return;

 // Beta ==0 is a code for processing only cdps
  if (beta==0) sinbeta=1;

  if ( xm && beta ){
    i=0;
    tc[i]=2*MAX(xm,h)/(sinbeta*velint[0]);  // First time to use for a ray coming horizontally
    //if (i==0 && tc[0]<0.1) fprintf(stderr,"xm=%f,h=%f\n",xm,h);
    itc=int (tc[i]/dt);
    it0=t2t0(t[MIN(itc,nt-1)],MAX(xm,h),velint[MIN(itc,nt-1)],dt);
    vel1=velint[MIN(it0,nt-1)];
    hec[i]=xm2+h2-4*xm2*h2/((tc[i]*vel1)*(tc[i]*vel1));

    if (hec[i] > 0) hec[i]=sqrt(hec[i]);
    else hec[i]=0;
    
    hecmax=sqrt(xm2+h2);
 
    //fprintf(stderr,"cdp=%f,xcsp=%f\n",cdp,xcsp);    
    do{
      i++;
      hec[i]=hec[i-1]+dhe;
      hec2=hec[i]*hec[i];
      // xm2+h2 would be enough for horizontal layers
      crossterm=xm2+h2-hec2;
      if (crossterm>0) {
	crossterm=sqrt(crossterm);
	tc[i]=2*xm*h/(vel1*crossterm);
	itc=(int) (tc[i]/dt+0.5);
	it0=t2t0(t[MIN(itc,nt-1)],hec[i-1],velint[MIN(itc,nt-1)],dt);
	vel1=velint[MIN(it0,nt-1)];
      }
      else tc[i]=t[nt-2];   
    }while((tc[i] < t[nt-2]) && (i<nhmax) && (hec[i]<hecmax));
    nhec=i++;  
  }
  else if (!xm){  // xm=0 is a CMP, hence he=h; h=0 is a zero offset trace

    tc[0]=0;//2*h/(sinbeta*velint[0]);
    hec[0]=h;
    nhec=1;
  }
  else return;

  tc[nhec]=t[nt-2];
  hec[nhec]=hecmax;
       
  it=int (tc[0]/dt+0.5);
  ihec=(int) ((hec[0]-he[0])/dhe);
  //int inithec=ihec; /* For printing */
  //fprintf(stderr,"ihec=%d\n",ihec);
  
  if (ihec>0)
    for (i=0;i<nhec;i++){
      itc=(int) (tc[i+1]/dt+0.5);
      if ((ihec < nhmax)&&(itc<nt)){    
	do{
	  it++;
	  //fprintf(stderr,"ihec=%d,it=%d,nhmax=%d,nt=%d\n",ihec,it,nhmax,nt);
	  csp[ihec][it]+=trace[it];
	  cspfold[ihec][it]=cspfold[ihec][it]+( int) fabs(scalefold*trace[it]);
	  //if (fabs(tr.data[it])>scalefold) cspfold[ihec][it]++;	    
	}while(it<=itc);
	ihec++;
      }
      else break;
    }
  //fprintf(stderr,"xm=%f, h=%f inithec=%d, nhmax=%d\n", xm, h, inithec, nhmax);
  return;
}  

#define MINTIME 10
  
void equiv_offset_test(segy tr,float **csp,unsigned short **cspfold,float *t,float xcsp,float *velint,float dt,float *he, float dhe, int nt, int nhmax, float beta, float cdpspace,
float scalefold)
{


  int it00;
  int it0;      /* Index */
  float hei;    /* equivalent offset for every trace */
  int   ihe ;   /* equivalent offset index */   
  float cdp;    /* midpoint coordinate */
  float xm;     /* midpoint distance from trace to scp = cdp_tr - x_scp */
  float xm2;             // xm*xm
  float h;               /* halfoffset */
  float h2;              // h*h
  float vel1;    // Temporal variable for velocity
  register int it;
  float t00;
  float sinbeta=sin(beta);
  int sign;
  int verbose2=0; 

  if (tr.offset > 0) sign=1;
  else sign=-1;
  if (verbose2) fprintf(stderr,"dhe=%f\n",dhe);
  
  h=fabs((float) tr.offset);
  h/=2.0;  // halfoffset
  //cdp=tr.sx+tr.offset/2.;

  cdp=tr.cdp;
  xm=fabs(cdp-xcsp);
  xm*=cdpspace;
  xm2=xm*xm;
  h2=h*h;
  // First time to use for a ray coming horizontally
  if (beta!=0) t00=2*MAX(xm,h)/(sinbeta*velint[0]);
  else t00=2*h/(velint[0]);
  
  it00=(int) (t00/dt+0.5);
  if (it00>=nt) it00=nt-1;
  it0=t2t0(t[MIN(it00,nt-1)],MAX(xm,h),velint[MIN(it00,nt-1)],dt);
  vel1=velint[MIN(it0,nt-1)];

  hei=xm2+h2-4*xm2*h2/((t[it00]*vel1)*(t[it00]*vel1));
  if (hei > 0) hei=sign*sqrt(hei);
  else hei=0;
  ihe=(int) ((hei-he[0])/dhe+sign*0.5);
  if (ihe>=0 &&  ihe < nhmax){
    it00=MAX(it00,MINTIME);
    for (it=it00;it< nt;it++){
      // For every t[it] estimate vel0 and equiv offset
      // The velocity must be estimated at the apex, i.e., t0
      it0=t2t0(t[MIN(it,nt-1)],hei,velint[MIN(it,nt-1)],dt);
      vel1=velint[MIN(it0,nt-1)];
      // xm2+h2 would be enough for horizontal layers
      // the crossterm corrects for xm != 0
      hei=xm2+h2-4*xm2*h2/((t[it]*vel1)*(t[it]*vel1));
      if (hei > 0) hei=sign*sqrt(hei);
      else hei=0;
      ihe=(int) ((hei-he[0])/dhe+sign*0.5);
      if (ihe >=0 && ihe < nhmax){
	csp[ihe][it]+=tr.data[it];
	//	cspfold[ihe][it]++;
	cspfold[ihe][it]=cspfold[ihe][it]+(int) fabs(scalefold*tr.data[it]);
      }      
    }
  }
  return;
}  

void equiv_offset_test(float *trace, int offset, int cdp, float **csp,unsigned short **cspfold,float *t,float xcsp,float *velint,float dt,float *he, float dhe, int nt, int nhmax, float beta, float cdpspace,
float scalefold)
{


  int it00;
  int it0;      /* Index */
  float hei;    /* equivalent offset for every trace */
  int   ihe ;   /* equivalent offset index */   
  //  float cdp;    /* midpoint coordinate */
  float xm;     /* midpoint distance from trace to scp = cdp_tr - x_scp */
  float xm2;             // xm*xm
  float h;               /* halfoffset */
  float h2;              // h*h
  float vel1;    // Temporal variable for velocity
  register int it;
  float t00;
  float sinbeta=sin(beta);
  int sign;
  int verbose2=0; 

  if (offset > 0) sign=1;
  else sign=-1;
  if (verbose2) fprintf(stderr,"dhe=%f\n",dhe);
  
  h=fabs((float) offset);
  h/=2.0;  // halfoffset
  //cdp=tr.sx+tr.offset/2.;

  xm=fabs(cdp-xcsp);
  xm*=cdpspace;
  xm2=xm*xm;
  h2=h*h;
  // First time to use for a ray coming horizontally
  if (beta!=0) t00=2*MAX(xm,h)/(sinbeta*velint[0]);
  else t00=2*h/(velint[0]);
  
  it00=(int) (t00/dt+0.5);
  if (it00>=nt) it00=nt-1;
  it0=t2t0(t[MIN(it00,nt-1)],MAX(xm,h),velint[MIN(it00,nt-1)],dt);
  vel1=velint[MIN(it0,nt-1)];

  hei=xm2+h2-4*xm2*h2/((t[it00]*vel1)*(t[it00]*vel1));
  if (hei > 0) hei=sign*sqrt(hei);
  else hei=0;
  ihe=(int) ((hei-he[0])/dhe+sign*0.5);
  if (ihe>=0 &&  ihe < nhmax){
    it00=MAX(it00,MINTIME);
    for (it=it00;it< nt;it++){
      // For every t[it] estimate vel0 and equiv offset
      // The velocity must be estimated at the apex, i.e., t0
      it0=t2t0(t[MIN(it,nt-1)],hei,velint[MIN(it,nt-1)],dt);
      vel1=velint[MIN(it0,nt-1)];
      // xm2+h2 would be enough for horizontal layers
      // the crossterm corrects for xm != 0
      hei=xm2+h2-4*xm2*h2/((t[it]*vel1)*(t[it]*vel1));
      if (hei > 0) hei=sign*sqrt(hei);
      else hei=0;
      ihe=(int) ((hei-he[0])/dhe+sign*0.5);
      if (ihe >=0 && ihe < nhmax){
	csp[ihe][it]+=trace[it];
	//	cspfold[ihe][it]++;
	cspfold[ihe][it]=cspfold[ihe][it]+(int) fabs(scalefold*trace[it]);
      }      
    }
  }
  return;
}  

    











