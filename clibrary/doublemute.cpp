#include "dan.h"
#include "su.h"

/* the minimum time to mute follows a slope give by it0-(iq-nmute)*slope
   This slope is a function of nt/nq. In general values between 2-5 must work fine */

#define SLOPE 2

void nmo(float *d,float *m,float *t,float h,float *vel,int adj,int nt,float dt,
	 float smute);

int taper(float **data, int nt, int nh, int ntaper,int flag);

void hrrti(float **d,float *pos,float dt,float **m,float *q, float fmax, 
	   int nt, int  nh, int nq, int rtmethod, float depth);

void doublemute(float **data, float **model, int nq, int nh, int nt, float *q, float *h, float *t, float parmute1, float parmute2, int mute1, int mute2, float fmax, int rtmethod, float depth, float t0)

{

  int side=0; // side 2 --> right mute ; side 1 left mute
  float **datarec=0;
  float dpd;
  float dpdp;
  float scale;
  float **modeltemp=0;
  int iq,ih,it;
  int nmute=0;
  float dt=t[1]-t[0];
  int it0=(int) (t0/dt); // define a minimum time for mute2
  int slope=SLOPE; // Defines a slope in the muting

  modeltemp=ealloc2float(nt,nq);
  datarec=ealloc2float(nt,nh);
  
  if (mute1){
    fprintf(stderr,"Muting artifacts before q=%f \n",parmute1);
    side=1;  // For artifacts use side 1 (filtering) 
    iq=0; while(q[iq]<parmute1) iq++; 
    nmute=iq;
    fprintf(stderr,"MUTING 1 at nmute=%d************************\n",nmute);
    for (iq=0;iq<nq;iq++)
      memcpy(modeltemp[iq],model[iq],nt*FSIZE);
    taper(modeltemp,nt,nq,nmute,side);
    if (0) plotgather_pipe(modeltemp,nq,nt,"mute1");
    hrrti(data,h,dt,modeltemp,q,fmax,nt,nh,nq,rtmethod,depth);
    if (0) plotgather_pipe(data,nh,nt,"mute1");
  }
  

  if (mute2){
    side=1;
    iq=0; while(q[iq]<parmute2) iq++; 
    if (side==2) nmute=nq-iq;
    else if (side==1) nmute=iq;
    fprintf(stderr,"MUTING 2 at nmute=%d************************\n",nmute);
    for (iq=0;iq<nq;iq++)
      memcpy(modeltemp[iq],model[iq],nt*FSIZE);
    taper(modeltemp,nt,nq,nmute,side);
    for (iq=0;iq<nq;iq++) for (it=0;it<(it0-(iq-nmute)*slope);it++) modeltemp[iq][it]=0; 
    if (1){
      save_gather(modeltemp,nq,q,nt,dt,"mutedmodel.su"); 
      system("cat mutedmodel.su >> radondata");      
      system("suxwigb key=f2 < mutedmodel.su title=mutedmodel perc=99 & \n");
    } 
  
  
    hrrti(datarec,h,dt,modeltemp,q,fmax,nt,nh,nq,rtmethod,depth);
    if (0){
      save_gather(datarec,nh,h,nt,dt,"multiples.su");      
      system("suxwigb key=offset < multiples.su title=nmo_multiples perc=97 & \n");
    } 
    if ((mute1)||(mute2)){
      fprintf(stderr,"substracting multiples \n");
      dpd=dot(nh*nt,data[0],datarec[0]);
      dpdp=dot(nh*nt,datarec[0],datarec[0]);
      scale=dpd/dpdp;
      fprintf(stderr,"scale 2===>%f\n",scale);
      for (ih=0;ih<nh;ih++) 
	for (it=0;it<nt;it++)
	     data[ih][it]=data[ih][it]-scale*datarec[ih][it];
    }
  }
  free2float(modeltemp);
  free2float(datarec);
  return;
}



 
 


